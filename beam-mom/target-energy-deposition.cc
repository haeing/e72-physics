#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TStyle.h"

#include "../basic-property.hh"

namespace {
constexpr double kElectronMassMeV = 0.51099895;
constexpr double kBetheKMeVCm2PerMol = 0.307075;

struct Material {
  const char* name;
  double density_g_cm3;
  double thickness_mm;
  int layers;
  double z_over_a;
  double mean_excitation_eV;
  int color;
};

struct ParticleSpecies {
  const char* name;
  const char* label;
  double mass_mev_c2;
  int color;
  int marker;
};

std::array<ParticleSpecies, 3> Particles()
{
  return {{{"pion", "#pi^{#pm}", mpi, kBlue + 1, 20},
           {"kaon", "K^{#pm}", mK, kRed + 1, 21},
           {"proton", "proton", mp, kGreen + 2, 22}}};
}

// Compound Z/A and I values are effective material values.  Densities and
// single-layer thicknesses are taken from basic-property.hh.  The assumed
// reaction point is at the center of the LH2 volume, so only the upstream
// target wall (one GFRP and one Kapton layer, plus two Mylar layers) and half
// of the LH2 length contribute to the incoming-particle energy loss.
std::array<Material, 4> TargetMaterials()
{
  return {{{"GFRP", gfrp_density, gfrp_thick, 1, 0.50, 100.0,
            kOrange + 7},
           {"Kapton", kapton_density, kapton_thick, 1, 0.5126, 79.6,
            kGreen + 2},
           {"Mylar", mylar_density, mylar_thick, 2, 0.5204, 78.7,
            kMagenta + 1},
           {"LH2 (to center)", lh2_density, 0.5 * lh2_thick, 1,
            1.0 / 1.008, 21.8,
            kBlue + 1}}};
}

double MeanStoppingPower(double momentum_mev_c, double mass_mev_c2,
                         const Material& material)
{
  if (momentum_mev_c <= 0.0) return 0.0;

  const double energy = std::hypot(momentum_mev_c, mass_mev_c2);
  const double beta2 = momentum_mev_c * momentum_mev_c / (energy * energy);
  const double gamma = energy / mass_mev_c2;
  const double mass_ratio = kElectronMassMeV / mass_mev_c2;
  const double wmax = 2.0 * kElectronMassMeV * beta2 * gamma * gamma /
      (1.0 + 2.0 * gamma * mass_ratio + mass_ratio * mass_ratio);
  const double excitation_mev = material.mean_excitation_eV * 1.0e-6;
  const double log_term = 0.5 * std::log(
      2.0 * kElectronMassMeV * beta2 * gamma * gamma * wmax /
      (excitation_mev * excitation_mev));

  // Mean collisional stopping power in MeV/cm.  All three particles have
  // |z|=1, so the charge-squared factor is one.
  return material.density_g_cm3 * kBetheKMeVCm2PerMol *
         material.z_over_a / beta2 * (log_term - beta2);
}

double DepositInMaterial(double& kinetic_energy_mev, double mass_mev_c2,
                         const Material& material)
{
  const double path_cm = material.thickness_mm * material.layers / 10.0;
  const int steps = std::max(1, static_cast<int>(std::ceil(path_cm / 0.02)));
  const double step_cm = path_cm / steps;
  double deposited = 0.0;

  for (int step = 0; step < steps && kinetic_energy_mev > 0.0; ++step) {
    const double total_energy = kinetic_energy_mev + mass_mev_c2;
    const double momentum = std::sqrt(std::max(
        0.0, total_energy * total_energy - mass_mev_c2 * mass_mev_c2));
    const double loss = std::min(
        kinetic_energy_mev,
        MeanStoppingPower(momentum, mass_mev_c2, material) * step_cm);
    kinetic_energy_mev -= loss;
    deposited += loss;
  }
  return deposited;
}
}  // namespace

void target_energy_deposition(double momentum_min_mev_c = 300.0,
                              double momentum_max_mev_c = 1000.0,
                              double momentum_step_mev_c = 50.0,
                              const char* output_prefix =
                                  "target-energy-deposition")
{
  if (momentum_min_mev_c <= 0.0 || momentum_max_mev_c < momentum_min_mev_c ||
      momentum_step_mev_c <= 0.0) {
    std::cerr << "Invalid momentum range or step." << std::endl;
    return;
  }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const auto materials = TargetMaterials();
  const auto particles = Particles();
  std::vector<double> momenta;
  std::array<std::array<std::vector<double>, 4>, 3> deposits;
  std::array<std::vector<double>, 3> total_deposits;
  std::array<std::vector<double>, 3> exit_momenta;
  std::array<std::vector<double>, 3> momentum_losses;

  for (double momentum = momentum_min_mev_c;
       momentum <= momentum_max_mev_c + 0.5 * momentum_step_mev_c;
       momentum += momentum_step_mev_c) {
    momenta.push_back(momentum);
    for (std::size_t particle = 0; particle < particles.size(); ++particle) {
      const double mass = particles[particle].mass_mev_c2;
      double kinetic_energy = std::hypot(momentum, mass) - mass;
      double total_deposit = 0.0;
      for (std::size_t material = 0; material < materials.size(); ++material) {
        const double deposit =
            DepositInMaterial(kinetic_energy, mass, materials[material]);
        deposits[particle][material].push_back(deposit);
        total_deposit += deposit;
      }
      total_deposits[particle].push_back(total_deposit);
      const double exit_total_energy = kinetic_energy + mass;
      const double exit_momentum = std::sqrt(std::max(
          0.0, exit_total_energy * exit_total_energy - mass * mass));
      exit_momenta[particle].push_back(exit_momentum);
      momentum_losses[particle].push_back(momentum - exit_momentum);
    }
  }

  const std::string prefix(output_prefix);
  std::ofstream csv(prefix + ".csv");
  csv << "momentum_MeV_c";
  for (const auto& particle : particles) {
    for (const auto& material : materials)
      csv << ',' << particle.name << '_' << material.name << "_MeV";
    csv << ',' << particle.name << "_total_MeV," << particle.name
        << "_exit_momentum_MeV_c," << particle.name
        << "_momentum_loss_MeV_c";
  }
  csv << '\n';
  csv << std::fixed << std::setprecision(6);

  std::cout << std::fixed << std::setprecision(3)
            << " p [MeV/c]   pion: dE/dp    kaon: dE/dp  proton: dE/dp"
            << std::endl;
  for (std::size_t point = 0; point < momenta.size(); ++point) {
    csv << momenta[point];
    std::cout << std::setw(10) << momenta[point];
    for (std::size_t particle = 0; particle < particles.size(); ++particle) {
      for (std::size_t material = 0; material < materials.size(); ++material)
        csv << ',' << deposits[particle][material][point];
      csv << ',' << total_deposits[particle][point] << ','
          << exit_momenta[particle][point] << ','
          << momentum_losses[particle][point];
      std::cout << std::setw(9) << total_deposits[particle][point] << '/'
                << std::setw(7) << momentum_losses[particle][point];
    }
    csv << '\n';
    std::cout << std::endl;
  }

  TFile output((prefix + ".root").c_str(), "RECREATE");
  TCanvas canvas("c_target_energy_deposition", "Target energy deposition",
                 900, 1000);
  canvas.Divide(1, 2);
  TMultiGraph comparison;
  comparison.SetTitle(
      "Energy loss to the LH2 center;Initial momentum [MeV/#it{c}];Mean energy loss [MeV]");
  TLegend legend(0.66, 0.62, 0.88, 0.87);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);

  std::array<TGraph, 3> total_graphs;
  std::array<TGraph, 3> momentum_loss_graphs;
  TMultiGraph momentum_comparison;
  momentum_comparison.SetTitle(
      "Momentum loss to the LH2 center;Initial momentum [MeV/#it{c}];#Deltap [MeV/#it{c}]");
  TLegend momentum_legend(0.66, 0.62, 0.88, 0.87);
  momentum_legend.SetBorderSize(0);
  momentum_legend.SetFillStyle(0);
  for (std::size_t particle = 0; particle < particles.size(); ++particle) {
    total_graphs[particle] = TGraph(momenta.size(), momenta.data(),
                                    total_deposits[particle].data());
    total_graphs[particle].SetName((std::string("g_deposit_total_") +
                                    particles[particle].name).c_str());
    total_graphs[particle].SetLineColor(particles[particle].color);
    total_graphs[particle].SetMarkerColor(particles[particle].color);
    total_graphs[particle].SetLineWidth(3);
    total_graphs[particle].SetMarkerStyle(particles[particle].marker);
    comparison.Add(&total_graphs[particle], "LP");
    legend.AddEntry(&total_graphs[particle], particles[particle].label, "lp");

    momentum_loss_graphs[particle] =
        TGraph(momenta.size(), momenta.data(), momentum_losses[particle].data());
    momentum_loss_graphs[particle].SetName(
        (std::string("g_momentum_loss_") + particles[particle].name).c_str());
    momentum_loss_graphs[particle].SetLineColor(particles[particle].color);
    momentum_loss_graphs[particle].SetMarkerColor(particles[particle].color);
    momentum_loss_graphs[particle].SetLineWidth(3);
    momentum_loss_graphs[particle].SetMarkerStyle(particles[particle].marker);
    momentum_comparison.Add(&momentum_loss_graphs[particle], "LP");
    momentum_legend.AddEntry(&momentum_loss_graphs[particle],
                             particles[particle].label, "lp");
  }

  const std::string pdf_name = prefix + ".pdf";
  canvas.cd(1);
  comparison.Draw("A");
  legend.Draw();
  canvas.cd(2);
  momentum_comparison.Draw("A");
  momentum_legend.Draw();
  canvas.Print(pdf_name.c_str());
  for (auto& graph : total_graphs) graph.Write();
  for (auto& graph : momentum_loss_graphs) graph.Write();
  output.Close();

  std::cout << "Wrote " << pdf_name << ", " << prefix << ".root, and "
            << prefix << ".csv" << std::endl;
}
