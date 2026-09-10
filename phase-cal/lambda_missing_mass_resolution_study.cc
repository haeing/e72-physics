// Toy missing-mass resolution study.
// Edit the configuration below, then run:
//   root -l -b -q lambda_missing_mass_resolution_study.cc
#include <array>
#include <cmath>
#include <memory>
#include <vector>

#include <TCanvas.h>
#include <TGenPhaseSpace.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TString.h>
#include <TStyle.h>
#include <TPaveText.h>

namespace {
constexpr Double_t kKaonMass = 0.493677;
constexpr Double_t kProtonMass = 0.938272;
constexpr Double_t kPionMass = 0.139570;
constexpr Double_t kPi0Mass = 0.134977;
constexpr Double_t kLambdaMass = 1.115683;
constexpr Double_t kEtaMass = 0.547862;

// ---- Study configuration ---------------------------------------------------
constexpr Double_t kBeamMomentumMeV = 735.;
constexpr Int_t kEventsPerReaction = 200000;
constexpr UInt_t kRandomSeed = 20260907;
// Gaussian relative sigma of the momentum magnitude. Directions are retained.
// Set any value to zero for a truth-level reference.
constexpr Double_t kBeamRelativeMomentumSigma = 1.e-3;
constexpr Double_t kProtonRelativeMomentumSigma = 0.040;
constexpr Double_t kPionRelativeMomentumSigma = 0.040;
// Gaussian sigma of the reconstructed Lambda direction (opening-angle error).
constexpr Double_t kLambdaAngularSigmaDeg = 1.0;
// Gaussian-core fit half-widths.  M_X^2 is the more robust resolution metric
// for #pi^0 because its negative tail is retained rather than threshold-folded.
constexpr Double_t kMassFitHalfWidth = 0.070;      // GeV/c^2
constexpr Double_t kMass2FitHalfWidth = 0.100;     // (GeV/c^2)^2
// ---------------------------------------------------------------------------

struct Reaction {
  TString label;
  Double_t missing_mass;
  Color_t color;
  Double_t weight; // relative yield in the inclusive mixture
};

const std::array<Reaction, 2> kReactions{{
  {"#Lambda#eta", kEtaMass, kRed + 1, 1.},
  {"#Lambda#pi^{0}", kPi0Mass, kBlue + 1, 1.},
}};

TLorentzVector MissingFourMomentum(const TVector3& lambda_momentum,
                                   const TVector3& beam_momentum)
{
  TLorentzVector beam;
  beam.SetXYZM(beam_momentum.X(), beam_momentum.Y(), beam_momentum.Z(), kKaonMass);
  TLorentzVector target;
  target.SetXYZM(0., 0., 0., kProtonMass);
  // This is the same Lambda mass constraint used in DstTPCLambdaEta.
  TLorentzVector lambda;
  lambda.SetXYZM(lambda_momentum.X(), lambda_momentum.Y(), lambda_momentum.Z(), kLambdaMass);
  return beam + target - lambda;
}

TVector3 SmearMagnitude(const TVector3& momentum, Double_t relative_sigma,
                        TRandom3& random)
{
  if (relative_sigma <= 0. || momentum.Mag2() <= 0.) return momentum;
  return momentum * random.Gaus(1., relative_sigma);
}

// Apply an independent Gaussian opening-angle error while preserving |p|.
// The deflection azimuth around the original direction is uniform.
TVector3 SmearDirection(const TVector3& momentum, Double_t angular_sigma_deg,
                        TRandom3& random)
{
  if (angular_sigma_deg <= 0. || momentum.Mag2() <= 0.) return momentum;
  const TVector3 direction = momentum.Unit();
  const TVector3 reference = std::abs(direction.Z()) < 0.9
    ? TVector3(0., 0., 1.) : TVector3(1., 0., 0.);
  const TVector3 transverse_1 = reference.Cross(direction).Unit();
  const TVector3 transverse_2 = direction.Cross(transverse_1).Unit();
  const Double_t delta = random.Gaus(0., angular_sigma_deg*TMath::DegToRad());
  const Double_t azimuth = random.Uniform(0., 2.*TMath::Pi());
  const TVector3 smeared_direction = direction*TMath::Cos(delta)
    + (transverse_1*TMath::Cos(azimuth) + transverse_2*TMath::Sin(azimuth))*TMath::Sin(delta);
  return momentum.Mag()*smeared_direction;
}

void Style(TH1D& histogram, Color_t color)
{
  histogram.SetLineColor(color);
  histogram.SetLineWidth(2);
  histogram.SetFillStyle(0);
}
}

void lambda_missing_mass_resolution_study()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const Double_t beam_momentum = kBeamMomentumMeV/1000.;
  const TVector3 beam_truth(0., 0., beam_momentum);
  TLorentzVector beam_four;
  beam_four.SetXYZM(0., 0., beam_momentum, kKaonMass);
  TLorentzVector target;
  target.SetXYZM(0., 0., 0., kProtonMass);
  TLorentzVector initial_state = beam_four + target;
  const std::array<Double_t, 2> lambda_decay_masses{{kProtonMass, kPionMass}};

  auto truth_mass = std::make_unique<TH1D>("h_missing_mass_truth_all", "", 650, 0., 0.65);
  auto smeared_mass = std::make_unique<TH1D>("h_missing_mass_smeared_all", "", 650, 0., 0.65);
  auto truth_mass2 = std::make_unique<TH1D>("h_missing_mass2_truth_all", "", 550, -0.10, 0.45);
  auto smeared_mass2 = std::make_unique<TH1D>("h_missing_mass2_smeared_all", "", 550, -0.10, 0.45);
  for (auto* histogram : {truth_mass.get(), smeared_mass.get(), truth_mass2.get(), smeared_mass2.get()}) {
    histogram->SetDirectory(nullptr); histogram->Sumw2(); Style(*histogram, kBlack);
  }

  std::array<std::unique_ptr<TH1D>, kReactions.size()> truth_mass_by_reaction;
  std::array<std::unique_ptr<TH1D>, kReactions.size()> smeared_mass_by_reaction;
  std::array<std::unique_ptr<TH1D>, kReactions.size()> truth_mass2_by_reaction;
  std::array<std::unique_ptr<TH1D>, kReactions.size()> smeared_mass2_by_reaction;
  // Reconstructed Lambda lab phase space, using the same daughter smearing
  // as the missing-mass resolution study.
  std::array<std::unique_ptr<TH2D>, kReactions.size()> smeared_lambda_lab_phase_space;
  for (std::size_t i=0; i<kReactions.size(); ++i) {
    truth_mass_by_reaction[i] = std::make_unique<TH1D>(Form("h_missing_mass_truth_%zu", i), "", 650, 0., 0.65);
    smeared_mass_by_reaction[i] = std::make_unique<TH1D>(Form("h_missing_mass_smeared_%zu", i), "", 650, 0., 0.65);
    truth_mass2_by_reaction[i] = std::make_unique<TH1D>(Form("h_missing_mass2_truth_%zu", i), "", 550, -0.10, 0.45);
    smeared_mass2_by_reaction[i] = std::make_unique<TH1D>(Form("h_missing_mass2_smeared_%zu", i), "", 550, -0.10, 0.45);
    smeared_lambda_lab_phase_space[i] = std::make_unique<TH2D>(
      Form("h_lambda_lab_phase_space_smeared_%zu", i), "", 200, -1., 1., 200, 0., 1.2);
    smeared_lambda_lab_phase_space[i]->SetDirectory(nullptr);
    smeared_lambda_lab_phase_space[i]->Sumw2();
    for (auto* histogram : {truth_mass_by_reaction[i].get(), smeared_mass_by_reaction[i].get(),
                            truth_mass2_by_reaction[i].get(), smeared_mass2_by_reaction[i].get()}) {
      histogram->SetDirectory(nullptr); histogram->Sumw2(); Style(*histogram, kReactions[i].color);
    }
  }

  TRandom3 random(kRandomSeed);
  for (std::size_t reaction_id=0; reaction_id<kReactions.size(); ++reaction_id) {
    const std::array<Double_t, 2> production_masses{{kLambdaMass, kReactions[reaction_id].missing_mass}};
    TGenPhaseSpace production;
    if (!production.SetDecay(initial_state, production_masses.size(), production_masses.data())) continue;
    TGenPhaseSpace lambda_decay;
    for (Int_t event=0; event<kEventsPerReaction; ++event) {
      production.Generate();
      auto* lambda_truth = production.GetDecay(0);
      if (!lambda_decay.SetDecay(*lambda_truth, lambda_decay_masses.size(), lambda_decay_masses.data())) continue;
      lambda_decay.Generate();
      auto* proton_truth = lambda_decay.GetDecay(0);
      auto* pion_truth = lambda_decay.GetDecay(1);

      const auto truth_missing = MissingFourMomentum(lambda_truth->Vect(), beam_truth);
      const Double_t weight = kReactions[reaction_id].weight;
      truth_mass->Fill(truth_missing.M(), weight);
      truth_mass2->Fill(truth_missing.M2(), weight);
      truth_mass_by_reaction[reaction_id]->Fill(truth_missing.M(), weight);
      truth_mass2_by_reaction[reaction_id]->Fill(truth_missing.M2(), weight);

      const TVector3 proton_smeared = SmearMagnitude(proton_truth->Vect(), kProtonRelativeMomentumSigma, random);
      const TVector3 pion_smeared = SmearMagnitude(pion_truth->Vect(), kPionRelativeMomentumSigma, random);
      TLorentzVector proton_four, pion_four;
      proton_four.SetXYZM(proton_smeared.X(), proton_smeared.Y(), proton_smeared.Z(), kProtonMass);
      pion_four.SetXYZM(pion_smeared.X(), pion_smeared.Y(), pion_smeared.Z(), kPionMass);
      const TVector3 lambda_smeared = SmearDirection(
        (proton_four+pion_four).Vect(), kLambdaAngularSigmaDeg, random);
      const Double_t lambda_smeared_momentum = lambda_smeared.Mag();
      if (lambda_smeared_momentum > 0.) {
        smeared_lambda_lab_phase_space[reaction_id]->Fill(
          lambda_smeared.Z()/lambda_smeared_momentum, lambda_smeared_momentum, weight);
      }
      const TVector3 beam_smeared = SmearMagnitude(beam_truth, kBeamRelativeMomentumSigma, random);
      const auto smeared_missing = MissingFourMomentum(lambda_smeared, beam_smeared);
      smeared_mass->Fill(smeared_missing.M(), weight);
      smeared_mass2->Fill(smeared_missing.M2(), weight);
      smeared_mass_by_reaction[reaction_id]->Fill(smeared_missing.M(), weight);
      smeared_mass2_by_reaction[reaction_id]->Fill(smeared_missing.M2(), weight);
    }
  }

  std::array<std::unique_ptr<TF1>, kReactions.size()> mass_fits;
  std::array<std::unique_ptr<TF1>, kReactions.size()> mass2_fits;
  for (std::size_t i=0; i<kReactions.size(); ++i) {
    const Double_t mass_center = kReactions[i].missing_mass;
    const Double_t mass2_center = mass_center*mass_center;
    mass_fits[i] = std::make_unique<TF1>(Form("fit_missing_mass_%zu", i), "gaus",
                                         std::max(0., mass_center-kMassFitHalfWidth),
                                         std::min(0.65, mass_center+kMassFitHalfWidth));
    mass2_fits[i] = std::make_unique<TF1>(Form("fit_missing_mass2_%zu", i), "gaus",
                                          mass2_center-kMass2FitHalfWidth,
                                          mass2_center+kMass2FitHalfWidth);
    mass_fits[i]->SetLineColor(kReactions[i].color); mass_fits[i]->SetLineStyle(2); mass_fits[i]->SetLineWidth(2);
    mass2_fits[i]->SetLineColor(kReactions[i].color); mass2_fits[i]->SetLineStyle(2); mass2_fits[i]->SetLineWidth(2);
    smeared_mass_by_reaction[i]->Fit(mass_fits[i].get(), "RQ0");
    smeared_mass2_by_reaction[i]->Fit(mass2_fits[i].get(), "RQ0");
    Info("lambda_missing_mass_resolution_study",
         "%s: sigma(M_X) = %.6f GeV/c^2, sigma(M_X^2) = %.6f (GeV/c^2)^2",
         kReactions[i].label.Data(), mass_fits[i]->GetParameter(2), mass2_fits[i]->GetParameter(2));
  }

  const TString pdf_name = Form("lambda_missing_mass_resolution_mom%.0fMeV.pdf", kBeamMomentumMeV);
  TCanvas canvas("c_missing_mass_resolution", "Missing mass resolution study", 1000, 800);
  canvas.Print(pdf_name + "[");
  using FitArray = std::array<std::unique_ptr<TF1>, kReactions.size()>;
  auto draw_overlay = [&](TH1D& total, const auto& components, const TString& title,
                          const FitArray* fits, Bool_t show_fit) {
    total.SetTitle(title);
    total.SetMaximum(1.15*total.GetMaximum());
    total.Draw("hist");
    TLegend legend(.53, .66, .88, .88);
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(.027);
    legend.AddEntry(&total, "inclusive (#Lambda#eta + #Lambda#pi^{0})", "l");
    for (std::size_t i=0; i<kReactions.size(); ++i) {
      components[i]->Draw("hist same");
      legend.AddEntry(components[i].get(), kReactions[i].label, "l");
      if (show_fit) { (*fits)[i]->Draw("same"); legend.AddEntry((*fits)[i].get(), kReactions[i].label + " Gaussian core", "l"); }
    }
    legend.Draw();
    if (show_fit) {
      TPaveText label(.14, .70, .50, .88, "NDC");
      label.SetBorderSize(0); label.SetFillStyle(0); label.SetTextAlign(12); label.SetTextSize(.027);
      label.AddText("Gaussian-core #sigma:");
      for (std::size_t i=0; i<kReactions.size(); ++i) label.AddText(Form("%s  %.5f", kReactions[i].label.Data(), (*fits)[i]->GetParameter(2)));
      label.Draw();
    }
    canvas.Print(pdf_name);
    canvas.Clear();
  };
  draw_overlay(*truth_mass, truth_mass_by_reaction,
               "Truth missing mass (no resolution);M_{X}(#Lambda) [GeV/c^{2}];Events / 1 MeV/c^{2}", nullptr, false);
  draw_overlay(*smeared_mass, smeared_mass_by_reaction,
               Form("Smeared missing mass: #sigma_{p}/p = %.1f%% (p), %.1f%% (#pi), %.1f%% (K^{-});M_{X}(#Lambda) [GeV/c^{2}];Events / 1 MeV/c^{2}",
                    100.*kProtonRelativeMomentumSigma, 100.*kPionRelativeMomentumSigma, 100.*kBeamRelativeMomentumSigma), &mass_fits, true);
  draw_overlay(*truth_mass2, truth_mass2_by_reaction,
               "Truth missing mass squared (no resolution);M_{X}^{2}(#Lambda) [(GeV/c^{2})^{2}];Events / 0.001 (GeV/c^{2})^{2}", nullptr, false);
  draw_overlay(*smeared_mass2, smeared_mass2_by_reaction,
               "Smeared missing mass squared;M_{X}^{2}(#Lambda) [(GeV/c^{2})^{2}];Events / 0.001 (GeV/c^{2})^{2}", &mass2_fits, true);

  auto smeared_lambda_lab_phase_space_all = std::unique_ptr<TH2D>(
    dynamic_cast<TH2D*>(smeared_lambda_lab_phase_space[0]->Clone(
      "h_lambda_lab_phase_space_smeared_all")));
  smeared_lambda_lab_phase_space_all->SetDirectory(nullptr);
  smeared_lambda_lab_phase_space_all->Add(smeared_lambda_lab_phase_space[1].get());
  canvas.Clear();
  smeared_lambda_lab_phase_space_all->SetTitle(
    Form("Smeared #Lambda lab phase space: #Lambda#eta + #Lambda#pi^{0} (#sigma_{#alpha}^{#Lambda}=%.1f^{#circ});cos #theta_{#Lambda}^{lab};|#vec{p}_{#Lambda}^{lab}| [GeV/c]",
         kLambdaAngularSigmaDeg));
  smeared_lambda_lab_phase_space_all->Draw("colz");
  canvas.Print(pdf_name);
  canvas.Print(pdf_name + "]");
}
