// Summarize E72Physics beam scalers for yield normalization.
//
// Run:
//   root -l -b -q kbeam_statistics.cc
//
// Outputs:
//   kbeam_statistics.txt  : run-by-run rows followed by the momentum summary
//
// The primary beam count is the scaler K-Beam counter. TRIG-D is retained
// only as a diagnostic because it matches K-Beam only for some run periods.
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TError.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TStyle.h>

namespace {
const std::string kRunListDirectory =
  "/home/had/haein/work/e72/ana/e72/runmanager/runlist/tpc_runinfo";
const std::string kScalerDirectory =
  "/gpfs/group/had/sks/E72/JPARC2025Nov/share/scaler";
const std::string kOutput = "kbeam_statistics.txt";
const std::string kPdfOutput = "kbeam_statistics.pdf";
const std::string kMomentumBinnedOutput = "kbeam_momentum_binned.txt";
const std::string kBeamProfileDirectory =
  "/gpfs/group/had/sks/E72/JPARC2025Nov/beam_simul";
// One representative beam profile per central-momentum setting.  For 735,
// use run 2270 as requested; edit this table if another reference is wanted.
const std::map<Int_t, std::string> kBeamProfileFiles{{
  {645, "beam_profile_run02978_-645.root"},
  {665, "beam_profile_run03006_-665.root"},
  {685, "beam_profile_run02883_-685.root"},
  {715, "beam_profile_run02629_-715.root"},
  {735, "beam_profile_run02270_-735.root"},
  {755, "beam_profile_run02723_-755.root"},
  {790, "beam_profile_run02800_-790.root"},
  {814, "beam_profile_run02912_-814.root"},
  {842, "beam_profile_run03031_-842.root"},
  {870, "beam_profile_run02992_-870.root"},
  {933, "beam_profile_run02856_-933.root"}
}};
constexpr Double_t kMomentumBinWidth = 0.004; // GeV/c
constexpr Double_t kMomentumBinMin = 0.50;    // GeV/c
constexpr Int_t kMomentumBins = 150;          // 0.50--1.10 GeV/c

struct ScalerCounts {
  Long64_t k_beam = -1;
  Long64_t trig_d = -1;
  double daq_efficiency = -1.;
};

struct RunRow {
  Int_t momentum_mev = 0;
  Int_t run = 0;
  std::string run_list_file;
  ScalerCounts scaler;
};

std::string Trim(const std::string& input)
{
  const auto first = input.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return "";
  const auto last = input.find_last_not_of(" \t\r\n");
  return input.substr(first, last-first+1);
}

double NormalizeEfficiency(double value)
{
  return value > 1. ? value/100. : value;
}

ScalerCounts ReadScaler(Int_t run)
{
  ScalerCounts result;
  const TString file_name = Form("%s/scaler_%05d.txt", kScalerDirectory.c_str(), run);
  std::ifstream input(file_name.Data());
  if (!input) {
    Warning("kbeam_statistics", "Cannot open %s", file_name.Data());
    return result;
  }

  std::string line;
  while (std::getline(input, line)) {
    std::istringstream stream(line);
    std::string name;
    double value = 0.;
    if (!(stream >> name >> value)) continue;
    if (name == "K-Beam") result.k_beam = static_cast<Long64_t>(value);
    if (name == "TRIG-D") result.trig_d = static_cast<Long64_t>(value);
    if (name == "DAQ-Eff") result.daq_efficiency = NormalizeEfficiency(value);
  }
  return result;
}

std::vector<Int_t> ReadE72PhysicsRuns(const std::filesystem::path& file)
{
  std::vector<Int_t> runs;
  std::ifstream input(file);
  std::string line;
  while (std::getline(input, line)) {
    const std::string trimmed = Trim(line);
    if (trimmed.empty() || trimmed.front() == '#') continue;
    // Keep only actual E72Physics entries; explicitly reject pedestal lines.
    if (trimmed.find("E72Physics") == std::string::npos ||
        trimmed.find("Pedestal") != std::string::npos) continue;
    std::istringstream stream(trimmed);
    Int_t run = 0;
    if (stream >> run) runs.push_back(run);
  }
  std::sort(runs.begin(), runs.end());
  runs.erase(std::unique(runs.begin(), runs.end()), runs.end());
  return runs;
}

bool AddScaledBeamProfile(const Int_t momentum_mev, const Double_t k_beam_total,
                          const Double_t effective_k_beam_total,
                          TH1D& k_beam_spectrum, TH1D& effective_spectrum)
{
  const auto profile_it = kBeamProfileFiles.find(momentum_mev);
  if (profile_it == kBeamProfileFiles.end()) {
    Warning("kbeam_statistics", "No beam profile configured for %d MeV/c", momentum_mev);
    return false;
  }
  const TString file_name = Form("%s/%s", kBeamProfileDirectory.c_str(),
                                 profile_it->second.c_str());
  TFile input(file_name, "READ");
  auto* tree = input.IsZombie() ? nullptr : dynamic_cast<TTree*>(input.Get("tr"));
  if (!tree || !tree->GetBranch("pInx") || !tree->GetBranch("pIny") ||
      !tree->GetBranch("pInz")) {
    Warning("kbeam_statistics", "Cannot read pInx/pIny/pInz from %s", file_name.Data());
    return false;
  }

  TH1D profile(Form("h_profile_%d", momentum_mev), "", kMomentumBins,
               kMomentumBinMin, kMomentumBinMin+kMomentumBins*kMomentumBinWidth);
  profile.SetDirectory(nullptr);
  TTreeReader reader(tree);
  TTreeReaderValue<Double_t> px(reader, "pInx");
  TTreeReaderValue<Double_t> py(reader, "pIny");
  TTreeReaderValue<Double_t> pz(reader, "pInz");
  while (reader.Next()) profile.Fill(std::sqrt((*px)*(*px) + (*py)*(*py) + (*pz)*(*pz)));
  const Double_t entries_in_range = profile.Integral(1, profile.GetNbinsX());
  Info("kbeam_statistics", "beam profile p=%d MeV/c (%s): mean=%.3f MeV/c, RMS=%.3f MeV/c, entries=%.0f",
       momentum_mev, profile_it->second.c_str(), 1000.*profile.GetMean(),
       1000.*profile.GetRMS(), entries_in_range);
  if (entries_in_range <= 0.) {
    Warning("kbeam_statistics", "No profile entries in the configured momentum range for %s", file_name.Data());
    return false;
  }
  k_beam_spectrum.Add(&profile, k_beam_total/entries_in_range);
  effective_spectrum.Add(&profile, effective_k_beam_total/entries_in_range);
  return true;
}
}

void kbeam_statistics()
{
  namespace fs = std::filesystem;
  const std::regex run_list_pattern(R"(E72Physics_([0-9]+)_[0-9]+\.txt)");
  std::map<Int_t, std::set<Int_t>> runs_by_momentum;
  std::map<Int_t, std::map<Int_t, std::string>> source_by_momentum_and_run;

  for (const auto& entry : fs::directory_iterator(kRunListDirectory)) {
    if (!entry.is_regular_file()) continue;
    std::smatch match;
    const std::string name = entry.path().filename().string();
    if (!std::regex_match(name, match, run_list_pattern)) continue;
    const Int_t momentum_mev = std::stoi(match[1].str());
    for (const Int_t run : ReadE72PhysicsRuns(entry.path())) {
      runs_by_momentum[momentum_mev].insert(run);
      source_by_momentum_and_run[momentum_mev][run] = name;
    }
  }

  if (runs_by_momentum.empty()) {
    Error("kbeam_statistics", "No E72Physics_<momentum>_<sequence>.txt files found in %s",
          kRunListDirectory.c_str());
    return;
  }

  std::vector<RunRow> rows;
  for (const auto& [momentum_mev, runs] : runs_by_momentum) {
    for (const Int_t run : runs) {
      RunRow row;
      row.momentum_mev = momentum_mev;
      row.run = run;
      row.run_list_file = source_by_momentum_and_run[momentum_mev][run];
      row.scaler = ReadScaler(run);
      rows.push_back(row);
    }
  }

  std::vector<RunRow> rows_by_run = rows;
  std::sort(rows_by_run.begin(), rows_by_run.end(),
            [](const RunRow& left, const RunRow& right) { return left.run < right.run; });

  std::ofstream run_output(kOutput);
  run_output << "# E72Physics physics runs only; source run lists: "
             << kRunListDirectory << "\n";
  run_output << "# K-beam normalization: scaler K-Beam; DAQ-Eff is K-Beam weighted.\n";
  run_output << "# TRIG-D is retained only as a diagnostic scaler counter.\n";
  run_output << "# rows are globally sorted by run number.\n";
  run_output << "# momentum_mev run nKbeam_K-Beam TRIG-D DAQ-Eff K-BeamxDAQ-Eff run_list\n";
  run_output << std::fixed << std::setprecision(6);
  for (const auto& row : rows_by_run) {
    const double effective_beam =
      (row.scaler.k_beam >= 0 && row.scaler.daq_efficiency >= 0.)
      ? row.scaler.k_beam*row.scaler.daq_efficiency : -1.;
    run_output << std::setw(5) << row.momentum_mev << " "
               << std::setw(5) << row.run << " "
               << std::setw(12) << row.scaler.k_beam << " "
               << std::setw(12) << row.scaler.trig_d << " "
               << std::setw(9) << row.scaler.daq_efficiency << " "
               << std::setw(14) << effective_beam << " "
               << row.run_list_file << "\n";
  }

  std::ofstream momentum_output(kOutput, std::ios::app);
  momentum_output << "\n[momentum_summary]\n";
  momentum_output << "# E72Physics trigger runs only; K-Beam sums and K-Beam-weighted DAQ efficiency.\n";
  momentum_output << "# momentum_mev n_runs nKbeam_sum_K-Beam TRIG-D_sum "
                  << "DAQ-Eff_weighted K-BeamxDAQ-Eff\n";
  momentum_output << std::fixed << std::setprecision(6);
  std::map<Int_t, Double_t> k_beam_total_by_momentum;
  std::map<Int_t, Double_t> effective_k_beam_total_by_momentum;
  for (const auto& [momentum_mev, runs] : runs_by_momentum) {
    Long64_t k_beam_sum = 0;
    Long64_t trig_d_sum = 0;
    double weighted_daq_sum = 0.;
    for (const auto& row : rows) {
      if (row.momentum_mev != momentum_mev) continue;
      if (row.scaler.trig_d >= 0) trig_d_sum += row.scaler.trig_d;
      if (row.scaler.k_beam >= 0 && row.scaler.daq_efficiency >= 0.) {
        k_beam_sum += row.scaler.k_beam;
        weighted_daq_sum += row.scaler.k_beam*row.scaler.daq_efficiency;
      }
    }
    const double daq_eff_weighted = k_beam_sum > 0 ? weighted_daq_sum/k_beam_sum : -1.;
    momentum_output << std::setw(5) << momentum_mev << " "
                    << std::setw(6) << runs.size() << " "
                    << std::setw(12) << k_beam_sum << " "
                    << std::setw(12) << trig_d_sum << " "
                    << std::setw(16) << daq_eff_weighted << " "
                    << std::setw(14) << weighted_daq_sum << "\n";
    k_beam_total_by_momentum[momentum_mev] = k_beam_sum;
    effective_k_beam_total_by_momentum[momentum_mev] = weighted_daq_sum;
    Info("kbeam_statistics",
         "p=%d MeV/c: %zu runs, K-Beam=%lld, TRIG-D=%lld, weighted DAQ-Eff=%.6f",
         momentum_mev, runs.size(), k_beam_sum, trig_d_sum, daq_eff_weighted);
  }

  TH1D k_beam_binned("h_kbeam_by_actual_momentum", "", kMomentumBins,
                     kMomentumBinMin, kMomentumBinMin+kMomentumBins*kMomentumBinWidth);
  TH1D effective_k_beam_binned("h_kbeam_daq_by_actual_momentum", "", kMomentumBins,
                               kMomentumBinMin, kMomentumBinMin+kMomentumBins*kMomentumBinWidth);
  k_beam_binned.SetDirectory(nullptr); effective_k_beam_binned.SetDirectory(nullptr);
  k_beam_binned.Sumw2(); effective_k_beam_binned.Sumw2();
  for (const auto& [momentum_mev, k_beam_total] : k_beam_total_by_momentum)
    AddScaledBeamProfile(momentum_mev, k_beam_total,
                         effective_k_beam_total_by_momentum[momentum_mev],
                         k_beam_binned, effective_k_beam_binned);

  std::ofstream binned_output(kMomentumBinnedOutput);
  binned_output << "# K-Beam totals distributed over actual K- momentum using beam_simul profiles.\n";
  binned_output << "# 735 MeV/c profile: beam_profile_run02270_-735.root\n";
  binned_output << "# momentum_low_GeV momentum_high_GeV momentum_center_GeV K-Beam_total K-BeamxDAQ-Eff_total\n";
  binned_output << std::fixed << std::setprecision(6);
  for (Int_t bin=1; bin<=kMomentumBins; ++bin) {
    binned_output << std::setw(10) << k_beam_binned.GetXaxis()->GetBinLowEdge(bin) << " "
                  << std::setw(10) << k_beam_binned.GetXaxis()->GetBinUpEdge(bin) << " "
                  << std::setw(10) << k_beam_binned.GetBinCenter(bin) << " "
                  << std::setw(16) << k_beam_binned.GetBinContent(bin) << " "
                  << std::setw(16) << effective_k_beam_binned.GetBinContent(bin) << "\n";
  }

  gStyle->SetOptStat(0);
  TGraph run_nkbeam;
  TGraph run_daq_efficiency;
  Int_t n_run_points = 0;
  Int_t n_daq_points = 0;
  for (const auto& row : rows_by_run) {
    if (row.scaler.k_beam >= 0)
      run_nkbeam.SetPoint(n_run_points++, row.run, row.scaler.k_beam);
    if (row.scaler.daq_efficiency >= 0.)
      run_daq_efficiency.SetPoint(n_daq_points++, row.run, row.scaler.daq_efficiency);
  }

  TGraph momentum_nkbeam;
  TGraph momentum_effective_nkbeam;
  TGraph momentum_daq_efficiency;
  Int_t n_momentum_points = 0;
  for (const auto& [momentum_mev, runs] : runs_by_momentum) {
    Long64_t n_kbeam = 0;
    double effective_nkbeam = 0.;
    for (const auto& row : rows) {
      if (row.momentum_mev != momentum_mev || row.scaler.k_beam < 0 ||
          row.scaler.daq_efficiency < 0.) continue;
      n_kbeam += row.scaler.k_beam;
      effective_nkbeam += row.scaler.k_beam*row.scaler.daq_efficiency;
    }
    if (n_kbeam <= 0) continue;
    momentum_nkbeam.SetPoint(n_momentum_points, momentum_mev, n_kbeam);
    momentum_effective_nkbeam.SetPoint(n_momentum_points, momentum_mev, effective_nkbeam);
    momentum_daq_efficiency.SetPoint(n_momentum_points, momentum_mev, effective_nkbeam/n_kbeam);
    ++n_momentum_points;
  }

  TCanvas canvas("c_kbeam_statistics", "K-beam statistics", 1100, 800);
  canvas.Print((kPdfOutput + "[").c_str());
  run_nkbeam.SetTitle("E72Physics run-by-run K-beam count;Run number;n_{K beam} = K-Beam");
  run_nkbeam.SetMarkerStyle(20); run_nkbeam.SetMarkerSize(.55);
  run_nkbeam.SetMarkerColor(kBlue + 1); run_nkbeam.Draw("AP");
  canvas.Print(kPdfOutput.c_str());

  canvas.Clear();
  run_daq_efficiency.SetTitle("E72Physics run-by-run DAQ efficiency;Run number;DAQ efficiency");
  run_daq_efficiency.SetMarkerStyle(20); run_daq_efficiency.SetMarkerSize(.55);
  run_daq_efficiency.SetMarkerColor(kRed + 1); run_daq_efficiency.SetMinimum(0.); run_daq_efficiency.SetMaximum(1.05);
  run_daq_efficiency.Draw("AP");
  canvas.Print(kPdfOutput.c_str());

  canvas.Clear();
  momentum_nkbeam.SetTitle("E72Physics K-beam statistics by central momentum;K^{-} beam momentum [MeV/c];n_{K beam} = #Sigma K-Beam");
  momentum_nkbeam.SetMarkerStyle(20); momentum_nkbeam.SetMarkerColor(kBlue + 1); momentum_nkbeam.SetLineColor(kBlue + 1);
  momentum_nkbeam.Draw("APL");
  momentum_effective_nkbeam.SetMarkerStyle(21); momentum_effective_nkbeam.SetMarkerColor(kRed + 1); momentum_effective_nkbeam.SetLineColor(kRed + 1);
  momentum_effective_nkbeam.Draw("PL same");
  TLegend beam_legend(.58, .70, .88, .87);
  beam_legend.SetBorderSize(0); beam_legend.SetFillStyle(0);
  beam_legend.AddEntry(&momentum_nkbeam, "#Sigma K-Beam (n_{K beam})", "pl");
  beam_legend.AddEntry(&momentum_effective_nkbeam, "#Sigma K-Beam #times DAQ-Eff", "pl");
  beam_legend.Draw();
  canvas.Print(kPdfOutput.c_str());

  canvas.Clear();
  momentum_daq_efficiency.SetTitle("K-Beam-weighted DAQ efficiency by central momentum;K^{-} beam momentum [MeV/c];DAQ efficiency");
  momentum_daq_efficiency.SetMarkerStyle(20); momentum_daq_efficiency.SetMarkerColor(kRed + 1); momentum_daq_efficiency.SetLineColor(kRed + 1);
  momentum_daq_efficiency.SetMinimum(0.); momentum_daq_efficiency.SetMaximum(1.05);
  momentum_daq_efficiency.Draw("APL");
  canvas.Print(kPdfOutput.c_str());

  canvas.Clear();
  k_beam_binned.SetTitle("K-Beam count versus actual K^{-} momentum (beam-profile weighted);|#vec{p}_{K^{-}}| [GeV/c];K-Beam count");
  k_beam_binned.SetLineColor(kBlue + 1); k_beam_binned.SetLineWidth(2);
  k_beam_binned.Draw("hist");
  canvas.Print(kPdfOutput.c_str());

  canvas.Clear();
  effective_k_beam_binned.SetTitle("DAQ-efficiency corrected K-Beam count versus actual K^{-} momentum;|#vec{p}_{K^{-}}| [GeV/c];K-Beam #times DAQ-Eff");
  effective_k_beam_binned.SetLineColor(kRed + 1); effective_k_beam_binned.SetLineWidth(2);
  effective_k_beam_binned.Draw("hist");
  canvas.Print(kPdfOutput.c_str());
  canvas.Print((kPdfOutput + "]").c_str());

  Info("kbeam_statistics", "Wrote %s, %s, and %s", kOutput.c_str(),
       kMomentumBinnedOutput.c_str(), kPdfOutput.c_str());
}
