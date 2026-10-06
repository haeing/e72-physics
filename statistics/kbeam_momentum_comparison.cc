// Compare the reconstructed K18 momentum distributions for 685, 715, and 735 MeV/c.
// Run from this directory with:
//   root -l -b -q kbeam_momentum_comparison.cc

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TDatime.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

namespace {
struct Setting {
  Int_t central_momentum_mev;
  Int_t input_run;
  const char* input_file;
};

struct MomentumStatistics {
  Long64_t k_beam = 0;
  Double_t effective_beam = 0.;
  std::vector<Int_t> runs;
};

const std::vector<Setting> kSettings = {
  {735, 2448, "/gpfs/home/had/haein/data/JPARC2025Nov_root/blc/run02448_K18Tracking.root"},
  {715, 2682, "/gpfs/home/had/haein/data/JPARC2025Nov_root/blc/run02682_K18Tracking.root"},
  {685, 2883, "/gpfs/home/had/haein/data/JPARC2025Nov_root/blc/run02883_K18Tracking.root"}
};
const std::string kStatisticsFile = "kbeam_statistics.txt";
const std::string kOutputDirectory = "result";
constexpr Int_t kMomentumBins = 220;
constexpr Double_t kMomentumMin = 0.67;
constexpr Double_t kMomentumMax = 0.80;

std::map<Int_t, MomentumStatistics> ReadMomentumStatistics()
{
  std::map<Int_t, MomentumStatistics> statistics;
  std::ifstream input(kStatisticsFile);
  if (!input) {
    Error("kbeam_momentum_comparison", "Cannot open %s", kStatisticsFile.c_str());
    return statistics;
  }

  std::string line;
  while (std::getline(input, line)) {
    std::istringstream row(line);
    Int_t momentum = 0, run = 0;
    Long64_t k_beam = 0, trig_d = 0;
    Double_t daq_eff = 0., effective_beam = 0.;
    std::string run_list;
    if (!(row >> momentum >> run >> k_beam >> trig_d >> daq_eff >> effective_beam >> run_list)) continue;
    if (k_beam <= 0 || daq_eff < 0. || effective_beam < 0.) continue;
    auto& total = statistics[momentum];
    total.k_beam += k_beam;
    total.effective_beam += effective_beam;
    total.runs.push_back(run);
  }
  return statistics;
}
}

void kbeam_momentum_comparison()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(kOutputDirectory.c_str(), kTRUE);
  const auto statistics = ReadMomentumStatistics();
  std::vector<std::unique_ptr<TH1D>> spectra;
  std::vector<Double_t> effective_beam;
  for (const auto& setting : kSettings) {
    auto histogram = std::make_unique<TH1D>(
      Form("h_pk18_%d_scaled", setting.central_momentum_mev), "",
      kMomentumBins, kMomentumMin, kMomentumMax);
    histogram->SetDirectory(nullptr);
    histogram->Sumw2();

    std::unique_ptr<TFile> input(TFile::Open(setting.input_file, "READ"));
    if (!input || input->IsZombie()) {
      Error("kbeam_momentum_comparison", "Cannot open %s", setting.input_file);
      spectra.push_back(std::move(histogram));
      effective_beam.push_back(0.);
      continue;
    }
    auto* tree = dynamic_cast<TTree*>(input->Get("k18"));
    if (!tree || !tree->GetBranch("pk18")) {
      Error("kbeam_momentum_comparison", "k18/pk18 is missing in %s", setting.input_file);
      spectra.push_back(std::move(histogram));
      effective_beam.push_back(0.);
      continue;
    }

    TTreeReader reader(tree);
    TTreeReaderValue<std::vector<Double_t>> pk18(reader, "pk18");
    Long64_t entries = 0;
    while (reader.Next()) {
      for (const Double_t momentum : *pk18) {
        if (!std::isfinite(momentum)) continue;
        histogram->Fill(momentum);
        ++entries;
      }
    }

    const auto stats_it = statistics.find(setting.central_momentum_mev);
    const Double_t effective = stats_it == statistics.end() ? 0. : stats_it->second.effective_beam;
    if (entries > 0 && effective > 0.) histogram->Scale(effective / entries);
    effective_beam.push_back(effective);
    Info("kbeam_momentum_comparison", "p=%d, shape run=%d, pk18 entries=%lld, K-Beam sum=%.0f, K-Beam*DAQ sum=%.3e",
         setting.central_momentum_mev, setting.input_run, entries,
         stats_it == statistics.end() ? 0. : static_cast<Double_t>(stats_it->second.k_beam), effective);
    spectra.push_back(std::move(histogram));
  }

  Int_t first_run = std::numeric_limits<Int_t>::max();
  Int_t last_run = 0;
  for (const auto& [momentum, stats] : statistics)
    for (const Int_t run : stats.runs) {
      if (std::find_if(kSettings.begin(), kSettings.end(), [momentum](const Setting& setting) {
            return setting.central_momentum_mev == momentum;
          }) == kSettings.end()) continue;
      first_run = std::min(first_run, run);
      last_run = std::max(last_run, run);
    }
  const std::string run_suffix = last_run > 0
    ? Form("_run%05d-run%05d", first_run, last_run) : "_no_valid_runs";
  const std::string pdf_name = kOutputDirectory + "/kbeam_momentum_comparison" + run_suffix + ".pdf";
  const std::string root_name = kOutputDirectory + "/kbeam_momentum_comparison" + run_suffix + ".root";

  TCanvas canvas("c_kbeam_momentum_comparison", "K-beam weighted K18 momentum", 1000, 800);
  canvas.Print((pdf_name + "[").c_str());
  TPaveText title(.07, .80, .93, .97, "NDC");
  title.SetBorderSize(0); title.SetFillStyle(0); title.SetTextAlign(12); title.SetTextSize(.025);
  title.AddText("kbeam_momentum_comparison.cc");
  title.AddText("K18 pk18 shape from each configured ROOT file; normalization from momentum totals in kbeam_statistics.txt");
  TDatime now;
  title.AddText(Form("Executed: %04d-%02d-%02d %02d:%02d:%02d; totals use K-Beam weighted DAQ efficiency",
                     now.GetYear(), now.GetMonth(), now.GetDay(), now.GetHour(), now.GetMinute(), now.GetSecond()));
  title.Draw();
  for (std::size_t i=0; i<kSettings.size(); ++i) {
    const auto stats_it = statistics.find(kSettings[i].central_momentum_mev);
    const std::vector<Int_t> empty_runs;
    const auto& runs = stats_it == statistics.end() ? empty_runs : stats_it->second.runs;
    const Double_t k_beam_sum = stats_it == statistics.end()
      ? 0. : static_cast<Double_t>(stats_it->second.k_beam);
    const std::string input_file = kSettings[i].input_file;
    const std::string input_basename = input_file.substr(input_file.find_last_of('/') + 1);
    const Double_t column_left = .065 + i * .30;
    TPaveText details(column_left, .05, column_left + .285, .77, "NDC");
    details.SetBorderSize(0); details.SetFillStyle(0); details.SetTextAlign(12); details.SetTextSize(.016);
    details.AddText(Form("%d MeV/c", kSettings[i].central_momentum_mev));
    details.AddText(Form("Shape run %05d", kSettings[i].input_run));
    details.AddText(input_basename.c_str());
    details.AddText(Form("K-Beam sum: %.0f", k_beam_sum));
    details.AddText(Form("K-Beam #times DAQ sum: %.3e", effective_beam[i]));
    details.AddText("Scaler runs:");
    for (std::size_t j=0; j<runs.size(); j+=5) {
      std::ostringstream run_line;
      const std::size_t stop = std::min(j+5, runs.size());
      for (std::size_t k=j; k<stop; ++k) {
        if (k != j) run_line << " ";
        run_line << Form("%05d", runs[k]);
      }
      details.AddText(run_line.str().c_str());
    }
    details.Draw();
  }
  canvas.Print(pdf_name.c_str());

  auto combined = std::make_unique<TH1D>("h_pk18_685_plus_715_plus_735_scaled",
                                           "", kMomentumBins, kMomentumMin, kMomentumMax);
  combined->SetDirectory(nullptr);
  combined->Sumw2();
  for (const auto& spectrum : spectra) combined->Add(spectrum.get());
  if (!spectra.empty()) spectra.front()->SetMaximum(combined->GetMaximum() * 1.12);

  const Int_t colors[] = {kGreen + 2, kBlue + 1, kRed + 1};
  for (std::size_t i=0; i<spectra.size(); ++i) {
    spectra[i]->SetLineColor(colors[i]);
    spectra[i]->SetLineWidth(i == 0 ? 3 : 2);
    spectra[i]->SetTitle("Individual and combined K18 momentum counts;|p_{K18}| [GeV/c];K-Beam#timesDAQ weighted counts");
    spectra[i]->Draw(i == 0 ? "hist" : "hist same");
  }
  combined->SetLineColor(kBlack);
  combined->SetLineStyle(2);
  combined->SetLineWidth(3);
  combined->Draw("hist same");
  TLegend legend(.20, .68, .40, .88);
  legend.SetBorderSize(0); legend.SetFillStyle(0);
  legend.SetTextSize(.04);
  for (std::size_t i=0; i<spectra.size(); ++i)
    legend.AddEntry(spectra[i].get(), Form("%d MeV/c", kSettings[i].central_momentum_mev), "l");
  legend.AddEntry(combined.get(), "Combined: 685 + 715 + 735", "l");
  legend.Draw();
  canvas.Print(pdf_name.c_str());
  canvas.Print((pdf_name + "]").c_str());

  TFile output(root_name.c_str(), "RECREATE");
  for (auto& histogram : spectra) histogram->Write();
  combined->Write();
  canvas.Write("c_kbeam_momentum_comparison");
  output.Write();
  output.Close();
  Info("kbeam_momentum_comparison", "Wrote %s and %s", pdf_name.c_str(), root_name.c_str());
}
