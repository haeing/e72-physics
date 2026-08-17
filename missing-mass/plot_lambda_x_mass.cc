#include <array>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

namespace {
struct Reaction {
  const char *fileName;
  const char *label;
  const char *crossSectionGraphName;
  Color_t color;
};

const std::array<Reaction, 7> kReactions = {{
    {"lambda_eta_DstTPCGeant4LambdaEta.root", "#Lambda#eta", "g_Lambda_eta_pK", kRed + 1},
    {"lambda_pi_DstTPCGeant4LambdaEta.root", "#Lambda#pi^{0}", "g_Lambda_pi0_pK", kBlue + 1},
    {"lambda_pi_pi_DstTPCGeant4LambdaEta.root", "#Lambda#pi^{+}#pi^{-}", "g_Lambda_pip_pim_pK", kGreen + 2},
    {"sigma0_pi0_DstTPCGeant4LambdaEta.root", "#Sigma^{0}#pi^{0}", "g_Sigma0_pi0_pK", kMagenta + 1},
    {"sigma0_pip_pim_DstTPCGeant4LambdaEta.root", "#Sigma^{0}#pi^{+}#pi^{-}", "g_Sigma0_pip_pim_pK", kOrange + 7},
    {"sigmap_pim_DstTPCGeant4LambdaEta.root", "#Sigma^{+}#pi^{-}", "g_Sigmap_pim_pK", kCyan + 2},
    {"sigmap_pim_pi0_DstTPCGeant4LambdaEta.root", "#Sigma^{+}#pi^{-}#pi^{0}", "g_Sigmap_pim_pi0_pK", kViolet + 1},
}};

bool GetCrossSectionMb(TFile *input, const char *graphName,
                       double beamMomentumMeV, double &crossSectionMb)
{
  auto *graph = dynamic_cast<TGraph *>(input->Get(graphName));
  if (!graph || graph->GetN() == 0) {
    std::cerr << "Cross-section graph not found: " << graphName << std::endl;
    return false;
  }

  double minMomentum;
  double unusedCrossSection;
  graph->GetPoint(0, minMomentum, unusedCrossSection);
  double maxMomentum = minMomentum;
  for (int point = 1; point < graph->GetN(); ++point) {
    double momentum;
    graph->GetPoint(point, momentum, unusedCrossSection);
    minMomentum = std::min(minMomentum, momentum);
    maxMomentum = std::max(maxMomentum, momentum);
  }
  if (beamMomentumMeV < minMomentum || beamMomentumMeV > maxMomentum) {
    std::cerr << graphName << " has no data at p_{K-} = "
              << beamMomentumMeV << " MeV/c (available range: "
              << minMomentum << "-" << maxMomentum << " MeV/c)."
              << std::endl;
    return false;
  }

  crossSectionMb = graph->Eval(beamMomentumMeV);
  return true;
}

void FillMassHistogram(TTree *tree, const char *branchName, TH1D *histogram)
{
  TTreeReader reader(tree);
  TTreeReaderValue<std::vector<double>> masses(reader, branchName);

  while (reader.Next()) {
    for (const double mass : *masses)
      histogram->Fill(mass);
  }
}

bool IsSelectedTriggerFlag(Int_t triggerFlag)
{
  return triggerFlag == 18 || triggerFlag == 20 || triggerFlag == 22;
}

void FillTriggeredMassHistogram(TTree *tree, const char *branchName,
                                TH1D *histogram)
{
  TTreeReader reader(tree);
  TTreeReaderValue<Int_t> triggerFlag(reader, "trig_flag_g4");
  TTreeReaderValue<std::vector<double>> masses(reader, branchName);

  while (reader.Next()) {
    if (!IsSelectedTriggerFlag(*triggerFlag))
      continue;
    for (const double mass : *masses)
      histogram->Fill(mass);
  }
}

void FillTriggeredConstrainedXMass(TTree *tree, TH1D *histogram)
{
  TTreeReader reader(tree);
  TTreeReaderValue<Int_t> triggerFlag(reader, "trig_flag_g4");
  TTreeReaderValue<std::vector<double>> lambdaMasses(reader, "lambda_mass");
  TTreeReaderValue<std::vector<double>> xMasses(reader, "X_mass");
  TTreeReaderValue<std::vector<double>> xMass2(reader, "X_mass2");

  while (reader.Next()) {
    if (!IsSelectedTriggerFlag(*triggerFlag))
      continue;
    const std::size_t candidates = std::min({lambdaMasses->size(), xMasses->size(), xMass2->size()});
    for (std::size_t i = 0; i < candidates; ++i) {
      if (1.10 <= lambdaMasses->at(i) && lambdaMasses->at(i) <= 1.12)
        histogram->Fill(xMasses->at(i));
    }
  }
}

void FillMassCorrelation(TTree *tree, TH2D *histogram, double weight)
{
  TTreeReader reader(tree);
  TTreeReaderValue<std::vector<double>> lambdaMasses(reader, "lambda_mass");
  TTreeReaderValue<std::vector<double>> xMasses(reader, "X_mass");
  TTreeReaderValue<std::vector<double>> xMass2(reader, "X_mass2");

  while (reader.Next()) {
    const std::size_t candidates = std::min({lambdaMasses->size(), xMasses->size(), xMass2->size()});
    for (std::size_t i = 0; i < candidates; ++i)
      histogram->Fill(xMasses->at(i), lambdaMasses->at(i), weight);
  }
}

void FillDataMassHistograms(TTree *tree, TH1D *lambdaAll,
                            TH1D *lambdaCut, TH1D *xAll, TH1D *xCut,
                            TH1D *xMass2All, TH1D *xMass2Cut, TH2D *correlation)
{
  TTreeReader reader(tree);
  TTreeReaderValue<std::vector<double>> lambdaMasses(reader, "lambda_mass");
  TTreeReaderValue<std::vector<double>> xMasses(reader, "X_mass");
  TTreeReaderValue<std::vector<double>> xMass2(reader, "X_mass2");

  while (reader.Next()) {
    const std::size_t candidates = std::min({lambdaMasses->size(), xMasses->size(), xMass2->size()});
    for (std::size_t i = 0; i < candidates; ++i) {
      const double lambdaMass = lambdaMasses->at(i);
      const double xMass = xMasses->at(i);
      const double xMassSquared = xMass2->at(i);
      lambdaAll->Fill(lambdaMass);
      xAll->Fill(xMass);
      xMass2All->Fill(xMassSquared);
      correlation->Fill(xMass, lambdaMass);
      if (1.10 <= lambdaMass && lambdaMass <= 1.12) {
        lambdaCut->Fill(lambdaMass);
        xCut->Fill(xMass);
        xMass2Cut->Fill(xMassSquared);
      }
    }
  }
}

void DrawMassPage(TCanvas *canvas, const char *pdfName, TH1D *total,
                  const std::vector<std::unique_ptr<TH1D>> &byReaction,
                  const char *xTitle, const char *pageTitle)
{
  total->SetTitle(Form("%s;%s;Candidates / bin [mb / generated event]",
                       pageTitle,
                       xTitle));
  total->SetLineColor(kBlack);
  total->SetLineWidth(3);
  total->SetFillStyle(0);
  total->Draw("hist");

  auto *legend = new TLegend(0.18, 0.55, 0.46, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(total, "All reactions", "l");

  for (std::size_t i = 0; i < byReaction.size(); ++i) {
    byReaction[i]->SetLineColor(kReactions[i].color);
    byReaction[i]->SetLineWidth(2);
    byReaction[i]->SetFillStyle(0);
    byReaction[i]->Draw("hist same");
    legend->AddEntry(byReaction[i].get(), kReactions[i].label, "l");
  }

  legend->Draw();
  canvas->Print(pdfName);
  canvas->Clear();
}

} // namespace

// Run with:
//   root -l -b -q 'plot_lambda_x_mass.cc()'
//   root -l -b -q 'plot_lambda_x_mass.cc(800)'
void plot_lambda_x_mass(int beamMomentumMeV = 735,
                        const char *crossSectionFileName = "../lambda-yield/cs-info/total_cross_section.root")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const std::string inputDir =
    gSystem->ExpandPathName(Form("~/simul-data/e72-k18ana/mom-%d",beamMomentumMeV));
  const std::string outputPdf = Form("lambda-x-mass-by-reaction-mom-%d.pdf", beamMomentumMeV);
  const std::string dataPath = gSystem->ExpandPathName(Form(
      "~/data/JPARC2025Nov_root/physics-%d/run02447_DstTPCLambdaEta.root",
      beamMomentumMeV));

  const std::string crossSectionPath = gSystem->ExpandPathName(crossSectionFileName);
  std::unique_ptr<TFile> crossSectionInput(
      TFile::Open(crossSectionPath.c_str(), "READ"));
  if (!crossSectionInput || crossSectionInput->IsZombie()) {
    std::cerr << "Cannot open cross-section file: " << crossSectionPath << std::endl;
    return;
  }

  std::array<double, kReactions.size()> crossSectionsMb;
  for (std::size_t i = 0; i < kReactions.size(); ++i) {
    if (!GetCrossSectionMb(crossSectionInput.get(),
                           kReactions[i].crossSectionGraphName,
                           beamMomentumMeV, crossSectionsMb[i]))
      return;
  }
  // Masses in the input trees are in GeV/c^{2}.
  auto lambdaTotal = std::make_unique<TH1D>(
      "h_lambda_mass_total", ";M_{p#pi^{-}} [GeV/#it{c}^{2}];Candidates / bin",
      110, 1.05, 1.16);
  auto xMass2Total = std::make_unique<TH1D>(
      "h_X_mass2_total", ";M_{X}^{2} [GeV^{2}/#it{c}^{4}];Candidates / bin",
      250, 0.0, 0.5);
  auto xTotal = std::make_unique<TH1D>(
      "h_X_mass_total", ";M_{X} [GeV/#it{c}^{2}];Candidates / bin",
      234, 0.0, 0.702);
  lambdaTotal->SetDirectory(nullptr);
  xTotal->SetDirectory(nullptr);
  xMass2Total->SetDirectory(nullptr);
  auto xLambdaTotal = std::make_unique<TH2D>(
      "h_X_vs_lambda_mass_total",
      ";M_{X} [GeV/#it{c}^{2}];M_{p#pi^{-}} [GeV/#it{c}^{2}];",
      234, 0.0, 0.702, 110, 1.05, 1.16);
  xLambdaTotal->SetDirectory(nullptr);

  std::vector<std::unique_ptr<TH1D>> lambdaByReaction;
  std::vector<std::unique_ptr<TH1D>> xByReaction;
  std::vector<std::unique_ptr<TH1D>> xMass2ByReaction;
  auto lambdaTriggeredTotal = std::make_unique<TH1D>(
      "h_lambda_mass_triggered_total", "", 110, 1.05, 1.16);
  auto xMass2TriggeredTotal = std::make_unique<TH1D>(
      "h_X_mass2_triggered_total", "", 250, 0.0, 0.5);
  auto xTriggeredTotal = std::make_unique<TH1D>(
      "h_X_mass_triggered_total", "", 234, 0.0, 0.702);
  lambdaTriggeredTotal->SetDirectory(nullptr);
  xTriggeredTotal->SetDirectory(nullptr);
  xMass2TriggeredTotal->SetDirectory(nullptr);
  auto etaSignalTemplate = std::make_unique<TH1D>(
      "h_eta_signal_template", "", 234, 0.0, 0.702);
  auto etaBackgroundTemplate = std::make_unique<TH1D>(
      "h_eta_background_template", "", 234, 0.0, 0.702);
  etaSignalTemplate->SetDirectory(nullptr);
  etaBackgroundTemplate->SetDirectory(nullptr);
  std::vector<std::unique_ptr<TH1D>> lambdaTriggeredByReaction;
  std::vector<std::unique_ptr<TH1D>> xTriggeredByReaction;
  std::vector<std::unique_ptr<TH1D>> xMass2TriggeredByReaction;

  for (std::size_t i = 0; i < kReactions.size(); ++i) {
    const std::string path = inputDir + "/" + kReactions[i].fileName;
    std::unique_ptr<TFile> input(TFile::Open(path.c_str(), "READ"));
    if (!input || input->IsZombie()) {
      std::cerr << "Cannot open " << path << std::endl;
      continue;
    }

    auto *tree = dynamic_cast<TTree *>(input->Get("tpc"));
    if (!tree) {
      std::cerr << "tpc tree not found in " << path << std::endl;
      continue;
    }

    const double crossSectionMb = crossSectionsMb[i];

    auto lambdaHist = std::make_unique<TH1D>(
        Form("h_lambda_mass_%zu", i), "", 110, 1.05, 1.16);
    auto xMass2Hist = std::make_unique<TH1D>(
        Form("h_X_mass2_%zu", i), "", 250, 0.0, 0.5);
    auto xHist = std::make_unique<TH1D>(
        Form("h_X_mass_%zu", i), "", 234, 0.0, 0.702);
    lambdaHist->SetDirectory(nullptr);
    xHist->SetDirectory(nullptr);
    xMass2Hist->SetDirectory(nullptr);
    auto lambdaTriggeredHist = std::make_unique<TH1D>(
        Form("h_lambda_mass_triggered_%zu", i), "", 110, 1.05, 1.16);
    auto xMass2TriggeredHist = std::make_unique<TH1D>(
        Form("h_X_mass2_triggered_%zu", i), "", 250, 0.0, 0.5);
    auto xTriggeredHist = std::make_unique<TH1D>(
        Form("h_X_mass_triggered_%zu", i), "", 234, 0.0, 0.702);
    lambdaTriggeredHist->SetDirectory(nullptr);
    xTriggeredHist->SetDirectory(nullptr);
    xMass2TriggeredHist->SetDirectory(nullptr);
    auto templateHist = std::make_unique<TH1D>(
        Form("h_eta_template_%zu", i), "", 234, 0.0, 0.702);
    templateHist->SetDirectory(nullptr);

    FillMassHistogram(tree, "lambda_mass", lambdaHist.get());
    FillMassHistogram(tree, "X_mass", xHist.get());
    FillMassHistogram(tree, "X_mass2", xMass2Hist.get());
    FillTriggeredMassHistogram(tree, "lambda_mass", lambdaTriggeredHist.get());
    FillTriggeredMassHistogram(tree, "X_mass", xTriggeredHist.get());
    FillTriggeredMassHistogram(tree, "X_mass2", xMass2TriggeredHist.get());
    FillTriggeredConstrainedXMass(tree, templateHist.get());
    const double weight = crossSectionMb / tree->GetEntries();
    lambdaHist->Scale(weight);
    xHist->Scale(weight);
    xMass2Hist->Scale(weight);
    lambdaTriggeredHist->Scale(weight);
    xTriggeredHist->Scale(weight);
    xMass2TriggeredHist->Scale(weight);
    templateHist->Scale(weight);
    FillMassCorrelation(tree, xLambdaTotal.get(), weight);
    std::cout << kReactions[i].label << ": #sigma("
              << beamMomentumMeV << " MeV/c) = "
              << crossSectionMb << " mb, weight = "
              << weight << " mb/event" << std::endl;
    lambdaTotal->Add(lambdaHist.get());
    xTotal->Add(xHist.get());
    xMass2Total->Add(xMass2Hist.get());
    lambdaTriggeredTotal->Add(lambdaTriggeredHist.get());
    xTriggeredTotal->Add(xTriggeredHist.get());
    xMass2TriggeredTotal->Add(xMass2TriggeredHist.get());
    if (i == 0)
      etaSignalTemplate->Add(templateHist.get());
    else
      etaBackgroundTemplate->Add(templateHist.get());
    lambdaByReaction.push_back(std::move(lambdaHist));
    xByReaction.push_back(std::move(xHist));
    xMass2ByReaction.push_back(std::move(xMass2Hist));
    lambdaTriggeredByReaction.push_back(std::move(lambdaTriggeredHist));
    xTriggeredByReaction.push_back(std::move(xTriggeredHist));
    xMass2TriggeredByReaction.push_back(std::move(xMass2TriggeredHist));
  }

  auto etaTemplateTotal = std::unique_ptr<TH1D>(
      static_cast<TH1D *>(etaSignalTemplate->Clone("h_eta_template_total")));
  etaTemplateTotal->SetDirectory(nullptr);
  etaTemplateTotal->Add(etaBackgroundTemplate.get());
  const int etaWindowBinMin = etaSignalTemplate->FindBin(0.53);
  const int etaWindowBinMax = etaSignalTemplate->FindBin(0.57);
  const double etaSignalWindow = etaSignalTemplate->Integral(etaWindowBinMin, etaWindowBinMax);
  const double etaBackgroundWindow = etaBackgroundTemplate->Integral(etaWindowBinMin, etaWindowBinMax);
  const double etaWindowTotal = etaSignalWindow + etaBackgroundWindow;
  const double etaSignalPercent = etaWindowTotal > 0.0
      ? 100.0 * etaSignalWindow / etaWindowTotal : 0.0;
  const double etaBackgroundPercent = etaWindowTotal > 0.0
      ? 100.0 * etaBackgroundWindow / etaWindowTotal : 0.0;
  auto lambdaData = std::make_unique<TH1D>(
      "h_lambda_mass_data",
      ";M_{p#pi^{-}} [GeV/#it{c}^{2}];Candidates / bin",
      110, 1.05, 1.16);
  auto lambdaDataCut = std::make_unique<TH1D>(
      "h_lambda_mass_data_lambda_constraint",
      ";M_{p#pi^{-}} [GeV/#it{c}^{2}];Counts",
      110, 1.05, 1.16);
  auto xData = std::make_unique<TH1D>(
      "h_X_mass_data_lambda_constraint",
      ";M_{X} [GeV/#it{c}^{2}];Candidates / bin",
      234, 0.0, 0.702);
  auto xMass2DataAll = std::make_unique<TH1D>(
      "h_X_mass2_data", ";M_{X}^{2} [GeV^{2}/#it{c}^{4}];Counts / bin", 250, 0.0, 0.5);
  auto xMass2Data = std::make_unique<TH1D>(
      "h_X_mass2_data_lambda_constraint", ";M_{X}^{2} [GeV^{2}/#it{c}^{4}];Candidates / bin", 250, 0.0, 0.5);
  auto xDataAll = std::make_unique<TH1D>(
      "h_X_mass_data",
      ";M_{X} [GeV/#it{c}^{2}];Counts",
      234, 0.0, 0.702);
  lambdaData->SetDirectory(nullptr);
  lambdaDataCut->SetDirectory(nullptr);
  xDataAll->SetDirectory(nullptr);
  xData->SetDirectory(nullptr);
  xMass2DataAll->SetDirectory(nullptr);
  xMass2Data->SetDirectory(nullptr);

  auto xLambdaData = std::make_unique<TH2D>(
      "h_X_vs_lambda_mass_data",
      ";M_{X} [GeV/#it{c}^{2}];M_{p#pi^{-}} [GeV/#it{c}^{2}];Candidates / bin",
      234, 0.0, 0.702, 110, 1.05, 1.16);
  xLambdaData->SetDirectory(nullptr);
  std::unique_ptr<TFile> dataFile(TFile::Open(dataPath.c_str(), "READ"));
  if (!dataFile || dataFile->IsZombie()) {
    std::cerr << "Cannot open data file " << dataPath << std::endl;
  } else if (auto *dataTree = dynamic_cast<TTree *>(dataFile->Get("tpc"))) {
    FillDataMassHistograms(dataTree, lambdaData.get(), lambdaDataCut.get(),
                           xDataAll.get(), xData.get(), xMass2DataAll.get(), xMass2Data.get(),
                           xLambdaData.get());
  } else {
    std::cerr << "tpc tree not found in " << dataPath << std::endl;
  }
  auto dataSignalEstimate = std::unique_ptr<TH1D>(
      static_cast<TH1D *>(xData->Clone("h_data_signal_mc_fraction")));
  auto dataBackgroundEstimate = std::unique_ptr<TH1D>(
      static_cast<TH1D *>(xData->Clone("h_data_background_mc_fraction")));
  dataSignalEstimate->SetDirectory(nullptr);
  dataBackgroundEstimate->SetDirectory(nullptr);
  dataSignalEstimate->Reset("ICES");
  dataBackgroundEstimate->Reset("ICES");
  for (int bin = 1; bin <= xData->GetNbinsX(); ++bin) {
    const double templateTotal = etaTemplateTotal->GetBinContent(bin);
    const double signalFraction = templateTotal > 0.0
        ? etaSignalTemplate->GetBinContent(bin) / templateTotal : 0.0;
    const double dataCount = xData->GetBinContent(bin);
    dataSignalEstimate->SetBinContent(bin, dataCount * signalFraction);
    dataBackgroundEstimate->SetBinContent(bin, dataCount * (1.0 - signalFraction));
  }
  const int dataEtaBinMin = xData->FindBin(0.53);
  const int dataEtaBinMax = xData->FindBin(0.57);
  const double dataSignalEstimateWindow =
      dataSignalEstimate->Integral(dataEtaBinMin, dataEtaBinMax);
  const double dataBackgroundEstimateWindow =
      dataBackgroundEstimate->Integral(dataEtaBinMin, dataEtaBinMax);
  const double etaBackgroundFitMin = 0.00;
  const double etaBackgroundFitMax = 0.55;
  const double etaFitMin = 0.53;
  const double etaFitMax = 0.57;
  auto etaBackgroundFit = std::make_unique<TF1>(
      "f_eta_background_fit", "pol2", etaBackgroundFitMin, etaBackgroundFitMax);
  const int etaBackgroundFitStatus =
      xData->Fit(etaBackgroundFit.get(), "Q0R");
  auto etaBackground = std::make_unique<TF1>(
      "f_eta_background", "pol2", 0.0, 0.702);
  etaBackground->SetParameters(etaBackgroundFit->GetParameter(0),
                              etaBackgroundFit->GetParameter(1),
                              etaBackgroundFit->GetParameter(2));
  auto etaFit = std::make_unique<TF1>(
      "f_eta_total", "gaus(0)+pol2(3)", etaFitMin, etaFitMax);
  const double etaPeakHeight = xData->GetBinContent(xData->FindBin(0.548));
  etaFit->SetParameters(std::max(1.0, etaPeakHeight - etaBackground->Eval(0.548)),
                        0.548, 0.008, etaBackground->GetParameter(0),
                        etaBackground->GetParameter(1), etaBackground->GetParameter(2));
  etaFit->SetParLimits(0, 0.0, 1.0e9);
  etaFit->SetParLimits(1, 0.530, 0.565);
  etaFit->SetParLimits(2, 0.002, 0.030);
  etaFit->FixParameter(3, etaBackground->GetParameter(0));
  etaFit->FixParameter(4, etaBackground->GetParameter(1));
  etaFit->FixParameter(5, etaBackground->GetParameter(2));
  const int etaFitStatus = xData->Fit(etaFit.get(), "Q0R");
  auto etaSignal = std::make_unique<TF1>(
      "f_eta_signal", "gaus", etaFitMin, etaFitMax);
  etaSignal->SetParameters(etaFit->GetParameter(0), etaFit->GetParameter(1),
                           etaFit->GetParameter(2));
  const double etaYield = etaFit->GetParameter(0) * etaFit->GetParameter(2) *
                          std::sqrt(2.0 * std::acos(-1.0)) / xData->GetBinWidth(1);
  std::cout << "eta background fit status = " << etaBackgroundFitStatus
            << ", eta signal fit status = " << etaFitStatus
            << ", raw eta yield = " << etaYield << std::endl;

  auto *canvas = new TCanvas("c_mass", "Mass by reaction", 1000, 750);
  canvas->SetTicks(1, 1);
  canvas->Print((std::string(outputPdf) + "[").c_str());
  DrawMassPage(canvas, outputPdf.c_str(), lambdaTotal.get(), lambdaByReaction,
               "M_{p#pi^{-}} [GeV/#it{c}^{2}]", "");
  DrawMassPage(canvas, outputPdf.c_str(), xTotal.get(), xByReaction,
               "M_{X} [GeV/#it{c}^{2}]", "");
  DrawMassPage(canvas, outputPdf.c_str(), xMass2Total.get(), xMass2ByReaction,
               "M_{X}^{2} [GeV^{2}/#it{c}^{4}]", "");
  DrawMassPage(canvas, outputPdf.c_str(), lambdaTriggeredTotal.get(),
               lambdaTriggeredByReaction, "M_{p#pi^{-}} [GeV/#it{c}^{2}]",
               "Simulation: trig_flag_g4 = 18, 20, 22");
  DrawMassPage(canvas, outputPdf.c_str(), xTriggeredTotal.get(), xTriggeredByReaction,
               "M_{X} [GeV/#it{c}^{2}]", "Simulation: trig_flag_g4 = 18, 20, 22");
  DrawMassPage(canvas, outputPdf.c_str(), xMass2TriggeredTotal.get(), xMass2TriggeredByReaction,
               "M_{X}^{2} [GeV^{2}/#it{c}^{4}]", "Simulation: trig_flag_g4 = 18, 20, 22");
  etaTemplateTotal->SetTitle("Simulation signal/background template;M_{X} [GeV/#it{c}^{2}];Counts / (3 MeV/#it{c}^{2}) [mb / generated event]");
  etaTemplateTotal->SetLineColor(kBlack);
  etaTemplateTotal->SetLineWidth(3);
  etaTemplateTotal->Draw("hist");
  etaBackgroundTemplate->SetLineColor(kBlue + 1);
  etaBackgroundTemplate->SetLineWidth(3);
  etaBackgroundTemplate->Draw("hist same");
  etaSignalTemplate->SetLineColor(kRed + 1);
  etaSignalTemplate->SetLineWidth(3);
  etaSignalTemplate->Draw("hist same");
  auto *templateLegend = new TLegend(0.18, 0.63, 0.53, 0.88);
  templateLegend->SetBorderSize(0);
  templateLegend->SetFillStyle(0);
  templateLegend->AddEntry(etaTemplateTotal.get(), "Signal + background", "l");
  templateLegend->AddEntry(etaSignalTemplate.get(), "#Lambda#eta signal", "l");
  templateLegend->AddEntry(etaBackgroundTemplate.get(), "Other reactions (background)", "l");
  templateLegend->Draw();
  auto *templateText = new TPaveText(0.18, 0.48, 0.53, 0.60, "NDC");
  templateText->SetBorderSize(0);
  templateText->SetFillStyle(0);
  templateText->AddText("530 < M_{X} < 570 MeV/#it{c}^{2}");
  templateText->AddText(Form("Signal: %.1f%%, Background: %.1f%%",
                             etaSignalPercent, etaBackgroundPercent));
  templateText->Draw();
  canvas->Print(outputPdf.c_str());
  canvas->Clear();
  xLambdaTotal->SetTitle("Simulation (cross-section weighted);M_{X} [GeV/#it{c}^{2}];M_{p#pi^{-}} [GeV/#it{c}^{2}]");
  xLambdaTotal->Draw("colz");
  canvas->Print(outputPdf.c_str());
  canvas->Clear();
  lambdaData->SetTitle("Run 02447 data: #Lambda mass;M_{p#pi^{-}} [GeV/#it{c}^{2}];Counts / (1 MeV/#it{c}^{2})");
  lambdaData->SetLineColor(kBlack);
  lambdaData->SetLineWidth(3);
  lambdaData->Draw("hist");
  lambdaDataCut->SetLineColor(kRed + 1);
  lambdaDataCut->SetLineWidth(3);
  lambdaDataCut->Draw("hist same");
  auto *lambdaLegend = new TLegend(0.18, 0.53, 0.50, 0.88);
  lambdaLegend->SetBorderSize(0);
  lambdaLegend->SetFillStyle(0);
  lambdaLegend->AddEntry(lambdaData.get(), "All candidates", "l");
  lambdaLegend->AddEntry(lambdaDataCut.get(), "1.10 #leq M_{p#pi^{-}} #leq 1.12", "l");
  lambdaLegend->Draw();
  canvas->Print(outputPdf.c_str());
  canvas->Clear();
  xMass2DataAll->SetTitle("Run 02447 data: X mass squared;M_{X}^{2} [GeV^{2}/#it{c}^{4}];Candidates / bin");
  xMass2DataAll->SetLineColor(kBlack);
  xMass2DataAll->SetLineWidth(3);
  xMass2DataAll->Draw("hist");
  xMass2Data->SetLineColor(kRed + 1);
  xMass2Data->SetLineWidth(3);
  xMass2Data->Draw("hist same");
  auto *xMass2Legend = new TLegend(0.18, 0.63, 0.50, 0.88);
  xMass2Legend->SetBorderSize(0);
  xMass2Legend->SetFillStyle(0);
  xMass2Legend->AddEntry(xMass2DataAll.get(), "All candidates", "l");
  xMass2Legend->AddEntry(xMass2Data.get(), "1.10 #leq M_{p#pi^{-}} #leq 1.12", "l");
  xMass2Legend->Draw();
  canvas->Print(outputPdf.c_str());
  canvas->Clear();
  xDataAll->SetTitle("Run 02447 data: X mass;M_{X} [GeV/#it{c}^{2}];Counts / (3 MeV/#it{c}^{2})");
  xDataAll->SetLineColor(kBlack);
  xDataAll->SetLineWidth(3);
  xDataAll->Draw("hist");
  xData->SetLineColor(kRed + 1);
  xData->SetLineWidth(3);
  xData->Draw("hist same");
  etaBackground->SetLineColor(kBlue + 1);
  etaBackground->SetLineStyle(2);
  etaBackground->SetLineWidth(3);
  etaBackground->Draw("same");
  etaSignal->SetLineColor(kGreen + 2);
  etaSignal->SetLineStyle(2);
  etaSignal->SetLineWidth(3);
  etaSignal->Draw("same");
  etaFit->SetLineColor(kMagenta + 1);
  etaFit->SetLineWidth(3);
  etaFit->Draw("same");
  auto *etaYieldText = new TPaveText(0.18, 0.34, 0.50, 0.50, "NDC");
  etaYieldText->SetBorderSize(0);
  etaYieldText->SetFillStyle(0);
  etaYieldText->AddText(Form("#eta fit: %.0f < M_{X} < %.0f MeV/#it{c}^{2}",
                             1000.0 * etaFitMin, 1000.0 * etaFitMax));
  etaYieldText->AddText(Form("Raw #eta yield = %.0f", etaYield));
  etaYieldText->AddText("Background fit: 0 < M_{X} < 550 MeV/#it{c}^{2}");
  etaYieldText->Draw();
  auto *xLegend = new TLegend(0.18, 0.53, 0.50, 0.88);
  xLegend->SetBorderSize(0);
  xLegend->SetFillStyle(0);
  xLegend->AddEntry(xDataAll.get(), "All candidates", "l");
  xLegend->AddEntry(xData.get(), "1.10 #leq M_{p#pi^{-}} #leq 1.12", "l");
  xLegend->AddEntry(etaBackground.get(), "Background (pol2)", "l");
  xLegend->AddEntry(etaSignal.get(), "#eta signal (Gaussian)", "l");
  xLegend->AddEntry(etaFit.get(), "Signal + background fit", "l");
  xLegend->Draw();
  canvas->Print(outputPdf.c_str());
  canvas->Clear();
  xData->SetTitle("Run 02447 data yield split by MC S/(S+B);M_{X} [GeV/#it{c}^{2}];Counts / (3 MeV/#it{c}^{2})");
  xData->SetLineColor(kBlack);
  xData->SetLineWidth(3);
  xData->Draw("hist");
  dataBackgroundEstimate->SetLineColor(kBlue + 1);
  dataBackgroundEstimate->SetLineWidth(3);
  dataBackgroundEstimate->Draw("hist same");
  dataSignalEstimate->SetLineColor(kRed + 1);
  dataSignalEstimate->SetLineWidth(3);
  dataSignalEstimate->Draw("hist same");
  auto *dataTemplateLegend = new TLegend(0.18, 0.63, 0.53, 0.88);
  dataTemplateLegend->SetBorderSize(0);
  dataTemplateLegend->SetFillStyle(0);
  dataTemplateLegend->AddEntry(xData.get(), "Lambda-cut data", "l");
  dataTemplateLegend->AddEntry(dataSignalEstimate.get(), "MC-fraction signal estimate", "l");
  dataTemplateLegend->AddEntry(dataBackgroundEstimate.get(), "MC-fraction background estimate", "l");
  dataTemplateLegend->Draw();
  auto *dataTemplateText = new TPaveText(0.18, 0.48, 0.53, 0.60, "NDC");
  dataTemplateText->SetBorderSize(0);
  dataTemplateText->SetFillStyle(0);
  dataTemplateText->AddText("530 < M_{X} < 570 MeV/#it{c}^{2}");
  dataTemplateText->AddText(Form("Estimated S = %.0f, B = %.0f",
                                 dataSignalEstimateWindow, dataBackgroundEstimateWindow));
  dataTemplateText->Draw();
  canvas->Print(outputPdf.c_str());
  canvas->Clear();
  xLambdaData->SetTitle("Run 02447 data;M_{X} [GeV/#it{c}^{2}];M_{p#pi^{-}} [GeV/#it{c}^{2}]");
  xLambdaData->Draw("colz");
  canvas->Print(outputPdf.c_str());
  canvas->Clear();
  canvas->Print((std::string(outputPdf) + "]").c_str());
}
