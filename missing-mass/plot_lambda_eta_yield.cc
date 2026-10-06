// Run with: root -l -b -q plot_lambda_eta_yield.cc
#include <algorithm>
#include <cmath>
#include <memory>
#include <map>
#include <vector>

#include <TCanvas.h>
#include <TDatime.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TParameter.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

namespace {
// Edit these settings before running.
constexpr Int_t kBeamMomentum = 715; // MeV/c; selects physics-<momentum> input directory
  //735
  /*
const std::vector<Int_t> kRunNumbers = {
  2447, 2449, 2450, 2451, 2452, 2453, 2454, 2456,
  2457, 2458, 2459, 2460, 2462, 2463, 2465, 2468
};
  */

  //715
  const std::vector<Int_t> kRunNumbers = {
    2682, 2683, 2684, 2686, 2687, 2689, 2690, 2691, 2692, 2693
};

const TString kDataDirectory = Form("/gpfs/home/had/haein/data/JPARC2025Nov_root/physics-%d", kBeamMomentum);
const TString kOutputDirectory = "result";

// LambdaEta selection and eta-yield extraction (GeV/c^2).
constexpr Bool_t kRequireProductionVertexInsideLH2 = true;
constexpr Double_t kTargetZ = -143.;
constexpr Double_t kLH2Radius = 38.;
constexpr Double_t kLH2HalfLengthY = 48.;
constexpr Double_t kLambdaMassMin = 1.09;
constexpr Double_t kLambdaMassMax = 1.13;
constexpr Double_t kLambdaMassDisplayMin = 1.04;
constexpr Double_t kLambdaMassDisplayMax = 1.20;
constexpr Double_t kLambdaPeakFitMin = 1.10;
constexpr Double_t kLambdaPeakFitMax = 1.13;
constexpr Double_t kLambdaPeakMassInitial = 1.115683;
constexpr Double_t kLambdaPeakSigmaInitial = 0.004;
constexpr Int_t kMissingMassBins = 300; // 2 MeV/c^2 bins
constexpr Double_t kMissingMassMin = 0.;
constexpr Double_t kMissingMassMax = 0.60;
constexpr Double_t kEtaMassInitial = 0.547862;
constexpr Double_t kEtaPeakFitMin = 0.50;
constexpr Double_t kEtaPeakFitMax = 0.60;
constexpr Double_t kEtaPeakMeanMin = 0.54;
constexpr Double_t kEtaPeakMeanMax = 0.56;
constexpr Double_t kEtaPeakSigmaMin = 0.001;
constexpr Double_t kEtaPeakSigmaMax = 0.005;
constexpr Double_t kEtaBackgroundFitMin = 0.0;
constexpr Double_t kEtaBackgroundFitMax = 0.53;
constexpr Double_t kKaonMass = 0.493677;
constexpr Double_t kProtonMass = 0.938272;
constexpr Double_t kLambdaMassNominal = 1.115683;
// Production-vertex beam momentum binning from plot_lambda_missing_mass_vertex.cc.
constexpr Double_t kBeamMomentumBinWidth = 0.004; // GeV/c (4 MeV/c)
constexpr Double_t kBeamMomentumMin = 0.700;
constexpr Int_t kBeamMomentumBins = 25; // 0.700--0.800 GeV/c
constexpr Double_t kBeamMomentumMax = kBeamMomentumMin + kBeamMomentumBins*kBeamMomentumBinWidth;
constexpr Double_t kMinEntriesPerBeamMassPanel = 20.;

bool IsInsideLH2(Double_t x, Double_t y, Double_t z)
{
  return std::hypot(x, z-kTargetZ) < kLH2Radius && std::abs(y) < kLH2HalfLengthY;
}

TString RunRangeTag()
{
  if (kRunNumbers.empty()) return "noruns";
  return Form("run%05d-run%05d", kRunNumbers.front(), kRunNumbers.back());
}

Double_t LambdaPeakWithLinearBackground(Double_t* x, Double_t* p)
{
  const Double_t pull = (x[0]-p[1])/p[2];
  return p[0]*std::exp(-0.5*pull*pull) + p[3] + p[4]*x[0];
}

Double_t EtaLinearBackground(Double_t* x, Double_t* p)
{
  return p[0] + p[1]*x[0];
}

// Linear continuum truncated at the kinematic endpoint mX_max = sqrt(s) - mLambda.
Double_t EtaGaussianPlusLinearBackground(Double_t* x, Double_t* p)
{
  const Double_t endpoint = p[5];
  const Double_t continuum = (x[0] >= 0. && x[0] < endpoint && endpoint > 0.)
    ? p[3] + p[4]*x[0] : 0.;
  const Double_t pull = (x[0]-p[1])/p[2];
  const Double_t signal = p[0]*std::exp(-0.5*pull*pull);
  return signal + continuum;
}

Double_t EtaGaussianComponent(Double_t* x, Double_t* p)
{
  const Double_t pull = (x[0]-p[1])/p[2];
  return p[0]*std::exp(-0.5*pull*pull);
}

Double_t EtaLinearBackgroundComponent(Double_t* x, Double_t* p)
{
  const Double_t endpoint = p[2];
  if (x[0] < 0. || x[0] >= endpoint || endpoint <= 0.) return 0.;
  return p[0] + p[1]*x[0];
}

Double_t KinematicMissingMassEndpoint(Double_t beam_momentum)
{
  if (!std::isfinite(beam_momentum) || beam_momentum <= 0.) return 0.;
  const Double_t beam_energy = std::sqrt(beam_momentum*beam_momentum + kKaonMass*kKaonMass);
  const Double_t s = kKaonMass*kKaonMass + kProtonMass*kProtonMass +
    2.*kProtonMass*beam_energy;
  return std::sqrt(s) - kLambdaMassNominal;
}

struct LambdaPeakFitResult {
  Double_t mass = kLambdaPeakMassInitial;
  Double_t mass_error = 0.;
  Double_t sigma = kLambdaPeakSigmaInitial;
  Double_t sigma_error = 0.;
  Int_t fit_status = -1;
};

LambdaPeakFitResult FitLambdaPeak(TH1D& histogram, TF1& fit_function)
{
  LambdaPeakFitResult result;
  fit_function.SetParameters(std::max(1., histogram.GetMaximum()), kLambdaPeakMassInitial,
                             kLambdaPeakSigmaInitial, 1., 0.);
  fit_function.SetParLimits(1, 1.10, 1.13);
  fit_function.SetParLimits(2, 0.0005, 0.02);
  TFitResultPtr fit = histogram.Fit(&fit_function, "RQS0");
  result.fit_status = static_cast<Int_t>(fit);
  result.mass = fit_function.GetParameter(1);
  result.mass_error = fit_function.GetParError(1);
  result.sigma = fit_function.GetParameter(2);
  result.sigma_error = fit_function.GetParError(2);
  return result;
}

struct EtaPeakFitResult {
  Double_t mass = kEtaMassInitial;
  Double_t mass_error = 0.;
  Double_t sigma = 0.;
  Double_t sigma_error = 0.;
  Double_t yield = 0.;
  Double_t yield_error = 0.;
  Int_t fit_status = -1;
};

EtaPeakFitResult FitEtaPeak(TH1D& mass, TF1& fit_function, Double_t endpoint)
{
  EtaPeakFitResult result;
  result.fit_status = -2; // too few entries
  if (mass.GetEntries() < kMinEntriesPerBeamMassPanel) return result;
  if (!std::isfinite(endpoint) || endpoint <= kEtaMassInitial) {
    result.fit_status = -3; // kinematically below the fitted eta-peak region
    return result;
  }

  Double_t amplitude = 1.;
  for (Int_t bin=mass.GetXaxis()->FindFixBin(kEtaPeakMeanMin);
       bin<=mass.GetXaxis()->FindFixBin(kEtaPeakMeanMax-1.e-9); ++bin) {
    const Double_t x = mass.GetBinCenter(bin);
    amplitude = std::max(amplitude, mass.GetBinContent(bin));
  }

  // Determine the linear continuum using only the non-eta sideband, then
  // freeze it while fitting the eta peak so peak bins cannot reshape it.
  TF1 sidebandFit(Form("eta_background_sideband_%s", mass.GetName()),
                  EtaLinearBackground, kEtaBackgroundFitMin, kEtaBackgroundFitMax, 2);
  sidebandFit.SetParameters(std::max(1., mass.GetMaximum()), 0.);
  mass.Fit(&sidebandFit, "RQS0N");
  const Double_t backgroundIntercept = sidebandFit.GetParameter(0);
  const Double_t backgroundSlope = sidebandFit.GetParameter(1);

  fit_function.SetParameters(amplitude, kEtaMassInitial, 0.01,
                             backgroundIntercept, backgroundSlope, endpoint);
  fit_function.SetParLimits(0, 0., std::max(10., 10.*mass.GetMaximum()));
  fit_function.SetParLimits(1, kEtaPeakMeanMin, kEtaPeakMeanMax);
  fit_function.SetParLimits(2, kEtaPeakSigmaMin, kEtaPeakSigmaMax);
  fit_function.SetParLimits(3, 0., std::max(10., 10.*mass.GetMaximum()));
  fit_function.SetParLimits(4, -20.*std::max(1., mass.GetMaximum()),
                            20.*std::max(1., mass.GetMaximum()));
  fit_function.FixParameter(3, backgroundIntercept);
  fit_function.FixParameter(4, backgroundSlope);
  fit_function.FixParameter(5, endpoint);
  fit_function.SetRange(kEtaPeakFitMin, std::min(kEtaPeakFitMax, endpoint));
  TFitResultPtr fit = mass.Fit(&fit_function, "RQS0");
  result.fit_status = static_cast<Int_t>(fit);
  result.mass = fit_function.GetParameter(1);
  result.mass_error = fit_function.GetParError(1);
  result.sigma = fit_function.GetParameter(2);
  result.sigma_error = fit_function.GetParError(2);

  const Double_t bin_width = mass.GetXaxis()->GetBinWidth(1);
  const Double_t area_factor = std::sqrt(2.*std::acos(-1.))/bin_width;
  const Double_t height = fit_function.GetParameter(0);
  result.yield = area_factor*height*result.sigma;
  if (fit.Get() && fit->CovMatrixStatus() > 0) {
    const Double_t var_height = fit->CovMatrix(0, 0);
    const Double_t var_sigma = fit->CovMatrix(2, 2);
    const Double_t cov_height_sigma = fit->CovMatrix(0, 2);
    const Double_t variance = area_factor*area_factor*
      (result.sigma*result.sigma*var_height + height*height*var_sigma +
       2.*height*result.sigma*cov_height_sigma);
    result.yield_error = std::sqrt(std::max(0., variance));
  } else {
    result.yield_error = area_factor*result.sigma*fit_function.GetParError(0);
  }
  fit_function.SetRange(kMissingMassMin, std::min(kEtaPeakFitMax, endpoint));
  return result;
}

struct EffectiveNtTpcHistograms {
  Int_t effectiveNtTpc;
  TH1D missingMass;
  TH2D missingMassVsProductionBeamMomentum;
  TH1D lambdaMassBeforeCut;
  std::vector<std::unique_ptr<TH1D>> missingMassByBeamMomentum;
  Double_t selectedBeamMomentumSum = 0.;
  Long64_t selectedBeamMomentumCount = 0;
  std::vector<Double_t> selectedBeamMomentumSumsByBin;
  std::vector<Long64_t> selectedBeamMomentumCountsByBin;
  Long64_t eventsRead = 0;
  Long64_t candidatesRead = 0;
  Long64_t lambdaCandidatesBeforeCut = 0;
  Long64_t candidatesSelected = 0;

  explicit EffectiveNtTpcHistograms(Int_t n)
    : effectiveNtTpc(n),
      missingMass(Form("h_lambda_eta_missing_mass_effnttpc_%d", n), "",
                  kMissingMassBins, kMissingMassMin, kMissingMassMax),
      missingMassVsProductionBeamMomentum(
        Form("h_missing_mass_vs_prod_k_momentum_effnttpc_%d", n), "",
        kMissingMassBins, kMissingMassMin, kMissingMassMax,
        100, kBeamMomentumMin, kBeamMomentumMax),
      lambdaMassBeforeCut(Form("h_lambda_mass_before_cut_effnttpc_%d", n), "", 320,
                          kLambdaMassDisplayMin, kLambdaMassDisplayMax),
      selectedBeamMomentumSumsByBin(kBeamMomentumBins, 0.),
      selectedBeamMomentumCountsByBin(kBeamMomentumBins, 0)
  {
    missingMass.SetDirectory(nullptr);
    missingMass.Sumw2();
    missingMassVsProductionBeamMomentum.SetDirectory(nullptr);
    missingMassVsProductionBeamMomentum.Sumw2();
    lambdaMassBeforeCut.SetDirectory(nullptr);
    lambdaMassBeforeCut.Sumw2();
    missingMassByBeamMomentum.reserve(kBeamMomentumBins);
    for (Int_t ibin=0; ibin<kBeamMomentumBins; ++ibin) {
      auto histogram = std::make_unique<TH1D>(
        Form("h_missing_mass_effnttpc_%d_beam_bin_%02d", n, ibin), "",
        kMissingMassBins, kMissingMassMin, kMissingMassMax);
      histogram->SetDirectory(nullptr);
      histogram->Sumw2();
      missingMassByBeamMomentum.push_back(std::move(histogram));
    }
  }
};

void FillRunOnce(Int_t run,
                 std::map<Int_t, std::unique_ptr<EffectiveNtTpcHistograms>>& categories)
{
  const TString path = Form("%s/run%05d_DstTPCLambdaEta.root", kDataDirectory.Data(), run);
  std::unique_ptr<TFile> input(TFile::Open(path, "READ"));
  if (!input || input->IsZombie()) {
    Warning("plot_lambda_eta_yield", "Cannot open %s", path.Data());
    return;
  }
  auto* tree = dynamic_cast<TTree*>(input->Get("tpc"));
  if (!tree) {
    Warning("plot_lambda_eta_yield", "Tree 'tpc' is absent in %s", path.Data());
    return;
  }
  const char* required[] = {"X_mass", "X_prod_found", "X_prod_vtx_x", "X_prod_vtx_y", "X_prod_vtx_z",
                            "X_prod_beam_mom", "lambda_mass", "lambda_pion_track_id", "effective_ntTpc"};
  for (const auto* name : required) {
    if (!tree->GetBranch(name)) {
      Warning("plot_lambda_eta_yield", "Branch %s is absent in %s", name, path.Data());
      return;
    }
  }
  const Bool_t hasVetoStatus = tree->GetBranch("pion_beam_window_status") != nullptr;
  if (!hasVetoStatus)
    Warning("plot_lambda_eta_yield", "pion_beam_window_status absent in %s; using the Lambda candidate tree selection", path.Data());

  TTreeReader reader(tree);
  TTreeReaderValue<Int_t> effectiveNtTpc(reader, "effective_ntTpc");
  TTreeReaderValue<std::vector<Double_t>> xMass(reader, "X_mass");
  TTreeReaderValue<std::vector<Int_t>> productionFound(reader, "X_prod_found");
  TTreeReaderValue<std::vector<Double_t>> prodX(reader, "X_prod_vtx_x");
  TTreeReaderValue<std::vector<Double_t>> prodY(reader, "X_prod_vtx_y");
  TTreeReaderValue<std::vector<Double_t>> prodZ(reader, "X_prod_vtx_z");
  TTreeReaderValue<std::vector<Double_t>> productionBeamMomentum(reader, "X_prod_beam_mom");
  TTreeReaderValue<std::vector<Double_t>> lambdaMass(reader, "lambda_mass");
  TTreeReaderValue<std::vector<Int_t>> pionTrackId(reader, "lambda_pion_track_id");
  std::unique_ptr<TTreeReaderValue<std::vector<Int_t>>> vetoStatus;
  if (hasVetoStatus)
    vetoStatus = std::make_unique<TTreeReaderValue<std::vector<Int_t>>>(reader, "pion_beam_window_status");

  while (reader.Next()) {
    const Int_t effectiveCount = *effectiveNtTpc;
    if (effectiveCount < 2) continue;
    std::vector<Int_t> categoryCodes = {3}; // effective nTPC >= 2
    if (effectiveCount == 2) categoryCodes.push_back(2); // exactly 2
    if (effectiveCount >= 3) categoryCodes.push_back(4); // effective nTPC >= 3
    const auto n = std::min({xMass->size(), productionFound->size(), prodX->size(), prodY->size(), prodZ->size(),
                             productionBeamMomentum->size(), lambdaMass->size(), pionTrackId->size()});
    for (const Int_t categoryCode : categoryCodes) {
      auto categoryIt = categories.find(categoryCode);
      if (categoryIt == categories.end()) {
        categoryIt = categories.emplace(categoryCode,
          std::make_unique<EffectiveNtTpcHistograms>(categoryCode)).first;
      }
      EffectiveNtTpcHistograms& category = *categoryIt->second;
      ++category.eventsRead;
      for (std::size_t i=0; i<n; ++i) {
      ++category.candidatesRead;
      if (!productionFound->at(i) || !std::isfinite(lambdaMass->at(i))) continue;
      if (kRequireProductionVertexInsideLH2 &&
          !IsInsideLH2(prodX->at(i), prodY->at(i), prodZ->at(i))) continue;
      if (hasVetoStatus) {
        const Int_t pionId = pionTrackId->at(i);
        if (pionId < 0 || static_cast<std::size_t>(pionId) >= (*vetoStatus)->size() ||
            (*vetoStatus)->at(pionId) != 0) continue;
      }
      // Keep this spectrum before applying the Lambda invariant-mass window.
      category.lambdaMassBeforeCut.Fill(lambdaMass->at(i));
      ++category.lambdaCandidatesBeforeCut;
      if (!std::isfinite(xMass->at(i)) ||
          !(lambdaMass->at(i) > kLambdaMassMin && lambdaMass->at(i) < kLambdaMassMax)) continue;
      category.missingMass.Fill(xMass->at(i));
      const Double_t beamMomentum = productionBeamMomentum->at(i);
      if (std::isfinite(beamMomentum) && beamMomentum >= kBeamMomentumMin &&
          beamMomentum < kBeamMomentumMin + kBeamMomentumBins*kBeamMomentumBinWidth) {
        const Int_t beamBin = static_cast<Int_t>((beamMomentum-kBeamMomentumMin)/kBeamMomentumBinWidth);
        if (beamBin >= 0 && beamBin < static_cast<Int_t>(category.missingMassByBeamMomentum.size())) {
          category.missingMassByBeamMomentum.at(beamBin)->Fill(xMass->at(i));
          category.selectedBeamMomentumSum += beamMomentum;
          ++category.selectedBeamMomentumCount;
          category.selectedBeamMomentumSumsByBin.at(beamBin) += beamMomentum;
          ++category.selectedBeamMomentumCountsByBin.at(beamBin);
          category.missingMassVsProductionBeamMomentum.Fill(xMass->at(i), beamMomentum);
        }
      }
      ++category.candidatesSelected;
      }
    }
  }
}
}

void plot_lambda_eta_yield()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(kOutputDirectory, kTRUE);
  std::map<Int_t, std::unique_ptr<EffectiveNtTpcHistograms>> categories;
  // One pass through each run's TTree. All yield inputs are filled together.
  for (const Int_t run : kRunNumbers)
    FillRunOnce(run, categories);

  if (categories.empty()) {
    Warning("plot_lambda_eta_yield", "No events with effective_ntTpc >= 2 were found");
    return;
  }

  for (auto& categoryEntry : categories) {
  EffectiveNtTpcHistograms& category = *categoryEntry.second;
  TH1D& missingMass = category.missingMass;
  TH1D& lambdaMassBeforeCut = category.lambdaMassBeforeCut;
  auto& missingMassByBeamMomentum = category.missingMassByBeamMomentum;
  Double_t& selectedBeamMomentumSum = category.selectedBeamMomentumSum;
  Long64_t& selectedBeamMomentumCount = category.selectedBeamMomentumCount;
  auto& selectedBeamMomentumSumsByBin = category.selectedBeamMomentumSumsByBin;
  auto& selectedBeamMomentumCountsByBin = category.selectedBeamMomentumCountsByBin;
  Long64_t eventsRead = category.eventsRead;
  Long64_t candidatesRead = category.candidatesRead;
  Long64_t lambdaCandidatesBeforeCut = category.lambdaCandidatesBeforeCut;
  Long64_t candidatesSelected = category.candidatesSelected;
  const Bool_t exactlyTwoEffectiveTracks = category.effectiveNtTpc == 2;
  const Bool_t atLeastTwoEffectiveTracks = category.effectiveNtTpc == 3;
  const TString categoryTag = exactlyTwoEffectiveTracks ? "effectiveNtTpc2" :
    (atLeastTwoEffectiveTracks ? "effectiveNtTpc2plus" : "effectiveNtTpc3plus");
  const TString pdfName = Form("%s/lambda_eta_yield_mom%d_%s_%s.pdf",
                               kOutputDirectory.Data(), kBeamMomentum,
                               RunRangeTag().Data(), categoryTag.Data());
  TString rootName = pdfName;
  rootName.ReplaceAll(".pdf", ".root");

  TF1 lambdaPeakFit("lambda_peak_with_linear_background", LambdaPeakWithLinearBackground,
                    kLambdaPeakFitMin, kLambdaPeakFitMax, 5);
  const LambdaPeakFitResult lambdaFit = FitLambdaPeak(lambdaMassBeforeCut, lambdaPeakFit);

  const Double_t averageBeamMomentum = selectedBeamMomentumCount > 0
    ? selectedBeamMomentumSum/selectedBeamMomentumCount : 0.001*kBeamMomentum;
  const Double_t totalKinematicEndpoint = KinematicMissingMassEndpoint(averageBeamMomentum);
  TF1 etaPeakFit("eta_gaussian_plus_linear_background", EtaGaussianPlusLinearBackground,
                 kMissingMassMin, std::min(kEtaPeakFitMax, totalKinematicEndpoint), 6);
  etaPeakFit.SetNpx(1200); // smooth rendering; does not change fitted parameters
  const EtaPeakFitResult totalEtaPeak = FitEtaPeak(
    missingMass, etaPeakFit, totalKinematicEndpoint);
  std::vector<EtaPeakFitResult> beamEtaPeakResults;
  beamEtaPeakResults.reserve(kBeamMomentumBins);
  std::vector<std::unique_ptr<TF1>> beamEtaPeakFunctions;
  beamEtaPeakFunctions.reserve(kBeamMomentumBins);
  std::vector<Double_t> beamKinematicEndpoints(kBeamMomentumBins, 0.);
  TH1D etaYieldByBeamMomentum("h_lambda_eta_yield_vs_production_beam_momentum", "",
                               kBeamMomentumBins, kBeamMomentumMin,
                               kBeamMomentumMin+kBeamMomentumBins*kBeamMomentumBinWidth);
  etaYieldByBeamMomentum.SetDirectory(nullptr);
  etaYieldByBeamMomentum.Sumw2();
  TH1D etaPeakMassByBeamMomentum("h_eta_peak_mass_vs_production_beam_momentum", "",
                                  kBeamMomentumBins, kBeamMomentumMin,
                                  kBeamMomentumMin+kBeamMomentumBins*kBeamMomentumBinWidth);
  etaPeakMassByBeamMomentum.SetDirectory(nullptr);
  etaPeakMassByBeamMomentum.Sumw2();
  for (Int_t ibin=0; ibin<kBeamMomentumBins; ++ibin) {
    beamEtaPeakFunctions.emplace_back(std::make_unique<TF1>(
      Form("eta_gaussian_plus_linear_background_beam_bin_%02d", ibin), EtaGaussianPlusLinearBackground,
      kMissingMassMin, kEtaPeakFitMax, 6));
    beamEtaPeakFunctions.back()->SetNpx(1200);
    const Double_t averageBeamMomentumInBin = selectedBeamMomentumCountsByBin.at(ibin) > 0
      ? selectedBeamMomentumSumsByBin.at(ibin)/selectedBeamMomentumCountsByBin.at(ibin)
      : kBeamMomentumMin + (ibin+0.5)*kBeamMomentumBinWidth;
    beamKinematicEndpoints.at(ibin) = KinematicMissingMassEndpoint(averageBeamMomentumInBin);
    const EtaPeakFitResult peak = FitEtaPeak(*missingMassByBeamMomentum.at(ibin),
      *beamEtaPeakFunctions.back(), beamKinematicEndpoints.at(ibin));
    beamEtaPeakResults.push_back(peak);
    if (peak.fit_status == 0) {
      // The momentum-dependent Lambda-eta yield is the fitted Gaussian signal area.
      etaYieldByBeamMomentum.SetBinContent(ibin+1, peak.yield);
      etaYieldByBeamMomentum.SetBinError(ibin+1, peak.yield_error);
      etaPeakMassByBeamMomentum.SetBinContent(ibin+1, peak.mass);
      etaPeakMassByBeamMomentum.SetBinError(ibin+1, peak.mass_error);
    }
  }

  TCanvas canvas("c_lambda_eta_yield", "LambdaEta yield", 1000, 800);
  canvas.Print(pdfName + "[");
  TPaveText title(.08, .08, .92, .92, "NDC");
  title.SetBorderSize(0);
  title.SetFillStyle(0);
  title.SetTextAlign(12);
  title.SetTextSize(.028);
  TDatime now;
  title.AddText("LambdaEta yield from production-vertex missing mass");
  title.AddText(Form("Macro: plot_lambda_eta_yield.cc | Executed: %04d-%02d-%02d %02d:%02d:%02d",
                     now.GetYear(), now.GetMonth(), now.GetDay(), now.GetHour(), now.GetMinute(), now.GetSecond()));
  title.AddText(Form("Beam momentum setting: %d MeV/c", kBeamMomentum));
  title.AddText(exactlyTwoEffectiveTracks ? "Effective nTPC = 2" :
                (atLeastTwoEffectiveTracks ? "Effective nTPC >= 2" : "Effective nTPC >= 3"));
  TString runList = "Runs used: ";
  for (const Int_t run : kRunNumbers) runList += Form("%s%05d", runList == "Runs used: " ? "" : ", ", run);
  title.AddText(runList);
  title.AddText(Form("Input: %s/runXXXXX_DstTPCLambdaEta.root", kDataDirectory.Data()));
  title.AddText("Selection: X_prod_found; Lambda mass window 1.09 < M(p pi-) < 1.13 GeV/c^2");
  title.AddText(Form("Lambda mass displayed: %.2f to %.2f GeV/c^2; peak fit: %.2f to %.2f GeV/c^2",
                     kLambdaMassDisplayMin, kLambdaMassDisplayMax, kLambdaPeakFitMin, kLambdaPeakFitMax));
  title.AddText(Form("Pre-cut Lambda peak fit model: Gaussian + linear background, %.3f to %.3f GeV/c^2",
                     kLambdaPeakFitMin, kLambdaPeakFitMax));
  title.AddText(Form("Production vertex LH2 cut: %s", kRequireProductionVertexInsideLH2 ? "on" : "off"));
  title.AddText("Pion selection: lambda_pion_track_id matched to pion_beam_window_status == 0 when branch exists");
  title.AddText(Form("Combined eta fit: Gaussian + linear background up to kinematic endpoint, %.2f to %.2f GeV/c^2",
                     kEtaPeakFitMin, kEtaPeakFitMax));
  title.AddText(Form("Linear background fit: %.2f--%.2f GeV/c^{2}; fixed in the eta-peak fit",
                     kEtaBackgroundFitMin, kEtaBackgroundFitMax));
  title.AddText(Form("Selected-candidate mean beam momentum = %.4f GeV/c; kinematic endpoint = %.4f GeV/c^2",
                     averageBeamMomentum, totalKinematicEndpoint));
  title.AddText("Momentum-dependent eta yield: Gaussian signal area from the combined fit");
  title.AddText(Form("Production beam momentum: %d bins x %.0f MeV/c from %.0f to %.0f MeV/c",
                     kBeamMomentumBins, 1000.*kBeamMomentumBinWidth, 1000.*kBeamMomentumMin,
                     1000.*(kBeamMomentumMin+kBeamMomentumBins*kBeamMomentumBinWidth)));
  title.AddText(Form("Events read: %lld | candidates checked: %lld | pre-cut Lambda entries: %lld | eta selection: %lld",
                     eventsRead, candidatesRead, lambdaCandidatesBeforeCut, candidatesSelected));
  title.Draw();
  canvas.Print(pdfName);

  canvas.Clear();
  lambdaMassBeforeCut.SetTitle("#Lambda invariant mass before Lambda-mass cut;M(p#pi^{-}) [GeV/c^{2}];Counts");
  lambdaMassBeforeCut.SetLineColor(kBlue+1);
  lambdaMassBeforeCut.SetLineWidth(2);
  lambdaMassBeforeCut.GetXaxis()->SetRangeUser(kLambdaMassDisplayMin, kLambdaMassDisplayMax);
  lambdaMassBeforeCut.Draw("hist");
  lambdaPeakFit.SetLineColor(kRed+1);
  lambdaPeakFit.SetLineWidth(2);
  lambdaPeakFit.Draw("same");
  TPaveText lambdaFitBox(.57, .68, .89, .88, "NDC");
  lambdaFitBox.SetBorderSize(0);
  lambdaFitBox.SetFillStyle(0);
  lambdaFitBox.SetTextAlign(12);
  lambdaFitBox.AddText(Form("M_{#Lambda} = %.5f #pm %.5f GeV/c^{2}", lambdaFit.mass, lambdaFit.mass_error));
  lambdaFitBox.AddText(Form("#sigma = %.5f #pm %.5f GeV/c^{2}", lambdaFit.sigma, lambdaFit.sigma_error));
  lambdaFitBox.AddText(Form("Fit status: %d", lambdaFit.fit_status));
  lambdaFitBox.Draw();
  canvas.Print(pdfName);

  canvas.Clear();
  missingMass.SetTitle(Form("LambdaEta missing mass;M_{X} [GeV/c^{2}];Counts / %.1f MeV/c^{2}",
                            1000.*(kMissingMassMax-kMissingMassMin)/kMissingMassBins));
  missingMass.SetLineColor(kBlue+1);
  missingMass.SetLineWidth(2);
  missingMass.GetXaxis()->SetRangeUser(0., kMissingMassMax);
  missingMass.Draw("hist");
  TLegend totalFitLegend(.19, .49, .52, .67);
  totalFitLegend.SetBorderSize(0);
  totalFitLegend.SetFillStyle(0);
  totalFitLegend.AddEntry(&missingMass, "Data", "l");
  TF1 totalEtaSignal("total_eta_signal_component", EtaGaussianComponent,
                     kEtaPeakFitMin, kEtaPeakFitMax, 3);
  TF1 totalEtaContinuum("total_eta_linear_background_component",
                        EtaLinearBackgroundComponent, kMissingMassMin,
                        totalKinematicEndpoint, 3);
  totalEtaSignal.SetNpx(1200);
  totalEtaContinuum.SetNpx(1200);
  if (totalEtaPeak.fit_status == 0) {
    totalEtaSignal.SetParameters(etaPeakFit.GetParameter(0), etaPeakFit.GetParameter(1),
                                 etaPeakFit.GetParameter(2));
    totalEtaSignal.SetLineColor(kMagenta+1);
    totalEtaSignal.SetLineStyle(2);
    totalEtaSignal.SetLineWidth(2);
    totalEtaSignal.Draw("same");
    totalEtaContinuum.SetParameters(etaPeakFit.GetParameter(3), etaPeakFit.GetParameter(4),
                                    etaPeakFit.GetParameter(5));
    totalEtaContinuum.SetLineColor(kGreen+2);
    totalEtaContinuum.SetLineStyle(2);
    totalEtaContinuum.SetLineWidth(2);
    totalEtaContinuum.Draw("same");
    etaPeakFit.SetLineColor(kBlack);
    etaPeakFit.SetLineWidth(2);
    etaPeakFit.Draw("same");
    totalFitLegend.AddEntry(&etaPeakFit, "Signal + endpoint background fit", "l");
    totalFitLegend.AddEntry(&totalEtaSignal, "#eta peak component", "l");
    totalFitLegend.AddEntry(&totalEtaContinuum, "Linear background component", "l");
  }
  totalFitLegend.Draw();
  TPaveText yieldBox(.19, .70, .52, .91, "NDC");
  yieldBox.SetBorderSize(0);
  yieldBox.SetFillStyle(0);
  yieldBox.SetTextAlign(12);
  if (totalEtaPeak.fit_status == 0) {
    yieldBox.AddText(Form("Combined #eta fit: M = %.5f #pm %.5f", totalEtaPeak.mass, totalEtaPeak.mass_error));
    yieldBox.AddText(Form("#sigma = %.5f; N_{#eta}^{fit} = %.1f #pm %.1f",
                          totalEtaPeak.sigma, totalEtaPeak.yield, totalEtaPeak.yield_error));
  }
  yieldBox.Draw();
  canvas.Print(pdfName);

  canvas.Clear();
  category.missingMassVsProductionBeamMomentum.SetTitle(
    "Production-vertex K^{-} momentum vs #Lambda#eta missing mass;M_{X} [GeV/c^{2}];"
    "|#vec{p}_{K^{-}}^{prod}| [GeV/c]");
  category.missingMassVsProductionBeamMomentum.GetXaxis()->SetRangeUser(kMissingMassMin,
                                                                         kMissingMassMax);
  category.missingMassVsProductionBeamMomentum.GetYaxis()->SetRangeUser(kBeamMomentumMin,
                                                                         kBeamMomentumMax);
  category.missingMassVsProductionBeamMomentum.Draw("COLZ");
  canvas.Print(pdfName);

  /* Removed daughter-scale scan; momentum-bin analysis replaces it.
  for (std::size_t i=0; i<kDaughterMomentumScales.size(); ++i) {
    const Double_t scale = kDaughterMomentumScales[i];
    canvas.Clear();
    TH1D& lambdaHistogram = *scaledLambdaMass[i];
    lambdaHistogram.SetTitle(Form("Scaled #Lambda mass, common daughter momentum scale = %.2f;M(p#pi^{-}) [GeV/c^{2}];Counts",
                                  scale));
    lambdaHistogram.SetLineColor(kBlue+1);
    lambdaHistogram.SetLineWidth(2);
    lambdaHistogram.GetXaxis()->SetRangeUser(kLambdaMassDisplayMin, kLambdaMassDisplayMax);
    lambdaHistogram.Draw("hist");
    if (scaledLambdaPeakResults[i].fit_status == 0) {
      scaledLambdaFits[i]->SetLineColor(kRed+1);
      scaledLambdaFits[i]->SetLineWidth(2);
      scaledLambdaFits[i]->Draw("same");
    }
    TPaveText lambdaScaleBox(.57, .68, .89, .88, "NDC");
    lambdaScaleBox.SetBorderSize(0);
    lambdaScaleBox.SetFillStyle(0);
    lambdaScaleBox.SetTextAlign(12);
    lambdaScaleBox.AddText(Form("Scale = %.2f", scale));
    if (scaledLambdaPeakResults[i].fit_status == 0) {
      lambdaScaleBox.AddText(Form("M_{#Lambda} = %.5f #pm %.5f GeV/c^{2}",
                                  scaledLambdaPeakResults[i].mass,
                                  scaledLambdaPeakResults[i].mass_error));
      lambdaScaleBox.AddText(Form("#sigma = %.5f GeV/c^{2}", scaledLambdaPeakResults[i].sigma));
    }
    else lambdaScaleBox.AddText(Form("Fit failed/status %d", scaledLambdaPeakResults[i].fit_status));
    lambdaScaleBox.Draw();
    canvas.Print(pdfName);

    canvas.Clear();
    TH1D& etaHistogram = *scaledMissingMass[i];
    etaHistogram.SetTitle(Form("Scaled missing mass, common daughter momentum scale = %.2f;M_{X} [GeV/c^{2}];Counts",
                               scale));
    etaHistogram.SetLineColor(kBlue+1);
    etaHistogram.SetLineWidth(2);
    etaHistogram.GetXaxis()->SetRangeUser(0.45, 0.60);
    etaHistogram.Draw("hist");
    if (scaledEtaPeakResults[i].fit_status == 0) {
      scaledEtaFits[i]->SetLineColor(kRed+1);
      scaledEtaFits[i]->SetLineWidth(2);
      scaledEtaFits[i]->Draw("same");
    }
    TPaveText etaScaleBox(.57, .68, .89, .88, "NDC");
    etaScaleBox.SetBorderSize(0);
    etaScaleBox.SetFillStyle(0);
    etaScaleBox.SetTextAlign(12);
    etaScaleBox.AddText(Form("Scale = %.2f", scale));
    etaScaleBox.AddText("Background from linear fit, 0-0.5 GeV/c^{2}");
    if (scaledEtaPeakResults[i].fit_status == 0) {
      etaScaleBox.AddText(Form("M_{#eta} = %.5f #pm %.5f GeV/c^{2}",
                               scaledEtaPeakResults[i].mass,
                               scaledEtaPeakResults[i].mass_error));
      etaScaleBox.AddText(Form("#sigma = %.5f GeV/c^{2}", scaledEtaPeakResults[i].sigma));
    }
    else etaScaleBox.AddText(Form("Fit failed/status %d", scaledEtaPeakResults[i].fit_status));
    etaScaleBox.Draw();
    canvas.Print(pdfName);
  }

  canvas.Clear();
  if (lambdaPeakVsScale.GetN() > 0) {
    lambdaPeakVsScale.SetMarkerStyle(20);
    lambdaPeakVsScale.SetMarkerSize(1.2);
    lambdaPeakVsScale.SetMinimum(kLambdaScalePeakFitMin);
    lambdaPeakVsScale.SetMaximum(kLambdaScalePeakFitMax);
    lambdaPeakVsScale.Draw("AP");
    lambdaPeakVsScale.GetXaxis()->SetLimits(kDaughterMomentumScales.front()-0.01,
                                             kDaughterMomentumScales.back()+0.01);
    TLine pdgLine(kDaughterMomentumScales.front()-0.01, kLambdaPeakMassInitial,
                  kDaughterMomentumScales.back()+0.01, kLambdaPeakMassInitial);
    pdgLine.SetLineColor(kRed+1);
    pdgLine.SetLineStyle(2);
    pdgLine.Draw("same");
  }
  else {
    TPaveText noLambdaFits(.15, .4, .85, .6, "NDC");
    noLambdaFits.SetBorderSize(0);
    noLambdaFits.AddText("No successful scaled Lambda peak fits");
    noLambdaFits.Draw();
  }
  canvas.Print(pdfName);

  canvas.Clear();
  if (etaPeakVsScale.GetN() > 0) {
    etaPeakVsScale.SetMarkerStyle(20);
    etaPeakVsScale.SetMarkerSize(1.2);
    etaPeakVsScale.SetMinimum(kEtaScalePeakFitMin);
    etaPeakVsScale.SetMaximum(kEtaScalePeakFitMax);
    etaPeakVsScale.Draw("AP");
    etaPeakVsScale.GetXaxis()->SetLimits(kDaughterMomentumScales.front()-0.01,
                                         kDaughterMomentumScales.back()+0.01);
    TLine etaPdgLine(kDaughterMomentumScales.front()-0.01, kEtaMassInitial,
                     kDaughterMomentumScales.back()+0.01, kEtaMassInitial);
    etaPdgLine.SetLineColor(kRed+1);
    etaPdgLine.SetLineStyle(2);
    etaPdgLine.Draw("same");
  }
  else {
    TPaveText noEtaFits(.15, .4, .85, .6, "NDC");
    noEtaFits.SetBorderSize(0);
    noEtaFits.AddText("No successful scaled eta peak fits");
    noEtaFits.Draw();
  }
  canvas.Print(pdfName);
  */
  for (Int_t ibin=0; ibin<kBeamMomentumBins; ++ibin) {
    TH1D& mass = *missingMassByBeamMomentum.at(ibin);
    const Double_t pLow = kBeamMomentumMin + ibin*kBeamMomentumBinWidth;
    const Double_t pHigh = pLow + kBeamMomentumBinWidth;
    canvas.Clear();
    mass.SetTitle(Form("Production-vertex K^{-} momentum %.0f--%.0f MeV/c;M_{X} [GeV/c^{2}];Counts / %.1f MeV/c^{2}",
                       1000.*pLow, 1000.*pHigh,
                       1000.*(kMissingMassMax-kMissingMassMin)/kMissingMassBins));
    mass.SetLineColor(kBlue+1);
    mass.SetLineWidth(2);
    mass.GetXaxis()->SetRangeUser(kMissingMassMin, kMissingMassMax);
    mass.Draw("hist");
    const EtaPeakFitResult& peak = beamEtaPeakResults.at(ibin);
    TLegend binFitLegend(.19, .49, .52, .67);
    binFitLegend.SetBorderSize(0);
    binFitLegend.SetFillStyle(0);
    binFitLegend.AddEntry(&mass, "Data", "l");
    TF1 binEtaSignal(Form("eta_signal_component_bin_%02d", ibin), EtaGaussianComponent,
                     kEtaPeakFitMin, kEtaPeakFitMax, 3);
    TF1 binEtaContinuum(Form("eta_linear_background_component_bin_%02d", ibin),
                        EtaLinearBackgroundComponent, kMissingMassMin,
                        beamKinematicEndpoints.at(ibin), 3);
    binEtaSignal.SetNpx(1200);
    binEtaContinuum.SetNpx(1200);
    if (peak.fit_status == 0) {
      binEtaSignal.SetParameters(beamEtaPeakFunctions.at(ibin)->GetParameter(0),
                                 beamEtaPeakFunctions.at(ibin)->GetParameter(1),
                                 beamEtaPeakFunctions.at(ibin)->GetParameter(2));
      binEtaSignal.SetLineColor(kMagenta+1);
      binEtaSignal.SetLineStyle(2);
      binEtaSignal.SetLineWidth(2);
      binEtaSignal.Draw("same");
      binEtaContinuum.SetParameters(beamEtaPeakFunctions.at(ibin)->GetParameter(3),
                                    beamEtaPeakFunctions.at(ibin)->GetParameter(4),
                                    beamEtaPeakFunctions.at(ibin)->GetParameter(5));
      binEtaContinuum.SetLineColor(kGreen+2);
      binEtaContinuum.SetLineStyle(2);
      binEtaContinuum.SetLineWidth(2);
      binEtaContinuum.Draw("same");
      beamEtaPeakFunctions.at(ibin)->SetLineColor(kBlack);
      beamEtaPeakFunctions.at(ibin)->SetLineWidth(2);
      beamEtaPeakFunctions.at(ibin)->Draw("same");
      binFitLegend.AddEntry(beamEtaPeakFunctions.at(ibin).get(), "Signal + endpoint background fit", "l");
      binFitLegend.AddEntry(&binEtaSignal, "#eta peak component", "l");
      binFitLegend.AddEntry(&binEtaContinuum, "Linear background component", "l");
    }
    binFitLegend.Draw();
    TPaveText binNote(.19, .70, .52, .91, "NDC");
    binNote.SetBorderSize(0);
    binNote.SetFillStyle(0);
    binNote.SetTextAlign(12);
    binNote.AddText(Form("Entries = %.0f", mass.GetEntries()));
    binNote.AddText(Form("Kinematic M_{X}^{max} = %.4f GeV/c^{2}", beamKinematicEndpoints.at(ibin)));
    if (mass.GetEntries() < kMinEntriesPerBeamMassPanel)
      binNote.AddText(Form("Low statistics (< %.0f entries)", kMinEntriesPerBeamMassPanel));
    if (peak.fit_status == 0) {
      binNote.AddText(Form("M_{#eta} = %.5f #pm %.5f", peak.mass, peak.mass_error));
      binNote.AddText(Form("#sigma = %.4f; N_{#eta}^{fit} = %.1f #pm %.1f",
                           peak.sigma, peak.yield, peak.yield_error));
    }
    binNote.Draw();
    canvas.Print(pdfName);
  }

  canvas.Clear();
  etaYieldByBeamMomentum.SetTitle(Form("Gaussian-fit #Lambda#eta yield vs production-vertex K^{-} momentum;|#vec{p}_{K^{-}}^{prod}| [GeV/c];Fitted #eta yield per %.0f MeV/c",
                                      1000.*kBeamMomentumBinWidth));
  etaYieldByBeamMomentum.SetMarkerStyle(20);
  etaYieldByBeamMomentum.SetMarkerColor(kBlue+1);
  etaYieldByBeamMomentum.SetLineColor(kBlue+1);
  etaYieldByBeamMomentum.Draw("E1");
  TPaveText momentumYieldNote(.14, .72, .75, .89, "NDC");
  momentumYieldNote.SetBorderSize(0);
  momentumYieldNote.SetFillStyle(0);
  momentumYieldNote.SetTextAlign(12);
  momentumYieldNote.AddText(Form("Combined-fit range %.2f--%.2f GeV/c^{2}; Lambda mass window %.2f--%.2f GeV/c^{2}",
                                 kEtaPeakFitMin, kEtaPeakFitMax, kLambdaMassMin, kLambdaMassMax));
  momentumYieldNote.AddText(Form("Linear background fit range: %.2f--%.2f GeV/c^{2}",
                                 kEtaBackgroundFitMin, kEtaBackgroundFitMax));
  momentumYieldNote.Draw();
  canvas.Print(pdfName);

  canvas.Clear();
  etaPeakMassByBeamMomentum.SetTitle("Fitted #eta peak position vs production-vertex K^{-} momentum;|#vec{p}_{K^{-}}^{prod}| [GeV/c];Fitted M_{X} [GeV/c^{2}]");
  etaPeakMassByBeamMomentum.SetMarkerStyle(20);
  etaPeakMassByBeamMomentum.SetMarkerColor(kMagenta+1);
  etaPeakMassByBeamMomentum.SetLineColor(kMagenta+1);
  etaPeakMassByBeamMomentum.Draw("E1");
  canvas.Print(pdfName);
  canvas.Print(pdfName + "]");

  TFile output(rootName, "RECREATE");
  TParameter<Int_t>("EffectiveNtTpcCategoryCode", category.effectiveNtTpc).Write();
  missingMass.Write();
  category.missingMassVsProductionBeamMomentum.Write();
  lambdaMassBeforeCut.Write();
  lambdaPeakFit.Write();
  TParameter<Double_t>("LambdaPeakMass", lambdaFit.mass).Write();
  TParameter<Double_t>("LambdaPeakMassError", lambdaFit.mass_error).Write();
  TParameter<Double_t>("LambdaPeakSigma", lambdaFit.sigma).Write();
  TParameter<Double_t>("LambdaPeakSigmaError", lambdaFit.sigma_error).Write();
  TParameter<Int_t>("LambdaPeakFitStatus", lambdaFit.fit_status).Write();
  TParameter<Double_t>("EtaPeakFitMass", totalEtaPeak.mass).Write();
  TParameter<Double_t>("EtaPeakFitMassError", totalEtaPeak.mass_error).Write();
  TParameter<Double_t>("EtaPeakFitSigma", totalEtaPeak.sigma).Write();
  TParameter<Double_t>("EtaPeakFitYield", totalEtaPeak.yield).Write();
  TParameter<Double_t>("EtaPeakFitYieldError", totalEtaPeak.yield_error).Write();
  TParameter<Int_t>("EtaPeakFitStatus", totalEtaPeak.fit_status).Write();
  TParameter<Double_t>("SelectedMeanBeamMomentum", averageBeamMomentum).Write();
  TParameter<Double_t>("KinematicMissingMassEndpoint", totalKinematicEndpoint).Write();
  TParameter<Long64_t>("EventsRead", eventsRead).Write();
  TParameter<Long64_t>("LambdaCandidatesBeforeMassCut", lambdaCandidatesBeforeCut).Write();
  TParameter<Long64_t>("LambdaCandidatesSelected", candidatesSelected).Write();
  /* Removed daughter-scale scan output.
  lambdaPeakVsScale.Write();
  etaPeakVsScale.Write();
  for (std::size_t i=0; i<kDaughterMomentumScales.size(); ++i) {
    const Int_t scalePercent = static_cast<Int_t>(std::lround(100.*kDaughterMomentumScales[i]));
    scaledLambdaMass[i]->Write(Form("ScaledLambdaMass_Scale%04d", scalePercent));
    scaledMissingMass[i]->Write(Form("ScaledMissingMass_Scale%04d", scalePercent));
    scaledLambdaFits[i]->Write();
    scaledEtaFits[i]->Write();
    scaledEtaSidebandFits[i]->Write();
    TParameter<Double_t>(Form("LambdaPeakMass_Scale%04d", scalePercent),
                         scaledLambdaPeakResults[i].mass).Write();
    TParameter<Double_t>(Form("LambdaPeakMassError_Scale%04d", scalePercent),
                         scaledLambdaPeakResults[i].mass_error).Write();
    TParameter<Int_t>(Form("LambdaPeakFitStatus_Scale%04d", scalePercent),
                      scaledLambdaPeakResults[i].fit_status).Write();
    TParameter<Double_t>(Form("EtaPeakMass_Scale%04d", scalePercent),
                         scaledEtaPeakResults[i].mass).Write();
    TParameter<Double_t>(Form("EtaPeakMassError_Scale%04d", scalePercent),
                         scaledEtaPeakResults[i].mass_error).Write();
    TParameter<Int_t>(Form("EtaPeakFitStatus_Scale%04d", scalePercent),
                      scaledEtaPeakResults[i].fit_status).Write();
  }
  */
  etaYieldByBeamMomentum.Write();
  etaPeakMassByBeamMomentum.Write();
  for (Int_t ibin=0; ibin<kBeamMomentumBins; ++ibin) {
    missingMassByBeamMomentum.at(ibin)->Write();
    beamEtaPeakFunctions.at(ibin)->Write();
    TParameter<Double_t>(Form("EtaPeakFitMass_beam_bin_%02d", ibin),
                         beamEtaPeakResults.at(ibin).mass).Write();
    TParameter<Double_t>(Form("EtaPeakFitMassError_beam_bin_%02d", ibin),
                         beamEtaPeakResults.at(ibin).mass_error).Write();
    TParameter<Double_t>(Form("EtaPeakFitSigma_beam_bin_%02d", ibin),
                         beamEtaPeakResults.at(ibin).sigma).Write();
    TParameter<Double_t>(Form("EtaPeakFitYield_beam_bin_%02d", ibin),
                         beamEtaPeakResults.at(ibin).yield).Write();
    TParameter<Double_t>(Form("EtaPeakFitYieldError_beam_bin_%02d", ibin),
                         beamEtaPeakResults.at(ibin).yield_error).Write();
    TParameter<Int_t>(Form("EtaPeakFitStatus_beam_bin_%02d", ibin),
                      beamEtaPeakResults.at(ibin).fit_status).Write();
    TParameter<Double_t>(Form("KinematicMissingMassEndpoint_beam_bin_%02d", ibin),
                         beamKinematicEndpoints.at(ibin)).Write();
  }
  output.Close();

  Info("plot_lambda_eta_yield", "Gaussian-fit eta yield = %.1f +/- %.1f; PDF: %s; ROOT: %s",
       totalEtaPeak.yield, totalEtaPeak.yield_error, pdfName.Data(), rootName.Data());
  }
}
