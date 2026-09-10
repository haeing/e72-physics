// Edit kRunNumbers here, then run: root -l -b -q plot_lambda_missing_mass_vertex.cc
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <memory>
#include <utility>
#include <sstream>
#include <vector>

#include <TCanvas.h>
#include <TError.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TPaveText.h>

namespace {
// Change only this value (e.g. 715 or 755): data, acceptance, and PDF
// names all follow it.  Add/remove run numbers below; no function arguments.
constexpr Int_t kBeamMomentum = 735; // MeV/c
           const std::vector<int> kRunNumbers = {2447, 2449, 2450, 2451, 2452, 2453,2454,2456,2457,2458,2459,2460,2462,2463,2465,2468}; //735
  //    const std::vector<int> kRunNumbers = {2682,2683,2684,2686,2687,2689,2690,2691};//715
const TString kDataDir = Form("/gpfs/home/had/haein/data/JPARC2025Nov_root/physics-%d", kBeamMomentum);
// Display and fit settings for the #eta missing-mass peak (GeV/c^2).
constexpr Double_t kMissingMassDisplayMax = 0.60;
constexpr Double_t kEtaMass = 0.547862;
constexpr Double_t kEtaPeakFitHalfWidth = 0.008;       // Gaussian core: +/- 8 MeV/c^2
constexpr Double_t kBackgroundFitMin = 0.25;             // GeV/c^2; excludes the pi0 peak
constexpr Double_t kBackgroundFitMax = 0.50;             // GeV/c^2; eta peak is outside this fit range
constexpr Double_t kMinEntriesPerBeamMassPanel = 20.;
constexpr Double_t kEtaYieldMin = 0.50;                  // GeV/c^2
constexpr Double_t kEtaYieldMax = 0.60;                  // GeV/c^2
// Production-angle data selection.
constexpr Double_t kMissingMassMin = 0.50; // GeV/c^2, production-vertex X_mass
// Lambda invariant-mass window for optional missing-mass background suppression.
constexpr Double_t kLambdaMassMin = 1.09; // GeV/c^2
constexpr Double_t kLambdaMassMax = 1.13; // GeV/c^2
// Toggle the reconstructed production-vertex fiducial cut for all selected plots.
constexpr Bool_t kRequireLH2Inside = true;
constexpr Double_t kTargetZ = -143.;       // mm
constexpr Double_t kLH2Radius = 38.;       // mm, strict interior
constexpr Double_t kLH2HalfLengthY = 48.;  // mm, strict interior
constexpr Int_t kCosThetaBins = 20;
// Beam-momentum study: 4 MeV/c bins from 0.70 to 0.80 GeV/c.
// The 735-MeV/c trigger acceptance is intentionally used for every bin.
constexpr Double_t kBeamMomentumBinWidth = 0.004; // GeV/c
constexpr Double_t kBeamMomentumMin = 0.700;      // GeV/c
constexpr Int_t kBeamMomentumBins = 25;
// K-Beam exposure table made by statistics/kbeam_statistics.cc.
const TString kKBeamStatisticsFile =
  "~/work/git/e72-physics/statistics/kbeam_statistics.txt";
TString AcceptanceFileName()
{
  return gSystem->ExpandPathName(Form("~/simul-data/e72-lambda-yield/mom-%d/e72_lambda_yield-%d_lambda_eta.root", kBeamMomentum,kBeamMomentum));
}

// The merged run list is part of the output name, so different selections do
// not overwrite each other's PDF (e.g. ..._runs02447_02449.pdf).
TString BuildPdfName()
{
  TString run_tag;
  for (const int run : kRunNumbers) {
    if (!run_tag.IsNull()) run_tag += "_";
    run_tag += Form("%05d", run);
  }
  return Form("lambda_missing_mass_vertex_mom%d_runs%s.pdf", kBeamMomentum, run_tag.Data());
}

// Fit the low-mass continuum below the eta peak, then extrapolate its
// linear shape beneath the eta signal for the background subtraction.
Double_t EtaBackgroundLinear(Double_t* x, Double_t* parameter)
{
  return parameter[0] + parameter[1]*x[0];
}

bool IsInsideLH2(Double_t x, Double_t y, Double_t z)
{
  return std::hypot(x, z-kTargetZ) < kLH2Radius && std::abs(y) < kLH2HalfLengthY;
}

const char* LH2SelectionLabel()
{
  return kRequireLH2Inside ? "LH2 inside" : "no LH2 cut";
}

struct BackgroundSubtractedYield {
  Double_t value = 0.;
  Double_t error = 0.;
};

BackgroundSubtractedYield ExtractEtaYield(TH1D& histogram, const TString& fit_name)
{
  TF1 background(fit_name, EtaBackgroundLinear, kBackgroundFitMin, kBackgroundFitMax, 2);
  background.SetParameters(std::max(1., histogram.GetBinContent(1)), 0.);
  histogram.Fit(&background, "RQ0");
  BackgroundSubtractedYield yield;
  const Int_t first_bin = histogram.GetXaxis()->FindFixBin(kEtaYieldMin);
  const Int_t last_bin = histogram.GetXaxis()->FindFixBin(kEtaYieldMax-1.e-9);
  Double_t variance = 0.;
  for (Int_t bin=first_bin; bin<=last_bin; ++bin) {
    yield.value += histogram.GetBinContent(bin)-background.Eval(histogram.GetBinCenter(bin));
    variance += TMath::Power(histogram.GetBinError(bin), 2);
  }
  yield.error = TMath::Sqrt(variance); // data statistical error; fit uncertainty omitted
  return yield;
}

// AnaManager.cc trigger bits: beam=16, KVC=8, HTOF multiplicity=4,
// HTOF forward=2, TPC condition=1.  Require beam and either HTOF trigger.
bool PassPhysicsTrigger(Int_t trig_flag)
{
  constexpr Int_t kBeam = 1 << 4;
  constexpr Int_t kHtofMp2 = 1 << 2;
  constexpr Int_t kHtofForward = 1 << 1;
  return (trig_flag & kBeam) && (trig_flag & (kHtofMp2 | kHtofForward));
}

struct KBeamExposure {
  Double_t selected = 0.;
  Double_t all_e72physics = 0.;
  Int_t selected_runs_found = 0;
};

// Read the run-by-run section only.  The exposure is K-Beam x DAQ-Eff so the
// recorded selected yield can be scaled to all E72Physics runs at this setting.
KBeamExposure ReadKBeamExposure()
{
  KBeamExposure exposure;
  const TString file_name = gSystem->ExpandPathName(kKBeamStatisticsFile);
  std::ifstream input(file_name.Data());
  if (!input) {
    Warning("plot_lambda_missing_mass_vertex", "Cannot open %s", file_name.Data());
    return exposure;
  }
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#' || line[0] == '[') continue;
    std::istringstream stream(line);
    Int_t momentum_mev = 0, run = 0;
    Double_t k_beam = -1., trig_d = -1., daq_efficiency = -1., effective_k_beam = -1.;
    if (!(stream >> momentum_mev >> run >> k_beam >> trig_d >> daq_efficiency >> effective_k_beam)) continue;
    if (momentum_mev != kBeamMomentum || effective_k_beam < 0.) continue;
    exposure.all_e72physics += effective_k_beam;
    if (std::find(kRunNumbers.begin(), kRunNumbers.end(), run) != kRunNumbers.end()) {
      exposure.selected += effective_k_beam;
      ++exposure.selected_runs_found;
    }
  }
  return exposure;
}

// c*tau is evaluated for every candidate with a reconstructed production
// vertex and Lambda decay vertex.  It intentionally has no eta, LH2, or
// Lambda invariant-mass selection.
void FillLambdaCTau(TFile& input, TH1D& ctau, TH1D& flight_length)
{
  auto* tree = dynamic_cast<TTree*>(input.Get("tpc"));
  if (!tree) return;
  const std::array<const char*, 10> needed{{"X_prod_found", "X_prod_vtx_x", "X_prod_vtx_y", "X_prod_vtx_z",
                                             "lambda_vtx_x", "lambda_vtx_y", "lambda_vtx_z",
                                             "lambda_mom_x", "lambda_mom_y", "lambda_mom_z"}};
  for (const auto* name : needed) if (!tree->GetBranch(name)) {
    Warning("plot_lambda_missing_mass_vertex", "%s is absent; skip c#tau for %s", name, input.GetName());
    return;
  }
  TTreeReader reader(tree);
  TTreeReaderValue<std::vector<Int_t>> production_found(reader, "X_prod_found");
  TTreeReaderValue<std::vector<Double_t>> prod_x(reader, "X_prod_vtx_x");
  TTreeReaderValue<std::vector<Double_t>> prod_y(reader, "X_prod_vtx_y");
  TTreeReaderValue<std::vector<Double_t>> prod_z(reader, "X_prod_vtx_z");
  TTreeReaderValue<std::vector<Double_t>> decay_x(reader, "lambda_vtx_x");
  TTreeReaderValue<std::vector<Double_t>> decay_y(reader, "lambda_vtx_y");
  TTreeReaderValue<std::vector<Double_t>> decay_z(reader, "lambda_vtx_z");
  TTreeReaderValue<std::vector<Double_t>> lambda_px(reader, "lambda_mom_x");
  TTreeReaderValue<std::vector<Double_t>> lambda_py(reader, "lambda_mom_y");
  TTreeReaderValue<std::vector<Double_t>> lambda_pz(reader, "lambda_mom_z");
  constexpr Double_t kLambdaMassPDG = 1.115683; // GeV/c^2
  while (reader.Next()) {
    const auto n = std::min({production_found->size(), prod_x->size(), prod_y->size(), prod_z->size(),
                             decay_x->size(), decay_y->size(), decay_z->size(),
                             lambda_px->size(), lambda_py->size(), lambda_pz->size()});
    for (std::size_t i=0; i<n; ++i) {
      if (!production_found->at(i)) continue;
      const Double_t dx = decay_x->at(i)-prod_x->at(i);
      const Double_t dy = decay_y->at(i)-prod_y->at(i);
      const Double_t dz = decay_z->at(i)-prod_z->at(i);
      const Double_t p2 = lambda_px->at(i)*lambda_px->at(i) + lambda_py->at(i)*lambda_py->at(i) + lambda_pz->at(i)*lambda_pz->at(i);
      if (!std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(dz) || !std::isfinite(p2) || p2 <= 0.) continue;
      const Double_t length = std::sqrt(dx*dx + dy*dy + dz*dz);
      if (!std::isfinite(length)) continue;
      flight_length.Fill(length);
      ctau.Fill(length*kLambdaMassPDG/std::sqrt(p2));
    }
  }
}

void FillSelectedDataAngles(TFile& input, TH1D& selected_lambda_mass,
                            TH1D& selected_lambda_mass_high_mm,
                            TH2D& production_zx_no_mm, TH2D& production_zx_high_mm,
                            TH2D& production_zy_no_mm, TH2D& production_zy_high_mm,
                            TH2D& production_xy_no_mm, TH2D& production_xy_high_mm, TH1D& lab, TH1D& angle_cm, TH1D& costheta_cm,
                            TH2D& costheta_lab_vs_lambda_mom_all,
                            TH2D& costheta_lab_vs_lambda_mom_low,
                            TH2D& costheta_lab_vs_lambda_mom_high,
                            TH2D& beam_momentum_vs_costheta,
                            const std::vector<TH1D*>& missing_mass_by_beam,
                            const std::vector<TH1D*>& missing_mass_by_beam_lambda_cut,
                            TH1D& missing_mass_no_lambda_cut,
                            TH1D& missing_mass_lambda_cut,
                            TH1D& missing_mass2_no_lambda_cut,
                            TH1D& missing_mass2_lambda_cut,
                            TH3D& missing_mass_by_beam_costheta)
{
  auto* tree = dynamic_cast<TTree*>(input.Get("tpc"));
  if (!tree) return;
  const std::array<const char*, 14> needed{{"X_mass", "X_mass2", "X_prod_found", "X_prod_vtx_x", "X_prod_vtx_y", "X_prod_vtx_z",
                                              "lambda_mass", "lambda_beam_opening_angle_lab", "lambda_production_angle_cm", "lambda_production_costheta_cm",
                                              "lambda_mom_x", "lambda_mom_y", "lambda_mom_z",
                                              "X_prod_beam_mom"}};
  for (const auto* name : needed) if (!tree->GetBranch(name)) {
    Warning("plot_lambda_missing_mass_vertex", "%s is absent; skip angle cuts for %s", name, input.GetName());
    return;
  }
  TTreeReader reader(tree);
  TTreeReaderValue<std::vector<Double_t>> x_mass(reader, "X_mass");
  TTreeReaderValue<std::vector<Double_t>> x_mass2(reader, "X_mass2");
  TTreeReaderValue<std::vector<Int_t>> production_found(reader, "X_prod_found");
  TTreeReaderValue<std::vector<Double_t>> prod_x(reader, "X_prod_vtx_x");
  TTreeReaderValue<std::vector<Double_t>> prod_y(reader, "X_prod_vtx_y");
  TTreeReaderValue<std::vector<Double_t>> prod_z(reader, "X_prod_vtx_z");
  TTreeReaderValue<std::vector<Double_t>> lambda_mass(reader, "lambda_mass");
  TTreeReaderValue<std::vector<Double_t>> lab_angle(reader, "lambda_beam_opening_angle_lab");
  TTreeReaderValue<std::vector<Double_t>> cm_angle(reader, "lambda_production_angle_cm");
  TTreeReaderValue<std::vector<Double_t>> cm_costheta(reader, "lambda_production_costheta_cm");
  TTreeReaderValue<std::vector<Double_t>> lambda_px(reader, "lambda_mom_x");
  TTreeReaderValue<std::vector<Double_t>> lambda_py(reader, "lambda_mom_y");
  TTreeReaderValue<std::vector<Double_t>> lambda_pz(reader, "lambda_mom_z");
  // Post-loss RK momentum at the beam--Lambda production closest approach.
  TTreeReaderValue<std::vector<Double_t>> prod_beam_mom(reader, "X_prod_beam_mom");
  while (reader.Next()) {
    const auto n = std::min({x_mass->size(), x_mass2->size(), production_found->size(), prod_x->size(), prod_y->size(), prod_z->size(),
                             lambda_mass->size(), lab_angle->size(), cm_angle->size(), cm_costheta->size(),
                             lambda_px->size(), lambda_py->size(), lambda_pz->size(),
                             prod_beam_mom->size()});
    for (std::size_t i=0; i<n; ++i) {
      if (!production_found->at(i) || !std::isfinite(x_mass->at(i)) ||
          (kRequireLH2Inside && !IsInsideLH2(prod_x->at(i), prod_y->at(i), prod_z->at(i)))) continue;
      const Bool_t pass_lambda_mass_window = std::isfinite(lambda_mass->at(i)) &&
        lambda_mass->at(i) > kLambdaMassMin && lambda_mass->at(i) < kLambdaMassMax;
      if (!pass_lambda_mass_window) continue;
      // All selected observables use the Lambda mass window above.
      if (std::isfinite(prod_z->at(i)) && std::isfinite(prod_x->at(i)))
        production_zx_no_mm.Fill(prod_z->at(i), prod_x->at(i));
      if (std::isfinite(prod_z->at(i)) && std::isfinite(prod_y->at(i)))
        production_zy_no_mm.Fill(prod_z->at(i), prod_y->at(i));
      if (std::isfinite(prod_x->at(i)) && std::isfinite(prod_y->at(i)))
        production_xy_no_mm.Fill(prod_x->at(i), prod_y->at(i));
      if (std::isfinite(lambda_mass->at(i))) selected_lambda_mass.Fill(lambda_mass->at(i));
      missing_mass_no_lambda_cut.Fill(x_mass->at(i));
      if (std::isfinite(x_mass2->at(i))) missing_mass2_no_lambda_cut.Fill(x_mass2->at(i));
      missing_mass_lambda_cut.Fill(x_mass->at(i));
      if (std::isfinite(x_mass2->at(i)) ) missing_mass2_lambda_cut.Fill(x_mass2->at(i));
      if (std::isfinite(prod_beam_mom->at(i))) {
        const Int_t beam_bin = static_cast<Int_t>((prod_beam_mom->at(i)-kBeamMomentumMin)/kBeamMomentumBinWidth);
        if (beam_bin >= 0 && beam_bin < static_cast<Int_t>(missing_mass_by_beam.size())) {
          missing_mass_by_beam.at(beam_bin)->Fill(x_mass->at(i));
          missing_mass_by_beam_lambda_cut.at(beam_bin)->Fill(x_mass->at(i));
          if (std::isfinite(cm_costheta->at(i)))
            missing_mass_by_beam_costheta.Fill(prod_beam_mom->at(i), cm_costheta->at(i), x_mass->at(i));
        }
      }
      const Bool_t is_high_missing_mass = x_mass->at(i) > kMissingMassMin;
      // Retain the existing 1D angular spectra as the M_X > 0.5 GeV/c^2 selection.
      if (is_high_missing_mass) {
        if (std::isfinite(prod_z->at(i)) && std::isfinite(prod_x->at(i)))
          production_zx_high_mm.Fill(prod_z->at(i), prod_x->at(i));
        if (std::isfinite(prod_z->at(i)) && std::isfinite(prod_y->at(i)))
          production_zy_high_mm.Fill(prod_z->at(i), prod_y->at(i));
        if (std::isfinite(prod_x->at(i)) && std::isfinite(prod_y->at(i)))
          production_xy_high_mm.Fill(prod_x->at(i), prod_y->at(i));
        if (std::isfinite(lambda_mass->at(i))) selected_lambda_mass_high_mm.Fill(lambda_mass->at(i));
        if (std::isfinite(lab_angle->at(i))) lab.Fill(lab_angle->at(i));
        if (std::isfinite(cm_angle->at(i))) angle_cm.Fill(cm_angle->at(i));
        if (std::isfinite(cm_costheta->at(i))) {
          costheta_cm.Fill(cm_costheta->at(i));
          if (std::isfinite(prod_beam_mom->at(i))) beam_momentum_vs_costheta.Fill(prod_beam_mom->at(i), cm_costheta->at(i));
        }
      }
      const Double_t lambda_mom = std::sqrt(lambda_px->at(i)*lambda_px->at(i) + lambda_py->at(i)*lambda_py->at(i) + lambda_pz->at(i)*lambda_pz->at(i));
      // Matches lambda-eta-phase.cc: TLorentzVector::CosTheta(), i.e. pz/|p| in lab.
      if (std::isfinite(lambda_mom) && lambda_mom > 0.) {
        const Double_t costheta_lab = lambda_pz->at(i)/lambda_mom;
        costheta_lab_vs_lambda_mom_all.Fill(costheta_lab, lambda_mom);
        if (is_high_missing_mass) costheta_lab_vs_lambda_mom_high.Fill(costheta_lab, lambda_mom);
        else costheta_lab_vs_lambda_mom_low.Fill(costheta_lab, lambda_mom);
      }
    }
  }
}

bool FillTriggerAcceptance(TH1D& generated, TH1D& triggered)
{
  const TString acceptance_file = AcceptanceFileName();
  std::unique_ptr<TFile> input(TFile::Open(acceptance_file, "READ"));
  if (!input || input->IsZombie()) { Error("plot_lambda_missing_mass_vertex", "Cannot open %s", acceptance_file.Data()); return false; }
  auto* tree = dynamic_cast<TTree*>(input->Get("g4hyptpc_light"));
  if (!tree || !tree->GetBranch("cos_theta_lambda") || !tree->GetBranch("trig_flag")) {
    Error("plot_lambda_missing_mass_vertex", "g4hyptpc_light cos_theta_lambda/trig_flag is absent"); return false;
  }
  TTreeReader reader(tree);
  TTreeReaderValue<Double_t> cos_theta_lambda(reader, "cos_theta_lambda");
  TTreeReaderValue<Int_t> trig_flag(reader, "trig_flag");
  while (reader.Next()) {
    if (!std::isfinite(*cos_theta_lambda) || *cos_theta_lambda < -1. || *cos_theta_lambda > 1.) continue;
    generated.Fill(*cos_theta_lambda);
    if (PassPhysicsTrigger(*trig_flag)) triggered.Fill(*cos_theta_lambda);
  }
  return true;
}

std::unique_ptr<TH1D> MergeMassSquaredBranch(const TString& branch_name, const TString& histogram_name)
{
  auto merged = std::make_unique<TH1D>(histogram_name,
    ";M_{X}^{2} [(GeV/c^{2})^{2}];Counts / 0.0014 (GeV/c^{2})^{2}",
    500, -0.20, 0.50);
  merged->SetDirectory(nullptr);
  for (const int run : kRunNumbers) {
    const TString path = Form("%s/run%05d_DstTPCLambdaEta.root", kDataDir.Data(), run);
    std::unique_ptr<TFile> input(TFile::Open(path, "READ"));
    if (!input || input->IsZombie()) { Warning("plot_lambda_missing_mass_vertex", "Cannot open %s", path.Data()); continue; }
    auto* tree = dynamic_cast<TTree*>(input->Get("tpc"));
    if (!tree || !tree->GetBranch(branch_name)) {
      Warning("plot_lambda_missing_mass_vertex", "%s is absent in run %05d", branch_name.Data(), run);
      continue;
    }
    TTreeReader reader(tree);
    TTreeReaderValue<std::vector<Double_t>> mass2(reader, branch_name);
    while (reader.Next()) for (const Double_t value : *mass2) if (std::isfinite(value)) merged->Fill(value);
  }
  return merged;
}

std::unique_ptr<TH1> MergeHistogram(const TString& histogram_name)
{
  std::unique_ptr<TH1> merged;
  for (const int run : kRunNumbers) {
    const TString path = Form("%s/run%05d_DstTPCLambdaEta.root", kDataDir.Data(), run);
    std::unique_ptr<TFile> input(TFile::Open(path, "READ"));
    if (!input || input->IsZombie()) {
      Warning("plot_lambda_missing_mass_vertex", "Cannot open %s", path.Data());
      continue;
    }
    auto* source = dynamic_cast<TH1*>(input->Get(histogram_name));
    if (!source) {
      Warning("plot_lambda_missing_mass_vertex", "%s is absent in run %05d", histogram_name.Data(), run);
      continue;
    }
    if (!merged) {
      merged.reset(dynamic_cast<TH1*>(source->Clone(histogram_name + "_merged")));
      merged->SetDirectory(nullptr);
    } else merged->Add(source);
  }
  return merged;
}
}

void plot_lambda_missing_mass_vertex()
{
  gStyle->SetOptStat(0);
  const TString pdf_name = BuildPdfName();
  const TString target_selection = LH2SelectionLabel();
  // Use only the RK beam momentum at the reconstructed Lambda production vertex.
  auto production = MergeHistogram("LambdaEta_MissingMass_ProductionVtx");
  // Read the corresponding production-vertex MM2 branch directly.
  auto production_mass2 = MergeMassSquaredBranch("X_mass2", "h_missing_mass2_production");
  auto dca = MergeHistogram("LambdaEta_BeamLambdaProductionDCA");
  auto production_zx = MergeHistogram("LambdaEta_ProductionVtxZX");
  auto opening_lab = MergeHistogram("LambdaEta_LambdaBeamOpeningAngleLab");
  auto angle_cm = MergeHistogram("LambdaEta_LambdaProductionAngleCM");
  auto costheta_cm = MergeHistogram("LambdaEta_LambdaProductionCosThetaCM");

  // Reconstructed flight observables: all candidates with a production and decay vertex.
  // Deliberately independent of M_X (eta), LH2, and Lambda-mass selections.
  auto lambda_flight_length_all = std::make_unique<TH1D>(
    "h_lambda_flight_length_all", "", 500, 0., 1000.);
  auto lambda_ctau_all = std::make_unique<TH1D>(
    "h_lambda_ctau_all", "", 500, 0., 500.);
  auto selected_lambda_mass = std::make_unique<TH1D>(
    "h_selected_lambda_invariant_mass", "", 240, 1.08, 1.20);
  auto selected_lambda_mass_high_mm = std::make_unique<TH1D>(
    "h_selected_lambda_invariant_mass_high_mm", "", 240, 1.08, 1.20);
  auto production_zx_no_mm = std::make_unique<TH2D>(
    "h_production_vtx_zx_no_mm", "", 500, -250., 250., 500, -250., 250.);
  auto production_zx_high_mm = std::make_unique<TH2D>(
    "h_production_vtx_zx_high_mm", "", 500, -250., 250., 500, -250., 250.);
  auto production_zy_no_mm = std::make_unique<TH2D>(
    "h_production_vtx_zy_no_mm", "", 500, -250., 250., 500, -250., 250.);
  auto production_zy_high_mm = std::make_unique<TH2D>(
    "h_production_vtx_zy_high_mm", "", 500, -250., 250., 500, -250., 250.);
  auto production_xy_no_mm = std::make_unique<TH2D>(
    "h_production_vtx_xy_no_mm", "", 500, -250., 250., 500, -250., 250.);
  auto production_xy_high_mm = std::make_unique<TH2D>(
    "h_production_vtx_xy_high_mm", "", 500, -250., 250., 500, -250., 250.);
  auto selected_lab = std::make_unique<TH1D>("h_selected_lambda_lab", "", 180, 0., 180.);
  auto selected_cm = std::make_unique<TH1D>("h_selected_lambda_cm", "", 180, 0., 180.);
  auto selected_costheta = std::make_unique<TH1D>("h_selected_lambda_costheta", "", kCosThetaBins, -1., 1.);
  auto selected_costheta_lab_vs_lambda_mom_all = std::make_unique<TH2D>("h_selected_lambda_costheta_lab_vs_lambda_mom_all", "", 100, -1., 1., 120, 0., 1.2);
  auto selected_costheta_lab_vs_lambda_mom_low = std::make_unique<TH2D>("h_selected_lambda_costheta_lab_vs_lambda_mom_low", "", 100, -1., 1., 120, 0., 1.2);
  auto selected_costheta_lab_vs_lambda_mom_high = std::make_unique<TH2D>("h_selected_lambda_costheta_lab_vs_lambda_mom_high", "", 100, -1., 1., 120, 0., 1.2);
  auto generated_costheta = std::make_unique<TH1D>("h_generated_lambda_costheta", "", kCosThetaBins, -1., 1.);
  auto beam_momentum_vs_costheta = std::make_unique<TH2D>(
    "h_selected_lambda_beam_momentum_vs_costheta", "", kBeamMomentumBins, kBeamMomentumMin,
    kBeamMomentumMin+kBeamMomentumBins*kBeamMomentumBinWidth, kCosThetaBins, -1., 1.);
  auto triggered_costheta = std::make_unique<TH1D>("h_triggered_lambda_costheta", "", kCosThetaBins, -1., 1.);
  auto missing_mass_by_beam_costheta = std::make_unique<TH3D>(
    "h_missing_mass_by_production_beam_costheta", "",
    kBeamMomentumBins, kBeamMomentumMin, kBeamMomentumMin+kBeamMomentumBins*kBeamMomentumBinWidth,
    kCosThetaBins, -1., 1., 300, 0., kMissingMassDisplayMax);
  missing_mass_by_beam_costheta->SetDirectory(nullptr);
  missing_mass_by_beam_costheta->Sumw2();
  auto missing_mass_no_lambda_cut = std::make_unique<TH1D>(
    "h_missing_mass_no_lambda_mass_cut", "", 300, 0., kMissingMassDisplayMax);
  auto missing_mass_lambda_cut = std::make_unique<TH1D>(
    "h_missing_mass_lambda_mass_cut", "", 300, 0., kMissingMassDisplayMax);
  auto missing_mass2_no_lambda_cut = std::make_unique<TH1D>(
    "h_missing_mass2_no_lambda_mass_cut", "", 500, -0.20, 0.50);
  auto missing_mass2_lambda_cut = std::make_unique<TH1D>(
    "h_missing_mass2_lambda_mass_cut", "", 500, -0.20, 0.50);
  for (auto* histogram : {missing_mass_no_lambda_cut.get(), missing_mass_lambda_cut.get(),
                          missing_mass2_no_lambda_cut.get(), missing_mass2_lambda_cut.get()}) {
    histogram->SetDirectory(nullptr);
    histogram->Sumw2();
  }
  std::vector<std::unique_ptr<TH1D>> missing_mass_by_beam_storage;
  std::vector<TH1D*> missing_mass_by_beam;
  missing_mass_by_beam_storage.reserve(kBeamMomentumBins);
  missing_mass_by_beam.reserve(kBeamMomentumBins);
  std::vector<std::unique_ptr<TH1D>> missing_mass_by_beam_lambda_cut_storage;
  std::vector<TH1D*> missing_mass_by_beam_lambda_cut;
  missing_mass_by_beam_lambda_cut_storage.reserve(kBeamMomentumBins);
  missing_mass_by_beam_lambda_cut.reserve(kBeamMomentumBins);
  for (Int_t ibeam=0; ibeam<kBeamMomentumBins; ++ibeam) {
    auto histogram = std::make_unique<TH1D>(Form("h_missing_mass_prod_beam_bin_%d", ibeam),
                                            "", 300, 0., kMissingMassDisplayMax);
    histogram->SetDirectory(nullptr);
    histogram->Sumw2();
    missing_mass_by_beam.push_back(histogram.get());
    missing_mass_by_beam_storage.push_back(std::move(histogram));
    auto cut_histogram = std::make_unique<TH1D>(Form("h_missing_mass_lambda_cut_prod_beam_bin_%d", ibeam),
                                                "", 300, 0., kMissingMassDisplayMax);
    cut_histogram->SetDirectory(nullptr);
    cut_histogram->Sumw2();
    missing_mass_by_beam_lambda_cut.push_back(cut_histogram.get());
    missing_mass_by_beam_lambda_cut_storage.push_back(std::move(cut_histogram));
  }
  for (auto* hist : {lambda_flight_length_all.get(), lambda_ctau_all.get(), selected_lambda_mass.get(), selected_lambda_mass_high_mm.get(), selected_lab.get(), selected_cm.get(), selected_costheta.get(), generated_costheta.get(), triggered_costheta.get()}) {
    hist->SetDirectory(nullptr); hist->Sumw2();
  }
  for (auto* hist : {selected_costheta_lab_vs_lambda_mom_all.get(), selected_costheta_lab_vs_lambda_mom_low.get(), selected_costheta_lab_vs_lambda_mom_high.get(), beam_momentum_vs_costheta.get(),
                     production_zx_no_mm.get(), production_zx_high_mm.get(), production_zy_no_mm.get(), production_zy_high_mm.get(), production_xy_no_mm.get(), production_xy_high_mm.get()}) {
    hist->SetDirectory(nullptr); hist->Sumw2();
  }
  for (const int run : kRunNumbers) {
    const TString path = Form("%s/run%05d_DstTPCLambdaEta.root", kDataDir.Data(), run);
    std::unique_ptr<TFile> input(TFile::Open(path, "READ"));
    if (!input || input->IsZombie()) { Warning("plot_lambda_missing_mass_vertex", "Cannot open %s", path.Data()); continue; }
    FillLambdaCTau(*input, *lambda_ctau_all, *lambda_flight_length_all);
    FillSelectedDataAngles(*input, *selected_lambda_mass, *selected_lambda_mass_high_mm,
                           *production_zx_no_mm, *production_zx_high_mm,
                           *production_zy_no_mm, *production_zy_high_mm, *production_xy_no_mm, *production_xy_high_mm, *selected_lab, *selected_cm, *selected_costheta,
                           *selected_costheta_lab_vs_lambda_mom_all,
                           *selected_costheta_lab_vs_lambda_mom_low,
                           *selected_costheta_lab_vs_lambda_mom_high, *beam_momentum_vs_costheta,
                           missing_mass_by_beam, missing_mass_by_beam_lambda_cut,
                           *missing_mass_no_lambda_cut, *missing_mass_lambda_cut,
                           *missing_mass2_no_lambda_cut, *missing_mass2_lambda_cut,
                           *missing_mass_by_beam_costheta);
  }
  // Use the selected (LH2-inside and Lambda mass-window) spectra for the
  // primary missing-mass and missing-mass-squared pages as well.
  production.reset(dynamic_cast<TH1*>(missing_mass_lambda_cut->Clone("h_missing_mass_selected_lambda_eta")));
  production->SetDirectory(nullptr);
  production_mass2.reset(dynamic_cast<TH1D*>(missing_mass2_lambda_cut->Clone("h_missing_mass2_selected_lambda_eta")));
  production_mass2->SetDirectory(nullptr);

  auto lambda_eta_yield_by_beam = std::make_unique<TH1D>(
    "h_lambda_eta_yield_background_subtracted", "", kBeamMomentumBins, kBeamMomentumMin,
    kBeamMomentumMin+kBeamMomentumBins*kBeamMomentumBinWidth);
  lambda_eta_yield_by_beam->SetDirectory(nullptr);
  lambda_eta_yield_by_beam->Sumw2();
  for (Int_t ibeam=0; ibeam<kBeamMomentumBins; ++ibeam) {
    const auto eta_yield = ExtractEtaYield(*missing_mass_by_beam_lambda_cut.at(ibeam),
                                           Form("eta_background_prod_beam_bin_%d", ibeam));
    // A negative background-subtracted fluctuation is not a physical yield.
    lambda_eta_yield_by_beam->SetBinContent(ibeam+1, std::max(0., eta_yield.value));
    lambda_eta_yield_by_beam->SetBinError(ibeam+1, eta_yield.error);
  }

  const Bool_t has_acceptance = FillTriggerAcceptance(*generated_costheta, *triggered_costheta);
  std::unique_ptr<TH1D> acceptance;
  std::unique_ptr<TH1D> selected_costheta_corrected;
  auto beam_momentum_yield_raw = std::make_unique<TH1D>(
    "h_lambda_yield_vs_beam_momentum_raw", "", kBeamMomentumBins, kBeamMomentumMin,
    kBeamMomentumMin+kBeamMomentumBins*kBeamMomentumBinWidth);
  auto beam_momentum_yield_corrected = std::make_unique<TH1D>(
    "h_lambda_yield_vs_beam_momentum_acceptance_corrected", "", kBeamMomentumBins, kBeamMomentumMin,
    kBeamMomentumMin+kBeamMomentumBins*kBeamMomentumBinWidth);
  beam_momentum_yield_raw->SetDirectory(nullptr);
  beam_momentum_yield_corrected->SetDirectory(nullptr);
  beam_momentum_yield_raw->Sumw2();
  beam_momentum_yield_corrected->Sumw2();
  std::unique_ptr<TH2D> lambda_eta_yield_by_beam_costheta;
  if (has_acceptance) {
    acceptance.reset(dynamic_cast<TH1D*>(triggered_costheta->Clone("h_trigger_acceptance_costheta")));
    acceptance->Divide(triggered_costheta.get(), generated_costheta.get(), 1., 1., "B");
    selected_costheta_corrected = std::make_unique<TH1D>(
      "h_lambda_eta_yield_costheta_background_subtracted", "", kCosThetaBins, -1., 1.);
    selected_costheta_corrected->SetDirectory(nullptr);
    selected_costheta_corrected->Sumw2();
    lambda_eta_yield_by_beam_costheta = std::make_unique<TH2D>(
      "h_lambda_eta_yield_by_beam_costheta_background_subtracted", "",
      kBeamMomentumBins, kBeamMomentumMin, kBeamMomentumMin+kBeamMomentumBins*kBeamMomentumBinWidth,
      kCosThetaBins, -1., 1.);
    lambda_eta_yield_by_beam_costheta->SetDirectory(nullptr);
    lambda_eta_yield_by_beam_costheta->Sumw2();

    // Fit M_X in each cos(theta) bin after integrating over production beam
    // momentum.  This keeps enough statistics for a stable background fit.
    for (Int_t icos=1; icos<=kCosThetaBins; ++icos) {
      std::unique_ptr<TH1D> mass_bin(missing_mass_by_beam_costheta->ProjectionZ(
        Form("h_missing_mass_costheta_%d", icos), 1, kBeamMomentumBins, icos, icos));
      mass_bin->SetDirectory(nullptr);
      const auto eta_yield = ExtractEtaYield(*mass_bin, Form("eta_background_costheta_%d", icos));
      const Double_t raw = std::max(0., eta_yield.value);
      const Double_t acc = acceptance->GetBinContent(icos);
      if (acc > 0.) {
        selected_costheta_corrected->SetBinContent(icos, raw/acc);
        selected_costheta_corrected->SetBinError(icos, eta_yield.error/acc);
      }
    }

    // The 4-MeV/c yield uses its own p_K-binned M_X fit.  Obtain its
    // effective acceptance from the selected angular distribution, rather
    // than attempting unstable fits in every (p_K, cos theta) cell.
    for (Int_t ibeam=1; ibeam<=kBeamMomentumBins; ++ibeam) {
      const Double_t raw = lambda_eta_yield_by_beam->GetBinContent(ibeam);
      const Double_t raw_error = lambda_eta_yield_by_beam->GetBinError(ibeam);
      Double_t observed = 0., acceptance_weighted = 0.;
      for (Int_t icos=1; icos<=kCosThetaBins; ++icos) {
        const Double_t count = beam_momentum_vs_costheta->GetBinContent(ibeam, icos);
        const Double_t acc = acceptance->GetBinContent(icos);
        observed += count;
        if (acc > 0.) acceptance_weighted += count/acc;
      }
      const Double_t inverse_effective_acceptance = observed > 0. ? acceptance_weighted/observed : 0.;
      const Double_t corrected = raw*inverse_effective_acceptance;
      const Double_t corrected_error = raw_error*inverse_effective_acceptance;
      beam_momentum_yield_raw->SetBinContent(ibeam, raw);
      beam_momentum_yield_raw->SetBinError(ibeam, raw_error);
      beam_momentum_yield_corrected->SetBinContent(ibeam, corrected);
      beam_momentum_yield_corrected->SetBinError(ibeam, corrected_error);
      if (acceptance_weighted > 0.) {
        for (Int_t icos=1; icos<=kCosThetaBins; ++icos) {
          const Double_t count = beam_momentum_vs_costheta->GetBinContent(ibeam, icos);
          const Double_t acc = acceptance->GetBinContent(icos);
          if (acc <= 0.) continue;
          const Double_t fraction = (count/acc)/acceptance_weighted;
          lambda_eta_yield_by_beam_costheta->SetBinContent(ibeam, icos, corrected*fraction);
          lambda_eta_yield_by_beam_costheta->SetBinError(ibeam, icos, corrected_error*fraction);
        }
      }
    }
  }

  const KBeamExposure kbeam_exposure = ReadKBeamExposure();
  const Double_t full_run_scale = kbeam_exposure.selected > 0.
    ? kbeam_exposure.all_e72physics/kbeam_exposure.selected : 0.;
  std::unique_ptr<TH1D> beam_momentum_yield_expected_all_runs;
  if (has_acceptance && full_run_scale > 0.) {
    beam_momentum_yield_expected_all_runs.reset(dynamic_cast<TH1D*>(
      beam_momentum_yield_corrected->Clone("h_lambda_eta_yield_expected_all_e72physics")));
    beam_momentum_yield_expected_all_runs->SetDirectory(nullptr);
    beam_momentum_yield_expected_all_runs->Scale(full_run_scale);
  } else if (has_acceptance) {
    Warning("plot_lambda_missing_mass_vertex", "No valid K-Beam exposure for the selected run list; skip all-run yield estimate.");
  }

  if (!production) {
    Error("plot_lambda_missing_mass_vertex", "Production-vertex missing-mass histogram not found. Re-run DstTPCLambdaEta with the production-vertex code.");
    return;
  }

  TCanvas canvas("c_missing_mass_vertex", "Lambda missing mass: production vertex", 1000, 800);
  canvas.Print(pdf_name + "[");

  canvas.Clear();
  lambda_flight_length_all->SetTitle("#Lambda flight length: production to decay vertex (no M_{X}, LH2, or M_{p#pi^{-}} cut);L [mm];Candidates");
  lambda_flight_length_all->SetLineColor(kBlue + 1);
  lambda_flight_length_all->SetLineWidth(2);
  lambda_flight_length_all->Draw("hist");
  canvas.Print(pdf_name);

  canvas.Clear();
  lambda_ctau_all->SetTitle("#Lambda c#tau from reconstructed vertices (no M_{X}, LH2, or M_{p#pi^{-}} cut);c#tau [mm];Candidates");
  lambda_ctau_all->SetLineColor(kBlue + 1);
  lambda_ctau_all->SetLineWidth(2);
  lambda_ctau_all->Draw("hist");
  canvas.Print(pdf_name);
  // Presentation spectrum. Numerical background subtraction for yields is kept separate and not drawn.
  production->SetTitle("Missing mass: LH2 inside, 1.09<M_{p#pi^{-}}<1.13;M_{X} [GeV/c^{2}];Counts / 2 MeV/c^{2}");
  // continuum with a straight line, then fit the background-subtracted eta core.
  production->SetTitle("#Lambda missing mass: LH2 inside, 1.09<M_{p#pi^{-}}<1.13;M_{X} [GeV/c^{2}];Counts / 2 MeV/c^{2}");
  production->SetLineColor(kRed + 1); production->SetLineWidth(2);
  production->GetXaxis()->SetRangeUser(0., kMissingMassDisplayMax);
  production->SetMaximum(1.15*production->GetMaximum());
  production->Draw("hist");
  canvas.Print(pdf_name);

  canvas.Clear();
  production_mass2->SetTitle("#Lambda missing mass squared: LH2 inside, 1.09<M_{p#pi^{-}}<1.13;M_{X}^{2} [(GeV/c^{2})^{2}];Counts / 0.0014 (GeV/c^{2})^{2}");
  production_mass2->SetLineColor(kRed + 1); production_mass2->SetLineWidth(2); production_mass2->Draw("hist");
  canvas.Print(pdf_name);

  // Missing-mass spectra in 4-MeV/c production-vertex K18 momentum bins.
  constexpr Int_t kMissingMassPanelsPerPage = 9;
  for (Int_t first_beam_bin=0; first_beam_bin<kBeamMomentumBins;
       first_beam_bin+=kMissingMassPanelsPerPage) {
    TCanvas missing_mass_canvas(Form("c_missing_mass_by_production_beam_%d", first_beam_bin),
                                "Missing mass by production beam momentum", 1200, 900);
    missing_mass_canvas.Divide(3, 3);
    for (Int_t panel=0; panel<kMissingMassPanelsPerPage; ++panel) {
      const Int_t ibeam = first_beam_bin + panel;
      if (ibeam >= kBeamMomentumBins) break;
      missing_mass_canvas.cd(panel+1);
      auto* mass_bin = missing_mass_by_beam_lambda_cut.at(ibeam);
      const Double_t p_low = 1000.*(kBeamMomentumMin + ibeam*kBeamMomentumBinWidth);
      const Double_t p_high = p_low + 1000.*kBeamMomentumBinWidth;
      mass_bin->SetTitle(Form("p_{K}^{prod}=%.0f--%.0f MeV/c (LH2 inside, 1.09<M_{p#pi^{-}}<1.13);M_{X} [GeV/c^{2}];Counts / 2 MeV/c^{2}", p_low, p_high));
      mass_bin->SetLineColor(kBlue + 1);
      mass_bin->SetLineWidth(2);
      mass_bin->Draw("hist");
      if (mass_bin->GetEntries() < kMinEntriesPerBeamMassPanel) {
        TPaveText low_stat_note(.18, .72, .82, .86, "NDC");
        low_stat_note.SetBorderSize(0); low_stat_note.SetFillStyle(0);
        low_stat_note.SetTextAlign(22); low_stat_note.SetTextSize(.065);
        low_stat_note.AddText(Form("low statistics: N < %.0f", kMinEntriesPerBeamMassPanel));
        low_stat_note.Draw();
        continue;
      }
      const auto eta_yield = ExtractEtaYield(*mass_bin, Form("eta_background_integral_prod_beam_%d", ibeam));
      TPaveText panel_note(.14, .76, .62, .88, "NDC");
      panel_note.SetBorderSize(0); panel_note.SetFillStyle(0);
      panel_note.SetTextAlign(12); panel_note.SetTextSize(.05);
      panel_note.AddText(Form("N_{#eta}(0.50--0.60) = %.0f #pm %.0f",
                              std::max(0., eta_yield.value), eta_yield.error));
      panel_note.Draw();
    }
    missing_mass_canvas.Print(pdf_name);
  }

  canvas.Clear();
  lambda_eta_yield_by_beam->SetTitle(Form("Background-subtracted #Lambda#eta yield vs production-vertex K18 momentum (%s, 1.09<M_{p#pi^{-}}<1.13);|#vec{p}_{K^{-}}^{prod}| [GeV/c];#Lambda#eta yield", target_selection.Data()));
  lambda_eta_yield_by_beam->SetMarkerStyle(20);
  lambda_eta_yield_by_beam->SetMarkerColor(kBlue + 1);
  lambda_eta_yield_by_beam->SetLineColor(kBlue + 1);
  lambda_eta_yield_by_beam->Draw("E1");
  TPaveText yield_note(.14, .72, .65, .88, "NDC");
  yield_note.SetBorderSize(0); yield_note.SetFillStyle(0); yield_note.SetTextAlign(12); yield_note.SetTextSize(.028);
  yield_note.AddText("M_{X} background-subtracted in 0.50--0.60 GeV/c^{2}");
  yield_note.AddText("#Lambda#eta yield: integral of data-background in 0.50--0.60 GeV/c^{2}");
  yield_note.AddText("Errors: data statistical only (fit uncertainty omitted)");
  yield_note.Draw();
  canvas.Print(pdf_name);

  if (dca) {
    canvas.Clear();
    dca->SetTitle("K18 RK beam--#Lambda production closest approach;DCA [mm];Counts");
    dca->SetLineColor(kBlack); dca->SetLineWidth(2); dca->Draw("hist");
    canvas.Print(pdf_name);
  }
  if (production_zx) {
    canvas.Clear();
    production_zx->SetTitle("Reconstructed #Lambda production vertex (DST inclusive);Z [mm];X [mm]");
    production_zx->Draw("colz");
    canvas.Print(pdf_name);
  }
  for (const auto& item : std::vector<std::pair<TH2D*, TString>>{
         {production_zx_no_mm.get(), Form("Reconstructed #Lambda production vertex (%s, no M_{X} cut);Z [mm];X [mm]", target_selection.Data())},
         {production_zx_high_mm.get(), Form("Reconstructed #Lambda production vertex (%s, M_{X}>%.2f GeV/c^{2});Z [mm];X [mm]", target_selection.Data(), kMissingMassMin)},
         {production_zy_no_mm.get(), Form("Reconstructed #Lambda production vertex (%s, no M_{X} cut);Z [mm];Y [mm]", target_selection.Data())},
         {production_zy_high_mm.get(), Form("Reconstructed #Lambda production vertex (%s, M_{X}>%.2f GeV/c^{2});Z [mm];Y [mm]", target_selection.Data(), kMissingMassMin)},
         {production_xy_no_mm.get(), Form("Reconstructed #Lambda production vertex (%s, no M_{X} cut);X [mm];Y [mm]", target_selection.Data())},
         {production_xy_high_mm.get(), Form("Reconstructed #Lambda production vertex (%s, M_{X}>%.2f GeV/c^{2});X [mm];Y [mm]", target_selection.Data(), kMissingMassMin)}}) {
    canvas.Clear();
    item.first->SetTitle(item.second);
    item.first->Draw("colz");
    canvas.Print(pdf_name);
  }
  for (const auto& item : std::vector<std::pair<TH1*, TString>>{
         {opening_lab.get(), "K^{-}-#Lambda opening angle (lab);#alpha_{K^{-}#Lambda}^{lab} [deg];Counts"},
         {angle_cm.get(), "#Lambda production angle in K^{-}p CM;#theta_{#Lambda}^{CM} [deg];Counts"},
         {costheta_cm.get(), "#Lambda production angle in K^{-}p CM;cos #theta_{#Lambda}^{CM};Counts"}}) {
    if (!item.first) continue;
    canvas.Clear();
    item.first->SetTitle(item.second);
    item.first->SetLineColor(kBlack); item.first->SetLineWidth(2); item.first->Draw("hist");
    canvas.Print(pdf_name);
  }
  for (const auto& item : std::vector<std::pair<TH1*, TString>>{
         {selected_lab.get(), Form("Selected #Lambda production angle (M_{X}>%.2f, %s);#alpha_{K^{-}#Lambda}^{lab} [deg];Candidates", kMissingMassMin, target_selection.Data())},
         {selected_cm.get(), Form("Selected #Lambda production angle (M_{X}>%.2f, %s);#theta_{#Lambda}^{CM} [deg];Candidates", kMissingMassMin, target_selection.Data())},
         {selected_costheta.get(), Form("Selected #Lambda production angle (M_{X}>%.2f, %s);cos #theta_{#Lambda}^{CM};Candidates", kMissingMassMin, target_selection.Data())}}) {
    canvas.Clear(); item.first->SetTitle(item.second); item.first->SetLineColor(kRed + 1); item.first->SetLineWidth(2); item.first->Draw("hist"); canvas.Print(pdf_name);
  }
  for (const auto& item : std::vector<std::pair<TH2D*, TString>>{
         {selected_costheta_lab_vs_lambda_mom_all.get(), "all M_{X}"},
         {selected_costheta_lab_vs_lambda_mom_low.get(), Form("M_{X} #leq %.2f GeV/c^{2}", kMissingMassMin)},
         {selected_costheta_lab_vs_lambda_mom_high.get(), Form("M_{X} > %.2f GeV/c^{2}", kMissingMassMin)}}) {
    canvas.Clear();
    item.first->SetTitle(Form("Selected #Lambda: lab cos#theta vs lab momentum (%s, %s);cos #theta_{#Lambda}^{lab}=p_{z}/|#vec{p}|;|#vec{p}_{#Lambda}|_{lab} [GeV/c]", target_selection.Data(), item.second.Data()));
    item.first->Draw("colz");
    canvas.Print(pdf_name);
  }

  canvas.Clear();
  selected_lambda_mass->SetTitle(
    Form("Selected #Lambda invariant mass (production vertex, %s);M_{p#pi^{-}} [GeV/c^{2}];Candidates / 0.5 MeV/c^{2}", target_selection.Data()));
  selected_lambda_mass->SetLineColor(kRed + 1);
  selected_lambda_mass->SetLineWidth(2);
  selected_lambda_mass->Draw("hist");
  canvas.Print(pdf_name);
  canvas.Clear();
  selected_lambda_mass_high_mm->SetTitle(
    Form("Selected #Lambda invariant mass (production vertex, %s, M_{X}>0.50);M_{p#pi^{-}} [GeV/c^{2}];Candidates / 0.5 MeV/c^{2}", target_selection.Data()));
  selected_lambda_mass_high_mm->SetLineColor(kBlue + 1);
  selected_lambda_mass_high_mm->SetLineWidth(2);
  selected_lambda_mass_high_mm->Draw("hist");
  canvas.Print(pdf_name);

  if (has_acceptance) {
    canvas.Clear();
    acceptance->SetTitle("Trigger acceptance for generated #Lambda#eta;cos #theta_{#Lambda}^{CM};Acceptance");
    acceptance->SetMarkerStyle(20); acceptance->SetMarkerColor(kBlue + 1); acceptance->SetLineColor(kBlue + 1); acceptance->SetMinimum(0.); acceptance->SetMaximum(1.05); acceptance->Draw("E1");
    canvas.Print(pdf_name);
    canvas.Clear();
    selected_costheta_corrected->SetTitle(Form("Background-subtracted #Lambda#eta production angle, trigger-acceptance corrected (%s, 1.09<M_{p#pi^{-}}<1.13);cos #theta_{#Lambda}^{CM};#Lambda#eta yield / Acceptance", target_selection.Data()));
    selected_costheta_corrected->SetLineColor(kMagenta + 1); selected_costheta_corrected->SetLineWidth(2); selected_costheta_corrected->Draw("E1");
    canvas.Print(pdf_name);
    canvas.Clear();
    Double_t selected_yield_error = 0.;
    const Double_t selected_yield = beam_momentum_yield_corrected->IntegralAndError(
      1, kBeamMomentumBins, selected_yield_error);
    beam_momentum_yield_corrected->SetTitle(
      Form("Background-subtracted #Lambda#eta yield vs production-vertex K18 momentum (%s, 1.09<M_{p#pi^{-}}<1.13);|#vec{p}_{K^{-}}^{prod}| [GeV/c];#Lambda#eta yield", target_selection.Data()));
    beam_momentum_yield_corrected->SetMarkerStyle(20);
    beam_momentum_yield_corrected->SetMarkerColor(kMagenta + 1);
    beam_momentum_yield_corrected->SetLineColor(kMagenta + 1);
    beam_momentum_yield_corrected->Draw("E1");
    TPaveText beam_note(.14, .72, .63, .88, "NDC");
    beam_note.SetBorderSize(0); beam_note.SetFillStyle(0); beam_note.SetTextAlign(12); beam_note.SetTextSize(.028);
    beam_note.AddText("4 MeV/c bins; production-vertex momentum");
    beam_note.AddText("M_{X} background-subtracted, then acceptance corrected (735-MeV/c simulation)");
    beam_note.AddText(Form("Selected integral = %.0f #pm %.0f", selected_yield, selected_yield_error));
    if (full_run_scale > 0.) beam_note.AddText(Form("All/selected K-Beam#timesDAQ factor = %.3f", full_run_scale));
    beam_note.Draw();
    canvas.Print(pdf_name);

    if (beam_momentum_yield_expected_all_runs) {
      canvas.Clear();
      Double_t expected_yield_error = 0.;
      const Double_t expected_yield = beam_momentum_yield_expected_all_runs->IntegralAndError(
        1, kBeamMomentumBins, expected_yield_error);
      beam_momentum_yield_expected_all_runs->SetTitle(
        "Expected #Lambda#eta yield for all E72Physics runs (K-Beam scaled);|#vec{p}_{K^{-}}^{prod}| [GeV/c];Expected #Lambda#eta yield");
      beam_momentum_yield_expected_all_runs->SetMarkerStyle(20);
      beam_momentum_yield_expected_all_runs->SetMarkerColor(kBlue + 1);
      beam_momentum_yield_expected_all_runs->SetLineColor(kBlue + 1);
      beam_momentum_yield_expected_all_runs->Draw("E1");
      TPaveText expected_note(.14, .72, .67, .88, "NDC");
      expected_note.SetBorderSize(0); expected_note.SetFillStyle(0); expected_note.SetTextAlign(12); expected_note.SetTextSize(.028);
      expected_note.AddText(Form("Selected K-Beam#timesDAQ = %.3e", kbeam_exposure.selected));
      expected_note.AddText(Form("All E72Physics K-Beam#timesDAQ = %.3e", kbeam_exposure.all_e72physics));
      expected_note.AddText(Form("Expected integral = %.0f #pm %.0f", expected_yield, expected_yield_error));
      expected_note.Draw();
      canvas.Print(pdf_name);
      Info("plot_lambda_missing_mass_vertex",
           "Selected #Lambda#eta yield = %.1f +/- %.1f; all-E72Physics expectation = %.1f +/- %.1f (K-Beam scale %.6f)",
           selected_yield, selected_yield_error, expected_yield, expected_yield_error, full_run_scale);
    }

    // One acceptance-corrected CM production-angle spectrum per 4-MeV/c
    // production-vertex beam-momentum bin; nine panels per PDF page.
    constexpr Int_t kAnglePanelsPerPage = 9;
    for (Int_t first_beam_bin=1; first_beam_bin<=kBeamMomentumBins;
         first_beam_bin+=kAnglePanelsPerPage) {
      TCanvas beam_angle_canvas(Form("c_lambda_angle_by_beam_momentum_%d", first_beam_bin),
                                "Lambda CM production angle by beam momentum", 1200, 900);
      beam_angle_canvas.Divide(3, 3);
      for (Int_t panel=0; panel<kAnglePanelsPerPage; ++panel) {
        const Int_t ibeam = first_beam_bin + panel;
        if (ibeam > kBeamMomentumBins) break;
        beam_angle_canvas.cd(panel+1);
        auto* angle_bin = lambda_eta_yield_by_beam_costheta->ProjectionY(
          Form("h_lambda_eta_costheta_prod_beam_bin_%d", ibeam), ibeam, ibeam);
        angle_bin->SetDirectory(nullptr);
        const Double_t p_low = 1000.*beam_momentum_vs_costheta->GetXaxis()->GetBinLowEdge(ibeam);
        const Double_t p_high = 1000.*beam_momentum_vs_costheta->GetXaxis()->GetBinUpEdge(ibeam);
        angle_bin->SetTitle(Form("p_{K}^{prod}=%.0f--%.0f MeV/c;cos #theta_{#Lambda}^{CM};Background-subtracted, acceptance-corrected #Lambda#eta yield",
                                 p_low, p_high));
        angle_bin->SetMarkerStyle(20);
        angle_bin->SetMarkerSize(.65);
        angle_bin->SetMarkerColor(kMagenta + 1);
        angle_bin->SetLineColor(kMagenta + 1);
        angle_bin->Draw("E1");
      }
      beam_angle_canvas.Print(pdf_name);
    }
  }
  canvas.Print(pdf_name + "]");
  Info("plot_lambda_missing_mass_vertex", "Wrote %s", pdf_name.Data());
}
