// Lambda yield in data (DstTPCLambdaEta) for each cut_flag selection, fitted
// as in yield-fit.cc, to compare with the simulation expectation (simul/yield.cc).
// cut_flag digits: 10000 ntrack (at most one good extra track), 1000 DCA, 100 chisqr,
// 10 kinematics, 1 LH2.
// The vertex ntrack (not in cut_flag) vetoes only good extra tracks that pass
// within extra_track::kVertexTrackDca of the production vertex (see ntrack-class.cc).
//   root -l -b -q -e '.L yield-compare.cc+' -e 'yield_compare()'
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>

#include <fstream>
#include <future>
#include <map>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

#include "extra-track.h"

using namespace std;

namespace {
const vector<int> kRuns = {2447, 2449, 2450, 2451, 2452};
const char *kDataDir = "/home/had/haein/data/JPARC2025Nov_root/physics-735";
const double kXmin = 1.07, kXmax = 1.20;     // fit range [GeV/c2]
const double kCutMin = 1.09, kCutMax = 1.13; // signal count window [GeV/c2]

// A selection lists the cut_flag digits it requires of the one Lambda the DST keeps.
// mm selections also require the missing mass X_mass >= kMissingMassCut (Lambda eta).
// vtxNtrack selections also require the vertex ntrack.
// kin selections also require the momentum part of the DST kinematics cut; the
// cut_flag digit 10 is not used since it also holds the mass window 1.09-1.12,
// which would remove the sidebands of the fit.
struct Selection {
  const char *name, *title;
  vector<int> places;
  bool mm = false;
  bool vtxNtrack = false;
  bool kin = false;
};
const double kMissingMassCut = 0.5; // GeV/c2
// PassSpecialKinematics of DstTPCLambdaEta.cc without the mass window.
const double kKinPionMomMax = 0.25;                         // GeV/c
const double kKinLambdaMomMin = 0.4, kKinLambdaMomMax = 0.6; // GeV/c
const vector<Selection> kSelections = {
    {"none", "no cut", {}},
    {"ntrack", "ntrack", {10000}},
    {"dca", "DCA", {1000}},
    {"chisqr", "#chi^{2}", {100}},
    {"lh2", "LH2", {1}},
    {"dca_chisqr", "DCA, #chi^{2}", {1000, 100}},
    {"ntrack_dca_chisqr", "ntrack, DCA, #chi^{2} (previous selection)", {10000, 1000, 100}},
    {"dca_chisqr_lh2", "DCA, #chi^{2}, LH2", {1000, 100, 1}},
    {"ntrack_dca_chisqr_lh2", "ntrack, DCA, #chi^{2}, LH2", {10000, 1000, 100, 1}},
    {"dca_chisqr_mm", "DCA, #chi^{2}, MM #geq 0.5", {1000, 100}, true},
    {"dca_chisqr_lh2_mm", "DCA, #chi^{2}, LH2, MM #geq 0.5", {1000, 100, 1}, true},
    {"ntrack_dca_chisqr_mm", "ntrack, DCA, #chi^{2}, MM #geq 0.5", {10000, 1000, 100}, true},
    {"ntrack_dca_chisqr_lh2_mm", "ntrack, DCA, #chi^{2}, LH2, MM #geq 0.5", {10000, 1000, 100, 1}, true},
    {"vntrack", "vertex ntrack", {}, false, true},
    {"vntrack_dca_chisqr", "vertex ntrack, DCA, #chi^{2}", {1000, 100}, false, true},
    {"vntrack_dca_chisqr_lh2", "vertex ntrack, DCA, #chi^{2}, LH2", {1000, 100, 1}, false, true},
    {"vntrack_dca_chisqr_mm", "vertex ntrack, DCA, #chi^{2}, MM #geq 0.5", {1000, 100}, true, true},
    {"vntrack_dca_chisqr_lh2_mm", "vertex ntrack, DCA, #chi^{2}, LH2, MM #geq 0.5", {1000, 100, 1}, true, true},
    {"dca_chisqr_kin", "DCA, #chi^{2}, kin.", {1000, 100}, false, false, true},
    {"ntrack_dca_chisqr_kin", "ntrack, DCA, #chi^{2}, kin.", {10000, 1000, 100}, false, false, true},
    {"vntrack_dca_chisqr_kin", "vertex ntrack, DCA, #chi^{2}, kin.", {1000, 100}, false, true, true},
    {"dca_chisqr_lh2_kin", "DCA, #chi^{2}, LH2, kin.", {1000, 100, 1}, false, false, true},
    {"ntrack_dca_chisqr_lh2_kin", "ntrack, DCA, #chi^{2}, LH2, kin.", {10000, 1000, 100, 1}, false, false, true},
    {"vntrack_dca_chisqr_lh2_kin", "vertex ntrack, DCA, #chi^{2}, LH2, kin.", {1000, 100, 1}, false, true, true},
    {"dca_chisqr_lh2_kin_mm", "DCA, #chi^{2}, LH2, kin., MM #geq 0.5", {1000, 100, 1}, true, false, true},
    {"ntrack_dca_chisqr_lh2_kin_mm", "ntrack, DCA, #chi^{2}, LH2, kin., MM #geq 0.5", {10000, 1000, 100, 1}, true,
     false, true},
    {"vntrack_dca_chisqr_lh2_kin_mm", "vertex ntrack, DCA, #chi^{2}, LH2, kin., MM #geq 0.5", {1000, 100, 1}, true,
     true, true},
};
// Missing-mass plot: X_mass of the Lambda candidates in the signal window.
const vector<const char *> kMmPlots = {"dca_chisqr",         "dca_chisqr_lh2",        "ntrack_dca_chisqr",
                                       "ntrack_dca_chisqr_lh2", "vntrack_dca_chisqr", "vntrack_dca_chisqr_lh2",
                                       "dca_chisqr_lh2_kin",    "ntrack_dca_chisqr_lh2_kin",
                                       "vntrack_dca_chisqr_lh2_kin"};
const int kMmBins = 90;
const double kMmMin = -0.1, kMmMax = 0.8; // GeV/c2
// Lambda eta: X_mass of the Lambda window minus the sidebands (same total width),
// fitted with a Gaussian eta peak and a free quadratic times a smeared square-root
// edge at the kinematic endpoint. The endpoint is set by the beam kinematics, so its
// position and smearing come from the first selection and are fixed in the others.
// Each selection is compared with the simulated Lambda eta of the same condition
// (yield-sim.root branch).
// No LH2 cut: the whole target region, as the simulation columns without _lh2.
const vector<pair<const char *, const char *>> kEtaFits = {
    {"dca_chisqr", "tpc"},
    {"ntrack_dca_chisqr", "ntrack"},
    {"vntrack_dca_chisqr", "ntrack0"},
    {"dca_chisqr_kin", "kin"},
    {"ntrack_dca_chisqr_kin", "ntrack_kin"},
    {"vntrack_dca_chisqr_kin", "ntrack0_kin"}};
const double kEtaSigMin = 1.105, kEtaSigMax = 1.125; // Lambda window [GeV/c2]
bool InEtaSideband(double m) { return (m > 1.085 && m < 1.095) || (m > 1.135 && m < 1.145); }
const int kEtaBins = 160;
const double kEtaMin = 0.30, kEtaMax = 0.70;       // GeV/c2, 2.5 MeV bins
const double kEtaFitMin = 0.45, kEtaFitMax = 0.60; // nominal fit range
const vector<double> kEtaFitMinAlt = {0.40, 0.48}; // lower fit edges for the systematic
const double kMassEta = 0.547862;
// chisqr check: p and pi- chisqr of the Lambda candidates passing the DCA cut, in the
// Lambda window (kEtaSigMin-kEtaSigMax) and in its sidebands (same total width), so
// window - sidebands is the Lambda signal and the sidebands its background. "max" is
// the larger of the two, which the cut chisqr <= kChi2Cut on both acts on.
const vector<const char *> kChi2Hists = {"chi2p_win", "chi2p_sb", "chi2pi_win", "chi2pi_sb", "chi2max_win", "chi2max_sb"};
const int kChi2Bins = 200;
const double kChi2Max = 20.; // larger values go to the last bin
const double kChi2Cut = 2.;  // DstTPCLambdaEta.cc kTightChiSqr
// Simulation expectation written by simul/yield.cc.
const char *kSimFile = "result/yield-sim.root";
// K-beam statistics of all physics runs: momentum, run, K-Beam, TRIG-D, DAQ-Eff,
// K-Beam x DAQ-Eff, run list. The yield of all kBeamMomentum runs is these runs
// scaled by the ratio of K-Beam x DAQ-Eff.
const char *kKbeamFile = "/home/had/haein/work/git/e72-physics/statistics/kbeam_statistics.txt";
const int kBeamMomentum = 735; // MeV/c
struct KbeamShare {
  double these = 0, all = 0;           // K-Beam x DAQ-Eff
  double kbeamThese = 0, kbeamAll = 0; // K-Beam
  int nThese = 0, nAll = 0;
};
KbeamShare ReadKbeamShare()
{
  KbeamShare share;
  ifstream in(kKbeamFile);
  string line;
  while (getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    istringstream ss(line);
    int momentum = 0, run = 0;
    double kbeam = 0, trigD = 0, daqEff = 0, kbeamDaq = 0;
    if (!(ss >> momentum >> run >> kbeam >> trigD >> daqEff >> kbeamDaq) || momentum != kBeamMomentum || kbeam <= 0)
      continue;
    share.all += kbeamDaq;
    share.kbeamAll += kbeam;
    ++share.nAll;
    if (find(kRuns.begin(), kRuns.end(), run) != kRuns.end()) {
      share.these += kbeamDaq;
      share.kbeamThese += kbeam;
      ++share.nThese;
    }
  }
  if (share.nAll == 0) cerr << "No " << kBeamMomentum << " MeV/c runs in " << kKbeamFile << endl;
  return share;
}
bool Pass(const Selection &sel, int flag)
{
  for (int place : sel.places)
    if ((flag / place) % 10 != 1) return false;
  return true;
}

// Same model as yield-fit.cc: asymmetric exponential background + Gaussian.
double asymExpBkgOnly(double *x, double *par)
{
  double t = x[0] - par[1];
  double bkg = par[4];
  if (t >= 0) bkg += par[0] * (TMath::Exp(-t / par[3]) - TMath::Exp(-t / par[2]));
  return bkg;
}
double gaussOnly(double *x, double *par) { return par[0] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[1]) / par[2], 2)); }
double asymExpBkgPlusGauss(double *x, double *par) { return asymExpBkgOnly(x, par) + gaussOnly(x, par + 5); }

// Lambda eta model. par: 0 N_eta, 1 mean, 2 sigma, 3-5 quadratic in (MM - 0.5),
// 6 endpoint, 7 endpoint smearing, 8 bin width (fixed).
double EtaEdge(double m, double endpoint, double smear)
{
  double sum = 0, wsum = 0;
  for (int k = -40; k <= 40; ++k) {
    const double d = k * smear / 8.;
    const double w = TMath::Exp(-0.5 * d * d / (smear * smear));
    const double u = endpoint - (m - d);
    sum += w * (u > 0 ? TMath::Sqrt(u) : 0.);
    wsum += w;
  }
  return sum / wsum;
}
double etaBkg(double *x, double *par)
{
  const double u = x[0] - 0.5;
  return (par[3] + par[4] * u + par[5] * u * u) * EtaEdge(x[0], par[6], par[7]);
}
double etaSig(double *x, double *par)
{
  return par[0] * par[8] / (TMath::Sqrt(2 * TMath::Pi()) * par[2]) *
         TMath::Exp(-0.5 * TMath::Power((x[0] - par[1]) / par[2], 2));
}
double etaTot(double *x, double *par) { return etaBkg(x, par) + etaSig(x, par); }
struct EtaFit {
  double n = 0, err = 0, mean = 0, sigma = 0, chi2 = 0;
  int ndf = 0;
  TF1 *tot = nullptr;
};
// edge: fit whose endpoint and smearing are fixed (none: both free).
EtaFit FitEta(TH1D *h, double lo, double hi, const char *name, const EtaFit *edge = nullptr)
{
  auto *f = new TF1(name, etaTot, lo, hi, 9);
  f->SetNpx(500);
  const double peak = h->GetBinContent(h->FindBin(0.555)), cont = h->GetBinContent(h->FindBin(0.50));
  f->SetParameters(std::max(1., (peak - cont) * 4), 0.555, 0.005, std::max(1., cont) * 4, 0, 0, 0.563, 0.002);
  f->FixParameter(8, h->GetBinWidth(1));
  f->SetParLimits(0, 0, 1e6);
  f->SetParLimits(1, 0.540, 0.565);
  f->SetParLimits(2, 0.002, 0.012);
  f->SetParLimits(6, 0.550, 0.580);
  f->SetParLimits(7, 0.001, 0.02);
  if (edge) {
    f->FixParameter(6, edge->tot->GetParameter(6));
    f->FixParameter(7, edge->tot->GetParameter(7));
  }
  h->Fit(f, "RQ0");
  h->Fit(f, "RQ0M");
  return {f->GetParameter(0), f->GetParError(0), f->GetParameter(1), f->GetParameter(2), f->GetChisquare(),
          f->GetNDF(), f};
}

// Name of histogram s of ReadRun: the mass of each selection, the X_mass of each
// kMmPlots selection, then the sideband-subtracted X_mass of each kEtaFits selection.
string HistName(size_t s)
{
  if (s < kSelections.size()) return kSelections[s].name;
  s -= kSelections.size();
  if (s < kMmPlots.size()) return Form("mm_%s", kMmPlots[s]);
  s -= kMmPlots.size();
  if (s < kEtaFits.size()) return Form("eta_%s", kEtaFits[s].first);
  return kChi2Hists[s - kEtaFits.size()];
}

// Mass histograms of one run, one per selection, and the X_mass histograms.
vector<TH1D *> ReadRun(int run)
{
  vector<TH1D *> hists;
  for (const auto &sel : kSelections) {
    auto *h = new TH1D(Form("h_%s_%d", sel.name, run), "", 130, kXmin, kXmax);
    h->SetDirectory(nullptr);
    hists.push_back(h);
  }
  vector<const Selection *> mmSel;
  for (const char *name : kMmPlots) {
    for (const auto &sel : kSelections)
      if (string(sel.name) == name) mmSel.push_back(&sel);
    auto *h = new TH1D(Form("h_%s_%d", HistName(hists.size()).c_str(), run), "", kMmBins, kMmMin, kMmMax);
    h->SetDirectory(nullptr);
    hists.push_back(h);
  }
  vector<const Selection *> etaSel;
  for (const auto &[name, column] : kEtaFits) {
    for (const auto &sel : kSelections)
      if (string(sel.name) == name) etaSel.push_back(&sel);
    auto *h = new TH1D(Form("h_%s_%d", HistName(hists.size()).c_str(), run), "", kEtaBins, kEtaMin, kEtaMax);
    h->SetDirectory(nullptr);
    h->Sumw2();
    hists.push_back(h);
  }
  const size_t chi2Index = hists.size();
  for (size_t j = 0; j < kChi2Hists.size(); ++j) {
    auto *h = new TH1D(Form("h_%s_%d", HistName(hists.size()).c_str(), run), "", kChi2Bins, 0, kChi2Max);
    h->SetDirectory(nullptr);
    hists.push_back(h);
  }
  TFile file(Form("%s/run%05d_DstTPCLambdaEta.root", kDataDir, run));
  TTree *tree = file.IsZombie() ? nullptr : file.Get<TTree>("tpc");
  if (!tree || !tree->GetBranch("cut_flag")) {
    cerr << "No tpc tree with cut_flag for run " << run << endl;
    return hists;
  }
  vector<double> *mass = nullptr, *missingMass = nullptr;
  vector<int> *flag = nullptr;
  Int_t ntTpc = 0;
  vector<int> *isBeam = nullptr, *isK18 = nullptr, *isAccidental = nullptr;
  vector<double> *chisqr = nullptr;
  // Vertex ntrack: Lambda tracks, extra-track quality, helices and production vertex.
  vector<int> *protonId = nullptr, *pionId = nullptr, *prodFound = nullptr, *fitFound = nullptr;
  vector<vector<double>> *hitZ = nullptr, *clX = nullptr, *clZ = nullptr;
  vector<double> *hcx = nullptr, *hcy = nullptr, *hz0 = nullptr, *hr = nullptr, *hdz = nullptr, *htmin = nullptr,
                 *htmax = nullptr;
  vector<double> *prodX = nullptr, *prodY = nullptr, *prodZ = nullptr, *fitX = nullptr, *fitY = nullptr, *fitZ = nullptr;
  // Kinematics: Lambda and pion momenta.
  vector<double> *lamPx = nullptr, *lamPy = nullptr, *lamPz = nullptr, *piPx = nullptr, *piPy = nullptr, *piPz = nullptr;
  tree->SetBranchStatus("*", 0);
  for (const char *b : {"lambda_mass", "X_mass", "cut_flag", "ntTpc", "is_beam", "is_k18", "is_accidental", "chisqr",
                        "lambda_proton_track_id", "lambda_pion_track_id", "hitpos_z", "track_cluster_x_center",
                        "track_cluster_z_center", "helix_cx", "helix_cy", "helix_z0", "helix_r", "helix_dz",
                        "helix_theta_min", "helix_theta_max", "X_prod_found", "X_prod_vtx_x", "X_prod_vtx_y",
                        "X_prod_vtx_z", "X_prod_fit_found", "X_prod_fit_vtx_x", "X_prod_fit_vtx_y", "X_prod_fit_vtx_z",
                        "lambda_mom_x", "lambda_mom_y", "lambda_mom_z", "lambda_pion_mom_x", "lambda_pion_mom_y",
                        "lambda_pion_mom_z"})
    tree->SetBranchStatus(b, 1);
  tree->SetBranchAddress("lambda_proton_track_id", &protonId);
  tree->SetBranchAddress("lambda_pion_track_id", &pionId);
  tree->SetBranchAddress("hitpos_z", &hitZ);
  tree->SetBranchAddress("track_cluster_x_center", &clX);
  tree->SetBranchAddress("track_cluster_z_center", &clZ);
  tree->SetBranchAddress("helix_cx", &hcx);
  tree->SetBranchAddress("helix_cy", &hcy);
  tree->SetBranchAddress("helix_z0", &hz0);
  tree->SetBranchAddress("helix_r", &hr);
  tree->SetBranchAddress("helix_dz", &hdz);
  tree->SetBranchAddress("helix_theta_min", &htmin);
  tree->SetBranchAddress("helix_theta_max", &htmax);
  tree->SetBranchAddress("X_prod_found", &prodFound);
  tree->SetBranchAddress("X_prod_vtx_x", &prodX);
  tree->SetBranchAddress("X_prod_vtx_y", &prodY);
  tree->SetBranchAddress("X_prod_vtx_z", &prodZ);
  tree->SetBranchAddress("X_prod_fit_found", &fitFound);
  tree->SetBranchAddress("X_prod_fit_vtx_x", &fitX);
  tree->SetBranchAddress("X_prod_fit_vtx_y", &fitY);
  tree->SetBranchAddress("X_prod_fit_vtx_z", &fitZ);
  tree->SetBranchAddress("lambda_mom_x", &lamPx);
  tree->SetBranchAddress("lambda_mom_y", &lamPy);
  tree->SetBranchAddress("lambda_mom_z", &lamPz);
  tree->SetBranchAddress("lambda_pion_mom_x", &piPx);
  tree->SetBranchAddress("lambda_pion_mom_y", &piPy);
  tree->SetBranchAddress("lambda_pion_mom_z", &piPz);
  tree->SetBranchAddress("lambda_mass", &mass);
  tree->SetBranchAddress("X_mass", &missingMass);
  tree->SetBranchAddress("cut_flag", &flag);
  tree->SetBranchAddress("ntTpc", &ntTpc);
  tree->SetBranchAddress("is_beam", &isBeam);
  tree->SetBranchAddress("is_k18", &isK18);
  tree->SetBranchAddress("is_accidental", &isAccidental);
  tree->SetBranchAddress("chisqr", &chisqr);
  tree->SetCacheSize(200 * 1024 * 1024);
  tree->AddBranchToCache("*", kTRUE);
  tree->StopCacheLearningPhase();
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    auto beamLike = [&](int it) { return isBeam->at(it) || isK18->at(it) || isAccidental->at(it); };
    // Vertex ntrack of Lambda k: no good extra track near its production vertex
    // (the likelihood fit, else the closest approach of beam and Lambda).
    auto vertexNtrack = [&](size_t k) {
      double vx, vy, vz;
      if (k < fitFound->size() && fitFound->at(k)) {
        vx = fitX->at(k), vy = fitY->at(k), vz = fitZ->at(k);
      } else if (k < prodFound->size() && prodFound->at(k)) {
        vx = prodX->at(k), vy = prodY->at(k), vz = prodZ->at(k);
      } else {
        return false;
      }
      for (int it = 0; it < ntTpc; ++it) {
        if (it == protonId->at(k) || it == pionId->at(k) || beamLike(it)) continue;
        if (hitZ->at(it).size() <= extra_track::kShortHits || extra_track::IsUpstream(clX->at(it), clZ->at(it))) continue;
        const double d = extra_track::HelixDistance(hcx->at(it), hcy->at(it), hz0->at(it), hr->at(it), hdz->at(it),
                                                    htmin->at(it), htmax->at(it), vx, vy, vz);
        if (d < extra_track::kVertexTrackDca) return false;
      }
      return true;
    };
    for (size_t k = 0; k < mass->size() && k < flag->size(); ++k) {
      const bool passMM = k < missingMass->size() && missingMass->at(k) >= kMissingMassCut;
      const bool passVtx = k < protonId->size() && k < pionId->size() && vertexNtrack(k);
      bool passKin = false;
      if (k < lamPx->size() && k < piPx->size()) {
        const double pLam = TMath::Sqrt(lamPx->at(k) * lamPx->at(k) + lamPy->at(k) * lamPy->at(k) + lamPz->at(k) * lamPz->at(k));
        const double pPi = TMath::Sqrt(piPx->at(k) * piPx->at(k) + piPy->at(k) * piPy->at(k) + piPz->at(k) * piPz->at(k));
        passKin = pPi < kKinPionMomMax && pLam >= kKinLambdaMomMin && pLam <= kKinLambdaMomMax;
      }
      auto accept = [&](const Selection &sel) {
        return Pass(sel, flag->at(k)) && (!sel.mm || passMM) && (!sel.vtxNtrack || passVtx) && (!sel.kin || passKin);
      };
      for (size_t s = 0; s < kSelections.size(); ++s)
        if (accept(kSelections[s])) hists[s]->Fill(mass->at(k));
      const double m = mass->at(k);
      // chisqr of the Lambda p and pi- after the DCA cut alone.
      if ((flag->at(k) / 1000) % 10 == 1 && k < protonId->size() && k < pionId->size()) {
        const int side = (m > kEtaSigMin && m < kEtaSigMax) ? 0 : InEtaSideband(m) ? 1 : -1;
        if (side >= 0) {
          const double cp = min(chisqr->at(protonId->at(k)), kChi2Max - 1e-6);
          const double cpi = min(chisqr->at(pionId->at(k)), kChi2Max - 1e-6);
          hists[chi2Index + side]->Fill(cp);
          hists[chi2Index + 2 + side]->Fill(cpi);
          hists[chi2Index + 4 + side]->Fill(max(cp, cpi));
        }
      }
      if (k >= missingMass->size()) continue;
      const double w = (m > kEtaSigMin && m < kEtaSigMax) ? 1. : InEtaSideband(m) ? -1. : 0.;
      if (w != 0.)
        for (size_t j = 0; j < etaSel.size(); ++j)
          if (accept(*etaSel[j])) hists[kSelections.size() + kMmPlots.size() + j]->Fill(missingMass->at(k), w);
      if (m <= kCutMin || m >= kCutMax) continue;
      for (size_t j = 0; j < mmSel.size(); ++j)
        if (accept(*mmSel[j])) hists[kSelections.size() + j]->Fill(missingMass->at(k));
    }
  }
  cout << "run" << run << ": " << tree->GetEntries() << " events" << endl;
  return hists;
}
} // namespace

void yield_compare(const char *outDir = "result")
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);
  // Reading is latency bound on GPFS, so the runs are read in parallel.
  ROOT::EnableThreadSafety();
  vector<future<vector<TH1D *>>> futures;
  for (int run : kRuns) futures.push_back(async(launch::async, ReadRun, run));
  vector<TH1D *> sum;
  for (auto &f : futures) {
    auto hists = f.get();
    if (sum.empty()) {
      for (size_t s = 0; s < hists.size(); ++s)
        sum.push_back((TH1D *)hists[s]->Clone(Form("h_%s", HistName(s).c_str())));
    } else {
      for (size_t s = 0; s < hists.size(); ++s) sum[s]->Add(hists[s]);
    }
  }

  const string runTag = Form("run%05d-%05d-n%zu", kRuns.front(), kRuns.back(), kRuns.size());
  const string pdf = Form("%s/yield-compare-%s.pdf", outDir, runTag.c_str());

  // Fit every selection first; the summary pages go in front of the fits.
  struct FitResult {
    double nSig, nBkg, mean, sigma;
    TF1 *tot, *bkg, *sig;
  };
  vector<FitResult> fits;
  map<string, size_t> index;
  cout << "\nselection" << string(36, ' ') << "  entries   N_sig(1.09-1.13)   N_bkg   mean[GeV]  sigma[GeV]" << endl;
  for (size_t s = 0; s < kSelections.size(); ++s) {
    TH1D *h = sum[s];
    const double binw = h->GetBinWidth(1);
    const double peak = h->GetBinContent(h->FindBin(1.1157));
    const double side = h->GetBinContent(h->FindBin(1.18));
    auto *fTot = new TF1(Form("fTot_%zu", s), asymExpBkgPlusGauss, kXmin, kXmax, 8);
    // Starting values follow yield-fit.cc, scaled to this histogram.
    fTot->SetParameters(2 * peak, 1.0775, 0.0225, 0.045, side, std::max(1.0, peak - side), 1.1157, 0.004);
    fTot->SetParLimits(7, 0.001, 0.02);
    h->Fit(fTot, "RQ0");
    auto *fBkg = new TF1(Form("fBkg_%zu", s), asymExpBkgOnly, kXmin, kXmax, 5);
    auto *fSig = new TF1(Form("fSig_%zu", s), gaussOnly, kXmin, kXmax, 3);
    for (int i = 0; i < 5; ++i) fBkg->SetParameter(i, fTot->GetParameter(i));
    for (int i = 0; i < 3; ++i) fSig->SetParameter(i, fTot->GetParameter(5 + i));
    FitResult r{fSig->Integral(kCutMin, kCutMax) / binw, fBkg->Integral(kCutMin, kCutMax) / binw,
                fTot->GetParameter(6), fTot->GetParameter(7), fTot, fBkg, fSig};
    printf("%-45s %8.0f   %14.0f   %8.0f   %.5f    %.5f\n", kSelections[s].title, h->GetEntries(), r.nSig, r.nBkg,
           r.mean, r.sigma);
    fits.push_back(r);
    index[kSelections[s].name] = s;
  }

  gStyle->SetPaperSize(27, 21); // landscape page matching the 900x700 canvas
  gStyle->SetPadBottomMargin(0.14); // room for the x-axis titles with subscripts
  gStyle->SetPadLeftMargin(0.13);
  TCanvas c("c_yield_compare", "", 900, 700);
  c.Print((pdf + "[").c_str());

  // Page 0: simulation scaled to these runs next to the data.
  {
    TFile simFile(kSimFile);
    TTree *runTree = simFile.IsZombie() ? nullptr : simFile.Get<TTree>("runs");
    TTree *yieldTree = simFile.IsZombie() ? nullptr : simFile.Get<TTree>("yield");
    if (!runTree || !yieldTree) {
      cerr << "No simulation yield in " << kSimFile << " (run simul/yield.cc)" << endl;
    } else {
      Int_t simRun = 0;
      Double_t trigDaq = 0;
      runTree->SetBranchAddress("run", &simRun);
      runTree->SetBranchAddress("trig_daq", &trigDaq);
      double simAll = 0, simUsed = 0;
      int nUsed = 0;
      for (Long64_t i = 0; i < runTree->GetEntries(); ++i) {
        runTree->GetEntry(i);
        simAll += trigDaq;
        if (find(kRuns.begin(), kRuns.end(), simRun) != kRuns.end()) {
          simUsed += trigDaq;
          ++nUsed;
        }
      }
      const double scale = simAll > 0 ? simUsed / simAll : 0;
      string *reaction = nullptr;
      Int_t lambdaPeak = 0;
      // Columns: all, LH2 for no ntrack, ntrack, MM, ntrack + MM (printed per reaction),
      // then ntrack0 and kinematics (only their Lambda peak sums are printed).
      const vector<const char *> simBranches = {"tpc",           "tpc_lh2",         "ntrack",         "ntrack_lh2",
                                                "mm",            "mm_lh2",          "ntrack_mm",      "ntrack_mm_lh2",
                                                "ntrack0",       "ntrack0_lh2",     "ntrack0_mm",     "ntrack0_mm_lh2",
                                                "kin_lh2",       "ntrack_kin_lh2",  "ntrack0_kin_lh2", "mm_kin_lh2",
                                                "ntrack_mm_kin_lh2", "ntrack0_mm_kin_lh2", "kin", "ntrack_kin",
                                                "ntrack0_kin"};
      const size_t nCol = simBranches.size(), nMainCol = 8;
      vector<Double_t> simValue(nCol, 0.);
      yieldTree->SetBranchAddress("reaction", &reaction);
      yieldTree->SetBranchAddress("lambda_peak", &lambdaPeak);
      for (size_t k = 0; k < nCol; ++k) yieldTree->SetBranchAddress(simBranches[k], &simValue[k]);
      const auto *effTpc = simFile.Get<TParameter<double>>("eff_tpc_L");
      const auto *thick = simFile.Get<TParameter<double>>("lh2_thick");

      c.Clear();
      auto *simPave = new TPaveText(0.02, 0.02, 0.98, 0.98, "NDC");
      simPave->SetFillColor(0);
      simPave->SetBorderSize(0);
      simPave->SetTextFont(82);
      simPave->SetTextAlign(12);
      simPave->SetTextSize(0.015);
      simPave->AddText(Form("Data vs simulation, %s", runTag.c_str()));
      simPave->AddText(Form("sim: %s x TRIG-D x DAQ share of these runs (%d) = %.4f", kSimFile, nUsed, scale));
      simPave->AddText(Form("     x eff_tpc_L %.2f, LH2 %.0f mm", effTpc ? effTpc->GetVal() : -1.,
                            thick ? thick->GetVal() : -1.));
      simPave->AddText("sim: physics trigger, p and pi- of the Lambda each hit >= 5 TPC layers;");
      simPave->AddText("     ntrack: at most one other track hits >= 5 TPC layers (primary K- beam aside), as data;");
      simPave->AddText("     ntrack0: no other track; MM: true missing mass (K- p -> Lambda X) >= 0.5 GeV/c2;");
      simPave->AddText("     kin: true p_Lambda 0.4-0.6 and pi- at its first TPC hit < 0.25 GeV/c");
      simPave->AddText(Form("data vertex ntrack: no good extra track within %.0f mm of the production vertex;",
                            extra_track::kVertexTrackDca));
      simPave->AddText("     the simulation has no pile-up, so it is compared with the sim ntrack0");
      simPave->AddText(" ");
      simPave->AddText("                              all             ntrack            MM          ntrack, MM");
      simPave->AddText("Simulation                all     LH2      all     LH2      all     LH2      all     LH2");
      vector<double> peak(nCol, 0.), sigmaP(nCol, 0.), etaSimRow(nCol, 0.);
      auto simRow = [&](const char *label, const vector<double> &v, const char *note) {
        TString line = Form("  %-22s", label);
        for (size_t k = 0; k < nMainCol; ++k) line += Form(" %7.0f%s", v[k], k % 2 ? " " : "");
        simPave->AddText(line + note);
      };
      for (Long64_t i = 0; i < yieldTree->GetEntries(); ++i) {
        yieldTree->GetEntry(i);
        vector<double> v(nCol);
        for (size_t k = 0; k < nCol; ++k) {
          v[k] = simValue[k] * scale;
          (lambdaPeak ? peak : sigmaP)[k] += v[k];
        }
        simRow(reaction->c_str(), v, lambdaPeak ? "" : " (no Lambda)");
        if (*reaction == "Lambda_eta") etaSimRow = v;
      }
      simRow("Lambda peak sum", peak, "");
      simRow("Sigma+ sum (no Lambda)", sigmaP, "");
      for (size_t k = nMainCol; k < nCol; k += 5) {
        TString line = "  Lambda peak:";
        for (size_t j = k; j < min(nCol, k + 5); ++j) line += Form("  %s %.0f", simBranches[j], peak[j]);
        simPave->AddText(line);
      }
      simPave->AddText(" ");
      simPave->AddText("Data N_sig (fit)                                 N_sig   Lambda peak sim          data/sim");
      // Each selection is compared with the simulation column of the same LH2, ntrack, MM and kinematics condition.
      const vector<pair<const char *, const char *>> dataRows = {
          {"none", "tpc"},
          {"dca_chisqr", "tpc"},
          {"dca_chisqr_lh2", "tpc_lh2"},
          {"ntrack_dca_chisqr", "ntrack"},
          {"ntrack_dca_chisqr_lh2", "ntrack_lh2"},
          {"dca_chisqr_mm", "mm"},
          {"dca_chisqr_lh2_mm", "mm_lh2"},
          {"ntrack_dca_chisqr_mm", "ntrack_mm"},
          {"ntrack_dca_chisqr_lh2_mm", "ntrack_mm_lh2"},
          {"vntrack_dca_chisqr", "ntrack0"},
          {"vntrack_dca_chisqr_lh2", "ntrack0_lh2"},
          {"vntrack_dca_chisqr_mm", "ntrack0_mm"},
          {"vntrack_dca_chisqr_lh2_mm", "ntrack0_mm_lh2"},
          {"dca_chisqr_kin", "kin"},
          {"ntrack_dca_chisqr_kin", "ntrack_kin"},
          {"vntrack_dca_chisqr_kin", "ntrack0_kin"},
          {"dca_chisqr_lh2_kin", "kin_lh2"},
          {"ntrack_dca_chisqr_lh2_kin", "ntrack_kin_lh2"},
          {"vntrack_dca_chisqr_lh2_kin", "ntrack0_kin_lh2"},
          {"dca_chisqr_lh2_kin_mm", "mm_kin_lh2"},
          {"ntrack_dca_chisqr_lh2_kin_mm", "ntrack_mm_kin_lh2"},
          {"vntrack_dca_chisqr_lh2_kin_mm", "ntrack0_mm_kin_lh2"}};
      for (const auto &[name, column] : dataRows) {
        const size_t is = index.at(name);
        const string want = column;
        const size_t k = find_if(simBranches.begin(), simBranches.end(), [&](const char *b) { return want == b; }) -
                         simBranches.begin();
        const double n = fits[is].nSig, sim = peak.at(k);
        TString title = kSelections[is].title;
        title.ReplaceAll("#chi^{2}", "chi2").ReplaceAll("p#pi^{-}", "p-pi").ReplaceAll(" (previous selection)", "");
        title.ReplaceAll("#geq", ">=");
        simPave->AddText(Form("  %-46s %8.0f   %7.0f %-20s %6.3f", title.Data(), n, sim, Form("(%s)", column),
                              sim > 0 ? n / sim : 0.));
      }
      simPave->Draw();
      c.Print(pdf.c_str());

      // Missing mass of the Lambda candidates in the signal window, one page per selection.
      for (size_t j = 0; j < kMmPlots.size(); ++j) {
        c.Clear();
        gPad->SetTopMargin(0.08);
        TH1D *h = sum[kSelections.size() + j];
        TString title = kSelections[index.at(kMmPlots[j])].title;
        title.ReplaceAll(" (previous selection)", "");
        h->SetTitle(Form("%s: %s, #it{m}_{p#pi^{-}} %.2f-%.2f;MM [GeV/#it{c}^{2}];counts / 10 MeV/#it{c}^{2}",
                         runTag.c_str(), title.Data(), kCutMin, kCutMax));
        h->SetLineColor(kBlack);
        h->SetLineWidth(2);
        h->SetMaximum(1.1 * h->GetMaximum());
        h->Draw("hist");
        TLine cut;
        cut.SetLineStyle(2);
        cut.SetLineColor(kMagenta + 1);
        cut.DrawLine(kMissingMassCut, 0, kMissingMassCut, h->GetMaximum());
        c.Print(pdf.c_str());
      }

      // Lambda eta: fit of the sideband-subtracted X_mass, one page per selection, then a table.
      struct EtaRow {
        string title, column;
        EtaFit fit;
        double syst, sim;
      };
      vector<EtaRow> etaRows;
      EtaFit edgeFit;           // nominal fit of the first selection
      vector<EtaFit> edgeFitAlt; // and with the other fit ranges
      for (size_t j = 0; j < kEtaFits.size(); ++j) {
        TH1D *h = sum[kSelections.size() + kMmPlots.size() + j];
        const auto &[name, column] = kEtaFits[j];
        EtaFit fit = FitEta(h, kEtaFitMin, kEtaFitMax, Form("fEta_%zu", j), j ? &edgeFit : nullptr);
        // The first fit sets the endpoint; refit with it fixed for a stable error matrix.
        if (j == 0) {
          edgeFit = fit;
          fit = FitEta(h, kEtaFitMin, kEtaFitMax, Form("fEta_%zu_fixed", j), &edgeFit);
        }
        double syst = 0;
        for (size_t a = 0; a < kEtaFitMinAlt.size(); ++a) {
          const EtaFit alt =
              FitEta(h, kEtaFitMinAlt[a], kEtaFitMax, Form("fEta_%zu_%zu", j, a), j ? &edgeFitAlt[a] : nullptr);
          if (j == 0) edgeFitAlt.push_back(alt);
          syst = max(syst, std::abs(alt.n - fit.n));
        }
        const string want = column;
        const size_t k = find_if(simBranches.begin(), simBranches.end(), [&](const char *b) { return want == b; }) -
                         simBranches.begin();
        TString title = kSelections[index.at(name)].title;
        title.ReplaceAll(" (previous selection)", "");
        etaRows.push_back({title.Data(), column, fit, syst, etaSimRow.at(k)});

        c.Clear();
        gPad->SetTopMargin(0.08);
        h->SetTitle(Form("%s: %s, #Lambda %.3f-%.3f minus sidebands;MM [GeV/#it{c}^{2}];counts / %.1f MeV/#it{c}^{2}",
                         runTag.c_str(), title.Data(), kEtaSigMin, kEtaSigMax, 1e3 * h->GetBinWidth(1)));
        h->GetXaxis()->SetRangeUser(0.40, 0.60);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.6);
        h->SetLineColor(kBlack);
        h->SetMinimum(min(0., h->GetMinimum()));
        h->SetMaximum(1.2 * h->GetMaximum());
        h->Draw("e");
        fit.tot->SetLineColor(kRed + 1);
        fit.tot->Draw("same");
        auto *bkg = new TF1(Form("fEtaBkg_%zu", j), etaBkg, kEtaFitMin, kEtaFitMax, 9);
        bkg->SetParameters(fit.tot->GetParameters());
        bkg->SetNpx(500);
        bkg->SetLineColor(kAzure + 2);
        bkg->SetLineStyle(2);
        bkg->Draw("same");
        TLatex t;
        t.SetNDC();
        t.SetTextSize(0.032);
        t.DrawLatex(0.16, 0.86, Form("N_{#Lambda#eta} = %.0f #pm %.0f (stat) #pm %.0f (fit range)", fit.n, fit.err, syst));
        t.DrawLatex(0.16, 0.81, Form("mean %.4f (m_{#eta} %.4f), #sigma %.4f GeV/#it{c}^{2}, #chi^{2}/ndf %.0f/%d", fit.mean,
                                     kMassEta, fit.sigma, fit.chi2, fit.ndf));
        t.DrawLatex(0.16, 0.76, Form("sim #Lambda#eta (%s) %.0f, data/sim %.3f", column, etaRows.back().sim,
                                     etaRows.back().sim > 0 ? fit.n / etaRows.back().sim : 0.));
        c.Print(pdf.c_str());
      }
      c.Clear();
      auto *etaPave = new TPaveText(0.03, 0.03, 0.97, 0.97, "NDC");
      etaPave->SetFillColor(0);
      etaPave->SetBorderSize(0);
      etaPave->SetTextFont(82);
      etaPave->SetTextAlign(12);
      etaPave->SetTextSize(0.0135);
      etaPave->AddText(Form("Lambda eta yield, %s", runTag.c_str()));
      etaPave->AddText(Form("fit of X_mass (Lambda %.3f-%.3f minus sidebands) in %.2f-%.2f GeV/c2:", kEtaSigMin,
                            kEtaSigMax, kEtaFitMin, kEtaFitMax));
      etaPave->AddText("  Gaussian eta + (quadratic) x smeared sqrt edge at the kinematic endpoint");
      etaPave->AddText(Form("  endpoint %.4f, smearing %.4f GeV/c2 from the first selection, fixed in the others",
                            edgeFit.tot->GetParameter(6), edgeFit.tot->GetParameter(7)));
      etaPave->AddText(Form("  systematic: largest change with the lower fit edge at %.2f or %.2f",
                            kEtaFitMinAlt[0], kEtaFitMinAlt[1]));
      etaPave->AddText("  sim: Lambda_eta row of yield-sim.root (true kinematics, x eff_tpc_L), same condition");
      const KbeamShare kbeam = ReadKbeamShare();
      const double all = kbeam.these > 0 ? kbeam.all / kbeam.these : 0.;
      etaPave->AddText(Form("  all %d MeV/c runs: these runs x (K-Beam x DAQ-Eff of all / these) = x %.2f;", kBeamMomentum,
                            all));
      etaPave->AddText("     stat x sqrt of that ratio, syst and sim x the ratio (runs from kbeam_statistics.txt)");
      etaPave->AddText(" ");
      // One table: the statistics, then N_eta (+- stat +- syst) and the simulation for these runs and all runs.
      const TString theseHead = Form("runs %d-%d", kRuns.front(), kRuns.back());
      const TString allHead = Form("all %d MeV/c", kBeamMomentum);
      auto tableRow = [&](const char *label, const TString &these, const TString &allRuns, const TString &simThese,
                          const TString &simAll, const TString &ratio) {
        etaPave->AddText(Form("  %-34s %-20s %-20s %9s %10s %9s", label, these.Data(), allRuns.Data(), simThese.Data(),
                              simAll.Data(), ratio.Data()));
      };
      tableRow("", theseHead, allHead, "sim", "sim all", "data/sim");
      tableRow("runs", Form("%d", kbeam.nThese), Form("%d", kbeam.nAll), "", "", "");
      tableRow("K-Beam", Form("%.4g", kbeam.kbeamThese), Form("%.4g", kbeam.kbeamAll), "", "", "");
      tableRow("K-Beam x DAQ-Eff", Form("%.4g", kbeam.these), Form("%.4g", kbeam.all), "", "", "");
      etaPave->AddText(" ");
      tableRow("N_Lambda_eta +- stat +- syst", "", "", "", "", "");
      for (const auto &r : etaRows) {
        TString title = r.title;
        title.ReplaceAll("#chi^{2}", "chi2");
        tableRow(Form("  %s", title.Data()), Form("%.0f +- %.0f +- %.0f", r.fit.n, r.fit.err, r.syst),
                 Form("%.0f +- %.0f +- %.0f", r.fit.n * all, r.fit.err * TMath::Sqrt(all), r.syst * all),
                 Form("%.0f", r.sim), Form("%.0f", r.sim * all), Form("%.3f", r.sim > 0 ? r.fit.n / r.sim : 0.));
        printf("Lambda eta %-42s %6.0f +- %4.0f +- %4.0f  sim %6.0f  data/sim %.3f  all runs %6.0f +- %4.0f +- %4.0f "
               "(sim %6.0f)\n",
               title.Data(), r.fit.n, r.fit.err, r.syst, r.sim, r.sim > 0 ? r.fit.n / r.sim : 0., r.fit.n * all,
               r.fit.err * TMath::Sqrt(all), r.syst * all, r.sim * all);
      }
      etaPave->AddText(" ");
      etaPave->AddText("  sim column (yield-sim.root): tpc, ntrack, ntrack0, kin, ntrack_kin, ntrack0_kin in this order");
      etaPave->Draw();
      c.Print(pdf.c_str());
      c.Clear();
    }
  }

  // chisqr check: signal (window - sidebands) and background (sidebands) after the DCA cut.
  {
    const size_t i0 = kSelections.size() + kMmPlots.size() + kEtaFits.size();
    const char *what[3] = {"p", "#pi^{-}", "max(p, #pi^{-})"};
    c.Clear();
    c.Divide(3, 1);
    vector<TH1D *> sig(3), bkg(3);
    for (int q = 0; q < 3; ++q) {
      c.cd(q + 1);
      gPad->SetLogy();
      gPad->SetLeftMargin(0.15);
      gPad->SetTopMargin(0.12);
      sig[q] = (TH1D *)sum[i0 + 2 * q]->Clone(Form("h_chi2sig_%d", q));
      sig[q]->Add(sum[i0 + 2 * q + 1], -1);
      bkg[q] = (TH1D *)sum[i0 + 2 * q + 1]->Clone(Form("h_chi2bkg_%d", q));
      for (TH1D *h : {sig[q], bkg[q]})
        if (h->Integral() > 0) h->Scale(1. / h->Integral());
      sig[q]->SetTitle(Form("%s: #chi^{2} of %s, DCA cut;#chi^{2};normalized", runTag.c_str(), what[q]));
      sig[q]->SetLineColor(kRed + 1);
      bkg[q]->SetLineColor(kAzure + 2);
      for (TH1D *h : {sig[q], bkg[q]}) {
        h->SetLineWidth(2);
        h->SetMinimum(1e-5);
      }
      sig[q]->SetMaximum(2 * max(sig[q]->GetMaximum(), bkg[q]->GetMaximum()));
      sig[q]->Draw("hist");
      bkg[q]->Draw("hist same");
      TLine cut;
      cut.SetLineStyle(2);
      cut.SetLineColor(kMagenta + 1);
      cut.DrawLine(kChi2Cut, 1e-5, kChi2Cut, sig[q]->GetMaximum());
      if (q == 0) {
        auto *leg = new TLegend(0.3, 0.74, 0.95, 0.86);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.045);
        leg->AddEntry(sig[q], "#Lambda signal (window - sidebands)", "l");
        leg->AddEntry(bkg[q], "background (sidebands)", "l");
        leg->Draw();
      }
    }
    c.Print(pdf.c_str());

    // Cut scan on both chisqr (max of the two): signal and background kept, S/sqrt(S+B).
    TH1D *win = sum[i0 + 4], *sb = sum[i0 + 5];
    const int n = win->GetNbinsX();
    const double sigAll = win->Integral() - sb->Integral(), bkgAll = sb->Integral();
    auto *gSig = new TGraph(), *gBkg = new TGraph(), *gSignif = new TGraph();
    double signifMax = 0;
    for (int b = 1; b <= n; ++b) {
      const double cutValue = win->GetXaxis()->GetBinUpEdge(b);
      const double S = win->Integral(1, b) - sb->Integral(1, b), B = sb->Integral(1, b);
      gSig->SetPoint(gSig->GetN(), cutValue, sigAll > 0 ? S / sigAll : 0);
      gBkg->SetPoint(gBkg->GetN(), cutValue, bkgAll > 0 ? B / bkgAll : 0);
      const double z = S + B > 0 ? S / TMath::Sqrt(S + B) : 0;
      gSignif->SetPoint(gSignif->GetN(), cutValue, z);
      signifMax = max(signifMax, z);
    }
    c.Clear();
    c.Divide(2, 1);
    c.cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetGridy();
    gSig->SetTitle(Form("%s: kept fraction vs #chi^{2} cut (both p and #pi^{-}), DCA cut;#chi^{2} cut;fraction kept",
                        runTag.c_str()));
    gSig->SetLineColor(kRed + 1);
    gSig->SetLineWidth(2);
    gBkg->SetLineColor(kAzure + 2);
    gBkg->SetLineWidth(2);
    gSig->SetMinimum(0);
    gSig->SetMaximum(1.05);
    gSig->Draw("al");
    gBkg->Draw("l same");
    TLine cut;
    cut.SetLineStyle(2);
    cut.SetLineColor(kMagenta + 1);
    cut.DrawLine(kChi2Cut, 0, kChi2Cut, 1.05);
    auto *leg = new TLegend(0.4, 0.2, 0.88, 0.34);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(gSig, Form("#Lambda signal (all: %.0f)", sigAll), "l");
    leg->AddEntry(gBkg, Form("background (all: %.0f)", bkgAll), "l");
    leg->Draw();
    c.cd(2);
    gPad->SetLeftMargin(0.15);
    gSignif->SetTitle(Form("%s: S/#sqrt{S+B} in #Lambda %.3f-%.3f;#chi^{2} cut;S/#sqrt{S+B}", runTag.c_str(), kEtaSigMin,
                           kEtaSigMax));
    gSignif->SetLineWidth(2);
    gSignif->SetMinimum(0);
    gSignif->SetMaximum(1.1 * signifMax);
    gSignif->Draw("al");
    cut.DrawLine(kChi2Cut, 0, kChi2Cut, 1.1 * signifMax);
    c.Print(pdf.c_str());
    c.Clear();
    for (double cv : {1., 1.5, 2., 3., 5., 10.}) {
      const int b = win->FindBin(cv - 1e-6);
      const double S = win->Integral(1, b) - sb->Integral(1, b), B = sb->Integral(1, b);
      printf("chisqr cut %4.1f: signal kept %.3f, background kept %.3f, S/sqrt(S+B) %.1f\n", cv,
             sigAll > 0 ? S / sigAll : 0, bkgAll > 0 ? B / bkgAll : 0, S + B > 0 ? S / TMath::Sqrt(S + B) : 0);
    }
  }

  // Page 1: each cut alone and the cumulative chain, relative to no cut.
  auto row = [&](TPaveText *pave, const char *label, const char *name, const char *ref) {
    const FitResult &r = fits[index.at(name)], &r0 = fits[index.at(ref)];
    pave->AddText(Form("  %-34s %8.0f   %6.3f   %8.0f   %6.3f", label, r.nSig, r.nSig / r0.nSig, r.nBkg,
                       r.nBkg / r0.nBkg));
  };
  c.Clear();
  auto *pave = new TPaveText(0.03, 0.03, 0.97, 0.97, "NDC");
  pave->SetFillColor(0);
  pave->SetBorderSize(0);
  pave->SetTextFont(82);
  pave->SetTextAlign(12);
  pave->SetTextSize(0.020);
  pave->AddText(Form("Lambda yield vs cut, %s (N_sig, N_bkg in %.2f-%.2f GeV/c2)", runTag.c_str(), kCutMin, kCutMax));
  pave->AddText(" ");
  pave->AddText("Each cut alone, one Lambda per event    N_sig   /nocut     N_bkg   /nocut");
  row(pave, "no cut", "none", "none");
  row(pave, "ntrack", "ntrack", "none");
  row(pave, "vertex ntrack", "vntrack", "none");
  row(pave, "DCA < 5 mm", "dca", "none");
  row(pave, "chi2 <= 2 (p and pi)", "chisqr", "none");
  row(pave, "LH2 production vertex", "lh2", "none");
  pave->AddText(" ");
  pave->AddText("Cumulative (step ratio)                 N_sig   /prev      N_bkg   /prev");
  const vector<pair<const char *, const char *>> chain = {{"no cut", "none"},
                                                           {"+ DCA, chi2", "dca_chisqr"},
                                                           {"+ ntrack", "ntrack_dca_chisqr"},
                                                           {"+ LH2", "ntrack_dca_chisqr_lh2"}};
  for (size_t k = 0; k < chain.size(); ++k)
    row(pave, chain[k].first, chain[k].second, chain[k == 0 ? 0 : k - 1].second);
  const FitResult &first = fits[index.at(chain.front().second)], &last = fits[index.at(chain.back().second)];
  pave->AddText(Form("  %-34s %8s   %6.3f   %8s   %6.3f", "total (last / first)", "", last.nSig / first.nSig, "",
                     last.nBkg / first.nBkg));
  pave->AddText(" ");
  pave->AddText("Vertex ntrack instead of ntrack (step ratio)");
  row(pave, "+ DCA, chi2", "dca_chisqr", "none");
  row(pave, "+ vertex ntrack", "vntrack_dca_chisqr", "dca_chisqr");
  row(pave, "+ LH2", "vntrack_dca_chisqr_lh2", "vntrack_dca_chisqr");
  row(pave, "vertex ntrack / ntrack", "vntrack_dca_chisqr", "ntrack_dca_chisqr");
  row(pave, "vertex ntrack / ntrack, LH2", "vntrack_dca_chisqr_lh2", "ntrack_dca_chisqr_lh2");
  row(pave, "vertex ntrack / ntrack, LH2, MM", "vntrack_dca_chisqr_lh2_mm", "ntrack_dca_chisqr_lh2_mm");
  pave->AddText(" ");
  pave->AddText(Form("Kinematics: p_pi < %.2f, %.1f <= p_Lambda <= %.1f GeV/c (step ratio)", kKinPionMomMax,
                     kKinLambdaMomMin, kKinLambdaMomMax));
  row(pave, "DCA, chi2, LH2 + kin.", "dca_chisqr_lh2_kin", "dca_chisqr_lh2");
  row(pave, "ntrack, DCA, chi2, LH2 + kin.", "ntrack_dca_chisqr_lh2_kin", "ntrack_dca_chisqr_lh2");
  row(pave, "vtx ntrack, DCA, chi2, LH2 + kin.", "vntrack_dca_chisqr_lh2_kin", "vntrack_dca_chisqr_lh2");
  row(pave, "DCA, chi2, LH2, MM + kin.", "dca_chisqr_lh2_kin_mm", "dca_chisqr_lh2_mm");
  row(pave, "ntrack, DCA, chi2, LH2, MM + kin.", "ntrack_dca_chisqr_lh2_kin_mm", "ntrack_dca_chisqr_lh2_mm");
  row(pave, "vtx ntrack, ..., LH2, MM + kin.", "vntrack_dca_chisqr_lh2_kin_mm", "vntrack_dca_chisqr_lh2_mm");
  pave->AddText("(the cut_flag kinematics digit also holds the mass window 1.09-1.12; not used)");
  pave->Draw();
  c.Print(pdf.c_str());

  // Page 2: the cumulative chain as a bar chart.
  c.Clear();
  auto *hChain = new TH1D("h_chain", Form("%s: #Lambda yield, cumulative cuts;;N_{sig}", runTag.c_str()),
                          chain.size(), 0, chain.size());
  hChain->SetDirectory(nullptr);
  for (size_t k = 0; k < chain.size(); ++k) {
    hChain->SetBinContent(k + 1, fits[index.at(chain[k].second)].nSig);
    hChain->GetXaxis()->SetBinLabel(k + 1, chain[k].first);
  }
  hChain->SetFillColor(kAzure - 9);
  hChain->SetBarWidth(0.7);
  hChain->SetBarOffset(0.15);
  hChain->SetMinimum(0);
  gPad->SetBottomMargin(0.15);
  hChain->Draw("bar");
  TLatex barText;
  barText.SetTextAlign(21);
  barText.SetTextSize(0.03);
  for (size_t k = 0; k < chain.size(); ++k) {
    const double n = hChain->GetBinContent(k + 1);
    barText.DrawLatex(k + 0.5, n * 1.02, Form("%.0f (%.0f%%)", n, 100. * n / hChain->GetBinContent(1)));
  }
  c.Print(pdf.c_str());
  gPad->SetBottomMargin(0.14);

  for (size_t s = 0; s < kSelections.size(); ++s) {
    TH1D *h = sum[s];
    const FitResult &r = fits[s];
    c.Clear();
    h->SetTitle(Form("%s: %s;#it{m}_{p#pi^{-}} [GeV/#it{c}^{2}];Counts/(0.001 GeV/#it{c}^{2})", runTag.c_str(),
                     kSelections[s].title));
    h->SetFillColorAlpha(kBlue, 0.1);
    h->Draw();
    r.tot->SetLineColor(kGreen + 2);
    r.bkg->SetLineColor(kBlack);
    r.bkg->SetLineStyle(2);
    r.sig->SetLineColor(kRed);
    for (TF1 *f : {r.tot, r.bkg, r.sig}) {
      f->SetLineWidth(3);
      f->Draw("same");
    }
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.035);
    lat.DrawLatex(0.55, 0.54, Form("Yield = %.0f", r.nSig));
    lat.DrawLatex(0.55, 0.49, Form("#mu = %.5f", r.mean));
    lat.DrawLatex(0.55, 0.44, Form("#sigma = %.5f", r.sigma));
    lat.DrawLatex(0.55, 0.39, Form("Bkg = %.0f", r.nBkg));
    c.Print(pdf.c_str());
  }
  c.Print((pdf + "]").c_str());
  cout << "Wrote " << pdf << endl;
}
