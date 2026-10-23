// Event classes of the extra tracks of Lambda eta candidates (DCA, chisqr and
// X_mass >= 0.5 GeV/c2), to see what the ntrack cut of DstTPCLambdaEta.cc
// removes. Extra tracks are those that are not beam-like and not the p/pi of
// the Lambda; "good" extra tracks are neither short (<= 6 hits) nor upstream
// beam-window tracks. The DCA of an extra track is the closest distance of its
// helix to the production vertex (X_prod_fit_vtx, or X_prod_vtx when the fit
// failed); DCAs beyond the histogram range are drawn in the last bin.
// Every histogram is filled with +1 in the Lambda window and
// -1 in the sidebands of the same total width, so it is sideband subtracted.
//   root -l -b -q -e '.L ntrack-class.cc+' -e 'ntrack_class()'        // LH2
//   root -l -b -q -e '.L ntrack-class.cc+' -e 'ntrack_class("result", false)' // all vertices
#include <TCanvas.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <future>
#include <iostream>
#include <vector>

#include "extra-track.h"

using namespace std;
using extra_track::HelixDistance;
using extra_track::IsUpstream;

namespace {
const vector<int> kRuns = {2447, 2449, 2450, 2451, 2452};
const char *kDataDir = "/home/had/haein/data/JPARC2025Nov_root/physics-735";
const double kSigMin = 1.105, kSigMax = 1.125;
bool InSideband(double m) { return (m > 1.085 && m < 1.095) || (m > 1.135 && m < 1.145); }
const double kMissingMassCut = 0.5;   // GeV/c2
const double kVertexDca = 10.;        // mm, extra track to the production vertex
const double kFarDca = extra_track::kVertexTrackDca; // mm, boundary of the far group in the Lambda MM split
const double kDcaMax = 60.;           // mm, histogram range; larger DCAs go to the last bin
// Lambda MM regions to compare the far single extra tracks: continuum and Lambda eta peak.
const double kContMin = 0.50, kContMax = 0.54, kPeakMin = 0.545, kPeakMax = 0.565;
const double kExclusiveMM2 = 0.04;    // GeV2/c4, |MM2(K- p -> Lambda pi+ pi- X)|
constexpr double kZTarget = -143.;
const double kMassK = 0.493677, kMassP = 0.93827208, kMassPi = 0.13957039, kMassLambda = 1.115683;

const vector<const char *> kClasses = {"no extra track",
                                       "only short/upstream extra",
                                       "1 good extra",
                                       "2 good, +-, MM^{2}(#Lambda#pi^{+}#pi^{-}) #approx 0",
                                       "2 good, +-, MM^{2} #neq 0",
                                       "2 good, same charge",
                                       "3+ good extra"};
const vector<const char *> kSpecies = {"#pi^{-}", "#pi^{+}", "p", "K^{#pm}", "e/other"};
int Species(int pid, int charge)
{
  if (pid & 0x4 && charge > 0) return 2;
  if (pid & 0x1) return charge < 0 ? 0 : 1;
  if (pid & 0x2) return 3;
  return 4;
}

struct Hists {
  TH1D *cls, *nAll, *nGood, *mm2, *mm2Same, *lambdaMM[7];
  TH1D *oneDca[5], *oneSpecies[2], *oneMom[2], *oneVtz[2], *oneMM[3];
  TH1D *twoDca, *threeN, *threeSpecies, *goodDca;
  // Far (DCA >= kFarDca) single good extra track, [0] continuum and [1] peak of the Lambda MM.
  TH1D *farSpecies[2], *farCharge[2], *farMom[2], *farInner[2], *farPileup[2], *farHits[2];
  TH2D *farDedx[2];
};
Hists MakeHists(int run)
{
  Hists h;
  auto H = [&](const char *name, const char *title, int n, double lo, double hi) {
    auto *x = new TH1D(Form("%s_%d", name, run), title, n, lo, hi);
    x->SetDirectory(nullptr);
    x->Sumw2();
    return x;
  };
  h.cls = H("cls", ";;#Lambda events", kClasses.size(), 0, kClasses.size());
  h.nAll = H("nAll", ";extra tracks (all);#Lambda events", 8, -0.5, 7.5);
  h.nGood = H("nGood", ";good extra tracks;#Lambda events", 8, -0.5, 7.5);
  h.mm2 = H("mm2", ";MM^{2}(K^{-}p#rightarrow#Lambda#pi^{+}#pi^{-}X) [GeV^{2}/#it{c}^{4}];#Lambda events / 0.005", 100, -0.2, 0.3);
  h.mm2Same = H("mm2Same", ";MM^{2} [GeV^{2}/#it{c}^{4}];#Lambda events / 0.005", 100, -0.2, 0.3);
  for (int c = 0; c < 7; ++c)
    h.lambdaMM[c] = H(Form("lambdaMM%d", c), ";MM(K^{-}p#rightarrow#LambdaX) [GeV/#it{c}^{2}];#Lambda events / 5 MeV", 14, 0.50, 0.57);
  for (int s = 0; s < 5; ++s)
    h.oneDca[s] = H(Form("oneDca%d", s), ";DCA to the production vertex [mm];#Lambda events / 1 mm", 60, 0, kDcaMax);
  for (int v = 0; v < 3; ++v)
    h.oneMM[v] = H(Form("oneMM%d", v), ";MM(K^{-}p#rightarrow#LambdaX) [GeV/#it{c}^{2}];#Lambda events / 5 MeV", 14, 0.50, 0.57);
  for (int v = 0; v < 2; ++v) {
    h.oneSpecies[v] = H(Form("oneSpecies%d", v), ";;#Lambda events", 5, 0, 5);
    h.oneMom[v] = H(Form("oneMom%d", v), ";extra track p [GeV/#it{c}];#Lambda events / 20 MeV", 50, 0, 1);
    h.oneVtz[v] = H(Form("oneVtz%d", v), ";production vertex z - Z_{T} [mm];#Lambda events / 10 mm", 60, -300, 300);
  }
  h.twoDca = H("twoDca", ";DCA to the production vertex [mm];tracks / 1 mm", 60, 0, kDcaMax);
  h.goodDca = H("goodDca", ";DCA to the production vertex [mm];tracks / 1 mm", 60, 0, kDcaMax);
  h.threeN = H("threeN", ";good extra tracks from the vertex (DCA < 10 mm);#Lambda events", 8, -0.5, 7.5);
  h.threeSpecies = H("threeSpecies", ";;tracks", 5, 0, 5);
  for (int r = 0; r < 2; ++r) {
    h.farSpecies[r] = H(Form("farSpecies%d", r), ";;fraction", 5, 0, 5);
    h.farCharge[r] = H(Form("farCharge%d", r), ";charge;fraction", 3, -1.5, 1.5);
    h.farMom[r] = H(Form("farMom%d", r), ";extra track p [GeV/#it{c}];fraction / 50 MeV", 20, 0, 1);
    h.farInner[r] = H(Form("farInner%d", r), ";innermost hit to the production vertex [mm];fraction / 10 mm", 30, 0, 300);
    h.farPileup[r] = H(Form("farPileup%d", r), ";DCA to another beam-like track [mm];fraction / 2 mm", 20, 0, 40);
    h.farHits[r] = H(Form("farHits%d", r), ";hits on the track;fraction / 2", 20, 0, 40);
    h.farDedx[r] = new TH2D(Form("farDedx%d_%d", r, run), ";charge #times p [GeV/#it{c}];dE/dx [arb.]", 50, -1, 1, 50, 0, 250);
    h.farDedx[r]->SetDirectory(nullptr);
  }
  return h;
}
vector<TH1 *> All(const Hists &h)
{
  vector<TH1 *> v = {h.cls, h.nAll, h.nGood, h.mm2, h.mm2Same, h.twoDca, h.goodDca, h.threeN, h.threeSpecies};
  for (auto *x : h.lambdaMM) v.push_back(x);
  for (auto *x : h.oneDca) v.push_back(x);
  for (auto *x : h.oneMM) v.push_back(x);
  for (int k = 0; k < 2; ++k) {
    v.push_back(h.oneSpecies[k]);
    v.push_back(h.oneMom[k]);
    v.push_back(h.oneVtz[k]);
    for (TH1 *x : {(TH1 *)h.farSpecies[k], (TH1 *)h.farCharge[k], (TH1 *)h.farMom[k], (TH1 *)h.farInner[k],
                   (TH1 *)h.farPileup[k], (TH1 *)h.farHits[k], (TH1 *)h.farDedx[k]})
      v.push_back(x);
  }
  return v;
}

Hists ReadRun(int run, bool lh2)
{
  Hists h = MakeHists(run);
  TFile file(Form("%s/run%05d_DstTPCLambdaEta.root", kDataDir, run));
  TTree *tree = file.IsZombie() ? nullptr : file.Get<TTree>("tpc");
  if (!tree || !tree->GetBranch("cut_flag")) {
    cerr << "No tpc tree with cut_flag for run " << run << endl;
    return h;
  }
  vector<double> *mass = nullptr, *missingMass = nullptr, *lx = nullptr, *ly = nullptr, *lz = nullptr;
  vector<double> *bx = nullptr, *by = nullptr, *bz = nullptr, *mom0 = nullptr;
  vector<int> *flag = nullptr, *protonId = nullptr, *pionId = nullptr;
  vector<double> *dedx = nullptr;
  vector<vector<double>> *closeDist = nullptr, *hitX = nullptr, *hitY = nullptr;
  vector<int> *isBeam = nullptr, *isK18 = nullptr, *isAccidental = nullptr, *pid = nullptr, *charge = nullptr;
  vector<vector<double>> *momX = nullptr, *momY = nullptr, *momZ = nullptr;
  vector<vector<double>> *hitZ = nullptr, *clX = nullptr, *clZ = nullptr;
  vector<double> *hcx = nullptr, *hcy = nullptr, *hz0 = nullptr, *hr = nullptr, *hdz = nullptr, *htmin = nullptr,
                 *htmax = nullptr;
  vector<int> *prodFound = nullptr, *fitFound = nullptr;
  vector<double> *prodX = nullptr, *prodY = nullptr, *prodZ = nullptr, *fitX = nullptr, *fitY = nullptr, *fitZ = nullptr;
  Int_t ntTpc = 0;
  tree->SetBranchStatus("*", 0);
  for (const char *b : {"lambda_mass", "X_mass", "lambda_mom_x", "lambda_mom_y", "lambda_mom_z", "X_prod_beam_px",
                        "X_prod_beam_py", "X_prod_beam_pz", "cut_flag", "lambda_proton_track_id",
                        "lambda_pion_track_id", "is_beam", "is_k18", "is_accidental", "pid", "charge", "mom0", "dEdx", "closeDistTpc", "hitpos_x", "hitpos_y",
                        "mom_vtx", "mom_vty", "mom_vtz", "ntTpc", "hitpos_z", "helix_cx", "helix_cy", "helix_z0", "helix_r",
                        "helix_dz", "helix_theta_min", "helix_theta_max", "X_prod_found", "X_prod_vtx_x",
                        "X_prod_vtx_y", "X_prod_vtx_z", "X_prod_fit_found", "X_prod_fit_vtx_x",
                        "X_prod_fit_vtx_y", "X_prod_fit_vtx_z",
                        "track_cluster_x_center", "track_cluster_z_center"})
    tree->SetBranchStatus(b, 1);
  tree->SetBranchAddress("lambda_mass", &mass);
  tree->SetBranchAddress("X_mass", &missingMass);
  tree->SetBranchAddress("lambda_mom_x", &lx);
  tree->SetBranchAddress("lambda_mom_y", &ly);
  tree->SetBranchAddress("lambda_mom_z", &lz);
  tree->SetBranchAddress("X_prod_beam_px", &bx);
  tree->SetBranchAddress("X_prod_beam_py", &by);
  tree->SetBranchAddress("X_prod_beam_pz", &bz);
  tree->SetBranchAddress("cut_flag", &flag);
  tree->SetBranchAddress("lambda_proton_track_id", &protonId);
  tree->SetBranchAddress("lambda_pion_track_id", &pionId);
  tree->SetBranchAddress("is_beam", &isBeam);
  tree->SetBranchAddress("is_k18", &isK18);
  tree->SetBranchAddress("is_accidental", &isAccidental);
  tree->SetBranchAddress("pid", &pid);
  tree->SetBranchAddress("charge", &charge);
  tree->SetBranchAddress("mom0", &mom0);
  tree->SetBranchAddress("dEdx", &dedx);
  tree->SetBranchAddress("closeDistTpc", &closeDist);
  tree->SetBranchAddress("hitpos_x", &hitX);
  tree->SetBranchAddress("hitpos_y", &hitY);
  tree->SetBranchAddress("mom_vtx", &momX);
  tree->SetBranchAddress("mom_vty", &momY);
  tree->SetBranchAddress("mom_vtz", &momZ);
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
  tree->SetBranchAddress("ntTpc", &ntTpc);
  tree->SetBranchAddress("hitpos_z", &hitZ);
  tree->SetBranchAddress("track_cluster_x_center", &clX);
  tree->SetBranchAddress("track_cluster_z_center", &clZ);
  tree->SetCacheSize(200 * 1024 * 1024);
  tree->AddBranchToCache("*", kTRUE);
  tree->StopCacheLearningPhase();
  long nNoVertex = 0;
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    if (flag->empty() || mass->empty() || protonId->empty() || pionId->empty() || missingMass->empty()) continue;
    const int f = flag->at(0);
    if ((f / 1000) % 10 != 1 || (f / 100) % 10 != 1) continue; // DCA and chisqr
    if (lh2 && f % 10 != 1) continue;
    if (!(missingMass->at(0) >= kMissingMassCut)) continue;
    const double m = mass->at(0);
    const double w = (m > kSigMin && m < kSigMax) ? 1. : InSideband(m) ? -1. : 0.;
    if (w == 0.) continue;
    const int ip = protonId->at(0), ipi = pionId->at(0);
    auto beamLike = [&](int it) { return isBeam->at(it) || isK18->at(it) || isAccidental->at(it); };
    // Production vertex: the likelihood fit, else the closest approach of beam and Lambda.
    double vx, vy, vz;
    if (!fitFound->empty() && fitFound->at(0)) {
      vx = fitX->at(0), vy = fitY->at(0), vz = fitZ->at(0);
    } else if (!prodFound->empty() && prodFound->at(0)) {
      vx = prodX->at(0), vy = prodY->at(0), vz = prodZ->at(0);
    } else {
      ++nNoVertex;
      continue;
    }
    auto dcaVertex = [&](int it) {
      const double d = HelixDistance(hcx->at(it), hcy->at(it), hz0->at(it), hr->at(it), hdz->at(it), htmin->at(it),
                                     htmax->at(it), vx, vy, vz);
      return std::isfinite(d) ? min(d, kDcaMax - 0.5) : kDcaMax - 0.5;
    };
    vector<int> good;
    int nExtra = 0;
    for (int it = 0; it < ntTpc; ++it) {
      if (it == ip || it == ipi || beamLike(it)) continue;
      ++nExtra;
      if (hitZ->at(it).size() <= extra_track::kShortHits || IsUpstream(clX->at(it), clZ->at(it))) continue;
      good.push_back(it);
      h.goodDca->Fill(dcaVertex(it), w);
    }
    h.nAll->Fill(nExtra, w);
    h.nGood->Fill(good.size(), w);

    // Missing mass squared of K- p -> Lambda pi+ pi- X with two good extra tracks.
    auto mm2 = [&](int a, int b) {
      TLorentzVector A, B, K, L;
      A.SetXYZM(momX->at(a).at(b), momY->at(a).at(b), momZ->at(a).at(b), kMassPi);
      B.SetXYZM(momX->at(b).at(a), momY->at(b).at(a), momZ->at(b).at(a), kMassPi);
      K.SetXYZM(bx->at(0), by->at(0), bz->at(0), kMassK);
      L.SetXYZM(lx->at(0), ly->at(0), lz->at(0), kMassLambda);
      return (K + TLorentzVector(0, 0, 0, kMassP) - L - A - B).M2();
    };
    int c;
    if (nExtra == 0) c = 0;
    else if (good.empty()) c = 1;
    else if (good.size() == 1) c = 2;
    else if (good.size() == 2 && charge->at(good[0]) * charge->at(good[1]) < 0) {
      const double x = mm2(good[0], good[1]);
      h.mm2->Fill(x, w);
      c = std::abs(x) < kExclusiveMM2 ? 3 : 4;
    } else if (good.size() == 2) c = 5;
    else c = 6;
    h.cls->Fill(c, w);
    h.lambdaMM[c]->Fill(missingMass->at(0), w);

    if (c == 2) {
      const int k = good[0];
      const int v = dcaVertex(k) < kVertexDca;
      const int s = Species(pid->at(k), charge->at(k));
      h.oneDca[s]->Fill(dcaVertex(k), w);
      h.oneMM[dcaVertex(k) < kVertexDca ? 0 : dcaVertex(k) < kFarDca ? 1 : 2]->Fill(missingMass->at(0), w);
      h.oneSpecies[v]->Fill(s, w);
      h.oneMom[v]->Fill(mom0->at(k), w);
      h.oneVtz[v]->Fill(vz - kZTarget, w);
      const double mm = missingMass->at(0);
      const int r = (mm > kContMin && mm < kContMax) ? 0 : (mm > kPeakMin && mm < kPeakMax) ? 1 : -1;
      if (dcaVertex(k) >= kFarDca && r >= 0) {
        h.farSpecies[r]->Fill(s, w);
        h.farCharge[r]->Fill(charge->at(k), w);
        h.farMom[r]->Fill(mom0->at(k), w);
        h.farHits[r]->Fill(min<double>(hitZ->at(k).size(), 39.5), w);
        h.farDedx[r]->Fill(charge->at(k) * mom0->at(k), dedx->at(k), w);
        double inner = 1e9;
        for (size_t j = 0; j < hitZ->at(k).size(); ++j)
          inner = min(inner, std::hypot(hitX->at(k)[j] - vx, hitY->at(k)[j] - vy, hitZ->at(k)[j] - vz));
        h.farInner[r]->Fill(min(inner, 295.), w);
        // Closest beam-like track other than the K18 one: a vertex of a pile-up beam (none or >= 30 mm: last bin).
        double pile = 1e9;
        for (int j = 0; j < ntTpc; ++j)
          if (j != k && beamLike(j) && !isK18->at(j) && std::isfinite(closeDist->at(k).at(j))) pile = min(pile, closeDist->at(k).at(j));
        h.farPileup[r]->Fill(min(pile, 39.), w);
      }
    }
    if (c == 5)
      for (int k : good) h.twoDca->Fill(dcaVertex(k), w);
    if (c == 6) {
      int nv = 0;
      for (int k : good) {
        nv += dcaVertex(k) < kVertexDca;
        h.threeSpecies->Fill(Species(pid->at(k), charge->at(k)), w);
      }
      h.threeN->Fill(nv, w);
    }
  }
  cout << "run" << run << " done, " << nNoVertex << " sideband-weighted entries without a production vertex" << endl;
  return h;
}
} // namespace

void ntrack_class(const char *outDir = "result", bool lh2 = true)
{
  ROOT::EnableThreadSafety();
  gStyle->SetOptStat(0);
  gStyle->SetPaperSize(27, 18); // landscape; larger sizes get clipped at the left edge of the PDF
  gSystem->mkdir(outDir, kTRUE);
  vector<future<Hists>> jobs;
  for (int run : kRuns) jobs.push_back(async(launch::async, ReadRun, run, lh2));
  Hists sum = MakeHists(0);
  vector<TH1 *> total = All(sum);
  for (auto &job : jobs) {
    vector<TH1 *> part = All(job.get());
    for (size_t i = 0; i < total.size(); ++i) total[i]->Add(part[i]);
  }

  const TString sel = Form("DCA, #chi^{2}%s, MM #geq 0.5, #Lambda %.3f-%.3f minus sidebands", lh2 ? ", LH2" : "",
                           kSigMin, kSigMax);
  const TString runTag = Form("run%05d-%05d-n%zu", kRuns.front(), kRuns.back(), kRuns.size());
  const TString pdf = Form("%s/ntrack-class-%s%s.pdf", outDir, runTag.Data(), lh2 ? "-lh2" : "");
  TCanvas canvas("c_ntrack_class", "", 1800, 1200);
  canvas.Print(pdf + "[");
  gStyle->SetTitleFontSize(0.05);
  TGaxis::SetMaxDigits(4);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleX(0.5);
  // Every page: a header line and a body pad divided into nx x ny pads.
  TPad *body = nullptr;
  auto page = [&](int nx, int ny, const char *what) {
    canvas.Clear();
    canvas.cd();
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.02);
    t.SetTextAlign(13);
    t.DrawLatex(0.01, 0.985, Form("%s: %s  |  %s", runTag.Data(), sel.Data(), what));
    body = new TPad("body", "", 0, 0, 1, 0.965);
    body->Draw();
    body->Divide(nx, ny, 0.002, 0.002);
  };
  // Go to pad i; a wide right margin keeps the legend outside the frame.
  auto pad = [&](int i, bool legendRoom) {
    body->cd(i);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.09);
    gPad->SetRightMargin(legendRoom ? 0.28 : 0.04);
  };
  auto legend = [&](int n) {
    const double x = 1 - gPad->GetRightMargin() + 0.01, y = 1 - gPad->GetTopMargin();
    auto *l = new TLegend(x, y - 0.06 * n, 0.995, y);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextSize(0.035);
    return l;
  };
  auto style = [](TH1 *h, int color, const char *title) {
    h->SetTitle(title);
    h->SetLineColor(color);
    h->SetLineWidth(2);
    h->SetFillStyle(0);
    for (auto *a : {h->GetXaxis(), h->GetYaxis()}) {
      a->SetTitleSize(0.05);
      a->SetLabelSize(0.045);
    }
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetNdivisions(507);
    h->GetYaxis()->SetTitleOffset(1.45);
  };
  // Draw histograms on one frame with a common maximum.
  auto overlay = [](const vector<TH1 *> &hs) {
    double hi = 0, lo = 0;
    for (auto *h : hs) {
      hi = max(hi, h->GetMaximum());
      lo = min(lo, h->GetMinimum());
    }
    for (size_t i = 0; i < hs.size(); ++i) {
      hs[i]->SetMaximum(1.1 * hi);
      hs[i]->SetMinimum(lo < 0 ? 1.1 * lo : 0);
      hs[i]->Draw(i ? "hist same" : "hist");
    }
  };
  TLine line;
  line.SetLineStyle(2);
  line.SetLineColor(kMagenta + 1);
  // DCA plots: mark the 10 mm vertex cut, the 30 mm boundary and the overflow bin.
  auto dcaNote = [&](TH1 *h) {
    const double ymax = h->GetMaximum();
    line.DrawLine(kVertexDca, 0, kVertexDca, ymax);
    TLine l;
    l.SetLineStyle(3);
    l.DrawLine(kFarDca, 0, kFarDca, ymax);
    TLatex t;
    t.SetTextSize(0.035);
    t.SetTextAngle(90);
    t.DrawLatex(kDcaMax - 1.5, 0.15 * ymax, Form("#geq %.0f mm", kDcaMax - 1));
  };
  auto speciesLabels = [](TH1 *h) {
    for (int s = 0; s < 5; ++s) h->GetXaxis()->SetBinLabel(s + 1, kSpecies[s]);
    h->GetXaxis()->SetLabelSize(0.07);
  };
  const double nLambda = sum.cls->Integral();
  const vector<const char *> kShort = {"no extra", "short/upstream only", "1 good", "2 good #pm, #Lambda#pi#pi",
                                       "2 good #pm, other", "2 good, same q", "3+ good"};
  const int colors[7] = {kBlack, kGray + 1, kAzure + 2, kRed + 1, kOrange + 1, kGreen + 2, kViolet};

  // Page 1: event classes.
  page(2, 1, "current ntrack keeps only \"no extra\"");
  pad(1, false);
  gPad->SetLeftMargin(0.45);
  gPad->SetRightMargin(0.06);
  auto *hCls = (TH1D *)sum.cls->Clone("hCls");
  style(hCls, kAzure + 2, Form("event classes, N_{#Lambda} = %.0f;;#Lambda events", nLambda));
  for (size_t k = 0; k < kClasses.size(); ++k) hCls->GetXaxis()->SetBinLabel(k + 1, kClasses[k]);
  hCls->GetXaxis()->SetLabelSize(0.04);
  hCls->GetYaxis()->SetLabelSize(0.035);
  hCls->GetYaxis()->SetNdivisions(505);
  hCls->GetYaxis()->SetTitleSize(0.04);
  hCls->SetFillColorAlpha(kAzure + 2, 0.35);
  hCls->SetFillStyle(1001);
  hCls->SetBarWidth(0.7);
  hCls->SetBarOffset(0.15);
  hCls->SetMaximum(1.45 * hCls->GetMaximum());
  hCls->Draw("hbar");
  {
    TLatex t;
    t.SetTextSize(0.03);
    t.SetTextAlign(12);
    for (size_t k = 0; k < kClasses.size(); ++k) {
      const double n = hCls->GetBinContent(k + 1);
      t.DrawLatex(max(0., n) + 0.02 * hCls->GetMaximum(), k + 0.5, Form("%.0f (%.1f%%)", n, 100 * n / nLambda));
    }
  }
  pad(2, true);
  gPad->SetRightMargin(0.36);
  style(sum.nAll, kBlack, "number of extra tracks;extra tracks;#Lambda events");
  style(sum.nGood, kRed + 1, "number of extra tracks;extra tracks;#Lambda events");
  overlay({sum.nAll, sum.nGood});
  auto *leg1 = legend(2);
  leg1->SetTextSize(0.028);
  leg1->AddEntry(sum.nAll, "all (current ntrack)", "l");
  leg1->AddEntry(sum.nGood, "good (no short/upstream)", "l");
  leg1->Draw();
  canvas.Print(pdf);

  // Page 2: two good opposite-charge extra tracks: exclusive Lambda pi+ pi-?
  page(2, 1, "2 good extra tracks with opposite charge");
  pad(1, false);
  style(sum.mm2, kBlack, "2 good extra tracks (#pm) as #pi^{+}#pi^{-}");
  overlay({sum.mm2});
  {
    const double ymax = sum.mm2->GetMaximum();
    line.SetLineColor(kRed + 1);
    for (double x : {-kExclusiveMM2, kExclusiveMM2}) line.DrawLine(x, 0, x, ymax);
    line.SetLineColor(kMagenta + 1);
    line.DrawLine(kMassPi * kMassPi * 0.94, 0, kMassPi * kMassPi * 0.94, 0.6 * ymax);
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.032);
    t.DrawLatex(0.55, 0.84, "0: #Lambda#pi^{+}#pi^{-} exclusive");
    t.DrawLatex(0.55, 0.79, "0.018: #eta#rightarrow#pi^{+}#pi^{-}#pi^{0}");
    t.DrawLatex(0.55, 0.74, Form("|MM^{2}| < %.2f: %.0f of %.0f", kExclusiveMM2,
                                 sum.mm2->Integral(sum.mm2->FindBin(-kExclusiveMM2 + 1e-6), sum.mm2->FindBin(kExclusiveMM2 - 1e-6)),
                                 sum.mm2->Integral()));
  }
  pad(2, true);
  gPad->SetRightMargin(0.36);
  auto *leg2 = legend(7);
  leg2->SetTextSize(0.028);
  vector<TH1 *> classMM;
  for (int c = 0; c < 7; ++c) {
    style(sum.lambdaMM[c], colors[c], "#Lambda missing mass per class");
    classMM.push_back(sum.lambdaMM[c]);
    leg2->AddEntry(sum.lambdaMM[c], kShort[c], "l");
  }
  sum.lambdaMM[0]->GetXaxis()->SetNdivisions(504);
  overlay(classMM);
  leg2->Draw();
  canvas.Print(pdf);

  // Page 3: one good extra track: from the production vertex or not?
  page(2, 2, "1 good extra track");
  pad(1, true);
  const int sColors[5] = {kBlue + 1, kRed + 1, kGreen + 2, kOrange + 1, kGray + 2};
  auto *leg3 = legend(5);
  vector<TH1 *> dcaSpecies;
  for (int s = 0; s < 5; ++s) {
    style(sum.oneDca[s], sColors[s], "DCA to the production vertex");
    dcaSpecies.push_back(sum.oneDca[s]);
    leg3->AddEntry(sum.oneDca[s], kSpecies[s], "l");
  }
  overlay(dcaSpecies);
  dcaNote(sum.oneDca[0]);
  leg3->Draw();
  auto pair2 = [&](int i, TH1D *a, TH1D *b, const char *title, bool labels) {
    pad(i, true);
    style(a, kRed + 1, title);
    style(b, kBlue + 1, title);
    if (labels) {
      speciesLabels(a);
      speciesLabels(b);
    }
    overlay({a, b});
    auto *l = legend(4);
    l->AddEntry(a, Form("DCA < %.0f mm", kVertexDca), "l");
    l->AddEntry(a, Form("  (%.0f)", a->Integral()), "");
    l->AddEntry(b, Form("DCA #geq %.0f mm", kVertexDca), "l");
    l->AddEntry(b, Form("  (%.0f)", b->Integral()), "");
    l->Draw();
  };
  pair2(2, sum.oneSpecies[1], sum.oneSpecies[0], "species", true);
  pair2(3, sum.oneMom[1], sum.oneMom[0], "momentum", false);
  pair2(4, sum.oneVtz[1], sum.oneVtz[0], "production vertex z", false);
  canvas.Print(pdf);

  // Page 4: two same-charge and 3+ good extra tracks, and all good extra tracks.
  page(2, 2, "2 same-charge and 3+ good extra tracks");
  pad(1, false);
  style(sum.twoDca, kGreen + 2, "2 good, same charge: DCA to the production vertex");
  overlay({sum.twoDca});
  dcaNote(sum.twoDca);
  pad(2, false);
  style(sum.threeN, kViolet, "3+ good: how many from the vertex");
  overlay({sum.threeN});
  pad(3, false);
  style(sum.threeSpecies, kViolet, "3+ good: species of the tracks");
  speciesLabels(sum.threeSpecies);
  overlay({sum.threeSpecies});
  pad(4, false);
  style(sum.goodDca, kBlack, "all good extra tracks: DCA to the production vertex");
  overlay({sum.goodDca});
  dcaNote(sum.goodDca);
  canvas.Print(pdf);

  // Page 5: Lambda MM of 1 good extra track by its DCA, compared with no extra track.
  page(2, 1, "1 good extra track: #Lambda missing mass by DCA");
  const char *oneLabel[3] = {"1 good, DCA < 10 mm", "1 good, DCA 10-30 mm", "1 good, DCA #geq 30 mm"};
  const int oneColor[3] = {kRed + 1, kOrange + 1, kBlue + 1};
  pad(1, true);
  gPad->SetRightMargin(0.36);
  auto *leg5 = legend(4);
  leg5->SetTextSize(0.028);
  style(sum.lambdaMM[0], kBlack, "#Lambda missing mass");
  vector<TH1 *> mmByDca = {sum.lambdaMM[0]};
  leg5->AddEntry(sum.lambdaMM[0], "no extra", "l");
  for (int v = 0; v < 3; ++v) {
    style(sum.oneMM[v], oneColor[v], "#Lambda missing mass");
    mmByDca.push_back(sum.oneMM[v]);
    leg5->AddEntry(sum.oneMM[v], oneLabel[v], "l");
  }
  sum.lambdaMM[0]->GetXaxis()->SetNdivisions(504);
  overlay(mmByDca);
  leg5->Draw();
  pad(2, true);
  gPad->SetRightMargin(0.36);
  auto *leg6 = legend(3);
  leg6->SetTextSize(0.028);
  vector<TH1 *> shapes;
  const int shapeOf[3] = {-1, 0, 2};
  for (int i = 0; i < 3; ++i) {
    TH1 *src = shapeOf[i] < 0 ? (TH1 *)sum.lambdaMM[0] : (TH1 *)sum.oneMM[shapeOf[i]];
    auto *x = (TH1 *)src->Clone(Form("shape%d", i));
    x->Scale(1. / x->Integral());
    style(x, src->GetLineColor(), "#Lambda missing mass shape;MM(K^{-}p#rightarrow#LambdaX) [GeV/#it{c}^{2}];unit area");
    shapes.push_back(x);
    leg6->AddEntry(x, shapeOf[i] < 0 ? "no extra" : oneLabel[shapeOf[i]], "l");
  }
  shapes[0]->GetXaxis()->SetNdivisions(504);
  overlay(shapes);
  leg6->Draw();
  canvas.Print(pdf);

  // Page 6: far single extra track in the Lambda MM continuum and Lambda eta peak.
  page(3, 2, Form("1 good extra track with DCA #geq %.0f mm: #Lambda MM %.2f-%.2f vs %.3f-%.3f (unit area)", kFarDca,
                  kContMin, kContMax, kPeakMin, kPeakMax));
  auto farPair = [&](int i, TH1D *const *x, const char *title, bool labels) {
    pad(i, true);
    vector<TH1 *> hs;
    for (int r = 0; r < 2; ++r) {
      auto *c = (TH1D *)x[r]->Clone(Form("%s_draw", x[r]->GetName()));
      if (c->Integral() != 0) c->Scale(1. / c->Integral());
      style(c, r ? kRed + 1 : kBlack, title);
      if (labels) speciesLabels(c);
      hs.push_back(c);
    }
    overlay(hs);
    auto *l = legend(4);
    l->SetTextSize(0.04);
    l->AddEntry(hs[0], "MM continuum", "l");
    l->AddEntry(hs[0], Form("  (%.0f)", x[0]->Integral()), "");
    l->AddEntry(hs[1], "MM #Lambda#eta peak", "l");
    l->AddEntry(hs[1], Form("  (%.0f)", x[1]->Integral()), "");
    l->Draw();
  };
  farPair(1, sum.farSpecies, "species", true);
  farPair(2, sum.farCharge, "charge", false);
  farPair(3, sum.farMom, "momentum", false);
  farPair(4, sum.farInner, "innermost hit to the vertex", false);
  farPair(5, sum.farPileup, "closest pile-up beam track (last bin: none)", false);
  farPair(6, sum.farHits, "number of hits", false);
  canvas.Print(pdf);
  page(2, 1, Form("1 good extra track with DCA #geq %.0f mm: dE/dx", kFarDca));
  for (int r = 0; r < 2; ++r) {
    pad(r + 1, false);
    gPad->SetRightMargin(0.12);
    style(sum.farDedx[r], kBlack, r ? "#Lambda MM #Lambda#eta peak" : "#Lambda MM continuum");
    sum.farDedx[r]->SetMinimum(0);
    sum.farDedx[r]->Draw("colz");
  }
  canvas.Print(pdf);

  canvas.Print(pdf + "]");
  for (int r = 0; r < 2; ++r) {
    const double n = sum.farSpecies[r]->Integral();
    printf("far single track, %s: N=%.0f", r ? "peak" : "continuum", n);
    for (int sp = 0; sp < 5; ++sp) printf("  %s %.2f", kSpecies[sp], sum.farSpecies[r]->GetBinContent(sp + 1) / n);
    printf("  q- %.2f q+ %.2f", sum.farCharge[r]->GetBinContent(1) / n, sum.farCharge[r]->GetBinContent(3) / n);
    printf("  pileup<10mm %.2f  inner<50mm %.2f  p<0.1 %.2f\n", sum.farPileup[r]->Integral(1, 5) / n,
           sum.farInner[r]->Integral(1, 5) / n, sum.farMom[r]->Integral(1, 2) / n);
  }
  auto window = [](TH1 *x) {
    return x->Integral(x->FindBin(0.5301), x->FindBin(0.5649));
  };
  printf("Lambda MM 0.53-0.565 / total(0.5-0.57):\n");
  printf("  no extra             %7.0f / %7.0f\n", window(sum.lambdaMM[0]), sum.lambdaMM[0]->Integral());
  for (int v = 0; v < 3; ++v) printf("  %-36s %7.0f / %7.0f\n", oneLabel[v], window(sum.oneMM[v]), sum.oneMM[v]->Integral());
  printf("Lambda MM bins (low edge: no extra, 1 good DCA <10, 10-30, >=30 mm):\n");
  for (int b = sum.lambdaMM[0]->FindBin(0.5001); b <= sum.lambdaMM[0]->GetNbinsX(); ++b)
    printf("  %.3f %7.1f %7.1f %7.1f %7.1f\n", sum.lambdaMM[0]->GetBinLowEdge(b), sum.lambdaMM[0]->GetBinContent(b),
           sum.oneMM[0]->GetBinContent(b), sum.oneMM[1]->GetBinContent(b), sum.oneMM[2]->GetBinContent(b));
  cout << "Lambda events " << nLambda << endl;
  for (size_t k = 0; k < kClasses.size(); ++k)
    printf("  %-45s %7.0f  %5.1f%%\n", kClasses[k], sum.cls->GetBinContent(k + 1), 100 * sum.cls->GetBinContent(k + 1) / nLambda);
  cout << "Wrote " << pdf << endl;
}
