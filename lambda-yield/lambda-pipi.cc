// Missing mass squared of K- p -> Lambda pi+ pi- X for Lambda candidates with
// exactly two good extra tracks of opposite charge (the "2 good +-" classes of
// ntrack-class.cc), to see what sits next to the exclusive Lambda pi+ pi- peak
// at MM2 = 0, e.g. X = pi0 (MM2 = 0.018) from eta -> pi+ pi- pi0.
// Selection: DCA, chisqr (and LH2), no X_mass cut; histograms are filled with +1
// in the Lambda window and -1 in the sidebands, so they are sideband subtracted.
//   root -l -b -q -e '.L lambda-pipi.cc+' -e 'lambda_pipi()'        // LH2
//   root -l -b -q -e '.L lambda-pipi.cc+' -e 'lambda_pipi("result", false)' // all vertices
#include <TCanvas.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLorentzVector.h>
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
using extra_track::IsUpstream;

namespace {
const vector<int> kRuns = {2447, 2449, 2450, 2451, 2452};
const char *kDataDir = "/home/had/haein/data/JPARC2025Nov_root/physics-735";
const double kSigMin = 1.105, kSigMax = 1.125;
bool InSideband(double m) { return (m > 1.085 && m < 1.095) || (m > 1.135 && m < 1.145); }
const double kMassK = 0.493677, kMassP = 0.93827208, kMassPi = 0.13957039, kMassPi0 = 0.1349768,
             kMassLambda = 1.115683, kMassEta = 0.547862;
// X_mass (Lambda missing mass) regions; the eta peak is the one of ntrack-class.cc.
const vector<double> kMmEdges = {-1., 0.5, 0.545, 0.565};
const vector<const char *> kMmNames = {"MM < 0.5", "0.5 #leq MM < 0.545", "#eta peak 0.545-0.565"};
// MM2(Lambda pi+ pi-) windows for the pi+ pi- mass: around 0 and around m_pi0^2.
const double kZeroMax = 0.008, kPi0Min = 0.012, kPi0Max = 0.028; // GeV2/c4, on bin edges

struct Hists {
  TH1D *mm2All, *mm2Pid, *mm2Mm[3], *mpipi[2], *mmLambda[2];
  TH2D *mm2VsMm;
};
TH1D *H(const char *name, const char *title, int n, double lo, double hi)
{
  auto *h = new TH1D(name, title, n, lo, hi);
  h->SetDirectory(nullptr);
  return h;
}
Hists MakeHists(int run)
{
  Hists h;
  const char *mm2Title = ";MM^{2}(K^{-}p#rightarrow#Lambda#pi^{+}#pi^{-}X) [GeV^{2}/#it{c}^{4}];#Lambda events / 0.004";
  h.mm2All = H(Form("mm2All_%d", run), mm2Title, 75, -0.1, 0.2);
  h.mm2Pid = H(Form("mm2Pid_%d", run), mm2Title, 75, -0.1, 0.2);
  for (int r = 0; r < 3; ++r) h.mm2Mm[r] = H(Form("mm2Mm%d_%d", r, run), mm2Title, 75, -0.1, 0.2);
  for (int r = 0; r < 2; ++r) {
    h.mpipi[r] = H(Form("mpipi%d_%d", r, run), ";#it{m}_{#pi^{+}#pi^{-}} [GeV/#it{c}^{2}];#Lambda events / 10 MeV/#it{c}^{2}",
                   50, 0.25, 0.75);
    h.mmLambda[r] = H(Form("mmLambda%d_%d", r, run), ";MM(K^{-}p#rightarrow#LambdaX) [GeV/#it{c}^{2}];#Lambda events / 10 MeV/#it{c}^{2}",
                      50, 0.1, 0.6);
  }
  h.mm2VsMm = new TH2D(Form("mm2VsMm_%d", run),
                       ";MM(K^{-}p#rightarrow#LambdaX) [GeV/#it{c}^{2}];MM^{2}(#Lambda#pi^{+}#pi^{-}X) [GeV^{2}/#it{c}^{4}]",
                       50, 0.1, 0.6, 50, -0.1, 0.2);
  h.mm2VsMm->SetDirectory(nullptr);
  return h;
}
vector<TH1 *> All(const Hists &h)
{
  vector<TH1 *> v = {h.mm2All, h.mm2Pid, h.mm2Mm[0], h.mm2Mm[1], h.mm2Mm[2], h.mpipi[0], h.mpipi[1],
                     h.mmLambda[0], h.mmLambda[1], h.mm2VsMm};
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
  vector<double> *bx = nullptr, *by = nullptr, *bz = nullptr;
  vector<int> *flag = nullptr, *protonId = nullptr, *pionId = nullptr;
  vector<int> *isBeam = nullptr, *isK18 = nullptr, *isAccidental = nullptr, *pid = nullptr, *charge = nullptr;
  vector<vector<double>> *momX = nullptr, *momY = nullptr, *momZ = nullptr;
  vector<vector<double>> *hitZ = nullptr, *clX = nullptr, *clZ = nullptr;
  Int_t ntTpc = 0;
  tree->SetBranchStatus("*", 0);
  for (const char *b : {"lambda_mass", "X_mass", "lambda_mom_x", "lambda_mom_y", "lambda_mom_z", "X_prod_beam_px",
                        "X_prod_beam_py", "X_prod_beam_pz", "cut_flag", "lambda_proton_track_id",
                        "lambda_pion_track_id", "is_beam", "is_k18", "is_accidental", "pid", "charge", "mom_vtx",
                        "mom_vty", "mom_vtz", "ntTpc", "hitpos_z", "track_cluster_x_center", "track_cluster_z_center"})
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
  tree->SetBranchAddress("mom_vtx", &momX);
  tree->SetBranchAddress("mom_vty", &momY);
  tree->SetBranchAddress("mom_vtz", &momZ);
  tree->SetBranchAddress("ntTpc", &ntTpc);
  tree->SetBranchAddress("hitpos_z", &hitZ);
  tree->SetBranchAddress("track_cluster_x_center", &clX);
  tree->SetBranchAddress("track_cluster_z_center", &clZ);
  tree->SetCacheSize(200 * 1024 * 1024);
  tree->AddBranchToCache("*", kTRUE);
  tree->StopCacheLearningPhase();
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    if (flag->empty() || mass->empty() || protonId->empty() || pionId->empty() || missingMass->empty()) continue;
    const int f = flag->at(0);
    if ((f / 1000) % 10 != 1 || (f / 100) % 10 != 1) continue; // DCA and chisqr
    if (lh2 && f % 10 != 1) continue;
    const double m = mass->at(0);
    const double w = (m > kSigMin && m < kSigMax) ? 1. : InSideband(m) ? -1. : 0.;
    if (w == 0.) continue;
    const int ip = protonId->at(0), ipi = pionId->at(0);
    vector<int> good;
    for (int it = 0; it < ntTpc; ++it) {
      if (it == ip || it == ipi || isBeam->at(it) || isK18->at(it) || isAccidental->at(it)) continue;
      if (hitZ->at(it).size() <= extra_track::kShortHits || IsUpstream(clX->at(it), clZ->at(it))) continue;
      good.push_back(it);
    }
    if (good.size() != 2 || charge->at(good[0]) * charge->at(good[1]) >= 0) continue;
    const int a = good[0], b = good[1];
    TLorentzVector A, B, K, L;
    A.SetXYZM(momX->at(a).at(b), momY->at(a).at(b), momZ->at(a).at(b), kMassPi);
    B.SetXYZM(momX->at(b).at(a), momY->at(b).at(a), momZ->at(b).at(a), kMassPi);
    K.SetXYZM(bx->at(0), by->at(0), bz->at(0), kMassK);
    L.SetXYZM(lx->at(0), ly->at(0), lz->at(0), kMassLambda);
    const double x = (K + TLorentzVector(0, 0, 0, kMassP) - L - A - B).M2();
    const double mm = missingMass->at(0);
    if (!std::isfinite(x) || !std::isfinite(mm)) continue;
    h.mm2All->Fill(x, w);
    if ((pid->at(a) & 0x1) && (pid->at(b) & 0x1)) h.mm2Pid->Fill(x, w);
    for (int r = 0; r < 3; ++r) {
      const double hi = r + 1 < (int)kMmEdges.size() ? kMmEdges[r + 1] : 1e9;
      if (mm >= kMmEdges[r] && mm < hi) h.mm2Mm[r]->Fill(x, w);
    }
    h.mm2VsMm->Fill(mm, x, w);
    const int z = std::abs(x) < kZeroMax ? 0 : (x > kPi0Min && x < kPi0Max) ? 1 : -1;
    if (z >= 0) {
      h.mpipi[z]->Fill((A + B).M(), w);
      h.mmLambda[z]->Fill(mm, w);
    }
  }
  cout << "run" << run << " done" << endl;
  return h;
}
} // namespace

void lambda_pipi(const char *outDir = "result", bool lh2 = true)
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

  const TString sel = Form("DCA, #chi^{2}%s, no MM cut, #Lambda %.3f-%.3f minus sidebands, 2 good extra tracks (#pm)",
                           lh2 ? ", LH2" : "", kSigMin, kSigMax);
  const TString runTag = Form("run%05d-%05d-n%zu", kRuns.front(), kRuns.back(), kRuns.size());
  const TString pdf = Form("%s/lambda-pipi-%s%s.pdf", outDir, runTag.Data(), lh2 ? "-lh2" : "");
  TCanvas canvas("c_lambda_pipi", "", 1800, 1200);
  canvas.Print(pdf + "[");
  gStyle->SetTitleFontSize(0.05);
  TGaxis::SetMaxDigits(4);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleX(0.5);
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
  auto pad = [&](int i, bool legendRoom) {
    body->cd(i);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.09);
    gPad->SetRightMargin(legendRoom ? 0.32 : 0.04);
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
  // Mark MM2 = 0 (Lambda pi+ pi-) and m_pi0^2 (Lambda pi+ pi- pi0).
  auto mm2Lines = [&](double ymax) {
    TLine l;
    l.SetLineStyle(2);
    l.SetLineColor(kMagenta + 1);
    l.DrawLine(0, 0, 0, ymax);
    l.SetLineColor(kGreen + 2);
    l.DrawLine(kMassPi0 * kMassPi0, 0, kMassPi0 * kMassPi0, ymax);
  };
  auto count = [](TH1 *h, double lo, double hi) { return h->Integral(h->FindBin(lo + 1e-6), h->FindBin(hi - 1e-6)); };

  // Page 1: MM2 of all, and with both extra tracks pion-PID.
  page(2, 1, "MM^{2}(#Lambda#pi^{+}#pi^{-}): magenta 0, green #it{m}_{#pi^{0}}^{2} = 0.018");
  pad(1, true);
  style(sum.mm2All, kBlack, "2 good extra tracks (#pm) as #pi^{+}#pi^{-}");
  style(sum.mm2Pid, kRed + 1, "2 good extra tracks (#pm) as #pi^{+}#pi^{-}");
  overlay({sum.mm2All, sum.mm2Pid});
  mm2Lines(sum.mm2All->GetMaximum());
  {
    auto *l = legend(4);
    l->SetTextSize(0.03);
    l->AddEntry(sum.mm2All, Form("all: %.0f", sum.mm2All->Integral()), "l");
    l->AddEntry(sum.mm2Pid, Form("both #pi PID: %.0f", sum.mm2Pid->Integral()), "l");
    l->AddEntry((TObject *)nullptr, Form("|MM^{2}| < %.3f: %.0f", kZeroMax, count(sum.mm2All, -kZeroMax, kZeroMax)), "");
    l->AddEntry((TObject *)nullptr, Form("%.3f-%.3f: %.0f", kPi0Min, kPi0Max, count(sum.mm2All, kPi0Min, kPi0Max)), "");
    l->Draw();
  }
  pad(2, false);
  gPad->SetRightMargin(0.12);
  sum.mm2VsMm->SetTitle("MM^{2}(#Lambda#pi^{+}#pi^{-}) vs MM(#Lambda)");
  for (auto *a : {sum.mm2VsMm->GetXaxis(), sum.mm2VsMm->GetYaxis()}) {
    a->SetTitleSize(0.05);
    a->SetLabelSize(0.045);
  }
  sum.mm2VsMm->GetYaxis()->SetTitleOffset(1.45);
  sum.mm2VsMm->SetMinimum(0);
  sum.mm2VsMm->Draw("colz");
  {
    TLine l;
    l.SetLineStyle(2);
    l.SetLineColor(kGreen + 2);
    l.DrawLine(0.1, kMassPi0 * kMassPi0, 0.6, kMassPi0 * kMassPi0);
    l.SetLineColor(kMagenta + 1);
    l.DrawLine(0.1, 0, 0.6, 0);
    l.DrawLine(kMassEta, -0.1, kMassEta, 0.2);
  }
  canvas.Print(pdf);

  // Page 2: MM2 in regions of the Lambda missing mass.
  page(3, 1, "MM^{2}(#Lambda#pi^{+}#pi^{-}) by MM(#Lambda) region");
  const int colors[3] = {kBlack, kAzure + 2, kRed + 1};
  for (int r = 0; r < 3; ++r) {
    pad(r + 1, false);
    style(sum.mm2Mm[r], colors[r], Form("%s: %.0f (0: %.0f, #pi^{0}: %.0f)", kMmNames[r], sum.mm2Mm[r]->Integral(),
                                        count(sum.mm2Mm[r], -kZeroMax, kZeroMax), count(sum.mm2Mm[r], kPi0Min, kPi0Max)));
    overlay({sum.mm2Mm[r]});
    mm2Lines(sum.mm2Mm[r]->GetMaximum());
  }
  canvas.Print(pdf);

  // Page 3: m(pi+ pi-) and MM(Lambda) near MM2 = 0 and near m_pi0^2.
  page(2, 1, Form("|MM^{2}| < %.3f vs MM^{2} %.3f-%.3f GeV^{2}/#it{c}^{4}", kZeroMax, kPi0Min, kPi0Max));
  pad(1, true);
  style(sum.mpipi[0], kMagenta + 1, "#it{m}_{#pi^{+}#pi^{-}}");
  style(sum.mpipi[1], kGreen + 2, "#it{m}_{#pi^{+}#pi^{-}}");
  overlay({sum.mpipi[0], sum.mpipi[1]});
  {
    TLine l;
    l.SetLineStyle(3);
    l.DrawLine(kMassEta - kMassPi0, 0, kMassEta - kMassPi0, sum.mpipi[0]->GetMaximum());
    auto *leg = legend(3);
    leg->SetTextSize(0.03);
    leg->AddEntry(sum.mpipi[0], Form("MM^{2} #approx 0: %.0f", sum.mpipi[0]->Integral()), "l");
    leg->AddEntry(sum.mpipi[1], Form("MM^{2} #approx m_{#pi^{0}}^{2}: %.0f", sum.mpipi[1]->Integral()), "l");
    leg->AddEntry((TObject *)nullptr, "dotted: m_{#eta} - m_{#pi^{0}}", "");
    leg->Draw();
  }
  pad(2, true);
  style(sum.mmLambda[0], kMagenta + 1, "MM(#Lambda)");
  style(sum.mmLambda[1], kGreen + 2, "MM(#Lambda)");
  overlay({sum.mmLambda[0], sum.mmLambda[1]});
  {
    TLine l;
    l.SetLineStyle(3);
    l.DrawLine(kMassEta, 0, kMassEta, sum.mmLambda[0]->GetMaximum());
    auto *leg = legend(3);
    leg->SetTextSize(0.03);
    leg->AddEntry(sum.mmLambda[0], "MM^{2} #approx 0", "l");
    leg->AddEntry(sum.mmLambda[1], "MM^{2} #approx m_{#pi^{0}}^{2}", "l");
    leg->AddEntry((TObject *)nullptr, "dotted: m_{#eta}", "");
    leg->Draw();
  }
  canvas.Print(pdf);
  canvas.Print(pdf + "]");

  cout << "N 2 good +-: " << sum.mm2All->Integral() << ", |MM2| < " << kZeroMax << ": "
       << count(sum.mm2All, -kZeroMax, kZeroMax) << ", " << kPi0Min << "-" << kPi0Max << ": "
       << count(sum.mm2All, kPi0Min, kPi0Max) << endl;
  for (int r = 0; r < 3; ++r)
    cout << kMmNames[r] << ": " << sum.mm2Mm[r]->Integral() << " (0: " << count(sum.mm2Mm[r], -kZeroMax, kZeroMax)
         << ", pi0: " << count(sum.mm2Mm[r], kPi0Min, kPi0Max) << ")" << endl;
  cout << pdf << endl;
}
