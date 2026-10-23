// What are the extra tracks that make the ntrack cut of DstTPCLambdaEta.cc fail?
// Lambda events passing DCA and chisqr are split into a signal window and
// sidebands of the same total width; the sideband-subtracted distributions show
// the extra (not beam-like, not p/pi of the Lambda) tracks of real Lambda events.
//   root -l -b -q -e '.L ntrack-extra.cc+' -e 'ntrack_extra()'
// ntrack_extra("result", true) keeps only the Lambda eta region, X_mass >= 0.5 GeV/c2.
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
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

using namespace std;

namespace {
const vector<int> kRuns = {2447, 2449, 2450, 2451};
const char *kDataDir = "/home/had/haein/data/JPARC2025Nov_root/physics-735";
// Signal window and two sidebands, 20 MeV/c2 in total each.
const double kSigMin = 1.105, kSigMax = 1.125;
const double kMissingMassCut = 0.5; // GeV/c2, Lambda eta
bool InSideband(double m) { return (m > 1.085 && m < 1.095) || (m > 1.135 && m < 1.145); }

// Extra track categories from the stored pid bits (1 pi, 2 K, 4 p) and charge.
const vector<const char *> kCategories = {"#pi^{-}", "#pi^{+}", "p", "K^{#pm}", "e/other (-)", "e/other (+)"};
int Category(int pid, int charge)
{
  if (pid & 0x4 && charge > 0) return 2;
  if (pid & 0x1) return charge < 0 ? 0 : 1;
  if (pid & 0x2) return 3;
  return charge < 0 ? 4 : 5;
}

// Extra track classes, checked in this order; ignoring classes 0..k-1 is an
// ntrack cut that does not count them.
constexpr double kZTarget = -143.;
const vector<const char *> kClasses = {"upstream beam window", "duplicate of p/#pi", "short (#leq 6 hits)",
                                       "from vertex (DCA_{beam} < 10 mm)", "other"};
constexpr int kNClasses = 5;
double Angle(double ax, double ay, double az, double bx, double by, double bz)
{
  const double c = (ax * bx + ay * by + az * bz) / sqrt((ax * ax + ay * ay + az * az) * (bx * bx + by * by + bz * bz));
  return acos(max(-1., min(1., c)));
}
// Same test as IsUpstreamBeamWindowTrack of DstTPCLambdaEta.cc.
bool IsUpstream(const vector<double> &x, const vector<double> &z)
{
  bool any = false;
  for (size_t i = 0; i < min(x.size(), z.size()); ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(z[i])) continue;
    any = true;
    if (z[i] >= kZTarget || x[i] <= -50. || x[i] >= 50.) return false;
  }
  return any;
}
// Classes 0-3; the caller splits class 3 by the DCA to the beam.
int Classify(int it, int ip, int ipi, int q, int qp, int qpi, const vector<vector<double>> &clX,
             const vector<vector<double>> &clZ, const vector<vector<double>> &hitZ,
             const vector<vector<double>> &dist, const vector<vector<double>> &px, const vector<vector<double>> &py,
             const vector<vector<double>> &pz)
{
  if (IsUpstream(clX[it], clZ[it])) return 0;
  // Same charge as a daughter, close to it and nearly parallel at the closest point.
  for (auto [d, qd] : {pair<int, int>{ip, qp}, pair<int, int>{ipi, qpi}})
    if (q == qd && dist[it][d] < 3. &&
        Angle(px[it][d], py[it][d], pz[it][d], px[d][it], py[d][it], pz[d][it]) < 0.2)
      return 1;
  if (hitZ[it].size() <= 6) return 2;
  return 3;
}

// Index 0: signal window, 1: sidebands.
struct Hists {
  TH1D *nExtra[2][2]; // [LH2 fail/pass][window]
  TH1D *category[2], *mom[2], *dcaBeam[2], *dcaLambda[2], *vtz[2], *chisqr[2], *nhit[2], *firstZ[2];
  TH2D *dedx[2];
  TH1D *trackClass[2], *zeroWhenIgnored[2], *events[2];
};
Hists MakeHists(int run)
{
  Hists h;
  for (int w = 0; w < 2; ++w) {
    const char *tag = w == 0 ? "sig" : "sb";
    for (int l = 0; l < 2; ++l)
      h.nExtra[l][w] = new TH1D(Form("nextra_%d_%s_%d", l, tag, run), "", 8, -0.5, 7.5);
    h.category[w] = new TH1D(Form("cat_%s_%d", tag, run), "", kCategories.size(), -0.5, kCategories.size() - 0.5);
    h.mom[w] = new TH1D(Form("mom_%s_%d", tag, run), "", 60, 0, 1.5);
    h.dcaBeam[w] = new TH1D(Form("dcab_%s_%d", tag, run), "", 50, 0, 100);
    h.dcaLambda[w] = new TH1D(Form("dcal_%s_%d", tag, run), "", 50, 0, 100);
    h.vtz[w] = new TH1D(Form("vtz_%s_%d", tag, run), "", 70, -300, 400);
    h.chisqr[w] = new TH1D(Form("chi_%s_%d", tag, run), "", 50, 0, 10);
    h.nhit[w] = new TH1D(Form("nhit_%s_%d", tag, run), "", 33, -0.5, 32.5);
    h.firstZ[w] = new TH1D(Form("fz_%s_%d", tag, run), "", 60, -300, 300);
    h.trackClass[w] = new TH1D(Form("class_%s_%d", tag, run), "", kNClasses, -0.5, kNClasses - 0.5);
    h.zeroWhenIgnored[w] = new TH1D(Form("zero_%s_%d", tag, run), "", kNClasses, -0.5, kNClasses - 0.5);
    h.events[w] = new TH1D(Form("events_%s_%d", tag, run), "", 1, -0.5, 0.5);
    h.dedx[w] = new TH2D(Form("dedx_%s_%d", tag, run), "", 100, -1.5, 1.5, 100, 0, 200);
  }
  return h;
}
vector<TH1 *> All(Hists &h)
{
  vector<TH1 *> v;
  for (int w = 0; w < 2; ++w) {
    v.push_back(h.nExtra[0][w]);
    v.push_back(h.nExtra[1][w]);
    for (TH1 *x : {(TH1 *)h.category[w], (TH1 *)h.mom[w], (TH1 *)h.dcaBeam[w], (TH1 *)h.dcaLambda[w], (TH1 *)h.vtz[w],
                   (TH1 *)h.chisqr[w], (TH1 *)h.nhit[w], (TH1 *)h.firstZ[w], (TH1 *)h.dedx[w],
                   (TH1 *)h.trackClass[w], (TH1 *)h.zeroWhenIgnored[w], (TH1 *)h.events[w]})
      v.push_back(x);
  }
  return v;
}

Hists ReadRun(int run, bool mmCut)
{
  Hists h = MakeHists(run);
  for (TH1 *x : All(h)) x->SetDirectory(nullptr);
  TFile file(Form("%s/run%05d_DstTPCLambdaEta.root", kDataDir, run));
  TTree *tree = file.IsZombie() ? nullptr : file.Get<TTree>("tpc");
  if (!tree || !tree->GetBranch("cut_flag")) {
    cerr << "No tpc tree with cut_flag for run " << run << endl;
    return h;
  }
  vector<double> *mass = nullptr, *missingMass = nullptr, *chisqr = nullptr, *mom0 = nullptr, *dEdx = nullptr;
  vector<int> *flag = nullptr, *protonId = nullptr, *pionId = nullptr;
  vector<int> *isBeam = nullptr, *isK18 = nullptr, *isAccidental = nullptr, *pid = nullptr, *charge = nullptr;
  vector<vector<double>> *closeDist = nullptr, *vtz = nullptr, *hitZ = nullptr, *clX = nullptr, *clZ = nullptr;
  vector<vector<double>> *momX = nullptr, *momY = nullptr, *momZ = nullptr;
  Int_t ntTpc = 0;
  tree->SetBranchStatus("*", 0);
  for (const char *b : {"lambda_mass", "X_mass", "cut_flag", "lambda_proton_track_id", "lambda_pion_track_id", "ntTpc", "is_beam",
                        "is_k18", "is_accidental", "pid", "charge", "chisqr", "mom0", "dEdx", "closeDistTpc", "vtzTpc",
                        "hitpos_z", "track_cluster_x_center", "track_cluster_z_center", "mom_vtx", "mom_vty",
                        "mom_vtz"})
    tree->SetBranchStatus(b, 1);
  tree->SetBranchAddress("lambda_mass", &mass);
  tree->SetBranchAddress("X_mass", &missingMass);
  tree->SetBranchAddress("cut_flag", &flag);
  tree->SetBranchAddress("lambda_proton_track_id", &protonId);
  tree->SetBranchAddress("lambda_pion_track_id", &pionId);
  tree->SetBranchAddress("ntTpc", &ntTpc);
  tree->SetBranchAddress("is_beam", &isBeam);
  tree->SetBranchAddress("is_k18", &isK18);
  tree->SetBranchAddress("is_accidental", &isAccidental);
  tree->SetBranchAddress("pid", &pid);
  tree->SetBranchAddress("charge", &charge);
  tree->SetBranchAddress("chisqr", &chisqr);
  tree->SetBranchAddress("mom0", &mom0);
  tree->SetBranchAddress("dEdx", &dEdx);
  tree->SetBranchAddress("closeDistTpc", &closeDist);
  tree->SetBranchAddress("vtzTpc", &vtz);
  tree->SetBranchAddress("hitpos_z", &hitZ);
  tree->SetBranchAddress("track_cluster_x_center", &clX);
  tree->SetBranchAddress("track_cluster_z_center", &clZ);
  tree->SetBranchAddress("mom_vtx", &momX);
  tree->SetBranchAddress("mom_vty", &momY);
  tree->SetBranchAddress("mom_vtz", &momZ);
  tree->SetCacheSize(200 * 1024 * 1024);
  tree->AddBranchToCache("*", kTRUE);
  tree->StopCacheLearningPhase();
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    if (flag->empty() || mass->empty() || protonId->empty() || pionId->empty()) continue;
    const int f = flag->at(0);
    if ((f / 1000) % 10 != 1 || (f / 100) % 10 != 1) continue; // DCA and chisqr
    if (mmCut && (missingMass->empty() || !(missingMass->at(0) >= kMissingMassCut))) continue;
    const double m = mass->at(0);
    const int w = (m > kSigMin && m < kSigMax) ? 0 : InSideband(m) ? 1 : -1;
    if (w < 0) continue;
    const int ip = protonId->at(0), ipi = pionId->at(0);
    auto beamLike = [&](int it) { return isBeam->at(it) || isK18->at(it) || isAccidental->at(it); };
    // Beam track for the DCA: the K18-matched track, else any beam track.
    int beam = -1;
    for (int it = 0; it < ntTpc && beam < 0; ++it)
      if (isK18->at(it)) beam = it;
    for (int it = 0; it < ntTpc && beam < 0; ++it)
      if (isBeam->at(it)) beam = it;
    int nExtra = 0;
    int nLeft[kNClasses] = {0};
    for (int it = 0; it < ntTpc; ++it) {
      if (it == ip || it == ipi || beamLike(it)) continue;
      ++nExtra;
      int c = Classify(it, ip, ipi, charge->at(it), charge->at(ip), charge->at(ipi), *clX, *clZ, *hitZ, *closeDist,
                       *momX, *momY, *momZ);
      if (c == 3 && (beam < 0 || closeDist->at(it).at(beam) >= 10.)) c = 4;
      h.trackClass[w]->Fill(c);
      // A class-c track is still counted when only the classes before c are ignored.
      for (int k = 0; k <= c; ++k) ++nLeft[k];
      h.category[w]->Fill(Category(pid->at(it), charge->at(it)));
      h.mom[w]->Fill(mom0->at(it));
      h.dedx[w]->Fill(charge->at(it) * mom0->at(it), dEdx->at(it));
      h.chisqr[w]->Fill(chisqr->at(it));
      h.nhit[w]->Fill(hitZ->at(it).size());
      if (!hitZ->at(it).empty()) h.firstZ[w]->Fill(*min_element(hitZ->at(it).begin(), hitZ->at(it).end()));
      h.dcaLambda[w]->Fill(min(closeDist->at(it).at(ip), closeDist->at(it).at(ipi)));
      if (beam >= 0) {
        h.dcaBeam[w]->Fill(closeDist->at(it).at(beam));
        h.vtz[w]->Fill(vtz->at(it).at(beam));
      }
    }
    h.nExtra[f % 10][w]->Fill(nExtra);
    for (int k = 0; k < kNClasses; ++k)
      if (nLeft[k] == 0) h.zeroWhenIgnored[w]->Fill(k);
    h.events[w]->Fill(0);
  }
  return h;
}

// Sideband-subtracted copy (the windows have the same width).
TH1 *Subtract(TH1 *sig, TH1 *sb, const char *title)
{
  TH1 *d = (TH1 *)sig->Clone(Form("%s_sub", sig->GetName()));
  d->Add(sb, -1);
  d->SetTitle(title);
  return d;
}
} // namespace

void ntrack_extra(const char *outDir = "result", bool mmCut = false)
{
  ROOT::EnableThreadSafety();
  gStyle->SetOptStat(0);
  gStyle->SetPaperSize(30, 20); // landscape page matching the 1800x1200 canvas
  gSystem->mkdir(outDir, kTRUE);
  vector<future<Hists>> jobs;
  for (int run : kRuns) jobs.push_back(async(launch::async, ReadRun, run, mmCut));
  Hists sum = MakeHists(0);
  vector<TH1 *> total = All(sum);
  for (auto &job : jobs) {
    Hists h = job.get();
    vector<TH1 *> part = All(h);
    for (size_t i = 0; i < total.size(); ++i) total[i]->Add(part[i]);
  }

  // Summary: fraction of Lambda events (sideband subtracted) with n extra tracks.
  TPaveText text(0.05, 0.05, 0.95, 0.95, "NDC");
  text.SetTextAlign(12);
  text.SetTextFont(42);
  text.SetFillColor(0);
  text.AddText(Form("runs %d-%d (n=%zu), DCA and #chi^{2}%s passed, signal %.3f-%.3f minus sidebands", kRuns.front(),
                    kRuns.back(), kRuns.size(), mmCut ? ", MM #geq 0.5" : "", kSigMin, kSigMax));
  TH1 *nAll = Subtract(sum.nExtra[0][0], sum.nExtra[0][1], "");
  nAll->Add(sum.nExtra[1][0]);
  nAll->Add(sum.nExtra[1][1], -1);
  for (int l = -1; l < 2; ++l) {
    TH1 *n = l < 0 ? nAll : Subtract(sum.nExtra[l][0], sum.nExtra[l][1], "");
    const double all = n->Integral();
    TString line = Form("%-12s N_{#Lambda} = %7.0f:", l < 0 ? "all" : l == 1 ? "LH2" : "outside LH2", all);
    for (int b = 1; b <= 5; ++b)
      line += Form("  n_{extra}=%d%s %.1f%%", b - 1, b == 5 ? "+" : "",
                   100 * (b == 5 ? n->Integral(5, n->GetNbinsX()) : n->GetBinContent(b)) / all);
    cout << line << endl;
    text.AddText(line);
  }
  TH1 *cat = Subtract(sum.category[0], sum.category[1], "");
  TString line = "extra tracks:";
  for (size_t c = 0; c < kCategories.size(); ++c)
    line += Form("  %s %.1f%%", kCategories[c], 100 * cat->GetBinContent(c + 1) / cat->Integral());
  cout << line << endl;
  text.AddText(line);
  TH1 *cls = Subtract(sum.trackClass[0], sum.trackClass[1], "");
  TH1 *zero = Subtract(sum.zeroWhenIgnored[0], sum.zeroWhenIgnored[1], "");
  const double nLambda = sum.events[0]->GetBinContent(1) - sum.events[1]->GetBinContent(1);
  for (int k = 0; k < kNClasses; ++k) {
    TString l = Form("class %d %-34s %5.1f%% of extra tracks", k, kClasses[k], 100 * cls->GetBinContent(k + 1) / cls->Integral());
    if (k > 0) l += Form("   |  ntrack pass when classes 0-%d ignored: %.1f%%", k - 1, 100 * zero->GetBinContent(k + 1) / nLambda);
    else l += Form("   |  ntrack pass (current): %.1f%%", 100 * zero->GetBinContent(1) / nLambda);
    cout << l << endl;
    text.AddText(l);
  }

  const TString pdf = Form("%s/ntrack-extra-run%05d-%05d-n%zu%s.pdf", outDir, kRuns.front(), kRuns.back(), kRuns.size(),
                          mmCut ? "-mm" : "");
  TCanvas canvas("c_ntrack_extra", "", 1800, 1200);
  canvas.Print(pdf + "[");
  text.Draw();
  canvas.Print(pdf);

  canvas.Clear();
  canvas.Divide(3, 2);
  canvas.cd(1);
  TH1 *nIn = Subtract(sum.nExtra[1][0], sum.nExtra[1][1], "n_{extra} (sideband subtracted);extra tracks;#Lambda events");
  TH1 *nOut = Subtract(sum.nExtra[0][0], sum.nExtra[0][1], "");
  nIn->Scale(1. / nIn->Integral());
  nOut->Scale(1. / nOut->Integral());
  nIn->SetLineColor(kBlue + 1);
  nOut->SetLineColor(kRed + 1);
  nIn->SetMaximum(1);
  nIn->GetYaxis()->SetTitle("fraction");
  nIn->Draw("HIST");
  nOut->Draw("HIST SAME");
  auto *legend = new TLegend(0.5, 0.7, 0.88, 0.88);
  legend->AddEntry(nIn, "production vtx in LH2", "l");
  legend->AddEntry(nOut, "outside LH2", "l");
  legend->Draw();
  canvas.cd(2);
  cat->SetTitle("extra track species (stored pid);;extra tracks");
  for (size_t c = 0; c < kCategories.size(); ++c) cat->GetXaxis()->SetBinLabel(c + 1, kCategories[c]);
  cat->Draw("HIST");
  canvas.cd(3);
  Subtract(sum.mom[0], sum.mom[1], "extra track momentum;mom0 [GeV/c];extra tracks")->Draw("HIST");
  canvas.cd(4);
  Subtract(sum.dcaBeam[0], sum.dcaBeam[1], "DCA to the beam track;DCA [mm];extra tracks")->Draw("HIST");
  canvas.cd(5);
  Subtract(sum.dcaLambda[0], sum.dcaLambda[1], "DCA to the #Lambda daughters (min);DCA [mm];extra tracks")->Draw("HIST");
  canvas.cd(6);
  Subtract(sum.vtz[0], sum.vtz[1], "z of the closest point to the beam;z [mm];extra tracks")->Draw("HIST");
  canvas.Print(pdf);

  canvas.Clear();
  canvas.Divide(2, 2);
  canvas.cd(1);
  sum.dedx[0]->SetTitle("extra tracks, signal window;q #times mom0 [GeV/c];dE/dx");
  sum.dedx[0]->Draw("COLZ");
  canvas.cd(2);
  Subtract(sum.chisqr[0], sum.chisqr[1], "extra track #chi^{2};#chi^{2};extra tracks")->Draw("HIST");
  canvas.cd(3);
  Subtract(sum.nhit[0], sum.nhit[1], "extra track hits;hits;extra tracks")->Draw("HIST");
  canvas.cd(4);
  Subtract(sum.firstZ[0], sum.firstZ[1], "most upstream hit z;z [mm];extra tracks")->Draw("HIST");
  canvas.Print(pdf);
  canvas.Print(pdf + "]");
  cout << "Wrote " << pdf << endl;
}
