#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TDatime.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TTree.h>

#include "../DstTPCBranches.hh"

using namespace std;

const vector<int> runnumbers = {2447, 2449, 2450, 2451, 2452, 2453, 2454, 2456, 2457, 2458, 2459, 2460};
//const vector<int> runnumbers = {2447};

const double lambda_peak_min = 1.095; // GeV/c2
const double lambda_peak_max = 1.125; // GeV/c2
const double tight_close_dist = 5.;   // mm
const double tight_chisqr = 2.;
const int lambda_bkg_poly_order = 1; // choose 1, 2, or 3

const double mass_pi = 0.13957039; // GeV/c2
const double mass_p = 0.93827208;  // GeV/c2

bool IsPion(int pid) { return (pid & 0x1) != 0; }
bool IsProton(int pid) { return (pid & 0x4) != 0; }

void DrawPeakWindowLines(double ymax)
{
  TLine *lmin = new TLine(lambda_peak_min, 0., lambda_peak_min, ymax);
  TLine *lmax = new TLine(lambda_peak_max, 0., lambda_peak_max, ymax);
  lmin->SetLineColor(kRed + 1);
  lmax->SetLineColor(kRed + 1);
  lmin->SetLineStyle(2);
  lmax->SetLineStyle(2);
  lmin->SetLineWidth(2);
  lmax->SetLineWidth(2);
  lmin->Draw("same");
  lmax->Draw("same");
}

void DrawOverlay(TH1D *hall, TH1D *hsel, const char *all_label, const char *sel_label)
{
  hall->SetLineColor(kBlack);
  hsel->SetLineColor(kRed + 1);
  hsel->SetLineWidth(2);
  hall->Draw("hist");
  hsel->Draw("hist same");

  TLegend *leg = new TLegend(0.55, 0.72, 0.80, 0.88);
  leg->AddEntry(hall, all_label, "l");
  leg->AddEntry(hsel, sel_label, "l");
  leg->Draw();
}

struct YieldResult {
  double mean = 0.;
  double sigma = 0.;
  double signal = 0.;
  double background = 0.;
  double total = 0.;
};

YieldResult FitLambdaYield(TH1D *hist, const char *tag)
{
  YieldResult result;
  if (!hist || hist->GetEntries() <= 0)
    return result;

  const double fit_min = 1.08;
  const double fit_max = 1.14;
  const int bkg_order = max(1, min(3, lambda_bkg_poly_order));

  string bkg_expr = "[3]";
  string bkg_only_expr = "[0]";
  for (int i = 1; i <= bkg_order; ++i) {
    bkg_expr += Form(" + [%d]*TMath::Power(x,%d)", 3 + i, i);
    bkg_only_expr += Form(" + [%d]*TMath::Power(x,%d)", i, i);
  }

  TF1 *fit = new TF1(Form("fit_%s", tag),
                     ("[0]*TMath::Gaus(x,[1],[2],0) + " + bkg_expr).c_str(),
                     fit_min, fit_max);
  fit->SetParName(0, "Amp");
  fit->SetParName(1, "Mean");
  fit->SetParName(2, "Sigma");
  fit->SetParameter(0, hist->GetMaximum());
  fit->SetParameter(1, 1.115);
  fit->SetParameter(2, 0.004);
  for (int i = 0; i <= bkg_order; ++i) {
    fit->SetParName(3 + i, Form("p%d", i));
    fit->SetParameter(3 + i, (i == 0) ? 100. : 0.);
  }
  fit->SetParLimits(1, 1.108, 1.122);
  fit->SetParLimits(2, 0.001, 0.015);

  hist->Fit(fit, "R");

  TF1 *bkg = new TF1(Form("bkg_%s", tag), bkg_only_expr.c_str(), fit_min, fit_max);
  for (int i = 0; i <= bkg_order; ++i)
    bkg->SetParameter(i, fit->GetParameter(3 + i));
  bkg->SetLineColor(kGreen + 2);
  bkg->SetLineStyle(2);
  bkg->Draw("same");

  TF1 *sig = new TF1(Form("sig_%s", tag), "[0]*TMath::Gaus(x,[1],[2],0)", fit_min, fit_max);
  sig->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2));
  sig->SetLineColor(kRed);
  sig->Draw("same");

  fit->SetLineColor(kMagenta);
  fit->Draw("same");

  result.mean = fit->GetParameter(1);
  result.sigma = fit->GetParameter(2);
  const double yield_min = result.mean - 3.0 * result.sigma;
  const double yield_max = result.mean + 3.0 * result.sigma;
  const double binw = hist->GetBinWidth(1);

  result.signal = sig->Integral(yield_min, yield_max) / binw;
  result.background = bkg->Integral(yield_min, yield_max) / binw;
  result.total = fit->Integral(yield_min, yield_max) / binw;

  TLegend *leg = new TLegend(0.55, 0.60, 0.80, 0.82);
  leg->AddEntry(hist, "Data", "l");
  leg->AddEntry(fit, Form("Gaussian + pol%d", bkg_order), "l");
  leg->AddEntry(sig, "Lambda peak", "l");
  leg->AddEntry(bkg, "Background", "l");
  leg->Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.55, 0.54, Form("Yield = %.0f", result.signal));
  lat.DrawLatex(0.55, 0.49, Form("#mu = %.5f", result.mean));
  lat.DrawLatex(0.55, 0.44, Form("#sigma = %.5f", result.sigma));
  lat.DrawLatex(0.55, 0.39, Form("Bkg = %.0f", result.background));

  cout << tag << " Background polynomial order = " << bkg_order << endl;
  cout << tag << " Lambda mean  = " << result.mean << endl;
  cout << tag << " Lambda sigma = " << result.sigma << endl;
  cout << tag << " Signal yield = " << result.signal << endl;
  cout << tag << " Bkg yield    = " << result.background << endl;
  cout << tag << " Total yield  = " << result.total << endl;

  return result;
}


void cut_condition()
{
  string dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root/physics-735";
  string outpdf;
  if (runnumbers.size() == 1)
    outpdf = Form("result/cut-condition-run%05d.pdf", runnumbers.front());
  else
    outpdf = Form("result/cut-condition-run%05d-%05d-n%zu.pdf",
                  runnumbers.front(), runnumbers.back(), runnumbers.size());

  TCanvas *c1 = new TCanvas("c1", "c1", 900, 700);
  TPaveText *p = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
  p->AddText("cut-condition.cc");
  p->AddText("Lambda cut study from p + pi^{-}");
  if (runnumbers.size() == 1)
    p->AddText(Form("run%05d", runnumbers.front()));
  else
    p->AddText(Form("runs %05d-%05d, n=%zu", runnumbers.front(), runnumbers.back(), runnumbers.size()));
  p->AddText("comparison: all p#pi^{-}, Lambda peak, and tight quality cut");
  p->AddText(Form("peak window: %.3f < M(p#pi^{-}) < %.3f GeV/c^{2}",
                  lambda_peak_min, lambda_peak_max));
  p->AddText(Form("tight cut: closeDist < %.1f mm, #chi^{2}_{p} <= %.1f, #chi^{2}_{#pi} <= %.1f",
                  tight_close_dist, tight_chisqr, tight_chisqr));
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",
                  now.GetYear(), now.GetMonth(), now.GetDay(),
                  now.GetHour(), now.GetMinute(), now.GetSecond()));
  p->Draw();
  c1->Print((outpdf + "(").c_str());
  c1->Clear();

  TH1D *hist_mass = new TH1D("hist_mass",
                             "p#pi^{-} invariant mass;M(p#pi^{-}) [GeV/c^{2}];Counts",
                             200, 1.05, 1.25);
  TH1D *hist_mass_tight = new TH1D("hist_mass_tight",
                                   "p#pi^{-} invariant mass, tight cut;M(p#pi^{-}) [GeV/c^{2}];Counts",
                                   200, 1.05, 1.25);

  TH1D *hist_dst_all = new TH1D("hist_dst_all",
                                "p#pi^{-} close distance;closeDistTpc [mm];Counts",
                                120, 0, 30);
  TH1D *hist_dst_peak = new TH1D("hist_dst_peak",
                                 "p#pi^{-} close distance;closeDistTpc [mm];Counts",
                                 120, 0, 30);
  TH1D *hist_dst_tight = new TH1D("hist_dst_tight",
                                  "p#pi^{-} close distance, tight cut;closeDistTpc [mm];Counts",
                                  120, 0, 30);

  TH1D *hist_pid = new TH1D("hist_pid", "track PID bit pattern;pid;Counts", 10, -0.5, 9.5);

  TH1D *hist_chisqr_p_all = new TH1D("hist_chisqr_p_all",
                                     "proton track #chi^{2};#chi^{2};Counts",
                                     120, 0, 10);
  TH1D *hist_chisqr_p_peak = new TH1D("hist_chisqr_p_peak",
                                      "proton track #chi^{2};#chi^{2};Counts",
                                      120, 0, 10);
  TH1D *hist_chisqr_p_tight = new TH1D("hist_chisqr_p_tight",
                                       "proton track #chi^{2}, tight cut;#chi^{2};Counts",
                                       120, 0, 10);
  TH1D *hist_chisqr_pi_all = new TH1D("hist_chisqr_pi_all",
                                      "#pi^{-} track #chi^{2};#chi^{2};Counts",
                                      120, 0, 10);
  TH1D *hist_chisqr_pi_peak = new TH1D("hist_chisqr_pi_peak",
                                       "#pi^{-} track #chi^{2};#chi^{2};Counts",
                                       120, 0, 10);
  TH1D *hist_chisqr_pi_tight = new TH1D("hist_chisqr_pi_tight",
                                        "#pi^{-} track #chi^{2}, tight cut;#chi^{2};Counts",
                                        120, 0, 10);
  TH2D *hist_chisqr_pair_all = new TH2D("hist_chisqr_pair_all",
                                        "track #chi^{2}, all p#pi^{-};proton #chi^{2};#pi^{-} #chi^{2}",
                                        120, 0, 10, 120, 0, 10);
  TH2D *hist_chisqr_pair_peak = new TH2D("hist_chisqr_pair_peak",
                                         "track #chi^{2}, Lambda peak;proton #chi^{2};#pi^{-} #chi^{2}",
                                         120, 0, 10, 120, 0, 10);
  TH2D *hist_chisqr_pair_tight = new TH2D("hist_chisqr_pair_tight",
                                          "track #chi^{2}, tight cut;proton #chi^{2};#pi^{-} #chi^{2}",
                                          120, 0, 10, 120, 0, 10);

  TH2D *hist_vtx_xz_all = new TH2D("hist_vtx_xz_all",
                                   "vertex, all p#pi^{-};vtx Z [mm];vtx X [mm]",
                                   240, -300, 300, 240, -300, 300);
  TH2D *hist_vtx_yz_all = new TH2D("hist_vtx_yz_all",
                                   "vertex, all p#pi^{-};vtx Z [mm];vtx Y [mm]",
                                   240, -300, 300, 240, -300, 300);
  TH2D *hist_vtx_xy_all = new TH2D("hist_vtx_xy_all",
                                   "vertex, all p#pi^{-};vtx X [mm];vtx Y [mm]",
                                   240, -300, 300, 240, -300, 300);
  TH2D *hist_vtx_xz_peak = new TH2D("hist_vtx_xz_peak",
                                    "vertex, Lambda peak;vtx Z [mm];vtx X [mm]",
                                    240, -300, 300, 240, -300, 300);
  TH2D *hist_vtx_yz_peak = new TH2D("hist_vtx_yz_peak",
                                    "vertex, Lambda peak;vtx Z [mm];vtx Y [mm]",
                                    240, -300, 300, 240, -300, 300);
  TH2D *hist_vtx_xy_peak = new TH2D("hist_vtx_xy_peak",
                                    "vertex, Lambda peak;vtx X [mm];vtx Y [mm]",
                                    240, -300, 300, 240, -300, 300);
  TH2D *hist_vtx_xz_tight = new TH2D("hist_vtx_xz_tight",
                                     "vertex, tight cut;vtx Z [mm];vtx X [mm]",
                                     240, -300, 300, 240, -300, 300);
  TH2D *hist_vtx_yz_tight = new TH2D("hist_vtx_yz_tight",
                                     "vertex, tight cut;vtx Z [mm];vtx Y [mm]",
                                     240, -300, 300, 240, -300, 300);
  TH2D *hist_vtx_xy_tight = new TH2D("hist_vtx_xy_tight",
                                     "vertex, tight cut;vtx X [mm];vtx Y [mm]",
                                     240, -300, 300, 240, -300, 300);

  TH1D *hist_p_mass_peak = new TH1D("hist_p_mass_peak",
                                    "proton candidate invariant mass, Lambda peak;m_{p} [GeV/c^{2}];Counts",
                                    120, 0.90, 0.98);
  TH1D *hist_pi_mass_peak = new TH1D("hist_pi_mass_peak",
                                     "#pi^{-} candidate invariant mass, Lambda peak;m_{#pi} [GeV/c^{2}];Counts",
                                     120, 0.10, 0.18);
  TH1D *hist_p_mom_peak = new TH1D("hist_p_mom_peak",
                                   "proton momentum at mom0, Lambda peak;p_{p} [GeV/c];Counts",
                                   160, 0, 2.0);
  TH1D *hist_pi_mom_peak = new TH1D("hist_pi_mom_peak",
                                    "#pi^{-} momentum at mom0, Lambda peak;p_{#pi} [GeV/c];Counts",
                                    160, 0, 2.0);

  TH1D *hist_p_mass_tight = new TH1D("hist_p_mass_tight",
                                     "proton candidate invariant mass, tight cut;m_{p} [GeV/c^{2}];Counts",
                                     120, 0.90, 0.98);
  TH1D *hist_pi_mass_tight = new TH1D("hist_pi_mass_tight",
                                      "#pi^{-} candidate invariant mass, tight cut;m_{#pi} [GeV/c^{2}];Counts",
                                      120, 0.10, 0.18);
  TH1D *hist_p_mom_tight = new TH1D("hist_p_mom_tight",
                                    "proton momentum at mom0, tight cut;p_{p} [GeV/c];Counts",
                                    160, 0, 2.0);
  TH1D *hist_pi_mom_tight = new TH1D("hist_pi_mom_tight",
                                     "#pi^{-} momentum at mom0, tight cut;p_{#pi} [GeV/c];Counts",
                                     160, 0, 2.0);

  TH2D *hist_dedx_vs_signed_mom_peak = new TH2D("hist_dedx_vs_signed_mom_peak",
                                                "dE/dx vs signed momentum, Lambda peak;#it{p/z} [GeV/#it{c}];TPC #LT#it{dE/dx}#GT (a.u.)",
                                                240, -2.0, 2.0, 240, 0, 500);
  TH2D *hist_dedx_vs_signed_mom_p_peak = new TH2D("hist_dedx_vs_signed_mom_p_peak",
                                                  "proton dE/dx vs signed momentum, Lambda peak;#it{p/z} [GeV/#it{c}];TPC #LT#it{dE/dx}#GT (a.u.)",
                                                  240, -2.0, 2.0, 240, 0, 500);
  TH2D *hist_dedx_vs_signed_mom_pi_peak = new TH2D("hist_dedx_vs_signed_mom_pi_peak",
                                                   "#pi^{-} dE/dx vs signed momentum, Lambda peak;#it{p/z} [GeV/#it{c}];TPC #LT#it{dE/dx}#GT (a.u.)",
                                                   240, -2.0, 2.0, 240, 0, 500);
  TH2D *hist_charge_vs_mom_peak = new TH2D("hist_charge_vs_mom_peak",
                                           "charge vs momentum, Lambda peak;p_{mom0} [GeV/c];charge",
                                           160, 0, 2.0, 5, -2.5, 2.5);
  TH2D *hist_dedx_vs_signed_mom_tight = new TH2D("hist_dedx_vs_signed_mom_tight",
                                                 "dE/dx vs signed momentum, tight cut;#it{p/z} [GeV/#it{c}];TPC #LT#it{dE/dx}#GT (a.u.)",
                                                 240, -2.0, 2.0, 240, 0, 500);
  TH2D *hist_dedx_vs_signed_mom_p_tight = new TH2D("hist_dedx_vs_signed_mom_p_tight",
                                                   "proton dE/dx vs signed momentum, tight cut;#it{p/z} [GeV/#it{c}];TPC #LT#it{dE/dx}#GT (a.u.)",
                                                   240, -2.0, 2.0, 240, 0, 500);
  TH2D *hist_dedx_vs_signed_mom_pi_tight = new TH2D("hist_dedx_vs_signed_mom_pi_tight",
                                                    "#pi^{-} dE/dx vs signed momentum, tight cut;#it{p/z} [GeV/#it{c}];TPC #LT#it{dE/dx}#GT (a.u.)",
                                                    240, -2.0, 2.0, 240, 0, 500);
  TH2D *hist_charge_vs_mom_tight = new TH2D("hist_charge_vs_mom_tight",
                                            "charge vs momentum, tight cut;p_{mom0} [GeV/c];charge",
                                            160, 0, 2.0, 5, -2.5, 2.5);

  Long64_t n_ppi_cand = 0;
  Long64_t n_peak_cand = 0;
  Long64_t n_tight_cand = 0;
  for (int runnumber : runnumbers) {
    TFile *file = new TFile(Form("%s/run%05d_DstTPCHelixTracking.root",
                                 dir.c_str(), runnumber));
    if (!file || file->IsZombie()) {
      cerr << "Cannot open input file for run " << runnumber << endl;
      if (file)
        file->Close();
      continue;
    }

    TTree *tree = GetTree(file, "tpc");
    if (!tree) {
      cerr << "No TTree found for run " << runnumber << endl;
      file->Close();
      continue;
    }

    dst_tpc::HelixTracks tracks;
    dst_tpc::HelixPairs pairs;
    bool ok = true;
    ok &= tracks.SetBranchAddresses(tree);
    ok &= pairs.SetBranchAddresses(tree);
    if (!ok) {
      cerr << "Some helix track/pair branches were not connected for run " << runnumber << endl;
      file->Close();
      continue;
    }

    cout << "Processing run " << runnumber << ", entries = " << tree->GetEntries() << endl;
    const Long64_t n_entries = tree->GetEntries();
    for (Long64_t i = 0; i < n_entries; ++i) {
      tree->GetEntry(i);

    for (Int_t ip = 0; ip < tracks.ntTpc; ++ip) {
      const int pid_p = tracks.pid->at(ip);
      const int charge_p = tracks.charge->at(ip);
      hist_pid->Fill(pid_p);

      if (!IsProton(pid_p) || charge_p <= 0)
        continue;

      for (Int_t ipi = 0; ipi < tracks.ntTpc; ++ipi) {
        if (ip == ipi)
          continue;

      const int pid_pi = tracks.pid->at(ipi);
        const int charge_pi = tracks.charge->at(ipi);
        if (!IsPion(pid_pi) || charge_pi >= 0)
          continue;

        const double close_dist = pairs.closeDistTpc->at(ip).at(ipi);
        if (!std::isfinite(close_dist))
          continue;

        const double vtx_x = pairs.vtxTpc->at(ip).at(ipi);
        const double vtx_y = pairs.vtyTpc->at(ip).at(ipi);
        const double vtx_z = pairs.vtzTpc->at(ip).at(ipi);
        if (!std::isfinite(vtx_x) || !std::isfinite(vtx_y) || !std::isfinite(vtx_z))
          continue;

        const double p_px = pairs.mom_vtx->at(ip).at(ipi);
        const double p_py = pairs.mom_vty->at(ip).at(ipi);
        const double p_pz = pairs.mom_vtz->at(ip).at(ipi);
        const double pi_px = pairs.mom_vtx->at(ipi).at(ip);
        const double pi_py = pairs.mom_vty->at(ipi).at(ip);
        const double pi_pz = pairs.mom_vtz->at(ipi).at(ip);
        if (!std::isfinite(p_px) || !std::isfinite(p_py) || !std::isfinite(p_pz) ||
            !std::isfinite(pi_px) || !std::isfinite(pi_py) || !std::isfinite(pi_pz))
          continue;

        TLorentzVector lv_p;
        TLorentzVector lv_pi;
        lv_p.SetXYZM(p_px, p_py, p_pz, mass_p);
        lv_pi.SetXYZM(pi_px, pi_py, pi_pz, mass_pi);

        const double inv_mass = (lv_p + lv_pi).M();
        const double chisqr_p = tracks.chisqr->at(ip);
        const double chisqr_pi = tracks.chisqr->at(ipi);
        const double mom0_p = tracks.mom0->at(ip);
        const double mom0_pi = tracks.mom0->at(ipi);
        ++n_ppi_cand;

        hist_mass->Fill(inv_mass);
        hist_dst_all->Fill(close_dist);
        hist_chisqr_p_all->Fill(chisqr_p);
        hist_chisqr_pi_all->Fill(chisqr_pi);
        hist_chisqr_pair_all->Fill(chisqr_p, chisqr_pi);
        hist_vtx_xz_all->Fill(vtx_z, vtx_x);
        hist_vtx_yz_all->Fill(vtx_z, vtx_y);
        hist_vtx_xy_all->Fill(vtx_x, vtx_y);

        const bool is_peak = (lambda_peak_min <= inv_mass && inv_mass <= lambda_peak_max);
        const bool is_tight = (close_dist < tight_close_dist &&
                               chisqr_p <= tight_chisqr && chisqr_pi <= tight_chisqr);

        if (is_peak) {
          hist_dst_peak->Fill(close_dist);
          hist_chisqr_p_peak->Fill(chisqr_p);
          hist_chisqr_pi_peak->Fill(chisqr_pi);
          hist_chisqr_pair_peak->Fill(chisqr_p, chisqr_pi);
          hist_vtx_xz_peak->Fill(vtx_z, vtx_x);
          hist_vtx_yz_peak->Fill(vtx_z, vtx_y);
          hist_vtx_xy_peak->Fill(vtx_x, vtx_y);
          hist_p_mass_peak->Fill(lv_p.M());
          hist_pi_mass_peak->Fill(lv_pi.M());
          hist_p_mom_peak->Fill(mom0_p);
          hist_pi_mom_peak->Fill(mom0_pi);
          hist_dedx_vs_signed_mom_peak->Fill(charge_p * mom0_p, tracks.dEdx->at(ip));
          hist_dedx_vs_signed_mom_peak->Fill(charge_pi * mom0_pi, tracks.dEdx->at(ipi));
          hist_dedx_vs_signed_mom_p_peak->Fill(charge_p * mom0_p, tracks.dEdx->at(ip));
          hist_dedx_vs_signed_mom_pi_peak->Fill(charge_pi * mom0_pi, tracks.dEdx->at(ipi));
          hist_charge_vs_mom_peak->Fill(mom0_p, charge_p);
          hist_charge_vs_mom_peak->Fill(mom0_pi, charge_pi);
          ++n_peak_cand;
        }

        if (is_tight) {
          hist_mass_tight->Fill(inv_mass);
          hist_dst_tight->Fill(close_dist);
          hist_chisqr_p_tight->Fill(chisqr_p);
          hist_chisqr_pi_tight->Fill(chisqr_pi);
          hist_chisqr_pair_tight->Fill(chisqr_p, chisqr_pi);
          hist_vtx_xz_tight->Fill(vtx_z, vtx_x);
          hist_vtx_yz_tight->Fill(vtx_z, vtx_y);
          hist_vtx_xy_tight->Fill(vtx_x, vtx_y);
          hist_p_mass_tight->Fill(lv_p.M());
          hist_pi_mass_tight->Fill(lv_pi.M());
          hist_p_mom_tight->Fill(mom0_p);
          hist_pi_mom_tight->Fill(mom0_pi);
          hist_dedx_vs_signed_mom_tight->Fill(charge_p * mom0_p, tracks.dEdx->at(ip));
          hist_dedx_vs_signed_mom_tight->Fill(charge_pi * mom0_pi, tracks.dEdx->at(ipi));
          hist_dedx_vs_signed_mom_p_tight->Fill(charge_p * mom0_p, tracks.dEdx->at(ip));
          hist_dedx_vs_signed_mom_pi_tight->Fill(charge_pi * mom0_pi, tracks.dEdx->at(ipi));
          hist_charge_vs_mom_tight->Fill(mom0_p, charge_p);
          hist_charge_vs_mom_tight->Fill(mom0_pi, charge_pi);
          ++n_tight_cand;
        }
      }
    }
    }
    file->Close();
  }

  cout << "p pi candidates = " << n_ppi_cand << endl;
  cout << "Lambda peak candidates = " << n_peak_cand << endl;
  cout << "tight cut candidates = " << n_tight_cand << endl;

  hist_mass->SetLineColor(kBlack);
  hist_mass->Draw("hist");
  DrawPeakWindowLines(hist_mass->GetMaximum() * 1.05);
  TLegend *leg_mass = new TLegend(0.55, 0.72, 0.80, 0.88);
  leg_mass->AddEntry(hist_mass, "all p#pi^{-}", "l");
  leg_mass->AddEntry((TObject*)0, Form("peak: %.3f-%.3f", lambda_peak_min, lambda_peak_max), "");
  leg_mass->Draw();
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_mass->SetLineColor(kBlack);
  hist_mass_tight->SetLineColor(kBlue + 1);
  hist_mass_tight->SetLineWidth(2);
  hist_mass->Draw("hist");
  hist_mass_tight->Draw("hist same");
  DrawPeakWindowLines(hist_mass->GetMaximum() * 1.05);
  TLegend *leg_mass_tight = new TLegend(0.50, 0.68, 0.80, 0.88);
  leg_mass_tight->AddEntry(hist_mass, "all p#pi^{-}", "l");
  leg_mass_tight->AddEntry(hist_mass_tight, Form("closeDist < %.1f mm, #chi^{2} <= %.1f", tight_close_dist, tight_chisqr), "l");
  leg_mass_tight->AddEntry((TObject*)0, Form("peak: %.3f-%.3f", lambda_peak_min, lambda_peak_max), "");
  leg_mass_tight->Draw();
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_mass_tight->Draw("hist");
  DrawPeakWindowLines(hist_mass_tight->GetMaximum() * 1.05);
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_mass->SetLineColor(kBlue);
  hist_mass->Draw("hist");
  FitLambdaYield(hist_mass, "all_ppi");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_mass_tight->SetLineColor(kBlue);
  hist_mass_tight->Draw("hist");
  FitLambdaYield(hist_mass_tight, "tight_cut");
  c1->Print(outpdf.c_str());
  c1->Clear();

  DrawOverlay(hist_dst_all, hist_dst_peak, "all p#pi^{-}", "Lambda peak selected");
  c1->Print(outpdf.c_str());
  c1->Clear();

  DrawOverlay(hist_dst_all, hist_dst_tight, "all p#pi^{-}", "tight cut");
  c1->Print(outpdf.c_str());
  c1->Clear();

  DrawOverlay(hist_chisqr_p_all, hist_chisqr_p_peak, "all p#pi^{-}", "Lambda peak selected");
  c1->Print(outpdf.c_str());
  c1->Clear();

  DrawOverlay(hist_chisqr_p_all, hist_chisqr_p_tight, "all p#pi^{-}", "tight cut");
  c1->Print(outpdf.c_str());
  c1->Clear();

  DrawOverlay(hist_chisqr_pi_all, hist_chisqr_pi_peak, "all p#pi^{-}", "Lambda peak selected");
  c1->Print(outpdf.c_str());
  c1->Clear();

  DrawOverlay(hist_chisqr_pi_all, hist_chisqr_pi_tight, "all p#pi^{-}", "tight cut");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_chisqr_pair_all->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_chisqr_pair_peak->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_chisqr_pair_tight->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_pid->Draw("hist");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_vtx_xz_all->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_vtx_yz_all->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_vtx_xy_all->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_vtx_xz_peak->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_vtx_yz_peak->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_vtx_xy_peak->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_vtx_xz_tight->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_vtx_yz_tight->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_vtx_xy_tight->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_dedx_vs_signed_mom_peak->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_dedx_vs_signed_mom_p_peak->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_dedx_vs_signed_mom_pi_peak->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_charge_vs_mom_peak->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_dedx_vs_signed_mom_tight->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_dedx_vs_signed_mom_p_tight->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_dedx_vs_signed_mom_pi_tight->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_charge_vs_mom_tight->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_p_mass_peak->SetLineColor(kBlue + 1);
  hist_p_mass_peak->Draw("hist");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_pi_mass_peak->SetLineColor(kRed + 1);
  hist_pi_mass_peak->Draw("hist");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_p_mom_peak->SetLineColor(kBlue + 1);
  hist_p_mom_peak->Draw("hist");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_pi_mom_peak->SetLineColor(kRed + 1);
  hist_pi_mom_peak->Draw("hist");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_p_mass_tight->SetLineColor(kBlue + 1);
  hist_p_mass_tight->Draw("hist");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_pi_mass_tight->SetLineColor(kRed + 1);
  hist_pi_mass_tight->Draw("hist");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_p_mom_tight->SetLineColor(kBlue + 1);
  hist_p_mom_tight->Draw("hist");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_pi_mom_tight->SetLineColor(kRed + 1);
  hist_pi_mom_tight->Draw("hist");
  c1->Print((outpdf + ")").c_str());

}
