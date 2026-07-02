#include <cmath>
#include <iostream>
#include <string>

#include <TCanvas.h>
#include <TDatime.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TTree.h>

#include "../DstTPCBranches.hh"

using namespace std;

const int runnumber = 2447;

const double lambda_peak_min = 1.095; // GeV/c2
const double lambda_peak_max = 1.125; // GeV/c2

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

void DrawOverlay(TH1D *hall, TH1D *hpeak, const char *all_label, const char *peak_label)
{
  hall->SetLineColor(kBlack);
  hpeak->SetLineColor(kRed + 1);
  hpeak->SetLineWidth(2);
  hall->Draw("hist");
  hpeak->Draw("hist same");

  TLegend *leg = new TLegend(0.55, 0.72, 0.80, 0.88);
  leg->AddEntry(hall, all_label, "l");
  leg->AddEntry(hpeak, peak_label, "l");
  leg->Draw();
}

void cut_condition()
{  string dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root/physics-735";
  TFile *file = new TFile(Form("%s/run%05d_DstTPCHelixTracking.root",
                               dir.c_str(), runnumber));
  string outpdf = Form("result/cut-condition-run%05d.pdf", runnumber);

  if (!file || file->IsZombie()) {
    cerr << "Cannot open input file" << endl;
    return;
  }

  TTree *tree = GetTree(file, "tpc");
  if (!tree) {
    cerr << "No TTree found in input file" << endl;
    file->Close();
    return;
  }

  dst_tpc::HelixTracks tracks;
  dst_tpc::HelixPairs pairs;
  bool ok = true;
  ok &= tracks.SetBranchAddresses(tree);
  ok &= pairs.SetBranchAddresses(tree);
  if (!ok) {
    cerr << "Some helix track/pair branches were not connected" << endl;
    file->Close();
    return;
  }


  TCanvas *c1 = new TCanvas("c1", "c1", 900, 700);
  TPaveText *p = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
  p->AddText("cut-condition.cc");
  p->AddText("Lambda cut study from p + pi^{-}");
  p->AddText(Form("run%05d", runnumber));
  p->AddText("comparison: all p#pi^{-} candidates vs Lambda peak selected");
  p->AddText(Form("peak window: %.3f < M(p#pi^{-}) < %.3f GeV/c^{2}",
                  lambda_peak_min, lambda_peak_max));
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

  TH1D *hist_dst_all = new TH1D("hist_dst_all",
                                "p#pi^{-} close distance;closeDistTpc [mm];Counts",
                                120, 0, 30);
  TH1D *hist_dst_peak = new TH1D("hist_dst_peak",
                                 "p#pi^{-} close distance;closeDistTpc [mm];Counts",
                                 120, 0, 30);

  TH1D *hist_pid = new TH1D("hist_pid", "track PID bit pattern;pid;Counts", 10, -0.5, 9.5);

  TH1D *hist_chisqr_p_all = new TH1D("hist_chisqr_p_all",
                                     "proton track #chi^{2};#chi^{2};Counts",
                                     120, 0, 10);
  TH1D *hist_chisqr_p_peak = new TH1D("hist_chisqr_p_peak",
                                      "proton track #chi^{2};#chi^{2};Counts",
                                      120, 0, 10);
  TH1D *hist_chisqr_pi_all = new TH1D("hist_chisqr_pi_all",
                                      "#pi^{-} track #chi^{2};#chi^{2};Counts",
                                      120, 0, 10);
  TH1D *hist_chisqr_pi_peak = new TH1D("hist_chisqr_pi_peak",
                                       "#pi^{-} track #chi^{2};#chi^{2};Counts",
                                       120, 0, 10);
  TH2D *hist_chisqr_pair_all = new TH2D("hist_chisqr_pair_all",
                                        "track #chi^{2}, all p#pi^{-};proton #chi^{2};#pi^{-} #chi^{2}",
                                        120, 0, 10, 120, 0, 10);
  TH2D *hist_chisqr_pair_peak = new TH2D("hist_chisqr_pair_peak",
                                         "track #chi^{2}, Lambda peak;proton #chi^{2};#pi^{-} #chi^{2}",
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

  TH2D *hist_dedx_vs_signed_mom_peak = new TH2D("hist_dedx_vs_signed_mom_peak",
                                                "dE/dx vs signed momentum, Lambda peak;q p_{mom0} [GeV/c];dE/dx",
                                                240, -2.0, 2.0, 240, 0, 500);
  TH2D *hist_dedx_vs_signed_mom_p_peak = new TH2D("hist_dedx_vs_signed_mom_p_peak",
                                                  "proton dE/dx vs signed momentum, Lambda peak;q p_{mom0} [GeV/c];dE/dx",
                                                  240, -2.0, 2.0, 240, 0, 500);
  TH2D *hist_dedx_vs_signed_mom_pi_peak = new TH2D("hist_dedx_vs_signed_mom_pi_peak",
                                                   "#pi^{-} dE/dx vs signed momentum, Lambda peak;q p_{mom0} [GeV/c];dE/dx",
                                                   240, -2.0, 2.0, 240, 0, 500);
  TH2D *hist_charge_vs_mom_peak = new TH2D("hist_charge_vs_mom_peak",
                                           "charge vs momentum, Lambda peak;p_{mom0} [GeV/c];charge",
                                           160, 0, 2.0, 5, -2.5, 2.5);

  Long64_t n_ppi_cand = 0;
  Long64_t n_peak_cand = 0;
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
        ++n_ppi_cand;

        hist_mass->Fill(inv_mass);
        hist_dst_all->Fill(close_dist);
        hist_chisqr_p_all->Fill(tracks.chisqr->at(ip));
        hist_chisqr_pi_all->Fill(tracks.chisqr->at(ipi));
        hist_chisqr_pair_all->Fill(tracks.chisqr->at(ip), tracks.chisqr->at(ipi));
        hist_vtx_xz_all->Fill(vtx_z, vtx_x);
        hist_vtx_yz_all->Fill(vtx_z, vtx_y);
        hist_vtx_xy_all->Fill(vtx_x, vtx_y);
        if (inv_mass < lambda_peak_min || inv_mass > lambda_peak_max)
          continue;

        hist_dst_peak->Fill(close_dist);
        hist_chisqr_p_peak->Fill(tracks.chisqr->at(ip));
        hist_chisqr_pi_peak->Fill(tracks.chisqr->at(ipi));
        hist_chisqr_pair_peak->Fill(tracks.chisqr->at(ip), tracks.chisqr->at(ipi));
        hist_vtx_xz_peak->Fill(vtx_z, vtx_x);
        hist_vtx_yz_peak->Fill(vtx_z, vtx_y);
        hist_vtx_xy_peak->Fill(vtx_x, vtx_y);
        hist_p_mass_peak->Fill(lv_p.M());
        hist_pi_mass_peak->Fill(lv_pi.M());
        hist_p_mom_peak->Fill(tracks.mom0->at(ip));
        hist_pi_mom_peak->Fill(tracks.mom0->at(ipi));
        hist_dedx_vs_signed_mom_peak->Fill(charge_p * tracks.mom0->at(ip), tracks.dEdx->at(ip));
        hist_dedx_vs_signed_mom_peak->Fill(charge_pi * tracks.mom0->at(ipi), tracks.dEdx->at(ipi));
        hist_dedx_vs_signed_mom_p_peak->Fill(charge_p * tracks.mom0->at(ip), tracks.dEdx->at(ip));
        hist_dedx_vs_signed_mom_pi_peak->Fill(charge_pi * tracks.mom0->at(ipi), tracks.dEdx->at(ipi));
        hist_charge_vs_mom_peak->Fill(tracks.mom0->at(ip), charge_p);
        hist_charge_vs_mom_peak->Fill(tracks.mom0->at(ipi), charge_pi);
        ++n_peak_cand;
      }
    }
  }

  cout << "p pi candidates = " << n_ppi_cand << endl;
  cout << "Lambda peak candidates = " << n_peak_cand << endl;

  hist_mass->SetLineColor(kBlack);
  hist_mass->Draw("hist");
  DrawPeakWindowLines(hist_mass->GetMaximum() * 1.05);
  TLegend *leg_mass = new TLegend(0.55, 0.72, 0.80, 0.88);
  leg_mass->AddEntry(hist_mass, "all p#pi^{-}", "l");
  leg_mass->AddEntry((TObject*)0, Form("peak: %.3f-%.3f", lambda_peak_min, lambda_peak_max), "");
  leg_mass->Draw();
  c1->Print(outpdf.c_str());
  c1->Clear();

  DrawOverlay(hist_dst_all, hist_dst_peak, "all p#pi^{-}", "Lambda peak selected");
  c1->Print(outpdf.c_str());
  c1->Clear();

  DrawOverlay(hist_chisqr_p_all, hist_chisqr_p_peak, "all p#pi^{-}", "Lambda peak selected");
  c1->Print(outpdf.c_str());
  c1->Clear();

  DrawOverlay(hist_chisqr_pi_all, hist_chisqr_pi_peak, "all p#pi^{-}", "Lambda peak selected");
  c1->Print(outpdf.c_str());
  c1->Clear();


  hist_chisqr_pair_all->Draw("colz");
  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_chisqr_pair_peak->Draw("colz");
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
  c1->Print((outpdf + ")").c_str());

  file->Close();
}
