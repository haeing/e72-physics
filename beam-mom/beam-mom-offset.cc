#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TInterpreter.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TParticle.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVector3.h"

namespace {
constexpr double kRadToDeg = 180.0 / TMath::Pi();
constexpr double kDegToRad = TMath::Pi() / 180.0;

struct BeamOffset {
  double rot_x_deg = 0.0;      // Delta atan2(py,pz), rotation around x axis.
  double rot_y_deg = 0.0;      // Delta atan2(px,pz), rotation around y axis.
  double rot_z_deg = 0.0;      // Delta atan2(py,px), rotation around z axis.
  double mom_scale = 1.0;      // Mean |p_TGT| / |p_BEAM|.
};

const TParticle* FindBeamParticle(const std::vector<TParticle>* particles)
{
  if (!particles || particles->empty()) {
    return nullptr;
  }

  const auto kaon = std::find_if(particles->begin(), particles->end(),
                                 [](const TParticle& p) {
                                   return p.GetPdgCode() == -321;
                                 });
  return kaon != particles->end() ? &(*kaon) : &particles->front();
}

TVector3 Momentum(const TParticle& p)
{
  return TVector3(p.Px(), p.Py(), p.Pz());
}

double WrapDeg(double angle)
{
  while (angle > 180.0) {
    angle -= 360.0;
  }
  while (angle <= -180.0) {
    angle += 360.0;
  }
  return angle;
}

double RotationAngleDeg(const TVector3& p, char axis)
{
  if (axis == 'x') {
    return std::atan2(p.Y(), p.Z()) * kRadToDeg;
  }
  if (axis == 'y') {
    return std::atan2(p.X(), p.Z()) * kRadToDeg;
  }
  return std::atan2(p.Y(), p.X()) * kRadToDeg;
}

double AxisDirectionAngleDeg(const TVector3& p, char axis)
{
  const double mag = p.Mag();
  if (mag <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double component = p.Z();
  if (axis == 'x') {
    component = p.X();
  } else if (axis == 'y') {
    component = p.Y();
  }

  const double cos_angle = std::clamp(component / mag, -1.0, 1.0);
  return std::acos(cos_angle) * kRadToDeg;
}

void StyleHist(TH1D& h, int color)
{
  h.SetLineColor(color);
  h.SetLineWidth(2);
}

void DrawPair(TCanvas& canvas, TH1D& h_before, TH1D& h_after,
              const char* pdf_name, const char* legend_before,
              const char* legend_after)
{
  canvas.Clear();
  const double ymax = 1.15 * std::max(h_before.GetMaximum(), h_after.GetMaximum());
  h_before.SetMaximum(ymax > 0.0 ? ymax : 1.0);
  h_before.Draw("hist");
  h_after.Draw("hist same");

  TLegend legend(0.58, 0.74, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.AddEntry(&h_before, legend_before, "l");
  legend.AddEntry(&h_after, legend_after, "l");
  legend.Draw();
  canvas.Print(pdf_name);
}
}  // namespace

// Apply the BEAM -> Target beam-momentum offset to any measured beam TVector3.
// The default mode applies x/y rotations and momentum scale.  The z rotation is
// available but disabled by default because it is poorly constrained for a beam
// nearly parallel to z; use full_xyz=true if you explicitly want it included.
TVector3 ApplyBeamMomentumOffset(const TVector3& p_beam,
                                 const BeamOffset& offset,
                                 bool full_xyz = false,
                                 bool apply_momentum_scale = true)
{
  TVector3 p = p_beam;
  p.RotateX(offset.rot_x_deg * kDegToRad);
  p.RotateY(offset.rot_y_deg * kDegToRad);
  if (full_xyz) {
    p.RotateZ(offset.rot_z_deg * kDegToRad);
  }
  if (apply_momentum_scale) {
    p *= offset.mom_scale;
  }
  return p;
}

void beam_mom_offset(
    const char* input_name =
        "/gpfs/group/had/sks/Users/haein/simul-data/e72-beam/e72_kbeam_735.root",
    const char* output_name = "beam-mom-offset.root",
    const char* pdf_name = "beam-mom-offset.pdf",
    Long64_t max_entries = -1)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1110);
  gInterpreter->GenerateDictionary("vector<TParticle>", "vector;TParticle.h");

  TFile input(input_name, "READ");
  if (input.IsZombie()) {
    std::cerr << "Failed to open input file: " << input_name << std::endl;
    return;
  }

  auto* tree = dynamic_cast<TTree*>(input.Get("g4hyptpc"));
  if (!tree) {
    std::cerr << "Cannot find TTree g4hyptpc in " << input_name << std::endl;
    return;
  }

  std::vector<TParticle>* beam = nullptr;
  std::vector<TParticle>* tgt = nullptr;
  tree->SetBranchAddress("BEAM", &beam);
  tree->SetBranchAddress("TGT", &tgt);

  TH1D h_dtheta_3d("h_dtheta_3d",
                   "BEAM to TGT 3D direction change;#Delta#theta_{3D} [deg];Events",
                   400, 0.0, 40.0);
  TH1D h_drot_x("h_drot_x",
                "Rotation around x axis: #Delta atan2(p_{y},p_{z});#Delta#theta_{x} [deg];Events",
                400, -20.0, 20.0);
  TH1D h_drot_y("h_drot_y",
                "Rotation around y axis: #Delta atan2(p_{x},p_{z});#Delta#theta_{y} [deg];Events",
                400, -20.0, 20.0);
  TH1D h_drot_z("h_drot_z",
                "Rotation around z axis: #Delta atan2(p_{y},p_{x});#Delta#theta_{z} [deg];Events",
                400, -180.0, 180.0);
  TH1D h_daxis_x("h_daxis_x",
                 "Change in angle to x axis;#Delta angle-to-x [deg];Events",
                 400, -20.0, 20.0);
  TH1D h_daxis_y("h_daxis_y",
                 "Change in angle to y axis;#Delta angle-to-y [deg];Events",
                 400, -20.0, 20.0);
  TH1D h_daxis_z("h_daxis_z",
                 "Change in angle to z axis;#Delta angle-to-z [deg];Events",
                 400, -20.0, 20.0);
  TH1D h_mom_beam("h_mom_beam", "Beam momentum at BEAM;|p| [same unit as file];Events",
                 400, 650.0, 800.0);
  TH1D h_mom_tgt("h_mom_tgt", "Beam momentum at TGT;|p| [same unit as file];Events",
                 400, 650.0, 800.0);
  TH1D h_mom_ratio("h_mom_ratio", "Momentum scale |p_{TGT}|/|p_{BEAM}|;scale;Events",
                   400, 0.90, 1.05);

  Long64_t nentries = tree->GetEntries();
  if (max_entries >= 0 && max_entries < nentries) {
    nentries = max_entries;
  }

  Long64_t used = 0;
  Long64_t skipped = 0;
  for (Long64_t entry = 0; entry < nentries; ++entry) {
    tree->GetEntry(entry);

    const TParticle* p_beam = FindBeamParticle(beam);
    const TParticle* p_tgt = FindBeamParticle(tgt);
    if (!p_beam || !p_tgt) {
      ++skipped;
      continue;
    }

    const TVector3 mom_beam = Momentum(*p_beam);
    const TVector3 mom_tgt = Momentum(*p_tgt);
    if (mom_beam.Mag() <= 0.0 || mom_tgt.Mag() <= 0.0) {
      ++skipped;
      continue;
    }

    h_dtheta_3d.Fill(mom_beam.Angle(mom_tgt) * kRadToDeg);
    h_drot_x.Fill(WrapDeg(RotationAngleDeg(mom_tgt, 'x') -
                          RotationAngleDeg(mom_beam, 'x')));
    h_drot_y.Fill(WrapDeg(RotationAngleDeg(mom_tgt, 'y') -
                          RotationAngleDeg(mom_beam, 'y')));
    h_drot_z.Fill(WrapDeg(RotationAngleDeg(mom_tgt, 'z') -
                          RotationAngleDeg(mom_beam, 'z')));
    h_daxis_x.Fill(AxisDirectionAngleDeg(mom_tgt, 'x') -
                   AxisDirectionAngleDeg(mom_beam, 'x'));
    h_daxis_y.Fill(AxisDirectionAngleDeg(mom_tgt, 'y') -
                   AxisDirectionAngleDeg(mom_beam, 'y'));
    h_daxis_z.Fill(AxisDirectionAngleDeg(mom_tgt, 'z') -
                   AxisDirectionAngleDeg(mom_beam, 'z'));
    h_mom_beam.Fill(mom_beam.Mag());
    h_mom_tgt.Fill(mom_tgt.Mag());
    h_mom_ratio.Fill(mom_tgt.Mag() / mom_beam.Mag());
    ++used;
  }

  BeamOffset offset;
  offset.rot_x_deg = h_drot_x.GetMean();
  offset.rot_y_deg = h_drot_y.GetMean();
  offset.rot_z_deg = h_drot_z.GetMean();
  offset.mom_scale = h_mom_ratio.GetMean();

  TH1D h_diff_px_raw("h_diff_px_raw", "TGT - BEAM p_{x};#Delta p_{x};Events", 400, -150.0, 150.0);
  TH1D h_diff_py_raw("h_diff_py_raw", "TGT - BEAM p_{y};#Delta p_{y};Events", 400, -80.0, 80.0);
  TH1D h_diff_pz_raw("h_diff_pz_raw", "TGT - BEAM p_{z};#Delta p_{z};Events", 400, -80.0, 80.0);
  TH1D h_diff_p_raw("h_diff_p_raw", "TGT - BEAM |p|;#Delta |p|;Events", 400, -80.0, 20.0);
  TH1D h_angle_raw("h_angle_raw", "Angle between TGT and BEAM;angle [deg];Events", 400, 0.0, 20.0);

  TH1D h_diff_px_cor("h_diff_px_cor", "TGT - corrected p_{x};#Delta p_{x};Events", 400, -150.0, 150.0);
  TH1D h_diff_py_cor("h_diff_py_cor", "TGT - corrected p_{y};#Delta p_{y};Events", 400, -80.0, 80.0);
  TH1D h_diff_pz_cor("h_diff_pz_cor", "TGT - corrected p_{z};#Delta p_{z};Events", 400, -80.0, 80.0);
  TH1D h_diff_p_cor("h_diff_p_cor", "TGT - corrected |p|;#Delta |p|;Events", 400, -80.0, 20.0);
  TH1D h_angle_cor("h_angle_cor", "Angle between TGT and corrected;angle [deg];Events", 400, 0.0, 20.0);

  TH2D h_px_actual_vs_cor("h_px_actual_vs_cor", "Target p_{x}: actual vs corrected;corrected p_{x};actual p_{x}",
                          400, -50.0, 150.0, 400, -50.0, 150.0);
  TH2D h_py_actual_vs_cor("h_py_actual_vs_cor", "Target p_{y}: actual vs corrected;corrected p_{y};actual p_{y}",
                          400, -50.0, 50.0, 400, -50.0, 50.0);
  TH2D h_pz_actual_vs_cor("h_pz_actual_vs_cor", "Target p_{z}: actual vs corrected;corrected p_{z};actual p_{z}",
                          400, 650.0, 780.0, 400, 650.0, 780.0);

  for (Long64_t entry = 0; entry < nentries; ++entry) {
    tree->GetEntry(entry);

    const TParticle* p_beam = FindBeamParticle(beam);
    const TParticle* p_tgt = FindBeamParticle(tgt);
    if (!p_beam || !p_tgt) {
      continue;
    }

    const TVector3 mom_beam = Momentum(*p_beam);
    const TVector3 mom_tgt = Momentum(*p_tgt);
    if (mom_beam.Mag() <= 0.0 || mom_tgt.Mag() <= 0.0) {
      continue;
    }

    const TVector3 mom_cor = ApplyBeamMomentumOffset(mom_beam, offset);

    h_diff_px_raw.Fill(mom_tgt.X() - mom_beam.X());
    h_diff_py_raw.Fill(mom_tgt.Y() - mom_beam.Y());
    h_diff_pz_raw.Fill(mom_tgt.Z() - mom_beam.Z());
    h_diff_p_raw.Fill(mom_tgt.Mag() - mom_beam.Mag());
    h_angle_raw.Fill(mom_tgt.Angle(mom_beam) * kRadToDeg);

    h_diff_px_cor.Fill(mom_tgt.X() - mom_cor.X());
    h_diff_py_cor.Fill(mom_tgt.Y() - mom_cor.Y());
    h_diff_pz_cor.Fill(mom_tgt.Z() - mom_cor.Z());
    h_diff_p_cor.Fill(mom_tgt.Mag() - mom_cor.Mag());
    h_angle_cor.Fill(mom_tgt.Angle(mom_cor) * kRadToDeg);

    h_px_actual_vs_cor.Fill(mom_cor.X(), mom_tgt.X());
    h_py_actual_vs_cor.Fill(mom_cor.Y(), mom_tgt.Y());
    h_pz_actual_vs_cor.Fill(mom_cor.Z(), mom_tgt.Z());
  }

  std::vector<TH1D*> hists1 = {
      &h_dtheta_3d, &h_drot_x, &h_drot_y, &h_drot_z,
      &h_daxis_x, &h_daxis_y, &h_daxis_z, &h_mom_beam,
      &h_mom_tgt, &h_mom_ratio,
      &h_diff_px_raw, &h_diff_py_raw, &h_diff_pz_raw, &h_diff_p_raw, &h_angle_raw,
      &h_diff_px_cor, &h_diff_py_cor, &h_diff_pz_cor, &h_diff_p_cor, &h_angle_cor};
  for (auto* h : hists1) {
    StyleHist(*h, kBlack);
  }
  StyleHist(h_diff_px_raw, kGray + 2);
  StyleHist(h_diff_py_raw, kGray + 2);
  StyleHist(h_diff_pz_raw, kGray + 2);
  StyleHist(h_diff_p_raw, kGray + 2);
  StyleHist(h_angle_raw, kGray + 2);
  StyleHist(h_diff_px_cor, kRed + 1);
  StyleHist(h_diff_py_cor, kRed + 1);
  StyleHist(h_diff_pz_cor, kRed + 1);
  StyleHist(h_diff_p_cor, kRed + 1);
  StyleHist(h_angle_cor, kRed + 1);

  TFile output(output_name, "RECREATE");
  for (auto* h : hists1) {
    h->Write();
  }
  h_px_actual_vs_cor.Write();
  h_py_actual_vs_cor.Write();
  h_pz_actual_vs_cor.Write();
  output.Close();

  TCanvas canvas("c_beam_mom_offset", "c_beam_mom_offset", 1100, 850);
  TPaveText title(0.10, 0.10, 0.90, 0.90, "NDC");
  title.AddText("beam-mom-offset.cc");
  title.AddText("BEAM momentum offset applied to Target");
  title.AddText(Form("Input: %s", input_name));
  title.AddText(Form("rot x/y/z = %.5f / %.5f / %.5f deg", offset.rot_x_deg, offset.rot_y_deg, offset.rot_z_deg));
  title.AddText(Form("momentum scale = %.8f", offset.mom_scale));
  TDatime now;
  title.AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",
                     now.GetYear(), now.GetMonth(), now.GetDay(),
                     now.GetHour(), now.GetMinute(), now.GetSecond()));
  title.Draw();
  canvas.Print((std::string(pdf_name) + "(").c_str());

  canvas.Clear();
  canvas.Divide(2, 2);
  canvas.cd(1); h_dtheta_3d.Draw("hist");
  canvas.cd(2); h_drot_x.Draw("hist");
  canvas.cd(3); h_drot_y.Draw("hist");
  canvas.cd(4); h_mom_ratio.Draw("hist");
  canvas.Print(pdf_name);

  DrawPair(canvas, h_diff_px_raw, h_diff_px_cor, pdf_name, "TGT - BEAM", "TGT - corrected");
  DrawPair(canvas, h_diff_py_raw, h_diff_py_cor, pdf_name, "TGT - BEAM", "TGT - corrected");
  DrawPair(canvas, h_diff_pz_raw, h_diff_pz_cor, pdf_name, "TGT - BEAM", "TGT - corrected");
  DrawPair(canvas, h_diff_p_raw, h_diff_p_cor, pdf_name, "TGT - BEAM", "TGT - corrected");
  DrawPair(canvas, h_angle_raw, h_angle_cor, pdf_name, "TGT vs BEAM", "TGT vs corrected");

  canvas.Clear();
  canvas.Divide(2, 2);
  canvas.cd(1); h_px_actual_vs_cor.Draw("colz");
  canvas.cd(2); h_py_actual_vs_cor.Draw("colz");
  canvas.cd(3); h_pz_actual_vs_cor.Draw("colz");
  canvas.Print((std::string(pdf_name) + ")").c_str());

  std::cout << "Input: " << input_name << '\n'
            << "Output ROOT: " << output_name << '\n'
            << "Output PDF: " << pdf_name << '\n'
            << "Entries read: " << nentries << '\n'
            << "Used events: " << used << '\n'
            << "Skipped events: " << skipped << '\n'
            << "Offset rot x/y/z [deg]: " << offset.rot_x_deg << " / "
            << offset.rot_y_deg << " / " << offset.rot_z_deg << '\n'
            << "Momentum scale |p_TGT|/|p_BEAM|: " << offset.mom_scale << '\n'
            << "Mean 3D direction change [deg]: " << h_dtheta_3d.GetMean()
            << " +- " << h_dtheta_3d.GetMeanError() << '\n'
            << "Mean angle after correction [deg]: " << h_angle_cor.GetMean()
            << " +- " << h_angle_cor.GetMeanError() << '\n'
            << "Mean TGT-corrected dpx/dpy/dpz: "
            << h_diff_px_cor.GetMean() << " / "
            << h_diff_py_cor.GetMean() << " / "
            << h_diff_pz_cor.GetMean() << std::endl;
}
