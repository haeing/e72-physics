// Estimate the K- beam energy loss before entering TargetLV.
//
// Run with CERN ROOT:
//   root -l -b -q 'beam_upstream_eloss.C()'
// or: root -l -b -q 'beam_upstream_eloss.C("input.root","output.root")'

#include <cmath>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TPaveText.h"
#include "TString.h"
#include "TTree.h"

namespace {
constexpr double kKaonMinusMassMeV = 493.677;

double KineticEnergyMeV(const TParticle& particle)
{
  const double p = particle.P();
  return std::sqrt(p*p + kKaonMinusMassMeV*kKaonMinusMassMeV)
         - kKaonMinusMassMeV;
}

const TParticle* FindPrimaryKaonAtTarget(const std::vector<TParticle>& hits)
{
  // VHitInfo stores the Geant4 track ID in TParticle daughter(0).
  for (const auto& hit : hits) {
    if (hit.GetPdgCode() == -321 && hit.GetDaughter(0) == 1)
      return &hit;
  }
  return nullptr;
}
} // namespace

void beam_upstream_eloss(
    const char* input = "/home/had/haein/simul-data/e72-k18ana/mom-735/e72_tpcon_beam_735.root",
    const char* output = "beam_upstream_eloss_735.root",
    const char* pdf = "beam_upstream_eloss_735.pdf")
{
  TFile inputFile(input, "READ");
  if (inputFile.IsZombie()) {
    std::cerr << "Cannot open input file: " << input << std::endl;
    return;
  }
  auto* inputTree = dynamic_cast<TTree*>(inputFile.Get("g4hyptpc"));
  if (!inputTree) {
    std::cerr << "Tree g4hyptpc was not found." << std::endl;
    return;
  }
  if (!inputTree->GetBranch("BEAM") || !inputTree->GetBranch("TGT")) {
    std::cerr << "Required BEAM and/or TGT branch is missing." << std::endl;
    return;
  }

  std::vector<TParticle>* beam = nullptr;
  std::vector<TParticle>* tgt = nullptr;
  inputTree->SetBranchAddress("BEAM", &beam);
  inputTree->SetBranchAddress("TGT", &tgt);

  TFile outputFile(output, "RECREATE");
  Long64_t nAccepted = 0;
  double p_initial = 0., p_target_in = 0.;
  double t_initial = 0., t_target_in = 0.;
  double delta_p = 0., delta_t = 0.;

  TH1D hDeltaT("h_delta_t",
               "K^{-} kinetic-energy loss before TargetLV;#DeltaT [MeV];Events",
               400, -1., 20.);
  TH1D hDeltaP("h_delta_p",
               "K^{-} momentum loss before TargetLV;#Delta p [MeV/#it{c}];Events",
               400, -1., 20.);

  Long64_t nNoBeam = 0, nNoTargetKaon = 0;
  for (Long64_t entry = 0, n = inputTree->GetEntries(); entry < n; ++entry) {
    inputTree->GetEntry(entry);
    if (!beam || beam->empty()) {
      ++nNoBeam;
      continue;
    }
    const TParticle* beamKaon = nullptr;
    for (const auto& particle : *beam) {
      if (particle.GetPdgCode() == -321) {
        beamKaon = &particle;
        break;
      }
    }
    const TParticle* targetKaon = tgt ? FindPrimaryKaonAtTarget(*tgt) : nullptr;
    if (!beamKaon || !targetKaon) {
      ++nNoTargetKaon;
      continue;
    }

    p_initial = beamKaon->P();
    p_target_in = targetKaon->P();
    t_initial = KineticEnergyMeV(*beamKaon);
    t_target_in = KineticEnergyMeV(*targetKaon);
    delta_p = p_initial - p_target_in;
    delta_t = t_initial - t_target_in;
    ++nAccepted;
    hDeltaT.Fill(delta_t);
    hDeltaP.Fill(delta_p);
  }

  // Fit the positive physical-loss region; the Landau location parameter is the MPV.
  if (hDeltaT.GetEntries() > 0) hDeltaT.Fit("landau", "Q", "", 3., 10.);
  if (hDeltaP.GetEntries() > 0) hDeltaP.Fit("landau", "Q", "", 3., 10.);
  hDeltaT.Write();
  hDeltaP.Write();
  TCanvas canvas("canvas", "Beam upstream energy loss", 800, 600);
  const TString pdfName(pdf);
  canvas.Print(pdfName + "[");
  auto drawWithLandau = [&canvas](TH1D& histogram) {
    canvas.Clear();
    histogram.Draw("hist");
    auto* fit = histogram.GetFunction("landau");
    if (!fit) return;
    fit->SetLineColor(kRed + 1);
    fit->SetLineWidth(2);
    fit->Draw("same");
    auto* label = new TPaveText(0.38, 0.68, 0.78, 0.88, "NDC");
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->SetTextAlign(12);
    label->AddText(Form("MPV = %.4f #pm %.4f", fit->GetParameter(1),
                        fit->GetParError(1)));
    label->AddText(Form("#sigma = %.4f #pm %.4f", fit->GetParameter(2),
                        fit->GetParError(2)));
    label->AddText(Form("#chi^{2}/NDF = %.1f / %d", fit->GetChisquare(),
                        fit->GetNDF()));
    label->Draw();
  };
  drawWithLandau(hDeltaT);
  canvas.Print(pdfName);
  drawWithLandau(hDeltaP);
  canvas.Print(pdfName);
  canvas.Print(pdfName + "]");
  outputFile.Close();
  std::cout << "Wrote histograms to " << pdf << std::endl;
  std::cout << "Wrote " << nAccepted << " events to " << output
            << "\nSkipped: no BEAM=" << nNoBeam
            << ", no primary K- TGT entry=" << nNoTargetKaon << std::endl;
}
