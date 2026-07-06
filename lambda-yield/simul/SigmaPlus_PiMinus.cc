#include <iostream>
#include <cmath>

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"

#include "../../basic-property.hh"

void SigmaPlus_PiMinus(int nEvent = 100000)
{
  // -----------------------------
  // Constants
  // Masses are assumed to be defined in basic-property.hh
  // mK, mp, mSigmap, mpi
  // Unit: MeV, cm
  // -----------------------------
  const double pK_lab    = 735.0;  // MeV/c
  const double cTauSigmap = 2.40;  // cm, Sigma+ c*tau

  TRandom3 rand(0);

  TH1D *hLmean = new TH1D("hLmean",
                          "Mean #Sigma^{+} decay length;#beta#gamma c#tau [cm];Counts",
                          200, 0, 5);

  TH1D *hVertexDistance = new TH1D("hVertexDistance",
                                   "#Sigma^{+} decay vertex distance;Distance from primary vertex [cm];Counts",
                                   200, 0, 15);

  TH2D *hPvsL = new TH2D("hPvsL",
                         "#Sigma^{+} momentum vs decay distance;Decay distance [cm];p_{#Sigma^{+}}^{lab} [MeV/c]",
                         200, 0, 15, 200, 0, 1200);

  // -----------------------------
  // Initial state in lab
  // K- beam + target proton at rest
  // -----------------------------
  TLorentzVector beamK(0.0, 0.0, pK_lab,
                       std::sqrt(pK_lab*pK_lab + mK*mK));

  TLorentzVector targetP(0.0, 0.0, 0.0, mp);

  TLorentzVector totalLab = beamK + targetP;

  const double W = totalLab.M();
  const double s = W * W;

  // -----------------------------
  // 2-body momentum in CM
  // K- p -> Sigma+ pi-
  // -----------------------------
  const double pCM =
    std::sqrt((s - std::pow(mSigmap + mpi, 2)) *
              (s - std::pow(mSigmap - mpi, 2))) / (2.0 * W);

  const double ESigmapCM = std::sqrt(pCM*pCM + mSigmap*mSigmap);
  const double EpiCM     = std::sqrt(pCM*pCM + mpi*mpi);

  TVector3 betaCM = totalLab.BoostVector();

  std::cout << "K- beam momentum = " << pK_lab << " MeV/c" << std::endl;
  std::cout << "W                = " << W      << " MeV"   << std::endl;
  std::cout << "pCM              = " << pCM    << " MeV/c" << std::endl;

  for (int i = 0; i < nEvent; i++) {

    // -----------------------------
    // Isotropic generation in CM
    // -----------------------------
    double cosTheta = rand.Uniform(-1.0, 1.0);
    double sinTheta = std::sqrt(1.0 - cosTheta*cosTheta);
    double phi      = rand.Uniform(0.0, 2.0 * M_PI);

    double px = pCM * sinTheta * std::cos(phi);
    double py = pCM * sinTheta * std::sin(phi);
    double pz = pCM * cosTheta;

    TLorentzVector sigmapCM(px, py, pz, ESigmapCM);
    TLorentzVector pimCM(-px, -py, -pz, EpiCM);

    // -----------------------------
    // Boost CM -> lab
    // -----------------------------
    TLorentzVector sigmapLab = sigmapCM;
    TLorentzVector pimLab    = pimCM;

    sigmapLab.Boost(betaCM);
    pimLab.Boost(betaCM);

    // -----------------------------
    // Primary vertex
    // pi- is produced here
    // -----------------------------
    TVector3 vtxPrimary(0.0, 0.0, 0.0);

    // -----------------------------
    // Sigma+ decay vertex
    //
    // Mean decay length:
    // Lmean = beta * gamma * cTau
    //       = p / m * cTau
    // -----------------------------
    double pSigmapLab = sigmapLab.P();
    double Lmean = (pSigmapLab / mSigmap) * cTauSigmap;

    // Actual decay length follows exponential distribution
    double Ldecay = rand.Exp(Lmean);

    TVector3 sigmapDir = sigmapLab.Vect().Unit();
    TVector3 vtxSigmapDecay = vtxPrimary + Ldecay * sigmapDir;

    double vertexDistance = (vtxSigmapDecay - vtxPrimary).Mag();

    hLmean->Fill(Lmean);
    hVertexDistance->Fill(vertexDistance);
    hPvsL->Fill(vertexDistance, pSigmapLab);
  }

  std::cout << "Generated events = " << nEvent << std::endl;
  std::cout << "Mean vertex distance = "
            << hVertexDistance->GetMean() << " cm" << std::endl;

  TCanvas *c1 = new TCanvas("c1", "Sigma+ vertex distance", 800, 600);
  hVertexDistance->Draw();

  TCanvas *c2 = new TCanvas("c2", "Sigma+ mean decay length", 800, 600);
  hLmean->Draw();

  TCanvas *c3 = new TCanvas("c3", "Sigma+ momentum vs vertex distance", 800, 600);
  hPvsL->Draw("colz");
}
