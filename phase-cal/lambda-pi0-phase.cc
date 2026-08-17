#include <array>
#include <iostream>
#include <string>

#include <TCanvas.h>
#include <TGenPhaseSpace.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>

void lambda_pi0_phase(Double_t beamMomentumMeV = 735.0)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  const Double_t kaonMass = 0.493;
  const Double_t protonMass = 0.938;
  const Double_t chargedPionMass = 0.1396;
  const Double_t neutralPionMass = 0.135;
  const Double_t lambdaMass = 1.115;
  const Double_t beamMomentum = beamMomentumMeV / 1000.0;

  TLorentzVector target(0.0, 0.0, 0.0, protonMass);
  TLorentzVector beam(0.0, 0.0, beamMomentum,
                      TMath::Sqrt(beamMomentum * beamMomentum +
                                  kaonMass * kaonMass));
  TLorentzVector centerOfMass = beam + target;
  const std::array<Double_t, 2> productionMasses = {
      lambdaMass, neutralPionMass};
  const std::array<Double_t, 2> lambdaDecayMasses = {
      protonMass, chargedPionMass};

  TGenPhaseSpace production;
  if (!production.SetDecay(centerOfMass, productionMasses.size(),
                           productionMasses.data())) {
    std::cerr << "K^{-}p -> #Lambda#pi^{0} is below threshold at p_{K-} = "
              << beamMomentumMeV << " MeV/c." << std::endl;
    return;
  }

  const std::array<const char *, 4> particleNames = {
      "lambda", "pi0", "proton_from_lambda", "pion_from_lambda"};
  const std::array<const char *, 4> particleTitles = {
      "#Lambda", "#pi^{0}", "p from #Lambda", "#pi^{-} from #Lambda"};
  std::array<TH2F *, 4> phaseSpaceHistograms = {
      new TH2F("h_lambda_phase_space", "#Lambda;cos#theta;momentum [GeV/#it{c}]",
               200, -1.0, 1.0, 200, 0.0, 1.2),
      new TH2F("h_pi0_phase_space", "#pi^{0};cos#theta;momentum [GeV/#it{c}]",
               200, -1.0, 1.0, 200, 0.0, 1.2),
      new TH2F("h_proton_from_lambda_phase_space",
               "p from #Lambda;cos#theta;momentum [GeV/#it{c}]",
               200, -1.0, 1.0, 200, 0.0, 1.2),
      new TH2F("h_pion_from_lambda_phase_space",
               "#pi^{-} from #Lambda;cos#theta;momentum [GeV/#it{c}]",
               200, -1.0, 1.0, 200, 0.0, 1.2)};
  const std::array<const char *, 4> componentNames = {"p", "px", "py", "pz"};
  const std::array<const char *, 4> componentTitles = {
      "|#vec{p}|", "p_{x}", "p_{y}", "p_{z}"};
  std::array<std::array<TH1F *, 4>, 4> momentumHistograms{};
  for (std::size_t particle = 0; particle < particleNames.size(); ++particle)
    for (std::size_t component = 0; component < componentNames.size(); ++component) {
      const Double_t minimum = component == 0 ? 0.0 : -1.2;
      momentumHistograms[particle][component] = new TH1F(
          Form("h_%s_%s", particleNames[particle], componentNames[component]),
          Form("%s: %s in lab;%s [GeV/#it{c}];Events",
               particleTitles[particle], componentTitles[component],
               componentTitles[component]), 240, minimum, 1.2);
    }

  TGenPhaseSpace lambdaDecay;
  constexpr Int_t generatedEvents = 100000;
  for (Int_t event = 0; event < generatedEvents; ++event) {
    production.Generate();
    TLorentzVector *lambda = production.GetDecay(0);
    TLorentzVector *pi0 = production.GetDecay(1);
    if (!lambdaDecay.SetDecay(*lambda, lambdaDecayMasses.size(),
                              lambdaDecayMasses.data())) {
      std::cerr << "Cannot decay #Lambda -> p#pi^{-}." << std::endl;
      return;
    }
    lambdaDecay.Generate();
    TLorentzVector *proton = lambdaDecay.GetDecay(0);
    TLorentzVector *pion = lambdaDecay.GetDecay(1);
    const std::array<TLorentzVector *, 4> particles = {
        lambda, pi0, proton, pion};
    for (std::size_t particle = 0; particle < particles.size(); ++particle) {
      phaseSpaceHistograms[particle]->Fill(particles[particle]->CosTheta(),
                                           particles[particle]->P());
      momentumHistograms[particle][0]->Fill(particles[particle]->P());
      momentumHistograms[particle][1]->Fill(particles[particle]->Px());
      momentumHistograms[particle][2]->Fill(particles[particle]->Py());
      momentumHistograms[particle][3]->Fill(particles[particle]->Pz());
    }
  }

  const std::string outputPdf =
      Form("lambda-pi0-phase-%.0fMeV.pdf", beamMomentumMeV);
  auto *canvas = new TCanvas("c_lambda_pi0", "lambda-pi0 phase space", 900, 700);
  canvas->Print((outputPdf + "(").c_str());
  for (TH2F *histogram : phaseSpaceHistograms) {
    histogram->Draw("colz");
    canvas->Print(outputPdf.c_str());
    canvas->Clear();
  }
  for (std::size_t particle = 0; particle < particleNames.size(); ++particle) {
    canvas->Divide(2, 2);
    for (std::size_t component = 0; component < componentNames.size(); ++component) {
      canvas->cd(component + 1);
      momentumHistograms[particle][component]->SetLineWidth(2);
      momentumHistograms[particle][component]->Draw("hist");
    }
    canvas->Print(particle + 1 == particleNames.size()
                      ? (outputPdf + ")").c_str() : outputPdf.c_str());
    canvas->Clear();
  }
}
