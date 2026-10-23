// Run: root -l -b -q lambda_eta_momentum_offset_study.cc
// Six 2D plots: Lambda invariant mass / eta missing mass vs p, pi, beam offset.
// Resolution prescription/values copied from lambda_missing_mass_resolution_study.cc:
//   kBeamRelativeMomentumSigma   = 0.001;
//   kProtonRelativeMomentumSigma = 0.040;
//   kPionRelativeMomentumSigma   = 0.040;
//   SmearMagnitude: p_true * Gaus(1, relative_sigma), with direction retained.
// Here p_used = p_true * Gaus(1, sigma_p/p) + offset in magnitude.
// Offset = p_used - p_true before random resolution; negative = underestimated p.
// Lambda mass plot: unconstrained M(p pi-).
// Eta missing mass: PDG mass constraint by default; last argument false uses Pp+Ppi.
// Unconstrained example:
// root -l -b -q 'lambda_eta_momentum_offset_study.cc(725.,100000,-20.,20.,81,.04,.04,.001,"lambda_eta_momentum_offset_725MeV_fourvector.pdf",false,false)'
// Only magnitude resolution is used, not its extra Lambda angular smearing.
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <limits>
#include <stdexcept>
#include <vector>
#include <TCanvas.h>
#include <TAxis.h>
#include <TGenPhaseSpace.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TStyle.h>

namespace eta_offset {
constexpr double kKaonMass = .493677;
constexpr double kProtonMass = .938272;
constexpr double kPionMass = .139570;
constexpr double kLambdaMass = 1.115683;
constexpr double kEtaMass = .547862;
// Keep these aligned with lambda_missing_mass_resolution_study.cc.
constexpr double kBeamRelativeMomentumSigma = .001;
constexpr double kProtonRelativeMomentumSigma = .040;
constexpr double kPionRelativeMomentumSigma = .040;
TLorentzVector Four(const TVector3& momentum, double mass) {
  TLorentzVector result;
  result.SetVectM(momentum, mass);
  return result;
}
}

void lambda_eta_momentum_offset_study(
    double beamMomentumMeV = 725., int nEvents = 100000,
    double offsetMinMeV = -20., double offsetMaxMeV = 20., int offsetPoints = 81,
    double protonRelativeSigma = eta_offset::kProtonRelativeMomentumSigma,
    double pionRelativeSigma = eta_offset::kPionRelativeMomentumSigma,
    double beamRelativeSigma = eta_offset::kBeamRelativeMomentumSigma,
    const char* outputPdf = "lambda_eta_momentum_offset_725MeV.pdf",
    bool allowedRegionOnly = false, bool constrainLambdaMass = true,
    double angleOffsetMinMrad = -20., double angleOffsetMaxMrad = 20.,
    int angleOffsetPoints = 81) {
  using namespace eta_offset;
  if (!std::isfinite(angleOffsetMinMrad) || !std::isfinite(angleOffsetMaxMrad) ||
      angleOffsetMaxMrad<=angleOffsetMinMrad || angleOffsetPoints<2)
    throw std::runtime_error("Invalid opening-angle offset scan");
  const char* lambdaLabel = constrainLambdaMass ? "#Lambda_{PDG}" : "(p#pi^{-})_{reco}";
  auto missingLambda = [&](const TLorentzVector& daughters) {
    return constrainLambdaMass ? Four(daughters.Vect(), kLambdaMass) : daughters;
  };
  auto formatLambdaAxis = [](TAxis* axis) {
    axis->SetNdivisions(504);
  };
  if (!std::isfinite(beamMomentumMeV) || beamMomentumMeV < 0. || nEvents <= 0 ||
      !std::isfinite(offsetMinMeV) || !std::isfinite(offsetMaxMeV) ||
      offsetMaxMeV <= offsetMinMeV || offsetPoints < 2 ||
      !std::isfinite(protonRelativeSigma) || protonRelativeSigma < 0. ||
      !std::isfinite(pionRelativeSigma) || pionRelativeSigma < 0. ||
      !std::isfinite(beamRelativeSigma) || beamRelativeSigma < 0.)
    throw std::runtime_error("Invalid beam, event count, offset scan, or resolution");

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  auto initial = Four(TVector3(0., 0., beamMomentumMeV / 1000.), kKaonMass)
               + Four(TVector3(), kProtonMass);
  double productionMasses[2] = {kLambdaMass, kEtaMass};
  double decayMasses[2] = {kProtonMass, kPionMass};
  TGenPhaseSpace production, decay;
  if (!production.SetDecay(initial, 2, productionMasses))
    throw std::runtime_error("True beam is at/below Lambda eta threshold");

  TRandom3 random(20261006);
  struct RandomGuard {
    TRandom* previous;
    ~RandomGuard() { gRandom = previous; }
  } guard{gRandom};
  gRandom = &random; // TGenPhaseSpace uses gRandom too.

  const double offsetStep = (offsetMaxMeV - offsetMinMeV) / (offsetPoints - 1);
  const double yMin = offsetMinMeV - offsetStep / 2.;
  const double yMax = offsetMaxMeV + offsetStep / 2.;
  const std::array<const char*, 3> modeNames{{"Proton", "Pion", "Beam"}};
  const std::array<const char*, 3> modeIds{{"proton", "pion", "beam"}};
  const std::array<const char*, 3> offsetLabels{{"#delta p_{p}", "#delta p_{#pi}", "#delta p_{K^{-}}"}};
  const std::array<double, 2> massMin{{1.09, .48}}, massMax{{1.14, .60}};
  std::array<std::array<std::unique_ptr<TH2D>, 2>, 3> histograms;
  std::array<Long64_t, 3> invalid{{0, 0, 0}};
  std::array<std::array<Long64_t, 2>, 3> outsideMassRange{};
  struct MissingStatus {
    Long64_t validTracks = 0, negativeM2 = 0, belowCut = 0;
    double minM2 = std::numeric_limits<double>::infinity();
  };
  std::array<std::vector<MissingStatus>, 3> missingStatus;
  for (auto& status : missingStatus) status.resize(offsetPoints);
  for (int mode = 0; mode < 3; ++mode) {
    for (int kind = 0; kind < 2; ++kind) {
      const char* particle = kind == 0 ? "lambda" : "eta";
      const char* axis = kind == 0 ? "M(p#pi^{-}) [GeV/c^{2}]" : "M_{miss}(K^{-}p-#Lambda) [GeV/c^{2}]";
      histograms[mode][kind].reset(new TH2D(
        Form("%s_mass_%s_offset", particle, modeIds[mode]),
        Form("%s offset: %s (p/#pi/K resolution %.1f/%.1f/%.1f%%);%s;%s [MeV/c]",
             modeNames[mode], kind == 0 ? "#Lambda" : "#eta missing mass",
             100.*protonRelativeSigma, 100.*pionRelativeSigma, 100.*beamRelativeSigma,
             axis, offsetLabels[mode]),
        400, massMin[kind], massMax[kind], offsetPoints, yMin, yMax));
      histograms[mode][kind]->SetDirectory(nullptr);
      if (kind == 0) formatLambdaAxis(histograms[mode][kind]->GetXaxis());
    }
  }
  const auto target = Four(TVector3(), kProtonMass);
  struct SmearedEvent {
    TVector3 protonDirection, pionDirection;
    double proton, pion, beam;
  };
  std::vector<SmearedEvent> smearedEvents;
  smearedEvents.reserve(nEvents);
  // Isotropic two-body production/decay toy, no acceptance or material model.
  // Each event uses the same independent proton/pion Gaussian draws at all
  // offsets and all modes, so scan differences reflect the offset.
  // Other particles have zero offset but retain their momentum resolution.
  for (int event = 0; event < nEvents; ++event) {
    production.Generate();
    auto lambda = *production.GetDecay(0);
    if (!decay.SetDecay(lambda, 2, decayMasses))
      throw std::runtime_error("Lambda decay failed");
    decay.Generate();
    const auto protonTrue = decay.GetDecay(0)->Vect();
    const auto pionTrue = decay.GetDecay(1)->Vect();
    const double truthMass = (Four(protonTrue, kProtonMass)
                            + Four(pionTrue, kPionMass)).M();
    if (std::abs(truthMass - kLambdaMass) > 1.e-9)
      throw std::runtime_error("Truth Lambda mass closure failed");
    const auto truthMissing = initial - Four(protonTrue + pionTrue, kLambdaMass);
    if (std::abs(truthMissing.M() - kEtaMass) > 1.e-9)
      throw std::runtime_error("Truth eta missing mass closure failed");
    const double beamSmeared = beamMomentumMeV / 1000. * random.Gaus(1., beamRelativeSigma);
    const auto protonDirection = protonTrue.Unit();
    const auto pionDirection = pionTrue.Unit();
    const double protonSmeared = protonTrue.Mag() * random.Gaus(1., protonRelativeSigma);
    const double pionSmeared = pionTrue.Mag() * random.Gaus(1., pionRelativeSigma);
    smearedEvents.push_back({protonDirection, pionDirection,
                             protonSmeared, pionSmeared, beamSmeared});
    for (int j = 0; j < offsetPoints; ++j) {
      const double offsetMeV = offsetMinMeV + j * offsetStep;
      for (int mode = 0; mode < 3; ++mode) {
        const double p = protonSmeared + (mode == 0 ? offsetMeV / 1000. : 0.);
        const double pi = pionSmeared + (mode == 1 ? offsetMeV / 1000. : 0.);
        const double beam = beamSmeared + (mode == 2 ? offsetMeV / 1000. : 0.);
        // Skip negative magnitudes instead of clamping or reversing direction.
        if (p < 0. || pi < 0. || beam < 0.) { ++invalid[mode]; continue; }
        const auto lambdaReco = Four(protonDirection * p, kProtonMass)
                              + Four(pionDirection * pi, kPionMass);
        const auto lambdaForMissing = missingLambda(lambdaReco);
        const auto missing = Four(TVector3(0., 0., beam), kKaonMass)
                           + target - lambdaForMissing;
        auto& status = missingStatus[mode][j];
        ++status.validTracks;
        const double missingM2 = missing.M2();
        status.minM2 = std::min(status.minM2, missingM2);
        if (missingM2 < 0.) ++status.negativeM2;
        else if (missingM2 < .5*.5) ++status.belowCut;
        // Negative Mmiss^2 has no real mass: count it, but do not fill the
        // eta mass plot using ROOT's signed square-root convention.
        const std::array<double, 2> masses{{lambdaReco.M(),
          missingM2 >= 0. ? std::sqrt(missingM2) : 0.}};
        for (int kind = 0; kind < 2; ++kind) {
          if (kind == 1 && missingM2 < 0.) continue;
          histograms[mode][kind]->Fill(masses[kind], offsetMeV);
          if (masses[kind] < massMin[kind] || masses[kind] >= massMax[kind])
            ++outsideMassRange[mode][kind];
        }
      }
    }
  }
  // Equal scales within each particle's three panels for comparison.
  std::array<double, 2> maximum{{0., 0.}};
  for (int mode = 0; mode < 3; ++mode)
    for (int kind = 0; kind < 2; ++kind)
      maximum[kind] = std::max(maximum[kind], histograms[mode][kind]->GetMaximum());
  std::array<double, 3> maxBelowCutPercent{};
  for (int mode = 0; mode < 3; ++mode) {
    Long64_t below = 0, negative = 0;
    for (const auto& status : missingStatus[mode]) {
      below += status.belowCut; negative += status.negativeM2;
      if (status.validTracks)
        maxBelowCutPercent[mode] = std::max(maxBelowCutPercent[mode],
          100. * status.belowCut / status.validTracks);
    }
    Info("eta_mass_cut", "%s, beam %.0f MeV/c, offset %.1f to %.1f: Mmiss < 0.5 events %lld; negative M2 %lld; maximum below-cut fraction %.4f%%",
         modeNames[mode], beamMomentumMeV, offsetMinMeV, offsetMaxMeV,
         below, negative, maxBelowCutPercent[mode]);
  }
  // Beam offset cannot affect M(p pi-); all its offset rows must be identical.
  for (int x = 0; x <= histograms[2][0]->GetNbinsX() + 1; ++x)
    for (int y = 2; y <= offsetPoints; ++y)
      if (invalid[2] == 0 &&
          histograms[2][0]->GetBinContent(x, y) != histograms[2][0]->GetBinContent(x, 1))
        throw std::runtime_error("Beam-only offset changed Lambda invariant mass");

  TCanvas canvas("c_lambda_eta_mass_offsets", "Mass vs momentum offset", 1600, 750);
  canvas.Print(TString(outputPdf) + "[");
  for (int mode = 0; mode < 3; ++mode) {
    canvas.Divide(2, 1);
    TLine cutLine(.5, yMin, .5, yMax);
    cutLine.SetLineColor(kRed); cutLine.SetLineStyle(2); cutLine.SetLineWidth(2);
    std::array<TLine, 2> massLines{{
      TLine(kLambdaMass, yMin, kLambdaMass, yMax),
      TLine(kEtaMass, yMin, kEtaMass, yMax)
    }};
    for (int kind = 0; kind < 2; ++kind) {
      canvas.cd(kind + 1);
      gPad->SetLeftMargin(.14);
      gPad->SetRightMargin(.17);
      gPad->SetBottomMargin(.14);
      auto& histogram = *histograms[mode][kind];
      histogram.SetMinimum(0.);
      histogram.SetMaximum(maximum[kind]);
      histogram.Draw("colz");
      massLines[kind].SetLineColor(kRed);
      massLines[kind].SetLineWidth(2);
      massLines[kind].Draw("same");
      if (kind == 1) cutLine.Draw("same");
      Info("lambda_eta_momentum_offset_study",
           "%s offset, %s: entries %.0f; invalid %lld; outside mass axis %lld",
           modeNames[mode], kind == 0 ? "Lambda" : "eta missing mass",
           histogram.GetEntries(), invalid[mode], outsideMassRange[mode][kind]);
    }
    canvas.Print(outputPdf);
    canvas.Clear();
  }
  // Wider offset scan focused on migration below the eta region Mmiss >= 0.5.
  // Same truth and resolution draws as above; no new events or extra files.
  constexpr double diagnosticMin = -200., diagnosticMax = 200., diagnosticStep = 2.;
  constexpr int diagnosticPoints = 201;
  const double diagnosticYMin = diagnosticMin - diagnosticStep / 2.;
  const double diagnosticYMax = diagnosticMax + diagnosticStep / 2.;
  TCanvas wideCanvas("c_eta_cut_wide", "Eta mass and below-cut fraction", 1600, 750);
  for (int mode = 0; mode < 3; ++mode) {
    TH2D mass(Form("eta_mass_wide_%s", modeIds[mode]),
      Form("%s offset: #eta (beam %.0f MeV/c);M_{miss} [GeV/c^{2}];%s [MeV/c]",
           modeNames[mode], beamMomentumMeV, offsetLabels[mode]),
      500, .35, .70, diagnosticPoints, diagnosticYMin, diagnosticYMax);
    mass.SetDirectory(nullptr);
    TGraph belowFraction;
    belowFraction.SetName(Form("below_cut_fraction_%s", modeIds[mode]));
    int firstNegative = -1, firstPositive = -1;
    Long64_t skipped = 0, negativeTotal = 0;
    for (int j = 0; j < diagnosticPoints; ++j) {
      const double offset = diagnosticMin + j * diagnosticStep;
      Long64_t below = 0, validTracks = 0;
      for (const auto& e : smearedEvents) {
        const double p = e.proton + (mode == 0 ? offset / 1000. : 0.);
        const double pi = e.pion + (mode == 1 ? offset / 1000. : 0.);
        const double beam = e.beam + (mode == 2 ? offset / 1000. : 0.);
        if (p < 0. || pi < 0. || beam < 0.) { ++skipped; continue; }
        ++validTracks;
        const auto daughters = Four(e.protonDirection*p, kProtonMass)
                             + Four(e.pionDirection*pi, kPionMass);
        const double m2 = (Four(TVector3(0., 0., beam), kKaonMass) + target
                        - missingLambda(daughters)).M2();
        if (m2 < 0.) { ++negativeTotal; continue; }
        const double missingMass = std::sqrt(m2);
        mass.Fill(missingMass, offset);
        if (missingMass < .5) ++below;
      }
      // Denominator = events with valid reconstructed track magnitudes.
      // Negative M2 is counted separately, not mislabeled as a real M < 0.5.
      if (validTracks) belowFraction.SetPoint(belowFraction.GetN(),
        100. * below / validTracks, offset);
      if (below && offset <= 0.) firstNegative = j;
      if (below && offset > 0. && firstPositive < 0) firstPositive = j;
      if (j == 0 || j == diagnosticPoints-1)
        Info("eta_mass_cut_wide", "%s offset %.0f: below 0.5 = %lld/%lld (%.4f%%)",
             modeNames[mode], offset, below, validTracks,
             validTracks ? 100.*below/validTracks : 0.);
    }
    wideCanvas.Divide(2, 1);
    wideCanvas.cd(1); gPad->SetLeftMargin(.14); gPad->SetRightMargin(.17);
    gPad->SetBottomMargin(.14);
    mass.Draw("colz");
    TLine cut(.5, diagnosticYMin, .5, diagnosticYMax);
    cut.SetLineColor(kRed); cut.SetLineStyle(2); cut.SetLineWidth(2); cut.Draw("same");
    TLine eta(kEtaMass, diagnosticYMin, kEtaMass, diagnosticYMax);
    eta.SetLineColor(kRed); eta.SetLineWidth(2); eta.Draw("same");
    wideCanvas.cd(2); gPad->SetLeftMargin(.14); gPad->SetRightMargin(.10);
    gPad->SetBottomMargin(.14);
    auto* frame = gPad->DrawFrame(-1., diagnosticYMin, 101., diagnosticYMax);
    frame->SetTitle(Form("%s offset: fraction below 0.5;Events with 0 #leq M_{miss} < 0.5 [%%];%s [MeV/c]",
                         modeNames[mode], offsetLabels[mode]));
    belowFraction.SetLineColor(kBlue+1); belowFraction.SetLineWidth(2);
    belowFraction.Draw("L same");
    for (int j : {firstNegative, firstPositive}) {
      if (j < 0) continue;
      const double offset = diagnosticMin + j * diagnosticStep;
      Info("eta_mass_cut_wide", "%s: closest sampled offset with Mmiss < 0.5: %.0f MeV/c",
           modeNames[mode], offset);
    }
    Info("eta_mass_cut_wide", "%s: invalid track magnitudes %lld; negative M2 %lld",
         modeNames[mode], skipped, negativeTotal);
    wideCanvas.Print(outputPdf);
    wideCanvas.Clear();
  }
  // x uses the selected Lambda energy prescription; y always uses M(p pi).
  // The boundary is a reference when the PDG mass constraint is used,
  // and a consistent joint-mass boundary when the daughter sum is used.
  const double sqrtS = initial.M();
  constexpr double jointXMin = .48, jointXMax = .60;
  constexpr double jointYMin = 1.08, jointYMax = 1.15;
  TH2D joint("lambda_eta_allowed_mass_region",
    Form("K^{-}p reconstructed mass pairs, beam %.0f MeV/c;M_{miss}(K^{-}p-%s) [GeV/c^{2}];M(p#pi^{-}) [GeV/c^{2}]", beamMomentumMeV, lambdaLabel),
    480, jointXMin, jointXMax, 400, jointYMin, jointYMax);
  joint.SetDirectory(nullptr);
  formatLambdaAxis(joint.GetYaxis());
  Long64_t jointRejected = 0, jointOutside = 0, jointOutsideRegion = 0;
  if (allowedRegionOnly) {
    // A geometric mask, not an event density or detector acceptance.
    for (int x = 1; x <= joint.GetNbinsX(); ++x)
      for (int y = 1; y <= joint.GetNbinsY(); ++y) {
        const double missingMass = joint.GetXaxis()->GetBinCenter(x);
        const double lambdaMass = joint.GetYaxis()->GetBinCenter(y);
        if (missingMass >= 0. && lambdaMass >= kProtonMass+kPionMass &&
            missingMass+lambdaMass <= sqrtS)
          joint.SetBinContent(x, y, 1.);
      }
  } else {
    for (const auto& e : smearedEvents) {
      if (e.proton < 0. || e.pion < 0.) { ++jointRejected; continue; }
      const auto lambda = Four(e.protonDirection*e.proton, kProtonMass)
                        + Four(e.pionDirection*e.pion, kPionMass);
      const auto missing = initial-missingLambda(lambda);
      if (missing.M2() < 0.) { ++jointRejected; continue; }
      const double x = std::sqrt(missing.M2()), y = lambda.M();
      if (x+y > sqrtS+1.e-10) ++jointOutsideRegion;
      joint.Fill(x, y);
      if (x < jointXMin || x >= jointXMax || y < jointYMin || y >= jointYMax)
        ++jointOutside;
    }
    Info("allowed_mass_events", "Generated %d; plotted %.0f; outside reference region %lld; undefined mass/track %lld; outside displayed axes %lld",
         nEvents, joint.GetEntries(), jointOutsideRegion, jointRejected, jointOutside);
  }
  TCanvas jointCanvas("c_allowed_mass_region", "Kinematically allowed mass region", 1000, 950);
  jointCanvas.SetLeftMargin(.14); jointCanvas.SetRightMargin(.15); jointCanvas.SetBottomMargin(.14);
  joint.SetMinimum(0.); joint.Draw(allowedRegionOnly ? "col" : "colz");
  const double boundaryLeft = std::max(jointXMin, sqrtS-jointYMax);
  const double boundaryRight = std::min(jointXMax, sqrtS-jointYMin);
  TLine boundary(boundaryLeft, sqrtS-boundaryLeft, boundaryRight, sqrtS-boundaryRight);
  boundary.SetLineColor(kBlack); boundary.SetLineWidth(3); boundary.Draw("same");
  TLine etaRef(kEtaMass, jointYMin, kEtaMass, jointYMax);
  TLine lambdaRef(jointXMin, kLambdaMass, jointXMax, kLambdaMass);
  etaRef.SetLineColor(kRed); etaRef.SetLineWidth(1); etaRef.Draw("same");
  lambdaRef.SetLineColor(kRed); lambdaRef.SetLineWidth(1); lambdaRef.Draw("same");
  jointCanvas.Print(outputPdf);

  // Offset changes the reconstruction, not the truth production at 725 MeV/c.
  // Keep the nominal-beam allowed region fixed for all comparison panels.
  const std::array<double, 3> jointOffsets{{-20., 0., 20.}}; // MeV/c
  auto scenarioOffsets = [&](int mode, int panel) {
    std::array<double, 3> d{{0., 0., 0.}}; // p, pi, beam in GeV/c
    if (panel == 3) d = {{-.020, -.020, .020}};
    else d[mode] = jointOffsets[panel]/1000.;
    return d;
  };
  auto scenarioTitle = [&](int mode, int panel) -> TString {
    if (panel == 3) return "Combined: p -20, #pi -20, K +20 MeV/c";
    return Form("%s offset %+.0f MeV/c", modeNames[mode], jointOffsets[panel]);
  };
  TCanvas offsetJointCanvas("c_joint_mass_offsets", "Mass correlation with offsets", 1600, 1400);
  for (int mode = 0; mode < 3; ++mode) {
    std::array<std::unique_ptr<TH2D>, 4> distributions;
    std::array<Long64_t, 4> outsideNominalRegion{};
    double colorMaximum = 0.;
    for (int panel = 0; panel < 4; ++panel) {
      const auto offsets = scenarioOffsets(mode, panel);
      const double usedBeam = beamMomentumMeV/1000. + offsets[2];
      const auto panelInitial = Four(TVector3(0., 0., usedBeam), kKaonMass)+target;
      distributions[panel].reset(new TH2D(Form("joint_%s_offset_%d", modeIds[mode], panel),
        Form("%s, beam %.0f MeV/c;M_{miss}(K^{-}p-%s) [GeV/c^{2}];M(p#pi^{-}) [GeV/c^{2}]",
             scenarioTitle(mode,panel).Data(), usedBeam*1000., lambdaLabel),
        300, jointXMin, jointXMax, 250, jointYMin, jointYMax));
      distributions[panel]->SetDirectory(nullptr);
      formatLambdaAxis(distributions[panel]->GetYaxis());
      Long64_t rejected = 0, outside = 0;
      for (const auto& e : smearedEvents) {
        const double p = e.proton+offsets[0];
        const double pi = e.pion+offsets[1];
        if (p < 0. || pi < 0. || usedBeam < 0.) { ++rejected; continue; }
        const auto lambda = Four(e.protonDirection*p, kProtonMass)
                          + Four(e.pionDirection*pi, kPionMass);
        const auto missing = panelInitial-missingLambda(lambda);
        if (missing.M2() < 0.) { ++rejected; continue; }
        const double x = std::sqrt(missing.M2()), y = lambda.M();
        if (x+y > sqrtS+1.e-10) ++outsideNominalRegion[panel];
        distributions[panel]->Fill(x, y);
        if (x < jointXMin || x >= jointXMax || y < jointYMin || y >= jointYMax)
          ++outside;
      }
      colorMaximum = std::max(colorMaximum, distributions[panel]->GetMaximum());
      Info("joint_mass_offset", "%s: plotted %.0f; outside nominal region %lld; invalid %lld; outside axes %lld",
           scenarioTitle(mode,panel).Data(), distributions[panel]->GetEntries(), outsideNominalRegion[panel], rejected, outside);
    }
    offsetJointCanvas.Divide(2, 2);
    std::array<TLine, 4> boundaries, etaLines, lambdaLines;
    for (int panel = 0; panel < 4; ++panel) {
      offsetJointCanvas.cd(panel+1);
      gPad->SetLeftMargin(.17); gPad->SetRightMargin(.18); gPad->SetBottomMargin(.15);
      distributions[panel]->SetMinimum(0.); distributions[panel]->SetMaximum(colorMaximum);
      distributions[panel]->Draw("colz");
      const double left = std::max(jointXMin, sqrtS-jointYMax);
      const double right = std::min(jointXMax, sqrtS-jointYMin);
      boundaries[panel] = TLine(left, sqrtS-left, right, sqrtS-right);
      boundaries[panel].SetLineWidth(3); boundaries[panel].Draw("same");
      etaLines[panel] = TLine(kEtaMass, jointYMin, kEtaMass, jointYMax);
      lambdaLines[panel] = TLine(jointXMin, kLambdaMass, jointXMax, kLambdaMass);
      for (auto* line : {&etaLines[panel], &lambdaLines[panel]}) {
        line->SetLineColor(kRed); line->SetLineWidth(1); line->Draw("same");
      }
    }
    offsetJointCanvas.Print(outputPdf);
    offsetJointCanvas.Clear();
  }

  // Fixed PDG Lambda mass, varying true beam momentum; geometric mask only.
  constexpr double beamPlotMin = .71, beamPlotMax = .77;
  constexpr double missingPlotMin = .50, missingPlotMax = .58;
  TH2D beamRegion("beam_missing_mass_allowed_region",
    "K^{-}p #rightarrow #Lambda X: allowed region;M_{miss} [GeV/c^{2}];p_{K^{-}} [GeV/c]",
    500, missingPlotMin, missingPlotMax, 500, beamPlotMin, beamPlotMax);
  beamRegion.SetDirectory(nullptr);
  auto maxMissingMass = [&](double beamMomentum) {
    const double s = kKaonMass*kKaonMass + kProtonMass*kProtonMass
                   + 2.*kProtonMass*std::sqrt(beamMomentum*beamMomentum+kKaonMass*kKaonMass);
    return std::sqrt(s)-kLambdaMass;
  };
  for (int y = 1; y <= beamRegion.GetNbinsY(); ++y) {
    const double allowedMax = maxMissingMass(beamRegion.GetYaxis()->GetBinCenter(y));
    for (int x = 1; x <= beamRegion.GetNbinsX(); ++x)
      if (beamRegion.GetXaxis()->GetBinCenter(x) <= allowedMax)
        beamRegion.SetBinContent(x, y, 1.);
  }
  TGraph beamBoundary;
  for (int j = 0; j <= 600; ++j) {
    const double momentum = beamPlotMin+(beamPlotMax-beamPlotMin)*j/600.;
    beamBoundary.SetPoint(j, maxMissingMass(momentum), momentum);
  }
  const double thresholdEnergy = ((kLambdaMass+kEtaMass)*(kLambdaMass+kEtaMass)
    - kProtonMass*kProtonMass-kKaonMass*kKaonMass)/(2.*kProtonMass);
  const double thresholdBeam = std::sqrt(thresholdEnergy*thresholdEnergy-kKaonMass*kKaonMass);
  TCanvas beamCanvas("c_beam_missing_mass_allowed", "Beam momentum vs allowed missing mass", 1000, 1000);
  beamCanvas.SetLeftMargin(.14); beamCanvas.SetRightMargin(.10); beamCanvas.SetBottomMargin(.14);
  beamRegion.SetMinimum(0.); beamRegion.Draw("col");
  beamBoundary.SetLineColor(kBlack); beamBoundary.SetLineWidth(3); beamBoundary.Draw("L same");
  TLine etaBeamLine(kEtaMass, beamPlotMin, kEtaMass, beamPlotMax);
  etaBeamLine.SetLineColor(kRed); etaBeamLine.SetLineWidth(2); etaBeamLine.Draw("same");
  TLine thresholdLine(missingPlotMin, thresholdBeam, missingPlotMax, thresholdBeam);
  thresholdLine.SetLineColor(kRed); thresholdLine.SetLineStyle(2); thresholdLine.SetLineWidth(2); thresholdLine.Draw("same");
  beamCanvas.Print(outputPdf);

  // Uniform random TRUE beam over the requested interval. Offsets are
  // reconstruction biases only; y is the same true pK in all twelve panels.
  // Below-threshold proposals cannot produce Lambda eta. No cross-section
  // or detector acceptance weighting is imposed on this isotropic toy.
  std::array<std::array<std::unique_ptr<TH2D>, 4>, 3> randomMaps;
  struct RandomStatus {
    Long64_t undefined = 0, outsideBoundary = 0, outsideAxes = 0;
    double massShiftSum = 0.;
  };
  std::array<std::array<RandomStatus, 4>, 3> randomStatus{};
  for (int mode = 0; mode < 3; ++mode) {
    for (int panel = 0; panel < 4; ++panel) {
      randomMaps[mode][panel].reset(new TH2D(Form("random_pK_%s_offset_%d", modeIds[mode], panel),
        Form("Uniform p_{K,true}: %s;M_{miss} (%s) [GeV/c^{2}];p_{K,true} [GeV/c]",
             scenarioTitle(mode,panel).Data(), lambdaLabel),
        400, .48, .60, 300, beamPlotMin, beamPlotMax));
      randomMaps[mode][panel]->SetDirectory(nullptr);
    }
  }
  TGenPhaseSpace randomProduction, randomDecay;
  Long64_t belowThreshold = 0, produced = 0;
  for (int attempt = 0; attempt < nEvents; ++attempt) {
    const double trueBeam = random.Uniform(beamPlotMin, beamPlotMax);
    auto randomInitial = Four(TVector3(0., 0., trueBeam), kKaonMass)+target;
    if (!randomProduction.SetDecay(randomInitial, 2, productionMasses)) {
      ++belowThreshold; continue;
    }
    randomProduction.Generate();
    auto lambda = *randomProduction.GetDecay(0);
    if (!randomDecay.SetDecay(lambda, 2, decayMasses))
      throw std::runtime_error("Random-beam Lambda decay failed");
    randomDecay.Generate();
    ++produced;
    const auto protonTrue = randomDecay.GetDecay(0)->Vect();
    const auto pionTrue = randomDecay.GetDecay(1)->Vect();
    const double protonMeasured = protonTrue.Mag()*random.Gaus(1., protonRelativeSigma);
    const double pionMeasured = pionTrue.Mag()*random.Gaus(1., pionRelativeSigma);
    const double measuredBeam = trueBeam*random.Gaus(1., beamRelativeSigma);
    if (std::abs((randomInitial-lambda).M()-kEtaMass) > 1.e-9)
      throw std::runtime_error("Random-beam truth missing mass closure failed");
    for (int mode = 0; mode < 3; ++mode) {
      for (int panel = 0; panel < 4; ++panel) {
        auto& status = randomStatus[mode][panel];
        const auto offsets = scenarioOffsets(mode,panel);
        const double p = protonMeasured+offsets[0];
        const double pi = pionMeasured+offsets[1];
        const double usedBeam = measuredBeam+offsets[2];
        if (p < 0. || pi < 0. || usedBeam < 0.) { ++status.undefined; continue; }
        const auto daughters = Four(protonTrue.Unit()*p, kProtonMass)
                             + Four(pionTrue.Unit()*pi, kPionMass);
        const auto reconstructedLambda = missingLambda(daughters);
        const auto missing = Four(TVector3(0., 0., usedBeam), kKaonMass)+target-reconstructedLambda;
        if (missing.M2() < 0.) { ++status.undefined; continue; }
        const double missingMass = std::sqrt(missing.M2());
        randomMaps[mode][panel]->Fill(missingMass, trueBeam);
        status.massShiftSum += (missingMass-kEtaMass)*1000.;
        if (missingMass > maxMissingMass(trueBeam)+1.e-10) ++status.outsideBoundary;
        if (missingMass < .48 || missingMass >= .60) ++status.outsideAxes;
      }
    }
  }
  Info("random_beam_events", "Uniform true pK %.2f to %.2f: proposals %d; below threshold %lld; produced %lld",
       beamPlotMin, beamPlotMax, nEvents, belowThreshold, produced);
  // All zero-offset panels must contain the exact same paired sample.
  for (int mode = 1; mode < 3; ++mode)
    for (int x = 0; x <= randomMaps[0][1]->GetNbinsX()+1; ++x)
      for (int y = 0; y <= randomMaps[0][1]->GetNbinsY()+1; ++y)
        if (randomMaps[mode][1]->GetBinContent(x,y) != randomMaps[0][1]->GetBinContent(x,y))
          throw std::runtime_error("Random-pK zero-offset panels do not match");
  double randomColorMaximum = 0.;
  for (auto& row : randomMaps)
    for (auto& h : row) randomColorMaximum = std::max(randomColorMaximum, h->GetMaximum());
  TCanvas randomBeamCanvas("c_random_pK_offsets", "Uniform beam with track and beam offsets", 1600, 1400);
  for (int mode = 0; mode < 3; ++mode) {
    randomBeamCanvas.Divide(2, 2);
    std::array<TLine, 4> randomEtaLines, randomThresholdLines;
    for (int panel = 0; panel < 4; ++panel) {
      randomBeamCanvas.cd(panel+1);
      gPad->SetLeftMargin(.16); gPad->SetRightMargin(.18); gPad->SetBottomMargin(.15);
      auto& h = *randomMaps[mode][panel];
      const auto& status = randomStatus[mode][panel];
      h.SetMinimum(0.); h.SetMaximum(randomColorMaximum); h.Draw("colz");
      beamBoundary.Draw("L same"); // No cut: points outside the boundary retained.
      randomEtaLines[panel] = TLine(kEtaMass, beamPlotMin, kEtaMass, beamPlotMax);
      randomEtaLines[panel].SetLineColor(kRed); randomEtaLines[panel].SetLineWidth(1); randomEtaLines[panel].Draw("same");
      randomThresholdLines[panel] = TLine(.48, thresholdBeam, .60, thresholdBeam);
      randomThresholdLines[panel].SetLineColor(kRed); randomThresholdLines[panel].SetLineStyle(2); randomThresholdLines[panel].SetLineWidth(1); randomThresholdLines[panel].Draw("same");
      Info("random_pK_offset", "%s: entries %.0f; outside true boundary %lld; undefined %lld; outside axes %lld; mean mass shift %.4f MeV/c^2",
           scenarioTitle(mode,panel).Data(), h.GetEntries(), status.outsideBoundary, status.undefined, status.outsideAxes,
           h.GetEntries() ? status.massShiftSum/h.GetEntries() : 0.);
    }
    randomBeamCanvas.Print(outputPdf);
    randomBeamCanvas.Clear();
  }
  // Beam-Lambda opening-angle bias only. Rotate both daughters together so
  // invariant p-pi mass and Lambda momentum magnitude stay unchanged.
  const double angleStep=(angleOffsetMaxMrad-angleOffsetMinMrad)/(angleOffsetPoints-1);
  const double angleYmin=angleOffsetMinMrad-angleStep/2.;
  const double angleYmax=angleOffsetMaxMrad+angleStep/2.;
  TH2D angleLambda("eta_angle_lambda","K^{-}-#Lambda angle offset: reconstructed #Lambda;M(p#pi^{-}) [GeV/c^{2}];#delta#theta_{K#Lambda} [mrad]",
    400,1.08,1.16,angleOffsetPoints,angleYmin,angleYmax);
  TH2D angleMissing("eta_angle_missing","K^{-}-#Lambda angle offset: missing mass;M_{X} [GeV/c^{2}];#delta#theta_{K#Lambda} [mrad]",
    400,.48,.60,angleOffsetPoints,angleYmin,angleYmax);
  angleLambda.SetDirectory(nullptr);angleMissing.SetDirectory(nullptr);
  angleLambda.GetXaxis()->SetNdivisions(504);
  std::vector<Long64_t> angleUndefined(angleOffsetPoints,0);
  for (const auto& e:smearedEvents) {
    if (e.proton<0. || e.pion<0. || e.beam<0.) continue;
    const auto proton=Four(e.protonDirection*e.proton,kProtonMass);
    const auto pion=Four(e.pionDirection*e.pion,kPionMass);
    const auto lambda=proton+pion;
    const auto incoming=Four(TVector3(0.,0.,e.beam),kKaonMass)+Four(TVector3(),kProtonMass);
    auto axis=TVector3(0.,0.,1.).Cross(lambda.Vect());
    if (axis.Mag2()<1.e-24) axis=TVector3(0.,1.,0.);
    axis=axis.Unit();
    for (int j=0;j<angleOffsetPoints;++j) {
      const double offset=angleOffsetMinMrad+j*angleStep;
      auto rotatedProton=proton;auto rotatedPion=pion;
      rotatedProton.Rotate(offset*.001,axis);rotatedPion.Rotate(offset*.001,axis);
      const auto rotated=rotatedProton+rotatedPion;
      if (std::abs(rotated.M2()-lambda.M2())>1.e-10)
        throw std::runtime_error("Opening-angle rotation changed Lambda invariant mass");
      angleLambda.Fill(rotated.M(),offset);
      const auto missing=incoming-missingLambda(rotated);
      if (missing.M2()<0.) {++angleUndefined[j];continue;}
      angleMissing.Fill(std::sqrt(missing.M2()),offset);
    }
  }
  // Common rotation must preserve the entire reconstructed mass histogram.
  for (int x=1;x<=angleLambda.GetNbinsX();++x)
    for (int y=2;y<=angleOffsetPoints;++y)
      if (angleLambda.GetBinContent(x,y)!=angleLambda.GetBinContent(x,1))
        throw std::runtime_error("Angle-only offset changed Lambda mass distribution");
  TCanvas angleCanvas("c_eta_angle_offset","Opening-angle offset",1600,1400);
  angleCanvas.Divide(2,1);
  TLine angleLambdaLine(kLambdaMass,angleYmin,kLambdaMass,angleYmax);
  TLine angleEtaLine(kEtaMass,angleYmin,kEtaMass,angleYmax);
  for (int panel=0;panel<2;++panel) {
    angleCanvas.cd(panel+1);gPad->SetLeftMargin(.16);gPad->SetRightMargin(.18);gPad->SetBottomMargin(.15);
    (panel==0 ? angleLambda:angleMissing).Draw("colz");
    auto& line=panel==0 ? angleLambdaLine:angleEtaLine;
    line.SetLineColor(kRed);line.SetLineWidth(1);line.Draw("same");
  }
  angleCanvas.Print(outputPdf);
  for (int j=0;j<angleOffsetPoints;++j)
    if (j==0 || j==angleOffsetPoints/2 || j==angleOffsetPoints-1)
      Info("eta_angle_offset","offset %+.2f mrad: undefined M_X %lld",angleOffsetMinMrad+j*angleStep,angleUndefined[j]);
  angleCanvas.Print(TString(outputPdf) + "]");
}
