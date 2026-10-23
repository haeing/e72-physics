// Run: root -l -b -q phase_lambda_X.cc
// Produces linear and logarithmic PDFs from the same generated events.
// Configure reactions/decays below. All final-state p-pi- combinations are used.
// Isotropic phase-space toy, weighted by measured production cross sections.
// Forced decay channels: weights are CONDITIONAL on these decays, without BRs.
// No acceptance, polarization or material transport; Sigma+ vertices propagated.
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGenPhaseSpace.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TStyle.h>

namespace lambda_x {
// ----------------------- Parameters: edit here -----------------------------
enum class Particle { Proton, Neutron, PiMinus, PiPlus, Pi0, Gamma,
                      Lambda, Eta, Sigma0, SigmaPlus };
using P = Particle;
struct Reaction {
  bool enabled;
  const char* name; // ROOT title/legend notation
  std::vector<Particle> finalState; // production particles, before their decays
  double relativeWeight;
  Color_t color;
  const char* crossSectionGraph; // pK in MeV/c, sigma in mb
};
// Add/remove rows or change enabled. Each active row gets the same MC sample size.
// Weight per event = interpolated sigma(pK) * relativeWeight / MC sample size.
// These production cross sections do not include the forced decay BR factors.
constexpr bool kUseCrossSections = true;
const char* const kCrossSectionFile = "../lambda-yield/cs-info/total_cross_section.root";
const std::vector<Reaction> kReactions = {
  {true, "#Lambda#eta",             {P::Lambda, P::Eta},                  1., kRed+1, "g_Lambda_eta_pK"},
  {true, "#Lambda#pi^{0}",          {P::Lambda, P::Pi0},                  1., kBlue+1, "g_Lambda_pi0_pK"},
  {true, "#Sigma^{0}#pi^{0}",       {P::Sigma0, P::Pi0},                  1., kGreen+2, "g_Sigma0_pi0_pK"},
  {true, "#Sigma^{+}#pi^{-}",       {P::SigmaPlus, P::PiMinus},           1., kMagenta+1, "g_Sigmap_pim_pK"},
  {true, "#Sigma^{+}#pi^{-}#pi^{0}",{P::SigmaPlus, P::PiMinus, P::Pi0},   1., kOrange+7, "g_Sigmap_pim_pi0_pK"},
  {true, "#Lambda#pi^{+}#pi^{-}", {P::Lambda, P::PiPlus, P::PiMinus}, 1., kCyan+2, "g_Lambda_pip_pim_pK"},
  {true, "#Sigma^{0}#pi^{+}#pi^{-}", {P::Sigma0, P::PiPlus, P::PiMinus}, 1., kViolet+1, "g_Sigma0_pip_pim_pK"},
};
struct Decay { Particle mother; std::vector<Particle> daughters; };
// These channels are FORCED, not sampled by physical branching fractions.
// Removing a row makes that species stable; adding a row enables a new decay.
// Sigma+ -> n pi+ does not provide a proton and hence no p-pi- candidate here.
const std::vector<Decay> kDecays = {
  {P::Lambda,    {P::Proton, P::PiMinus}},
  {P::Sigma0,    {P::Lambda, P::Gamma}},
  {P::SigmaPlus, {P::Proton, P::Pi0}},
  {P::Eta,       {P::Pi0, P::Pi0, P::Pi0}}, // neutral mode only: eta -> 3 pi0 -> 6 gamma
  {P::Pi0,       {P::Gamma, P::Gamma}},
};
constexpr double kBeamMomentumMeV = 735.;
constexpr int kEventsPerReaction = 100000;
constexpr UInt_t kSeed = 20261009;
// PDG 2025: https://pdg.lbl.gov/2025/tables/rpp2025-tab-baryons-Sigma.pdf
constexpr double kSigmaPlusLifetimeSec = .8018e-10;
constexpr double kSpeedOfLightMmSec = 2.99792458e11;
constexpr double kVertexRangeMm = 10.; // displayed range; tails retained in overflow
constexpr bool kApplyResolution = true;
constexpr bool kConstrainLambdaMass = true; // false: missing uses Pp + Ppi
// Same magnitude-resolution prescription as lambda_missing_mass_resolution_study.cc.
constexpr double kBeamSigma = .001, kProtonSigma = .040, kPionSigma = .040;
// Offset = p_used - p_true, in MeV/c; directions are unchanged.
constexpr double kBeamOffsetMeV = 0., kProtonOffsetMeV = 0., kPionOffsetMeV = 0.;
constexpr double kLambdaMass = 1.115683, kKaonMass = .493677;
constexpr double kProtonMass = .938272, kPionMass = .139570;
// Lab-frame momentum axes; z is the incident kaon direction (GeV/c).
constexpr int kMomentumBins = 240;
constexpr double kMomentumMax = 1.2;
constexpr double kPionMomentumMax = .7; // pion momentum magnitude display limit
// Simultaneous selections on measured lab momenta and reconstructed p-pi mass.
constexpr double kSelectedPionMax = .25;
constexpr double kSelectedLambdaMomentumMin = .4, kSelectedLambdaMomentumMax = .6;
constexpr double kSelectedLambdaMassCenter = kLambdaMass, kSelectedLambdaMassHalfWidth = .015;
constexpr int kMassBins = 500;
constexpr double kLambdaAxisMin = 1.06, kLambdaAxisMax = 1.70;
constexpr double kMissingAxisMin = 0., kMissingAxisMax = .65;
// ---------------------------------------------------------------------------
// Masses in GeV/c^2; Sigma masses from PDG 2025 baryon summary.
// https://pdg.lbl.gov/2025/tables/rpp2025-sum-baryons.pdf
// https://pdg.lbl.gov/2025/tables/rpp2025-sum-mesons.pdf

double Mass(Particle particle) {
  switch (particle) {
    case P::Proton: return kProtonMass;
    case P::Neutron: return .939565;
    case P::PiMinus: case P::PiPlus: return kPionMass;
    case P::Pi0: return .134977;
    case P::Gamma: return 0.;
    case P::Lambda: return kLambdaMass;
    case P::Eta: return .547862;
    case P::Sigma0: return 1.192642;
    case P::SigmaPlus: return 1.189370;
  }
  throw std::runtime_error("Unknown particle");
}
// Use the TRUE beam momentum, never the reconstruction offset/smeared momentum.
// Linear interpolation only inside the supplied measurement range.
double CrossSection(TFile* file, const Reaction& reaction, double momentumMeV) {
  if (!kUseCrossSections) return 1.;
  auto* graph = dynamic_cast<TGraph*>(file->Get(reaction.crossSectionGraph));
  if (!graph || graph->GetN()<2)
    throw std::runtime_error(Form("Missing cross-section graph: %s",reaction.crossSectionGraph));
  double lower=graph->GetX()[0], upper=lower;
  for (int i=1; i<graph->GetN(); ++i) {
    lower=std::min(lower,graph->GetX()[i]);
    upper=std::max(upper,graph->GetX()[i]);
  }
  if (momentumMeV<lower || momentumMeV>upper)
    throw std::runtime_error(Form("%s: beam %.3f MeV/c outside cross-section data [%.3f, %.3f]; no extrapolation",reaction.name,momentumMeV,lower,upper));
  const double sigma=graph->Eval(momentumMeV);
  if (!std::isfinite(sigma) || sigma<0.)
    throw std::runtime_error("Invalid cross-section value");
  return sigma;
}
// Two-body Lambda eta reconstruction using ONLY measured beam magnitude,
// opening-angle cosine and fixed PDG masses; no daughter momentum magnitude.
std::vector<double> EtaMomentumSolutions(double beam,double cosine) {
  const double energy=std::sqrt(beam*beam+kKaonMass*kKaonMass)+kProtonMass;
  const double invariant=energy*energy-beam*beam;
  if (invariant<std::pow(kLambdaMass+Mass(P::Eta),2)) return {};
  const double B=beam*cosine;
  const double C=(invariant+kLambdaMass*kLambdaMass-Mass(P::Eta)*Mass(P::Eta))/2.;
  const double a=energy*energy-B*B,b=-2.*B*C,c=energy*energy*kLambdaMass*kLambdaMass-C*C;
  double discriminant=b*b-4.*a*c;
  if (discriminant < -1.e-12) return {};
  discriminant=std::max(0.,discriminant);
  std::vector<double> roots;
  for (double sign:{-1.,1.}) {
    const double momentum=(-b+sign*std::sqrt(discriminant))/(2.*a);
    if(momentum<0. || std::abs(energy*std::sqrt(momentum*momentum+kLambdaMass*kLambdaMass)-B*momentum-C)>1.e-8) continue;
    if(roots.empty() || std::abs(momentum-roots.back())>1.e-10) roots.push_back(momentum);
  }
  return roots;
}
const char* YieldAxis() {
  return "Counts"; // Bin contents retain the configured cross-section weights.
}
TLorentzVector Four(const TVector3& p, double mass) {
  TLorentzVector result; result.SetVectM(p, mass); return result;
}

// TGenPhaseSpace returns weights normalized to an upper bound of one.
// Rejection-unweight every two-/three-body step to avoid biasing three-body shapes.
void GenerateUnweighted(TGenPhaseSpace& generator, TRandom3& random) {
  for (int attempt = 0; attempt < 1000000; ++attempt) {
    const double weight = generator.Generate();
    if (!std::isfinite(weight) || weight < 0. || weight > 1.+1.e-9)
      throw std::runtime_error("Invalid normalized TGenPhaseSpace weight");
    if (random.Uniform() <= weight) return;
  }
  throw std::runtime_error("Phase-space unweighting did not converge");
}
struct Track {
  Particle species;
  TLorentzVector truth;
  int lambdaOrigin; // genuine pair only if p/pi share the same positive origin
  TVector3 vertex; // true production position in lab, mm
};
void DecayToFinal(Particle species, TLorentzVector four, int origin,
                  int& nextOrigin, std::vector<Track>& final,
                  TRandom3& random, TRandom3& vertexRandom, TVector3 vertex=TVector3(), int depth=0) {
  if (depth > 12) throw std::runtime_error("Decay configuration cycle/depth exceeded");
  const auto rule = std::find_if(kDecays.begin(), kDecays.end(),
    [&](const Decay& d) { return d.mother == species; });
  if (rule == kDecays.end()) { final.push_back({species, four, origin, vertex}); return; }
  if (rule->daughters.size() < 2 || rule->daughters.size() > 18)
    throw std::runtime_error("Decay must have 2 to 18 daughters");
  if (species == P::Lambda) origin = ++nextOrigin;
  // Origin fixed at (0,0,0); straight flight, no field/material or vertex smearing.
  // Only Sigma+ displacement is modeled, for the requested Sigma+ reaction pages.
  if (species == P::SigmaPlus)
    vertex += four.Vect() *
      (kSpeedOfLightMmSec * vertexRandom.Exp(kSigmaPlusLifetimeSec) / Mass(species));
  std::vector<double> masses;
  for (auto daughter : rule->daughters) masses.push_back(Mass(daughter));
  TGenPhaseSpace decay;
  if (!decay.SetDecay(four, masses.size(), masses.data()))
    throw std::runtime_error("Configured decay is kinematically impossible");
  GenerateUnweighted(decay, random);
  for (std::size_t i=0; i<rule->daughters.size(); ++i)
    DecayToFinal(rule->daughters[i], *decay.GetDecay(i), origin,
                 nextOrigin, final, random, vertexRandom, vertex, depth+1);
}
struct MeasuredTrack { const Track* track; TLorentzVector four; };
struct Histograms {
  Reaction reaction;
  std::unique_ptr<TH1D> lambda, missing, selectedMissing, genuineLambda, wrongLambda;
  std::unique_ptr<TH2D> phase, correlation, massMomentum, massPionMomentum, missingPionMomentum;
  std::unique_ptr<TH1D> separation;
  std::unique_ptr<TH2D> protonZX, pionZX, protonXY, pionXY, displacement;
  std::array<std::unique_ptr<TH1D>,4> momentum; // p, px, py, pz
  std::array<std::unique_ptr<TH1D>,4> protonMomentum; // measured tracks, once per proton
  std::array<std::unique_ptr<TH1D>,4> pionMomentum; // measured pi- tracks, once per pion
  Long64_t generated=0, candidates=0, genuine=0, negativeM2=0, invalidTracks=0;
  Histograms(const Reaction& r, int index): reaction(r) {
    lambda.reset(new TH1D(Form("lx_mass_%d",index),Form("%s: reconstructed #Lambda;M(p#pi^{-}) [GeV/c^{2}];%s",r.name,YieldAxis()),kMassBins,kLambdaAxisMin,kLambdaAxisMax));
    missing.reset(new TH1D(Form("lx_missing_%d",index),Form("%s: missing mass;M_{X} [GeV/c^{2}];%s",r.name,YieldAxis()),kMassBins,kMissingAxisMin,kMissingAxisMax));
    selectedMissing.reset(new TH1D(Form("lx_selected_missing_%d",index),Form("%s: selected missing mass;M_{X} [GeV/c^{2}];Counts",r.name),kMassBins,kMissingAxisMin,kMissingAxisMax));
    selectedMissing->SetLineColor(r.color);selectedMissing->SetLineWidth(2);
    selectedMissing->GetYaxis()->SetTitleOffset(1.6);
    genuineLambda.reset(new TH1D(Form("lx_genuine_%d",index),"",kMassBins,kLambdaAxisMin,kLambdaAxisMax));
    wrongLambda.reset(new TH1D(Form("lx_wrong_%d",index),"",kMassBins,kLambdaAxisMin,kLambdaAxisMax));
    phase.reset(new TH2D(Form("lx_phase_%d",index),Form("%s: p#pi^{-} candidate lab phase;#theta_{#Lambda}^{lab} [degree];p_{#Lambda}^{lab} [GeV/c]",r.name),180,0.,180.,200,0.,1.2));
    correlation.reset(new TH2D(Form("lx_correlation_%d",index),Form("%s: reconstructed masses;M_{X} [GeV/c^{2}];M(p#pi^{-}) [GeV/c^{2}]",r.name),300,kMissingAxisMin,kMissingAxisMax,300,kLambdaAxisMin,kLambdaAxisMax));
    massMomentum.reset(new TH2D(Form("lx_mass_momentum_%d",index),
      Form("%s: reconstructed #Lambda mass vs momentum;M(p#pi^{-}) [GeV/c^{2}];p_{#Lambda} [GeV/c]",r.name),
      300,kLambdaAxisMin,kLambdaAxisMax,kMomentumBins,0.,kMomentumMax));
    massMomentum->GetXaxis()->SetNdivisions(504);
    massMomentum->GetZaxis()->SetTitle("");
    massPionMomentum.reset(new TH2D(Form("lx_mass_pion_momentum_%d",index),
      Form("%s: reconstructed mass vs #pi^{-} momentum;M(p#pi^{-}) [GeV/c^{2}];p_{#pi^{-}} [GeV/c]",r.name),
      300,1.05,1.50,kMomentumBins,0.,kPionMomentumMax));
    massPionMomentum->GetXaxis()->SetNdivisions(504);
    massPionMomentum->GetZaxis()->SetTitle("");
    missingPionMomentum.reset(new TH2D(Form("lx_missing_pion_momentum_%d",index),
      Form("%s: missing mass vs #pi^{-} momentum;M_{X} [GeV/c^{2}];p_{#pi^{-}} [GeV/c]",r.name),
      300,kMissingAxisMin,kMissingAxisMax,kMomentumBins,0.,kPionMomentumMax));
    missingPionMomentum->GetZaxis()->SetTitle("");
    const std::vector<TH1*> all{lambda.get(),missing.get(),selectedMissing.get(),genuineLambda.get(),wrongLambda.get(),phase.get(),correlation.get(),massMomentum.get(),massPionMomentum.get(),missingPionMomentum.get()};
    for (auto* h:all) { h->SetDirectory(nullptr); h->Sumw2(); }
    separation.reset(new TH1D(Form("lx_separation_%d",index),Form("%s: p-#pi^{-} production separation;|r_{p}-r_{#pi^{-}}| [mm];Counts",r.name),240,0.,kVertexRangeMm));
    auto vertexPlot=[&](const char* name,const char* title) {
      auto result=std::unique_ptr<TH2D>(new TH2D(Form("lx_%s_%d",name,index),title,160,-kVertexRangeMm,kVertexRangeMm,160,-kVertexRangeMm,kVertexRangeMm));
      result->SetDirectory(nullptr);result->Sumw2();return result;
    };
    protonZX=vertexPlot("protonZX",Form("%s: lab production positions;z [mm];x [mm]",r.name));
    pionZX=vertexPlot("pionZX","");
    protonXY=vertexPlot("protonXY",Form("%s: lab production positions;x [mm];y [mm]",r.name));
    pionXY=vertexPlot("pionXY","");
    displacement.reset(new TH2D(Form("lx_displacement_%d",index),Form("%s: lab p-#pi^{-} separation;z_{p}-z_{#pi^{-}} [mm];Transverse separation [mm]",r.name),200,-kVertexRangeMm,kVertexRangeMm,160,0.,kVertexRangeMm));
    for (TH1* plot : {static_cast<TH1*>(separation.get()),static_cast<TH1*>(displacement.get())}) {plot->SetDirectory(nullptr);plot->Sumw2();}
    separation->SetLineColor(r.color);separation->SetLineWidth(2);
    separation->GetYaxis()->SetTitleOffset(1.6);
    const char* symbols[] = {"p_{#Lambda}","p_{x,#Lambda}","p_{y,#Lambda}","p_{z,#Lambda}"};
    for (int component=0; component<4; ++component) {
      momentum[component].reset(new TH1D(Form("lx_momentum_%d_%d",index,component),
        Form("%s: reconstructed #Lambda momentum;%s [GeV/c];Counts",r.name,symbols[component]),
        kMomentumBins,component==0 ? 0. : -kMomentumMax,kMomentumMax));
      auto* plot=momentum[component].get();
      plot->SetDirectory(nullptr);plot->Sumw2();plot->SetLineColor(r.color);plot->SetLineWidth(2);
      plot->GetYaxis()->SetTitleOffset(1.6);plot->GetXaxis()->SetNdivisions(504);
    }
    const char* protonSymbols[] = {"p_{p}","p_{x,p}","p_{y,p}","p_{z,p}"};
    for (int component=0; component<4; ++component) {
      protonMomentum[component].reset(new TH1D(Form("lx_proton_momentum_%d_%d",index,component),
        Form("%s: proton lab momentum;%s [GeV/c];Counts",r.name,protonSymbols[component]),
        kMomentumBins,component==0 ? 0. : -kMomentumMax,kMomentumMax));
      auto* plot=protonMomentum[component].get();
      plot->SetDirectory(nullptr);plot->Sumw2();plot->SetLineColor(r.color);plot->SetLineWidth(2);
      plot->GetYaxis()->SetTitleOffset(1.6);plot->GetXaxis()->SetNdivisions(504);
    }
    const char* pionSymbols[] = {"p_{#pi^{-}}","p_{x,#pi^{-}}","p_{y,#pi^{-}}","p_{z,#pi^{-}}"};
    for (int component=0; component<4; ++component) {
      pionMomentum[component].reset(new TH1D(Form("lx_pion_momentum_%d_%d",index,component),
        Form("%s: #pi^{-} lab momentum;%s [GeV/c];Counts",r.name,pionSymbols[component]),
        kMomentumBins,component==0 ? 0. : -kMomentumMax,
        component==0 ? kPionMomentumMax : kMomentumMax));
      auto* plot=pionMomentum[component].get();
      plot->SetDirectory(nullptr);plot->Sumw2();plot->SetLineColor(r.color);plot->SetLineWidth(2);
      plot->GetYaxis()->SetTitleOffset(1.6);plot->GetXaxis()->SetNdivisions(504);
    }
    phase->GetZaxis()->SetTitle(kUseCrossSections ? "Conditional mb/bin":"Candidates / event / bin");
    correlation->GetZaxis()->SetTitle(kUseCrossSections ? "Conditional mb/bin":"Candidates / event / bin");
    lambda->SetLineColor(r.color); lambda->SetLineWidth(2);
    missing->SetLineColor(r.color); missing->SetLineWidth(2);
    lambda->GetYaxis()->SetTitleOffset(1.6);
    missing->GetYaxis()->SetTitleOffset(1.6);
    genuineLambda->SetLineColor(kGreen+2); genuineLambda->SetLineStyle(2);
    wrongLambda->SetLineColor(kMagenta+1); wrongLambda->SetLineStyle(3);
    correlation->GetYaxis()->SetNdivisions(504);
    lambda->GetXaxis()->SetNdivisions(504);
  }
};
}

void phase_lambda_X(double beamMomentumMeV=lambda_x::kBeamMomentumMeV,
    int eventsPerReaction=lambda_x::kEventsPerReaction,
    bool constrainLambdaMass=lambda_x::kConstrainLambdaMass,
    bool applyResolution=lambda_x::kApplyResolution,
    const char* outputPdf="phase_lambda_X.pdf") {
  using namespace lambda_x;
  if (!std::isfinite(beamMomentumMeV) || beamMomentumMeV<0. || eventsPerReaction<=0)
    throw std::runtime_error("Invalid beam momentum or event count");
  gROOT->SetBatch(kTRUE); gStyle->SetOptStat(0); gStyle->SetOptTitle(1);
  TRandom3 random(kSeed);
  TRandom3 vertexRandom(kSeed+1); // independent stream preserves original momentum samples
  struct RandomGuard { TRandom* old; ~RandomGuard(){gRandom=old;} } guard{gRandom};
  gRandom=&random;
  std::unique_ptr<TFile> crossSections;
  if (kUseCrossSections) {
    crossSections.reset(TFile::Open(kCrossSectionFile,"READ"));
    if (!crossSections || crossSections->IsZombie())
      throw std::runtime_error(Form("Cannot read cross sections: %s",kCrossSectionFile));
  }
  const double beam=beamMomentumMeV/1000.;
  const auto target=Four(TVector3(),kProtonMass);
  auto initial=Four(TVector3(0.,0.,beam),kKaonMass)+target;
  // Each physical angular solution contributes once; this is solution multiplicity.
  auto etaMomentumComparison=std::unique_ptr<TH2D>(new TH2D("lx_eta_momentum_comparison",
    "#Lambda#eta: both momentum solutions;p_{#Lambda}^{beam,#theta} [GeV/c];p_{#Lambda}^{p#pi^{-}} [GeV/c]",
    300,.3,.8,300,.3,.8));
  etaMomentumComparison->SetDirectory(nullptr);etaMomentumComparison->Sumw2();
  Long64_t etaAngleCandidates=0,etaAngleNoSolution=0;
  std::vector<std::unique_ptr<Histograms>> histograms;
  for (std::size_t id=0; id<kReactions.size(); ++id) {
    const auto& reaction=kReactions[id];
    if (!reaction.enabled) continue;
    if (reaction.finalState.size()<2 || reaction.finalState.size()>18 ||
        !std::isfinite(reaction.relativeWeight) || reaction.relativeWeight<=0.)
      throw std::runtime_error("Invalid reaction configuration");
    std::vector<double> masses;
    for (auto particle:reaction.finalState) masses.push_back(Mass(particle));
    TGenPhaseSpace production;
    if (!production.SetDecay(initial,masses.size(),masses.data())) {
      Warning("phase_lambda_X","%s below threshold at %.1f MeV/c; skipped",reaction.name,beamMomentumMeV);
      continue;
    }
    const double sigma=CrossSection(crossSections.get(),reaction,beamMomentumMeV);
    const double eventWeight=sigma*reaction.relativeWeight/eventsPerReaction;
    Info("phase_lambda_X","%s at %.1f MeV/c: sigma %.6f mb, conditional event weight %.9g",reaction.name,beamMomentumMeV,sigma,eventWeight);
    histograms.emplace_back(new Histograms(reaction,id));
    auto& h=*histograms.back();
    for (int event=0; event<eventsPerReaction; ++event) {
      GenerateUnweighted(production,random);
      std::vector<Track> final;
      int nextOrigin=0;
      for (std::size_t i=0; i<reaction.finalState.size(); ++i)
        DecayToFinal(reaction.finalState[i],*production.GetDecay(i),-1,nextOrigin,final,random,vertexRandom);
      TLorentzVector finalSum;
      for (const auto& track:final) finalSum+=track.truth;
      const auto difference=finalSum-initial;
      if (std::abs(difference.E())>1.e-8 || difference.P()>1.e-8)
        throw std::runtime_error("Final-state four-momentum closure failed");
      ++h.generated;
      const double measuredBeam=beam*(applyResolution ? random.Gaus(1.,kBeamSigma):1.)+kBeamOffsetMeV/1000.;
      if (measuredBeam<0.) { ++h.invalidTracks; continue; }
      const auto initialUsed=Four(TVector3(0.,0.,measuredBeam),kKaonMass)+target;
      std::vector<MeasuredTrack> protons,pions;
      for (const auto& track:final) {
        if (track.species!=P::Proton && track.species!=P::PiMinus) continue;
        const bool proton=track.species==P::Proton;
        const double sigma=proton ? kProtonSigma:kPionSigma;
        const double offset=proton ? kProtonOffsetMeV:kPionOffsetMeV;
        const double magnitude=track.truth.P()*(applyResolution ? random.Gaus(1.,sigma):1.)+offset/1000.;
        if (magnitude<0.) { ++h.invalidTracks; continue; }
        MeasuredTrack measured{&track,Four(track.truth.Vect().Unit()*magnitude,Mass(track.species))};
        (proton ? protons:pions).push_back(measured);
      }
      // Fill once per measured proton, irrespective of the number of pion partners.
      for (const auto& proton:protons) {
        const double components[] = {proton.four.P(),proton.four.Px(),proton.four.Py(),proton.four.Pz()};
        for (int component=0; component<4; ++component)
          h.protonMomentum[component]->Fill(components[component],eventWeight);
      }
      for (const auto& pion:pions) {
        const double components[] = {pion.four.P(),pion.four.Px(),pion.four.Py(),pion.four.Pz()};
        for (int component=0; component<4; ++component)
          h.pionMomentum[component]->Fill(components[component],eventWeight);
      }
      for (const auto& proton:protons) for (const auto& pion:pions) {
        ++h.candidates;
        const bool genuine=proton.track->lambdaOrigin>0 &&
          proton.track->lambdaOrigin==pion.track->lambdaOrigin;
        h.genuine+=genuine;
        const auto lambda=proton.four+pion.four;
        const auto lambdaUsed=constrainLambdaMass ? Four(lambda.Vect(),kLambdaMass):lambda;
        const auto missing=initialUsed-lambdaUsed;
        const double weight=eventWeight;
        const bool lambdaEta=reaction.finalState.size()==2 &&
          std::find(reaction.finalState.begin(),reaction.finalState.end(),P::Lambda)!=reaction.finalState.end() &&
          std::find(reaction.finalState.begin(),reaction.finalState.end(),P::Eta)!=reaction.finalState.end();
        if(lambdaEta) {
          ++etaAngleCandidates;
          const auto solutions=EtaMomentumSolutions(measuredBeam,lambda.CosTheta());
          if(solutions.empty()) ++etaAngleNoSolution;
          else {
            for(double solution:solutions)etaMomentumComparison->Fill(solution,lambda.P(),weight);
            if(!applyResolution && kBeamOffsetMeV==0. && kProtonOffsetMeV==0. && kPionOffsetMeV==0.) {
              double closest=1.e100;for(double solution:solutions)closest=std::min(closest,std::abs(solution-lambda.P()));
              if(closest>1.e-7) throw std::runtime_error("Truth Lambda eta momentum matches neither angular solution");
            }
          }
        }
        const auto& rp=proton.track->vertex;
        const auto& rpi=pion.track->vertex;
        const auto delta=rp-rpi;
        h.protonZX->Fill(rp.Z(),rp.X(),weight);h.pionZX->Fill(rpi.Z(),rpi.X(),weight);
        h.protonXY->Fill(rp.X(),rp.Y(),weight);h.pionXY->Fill(rpi.X(),rpi.Y(),weight);
        h.separation->Fill(delta.Mag(),weight);h.displacement->Fill(delta.Z(),delta.Perp(),weight);
        const double components[] = {lambda.P(),lambda.Px(),lambda.Py(),lambda.Pz()};
        for (int component=0; component<4; ++component)
          h.momentum[component]->Fill(components[component],weight);
        h.lambda->Fill(lambda.M(),weight);
        h.massMomentum->Fill(lambda.M(),lambda.P(),weight);
        h.massPionMomentum->Fill(lambda.M(),pion.four.P(),weight);
        (genuine ? h.genuineLambda:h.wrongLambda)->Fill(lambda.M(),weight);
        h.phase->Fill(std::min(lambda.Theta()*180./std::acos(-1.),180.-1.e-12),lambda.P(),weight);
        if (!applyResolution && kProtonOffsetMeV==0. && kPionOffsetMeV==0. && genuine &&
            std::abs(lambda.M()-kLambdaMass)>1.e-8)
          throw std::runtime_error("Truth genuine Lambda mass closure failed");
        // Retain M^2 information through the count; no real M exists for M^2<0.
        // No Lambda mass window, MM>0.5 selection, or kinematic-region cut.
        if (missing.M2()<0.) { ++h.negativeM2; continue; }
        const double mass=std::sqrt(missing.M2());
        h.missing->Fill(mass,weight);
        if (pion.four.P()<kSelectedPionMax &&
            lambda.P()>kSelectedLambdaMomentumMin && lambda.P()<kSelectedLambdaMomentumMax &&
            std::abs(lambda.M()-kSelectedLambdaMassCenter)<kSelectedLambdaMassHalfWidth)
          h.selectedMissing->Fill(mass,weight);
        h.missingPionMomentum->Fill(mass,pion.four.P(),weight);
        h.correlation->Fill(mass,lambda.M(),weight);
      }
    }
    Info("phase_lambda_X","%s: generated %lld, all p-pi- candidates %lld, genuine Lambda %lld, non-Lambda pairs %lld, negative Mmiss^2 %lld, invalid momenta %lld",
         reaction.name,h.generated,h.candidates,h.genuine,h.candidates-h.genuine,h.negativeM2,h.invalidTracks);
    if (!h.candidates) Warning("phase_lambda_X","%s has no final p-pi- pairs with current decay configuration",reaction.name);
  }
  if (histograms.empty()) throw std::runtime_error("No enabled reaction is above threshold");

  auto totalLambda=std::unique_ptr<TH1D>(static_cast<TH1D*>(histograms[0]->lambda->Clone("lx_total_lambda")));
  auto totalMissing=std::unique_ptr<TH1D>(static_cast<TH1D*>(histograms[0]->missing->Clone("lx_total_missing")));
  auto totalPhase=std::unique_ptr<TH2D>(static_cast<TH2D*>(histograms[0]->phase->Clone("lx_total_phase")));
  auto totalCorrelation=std::unique_ptr<TH2D>(static_cast<TH2D*>(histograms[0]->correlation->Clone("lx_total_correlation")));
  for (std::size_t i=1;i<histograms.size();++i) {
    totalLambda->Add(histograms[i]->lambda.get()); totalMissing->Add(histograms[i]->missing.get());
    totalPhase->Add(histograms[i]->phase.get()); totalCorrelation->Add(histograms[i]->correlation.get());
  }
  auto totalMassMomentum=std::unique_ptr<TH2D>(static_cast<TH2D*>(histograms[0]->massMomentum->Clone("lx_total_mass_momentum")));
  for (std::size_t i=1;i<histograms.size();++i)
    totalMassMomentum->Add(histograms[i]->massMomentum.get());
  totalMassMomentum->SetTitle("All reactions: reconstructed #Lambda mass vs momentum;M(p#pi^{-}) [GeV/c^{2}];p_{#Lambda} [GeV/c]");
  auto totalMassPionMomentum=std::unique_ptr<TH2D>(static_cast<TH2D*>(histograms[0]->massPionMomentum->Clone("lx_total_mass_pion_momentum")));
  for (std::size_t i=1;i<histograms.size();++i)
    totalMassPionMomentum->Add(histograms[i]->massPionMomentum.get());
  totalMassPionMomentum->SetTitle("All reactions: reconstructed mass vs #pi^{-} momentum;M(p#pi^{-}) [GeV/c^{2}];p_{#pi^{-}} [GeV/c]");
  auto totalMissingPionMomentum=std::unique_ptr<TH2D>(static_cast<TH2D*>(histograms[0]->missingPionMomentum->Clone("lx_total_missing_pion_momentum")));
  for (std::size_t i=1;i<histograms.size();++i)
    totalMissingPionMomentum->Add(histograms[i]->missingPionMomentum.get());
  totalMissingPionMomentum->SetTitle("All reactions: missing mass vs #pi^{-} momentum;M_{X} [GeV/c^{2}];p_{#pi^{-}} [GeV/c]");
  std::array<std::unique_ptr<TH1D>,4> totalMomentum;
  for (int component=0; component<4; ++component) {
    totalMomentum[component].reset(static_cast<TH1D*>(histograms[0]->momentum[component]->Clone(Form("lx_total_momentum_%d",component))));
    for (std::size_t i=1; i<histograms.size(); ++i)
      totalMomentum[component]->Add(histograms[i]->momentum[component].get());
    totalMomentum[component]->SetLineColor(kBlack);
  }
  std::array<std::unique_ptr<TH1D>,4> totalProtonMomentum;
  for (int component=0; component<4; ++component) {
    totalProtonMomentum[component].reset(static_cast<TH1D*>(histograms[0]->protonMomentum[component]->Clone(Form("lx_total_proton_momentum_%d",component))));
    for (std::size_t i=1; i<histograms.size(); ++i)
      totalProtonMomentum[component]->Add(histograms[i]->protonMomentum[component].get());
    totalProtonMomentum[component]->SetLineColor(kBlack);
  }
  std::array<std::unique_ptr<TH1D>,4> totalPionMomentum;
  for (int component=0; component<4; ++component) {
    totalPionMomentum[component].reset(static_cast<TH1D*>(histograms[0]->pionMomentum[component]->Clone(Form("lx_total_pion_momentum_%d",component))));
    for (std::size_t i=1; i<histograms.size(); ++i)
      totalPionMomentum[component]->Add(histograms[i]->pionMomentum[component].get());
    totalPionMomentum[component]->SetLineColor(kBlack);
  }
  auto totalSelectedMissing=std::unique_ptr<TH1D>(static_cast<TH1D*>(histograms[0]->selectedMissing->Clone("lx_total_selected_missing")));
  for (std::size_t i=1;i<histograms.size();++i)
    totalSelectedMissing->Add(histograms[i]->selectedMissing.get());
  totalSelectedMissing->SetLineColor(kBlack);
  totalSelectedMissing->SetTitle("All reactions: selected missing mass;M_{X} [GeV/c^{2}];Counts");
  for (auto& h:histograms)
    Info("phase_lambda_X","%s: selected missing-mass weighted integral %.8g",h->reaction.name,h->selectedMissing->Integral(0,h->selectedMissing->GetNbinsX()+1));
  Info("eta_angle_momentum","Lambda eta candidates %lld; no physical angular solution %lld",etaAngleCandidates,etaAngleNoSolution);
  TCanvas canvas("c_phase_lambda_X","Multiple reaction reconstruction",1600,1400);
  TString logPdf(outputPdf);
  if (logPdf.EndsWith(".pdf")) logPdf.Remove(logPdf.Length()-4);
  logPdf += "_log.pdf";
  // Use a shared positive floor so faint reaction components remain visible.
  double smallestWeight=1.e100;
  for (auto& h:histograms) {
    for (TH1* plot : {static_cast<TH1*>(h->lambda.get()),static_cast<TH1*>(h->missing.get()),
                      static_cast<TH1*>(h->phase.get()),static_cast<TH1*>(h->correlation.get())})
      for (int bin=0; bin<plot->GetNcells(); ++bin)
        if (plot->GetBinContent(bin)>0.) smallestWeight=std::min(smallestWeight,plot->GetBinContent(bin));
  }
  if (smallestWeight==1.e100) smallestWeight=1.;
  for (bool logarithmic : {false,true}) {
  const TString pdf = logarithmic ? logPdf : TString(outputPdf);
  auto setScale=[&](TH1* plot) {
    plot->SetMinimum(logarithmic ? smallestWeight*.5 : 0.);
    plot->SetMaximum(-1111.); // Let ROOT select the upper range for each plot.
  };
  for (TH1* plot : {static_cast<TH1*>(totalLambda.get()),static_cast<TH1*>(totalMissing.get()),
                    static_cast<TH1*>(totalPhase.get()),static_cast<TH1*>(totalCorrelation.get())}) setScale(plot);
  for (auto& h:histograms)
    for (TH1* plot : {static_cast<TH1*>(h->lambda.get()),static_cast<TH1*>(h->missing.get()),
                      static_cast<TH1*>(h->phase.get()),static_cast<TH1*>(h->correlation.get())}) setScale(plot);
  canvas.Print(pdf+"[");
  auto margins=[&](){
    gPad->SetLeftMargin(.17);gPad->SetRightMargin(.18);gPad->SetBottomMargin(.15);
    gPad->SetLogy(logarithmic && gPad->GetNumber()<=2);
    gPad->SetLogz(logarithmic && gPad->GetNumber()>=3);
  };
  canvas.Divide(2,2);
  TLegend lambdaLegend(.53,.52,.87,.88),missingLegend(.48,.52,.82,.88);
  for (auto* legend:{&lambdaLegend,&missingLegend}) {legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.030);}
  canvas.cd(1); margins(); totalLambda->SetTitle(Form("All reactions, beam %.0f MeV/c;M(p#pi^{-}) [GeV/c^{2}];%s",beamMomentumMeV,YieldAxis()));
  totalLambda->SetLineColor(kBlack); totalLambda->Draw("hist"); lambdaLegend.AddEntry(totalLambda.get(),"Total","l");
  for (auto& h:histograms) {h->lambda->Draw("hist same");lambdaLegend.AddEntry(h->lambda.get(),h->reaction.name,"l");}
  lambdaLegend.Draw();
  canvas.cd(2); margins(); totalMissing->SetTitle(Form("All reactions, #Lambda %s;M_{X} [GeV/c^{2}];%s",constrainLambdaMass ? "PDG mass":"p#pi^{-} four-vector",YieldAxis()));
  totalMissing->SetLineColor(kBlack);totalMissing->Draw("hist");missingLegend.AddEntry(totalMissing.get(),"Total","l");
  for (auto& h:histograms) {h->missing->Draw("hist same");missingLegend.AddEntry(h->missing.get(),h->reaction.name,"l");}
  missingLegend.Draw();
  canvas.cd(3);margins();totalPhase->SetTitle("All reactions: p#pi^{-} lab phase;#theta_{#Lambda}^{lab} [degree];p_{#Lambda}^{lab} [GeV/c]");totalPhase->Draw("colz");
  canvas.cd(4);margins();totalCorrelation->SetTitle("All reactions: reconstructed masses;M_{X} [GeV/c^{2}];M(p#pi^{-}) [GeV/c^{2}]");totalCorrelation->Draw("colz");
  canvas.Print(pdf);canvas.Clear();
  for (auto& h:histograms) {
    canvas.Divide(2,2);
    TLegend pairLegend(.48,.67,.82,.88);pairLegend.SetBorderSize(0);pairLegend.SetFillStyle(0);pairLegend.SetTextSize(.032);
    canvas.cd(1);margins();h->lambda->Draw("hist");h->genuineLambda->Draw("hist same");h->wrongLambda->Draw("hist same");
    pairLegend.AddEntry(h->lambda.get(),"All p#pi^{-}","l");pairLegend.AddEntry(h->genuineLambda.get(),"Same #Lambda decay","l");pairLegend.AddEntry(h->wrongLambda.get(),"Other pairing","l");pairLegend.Draw();
    canvas.cd(2);margins();h->missing->Draw("hist");
    canvas.cd(3);margins();h->phase->Draw("colz");
    canvas.cd(4);margins();h->correlation->Draw("colz");
    canvas.Print(pdf);canvas.Clear();
  }
  // Momentum distributions include every reconstructed p-pi candidate,
  // including candidates with negative missing M^2, before any missing-mass cut.
  const char* momentumSymbols[] = {"p_{#Lambda}","p_{x,#Lambda}","p_{y,#Lambda}","p_{z,#Lambda}"};
  canvas.Divide(2,2);
  TLegend momentumLegend(.60,.52,.94,.88);
  momentumLegend.SetBorderSize(0);momentumLegend.SetFillStyle(0);momentumLegend.SetTextSize(.030);
  momentumLegend.AddEntry(totalMomentum[0].get(),"Total","l");
  for (auto& h:histograms) momentumLegend.AddEntry(h->momentum[0].get(),h->reaction.name,"l");
  for (int component=0; component<4; ++component) {
    canvas.cd(component+1);margins();gPad->SetLogz(false);gPad->SetLogy(logarithmic);
    auto* plot=totalMomentum[component].get();setScale(plot);
    plot->SetTitle(Form("All reactions, beam %.0f MeV/c;%s [GeV/c];Counts",beamMomentumMeV,momentumSymbols[component]));
    plot->Draw("hist");
    for (auto& h:histograms) {
      setScale(h->momentum[component].get());h->momentum[component]->Draw("hist same");
    }
    momentumLegend.Draw();
  }
  canvas.Print(pdf);canvas.Clear();
  for (auto& h:histograms) {
    canvas.Divide(2,2);
    for (int component=0; component<4; ++component) {
      canvas.cd(component+1);margins();gPad->SetLogz(false);gPad->SetLogy(logarithmic);
      setScale(h->momentum[component].get());h->momentum[component]->Draw("hist");
    }
    canvas.Print(pdf);canvas.Clear();
  }
  for (auto& h:histograms) {
    if (std::find(h->reaction.finalState.begin(),h->reaction.finalState.end(),P::SigmaPlus)==h->reaction.finalState.end()) continue;
    canvas.Divide(2,2);
    TLegend vertexLegend(.52,.70,.81,.87);
    vertexLegend.SetBorderSize(0);vertexLegend.SetFillStyle(0);vertexLegend.SetTextSize(.032);
    for (auto* plot:{h->protonZX.get(),h->protonXY.get()}) {plot->SetMarkerColor(kBlue+1);plot->SetMarkerStyle(7);}
    for (auto* plot:{h->pionZX.get(),h->pionXY.get()}) {plot->SetMarkerColor(kRed+1);plot->SetMarkerStyle(29);plot->SetMarkerSize(1.4);}
    vertexLegend.AddEntry(h->protonZX.get(),"Proton production","p");vertexLegend.AddEntry(h->pionZX.get(),"#pi^{-} production","p");
    canvas.cd(1);margins();gPad->SetLogy(false);gPad->SetLogz(false);
    h->protonZX->Draw("SCAT");h->pionZX->Draw("SCAT SAME");vertexLegend.Draw();
    canvas.cd(2);margins();gPad->SetLogy(false);gPad->SetLogz(false);
    h->protonXY->Draw("SCAT");h->pionXY->Draw("SCAT SAME");vertexLegend.Draw();
    canvas.cd(3);margins();gPad->SetLogz(false);gPad->SetLogy(logarithmic);setScale(h->separation.get());h->separation->Draw("hist");
    canvas.cd(4);margins();gPad->SetLogy(false);gPad->SetLogz(logarithmic);setScale(h->displacement.get());h->displacement->Draw("colz");
    canvas.Print(pdf);canvas.Clear();
  }
  const char* protonMomentumSymbols[] = {"p_{p}","p_{x,p}","p_{y,p}","p_{z,p}"};
  canvas.Divide(2,2);
  TLegend protonMomentumLegend(.48,.52,.82,.88);
  protonMomentumLegend.SetBorderSize(0);protonMomentumLegend.SetFillStyle(0);protonMomentumLegend.SetTextSize(.030);
  protonMomentumLegend.AddEntry(totalProtonMomentum[0].get(),"Total","l");
  for (auto& h:histograms) protonMomentumLegend.AddEntry(h->protonMomentum[0].get(),h->reaction.name,"l");
  for (int component=0; component<4; ++component) {
    canvas.cd(component+1);margins();gPad->SetLogz(false);gPad->SetLogy(logarithmic);
    auto* plot=totalProtonMomentum[component].get();setScale(plot);
    plot->SetTitle(Form("All reactions: proton, beam %.0f MeV/c;%s [GeV/c];Counts",beamMomentumMeV,protonMomentumSymbols[component]));
    plot->Draw("hist");
    for (auto& h:histograms) {
      setScale(h->protonMomentum[component].get());h->protonMomentum[component]->Draw("hist same");
    }
    protonMomentumLegend.Draw();
  }
  canvas.Print(pdf);canvas.Clear();
  for (auto& h:histograms) {
    canvas.Divide(2,2);
    for (int component=0; component<4; ++component) {
      canvas.cd(component+1);margins();gPad->SetLogz(false);gPad->SetLogy(logarithmic);
      setScale(h->protonMomentum[component].get());h->protonMomentum[component]->Draw("hist");
    }
    canvas.Print(pdf);canvas.Clear();
  }
  const char* pionMomentumSymbols[] = {"p_{#pi^{-}}","p_{x,#pi^{-}}","p_{y,#pi^{-}}","p_{z,#pi^{-}}"};
  canvas.Divide(2,2);
  TLegend pionMomentumLegend(.48,.52,.82,.88);
  pionMomentumLegend.SetBorderSize(0);pionMomentumLegend.SetFillStyle(0);pionMomentumLegend.SetTextSize(.030);
  pionMomentumLegend.AddEntry(totalPionMomentum[0].get(),"Total","l");
  for (auto& h:histograms) pionMomentumLegend.AddEntry(h->pionMomentum[0].get(),h->reaction.name,"l");
  for (int component=0; component<4; ++component) {
    canvas.cd(component+1);margins();gPad->SetLogz(false);gPad->SetLogy(logarithmic);
    auto* plot=totalPionMomentum[component].get();setScale(plot);
    plot->SetTitle(Form("All reactions: #pi^{-}, beam %.0f MeV/c;%s [GeV/c];Counts",beamMomentumMeV,pionMomentumSymbols[component]));
    plot->Draw("hist");
    for (auto& h:histograms) {
      setScale(h->pionMomentum[component].get());h->pionMomentum[component]->Draw("hist same");
    }
    pionMomentumLegend.Draw();
  }
  canvas.Print(pdf);canvas.Clear();
  for (auto& h:histograms) {
    canvas.Divide(2,2);
    for (int component=0; component<4; ++component) {
      canvas.cd(component+1);margins();gPad->SetLogz(false);gPad->SetLogy(logarithmic);
      setScale(h->pionMomentum[component].get());h->pionMomentum[component]->Draw("hist");
    }
    canvas.Print(pdf);canvas.Clear();
  }
  // All candidates, before any missing-mass selection; shared momentum smearing.
  auto drawMassMomentum=[&](TH2D* plot) {
    canvas.Divide(1,1);canvas.cd(1);margins();
    gPad->SetLogy(false);gPad->SetLogz(logarithmic);
    setScale(plot);plot->Draw("colz");canvas.Print(pdf);canvas.Clear();
  };
  drawMassMomentum(totalMassMomentum.get());
  for (auto& h:histograms) drawMassMomentum(h->massMomentum.get());
  drawMassMomentum(totalMassPionMomentum.get());
  for (auto& h:histograms) drawMassMomentum(h->massPionMomentum.get());
  drawMassMomentum(totalMissingPionMomentum.get());
  for (auto& h:histograms) drawMassMomentum(h->missingPionMomentum.get());
  canvas.Divide(1,1);canvas.cd(1);margins();gPad->SetLogy(logarithmic);gPad->SetLogz(false);
  TLegend selectedLegend(.60,.52,.94,.88);
  selectedLegend.SetBorderSize(0);selectedLegend.SetFillStyle(0);selectedLegend.SetTextSize(.030);
  setScale(totalSelectedMissing.get());totalSelectedMissing->Draw("hist");
  selectedLegend.AddEntry(totalSelectedMissing.get(),"Total","l");
  for (auto& h:histograms) {
    setScale(h->selectedMissing.get());h->selectedMissing->Draw("hist same");
    selectedLegend.AddEntry(h->selectedMissing.get(),h->reaction.name,"l");
  }
  selectedLegend.Draw();canvas.Print(pdf);canvas.Clear();
  for (auto& h:histograms) {
    canvas.Divide(1,1);canvas.cd(1);margins();gPad->SetLogy(logarithmic);gPad->SetLogz(false);
    setScale(h->selectedMissing.get());h->selectedMissing->Draw("hist");canvas.Print(pdf);canvas.Clear();
  }
  if(etaAngleCandidates) {
    // Never select a branch using the reconstructed momentum being compared.
    // Events beyond the angular kinematic boundary are counted and not forced to a root.
    {
      auto& plot=etaMomentumComparison;
      canvas.Divide(1,1);canvas.cd(1);margins();gPad->SetLogy(false);gPad->SetLogz(logarithmic);
      setScale(plot.get());plot->Draw("colz");
      TLine diagonal(.3,.3,.8,.8);
      diagonal.SetLineColor(kRed);diagonal.SetLineWidth(1);diagonal.SetLineStyle(2);diagonal.Draw("same");
      canvas.Print(pdf);canvas.Clear();
    }
  }
  canvas.Print(pdf+"]");
  }
}
