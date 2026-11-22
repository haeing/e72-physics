// Run: root -l -b -q lambda-pipi-phase.cc
// K- p -> Lambda pi+ pi-, Lambda -> p pi-; lab frame, no resolution.
// Isotropic three-body phase space; no resonance model or acceptance.
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <TCanvas.h>
#include <TGenPhaseSpace.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>

namespace lambda_pipi_phase_config {
constexpr double kBeamMomentumGeV = .735;
constexpr int kEvents = 100000;
constexpr UInt_t kSeed = 20261026;
constexpr double kKaonMass = .493677, kProtonMass = .938272;
constexpr double kLambdaMass = 1.115683, kPionMass = .139570;
constexpr double kMomentumMax = 1.2;
void Generate(TGenPhaseSpace& generator,TRandom3& random) {
  // Rejection unweighting is essential for a three-body phase-space shape.
  for(int attempt=0;attempt<1000000;++attempt) {
    const double weight=generator.Generate();
    if(!std::isfinite(weight) || weight<0. || weight>1.+1.e-9)
      throw std::runtime_error("Invalid phase-space weight");
    if(random.Uniform()<=weight)return;
  }
  throw std::runtime_error("Phase-space unweighting failed");
}
}

void lambda_pipi_phase(double beamMomentumGeV=lambda_pipi_phase_config::kBeamMomentumGeV,
    int events=lambda_pipi_phase_config::kEvents,
    const char* outputPdf="lambda-pipi-phase.pdf") {
  using namespace lambda_pipi_phase_config;
  if(!std::isfinite(beamMomentumGeV) || beamMomentumGeV<0. || events<=0)
    throw std::runtime_error("Invalid beam momentum or event count");
  gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);gStyle->SetPalette(1);
  TRandom3 random(kSeed);
  struct Guard {TRandom* old;~Guard(){gRandom=old;}} guard{gRandom};gRandom=&random;
  const TLorentzVector target(0.,0.,0.,kProtonMass);
  const TLorentzVector beam(0.,0.,beamMomentumGeV,std::hypot(beamMomentumGeV,kKaonMass));
  auto initial=beam+target;
  double productionMasses[]={kLambdaMass,kPionMass,kPionMass};
  double decayMasses[]={kProtonMass,kPionMass};
  TGenPhaseSpace production;
  if(!production.SetDecay(initial,3,productionMasses))
    throw std::runtime_error("Beam below Lambda pi+ pi- threshold");
  const char* names[]={"lambda","pip_primary","pim_primary","proton_lambda","pim_lambda"};
  const char* titles[]={"#Lambda","#pi^{+} from reaction","#pi^{-} from reaction","p from #Lambda","#pi^{-} from #Lambda"};
  const char* components[]={"|#vec{p}|","p_{x}","p_{y}","p_{z}"};
  std::array<std::unique_ptr<TH2D>,5> phase;
  std::array<std::array<std::unique_ptr<TH1D>,4>,5> momentum;
  for(int i=0;i<5;++i) {
    phase[i].reset(new TH2D(Form("lpp_phase_%s",names[i]),
      Form("%s in lab;cos#theta;momentum [GeV/#it{c}]",titles[i]),200,-1.,1.,240,0.,kMomentumMax));
    phase[i]->SetDirectory(nullptr);
    for(int j=0;j<4;++j) {
      momentum[i][j].reset(new TH1D(Form("lpp_momentum_%s_%d",names[i],j),
        Form("%s: %s in lab;%s [GeV/#it{c}];Counts",titles[i],components[j],components[j]),
        240,j==0 ? 0. : -kMomentumMax,kMomentumMax));
      momentum[i][j]->SetDirectory(nullptr);momentum[i][j]->SetLineWidth(2);
    }
  }
  for(int event=0;event<events;++event) {
    Generate(production,random);
    auto lambda=*production.GetDecay(0);
    const auto pip=*production.GetDecay(1),pim=*production.GetDecay(2);
    TGenPhaseSpace decay;
    if(!decay.SetDecay(lambda,2,decayMasses))throw std::runtime_error("Lambda decay failed");
    Generate(decay,random);
    const auto proton=*decay.GetDecay(0),daughterPion=*decay.GetDecay(1);
    const auto closure=proton+daughterPion+pip+pim-initial;
    if(std::abs(closure.E())>1.e-9 || closure.P()>1.e-9)
      throw std::runtime_error("Final four-momentum closure failed");
    const std::array<TLorentzVector,5> particles{{lambda,pip,pim,proton,daughterPion}};
    for(int i=0;i<5;++i) {
      const auto& p=particles[i];phase[i]->Fill(p.CosTheta(),p.P());
      const double values[]={p.P(),p.Px(),p.Py(),p.Pz()};
      for(int j=0;j<4;++j)momentum[i][j]->Fill(values[j]);
    }
  }
  TCanvas canvas("c_lambda_pipi_phase","Lambda pi+ pi- lab phase",1100,1000);
  canvas.Print(TString(outputPdf)+"[");
  for(int i=0;i<5;++i) {
    canvas.SetLeftMargin(.15);canvas.SetRightMargin(.16);canvas.SetBottomMargin(.14);
    phase[i]->Draw("colz");canvas.Print(outputPdf);canvas.Clear();
    Info("lambda_pipi_phase","%s: generated %.0f, plotted %.0f",names[i],phase[i]->GetEntries(),phase[i]->Integral());
  }
  for(int i=0;i<5;++i) {
    canvas.Divide(2,2);
    for(int j=0;j<4;++j) {
      canvas.cd(j+1);gPad->SetLeftMargin(.15);gPad->SetBottomMargin(.14);
      momentum[i][j]->Draw("hist");
    }
    canvas.Print(outputPdf);canvas.Clear();
  }
  canvas.Print(TString(outputPdf)+"]");
}
