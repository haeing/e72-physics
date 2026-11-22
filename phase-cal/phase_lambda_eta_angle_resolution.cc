// Run: root -l -b -q phase_lambda_eta_angle_resolution.cc
// Reuse phase_lambda_X settings and helpers; do not run its other plots.
#include "phase_lambda_X.cc"

namespace eta_angle_resolution {
// Each row: relative proton, pion, beam momentum sigma.
const std::vector<std::array<double,3>> kResolutions = {
  {{0.,0.,0.}}, {{0.,0.,.001}}, {{.01,.01,.001}},
  {{.02,.02,.001}}, {{.04,.04,.001}}, {{.06,.06,.001}}
};
}

void phase_lambda_eta_angle_resolution(double beamMomentumMeV=lambda_x::kBeamMomentumMeV,
    int events=lambda_x::kEventsPerReaction,
    const char* outputPdf="phase_lambda_eta_angle_resolution.pdf") {
  using namespace lambda_x;
  using namespace eta_angle_resolution;
  if(!std::isfinite(beamMomentumMeV) || beamMomentumMeV<=0. || events<=0 || kResolutions.empty())
    throw std::runtime_error("Invalid resolution scan configuration");
  gROOT->SetBatch(true);gStyle->SetOptStat(0);
  TRandom3 random(kSeed),noise(kSeed+2);
  struct Guard {TRandom* old;~Guard(){gRandom=old;}} guard{gRandom};gRandom=&random;
  const double beam=beamMomentumMeV/1000.;
  const auto target=Four(TVector3(),kProtonMass);
  auto initial=Four(TVector3(0.,0.,beam),kKaonMass)+target;
  double productionMasses[]={kLambdaMass,Mass(P::Eta)},decayMasses[]={kProtonMass,kPionMass};
  TGenPhaseSpace production;
  if(!production.SetDecay(initial,2,productionMasses))throw std::runtime_error("Lambda eta below threshold");
  std::vector<std::unique_ptr<TH2D>> plots;
  std::vector<Long64_t> invalid(kResolutions.size(),0);
  for(std::size_t i=0;i<kResolutions.size();++i) {
    const auto& r=kResolutions[i];
    for(double sigma:r)if(!std::isfinite(sigma)||sigma<0.)throw std::runtime_error("Invalid relative momentum sigma");
    plots.emplace_back(new TH2D(Form("eta_angle_resolution_%zu",i),
      Form("#sigma_{p}/p: p %.1f%%, #pi^{-} %.1f%%, K^{-} %.1f%%;#theta_{#Lambda}^{lab} [degree];cos#theta_{#Lambda}^{CM}",100.*r[0],100.*r[1],100.*r[2]),
      360,0.,10.,300,-1.,1.));
    plots.back()->SetDirectory(nullptr);
  }
  // Same true events and standard-normal draws at every resolution.
  // No momentum offsets; genuine Lambda daughters only. Lambda energy uses PDG mass.
  for(int event=0;event<events;++event) {
    GenerateUnweighted(production,random);auto lambda=*production.GetDecay(0);
    TGenPhaseSpace decay;
    if(!decay.SetDecay(lambda,2,decayMasses))throw std::runtime_error("Lambda decay failed");
    GenerateUnweighted(decay,random);
    const auto proton=*decay.GetDecay(0),pion=*decay.GetDecay(1);
    const double zp=noise.Gaus(),zpi=noise.Gaus(),zk=noise.Gaus();
    for(std::size_t i=0;i<kResolutions.size();++i) {
      const auto& r=kResolutions[i];
      const double pp=proton.P()*(1.+r[0]*zp),ppi=pion.P()*(1.+r[1]*zpi),pk=beam*(1.+r[2]*zk);
      if(pp<0.||ppi<0.||pk<0.){++invalid[i];continue;}
      const auto reconstructed=Four(proton.Vect().Unit()*pp,kProtonMass)+Four(pion.Vect().Unit()*ppi,kPionMass);
      auto usedLambda=Four(reconstructed.Vect(),kLambdaMass);
      const auto incoming=Four(TVector3(0.,0.,pk),kKaonMass)+target;
      usedLambda.Boost(-incoming.BoostVector());
      if(r[0]==0. && r[1]==0. && r[2]==0.) {
        auto truth=lambda;truth.Boost(-initial.BoostVector());
        if(std::abs(usedLambda.CosTheta()-truth.CosTheta())>1.e-8)
          throw std::runtime_error("Zero-resolution CM angle closure failed");
      }
      plots[i]->Fill(reconstructed.Theta()*180./std::acos(-1.),std::min(usedLambda.CosTheta(),1.-1.e-12));
    }
  }
  double maximum=0.;for(auto& h:plots)maximum=std::max(maximum,h->GetMaximum());
  TString logPdf(outputPdf);if(logPdf.EndsWith(".pdf"))logPdf.Remove(logPdf.Length()-4);logPdf+="_log.pdf";
  TCanvas canvas("c_eta_angle_resolution","Lambda eta angle resolution comparison",1800,1200);
  for(bool logarithmic:{false,true}) {
    canvas.Clear();canvas.Divide(3,(kResolutions.size()+2)/3);
    for(std::size_t i=0;i<plots.size();++i) {
      canvas.cd(i+1);gPad->SetLeftMargin(.16);gPad->SetRightMargin(.18);gPad->SetBottomMargin(.16);
      gPad->SetLogz(logarithmic);plots[i]->SetMinimum(logarithmic ? .5:0.);plots[i]->SetMaximum(maximum);
      plots[i]->Draw("colz");
    }
    canvas.Print(logarithmic ? logPdf:TString(outputPdf));
  }
  for(std::size_t i=0;i<plots.size();++i)
    Info("eta_angle_resolution","panel %zu: entries %.0f, inside axes %.0f, invalid momenta %lld",i+1,plots[i]->GetEntries(),plots[i]->Integral(),invalid[i]);
}
