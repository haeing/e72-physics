// Run: root -l -b -q phase_lambda_pip_pim.cc
// K- p -> Lambda pi+ pi-: unweighted, isotropic three-body phase space.
// Lab theta is relative to the incident beam (+z). No acceptance or dynamics.
// The plotted pi- is the production pi-, not the Lambda-decay pi-.
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <TCanvas.h>
#include <TGenPhaseSpace.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TLine.h>
#include <vector>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>

namespace lambda_pip_pim {
// ---------------- Parameters ----------------
constexpr double kBeamMomentumMeV = 735.;
constexpr int kEvents = 100000;
constexpr UInt_t kSeed = 20261012;
constexpr bool kApplyResolution = false;
constexpr double kProtonSigma = .04, kPionSigma = .04;
constexpr double kKaonMass = .493677, kProtonMass = .938272;
constexpr double kPionMass = .139570, kLambdaMass = 1.115683;
constexpr double kEtaMass = .547862;
const std::vector<double> kEtaLabAnglesDeg = {2.,4.,6.,7.};
constexpr double kEtaAngleHalfWidthDeg = .15;
constexpr double kEtaThetaMaxDeg = 15.;
constexpr int kThetaBins = 180, kMomentumBins = 240;
constexpr double kLambdaMomentumMax = 1.2, kPionMomentumMax = .7;
// Lambda is reconstructed from its p pi- decay when resolution is enabled.
// Each pion magnitude is smeared by 4%; track directions remain unchanged.
// Incoming beam is fixed for production; no reconstructed beam enters these plots.
// --------------------------------------------
TLorentzVector Four(const TVector3& p, double mass) {
  TLorentzVector out;out.SetVectM(p,mass);return out;
}
void Unweight(TGenPhaseSpace& generator,TRandom3& random) {
  for (int attempt=0;attempt<1000000;++attempt) {
    const double weight=generator.Generate();
    if (!std::isfinite(weight) || weight<0. || weight>1.+1.e-9)
      throw std::runtime_error("Invalid phase-space weight");
    if (random.Uniform()<=weight) return;
  }
  throw std::runtime_error("Unweighting failed");
}
// Solve E_initial E_Lambda - p_beam p_Lambda cos(theta) = (s+mL^2-meta^2)/2.
std::vector<double> EtaLabSolutions(const TLorentzVector& initial,double thetaDeg) {
  const double A=initial.E(), B=initial.P()*std::cos(thetaDeg*std::acos(-1.)/180.);
  const double C=(initial.M2()+kLambdaMass*kLambdaMass-kEtaMass*kEtaMass)/2.;
  const double a=A*A-B*B,b=-2.*B*C,c=A*A*kLambdaMass*kLambdaMass-C*C;
  const double discriminant=b*b-4.*a*c;
  std::vector<double> solutions;
  if(discriminant<0.)return solutions;
  for(double sign:{-1.,1.}) {
    const double p=(-b+sign*std::sqrt(discriminant))/(2.*a);
    if(p>=0. && std::abs(A*std::sqrt(p*p+kLambdaMass*kLambdaMass)-B*p-C)<1.e-8)
      solutions.push_back(p);
  }
  return solutions;
}
TLorentzVector Measure(const TLorentzVector& truth,double mass,double sigma,
                      bool apply,TRandom3& random) {
  const double p=truth.P()*(apply ? random.Gaus(1.,sigma):1.);
  if (p<0.) throw std::runtime_error("Negative measured momentum");
  return Four(truth.Vect().Unit()*p,mass);
}
}

void phase_lambda_pip_pim(double beamMomentumMeV=lambda_pip_pim::kBeamMomentumMeV,
    int events=lambda_pip_pim::kEvents,
    bool applyResolution=lambda_pip_pim::kApplyResolution,
    const char* outputPdf="phase_lambda_pip_pim.pdf") {
  using namespace lambda_pip_pim;
  if (!std::isfinite(beamMomentumMeV) || beamMomentumMeV<0. || events<=0)
    throw std::runtime_error("Invalid beam or event count");
  gROOT->SetBatch(true);gStyle->SetOptStat(0);gStyle->SetOptTitle(1);
  TRandom3 random(kSeed);
  struct Guard {TRandom* old;~Guard(){gRandom=old;}} guard{gRandom};gRandom=&random;
  auto initial=Four(TVector3(0.,0.,beamMomentumMeV/1000.),kKaonMass)+Four(TVector3(),kProtonMass);
  double masses[]={kLambdaMass,kPionMass,kPionMass};
  TGenPhaseSpace production;
  if (!production.SetDecay(initial,3,masses))
    throw std::runtime_error("Beam below Lambda pi+ pi- production threshold");
  const char* particles[]={"#Lambda","#pi^{+}","#pi^{-}"};
  std::array<std::unique_ptr<TH2D>,3> plots;
  for (int i=0;i<3;++i) {
    plots[i].reset(new TH2D(Form("lpp_theta_p_%d",i),
      Form("K^{-}p #rightarrow #Lambda#pi^{+}#pi^{-}: %s, beam %.0f MeV/c;#theta_{%s}^{lab} [degree];p_{%s}^{lab} [GeV/c]",particles[i],beamMomentumMeV,particles[i],particles[i]),
      kThetaBins,0.,180.,kMomentumBins,0.,i==0 ? kLambdaMomentumMax:kPionMomentumMax));
    plots[i]->SetDirectory(nullptr);plots[i]->GetZaxis()->SetTitle("");
  }
  double decayMasses[]={kProtonMass,kPionMass};
  for (int event=0;event<events;++event) {
    Unweight(production,random);
    const auto lambda=*production.GetDecay(0);
    const auto pip=*production.GetDecay(1),pim=*production.GetDecay(2);
    const auto closure=lambda+pip+pim-initial;
    if (std::abs(closure.E())>1.e-9 || closure.P()>1.e-9)
      throw std::runtime_error("Production four-momentum closure failed");
    TGenPhaseSpace decay;auto mother=lambda;
    if (!decay.SetDecay(mother,2,decayMasses)) throw std::runtime_error("Lambda decay failed");
    Unweight(decay,random);
    const auto reconstructed=Measure(*decay.GetDecay(0),kProtonMass,kProtonSigma,applyResolution,random)
      +Measure(*decay.GetDecay(1),kPionMass,kPionSigma,applyResolution,random);
    const std::array<TLorentzVector,3> used{{reconstructed,
      Measure(pip,kPionMass,kPionSigma,applyResolution,random),
      Measure(pim,kPionMass,kPionSigma,applyResolution,random)}};
    if (!applyResolution && (used[0]-lambda).P()>1.e-9)
      throw std::runtime_error("Truth Lambda reconstruction closure failed");
    for (int i=0;i<3;++i) plots[i]->Fill(used[i].Theta()*180./std::acos(-1.),used[i].P());
  }
  // Extra Lambda eta pages: truth lab-angle selection, with no detector resolution.
  TH2D etaMap("lpp_eta_truth","K^{-}p #rightarrow #Lambda#eta: lab kinematics;#theta_{#Lambda}^{lab} [degree];p_{#Lambda}^{lab} [GeV/c]",
    300,0.,kEtaThetaMaxDeg,400,0.,kLambdaMomentumMax);
  etaMap.SetDirectory(nullptr);
  std::vector<std::unique_ptr<TH1D>> etaSlices;
  for(std::size_t i=0;i<kEtaLabAnglesDeg.size();++i) {
    etaSlices.emplace_back(new TH1D(Form("lpp_eta_slice_%zu",i),
      Form("#Lambda#eta: lab #theta = %.1f #pm %.2f degree;p_{#Lambda}^{lab} [GeV/c];Counts",kEtaLabAnglesDeg[i],kEtaAngleHalfWidthDeg),
      600,0.,kLambdaMomentumMax));
    etaSlices.back()->SetDirectory(nullptr);
  }
  TGraph lowerBranch,upperBranch;
  for(int i=0;i<=3000;++i) {
    const double theta=kEtaThetaMaxDeg*i/3000.;const auto roots=EtaLabSolutions(initial,theta);
    if(roots.size()==2) {
      lowerBranch.SetPoint(lowerBranch.GetN(),theta,roots[0]);upperBranch.SetPoint(upperBranch.GetN(),theta,roots[1]);
    }
  }
  double etaMasses[]={kLambdaMass,kEtaMass};TGenPhaseSpace etaProduction;
  const bool etaAllowed=etaProduction.SetDecay(initial,2,etaMasses);
  if(etaAllowed) {
    TRandom3 etaRandom(kSeed+1);gRandom=&etaRandom;
    for(int event=0;event<events;++event) {
      Unweight(etaProduction,etaRandom);const auto& lambda=*etaProduction.GetDecay(0);
      const double theta=lambda.Theta()*180./std::acos(-1.);
      etaMap.Fill(theta,lambda.P());
      for(std::size_t i=0;i<etaSlices.size();++i)
        if(std::abs(theta-kEtaLabAnglesDeg[i])<kEtaAngleHalfWidthDeg)etaSlices[i]->Fill(lambda.P());
    }
    gRandom=&random;
  } else Warning("phase_lambda_pip_pim","Lambda eta below threshold: extra pages skipped");
  TString logPdf(outputPdf);if(logPdf.EndsWith(".pdf"))logPdf.Remove(logPdf.Length()-4);logPdf+="_log.pdf";
  TCanvas canvas("c_lambda_pip_pim","Lab momentum vs angle",1100,1000);
  for (bool logarithmic:{false,true}) {
    const TString pdf=logarithmic ? logPdf:TString(outputPdf);canvas.Print(pdf+"[");
    for (int i=0;i<3;++i) {
      canvas.Clear();canvas.SetLeftMargin(.15);canvas.SetRightMargin(.16);canvas.SetBottomMargin(.14);
      canvas.SetLogy(false);canvas.SetLogz(logarithmic);plots[i]->SetMinimum(logarithmic ? .5:0.);
      plots[i]->Draw("colz");canvas.Print(pdf);
      if (!logarithmic) Info("phase_lambda_pip_pim","%s: entries %.0f, plotted integral %.0f",particles[i],plots[i]->GetEntries(),plots[i]->Integral());
    }
    if(etaAllowed) {
      canvas.Clear();canvas.SetLogy(false);canvas.SetLogz(logarithmic);
      etaMap.SetMinimum(logarithmic ? .5:0.);etaMap.Draw("colz");
      for(auto* graph:{&lowerBranch,&upperBranch}) {graph->SetLineColor(kRed);graph->SetLineWidth(1);graph->Draw("L same");}
      canvas.Print(pdf);
      for(std::size_t i=0;i<etaSlices.size();++i) {
        canvas.Clear();canvas.SetLogz(false);canvas.SetLogy(logarithmic);
        auto& h=*etaSlices[i];h.SetMinimum(logarithmic ? .5:0.);
        h.SetMaximum(std::max(1.,h.GetMaximum()*(logarithmic ? 2.:1.15)));h.Draw("hist");
        std::vector<TLine> lines;
        const auto roots=EtaLabSolutions(initial,kEtaLabAnglesDeg[i]);
        lines.reserve(roots.size());
        for(double p:roots) {
          lines.emplace_back(p,logarithmic ? .5:0.,p,h.GetMaximum());
          lines.back().SetLineColor(kRed);lines.back().SetLineWidth(1);lines.back().Draw("same");
          if(!logarithmic)Info("eta_lab_angle","theta %.2f degree: exact p_Lambda %.6f GeV/c",kEtaLabAnglesDeg[i],p);
        }
        canvas.Print(pdf);
      }
    }
    canvas.Print(pdf+"]");
  }
}
