// Run: root -l -b -q lambda_max_scattering_angle.cc
// Exact lab kinematic boundaries for direct K- p -> Lambda X production.
// Two-body X=eta/pi0: the theta(pLambda) curve is the allowed locus itself.
// Three-body X=pi+pi-: minimum recoil mass 2 mPi gives the upper-angle envelope.
// No resolution, acceptance or resonance model. All boundaries include endpoints.
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
namespace lambda_angle_limit {
// ---------------- Parameters ----------------
constexpr double kBeamMomentumMeV = 735.;
constexpr double kBeamScanMinMeV = 710., kBeamScanMaxMeV = 770.;
constexpr double kLambdaMomentumMaxGeV = 1.;
constexpr double kKaonMass=.493677,kProtonMass=.938272,kLambdaMass=1.115683;
struct Reaction {bool enabled;const char* name;double recoilMinMass;Color_t color;};
const std::vector<Reaction> kReactions={
  {true,"#Lambda#eta",.547862,kRed+1},
  {true,"#Lambda#pi^{0}",.134977,kBlue+1},
  {true,"#Lambda#pi^{+}#pi^{-}",2.*.139570,kGreen+2}
};
constexpr int kPoints=3000;
// --------------------------------------------
struct Kinematics {double beam,energy,sqrtS,q,lambdaEnergyCm;bool allowed;};
Kinematics Kin(double beam,double recoil) {
  Kinematics k{};k.beam=beam;k.energy=std::hypot(beam,kKaonMass)+kProtonMass;
  k.sqrtS=std::sqrt(k.energy*k.energy-beam*beam);
  k.allowed=k.sqrtS>=kLambdaMass+recoil;
  if(!k.allowed)return k;
  const double s=k.sqrtS*k.sqrtS;
  k.q=std::sqrt(std::max(0.,(s-std::pow(kLambdaMass+recoil,2))*(s-std::pow(kLambdaMass-recoil,2))))/(2.*k.sqrtS);
  k.lambdaEnergyCm=(s+kLambdaMass*kLambdaMass-recoil*recoil)/(2.*k.sqrtS);return k;
}
double MaxAngle(const Kinematics& k) {
  if(!k.allowed)return -1.;
  if(k.beam==0.)return k.q>0. ? 180.:0.;
  const double beta=k.beam/k.energy;
  if(k.q>beta*k.lambdaEnergyCm+1.e-12)return 180.;
  if(std::abs(k.q-beta*k.lambdaEnergyCm)<1.e-12)return 90.;
  return std::asin(std::clamp(k.q*k.sqrtS/(k.beam*kLambdaMass),0.,1.))*180./std::acos(-1.);
}
}
void lambda_max_scattering_angle(double beamMomentumMeV=lambda_angle_limit::kBeamMomentumMeV,
    const char* outputPdf="lambda_max_scattering_angle.pdf") {
  using namespace lambda_angle_limit;
  if(!std::isfinite(beamMomentumMeV)||beamMomentumMeV<=0.||kBeamScanMaxMeV<=kBeamScanMinMeV)
    throw std::runtime_error("Invalid beam momentum or scan range");
  gROOT->SetBatch(true);gStyle->SetOptStat(0);
  std::vector<std::unique_ptr<TGraph>> fixed,scan;std::vector<const Reaction*> active;
  for(const auto& reaction:kReactions) {
    if(!reaction.enabled)continue;active.push_back(&reaction);
    fixed.emplace_back(new TGraph);scan.emplace_back(new TGraph);
    const auto k=Kin(beamMomentumMeV/1000.,reaction.recoilMinMass);
    if(k.allowed) {
      // Generate the recoil-minimum two-body boundary parametrically in CM.
      // For the three-body reaction this bounds all larger dipion masses.
      const double beta=k.beam/k.energy,gamma=k.energy/k.sqrtS;
      for(int i=0;i<=kPoints;++i) {
        const double u=-1.+2.*i/kPoints;
        const double transverse=k.q*std::sqrt(std::max(0.,1.-u*u));
        const double longitudinal=gamma*(k.q*u+beta*k.lambdaEnergyCm);
        const double momentum=std::hypot(transverse,longitudinal);
        if(momentum<1.e-12)continue;
        const double theta=std::atan2(transverse,longitudinal)*180./std::acos(-1.);
        // Independent invariant-mass check of every boundary point.
        const double missing2=std::pow(k.energy-std::hypot(momentum,kLambdaMass),2)-
          (transverse*transverse+std::pow(k.beam-longitudinal,2));
        if(std::abs(missing2-reaction.recoilMinMass*reaction.recoilMinMass)>1.e-9)
          throw std::runtime_error("Boundary recoil mass closure failed");
        fixed.back()->SetPoint(fixed.back()->GetN(),momentum,theta);
      }
    }
    for(int i=0;i<=kPoints;++i) {
      const double beam=kBeamScanMinMeV+(kBeamScanMaxMeV-kBeamScanMinMeV)*i/kPoints;
      const double angle=MaxAngle(Kin(beam/1000.,reaction.recoilMinMass));
      if(angle>=0.)scan.back()->SetPoint(scan.back()->GetN(),beam,angle);
    }
    for(auto* graph:{fixed.back().get(),scan.back().get()}){graph->SetLineColor(reaction.color);graph->SetLineWidth(2);}
    Info("lambda_angle_limit","%s at %.1f MeV/c: maximum lab angle %.6f degree",reaction.name,beamMomentumMeV,MaxAngle(k));
  }
  if(active.empty())throw std::runtime_error("No enabled reaction");
  TCanvas canvas("c_lambda_angle_limit","Maximum Lambda lab scattering angle",1100,1000);
  canvas.SetLeftMargin(.15);canvas.SetBottomMargin(.14);canvas.SetRightMargin(.08);
  canvas.Print(TString(outputPdf)+"[");
  auto legend=[&](const std::vector<std::unique_ptr<TGraph>>& graphs){
    auto l=std::make_unique<TLegend>(.50,.69,.89,.88);l->SetFillStyle(0);l->SetBorderSize(0);l->SetTextSize(.030);
    for(std::size_t i=0;i<active.size();++i)l->AddEntry(graphs[i].get(),active[i]->name,"l");return l;
  };
  double ymax=0.;for(auto* r:active)ymax=std::max(ymax,MaxAngle(Kin(beamMomentumMeV/1000.,r->recoilMinMass)));
  auto* frame=canvas.DrawFrame(0.,0.,kLambdaMomentumMaxGeV,std::max(10.,std::min(180.,ymax*1.15)));
  frame->SetTitle(Form("Lab angular boundary, beam %.0f MeV/c;p_{#Lambda}^{lab} [GeV/c];#theta_{#Lambda,max}^{lab} [degree]",beamMomentumMeV));
  for(auto& graph:fixed)if(graph->GetN())graph->Draw("L same");
  auto l1=legend(fixed);l1->Draw();canvas.Print(outputPdf);canvas.Clear();
  ymax=0.;for(auto& graph:scan)for(int i=0;i<graph->GetN();++i)ymax=std::max(ymax,graph->GetY()[i]);
  frame=canvas.DrawFrame(kBeamScanMinMeV,0.,kBeamScanMaxMeV,std::max(10.,std::min(180.,ymax*1.15)));
  frame->SetTitle("Maximum Lambda lab angle;p_{K^{-}}^{lab} [MeV/c];#theta_{#Lambda,max}^{lab} [degree]");
  for(auto& graph:scan)if(graph->GetN())graph->Draw("L same");
  auto l2=legend(scan);l2->Draw();canvas.Print(outputPdf);canvas.Clear();
  // Dedicated Lambda-eta zoom, with exact maximum-angle point.
  for(std::size_t i=0;i<active.size();++i) {
    if(std::abs(active[i]->recoilMinMass-.547862)>1.e-9||fixed[i]->GetN()==0)continue;
    const auto k=Kin(beamMomentumMeV/1000.,active[i]->recoilMinMass);
    double pmin=1.e9,pmax=0.;for(int j=0;j<fixed[i]->GetN();++j){pmin=std::min(pmin,fixed[i]->GetX()[j]);pmax=std::max(pmax,fixed[i]->GetX()[j]);}
    frame=canvas.DrawFrame(std::max(0.,pmin-.02),0.,pmax+.02,std::max(1.,MaxAngle(k)*1.2));
    frame->SetTitle(Form("#Lambda#eta lab angular boundary, beam %.0f MeV/c;p_{#Lambda}^{lab} [GeV/c];#theta_{#Lambda}^{lab} [degree]",beamMomentumMeV));
    fixed[i]->Draw("L same");canvas.Print(outputPdf);canvas.Clear();
  }
  canvas.Print(TString(outputPdf)+"]");
}
