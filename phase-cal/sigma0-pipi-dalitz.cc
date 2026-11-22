// Run: root -l -b -q sigma0-pipi-dalitz.cc
// Phase-space events reweighted by an incoherent Breit-Wigner toy model.
// Sigma0 -> Lambda gamma, Lambda -> p pi-. No detector resolution or cuts.
// Toy fractions are configurable assumptions, not measured resonance yields.
// Normalize the summed model to the measured production cross section.
#include "lambda-pipi-phase.cc"
#include <TGraph.h>
#include <TFile.h>
#include <TLine.h>
#include <TLegend.h>
namespace sigma0_pipi {
constexpr double kBeamMomentumMeV=735.;
constexpr int kEvents=100000;
constexpr double kSigmaMass=1.192642;
constexpr double kSigmaStarPlus=1.38283,kSigmaStarMinus=1.3872;
constexpr bool kShowResonanceGuides=true;
constexpr double kSigmaStarPlusWidth=.0362,kSigmaStarMinusWidth=.0394;
// Illustrative event fractions AFTER integration over the allowed phase space.
constexpr double kNonresonantFraction=.20,kPlusFraction=.40,kMinusFraction=.40;
constexpr bool kUseCrossSection=true;
const char* const kCrossSectionFile="../lambda-yield/cs-info/total_cross_section.root";
// Constant-width relativistic |BW(s)|^2 with unit height at s=m0^2.
// No spin/helicity amplitudes, P-wave factors, interference or detector acceptance.
double BW(double invariantSquared,double mass,double width) {
  const double scale=mass*width;
  return scale*scale/(std::pow(invariantSquared-mass*mass,2)+scale*scale);
}
double CrossSection(double beamMeV) {
  if(!kUseCrossSection)return 1.;
  TFile input(kCrossSectionFile,"READ");
  if(input.IsZombie())throw std::runtime_error("Cannot open cross-section file");
  auto* graph=dynamic_cast<TGraph*>(input.Get("g_Sigma0_pip_pim_pK"));
  if(!graph||graph->GetN()<2)throw std::runtime_error("Missing Sigma0 pi+ pi- cross section");
  double lo=graph->GetX()[0],hi=lo;
  for(int i=1;i<graph->GetN();++i){lo=std::min(lo,graph->GetX()[i]);hi=std::max(hi,graph->GetX()[i]);}
  if(beamMeV<lo||beamMeV>hi)throw std::runtime_error("Beam outside cross-section data range; no extrapolation");
  const double sigma=graph->Eval(beamMeV);
  if(!std::isfinite(sigma)||sigma<0.)throw std::runtime_error("Invalid cross section");return sigma;
}
void SetFloor(TH1* plot,bool logarithmic) {
  double minimum=1.e100;
  for(int bin=0;bin<plot->GetNcells();++bin)if(plot->GetBinContent(bin)>0.)minimum=std::min(minimum,plot->GetBinContent(bin));
  plot->SetMinimum(logarithmic ? (minimum==1.e100 ? 1.e-12:minimum*.5):0.);
}
// PDG: https://pdg.lbl.gov/2025/listings/rpp2025-list-sigma-1385.pdf
TGraph Boundary(double W,double m1,double m2,double m3) {
  TGraph graph;std::vector<std::array<double,3>> points;
  for(int i=0;i<=1000;++i){
    double x=std::pow(m1+m2,2)+(std::pow(W-m3,2)-std::pow(m1+m2,2))*i/1000.;
    double E1=(x+m1*m1-m2*m2)/(2.*std::sqrt(x)),E3=(W*W-x-m3*m3)/(2.*std::sqrt(x));
    double delta=2.*std::sqrt(std::max(0.,E1*E1-m1*m1))*std::sqrt(std::max(0.,E3*E3-m3*m3));
    double center=m1*m1+m3*m3+2.*E1*E3;points.push_back({x,center-delta,center+delta});
  }
  for(auto&p:points)graph.SetPoint(graph.GetN(),p[0],p[1]);
  for(auto it=points.rbegin();it!=points.rend();++it)graph.SetPoint(graph.GetN(),(*it)[0],(*it)[2]);
  graph.SetPoint(graph.GetN(),points[0][0],points[0][1]);graph.SetLineWidth(2);return graph;
}
}
void sigma0_pipi_dalitz(double beamMomentumMeV=sigma0_pipi::kBeamMomentumMeV,
    int events=sigma0_pipi::kEvents,const char* outputPdf="sigma0-pipi-dalitz.pdf") {
  using namespace lambda_pipi_phase_config;using namespace sigma0_pipi;
  if(!std::isfinite(beamMomentumMeV)||beamMomentumMeV<=0.||events<=0)throw std::runtime_error("Invalid beam or event count");
  const double fractionSum=kNonresonantFraction+kPlusFraction+kMinusFraction;
  if(kNonresonantFraction<0.||kPlusFraction<0.||kMinusFraction<0.||std::abs(fractionSum-1.)>1.e-10)
    throw std::runtime_error("Nonnegative toy fractions must sum to one");
  const double crossSection=CrossSection(beamMomentumMeV);
  gROOT->SetBatch(true);gStyle->SetOptStat(0);
  TRandom3 random(20261102);struct Guard{TRandom* old;~Guard(){gRandom=old;}}guard{gRandom};gRandom=&random;
  TLorentzVector initial(0.,0.,beamMomentumMeV/1000.,std::hypot(beamMomentumMeV/1000.,kKaonMass)+kProtonMass);
  double masses[]={kSigmaMass,kPionMass,kPionMass},sigmaDaughters[]={kLambdaMass,0.},lambdaDaughters[]={kProtonMass,kPionMass};
  TGenPhaseSpace production;if(!production.SetDecay(initial,3,masses))throw std::runtime_error("Sigma0 pi+ pi- below threshold");
  const double W=initial.M(),smin=std::pow(kSigmaMass+kPionMass,2),smax=std::pow(W-kPionMass,2);
  const double ppmin=4.*kPionMass*kPionMass,ppmax=std::pow(W-kSigmaMass,2);
  TH2D dalitz("spp_dalitz","K^{-}p #rightarrow #Sigma^{0}#pi^{+}#pi^{-}: incoherent BW toy;M^{2}(#Sigma^{0}#pi^{+}) [GeV^{2}/c^{4}];M^{2}(#Sigma^{0}#pi^{-}) [GeV^{2}/c^{4}]",220,smin-.015,smax+.015,220,smin-.015,smax+.015);
  TH2D alternate("spp_alternate","#Sigma^{0}#pi^{+}#pi^{-}: incoherent BW toy;M^{2}(#pi^{+}#pi^{-}) [GeV^{2}/c^{4}];M^{2}(#Sigma^{0}#pi^{+}) [GeV^{2}/c^{4}]",220,ppmin-.005,ppmax+.005,220,smin-.015,smax+.015);
  TH1D mPlus("spp_mplus","#Sigma^{0}#pi^{+};M(#Sigma^{0}#pi^{+}) [GeV/c^{2}];Counts",240,std::sqrt(smin),std::sqrt(smax));
  TH1D mMinus("spp_mminus","#Sigma^{0}#pi^{-};M(#Sigma^{0}#pi^{-}) [GeV/c^{2}];Counts",240,std::sqrt(smin),std::sqrt(smax));
  TH1D mPipi("spp_mpipi","#pi^{+}#pi^{-};M(#pi^{+}#pi^{-}) [GeV/c^{2}];Counts",240,2.*kPionMass,W-kSigmaMass);
  TH1D photon("spp_photon","Photon from #Sigma^{0};E_{#gamma}^{lab} [GeV];Counts",200,0.,.3);
  const char* names[]={"#Sigma^{0}","#Lambda from #Sigma^{0}","#pi^{+} from reaction","#pi^{-} from reaction","p from #Lambda","#pi^{-} from #Lambda"};
  std::array<std::unique_ptr<TH2D>,6> phases;
  for(int i=0;i<6;++i)phases[i].reset(new TH2D(Form("spp_phase_%d",i),Form("%s;#theta^{lab} [degree];p^{lab} [GeV/c]",names[i]),180,0.,180.,200,0.,1.2));
  auto phaseSpace=std::unique_ptr<TH2D>(static_cast<TH2D*>(dalitz.Clone("spp_pure_phase_space")));
  phaseSpace->SetDirectory(nullptr);
  phaseSpace->SetTitle("#Sigma^{0}#pi^{+}#pi^{-}: pure phase space;M^{2}(#Sigma^{0}#pi^{+}) [GeV^{2}/c^{4}];M^{2}(#Sigma^{0}#pi^{-}) [GeV^{2}/c^{4}]");
  const std::vector<TH1*> allPlots{&dalitz,&alternate,&mPlus,&mMinus,&mPipi,&photon,phaseSpace.get()};
  for(auto* plot:allPlots)plot->Sumw2();
  for(auto& plot:phases)plot->Sumw2();
  struct Event {double a,b,c,photonEnergy;std::array<TLorentzVector,6> particles;};
  std::vector<Event> sample;sample.reserve(events);
  double plusIntegral=0.,minusIntegral=0.;
  for(int event=0;event<events;++event){
    Generate(production,random);auto sigma=*production.GetDecay(0);const auto pip=*production.GetDecay(1),pim=*production.GetDecay(2);
    TGenPhaseSpace sd;if(!sd.SetDecay(sigma,2,sigmaDaughters))throw std::runtime_error("Sigma decay failed");Generate(sd,random);
    auto lambda=*sd.GetDecay(0);const auto gamma=*sd.GetDecay(1);
    TGenPhaseSpace ld;if(!ld.SetDecay(lambda,2,lambdaDaughters))throw std::runtime_error("Lambda decay failed");Generate(ld,random);
    const auto p=*ld.GetDecay(0),pi=*ld.GetDecay(1);
    const auto reconstructed=p+pi+gamma;
    const auto closure=p+pi+gamma+pip+pim-initial;
    if(std::abs(closure.E())>1.e-9||closure.P()>1.e-9||std::abs(reconstructed.M()-kSigmaMass)>1.e-9)throw std::runtime_error("Decay closure failed");
    const double a=(sigma+pip).M2(),b=(sigma+pim).M2(),c=(pip+pim).M2();
    if(std::abs(a+b+c-(initial.M2()+kSigmaMass*kSigmaMass+2.*kPionMass*kPionMass))>1.e-9)throw std::runtime_error("Dalitz identity failed");
    sample.push_back({a,b,c,gamma.E(),{{sigma,lambda,pip,pim,p,pi}}});
    plusIntegral+=BW(a,kSigmaStarPlus,kSigmaStarPlusWidth);
    minusIntegral+=BW(b,kSigmaStarMinus,kSigmaStarMinusWidth);
  }
  if(plusIntegral<=0.||minusIntegral<=0.)throw std::runtime_error("Empty Breit-Wigner normalization");
  double totalWeight=0.,plusWeight=0.,minusWeight=0.,nonresWeight=0.;
  for(const auto& event:sample){
    const double nonres=crossSection*kNonresonantFraction/events;
    const double plus=crossSection*kPlusFraction*BW(event.a,kSigmaStarPlus,kSigmaStarPlusWidth)/plusIntegral;
    const double minus=crossSection*kMinusFraction*BW(event.b,kSigmaStarMinus,kSigmaStarMinusWidth)/minusIntegral;
    const double weight=nonres+plus+minus;
    totalWeight+=weight;plusWeight+=plus;minusWeight+=minus;nonresWeight+=nonres;
    dalitz.Fill(event.a,event.b,weight);alternate.Fill(event.c,event.a,weight);
    mPlus.Fill(std::sqrt(event.a),weight);mMinus.Fill(std::sqrt(event.b),weight);mPipi.Fill(std::sqrt(event.c),weight);photon.Fill(event.photonEnergy,weight);
    phaseSpace->Fill(event.a,event.b,crossSection/events);
    for(int i=0;i<6;++i)phases[i]->Fill(event.particles[i].Theta()*180./std::acos(-1.),event.particles[i].P(),weight);
  }
  if(std::abs(totalWeight-crossSection)>1.e-9*std::max(1.,crossSection))throw std::runtime_error("Cross-section normalization failed");
  Info("sigma0_pipi_model","sigma %.6f mb; nonres %.6f, Sigma*+ %.6f, Sigma*- %.6f mb; no interference or decay BR correction",crossSection,nonresWeight,plusWeight,minusWeight);
  // Weighted integrals are conditional production cross sections, not raw counts.
  for(auto* plot:{&mPlus,&mMinus,&mPipi,&photon})plot->GetYaxis()->SetTitle("Weighted counts");
  auto boundary=Boundary(W,kSigmaMass,kPionMass,kPionMass);
  auto alternativeBoundary=Boundary(W,kPionMass,kPionMass,kSigmaMass);
  TString logPdf(outputPdf);if(logPdf.EndsWith(".pdf"))logPdf.Remove(logPdf.Length()-4);logPdf+="_log.pdf";
  TCanvas canvas("c_sigma0_pipi_dalitz","Sigma0 pi+ pi-",1200,1100);
  for(bool logarithmic:{false,true}){
    const TString pdf=logarithmic?logPdf:TString(outputPdf);canvas.Print(pdf+"[");
    canvas.SetLeftMargin(.16);canvas.SetRightMargin(.17);canvas.SetBottomMargin(.15);canvas.SetLogz(logarithmic);
    SetFloor(&dalitz,logarithmic);dalitz.Draw("colz");boundary.Draw("L same");
    TLine plus(kSigmaStarPlus*kSigmaStarPlus,smin,kSigmaStarPlus*kSigmaStarPlus,smax),minus(smin,kSigmaStarMinus*kSigmaStarMinus,smax,kSigmaStarMinus*kSigmaStarMinus);
    TLegend legend(.20,.72,.61,.87);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.025);
    if(kShowResonanceGuides){for(auto* l:{&plus,&minus}){l->SetLineStyle(2);l->SetLineWidth(1);}plus.SetLineColor(kRed);minus.SetLineColor(kMagenta+1);plus.Draw("same");minus.Draw("same");legend.AddEntry(&plus,"#Sigma(1385)^{+} mass reference","l");legend.AddEntry(&minus,"#Sigma(1385)^{-} mass reference","l");legend.Draw();}
    canvas.Print(pdf);canvas.Clear();
    SetFloor(&alternate,logarithmic);alternate.Draw("colz");alternativeBoundary.Draw("L same");canvas.Print(pdf);canvas.Clear();canvas.SetLogz(false);
    canvas.Divide(2,2);const std::array<TH1D*,4> projections{{&mPlus,&mMinus,&mPipi,&photon}};
    for(int i=0;i<4;++i){canvas.cd(i+1);gPad->SetLeftMargin(.16);gPad->SetBottomMargin(.15);gPad->SetLogy(logarithmic);SetFloor(projections[i],logarithmic);projections[i]->Draw("hist");}
    canvas.Print(pdf);canvas.Clear();canvas.Divide(3,2);
    for(int i=0;i<6;++i){canvas.cd(i+1);gPad->SetLeftMargin(.16);gPad->SetBottomMargin(.15);gPad->SetRightMargin(.17);gPad->SetLogz(logarithmic);SetFloor(phases[i].get(),logarithmic);phases[i]->Draw("colz");}
    canvas.Print(pdf);canvas.Clear();
    canvas.SetLogy(false);canvas.SetLogz(logarithmic);
    SetFloor(phaseSpace.get(),logarithmic);phaseSpace->Draw("colz");boundary.Draw("L same");
    canvas.Print(pdf);canvas.Clear();canvas.Print(pdf+"]");
  }
  Info("sigma0_pipi_dalitz","beam %.1f MeV/c; sqrt(s) %.6f GeV; M(Sigma pi) %.6f..%.6f; M(pi pi) %.6f..%.6f GeV; entries %.0f",beamMomentumMeV,W,std::sqrt(smin),std::sqrt(smax),std::sqrt(ppmin),std::sqrt(ppmax),dalitz.GetEntries());
}
