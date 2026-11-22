// root -l -b -q 'hyperon_pipi_kinematic_comparison.cc'
// Pure unweighted phase space, correct primary pion assignment; no acceptance.
// Reconstruct Lambda from p pi- without mass constraint. Sigma0 -> Lambda gamma.
// Resolution cases: ideal and independent sigma(p)/p on ALL FOUR charged tracks.
// Exact beam and track directions. Default mixture is specific to 735 MeV/c.
#include "lambda-pipi-phase.cc"
#include <TFile.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TNamed.h>
#include <fstream>
#include <vector>
namespace hyperon_compare {
struct Var{const char* name;const char* title;double lo,hi;};
const std::vector<Var> vars={
 {"mm_pipi","MM(#pi#pi) [GeV/c^{2}]",.9,1.4},
 {"p_lambda","p_{#Lambda}^{lab} [GeV/c]",0,1.1},
 {"theta_lambda","#theta_{#Lambda}^{lab} [degree]",0,90},
 {"m_pipi","M(#pi#pi) [GeV/c^{2}]",.27,.65},
 {"p_pipi","p_{#pi#pi}^{lab} [GeV/c]",0,.95},
 {"theta_pipi","#theta_{#pi#pi}^{lab} [degree]",0,180},
 {"opening","#angle(#Lambda,#pi#pi)^{lab} [degree]",0,180},
 {"missing_energy","E_{miss}^{lab} [GeV]",-.3,.4},
 {"missing_p","|p_{miss}^{lab}| [GeV/c]",0,.45},
 {"cm_deficit","#sqrt{s}-E_{visible}^{CM} [GeV]",-.25,.3},
 {"visible_mass","M(#Lambda#pi#pi) [GeV/c^{2}]",1.25,1.95},
 {"missing_m2","M_{miss}^{2}(#Lambda#pi#pi) [GeV^{2}/c^{4}]",-.12,.12},
 {"p_lambda_cm","p_{#Lambda}^{CM} [GeV/c]",0,.5},
 {"m_lpip","M(#Lambda#pi^{+}) [GeV/c^{2}]",1.2,1.7},
 {"m_lpim","M(#Lambda#pi^{-}) [GeV/c^{2}]",1.2,1.7},
 {"pt_balance","|p_{T,visible}| [GeV/c]",0,.4},
 {"lambda_mass","M(p#pi^{-}) [GeV/c^{2}]",1.04,1.2},
 {"p_pip","p_{#pi^{+}}^{lab} [GeV/c]",0,.7},
 {"p_pim","p_{#pi^{-}}^{lab} [GeV/c]",0,.7},
 {"pion_opening","#angle(#pi^{+},#pi^{-})^{lab} [degree]",0,180}
};
}
void hyperon_pipi_kinematic_comparison(double beamMeV=735.,int events=100000,
 double resolution=.10,double sigmaFraction=.17075){
 using namespace lambda_pipi_phase_config;using namespace hyperon_compare;
 if(events<=0||!std::isfinite(beamMeV)||beamMeV<=0||!std::isfinite(resolution)||resolution<0||
 !std::isfinite(sigmaFraction)||sigmaFraction<0||sigmaFraction>1)throw std::runtime_error("Invalid arguments");
 gROOT->SetBatch(true);gStyle->SetOptStat(0);gStyle->SetPalette(kViridis);gStyle->SetTitleFontSize(.035);
 gStyle->SetPaperSize(30.,14.);
 TRandom3 random(20261117);struct Guard{TRandom* old;~Guard(){gRandom=old;}}guard{gRandom};gRandom=&random;
 TLorentzVector initial(0,0,beamMeV/1000.,std::hypot(beamMeV/1000.,kKaonMass)+kProtonMass);
 const double mass[]={kLambdaMass,1.192642},fraction[]={1-sigmaFraction,sigmaFraction};
 const char* channel[]={"Lambda pipi","Sigma0 pipi"};
 const int nv=vars.size();
 // Exhaustive pairwise atlas of the twelve core observables, plus Dalitz and tracks.
 std::vector<std::pair<int,int>> pairs;
 for(int a=0;a<12;++a)for(int b=a+1;b<12;++b)pairs.emplace_back(a,b);
 for(auto p:std::vector<std::pair<int,int>>{{12,3},{12,9},{13,14},{13,3},{14,3},{15,7},{15,9},{16,0},{16,9},{17,18},{19,3},{19,0}})pairs.push_back(p);
 struct Plots{std::vector<std::unique_ptr<TH1D>> h;std::vector<std::unique_ptr<TH2D>> h2;};
 Plots plots[2][2];
 for(int r=0;r<2;++r)for(int c=0;c<2;++c){
  for(int v=0;v<nv;++v){const auto& x=vars[v];
   auto h=std::make_unique<TH1D>(Form("r%d_c%d_%s",r,c,x.name),Form(";%s;Probability / bin",x.title),120,x.lo,x.hi);
   h->SetDirectory(nullptr);h->Sumw2();plots[r][c].h.push_back(std::move(h));}
  for(size_t j=0;j<pairs.size();++j){auto a=pairs[j].first,b=pairs[j].second;
   auto h=std::make_unique<TH2D>(Form("r%d_c%d_%s_vs_%s",r,c,vars[b].name,vars[a].name),
    Form(";%s;%s",vars[a].title,vars[b].title),65,vars[a].lo,vars[a].hi,65,vars[b].lo,vars[b].hi);
   h->SetDirectory(nullptr);h->Sumw2();plots[r][c].h2.push_back(std::move(h));}
 }
 std::ofstream report("hyperon_pipi_kinematic_comparison.txt");
 report<<"Pure phase space; correct pion assignment; no acceptance or resonance model.\n"
 <<"All FOUR charged tracks smeared independently in magnitude; directions/beam exact. No Lambda mass constraint.\n";
 double selected[2][2][3]={};
 auto smear=[&](const TLorentzVector& p,double m,double z,double res){TLorentzVector q;q.SetVectM(p.Vect()*(1+z*res),m);return q;};
 for(int c=0;c<2;++c){
  double pm[]={mass[c],kPionMass,kPionMass};TGenPhaseSpace production;
  if(!production.SetDecay(initial,3,pm))throw std::runtime_error("Production threshold");
  for(int e=0;e<events;++e){
   Generate(production,random);auto y=*production.GetDecay(0),pip=*production.GetDecay(1),pim=*production.GetDecay(2),lambda=y;
   TLorentzVector gamma;
   if(c){double dm[]={kLambdaMass,0.};TGenPhaseSpace d;if(!d.SetDecay(y,2,dm))throw std::runtime_error("Sigma decay");Generate(d,random);lambda=*d.GetDecay(0);gamma=*d.GetDecay(1);}
   double dm[]={kProtonMass,kPionMass};TGenPhaseSpace d;if(!d.SetDecay(lambda,2,dm))throw std::runtime_error("Lambda decay");Generate(d,random);
   auto proton=*d.GetDecay(0),daughter=*d.GetDecay(1);
   auto closure=proton+daughter+pip+pim+gamma-initial;
   if(std::abs(closure.E())>1e-9||closure.P()>1e-9||std::abs((initial-pip-pim).M()-mass[c])>1e-9)throw std::runtime_error("Closure failed");
   double z[4];for(auto& v:z)do{v=random.Gaus();}while(1+resolution*v<0);
   for(int r=0;r<2;++r){double res=r?resolution:0.;
    auto a=smear(pip,kPionMass,z[0],res),b=smear(pim,kPionMass,z[1],res);
    auto l=smear(proton,kProtonMass,z[2],res)+smear(daughter,kPionMass,z[3],res);
    auto pp=a+b,visible=l+pp,missing=initial-visible,cm=visible,lcm=l;
    cm.Boost(-initial.BoostVector());lcm.Boost(-initial.BoostVector());
    double val[]={ (initial-pp).M(),l.P(),l.Theta()*TMath::RadToDeg(),pp.M(),pp.P(),pp.Theta()*TMath::RadToDeg(),
     l.Vect().Angle(pp.Vect())*TMath::RadToDeg(),missing.E(),missing.P(),initial.M()-cm.E(),visible.M(),missing.M2(),
     lcm.P(),(l+a).M(),(l+b).M(),visible.Pt(),l.M(),a.P(),b.P(),a.Vect().Angle(b.Vect())*TMath::RadToDeg()};
    for(int v=0;v<nv;++v)plots[r][c].h[v]->Fill(val[v],1./events);
    for(size_t j=0;j<pairs.size();++j)plots[r][c].h2[j]->Fill(val[pairs[j].first],val[pairs[j].second],1./events);
    // Fixed illustrative Sigma-enrichment cuts, not optimized on this sample.
    if(val[0]>(mass[0]+mass[1])/2)selected[r][c][0]++;
    if(val[9]>.04)selected[r][c][1]++;
    if(val[7]>.04)selected[r][c][2]++;
   }
  }
 }
 TFile output("hyperon_pipi_kinematic_comparison.root","RECREATE");
 if(output.IsZombie())throw std::runtime_error("Cannot write ROOT output");
 TNamed metadata("assumptions",Form("beam %.1f MeV/c; sigma fraction %.5f; all-track sigma_p/p %.4f; exact directions/beam; pure phase space; correct pion assignment",beamMeV,sigmaFraction,resolution));metadata.Write();
 TCanvas canvas("c_hyperon_compare","Kinematic comparison",1500,650);
 const char* pdf="hyperon_pipi_kinematic_comparison.pdf";canvas.Print(TString(pdf)+"[");
 canvas.cd();TLatex intro;intro.SetNDC();intro.SetTextSize(.033);
 const std::vector<TString> notes={
  "Lambda pipi vs Sigma0 pipi: kinematic atlas",
  Form("Beam %.0f MeV/c; %d events/channel; Sigma0 production fraction %.3f%%",beamMeV,events,100*sigmaFraction),
  "Sigma0 -> Lambda gamma; Lambda reconstructed from p pi-; photon omitted",
  Form("Smeared case: independent %.1f%% momentum errors on ALL FOUR charged tracks",100*resolution),
  "Exact beam/directions; correct pion assignment; no Lambda mass constraint",
  "Pure phase space; no resonance model, detector acceptance, or selection cuts",
  "Smeared case first; ideal case second. See .txt for page index and cut results.",
  "1D: unit-area shapes and production-weighted shapes; 2D: Lambda / Sigma0 / mixture",
  "Core variables: MM(pipi), Lambda and pipi momenta/angles, M(pipi), opening angle,",
  "missing energy/momentum, CM energy deficit, visible mass, missing mass squared.",
  "78 2D pairs/case: all 66 pairs among 12 core variables + 12 further correlations.",
  "CM deficit = sqrt(s) - E_CM(Lambda + pi+ + pi-). No mass constraint used.",
  "Missing mass squared of the full visible system is zero for BOTH ideal reactions.",
  "The missing energy and momentum distinguish no missing particle from a photon."
 };
 for(size_t i=0;i<notes.size();++i)intro.DrawLatex(.05,.92-.06*i,notes[i]);
 canvas.Print(pdf);canvas.Clear();int page=1;
 auto header=[&](const char* title,int r){canvas.cd();TLatex text;text.SetNDC();text.SetTextAlign(13);text.SetTextSize(.025);text.DrawLatex(.08,.985,Form("%s | p_{K} = %.0f MeV/c | all-track #sigma_{p}/p = %.1f%%",title,beamMeV,100.*(r?resolution:0.)));};
 for(int r:{1,0}){
  for(int v=0;v<nv;++v){canvas.Clear();canvas.Divide(2,1,.005,.07);
   std::vector<std::unique_ptr<TH1D>> displayed;
   for(int mode=0;mode<2;++mode){canvas.cd(mode+1);gPad->SetLeftMargin(.14);gPad->SetBottomMargin(.16);gPad->SetTopMargin(.10);gPad->SetLogy(mode);
    double peak=0;for(int c=0;c<2;++c)peak=std::max(peak,plots[r][c].h[v]->GetMaximum()*(mode?fraction[c]:1.));
    std::unique_ptr<TH1D> shown[2];for(int c=0;c<2;++c){shown[c].reset(static_cast<TH1D*>(plots[r][c].h[v]->Clone(Form("display_%d",c))));shown[c]->SetDirectory(nullptr);if(mode)shown[c]->Scale(fraction[c]);shown[c]->SetLineColor(c?kRed+1:kBlue+1);shown[c]->SetLineWidth(2);shown[c]->SetTitle(mode?"Cross-section fractions (log)":"Unit-area shapes");shown[c]->SetMinimum(mode?1e-6:0);shown[c]->SetMaximum(peak*(mode?6:1.3));shown[c]->Draw(c?"hist same":"hist");}
    TLegend leg(.59,.72,.94,.89);leg.SetTextSize(.032);leg.SetBorderSize(0);leg.SetFillStyle(0);for(int c=0;c<2;++c)leg.AddEntry(shown[c].get(),channel[c],"l");leg.Draw();
    for(auto& h:shown)displayed.push_back(std::move(h));leg.DrawClone();
   }header(vars[v].name,r);canvas.Print(pdf);
   report<<"page "<<++page<<" resolution="<<(r?resolution:0)<<" 1D "<<vars[v].name<<"\n";
  }
  for(size_t j=0;j<pairs.size();++j){canvas.Clear();canvas.Divide(3,1,.005,.07);
   std::unique_ptr<TH2D> mix(static_cast<TH2D*>(plots[r][0].h2[j]->Clone("display_mix")));mix->SetDirectory(nullptr);mix->Scale(fraction[0]);mix->Add(plots[r][1].h2[j].get(),fraction[1]);
   double maximum=std::max(plots[r][0].h2[j]->GetMaximum(),plots[r][1].h2[j]->GetMaximum());
   for(int c=0;c<3;++c){canvas.cd(c+1);gPad->SetLeftMargin(.17);gPad->SetRightMargin(.19);gPad->SetBottomMargin(.18);gPad->SetTopMargin(.12);
    auto* h=c<2?plots[r][c].h2[j].get():mix.get();h->SetTitle(c<2?channel[c]:Form("Mixture: %.1f%% / %.1f%%",100*fraction[0],100*fraction[1]));h->SetMinimum(0);h->SetMaximum(maximum);h->GetXaxis()->SetTitleSize(.045);h->GetYaxis()->SetTitleSize(.045);h->Draw("colz");}
   header(Form("%s vs %s",vars[pairs[j].second].name,vars[pairs[j].first].name),r);canvas.Print(pdf);canvas.Clear();
   report<<"page "<<++page<<" resolution="<<(r?resolution:0)<<" 2D "<<vars[pairs[j].first].name<<" vs "<<vars[pairs[j].second].name<<"\n";
   output.cd();mix->Write(Form("r%d_mix_pair%zu",r,j));
  }
  const char* cuts[]={"MM(pipi) > midpoint masses","CM missing energy > 40 MeV","Lab missing energy > 40 MeV"};
  for(int k=0;k<3;++k){double el=selected[r][0][k]/events,es=selected[r][1][k]/events;
   double denom=fraction[0]*el+fraction[1]*es;
   report<<"resolution="<<(r?resolution:0)<<" cut="<<cuts[k]<<" Lambda efficiency="<<el<<" Sigma efficiency="<<es<<" Sigma purity="<<(denom?fraction[1]*es/denom:0)<<"\n";
  }
  output.cd();for(int c=0;c<2;++c){for(auto& h:plots[r][c].h)h->Write();for(auto& h:plots[r][c].h2)h->Write();}
 }
 canvas.Print(TString(pdf)+"]");output.Close();report.close();
}
