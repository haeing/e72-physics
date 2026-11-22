// Run: root -l -b -q 'pipi_missing_mass_mixture.cc'
// K- p -> Y pi+ pi-, Y=Lambda or Sigma0, unweighted three-body phase space.
// MM(pi+ pi-) includes the whole hyperon; its decay need not be generated.
// Default production fractions at 735 MeV/c: Lambda=82.925%, Sigma0=17.075%.
// From linear interpolation of g_*_pip_pim_pK in
// ../lambda-yield/cs-info/total_cross_section.root (2.90898 and 0.59898 mb).
// The fixed fraction must be overridden for other beam momenta.
// Independent Gaussian pion momentum-magnitude errors; exact directions/beam.
// No acceptance, resonance model, or branching-ratio corrections.
// Fourth argument is now the MAXIMUM sigma(p)/p in the resolution scan.
#include "lambda-pipi-phase.cc"
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <vector>

void pipi_missing_mass_mixture(double beamMomentumMeV=735., int eventsPerChannel=100000,
    double sigmaFraction=0.17075, double maximumMomentumSigma=0.10,
    const char* outputPdf="pipi_missing_mass_mixture.pdf",
    const char* outputRoot="pipi_missing_mass_mixture.root") {
  using namespace lambda_pipi_phase_config;
  constexpr double sigmaMass=1.192642;
  if(!std::isfinite(beamMomentumMeV)||beamMomentumMeV<=0.||eventsPerChannel<=0||
     !std::isfinite(sigmaFraction)||sigmaFraction<0.||sigmaFraction>1.||
     !std::isfinite(maximumMomentumSigma)||maximumMomentumSigma<0.)
    throw std::runtime_error("Invalid input");
  gROOT->SetBatch(true);gStyle->SetOptStat(0);
  TRandom3 random(20261116);
  struct Guard {TRandom* old;~Guard(){gRandom=old;}} guard{gRandom};gRandom=&random;
  TLorentzVector initial(0.,0.,beamMomentumMeV/1000.,
    std::hypot(beamMomentumMeV/1000.,kKaonMass)+kProtonMass);
  const double masses[]={kLambdaMass,sigmaMass},fractions[]={1.-sigmaFraction,sigmaFraction};
  const char* labels[]={"#Lambda#pi^{+}#pi^{-}","#Sigma^{0}#pi^{+}#pi^{-}"};
  const double scales[]={0.,.1,.2,.3,.5,.7,1.};
  constexpr int n=7;
  std::array<std::array<std::unique_ptr<TH1D>,2>,n> hist;
  std::array<std::unique_ptr<TH1D>,n> sums;
  TGraph widths[2];
  double residualSum[n][2]={},residualSquaredSum[n][2]={};
  for(int r=0;r<n;++r)for(int c=0;c<2;++c){
    hist[r][c].reset(new TH1D(Form("mm_r%d_c%d",r,c),
      Form("#sigma_{p}/p = %.2f%%;MM(#pi^{+}#pi^{-}) [GeV/c^{2}];Weighted entries / 0.5 MeV/c^{2}",100.*maximumMomentumSigma*scales[r]),
      1400,.8,1.5));
    hist[r][c]->SetDirectory(nullptr);hist[r][c]->Sumw2();
    hist[r][c]->SetLineColor(c==0?kBlue+1:kRed+1);hist[r][c]->SetLineWidth(2);
  }
  auto smear=[](const TLorentzVector& p,double z,double resolution){
    TLorentzVector result;
    result.SetVectM(p.Vect()*(1.+resolution*z),kPionMass);return result;
  };
  for(int c=0;c<2;++c){
    double daughters[]={masses[c],kPionMass,kPionMass};TGenPhaseSpace production;
    if(!production.SetDecay(initial,3,daughters))throw std::runtime_error("Below production threshold");
    for(int e=0;e<eventsPerChannel;++e){
      Generate(production,random);
      const auto hyperon=*production.GetDecay(0),pip=*production.GetDecay(1),pim=*production.GetDecay(2);
      const auto closure=hyperon+pip+pim-initial;
      if(std::abs(closure.E())>1.e-9||closure.P()>1.e-9||
         std::abs((initial-pip-pim).M()-masses[c])>1.e-9)
        throw std::runtime_error("Truth closure failed");
      // Common Gaussian draws across scan points isolate the resolution effect.
      double a,b;do{a=random.Gaus();}while(1.+maximumMomentumSigma*a<0.);
      do{b=random.Gaus();}while(1.+maximumMomentumSigma*b<0.);
      for(int r=0;r<n;++r){
        const double resolution=maximumMomentumSigma*scales[r];
        const double mm=(initial-smear(pip,a,resolution)-smear(pim,b,resolution)).M();
        const double residual=mm-masses[c];
        residualSum[r][c]+=residual;residualSquaredSum[r][c]+=residual*residual;
        hist[r][c]->Fill(mm,fractions[c]/eventsPerChannel);
      }
    }
  }
  TFile output(outputRoot,"RECREATE");if(output.IsZombie())throw std::runtime_error("Cannot create ROOT file");
  for(int r=0;r<n;++r){
    sums[r].reset(static_cast<TH1D*>(hist[r][0]->Clone(Form("mm_r%d_sum",r))));
    sums[r]->SetDirectory(nullptr);sums[r]->Add(hist[r][1].get());sums[r]->SetLineColor(kBlack);
    if(std::abs(sums[r]->Integral(0,sums[r]->GetNbinsX()+1)-1.)>1.e-8)throw std::runtime_error("Normalization failed");
    for(int c=0;c<2;++c){
      const double mean=residualSum[r][c]/eventsPerChannel;
      const double rms=1000.*std::sqrt(std::max(0.,residualSquaredSum[r][c]/eventsPerChannel-mean*mean));
      widths[c].SetPoint(r,100.*maximumMomentumSigma*scales[r],rms);
      Info("pipi_mm_resolution","sigma_p/p %.2f%%, channel %d: MM RMS %.4f MeV/c^2, overflow weight %.3g",
        100.*maximumMomentumSigma*scales[r],c,rms,
        hist[r][c]->GetBinContent(0)+hist[r][c]->GetBinContent(hist[r][c]->GetNbinsX()+1));
      hist[r][c]->Write();
    }
    sums[r]->Write();
  }
  TCanvas canvas("c_pipi_mm_mix","Missing mass resolution scan",1100,800);
  canvas.SetLeftMargin(.14);canvas.SetBottomMargin(.14);canvas.SetRightMargin(.05);
  canvas.Print(TString(outputPdf)+"[");
  for(int r=0;r<n;++r)for(bool log:{false,true}){
    canvas.SetLogy(log);
    auto& sum=*sums[r];const double peak=sum.GetBinContent(sum.GetMaximumBin());
    sum.SetMinimum(log?1.e-7:0.);sum.SetMaximum(peak*(log?5.:1.25));sum.Draw("hist");
    for(int c=0;c<2;++c)hist[r][c]->Draw("hist same");
    TLegend legend(.64,.67,.93,.87);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.028);
    legend.AddEntry(&sum,"Mixture","l");
    for(int c=0;c<2;++c)legend.AddEntry(hist[r][c].get(),Form("%s (%.1f%%)",labels[c],100.*fractions[c]),"l");
    legend.Draw();TLatex note;note.SetNDC();note.SetTextSize(.025);
    note.DrawLatex(.16,.93,Form("K^{-}p #rightarrow #pi^{+}#pi^{-}X; p_{K^{-}} = %.0f MeV/c; pion momentum errors only",beamMomentumMeV));
    canvas.Print(outputPdf);canvas.Clear();
  }
  canvas.SetLogy(false);
  widths[0].SetTitle("Missing-mass width vs momentum resolution;#sigma_{p}/p [%];RMS of MM [MeV/c^{2}]");
  for(int c=0;c<2;++c){widths[c].SetLineColor(c==0?kBlue+1:kRed+1);widths[c].SetMarkerColor(c==0?kBlue+1:kRed+1);widths[c].SetMarkerStyle(20+c);widths[c].SetLineWidth(2);}
  widths[0].SetMinimum(0.);widths[0].Draw("ALP");widths[1].Draw("LP same");
  TLegend legend(.16,.70,.49,.86);legend.SetBorderSize(0);legend.SetFillStyle(0);
  for(int c=0;c<2;++c)legend.AddEntry(&widths[c],labels[c],"lp");legend.Draw();
  canvas.Print(outputPdf);canvas.Print(TString(outputPdf)+"]");
  output.cd();widths[0].Write("mm_rms_lambda");widths[1].Write("mm_rms_sigma0");output.Close();
}
