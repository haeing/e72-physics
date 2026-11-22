// Run: root -l -b -q lambda_pipi_momentum_offset_study.cc
// Exclusive K- p -> Lambda pi+ pi-, Lambda -> p pi-.
// Reconstruct beam using (1) final pz sum, (2) final energy sum, separately.
// Truth-matched Lambda daughters; no combinatorial pairs, acceptance or resonances.
#include "lambda-pipi-phase.cc"
#include <TLegend.h>
#include <TLine.h>
#include <vector>
#include <algorithm>

namespace pipi_offset {
// ---------------- Parameters ----------------
constexpr double kBeamMomentumMeV = 735.;
constexpr int kEvents = 100000;
constexpr double kOffsetMinMeV = -20., kOffsetMaxMeV = 20.;
constexpr int kOffsetPoints = 81;
constexpr bool kApplyResolution = true;
constexpr double kProtonRelativeSigma = .04, kPionRelativeSigma = .04;
constexpr bool kConstrainLambdaMass = true;
constexpr double kDeltaBeamMinMeV = -150., kDeltaBeamMaxMeV = 150.;
// Four distinct daughter/production tracks plus a composite Lambda offset.
// Offset adds to the lab momentum MAGNITUDE; directions remain unchanged.
// Composite Lambda offset is applied after summing its measured daughter momenta.
const std::array<const char*,5> kModeNames{{"p from #Lambda", "#pi^{-} from #Lambda",
  "production #pi^{+}", "production #pi^{-}", "reconstructed #Lambda"}};
const std::array<Color_t,3> kColors{{kBlue+1,kBlack,kRed+1}};
struct Result {
  std::array<std::unique_ptr<TH2D>,2> scan;
  std::array<std::array<std::unique_ptr<TH1D>,3>,2> distribution;
  Long64_t invalidTracks=0, invalidEnergy=0, invalidLongitudinal=0;
  std::vector<std::array<double,2>> shifts;
  std::vector<std::array<Long64_t,2>> counts;
};
}

void lambda_pipi_momentum_offset_study(double beamMomentumMeV=pipi_offset::kBeamMomentumMeV,
    int events=pipi_offset::kEvents,
    double offsetMinMeV=pipi_offset::kOffsetMinMeV,
    double offsetMaxMeV=pipi_offset::kOffsetMaxMeV,
    int offsetPoints=pipi_offset::kOffsetPoints,
    bool applyResolution=pipi_offset::kApplyResolution,
    bool constrainLambdaMass=pipi_offset::kConstrainLambdaMass,
    const char* outputPdf="lambda_pipi_momentum_offset_study.pdf") {
  using namespace pipi_offset;
  using namespace lambda_pipi_phase_config;
  if(!std::isfinite(beamMomentumMeV)||beamMomentumMeV<=0.||events<=0||
     !std::isfinite(offsetMinMeV)||!std::isfinite(offsetMaxMeV)||
     offsetMaxMeV<=offsetMinMeV||offsetPoints<3)
    throw std::runtime_error("Invalid beam, event count or offset scan");
  gROOT->SetBatch(true);gStyle->SetOptStat(0);gStyle->SetPalette(1);
  TRandom3 random(20261026),noise(20261027);
  struct Guard {TRandom* old;~Guard(){gRandom=old;}} guard{gRandom};gRandom=&random;
  auto four=[](const TVector3& p,double m){TLorentzVector v;v.SetVectM(p,m);return v;};
  const auto target=four(TVector3(),kProtonMass);
  auto initial=four(TVector3(0.,0.,beamMomentumMeV/1000.),kKaonMass)+target;
  double masses[]={kLambdaMass,kPionMass,kPionMass},decayMasses[]={kProtonMass,kPionMass};
  TGenPhaseSpace production;
  if(!production.SetDecay(initial,3,masses))throw std::runtime_error("Lambda pi+ pi- production below threshold");
  const double step=(offsetMaxMeV-offsetMinMeV)/(offsetPoints-1);
  // Include exactly zero for a paired baseline even if the requested grid omits it.
  std::vector<double> offsets;for(int j=0;j<offsetPoints;++j)offsets.push_back(offsetMinMeV+j*step);
  bool hasZero=false;for(double v:offsets)if(std::abs(v)<1.e-10)hasZero=true;
  if(!hasZero)offsets.push_back(0.);
  const std::array<double,3> comparisonOffsets{{offsetMinMeV,0.,offsetMaxMeV}};
  const char* methods[]={"final p_{z} sum","final energy sum"};
  std::array<Result,5> results;
  for(int mode=0;mode<5;++mode) {
    auto& r=results[mode];r.shifts.resize(offsets.size());r.counts.resize(offsets.size());
    for(int method=0;method<2;++method) {
      r.scan[method].reset(new TH2D(Form("pipi_beam_scan_%d_%d",mode,method),
        Form("%s offset: %s; p_{K}^{reco}-p_{K}^{truth} [MeV/c];#delta p [MeV/c]",kModeNames[mode],methods[method]),
        400,kDeltaBeamMinMeV,kDeltaBeamMaxMeV,offsetPoints,offsetMinMeV-step/2.,offsetMaxMeV+step/2.));
      r.scan[method]->SetDirectory(nullptr);
      for(int panel=0;panel<3;++panel) {
        r.distribution[method][panel].reset(new TH1D(Form("pipi_beam_dist_%d_%d_%d",mode,method,panel),
          Form("%s offset: %s;p_{K}^{reco} [MeV/c];Counts",kModeNames[mode],methods[method]),
          400,beamMomentumMeV+kDeltaBeamMinMeV,beamMomentumMeV+kDeltaBeamMaxMeV));
        auto* h=r.distribution[method][panel].get();h->SetDirectory(nullptr);h->SetLineColor(kColors[panel]);h->SetLineWidth(2);
      }
    }
  }
  for(int event=0;event<events;++event) {
    Generate(production,random);auto lambda=*production.GetDecay(0);
    const auto pip=*production.GetDecay(1),pim=*production.GetDecay(2);
    TGenPhaseSpace decay;if(!decay.SetDecay(lambda,2,decayMasses))throw std::runtime_error("Lambda decay failed");
    Generate(decay,random);
    const std::array<TLorentzVector,4> truth{{*decay.GetDecay(0),*decay.GetDecay(1),pip,pim}};
    std::array<TVector3,4> directions;
    std::array<double,4> measured;
    for(int t=0;t<4;++t) {
      directions[t]=truth[t].Vect().Unit();
      const double sigma=t==0 ? kProtonRelativeSigma:kPionRelativeSigma;
      measured[t]=truth[t].P()*(applyResolution ? 1.+sigma*noise.Gaus():1.);
    }
    // Closure at truth before detector smearing and every systematic variation.
    const auto trueSum=truth[0]+truth[1]+truth[2]+truth[3];
    if(std::abs(trueSum.E()-initial.E())>1.e-9 || (trueSum-initial).P()>1.e-9)
      throw std::runtime_error("Truth four-momentum closure failed");
    for(int mode=0;mode<5;++mode)for(std::size_t j=0;j<offsets.size();++j) {
      const double offset=offsets[j];auto magnitudes=measured;
      if(mode<4)magnitudes[mode]+=offset/1000.;
      bool invalid=false;for(double p:magnitudes)if(p<0.)invalid=true;
      auto& r=results[mode];if(invalid){++r.invalidTracks;continue;}
      const auto proton=four(directions[0]*magnitudes[0],kProtonMass);
      const auto daughterPion=four(directions[1]*magnitudes[1],kPionMass);
      const auto lambdaSum=proton+daughterPion;
      auto usedLambda=constrainLambdaMass ? four(lambdaSum.Vect(),kLambdaMass):lambdaSum;
      if(mode==4) {
        const double p=usedLambda.P()+offset/1000.;
        if(p<0.){++r.invalidTracks;continue;}
        usedLambda=four(usedLambda.Vect().Unit()*p,constrainLambdaMass ? kLambdaMass:lambdaSum.M());
      }
      const auto finalSum=usedLambda+four(directions[2]*magnitudes[2],kPionMass)+four(directions[3]*magnitudes[3],kPionMass);
      const double energyK=finalSum.E()-kProtonMass;
      const std::array<double,2> beamRecon{{finalSum.Pz()*1000.,
        energyK>=kKaonMass ? std::sqrt(std::max(0.,energyK*energyK-kKaonMass*kKaonMass))*1000.:-1.}};
      for(int method=0;method<2;++method) {
        if(beamRecon[method]<0.){if(method==0)++r.invalidLongitudinal;else ++r.invalidEnergy;continue;}
        const double delta=beamRecon[method]-beamMomentumMeV;
        if(!applyResolution && std::abs(offset)<1.e-10 && std::abs(delta)>1.e-5)
          throw std::runtime_error("Zero-resolution/zero-offset beam reconstruction closure failed");
        if(j<static_cast<std::size_t>(offsetPoints))r.scan[method]->Fill(delta,offset);
        r.shifts[j][method]+=delta;++r.counts[j][method];
        for(int panel=0;panel<3;++panel)
          if(std::abs(offset-comparisonOffsets[panel])<1.e-10)r.distribution[method][panel]->Fill(beamRecon[method]);
      }
    }
  }
  TString logPdf(outputPdf);if(logPdf.EndsWith(".pdf"))logPdf.Remove(logPdf.Length()-4);logPdf+="_log.pdf";
  TCanvas canvas("c_pipi_beam_offsets","Beam reconstruction vs momentum offsets",1600,1400);
  for(bool logarithmic:{false,true}) {
    const TString pdf=logarithmic ? logPdf:TString(outputPdf);canvas.Print(pdf+"[");
    for(int mode=0;mode<5;++mode) {
      canvas.Divide(2,2);std::array<TLine,4> references;
      std::array<TLegend,2> legends{{TLegend(.57,.66,.88,.88),TLegend(.57,.66,.88,.88)}};
      for(int method=0;method<2;++method) {
        canvas.cd(method+1);gPad->SetLeftMargin(.17);gPad->SetRightMargin(.18);gPad->SetBottomMargin(.15);gPad->SetLogz(logarithmic);
        auto* h=results[mode].scan[method].get();h->SetMinimum(logarithmic ? .5:0.);h->Draw("colz");
        references[method]=TLine(0.,offsetMinMeV-step/2.,0.,offsetMaxMeV+step/2.);
        references[method].SetLineColor(kRed);references[method].SetLineWidth(1);references[method].Draw("same");
        canvas.cd(method+3);gPad->SetLeftMargin(.17);gPad->SetRightMargin(.12);gPad->SetBottomMargin(.15);gPad->SetLogy(logarithmic);
        double maximum=0.;for(auto& plot:results[mode].distribution[method])maximum=std::max(maximum,plot->GetMaximum());
        auto& legend=legends[method];legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.030);
        for(int panel=0;panel<3;++panel) {
          auto* plot=results[mode].distribution[method][panel].get();plot->SetMinimum(logarithmic ? .5:0.);plot->SetMaximum(std::max(1.,maximum*(logarithmic ? 2.:1.2)));
          plot->Draw(panel==0 ? "hist":"hist same");legend.AddEntry(plot,Form("%+.0f MeV/c",comparisonOffsets[panel]),"l");
        }
        legend.Draw();references[method+2]=TLine(beamMomentumMeV,logarithmic ? .5:0.,beamMomentumMeV,std::max(1.,maximum));
        references[method+2].SetLineColor(kRed);references[method+2].SetLineStyle(2);references[method+2].Draw("same");
      }
      canvas.Print(pdf);canvas.Clear();
    }
    canvas.Print(pdf+"]");
  }
  for(int mode=0;mode<5;++mode) {
    const auto& r=results[mode];
    for(std::size_t j=0;j<offsets.size();++j)
      if(j==0 || j==static_cast<std::size_t>(offsetPoints-1) || std::abs(offsets[j])<1.e-10)
        Info("pipi_beam_offset","%s offset %+.1f: mean beam shift pz %+.4f, energy %+.4f MeV/c",kModeNames[mode],offsets[j],r.counts[j][0] ? r.shifts[j][0]/r.counts[j][0]:0.,r.counts[j][1] ? r.shifts[j][1]/r.counts[j][1]:0.);
    Info("pipi_beam_offset","%s: invalid tracks %lld; invalid energy %lld; negative reconstructed pz %lld",kModeNames[mode],r.invalidTracks,r.invalidEnergy,r.invalidLongitudinal);
  }
}
