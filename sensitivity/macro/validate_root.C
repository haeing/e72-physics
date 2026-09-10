#include "RootData.h"
#include <fstream>
#include <sstream>
#include <TCanvas.h>
#include <TKey.h>
void validate_root() {
  const double momenta[]={724,726,728,730,732,734,738,742,746,750,754,758,762,766,770};
  TFile y("results/yield.root"), a("results/acceptance.root"), d("results/diff_cross.root"),
        p("results/polar_th.root"), c("results/crystalball.root");
  auto *entry=RootData::graph(y,"entry_mom",15);
  auto *lumi=RootData::graph(y,"luminosity_mom",15);
  int graphs=0;
  for (int i=0;i<15;++i) {
    auto *acc=RootData::graph(a,Form("momentum/p%03.0f/acceptance",momenta[i]),15);
    auto *yield=RootData::graph(d,Form("momentum/p%03.0f/yield",momenta[i]),15);
    auto *pol=RootData::graph(p,Form("momentum/p%03.0f/polarization_uncertainty",momenta[i]),15);
    auto *cs=RootData::graph(c,Form("momentum/p%03.0f/differential_cs",momenta[i]),15);
    graphs+=4;
    double sum=0;
    for(int j=0;j<15;++j) {
      double n=yield->GetY()[j], efficiency=acc->GetY()[j];
      if(!(n>0 && efficiency>0 && efficiency<=1)) throw std::runtime_error("Invalid yield/acceptance");
      sum+=n;
      if(std::abs(pol->GetY()[j]-sqrt(3)/(0.65*sqrt(n)))>1e-12)
        throw std::runtime_error("Polarization uncertainty mismatch");
      double error=sqrt(n)/(lumi->GetY()[i]*(2./15)*2*TMath::Pi()*efficiency);
      if(!std::isfinite(cs->GetY()[j]) || std::abs(cs->GetEY()[j]-error)>1e-12)
        throw std::runtime_error("Differential CS error mismatch");
    }
    if(std::abs(sum-entry->GetY()[i])>1e-8) throw std::runtime_error("Yield conservation failed");
  }
  // Compare the stored full-precision values against the previous TXT pipeline.
  for (const char *name : {"entry_mom", "luminosity_mom", "entry_mom_cos", "acceptance_mom_cos"}) {
    std::ifstream in(Form("legacy/last_txt_results/%s.txt",name));
    if(!in) throw std::runtime_error("Missing comparison snapshot");
    std::string line; int row=0;
    while(std::getline(in,line)) {
      if(line.empty() || line[0]=='#') continue;
      std::istringstream stream(line); double v; std::vector<double> values;
      while(stream>>v) values.push_back(v);
      double expected=values.back(), actual;
      if(values.size()==4) actual=(std::string(name)=="entry_mom"?entry:lumi)->GetY()[row];
      else {
        TFile &file=std::string(name)=="entry_mom_cos"?d:a;
        const char *key=std::string(name)=="entry_mom_cos"?"yield":"acceptance";
        actual=RootData::graph(file,Form("momentum/p%03.0f/%s",values[0],key),15)->GetY()[row%15];
      }
      if(std::abs(actual-expected)>1e-5*std::max(1.,std::abs(expected)))
        throw std::runtime_error("Regression against previous TXT result");
      ++row;
    }
  }
  Printf("PASS: %d momentum-specific graphs, yield conservation, uncertainty formulas, TXT regression",graphs);
}
