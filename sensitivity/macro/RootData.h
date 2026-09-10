#ifndef SENSITIVITY_ROOT_DATA_H
#define SENSITIVITY_ROOT_DATA_H
#include <TFile.h>
#include <TDirectory.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TParameter.h>
#include <TSystem.h>
#include <cmath>
#include <stdexcept>
namespace RootData {
inline TGraphErrors *graph(TFile &file, const char *key, int n) {
  auto *g = dynamic_cast<TGraphErrors *>(file.Get(key));
  if (file.IsZombie() || !g || g->GetN() != n)
    throw std::runtime_error(std::string("Missing or incompatible ROOT graph: ") + file.GetName() + ":" + key);
  return g;
}
inline void readMomentum(const char *path, const char *key, int n, const double *mom, double *values) {
  TFile file(path, "READ");
  auto *g = graph(file, key, n);
  for (int i=0; i<n; ++i) {
    if (std::abs(g->GetX()[i]-mom[i]) > 1e-6 || !std::isfinite(g->GetY()[i]))
      throw std::runtime_error("Invalid momentum grid or value");
    values[i] = g->GetY()[i];
  }
}
template<int M, int C>
void readAngular(const char *path, const char *key, const double *mom, double (&values)[M][C]) {
  TFile file(path, "READ");
  for (int i=0; i<M; ++i) {
    auto *g = graph(file, Form("momentum/p%03.0f/%s",mom[i],key), C);
    for (int j=0; j<C; ++j) {
      if (std::abs(g->GetX()[j]-(-1.+(j+0.5)*2./C)) > 1e-6 || !std::isfinite(g->GetY()[j]))
        throw std::runtime_error("Invalid angular grid or value");
      values[i][j] = g->GetY()[j];
    }
  }
}
inline TDirectory *momentum(TFile &file, double p) {
  auto *base=file.GetDirectory("momentum");
  if (!base) base=file.mkdir("momentum");
  auto *dir=base->mkdir(Form("p%03.0f",p));
  dir->cd();
  TParameter<double>("momentum_MeV_c",p).Write();
  TParameter<double>("half_width_MeV_c",1.).Write();
  return dir;
}
inline void angular(TGraphErrors *g, const char *name) {
  g->Write(name);
  TH1D h(Form("h_%s",name),g->GetTitle(),g->GetN(),-1.,1.);
  h.SetDirectory(nullptr);
  for (int j=0; j<g->GetN(); ++j) {
    h.SetBinContent(j+1,g->GetY()[j]);
    h.SetBinError(j+1,g->GetEY()[j]);
  }
  h.Write();
}
}
#endif
