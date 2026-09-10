#include <TCanvas.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TCollection.h>
#include <TString.h>

// Run each stage in its own ROOT process: legacy globals share names.
void run_stage(const char *stage) {
  gROOT->SetBatch(kTRUE);
  gSystem->mkdir("results", kTRUE);
  Int_t error = 0;
  gROOT->ProcessLine(Form(".x macro/%s.cc", stage), &error);
  if (error) gSystem->Exit(1);
  TFile output(Form("results/%s.root", stage), "UPDATE");
  TIter next(gROOT->GetListOfCanvases());
  while (auto *object = next()) {
    auto *canvas = dynamic_cast<TCanvas *>(object);
    if (!canvas) continue;
    canvas->Write(canvas->GetName(), TObject::kOverwrite);
  }
  output.Close();
}
