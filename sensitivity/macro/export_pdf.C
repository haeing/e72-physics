#include <TCanvas.h>
#include <TFile.h>
#include <TKey.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TString.h>
#include <memory>

void export_pdf() {
  gROOT->SetBatch(kTRUE);
  const char *pdf = "results/sensitivity.pdf";
  const char *stages[] = {"yield", "acceptance", "diff_cross", "polar_th", "crystalball"};
  TCanvas book("pdf_book", "Sensitivity");
  book.Print(Form("%s[", pdf));
  int pages = 0;
  for (const char *stage : stages) {
    TFile input(Form("results/%s.root", stage), "READ");
    if (input.IsZombie()) gSystem->Exit(1);
    TIter next(input.GetListOfKeys());
    while (auto *key = dynamic_cast<TKey *>(next())) {
      std::unique_ptr<TObject> object(key->ReadObj());
      auto *canvas = dynamic_cast<TCanvas *>(object.get());
      if (!canvas) continue;
      canvas->Print(pdf, Form("Title:%s / %s", stage, canvas->GetName()));
      ++pages;
    }
  }
  book.Print(Form("%s]", pdf));
  if (!pages) gSystem->Exit(1);
  Printf("Saved %d pages to %s", pages, pdf);
}
