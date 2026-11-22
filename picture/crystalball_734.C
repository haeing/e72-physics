#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <stdexcept>

// Run: root -l -b -q crystalball_734.C
// Read the original Crystal Ball points only (no projected E72 errors or fit).
void crystalball_734()
{
    const TString directory = gSystem->DirName(__FILE__);
    TFile input(directory + "/../sensitivity/results/crystalball.root", "READ");
    if (input.IsZombie())
        throw std::runtime_error("Cannot open sensitivity/results/crystalball.root");

    TGraphErrors *source = nullptr;
    input.GetObject("momentum/p734/crystal_ball", source);
    if (!source)
        throw std::runtime_error("Missing momentum/p734/crystal_ball");
    auto graph = static_cast<TGraphErrors *>(source->Clone("crystal_ball_734"));
    input.Close();

    gStyle->SetOptStat(0);
    auto canvas = new TCanvas("crystalball_734_canvas", "Crystal Ball: 734 MeV/c", 850, 700);
    canvas->SetLeftMargin(0.15);
    canvas->SetBottomMargin(0.14);
    graph->SetTitle("Crystal Ball, p_{K^{-}} = 734 MeV/c;cos#theta_{#eta}^{CM};d#sigma/d#Omega [mb/sr]");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.2);
    graph->SetMarkerColor(kBlue + 2);
    graph->SetLineColor(kBlue + 2);
    graph->SetLineWidth(2);
    graph->Draw("APZ");
    graph->GetXaxis()->SetLimits(-1., 1.);
    graph->GetXaxis()->SetRange(0, 0); // Clear the persisted source-axis bin selection.
    graph->SetMinimum(0.);
    for (auto axis : {graph->GetXaxis(), graph->GetYaxis()}) {
        axis->SetTitleSize(0.05);
        axis->SetLabelSize(0.045);
    }
    graph->GetYaxis()->SetTitleOffset(1.35);
    canvas->Modified();
    canvas->Update();
    canvas->SaveAs(directory + "/crystalball_734.png");
    canvas->SaveAs(directory + "/crystalball_734.pdf");
    canvas->SaveAs(directory + "/crystalball_734.svg");
}
