#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TLegend.h>
#include <TColor.h>
#include <vector>
#include <TStyle.h>
#include <TPad.h>
#include <TString.h>
#include <cmath>
#include <iostream>
#include <stdexcept>

// Run: root -l -b -q 'RunTimeStat.C()'
// Existing full yield is BEFORE the vertex cut, default: table total 70379.
// Use RunTimeStat(9341) for the analyzed subset instead.
// Optional existingFullVariance is AFTER-cut variance; -1 uses Poisson.
// Existing and future data are assumed to have the same rates/composition.
void RunTimeStat(double existingFullYieldBeforeCut = 70379.,
                 double existingEmptyHours = 0.,
                 double frameFractionBeforeCut = 30289. / 70379.,
                 double fullRateBeforeCut = 936.,
                 double emptyRateBeforeCut = 400.,
                 double hydrogenCutEfficiency = 0.80,
                 double frameCutEfficiency = 0.20,
                 double existingFullVariance = -1.)
{
    if (existingFullYieldBeforeCut < 0 || existingEmptyHours < 0 ||
        !std::isfinite(existingFullVariance) ||
        (existingFullVariance < 0 && existingFullVariance != -1.) ||
        frameFractionBeforeCut <= 0 || frameFractionBeforeCut >= 1 ||
        fullRateBeforeCut <= 0 || emptyRateBeforeCut <= 0 ||
        hydrogenCutEfficiency <= 0 || hydrogenCutEfficiency > 1 ||
        frameCutEfficiency <= 0 || frameCutEfficiency > 1)
        throw std::invalid_argument("Invalid time, rate, fraction or efficiency.");

    const double existingFullHours = existingFullYieldBeforeCut / fullRateBeforeCut;
    const double initialSignal = existingFullYieldBeforeCut *
        (1. - frameFractionBeforeCut) * hydrogenCutEfficiency;
    const double initialBackground = existingFullYieldBeforeCut *
        frameFractionBeforeCut * frameCutEfficiency;
    const double initialFull = initialSignal + initialBackground;
    const double initialVariance = existingFullVariance >= 0 ? existingFullVariance : initialFull;

    // All rates below are selected events per hour AFTER the vertex cut.
    const double signalRate = fullRateBeforeCut *
        (1. - frameFractionBeforeCut) * hydrogenCutEfficiency;
    const double backgroundRate = fullRateBeforeCut *
        frameFractionBeforeCut * frameCutEfficiency;
    const double fullRate = signalRate + backgroundRate;
    const double emptyRate = emptyRateBeforeCut * frameCutEfficiency;

    // Background-rate transfer factor, distinct from exposure/time scaling.
    // In an analysis this must be validated with beam normalization and
    // full/empty transfer corrections; it is treated as exact in this forecast.
    const double k = backgroundRate / emptyRate;

    // Use 0.1-hour bins over 0..24 physics hours and 3..12 empty hours.
    // Zero total exposure is excluded because its uncertainty is undefined.
    const int nx = 240, ny = 90;
    auto total = new TH2D("hTotal", "Total statistical uncertainty;Additional Physics run time [h];Empty-target run time [h];#sigma_{S}/S [%]",
                          nx, 0., 24., ny, 3., 12.);
    auto full = new TH2D("hFull", "Physics statistical contribution;Additional Physics run time [h];Empty-target run time [h];#sigma_{S,F}/S [%]",
                         nx, 0., 24., ny, 3., 12.);
    auto empty = new TH2D("hEmpty", "Empty statistical contribution;Additional Physics run time [h];Empty-target run time [h];#sigma_{S,E}/S [%]",
                          nx, 0., 24., ny, 3., 12.);
    TH2D* maps[] = {total, full, empty};
    for (auto h : maps) { h->SetDirectory(nullptr); h->SetStats(false); }

    total->SetMinimum(0.7);
    total->SetMaximum(1.5);

    for (int ix = 1; ix <= nx; ++ix) {
        const double addF = total->GetXaxis()->GetBinCenter(ix);
        const double tF = existingFullHours + addF;
        for (int iy = 1; iy <= ny; ++iy) {
            const double tE = existingEmptyHours + total->GetYaxis()->GetBinCenter(iy);
            const double expectedSignal = initialSignal + signalRate * addF;
            const double alpha = k * tF / tE;
            const double varianceFull = initialVariance + fullRate * addF;
            const double varianceEmpty = alpha * alpha * emptyRate * tE;
            const double dFull = std::sqrt(varianceFull) / expectedSignal;
            const double dEmpty = std::sqrt(varianceEmpty) / expectedSignal;
            total->SetBinContent(ix, iy, 100. * std::hypot(dFull, dEmpty));
            full->SetBinContent(ix, iy, 100. * dFull);
            empty->SetBinContent(ix, iy, 100. * dEmpty);
        }
    }

    gStyle->SetPalette(kRainBow);
    gStyle->SetNumberContours(100);
    auto canvas = new TCanvas("cRunTimeStat", "Run-time statistical forecast", 1800, 650);
    canvas->Divide(3, 1);
    for (int i = 0; i < 3; ++i) {
        canvas->cd(i + 1);
        gPad->SetRightMargin(.19);
        gPad->SetLeftMargin(.12);
        gPad->SetBottomMargin(.13);
        maps[i]->GetZaxis()->SetTitleOffset(1.3);
        maps[i]->Draw("COLZ");

    }
    canvas->SaveAs("RunTimeStat.png");
    canvas->SaveAs("RunTimeStat.pdf");

    // Fixed additional budget; retain the existing data in the forecast.
    auto allocation = new TGraph();
    allocation->SetName("gFixed24h");
    allocation->SetTitle("24 h additional run time: Physics + Empty = 24 h;Empty-target run time [h];Total relative statistical uncertainty [%]");
    double bestEmpty = 0., bestError = 1.e100;
    for (int j = 1; j <= 2400; ++j) {
        const double addE = j * 0.01;
        const double addF = 24. - addE;
        const double signal = initialSignal + signalRate * addF;
        if (signal <= 0) continue;
        const double tE = existingEmptyHours + addE;
        const double alpha = k * (existingFullHours + addF) / tE;
        const double variance = initialVariance + fullRate * addF
                              + alpha * alpha * emptyRate * tE;
        const double error = 100. * std::sqrt(variance) / signal;
        allocation->SetPoint(allocation->GetN(), addE, error);
        if (error < bestError) { bestError = error; bestEmpty = addE; }
    }
    auto allocationCanvas = new TCanvas("cFixed24h", "24-hour run allocation", 1000, 700);
    allocationCanvas->SetLeftMargin(.13);
    allocationCanvas->SetBottomMargin(.13);
    allocationCanvas->SetGrid();
    allocation->SetLineColor(kBlue + 1);
    allocation->SetLineWidth(3);
    allocation->Draw("AL");
    allocation->GetXaxis()->SetLimits(0., 24.);
    // Show the useful allocation range; the curve diverges at zero empty time.
    allocation->SetMinimum(0.);
    allocation->SetMaximum(3.);
    auto label = new TLatex();
    label->SetNDC();
    label->SetTextSize(.03);
    allocationCanvas->SaveAs("RunTimeStat_24h.png");
    allocationCanvas->SaveAs("RunTimeStat_24h.pdf");
    std::cout << "24 h allocation optimum: empty=" << bestEmpty
              << " h, physics=" << 24. - bestEmpty
              << " h, total uncertainty=" << bestError << "%\n";
    auto comparisonCanvas = new TCanvas("cAllocation12to24h", "12-24 hour run allocation", 1100, 750);
    comparisonCanvas->SetLeftMargin(.13);
    comparisonCanvas->SetBottomMargin(.13);
    comparisonCanvas->SetGrid();
    auto legend = new TLegend(.57, .48, .88, .86);
    legend->SetHeader("Additional Physics + Empty");
    legend->SetNColumns(2);
    legend->SetBorderSize(0);
    legend->SetTextSize(.025);
    std::vector<TGraph*> allocationCurves;
    for (int budget = 12; budget <= 24; budget += 6) {
        auto curve = new TGraph();
        curve->SetName(TString::Format("gAllocation%dh", budget));
        curve->SetTitle("12-24 h additional run time;Empty-target run time [h];Total relative statistical uncertainty [%]");
        for (int j = 1; j <= budget * 100; ++j) {
            const double addE = j * .01;
            const double addF = budget - addE;
            const double signal = initialSignal + signalRate * addF;
            if (signal <= 0) continue;
            const double tE = existingEmptyHours + addE;
            const double alpha = k * (existingFullHours + addF) / tE;
            const double variance = initialVariance + fullRate * addF
                                  + alpha * alpha * emptyRate * tE;
            curve->SetPoint(curve->GetN(), addE, 100. * std::sqrt(variance) / signal);
        }
        const int color = TColor::GetColorPalette((budget - 12) * (TColor::GetNumberOfColors() - 1) / 12);
        curve->SetLineColor(color);
        curve->SetLineWidth(2);
        curve->SetLineStyle(1 + (budget - 12) / 6);
        curve->Draw(budget == 12 ? "AL" : "L SAME");
        if (budget == 12) {
            curve->GetXaxis()->SetLimits(0., 24.);
            curve->SetMinimum(0.);
            curve->SetMaximum(3.);
        }

        legend->AddEntry(curve, TString::Format("%d h", budget), "l");
        allocationCurves.push_back(curve);
    }
    legend->Draw();
    comparisonCanvas->SaveAs("RunTimeStat_12to24h.png");
    comparisonCanvas->SaveAs("RunTimeStat_12to24h.pdf");
    TFile output("RunTimeStat.root", "RECREATE");
    for (auto h : maps) h->Write();
    canvas->Write();
    allocation->Write();
    allocationCanvas->Write();
    for (auto curve : allocationCurves) curve->Write();
    comparisonCanvas->Write();
    output.Close();

    std::cout << "Post-cut rates [events/h]: signal=" << signalRate
              << ", frame=" << backgroundRate << ", full=" << fullRate
              << ", empty=" << emptyRate << "; k=" << k << '\n';
    std::cout << "Existing equivalent hours: full=" << existingFullHours
              << ", empty=" << existingEmptyHours << '\n';
    std::cout << "Existing yields: pre-cut=" << existingFullYieldBeforeCut
              << ", post-cut signal=" << initialSignal << ", frame=" << initialBackground
              << "; post-cut full variance=" << initialVariance << '\n';
    std::cout << "ADDITIONAL-time examples:\n";
    for (double addF : {0., 20., 50.}) {
        for (double addE : {6., 12., 24.}) {
            const double tE = existingEmptyHours + addE;
            const double alpha = k * (existingFullHours + addF) / tE;
            const double signal = initialSignal + signalRate * addF;
            const double variance = initialVariance + fullRate * addF
                                  + alpha * alpha * emptyRate * tE;
            std::cout << "  +physics=" << addF << " h, +empty=" << addE
                      << " h: relative statistical error=" << 100.*std::sqrt(variance)/signal << "%\n";
        }
    }
}
