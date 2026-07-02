#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <cmath>
#include <iomanip>


#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <TObject.h>
#include <TROOT.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>

#include "../../basic-property.hh"

using namespace std;

const double br_Sigma0_to_Lambda = 1.0;

const vector<int> runnumbers = {2447, 2449, 2450, 2451, 2452, 2453, 2454, 2456, 2457, 2458, 2459, 2460};

struct ScalerInfo {
    Long64_t trigD;
    double daqEff;
};

struct CombinedScalerInfo {
    vector<int> runnumbers;
    Long64_t trigD;
    double daqEff;
};

struct TriggerAcceptanceInfo {
    string reactionName;
    string fileName;
    string graphName;
    Long64_t nEntry;
    Long64_t nBeam;
    Long64_t nFlag1000;
    Long64_t nFlag1101;
    Long64_t nFlag1011;
    Long64_t nFlag1111;
    double pctFlag1101;
    double pctFlag1011;
    double pctFlag1111;
    vector<double> binLow;
    vector<double> binHigh;
    vector<double> binCenter;
    vector<Long64_t> binEntry;
    vector<Long64_t> binDen;
    vector<Long64_t> binNum;
    vector<double> binAcc;
    vector<double> binAccErr;
    vector<double> binCrossSection;
};

struct BeamSpectrum {
    string fileName;
    Long64_t nEntry;
    double scale;
    vector<double> binLow;
    vector<double> binHigh;
    vector<double> binCenter;
    vector<double> binRaw;
    vector<double> binScaled;
};

struct YieldEstimateInfo {
    string reactionName;
    vector<double> binLow;
    vector<double> binHigh;
    vector<double> binCenter;
    vector<double> nBeam;
    vector<double> crossSectionMb;
    vector<double> geomEff;
    vector<double> yield;
    double totalYield;
};

ScalerInfo GetScalerInfo(int runnumber)
{
    string scaler_dir =
        "/gpfs/group/had/sks/E72/JPARC2025Nov/share/scaler";

    string scaler_file =
        Form("%s/scaler_%05d.txt",
             scaler_dir.c_str(), runnumber);

    ifstream fin(scaler_file);

    if (!fin.is_open()) {
        cerr << "Cannot open " << scaler_file << endl;
        return {-1, -1.0};
    }

    string line;
    Long64_t trigD = -1;
    double daqEff = -1.0;

    while (getline(fin, line)) {

        stringstream ss(line);

        string name;
        double value;

        if (!(ss >> name))
            continue;

        if (!(ss >> value))
            continue;

        if (name == "TRIG-D")
            trigD = static_cast<Long64_t>(value);
        else if (name == "DAQ-Eff")
            daqEff = value;

        if (trigD >= 0 && daqEff >= 0.0)
            return {trigD, daqEff};
    }

    if (trigD < 0)
        cerr << "TRIG-D not found in " << scaler_file << endl;
    if (daqEff < 0.0)
        cerr << "DAQ-Eff not found in " << scaler_file << endl;

    return {trigD, daqEff};
}

double NormalizeScalerEfficiency(double eff)
{
    if (eff > 1.0)
        return eff / 100.0;
    return eff;
}

CombinedScalerInfo GetCombinedScalerInfo(const vector<int>& runs)
{
    CombinedScalerInfo combined;
    combined.runnumbers = runs;
    combined.trigD = 0;
    combined.daqEff = 0.0;

    double weightedDaqEff = 0.0;

    cout << "runs : ";
    for (size_t i = 0; i < runs.size(); ++i) {
        if (i > 0)
            cout << ", ";
        cout << runs[i];
    }
    cout << endl;

    for (int run : runs) {
        ScalerInfo scaler = GetScalerInfo(run);
        if (scaler.trigD <= 0 || scaler.daqEff < 0.0) {
            cerr << "Skip run " << run << " because scaler info is invalid" << endl;
            continue;
        }

        double daqEff = NormalizeScalerEfficiency(scaler.daqEff);
        combined.trigD += scaler.trigD;
        weightedDaqEff += scaler.trigD * daqEff;

        cout << "  run " << run
             << " TRIG-D : " << scaler.trigD
             << ", DAQ-Eff : " << daqEff << endl;
    }

    if (combined.trigD > 0)
        combined.daqEff = weightedDaqEff / combined.trigD;

    return combined;
}
TGraph *LoadCrossSectionGraph(TFile *file, const string& graph_name)
{
    if (!file || file->IsZombie())
        return nullptr;

    TObject *obj = file->Get(graph_name.c_str());
    TGraph *graph = dynamic_cast<TGraph*>(obj);
    if (!graph)
        cerr << graph_name << " not found in cross-section file" << endl;

    return graph;
}

double InterpolateGraph(TGraph *graph, double x)
{
    if (!graph || graph->GetN() <= 0)
        return 0.0;

    vector<pair<double, double>> points;
    for (int i = 0; i < graph->GetN(); ++i) {
        double gx = 0.0;
        double gy = 0.0;
        graph->GetPoint(i, gx, gy);
        points.push_back({gx, gy});
    }

    sort(points.begin(), points.end());

    if (x < points.front().first || x > points.back().first)
        return 0.0;

    for (size_t i = 1; i < points.size(); ++i) {
        if (x > points[i].first)
            continue;

        double x0 = points[i - 1].first;
        double y0 = points[i - 1].second;
        double x1 = points[i].first;
        double y1 = points[i].second;

        if (x1 == x0)
            return y0;

        double t = (x - x0) / (x1 - x0);
        return y0 + t * (y1 - y0);
    }

    return points.back().second;
}

int FindMomentumBin(double mom, const vector<double>& bins)
{
    if (bins.size() < 2 || mom < bins.front() || mom >= bins.back())
        return -1;

    auto it = upper_bound(bins.begin(), bins.end(), mom);
    return static_cast<int>(it - bins.begin()) - 1;
}

BeamSpectrum GetBeamSpectrum(const string& beam_file,
                             const vector<double>& mom_bins,
                             double nKbeam)
{
    BeamSpectrum beam;
    beam.fileName = beam_file;
    beam.nEntry = 0;
    beam.scale = 0.0;

    const int nBins = static_cast<int>(mom_bins.size()) - 1;
    for (int i = 0; i < nBins; ++i) {
        double low = mom_bins[i];
        double high = mom_bins[i + 1];
        beam.binLow.push_back(low);
        beam.binHigh.push_back(high);
        beam.binCenter.push_back(0.5 * (low + high));
        beam.binRaw.push_back(0.0);
        beam.binScaled.push_back(0.0);
    }

    TFile *file = TFile::Open(beam_file.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cerr << "Cannot open " << beam_file << endl;
        if (file)
            file->Close();
        return beam;
    }

    TTree *tree = dynamic_cast<TTree*>(file->Get("g4hyptpc_light"));
    if (!tree)
        tree = dynamic_cast<TTree*>(file->Get("g4hyptpc"));

    if (!tree) {
        cerr << "g4hyptpc_light or g4hyptpc not found in " << beam_file << endl;
        file->Close();
        return beam;
    }

    if (!tree->GetBranch("mom_kaon_lab")) {
        cerr << "mom_kaon_lab not found in " << beam_file << endl;
        file->Close();
        return beam;
    }

    Double_t momKaonLab = 0.0;
    tree->SetBranchAddress("mom_kaon_lab", &momKaonLab);

    beam.nEntry = tree->GetEntries();
    for (Long64_t i = 0; i < beam.nEntry; ++i) {
        tree->GetEntry(i);
        int bin = FindMomentumBin(momKaonLab, mom_bins);
        if (bin >= 0)
            beam.binRaw[bin] += 1.0;
    }

    double rawIntegral = 0.0;
    for (double count : beam.binRaw)
        rawIntegral += count;

    if (rawIntegral > 0.0 && nKbeam > 0.0) {
        beam.scale = nKbeam / rawIntegral;
        for (size_t i = 0; i < beam.binRaw.size(); ++i)
            beam.binScaled[i] = beam.binRaw[i] * beam.scale;
    }

    file->Close();
    return beam;
}

TriggerAcceptanceInfo GetTriggerAcceptance(const string& reaction_name,
                                           const string& sim_file,
                                           const string& graph_name,
                                           TGraph *cross_section_graph,
                                           const vector<double>& mom_bins)
{
    TriggerAcceptanceInfo info;
    info.reactionName = reaction_name;
    info.fileName = sim_file;
    info.graphName = graph_name;
    info.nEntry = 0;
    info.nBeam = 0;
    info.nFlag1000 = 0;
    info.nFlag1101 = 0;
    info.nFlag1011 = 0;
    info.nFlag1111 = 0;
    info.pctFlag1101 = 0.0;
    info.pctFlag1011 = 0.0;
    info.pctFlag1111 = 0.0;

    const int nBins = static_cast<int>(mom_bins.size()) - 1;
    for (int i = 0; i < nBins; ++i) {
        double low = mom_bins[i];
        double high = mom_bins[i + 1];
        double center = 0.5 * (low + high);
        info.binLow.push_back(low);
        info.binHigh.push_back(high);
        info.binCenter.push_back(center);
        info.binEntry.push_back(0);
        info.binDen.push_back(0);
        info.binNum.push_back(0);
        info.binAcc.push_back(0.0);
        info.binAccErr.push_back(0.0);
        info.binCrossSection.push_back(InterpolateGraph(cross_section_graph, center));
    }

    TFile *file = TFile::Open(sim_file.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cerr << "Cannot open " << sim_file << endl;
        if (file)
            file->Close();
        return info;
    }

    TTree *tree = dynamic_cast<TTree*>(file->Get("g4hyptpc_light"));
    if (!tree) {
        cerr << "g4hyptpc_light not found in " << sim_file << endl;
        file->Close();
        return info;
    }

    if (!tree->GetBranch("trig_flag") || !tree->GetBranch("mom_kaon_lab")) {
        cerr << "trig_flag or mom_kaon_lab not found in " << sim_file << endl;
        file->Close();
        return info;
    }

    Int_t trigFlag = 0;
    Double_t momKaonLab = 0.0;
    tree->SetBranchAddress("trig_flag", &trigFlag);
    tree->SetBranchAddress("mom_kaon_lab", &momKaonLab);

    const Int_t kBeamFlag = 0x8;  // 1000
    const Int_t kFlag1000 = 0x8;
    const Int_t kFlag1101 = 0xd;
    const Int_t kFlag1011 = 0xb;
    const Int_t kFlag1111 = 0xf;

    info.nEntry = tree->GetEntries();
    for (Long64_t i = 0; i < info.nEntry; ++i) {
        tree->GetEntry(i);

        int bin = FindMomentumBin(momKaonLab, mom_bins);
        if (bin >= 0)
            ++info.binEntry[bin];

        if ((trigFlag & kBeamFlag) != 0)
            ++info.nBeam;

        if (trigFlag == kFlag1000)
            ++info.nFlag1000;
        else if (trigFlag == kFlag1101)
            ++info.nFlag1101;
        else if (trigFlag == kFlag1011)
            ++info.nFlag1011;
        else if (trigFlag == kFlag1111)
            ++info.nFlag1111;

        if (bin < 0)
            continue;

        bool isDen = (trigFlag == kFlag1000 || trigFlag == kFlag1101 ||
                      trigFlag == kFlag1011 || trigFlag == kFlag1111);
        bool isNum = (trigFlag == kFlag1101 || trigFlag == kFlag1011 ||
                      trigFlag == kFlag1111);

        if (isDen)
            ++info.binDen[bin];
        if (isNum)
            ++info.binNum[bin];
    }

    Long64_t nDenTotal = info.nFlag1000 + info.nFlag1101 + info.nFlag1011 + info.nFlag1111;
    if (nDenTotal > 0) {
        info.pctFlag1101 = 100.0 * info.nFlag1101 / nDenTotal;
        info.pctFlag1011 = 100.0 * info.nFlag1011 / nDenTotal;
        info.pctFlag1111 = 100.0 * info.nFlag1111 / nDenTotal;
    }

    for (size_t i = 0; i < info.binAcc.size(); ++i) {
        if (info.binDen[i] > 0) {
            double acc = static_cast<double>(info.binNum[i]) / info.binDen[i];
            info.binAcc[i] = 100.0 * acc;
            info.binAccErr[i] = 100.0 * sqrt(acc * (1.0 - acc) / info.binDen[i]);
        }
    }

    file->Close();
    return info;
}

TH1D *MakeCountHist(const TriggerAcceptanceInfo& info, const char *name, const char *title, bool use_denominator)
{
    const int nBins = info.binLow.size();
    vector<double> edges;
    for (int i = 0; i < nBins; ++i)
        edges.push_back(info.binLow[i]);
    edges.push_back(info.binHigh.back());

    TH1D *hist = new TH1D(name, title, nBins, edges.data());
    hist->SetDirectory(nullptr);
    for (int i = 0; i < nBins; ++i) {
        Long64_t count = use_denominator ? info.binDen[i] : info.binEntry[i];
        hist->SetBinContent(i + 1, count);
    }
    return hist;
}

TH1D *MakeBeamHist(const BeamSpectrum& beam, const char *name, const char *title)
{
    const int nBins = beam.binLow.size();
    vector<double> edges;
    for (int i = 0; i < nBins; ++i)
        edges.push_back(beam.binLow[i]);
    edges.push_back(beam.binHigh.back());

    TH1D *hist = new TH1D(name, title, nBins, edges.data());
    hist->SetDirectory(nullptr);
    for (int i = 0; i < nBins; ++i)
        hist->SetBinContent(i + 1, beam.binScaled[i]);
    return hist;
}

TGraphErrors *MakeAcceptanceGraph(const TriggerAcceptanceInfo& info, const char *name)
{
    TGraphErrors *graph = new TGraphErrors(info.binCenter.size());
    graph->SetName(name);
    for (size_t i = 0; i < info.binCenter.size(); ++i) {
        double xerr = 0.5 * (info.binHigh[i] - info.binLow[i]);
        graph->SetPoint(i, info.binCenter[i], info.binAcc[i]);
        graph->SetPointError(i, xerr, info.binAccErr[i]);
    }
    return graph;
}

TGraphErrors *MakeCrossSectionGraph(const TriggerAcceptanceInfo& info, const char *name)
{
    TGraphErrors *graph = new TGraphErrors(info.binCenter.size());
    graph->SetName(name);
    for (size_t i = 0; i < info.binCenter.size(); ++i) {
        double xerr = 0.5 * (info.binHigh[i] - info.binLow[i]);
        graph->SetPoint(i, info.binCenter[i], info.binCrossSection[i]);
        graph->SetPointError(i, xerr, 0.0);
    }
    return graph;
}

double NormalizeEfficiency(double eff)
{
    if (eff > 1.0)
        return eff / 100.0;
    return eff;
}

double GetLH2TargetArealDensity()
{
    const double mmtocm = 0.1;
    return lh2_density * N_A * lh2_thick * mmtocm * lh2_layer / lh2_W;
}

YieldEstimateInfo EstimateLambdaYield(const TriggerAcceptanceInfo& info,
                                      const BeamSpectrum& beam,
                                      double targetArealDensity,
                                      double daqEff,
                                      double branchingFraction)
{
    YieldEstimateInfo result;
    result.reactionName = info.reactionName;
    result.binLow = info.binLow;
    result.binHigh = info.binHigh;
    result.binCenter = info.binCenter;
    result.totalYield = 0.0;

    const double mbToCm2 = 1.0e-27;
    const double daqEffNorm = NormalizeEfficiency(daqEff);

    for (size_t i = 0; i < info.binCenter.size(); ++i) {
        double nBeam = (i < beam.binScaled.size()) ? beam.binScaled[i] : 0.0;
        double sigmaMb = info.binCrossSection[i];
        double geomEff = info.binAcc[i] / 100.0;
        double binYield = nBeam * targetArealDensity * sigmaMb * mbToCm2 *
                          geomEff * daqEffNorm * branchingFraction;

        result.nBeam.push_back(nBeam);
        result.crossSectionMb.push_back(sigmaMb);
        result.geomEff.push_back(geomEff);
        result.yield.push_back(binYield);
        result.totalYield += binYield;
    }

    return result;
}

YieldEstimateInfo SumYieldInfos(const vector<YieldEstimateInfo>& yieldInfos,
                                const string& name)
{
    YieldEstimateInfo total;
    total.reactionName = name;
    total.totalYield = 0.0;

    if (yieldInfos.empty())
        return total;

    const auto& first = yieldInfos.front();
    total.binLow = first.binLow;
    total.binHigh = first.binHigh;
    total.binCenter = first.binCenter;
    total.nBeam = first.nBeam;
    total.crossSectionMb.assign(first.binCenter.size(), 0.0);
    total.geomEff.assign(first.binCenter.size(), 0.0);
    total.yield.assign(first.binCenter.size(), 0.0);

    for (const auto& info : yieldInfos) {
        for (size_t i = 0; i < total.yield.size() && i < info.yield.size(); ++i) {
            total.yield[i] += info.yield[i];
            total.totalYield += info.yield[i];
        }
    }

    return total;
}

TH1D *MakeYieldHist(const YieldEstimateInfo& info, const char *name, const char *title)
{
    const int nBins = info.binLow.size();
    vector<double> edges;
    for (int i = 0; i < nBins; ++i)
        edges.push_back(info.binLow[i]);
    edges.push_back(info.binHigh.back());

    TH1D *hist = new TH1D(name, title, nBins, edges.data());
    hist->SetDirectory(nullptr);
    for (int i = 0; i < nBins; ++i)
        hist->SetBinContent(i + 1, info.yield[i]);
    return hist;
}

void PrintYieldTable(const vector<YieldEstimateInfo>& yieldInfos,
                     double targetArealDensity,
                     double daqEff,
                     double branchingFraction)
{
    cout << fixed << setprecision(6);
    cout << "target areal density used [/cm2] : " << scientific
         << targetArealDensity << fixed << endl;
    cout << "DAQ efficiency : " << NormalizeEfficiency(daqEff)
         << ", Lambda branching fraction : " << branchingFraction << endl;

    for (const auto& info : yieldInfos) {
        cout << "\n[" << info.reactionName << "] expected Lambda yield" << endl;
        bool isTotal = info.reactionName == "Total";
        cout << "bin [MeV/c]"
             << setw(16) << "Nbeam"
             << setw(14) << "sigma[mb]"
             << setw(14) << "geomEff"
             << setw(16) << "yield" << endl;

        for (size_t i = 0; i < info.binCenter.size(); ++i) {
            cout << setw(4) << static_cast<int>(info.binLow[i]) << "-"
                 << setw(3) << static_cast<int>(info.binHigh[i])
                 << setw(16) << setprecision(2) << info.nBeam[i];
            if (isTotal) {
                cout << setw(14) << "-"
                     << setw(14) << "-";
            } else {
                cout << setw(14) << setprecision(6) << info.crossSectionMb[i]
                     << setw(14) << setprecision(6) << info.geomEff[i];
            }
            cout << setw(16) << setprecision(6) << info.yield[i] << endl;
        }
        cout << "total yield : " << setprecision(6) << info.totalYield << endl;
    }
}



void DrawBeamPage(const BeamSpectrum& beam, TCanvas *c, const string& pdf_file)
{
    c->SetCanvasSize(900, 1100);
    c->Clear();
    c->Divide(1, 1);
    c->cd(1);

    TH1D *hBeam = MakeBeamHist(beam, "h_beam_scaled",
                               "Beam spectrum scaled to TRIG-D;K^{-} momentum [MeV/#it{c}];Counts");
    hBeam->SetLineColor(kMagenta + 2);
    hBeam->SetLineWidth(2);
    hBeam->SetFillColor(0);
    hBeam->Draw("hist");

    TLegend *leg = new TLegend(0.62, 0.78, 0.80, 0.88);
    leg->AddEntry(hBeam, "beam spectrum scaled to TRIG-D", "l");
    leg->Draw();

    c->Print(pdf_file.c_str());
    delete hBeam;
    delete leg;
}

void DrawYieldPage(const YieldEstimateInfo& yieldInfo, TCanvas *c, const string& pdf_file)
{
    c->SetCanvasSize(850, 850);
    c->Clear();
    c->Divide(1, 1);
    c->cd(1);

    TH1D *hYield = MakeYieldHist(yieldInfo,
                                 Form("h_yield_%s", yieldInfo.reactionName.c_str()),
                                 Form("%s;K^{-} momentum [MeV/#it{c}];Expected #Lambda yield",
                                      yieldInfo.reactionName.c_str()));
    hYield->SetLineColor(kViolet + 2);
    hYield->SetLineWidth(2);
    hYield->SetFillColor(kViolet - 9);
    TGaxis::SetMaxDigits(3);
    hYield->Draw("hist");

    TLegend *leg = new TLegend(0.62, 0.78, 0.80, 0.88);
    leg->AddEntry(hYield, "beam #times target #times #sigma #times eff. #times BR", "f");
    leg->Draw();

    c->Print(pdf_file.c_str());
    TGaxis::SetMaxDigits(5);
    delete hYield;
    delete leg;
}

void DrawTriggerAcceptancePdf(const vector<TriggerAcceptanceInfo>& infos,
                              const BeamSpectrum& beam,
                              const vector<YieldEstimateInfo>& yieldInfos,
                              const string& pdf_file)
{
    TCanvas *c = new TCanvas("c_yield", "yield", 900, 1100);
    c->Print((pdf_file + "[").c_str());

    DrawBeamPage(beam, c, pdf_file);

    for (size_t i = 0; i < infos.size(); ++i) {
        const auto& info = infos[i];
        c->SetCanvasSize(900, 1100);
        c->Clear();
        c->Divide(1, 3);

        c->cd(1);
        TH1D *hEntry = MakeCountHist(info, Form("h_entry_%s", info.reactionName.c_str()),
                                     Form("%s;K^{-} momentum [MeV/#it{c}];Counts", info.reactionName.c_str()),
                                     false);
        TH1D *hDen = MakeCountHist(info, Form("h_den_%s", info.reactionName.c_str()),
                                   Form("%s;K^{-} momentum [MeV/#it{c}];Counts", info.reactionName.c_str()),
                                   true);
        hEntry->SetLineColor(kBlack);
        hEntry->SetFillColor(0);
        hEntry->Draw("hist");
        hDen->SetLineColor(kBlue + 1);
        hDen->SetFillColor(kAzure - 9);
        hDen->Draw("hist same");
        TLegend *leg = new TLegend(0.62, 0.72, 0.80, 0.88);
        leg->AddEntry(hEntry, "simulation entries", "l");
        leg->AddEntry(hDen, "physics trigger denominator", "f");
        leg->Draw();

        c->cd(2);
        TGraphErrors *gAcc = MakeAcceptanceGraph(info, Form("g_acc_%s", info.reactionName.c_str()));
        gAcc->SetTitle(Form("%s;K^{-} momentum [MeV/c];geometrical acceptance [%%]", info.reactionName.c_str()));
        gAcc->SetMarkerStyle(20);
        gAcc->SetMarkerColor(kRed + 1);
        gAcc->SetLineColor(kRed + 1);
        gAcc->SetMinimum(0.0);
        gAcc->SetMaximum(100.0);
        gAcc->Draw("APL");

        c->cd(3);
        TGraphErrors *gXsec = MakeCrossSectionGraph(info, Form("g_xsec_%s", info.reactionName.c_str()));
        gXsec->SetTitle(Form("%s;K^{-} momentum [MeV/c];total cross section [mb]", info.reactionName.c_str()));
        gXsec->SetMarkerStyle(21);
        gXsec->SetMarkerColor(kGreen + 2);
        gXsec->SetLineColor(kGreen + 2);
        gXsec->Draw("APL");

        c->Print(pdf_file.c_str());
        delete hEntry;
        delete hDen;
        delete leg;
        delete gAcc;
        delete gXsec;

        if (i < yieldInfos.size())
            DrawYieldPage(yieldInfos[i], c, pdf_file);
    }

    if (yieldInfos.size() > infos.size())
        DrawYieldPage(yieldInfos.back(), c, pdf_file);

    c->Print((pdf_file + "]").c_str());
    delete c;
}

void yield(){
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  
  CombinedScalerInfo scaler = GetCombinedScalerInfo(runnumbers);
  Long64_t nKbeam = scaler.trigD;
  double daqEff = scaler.daqEff;

  if (nKbeam <= 0 || daqEff <= 0.0) {
      cerr << "No valid scaler information for requested runs" << endl;
      return;
  }

  cout << "combined nkbeam : " << nKbeam
       << ", weighted daqEff : " << daqEff << endl;

  string sim_dir = "/home/had/haein/simul-data/e72-lambda-yield";
  string beam_file = sim_dir + "/beam_735_beam.root";
  string xsec_file = "/home/had/haein/work/git/e72-physics/lambda-yield/cs-info/total_cross_section.root";

  // Edit these edges after checking the PDF. Units are MeV/c.
  std::vector<double> momBins;
  for (double mom = 700; mom <= 780; mom += 5) {
    momBins.push_back(mom);
  }

  BeamSpectrum beamSpectrum = GetBeamSpectrum(beam_file, momBins, nKbeam);

  double nKbeamInBins = 0.0;
  for (double count : beamSpectrum.binScaled)
      nKbeamInBins += count;

  cout << "beam spectrum entries : " << beamSpectrum.nEntry
       << ", scale : " << beamSpectrum.scale
       << ", scaled integral : " << nKbeamInBins << endl;

  TFile *fxsec = TFile::Open(xsec_file.c_str(), "READ");
  if (!fxsec || fxsec->IsZombie()) {
      cerr << "Cannot open " << xsec_file << endl;
      if (fxsec)
          fxsec->Close();
      return;
  }

  TGraph *gLambdaEtaXsec = LoadCrossSectionGraph(fxsec, "g_Lambda_eta_pK");
  TGraph *gLambdaPiXsec = LoadCrossSectionGraph(fxsec, "g_Lambda_pi0_pK");
  TGraph *gSigmaPiXsec = LoadCrossSectionGraph(fxsec, "g_Sigma0_pi0_pK");

  TriggerAcceptanceInfo lambdaEtaAcc = GetTriggerAcceptance("Lambda_eta",
                                                            sim_dir + "/e72_lambda_yield_lambda_eta.root",
                                                            "g_Lambda_eta_pK",
                                                            gLambdaEtaXsec,
                                                            momBins);
  TriggerAcceptanceInfo lambdaPiAcc = GetTriggerAcceptance("Lambda_pi0",
                                                           sim_dir + "/e72_lambda_yield_lambda_pi.root",
                                                           "g_Lambda_pi0_pK",
                                                           gLambdaPiXsec,
                                                           momBins);
  TriggerAcceptanceInfo sigmaPiAcc = GetTriggerAcceptance("Sigma0_pi0",
                                                          sim_dir + "/e72_lambda_yield_sigma_pi.root",
                                                          "g_Sigma0_pi0_pK",
                                                          gSigmaPiXsec,
                                                          momBins);

  vector<double> xsecLambdaEta = lambdaEtaAcc.binCrossSection;
  vector<double> xsecLambdaPi = lambdaPiAcc.binCrossSection;
  vector<double> xsecSigmaPi = sigmaPiAcc.binCrossSection;

  Long64_t nBeamLambdaEta = lambdaEtaAcc.nBeam;
  Long64_t nTrigDenLambdaEta = lambdaEtaAcc.nFlag1000 + lambdaEtaAcc.nFlag1101 + lambdaEtaAcc.nFlag1011 + lambdaEtaAcc.nFlag1111;
  Long64_t nTrigNumLambdaEta = lambdaEtaAcc.nFlag1101 + lambdaEtaAcc.nFlag1011 + lambdaEtaAcc.nFlag1111;

  Long64_t nBeamLambdaPi = lambdaPiAcc.nBeam;
  Long64_t nTrigDenLambdaPi = lambdaPiAcc.nFlag1000 + lambdaPiAcc.nFlag1101 + lambdaPiAcc.nFlag1011 + lambdaPiAcc.nFlag1111;
  Long64_t nTrigNumLambdaPi = lambdaPiAcc.nFlag1101 + lambdaPiAcc.nFlag1011 + lambdaPiAcc.nFlag1111;

  Long64_t nBeamSigmaPi = sigmaPiAcc.nBeam;
  Long64_t nTrigDenSigmaPi = sigmaPiAcc.nFlag1000 + sigmaPiAcc.nFlag1101 + sigmaPiAcc.nFlag1011 + sigmaPiAcc.nFlag1111;
  Long64_t nTrigNumSigmaPi = sigmaPiAcc.nFlag1101 + sigmaPiAcc.nFlag1011 + sigmaPiAcc.nFlag1111;

  double geomAccLambdaEta = nTrigDenLambdaEta > 0 ? 100.0 * nTrigNumLambdaEta / nTrigDenLambdaEta : 0.0;
  double geomAccLambdaPi = nTrigDenLambdaPi > 0 ? 100.0 * nTrigNumLambdaPi / nTrigDenLambdaPi : 0.0;
  double geomAccSigmaPi = nTrigDenSigmaPi > 0 ? 100.0 * nTrigNumSigmaPi / nTrigDenSigmaPi : 0.0;

  
  double mmtocm = 0.1;
  double N_lh2 = GetLH2TargetArealDensity();
  double N_gfrp = gfrp_density*N_A*gfrp_thick*mmtocm*gfrp_layer/gfrp_W;
  double N_kapton = kapton_density*N_A*kapton_thick*mmtocm*kapton_layer/kapton_W;
  double N_mylar = mylar_density*N_A*mylar_thick*mmtocm*mylar_layer/mylar_W;

  double N_target_total = N_lh2+N_gfrp+N_kapton+N_mylar;
  cout << "LH2 target : " << N_lh2
       << ", total target used for yield : " << N_target_total << endl;

  vector<TriggerAcceptanceInfo> accInfos = {lambdaEtaAcc, lambdaPiAcc, sigmaPiAcc};
  vector<YieldEstimateInfo> yieldInfos;
  yieldInfos.push_back(EstimateLambdaYield(lambdaEtaAcc, beamSpectrum, N_target_total, daqEff, br_L));
  yieldInfos.push_back(EstimateLambdaYield(lambdaPiAcc, beamSpectrum, N_target_total, daqEff, br_L));
  yieldInfos.push_back(EstimateLambdaYield(sigmaPiAcc, beamSpectrum, N_target_total, daqEff,
                                           br_L * br_Sigma0_to_Lambda));
  yieldInfos.push_back(SumYieldInfos(yieldInfos, "Total"));

  cout << "Sigma0 -> Lambda gamma branching fraction : "
       << br_Sigma0_to_Lambda << endl;
  PrintYieldTable(yieldInfos, N_target_total, daqEff, br_L);
  DrawTriggerAcceptancePdf(accInfos, beamSpectrum, yieldInfos, "simul/yield.pdf");

  fxsec->Close();


  
}
