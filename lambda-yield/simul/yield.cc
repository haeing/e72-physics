#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <cmath>
#include <iomanip>
#include <map>
#include <set>
#include <future>


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
#include <TPaveText.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TParameter.h>
#include <TParticle.h>
#include <TLorentzVector.h>

#include "../../basic-property.hh"

// Dictionaries for the vertex and BEAM branches read by ACLiC.
#ifdef __ROOTCLING__
#pragma link C++ class std::vector<std::vector<int>>+;
#pragma link C++ class std::vector<std::vector<double>>+;
#pragma link C++ class std::vector<TParticle>+;
#endif

using namespace std;

const vector<int> runnumbers = {2447, 2449, 2450, 2451, 2452, 2453, 2454, 2456, 2457, 2458, 2459, 2460, 2462, 2463, 2465, 2466, 2468};
//const vector<int> runnumbers = {2447};

struct ScalerInfo {
    Long64_t trigD;
    double daqEff;
};

struct CombinedScalerInfo {
    vector<int> runnumbers;
    Long64_t trigD;
    double daqEff;
};

// Combine simulation (beam generator 7201 + reaction generator): an accepted
// beam event is followed by one reaction event, whose vertex is sampled along
// the beam path in LH2 and the target frame weighted by the areal density.
// mom_kaon_lab of a reaction event is the K- momentum at the reaction vertex.
struct TriggerAcceptanceInfo {
    string reactionName;
    string fileName;
    string graphName;
    Long64_t nEntry;
    Long64_t nBeam;      // reaction events with the beam trigger bit
    Long64_t nPhysTrig;  // reaction events passing the physics trigger
    double nGen;         // generated reactions with the beam trigger bit (accepted beams)
    double nGenLH2;      // expected number of them in LH2 (path-weight fraction)
    double sigmaSumPass;           // sum of sigma(p_vtx) [mb] over physics-trigger events
    map<int, double> sigmaSumVolume; // same, per reaction vertex volume
    map<int, Long64_t> nPhysTrigVolume;
    // Physics-trigger events passing the TPC bit rule (see CheckTpcPair).
    Long64_t nTpcPair;
    double sigmaSumTpc;
    map<int, double> sigmaSumVolumeTpc;
    map<int, Long64_t> nTpcPairVolume;
    // Of them, no other track in the TPC (ntrack).
    Long64_t nTpcNtrack = 0;
    double sigmaSumTpcNtrack = 0.0;
    map<int, double> sigmaSumVolumeTpcNtrack;
    // Same two with the missing mass cut (MM >= kMissingMassCut).
    Long64_t nTpcMM = 0, nTpcNtrackMM = 0;
    double sigmaSumTpcMM = 0.0, sigmaSumTpcNtrackMM = 0.0;
    map<int, double> sigmaSumVolumeTpcMM, sigmaSumVolumeTpcNtrackMM;
    // Further selections on top of the TPC pair, by name (kExtraSelections).
    map<string, Long64_t> nSel;
    map<string, double> sigmaSumSel;
    map<string, map<int, double>> sigmaSumVolumeSel;
    // Reaction-event counts per trigger condition, each on top of the beam bit.
    Long64_t nReaction;  // saved reaction events
    Long64_t nKvcVeto;   // beam x KVC veto
    Long64_t nHtofMp;    // beam x KVC veto x HTOF Mp >= 2
    Long64_t nHtofFwd;   // beam x KVC veto x HTOF forward
    Long64_t nHtofBoth;  // beam x KVC veto x HTOF Mp >= 2 x HTOF forward
    map<int, Long64_t> nBeamVolume;  // beam-trigger reactions per vertex volume
    vector<double> binLow;
    vector<double> binHigh;
    vector<double> binCenter;
    vector<Long64_t> binEntry;
    vector<Long64_t> binBeam;
    vector<Long64_t> binPhysTrig;
    vector<double> binSigmaSum;
    vector<double> binAcc;
    vector<double> binAccErr;
    vector<double> binCrossSection;
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
    map<int, double> volumeYield;
    // Same, for events whose p and pi- pass the TPC bit rule (also x eff_tpc_L).
    double totalYieldTpc = 0.0;
    map<int, double> volumeYieldTpc;
    double totalYieldNtrack = 0.0; // also ntrack
    map<int, double> volumeYieldNtrack;
    double totalYieldMM = 0.0, totalYieldNtrackMM = 0.0; // also MM >= kMissingMassCut
    map<int, double> volumeYieldMM, volumeYieldNtrackMM;
    map<string, double> totalYieldSel; // kExtraSelections
    map<string, map<int, double>> volumeYieldSel;
};

struct ReactionConfig {
    string reactionName;
    string simFileName;
    string graphName;
    double branchingFraction;
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

const int kLH2VolumeId = 1000;

string VolumeName(int volume_id)
{
    switch (volume_id) {
    case 1000: return "LH2";
    case 1001: return "Kapton";
    case 1002: return "Al";
    case 1003: return "Mylar";
    case 1004: return "GFRP";
    case 1005: return "G10 top";
    case 1006: return "G10 bottom";
    default: return Form("volume %d", volume_id);
    }
}

// TPC bit of trig_flag, recomputed from the hit vectors: the kTPCBit of
// AnaManager.cc reads the TPC TParticle mothers as parent ID and layer, which
// hold 1 or 1000 with TPCPadOn, so the bit is never set. Same rule here: the
// p and pi- of the Lambda each hit >= 5 different layers. The Lambda daughters
// are a p and a pi- track with the same parent track (parentpidtpc is -9999 when
// the Lambda comes from Sigma0); for Sigma+ -> p pi0 it is the p of the Sigma+
// with the pi- of most layers.
// ntrack (as cut_flag of DstTPCLambdaEta): at most kMaxOtherTracks other tracks
// hit >= 5 layers, apart from a primary K- beam; ntrack0: none (compared with the
// data vertex ntrack, as the simulation has no pile-up).
const int kTrigTpcLayerMulti = 5;
const int kMaxOtherTracks = 1;
// Kinematics: PassSpecialKinematics of DstTPCLambdaEta without the mass window,
// from the true Lambda at its production and the pi- momentum at its first TPC hit.
const double kKinPionMomMax = 0.25;                         // GeV/c
const double kKinLambdaMomMin = 0.4, kKinLambdaMomMax = 0.6; // GeV/c
// Selections beyond tpc, ntrack, mm and ntrack_mm, written as <name> and <name>_lh2.
const vector<string> kExtraSelections = {"ntrack0", "ntrack0_mm", "kin", "ntrack_kin", "ntrack0_kin",
                                         "mm_kin", "ntrack_mm_kin", "ntrack0_mm_kin"};
// Missing mass of K- p -> Lambda X as in DstTPCLambdaEta (target p at rest,
// Lambda mass fixed), from the true K- at the vertex (BEAM, MeV/c) and the
// true Lambda (vtx_trackpid 3122, GeV/c). Lambda eta is selected with MM >= 0.5.
const double kMissingMassCut = 0.5; // GeV/c2
double TruthMissingMass(const vector<TParticle>& beam, const vector<vector<int>>& vtxPid,
                        const vector<vector<double>>& vtxPx, const vector<vector<double>>& vtxPy,
                        const vector<vector<double>>& vtxPz)
{
    const double kMassK = 0.493677, kMassP = 0.93827208, kMassLambda = 1.115683;
    const TParticle* kaon = nullptr;
    for (const auto& b : beam)
        if (b.GetPdgCode() == -321) { kaon = &b; break; }
    if (!kaon) return -1.0;
    const TLorentzVector k(kaon->Px() * 1e-3, kaon->Py() * 1e-3, kaon->Pz() * 1e-3,
                           sqrt(kaon->P() * kaon->P() * 1e-6 + kMassK * kMassK));
    for (size_t iv = 0; iv < vtxPid.size() && iv < vtxPx.size() && iv < vtxPy.size() && iv < vtxPz.size(); ++iv)
        for (size_t j = 0; j < vtxPid[iv].size() && j < vtxPx[iv].size() && j < vtxPy[iv].size() && j < vtxPz[iv].size(); ++j) {
            if (vtxPid[iv][j] != 3122) continue;
            TLorentzVector lambda;
            lambda.SetXYZM(vtxPx[iv][j], vtxPy[iv][j], vtxPz[iv][j], kMassLambda);
            const TLorentzVector missing = k + TLorentzVector(0, 0, 0, kMassP) - lambda;
            return missing.M2() > 0 ? missing.M() : 0.0;
        }
    return -1.0;
}
// Momentum [GeV/c] of the true Lambda (vtx_trackpid 3122), or -1.
double TruthLambdaMomentum(const vector<vector<int>>& vtxPid, const vector<vector<double>>& vtxPx,
                           const vector<vector<double>>& vtxPy, const vector<vector<double>>& vtxPz)
{
    for (size_t iv = 0; iv < vtxPid.size() && iv < vtxPx.size() && iv < vtxPy.size() && iv < vtxPz.size(); ++iv)
        for (size_t j = 0; j < vtxPid[iv].size() && j < vtxPx[iv].size() && j < vtxPy[iv].size() && j < vtxPz[iv].size(); ++j)
            if (vtxPid[iv][j] == 3122) return sqrt(vtxPx[iv][j] * vtxPx[iv][j] + vtxPy[iv][j] * vtxPy[iv][j] + vtxPz[iv][j] * vtxPz[iv][j]);
    return -1.0;
}
struct TpcPairResult {
    bool pair = false;
    bool ntrack = false;  // at most kMaxOtherTracks other tracks
    bool ntrack0 = false; // no other track
    double pionMom = -1.0; // GeV/c, at the first TPC hit of the pi-
};
TpcPairResult CheckTpcPair(const vector<int>& trackId, const vector<int>& pid, const vector<int>& parentId,
                           const vector<int>& parentPid, const vector<int>& layer, const vector<double>& mom,
                           const vector<double>& time)
{
    const int kProton = 2212, kPiMinus = -211, kKMinus = -321, kSigmaPlus = 3222;
    struct Track {
        int pid, parentId, parentPid;
        set<int> layers;
        double firstTime, firstMom;
    };
    map<int, Track> tracks;
    const size_t n = min({trackId.size(), pid.size(), parentId.size(), parentPid.size(), layer.size(), mom.size(), time.size()});
    for (size_t i = 0; i < n; ++i) {
        if (layer[i] < 0 || layer[i] >= 32) continue;
        auto it = tracks.find(trackId[i]);
        if (it == tracks.end())
            it = tracks.emplace(trackId[i], Track{pid[i], parentId[i], parentPid[i], {}, time[i], mom[i]}).first;
        it->second.layers.insert(layer[i]);
        if (time[i] < it->second.firstTime) {
            it->second.firstTime = time[i];
            it->second.firstMom = mom[i];
        }
    }
    auto ok = [](const Track& t) { return static_cast<int>(t.layers.size()) >= kTrigTpcLayerMulti; };
    int proton = -1, pion = -1;
    for (const auto& [ip, p] : tracks) {
        if (p.pid != kProton || p.parentId <= 0 || !ok(p)) continue;
        for (const auto& [ipi, pi] : tracks)
            if (pi.pid == kPiMinus && pi.parentId == p.parentId && ok(pi)) { proton = ip; pion = ipi; break; }
        if (proton >= 0) break;
    }
    if (proton < 0) {
        for (const auto& [ip, p] : tracks)
            if (p.pid == kProton && p.parentPid == kSigmaPlus && ok(p)) { proton = ip; break; }
        size_t most = 0;
        for (const auto& [ipi, pi] : tracks)
            if (proton >= 0 && pi.pid == kPiMinus && ok(pi) && pi.layers.size() > most) { pion = ipi; most = pi.layers.size(); }
        if (pion < 0) proton = -1;
    }
    TpcPairResult result;
    result.pair = proton >= 0;
    if (!result.pair) return result;
    result.pionMom = tracks.at(pion).firstMom;
    int nOther = 0;
    for (const auto& [it, t] : tracks) {
        if (it == proton || it == pion || !ok(t)) continue;
        if (t.pid == kKMinus && t.parentId == 0) continue;
        ++nOther;
    }
    result.ntrack = nOther <= kMaxOtherTracks;
    result.ntrack0 = nOther == 0;
    return result;
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
    info.nReaction = 0;
    info.nKvcVeto = 0;
    info.nHtofMp = 0;
    info.nHtofFwd = 0;
    info.nHtofBoth = 0;
    info.nPhysTrig = 0;
    info.nGen = 0.0;
    info.nGenLH2 = 0.0;
    info.sigmaSumPass = 0.0;
    info.nTpcPair = 0;
    info.sigmaSumTpc = 0.0;

    const int nBins = static_cast<int>(mom_bins.size()) - 1;
    for (int i = 0; i < nBins; ++i) {
        double low = mom_bins[i];
        double high = mom_bins[i + 1];
        double center = 0.5 * (low + high);
        info.binLow.push_back(low);
        info.binHigh.push_back(high);
        info.binCenter.push_back(center);
        info.binEntry.push_back(0);
        info.binBeam.push_back(0);
        info.binPhysTrig.push_back(0);
        info.binSigmaSum.push_back(0.0);
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

    TTree *tree = dynamic_cast<TTree*>(file->Get("g4hyptpc"));
    if (!tree) {
        cerr << "g4hyptpc not found in " << sim_file << endl;
        file->Close();
        return info;
    }

    for (const char *name : {"generator", "trig_flag", "mom_kaon_lab", "vtx_volume_id",
                             "reaction_path_volume_id", "reaction_path_weight", "trackidtpc", "pidtpc", "parentidtpc", "parentpidtpc", "layertpc", "pptpc", "timetpc",
                             "BEAM", "vtx_trackpid", "vtx_px", "vtx_py", "vtx_pz"}) {
        if (!tree->GetBranch(name)) {
            cerr << name << " not found in " << sim_file << endl;
            file->Close();
            return info;
        }
    }

    Int_t generator = 0;
    Int_t trigFlag = 0;
    Double_t momKaonLab = 0.0;
    Int_t vtxVolumeId = -1;
    vector<int> *pathVolumeId = nullptr;
    vector<double> *pathWeight = nullptr;
    vector<int> *trackIdTpc = nullptr, *pidTpc = nullptr, *parentIdTpc = nullptr, *parentPidTpc = nullptr, *layerTpc = nullptr;
    vector<double> *ppTpc = nullptr, *timeTpc = nullptr;
    tree->SetBranchStatus("*", 0);
    for (const char *name : {"generator", "trig_flag", "mom_kaon_lab", "vtx_volume_id",
                             "reaction_path_volume_id", "reaction_path_weight", "trackidtpc", "pidtpc", "parentidtpc", "parentpidtpc", "layertpc", "pptpc", "timetpc",
                             "BEAM", "vtx_trackpid", "vtx_px", "vtx_py", "vtx_pz"})
        tree->SetBranchStatus(name, 1);
    tree->SetBranchAddress("generator", &generator);
    tree->SetBranchAddress("trig_flag", &trigFlag);
    tree->SetBranchAddress("mom_kaon_lab", &momKaonLab);
    tree->SetBranchAddress("vtx_volume_id", &vtxVolumeId);
    tree->SetBranchAddress("reaction_path_volume_id", &pathVolumeId);
    tree->SetBranchAddress("reaction_path_weight", &pathWeight);
    vector<TParticle> *beamParticle = nullptr;
    vector<vector<int>> *vtxTrackPid = nullptr;
    vector<vector<double>> *vtxPx = nullptr, *vtxPy = nullptr, *vtxPz = nullptr;
    tree->SetBranchAddress("BEAM", &beamParticle);
    tree->SetBranchAddress("vtx_trackpid", &vtxTrackPid);
    tree->SetBranchAddress("vtx_px", &vtxPx);
    tree->SetBranchAddress("vtx_py", &vtxPy);
    tree->SetBranchAddress("vtx_pz", &vtxPz);
    tree->SetBranchAddress("trackidtpc", &trackIdTpc);
    tree->SetBranchAddress("pidtpc", &pidTpc);
    tree->SetBranchAddress("parentidtpc", &parentIdTpc);
    tree->SetBranchAddress("parentpidtpc", &parentPidTpc);
    tree->SetBranchAddress("layertpc", &layerTpc);
    tree->SetBranchAddress("pptpc", &ppTpc);
    tree->SetBranchAddress("timetpc", &timeTpc);
    // The baskets of these small branches are scattered over a ~80 GB file;
    // prefetch them with a tree cache instead of one read per basket.
    tree->SetCacheSize(200 * 1024 * 1024);
    tree->AddBranchToCache("*", kTRUE);
    tree->StopCacheLearningPhase();

    const Int_t kBeamGenerator = 7201;
    // trig_flag bits (AnaManager.cc): beam, KVC, HTOF Mp >= 2, HTOF forward, TPC
    const Int_t kBeamFlag = 0x10;  // 10000: kaon beam
    const Int_t kKvcSignalFlag = 0x8;  // 01000: KVC signal present
    const Int_t kHtofMpFlag = 0x4;  // 00100: HTOF multiplicity >= 2
    const Int_t kHtofFwdFlag = 0x2;  // 00010: HTOF forward segment
    Long64_t nOutside = 0;

    info.nEntry = tree->GetEntries();
    for (Long64_t i = 0; i < info.nEntry; ++i) {
        tree->GetEntry(i);
        bool isBeam = ((trigFlag & kBeamFlag) != 0);

        if (generator == kBeamGenerator) {
            // Every accepted beam event leads to one reaction (also below threshold,
            // which is not saved); its LH2 share is the LH2 fraction of the path weight.
            if (!isBeam)
                continue;
            double total = 0.0;
            double lh2 = 0.0;
            for (size_t k = 0; k < pathWeight->size() && k < pathVolumeId->size(); ++k) {
                total += pathWeight->at(k);
                if (pathVolumeId->at(k) == kLH2VolumeId)
                    lh2 += pathWeight->at(k);
            }
            if (total <= 0.0)
                continue;
            info.nGen += 1.0;
            info.nGenLH2 += lh2 / total;
            continue;
        }

        int bin = FindMomentumBin(momKaonLab, mom_bins);
        if (bin >= 0)
            ++info.binEntry[bin];

        // Physics trigger: kaon beam x KVC veto x (HTOF Mp >= 2 or HTOF forward).
        bool hasKvcSignal = ((trigFlag & kKvcSignalFlag) != 0);
        bool isPhysTrig = isBeam && !hasKvcSignal &&
            (trigFlag & (kHtofMpFlag | kHtofFwdFlag)) != 0;

        const bool hasMp = (trigFlag & kHtofMpFlag) != 0;
        const bool hasFwd = (trigFlag & kHtofFwdFlag) != 0;
        ++info.nReaction;
        if (isBeam) {
            ++info.nBeam;
            ++info.nBeamVolume[vtxVolumeId];
        }
        if (isBeam && !hasKvcSignal) {
            ++info.nKvcVeto;
            if (hasMp)
                ++info.nHtofMp;
            if (hasFwd)
                ++info.nHtofFwd;
            if (hasMp && hasFwd)
                ++info.nHtofBoth;
        }
        if (isPhysTrig) {
            const double sigma = InterpolateGraph(cross_section_graph, momKaonLab);
            ++info.nPhysTrig;
            info.sigmaSumPass += sigma;
            info.sigmaSumVolume[vtxVolumeId] += sigma;
            ++info.nPhysTrigVolume[vtxVolumeId];
            const TpcPairResult tpcPair = CheckTpcPair(*trackIdTpc, *pidTpc, *parentIdTpc, *parentPidTpc, *layerTpc,
                                                       *ppTpc, *timeTpc);
            const bool passMM = tpcPair.pair &&
                TruthMissingMass(*beamParticle, *vtxTrackPid, *vtxPx, *vtxPy, *vtxPz) >= kMissingMassCut;
            if (tpcPair.pair) {
                const double pLambda = TruthLambdaMomentum(*vtxTrackPid, *vtxPx, *vtxPy, *vtxPz);
                const bool passKin = tpcPair.pionMom >= 0 && tpcPair.pionMom < kKinPionMomMax &&
                    pLambda >= kKinLambdaMomMin && pLambda <= kKinLambdaMomMax;
                const map<string, bool> pass = {
                    {"ntrack0", tpcPair.ntrack0},
                    {"ntrack0_mm", tpcPair.ntrack0 && passMM},
                    {"kin", passKin},
                    {"ntrack_kin", tpcPair.ntrack && passKin},
                    {"ntrack0_kin", tpcPair.ntrack0 && passKin},
                    {"mm_kin", passMM && passKin},
                    {"ntrack_mm_kin", tpcPair.ntrack && passMM && passKin},
                    {"ntrack0_mm_kin", tpcPair.ntrack0 && passMM && passKin}};
                for (const string& name : kExtraSelections) {
                    if (!pass.at(name)) continue;
                    ++info.nSel[name];
                    info.sigmaSumSel[name] += sigma;
                    info.sigmaSumVolumeSel[name][vtxVolumeId] += sigma;
                }
            }
            if (tpcPair.ntrack) {
                ++info.nTpcNtrack;
                info.sigmaSumTpcNtrack += sigma;
                info.sigmaSumVolumeTpcNtrack[vtxVolumeId] += sigma;
            }
            if (passMM) {
                ++info.nTpcMM;
                info.sigmaSumTpcMM += sigma;
                info.sigmaSumVolumeTpcMM[vtxVolumeId] += sigma;
            }
            if (passMM && tpcPair.ntrack) {
                ++info.nTpcNtrackMM;
                info.sigmaSumTpcNtrackMM += sigma;
                info.sigmaSumVolumeTpcNtrackMM[vtxVolumeId] += sigma;
            }
            if (tpcPair.pair) {
                ++info.nTpcPair;
                info.sigmaSumTpc += sigma;
                info.sigmaSumVolumeTpc[vtxVolumeId] += sigma;
                ++info.nTpcPairVolume[vtxVolumeId];
            }
            if (bin < 0)
                ++nOutside;
            else
                info.binSigmaSum[bin] += sigma;
        }

        if (bin < 0)
            continue;

        if (isBeam)
            ++info.binBeam[bin];
        if (isPhysTrig)
            ++info.binPhysTrig[bin];
    }

    if (nOutside > 0)
        cerr << reaction_name << ": " << nOutside
             << " physics-trigger events outside the momentum bins (kept in the total)" << endl;

    for (size_t i = 0; i < info.binAcc.size(); ++i) {
        if (info.binBeam[i] > 0) {
            double acc = static_cast<double>(info.binPhysTrig[i]) / info.binBeam[i];
            info.binAcc[i] = 100.0 * acc;
            info.binAccErr[i] = 100.0 * sqrt(acc * (1.0 - acc) / info.binBeam[i]);
        }
    }

    file->Close();
    return info;
}

TH1D *MakeCountHist(const TriggerAcceptanceInfo& info, const char *name, const char *title, int count_type)
{
    const int nBins = info.binLow.size();
    vector<double> edges;
    for (int i = 0; i < nBins; ++i)
        edges.push_back(info.binLow[i]);
    edges.push_back(info.binHigh.back());

    TH1D *hist = new TH1D(name, title, nBins, edges.data());
    hist->SetDirectory(nullptr);
    for (int i = 0; i < nBins; ++i) {
        Long64_t count = info.binEntry[i];
        if (count_type == 1)
            count = info.binBeam[i];
        else if (count_type == 2)
            count = info.binPhysTrig[i];
        hist->SetBinContent(i + 1, count);
    }
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

void PrintCountWithUnits(const string& label, double value)
{
    ios::fmtflags oldFlags = cout.flags();
    streamsize oldPrecision = cout.precision();

    cout << label << " : "
         << scientific << setprecision(6) << value
         << fixed << setprecision(6) << " (" << value / 1.0e6 << " M)"
         << endl;

    cout.flags(oldFlags);
    cout.precision(oldPrecision);
}

double GetLH2TargetArealDensity()
{
    const double mmtocm = 0.1;
    return lh2_density * N_A * lh2_thick * mmtocm * lh2_layer / lh2_W;
}

// The LH2 target density normalizes the LH2 share of the generated reactions;
// the frame reactions come with the same per-areal-density rate as in the
// simulation, so no frame target is added here:
//   yield = N_K * n_LH2 * sum_pass sigma(p_vtx) / N_gen(LH2) * DAQ * BR * eff
YieldEstimateInfo EstimateLambdaYield(const TriggerAcceptanceInfo& info,
                                      double nKbeam,
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

    if (info.nGenLH2 <= 0.0 || info.nBeam <= 0) {
        cerr << info.reactionName << ": no generated LH2 reaction, yield set to 0" << endl;
        result.nBeam.assign(info.binCenter.size(), 0.0);
        result.crossSectionMb = info.binCrossSection;
        result.geomEff.assign(info.binCenter.size(), 0.0);
        result.yield.assign(info.binCenter.size(), 0.0);
        return result;
    }

    const double mbToCm2 = 1.0e-27;
    const double daqEffNorm = NormalizeEfficiency(daqEff);
    const double scale = nKbeam * targetArealDensity * mbToCm2 / info.nGenLH2 *
        daqEffNorm * branchingFraction * eff_tpc_L;

    for (size_t i = 0; i < info.binCenter.size(); ++i) {
        // Data beams in this vertex-momentum bin, from the simulated spectrum.
        double nBeam = nKbeam * static_cast<double>(info.binBeam[i]) / info.nBeam;
        double geomEff = info.binAcc[i] / 100.0;
        result.nBeam.push_back(nBeam);
        result.crossSectionMb.push_back(info.binCrossSection[i]);
        result.geomEff.push_back(geomEff);
        result.yield.push_back(scale * info.binSigmaSum[i]);
    }
    result.totalYield = scale * info.sigmaSumPass;
    for (const auto& [volume, sigmaSum] : info.sigmaSumVolume)
        result.volumeYield[volume] = scale * sigmaSum;
    result.totalYieldTpc = scale * info.sigmaSumTpc;
    for (const auto& [volume, sigmaSum] : info.sigmaSumVolumeTpc)
        result.volumeYieldTpc[volume] = scale * sigmaSum;
    result.totalYieldNtrack = scale * info.sigmaSumTpcNtrack;
    for (const auto& [volume, sigmaSum] : info.sigmaSumVolumeTpcNtrack)
        result.volumeYieldNtrack[volume] = scale * sigmaSum;
    result.totalYieldMM = scale * info.sigmaSumTpcMM;
    for (const auto& [volume, sigmaSum] : info.sigmaSumVolumeTpcMM)
        result.volumeYieldMM[volume] = scale * sigmaSum;
    result.totalYieldNtrackMM = scale * info.sigmaSumTpcNtrackMM;
    for (const auto& [volume, sigmaSum] : info.sigmaSumVolumeTpcNtrackMM)
        result.volumeYieldNtrackMM[volume] = scale * sigmaSum;
    for (const auto& [name, sigmaSum] : info.sigmaSumSel)
        result.totalYieldSel[name] = scale * sigmaSum;
    for (const auto& [name, volumes] : info.sigmaSumVolumeSel)
        for (const auto& [volume, sigmaSum] : volumes)
            result.volumeYieldSel[name][volume] = scale * sigmaSum;

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
        for (size_t i = 0; i < total.yield.size() && i < info.yield.size(); ++i)
            total.yield[i] += info.yield[i];
        total.totalYield += info.totalYield;
        for (const auto& [volume, yield] : info.volumeYield)
            total.volumeYield[volume] += yield;
        total.totalYieldTpc += info.totalYieldTpc;
        for (const auto& [volume, yield] : info.volumeYieldTpc)
            total.volumeYieldTpc[volume] += yield;
        total.totalYieldNtrack += info.totalYieldNtrack;
        for (const auto& [volume, yield] : info.volumeYieldNtrack)
            total.volumeYieldNtrack[volume] += yield;
        total.totalYieldMM += info.totalYieldMM;
        for (const auto& [volume, yield] : info.volumeYieldMM)
            total.volumeYieldMM[volume] += yield;
        total.totalYieldNtrackMM += info.totalYieldNtrackMM;
        for (const auto& [volume, yield] : info.volumeYieldNtrackMM)
            total.volumeYieldNtrackMM[volume] += yield;
        for (const auto& [name, yield] : info.totalYieldSel)
            total.totalYieldSel[name] += yield;
        for (const auto& [name, volumes] : info.volumeYieldSel)
            for (const auto& [volume, yield] : volumes)
                total.volumeYieldSel[name][volume] += yield;
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
                     double daqEff)
{
    cout << fixed << setprecision(6);
    cout << "target areal density used [/cm2] : " << scientific
         << targetArealDensity << fixed << endl;
    cout << "DAQ efficiency : " << NormalizeEfficiency(daqEff) << endl;

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
        PrintCountWithUnits("total yield", info.totalYield);
        for (const auto& [volume, yield] : info.volumeYield)
            cout << "  vertex in " << setw(10) << VolumeName(volume) << " : "
                 << scientific << setprecision(3) << yield << fixed
                 << " (" << setprecision(1) << (info.totalYield > 0 ? 100.0 * yield / info.totalYield : 0.0)
                 << " %)" << endl;
        PrintCountWithUnits("total yield with p, pi- TPC hits in >= 5 layers", info.totalYieldTpc);
        for (const auto& [volume, yield] : info.volumeYieldTpc)
            cout << "  vertex in " << setw(10) << VolumeName(volume) << " : "
                 << scientific << setprecision(3) << yield << fixed << endl;
    }
}

void PrintTriggerAcceptanceCounts(const vector<TriggerAcceptanceInfo>& infos)
{
    cout << fixed << setprecision(6);

    for (const auto& info : infos) {
        cout << "\n[" << info.reactionName << "] simulation trigger acceptance counts" << endl;
        cout << "bin [MeV/c]"
             << setw(14) << "all"
             << setw(14) << "beam"
             << setw(14) << "physTrig"
             << setw(14) << "acc[%]" << endl;

        Long64_t totalEntry = 0;
        Long64_t totalBeam = 0;
        Long64_t totalPhysTrig = 0;

        for (size_t i = 0; i < info.binCenter.size(); ++i) {
            totalEntry += info.binEntry[i];
            totalBeam += info.binBeam[i];
            totalPhysTrig += info.binPhysTrig[i];

            cout << setw(4) << static_cast<int>(info.binLow[i]) << "-"
                 << setw(3) << static_cast<int>(info.binHigh[i])
                 << setw(14) << info.binEntry[i]
                 << setw(14) << info.binBeam[i]
                 << setw(14) << info.binPhysTrig[i]
                 << setw(14) << info.binAcc[i] << endl;
        }

        double totalAcc = (totalBeam > 0) ? 100.0 * totalPhysTrig / totalBeam : 0.0;
        cout << "total"
             << setw(17) << totalEntry
             << setw(14) << totalBeam
             << setw(14) << totalPhysTrig
             << setw(14) << totalAcc << endl;
        cout << "generated reactions with beam trigger : " << info.nGen
             << ", in LH2 (expected) : " << info.nGenLH2
             << " (" << (info.nGen > 0 ? 100.0 * info.nGenLH2 / info.nGen : 0.0) << " %)" << endl;
        for (const auto& [volume, count] : info.nPhysTrigVolume)
            cout << "  physTrig vertex in " << setw(10) << VolumeName(volume) << " : " << count << endl;
    }
}



void DrawYieldPage(const YieldEstimateInfo& yieldInfo, TCanvas *c, const string& pdf_file)
{
    c->SetCanvasSize(850, 850);
    c->Clear();
    c->Divide(1, 1);
    c->cd(1);

    TH1D *hYield = MakeYieldHist(yieldInfo,
                                 Form("h_yield_%s", yieldInfo.reactionName.c_str()),
                                 Form("%s;K^{-} momentum at vertex [MeV/#it{c}];Expected #Lambda yield",
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

double Percent(double num, double den)
{
    return den > 0 ? 100.0 * num / den : 0.0;
}

// One page with the counts of each trigger condition, the vertex volumes and the yield.
void DrawConditionPage(const TriggerAcceptanceInfo& info, const YieldEstimateInfo& yieldInfo,
                       TCanvas *c, const string& pdf_file)
{
    c->SetCanvasSize(900, 1100);
    c->Clear();
    TPaveText *pave = new TPaveText(0.04, 0.03, 0.96, 0.97, "NDC");
    pave->SetFillColor(0);
    pave->SetBorderSize(0);
    pave->SetTextFont(82);
    pave->SetTextAlign(12);
    pave->SetTextSize(0.019);
    pave->AddText(Form("%s : %s", info.reactionName.c_str(), gSystem->BaseName(info.fileName.c_str())));
    pave->AddText(" ");
    pave->AddText("Trigger condition (reaction events)        count     % of beam");
    const double beam = static_cast<double>(info.nBeam);
    auto row = [&](const char *label, double count) {
        pave->AddText(Form("  %-40s %10.0f   %6.2f", label, count, Percent(count, beam)));
    };
    row("saved reaction events", static_cast<double>(info.nReaction));
    row("kaon beam", beam);
    row("beam x KVC veto", static_cast<double>(info.nKvcVeto));
    row("  x HTOF Mp >= 2", static_cast<double>(info.nHtofMp));
    row("  x HTOF forward", static_cast<double>(info.nHtofFwd));
    row("  x HTOF Mp >= 2 and forward", static_cast<double>(info.nHtofBoth));
    row("physics trigger (Mp >= 2 or forward)", static_cast<double>(info.nPhysTrig));
    row("  x p and pi- TPC hits", static_cast<double>(info.nTpcPair));
    row(Form("    x <= %d other TPC track (ntrack)", kMaxOtherTracks), static_cast<double>(info.nTpcNtrack));
    row("    x MM >= 0.5 GeV/c2", static_cast<double>(info.nTpcMM));
    row("    x ntrack x MM >= 0.5 GeV/c2", static_cast<double>(info.nTpcNtrackMM));
    auto nSel = [&](const char *name) {
        const auto it = info.nSel.find(name);
        return it == info.nSel.end() ? 0.0 : static_cast<double>(it->second);
    };
    row("    x no other TPC track (ntrack0)", nSel("ntrack0"));
    row("    x kinematics", nSel("kin"));
    row("    x ntrack x MM x kinematics", nSel("ntrack_mm_kin"));
    pave->AddText(" ");
    pave->AddText(Form("Generated reactions with beam trigger : %.0f (LH2 share %.0f, %.2f %%)",
                       info.nGen, info.nGenLH2, Percent(info.nGenLH2, info.nGen)));
    pave->AddText(" ");
    pave->AddText("Vertex volume        beam trig.   phys. trig.   acc %      yield");
    for (const auto& [volume, nBeamVol] : info.nBeamVolume) {
        const auto itPhys = info.nPhysTrigVolume.find(volume);
        const Long64_t nPhys = itPhys == info.nPhysTrigVolume.end() ? 0 : itPhys->second;
        const auto itYield = yieldInfo.volumeYield.find(volume);
        const double yield = itYield == yieldInfo.volumeYield.end() ? 0.0 : itYield->second;
        pave->AddText(Form("  %-16s %12lld  %12lld   %6.2f   %10.3e", VolumeName(volume).c_str(),
                           nBeamVol, nPhys, Percent(nPhys, nBeamVol), yield));
    }
    pave->AddText(" ");
    pave->AddText(Form("Expected #Lambda yield : %.3e (eff_tpc_L = %.2f)", yieldInfo.totalYield, eff_tpc_L));
    pave->AddText(Form("  with p, pi- in >= 5 TPC layers : %.3e (also x eff_tpc_L)", yieldInfo.totalYieldTpc));
    pave->AddText(Form("  + ntrack : %.3e,  + MM >= 0.5 : %.3e,  + ntrack + MM : %.3e", yieldInfo.totalYieldNtrack,
                       yieldInfo.totalYieldMM, yieldInfo.totalYieldNtrackMM));
    for (const auto& [volume, yield] : yieldInfo.volumeYieldTpc)
        pave->AddText(Form("    %-16s %10.3e", VolumeName(volume).c_str(), yield));
    pave->Draw();
    c->Print(pdf_file.c_str());
    delete pave;
}

// Summary of all reactions on one page.
void DrawSummaryPage(const vector<TriggerAcceptanceInfo>& infos,
                     const vector<YieldEstimateInfo>& yieldInfos,
                     TCanvas *c, const string& pdf_file)
{
    // The PDF page keeps the size of its first page, so stay in portrait.
    c->SetCanvasSize(900, 1100);
    c->Clear();
    TPaveText *pave = new TPaveText(0.01, 0.30, 0.99, 0.97, "NDC");
    pave->SetFillColor(0);
    pave->SetBorderSize(0);
    pave->SetTextFont(82);
    pave->SetTextAlign(12);
    pave->SetTextSize(0.0115);
    pave->AddText("Simulation Lambda yield summary");
    pave->AddText("trigger columns: % of beam-trigger reactions, each with KVC veto");
    pave->AddText(Form("yields x eff_tpc_L = %.2f;  TPC: p and pi- also hit >= 5 TPC layers", eff_tpc_L));
    pave->AddText(" ");
    pave->AddText("reaction           beam  KVC veto   HTOF Mp  HTOF Fwd phys.trig  TPC hit      yield  LH2 yield  TPC yield    TPC LH2");
    for (size_t i = 0; i < infos.size() && i < yieldInfos.size(); ++i) {
        const auto& info = infos[i];
        const double beam = static_cast<double>(info.nBeam);
        const auto itLH2 = yieldInfos[i].volumeYield.find(kLH2VolumeId);
        const double lh2 = itLH2 == yieldInfos[i].volumeYield.end() ? 0.0 : itLH2->second;
        const auto itTpcLH2 = yieldInfos[i].volumeYieldTpc.find(kLH2VolumeId);
        const double tpcLH2 = itTpcLH2 == yieldInfos[i].volumeYieldTpc.end() ? 0.0 : itTpcLH2->second;
        pave->AddText(Form("%-14s %9lld %9.2f %9.2f %9.2f %9.2f %8.2f %10.3e %10.3e %10.3e %10.3e",
                           info.reactionName.c_str(), info.nBeam, Percent(info.nKvcVeto, beam),
                           Percent(info.nHtofMp, beam), Percent(info.nHtofFwd, beam),
                           Percent(info.nPhysTrig, beam), Percent(info.nTpcPair, beam),
                           yieldInfos[i].totalYield, lh2, yieldInfos[i].totalYieldTpc, tpcLH2));
    }
    if (yieldInfos.size() > infos.size()) {
        const auto& total = yieldInfos.back();
        const auto itLH2 = total.volumeYield.find(kLH2VolumeId);
        const auto itTpcLH2 = total.volumeYieldTpc.find(kLH2VolumeId);
        pave->AddText(Form("%-14s %9s %9s %9s %9s %9s %8s %10.3e %10.3e %10.3e %10.3e", "Total", "", "", "", "", "", "",
                           total.totalYield, itLH2 == total.volumeYield.end() ? 0.0 : itLH2->second,
                           total.totalYieldTpc,
                           itTpcLH2 == total.volumeYieldTpc.end() ? 0.0 : itTpcLH2->second));
    }
    pave->Draw();
    c->Print(pdf_file.c_str());
    delete pave;
}

void DrawTriggerAcceptancePdf(const vector<TriggerAcceptanceInfo>& infos,
                              const vector<YieldEstimateInfo>& yieldInfos,
                              const string& pdf_file)
{
    TCanvas *c = new TCanvas("c_yield", "yield", 900, 1100);
    c->Print((pdf_file + "[").c_str());

    DrawSummaryPage(infos, yieldInfos, c, pdf_file);

    for (size_t i = 0; i < infos.size(); ++i) {
        const auto& info = infos[i];
        if (i < yieldInfos.size())
            DrawConditionPage(info, yieldInfos[i], c, pdf_file);
        c->SetCanvasSize(900, 1100);
        c->Clear();
        c->Divide(1, 3);

        c->cd(1);
        TH1D *hEntry = MakeCountHist(info, Form("h_entry_%s", info.reactionName.c_str()),
                                     Form("%s;K^{-} momentum at vertex [MeV/#it{c}];Counts", info.reactionName.c_str()),
                                     0);
        TH1D *hDen = MakeCountHist(info, Form("h_den_%s", info.reactionName.c_str()),
                                   Form("%s;K^{-} momentum at vertex [MeV/#it{c}];Counts", info.reactionName.c_str()),
                                   1);
        TH1D *hPhysTrig = MakeCountHist(info, Form("h_phys_trig_%s", info.reactionName.c_str()),
                                        Form("%s;K^{-} momentum at vertex [MeV/#it{c}];Counts", info.reactionName.c_str()),
                                        2);
        hEntry->SetLineColor(kBlack);
        hEntry->SetFillColor(0);
        hEntry->Draw("hist");
        hDen->SetLineColor(kBlue + 1);
        hDen->SetFillColor(kAzure - 9);
        hDen->Draw("hist same");
        hPhysTrig->SetLineColor(kRed + 1);
        hPhysTrig->SetLineWidth(2);
        hPhysTrig->SetFillColor(0);
        hPhysTrig->Draw("hist same");
        TLegend *leg = new TLegend(0.58, 0.68, 0.84, 0.88);
        leg->AddEntry(hEntry, "reaction events", "l");
        leg->AddEntry(hDen, "beam trigger total", "f");
        leg->AddEntry(hPhysTrig, "physics trigger passed", "l");
        leg->Draw();

        c->cd(2);
        TGraphErrors *gAcc = MakeAcceptanceGraph(info, Form("g_acc_%s", info.reactionName.c_str()));
        gAcc->SetTitle(Form("%s;K^{-} momentum at vertex [MeV/c];geometrical acceptance [%%]", info.reactionName.c_str()));
        gAcc->SetMarkerStyle(20);
        gAcc->SetMarkerColor(kRed + 1);
        gAcc->SetLineColor(kRed + 1);
        gAcc->SetMinimum(0.0);
        gAcc->SetMaximum(100.0);
        gAcc->Draw("APL");

        c->cd(3);
        TGraphErrors *gXsec = MakeCrossSectionGraph(info, Form("g_xsec_%s", info.reactionName.c_str()));
        gXsec->SetTitle(Form("%s;K^{-} momentum at vertex [MeV/c];total cross section [mb]", info.reactionName.c_str()));
        gXsec->SetMarkerStyle(21);
        gXsec->SetMarkerColor(kGreen + 2);
        gXsec->SetLineColor(kGreen + 2);
        gXsec->Draw("APL");

        c->Print(pdf_file.c_str());
        delete hEntry;
        delete hDen;
        delete hPhysTrig;
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

// Yields per reaction and per-run TRIG-D x DAQ-Eff, for yield-compare.cc to
// scale the simulation to the data runs.
void WriteYieldRoot(const vector<ReactionConfig>& reactions, const vector<YieldEstimateInfo>& yieldInfos,
                    const vector<int>& runs, const string& out_file)
{
    TFile out(out_file.c_str(), "RECREATE");
    TTree runTree("runs", "TRIG-D x DAQ-Eff of the runs used");
    Int_t run = 0;
    Double_t trigDaq = 0.0;
    runTree.Branch("run", &run);
    runTree.Branch("trig_daq", &trigDaq);
    for (int r : runs) {
        ScalerInfo scaler = GetScalerInfo(r);
        if (scaler.trigD <= 0 || scaler.daqEff < 0.0)
            continue;
        run = r;
        trigDaq = scaler.trigD * NormalizeScalerEfficiency(scaler.daqEff);
        runTree.Fill();
    }
    TTree yieldTree("yield", "Expected Lambda yield per reaction (x eff_tpc_L)");
    string reaction;
    Int_t lambdaPeak = 0;
    Double_t total = 0.0, lh2 = 0.0, tpc = 0.0, tpcLH2 = 0.0, ntrack = 0.0, ntrackLH2 = 0.0;
    Double_t mm = 0.0, mmLH2 = 0.0, ntrackMM = 0.0, ntrackMMLH2 = 0.0;
    yieldTree.Branch("reaction", &reaction);
    yieldTree.Branch("lambda_peak", &lambdaPeak); // 0 for Sigma+ -> p pi0 (no Lambda)
    yieldTree.Branch("total", &total);
    yieldTree.Branch("lh2", &lh2);
    yieldTree.Branch("tpc", &tpc);
    yieldTree.Branch("tpc_lh2", &tpcLH2);
    yieldTree.Branch("ntrack", &ntrack);       // tpc and at most kMaxOtherTracks other TPC tracks
    yieldTree.Branch("ntrack_lh2", &ntrackLH2);
    yieldTree.Branch("mm", &mm);               // tpc and MM >= 0.5
    yieldTree.Branch("mm_lh2", &mmLH2);
    yieldTree.Branch("ntrack_mm", &ntrackMM);  // tpc, ntrack and MM >= 0.5
    yieldTree.Branch("ntrack_mm_lh2", &ntrackMMLH2);
    map<string, pair<Double_t, Double_t>> extra; // kExtraSelections: all, LH2
    for (const string& name : kExtraSelections) {
        yieldTree.Branch(name.c_str(), &extra[name].first);
        yieldTree.Branch((name + "_lh2").c_str(), &extra[name].second);
    }
    for (size_t i = 0; i < reactions.size() && i < yieldInfos.size(); ++i) {
        const auto& y = yieldInfos[i];
        reaction = reactions[i].reactionName;
        lambdaPeak = reaction.rfind("Sigmap", 0) == 0 ? 0 : 1;
        total = y.totalYield;
        tpc = y.totalYieldTpc;
        const auto itLH2 = y.volumeYield.find(kLH2VolumeId);
        lh2 = itLH2 == y.volumeYield.end() ? 0.0 : itLH2->second;
        const auto itTpcLH2 = y.volumeYieldTpc.find(kLH2VolumeId);
        tpcLH2 = itTpcLH2 == y.volumeYieldTpc.end() ? 0.0 : itTpcLH2->second;
        auto lh2Of = [](const map<int, double>& m) {
            const auto it = m.find(kLH2VolumeId);
            return it == m.end() ? 0.0 : it->second;
        };
        ntrack = y.totalYieldNtrack;
        ntrackLH2 = lh2Of(y.volumeYieldNtrack);
        mm = y.totalYieldMM;
        mmLH2 = lh2Of(y.volumeYieldMM);
        ntrackMM = y.totalYieldNtrackMM;
        ntrackMMLH2 = lh2Of(y.volumeYieldNtrackMM);
        for (const string& name : kExtraSelections) {
            const auto itTotal = y.totalYieldSel.find(name);
            const auto itVolume = y.volumeYieldSel.find(name);
            extra[name].first = itTotal == y.totalYieldSel.end() ? 0.0 : itTotal->second;
            extra[name].second = itVolume == y.volumeYieldSel.end() ? 0.0 : lh2Of(itVolume->second);
        }
        yieldTree.Fill();
    }
    TParameter<double>("eff_tpc_L", eff_tpc_L).Write();
    TParameter<double>("lh2_thick", lh2_thick).Write();
    out.Write();
    cout << "Wrote " << out_file << endl;
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
  PrintCountWithUnits("number of total K- beam", static_cast<double>(nKbeam));

  string sim_dir = "/home/had/haein/simul-data/e72-k18ana/mom-735";
  string xsec_file = "/home/had/haein/work/git/e72-physics/lambda-yield/cs-info/total_cross_section.root";

  // Edit these edges after checking the PDF. Units are MeV/c.
  std::vector<double> momBins;
  for (double mom = 690; mom <= 790; mom += 2) {
    momBins.push_back(mom);
  }

  TFile *fxsec = TFile::Open(xsec_file.c_str(), "READ");
  if (!fxsec || fxsec->IsZombie()) {
      cerr << "Cannot open " << xsec_file << endl;
      if (fxsec)
          fxsec->Close();
      return;
  }

  vector<ReactionConfig> reactions = {
      // Add a new reaction here: {label, simulation ROOT file, total-cross-section graph, Lambda-visible branching fraction}
      {"Lambda_eta", "e72_tpcon-735_lambda_eta.root", "g_Lambda_eta_pK", br_L_ppi},
      {"Lambda_pi0", "e72_tpcon-735_lambda_pi.root", "g_Lambda_pi0_pK", br_L_ppi},
      {"Sigma0_pi0", "e72_tpcon-735_sigma0_pi0.root", "g_Sigma0_pi0_pK", br_Sigma0_L * br_L_ppi},
      {"Sigmap_pim", "e72_tpcon-735_sigmap_pim.root", "g_Sigmap_pim_pK", br_Sigmap_ppi0},
      {"Lambda_pip_pim", "e72_tpcon-735_lambda_pi_pi.root", "g_Lambda_pip_pim_pK", br_L_ppi},
      {"Sigma0_pip_pim", "e72_tpcon-735_sigma0_pip_pim.root", "g_Sigma0_pip_pim_pK", br_Sigma0_L * br_L_ppi},
      {"Sigmap_pim_pi0", "e72_tpcon-735_sigmap_pim_pi0.root", "g_Sigmap_pim_pi0_pK", br_Sigmap_ppi0}
      
  };

  // The target-frame reactions are in the simulation; only LH2 normalizes.
  double N_lh2 = GetLH2TargetArealDensity();
  PrintCountWithUnits("number of LH2 target", N_lh2);

  vector<TriggerAcceptanceInfo> accInfos;
  vector<YieldEstimateInfo> yieldInfos;

  // Reading is latency bound on GPFS, so the reaction files are read in parallel.
  ROOT::EnableThreadSafety();
  vector<future<TriggerAcceptanceInfo>> accFutures;
  for (const auto& reaction : reactions) {
      TGraph *xsecGraph = LoadCrossSectionGraph(fxsec, reaction.graphName);
      accFutures.push_back(async(launch::async, GetTriggerAcceptance,
                                 reaction.reactionName,
                                 sim_dir + "/" + reaction.simFileName,
                                 reaction.graphName,
                                 xsecGraph,
                                 momBins));
  }

  for (size_t ir = 0; ir < reactions.size(); ++ir) {
      const auto& reaction = reactions[ir];
      TriggerAcceptanceInfo acc = accFutures[ir].get();
      yieldInfos.push_back(EstimateLambdaYield(acc,
                                               static_cast<double>(nKbeam),
                                               N_lh2,
                                               daqEff,
                                               reaction.branchingFraction));
      accInfos.push_back(acc);
  }

  yieldInfos.push_back(SumYieldInfos(yieldInfos, "Total"));
  PrintTriggerAcceptanceCounts(accInfos);
  PrintYieldTable(yieldInfos, N_lh2, daqEff);
  DrawTriggerAcceptancePdf(accInfos, yieldInfos, "yield.pdf");
  WriteYieldRoot(reactions, yieldInfos, runnumbers, "../result/yield-sim.root");
  
  fxsec->Close();


  
}
