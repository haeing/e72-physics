#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>
class TCanvas;
void compare_beam_truth_at_vertex(TCanvas&, const std::string&);

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

namespace {

constexpr double kProtonMassGeV = 0.93827208816;
constexpr double kPionMassGeV = 0.13957039;
constexpr double kLambdaMassGeV = 1.115683;

struct TpcTrack {
  int id;
  int pdg;
  double px;
  double py;
  double pz;
  double absY;
};

double Energy(double px, double py, double pz, double mass)
{
  return std::sqrt(px * px + py * py + pz * pz + mass * mass);
}

double InvariantMass(const TpcTrack &proton, const TpcTrack &pion)
{
  const double e = Energy(proton.px, proton.py, proton.pz, kProtonMassGeV) +
                   Energy(pion.px, pion.py, pion.pz, kPionMassGeV);
  const double px = proton.px + pion.px;
  const double py = proton.py + pion.py;
  const double pz = proton.pz + pion.pz;
  const double mass2 = e * e - px * px - py * py - pz * pz;
  return std::sqrt(std::max(0., mass2));
}

} // namespace

void analyze_sim_lambda_truth()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const std::string inputPath = gSystem->ExpandPathName("~/simul-data/e72-k18ana/mom-735/e72_tpcon_lambda_pi.root");
  std::unique_ptr<TFile> input(TFile::Open(inputPath.c_str(), "READ"));
  if (!input || input->IsZombie()) {
    std::cerr << "Cannot open input file: " << inputPath << std::endl;
    return;
  }
  auto *tree = dynamic_cast<TTree *>(input->Get("g4hyptpc"));
  if (!tree) {
    std::cerr << "TTree 'g4hyptpc' was not found in " << inputPath << std::endl;
    return;
  }

  auto hAll = std::make_unique<TH1D>("h_lambda_p_pim_all",
      ";M_{p#pi^{-}} [GeV/#it{c}^{2}];Candidates", 240, 1.05, 1.25);
  auto hTrue = std::make_unique<TH1D>("h_lambda_p_pim_true",
      ";M_{p#pi^{-}} [GeV/#it{c}^{2}];Truth-matched candidates", 240, 1.05, 1.25);
  auto hMissing = std::make_unique<TH1D>("h_missing_mass_lambda_truth",
      ";M_{X}(K^{-}p #rightarrow #Lambda X) [GeV/#it{c}^{2}];Events", 240, 0.0, 0.70);
  hAll->SetDirectory(nullptr);
  hTrue->SetDirectory(nullptr);
  hMissing->SetDirectory(nullptr);

  TTreeReader reader(tree);
  TTreeReaderValue<Int_t> eventNumber(reader, "evnum");
  TTreeReaderValue<std::vector<Int_t>> hitTrackId(reader, "trackidtpc");
  TTreeReaderValue<std::vector<Int_t>> hitPid(reader, "pidtpc");
  TTreeReaderValue<std::vector<Double_t>> hitY(reader, "y0tpc");
  TTreeReaderValue<std::vector<Double_t>> hitPx(reader, "pxtpc");
  TTreeReaderValue<std::vector<Double_t>> hitPy(reader, "pytpc");
  TTreeReaderValue<std::vector<Double_t>> hitPz(reader, "pztpc");
  TTreeReaderValue<std::vector<Int_t>> vtxType(reader, "vtx_type");
  TTreeReaderValue<std::vector<Int_t>> vtxMotherPdg(reader, "vtx_motherpid");
  TTreeReaderValue<std::vector<std::vector<Int_t>>> vtxTrackId(reader, "vtx_trackid");
  TTreeReaderValue<std::vector<std::vector<Int_t>>> vtxTrackPdg(reader, "vtx_trackpid");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> vtxPx(reader, "vtx_px");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> vtxPy(reader, "vtx_py");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> vtxPz(reader, "vtx_pz");
  TTreeReaderValue<std::vector<TParticle>> beam(reader, "BEAM");

  const Long64_t availableEvents = tree->GetEntries();
  const Long64_t eventsToRead = availableEvents;
  constexpr int kPrintEvents = 10;
  Long64_t readEvents = 0;
  Long64_t candidates = 0;
  Long64_t truthMatched = 0;

  while (readEvents < eventsToRead && reader.Next()) {
    const std::size_t hitCount = std::min({hitTrackId->size(), hitPid->size(), hitY->size(),
                                           hitPx->size(), hitPy->size(), hitPz->size()});
    std::map<int, TpcTrack> tracks;
    for (std::size_t i = 0; i < hitCount; ++i) {
      const int id = hitTrackId->at(i);
      const double absY = std::abs(hitY->at(i));
      auto it = tracks.find(id);
      if (it == tracks.end() || absY < it->second.absY) {
        tracks[id] = {id, hitPid->at(i), hitPx->at(i), hitPy->at(i), hitPz->at(i), absY};
      }
    }

    std::set<std::pair<int, int>> lambdaDaughterPairs;
    const std::size_t vertexCount = std::min({vtxType->size(), vtxMotherPdg->size(),
                                              vtxTrackId->size(), vtxTrackPdg->size(), vtxPx->size(), vtxPy->size(), vtxPz->size()});
    const TParticle *beamKaon = nullptr;
    for (const auto &particle : *beam) if (particle.GetPdgCode() == -321) { beamKaon = &particle; break; }
    if (beamKaon) for (std::size_t ivtx = 0; ivtx < vertexCount; ++ivtx) if (vtxType->at(ivtx) == 0) {
      const auto &pdgs = vtxTrackPdg->at(ivtx); const auto &pxs = vtxPx->at(ivtx);
      const auto &pys = vtxPy->at(ivtx); const auto &pzs = vtxPz->at(ivtx);
      const std::size_t n = std::min({pdgs.size(), pxs.size(), pys.size(), pzs.size()});
      for (std::size_t i = 0; i < n; ++i) if (pdgs[i] == 3122) {
        const double le = std::sqrt(pxs[i]*pxs[i]+pys[i]*pys[i]+pzs[i]*pzs[i]+kLambdaMassGeV*kLambdaMassGeV);
        const double e = beamKaon->Energy()/1000. + kProtonMassGeV - le;
        const double x = beamKaon->Px()/1000. - pxs[i], y = beamKaon->Py()/1000. - pys[i], z = beamKaon->Pz()/1000. - pzs[i];
        hMissing->Fill(std::sqrt(std::max(0., e*e-x*x-y*y-z*z)));
      }
    }
    for (std::size_t ivtx = 0; ivtx < vertexCount; ++ivtx) {
      if (vtxType->at(ivtx) != 1 || vtxMotherPdg->at(ivtx) != 3122)
        continue;
      const auto &ids = vtxTrackId->at(ivtx);
      const auto &pdgs = vtxTrackPdg->at(ivtx);
      int protonId = -1;
      int pionId = -1;
      for (std::size_t idau = 0; idau < std::min(ids.size(), pdgs.size()); ++idau) {
        if (pdgs[idau] == 2212) protonId = ids[idau];
        if (pdgs[idau] == -211) pionId = ids[idau];
      }
      if (protonId >= 0 && pionId >= 0)
        lambdaDaughterPairs.emplace(protonId, pionId);
    }

    std::vector<TpcTrack> protons;
    std::vector<TpcTrack> pions;
    for (const auto &trackEntry : tracks) {
      const TpcTrack &track = trackEntry.second;
      if (track.pdg == 2212) protons.push_back(track);
      if (track.pdg == -211) pions.push_back(track);
    }

    Long64_t eventMatches = 0;
    for (const auto &proton : protons) {
      for (const auto &pion : pions) {
        const double mass = InvariantMass(proton, pion);
        const bool isTrueLambda = lambdaDaughterPairs.count({proton.id, pion.id}) != 0;
        hAll->Fill(mass);
        ++candidates;
        if (isTrueLambda) {
          hTrue->Fill(mass);
          ++truthMatched;
          ++eventMatches;
        }
      }
    }
    if (readEvents < kPrintEvents) {
      std::cout << "Event " << *eventNumber << ": " << protons.size() << " p, "
                << pions.size() << " pi-, " << eventMatches
                << " truth-matched Lambda -> p pi- candidates" << std::endl;
    }
    ++readEvents;
  }

  hAll->SetLineColor(kBlack);
  hTrue->SetLineColor(kRed + 1);
  hTrue->SetLineWidth(2);
  const std::string outputPdf = "sim-lambda-p-pim-truth-mom-735.pdf";
  TCanvas canvas("c_lambda_truth", "Lambda truth matching", 900, 700);
  canvas.Print((outputPdf + "[").c_str());
  hAll->Draw("hist");
  canvas.Print(outputPdf.c_str());
  canvas.Clear();
  hTrue->Draw("hist");
  canvas.Print(outputPdf.c_str());
  canvas.Clear();
  hMissing->Draw("hist");
  canvas.Print(outputPdf.c_str());
  compare_beam_truth_at_vertex(canvas, outputPdf);
  canvas.Print((outputPdf + "]").c_str());
  std::cout << "Read " << readEvents << " / " << availableEvents << " events; built "
            << candidates << " p pi- candidates, of which " << truthMatched
            << " are matched to Lambda -> p pi-." << std::endl;
  std::cout << "Wrote " << outputPdf << std::endl;
}
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

namespace {
struct ReactionVertex { double z, px, py, pz, lx, ly, lz; };
double Magnitude(double x, double y, double z) { return std::sqrt(x*x + y*y + z*z); }
}

// Compare the truth beam momentum at the reaction vertex in the original
// Geant4 lambda-pi file with the closest truth-beam point retained by
// DstTPCGeant4HelixTracking.  Entries are paired by effective_evnum.
void compare_beam_truth_at_vertex(TCanvas &canvas, const std::string &outputPdf)
{
  gROOT->SetBatch(kTRUE); gStyle->SetOptStat(0);
  const std::string dir = gSystem->ExpandPathName("~/simul-data/e72-k18ana/mom-735");
  std::unique_ptr<TFile> g4(TFile::Open((dir + "/e72_tpcon_lambda_pi.root").c_str()));
  std::unique_ptr<TFile> dst(TFile::Open((dir + "/lambda_pi_DstTPCGeant4HelixTracking.root").c_str()));
  if (!g4 || g4->IsZombie() || !dst || dst->IsZombie()) { std::cerr << "Cannot open input ROOT file." << std::endl; return; }
  auto *g4tree = dynamic_cast<TTree *>(g4->Get("g4hyptpc"));
  auto *dsttree = dynamic_cast<TTree *>(dst->Get("tpc"));
  if (!g4tree || !dsttree) { std::cerr << "Required tree was not found." << std::endl; return; }

  std::unordered_map<int, ReactionVertex> vertices;
  TTreeReader g4reader(g4tree);
  TTreeReaderValue<Int_t> effective(g4reader, "effective_evnum");
  TTreeReaderValue<std::vector<TParticle>> beam(g4reader, "BEAM");
  TTreeReaderValue<std::vector<Int_t>> vtxType(g4reader, "vtx_type");
  TTreeReaderValue<std::vector<Double_t>> vtxZ(g4reader, "vtx_z");
  TTreeReaderValue<std::vector<std::vector<Int_t>>> vtxPdg(g4reader, "vtx_trackpid");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> vtxPx(g4reader, "vtx_px"), vtxPy(g4reader, "vtx_py"), vtxPz(g4reader, "vtx_pz");
  while (g4reader.Next()) {
    const TParticle *kaon = nullptr;
    for (const auto &p : *beam) if (p.GetPdgCode() == -321) { kaon = &p; break; }
    if (!kaon || vtxType->size()!=vtxZ->size()) continue;
    for (std::size_t i=0;i<vtxType->size();++i) if(vtxType->at(i)==0) for(std::size_t j=0;j<vtxPdg->at(i).size();++j) if(vtxPdg->at(i).at(j)==3122){ vertices[*effective]={vtxZ->at(i),kaon->Px()/1000.,kaon->Py()/1000.,kaon->Pz()/1000.,vtxPx->at(i).at(j),vtxPy->at(i).at(j),vtxPz->at(i).at(j)}; break; }
  }

  TH1D hDz("h_beam_truth_point_dz", ";z_{beam point}-z_{reaction vertex} [mm];Events", 200, -100., 100.);
  TH1D hDp("h_beam_truth_point_dp", ";|p|_{nearest truth point}-|p|_{vertex} [GeV/#it{c}];Events", 200, -0.05, 0.05);
  TH1D hDpx("h_beam_truth_point_dpx", ";p_{x}^{nearest}-p_{x}^{vertex} [GeV/#it{c}];Events", 200, -0.02, 0.02);
  TTreeReader dstreader(dsttree);
  TTreeReaderValue<Int_t> event(dstreader, "event_number");
  TTreeReaderValue<std::vector<Double_t>> z(dstreader, "beam_z_g4");
  TH1D hBeam("h_mm_beam", ";M_{X} [GeV/#it{c}^{2}];Events",240,0.,0.4), hNear("h_mm_near", ";M_{X} [GeV/#it{c}^{2}];Events",240,0.,0.4);
  TTreeReaderValue<std::vector<Double_t>> px(dstreader, "beam_px_g4");
  TTreeReaderValue<std::vector<Double_t>> py(dstreader, "beam_py_g4");
  TTreeReaderValue<std::vector<Double_t>> pz(dstreader, "beam_pz_g4");
  Long64_t matched = 0;
  while (dstreader.Next()) {
    const auto found = vertices.find(*event); if (found == vertices.end()) continue;
    if (z->size() != px->size() || z->size() != py->size() || z->size() != pz->size()) continue;
    std::size_t best = z->size(); double bestDz = 1.e99;
    for (std::size_t i = 0; i < z->size(); ++i) if (std::isfinite(z->at(i)) && std::abs(z->at(i)-found->second.z) < bestDz) { bestDz = std::abs(z->at(i)-found->second.z); best = i; }
    if (best == z->size()) continue;
    const auto &v = found->second;
    hDz.Fill(z->at(best)-v.z); hDp.Fill(Magnitude(px->at(best),py->at(best),pz->at(best))-Magnitude(v.px,v.py,v.pz)); hDpx.Fill(px->at(best)-v.px); ++matched;
    auto mm=[&](double a,double b,double c){double pb=Magnitude(a,b,c),pl=Magnitude(v.lx,v.ly,v.lz),e=std::sqrt(pb*pb+.493677*.493677)+.93827208-std::sqrt(pl*pl+1.115683*1.115683); return std::sqrt(std::max(0.,e*e-(a-v.lx)*(a-v.lx)-(b-v.ly)*(b-v.ly)-(c-v.lz)*(c-v.lz)));}; hBeam.Fill(mm(v.px,v.py,v.pz)); hNear.Fill(mm(px->at(best),py->at(best),pz->at(best)));
  }
  canvas.Clear(); hBeam.Draw(); canvas.Print(outputPdf.c_str()); canvas.Clear(); hNear.Draw(); canvas.Print(outputPdf.c_str());
  //c.Print((outputPdf+"[").c_str()); hBeam.Draw(); c.Print(outputPdf.c_str()); c.Clear(); hNear.Draw(); c.Print(outputPdf.c_str()); c.Print((outputPdf+"]").c_str());
}
