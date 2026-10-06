// PlotReactionBeamMomentumVsTargetThickness.cc

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TChain.h>
#include <TDatime.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TParticle.h>
#include <TTree.h>
#include <TStyle.h>

void PlotReactionBeamMomentumVsTargetThickness()
{
  // ================= User-configurable parameters =================
  const std::vector<std::string> inputFiles = {
    "/home/had/haein/simul-data/e72-k18ana/mom-735/"
    "e72_tpcon-735_beam_735.root"
  };
  const std::string outputPdf =
    "/home/had/haein/work/git/e72-physics/target/"
    "ReactionBeamMomentumVsTargetThickness_e72_tpcon-735_beam_735.pdf";

  // true: include all recorded target/frame volumes.
  // false: use only LH2 target volume copy number 1000.
  const bool includeTargetFrame = true;
  const int targetVolumeId = 1000;
  const int firstGenerator = 7201;
  const int reactionGenerator = 7202;

  const int nMomentumBins = 240;
  const double momentumMin = 0.7;   // GeV/c
  const double momentumMax = 0.8;   // GeV/c
  const int nThicknessBins = 240;
  const double thicknessMin = 0.0; // g/cm^2
  const double thicknessMax = 5.0; // g/cm^2
  const double pathLengthMin = 0.0; // mm
  const double pathLengthMax = 100.0; // mm
  // ================================================================

  TChain chain("g4hyptpc");
  std::ostringstream inputSummary;
  int nFilesAdded = 0;
  for (const auto& fileName : inputFiles) {
    TFile input(fileName.c_str(), "READ");
    if (input.IsZombie()) {
      Warning("PlotReactionBeamMomentumVsTargetThickness",
              "Cannot open input file: %s", fileName.c_str());
      continue;
    }
    auto* tree = dynamic_cast<TTree*>(input.Get("g4hyptpc"));
    if (!tree) {
      Warning("PlotReactionBeamMomentumVsTargetThickness",
              "Tree g4hyptpc not found in: %s", fileName.c_str());
      continue;
    }
    chain.Add(fileName.c_str());
    inputSummary << fileName << "\n";
    ++nFilesAdded;
  }
  if (nFilesAdded == 0) {
    Error("PlotReactionBeamMomentumVsTargetThickness",
          "No usable input ROOT file was found.");
    return;
  }

  std::vector<TParticle>* beam = nullptr;
  std::vector<int>* volumeId = nullptr;
  std::vector<double>* pathWeight = nullptr;
  std::vector<double>* pathLength = nullptr;
  int generator = -1;
  if (chain.SetBranchAddress("generator", &generator) < 0 ||
      chain.SetBranchAddress("BEAM", &beam) < 0 ||
      chain.SetBranchAddress("reaction_path_volume_id", &volumeId) < 0 ||
      chain.SetBranchAddress("reaction_path_weight", &pathWeight) < 0 ||
      chain.SetBranchAddress("reaction_path_length", &pathLength) < 0) {
    Error("PlotReactionBeamMomentumVsTargetThickness",
          "Required branch is missing.");
    return;
  }

  auto* h2 = new TH2D(
    "hMomentumVsThickness",
    "Reaction-vertex beam momentum vs target thickness;"
    "beam momentum at reaction vertex [GeV/c];"
    "density-weighted path length [g/cm^{2}]",
    nMomentumBins, momentumMin, momentumMax,
    nThicknessBins, thicknessMin, thicknessMax);
  h2->SetDirectory(nullptr);
  h2->SetStats(false);
  h2->SetContour(100);

  auto* h2PathLength = new TH2D(
    "hMomentumVsPathLength",
    "Reaction-vertex beam momentum vs geometric path length;"
    "beam momentum at reaction vertex [GeV/c];"
    "geometric path length [mm]",
    nMomentumBins, momentumMin, momentumMax,
    nThicknessBins, pathLengthMin, pathLengthMax);
  h2PathLength->SetDirectory(nullptr);
  h2PathLength->SetStats(false);
  h2PathLength->SetContour(100);

  auto* h2Initial = new TH2D(
    "hInitialMomentumVsThickness",
    "Initial beam momentum vs target thickness;"
    "initial beam momentum [GeV/c];"
    "density-weighted path length [g/cm^{2}]",
    nMomentumBins, momentumMin, momentumMax,
    nThicknessBins, thicknessMin, thicknessMax);
  h2Initial->SetDirectory(nullptr);
  h2Initial->SetStats(false);
  h2Initial->SetContour(100);

  auto* h2InitialPathLength = new TH2D(
    "hInitialMomentumVsPathLength",
    "Initial beam momentum vs geometric path length;"
    "initial beam momentum [GeV/c];"
    "geometric path length [mm]",
    nMomentumBins, momentumMin, momentumMax,
    nThicknessBins, pathLengthMin, pathLengthMax);
  h2InitialPathLength->SetDirectory(nullptr);
  h2InitialPathLength->SetStats(false);
  h2InitialPathLength->SetContour(100);

  auto* h2InitialVsTgtTotalLength = new TH2D(
    "hInitialMomentumVsTgtTotalLength",
    "Initial beam momentum vs total TGT path length;"
    "initial beam momentum [GeV/c];"
    "total TGT path length [mm]",
    nMomentumBins, momentumMin, momentumMax,
    nThicknessBins, pathLengthMin, pathLengthMax);
  h2InitialVsTgtTotalLength->SetDirectory(nullptr);
  h2InitialVsTgtTotalLength->SetStats(false);
  h2InitialVsTgtTotalLength->SetContour(100);

  //const Long64_t nEntries = chain.GetEntries();
  const Long64_t nEntries = 10000;
  Long64_t nFilled = 0;
  Long64_t nReaction = 0;
  Long64_t nMissingInitialMomentum = 0;
  Long64_t nPathMismatch = 0;
  std::size_t maxSegments = 0;
  double maxThickness = 0.0;
  double maxGeometricLength = 0.0;
  double initialMomentum = -1.0;
  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    chain.GetEntry(entry);
    if (!beam || beam->empty() || !volumeId || !pathWeight || !pathLength)
      continue;

    const double beamMomentum = beam->front().P() / 1000.0; // GeV/c
    if (generator == firstGenerator) {
      initialMomentum = beamMomentum;
      continue;
    }
    if (generator != reactionGenerator) continue;
    ++nReaction;
    if (initialMomentum < 0.0) {
      ++nMissingInitialMomentum;
      continue;
    }
    const std::size_t nSegments =
      std::min({volumeId->size(), pathWeight->size(), pathLength->size()});
    if (volumeId->size() != pathWeight->size() ||
        volumeId->size() != pathLength->size()) {
      ++nPathMismatch;
    }
    double thickness = 0.0;
    double geometricLength = 0.0;
    double tgtGeometricLength = 0.0;
    for (std::size_t i = 0; i < nSegments; ++i) {
      if (includeTargetFrame || (*volumeId)[i] == targetVolumeId) {
        thickness += (*pathWeight)[i];
        geometricLength += (*pathLength)[i];
      }
      if ((*volumeId)[i] == targetVolumeId)
        tgtGeometricLength += (*pathLength)[i];
    }
    maxSegments = std::max(maxSegments, nSegments);
    maxThickness = std::max(maxThickness, thickness);
    maxGeometricLength = std::max(maxGeometricLength, geometricLength);
    h2->Fill(beamMomentum, thickness);
    h2PathLength->Fill(beamMomentum, geometricLength);
    h2Initial->Fill(initialMomentum, thickness);
    h2InitialPathLength->Fill(initialMomentum, geometricLength);
    h2InitialVsTgtTotalLength->Fill(initialMomentum, tgtGeometricLength);
    ++nFilled;
  }

  gStyle->SetOptStat(0);
  TCanvas canvas("canvas", "canvas", 1000, 800);
  canvas.SetRightMargin(0.14);
  canvas.Print((outputPdf + "[").c_str());

  // Information/title page.
  canvas.Clear();
  TLatex text;
  text.SetNDC(true);
  text.SetTextFont(42);
  text.SetTextSize(0.035);
  text.DrawLatex(0.08, 0.90,
                 "Reaction beam momentum vs target thickness");
  text.SetTextSize(0.026);
  TDatime now;
  std::ostringstream dateLine;
  dateLine << "Executed: " << now.AsSQLString();
  text.DrawLatex(0.08, 0.84, dateLine.str().c_str());
  text.DrawLatex(0.08, 0.79,
                 "Macro: PlotReactionBeamMomentumVsTargetThickness.cc");
  text.DrawLatex(0.08, 0.74,
                 includeTargetFrame
                   ? "Thickness: all recorded target/frame volumes"
                   : "Thickness: LH2 target only (volume ID 1000)");
  text.DrawLatex(0.08, 0.69, "Thickness unit: g/cm^{2}");
  text.DrawLatex(0.08, 0.64,
                 "Reaction momentum: BEAM branch on reaction event");
  text.DrawLatex(0.08, 0.60,
                 "Initial momentum: preceding first-generator beam event");

  std::ostringstream countLine;
  countLine << "Entries: " << nEntries << "   Reaction: " << nReaction
            << "   Filled: " << nFilled;
  text.DrawLatex(0.08, 0.54, countLine.str().c_str());
  std::ostringstream diagnosticLine;
  diagnosticLine << "Missing initial: " << nMissingInitialMomentum
                 << "   Path-size mismatch: " << nPathMismatch;
  text.DrawLatex(0.08, 0.50, diagnosticLine.str().c_str());
  std::ostringstream rangeLine;
  rangeLine << "Max path segments: " << maxSegments
            << "   Max thickness: " << maxThickness
            << " g/cm^{2}   Max length: " << maxGeometricLength << " mm";
  text.DrawLatex(0.08, 0.46, rangeLine.str().c_str());
  text.DrawLatex(0.08, 0.40, "Input files:");
  text.SetTextSize(0.021);
  double y = 0.36;
  std::istringstream files(inputSummary.str());
  std::string line;
  while (std::getline(files, line)) {
    text.DrawLatex(0.11, y, line.c_str());
    y -= 0.032;
  }
  canvas.Modified();
  canvas.Update();
  canvas.Print(outputPdf.c_str());

  // One plot per page.
  canvas.Clear();
  h2->Draw("COLZ");
  canvas.Modified();
  canvas.Update();
  canvas.Print(outputPdf.c_str());
  canvas.Clear();
  h2PathLength->Draw("COLZ");
  canvas.Modified();
  canvas.Update();
  canvas.Print(outputPdf.c_str());
  canvas.Clear();
  h2Initial->Draw("COLZ");
  canvas.Modified();
  canvas.Update();
  canvas.Print(outputPdf.c_str());
  canvas.Clear();
  h2InitialPathLength->Draw("COLZ");
  canvas.Modified();
  canvas.Update();
  canvas.Print(outputPdf.c_str());
  canvas.Clear();
  h2InitialVsTgtTotalLength->Draw("COLZ");
  canvas.Modified();
  canvas.Update();
  canvas.Print(outputPdf.c_str());
  canvas.Print((outputPdf + "]").c_str());

  delete h2;
  delete h2PathLength;
  delete h2Initial;
  delete h2InitialPathLength;
  delete h2InitialVsTgtTotalLength;
  std::cout << "Wrote " << outputPdf << " using " << nFilled
            << " entries." << std::endl;
}
