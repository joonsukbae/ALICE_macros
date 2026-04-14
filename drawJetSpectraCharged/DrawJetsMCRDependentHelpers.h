#ifndef DRAWJETSMCRDEPENDENTHELPERS_H
#define DRAWJETSMCRDEPENDENTHELPERS_H

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TSystem.h"
#include "Rtypes.h"
#include "Filipad2.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldSvd.h"
#include "TSVDUnfold_local.h"
#include "TMatrixD.h"
#include "TLatex.h"
#include "TLine.h"
#include "TF1.h"
#include "TStyle.h"
#include "TList.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <cmath>
#include <algorithm>
using namespace std;

// ============================================
// Logging: tee stderr to file for debugging
// ============================================
class TeeBuf : public std::streambuf {
  std::streambuf* fSb1;  // original stderr
  std::streambuf* fSb2;  // log file
public:
  TeeBuf(std::streambuf* sb1, std::streambuf* sb2) : fSb1(sb1), fSb2(sb2) {}
  int overflow(int c) override {
    if (c == EOF) return !EOF;
    if (fSb1->sputc(c) == EOF) return EOF;
    if (fSb2->sputc(c) == EOF) return EOF;
    return c;
  }
  int sync() override {
    fSb1->pubsync();
    fSb2->pubsync();
    return 0;
  }
};

struct LogFileState {
  std::ofstream file;
  TeeBuf* tee = nullptr;
  std::streambuf* origBuf = nullptr;
  bool active = false;
};

inline LogFileState& GetLogState() {
  static LogFileState s;
  return s;
}

// Call at the start of a config set run. Creates/appends analysis.log in outputDir.
inline void OpenLogFile(const char* outputDir) {
  auto& s = GetLogState();
  if (s.active) return;
  gSystem->mkdir(outputDir, true);
  s.file.open(Form("%s/analysis.log", outputDir), std::ios::app);
  if (s.file.is_open()) {
    s.origBuf = std::cerr.rdbuf();
    s.tee = new TeeBuf(s.origBuf, s.file.rdbuf());
    std::cerr.rdbuf(s.tee);
    s.active = true;
    // Timestamp
    time_t now = time(nullptr);
    char timeBuf[64];
    strftime(timeBuf, sizeof(timeBuf), "%Y-%m-%d %H:%M:%S", localtime(&now));
    std::cerr << "======== Log opened: " << timeBuf << " ========" << std::endl;
  }
}

// Call at the end of a config set run. Restores stderr.
inline void CloseLogFile() {
  auto& s = GetLogState();
  if (!s.active) return;
  std::cerr << "======== Log closed ========" << std::endl;
  std::cerr.rdbuf(s.origBuf);
  delete s.tee;
  s.tee = nullptr;
  s.file.close();
  s.active = false;
}

// Forward declarations
class TFile;
class TH1;
class TH2;
class TH3;
class TLegend;
class TPad;
class TCanvas;

// ============================================
// Configuration: R-dependent data/MC inputs
// ============================================

struct RConfig {
  Double_t R;
  TString dataRunNumber;
  TString mcRunNumber;
  TString dataDir;
  TString mcDir;
  TString label;
  bool isJJ;  // true = JJ MC (requires weighted histograms), false = MB MC

  RConfig() : R(0.0), isJJ(false) {}

  RConfig(Double_t r, const char* dataRun, const char* mcRun,
           const char* dDir, const char* mDir, const char* lbl, bool jj = false)
    : R(r), dataRunNumber(dataRun), mcRunNumber(mcRun),
      dataDir(dDir), mcDir(mDir), label(lbl), isJJ(jj) {}
  
  TString GetDataFile() const {
    static const TString baseDir = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/";
    static const TString fileSuffix = "_AnalysisResults.root";
    return baseDir + dataRunNumber + fileSuffix;
  }
  TString GetMcFile() const {
    static const TString baseDir = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/";
    static const TString fileSuffix = "_AnalysisResults.root";
    return baseDir + mcRunNumber + fileSuffix;
  }
};

// ============================================
// Global variables as functions (to avoid SIOF)
// ============================================

// Cross section mode: NonUE (inclusive) or UE (UE-subtracted)
enum EXsecMode {
  kNonUE = 0,
  kUE = 1
};

// Global mode getter/setter
inline EXsecMode& GetCurrentMode() {
  static EXsecMode gMode = kNonUE;
  return gMode;
}

inline void SetCurrentMode(EXsecMode mode) {
  GetCurrentMode() = mode;
}

// Cross-section normalization mode
enum EXsecNormMode {
  kNormSoft = 0,    // Framework luminosity (per-run CCDB σ_vis via hLumiTVX)
  kNormDirect = 1   // Hardcoded σ_vis (53.0 mb for 2022, 53.4 mb for 2023)
};

inline EXsecNormMode& GetXsecNormMode() {
  static EXsecNormMode mode = kNormSoft;
  return mode;
}
inline void SetXsecNormMode(EXsecNormMode m) { GetXsecNormMode() = m; }

// Unfolding method selection
enum EUnfoldMethod {
  kUnfoldSVD = 0,
  kUnfoldBayes = 1
};

inline EUnfoldMethod& GetUnfoldMethod() {
  static EUnfoldMethod method = kUnfoldSVD;
  // static EUnfoldMethod method = kUnfoldBayes;
  return method;
}
inline void SetUnfoldMethod(EUnfoldMethod m) { GetUnfoldMethod() = m; }

// ============================================
// Multi-set configuration
// ============================================

struct RConfigSet {
  TString name;                     // e.g. "LHC22o_LHC25a2b_570220"
  TString outputDir;                // e.g. "../plots/RDependentComparison/LHC22o..."
  EXsecMode mode;                   // kNonUE or kUE
  std::vector<RConfig> configs;     // the 7 R-value configs
  double jetEtaCut;                 // <0: R-dependent |eta|<0.9-R (default)
                                    // >=0: fixed |eta|<jetEtaCut for all R
  double productionEtaCut;          // Eta cut used in O2Physics task/jet finder
                                    // <0: R-dependent |eta|<0.9-R (standard)
                                    // >=0: fixed |eta|<productionEtaCut
                                    // Used for histograms without eta axis (constituent pT,
                                    // track pT, 1D truth) that can't be re-projected.
  bool runBothModes;                // If true, auto-generate both kNonUE and kUE sets
                                    // s.mode is used as the "primary" mode; the other is
                                    // created with "_UEsub" appended (or removed) in name/outputDir.

  RConfigSet() : mode(kNonUE), jetEtaCut(-1), productionEtaCut(-1), runBothModes(false) {}
};

inline RConfigSet& GetCurrentConfigSet() {
  static RConfigSet gCurrentSet;
  return gCurrentSet;
}

inline void SetCurrentConfigSet(const RConfigSet& set) {
  GetCurrentConfigSet() = set;
  SetCurrentMode(set.mode);
}

// Get jet eta acceptance for a given R, using current config set's jetEtaCut
// jetEtaCut < 0 (default): R-dependent |eta| < 0.9 - R
// jetEtaCut >= 0: fixed |eta| < jetEtaCut for all R
inline Double_t GetJetEtaMax(Double_t R) {
  double cut = GetCurrentConfigSet().jetEtaCut;
  return (cut >= 0) ? cut : (0.9 - R);
}

inline Double_t GetDeltaEta(Double_t R) {
  return 2.0 * TMath::Abs(GetJetEtaMax(R));
}

// Production-level eta: what the O2Physics task/jet finder actually used.
// Histograms without an eta axis (constituent pT, track pT, 1D truth, response matrix)
// were filled with this cut and cannot be re-projected.
inline Double_t GetProductionEtaMax(Double_t R) {
  double cut = GetCurrentConfigSet().productionEtaCut;
  return (cut >= 0) ? cut : (0.9 - R);
}

inline Double_t GetProductionDeltaEta(Double_t R) {
  return 2.0 * TMath::Abs(GetProductionEtaMax(R));
}

// Load reco jet pT with analysis-level eta cut applied.
// If jetEtaCut >= 0, projects from 3D histogram (h3_jet_pt_jet_eta_jet_phi) with eta cut.
// Falls back to 1D histogram (h_jet_pt) if 3D not available or jetEtaCut < 0.
// Works for both data and MC reco files. UE-sub uses the rho-area-subtracted version.
inline TH1* LoadJetPtWithEtaCut(TFile* file, const char* dir, bool useUEsub, const char* cloneName) {
  double cut = GetCurrentConfigSet().jetEtaCut;
  // Try 3D first when a custom eta cut is set
  if (cut >= 0) {
    const char* h3Name = useUEsub ? "h3_jet_pt_jet_eta_jet_phi_rhoareasubtracted"
                                  : "h3_jet_pt_jet_eta_jet_phi";
    TH3* h3 = (TH3*)file->Get(Form("%s/%s", dir, h3Name));
    if (h3) {
      int etaBinLo = h3->GetYaxis()->FindBin(-cut + 1e-6);
      int etaBinHi = h3->GetYaxis()->FindBin( cut - 1e-6);
      TH1* hProj = h3->ProjectionX(cloneName, etaBinLo, etaBinHi, 1, h3->GetNbinsZ());
      if (hProj) {
        hProj->SetDirectory(nullptr);
        std::cerr << "[Info] LoadJetPtWithEtaCut: projected 3D with |eta| < " << cut
                  << " [bins " << etaBinLo << "-" << etaBinHi << "] -> " << cloneName << std::endl;
        return hProj;
      }
    }
    std::cerr << "[Warning] LoadJetPtWithEtaCut: 3D histogram not found, falling back to 1D" << std::endl;
  }
  // Fallback: 1D histogram
  const char* h1Name = useUEsub ? "h_jet_pt_rhoareasubtracted" : "h_jet_pt";
  TH1* h1 = (TH1*)file->Get(Form("%s/%s", dir, h1Name));
  if (h1) {
    TH1* hClone = (TH1*)h1->Clone(cloneName);
    hClone->SetDirectory(nullptr);
    return hClone;
  }
  return nullptr;
}

// Auto-detect dataset year from config set name ("LHC22..." → 2022, "LHC23..." → 2023)
inline int DetectDatasetYear(const TString& configSetName) {
  if (configSetName.Contains("LHC23") || configSetName.Contains("2023")) return 2023;
  if (configSetName.Contains("LHC24") || configSetName.Contains("2024")) return 2024;
  return 2022;  // default
}

inline bool GetNORMEVENTS() { return true; }
inline const char* GetEventObj() { return "h_collisions"; }
inline const char* GetEventWObj() { return "h_collisions_weighted"; }
inline bool GetDRAWPLOTS() { return true; }
inline int& GetNN() { static int nn = 0; return nn; }
inline bool GetREBINON() { return true; }
inline double GetPlotPtMin() { return 5.0; }
inline double GetPlotPtMax() { return 140.0; }

// Bin arrays
inline const Double_t* GetTrackptbin() {
  static const Double_t Trackptbin[21] = {-0.5, 1.5, 3.5, 5.5, 7.5, 9.5, 14.5, 19.5, 24.5, 29.5, 39.5, 49.5, 59.5, 69.5, 79.5, 89.5, 99.5, 119.5, 139.5, 169.5, 199.5};
  return Trackptbin;
}
// Master bin edges — superset covering all possible histogram ranges.
// Clipped to the actual histogram axis range by BuildPtBinning().
// Includes negative bins for UE-subtracted jets (pT_reco = pT_raw - rho*A can be < 0).
inline const std::vector<Double_t>& GetMasterRecoBins() {
  static const std::vector<Double_t> v = {
    -200, -100, -50, -20, -10, -5, -2, -1,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20,
    25, 30, 40, 50, 60, 70, 85, 100, 140, 200, 300};
  return v;
}
inline const std::vector<Double_t>& GetMasterTruthBins() {
  static const std::vector<Double_t> v = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20,
    25, 30, 40, 50, 60, 70, 85, 100, 140, 200, 300};
  return v;
}

inline int GetNTrackptbin() { return 20; }

// Backward-compatible fixed bin arrays (used by DrawDataCrossSection.C etc.)
inline const Double_t* GetPtbin() {
  static const Double_t ptbin[21] = {5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
  return ptbin;
}
inline const Double_t* GetPtbinGen() {
  static const Double_t ptbinGen[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
  return ptbinGen;
}
inline int GetNptBins() { return 20; }
inline int GetNptBinsGen() { return 25; }

// Dynamic pT binning: bins are clipped from master array to match the histogram range.
struct PtBinning {
  std::vector<Double_t> edges; // owned storage
  const Double_t* bins;        // pointer to edges[0] for TH1::Rebin compatibility
  int nBins;                   // number of bins = edges.size() - 1

  PtBinning() : bins(nullptr), nBins(0) {}
  PtBinning(const std::vector<Double_t>& e) : edges(e), bins(edges.data()), nBins((int)edges.size() - 1) {}
  PtBinning(const PtBinning& o) : edges(o.edges), bins(edges.data()), nBins(o.nBins) {}
  PtBinning& operator=(const PtBinning& o) {
    if (this != &o) { edges = o.edges; bins = edges.data(); nBins = o.nBins; }
    return *this;
  }
};

// Build a PtBinning by selecting master edges that fall within [axisMin, axisMax].
// Safety: always includes axisMin and axisMax as edges to prevent data loss at boundaries.
// All O2Physics histograms use 1 GeV uniform bins with integer boundaries, so axis
// boundaries are always valid bin edges for TH1::Rebin.
inline PtBinning BuildPtBinning(const std::vector<Double_t>& master, double axisMin, double axisMax) {
  const double eps = 0.01;  // tolerance for floating-point axis limits
  std::vector<Double_t> sel;
  for (auto e : master) {
    if (e >= axisMin - eps && e <= axisMax + eps) sel.push_back(e);
  }
  if (sel.size() < 2) {
    std::cerr << "[Binning] WARNING: fewer than 2 edges in [" << axisMin << ", " << axisMax
              << "], falling back to first/last" << std::endl;
    sel = {axisMin, axisMax};
  }
  // Safety: ensure axis boundaries are included to prevent data loss.
  // If the first selected edge is above axisMin, prepend axisMin.
  // If the last selected edge is below axisMax, append axisMax.
  if (sel.front() > axisMin + eps) {
    sel.insert(sel.begin(), axisMin);
    std::cerr << "[Binning] Safety: prepended axis min " << axisMin << std::endl;
  }
  if (sel.back() < axisMax - eps) {
    sel.push_back(axisMax);
    std::cerr << "[Binning] Safety: appended axis max " << axisMax << std::endl;
  }
  PtBinning b(sel);
  std::cerr << "[Binning] Built " << b.nBins << " bins in [" << sel.front() << ", " << sel.back() << "]" << std::endl;
  return b;
}

// Auto-detect from a 1D histogram axis
inline PtBinning DetectRecoBinning(TH1* h) {
  return BuildPtBinning(GetMasterRecoBins(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
}
inline PtBinning DetectTruthBinning(TH1* h) {
  double xmin = TMath::Max(0.0, h->GetXaxis()->GetXmin());  // particle-level pT >= 0
  return BuildPtBinning(GetMasterTruthBins(), xmin, h->GetXaxis()->GetXmax());
}

// Auto-detect from a 2D response matrix (both axes)
// Truth (particle-level) pT is always >= 0.  UE-subtracted response matrices use
// jetPtAxisRhoAreaSub [-200, 200] for both axes, but negative truth pT has no content.
// Clipping truth at 0 ensures the truth binning {0,1,...,200} matches POWHEG and MC truth
// histograms — otherwise TH1::Divide silently fails due to 26 vs 25 bin mismatch.
inline void DetectResponseBinning(TH2* h2, PtBinning& recoBin, PtBinning& truthBin) {
  recoBin  = BuildPtBinning(GetMasterRecoBins(),  h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax());
  double truthYmin = TMath::Max(0.0, h2->GetYaxis()->GetXmin());
  if (truthYmin != h2->GetYaxis()->GetXmin()) {
    std::cerr << "[Binning] DetectResponseBinning: clipping truth Y from "
              << h2->GetYaxis()->GetXmin() << " to " << truthYmin
              << " (particle-level pT >= 0)" << std::endl;
  }
  truthBin = BuildPtBinning(GetMasterTruthBins(), truthYmin, h2->GetYaxis()->GetXmax());
}

// Convenience: first/last reco bin edges (after detection)
inline Double_t GetRecoPtMin(const PtBinning& b) { return b.edges.front(); }
inline Double_t GetRecoPtMax(const PtBinning& b) { return b.edges.back(); }

// Safe rebin: rebin hSrc to target binning even if hSrc has a narrower axis range.
// - If source covers target: uses standard TH1::Rebin (fast, exact).
// - If source is narrower: creates new histogram, fills from source, zero-pads outside.
// This is essential when data (e.g. [0,200]) must match response matrix reco (e.g. [0,300]).
inline TH1* RebinToTarget(TH1* hSrc, const PtBinning& target, const char* name) {
  if (!hSrc) return nullptr;
  double srcMin = hSrc->GetXaxis()->GetXmin();
  double srcMax = hSrc->GetXaxis()->GetXmax();
  double tgtMin = target.edges.front();
  double tgtMax = target.edges.back();

  // Standard path: source covers target range
  if (srcMin <= tgtMin + 0.01 && srcMax >= tgtMax - 0.01) {
    TH1* h = hSrc->Rebin(target.nBins, name, target.bins);
    if (h) h->SetDirectory(nullptr);
    return h;
  }

  // Extended path: source is narrower, create target and transfer content
  std::cerr << "[Binning] RebinToTarget: source [" << srcMin << ", " << srcMax
            << "] narrower than target [" << tgtMin << ", " << tgtMax
            << "], zero-padding outside" << std::endl;
  TH1D* hOut = new TH1D(name, hSrc->GetTitle(), target.nBins, target.bins);
  hOut->SetDirectory(nullptr);
  hOut->Sumw2();
  for (int i = 1; i <= hSrc->GetNbinsX(); i++) {
    double pt = hSrc->GetXaxis()->GetBinCenter(i);
    int tBin = hOut->FindBin(pt);
    if (tBin >= 1 && tBin <= target.nBins) {
      hOut->SetBinContent(tBin, hOut->GetBinContent(tBin) + hSrc->GetBinContent(i));
      double e1 = hOut->GetBinError(tBin);
      double e2 = hSrc->GetBinError(i);
      hOut->SetBinError(tBin, TMath::Sqrt(e1 * e1 + e2 * e2));
    }
  }
  return hOut;
}

// Object names (mode-dependent for jet histograms)
inline const char* GetTrackPtObj() { return "h_track_pt"; }  // Track names are same for both modes
inline const char* GetTrackEtaObj() { return "h2_track_eta_track_phi"; }
inline const char* GetTrackPhiObj() { return "h2_track_eta_track_phi"; }
inline const char* GetJetPtObj() {
  return (GetCurrentMode() == kUE) ? "h3_jet_pt_jet_eta_jet_phi_rhoareasubtracted" : "h3_jet_pt_jet_eta_jet_phi";
}
inline const char* GetJetEtaObj() {
  return (GetCurrentMode() == kUE) ? "h3_jet_pt_jet_eta_jet_phi_rhoareasubtracted" : "h3_jet_pt_jet_eta_jet_phi";
}
inline const char* GetJetPhiObj() {
  return (GetCurrentMode() == kUE) ? "h3_jet_pt_jet_eta_jet_phi_rhoareasubtracted" : "h3_jet_pt_jet_eta_jet_phi";
}
inline const char* GetJetPtMCPObj() {
  return (GetCurrentMode() == kUE) ? "h_jet_pt_part_rhoareasubtracted" : "h_jet_pt_part";
}
inline const char* GetJetAreaObj() {
  return (GetCurrentMode() == kUE) ? "h2_jet_pt_jet_area_rhoareasubtracted" : "h2_jet_pt_jet_area";
}
inline const char* GetConstituentPtObj() {
  return (GetCurrentMode() == kUE) ? "h2_jet_pt_track_pt_rhoareasubtracted" : "h2_jet_pt_track_pt";
}
inline const char* GetNtracksObj() {
  return (GetCurrentMode() == kUE) ? "h2_jet_pt_jet_ntracks_rhoareasubtracted" : "h2_jet_pt_jet_ntracks";
}
inline const char* GetJetResolutionObj() {
  return (GetCurrentMode() == kUE) ? "h2_jet_pt_mcp_jet_pt_diff_matchedgeo_rhoareasubtracted" : "h2_jet_pt_mcp_jet_pt_diff_matchedgeo";
}
inline const char* GetLumiTVXObj() { return "eventselection-run3/luminosity/hLumiTVX"; }

// Get 1D jet pT from a ROOT file, applying eta cut when jetEtaCut >= 0.
// When jetEtaCut < 0 (default): returns 1D h_jet_pt directly (R-dependent eta already applied at production).
// When jetEtaCut >= 0: projects 3D h3_jet_pt_jet_eta_jet_phi within |eta| < jetEtaCut.
// Caller owns the returned histogram (SetDirectory(0) applied).
inline TH1* GetJetPt1D(TFile* file, const char* dir, const char* uniqueName) {
  const auto& set = GetCurrentConfigSet();
  if (set.jetEtaCut >= 0) {
    // Fixed eta: project 3D within |eta| < jetEtaCut
    TH3* h3 = (TH3*)file->Get(Form("%s/%s", dir, GetJetPtObj()));
    if (!h3) return nullptr;
    int binLo = h3->GetYaxis()->FindBin(-set.jetEtaCut + 1e-6);
    int binHi = h3->GetYaxis()->FindBin( set.jetEtaCut - 1e-6);
    TH1* h1 = h3->ProjectionX(uniqueName, binLo, binHi);
    h1->SetDirectory(0);
    return h1;
  } else {
    // R-dependent: use pre-computed 1D histogram
    const char* objName = (GetCurrentMode() == kUE) ? "h_jet_pt_rhoareasubtracted" : "h_jet_pt";
    TH1* h1 = (TH1*)file->Get(Form("%s/%s", dir, objName));
    if (!h1) return nullptr;
    h1 = (TH1*)h1->Clone(uniqueName);
    h1->SetDirectory(0);
    return h1;
  }
}

// Helper function to add "UE subtracted" label to jet plots (not track plots)
// Default position: left-top corner, below the frame top line
inline void AddUESubtractedLabel(TPad* pad, Double_t x = 0.18, Double_t y = 0.82, Double_t size = 0.045) {
  if (GetCurrentMode() == kUE && pad) {
    pad->cd();
    TLatex* ueText = new TLatex();
    ueText->SetNDC();
    ueText->SetTextFont(42);
    ueText->SetTextSize(size);
    ueText->DrawLatex(x, y, "UE subtracted");
  }
}

// Title strings
inline TString GetTrackPtTitleX() { return "#it{p}_{T, track}^{reco} (GeV/#it{c})"; }
inline TString GetTrackEtaTitleX() { return "#it{#eta}_{track}"; }
inline TString GetTrackPhiTitleX() { return "#it{#varphi}_{track}"; }
inline TString GetJetPtTitleX() { return "#it{p}_{T, jet}^{reco} (GeV/#it{c})"; }
inline TString GetJetPtGenTitleX() { return "#it{p}_{T, jet}^{true} (GeV/#it{c})"; }
inline TString GetJetEtaTitleX() { return "#it{#eta}_{jet}"; }
inline TString GetJetPhiTitleX() { return "#it{#varphi}_{jet}"; }
inline TString GetJRETitleX() { return "#it{p}_{T, jet}^{true} (GeV/#it{c})"; }
inline TString GetJRETitleY() { return "1/#it{N}_{evt} d#it{N}/d#it{p}_{T}"; }
inline TString GetJRPTitleX() { return "#it{p}_{T, jet}^{reco} (GeV/#it{c})"; }
inline TString GetJRPTitleY() { return "1/#it{N}_{evt} d#it{N}/d#it{p}_{T}"; }

// Load POWHEG NLO + PYTHIA8 (CT18nlo, Monash tune, dijet, 13.6 TeV) cross section
// Returns d²σ/dpT/dη [mb/(GeV/c)] for the given R and current UE mode.
// Available R: 0.2–0.6.  Returns nullptr for R outside this range.
// Histogram naming: InclusiveJetXSection_R{N} (non-UE), BkgSubtractedJetXSection_R{N} (UE-sub)
// POWHEG eta convention: |η_jet| < 0.9-R, matching our analysis.
inline TH1* LoadPowhegForR(double R) {
  int iR = (int)(R * 10 + 0.5);  // 0.2→2, ..., 0.6→6
  if (iR < 2 || iR > 6) {
    return nullptr;  // No POWHEG for R=0.1 or R=0.7
  }

  TString powhegFile = Form("%s/cernbox/workspace/O2Physics/jets/POWHEG/Pythia8bJetSpectra_dijet_CT18nlo_MonashTune.root",
                            gSystem->HomeDirectory());
  TFile* fPOWHEG = TFile::Open(powhegFile.Data(), "READ");
  if (!fPOWHEG || fPOWHEG->IsZombie()) {
    std::cerr << "[Warning] Cannot open POWHEG file: " << powhegFile.Data() << std::endl;
    if (fPOWHEG) { fPOWHEG->Close(); delete fPOWHEG; }
    return nullptr;
  }

  // Select histogram based on UE mode
  TString hName;
  if (GetCurrentMode() == kUE) {
    hName = Form("BkgSubtractedJetXSection_R%d", iR);
  } else {
    hName = Form("InclusiveJetXSection_R%d", iR);
  }

  TH1D* hRaw = (TH1D*)fPOWHEG->Get(hName.Data());
  if (!hRaw) {
    std::cerr << "[Warning] Cannot find POWHEG histogram: " << hName.Data() << std::endl;
    fPOWHEG->Close(); delete fPOWHEG;
    return nullptr;
  }

  // Get number of events for normalization
  TH1* hNEvent = (TH1*)fPOWHEG->Get("hNEvent");
  Double_t nEvents = 20000000.0;  // fallback
  if (hNEvent) {
    nEvents = hNEvent->GetEntries();
  }

  // Clone and detach from file
  TH1* hClone = (TH1*)hRaw->Clone(Form("hPOWHEG_%s_R%.1f_raw", (GetCurrentMode() == kUE) ? "UE" : "NonUE", R));
  hClone->SetDirectory(0);
  fPOWHEG->Close(); delete fPOWHEG;

  // Rebin to analysis truth binning BEFORE normalization (content is raw weighted counts).
  // POWHEG has uniform 1 GeV bins; analysis uses non-uniform bins {0,1,...,10,12,14,...,200}.
  // All analysis bin edges are at integer GeV, so they align with POWHEG's 1 GeV grid.
  const Double_t* truthEdges = GetPtbinGen();
  Int_t nTruthBins = GetNptBinsGen();  // 25
  TH1* hRebinned = hClone->Rebin(nTruthBins, Form("hPOWHEG_%s_R%.1f", (GetCurrentMode() == kUE) ? "UE" : "NonUE", R), truthEdges);
  hRebinned->SetDirectory(0);
  delete hClone;

  // Normalize: d²σ/dpT/dη = histogram / nEvents / binWidth / Δη
  // Histograms are filled with event weight (already in mb via ×1e-9 in RunPythia8.C).
  // Scale(1/nEvents, "width") gives dσ/dpT [mb/(GeV/c)] integrated over η window.
  // Then divide by Δη = 2×(0.9-R) for d²σ/dpT/dη.
  hRebinned->Scale(1.0 / nEvents, "width");
  Double_t deltaEta = 2.0 * (0.9 - R);  // POWHEG uses same 0.9-R convention
  if (deltaEta > 0) {
    hRebinned->Scale(1.0 / deltaEta);
  }

  std::cerr << "[Info] Loaded POWHEG NLO for R=" << R
            << " (" << hName.Data() << "): nEvents=" << Form("%.0f", nEvents)
            << ", deltaEta=" << Form("%.1f", deltaEta)
            << ", rebinned to " << nTruthBins << " analysis bins" << std::endl;

  return hRebinned;
}

// ============================================
// Standalone Model Loading (PYTHIA8 variants, Herwig7)
// ============================================

struct StandaloneModel {
  TString name;      // "Monash", "Rope", "Shoving", "Herwig7"
  TString fileName;  // just filename
  TString legLabel;  // legend label
  Int_t lineStyle;   // distinct per model (for AllR canvases where color=R)
  Int_t lineWidth;   // line width
  Color_t lineColor; // model-specific color (for per-R canvases where R is fixed)
};

inline std::vector<StandaloneModel> GetStandaloneModels() {
  std::vector<StandaloneModel> models;
  models.push_back({"Monash",  "merged_JJMC_PYTHIA8_Monash_pp_13600GeV.root",  "PYTHIA8 Monash",  2, 2, kBlue+1});
  models.push_back({"Rope",    "merged_JJMC_PYTHIA8_Rope_pp_13600GeV.root",    "PYTHIA8 Rope",    7, 2, kRed+1});
  models.push_back({"Shoving", "merged_JJMC_PYTHIA8_Shoving_pp_13600GeV.root", "PYTHIA8 Shoving", 9, 2, kGreen+2});
  // Herwig7 disabled: pT-hat stitching bug causes saw-tooth artifacts. Re-enable after KIAF re-production.
  // models.push_back({"Herwig7", "merged_JJMC_Herwig7_Herwig_pp_13600GeV.root",  "Herwig7",         4, 2, kMagenta+1});
  return models;
}

static const Color_t kPowhegColor = kOrange+1;

inline TString GetStandaloneModelDir() {
  return Form("%s/cernbox/workspace/O2Physics/jets/PYTHIA_standalone_KIAF", gSystem->HomeDirectory());
}

inline TH1* LoadStandaloneModelXsec(const StandaloneModel& model, double R) {
  int iR = (int)(R * 10 + 0.5);
  if (iR < 1 || iR > 7) {
    std::cerr << "[Model] Invalid R=" << R << " for " << model.name << std::endl;
    return nullptr;
  }

  TString filePath = Form("%s/%s", GetStandaloneModelDir().Data(), model.fileName.Data());
  TFile* f = TFile::Open(filePath.Data(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "[Model] Cannot open " << filePath << std::endl;
    if (f) { f->Close(); delete f; }
    return nullptr;
  }

  // Select histogram based on UE mode
  TString hName;
  if (GetCurrentMode() == kUE) {
    hName = Form("hJetXsec_UEsub_R0%d", iR);
  } else {
    hName = Form("hJetXsec_R0%d", iR);
  }

  TH1* hRaw = (TH1*)f->Get(hName.Data());
  if (!hRaw) {
    std::cerr << "[Model] Cannot find " << hName << " in " << model.fileName << std::endl;
    f->Close(); delete f;
    return nullptr;
  }

  // Get number of events
  TH1* hNev = (TH1*)f->Get("hNeventTotal");
  Double_t nEvents = 1.0;
  if (hNev) {
    nEvents = hNev->GetBinContent(1);
  } else {
    std::cerr << "[Model] WARNING: hNeventTotal not found in " << model.fileName << std::endl;
  }

  // Clone and detach
  TH1* hClone = (TH1*)hRaw->Clone(Form("hModel_%s_%s_R%.1f_raw",
    model.name.Data(), (GetCurrentMode() == kUE) ? "UE" : "NonUE", R));
  hClone->SetDirectory(0);
  f->Close(); delete f;

  // Rebin to analysis truth binning
  const Double_t* truthEdges = GetPtbinGen();
  Int_t nTruthBins = GetNptBinsGen();
  TH1* hRebinned = hClone->Rebin(nTruthBins,
    Form("hModel_%s_%s_R%.1f", model.name.Data(), (GetCurrentMode() == kUE) ? "UE" : "NonUE", R),
    truthEdges);
  hRebinned->SetDirectory(0);
  delete hClone;

  // Normalize: d²σ/dpT/dη
  // Histogram content is ALREADY dσ [mb] per pT bin (weighted by σ_gen/N_events in production).
  // Do NOT divide by nEvents — just divide by bin width and eta window.
  hRebinned->Scale(1.0, "width");
  Double_t deltaEta = 2.0 * (0.9 - R);
  if (deltaEta > 0) {
    hRebinned->Scale(1.0 / deltaEta);
  }

  return hRebinned;
}

inline TH1* LoadStandaloneModelRatio(const StandaloneModel& model, double numR, double denR) {
  int iNum = (int)(numR * 10 + 0.5);
  int iDen = (int)(denR * 10 + 0.5);

  TString filePath = Form("%s/%s", GetStandaloneModelDir().Data(), model.fileName.Data());
  TFile* f = TFile::Open(filePath.Data(), "READ");
  if (!f || f->IsZombie()) {
    if (f) { f->Close(); delete f; }
    return nullptr;
  }

  TString hName;
  if (GetCurrentMode() == kUE) {
    hName = Form("hJetXsecRatio_UEsub_R0%d_R0%d", iNum, iDen);
  } else {
    hName = Form("hJetXsecRatio_R0%d_R0%d", iNum, iDen);
  }

  TH1* hRaw = (TH1*)f->Get(hName.Data());
  if (!hRaw) {
    std::cerr << "[Model] Cannot find " << hName << " in " << model.fileName << std::endl;
    f->Close(); delete f;
    return nullptr;
  }

  TH1* hClone = (TH1*)hRaw->Clone(Form("hModelRatio_%s_R%.1f_R%.1f", model.name.Data(), numR, denR));
  hClone->SetDirectory(0);
  f->Close(); delete f;

  return hClone;
}

// Load standalone model per-event normalized yield: (1/N_evt) dN/dpT
// This matches the MC particle-level normalization in DrawJetPtPart
inline TH1* LoadStandaloneModelNormYield(const StandaloneModel& model, double R) {
  int iR = (int)(R * 10 + 0.5);
  if (iR < 1 || iR > 7) {
    std::cerr << "[Model] Invalid R=" << R << " for NormYield " << model.name << std::endl;
    return nullptr;
  }

  TString filePath = Form("%s/%s", GetStandaloneModelDir().Data(), model.fileName.Data());
  TFile* f = TFile::Open(filePath.Data(), "READ");
  if (!f || f->IsZombie()) {
    if (f) { f->Close(); delete f; }
    return nullptr;
  }

  TString hName;
  if (GetCurrentMode() == kUE) {
    hName = Form("hJetNormYield_UEsub_R0%d", iR);
  } else {
    hName = Form("hJetNormYield_R0%d", iR);
  }

  TH1* hRaw = (TH1*)f->Get(hName.Data());
  if (!hRaw) {
    std::cerr << "[Model] Cannot find " << hName << " in " << model.fileName << std::endl;
    f->Close(); delete f;
    return nullptr;
  }

  TH1* hClone = (TH1*)hRaw->Clone(Form("hModelNormYield_%s_%s_R%.1f_raw",
    model.name.Data(), (GetCurrentMode() == kUE) ? "UE" : "NonUE", R));
  hClone->SetDirectory(0);
  f->Close(); delete f;

  // Rebin to analysis truth binning
  const Double_t* truthEdges = GetPtbinGen();
  Int_t nTruthBins = GetNptBinsGen();
  TH1* hRebinned = hClone->Rebin(nTruthBins,
    Form("hModelNormYield_%s_%s_R%.1f", model.name.Data(), (GetCurrentMode() == kUE) ? "UE" : "NonUE", R),
    truthEdges);
  hRebinned->SetDirectory(0);
  delete hClone;

  // hJetNormYield is (1/N_evt) dN per pT bin per unit η
  // Scale by bin width → (1/N_evt) d²N/(dη dpT)
  // Then multiply by Δη → (1/N_evt) dN/dpT (matching O2Physics MC truth format)
  Double_t deltaEta = 2.0 * (0.9 - R);
  hRebinned->Scale(1.0, "width");
  if (deltaEta > 0) {
    hRebinned->Scale(deltaEta);  // multiply by Δη to get total per-event yield
  }

  return hRebinned;
}

// Forward declarations (defined below GetAllRConfigSets)
inline Double_t Nevents(const char *fileName, const char *eventDir, const char *eventObj, Int_t ifMCP = 0, bool isJJ = false);
inline Double_t GetRCTPassFraction(const char *MCfileName, const char *Dir, bool isJJ = false);
inline const char* GetLumiCalcDir();
inline Double_t GetSel8Efficiency(const char *fileName, bool isMC = true);

// Load MC particle-level truth as cross-section d²σ/dpT/dη [mb/GeV]
// For JJ MC: histogram is already cross-section weighted (σ_gen/N_gen per event)
//   → just rebin, divide by bin width and Δη
// For MB MC: raw counts → (N_jets / N_INEL) × σ_INEL / ΔpT / Δη
//   σ_INEL ≈ 79 mb (PYTHIA8 total inelastic at 13.6 TeV)
inline TH1* LoadMCTruthXsec(const RConfig& config) {
  const Double_t kSigmaINEL = 79.0;  // mb, total inelastic (ND+SD+DD) at 13.6 TeV

  TFile* f = TFile::Open(config.GetMcFile().Data(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "[MCTruth] Cannot open MC file: " << config.GetMcFile() << std::endl;
    if (f) { f->Close(); delete f; }
    return nullptr;
  }

  // Get particle-level jet pT (mode-dependent)
  const char* truthObj = GetJetPtMCPObj();
  TH1* hRaw = (TH1*)f->Get(Form("%s/%s", config.mcDir.Data(), truthObj));
  if (!hRaw) {
    std::cerr << "[MCTruth] Cannot find " << truthObj << " in " << config.mcDir << std::endl;
    f->Close(); delete f;
    return nullptr;
  }

  TH1* hClone = (TH1*)hRaw->Clone(Form("hMCTruth_raw_R%.1f_%s", config.R,
    (GetCurrentMode() == kUE) ? "UE" : "NonUE"));
  hClone->SetDirectory(0);

  // For MB MC: get N_INEL from h_mccollisions bin 0.5 (all MC collisions)
  Double_t nINEL = 1.0;
  if (!config.isJJ) {
    TH1* hMcColl = (TH1*)f->Get(Form("%s/h_mccollisions", config.mcDir.Data()));
    if (hMcColl) {
      nINEL = hMcColl->GetBinContent(hMcColl->FindBin(0.5));
      std::cerr << "[MCTruth] MB MC N_INEL = " << nINEL << std::endl;
    } else {
      std::cerr << "[MCTruth] WARNING: h_mccollisions not found, cannot normalize MB MC truth" << std::endl;
    }
  }
  f->Close(); delete f;

  // Rebin to analysis truth binning
  const Double_t* truthEdges = GetPtbinGen();
  Int_t nTruthBins = GetNptBinsGen();
  TH1* hRebinned = hClone->Rebin(nTruthBins,
    Form("hMCTruthXsec_R%.1f_%s", config.R, (GetCurrentMode() == kUE) ? "UE" : "NonUE"),
    truthEdges);
  hRebinned->SetDirectory(0);
  delete hClone;

  // Normalize to d²σ/dpT/dη [mb/GeV]
  Double_t deltaEta = GetDeltaEta(config.R);
  if (config.isJJ) {
    // JJ MC: histogram content = Σ(σ_gen/N_gen × n_jets) = dσ per bin [mb]
    hRebinned->Scale(1.0, "width");
    if (deltaEta > 0) hRebinned->Scale(1.0 / deltaEta);
  } else {
    // MB MC: raw counts → cross-section
    if (nINEL > 0 && deltaEta > 0) {
      hRebinned->Scale(kSigmaINEL / (nINEL * deltaEta), "width");
    }
  }

  std::cerr << "[MCTruth] Loaded R=" << config.R << " (isJJ=" << config.isJJ
            << "), integral=" << Form("%.4e", hRebinned->Integral("width"))
            << " mb" << std::endl;

  return hRebinned;
}

// Load MC particle-level truth as invariant yield (1/Nevt) d²N/(dpT dη)
// Uses nevtsMCP (weighted selected events for JJ, selected events for MB)
// Consistent with DrawInvariantYield normalization
inline TH1* LoadMCTruthNormYield(const RConfig& config) {
  TFile* f = TFile::Open(config.GetMcFile().Data(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "[MCTruthYield] Cannot open MC file: " << config.GetMcFile() << std::endl;
    if (f) { f->Close(); delete f; }
    return nullptr;
  }

  const char* truthObj = GetJetPtMCPObj();
  TH1* hRaw = (TH1*)f->Get(Form("%s/%s", config.mcDir.Data(), truthObj));
  if (!hRaw) {
    std::cerr << "[MCTruthYield] Cannot find " << truthObj << " in " << config.mcDir << std::endl;
    f->Close(); delete f;
    return nullptr;
  }

  TH1* hClone = (TH1*)hRaw->Clone(Form("hMCTruthYield_raw_R%.1f_%s", config.R,
    (GetCurrentMode() == kUE) ? "UE" : "NonUE"));
  hClone->SetDirectory(0);

  // Close file BEFORE calling Nevents() to avoid TFile::Open cache collision
  f->Close(); delete f;

  // Get number of selected events (bin 5.5 = all cuts)
  // For JJ MC: returns weighted count from h_mccollisions_weighted
  Double_t nevtsMCP = Nevents(config.GetMcFile().Data(), config.mcDir.Data(),
                               GetEventObj(), 1, config.isJJ);

  // Rebin to analysis truth binning
  const Double_t* truthEdges = GetPtbinGen();
  Int_t nTruthBins = GetNptBinsGen();
  TH1* hRebinned = hClone->Rebin(nTruthBins,
    Form("hMCTruthYield_R%.1f_%s", config.R, (GetCurrentMode() == kUE) ? "UE" : "NonUE"),
    truthEdges);
  hRebinned->SetDirectory(0);
  delete hClone;

  // Normalize to (1/Nevt) d²N/(dpT dη)
  Double_t deltaEta = GetDeltaEta(config.R);
  if (nevtsMCP > 0 && deltaEta > 0) {
    hRebinned->Scale(1.0 / (nevtsMCP * deltaEta), "width");
  }

  std::cerr << "[MCTruthYield] Loaded R=" << config.R << " (isJJ=" << config.isJJ
            << ", nevtsMCP=" << Form("%.0f", nevtsMCP)
            << "), integral=" << Form("%.4e", hRebinned->Integral("width")) << std::endl;

  return hRebinned;
}

// Configuration functions
inline std::vector<Color_t> GetRColors() {
  std::vector<Color_t> colors;
  colors.push_back(kBlack);
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kGreen+2);
  colors.push_back(kMagenta+2);
  colors.push_back(kCyan+2);
  colors.push_back(kOrange+7);
  return colors;
}

// Map R value to a fixed color: R=0.1→black, 0.2→red, ..., 0.7→orange
// Use this on combined canvases so colors are stable regardless of which R values are present.
inline Color_t GetColorForR(double R) {
  int idx = (int)(R * 10 + 0.5) - 1;  // 0.1→0, 0.2→1, ..., 0.7→6
  static const Color_t map[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta+2, kCyan+2, kOrange+7};
  if (idx >= 0 && idx < 7) return map[idx];
  return kBlack;
}

// Model-specific colors (distinct from per-R data colors)
inline Color_t GetPYTHIAColor()  { return kGray+2; }    // dark gray
inline Color_t GetPOWHEGColor()  { return kViolet+1; }  // purple
inline Color_t GetMCTruthColor() { return kGray+2; }    // dark gray (same role as PYTHIA in UE mode)

// Build color vector from R configs (each R gets its fixed color)
inline std::vector<Color_t> GetRColorsForConfigs(const std::vector<RConfig>& configs) {
  std::vector<Color_t> colors;
  for (const auto& cfg : configs) {
    colors.push_back(GetColorForR(cfg.R));
  }
  return colors;
}

// Returns the R configs for the currently active config set
inline std::vector<RConfig> GetRConfigs() {
  return GetCurrentConfigSet().configs;
}

// All config sets: uncomment/comment sets to control which run
inline std::vector<RConfigSet> GetAllRConfigSets() {
  std::vector<RConfigSet> sets;

  // ---- NonUE sets ----

  // // 2022 small data + MB MC (LHC24f3c) tunerA
  // { RConfigSet s;
  //   s.name = "LHC22o-pass7-small_LHC24f3c";
  //   s.outputDir = "../plots/RDependentComparison/LHC22o-pass7-small_LHC24f3c";
  //   s.mode = kNonUE;
  //   s.configs.push_back(RConfig(0.1, "562004", "561898", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.1"));
  //   s.configs.push_back(RConfig(0.2, "562003", "561899", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.2"));
  //   s.configs.push_back(RConfig(0.3, "562005", "561900", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.3"));
  //   s.configs.push_back(RConfig(0.4, "562006", "561901", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.4"));
  //   s.configs.push_back(RConfig(0.5, "562007", "561902", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.5"));
  //   s.configs.push_back(RConfig(0.6, "562008", "561903", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.6"));
  //   s.configs.push_back(RConfig(0.7, "562009", "561904", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.7"));
  //   sets.push_back(s); }

  // // 2022 small data + JJ MC (LHC25a2b) tunerA
  // { RConfigSet s;
  //   s.name = "LHC22o-pass7-small_LHC25a2b_570220";
  //   s.outputDir = "../plots/RDependentComparison/LHC22o-pass7-small_LHC25a2b_570220";
  //   s.mode = kNonUE;
  //   s.configs.push_back(RConfig(0.1, "562004", "570218", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.1", true));
  //   s.configs.push_back(RConfig(0.2, "562003", "570219", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.2", true));
  //   s.configs.push_back(RConfig(0.3, "562005", "574218", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.3", true));
  //   s.configs.push_back(RConfig(0.4, "562006", "570220", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.4", true));
  //   s.configs.push_back(RConfig(0.5, "562007", "570221", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.5", true));
  //   s.configs.push_back(RConfig(0.6, "562008", "570222", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.6", true));
  //   s.configs.push_back(RConfig(0.7, "562009", "570223", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.7", true));
  //   sets.push_back(s); }

  // // // ---- UE sets ----

  // // 2022 small data + MB MC (LHC24f3c) pTsmearing1p5
  // { RConfigSet s;
  //   s.name = "LHC22o-pass7-small_LHC24f3c_pTsmearing1p5";
  //   s.outputDir = "../plots/RDependentComparison/LHC22o-pass7-small_LHC24f3c_pTsmearing1p5";
  //   s.runBothModes = true;
  // //   s.configs.push_back(RConfig(0.1, "592451", "598076", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.1"));
  //   s.configs.push_back(RConfig(0.2, "592377", "594013", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.2"));
  // //   s.configs.push_back(RConfig(0.3, "592379", "598075", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.3"));
  //   s.configs.push_back(RConfig(0.4, "592378", "594014", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.4"));
  // //   s.configs.push_back(RConfig(0.5, "592380", "598074", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.5"));
  //   s.configs.push_back(RConfig(0.6, "592381", "598073", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.6"));
  // //   s.configs.push_back(RConfig(0.7, "592382", "594015", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.7"));
  //   sets.push_back(s); }


  // // 2022 full data + JJ MC (LHC25a2b) Tuner A
  // { RConfigSet s;
  //   s.name = "LHC22o-pass7-small_LHC25a2b_tunerA";
  //   s.outputDir = "../plots/RDependentComparison/LHC22o-pass7-small_LHC25a2b_tunerA";
  //   s.runBothModes = true;
  //   s.configs.push_back(RConfig(0.1, "592451", "615372"/*LHC25a2 w/ myTune*/, "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.1", true));
  //   s.configs.push_back(RConfig(0.2, "592377", "599366", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.2", true));
  //   s.configs.push_back(RConfig(0.3, "592379", "599367", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.3", true));
  //   s.configs.push_back(RConfig(0.4, "592378", "605663", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.4", true));
  //   s.configs.push_back(RConfig(0.5, "592380", "599368", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.5", true));
  //   s.configs.push_back(RConfig(0.6, "592381", "599369", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.6", true));
  //   s.configs.push_back(RConfig(0.7, "592382", "599370", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.7", true));
  //   sets.push_back(s); 
  //  }

  // // 2022 small data + JJ MC (LHC25a2b) pTsmearing1p5
  // { RConfigSet s;
  //   s.name = "LHC22o-pass7-small_LHC25a2b_pTsmearing1p5_woRdepMatching";
  //   s.outputDir = "../plots/RDependentComparison/LHC22o-pass7-small_LHC25a2b_pTsmearing1p5_woRdepMatching";
  //   s.runBothModes = true;
  // //   s.configs.push_back(RConfig(0.1, "592451", "594021", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.1", true));
  //   s.configs.push_back(RConfig(0.2, "592377", "594022"/*LHC25a2*/, "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.2", true));
  // //   s.configs.push_back(RConfig(0.3, "592379", "594023"/*LHC25a2*/, "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.3", true));
  //   s.configs.push_back(RConfig(0.4, "592378", "593757", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.4", true));
  //   // s.configs.push_back(RConfig(0.5, "592380", "593758", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.5", true));
  //   s.configs.push_back(RConfig(0.6, "592381", "593759", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.6", true));
  // //   s.configs.push_back(RConfig(0.7, "592382", "593756", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.7", true));
  //   sets.push_back(s);
  // }

  // // 2022 full data + JJ MC (LHC25a2) MyTuner
  // { RConfigSet s;
  //   s.name = "LHC22o-pass7-small_LHC25a2_MyTuner";
  //   s.outputDir = "../plots/RDependentComparison/LHC22o-pass7-small_LHC25a2_MyTuner";
  //   s.runBothModes = true;
  // //   s.configs.push_back(RConfig(0.1, "592451", /*"616540",*/ "617213",/*LHC252b*/ "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.1", true));
  //   s.configs.push_back(RConfig(0.2, "592377", /*"616541",*/ "617214",/*LHC252b*/ "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.2", true));
  // //   s.configs.push_back(RConfig(0.3, "592379", "616542", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.3", true));
  //   s.configs.push_back(RConfig(0.4, "592378", "616543", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.4", true));
  // //   s.configs.push_back(RConfig(0.5, "592380", "616544", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.5", true));
  //   s.configs.push_back(RConfig(0.6, "592381", "619926", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.6", true));
  // //   s.configs.push_back(RConfig(0.7, "592382", "619927", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.7", true));
  //   sets.push_back(s);
  // }

  // // 2022 full data + JJ MC (LHC25a2) MyTuner, fixed |eta|<0.2
  // { RConfigSet s;
  //   s.name = "LHC22o-pass7-small_LHC25a2_MyTuner_deta0p2";
  //   s.outputDir = "../plots/RDependentComparison/LHC22o-pass7-small_LHC25a2_MyTuner_deta0p2";
  //   s.runBothModes = true;
  //   s.jetEtaCut = 0.2;  // fixed |eta_jet| < 0.2 for all R
  // //   s.configs.push_back(RConfig(0.1, "592451", "618677", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.1", true));
  //   s.configs.push_back(RConfig(0.2, "592377", "618678", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.2", true));
  // //   s.configs.push_back(RConfig(0.3, "592379", "618680", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.3", true));
  //   s.configs.push_back(RConfig(0.4, "592378", "618681", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.4", true));
  // //   s.configs.push_back(RConfig(0.5, "592380", "618682", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.5", true));
  //   s.configs.push_back(RConfig(0.6, "592381", "618683", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.6", true));
  // //   s.configs.push_back(RConfig(0.7, "592382", "618679", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.7", true));
  //   sets.push_back(s);
  // }

  // // 2023 small data + MB MC (LHC23k4h) MyTuner
  // { RConfigSet s;
  //   s.name = "LHC23-pass4-Thin_small_LHC23k4h_MyTuner";
  //   s.outputDir = "../plots/RDependentComparison/LHC23-pass4-Thin_small_LHC23k4h_MyTuner";
  //   // s.mode = kUE;
  //   s.runBothModes = true;
  //   s.configs.push_back(RConfig(0.1, "608780", "606091", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.1"));
  //   s.configs.push_back(RConfig(0.2, "R02oneTrain", "606086", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.2"));
  //   s.configs.push_back(RConfig(0.3, "608779", "606087", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.3"));
  //   s.configs.push_back(RConfig(0.4, "608785", "606088", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.4"));
  //   s.configs.push_back(RConfig(0.5, "608782", "606089", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.5"));
  //   s.configs.push_back(RConfig(0.6, "608783", "606090", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.6"));
  //   s.configs.push_back(RConfig(0.7, "608784", "606092", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.7"));
  //   sets.push_back(s);
  // }

  // // 2023 small data + MB MC (LHC23k4h) pTsmearing1p5
  // { RConfigSet s;
  //   s.name = "LHC23-pass4-Thin_small_LHC23k4h_pTsmearing1p5_woRdepMatching";
  //   s.outputDir = "../plots/RDependentComparison/LHC23-pass4-Thin_small_LHC23k4h_pTsmearing1p5_woRdepMatching";
  //   s.runBothModes = true;
  //   s.configs.push_back(RConfig(0.1, "608780", "594176", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.1"));
  //   s.configs.push_back(RConfig(0.2, "608781"/*"594032"*/, "594177", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.2"));
  //   s.configs.push_back(RConfig(0.3, "608779", "594178", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.3"));
  //   s.configs.push_back(RConfig(0.4, "608785", "594179", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.4"));
  //   s.configs.push_back(RConfig(0.5, "608782", "594180", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.5"));
  //   s.configs.push_back(RConfig(0.6, "608783", "594181", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.6"));
  //   s.configs.push_back(RConfig(0.7, "608784", "594182", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.7"));
  //   sets.push_back(s);
  // }

  // 2023 small data + MB MC (LHC23k4i) custom tuner
  { RConfigSet s;
    s.name = "LHC23-pass4-Thin_small_LHC23k4i_2022customTuner";
    s.outputDir = "../plots/RDependentComparison/LHC23-pass4-Thin_small_LHC23k4i_2022customTuner";
    s.runBothModes = true;
    s.configs.push_back(RConfig(0.1, "608780", "626305", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.1"));
    s.configs.push_back(RConfig(0.2, "608781"/*"626306"*/, "594177", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.2"));
    s.configs.push_back(RConfig(0.3, "608779", "626307", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.3"));
    s.configs.push_back(RConfig(0.4, "608785", "626304", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.4"));
    s.configs.push_back(RConfig(0.5, "608782", "626308", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.5"));
    s.configs.push_back(RConfig(0.6, "608783", "626309", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.6"));
    s.configs.push_back(RConfig(0.7, "608784", "626310", "jet-spectra-charged", "jet-spectra-charged", "#it{R} = 0.7"));
    sets.push_back(s);
  }

  // Expand runBothModes: duplicate each such set with the other mode
  std::vector<RConfigSet> expanded;
  for (auto& s : sets) {
    if (s.runBothModes) {
      // Ensure primary is kNonUE (no suffix), secondary is kUE (with _UEsub)
      RConfigSet nonUE = s;
      nonUE.runBothModes = false;
      nonUE.mode = kNonUE;
      // Strip _UEsub suffix if present in the original name/outputDir
      nonUE.name.ReplaceAll("_UEsub", "");
      nonUE.outputDir.ReplaceAll("_UEsub", "");

      RConfigSet ue = s;
      ue.runBothModes = false;
      ue.mode = kUE;
      // Ensure _UEsub suffix is present
      if (!ue.name.Contains("_UEsub")) ue.name += "_UEsub";
      if (!ue.outputDir.Contains("_UEsub")) ue.outputDir += "_UEsub";

      expanded.push_back(nonUE);
      expanded.push_back(ue);
    } else {
      expanded.push_back(s);
    }
  }

  return expanded;
}

// Legacy compatibility: GetRConfigsUE() removed — all sets are in GetAllRConfigSets()

inline bool GetDrawTrack() { return true; }
inline bool GetDrawJet() { return true; }
inline bool GetDrawNtracks() { return true; }
inline bool GetDrawConstituentPt() { return true; }
inline bool GetDrawJetPtPart() { return true; }
inline bool GetDrawJetArea() { return true; }
inline bool GetDrawResponseMatrix() { return true; }
inline bool GetDrawPurityEfficiency() { return true; }
inline bool GetDrawKinematicEfficiency() { return true; }
inline bool GetDrawInvariantYield() { return true; }
inline bool GetDrawCrossSection() { return true; }
inline bool GetDrawJetResolution() { return true; }
inline bool GetDrawCustomRatios() { return true; }  // Set to true to enable custom ratio plots

// Custom ratio pairs: {numerator R, denominator R}
// Example: {{0.2, 0.3}, {0.2, 0.4}} means R=0.2/R=0.3 and R=0.2/R=0.4
// Each pair will generate 2 plots: Data ratio (solid) and MC ratio (dashed)
// Usage: Uncomment and modify the pairs below, then set GetDrawCustomRatios() to true
inline std::vector<std::pair<double, double>> GetCustomRatioPairs() {
  std::vector<std::pair<double, double>> pairs;
  // Add your custom ratio pairs here
  // Example (uncomment to use):
  pairs.push_back(std::make_pair(0.2, 0.3));
  pairs.push_back(std::make_pair(0.2, 0.4));
  pairs.push_back(std::make_pair(0.2, 0.5));
  pairs.push_back(std::make_pair(0.2, 0.6));
  pairs.push_back(std::make_pair(0.2, 0.7));
  // Add more pairs as needed:
  // pairs.push_back(std::make_pair(0.3, 0.4));
  // pairs.push_back(std::make_pair(0.1, 0.2));
  return pairs;
}

// Returns the output directory for the currently active config set
inline TString GetOutputDir() {
  return GetCurrentConfigSet().outputDir;
}

// ============================================
// Helper functions
// ============================================

inline Double_t Nevents(const char *fileName, const char *eventDir, const char *eventObj, Int_t ifMCP = 0, bool isJJ = false) {
  auto file = TFile::Open(fileName, "open");
  if (!file || file->IsZombie()) {
    std::cerr << "[Error] Cannot open file: " << fileName << std::endl;
    return 1.0;
  }

  TH1D *hNevents = nullptr;
  auto hNeventsTemp1 = (TH1D *)file->Get(Form("%s/%s", eventDir, GetEventWObj()));
  auto hNeventsTemp2 = (TH1D *)file->Get(Form("%s/%s", eventDir, eventObj));

  TH1D *hMcPNevts = nullptr;
  auto hNeventsTemp3 = (TH1D *)file->Get(Form("%s/%s", eventDir, "h_mcColl_counts_weight"));
  auto hNeventsTemp4 = (TH1D *)file->Get(Form("%s/%s", eventDir, "h_mcColl_counts"));
  auto hNeventsTemp5 = (TH1D *)file->Get(Form("%s/%s", eventDir, "h_mccollisions"));
  auto hNeventsTemp6 = (TH1D *)file->Get(Form("%s/%s", eventDir, "h_mccollisions_weighted"));

  // Data and MB MC: always use unweighted h_collisions
  // JJ MC (isJJ=true): prefer h_collisions_weighted, fall back to unweighted
  bool usedWeightedData = false;
  if (isJJ && hNeventsTemp1) {
    hNevents = hNeventsTemp1;
    usedWeightedData = true;
  } else if (hNeventsTemp2) {
    hNevents = hNeventsTemp2;
  }

  // MCP event count: JJ MC needs weighted denominator, MB/data use unweighted
  // For h_mccollisions* (new format, 10 bins [0,10]): use FindBin(5.5) = all cuts passed
  // For h_mcColl_counts* (old format): use bin 4
  Int_t mcBinIndex = 0;
  bool usedWeightedMC = false;
  if (isJJ) {
    // JJ MC: prefer weighted histograms
    if (hNeventsTemp3) {
      hMcPNevts = hNeventsTemp3;
      mcBinIndex = 4;
      usedWeightedMC = true;
    } else if (hNeventsTemp6) {
      hMcPNevts = hNeventsTemp6;
      mcBinIndex = hNeventsTemp6->FindBin(5.5);
      usedWeightedMC = true;
    } else if (hNeventsTemp4) {
      hMcPNevts = hNeventsTemp4;
      mcBinIndex = 4;
    } else if (hNeventsTemp5) {
      hMcPNevts = hNeventsTemp5;
      mcBinIndex = hNeventsTemp5->FindBin(5.5);
    }
  } else {
    // MB MC or data: always use unweighted
    if (hNeventsTemp4) {
      hMcPNevts = hNeventsTemp4;
      mcBinIndex = 4;
    } else if (hNeventsTemp5) {
      hMcPNevts = hNeventsTemp5;
      mcBinIndex = hNeventsTemp5->FindBin(5.5);
    }
  }

  // JJ MC validation: warn if weighted histograms are missing
  if (isJJ) {
    if (ifMCP == 0 && !usedWeightedData) {
      std::cerr << "[Warning] JJ MC file missing h_collisions_weighted: " << fileName << std::endl;
      std::cerr << "  -> Enable processCollisionsWeighted=true in O2Physics config" << std::endl;
    }
    if (ifMCP == 1 && !usedWeightedMC) {
      std::cerr << "[Warning] JJ MC file missing weighted MC collision histograms: " << fileName << std::endl;
      std::cerr << "  -> Enable processMCCollisionsWeighted=true and processSpectraMCPWeighted=true" << std::endl;
      std::cerr << "  -> Also enable processJetsMatchedWeighted=true for correct response matrix" << std::endl;
    }
  }

  Double_t nevents = 0;
  if (ifMCP == 0) {
    if (hNevents) {
      // bin 3.5 = "occupancycut" = events passing ALL cuts (sel8 + centrality + occupancy)
      // Must match the selection level at which jets are actually filled
      Int_t binIndex = hNevents->FindBin(3.5);
      nevents = hNevents->GetBinContent(binIndex);
    } else {
      return 1.0;
    }
  } else if (ifMCP == 1) {
    if (hMcPNevts) {
      nevents = hMcPNevts->GetBinContent(mcBinIndex);
    } else {
      return 1.0;
    }
  }

  file->Close();
  return nevents;
}

// Get RCT+sel8 pass fraction from h_mccollisions histogram
// h_mccollisions bins (from jetSpectraCharged.cxx applyMCCollisionCuts):
//   bin 1 (0.5): all MC collisions
//   bin 2 (1.5): has reco collision
//   bin 3 (2.5): after split collision cut
//   bin 4 (3.5): after sel8+RCT
//   bin 5 (4.5): after centrality
//   bin 6 (5.5): after occupancy
// Returns bin4/bin3 = fraction of collisions passing sel8+RCT after prior cuts.
// This is the COMBINED sel8+RCT efficiency — cannot be separated from this histogram.
inline Double_t GetRCTPassFraction(const char *MCfileName, const char *Dir, bool isJJ = false) {
  auto file = TFile::Open(MCfileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "[RCT] Cannot open MC file: " << MCfileName << std::endl;
    return 1.0;
  }

  // JJ MC: try weighted first; MB MC: unweighted only
  TH1 *hMcColl = nullptr;
  const char* usedName = nullptr;
  if (isJJ) {
    const char* histNames[] = {"h_mccollisions_weighted", "h_mccollisions",
                                "h_mcColl_counts_weight", "h_mcColl_counts"};
    for (int i = 0; i < 4; ++i) {
      hMcColl = (TH1 *)file->Get(Form("%s/%s", Dir, histNames[i]));
      if (hMcColl) { usedName = histNames[i]; break; }
    }
  } else {
    const char* histNames[] = {"h_mccollisions", "h_mcColl_counts"};
    for (int i = 0; i < 2; ++i) {
      hMcColl = (TH1 *)file->Get(Form("%s/%s", Dir, histNames[i]));
      if (hMcColl) { usedName = histNames[i]; break; }
    }
  }

  if (!hMcColl) {
    std::cerr << "[RCT] No h_mccollisions histogram found in " << Dir
              << " — cannot compute RCT fraction. Returning 1.0" << std::endl;
    file->Close();
    return 1.0;
  }

  // bin 3 = after split collision cut (at 2.5), bin 4 = after sel8+RCT (at 3.5)
  Int_t binBeforeRCT = hMcColl->FindBin(2.5);
  Int_t binAfterRCT  = hMcColl->FindBin(3.5);
  Double_t nBefore = hMcColl->GetBinContent(binBeforeRCT);
  Double_t nAfter  = hMcColl->GetBinContent(binAfterRCT);

  Double_t fraction = (nBefore > 0) ? (nAfter / nBefore) : 1.0;

  std::cerr << "[RCT] " << Dir << " (" << usedName << "): "
            << "before sel8+RCT=" << nBefore
            << ", after=" << nAfter
            << ", pass fraction=" << Form("%.4f", fraction) << std::endl;

  file->Close();
  return fraction;
}

// Get sel8 (event selection only, no RCT) efficiency from jet-luminosity-calculator
// The luminosity calculator "counter" histogram bins (from luminosityCalculator.cxx):
//   bin  1: BC
//   bin  2: BC+TVX
//   bin  3: BC+TVX+NoTFB
//   bin  4: BC+TVX+NoTFB+NoITSROFB
//   bin  5: Coll
//   bin  6: Coll+TVX
//   bin  7: Coll+TVX+VtxZ+Sel8              (data)
//   bin  8: Coll+TVX+VtxZ+Sel8Full          (data, full)
//   bin  9: Coll+TVX+VtxZ+Sel8FullPbPb      (data, PbPb)
//   bin 10: Coll+TVX+VtxZ+SelMC             (MC)
//   bin 11: Coll+TVX+VtxZ+SelMCFull         (MC, full)
//   bin 12: Coll+TVX+VtxZ+SelMCFullPbPb     (MC, PbPb)
//   bin 13: Coll+TVX+VtxZ+SelUnanchoredMC
//   bin 14: Coll+TVX+VtxZ+SelTVX
//   bin 15: Coll+TVX+VtxZ+Sel7
//   bin 16: Coll+TVX+VtxZ+Sel7KINT7
//   bin 17: custom
//
// For MC: sel8 efficiency = bin(10)/bin(5) = (Coll+TVX+VtxZ+SelMC) / Coll
// For data: sel8 efficiency = bin(7)/bin(5) = (Coll+TVX+VtxZ+Sel8) / Coll
// This is event-level: same for all jet R values in the same dataset.
//
// Directory in ROOT file: "jet-luminosity-calculator"
// Set kLumiCalcDir below if your task has a different directory name.

inline const char* GetLumiCalcDir() { return "jet-luminosity-calculator"; }

inline Double_t GetSel8Efficiency(const char *fileName, bool isMC = true) {
  auto file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "[Sel8] Cannot open file: " << fileName << std::endl;
    return -1.0;  // negative = not available
  }

  TH1 *hCounter = (TH1 *)file->Get(Form("%s/counter", GetLumiCalcDir()));
  if (!hCounter) {
    // Silently return — luminosityProducer is optional (cross-check only)
    file->Close();
    return -1.0;
  }

  // bin 5 = Coll (all collisions)
  // bin 7 = Coll+TVX+VtxZ+Sel8  (data)
  // bin 10 = Coll+TVX+VtxZ+SelMC (MC)
  Double_t nColl = hCounter->GetBinContent(5);
  Int_t sel8Bin = isMC ? 10 : 7;
  Double_t nSel8 = hCounter->GetBinContent(sel8Bin);
  const char* sel8Label = isMC ? "SelMC" : "Sel8";

  Double_t eff = (nColl > 0) ? (nSel8 / nColl) : -1.0;

  std::cerr << "[Sel8] " << fileName << ": "
            << "Coll=" << nColl
            << ", Coll+TVX+VtxZ+" << sel8Label << "=" << nSel8
            << ", sel8 eff=" << Form("%.4f", eff) << std::endl;

  file->Close();
  return eff;
}

// ============================================
// Decomposed normalization: per-run TVX counting
// ============================================
// Uses event-selection luminosity histograms to compute:
//   N_TVX_eff = Σ_{good RCT runs} hCounterTVXafterBCcuts[run]
//
// This is the most accurate approach:
//   - hCounterTVXafterBCcuts = TVX after sel8 BC cuts (TFBorder + ITSROFBorder)
//   - hLumiTVXafterBCcutsRCT(run, rctLabel) > 0 → run passes RCT
//   - Sum only over good runs → exact effective TVX count, properly weighted
//
// Also computes diagnostic decomposition: ε_sel8, ε_RCT, ε_sel8×RCT
//
// rctYBin: Y-axis bin in hLumiTVXafterBCcutsRCT for the RCT label
//          (CBT_hadronPID = 5, check labels in your data)

struct NormDecomposition {
  Double_t nTVX;             // hCounterTVX total (all runs)
  Double_t nTVXafterBC;      // hCounterTVXafterBCcuts total (all runs, sel8 applied)
  Double_t nTVXafterBCRCT;   // hCounterTVXafterBCcuts summed for good RCT runs only
  Double_t epsSel8;          // = nTVXafterBC / nTVX (sel8 efficiency)
  Double_t epsRCT;           // = nTVXafterBCRCT / nTVXafterBC (RCT pass fraction)
  Double_t epsSel8RCT;       // = nTVXafterBCRCT / nTVX (combined)
  // Luminosity from framework (uses σ_vis from CCDB, per-run)  [μb⁻¹]
  Double_t lumiTVX;          // hLumiTVX total (all runs, all TVX BCs)
  Double_t lumiTVXafterBC;   // hLumiTVXafterBCcuts total (all runs, after sel8 BC cuts)
  Double_t lumiTVXafterBCRCT; // hLumiTVXafterBCcutsRCT summed for good RCT runs (CBT_hadronPID)
  Double_t lumiEpsSel8;      // = lumiTVXafterBC / lumiTVX (luminosity-based sel8 eff)
  Double_t lumiEpsRCT;       // = lumiTVXafterBCRCT / lumiTVXafterBC (luminosity-based RCT eff)
  int nRunsTotal;
  int nRunsGood;
  bool valid;
};

inline NormDecomposition GetNormDecomposition(const char *dataFileName, int rctYBin = 5) {
  NormDecomposition result = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false};

  auto file = TFile::Open(dataFileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "[NormDecomp] Cannot open: " << dataFileName << std::endl;
    return result;
  }

  const char* lumiDir = "eventselection-run3/luminosity";
  TH1D *hTVX = (TH1D *)file->Get(Form("%s/hCounterTVX", lumiDir));
  TH1D *hBC  = (TH1D *)file->Get(Form("%s/hCounterTVXafterBCcuts", lumiDir));
  TH2D *hRCT = (TH2D *)file->Get(Form("%s/hLumiTVXafterBCcutsRCT", lumiDir));

  // Luminosity histograms (framework-computed, using σ_vis from CCDB) [μb⁻¹]
  TH1D *hLumiTVX   = (TH1D *)file->Get(Form("%s/hLumiTVX", lumiDir));
  TH1D *hLumiBC    = (TH1D *)file->Get(Form("%s/hLumiTVXafterBCcuts", lumiDir));

  if (!hTVX || !hBC) {
    std::cerr << "[NormDecomp] Missing hCounterTVX or hCounterTVXafterBCcuts" << std::endl;
    file->Close();
    return result;
  }

  // Print RCT label for verification
  if (hRCT) {
    TString rctLabel = hRCT->GetYaxis()->GetBinLabel(rctYBin);
    std::cerr << "[NormDecomp] RCT label for ybin " << rctYBin << ": \"" << rctLabel << "\"" << std::endl;
  }

  // Accumulate total luminosities (all runs)
  if (hLumiTVX) result.lumiTVX = hLumiTVX->Integral(1, hLumiTVX->GetNbinsX());
  if (hLumiBC)  result.lumiTVXafterBC = hLumiBC->Integral(1, hLumiBC->GetNbinsX());

  // Per-run analysis
  for (int ix = 1; ix <= hTVX->GetNbinsX(); ix++) {
    Double_t cTVX = hTVX->GetBinContent(ix);
    if (cTVX < 1) continue;

    Double_t cBC = hBC->GetBinContent(ix);
    result.nTVX += cTVX;
    result.nTVXafterBC += cBC;
    result.nRunsTotal++;

    // Check if this run passes RCT
    bool passesRCT = false;
    if (hRCT) {
      passesRCT = (hRCT->GetBinContent(ix, rctYBin) > 0);
    } else {
      passesRCT = true;  // no RCT histogram → assume all pass
    }

    if (passesRCT) {
      result.nTVXafterBCRCT += cBC;
      result.nRunsGood++;
      // Accumulate luminosity for good-RCT runs after BC cuts
      if (hRCT) {
        result.lumiTVXafterBCRCT += hRCT->GetBinContent(ix, rctYBin);
      }
    }
  }

  // Compute efficiencies (counter-based)
  if (result.nTVX > 0 && result.nTVXafterBC > 0) {
    result.epsSel8 = result.nTVXafterBC / result.nTVX;
    result.epsRCT = result.nTVXafterBCRCT / result.nTVXafterBC;
    result.epsSel8RCT = result.nTVXafterBCRCT / result.nTVX;
    result.valid = true;
  }

  // Compute luminosity-based efficiencies (for kNormSoft mode)
  if (result.lumiTVX > 0) {
    result.lumiEpsSel8 = result.lumiTVXafterBC / result.lumiTVX;
    if (result.lumiTVXafterBC > 0) {
      result.lumiEpsRCT = result.lumiTVXafterBCRCT / result.lumiTVXafterBC;
    }
  }

  std::cerr << "============================================" << std::endl;
  std::cerr << "[NormDecomp] Per-run TVX decomposition:" << std::endl;
  std::cerr << "  N_TVX (all)         = " << Form("%.0f", result.nTVX) << std::endl;
  std::cerr << "  N_TVX_sel8          = " << Form("%.0f", result.nTVXafterBC) << std::endl;
  std::cerr << "  N_TVX_sel8_RCT      = " << Form("%.0f", result.nTVXafterBCRCT) << std::endl;
  std::cerr << "  eps_sel8            = " << Form("%.5f", result.epsSel8) << std::endl;
  std::cerr << "  eps_RCT             = " << Form("%.5f", result.epsRCT) << std::endl;
  std::cerr << "  eps_sel8*RCT        = " << Form("%.5f", result.epsSel8RCT) << std::endl;
  std::cerr << "  Runs: " << result.nRunsGood << " good / " << result.nRunsTotal << " total" << std::endl;
  std::cerr << "  --- Luminosity [ub^-1] ---" << std::endl;
  std::cerr << "  L_TVX (all)         = " << Form("%.2f", result.lumiTVX) << std::endl;
  std::cerr << "  L_TVX_sel8          = " << Form("%.2f", result.lumiTVXafterBC) << std::endl;
  std::cerr << "  L_TVX_sel8_RCT      = " << Form("%.2f", result.lumiTVXafterBCRCT)
            << (result.lumiTVXafterBCRCT > 0 ? "" : "  ** MISSING **") << std::endl;
  std::cerr << "  L_eps_sel8          = " << Form("%.5f", result.lumiEpsSel8) << std::endl;
  std::cerr << "  L_eps_RCT           = " << Form("%.5f", result.lumiEpsRCT) << std::endl;
  std::cerr << "============================================" << std::endl;

  file->Close();
  return result;
}

// Get coll→sel8+RCT efficiency from DATA h_collisions histogram (for cross-check).
// h_collisions bins (from jetSpectraCharged.cxx applyCollisionCuts):
//   bin 1 (0.5): all collisions
//   bin 2 (1.5): after selectCollision (sel8+RCT+mbGap)
//   bin 3 (2.5): after centrality
//   bin 4 (3.5): after occupancy
// Returns bin(2)/bin(1) = fraction of reconstructed collisions passing sel8+RCT.
inline Double_t GetCollToSel8RCTEff(const char *dataFileName, const char *Dir) {
  auto file = TFile::Open(dataFileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "[Sel8RCT] Cannot open data file: " << dataFileName << std::endl;
    return -1.0;
  }

  // Data always uses unweighted h_collisions
  TH1 *hColl = (TH1 *)file->Get(Form("%s/h_collisions", Dir));

  if (!hColl) {
    std::cerr << "[Sel8RCT] No h_collisions histogram found in " << Dir << std::endl;
    file->Close();
    return -1.0;
  }

  Double_t nAll    = hColl->GetBinContent(hColl->FindBin(0.5));  // all collisions
  Double_t nSel8   = hColl->GetBinContent(hColl->FindBin(1.5));  // after sel8+RCT

  Double_t eff = (nAll > 0) ? (nSel8 / nAll) : -1.0;

  std::cerr << "[Sel8RCT] " << Dir << ": allColl=" << nAll
            << ", sel8+RCT=" << nSel8
            << ", eff_coll->sel8RCT=" << Form("%.5f", eff) << std::endl;

  file->Close();
  return eff;
}

// Get z-vertex ±10cm efficiency from h_collisions_Zvertex using Gaussian extrapolation.
// Fits a Gaussian in [-10,10] cm (with mean fixed to histogram mean),
// then returns Int(-10,10)/Int(-inf,+inf) of the fitted function.
// Same method as DrawZVertexExtrapolation.C.
inline Double_t GetZvtx10Efficiency(const char *dataFileName, const char *Dir,
                                    const char *plotDir = nullptr, const char *plotLabel = nullptr) {
  auto file = TFile::Open(dataFileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "[Zvtx] Cannot open data file: " << dataFileName << std::endl;
    return -1.0;
  }

  // Try both capitalizations (older tasks use Zvertex, newer use zvertex)
  TH1 *hZvtx = (TH1 *)file->Get(Form("%s/h_collisions_Zvertex", Dir));
  if (!hZvtx) hZvtx = (TH1 *)file->Get(Form("%s/h_collisions_zvertex", Dir));
  if (!hZvtx) {
    std::cerr << "[Zvtx] No h_collisions_Zvertex/zvertex in " << Dir << std::endl;
    file->Close();
    return -1.0;
  }
  hZvtx->SetDirectory(nullptr);

  // Fit Gaussian: mean initialized to histogram mean, free to float
  double histMean = hZvtx->GetMean();
  TF1 *fitFunc = new TF1("zvtxFit", "gaus", -10, 10);
  fitFunc->SetParameter(1, histMean);
  hZvtx->Fit(fitFunc, "RQ0");

  // Compute ratio: Int[-10,10] / Int[-inf,+inf] using fitted Gaussian
  double num = fitFunc->Integral(-10, 10);
  double den = fitFunc->Integral(-1e6, 1e6);  // effectively (-inf, +inf) for Gaussian
  double eff = (den > 0) ? (num / den) : -1.0;

  double fittedMean = fitFunc->GetParameter(1);
  double sigma = fitFunc->GetParameter(2);
  std::cerr << "[Zvtx] " << Dir << ": Gaussian mean=" << Form("%.3f", fittedMean)
            << " cm, sigma=" << Form("%.3f", sigma)
            << " cm, zvtx10 eff=" << Form("%.5f", eff) << std::endl;

  // Draw QA plot if output directory is provided
  if (plotDir && plotLabel) {
    gSystem->mkdir(plotDir, true);
    TString label = plotLabel;

    TCanvas *cZvtx = new TCanvas(Form("cZvtx_%s", label.Data()), "", 700, 500);
    cZvtx->SetGrid();

    hZvtx->SetStats(0);
    hZvtx->SetTitle("");
    hZvtx->GetXaxis()->SetTitle("z_{vtx} (cm)");
    hZvtx->GetYaxis()->SetTitle("Events");
    hZvtx->GetXaxis()->SetRangeUser(-20, 20);
    hZvtx->SetMarkerStyle(20);
    hZvtx->SetMarkerSize(0.5);
    hZvtx->SetLineColor(kBlack);
    hZvtx->Draw("PE");

    // Draw full Gaussian fit extended beyond ±10
    TF1 *fitFull = new TF1("zvtxFitFull", "gaus", -20, 20);
    fitFull->SetParameters(fitFunc->GetParameter(0), fitFunc->GetParameter(1), sigma);
    fitFull->SetLineColor(kRed);
    fitFull->SetLineWidth(2);
    fitFull->Draw("same");

    // Shade the |z| < 10 cm region
    TF1 *fitShade = new TF1("zvtxFitShade", "gaus", -10, 10);
    fitShade->SetParameters(fitFunc->GetParameter(0), fitFunc->GetParameter(1), sigma);
    fitShade->SetLineColor(kRed);
    fitShade->SetLineWidth(0);
    fitShade->SetFillColorAlpha(kRed, 0.2);
    fitShade->Draw("FC same");

    // |z| = 10 cm boundary lines
    Double_t yMax = hZvtx->GetMaximum();
    TLine *lineL = new TLine(-10, 0, -10, yMax * 0.8);
    TLine *lineR = new TLine(10, 0, 10, yMax * 0.8);
    lineL->SetLineColor(kBlue); lineL->SetLineWidth(2); lineL->SetLineStyle(2);
    lineR->SetLineColor(kBlue); lineR->SetLineWidth(2); lineR->SetLineStyle(2);
    lineL->Draw();
    lineR->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.15, 0.85, Form("#mu_{z} = %.3f cm", fittedMean));
    latex.DrawLatex(0.15, 0.80, Form("#sigma_{z} = %.3f cm", sigma));
    latex.DrawLatex(0.15, 0.75, Form("#varepsilon_{|z|<10} = %.5f", eff));
    latex.DrawLatex(0.15, 0.70, label.Data());

    TLegend *leg = new TLegend(0.55, 0.75, 0.92, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0, 0);
    leg->AddEntry(hZvtx, "Data z_{vtx}", "pe");
    leg->AddEntry(fitFull, "Gaussian fit", "l");
    leg->AddEntry(fitShade, "|z_{vtx}| < 10 cm", "f");
    leg->Draw();

    cZvtx->Print(Form("%s/ZvtxEfficiency_%s.pdf", plotDir, label.Data()));
    delete cZvtx;
    delete fitFull;
    delete fitShade;
  }

  delete fitFunc;
  file->Close();
  return eff;
}

// ============================================
// Unified cross-section normalization
// ============================================
// Two modes:
//   kNormSoft (default):
//     L_correct = hLumiTVXafterBCcutsRCT × (σ_vis_CCDB / σ_vis_correct)
//     normFactor = 1 / L_correct[mb⁻¹] / ε_zvtx10
//     - hLumiTVXafterBCcutsRCT is BC-level (no collision reconstruction bias)
//     - Already contains: TVX trigger, sel8 BC cuts, RCT, framework pile-up correction
//     - σ_vis_CCDB estimated from data: nTVX / hLumiTVX (≈ σ_vis_CCDB/<pileupCorr>)
//     - No ε_INEL→TVX: for jets, TVX captures all INEL events with jets
//       (SD/DD events without TVX have no hard scattering → no jets)
//     - No ε_sel8: collision-level sel8 cuts (vertex quality, ITS-TPC vertex)
//       are ~100% efficient for jet events (many tracks → good vertex)
//
//   kNormDirect:
//     L_eff = N_TVX × pileupCorr / σ_vis × ε_sel8 × ε_RCT
//     normFactor = 1 / L_eff[mb⁻¹] / ε_zvtx10
//     - Uses hardcoded σ_vis and conservative pile-up (μ=0.05)
//     - Counter-based sel8/RCT efficiencies
//     - Same physics: no ε_INEL→TVX needed for jets
//
// All efficiency routes are printed as diagnostics regardless of mode.

struct XsecNormResult {
  Double_t normFactor;       // final: multiply unfolded spectrum by this
  Double_t lumiEff;          // L_correct [μb⁻¹]
  Double_t epsZvtx10;        // |z| < 10 cm fraction
  Double_t sigmaVisCCDB_est; // estimated σ_vis_CCDB from data [mb] (kNormSoft only)
  bool valid;
};

inline XsecNormResult ComputeXsecNormFactor(
    const NormDecomposition& decomp,
    Double_t epsZvtx10,
    int datasetYear,
    EXsecNormMode normMode)
{
  // Year-dependent constants
  const Double_t kSigmaVisTVX2022  = 53.0;   // mb (temporary, no precise reference yet)
  const Double_t kSigmaVisTVX2023  = 53.4;   // mb
  // Reference ε_INEL→TVX from CrossSectionEfficiency.cxx (diagnostic only, NOT applied)
  const Double_t kEpsINELtoTVX2022 = 0.9661;
  const Double_t kEpsINELtoTVX2023 = 0.9227;
  // Pile-up correction for kNormDirect: μ/(1-exp(-μ)), conservative μ=0.05
  const Double_t kMuPileup = 0.05;
  const Double_t kPileupCorr = kMuPileup / (1.0 - TMath::Exp(-kMuPileup));  // ≈ 1.025

  Double_t sigmaVis = (datasetYear >= 2023) ? kSigmaVisTVX2023 : kSigmaVisTVX2022;

  XsecNormResult r = {};
  r.epsZvtx10 = epsZvtx10;
  r.sigmaVisCCDB_est = 0;
  r.valid = false;

  if (normMode == kNormSoft) {
    // ---- kNormSoft: rescale framework luminosity for correct σ_vis ----
    // hLumiTVXafterBCcutsRCT is BC-level: TVX trigger + sel8 BC cuts + RCT + pile-up correction.
    // No collision reconstruction bias. Already contains everything except zvtx.
    //
    // Problem: uses σ_vis_CCDB which may be outdated (e.g., ~59 mb instead of 53.4 mb).
    // Solution: estimate σ_vis_CCDB from data, then rescale.
    //
    // hLumiTVX = Σ (pileupCorr_i / σ_vis_CCDB)  [μb⁻¹]
    // hCounterTVX = N_TVX  (raw BC count, no pile-up)
    // σ_vis_CCDB_est ≈ nTVX / hLumiTVX  (≈ σ_vis_CCDB / <pileupCorr>, ~2% bias)
    //
    // L_correct = hLumiTVXafterBCcutsRCT × (σ_vis_CCDB_est / σ_vis_correct)
    //           = hLumiTVXafterBCcutsRCT × (nTVX / hLumiTVX) / (σ_vis_correct × 1e3)
    //
    // Note: the ~2% pile-up in σ_vis_CCDB_est approximately cancels, but the
    // luminosity-weighted efficiency structure (per-run RCT quality) is preserved.
    //
    // No ε_INEL→TVX: TVX captures all NSD events, and no jets exist in SD events.
    // No ε_sel8_coll: collision-level sel8 cuts are ~100% efficient for jet events.
    if (decomp.lumiTVX > 0 && decomp.nTVX > 0) {
      // σ_vis_CCDB_est [μb] = nTVX / lumiTVX [μb⁻¹];  convert to mb: / 1e3
      r.sigmaVisCCDB_est = decomp.nTVX / decomp.lumiTVX / 1e3;  // mb
      Double_t rescale = r.sigmaVisCCDB_est / sigmaVis;
      r.lumiEff = decomp.lumiTVXafterBCRCT * rescale;
    } else {
      std::cerr << "[XsecNorm] Warning: hLumiTVX or hCounterTVX missing, cannot estimate sigma_vis_CCDB" << std::endl;
      r.lumiEff = decomp.lumiTVXafterBCRCT;  // fallback: use as-is (wrong σ_vis)
    }
  } else {
    // ---- kNormDirect: counter-based with hardcoded pile-up ----
    // L_TVX = N_TVX × pileupCorr / σ_vis  (BC-level, conservative μ=0.05)
    // L_eff = L_TVX × ε_sel8 × ε_RCT  (counter-based efficiencies)
    Double_t lumiTVX_ub = decomp.nTVX * kPileupCorr / sigmaVis / 1e3;
    r.lumiEff = lumiTVX_ub * decomp.epsSel8 * decomp.epsRCT;
  }

  if (r.lumiEff > 0 && r.epsZvtx10 > 0) {
    Double_t lumiEffMb = r.lumiEff * 1e3;  // μb⁻¹ → mb⁻¹
    r.normFactor = 1.0 / lumiEffMb / r.epsZvtx10;
    r.valid = true;
  }

  // ---- Diagnostics: print exact equation and all parameters ----
  Double_t epsINELtoTVX_ref = (datasetYear >= 2023) ? kEpsINELtoTVX2023 : kEpsINELtoTVX2022;
  std::cerr << "============================================" << std::endl;
  std::cerr << "[XsecNorm] Mode: " << (normMode == kNormSoft ? "SOFT" : "DIRECT") << std::endl;
  std::cerr << "[XsecNorm] Dataset year: " << datasetYear
            << "  sigma_vis_correct = " << Form("%.1f", sigmaVis) << " mb" << std::endl;
  std::cerr << std::endl;

  if (normMode == kNormSoft) {
    std::cerr << "[XsecNorm] Equation (kNormSoft):" << std::endl;
    std::cerr << "[XsecNorm]   dsigma/dpT = N_jets * normFactor / Delta_eta / Delta_pT" << std::endl;
    std::cerr << "[XsecNorm]   normFactor = 1 / L_correct[mb] / eps_zvtx" << std::endl;
    std::cerr << "[XsecNorm]   L_correct  = hLumiTVXafterBCcutsRCT * (sigma_vis_CCDB / sigma_vis_correct)" << std::endl;
    std::cerr << "[XsecNorm]   sigma_vis_CCDB estimated from hCounterTVX / hLumiTVX" << std::endl;
    std::cerr << "[XsecNorm]   (BC-level: no collision reco bias; includes pile-up correction + sel8 BC + RCT)" << std::endl;
    std::cerr << "[XsecNorm]   (No eps_INEL->TVX: TVX captures all jet-producing INEL events)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "[XsecNorm] Parameters:" << std::endl;
    std::cerr << "[XsecNorm]   hLumiTVXafterBCcutsRCT = " << Form("%.2f", decomp.lumiTVXafterBCRCT) << " ub^-1" << std::endl;
    std::cerr << "[XsecNorm]   hLumiTVX               = " << Form("%.2f", decomp.lumiTVX) << " ub^-1" << std::endl;
    std::cerr << "[XsecNorm]   hCounterTVX (N_TVX)    = " << Form("%.0f", decomp.nTVX) << std::endl;
    std::cerr << "[XsecNorm]   sigma_vis_CCDB (est)   = " << Form("%.1f", r.sigmaVisCCDB_est) << " mb"
              << "  (= N_TVX/hLumiTVX, includes ~2%% pileup bias)" << std::endl;
    std::cerr << "[XsecNorm]   sigma_vis_correct      = " << Form("%.1f", sigmaVis) << " mb" << std::endl;
    std::cerr << "[XsecNorm]   rescale factor          = " << Form("%.4f", r.sigmaVisCCDB_est / sigmaVis)
              << "  (sigma_vis_CCDB / sigma_vis_correct)" << std::endl;
    std::cerr << "[XsecNorm]   L_correct              = " << Form("%.2f", r.lumiEff) << " ub^-1 = "
              << Form("%.2f", r.lumiEff*1e3) << " mb^-1" << std::endl;
    std::cerr << "[XsecNorm]   Runs: " << decomp.nRunsGood << " good / " << decomp.nRunsTotal << " total" << std::endl;
  } else {
    Double_t lumiTVX_mb = decomp.nTVX * kPileupCorr / sigmaVis;
    std::cerr << "[XsecNorm] Equation (kNormDirect):" << std::endl;
    std::cerr << "[XsecNorm]   dsigma/dpT = N_jets * normFactor / Delta_eta / Delta_pT" << std::endl;
    std::cerr << "[XsecNorm]   normFactor = 1 / L_eff[mb] / eps_zvtx" << std::endl;
    std::cerr << "[XsecNorm]   L_eff = L_TVX * eps_sel8 * eps_RCT" << std::endl;
    std::cerr << "[XsecNorm]   L_TVX = N_TVX * pileupCorr / sigma_vis" << std::endl;
    std::cerr << "[XsecNorm]   pileupCorr = mu/(1-exp(-mu)), mu = " << Form("%.3f", kMuPileup) << std::endl;
    std::cerr << std::endl;
    std::cerr << "[XsecNorm] Parameters:" << std::endl;
    std::cerr << "[XsecNorm]   N_TVX (hCounterTVX) = " << Form("%.0f", decomp.nTVX) << std::endl;
    std::cerr << "[XsecNorm]   sigma_vis           = " << Form("%.1f", sigmaVis) << " mb" << std::endl;
    std::cerr << "[XsecNorm]   pileupCorr (mu=" << Form("%.3f", kMuPileup) << ") = " << Form("%.4f", kPileupCorr) << std::endl;
    std::cerr << "[XsecNorm]   L_TVX               = " << Form("%.2f", lumiTVX_mb) << " mb^-1 = "
              << Form("%.2f", lumiTVX_mb/1e3) << " ub^-1" << std::endl;
    std::cerr << "[XsecNorm]   eps_sel8 (counter)   = " << Form("%.5f", decomp.epsSel8) << std::endl;
    std::cerr << "[XsecNorm]   eps_RCT  (counter)   = " << Form("%.5f", decomp.epsRCT)
              << "  (" << decomp.nRunsGood << "/" << decomp.nRunsTotal << " runs pass)" << std::endl;
    std::cerr << "[XsecNorm]   L_eff               = " << Form("%.2f", r.lumiEff*1e3) << " mb^-1 = "
              << Form("%.2f", r.lumiEff) << " ub^-1" << std::endl;
  }
  std::cerr << "[XsecNorm]   eps_zvtx10          = " << Form("%.5f", r.epsZvtx10) << std::endl;
  std::cerr << std::endl;
  if (r.valid) {
    std::cerr << "[XsecNorm] Result:" << std::endl;
    std::cerr << "[XsecNorm]   normFactor = " << Form("%.4e", r.normFactor) << std::endl;
    std::cerr << "[XsecNorm]   = 1 / " << Form("%.2f", r.lumiEff*1e3) << " mb^-1"
              << " / " << Form("%.5f", r.epsZvtx10) << std::endl;
  } else {
    std::cerr << "[XsecNorm] ** FAILED ** -- normFactor not computed" << std::endl;
  }
  std::cerr << std::endl;

  // ---- Reference: eps_INEL->TVX (NOT applied, diagnostic only) ----
  std::cerr << "  --- Reference (NOT applied): eps_INEL->TVX = " << Form("%.4f", epsINELtoTVX_ref)
            << " (from CrossSectionEfficiency MC)" << std::endl;
  std::cerr << "  --- Omitted because: TVX captures all jet-producing INEL events" << std::endl;
  std::cerr << "  --- (includes ~" << Form("%.1f", (1-epsINELtoTVX_ref)*100)
            << "%% from reco+sel8 collision cuts, ~100%% efficient for jets)" << std::endl;

  // Diagnostic: compare sel8/RCT efficiency from counter vs luminosity routes
  std::cerr << "  --- Cross-check: sel8/RCT efficiency routes ---" << std::endl;
  std::cerr << "  Counter-based:  eps_sel8=" << Form("%.5f", decomp.epsSel8)
            << "  eps_RCT=" << Form("%.5f", decomp.epsRCT)
            << "  product=" << Form("%.5f", decomp.epsSel8RCT) << std::endl;
  std::cerr << "  Lumi-based:     eps_sel8=" << Form("%.5f", decomp.lumiEpsSel8)
            << "  eps_RCT=" << Form("%.5f", decomp.lumiEpsRCT)
            << "  product=" << Form("%.5f", decomp.lumiEpsSel8 * decomp.lumiEpsRCT) << std::endl;
  std::cerr << "  L_eff check:    L_TVXafterBCRCT/L_TVX = "
            << Form("%.5f", decomp.lumiTVX > 0 ? decomp.lumiTVXafterBCRCT / decomp.lumiTVX : 0)
            << "  (should ~ counter product)" << std::endl;
  std::cerr << "============================================" << std::endl;

  return r;
}

// Find optimal SVD k using d-vector method + stability validation.
// Step 1: d-vector finds the signal/noise boundary (standard ALICE method,
//         Höcker & Kartvelishvili, NIM A372, 469, 1996).
//         d_i = (U^T b~)_i: data projected onto SVD basis.
//         Signal components: |d_i| >> 1; Noise: |d_i| ~ sqrt(2/pi) ~ 0.8.
//         Consecutive-below criterion: two consecutive |d_i| < 1.0.
// Step 2: Starting from the d-vector boundary, find the first k where BOTH
//         |y(k-1)/y(k)-1| < ε AND |y(k)/y(k+1)-1| < ε (bilateral stability).
//         This ensures k is on the convergence plateau (not at the knee).
inline Int_t FindOptimalSvdK_Dvector(RooUnfoldResponse *Response, TH1 *hRecoData,
                                      Double_t stabilityThreshold = 0.03,
                                      const char *plotLabel = nullptr,
                                      const char *plotDir = nullptr,
                                      Double_t chi2Threshold = 2.0) {
  if (!Response || !hRecoData)
    return 6;

  const Int_t nBins = hRecoData->GetNbinsX();
  const Int_t maxK = TMath::Min(nBins, 25);
  const Double_t dvecThreshold = 1.0;

  std::cerr << "[Info] SVD k optimization: d-vector (upper bound) + refolding chi2 + bilateral stability" << std::endl;

  // Suppress TDecompSVD convergence warnings during k-scan
  Int_t savedErrorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;

  // --- Step 1: Unfold at each k, compute refolding chi2 ---
  std::vector<TH1*> results;
  std::vector<double> chi2Values;

  for (Int_t k = 1; k <= maxK; ++k) {
    RooUnfoldSvd unfold(Response, hRecoData, k);
    TH1* h = (TH1*)unfold.Hreco()->Clone(Form("hSVD_stab_k%d", k));
    h->SetDirectory(nullptr);
    results.push_back(h);

    // Refolding chi2
    TH1* hRefold = (TH1*)Response->ApplyToTruth(h, Form("hRefold_SVD_k%d", k));
    double chi2 = 0.0;
    int ndf = 0;
    if (hRefold) {
      hRefold->SetDirectory(nullptr);
      for (int i = 1; i <= hRecoData->GetNbinsX(); ++i) {
        double ptLow = hRecoData->GetXaxis()->GetBinLowEdge(i);
        double ptHigh = hRecoData->GetXaxis()->GetBinUpEdge(i);
        if (ptLow < GetPlotPtMin() || ptHigh > GetPlotPtMax()) continue;
        const double measured = hRecoData->GetBinContent(i);
        const double measuredErr = hRecoData->GetBinError(i);
        const double refolded = hRefold->GetBinContent(i);
        if (measuredErr > 0 && measured > 0) {
          double pull = (measured - refolded) / measuredErr;
          chi2 += pull * pull;
          ndf++;
        }
      }
      delete hRefold;
    }
    double reducedChi2 = (ndf > 0) ? chi2 / ndf : 1e10;
    chi2Values.push_back(reducedChi2);

    std::cerr << "[Info]   SVD k=" << k
              << ": refolding chi2/ndf = " << reducedChi2 << " (ndf=" << ndf << ")"
              << (reducedChi2 < chi2Threshold ? "  [ok]" : "") << std::endl;
  }

  // --- Step 2: d-vector (upper bound) ---
  RooUnfoldSvd unfoldFull(Response, hRecoData, nBins);
  unfoldFull.Hreco();

  TSVDUnfold_local *svdImpl = unfoldFull.Impl();
  Int_t dvectorK = maxK;
  TH1D *hD = nullptr;
  if (svdImpl) {
    hD = svdImpl->GetD();
    if (hD) {
      const Int_t nD = hD->GetNbinsX();
      std::cerr << "[Info] d-vector values:" << std::endl;
      for (Int_t i = 1; i <= nD; ++i) {
        Double_t absDi = TMath::Abs(hD->GetBinContent(i));
        std::cerr << "[Info]   d_" << i << " = " << hD->GetBinContent(i)
                  << "  |d_" << i << "| = " << absDi
                  << (absDi > dvecThreshold ? "  > 1 (signal)" : "  <= 1 (noise)") << std::endl;
      }
      // Robust boundary: require 3+ out of 4 consecutive values below threshold
      // (avoids premature capping from isolated noise pockets, e.g. R=0.2 UE-sub)
      const Int_t windowSize = 4;
      const Int_t noiseRequired = 3;
      for (Int_t i = 1; i <= nD - windowSize + 1; ++i) {
        Int_t noiseCount = 0;
        for (Int_t j = 0; j < windowSize; ++j) {
          if (TMath::Abs(hD->GetBinContent(i + j)) < dvecThreshold)
            noiseCount++;
        }
        if (noiseCount >= noiseRequired) {
          dvectorK = i - 1;
          break;
        }
      }
      if (dvectorK < 1) dvectorK = 1;
    }
  }
  std::cerr << "[Info] d-vector upper bound: k=" << dvectorK << std::endl;

  // --- Step 3: Find first k where chi2/ndf < threshold ---
  Int_t chi2K = -1;
  for (int i = 0; i < maxK; ++i) {
    if (chi2Values[i] < chi2Threshold) {
      chi2K = i + 1;
      break;
    }
  }

  Int_t candidateK;
  if (chi2K > 0) {
    candidateK = chi2K;
    std::cerr << "[Info] Refolding chi2 converged at k=" << chi2K << std::endl;
  } else {
    candidateK = 2;
    std::cerr << "[Warning] Refolding chi2/ndf never dropped below " << chi2Threshold
              << ", starting bilateral stability from k=2." << std::endl;
  }

  // Cap by d-vector: don't search past signal/noise boundary
  Int_t searchMax = TMath::Min(dvectorK + 1, maxK);
  std::cerr << "[Info] Bilateral stability scan range: k=" << candidateK
            << " to " << searchMax << " (d-vector cap)" << std::endl;

  // --- Step 4: Bilateral stability within [candidateK, searchMax] ---
  std::cerr << "[Info] Bilateral stability validation (threshold = "
            << stabilityThreshold * 100 << "%):" << std::endl;

  auto maxStepChange = [&](Int_t kA, Int_t kB) -> Double_t {
    Int_t idxA = kA - 1, idxB = kB - 1;
    if (idxA < 0 || idxA >= maxK || idxB < 0 || idxB >= maxK) return 999.;
    Double_t maxRel = 0.0;
    for (Int_t bin = 1; bin <= results[idxA]->GetNbinsX(); ++bin) {
      Double_t ptLow = results[idxA]->GetXaxis()->GetBinLowEdge(bin);
      Double_t ptHigh = results[idxA]->GetXaxis()->GetBinUpEdge(bin);
      if (ptLow < GetPlotPtMin() || ptHigh > GetPlotPtMax()) continue;
      Double_t a = results[idxA]->GetBinContent(bin);
      Double_t b = results[idxB]->GetBinContent(bin);
      if (a > 0 && b > 0) {
        Double_t rel = TMath::Abs(b / a - 1.0);
        if (rel > maxRel) maxRel = rel;
      }
    }
    return maxRel;
  };

  // Collect stability data for plotting (full range for QA)
  std::vector<Double_t> backSteps, fwdSteps;
  for (Int_t k = 2; k < maxK; ++k) {
    backSteps.push_back(maxStepChange(k - 1, k));
    fwdSteps.push_back(maxStepChange(k, k + 1));
  }

  // Scan for bilateral stability plateau: find ALL stable k, pick the one with best chi2.
  // This avoids stopping too early (e.g. R=0.2 UE-sub: stable at k=5 but chi2 keeps improving).
  Int_t firstStableK = -1;
  Int_t bestStableK = -1;
  Double_t bestChi2 = 1e10;
  for (Int_t k = candidateK; k < searchMax; ++k) {
    Double_t backStep = maxStepChange(k - 1, k);
    Double_t fwdStep  = maxStepChange(k, k + 1);
    Bool_t backOK = (backStep < stabilityThreshold);
    Bool_t fwdOK  = (fwdStep  < stabilityThreshold);

    std::cerr << "[Info]   k=" << k
              << ": k-1->k = " << backStep * 100 << "%" << (backOK ? " [ok]" : "")
              << ",  k->k+1 = " << fwdStep * 100 << "%" << (fwdOK ? " [ok]" : "")
              << (backOK && fwdOK ? "  => STABLE" : "") << std::endl;

    if (backOK && fwdOK) {
      if (firstStableK < 0) firstStableK = k;
      // Among all stable k, pick the one with best refolding chi2
      if (k - 1 < (Int_t)chi2Values.size() && chi2Values[k - 1] < bestChi2) {
        bestChi2 = chi2Values[k - 1];
        bestStableK = k;
      }
    }
  }

  Int_t finalK = candidateK;
  Bool_t foundStable = (bestStableK > 0);
  if (foundStable) {
    finalK = bestStableK;
  }

  if (!foundStable) {
    if (chi2K > 0) {
      finalK = chi2K;
      std::cerr << "[Warning] No bilateral stability within d-vector range. "
                << "Using chi2 convergence k=" << chi2K << std::endl;
    } else {
      std::cerr << "[Warning] SVD optimization: neither chi2 convergence nor bilateral stability found. "
                << "Using k=" << finalK << ". Consider relaxing thresholds." << std::endl;
    }
  } else {
    std::cerr << "[Info] Optimal k=" << finalK
              << " (first stable at k=" << firstStableK
              << ", best chi2=" << Form("%.3f", bestChi2)
              << ", chi2 converged at k=" << chi2K
              << ", d-vector cap at k=" << dvectorK
              << ", stable at k=" << finalK << ")" << std::endl;
  }

  if (finalK < 2) finalK = 2;
  if (finalK > 25) finalK = 25;

  std::cerr << "[Info] Optimal SVD k (chi2 + d-vector + stability): " << finalK << std::endl;

  // Restore ROOT error level after k-scan
  gErrorIgnoreLevel = savedErrorLevel;

  // --- QA Plots ---
  if (plotDir && plotLabel) {
    gSystem->mkdir(plotDir, true);

    // (A) D-vector plot: |d_i| vs i with threshold line
    if (hD) {
      const Int_t nD = hD->GetNbinsX();
      TCanvas *cDvec = new TCanvas(Form("cDvec_%s", plotLabel), "", 700, 500);
      cDvec->SetLogy();
      cDvec->SetGrid();

      TH1D *hAbsD = new TH1D(Form("hAbsD_%s", plotLabel), Form("|d_{i}| vs i  (%s);Singular value index i;|d_{i}|", plotLabel),
                              nD, 0.5, nD + 0.5);
      hAbsD->SetDirectory(nullptr);
      for (Int_t i = 1; i <= nD; ++i) {
        hAbsD->SetBinContent(i, TMath::Abs(hD->GetBinContent(i)));
      }
      hAbsD->SetMarkerStyle(20);
      hAbsD->SetMarkerSize(1.0);
      hAbsD->SetLineWidth(2);
      hAbsD->SetMinimum(0.05);
      hAbsD->SetMaximum(hAbsD->GetMaximum() * 5);
      hAbsD->Draw("P");

      TLine *line1 = new TLine(0.5, 1.0, nD + 0.5, 1.0);
      line1->SetLineColor(kRed);  line1->SetLineWidth(2);  line1->SetLineStyle(2);
      line1->Draw();

      TLine *lineK = new TLine(finalK + 0.5, hAbsD->GetMinimum(), finalK + 0.5, hAbsD->GetMaximum());
      lineK->SetLineColor(kBlue);  lineK->SetLineWidth(2);  lineK->SetLineStyle(1);
      lineK->Draw();

      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.04);
      latex.DrawLatex(0.55, 0.85, Form("Optimal k = %d", finalK));
      latex.DrawLatex(0.55, 0.80, Form("d-vector boundary = %d", dvectorK));
      latex.DrawLatex(0.55, 0.75, Form("#chi^{2} converged at k = %d", chi2K));

      TLegend *leg = new TLegend(0.55, 0.58, 0.88, 0.73);
      leg->SetBorderSize(0);  leg->SetFillColorAlpha(0, 0);
      leg->AddEntry(hAbsD, "|d_{i}| from data", "p");
      leg->AddEntry(line1, "|d_{i}| = 1 threshold", "l");
      leg->AddEntry(lineK, Form("k = %d (chosen)", finalK), "l");
      leg->Draw();

      cDvec->Print(Form("%s/SVD_Dvector_%s.pdf", plotDir, plotLabel));
      delete cDvec;
    }

    // (B) Refolding chi2/ndf vs k (start from k=2 to avoid k=1 dominating y-range)
    {
      TCanvas *cChi2 = new TCanvas(Form("cChi2SVD_%s", plotLabel), "", 700, 500);
      cChi2->SetGrid();

      // Skip k=1 (huge chi2 compresses the interesting region)
      Int_t nPtsK = maxK - 1;  // k=2..maxK
      std::vector<Double_t> kVals;
      std::vector<Double_t> chi2Vals;
      for (int i = 1; i < maxK; ++i) {
        kVals.push_back(i + 1);  // k=2,3,...,maxK
        chi2Vals.push_back(chi2Values[i]);
      }

      TGraph *gChi2 = new TGraph(nPtsK, kVals.data(), chi2Vals.data());
      gChi2->SetMarkerStyle(20);  gChi2->SetMarkerSize(1.0);
      gChi2->SetLineWidth(2);
      gChi2->SetTitle(Form("SVD refolding #chi^{2}/ndf (%s);k;#chi^{2}/ndf", plotLabel));
      Double_t yMaxChi2 = 0;
      for (int i = 0; i < nPtsK; ++i) yMaxChi2 = TMath::Max(yMaxChi2, chi2Vals[i]);
      gChi2->GetYaxis()->SetRangeUser(0, TMath::Max(yMaxChi2 * 1.2, chi2Threshold * 3));
      gChi2->Draw("ALP");

      TLine *lineThr = new TLine(2, chi2Threshold, maxK, chi2Threshold);
      lineThr->SetLineColor(kRed);  lineThr->SetLineWidth(2);  lineThr->SetLineStyle(2);
      lineThr->Draw();

      TLine *lineK = new TLine(finalK, 0, finalK, gChi2->GetYaxis()->GetXmax() * 0.9);
      lineK->SetLineColor(kBlue);  lineK->SetLineWidth(2);  lineK->SetLineStyle(1);
      lineK->Draw();

      TLegend *leg = new TLegend(0.50, 0.70, 0.88, 0.88);
      leg->SetBorderSize(0);  leg->SetFillColorAlpha(0, 0);
      leg->AddEntry(gChi2, "#chi^{2}/ndf", "lp");
      leg->AddEntry(lineThr, Form("threshold = %.1f", chi2Threshold), "l");
      leg->AddEntry(lineK, Form("k = %d (chosen)", finalK), "l");
      leg->Draw();

      cChi2->Print(Form("%s/SVD_Chi2_%s.pdf", plotDir, plotLabel));
      delete cChi2;
    }

    // (C) Bilateral stability: max |y(k±1)/y(k) - 1| vs k
    {
      TCanvas *cStab = new TCanvas(Form("cStab_%s", plotLabel), "", 700, 500);
      cStab->SetLogy();
      cStab->SetGrid();

      Int_t nPts = (Int_t)backSteps.size();
      TH1D *hBack = new TH1D(Form("hBack_%s", plotLabel),
        Form("Bilateral stability (%s);k;max |y(k#pm1)/y(k) #minus 1|", plotLabel),
        nPts, 1.5, nPts + 1.5);
      TH1D *hFwd = new TH1D(Form("hFwd_%s", plotLabel), "", nPts, 1.5, nPts + 1.5);
      hBack->SetDirectory(nullptr);
      hFwd->SetDirectory(nullptr);
      for (Int_t i = 0; i < nPts; ++i) {
        hBack->SetBinContent(i + 1, backSteps[i]);
        hFwd->SetBinContent(i + 1, fwdSteps[i]);
      }
      hBack->SetMarkerStyle(20); hBack->SetMarkerColor(kBlue); hBack->SetLineColor(kBlue);
      hFwd->SetMarkerStyle(24);  hFwd->SetMarkerColor(kRed);  hFwd->SetLineColor(kRed);
      hBack->SetMinimum(1e-4);
      hBack->SetMaximum(1.0);
      hBack->Draw("P");
      hFwd->Draw("P same");

      TLine *lineThresh = new TLine(1.5, stabilityThreshold, nPts + 1.5, stabilityThreshold);
      lineThresh->SetLineColor(kGreen + 2);  lineThresh->SetLineWidth(2);  lineThresh->SetLineStyle(2);
      lineThresh->Draw();

      TLine *lineK = new TLine(finalK, hBack->GetMinimum(), finalK, hBack->GetMaximum());
      lineK->SetLineColor(kBlue);  lineK->SetLineWidth(2);
      lineK->Draw();

      // Show d-vector cap as dashed vertical line
      if (dvectorK < maxK) {
        TLine *lineDvec = new TLine(dvectorK, hBack->GetMinimum(), dvectorK, hBack->GetMaximum());
        lineDvec->SetLineColor(kGray + 2);  lineDvec->SetLineWidth(2);  lineDvec->SetLineStyle(3);
        lineDvec->Draw();
      }

      TLegend *leg = new TLegend(0.45, 0.68, 0.88, 0.88);
      leg->SetBorderSize(0);  leg->SetFillColorAlpha(0, 0);
      leg->AddEntry(hBack, "|y(k#minus1)/y(k) #minus 1|", "p");
      leg->AddEntry(hFwd, "|y(k)/y(k+1) #minus 1|", "p");
      leg->AddEntry(lineThresh, Form("%.0f%% threshold", stabilityThreshold * 100), "l");
      leg->Draw();

      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.04);
      latex.DrawLatex(0.50, 0.63, Form("Optimal k = %d", finalK));

      cStab->Print(Form("%s/SVD_BilateralStability_%s.pdf", plotDir, plotLabel));
      delete cStab;
    }
  }

  // Cleanup
  for (auto* h : results) delete h;

  return finalK;
}

// Plot Bayesian refolding chi2/ndf vs iteration QA (self-contained: unfolds internally)
inline void PlotBayesChi2Diagnostic(RooUnfoldResponse *Response, TH1 *hRecoData,
                                     Int_t optimalIter, Double_t chi2Threshold,
                                     Int_t maxIter,
                                     const char *plotLabel = nullptr,
                                     const char *plotDir = nullptr) {
  if (!plotDir || !plotLabel) return;
  const Int_t nBins = hRecoData->GetNbinsX();
  std::vector<Double_t> iterVals, chi2Vals;

  for (int iter = 1; iter <= maxIter; ++iter) {
    RooUnfoldBayes unfold(Response, hRecoData, iter);
    TH1 *hUnf = (TH1 *)unfold.Hreco()->Clone(Form("hBayes_chi2plot_%s_iter%d", plotLabel, iter));
    hUnf->SetDirectory(nullptr);
    TH1 *hRefold = (TH1 *)Response->ApplyToTruth(hUnf, Form("hRefold_chi2plot_%s_iter%d", plotLabel, iter));
    double chi2 = 0.0;
    int ndf = 0;
    if (hRefold) {
      hRefold->SetDirectory(nullptr);
      for (int i = 1; i <= nBins; ++i) {
        double ptLow = hRecoData->GetXaxis()->GetBinLowEdge(i);
        double ptHigh = hRecoData->GetXaxis()->GetBinUpEdge(i);
        if (ptLow < GetPlotPtMin() || ptHigh > GetPlotPtMax()) continue;
        double measured = hRecoData->GetBinContent(i);
        double measuredErr = hRecoData->GetBinError(i);
        double refolded = hRefold->GetBinContent(i);
        if (measuredErr > 0 && measured > 0) {
          double pull = (measured - refolded) / measuredErr;
          chi2 += pull * pull;
          ndf++;
        }
      }
      delete hRefold;
    }
    delete hUnf;
    double reducedChi2 = (ndf > 0) ? chi2 / ndf : 0.0;
    iterVals.push_back(iter);
    chi2Vals.push_back(reducedChi2);
  }

  if (iterVals.empty()) return;

  TGraph *gChi2 = new TGraph(iterVals.size(), iterVals.data(), chi2Vals.data());
  gChi2->SetMarkerStyle(20);  gChi2->SetMarkerSize(1.0);
  gChi2->SetMarkerColor(kBlack);  gChi2->SetLineColor(kBlack);  gChi2->SetLineWidth(2);

  TCanvas *c = new TCanvas(Form("cBayesChi2_%s", plotLabel), "Bayes chi2/ndf", 600, 500);
  c->SetGridy();

  Double_t yMax = 0;
  for (auto v : chi2Vals) yMax = TMath::Max(yMax, v);
  yMax = TMath::Max(yMax * 1.3, chi2Threshold * 3.0);

  gChi2->SetTitle("");
  gChi2->GetXaxis()->SetTitle("Bayesian iteration");
  gChi2->GetYaxis()->SetTitle("Refolding #chi^{2}/ndf");
  gChi2->GetXaxis()->CenterTitle();  gChi2->GetYaxis()->CenterTitle();
  gChi2->GetYaxis()->SetRangeUser(0, yMax);
  gChi2->GetXaxis()->SetLimits(0.5, maxIter + 0.5);
  gChi2->Draw("ALP");

  TLine *lineThresh = new TLine(0.5, chi2Threshold, maxIter + 0.5, chi2Threshold);
  lineThresh->SetLineColor(kOrange + 1);
  lineThresh->SetLineStyle(7);  lineThresh->SetLineWidth(2);
  lineThresh->Draw("same");

  TLine *lineOne = new TLine(0.5, 1.0, maxIter + 0.5, 1.0);
  lineOne->SetLineColor(kRed);
  lineOne->SetLineStyle(7);  lineOne->SetLineWidth(2);
  lineOne->Draw("same");

  TLine *lineOpt = new TLine(optimalIter, 0, optimalIter, yMax * 0.9);
  lineOpt->SetLineColor(kBlue);
  lineOpt->SetLineStyle(7);  lineOpt->SetLineWidth(2);
  lineOpt->Draw("same");

  TLegend *leg = new TLegend(0.45, 0.65, 0.88, 0.88);
  leg->SetBorderSize(0);  leg->SetFillStyle(0);  leg->SetTextSize(0.035);
  leg->AddEntry(gChi2, "Refolding #chi^{2}/ndf", "lp");
  leg->AddEntry(lineThresh, Form("Threshold = %.1f", chi2Threshold), "l");
  leg->AddEntry(lineOne, "#chi^{2}/ndf = 1", "l");
  leg->AddEntry(lineOpt, Form("Optimal iter = %d", optimalIter), "l");
  leg->Draw();

  gSystem->mkdir(plotDir, kTRUE);
  c->SaveAs(Form("%s/Bayes_Chi2_%s.pdf", plotDir, plotLabel));
  std::cerr << "[Info] Bayes chi2/ndf plot saved to " << plotDir << "/Bayes_Chi2_" << plotLabel << ".pdf" << std::endl;
  delete c;
}

// Plot Bayesian bilateral stability QA (self-contained: unfolds internally)
inline void PlotBayesBilateralStability(RooUnfoldResponse *Response, TH1 *hRecoData,
                                         Int_t optimalIter, Double_t threshold,
                                         const char *plotLabel = nullptr,
                                         const char *plotDir = nullptr) {
  if (!plotDir || !plotLabel) return;
  const Int_t nBins = hRecoData->GetNbinsX();
  const Int_t maxIter = TMath::Min(nBins, 20);
  if (maxIter < 3) return;

  std::vector<TH1*> results;
  for (Int_t iter = 1; iter <= maxIter; ++iter) {
    RooUnfoldBayes unfold(Response, hRecoData, iter);
    TH1* h = (TH1*)unfold.Hreco()->Clone(Form("hBayesStabPlot_%s_iter%d", plotLabel, iter));
    h->SetDirectory(nullptr);
    results.push_back(h);
  }

  auto maxStepChange = [&](Int_t iterA, Int_t iterB) -> Double_t {
    Int_t idxA = iterA - 1, idxB = iterB - 1;
    if (idxA < 0 || idxA >= maxIter || idxB < 0 || idxB >= maxIter) return -1.;
    Double_t maxRel = 0.0;
    for (Int_t bin = 1; bin <= results[idxA]->GetNbinsX(); ++bin) {
      Double_t ptLow = results[idxA]->GetXaxis()->GetBinLowEdge(bin);
      Double_t ptHigh = results[idxA]->GetXaxis()->GetBinUpEdge(bin);
      if (ptLow < GetPlotPtMin() || ptHigh > GetPlotPtMax()) continue;
      Double_t a = results[idxA]->GetBinContent(bin);
      Double_t b = results[idxB]->GetBinContent(bin);
      if (a > 0 && b > 0) {
        Double_t rel = TMath::Abs(b / a - 1.0);
        if (rel > maxRel) maxRel = rel;
      }
    }
    return maxRel;
  };

  std::vector<Double_t> iterVals, fwdVals, backVals;
  for (Int_t it = 2; it < maxIter; ++it) {
    Double_t fwd  = maxStepChange(it, it + 1);
    Double_t back = maxStepChange(it - 1, it);
    if (fwd < 0 || back < 0) continue;
    iterVals.push_back(it);
    fwdVals.push_back(fwd * 100.0);
    backVals.push_back(back * 100.0);
  }

  for (auto* h : results) delete h;
  if (iterVals.empty()) return;

  TGraph *gFwd  = new TGraph(iterVals.size(), iterVals.data(), fwdVals.data());
  TGraph *gBack = new TGraph(iterVals.size(), iterVals.data(), backVals.data());
  gFwd->SetMarkerStyle(20);  gFwd->SetMarkerSize(1.0);
  gFwd->SetMarkerColor(kRed);  gFwd->SetLineColor(kRed);  gFwd->SetLineWidth(2);
  gBack->SetMarkerStyle(21); gBack->SetMarkerSize(1.0);
  gBack->SetMarkerColor(kBlue); gBack->SetLineColor(kBlue); gBack->SetLineWidth(2);

  TCanvas *c = new TCanvas(Form("cBayesStab_%s", plotLabel), "Bayes Bilateral Stability", 700, 500);
  c->SetGridy();

  Double_t yMax = 0;
  for (size_t i = 0; i < fwdVals.size(); ++i)
    yMax = TMath::Max(yMax, TMath::Max(fwdVals[i], backVals[i]));
  yMax = TMath::Max(yMax * 1.3, threshold * 100 * 2.0);

  gBack->SetTitle("");
  gBack->GetXaxis()->SetTitle("Bayesian iteration");
  gBack->GetYaxis()->SetTitle("max |#it{y}(#it{n}#pm1)/#it{y}(#it{n}) #minus 1| (%)");
  gBack->GetXaxis()->CenterTitle();  gBack->GetYaxis()->CenterTitle();
  gBack->GetYaxis()->SetRangeUser(0, yMax);
  gBack->GetXaxis()->SetLimits(iterVals.front() - 0.5, iterVals.back() + 0.5);
  gBack->Draw("ALP");
  gFwd->Draw("LP same");

  TLine *lineThresh = new TLine(iterVals.front() - 0.5, threshold * 100,
                                 iterVals.back() + 0.5, threshold * 100);
  lineThresh->SetLineColor(kGreen + 2);
  lineThresh->SetLineStyle(7);  lineThresh->SetLineWidth(2);
  lineThresh->Draw("same");

  TLine *lineOpt = new TLine(optimalIter, 0, optimalIter, yMax * 0.9);
  lineOpt->SetLineColor(kBlack);
  lineOpt->SetLineStyle(2);  lineOpt->SetLineWidth(2);
  lineOpt->Draw("same");

  TLegend *leg = new TLegend(0.45, 0.65, 0.88, 0.88);
  leg->SetBorderSize(0);  leg->SetFillStyle(0);  leg->SetTextSize(0.035);
  leg->AddEntry(gBack, "max |#it{y}(#it{n}#minus1)/#it{y}(#it{n}) #minus 1|", "lp");
  leg->AddEntry(gFwd,  "max |#it{y}(#it{n})/#it{y}(#it{n}+1) #minus 1|", "lp");
  leg->AddEntry(lineThresh, Form("#varepsilon = %.0f%%", threshold * 100), "l");
  leg->AddEntry(lineOpt, Form("Optimal iter = %d", optimalIter), "l");
  leg->Draw();

  gSystem->mkdir(plotDir, kTRUE);
  c->SaveAs(Form("%s/Bayes_BilateralStability_%s.pdf", plotDir, plotLabel));
  std::cerr << "[Info] Bayes bilateral stability plot saved to " << plotDir
            << "/Bayes_BilateralStability_" << plotLabel << ".pdf" << std::endl;
  delete c;
}

// Find optimal Bayesian iteration using refolding chi2/ndf + bilateral stability.
// Analogous to SVD d-vector + bilateral stability:
//   Step 1: Refolding chi2/ndf sets lower bound (first iter where chi2/ndf < threshold).
//   Step 2: Bilateral stability finds the convergence plateau.
inline Int_t FindOptimalBayesIter(RooUnfoldResponse *Response, TH1 *hRecoData,
                                   Int_t maxIter = 20,
                                   Double_t chi2Threshold = 2.0,
                                   Double_t stabilityThreshold = 0.03,
                                   const char *plotLabel = nullptr,
                                   const char *plotDir = nullptr) {
  std::cerr << "[Info] Bayesian optimization: refolding chi2 + bilateral stability..." << std::endl;

  // Unfold at each iteration and store results
  std::vector<TH1*> results;
  std::vector<double> chi2Values;

  for (int iter = 1; iter <= maxIter; ++iter) {
    RooUnfoldBayes unfoldBayes(Response, hRecoData, iter);
    TH1 *hUnf = (TH1 *)unfoldBayes.Hreco()->Clone(Form("hBayes_opt_%s_iter%d", plotLabel ? plotLabel : "x", iter));
    hUnf->SetDirectory(nullptr);
    results.push_back(hUnf);

    // Refolding chi2
    TH1 *hRefold = (TH1 *)Response->ApplyToTruth(hUnf, Form("hRefold_Bayes_opt_%s_iter%d", plotLabel ? plotLabel : "x", iter));
    double chi2 = 0.0;
    int ndf = 0;
    if (hRefold) {
      hRefold->SetDirectory(nullptr);
      for (int i = 1; i <= hRecoData->GetNbinsX(); ++i) {
        double ptLow = hRecoData->GetXaxis()->GetBinLowEdge(i);
        double ptHigh = hRecoData->GetXaxis()->GetBinUpEdge(i);
        if (ptLow < GetPlotPtMin() || ptHigh > GetPlotPtMax()) continue;
        const double measured = hRecoData->GetBinContent(i);
        const double measuredErr = hRecoData->GetBinError(i);
        const double refolded = hRefold->GetBinContent(i);
        if (measuredErr > 0 && measured > 0) {
          double pull = (measured - refolded) / measuredErr;
          chi2 += pull * pull;
          ndf++;
        }
      }
      delete hRefold;
    }
    double reducedChi2 = (ndf > 0) ? chi2 / ndf : 1e10;
    chi2Values.push_back(reducedChi2);

    std::cerr << "[Info]   Bayes iter " << iter
              << ": refolding chi2/ndf = " << reducedChi2 << " (ndf=" << ndf << ")"
              << (reducedChi2 < chi2Threshold ? "  [ok]" : "")
              << std::endl;
  }

  // Step 1: Find first iteration where chi2/ndf < threshold
  Int_t chi2Iter = -1;
  for (int i = 0; i < maxIter; ++i) {
    if (chi2Values[i] < chi2Threshold) {
      chi2Iter = i + 1;
      break;
    }
  }

  Int_t candidateIter = chi2Iter;
  if (chi2Iter < 0) {
    std::cerr << "[Warning] Refolding chi2/ndf never dropped below " << chi2Threshold
              << ", starting bilateral stability from iter=2." << std::endl;
    candidateIter = 2;
  }
  std::cerr << "[Info] chi2/ndf boundary: iter=" << (chi2Iter > 0 ? chi2Iter : -1)
            << ", starting bilateral scan from iter=" << candidateIter << std::endl;

  // Step 2: Bilateral stability (same criterion as SVD)
  std::cerr << "[Info] Bilateral stability validation (threshold = "
            << stabilityThreshold * 100 << "%):" << std::endl;

  auto maxStepChange = [&](Int_t iterA, Int_t iterB) -> Double_t {
    Int_t idxA = iterA - 1, idxB = iterB - 1;
    if (idxA < 0 || idxA >= maxIter || idxB < 0 || idxB >= maxIter) return 999.;
    Double_t maxRel = 0.0;
    for (Int_t bin = 1; bin <= results[idxA]->GetNbinsX(); ++bin) {
      Double_t ptLow = results[idxA]->GetXaxis()->GetBinLowEdge(bin);
      Double_t ptHigh = results[idxA]->GetXaxis()->GetBinUpEdge(bin);
      if (ptLow < GetPlotPtMin() || ptHigh > GetPlotPtMax()) continue;
      Double_t a = results[idxA]->GetBinContent(bin);
      Double_t b = results[idxB]->GetBinContent(bin);
      if (a > 0 && b > 0) {
        Double_t rel = TMath::Abs(b / a - 1.0);
        if (rel > maxRel) maxRel = rel;
      }
    }
    return maxRel;
  };

  // Scan for bilateral stability plateau: find ALL stable iter, pick the one with best chi2.
  // This avoids stopping too early (same fix as SVD: R=0.2 UE-sub first stable at low iter
  // but chi2 keeps improving at higher iter).
  Int_t firstStableIter = -1;
  Int_t bestStableIter = -1;
  Double_t bestChi2 = 1e10;
  for (Int_t it = candidateIter; it < maxIter; ++it) {
    Double_t backStep = maxStepChange(it - 1, it);
    Double_t fwdStep  = maxStepChange(it, it + 1);
    Bool_t backOK = (backStep < stabilityThreshold);
    Bool_t fwdOK  = (fwdStep  < stabilityThreshold);

    std::cerr << "[Info]   iter=" << it
              << ": iter-1->iter = " << backStep * 100 << "%" << (backOK ? " [ok]" : "")
              << ",  iter->iter+1 = " << fwdStep * 100 << "%" << (fwdOK ? " [ok]" : "")
              << (backOK && fwdOK ? "  => STABLE" : "") << std::endl;

    if (backOK && fwdOK) {
      if (firstStableIter < 0) firstStableIter = it;
      // Among all stable iter, pick the one with best refolding chi2
      if (it - 1 < (Int_t)chi2Values.size() && chi2Values[it - 1] < bestChi2) {
        bestChi2 = chi2Values[it - 1];
        bestStableIter = it;
      }
    }
  }

  Int_t finalIter = candidateIter;
  Bool_t foundStable = (bestStableIter > 0);
  if (foundStable) {
    finalIter = bestStableIter;
  }

  if (!foundStable) {
    std::cerr << "[Warning] Bayesian bilateral stability never satisfied (threshold "
              << stabilityThreshold * 100 << "%). Using iter=" << finalIter
              << ". Consider increasing maxIter or relaxing threshold." << std::endl;
  } else {
    std::cerr << "[Info] Bayesian optimal iter=" << finalIter
              << " (first stable=" << firstStableIter
              << ", best chi2=" << bestChi2
              << " at iter=" << bestStableIter << ")" << std::endl;
  }

  // Cleanup
  for (auto* h : results) delete h;

  if (finalIter < 2) finalIter = 2;
  if (finalIter > maxIter) finalIter = maxIter;

  std::cerr << "[Info] Optimal Bayes iter (chi2 + bilateral stability): " << finalIter << std::endl;

  // QA plots (if output location provided)
  if (plotLabel && plotDir) {
    PlotBayesChi2Diagnostic(Response, hRecoData, finalIter, chi2Threshold, maxIter, plotLabel, plotDir);
    PlotBayesBilateralStability(Response, hRecoData, finalIter, stabilityThreshold, plotLabel, plotDir);
  }

  return finalIter;
}

// Refolding test: unfold data, fold back through response, compare to measured
inline void PlotRefoldingTest(RooUnfoldResponse *Response, TH1 *hRecoData,
                               TH1 *hUnfolded, Int_t optK,
                               const char *plotLabel, const char *plotDir,
                               Double_t nEvents = 1.0) {
  if (!Response || !hRecoData || !hUnfolded || !plotDir || !plotLabel) return;
  gSystem->mkdir(plotDir, true);

  TH1 *hRefold = (TH1 *)Response->ApplyToTruth(hUnfolded, Form("hRefold_%s", plotLabel));
  if (!hRefold) return;
  hRefold->SetDirectory(nullptr);

  // Upper pad: measured vs refolded
  TCanvas *cRef = new TCanvas(Form("cRefold_%s", plotLabel), "", 700, 700);
  TPad *padUp = new TPad("padUp", "", 0, 0.35, 1, 1.0);
  padUp->SetBottomMargin(0.02);
  padUp->SetLogy();
  padUp->SetTicks(1, 1);
  padUp->Draw();
  TPad *padLo = new TPad("padLo", "", 0, 0.0, 1, 0.35);
  padLo->SetTopMargin(0.02);
  padLo->SetBottomMargin(0.3);
  padLo->SetTicks(1, 1);
  padLo->Draw();

  padUp->cd();
  TH1 *hMeas = (TH1 *)hRecoData->Clone(Form("hMeas_%s", plotLabel));
  hMeas->SetDirectory(nullptr);
  // Normalize by events and bin width for proper density display
  if (nEvents > 0) hMeas->Scale(1.0 / nEvents, "width");
  hMeas->SetTitle(Form("Refolding test (%s)", plotLabel));
  hMeas->GetYaxis()->SetTitle("1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
  hMeas->GetXaxis()->SetLabelSize(0);
  hMeas->SetMarkerStyle(20);
  hMeas->SetMarkerSize(0.8);
  hMeas->SetLineColor(kBlack);
  hMeas->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
  hMeas->Draw("pe");

  TH1 *hRef = (TH1 *)hRefold->Clone(Form("hRefDraw_%s", plotLabel));
  hRef->SetDirectory(nullptr);
  // Same normalization as measured
  if (nEvents > 0) hRef->Scale(1.0 / nEvents, "width");
  hRef->SetLineColor(kRed);
  hRef->SetLineWidth(2);
  hRef->SetMarkerStyle(24);
  hRef->SetMarkerColor(kRed);
  hRef->SetMarkerSize(0.8);
  hRef->Draw("pe same");

  TLegend *leg = new TLegend(0.50, 0.70, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);
  leg->AddEntry(hMeas, "Measured (purity-corrected)", "pe");
  leg->AddEntry(hRef, Form("Refolded (SVD k=%d)", optK), "pe");
  leg->Draw();

  // Lower pad: ratio
  padLo->cd();
  TH1 *hRat = (TH1 *)hRefold->Clone(Form("hRefRat_%s", plotLabel));
  hRat->SetDirectory(nullptr);
  hRat->Divide(hRecoData);
  hRat->SetTitle("");
  hRat->GetYaxis()->SetTitle("Refolded / Measured");
  hRat->GetXaxis()->SetTitle("#it{p}_{T,jet}^{ch} (GeV/#it{c})");
  hRat->GetYaxis()->SetRangeUser(0.95, 1.05);
  hRat->GetYaxis()->SetNdivisions(505);
  hRat->GetYaxis()->SetTitleSize(0.10);
  hRat->GetYaxis()->SetTitleOffset(0.5);
  hRat->GetYaxis()->SetLabelSize(0.08);
  hRat->GetXaxis()->SetTitleSize(0.10);
  hRat->GetXaxis()->SetLabelSize(0.08);
  hRat->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
  hRat->SetLineColor(kRed);
  hRat->SetMarkerColor(kRed);
  hRat->SetMarkerStyle(20);
  hRat->SetMarkerSize(0.8);
  hRat->Draw("pe");

  TLine *unity = new TLine(5, 1.0, 200, 1.0);
  unity->SetLineColor(kBlack);
  unity->SetLineStyle(2);
  unity->Draw();

  cRef->Print(Form("%s/RefoldingTest_%s.pdf", plotDir, plotLabel));
  delete cRef;
  delete hRefold;
}

template <typename T>
inline void hset(T &hid, TString xtit = "", TString ytit = "", double titoffx = 0.9,
          double titoffy = 1.2, double titsizex = 0.06, double titsizey = 0.06,
          double labeloffx = 0.01, double labeloffy = 0.001,
          double labelsizex = 0.05, double labelsizey = 0.05, int divx = 510,
          int divy = 510) {
  hid.GetXaxis()->CenterTitle(1);
  hid.GetYaxis()->CenterTitle(1);
  hid.GetXaxis()->SetTitle(xtit);
  hid.GetYaxis()->SetTitle(ytit);
  hid.GetXaxis()->SetTitleOffset(titoffx);
  hid.GetYaxis()->SetTitleOffset(titoffy);
  hid.GetXaxis()->SetTitleSize(titsizex);
  hid.GetYaxis()->SetTitleSize(titsizey);
  hid.GetXaxis()->SetLabelOffset(labeloffx);
  hid.GetYaxis()->SetLabelOffset(labeloffy);
  hid.GetXaxis()->SetLabelSize(labelsizex);
  hid.GetYaxis()->SetLabelSize(labelsizey);
  hid.GetXaxis()->SetNdivisions(divx);
  hid.GetYaxis()->SetNdivisions(divy);
}

template <typename T>
inline void hoptset(T &hid, Double_t Nevts, Color_t colorID, Double_t Xmin, Double_t Xmax,
           Double_t Ymin, Double_t Ymax, Double_t MarkerSize = 1.0, Int_t LineStyle = 1,
           Int_t LineWidth = 2, Int_t MarkerStyle = 20) {
  // Only scale if Nevts > 0 to avoid division by zero (for efficiency/purity plots, Nevts=0 means no normalization needed)
  if (Nevts > 0) {
    hid.Scale(1.0 / Nevts, "width");
  }
  hid.SetLineColor(colorID);
  hid.SetMarkerColor(colorID);
  hid.SetMarkerStyle(MarkerStyle);
  hid.SetMarkerSize(MarkerSize);
  hid.SetLineStyle(LineStyle);
  hid.SetLineWidth(LineWidth);
  hid.GetXaxis()->SetRangeUser(Xmin, Xmax);
  hid.GetYaxis()->SetRangeUser(Ymin, Ymax);
  hid.SetStats(0);  // Disable stat box
}

inline void optFili(TPad &pad, Int_t logx = 0, Int_t logy = 0, Int_t logz = 0, Int_t grid = 0) {
  pad.SetLogx(logx);
  pad.SetLogy(logy);
  pad.SetLogz(logz);
  pad.SetGrid(grid);
}

inline void setpad(TVirtualPad *pad, Double_t tmargin = 0.02, Double_t bmargin = 0.15,
            Double_t lmargin = 0.13, Double_t rmargin = 0.05) {
  pad->SetTopMargin(tmargin);
  pad->SetLeftMargin(lmargin);
  pad->SetRightMargin(rmargin);
  pad->SetBottomMargin(bmargin);
}

inline TH1 *DrawRatioTH1(TH1 *hNum, TH1 *hDenom) {
  auto axN = hNum->GetXaxis();
  auto axD = hDenom->GetXaxis();
  const int nN = axN->GetNbins();
  const int nD = axD->GetNbins();
  
  auto edgesEqual = [&](int nA, TAxis* a, int nB, TAxis* b) -> bool {
    if (nA != nB) return false;
    for (int i = 1; i <= nA; ++i) {
      double loA = a->GetBinLowEdge(i);
      double hiA = a->GetBinUpEdge(i);
      double loB = b->GetBinLowEdge(i);
      double hiB = b->GetBinUpEdge(i);
      if (TMath::Abs(loA - loB) > 1e-6 || TMath::Abs(hiA - hiB) > 1e-6) return false;
    }
    return true;
  };
  
  TH1 *hRatio = nullptr;
  if (edgesEqual(nN, axN, nD, axD)) {
    hRatio = (TH1 *)hDenom->Clone("hRatio");
    hRatio->Reset();
    hRatio->SetStats(0);  // Disable stat box
    for (int i = 1; i <= nN; ++i) {
      double yNum = hNum->GetBinContent(i);
      double eNum = hNum->GetBinError(i);
      double yDen = hDenom->GetBinContent(i);
      double eDen = hDenom->GetBinError(i);
      if (TMath::Abs(yDen) > 1e-10) {  // Avoid division by zero or very small numbers
        double ratio = yNum / yDen;
        // Proper error propagation: sqrt((eNum/yDen)^2 + (yNum*eDen/yDen^2)^2)
        double term1 = eNum / yDen;
        double term2 = (yNum * eDen) / (yDen * yDen);
        double ratioErr = TMath::Sqrt(term1 * term1 + term2 * term2);
        hRatio->SetBinContent(i, ratio);
        hRatio->SetBinError(i, ratioErr);
      } else {
        hRatio->SetBinContent(i, 0);
        hRatio->SetBinError(i, 0);
      }
    }
  } else {
    hRatio = (TH1 *)hDenom->Clone("hRatio");
    hRatio->Reset();
    hRatio->SetStats(0);  // Disable stat box
    for (int i = 1; i <= hDenom->GetNbinsX(); ++i) {
      double xCenter = hDenom->GetXaxis()->GetBinCenter(i);
      int binNum = hNum->GetXaxis()->FindBin(xCenter);
      if (binNum >= 1 && binNum <= hNum->GetNbinsX()) {
        double yNum = hNum->GetBinContent(binNum);
        double eNum = hNum->GetBinError(binNum);
        double yDen = hDenom->GetBinContent(i);
        double eDen = hDenom->GetBinError(i);
        if (TMath::Abs(yDen) > 1e-10) {  // Avoid division by zero or very small numbers
          double ratio = yNum / yDen;
          // Proper error propagation: sqrt((eNum/yDen)^2 + (yNum*eDen/yDen^2)^2)
          double term1 = eNum / yDen;
          double term2 = (yNum * eDen) / (yDen * yDen);
          double ratioErr = TMath::Sqrt(term1 * term1 + term2 * term2);
          hRatio->SetBinContent(i, ratio);
          hRatio->SetBinError(i, ratioErr);
        } else {
          hRatio->SetBinContent(i, 0);
          hRatio->SetBinError(i, 0);
        }
      }
    }
  }
  return hRatio;
}

inline TH1 *DrawRatio(const char *ratioName, TH1 *refHist, TH1 *testHist,
               TString AxisTitleX, TString AxisTitleY, Color_t colorID,
               Double_t Ymin, Double_t Ymax, Double_t MarkerSize = 1) {
  TH1 *ratioHist = DrawRatioTH1(testHist, refHist);
  hset(*ratioHist, AxisTitleX, AxisTitleY, 1.2, 1.0, 0.07, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 505);
  ratioHist->GetYaxis()->SetRangeUser(Ymin, Ymax);
  ratioHist->SetMarkerColor(colorID);
  ratioHist->SetLineColor(colorID);
  ratioHist->SetMarkerStyle(testHist->GetMarkerStyle());  // follow numerator
  ratioHist->SetMarkerSize(MarkerSize);
  ratioHist->SetStats(0);  // Disable stat box
  ratioHist->Draw("esame");
  return ratioHist;
}

void ALICEfigureLegend(const char *FigureLabel, double x11, double y11, double x12, double y12, double x21, double y21, double x22, double y22, double textsize=0.043, double Rvalue = -1) {
    TLegend *ALICEleg2 =
        new TLegend(x21, y21, x22, y22, NULL,"brNDC");
    ALICEleg2->SetTextSize(textsize);
    ALICEleg2->SetBorderSize(0);
    ALICEleg2->SetTextAlign(12);
    ALICEleg2->SetFillColorAlpha(0,0); 
    ALICEleg2->AddEntry("", "Anti-#it{k}_{T}, charged-particle jets", "");
    if (Rvalue != -1) {
    ALICEleg2->AddEntry("", Form("#it{R}=%.1f, |#it{#eta}_{jet}| < %.1f", Rvalue, GetJetEtaMax(Rvalue)), "");
    } else {
    ALICEleg2->AddEntry("", (GetCurrentConfigSet().jetEtaCut >= 0)
      ? Form("|#it{#eta}_{jet}| < %.1f", GetCurrentConfigSet().jetEtaCut)
      : "|#it{#eta}_{jet}| < 0.9-#it{R}", "");
    }
    TLegend *ALICEleg1 =
        new TLegend(x11, y11, x12, y12, NULL,"brNDC");
    ALICEleg1->SetTextSize(textsize);
    ALICEleg1->SetBorderSize(0);
    ALICEleg1->SetTextAlign(12);
    ALICEleg1->SetFillColorAlpha(0,0); 
    ALICEleg1->AddEntry("", FigureLabel, "");
    ALICEleg1->AddEntry("", "pp #sqrt{#it{s}} = 13.6 TeV", "");
    ALICEleg1->AddEntry("", "#it{p}_{T, track} > 0.15 GeV/#it{c}", "");
    ALICEleg1->AddEntry("", "|#it{#eta}_{track}| < 0.9", "");

    ALICEleg1->Draw();
    ALICEleg2->Draw();
}

// Drawing functions - simplified versions
inline TH1 *DrawTrackPt(const char *fileName, const char *histName, Double_t Nevts,
                 TLegend *legend, Color_t colorID, const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "OPEN");
  if (!file || file->IsZombie()) return nullptr;
  
  TH1* TrackPt = (TH1*) file->Get(Form("%s/%s", Dir, GetTrackPtObj()));
  if (!TrackPt) {
    file->Close();
    return nullptr;
  }
  
  // Clone and detach from file before rebinning
  TrackPt = (TH1*)TrackPt->Clone(Form("TrackPt_temp_%s", histName));
  TrackPt->SetDirectory(0);
  
  if (GetREBINON()) {
    TH1* rebinned = TrackPt->Rebin(GetNTrackptbin(), Form("TrackPt_%s", histName), GetTrackptbin());
    if (rebinned) {
      rebinned->SetDirectory(0);
      delete TrackPt;
      TrackPt = rebinned;
    }
  }
  
  legend->AddEntry(TrackPt, histName);
  hset(*TrackPt, GetTrackPtTitleX(), "1/#it{N}_{evt} d#it{N}/d#it{p}_{T}", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
  if (GetNORMEVENTS()) {
    hoptset(*TrackPt, Nevts, colorID, 0, 200, 2e-10, 1e0 + 0.05, 1, 1, 2, colorID==kBlack? 21 : 20);
  } else {
    hoptset(*TrackPt, Nevts, colorID, 0, 200, 2e-10, 1e0 + 0.05, 0.6, 1, 2, 24);
  }
  TrackPt->Draw("pesame");
  file->Close();
  return TrackPt;
}

inline TH1 *DrawTrackEta(const char *fileName, const char *histName, Double_t Nevts,
                  TLegend *legend, Color_t colorID, const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  if (!file || file->IsZombie()) return nullptr;
  
  TH2 *h2TrackEta = (TH2 *)file->Get(Form("%s/%s", Dir, GetTrackEtaObj()));
  if (!h2TrackEta) {
    file->Close();
    return nullptr;
  }
  
  TH1* TrackEta = (TH1*) h2TrackEta->ProjectionX(Form("TrackEta_%s", histName), 1, h2TrackEta->GetNbinsY(), "e");
  if (!TrackEta) {
    file->Close();
    return nullptr;
  }
  
  TrackEta->SetDirectory(0);
  // Rebin by 5: ensures bin edges at multiples of 0.1, aligning with track |eta|<0.9 boundary
  TH1 *h1TrackEta = TrackEta->Rebin(5, Form("TrackEta_rebinned_%s", histName));
  if (h1TrackEta) {
    h1TrackEta->SetDirectory(0);
    delete TrackEta;
  } else {
    h1TrackEta = TrackEta;
  }
  
  legend->AddEntry(h1TrackEta, histName);
  hset(*h1TrackEta, GetTrackEtaTitleX(), "1/#it{N}_{evt} d#it{N}/d#it{#eta}", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
  
  if (GetNORMEVENTS()) {
    hoptset(*h1TrackEta, Nevts, colorID, -0.89, 0.89, 0, 16);
  } else {
    hoptset(*h1TrackEta, Nevts, colorID, -0.89, 0.89, 0., 1.5);
  }
  h1TrackEta->Draw("esame");
  file->Close();
  return h1TrackEta;
}

inline TH1 *DrawTrackPhi(const char *fileName, const char *histName, Double_t Nevts,
                  TLegend *legend, Color_t colorID, const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  if (!file || file->IsZombie()) return nullptr;
  
  TH2 *h2TrackPhi = (TH2 *)file->Get(Form("%s/%s", Dir, GetTrackPhiObj()));
  if (!h2TrackPhi) {
    file->Close();
    return nullptr;
  }
  
  TH1* TrackPhi = (TH1*) h2TrackPhi->ProjectionY(Form("TrackPhi_%s", histName), 1, h2TrackPhi->GetNbinsX(), "e");
  if (!TrackPhi) {
    file->Close();
    return nullptr;
  }
  
  TrackPhi->SetDirectory(0);
  legend->AddEntry(TrackPhi, histName);
  hset(*TrackPhi, GetTrackPhiTitleX(), "1/#it{N}_{evt} d#it{N}/d#it{#varphi}", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
  
  // Display range clipped to [0.01, 6.19] to avoid partially populated edge bins
  // (phiAxis = {160, -1, 7}, physics range [0, 2pi~6.28], last full bin edge ~6.2)
  if (GetNORMEVENTS()) {
    hoptset(*TrackPhi, Nevts, colorID, 0.01, 6.19, 0, 5);
  } else {
    hoptset(*TrackPhi, Nevts, colorID, 0.01, 6.19, 0, 0.3);
  }
  TrackPhi->Draw("esame");
  file->Close();
  return TrackPhi;
}

inline TH1 *DrawConstituentPt(const char *fileName, const char *histName, Double_t Nevts,
                       TLegend *legend, Color_t colorID, const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  if (!file || file->IsZombie()) return nullptr;
  
  TH2 *h2ConstituentPt = (TH2 *)file->Get(Form("%s/%s", Dir, GetConstituentPtObj()));
  if (!h2ConstituentPt) {
    file->Close();
    return nullptr;
  }
  
  // Project Y-axis (track pT) for all jet pT bins
  TH1* ConstituentPt = (TH1*) h2ConstituentPt->ProjectionY(Form("ConstituentPt_%s", histName), 1, h2ConstituentPt->GetNbinsX(), "e");
  if (!ConstituentPt) {
    file->Close();
    return nullptr;
  }
  
  ConstituentPt->SetDirectory(0);
  
  // Rebin if needed (same as DrawTrackPt)
  if (GetREBINON()) {
    TH1* rebinned = ConstituentPt->Rebin(GetNTrackptbin(), Form("ConstituentPt_%s", histName), GetTrackptbin());
    if (rebinned) {
      rebinned->SetDirectory(0);
      delete ConstituentPt;
      ConstituentPt = rebinned;
    }
  }
  
  legend->AddEntry(ConstituentPt, histName);
  hset(*ConstituentPt, "#it{p}_{T, constituent}^{reco} (GeV/#it{c})", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T}", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
  if (GetNORMEVENTS()) {
    hoptset(*ConstituentPt, Nevts, colorID, 0, 99, 2e-10, 1e0 + 0.05, 1, 1, 2, colorID==kBlack? 21 : 20);
  } else {
    hoptset(*ConstituentPt, Nevts, colorID, 0, 99, 2e-10, 1e0 + 0.05, 0.6, 1, 2, 24);
  }
  ConstituentPt->Draw("pesame");
  file->Close();
  return ConstituentPt;
}

inline TH1 *DrawNtracksDistribution(const char *fileName, const char *histName, Double_t Nevts,
                             TLegend *legend, Color_t colorID, const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  if (!file || file->IsZombie()) return nullptr;
  
  TH2 *h2Ntracks = (TH2 *)file->Get(Form("%s/%s", Dir, GetNtracksObj()));
  if (!h2Ntracks) {
    file->Close();
    return nullptr;
  }
  
  // Project Y-axis (Ntracks) for all jet pT bins
  TH1* Ntracks = (TH1*) h2Ntracks->ProjectionY(Form("Ntracks_%s", histName), 1, h2Ntracks->GetNbinsX(), "e");
  if (!Ntracks) {
    file->Close();
    return nullptr;
  }
  
  Ntracks->SetDirectory(0);
  
  legend->AddEntry(Ntracks, histName);
  hset(*Ntracks, "N_{tracks}", "1/#it{N}_{evt} d#it{N}/dN_{tracks}", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
  if (GetNORMEVENTS()) {
    hoptset(*Ntracks, Nevts, colorID, 0, 40, 1e-9, 1e-1, 1, 1, 2, colorID==kBlack? 21 : 20);
  } else {
    hoptset(*Ntracks, Nevts, colorID, 0, 40, 1e-9, 1e-1, 0.6, 1, 2, 24);
  }
  Ntracks->Draw("pesame");
  file->Close();
  return Ntracks;
}

// Helper function for Y-axis range
inline std::pair<double, double> getYAxisRange(TH1* histogram, Double_t Nevts, double Xmin, double Xmax, double UpMarginFactor = 5, double DownMarginFactor = 15, double extraScale = 1.0) {
    if (!histogram || Nevts <= 0) {
        return {0, 1};
    }
    int binLow = histogram->FindBin(Xmin);
    int binHigh = histogram->FindBin(Xmax);
    if (binLow < 1) binLow = 1;
    if (binHigh > histogram->GetNbinsX()) binHigh = histogram->GetNbinsX();
    double widthLow = histogram->GetBinWidth(binLow);
    double widthHigh = histogram->GetBinWidth(binHigh);
    // Estimate Y values after Scale(1/Nevts, "width") and optional extraScale (e.g. 1/deltaEta)
    double yMax = histogram->GetBinContent(binLow) / Nevts / widthLow * extraScale;
    double yMin = histogram->GetBinContent(binHigh) / Nevts / widthHigh * extraScale;
    yMin = yMin / DownMarginFactor;
    yMax = yMax * UpMarginFactor;
    return {yMin, yMax};
}

// Jet drawing functions - implementations
inline TH1 *DrawJetPt(const char *fileName, const char *histName, const char *Obj, Double_t Nevts,
               TLegend *legend, Color_t colorID, Int_t i = 0, const char *Dir = nullptr, Double_t R = 0.4) {
  auto file = TFile::Open(fileName, "open");
  if (!file || file->IsZombie()) return nullptr;
  TH3 *h3JetPt = (TH3 *)file->Get(Form("%s/%s", Dir, Obj));
  if (!h3JetPt) {
    file->Close();
    return nullptr;
  }
  // Apply eta cut: fixed |eta|<jetEtaCut or full range (R-dependent already applied at production)
  int etaBinLo = 1, etaBinHi = h3JetPt->GetNbinsY();
  double cut = GetCurrentConfigSet().jetEtaCut;
  if (cut >= 0) {
    etaBinLo = h3JetPt->GetYaxis()->FindBin(-cut + 1e-6);
    etaBinHi = h3JetPt->GetYaxis()->FindBin( cut - 1e-6);
  }
  TH1 *JetPt = h3JetPt->ProjectionX(Form("JetPt_%s", histName), etaBinLo, etaBinHi, 1, h3JetPt->GetNbinsZ());
  if (!JetPt) {
    file->Close();
    return nullptr;
  }
  JetPt->SetDirectory(0);
  if (GetREBINON()) {
    PtBinning recoBin = DetectRecoBinning(JetPt);
    TH1* rebinned = RebinToTarget(JetPt, recoBin, Form("JetPt_%s", histName));
    if (rebinned) {
      delete JetPt;
      JetPt = rebinned;
    }
  }
  legend->AddEntry(JetPt, histName);
  // Apply dEta normalization (R-dependent or fixed, from config)
  Double_t deltaEta = GetDeltaEta(R);
  hset(*JetPt, GetJetPtTitleX(), Form("1/#it{N}_{evt} d#it{N}/(d#it{p}_{T}d#it{#eta})"), 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
  double extraScale = (deltaEta > 0) ? (1.0 / deltaEta) : 1.0;
  auto [yMin, yMax] = getYAxisRange(JetPt, Nevts, GetPlotPtMin(), GetPlotPtMax(), 5, 15, extraScale);
  hoptset(*JetPt, Nevts, colorID, GetPlotPtMin(), GetPlotPtMax(), yMin, yMax, colorID==kBlack? 1.2 : 1, 1, 2, colorID==kBlack? 21 : 20);
  // Normalize by dEta
  if (deltaEta > 0) {
    JetPt->Scale(1.0 / deltaEta);
  }
  JetPt->Draw("esame");
  file->Close();
  return JetPt;
}

inline TH1 *DrawJetEta(const char *fileName, const char *histName, const char *Obj, Double_t Nevts,
                TLegend *legend, Color_t colorID, const char *Dir = nullptr, Double_t R = 0.4) {
  auto file = TFile::Open(fileName, "open");
  if (!file || file->IsZombie()) return nullptr;
  TH3D *H3JetEta = (TH3D *)file->Get(Form("%s/%s", Dir, Obj));
  if (!H3JetEta) {
    file->Close();
    return nullptr;
  }
  // ProjectionY(name, Xmin, Xmax, Zmin, Zmax) — integrate over pT (X) and phi (Z)
  TH1 *JetEta = H3JetEta->ProjectionY(Form("JetEta_%s", histName), 1, H3JetEta->GetNbinsX(), 1, H3JetEta->GetNbinsZ());
  if (!JetEta) {
    file->Close();
    return nullptr;
  }
  JetEta->SetDirectory(0);
  // Rebin by 5: ensures bin edges at multiples of 0.1, aligning with fiducial |eta|<0.9-R
  TH1 *h1JetEta = JetEta->Rebin(5, Form("JetEta_rebinned_%s", histName));
  if (h1JetEta) {
    h1JetEta->SetDirectory(0);
    delete JetEta;
  } else {
    h1JetEta = JetEta;
  }
  legend->AddEntry(h1JetEta, histName);
  // Apply dEta normalization (R-dependent or fixed, from config)
  Double_t deltaEta = GetDeltaEta(R);
  hset(*h1JetEta, GetJetEtaTitleX(), Form("1/#it{N}_{evt} d#it{N}/(d#it{#eta})"), 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
  // R-dependent display range: show only fiducial |eta| < 0.9-R, offset by 0.01 to avoid edge bins
  Double_t etaDisplay = GetJetEtaMax(R) - 0.01;
  if (GetNORMEVENTS()) {
    hoptset(*h1JetEta, Nevts, colorID, -etaDisplay, etaDisplay, 0., 2.0);
  } else {
    hoptset(*h1JetEta, Nevts, colorID, -etaDisplay, etaDisplay, 0., 1.5);
  }
  // Normalize by dEta
  if (deltaEta > 0) {
    h1JetEta->Scale(1.0 / deltaEta);
  }
  h1JetEta->Draw("esame");
  file->Close();
  return h1JetEta;
}

inline TH1 *DrawJetPhi(const char *fileName, const char *histName, const char *Obj, Double_t Nevts,
                TLegend *legend, Color_t colorID, const char *Dir = nullptr, Double_t R = 0.4) {
  auto file = TFile::Open(fileName, "open");
  if (!file || file->IsZombie()) return nullptr;
  TH3D *H3JetPhi = (TH3D *)file->Get(Form("%s/%s", Dir, Obj));
  if (!H3JetPhi) {
    file->Close();
    return nullptr;
  }
  // ProjectionZ(name, Xmin, Xmax, Ymin, Ymax) — integrate over pT (X) and eta (Y)
  int etaBinLo = 1, etaBinHi = H3JetPhi->GetNbinsY();
  double cut = GetCurrentConfigSet().jetEtaCut;
  if (cut >= 0) {
    etaBinLo = H3JetPhi->GetYaxis()->FindBin(-cut + 1e-6);
    etaBinHi = H3JetPhi->GetYaxis()->FindBin( cut - 1e-6);
  }
  TH1D *JetPhi = H3JetPhi->ProjectionZ(Form("JetPhi_%s", histName), 1, H3JetPhi->GetNbinsX(), etaBinLo, etaBinHi);
  if (!JetPhi) {
    file->Close();
    return nullptr;
  }
  JetPhi->SetDirectory(0);
  // Rebin by 4: finer bins (0.2 width) to reduce edge effects at phi~2pi boundary
  TH1 *h1JetPhi = JetPhi->Rebin(4, Form("JetPhi_rebinned_%s", histName));
  if (h1JetPhi) {
    h1JetPhi->SetDirectory(0);
    delete JetPhi;
    JetPhi = (TH1D*)h1JetPhi;
  }
  legend->AddEntry(JetPhi, histName, "p");
  // Apply dEta normalization (R-dependent or fixed, from config)
  Double_t deltaEta = GetDeltaEta(R);
  hset(*JetPhi, GetJetPhiTitleX(), Form("1/#it{N}_{evt} d#it{N}/(d#it{#varphi}d#it{#eta})"), 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
  // Display range clipped to [0.01, 6.19] to avoid partially populated edge bins
  // (phiAxis = {160, -1, 7}, physics range [0, 2pi~6.28], last full bin edge ~6.2)
  if (GetNORMEVENTS()) {
    hoptset(*JetPhi, Nevts, colorID, 0.01, 6.19, 0., 0.3);
  } else {
    hoptset(*JetPhi, Nevts, colorID, 0.01, 6.19, 0., 0.3);
  }
  // Normalize by dEta
  if (deltaEta > 0) {
    JetPhi->Scale(1.0 / deltaEta);
  }
  JetPhi->SetMarkerStyle(20);
  JetPhi->SetMarkerColor(colorID);
  JetPhi->SetMarkerSize(0.8);
  JetPhi->Draw("esame");
  file->Close();
  return JetPhi;
}

inline TH1 *DrawJetPtMCP(const char *fileName, const char *histName, const char *Obj, Double_t Nevts,
                  TLegend *legend, Color_t colorID, Int_t i = 0, const char *Dir = nullptr, const char *legendLabel = nullptr) {
  TFile *file = TFile::Open(fileName, "open");
  if (!file || file->IsZombie()) {
    if (file) file->Close();
    return nullptr;
  }
  
  TH1 *JetPtMCP = nullptr;
  TString histPath = Form("%s/%s", Dir, Obj);
  
  TH1 *tempHist = (TH1 *)file->Get(histPath.Data());
  if (tempHist) {
    JetPtMCP = (TH1 *)tempHist->Clone(Form("JetPtMCP_%s", histName));
    JetPtMCP->SetDirectory(0);
  } else {
    TH3 *h3JetPtMCP = (TH3 *)file->Get(histPath.Data());
    if (h3JetPtMCP) {
      JetPtMCP = h3JetPtMCP->ProjectionX(Form("JetPtMCP_%s", histName), 1, h3JetPtMCP->GetNbinsY(), 1, h3JetPtMCP->GetNbinsZ());
      if (JetPtMCP) {
        JetPtMCP->SetDirectory(0);
      }
    } else {
      file->Close();
      return nullptr;
    }
  }
  
  file->Close();
  
  if (!JetPtMCP) return nullptr;
  
  if (GetREBINON()) {
    PtBinning truthBin = DetectTruthBinning(JetPtMCP);
    TH1 *rebinnedHist = RebinToTarget(JetPtMCP, truthBin, Form("JetPtMCP_%s_rebinned", histName));
    if (rebinnedHist) {
      delete JetPtMCP;
      JetPtMCP = rebinnedHist;
    }
  }
  
  // Use provided legend label or fall back to histName
  const char* labelToUse = (legendLabel && strlen(legendLabel) > 0) ? legendLabel : histName;
  legend->AddEntry(JetPtMCP, labelToUse);
  JetPtMCP->SetTitle("");  // Remove histogram title to avoid "partvjet pT" text
  hset(*JetPtMCP, GetJetPtGenTitleX(), "1/#it{N}_{evt} d#it{N}/d#it{p}_{T}", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510); 
  hoptset(*JetPtMCP, Nevts, colorID, 0, GetPlotPtMax(), 1e-8, 1e-0, 1.2 - 0.2 * i);
  // Note: Don't draw here - let the caller decide when to draw
  // JetPtMCP->Draw("esame");
  
  return JetPtMCP;
}

inline TH2 *DrawJetArea(const char *fileName, const char *histName, const char *Obj, Double_t Nevts,
                 TLegend *legend, Color_t colorID, const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) return nullptr;
  
  TH2 *JetArea = (TH2 *)file->Get(Form("%s/%s", Dir, Obj));
  if (!JetArea) {
    file->Close();
    return nullptr;
  }
  
  // Clone and detach from file before closing
  JetArea = (TH2*)JetArea->Clone(Form("JetArea_%s", histName));
  JetArea->SetDirectory(0);
  
  file->Close();
  
  legend->AddEntry(JetArea, histName, "l");
  hset(*JetArea, GetJetPtTitleX(), "#it{A}_{jet}", 1.2, 1.0, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
  JetArea->Scale(1. / JetArea->Integral(), "width");
  JetArea->GetXaxis()->SetRangeUser(0., GetPlotPtMax());
  JetArea->GetYaxis()->SetRangeUser(0., 2.);
  JetArea->GetZaxis()->SetRangeUser(1e-7, 1e0);
  JetArea->SetStats(0);  // Disable stat box
  JetArea->Draw("colz");
  
  return JetArea;
}

// DrawJetMatching - R-dependent unfolding with result caching
// Caches unfolded histograms by (MCfile, Datafile, Dir, RefDir) to avoid
// redundant unfolding when called from multiple functions for the same R/MC pair.
// Call ClearUnfoldCache() between config sets to reset.
static std::map<TString, TH1F*> gUnfoldCache;
inline void ClearUnfoldCache() {
  for (auto &p : gUnfoldCache) delete p.second;
  gUnfoldCache.clear();
  std::cerr << "[Cache] Unfolding cache cleared" << std::endl;
}

inline TH1 *DrawJetMatching(const char *MCfileName, const char *DfileName, const char *histName,
                     Double_t NevtsData, Double_t NevtsDataUnTrig, Double_t NevtsMCD, Double_t NevtsMCP,
                     TPad *jrep, TPad *ratiojrep, TLegend *JRElegend, TPad *consistencyp,
                     TPad *ratioconsistp, TLegend *CONSISTlegend, TPad *unfoldp, TPad *ratunfoldp,
                     TLegend *unfoldlegend, Color_t colorID, TPad *jrpp = nullptr,
                     TPad *ratiojrpp = nullptr, TLegend *JRPlegend = nullptr, TFile *savefile = nullptr,
                     const char *RefDir = nullptr, const char *Dir = nullptr, const char *mccollcounterfile = nullptr,
                     const char *RecoJetObj = nullptr,
                              const char *TrueJetObj = nullptr,
                              const char *Resp2DObj = nullptr,
                              Int_t svdKOverride = -1) {
  // Use mode-dependent defaults if not provided
  if (!RecoJetObj) {
    RecoJetObj = (GetCurrentMode() == kUE) ? "h_jet_pt_rhoareasubtracted" : "h_jet_pt";
  }
  if (!TrueJetObj) {
    TrueJetObj = (GetCurrentMode() == kUE) ? "h_jet_pt_part_rhoareasubtracted" : "h_jet_pt_part";
  }
  if (!Resp2DObj) {
    Resp2DObj = (GetCurrentMode() == kUE) ? "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted_mcdetaconstraint" : "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint";
  }

  // --- Cache lookup ---
  TString cacheKey = Form("%s|%s|%s|%s|%s", MCfileName, DfileName, Dir ? Dir : "", RefDir ? RefDir : "",
                          GetUnfoldMethod() == kUnfoldBayes ? "Bayes" : "SVD");
  {
    auto cacheIt = gUnfoldCache.find(cacheKey);
    if (cacheIt != gUnfoldCache.end() && cacheIt->second) {
      std::cerr << "[Cache] Returning cached unfolded result for " << histName << std::endl;
      TH1F *hDcorrected = (TH1F *)cacheIt->second->Clone(Form("hDcorrected_%s", histName));
      hDcorrected->SetDirectory(0);
      hDcorrected->SetStats(0);
      return hDcorrected;
    }
  }

  // --- Full unfolding (runs once per unique R/MC/mode) ---
  static int n = 0;

  // Open data file
  auto Dfile = TFile::Open(DfileName, "READ");
  if (!Dfile || Dfile->IsZombie()) {
    std::cerr << "[Error] Cannot open data file: " << DfileName << std::endl;
    return nullptr;
  }

  TH1 *DJetpt = GetJetPt1D(Dfile, RefDir, Form("DJetpt_%d", n));
  if (!DJetpt) {
    std::cerr << "[Error] Cannot find reco jet pt histogram in data file (jetEtaCut="
              << GetCurrentConfigSet().jetEtaCut << ")" << std::endl;
    Dfile->Close();
    return nullptr;
  }

  // Open MC file
  auto MCfile = TFile::Open(MCfileName, "READ");
  if (!MCfile || MCfile->IsZombie()) {
    std::cerr << "[Error] Cannot open MC file: " << MCfileName << std::endl;
    Dfile->Close();
    return nullptr;
  }

  // MC truth: only 1D available (no eta axis), uses production-level eta as-is
  TH1 *JetMCPPt = (TH1 *)MCfile->Get(Form("%s/%s", Dir, TrueJetObj));
  if (!JetMCPPt) {
    std::cerr << "[Error] Cannot find true jet pt histogram (\""<< TrueJetObj <<"\") in MC file" << std::endl;
    MCfile->Close();
    Dfile->Close();
    return nullptr;
  }

  // MC reco: use 3D projection with jetEtaCut when set (consistent with data)
  TH1 *JetMCDPt = GetJetPt1D(MCfile, Dir, Form("JetMCDPt_%d", n));
  if (!JetMCDPt) {
    std::cerr << "[Error] Cannot find reco jet pt histogram in MC file" << std::endl;
    MCfile->Close();
    Dfile->Close();
    return nullptr;
  }

  // Load 2D correlation (response) histogram
  TH2 *HCorrelate2D = (TH2 *)MCfile->Get(Form("%s/%s", Dir, Resp2DObj));
  if (!HCorrelate2D) {
    std::cerr << "[Error] Cannot find 2D correlation (response) histogram (\""<< Resp2DObj <<"\")" << std::endl;
    MCfile->Close();
    Dfile->Close();
    return nullptr;
  }

  // Detect binning from response matrix (authoritative source for consistent rebinning)
  PtBinning recoBin, truthBin;
  DetectResponseBinning(HCorrelate2D, recoBin, truthBin);

  // Rebin all 1D histograms with detected binning (safe: handles narrower data range)
  if (GetREBINON()) {
    DJetpt   = RebinToTarget(DJetpt,   recoBin,  Form("DJetpt_%s", histName));
    JetMCPPt = RebinToTarget(JetMCPPt, truthBin, Form("JetMCPPt_%s", histName));
    JetMCDPt = RebinToTarget(JetMCDPt, recoBin,  Form("JetMCDPt_%s", histName));
  }

  // Rebin 2D response
  TH2F *h2HCorrelate = nullptr;
  if (GetREBINON()) {
    h2HCorrelate = new TH2F(Form("hcorrelate_%s", histName), Form("correlate_%s", histName),
                            recoBin.nBins, recoBin.bins, truthBin.nBins, truthBin.bins);
    for (Int_t ix = 1; ix <= HCorrelate2D->GetNbinsX(); ++ix) {
      for (Int_t iy = 1; iy <= HCorrelate2D->GetNbinsY(); ++iy) {
        Double_t content = HCorrelate2D->GetBinContent(ix, iy);
        Double_t error = HCorrelate2D->GetBinError(ix, iy);
        Int_t xbin = h2HCorrelate->GetXaxis()->FindBin(HCorrelate2D->GetXaxis()->GetBinCenter(ix));
        Int_t ybin = h2HCorrelate->GetYaxis()->FindBin(HCorrelate2D->GetYaxis()->GetBinCenter(iy));
        Double_t currentContent = h2HCorrelate->GetBinContent(xbin, ybin);
        Double_t currentError = h2HCorrelate->GetBinError(xbin, ybin);
        h2HCorrelate->SetBinContent(xbin, ybin, currentContent + content);
        h2HCorrelate->SetBinError(xbin, ybin, std::sqrt(currentError * currentError + error * error));
      }
    }
  } else {
    h2HCorrelate = (TH2F *)HCorrelate2D->Clone(Form("hcorrelate_%s", histName));
  }

  // Project to get matched distributions
  TH1 *MCDMatchedpt = (TH1 *)h2HCorrelate->ProjectionX(Form("hMCDMatched_%s", histName), 1, h2HCorrelate->GetNbinsY(), "e");
  TH1 *MCPMatchedpt = (TH1 *)h2HCorrelate->ProjectionY(Form("hMCPMatched_%s", histName), 1, h2HCorrelate->GetNbinsX(), "e");

  // Clone rebinned 2D response (raw matched counts)
  TH2F *Respt = (TH2F *)h2HCorrelate->Clone(Form("hist_%i", ++n));
  if (Respt->GetSumw2N() == 0) Respt->Sumw2();

  // Kernel-style: 3-arg constructor (matched pairs only, no Fill/Miss/Fake)
  RooUnfoldResponse *Response = new RooUnfoldResponse(MCDMatchedpt, MCPMatchedpt, Respt);

  // Purity from ORIGINAL counts
  TH1F *hPurity = (TH1F *)JetMCDPt->Clone(Form("hPurity_%s", histName));
  hPurity->SetDirectory(nullptr);
  hPurity->Reset();
  for (int i = 1; i <= hPurity->GetNbinsX(); ++i) {
    const double denom = JetMCDPt->GetBinContent(i);
    const double numer = MCDMatchedpt->GetBinContent(i);
    double p = (denom > 0) ? (numer / denom) : 0.0;
    hPurity->SetBinContent(i, TMath::Max(0.0, TMath::Min(1.0, p)));
  }

  // Efficiency from ORIGINAL counts
  TH1F *hEfficiency = (TH1F *)JetMCPPt->Clone(Form("hEfficiency_%s", histName));
  hEfficiency->SetDirectory(nullptr);
  hEfficiency->Reset();
  for (int i = 1; i <= hEfficiency->GetNbinsX(); ++i) {
    const double denom = JetMCPPt->GetBinContent(i);
    const double numer = MCPMatchedpt->GetBinContent(i);
    double eff = (denom > 0) ? (numer / denom) : 0.0;
    hEfficiency->SetBinContent(i, TMath::Max(0.0, TMath::Min(1.0, eff)));
  }

  // Purity correction BEFORE unfolding
  TH1F *DJetptMatched = (TH1F *)DJetpt->Clone(Form("DJetptMatched_%s", histName));
  DJetptMatched->SetDirectory(nullptr);
  for (int i = 1; i <= DJetptMatched->GetNbinsX(); ++i) {
    const double x = DJetpt->GetBinContent(i);
    const double ex = DJetpt->GetBinError(i);
    const double p = hPurity->GetBinContent(i);
    DJetptMatched->SetBinContent(i, x * p);
    DJetptMatched->SetBinError(i, ex * p);
  }

  // --- Determine regularization and unfold ---
  TString outputDir = GetOutputDir();
  Int_t savedErrLvl = gErrorIgnoreLevel;

  TH1F *hDcorrectedTemp = nullptr;
  RooUnfold *unfoldPtr = nullptr;
  Int_t optimalReg = -1;
  TString methodLabel;

  if (GetUnfoldMethod() == kUnfoldBayes) {
    // Bayesian: chi2/ndf + bilateral stability
    Int_t optIter = FindOptimalBayesIter(Response, DJetptMatched, 20, 2.0, 0.03,
                                          histName, outputDir.Data());
    if (optIter < 2) optIter = 2;
    optimalReg = optIter;
    methodLabel = Form("Bayes(iter=%d)", optIter);

    std::cerr << "[Info] Performing Bayesian unfolding with iter=" << optIter << std::endl;
    gErrorIgnoreLevel = kFatal;
    auto *unfoldBayes = new RooUnfoldBayes(Response, DJetptMatched, optIter);
    hDcorrectedTemp = (TH1F *)unfoldBayes->Hreco();
    unfoldPtr = unfoldBayes;
  } else {
    // SVD: d-vector + bilateral stability
    Int_t optK = FindOptimalSvdK_Dvector(Response, DJetptMatched, 0.03,
                                          histName, outputDir.Data());
    if (svdKOverride > 0 && svdKOverride != optK) {
      std::cerr << "[Info] MC closure k=" << svdKOverride
                << " vs data d-vector k=" << optK
                << " (using data, standard practice)" << std::endl;
    }
    if (optK < 2) optK = 2;
    if (optK > 25) optK = 25;
    optimalReg = optK;
    methodLabel = Form("SVD(k=%d)", optK);

    std::cerr << "[Info] Performing SVD unfolding with k=" << optK << std::endl;
    gErrorIgnoreLevel = kFatal;
    auto *unfoldSvd = new RooUnfoldSvd(Response, DJetptMatched, optK);
    hDcorrectedTemp = (TH1F *)unfoldSvd->Hreco();
    unfoldPtr = unfoldSvd;
  }

  if (!hDcorrectedTemp) {
    gErrorIgnoreLevel = savedErrLvl;
    std::cerr << "[Error] Unfolding failed (" << methodLabel << ")" << std::endl;
    delete unfoldPtr;
    MCfile->Close(); Dfile->Close(); delete Response;
    return nullptr;
  }

  TH1F *hDcorrected = (TH1F *)hDcorrectedTemp->Clone(Form("hDcorrected_%s", histName));
  hDcorrected->SetDirectory(0);
  hDcorrected->SetStats(0);

  // Statistical errors from covariance matrix
  TMatrixD covMatrix = unfoldPtr->Ereco();
  gErrorIgnoreLevel = savedErrLvl;
  int nBins = hDcorrected->GetNbinsX();
  bool useCovMatrix = (covMatrix.GetNrows() >= nBins && covMatrix.GetNcols() >= nBins);
  for (int i = 1; i <= nBins; ++i) {
    int idx = i - 1;
    if (useCovMatrix && idx < covMatrix.GetNrows()) {
      double covDiag = covMatrix(idx, idx);
      if (covDiag > 0) hDcorrected->SetBinError(i, TMath::Sqrt(covDiag));
    }
  }

  // Refolding test QA plot (before efficiency correction)
  PlotRefoldingTest(Response, DJetptMatched, hDcorrectedTemp, optimalReg, histName, outputDir.Data(), NevtsData);

  // Efficiency correction AFTER unfolding
  for (int i = 1; i <= hDcorrected->GetNbinsX(); ++i) {
    const int effBin = hEfficiency->GetXaxis()->FindBin(hDcorrected->GetXaxis()->GetBinCenter(i));
    Double_t eff = (effBin >= 1 && effBin <= hEfficiency->GetNbinsX()) ? hEfficiency->GetBinContent(effBin) : 0.0;
    if (eff > 0) {
      hDcorrected->SetBinContent(i, hDcorrected->GetBinContent(i) / eff);
      hDcorrected->SetBinError(i, hDcorrected->GetBinError(i) / eff);
    } else {
      hDcorrected->SetBinContent(i, 0.0);
      hDcorrected->SetBinError(i, 0.0);
    }
  }

  delete unfoldPtr;
  std::cerr << "[Info] " << methodLabel << " unfolding completed for " << histName << std::endl;

  // Draw Jet Reconstruction Efficiency if pads provided
  if (jrep && ratiojrep && JRElegend) {
    jrep->cd();
    auto JREp = (TH1 *)JetMCPPt->Clone(Form("hist_%i", ++n));
    JREp->GetXaxis()->SetRangeUser(5, GetPlotPtMax());
    JREp->SetStats(0);
    JRElegend->AddEntry("", histName, "");
    JRElegend->AddEntry(JREp, "Generated jets");
    hset(*JREp, GetJRETitleX(), GetJRETitleY(), 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    hoptset(*JREp, NevtsMCP, kBlack, 10, GetPlotPtMax(), 3e-10, 1e-1);

    auto mcpmatchedpt = (TH1 *)MCPMatchedpt->Clone(Form("hist_%i", ++n));
    mcpmatchedpt->GetXaxis()->SetRangeUser(0, GetPlotPtMax());
    mcpmatchedpt->SetStats(0);
    JRElegend->AddEntry(mcpmatchedpt, "Matched jets");
    hset(*mcpmatchedpt, GetJRETitleX(), GetJRETitleY(), 0.8, 1.4, 0.04, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    hoptset(*mcpmatchedpt, NevtsMCP, kRed, 10, GetPlotPtMax(), 3e-10, 1e-1);
    JREp->GetYaxis()->SetNdivisions(505);
    JREp->Draw("PE");
    mcpmatchedpt->Draw("PEsame");

    ratiojrep->cd();
    TH1F *__jetpt_ratio_frame = new TH1F(Form("frame_ratio_%s", histName), "", 100, GetPlotPtMin(), GetPlotPtMax());
    __jetpt_ratio_frame->SetDirectory(nullptr);
    __jetpt_ratio_frame->SetStats(0);
    __jetpt_ratio_frame->SetMinimum(0.5);
    __jetpt_ratio_frame->SetMaximum(1.01);
    __jetpt_ratio_frame->GetXaxis()->SetTitle(GetJRETitleX());
    __jetpt_ratio_frame->GetYaxis()->SetTitle("#it{#varepsilon}_{reco}^{jet}");
    __jetpt_ratio_frame->SetLineColor(0);
    __jetpt_ratio_frame->SetMarkerSize(0);
    __jetpt_ratio_frame->Draw("AXIS");

    auto jre = (TH1 *)mcpmatchedpt->Clone(Form("hist_%i", ++n));
    jre->Divide(mcpmatchedpt, JREp, 1., 1., "B");
    jre->SetStats(0);
    hset(*jre, GetJRETitleX(), "#it{#varepsilon}_{reco}^{jet}", 1.2, 1.0, 0.07, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 510);
    jre->SetMarkerColor(kRed);
    jre->SetLineColor(kRed);
    jre->SetMarkerSize(.7);
    jre->SetMarkerStyle(22);
    jre->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
    jre->GetYaxis()->SetRangeUser(0.5, 1.01);
    jre->SetFillColorAlpha(kRed, 0.3);
    jre->Draw("pe");
  }

  // Draw Jet Reconstruction Purity if pads provided
  if (jrpp && ratiojrpp && JRPlegend) {
    jrpp->cd();
    auto JRPp = (TH1F *)JetMCDPt->Clone(Form("hist_%i", ++n));
    JRPp->SetStats(0);
    JRPlegend->AddEntry("", histName, "");
    JRPlegend->AddEntry(JRPp, "Detector level jets");
    hset(*JRPp, GetJRPTitleX(), GetJRPTitleY(), 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    hoptset(*JRPp, NevtsMCD, kBlack, 5, GetPlotPtMax(), 5e-9, 1e-3, 1, 1, 1, 20);

    auto mcdmatchedpt = (TH1F *)MCDMatchedpt->Clone(Form("hist_%i", ++n));
    mcdmatchedpt->SetStats(0);
    JRPlegend->AddEntry(mcdmatchedpt, "Matched jets in Detector level");
    hset(*mcdmatchedpt, GetJRPTitleX(), GetJRPTitleY(), 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    hoptset(*mcdmatchedpt, NevtsMCD, kRed, 5, GetPlotPtMax(), 3e-10, 1e-3, 1, 1, 1, 24);
    JRPp->GetYaxis()->SetNdivisions(505);
    JRPp->Draw("pe");
    mcdmatchedpt->Draw("pesame");

    ratiojrpp->cd();
    auto jrp = (TH1F *)MCDMatchedpt->Clone(Form("hist_%i", ++n));
    jrp->Divide(jrp, JetMCDPt, 1., 1., "B");
    jrp->SetStats(0);
    hset(*jrp, GetJRPTitleX(), "Jet purity", 1.2, 1.0, 0.07, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 505);
    jrp->SetMarkerColor(kRed);
    jrp->SetLineColor(kRed);
    jrp->SetMarkerSize(.7);
    jrp->SetMarkerStyle(24);
    jrp->GetXaxis()->SetRangeUser(0., GetPlotPtMax());
    jrp->GetYaxis()->SetRangeUser(0, 1.01);
    jrp->SetFillColorAlpha(kRed, 0.3);
    jrp->Draw("pe");
  }

  // Draw unfolded result if pad provided
  if (unfoldp && unfoldlegend) {
    unfoldp->cd();
    gPad->SetTicks(1, 1);
    unfoldlegend->AddEntry(hDcorrected, histName);
    hset(*hDcorrected, GetJetPtGenTitleX(), "1/#it{N}_{evt} d#it{N}/d#it{p}_{T}", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    hoptset(*hDcorrected, 0, colorID, GetPlotPtMin(), GetPlotPtMax(), 1e-9, 1e-1, 0.7, 1, 2, 24);
    hDcorrected->Draw("pe");
  }

  // Store in cache and clean up
  gUnfoldCache[cacheKey] = (TH1F *)hDcorrected->Clone(Form("hCache_%s", histName));
  gUnfoldCache[cacheKey]->SetDirectory(0);

  MCfile->Close();
  Dfile->Close();
  delete Response;

  return hDcorrected;
}

#endif // DRAWJETSMCRDEPENDENTHELPERS_H
