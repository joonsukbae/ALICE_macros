// =============================================================================
// DrawMcClosureTest.C
// MC Closure Test using SVD and Bayesian unfolding methods
// =============================================================================
#include "Filipad2.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TParameter.h>
#include <TString.h>
#include <TSystem.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

// ---- Tee stream buffer: duplicates output to both original stream and log file ----
class TeeBuf : public std::streambuf {
public:
  TeeBuf(std::streambuf *orig, std::streambuf *log) : fOrig(orig), fLog(log) {}
protected:
  int overflow(int c) override {
    if (c == EOF) return !EOF;
    int r1 = fOrig->sputc(c);
    int r2 = fLog->sputc(c);
    return (r1 == EOF || r2 == EOF) ? EOF : c;
  }
  int sync() override {
    fOrig->pubsync();
    fLog->pubsync();
    return 0;
  }
private:
  std::streambuf *fOrig, *fLog;
};

#include <RooUnfoldBayes.h>
#include <RooUnfoldResponse.h>
#include <RooUnfoldSvd.h>
#include <TSVDUnfold_local.h>

// =============================================================================
// Configuration
// =============================================================================

// UE mode configuration
struct UEMode {
  const char *label;       // "NoUE" or "UESub"
  const char *outputDir;   // top-level output directory
  const char *histTrue;    // truth histogram name
  const char *histReco;    // reco histogram name
  const char *histResp2D;  // response matrix histogram name
  bool useUEDir;           // true = use dirNameUE, false = use dirName
};

const UEMode kUEModes[2] = {
  {"NoUE", "plots/MCClosureTest",
   "h_jet_pt_part", "h_jet_pt",
   "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint",
   false},
  {"UESub", "plots/MCClosureTest_UESub",
   "h_jet_pt_part_rhoareasubtracted", "h_jet_pt_rhoareasubtracted",
   "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted_mcdetaconstraint",
   true}
};

// Input MC set structure
struct MCInputSet {
  TString responsePath; // Response file path (contains response matrix + truth)
  TString recoPath;     // Reco file path (contains pseudo-data)
  TString name;         // Set name for output files
  TString dirName;      // Directory name inside ROOT file (non-UE)
  TString runId;        // Run/train number for unique identification
  TString dirNameUE;    // Directory for UE-sub histograms (empty = same as dirName)

  // Returns the directory name for a given UE mode
  TString Dir(const UEMode &mode) const {
    if (mode.useUEDir && dirNameUE.Length() > 0)
      return dirNameUE;
    return dirName;
  }
};

// Analysis binning (reco-level)
Double_t ptbin[21] = {5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20,
                      25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
const Int_t nptBins = 20;

// Generator-level binning (finer, starts from 0)
Double_t ptbinGen[26] = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10,  12,  14,
                         16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
const Int_t nptBinsGen = 25;

// Axis titles
TString JetPtTitleX = "#it{p}_{T, jet} (GeV/#it{c})";
TString JetPtTitleY = "1/N_{evt} dN/d#it{p}_{T}";

// =============================================================================
// Normalization summary (per file)
// =============================================================================
// RESPONSE FILE (input.responsePath):
//   - h_jet_pt_part, h_jet_pt, h2_jet_pt_mcd_jet_pt_mcp_... : response matrix
//     and truth/reco spectra. Used only to BUILD the response; no nEvent from
//     this file is used for closure normalization.
//
// RECO FILE (input.recoPath):
//   - h_jet_pt (MCD): unfolded with response → UNFOLDED spectrum
//   - h_jet_pt_part (MCP): truth spectrum for closure comparison
//   - nEventsUnfolded: from h_collisions_weight bin 3 (JJ) or h_collisions bin
//   3 (MB)
//     → used to normalize UNFOLDED: unfolded / nEventsUnfolded
//   - nEventsTruth: from h_mcColl_counts_weight bin 4 (JJ) or
//   h_mccollisions_weight
//     bin 4 or h_mcColl_counts/h_mccollisions bin 4 (MB)
//     → used to normalize TRUTH: hPseudoTruth / nEventsTruth
//
// Closure comparison: (unfolded/nEventsUnfolded) vs (truth/nEventsTruth)
// =============================================================================

// =============================================================================
// Helper functions
// =============================================================================

template <typename T>
void hset(T &hid, TString xtit = "", TString ytit = "", double titoffx = 0.9,
          double titoffy = 1.2, double titsizex = 0.06, double titsizey = 0.06,
          double labeloffx = 0.01, double labeloffy = 0.001,
          double labelsizex = 0.05, double labelsizey = 0.05, int divx = 510,
          int divy = 510) {
  hid.GetXaxis()->CenterTitle(true);
  hid.GetYaxis()->CenterTitle(true);
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
  hid.GetXaxis()->SetTitle(xtit);
  hid.GetYaxis()->SetTitle(ytit);
}

void optFili(TPad &pid, Int_t gridx, Int_t gridy, Int_t logx, Int_t logy) {
  pid.SetGridx(gridx);
  pid.SetGridy(gridy);
  pid.SetLogx(logx);
  pid.SetLogy(logy);
}

// Get event counts for both reco and truth CONSISTENTLY
// Both must use weighted OR both must use non-weighted
// Returns: nEventsReco, nEventsTruth via reference; returns true if found
bool GetEventCounts(TFile *f, const char *dirName, Double_t &nEventsReco,
                    Double_t &nEventsTruth) {
  nEventsReco = -1;
  nEventsTruth = -1;

  if (!f || f->IsZombie())
    return false;

  TDirectory *dir = dynamic_cast<TDirectory *>(f->Get(dirName));
  if (!dir) {
    std::cerr << "[Warning] Directory not found: " << dirName << std::endl;
    return false;
  }

  // MCD (reco) priority: h_collisions_weighted -> h_collisions
  // MCP (truth) priority: h_mcColl_counts_weight -> h_mcColl_counts ->
  // h_mccollisions_weighted -> h_mccollisions
  TH1 *hCollWeighted = dynamic_cast<TH1 *>(dir->Get("h_collisions_weighted"));
  TH1 *hColl = dynamic_cast<TH1 *>(dir->Get("h_collisions"));
  TH1 *hMcCollWeight = dynamic_cast<TH1 *>(dir->Get("h_mcColl_counts_weight"));
  TH1 *hMcColl = dynamic_cast<TH1 *>(dir->Get("h_mcColl_counts"));
  TH1 *hMcCollisionsWeighted =
      dynamic_cast<TH1 *>(dir->Get("h_mccollisions_weighted"));
  TH1 *hMcCollisions = dynamic_cast<TH1 *>(dir->Get("h_mccollisions"));

  std::cerr << "[Debug] Available histograms:" << std::endl;
  std::cerr << "  h_collisions_weighted: " << (hCollWeighted ? "YES" : "NO")
            << std::endl;
  std::cerr << "  h_collisions: " << (hColl ? "YES" : "NO") << std::endl;
  std::cerr << "  h_mcColl_counts_weight: " << (hMcCollWeight ? "YES" : "NO")
            << std::endl;
  std::cerr << "  h_mcColl_counts: " << (hMcColl ? "YES" : "NO") << std::endl;
  std::cerr << "  h_mccollisions_weighted: "
            << (hMcCollisionsWeighted ? "YES" : "NO") << std::endl;
  std::cerr << "  h_mccollisions: " << (hMcCollisions ? "YES" : "NO")
            << std::endl;

  // Strategy: Use CONSISTENT pair (both weighted or both non-weighted)
  // Weighted pair: h_collisions_weighted + (h_mcColl_counts_weight or
  // h_mccollisions_weighted) Non-weighted pair: h_collisions + (h_mcColl_counts
  // or h_mccollisions)

  bool useWeighted = false;
  if (hCollWeighted && (hMcCollWeight || hMcCollisionsWeighted)) {
    useWeighted = true;
  }

  if (useWeighted) {
    // MCD: h_collisions_weighted bin 3
    nEventsReco = hCollWeighted->GetBinContent(3);
    // MCP: h_mcColl_counts_weight -> h_mccollisions_weighted (priority)
    if (hMcCollWeight) {
      nEventsTruth = hMcCollWeight->GetBinContent(4);
      std::cerr << "[Info] Using WEIGHTED pair: h_collisions_weighted bin 3 = "
                << nEventsReco
                << ", h_mcColl_counts_weight bin 4 = " << nEventsTruth
                << std::endl;
    } else {
      nEventsTruth = hMcCollisionsWeighted->GetBinContent(4);
      std::cerr << "[Info] Using WEIGHTED pair: h_collisions_weighted bin 3 = "
                << nEventsReco
                << ", h_mccollisions_weighted bin 4 = " << nEventsTruth
                << std::endl;
    }
  } else {
    // MCD: h_collisions bin 3
    if (hColl)
      nEventsReco = hColl->GetBinContent(3);
    // MCP: h_mcColl_counts -> h_mccollisions (priority)
    if (hMcColl) {
      nEventsTruth = hMcColl->GetBinContent(4);
      std::cerr << "[Info] Using NON-WEIGHTED pair: h_collisions bin 3 = "
                << nEventsReco << ", h_mcColl_counts bin 4 = " << nEventsTruth
                << std::endl;
    } else if (hMcCollisions) {
      nEventsTruth = hMcCollisions->GetBinContent(4);
      std::cerr << "[Info] Using NON-WEIGHTED pair: h_collisions bin 3 = "
                << nEventsReco << ", h_mccollisions bin 4 = " << nEventsTruth
                << std::endl;
    }
  }

  // Verify we got both
  if (nEventsReco <= 0 || nEventsTruth <= 0) {
    std::cerr << "[Error] Could not find consistent event count pair!"
              << std::endl;
    return false;
  }

  std::cerr << "[Info] Event count ratio (reco/truth): "
            << nEventsReco / nEventsTruth << std::endl;
  return true;
}

// Legacy function for compatibility (deprecated)
Double_t GetEventCountReco(TFile *f, const char *dirName) {
  Double_t nReco = -1, nTruth = -1;
  GetEventCounts(f, dirName, nReco, nTruth);
  return nReco;
}

// Legacy function for compatibility (deprecated)
Double_t GetEventCountTruth(TFile *f, const char *dirName) {
  if (!f || f->IsZombie())
    return -1;

  TDirectory *dir = dynamic_cast<TDirectory *>(f->Get(dirName));
  if (!dir) {
    std::cerr << "[Warning] Directory not found: " << dirName << std::endl;
    return -1;
  }

  // MCP priority: h_mcColl_counts_weight -> h_mcColl_counts ->
  // h_mccollisions_weighted -> h_mccollisions
  TH1 *h = dynamic_cast<TH1 *>(dir->Get("h_mcColl_counts_weight"));
  if (h) {
    Double_t n = h->GetBinContent(4);
    std::cerr << "[Info] GetEventCountTruth: h_mcColl_counts_weight bin 4 = "
              << n << std::endl;
    return n;
  }
  h = dynamic_cast<TH1 *>(dir->Get("h_mcColl_counts"));
  if (h) {
    Double_t n = h->GetBinContent(4);
    std::cerr << "[Info] GetEventCountTruth: h_mcColl_counts bin 4 = " << n
              << std::endl;
    return n;
  }
  h = dynamic_cast<TH1 *>(dir->Get("h_mccollisions_weighted"));
  if (h) {
    Double_t n = h->GetBinContent(4);
    std::cerr << "[Info] GetEventCountTruth: h_mccollisions_weighted bin 4 = "
              << n << std::endl;
    return n;
  }
  h = dynamic_cast<TH1 *>(dir->Get("h_mccollisions"));
  if (h) {
    Double_t n = h->GetBinContent(4);
    std::cerr << "[Info] GetEventCountTruth: h_mccollisions bin 4 = " << n
              << std::endl;
    return n;
  }
  return -1;
}

// =============================================================================
// Histogram loading
// =============================================================================

// Load histograms from MC input set
// nEventsUnfolded: for normalizing unfolded result (from h_collisions in reco
// file) nEventsTruth: for normalizing truth (from h_mcColl_counts in reco file)
bool LoadMCInputSet(const MCInputSet &input, const UEMode &mode,
                    TH1D *&hTruth, TH1D *&hReco,
                    TH2D *&hResp2D, TH1D *&hPseudoData, TH1D *&hPseudoTruth,
                    Double_t &nEventsUnfolded, Double_t &nEventsTruth) {

  TString activeDir = input.Dir(mode);
  std::cerr << "[Info] Loading MC set: " << input.name
            << " [" << mode.label << "]" << std::endl;
  std::cerr << "       Response: " << input.responsePath << std::endl;
  std::cerr << "       Reco:     " << input.recoPath << std::endl;
  std::cerr << "       Dir:      " << activeDir << std::endl;

  // ========================================
  // Load from response file
  // ========================================
  TFile *fResp = TFile::Open(input.responsePath, "READ");
  if (!fResp || fResp->IsZombie()) {
    std::cerr << "[Error] Cannot open response file: " << input.responsePath
              << std::endl;
    return false;
  }

  TDirectory *dirResp = dynamic_cast<TDirectory *>(fResp->Get(activeDir));
  if (!dirResp) {
    std::cerr << "[Skip] Directory not found in response file: "
              << activeDir << std::endl;
    return false;
  }

  // Get raw histograms
  TH1 *hTruthRaw = dynamic_cast<TH1 *>(dirResp->Get(mode.histTrue));
  TH2 *hResp2DRaw = dynamic_cast<TH2 *>(dirResp->Get(mode.histResp2D));
  TH1 *hRecoRaw = dynamic_cast<TH1 *>(dirResp->Get(mode.histReco));

  if (!hTruthRaw || !hResp2DRaw || !hRecoRaw) {
    std::cerr << "[Skip] Missing histograms in response file for " << mode.label << ":" << std::endl;
    std::cerr << "  " << mode.histTrue << ": " << (hTruthRaw ? "found" : "NOT FOUND") << std::endl;
    std::cerr << "  " << mode.histResp2D << ": " << (hResp2DRaw ? "found" : "NOT FOUND") << std::endl;
    std::cerr << "  " << mode.histReco << ": " << (hRecoRaw ? "found" : "NOT FOUND") << std::endl;
    return false;
  }

  // Response file nEvents not needed for closure test (using reco file for
  // normalization)

  // Debug: print histogram types
  std::cerr << "[Debug] hTruthRaw type: " << hTruthRaw->ClassName()
            << ", nBins: " << hTruthRaw->GetNbinsX() << std::endl;
  std::cerr << "[Debug] hRecoRaw type: " << hRecoRaw->ClassName()
            << ", nBins: " << hRecoRaw->GetNbinsX() << std::endl;
  std::cerr << "[Debug] hResp2DRaw type: " << hResp2DRaw->ClassName()
            << ", nBinsX: " << hResp2DRaw->GetNbinsX()
            << ", nBinsY: " << hResp2DRaw->GetNbinsY() << std::endl;

  // Rebin to analysis binning using ROOT's Rebin
  // TH1F::Rebin returns TH1F*, so we must copy to TH1D (dynamic_cast<TH1D*>
  // would be nullptr)
  TH1 *hTruthRebinned = hTruthRaw->Rebin(
      nptBinsGen, Form("hTruth_%s_%s", input.name.Data(), mode.label), ptbinGen);
  hTruth = dynamic_cast<TH1D *>(hTruthRebinned);
  if (!hTruth) {
    hTruth = new TH1D(Form("hTruth_%s_%s", input.name.Data(), mode.label),
                      hTruthRaw->GetTitle(), nptBinsGen, ptbinGen);
    hTruth->SetDirectory(nullptr);
    for (int i = 1; i <= nptBinsGen; ++i) {
      hTruth->SetBinContent(i, hTruthRebinned->GetBinContent(i));
      hTruth->SetBinError(i, hTruthRebinned->GetBinError(i));
    }
  } else {
    hTruth->SetDirectory(nullptr);
  }
  std::cerr << "[Debug] hTruth rebinned (ROOT Rebin)" << std::endl;

  TH1 *hRecoRebinned =
      hRecoRaw->Rebin(nptBins, Form("hReco_%s_%s", input.name.Data(), mode.label), ptbin);
  hReco = dynamic_cast<TH1D *>(hRecoRebinned);
  if (!hReco) {
    hReco = new TH1D(Form("hReco_%s_%s", input.name.Data(), mode.label), hRecoRaw->GetTitle(),
                     nptBins, ptbin);
    hReco->SetDirectory(nullptr);
    for (int i = 1; i <= nptBins; ++i) {
      hReco->SetBinContent(i, hRecoRebinned->GetBinContent(i));
      hReco->SetBinError(i, hRecoRebinned->GetBinError(i));
    }
  } else {
    hReco->SetDirectory(nullptr);
  }
  std::cerr << "[Debug] hReco rebinned (ROOT Rebin)" << std::endl;

  // Build rebinned 2D matrix: X=reco (mcd), Y=true (mcp)
  // Use bin center method as in DrawJetsMCfJetQA.h (line 337-352) and backup
  // code
  hResp2D = new TH2D(Form("hResp2D_%s_%s", input.name.Data(), mode.label),
                     "Response;#it{p}_{T}^{reco};#it{p}_{T}^{true}", nptBins,
                     ptbin, nptBinsGen, ptbinGen);
  hResp2D->SetDirectory(nullptr);
  if (hResp2D->GetSumw2N() == 0)
    hResp2D->Sumw2();

  for (int ix = 1; ix <= hResp2DRaw->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hResp2DRaw->GetNbinsY(); ++iy) {
      double content = hResp2DRaw->GetBinContent(ix, iy);
      double error = hResp2DRaw->GetBinError(ix, iy);
      if (content == 0. && error == 0.)
        continue;

      // Use bin center to find target bin (same as DrawJetsMCfJetQA.h)
      double xCenter = hResp2DRaw->GetXaxis()->GetBinCenter(ix);
      double yCenter = hResp2DRaw->GetYaxis()->GetBinCenter(iy);
      int xbin = hResp2D->GetXaxis()->FindBin(xCenter);
      int ybin = hResp2D->GetYaxis()->FindBin(yCenter);

      // Skip overflow/underflow
      if (xbin < 1 || xbin > nptBins || ybin < 1 || ybin > nptBinsGen)
        continue;

      double currentContent = hResp2D->GetBinContent(xbin, ybin);
      double currentError = hResp2D->GetBinError(xbin, ybin);
      hResp2D->SetBinContent(xbin, ybin, currentContent + content);
      hResp2D->SetBinError(
          xbin, ybin, std::sqrt(currentError * currentError + error * error));
    }
  }
  std::cerr << "[Debug] hResp2D rebinned (bin center method)" << std::endl;

  // fResp->Close();  // Keep open for ROOT browser QA

  // ========================================
  // Load from reco file (pseudo-data)
  // ========================================
  TFile *fReco = TFile::Open(input.recoPath, "READ");
  if (!fReco || fReco->IsZombie()) {
    std::cerr << "[Error] Cannot open reco file: " << input.recoPath
              << std::endl;
    return false;
  }

  TDirectory *dirReco = dynamic_cast<TDirectory *>(fReco->Get(activeDir));
  if (!dirReco) {
    std::cerr << "[Skip] Directory not found in reco file: " << activeDir
              << std::endl;
    return false;
  }

  TH1 *hPseudoDataRaw = dynamic_cast<TH1 *>(dirReco->Get(mode.histReco));
  TH1 *hPseudoTruthRaw = dynamic_cast<TH1 *>(dirReco->Get(mode.histTrue));

  if (!hPseudoDataRaw) {
    std::cerr << "[Error] Cannot find pseudo-data histogram" << std::endl;
    // fReco->Close();  // Keep open for ROOT browser QA
    return false;
  }

  // Get nEvents from reco file for normalization
  // IMPORTANT: Both must use consistent pair (weighted or non-weighted)
  if (!GetEventCounts(fReco, activeDir, nEventsUnfolded, nEventsTruth)) {
    std::cerr << "[Error] Failed to get event counts" << std::endl;
    return false;
  }

  // Use ROOT's Rebin; copy to TH1D if Rebin returns TH1F
  TH1 *hPseudoDataRebinned = hPseudoDataRaw->Rebin(
      nptBins, Form("hPseudoData_%s_%s", input.name.Data(), mode.label), ptbin);
  hPseudoData = dynamic_cast<TH1D *>(hPseudoDataRebinned);
  if (!hPseudoData) {
    hPseudoData = new TH1D(Form("hPseudoData_%s_%s", input.name.Data(), mode.label),
                           hPseudoDataRaw->GetTitle(), nptBins, ptbin);
    hPseudoData->SetDirectory(nullptr);
    for (int i = 1; i <= nptBins; ++i) {
      hPseudoData->SetBinContent(i, hPseudoDataRebinned->GetBinContent(i));
      hPseudoData->SetBinError(i, hPseudoDataRebinned->GetBinError(i));
    }
  } else {
    hPseudoData->SetDirectory(nullptr);
  }

  if (hPseudoTruthRaw) {
    TH1 *hPseudoTruthRebinned = hPseudoTruthRaw->Rebin(
        nptBinsGen, Form("hPseudoTruth_%s_%s", input.name.Data(), mode.label), ptbinGen);
    hPseudoTruth = dynamic_cast<TH1D *>(hPseudoTruthRebinned);
    if (!hPseudoTruth) {
      hPseudoTruth =
          new TH1D(Form("hPseudoTruth_%s_%s", input.name.Data(), mode.label),
                   hPseudoTruthRaw->GetTitle(), nptBinsGen, ptbinGen);
      hPseudoTruth->SetDirectory(nullptr);
      for (int i = 1; i <= nptBinsGen; ++i) {
        hPseudoTruth->SetBinContent(i, hPseudoTruthRebinned->GetBinContent(i));
        hPseudoTruth->SetBinError(i, hPseudoTruthRebinned->GetBinError(i));
      }
    } else {
      hPseudoTruth->SetDirectory(nullptr);
    }
  } else {
    hPseudoTruth = nullptr;
  }

  // fReco->Close();  // Keep open for ROOT browser QA

  std::cerr << "[Info] Histograms loaded successfully" << std::endl;
  return true;
}

// =============================================================================
// Build RooUnfoldResponse
// =============================================================================

// Build RooUnfoldResponse for kernel-style unfolding (matched pairs only)
// - Purity (fake) and efficiency (miss) are handled externally.
RooUnfoldResponse *BuildResponse(TH1D *hTrue, TH1D *hReco, TH2D *hResp2D,
                                 TH1D *&hPurity, TH1D *&hEfficiency) {
  hPurity = nullptr;
  hEfficiency = nullptr;
  if (!hTrue || !hReco || !hResp2D) {
    std::cerr << "[Error] BuildResponse: null input" << std::endl;
    return nullptr;
  }

  // Clone the already rebinned 2D matrix (same as backup code)
  TH2D *Respt = dynamic_cast<TH2D *>(hResp2D->Clone("Respt"));
  Respt->SetDirectory(nullptr);

  // Get matched projections from 2D matrix
  TH1D *hRecoMatched = dynamic_cast<TH1D *>(
      hResp2D->ProjectionX("hRecoMatched", 1, hResp2D->GetNbinsY(), "e"));
  TH1D *hTrueMatched = dynamic_cast<TH1D *>(
      hResp2D->ProjectionY("hTrueMatched", 1, hResp2D->GetNbinsX(), "e"));
  hRecoMatched->SetDirectory(nullptr);
  hTrueMatched->SetDirectory(nullptr);

  // Kernel-style correction factors
  // Purity(pT_reco) = matchedReco / reco
  // Efficiency(pT_true) = matchedTruth / truth
  hPurity = dynamic_cast<TH1D *>(hReco->Clone("hPurity"));
  if (hPurity) {
    hPurity->SetDirectory(nullptr);
    hPurity->Reset();
    for (int i = 1; i <= hPurity->GetNbinsX(); ++i) {
      const double denom = hReco->GetBinContent(i);
      const double numer = hRecoMatched->GetBinContent(i);
      double p = (denom > 0) ? (numer / denom) : 0.0;
      if (p < 0)
        p = 0.0;
      if (p > 1)
        p = 1.0;
      hPurity->SetBinContent(i, p);
    }
  }

  hEfficiency = dynamic_cast<TH1D *>(hTrue->Clone("hEfficiency"));
  if (hEfficiency) {
    hEfficiency->SetDirectory(nullptr);
    hEfficiency->Reset();
    for (int i = 1; i <= hEfficiency->GetNbinsX(); ++i) {
      const double denom = hTrue->GetBinContent(i);
      const double numer = hTrueMatched->GetBinContent(i);
      double eff = (denom > 0) ? (numer / denom) : 0.0;
      if (eff < 0)
        eff = 0.0;
      if (eff > 1)
        eff = 1.0;
      hEfficiency->SetBinContent(i, eff);
    }
  }

  // ========================================
  // Step 3: Calculate fake and miss
  // fake = reco - matched_reco
  // miss = true - matched_true
  // ========================================
  TH1D *fake = dynamic_cast<TH1D *>(hReco->Clone("hFake_tmp"));
  fake->SetDirectory(nullptr);
  fake->Add(hRecoMatched, -1);

  TH1D *miss = dynamic_cast<TH1D *>(hTrue->Clone("hMiss_tmp"));
  miss->SetDirectory(nullptr);
  miss->Add(hTrueMatched, -1);

  // Debug output
  std::cerr << "[Info] Response matrix stats:" << std::endl;
  std::cerr << "  Reco total: " << hReco->Integral()
            << ", Matched: " << hRecoMatched->Integral() << " (ratio: "
            << (hReco->Integral() > 0
                    ? hRecoMatched->Integral() / hReco->Integral()
                    : 0)
            << ")" << std::endl;
  std::cerr << "  True total: " << hTrue->Integral()
            << ", Matched: " << hTrueMatched->Integral() << " (ratio: "
            << (hTrue->Integral() > 0
                    ? hTrueMatched->Integral() / hTrue->Integral()
                    : 0)
            << ")" << std::endl;
  std::cerr << "  Fake: " << fake->Integral() << ", Miss: " << miss->Integral()
            << std::endl;

  // ========================================
  // Step 4: Create RooUnfoldResponse (matched pairs only)
  // ========================================
  // 3-arg constructor: directly sets _res=Respt, _mes=hRecoMatched, _tru=hTrueMatched
  // No Fill() loop — avoids doubling _mes/_tru that occurs with 2-arg + Fill pattern
  // Preserves original bin errors (Sumw2) from the 2D matrix
  RooUnfoldResponse *resp = new RooUnfoldResponse(hRecoMatched, hTrueMatched, Respt);
  std::cerr << "[Info] RooUnfoldResponse built with 3-arg constructor (matched pairs, raw counts)" << std::endl;

  // NOTE: Do not add Miss/Fake into RooUnfoldResponse (handled externally).

  return resp;
}

// =============================================================================
// Calculate ratio by pT (same binning)
// =============================================================================

// Copy truth-space histogram onto the exact axis of response->Htruth()
// ApplyToTruth() uses BIN INDEX: refolded[i] = sum_j R(i,j)*truth(j). Truth bin
// j must match response Y-axis (truth) bin j. Copy BY INDEX when nbins match.
TH1D *CopyOntoTruthAxis(TH1D *hSource, const TH1 *hTruthTemplate) {
  if (!hSource || !hTruthTemplate)
    return nullptr;
  const TAxis *ax = hTruthTemplate->GetXaxis();
  int n = ax->GetNbins();
  TH1D *hOut = nullptr;
  if (ax->GetXbins() && ax->GetXbins()->GetArray()) {
    hOut = new TH1D(Form("%s_onTruthAxis", hSource->GetName()),
                    hSource->GetTitle(), n, ax->GetXbins()->GetArray());
  } else {
    hOut = new TH1D(Form("%s_onTruthAxis", hSource->GetName()),
                    hSource->GetTitle(), n, ax->GetXmin(), ax->GetXmax());
  }
  hOut->SetDirectory(nullptr);
  int nSrc = hSource->GetNbinsX();
  if (nSrc == n) {
    for (int j = 1; j <= n; ++j) {
      hOut->SetBinContent(j, hSource->GetBinContent(j));
      hOut->SetBinError(j, hSource->GetBinError(j));
    }
  } else {
    for (int j = 1; j <= n; ++j) {
      double xCenter = ax->GetBinCenter(j);
      int binSrc = hSource->GetXaxis()->FindBin(xCenter);
      if (binSrc >= 1 && binSrc <= nSrc) {
        hOut->SetBinContent(j, hSource->GetBinContent(binSrc));
        hOut->SetBinError(j, hSource->GetBinError(binSrc));
      }
    }
  }
  return hOut;
}

// Manual refolding using the migration kernel:
// refolded_reco[i] = sum_j K(i|j) * truth_detected[j]
// where K(i|j) = N(i,j) / Ntruth_train(j).
// NOTE: Hresponse() stores COUNTS N(i,j). We must normalise by
// response->Htruth() to get K.
TH1D *RefoldManual(RooUnfoldResponse *response, TH1D *hTruthAligned) {
  if (!response || !hTruthAligned)
    return nullptr;
  TH2D *hR = dynamic_cast<TH2D *>(response->Hresponse());
  TH1 *hMeasured = response->Hmeasured();
  const TH1 *hTruthTrain = response->Htruth();
  if (!hR || !hMeasured || !hTruthTrain)
    return nullptr;
  int nReco = hR->GetNbinsX();
  int nTruth = hR->GetNbinsY();
  if (hTruthAligned->GetNbinsX() != nTruth) {
    std::cerr << "[Error] RefoldManual: truth bins "
              << hTruthAligned->GetNbinsX() << " != response Y bins " << nTruth
              << std::endl;
    return nullptr;
  }
  const TAxis *axMeas = hMeasured->GetXaxis();
  TH1D *hRefolded = nullptr;
  if (axMeas->GetXbins() && axMeas->GetXbins()->GetArray()) {
    hRefolded = new TH1D("hRefolded_manual", "Refolded (manual)", nReco,
                         axMeas->GetXbins()->GetArray());
  } else {
    hRefolded = new TH1D("hRefolded_manual", "Refolded (manual)", nReco,
                         axMeas->GetXmin(), axMeas->GetXmax());
  }
  hRefolded->SetDirectory(nullptr);
  for (int i = 1; i <= nReco; ++i) {
    double sum = 0.0;
    for (int j = 1; j <= nTruth; ++j) {
      const double nTrain = hTruthTrain->GetBinContent(j);
      if (nTrain <= 0) {
        continue;
      }
      const double kij = hR->GetBinContent(i, j) / nTrain;
      sum += kij * hTruthAligned->GetBinContent(j);
    }
    hRefolded->SetBinContent(i, sum);
    hRefolded->SetBinError(i,
                           0.); // optional: propagate from matrix/truth errors
  }
  return hRefolded;
}

// Copy reco-space histogram onto the exact axis of response->Hmeasured() (RM
// X-axis) Refolded output may have different axis; for ratio with measured use
// SAME reco binning.
TH1D *CopyOntoMeasuredAxis(TH1D *hSource, const TH1 *hMeasuredTemplate) {
  if (!hSource || !hMeasuredTemplate)
    return nullptr;
  const TAxis *ax = hMeasuredTemplate->GetXaxis();
  int n = ax->GetNbins();
  TH1D *hOut = nullptr;
  if (ax->GetXbins() && ax->GetXbins()->GetArray()) {
    hOut = new TH1D(Form("%s_onMeasuredAxis", hSource->GetName()),
                    hSource->GetTitle(), n, ax->GetXbins()->GetArray());
  } else {
    hOut = new TH1D(Form("%s_onMeasuredAxis", hSource->GetName()),
                    hSource->GetTitle(), n, ax->GetXmin(), ax->GetXmax());
  }
  hOut->SetDirectory(nullptr);
  int nSrc = hSource->GetNbinsX();
  if (nSrc == n) {
    for (int i = 1; i <= n; ++i) {
      hOut->SetBinContent(i, hSource->GetBinContent(i));
      hOut->SetBinError(i, hSource->GetBinError(i));
    }
  } else {
    for (int i = 1; i <= n; ++i) {
      double xCenter = ax->GetBinCenter(i);
      int binSrc = hSource->GetXaxis()->FindBin(xCenter);
      if (binSrc >= 1 && binSrc <= nSrc) {
        hOut->SetBinContent(i, hSource->GetBinContent(binSrc));
        hOut->SetBinError(i, hSource->GetBinError(binSrc));
      }
    }
  }
  return hOut;
}

TH1D *CalculateRatioByPt(TH1D *hNum, TH1D *hDenom) {
  if (!hNum || !hDenom)
    return nullptr;

  TH1D *hRatio =
      dynamic_cast<TH1D *>(hNum->Clone(Form("%s_ratio", hNum->GetName())));
  hRatio->SetDirectory(nullptr);

  for (int i = 1; i <= hRatio->GetNbinsX(); ++i) {
    double num = hNum->GetBinContent(i);
    double denom = hDenom->GetBinContent(i);
    double numErr = hNum->GetBinError(i);
    double denomErr = hDenom->GetBinError(i);

    if (denom > 0 && num > 0) {
      double ratio = num / denom;
      double relErrNum = numErr / num;
      double relErrDenom = denomErr / denom;
      double ratioErr = ratio * TMath::Sqrt(relErrNum * relErrNum +
                                            relErrDenom * relErrDenom);

      hRatio->SetBinContent(i, ratio);
      hRatio->SetBinError(i, ratioErr);
    } else {
      hRatio->SetBinContent(i, 0);
      hRatio->SetBinError(i, 0);
    }
  }

  return hRatio;
}

// =============================================================================
// Optimization functions
// =============================================================================

// Find optimal Bayesian iteration using convergence criterion.
// Primary: first iteration where max relative bin-by-bin change < threshold (1%).
// This is robust for high-statistics MC where chi2/ndf can remain >> 1
// due to tiny statistical errors.
// Fallback: if no convergence, use chi2/ndf closest to 1.0.
// Find optimal Bayesian iteration using chi2/ndf + bilateral stability.
// Analogous to SVD d-vector + bilateral stability:
//   Step 1: chi2/ndf vs truth sets lower bound (first iter where chi2/ndf < threshold).
//   Step 2: Bilateral stability finds the convergence plateau.
Int_t FindOptimalBayesIter(TH1D *hTruth, const std::vector<TH1D *> &unfolded,
                           Double_t nEventsTruth, Double_t nEventsUnfolded,
                           Double_t chi2Threshold = 2.0,
                           Double_t stabilityThreshold = 0.03) {
  if (unfolded.size() < 2)
    return 4;

  const Int_t maxIter = static_cast<Int_t>(unfolded.size());

  // ---- Step 1: chi2/ndf vs truth ----
  // Find first iteration where chi2/ndf < threshold (analogous to d-vector boundary)
  TH1D *hTruthNorm = dynamic_cast<TH1D *>(hTruth->Clone("hTruthNorm_tmp"));
  hTruthNorm->SetDirectory(nullptr);
  if (nEventsTruth > 0)
    hTruthNorm->Scale(1.0 / nEventsTruth, "width");

  Int_t chi2Iter = -1;
  std::cerr << "[Info] Bayesian chi2/ndf scan (threshold = " << chi2Threshold << "):" << std::endl;

  for (size_t i = 0; i < unfolded.size(); ++i) {
    int iter = static_cast<int>(i) + 1;

    TH1D *hUnfNorm = dynamic_cast<TH1D *>(unfolded[i]->Clone("hUnfNorm_tmp"));
    hUnfNorm->SetDirectory(nullptr);
    if (nEventsUnfolded > 0)
      hUnfNorm->Scale(1.0 / nEventsUnfolded, "width");

    double chi2 = 0.0;
    int ndf = 0;

    for (int bin = 1; bin <= hTruthNorm->GetNbinsX(); ++bin) {
      double ptLow = hTruthNorm->GetXaxis()->GetBinLowEdge(bin);
      double ptHigh = hTruthNorm->GetXaxis()->GetBinUpEdge(bin);
      if (ptLow < 5.0 || ptHigh > 200.0) continue;

      double truth = hTruthNorm->GetBinContent(bin);
      double unfoldedVal = hUnfNorm->GetBinContent(bin);
      double err = hUnfNorm->GetBinError(bin);
      double truthErr = hTruthNorm->GetBinError(bin);
      double totalErr = TMath::Sqrt(err * err + truthErr * truthErr);

      if (truth > 0 && unfoldedVal > 0 && totalErr > 0) {
        double diff = unfoldedVal - truth;
        chi2 += (diff * diff) / (totalErr * totalErr);
        ndf++;
      }
    }

    double reducedChi2 = (ndf > 0) ? chi2 / ndf : 1e10;
    std::cerr << "[Info]   Bayes iter " << iter
              << ": chi2/ndf = " << reducedChi2 << " (ndf=" << ndf << ")"
              << (reducedChi2 < chi2Threshold ? "  [ok]" : "")
              << std::endl;

    if (chi2Iter < 0 && reducedChi2 < chi2Threshold) {
      chi2Iter = iter;
    }

    delete hUnfNorm;
  }

  delete hTruthNorm;

  // If chi2 never drops below threshold, fall back to bilateral stability from iter=2
  Int_t candidateIter = chi2Iter;
  if (chi2Iter < 0) {
    std::cerr << "[Warning] chi2/ndf never dropped below " << chi2Threshold
              << ", starting bilateral stability from iter=2." << std::endl;
    candidateIter = 2;
  }
  std::cerr << "[Info] chi2/ndf boundary: iter=" << (chi2Iter > 0 ? chi2Iter : -1)
            << ", starting bilateral scan from iter=" << candidateIter << std::endl;

  // ---- Step 2: Bilateral stability validation ----
  // Same criterion as SVD: first iter where BOTH |y(iter-1)/y(iter)-1| < ε
  // AND |y(iter)/y(iter+1)-1| < ε across all pT bins in [5, 100] GeV.
  std::cerr << "[Info] Bilateral stability validation (threshold = "
            << stabilityThreshold * 100 << "%):" << std::endl;

  auto maxStepChange = [&](Int_t iterA, Int_t iterB) -> Double_t {
    Int_t idxA = iterA - 1, idxB = iterB - 1;
    if (idxA < 0 || idxA >= maxIter || idxB < 0 || idxB >= maxIter) return 999.;
    Double_t maxRel = 0.0;
    for (Int_t bin = 1; bin <= unfolded[idxA]->GetNbinsX(); ++bin) {
      Double_t ptLow = unfolded[idxA]->GetXaxis()->GetBinLowEdge(bin);
      Double_t ptHigh = unfolded[idxA]->GetXaxis()->GetBinUpEdge(bin);
      if (ptLow < 5.0 || ptHigh > 200.0) continue;
      Double_t a = unfolded[idxA]->GetBinContent(bin);
      Double_t b = unfolded[idxB]->GetBinContent(bin);
      if (a > 0 && b > 0) {
        Double_t rel = TMath::Abs(b / a - 1.0);
        if (rel > maxRel) maxRel = rel;
      }
    }
    return maxRel;
  };

  // Scan for bilateral stability plateau: find ALL stable iter, pick the LAST one.
  // MC closure doesn't have refolding chi2, so we pick the highest stable iter
  // (least regularization bias). Same reasoning as the SVD "last stable k" fix.
  Int_t firstStableIter = -1;
  Int_t lastStableIter = -1;
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
      lastStableIter = it;
    }
  }

  Int_t finalIter = candidateIter;
  Bool_t foundStable = (lastStableIter > 0);
  if (foundStable) {
    finalIter = lastStableIter;
  }

  if (!foundStable) {
    std::cerr << "[Warning] Bayesian bilateral stability never satisfied (threshold "
              << stabilityThreshold * 100 << "%). Using iter=" << finalIter
              << ". Consider increasing maxIter or relaxing threshold." << std::endl;
  } else if (lastStableIter != firstStableIter) {
    std::cerr << "[Info] Bayesian stability: first stable iter=" << firstStableIter
              << ", last stable iter=" << lastStableIter
              << " — using last (least regularization bias)" << std::endl;
  }

  if (finalIter < 2) finalIter = 2;
  if (finalIter > maxIter) finalIter = maxIter;

  std::cerr << "[Info] Optimal Bayes iter (chi2 + bilateral stability): " << finalIter << std::endl;
  return finalIter;
}

// Find optimal SVD k using d-vector method (Höcker & Kartvelishvili, NIM A372, 469, 1996),
// combined with bilateral stability criterion.
//
// Step 1 (d-vector): Find the signal/noise boundary in the SVD basis.
//   k_min = last index where |d_i| > 1.0 (two-consecutive-below criterion).
// Step 2 (bilateral stability): Starting from k_min, find the first k where BOTH
//   |y(k-1)/y(k) - 1| < ε AND |y(k)/y(k+1) - 1| < ε across all pT bins.
//   This ensures k is on the convergence plateau (not at the knee).
//
// This ensures the d-vector result is validated against actual result stability,
// avoiding edge cases where k lands at the convergence knee rather than on the plateau.
//
// Inputs:
//   response      - RooUnfold response matrix
//   hRecoData     - measured (pseudo-)data histogram
//   svdUnfolded   - vector of unfolded histograms for k = 1, 2, ..., N
//                   (index 0 = k=1, etc.), used for stability check
//   stabilityThreshold - max relative bin-by-bin change for bilateral stability (default 3%)
Int_t FindOptimalSvdK_Dvector(RooUnfoldResponse *response, TH1D *hRecoData,
                               const std::vector<TH1D *> &svdUnfolded,
                               Double_t stabilityThreshold = 0.03) {
  if (!response || !hRecoData)
    return 6;

  const Int_t nBins = hRecoData->GetNbinsX();
  const Double_t threshold = 1.0;

  std::cerr << "[Info] Computing d-vector for SVD k optimization..." << std::endl;

  RooUnfoldSvd unfoldFull(response, hRecoData, nBins);
  unfoldFull.Hreco(); // trigger unfolding to populate TSVDUnfold_local

  TSVDUnfold_local *svdImpl = unfoldFull.Impl();
  if (!svdImpl) {
    std::cerr << "[Error] FindOptimalSvdK_Dvector: Could not access TSVDUnfold_local" << std::endl;
    return 6;
  }

  TH1D *hD = svdImpl->GetD();
  if (!hD) {
    std::cerr << "[Error] FindOptimalSvdK_Dvector: GetD() returned null" << std::endl;
    return 6;
  }

  std::cerr << "[Info] d-vector values:" << std::endl;
  const Int_t nD = hD->GetNbinsX();
  for (Int_t i = 1; i <= nD; ++i) {
    Double_t absDi = TMath::Abs(hD->GetBinContent(i));
    std::cerr << "[Info]   d_" << i << " = " << hD->GetBinContent(i)
              << "  |d_" << i << "| = " << absDi
              << (absDi > threshold ? "  > 1 (signal)" : "  <= 1 (noise)") << std::endl;
  }

  // Step 1: d-vector boundary — robust criterion: 3+ out of 4 consecutive below threshold
  // (avoids premature capping from isolated noise pockets, e.g. R=0.2 UE-sub)
  Int_t dvectorK = nD;
  const Int_t windowSize = 4;
  const Int_t noiseRequired = 3;
  for (Int_t i = 1; i <= nD - windowSize + 1; ++i) {
    Int_t noiseCount = 0;
    for (Int_t j = 0; j < windowSize; ++j) {
      if (TMath::Abs(hD->GetBinContent(i + j)) < threshold)
        noiseCount++;
    }
    if (noiseCount >= noiseRequired) {
      dvectorK = i - 1;
      break;
    }
  }
  if (dvectorK < 1) dvectorK = 1;

  const Int_t maxK = static_cast<Int_t>(svdUnfolded.size());

  // If d-vector found no clear boundary (dvectorK == nD), fall back to bilateral stability from k=2
  Int_t candidateK = dvectorK;
  if (dvectorK >= maxK) {
    std::cerr << "[Warning] d-vector found no clear signal/noise boundary (all |d_i| ~ 1)." << std::endl;
    std::cerr << "[Info] Falling back to bilateral stability scan from k=2." << std::endl;
    candidateK = 2;
  }
  std::cerr << "[Info] d-vector boundary: k=" << dvectorK
            << ", starting bilateral scan from k=" << candidateK << std::endl;

  // Step 2: Bilateral stability validation
  // k is on the plateau if BOTH |y(k-1)/y(k) - 1| < ε AND |y(k)/y(k+1) - 1| < ε.
  // This ensures the result is locally converged (not at the knee).
  std::cerr << "[Info] Bilateral stability validation (threshold = "
            << stabilityThreshold * 100 << "%):" << std::endl;

  // Helper: compute max |y(kA)/y(kB) - 1| across pT bins
  auto maxStepChange = [&](Int_t kA, Int_t kB) -> Double_t {
    Int_t idxA = kA - 1, idxB = kB - 1;
    if (idxA < 0 || idxA >= maxK || idxB < 0 || idxB >= maxK) return 999.;
    Double_t maxRel = 0.0;
    for (Int_t bin = 1; bin <= svdUnfolded[idxA]->GetNbinsX(); ++bin) {
      Double_t ptLow = svdUnfolded[idxA]->GetXaxis()->GetBinLowEdge(bin);
      Double_t ptHigh = svdUnfolded[idxA]->GetXaxis()->GetBinUpEdge(bin);
      if (ptLow < 5.0 || ptHigh > 200.0) continue; // exclude unreliable high-pT bins
      Double_t a = svdUnfolded[idxA]->GetBinContent(bin);
      Double_t b = svdUnfolded[idxB]->GetBinContent(bin);
      if (a > 0 && b > 0) {
        Double_t rel = TMath::Abs(b / a - 1.0);
        if (rel > maxRel) maxRel = rel;
      }
    }
    return maxRel;
  };

  // Scan for bilateral stability plateau: find ALL stable k, pick the LAST one.
  // MC closure doesn't have refolding chi2, so we pick the highest stable k
  // (least regularization bias) as a proxy for the best closure quality.
  // This avoids stopping too early (e.g. R=0.2 UE-sub: first stable at k=5
  // but closure quality keeps improving at higher k).
  Int_t firstStableK = -1;
  Int_t lastStableK = -1;
  for (Int_t k = candidateK; k < maxK; ++k) {
    Double_t backStep = maxStepChange(k - 1, k);  // |y(k-1)/y(k) - 1|
    Double_t fwdStep  = maxStepChange(k, k + 1);  // |y(k)/y(k+1) - 1|
    Bool_t backOK = (backStep < stabilityThreshold);
    Bool_t fwdOK  = (fwdStep  < stabilityThreshold);

    std::cerr << "[Info]   k=" << k
              << ": k-1->k = " << backStep * 100 << "%" << (backOK ? " [ok]" : "")
              << ",  k->k+1 = " << fwdStep * 100 << "%" << (fwdOK ? " [ok]" : "")
              << (backOK && fwdOK ? "  => STABLE" : "") << std::endl;

    if (backOK && fwdOK) {
      if (firstStableK < 0) firstStableK = k;
      lastStableK = k;
    }
  }

  Int_t finalK = candidateK;
  Bool_t foundStable = (lastStableK > 0);
  if (foundStable) {
    finalK = lastStableK;
  }

  if (!foundStable) {
    std::cerr << "[Warning] SVD bilateral stability never satisfied (threshold "
              << stabilityThreshold * 100 << "%). Using k=" << finalK
              << ". Consider relaxing threshold." << std::endl;
  } else if (lastStableK != firstStableK) {
    std::cerr << "[Info] Bilateral stability: first stable k=" << firstStableK
              << ", last stable k=" << lastStableK
              << " — using last (least regularization bias)" << std::endl;
  } else {
    std::cerr << "[Info] d-vector k=" << candidateK
              << " confirmed on plateau at k=" << finalK << std::endl;
  }

  // Clamp to [2, 20]
  if (finalK < 2) finalK = 2;
  if (finalK > 25) finalK = 25;

  std::cerr << "[Info] Optimal SVD k (d-vector + bilateral stability): " << finalK << std::endl;
  return finalK;
}


// =============================================================================
// Plot functions
// =============================================================================

// Plot d-vector |d_i| vs index (Höcker & Kartvelishvili Fig. 1c/2c)
void PlotDvector(RooUnfoldResponse *response, TH1D *hRecoData, Int_t optimalK,
                 const char *outDir, const char *outPdf, const char *setName = "") {
  if (!response || !hRecoData)
    return;

  const Double_t noiseExpectation = TMath::Sqrt(2.0 / TMath::Pi()); // ~0.798

  // Run SVD with k=nBins to get the full d-vector
  const Int_t nBins = hRecoData->GetNbinsX();
  RooUnfoldSvd unfoldFull(response, hRecoData, nBins);
  unfoldFull.Hreco();

  TSVDUnfold_local *svdImpl = unfoldFull.Impl();
  if (!svdImpl) return;
  TH1D *hD = svdImpl->GetD();
  if (!hD) return;

  // Build TGraph with integer x-values (hD bin centers may not be integers)
  const Int_t nD = hD->GetNbinsX();
  std::vector<Double_t> xVals(nD), yVals(nD);
  Double_t yMin = 1e30, yMax = 0;
  for (Int_t i = 0; i < nD; ++i) {
    xVals[i] = i + 1;  // SVD index 1, 2, 3, ...
    yVals[i] = TMath::Abs(hD->GetBinContent(i + 1));
    if (yVals[i] > 0 && yVals[i] < yMin) yMin = yVals[i];
    if (yVals[i] > yMax) yMax = yVals[i];
  }
  // Ensure all points visible: set minimum to half the smallest nonzero value
  if (yMin > 0) yMin *= 0.5;
  else yMin = 0.01;

  TGraph *gD = new TGraph(nD, xVals.data(), yVals.data());
  gD->SetMarkerStyle(20);
  gD->SetMarkerSize(1.0);
  gD->SetMarkerColor(kBlack);
  gD->SetLineColor(kBlack);

  TCanvas *c = new TCanvas(Form("cDvector_%s", setName), "d-vector", 600, 500);
  c->SetLogy();
  c->SetGridy();

  gD->SetTitle("");
  gD->GetXaxis()->SetTitle("SVD index #it{i}");
  gD->GetYaxis()->SetTitle("|#it{d}_{i}|");
  gD->GetXaxis()->CenterTitle();
  gD->GetYaxis()->CenterTitle();
  gD->GetYaxis()->SetRangeUser(yMin, yMax * 5.0);
  gD->GetXaxis()->SetLimits(0.5, nD + 0.5);
  gD->Draw("AP");

  // Horizontal line at |d_i| = 1 (signal/noise boundary)
  TLine *lineThresh = new TLine(0.5, 1.0, nD + 0.5, 1.0);
  lineThresh->SetLineColor(kOrange + 1);
  lineThresh->SetLineStyle(7);
  lineThresh->SetLineWidth(2);
  lineThresh->Draw("same");

  // Horizontal line at noise expectation sqrt(2/pi)
  TLine *lineNoise = new TLine(0.5, noiseExpectation, nD + 0.5, noiseExpectation);
  lineNoise->SetLineColor(kRed);
  lineNoise->SetLineStyle(7);
  lineNoise->SetLineWidth(2);
  lineNoise->Draw("same");

  // Vertical line at chosen k
  TLine *lineK = new TLine(optimalK, yMin, optimalK, yMax * 3.0);
  lineK->SetLineColor(kBlue);
  lineK->SetLineStyle(7);
  lineK->SetLineWidth(2);
  lineK->Draw("same");

  TLegend *leg = new TLegend(0.45, 0.65, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->AddEntry(gD, "|d_{i}|", "p");
  leg->AddEntry(lineThresh, "Signal/noise boundary (|d_{i}| = 1)", "l");
  leg->AddEntry(lineNoise, Form("#sqrt{2/#pi} #approx %.2f", noiseExpectation), "l");
  leg->AddEntry(lineK, Form("Optimal k = %d", optimalK), "l");
  leg->Draw();

  gSystem->Exec(Form("mkdir -p %s", outDir));
  c->SaveAs(Form("%s/%s", outDir, outPdf));
}

// Plot bilateral stability QA: max step-change vs k for both directions
// Shows forward (k→k+1) and backward (k-1→k) max relative changes,
// with horizontal threshold line and vertical line at chosen k.
void PlotBilateralStability(const std::vector<TH1D *> &svdUnfolded, Int_t optimalK,
                             Double_t threshold, const char *outDir,
                             const char *outPdf, const char *setName = "",
                             const char *xTitle = "SVD regularization parameter #it{k}",
                             const char *optLabel = "k") {
  const Int_t maxK = static_cast<Int_t>(svdUnfolded.size());
  if (maxK < 3) return;

  // Helper: max |y(kA)/y(kB) - 1| across pT bins in [5, 100] GeV
  auto maxStepChange = [&](Int_t kA, Int_t kB) -> Double_t {
    Int_t idxA = kA - 1, idxB = kB - 1;
    if (idxA < 0 || idxA >= maxK || idxB < 0 || idxB >= maxK) return -1.;
    Double_t maxRel = 0.0;
    for (Int_t bin = 1; bin <= svdUnfolded[idxA]->GetNbinsX(); ++bin) {
      Double_t ptLow = svdUnfolded[idxA]->GetXaxis()->GetBinLowEdge(bin);
      Double_t ptHigh = svdUnfolded[idxA]->GetXaxis()->GetBinUpEdge(bin);
      if (ptLow < 5.0 || ptHigh > 200.0) continue;
      Double_t a = svdUnfolded[idxA]->GetBinContent(bin);
      Double_t b = svdUnfolded[idxB]->GetBinContent(bin);
      if (a > 0 && b > 0) {
        Double_t rel = TMath::Abs(b / a - 1.0);
        if (rel > maxRel) maxRel = rel;
      }
    }
    return maxRel;
  };

  // Build graphs: forward step (k→k+1) and backward step (k-1→k)
  std::vector<Double_t> kVals, fwdVals, backVals;
  for (Int_t k = 2; k < maxK; ++k) {
    Double_t fwd  = maxStepChange(k, k + 1);
    Double_t back = maxStepChange(k - 1, k);
    if (fwd < 0 || back < 0) continue;
    kVals.push_back(k);
    fwdVals.push_back(fwd * 100.0);   // percent
    backVals.push_back(back * 100.0);  // percent
  }
  if (kVals.empty()) return;

  TGraph *gFwd  = new TGraph(kVals.size(), kVals.data(), fwdVals.data());
  TGraph *gBack = new TGraph(kVals.size(), kVals.data(), backVals.data());

  gFwd->SetMarkerStyle(20);
  gFwd->SetMarkerSize(1.0);
  gFwd->SetMarkerColor(kRed);
  gFwd->SetLineColor(kRed);
  gFwd->SetLineWidth(2);

  gBack->SetMarkerStyle(21);
  gBack->SetMarkerSize(1.0);
  gBack->SetMarkerColor(kBlue);
  gBack->SetLineColor(kBlue);
  gBack->SetLineWidth(2);

  TCanvas *c = new TCanvas(Form("cStability_%s", setName), "Bilateral Stability", 700, 500);
  c->SetGridy();

  // Find y range
  Double_t yMax = 0;
  for (size_t i = 0; i < fwdVals.size(); ++i)
    yMax = TMath::Max(yMax, TMath::Max(fwdVals[i], backVals[i]));
  yMax = TMath::Max(yMax * 1.3, threshold * 100 * 2.0);

  gBack->SetTitle("");
  gBack->GetXaxis()->SetTitle(xTitle);
  gBack->GetYaxis()->SetTitle("max |#it{y}(#it{n}#pm1)/#it{y}(#it{n}) #minus 1| (%)");
  gBack->GetXaxis()->CenterTitle();
  gBack->GetYaxis()->CenterTitle();
  gBack->GetYaxis()->SetRangeUser(0, yMax);
  gBack->GetXaxis()->SetLimits(kVals.front() - 0.5, kVals.back() + 0.5);
  gBack->Draw("ALP");
  gFwd->Draw("LP same");

  // Threshold line
  TLine *lineThresh = new TLine(kVals.front() - 0.5, threshold * 100,
                                 kVals.back() + 0.5, threshold * 100);
  lineThresh->SetLineColor(kGreen + 2);
  lineThresh->SetLineStyle(7);
  lineThresh->SetLineWidth(2);
  lineThresh->Draw("same");

  // Vertical line at chosen k/iter
  TLine *lineK = new TLine(optimalK, 0, optimalK, yMax * 0.9);
  lineK->SetLineColor(kBlack);
  lineK->SetLineStyle(2);
  lineK->SetLineWidth(2);
  lineK->Draw("same");

  TLegend *leg = new TLegend(0.45, 0.65, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->AddEntry(gBack, "max |#it{y}(#it{n}#minus1)/#it{y}(#it{n}) #minus 1|", "lp");
  leg->AddEntry(gFwd,  "max |#it{y}(#it{n})/#it{y}(#it{n}+1) #minus 1|", "lp");
  leg->AddEntry(lineThresh, Form("#varepsilon = %.0f%%", threshold * 100), "l");
  leg->AddEntry(lineK, Form("Optimal %s = %d", optLabel, optimalK), "l");
  leg->Draw();

  gSystem->Exec(Form("mkdir -p %s", outDir));
  c->SaveAs(Form("%s/%s", outDir, outPdf));
}

// Plot Bayesian chi2/ndf vs iteration QA (analogous to SVD d-vector plot)
// Shows chi2/ndf vs truth at each iteration, with threshold line and optimal iter marker.
void PlotBayesChi2(TH1D *hTruth, const std::vector<TH1D *> &unfolded,
                    Double_t nEventsTruth, Double_t nEventsUnfolded,
                    Int_t optimalIter, Double_t chi2Threshold,
                    const char *outDir, const char *outPdf, const char *setName = "") {
  if (!hTruth || unfolded.empty()) return;

  TH1D *hTruthNorm = dynamic_cast<TH1D *>(hTruth->Clone("hTruthNorm_plot"));
  hTruthNorm->SetDirectory(nullptr);
  if (nEventsTruth > 0) hTruthNorm->Scale(1.0 / nEventsTruth, "width");

  const Int_t maxIter = static_cast<Int_t>(unfolded.size());
  std::vector<Double_t> iterVals, chi2Vals;

  for (Int_t i = 0; i < maxIter; ++i) {
    Int_t iter = i + 1;
    TH1D *hUnfNorm = dynamic_cast<TH1D *>(unfolded[i]->Clone("hUnfNorm_plot"));
    hUnfNorm->SetDirectory(nullptr);
    if (nEventsUnfolded > 0) hUnfNorm->Scale(1.0 / nEventsUnfolded, "width");

    double chi2 = 0.0;
    int ndf = 0;
    for (int bin = 1; bin <= hTruthNorm->GetNbinsX(); ++bin) {
      double ptLow = hTruthNorm->GetXaxis()->GetBinLowEdge(bin);
      double ptHigh = hTruthNorm->GetXaxis()->GetBinUpEdge(bin);
      if (ptLow < 5.0 || ptHigh > 200.0) continue;
      double truth = hTruthNorm->GetBinContent(bin);
      double unfVal = hUnfNorm->GetBinContent(bin);
      double err = hUnfNorm->GetBinError(bin);
      double truthErr = hTruthNorm->GetBinError(bin);
      double totalErr = TMath::Sqrt(err * err + truthErr * truthErr);
      if (truth > 0 && unfVal > 0 && totalErr > 0) {
        double diff = unfVal - truth;
        chi2 += (diff * diff) / (totalErr * totalErr);
        ndf++;
      }
    }
    double reducedChi2 = (ndf > 0) ? chi2 / ndf : 0.0;
    iterVals.push_back(iter);
    chi2Vals.push_back(reducedChi2);
    delete hUnfNorm;
  }
  delete hTruthNorm;

  if (iterVals.empty()) return;

  TGraph *gChi2 = new TGraph(iterVals.size(), iterVals.data(), chi2Vals.data());
  gChi2->SetMarkerStyle(20);
  gChi2->SetMarkerSize(1.0);
  gChi2->SetMarkerColor(kBlack);
  gChi2->SetLineColor(kBlack);
  gChi2->SetLineWidth(2);

  TCanvas *c = new TCanvas(Form("cBayesChi2_%s", setName), "Bayes chi2/ndf", 600, 500);
  c->SetGridy();

  Double_t yMax = 0;
  for (auto v : chi2Vals) yMax = TMath::Max(yMax, v);
  yMax = TMath::Max(yMax * 1.3, chi2Threshold * 3.0);

  gChi2->SetTitle("");
  gChi2->GetXaxis()->SetTitle("Bayesian iteration");
  gChi2->GetYaxis()->SetTitle("#chi^{2}/ndf (vs truth)");
  gChi2->GetXaxis()->CenterTitle();
  gChi2->GetYaxis()->CenterTitle();
  gChi2->GetYaxis()->SetRangeUser(0, yMax);
  gChi2->GetXaxis()->SetLimits(0.5, maxIter + 0.5);
  gChi2->Draw("ALP");

  // Threshold line
  TLine *lineThresh = new TLine(0.5, chi2Threshold, maxIter + 0.5, chi2Threshold);
  lineThresh->SetLineColor(kOrange + 1);
  lineThresh->SetLineStyle(7);
  lineThresh->SetLineWidth(2);
  lineThresh->Draw("same");

  // chi2/ndf = 1 reference line
  TLine *lineOne = new TLine(0.5, 1.0, maxIter + 0.5, 1.0);
  lineOne->SetLineColor(kRed);
  lineOne->SetLineStyle(7);
  lineOne->SetLineWidth(2);
  lineOne->Draw("same");

  // Vertical line at optimal iter
  TLine *lineOpt = new TLine(optimalIter, 0, optimalIter, yMax * 0.9);
  lineOpt->SetLineColor(kBlue);
  lineOpt->SetLineStyle(7);
  lineOpt->SetLineWidth(2);
  lineOpt->Draw("same");

  TLegend *leg = new TLegend(0.45, 0.65, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->AddEntry(gChi2, "#chi^{2}/ndf", "lp");
  leg->AddEntry(lineThresh, Form("Threshold = %.1f", chi2Threshold), "l");
  leg->AddEntry(lineOne, "#chi^{2}/ndf = 1", "l");
  leg->AddEntry(lineOpt, Form("Optimal iter = %d", optimalIter), "l");
  leg->Draw();

  gSystem->Exec(Form("mkdir -p %s", outDir));
  c->SaveAs(Form("%s/%s", outDir, outPdf));
}

// Plot deviation ratio (optimal ± 1) using Filipad2
// IMPORTANT: Use separate nEvents for truth and unfolded
// setName: unique canvas name per input set
void PlotDeviationRatio(TH1D *hTruth, const std::vector<TH1D *> &unfolded,
                        Int_t optimalParam, const char *method,
                        Double_t nEventsTruth, Double_t nEventsUnfolded,
                        const char *outDir, const char *outPdf,
                        const char *setName = "") {
  // Clone and normalize truth
  TH1D *hTruthPlot = dynamic_cast<TH1D *>(hTruth->Clone("hTruth_plot"));
  hTruthPlot->SetDirectory(nullptr);
  if (nEventsTruth > 0)
    hTruthPlot->Scale(1.0 / nEventsTruth, "width");

  Filipad2 *pad = new Filipad2(Form("%s_deviation_%s", method, setName), 1, 2,
                               0.4, 100, 50, 0.7, 1, 1);
  pad->Draw();
  TPad *top = pad->GetPad(1);
  TPad *bottom = pad->GetPad(2);
  optFili(*top, 1, 1, 0, 1);
  optFili(*bottom, 1, 1, 0, 0);

  // Top pad: yields
  top->cd();
  hset(*hTruthPlot, JetPtTitleX, JetPtTitleY, 0.9, 1.3, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hTruthPlot->SetMarkerColor(kBlack);
  hTruthPlot->SetLineColor(kBlack);
  hTruthPlot->SetMarkerStyle(20);
  hTruthPlot->GetXaxis()->SetRangeUser(5., 200.);
  hTruthPlot->Draw("pe");

  TLegend *leg =
      new TLegend(0.55, 0.55, 0.88, 0.88,
                  Form("%s = %d #pm 1", method, optimalParam), "brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hTruthPlot, "MC truth", "pl");

  int colors[3] = {kRed + 1, kBlue + 1, kGreen + 2};
  std::vector<int> params;
  if (optimalParam > 2)
    params.push_back(optimalParam - 1);
  params.push_back(optimalParam);
  int maxParam = static_cast<int>(unfolded.size());
  if (optimalParam < maxParam)
    params.push_back(optimalParam + 1);

  std::vector<TH1D *> unfoldedPlot;
  for (size_t idx = 0; idx < params.size(); ++idx) {
    int param = params[idx];
    int i = param - 1; // 0-based index
    if (i < 0 || i >= static_cast<int>(unfolded.size()))
      continue;

    TH1D *h =
        dynamic_cast<TH1D *>(unfolded[i]->Clone(Form("hUnfolded_%d", param)));
    h->SetDirectory(nullptr);
    // Normalize unfolded with RECO nEvents
    if (nEventsUnfolded > 0)
      h->Scale(1.0 / nEventsUnfolded, "width");

    h->SetMarkerColor(colors[idx % 3]);
    h->SetLineColor(colors[idx % 3]);
    h->SetMarkerStyle(20);
    h->GetXaxis()->SetRangeUser(5., 200.);
    h->Draw("pe same");

    const char *label = (param == optimalParam)
                            ? Form("%s = %d (optimal)", method, param)
                            : Form("%s = %d", method, param);
    leg->AddEntry(h, label, "pl");
    unfoldedPlot.push_back(h);
  }
  leg->Draw();

  // Bottom pad: ratios
  bottom->cd();
  bottom->SetTicks(1, 1);

  std::vector<TH1D *> ratios;
  for (size_t idx = 0; idx < unfoldedPlot.size(); ++idx) {
    TH1D *ratio = CalculateRatioByPt(unfoldedPlot[idx], hTruthPlot);
    ratio->SetName(Form("hRatio_%d", params[idx]));
    ratios.push_back(ratio);
  }

  if (!ratios.empty()) {
    hset(*ratios[0], JetPtTitleX, "unfolded / truth", 1.1, 0.9, 0.06, 0.06,
         0.01, 0.01, 0.06, 0.06, 510, 505);

    for (size_t idx = 0; idx < ratios.size(); ++idx) {
      ratios[idx]->SetMarkerColor(colors[idx % 3]);
      ratios[idx]->SetLineColor(colors[idx % 3]);
      ratios[idx]->SetMarkerStyle(20);
      ratios[idx]->GetXaxis()->SetRangeUser(5., 200.);
      ratios[idx]->GetYaxis()->SetRangeUser(0.48, 1.52);
      ratios[idx]->Draw(idx == 0 ? "pe" : "pe same");
    }

    TLine *line = new TLine(5., 1., 200., 1.);
    line->SetLineStyle(2);
    line->Draw("same");
  }

  gSystem->Exec(Form("mkdir -p %s", outDir));
  pad->C->SaveAs(Form("%s/%s", outDir, outPdf));

  // Cleanup - disabled for ROOT browser QA
  // delete hTruthPlot;
  // for (auto *h : unfoldedPlot) delete h;
  // for (auto *r : ratios) delete r;
}

// Plot Refolding test (setName for unique canvas)
void PlotRefoldingTest(RooUnfoldResponse *response, TH1D *hReco,
                       TH1D *hUnfolded, Int_t optimalK, Double_t nEvents,
                       const char *outDir, const char *outPdf,
                       const char *setName = "") {
  // Clone for plotting
  TH1D *hRecoPlot = dynamic_cast<TH1D *>(hReco->Clone("hReco_refold"));
  hRecoPlot->SetDirectory(nullptr);

  // Debug: print binning info
  std::cerr << "[Debug] Refolding Test:" << std::endl;
  std::cerr << "  hReco: nbins=" << hReco->GetNbinsX()
            << ", xmin=" << hReco->GetXaxis()->GetXmin()
            << ", xmax=" << hReco->GetXaxis()->GetXmax() << std::endl;
  std::cerr << "  hUnfolded: nbins=" << hUnfolded->GetNbinsX()
            << ", xmin=" << hUnfolded->GetXaxis()->GetXmin()
            << ", xmax=" << hUnfolded->GetXaxis()->GetXmax() << std::endl;
  std::cerr << "  Response Htruth: nbins=" << response->Htruth()->GetNbinsX()
            << ", xmin=" << response->Htruth()->GetXaxis()->GetXmin()
            << ", xmax=" << response->Htruth()->GetXaxis()->GetXmax()
            << std::endl;
  std::cerr << "  Response Hmeasured: nbins="
            << response->Hmeasured()->GetNbinsX()
            << ", xmin=" << response->Hmeasured()->GetXaxis()->GetXmin()
            << ", xmax=" << response->Hmeasured()->GetXaxis()->GetXmax()
            << std::endl;
  std::cerr << "  Response Hfakes: nbins=" << response->Hfakes()->GetNbinsX()
            << ", xmin=" << response->Hfakes()->GetXaxis()->GetXmin()
            << ", xmax=" << response->Hfakes()->GetXaxis()->GetXmax()
            << std::endl;

  if (nEvents > 0)
    hRecoPlot->Scale(1.0 / nEvents, "width");

  // CRITICAL: ApplyToTruth expects truth-space histogram with EXACT same axis
  // as response->Htruth() RooUnfold Hreco() may return different axis; mismatch
  // causes wrong refolded shape
  TH1D *hUnfoldedAligned = CopyOntoTruthAxis(hUnfolded, response->Htruth());
  if (!hUnfoldedAligned) {
    std::cerr << "[Error] CopyOntoTruthAxis failed" << std::endl;
    return;
  }
  std::cerr << "  hUnfoldedAligned: nbins=" << hUnfoldedAligned->GetNbinsX()
            << " (same as Htruth)" << std::endl;

  // Refold by explicit matrix multiplication: refolded_reco[i] = sum_j
  // R(i,j)*truth[j] RM: Hresponse() X=reco (20 bins), Y=truth (25 bins). No
  // ApplyToTruth (convention/axis issues).
  TH1D *hRefolded = RefoldManual(response, hUnfoldedAligned);
  if (!hRefolded) {
    std::cerr << "[Error] RefoldManual failed" << std::endl;
    return;
  }
  std::cerr << "  hRefolded (manual, reco axis): nbins="
            << hRefolded->GetNbinsX()
            << ", xmin=" << hRefolded->GetXaxis()->GetXmin()
            << ", xmax=" << hRefolded->GetXaxis()->GetXmax() << std::endl;

  if (nEvents > 0)
    hRefolded->Scale(1.0 / nEvents, "width");

  // Plot
  Filipad2 *pad = new Filipad2(Form("RefoldingTest_%s", setName), 1, 2, 0.4,
                               100, 50, 0.7, 1, 1);
  pad->Draw();
  TPad *top = pad->GetPad(1);
  TPad *bottom = pad->GetPad(2);
  optFili(*top, 1, 1, 0, 1);
  optFili(*bottom, 1, 1, 0, 0);

  top->cd();
  hset(*hRecoPlot, JetPtTitleX, JetPtTitleY, 0.9, 1.3, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hRecoPlot->SetMarkerColor(kBlack);
  hRecoPlot->SetLineColor(kBlack);
  hRecoPlot->SetMarkerStyle(20);
  hRecoPlot->GetXaxis()->SetRangeUser(5., 200.);
  hRecoPlot->Draw("pe");

  hRefolded->SetMarkerColor(kRed + 1);
  hRefolded->SetLineColor(kRed + 1);
  hRefolded->SetMarkerStyle(24);
  hRefolded->Draw("pe same");

  TLegend *leg =
      new TLegend(0.55, 0.65, 0.88, 0.88,
                  Form("Refolding test (k = %d)", optimalK), "brNDC");
  leg->SetBorderSize(0);
  leg->AddEntry(hRecoPlot, "Measured reco", "pl");
  leg->AddEntry(hRefolded, "Refolded", "pl");
  leg->Draw();

  bottom->cd();
  bottom->SetTicks(1, 1);

  TH1D *ratio = CalculateRatioByPt(hRefolded, hRecoPlot);
  hset(*ratio, JetPtTitleX, "refolded / measured", 1.1, 0.9, 0.06, 0.06, 0.01,
       0.01, 0.06, 0.06, 510, 505);
  ratio->SetMarkerColor(kRed + 1);
  ratio->SetLineColor(kRed + 1);
  ratio->SetMarkerStyle(20);
  ratio->GetXaxis()->SetRangeUser(5., 200.);
  ratio->GetYaxis()->SetRangeUser(0.48, 1.52);
  ratio->Draw("pe");

  TLine *line = new TLine(5., 1., 200., 1.);
  line->SetLineStyle(2);
  line->Draw("same");

  gSystem->Exec(Form("mkdir -p %s", outDir));
  pad->C->SaveAs(Form("%s/%s", outDir, outPdf));

  // Cleanup - disabled for ROOT browser QA
  // delete hRecoPlot;
  // delete hRefolded;
  // delete ratio;
}

// =============================================================================
// Main function
// =============================================================================

void DrawMcClosureTest() {
  // Set up log file: tee all cerr and cout to plots/MCClosureTest/log.txt
  gSystem->Exec("mkdir -p plots/MCClosureTest");
  std::ofstream logFile("plots/MCClosureTest/log.txt");
  TeeBuf teeCerr(std::cerr.rdbuf(), logFile.rdbuf());
  TeeBuf teeCout(std::cout.rdbuf(), logFile.rdbuf());
  std::streambuf *origCerr = std::cerr.rdbuf(&teeCerr);
  std::streambuf *origCout = std::cout.rdbuf(&teeCout);

  std::cerr << "========================================" << std::endl;
  std::cerr << "MC Closure Test" << std::endl;
  std::cerr << "========================================" << std::endl;

  // Define input MC sets
  std::vector<MCInputSet> inputSets = {
      // LHC25a2b
      {"../../../jets/AnalysisResults/McClosure/"
        "merged_Response_LHC24f3c_515446.root",
        "../../../jets/AnalysisResults/McClosure/"
        "merged_MCreco_LHC24f3c_515446.root",
        "LHC24f3c_tunerA", "jet-spectra-charged_id34413", "515446", "jet-spectra-charged_id37200"},
      {"../../../jets/AnalysisResults/McClosure/"
       "merged_Response_LHC25a2b_516969.root",
       "../../../jets/AnalysisResults/McClosure/"
       "merged_MCreco_LHC25a2b_516969.root",
       "LHC25a2b_tunerA", "jet-spectra-charged_Nmax1p5", "516969_Nmax1p5"},
      {"../../../jets/AnalysisResults/McClosure/"
       "merged_Response_LHC25a2b_516969.root",
       "../../../jets/AnalysisResults/McClosure/"
       "merged_MCreco_LHC25a2b_516969.root",
       "LHC25a2b_tunerA_Nmax4", "jet-spectra-charged", "516969_Nmax4"},
      {"../../../jets/AnalysisResults/McClosure/"
       "merged_Response_LHC25a2b_Btune1p2.root",
       "../../../jets/AnalysisResults/McClosure/"
       "merged_MCreco_LHC25a2b_Btune1p2.root",
       "LHC25a2b_tunerB1p2", "jet-spectra-charged", "596838"},
      {"../../../jets/AnalysisResults/McClosure/"
       "merged_Response_LHC25a2b_Btune1p5.root",
       "../../../jets/AnalysisResults/McClosure/"
       "merged_MCreco_LHC25a2b_Btune1p5.root",
       "LHC25a2b_tunerB1p5", "jet-spectra-charged", "596836"},
      {"../../../jets/AnalysisResults/McClosure/"
       "merged_Response_LHC25a2b_Btune1p8.root",
       "../../../jets/AnalysisResults/McClosure/"
       "merged_MCreco_LHC25a2b_Btune1p8.root",
       "LHC25a2b_tunerB1p8", "jet-spectra-charged", "596837"},

      // LHC23k4h
      {"../../../jets/AnalysisResults/McClosure/"
       "merged_Response_LHC23k4h_606091_R01.root",
       "../../../jets/AnalysisResults/McClosure/"
       "merged_MCreco_LHC23k4h_606091_R01.root",
       "LHC23k4h_606091_R01", "jet-spectra-charged", "606091"},
      {"../../../jets/AnalysisResults/McClosure/"
       "merged_Response_LHC23k4h_606086_R02.root",
       "../../../jets/AnalysisResults/McClosure/"
       "merged_MCreco_LHC23k4h_606086_R02.root",
       "LHC23k4h_606086_R02", "jet-spectra-charged", "606086"},
      {"../../../jets/AnalysisResults/McClosure/"
       "merged_Response_LHC23k4h_606087_R03.root",
       "../../../jets/AnalysisResults/McClosure/"
       "merged_MCreco_LHC23k4h_606087_R03.root",
       "LHC23k4h_606087_R03", "jet-spectra-charged", "606087"},
      {"../../../jets/AnalysisResults/McClosure/"
       "merged_Response_LHC23k4h_606088_R04.root",
       "../../../jets/AnalysisResults/McClosure/"
       "merged_MCreco_LHC23k4h_606088_R04.root",
       "LHC23k4h_606088_R04", "jet-spectra-charged", "606088"},
      {"../../../jets/AnalysisResults/McClosure/"
       "merged_Response_LHC23k4h_606089_R05.root",
       "../../../jets/AnalysisResults/McClosure/"
       "merged_MCreco_LHC23k4h_606089_R05.root",
       "LHC23k4h_606089_R05", "jet-spectra-charged", "606089"},
      {"../../../jets/AnalysisResults/McClosure/"
       "merged_Response_LHC23k4h_606090_R06.root",
       "../../../jets/AnalysisResults/McClosure/"
       "merged_MCreco_LHC23k4h_606090_R06.root",
       "LHC23k4h_606090_R06", "jet-spectra-charged", "606090"},
      {"../../../jets/AnalysisResults/McClosure/"
       "merged_Response_LHC23k4h_606092_R07.root",
       "../../../jets/AnalysisResults/McClosure/"
       "merged_MCreco_LHC23k4h_606092_R07.root",
       "LHC23k4h_606092_R07", "jet-spectra-charged", "606092"},
  };

  // Collect optimal regularization parameters for output
  // Key = name, value = {runId, kSVD, iterBayes}
  struct OptResult {
    TString runId;
    int kSVD;
    int iterBayes;
  };
  std::map<TString, OptResult> optimalParams[2]; // indexed by iMode (0=NoUE, 1=UESub)

  // Process each input set
  for (const auto &input : inputSets) {
    std::cerr << "\n======== Processing: " << input.name
              << " ========" << std::endl;

    // Loop over both UE modes (non-UE and UE-subtracted)
    for (int iMode = 0; iMode < 2; ++iMode) {
    const UEMode &mode = kUEModes[iMode];

    std::cerr << "\n---- Mode: " << mode.label << " ----" << std::endl;

    // Per-dataset per-mode output directory
    TString datasetDir = Form("%s/%s", mode.outputDir, input.name.Data());
    gSystem->Exec(Form("mkdir -p %s", datasetDir.Data()));

    // Load histograms
    TH1D *hTruth = nullptr, *hReco = nullptr, *hPseudoData = nullptr,
         *hPseudoTruth = nullptr;
    TH2D *hResp2D = nullptr;
    Double_t nEventsUnfolded = 0, nEventsTruth = 0;

    if (!LoadMCInputSet(input, mode, hTruth, hReco, hResp2D, hPseudoData,
                        hPseudoTruth, nEventsUnfolded, nEventsTruth)) {
      std::cerr << "[Skip] " << input.name << " / " << mode.label
                << ": histograms not available" << std::endl;
      continue;
    }

    // MC Closure Test:
    // - Unfold reco file's MCD jet (hPseudoData) using response file's response
    // matrix
    // - Compare unfolded/nEventsUnfolded with hPseudoTruth/nEventsTruth (both
    // from reco file)
    TH1D *hTruthForClosure = hPseudoTruth ? hPseudoTruth : hTruth;

    std::cerr << "[Info] nEventsUnfolded (h_collisions): " << nEventsUnfolded
              << std::endl;
    std::cerr << "[Info] nEventsTruth (h_mcColl_counts): " << nEventsTruth
              << std::endl;

    // Build response (kernel-style) + correction factors
    TH1D *hPurity = nullptr;
    TH1D *hEfficiency = nullptr;
    RooUnfoldResponse *response =
        BuildResponse(hTruth, hReco, hResp2D, hPurity, hEfficiency);
    if (!response) {
      std::cerr << "[Error] Failed to build response" << std::endl;
      continue;
    }

    if (!hPurity || !hEfficiency) {
      std::cerr << "[Error] Failed to build purity/efficiency histograms"
                << std::endl;
      continue;
    }

    // Purity correction (fake subtraction) BEFORE unfolding
    TH1D *hPseudoDataMatched =
        dynamic_cast<TH1D *>(hPseudoData->Clone("hPseudoDataMatched"));
    hPseudoDataMatched->SetDirectory(nullptr);
    for (int i = 1; i <= hPseudoDataMatched->GetNbinsX(); ++i) {
      const double x = hPseudoData->GetBinContent(i);
      const double ex = hPseudoData->GetBinError(i);
      const double p = hPurity->GetBinContent(i);
      hPseudoDataMatched->SetBinContent(i, x * p);
      hPseudoDataMatched->SetBinError(i, ex * p);
    }

    // ============================================
    // Bayesian unfolding (iter = 1 to 20)
    // ============================================
    std::cerr << "\n[Step] Bayesian unfolding (iter = 1-20)..." << std::endl;
    std::vector<TH1D *> bayesUnfolded;

    for (int iter = 1; iter <= 20; ++iter) {
      RooUnfoldBayes unfold(response, hPseudoDataMatched, iter);
      unfold.SetVerbose(0);
      TH1D *h = dynamic_cast<TH1D *>(unfold.Hreco());
      if (h) {
        TH1D *hClone = dynamic_cast<TH1D *>(
            h->Clone(Form("hBayes_iter%d_%s_%s", iter, input.name.Data(), mode.label)));
        hClone->SetDirectory(nullptr);

        // Efficiency correction (miss) AFTER unfolding
        for (int b = 1; b <= hClone->GetNbinsX(); ++b) {
          const double pt = hClone->GetXaxis()->GetBinCenter(b);
          const int effBin = hEfficiency->GetXaxis()->FindBin(pt);
          double eff = 0.0;
          if (effBin >= 1 && effBin <= hEfficiency->GetNbinsX()) {
            eff = hEfficiency->GetBinContent(effBin);
          }
          if (eff > 0) {
            hClone->SetBinContent(b, hClone->GetBinContent(b) / eff);
            hClone->SetBinError(b, hClone->GetBinError(b) / eff);
          } else {
            hClone->SetBinContent(b, 0.0);
            hClone->SetBinError(b, 0.0);
          }
        }

        bayesUnfolded.push_back(hClone);
      }
    }

    // Find optimal iteration using chi2/ndf ≈ 1.0 criterion (closure truth available)
    Int_t optimalIter = FindOptimalBayesIter(
        hTruthForClosure, bayesUnfolded, nEventsTruth, nEventsUnfolded);

    // QA plots: Bayesian chi2/ndf and bilateral stability
    TString setTag = Form("%s_%s", input.name.Data(), mode.label);
    PlotBayesChi2(hTruthForClosure, bayesUnfolded, nEventsTruth, nEventsUnfolded,
                  optimalIter, 2.0, datasetDir.Data(), "Bayes_Chi2.pdf", setTag.Data());
    PlotBilateralStability(bayesUnfolded, optimalIter, 0.03,
                           datasetDir.Data(), "Bayes_BilateralStability.pdf", setTag.Data(),
                           "Bayesian iteration", "iter");

    // Plot Bayesian deviation ratio
    PlotDeviationRatio(hTruthForClosure, bayesUnfolded, optimalIter, "iter",
                       nEventsTruth, nEventsUnfolded,
                       datasetDir.Data(), "Bayes_deviation.pdf",
                       setTag.Data());

    // ============================================
    // SVD unfolding (k = 1 to 25)
    // ============================================
    std::cerr << "\n[Step] SVD unfolding (k = 1-25)..." << std::endl;
    std::vector<TH1D *> svdUnfolded;

    for (int k = 1; k <= 25; ++k) {
      RooUnfoldSvd unfoldSvd(response, hPseudoDataMatched, k);
      unfoldSvd.SetVerbose(0);
      TH1D *h = dynamic_cast<TH1D *>(unfoldSvd.Hreco());
      if (h) {
        TH1D *hClone = dynamic_cast<TH1D *>(
            h->Clone(Form("hSVD_k%d_%s_%s", k, input.name.Data(), mode.label)));
        hClone->SetDirectory(nullptr);

        // CRITICAL: Get statistical errors from covariance matrix for SVD
        TMatrixD covMatrix = unfoldSvd.Ereco();
        int nBins = hClone->GetNbinsX();
        bool useCovMatrix =
            (covMatrix.GetNrows() >= nBins && covMatrix.GetNcols() >= nBins);

        for (int i = 1; i <= nBins; ++i) {
          int matrixIdx = i - 1;
          if (useCovMatrix && matrixIdx < covMatrix.GetNrows()) {
            double covDiag = covMatrix(matrixIdx, matrixIdx);
            if (covDiag > 0) {
              hClone->SetBinError(i, TMath::Sqrt(covDiag));
            }
          }
        }

        // Efficiency correction (miss) AFTER unfolding
        for (int b = 1; b <= hClone->GetNbinsX(); ++b) {
          const double pt = hClone->GetXaxis()->GetBinCenter(b);
          const int effBin = hEfficiency->GetXaxis()->FindBin(pt);
          double eff = 0.0;
          if (effBin >= 1 && effBin <= hEfficiency->GetNbinsX()) {
            eff = hEfficiency->GetBinContent(effBin);
          }
          if (eff > 0) {
            hClone->SetBinContent(b, hClone->GetBinContent(b) / eff);
            hClone->SetBinError(b, hClone->GetBinError(b) / eff);
          } else {
            hClone->SetBinContent(b, 0.0);
            hClone->SetBinError(b, 0.0);
          }
        }

        svdUnfolded.push_back(hClone);
      }
    }

    // Find optimal k using d-vector + bilateral stability
    Int_t optimalK = FindOptimalSvdK_Dvector(response, hPseudoDataMatched, svdUnfolded);
    std::cerr << "[Info] Using optimal k = " << optimalK
              << " (d-vector + bilateral stability)" << std::endl;

    // QA plots: SVD d-vector and bilateral stability
    PlotDvector(response, hPseudoDataMatched, optimalK,
                datasetDir.Data(), "SVD_Dvector.pdf", setTag.Data());
    PlotBilateralStability(svdUnfolded, optimalK, 0.03,
                           datasetDir.Data(), "SVD_BilateralStability.pdf", setTag.Data());

    // Plot SVD deviation ratio
    PlotDeviationRatio(hTruthForClosure, svdUnfolded, optimalK, "k",
                       nEventsTruth, nEventsUnfolded,
                       datasetDir.Data(), "SVD_deviation.pdf",
                       setTag.Data());

    // ============================================
    // Refolding test
    // ============================================
    std::cerr << "\n[Step] Refolding test..." << std::endl;
    if (optimalK >= 1 && optimalK <= static_cast<int>(svdUnfolded.size())) {
      // Refold in matched-reco space: need detected-truth = efficiency * truth
      TH1D *hUnfoldedDetectedForRefold =
          dynamic_cast<TH1D *>(svdUnfolded[optimalK - 1]->Clone(
              Form("hUnfoldedDetectedForRefold_%s_%s", input.name.Data(), mode.label)));
      hUnfoldedDetectedForRefold->SetDirectory(nullptr);
      for (int b = 1; b <= hUnfoldedDetectedForRefold->GetNbinsX(); ++b) {
        const double pt =
            hUnfoldedDetectedForRefold->GetXaxis()->GetBinCenter(b);
        const int effBin = hEfficiency->GetXaxis()->FindBin(pt);
        double eff = 0.0;
        if (effBin >= 1 && effBin <= hEfficiency->GetNbinsX()) {
          eff = hEfficiency->GetBinContent(effBin);
        }
        hUnfoldedDetectedForRefold->SetBinContent(
            b, hUnfoldedDetectedForRefold->GetBinContent(b) * eff);
        hUnfoldedDetectedForRefold->SetBinError(
            b, hUnfoldedDetectedForRefold->GetBinError(b) * eff);
      }

      PlotRefoldingTest(response, hPseudoDataMatched,
                        hUnfoldedDetectedForRefold, optimalK, nEventsUnfolded,
                        datasetDir.Data(), "Refolding.pdf",
                        setTag.Data());
    }

    // ============================================
    // Save MC closure uncertainty
    // ============================================
    std::cerr << "\n[Step] Saving MC closure uncertainty..." << std::endl;
    TFile *fOut =
        new TFile(Form("%s/MCClosure.root", datasetDir.Data()),
                  "RECREATE");

    if (optimalK >= 1 && optimalK <= static_cast<int>(svdUnfolded.size())) {
      // Save normalized unfolded
      TH1D *hUnfoldedNorm =
          dynamic_cast<TH1D *>(svdUnfolded[optimalK - 1]->Clone("hUnfolded"));
      hUnfoldedNorm->SetDirectory(fOut);
      if (nEventsUnfolded > 0)
        hUnfoldedNorm->Scale(1.0 / nEventsUnfolded, "width");

      // Save normalized truth
      TH1D *hTruthNorm =
          dynamic_cast<TH1D *>(hTruthForClosure->Clone("hTruth"));
      hTruthNorm->SetDirectory(fOut);
      if (nEventsTruth > 0)
        hTruthNorm->Scale(1.0 / nEventsTruth, "width");

      // Calculate MC closure uncertainty: |ratio - 1| in percent
      TH1D *hUncertainty =
          dynamic_cast<TH1D *>(hTruthNorm->Clone("hMCClosureUncertainty"));
      hUncertainty->SetDirectory(fOut);
      hUncertainty->Reset();

      for (int i = 1; i <= hUncertainty->GetNbinsX(); ++i) {
        double truth = hTruthNorm->GetBinContent(i);
        double unfolded = hUnfoldedNorm->GetBinContent(i);
        if (truth > 0 && unfolded > 0) {
          double ratio = unfolded / truth;
          double deviation = TMath::Abs(ratio - 1.0) * 100.0; // in percent
          hUncertainty->SetBinContent(i, deviation);
        }
      }

      // Save optimal parameters for use in DrawUnfoldingSystematicUncertainty.C
      TParameter<int> paramSVDk("optimalSVDk", optimalK);
      paramSVDk.Write();
      TParameter<int> paramBayesIter("optimalBayesIter", optimalIter);
      paramBayesIter.Write();
      std::cerr << "[Info] Saved optimal parameters: SVD k=" << optimalK
                << ", Bayes iter=" << optimalIter << std::endl;

      // Collect for accumulated text file
      optimalParams[iMode][input.name] = {input.runId, optimalK, optimalIter};

      fOut->Write();
    }
    // fOut->Close();  // Keep open for ROOT browser QA

    std::cerr << "[Info] Results saved to: " << datasetDir << "/MCClosure.root" << std::endl;

    // Cleanup - disabled for ROOT browser QA
    // delete response;
    // delete hTruth;
    // delete hReco;
    // delete hResp2D;
    // delete hPseudoData;
    // if (hPseudoTruth) delete hPseudoTruth;
    // for (auto *h : bayesUnfolded) delete h;
    // for (auto *h : svdUnfolded) delete h;
    } // end iMode loop
  } // end inputSets loop

  // ============================================
  // Write accumulated optimal parameters to text file (one per UE mode)
  // Reads existing file first, merges (overwrites duplicates), writes back.
  // Format: tab-separated, easily parseable by other macros via
  //   ifstream >> name >> kSVD >> iterBayes
  // ============================================
  for (int iMode = 0; iMode < 2; ++iMode) {
    if (optimalParams[iMode].empty())
      continue;

    const UEMode &mode = kUEModes[iMode];
    gSystem->Exec(Form("mkdir -p %s", mode.outputDir));
    TString paramFile = Form("%s/OptimalRegularization.txt", mode.outputDir);

    // Read existing entries (if file exists)
    std::map<TString, OptResult> existing;
    {
      std::ifstream fin(paramFile.Data());
      if (fin.is_open()) {
        std::string line;
        while (std::getline(fin, line)) {
          if (line.empty() || line[0] == '#')
            continue;
          std::istringstream iss(line);
          std::string name, runId;
          int k, iter;
          if (iss >> name >> runId >> k >> iter) {
            existing[name.c_str()] = {runId.c_str(), k, iter};
          }
        }
        fin.close();
        std::cerr << "[Info] Read " << existing.size()
                  << " existing entries from " << paramFile << std::endl;
      }
    }

    // Merge: new results overwrite existing entries for the same name
    int nUpdated = 0, nNew = 0;
    for (const auto &p : optimalParams[iMode]) {
      if (existing.count(p.first)) {
        nUpdated++;
      } else {
        nNew++;
      }
      existing[p.first] = p.second;
    }

    // Write back the full merged file
    {
      std::ofstream fout(paramFile.Data());
      fout << "# OptimalRegularization.txt (" << mode.label << ")" << std::endl;
      fout << "# Generated by DrawMcClosureTest.C" << std::endl;
      fout << "# Name\tRunId\tkSVD\titerBayes" << std::endl;
      for (const auto &e : existing) {
        TString rid = e.second.runId.Length() > 0 ? e.second.runId : "-";
        fout << e.first << "\t" << rid << "\t" << e.second.kSVD << "\t"
             << e.second.iterBayes << std::endl;
      }
      fout.close();
    }

    std::cerr << "[Info] Wrote " << paramFile << " (" << existing.size()
              << " entries: " << nNew << " new, " << nUpdated << " updated)"
              << std::endl;
  }

  std::cerr << "\n========================================" << std::endl;
  std::cerr << "MC Closure Test completed!" << std::endl;
  std::cerr << "========================================" << std::endl;

  // Restore original streams
  std::cerr.rdbuf(origCerr);
  std::cout.rdbuf(origCout);
  logFile.close();
  std::cerr << "[Info] Log saved to plots/MCClosureTest/log.txt" << std::endl;
}
