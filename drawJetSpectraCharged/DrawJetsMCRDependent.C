#include "DrawJetsMCRDependentHelpers.h"
#include "TLegendEntry.h"
#include "TParameter.h"
#include <set>

// ============================================
// User Configuration (edit here)
// ============================================

// Unfolding regularization: DATA d-vector determines nominal k (standard ALICE practice)
// MC closure values from OptimalRegularization.txt are logged for cross-check only.
// Format: Name  RunId  kSVD  iterBayes
// Separate files for non-UE and UE-subtracted (DrawMcClosureTest.C writes to different dirs)
const char* kOptimalRegFileNonUE = "plots/MCClosureTest/OptimalRegularization.txt";
const char* kOptimalRegFileUE    = "plots/MCClosureTest_UESub/OptimalRegularization.txt";

// Temporary fallback k when mcRunNumber is NOT found in the txt file (loud warning)
// kSvdKFallback removed: d-vector auto-optimization replaces fixed fallback

// MB/JJ stitching for final Run 3 cross section
// - If enabled, build a stitched spectrum per R on the unfolded axis (true jet pT)
// - Hard switch: pT < switch -> MB, pT >= switch -> JJ
// - switch must match an edge in ptbinGen, otherwise stitching is aborted with an error
const bool kEnableMBJJStitching = true;
const Double_t kPtSwitchStitch = 20.0;

// When enabled, stitched spectra are written as additional UE plots
const bool kDrawStitchedCrossSection = true;

// R-ratio reference: sigma(R) / sigma(R_ref)
const Double_t kRRatioRef = 0.4;

// ============================================
// Cross-section normalization constants
// ============================================
// Year-dependent constants (ε_INEL→TVX, σ_vis) are now inside ComputeXsecNormFactor().
// Dataset year is auto-detected from RConfigSet name via DetectDatasetYear().
// Norm mode (kNormSoft / kNormDirect) is set via SetXsecNormMode().

// RCT label Y-bin in hLumiTVXafterBCcutsRCT (CBT_hadronPID = 5)
const Int_t kRCTYBin = 5;

// ============================================
// Internal caches (filled by DrawCrossSection)
// ============================================

// UE mode: store normalized Run 3 cross sections per R for MB/JJ and stitched.
// These are used for stitching plots and for R-ratio plots.
static std::map<double, TH1*> gUE_Xsec_MB;
static std::map<double, TH1*> gUE_Xsec_JJ;
static std::map<double, TH1*> gUE_Xsec_Stitched;

static std::map<double, TString> gUE_McRun_MB;
static std::map<double, TString> gUE_McRun_JJ;
static std::map<double, TH1*> gUE_POWHEG;  // POWHEG NLO per R (for UE AllR canvas)
static std::map<std::string, std::map<double, TH1*>> gUE_Models;  // standalone models per R (for UE AllR canvas)
static std::map<double, TH1*> gUE_MCTruth;  // MC particle-level cross-section per R (for UE AllR canvas)

static void ClearXsecCachesUE()
{
  for (auto& kv : gUE_Xsec_MB) {
    delete kv.second;
  }
  for (auto& kv : gUE_Xsec_JJ) {
    delete kv.second;
  }
  for (auto& kv : gUE_Xsec_Stitched) {
    delete kv.second;
  }
  for (auto& kv : gUE_POWHEG) {
    delete kv.second;
  }
  for (auto& mkv : gUE_Models) {
    for (auto& kv : mkv.second) delete kv.second;
  }
  for (auto& kv : gUE_MCTruth) {
    delete kv.second;
  }
  gUE_Xsec_MB.clear();
  gUE_Xsec_JJ.clear();
  gUE_Xsec_Stitched.clear();
  gUE_McRun_MB.clear();
  gUE_McRun_JJ.clear();
  gUE_POWHEG.clear();
  gUE_Models.clear();
  gUE_MCTruth.clear();
}

// ============================================
// Cross-section computation cache
// ============================================
// Bump kXsecCacheVersion when the unfolding/normalization algorithm changes.
// Set kForceRecompute = true to ignore existing cache and recompute from scratch.
static const int kXsecCacheVersion = 3;
static bool kForceRecompute = false;

// ============================================
// MC closure regularization (for cross-check logging only)
// ============================================
struct RegParams { int kSVD; int iterBayes; };

static std::map<TString, RegParams> LoadOptimalRegularization(const char* path)
{
  std::map<TString, RegParams> result;
  std::ifstream fin(path);
  if (!fin.is_open()) {
    std::cerr << "[Warning] Cannot open " << path
              << "; will use d-vector auto-optimization for all R values" << std::endl;
    return result;
  }
  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    std::string name, runId;
    int k, iter;
    if (iss >> name >> runId >> k >> iter) {
      result[runId.c_str()] = {k, iter};
    }
  }
  fin.close();
  std::cerr << "[Info] Loaded " << result.size()
            << " regularization entries from " << path << std::endl;
  return result;
}

static Int_t GetSvdKForRun(const TString& mcRunNumber,
                           const std::map<TString, RegParams>& regMap)
{
  auto it = regMap.find(mcRunNumber);
  if (it != regMap.end()) {
    std::cerr << "[Info] MC closure SVD k=" << it->second.kSVD
              << " for mcRun=" << mcRunNumber << " (cross-check only)" << std::endl;
    return it->second.kSVD;
  }
  // NOT FOUND — return -1 to trigger d-vector auto-optimization
  std::cerr << "[Info] mcRunNumber=" << mcRunNumber
            << " not in OptimalRegularization.txt; will use d-vector optimization" << std::endl;
  return -1;
}

static bool IsAxisEdge(const TAxis* ax, double x, double tol = 1e-6)
{
  if (!ax) {
    return false;
  }
  const int n = ax->GetNbins();
  for (int i = 1; i <= n; ++i) {
    const double low = ax->GetBinLowEdge(i);
    if (std::abs(low - x) < tol) {
      return true;
    }
  }
  // Also allow the upper edge
  const double up = ax->GetBinUpEdge(n);
  return (std::abs(up - x) < tol);
}

static bool HaveSameBinning(const TH1* a, const TH1* b, double tol = 1e-9)
{
  if (!a || !b) {
    return false;
  }
  const TAxis* axA = a->GetXaxis();
  const TAxis* axB = b->GetXaxis();
  if (!axA || !axB) {
    return false;
  }
  if (axA->GetNbins() != axB->GetNbins()) {
    return false;
  }
  const int n = axA->GetNbins();
  for (int i = 1; i <= n; ++i) {
    const double aLow = axA->GetBinLowEdge(i);
    const double bLow = axB->GetBinLowEdge(i);
    if (std::abs(aLow - bLow) > tol) {
      return false;
    }
  }
  const double aUp = axA->GetBinUpEdge(n);
  const double bUp = axB->GetBinUpEdge(n);
  return (std::abs(aUp - bUp) <= tol);
}

static TH1* StitchMBJJByPtSwitch(const TH1* hMB, const TH1* hJJ, double ptSwitch, const char* name)
{
  if (!hMB || !hJJ) {
    return nullptr;
  }
  if (!HaveSameBinning(hMB, hJJ)) {
    std::cerr << "[Error] Cannot stitch: MB/JJ histograms have different binning." << std::endl;
    return nullptr;
  }
  if (!IsAxisEdge(hMB->GetXaxis(), ptSwitch)) {
    std::cerr << "[Error] Cannot stitch: pT_switch=" << ptSwitch << " is not an edge of the unfolded axis (ptbinGen)." << std::endl;
    std::cerr << "[Error] Please set pT_switch to one of the ptbinGen edges (GetPtbinGen())." << std::endl;
    return nullptr;
  }

  TH1* out = (TH1*)hMB->Clone(name);
  out->SetDirectory(0);
  out->SetStats(0);

  for (int i = 1; i <= out->GetNbinsX(); ++i) {
    const double low = out->GetXaxis()->GetBinLowEdge(i);
    const TH1* src = (low < ptSwitch) ? hMB : hJJ;
    out->SetBinContent(i, src->GetBinContent(i));
    out->SetBinError(i, src->GetBinError(i));
  }
  return out;
}

///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// R-dependent comparison macro  /////////
////////// author: Joonsuk Bae           /////////
////////// E-mail: jbae@cern.ch          /////////
////////// Last Modified: 2026           /////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

// Draw histogram with custom style (MC: empty square, Data: filled circle)
void DrawHistWithStyle(TH1* hist, Color_t color, bool isData, bool isFirst = false) {
  if (!hist) return;
  
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetLineWidth(2);
  
  if (isData) {
    hist->SetMarkerStyle(20);  // Filled circle
    hist->SetFillStyle(0);
  } else {
    hist->SetMarkerStyle(25);  // Empty square
    hist->SetFillStyle(0);
  }
  
  hist->SetMarkerSize(0.8);
  
  if (isFirst) {
    hist->Draw("pe");
  } else {
    hist->Draw("pesame");
  }
}

// ============================================
// Main plotting functions
// ============================================

void DrawTrackObservables() {
  if (!GetDrawTrack()) {
    std::cout << "[DrawTrackObservables] DRAW_TRACK is false, skipping..." << std::endl;
    return;
  }

  std::vector<RConfig> rConfigs = GetRConfigs();
  TString outputDir = GetOutputDir();

  // Track properties are R-independent — use R=0.4 files as representative.
  const RConfig* repConfig = nullptr;
  for (const auto& cfg : rConfigs) {
    if (TMath::Abs(cfg.R - 0.4) < 0.01) { repConfig = &cfg; break; }
  }
  if (!repConfig) {
    // Fallback: use the first available config
    if (rConfigs.empty()) {
      std::cerr << "[DrawTrackObservables] No R configs available. Skipping." << std::endl;
      return;
    }
    repConfig = &rConfigs[0];
  }

  std::cout << "[DrawTrackObservables] Starting track observables (using files from R="
            << repConfig->R << " as representative)..." << std::endl;
  std::vector<TString> observables = {"pt", "eta", "phi"};

  for (const auto& obs : observables) {
    std::cout << "[DrawTrackObservables] Processing " << obs << "..." << std::endl;

    // Create pad
    Filipad2* pad = new Filipad2(Form("Track%s", obs.Data()), ++GetNN(), 2, 0.4, 100, 50, 0.7, 1, 1);
    pad->Draw();
    TPad* mainPad = pad->GetPad(1);
    TPad* ratioPad = pad->GetPad(2);
    if (!mainPad || !ratioPad) {
      std::cerr << "[Error] Failed to get pads for " << obs << std::endl;
      continue;
    }
    int logy = (obs == "pt") ? 1 : 0;
    optFili(*mainPad, 0, logy, 0, 1);
    optFili(*ratioPad, 0, 0, 0, 0);

    TLegend* leg = new TLegend(0.55, 0.7, 0.9, 0.95, NULL, "brNDC");
    leg->SetTextSize(0.045);
    leg->SetBorderSize(0);

    mainPad->cd();

    // --- MC ---
    TString mcFile = repConfig->GetMcFile();
    Double_t nevtsMC = GetNORMEVENTS() ? Nevents(mcFile.Data(), repConfig->mcDir.Data(), GetEventObj(), 0, repConfig->isJJ) : 1.0;

    TH1* hMC = nullptr;
    if (obs == "pt") {
      hMC = DrawTrackPt(mcFile.Data(), "MC", nevtsMC, leg, kRed+1, repConfig->mcDir.Data());
    } else if (obs == "eta") {
      hMC = DrawTrackEta(mcFile.Data(), "MC", nevtsMC, leg, kRed+1, repConfig->mcDir.Data());
    } else if (obs == "phi") {
      hMC = DrawTrackPhi(mcFile.Data(), "MC", nevtsMC, leg, kRed+1, repConfig->mcDir.Data());
    }
    if (hMC) {
      hMC->SetMarkerStyle(25);  // Empty square for MC
      hMC->SetFillStyle(0);
    }

    // --- Data ---
    TString dataFile = repConfig->GetDataFile();
    Double_t nevtsData = GetNORMEVENTS() ? Nevents(dataFile.Data(), repConfig->dataDir.Data(), GetEventObj(), 0) : 1.0;

    TH1* hData = nullptr;
    if (obs == "pt") {
      hData = DrawTrackPt(dataFile.Data(), "Data", nevtsData, leg, kBlack, repConfig->dataDir.Data());
    } else if (obs == "eta") {
      hData = DrawTrackEta(dataFile.Data(), "Data", nevtsData, leg, kBlack, repConfig->dataDir.Data());
    } else if (obs == "phi") {
      hData = DrawTrackPhi(dataFile.Data(), "Data", nevtsData, leg, kBlack, repConfig->dataDir.Data());
    }
    if (hData) {
      hData->SetMarkerStyle(20);  // Filled circle for Data
      hData->SetFillStyle(1001);
    }

    mainPad->cd();
    leg->Draw();

    // --- Ratio pad: MC / Data ---
    ratioPad->cd();
    if (hMC && hData) {
      double ratioYMin = (obs == "pt") ? 0.0 : 0.5;
      double ratioYMax = (obs == "pt") ? 3.0 : 2.0;
      TString xAxisTitle;
      if (obs == "pt")       xAxisTitle = GetTrackPtTitleX();
      else if (obs == "eta") xAxisTitle = GetTrackEtaTitleX();
      else if (obs == "phi") xAxisTitle = GetTrackPhiTitleX();

      DrawRatio(Form("Ratio_%s", obs.Data()),
                hData, hMC, xAxisTitle, "MC / Data", kRed+1, ratioYMin, ratioYMax);
    }

    if (GetDRAWPLOTS()) {
      if (pad->C) {
        TString outputPath = Form("%s/Track%s.pdf", outputDir.Data(), obs.Data());
        pad->C->Print(outputPath.Data());
      }
    }

    std::cout << "[DrawTrackObservables] Completed " << obs << std::endl;
  }

  std::cerr << "[DrawTrackObservables] All track observables completed!" << std::endl;

  // ========== Track pT per R: Data vs MC with MC/Data ratio ==========
  // Each R config draws its own Data + MC pair; ratio pad shows MC/Data per R.
  {
    Filipad2* padCons = new Filipad2("TrackPt_AllR", ++GetNN(), 2, 0.4, 100, 50, 0.7, 1, 1);
    padCons->Draw();
    TPad* mainPadCons = padCons->GetPad(1);
    TPad* ratioPadCons = padCons->GetPad(2);
    optFili(*mainPadCons, 0, 1, 0, 1);
    optFili(*ratioPadCons, 0, 0, 0, 0);

    TLegend* legCons = new TLegend(0.45, 0.50, 0.92, 0.95, NULL, "brNDC");
    legCons->SetTextSize(0.035);
    legCons->SetBorderSize(0);

    mainPadCons->cd();

    std::vector<TH1*> dataHists;
    std::vector<TH1*> mcHists;
    bool firstDraw = true;

    for (size_t i = 0; i < rConfigs.size(); ++i) {
      const auto& cfg = rConfigs[i];
      int rIdx = (int)(cfg.R * 10 + 0.5);
      Color_t col = GetColorForR(cfg.R);

      // --- Data ---
      TH1* hData = nullptr;
      {
        auto file = TFile::Open(cfg.GetDataFile().Data(), "READ");
        if (file && !file->IsZombie()) {
          hData = (TH1*)file->Get(Form("%s/%s", cfg.dataDir.Data(), GetTrackPtObj()));
          if (hData) {
            hData = (TH1*)hData->Clone(Form("TrackPtData_R%02d", rIdx));
            hData->SetDirectory(0);
          }
          file->Close();
        }
      }
      if (hData && GetREBINON()) {
        TH1* rb = hData->Rebin(GetNTrackptbin(), Form("TrackPtDataReb_R%02d", rIdx), GetTrackptbin());
        if (rb) { rb->SetDirectory(0); delete hData; hData = rb; }
      }

      // --- MC ---
      TH1* hMC = nullptr;
      {
        auto file = TFile::Open(cfg.GetMcFile().Data(), "READ");
        if (file && !file->IsZombie()) {
          hMC = (TH1*)file->Get(Form("%s/%s", cfg.mcDir.Data(), GetTrackPtObj()));
          if (hMC) {
            hMC = (TH1*)hMC->Clone(Form("TrackPtMC_R%02d", rIdx));
            hMC->SetDirectory(0);
          }
          file->Close();
        }
      }
      if (hMC && GetREBINON()) {
        TH1* rb = hMC->Rebin(GetNTrackptbin(), Form("TrackPtMCReb_R%02d", rIdx), GetTrackptbin());
        if (rb) { rb->SetDirectory(0); delete hMC; hMC = rb; }
      }

      // Normalize
      if (hData) {
        Double_t nevts = GetNORMEVENTS() ? Nevents(cfg.GetDataFile().Data(), cfg.dataDir.Data(), GetEventObj(), 0) : 1.0;
        hset(*hData, GetTrackPtTitleX(), "1/#it{N}_{evt} d#it{N}/d#it{p}_{T}", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
        hoptset(*hData, nevts, col, 0, 200, 2e-10, 1e0 + 0.05, 1, 1, 2, 20);
        hData->SetStats(0);
      }
      if (hMC) {
        Double_t nevtsMC = GetNORMEVENTS() ? Nevents(cfg.GetMcFile().Data(), cfg.mcDir.Data(), GetEventObj(), 0, cfg.isJJ) : 1.0;
        hset(*hMC, GetTrackPtTitleX(), "1/#it{N}_{evt} d#it{N}/d#it{p}_{T}", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
        hoptset(*hMC, nevtsMC, col, 0, 200, 2e-10, 1e0 + 0.05, 1, 1, 2, 25);
        hMC->SetMarkerSize(0.8);
        hMC->SetStats(0);
      }

      // Draw on main pad
      mainPadCons->cd();
      if (hData) {
        if (firstDraw) {
          hData->GetXaxis()->SetLabelSize(0);
          hData->GetXaxis()->SetTitleSize(0);
          hData->Draw("pe");
          firstDraw = false;
        } else {
          hData->Draw("pesame");
        }
        legCons->AddEntry(hData, Form("Data, %s", cfg.label.Data()), "pe");
      }
      if (hMC) {
        hMC->Draw("pesame");
        if (i == 0) legCons->AddEntry(hMC, "MC (open markers)", "pe");
      }

      dataHists.push_back(hData);
      mcHists.push_back(hMC);
    }

    mainPadCons->cd();
    legCons->Draw();

    // MC/Data ratio per R
    ratioPadCons->cd();
    bool firstRatio = true;
    for (size_t i = 0; i < rConfigs.size(); ++i) {
      if (!mcHists[i] || !dataHists[i]) continue;
      TH1* hRat = (TH1*)mcHists[i]->Clone(Form("TrackPtMCDataRatio_R%02d", (int)(rConfigs[i].R * 10 + 0.5)));
      hRat->SetDirectory(0);
      hRat->Divide(dataHists[i]);
      Color_t col = GetColorForR(rConfigs[i].R);
      hRat->SetLineColor(col);
      hRat->SetMarkerColor(col);
      hRat->SetMarkerStyle(20);
      hRat->SetMarkerSize(0.8);
      hRat->SetStats(0);
      if (firstRatio) {
        hset(*hRat, GetTrackPtTitleX(), "MC / Data",
             1.2, 1.0, 0.07, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 505);
        hRat->GetYaxis()->SetRangeUser(0.0, 3.0);
        hRat->Draw("pe");
        firstRatio = false;
      } else {
        hRat->Draw("pesame");
      }
    }
    if (!firstRatio) {
      TLine* line1 = new TLine(0, 1.0, 200, 1.0);
      line1->SetLineColor(kBlack);
      line1->SetLineStyle(2);
      line1->Draw();
    }

    if (GetDRAWPLOTS() && padCons->C) {
      padCons->C->Print(Form("%s/TrackPt_AllR.pdf", outputDir.Data()));
    }
  }
}

void DrawConstituentPtObservables() {
  if (!GetDrawConstituentPt()) {
    std::cout << "[DrawConstituentPtObservables] DRAW_CONSTITUENT_PT is false, skipping..." << std::endl;
    return;
  }
  
  // Get R configurations and colors
  std::vector<RConfig> rConfigs = GetRConfigs();
  std::vector<Color_t> rColors = GetRColorsForConfigs(rConfigs);
  TString outputDir = GetOutputDir();
  
  std::cout << "[DrawConstituentPtObservables] Starting constituent pT drawing..." << std::endl;
  
  // Create pad
  Filipad2* pad = new Filipad2("ConstituentPt", ++GetNN(), 2, 0.4, 100, 50, 0.7, 1, 1);
  pad->Draw();
  TPad* mainPad = pad->GetPad(1);
  TPad* ratioPad = pad->GetPad(2);
  if (!mainPad || !ratioPad) {
    std::cerr << "[Error] Failed to get pads for constituent pT" << std::endl;
    return;
  }
  
  optFili(*mainPad, 0, 1, 0, 1);  // logx=0, logy=1 (y축 log)
  optFili(*ratioPad, 0, 0, 0, 0);  // logx=0, logy=0 (ratio y축 linear)
  
  // Legend shifted left by 0.5 in x to stay inside pad
  TLegend* leg = new TLegend(0.7, 0.425, 0.95, 0.95, NULL, "brNDC");
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  
  mainPad->cd();
  
  std::vector<TH1*> mcHists, dataHists;
  
  // Draw all R values
  for (size_t i = 0; i < rConfigs.size(); ++i) {
    const auto& config = rConfigs[i];
    Color_t color = rColors[i];
    
    std::cout << "[DrawConstituentPtObservables] Processing R=" << config.R << std::endl;
    
    // Get MC histogram
    TString mcFile = config.GetMcFile();
    Double_t nevtsMC = GetNORMEVENTS() ? Nevents(mcFile.Data(), config.mcDir.Data(), GetEventObj(), 0, config.isJJ) : 1.0;
    
    TH1* hMC = DrawConstituentPt(mcFile.Data(), Form("%s MC", config.label.Data()), nevtsMC, leg, color, config.mcDir.Data());
    
    if (hMC) {
      hMC->SetMarkerStyle(25);  // Empty square for MC
      hMC->SetFillStyle(0);
    }
    mcHists.push_back(hMC);  // always push (even nullptr) to keep aligned with rConfigs

    // Get Data histogram
    TString dataFile = config.GetDataFile();
    Double_t nevtsData = GetNORMEVENTS() ? Nevents(dataFile.Data(), config.dataDir.Data(), GetEventObj(), 0) : 1.0;

    TH1* hData = DrawConstituentPt(dataFile.Data(), Form("%s Data", config.label.Data()), nevtsData, leg, color, config.dataDir.Data());

    if (hData) {
      hData->SetMarkerStyle(20);  // Filled circle for Data
      hData->SetFillStyle(1001);
    }
    dataHists.push_back(hData);  // always push (even nullptr) to keep aligned with rConfigs
  }
  
  mainPad->cd();
  leg->Draw();
  
  // Draw ratio plots (Data/MC for each R)
  ratioPad->cd();
  for (size_t i = 0; i < std::min(mcHists.size(), dataHists.size()); ++i) {
    if (!mcHists[i] || !dataHists[i]) {
      std::cerr << "[Warning] Missing histogram for ratio plot, R=" << rConfigs[i].R << std::endl;
      continue;
    }
    
    // Draw MC/Data = MC / Data
    TH1* ratio = DrawRatio(Form("Ratio_ConstituentPt_R%.1f", rConfigs[i].R), 
                          dataHists[i], mcHists[i],
                          "#it{p}_{T, constituent}^{reco} (GeV/#it{c})", "MC / Data", rColors[i], 0.0, 5.0);
    if (!ratio) {
      std::cerr << "[Error] Failed to create ratio plot for R=" << rConfigs[i].R << std::endl;
    }
  }
  
  if (GetDRAWPLOTS()) {
    if (!pad->C) {
      std::cerr << "[Error] Canvas is null, cannot print!" << std::endl;
    } else {
      TString outputPath = Form("%s/ConstituentPt_RDependent.pdf", outputDir.Data());
      std::cout << "[DrawConstituentPtObservables] Printing to " << outputPath.Data() << std::endl;
      pad->C->Print(outputPath.Data());
    }
  }
  
  std::cout << "[DrawConstituentPtObservables] Completed!" << std::endl;
}

void DrawNtracksObservables() {
  if (!GetDrawNtracks()) {
    std::cout << "[DrawNtracksObservables] DRAW_NTRACKS is false, skipping..." << std::endl;
    return;
  }
  
  // Get R configurations and colors
  std::vector<RConfig> rConfigs = GetRConfigs();
  std::vector<Color_t> rColors = GetRColorsForConfigs(rConfigs);
  TString outputDir = GetOutputDir();
  
  std::cout << "[DrawNtracksObservables] Starting Ntracks distribution drawing..." << std::endl;
  
  // Create pad
  Filipad2* pad = new Filipad2("Ntracks", ++GetNN(), 2, 0.4, 100, 50, 0.7, 1, 1);
  pad->Draw();
  TPad* mainPad = pad->GetPad(1);
  TPad* ratioPad = pad->GetPad(2);
  if (!mainPad || !ratioPad) {
    std::cerr << "[Error] Failed to get pads for Ntracks" << std::endl;
    return;
  }
  
  optFili(*mainPad, 0, 1, 0, 1);  // logx=0, logy=1 (y축 log)
  optFili(*ratioPad, 0, 0, 0, 0);  // logx=0, logy=0 (ratio y축 linear)
  
  // Legend shifted left by 0.5 in x to stay inside pad
  TLegend* leg = new TLegend(0.7, 0.425, 0.95, 0.95, NULL, "brNDC");
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  
  mainPad->cd();
  
  std::vector<TH1*> mcHists, dataHists;
  
  // Draw all R values
  for (size_t i = 0; i < rConfigs.size(); ++i) {
    const auto& config = rConfigs[i];
    Color_t color = rColors[i];
    
    std::cout << "[DrawNtracksObservables] Processing R=" << config.R << std::endl;
    
    // Get MC histogram
    TString mcFile = config.GetMcFile();
    Double_t nevtsMC = GetNORMEVENTS() ? Nevents(mcFile.Data(), config.mcDir.Data(), GetEventObj(), 0, config.isJJ) : 1.0;
    
    TH1* hMC = DrawNtracksDistribution(mcFile.Data(), Form("%s MC", config.label.Data()), nevtsMC, leg, color, config.mcDir.Data());
    
    if (hMC) {
      hMC->SetMarkerStyle(25);  // Empty square for MC
      hMC->SetFillStyle(0);
    }
    mcHists.push_back(hMC);  // always push (even nullptr) to keep aligned with rConfigs

    // Get Data histogram
    TString dataFile = config.GetDataFile();
    Double_t nevtsData = GetNORMEVENTS() ? Nevents(dataFile.Data(), config.dataDir.Data(), GetEventObj(), 0) : 1.0;

    TH1* hData = DrawNtracksDistribution(dataFile.Data(), Form("%s Data", config.label.Data()), nevtsData, leg, color, config.dataDir.Data());

    if (hData) {
      hData->SetMarkerStyle(20);  // Filled circle for Data
      hData->SetFillStyle(1001);
    }
    dataHists.push_back(hData);  // always push (even nullptr) to keep aligned with rConfigs
  }
  
  mainPad->cd();
  leg->Draw();
  
  // Draw ratio plots (Data/MC for each R)
  ratioPad->cd();
  for (size_t i = 0; i < std::min(mcHists.size(), dataHists.size()); ++i) {
    if (!mcHists[i] || !dataHists[i]) {
      std::cerr << "[Warning] Missing histogram for ratio plot, R=" << rConfigs[i].R << std::endl;
      continue;
    }
    
    // Draw MC/Data = MC / Data
    TH1* ratio = DrawRatio(Form("Ratio_Ntracks_R%.1f", rConfigs[i].R), 
                          dataHists[i], mcHists[i],
                          "N_{tracks}", "MC / Data", rColors[i], 0.0, 3.5);
    if (!ratio) {
      std::cerr << "[Error] Failed to create ratio plot for R=" << rConfigs[i].R << std::endl;
    }
  }
  
  if (GetDRAWPLOTS()) {
    if (!pad->C) {
      std::cerr << "[Error] Canvas is null, cannot print!" << std::endl;
    } else {
      TString outputPath = Form("%s/Ntracks_RDependent.pdf", outputDir.Data());
      std::cout << "[DrawNtracksObservables] Printing to " << outputPath.Data() << std::endl;
      pad->C->Print(outputPath.Data());
    }
  }
  
  std::cout << "[DrawNtracksObservables] Completed!" << std::endl;
}

void DrawJetObservables() {
  if (!GetDrawJet()) return;
  
  // Get R configurations and colors
  std::vector<RConfig> rConfigs = GetRConfigs();
  std::vector<Color_t> rColors = GetRColorsForConfigs(rConfigs);
  TString outputDir = GetOutputDir();
  
  std::vector<TString> observables = {"pt", "eta", "phi"};
  
  for (const auto& obs : observables) {
    Filipad2* pad = new Filipad2(Form("Jet%s", obs.Data()), ++GetNN(), 2, 0.4, 100, 50, 0.7, 1, 1);
    pad->Draw();
    TPad* mainPad = pad->GetPad(1);
    TPad* ratioPad = pad->GetPad(2);
    // pt는 y축 log scale, eta/phi는 linear scale
    int logy = (obs == "pt") ? 1 : 0;
    optFili(*mainPad, 0, logy, 0, 1);  // logx=0 (x축 linear), logy는 obs에 따라 결정
    optFili(*ratioPad, 0, 0, 0, 0);  // logx=0, logy=0 (ratio y축도 linear)
    
    TLegend* leg = new TLegend(0.5, 0.425, 0.9, 0.95, NULL, "brNDC");  // 높이 1.5배 (0.6->0.425)
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    
    mainPad->cd();
    
    std::vector<TH1*> mcHists, dataHists;
    
    for (size_t i = 0; i < rConfigs.size(); ++i) {
      const auto& config = rConfigs[i];
      Color_t color = rColors[i];
      
      Double_t nevtsMC = GetNORMEVENTS() ? Nevents(config.GetMcFile().Data(), config.mcDir.Data(), GetEventObj(), 0, config.isJJ) : 1.0;
      Double_t nevtsData = GetNORMEVENTS() ? Nevents(config.GetDataFile().Data(), config.dataDir.Data(), GetEventObj(), 0) : 1.0;
      
      TH1* hMC = nullptr;
      TH1* hData = nullptr;
      
      if (obs == "pt") {
        hMC = DrawJetPt(config.GetMcFile().Data(), Form("%s MC", config.label.Data()), GetJetPtObj(), nevtsMC, leg, color, 0, config.mcDir.Data(), config.R);
        hData = DrawJetPt(config.GetDataFile().Data(), Form("%s Data", config.label.Data()), GetJetPtObj(), nevtsData, leg, color, 0, config.dataDir.Data(), config.R);
      } else if (obs == "eta") {
        hMC = DrawJetEta(config.GetMcFile().Data(), Form("%s MC", config.label.Data()), GetJetEtaObj(), nevtsMC, leg, color, config.mcDir.Data(), config.R);
        hData = DrawJetEta(config.GetDataFile().Data(), Form("%s Data", config.label.Data()), GetJetEtaObj(), nevtsData, leg, color, config.dataDir.Data(), config.R);
      } else if (obs == "phi") {
        hMC = DrawJetPhi(config.GetMcFile().Data(), Form("%s MC", config.label.Data()), GetJetPhiObj(), nevtsMC, leg, color, config.mcDir.Data(), config.R);
        hData = DrawJetPhi(config.GetDataFile().Data(), Form("%s Data", config.label.Data()), GetJetPhiObj(), nevtsData, leg, color, config.dataDir.Data(), config.R);
      }
      
      if (hMC) {
        hMC->SetMarkerStyle(25);  // Empty square for MC
        hMC->SetFillStyle(0);
      }
      mcHists.push_back(hMC);  // always push (even nullptr) to keep aligned with rConfigs
      if (hData) {
        hData->SetMarkerStyle(20);  // Filled circle for Data
        hData->SetFillStyle(1001);
      }
      dataHists.push_back(hData);  // always push (even nullptr) to keep aligned with rConfigs
    }
    
    // Auto-range y-axis for eta/phi to show all R values
    mainPad->cd();
    if (obs == "eta" || obs == "phi") {
      Double_t yMaxAll = 0;
      for (auto* h : mcHists)   if (h) yMaxAll = TMath::Max(yMaxAll, h->GetMaximum());
      for (auto* h : dataHists) if (h) yMaxAll = TMath::Max(yMaxAll, h->GetMaximum());
      if (yMaxAll > 0) {
        // Set generous y-range on the first valid histogram (controls axes)
        for (auto* h : mcHists) {
          if (h) { h->GetYaxis()->SetRangeUser(0, yMaxAll * 1.4); break; }
        }
        for (auto* h : dataHists) {
          if (h) { h->GetYaxis()->SetRangeUser(0, yMaxAll * 1.4); break; }
        }
      }
    }
    leg->Draw();

    // Add "UE subtracted" label for UE mode (jet plots only, not track plots)
    AddUESubtractedLabel(mainPad);
    
    // Ratio plots
    ratioPad->cd();
    TString xAxisTitle;
    if (obs == "pt") {
      xAxisTitle = GetJetPtTitleX();
    } else if (obs == "eta") {
      xAxisTitle = GetJetEtaTitleX();
    } else if (obs == "phi") {
      xAxisTitle = GetJetPhiTitleX();
    }
    for (size_t i = 0; i < std::min(mcHists.size(), dataHists.size()); ++i) {
      if (!mcHists[i] || !dataHists[i]) continue;
      // Draw MC/Data = MC / Data
      DrawRatio(Form("Ratio_Jet%s_R%.1f", obs.Data(), rConfigs[i].R),
               dataHists[i], mcHists[i],
               xAxisTitle, "MC / Data", rColors[i], 0.5, 2.5);
    }
    
    if (GetDRAWPLOTS()) {
      pad->C->Print(Form("%s/Jet%s_RDependent.pdf", outputDir.Data(), obs.Data()));
    }
  }
}

void DrawJetPtPart() {
  if (!GetDrawJetPtPart()) return;

  std::cout << "[DrawJetPtPart] Starting..." << std::endl;
  
  // Get R configurations and colors
  std::vector<RConfig> rConfigs = GetRConfigs();
  std::vector<Color_t> rColors = GetRColorsForConfigs(rConfigs);
  TString outputDir = GetOutputDir();
  
  // Create canvas with 4 pads: main (0.4) + 3 ratio pads (each 0.2)
  TString padName = Form("JetPtPart_RDependent_%d", ++GetNN());
  TCanvas* can = new TCanvas(padName.Data(), "", 800, 1000);
  can->SetTopMargin(0.);
  can->SetBottomMargin(0.);
  can->SetLeftMargin(0.);
  can->SetRightMargin(0.);
  can->Draw();
  
  // Pad sizes: main=0.4, ratio1=0.2, ratio2=0.2, ratio3=0.2
  // From bottom to top: ratio3 (0-0.2), ratio2 (0.2-0.4), ratio1 (0.4-0.6), main (0.6-1.0)
  float mainRatio = 0.4f;
  float ratioSize = 0.2f;
  
  TPad* mainPad = new TPad("mainPad", "mainPad", 0.0, 0.6, 1.0, 1.0, 0);
  mainPad->SetTopMargin(0.02 / mainRatio);
  mainPad->SetBottomMargin(0.0015);
  mainPad->SetLeftMargin(0.15);
  mainPad->SetRightMargin(0.03);
  mainPad->Draw();
  
  TPad* ratioPad1 = new TPad("ratioPad1", "ratioPad1", 0.0, 0.43, 1.0, 0.6, 0);
  ratioPad1->SetTopMargin(0.0015);
  ratioPad1->SetBottomMargin(0.0015);
  ratioPad1->SetLeftMargin(0.15);
  ratioPad1->SetRightMargin(0.03);
  ratioPad1->Draw();
  
  TPad* ratioPad2 = new TPad("ratioPad2", "ratioPad2", 0.0, 0.26, 1.0, 0.43, 0);
  ratioPad2->SetTopMargin(0.0015);
  ratioPad2->SetBottomMargin(0.0015);
  ratioPad2->SetLeftMargin(0.15);
  ratioPad2->SetRightMargin(0.03);
  ratioPad2->Draw();
  
  TPad* ratioPad3 = new TPad("ratioPad3", "ratioPad3", 0.0, 0.0, 1.0, 0.26, 0);
  ratioPad3->SetTopMargin(0.0015);
  ratioPad3->SetBottomMargin(0.08 / ratioSize);
  ratioPad3->SetLeftMargin(0.15);
  ratioPad3->SetRightMargin(0.03);
  ratioPad3->Draw();
  
  optFili(*mainPad, 0, 1, 0, 1);  // logx=0 (x축 linear), logy=1 (y축 log)
  optFili(*ratioPad1, 0, 0, 0, 0);  // logx=0, logy=0 (ratio y축도 linear)
  optFili(*ratioPad2, 0, 0, 0, 0);  // logx=0, logy=0 (ratio y축도 linear)
  optFili(*ratioPad3, 0, 0, 0, 0);  // logx=0, logy=0 (ratio y축도 linear)
  
  TLegend* leg = new TLegend(0.65, 0.3, 0.9, 0.9, NULL, "brNDC");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  
  // Add header to legend
  leg->AddEntry((TObject*)0, "MC particle-level jets", "");
  
  mainPad->cd();
  
  std::vector<TH1*> mcpHists;
  bool firstHist = true;
  
  // Draw all R values
  for (size_t i = 0; i < rConfigs.size(); ++i) {
    const auto& config = rConfigs[i];
    Double_t nevtsMCP = GetNORMEVENTS() ? Nevents(config.GetMcFile().Data(), config.mcDir.Data(), GetEventObj(), 1, config.isJJ) : 1.0;
    
    // Create unique histogram name for each R value
    TString uniqueName = Form("JetPtMCP_R%.1f", config.R);
    TString legendLabel = Form("R=%.1f", config.R);
    
    TH1* hMCP = DrawJetPtMCP(config.GetMcFile().Data(), uniqueName.Data(), GetJetPtMCPObj(), nevtsMCP, leg, rColors[i], i, config.mcDir.Data(), legendLabel.Data());
    if (hMCP) {
      // Redraw with correct option (first one uses "pe", others use "pesame")
      mainPad->cd();
      if (firstHist) {
        hMCP->Draw("pe");
        firstHist = false;
      } else {
        hMCP->Draw("pesame");
      }
    }
    mcpHists.push_back(hMCP);  // always push (even nullptr) to keep aligned with rConfigs
  }

  // Overlay standalone model predictions
  auto models = GetStandaloneModels();
  for (auto& model : models) {
    bool legendAdded = false;
    for (size_t i = 0; i < rConfigs.size(); ++i) {
      const auto& config = rConfigs[i];
      TH1* hModel = LoadStandaloneModelNormYield(model, config.R);
      if (!hModel) continue;
      mainPad->cd();
      hModel->SetLineColor(rColors[i]);
      hModel->SetMarkerColor(rColors[i]);
      hModel->SetLineStyle(model.lineStyle);
      hModel->SetLineWidth(model.lineWidth);
      hModel->SetMarkerStyle(0);
      hModel->GetXaxis()->SetRangeUser(0, GetPlotPtMax());
      hModel->Draw("l same");
      if (!legendAdded) {
        leg->AddEntry(hModel, model.legLabel.Data(), "l");
        legendAdded = true;
      }
    }
  }

  leg->Draw();

  // Add "UE subtracted" label for UE mode (jet plots only)
  AddUESubtractedLabel(mainPad);
  
  // Find indices for R=0.2, 0.4, 0.6
  int idx02 = -1, idx04 = -1, idx06 = -1;
  for (size_t i = 0; i < rConfigs.size(); ++i) {
    if (TMath::Abs(rConfigs[i].R - 0.2) < 0.01) idx02 = i;
    if (TMath::Abs(rConfigs[i].R - 0.4) < 0.01) idx04 = i;
    if (TMath::Abs(rConfigs[i].R - 0.6) < 0.01) idx06 = i;
  }
  
  // Ratio plot 1: R=0.2/R=0.4 (top ratio pad)
  if (idx02 >= 0 && idx04 >= 0 && idx02 < mcpHists.size() && idx04 < mcpHists.size() && mcpHists[idx02] && mcpHists[idx04]) {
    ratioPad1->cd();
    TH1* ratio1 = DrawRatio(Form("Ratio_MCP_R02_R04"), mcpHists[idx04], mcpHists[idx02],
                            GetJetPtGenTitleX(), "R=0.2 / R=0.4", rColors[idx02], 0.0, 3.0);
    if (ratio1) {
      // Increase title sizes for ratio plots
      ratio1->GetXaxis()->SetTitleSize(0.14);
      ratio1->GetYaxis()->SetTitleSize(0.14);
      ratio1->GetXaxis()->SetLabelSize(0.10);
      ratio1->GetYaxis()->SetLabelSize(0.10);
      ratio1->GetYaxis()->SetTitleOffset(0.4);
    }
  }
  
  // Ratio plot 2: R=0.4/R=0.6 (middle ratio pad)
  if (idx04 >= 0 && idx06 >= 0 && idx04 < mcpHists.size() && idx06 < mcpHists.size() && mcpHists[idx04] && mcpHists[idx06]) {
    ratioPad2->cd();
    TH1* ratio2 = DrawRatio(Form("Ratio_MCP_R04_R06"), mcpHists[idx06], mcpHists[idx04],
                            GetJetPtGenTitleX(), "R=0.4 / R=0.6", rColors[idx04], 0.0, 3.0);
    if (ratio2) {
      // Increase title sizes for ratio plots
      ratio2->GetXaxis()->SetTitleSize(0.14);
      ratio2->GetYaxis()->SetTitleSize(0.14);
      ratio2->GetXaxis()->SetLabelSize(0.10);
      ratio2->GetYaxis()->SetLabelSize(0.10);
      ratio2->GetYaxis()->SetTitleOffset(0.4);
    }
  }
  
  // Ratio plot 3: R=0.2/R=0.6 (bottom ratio pad)
  if (idx02 >= 0 && idx06 >= 0 && idx02 < mcpHists.size() && idx06 < mcpHists.size() && mcpHists[idx02] && mcpHists[idx06]) {
    ratioPad3->cd();
    TH1* ratio3 = DrawRatio(Form("Ratio_MCP_R02_R06"), mcpHists[idx06], mcpHists[idx02],
                            GetJetPtGenTitleX(), "R=0.2 / R=0.6", rColors[idx06], 0.0, 3.0);
    if (ratio3) {
      // Increase title sizes for ratio plots
      ratio3->GetXaxis()->SetTitleSize(0.09);
      ratio3->GetYaxis()->SetTitleSize(0.09);
      ratio3->GetXaxis()->SetLabelSize(0.1);
      ratio3->GetYaxis()->SetLabelSize(0.1);
      ratio3->GetYaxis()->SetTitleOffset(0.6);
    }
  }
  
  if (GetDRAWPLOTS()) {
    can->Print(Form("%s/JetPtPart_RDependent.pdf", outputDir.Data()));
  }
}

void DrawJetAreaPerR() {
  if (!GetDrawJetArea()) return;
  
  // Get R configurations and colors
  std::vector<RConfig> rConfigs = GetRConfigs();
  std::vector<Color_t> rColors = GetRColorsForConfigs(rConfigs);
  TString outputDir = GetOutputDir();
  
  for (size_t i = 0; i < rConfigs.size(); ++i) {
    const auto& config = rConfigs[i];
    
    TCanvas* can = new TCanvas(Form("JetArea_R%.1f", config.R), Form("Jet Area R=%.1f", config.R), 800, 600);
    can->SetLogz(1);
    setpad(can, 0.02, 0.12, 0.12, 0.01);
    can->Draw();
    
    TLegend* leg = new TLegend(0.6, 0.6, 0.9, 0.9, NULL, "brNDC");  // 높이 1.5배 (0.7->0.6)
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    
    Double_t nevtsMC = GetNORMEVENTS() ? Nevents(config.GetMcFile().Data(), config.mcDir.Data(), GetEventObj(), 0, config.isJJ) : 1.0;
    
    TH2* hArea = DrawJetArea(config.GetMcFile().Data(), config.label.Data(), GetJetAreaObj(), nevtsMC, leg, rColors[i], config.mcDir.Data());
    
    leg->Draw();
    
    // Add "UE subtracted" label for UE mode (jet plots only)
    // Here we draw directly on the canvas (no Filipad2 mainPad in this function).
    AddUESubtractedLabel(can);
    
    if (GetDRAWPLOTS()) {
      can->Print(Form("%s/JetArea_R%.1f.pdf", outputDir.Data(), config.R));
    }
  }
}

void DrawResponseMatrixPerR() {
  if (!GetDrawResponseMatrix()) return;
  
  // Get R configurations
  std::vector<RConfig> rConfigs = GetRConfigs();
  TString outputDir = GetOutputDir();
  
  for (size_t i = 0; i < rConfigs.size(); ++i) {
    const auto& config = rConfigs[i];
    
    TCanvas* can = new TCanvas(Form("ResponseMatrix_R%.1f", config.R), Form("Response Matrix R=%.1f", config.R), 800, 600);
    can->SetLogz(1);
    setpad(can, 0.02, 0.15, 0.15, 0.15);  // Reduced left margin, increased right margin for color scale
    can->Draw();
    
    TFile* mcFile = TFile::Open(config.GetMcFile().Data(), "READ");
    if (!mcFile || mcFile->IsZombie()) {
      std::cerr << "[Error] Cannot open MC file: " << config.GetMcFile().Data() << " for R=" << config.R << std::endl;
      continue;
    }
    
    // Get response matrix (2D correlation plot) - use mode-dependent name
    const char* respMatrixName = (GetCurrentMode() == kUE) ? 
      "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted_mcdetaconstraint" : 
      "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint";
    TH2* h2Correlate = (TH2*)mcFile->Get(Form("%s/%s", config.mcDir.Data(), respMatrixName));
    if (!h2Correlate) {
      std::cerr << "[Error] Response matrix not found for R=" << config.R << " at path: " << Form("%s/%s", config.mcDir.Data(), respMatrixName) << std::endl;
      mcFile->Close();
      continue;
    }
    
    // Clone and detach from file
    TH2* h2CorrelateClone = (TH2*)h2Correlate->Clone(Form("h2Correlate_R%.1f", config.R));
    h2CorrelateClone->SetDirectory(0);
    h2CorrelateClone->SetStats(0);  // Disable stat box
    h2CorrelateClone->SetTitle("");  // Remove histogram title

    mcFile->Close();

    // Normalize each row (fixed true pT bin) to probability: P(reco pT | true pT)
    // X axis = reco pT, Y axis = true pT
    for (Int_t jy = 1; jy <= h2CorrelateClone->GetNbinsY(); ++jy) {
      Double_t rowSum = 0;
      for (Int_t ix = 1; ix <= h2CorrelateClone->GetNbinsX(); ++ix) {
        rowSum += h2CorrelateClone->GetBinContent(ix, jy);
      }
      if (rowSum > 0) {
        for (Int_t ix = 1; ix <= h2CorrelateClone->GetNbinsX(); ++ix) {
          h2CorrelateClone->SetBinContent(ix, jy, h2CorrelateClone->GetBinContent(ix, jy) / rowSum);
          h2CorrelateClone->SetBinError(ix, jy, h2CorrelateClone->GetBinError(ix, jy) / rowSum);
        }
      }
    }

    // Set up histogram
    can->cd();
    can->SetLogz(1);  // Log z for probability distribution
    hset(*h2CorrelateClone, "#it{p}_{T, jet}^{reco} (GeV/#it{c})", "#it{p}_{T, jet}^{true} (GeV/#it{c})",
         1.2, 1.3, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    h2CorrelateClone->GetXaxis()->SetRangeUser(0, 200);
    h2CorrelateClone->GetYaxis()->SetRangeUser(0, 200);
    h2CorrelateClone->GetZaxis()->SetRangeUser(1e-10, 1);
    h2CorrelateClone->Draw("colz");
    
    // Add "UE subtracted" label for UE mode (jet plots only)
    AddUESubtractedLabel(can);
    
    if (GetDRAWPLOTS()) {
      can->Print(Form("%s/ResponseMatrix_R%.1f.pdf", outputDir.Data(), config.R));
    }
  }
}

void DrawPurityEfficiency() {
  if (!GetDrawPurityEfficiency()) return;
  
  // Get R configurations and colors
  std::vector<RConfig> rConfigs = GetRConfigs();
  std::vector<Color_t> rColors = GetRColorsForConfigs(rConfigs);
  
  // Create static canvases (similar to DrawJetMatching implementation)
  static TCanvas* sJREcanvas = nullptr;
  static TLegend* sJREcanLegend = nullptr;
  static bool sJREfirst = true;
  
  static TCanvas* sJRPcanvas = nullptr;
  static TLegend* sJRPcanLegend = nullptr;
  static bool sJRPfirst = true;
  
  // Reset and clear if canvases already exist (to avoid duplicate drawing)
  if (sJREcanvas && !sJREfirst) {
    sJREcanvas->Clear();
    sJREfirst = true;  // Reset to allow redrawing
    if (sJREcanLegend) {
      sJREcanLegend->Clear();
    }
  }
  
  if (sJRPcanvas && !sJRPfirst) {
    sJRPcanvas->Clear();
    sJRPfirst = true;  // Reset to allow redrawing
    if (sJRPcanLegend) {
      sJRPcanLegend->Clear();
    }
  }
  
  if (!sJREcanvas) {
    sJREcanvas = new TCanvas("JREcanvas_AllR", "", 800, 600);  // Remove title
    setpad(sJREcanvas, 0.02, 0.15, 0.15, 0.05);
    sJREcanvas->Draw();
    // Legend: height increased by 20%, positioned at bottom right corner
    // Current height: 0.233, new height: 0.233 * 1.2 = 0.28
    // y_min = 0.1 (10% from pad bottom), y_max = 0.1 + 0.28 = 0.38
    sJREcanLegend = new TLegend(0.75, 0.20, 0.95, 0.45, NULL, "brNDC");  // Bottom right, 20% taller, 10% from bottom
    sJREcanLegend->SetBorderSize(0);
    sJREcanLegend->SetTextSize(0.04);
    sJREcanLegend->SetFillColorAlpha(0, 0);
  }
  
  if (!sJRPcanvas) {
    sJRPcanvas = new TCanvas("JRPcanvas_AllR", "", 800, 600);  // Remove title
    setpad(sJRPcanvas, 0.02, 0.15, 0.15, 0.05);
    sJRPcanvas->Draw();
    // Legend: height increased by 20%, positioned at bottom right corner
    sJRPcanLegend = new TLegend(0.75, 0.20, 0.95, 0.45, NULL, "brNDC");  // Bottom right, 20% taller, 10% from bottom
    sJRPcanLegend->SetBorderSize(0);
    sJRPcanLegend->SetTextSize(0.04);
    sJRPcanLegend->SetFillColorAlpha(0, 0);
  }
  
  for (size_t i = 0; i < rConfigs.size(); ++i) {
    const auto& config = rConfigs[i];
    
    Double_t nevtsData = GetNORMEVENTS() ? Nevents(config.GetDataFile().Data(), config.dataDir.Data(), GetEventObj(), 0) : 1.0;
    Double_t nevtsDataUnTrig = GetNORMEVENTS() ? Nevents(config.GetDataFile().Data(), config.dataDir.Data(), GetEventObj(), 1) : 1.0;
    Double_t nevtsMCD = GetNORMEVENTS() ? Nevents(config.GetMcFile().Data(), config.mcDir.Data(), GetEventObj(), 0, config.isJJ) : 1.0;
    Double_t nevtsMCP = GetNORMEVENTS() ? Nevents(config.GetMcFile().Data(), config.mcDir.Data(), GetEventObj(), 1, config.isJJ) : 1.0;
    
    // Open MC file and calculate efficiency/purity directly (following DrawJetsMCfJetQA.h approach)
    TFile* mcFile = TFile::Open(config.GetMcFile().Data(), "READ");
    if (!mcFile || mcFile->IsZombie()) {
      std::cerr << "[Error] Cannot open MC file: " << config.GetMcFile().Data() << " for R=" << config.R << std::endl;
      continue;
    }
    
    // Get histograms
    TH1* JetMCPPt = (TH1*)mcFile->Get(Form("%s/%s", config.mcDir.Data(), GetJetPtMCPObj()));
    if (!JetMCPPt) {
      std::cerr << "[Error] Cannot find " << GetJetPtMCPObj() << " for R=" << config.R << std::endl;
      mcFile->Close();
      continue;
    }
    
    // Use mode-dependent histogram: h_jet_pt_rhoareasubtracted for UE, h_jet_pt for NonUE
    const char* jetPtMCDName = (GetCurrentMode() == kUE) ? "h_jet_pt_rhoareasubtracted" : "h_jet_pt";
    TH1* JetMCDPt = (TH1*)mcFile->Get(Form("%s/%s", config.mcDir.Data(), jetPtMCDName));
    if (!JetMCDPt) {
      std::cerr << "[Error] Cannot find " << jetPtMCDName << " for R=" << config.R << std::endl;
      mcFile->Close();
      continue;
    }
    
    // Use mode-dependent response matrix name
    const char* respMatrixName = (GetCurrentMode() == kUE) ? 
      "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted_mcdetaconstraint" : 
      "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint";
    TH2* HCorrelate2D = (TH2*)mcFile->Get(Form("%s/%s", config.mcDir.Data(), respMatrixName));
    if (!HCorrelate2D) {
      std::cerr << "[Error] Cannot find 2D correlation histogram for R=" << config.R << " at path: " << Form("%s/%s", config.mcDir.Data(), respMatrixName) << std::endl;
      mcFile->Close();
      continue;
    }
    
    // Detect binning from response matrix (authoritative for consistent rebinning)
    PtBinning recoBin, truthBin;
    DetectResponseBinning(HCorrelate2D, recoBin, truthBin);

    // Rebin 1D histograms with detected binning (safe: handles narrower source range)
    if (GetREBINON()) {
      JetMCPPt = RebinToTarget(JetMCPPt, truthBin, Form("JetMCPPt_R%.1f", config.R));
      JetMCDPt = RebinToTarget(JetMCDPt, recoBin, Form("JetMCDPt_R%.1f", config.R));
    }

    // Create rebinned 2D correlation if needed
    TH2F* h2HCorrelate = nullptr;
    if (GetREBINON()) {
      h2HCorrelate = new TH2F(Form("hcorrelate_R%.1f", config.R), "",
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
      h2HCorrelate = (TH2F*)HCorrelate2D->Clone(Form("hcorrelate_R%.1f", config.R));
      h2HCorrelate->SetTitle("");
    }

    // Project to get matched distributions
    TH1* MCDMatchedpt = (TH1*)h2HCorrelate->ProjectionX(Form("hMCDMatched_R%.1f", config.R), 1, h2HCorrelate->GetNbinsY(), "e");
    TH1* MCPMatchedpt = (TH1*)h2HCorrelate->ProjectionY(Form("hMCPMatched_R%.1f", config.R), 1, h2HCorrelate->GetNbinsX(), "e");

    MCDMatchedpt->SetTitle("");
    MCPMatchedpt->SetTitle("");

    // Rebin projections to match 1D histograms (safety: should already match from h2HCorrelate)
    if (GetREBINON()) {
      TH1* tempMCPMatchedpt = RebinToTarget(MCPMatchedpt, truthBin, Form("tempMCPMatchedpt_R%.1f", config.R));
      delete MCPMatchedpt;
      MCPMatchedpt = tempMCPMatchedpt;
    }

    if (GetREBINON()) {
      TH1* tempMCDMatchedpt = RebinToTarget(MCDMatchedpt, recoBin, Form("tempMCDMatchedpt_R%.1f", config.R));
      delete MCDMatchedpt;
      MCDMatchedpt = tempMCDMatchedpt;
    }
    
    MCDMatchedpt->SetDirectory(0);
    MCPMatchedpt->SetDirectory(0);
    JetMCPPt->SetDirectory(0);
    JetMCDPt->SetDirectory(0);
    
    mcFile->Close();
    
    // Calculate JRE: mcpmatchedpt / JetMCPPt
    // Ensure both histograms have the same binning
    auto JREp = (TH1*)JetMCPPt->Clone(Form("JREp_R%.1f", config.R));
    JREp->SetDirectory(0);
    JREp->GetXaxis()->SetRangeUser(5, GetPlotPtMax());
    
    auto mcpmatchedpt = (TH1*)MCPMatchedpt->Clone(Form("mcpmatchedpt_R%.1f", config.R));
    mcpmatchedpt->SetDirectory(0);
    mcpmatchedpt->GetXaxis()->SetRangeUser(0, GetPlotPtMax());
    
    // Verify binning matches before division
    if (mcpmatchedpt->GetNbinsX() != JREp->GetNbinsX()) {
      std::cerr << "[Warning] Binning mismatch for JRE calculation: mcpmatchedpt has " 
                << mcpmatchedpt->GetNbinsX() << " bins, JREp has " << JREp->GetNbinsX() << " bins for R=" << config.R << std::endl;
    }
    
    TH1* jre = (TH1*)mcpmatchedpt->Clone(Form("jre_R%.1f", config.R));
    jre->Divide(mcpmatchedpt, JREp, 1., 1., "B");
          jre->SetDirectory(0);
    jre->SetStats(0);
    jre->SetTitle("");  // Remove histogram title
    
    // Draw JRE on static canvas
          sJREcanvas->cd();
    gPad->SetGrid(0);
    hoptset(*jre, 0, rColors[i], 5, GetPlotPtMax(), 0.0, 1.1, 1, 1, 2, 20);  // Y range: 0 to 1.1
          hset(*jre, GetJRETitleX(), "#it{#varepsilon}_{reco}^{jet}", 1.3, 1.0, 0.05, 0.07, 0.01, 0.01, 0.05, 0.05, 510, 510);
          jre->SetMarkerStyle(20);
          jre->SetMarkerColor(rColors[i]);
          jre->SetLineColor(rColors[i]);
    if (sJREfirst) {
      jre->Draw("pe");
          sJREfirst = false;
    } else {
      jre->Draw("pesame");
    }
    // Only add entry once per R value
    bool jreEntryExists = false;
    TList* jreEntries = sJREcanLegend->GetListOfPrimitives();
    if (jreEntries) {
      TIter next(jreEntries);
      TObject* obj;
      while ((obj = next())) {
        TLegendEntry* entry = (TLegendEntry*)obj;
        if (entry && TString(entry->GetLabel()) == config.label) {
          jreEntryExists = true;
          break;
        }
      }
    }
    if (!jreEntryExists) {
      sJREcanLegend->AddEntry(jre, config.label.Data(), "pe");
    }
    
    // Calculate JRP: mcdmatchedpt / JetMCDPt
    auto JRPp = (TH1F*)JetMCDPt->Clone(Form("JRPp_R%.1f", config.R));
    JRPp->SetDirectory(0);
    
    auto mcdmatchedpt = (TH1F*)MCDMatchedpt->Clone(Form("mcdmatchedpt_R%.1f", config.R));
    mcdmatchedpt->SetDirectory(0);
    
    // Verify binning matches before division
    if (mcdmatchedpt->GetNbinsX() != JetMCDPt->GetNbinsX()) {
      std::cerr << "[Warning] Binning mismatch for JRP calculation: mcdmatchedpt has " 
                << mcdmatchedpt->GetNbinsX() << " bins, JetMCDPt has " << JetMCDPt->GetNbinsX() << " bins for R=" << config.R << std::endl;
    }
    
    TH1* jrp = (TH1F*)mcdmatchedpt->Clone(Form("jrp_R%.1f", config.R));
    jrp->Divide(mcdmatchedpt, JetMCDPt, 1., 1., "B");
          jrp->SetDirectory(0);
    jrp->SetStats(0);
    jrp->SetTitle("");  // Remove histogram title
    
    // Draw JRP on static canvas
          sJRPcanvas->cd();
    gPad->SetGrid(0);
    hoptset(*jrp, 0, rColors[i], 5, GetPlotPtMax(), 0.0, 1.1, 1, 1, 2, 20);  // Y range: 0 to 1.1
          hset(*jrp, GetJRPTitleX(), "Jet purity", 1.3, 1.0, 0.05, 0.07, 0.01, 0.01, 0.05, 0.05, 510, 510);
    jrp->SetMarkerStyle(20);  // Same marker style as efficiency (filled circle)
          jrp->SetMarkerColor(rColors[i]);
          jrp->SetLineColor(rColors[i]);
    if (sJRPfirst) {
      jrp->Draw("pe");
          sJRPfirst = false;
    } else {
      jrp->Draw("pesame");
    }
    // Only add entry once per R value
    bool jrpEntryExists = false;
    TList* jrpEntries = sJRPcanLegend->GetListOfPrimitives();
    if (jrpEntries) {
      TIter next(jrpEntries);
      TObject* obj;
      while ((obj = next())) {
        TLegendEntry* entry = (TLegendEntry*)obj;
        if (entry && TString(entry->GetLabel()) == config.label) {
          jrpEntryExists = true;
          break;
        }
      }
    }
    if (!jrpEntryExists) {
      sJRPcanLegend->AddEntry(jrp, config.label.Data(), "pe");
    }
  }
  
  TString outputDir = GetOutputDir();
  
  sJREcanvas->cd();
  sJREcanLegend->Draw();
  // Add "UE subtracted" label for UE mode (jet plots only)
  AddUESubtractedLabel(sJREcanvas, 0.4, 0.6);
  if (GetDRAWPLOTS()) {
    sJREcanvas->Print(Form("%s/JetReconstructionEfficiency_AllR.pdf", outputDir.Data()));
  }
  
  sJRPcanvas->cd();
  sJRPcanLegend->Draw();
  // Add "UE subtracted" label for UE mode (jet plots only)
  AddUESubtractedLabel(sJRPcanvas, 0.4, 0.6);
  if (GetDRAWPLOTS()) {
    sJRPcanvas->Print(Form("%s/JetReconstructionPurity_AllR.pdf", outputDir.Data()));
  }
}

void DrawKinematicEfficiency() {
  if (!GetDrawKinematicEfficiency()) return;
  
  // Get R configurations and colors
  std::vector<RConfig> rConfigs = GetRConfigs();
  std::vector<Color_t> rColors = GetRColorsForConfigs(rConfigs);
  TString outputDir = GetOutputDir();
  
  static TCanvas* sJKEcanvas = nullptr;
  static TLegend* sJKEcanLegend = nullptr;
  static bool sJKEfirst = true;
  
  // Reset and clear if canvas already exists (to avoid duplicate drawing)
  if (sJKEcanvas && !sJKEfirst) {
    sJKEcanvas->Clear();
    sJKEfirst = true;  // Reset to allow redrawing
    if (sJKEcanLegend) {
      sJKEcanLegend->Clear();
    }
  }
  
  if (!sJKEcanvas) {
    sJKEcanvas = new TCanvas("JKEcanvas_AllR", "", 800, 600);  // Remove title
    setpad(sJKEcanvas, 0.02, 0.15, 0.15, 0.05);
    sJKEcanvas->Draw();
    // Legend: height increased by 20%, positioned at bottom right corner
    sJKEcanLegend = new TLegend(0.75, 0.20, 0.95, 0.45, NULL, "brNDC");  // Bottom right, 20% taller
    sJKEcanLegend->SetBorderSize(0);
    sJKEcanLegend->SetTextSize(0.04);
    sJKEcanLegend->SetFillColorAlpha(0, 0);
  }
  
  for (size_t i = 0; i < rConfigs.size(); ++i) {
    const auto& config = rConfigs[i];
    
    TFile* mcFile = TFile::Open(config.GetMcFile().Data(), "READ");
    if (!mcFile || mcFile->IsZombie()) continue;
    
    // Get matching 2D histogram
    const char* respMatrixName = (GetCurrentMode() == kUE) ?
      "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted_mcdetaconstraint" :
      "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint";
    TH2* HCorrelate2D = (TH2*)mcFile->Get(Form("%s/%s", config.mcDir.Data(), respMatrixName));
    if (!HCorrelate2D) {
      std::cerr << "[Error] Matching histogram not found for R=" << config.R << std::endl;
      mcFile->Close();
      continue;
    }
    
    // Detect binning from response matrix
    PtBinning recoBin, truthBin;
    DetectResponseBinning(HCorrelate2D, recoBin, truthBin);

    // Create rebinned 2D correlation if needed
    TH2F* h2HCorrelate = nullptr;
    if (GetREBINON()) {
      h2HCorrelate = new TH2F(Form("hcorrelate_JKE_R%.1f", config.R), "",
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
      h2HCorrelate = (TH2F*)HCorrelate2D->Clone(Form("hcorrelate_JKE_R%.1f", config.R));
      h2HCorrelate->SetTitle("");
    }

    // Get total MCP jet distribution (denominator)
    TH1* rec_total_window = (TH1*)h2HCorrelate->ProjectionY(Form("rec_total_R%.1f", config.R), 0, -1, "e");
    rec_total_window->SetTitle("");

    // Project y-axis (true) for x-axis (reco) range: GetPlotPtMin() to GetPlotPtMax()
    Double_t jetPtMin = GetPlotPtMin();
    Double_t jetPtMax = GetPlotPtMax();
    int binLow = h2HCorrelate->GetXaxis()->FindBin(jetPtMin + 1e-5);
    int binHigh = h2HCorrelate->GetXaxis()->FindBin(jetPtMax - 1e-5);
    TH1* rec_selected_window = h2HCorrelate->ProjectionY(Form("rec_selected_R%.1f", config.R), binLow, binHigh, "e");
    rec_selected_window->SetTitle("");

    // Rebin projections with detected truth binning
    if (GetREBINON()) {
      rec_total_window = RebinToTarget(rec_total_window, truthBin, Form("rec_total_R%.1f_rebinned", config.R));
      rec_selected_window = RebinToTarget(rec_selected_window, truthBin, Form("rec_selected_R%.1f_rebinned", config.R));
    }
    
    rec_total_window->SetDirectory(0);
    rec_selected_window->SetDirectory(0);
    
    // Calculate kinematic efficiency: rec_selected_window / rec_total_window
    // This is a ratio, so no normalization by events needed
    TH1* jke = (TH1*)rec_selected_window->Clone(Form("jke_R%.1f", config.R));
    jke->Divide(rec_selected_window, rec_total_window, 1., 1., "B");
    jke->SetDirectory(0);
    jke->SetStats(0);
    jke->SetTitle("");
    
    // Get true jet pT range from histogram
    Double_t truePtMin = rec_total_window->GetXaxis()->GetXmin();
    Double_t truePtMax = rec_total_window->GetXaxis()->GetXmax();
    
    sJKEcanvas->cd();
    gPad->SetGrid(0);
    hoptset(*jke, 0, rColors[i], 0, 200, 0.0, 1.1, 1, 1, 2, 20);  // X range: [0, 200], Y range: 0 to 1.1
    hset(*jke, GetJRETitleX(), "Kinematic efficiency", 1.3, 1.0, 0.05, 0.07, 0.01, 0.01, 0.05, 0.05, 510, 510);
    jke->SetMarkerStyle(20);  // Same marker style as efficiency/purity
    jke->SetMarkerColor(rColors[i]);
    jke->SetLineColor(rColors[i]);
    if (sJKEfirst) {
      jke->Draw("pe");
      sJKEfirst = false;
    } else {
      jke->Draw("pesame");
    }
    
    // Draw vertical lines at the analysis pT range boundaries [5, 140]
    const double kAnalysisPtMin = 5.0;
    const double kAnalysisPtMax = 140.0;
    TLine* lineKineL = new TLine(kAnalysisPtMin, 0.00, kAnalysisPtMin, 1.1);
    lineKineL->SetLineColor(kGray + 3);
    lineKineL->SetLineStyle(2);
    lineKineL->SetLineWidth(2);
    lineKineL->Draw("lsame");

    TLine* lineKineR = new TLine(kAnalysisPtMax, 0.00, kAnalysisPtMax, 1.1);
    lineKineR->SetLineColor(kGray + 3);
    lineKineR->SetLineStyle(2);
    lineKineR->SetLineWidth(2);
    lineKineR->Draw("lsame");
    
    // Only add entry once per R value
    bool jkeEntryExists = false;
    TList* jkeEntries = sJKEcanLegend->GetListOfPrimitives();
    if (jkeEntries) {
      TIter next(jkeEntries);
      TObject* obj;
      while ((obj = next())) {
        TLegendEntry* entry = (TLegendEntry*)obj;
        if (entry && TString(entry->GetLabel()) == config.label) {
          jkeEntryExists = true;
          break;
        }
      }
    }
    if (!jkeEntryExists) {
    sJKEcanLegend->AddEntry(jke, config.label.Data(), "pe");
    }
    
    mcFile->Close();
  }
  
  sJKEcanvas->cd();
  sJKEcanLegend->Draw();
  // Add "UE subtracted" label for UE mode (jet plots only)
  AddUESubtractedLabel(sJKEcanvas, 0.4, 0.6);
  if (GetDRAWPLOTS()) {
    sJKEcanvas->Print(Form("%s/KinematicEfficiency_AllR.pdf", outputDir.Data()));
  }
}

void DrawInvariantYield() {
  if (!GetDrawInvariantYield()) return;

  // Get R configurations and colors
  std::vector<RConfig> rConfigs = GetRConfigs();
  std::vector<Color_t> rColors = GetRColorsForConfigs(rConfigs);
  TString outputDir = GetOutputDir();

  // Load closure-derived regularization parameters
  const char* regFile = (GetCurrentMode() == kUE) ? kOptimalRegFileUE : kOptimalRegFileNonUE;
  auto regMap = LoadOptimalRegularization(regFile);

  // Create canvas with main pad + 1 ratio pad (MC/Data only)
  TString modeTag = (GetCurrentMode() == kUE) ? "UE" : "NonUE";
  int uniqueId = ++GetNN();
  TString padName = Form("InvariantYield_%s_AllR_%d", modeTag.Data(), uniqueId);
  TCanvas* yieldCanvas = new TCanvas(padName.Data(), "Invariant Yield", 800, 700);
  yieldCanvas->SetTopMargin(0.);
  yieldCanvas->SetBottomMargin(0.);
  yieldCanvas->SetLeftMargin(0.);
  yieldCanvas->SetRightMargin(0.);
  yieldCanvas->Draw();

  // Pad sizes: main=65%, ratio=35% (ratio pad has bottom margin for x-axis)
  const float mainFrac = 0.65f;
  const float ratioFrac = 0.35f;

  TPad* mainPad = new TPad(Form("invY_main_%d", uniqueId), "", 0.0, ratioFrac, 1.0, 1.0, 0);
  mainPad->SetTopMargin(0.03 / mainFrac);
  mainPad->SetBottomMargin(0.001);
  mainPad->SetLeftMargin(0.15);
  mainPad->SetRightMargin(0.03);
  mainPad->Draw();

  TPad* ratioPad1 = new TPad(Form("invY_rat1_%d", uniqueId), "", 0.0, 0.0, 1.0, ratioFrac, 0);
  ratioPad1->SetTopMargin(0.001);
  ratioPad1->SetBottomMargin(0.12 / ratioFrac);
  ratioPad1->SetLeftMargin(0.15);
  ratioPad1->SetRightMargin(0.03);
  ratioPad1->Draw();

  // Configure pads
  optFili(*mainPad, 0, 1, 0, 1);  // logx=0, logy=1
  optFili(*ratioPad1, 0, 0, 0, 0);  // logx=0, logy=0
  
  TLegend* legDataset = new TLegend(0.55, 0.30, 0.75, 0.45, NULL, "brNDC");
  legDataset->SetTextSize(0.06);
  legDataset->SetBorderSize(0);
   
  TLegend* leg = new TLegend(0.82, 0.40, 0.95, 0.90, NULL, "brNDC");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  
  // Store histograms for ratio plots
  std::vector<TH1*> invYieldDataHistos;  // Data invariant yield histograms
  std::vector<TH1*> invYieldMCHistos;    // MC invariant yield histograms
  bool firstDraw = true;  // track first successful draw for axis setup

  for (size_t i = 0; i < rConfigs.size(); ++i) {
    const auto& config = rConfigs[i];
    
    // Create temporary pad for DrawJetMatching to get unfolded result
    Filipad2* tempPad = new Filipad2(Form("tempUnfold_InvYield_R%.1f", config.R), ++GetNN(), 2, 0.3, 100, 50, 0.7, 1, 1);
    tempPad->Draw();
    TPad* unfoldpad = tempPad->GetPad(1);
    TPad* ratunfoldpad = tempPad->GetPad(2);
    TLegend* tempLeg = new TLegend(0.75, 0.3, 0.95, 0.65, NULL, "brNDC");
    
    Double_t nevtsData = GetNORMEVENTS() ? Nevents(config.GetDataFile().Data(), config.dataDir.Data(), GetEventObj(), 0) : 1.0;
    Double_t nevtsDataUnTrig = GetNORMEVENTS() ? Nevents(config.GetDataFile().Data(), config.dataDir.Data(), GetEventObj(), 1) : 1.0;
    Double_t nevtsMCD = GetNORMEVENTS() ? Nevents(config.GetMcFile().Data(), config.mcDir.Data(), GetEventObj(), 0, config.isJJ) : 1.0;
    Double_t nevtsMCP = GetNORMEVENTS() ? Nevents(config.GetMcFile().Data(), config.mcDir.Data(), GetEventObj(), 1, config.isJJ) : 1.0;
    
    // MC closure k for cross-check logging (data d-vector determines actual k)
    const Int_t svdKOverride = GetSvdKForRun(config.mcRunNumber, regMap);

    // Call DrawJetMatching to get unfolded result
    TH1* hDcorrected = DrawJetMatching(config.GetMcFile().Data(), config.GetDataFile().Data(), config.label.Data(),
                   nevtsData, nevtsDataUnTrig, nevtsMCD, nevtsMCP,
                   nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                   unfoldpad, ratunfoldpad, tempLeg, rColors[i],
                   nullptr, nullptr, nullptr, nullptr, config.dataDir.Data(),
                   config.mcDir.Data(), "",
                   nullptr, nullptr, nullptr,
                   svdKOverride);

    // Close temporary canvas to avoid cluttering the screen
    if (tempPad->C) {
      tempPad->C->Close();
    }
    
    if (hDcorrected) {
      // Create InvYData following DrawJetsMCfJetQA.h approach (line 906-908)
      // hDcorrected is unfolded but not normalized yet (DrawJetMatching doesn't normalize)
      // Clone and normalize using hoptset (same as DrawJetsMCfJetQA.h)
      TH1* InvYData = (TH1F*)hDcorrected->Clone(Form("invYield_Data_R%.1f", config.R));
      InvYData->SetDirectory(0);
      InvYData->SetStats(0);
      
      // Apply normalization and styling (same as DrawJetsMCfJetQA.h line 908)
      // Normalize by events and bin width
      if (nevtsData > 0) {
        InvYData->Scale(1.0 / nevtsData, "width");
      }
      InvYData->SetLineColor(rColors[i]);
      InvYData->SetMarkerColor(rColors[i]);
      InvYData->SetMarkerStyle(20);
      InvYData->SetMarkerSize(0.6);
      InvYData->SetLineStyle(1);
      InvYData->SetLineWidth(2);
      InvYData->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
      InvYData->SetStats(0);
      
      Double_t deltaEta = GetDeltaEta(config.R);
      if (deltaEta > 0) {
        InvYData->Scale(1.0 / deltaEta);
      }
      
      // Get MC histogram for invariant yield
      TFile* mcFile = TFile::Open(config.GetMcFile().Data(), "READ");
      if (mcFile && !mcFile->IsZombie()) {
        TH1* JetMCPPt = (TH1*)mcFile->Get(Form("%s/%s", config.mcDir.Data(), GetJetPtMCPObj()));
        if (JetMCPPt) {
          if (GetREBINON()) {
            PtBinning truthBin = DetectTruthBinning(JetMCPPt);
            JetMCPPt = RebinToTarget(JetMCPPt, truthBin, Form("JetMCPPt_R%.1f", config.R));
          }
          JetMCPPt->SetDirectory(0);

          // Create InvYMC following DrawJetsMCfJetQA.h approach (line 913-915)
          TH1* InvYMC = (TH1*)JetMCPPt->Clone(Form("invYield_MC_R%.1f", config.R));
          InvYMC->SetDirectory(0);
          InvYMC->SetStats(0);
          
          // Use nevtsMCP (from Nevents(...,1)) which correctly prioritizes
          // weighted MC collision counts for JJ MC (pT-hat weighted spectra)
          Double_t Nmccollcount = nevtsMCP;

          // Apply normalization and styling
          // Normalize by events and bin width
          if (Nmccollcount > 0) {
            InvYMC->Scale(1.0 / Nmccollcount, "width");
          }
          InvYMC->SetLineColor(rColors[i]);
          InvYMC->SetMarkerColor(rColors[i]);
          InvYMC->SetMarkerStyle(33);
          InvYMC->SetMarkerSize(0.75);
          InvYMC->SetLineStyle(1);
          InvYMC->SetLineWidth(2);
          InvYMC->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
          InvYMC->SetStats(0);
          
          Double_t deltaEta = GetDeltaEta(config.R);
          if (deltaEta > 0) {
            InvYMC->Scale(1.0 / deltaEta);
          }
          
          // Store MC histogram
          invYieldMCHistos.push_back(InvYMC);
        } else {
          invYieldMCHistos.push_back(nullptr);  // keep aligned with invYieldDataHistos
        }
        mcFile->Close();
      } else {
        invYieldMCHistos.push_back(nullptr);  // keep aligned with invYieldDataHistos
      }
      
      // Draw on main pad
      mainPad->cd();
      gPad->SetGrid(0);  // Disable grid on pad
      gPad->SetGridx(0);  // Disable x grid on pad
      gPad->SetGridy(0);  // Disable y grid on pad
      if (firstDraw) {
        // Set y-axis title for invariant yield: dN/(dpT dη)
        hset(*InvYData, GetJetPtGenTitleX(), "1/N_{evt} dN_{jet}/d#it{p}_{T}d#it{#eta}", 0.9, 0.9, 0.05, 0.07, 0.01, 0.01, 0.05, 0.05, 510, 505);
        InvYData->Draw("pe");
      } else {
        InvYData->Draw("pesame");
      }
      
      // Draw MC particle level (dashed lines) on main pad immediately after Data
      if (invYieldMCHistos.size() > i && invYieldMCHistos[i]) {
        TH1* InvYMC = invYieldMCHistos[i];
        InvYMC->SetLineStyle(2);  // Dashed line
        InvYMC->SetLineWidth(2);
        InvYMC->SetMarkerStyle(0);  // No markers for MC
        InvYMC->Draw("plsame");
      }
      
      // Add both Data and MC to legend together for each R
      leg->AddEntry(InvYData, config.label.Data(), "lpe");
      if (firstDraw && invYieldMCHistos.size() > i && invYieldMCHistos[i]) {
        legDataset->AddEntry(InvYData, "ALICE data", "lpe");
        TString mcLegLabel = (GetCurrentMode() == kUE) ? "MC truth" : "PYTHIA8 Monash";
        legDataset->AddEntry(invYieldMCHistos[i], mcLegLabel.Data(), "lpe");
      }
      
      // Store for ratio plots
      invYieldDataHistos.push_back(InvYData);
      firstDraw = false;
    } else {
      std::cerr << "[Error] Failed to get unfolded histogram for R=" << config.R << std::endl;
      invYieldDataHistos.push_back(nullptr);  // keep aligned with rConfigs
      invYieldMCHistos.push_back(nullptr);    // keep aligned with rConfigs
    }
  }

  // Draw legend on main pad
  mainPad->cd();
  leg->Draw();
  legDataset->Draw();
  ALICEfigureLegend("ALICE WIP", 0.20, 0.65, 0.42, 0.92, 0.13, 0.03, 0.35, 0.2, 0.06);
  
  // Set y-axis range after all histograms are drawn
  // Find the minimum and maximum values across all Data and MC histograms
  Double_t yMin = 1e10;
  Double_t yMax = 1e-10;
  for (size_t i = 0; i < invYieldDataHistos.size(); ++i) {
    if (invYieldDataHistos[i]) {
      Int_t binLow = invYieldDataHistos[i]->GetXaxis()->FindBin(GetPlotPtMin());
      Int_t binHigh = invYieldDataHistos[i]->GetXaxis()->FindBin(GetPlotPtMax());
      for (Int_t bin = binLow; bin <= binHigh; ++bin) {
        Double_t content = invYieldDataHistos[i]->GetBinContent(bin);
        if (content > 0) {
          yMin = TMath::Min(yMin, content);
          yMax = TMath::Max(yMax, content);
        }
      }
    }
    if (i < invYieldMCHistos.size() && invYieldMCHistos[i]) {
      Int_t binLow = invYieldMCHistos[i]->GetXaxis()->FindBin(GetPlotPtMin());
      Int_t binHigh = invYieldMCHistos[i]->GetXaxis()->FindBin(GetPlotPtMax());
      for (Int_t bin = binLow; bin <= binHigh; ++bin) {
        Double_t content = invYieldMCHistos[i]->GetBinContent(bin);
        if (content > 0) {
          yMin = TMath::Min(yMin, content);
          yMax = TMath::Max(yMax, content);
        }
      }
    }
  }
  
  // Apply some margin (factor of 0.5 for min, 2.0 for max)
  if (yMin < 1e9 && yMax > 1e-9) {
    yMin = yMin * 0.5;
    yMax = yMax * 2.0;
    // Set minimum y-axis value to 1e-8
    yMin = TMath::Max(yMin, 1e-8);  // Minimum y-axis value set to 1e-8
    yMax = TMath::Min(yMax, 1e0);
    
    // Set range on the pad
    mainPad->cd();
    gPad->Update();
    if (gPad->GetLogy()) {
      gPad->SetLogy(0);  // Temporarily turn off log scale
      gPad->Update();
    }
    // Set range on the first drawn histogram (which controls the axis)
    for (size_t ih = 0; ih < invYieldDataHistos.size(); ++ih) {
      if (invYieldDataHistos[ih]) {
        invYieldDataHistos[ih]->GetYaxis()->SetRangeUser(yMin, yMax);
        gPad->SetLogy(1);  // Turn log scale back on
        gPad->Update();
        break;
      }
    }
  }
  
  // Ratio Pad 1: MC/Data for all R values
  if (invYieldMCHistos.size() == invYieldDataHistos.size() && invYieldMCHistos.size() > 0) {
    ratioPad1->cd();
    gPad->SetTicks(1, 1);
    gPad->SetGrid(0);
    gPad->SetGridx(0);
    gPad->SetGridy(0);
    
    // Create frame for ratio plot
    TH1F* ratioFrame1 = new TH1F(Form("ratioFrame1_Data_MC_%d", uniqueId), "", 100, GetPlotPtMin(), GetPlotPtMax());
    ratioFrame1->SetDirectory(nullptr);
    ratioFrame1->SetStats(0);
    ratioFrame1->GetXaxis()->SetTitle(GetJetPtGenTitleX());
    ratioFrame1->GetYaxis()->SetTitle("MC / Data");
    ratioFrame1->GetXaxis()->SetTitleSize(0.1);
    ratioFrame1->GetYaxis()->SetTitleSize(0.1);
    ratioFrame1->GetXaxis()->SetLabelSize(0.1);
    ratioFrame1->GetYaxis()->SetLabelSize(0.1);
    ratioFrame1->GetYaxis()->SetTitleOffset(0.675);
    ratioFrame1->SetMinimum(0.5);
    ratioFrame1->SetMaximum(3.1);
    ratioFrame1->SetLineColor(0);
    ratioFrame1->SetMarkerSize(0);
    ratioFrame1->Draw("AXIS");

    // Draw MC/Data ratio for all R values
    for (size_t i = 0; i < invYieldMCHistos.size(); ++i) {
      if (!invYieldMCHistos[i] || !invYieldDataHistos[i]) continue;
      TH1* ratioHist = (TH1*)invYieldMCHistos[i]->Clone(Form("ratio_MC_Data_R%.1f", rConfigs[i].R));
      ratioHist->Divide(ratioHist, invYieldDataHistos[i], 1., 1., "");
      ratioHist->SetStats(0);
      hset(*ratioHist, GetJetPtGenTitleX(), "MC / Data", 1.3, 0.675, 0.11, 0.1, 0.01, 0.01, 0.1, 0.1, 510, 503);
      ratioHist->SetLineColor(rColors[i]);
      ratioHist->SetMarkerColor(rColors[i]);
      ratioHist->SetMarkerStyle(20);
      ratioHist->SetLineWidth(2);
      ratioHist->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
      ratioHist->GetYaxis()->SetRangeUser(0.5, 3.1);
      ratioHist->Draw("pesame");
    }
    
    // Draw horizontal line at 1.0
    TLine* line1 = new TLine(GetPlotPtMin(), 1.0, GetPlotPtMax(), 1.0);
    line1->SetLineColor(kBlack);
    line1->SetLineStyle(2);
    line1->SetLineWidth(1);
    line1->Draw("lsame");
  }
  
  // Add "UE subtracted" label for UE mode (between ALICE WIP and R labels)
  mainPad->cd();
  AddUESubtractedLabel(mainPad, 0.55, 0.85, 0.05);
  
  if (GetDRAWPLOTS()) {
    yieldCanvas->Print(Form("%s/InvariantYield_%s_AllR.pdf", outputDir.Data(), modeTag.Data()));
  }
}

void DrawCustomRatios() {
  if (!GetDrawCustomRatios()) {
    std::cout << "[DrawCustomRatios] GetDrawCustomRatios() is false. Returning." << std::endl;
    return;
  }

  // Track per-mode calls: allow once per mode (NonUE and UE independently)
  // Get custom ratio pairs
  std::vector<std::pair<double, double>> ratioPairs = GetCustomRatioPairs();
  std::cout << "[DrawCustomRatios] Number of ratio pairs: " << ratioPairs.size() << std::endl;
  if (ratioPairs.empty()) {
    std::cout << "[DrawCustomRatios] No ratio pairs specified. Returning." << std::endl;
    return;
  }
  
  for (const auto& pair : ratioPairs) {
    std::cout << "[DrawCustomRatios] Ratio pair: R=" << pair.first << " / R=" << pair.second << std::endl;
  }
  
  // Get R configurations and colors
  std::vector<RConfig> rConfigs = GetRConfigs();
  std::vector<Color_t> rColors = GetRColorsForConfigs(rConfigs);
  TString outputDir = GetOutputDir();

  const EXsecMode mode = GetCurrentMode();

  // Load closure-derived regularization parameters
  const char* regFileCustom = (mode == kUE) ? kOptimalRegFileUE : kOptimalRegFileNonUE;
  auto regMapCustom = LoadOptimalRegularization(regFileCustom);

  // Create canvas with single pad for ratio plot
  TString padName = Form("CustomRatios_%d", ++GetNN());
  TCanvas* ratioCanvas = new TCanvas(padName.Data(), "Custom Ratios", 800, 600);
  setpad(ratioCanvas, 0.02, 0.12, 0.12, 0.01);
  ratioCanvas->SetGrid(0);
  ratioCanvas->SetGridx(0);
  ratioCanvas->SetGridy(0);
  ratioCanvas->Draw();
  ratioCanvas->cd();
  gPad->SetGrid(0);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  
  // Store histograms for all R values
  std::vector<TH1*> invYieldDataHistos;
  std::vector<double> rValues;

  if (mode == kUE) {
    // Try stitched first, then MB, then JJ as fallback
    const std::map<double, TH1*>* useMap = nullptr;
    TString sourceLabel;
    if (!gUE_Xsec_Stitched.empty()) {
      useMap = &gUE_Xsec_Stitched;
      sourceLabel = "stitched (MB+JJ)";
    } else if (!gUE_Xsec_MB.empty()) {
      useMap = &gUE_Xsec_MB;
      sourceLabel = "MB only (no JJ for stitching)";
    } else if (!gUE_Xsec_JJ.empty()) {
      useMap = &gUE_Xsec_JJ;
      sourceLabel = "JJ only";
    } else {
      std::cerr << "[DrawCustomRatios] UE mode: no cached spectra available. Run DrawCrossSection() first." << std::endl;
      return;
    }
    std::cout << "[DrawCustomRatios] UE mode: using " << sourceLabel.Data() << " spectra for ratios." << std::endl;

    // Build map from R → best MC config (prefer JJ for high-pT statistics)
    std::vector<RConfig> rConfigs_ue = GetRConfigs();
    std::map<double, const RConfig*> bestMCByR;
    for (const auto& cfg : rConfigs_ue) {
      double rKey = TMath::Nint(cfg.R * 10) / 10.0;
      bool isJJ = cfg.isJJ;
      if (bestMCByR.find(rKey) == bestMCByR.end() || isJJ) {
        bestMCByR[rKey] = &cfg;
      }
    }

    for (const auto& kv : *useMap) {
      const double R = kv.first;
      rValues.push_back(R);
      TH1* h = (TH1*)kv.second->Clone(Form("run3_xsec_R%.1f", R));
      h->SetDirectory(0);
      h->SetStats(0);
      invYieldDataHistos.push_back(h);
    }
  } else {
    // Load histograms for all R values (legacy path: unfold inside this function)
    std::cout << "[DrawCustomRatios] Loading histograms for " << rConfigs.size() << " R values..." << std::endl;
    for (size_t i = 0; i < rConfigs.size(); ++i) {
      const auto& config = rConfigs[i];
      rValues.push_back(config.R);
      std::cout << "[DrawCustomRatios] Processing R=" << config.R << "..." << std::endl;

      // Create temporary pad for DrawJetMatching to get unfolded result
      Filipad2* tempPad = new Filipad2(Form("tempUnfold_CustomRatios_R%.1f", config.R), ++GetNN(), 2, 0.3, 100, 50, 0.7, 1, 1);
      tempPad->Draw();
      TPad* unfoldpad = tempPad->GetPad(1);
      TPad* ratunfoldpad = tempPad->GetPad(2);
      TLegend* tempLeg = new TLegend(0.7, 0.9, 0.9, 0.95, NULL, "brNDC");

      Double_t nevtsData = GetNORMEVENTS() ? Nevents(config.GetDataFile().Data(), config.dataDir.Data(), GetEventObj(), 0) : 1.0;
      Double_t nevtsDataUnTrig = GetNORMEVENTS() ? Nevents(config.GetDataFile().Data(), config.dataDir.Data(), GetEventObj(), 1) : 1.0;
      Double_t nevtsMCD = GetNORMEVENTS() ? Nevents(config.GetMcFile().Data(), config.mcDir.Data(), GetEventObj(), 0, config.isJJ) : 1.0;
      Double_t nevtsMCP = GetNORMEVENTS() ? Nevents(config.GetMcFile().Data(), config.mcDir.Data(), GetEventObj(), 1, config.isJJ) : 1.0;

      // MC closure k for cross-check logging (data d-vector determines actual k)
      const Int_t svdKOverrideCustom = GetSvdKForRun(config.mcRunNumber, regMapCustom);

      // Call DrawJetMatching to get unfolded result
      TH1* hDcorrected = DrawJetMatching(config.GetMcFile().Data(), config.GetDataFile().Data(), config.label.Data(),
                                        nevtsData, nevtsDataUnTrig, nevtsMCD, nevtsMCP,
                                        nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                                        unfoldpad, ratunfoldpad, tempLeg, rColors[i],
                                        nullptr, nullptr, nullptr, nullptr, config.dataDir.Data(),
                                        config.mcDir.Data(), "",
                                        nullptr, nullptr, nullptr,
                                        svdKOverrideCustom);

      // Close temporary canvas to avoid cluttering the screen
      if (tempPad->C) {
        tempPad->C->Close();
      }

      if (hDcorrected) {
        // Create InvYData
        TH1* InvYData = (TH1F*)hDcorrected->Clone(Form("invYield_Data_R%.1f", config.R));
        InvYData->SetDirectory(0);
        InvYData->SetStats(0);
        hoptset(*InvYData, nevtsData, rColors[i], GetPlotPtMin(), GetPlotPtMax(), 6e-9, 5e-3, 0.6, 1, 2, 20);

        Double_t deltaEta = GetDeltaEta(config.R);
        if (deltaEta > 0) {
          InvYData->Scale(1.0 / deltaEta);
        }

        invYieldDataHistos.push_back(InvYData);
      } else {
        std::cerr << "[DrawCustomRatios] Failed to get unfolded histogram for R=" << config.R << std::endl;
        invYieldDataHistos.push_back(nullptr);
      }
    }
  }
  
  std::cout << "[DrawCustomRatios] Loaded " << invYieldDataHistos.size() << " Data histograms" << std::endl;
  
  // Switch to ratio canvas
  ratioCanvas->cd();
  gPad->SetGrid(0);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  
  // Create frame for ratio plot
  TH1F* ratioFrame = new TH1F(Form("ratioFrame_Custom_%d", (int)mode), "", 100, GetPlotPtMin(), GetPlotPtMax());
  ratioFrame->SetDirectory(nullptr);
  ratioFrame->SetStats(0);
  ratioFrame->SetMinimum(0.5);
  ratioFrame->SetMaximum(1.5);
  ratioFrame->GetXaxis()->SetTitle(GetJetPtGenTitleX());
  ratioFrame->GetYaxis()->SetTitle("Cross section ratio");
  hset(*ratioFrame, GetJetPtGenTitleX(), "Cross section ratio", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
  ratioFrame->Draw("AXIS");
  
  // Helper function to find R index
  auto findRIndex = [&](double r) -> int {
    for (size_t i = 0; i < rValues.size(); ++i) {
      if (TMath::Abs(rValues[i] - r) < 0.01) return i;
    }
    return -1;
  };
  
  // Hard-coded Run2 reference data from CSV file
  // Structure: {pT, ratio, stat_err+, stat_err-, sys_err+, sys_err-}
  // UE subtracted data (backup)
  struct Run2DataPointUEsub {
    double pT, ratio, statErrP, statErrM, sysErrP, sysErrM;
  };
  
  std::map<std::pair<double, double>, std::vector<Run2DataPointUEsub>> run2DataUEsub;
  
  // Without UE subtraction data (from HEPData-ins2026265-v1-Figure_A3.csv)
  struct Run2DataPoint {
    double pT, ratio, statErrP, statErrM, sysErrP, sysErrM;
  };
  
  std::map<std::pair<double, double>, std::vector<Run2DataPoint>> run2Data;
  
  // R=0.2/R=0.3 (UE subtracted)
  run2DataUEsub[{0.2, 0.3}] = {
    {5.5, 0.68314, 0.00046278, 0.00046278, 0.0502, 0.0502},
    {6.5, 0.64336, 0.00054781, 0.00054781, 0.043892, 0.043892},
    {7.5, 0.62269, 0.00064853, 0.00064853, 0.039569, 0.039569},
    {8.5, 0.60581, 0.00073851, 0.00073851, 0.048042, 0.048042},
    {9.5, 0.59767, 0.00080345, 0.00080345, 0.042282, 0.042282},
    {11.0, 0.59766, 0.00088933, 0.00088933, 0.038594, 0.038594},
    {13.0, 0.60601, 0.0010341, 0.0010341, 0.037348, 0.037348},
    {15.0, 0.61398, 0.0012314, 0.0012314, 0.036927, 0.036927},
    {17.0, 0.61263, 0.0014032, 0.0014032, 0.035343, 0.035343},
    {19.0, 0.63468, 0.0016345, 0.0016345, 0.034681, 0.034681},
    {22.5, 0.65715, 0.0018511, 0.0018511, 0.033809, 0.033809},
    {27.5, 0.684, 0.0023689, 0.0023689, 0.032576, 0.032576},
    {35.0, 0.71333, 0.0030635, 0.0030635, 0.030877, 0.030877},
    {45.0, 0.74709, 0.0038978, 0.0038978, 0.029764, 0.029764},
    {55.0, 0.76908, 0.0047218, 0.0047218, 0.029595, 0.029595},
    {65.0, 0.77983, 0.0054505, 0.0054505, 0.029515, 0.029515},
    {77.5, 0.79118, 0.006104, 0.006104, 0.030538, 0.030538},
    {92.5, 0.79969, 0.0066222, 0.0066222, 0.030401, 0.030401},
    {120.0, 0.81039, 0.0070234, 0.0070234, 0.030919, 0.030919}
  };
  
  // R=0.2/R=0.4 (UE subtracted)
  run2DataUEsub[{0.2, 0.4}] = {
    {5.5, 0.52233, 0.00035637, 0.00035637, 0.038723, 0.038723},
    {6.5, 0.47243, 0.00040424, 0.00040424, 0.033965, 0.033965},
    {7.5, 0.44482, 0.00046468, 0.00046468, 0.031951, 0.031951},
    {8.5, 0.4256, 0.0005201, 0.0005201, 0.031176, 0.031176},
    {9.5, 0.418, 0.00056245, 0.00056245, 0.031028, 0.031028},
    {11.0, 0.4205, 0.00062535, 0.00062535, 0.030845, 0.030845},
    {13.0, 0.43513, 0.00074214, 0.00074214, 0.03063, 0.03063},
    {15.0, 0.45195, 0.00090707, 0.00090707, 0.030411, 0.030411},
    {17.0, 0.45731, 0.0010482, 0.0010482, 0.029158, 0.029158},
    {19.0, 0.47998, 0.001236, 0.001236, 0.028849, 0.028849},
    {22.5, 0.51128, 0.0014389, 0.0014389, 0.029133, 0.029133},
    {27.5, 0.55206, 0.00191, 0.00191, 0.030007, 0.030007},
    {35.0, 0.59232, 0.0025402, 0.0025402, 0.030901, 0.030901},
    {45.0, 0.6374, 0.0033183, 0.0033183, 0.032726, 0.032726},
    {55.0, 0.66496, 0.0040705, 0.0040705, 0.033947, 0.033947},
    {65.0, 0.67877, 0.004727, 0.004727, 0.034348, 0.034348},
    {77.5, 0.69286, 0.0053231, 0.0053231, 0.035349, 0.035349},
    {92.5, 0.70625, 0.0058216, 0.0058216, 0.035199, 0.035199},
    {120.0, 0.7198, 0.006208, 0.006208, 0.035676, 0.035676}
  };
  
  // R=0.2/R=0.5 (UE subtracted)
  run2DataUEsub[{0.2, 0.5}] = {
    {5.5, 0.42858, 0.00029716, 0.00029716, 0.031986, 0.031986},
    {6.5, 0.37563, 0.00032613, 0.00032613, 0.025796, 0.025796},
    {7.5, 0.34439, 0.00036425, 0.00036425, 0.021359, 0.021359},
    {8.5, 0.32548, 0.00040274, 0.00040274, 0.025197, 0.025197},
    {9.5, 0.31784, 0.00043227, 0.00043227, 0.022074, 0.022074},
    {11.0, 0.31981, 0.00047938, 0.00047938, 0.020959, 0.020959},
    {13.0, 0.33588, 0.00057744, 0.00057744, 0.021662, 0.021662},
    {15.0, 0.35634, 0.0007222, 0.0007222, 0.022075, 0.022075},
    {17.0, 0.36584, 0.00084658, 0.00084658, 0.021622, 0.021622},
    {19.0, 0.38616, 0.001002, 0.001002, 0.022465, 0.022465},
    {22.5, 0.42311, 0.0011973, 0.0011973, 0.025147, 0.025147},
    {27.5, 0.47237, 0.0016427, 0.0016427, 0.028537, 0.028537},
    {35.0, 0.51944, 0.0022383, 0.0022383, 0.031298, 0.031298},
    {45.0, 0.57154, 0.0029879, 0.0029879, 0.035063, 0.035063},
    {55.0, 0.60101, 0.0036917, 0.0036917, 0.037512, 0.037512},
    {65.0, 0.61706, 0.004309, 0.004309, 0.03874, 0.03874},
    {77.5, 0.63189, 0.0048653, 0.0048653, 0.040602, 0.040602},
    {92.5, 0.64751, 0.0053469, 0.0053469, 0.041245, 0.041245},
    {120.0, 0.66338, 0.00573, 0.00573, 0.042205, 0.042205}
  };
  
  // R=0.2/R=0.6 (UE subtracted)
  run2DataUEsub[{0.2, 0.6}] = {
    {5.5, 0.37484, 0.0002672, 0.0002672, 0.029153, 0.029153},
    {6.5, 0.32346, 0.00028925, 0.00028925, 0.023138, 0.023138},
    {7.5, 0.28829, 0.00031284, 0.00031284, 0.018172, 0.018172},
    {8.5, 0.26638, 0.00033797, 0.00033797, 0.017488, 0.017488},
    {9.5, 0.25821, 0.00035935, 0.00035935, 0.015305, 0.015305},
    {11.0, 0.25945, 0.00039647, 0.00039647, 0.014838, 0.014838},
    {13.0, 0.27488, 0.00048215, 0.00048215, 0.0159, 0.0159},
    {15.0, 0.29531, 0.0006128, 0.0006128, 0.017224, 0.017224},
    {17.0, 0.30738, 0.00072781, 0.00072781, 0.017067, 0.017067},
    {19.0, 0.326, 0.00086189, 0.00086189, 0.017642, 0.017642},
    {22.5, 0.36001, 0.0010333, 0.0010333, 0.019401, 0.019401},
    {27.5, 0.41752, 0.001472, 0.001472, 0.023368, 0.023368},
    {35.0, 0.47105, 0.002058, 0.002058, 0.026521, 0.026521},
    {45.0, 0.53108, 0.0028135, 0.0028135, 0.030766, 0.030766},
    {55.0, 0.56257, 0.0034988, 0.0034988, 0.033679, 0.033679},
    {65.0, 0.57829, 0.004085, 0.004085, 0.034455, 0.034455},
    {77.5, 0.59367, 0.0046203, 0.0046203, 0.036391, 0.036391},
    {92.5, 0.61146, 0.0051004, 0.0051004, 0.037999, 0.037999},
    {120.0, 0.62981, 0.005493, 0.005493, 0.039001, 0.039001}
  };
  
  // R=0.2/R=0.7 (UE subtracted)
  run2DataUEsub[{0.2, 0.7}] = {
    {5.5, 0.33693, 0.00025076, 0.00025076, 0.020893, 0.020893},
    {6.5, 0.29193, 0.00027516, 0.00027516, 0.020688, 0.020688},
    {7.5, 0.25396, 0.00028926, 0.00028926, 0.016758, 0.016758},
    {8.5, 0.22937, 0.00030532, 0.00030532, 0.013686, 0.013686},
    {9.5, 0.21992, 0.00032026, 0.00032026, 0.011708, 0.011708},
    {11.0, 0.21997, 0.00034975, 0.00034975, 0.011464, 0.011464},
    {13.0, 0.23329, 0.00042721, 0.00042721, 0.012311, 0.012311},
    {15.0, 0.25421, 0.00055461, 0.00055461, 0.013596, 0.013596},
    {17.0, 0.26515, 0.00066001, 0.00066001, 0.014422, 0.014422},
    {19.0, 0.28051, 0.00077301, 0.00077301, 0.015018, 0.015018},
    {22.5, 0.30858, 0.0009154, 0.0009154, 0.015877, 0.015877},
    {27.5, 0.37244, 0.0013556, 0.0013556, 0.020495, 0.020495},
    {35.0, 0.4354, 0.0019672, 0.0019672, 0.023542, 0.023542},
    {45.0, 0.50553, 0.0027705, 0.0027705, 0.028148, 0.028148},
    {55.0, 0.53954, 0.0034677, 0.0034677, 0.03112, 0.03112},
    {65.0, 0.55773, 0.0040657, 0.0040657, 0.033802, 0.033802},
    {77.5, 0.57278, 0.004594, 0.004594, 0.034966, 0.034966},
    {92.5, 0.59171, 0.0050811, 0.0050811, 0.037461, 0.037461},
    {120.0, 0.61323, 0.0055019, 0.0055019, 0.038474, 0.038474}
  };
  
  // Without UE subtraction data (from HEPData-ins2026265-v1-Figure_A3.csv)
  // R=0.2/R=0.3
  run2Data[{0.2, 0.3}] = {
    {5.5, 0.58656, 0.00031472, 0.00031472, 0.033157, 0.033157},
    {6.5, 0.55959, 0.00038644, 0.00038644, 0.030232, 0.030232},
    {7.5, 0.54908, 0.0004718, 0.0004718, 0.029003, 0.029003},
    {8.5, 0.5409, 0.00054813, 0.00054813, 0.027684, 0.027684},
    {9.5, 0.54155, 0.00062037, 0.00062037, 0.02728, 0.02728},
    {11.0, 0.54736, 0.00071379, 0.00071379, 0.027315, 0.027315},
    {13.0, 0.55929, 0.00083888, 0.00083888, 0.027866, 0.027866},
    {15.0, 0.57141, 0.00099166, 0.00099166, 0.028928, 0.028928},
    {17.0, 0.57536, 0.0011451, 0.0011451, 0.028821, 0.028821},
    {19.0, 0.6016, 0.0013892, 0.0013892, 0.029381, 0.029381},
    {22.5, 0.62642, 0.0017412, 0.0017412, 0.029574, 0.029574},
    {27.5, 0.65066, 0.0022383, 0.0022383, 0.02896, 0.02896},
    {35.0, 0.687, 0.0029346, 0.0029346, 0.028064, 0.028064},
    {45.0, 0.72693, 0.0031631, 0.0031631, 0.027676, 0.027676},
    {55.0, 0.75096, 0.0037066, 0.0037066, 0.027971, 0.027971},
    {65.0, 0.76363, 0.0041662, 0.0041662, 0.028229, 0.028229},
    {77.5, 0.77648, 0.0045731, 0.0045731, 0.029431, 0.029431},
    {92.5, 0.78611, 0.0048923, 0.0048923, 0.029374, 0.029374},
    {120.0, 0.79747, 0.0051431, 0.0051431, 0.029896, 0.029896}
  };
  
  // R=0.2/R=0.4
  run2Data[{0.2, 0.4}] = {
    {5.5, 0.37799, 0.00020128, 0.00020128, 0.020695, 0.020695},
    {6.5, 0.34545, 0.00023656, 0.00023656, 0.017894, 0.017894},
    {7.5, 0.33112, 0.00028205, 0.00028205, 0.016664, 0.016664},
    {8.5, 0.3232, 0.00032473, 0.00032473, 0.015716, 0.015716},
    {9.5, 0.32528, 0.00036941, 0.00036941, 0.015314, 0.015314},
    {11.0, 0.33425, 0.00043216, 0.00043216, 0.015571, 0.015571},
    {13.0, 0.35265, 0.00052486, 0.00052486, 0.016179, 0.016179},
    {15.0, 0.37449, 0.00064574, 0.00064574, 0.017144, 0.017144},
    {17.0, 0.38793, 0.00076776, 0.00076776, 0.017667, 0.017667},
    {19.0, 0.41673, 0.0009574, 0.0009574, 0.01907, 0.01907},
    {22.5, 0.4516, 0.0012497, 0.0012497, 0.020995, 0.020995},
    {27.5, 0.49188, 0.0016857, 0.0016857, 0.023121, 0.023121},
    {35.0, 0.5387, 0.0022939, 0.0022939, 0.025414, 0.025414},
    {45.0, 0.59231, 0.0025651, 0.0025651, 0.0285, 0.0285},
    {55.0, 0.62381, 0.0030636, 0.0030636, 0.030586, 0.030586},
    {65.0, 0.64158, 0.0034818, 0.0034818, 0.031755, 0.031755},
    {77.5, 0.65867, 0.0038581, 0.0038581, 0.033331, 0.033331},
    {92.5, 0.67419, 0.0041721, 0.0041721, 0.033655, 0.033655},
    {120.0, 0.68899, 0.0044179, 0.0044179, 0.034508, 0.034508}
  };
  
  // R=0.2/R=0.5
  run2Data[{0.2, 0.5}] = {
    {5.5, 0.26441, 0.0001405, 0.0001405, 0.014458, 0.014458},
    {6.5, 0.22951, 0.00015665, 0.00015665, 0.011785, 0.011785},
    {7.5, 0.21233, 0.00018018, 0.00018018, 0.010458, 0.010458},
    {8.5, 0.20361, 0.00020378, 0.00020378, 0.0098259, 0.0098259},
    {9.5, 0.20368, 0.00023033, 0.00023033, 0.0096365, 0.0096365},
    {11.0, 0.21017, 0.0002705, 0.0002705, 0.009819, 0.009819},
    {13.0, 0.22658, 0.00033586, 0.00033586, 0.010297, 0.010297},
    {15.0, 0.24899, 0.00042814, 0.00042814, 0.010864, 0.010864},
    {17.0, 0.26657, 0.00052653, 0.00052653, 0.011699, 0.011699},
    {19.0, 0.29415, 0.00067468, 0.00067468, 0.01398, 0.01398},
    {22.5, 0.33354, 0.00092189, 0.00092189, 0.017651, 0.017651},
    {27.5, 0.38112, 0.0013058, 0.0013058, 0.021564, 0.021564},
    {35.0, 0.43217, 0.0018417, 0.0018417, 0.02513, 0.02513},
    {45.0, 0.49669, 0.0021516, 0.0021516, 0.030004, 0.030004},
    {55.0, 0.53181, 0.0026127, 0.0026127, 0.033083, 0.033083},
    {65.0, 0.55377, 0.0030063, 0.0030063, 0.034908, 0.034908},
    {77.5, 0.57316, 0.0033581, 0.0033581, 0.037143, 0.037143},
    {92.5, 0.59184, 0.0036633, 0.0036633, 0.038141, 0.038141},
    {120.0, 0.60935, 0.0039079, 0.0039079, 0.039285, 0.039285}
  };
  
  // R=0.2/R=0.6
  run2Data[{0.2, 0.6}] = {
    {5.5, 0.2042, 0.00010887, 0.00010887, 0.011876, 0.011876},
    {6.5, 0.16738, 0.00011445, 0.00011445, 0.0091042, 0.0091042},
    {7.5, 0.14844, 0.00012605, 0.00012605, 0.0075066, 0.0075066},
    {8.5, 0.13764, 0.00013779, 0.00013779, 0.0066301, 0.0066301},
    {9.5, 0.13479, 0.00015234, 0.00015234, 0.0063378, 0.0063378},
    {11.0, 0.13679, 0.00017581, 0.00017581, 0.006371, 0.006371},
    {13.0, 0.14719, 0.00021791, 0.00021791, 0.0067901, 0.0067901},
    {15.0, 0.16429, 0.00028249, 0.00028249, 0.0077256, 0.0077256},
    {17.0, 0.18118, 0.00035814, 0.00035814, 0.0083638, 0.0083638},
    {19.0, 0.20592, 0.00047276, 0.00047276, 0.0098578, 0.0098578},
    {22.5, 0.24215, 0.00067032, 0.00067032, 0.012209, 0.012209},
    {27.5, 0.29302, 0.0010068, 0.0010068, 0.015932, 0.015932},
    {35.0, 0.34791, 0.0014897, 0.0014897, 0.019421, 0.019421},
    {45.0, 0.42042, 0.0018301, 0.0018301, 0.024442, 0.024442},
    {55.0, 0.45915, 0.0022675, 0.0022675, 0.027763, 0.027763},
    {65.0, 0.48317, 0.002637, 0.002637, 0.029199, 0.029199},
    {77.5, 0.50479, 0.0029734, 0.0029734, 0.031438, 0.031438},
    {92.5, 0.52642, 0.0032759, 0.0032759, 0.033271, 0.033271},
    {120.0, 0.54656, 0.003524, 0.003524, 0.034429, 0.034429}
  };
  
  // R=0.2/R=0.7
  run2Data[{0.2, 0.7}] = {
    {5.5, 0.17087, 9.2302e-05, 9.2302e-05, 0.0095031, 0.0095031},
    {6.5, 0.13236, 9.1421e-05, 9.1421e-05, 0.006995, 0.006995},
    {7.5, 0.11175, 9.5667e-05, 9.5667e-05, 0.0060683, 0.0060683},
    {8.5, 0.099501, 0.00010033, 0.00010033, 0.0048926, 0.0048926},
    {9.5, 0.094428, 0.00010731, 0.00010731, 0.0044446, 0.0044446},
    {11.0, 0.09276, 0.00011968, 0.00011968, 0.004383, 0.004383},
    {13.0, 0.096838, 0.00014385, 0.00014385, 0.0045254, 0.0045254},
    {15.0, 0.10775, 0.00018617, 0.00018617, 0.0050768, 0.0050768},
    {17.0, 0.11938, 0.00023734, 0.00023734, 0.0059089, 0.0059089},
    {19.0, 0.13817, 0.00031898, 0.00031898, 0.0069896, 0.0069896},
    {22.5, 0.16889, 0.00047038, 0.00047038, 0.0084684, 0.0084684},
    {27.5, 0.21768, 0.00075432, 0.00075432, 0.011935, 0.011935},
    {35.0, 0.27205, 0.0011792, 0.0011792, 0.014829, 0.014829},
    {45.0, 0.35415, 0.0015626, 0.0015626, 0.020004, 0.020004},
    {55.0, 0.39629, 0.0019854, 0.0019854, 0.02325, 0.02325},
    {65.0, 0.42454, 0.0023518, 0.0023518, 0.026182, 0.026182},
    {77.5, 0.44763, 0.002677, 0.002677, 0.027832, 0.027832},
    {92.5, 0.47108, 0.0029766, 0.0029766, 0.030353, 0.030353},
    {120.0, 0.49401, 0.0032344, 0.0032344, 0.031527, 0.031527}
  };
  
  // Draw ratios for each pair
  bool firstDraw = true;
  std::cout << "[DrawCustomRatios] Drawing " << ratioPairs.size() << " ratio pairs..." << std::endl;
  
  // Color array for ratio pairs (different from R colors to distinguish ratios)
  std::vector<Color_t> ratioColors;
  ratioColors.push_back(kRed);
  ratioColors.push_back(kBlue);
  ratioColors.push_back(kGreen+2);
  ratioColors.push_back(kMagenta+2);
  ratioColors.push_back(kCyan+2);
  ratioColors.push_back(kOrange+7);
  ratioColors.push_back(kYellow+2);
  ratioColors.push_back(kPink+2);
  ratioColors.push_back(kSpring+5);
  ratioColors.push_back(kTeal+5);
  
  // Store ratio histograms for double ratio plot
  std::vector<TH1*> ratioHistDataList;
  std::map<std::string, std::vector<TH1*>> ratioHistModelLists;  // per model name
  std::vector<TH1*> ratioHistPOWHEGList;  // Store POWHEG ratio histograms for POWHEG/Data double ratio
  std::vector<std::pair<double, double>> ratioPairsList;
  std::map<std::pair<double, double>, TGraphErrors*> run2Graphs;  // Store Run2 graphs for double ratio
  auto crModels = GetStandaloneModels();

  // Create legend for Canvas A (Run3 vs Run2 only)
  TLegend* commonLeg = new TLegend(0.35, 0.75, 0.60, 0.95, NULL, "brNDC");
  commonLeg->SetTextSize(0.035);
  commonLeg->SetBorderSize(0);

  TH1F* dummyRun3 = new TH1F("dummyRun3", "", 1, 0, 1);
  dummyRun3->SetLineColor(kBlack); dummyRun3->SetMarkerColor(kBlack);
  dummyRun3->SetMarkerStyle(20); dummyRun3->SetLineWidth(2);
  commonLeg->AddEntry(dummyRun3, "Run 3", "lpe");

  TH1F* dummyRun2 = new TH1F("dummyRun2", "", 1, 0, 1);
  dummyRun2->SetFillColorAlpha(kBlack, 0.35); dummyRun2->SetFillStyle(1001);
  dummyRun2->SetMarkerColor(kBlack); dummyRun2->SetMarkerStyle(29);
  dummyRun2->SetLineColor(kBlack); dummyRun2->SetLineStyle(3); dummyRun2->SetLineWidth(2);
  commonLeg->AddEntry(dummyRun2, "Run 2", "lpef");

  // Load POWHEG cross sections for all unique R values needed in ratio pairs
  std::map<double, TH1*> powhegByR;
  {
    std::set<double> neededR;
    for (const auto& pair : ratioPairs) {
      neededR.insert(pair.first);
      neededR.insert(pair.second);
    }
    for (double R : neededR) {
      TH1* hPow = LoadPowhegForR(R);
      if (hPow) {
        powhegByR[R] = hPow;
      }
    }
  }

  // Create separate legend for R ratio pairs (right side)
  TLegend* ratioPairLeg = new TLegend(0.70, 0.744348, 0.90, 0.95, NULL, "brNDC");
  ratioPairLeg->SetTextSize(0.04);
  ratioPairLeg->SetBorderSize(0);
  
  size_t pairIdx = 0;
  for (const auto& pair : ratioPairs) {
    double numR = pair.first;
    double denR = pair.second;
    
    std::cout << "[DrawCustomRatios] Processing ratio R=" << numR << " / R=" << denR << std::endl;
    
    int numIdx = findRIndex(numR);
    int denIdx = findRIndex(denR);
    
    std::cout << "[DrawCustomRatios] Found indices: numIdx=" << numIdx << ", denIdx=" << denIdx << std::endl;
    
    if (numIdx < 0 || denIdx < 0) {
      std::cerr << "[Warning] R value not found: " << numR << " or " << denR << std::endl;
      continue;
    }
    
    if (numIdx >= invYieldDataHistos.size() || denIdx >= invYieldDataHistos.size() ||
        !invYieldDataHistos[numIdx] || !invYieldDataHistos[denIdx]) {
      std::cerr << "[Warning] Histogram not available for R=" << numR << " (idx=" << numIdx << ") or R=" << denR << " (idx=" << denIdx << ")" << std::endl;
      std::cerr << "[Warning] Data histograms size: " << invYieldDataHistos.size() << std::endl;
      continue;
    }
    
    std::cout << "[DrawCustomRatios] Drawing Data ratio..." << std::endl;
    
    // Make sure we're on the ratio canvas
    ratioCanvas->cd();
    
    // Get color for this ratio pair: use denominator R's color
    // e.g. R=0.2/R=0.7 uses R=0.7's color (orange)
    Color_t pairColor = GetColorForR(denR);
    
    // Draw Data ratio (smaller axis labels/titles)
    // Use pairColor for actual graph (legend will use black dummy)
    TH1* ratioHistData = (TH1*)invYieldDataHistos[numIdx]->Clone(Form("ratio_Data_R%.1f_R%.1f", numR, denR));
    // Use standard error propagation for A/B (do not use binomial option)
    ratioHistData->Divide(ratioHistData, invYieldDataHistos[denIdx], 1., 1., "");
    ratioHistData->SetStats(0);
    hset(*ratioHistData, "#it{p}_{T, jet}^{ch} (GeV/#it{c})", Form("#sigma(#it{R}=%.1f) / #sigma(#it{R}=X)", numR), 1.0, 1.0, 0.05, 0.06, 0.01, 0.01, 0.045, 0.045, 510, 505);
    ratioHistData->SetLineColor(pairColor);
    ratioHistData->SetMarkerColor(pairColor);
    ratioHistData->SetMarkerStyle(20);  // Full circle
    ratioHistData->SetLineWidth(2);
    ratioHistData->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
    ratioHistData->GetYaxis()->SetRangeUser(0.0, 1.5);
    ratioHistData->Draw(firstDraw ? "pe" : "pesame");
    
    // Store for double ratio plot
    ratioHistDataList.push_back((TH1*)ratioHistData->Clone(Form("ratio_Data_R%.1f_R%.1f_stored", numR, denR)));
    ratioPairsList.push_back(pair);
    
    // Compute standalone model R-ratios (store for Canvas B and combined canvas)
    for (auto& m : crModels) {
      TH1* hModelNum = LoadStandaloneModelXsec(m, numR);
      TH1* hModelDen = LoadStandaloneModelXsec(m, denR);
      TH1* hModelRatio = nullptr;
      if (hModelNum && hModelDen) {
        hModelRatio = (TH1*)hModelNum->Clone(Form("hModelRatio_%s_R%.1f_R%.1f", m.name.Data(), numR, denR));
        hModelRatio->Divide(hModelNum, hModelDen, 1., 1., "");
        hModelRatio->SetDirectory(0);
      }
      if (hModelNum) delete hModelNum;
      if (hModelDen) delete hModelDen;
      if (hModelRatio) {
        ratioHistModelLists[m.name.Data()].push_back(hModelRatio);
      } else {
        ratioHistModelLists[m.name.Data()].push_back(nullptr);
      }
    }

    // Compute POWHEG NLO ratio (store for Canvas B and combined canvas)
    if (powhegByR.count(numR) && powhegByR.count(denR)) {
      TH1* hPowNum = powhegByR[numR];
      TH1* hPowDen = powhegByR[denR];
      TH1* ratioHistPOWHEG = (TH1*)hPowNum->Clone(Form("ratio_POWHEG_R%.1f_R%.1f", numR, denR));
      ratioHistPOWHEG->Divide(hPowNum, hPowDen, 1., 1., "");
      ratioHistPOWHEG->SetDirectory(0);
      ratioHistPOWHEGList.push_back(ratioHistPOWHEG);
    } else {
      ratioHistPOWHEGList.push_back(nullptr);
    }

    // Draw Run2 reference data on Canvas A
    TGraphErrors* run2Graph = nullptr;

    auto drawRun2Reference = [&](const auto& dataPoints) {
      int nPoints = dataPoints.size();
      
      // Create TGraphErrors for statistical errors (with data points)
      // Use pairColor for actual graph (legend will use black dummy)
      run2Graph = new TGraphErrors(nPoints);
      run2Graph->SetName(Form("run2Graph_R%.1f_R%.1f", numR, denR));
      run2Graph->SetLineColor(pairColor);
      run2Graph->SetMarkerColor(pairColor);
      run2Graph->SetMarkerStyle(29);  // Star
      run2Graph->SetLineStyle(3);     // Dotted line
      run2Graph->SetLineWidth(2);
      run2Graph->SetMarkerSize(1.2);
      // Also give it the same fill style as the syst box so that the legend
      // entry shows a box + star (even though the actual boxes are drawn by run2SysBox)
      run2Graph->SetFillColorAlpha(pairColor, 0.35);
      run2Graph->SetFillStyle(1001);
      
      // Create TGraphAsymmErrors for systematic uncertainty boxes (without points)
      TGraphAsymmErrors* run2SysBox = new TGraphAsymmErrors(nPoints);
      run2SysBox->SetName(Form("run2SysBox_R%.1f_R%.1f", numR, denR));
      run2SysBox->SetFillColorAlpha(pairColor, 0.35);  // Slightly higher opacity
      run2SysBox->SetFillStyle(1001);
      run2SysBox->SetLineWidth(0);
      run2SysBox->SetLineColor(0);
      
      // Hard-coded bin edges matching Data/MC histogram binning (ptbin from GetPtbin())
      // ptbin = {5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200}
      // Run 2 data points: 5.5, 6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 17.0, 19.0, 22.5, 27.5, 35.0, 45.0, 55.0, 65.0, 77.5, 92.5, 120.0
      // Each point corresponds to a bin center, so bin edges are from ptbin array
      static const double ptbinEdges[] = {5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
      static const double run2PtValues[] = {5.5, 6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 17.0, 19.0, 22.5, 27.5, 35.0, 45.0, 55.0, 65.0, 77.5, 92.5, 120.0};
      
      for (int i = 0; i < nPoints; ++i) {
        const auto& pt = dataPoints[i];
        run2Graph->SetPoint(i, pt.pT, pt.ratio);
        run2Graph->SetPointError(i, 0, pt.statErrP);  // Statistical error (y only)
        
        // Find corresponding bin edges from ptbin array
        // Match pt.pT to run2PtValues to find bin index
        int binIdx = -1;
        for (int j = 0; j < 19; ++j) {
          if (TMath::Abs(pt.pT - run2PtValues[j]) < 0.1) {
            binIdx = j;
            break;
          }
        }
        
        double exL = 0.5, exR = 0.5;  // Default values
        if (binIdx >= 0 && binIdx < 19) {
          // Find bin edges: binIdx corresponds to bin [ptbinEdges[binIdx], ptbinEdges[binIdx+1])
          if (binIdx < 20) {
            double binLeft = ptbinEdges[binIdx];
            double binRight = ptbinEdges[binIdx + 1];
            exL = pt.pT - binLeft;
            exR = binRight - pt.pT;
          }
        }
        
        // Systematic uncertainty box
        run2SysBox->SetPoint(i, pt.pT, pt.ratio);
        run2SysBox->SetPointError(i, exL, exR, pt.sysErrM, pt.sysErrP);
      }
      
      // Draw systematic uncertainty boxes first (behind everything)
      run2SysBox->Draw("E2 same");
      
      // Draw statistical error with data points
      run2Graph->Draw("pz same");
      
      // Store for double ratio plot
      run2Graphs[pair] = run2Graph;
    };

    if (mode == kUE) {
      auto it = run2DataUEsub.find(pair);
      if (it != run2DataUEsub.end()) {
        drawRun2Reference(it->second);
      }
    } else {
      auto it = run2Data.find(pair);
      if (it != run2Data.end()) {
        drawRun2Reference(it->second);
      }
    }
    
    // Add R ratio pair to separate legend with colored marker
    TH1F* dummyPair = new TH1F(Form("dummyPair_R%.1f_R%.1f", numR, denR), "", 1, 0, 1);
    dummyPair->SetLineColor(pairColor);
    dummyPair->SetMarkerColor(pairColor);
    dummyPair->SetMarkerStyle(20);
    dummyPair->SetLineWidth(2);
    ratioPairLeg->AddEntry(dummyPair, Form("#it{R}=%.1f / #it{R}=%.1f", numR, denR), "LEP");
    
    firstDraw = false;
    pairIdx++;
  }
  
  // Draw both legends after all pairs are processed
  ratioCanvas->cd();
  commonLeg->Draw();  // Left side: fixed entries
  ratioPairLeg->Draw();  // Right side: R ratio pairs
  
  // Draw ALICE figure legend on ratio canvas
  ratioCanvas->cd();
  ALICEfigureLegend("ALICE WIP", 0.15, 0.60, 0.35, 0.90, 0.70, 0.05, 0.90, 0.15, 0.04);
  
  // Update canvas to ensure all plots are visible
  ratioCanvas->Update();
  ratioCanvas->Modified();
  ratioCanvas->Draw();
  
  if (GetDRAWPLOTS()) {
    ratioCanvas->Print(Form("%s/CustomRatios_Run2.pdf", outputDir.Data()));
  }

  // ========== Canvas A2: Double Ratio Run3/Run2 ==========
  if (!run2Graphs.empty() && !ratioHistDataList.empty()) {
    TCanvas* drCanvas = new TCanvas(Form("DoubleRatio_Run3Run2_%d", ++GetNN()),
                                    "Double Ratio Run3/Run2", 800, 500);
    setpad(drCanvas, 0.02, 0.14, 0.14, 0.01);
    drCanvas->SetGrid(0);
    drCanvas->Draw();
    drCanvas->cd();

    TH1F* drFrame = new TH1F(Form("drFrame_%d", GetNN()), "", 100,
                               GetPlotPtMin(), GetPlotPtMax());
    drFrame->SetMinimum(0.7);
    drFrame->SetMaximum(1.3);
    hset(*drFrame, "#it{p}_{T, jet}^{ch} (GeV/#it{c})",
         "R-ratio Run 3 / Run 2", 1.0, 1.0, 0.05, 0.06, 0.01, 0.01, 0.045, 0.045, 510, 505);
    drFrame->Draw("AXIS");

    TLine* drUnity = new TLine(GetPlotPtMin(), 1.0, GetPlotPtMax(), 1.0);
    drUnity->SetLineStyle(2); drUnity->SetLineWidth(1);
    drUnity->Draw();

    TLegend* drLeg = new TLegend(0.60, 0.70, 0.90, 0.92, NULL, "brNDC");
    drLeg->SetTextSize(0.035);
    drLeg->SetBorderSize(0);

    for (size_t pi = 0; pi < ratioPairsList.size(); ++pi) {
      auto pair = ratioPairsList[pi];
      TH1* hR3 = (pi < ratioHistDataList.size()) ? ratioHistDataList[pi] : nullptr;
      auto itR2 = run2Graphs.find(pair);
      if (!hR3 || itR2 == run2Graphs.end()) continue;

      TGraphErrors* gR2 = itR2->second;
      Color_t pairColor = GetColorForR(pair.second);

      // Compute double ratio: for each Run2 pT point, find matching Run3 bin
      int nR2 = gR2->GetN();
      std::vector<double> drPt, drVal, drErr;
      for (int i = 0; i < nR2; i++) {
        double r2x, r2y;
        gR2->GetPoint(i, r2x, r2y);
        double r2e = gR2->GetErrorY(i);
        if (r2x > GetPlotPtMax() || r2y <= 0) continue;

        int bin = hR3->FindBin(r2x);
        double r3y = hR3->GetBinContent(bin);
        double r3e = hR3->GetBinError(bin);
        if (r3y <= 0) continue;

        double dr = r3y / r2y;
        double err = dr * TMath::Sqrt(TMath::Power(r3e/r3y, 2) + TMath::Power(r2e/r2y, 2));
        drPt.push_back(r2x);
        drVal.push_back(dr);
        drErr.push_back(err);
      }

      if (drPt.empty()) continue;
      TGraphErrors* gDR = new TGraphErrors(drPt.size(), drPt.data(), drVal.data(), nullptr, drErr.data());
      gDR->SetMarkerStyle(20);
      gDR->SetMarkerSize(0.7);
      gDR->SetMarkerColor(pairColor);
      gDR->SetLineColor(pairColor);
      gDR->Draw("pz same");

      drLeg->AddEntry(gDR, Form("R=%.1f / R=%.1f", pair.first, pair.second), "lpe");
    }

    drLeg->Draw();

    TString modeLabel = (mode == kUE) ? "UE subtracted" : "Inclusive";
    ALICEfigureLegend("ALICE WIP", 0.16, 0.65, 0.30, 0.90, 0.78, 0.05, 0.90, 0.16, 0.035);

    drCanvas->Update();
    if (GetDRAWPLOTS()) {
      TString suffix = (mode == kUE) ? "_UEsub" : "";
      drCanvas->Print(Form("%s/DoubleRatio_Run3vsRun2%s.pdf", outputDir.Data(), suffix.Data()));
    }
  }

  // ========== Canvas B: CustomRatios_Models — R-ratios with MC-truth + models ==========
  {
    TCanvas* modelRatioCanvas = new TCanvas(Form("CustomRatios_Models_%d", ++GetNN()),
                                            "Custom Ratios: Models", 800, 600);
    setpad(modelRatioCanvas, 0.02, 0.12, 0.12, 0.01);
    modelRatioCanvas->SetGrid(0);
    modelRatioCanvas->Draw();

    // Legend for models
    TLegend* modLeg = new TLegend(0.35, 0.55, 0.60, 0.95, NULL, "brNDC");
    modLeg->SetTextSize(0.035);
    modLeg->SetBorderSize(0);
    TH1F* dumR3m = new TH1F(Form("dumR3m_%d", GetNN()), "", 1, 0, 1);
    dumR3m->SetLineColor(kBlack); dumR3m->SetMarkerStyle(20); dumR3m->SetLineWidth(2);
    modLeg->AddEntry(dumR3m, "Run 3 Data", "lpe");
    TH1F* dumMCTm = new TH1F(Form("dumMCTm_%d", GetNN()), "", 1, 0, 1);
    dumMCTm->SetLineColor(kBlack); dumMCTm->SetLineStyle(7); dumMCTm->SetLineWidth(2);
    modLeg->AddEntry(dumMCTm, "MC particle-level", "l");
    for (auto& m : crModels) {
      TH1F* dumM = new TH1F(Form("dumCRM_%s_%d", m.name.Data(), GetNN()), "", 1, 0, 1);
      dumM->SetLineColor(kBlack); dumM->SetLineStyle(m.lineStyle); dumM->SetLineWidth(m.lineWidth);
      modLeg->AddEntry(dumM, m.legLabel.Data(), "l");
    }
    TH1F* dumPm = new TH1F(Form("dumPm_%d", GetNN()), "", 1, 0, 1);
    dumPm->SetLineColor(kBlack); dumPm->SetLineStyle(3); dumPm->SetLineWidth(3);
    modLeg->AddEntry(dumPm, "POWHEG NLO", "l");

    // R-pair legend
    TLegend* modPairLeg = new TLegend(0.70, 0.744, 0.90, 0.95, NULL, "brNDC");
    modPairLeg->SetTextSize(0.04);
    modPairLeg->SetBorderSize(0);

    bool firstModDraw = true;
    for (size_t pi = 0; pi < ratioPairsList.size(); ++pi) {
      double numR = ratioPairsList[pi].first;
      double denR = ratioPairsList[pi].second;
      Color_t pairColor = GetColorForR(denR);

      // Data R-ratio
      if (pi < ratioHistDataList.size() && ratioHistDataList[pi]) {
        TH1* hDR = (TH1*)ratioHistDataList[pi]->Clone(Form("modRatData_%zu", pi));
        hDR->SetLineColor(pairColor); hDR->SetMarkerColor(pairColor);
        hDR->SetMarkerStyle(20); hDR->SetLineWidth(2);
        if (firstModDraw) {
          hset(*hDR, "#it{p}_{T, jet}^{ch} (GeV/#it{c})",
               Form("#sigma(#it{R}=%.1f) / #sigma(#it{R}=X)", numR),
               1.0, 1.0, 0.05, 0.06, 0.01, 0.01, 0.045, 0.045, 510, 505);
          hDR->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
          hDR->GetYaxis()->SetRangeUser(0.0, 1.5);
          hDR->Draw("pe");
          firstModDraw = false;
        } else {
          hDR->Draw("pe same");
        }
      }

      // MC truth R-ratio
      {
        int numIdx = findRIndex(numR);
        int denIdx = findRIndex(denR);
        TH1* hMCTnum = nullptr;
        TH1* hMCTden = nullptr;
        if (mode == kUE) {
          if (gUE_MCTruth.count(numR)) hMCTnum = gUE_MCTruth[numR];
          if (gUE_MCTruth.count(denR)) hMCTden = gUE_MCTruth[denR];
        } else {
          // Non-UE: load from configs
          for (size_t ci = 0; ci < rConfigs.size(); ++ci) {
            if (TMath::Abs(rConfigs[ci].R - numR) < 0.01 && !hMCTnum)
              hMCTnum = LoadMCTruthXsec(rConfigs[ci]);
            if (TMath::Abs(rConfigs[ci].R - denR) < 0.01 && !hMCTden)
              hMCTden = LoadMCTruthXsec(rConfigs[ci]);
          }
        }
        if (hMCTnum && hMCTden) {
          TH1* hMCTratio = (TH1*)hMCTnum->Clone(Form("MCTratio_R%.1f_R%.1f", numR, denR));
          hMCTratio->Divide(hMCTnum, hMCTden, 1., 1., "");
          hMCTratio->SetStats(0);
          hMCTratio->SetLineColor(pairColor); hMCTratio->SetMarkerColor(pairColor);
          hMCTratio->SetLineStyle(7); hMCTratio->SetLineWidth(2); hMCTratio->SetMarkerStyle(0);
          hMCTratio->Draw("l same");
        }
        // Clean up non-UE loaded histograms
        if (mode != kUE) { delete hMCTnum; delete hMCTden; }
      }

      // Standalone model R-ratios
      for (auto& m : crModels) {
        if (ratioHistModelLists.count(m.name.Data()) &&
            pi < ratioHistModelLists[m.name.Data()].size() &&
            ratioHistModelLists[m.name.Data()][pi]) {
          TH1* hMR = (TH1*)ratioHistModelLists[m.name.Data()][pi]->Clone(
            Form("modRat_%s_%zu", m.name.Data(), pi));
          hMR->SetStats(0);
          hMR->SetLineColor(pairColor); hMR->SetMarkerColor(pairColor);
          hMR->SetLineStyle(m.lineStyle); hMR->SetLineWidth(m.lineWidth);
          hMR->SetMarkerStyle(0);
          hMR->Draw("l same");
        }
      }

      // POWHEG R-ratio
      if (pi < ratioHistPOWHEGList.size() && ratioHistPOWHEGList[pi]) {
        TH1* hPR = (TH1*)ratioHistPOWHEGList[pi]->Clone(Form("modRatPow_%zu", pi));
        hPR->SetStats(0);
        hPR->SetLineColor(pairColor); hPR->SetMarkerColor(pairColor);
        hPR->SetLineStyle(3); hPR->SetLineWidth(3); hPR->SetMarkerStyle(0);
        hPR->Draw("l same");
      }

      // R-pair legend entry
      TH1F* dumPair = new TH1F(Form("dumModPair_%zu_%d", pi, GetNN()), "", 1, 0, 1);
      dumPair->SetLineColor(pairColor); dumPair->SetMarkerColor(pairColor);
      dumPair->SetMarkerStyle(20); dumPair->SetLineWidth(2);
      modPairLeg->AddEntry(dumPair, Form("#it{R}=%.1f / #it{R}=%.1f", numR, denR), "LEP");
    }

    modelRatioCanvas->cd();
    modLeg->Draw();
    modPairLeg->Draw();
    ALICEfigureLegend("ALICE WIP", 0.15, 0.50, 0.35, 0.80, 0.70, 0.05, 0.90, 0.15, 0.04);

    if (GetDRAWPLOTS()) {
      modelRatioCanvas->Update(); modelRatioCanvas->Modified();
      modelRatioCanvas->Print(Form("%s/CustomRatios_Models.pdf", outputDir.Data()));
    }
  }

  // ========== CustomRatios_Combined: Horizontal per-R layout, double ratios ==========
  // Each column = one R-ratio pair, showing Model R-ratio / Data R-ratio for all models
  if (ratioHistDataList.size() > 0) {
    int nPairs = (int)ratioPairsList.size();
    // Models + POWHEG + MC-truth
    int nModelRows = (int)crModels.size() + 2;  // +1 POWHEG, +1 MC-truth

    int canW = 400 * nPairs;
    int canH = 400;
    TCanvas* combinedCanvas = new TCanvas(Form("CombinedRatios_%d", ++GetNN()),
                                          "Combined Double Ratios", canW, canH);
    combinedCanvas->Divide(nPairs, 1, 0.001, 0.001);

    for (int ip = 0; ip < nPairs; ++ip) {
      combinedCanvas->cd(ip + 1);
      gPad->SetGrid(0);
      gPad->SetLeftMargin(ip == 0 ? 0.18 : 0.05);
      gPad->SetRightMargin(ip == nPairs - 1 ? 0.05 : 0.01);
      gPad->SetBottomMargin(0.15);
      gPad->SetTopMargin(0.05);

      double numR = ratioPairsList[ip].first;
      double denR = ratioPairsList[ip].second;
      TString panelTitle = Form("#sigma(#it{R}=%.1f)/#sigma(#it{R}=%.1f)", numR, denR);

      TH1F* frame = new TH1F(Form("combFrame_%d_%d", ip, GetNN()), "", 100, 5, 200);
      frame->SetDirectory(nullptr); frame->SetStats(0);
      frame->SetMinimum(0.5); frame->SetMaximum(3.5);
      frame->GetXaxis()->SetTitle(GetJetPtGenTitleX());
      frame->GetYaxis()->SetTitle(ip == 0 ? "Model / Data (R-ratio)" : "");
      frame->GetXaxis()->SetTitleSize(0.06); frame->GetYaxis()->SetTitleSize(0.06);
      frame->GetXaxis()->SetLabelSize(0.05);
      frame->GetYaxis()->SetLabelSize(ip == 0 ? 0.05 : 0.0);
      frame->GetXaxis()->SetTitleOffset(1.1);
      frame->GetYaxis()->SetTitleOffset(1.3);
      frame->GetXaxis()->SetNdivisions(510); frame->GetYaxis()->SetNdivisions(505);
      frame->Draw("AXIS");

      // Panel title
      TLatex* lat = new TLatex(0.5, 0.92, panelTitle.Data());
      lat->SetNDC(); lat->SetTextSize(0.06); lat->SetTextAlign(22);
      lat->Draw();

      // Unity line
      TLine* uLine = new TLine(5, 1.0, 200, 1.0);
      uLine->SetLineColor(kBlack); uLine->SetLineStyle(2); uLine->SetLineWidth(1);
      uLine->Draw("lsame");

      if (!ratioHistDataList[ip]) continue;

      // MC-truth double ratio
      {
        TH1* hMCTnum = nullptr;
        TH1* hMCTden = nullptr;
        if (mode == kUE) {
          if (gUE_MCTruth.count(numR)) hMCTnum = gUE_MCTruth[numR];
          if (gUE_MCTruth.count(denR)) hMCTden = gUE_MCTruth[denR];
        } else {
          for (size_t ci = 0; ci < rConfigs.size(); ++ci) {
            if (TMath::Abs(rConfigs[ci].R - numR) < 0.01 && !hMCTnum)
              hMCTnum = LoadMCTruthXsec(rConfigs[ci]);
            if (TMath::Abs(rConfigs[ci].R - denR) < 0.01 && !hMCTden)
              hMCTden = LoadMCTruthXsec(rConfigs[ci]);
          }
        }
        if (hMCTnum && hMCTden) {
          TH1* hMCTratio = (TH1*)hMCTnum->Clone(Form("combMCT_%d_%d", ip, GetNN()));
          hMCTratio->Divide(hMCTnum, hMCTden, 1., 1., "");
          TH1* hDR = (TH1*)hMCTratio->Clone(Form("combMCTdr_%d_%d", ip, GetNN()));
          hDR->Divide(hMCTratio, ratioHistDataList[ip], 1., 1., "");
          hDR->SetStats(0);
          hDR->SetLineColor(kGray+2); hDR->SetLineStyle(7); hDR->SetLineWidth(2); hDR->SetMarkerStyle(0);
          hDR->Draw("l same");
          delete hMCTratio;
        }
        if (mode != kUE) { delete hMCTnum; delete hMCTden; }
      }

      // Standalone model double ratios (use model-specific colors)
      for (size_t im = 0; im < crModels.size(); ++im) {
        if (ratioHistModelLists.count(crModels[im].name.Data()) &&
            ip < (int)ratioHistModelLists[crModels[im].name.Data()].size() &&
            ratioHistModelLists[crModels[im].name.Data()][ip]) {
          TH1* hMR = ratioHistModelLists[crModels[im].name.Data()][ip];
          TH1* hDR = (TH1*)hMR->Clone(Form("combDR_%s_%d_%d", crModels[im].name.Data(), ip, GetNN()));
          hDR->Divide(hMR, ratioHistDataList[ip], 1., 1., "");
          hDR->SetStats(0);
          hDR->SetLineColor(crModels[im].lineColor);
          hDR->SetLineStyle(crModels[im].lineStyle); hDR->SetLineWidth(crModels[im].lineWidth);
          hDR->SetMarkerStyle(0);
          hDR->Draw("l same");
        }
      }

      // POWHEG double ratio
      if (ip < (int)ratioHistPOWHEGList.size() && ratioHistPOWHEGList[ip]) {
        TH1* hDR = (TH1*)ratioHistPOWHEGList[ip]->Clone(Form("combDRpow_%d_%d", ip, GetNN()));
        hDR->Divide(ratioHistPOWHEGList[ip], ratioHistDataList[ip], 1., 1., "");
        hDR->SetStats(0);
        hDR->SetLineColor(kPowhegColor);
        hDR->SetLineStyle(3); hDR->SetLineWidth(3); hDR->SetMarkerStyle(0);
        hDR->Draw("l same");
      }

      // Legend in first panel
      if (ip == 0) {
        TLegend* combLeg = new TLegend(0.22, 0.45, 0.65, 0.88, NULL, "brNDC");
        combLeg->SetTextSize(0.045); combLeg->SetBorderSize(0); combLeg->SetFillStyle(0);
        TH1F* dmct = new TH1F(Form("dmct_c_%d", GetNN()), "", 1, 0, 1);
        dmct->SetLineColor(kGray+2); dmct->SetLineStyle(7); dmct->SetLineWidth(2);
        combLeg->AddEntry(dmct, "MC part-level", "l");
        for (auto& m : crModels) {
          TH1F* dm = new TH1F(Form("dm_%s_c_%d", m.name.Data(), GetNN()), "", 1, 0, 1);
          dm->SetLineColor(m.lineColor); dm->SetLineStyle(m.lineStyle); dm->SetLineWidth(m.lineWidth);
          combLeg->AddEntry(dm, m.legLabel.Data(), "l");
        }
        TH1F* dpw = new TH1F(Form("dpw_c_%d", GetNN()), "", 1, 0, 1);
        dpw->SetLineColor(kPowhegColor); dpw->SetLineStyle(3); dpw->SetLineWidth(3);
        combLeg->AddEntry(dpw, "POWHEG NLO", "l");
        combLeg->Draw();
      }
    }

    if (GetDRAWPLOTS()) {
      combinedCanvas->Update(); combinedCanvas->Modified();
      combinedCanvas->Print(Form("%s/CustomRatios_Combined.pdf", outputDir.Data()));
    }
  }
}

// ============================================
// Helper Functions for Cross Section Comparison
// ============================================

// Run 2 UE-subtracted reference (HEPData Figure 3)
struct Run2UEXsecPoint {
  double pt;
  double sigma;
  double statErr;
  double sysErr;
};

// Hard-coded Run 2 UE-subtracted cross sections
// Source: HEPData-ins2026265-v1-Figure_3.csv (Yongzhen Hou, pp 13 TeV)
inline std::vector<Run2UEXsecPoint> GetRun2UEXsecData(double R) {
  std::vector<Run2UEXsecPoint> points;

  // R = 0.2
  if (TMath::Abs(R - 0.2) < 0.01) {
    points = {
      {5.5, 0.39279, 7.6253e-05, 0.025213}, {6.5, 0.19321, 4.7465e-05, 0.012373},
      {7.5, 0.10585, 3.1971e-05, 0.006898}, {8.5, 0.063039, 2.245e-05, 0.0040825},
      {9.5, 0.040221, 1.5825e-05, 0.0026695},{11.0, 0.022846, 1.0016e-05, 0.0014934},
      {13.0, 0.011834, 5.9707e-06, 0.00080998},{15.0, 0.0068225, 4.0547e-06, 0.00045105},
      {17.0, 0.0041938, 2.8491e-06, 0.00030587},{19.0, 0.0027479, 2.0939e-06, 0.00020752},
      {22.5, 0.0014816, 1.2336e-06, 0.00011434},{27.5, 0.00065555, 6.6925e-07, 5.136e-05},
      {35.0, 0.00025239, 3.2027e-07, 2.037e-05},{45.0, 8.2464e-05, 1.277e-07, 7.0804e-06},
      {55.0, 3.2965e-05, 6.0345e-08, 3.0121e-06},{65.0, 1.5071e-05, 3.1532e-08, 1.4066e-06},
      {77.5, 6.6585e-06, 1.5427e-08, 6.2289e-07},{92.5, 2.7858e-06, 6.9443e-09, 2.5918e-07},
      {120.0, 8.3394e-07, 2.1792e-09, 7.9439e-08}
    };
  }
  // R = 0.3
  else if (TMath::Abs(R - 0.3) < 0.01) {
    points = {
      {5.5, 0.57284, 0.00010932, 0.048655}, {6.5, 0.29811, 7.0364e-05, 0.025192},
      {7.5, 0.16905, 4.8196e-05, 0.013967}, {8.5, 0.10297, 3.4136e-05, 0.010109},
      {9.5, 0.066589, 2.3955e-05, 0.0061598},{11.0, 0.037983, 1.4991e-05, 0.0033773},
      {13.0, 0.019481, 8.8225e-06, 0.0017249},{15.0, 0.01105, 5.9373e-06, 0.00097616},
      {17.0, 0.0066902, 4.1255e-06, 0.00060658},{19.0, 0.0043035, 2.9769e-06, 0.00039429},
      {22.5, 0.0022462, 1.6924e-06, 0.00020498},{27.5, 0.00095148, 8.839e-07, 8.8885e-05},
      {35.0, 0.00034941, 4.0342e-07, 3.2996e-05},{45.0, 0.00010859, 1.5224e-07, 1.0481e-05},
      {55.0, 4.2118e-05, 6.9303e-08, 4.1755e-06},{65.0, 1.8936e-05, 3.5354e-08, 1.8735e-06},
      {77.5, 8.2351e-06, 1.6917e-08, 8.588e-07},{92.5, 3.3971e-06, 7.4713e-09, 3.5434e-07},
      {120.0, 1.0038e-06, 2.3065e-09, 1.0718e-07}
    };
  }
  // R = 0.4
  else if (TMath::Abs(R - 0.4) < 0.01) {
    points = {
      {5.5, 0.73271, 0.00014903, 0.062585}, {6.5, 0.39868, 9.8507e-05, 0.035592},
      {7.5, 0.2329, 6.8456e-05, 0.021378},  {8.5, 0.14445, 4.91e-05, 0.013835},
      {9.5, 0.094007, 3.4181e-05, 0.0094001},{11.0, 0.053479, 2.0988e-05, 0.0053987},
      {13.0, 0.026927, 1.2121e-05, 0.0027346},{15.0, 0.014953, 8.0655e-06, 0.0014948},
      {17.0, 0.0089322, 5.5281e-06, 0.00089048},{19.0, 0.0056478, 3.8943e-06, 0.00057549},
      {22.5, 0.0028768, 2.1461e-06, 0.00028901},{27.5, 0.0011753, 1.079e-06, 0.00011904},
      {35.0, 0.00041953, 4.7569e-07, 4.268e-05},{45.0, 0.00012699, 1.7316e-07, 1.3475e-05},
      {55.0, 4.8626e-05, 7.7014e-08, 5.2298e-06},{65.0, 2.1715e-05, 3.8654e-08, 2.3392e-06},
      {77.5, 9.3914e-06, 1.825e-08, 1.049e-06},{92.5, 3.8434e-06, 7.9487e-09, 4.2323e-07},
      {120.0, 1.1284e-06, 2.4285e-09, 1.2455e-07}
    };
  }
  // R = 0.5
  else if (TMath::Abs(R - 0.5) < 0.01) {
    points = {
      {5.5, 0.86828, 0.00020043, 0.07479}, {6.5, 0.48908, 0.0001365, 0.042971},
      {7.5, 0.29388, 9.6539e-05, 0.025756}, {8.5, 0.18504, 7.0388e-05, 0.018885},
      {9.5, 0.12152, 4.8926e-05, 0.012162},{11.0, 0.069247, 2.946e-05, 0.0070507},
      {13.0, 0.034442, 1.6789e-05, 0.0035605},{15.0, 0.01877, 1.1101e-05, 0.0019294},
      {17.0, 0.011047, 7.4724e-06, 0.0011449},{19.0, 0.0069536, 5.1588e-06, 0.00072503},
      {22.5, 0.0034532, 2.7229e-06, 0.00037111},{27.5, 0.001365, 1.3196e-06, 0.00015309},
      {35.0, 0.00047571, 5.6588e-07, 5.4436e-05},{45.0, 0.0001409, 2.0032e-07, 1.66e-05},
      {55.0, 5.3565e-05, 8.7824e-08, 6.4327e-06},{65.0, 2.3801e-05, 4.3555e-08, 2.8644e-06},
      {77.5, 1.026e-05, 2.0376e-08, 1.2729e-06},{92.5, 4.1789e-06, 8.7922e-09, 5.1604e-07},
      {120.0, 1.2205e-06, 2.6639e-09, 1.5134e-07}
    };
  }
  // R = 0.6
  else if (TMath::Abs(R - 0.6) < 0.01) {
    points = {
      {5.5, 0.97671, 0.00026394, 0.082684}, {6.5, 0.56118, 0.0001862, 0.048606},
      {7.5, 0.34618, 0.00013431, 0.029971}, {8.5, 0.22249, 0.00010011, 0.021061},
      {9.5, 0.14741, 7.0049e-05, 0.014369},{11.0, 0.084417, 4.1758e-05, 0.0083738},
      {13.0, 0.0417, 2.3739e-05, 0.004234},{15.0, 0.022456, 1.5725e-05, 0.0024118},
      {17.0, 0.013059, 1.0411e-05, 0.001341},{19.0, 0.0081984, 7.021e-06, 0.00088122},
      {22.5, 0.0040321, 3.5788e-06, 0.00044384},{27.5, 0.0015348, 1.6647e-06, 0.00017747},
      {35.0, 0.0005211, 6.9559e-07, 6.1686e-05},{45.0, 0.00015077, 2.3989e-07, 1.8672e-05},
      {55.0, 5.687e-05, 1.0386e-07, 7.3169e-06},{65.0, 2.5243e-05, 5.1188e-08, 3.2319e-06},
      {77.5, 1.0845e-05, 2.3756e-08, 1.4131e-06},{92.5, 4.3994e-06, 1.0173e-08, 5.7709e-07},
      {120.0, 1.2777e-06, 3.0575e-09, 1.6651e-07}
    };
  }
  // R = 0.7
  else if (TMath::Abs(R - 0.7) < 0.01) {
    points = {
      {5.5, 1.0815, 0.00035038, 0.082193}, {6.5, 0.62221, 0.00025633, 0.055258},
      {7.5, 0.39121, 0.00018964, 0.033611}, {8.5, 0.25658, 0.00014485, 0.022893},
      {9.5, 0.17185, 0.00010275, 0.01618},  {11.0, 0.099064, 6.12e-05, 0.0099809},
      {13.0, 0.048885, 3.5214e-05, 0.004959},{15.0, 0.026002, 2.3401e-05, 0.0027029},
      {17.0, 0.015078, 1.5457e-05, 0.0016587},{19.0, 0.0094971, 1.023e-05, 0.0010626},
      {22.5, 0.004699, 5.4377e-06, 0.00052095},{27.5, 0.0017396, 2.4159e-06, 0.00020834},
      {35.0, 0.00057631, 1.0265e-06, 6.9532e-05},{45.0, 0.00016294, 3.6086e-07, 2.0523e-05},
      {55.0, 6.1234e-05, 1.607e-07, 8.1226e-06},{65.0, 2.7226e-05, 8.117e-08, 3.7395e-06},
      {77.5, 1.1687e-05, 3.824e-08, 1.5777e-06},{92.5, 4.7354e-06, 1.6527e-08, 6.3909e-07},
      {120.0, 1.3679e-06, 4.9716e-09, 1.7564e-07}
    };
  }

  return points;
}

// Create Run2 UE-sub TGraphErrors (stat) and TGraphAsymmErrors (sys)
// with x-errors matched to Run 3 binning
inline std::pair<TGraphErrors*, TGraphAsymmErrors*>
CreateRun2UEGraphs(double R) {
  std::vector<Run2UEXsecPoint> points = GetRun2UEXsecData(R);
  if (points.empty()) return {nullptr, nullptr};

  const Double_t* ptbinGen = GetPtbinGen();
  int nBinsGen = GetNptBinsGen();

  int nPoints = points.size();
  TGraphErrors* grStat = new TGraphErrors(nPoints);
  TGraphAsymmErrors* grSys = new TGraphAsymmErrors(nPoints);

  grStat->SetName(Form("Run2UE_Stat_R%.1f", R));
  grSys->SetName(Form("Run2UE_Sys_R%.1f", R));

  // Run2 styling: black dashed line (style 2), gray box, open square marker
  grStat->SetLineColor(kBlack);
  grStat->SetMarkerColor(kBlack);
  grStat->SetMarkerStyle(25);  // Open square
  grStat->SetLineStyle(2);     // Dashed
  grStat->SetLineWidth(2);

  grSys->SetFillColorAlpha(kGray + 1, 0.5);
  grSys->SetFillStyle(1001);
  grSys->SetLineWidth(0);
  grSys->SetLineColor(0);

  for (int i = 0; i < nPoints; ++i) {
    const auto& pt = points[i];

    // Match Run 3 bin center
    int binIdx = -1;
    for (int j = 0; j < nBinsGen; ++j) {
      double binCenter = 0.5 * (ptbinGen[j] + ptbinGen[j+1]);
      if (TMath::Abs(pt.pt - binCenter) < 0.11) {
        binIdx = j;
        break;
      }
    }

    double exLow = 0.0, exHigh = 0.0;
    if (binIdx >= 0 && binIdx < nBinsGen) {
      double binLow = ptbinGen[binIdx];
      double binHigh = ptbinGen[binIdx+1];
      exLow  = pt.pt - binLow;
      exHigh = binHigh - pt.pt;
    } else {
      // Fallback: symmetric errors from neighbor spacing
      if (i > 0 && i < nPoints - 1) {
        exLow  = 0.5 * (pt.pt - points[i-1].pt);
        exHigh = 0.5 * (points[i+1].pt - pt.pt);
      } else if (i == 0 && nPoints > 1) {
        exLow = exHigh = 0.5 * (points[i+1].pt - pt.pt);
      } else if (i == nPoints - 1 && nPoints > 1) {
        exLow = exHigh = 0.5 * (pt.pt - points[i-1].pt);
      }
    }

    grStat->SetPoint(i, pt.pt, pt.sigma);
    grStat->SetPointError(i, 0.5 * (exLow + exHigh), pt.statErr);

    grSys->SetPoint(i, pt.pt, pt.sigma);
    grSys->SetPointError(i, exLow, exHigh, pt.sysErr, pt.sysErr);
  }

  return {grStat, grSys};
}

// Helper: create ratio histogram (Run3 UE-sub / Run2 UE-sub)
inline TH1* CreateRatioTH1vsTGraphUE(TH1* hNum, TGraphErrors* grDenom,
                                     const char* name) {
  if (!hNum || !grDenom) return nullptr;

  TH1* hRatio = (TH1*)hNum->Clone(name);
  hRatio->Reset();
  hRatio->SetDirectory(0);
  hRatio->SetStats(0);

  int nPoints = grDenom->GetN();
  double* xDen = grDenom->GetX();
  double* yDen = grDenom->GetY();
  double* eDen = grDenom->GetEY();

  for (int i = 1; i <= hNum->GetNbinsX(); ++i) {
    double xLow = hNum->GetBinLowEdge(i);
    double xUp  = xLow + hNum->GetBinWidth(i);
    double yDenVal = 0.0;
    double eDenVal = 0.0;
    bool found = false;

    for (int j = 0; j < nPoints; ++j) {
      if (xDen[j] >= xLow && xDen[j] < xUp) {
        yDenVal = yDen[j];
        eDenVal = eDen[j];
        found = true;
        break;
      }
    }

    // Diagnostic: print for high-pT bins (potential ratio issues)
    if (xLow >= 85) {
      double yNum = hNum->GetBinContent(i);
      std::cerr << "[RatioUE-diag] " << name << " bin[" << xLow << "," << xUp
                << "] found=" << found << " yDen=" << yDenVal << " yNum=" << yNum << std::endl;
    }

    if (!found || yDenVal <= 0) continue;

    double yNum = hNum->GetBinContent(i);
    double eNum = hNum->GetBinError(i);
    if (yNum <= 0) continue;

    double ratio = yNum / yDenVal;
    double relNum = eNum / yNum;
    double relDen = eDenVal / yDenVal;
    double eRatio = ratio * TMath::Sqrt(relNum * relNum + relDen * relDen);

    hRatio->SetBinContent(i, ratio);
    hRatio->SetBinError(i, eRatio);
  }

  return hRatio;
}

// Load PYTHIA8 13.6 TeV histogram for given R
// Get Run2 cross section data for given R (hard-coded from CSV)
struct Run2XsecPoint {
  double pt, sigma, statErr, sysErr;
};

std::vector<Run2XsecPoint> GetRun2XsecData(double R) {
  std::vector<Run2XsecPoint> points;
  
  // R=0.2
  if (TMath::Abs(R - 0.2) < 0.01) {
    points = {
      {5.5, 0.41844, 7.5615e-05, 0.025398}, {6.5, 0.20398, 4.7258e-05, 0.012547},
      {7.5, 0.11078, 3.1837e-05, 0.0070092}, {8.5, 0.065436, 2.2311e-05, 0.0041402},
      {9.5, 0.041482, 1.5753e-05, 0.0027067}, {11.0, 0.023514, 9.9869e-06, 0.0015188},
      {13.0, 0.012152, 5.9537e-06, 0.00082559}, {15.0, 0.0069859, 4.0296e-06, 0.00045924},
      {17.0, 0.0042814, 2.8277e-06, 0.00031121}, {19.0, 0.002798, 2.0827e-06, 0.00021082},
      {22.5, 0.0015066, 1.2322e-06, 0.00011611}, {27.5, 0.00066571, 6.683e-07, 5.2124e-05},
      {35.0, 0.00025605, 3.1841e-07, 2.0669e-05}, {45.0, 8.3592e-05, 1.2624e-07, 7.1843e-06},
      {55.0, 3.3397e-05, 5.9352e-08, 3.0558e-06}, {65.0, 1.5262e-05, 3.0882e-08, 1.427e-06},
      {77.5, 6.7409e-06, 1.506e-08, 6.3208e-07}, {92.5, 2.8197e-06, 6.7636e-09, 2.6313e-07},
      {120.0, 8.4402e-07, 2.1194e-09, 8.0712e-08}
    };
  }
  // R=0.3
  else if (TMath::Abs(R - 0.3) < 0.01) {
    points = {
      {5.5, 0.70636, 0.00011433, 0.047252}, {6.5, 0.36018, 7.3718e-05, 0.025187},
      {7.5, 0.19998, 5.0381e-05, 0.014265}, {8.5, 0.11935, 3.549e-05, 0.0087693},
      {9.5, 0.075733, 2.4876e-05, 0.0056478}, {11.0, 0.042625, 1.5549e-05, 0.003231},
      {13.0, 0.02159, 9.1255e-06, 0.0016828}, {15.0, 0.012099, 6.1023e-06, 0.00096227},
      {17.0, 0.0072417, 4.2244e-06, 0.00060607}, {19.0, 0.0046109, 3.0462e-06, 0.00039733},
      {22.5, 0.0023891, 1.7317e-06, 0.00020806}, {27.5, 0.0010047, 9.0217e-07, 9.0813e-05},
      {35.0, 0.00036671, 4.1001e-07, 3.3842e-05}, {45.0, 0.00011338, 1.5408e-07, 1.0764e-05},
      {55.0, 4.3788e-05, 6.9903e-08, 4.2905e-06}, {65.0, 1.962e-05, 3.5568e-08, 1.9241e-06},
      {77.5, 8.51e-06, 1.6988e-08, 8.8232e-07}, {92.5, 3.5038e-06, 7.4925e-09, 3.6405e-07},
      {120.0, 1.0341e-06, 2.3111e-09, 1.103e-07}
    };
  }
  // R=0.4
  else if (TMath::Abs(R - 0.4) < 0.01) {
    points = {
      {5.5, 1.0972, 0.00016396, 0.074335}, {6.5, 0.58418, 0.00010887, 0.042393},
      {7.5, 0.33224, 7.5637e-05, 0.024799}, {8.5, 0.20001, 5.3822e-05, 0.015407},
      {9.5, 0.12614, 3.7416e-05, 0.010138}, {11.0, 0.069825, 2.2965e-05, 0.0057313},
      {13.0, 0.034225, 1.3182e-05, 0.0029103}, {15.0, 0.018498, 8.6563e-06, 0.001585},
      {17.0, 0.010764, 5.8873e-06, 0.00094719}, {19.0, 0.0066421, 4.1364e-06, 0.00061878},
      {22.5, 0.0033177, 2.2742e-06, 0.00031195}, {27.5, 0.0013313, 1.136e-06, 0.00012884},
      {35.0, 0.00046794, 4.9698e-07, 4.6138e-05}, {45.0, 0.00013985, 1.7976e-07, 1.4546e-05},
      {55.0, 5.2994e-05, 7.9569e-08, 5.6228e-06}, {65.0, 2.3472e-05, 3.9801e-08, 2.5033e-06},
      {77.5, 1.0087e-05, 1.8747e-08, 1.1183e-06}, {92.5, 4.1091e-06, 8.1518e-09, 4.4943e-07},
      {120.0, 1.2028e-06, 2.4881e-09, 1.3178e-07}
    };
  }
  // R=0.5
  else if (TMath::Abs(R - 0.5) < 0.01) {
    points = {
      {5.5, 1.57, 0.00023005, 0.10651}, {6.5, 0.88103, 0.0001582, 0.064779},
      {7.5, 0.51888, 0.00011259, 0.040347}, {8.5, 0.31783, 8.148e-05, 0.02555},
      {9.5, 0.20177, 5.6703e-05, 0.01709}, {11.0, 0.11108, 3.4358e-05, 0.0099431},
      {13.0, 0.053294, 1.9453e-05, 0.004869}, {15.0, 0.027861, 1.2594e-05, 0.0025515},
      {17.0, 0.015669, 8.3744e-06, 0.0014813}, {19.0, 0.0094169, 5.757e-06, 0.00091899},
      {22.5, 0.0044912, 3.0312e-06, 0.00046255}, {27.5, 0.0017143, 1.4494e-06, 0.00018721},
      {35.0, 0.00058007, 6.1253e-07, 6.5211e-05}, {45.0, 0.00016772, 2.1438e-07, 1.9522e-05},
      {55.0, 6.2547e-05, 9.323e-08, 7.4437e-06}, {65.0, 2.7375e-05, 4.5977e-08, 3.2697e-06},
      {77.5, 1.1666e-05, 2.1427e-08, 1.4378e-06}, {92.5, 4.7123e-06, 9.2221e-09, 5.779e-07},
      {120.0, 1.3689e-06, 2.7899e-09, 1.6831e-07}
    };
  }
  // R=0.6
  else if (TMath::Abs(R - 0.6) < 0.01) {
    points = {
      {5.5, 2.0341, 0.00031379, 0.13237}, {6.5, 1.2093, 0.00022468, 0.086216},
      {7.5, 0.74324, 0.00016445, 0.056814}, {8.5, 0.47063, 0.00012216, 0.038551},
      {9.5, 0.30534, 8.5923e-05, 0.027247}, {11.0, 0.17086, 5.2164e-05, 0.01576},
      {13.0, 0.08204, 2.9589e-05, 0.0077182}, {15.0, 0.042196, 1.9148e-05, 0.0042407},
      {17.0, 0.023068, 1.2499e-05, 0.0022311}, {19.0, 0.013484, 8.3817e-06, 0.0013943},
      {22.5, 0.0061779, 4.255e-06, 0.00066363}, {27.5, 0.0022162, 1.9323e-06, 0.00025252},
      {35.0, 0.00071612, 7.8683e-07, 8.3948e-05}, {45.0, 0.00019915, 2.6612e-07, 2.4497e-05},
      {55.0, 7.282e-05, 1.1368e-07, 9.3194e-06}, {65.0, 3.1557e-05, 5.5526e-08, 4.0201e-06},
      {77.5, 1.3315e-05, 2.5612e-08, 1.7263e-06}, {92.5, 5.3326e-06, 1.0924e-08, 6.9561e-07},
      {120.0, 1.536e-06, 3.2751e-09, 1.9806e-07}
    };
  }
  // R=0.7
  else if (TMath::Abs(R - 0.7) < 0.01) {
    points = {
      {5.5, 2.4314, 0.00043375, 0.16635}, {6.5, 1.53, 0.00032329, 0.11086},
      {7.5, 0.98795, 0.00024463, 0.074637}, {8.5, 0.65147, 0.00018752, 0.052778},
      {9.5, 0.43637, 0.00013418, 0.039069}, {11.0, 0.25245, 8.2418e-05, 0.024531},
      {13.0, 0.12473, 4.7667e-05, 0.012127}, {15.0, 0.064384, 3.1206e-05, 0.0064126},
      {17.0, 0.034976, 2.0356e-05, 0.0037265}, {19.0, 0.020087, 1.3397e-05, 0.0021984},
      {22.5, 0.0088944, 6.5971e-06, 0.0009716}, {27.5, 0.0029685, 2.832e-06, 0.00035246},
      {35.0, 0.00090089, 1.1035e-06, 0.000108}, {45.0, 0.00023787, 3.5907e-07, 2.9811e-05},
      {55.0, 8.485e-05, 1.5086e-07, 1.1207e-05}, {65.0, 3.6275e-05, 7.3054e-08, 4.9618e-06},
      {77.5, 1.5128e-05, 3.3403e-08, 2.0325e-06}, {92.5, 6.0063e-06, 1.4148e-08, 8.061e-07},
      {120.0, 1.7127e-06, 4.2033e-09, 2.1739e-07}
    };
  }
  
  return points;
}

// Create Run2 TGraphErrors (stat) and TGraphAsymmErrors (sys) from data points
std::pair<TGraphErrors*, TGraphAsymmErrors*> CreateRun2Graphs(double R) {
  std::vector<Run2XsecPoint> points = GetRun2XsecData(R);
  if (points.empty()) {
    return {nullptr, nullptr};
  }
  
  // Get Run3 binning to match bin edges
  const Double_t* ptbinGen = GetPtbinGen();
  int nBinsGen = GetNptBinsGen();
  
  int nPoints = points.size();
  TGraphErrors* grStat = new TGraphErrors(nPoints);
  TGraphAsymmErrors* grSys = new TGraphAsymmErrors(nPoints);
  
  grStat->SetName(Form("Run2_Stat_R%.1f", R));
  grSys->SetName(Form("Run2_Sys_R%.1f", R));
  
  // Run2 styling: black points with stat errors, gray box for sys
  grStat->SetLineColor(kBlack);
  grStat->SetMarkerColor(kBlack);
  grStat->SetMarkerStyle(20);
  grStat->SetLineWidth(2);
  
  grSys->SetFillColorAlpha(kGray+1, 0.5);
  grSys->SetFillStyle(1001);
  grSys->SetLineWidth(0);
  grSys->SetLineColor(0);
  
  for (int i = 0; i < nPoints; ++i) {
    const auto& pt = points[i];
    
    // Find which Run3 bin this point belongs to (pt.pt is bin center)
    int binIdx = -1;
    for (int j = 0; j < nBinsGen; ++j) {
      double binCenter = (ptbinGen[j] + ptbinGen[j+1]) / 2.0;
      if (TMath::Abs(pt.pt - binCenter) < 0.1) {  // Match bin center within 0.1 GeV
        binIdx = j;
        break;
      }
    }
    
    // Calculate x-axis error from Run3 bin edges
    double exLow = 0, exHigh = 0;
    if (binIdx >= 0 && binIdx < nBinsGen) {
      double binLow = ptbinGen[binIdx];
      double binHigh = ptbinGen[binIdx + 1];
      exLow = pt.pt - binLow;
      exHigh = binHigh - pt.pt;
    } else {
      // Fallback: use half distance to neighbors
      if (i > 0 && i < nPoints - 1) {
        exLow = (pt.pt - points[i-1].pt) / 2.0;
        exHigh = (points[i+1].pt - pt.pt) / 2.0;
      } else if (i == 0) {
        exLow = (points[i+1].pt - pt.pt) / 2.0;
        exHigh = exLow;
      } else {
        exLow = (pt.pt - points[i-1].pt) / 2.0;
        exHigh = exLow;
      }
    }
    
    grStat->SetPoint(i, pt.pt, pt.sigma);
    grStat->SetPointError(i, (exLow + exHigh) / 2.0, pt.statErr);
    
    grSys->SetPoint(i, pt.pt, pt.sigma);
    grSys->SetPointError(i, exLow, exHigh, pt.sysErr, pt.sysErr);
  }
  
  return {grStat, grSys};
}

// Helper function to create ratio histogram from TH1 and TGraphErrors
// Only uses overlapping pT range and matches binning
TH1* CreateRatioTH1vsTGraph(TH1* hNum, TGraphErrors* grDenom, const char* name) {
  if (!hNum || !grDenom) return nullptr;
  
  // Create ratio histogram with same binning as numerator
  TH1* hRatio = (TH1*)hNum->Clone(name);
  hRatio->Reset();
  hRatio->SetDirectory(0);
  hRatio->SetStats(0);
  
  int nPoints = grDenom->GetN();
  double* xDenom = grDenom->GetX();
  double* yDenom = grDenom->GetY();
  double* eyDenom = grDenom->GetEY();
  
  for (int i = 1; i <= hNum->GetNbinsX(); ++i) {
    double xCenter = hNum->GetBinCenter(i);
    double xLow = hNum->GetBinLowEdge(i);
    double xUp = hNum->GetBinLowEdge(i) + hNum->GetBinWidth(i);
    
    // Find matching point in TGraph (within bin range)
    double yDenomVal = 0;
    double eyDenomVal = 0;
    bool found = false;
    
    for (int j = 0; j < nPoints; ++j) {
      if (xDenom[j] >= xLow && xDenom[j] < xUp) {
        yDenomVal = yDenom[j];
        eyDenomVal = eyDenom[j];
        found = true;
        break;
      }
    }
    
    if (found && yDenomVal > 0) {
      double yNum = hNum->GetBinContent(i);
      double eyNum = hNum->GetBinError(i);
      
      if (yNum > 0) {
        double ratio = yNum / yDenomVal;
        double relErrNum = eyNum / yNum;
        double relErrDenom = eyDenomVal / yDenomVal;
        double ratioErr = ratio * TMath::Sqrt(relErrNum * relErrNum + relErrDenom * relErrDenom);
        
        hRatio->SetBinContent(i, ratio);
        hRatio->SetBinError(i, ratioErr);
      }
    }
  }
  
  return hRatio;
}

// ============================================
// Run 2 R-ratio from cross-section data
// ============================================
// Computes sigma(R)/sigma(Rref) from the hardcoded Run 2 cross sections.
// The ratio values are exact; error bars are approximate (uncorrelated propagation).

TGraphErrors* CreateRun2RRatioGraph(double R, double Rref, bool isUE) {
  // Skip self-ratio
  if (TMath::Abs(R - Rref) < 0.01) return nullptr;

  // Get cross-section data for numerator and denominator
  auto getXsec = [&](double r) {
    if (isUE) {
      std::vector<Run2UEXsecPoint> pts = GetRun2UEXsecData(r);
      // Convert to common format {pt, sigma, statErr}
      std::vector<std::tuple<double, double, double>> out;
      for (const auto& p : pts) out.push_back({p.pt, p.sigma, p.statErr});
      return out;
    } else {
      std::vector<Run2XsecPoint> pts = GetRun2XsecData(r);
      std::vector<std::tuple<double, double, double>> out;
      for (const auto& p : pts) out.push_back({p.pt, p.sigma, p.statErr});
      return out;
    }
  };

  auto numData = getXsec(R);
  auto denData = getXsec(Rref);
  if (numData.empty() || denData.empty()) return nullptr;

  // Build map from pT → {sigma, err} for denominator
  std::map<double, std::pair<double, double>> denMap;
  for (const auto& d : denData)
    denMap[std::get<0>(d)] = {std::get<1>(d), std::get<2>(d)};

  // Compute ratio at matching pT points
  std::vector<double> xVals, yVals, yErrs;
  for (const auto& n : numData) {
    double pt = std::get<0>(n);
    auto it = denMap.find(pt);
    if (it == denMap.end()) continue;
    double num = std::get<1>(n), numE = std::get<2>(n);
    double den = it->second.first, denE = it->second.second;
    if (den <= 0 || num <= 0) continue;

    double ratio = num / den;
    double relErr = TMath::Sqrt((numE / num) * (numE / num) + (denE / den) * (denE / den));
    xVals.push_back(pt);
    yVals.push_back(ratio);
    yErrs.push_back(ratio * relErr);
  }

  if (xVals.empty()) return nullptr;

  int nPts = xVals.size();
  TGraphErrors* gr = new TGraphErrors(nPts);
  gr->SetName(Form("Run2_RRatio_R%.1f_over_R%.1f%s", R, Rref, isUE ? "_UE" : ""));
  for (int i = 0; i < nPts; ++i) {
    gr->SetPoint(i, xVals[i], yVals[i]);
    gr->SetPointError(i, 0, yErrs[i]);
  }

  // Styling: open marker, dashed line
  gr->SetMarkerStyle(24);  // open circle
  gr->SetMarkerSize(0.8);
  gr->SetLineStyle(2);     // dashed
  gr->SetLineWidth(2);
  return gr;
}

// ============================================
// R-Ratio Plots: sigma(R) / sigma(R_ref)
// ============================================
// Works for both NonUE (vector-based) and UE (map-based) inputs.
// Computes bin-by-bin ratio with error propagation (uncorrelated).

void DrawRRatioPlots(const std::map<double, TH1*>& spectraByR,
                     const std::vector<Color_t>& colors,
                     const TString& outputDir,
                     const TString& tag = "")
{
  // Find reference spectrum
  TH1* hRef = nullptr;
  for (const auto& kv : spectraByR) {
    if (TMath::Abs(kv.first - kRRatioRef) < 0.01) {
      hRef = kv.second;
      break;
    }
  }
  if (!hRef) {
    std::cerr << "[Warning] R-ratio: reference R=" << kRRatioRef
              << " not found. Skipping R-ratio plot." << std::endl;
    return;
  }

  bool isUE = tag.Contains("UE");

  TString canvasName = tag.IsNull() ? "RRatio" : Form("RRatio_%s", tag.Data());
  Filipad2* pad = new Filipad2(canvasName.Data(), ++GetNN(), 2.2, 0.0, 0.0, 100, 50, 0.7, 5, 1);
  pad->Draw();
  TPad* mainPad = pad->GetPad(1);
  optFili(*mainPad, 0, 0, 0, 0);
  mainPad->cd();
  gPad->SetGrid(0);

  TLegend* leg = new TLegend(0.50, 0.55, 0.92, 0.95, NULL, "brNDC");
  leg->SetTextSize(0.032);
  leg->SetBorderSize(0);
  leg->SetNColumns(2);
  leg->SetHeader(Form("#sigma(#it{R}) / #sigma(#it{R} = %.1f)", kRRatioRef));

  bool firstDraw = true;
  bool run2LegendAdded = false;

  for (const auto& kv : spectraByR) {
    const double R = kv.first;
    TH1* hSpec = kv.second;
    if (!hSpec) continue;
    if (TMath::Abs(R - kRRatioRef) < 0.01) continue;  // skip self-ratio

    TH1* hRatio = (TH1*)hSpec->Clone(Form("RRatio_R%.1f_over_R%.1f%s", R, kRRatioRef, tag.Data()));
    hRatio->Reset();

    int nBins = hRatio->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
      double num = hSpec->GetBinContent(i);
      double den = hRef->GetBinContent(i);
      double numErr = hSpec->GetBinError(i);
      double denErr = hRef->GetBinError(i);

      if (den > 0 && num > 0) {
        double ratio = num / den;
        double relErr = TMath::Sqrt((numErr / num) * (numErr / num) + (denErr / den) * (denErr / den));
        hRatio->SetBinContent(i, ratio);
        hRatio->SetBinError(i, ratio * relErr);
      }
    }

    Color_t col = GetColorForR(R);
    hRatio->SetLineColor(col);
    hRatio->SetMarkerColor(col);
    hRatio->SetMarkerStyle(20);
    hRatio->SetLineWidth(2);

    if (firstDraw) {
      hset(*hRatio, GetJetPtGenTitleX(),
           Form("#sigma(#it{R}) / #sigma(#it{R} = %.1f)", kRRatioRef),
           1.0, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
      hRatio->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
      hRatio->GetYaxis()->SetRangeUser(0.0, 2.5);
      hRatio->Draw("pe");
      firstDraw = false;
    } else {
      hRatio->Draw("pe same");
    }

    leg->AddEntry(hRatio, Form("Run 3, #it{R} = %.1f", R), "pe");

    // Overlay Run 2 R-ratio computed from cross-section data
    TGraphErrors* grRun2 = CreateRun2RRatioGraph(R, kRRatioRef, isUE);
    if (grRun2) {
      grRun2->SetMarkerColor(col);
      grRun2->SetLineColor(col);
      grRun2->Draw("pz same");
      if (!run2LegendAdded) {
        // Add a single "Run 2" entry using first graph as template
        leg->AddEntry(grRun2, "Run 2 (open)", "pe");
        run2LegendAdded = true;
      }
    }

  }

  if (firstDraw) {
    std::cerr << "[Warning] R-ratio: no spectra to plot (only reference?)." << std::endl;
    return;
  }

  // Unity line
  TLine* lineUnity = new TLine(GetPlotPtMin(), 1.0, GetPlotPtMax(), 1.0);
  lineUnity->SetLineColor(kBlack);
  lineUnity->SetLineStyle(2);
  lineUnity->SetLineWidth(1);
  lineUnity->Draw("lsame");

  leg->Draw();

  if (GetDRAWPLOTS()) {
    TString fname = tag.IsNull()
      ? Form("%s/RRatio_ref_R%.1f.pdf", outputDir.Data(), kRRatioRef)
      : Form("%s/RRatio_%s_ref_R%.1f.pdf", outputDir.Data(), tag.Data(), kRRatioRef);
    pad->C->Print(fname.Data());
  }
}

// ============================================
// Event Selection (sel8 + RCT) Pass Fraction Plot
// ============================================
// Uses BC-level decomposition from GetNormDecomposition() (DATA file):
//   1. sel8 (BC-level)     = nTVXafterBC / nTVX
//   2. RCT (BC-level)      = nTVXafterBCRCT / nTVXafterBC
//   3. sel8×RCT combined   = nTVXafterBCRCT / nTVX
// All three use consistent denominators within the same BC chain.
// Per-R MC h_mccollisions pass fraction shown as cross-check markers.

void DrawRCTEfficiency() {
  std::vector<RConfig> rConfigs = GetRConfigs();
  std::vector<Color_t> rColors = GetRColorsForConfigs(rConfigs);
  TString outputDir = GetOutputDir();

  // --- Step 1: Get BC-level decomposition from DATA files ---
  // NormDecomposition is per-dataset (same for all R sharing a data file).
  // Collect unique data files and compute once per file.
  struct DatasetDecomp {
    TString dataFile;
    NormDecomposition nd;
  };
  std::vector<DatasetDecomp> datasets;
  std::set<TString> seenDataFiles;

  for (const auto& config : rConfigs) {
    TString df = config.GetDataFile();
    if (seenDataFiles.count(df)) continue;
    seenDataFiles.insert(df);

    NormDecomposition nd = GetNormDecomposition(df.Data());
    datasets.push_back({df, nd});
  }

  // Use the first valid NormDecomposition
  NormDecomposition bestND = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false};
  for (const auto& ds : datasets) {
    if (ds.nd.valid) { bestND = ds.nd; break; }
  }

  // --- Step 2: Get per-R MC-level pass fraction as cross-check ---
  struct EvtSelInfo {
    double R;
    TString label;
    double mcCombined;   // h_mccollisions bin(3.5)/bin(2.5) — MC-level cross-check
  };
  std::vector<EvtSelInfo> infos;
  std::set<TString> seenMcDirs;

  for (const auto& config : rConfigs) {
    TString key = config.GetMcFile() + ":" + config.mcDir;
    if (seenMcDirs.count(key)) continue;
    seenMcDirs.insert(key);

    double mc = GetRCTPassFraction(config.GetMcFile().Data(), config.mcDir.Data(), config.isJJ);
    infos.push_back({config.R, config.label, mc});
  }

  if (infos.empty()) {
    std::cerr << "[EvtSel] No configs to plot." << std::endl;
    return;
  }

  // --- Step 3: Print summary ---
  std::cerr << std::endl;
  std::cerr << "============================================================" << std::endl;
  std::cerr << "  Event Selection Efficiency Decomposition (BC-level)" << std::endl;
  std::cerr << "============================================================" << std::endl;
  if (bestND.valid) {
    std::cerr << Form("  eps_sel8    (BC) = %.5f  (nTVXafterBC / nTVX)", bestND.epsSel8) << std::endl;
    std::cerr << Form("  eps_RCT     (BC) = %.5f  (nTVXafterBCRCT / nTVXafterBC)", bestND.epsRCT) << std::endl;
    std::cerr << Form("  eps_sel8xRCT(BC) = %.5f  (nTVXafterBCRCT / nTVX)", bestND.epsSel8RCT) << std::endl;
    std::cerr << Form("  Runs: %d good / %d total", bestND.nRunsGood, bestND.nRunsTotal) << std::endl;
  } else {
    std::cerr << "  [NormDecomposition NOT AVAILABLE — no luminosity histograms in data]" << std::endl;
  }
  std::cerr << "------------------------------------------------------------" << std::endl;
  std::cerr << "  Per-R MC cross-check (h_mccollisions bin 3.5/2.5):" << std::endl;
  for (const auto& info : infos) {
    std::cerr << Form("    R=%-4.1f  sel8 pass = %.4f", info.R, info.mcCombined) << std::endl;
  }
  std::cerr << "============================================================" << std::endl << std::endl;

  // --- Step 4: Draw bar chart ---
  int nR = (int)infos.size();
  bool hasDecomp = bestND.valid;

  TCanvas* cRCT = new TCanvas(Form("cEvtSelEff_%d", ++GetNN()),
                               "Event Selection Efficiency Decomposition", 900, 550);
  cRCT->SetGrid(0, 1);
  cRCT->SetLeftMargin(0.10);
  cRCT->SetBottomMargin(0.15);
  cRCT->SetRightMargin(0.02);

  if (hasDecomp) {
    // --- BC-level decomposition available: show 3 bars per R ---
    // (values are same for all R sharing a dataset, but plot per R for layout)

    TH1F* hSel8RCT = new TH1F(Form("hSel8RCT_%d", GetNN()), "", nR, 0, nR);
    hSel8RCT->SetDirectory(nullptr);
    hSel8RCT->SetStats(0);
    hSel8RCT->GetYaxis()->SetTitle("BC-level efficiency");
    hSel8RCT->GetYaxis()->SetRangeUser(0.0, 1.15);
    hSel8RCT->GetYaxis()->SetTitleOffset(0.9);
    hSel8RCT->SetFillColor(kAzure - 4);
    hSel8RCT->SetLineColor(kAzure - 6);
    hSel8RCT->SetBarWidth(0.25);
    hSel8RCT->SetBarOffset(0.08);

    TH1F* hSel8 = new TH1F(Form("hSel8_%d", GetNN()), "", nR, 0, nR);
    hSel8->SetDirectory(nullptr);
    hSel8->SetFillColor(kOrange - 3);
    hSel8->SetLineColor(kOrange - 1);
    hSel8->SetBarWidth(0.25);
    hSel8->SetBarOffset(0.36);

    TH1F* hRCT = new TH1F(Form("hRCT_%d", GetNN()), "", nR, 0, nR);
    hRCT->SetDirectory(nullptr);
    hRCT->SetFillColor(kGreen - 6);
    hRCT->SetLineColor(kGreen - 3);
    hRCT->SetBarWidth(0.25);
    hRCT->SetBarOffset(0.64);

    for (int i = 0; i < nR; ++i) {
      // Find dataset for this R's data file
      NormDecomposition nd = bestND;  // fallback to best
      for (const auto& config : rConfigs) {
        if (TMath::Abs(config.R - infos[i].R) < 0.01) {
          for (const auto& ds : datasets) {
            if (ds.dataFile == config.GetDataFile() && ds.nd.valid) { nd = ds.nd; break; }
          }
          break;
        }
      }

      hSel8RCT->SetBinContent(i + 1, nd.epsSel8RCT);
      hSel8->SetBinContent(i + 1, nd.epsSel8);
      hRCT->SetBinContent(i + 1, nd.epsRCT);
      hSel8RCT->GetXaxis()->SetBinLabel(i + 1, Form("R = %.1f", infos[i].R));
    }
    hSel8RCT->GetXaxis()->SetLabelSize(0.06);
    hSel8RCT->Draw("bar");
    hSel8->Draw("bar same");
    hRCT->Draw("bar same");

    // Value labels
    TLatex latex;
    latex.SetTextAlign(21);
    latex.SetTextSize(0.028);
    for (int i = 0; i < nR; ++i) {
      double xBase = i + 0.5;
      NormDecomposition nd = bestND;
      for (const auto& config : rConfigs) {
        if (TMath::Abs(config.R - infos[i].R) < 0.01) {
          for (const auto& ds : datasets) {
            if (ds.dataFile == config.GetDataFile() && ds.nd.valid) { nd = ds.nd; break; }
          }
          break;
        }
      }
      latex.SetTextColor(kAzure - 6);
      latex.DrawLatex(xBase - 0.18, nd.epsSel8RCT + 0.02, Form("%.4f", nd.epsSel8RCT));
      latex.SetTextColor(kOrange - 1);
      latex.DrawLatex(xBase + 0.0, nd.epsSel8 + 0.02, Form("%.4f", nd.epsSel8));
      latex.SetTextColor(kGreen + 2);
      latex.DrawLatex(xBase + 0.18, nd.epsRCT + 0.02, Form("%.4f", nd.epsRCT));

      // MC cross-check marker
      latex.SetTextColor(kRed + 1);
      latex.SetTextSize(0.022);
      latex.DrawLatex(xBase, infos[i].mcCombined + 0.04, Form("MC: %.3f", infos[i].mcCombined));
      latex.SetTextSize(0.028);
    }
    latex.SetTextColor(kBlack);

    // MC cross-check markers (triangles)
    for (int i = 0; i < nR; ++i) {
      TMarker* mk = new TMarker(i + 0.22, infos[i].mcCombined, 26);  // open triangle
      mk->SetMarkerColor(kRed + 1);
      mk->SetMarkerSize(1.2);
      mk->Draw();
    }

    // Unity line
    TLine* lineUnity = new TLine(0, 1.0, nR, 1.0);
    lineUnity->SetLineColor(kBlack);
    lineUnity->SetLineStyle(2);
    lineUnity->SetLineWidth(1);
    lineUnity->Draw("lsame");

    // Legend
    TLegend* leg = new TLegend(0.12, 0.78, 0.92, 0.95);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.028);
    leg->SetNColumns(4);
    leg->AddEntry(hSel8RCT, "sel8#timesRCT (BC-level)", "f");
    leg->AddEntry(hSel8, "sel8 only (BC-level)", "f");
    leg->AddEntry(hRCT, "RCT only (BC-level)", "f");
    // Dummy marker for MC cross-check
    TH1F* dummyMC = new TH1F("dummyMCxcheck", "", 1, 0, 1);
    dummyMC->SetMarkerStyle(26);
    dummyMC->SetMarkerColor(kRed + 1);
    leg->AddEntry(dummyMC, "MC sel8 cross-check", "p");
    leg->Draw();

  } else {
    // --- No decomposition: show only MC-level combined ---
    TH1F* hCombined = new TH1F(Form("hCombined_%d", GetNN()), "", nR, 0, nR);
    hCombined->SetDirectory(nullptr);
    hCombined->SetStats(0);
    hCombined->GetYaxis()->SetTitle("Pass fraction");
    hCombined->GetYaxis()->SetRangeUser(0.0, 1.15);
    hCombined->GetYaxis()->SetTitleOffset(0.9);
    hCombined->SetFillColor(kAzure - 4);
    hCombined->SetLineColor(kAzure - 6);
    hCombined->SetBarWidth(0.5);
    hCombined->SetBarOffset(0.25);

    for (int i = 0; i < nR; ++i) {
      hCombined->SetBinContent(i + 1, infos[i].mcCombined);
      hCombined->GetXaxis()->SetBinLabel(i + 1, Form("R = %.1f", infos[i].R));
    }
    hCombined->GetXaxis()->SetLabelSize(0.06);
    hCombined->Draw("bar");

    // Value labels
    TLatex latex;
    latex.SetTextAlign(21);
    latex.SetTextSize(0.030);
    for (int i = 0; i < nR; ++i) {
      latex.SetTextColor(kAzure - 6);
      latex.DrawLatex(i + 0.5, infos[i].mcCombined + 0.02, Form("%.3f", infos[i].mcCombined));
    }
    latex.SetTextColor(kBlack);

    // Unity line
    TLine* lineUnity = new TLine(0, 1.0, nR, 1.0);
    lineUnity->SetLineColor(kBlack);
    lineUnity->SetLineStyle(2);
    lineUnity->SetLineWidth(1);
    lineUnity->Draw("lsame");

    // Legend
    TLegend* leg = new TLegend(0.12, 0.85, 0.65, 0.95);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.030);
    leg->AddEntry(hCombined, "sel8 pass fraction (MC, h_mccollisions)", "f");
    leg->Draw();

    // Warning text
    TLatex warn;
    warn.SetNDC();
    warn.SetTextFont(42);
    warn.SetTextSize(0.030);
    warn.SetTextColor(kRed + 1);
    warn.DrawLatex(0.15, 0.75, "No BC-level decomposition available");
    warn.DrawLatex(0.15, 0.70, "(missing eventselection-run3/luminosity histograms in data file)");
  }

  if (GetDRAWPLOTS()) {
    cRCT->Print(Form("%s/EvtSel_PassFraction.pdf", outputDir.Data()));
  }
}

// ============================================
// Cross-section cache: save/load helpers
// ============================================

static bool SaveXsecCache(const TString& path, int version, EXsecMode mode,
                           const std::vector<RConfig>& rConfigs,
                           const std::vector<TH1*>& vecUnfolded,
                           const std::vector<TH1*>& vecPOWHEG,
                           const std::map<std::string, std::vector<TH1*>>& vecModels)
{
  TFile* f = TFile::Open(path, "RECREATE");
  if (!f || f->IsZombie()) {
    std::cerr << "[XsecCache] ERROR: cannot create " << path << std::endl;
    delete f;
    return false;
  }

  TParameter<int> pVer("cacheVersion", version);   pVer.Write();
  TParameter<int> pMode("mode", (int)mode);         pMode.Write();
  TParameter<int> pN("nConfigs", (int)rConfigs.size()); pN.Write();

  for (size_t i = 0; i < rConfigs.size(); ++i) {
    TString dirName = Form("config_%zu", i);
    TDirectory* dir = f->mkdir(dirName);
    dir->cd();

    TParameter<double> pR("R", rConfigs[i].R);              pR.Write();
    TParameter<int> pJJ("isJJ", rConfigs[i].isJJ ? 1 : 0);  pJJ.Write();
    TNamed nMc("mcRunNumber", rConfigs[i].mcRunNumber.Data());   nMc.Write();
    TNamed nData("dataRunNumber", rConfigs[i].dataRunNumber.Data()); nData.Write();

    if (i < vecUnfolded.size() && vecUnfolded[i]) vecUnfolded[i]->Write("hUnfolded");
    if (i < vecPOWHEG.size()   && vecPOWHEG[i])   vecPOWHEG[i]->Write("hPOWHEG");
    for (const auto& mkv : vecModels) {
      if (i < mkv.second.size() && mkv.second[i]) {
        mkv.second[i]->Write(Form("hModel_%s", mkv.first.c_str()));
      }
    }
    f->cd();
  }

  f->Close();
  delete f;
  std::cerr << "[XsecCache] Saved " << rConfigs.size() << " configs to " << path << std::endl;
  return true;
}

static bool LoadXsecCache(const TString& path, int version, EXsecMode mode,
                           const std::vector<RConfig>& rConfigs,
                           std::vector<TH1*>& vecUnfolded,
                           std::vector<TH1*>& vecPOWHEG,
                           std::map<std::string, std::vector<TH1*>>& vecModels)
{
  if (gSystem->AccessPathName(path)) return false;  // file does not exist

  TFile* f = TFile::Open(path, "READ");
  if (!f || f->IsZombie()) { delete f; return false; }

  // Validate version
  TParameter<int>* pVer = (TParameter<int>*)f->Get("cacheVersion");
  if (!pVer || pVer->GetVal() != version) {
    std::cerr << "[XsecCache] Version mismatch (file=" << (pVer ? pVer->GetVal() : -1)
              << ", expected=" << version << "), recomputing." << std::endl;
    f->Close(); delete f; return false;
  }

  // Validate mode
  TParameter<int>* pMode = (TParameter<int>*)f->Get("mode");
  if (!pMode || pMode->GetVal() != (int)mode) {
    std::cerr << "[XsecCache] Mode mismatch, recomputing." << std::endl;
    f->Close(); delete f; return false;
  }

  // Validate nConfigs
  TParameter<int>* pN = (TParameter<int>*)f->Get("nConfigs");
  if (!pN || pN->GetVal() != (int)rConfigs.size()) {
    std::cerr << "[XsecCache] Config count mismatch (file=" << (pN ? pN->GetVal() : -1)
              << ", expected=" << rConfigs.size() << "), recomputing." << std::endl;
    f->Close(); delete f; return false;
  }

  vecUnfolded.clear(); vecPOWHEG.clear(); vecModels.clear();

  // Get model names for loading
  std::vector<StandaloneModel> models = GetStandaloneModels();
  for (auto& m : models) vecModels[m.name.Data()].clear();

  for (size_t i = 0; i < rConfigs.size(); ++i) {
    TString dirName = Form("config_%zu", i);
    TDirectory* dir = f->GetDirectory(dirName);
    if (!dir) {
      std::cerr << "[XsecCache] Missing directory " << dirName << ", recomputing." << std::endl;
      f->Close(); delete f; return false;
    }

    // Validate per-config metadata
    TParameter<double>* pR = (TParameter<double>*)dir->Get("R");
    TParameter<int>* pJJ = (TParameter<int>*)dir->Get("isJJ");
    TNamed* nMc = (TNamed*)dir->Get("mcRunNumber");
    TNamed* nData = (TNamed*)dir->Get("dataRunNumber");

    if (!pR || fabs(pR->GetVal() - rConfigs[i].R) > 0.001) {
      std::cerr << "[XsecCache] R mismatch at config " << i << ", recomputing." << std::endl;
      f->Close(); delete f; return false;
    }
    if (!pJJ || (pJJ->GetVal() != 0) != rConfigs[i].isJJ) {
      std::cerr << "[XsecCache] isJJ mismatch at config " << i << ", recomputing." << std::endl;
      f->Close(); delete f; return false;
    }
    if (!nMc || TString(nMc->GetTitle()) != rConfigs[i].mcRunNumber) {
      std::cerr << "[XsecCache] mcRunNumber mismatch at config " << i << ", recomputing." << std::endl;
      f->Close(); delete f; return false;
    }
    if (!nData || TString(nData->GetTitle()) != rConfigs[i].dataRunNumber) {
      std::cerr << "[XsecCache] dataRunNumber mismatch at config " << i << ", recomputing." << std::endl;
      f->Close(); delete f; return false;
    }

    // Load histograms (detach from file)
    auto loadH = [&](const char* name) -> TH1* {
      TH1* h = (TH1*)dir->Get(name);
      if (h) h->SetDirectory(0);
      return h;
    };

    vecUnfolded.push_back(loadH("hUnfolded"));
    vecPOWHEG.push_back(loadH("hPOWHEG"));
    for (auto& m : models) {
      vecModels[m.name.Data()].push_back(loadH(Form("hModel_%s", m.name.Data())));
    }
  }

  f->Close();
  delete f;
  return true;
}

void DrawCrossSection() {
  if (!GetDrawCrossSection()) return;

  // Mode is set by SetCurrentConfigSet() before calling this function
  const EXsecMode mode = GetCurrentMode();

  // Get R configurations and colors based on mode
  std::vector<RConfig> rConfigs = GetRConfigs();
  std::vector<Color_t> rColors = GetRColorsForConfigs(rConfigs);
  TString outputDir = GetOutputDir();

  bool includePOWHEG = true;             // POWHEG NLO for both modes (R=0.2–0.6)
  std::vector<StandaloneModel> models = GetStandaloneModels();

  // Load closure-derived regularization parameters (once per DrawCrossSection call)
  const char* regFile = (mode == kUE) ? kOptimalRegFileUE : kOptimalRegFileNonUE;
  auto regMap = LoadOptimalRegularization(regFile);

  if (mode == kUE) {
    ClearXsecCachesUE();
  }

  // ---- Normalization setup (per-R computation below, since data files can differ per train) ----
  int datasetYear = DetectDatasetYear(GetCurrentConfigSet().name);
  EXsecNormMode normMode = GetXsecNormMode();

  // ====================================================================
  // Phase 1: Computation (or load from cache)
  // ====================================================================
  TString cacheFile = Form("%s/XsecCache.root", outputDir.Data());
  std::vector<TH1*> vecUnfolded, vecPOWHEG;
  std::map<std::string, std::vector<TH1*>> vecModels;

  bool cacheLoaded = !kForceRecompute &&
    LoadXsecCache(cacheFile, kXsecCacheVersion, mode, rConfigs,
                  vecUnfolded, vecPOWHEG, vecModels);

  if (cacheLoaded) {
    std::cerr << "[XsecCache] Fast path: loaded " << rConfigs.size()
              << " configs from " << cacheFile << std::endl;
  } else {
    // Compute from scratch
    vecUnfolded.clear(); vecPOWHEG.clear(); vecModels.clear();
    for (auto& m : models) vecModels[m.name.Data()].clear();

    for (size_t i = 0; i < rConfigs.size(); ++i) {
      const auto& config = rConfigs[i];
      const bool isJJ = config.isJJ;
      const Int_t svdKOverride = GetSvdKForRun(config.mcRunNumber, regMap);

      Double_t nevtsData = GetNORMEVENTS() ? Nevents(config.GetDataFile().Data(), config.dataDir.Data(), GetEventObj(), 0) : 1.0;
      Double_t nevtsDataUnTrig = GetNORMEVENTS() ? Nevents(config.GetDataFile().Data(), config.dataDir.Data(), GetEventObj(), 1) : 1.0;
      Double_t nevtsMCD = GetNORMEVENTS() ? Nevents(config.GetMcFile().Data(), config.mcDir.Data(), GetEventObj(), 0, config.isJJ) : 1.0;
      Double_t nevtsMCP = GetNORMEVENTS() ? Nevents(config.GetMcFile().Data(), config.mcDir.Data(), GetEventObj(), 1, config.isJJ) : 1.0;

      TH1* hUnfolded = DrawJetMatching(config.GetMcFile().Data(), config.GetDataFile().Data(), config.label.Data(),
                                       nevtsData, nevtsDataUnTrig, nevtsMCD, nevtsMCP,
                                       nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                                       nullptr, nullptr, nullptr, rColors[i],
                                       nullptr, nullptr, nullptr, nullptr,
                                       config.dataDir.Data(), config.mcDir.Data(), "",
                                       nullptr, nullptr, nullptr,
                                       svdKOverride);

      if (!hUnfolded) {
        vecUnfolded.push_back(nullptr);
        vecPOWHEG.push_back(nullptr);
        for (auto& m : models) vecModels[m.name.Data()].push_back(nullptr);
        continue;
      }

      // ---- Compute per-R normalization from this config's data file ----
      NormDecomposition normDecomp = GetNormDecomposition(config.GetDataFile().Data(), kRCTYBin);
      TString zvtxLabel = Form("R%02d", (int)(config.R * 10 + 0.5));
      Double_t epsZvtx10 = GetZvtx10Efficiency(config.GetDataFile().Data(), config.dataDir.Data(),
                                                GetOutputDir().Data(), zvtxLabel.Data());
      Double_t epsCollCheck = GetCollToSel8RCTEff(config.GetDataFile().Data(), config.dataDir.Data());
      std::cerr << "[Normalization] R=" << config.R << " h_collisions check: eps_coll->sel8RCT = "
                << Form("%.5f", epsCollCheck) << std::endl;
      if (epsZvtx10 <= 0) {
        epsZvtx10 = 0.956;
        std::cerr << "[XsecNorm] Warning: zvtx10 fit failed for R=" << config.R
                  << ", using fallback=" << epsZvtx10 << std::endl;
      }
      XsecNormResult xsecNorm = ComputeXsecNormFactor(normDecomp, epsZvtx10, datasetYear, normMode);

      hUnfolded->SetStats(0);
      Double_t normFactor = xsecNorm.valid ? xsecNorm.normFactor : 0;
      if (!xsecNorm.valid) {
        std::cerr << "[Error] Cross-section normalization failed for R=" << config.R << std::endl;
        vecUnfolded.push_back(nullptr);
        vecPOWHEG.push_back(nullptr);
        for (auto& m : models) vecModels[m.name.Data()].push_back(nullptr);
        continue;
      }
      hUnfolded->Scale(normFactor, "width");
      Double_t deltaEta = GetDeltaEta(config.R);
      if (deltaEta > 0) {
        hUnfolded->Scale(1.0 / deltaEta);
      }
      hUnfolded->SetDirectory(0);

      // Load POWHEG NLO (both modes, R=0.2–0.6)
      TH1* hPOWHEG = nullptr;
      if (includePOWHEG) {
        hPOWHEG = LoadPowhegForR(config.R);
      }

      // Load standalone models (PYTHIA8 Monash/Rope/Shoving, Herwig7)
      for (auto& model : models) {
        TH1* hModel = LoadStandaloneModelXsec(model, config.R);
        vecModels[model.name.Data()].push_back(hModel);
      }

      vecUnfolded.push_back(hUnfolded);
      vecPOWHEG.push_back(hPOWHEG);
    }

    // Save cache
    SaveXsecCache(cacheFile, kXsecCacheVersion, mode, rConfigs,
                  vecUnfolded, vecPOWHEG, vecModels);
  }

  // ====================================================================
  // Phase 2: Drawing (always runs — applies styles, populates UE caches,
  //          draws individual + combined canvases)
  // ====================================================================

  std::vector<TH1*> allUnfoldedHists;
  std::map<std::string, std::vector<TH1*>> allModelHists;
  std::vector<TH1*> allPowhegHists;
  std::vector<TH1*> allMCTruthHists;
  std::vector<std::pair<TGraphErrors*, TGraphAsymmErrors*>> allRun2Graphs;

  for (size_t i = 0; i < rConfigs.size(); ++i) {
    const auto& config = rConfigs[i];
    const bool isJJ = config.isJJ;
    const TString ueTag = (mode == kUE) ? (isJJ ? "JJ" : "MB") : "";

    TH1* hUnfolded = (i < vecUnfolded.size()) ? vecUnfolded[i] : nullptr;
    if (!hUnfolded) continue;
    TH1* hPOWHEG  = (i < vecPOWHEG.size())  ? vecPOWHEG[i]  : nullptr;

    // Apply styles (always, so cosmetic changes work from cache)
    hUnfolded->SetStats(0);
    hUnfolded->SetLineColor(rColors[i]);
    hUnfolded->SetMarkerColor(rColors[i]);
    hUnfolded->SetMarkerStyle((mode == kUE && isJJ) ? 25 : 20);
    hUnfolded->SetLineWidth(2);

    // Style standalone models
    for (auto& model : models) {
      TH1* hM = vecModels.count(model.name.Data()) && i < vecModels[model.name.Data()].size()
                 ? vecModels[model.name.Data()][i] : nullptr;
      if (hM) {
        hM->SetLineColor(GetColorForR(config.R));
        hM->SetMarkerColor(GetColorForR(config.R));
        hM->SetLineStyle(model.lineStyle);
        hM->SetLineWidth(model.lineWidth);
        hM->SetMarkerStyle(0);
        hM->SetStats(0);
      }
    }

    if (hPOWHEG) {
      hPOWHEG->SetLineColor(GetColorForR(config.R));
      hPOWHEG->SetMarkerColor(GetColorForR(config.R));
      hPOWHEG->SetLineStyle(3);
      hPOWHEG->SetLineWidth(3);
      hPOWHEG->SetMarkerStyle(0);
      hPOWHEG->SetStats(0);
    }

    // Populate UE global caches (needed for stitching and AllR canvas)
    if (mode == kUE) {
      std::map<double, TH1*>& xsecCache = isJJ ? gUE_Xsec_JJ : gUE_Xsec_MB;
      std::map<double, TString>& cacheRun = isJJ ? gUE_McRun_JJ : gUE_McRun_MB;
      if (xsecCache.find(config.R) != xsecCache.end()) {
        std::cerr << "[Error] Duplicate UE " << (isJJ ? "JJ" : "MB") << " entry for R=" << config.R
                  << ". Stitching would be ambiguous. Not caching this one." << std::endl;
      } else {
        TH1* hStored = (TH1*)hUnfolded->Clone(Form("xSection_UE_%s_R%.1f", ueTag.Data(), config.R));
        hStored->SetDirectory(0);
        hStored->SetStats(0);
        xsecCache[config.R] = hStored;
        cacheRun[config.R] = config.mcRunNumber;
      }

      // Cache models for UE AllR canvas
      for (auto& model : models) {
        TH1* hM = vecModels.count(model.name.Data()) && i < vecModels[model.name.Data()].size()
                   ? vecModels[model.name.Data()][i] : nullptr;
        if (hM) {
          if (gUE_Models[model.name.Data()].count(config.R))
            delete gUE_Models[model.name.Data()][config.R];
          gUE_Models[model.name.Data()][config.R] = (TH1*)hM->Clone(
            Form("Model_%s_stored_R%.1f", model.name.Data(), config.R));
          gUE_Models[model.name.Data()][config.R]->SetDirectory(0);
        }
      }
    }

    if (mode == kUE && hPOWHEG) {
      if (gUE_POWHEG.count(config.R)) delete gUE_POWHEG[config.R];
      gUE_POWHEG[config.R] = (TH1*)hPOWHEG->Clone(Form("POWHEG_stored_R%.1f", config.R));
      gUE_POWHEG[config.R]->SetDirectory(0);
      std::cerr << "[POWHEG-UE] Stored POWHEG for R=" << config.R
                << ", integral=" << Form("%.4e", gUE_POWHEG[config.R]->Integral())
                << ", entries=" << gUE_POWHEG[config.R]->GetEntries() << std::endl;
    } else if (mode == kUE && !hPOWHEG) {
      std::cerr << "[POWHEG-UE] WARNING: hPOWHEG is null for R=" << config.R << std::endl;
    }

    // Cache MC truth for UE AllR canvas
    if (mode == kUE) {
      TH1* hMCTruthUE = LoadMCTruthXsec(config);
      if (hMCTruthUE) {
        if (gUE_MCTruth.count(config.R)) delete gUE_MCTruth[config.R];
        gUE_MCTruth[config.R] = hMCTruthUE;
        hMCTruthUE->SetDirectory(0);
      }
    }

    // Load Run2 data (cheap, always recomputed from hardcoded arrays)
    TGraphErrors* grRun2Stat = nullptr;
    TGraphAsymmErrors* grRun2Sys = nullptr;
    if (config.R >= 0.2) {
      if (mode == kUE) {
        auto run2Graphs = CreateRun2UEGraphs(config.R);
        grRun2Stat = run2Graphs.first;
        grRun2Sys = run2Graphs.second;
      } else {
        auto run2Graphs = CreateRun2Graphs(config.R);
        grRun2Stat = run2Graphs.first;
        grRun2Sys = run2Graphs.second;
      }
    }

    // ========== Create Filipad2 canvas for this R ==========
    int nRatioPads = 2;
    float ratio2Val = 0.25;

    TString padName;
    if (mode == kUE) {
      padName = Form("UE_%s_CrossSection_R%.1f_MC%s", ueTag.Data(), config.R, config.mcRunNumber.Data());
    } else {
      padName = Form("CrossSection_R%.1f", config.R);
    }

    Filipad2* padR = new Filipad2(padName.Data(), ++GetNN(), 2.2, 0.5, ratio2Val, 100, 50, 0.7, 5, 1);
    padR->Draw();
    TPad* mainPad = padR->GetPad(1);
    TPad* ratioPad1 = padR->GetPad(2);
    TPad* ratioPad2 = padR->GetPad(3);

    optFili(*mainPad, 0, 1, 0, 1);
    optFili(*ratioPad1, 0, 0, 0, 0);
    if (ratioPad2) {
      optFili(*ratioPad2, 0, 0, 0, 0);
    }

    TLegend* legR = new TLegend(0.5, 0.6, 0.75, 0.95, Form("#it{R} = %.1f", config.R), "brNDC");
    legR->SetTextSize(0.04);
    legR->SetBorderSize(0);

    // ========== Main pad: yield plot ==========
    mainPad->cd();
    gPad->SetGrid(0);

    bool firstDraw = true;
    TH1* hFirst = hUnfolded;
    if (!hFirst) {
      for (auto& model : models) {
        TH1* hM = vecModels.count(model.name.Data()) && i < vecModels[model.name.Data()].size()
                   ? vecModels[model.name.Data()][i] : nullptr;
        if (hM) { hFirst = hM; break; }
      }
    }

    // Draw Run2 systematic box first
    if (grRun2Sys) {
      grRun2Sys->Draw("E2");
      if (firstDraw && hFirst) {
        hset(*hFirst, GetJetPtGenTitleX(), "d^{2}#sigma/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
        hFirst->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
        hFirst->GetYaxis()->SetRangeUser(1e-6, 1e2);
        hFirst->GetXaxis()->SetLabelSize(0);
        hFirst->GetXaxis()->SetTitleSize(0);
        hFirst->Draw("AXIS");
        firstDraw = false;
      }
    }

    // Draw Run2 statistical errors
    if (grRun2Stat) {
      if (firstDraw) {
        double xMin = TMath::MinElement(grRun2Stat->GetN(), grRun2Stat->GetX());
        double xMax = TMath::MaxElement(grRun2Stat->GetN(), grRun2Stat->GetX());
        grRun2Stat->GetXaxis()->SetLimits(xMin, xMax);
        grRun2Stat->GetYaxis()->SetRangeUser(1e-6, 1e2);
        grRun2Stat->Draw("APZ");
        firstDraw = false;
      } else {
        grRun2Stat->Draw("PZ same");
      }
      legR->AddEntry(grRun2Stat, "Run 2 data", "pe");
    }

    // Standalone models and POWHEG only shown in ratio pad (pad 2), not in main yield pad
    // // Draw standalone models (PYTHIA8 Monash/Rope/Shoving, Herwig7)
    // for (size_t im = 0; im < models.size(); ++im) {
    //   TH1* hM = vecModels.count(models[im].name.Data()) && i < vecModels[models[im].name.Data()].size()
    //              ? vecModels[models[im].name.Data()][i] : nullptr;
    //   if (!hM) continue;
    //   if (firstDraw) {
    //     hset(*hM, GetJetPtGenTitleX(), "d^{2}#sigma/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    //     hM->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
    //     hM->GetXaxis()->SetLabelSize(0);
    //     hM->GetXaxis()->SetTitleSize(0);
    //     hM->Draw("l");
    //     firstDraw = false;
    //   } else {
    //     hM->Draw("l same");
    //   }
    //   legR->AddEntry(hM, models[im].legLabel.Data(), "l");
    // }

    // // Draw POWHEG NLO (both modes)
    // if (hPOWHEG) {
    //   if (firstDraw) {
    //     hset(*hPOWHEG, GetJetPtGenTitleX(), "d^{2}#sigma/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    //     hPOWHEG->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
    //     hPOWHEG->GetXaxis()->SetLabelSize(0);
    //     hPOWHEG->GetXaxis()->SetTitleSize(0);
    //     hPOWHEG->Draw("l");
    //     firstDraw = false;
    //   } else {
    //     hPOWHEG->Draw("l same");
    //   }
    //   legR->AddEntry(hPOWHEG, "POWHEG NLO", "l");
    // }

    // Draw Run3 unfolded data
    if (firstDraw) {
      hset(*hUnfolded, GetJetPtGenTitleX(), "d^{2}#sigma/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
      hUnfolded->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
      hUnfolded->GetXaxis()->SetLabelSize(0);
      hUnfolded->GetXaxis()->SetTitleSize(0);
      hUnfolded->Draw("pe");
    } else {
      hUnfolded->Draw("pe same");
    }

    if (mode == kUE) {
      legR->AddEntry(hUnfolded, Form("Run 3 data (%s response)", ueTag.Data()), "pe");
    } else {
      legR->AddEntry(hUnfolded, "Run 3 data", "pe");
    }

    legR->Draw();

    // ALICE figure legend
    ALICEfigureLegend("ALICE WIP", 0.2, 0.70, 0.4, 0.95, 0.2, 0.05, 0.4, 0.15, 0.04, config.R);

    // ========== Ratio pad 1: Run3 / Run2 ==========
    if (grRun2Stat && config.R >= 0.2) {
      ratioPad1->cd();
      gPad->SetGrid(0);

      TH1* hRatio1 = nullptr;
      if (mode == kUE) {
        hRatio1 = CreateRatioTH1vsTGraphUE(hUnfolded, grRun2Stat, Form("ratio_Run3_Run2_%s_R%.1f", ueTag.Data(), config.R));
      } else {
        hRatio1 = CreateRatioTH1vsTGraph(hUnfolded, grRun2Stat, Form("ratio_Run3_Run2_R%.1f", config.R));
      }

      if (hRatio1) {
        hRatio1->SetStats(0);
        hset(*hRatio1, GetJetPtGenTitleX(), "Run 3 / Run 2", 1.2, 0.45, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 503);
        hRatio1->SetLineColor(rColors[i]);
        hRatio1->SetMarkerColor(rColors[i]);
        hRatio1->SetMarkerStyle((mode == kUE && isJJ) ? 25 : 20);
        hRatio1->SetLineWidth(2);
        hRatio1->GetYaxis()->SetRangeUser(0.8, 1.5);
        hRatio1->Draw("pe");

        TLine* line1 = new TLine(GetPlotPtMin(), 1.0, GetPlotPtMax(), 1.0);
        line1->SetLineColor(kBlack);
        line1->SetLineStyle(2);
        line1->SetLineWidth(1);
        line1->Draw("lsame");
      }
    }

    // ========== Ratio pad 2: Model/Run3 ==========
    {
      bool hasAnyModel = (hPOWHEG != nullptr);
      for (auto& model : models) {
        if (vecModels.count(model.name.Data()) && i < vecModels[model.name.Data()].size()
            && vecModels[model.name.Data()][i]) { hasAnyModel = true; break; }
      }
      if (ratioPad2 && hasAnyModel) {
        ratioPad2->cd();
        gPad->SetGrid(0);

        bool firstRat2 = true;

        // Draw standalone model / Run3 ratios (use model-specific colors for per-R)
        for (size_t im = 0; im < models.size(); ++im) {
          TH1* hM = vecModels.count(models[im].name.Data()) && i < vecModels[models[im].name.Data()].size()
                     ? vecModels[models[im].name.Data()][i] : nullptr;
          if (!hM || !hUnfolded) continue;
          TH1* hRat = (TH1*)hM->Clone(Form("ratio_%s_Run3_R%.1f", models[im].name.Data(), config.R));
          hRat->Divide(hM, hUnfolded, 1., 1., "");
          hRat->SetStats(0);
          hRat->SetLineColor(models[im].lineColor);
          hRat->SetMarkerColor(models[im].lineColor);
          hRat->SetLineStyle(1);
          hRat->SetLineWidth(2);
          hRat->SetMarkerStyle(0);
          if (firstRat2) {
            hset(*hRat, GetJetPtGenTitleX(), "Model / Run 3", 1.2, 0.45, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 503);
            hRat->GetXaxis()->SetRangeUser(5, 140);
            hRat->GetYaxis()->SetRangeUser(0.5, 2.0);
            hRat->Draw("l");
            firstRat2 = false;
          } else {
            hRat->Draw("l same");
          }
        }

        // Draw POWHEG/Run3
        if (hPOWHEG) {
          if (hPOWHEG->GetNbinsX() != hUnfolded->GetNbinsX()) {
            std::cerr << "[POWHEG-Ratio] WARNING: bin mismatch R=" << config.R
                      << " POWHEG=" << hPOWHEG->GetNbinsX() << " vs hUnfolded=" << hUnfolded->GetNbinsX()
                      << " — Divide will fail silently!" << std::endl;
          }
          TH1* hRatioPow = (TH1*)hPOWHEG->Clone(Form("ratio_POWHEG_Run3_R%.1f", config.R));
          hRatioPow->Divide(hPOWHEG, hUnfolded, 1., 1., "");
          hRatioPow->SetStats(0);
          hRatioPow->SetLineColor(kPowhegColor);
          hRatioPow->SetMarkerColor(kPowhegColor);
          hRatioPow->SetLineStyle(1);
          hRatioPow->SetLineWidth(2);
          hRatioPow->SetMarkerStyle(0);
          if (firstRat2) {
            hset(*hRatioPow, GetJetPtGenTitleX(), "Model / Run 3", 1.2, 0.45, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 503);
            hRatioPow->GetXaxis()->SetRangeUser(5, 140);
            hRatioPow->GetYaxis()->SetRangeUser(0.5, 2.0);
            hRatioPow->Draw("l");
            firstRat2 = false;
          } else {
            hRatioPow->Draw("l same");
          }
        }

        TLine* line2 = new TLine(5, 1.0, 140, 1.0);
        line2->SetLineColor(kBlack);
        line2->SetLineStyle(2);
        line2->SetLineWidth(1);
        line2->Draw("lsame");

        // Small legend with model-specific colors in ratio pad
        TLegend* legRat2 = new TLegend(0.55, 0.55, 0.93, 0.95, NULL, "brNDC");
        legRat2->SetTextSize(0.06);
        legRat2->SetBorderSize(0);
        legRat2->SetFillStyle(0);
        for (auto& model : models) {
          TH1F* dumR = new TH1F(Form("dumRat_%s_R%.1f", model.name.Data(), config.R), "", 1, 0, 1);
          dumR->SetLineColor(model.lineColor);
          dumR->SetLineStyle(1);
          dumR->SetLineWidth(2);
          legRat2->AddEntry(dumR, model.legLabel.Data(), "l");
        }
        TH1F* dumPow = new TH1F(Form("dumRat_POWHEG_R%.1f", config.R), "", 1, 0, 1);
        dumPow->SetLineColor(kPowhegColor);
        dumPow->SetLineStyle(1);
        dumPow->SetLineWidth(2);
        legRat2->AddEntry(dumPow, "POWHEG NLO", "l");
        legRat2->Draw();
      }
    }

    if (mode == kUE) {
      mainPad->cd();
      AddUESubtractedLabel(mainPad, 0.70, 0.30, 0.045);
    }

    if (GetDRAWPLOTS()) {
      if (mode == kUE) {
        padR->C->Print(Form("%s/UE_CrossSection_%s_MC%s_R%.1f.pdf", outputDir.Data(), ueTag.Data(), config.mcRunNumber.Data(), config.R));
      } else {
        padR->C->Print(Form("%s/CrossSection_R%.1f.pdf", outputDir.Data(), config.R));
      }
    }

    // Store for combined canvas (NonUE only)
    if (mode == kNonUE) {
      allUnfoldedHists.push_back((TH1*)hUnfolded->Clone(Form("xSection_R%.1f_stored", config.R)));
      for (auto& model : models) {
        TH1* hM = vecModels.count(model.name.Data()) && i < vecModels[model.name.Data()].size()
                   ? vecModels[model.name.Data()][i] : nullptr;
        if (hM) {
          allModelHists[model.name.Data()].push_back((TH1*)hM->Clone(Form("model_%s_R%.1f_stored", model.name.Data(), config.R)));
        } else {
          allModelHists[model.name.Data()].push_back(nullptr);
        }
      }
      if (hPOWHEG) {
        allPowhegHists.push_back((TH1*)hPOWHEG->Clone(Form("powheg_R%.1f_stored", config.R)));
      } else {
        allPowhegHists.push_back(nullptr);
      }
      allRun2Graphs.push_back({grRun2Stat, grRun2Sys});

      // MC particle-level truth as cross-section
      TH1* hMCTruth = LoadMCTruthXsec(config);
      allMCTruthHists.push_back(hMCTruth);
    }
  }

  // UE mode: build stitched spectra (MB low pT + JJ high pT) and draw them.
  // Then return (no combined all-R canvas for UE at the moment).
  if (mode == kUE) {
    if (kEnableMBJJStitching) {
      // Build stitched spectra per R (requires both MB and JJ for that R)
      for (const auto& kv : gUE_Xsec_MB) {
        const double R = kv.first;
        const TH1* hMB = kv.second;
        auto itJJ = gUE_Xsec_JJ.find(R);
        if (itJJ == gUE_Xsec_JJ.end()) {
          std::cerr << "[Warning] Stitching skipped for R=" << R << ": missing JJ unfolded spectrum." << std::endl;
          continue;
        }
        const TH1* hJJ = itJJ->second;

        TH1* hStitched = StitchMBJJByPtSwitch(hMB, hJJ, kPtSwitchStitch, Form("xSection_UE_Stitched_R%.1f", R));
        if (!hStitched) {
          std::cerr << "[Error] Failed to stitch spectra for R=" << R << std::endl;
          continue;
        }
        gUE_Xsec_Stitched[R] = hStitched;

        if (!kDrawStitchedCrossSection) {
          continue;
        }

        // --- Plot stitched cross section for this R ---
        TGraphErrors* grRun2Stat = nullptr;
        TGraphAsymmErrors* grRun2Sys = nullptr;
        if (R >= 0.2) {
          auto run2Graphs = CreateRun2UEGraphs(R);
          grRun2Stat = run2Graphs.first;
          grRun2Sys = run2Graphs.second;
        }

        Filipad2* padR = new Filipad2(Form("UE_Stitched_CrossSection_R%.1f", R), ++GetNN(), 2.2, 0.3, 0.0, 100, 50, 0.7, 5, 1);
        padR->Draw();
        TPad* mainPad = padR->GetPad(1);
        TPad* ratioPad1 = padR->GetPad(2);
        optFili(*mainPad, 0, 1, 0, 1);
        optFili(*ratioPad1, 0, 0, 0, 0);

        TLegend* legR = new TLegend(0.5, 0.6, 0.9, 0.95, Form("#it{R} = %.1f", R), "brNDC");
        legR->SetTextSize(0.04);
        legR->SetBorderSize(0);

        mainPad->cd();
        gPad->SetGrid(0);

        bool firstDraw = true;
        TH1* hFirst = hStitched;

        if (grRun2Sys) {
          grRun2Sys->Draw("E2");
          if (firstDraw && hFirst) {
            hset(*hFirst, GetJetPtGenTitleX(), "d^{2}#sigma/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
            hFirst->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
            hFirst->GetYaxis()->SetRangeUser(1e-6, 1e2);
            hFirst->GetXaxis()->SetLabelSize(0);
            hFirst->GetXaxis()->SetTitleSize(0);
            hFirst->Draw("AXIS");
            firstDraw = false;
          }
        }

        if (grRun2Stat) {
          if (firstDraw) {
            double xMin = TMath::MinElement(grRun2Stat->GetN(), grRun2Stat->GetX());
            double xMax = TMath::MaxElement(grRun2Stat->GetN(), grRun2Stat->GetX());
            grRun2Stat->GetXaxis()->SetLimits(xMin, xMax);
            grRun2Stat->GetYaxis()->SetRangeUser(1e-6, 1e2);
            grRun2Stat->Draw("APZ");
            firstDraw = false;
          } else {
            grRun2Stat->Draw("PZ same");
          }
          legR->AddEntry(grRun2Stat, "Run 2 data", "pe");
        }

        if (firstDraw) {
          hset(*hStitched, GetJetPtGenTitleX(), "d^{2}#sigma/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
          hStitched->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
          hStitched->GetXaxis()->SetLabelSize(0);
          hStitched->GetXaxis()->SetTitleSize(0);
          hStitched->SetMarkerStyle(21);
          hStitched->Draw("pe");
        } else {
          hStitched->SetMarkerStyle(21);
          hStitched->Draw("pe same");
        }

        TString mbRun = gUE_McRun_MB.count(R) ? gUE_McRun_MB[R] : "";
        TString jjRun = gUE_McRun_JJ.count(R) ? gUE_McRun_JJ[R] : "";
        if (!mbRun.IsNull() && !jjRun.IsNull()) {
          legR->AddEntry(hStitched, Form("Run 3 data (stitched: MB MC%s < %.0f, JJ MC%s #geq %.0f)", mbRun.Data(), kPtSwitchStitch, jjRun.Data(), kPtSwitchStitch), "pe");
        } else {
          legR->AddEntry(hStitched, Form("Run 3 data (stitched at %.0f GeV)", kPtSwitchStitch), "pe");
        }

        legR->Draw();

        // Ratio pad: Run3 / Run2
        if (grRun2Stat && R >= 0.2) {
          ratioPad1->cd();
          gPad->SetGrid(0);
          TH1* hRatio1 = CreateRatioTH1vsTGraphUE(hStitched, grRun2Stat, Form("ratio_Run3_Run2_Stitched_R%.1f", R));
          if (hRatio1) {
            hRatio1->SetStats(0);
            hset(*hRatio1, GetJetPtGenTitleX(), "Run 3 / Run 2", 1.2, 0.45, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 503);
            hRatio1->SetLineColor(kBlack);
            hRatio1->SetMarkerColor(kBlack);
            hRatio1->SetMarkerStyle(21);
            hRatio1->SetLineWidth(2);
            hRatio1->GetYaxis()->SetRangeUser(0.8, 1.5);
            hRatio1->Draw("pe");
            TLine* line1 = new TLine(GetPlotPtMin(), 1.0, GetPlotPtMax(), 1.0);
            line1->SetLineColor(kBlack);
            line1->SetLineStyle(2);
            line1->SetLineWidth(1);
            line1->Draw("lsame");
          }
        }

        // UE label
        mainPad->cd();
        AddUESubtractedLabel(mainPad, 0.50, 0.90, 0.045);

        if (GetDRAWPLOTS()) {
          padR->C->Print(Form("%s/UE_CrossSection_Stitched_switch%.0f_R%.1f.pdf", outputDir.Data(), kPtSwitchStitch, R));
        }
      }

      // Warn about JJ-only R values
      for (const auto& kv : gUE_Xsec_JJ) {
        const double R = kv.first;
        if (gUE_Xsec_MB.find(R) == gUE_Xsec_MB.end()) {
          std::cerr << "[Warning] Stitching skipped for R=" << R << ": missing MB unfolded spectrum." << std::endl;
        }
      }
    }

    // R-ratio from stitched (or MB fallback) UE spectra
    if (!gUE_Xsec_Stitched.empty()) {
      DrawRRatioPlots(gUE_Xsec_Stitched, rColors, outputDir, "UE_Stitched");
    } else if (!gUE_Xsec_MB.empty()) {
      std::cerr << "[Info] Using MB spectra for R-ratio (no stitched available)" << std::endl;
      DrawRRatioPlots(gUE_Xsec_MB, rColors, outputDir, "UE_MB");
    }

    // ========== UE Combined canvas with all R values (Filipad2 with 3 pads) ==========
    {
      // Choose best available spectra: stitched > MB-only
      const std::map<double, TH1*>* useMap = nullptr;
      TString ueTag;
      if (!gUE_Xsec_Stitched.empty()) {
        useMap = &gUE_Xsec_Stitched;
        ueTag = "Stitched";
      } else if (!gUE_Xsec_MB.empty()) {
        useMap = &gUE_Xsec_MB;
        ueTag = "MB";
      }

      if (useMap && !useMap->empty()) {
        // Build ordered vectors from map (std::map is sorted by R)
        std::vector<double> ueRvals;
        std::vector<TH1*> ueHists;
        std::vector<TH1*> uePowhegHists;
        std::vector<TH1*> ueMCTruthHists;
        std::map<std::string, std::vector<TH1*>> ueModelHists;  // per model name
        std::vector<std::pair<TGraphErrors*, TGraphAsymmErrors*>> ueRun2Graphs;
        auto ueModels = GetStandaloneModels();
        std::cerr << "[AllR-UE] gUE_POWHEG has " << gUE_POWHEG.size() << " entries:";
        for (const auto& kv : gUE_POWHEG) std::cerr << " R=" << kv.first;
        std::cerr << std::endl;
        for (const auto& kv : *useMap) {
          double R = kv.first;
          ueRvals.push_back(R);
          ueHists.push_back(kv.second);
          bool hasPow = gUE_POWHEG.count(R) > 0;
          uePowhegHists.push_back(hasPow ? gUE_POWHEG[R] : nullptr);
          ueMCTruthHists.push_back(gUE_MCTruth.count(R) ? gUE_MCTruth[R] : nullptr);
          for (auto& m : ueModels) {
            bool hasModel = gUE_Models.count(m.name.Data()) && gUE_Models[m.name.Data()].count(R);
            ueModelHists[m.name.Data()].push_back(hasModel ? gUE_Models[m.name.Data()][R] : nullptr);
          }
          std::cerr << "[AllR-UE] R=" << R << " POWHEG=" << (hasPow ? "YES" : "NO") << std::endl;
          if (R >= 0.2) {
            auto run2G = CreateRun2UEGraphs(R);
            ueRun2Graphs.push_back(run2G);
          } else {
            ueRun2Graphs.push_back({nullptr, nullptr});
          }
        }

        bool hasTheoryOverlay = false;
        for (size_t i = 0; i < uePowhegHists.size(); ++i) {
          if (uePowhegHists[i]) { hasTheoryOverlay = true; break; }
        }
        if (!hasTheoryOverlay) {
          for (auto& m : ueModels) {
            for (auto& h : ueModelHists[m.name.Data()]) {
              if (h) { hasTheoryOverlay = true; break; }
            }
            if (hasTheoryOverlay) break;
          }
        }

        // ========== UE Canvas A: Run3 vs Run2 only ==========
        {
          const int uTitSz = 18, uLabSz = 16, uLegSz = 14;
          double umL = 0.15, umR = 0.03;
          double mainFracUA = 0.65, ratioFracUA = 0.35;

          TCanvas* canUA = new TCanvas(Form("CrossSection_AllR_UE_Run2_%s_%d", ueTag.Data(), ++GetNN()),
                                       "Cross Section AllR UE: Run3 vs Run2", 800, 700);
          canUA->SetFillStyle(4000); canUA->SetFillColor(10);
          gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
          canUA->SetTopMargin(0.); canUA->SetBottomMargin(0.);
          canUA->SetLeftMargin(0.); canUA->SetRightMargin(0.);
          canUA->Draw();

          TPad* mainPadUA = new TPad("upadMainUA", "Main", 0, ratioFracUA, 1, 1.0, 0);
          mainPadUA->SetTopMargin(0.04 / mainFracUA);
          mainPadUA->SetBottomMargin(0.001);
          mainPadUA->SetLeftMargin(umL); mainPadUA->SetRightMargin(umR);
          mainPadUA->Draw();

          TPad* ratioPadUA = new TPad("upadRatioUA", "Run3/Run2", 0, 0, 1, ratioFracUA, 0);
          ratioPadUA->SetTopMargin(0.001);
          ratioPadUA->SetBottomMargin(0.12 / ratioFracUA);
          ratioPadUA->SetLeftMargin(umL); ratioPadUA->SetRightMargin(umR);
          ratioPadUA->Draw();

          optFili(*mainPadUA, 0, 1, 0, 1);
          optFili(*ratioPadUA, 0, 0, 0, 0);

          // Main pad: Run3 + Run2 spectra
          mainPadUA->cd();
          gPad->SetGrid(0); gPad->SetLogy();

          TH1F* uframeA = new TH1F(Form("uframeA_%s", ueTag.Data()), "", 100, GetPlotPtMin(), GetPlotPtMax());
          uframeA->SetDirectory(nullptr); uframeA->SetStats(0);
          uframeA->SetMinimum(1e-6); uframeA->SetMaximum(1e2);
          uframeA->GetXaxis()->SetTitle(""); uframeA->GetXaxis()->SetLabelSize(0);
          uframeA->GetYaxis()->SetTitle("d^{2}#sigma/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]");
          uframeA->GetXaxis()->SetTitleFont(43); uframeA->GetYaxis()->SetTitleFont(43);
          uframeA->GetXaxis()->SetLabelFont(43); uframeA->GetYaxis()->SetLabelFont(43);
          uframeA->GetXaxis()->SetTitleSize(uTitSz); uframeA->GetYaxis()->SetTitleSize(uTitSz);
          uframeA->GetYaxis()->SetLabelSize(uLabSz);
          uframeA->GetYaxis()->SetTitleOffset(2.0);
          uframeA->GetXaxis()->SetNdivisions(510); uframeA->GetYaxis()->SetNdivisions(505);
          uframeA->Draw("AXIS");

          for (size_t i = 0; i < ueRun2Graphs.size(); ++i) {
            if (ueRun2Graphs[i].second) ueRun2Graphs[i].second->Draw("E2 same");
          }
          for (size_t i = 0; i < ueRun2Graphs.size(); ++i) {
            if (ueRun2Graphs[i].first) ueRun2Graphs[i].first->Draw("PZ same");
          }
          for (size_t i = 0; i < ueHists.size(); ++i) {
            if (ueHists[i]) {
              Color_t col = GetColorForR(ueRvals[i]);
              ueHists[i]->SetMarkerColor(col); ueHists[i]->SetLineColor(col);
              ueHists[i]->SetMarkerStyle(21);
              ueHists[i]->Draw("pe same");
            }
          }

          // Legends
          TLegend* ulegA = new TLegend(0.50, 0.60, 0.92, 0.90, NULL, "brNDC");
          ulegA->SetTextFont(43); ulegA->SetTextSize(uLegSz);
          ulegA->SetBorderSize(0); ulegA->SetFillStyle(0);
          TH1F* udA_r3 = new TH1F(Form("udA_r3_%s", ueTag.Data()), "", 1, 0, 1);
          udA_r3->SetLineColor(kBlack); udA_r3->SetMarkerStyle(21);
          ulegA->AddEntry(udA_r3, "Run 3", "lpe");
          TH1F* udA_r2 = new TH1F(Form("udA_r2_%s", ueTag.Data()), "", 1, 0, 1);
          udA_r2->SetFillColorAlpha(kBlack, 0.35); udA_r2->SetMarkerStyle(29);
          ulegA->AddEntry(udA_r2, "Run 2", "lpef");
          ulegA->Draw();

          TLegend* uRlegA = new TLegend(0.70, 0.35, 0.92, 0.58, NULL, "brNDC");
          uRlegA->SetTextFont(43); uRlegA->SetTextSize(uLegSz);
          uRlegA->SetBorderSize(0); uRlegA->SetFillStyle(0);
          for (size_t i = 0; i < ueRvals.size(); ++i) {
            TH1F* udR = new TH1F(Form("udRA_%s_%d", ueTag.Data(), (int)i), "", 1, 0, 1);
            udR->SetLineColor(GetColorForR(ueRvals[i])); udR->SetMarkerColor(GetColorForR(ueRvals[i]));
            udR->SetMarkerStyle(21); udR->SetLineWidth(2);
            uRlegA->AddEntry(udR, Form("#it{R} = %.1f", ueRvals[i]), "lpe");
          }
          uRlegA->Draw();
          ALICEfigureLegend("ALICE WIP", 0.2, 0.70, 0.4, 0.95, 0.2, 0.05, 0.4, 0.15, 0.04);
          AddUESubtractedLabel(mainPadUA, 0.70, 0.30, 0.045);

          // Ratio pad: Run3/Run2
          ratioPadUA->cd();
          gPad->SetGrid(0);
          TH1F* uframeRA = new TH1F(Form("uframeRA_%s", ueTag.Data()), "", 100, GetPlotPtMin(), GetPlotPtMax());
          uframeRA->SetDirectory(nullptr); uframeRA->SetStats(0);
          uframeRA->SetMinimum(0.5); uframeRA->SetMaximum(2.0);
          uframeRA->GetXaxis()->SetTitle(GetJetPtGenTitleX());
          uframeRA->GetYaxis()->SetTitle("Run 3 / Run 2");
          uframeRA->GetXaxis()->SetTitleFont(43); uframeRA->GetYaxis()->SetTitleFont(43);
          uframeRA->GetXaxis()->SetLabelFont(43); uframeRA->GetYaxis()->SetLabelFont(43);
          uframeRA->GetXaxis()->SetTitleSize(uTitSz); uframeRA->GetYaxis()->SetTitleSize(uTitSz);
          uframeRA->GetXaxis()->SetLabelSize(uLabSz); uframeRA->GetYaxis()->SetLabelSize(uLabSz);
          uframeRA->GetYaxis()->SetTitleOffset(2.0);
          uframeRA->GetXaxis()->SetTitleOffset(3.0);
          uframeRA->GetXaxis()->SetNdivisions(510); uframeRA->GetYaxis()->SetNdivisions(505);
          uframeRA->Draw("AXIS");

          for (size_t i = 0; i < ueHists.size(); ++i) {
            if (ueHists[i] && i < ueRun2Graphs.size() && ueRun2Graphs[i].first && ueRvals[i] >= 0.2) {
              TH1* hRatUE = CreateRatioTH1vsTGraphUE(ueHists[i], ueRun2Graphs[i].first, Form("ratio_Run3_Run2_AllR_UE_A_%zu", i));
              if (hRatUE) {
                Color_t col = GetColorForR(ueRvals[i]);
                hRatUE->SetStats(0); hRatUE->SetLineColor(col); hRatUE->SetMarkerColor(col);
                hRatUE->SetMarkerStyle(21); hRatUE->SetLineWidth(2);
                hRatUE->Draw("pe same");
              }
            }
          }
          TLine* lineUA = new TLine(GetPlotPtMin(), 1.0, GetPlotPtMax(), 1.0);
          lineUA->SetLineColor(kBlack); lineUA->SetLineStyle(2); lineUA->SetLineWidth(1);
          lineUA->Draw("lsame");

          if (GetDRAWPLOTS()) {
            canUA->Update(); canUA->Modified();
            canUA->Print(Form("%s/CrossSection_AllR_UE_Run2_%s.pdf", outputDir.Data(), ueTag.Data()));
          }
        }

        // ========== UE Canvas B: Run3 vs MC-truth vs Models ==========
        {
          const int uTitSz = 18, uLabSz = 16, uLegSz = 14;
          double umL = 0.15, umR = 0.03;
          // Ratio pads: MC-truth + standalone models + POWHEG
          int ueNRatioPads = (int)ueModels.size() + 2;  // +1 MC-truth, +1 POWHEG
          double mainFracUB = 0.35;
          double ratFracUB = (1.0 - mainFracUB) / ueNRatioPads;

          TCanvas* canUB = new TCanvas(Form("CrossSection_AllR_UE_Models_%s_%d", ueTag.Data(), ++GetNN()),
                                       "Cross Section AllR UE: Models", 800, 1600);
          canUB->SetFillStyle(4000); canUB->SetFillColor(10);
          gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
          canUB->SetTopMargin(0.); canUB->SetBottomMargin(0.);
          canUB->SetLeftMargin(0.); canUB->SetRightMargin(0.);
          canUB->Draw();

          TPad* mainPadUB = new TPad("upadMainUB", "Main", 0, 1.0 - mainFracUB, 1, 1.0, 0);
          mainPadUB->SetTopMargin(0.04 / mainFracUB);
          mainPadUB->SetBottomMargin(0.001);
          mainPadUB->SetLeftMargin(umL); mainPadUB->SetRightMargin(umR);
          mainPadUB->Draw();

          std::vector<TPad*> uePadsB;
          for (int mp = 0; mp < ueNRatioPads; ++mp) {
            double top = (1.0 - mainFracUB) - mp * ratFracUB;
            double bot = top - ratFracUB;
            bool isLast = (mp == ueNRatioPads - 1);
            TString mName;
            if (mp == 0) mName = "MCTruth";
            else if (mp <= (int)ueModels.size()) mName = ueModels[mp - 1].name;
            else mName = "POWHEG";
            TPad* pad = new TPad(Form("upadB_%s_%s", mName.Data(), ueTag.Data()), mName.Data(), 0, bot, 1, top, 0);
            pad->SetTopMargin(0.001);
            pad->SetBottomMargin(isLast ? 0.12 / ratFracUB : 0.001);
            pad->SetLeftMargin(umL); pad->SetRightMargin(umR);
            pad->Draw();
            uePadsB.push_back(pad);
          }

          optFili(*mainPadUB, 0, 1, 0, 1);
          for (auto* p : uePadsB) optFili(*p, 0, 0, 0, 0);

          // Main pad: yield with models and MC truth
          mainPadUB->cd();
          gPad->SetGrid(0); gPad->SetLogy();

          TH1F* uframeB = new TH1F(Form("uframeB_%s", ueTag.Data()), "", 100, GetPlotPtMin(), GetPlotPtMax());
          uframeB->SetDirectory(nullptr); uframeB->SetStats(0);
          uframeB->SetMinimum(1e-6); uframeB->SetMaximum(1e2);
          uframeB->GetXaxis()->SetTitle(""); uframeB->GetXaxis()->SetLabelSize(0);
          uframeB->GetYaxis()->SetTitle("d^{2}#sigma/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]");
          uframeB->GetXaxis()->SetTitleFont(43); uframeB->GetYaxis()->SetTitleFont(43);
          uframeB->GetXaxis()->SetLabelFont(43); uframeB->GetYaxis()->SetLabelFont(43);
          uframeB->GetXaxis()->SetTitleSize(uTitSz); uframeB->GetYaxis()->SetTitleSize(uTitSz);
          uframeB->GetYaxis()->SetLabelSize(uLabSz);
          uframeB->GetYaxis()->SetTitleOffset(2.0);
          uframeB->GetXaxis()->SetNdivisions(510); uframeB->GetYaxis()->SetNdivisions(505);
          uframeB->Draw("AXIS");

          // MC truth (dashed)
          for (size_t i = 0; i < ueMCTruthHists.size(); ++i) {
            if (!ueMCTruthHists[i]) continue;
            Color_t col = GetColorForR(ueRvals[i]);
            ueMCTruthHists[i]->SetLineColor(col); ueMCTruthHists[i]->SetMarkerColor(col);
            ueMCTruthHists[i]->SetLineStyle(7); ueMCTruthHists[i]->SetLineWidth(2);
            ueMCTruthHists[i]->SetMarkerStyle(0); ueMCTruthHists[i]->SetStats(0);
            ueMCTruthHists[i]->Draw("l same");
          }

          // Standalone models
          for (auto& m : ueModels) {
            for (size_t i = 0; i < ueModelHists[m.name.Data()].size(); ++i) {
              TH1* hM = ueModelHists[m.name.Data()][i];
              if (!hM) continue;
              Color_t col = GetColorForR(ueRvals[i]);
              hM->SetLineColor(col); hM->SetMarkerColor(col);
              hM->SetLineStyle(m.lineStyle); hM->SetLineWidth(m.lineWidth);
              hM->SetMarkerStyle(0);
              hM->Draw("l same");
            }
          }

          // POWHEG
          for (size_t i = 0; i < uePowhegHists.size(); ++i) {
            if (uePowhegHists[i]) {
              Color_t col = GetColorForR(ueRvals[i]);
              uePowhegHists[i]->SetLineColor(col); uePowhegHists[i]->SetMarkerColor(col);
              uePowhegHists[i]->SetLineStyle(3); uePowhegHists[i]->SetLineWidth(3);
              uePowhegHists[i]->SetMarkerStyle(0);
              uePowhegHists[i]->Draw("l same");
            }
          }

          // Run3 data (on top)
          for (size_t i = 0; i < ueHists.size(); ++i) {
            if (ueHists[i]) {
              Color_t col = GetColorForR(ueRvals[i]);
              ueHists[i]->SetMarkerColor(col); ueHists[i]->SetLineColor(col);
              ueHists[i]->SetMarkerStyle(21);
              ueHists[i]->Draw("pe same");
            }
          }

          // Legends
          TLegend* ulegB = new TLegend(0.35, 0.40, 0.62, 0.90, NULL, "brNDC");
          ulegB->SetTextFont(43); ulegB->SetTextSize(uLegSz);
          ulegB->SetBorderSize(0); ulegB->SetFillStyle(0);
          TH1F* udB_r3 = new TH1F(Form("udB_r3_%s", ueTag.Data()), "", 1, 0, 1);
          udB_r3->SetLineColor(kBlack); udB_r3->SetMarkerStyle(21);
          ulegB->AddEntry(udB_r3, "Run 3", "lpe");
          TH1F* udB_mct = new TH1F(Form("udB_mct_%s", ueTag.Data()), "", 1, 0, 1);
          udB_mct->SetLineColor(kBlack); udB_mct->SetLineStyle(7); udB_mct->SetLineWidth(2);
          ulegB->AddEntry(udB_mct, "MC particle-level", "l");
          for (auto& m : ueModels) {
            TH1F* udB_m = new TH1F(Form("udB_%s_%s", m.name.Data(), ueTag.Data()), "", 1, 0, 1);
            udB_m->SetLineColor(kBlack); udB_m->SetLineStyle(m.lineStyle); udB_m->SetLineWidth(m.lineWidth);
            ulegB->AddEntry(udB_m, m.legLabel.Data(), "l");
          }
          TH1F* udB_pw = new TH1F(Form("udB_pw_%s", ueTag.Data()), "", 1, 0, 1);
          udB_pw->SetLineColor(kBlack); udB_pw->SetLineStyle(3); udB_pw->SetLineWidth(3);
          ulegB->AddEntry(udB_pw, "POWHEG NLO", "l");
          ulegB->Draw();

          TLegend* uRlegB = new TLegend(0.70, 0.52, 0.92, 0.90, NULL, "brNDC");
          uRlegB->SetTextFont(43); uRlegB->SetTextSize(uLegSz);
          uRlegB->SetBorderSize(0); uRlegB->SetFillStyle(0);
          for (size_t i = 0; i < ueRvals.size(); ++i) {
            TH1F* udR = new TH1F(Form("udRB_%s_%d", ueTag.Data(), (int)i), "", 1, 0, 1);
            udR->SetLineColor(GetColorForR(ueRvals[i])); udR->SetMarkerColor(GetColorForR(ueRvals[i]));
            udR->SetMarkerStyle(21); udR->SetLineWidth(2);
            uRlegB->AddEntry(udR, Form("#it{R} = %.1f", ueRvals[i]), "lpe");
          }
          uRlegB->Draw();
          ALICEfigureLegend("ALICE WIP", 0.2, 0.70, 0.4, 0.95, 0.2, 0.05, 0.4, 0.15, 0.04);
          AddUESubtractedLabel(mainPadUB, 0.70, 0.30, 0.045);

          // Ratio pads: [0]=MC-truth, [1..nModels-1]=standalone, [nModels]=POWHEG
          for (int mp = 0; mp < ueNRatioPads; ++mp) {
            uePadsB[mp]->cd();
            gPad->SetGrid(0);

            bool isMCTruth = (mp == 0);
            bool isPowheg = (mp == ueNRatioPads - 1);
            int modelIdx = mp - 1;
            TString mLabel = isMCTruth ? "MC part-level" : (isPowheg ? "POWHEG NLO" : ueModels[modelIdx].legLabel);
            TString yTitle = Form("%s / Data", mLabel.Data());
            bool isLast = (mp == ueNRatioPads - 1);

            double yMin = (isMCTruth || isPowheg) ? 0.5 : 0.5;
            double yMax = (isMCTruth || isPowheg) ? 2.0 : 3.5;

            TH1F* umframe = new TH1F(Form("umframeB_%s_%d", ueTag.Data(), mp), "", 100, GetPlotPtMin(), GetPlotPtMax());
            umframe->SetDirectory(nullptr); umframe->SetStats(0);
            umframe->SetMinimum(yMin); umframe->SetMaximum(yMax);
            umframe->GetXaxis()->SetTitle(isLast ? GetJetPtGenTitleX() : "");
            umframe->GetYaxis()->SetTitle(yTitle.Data());
            umframe->GetXaxis()->SetTitleFont(43); umframe->GetYaxis()->SetTitleFont(43);
            umframe->GetXaxis()->SetLabelFont(43); umframe->GetYaxis()->SetLabelFont(43);
            umframe->GetXaxis()->SetTitleSize(uTitSz); umframe->GetYaxis()->SetTitleSize(uTitSz);
            umframe->GetXaxis()->SetLabelSize(isLast ? uLabSz : 0);
            umframe->GetYaxis()->SetLabelSize(uLabSz);
            umframe->GetYaxis()->SetTitleOffset(2.0);
            if (isLast) umframe->GetXaxis()->SetTitleOffset(3.5);
            umframe->GetXaxis()->SetNdivisions(510); umframe->GetYaxis()->SetNdivisions(505);
            umframe->Draw("AXIS");

            for (size_t i = 0; i < ueHists.size(); ++i) {
              if (!ueHists[i]) continue;
              TH1* hModel = nullptr;
              if (isMCTruth) {
                hModel = (i < ueMCTruthHists.size()) ? ueMCTruthHists[i] : nullptr;
              } else if (isPowheg) {
                hModel = (i < uePowhegHists.size()) ? uePowhegHists[i] : nullptr;
              } else {
                auto& mhists = ueModelHists[ueModels[modelIdx].name.Data()];
                hModel = (i < mhists.size()) ? mhists[i] : nullptr;
              }
              if (!hModel) continue;
              TH1* hRatUE = (TH1*)hModel->Clone(Form("ueRatB_%d_%s_%zu", mp, ueTag.Data(), i));
              hRatUE->Divide(hModel, ueHists[i], 1., 1., "");
              hRatUE->SetStats(0);
              Color_t col = GetColorForR(ueRvals[i]);
              hRatUE->SetLineColor(col); hRatUE->SetMarkerColor(col);
              hRatUE->SetLineStyle(1); hRatUE->SetLineWidth(2); hRatUE->SetMarkerStyle(0);
              hRatUE->Draw("l same");
            }

            TLine* umLine = new TLine(GetPlotPtMin(), 1.0, GetPlotPtMax(), 1.0);
            umLine->SetLineColor(kBlack); umLine->SetLineStyle(2); umLine->SetLineWidth(1);
            umLine->Draw("lsame");
          }

          if (GetDRAWPLOTS()) {
            canUB->Update(); canUB->Modified();
            canUB->Print(Form("%s/CrossSection_AllR_UE_Models_%s.pdf", outputDir.Data(), ueTag.Data()));
          }
        }
      }
    }

    return;
  }

  // ========== Canvas A: CrossSection_AllR_Run2 — Run3 vs Run2 only ==========
  {
    const int arTitSz = 18, arLabSz = 14, arLegSz = 13;
    double arL = 0.15, arR = 0.03;
    double mainFracA = 0.65, ratioFracA = 0.35;

    TCanvas* canA = new TCanvas(Form("CrossSection_AllR_Run2_%d", ++GetNN()),
                                "Cross Section AllR: Run3 vs Run2", 800, 700);
    canA->SetFillStyle(4000); canA->SetFillColor(10);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    canA->SetTopMargin(0.); canA->SetBottomMargin(0.);
    canA->SetLeftMargin(0.); canA->SetRightMargin(0.);
    canA->Draw();

    TPad* mainPadA = new TPad("arMainA", "Main", 0, ratioFracA, 1, 1.0, 0);
    mainPadA->SetTopMargin(0.04 / mainFracA);
    mainPadA->SetBottomMargin(0.001);
    mainPadA->SetLeftMargin(arL); mainPadA->SetRightMargin(arR);
    mainPadA->Draw();

    TPad* ratioPadA = new TPad("arRatioA", "Run3/Run2", 0, 0, 1, ratioFracA, 0);
    ratioPadA->SetTopMargin(0.001);
    ratioPadA->SetBottomMargin(0.12 / ratioFracA);
    ratioPadA->SetLeftMargin(arL); ratioPadA->SetRightMargin(arR);
    ratioPadA->Draw();

    optFili(*mainPadA, 0, 1, 0, 1);
    optFili(*ratioPadA, 0, 0, 0, 0);

    // Main pad: Run3 + Run2 spectra
    mainPadA->cd();
    gPad->SetGrid(0);

    TLegend* legA = new TLegend(0.50, 0.55, 0.92, 0.95, NULL, "brNDC");
    legA->SetTextFont(43); legA->SetTextSize(arLegSz);
    legA->SetBorderSize(0);

    bool firstA = true;
    TH1* hFirstA = (!allUnfoldedHists.empty() && allUnfoldedHists[0]) ? allUnfoldedHists[0] : nullptr;

    // Run2 systematic boxes
    for (size_t i = 0; i < allRun2Graphs.size(); ++i) {
      if (allRun2Graphs[i].second) {
        if (firstA && hFirstA) {
          hset(*hFirstA, "", "d^{2}#sigma/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
          hFirstA->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
          hFirstA->GetYaxis()->SetRangeUser(1e-6, 1e2);
          hFirstA->GetXaxis()->SetLabelSize(0); hFirstA->GetXaxis()->SetTitleSize(0);
          hFirstA->GetYaxis()->SetTitleFont(43); hFirstA->GetYaxis()->SetTitleSize(arTitSz);
          hFirstA->GetYaxis()->SetLabelFont(43); hFirstA->GetYaxis()->SetLabelSize(arLabSz);
          hFirstA->GetYaxis()->SetTitleOffset(2.5);
          hFirstA->Draw("AXIS"); firstA = false;
        }
        allRun2Graphs[i].second->Draw(firstA ? "E2" : "E2 same");
        if (firstA) firstA = false;
      }
    }
    for (size_t i = 0; i < allRun2Graphs.size(); ++i) {
      if (allRun2Graphs[i].first) {
        if (firstA) {
          allRun2Graphs[i].first->GetYaxis()->SetRangeUser(1e-6, 1e2);
          allRun2Graphs[i].first->Draw("APZ"); firstA = false;
        } else {
          allRun2Graphs[i].first->Draw("PZ same");
        }
        if (i == 0) legA->AddEntry(allRun2Graphs[i].first, "Run 2 data", "pe");
      }
    }
    for (size_t i = 0; i < allUnfoldedHists.size(); ++i) {
      if (allUnfoldedHists[i]) {
        allUnfoldedHists[i]->Draw(firstA ? "pe" : "pe same");
        if (firstA) firstA = false;
        legA->AddEntry(allUnfoldedHists[i], rConfigs[i].label.Data(), "pe");
      }
    }
    legA->Draw();
    ALICEfigureLegend("ALICE WIP", 0.2, 0.70, 0.4, 0.95, 0.2, 0.05, 0.4, 0.15, 0.04);

    // Ratio pad: Run3/Run2
    ratioPadA->cd();
    gPad->SetGrid(0);
    TH1F* frameA = new TH1F(Form("arFrameA_%d", GetNN()), "", 100, GetPlotPtMin(), GetPlotPtMax());
    frameA->SetDirectory(nullptr); frameA->SetStats(0);
    frameA->SetMinimum(0.5); frameA->SetMaximum(2.0);
    frameA->GetXaxis()->SetTitle(GetJetPtGenTitleX());
    frameA->GetYaxis()->SetTitle("Run 3 / Run 2");
    frameA->GetXaxis()->SetTitleFont(43); frameA->GetYaxis()->SetTitleFont(43);
    frameA->GetXaxis()->SetLabelFont(43); frameA->GetYaxis()->SetLabelFont(43);
    frameA->GetXaxis()->SetTitleSize(arTitSz); frameA->GetYaxis()->SetTitleSize(arTitSz);
    frameA->GetXaxis()->SetLabelSize(arLabSz); frameA->GetYaxis()->SetLabelSize(arLabSz);
    frameA->GetYaxis()->SetTitleOffset(2.5);
    frameA->GetXaxis()->SetTitleOffset(3.0);
    frameA->GetXaxis()->SetNdivisions(510); frameA->GetYaxis()->SetNdivisions(505);
    frameA->GetXaxis()->CenterTitle(1); frameA->GetYaxis()->CenterTitle(1);
    frameA->Draw("AXIS");

    for (size_t i = 0; i < allUnfoldedHists.size(); ++i) {
      if (i < allRun2Graphs.size() && allRun2Graphs[i].first && rConfigs[i].R >= 0.2) {
        TH1* hRatio1 = CreateRatioTH1vsTGraph(allUnfoldedHists[i], allRun2Graphs[i].first, Form("ratio_Run3_Run2_AllR_%zu", i));
        if (hRatio1) {
          hRatio1->SetStats(0);
          hRatio1->SetLineColor(rColors[i]); hRatio1->SetMarkerColor(rColors[i]);
          hRatio1->SetMarkerStyle(20); hRatio1->SetLineWidth(2);
          hRatio1->Draw("pe same");
        }
      }
    }
    TLine* lineA = new TLine(GetPlotPtMin(), 1.0, GetPlotPtMax(), 1.0);
    lineA->SetLineColor(kBlack); lineA->SetLineStyle(2); lineA->SetLineWidth(1);
    lineA->Draw("lsame");

    if (GetDRAWPLOTS()) {
      canA->Update(); canA->Modified();
      canA->Print(Form("%s/CrossSection_AllR_Run2.pdf", outputDir.Data()));
    }
  }

  // ========== Canvas B: CrossSection_AllR_Models — Run3 vs MC-truth vs Models ==========
  {
    auto allModels = GetStandaloneModels();
    // Ratio pads: MC-truth + standalone models + POWHEG
    int nModelPads = (int)allModels.size() + 2;  // +1 for MC-truth, +1 for POWHEG
    const int arTitSz = 18, arLabSz = 14, arLegSz = 13;
    double arL = 0.15, arR = 0.03;
    double mainFracB = 0.35;
    double modelFracB = (1.0 - mainFracB) / nModelPads;

    TCanvas* canB = new TCanvas(Form("CrossSection_AllR_Models_%d", ++GetNN()),
                                "Cross Section AllR: Run3 vs Models", 800, 1600);
    canB->SetFillStyle(4000); canB->SetFillColor(10);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    canB->SetTopMargin(0.); canB->SetBottomMargin(0.);
    canB->SetLeftMargin(0.); canB->SetRightMargin(0.);
    canB->Draw();

    // Main pad
    TPad* mainPadB = new TPad("arMainB", "Main", 0, 1.0 - mainFracB, 1, 1.0, 0);
    mainPadB->SetTopMargin(0.04 / mainFracB);
    mainPadB->SetBottomMargin(0.001);
    mainPadB->SetLeftMargin(arL); mainPadB->SetRightMargin(arR);
    mainPadB->Draw();

    // Ratio pads: [0]=MC-truth, [1..nModels-1]=standalone models, [nModels]=POWHEG
    std::vector<TPad*> modelPadsB;
    for (int mp = 0; mp < nModelPads; ++mp) {
      double top = (1.0 - mainFracB) - mp * modelFracB;
      double bot = top - modelFracB;
      bool isLast = (mp == nModelPads - 1);
      TString mName;
      if (mp == 0) mName = "MCTruth";
      else if (mp <= (int)allModels.size()) mName = allModels[mp - 1].name;
      else mName = "POWHEG";
      TPad* pad = new TPad(Form("arModelB_%s", mName.Data()), mName.Data(), 0, bot, 1, top, 0);
      pad->SetTopMargin(0.001);
      pad->SetBottomMargin(isLast ? 0.12 / modelFracB : 0.001);
      pad->SetLeftMargin(arL); pad->SetRightMargin(arR);
      pad->Draw();
      modelPadsB.push_back(pad);
    }

    optFili(*mainPadB, 0, 1, 0, 1);
    for (auto* p : modelPadsB) optFili(*p, 0, 0, 0, 0);

    // Main pad: yield with models
    mainPadB->cd();
    gPad->SetGrid(0);

    TLegend* legB = new TLegend(0.40, 0.35, 0.92, 0.95, NULL, "brNDC");
    legB->SetTextFont(43); legB->SetTextSize(arLegSz);
    legB->SetBorderSize(0);

    bool firstB = true;
    TH1* hFirstB = (!allUnfoldedHists.empty() && allUnfoldedHists[0]) ? allUnfoldedHists[0] : nullptr;

    auto drawFrameB = [&]() {
      if (firstB && hFirstB) {
        hset(*hFirstB, "", "d^{2}#sigma/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
        hFirstB->GetXaxis()->SetRangeUser(GetPlotPtMin(), GetPlotPtMax());
        hFirstB->GetYaxis()->SetRangeUser(1e-6, 1e2);
        hFirstB->GetXaxis()->SetLabelSize(0); hFirstB->GetXaxis()->SetTitleSize(0);
        hFirstB->GetYaxis()->SetTitleFont(43); hFirstB->GetYaxis()->SetTitleSize(arTitSz);
        hFirstB->GetYaxis()->SetLabelFont(43); hFirstB->GetYaxis()->SetLabelSize(arLabSz);
        hFirstB->GetYaxis()->SetTitleOffset(2.5);
        hFirstB->Draw("AXIS"); firstB = false;
      }
    };

    // MC truth lines (dashed gray)
    for (size_t i = 0; i < allMCTruthHists.size(); ++i) {
      if (!allMCTruthHists[i]) continue;
      drawFrameB();
      allMCTruthHists[i]->SetLineColor(rColors[i]);
      allMCTruthHists[i]->SetMarkerColor(rColors[i]);
      allMCTruthHists[i]->SetLineStyle(7);  // dashed
      allMCTruthHists[i]->SetLineWidth(2);
      allMCTruthHists[i]->SetMarkerStyle(0);
      allMCTruthHists[i]->SetStats(0);
      allMCTruthHists[i]->Draw("l same");
      if (i == 0) {
        TH1F* dumMCT = new TH1F(Form("dumMCT_%d", GetNN()), "", 1, 0, 1);
        dumMCT->SetLineColor(kBlack); dumMCT->SetLineStyle(7); dumMCT->SetLineWidth(2);
        legB->AddEntry(dumMCT, "MC particle-level", "l");
      }
    }

    // Standalone models
    for (auto& m : allModels) {
      if (!allModelHists.count(m.name.Data())) continue;
      auto& hists = allModelHists[m.name.Data()];
      bool legendAdded = false;
      for (size_t i = 0; i < hists.size(); ++i) {
        if (!hists[i]) continue;
        drawFrameB();
        hists[i]->SetLineStyle(m.lineStyle);
        hists[i]->SetLineWidth(m.lineWidth);
        hists[i]->SetMarkerStyle(0);
        hists[i]->Draw("l same");
        if (!legendAdded) { legB->AddEntry(hists[i], m.legLabel.Data(), "l"); legendAdded = true; }
      }
    }

    // POWHEG NLO
    for (size_t i = 0; i < allPowhegHists.size(); ++i) {
      if (allPowhegHists[i]) {
        drawFrameB();
        allPowhegHists[i]->SetLineStyle(3);
        allPowhegHists[i]->SetLineWidth(3);
        allPowhegHists[i]->Draw("l same");
        if (i == 0) legB->AddEntry(allPowhegHists[i], "POWHEG NLO", "l");
      }
    }

    // Run3 unfolded data (on top)
    for (size_t i = 0; i < allUnfoldedHists.size(); ++i) {
      if (allUnfoldedHists[i]) {
        drawFrameB();
        allUnfoldedHists[i]->Draw("pe same");
        legB->AddEntry(allUnfoldedHists[i], rConfigs[i].label.Data(), "pe");
      }
    }
    legB->Draw();
    ALICEfigureLegend("ALICE WIP", 0.2, 0.70, 0.4, 0.95, 0.2, 0.05, 0.4, 0.15, 0.04);

    // Ratio pads
    for (int mp = 0; mp < nModelPads; ++mp) {
      modelPadsB[mp]->cd();
      gPad->SetGrid(0);

      bool isMCTruth = (mp == 0);
      bool isPowheg = (mp == nModelPads - 1);
      int modelIdx = mp - 1;  // index into allModels (valid when !isMCTruth && !isPowheg)
      TString mLabel = isMCTruth ? "MC part-level" : (isPowheg ? "POWHEG NLO" : allModels[modelIdx].legLabel);
      TString yTitle = Form("%s / Data", mLabel.Data());
      bool isLast = (mp == nModelPads - 1);

      // MC-truth and POWHEG get tighter range; standalone models get wider range
      double yMin = (isMCTruth || isPowheg) ? 0.5 : 0.5;
      double yMax = (isMCTruth || isPowheg) ? 2.0 : 3.5;

      TH1F* mframe = new TH1F(Form("arMFrameB_%d_%d", mp, GetNN()), "", 100, GetPlotPtMin(), GetPlotPtMax());
      mframe->SetDirectory(nullptr); mframe->SetStats(0);
      mframe->SetMinimum(yMin); mframe->SetMaximum(yMax);
      mframe->GetXaxis()->SetTitle(isLast ? GetJetPtGenTitleX() : "");
      mframe->GetYaxis()->SetTitle(yTitle.Data());
      mframe->GetXaxis()->SetTitleFont(43); mframe->GetYaxis()->SetTitleFont(43);
      mframe->GetXaxis()->SetLabelFont(43); mframe->GetYaxis()->SetLabelFont(43);
      mframe->GetXaxis()->SetTitleSize(arTitSz); mframe->GetYaxis()->SetTitleSize(arTitSz);
      mframe->GetXaxis()->SetLabelSize(isLast ? arLabSz : 0);
      mframe->GetYaxis()->SetLabelSize(arLabSz);
      mframe->GetYaxis()->SetTitleOffset(2.5);
      if (isLast) mframe->GetXaxis()->SetTitleOffset(3.5);
      mframe->GetXaxis()->SetNdivisions(510); mframe->GetYaxis()->SetNdivisions(505);
      mframe->GetXaxis()->CenterTitle(1); mframe->GetYaxis()->CenterTitle(1);
      mframe->Draw("AXIS");

      for (size_t i = 0; i < allUnfoldedHists.size(); ++i) {
        if (!allUnfoldedHists[i]) continue;
        TH1* hModel = nullptr;
        if (isMCTruth) {
          hModel = (i < allMCTruthHists.size()) ? allMCTruthHists[i] : nullptr;
        } else if (isPowheg) {
          hModel = (i < allPowhegHists.size()) ? allPowhegHists[i] : nullptr;
        } else {
          auto& mhists = allModelHists[allModels[modelIdx].name.Data()];
          hModel = (i < mhists.size()) ? mhists[i] : nullptr;
        }
        if (!hModel) continue;
        TH1* hRat = (TH1*)hModel->Clone(Form("arRatB_%d_%zu_%d", mp, i, GetNN()));
        hRat->Divide(hModel, allUnfoldedHists[i], 1., 1., "");
        hRat->SetStats(0);
        hRat->SetLineColor(rColors[i]); hRat->SetMarkerColor(rColors[i]);
        hRat->SetLineStyle(1); hRat->SetLineWidth(2); hRat->SetMarkerStyle(0);
        hRat->Draw("l same");
      }

      TLine* mLine = new TLine(GetPlotPtMin(), 1.0, GetPlotPtMax(), 1.0);
      mLine->SetLineColor(kBlack); mLine->SetLineStyle(2); mLine->SetLineWidth(1);
      mLine->Draw("lsame");
    }

    if (GetDRAWPLOTS()) {
      canB->Update(); canB->Modified();
      canB->Print(Form("%s/CrossSection_AllR_Models.pdf", outputDir.Data()));
    }
  }

  // R-ratio from NonUE unfolded spectra
  if (!allUnfoldedHists.empty()) {
    std::map<double, TH1*> spectraByR;
    for (size_t i = 0; i < allUnfoldedHists.size(); ++i) {
      if (allUnfoldedHists[i]) {
        spectraByR[rConfigs[i].R] = allUnfoldedHists[i];
      }
    }
    DrawRRatioPlots(spectraByR, rColors, outputDir, "NonUE");
  }
}

// ============================================
// Save Data Results for standalone comparison
// ============================================

void SaveDataResults() {
  std::vector<RConfig> rConfigs = GetRConfigs();
  TString outputDir = GetOutputDir();
  EXsecMode mode = GetCurrentMode();
  TString suffix = (mode == kUE) ? "_UE" : "";

  // Load cross-sections from cache
  TString cachePath = Form("%s/XsecCache.root", outputDir.Data());
  std::vector<TH1*> vecUnfolded, vecPOWHEG;
  std::map<std::string, std::vector<TH1*>> vecModels;
  if (!LoadXsecCache(cachePath, kXsecCacheVersion, mode, rConfigs, vecUnfolded, vecPOWHEG, vecModels)) {
    std::cerr << "[SaveDataResults] No XsecCache found. Skipping." << std::endl;
    return;
  }

  // Compute invariant yield per R (from cached unfolding)
  const char* regFile = (mode == kUE) ? kOptimalRegFileUE : kOptimalRegFileNonUE;
  auto regMap = LoadOptimalRegularization(regFile);

  std::vector<TH1*> vecInvYield;
  for (int iR = 0; iR < (int)rConfigs.size(); ++iR) {
    const RConfig& config = rConfigs[iR];
    Filipad2* tempPad = new Filipad2(Form("tempSave_R%.1f", config.R), ++GetNN(), 2, 0.3, 100, 50, 0.7, 1, 1);
    tempPad->Draw();
    TPad* uf = tempPad->GetPad(1);
    TPad* rf = tempPad->GetPad(2);
    TLegend* tl = new TLegend(0.75, 0.3, 0.95, 0.65, NULL, "brNDC");

    Double_t nD  = GetNORMEVENTS() ? Nevents(config.GetDataFile().Data(), config.dataDir.Data(), GetEventObj(), 0) : 1.0;
    Double_t nDU = GetNORMEVENTS() ? Nevents(config.GetDataFile().Data(), config.dataDir.Data(), GetEventObj(), 1) : 1.0;
    Double_t nMD = GetNORMEVENTS() ? Nevents(config.GetMcFile().Data(), config.mcDir.Data(), GetEventObj(), 0, config.isJJ) : 1.0;
    Double_t nMP = GetNORMEVENTS() ? Nevents(config.GetMcFile().Data(), config.mcDir.Data(), GetEventObj(), 1, config.isJJ) : 1.0;
    const Int_t kOvr = GetSvdKForRun(config.mcRunNumber, regMap);

    TH1* hDcorr = DrawJetMatching(config.GetMcFile().Data(), config.GetDataFile().Data(), config.label.Data(),
                   nD, nDU, nMD, nMP, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                   uf, rf, tl, kBlack, nullptr, nullptr, nullptr, nullptr,
                   config.dataDir.Data(), config.mcDir.Data(), "", nullptr, nullptr, nullptr, kOvr);
    if (tempPad->C) tempPad->C->Close();

    if (hDcorr) {
      TH1* hY = (TH1*)hDcorr->Clone(Form("hDataInvYield_R%02d%s", (int)(config.R*10+0.5), suffix.Data()));
      hY->SetDirectory(0);
      Double_t dEta = GetDeltaEta(config.R);
      if (nD > 0) hY->Scale(1.0 / nD, "width");
      if (dEta > 0) hY->Scale(1.0 / dEta);
      vecInvYield.push_back(hY);
    } else {
      vecInvYield.push_back(nullptr);
    }
  }

  // Save to DataResults file
  TString outPath = Form("%s/DataResults%s.root", outputDir.Data(), suffix.Data());
  TFile* fout = TFile::Open(outPath, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "[SaveDataResults] Cannot create " << outPath << std::endl;
    if (fout) delete fout;
    return;
  }

  // Metadata
  TParameter<int> pMode("mode", (int)mode); pMode.Write();
  TParameter<int> pNR("nR", (int)rConfigs.size()); pNR.Write();

  // Per-R: cross-section, invariant yield
  for (int iR = 0; iR < (int)rConfigs.size(); ++iR) {
    int iRtag = (int)(rConfigs[iR].R * 10 + 0.5);
    TString tag = Form("R%02d", iRtag);

    TParameter<double> pR(Form("R_%d", iR), rConfigs[iR].R); pR.Write();

    if (iR < (int)vecUnfolded.size() && vecUnfolded[iR]) {
      vecUnfolded[iR]->Write(Form("hDataXsec_%s", tag.Data()));
    }
    if (iR < (int)vecInvYield.size() && vecInvYield[iR]) {
      vecInvYield[iR]->Write(Form("hDataInvYield_%s", tag.Data()));
    }
  }

  // R-ratios (all pairs: smaller R / larger R)
  for (int i = 0; i < (int)rConfigs.size(); ++i) {
    for (int j = i + 1; j < (int)rConfigs.size(); ++j) {
      if (!vecUnfolded[i] || !vecUnfolded[j]) continue;
      int iRnum = (int)(rConfigs[i].R * 10 + 0.5);
      int iRden = (int)(rConfigs[j].R * 10 + 0.5);
      TH1* hRat = (TH1*)vecUnfolded[i]->Clone(Form("hDataXsecRatio_R%02d_R%02d", iRnum, iRden));
      hRat->SetDirectory(0);
      // Account for different Δη: ratio = (σ_i/Δη_i) / (σ_j/Δη_j)
      // vecUnfolded already has d²σ/(dpT dη), so just divide
      hRat->Divide(vecUnfolded[i], vecUnfolded[j], 1., 1., "");
      hRat->Write();
      delete hRat;
    }
  }

  fout->Close();
  delete fout;
  std::cerr << "[SaveDataResults] Saved to " << outPath << std::endl;
}

// ============================================
// Model Comparison: per-model Model/Data ratio
// ============================================

void DrawModelComparison() {
  if (!GetDrawCrossSection()) return;

  std::vector<RConfig> rConfigs = GetRConfigs();
  TString outputDir = GetOutputDir();
  EXsecMode mode = GetCurrentMode();
  TString suffix = (mode == kUE) ? "_UE" : "";
  const Double_t kSigmaINEL = 79.0;  // mb

  std::vector<StandaloneModel> models = GetStandaloneModels();

  // Load cross-section data from cache
  TString cachePath = Form("%s/XsecCache.root", outputDir.Data());
  std::vector<TH1*> vecUnfolded, vecPOWHEG;
  std::map<std::string, std::vector<TH1*>> vecModels;

  if (!LoadXsecCache(cachePath, kXsecCacheVersion, mode, rConfigs, vecUnfolded, vecPOWHEG, vecModels)) {
    std::cerr << "[DrawModelComparison] No XsecCache found. Run DrawCrossSection() first." << std::endl;
    return;
  }

  int nR = (int)rConfigs.size();
  int canW = 400 * nR;
  int canH = 500;

  // Load regularization for unfolding (needed for invariant yield data)
  const char* regFile = (mode == kUE) ? kOptimalRegFileUE : kOptimalRegFileNonUE;
  auto regMap = LoadOptimalRegularization(regFile);

  // ====================================================================
  // Pre-compute data invariant yield per R (uses cached unfolding)
  // ====================================================================
  std::vector<TH1*> vecDataYield;
  for (int iR = 0; iR < nR; ++iR) {
    const RConfig& config = rConfigs[iR];

    // Unfolding (uses gUnfoldCache if already computed by DrawInvariantYield)
    Filipad2* tempPad = new Filipad2(Form("tempUnfold_MC_R%.1f", config.R), ++GetNN(), 2, 0.3, 100, 50, 0.7, 1, 1);
    tempPad->Draw();
    TPad* unfoldpad = tempPad->GetPad(1);
    TPad* ratunfoldpad = tempPad->GetPad(2);
    TLegend* tempLeg = new TLegend(0.75, 0.3, 0.95, 0.65, NULL, "brNDC");

    Double_t nevtsData = GetNORMEVENTS() ? Nevents(config.GetDataFile().Data(), config.dataDir.Data(), GetEventObj(), 0) : 1.0;
    Double_t nevtsDataUnTrig = GetNORMEVENTS() ? Nevents(config.GetDataFile().Data(), config.dataDir.Data(), GetEventObj(), 1) : 1.0;
    Double_t nevtsMCD = GetNORMEVENTS() ? Nevents(config.GetMcFile().Data(), config.mcDir.Data(), GetEventObj(), 0, config.isJJ) : 1.0;
    Double_t nevtsMCP = GetNORMEVENTS() ? Nevents(config.GetMcFile().Data(), config.mcDir.Data(), GetEventObj(), 1, config.isJJ) : 1.0;
    const Int_t svdKOverride = GetSvdKForRun(config.mcRunNumber, regMap);

    TH1* hDcorrected = DrawJetMatching(config.GetMcFile().Data(), config.GetDataFile().Data(), config.label.Data(),
                   nevtsData, nevtsDataUnTrig, nevtsMCD, nevtsMCP,
                   nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                   unfoldpad, ratunfoldpad, tempLeg, kBlack,
                   nullptr, nullptr, nullptr, nullptr, config.dataDir.Data(),
                   config.mcDir.Data(), "",
                   nullptr, nullptr, nullptr,
                   svdKOverride);

    if (tempPad->C) tempPad->C->Close();

    if (hDcorrected) {
      TH1* hYield = (TH1*)hDcorrected->Clone(Form("invYieldData_MC_%d_%d", iR, GetNN()));
      hYield->SetDirectory(0);
      Double_t deltaEta = GetDeltaEta(config.R);
      if (nevtsData > 0) hYield->Scale(1.0 / nevtsData, "width");
      if (deltaEta > 0) hYield->Scale(1.0 / deltaEta);
      vecDataYield.push_back(hYield);
    } else {
      vecDataYield.push_back(nullptr);
    }
  }

  // ====================================================================
  // Canvas A: Invariant yield per-R Model/Data ratios
  // ====================================================================
  TCanvas* canA = new TCanvas(Form("ModelComp_InvY_%d", ++GetNN()), "Model Comparison InvYield", canW, canH);
  canA->Divide(nR, 1, 0.001, 0.001);

  for (int iR = 0; iR < nR; ++iR) {
    canA->cd(iR + 1);
    gPad->SetGrid(0);
    gPad->SetLeftMargin(iR == 0 ? 0.18 : 0.05);
    gPad->SetRightMargin(iR == nR - 1 ? 0.05 : 0.01);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);

    TH1* hData = vecDataYield[iR];
    double R = rConfigs[iR].R;

    TH1F* frame = new TH1F(Form("frameMC_IY_%d_%d", iR, GetNN()), "", 100, 5, 200);
    frame->SetDirectory(nullptr);
    frame->SetStats(0);
    frame->SetMinimum(0.5);
    frame->SetMaximum(3.5);
    frame->GetXaxis()->SetTitle(GetJetPtGenTitleX());
    frame->GetYaxis()->SetTitle(iR == 0 ? "Model / Data" : "");
    frame->GetXaxis()->SetTitleSize(0.06);
    frame->GetYaxis()->SetTitleSize(0.06);
    frame->GetXaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetLabelSize(iR == 0 ? 0.05 : 0.0);
    frame->GetXaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetTitleOffset(1.3);
    frame->GetXaxis()->SetNdivisions(510);
    frame->GetYaxis()->SetNdivisions(505);
    frame->Draw("AXIS");

    TLatex* latR = new TLatex(0.5, 0.92, Form("#it{R} = %.1f", R));
    latR->SetNDC(); latR->SetTextSize(0.06); latR->SetTextAlign(22); latR->Draw();

    TLine* lineA = new TLine(5, 1.0, 200, 1.0);
    lineA->SetLineColor(kBlack); lineA->SetLineStyle(2); lineA->SetLineWidth(1);
    lineA->Draw("lsame");

    if (!hData) continue;

    // MC truth ratio (invariant yield)
    TH1* hMCTruth = LoadMCTruthNormYield(rConfigs[iR]);
    if (hMCTruth) {
      TH1* hRatMC = (TH1*)hMCTruth->Clone(Form("ratMCTruthIY_%d_%d", iR, GetNN()));
      hRatMC->Divide(hMCTruth, hData, 1., 1., "");
      hRatMC->SetStats(0);
      hRatMC->SetLineColor(kGray+2); hRatMC->SetLineStyle(7);
      hRatMC->SetLineWidth(2); hRatMC->SetMarkerStyle(0);
      hRatMC->Draw("l same");
      delete hMCTruth;
    }

    // Standalone model ratios (invariant yield)
    Double_t deltaEta = GetDeltaEta(R);
    for (int im = 0; im < (int)models.size(); ++im) {
      TH1* hModel = LoadStandaloneModelNormYield(models[im], R);
      if (!hModel) continue;
      // LoadStandaloneModelNormYield returns (1/Nevt) dN/dpT (×Δη format)
      // Convert to (1/Nevt) d²N/(dpT dη) by dividing by Δη
      if (deltaEta > 0) hModel->Scale(1.0 / deltaEta);

      TH1* hRat = (TH1*)hModel->Clone(Form("ratModelIY_%d_%d_%d", im, iR, GetNN()));
      hRat->Divide(hModel, hData, 1., 1., "");
      hRat->SetStats(0);
      hRat->SetLineColor(models[im].lineColor);
      hRat->SetLineStyle(models[im].lineStyle);
      hRat->SetLineWidth(models[im].lineWidth);
      hRat->SetMarkerStyle(0);
      hRat->Draw("l same");
      delete hModel;
    }

    // POWHEG ratio (cross-section / σ_INEL → per-INEL-event yield)
    TH1* hPOWHEGxs = LoadPowhegForR(R);
    if (hPOWHEGxs) {
      hPOWHEGxs->Scale(1.0 / kSigmaINEL);  // xsec → yield
      TH1* hRatP = (TH1*)hPOWHEGxs->Clone(Form("ratPOWHEGIY_%d_%d", iR, GetNN()));
      hRatP->Divide(hPOWHEGxs, hData, 1., 1., "");
      hRatP->SetStats(0);
      hRatP->SetLineColor(kRed+1); hRatP->SetLineStyle(3);
      hRatP->SetLineWidth(3); hRatP->SetMarkerStyle(0);
      hRatP->Draw("l same");
      delete hPOWHEGxs;
    }

    // Legend in first panel
    if (iR == 0) {
      TLegend* leg = new TLegend(0.22, 0.45, 0.75, 0.88, NULL, "brNDC");
      leg->SetTextSize(0.04); leg->SetBorderSize(0); leg->SetFillStyle(0);

      TH1F* dumMC = new TH1F(Form("dumMCTIY_%d", GetNN()), "", 1, 0, 1);
      dumMC->SetLineColor(kGray+2); dumMC->SetLineStyle(7); dumMC->SetLineWidth(2);
      leg->AddEntry(dumMC, "MC truth", "l");

      for (int im = 0; im < (int)models.size(); ++im) {
        TH1F* dum = new TH1F(Form("dumMIY_%d_%d", im, GetNN()), "", 1, 0, 1);
        dum->SetLineColor(models[im].lineColor);
        dum->SetLineStyle(models[im].lineStyle);
        dum->SetLineWidth(models[im].lineWidth);
        leg->AddEntry(dum, models[im].legLabel.Data(), "l");
      }

      TH1F* dumP = new TH1F(Form("dumPIY_%d", GetNN()), "", 1, 0, 1);
      dumP->SetLineColor(kRed+1); dumP->SetLineStyle(3); dumP->SetLineWidth(3);
      leg->AddEntry(dumP, "POWHEG NLO", "l");

      leg->Draw();
    }
  }

  TString pdfNameA = Form("%s/ModelComparison_InvYield_AllR%s.pdf", outputDir.Data(), suffix.Data());
  if (GetDRAWPLOTS()) {
    canA->Print(pdfNameA.Data());
  }

  // ====================================================================
  // Canvas B: Cross section per-R Model/Data ratios
  // ====================================================================
  TCanvas* canB = new TCanvas(Form("ModelComp_Xsec_%d", ++GetNN()), "Model Comparison Xsec", canW, canH);
  canB->Divide(nR, 1, 0.001, 0.001);

  for (int iR = 0; iR < nR; ++iR) {
    canB->cd(iR + 1);
    gPad->SetGrid(0);
    gPad->SetLeftMargin(iR == 0 ? 0.18 : 0.05);
    gPad->SetRightMargin(iR == nR - 1 ? 0.05 : 0.01);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);

    TH1* hData = (iR < (int)vecUnfolded.size()) ? vecUnfolded[iR] : nullptr;
    double R = rConfigs[iR].R;

    TH1F* frame = new TH1F(Form("frameMC_XS_%d_%d", iR, GetNN()), "", 100, 5, 200);
    frame->SetDirectory(nullptr);
    frame->SetStats(0);
    frame->SetMinimum(0.5);
    frame->SetMaximum(3.5);
    frame->GetXaxis()->SetTitle(GetJetPtGenTitleX());
    frame->GetYaxis()->SetTitle(iR == 0 ? "Model / Data" : "");
    frame->GetXaxis()->SetTitleSize(0.06);
    frame->GetYaxis()->SetTitleSize(0.06);
    frame->GetXaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetLabelSize(iR == 0 ? 0.05 : 0.0);
    frame->GetXaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetTitleOffset(1.3);
    frame->GetXaxis()->SetNdivisions(510);
    frame->GetYaxis()->SetNdivisions(505);
    frame->Draw("AXIS");

    TLatex* latR = new TLatex(0.5, 0.92, Form("#it{R} = %.1f", R));
    latR->SetNDC(); latR->SetTextSize(0.06); latR->SetTextAlign(22); latR->Draw();

    TLine* lineB = new TLine(5, 1.0, 200, 1.0);
    lineB->SetLineColor(kBlack); lineB->SetLineStyle(2); lineB->SetLineWidth(1);
    lineB->Draw("lsame");

    if (!hData) continue;

    // MC truth ratio (cross-section)
    TH1* hMCTruth = LoadMCTruthXsec(rConfigs[iR]);
    if (hMCTruth) {
      TH1* hRatMC = (TH1*)hMCTruth->Clone(Form("ratMCTruthXS_%d_%d", iR, GetNN()));
      hRatMC->Divide(hMCTruth, hData, 1., 1., "");
      hRatMC->SetStats(0);
      hRatMC->SetLineColor(kGray+2); hRatMC->SetLineStyle(7);
      hRatMC->SetLineWidth(2); hRatMC->SetMarkerStyle(0);
      hRatMC->Draw("l same");
      delete hMCTruth;
    }

    // Standalone model ratios (cross-section from cache)
    for (int im = 0; im < (int)models.size(); ++im) {
      TString mName = models[im].name;
      TH1* hModel = nullptr;
      if (vecModels.count(mName.Data()) && iR < (int)vecModels[mName.Data()].size())
        hModel = vecModels[mName.Data()][iR];
      if (!hModel) continue;

      TH1* hRat = (TH1*)hModel->Clone(Form("ratModelXS_%d_%d_%d", im, iR, GetNN()));
      hRat->Divide(hModel, hData, 1., 1., "");
      hRat->SetStats(0);
      hRat->SetLineColor(models[im].lineColor);
      hRat->SetLineStyle(models[im].lineStyle);
      hRat->SetLineWidth(models[im].lineWidth);
      hRat->SetMarkerStyle(0);
      hRat->Draw("l same");
    }

    // POWHEG ratio (cross-section from cache)
    TH1* hPOWHEG = (iR < (int)vecPOWHEG.size()) ? vecPOWHEG[iR] : nullptr;
    if (hPOWHEG) {
      TH1* hRatP = (TH1*)hPOWHEG->Clone(Form("ratPOWHEGXS_%d_%d", iR, GetNN()));
      hRatP->Divide(hPOWHEG, hData, 1., 1., "");
      hRatP->SetStats(0);
      hRatP->SetLineColor(kRed+1); hRatP->SetLineStyle(3);
      hRatP->SetLineWidth(3); hRatP->SetMarkerStyle(0);
      hRatP->Draw("l same");
    }

    // Legend in first panel
    if (iR == 0) {
      TLegend* leg = new TLegend(0.22, 0.45, 0.75, 0.88, NULL, "brNDC");
      leg->SetTextSize(0.04); leg->SetBorderSize(0); leg->SetFillStyle(0);

      TH1F* dumMC = new TH1F(Form("dumMCTXS_%d", GetNN()), "", 1, 0, 1);
      dumMC->SetLineColor(kGray+2); dumMC->SetLineStyle(7); dumMC->SetLineWidth(2);
      leg->AddEntry(dumMC, "MC truth", "l");

      for (int im = 0; im < (int)models.size(); ++im) {
        TH1F* dum = new TH1F(Form("dumMXS_%d_%d", im, GetNN()), "", 1, 0, 1);
        dum->SetLineColor(models[im].lineColor);
        dum->SetLineStyle(models[im].lineStyle);
        dum->SetLineWidth(models[im].lineWidth);
        leg->AddEntry(dum, models[im].legLabel.Data(), "l");
      }

      TH1F* dumP = new TH1F(Form("dumPXS_%d", GetNN()), "", 1, 0, 1);
      dumP->SetLineColor(kRed+1); dumP->SetLineStyle(3); dumP->SetLineWidth(3);
      leg->AddEntry(dumP, "POWHEG NLO", "l");

      leg->Draw();
    }
  }

  TString pdfNameB = Form("%s/ModelComparison_Xsec_AllR%s.pdf", outputDir.Data(), suffix.Data());
  if (GetDRAWPLOTS()) {
    canB->Print(pdfNameB.Data());
  }
}

// ============================================
// Jet pT Resolution
// ============================================

void DrawJetResolution() {
  if (!GetDrawJetResolution()) return;

  std::cout << "[DrawJetResolution] Starting..." << std::endl;

  std::vector<RConfig> rConfigs = GetRConfigs();
  std::vector<Color_t> rColors = GetRColorsForConfigs(rConfigs);
  TString outputDir = GetOutputDir();

  // pT classes for resolution distributions
  std::vector<std::pair<int, int>> pTClasses = {{10, 15}, {15, 20}, {20, 40}, {40, 60}};
  Color_t ptColors[] = {kBlack, kRed+1, kBlue+1, kGreen+2};

  // Vectors for resolution vs pT (all R overlaid)
  std::vector<TGraphErrors*> grMeans, grSigmas;
  std::vector<TString> grLabels;

  for (size_t i = 0; i < rConfigs.size(); ++i) {
    const auto& config = rConfigs[i];

    TFile* mcFile = TFile::Open(config.GetMcFile().Data(), "READ");
    if (!mcFile || mcFile->IsZombie()) {
      std::cerr << "[DrawJetResolution] Cannot open MC file for R=" << config.R << std::endl;
      grMeans.push_back(nullptr);
      grSigmas.push_back(nullptr);
      grLabels.push_back("");
      continue;
    }

    TH2* h2Res = (TH2*)mcFile->Get(Form("%s/%s", config.mcDir.Data(), GetJetResolutionObj()));
    if (!h2Res) {
      std::cerr << "[DrawJetResolution] Resolution histogram not found for R=" << config.R
                << " (" << GetJetResolutionObj() << ")" << std::endl;
      mcFile->Close();
      delete mcFile;
      grMeans.push_back(nullptr);
      grSigmas.push_back(nullptr);
      grLabels.push_back("");
      continue;
    }

    TH2* h2ResClone = (TH2*)h2Res->Clone(Form("h2JetRes_R%.1f_%d", config.R, (int)i));
    h2ResClone->SetDirectory(0);
    mcFile->Close();
    delete mcFile;

    // --- Part 1: Resolution distribution per pT class ---
    TCanvas* can = new TCanvas(Form("JetResolution_R%.1f_%d", config.R, ++GetNN()),
                               "", 700, 500);
    setpad(can, 0.03, 0.15, 0.15, 0.05);
    can->SetLogy(1);
    can->Draw();

    TLegend* leg = new TLegend(0.52, 0.52, 0.93, 0.93, NULL, "brNDC");
    leg->SetTextSize(0.035);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry((TObject*)0, Form("R = %.1f", config.R), "");

    bool first = true;
    for (size_t j = 0; j < pTClasses.size(); ++j) {
      int ptLow = pTClasses[j].first;
      int ptHigh = pTClasses[j].second;

      int binLow = h2ResClone->GetXaxis()->FindBin(ptLow);
      int binHigh = h2ResClone->GetXaxis()->FindBin(ptHigh - 1e-6);

      TH1D* hProj = h2ResClone->ProjectionY(
          Form("JetRes_R%.1f_pt%d_%d_%d", config.R, ptLow, ptHigh, (int)i), binLow, binHigh);
      hProj->SetDirectory(0);

      // Normalize to unit area (probability density)
      if (hProj->Integral() > 0) hProj->Scale(1.0 / hProj->Integral(), "width");

      hProj->SetLineColor(ptColors[j]);
      hProj->SetMarkerColor(ptColors[j]);
      hProj->SetMarkerStyle(20);
      hProj->SetMarkerSize(0.5);
      hProj->SetLineWidth(2);

      if (first) {
        hset(*hProj, "(#it{p}_{T}^{true} #minus #it{p}_{T}^{reco}) / #it{p}_{T}^{true}",
             "Probability density", 1.1, 1.3, 0.05, 0.05, 0.01, 0.01, 0.045, 0.045);
        hProj->GetXaxis()->SetRangeUser(-0.5, 1.0);
        hProj->GetYaxis()->SetRangeUser(1e-3, 50);
        hProj->Draw("hist e");
        first = false;
      } else {
        hProj->Draw("hist e same");
      }

      leg->AddEntry(hProj, Form("#it{p}_{T}^{true} [%d, %d] GeV/#it{c}, #mu=%.3f, #sigma=%.3f",
                                ptLow, ptHigh, hProj->GetMean(), hProj->GetRMS()), "l");
    }

    leg->Draw();
    AddUESubtractedLabel(can);

    if (GetDRAWPLOTS()) {
      can->Print(Form("%s/JetResolution_R%.1f.pdf", outputDir.Data(), config.R));
    }

    // --- Part 2: Compute mean and sigma vs pT for this R ---
    std::vector<double> ptCenters, ptErrors, means, meanErrors, sigmas, sigmaErrors;

    for (int ib = 1; ib <= h2ResClone->GetNbinsX(); ++ib) {
      double ptCenter = h2ResClone->GetXaxis()->GetBinCenter(ib);
      double ptHalfW = h2ResClone->GetXaxis()->GetBinWidth(ib) / 2.0;

      if (ptCenter < 5 || ptCenter > 150) continue;

      TH1D* hSlice = h2ResClone->ProjectionY(Form("slice_R%.1f_b%d_%d", config.R, ib, (int)i), ib, ib);
      if (hSlice->GetEntries() < 30) { delete hSlice; continue; }

      ptCenters.push_back(ptCenter);
      ptErrors.push_back(ptHalfW);
      means.push_back(hSlice->GetMean());
      meanErrors.push_back(hSlice->GetMeanError());
      sigmas.push_back(hSlice->GetRMS());
      sigmaErrors.push_back(hSlice->GetRMSError());
      delete hSlice;
    }

    if (ptCenters.empty()) {
      grMeans.push_back(nullptr);
      grSigmas.push_back(nullptr);
      grLabels.push_back("");
      delete h2ResClone;
      continue;
    }

    int nPts = ptCenters.size();
    TGraphErrors* grMean = new TGraphErrors(nPts, ptCenters.data(), means.data(),
                                            ptErrors.data(), meanErrors.data());
    TGraphErrors* grSigma = new TGraphErrors(nPts, ptCenters.data(), sigmas.data(),
                                             ptErrors.data(), sigmaErrors.data());

    grMean->SetLineColor(rColors[i]);
    grMean->SetMarkerColor(rColors[i]);
    grMean->SetMarkerStyle(20 + i);
    grMean->SetMarkerSize(0.8);

    grSigma->SetLineColor(rColors[i]);
    grSigma->SetMarkerColor(rColors[i]);
    grSigma->SetMarkerStyle(20 + i);
    grSigma->SetMarkerSize(0.8);

    grMeans.push_back(grMean);
    grSigmas.push_back(grSigma);
    grLabels.push_back(Form("R = %.1f", config.R));

    delete h2ResClone;
  }

  // --- Combined: Resolution mean & sigma vs jet pT ---
  // Top pad: mean (JES), Bottom pad: sigma (JER)
  Filipad2* padVsPt = new Filipad2("JetResVsPt", ++GetNN(), 2, 0.4, 100, 50, 0.7, 1, 1);
  padVsPt->Draw();

  TPad* pMean = padVsPt->GetPad(1);
  TPad* pSigma = padVsPt->GetPad(2);

  // --- Compute auto y-ranges from data ---
  double meanYMin = 1e9, meanYMax = -1e9;
  double sigmaYMin = 1e9, sigmaYMax = -1e9;
  for (size_t ig = 0; ig < grMeans.size(); ++ig) {
    if (grMeans[ig]) {
      for (int ip = 0; ip < grMeans[ig]->GetN(); ++ip) {
        double y = grMeans[ig]->GetY()[ip];
        double ey = grMeans[ig]->GetEY()[ip];
        if (y - ey < meanYMin) meanYMin = y - ey;
        if (y + ey > meanYMax) meanYMax = y + ey;
      }
    }
    if (grSigmas[ig]) {
      for (int ip = 0; ip < grSigmas[ig]->GetN(); ++ip) {
        double y = grSigmas[ig]->GetY()[ip];
        double ey = grSigmas[ig]->GetEY()[ip];
        if (y - ey < sigmaYMin) sigmaYMin = y - ey;
        if (y + ey > sigmaYMax) sigmaYMax = y + ey;
      }
    }
  }
  // Add 30% padding
  double meanRange = meanYMax - meanYMin;
  if (meanRange < 0.01) meanRange = 0.1; // fallback
  meanYMin -= 0.3 * meanRange;
  meanYMax += 0.3 * meanRange;
  double sigmaRange = sigmaYMax - sigmaYMin;
  if (sigmaRange < 0.01) sigmaRange = 0.1;
  sigmaYMin = 0; // sigma is always positive, start from 0
  sigmaYMax += 0.3 * sigmaRange;

  // --- Mean (JES) pad ---
  pMean->cd();
  pMean->SetGridy(1);
  TH1* hFrameMean = pMean->DrawFrame(0, meanYMin, 160, meanYMax);
  padVsPt->Hset(hFrameMean, "",
    "#LT(#it{p}_{T}^{true} #minus #it{p}_{T}^{reco}) / #it{p}_{T}^{true}#GT");

  TLegend* legVsPt = new TLegend(0.6, 0.15, 0.93, 0.55, NULL, "brNDC");
  legVsPt->SetTextSize(0.06);
  legVsPt->SetBorderSize(0);
  legVsPt->SetFillStyle(0);
  legVsPt->AddEntry((TObject*)0, "Jet Energy Scale", "");

  for (size_t i = 0; i < grMeans.size(); ++i) {
    if (!grMeans[i]) continue;
    pMean->cd();
    grMeans[i]->Draw("P same");
    legVsPt->AddEntry(grMeans[i], grLabels[i].Data(), "lp");
  }

  pMean->cd();
  TLine* zeroLine = new TLine(0, 0, 160, 0);
  zeroLine->SetLineStyle(2);
  zeroLine->SetLineColor(kGray+1);
  zeroLine->Draw();
  legVsPt->Draw();
  AddUESubtractedLabel(pMean);

  // --- Sigma (JER) pad ---
  pSigma->cd();
  pSigma->SetGridy(1);
  TH1* hFrameSigma = pSigma->DrawFrame(0, sigmaYMin, 160, sigmaYMax);
  padVsPt->Hset(hFrameSigma,
    "#it{p}_{T, jet}^{true} (GeV/#it{c})",
    "#sigma(#Deltap_{T} / p_{T}^{true})");

  for (size_t i = 0; i < grSigmas.size(); ++i) {
    if (!grSigmas[i]) continue;
    pSigma->cd();
    grSigmas[i]->Draw("P same");
  }

  if (GetDRAWPLOTS()) {
    padVsPt->C->Print(Form("%s/JetResolution_VsPt.pdf", outputDir.Data()));
  }

  std::cout << "[DrawJetResolution] Completed!" << std::endl;
}

// ============================================
// Main function
// ============================================

// ============================================
// Common draw sequence for one config set
// ============================================

void RunAnalysisForCurrentSet() {
  ClearUnfoldCache();  // Reset between config sets to free memory

  std::vector<RConfig> rConfigs = GetRConfigs();
  TString outputDir = GetOutputDir();
  gSystem->MakeDirectory(outputDir.Data());

  // Open log file (tee stderr → terminal + outputDir/analysis.log)
  OpenLogFile(outputDir.Data());

  std::cerr << "========================================" << std::endl;
  std::cerr << "Running: " << GetCurrentConfigSet().name << std::endl;
  std::cerr << "Mode: " << (GetCurrentMode() == kUE ? "UE-sub" : "NonUE") << std::endl;
  std::cerr << "Unfolding: " << (GetUnfoldMethod() == kUnfoldBayes ? "Bayesian" : "SVD") << std::endl;
  std::cerr << "Jet eta: " << (GetCurrentConfigSet().jetEtaCut >= 0
    ? Form("fixed |eta|<%.1f (deta=%.1f)", GetCurrentConfigSet().jetEtaCut, 2.0*GetCurrentConfigSet().jetEtaCut)
    : "R-dependent |eta|<0.9-R") << std::endl;
  std::cerr << "R values: " << rConfigs.size() << std::endl;
  std::cerr << "Output: " << outputDir << std::endl;
  std::cerr << "========================================" << std::endl;

  if (GetDrawTrack()) DrawTrackObservables();
  if (GetDrawJet()) DrawJetObservables();
  if (GetDrawJetPtPart()) DrawJetPtPart();
  if (GetDrawJetArea()) DrawJetAreaPerR();
  if (GetDrawResponseMatrix()) DrawResponseMatrixPerR();
  if (GetDrawJetResolution()) DrawJetResolution();
  if (GetDrawPurityEfficiency()) DrawPurityEfficiency();
  DrawRCTEfficiency();
  if (GetDrawKinematicEfficiency()) DrawKinematicEfficiency();
  if (GetDrawInvariantYield()) DrawInvariantYield();
  if (GetDrawCrossSection()) DrawCrossSection();
  if (GetDrawCrossSection()) DrawModelComparison();
  if (GetDrawCrossSection()) SaveDataResults();
  if (GetDrawCustomRatios()) DrawCustomRatios();
  if (GetDrawConstituentPt()) DrawConstituentPtObservables();
  if (GetDrawNtracks()) DrawNtracksObservables();

  std::cerr << "========================================" << std::endl;
  std::cerr << "Set complete: " << GetCurrentConfigSet().name << std::endl;
  std::cerr << "========================================" << std::endl;

  CloseLogFile();
}

// ============================================
// Entry points
// ============================================

void DrawJetsMCRDependentNonUE() {
  for (const auto& set : GetAllRConfigSets()) {
    if (set.mode == kNonUE) {
      SetCurrentConfigSet(set);
      RunAnalysisForCurrentSet();
    }
  }
}

void DrawJetsMCRDependentUE() {
  for (const auto& set : GetAllRConfigSets()) {
    if (set.mode == kUE) {
      SetCurrentConfigSet(set);
      RunAnalysisForCurrentSet();
    }
  }
}

// Default entry point: run all uncommented config sets
// Usage: root -l DrawJetsMCRDependent.C
void DrawJetsMCRDependent() {
  for (const auto& set : GetAllRConfigSets()) {
    SetCurrentConfigSet(set);
    RunAnalysisForCurrentSet();
  }
}

// NOTE:
// Run with:
//   root -l DrawJetsMCRDependent.C          (runs all uncommented sets)
// Or inside a ROOT session:
//   .L DrawJetsMCRDependent.C+
//   DrawJetsMCRDependentNonUE();  // non-UE sets only
//   DrawJetsMCRDependentUE();     // UE-sub sets only
//   DrawJetsMCRDependent();       // all sets
