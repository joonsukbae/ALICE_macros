// ============================================================
// XsectionModelComparison.C
// Standalone macro for comparing jet cross-section measurements
// with theoretical predictions (PYTHIA8, POWHEG, Herwig7).
//
// Reads pre-computed data results from DataResults.root
// (produced by DrawJetsMCRDependent.C) — no unfolding needed.
//
// Usage:
//   root -l XsectionModelComparison.C
//   root -l 'XsectionModelComparison.C("path/to/DataResults.root")'
// ============================================================

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TParameter.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <map>

// ============================================================
// Configuration — edit paths here
// ============================================================

// Standalone model directory (KIAF productions)
const char* kStandaloneDir = "/Users/js/cernbox/workspace/O2Physics/jets/PYTHIA_standalone_KIAF";

// POWHEG file
const char* kPOWHEGFile = "/Users/js/cernbox/workspace/O2Physics/jets/POWHEG/Pythia8bJetSpectra_dijet_CT18nlo_MonashTune.root";

// Generator-level MC (jet-finder-mcp-charged output, THnSparse format)
// These are pure particle-level runs (no detector sim) with large statistics
struct GenLevelMC {
  TString name;       // display name
  TString file;       // ROOT file path
  TString dir;        // TDirectory inside file (e.g. "jet-finder-mcp-charged")
  TString evtCountDir;  // directory for event counter (e.g. "jet-hadron-recoil")
  TString evtCountHist; // histogram name for N_events (e.g. "hZvtxSelected")
  double R;           // jet radius
  double etaMax;      // jet |eta| < etaMax
  Color_t color;
  Int_t lineStyle;
  Int_t lineWidth;
  bool enabled;
};

std::vector<GenLevelMC> GetGenLevelMCs() {
  std::vector<GenLevelMC> v;
  v.push_back({
    "PYTHIA8 Monash (HY)",
    "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/601994_AnalysisResults.root",
    "jet-finder-mcp-charged", "jet-hadron-recoil", "hZvtxSelected",
    0.4, 0.5, kBlue, 1, 2, true
  });
  // v.push_back({
  //   "PYTHIA8 Rope (HY)",
  //   "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/606723_AnalysisResults.root",
  //   "jet-finder-mcp-charged", "jet-hadron-recoil_id46506", "hZvtxSelected",
  //   0.4, 0.5, kRed-4, 1, 2, true
  // });
  // v.push_back({
  //   "PYTHIA8 Shoving (HY)",
  //   "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/602260_AnalysisResults.root",
  //   "jet-finder-mcp-charged", "jet-hadron-recoil", "hZvtxSelected",
  //   0.4, 0.5, kGreen+3, 1, 2, true
  // });
  return v;
}

// ALICE MC truth files (one per R; set to "" to skip)
// These are the O2Physics AnalysisResults files used for the response matrix
// Format: {R, filePath, directory, isJJ}
struct MCTruthConfig {
  double R;
  TString file;
  TString dir;
  bool isJJ;
  // Cross-section efficiency: file/dir for h_mccollisions_eventselection + h2_jet_pt_part_eventselection
  // If xsecEffFile is empty, looks in the same file/dir
  TString xsecEffFile;
  TString xsecEffDir;
  int zvtxBin;  // which Y bin in h2 is the zvtx step (for ε_zvtx correction). 0 = auto-detect
};

// MC truth configurations (match your active RConfigSet)
std::vector<MCTruthConfig> GetMCTruthConfigs() {
  std::vector<MCTruthConfig> v;
  // LHC23k4i MB MC (2023 data) — xsec efficiency from 605134 (LHC23k4h, same generator)
  TString xsecEffFile = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/605134_AnalysisResults.root";
  TString xsecEffDir  = "jet-cross-section-efficiency_selMC";
  v.push_back({0.2, "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/594177_AnalysisResults.root", "jet-spectra-charged", false, xsecEffFile, xsecEffDir, 0});
  v.push_back({0.4, "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/626304_AnalysisResults.root", "jet-spectra-charged", false, xsecEffFile, xsecEffDir, 0});
  v.push_back({0.6, "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/626309_AnalysisResults.root", "jet-spectra-charged", false, xsecEffFile, xsecEffDir, 0});

  // // LHC24f3c MB MC (2022 data)
  // v.push_back({0.2, "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/594013_AnalysisResults.root", "jet-spectra-charged", false});
  // v.push_back({0.4, "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/594014_AnalysisResults.root", "jet-spectra-charged", false});
  // v.push_back({0.6, "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/598073_AnalysisResults.root", "jet-spectra-charged", false});
  return v;
}

// Output directory
const char* kOutputDir = "plots/ModelComparisons";

// Plot range
const double kPtMin = 5.0;
const double kPtMax = 140.0;
const double kRatioMin = 0.3;
const double kRatioMax = 2.5;

// σ_INEL for POWHEG yield conversion [mb]
const double kSigmaINEL = 79.0;

// ============================================================
// Model definition
// ============================================================

struct Model {
  TString name;       // internal key
  TString fileName;   // ROOT file name (in kStandaloneDir)
  TString legLabel;   // legend text
  Color_t color;
  Int_t lineStyle;
  Int_t lineWidth;
  bool enabled;
  bool isMB;          // true = MB production (correct), false = HardQCD (MPI bias)

  // Histogram name pattern for cross-section
  TString GetXsecName(int iR, bool isUE) const {
    TString ue = isUE ? "_UEsub" : "";
    return Form("hJetXsec%s_R%02d", ue.Data(), iR);
  }
  TString GetYieldName(int iR, bool isUE) const {
    TString ue = isUE ? "_UEsub" : "";
    return Form("hJetNormYield%s_R%02d", ue.Data(), iR);
  }
};

std::vector<Model> GetModels() {
  std::vector<Model> v;
  //                name           fileName                                              legend          color        ls lw  on   isMB
  v.push_back({"MonashMB",  "merged_MBMC_10B_PYTHIA8_Monash_pp_13600GeV.root",   "PYTHIA8 Monash (MB, KIAF)",  kCyan+1,     3, 2, true,  true});
  v.push_back({"Monash",    "merged_JJMC_PYTHIA8_Monash_pp_13600GeV.root",   "PYTHIA8 Monash (JJ, KIAF)",  kRed,        9, 2, true, false});
  // v.push_back({"Rope",      "merged_JJMC_PYTHIA8_Rope_pp_13600GeV.root",     "PYTHIA8 Rope (JJ, KIAF)",    kRed+1,      10, 2, true, false});
  // v.push_back({"Shoving",   "merged_JJMC_PYTHIA8_Shoving_pp_13600GeV.root",  "PYTHIA8 Shoving (JJ, KIAF)", kGreen+3,    10, 2, true, false});
  // v.push_back({"Herwig7",   "merged_JJMC_Herwig7_Herwig_pp_13600GeV.root",   "Herwig 7 (JJ, KIAF)",        kOrange+7,   10, 2, true, false});  // pT-hat stitching bug
  // Add MB variants for Rope/Shoving when available:
  // v.push_back({"RopeMB",   "merged_MBMC_PYTHIA8_Rope_pp_13600GeV.root",    "PYTHIA8 Rope (MB)",   kRed,        1, 2, true,  true});
  // v.push_back({"ShovingMB","merged_MBMC_PYTHIA8_Shoving_pp_13600GeV.root", "PYTHIA8 Shoving (MB)",kGreen+2,    1, 2, true,  true});
  return v;
}

// ============================================================
// Analysis truth binning (must match DrawJetsMCRDependent)
// ============================================================

const int kNTruthBins = 25;
const double kTruthEdges[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};

double DeltaEta(double R) { return 2.0 * (0.9 - R); }

// ============================================================
// Loading functions
// ============================================================

// Load data histogram from DataResults.root
TH1* LoadData(TFile* fData, const char* prefix, int iR) {
  TString name = Form("%s_R%02d", prefix, iR);
  TH1* h = (TH1*)fData->Get(name);
  if (h) {
    h = (TH1*)h->Clone(Form("%s_clone", name.Data()));
    h->SetDirectory(0);
  }
  return h;
}

// Load standalone model cross-section d²σ/(dpT dη)
TH1* LoadModelXsec(const Model& model, double R, bool isUE) {
  int iR = (int)(R * 10 + 0.5);
  TString path = Form("%s/%s", kStandaloneDir, model.fileName.Data());
  TFile* f = TFile::Open(path, "READ");
  if (!f || f->IsZombie()) { delete f; return nullptr; }

  TH1* hRaw = (TH1*)f->Get(model.GetXsecName(iR, isUE));
  if (!hRaw) { f->Close(); delete f; return nullptr; }

  TH1* hClone = (TH1*)hRaw->Clone(Form("model_%s_xsec_R%02d", model.name.Data(), iR));
  hClone->SetDirectory(0);
  f->Close(); delete f;

  // Rebin to analysis binning
  TH1* hRebin = hClone->Rebin(kNTruthBins, Form("model_%s_xsec_R%02d_rb", model.name.Data(), iR), kTruthEdges);
  hRebin->SetDirectory(0);
  delete hClone;

  // Content is dσ per bin [mb] → d²σ/(dpT dη)
  hRebin->Scale(1.0, "width");
  double dEta = DeltaEta(R);
  if (dEta > 0) hRebin->Scale(1.0 / dEta);

  return hRebin;
}

// Load standalone model invariant yield (1/Nevt) d²N/(dpT dη)
TH1* LoadModelYield(const Model& model, double R, bool isUE) {
  int iR = (int)(R * 10 + 0.5);
  TString path = Form("%s/%s", kStandaloneDir, model.fileName.Data());
  TFile* f = TFile::Open(path, "READ");
  if (!f || f->IsZombie()) { delete f; return nullptr; }

  TH1* hRaw = (TH1*)f->Get(model.GetYieldName(iR, isUE));
  if (!hRaw) { f->Close(); delete f; return nullptr; }

  TH1* hClone = (TH1*)hRaw->Clone(Form("model_%s_yield_R%02d", model.name.Data(), iR));
  hClone->SetDirectory(0);
  f->Close(); delete f;

  TH1* hRebin = hClone->Rebin(kNTruthBins, Form("model_%s_yield_R%02d_rb", model.name.Data(), iR), kTruthEdges);
  hRebin->SetDirectory(0);
  delete hClone;

  // Content is (1/Nevt) dN per bin per η → d²N/(dpT dη)
  hRebin->Scale(1.0, "width");
  // Already per-η, no further deltaEta scaling

  return hRebin;
}

// Load POWHEG cross-section d²σ/(dpT dη)
TH1* LoadPOWHEG(double R, bool isUE) {
  int iR = (int)(R * 10 + 0.5);
  if (iR < 2 || iR > 6) return nullptr;

  TFile* f = TFile::Open(kPOWHEGFile, "READ");
  if (!f || f->IsZombie()) { delete f; return nullptr; }

  TString hName = isUE ? Form("BkgSubtractedJetXSection_R%d", iR)
                       : Form("InclusiveJetXSection_R%d", iR);
  TH1* hRaw = (TH1*)f->Get(hName);
  if (!hRaw) { f->Close(); delete f; return nullptr; }

  TH1* hNev = (TH1*)f->Get("hNEvent");
  double nEvents = hNev ? hNev->GetEntries() : 20000000.0;

  TH1* hClone = (TH1*)hRaw->Clone(Form("powheg_xsec_R%02d_raw", iR));
  hClone->SetDirectory(0);
  f->Close(); delete f;

  TH1* hRebin = hClone->Rebin(kNTruthBins, Form("powheg_xsec_R%02d", iR), kTruthEdges);
  hRebin->SetDirectory(0);
  delete hClone;

  hRebin->Scale(1.0 / nEvents, "width");
  double dEta = DeltaEta(R);
  if (dEta > 0) hRebin->Scale(1.0 / dEta);

  return hRebin;
}

// Load ALICE MC truth as per-event yield: (1/Nevt) d²N/(dpT dη)
// For MB MC: counts / N_selected / dpT / Δη  (bin 5.5 = all cuts passed)
// For JJ MC: weighted counts / nevtsMCP_weighted / dpT / Δη
TH1* LoadMCTruthYield(const MCTruthConfig& cfg, bool isUE) {
  TFile* f = TFile::Open(cfg.file, "READ");
  if (!f || f->IsZombie()) { delete f; return nullptr; }

  const char* truthObj = isUE ? "h_jet_pt_part_rhoareasubtracted" : "h_jet_pt_part";
  TH1* hRaw = (TH1*)f->Get(Form("%s/%s", cfg.dir.Data(), truthObj));
  if (!hRaw) { f->Close(); delete f; return nullptr; }

  TH1* hClone = (TH1*)hRaw->Clone(Form("mctruth_R%02d_raw", (int)(cfg.R*10+0.5)));
  hClone->SetDirectory(0);

  // Get N_selected (bin 5.5 = events passing all cuts, where h_jet_pt_part is filled)
  double nEvt = 1.0;
  const char* collHist = cfg.isJJ ? "h_mccollisions_weighted" : "h_mccollisions";
  TH1* hMC = (TH1*)f->Get(Form("%s/%s", cfg.dir.Data(), collHist));
  if (hMC) nEvt = hMC->GetBinContent(hMC->FindBin(5.5));
  f->Close(); delete f;

  TH1* hRebin = hClone->Rebin(kNTruthBins,
    Form("mctruth_yield_R%02d", (int)(cfg.R*10+0.5)), kTruthEdges);
  hRebin->SetDirectory(0);
  delete hClone;

  double dEta = DeltaEta(cfg.R);
  if (nEvt > 0 && dEta > 0) {
    hRebin->Scale(1.0 / (nEvt * dEta), "width");
  }

  std::cout << "[MCTruth] R=" << cfg.R << " Nevt(bin5.5)=" << Form("%.3e", nEvt)
            << " isJJ=" << cfg.isJJ << std::endl;

  return hRebin;
}

// Load ALICE MC truth as cross-section: d²σ/(dpT dη) [mb/GeV]
// Formula: σ = σ_INEL / N_INEL(bin1) × jets(INEL bin1) / dpT / dη × 1/ε_zvtx
// Uses cross-section efficiency histograms for proper INEL normalization.
// ε_zvtx corrects for jet finder z-vertex cut on the numerator.
// Falls back to h_jet_pt_part + h_mccollisions bin 0.5 if xsec eff not available.
TH1* LoadMCTruthXsec(const MCTruthConfig& cfg, bool isUE) {
  // --- Try cross-section efficiency method first ---
  // Only use if the xsecEff file has the same R (check selectedJetsRadius in config)
  // For now: xsecEff file is R=0.4 only (605134). Skip for other R values.
  bool xsecEffRmatch = (TMath::Abs(cfg.R - 0.4) < 0.01);  // TODO: make configurable
  if (!isUE && cfg.xsecEffFile.Length() > 0 && xsecEffRmatch) {
    TFile* fEff = TFile::Open(cfg.xsecEffFile, "READ");
    if (fEff && !fEff->IsZombie()) {
      TDirectory* dEff = (TDirectory*)fEff->Get(cfg.xsecEffDir);
      if (dEff) {
        TH2* h2 = (TH2*)dEff->Get("h2_jet_pt_part_eventselection");
        TH1* hNevt = (TH1*)dEff->Get("h_mccollisions_eventselection");
        if (h2 && hNevt) {
          // N_INEL = bin 1 of h_mccollisions_eventselection (true INEL, no cuts)
          double nINEL = hNevt->GetBinContent(1);

          // INEL jets (bin 1) — jet finder zvtx is baked in
          TH1D* hJets = (TH1D*)h2->ProjectionX(
            Form("mctruth_xsec_R%02d_raw", (int)(cfg.R*10+0.5)), 1, 1);
          hJets->SetDirectory(0);

          // Find ε_zvtx: step where "zvtx" label appears
          double eps_zvtx = 1.0;
          for (int iy = 2; iy <= h2->GetNbinsY(); iy++) {
            TString label = h2->GetYaxis()->GetBinLabel(iy);
            if (label.Contains("zvtx") || label.Contains("Zvtx")) {
              double nBefore = hNevt->GetBinContent(iy - 1);
              double nAfter = hNevt->GetBinContent(iy);
              if (nBefore > 0) eps_zvtx = nAfter / nBefore;
              break;
            }
          }

          fEff->Close(); delete fEff;

          TH1* hRebin = hJets->Rebin(kNTruthBins,
            Form("mctruth_xsec_R%02d", (int)(cfg.R*10+0.5)), kTruthEdges);
          hRebin->SetDirectory(0);
          delete hJets;

          double dEta = DeltaEta(cfg.R);
          if (nINEL > 0 && dEta > 0 && eps_zvtx > 0) {
            hRebin->Scale(kSigmaINEL / (nINEL * dEta * eps_zvtx), "width");
          }

          std::cout << "[MCTruthXsec] R=" << cfg.R
                    << " N_INEL=" << Form("%.3e", nINEL)
                    << " ε_zvtx=" << Form("%.4f", eps_zvtx)
                    << " (from " << cfg.xsecEffDir << ")" << std::endl;
          return hRebin;
        }
      }
      fEff->Close(); delete fEff;
    }
  }

  // --- Fallback: h_jet_pt_part + h_mccollisions bin 0.5 (biased by selection ε) ---
  TFile* f = TFile::Open(cfg.file, "READ");
  if (!f || f->IsZombie()) { delete f; return nullptr; }

  const char* truthObj = isUE ? "h_jet_pt_part_rhoareasubtracted" : "h_jet_pt_part";
  TH1* hRaw = (TH1*)f->Get(Form("%s/%s", cfg.dir.Data(), truthObj));
  if (!hRaw) { f->Close(); delete f; return nullptr; }
  TH1* hJets = (TH1*)hRaw->Clone(Form("mctruth_xsec_R%02d_raw", (int)(cfg.R*10+0.5)));
  hJets->SetDirectory(0);

  double nINEL = 1.0;
  const char* collHist = cfg.isJJ ? "h_mccollisions_weighted" : "h_mccollisions";
  TH1* hMC = (TH1*)f->Get(Form("%s/%s", cfg.dir.Data(), collHist));
  if (hMC) nINEL = hMC->GetBinContent(hMC->FindBin(0.5));
  f->Close(); delete f;

  TH1* hRebin = hJets->Rebin(kNTruthBins,
    Form("mctruth_xsec_R%02d", (int)(cfg.R*10+0.5)), kTruthEdges);
  hRebin->SetDirectory(0);
  delete hJets;

  double dEta = DeltaEta(cfg.R);
  if (cfg.isJJ) {
    hRebin->Scale(1.0, "width");
    if (dEta > 0) hRebin->Scale(1.0 / dEta);
  } else {
    if (nINEL > 0 && dEta > 0) {
      hRebin->Scale(kSigmaINEL / (nINEL * dEta), "width");
    }
  }

  std::cout << "[MCTruthXsec] FALLBACK R=" << cfg.R << " N(bin0.5)=" << Form("%.3e", nINEL)
            << " (biased, no ε_zvtx correction)" << std::endl;
  return hRebin;
}

// Load generator-level MC cross-section from THnSparse (jet-finder-mcp-charged output)
// d²σ/(dpT dη) = counts × σ_INEL / N_events / dpT / Δη
TH1* LoadGenLevelXsec(const GenLevelMC& cfg) {
  TFile* f = TFile::Open(cfg.file, "READ");
  if (!f || f->IsZombie()) { delete f; return nullptr; }

  // Get event count
  TDirectory* evtDir = (TDirectory*)f->Get(cfg.evtCountDir);
  if (!evtDir) { f->Close(); delete f; return nullptr; }
  TH1* hEvt = (TH1*)evtDir->Get(cfg.evtCountHist);
  double nEvents = hEvt ? hEvt->GetEntries() : 0;
  if (nEvents <= 0) { f->Close(); delete f; return nullptr; }

  // Get jet THnSparse
  TDirectory* jetDir = (TDirectory*)f->Get(cfg.dir);
  if (!jetDir) { f->Close(); delete f; return nullptr; }
  THnBase* hn = (THnBase*)jetDir->Get("hJetMCP");
  if (!hn) { f->Close(); delete f; return nullptr; }

  // Apply eta cut: |eta| < etaMax (axis 2)
  int etaLo = hn->GetAxis(2)->FindBin(-cfg.etaMax + 0.001);
  int etaHi = hn->GetAxis(2)->FindBin(cfg.etaMax - 0.001);
  hn->GetAxis(2)->SetRange(etaLo, etaHi);

  // Project pT (axis 1)
  TH1D* hPt = (TH1D*)hn->Projection(1);
  hPt->SetDirectory(0);
  f->Close(); delete f;

  // Rebin: use only bin edges compatible with the source binning (5 GeV grid)
  // Filter kTruthEdges to keep only multiples of the source bin width
  double srcBinW = hPt->GetBinWidth(1);
  std::vector<double> compatEdges;
  for (int i = 0; i <= kNTruthBins; ++i) {
    double edge = kTruthEdges[i];
    // Accept if edge is a multiple of srcBinW (within tolerance)
    if (TMath::Abs(TMath::Nint(edge / srcBinW) * srcBinW - edge) < 0.01) {
      compatEdges.push_back(edge);
    }
  }
  int nCompatBins = (int)compatEdges.size() - 1;
  TH1* hRebin = hPt->Rebin(nCompatBins,
    Form("genmc_%s_xsec_R%02d", cfg.name.Data(), (int)(cfg.R*10+0.5)), compatEdges.data());
  hRebin->SetDirectory(0);
  delete hPt;

  // Normalize to per-event yield: (1/Nevt) d²N/(dpT dη)
  double dEta = 2.0 * cfg.etaMax;
  hRebin->Scale(1.0 / (nEvents * dEta), "width");

  std::cout << "[GenLevelMC] " << cfg.name << " R=" << cfg.R
            << " Nevents=" << Form("%.3e", nEvents)
            << " dEta=" << dEta << std::endl;

  return hRebin;
}

// ============================================================
// Drawing helpers
// ============================================================

static int gNN = 0;

// Make ratio histogram. If binnings differ, interpolate hNum into hDen's binning.
TH1* MakeRatio(TH1* hNum, TH1* hDen, const char* name) {
  if (!hNum || !hDen) return nullptr;
  if (hNum->GetNbinsX() == hDen->GetNbinsX()) {
    // Same binning — direct divide
    TH1* hRat = (TH1*)hNum->Clone(name);
    hRat->SetDirectory(0);
    hRat->Divide(hNum, hDen, 1., 1., "");
    return hRat;
  }
  // Different binning — build ratio in hDen's binning by matching bin centers
  TH1* hRat = (TH1*)hDen->Clone(name);
  hRat->SetDirectory(0);
  hRat->Reset();
  for (int i = 1; i <= hDen->GetNbinsX(); ++i) {
    double pt = hDen->GetBinCenter(i);
    int jBin = hNum->FindBin(pt);
    if (jBin < 1 || jBin > hNum->GetNbinsX()) continue;
    double num = hNum->GetBinContent(jBin);
    double den = hDen->GetBinContent(i);
    if (den > 0 && num > 0) {
      hRat->SetBinContent(i, num / den);
      // Simple error propagation (relative errors added in quadrature)
      double eNum = hNum->GetBinError(jBin) / num;
      double eDen = hDen->GetBinError(i) / den;
      hRat->SetBinError(i, (num/den) * TMath::Sqrt(eNum*eNum + eDen*eDen));
    }
  }
  return hRat;
}

void StyleRatio(TH1* h, Color_t color, Int_t lstyle, Int_t lwidth) {
  h->SetStats(0);
  h->SetLineColor(color);
  h->SetLineStyle(lstyle);
  h->SetLineWidth(lwidth);
  h->SetMarkerStyle(0);
}

// ============================================================
// Main comparison function
// ============================================================

void XsectionModelComparison(
  const char* dataFile = "../plots/RDependentComparison/LHC23-pass4-Thin_small_LHC23k4i_2022customTuner/DataResults.root",
  const char* outputDir = "",
  bool isUE = false)
{
  gStyle->SetOptStat(0);

  // Default output dir = dataFile directory + /ModelComparisons/
  TString outDir = outputDir;
  if (outDir.IsNull()) {
    outDir = Form("%s/ModelComparisons", gSystem->DirName(dataFile));
  }
  gSystem->mkdir(outDir, true);
  outputDir = outDir.Data();

  TString suffix = isUE ? "_UE" : "";

  // ---- Load data ----
  TFile* fData = TFile::Open(dataFile, "READ");
  if (!fData || fData->IsZombie()) {
    std::cerr << "Cannot open data file: " << dataFile << std::endl;
    return;
  }

  TParameter<int>* pNR = (TParameter<int>*)fData->Get("nR");
  int nR = pNR ? pNR->GetVal() : 0;
  if (nR == 0) { std::cerr << "No R values in data file" << std::endl; return; }

  std::vector<double> rValues;
  for (int i = 0; i < nR; ++i) {
    TParameter<double>* pR = (TParameter<double>*)fData->Get(Form("R_%d", i));
    if (pR) rValues.push_back(pR->GetVal());
  }
  std::cout << "Loaded " << rValues.size() << " R values: ";
  for (double r : rValues) std::cout << r << " ";
  std::cout << std::endl;

  // Load data cross-sections and invariant yields
  std::vector<TH1*> dataXsec, dataYield;
  for (double R : rValues) {
    int iR = (int)(R * 10 + 0.5);
    dataXsec.push_back(LoadData(fData, "hDataXsec", iR));
    dataYield.push_back(LoadData(fData, "hDataInvYield", iR));
  }

  // ---- Load models ----
  std::vector<Model> models = GetModels();
  std::vector<MCTruthConfig> mcConfigs = GetMCTruthConfigs();

  // ---- Load POWHEG ----
  std::vector<TH1*> powhegXsec;
  for (double R : rValues) powhegXsec.push_back(LoadPOWHEG(R, isUE));

  // ---- Load MC truth (per-event yield, bin 5.5) ----
  std::vector<TH1*> mcTruthYield;
  for (double R : rValues) {
    TH1* h = nullptr;
    for (const auto& cfg : mcConfigs) {
      if (TMath::Abs(cfg.R - R) < 0.01) {
        h = LoadMCTruthYield(cfg, isUE);
        break;
      }
    }
    mcTruthYield.push_back(h);
  }

  // ---- Load MC truth (cross-section) ----
  // Only R values with xsecEff match get true INEL; others are nullptr
  std::vector<TH1*> mcTruthXsec;
  for (double R : rValues) {
    TH1* h = nullptr;
    bool isINEL = (TMath::Abs(R - 0.4) < 0.01);  // only R=0.4 has proper INEL from 605134
    if (isINEL) {
      for (const auto& cfg : mcConfigs) {
        if (TMath::Abs(cfg.R - R) < 0.01) {
          h = LoadMCTruthXsec(cfg, isUE);
          break;
        }
      }
    }
    mcTruthXsec.push_back(h);
  }

  // ---- Derive MC truth INEL yield = xsec / σ_INEL ----
  std::vector<TH1*> mcTruthYieldINEL;
  for (int iR = 0; iR < nR; ++iR) {
    if (mcTruthXsec[iR]) {
      TH1* h = (TH1*)mcTruthXsec[iR]->Clone(Form("mctruth_yield_INEL_R%02d", (int)(rValues[iR]*10+0.5)));
      h->SetDirectory(0);
      h->Scale(1.0 / kSigmaINEL);
      mcTruthYieldINEL.push_back(h);
    } else {
      mcTruthYieldINEL.push_back(nullptr);
    }
  }

  // ---- Load generator-level MC (per-event yield, THnSparse) ----
  std::vector<GenLevelMC> genMCs = GetGenLevelMCs();
  std::map<int, std::vector<std::pair<const GenLevelMC*, TH1*>>> genMCYield;
  for (int iR = 0; iR < nR; ++iR) {
    for (const auto& gmc : genMCs) {
      if (!gmc.enabled) continue;
      if (TMath::Abs(gmc.R - rValues[iR]) < 0.01) {
        TH1* h = LoadGenLevelXsec(gmc);  // returns per-event yield now
        if (h) genMCYield[iR].push_back({&gmc, h});
      }
    }
  }

  // ============================================================
  // Canvas 1: Per-event yield Model/Data ratio (per-R horizontal)
  // All models normalized as (1/Nevt) d²N/(dpT dη) for consistent comparison
  // ============================================================
  {
    int canW = 450 * nR;
    TCanvas* can = new TCanvas(Form("yieldRatioAll_%d", ++gNN), "Yield Model/Data", canW, 500);
    can->Divide(nR, 1, 0.001, 0.001);

    for (int iR = 0; iR < nR; ++iR) {
      can->cd(iR + 1);
      gPad->SetGrid(0);
      gPad->SetLeftMargin(iR == 0 ? 0.18 : 0.05);
      gPad->SetRightMargin(iR == nR - 1 ? 0.05 : 0.01);
      gPad->SetBottomMargin(0.15);
      gPad->SetTopMargin(0.05);

      TH1F* frame = new TH1F(Form("frXR_%d_%d", iR, gNN), "", 100, kPtMin, kPtMax);
      frame->SetMinimum(kRatioMin);
      frame->SetMaximum(kRatioMax);
      frame->GetXaxis()->SetTitle("#it{p}_{T}^{true} (GeV/#it{c})");
      frame->GetYaxis()->SetTitle(iR == 0 ? "Model / Data" : "");
      frame->GetXaxis()->SetTitleSize(0.06);
      frame->GetYaxis()->SetTitleSize(0.06);
      frame->GetXaxis()->SetLabelSize(0.05);
      frame->GetYaxis()->SetLabelSize(iR == 0 ? 0.05 : 0.0);
      frame->GetXaxis()->SetTitleOffset(1.1);
      frame->GetYaxis()->SetTitleOffset(1.3);
      frame->Draw("AXIS");

      TLatex* lat = new TLatex(0.5, 0.92, Form("#it{R} = %.1f", rValues[iR]));
      lat->SetNDC(); lat->SetTextSize(0.06); lat->SetTextAlign(22); lat->Draw();

      TLine* line = new TLine(kPtMin, 1.0, kPtMax, 1.0);
      line->SetLineStyle(2); line->Draw();

      if (!dataYield[iR]) continue;

      // MC truth ratio — selected (per-selected-event yield, bin 5.5)
      if (mcTruthYield[iR]) {
        TH1* hR = MakeRatio(mcTruthYield[iR], dataYield[iR], Form("ratMCTsel_%d_%d", iR, gNN));
        StyleRatio(hR, kGray+1, 2, 2);
        hR->Draw("hist e same");
      }

      // MC truth ratio — INEL (corrected, xsec/σ_INEL)
      if (mcTruthYieldINEL[iR]) {
        TH1* hR = MakeRatio(mcTruthYieldINEL[iR], dataYield[iR], Form("ratMCTinel_%d_%d", iR, gNN));
        StyleRatio(hR, kGray+1, 1, 2);  // same color, solid line
        hR->Draw("hist e same");
      }

      // Standalone model ratios (per-event yield)
      for (const auto& model : models) {
        if (!model.enabled) continue;
        TH1* hModel = LoadModelYield(model, rValues[iR], isUE);
        if (!hModel) continue;
        TH1* hR = MakeRatio(hModel, dataYield[iR], Form("ratM_%s_%d_%d", model.name.Data(), iR, gNN));
        StyleRatio(hR, model.color, model.lineStyle, model.lineWidth);
        hR->Draw("l same");
        delete hModel;
      }

      // POWHEG ratio (xsec / σ_INEL → per-event yield)
      if (powhegXsec[iR]) {
        TH1* hPY = (TH1*)powhegXsec[iR]->Clone(Form("pwYieldC1_%d_%d", iR, gNN));
        hPY->Scale(1.0 / kSigmaINEL);
        TH1* hR = MakeRatio(hPY, dataYield[iR], Form("ratPW_%d_%d", iR, gNN));
        StyleRatio(hR, kMagenta+1, 7, 2);
        hR->Draw("l same");
      }

      // Generator-level MC ratios (already per-event yield)
      if (genMCYield.count(iR)) {
        for (const auto& p : genMCYield[iR]) {
          TH1* hR = MakeRatio(p.second, dataYield[iR], Form("ratGMC_%s_%d_%d", p.first->name.Data(), iR, gNN));
          StyleRatio(hR, p.first->color, p.first->lineStyle, p.first->lineWidth);
          hR->Draw("l same");
        }
      }

      // Data/Data = 1 with error bars
      TH1* hDataSelf = (TH1*)dataYield[iR]->Clone(Form("dataSelf_%d_%d", iR, gNN));
      hDataSelf->Divide(dataXsec[iR], dataXsec[iR], 1., 1., "");
      hDataSelf->SetLineColor(kBlack);
      hDataSelf->SetMarkerStyle(20);
      hDataSelf->SetMarkerSize(0.5);
      hDataSelf->Draw("pe same");

      // Legend (first panel only) — built from ALL sources, not just this R
      if (iR == 0) {
        TLegend* leg = new TLegend(0.20976, 0.612003, 0.792212, 0.88709);
        leg->SetTextSize(0.033); leg->SetBorderSize(0); leg->SetFillStyle(0);

        leg->AddEntry(hDataSelf, "ALICE data", "lpe");

        // MC truth selected (dashed)
        bool hasMCTsel = false;
        for (auto h : mcTruthYield) if (h) { hasMCTsel = true; break; }
        if (hasMCTsel) {
          TH1F* d = new TH1F(Form("dMCTsel_%d", gNN), "", 1, 0, 1);
          d->SetLineColor(kGray+1); d->SetLineStyle(2); d->SetLineWidth(2);
          leg->AddEntry(d, "MC truth sel (PYTHIA8)", "l");
        }

        // MC truth INEL (solid)
        bool hasMCTinel = false;
        for (auto h : mcTruthYieldINEL) if (h) { hasMCTinel = true; break; }
        if (hasMCTinel) {
          TH1F* d = new TH1F(Form("dMCTinel_%d", gNN), "", 1, 0, 1);
          d->SetLineColor(kGray+1); d->SetLineStyle(1); d->SetLineWidth(2);
          leg->AddEntry(d, "MC truth INEL (PYTHIA8)", "l");
        }

        // Generator-level MC: check all R
        for (const auto& gmc : genMCs) {
          if (!gmc.enabled) continue;
          TH1F* d = new TH1F(Form("dGMC_%s_%d", gmc.name.Data(), gNN), "", 1, 0, 1);
          d->SetLineColor(gmc.color); d->SetLineStyle(gmc.lineStyle); d->SetLineWidth(gmc.lineWidth);
          leg->AddEntry(d, gmc.name.Data(), "l");
        }

        for (const auto& model : models) {
          if (!model.enabled) continue;
          TH1F* d = new TH1F(Form("d_%s_%d", model.name.Data(), gNN), "", 1, 0, 1);
          d->SetLineColor(model.color); d->SetLineStyle(model.lineStyle); d->SetLineWidth(model.lineWidth);
          leg->AddEntry(d, model.legLabel.Data(), "l");
        }

        // POWHEG: check any R
        bool hasPW = false;
        for (auto h : powhegXsec) if (h) { hasPW = true; break; }
        if (hasPW) {
          TH1F* d = new TH1F(Form("dPW_%d", gNN), "", 1, 0, 1);
          d->SetLineColor(kMagenta+1); d->SetLineStyle(7); d->SetLineWidth(2);
          leg->AddEntry(d, "POWHEG NLO (Hadi)", "l");
        }

        leg->Draw();
      }
    }

    TString pdf = Form("%s/ModelComparison_Yield_AllR%s.pdf", outputDir, suffix.Data());
    can->Print(pdf);
    std::cout << "Saved: " << pdf << std::endl;
  }

  // ============================================================
  // Canvas 2: Cross-section Model/Data ratio (per-R horizontal)
  // Models normalized as d²σ/(dpT dη): yield × σ_INEL
  // ============================================================
  {
    int canW = 450 * nR;
    TCanvas* can = new TCanvas(Form("xsecRatio_%d", ++gNN), "Xsec Model/Data", canW, 500);
    can->Divide(nR, 1, 0.001, 0.001);

    for (int iR = 0; iR < nR; ++iR) {
      can->cd(iR + 1);
      gPad->SetGrid(0);
      gPad->SetLeftMargin(iR == 0 ? 0.18 : 0.05);
      gPad->SetRightMargin(iR == nR - 1 ? 0.05 : 0.01);
      gPad->SetBottomMargin(0.15);
      gPad->SetTopMargin(0.05);

      TH1F* frame = new TH1F(Form("frXS_%d_%d", iR, gNN), "", 100, kPtMin, kPtMax);
      frame->SetMinimum(kRatioMin);
      frame->SetMaximum(kRatioMax);
      frame->GetXaxis()->SetTitle("#it{p}_{T}^{true} (GeV/#it{c})");
      frame->GetYaxis()->SetTitle(iR == 0 ? "Model / Data" : "");
      frame->GetXaxis()->SetTitleSize(0.06);
      frame->GetYaxis()->SetTitleSize(0.06);
      frame->GetXaxis()->SetLabelSize(0.05);
      frame->GetYaxis()->SetLabelSize(iR == 0 ? 0.05 : 0.0);
      frame->GetXaxis()->SetTitleOffset(1.1);
      frame->GetYaxis()->SetTitleOffset(1.3);
      frame->Draw("AXIS");

      TLatex* lat = new TLatex(0.5, 0.92, Form("#it{R} = %.1f", rValues[iR]));
      lat->SetNDC(); lat->SetTextSize(0.06); lat->SetTextAlign(22); lat->Draw();

      TLine* line = new TLine(kPtMin, 1.0, kPtMax, 1.0);
      line->SetLineStyle(2); line->Draw();

      if (!dataXsec[iR]) continue;

      // MC truth INEL ratio (cross-section, corrected)
      if (mcTruthXsec[iR]) {
        TH1* hR = MakeRatio(mcTruthXsec[iR], dataXsec[iR], Form("ratXsMCT_%d_%d", iR, gNN));
        StyleRatio(hR, kGray+1, 1, 2);
        hR->Draw("l same");
      }

      // Standalone model ratios (cross-section directly)
      for (const auto& model : models) {
        if (!model.enabled) continue;
        TH1* hModel = LoadModelXsec(model, rValues[iR], isUE);
        if (!hModel) continue;
        TH1* hR = MakeRatio(hModel, dataXsec[iR], Form("ratXsM_%s_%d_%d", model.name.Data(), iR, gNN));
        StyleRatio(hR, model.color, model.lineStyle, model.lineWidth);
        hR->Draw("l same");
        delete hModel;
      }

      // POWHEG ratio (already xsec)
      if (powhegXsec[iR]) {
        TH1* hR = MakeRatio(powhegXsec[iR], dataXsec[iR], Form("ratXsPW_%d_%d", iR, gNN));
        StyleRatio(hR, kMagenta+1, 7, 2);
        hR->Draw("l same");
      }

      // Generator-level MC ratios (yield × σ_INEL → xsec)
      if (genMCYield.count(iR)) {
        for (const auto& p : genMCYield[iR]) {
          TH1* hXs = (TH1*)p.second->Clone(Form("gmcXs_%s_%d_%d", p.first->name.Data(), iR, gNN));
          hXs->Scale(kSigmaINEL);
          TH1* hR = MakeRatio(hXs, dataXsec[iR], Form("ratXsGMC_%s_%d_%d", p.first->name.Data(), iR, gNN));
          StyleRatio(hR, p.first->color, p.first->lineStyle, p.first->lineWidth);
          hR->Draw("l same");
        }
      }

      // Data with errors
      TH1* hDataSelf = (TH1*)dataXsec[iR]->Clone(Form("dataSelfXs_%d_%d", iR, gNN));
      hDataSelf->Divide(dataXsec[iR], dataXsec[iR], 1., 1., "");
      hDataSelf->SetLineColor(kBlack); hDataSelf->SetMarkerStyle(20); hDataSelf->SetMarkerSize(0.5);
      hDataSelf->Draw("pe same");

      // Legend (first panel)
      if (iR == 0) {
        TLegend* leg = new TLegend(0.20976, 0.612003, 0.792212, 0.88709);
        leg->SetTextSize(0.033); leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->AddEntry(hDataSelf, "ALICE data", "lpe");

        bool hasMCT = false;
        for (auto h : mcTruthYield) if (h) { hasMCT = true; break; }
        if (hasMCT) {
          TH1F* d = new TH1F(Form("dXsMCT_%d", gNN), "", 1, 0, 1);
          d->SetLineColor(kGray+1); d->SetLineStyle(1); d->SetLineWidth(2);
          leg->AddEntry(d, "MC truth INEL (PYTHIA8)", "l");
        }
        for (const auto& gmc : genMCs) {
          if (!gmc.enabled) continue;
          TH1F* d = new TH1F(Form("dXsGMC_%s_%d", gmc.name.Data(), gNN), "", 1, 0, 1);
          d->SetLineColor(gmc.color); d->SetLineStyle(gmc.lineStyle); d->SetLineWidth(gmc.lineWidth);
          leg->AddEntry(d, gmc.name.Data(), "l");
        }
        for (const auto& model : models) {
          if (!model.enabled) continue;
          TH1F* d = new TH1F(Form("dXs_%s_%d", model.name.Data(), gNN), "", 1, 0, 1);
          d->SetLineColor(model.color); d->SetLineStyle(model.lineStyle); d->SetLineWidth(model.lineWidth);
          leg->AddEntry(d, model.legLabel.Data(), "l");
        }
        bool hasPW = false;
        for (auto h : powhegXsec) if (h) { hasPW = true; break; }
        if (hasPW) {
          TH1F* d = new TH1F(Form("dXsPW_%d", gNN), "", 1, 0, 1);
          d->SetLineColor(kMagenta+1); d->SetLineStyle(7); d->SetLineWidth(2);
          leg->AddEntry(d, "POWHEG NLO (Hadi)", "l");
        }
        leg->Draw();
      }
    }

    TString pdf = Form("%s/ModelComparison_Xsec_AllR%s.pdf", outputDir, suffix.Data());
    can->Print(pdf);
    std::cout << "Saved: " << pdf << std::endl;
  }

  // ============================================================
  // Canvas 3: Cross-section R-ratio comparison (R=0.2 / R=X)
  // ============================================================
  if (nR >= 2) {
    // Find R=0.2 index and collect denominator R values
    int iR02 = -1;
    for (int i = 0; i < nR; ++i) {
      if (TMath::Abs(rValues[i] - 0.2) < 0.01) { iR02 = i; break; }
    }
    if (iR02 < 0) { std::cout << "[RRatio] R=0.2 not found, skipping" << std::endl; }
    else {
    // Collect denominator R indices (all R > 0.2)
    std::vector<int> denIndices;
    for (int i = 0; i < nR; ++i) {
      if (rValues[i] > 0.25) denIndices.push_back(i);
    }
    int nPairs = (int)denIndices.size();
    int canW = 400 * nPairs;
    TCanvas* can = new TCanvas(Form("rratio_%d", ++gNN), "R-ratio Comparison", canW, 700);

    const float splitY = 0.38;

    for (int ip = 0; ip < nPairs; ++ip) {
      int j = denIndices[ip];
      double Rnum = 0.2, Rden = rValues[j];
      int iRnum = (int)(Rnum * 10 + 0.5);
      int iRden = (int)(Rden * 10 + 0.5);

      float padL = (float)ip / nPairs;
      float padR = (float)(ip + 1) / nPairs;
      float lm = (ip == 0) ? 0.18 : 0.04;
      float rm = (ip == nPairs - 1) ? 0.04 : 0.01;

      // Upper pad: R-ratio overlay
      TPad* pUp = new TPad(Form("pUpRR_%d_%d", ip, gNN), "", padL, splitY, padR, 1.0);
      pUp->SetTopMargin(0.06); pUp->SetBottomMargin(0.0);
      pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
      pUp->SetTicks(1, 1); pUp->Draw(); pUp->cd();

      TH1F* frUp = new TH1F(Form("frRRup_%d_%d", ip, gNN), "", 100, kPtMin, kPtMax);
      frUp->SetMinimum(0.0); frUp->SetMaximum(1.0);
      frUp->GetYaxis()->SetTitle(ip == 0 ? "#sigma(R=0.2) / #sigma(R=X)" : "");
      frUp->GetYaxis()->SetTitleSize(0.06); frUp->GetYaxis()->SetLabelSize(0.05);
      frUp->GetXaxis()->SetLabelSize(0); frUp->GetXaxis()->SetTitleSize(0);
      frUp->GetYaxis()->SetTitleOffset(1.2);
      frUp->Draw("AXIS");

      TLatex* lat = new TLatex(0.5, 0.92, Form("R=0.2 / R=%.1f", Rden));
      lat->SetNDC(); lat->SetTextSize(0.06); lat->SetTextAlign(22); lat->Draw();

      // Data
      TH1* hDataRR = (TH1*)fData->Get(Form("hDataXsecRatio_R%02d_R%02d", iRnum, iRden));
      if (hDataRR) {
        hDataRR = (TH1*)hDataRR->Clone(Form("dataRR_%d_%d", ip, gNN));
        hDataRR->SetDirectory(0);
        hDataRR->SetLineColor(kBlack); hDataRR->SetMarkerStyle(20); hDataRR->SetMarkerSize(0.6);
        hDataRR->Draw("pe same");
      }

      // Collect model R-ratios for double ratio
      struct ModelRR { TH1* hRatio; Color_t col; Int_t ls; Int_t lw; TString name; };
      std::vector<ModelRR> modelRRs;

      // Standalone models
      for (const auto& model : models) {
        if (!model.enabled) continue;
        TH1* hN = LoadModelYield(model, Rnum, isUE);
        TH1* hD = LoadModelYield(model, Rden, isUE);
        if (hN && hD) {
          TH1* hR = MakeRatio(hN, hD, Form("rrM_%s_%d_%d", model.name.Data(), ip, gNN));
          StyleRatio(hR, model.color, model.lineStyle, model.lineWidth);
          hR->Draw("hist e same");
          modelRRs.push_back({hR, model.color, model.lineStyle, model.lineWidth, model.legLabel});
        }
        delete hN; delete hD;
      }

      // POWHEG
      TH1* hPowRR = nullptr;
      if (iR02 < (int)powhegXsec.size() && j < (int)powhegXsec.size() &&
          powhegXsec[iR02] && powhegXsec[j]) {
        hPowRR = MakeRatio(powhegXsec[iR02], powhegXsec[j], Form("rrPW_%d_%d", ip, gNN));
        StyleRatio(hPowRR, kMagenta+1, 7, 2);
        hPowRR->Draw("hist e same");
        modelRRs.push_back({hPowRR, kMagenta+1, 7, 2, "POWHEG NLO"});
      }

      // MC truth INEL
      TH1* hMCTrr = nullptr;
      if (iR02 < (int)mcTruthYieldINEL.size() && j < (int)mcTruthYieldINEL.size() &&
          mcTruthYieldINEL[iR02] && mcTruthYieldINEL[j]) {
        hMCTrr = MakeRatio(mcTruthYieldINEL[iR02], mcTruthYieldINEL[j], Form("rrMCT_%d_%d", ip, gNN));
        StyleRatio(hMCTrr, kGray+1, 1, 2);
        hMCTrr->Draw("hist e same");
        modelRRs.push_back({hMCTrr, kGray+1, 1, 2, "MC truth INEL"});
      }

      if (hDataRR) hDataRR->Draw("pe same");  // redraw on top

      // Legend (first panel)
      if (ip == 0) {
        TLegend* leg = new TLegend(0.22, 0.55, 0.80, 0.88);
        leg->SetTextSize(0.04); leg->SetBorderSize(0); leg->SetFillStyle(0);
        if (hDataRR) leg->AddEntry(hDataRR, "ALICE data", "lpe");
        for (const auto& mr : modelRRs) {
          TH1F* dd = new TH1F(Form("dRR_%s_%d", mr.name.Data(), gNN), "", 1, 0, 1);
          dd->SetLineColor(mr.col); dd->SetLineStyle(mr.ls); dd->SetLineWidth(mr.lw);
          leg->AddEntry(dd, mr.name.Data(), "l");
        }
        leg->Draw();
      }

      // ---- Lower pad: Model/Data double ratio ----
      can->cd();
      TPad* pLo = new TPad(Form("pLoRR_%d_%d", ip, gNN), "", padL, 0, padR, splitY);
      pLo->SetTopMargin(0.0); pLo->SetBottomMargin(0.25);
      pLo->SetLeftMargin(lm); pLo->SetRightMargin(rm);
      pLo->SetTicks(1, 1); pLo->Draw(); pLo->cd();

      TH1F* frLo = new TH1F(Form("frRRlo_%d_%d", ip, gNN), "", 100, kPtMin, kPtMax);
      frLo->SetMinimum(0.75); frLo->SetMaximum(1.25);
      frLo->GetXaxis()->SetTitle("#it{p}_{T}^{ch} (GeV/#it{c})");
      frLo->GetXaxis()->SetTitleSize(0.08); frLo->GetXaxis()->SetLabelSize(0.07);
      frLo->GetXaxis()->SetTitleOffset(1.0);
      frLo->GetYaxis()->SetTitle(ip == 0 ? "Model / Data" : "");
      frLo->GetYaxis()->SetTitleSize(0.08); frLo->GetYaxis()->SetLabelSize(0.07);
      frLo->GetYaxis()->SetTitleOffset(0.7);
      frLo->GetYaxis()->SetNdivisions(505);
      frLo->Draw("AXIS");

      TLine* unity = new TLine(kPtMin, 1.0, kPtMax, 1.0);
      unity->SetLineStyle(2); unity->Draw();

      // Data at unity with errors
      if (hDataRR) {
        TH1* hDself = (TH1*)hDataRR->Clone(Form("dselfRR_%d_%d", ip, gNN));
        for (int ib = 1; ib <= hDself->GetNbinsX(); ++ib) {
          double v = hDself->GetBinContent(ib);
          double e = hDself->GetBinError(ib);
          if (v > 0) { hDself->SetBinContent(ib, 1.0); hDself->SetBinError(ib, e/v); }
          else { hDself->SetBinContent(ib, 0); hDself->SetBinError(ib, 0); }
        }
        hDself->SetMarkerStyle(20); hDself->SetMarkerSize(0.5);
        hDself->Draw("pe same");
      }

      // Model/Data double ratios
      for (const auto& mr : modelRRs) {
        if (!mr.hRatio || !hDataRR) continue;
        TH1* hDR = MakeRatio(mr.hRatio, hDataRR, Form("drRR_%s_%d_%d", mr.name.Data(), ip, gNN));
        if (!hDR) continue;
        hDR->SetLineColor(mr.col); hDR->SetLineStyle(mr.ls); hDR->SetLineWidth(mr.lw);
        hDR->SetMarkerStyle(0);
        hDR->Draw("hist e same");
      }

      can->cd();
    }

    TString pdf = Form("%s/RRatioComparison%s.pdf", outputDir, suffix.Data());
    can->Print(pdf);
    std::cout << "Saved: " << pdf << std::endl;
    } // end iR02 found
  }

  // ============================================================
  // Canvas 5: R=0.4 detailed comparison — CMS-style publication quality
  // Upper pad: per-event yield spectra (log scale)
  // Lower pad: Model / Data ratio
  // ============================================================
  {
    int iR04 = -1;
    for (int i = 0; i < nR; ++i) {
      if (TMath::Abs(rValues[i] - 0.4) < 0.01) { iR04 = i; break; }
    }
    if (iR04 < 0) {
      std::cout << "[Canvas5] R=0.4 not found in data, skipping" << std::endl;
    } else if (!dataYield[iR04]) {
      std::cout << "[Canvas5] No data yield for R=0.4, skipping" << std::endl;
    } else {

      // ---- Style ----
      gStyle->SetOptStat(0);
      gStyle->SetOptTitle(0);

      const int   kFont      = 42;
      const float kMarkerSz  = 0.8;
      const float splitY     = 0.40;    // ratio pad fraction (includes bottom margin for x-axis)
      const float lMargin    = 0.14;
      const float rMargin    = 0.035;

      // Effective text sizes (absolute on canvas, independent of pad height)
      const float absTitle   = 0.040;
      const float absLabel   = 0.035;
      const float absTex     = 0.033;
      const float absLeg     = 0.028;

      TCanvas* can = new TCanvas(Form("detailR04_%d", ++gNN), "R=0.4 Yield Comparison", 600, 700);
      can->SetFillColor(kWhite);
      can->cd();

      TPad* padUp = new TPad(Form("padUp_%d", gNN), "", 0, splitY, 1, 1);
      padUp->SetFillColor(kWhite);
      padUp->SetTopMargin(0.07);
      padUp->SetBottomMargin(0.0);
      padUp->SetLeftMargin(lMargin);
      padUp->SetRightMargin(rMargin);
      padUp->SetLogy();
      padUp->SetTicks(1, 1);
      padUp->Draw();

      TPad* padLo = new TPad(Form("padLo_%d", gNN), "", 0, 0, 1, splitY);
      padLo->SetFillColor(kWhite);
      padLo->SetTopMargin(0.0);
      padLo->SetBottomMargin(0.30);
      padLo->SetLeftMargin(lMargin);
      padLo->SetRightMargin(rMargin);
      padLo->SetTicks(1, 1);
      padLo->Draw();

      TH1* hData = dataYield[iR04];
      const double xMin = kPtMin;
      const double xMax = kPtMax;
      const float upH = 1.0 - splitY;
      const float loH = splitY;

      // ==== UPPER PAD ====
      padUp->cd();

      TH1F* frUp = new TH1F(Form("frUp_%d", gNN), "", 1, xMin, xMax);
      frUp->SetMinimum(1e-9);
      frUp->SetMaximum(2e-2);
      frUp->GetYaxis()->SetTitle("#frac{1}{#it{N}_{evt}} #frac{d^{2}#it{N}}{d#it{p}_{T}d#it{#eta}} (GeV/#it{c})^{#minus1}");
      frUp->GetYaxis()->SetTitleFont(kFont);
      frUp->GetYaxis()->SetTitleSize(absTitle / upH / 1.5);
      frUp->GetYaxis()->SetTitleOffset(2.5 * upH);
      frUp->GetYaxis()->SetLabelFont(kFont);
      frUp->GetYaxis()->SetLabelSize(absLabel / upH);
      frUp->GetXaxis()->SetLabelSize(0);
      frUp->GetXaxis()->SetTitleSize(0);
      frUp->GetYaxis()->SetNdivisions(510);
      frUp->Draw("AXIS");

      // Collect models
      struct DrawnModel { TH1* hYield; TString name; Color_t col; Int_t ls; Int_t lw; };
      std::vector<DrawnModel> drawn;

      // MC truth selected (dashed)
      if (mcTruthYield[iR04]) {
        TH1* h = (TH1*)mcTruthYield[iR04]->Clone(Form("mctSelR04_%d", gNN));
        h->SetLineColor(kGray+1); h->SetLineStyle(2); h->SetLineWidth(2); h->SetMarkerStyle(0);
        h->Draw("hist e same");
        drawn.push_back({h, "MC truth sel", kGray+1, 2, 2});
      }

      // MC truth INEL (solid)
      if (mcTruthYieldINEL[iR04]) {
        TH1* h = (TH1*)mcTruthYieldINEL[iR04]->Clone(Form("mctInelR04_%d", gNN));
        h->SetLineColor(kGray+1); h->SetLineStyle(1); h->SetLineWidth(2); h->SetMarkerStyle(0);
        h->Draw("hist e same");
        drawn.push_back({h, "MC truth INEL", kGray+1, 1, 2});
      }

      // GenMC (HY)
      if (genMCYield.count(iR04)) {
        for (const auto& p : genMCYield[iR04]) {
          TH1* h = (TH1*)p.second->Clone(Form("gmcR04_%s_%d", p.first->name.Data(), gNN));
          h->SetLineColor(p.first->color); h->SetLineStyle(p.first->lineStyle);
          h->SetLineWidth(p.first->lineWidth); h->SetMarkerStyle(0);
          h->Draw("hist e same");
          drawn.push_back({h, p.first->name, p.first->color, p.first->lineStyle, p.first->lineWidth});
        }
      }

      // KIAF standalone
      for (const auto& model : models) {
        if (!model.enabled) continue;
        TH1* h = LoadModelYield(model, 0.4, isUE);
        if (!h) continue;
        h->SetLineColor(model.color); h->SetLineStyle(model.lineStyle);
        h->SetLineWidth(model.lineWidth); h->SetMarkerStyle(0);
        h->Draw("hist e same");
        drawn.push_back({h, model.legLabel, model.color, model.lineStyle, model.lineWidth});
      }

      // POWHEG
      if (powhegXsec[iR04]) {
        TH1* h = (TH1*)powhegXsec[iR04]->Clone(Form("pwR04_%d", gNN));
        h->Scale(1.0 / kSigmaINEL);
        h->SetLineColor(kMagenta+1); h->SetLineStyle(7); h->SetLineWidth(2); h->SetMarkerStyle(0);
        h->Draw("hist e same");
        drawn.push_back({h, "POWHEG NLO (Hadi)", kMagenta+1, 7, 2});
      }

      // Data on top
      TH1* hDataDraw = (TH1*)hData->Clone(Form("dataR04_%d", gNN));
      hDataDraw->SetLineColor(kBlack); hDataDraw->SetMarkerColor(kBlack);
      hDataDraw->SetMarkerStyle(20); hDataDraw->SetMarkerSize(kMarkerSz);
      hDataDraw->SetLineWidth(1);
      hDataDraw->GetXaxis()->SetRangeUser(xMin, xMax);
      hDataDraw->Draw("pe same");

      // Physics info — top-right, away from data
      TLatex tex;
      tex.SetNDC(); tex.SetTextFont(42); tex.SetTextAlign(11);  // right-aligned

      tex.SetTextSize(absTex / upH * 1.1); tex.SetTextFont(62);
      tex.DrawLatex(0.20, 0.85, "ALICE WIP");
      tex.SetTextFont(42);

      tex.SetTextSize(absTex / upH * 0.9);
      tex.DrawLatex(0.16, 0.17, "pp, #sqrt{#it{s}} = 13.6 TeV");
      tex.DrawLatex(0.16, 0.11, "#it{p}_{T,track} > 0.15 GeV/#it{c}, |#it{#eta}_{track}| < 0.9");
      tex.DrawLatex(0.16, 0.05, "Charged jets, anti-#it{k}_{T}, #it{R} = 0.4, |#it{#eta}_{jet}| < 0.5");

      // Legend — position and style from markerLegend.c (user-tuned)
      TLegend* leg = new TLegend(0.45, 0.415, 0.95, 0.896);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.042);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);

      TLegendEntry* le;
      le = leg->AddEntry(hDataDraw, "ALICE data", "pe");
      le->SetMarkerStyle(20); le->SetMarkerSize(0.8); le->SetTextFont(42);
      for (const auto& dm : drawn) {
        TH1F* dum = new TH1F(Form("dR04_%s_%d", dm.name.Data(), gNN), "", 1, 0, 1);
        dum->SetLineColor(dm.col); dum->SetLineStyle(dm.ls); dum->SetLineWidth(dm.lw);
        le = leg->AddEntry(dum, dm.name.Data(), "l");
        le->SetLineColor(dm.col); le->SetLineStyle(dm.ls); le->SetLineWidth(dm.lw);
        le->SetTextFont(42);
      }
      leg->Draw();

      // ==== LOWER PAD ====
      padLo->cd();

      TH1F* frLo = new TH1F(Form("frLo_%d", gNN), "", 1, xMin, xMax);
      frLo->SetMinimum(0.45);
      frLo->SetMaximum(1.95);
      frLo->GetXaxis()->SetTitle("#it{p}_{T,jet}^{ch} (GeV/#it{c})");
      frLo->GetXaxis()->SetTitleFont(kFont);
      frLo->GetXaxis()->SetTitleSize(absTitle / loH);
      frLo->GetXaxis()->SetTitleOffset(1.0);
      frLo->GetXaxis()->SetLabelFont(kFont);
      frLo->GetXaxis()->SetLabelSize(absLabel / loH);
      frLo->GetYaxis()->SetTitle("Model / Data");
      frLo->GetYaxis()->SetTitleFont(kFont);
      frLo->GetYaxis()->SetTitleSize(absTitle / loH);
      frLo->GetYaxis()->SetTitleOffset(1.25 * loH);
      frLo->GetYaxis()->SetLabelFont(kFont);
      frLo->GetYaxis()->SetLabelSize(absLabel / loH);
      frLo->GetYaxis()->SetNdivisions(505);
      frLo->GetYaxis()->CenterTitle(true);
      frLo->Draw("AXIS");

      // Unity line
      TLine* unity = new TLine(xMin, 1.0, xMax, 1.0);
      unity->SetLineColor(kBlack); unity->SetLineStyle(2); unity->SetLineWidth(1);
      unity->Draw();

      // Data stat errors at unity
      TH1* hDself = (TH1*)hData->Clone(Form("dselfR04_%d", gNN));
      for (int ib = 1; ib <= hDself->GetNbinsX(); ++ib) {
        double val = hDself->GetBinContent(ib);
        double err = hDself->GetBinError(ib);
        if (val > 0) { hDself->SetBinContent(ib, 1.0); hDself->SetBinError(ib, err / val); }
        else { hDself->SetBinContent(ib, 0); hDself->SetBinError(ib, 0); }
      }
      hDself->SetLineColor(kBlack); hDself->SetMarkerColor(kBlack);
      hDself->SetMarkerStyle(20); hDself->SetMarkerSize(kMarkerSz);
      hDself->GetXaxis()->SetRangeUser(xMin, xMax);
      hDself->Draw("pe same");

      // Model ratios
      for (const auto& dm : drawn) {
        TH1* hR = MakeRatio(dm.hYield, hData, Form("ratR04_%s_%d", dm.name.Data(), gNN));
        if (!hR) continue;
        hR->SetLineColor(dm.col); hR->SetLineStyle(dm.ls); hR->SetLineWidth(dm.lw);
        hR->SetMarkerStyle(0);
        hR->Draw("hist e same");
      }
      hDself->Draw("pe same");

      TString pdf = Form("%s/DetailComparison_Yield_R04%s.pdf", outputDir, suffix.Data());
      can->Print(pdf);
      std::cout << "Saved: " << pdf << std::endl;
    }
  }

  // ============================================================
  // Canvas 6: R=0.4 cross-section detail — same layout as Canvas 5
  // Upper pad: d²σ/(dpT dη) spectra, Lower pad: Model/Data ratio
  // ============================================================
  {
    int iR04 = -1;
    for (int i = 0; i < nR; ++i) {
      if (TMath::Abs(rValues[i] - 0.4) < 0.01) { iR04 = i; break; }
    }
    if (iR04 >= 0 && dataXsec[iR04]) {

      gStyle->SetOptStat(0);
      gStyle->SetOptTitle(0);

      const int   kFont      = 42;
      const float kMarkerSz  = 0.8;
      const float splitY     = 0.40;
      const float lMargin    = 0.14;
      const float rMargin    = 0.035;
      const float absTitle   = 0.040;
      const float absLabel   = 0.035;
      const float absTex     = 0.033;

      TCanvas* can = new TCanvas(Form("detailXsR04_%d", ++gNN), "R=0.4 Xsec Comparison", 600, 700);
      can->SetFillColor(kWhite);
      can->cd();

      TPad* padUp = new TPad(Form("padUpXs_%d", gNN), "", 0, splitY, 1, 1);
      padUp->SetFillColor(kWhite);
      padUp->SetTopMargin(0.07);
      padUp->SetBottomMargin(0.0);
      padUp->SetLeftMargin(lMargin);
      padUp->SetRightMargin(rMargin);
      padUp->SetLogy();
      padUp->SetTicks(1, 1);
      padUp->Draw();

      TPad* padLo = new TPad(Form("padLoXs_%d", gNN), "", 0, 0, 1, splitY);
      padLo->SetFillColor(kWhite);
      padLo->SetTopMargin(0.0);
      padLo->SetBottomMargin(0.30);
      padLo->SetLeftMargin(lMargin);
      padLo->SetRightMargin(rMargin);
      padLo->SetTicks(1, 1);
      padLo->Draw();

      TH1* hData = dataXsec[iR04];
      const double xMin = kPtMin;
      const double xMax = kPtMax;
      const float upH = 1.0 - splitY;
      const float loH = splitY;

      // ==== UPPER PAD ====
      padUp->cd();

      TH1F* frUp = new TH1F(Form("frUpXs_%d", gNN), "", 1, xMin, xMax);
      frUp->SetMinimum(1e-7);
      frUp->SetMaximum(1e1);
      frUp->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#it{#eta}} (mb GeV^{#minus1}#it{c})");
      frUp->GetYaxis()->SetTitleFont(kFont);
      frUp->GetYaxis()->SetTitleSize(absTitle / upH / 1.5);
      frUp->GetYaxis()->SetTitleOffset(2.5 * upH);
      frUp->GetYaxis()->SetLabelFont(kFont);
      frUp->GetYaxis()->SetLabelSize(absLabel / upH);
      frUp->GetXaxis()->SetLabelSize(0);
      frUp->GetXaxis()->SetTitleSize(0);
      frUp->GetYaxis()->SetNdivisions(510);
      frUp->Draw("AXIS");

      // Models: yield × σ_INEL → cross-section
      struct DrawnModel { TH1* hXsec; TString name; Color_t col; Int_t ls; Int_t lw; };
      std::vector<DrawnModel> drawnXs;

      // MC truth (cross-section, bin 0.5 normalization)
      if (mcTruthXsec[iR04]) {
        TH1* h = (TH1*)mcTruthXsec[iR04]->Clone(Form("mctXsR04_%d", gNN));
        h->SetLineColor(kGray+1); h->SetLineStyle(1); h->SetLineWidth(2); h->SetMarkerStyle(0);
        h->Draw("hist e same");
        drawnXs.push_back({h, "MC truth INEL (PYTHIA8 Monash)", kGray+1, 1, 2});
      }

      // GenMC (HY): yield × σ_INEL
      if (genMCYield.count(iR04)) {
        for (const auto& p : genMCYield[iR04]) {
          TH1* h = (TH1*)p.second->Clone(Form("gmcXsR04_%s_%d", p.first->name.Data(), gNN));
          h->Scale(kSigmaINEL);
          h->SetLineColor(p.first->color); h->SetLineStyle(p.first->lineStyle);
          h->SetLineWidth(p.first->lineWidth); h->SetMarkerStyle(0);
          h->Draw("hist e same");
          drawnXs.push_back({h, p.first->name, p.first->color, p.first->lineStyle, p.first->lineWidth});
        }
      }

      // KIAF standalone (cross-section directly)
      for (const auto& model : models) {
        if (!model.enabled) continue;
        TH1* h = LoadModelXsec(model, 0.4, isUE);
        if (!h) continue;
        h->SetLineColor(model.color); h->SetLineStyle(model.lineStyle);
        h->SetLineWidth(model.lineWidth); h->SetMarkerStyle(0);
        h->Draw("hist e same");
        drawnXs.push_back({h, model.legLabel, model.color, model.lineStyle, model.lineWidth});
      }

      // POWHEG (already cross-section)
      if (powhegXsec[iR04]) {
        TH1* h = (TH1*)powhegXsec[iR04]->Clone(Form("pwXsR04_%d", gNN));
        h->SetLineColor(kMagenta+1); h->SetLineStyle(7); h->SetLineWidth(2); h->SetMarkerStyle(0);
        h->Draw("hist e same");
        drawnXs.push_back({h, "POWHEG NLO (Hadi)", kMagenta+1, 7, 2});
      }

      // Data on top
      TH1* hDataDraw = (TH1*)hData->Clone(Form("dataXsR04_%d", gNN));
      hDataDraw->SetLineColor(kBlack); hDataDraw->SetMarkerColor(kBlack);
      hDataDraw->SetMarkerStyle(20); hDataDraw->SetMarkerSize(kMarkerSz);
      hDataDraw->SetLineWidth(1);
      hDataDraw->GetXaxis()->SetRangeUser(xMin, xMax);
      hDataDraw->Draw("pe same");

      // Physics info
      TLatex tex;
      tex.SetNDC(); tex.SetTextFont(42); tex.SetTextAlign(11);
      tex.SetTextSize(absTex / upH * 1.1); tex.SetTextFont(62);
      tex.DrawLatex(0.20, 0.85, "ALICE WIP");
      tex.SetTextFont(42);
      tex.SetTextSize(absTex / upH * 0.9);
      tex.DrawLatex(0.16, 0.17, "pp, #sqrt{#it{s}} = 13.6 TeV");
      tex.DrawLatex(0.16, 0.11, "#it{p}_{T,track} > 0.15 GeV/#it{c}, |#it{#eta}_{track}| < 0.9");
      tex.DrawLatex(0.16, 0.05, "Charged jets, anti-#it{k}_{T}, #it{R} = 0.4, |#it{#eta}_{jet}| < 0.5");

      // Legend
      TLegend* leg = new TLegend(0.45, 0.415, 0.95, 0.896);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.042);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);

      TLegendEntry* le;
      le = leg->AddEntry(hDataDraw, "ALICE data", "pe");
      le->SetMarkerStyle(20); le->SetMarkerSize(0.8); le->SetTextFont(42);
      for (const auto& dm : drawnXs) {
        TH1F* dum = new TH1F(Form("dXsR04_%s_%d", dm.name.Data(), gNN), "", 1, 0, 1);
        dum->SetLineColor(dm.col); dum->SetLineStyle(dm.ls); dum->SetLineWidth(dm.lw);
        le = leg->AddEntry(dum, dm.name.Data(), "l");
        le->SetLineColor(dm.col); le->SetLineStyle(dm.ls); le->SetLineWidth(dm.lw);
        le->SetTextFont(42);
      }
      leg->Draw();

      // ==== LOWER PAD ====
      padLo->cd();

      TH1F* frLo = new TH1F(Form("frLoXs_%d", gNN), "", 1, xMin, xMax);
      frLo->SetMinimum(0.45);
      frLo->SetMaximum(1.95);
      frLo->GetXaxis()->SetTitle("#it{p}_{T,jet}^{ch} (GeV/#it{c})");
      frLo->GetXaxis()->SetTitleFont(kFont);
      frLo->GetXaxis()->SetTitleSize(absTitle / loH);
      frLo->GetXaxis()->SetTitleOffset(1.0);
      frLo->GetXaxis()->SetLabelFont(kFont);
      frLo->GetXaxis()->SetLabelSize(absLabel / loH);
      frLo->GetYaxis()->SetTitle("Model / Data");
      frLo->GetYaxis()->SetTitleFont(kFont);
      frLo->GetYaxis()->SetTitleSize(absTitle / loH);
      frLo->GetYaxis()->SetTitleOffset(1.25 * loH);
      frLo->GetYaxis()->SetLabelFont(kFont);
      frLo->GetYaxis()->SetLabelSize(absLabel / loH);
      frLo->GetYaxis()->SetNdivisions(505);
      frLo->GetYaxis()->CenterTitle(true);
      frLo->Draw("AXIS");

      TLine* unity = new TLine(xMin, 1.0, xMax, 1.0);
      unity->SetLineColor(kBlack); unity->SetLineStyle(2); unity->SetLineWidth(1);
      unity->Draw();

      // (no ±10% band — systematic uncertainties not yet included)

      // Data at unity
      TH1* hDself = (TH1*)hData->Clone(Form("dselfXsR04_%d", gNN));
      for (int ib = 1; ib <= hDself->GetNbinsX(); ++ib) {
        double val = hDself->GetBinContent(ib);
        double err = hDself->GetBinError(ib);
        if (val > 0) { hDself->SetBinContent(ib, 1.0); hDself->SetBinError(ib, err / val); }
        else { hDself->SetBinContent(ib, 0); hDself->SetBinError(ib, 0); }
      }
      hDself->SetLineColor(kBlack); hDself->SetMarkerColor(kBlack);
      hDself->SetMarkerStyle(20); hDself->SetMarkerSize(kMarkerSz);
      hDself->GetXaxis()->SetRangeUser(xMin, xMax);
      hDself->Draw("pe same");

      // Model ratios
      for (const auto& dm : drawnXs) {
        TH1* hR = MakeRatio(dm.hXsec, hData, Form("ratXsR04_%s_%d", dm.name.Data(), gNN));
        if (!hR) continue;
        hR->SetLineColor(dm.col); hR->SetLineStyle(dm.ls); hR->SetLineWidth(dm.lw);
        hR->SetMarkerStyle(0);
        hR->Draw("hist e same");
      }
      hDself->Draw("pe same");

      TString pdf = Form("%s/DetailComparison_Xsec_R04%s.pdf", outputDir, suffix.Data());
      can->Print(pdf);
      std::cout << "Saved: " << pdf << std::endl;
    }
  }

  // ============================================================
  // Canvas 7: R=0.4 Model / MC truth ratio (cross-section, ratio pad only)
  // ============================================================
  {
    int iR04 = -1;
    for (int i = 0; i < nR; ++i) {
      if (TMath::Abs(rValues[i] - 0.4) < 0.01) { iR04 = i; break; }
    }
    if (iR04 >= 0 && mcTruthXsec[iR04]) {

      const int kFont = 42;
      const float absTitle = 0.045;
      const float absLabel = 0.040;
      const float absTex   = 0.038;

      TCanvas* can = new TCanvas(Form("modelVsMCT_%d", ++gNN), "R=0.4 Model/MCtruth", 600, 400);
      can->SetFillColor(kWhite);
      can->SetLeftMargin(0.14);
      can->SetRightMargin(0.035);
      can->SetTopMargin(0.06);
      can->SetBottomMargin(0.15);
      can->SetTicks(1, 1);

      TH1* hDenom = mcTruthXsec[iR04];
      const double xMin = kPtMin;
      const double xMax = kPtMax;

      TH1F* frame = new TH1F(Form("frMvMCT_%d", gNN), "", 1, xMin, xMax);
      frame->SetMinimum(0.3);
      frame->SetMaximum(2.2);
      frame->GetXaxis()->SetTitle("#it{p}_{T,jet}^{ch} (GeV/#it{c})");
      frame->GetXaxis()->SetTitleFont(kFont);
      frame->GetXaxis()->SetTitleSize(absTitle);
      frame->GetXaxis()->SetLabelFont(kFont);
      frame->GetXaxis()->SetLabelSize(absLabel);
      frame->GetYaxis()->SetTitle("Model / MC truth");
      frame->GetYaxis()->SetTitleFont(kFont);
      frame->GetYaxis()->SetTitleSize(absTitle);
      frame->GetYaxis()->SetTitleOffset(1.1);
      frame->GetYaxis()->SetLabelFont(kFont);
      frame->GetYaxis()->SetLabelSize(absLabel);
      frame->GetYaxis()->SetNdivisions(505);
      frame->GetYaxis()->CenterTitle(true);
      frame->Draw("AXIS");

      TLine* unity = new TLine(xMin, 1.0, xMax, 1.0);
      unity->SetLineColor(kBlack); unity->SetLineStyle(2); unity->SetLineWidth(1);
      unity->Draw();

      // MC truth / MC truth = 1 with stat error bars (no band — bands = syst only in ALICE)
      TH1* hMCTself = (TH1*)hDenom->Clone(Form("mctSelf_%d", gNN));
      for (int ib = 1; ib <= hMCTself->GetNbinsX(); ++ib) {
        double val = hMCTself->GetBinContent(ib);
        double err = hMCTself->GetBinError(ib);
        if (val > 0) { hMCTself->SetBinContent(ib, 1.0); hMCTself->SetBinError(ib, err / val); }
        else { hMCTself->SetBinContent(ib, 0); hMCTself->SetBinError(ib, 0); }
      }
      hMCTself->SetLineColor(kGray+2); hMCTself->SetLineStyle(7); hMCTself->SetLineWidth(2);
      hMCTself->SetMarkerStyle(0);
      hMCTself->Draw("hist e same");

      struct DrawnRat { TString name; Color_t col; Int_t ls; Int_t lw; };
      std::vector<DrawnRat> legEntries;

      // GenMC (HY) / MC truth
      if (genMCYield.count(iR04)) {
        for (const auto& p : genMCYield[iR04]) {
          TH1* hXs = (TH1*)p.second->Clone(Form("gmcXsMCT_%s_%d", p.first->name.Data(), gNN));
          hXs->Scale(kSigmaINEL);
          TH1* hR = MakeRatio(hXs, hDenom, Form("ratMCT_gmc_%s_%d", p.first->name.Data(), gNN));
          if (hR) {
            hR->SetLineColor(p.first->color); hR->SetLineStyle(p.first->lineStyle);
            hR->SetLineWidth(p.first->lineWidth); hR->SetMarkerStyle(0);
            hR->Draw("hist e same");
            legEntries.push_back({p.first->name, p.first->color, p.first->lineStyle, p.first->lineWidth});
          }
        }
      }

      // KIAF standalone / MC truth
      for (const auto& model : models) {
        if (!model.enabled) continue;
        TH1* hModel = LoadModelXsec(model, 0.4, isUE);
        if (!hModel) continue;
        TH1* hR = MakeRatio(hModel, hDenom, Form("ratMCT_m_%s_%d", model.name.Data(), gNN));
        if (hR) {
          hR->SetLineColor(model.color); hR->SetLineStyle(model.lineStyle);
          hR->SetLineWidth(model.lineWidth); hR->SetMarkerStyle(0);
          hR->Draw("hist e same");
          legEntries.push_back({model.legLabel, model.color, model.lineStyle, model.lineWidth});
        }
        delete hModel;
      }

      // POWHEG / MC truth
      if (powhegXsec[iR04]) {
        TH1* hR = MakeRatio(powhegXsec[iR04], hDenom, Form("ratMCT_pw_%d", gNN));
        if (hR) {
          hR->SetLineColor(kMagenta+1); hR->SetLineStyle(7); hR->SetLineWidth(2);
          hR->SetMarkerStyle(0);
          hR->Draw("hist e same");
          legEntries.push_back({"POWHEG NLO (Hadi)", kMagenta+1, 7, 2});
        }
      }

      // Data / MC truth
      if (dataXsec[iR04]) {
        TH1* hR = MakeRatio(dataXsec[iR04], hDenom, Form("ratMCT_data_%d", gNN));
        if (hR) {
          hR->SetLineColor(kBlack); hR->SetMarkerColor(kBlack);
          hR->SetMarkerStyle(20); hR->SetMarkerSize(0.8);
          hR->Draw("pe same");
        }
      }

      // Physics info
      TLatex tex;
      tex.SetNDC(); tex.SetTextFont(62); tex.SetTextAlign(11);
      tex.SetTextSize(absTex * 1.1);
      tex.DrawLatex(0.16, 0.90, "ALICE WIP");
      tex.SetTextFont(42); tex.SetTextSize(absTex * 0.85);
      tex.DrawLatex(0.16, 0.84, "pp #sqrt{#it{s}} = 13.6 TeV, charged jets, #it{R} = 0.4");

      // Legend
      TLegend* leg = new TLegend(0.45, 0.55, 0.96, 0.93);
      leg->SetBorderSize(0); leg->SetTextFont(42); leg->SetTextSize(0.035);
      leg->SetFillStyle(0);

      TH1F* dumD = new TH1F(Form("dumD_%d", gNN), "", 1, 0, 1);
      dumD->SetMarkerStyle(20); dumD->SetMarkerSize(0.8);
      leg->AddEntry(dumD, "ALICE data", "pe");

      TH1F* dumMCT = new TH1F(Form("dumMCT_%d", gNN), "", 1, 0, 1);
      dumMCT->SetLineColor(kGray+2); dumMCT->SetLineStyle(7); dumMCT->SetLineWidth(2);
      leg->AddEntry(dumMCT, "MC truth INEL (PYTHIA8 Monash)", "l");

      for (const auto& le : legEntries) {
        TH1F* dum = new TH1F(Form("dumMvMCT_%s_%d", le.name.Data(), gNN), "", 1, 0, 1);
        dum->SetLineColor(le.col); dum->SetLineStyle(le.ls); dum->SetLineWidth(le.lw);
        leg->AddEntry(dum, le.name.Data(), "l");
      }
      leg->Draw();

      TString pdf = Form("%s/ModelVsMCtruth_Xsec_R04%s.pdf", outputDir, suffix.Data());
      can->Print(pdf);
      std::cout << "Saved: " << pdf << std::endl;
    }
  }

  fData->Close();
  delete fData;

  std::cout << "\n=== ModelComparisons complete ===" << std::endl;
  std::cout << "Output: " << outputDir << std::endl;
}
