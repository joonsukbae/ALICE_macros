#include "Filipad2.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TString.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

// ============================================================
//  TrackTuner QA: Compare 3 tuner variation output files
//  - DCAonly: DCA correction only (no q/pT smearing)
//  - constPtSmear1p5: Constant pT smearing ratio = 1.5
//  - customTuner: Graph-based tuner (Sigma1overPt file)
// ============================================================

// --- Plot toggles ---
const bool kDrawSigmaPtOverPt       = true;
const bool kDrawSigma1OverPt        = true;
const bool kDrawPtResidualSigma     = true;
const bool kDrawDCAxyWidth          = false;
const bool kDrawDCAzWidth           = false;
const bool kDrawDCAxyMean           = false;
const bool kDrawDCAzMean            = false;
const bool kDrawJetPtDet            = true;
const bool kDrawJetPtUESub          = false;
const bool kDrawJetPtParticle       = true;
const bool kDrawDeltaPtSigma        = false;
const bool kDrawPtResidualDist     = true;
const bool kDrawTrackPt            = true;
const bool kDrawTrackEta           = true;
const bool kDrawTrackPhi           = true;

// --- Source definition ---
struct Source {
  TString file;
  TString label;
  Color_t color;
  int marker;
};

const TString kBaseDir = "../../jets/AnalysisResults";
const TString kOutputDir = "plots/tunerComparison_2023MB";

std::vector<Source> kSources = {
  {kBaseDir + "/633426_AnalysisResults.root",  "2023 MB MC, tuner off",              kBlue,     20},
  {kBaseDir + "/639756_AnalysisResults.root",  "2023 MB MC, DCA only",              kCyan+2,   21},
  {kBaseDir + "/637148_AnalysisResults.root",  "2023 MB MC, DCA + map p_{T} smear", kRed,      22},
  {kBaseDir + "/639655_AnalysisResults.root",  "2023 MB MC, DCA + const 1.5 smear", kGreen+2, 33},
};

const int kRefIndex = 0; // ratio reference: DCAonly

// --- Data overlay (set kDrawData=true and file path to overlay data) ---
const bool kDrawData = true;
Source kDataSource = {
  kBaseDir + "/616559_AnalysisResults.root",
  "2023 Data",
  kBlack,
  24  // open circle
};

// --- Jet pT binning (from DrawJetsMCRDependentHelpers.h) ---
const Double_t kJetPtBinsReco[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
const int kNJetPtBinsReco = 25;
const Double_t kJetPtBinsRecoUE[] = {-200, -100, -50, -20, -10, -5, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
const int kNJetPtBinsRecoUE = 33;
const Double_t kJetPtBinsTruth[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
const int kNJetPtBinsTruth = 25;

// --- Histogram directories ---
const char* kDirTrackEff    = "track-efficiency";
const char* kDirPropagation = "propagation-service";
const char* kDirJetSpectra  = "jet-spectra-charged";

// --- Gaussian fit helper ---
struct FitResult {
  double mean    = 0;
  double meanErr = 0;
  double sigma   = 0;
  double sigmaErr= 0;
  bool ok = false;
};

FitResult FitCoreGaussian(TH1* h, double nSigma = 2.0) {
  FitResult r;
  if (!h || h->GetEffectiveEntries() < 20) return r;
  double mu0  = h->GetMean();
  double rms0 = h->GetRMS();
  if (rms0 <= 0) return r;
  for (int iter = 0; iter < 3; iter++) {
    double lo = mu0 - nSigma * rms0;
    double hi = mu0 + nSigma * rms0;
    TF1* f = new TF1("ftmp", "gaus", lo, hi);
    f->SetParameters(h->GetMaximum(), mu0, rms0);
    int st = h->Fit(f, "QNR0");
    if (st == 0) {
      mu0  = f->GetParameter(1);
      rms0 = std::abs(f->GetParameter(2));
      r.mean     = mu0;
      r.meanErr  = f->GetParError(1);
      r.sigma    = rms0;
      r.sigmaErr = f->GetParError(2);
      r.ok = true;
    }
    delete f;
    if (!r.ok) break;
  }
  return r;
}

// --- Extract mean profile from TH2 (Y mean per X bin) ---
// minEntries=1 for covariance matrix (each entry is one track's stored sigma)
TGraphErrors* ExtractMeanProfile(TH2* h2, const TString& name, double minEntries = 1) {
  if (!h2) return nullptr;
  TGraphErrors* g = new TGraphErrors();
  g->SetName(name);
  int np = 0;
  for (int ix = 1; ix <= h2->GetNbinsX(); ix++) {
    TH1D* py = h2->ProjectionY(Form("%s_py%d", name.Data(), ix), ix, ix, "e");
    if (py->GetEffectiveEntries() < minEntries) { delete py; continue; }
    double x  = h2->GetXaxis()->GetBinCenter(ix);
    double ex = 0.5 * h2->GetXaxis()->GetBinWidth(ix);
    double y  = py->GetMean();
    double ey = py->GetMeanError();
    g->SetPoint(np, x, y);
    g->SetPointError(np, ex, ey);
    np++;
    delete py;
  }
  cout << "  " << name << ": " << np << " points, x=["
       << (np>0 ? Form("%.1f", g->GetX()[0]) : "N/A") << ", "
       << (np>0 ? Form("%.1f", g->GetX()[np-1]) : "N/A") << "]" << endl;
  return g;
}

// --- Extract mean profile with variable bin merging ---
// mergeScheme: list of {ptThreshold, nMerge} — below first threshold no merge
TGraphErrors* ExtractMeanProfileVarBins(TH2* h2, const TString& name,
    const std::vector<std::pair<double,int>>& mergeScheme, double minEntries = 1) {
  if (!h2) return nullptr;
  TGraphErrors* g = new TGraphErrors();
  g->SetName(name);
  int np = 0;
  int ix = 1;
  int nbins = h2->GetNbinsX();
  while (ix <= nbins) {
    double ptLow = h2->GetXaxis()->GetBinLowEdge(ix);
    int nMerge = 1;
    for (const auto& ms : mergeScheme) {
      if (ptLow >= ms.first) nMerge = ms.second;
    }
    int ixEnd = std::min(ix + nMerge - 1, nbins);
    TH1D* py = h2->ProjectionY(Form("%s_vb%d", name.Data(), ix), ix, ixEnd, "e");
    if (py->GetEffectiveEntries() >= minEntries) {
      double xLo = h2->GetXaxis()->GetBinLowEdge(ix);
      double xHi = h2->GetXaxis()->GetBinUpEdge(ixEnd);
      double x = 0.5 * (xLo + xHi);
      double ex = 0.5 * (xHi - xLo);
      g->SetPoint(np, x, py->GetMean());
      g->SetPointError(np, ex, py->GetMeanError());
      np++;
    }
    delete py;
    ix = ixEnd + 1;
  }
  cout << "  " << name << ": " << np << " points, x=["
       << (np>0 ? Form("%.1f", g->GetX()[0]) : "N/A") << ", "
       << (np>0 ? Form("%.1f", g->GetX()[np-1]) : "N/A") << "]" << endl;
  return g;
}

// --- Extract sigma from Gaussian fit per pT bin of TH2 ---
TGraphErrors* ExtractSigmaVsPt(TH2* h2, const TString& name, double minEntries = 10) {
  if (!h2) return nullptr;
  TGraphErrors* g = new TGraphErrors();
  g->SetName(name);
  int np = 0;
  for (int ix = 1; ix <= h2->GetNbinsX(); ix++) {
    TH1D* py = h2->ProjectionY(Form("%s_sig%d", name.Data(), ix), ix, ix, "e");
    if (py->GetEffectiveEntries() < minEntries) { delete py; continue; }
    FitResult fr = FitCoreGaussian(py);
    if (!fr.ok) { delete py; continue; }
    double x  = h2->GetXaxis()->GetBinCenter(ix);
    double ex = 0.5 * h2->GetXaxis()->GetBinWidth(ix);
    g->SetPoint(np, x, fr.sigma);
    g->SetPointError(np, ex, fr.sigmaErr);
    np++;
    delete py;
  }
  cout << "  " << name << ": " << np << " points, x=["
       << (np>0 ? Form("%.1f", g->GetX()[0]) : "N/A") << ", "
       << (np>0 ? Form("%.1f", g->GetX()[np-1]) : "N/A") << "]" << endl;
  return g;
}

// --- DCA histograms: x=DCA, y=pT ---
// Iterate over Y-bins (pT), project X (DCA), extract mean
TGraphErrors* ExtractDCAMeanVsPt(TH2* h2, const TString& name, double minEntries = 1) {
  if (!h2) return nullptr;
  TGraphErrors* g = new TGraphErrors();
  g->SetName(name);
  int np = 0;
  for (int iy = 1; iy <= h2->GetNbinsY(); iy++) {
    TH1D* px = h2->ProjectionX(Form("%s_dcam%d", name.Data(), iy), iy, iy, "e");
    if (px->GetEffectiveEntries() < minEntries) { delete px; continue; }
    double x  = h2->GetYaxis()->GetBinCenter(iy);  // pT
    double ex = 0.5 * h2->GetYaxis()->GetBinWidth(iy);
    double y  = px->GetMean();   // mean DCA
    double ey = px->GetMeanError();
    g->SetPoint(np, x, y);
    g->SetPointError(np, ex, ey);
    np++;
    delete px;
  }
  cout << "  " << name << ": " << np << " points, x=["
       << (np>0 ? Form("%.1f", g->GetX()[0]) : "N/A") << ", "
       << (np>0 ? Form("%.1f", g->GetX()[np-1]) : "N/A") << "]" << endl;
  return g;
}

// Iterate over Y-bins (pT), project X (DCA), fit Gaussian for sigma
TGraphErrors* ExtractDCASigmaVsPt(TH2* h2, const TString& name, double minEntries = 10) {
  if (!h2) return nullptr;
  TGraphErrors* g = new TGraphErrors();
  g->SetName(name);
  int np = 0;
  for (int iy = 1; iy <= h2->GetNbinsY(); iy++) {
    TH1D* px = h2->ProjectionX(Form("%s_dcas%d", name.Data(), iy), iy, iy, "e");
    if (px->GetEffectiveEntries() < minEntries) { delete px; continue; }
    FitResult fr = FitCoreGaussian(px);
    if (!fr.ok || fr.sigma <= 0) { delete px; continue; }
    double x  = h2->GetYaxis()->GetBinCenter(iy);  // pT
    double ex = 0.5 * h2->GetYaxis()->GetBinWidth(iy);
    g->SetPoint(np, x, fr.sigma);
    g->SetPointError(np, ex, fr.sigmaErr);
    np++;
    delete px;
  }
  cout << "  " << name << ": " << np << " points, x=["
       << (np>0 ? Form("%.1f", g->GetX()[0]) : "N/A") << ", "
       << (np>0 ? Form("%.1f", g->GetX()[np-1]) : "N/A") << "]" << endl;
  return g;
}

// --- Merge two TGraphErrors (low + high pT) into one ---
TGraphErrors* MergeGraphs(TGraphErrors* gLow, TGraphErrors* gHigh, const TString& name) {
  TGraphErrors* g = new TGraphErrors();
  g->SetName(name);
  int np = 0;
  if (gLow) {
    for (int i = 0; i < gLow->GetN(); i++) {
      double x, y; gLow->GetPoint(i, x, y);
      g->SetPoint(np, x, y);
      g->SetPointError(np, gLow->GetErrorX(i), gLow->GetErrorY(i));
      np++;
    }
  }
  if (gHigh) {
    for (int i = 0; i < gHigh->GetN(); i++) {
      double x, y; gHigh->GetPoint(i, x, y);
      g->SetPoint(np, x, y);
      g->SetPointError(np, gHigh->GetErrorX(i), gHigh->GetErrorY(i));
      np++;
    }
  }
  return g;
}

// --- Make ratio graph: g / gRef ---
TGraphErrors* MakeRatioGraph(TGraphErrors* g, TGraphErrors* gRef, const TString& name) {
  if (!g || !gRef || gRef->GetN() == 0) return nullptr;
  TGraphErrors* gr = new TGraphErrors();
  gr->SetName(name);
  int np = 0;
  for (int i = 0; i < g->GetN(); i++) {
    double x, y; g->GetPoint(i, x, y);
    // Find closest reference point
    double bestDist = 1e9, yRef = 0, eyRef = 0;
    int bestJ = -1;
    for (int j = 0; j < gRef->GetN(); j++) {
      double xr, yr; gRef->GetPoint(j, xr, yr);
      if (std::abs(xr - x) < bestDist) {
        bestDist = std::abs(xr - x);
        yRef = yr;
        eyRef = gRef->GetErrorY(j);
        bestJ = j;
      }
    }
    if (bestJ < 0 || bestDist > 0.5 || yRef == 0) continue;
    double ratio = y / yRef;
    double ey = g->GetErrorY(i);
    double relErr = 0;
    if (y != 0) relErr += (ey / y) * (ey / y);
    if (yRef != 0) relErr += (eyRef / yRef) * (eyRef / yRef);
    double eRatio = ratio * std::sqrt(relErr);
    gr->SetPoint(np, x, ratio);
    gr->SetPointError(np, g->GetErrorX(i), eRatio);
    np++;
  }
  return gr;
}

// --- Ensure output directory exists ---
void EnsureDir(const TString& dir) {
  gSystem->mkdir(dir, true);
}

// --- Open file and get directory ---
TDirectory* OpenDir(const TString& filePath, const char* dirName) {
  TString fn = filePath;
  gSystem->ExpandPathName(fn);
  TFile* f = TFile::Open(fn, "READ");
  if (!f || f->IsZombie()) {
    cout << "Error: Cannot open " << fn << endl;
    return nullptr;
  }
  TDirectory* d = (TDirectory*)f->Get(dirName);
  if (!d) {
    cout << "Warning: No directory " << dirName << " in " << fn << endl;
  }
  return d;
}

// --- Style a graph ---
void StyleGraph(TGraphErrors* g, const Source& src) {
  if (!g) return;
  g->SetMarkerColor(src.color);
  g->SetLineColor(src.color);
  g->SetMarkerStyle(src.marker);
  g->SetMarkerSize(0.9);
  g->SetLineWidth(2);
}

// --- Generic graph comparison with Filipad2 ratio panel ---
static int sPadCount = 0;
void DrawGraphComparison(
    std::vector<TGraphErrors*>& graphs,
    const std::vector<Source>& sources, int refIdx,
    const TString& xtitle, const TString& ytitle,
    double xmin, double xmax, double ymin, double ymax,
    double ratioMin, double ratioMax,
    const TString& plotTitle, const TString& outName,
    bool logX = false, bool logY = false)
{
  if (graphs.empty()) return;
  EnsureDir(kOutputDir);

  sPadCount++;
  Filipad2* pad = new Filipad2(outName, sPadCount, 4.0, 0.4, 120, 80, 0.7 * 1.5, 1, 1);
  pad->Draw();

  TPad* top = pad->GetPad(1);
  top->SetLeftMargin(0.14);
  top->SetRightMargin(0.03);
  top->SetTopMargin(0.04);
  top->SetBottomMargin(0.001);
  top->SetGridx(1); top->SetGridy(1);
  if (logX) top->SetLogx(1);
  if (logY) top->SetLogy(1);
  top->cd();

  TH1F* hFrame = top->DrawFrame(xmin, ymin, xmax, ymax);
  hFrame->GetXaxis()->SetTitle("");
  hFrame->GetYaxis()->SetTitle(ytitle);
  hFrame->GetYaxis()->SetTitleOffset(1.0);
  hFrame->GetYaxis()->SetTitleSize(0.065);
  hFrame->GetYaxis()->SetLabelSize(0.055);

  TLatex* lat = new TLatex();
  lat->SetNDC(); lat->SetTextSize(0.055);
  lat->DrawLatex(0.17, 0.90, plotTitle);

  TLegend* leg = new TLegend(0.17, 0.62, 0.55, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);
  leg->SetTextSize(0.050);

  for (size_t i = 0; i < graphs.size(); i++) {
    if (!graphs[i]) continue;
    graphs[i]->Draw("PESAME");
    leg->AddEntry(graphs[i], sources[i].label, "pe");
  }
  leg->Draw();

  // Ratio panel
  TPad* bot = pad->GetPad(2);
  bot->SetLeftMargin(0.14);
  bot->SetRightMargin(0.03);
  bot->SetTopMargin(0.001);
  bot->SetBottomMargin(0.28);
  bot->SetGridx(1); bot->SetGridy(1);
  if (logX) bot->SetLogx(1);
  bot->cd();

  TH1F* hFrameR = bot->DrawFrame(xmin, ratioMin, xmax, ratioMax);
  hFrameR->GetXaxis()->SetTitle(xtitle);
  hFrameR->GetYaxis()->SetTitle(Form("Ratio to %s", sources[refIdx].label.Data()));
  hFrameR->GetXaxis()->SetTitleOffset(1.1);
  hFrameR->GetYaxis()->SetTitleOffset(0.55);
  hFrameR->GetXaxis()->SetTitleSize(0.11);
  hFrameR->GetYaxis()->SetTitleSize(0.09);
  hFrameR->GetXaxis()->SetLabelSize(0.10);
  hFrameR->GetYaxis()->SetLabelSize(0.08);
  hFrameR->GetYaxis()->SetNdivisions(505);

  TLine* line1 = new TLine(xmin, 1.0, xmax, 1.0);
  line1->SetLineStyle(2); line1->SetLineWidth(1);
  line1->Draw("SAME");

  TGraphErrors* gRef = (refIdx < (int)graphs.size()) ? graphs[refIdx] : nullptr;
  for (size_t i = 0; i < graphs.size(); i++) {
    if (!graphs[i] || (int)i == refIdx) continue;
    TGraphErrors* gRat = MakeRatioGraph(graphs[i], gRef, Form("ratio_%s_%zu", outName.Data(), i));
    if (gRat) {
      StyleGraph(gRat, sources[i]);
      gRat->Draw("PESAME");
    }
  }
  // Reference as flat line at 1
  if (gRef) {
    TGraphErrors* gRefRat = MakeRatioGraph(gRef, gRef, Form("ratio_ref_%s", outName.Data()));
    if (gRefRat) {
      StyleGraph(gRefRat, sources[refIdx]);
      gRefRat->SetMarkerSize(0.5);
      gRefRat->Draw("PESAME");
    }
  }

  pad->C->SaveAs(Form("%s/%s.pdf", kOutputDir.Data(), outName.Data()));
  cout << "Saved: " << kOutputDir << "/" << outName << ".pdf" << endl;
}

// --- Generic histogram comparison with Filipad2 ratio panel ---
void DrawHistComparison(
    std::vector<TH1*>& hists,
    const std::vector<Source>& sources, int refIdx,
    const TString& xtitle, const TString& ytitle,
    double xmin, double xmax, double ymin, double ymax,
    double ratioMin, double ratioMax,
    const TString& plotTitle, const TString& outName,
    bool logX = false, bool logY = true)
{
  if (hists.empty()) return;
  EnsureDir(kOutputDir);

  // Auto y-range if ymin==ymax==0
  if (ymin == 0 && ymax == 0) {
    double lo = 1e30, hi = -1e30;
    for (auto* h : hists) {
      if (!h) continue;
      for (int b = 1; b <= h->GetNbinsX(); b++) {
        double c = h->GetBinContent(b);
        if (c != 0) { lo = std::min(lo, c); hi = std::max(hi, c); }
      }
    }
    double margin = (hi - lo) * 0.15;
    ymin = lo - margin;
    ymax = hi + margin;
    if (logY) { ymin = lo * 0.3; ymax = hi * 5.0; }
  }

  sPadCount++;
  Filipad2* pad = new Filipad2(outName, sPadCount, 4.0, 0.4, 120, 80, 0.7 * 1.5, 1, 1);
  pad->Draw();

  TPad* top = pad->GetPad(1);
  top->SetLeftMargin(0.14);
  top->SetRightMargin(0.03);
  top->SetTopMargin(0.04);
  top->SetBottomMargin(0.001);
  top->SetGridx(1); top->SetGridy(1);
  if (logX) top->SetLogx(1);
  if (logY) top->SetLogy(1);
  top->cd();

  TH1F* hFrame = top->DrawFrame(xmin, ymin, xmax, ymax);
  hFrame->GetXaxis()->SetTitle("");
  hFrame->GetYaxis()->SetTitle(ytitle);
  hFrame->GetYaxis()->SetTitleOffset(1.0);
  hFrame->GetYaxis()->SetTitleSize(0.065);
  hFrame->GetYaxis()->SetLabelSize(0.055);

  TLatex* lat = new TLatex();
  lat->SetNDC(); lat->SetTextSize(0.055);
  lat->DrawLatex(0.17, 0.90, plotTitle);

  TLegend* leg = new TLegend(0.50, 0.65, 0.95, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);
  leg->SetTextSize(0.050);

  for (size_t i = 0; i < hists.size(); i++) {
    if (!hists[i]) continue;
    hists[i]->SetMarkerColor(sources[i].color);
    hists[i]->SetLineColor(sources[i].color);
    hists[i]->SetMarkerStyle(sources[i].marker);
    hists[i]->SetMarkerSize(0.7);
    hists[i]->Draw("PE SAME");
    leg->AddEntry(hists[i], sources[i].label, "pe");
  }
  leg->Draw();

  // Ratio panel
  TPad* bot = pad->GetPad(2);
  bot->SetLeftMargin(0.14);
  bot->SetRightMargin(0.03);
  bot->SetTopMargin(0.001);
  bot->SetBottomMargin(0.28);
  bot->SetGridx(1); bot->SetGridy(1);
  if (logX) bot->SetLogx(1);
  bot->cd();

  TH1F* hFrameR = bot->DrawFrame(xmin, ratioMin, xmax, ratioMax);
  hFrameR->GetXaxis()->SetTitle(xtitle);
  hFrameR->GetYaxis()->SetTitle(Form("Ratio to %s", sources[refIdx].label.Data()));
  hFrameR->GetXaxis()->SetTitleOffset(1.1);
  hFrameR->GetYaxis()->SetTitleOffset(0.55);
  hFrameR->GetXaxis()->SetTitleSize(0.11);
  hFrameR->GetYaxis()->SetTitleSize(0.09);
  hFrameR->GetXaxis()->SetLabelSize(0.10);
  hFrameR->GetYaxis()->SetLabelSize(0.08);
  hFrameR->GetYaxis()->SetNdivisions(505);

  TLine* line1 = new TLine(xmin, 1.0, xmax, 1.0);
  line1->SetLineStyle(2); line1->SetLineWidth(1);
  line1->Draw("SAME");

  TH1* hRef = (refIdx < (int)hists.size()) ? hists[refIdx] : nullptr;
  for (size_t i = 0; i < hists.size(); i++) {
    if (!hists[i]) continue;
    TH1* hRat = (TH1*)hists[i]->Clone(Form("hRat_%s_%zu", outName.Data(), i));
    hRat->SetDirectory(nullptr);
    if (hRef) hRat->Divide(hRef);
    hRat->SetMarkerColor(sources[i].color);
    hRat->SetLineColor(sources[i].color);
    hRat->SetMarkerStyle(sources[i].marker);
    hRat->SetMarkerSize(0.7);
    hRat->Draw("PE SAME");
  }

  pad->C->SaveAs(Form("%s/%s.pdf", kOutputDir.Data(), outName.Data()));
  cout << "Saved: " << kOutputDir << "/" << outName << ".pdf" << endl;
}

// ============================================================
//  PLOT FUNCTIONS
// ============================================================

// 1. sigma(pT)/pT vs pT from covariance matrix + MC residual overlay
void PlotSigmaPtOverPt() {
  if (!kDrawSigmaPtOverPt) return;
  cout << "\n=== sigma(pT)/pT ===" << endl;
  std::vector<Source> sources = kSources;
  if (kDrawData) sources.push_back(kDataSource);
  int refIdx = kDrawData ? (int)sources.size() - 1 : kRefIndex;

  // --- Extract covariance matrix sigma(pT)/pT ---
  std::vector<TGraphErrors*> graphsCov;
  for (size_t i = 0; i < sources.size(); i++) {
    cout << " Source " << i << ": " << sources[i].label << endl;
    TDirectory* d = OpenDir(sources[i].file, kDirTrackEff);
    if (!d) { graphsCov.push_back(nullptr); continue; }
    TH2* h2Lo = (TH2*)d->Get("h2_track_pt_track_sigmapt");
    TH2* h2Hi = (TH2*)d->Get("h2_track_pt_high_track_sigmapt");
    if (h2Lo) h2Lo->RebinX(5);
    cout << "  h2Lo: " << (h2Lo ? Form("entries=%.0f, xRange=[%.1f,%.1f], nbinsX=%d", h2Lo->GetEntries(), h2Lo->GetXaxis()->GetXmin(), h2Lo->GetXaxis()->GetXmax(), h2Lo->GetNbinsX()) : "NULL") << endl;
    cout << "  h2Hi: " << (h2Hi ? Form("entries=%.0f, xRange=[%.1f,%.1f], nbinsX=%d", h2Hi->GetEntries(), h2Hi->GetXaxis()->GetXmin(), h2Hi->GetXaxis()->GetXmax(), h2Hi->GetNbinsX()) : "NULL") << endl;
    TGraphErrors* gLo = ExtractMeanProfile(h2Lo, Form("sigmaPt_lo_%zu", i));
    TGraphErrors* gHi = ExtractMeanProfileVarBins(h2Hi, Form("sigmaPt_hi_%zu", i), {{40, 2}, {50, 5}, {70, 10}});
    TGraphErrors* gMerge = MergeGraphs(gLo, gHi, Form("sigmaPt_%zu", i));
    StyleGraph(gMerge, sources[i]);
    graphsCov.push_back(gMerge);
    delete gLo; delete gHi;
  }

  // --- Extract MC residual sigma(pT)/pT (MC only, no Data) ---
  // From h2_particle_pt_track_pt_deltaptoverparticlept: x=pT_truth, y=(pT_reco-pT_truth)/pT_truth
  std::vector<TGraphErrors*> graphsRes;
  for (size_t i = 0; i < kSources.size(); i++) {
    cout << " MC residual source " << i << ": " << kSources[i].label << endl;
    TDirectory* d = OpenDir(kSources[i].file, kDirTrackEff);
    if (!d) { graphsRes.push_back(nullptr); continue; }
    TH2* h2Res = (TH2*)d->Get("h2_particle_pt_track_pt_deltaptoverparticlept");
    if (!h2Res) { cout << "  deltaptoverparticlept: NOT FOUND" << endl; graphsRes.push_back(nullptr); continue; }
    cout << "  deltaptoverparticlept: entries=" << h2Res->GetEntries()
         << ", x=[" << h2Res->GetXaxis()->GetXmin() << "," << h2Res->GetXaxis()->GetXmax()
         << "], nbinsX=" << h2Res->GetNbinsX() << endl;
    // Use same variable binning as high-pT covariance histogram for consistency
    TGraphErrors* gRes = ExtractSigmaVsPt(h2Res, Form("sigmaPtRes_%zu", i), 50);
    if (gRes) {
      gRes->SetMarkerColor(kSources[i].color);
      gRes->SetLineColor(kSources[i].color);
      gRes->SetMarkerStyle(kSources[i].marker + 4); // open marker variant
      gRes->SetMarkerSize(0.7);
      gRes->SetLineWidth(2);
      gRes->SetLineStyle(2); // dashed
    }
    graphsRes.push_back(gRes);
  }

  // --- Draw combined plot ---
  EnsureDir(kOutputDir);
  sPadCount++;
  TString outName = "sigmaPtOverPt";
  Filipad2* pad = new Filipad2(outName, sPadCount, 4.0, 0.4, 120, 80, 0.7 * 1.5, 1, 1);
  pad->Draw();

  TPad* top = pad->GetPad(1);
  top->SetLeftMargin(0.14);
  top->SetRightMargin(0.03);
  top->SetTopMargin(0.04);
  top->SetBottomMargin(0.001);
  top->SetGridx(1); top->SetGridy(1);
  top->cd();

  TH1F* hFrame = top->DrawFrame(0, 0, 100, 0.4);
  hFrame->GetXaxis()->SetTitle("");
  hFrame->GetYaxis()->SetTitle("#sigma(#it{p}_{T})/#it{p}_{T}");
  hFrame->GetYaxis()->SetTitleOffset(1.0);
  hFrame->GetYaxis()->SetTitleSize(0.065);
  hFrame->GetYaxis()->SetLabelSize(0.055);

  TLatex* lat = new TLatex();
  lat->SetNDC(); lat->SetTextSize(0.050);
  lat->DrawLatex(0.17, 0.90, "#sigma(#it{p}_{T})/#it{p}_{T}");

  TLegend* leg = new TLegend(0.17, 0.45, 0.60, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);
  leg->SetTextSize(0.042);

  // Draw covariance matrix graphs (solid)
  for (size_t i = 0; i < graphsCov.size(); i++) {
    if (!graphsCov[i]) continue;
    graphsCov[i]->Draw("PESAME");
    leg->AddEntry(graphsCov[i], Form("%s (cov. mat.)", sources[i].label.Data()), "pe");
  }
  // Draw MC residual graphs (dashed, open markers)
  for (size_t i = 0; i < graphsRes.size(); i++) {
    if (!graphsRes[i]) continue;
    graphsRes[i]->Draw("PESAME");
    leg->AddEntry(graphsRes[i], Form("%s (MC residual)", kSources[i].label.Data()), "pe");
  }
  leg->Draw();

  // --- Ratio panel: all to Data (cov matrix) ---
  TPad* bot = pad->GetPad(2);
  bot->SetLeftMargin(0.14);
  bot->SetRightMargin(0.03);
  bot->SetTopMargin(0.001);
  bot->SetBottomMargin(0.28);
  bot->SetGridx(1); bot->SetGridy(1);
  bot->cd();

  TH1F* hFrameR = bot->DrawFrame(0, 0.0, 100, 2.0);
  hFrameR->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hFrameR->GetYaxis()->SetTitle(Form("Ratio to %s", sources[refIdx].label.Data()));
  hFrameR->GetXaxis()->SetTitleOffset(1.1);
  hFrameR->GetYaxis()->SetTitleOffset(0.55);
  hFrameR->GetXaxis()->SetTitleSize(0.11);
  hFrameR->GetYaxis()->SetTitleSize(0.09);
  hFrameR->GetXaxis()->SetLabelSize(0.10);
  hFrameR->GetYaxis()->SetLabelSize(0.08);
  hFrameR->GetYaxis()->SetNdivisions(505);

  TLine* line1 = new TLine(0, 1.0, 100, 1.0);
  line1->SetLineStyle(2); line1->SetLineWidth(1);
  line1->Draw("SAME");

  TGraphErrors* gRef = (refIdx < (int)graphsCov.size()) ? graphsCov[refIdx] : nullptr;
  // Covariance matrix ratios
  for (size_t i = 0; i < graphsCov.size(); i++) {
    if (!graphsCov[i] || (int)i == refIdx) continue;
    TGraphErrors* gRat = MakeRatioGraph(graphsCov[i], gRef, Form("ratio_cov_%zu", i));
    if (gRat) { StyleGraph(gRat, sources[i]); gRat->Draw("PESAME"); }
  }
  // MC residual ratios (to same Data cov matrix reference)
  for (size_t i = 0; i < graphsRes.size(); i++) {
    if (!graphsRes[i] || !gRef) continue;
    TGraphErrors* gRat = MakeRatioGraph(graphsRes[i], gRef, Form("ratio_res_%zu", i));
    if (gRat) {
      gRat->SetMarkerColor(kSources[i].color);
      gRat->SetLineColor(kSources[i].color);
      gRat->SetMarkerStyle(kSources[i].marker + 4);
      gRat->SetMarkerSize(0.7);
      gRat->SetLineStyle(2);
      gRat->Draw("PESAME");
    }
  }
  // Reference flat line
  if (gRef) {
    TGraphErrors* gRefRat = MakeRatioGraph(gRef, gRef, "ratio_ref_cov");
    if (gRefRat) { StyleGraph(gRefRat, sources[refIdx]); gRefRat->SetMarkerSize(0.5); gRefRat->Draw("PESAME"); }
  }

  pad->C->SaveAs(Form("%s/%s.pdf", kOutputDir.Data(), outName.Data()));
  cout << "Saved: " << kOutputDir << "/" << outName << ".pdf" << endl;
}

// 2. sigma(1/pT) vs pT from covariance matrix
void PlotSigma1OverPt() {
  if (!kDrawSigma1OverPt) return;
  cout << "\n=== sigma(1/pT) ===" << endl;
  std::vector<Source> sources = kSources;
  if (kDrawData) sources.push_back(kDataSource);
  int refIdx = kDrawData ? (int)sources.size() - 1 : kRefIndex;
  std::vector<TGraphErrors*> graphs;
  for (size_t i = 0; i < sources.size(); i++) {
    cout << " Source " << i << ": " << sources[i].label << endl;
    TDirectory* d = OpenDir(sources[i].file, kDirTrackEff);
    if (!d) { graphs.push_back(nullptr); continue; }
    TH2* h2Lo = (TH2*)d->Get("h2_track_pt_track_sigma1overpt");
    TH2* h2Hi = (TH2*)d->Get("h2_track_pt_high_track_sigma1overpt");
    if (h2Lo) h2Lo->RebinX(5);
    cout << "  h2Lo: " << (h2Lo ? Form("entries=%.0f, xRange=[%.1f,%.1f], nbinsX=%d", h2Lo->GetEntries(), h2Lo->GetXaxis()->GetXmin(), h2Lo->GetXaxis()->GetXmax(), h2Lo->GetNbinsX()) : "NULL") << endl;
    cout << "  h2Hi: " << (h2Hi ? Form("entries=%.0f, xRange=[%.1f,%.1f], nbinsX=%d", h2Hi->GetEntries(), h2Hi->GetXaxis()->GetXmin(), h2Hi->GetXaxis()->GetXmax(), h2Hi->GetNbinsX()) : "NULL") << endl;
    TGraphErrors* gLo = ExtractMeanProfile(h2Lo, Form("sigma1pt_lo_%zu", i));
    TGraphErrors* gHi = ExtractMeanProfileVarBins(h2Hi, Form("sigma1pt_hi_%zu", i), {{40, 2}, {50, 5}, {70, 10}});
    TGraphErrors* gMerge = MergeGraphs(gLo, gHi, Form("sigma1pt_%zu", i));
    StyleGraph(gMerge, sources[i]);
    graphs.push_back(gMerge);
    delete gLo; delete gHi;
  }
  DrawGraphComparison(graphs, sources, refIdx,
    "#it{p}_{T} (GeV/#it{c})", "#sigma(1/#it{p}_{T}) ((GeV/#it{c})^{-1})",
    0, 100, 0, 0.02, 0.0, 2.0,
    "#sigma(1/#it{p}_{T}) (cov. matrix)", "sigma1OverPt");
}

// 3. pT residual sigma vs pT (Gaussian fit to (pT_reco - pT_true)/pT_true)
void PlotPtResidualSigma() {
  if (!kDrawPtResidualSigma) return;
  cout << "\n=== pT residual sigma ===" << endl;
  std::vector<TGraphErrors*> graphs;
  for (size_t i = 0; i < kSources.size(); i++) {
    cout << " Source " << i << ": " << kSources[i].label << endl;
    TDirectory* d = OpenDir(kSources[i].file, kDirTrackEff);
    if (!d) { graphs.push_back(nullptr); continue; }
    TH2* h2Lo = (TH2*)d->Get("h2_particle_pt_track_pt_residual_associatedtrack_primary");
    TH2* h2Hi = (TH2*)d->Get("h2_particle_pt_high_track_pt_high_residual_associatedtrack_primary");
    if (h2Lo) h2Lo->RebinX(10);
    cout << "  h2Lo: " << (h2Lo ? Form("entries=%.0f, xRange=[%.1f,%.1f], nbinsX=%d", h2Lo->GetEntries(), h2Lo->GetXaxis()->GetXmin(), h2Lo->GetXaxis()->GetXmax(), h2Lo->GetNbinsX()) : "NULL") << endl;
    cout << "  h2Hi: " << (h2Hi ? Form("entries=%.0f, xRange=[%.1f,%.1f], nbinsX=%d", h2Hi->GetEntries(), h2Hi->GetXaxis()->GetXmin(), h2Hi->GetXaxis()->GetXmax(), h2Hi->GetNbinsX()) : "NULL") << endl;
    TGraphErrors* gLo = ExtractSigmaVsPt(h2Lo, Form("ptres_lo_%zu", i));
    TGraphErrors* gHi = ExtractSigmaVsPt(h2Hi, Form("ptres_hi_%zu", i));
    TGraphErrors* gMerge = MergeGraphs(gLo, gHi, Form("ptres_%zu", i));
    StyleGraph(gMerge, kSources[i]);
    graphs.push_back(gMerge);
    delete gLo; delete gHi;
  }
  DrawGraphComparison(graphs, kSources, kRefIndex,
    "#it{p}_{T}^{true} (GeV/#it{c})", "#sigma((#it{p}_{T}^{reco}-#it{p}_{T}^{true})/#it{p}_{T}^{true})",
    0, 100, 0, 0.25, 0.85, 1.15,
    "#it{p}_{T} residual #sigma (Gaussian fit)", "ptResidualSigma");
}

// 4. DCA_xy width vs pT (DCA histo: x=DCA, y=pT)
void PlotDCAxyWidth() {
  if (!kDrawDCAxyWidth) return;
  cout << "\n=== DCA_xy width ===" << endl;
  std::vector<Source> sources = kSources;
  if (kDrawData) sources.push_back(kDataSource);
  std::vector<TGraphErrors*> graphs;
  for (size_t i = 0; i < sources.size(); i++) {
    cout << " Source " << i << ": " << sources[i].label << endl;
    TDirectory* d = OpenDir(sources[i].file, kDirPropagation);
    if (!d) { graphs.push_back(nullptr); continue; }
    TH2* h2 = (TH2*)d->Get("hDCAxyVsPtRec");
    if (h2) cout << Form("  hDCAxyVsPtRec: entries=%.0f, x=[%.3f,%.3f] nbX=%d, y=[%.1f,%.1f] nbY=%d",
      h2->GetEntries(), h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax(), h2->GetNbinsX(),
      h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax(), h2->GetNbinsY()) << endl;
    TGraphErrors* g = ExtractDCASigmaVsPt(h2, Form("dcaxy_w_%zu", i));
    StyleGraph(g, sources[i]);
    graphs.push_back(g);
  }
  DrawGraphComparison(graphs, sources, kRefIndex,
    "#it{p}_{T} (GeV/#it{c})", "#sigma(DCA_{xy}) (cm)",
    0, 20, 0, 0.06, 0.90, 1.10,
    "DCA_{xy} width (Gaussian fit)", "DCAxy_width");
}

// 5. DCA_z width vs pT
void PlotDCAzWidth() {
  if (!kDrawDCAzWidth) return;
  cout << "\n=== DCA_z width ===" << endl;
  std::vector<Source> sources = kSources;
  if (kDrawData) sources.push_back(kDataSource);
  std::vector<TGraphErrors*> graphs;
  for (size_t i = 0; i < sources.size(); i++) {
    cout << " Source " << i << ": " << sources[i].label << endl;
    TDirectory* d = OpenDir(sources[i].file, kDirPropagation);
    if (!d) { graphs.push_back(nullptr); continue; }
    TH2* h2 = (TH2*)d->Get("hDCAzVsPtRec");
    if (h2) cout << Form("  hDCAzVsPtRec: entries=%.0f, x=[%.3f,%.3f] nbX=%d, y=[%.1f,%.1f] nbY=%d",
      h2->GetEntries(), h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax(), h2->GetNbinsX(),
      h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax(), h2->GetNbinsY()) << endl;
    TGraphErrors* g = ExtractDCASigmaVsPt(h2, Form("dcaz_w_%zu", i));
    StyleGraph(g, sources[i]);
    graphs.push_back(g);
  }
  DrawGraphComparison(graphs, sources, kRefIndex,
    "#it{p}_{T} (GeV/#it{c})", "#sigma(DCA_{z}) (cm)",
    0, 20, 0, 0.06, 0.90, 1.10,
    "DCA_{z} width (Gaussian fit)", "DCAz_width");
}

// 6. DCA_xy mean (bias) vs pT (DCA histo: x=DCA, y=pT)
void PlotDCAxyMean() {
  if (!kDrawDCAxyMean) return;
  cout << "\n=== DCA_xy mean ===" << endl;
  std::vector<Source> sources = kSources;
  if (kDrawData) sources.push_back(kDataSource);
  std::vector<TGraphErrors*> graphs;
  for (size_t i = 0; i < sources.size(); i++) {
    TDirectory* d = OpenDir(sources[i].file, kDirPropagation);
    if (!d) { graphs.push_back(nullptr); continue; }
    TH2* h2 = (TH2*)d->Get("hDCAxyVsPtRec");
    TGraphErrors* g = ExtractDCAMeanVsPt(h2, Form("dcaxy_m_%zu", i));
    StyleGraph(g, sources[i]);
    graphs.push_back(g);
  }
  // For mean: no ratio, just overlay
  EnsureDir(kOutputDir);
  TCanvas* c = new TCanvas("cDCAxyMean", "", 900, 600);
  c->SetLeftMargin(0.14); c->SetRightMargin(0.03);
  c->SetGridx(1); c->SetGridy(1);
  TH1F* hf = c->DrawFrame(0, -0.005, 20, 0.005);
  hf->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hf->GetYaxis()->SetTitle("#LT DCA_{xy} #GT (cm)");
  hf->GetYaxis()->SetTitleOffset(1.2);
  TLegend* leg = new TLegend(0.55, 0.70, 0.95, 0.90);
  leg->SetBorderSize(0); leg->SetFillColorAlpha(0,0); leg->SetTextSize(0.040);
  for (size_t i = 0; i < graphs.size(); i++) {
    if (!graphs[i]) continue;
    graphs[i]->Draw("PESAME");
    leg->AddEntry(graphs[i], sources[i].label, "pe");
  }
  TLine* l0 = new TLine(0, 0, 20, 0);
  l0->SetLineStyle(2); l0->Draw("SAME");
  TLatex* lat = new TLatex(); lat->SetNDC(); lat->SetTextSize(0.045);
  lat->DrawLatex(0.17, 0.88, "DCA_{xy} mean (bias)");
  leg->Draw();
  c->SaveAs(Form("%s/DCAxy_mean.pdf", kOutputDir.Data()));
  cout << "Saved: " << kOutputDir << "/DCAxy_mean.pdf" << endl;
}

// 7. DCA_z mean (bias) vs pT
void PlotDCAzMean() {
  if (!kDrawDCAzMean) return;
  cout << "\n=== DCA_z mean ===" << endl;
  std::vector<Source> sources = kSources;
  if (kDrawData) sources.push_back(kDataSource);
  std::vector<TGraphErrors*> graphs;
  for (size_t i = 0; i < sources.size(); i++) {
    TDirectory* d = OpenDir(sources[i].file, kDirPropagation);
    if (!d) { graphs.push_back(nullptr); continue; }
    TH2* h2 = (TH2*)d->Get("hDCAzVsPtRec");
    TGraphErrors* g = ExtractDCAMeanVsPt(h2, Form("dcaz_m_%zu", i));
    StyleGraph(g, sources[i]);
    graphs.push_back(g);
  }
  EnsureDir(kOutputDir);
  TCanvas* c = new TCanvas("cDCAzMean", "", 900, 600);
  c->SetLeftMargin(0.14); c->SetRightMargin(0.03);
  c->SetGridx(1); c->SetGridy(1);
  TH1F* hf = c->DrawFrame(0, -0.005, 20, 0.005);
  hf->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hf->GetYaxis()->SetTitle("#LT DCA_{z} #GT (cm)");
  hf->GetYaxis()->SetTitleOffset(1.2);
  TLegend* leg = new TLegend(0.55, 0.70, 0.95, 0.90);
  leg->SetBorderSize(0); leg->SetFillColorAlpha(0,0); leg->SetTextSize(0.040);
  for (size_t i = 0; i < graphs.size(); i++) {
    if (!graphs[i]) continue;
    graphs[i]->Draw("PESAME");
    leg->AddEntry(graphs[i], sources[i].label, "pe");
  }
  TLine* l0 = new TLine(0, 0, 20, 0);
  l0->SetLineStyle(2); l0->Draw("SAME");
  TLatex* lat = new TLatex(); lat->SetNDC(); lat->SetTextSize(0.045);
  lat->DrawLatex(0.17, 0.88, "DCA_{z} mean (bias)");
  leg->Draw();
  c->SaveAs(Form("%s/DCAz_mean.pdf", kOutputDir.Data()));
  cout << "Saved: " << kOutputDir << "/DCAz_mean.pdf" << endl;
}

// 8. Jet pT detector level
void PlotJetPtDet() {
  if (!kDrawJetPtDet) return;
  cout << "\n=== Jet pT (detector level) ===" << endl;
  std::vector<TH1*> hists;
  for (size_t i = 0; i < kSources.size(); i++) {
    TDirectory* d = OpenDir(kSources[i].file, kDirJetSpectra);
    if (!d) { hists.push_back(nullptr); continue; }
    TH1* h = (TH1*)d->Get("h_jet_pt");
    if (h) {
      h = h->Rebin(kNJetPtBinsReco, Form("hJetPtDet_%zu", i), kJetPtBinsReco);
      h->SetDirectory(nullptr);
      if (h->Integral() > 0) h->Scale(1.0 / h->Integral(), "width");
    }
    hists.push_back(h);
  }
  DrawHistComparison(hists, kSources, kRefIndex,
    "#it{p}_{T,jet}^{det} (GeV/#it{c})", "1/N dN/d#it{p}_{T}",
    0, 150, 1e-7, 1, 0.85, 1.15,
    "Jet #it{p}_{T} (detector level)", "jetPt_det");
}

// 9. Jet pT UE-subtracted
void PlotJetPtUESub() {
  if (!kDrawJetPtUESub) return;
  cout << "\n=== Jet pT (UE-subtracted) ===" << endl;
  std::vector<TH1*> hists;
  for (size_t i = 0; i < kSources.size(); i++) {
    TDirectory* d = OpenDir(kSources[i].file, kDirJetSpectra);
    if (!d) { hists.push_back(nullptr); continue; }
    TH1* h = (TH1*)d->Get("h_jet_pt_rhoareasubtracted");
    if (h) {
      h = h->Rebin(kNJetPtBinsRecoUE, Form("hJetPtUE_%zu", i), kJetPtBinsRecoUE);
      h->SetDirectory(nullptr);
      if (h->Integral() > 0) h->Scale(1.0 / h->Integral(), "width");
    }
    hists.push_back(h);
  }
  DrawHistComparison(hists, kSources, kRefIndex,
    "#it{p}_{T,jet}^{det} - #rho A (GeV/#it{c})", "1/N dN/d#it{p}_{T}",
    -20, 150, 1e-7, 1, 0.85, 1.15,
    "Jet #it{p}_{T} (UE-subtracted)", "jetPt_UEsub");
}

// 10. Jet pT particle level (sanity check - should be identical)
void PlotJetPtParticle() {
  if (!kDrawJetPtParticle) return;
  cout << "\n=== Jet pT (particle level) ===" << endl;
  std::vector<TH1*> hists;
  for (size_t i = 0; i < kSources.size(); i++) {
    TDirectory* d = OpenDir(kSources[i].file, kDirJetSpectra);
    if (!d) { hists.push_back(nullptr); continue; }
    TH1* h = (TH1*)d->Get("h_jet_pt_part");
    if (h) {
      h = h->Rebin(kNJetPtBinsTruth, Form("hJetPtPart_%zu", i), kJetPtBinsTruth);
      h->SetDirectory(nullptr);
      if (h->Integral() > 0) h->Scale(1.0 / h->Integral(), "width");
    }
    hists.push_back(h);
  }
  DrawHistComparison(hists, kSources, kRefIndex,
    "#it{p}_{T,jet}^{part} (GeV/#it{c})", "1/N dN/d#it{p}_{T}",
    0, 200, 1e-8, 1, 0.95, 1.05,
    "Jet #it{p}_{T} (particle level, sanity)", "jetPt_particle");
}

// 11. delta pT sigma from residual (pT_reco - pT_true) absolute
void PlotDeltaPtSigma() {
  if (!kDrawDeltaPtSigma) return;
  cout << "\n=== delta pT sigma ===" << endl;
  std::vector<TGraphErrors*> graphs;
  for (size_t i = 0; i < kSources.size(); i++) {
    TDirectory* d = OpenDir(kSources[i].file, kDirTrackEff);
    if (!d) { graphs.push_back(nullptr); continue; }
    TH2* h2Lo = (TH2*)d->Get("h2_particle_pt_track_pt_deltapt");
    // No high-pT version for deltapt; use low only
    TGraphErrors* g = ExtractSigmaVsPt(h2Lo, Form("dpt_%zu", i));
    StyleGraph(g, kSources[i]);
    graphs.push_back(g);
  }
  DrawGraphComparison(graphs, kSources, kRefIndex,
    "#it{p}_{T}^{true} (GeV/#it{c})", "#sigma(#it{p}_{T}^{reco} - #it{p}_{T}^{true}) (GeV/#it{c})",
    0, 10, 0, 0.5, 0.85, 1.15,
    "#Delta#it{p}_{T} = #it{p}_{T}^{reco} - #it{p}_{T}^{true} (absolute)", "deltaPt_sigma");
}

// 12. Track pT residual distributions per pT class (like DrawJetResolution)
void PlotPtResidualDistributions() {
  if (!kDrawPtResidualDist) return;
  cout << "\n=== pT residual distributions ===" << endl;

  // pT classes: {ptLow, ptHigh, useHighPtHisto}
  struct PtClass { double lo, hi; bool useHigh; Color_t color; };
  std::vector<PtClass> ptClasses = {
    {1,   2,   false, kBlack},
    {2,   5,   false, kRed+1},
    {5,   10,  false, kBlue+1},
    {10,  20,  true,  kGreen+2},
    {20,  50,  true,  kMagenta+1},
  };

  for (size_t i = 0; i < kSources.size(); i++) {
    cout << " Source " << i << ": " << kSources[i].label << endl;
    TDirectory* d = OpenDir(kSources[i].file, kDirTrackEff);
    if (!d) continue;

    TH2* h2Lo = (TH2*)d->Get("h2_particle_pt_track_pt_residual_associatedtrack_primary");
    TH2* h2Hi = (TH2*)d->Get("h2_particle_pt_high_track_pt_high_residual_associatedtrack_primary");
    if (!h2Lo && !h2Hi) continue;

    EnsureDir(kOutputDir);
    sPadCount++;
    TCanvas* can = new TCanvas(Form("cPtResDist_%zu", i), "", 700, 500);
    can->SetLeftMargin(0.15); can->SetRightMargin(0.05);
    can->SetTopMargin(0.03); can->SetBottomMargin(0.15);
    can->SetLogy(1);

    TH1F* hFrame = can->DrawFrame(-0.5, 1e-3, 0.5, 50);
    hFrame->GetXaxis()->SetTitle("(#it{p}_{T}^{true} #minus #it{p}_{T}^{reco}) / #it{p}_{T}^{true}");
    hFrame->GetYaxis()->SetTitle("Probability density");
    hFrame->GetXaxis()->SetTitleSize(0.05);
    hFrame->GetYaxis()->SetTitleSize(0.05);
    hFrame->GetXaxis()->SetLabelSize(0.04);
    hFrame->GetYaxis()->SetLabelSize(0.04);
    hFrame->GetYaxis()->SetTitleOffset(1.3);

    TLegend* leg = new TLegend(0.52, 0.50, 0.93, 0.93);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry((TObject*)0, kSources[i].label, "");

    for (size_t j = 0; j < ptClasses.size(); j++) {
      TH2* h2 = ptClasses[j].useHigh ? h2Hi : h2Lo;
      if (!h2) continue;

      int binLo = h2->GetXaxis()->FindBin(ptClasses[j].lo);
      int binHi = h2->GetXaxis()->FindBin(ptClasses[j].hi - 1e-6);

      TH1D* hProj = h2->ProjectionY(
          Form("ptres_dist_%zu_%zu", i, j), binLo, binHi);
      hProj->SetDirectory(0);

      if (hProj->Integral() > 0) hProj->Scale(1.0 / hProj->Integral(), "width");

      hProj->SetLineColor(ptClasses[j].color);
      hProj->SetLineWidth(2);
      hProj->Draw("hist e same");

      // Gaussian fit and overlay
      FitResult fr = FitCoreGaussian(hProj);
      if (fr.ok) {
        TF1* fGaus = new TF1(Form("fGaus_%zu_%zu", i, j), "gaus",
                             fr.mean - 3.0 * fr.sigma, fr.mean + 3.0 * fr.sigma);
        fGaus->SetParameters(hProj->GetBinContent(hProj->FindBin(fr.mean)), fr.mean, fr.sigma);
        fGaus->SetLineColor(ptClasses[j].color);
        fGaus->SetLineStyle(2); // dashed
        fGaus->SetLineWidth(2);
        fGaus->Draw("same");
        leg->AddEntry(hProj,
          Form("[%.0f, %.0f] GeV/#it{c}, #sigma_{fit}=%.4f",
               ptClasses[j].lo, ptClasses[j].hi, fr.sigma), "l");
        leg->AddEntry(fGaus, Form("  Gauss fit (#mu=%.4f)", fr.mean), "l");
      } else {
        leg->AddEntry(hProj,
          Form("[%.0f, %.0f] GeV/#it{c}, RMS=%.4f (fit failed)",
               ptClasses[j].lo, ptClasses[j].hi, hProj->GetRMS()), "l");
      }
    }

    leg->Draw();
    TString tag = kSources[i].label;
    tag.ReplaceAll(" ", "_").ReplaceAll(",", "").ReplaceAll("#times", "x").ReplaceAll("#", "");
    can->SaveAs(Form("%s/ptResidualDist_%s.pdf", kOutputDir.Data(), tag.Data()));
    cout << "Saved: " << kOutputDir << "/ptResidualDist_" << tag << ".pdf" << endl;
  }
}

// Adaptive rebinning: merge bins from left to right until relative stat error < maxRelErr.
// Returns a new histogram with variable bin widths, divided by bin width (density).
// The reference histogram (typically Data with lowest stats at high pT) determines the bin edges.
std::vector<double> BuildAdaptiveBinEdges(TH1* hRef, double maxRelErr = 0.03, double ptMax = 200.0) {
  std::vector<double> edges;
  edges.push_back(hRef->GetXaxis()->GetBinLowEdge(1));
  double sumW = 0, sumW2 = 0;
  for (int i = 1; i <= hRef->GetNbinsX(); i++) {
    double lo = hRef->GetXaxis()->GetBinLowEdge(i);
    if (lo >= ptMax) break;
    double c = hRef->GetBinContent(i);
    double e = hRef->GetBinError(i);
    sumW += c;
    sumW2 += e * e;
    bool lastBin = (i == hRef->GetNbinsX() || hRef->GetXaxis()->GetBinUpEdge(i) >= ptMax);
    double relErr = (sumW > 0) ? sqrt(sumW2) / sumW : 1.0;
    if (relErr < maxRelErr || lastBin) {
      edges.push_back(std::min(hRef->GetXaxis()->GetBinUpEdge(i), ptMax));
      sumW = 0; sumW2 = 0;
    }
  }
  return edges;
}

TH1* RebinToEdges(TH1* hOrig, const std::vector<double>& edges, const char* name) {
  int nBins = edges.size() - 1;
  TH1D* hNew = new TH1D(name, hOrig->GetTitle(), nBins, edges.data());
  hNew->SetDirectory(nullptr);
  hNew->Sumw2();
  for (int i = 1; i <= hOrig->GetNbinsX(); i++) {
    double x = hOrig->GetXaxis()->GetBinCenter(i);
    double c = hOrig->GetBinContent(i);
    double e = hOrig->GetBinError(i);
    int newBin = hNew->FindBin(x);
    if (newBin >= 1 && newBin <= nBins) {
      hNew->SetBinContent(newBin, hNew->GetBinContent(newBin) + c);
      hNew->SetBinError(newBin, sqrt(hNew->GetBinError(newBin) * hNew->GetBinError(newBin) + e * e));
    }
  }
  // Divide by bin width for density
  for (int i = 1; i <= nBins; i++) {
    double w = hNew->GetBinWidth(i);
    if (w > 0) {
      hNew->SetBinContent(i, hNew->GetBinContent(i) / w);
      hNew->SetBinError(i, hNew->GetBinError(i) / w);
    }
  }
  return hNew;
}

// 13. Track pT distribution (event-normalized)
void PlotTrackPt() {
  if (!kDrawTrackPt) return;
  cout << "\n=== Track pT distribution ===" << endl;
  std::vector<Source> sources = kSources;
  if (kDrawData) sources.push_back(kDataSource);
  int refIdx = kDrawData ? (int)sources.size() - 1 : kRefIndex;

  // First pass: load raw histograms and event-normalize (no bin width division yet)
  std::vector<TH1*> rawHists;
  for (size_t i = 0; i < sources.size(); i++) {
    cout << " Source " << i << ": " << sources[i].label << endl;
    TDirectory* d = OpenDir(sources[i].file, kDirTrackEff);
    TH1* hPt = nullptr;
    if (d) {
      TH2* h2 = (TH2*)d->Get("h2_centrality_track_pt");
      if (h2) {
        hPt = h2->ProjectionY(Form("hTrackPtRaw_%zu", i));
        hPt->SetDirectory(nullptr);
        cout << "  h2_centrality_track_pt → projected, entries=" << hPt->GetEntries() << endl;
      }
      if (!hPt) {
        TH3* h3 = (TH3*)d->Get("h3_track_pt_track_eta_track_phi_associatedtrack_primary");
        if (h3) {
          hPt = h3->ProjectionX(Form("hTrackPtRaw_%zu", i));
          hPt->SetDirectory(nullptr);
          cout << "  h3_associatedtrack_primary → projected X, entries=" << hPt->GetEntries() << endl;
        }
      }
      if (!hPt) { cout << "  No track pT histogram found" << endl; rawHists.push_back(nullptr); continue; }

      double nEvents = 1.0;
      TH1* hColl = (TH1*)d->Get("h_collisions");
      if (hColl) {
        for (int b = hColl->GetNbinsX(); b >= 1; b--) {
          if (hColl->GetBinContent(b) > 0) { nEvents = hColl->GetBinContent(b); break; }
        }
        cout << "  Nevents (h_collisions) = " << nEvents << endl;
      } else {
        TH1* hMcColl = (TH1*)d->Get("hMcCollCutsCounts");
        if (hMcColl) {
          for (int b = hMcColl->GetNbinsX(); b >= 1; b--) {
            if (hMcColl->GetBinContent(b) > 0) { nEvents = hMcColl->GetBinContent(b); break; }
          }
          cout << "  Nevents (hMcCollCutsCounts) = " << nEvents << endl;
        } else {
          cout << "  WARNING: No event counter found, using raw counts" << endl;
        }
      }
      if (nEvents > 0) hPt->Scale(1.0 / nEvents);
    }
    rawHists.push_back(hPt);
  }

  // Build adaptive bin edges from the reference (Data or first MC with lowest stats)
  TH1* hRef = (refIdx < (int)rawHists.size()) ? rawHists[refIdx] : nullptr;
  if (!hRef) { cout << "ERROR: reference histogram not found" << endl; return; }
  std::vector<double> edges = BuildAdaptiveBinEdges(hRef, 0.015, 200.0);
  cout << "  Adaptive binning: " << edges.size() - 1 << " bins (from " << edges.front()
       << " to " << edges.back() << " GeV)" << endl;

  // Rebin all histograms to the adaptive bin edges (includes /binWidth for density)
  std::vector<TH1*> hists;
  for (size_t i = 0; i < rawHists.size(); i++) {
    if (!rawHists[i]) { hists.push_back(nullptr); continue; }
    hists.push_back(RebinToEdges(rawHists[i], edges, Form("hTrackPt_%zu", i)));
  }

  // Find y-range from data
  double yminAuto = 1e30, ymaxAuto = 0;
  for (size_t i = 0; i < hists.size(); i++) {
    if (!hists[i]) continue;
    for (int b = 1; b <= hists[i]->GetNbinsX(); b++) {
      double c = hists[i]->GetBinContent(b);
      if (c > 0) {
        yminAuto = std::min(yminAuto, c);
        ymaxAuto = std::max(ymaxAuto, c);
      }
    }
  }
  double ylo = yminAuto * 0.3;
  double yhi = ymaxAuto * 5.0;

  DrawHistComparison(hists, sources, refIdx,
    "#it{p}_{T} (GeV/#it{c})", "(1/#it{N}_{evt}) d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}",
    0, 200, ylo, yhi,
    0.5, 1.5,
    "Track #it{p}_{T} distribution", "trackPt",
    false, true);
}

// 14. Track eta distribution (event-normalized)
void PlotTrackEta() {
  if (!kDrawTrackEta) return;
  cout << "\n=== Track eta distribution ===" << endl;
  std::vector<Source> sources = kSources;
  if (kDrawData) sources.push_back(kDataSource);
  int refIdx = kDrawData ? (int)sources.size() - 1 : kRefIndex;
  std::vector<TH1*> hists;
  for (size_t i = 0; i < sources.size(); i++) {
    cout << " Source " << i << ": " << sources[i].label << endl;
    TDirectory* d = OpenDir(sources[i].file, kDirTrackEff);
    TH1* h = nullptr;
    if (d) {
      TH2* h2 = (TH2*)d->Get("h2_centrality_track_eta");
      if (h2) {
        h = h2->ProjectionY(Form("hTrackEta_%zu", i));
        h->SetDirectory(nullptr);
      }
      if (!h) { cout << "  No track eta histogram found" << endl; hists.push_back(nullptr); continue; }
      double nEvents = 1.0;
      TH1* hColl = (TH1*)d->Get("h_collisions");
      if (hColl) {
        for (int b = hColl->GetNbinsX(); b >= 1; b--) {
          if (hColl->GetBinContent(b) > 0) { nEvents = hColl->GetBinContent(b); break; }
        }
      }
      if (nEvents > 0) h->Scale(1.0 / nEvents);
    }
    hists.push_back(h);
  }
  DrawHistComparison(hists, sources, refIdx,
    "#eta", "(1/#it{N}_{evt}) d#it{N}/d#eta",
    -1.0, 1.0, 0, 0,
    0.95, 1.05,
    "Track #eta distribution", "trackEta",
    false, false);
}

// 15. Track phi distribution (event-normalized)
void PlotTrackPhi() {
  if (!kDrawTrackPhi) return;
  cout << "\n=== Track phi distribution ===" << endl;
  std::vector<Source> sources = kSources;
  if (kDrawData) sources.push_back(kDataSource);
  int refIdx = kDrawData ? (int)sources.size() - 1 : kRefIndex;
  std::vector<TH1*> hists;
  for (size_t i = 0; i < sources.size(); i++) {
    cout << " Source " << i << ": " << sources[i].label << endl;
    TDirectory* d = OpenDir(sources[i].file, kDirTrackEff);
    TH1* h = nullptr;
    if (d) {
      TH2* h2 = (TH2*)d->Get("h2_centrality_track_phi");
      if (h2) {
        h = h2->ProjectionY(Form("hTrackPhi_%zu", i));
        h->SetDirectory(nullptr);
      }
      if (!h) { cout << "  No track phi histogram found" << endl; hists.push_back(nullptr); continue; }
      double nEvents = 1.0;
      TH1* hColl = (TH1*)d->Get("h_collisions");
      if (hColl) {
        for (int b = hColl->GetNbinsX(); b >= 1; b--) {
          if (hColl->GetBinContent(b) > 0) { nEvents = hColl->GetBinContent(b); break; }
        }
      }
      if (nEvents > 0) h->Scale(1.0 / nEvents);
    }
    hists.push_back(h);
  }
  DrawHistComparison(hists, sources, refIdx,
    "#varphi (rad)", "(1/#it{N}_{evt}) d#it{N}/d#varphi",
    0, TMath::TwoPi(), 0, 0,
    0.95, 1.05,
    "Track #varphi distribution", "trackPhi",
    false, false);
}

// ============================================================
//  MAIN
// ============================================================
void DrawTrackResolutionsNew() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  EnsureDir(kOutputDir);

  cout << "=== TrackTuner QA: Comparing " << kSources.size() << " MC sources ===" << endl;
  for (size_t i = 0; i < kSources.size(); i++) {
    cout << "  [" << i << "] " << kSources[i].label << " : " << kSources[i].file << endl;
  }
  cout << "  Reference: [" << kRefIndex << "] " << kSources[kRefIndex].label << endl;
  cout << "  Output: " << kOutputDir << "/" << endl;

  PlotSigmaPtOverPt();
  PlotSigma1OverPt();
  PlotPtResidualSigma();
  PlotDCAxyWidth();
  PlotDCAzWidth();
  PlotDCAxyMean();
  PlotDCAzMean();
  PlotJetPtDet();
  PlotJetPtUESub();
  PlotJetPtParticle();
  PlotDeltaPtSigma();
  PlotPtResidualDistributions();
  PlotTrackPt();
  PlotTrackEta();
  PlotTrackPhi();

  cout << "\n=== All plots done ===" << endl;
}
