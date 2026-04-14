///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Draw Track Selection QA        //////////
////////// Data vs MC comparison          //////////
////////// author: Joonsuk Bae            //////////
////////// E-mail: jbae@cern.ch           //////////
////////// Last Modified: 2026-02-25      //////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

#include "../Filipad2.h"
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>
#include <TMath.h>
#include <TStyle.h>
#include <TString.h>
#include <iostream>
#include <vector>

// ============================================================
// Configuration
// ============================================================

// Files with processTrackSelectionHistograms enabled
TString kMainDir = "~/cernbox/workspace/O2Physics/jets/AnalysisResults/";
TString kDataFile = "619713_AnalysisResults.root";  // Data: DCAz 0.2 cm variation (has trackselplot)
TString kMCFile   = "619512_AnalysisResults.root";  // MC: DCAz 0.2 cm variation (has trackselplot)
// Alternative: use LHC26b5 local MC
// TString kMCFile = "LHC26b5_local_pTsmearing1p5_AnalysisResults.root";

// File for efficiency (has 3D eff histograms but no trackselplot)
TString kMCEffFile = "596836_AnalysisResults.root";  // LHC25a2b, tuner 1.5

TString kHistDir = "track-efficiency";
TString kOutputDir = "plots/TrackQA";

// Plot ranges for pT axis
double kPtMin = 0.0;
double kPtMax = 200.0;

// Data/MC label
TString kDataLabel = "Data (619713)";
TString kMCLabel   = "MC (619512)";

// ============================================================
// Track QA variable definition
// ============================================================
struct TrackQAVar {
  TString hist1D;
  TString hist2D;
  TString xTitle;
  TString shortName;  // for filenames
  double xMin, xMax;
  double ratioMin, ratioMax;
  bool logY1D;  // log scale for 1D overlay
};

std::vector<TrackQAVar> GetTrackQAVariables() {
  std::vector<TrackQAVar> vars;
  vars.push_back({"h_trackselplot_tpccrossedrows",
                   "h2_trackselplot_pt_tpccrossedrows",
                   "TPC crossed rows", "TPCCrossedRows",
                   -0.5, 164.5, 0.5, 1.5, false});
  vars.push_back({"h_trackselplot_tpccrossedrowsoverfindable",
                   "h2_trackselplot_pt_tpccrossedrowsoverfindable",
                   "TPC crossed rows / findable", "TPCCrossedRowsOverFindable",
                   0.0, 1.2, 0.5, 1.5, true});
  vars.push_back({"h_trackselplot_chi2ncls_tpc",
                   "h2_trackselplot_pt_chi2ncls_tpc",
                   "TPC #chi^{2}/cluster", "Chi2TPC",
                   0.0, 10.0, 0.5, 1.5, false});
  vars.push_back({"h_trackselplot_chi2ncls_its",
                   "h2_trackselplot_pt_chi2ncls_its",
                   "ITS #chi^{2}/cluster", "Chi2ITS",
                   0.0, 40.0, 0.5, 1.5, true});
  vars.push_back({"h_trackselplot_dcaxy",
                   "h2_trackselplot_pt_dcaxy",
                   "DCA_{xy} (cm)", "DCAxy",
                   -1.0, 1.0, 0.0, 2.5, true});
  vars.push_back({"h_trackselplot_dcaz",
                   "h2_trackselplot_pt_dcaz",
                   "DCA_{z} (cm)", "DCAz",
                   -4.0, 4.0, 0.0, 2.5, true});
  return vars;
}

// ============================================================
// Helper Functions
// ============================================================

void hset(TH1& h, TString xtitle, TString ytitle,
          double xTitleOffset, double yTitleOffset,
          double xLabelSize, double yLabelSize,
          double xLabelOffset, double yLabelOffset,
          double xTickLength, double yTickLength,
          int xNDiv, int yNDiv) {
  h.GetXaxis()->SetTitle(xtitle);
  h.GetYaxis()->SetTitle(ytitle);
  h.GetXaxis()->SetTitleOffset(xTitleOffset);
  h.GetYaxis()->SetTitleOffset(yTitleOffset);
  h.GetXaxis()->SetLabelSize(xLabelSize);
  h.GetYaxis()->SetLabelSize(yLabelSize);
  h.GetXaxis()->SetLabelOffset(xLabelOffset);
  h.GetYaxis()->SetLabelOffset(yLabelOffset);
  h.GetXaxis()->SetTickLength(xTickLength);
  h.GetYaxis()->SetTickLength(yTickLength);
  h.GetXaxis()->SetNdivisions(xNDiv);
  h.GetYaxis()->SetNdivisions(yNDiv);
  h.GetXaxis()->CenterTitle(1);
  h.GetYaxis()->CenterTitle(1);
}

void optFili(TPad& pad, int gridx, int gridy, int logx, int logy) {
  pad.SetGridx(gridx);
  pad.SetGridy(gridy);
  pad.SetLogx(logx);
  pad.SetLogy(logy);
}

// Safe histogram retrieval: get from directory, detach from file
TH1* GetHist1D(TDirectory* dir, TString name) {
  TH1* h = (TH1*)dir->Get(name);
  if (!h) { std::cerr << "[Warning] Histogram not found: " << name << std::endl; return nullptr; }
  h->SetDirectory(0);
  return h;
}

TH2* GetHist2D(TDirectory* dir, TString name) {
  TH2* h = (TH2*)dir->Get(name);
  if (!h) { std::cerr << "[Warning] 2D Histogram not found: " << name << std::endl; return nullptr; }
  h->SetDirectory(0);
  return h;
}

TH3* GetHist3D(TDirectory* dir, TString name) {
  TH3* h = (TH3*)dir->Get(name);
  if (!h) { std::cerr << "[Warning] 3D Histogram not found: " << name << std::endl; return nullptr; }
  h->SetDirectory(0);
  return h;
}

// ============================================================
// Draw 1D overlay: Data vs MC (normalized to unit area) + ratio
// ============================================================
void Draw1DOverlays(TDirectory* dirData, TDirectory* dirMC) {
  int nn = 0;
  std::vector<TrackQAVar> vars = GetTrackQAVariables();

  for (int iv = 0; iv < (int)vars.size(); iv++) {
    TrackQAVar& v = vars[iv];
    TH1* hData = GetHist1D(dirData, v.hist1D);
    TH1* hMC   = GetHist1D(dirMC,   v.hist1D);
    if (!hData || !hMC) continue;

    // Normalize to unit area
    if (hData->Integral() > 0) hData->Scale(1.0 / hData->Integral());
    if (hMC->Integral() > 0)   hMC->Scale(1.0 / hMC->Integral());

    // Create Filipad2 canvas
    Filipad2* fpad = new Filipad2(Form("c1D_%s", v.shortName.Data()), ++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    fpad->Draw();

    // Upper pad: overlay
    TPad* p1 = fpad->GetPad(1);
    p1->cd();
    p1->SetTickx(); p1->SetTicky();
    if (v.logY1D) p1->SetLogy();

    hset(*hData, "", "Normalized",
         2.5, 1.5, 0.05, 0.05, 0.01, 0.001, 0.03, 0.03, 510, 510);
    hData->GetXaxis()->SetRangeUser(v.xMin, v.xMax);
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(0.6);
    hData->SetMarkerColor(kBlack);
    hData->SetLineColor(kBlack);
    hData->SetLineWidth(1);
    hData->SetStats(0);

    hMC->SetMarkerStyle(24);
    hMC->SetMarkerSize(0.6);
    hMC->SetMarkerColor(kRed);
    hMC->SetLineColor(kRed);
    hMC->SetLineWidth(1);
    hMC->SetStats(0);

    // Find Y range
    double ymax = TMath::Max(hData->GetMaximum(), hMC->GetMaximum()) * 1.5;
    if (v.logY1D) {
      hData->SetMinimum(1e-6);
      hData->SetMaximum(ymax * 5);
    } else {
      hData->SetMinimum(0);
      hData->SetMaximum(ymax);
    }

    hData->DrawCopy("PE");
    hMC->DrawCopy("PE SAME");

    TLegend* leg = new TLegend(0.55, 0.70, 0.92, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
    leg->AddEntry(hData, kDataLabel, "lp");
    leg->AddEntry(hMC, kMCLabel, "lp");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.DrawLatex(0.20, 0.85, v.xTitle);

    // Lower pad: ratio
    TPad* p2 = fpad->GetPad(2);
    p2->cd();
    p2->SetTickx(); p2->SetTicky();
    p2->SetGridy();

    TH1* hRatio = (TH1*)hData->Clone(Form("hRatio_%s", v.shortName.Data()));
    hRatio->Divide(hMC);
    hRatio->SetDirectory(0);

    hset(*hRatio, v.xTitle, "Data/MC",
         2.5, 0.8, 0.10, 0.10, 0.01, 0.001, 0.05, 0.03, 510, 510);
    hRatio->GetXaxis()->SetRangeUser(v.xMin, v.xMax);
    hRatio->SetMinimum(v.ratioMin);
    hRatio->SetMaximum(v.ratioMax);
    hRatio->DrawCopy("PE");

    TLine* line = new TLine(v.xMin, 1.0, v.xMax, 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kGray + 2);
    line->Draw();

    // Save
    fpad->C->SaveAs(Form("%s/1D_%s.pdf", kOutputDir.Data(), v.shortName.Data()));

    delete hData; delete hMC; delete hRatio;
  }
}

// ============================================================
// Draw 2D (pT vs property) COLZ: Data and MC side by side
// ============================================================
void Draw2DComparisons(TDirectory* dirData, TDirectory* dirMC) {
  std::vector<TrackQAVar> vars = GetTrackQAVariables();

  for (int iv = 0; iv < (int)vars.size(); iv++) {
    TrackQAVar& v = vars[iv];
    TH2* hData = GetHist2D(dirData, v.hist2D);
    TH2* hMC   = GetHist2D(dirMC,   v.hist2D);
    if (!hData || !hMC) continue;

    TCanvas* c = new TCanvas(Form("c2D_%s", v.shortName.Data()),
                              Form("2D %s", v.shortName.Data()),
                              1200, 500);
    c->Divide(2, 1);

    // Data
    c->cd(1);
    gPad->SetLogx();
    gPad->SetLogz();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.15);
    hData->SetTitle(Form("Data: %s", v.xTitle.Data()));
    hData->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hData->GetYaxis()->SetTitle(v.xTitle);
    hData->GetXaxis()->SetRangeUser(0.15, kPtMax);
    hData->GetYaxis()->SetRangeUser(v.xMin, v.xMax);
    hData->SetStats(0);
    hData->DrawCopy("COLZ");

    // MC
    c->cd(2);
    gPad->SetLogx();
    gPad->SetLogz();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.15);
    hMC->SetTitle(Form("MC: %s", v.xTitle.Data()));
    hMC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hMC->GetYaxis()->SetTitle(v.xTitle);
    hMC->GetXaxis()->SetRangeUser(0.15, kPtMax);
    hMC->GetYaxis()->SetRangeUser(v.xMin, v.xMax);
    hMC->SetStats(0);
    hMC->DrawCopy("COLZ");

    c->SaveAs(Form("%s/2D_%s.pdf", kOutputDir.Data(), v.shortName.Data()));
    delete hData; delete hMC;
    delete c;
  }
}

// ============================================================
// Draw profile plots: mean of property vs pT, Data vs MC + ratio
// ============================================================
void DrawProfiles(TDirectory* dirData, TDirectory* dirMC) {
  int nn = 100;
  std::vector<TrackQAVar> vars = GetTrackQAVariables();

  for (int iv = 0; iv < (int)vars.size(); iv++) {
    TrackQAVar& v = vars[iv];
    TH2* hData2D = GetHist2D(dirData, v.hist2D);
    TH2* hMC2D   = GetHist2D(dirMC,   v.hist2D);
    if (!hData2D || !hMC2D) continue;

    TProfile* profData = hData2D->ProfileX(Form("profData_%s", v.shortName.Data()));
    TProfile* profMC   = hMC2D->ProfileX(Form("profMC_%s", v.shortName.Data()));
    profData->SetDirectory(0);
    profMC->SetDirectory(0);

    // Create Filipad2 canvas
    Filipad2* fpad = new Filipad2(Form("cProf_%s", v.shortName.Data()), ++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    fpad->Draw();

    // Upper pad: overlay profiles
    TPad* p1 = fpad->GetPad(1);
    p1->cd();
    p1->SetTickx(); p1->SetTicky();
    p1->SetLogx();

    hset(*profData, "", Form("#LT%s#GT", v.xTitle.Data()),
         2.5, 1.5, 0.05, 0.05, 0.01, 0.001, 0.03, 0.03, 510, 510);
    profData->GetXaxis()->SetRangeUser(0.15, kPtMax);
    profData->SetMarkerStyle(20);
    profData->SetMarkerSize(0.6);
    profData->SetMarkerColor(kBlack);
    profData->SetLineColor(kBlack);
    profData->SetStats(0);

    profMC->SetMarkerStyle(24);
    profMC->SetMarkerSize(0.6);
    profMC->SetMarkerColor(kRed);
    profMC->SetLineColor(kRed);
    profMC->SetStats(0);

    profData->DrawCopy("PE");
    profMC->DrawCopy("PE SAME");

    TLegend* leg = new TLegend(0.55, 0.70, 0.92, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
    leg->AddEntry(profData, kDataLabel, "lp");
    leg->AddEntry(profMC, kMCLabel, "lp");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.DrawLatex(0.20, 0.85, Form("#LT%s#GT vs #it{p}_{T}", v.xTitle.Data()));

    // Lower pad: ratio
    TPad* p2 = fpad->GetPad(2);
    p2->cd();
    p2->SetTickx(); p2->SetTicky();
    p2->SetLogx();
    p2->SetGridy();

    // Build ratio from profiles: clone Data profile, divide by MC
    TH1* hProfRatio = (TH1*)profData->Clone(Form("hProfRatio_%s", v.shortName.Data()));
    hProfRatio->Divide(profMC);
    hProfRatio->SetDirectory(0);

    hset(*hProfRatio, "#it{p}_{T} (GeV/#it{c})", "Data/MC",
         2.5, 0.8, 0.10, 0.10, 0.01, 0.001, 0.05, 0.03, 510, 510);
    hProfRatio->GetXaxis()->SetRangeUser(0.15, kPtMax);
    hProfRatio->SetMinimum(0.9);
    hProfRatio->SetMaximum(1.1);
    hProfRatio->DrawCopy("PE");

    TLine* line = new TLine(0.15, 1.0, kPtMax, 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kGray + 2);
    line->Draw();

    fpad->C->SaveAs(Form("%s/Profile_%s.pdf", kOutputDir.Data(), v.shortName.Data()));

    delete hData2D; delete hMC2D; delete profData; delete profMC; delete hProfRatio;
  }
}

// ============================================================
// Draw single-track efficiency vs pT from 3D histograms
// ============================================================
void DrawEfficiency(TDirectory* dirMCEff) {
  // Low pT range [0, 10] GeV
  TH3* h3Truth_lo = GetHist3D(dirMCEff, "h3_particle_pt_particle_eta_particle_phi_mcpartofinterest");
  TH3* h3Match_lo = GetHist3D(dirMCEff, "h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary");
  // High pT range [10, 100] GeV
  TH3* h3Truth_hi = GetHist3D(dirMCEff, "h3_particle_pt_high_particle_eta_particle_phi_mcpartofinterest");
  TH3* h3Match_hi = GetHist3D(dirMCEff, "h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_primary");

  if (!h3Truth_lo || !h3Match_lo) {
    std::cerr << "[Warning] Cannot draw efficiency: missing 3D histograms in MC file" << std::endl;
    return;
  }

  // Restrict eta range to |eta| < 0.9 (track acceptance)
  // Eta axis: 200 bins, [-1, 1] -> bin width 0.01
  // |eta| < 0.9 -> bins [11, 190] (0-indexed: bin 1=-1.0, bin 10=-0.91, bin 11=-0.90, ..., bin 190=0.89, bin 191=0.90)
  int etaBinLo = h3Truth_lo->GetYaxis()->FindBin(-0.899);
  int etaBinHi = h3Truth_lo->GetYaxis()->FindBin(0.899);

  // Project to pT axis with eta cut
  TH1* hTruth_lo = h3Truth_lo->ProjectionX("hTruth_lo", etaBinLo, etaBinHi, 0, -1);
  TH1* hMatch_lo = h3Match_lo->ProjectionX("hMatch_lo", etaBinLo, etaBinHi, 0, -1);
  hTruth_lo->SetDirectory(0);
  hMatch_lo->SetDirectory(0);

  // Compute efficiency (low pT)
  TH1* hEff_lo = (TH1*)hMatch_lo->Clone("hEff_lo");
  hEff_lo->SetDirectory(0);
  hEff_lo->Divide(hMatch_lo, hTruth_lo, 1., 1., "B");  // Binomial errors

  TH1* hEff_hi = nullptr;
  if (h3Truth_hi && h3Match_hi) {
    int etaBinLo_hi = h3Truth_hi->GetYaxis()->FindBin(-0.899);
    int etaBinHi_hi = h3Truth_hi->GetYaxis()->FindBin(0.899);
    TH1* hTruth_hi = h3Truth_hi->ProjectionX("hTruth_hi", etaBinLo_hi, etaBinHi_hi, 0, -1);
    TH1* hMatch_hi = h3Match_hi->ProjectionX("hMatch_hi", etaBinLo_hi, etaBinHi_hi, 0, -1);
    hTruth_hi->SetDirectory(0);
    hMatch_hi->SetDirectory(0);

    hEff_hi = (TH1*)hMatch_hi->Clone("hEff_hi");
    hEff_hi->SetDirectory(0);
    hEff_hi->Divide(hMatch_hi, hTruth_hi, 1., 1., "B");

    delete hTruth_hi; delete hMatch_hi;
  }

  // Draw efficiency
  TCanvas* cEff = new TCanvas("cEff", "Track Efficiency", 800, 600);
  cEff->SetTickx(); cEff->SetTicky();
  cEff->SetLogx();
  cEff->SetLeftMargin(0.12);
  cEff->SetBottomMargin(0.12);

  hset(*hEff_lo, "#it{p}_{T} (GeV/#it{c})", "Tracking efficiency",
       1.2, 1.4, 0.04, 0.04, 0.01, 0.001, 0.03, 0.03, 510, 510);
  hEff_lo->GetXaxis()->SetRangeUser(0.15, 10.0);
  hEff_lo->SetMinimum(0.0);
  hEff_lo->SetMaximum(1.1);
  hEff_lo->SetMarkerStyle(20);
  hEff_lo->SetMarkerSize(0.5);
  hEff_lo->SetMarkerColor(kBlue + 1);
  hEff_lo->SetLineColor(kBlue + 1);
  hEff_lo->SetStats(0);
  hEff_lo->DrawCopy("PE");

  if (hEff_hi) {
    hEff_hi->SetMarkerStyle(20);
    hEff_hi->SetMarkerSize(0.7);
    hEff_hi->SetMarkerColor(kRed + 1);
    hEff_hi->SetLineColor(kRed + 1);
    // Draw high pT on same canvas — shift X range
    // High pT: [10, 100] GeV, 18 bins of 5 GeV
    hEff_hi->DrawCopy("PE SAME");
  }

  // Draw combined on wider range
  TCanvas* cEffWide = new TCanvas("cEffWide", "Track Efficiency (full pT)", 800, 600);
  cEffWide->SetTickx(); cEffWide->SetTicky();
  cEffWide->SetLogx();
  cEffWide->SetLeftMargin(0.12);
  cEffWide->SetBottomMargin(0.12);

  // Create combined histogram
  // Low: 200 bins [0, 10] = 50 MeV; High: 18 bins [10, 100] = 5 GeV
  // Use variable binning
  const int nBinsCombined = 200 + 18;
  double xbins[nBinsCombined + 1];
  // Low pT: 200 bins, 0 to 10
  for (int i = 0; i <= 200; i++) xbins[i] = 0.0 + i * 0.05;
  // High pT: 18 bins, 10 to 100
  for (int i = 1; i <= 18; i++) xbins[200 + i] = 10.0 + i * 5.0;

  TH1D* hEffCombined = new TH1D("hEffCombined", "", nBinsCombined, xbins);
  hEffCombined->SetDirectory(0);

  // Fill from low pT
  for (int ib = 1; ib <= hEff_lo->GetNbinsX(); ib++) {
    double pt = hEff_lo->GetBinCenter(ib);
    int binC = hEffCombined->FindBin(pt);
    hEffCombined->SetBinContent(binC, hEff_lo->GetBinContent(ib));
    hEffCombined->SetBinError(binC, hEff_lo->GetBinError(ib));
  }
  // Fill from high pT
  if (hEff_hi) {
    for (int ib = 1; ib <= hEff_hi->GetNbinsX(); ib++) {
      double pt = hEff_hi->GetBinCenter(ib);
      int binC = hEffCombined->FindBin(pt);
      hEffCombined->SetBinContent(binC, hEff_hi->GetBinContent(ib));
      hEffCombined->SetBinError(binC, hEff_hi->GetBinError(ib));
    }
  }

  hset(*hEffCombined, "#it{p}_{T} (GeV/#it{c})", "Tracking efficiency",
       1.2, 1.4, 0.04, 0.04, 0.01, 0.001, 0.03, 0.03, 510, 510);
  hEffCombined->GetXaxis()->SetRangeUser(0.15, 100.0);
  hEffCombined->SetMinimum(0.0);
  hEffCombined->SetMaximum(1.1);
  hEffCombined->SetMarkerStyle(20);
  hEffCombined->SetMarkerSize(0.5);
  hEffCombined->SetMarkerColor(kBlack);
  hEffCombined->SetLineColor(kBlack);
  hEffCombined->SetStats(0);
  hEffCombined->DrawCopy("PE");

  // Reference lines at 1.0 and at +-3% variation
  TLine* l1 = new TLine(0.15, 1.0, 100.0, 1.0);
  l1->SetLineStyle(2); l1->SetLineColor(kGray + 1); l1->Draw();
  TLine* l097 = new TLine(0.15, 0.97, 100.0, 0.97);
  l097->SetLineStyle(3); l097->SetLineColor(kGray + 1);
  // Don't draw 0.97 line — it would be too low on efficiency plot

  TLegend* legEff = new TLegend(0.15, 0.15, 0.55, 0.30);
  legEff->SetBorderSize(0);
  legEff->SetFillStyle(0);
  legEff->SetTextSize(0.035);
  legEff->AddEntry(hEffCombined, "Primary track eff. (|#eta| < 0.9)", "lp");
  legEff->Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04);
  latex.DrawLatex(0.15, 0.32, Form("MC: %s", kMCEffFile.Data()));

  cEffWide->SaveAs(Form("%s/Efficiency_vs_pT.pdf", kOutputDir.Data()));
  cEff->SaveAs(Form("%s/Efficiency_vs_pT_lowHigh.pdf", kOutputDir.Data()));

  // Also draw fake rate
  TH3* h3Fake_lo = GetHist3D(dirMCEff, "h3_track_pt_track_eta_track_phi_nonassociatedtrack");
  TH3* h3All_lo  = GetHist3D(dirMCEff, "h3_track_pt_track_eta_track_phi_associatedtrack_primary");
  if (h3Fake_lo && h3All_lo) {
    // Fake rate = non-associated / (associated_primary + non-associated)
    int etaLo = h3Fake_lo->GetYaxis()->FindBin(-0.899);
    int etaHi = h3Fake_lo->GetYaxis()->FindBin(0.899);
    TH1* hFake = h3Fake_lo->ProjectionX("hFake", etaLo, etaHi, 0, -1);
    TH1* hAllTracks = h3All_lo->ProjectionX("hAllTracks", etaLo, etaHi, 0, -1);
    hFake->SetDirectory(0);
    hAllTracks->SetDirectory(0);

    // Denominator = fake + primary associated
    TH1* hDenom = (TH1*)hFake->Clone("hDenom");
    hDenom->Add(hAllTracks);
    hDenom->SetDirectory(0);

    TH1* hFakeRate = (TH1*)hFake->Clone("hFakeRate");
    hFakeRate->SetDirectory(0);
    hFakeRate->Divide(hFake, hDenom, 1., 1., "B");

    TCanvas* cFake = new TCanvas("cFake", "Fake Rate", 800, 600);
    cFake->SetLogx(); cFake->SetTickx(); cFake->SetTicky();
    cFake->SetLeftMargin(0.12); cFake->SetBottomMargin(0.12);

    hset(*hFakeRate, "#it{p}_{T} (GeV/#it{c})", "Fake track rate",
         1.2, 1.4, 0.04, 0.04, 0.01, 0.001, 0.03, 0.03, 510, 510);
    hFakeRate->GetXaxis()->SetRangeUser(0.15, 10.0);
    hFakeRate->SetMinimum(0.0);
    hFakeRate->SetMaximum(0.15);
    hFakeRate->SetMarkerStyle(20);
    hFakeRate->SetMarkerSize(0.5);
    hFakeRate->SetMarkerColor(kRed + 1);
    hFakeRate->SetLineColor(kRed + 1);
    hFakeRate->SetStats(0);
    hFakeRate->DrawCopy("PE");

    TLegend* legF = new TLegend(0.50, 0.75, 0.88, 0.88);
    legF->SetBorderSize(0); legF->SetFillStyle(0); legF->SetTextSize(0.035);
    legF->AddEntry(hFakeRate, "Fake / (Fake + Primary)", "lp");
    legF->Draw();

    cFake->SaveAs(Form("%s/FakeRate_vs_pT.pdf", kOutputDir.Data()));

    delete hFake; delete hAllTracks; delete hDenom; delete hFakeRate; delete cFake;
  }

  // Also draw secondary contamination
  TH3* h3Sec_lo = GetHist3D(dirMCEff, "h3_track_pt_track_eta_track_phi_associatedtrack_nonprimary");
  if (h3Sec_lo && h3All_lo) {
    int etaLo = h3Sec_lo->GetYaxis()->FindBin(-0.899);
    int etaHi = h3Sec_lo->GetYaxis()->FindBin(0.899);
    TH1* hSec = h3Sec_lo->ProjectionX("hSec", etaLo, etaHi, 0, -1);
    TH1* hPrim = h3All_lo->ProjectionX("hPrimForSec", etaLo, etaHi, 0, -1);
    hSec->SetDirectory(0);
    hPrim->SetDirectory(0);

    TH1* hDenom = (TH1*)hSec->Clone("hDenomSec");
    hDenom->Add(hPrim);
    hDenom->SetDirectory(0);

    TH1* hSecFrac = (TH1*)hSec->Clone("hSecFrac");
    hSecFrac->SetDirectory(0);
    hSecFrac->Divide(hSec, hDenom, 1., 1., "B");

    TCanvas* cSec = new TCanvas("cSec", "Secondary Fraction", 800, 600);
    cSec->SetLogx(); cSec->SetTickx(); cSec->SetTicky();
    cSec->SetLeftMargin(0.12); cSec->SetBottomMargin(0.12);

    hset(*hSecFrac, "#it{p}_{T} (GeV/#it{c})", "Secondary fraction",
         1.2, 1.4, 0.04, 0.04, 0.01, 0.001, 0.03, 0.03, 510, 510);
    hSecFrac->GetXaxis()->SetRangeUser(0.15, 10.0);
    hSecFrac->SetMinimum(0.0);
    hSecFrac->SetMaximum(0.15);
    hSecFrac->SetMarkerStyle(20);
    hSecFrac->SetMarkerSize(0.5);
    hSecFrac->SetMarkerColor(kOrange + 7);
    hSecFrac->SetLineColor(kOrange + 7);
    hSecFrac->SetStats(0);
    hSecFrac->DrawCopy("PE");

    TLegend* legS = new TLegend(0.50, 0.75, 0.88, 0.88);
    legS->SetBorderSize(0); legS->SetFillStyle(0); legS->SetTextSize(0.035);
    legS->AddEntry(hSecFrac, "Secondary / (Secondary + Primary)", "lp");
    legS->Draw();

    cSec->SaveAs(Form("%s/SecondaryFraction_vs_pT.pdf", kOutputDir.Data()));

    delete hSec; delete hPrim; delete hDenom; delete hSecFrac; delete cSec;
  }

  // Cleanup
  delete h3Truth_lo; delete h3Match_lo;
  delete hTruth_lo; delete hMatch_lo; delete hEff_lo;
  if (h3Truth_hi) delete h3Truth_hi;
  if (h3Match_hi) delete h3Match_hi;
  if (hEff_hi) delete hEff_hi;
  delete hEffCombined;
  delete cEff; delete cEffWide;
  if (h3Fake_lo) delete h3Fake_lo;
  if (h3All_lo) delete h3All_lo;
  if (h3Sec_lo) delete h3Sec_lo;
}

// ============================================================
// Draw pT resolution and sigma(pT) from 2D histograms
// ============================================================
void DrawPtResolution(TDirectory* dirMCEff) {
  // pT resolution: (pT_reco - pT_truth) / pT_truth vs pT_truth
  TH2* h2Res = GetHist2D(dirMCEff, "h2_particle_pt_track_pt_deltaptoverparticlept");
  if (!h2Res) return;

  // Profile: mean resolution vs pT
  TProfile* profRes = h2Res->ProfileX("profRes");
  profRes->SetDirectory(0);

  TCanvas* cRes = new TCanvas("cRes", "pT Resolution", 800, 600);
  cRes->SetLogx(); cRes->SetTickx(); cRes->SetTicky();
  cRes->SetLeftMargin(0.12); cRes->SetBottomMargin(0.12);

  hset(*profRes, "#it{p}_{T}^{truth} (GeV/#it{c})", "#LT#Delta#it{p}_{T}/#it{p}_{T}^{truth}#GT",
       1.2, 1.4, 0.04, 0.04, 0.01, 0.001, 0.03, 0.03, 510, 510);
  profRes->GetXaxis()->SetRangeUser(0.15, 200.0);
  profRes->SetMinimum(-0.05);
  profRes->SetMaximum(0.05);
  profRes->SetMarkerStyle(20);
  profRes->SetMarkerSize(0.5);
  profRes->SetMarkerColor(kBlue + 1);
  profRes->SetLineColor(kBlue + 1);
  profRes->SetStats(0);
  profRes->DrawCopy("PE");

  TLine* l0 = new TLine(0.15, 0.0, 200.0, 0.0);
  l0->SetLineStyle(2); l0->SetLineColor(kGray + 1); l0->Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.15, 0.85, Form("MC: %s", kMCEffFile.Data()));

  cRes->SaveAs(Form("%s/PtResolution_vs_pT.pdf", kOutputDir.Data()));

  // 2D COLZ
  TCanvas* cRes2D = new TCanvas("cRes2D", "pT Resolution 2D", 700, 600);
  cRes2D->SetLogx(); cRes2D->SetLogz();
  cRes2D->SetLeftMargin(0.12); cRes2D->SetRightMargin(0.15); cRes2D->SetBottomMargin(0.12);
  hset(*h2Res, "#it{p}_{T}^{truth} (GeV/#it{c})", "#Delta#it{p}_{T}/#it{p}_{T}^{truth}",
       1.2, 1.4, 0.04, 0.04, 0.01, 0.001, 0.03, 0.03, 510, 510);
  h2Res->GetXaxis()->SetRangeUser(0.15, 200.0);
  h2Res->SetStats(0);
  h2Res->DrawCopy("COLZ");

  cRes2D->SaveAs(Form("%s/PtResolution_2D.pdf", kOutputDir.Data()));

  // sigma(pT)/pT vs pT
  TH2* h2Sigma = GetHist2D(dirMCEff, "h2_track_pt_track_sigmapt");
  if (h2Sigma) {
    TProfile* profSigma = h2Sigma->ProfileX("profSigma");
    profSigma->SetDirectory(0);

    TCanvas* cSigma = new TCanvas("cSigma", "sigma pT", 800, 600);
    cSigma->SetLogx(); cSigma->SetTickx(); cSigma->SetTicky();
    cSigma->SetLeftMargin(0.12); cSigma->SetBottomMargin(0.12);

    hset(*profSigma, "#it{p}_{T} (GeV/#it{c})", "#LT#sigma(#it{p}_{T})/#it{p}_{T}#GT (%)",
         1.2, 1.4, 0.04, 0.04, 0.01, 0.001, 0.03, 0.03, 510, 510);
    profSigma->GetXaxis()->SetRangeUser(0.15, 10.0);
    profSigma->SetMinimum(0.0);
    profSigma->SetMaximum(10.0);
    profSigma->SetMarkerStyle(20);
    profSigma->SetMarkerSize(0.5);
    profSigma->SetMarkerColor(kGreen + 2);
    profSigma->SetLineColor(kGreen + 2);
    profSigma->SetStats(0);
    profSigma->DrawCopy("PE");

    cSigma->SaveAs(Form("%s/SigmaPt_vs_pT.pdf", kOutputDir.Data()));

    delete profSigma; delete h2Sigma; delete cSigma;
  }

  delete h2Res; delete profRes; delete cRes; delete cRes2D;
}

// ============================================================
// Draw cut flow from hTrackCutsCounts
// ============================================================
void DrawCutFlow(TDirectory* dirData, TDirectory* dirMC) {
  TH1* hCutsData = dirData ? GetHist1D(dirData, "hTrackCutsCounts") : nullptr;
  TH1* hCutsMC   = dirMC   ? GetHist1D(dirMC,   "hTrackCutsCounts") : nullptr;

  // Draw whichever is available
  TH1* hCuts = hCutsMC ? hCutsMC : hCutsData;
  if (!hCuts) {
    std::cerr << "[Info] No hTrackCutsCounts found, skipping cut flow" << std::endl;
    return;
  }

  TCanvas* cCuts = new TCanvas("cCuts", "Track Cut Flow", 900, 600);
  cCuts->SetLeftMargin(0.12); cCuts->SetBottomMargin(0.20);
  cCuts->SetLogy();

  // Normalize to first bin (all tracks)
  TH1* hCutsNorm = (TH1*)hCuts->Clone("hCutsNorm");
  hCutsNorm->SetDirectory(0);
  double nAll = hCutsNorm->GetBinContent(1);
  if (nAll > 0) hCutsNorm->Scale(1.0 / nAll);

  hCutsNorm->SetMarkerStyle(20);
  hCutsNorm->SetMarkerSize(0.8);
  hCutsNorm->SetMarkerColor(kBlue + 1);
  hCutsNorm->SetLineColor(kBlue + 1);
  hCutsNorm->SetStats(0);
  hCutsNorm->GetYaxis()->SetTitle("Fraction of tracks");
  hCutsNorm->SetMinimum(1e-3);
  hCutsNorm->SetMaximum(2.0);
  hCutsNorm->DrawCopy("PE TEXT45");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  TString label = hCutsMC ? "MC" : "Data";
  latex.DrawLatex(0.15, 0.85, Form("Track cut flow (%s)", label.Data()));

  cCuts->SaveAs(Form("%s/CutFlow.pdf", kOutputDir.Data()));

  delete hCutsNorm;
  if (hCutsData) delete hCutsData;
  if (hCutsMC) delete hCutsMC;
  delete cCuts;
}

// ============================================================
// Save all histograms to a ROOT file
// ============================================================
void SaveToRootFile(TDirectory* dirData, TDirectory* dirMC, TDirectory* dirMCEff) {
  TFile* fout = TFile::Open(Form("%s/TrackQA.root", kOutputDir.Data()), "RECREATE");
  if (!fout) return;

  std::vector<TrackQAVar> vars = GetTrackQAVariables();

  // 1D histograms
  TDirectory* dir1D = fout->mkdir("1D");
  for (auto& v : vars) {
    if (dirData) {
      TH1* h = GetHist1D(dirData, v.hist1D);
      if (h) { dir1D->cd(); h->Write(Form("Data_%s", v.shortName.Data())); delete h; }
    }
    if (dirMC) {
      TH1* h = GetHist1D(dirMC, v.hist1D);
      if (h) { dir1D->cd(); h->Write(Form("MC_%s", v.shortName.Data())); delete h; }
    }
  }

  // 2D histograms
  TDirectory* dir2D = fout->mkdir("2D");
  for (auto& v : vars) {
    if (dirData) {
      TH2* h = GetHist2D(dirData, v.hist2D);
      if (h) { dir2D->cd(); h->Write(Form("Data_%s", v.shortName.Data())); delete h; }
    }
    if (dirMC) {
      TH2* h = GetHist2D(dirMC, v.hist2D);
      if (h) { dir2D->cd(); h->Write(Form("MC_%s", v.shortName.Data())); delete h; }
    }
  }

  fout->Close();
  delete fout;
  std::cerr << "[Info] Saved TrackQA.root" << std::endl;
}

// ============================================================
// Main entry point
// ============================================================
void DrawTrackQA() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // Create output directory
  gSystem->mkdir(kOutputDir, true);

  // --- Part 1: Track selection QA (Data vs MC) ---
  TFile* fData = TFile::Open(kMainDir + kDataFile, "READ");
  TFile* fMC   = TFile::Open(kMainDir + kMCFile, "READ");

  TDirectory* dirData = nullptr;
  TDirectory* dirMC   = nullptr;

  if (fData && !fData->IsZombie()) {
    dirData = (TDirectory*)fData->Get(kHistDir);
    if (!dirData) std::cerr << "[Warning] No " << kHistDir << " directory in Data file" << std::endl;
  } else {
    std::cerr << "[Warning] Cannot open Data file: " << kDataFile << std::endl;
  }

  if (fMC && !fMC->IsZombie()) {
    dirMC = (TDirectory*)fMC->Get(kHistDir);
    if (!dirMC) std::cerr << "[Warning] No " << kHistDir << " directory in MC file" << std::endl;
  } else {
    std::cerr << "[Warning] Cannot open MC file: " << kMCFile << std::endl;
  }

  // Check if trackselplot histograms exist
  bool hasTrackSel = false;
  if (dirData) {
    TH1* hTest = (TH1*)dirData->Get("h_trackselplot_tpccrossedrows");
    if (hTest) { hasTrackSel = true; delete hTest; }
  }
  if (!hasTrackSel && dirMC) {
    TH1* hTest = (TH1*)dirMC->Get("h_trackselplot_tpccrossedrows");
    if (hTest) { hasTrackSel = true; delete hTest; }
  }

  if (hasTrackSel && dirData && dirMC) {
    std::cerr << "[Info] Drawing 1D overlays..." << std::endl;
    Draw1DOverlays(dirData, dirMC);

    std::cerr << "[Info] Drawing 2D comparisons..." << std::endl;
    Draw2DComparisons(dirData, dirMC);

    std::cerr << "[Info] Drawing profile plots..." << std::endl;
    DrawProfiles(dirData, dirMC);

    // Save raw histograms
    SaveToRootFile(dirData, dirMC, nullptr);
  } else {
    std::cerr << "[Warning] Track selection histograms (h_trackselplot_*) not available." << std::endl;
    std::cerr << "  processTrackSelectionHistograms was not enabled in these train configs." << std::endl;
    std::cerr << "  Files with trackselplot: 610951, 611924, 615987, 616559, 619511-619514, 619522, 619713, 619714" << std::endl;
    std::cerr << "  LHC26b5_local_* files also have them." << std::endl;
  }

  // Cut flow (works with any track-efficiency file)
  std::cerr << "[Info] Drawing cut flow..." << std::endl;
  DrawCutFlow(dirData, dirMC);

  // --- Part 2: Efficiency & pT resolution (MC only) ---
  TFile* fMCEff = TFile::Open(kMainDir + kMCEffFile, "READ");
  TDirectory* dirMCEff = nullptr;
  if (fMCEff && !fMCEff->IsZombie()) {
    dirMCEff = (TDirectory*)fMCEff->Get(kHistDir);
  }

  // If the main MC file has efficiency histograms, use it; otherwise use kMCEffFile
  TDirectory* dirEff = dirMCEff;
  if (!dirEff && dirMC) {
    // Check if main MC file has 3D efficiency histograms
    TH3* hTest = (TH3*)dirMC->Get("h3_particle_pt_particle_eta_particle_phi_mcpartofinterest");
    if (hTest) { dirEff = dirMC; delete hTest; }
  }

  if (dirEff) {
    std::cerr << "[Info] Drawing efficiency plots..." << std::endl;
    DrawEfficiency(dirEff);

    std::cerr << "[Info] Drawing pT resolution..." << std::endl;
    DrawPtResolution(dirEff);
  } else {
    std::cerr << "[Warning] No efficiency histograms found in any MC file" << std::endl;
  }

  // Cleanup files
  if (fData) { fData->Close(); delete fData; }
  if (fMC)   { fMC->Close(); delete fMC; }
  if (fMCEff) { fMCEff->Close(); delete fMCEff; }

  std::cerr << "[Info] All done! Output in: " << kOutputDir << std::endl;
}
