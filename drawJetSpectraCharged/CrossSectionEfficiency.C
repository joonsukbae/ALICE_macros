// Compare MCP histograms from two AnalysisResults files
// Compares: 496218 vs 498287
// Histogram: jet-finder-charged-qa/h_jet_pt_part

#if defined(__CLING__) || defined(__CINT__) || defined(__ROOTCLING__)

#include "DrawJetsMCRDependentHelpers.h"
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLatex.h>
#include <TMath.h>
#include <algorithm>
#include <vector>
#include <iostream>

// Constants
const double PlotPtMin = 0;
const double PlotPtMax = 200;
// Fit range: 10 to 200 GeV/c
const double FitPtMin = 10.0;
const double FitPtMax = 200.0;
TString JetPtGenTitleX = "#it{p}_{T, jet}^{true} (GeV/#it{c})";
TString JetPtMCFinalTitleY = "d#it{N}/d#it{p}_{T}";
Double_t ptbinGen[27] = {0, 1, 2, 3, 4, 5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200, 300};


void CrossSectionEfficiencyINELsel8jetSpectraCharged() {
  // File path
  const char* fileName = "~/cernbox/workspace/O2Physics/jets/AnalysisResults/560874_AnalysisResults.root";
  const char* datasetName = "LHC25a2b";
  const char* histPath2D = "jet-spectra-charged/h2_jet_pt_part_eventselection";
  const char* histPathMcCollisions = "jet-spectra-charged/h_mccollisions";
  const char* histPathMcCollisionsZvtx = "jet-spectra-charged/h_mccollisions_zvertex";
  
  // Open file
  TFile *file = TFile::Open(fileName, "read");
  
  if (!file || file->IsZombie()) {
    std::cout << "Error: Cannot open file " << fileName << std::endl;
    return;
  }

  // Get 2D histogram
  TH2 *h2EventSelection = (TH2*)file->Get(histPath2D);
  if (!h2EventSelection) {
    std::cout << "Error: Cannot find histogram " << histPath2D << " in file " << fileName << std::endl;
    file->Close();
    return;
  }

  // Get h_mccollisions histogram
  TH1 *hMcCollisions = (TH1*)file->Get(histPathMcCollisions);
  if (!hMcCollisions) {
    std::cout << "Error: Cannot find histogram " << histPathMcCollisions << " in file " << fileName << std::endl;
    file->Close();
    return;
  }

  // Get h_mccollisions_zvertex for zvtx cut calculation
  TH1 *hMcCollisionsZvtx = (TH1*)file->Get(histPathMcCollisionsZvtx);
  if (!hMcCollisionsZvtx) {
    std::cout << "Error: Cannot find histogram " << histPathMcCollisionsZvtx << " in file " << fileName << std::endl;
    file->Close();
    return;
  }

  // Get event counts from h_mccollisions
  // Index mapping: y index 1->0.5, 3->1.5, 4->2.5, 5->3.5, 6->4.5, 7->5.5
  Double_t NmccollINEL = hMcCollisions->GetBinContent(hMcCollisions->FindBin(0.5));  // allMcColl
  Double_t NmccollNoRecoColl = hMcCollisions->GetBinContent(hMcCollisions->FindBin(1.5));
  Double_t NmccollSplitColl = hMcCollisions->GetBinContent(hMcCollisions->FindBin(2.5));
  Double_t NmccollRecoEvtSel = hMcCollisions->GetBinContent(hMcCollisions->FindBin(3.5));
  Double_t NmccollCentralityCut = hMcCollisions->GetBinContent(hMcCollisions->FindBin(4.5));
  Double_t NmccollOccupancyCut = hMcCollisions->GetBinContent(hMcCollisions->FindBin(5.5));  // This is SelMC

  // Calculate zvtx cut efficiency from h_mccollisions_zvertex
  // Gaussian fit with mean fixed to raw distribution mean
  double rawMean = hMcCollisionsZvtx->GetMean();
  double rawStd = hMcCollisionsZvtx->GetStdDev();
  
  // Create Gaussian fit function
  TF1 *gausFit = new TF1("gausFit", "gaus", -15, 15);
  gausFit->SetParameter(0, hMcCollisionsZvtx->GetMaximum());
  gausFit->SetParameter(1, rawMean);  // Set mean to raw distribution mean
  gausFit->SetParameter(2, rawStd);
  gausFit->FixParameter(1, rawMean);  // Fix the mean parameter
  hMcCollisionsZvtx->Fit(gausFit, "Q0", "", -15, 15);
  
  // Calculate integral ratio: [-10, +10] / total
  // This gives the efficiency of passing the zvtx cut
  double totalIntegral = gausFit->Integral(-15, 15);
  double zvtxRangeIntegral = gausFit->Integral(-10, 10);
  double zvtxEfficiency = (totalIntegral > 0) ? zvtxRangeIntegral / totalIntegral : 0.96;
  
  std::cout << "Zvtx cut calculation:" << std::endl;
  std::cout << "  Raw mean: " << rawMean << ", std: " << rawStd << std::endl;
  std::cout << "  Total integral: " << totalIntegral << std::endl;
  std::cout << "  [-10, +10] integral: " << zvtxRangeIntegral << std::endl;
  std::cout << "  Efficiency: " << zvtxEfficiency << std::endl;
  
  Double_t NmccollZvtx = NmccollINEL * zvtxEfficiency;
  
  std::cout << "Event counts:" << std::endl;
  std::cout << "  NmccollINEL: " << NmccollINEL << std::endl;
  std::cout << "  NmccollZvtx: " << NmccollZvtx << " (efficiency: " << zvtxEfficiency << ")" << std::endl;
  std::cout << "  NmccollNoRecoColl: " << NmccollNoRecoColl << std::endl;
  std::cout << "  NmccollSplitColl: " << NmccollSplitColl << std::endl;
  std::cout << "  NmccollRecoEvtSel: " << NmccollRecoEvtSel << std::endl;
  std::cout << "  NmccollCentralityCut: " << NmccollCentralityCut << std::endl;
  std::cout << "  NmccollOccupancyCut (SelMC): " << NmccollOccupancyCut << std::endl;

  // Project 2D histogram to get 1D histograms for each step
  // y index: 1=INEL, 2=zvtx, 3=noRecoColl, 4=splitColl, 5=recoEvtSel, 6=centralitycut, 7=occupancycut
  TH1 *hInel = h2EventSelection->ProjectionX("hInel", 1, 1);  // y index 1
  TH1 *hZvtx = h2EventSelection->ProjectionX("hZvtx", 2, 2);  // y index 2
  TH1 *hNoRecoColl = h2EventSelection->ProjectionX("hNoRecoColl", 3, 3);  // y index 3
  TH1 *hSplitColl = h2EventSelection->ProjectionX("hSplitColl", 4, 4);  // y index 4
  TH1 *hRecoEvtSel = h2EventSelection->ProjectionX("hRecoEvtSel", 5, 5);  // y index 5
  TH1 *hCentralityCut = h2EventSelection->ProjectionX("hCentralityCut", 6, 6);  // y index 6
  TH1 *hSelMC = h2EventSelection->ProjectionX("hSelMC", 7, 7);  // y index 7 (occupancycut = SelMC)
  
  if (!hInel || !hZvtx || !hNoRecoColl || !hSplitColl || !hRecoEvtSel || !hCentralityCut || !hSelMC) {
    std::cout << "Error: Failed to project histograms from 2D histogram" << std::endl;
    file->Close();
    return;
  }
  
  // Rebin histograms to ptbinGen (26 bins)
  TH1* hInelRebinned = (TH1*)hInel->Rebin(26, "hInel_rebinned", ptbinGen);
  TH1* hZvtxRebinned = (TH1*)hZvtx->Rebin(26, "hZvtx_rebinned", ptbinGen);
  TH1* hNoRecoCollRebinned = (TH1*)hNoRecoColl->Rebin(26, "hNoRecoColl_rebinned", ptbinGen);
  TH1* hSplitCollRebinned = (TH1*)hSplitColl->Rebin(26, "hSplitColl_rebinned", ptbinGen);
  TH1* hRecoEvtSelRebinned = (TH1*)hRecoEvtSel->Rebin(26, "hRecoEvtSel_rebinned", ptbinGen);
  TH1* hCentralityCutRebinned = (TH1*)hCentralityCut->Rebin(26, "hCentralityCut_rebinned", ptbinGen);
  TH1* hSelMCRebinned = (TH1*)hSelMC->Rebin(26, "hSelMC_rebinned", ptbinGen);
  
  if (!hInelRebinned || !hZvtxRebinned || !hNoRecoCollRebinned || !hSplitCollRebinned || !hRecoEvtSelRebinned || !hCentralityCutRebinned || !hSelMCRebinned) {
    std::cout << "Error: Rebin failed" << std::endl;
    file->Close();
    return;
  }
  
  // Replace original histograms with rebinned ones
  delete hInel;
  delete hZvtx;
  delete hNoRecoColl;
  delete hSplitColl;
  delete hRecoEvtSel;
  delete hCentralityCut;
  delete hSelMC;
  
  hInel = hInelRebinned;
  hZvtx = hZvtxRebinned;
  hNoRecoColl = hNoRecoCollRebinned;
  hSplitColl = hSplitCollRebinned;
  hRecoEvtSel = hRecoEvtSelRebinned;
  hCentralityCut = hCentralityCutRebinned;
  hSelMC = hSelMCRebinned;
  
  // Set directory to 0 to avoid issues
  hInel->SetDirectory(0);
  hZvtx->SetDirectory(0);
  hNoRecoColl->SetDirectory(0);
  hSplitColl->SetDirectory(0);
  hRecoEvtSel->SetDirectory(0);
  hCentralityCut->SetDirectory(0);
  hSelMC->SetDirectory(0);
  
  // Calculate jet counts in 10-200 GeV range BEFORE normalization
  int binLow = hInel->GetXaxis()->FindBin(10.0);
  int binHigh = hInel->GetXaxis()->FindBin(200.0);
  
  double nJetsInel = hInel->Integral(binLow, binHigh);
  double nJetsZvtx = hZvtx->Integral(binLow, binHigh);
  double nJetsNoRecoColl = hNoRecoColl->Integral(binLow, binHigh);
  double nJetsSplitColl = hSplitColl->Integral(binLow, binHigh);
  double nJetsRecoEvtSel = hRecoEvtSel->Integral(binLow, binHigh);
  double nJetsCentralityCut = hCentralityCut->Integral(binLow, binHigh);
  double nJetsSelMC = hSelMC->Integral(binLow, binHigh);
  
  // Clone histograms for raw jet count version (before normalization)
  // Only keep: Zvtx10, hasColl, selMC (will be labeled as INEL, Coll, TVX)
  TH1 *hZvtxRaw = (TH1*)hZvtx->Clone("hZvtxRaw");      // Will be labeled as "INEL"
  TH1 *hNoRecoCollRaw = (TH1*)hNoRecoColl->Clone("hNoRecoCollRaw");  // Will be labeled as "Coll"
  TH1 *hSelMCRaw = (TH1*)hSelMC->Clone("hSelMCRaw");   // Will be labeled as "TVX"
  
  // Create Filipad2 for event normalized plot with ratio
  static int nn = 0;
  Filipad2 *MCPPad = new Filipad2("MCP_Comparison_JetSpectraCharged", ++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
  MCPPad->Draw();
  
  // Get pads
  TPad *mainPad = MCPPad->GetPad(1);  // Upper pad for event normalized plot
  TPad *ratioPad = MCPPad->GetPad(2); // Lower pad for ratio
  
  // Set pad options (gridx, gridy, logx, logy)
  optFili(*mainPad, 0, 0, 0, 1);   // no grid, logy = 1
  optFili(*ratioPad, 0, 0, 0, 0);  // no grid, logy = 0
  
  // Draw main plot
  mainPad->cd();
  
  // Setup first histogram (reference) - normalize by width
  // Use better color scheme for visibility
  hset(*hInel, JetPtGenTitleX, "1/N_{evt} dN_{ch-jet}^{true}/dp_{T}", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
  hoptset(*hInel, 0, kRed, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 22);
  hInel->Scale(1./NmccollINEL, "width");
  hInel->SetLineStyle(1);
  hInel->SetLineWidth(2);
  hInel->Draw("pe");

  hoptset(*hZvtx, 0, kBlue, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 21);
  hZvtx->Scale(1./NmccollZvtx, "width");
  hZvtx->SetLineStyle(1);
  hZvtx->SetLineWidth(2);
  hZvtx->Draw("pesame");

  hoptset(*hNoRecoColl, 0, kGreen+2, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 23);
  hNoRecoColl->Scale(1./NmccollNoRecoColl, "width");
  hNoRecoColl->SetLineStyle(1);
  hNoRecoColl->SetLineWidth(2);
  hNoRecoColl->Draw("pesame");

  hoptset(*hSplitColl, 0, kMagenta+2, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 24);
  hSplitColl->Scale(1./NmccollSplitColl, "width");
  hSplitColl->SetLineStyle(1);
  hSplitColl->SetLineWidth(2);
  hSplitColl->Draw("pesame");

  hoptset(*hRecoEvtSel, 0, kBlack, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 20);
  hRecoEvtSel->Scale(1./NmccollRecoEvtSel, "width");
  hRecoEvtSel->SetLineStyle(1);
  hRecoEvtSel->SetLineWidth(3);
  hRecoEvtSel->Draw("pesame");
  
  // Create legend with integrated jet and event counts
  TLegend *leg = new TLegend(0.40, 0.40, 0.65, 0.95, NULL, "brNDC");
  leg->SetTextSize(0.045);
  leg->SetBorderSize(0);
  leg->AddEntry("", "MC particle-level", "");
  leg->AddEntry("", "(N_{jets} [10-200 GeV], N_{evt})", "");
  leg->AddEntry(hInel, Form("INEL (%.2e, %.1e)", nJetsInel, NmccollINEL), "pe");
  leg->AddEntry(hZvtx, Form("+Zvtx10 (%.2e, %.1e)", nJetsZvtx, NmccollZvtx), "pe");
  leg->AddEntry(hNoRecoColl, Form("+hasColl (%.2e, %.1e)", nJetsNoRecoColl, NmccollNoRecoColl), "pe");
  leg->AddEntry(hSplitColl, Form("+noSplit (%.2e, %.1e)", nJetsSplitColl, NmccollSplitColl), "pe");
  leg->AddEntry(hRecoEvtSel, Form("+selMC (%.2e, %.1e)", nJetsRecoEvtSel, NmccollRecoEvtSel), "pe");
  leg->Draw();
  
  // Update main pad to ensure it's visible in ROOT window
  mainPad->Update();
  
  // Draw ratio plot for event normalized
  ratioPad->cd();
  
  // Arrays for looping - all steps vs selMC (recoEvtSel)
  // Use same color scheme as main plot for consistency
  TH1* hists[] = {hInel, hZvtx, hNoRecoColl, hSplitColl};
  Color_t colors[] = {kRed, kBlue, kGreen+2, kMagenta+2};
  const char* labels[] = {"INEL", "Zvtx", "hasColl", "noSplit"};
  int nHists = 4;

  TString ratioName = "MCP_Ratio_Steps_vs_selMC";
  for (int i = 0; i < nHists; i++) {
    TString ratioName_i = Form("%s_%s", ratioName.Data(), labels[i]);
    TH1 *ratioHist = DrawRatioTH1(hists[i], hRecoEvtSel);
    ratioHist->SetName(ratioName_i.Data());
    ratioHist->SetDirectory(0);
    
    // Setup histogram properties
    if (i == 0) {
      // First histogram: setup axis and draw
      hset(*ratioHist, JetPtGenTitleX, "Ref. / (selMC)", 1.2, 1.0, 0.07, 0.06, 0.01, 0.01,
           0.07, 0.07, 510, 505);
      ratioHist->GetXaxis()->SetRangeUser(PlotPtMin, PlotPtMax);
      ratioHist->GetYaxis()->SetRangeUser(0.7, 1.3);
    }
    
    ratioHist->SetMarkerColor(colors[i]);
    ratioHist->SetLineColor(colors[i]);
    ratioHist->SetMarkerSize(0.8);
    ratioHist->SetLineStyle(1);
    ratioHist->SetLineWidth(2);
    ratioHist->SetFillColor(0);
    
    // Draw: first with "e", others with "esame"
    if (i == 0) {
      ratioHist->Draw("e");
      // Draw horizontal line at y=1
      TLine *line1 = new TLine(PlotPtMin, 1.0, PlotPtMax, 1.0);
      line1->SetLineColor(kBlack);
      line1->SetLineStyle(1);
      line1->SetLineWidth(1);
      line1->Draw("same");
      
      // Draw infoLatex header after first histogram is drawn
      TLatex *infoLatex = new TLatex();
      infoLatex->SetNDC();
      infoLatex->SetTextSize(0.06);
      infoLatex->SetTextColor(kBlack);
      infoLatex->SetTextAlign(12);  // Left and bottom aligned
      infoLatex->DrawLatex(0.18, 0.95, "Constant fit (10-200 GeV):");
    } else {
      ratioHist->Draw("esame");
    }
    
    TF1 *fConst = new TF1(Form("fConst_%d", i), "[0]", FitPtMin, FitPtMax);
    fConst->SetParameter(0, 1.0);  // Initial value
    fConst->SetLineColor(colors[i]);
    fConst->SetLineStyle(2);  // Dashed line
    fConst->SetLineWidth(3);
    ratioHist->Fit(fConst, "Q0", "", FitPtMin, FitPtMax);
    // Draw the fit function
    fConst->Draw("same");
    
    // Display fit result
    double fitValue = fConst->GetParameter(0);
    double fitError = fConst->GetParError(0);
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.07);
    latex->SetTextColor(colors[i]);
    latex->SetTextAlign(12);  // Left and bottom aligned
    double yPos = 0.95 - i * 0.08;  // Stack text vertically
    latex->DrawLatex(0.55, yPos, Form("%s: %.4f #pm %.4f", labels[i], fitValue, fitError));
  }
  
  // Update canvas to ensure everything is visible in ROOT window
  MCPPad->C->Update();

  // Print first canvas (event normalized + ratio)
  MCPPad->C->Print(Form("plots/MCP_jet_part_spectra_steps_vs_SelMC_%s.pdf", datasetName));
  
  // ========== Second plot: Raw jet counts (not normalized) ==========
  // Create Filipad2 for raw jet count plot with ratio
  static int nn2 = 0;
  Filipad2 *MCPPadRaw = new Filipad2("MCP_Comparison_JetCountsCharged", ++nn2, 2, 0.4, 100, 50, 0.7, 1, 1);
  MCPPadRaw->Draw();
  
  // Get pads
  TPad *mainPadRaw = MCPPadRaw->GetPad(1);  // Upper pad for raw jet count plot
  TPad *ratioPadRaw = MCPPadRaw->GetPad(2); // Lower pad for ratio
  
  // Set pad options (gridx, gridy, logx, logy)
  optFili(*mainPadRaw, 0, 0, 0, 1);   // no grid, logy = 1
  optFili(*ratioPadRaw, 0, 0, 0, 0);  // no grid, logy = 0
  
  // Draw main plot (raw jet counts, no normalization)
  // Only plot: Zvtx10 (labeled as INEL), hasColl (labeled as Coll), selMC (labeled as selMC)
  mainPadRaw->cd();
  
  // First histogram: Zvtx10 (labeled as "INEL")
  hset(*hZvtxRaw, JetPtGenTitleX, "dN_{ch-jet}^{true}/dp_{T}", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
  hoptset(*hZvtxRaw, 0, kRed, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 22);
  hZvtxRaw->Scale(1.0, "width");  // Only scale by width, not by events
  hZvtxRaw->SetLineStyle(1);
  hZvtxRaw->SetLineWidth(2);
  hZvtxRaw->Draw("pe");

  // Second histogram: hasColl (labeled as "Coll")
  hoptset(*hNoRecoCollRaw, 0, kBlue, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 21);
  hNoRecoCollRaw->Scale(1.0, "width");
  hNoRecoCollRaw->SetLineStyle(1);
  hNoRecoCollRaw->SetLineWidth(2);
  hNoRecoCollRaw->Draw("pesame");

  // Third histogram: selMC (labeled as "selMC")
  hoptset(*hSelMCRaw, 0, kGreen+2, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 23);
  hSelMCRaw->Scale(1.0, "width");
  hSelMCRaw->SetLineStyle(1);
  hSelMCRaw->SetLineWidth(2);
  hSelMCRaw->Draw("pesame");
  
  // Create legend with jet counts only (no event counts)
  TLegend *legRaw = new TLegend(0.40, 0.40, 0.65, 0.95, NULL, "brNDC");
  legRaw->SetTextSize(0.045);
  legRaw->SetBorderSize(0);
  legRaw->AddEntry("", "MC particle-level", "");
  legRaw->AddEntry("", "N_{jets} [10-200 GeV]", "");
  legRaw->AddEntry(hZvtxRaw, Form("INEL (%.2e)", nJetsZvtx), "pe");
  legRaw->AddEntry(hNoRecoCollRaw, Form("Coll (%.2e)", nJetsNoRecoColl), "pe");
  legRaw->AddEntry(hSelMCRaw, Form("selMC (%.2e)", nJetsSelMC), "pe");
  legRaw->Draw();
  
  // Update main pad to ensure it's visible in ROOT window
  mainPadRaw->Update();
  
  // Draw ratio plot for raw jet counts
  // Ratio: Coll and selMC vs INEL (reference)
  ratioPadRaw->cd();
  
  // Arrays for looping - Coll and selMC vs INEL
  TH1* histsRaw[] = {hNoRecoCollRaw, hSelMCRaw};
  Color_t colorsRaw[] = {kBlue, kGreen+2};
  const char* labelsRaw[] = {"Coll", "selMC"};
  int nHistsRaw = 2;

  TString ratioNameRaw = "MCP_Ratio_Coll_selMC_vs_INEL_Raw";
  for (int i = 0; i < nHistsRaw; i++) {
    TString ratioName_i = Form("%s_%s", ratioNameRaw.Data(), labelsRaw[i]);
    TH1 *ratioHist = DrawRatioTH1(histsRaw[i], hZvtxRaw);
    ratioHist->SetName(ratioName_i.Data());
    ratioHist->SetDirectory(0);
    
    // Setup histogram properties
    if (i == 0) {
      // First histogram: setup axis and draw
      hset(*ratioHist, JetPtGenTitleX, "Ref. / (INEL)", 1.2, 1.0, 0.07, 0.06, 0.01, 0.01,
           0.07, 0.07, 510, 505);
      ratioHist->GetXaxis()->SetRangeUser(PlotPtMin, PlotPtMax);
      ratioHist->GetYaxis()->SetRangeUser(0.8, 1.1);
    }
    
    ratioHist->SetMarkerColor(colorsRaw[i]);
    ratioHist->SetLineColor(colorsRaw[i]);
    ratioHist->SetMarkerSize(0.8);
    ratioHist->SetLineStyle(1);
    ratioHist->SetLineWidth(2);
    ratioHist->SetFillColor(0);
    
    // Draw: first with "e", others with "esame"
    if (i == 0) {
      ratioHist->Draw("e");
      // Draw horizontal line at y=1
      TLine *line1 = new TLine(PlotPtMin, 1.0, PlotPtMax, 1.0);
      line1->SetLineColor(kBlack);
      line1->SetLineStyle(1);
      line1->SetLineWidth(1);
      line1->Draw("same");
      
      // Draw infoLatex header after first histogram is drawn
      TLatex *infoLatex = new TLatex();
      infoLatex->SetNDC();
      infoLatex->SetTextSize(0.07);
      infoLatex->SetTextColor(kBlack);
      infoLatex->SetTextAlign(12);  // Left and bottom aligned
      infoLatex->DrawLatex(0.25, 0.95, "10-200 GeV fit:");
    } else {
      ratioHist->Draw("esame");
    }
    
    TF1 *fConst = new TF1(Form("fConstRaw_%d", i), "[0]", FitPtMin, FitPtMax);
    fConst->SetParameter(0, 1.0);  // Initial value
    fConst->SetLineColor(colorsRaw[i]);
    fConst->SetLineStyle(2);  // Dashed line
    fConst->SetLineWidth(3);
    ratioHist->Fit(fConst, "Q0", "", FitPtMin, FitPtMax);
    // Draw the fit function
    fConst->Draw("same");
    
    // Display fit result
    double fitValue = fConst->GetParameter(0);
    double fitError = fConst->GetParError(0);
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.07);
    latex->SetTextColor(colorsRaw[i]);
    latex->SetTextAlign(12);  // Left and bottom aligned
    double yPos = 0.95 - i * 0.08;  // Stack text vertically
    latex->DrawLatex(0.55, yPos, Form("%s: %.4f #pm %.4f", labelsRaw[i], fitValue, fitError));
  }
  
  // Update canvas to ensure everything is visible in ROOT window
  MCPPadRaw->C->Update();

  // Print second canvas (raw jet counts + ratio)
  MCPPadRaw->C->Print("plots/MCP_jet_part_counts_steps_vs_SelMC.pdf");
  
  // Cleanup
  // delete hInelRaw;
  // delete hZvtxRaw;
  // delete hNoRecoCollRaw;
  // delete hSplitCollRaw;
  // delete hRecoEvtSelRaw;
  // file->Close();
}

void CrossSectionEfficiencyINELsel8jetFinderQA() {
  // File paths
  const char* file1 = "../../../jets/AnalysisResults/496218_AnalysisResults.root";
  const char* file2 = "../../../jets/AnalysisResults/498287_AnalysisResults.root";
  const char* histPath = "jet-finder-charged-qa/h_jet_pt_part";
  
  // Open files
  TFile *f1 = TFile::Open(file1, "read");
  TFile *f2 = TFile::Open(file2, "read");
  
  if (!f1 || f1->IsZombie()) {
    std::cout << "Error: Cannot open file " << file1 << std::endl;
    return;
  }
  if (!f2 || f2->IsZombie()) {
    std::cout << "Error: Cannot open file " << file2 << std::endl;
    if (f1) f1->Close();
    return;
  }
  
  // Get histograms
  TH1 *h1 = (TH1*)f1->Get(histPath);
  TH1 *h2 = (TH1*)f2->Get(histPath);
  
  if (!h1) {
    std::cout << "Error: Cannot find histogram " << histPath << " in file " << file1 << std::endl;
    f1->Close();
    f2->Close();
    return;
  }
  if (!h2) {
    std::cout << "Error: Cannot find histogram " << histPath << " in file " << file2 << std::endl;
    f1->Close();
    f2->Close();
    return;
  }
  
  // Debug: Print original histogram info
  std::cout << "Original h1 entries: " << h1->GetEntries() << ", integral: " << h1->Integral() << std::endl;
  std::cout << "Original h1 range: [" << h1->GetXaxis()->GetXmin() << ", " << h1->GetXaxis()->GetXmax() << "]" << std::endl;
  std::cout << "Original h2 entries: " << h2->GetEntries() << ", integral: " << h2->Integral() << std::endl;
  std::cout << "Original h2 range: [" << h2->GetXaxis()->GetXmin() << ", " << h2->GetXaxis()->GetXmax() << "]" << std::endl;
  
  // Clone histograms to avoid modifying originals
  // ptbinGen has 27 elements = 26 bins, so first parameter should be 26
  TH1 *h1Clone = (TH1*)h1->Clone("h1_MCP");
  h1Clone = (TH1*)h1Clone->Rebin(26, "h1_MCP_rebinned", ptbinGen);
  if (!h1Clone) {
    std::cout << "Error: Rebin failed for h1" << std::endl;
    f1->Close();
    f2->Close();
    return;
  }
  
  TH1 *h2Clone = (TH1*)h2->Clone("h2_MCP");
  h2Clone = (TH1*)h2Clone->Rebin(26, "h2_MCP_rebinned", ptbinGen);
  if (!h2Clone) {
    std::cout << "Error: Rebin failed for h2" << std::endl;
    f1->Close();
    f2->Close();
    return;
  }
  
  // Debug: Print rebinned histogram info
  std::cout << "Rebinned h1 entries: " << h1Clone->GetEntries() << ", integral: " << h1Clone->Integral() << std::endl;
  std::cout << "Rebinned h1 range: [" << h1Clone->GetXaxis()->GetXmin() << ", " << h1Clone->GetXaxis()->GetXmax() << "]" << std::endl;
  std::cout << "Rebinned h2 entries: " << h2Clone->GetEntries() << ", integral: " << h2Clone->Integral() << std::endl;
  std::cout << "Rebinned h2 range: [" << h2Clone->GetXaxis()->GetXmin() << ", " << h2Clone->GetXaxis()->GetXmax() << "]" << std::endl;
  
  // Set directory to 0 to avoid issues
  h1Clone->SetDirectory(0);
  h2Clone->SetDirectory(0);
  
  // Create Filipad2 for ratio plot
  static int nn = 0;
  Filipad2 *MCPPad = new Filipad2("MCP_Comparison", ++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
  MCPPad->Draw();
  
  // Get pads
  TPad *mainPad = MCPPad->GetPad(1);  // Upper pad for main plot
  TPad *ratioPad = MCPPad->GetPad(2); // Lower pad for ratio
  
  // Set pad options (gridx, gridy, logx, logy)
  optFili(*mainPad, 1, 1, 0, 1);   // logy = 1
  optFili(*ratioPad, 1, 1, 0, 0);  // logy = 0
  
  // Draw main plot
  mainPad->cd();
  
  // Setup first histogram (reference) - normalize by width
  // 검정 네모 (black square)
  hset(*h1Clone, JetPtGenTitleX, JetPtMCFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
  hoptset(*h1Clone, 0, kBlack, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 22);
  h1Clone->Scale(1., "width");
  h1Clone->SetLineStyle(1);  // Solid line
  h1Clone->SetLineWidth(1);
  h1Clone->Draw("pe");
  
  // Setup second histogram - normalize by width
  // 빨간 동그라미 (red circle)
  hoptset(*h2Clone, 0, kRed, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 20);
  h2Clone->Scale(1., "width");
  h2Clone->SetLineStyle(1);  // Solid line
  h2Clone->SetLineWidth(1);
  h2Clone->Draw("pesame");
  
  // Create legend
  TLegend *leg = new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->AddEntry(h1Clone, "MCP selMC+zvtx", "pe");
  leg->AddEntry(h2Clone, "MCP INEL", "pe");
  leg->Draw();
  
  // Draw ratio plot
  ratioPad->cd();
  
  TString ratioName = "MCP_Ratio_496218_498287";
  TH1 *ratioHist = DrawRatio(ratioName.Data(), h1Clone, h2Clone, 
                             JetPtGenTitleX, "MCP (INEL / selMC+zvtx)", 
                             kRed, 0.5, 1.5, 0.75);  
  // Draw horizontal line at y=1
  TLine *line1 = new TLine(PlotPtMin, 1.0, PlotPtMax, 1.0);
  line1->SetLineColor(kBlack);
  line1->SetLineStyle(1);
  line1->SetLineWidth(1);
  line1->Draw("same");
  
  // Fit constant function to ratio histogram
  TF1 *fConst = new TF1("fConst", "[0]", FitPtMin, FitPtMax);
  fConst->SetParameter(0, 1.0);  // Initial value
  fConst->SetLineColor(kBlack);
  fConst->SetLineStyle(1);  // Solid line (실선)
  fConst->SetLineWidth(3);
  
  // Perform fit in range 10-140 GeV/c
  ratioHist->Fit(fConst, "Q0", "", FitPtMin, FitPtMax);
  
  // Draw the fit function (검정 실선으로)
  fConst->Draw("same");
  
  // Display fit result
  double fitValue = fConst->GetParameter(0);
  double fitError = fConst->GetParError(0);
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.06);
  latex->SetTextColor(kBlack);
  latex->SetTextAlign(12);  // Left and bottom aligned
  latex->DrawLatex(0.2, 0.25, Form("Constant fit: %.4f #pm %.4f", fitValue, fitError));
  
  // Debug: Print ratio histogram info
  std::cout << "Ratio histogram entries: " << ratioHist->GetEntries() << ", non-zero bins: ";
  int nNonZero = 0;
  for (int i = 1; i <= ratioHist->GetNbinsX(); i++) {
    if (ratioHist->GetBinContent(i) > 0) nNonZero++;
  }
  std::cout << nNonZero << std::endl;
  std::cout << "Constant fit value: " << fitValue << " +/- " << fitError << std::endl;
  
  // Print canvas
  MCPPad->C->Print("plots/MCP_jet_part_496218_vs_498287.pdf");
  
  // Cleanup
  f1->Close();
  f2->Close();
  
  std::cout << "MCP comparison plot saved as MCP_Comparison_496218_vs_498287.pdf" << std::endl;
}

void CrossSectionEfficiencySelMCSel8() {
  // File paths
  const char* fileSelMC = "~/cernbox/workspace/O2Physics/jets/AnalysisResults/553776_AnalysisResults.root";
  const char* fileSel8 = "~/cernbox/workspace/O2Physics/jets/AnalysisResults/496287_AnalysisResults.root";
  const char* histPath = "jet-spectra-charged/h_jet_pt";
  
  // Open files
  TFile *fSelMC = TFile::Open(fileSelMC, "read");
  TFile *fSel8 = TFile::Open(fileSel8, "read");
  
  if (!fSelMC || fSelMC->IsZombie()) {
    std::cout << "Error: Cannot open file " << fileSelMC << std::endl;
    return;
  }
  if (!fSel8 || fSel8->IsZombie()) {
    std::cout << "Error: Cannot open file " << fileSel8 << std::endl;
    if (fSelMC) fSelMC->Close();
    return;
  }
  
  // Get histograms
  TH1 *hSelMC = (TH1*)fSelMC->Get(histPath);
  TH1 *hSel8 = (TH1*)fSel8->Get(histPath);
  
  if (!hSelMC) {
    std::cout << "Error: Cannot find histogram " << histPath << " in file " << fileSelMC << std::endl;
    fSelMC->Close();
    fSel8->Close();
    return;
  }
  if (!hSel8) {
    std::cout << "Error: Cannot find histogram " << histPath << " in file " << fileSel8 << std::endl;
    fSelMC->Close();
    fSel8->Close();
    return;
  }
  
  // Clone and rebin histograms
  TH1 *hSelMCClone = (TH1*)hSelMC->Clone("hSelMC_CrossSection");
  hSelMCClone = (TH1*)hSelMCClone->Rebin(26, "hSelMC_CrossSection_rebinned", ptbinGen);
  
  TH1 *hSel8Clone = (TH1*)hSel8->Clone("hSel8_CrossSection");
  hSel8Clone = (TH1*)hSel8Clone->Rebin(26, "hSel8_CrossSection_rebinned", ptbinGen);
  
  if (!hSelMCClone || !hSel8Clone) {
    std::cout << "Error: Rebin failed" << std::endl;
    fSelMC->Close();
    fSel8->Close();
    return;
  }
  
  // Set directory to 0
  hSelMCClone->SetDirectory(0);
  hSel8Clone->SetDirectory(0);
  
  // Calculate jet counts in 10-140 GeV range (no normalization)
  int binLow = hSelMCClone->GetXaxis()->FindBin(10.0);
  int binHigh = hSelMCClone->GetXaxis()->FindBin(140.0);
  
  double nJetsSelMC = hSelMCClone->Integral(binLow, binHigh);
  double nJetsSel8 = hSel8Clone->Integral(binLow, binHigh);
  double ratio = (nJetsSel8 > 0) ? nJetsSelMC / nJetsSel8 : 0.0;
  
  std::cout << "Jet counts (10-140 GeV):" << std::endl;
  std::cout << "  selMC: " << nJetsSelMC << std::endl;
  std::cout << "  sel8: " << nJetsSel8 << std::endl;
  std::cout << "  Ratio (selMC/sel8): " << ratio << std::endl;
  
  // Create Filipad2 for ratio plot
  static int nn = 0;
  Filipad2 *CrossSectionPad = new Filipad2("CrossSection_SelMC_vs_Sel8", ++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
  CrossSectionPad->Draw();
  
  // Get pads
  TPad *mainPad = CrossSectionPad->GetPad(1);  // Upper pad for main plot
  TPad *ratioPad = CrossSectionPad->GetPad(2); // Lower pad for ratio
  
  // Set pad options
  optFili(*mainPad, 0, 0, 0, 1);   // no grid, logy = 1
  optFili(*ratioPad, 0, 0, 0, 0);  // no grid, logy = 0
  
  // Draw main plot (raw jet counts, no normalization)
  mainPad->cd();
  
  hset(*hSelMCClone, "#it{p}_{T, jet}^{raw} [GeV/c]", "dN_{jet}^{raw}/dp_{T}", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
  hoptset(*hSelMCClone, 0, kRed, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 20);
  hSelMCClone->Scale(1.0, "width");  // Only scale by width, not by events
  hSelMCClone->SetLineStyle(1);
  hSelMCClone->SetLineWidth(2);
  hSelMCClone->Draw("pe");
  
  hoptset(*hSel8Clone, 0, kBlue, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 21);
  hSel8Clone->Scale(1.0, "width");
  hSel8Clone->SetLineStyle(1);
  hSel8Clone->SetLineWidth(2);
  hSel8Clone->Draw("pesame");
  
  // Create legend
  TLegend *leg = new TLegend(0.40, 0.60, 0.65, 0.95, NULL, "brNDC");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->AddEntry("", "Raw data (N_{jets} [10-140 GeV])", "");
  leg->AddEntry(hSelMCClone, Form("selMC (%.2e)", nJetsSelMC), "pe");
  leg->AddEntry(hSel8Clone, Form("sel8 (%.2e)", nJetsSel8), "pe");
  leg->Draw();
  
  // Update main pad
  mainPad->Update();
  
  // Draw ratio plot (selMC / sel8)
  ratioPad->cd();
  
  TH1 *ratioHist = DrawRatioTH1(hSelMCClone, hSel8Clone);
  ratioHist->SetName("CrossSection_Ratio_SelMC_Sel8");
  ratioHist->SetDirectory(0);
  
  hset(*ratioHist, "#it{p}_{T, jet}^{raw} [GeV/c]", "selMC / sel8", 1.2, 1.0, 0.07, 0.06, 0.01, 0.01, 0.07, 0.07, 510, 505);
  ratioHist->GetXaxis()->SetRangeUser(PlotPtMin, PlotPtMax);
  ratioHist->GetYaxis()->SetRangeUser(0.8, 1.2);
  ratioHist->SetMarkerColor(kRed);
  ratioHist->SetLineColor(kRed);
  ratioHist->SetMarkerSize(0.8);
  ratioHist->SetLineStyle(1);
  ratioHist->SetLineWidth(2);
  ratioHist->SetFillColor(0);
  ratioHist->Draw("e");
  
  // Draw horizontal line at y=1
  TLine *line1 = new TLine(PlotPtMin, 1.0, PlotPtMax, 1.0);
  line1->SetLineColor(kBlack);
  line1->SetLineStyle(1);
  line1->SetLineWidth(1);
  line1->Draw("same");
  
  // Fit constant function in 10-140 GeV range
  TF1 *fConst = new TF1("fConst_CrossSection", "[0]", 10.0, 140.0);
  fConst->SetParameter(0, 1.0);
  fConst->SetLineColor(kRed);
  fConst->SetLineStyle(2);  // Dashed line
  fConst->SetLineWidth(3);
  ratioHist->Fit(fConst, "Q0", "", 10.0, 140.0);
  fConst->Draw("same");
  
  // Display fit result
  double fitValue = fConst->GetParameter(0);
  double fitError = fConst->GetParError(0);
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.07);
  latex->SetTextColor(kRed);
  latex->SetTextAlign(12);
  latex->DrawLatex(0.2, 0.25, Form("Constant fit (10-140 GeV): %.4f #pm %.4f", fitValue, fitError));
  
  // Update canvas
  CrossSectionPad->C->Update();
  
  // Print canvas
  CrossSectionPad->C->Print("plots/CrossSection_SelMC_vs_Sel8.pdf");
  
  // Cleanup
  fSelMC->Close();
  fSel8->Close();
  
  std::cout << "Cross section efficiency plot saved as CrossSection_SelMC_vs_Sel8.pdf" << std::endl;
}

//====================================================================
// New cross section / efficiency plots using updated jetSpectraCharged.cxx
// event-selection definition:
//  y-bin 1 -> INEL
//  y-bin 2 -> noRecoColl
//  y-bin 3 -> splitColl
//  y-bin 4 -> kTVX
//  y-bin 5 -> kTFBorder
//  y-bin 6 -> kITSROFBorder
//  y-bin 7 -> zvtx
//  y-bin 8 -> centrality cut
//
// Event counts are taken from h_mccollisions_eventselection with the same bin indices.
// This macro produces:
//  (1) INEL, noRecoColl, kTVX: raw jet-count spectra + ratios (selection / INEL)
//  (2) kTVX, kTFBorder, kITSROFBorder: raw jet-count spectra + ratios (selection / TVX)
//  (3) kTVX, kTFBorder, kITSROFBorder: invariant yields + ratios (selection / TVX)
//====================================================================
void CrossSectionEfficiencyINELCollNonSplitTVX()
{
  // Local pT range for this macro (0–200 GeV/c)
  const double PtMinLocal = 0.0;
  const double PtMaxLocal = 140.0;
  // File / histogram paths (MC particle-level jets from jetSpectraCharged.cxx)
  const char *fileName = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/569879_AnalysisResults.root"; // LHC25a2b w Mytuner new
  // const char *fileName = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/569879_AnalysisResults.root";
  // const char *fileName = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/570220_AnalysisResults.root"; // LHC25a2b w/o tuner old
  // const char *fileName = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/605134_AnalysisResults.root"; // LHC23k4h w Mytuner new
  // const char *fileName = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/615791_AnalysisResults.root"; // LHC25a2b w Mytuner new
  const char *datasetName = "LHC24f3c";
  // const char *datasetName = "LHC25a2b_new";
  // const char *datasetName = "LHC23k4h";
  const char *histPath2D = "jet-cross-section-efficiency/h2_jet_pt_part_eventselection";
  const char *histPathMcCollisions = "jet-cross-section-efficiency/h_mccollisions_eventselection";
  const char *histPathMcCollisionsWeight = "jet-cross-section-efficiency/h_mccollisions_eventselection_weighted";

  // Open file
  TFile *file = TFile::Open(fileName, "read");
  if (!file || file->IsZombie()) {
    std::cout << "Error: Cannot open file " << fileName << std::endl;
    return;
  }

  // Get 2D histogram (jet pT vs event selection)
  TH2 *h2EventSelection = (TH2 *)file->Get(histPath2D);
  if (!h2EventSelection) {
    std::cout << "Error: Cannot find histogram " << histPath2D << " in file " << fileName << std::endl;
    file->Close();
    return;
  }

  // Get mc-collision event counts
  TH1 *hMcCollisions = (TH1 *)file->Get(histPathMcCollisions);
  if (!hMcCollisions) {
    hMcCollisions = (TH1 *)file->Get(histPathMcCollisionsWeight);
    if (!hMcCollisions) {
      std::cout << "Error: Cannot find histogram " << histPathMcCollisions << " or " << histPathMcCollisionsWeight << " in file " << fileName << std::endl;
      file->Close();
      return;
    }
  }

  // Event counts per selection (new bin indices)
  Double_t NmccollINEL         = hMcCollisions->GetBinContent(1); // INEL
  Double_t NmccollNoRecoColl   = hMcCollisions->GetBinContent(2); // noRecoColl
  Double_t NmccollTVX          = hMcCollisions->GetBinContent(4); // kTVX
  Double_t NmccollTFBorder     = hMcCollisions->GetBinContent(5); // kTFBorder
  Double_t NmccollITSROFBorder = hMcCollisions->GetBinContent(6); // kITSROFBorder

  std::cout << "[CrossSectionEfficiencyINELCollNonSplitTVX] Event counts:" << std::endl;
  std::cout << "  INEL         (bin 1): " << NmccollINEL << std::endl;
  std::cout << "  noRecoColl   (bin 2): " << NmccollNoRecoColl << std::endl;
  std::cout << "  kTVX         (bin 4): " << NmccollTVX << std::endl;
  std::cout << "  kTFBorder    (bin 5): " << NmccollTFBorder << std::endl;
  std::cout << "  kITSROFBorder (bin 6): " << NmccollITSROFBorder << std::endl;

  // Project 2D histogram to pT for each selection
  // Plot 1: INEL, noRecoColl, kTVX
  TH1 *hINEL         = h2EventSelection->ProjectionX("hINEL",         1, 1);  // INEL
  TH1 *hNoRecoColl   = h2EventSelection->ProjectionX("hNoRecoColl",   2, 2);  // noRecoColl
  TH1 *hTVX          = h2EventSelection->ProjectionX("hTVX",          4, 4);  // kTVX
  // Plot 2: kTVX, kTFBorder, kITSROFBorder
  TH1 *hTFBorder     = h2EventSelection->ProjectionX("hTFBorder",     5, 5);  // kTFBorder
  TH1 *hITSROFBorder = h2EventSelection->ProjectionX("hITSROFBorder", 6, 6);  // kITSROFBorder

  if (!hINEL || !hNoRecoColl || !hTVX || !hTFBorder || !hITSROFBorder) {
    std::cout << "Error: Failed to project histograms from 2D histogram" << std::endl;
    file->Close();
    return;
  }

  // Rebin to common jet pT binning used elsewhere (ptbinGen, 26 bins)
  TH1 *hINELRebinned         = (TH1 *)hINEL->Rebin(26, "hINEL_rebinned", ptbinGen);
  TH1 *hNoRecoCollRebinned   = (TH1 *)hNoRecoColl->Rebin(26, "hNoRecoColl_rebinned", ptbinGen);
  TH1 *hTVXRebinned          = (TH1 *)hTVX->Rebin(26, "hTVX_rebinned", ptbinGen);
  TH1 *hTFBorderRebinned     = (TH1 *)hTFBorder->Rebin(26, "hTFBorder_rebinned", ptbinGen);
  TH1 *hITSROFBorderRebinned = (TH1 *)hITSROFBorder->Rebin(26, "hITSROFBorder_rebinned", ptbinGen);

  if (!hINELRebinned || !hNoRecoCollRebinned || !hTVXRebinned || !hTFBorderRebinned || !hITSROFBorderRebinned) {
    std::cout << "Error: Rebin failed" << std::endl;
    file->Close();
    return;
  }

  // Replace originals with rebinned and detach from file
  delete hINEL;
  delete hNoRecoColl;
  delete hTVX;
  delete hTFBorder;
  delete hITSROFBorder;

  hINEL         = hINELRebinned;
  hNoRecoColl   = hNoRecoCollRebinned;
  hTVX          = hTVXRebinned;
  hTFBorder     = hTFBorderRebinned;
  hITSROFBorder = hITSROFBorderRebinned;

  hINEL->SetDirectory(0);
  hNoRecoColl->SetDirectory(0);
  hTVX->SetDirectory(0);
  hTFBorder->SetDirectory(0);
  hITSROFBorder->SetDirectory(0);

  // Jet counts (integrated) in 0–50 GeV/c BEFORE any normalization
  int binLow = hTVX->GetXaxis()->FindBin(0.0);
  int binHigh = hTVX->GetXaxis()->FindBin(50.0);

  double nJetsINEL         = hINEL->Integral(binLow, binHigh);
  double nJetsNoRecoColl   = hNoRecoColl->Integral(binLow, binHigh);
  double nJetsTVX          = hTVX->Integral(binLow, binHigh);
  double nJetsTFBorder     = hTFBorder->Integral(binLow, binHigh);
  double nJetsITSROFBorder = hITSROFBorder->Integral(binLow, binHigh);

  // Jet counts for analysis range (10-140 GeV) for legend
  int binLowAnalysis = hTVX->GetXaxis()->FindBin(10.0);
  int binHighAnalysis = hTVX->GetXaxis()->FindBin(140.0);

  double nJetsTVXAnalysis      = hTVX->Integral(binLowAnalysis, binHighAnalysis);
  double nJetsTFBorderAnalysis = hTFBorder->Integral(binLowAnalysis, binHighAnalysis);
  double nJetsITSROFBorderAnalysis = hITSROFBorder->Integral(binLowAnalysis, binHighAnalysis);

  std::cout << "[CrossSectionEfficiencyINELCollNonSplitTVX] Jet counts (0–50 GeV/c):" << std::endl;
  std::cout << "  INEL         : " << nJetsINEL << std::endl;
  std::cout << "  noRecoColl   : " << nJetsNoRecoColl << std::endl;
  std::cout << "  kTVX         : " << nJetsTVX << std::endl;
  std::cout << "  kTFBorder    : " << nJetsTFBorder << std::endl;
  std::cout << "  kITSROFBorder : " << nJetsITSROFBorder << std::endl;

  // Clone for raw jet-count spectra (only scale by bin width when plotting)
  TH1 *hINELRaw         = (TH1 *)hINEL->Clone("hINELRaw");
  TH1 *hNoRecoCollRaw   = (TH1 *)hNoRecoColl->Clone("hNoRecoCollRaw");
  TH1 *hTVXRaw          = (TH1 *)hTVX->Clone("hTVXRaw");
  TH1 *hTFBorderRaw     = (TH1 *)hTFBorder->Clone("hTFBorderRaw");
  TH1 *hITSROFBorderRaw = (TH1 *)hITSROFBorder->Clone("hITSROFBorderRaw");

  // ===================== (1) Plot 1: INEL, noRecoColl, kTVX - ratios (selection / INEL) =====================
  {
    static int nn = 0;
    TCanvas *c1 = new TCanvas(Form("MCP_JetCounts_INEL_noRecoColl_TVX_%d", ++nn), "Jet Counts: INEL vs +Collisions vs +TVX", 800, 600);
    c1->cd();
    gStyle->SetOptStat(0);
    gPad->SetMargin(0.12, 0.02, 0.15, 0.06);  // left, right, bottom, top - small top for title, large bottom for x title

    // Ratios: selection / INEL
    TH1 *histsRaw[] = {hNoRecoCollRaw, hTVXRaw};
    Color_t colorsRaw[] = {kBlack, kBlue, kRed};  // Black for INEL (denominator)
    const char *labelsRaw[] = {"INEL", "+Collisions", "+TVX"};

    // Fit range: 10-140 GeV
    const double FitPtMin1 = 10.0;
    const double FitPtMax1 = 140.0;
    double fitValues[3] = {0, 0, 0};
    double fitErrors[3] = {0, 0, 0};
    TH1 *ratioHists[3] = {nullptr, nullptr, nullptr};  // Store ratio histograms for legend

    // First draw INEL/INEL = 1 (denominator reference)
    TH1 *ratioHistINEL = DrawRatioTH1(hINELRaw, hINELRaw);
    ratioHistINEL->SetDirectory(0);
    ratioHistINEL->SetTitle("");
    ratioHists[0] = ratioHistINEL;
    hset(*ratioHistINEL, JetPtGenTitleX, "selection / INEL", 1.2, 1.0,
         0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    ratioHistINEL->GetXaxis()->SetRangeUser(PtMinLocal, PtMaxLocal);
    ratioHistINEL->GetYaxis()->SetRangeUser(0.8, 1.2);
    ratioHistINEL->SetMarkerColor(colorsRaw[0]);
    ratioHistINEL->SetLineColor(colorsRaw[0]);
    ratioHistINEL->SetMarkerSize(0.8);
    ratioHistINEL->SetLineWidth(2);
    ratioHistINEL->SetFillColor(0);
    ratioHistINEL->Draw("pe");

    // Fit for INEL/INEL
    TF1 *fConstINEL = new TF1("fConst_INEL", "[0]", FitPtMin1, FitPtMax1);
    fConstINEL->SetParameter(0, 1.0);
    fConstINEL->SetLineColor(colorsRaw[0]);
    fConstINEL->SetLineStyle(2);
    fConstINEL->SetLineWidth(4);
    ratioHistINEL->Fit(fConstINEL, "Q0", "", FitPtMin1, FitPtMax1);
    fConstINEL->Draw("same");
    fitValues[0] = fConstINEL->GetParameter(0);
    fitErrors[0] = fConstINEL->GetParError(0);

    // Then draw other ratios
    for (int i = 0; i < 2; ++i) {
      TH1 *ratioHist = DrawRatioTH1(histsRaw[i], hINELRaw);
      ratioHist->SetDirectory(0);
      ratioHist->SetTitle("");  // Remove histogram title
      ratioHists[i+1] = ratioHist;  // Store for legend
      ratioHist->SetMarkerColor(colorsRaw[i+1]);
      ratioHist->SetLineColor(colorsRaw[i+1]);
      ratioHist->SetMarkerSize(0.8);
      ratioHist->SetLineWidth(2);
      ratioHist->SetFillColor(0);
      ratioHist->Draw("pesame");

      // Fit constant in 10-140 GeV range
      TF1 *fConst = new TF1(Form("fConst_%d", i), "[0]", FitPtMin1, FitPtMax1);
      fConst->SetParameter(0, 1.0);
      fConst->SetLineColor(colorsRaw[i+1]);
      fConst->SetLineStyle(2);
      fConst->SetLineWidth(4);
      ratioHist->Fit(fConst, "Q0", "", FitPtMin1, FitPtMax1);
      fConst->Draw("same");

      // Store fit results
      fitValues[i+1] = fConst->GetParameter(0);
      fitErrors[i+1] = fConst->GetParError(0);
    }

    // Title
    TLatex *title1 = new TLatex(0.5, 0.98, "Trigger efficiency for jets");
    title1->SetNDC();
    title1->SetTextAlign(22);  // Center align
    title1->SetTextSize(0.04);
    title1->SetTextFont(42);
    title1->Draw();

    // ALICE figure label (left top) - moved closer to corner
    ALICEfigureLegend("ALICE simulation WIP", 0.075188, 0.706957, 0.374687, 0.916522, 0.066416, 0.175652, 0.365915, 0.276522, 0.04, 0.4);

    // Create legend with fit results (right top) - use drawn histograms for LPE
    TLegend *legRaw = new TLegend(0.570175, 0.726957, 0.919799, 0.926957, NULL, "brNDC");
    legRaw->SetTextSize(0.036);
    legRaw->SetBorderSize(0);
    legRaw->SetHeader("MC particle-level (Fit [10-140 GeV])");
    legRaw->AddEntry(ratioHists[0], Form("%s (%.4f #pm %.4f)", labelsRaw[0], fitValues[0], fitErrors[0]), "lpe");
    legRaw->AddEntry(ratioHists[1], Form("%s (%.4f #pm %.4f)", labelsRaw[1], fitValues[1], fitErrors[1]), "lpe");
    legRaw->AddEntry(ratioHists[2], Form("%s (%.4f #pm %.4f)", labelsRaw[2], fitValues[2], fitErrors[2]), "lpe");
    legRaw->Draw();

    c1->Modified();
    c1->Update();
    c1->Print(Form("plots/MCP_jet_part_spectra_INEL_vs_noRecoColl_vs_kTVX_counts_%s.pdf",datasetName));
  }

  // ===================== (2) Plot 2: kTVX, kTFBorder, kITSROFBorder - ratios (selection / TVX) =====================
  {
  static int nn2 = 0;
    TCanvas *c2 = new TCanvas(Form("MCP_JetCounts_TVX_TFBorder_ITSROFBorder_%d", ++nn2), "Jet Counts: TVX vs +NoTFBorder vs +NoITSROFBorder", 800, 600);
    c2->cd();
    gStyle->SetOptStat(0);
    gPad->SetMargin(0.12, 0.02, 0.15, 0.06);  // left, right, bottom, top - small top for title, large bottom for x title

    // Ratios: selection / TVX
    // Plot range: 0-200 GeV
    const double PlotPtMax2 = 200.0;

    TH1 *histsRaw2[] = {hTFBorderRaw, hITSROFBorderRaw};
    Color_t colorsRaw2[] = {kBlack, kBlue, kRed};  // Black for TVX (denominator)
    const char *labelsRaw2[] = {"TVX", "+NoTFBorder", "+NoITSROFBorder"};
    TH1 *ratioHists2[3] = {nullptr, nullptr, nullptr};  // Store ratio histograms for legend

    // Fit range: 10-140 GeV
    const double FitPtMin2 = 10.0;
    const double FitPtMax2 = 140.0;
    double fitValues2[3] = {0, 0, 0};
    double fitErrors2[3] = {0, 0, 0};

    // First draw TVX/TVX = 1 (denominator reference)
    TH1 *ratioHistTVX = DrawRatioTH1(hTVXRaw, hTVXRaw);
    ratioHistTVX->SetDirectory(0);
    ratioHistTVX->SetTitle("");
    ratioHists2[0] = ratioHistTVX;
    hset(*ratioHistTVX, JetPtGenTitleX, "selection / TVX", 1.2, 1.0,
         0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    ratioHistTVX->GetXaxis()->SetRangeUser(PtMinLocal, PlotPtMax2);
    ratioHistTVX->GetYaxis()->SetRangeUser(0.5, 1.5);
    ratioHistTVX->SetMarkerColor(colorsRaw2[0]);
    ratioHistTVX->SetLineColor(colorsRaw2[0]);
    ratioHistTVX->SetMarkerSize(0.8);
    ratioHistTVX->SetLineWidth(2);
    ratioHistTVX->SetFillColor(0);
    ratioHistTVX->Draw("pe");

    // Fit for TVX/TVX
    TF1 *fConstTVX = new TF1("fConst_TVX_Plot2", "[0]", FitPtMin2, FitPtMax2);
    fConstTVX->SetParameter(0, 1.0);
    fConstTVX->SetLineColor(colorsRaw2[0]);
    fConstTVX->SetLineStyle(2);
    fConstTVX->SetLineWidth(4);
    ratioHistTVX->Fit(fConstTVX, "Q0", "", FitPtMin2, FitPtMax2);
    fConstTVX->Draw("same");
    fitValues2[0] = fConstTVX->GetParameter(0);
    fitErrors2[0] = fConstTVX->GetParError(0);

    // Then draw other ratios
    for (int i = 0; i < 2; ++i) {
      TH1 *ratioHist = DrawRatioTH1(histsRaw2[i], hTVXRaw);
      ratioHist->SetDirectory(0);
      ratioHist->SetTitle("");  // Remove histogram title
      ratioHists2[i+1] = ratioHist;  // Store for legend
      ratioHist->SetMarkerColor(colorsRaw2[i+1]);
      ratioHist->SetLineColor(colorsRaw2[i+1]);
      ratioHist->SetMarkerSize(0.8);
      ratioHist->SetLineWidth(2);
      ratioHist->SetFillColor(0);
      ratioHist->Draw("pesame");

      // Fit constant in 10-140 GeV range
      TF1 *fConst = new TF1(Form("fConst_Plot2_%d", i), "[0]", FitPtMin2, FitPtMax2);
      fConst->SetParameter(0, 1.0);
      fConst->SetLineColor(colorsRaw2[i+1]);
      fConst->SetLineStyle(2);
      fConst->SetLineWidth(4);
      ratioHist->Fit(fConst, "Q0", "", FitPtMin2, FitPtMax2);
      fConst->Draw("same");
      fitValues2[i+1] = fConst->GetParameter(0);
      fitErrors2[i+1] = fConst->GetParError(0);
    }

    // Title
    TLatex *title2 = new TLatex(0.5, 0.98, "Selection efficiency for jets");
    title2->SetNDC();
    title2->SetTextAlign(22);  // Center align
    title2->SetTextSize(0.04);
    title2->SetTextFont(42);
    title2->Draw();

    // ALICE figure label (left top) - moved closer to corner
    ALICEfigureLegend("ALICE simulation WIP", 0.075188, 0.706957, 0.374687, 0.916522, 0.066416, 0.175652, 0.365915, 0.276522, 0.04, 0.4);

    // Legend (right top) - Fit [10-140 GeV] format
    TLegend *legRaw2 = new TLegend(0.486216, 0.712174, 0.83584, 0.912174, NULL, "brNDC");
    legRaw2->SetTextSize(0.036);
    legRaw2->SetBorderSize(0);
    legRaw2->SetHeader("MC particle-level (Fit [10-140 GeV])");
    legRaw2->AddEntry(ratioHists2[0], Form("%s (%.4f #pm %.4f)", labelsRaw2[0], fitValues2[0], fitErrors2[0]), "lpe");
    legRaw2->AddEntry(ratioHists2[1], Form("%s (%.4f #pm %.4f)", labelsRaw2[1], fitValues2[1], fitErrors2[1]), "lpe");
    legRaw2->AddEntry(ratioHists2[2], Form("%s (%.4f #pm %.4f)", labelsRaw2[2], fitValues2[2], fitErrors2[2]), "lpe");
    legRaw2->Draw();

    c2->Modified();
    c2->Update();
    c2->Print(Form("plots/MCP_jet_part_spectra_kTVX_vs_kTFBorder_vs_kITSROFBorder_counts_%s.pdf",datasetName));
  }

  // ===================== (3) Invariant yields: kTVX, kTFBorder, kITSROFBorder - ratios (selection / TVX) =====================
  {
    // Clone spectra for yield (will be normalised by N_evt and width)
    TH1 *hTVXYield          = (TH1 *)hTVX->Clone("hTVXYield");
    TH1 *hTFBorderYield     = (TH1 *)hTFBorder->Clone("hTFBorderYield");
    TH1 *hITSROFBorderYield = (TH1 *)hITSROFBorder->Clone("hITSROFBorderYield");

    // Normalize by N_evt and width
    if (NmccollTVX > 0) {
      hTVXYield->Scale(1.0 / NmccollTVX, "width");
    }
    if (NmccollTFBorder > 0) {
      hTFBorderYield->Scale(1.0 / NmccollTFBorder, "width");
    }
    if (NmccollITSROFBorder > 0) {
      hITSROFBorderYield->Scale(1.0 / NmccollITSROFBorder, "width");
    }

    static int nn3 = 0;
    TCanvas *c3 = new TCanvas(Form("MCP_JetYield_TVX_TFBorder_ITSROFBorder_%d", ++nn3), "Invariant Yield: TVX vs +NoTFBorder vs +NoITSROFBorder", 800, 600);
    c3->cd();
    gStyle->SetOptStat(0);
    gPad->SetMargin(0.12, 0.02, 0.15, 0.06);  // left, right, bottom, top - small top for title, large bottom for x title

    // Ratios of invariant yield: selection / kTVX
    TH1 *histsYield[] = {hTFBorderYield, hITSROFBorderYield};
    Color_t colorsYield[] = {kBlack, kBlue, kRed};  // Black for TVX (denominator)
    const char *labelsYield[] = {"TVX", "+NoTFBorder", "+NoITSROFBorder"};
    TH1 *ratioHistsYield[3] = {nullptr, nullptr, nullptr};  // Store ratio histograms for legend

    // Fit range: 10-140 GeV
    const double FitPtMin3 = 10.0;
    const double FitPtMax3 = 140.0;
    double fitValues3[3] = {0, 0, 0};
    double fitErrors3[3] = {0, 0, 0};

    // First draw TVX/TVX = 1 (denominator reference)
    TH1 *ratioHistTVXYield = DrawRatioTH1(hTVXYield, hTVXYield);
    ratioHistTVXYield->SetDirectory(0);
    ratioHistTVXYield->SetTitle("");
    ratioHistsYield[0] = ratioHistTVXYield;
    hset(*ratioHistTVXYield, JetPtGenTitleX, "yield (sel.) / yield (TVX)", 1.2, 1.0,
         0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    ratioHistTVXYield->GetXaxis()->SetRangeUser(PtMinLocal, PtMaxLocal);
    ratioHistTVXYield->GetYaxis()->SetRangeUser(0.945, 1.055);
    ratioHistTVXYield->SetMarkerColor(colorsYield[0]);
    ratioHistTVXYield->SetLineColor(colorsYield[0]);
    ratioHistTVXYield->SetMarkerSize(0.8);
    ratioHistTVXYield->SetLineWidth(2);
    ratioHistTVXYield->SetFillColor(0);
    ratioHistTVXYield->Draw("pe");

    // Fit for TVX/TVX
    TF1 *fConstTVXYield = new TF1("fConst_TVXYield_Plot3", "[0]", FitPtMin3, FitPtMax3);
    fConstTVXYield->SetParameter(0, 1.0);
    fConstTVXYield->SetLineColor(colorsYield[0]);
    fConstTVXYield->SetLineStyle(2);
    fConstTVXYield->SetLineWidth(4);
    ratioHistTVXYield->Fit(fConstTVXYield, "Q0", "", FitPtMin3, FitPtMax3);
    fConstTVXYield->Draw("same");
    fitValues3[0] = fConstTVXYield->GetParameter(0);
    fitErrors3[0] = fConstTVXYield->GetParError(0);

    // Then draw other ratios
    for (int i = 0; i < 2; ++i) {
      TH1 *ratioHist = DrawRatioTH1(histsYield[i], hTVXYield);
      ratioHist->SetDirectory(0);
      ratioHist->SetTitle("");  // Remove histogram title
      ratioHistsYield[i+1] = ratioHist;  // Store for legend
      ratioHist->SetMarkerColor(colorsYield[i+1]);
      ratioHist->SetLineColor(colorsYield[i+1]);
      ratioHist->SetMarkerSize(0.8);
      ratioHist->SetLineWidth(2);
      ratioHist->SetFillColor(0);
      ratioHist->Draw("pesame");

      // Fit constant in 10-140 GeV range
      TF1 *fConst = new TF1(Form("fConst_Plot3_%d", i), "[0]", FitPtMin3, FitPtMax3);
      fConst->SetParameter(0, 1.0);
      fConst->SetLineColor(colorsYield[i+1]);
      fConst->SetLineStyle(2);
      fConst->SetLineWidth(4);
      ratioHist->Fit(fConst, "Q0", "", FitPtMin3, FitPtMax3);
      fConst->Draw("same");
      fitValues3[i+1] = fConst->GetParameter(0);
      fitErrors3[i+1] = fConst->GetParError(0);
    }

    // Title
    TLatex *title3 = new TLatex(0.5, 0.98, "Selection efficiency for jet yields");
    title3->SetNDC();
    title3->SetTextAlign(22);  // Center align
    title3->SetTextSize(0.04);
    title3->SetTextFont(42);
    title3->Draw();

    // ALICE figure label (left top) - moved closer to corner
    ALICEfigureLegend("ALICE simulation WIP", 0.075188, 0.706957, 0.374687, 0.916522, 0.066416, 0.175652, 0.365915, 0.276522, 0.04, 0.4);;

    // Legend (right top) - Fit [10-140 GeV] format
    TLegend *legYield = new TLegend(0.486216, 0.712174, 0.83584, 0.912174, NULL, "brNDC");
    legYield->SetTextSize(0.036);
  legYield->SetBorderSize(0);
    legYield->SetHeader("MC particle-level (Fit [10-140 GeV])");
    legYield->AddEntry(ratioHistsYield[0], Form("%s (%.4f #pm %.4f)", labelsYield[0], fitValues3[0], fitErrors3[0]), "lpe");
    legYield->AddEntry(ratioHistsYield[1], Form("%s (%.4f #pm %.4f)", labelsYield[1], fitValues3[1], fitErrors3[1]), "lpe");
    legYield->AddEntry(ratioHistsYield[2], Form("%s (%.4f #pm %.4f)", labelsYield[2], fitValues3[2], fitErrors3[2]), "lpe");
  legYield->Draw();
  
    c3->Modified();
    c3->Update();
    c3->Print(Form("plots/MCP_jet_part_spectra_kTVX_vs_kTFBorder_vs_kITSROFBorder_invariantYield_%s.pdf",datasetName));
  }

  // ===================== (4) Yield ratio plot: INEL, noRecoColl, kTVX - ratios (selection / INEL) =====================
  {
    // Clone spectra for yield (will be normalised by N_evt and width)
    TH1 *hINELYield         = (TH1 *)hINEL->Clone("hINELYield");
    TH1 *hNoRecoCollYield   = (TH1 *)hNoRecoColl->Clone("hNoRecoCollYield");
    TH1 *hTVXYield          = (TH1 *)hTVX->Clone("hTVXYield2");

    // Normalize by N_evt and width
    if (NmccollINEL > 0) {
      hINELYield->Scale(1.0 / NmccollINEL, "width");
    }
    if (NmccollNoRecoColl > 0) {
      hNoRecoCollYield->Scale(1.0 / NmccollNoRecoColl, "width");
    }
    if (NmccollTVX > 0) {
      hTVXYield->Scale(1.0 / NmccollTVX, "width");
    }

    static int nn4 = 0;
    TCanvas *c4 = new TCanvas(Form("MCP_JetYield_INEL_noRecoColl_TVX_%d", ++nn4), "Invariant Yield: INEL vs +Collisions vs +TVX", 800, 600);
    c4->cd();
    gStyle->SetOptStat(0);
    gPad->SetMargin(0.12, 0.02, 0.15, 0.06);  // left, right, bottom, top - small top for title, large bottom for x title

    // Ratios of invariant yield: selection / INEL
    TH1 *histsYield2[] = {hNoRecoCollYield, hTVXYield};
    Color_t colorsYield2[] = {kBlack, kBlue, kRed};  // Black for INEL (denominator)
    const char *labelsYield2[] = {"INEL", "+Collisions", "+TVX"};
    TH1 *ratioHistsYield2[3] = {nullptr, nullptr, nullptr};  // Store ratio histograms for legend

    // Fit range: 10-140 GeV
    const double FitPtMin4 = 10.0;
    const double FitPtMax4 = 140.0;
    double fitValues4[3] = {0, 0, 0};
    double fitErrors4[3] = {0, 0, 0};

    // First draw INEL/INEL = 1 (denominator reference)
    TH1 *ratioHistINELYield = DrawRatioTH1(hINELYield, hINELYield);
    ratioHistINELYield->SetDirectory(0);
    ratioHistINELYield->SetTitle("");
    ratioHistsYield2[0] = ratioHistINELYield;
    hset(*ratioHistINELYield, JetPtGenTitleX, "yield (sel.) / yield (INEL)", 1.2, 1.0,
         0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    ratioHistINELYield->GetXaxis()->SetRangeUser(PtMinLocal, PtMaxLocal);
    ratioHistINELYield->GetYaxis()->SetRangeUser(0.5, 2.0);
    ratioHistINELYield->SetMarkerColor(colorsYield2[0]);
    ratioHistINELYield->SetLineColor(colorsYield2[0]);
    ratioHistINELYield->SetMarkerSize(0.8);
    ratioHistINELYield->SetLineWidth(2);
    ratioHistINELYield->SetFillColor(0);
    ratioHistINELYield->Draw("pe");

    // Fit for INEL/INEL
    TF1 *fConstINELYield = new TF1("fConst_INELYield_Plot4", "[0]", FitPtMin4, FitPtMax4);
    fConstINELYield->SetParameter(0, 1.0);
    fConstINELYield->SetLineColor(colorsYield2[0]);
    fConstINELYield->SetLineStyle(2);
    fConstINELYield->SetLineWidth(4);
    ratioHistINELYield->Fit(fConstINELYield, "Q0", "", FitPtMin4, FitPtMax4);
    fConstINELYield->Draw("same");
    fitValues4[0] = fConstINELYield->GetParameter(0);
    fitErrors4[0] = fConstINELYield->GetParError(0);

    // Then draw other ratios
    for (int i = 0; i < 2; ++i) {
      TH1 *ratioHist = DrawRatioTH1(histsYield2[i], hINELYield);
      ratioHist->SetDirectory(0);
      ratioHist->SetTitle("");  // Remove histogram title
      ratioHistsYield2[i+1] = ratioHist;  // Store for legend
      ratioHist->SetMarkerColor(colorsYield2[i+1]);
      ratioHist->SetLineColor(colorsYield2[i+1]);
      ratioHist->SetMarkerSize(0.8);
      ratioHist->SetLineWidth(2);
      ratioHist->SetFillColor(0);
      ratioHist->Draw("pesame");

      // Fit constant in 10-140 GeV range
      TF1 *fConst = new TF1(Form("fConst_Plot4_%d", i), "[0]", FitPtMin4, FitPtMax4);
      fConst->SetParameter(0, 1.0);
      fConst->SetLineColor(colorsYield2[i+1]);
      fConst->SetLineStyle(2);
      fConst->SetLineWidth(4);
      ratioHist->Fit(fConst, "Q0", "", FitPtMin4, FitPtMax4);
      fConst->Draw("same");
      fitValues4[i+1] = fConst->GetParameter(0);
      fitErrors4[i+1] = fConst->GetParError(0);
    }

    // Title
    TLatex *title4 = new TLatex(0.5, 0.98, "Trigger efficiency for jet yields");
    title4->SetNDC();
    title4->SetTextAlign(22);  // Center align
    title4->SetTextSize(0.04);
    title4->SetTextFont(42);
    title4->Draw();

    // ALICE figure label (left bottom) - special for this plot
    ALICEfigureLegend("ALICE simulation WIP", 0.075188, 0.706957, 0.374687, 0.916522, 0.066416, 0.175652, 0.365915, 0.276522, 0.04, 0.4);

    // Legend (left bottom) - special for this plot - Fit [10-140 GeV] format
    TLegend *legYield2 = new TLegend(0.570175, 0.726957, 0.919799, 0.926957, NULL, "brNDC");
    legYield2->SetTextSize(0.036);
    legYield2->SetBorderSize(0);
    legYield2->SetHeader("MC particle-level (Fit [10-140 GeV])");
    legYield2->AddEntry(ratioHistsYield2[0], Form("%s (%.4f #pm %.4f)", labelsYield2[0], fitValues4[0], fitErrors4[0]), "lpe");
    legYield2->AddEntry(ratioHistsYield2[1], Form("%s (%.4f #pm %.4f)", labelsYield2[1], fitValues4[1], fitErrors4[1]), "lpe");
    legYield2->AddEntry(ratioHistsYield2[2], Form("%s (%.4f #pm %.4f)", labelsYield2[2], fitValues4[2], fitErrors4[2]), "lpe");
    legYield2->Draw();

    c4->Modified();
    c4->Update();
    c4->Print(Form("plots/MCP_jet_part_spectra_INEL_vs_noRecoColl_vs_kTVX_invariantYield_%s.pdf",datasetName));
  }

  file->Close();
}

void CrossSectionEfficiency() {
  // New updated selection breakdown using jetSpectraCharged.cxx
  CrossSectionEfficiencyINELCollNonSplitTVX();
  // Legacy plots (can be kept for comparison if needed)
  // CrossSectionEfficiencyINELsel8jetSpectraCharged();
  // CrossSectionEfficiencyINELsel8jetFinderQA();
  // CrossSectionEfficiencySelMCSel8();
}

#else
// This file is intended to be executed as a ROOT macro (cling). To run:
// root -l 'CrossSectionEfficiency.C()'
#endif