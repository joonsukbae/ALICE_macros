// Compare MCP histograms from two AnalysisResults files
// Compares: 496218 vs 498287
// Histogram: jet-finder-charged-qa/h_jet_pt_part

#if defined(__CLING__) || defined(__CINT__) || defined(__ROOTCLING__)

#include "Filipad2.h"
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

// Helper functions (simplified versions)
void optFili(TPad &pid, Int_t gridx, Int_t gridy, Int_t logx, Int_t logy) {
  pid.SetGridy(gridx);
  pid.SetGridx(gridy);
  pid.SetLogx(logx);
  pid.SetLogy(logy);
}

template <typename T>
void hset(T &hid, TString xtit = "", TString ytit = "", double titoffx = 0.9,
          double titoffy = 1.2, double titsizex = 0.06, double titsizey = 0.06,
          double labeloffx = 0.01, double labeloffy = 0.001,
          double labelsizex = 0.05, double labelsizey = 0.05, int divx = 510,
          int divy = 510) {
  hid.GetXaxis()->CenterTitle(1);
  hid.GetYaxis()->CenterTitle(1);
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

template <typename T>
void hoptset(T &hid, Double_t N = 1, Color_t color = kBlack, Double_t minX = 0,
             Double_t maxX = 100, Double_t minY = 0, Double_t maxY = 1,
             Double_t MarkerSize = .75, Int_t logy = 1, Int_t LineStyle = 1, 
             Int_t LineWidth = 1, Int_t MarkerStyle = 20) {
  if (N != 0) {
    double integral = hid.Integral();
    if (integral > 0) {
      if (N == 1) {
        hid.Scale(1. / integral, "width");
      } else if (N == 2) {
        hid.Scale(1. / integral, "");
      } else {
        hid.Scale(1. / N, "width");
      }
    }
  }
  hid.SetMarkerColor(color);
  hid.SetLineColor(color);
  hid.SetMarkerSize(MarkerSize);
  hid.SetMarkerStyle(MarkerStyle); 
  hid.SetLineStyle(LineStyle); 
  hid.SetLineWidth(LineWidth);
  hid.GetXaxis()->SetRangeUser(minX, maxX);
  if (minY > 0 && maxY > minY) {
    hid.GetYaxis()->SetRangeUser(minY, maxY);
  }
  hid.SetFillColor(0);  // No fill for point plots
}

TH1 *DrawRatioTH1(TH1 *hNum, TH1 *hDenom) {
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
            if (TMath::Abs(loA - loB) > 1e-9 || TMath::Abs(hiA - hiB) > 1e-9) return false;
        }
        return true;
    };

    if (edgesEqual(nN, axN, nD, axD)) {
        TH1 *hRatio = (TH1 *)hNum->Clone("hRatio");
        hRatio->Reset();
        for (int i = 1; i <= nN; i++) {
            double yNum = hNum->GetBinContent(i);
            double eNum = hNum->GetBinError(i);
            double yDen = hDenom->GetBinContent(i);
            double eDen = hDenom->GetBinError(i);
            if (yDen != 0) {
                double ratio = yNum / yDen;
                // Proper error propagation: sqrt((eNum/yDen)^2 + (yNum*eDen/yDen^2)^2)
                double ratioErr = 0.0;
                if (yNum > 0 && yDen > 0) {
                    double term1 = eNum / yDen;
                    double term2 = (yNum * eDen) / (yDen * yDen);
                    ratioErr = TMath::Sqrt(term1*term1 + term2*term2);
                }
                hRatio->SetBinContent(i, ratio);
                hRatio->SetBinError(i, ratioErr);
            } else {
                hRatio->SetBinContent(i, 0);
                hRatio->SetBinError(i, 0);
            }
        }
        return hRatio;
    }

    // Build common bin edges
    std::vector<double> eN(nN + 1), eD(nD + 1);
    for (int i = 0; i <= nN; ++i) eN[i] = (i < nN ? axN->GetBinLowEdge(i + 1) : axN->GetXmax());
    for (int i = 0; i <= nD; ++i) eD[i] = (i < nD ? axD->GetBinLowEdge(i + 1) : axD->GetXmax());
    eN[0] = axN->GetXmin(); eD[0] = axD->GetXmin();

    auto almostEqual = [](double a, double b) { return TMath::Abs(a - b) < 1e-9; };
    std::vector<double> common;
    for (double en : eN) {
        for (double ed : eD) {
            if (almostEqual(en, ed)) { common.push_back(en); break; }
        }
    }
    std::sort(common.begin(), common.end());
    common.erase(std::unique(common.begin(), common.end(), almostEqual), common.end());

    if (common.size() < 2) {
        TH1 *hRatio = (TH1 *)hNum->Clone("hRatio");
        hRatio->Reset();
        for (int i = 1; i <= nN; i++) {
            double xCenter = axN->GetBinCenter(i);
            double yNum = hNum->GetBinContent(i);
            double eNum = hNum->GetBinError(i);
            int ib = axD->FindBin(xCenter);
            double yDen = hDenom->GetBinContent(ib);
            double eDen = hDenom->GetBinError(ib);
            if (yDen != 0) {
                double ratio = yNum / yDen;
                // Proper error propagation
                double ratioErr = 0.0;
                if (yNum > 0 && yDen > 0) {
                    double term1 = eNum / yDen;
                    double term2 = (yNum * eDen) / (yDen * yDen);
                    ratioErr = TMath::Sqrt(term1*term1 + term2*term2);
                }
                hRatio->SetBinContent(i, ratio);
                hRatio->SetBinError(i, ratioErr);
            }
        }
        return hRatio;
    }

    int nC = static_cast<int>(common.size()) - 1;
    std::vector<double> edges = common;
    TH1 *numC = (TH1 *)hNum->Rebin(nC, Form("%s_rebinned_num", hNum->GetName()), edges.data());
    TH1 *denC = (TH1 *)hDenom->Rebin(nC, Form("%s_rebinned_den", hDenom->GetName()), edges.data());

    TH1 *hRatio = (TH1 *)numC->Clone("hRatio");
    hRatio->Reset();
    for (int i = 1; i <= nC; ++i) {
        double yNum = numC->GetBinContent(i);
        double eNum = numC->GetBinError(i);
        double yDen = denC->GetBinContent(i);
        double eDen = denC->GetBinError(i);
        if (yDen != 0) {
            double ratio = yNum / yDen;
            // Proper error propagation
            double ratioErr = 0.0;
            if (yNum > 0 && yDen > 0) {
                double term1 = eNum / yDen;
                double term2 = (yNum * eDen) / (yDen * yDen);
                ratioErr = TMath::Sqrt(term1*term1 + term2*term2);
            }
            hRatio->SetBinContent(i, ratio);
            hRatio->SetBinError(i, ratioErr);
        } else {
            hRatio->SetBinContent(i, 0);
            hRatio->SetBinError(i, 0);
        }
    }
    return hRatio;
}

TH1 *DrawRatio(const char *ratioName, TH1 *refHist, TH1 *testHist,
               TString AxisTitleX, TString AxisTitleY, Color_t colorID,
               Double_t Ymin, Double_t Ymax, Double_t MarkerSize = 1) {
  TH1 *ratioHist = DrawRatioTH1(testHist, refHist);
  ratioHist->SetName(ratioName);
  ratioHist->SetDirectory(0);
  hset(*ratioHist, AxisTitleX, AxisTitleY, 1.2, 1.0, 0.07, 0.07, 0.01, 0.01,
       0.07, 0.07, 510, 505);
  ratioHist->GetXaxis()->SetRangeUser(PlotPtMin, PlotPtMax);
  ratioHist->GetYaxis()->SetRangeUser(Ymin, Ymax);
  ratioHist->SetMarkerColor(colorID);
  ratioHist->SetLineColor(colorID);
  ratioHist->SetMarkerSize(MarkerSize);
  ratioHist->SetLineStyle(1);  // Solid line
  ratioHist->SetLineWidth(1);
  ratioHist->SetFillColor(0);
  // Draw with "e" first (not "same") to create the axis
  ratioHist->Draw("e");
  return ratioHist;
}

void MCPTestJetSpectraCharged() {
  // File path
  const char* fileName = "~/cernbox/workspace/O2Physics/jets/AnalysisResults/560874_AnalysisResults.root";
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
  TH1 *hInelRaw = (TH1*)hInel->Clone("hInelRaw");
  TH1 *hZvtxRaw = (TH1*)hZvtx->Clone("hZvtxRaw");
  TH1 *hNoRecoCollRaw = (TH1*)hNoRecoColl->Clone("hNoRecoCollRaw");
  TH1 *hSplitCollRaw = (TH1*)hSplitColl->Clone("hSplitCollRaw");
  TH1 *hRecoEvtSelRaw = (TH1*)hRecoEvtSel->Clone("hRecoEvtSelRaw");
  
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
      infoLatex->SetTextSize(0.07);
      infoLatex->SetTextColor(kBlack);
      infoLatex->SetTextAlign(12);  // Left and bottom aligned
      infoLatex->DrawLatex(0.25, 0.95, "10-200 GeV fit:");
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
  MCPPad->C->Print("plots/MCP_jet_part_spectra_steps_vs_SelMC.pdf");
  
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
  mainPadRaw->cd();
  
  // Setup first histogram (reference) - no normalization, just scale by width for display
  hset(*hInelRaw, JetPtGenTitleX, "dN_{ch-jet}^{true}/dp_{T}", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
  hoptset(*hInelRaw, 0, kRed, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 22);
  hInelRaw->Scale(1.0, "width");  // Only scale by width, not by events
  hInelRaw->SetLineStyle(1);
  hInelRaw->SetLineWidth(2);
  hInelRaw->Draw("pe");

  hoptset(*hZvtxRaw, 0, kBlue, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 21);
  hZvtxRaw->Scale(1.0, "width");
  hZvtxRaw->SetLineStyle(1);
  hZvtxRaw->SetLineWidth(2);
  hZvtxRaw->Draw("pesame");

  hoptset(*hNoRecoCollRaw, 0, kGreen+2, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 23);
  hNoRecoCollRaw->Scale(1.0, "width");
  hNoRecoCollRaw->SetLineStyle(1);
  hNoRecoCollRaw->SetLineWidth(2);
  hNoRecoCollRaw->Draw("pesame");

  hoptset(*hSplitCollRaw, 0, kMagenta+2, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 24);
  hSplitCollRaw->Scale(1.0, "width");
  hSplitCollRaw->SetLineStyle(1);
  hSplitCollRaw->SetLineWidth(2);
  hSplitCollRaw->Draw("pesame");

  hoptset(*hRecoEvtSelRaw, 0, kBlack, PlotPtMin, PlotPtMax, 1e-12, 1e+1, 0.75, 1, 1, 20);
  hRecoEvtSelRaw->Scale(1.0, "width");
  hRecoEvtSelRaw->SetLineStyle(1);
  hRecoEvtSelRaw->SetLineWidth(3);
  hRecoEvtSelRaw->Draw("pesame");
  
  // Create legend with jet counts only (no event counts)
  TLegend *legRaw = new TLegend(0.70, 0.40, 0.65, 0.95, NULL, "brNDC");
  legRaw->SetTextSize(0.045);
  legRaw->SetBorderSize(0);
  legRaw->AddEntry("", "MC particle-level", "");
  legRaw->AddEntry("", "N_{jets} [10-200 GeV]", "");
  legRaw->AddEntry(hInelRaw, Form("INEL (%.2e)", nJetsInel), "pe");
  legRaw->AddEntry(hZvtxRaw, Form("+Zvtx10 (%.2e)", nJetsZvtx), "pe");
  legRaw->AddEntry(hNoRecoCollRaw, Form("+hasColl (%.2e)", nJetsNoRecoColl), "pe");
  legRaw->AddEntry(hSplitCollRaw, Form("+noSplit (%.2e)", nJetsSplitColl), "pe");
  legRaw->AddEntry(hRecoEvtSelRaw, Form("+selMC (%.2e)", nJetsRecoEvtSel), "pe");
  legRaw->Draw();
  
  // Update main pad to ensure it's visible in ROOT window
  mainPadRaw->Update();
  
  // Draw ratio plot for raw jet counts
  ratioPadRaw->cd();
  
  // Arrays for looping - all steps vs selMC (recoEvtSel)
  TH1* histsRaw[] = {hInelRaw, hZvtxRaw, hNoRecoCollRaw, hSplitCollRaw};
  Color_t colorsRaw[] = {kRed, kBlue, kGreen+2, kMagenta+2};
  const char* labelsRaw[] = {"INEL", "Zvtx", "hasColl", "noSplit"};
  int nHistsRaw = 4;

  TString ratioNameRaw = "MCP_Ratio_Steps_vs_selMC_Raw";
  for (int i = 0; i < nHistsRaw; i++) {
    TString ratioName_i = Form("%s_%s", ratioNameRaw.Data(), labelsRaw[i]);
    TH1 *ratioHist = DrawRatioTH1(histsRaw[i], hRecoEvtSelRaw);
    ratioHist->SetName(ratioName_i.Data());
    ratioHist->SetDirectory(0);
    
    // Setup histogram properties
    if (i == 0) {
      // First histogram: setup axis and draw
      hset(*ratioHist, JetPtGenTitleX, "Ref. / (selMC)", 1.2, 1.0, 0.07, 0.06, 0.01, 0.01,
           0.07, 0.07, 510, 505);
      ratioHist->GetXaxis()->SetRangeUser(PlotPtMin, PlotPtMax);
      ratioHist->GetYaxis()->SetRangeUser(0.97, 1.35);
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

void MCPTestJetFinderChargedQA() {
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
  
  hset(*hSelMCClone, JetPtGenTitleX, "dN_{jet}/dp_{T}", 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
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
  leg->AddEntry(hSelMCClone, Form("selMC (N_{jets}=%.2e)", nJetsSelMC), "pe");
  leg->AddEntry(hSel8Clone, Form("sel8 (N_{jets}=%.2e)", nJetsSel8), "pe");
  leg->Draw();
  
  // Update main pad
  mainPad->Update();
  
  // Draw ratio plot (selMC / sel8)
  ratioPad->cd();
  
  TH1 *ratioHist = DrawRatioTH1(hSelMCClone, hSel8Clone);
  ratioHist->SetName("CrossSection_Ratio_SelMC_Sel8");
  ratioHist->SetDirectory(0);
  
  hset(*ratioHist, JetPtGenTitleX, "selMC / sel8", 1.2, 1.0, 0.07, 0.06, 0.01, 0.01, 0.07, 0.07, 510, 505);
  ratioHist->GetXaxis()->SetRangeUser(PlotPtMin, PlotPtMax);
  ratioHist->GetYaxis()->SetRangeUser(0.7, 1.3);
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

void CrossSectionEfficiency() {
  MCPTestJetSpectraCharged();
  MCPTestJetFinderChargedQA();
  CrossSectionEfficiencySelMCSel8();
}

#else
// This file is intended to be executed as a ROOT macro (cling). To run:
// root -l 'CrossSectionEfficiency.C()'
#endif