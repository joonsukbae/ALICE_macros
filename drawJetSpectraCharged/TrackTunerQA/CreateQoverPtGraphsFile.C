// Create Q/pT Correction Graphs File
// Purpose: Create ROOT file with sigmaVsPtMc and sigmaVsPtData graphs
//          in the same format as reference file for TrackTuner
// Author: Auto-generated
// Date: 2026-01-29

#include "../BSHelper.cxx"
#include "../Filipad2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TList.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstring>

// Histogram names (from DrawTrackResolutionsNew.C)
// IMPORTANT:
// - TrackTuner applies smearing on q/pT (a.k.a. 1/pT). The CCDB object that
//   TrackTuner consumes (graphs named sigmaVsPtMc / sigmaVsPtData) must contain
//   sigma(1/pT) values.
// - For QA / Run-1/2 comparisons it is often more intuitive to plot
//   sigma(pT)/pT. This macro therefore ALSO extracts sigma(pT)/pT and:
//     * stores it in the same ROOT file as extra graphs, and
//     * uses it for the QA PDF.
const char *DataTrackSigma1OverPtObj = "h2_track_pt_track_sigma1overpt";
const char *DataTrackSigma1OverPtHighObj = "h2_track_pt_high_track_sigma1overpt";

const char *DataTrackSigmaPtOverPtObj = "h2_track_pt_track_sigmapt";
const char *DataTrackSigmaPtOverPtHighObj = "h2_track_pt_high_track_sigmapt";

// Custom pT binning: use native TH2 bin width where statistics allow.
// TrackTuner.h uses TGraphErrors::Eval() (linear interpolation) so
// non-uniform spacing is perfectly fine.
// Native TH2 bins: low pT = 0.1 GeV (0-10), high pT = 1 GeV (10-100).
// NOTE: 10.0 GeV must remain an edge (boundary between low/high TH2).
static const std::vector<Double_t> kPtBinEdges = {
    // 0-10 GeV: native 0.1 GeV bins
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
    2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
    3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0,
    4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0,
    5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0,
    6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0,
    7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0,
    8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0,
    9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0,
    // 10-20 GeV: native 1 GeV bins
    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
    // 20-30 GeV: 2 GeV bins
    22.0, 24.0, 26.0, 28.0, 30.0,
    // 30-50 GeV: 5 GeV bins
    35.0, 40.0, 45.0, 50.0,
    // 50-100 GeV: 10 GeV bins
    60.0, 70.0, 80.0, 100.0
};

// Helper: project custom pT bins from a single TH2, appending to graph.
// Only processes custom bins that overlap with this TH2's x-range.
static void AddPointsFromTH2WithBinning(TH2* h2, TGraphErrors* graph,
                                         int& nPoints,
                                         const std::vector<Double_t>& binEdges)
{
  const Double_t h2Xmin = h2->GetXaxis()->GetXmin();
  const Double_t h2Xmax = h2->GetXaxis()->GetXmax();

  for (size_t iBin = 0; iBin + 1 < binEdges.size(); iBin++) {
    Double_t binLow  = binEdges[iBin];
    Double_t binHigh = binEdges[iBin + 1];

    // Skip custom bins outside this TH2's range
    if (binHigh <= h2Xmin || binLow >= h2Xmax) continue;

    // Clamp to TH2 range
    Double_t effLow  = std::max(binLow,  h2Xmin);
    Double_t effHigh = std::min(binHigh, h2Xmax);

    // Map to TH2 bin indices (small offset avoids landing on the edge)
    Int_t firstBin = h2->GetXaxis()->FindBin(effLow  + 1e-6);
    Int_t lastBin  = h2->GetXaxis()->FindBin(effHigh - 1e-6);
    if (firstBin > lastBin) continue;

    TH1* hProj = h2->ProjectionY(Form("proj_%d_%zu", nPoints, iBin),
                                  firstBin, lastBin, "e");
    if (hProj->GetEntries() < 10) { delete hProj; continue; }

    Double_t ptCenter  = (effLow + effHigh) / 2.0;
    Double_t meanSigma = hProj->GetMean();
    Double_t errSigma  = hProj->GetMeanError();

    if (meanSigma > 0 && ptCenter > 0) {
      graph->SetPoint(nPoints, ptCenter, meanSigma);
      graph->SetPointError(nPoints, (effHigh - effLow) / 2.0, errSigma);
      nPoints++;
    }
    delete hProj;
  }
}

static TGraphErrors* ExtractSigmaVsPtFromCovMatGeneric(const char* fileName,
                                                       const char* directory,
                                                       const char* histLowName,
                                                       const char* histHighName,
                                                       const char* graphName)
{
  TFile* file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open file " << fileName << std::endl;
    return nullptr;
  }

  TDirectory* dir = (TDirectory*)file->Get(directory);
  if (!dir) {
    std::cerr << "Error: Cannot find directory " << directory << std::endl;
    file->Close();
    return nullptr;
  }

  // Low pT
  TH2* h2Low = (TH2*)dir->Get(histLowName);
  if (!h2Low) {
    std::cerr << "Error: Cannot find " << histLowName << std::endl;
    file->Close();
    return nullptr;
  }

  // High pT
  TH2* h2High = (TH2*)dir->Get(histHighName);
  if (!h2High) {
    std::cerr << "Warning: Cannot find " << histHighName << std::endl;
  }

  TGraphErrors* graph = new TGraphErrors();
  graph->SetName(graphName);
  int nPoints = 0;

  // Process both TH2s with the same custom binning.
  // Each custom bin is handled by whichever TH2 covers it.
  AddPointsFromTH2WithBinning(h2Low, graph, nPoints, kPtBinEdges);
  if (h2High) {
    AddPointsFromTH2WithBinning(h2High, graph, nPoints, kPtBinEdges);
  }

  std::cout << "Extracted " << nPoints << " points from " << fileName
            << " (custom binning, " << kPtBinEdges.size() - 1 << " bins defined)"
            << std::endl;

  file->Close();
  return graph;
}

// Extract sigma(1/pT) vs pT
TGraphErrors* ExtractSigma1OverPtVsPtFromCovMat(const char* fileName, const char* directory)
{
  return ExtractSigmaVsPtFromCovMatGeneric(fileName,
                                           directory,
                                           DataTrackSigma1OverPtObj,
                                           DataTrackSigma1OverPtHighObj,
                                           "sigma1overptVsPt");
}

// Extract sigma(pT)/pT vs pT
TGraphErrors* ExtractSigmaPtOverPtVsPtFromCovMat(const char* fileName, const char* directory)
{
  return ExtractSigmaVsPtFromCovMatGeneric(fileName,
                                           directory,
                                           DataTrackSigmaPtOverPtObj,
                                           DataTrackSigmaPtOverPtHighObj,
                                           "sigmaptoverptVsPt");
}

// MC uses the same histogram names (covariance matrix)
TGraphErrors* ExtractSigma1OverPtVsPtFromCovMatMC(const char* fileName, const char* directory)
{
  return ExtractSigma1OverPtVsPtFromCovMat(fileName, directory);
}

TGraphErrors* ExtractSigmaPtOverPtVsPtFromCovMatMC(const char* fileName, const char* directory)
{
  return ExtractSigmaPtOverPtVsPtFromCovMat(fileName, directory);
}

static void PrintDataOverMcRatioSummary(const TGraphErrors* gData, const TGraphErrors* gMC, const char* label)
{
  if (!gData || !gMC) {
    return;
  }
  if (gData->GetN() < 2 || gMC->GetN() < 2) {
    return;
  }

  Double_t x0d = 0, y0d = 0, x1d = 0, y1d = 0;
  Double_t x0m = 0, y0m = 0, x1m = 0, y1m = 0;
  gData->GetPoint(0, x0d, y0d);
  gData->GetPoint(gData->GetN() - 1, x1d, y1d);
  gMC->GetPoint(0, x0m, y0m);
  gMC->GetPoint(gMC->GetN() - 1, x1m, y1m);
  const double xmin = std::max((double)x0d, (double)x0m);
  const double xmax = std::min((double)x1d, (double)x1m);

  const double testPts[] = {0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 80.0, 99.5};
  std::cout << "\n=== Data/MC ratio summary: " << label << " ===" << std::endl;
  std::cout << "  overlap pT range: [" << xmin << ", " << xmax << "]" << std::endl;

  for (double pt : testPts) {
    if (pt < xmin || pt > xmax) {
      continue;
    }
    const double yData = gData->Eval(pt);
    const double yMc = gMC->Eval(pt);
    if (!(yData > 0) || !(yMc > 0)) {
      continue;
    }
    std::cout << Form("  pT=%6.2f : data=%.6g  mc=%.6g  data/mc=%.4f", pt, yData, yMc, yData / yMc) << std::endl;
  }

  double rMin = 1e9;
  double rMax = 0.0;
  int nUsed = 0;
  for (int i = 0; i < gData->GetN(); ++i) {
    double x = 0.0, y = 0.0;
    gData->GetPoint(i, x, y);
    if (x < xmin || x > xmax) {
      continue;
    }
    const double yMc = gMC->Eval(x);
    if (!(y > 0) || !(yMc > 0)) {
      continue;
    }
    const double r = y / yMc;
    rMin = std::min(rMin, r);
    rMax = std::max(rMax, r);
    ++nUsed;
  }
  if (nUsed > 0) {
    std::cout << "  min/max data/mc over Data points (Eval MC): " << rMin << " / " << rMax << " (n=" << nUsed << ")" << std::endl;
  }
}

// Function to create and save plot (same format as AnalyzeReferenceQoverPtFile.C)
// Note: the ratio panel is always Data/MC, independent of y-units.
void CreateQoverPtPlot(TGraphErrors* gData,
                       TGraphErrors* gMC,
                       const char* outputFileName,
                       const char* yAxisTitle,
                       const char* tag = "")
{
  // Calculate ratio using common bin edges (same logic as AnalyzeReferenceQoverPtFile.C)
  std::vector<Double_t> edgesData, edgesMC;
  
  // Extract bin edges from Data graph
  for (int i = 0; i < gData->GetN(); i++) {
    Double_t x, y;
    gData->GetPoint(i, x, y);
    Double_t errX = gData->GetErrorX(i);
    edgesData.push_back(x - errX);
    edgesData.push_back(x + errX);
  }
  
  // Extract bin edges from MC graph
  for (int i = 0; i < gMC->GetN(); i++) {
    Double_t x, y;
    gMC->GetPoint(i, x, y);
    Double_t errX = gMC->GetErrorX(i);
    edgesMC.push_back(x - errX);
    edgesMC.push_back(x + errX);
  }
  
  // Sort and remove duplicates
  std::sort(edgesData.begin(), edgesData.end());
  edgesData.erase(std::unique(edgesData.begin(), edgesData.end()), edgesData.end());
  
  std::sort(edgesMC.begin(), edgesMC.end());
  edgesMC.erase(std::unique(edgesMC.begin(), edgesMC.end()), edgesMC.end());
  
  // Find common bin edges
  auto almostEqual = [](Double_t a, Double_t b) { return TMath::Abs(a - b) < 1e-6; };
  std::vector<Double_t> commonEdges;
  
  for (Double_t eData : edgesData) {
    for (Double_t eMC : edgesMC) {
      if (almostEqual(eData, eMC)) {
        commonEdges.push_back(eData);
        break;
      }
    }
  }
  
  // Also add edges that are close enough
  for (Double_t eData : edgesData) {
    bool found = false;
    for (Double_t eCommon : commonEdges) {
      if (almostEqual(eData, eCommon)) {
        found = true;
        break;
      }
    }
    if (!found) {
      for (Double_t eMC : edgesMC) {
        Double_t avgEdge = (eData + eMC) / 2.0;
        Double_t tolerance = TMath::Max(TMath::Abs(eData - eMC) * 0.5, 0.01);
        if (TMath::Abs(eData - eMC) < tolerance) {
          commonEdges.push_back(avgEdge);
          found = true;
          break;
        }
      }
    }
  }
  
  std::sort(commonEdges.begin(), commonEdges.end());
  commonEdges.erase(std::unique(commonEdges.begin(), commonEdges.end(), almostEqual), commonEdges.end());
  
  if (commonEdges.size() < 2) {
    commonEdges.clear();
    commonEdges.insert(commonEdges.end(), edgesData.begin(), edgesData.end());
    commonEdges.insert(commonEdges.end(), edgesMC.begin(), edgesMC.end());
    std::sort(commonEdges.begin(), commonEdges.end());
    commonEdges.erase(std::unique(commonEdges.begin(), commonEdges.end(), almostEqual), commonEdges.end());
  }
  
  // Convert graphs to histograms for rebinning
  int nCommonBins = commonEdges.size() - 1;
  if (nCommonBins < 1) {
    std::cerr << "Error: Not enough common bin edges!" << std::endl;
    return;
  }
  
  TH1D* hData = new TH1D("hData_rebin", "Data", nCommonBins, commonEdges.data());
  TH1D* hMC = new TH1D("hMC_rebin", "MC", nCommonBins, commonEdges.data());
  hData->SetDirectory(0);
  hMC->SetDirectory(0);
  
  // Fill histograms from graphs
  for (int i = 0; i < gData->GetN(); i++) {
    Double_t x, y;
    gData->GetPoint(i, x, y);
    Double_t errY = gData->GetErrorY(i);
    int bin = hData->FindBin(x);
    if (bin > 0 && bin <= hData->GetNbinsX()) {
      Double_t oldContent = hData->GetBinContent(bin);
      Double_t oldError = hData->GetBinError(bin);
      if (oldContent > 0) {
        Double_t weight1 = 1.0 / (oldError * oldError);
        Double_t weight2 = 1.0 / (errY * errY);
        Double_t newContent = (oldContent * weight1 + y * weight2) / (weight1 + weight2);
        Double_t newError = 1.0 / TMath::Sqrt(weight1 + weight2);
        hData->SetBinContent(bin, newContent);
        hData->SetBinError(bin, newError);
      } else {
        hData->SetBinContent(bin, y);
        hData->SetBinError(bin, errY);
      }
    }
  }
  
  for (int i = 0; i < gMC->GetN(); i++) {
    Double_t x, y;
    gMC->GetPoint(i, x, y);
    Double_t errY = gMC->GetErrorY(i);
    int bin = hMC->FindBin(x);
    if (bin > 0 && bin <= hMC->GetNbinsX()) {
      Double_t oldContent = hMC->GetBinContent(bin);
      Double_t oldError = hMC->GetBinError(bin);
      if (oldContent > 0) {
        Double_t weight1 = 1.0 / (oldError * oldError);
        Double_t weight2 = 1.0 / (errY * errY);
        Double_t newContent = (oldContent * weight1 + y * weight2) / (weight1 + weight2);
        Double_t newError = 1.0 / TMath::Sqrt(weight1 + weight2);
        hMC->SetBinContent(bin, newContent);
        hMC->SetBinError(bin, newError);
      } else {
        hMC->SetBinContent(bin, y);
        hMC->SetBinError(bin, errY);
      }
    }
  }
  
  // Create ratio graph
  TGraphErrors* gRatio = new TGraphErrors();
  gRatio->SetName("sigmaVsPtRatio");
  gRatio->SetTitle("Data/MC Ratio; p_{T} (GeV/c); #sigma_{Data} / #sigma_{MC}");
  
  int nRatioPoints = 0;
  for (int i = 1; i <= nCommonBins; i++) {
    Double_t yData = hData->GetBinContent(i);
    Double_t errData = hData->GetBinError(i);
    Double_t yMC = hMC->GetBinContent(i);
    Double_t errMC = hMC->GetBinError(i);
    
    if (yMC > 0 && yData > 0) {
      Double_t xCenter = hData->GetXaxis()->GetBinCenter(i);
      Double_t binWidth = hData->GetXaxis()->GetBinWidth(i);
      
      Double_t ratio = yData / yMC;
      Double_t errRatio = ratio * TMath::Sqrt(
        TMath::Power(errData / yData, 2) + TMath::Power(errMC / yMC, 2)
      );
      
      gRatio->SetPoint(nRatioPoints, xCenter, ratio);
      gRatio->SetPointError(nRatioPoints, binWidth / 2.0, errRatio);
      nRatioPoints++;
    }
  }
  
  delete hData;
  delete hMC;
  
  // Create Filipad2 plot
  Int_t nn = 0;
  TString padName = TString(outputFileName);
  padName.ReplaceAll(".root", "");
  if (tag && std::strlen(tag) > 0) {
    padName += tag;
  }
  Filipad2* pad = new Filipad2(padName, ++nn, 2.0, 0.4, 100, 50, 0.7, 1, 1);
  pad->Draw();
  
  TPad* toppad = pad->GetPad(1);
  toppad->SetLeftMargin(0.15);
  toppad->SetRightMargin(0.03);
  toppad->SetTopMargin(0.02);
  toppad->SetBottomMargin(0.0015);
  
  TPad* ratiopad = pad->GetPad(2);
  ratiopad->SetLeftMargin(0.15);
  ratiopad->SetRightMargin(0.03);
  ratiopad->SetTopMargin(0.0015);
  ratiopad->SetBottomMargin(0.24);
  
  // Top pad: Data and MC overlay
  toppad->cd();
  gMC->SetMarkerColor(kBlue);
  gMC->SetLineColor(kBlue);
  gMC->SetMarkerStyle(20);
  gMC->SetMarkerSize(0.8);
  gMC->SetTitle(Form(";#it{p}_{T} (GeV/#it{c});%s", yAxisTitle));
  
  Double_t maxMC = 0, maxData = 0;
  for (int i = 0; i < gMC->GetN(); i++) {
    Double_t x, y;
    gMC->GetPoint(i, x, y);
    if (y > maxMC) maxMC = y;
  }
  for (int i = 0; i < gData->GetN(); i++) {
    Double_t x, y;
    gData->GetPoint(i, x, y);
    if (y > maxData) maxData = y;
  }
  gMC->GetYaxis()->SetRangeUser(0, TMath::Max(maxMC, maxData) * 1.2);
  gMC->Draw("AP");
  
  gData->SetMarkerColor(kRed);
  gData->SetLineColor(kRed);
  gData->SetMarkerStyle(20);
  gData->SetMarkerSize(0.8);
  gData->Draw("P SAME");
  
  TLegend* leg = new TLegend(0.2, 0.7, 0.35, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gData, "Data", "lp");
  leg->AddEntry(gMC, "MC", "lp");
  leg->Draw();
  
  // Set axis labels and sizes (yield pad)
  gMC->GetXaxis()->SetTitleSize(0.068);
  gMC->GetYaxis()->SetTitleSize(0.055);
  gMC->GetXaxis()->SetLabelSize(0.058);
  gMC->GetYaxis()->SetLabelSize(0.050);
  gMC->GetXaxis()->SetTitleOffset(1.12);
  gMC->GetYaxis()->SetTitleOffset(1.3);
  
  // Bottom pad: Ratio
  ratiopad->cd();
  gRatio->SetMarkerColor(kBlack);
  gRatio->SetLineColor(kBlack);
  gRatio->SetMarkerStyle(20);
  gRatio->SetMarkerSize(0.8);
  gRatio->SetTitle(";#it{p}_{T} (GeV/#it{c});Data / MC");
  gRatio->GetYaxis()->SetRangeUser(0.45, 2.05);
  gRatio->Draw("AP");
  
  TLine* line1 = new TLine(gRatio->GetXaxis()->GetXmin(), 1.0, 
                           gRatio->GetXaxis()->GetXmax(), 1.0);
  line1->SetLineColor(kRed);
  line1->SetLineStyle(2);
  line1->SetLineWidth(1);
  line1->Draw("SAME");
  
  // Set axis labels and sizes for ratio pad
  gRatio->GetXaxis()->SetTitleSize(0.068);
  gRatio->GetYaxis()->SetTitleSize(0.075);
  gRatio->GetXaxis()->SetLabelSize(0.08);
  gRatio->GetYaxis()->SetLabelSize(0.075);
  gRatio->GetXaxis()->SetTitleOffset(1.12);
  gRatio->GetYaxis()->SetTitleOffset(0.90);
  
  // Add text info on top pad
  TLatex* lat = new TLatex();
  lat->SetNDC();
  lat->SetTextSize(0.04);
  TString fileName = TString(outputFileName);
  fileName = fileName(fileName.Last('/')+1, fileName.Length());
  toppad->cd();
  lat->DrawLatex(0.17, 0.92, Form("%s", fileName.Data()));
  
  TString pdfName = TString(outputFileName);
  pdfName.ReplaceAll(".root", "");
  if (tag && std::strlen(tag) > 0) {
    pdfName += tag;
  }
  pdfName += ".pdf";
  pad->C->SaveAs(Form("plots/%s", pdfName.Data()));
  std::cout << "Plot saved to: plots/" << pdfName.Data() << std::endl;
  
  delete gRatio;
}

  // Main function to create Q/pT correction file
void CreateQoverPtGraphsFile(
  const char* dataFile = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/610951_AnalysisResults.root",
  const char* mbMCFile = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/485159_AnalysisResults.root",
  const char* jjMCFile = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/611764_AnalysisResults.root",
  // const char* outputFile = "Sigma1overPt_Data_LHC22o_pass7_globalTracks_MC_LHC24f3c.root",
  const char* MBMC_outputFile = "Sigma1overPt_Data_LHC22o_pass7_globalTracks_MC_LHC24f3c.root",
  const char* JJMC_outputFile = "Sigma1overPt_Data_LHC22o_pass7_globalTracks_MC_LHC25a2b.root",
  const char* directory = "track-efficiency"
) {
  
  std::cout << "=== Creating Q/pT Correction Graphs File ===" << std::endl;
  std::cout << "Data file: " << dataFile << std::endl;
  std::cout << "MB MC file: " << mbMCFile << std::endl;
  std::cout << "JJ MC file: " << jjMCFile << std::endl;
  std::cout << "MB MC output file: " << MBMC_outputFile << std::endl;
  std::cout << "JJ MC output file: " << JJMC_outputFile << std::endl;
  
  // Extract graphs for CCDB (sigma(1/pT))
  std::cout << "\n=== Extracting Data graph (sigma(1/pT)) ===" << std::endl;
  TGraphErrors* gData1OverPt = ExtractSigma1OverPtVsPtFromCovMat(dataFile, directory);
  if (!gData1OverPt) {
    std::cerr << "Error: Failed to extract Data graph (sigma(1/pT))" << std::endl;
    return;
  }
  gData1OverPt->SetName("sigmaVsPtData");
  gData1OverPt->SetTitle("Data sigma(1/pT) vs p_{T}");

  std::cout << "\n=== Extracting MB MC graph (sigma(1/pT)) ===" << std::endl;
  TGraphErrors* gMBMC1OverPt = ExtractSigma1OverPtVsPtFromCovMatMC(mbMCFile, directory);
  if (!gMBMC1OverPt) {
    std::cerr << "Error: Failed to extract MB MC graph (sigma(1/pT))" << std::endl;
    delete gData1OverPt;
    return;
  }

  std::cout << "\n=== Extracting JJ MC graph (sigma(1/pT)) ===" << std::endl;
  TGraphErrors* gJJMC1OverPt = ExtractSigma1OverPtVsPtFromCovMatMC(jjMCFile, directory);
  if (!gJJMC1OverPt) {
    std::cerr << "Error: Failed to extract JJ MC graph (sigma(1/pT))" << std::endl;
    delete gData1OverPt;
    delete gMBMC1OverPt;
    return;
  }

  // Extract graphs for QA (sigma(pT)/pT)
  std::cout << "\n=== Extracting Data graph (sigma(pT)/pT) ===" << std::endl;
  TGraphErrors* gDataPtOverPt = ExtractSigmaPtOverPtVsPtFromCovMat(dataFile, directory);
  if (!gDataPtOverPt) {
    std::cerr << "Error: Failed to extract Data graph (sigma(pT)/pT)" << std::endl;
    delete gData1OverPt;
    delete gMBMC1OverPt;
    delete gJJMC1OverPt;
    return;
  }
  gDataPtOverPt->SetName("sigmaPtOverPtVsPtData");
  gDataPtOverPt->SetTitle("Data sigma(pT)/pT vs p_{T}");

  std::cout << "\n=== Extracting MB MC graph (sigma(pT)/pT) ===" << std::endl;
  TGraphErrors* gMBMCPtOverPt = ExtractSigmaPtOverPtVsPtFromCovMatMC(mbMCFile, directory);
  if (!gMBMCPtOverPt) {
    std::cerr << "Error: Failed to extract MB MC graph (sigma(pT)/pT)" << std::endl;
    delete gData1OverPt;
    delete gMBMC1OverPt;
    delete gJJMC1OverPt;
    delete gDataPtOverPt;
    return;
  }
  gMBMCPtOverPt->SetName("sigmaPtOverPtVsPtMc");

  std::cout << "\n=== Extracting JJ MC graph (sigma(pT)/pT) ===" << std::endl;
  TGraphErrors* gJJMCPtOverPt = ExtractSigmaPtOverPtVsPtFromCovMatMC(jjMCFile, directory);
  if (!gJJMCPtOverPt) {
    std::cerr << "Error: Failed to extract JJ MC graph (sigma(pT)/pT)" << std::endl;
    delete gData1OverPt;
    delete gMBMC1OverPt;
    delete gJJMC1OverPt;
    delete gDataPtOverPt;
    delete gMBMCPtOverPt;
    return;
  }
  gJJMCPtOverPt->SetName("sigmaPtOverPtVsPtMc");
  // Keep original name for plotting, but will rename when saving to file
  
  // Combine MB and JJ MC graphs (or use one based on pT range)
  // For now, we'll create separate files or combine them
  // Option 1: Use MB for low pT, JJ for high pT
  // Option 2: Create weighted average
  // Option 3: Create separate files
  
  // Create output directory
  gSystem->MakeDirectory("plots");

  // Quick ratio printouts before writing files
  PrintDataOverMcRatioSummary(gData1OverPt, gMBMC1OverPt, "Data / MB MC (sigma(1/pT), CCDB inputs)");
  PrintDataOverMcRatioSummary(gData1OverPt, gJJMC1OverPt, "Data / JJ MC (sigma(1/pT), CCDB inputs)");
  PrintDataOverMcRatioSummary(gDataPtOverPt, gMBMCPtOverPt, "Data / MB MC (sigma(pT)/pT, QA)");
  PrintDataOverMcRatioSummary(gDataPtOverPt, gJJMCPtOverPt, "Data / JJ MC (sigma(pT)/pT, QA)");
  
  // Create MB MC + Data file (only sigmaVsPtMc and sigmaVsPtData with sigma(1/pT))
  std::cout << "\n=== Creating MB MC file ===" << std::endl;
  TGraphErrors* gMBMC_renamed = (TGraphErrors*)gMBMC1OverPt->Clone("sigmaVsPtMc");
  TGraphErrors* gData_renamed_MB = (TGraphErrors*)gData1OverPt->Clone("sigmaVsPtData");
  TFile* outFileMB = TFile::Open(MBMC_outputFile, "RECREATE");
  TList* ccdb_object_MB = new TList();
  ccdb_object_MB->SetName("ccdb_object");
  ccdb_object_MB->Add(gMBMC_renamed);
  ccdb_object_MB->Add(gData_renamed_MB);
  outFileMB->WriteObject(ccdb_object_MB, "ccdb_object");
  outFileMB->Close();
  std::cout << "MB MC file created: " << MBMC_outputFile << std::endl;
  
  // Create plot for MB MC - sigma(pT)/pT
  TGraphErrors* gData_MB = (TGraphErrors*)gDataPtOverPt->Clone("gData_MB");
  TGraphErrors* gMBMC_plot = (TGraphErrors*)gMBMCPtOverPt->Clone("gMBMC_plot");
  CreateQoverPtPlot(gData_MB, gMBMC_plot, MBMC_outputFile, "#sigma(#it{p}_{T})/#it{p}_{T}", "_sigmaptoverpt");
  delete gData_MB;
  delete gMBMC_plot;

  // Create plot for MB MC - sigma(1/pT)
  TGraphErrors* gData_MB_1overpt = (TGraphErrors*)gData1OverPt->Clone("gData_MB_1overpt");
  TGraphErrors* gMBMC_plot_1overpt = (TGraphErrors*)gMBMC1OverPt->Clone("gMBMC_plot_1overpt");
  CreateQoverPtPlot(gData_MB_1overpt, gMBMC_plot_1overpt, MBMC_outputFile, "#sigma(1/#it{p}_{T}) (GeV^{-1}#it{c})", "_sigma1overpt");
  delete gData_MB_1overpt;
  delete gMBMC_plot_1overpt;
  
  // Create JJ MC + Data file (only sigmaVsPtMc and sigmaVsPtData with sigma(1/pT))
  std::cout << "\n=== Creating JJ MC file ===" << std::endl;
  TGraphErrors* gJJMC_renamed = (TGraphErrors*)gJJMC1OverPt->Clone("sigmaVsPtMc");
  TGraphErrors* gData_renamed_JJ = (TGraphErrors*)gData1OverPt->Clone("sigmaVsPtData");
  TFile* outFileJJ = TFile::Open(JJMC_outputFile, "RECREATE");
  TList* ccdb_object_JJ = new TList();
  ccdb_object_JJ->SetName("ccdb_object");
  ccdb_object_JJ->Add(gJJMC_renamed);
  ccdb_object_JJ->Add(gData_renamed_JJ);
  outFileJJ->WriteObject(ccdb_object_JJ, "ccdb_object");
  outFileJJ->Close();
  std::cout << "JJ MC file created: " << JJMC_outputFile << std::endl;
  
  // Create plot for JJ MC - sigma(pT)/pT
  TGraphErrors* gData_JJ = (TGraphErrors*)gDataPtOverPt->Clone("gData_JJ");
  TGraphErrors* gJJMC_plot = (TGraphErrors*)gJJMCPtOverPt->Clone("gJJMC_plot");
  CreateQoverPtPlot(gData_JJ, gJJMC_plot, JJMC_outputFile, "#sigma(#it{p}_{T})/#it{p}_{T}", "_sigmaptoverpt");
  delete gData_JJ;
  delete gJJMC_plot;

  // Create plot for JJ MC - sigma(1/pT)
  TGraphErrors* gData_JJ_1overpt = (TGraphErrors*)gData1OverPt->Clone("gData_JJ_1overpt");
  TGraphErrors* gJJMC_plot_1overpt = (TGraphErrors*)gJJMC1OverPt->Clone("gJJMC_plot_1overpt");
  CreateQoverPtPlot(gData_JJ_1overpt, gJJMC_plot_1overpt, JJMC_outputFile, "#sigma(1/#it{p}_{T}) (GeV^{-1}#it{c})", "_sigma1overpt");
  delete gData_JJ_1overpt;
  delete gJJMC_plot_1overpt;
  
  // Cleanup
  delete gMBMC1OverPt;
  delete gJJMC1OverPt;
  delete gData1OverPt;
  delete gMBMCPtOverPt;
  delete gJJMCPtOverPt;
  delete gDataPtOverPt;
  delete gMBMC_renamed;
  delete gJJMC_renamed;
  delete gData_renamed_MB;
  delete gData_renamed_JJ;
}
