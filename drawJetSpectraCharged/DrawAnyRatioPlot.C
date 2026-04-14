///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Draw Any Ratio Plot //////////
////////// author: Joonsuk Bae           ////////// 
////////// E-mail: jbae@cern.ch          //////////
////////// Last Modified: 2024           //////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

#include "../Filipad2.h"
#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TObject.h>
#include <TLegend.h>
#include <TLine.h>
#include <TSystem.h>
#include <TMath.h>
#include <TAxis.h>
#include <TString.h>
#include <iostream>
#include <vector>

// ============================================================
// Configuration Variables - Modify these to customize the plot
// ============================================================

// Input files and histograms
// Format: "directory/subdirectory/file.root" for file path
// Format: "directory/subdirectory/histogram_name" for histogram path in ROOT file
TString kNumeratorFile = "XsectionResults/Run3_CrossSection_data_498133_MC_516969.root";  // Numerator ROOT file path
TString kNumeratorHistPath = "Run3_CrossSection";  // Histogram path in numerator file (can include subdirectories: "dir/subdir/hist")

TString kDenominatorFile = "XsectionResults/Run3_CrossSection_data_498133_MC_515446.root";  // Denominator ROOT file path
TString kDenominatorHistPath = "Run3_CrossSection";  // Histogram path in denominator file (can include subdirectories: "dir/subdir/hist")

// Plot option: 1 = yield pad + ratio pad, 2 = ratio pad only
int kPlotOption = 1;

// Legend entries
TString kNumeratorLegend = "Run 3 data, (JJ MC, A tune)";      // Legend entry for numerator
TString kDenominatorLegend = "Run 3 data, (MB MC, A tune)";  // Legend entry for denominator

// Additional legends (optional - set empty string to disable)
// Yield pad additional legend
bool kUseYieldPadLegend = true;
double kYieldPadLegendX1 = 0.25;
double kYieldPadLegendY1 = 0.80;
double kYieldPadLegendX2 = 0.40;
double kYieldPadLegendY2 = 0.90;
TString kYieldPadLegendText = "ALICE WIP";  // Set text here if you want additional legend on yield pad

// Ratio pad additional legend
bool kUseRatioPadLegend = false;
double kRatioPadLegendX1 = 0.15;
double kRatioPadLegendY1 = 0.15;
double kRatioPadLegendX2 = 0.35;
double kRatioPadLegendY2 = 0.25;
TString kRatioPadLegendText = "";  // Set text here if you want additional legend on ratio pad

// Axis titles
TString kYieldPadXTitle = "#it{p}_{T} (GeV/#it{c})";
TString kYieldPadYTitle = "d^{2}#sigma/d#it{p}_{T}d#it{#eta}  [mb (GeV/#it{c})^{-1}]";
TString kRatioPadXTitle = "#it{p}_{T} (GeV/#it{c})";
TString kRatioPadYTitle = "Unfolded (JJ MC / MB MC)";

// Plot ranges
double kXMin = 0.0;
double kXMax = 140.0;
double kYieldYMin = 3e-7;
double kYieldYMax = 1e+3;
double kRatioYMin = 0.75;
double kRatioYMax = 1.05;

// Colors and styles
Color_t kNumeratorColor = kRed;
Color_t kDenominatorColor = kBlue;
int kNumeratorMarkerStyle = 20;
int kDenominatorMarkerStyle = 24;
double kMarkerSize = 0.8;
double kLineWidth = 2;

// Log scale options
bool kYieldLogX = false;
bool kYieldLogY = true;
bool kRatioLogX = false;
bool kRatioLogY = false;

// Output
TString kOutputName = "RatioPlot";  // Output file name (without extension)

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

void hoptset(TH1& h, double norm, Color_t color, 
             double xMin, double xMax, double yMin, double yMax,
             double markerSize, double lineWidth, int lineStyle, int markerStyle) {
  if (norm > 0) {
    h.Scale(norm, "width");
  }
  h.SetMarkerColor(color);
  h.SetLineColor(color);
  h.SetMarkerSize(markerSize);
  h.SetLineWidth(lineWidth);
  h.SetLineStyle(lineStyle);
  h.SetMarkerStyle(markerStyle);
  h.GetXaxis()->SetRangeUser(xMin, xMax);
  h.GetYaxis()->SetRangeUser(yMin, yMax);
  h.SetFillStyle(0);
  h.SetFillColor(0);
  h.SetStats(0);
}

void optFili(TPad& pad, int gridx, int gridy, int logx, int logy) {
  pad.SetGridx(gridx);
  pad.SetGridy(gridy);
  pad.SetLogx(logx);
  pad.SetLogy(logy);
}

// Helper function to get object type name
TString GetObjectType(TObject* obj) {
  if (!obj) return "Unknown";
  TString className = obj->ClassName();
  if (className.Contains("TH1D")) return "TH1D";
  if (className.Contains("TH1F")) return "TH1F";
  if (className.Contains("TH1I")) return "TH1I";
  if (className.Contains("TH1S")) return "TH1S";
  if (className.Contains("TH1C")) return "TH1C";
  if (className.Contains("TH1")) return "TH1";
  if (className.Contains("TGraphErrors")) return "TGraphErrors";
  if (className.Contains("TGraph")) return "TGraph";
  return className;
}

// Convert TGraph or TGraphErrors to TH1D
TH1D* GraphToHistogram(TGraph* graph) {
  if (!graph) return nullptr;
  
  // Find x-range from graph
  double xMin = 1e10;
  double xMax = -1e10;
  for (int i = 0; i < graph->GetN(); ++i) {
    double x, y;
    graph->GetPoint(i, x, y);
    if (x < xMin) xMin = x;
    if (x > xMax) xMax = x;
  }
  
  if (xMin >= xMax) return nullptr;
  
  // Create histogram with reasonable binning
  // Use 100 bins as default, or match graph points
  int nBins = TMath::Max(50, graph->GetN());
  TH1D* hist = new TH1D(Form("hist_from_graph_%p", graph), "", nBins, xMin, xMax);
  hist->SetDirectory(nullptr);
  
  for (int i = 0; i < graph->GetN(); ++i) {
    double x, y;
    graph->GetPoint(i, x, y);
    
    int bin = hist->FindBin(x);
    if (bin > 0 && bin <= hist->GetNbinsX()) {
      hist->SetBinContent(bin, y);
      
      // Set error if TGraphErrors
      TGraphErrors* grErr = dynamic_cast<TGraphErrors*>(graph);
      if (grErr) {
        double ey = grErr->GetErrorY(i);
        hist->SetBinError(bin, ey);
      } else {
        // For TGraph, use sqrt(y) as error estimate
        hist->SetBinError(bin, TMath::Sqrt(TMath::Abs(y)));
      }
    }
  }
  
  return hist;
}

// Check if two histograms have overlapping binning
struct BinningOverlap {
  bool hasOverlap;
  double xMin;
  double xMax;
  std::vector<int> numBins;  // Numerator bin indices in overlap region
  std::vector<int> denomBins; // Denominator bin indices in overlap region
};

BinningOverlap FindBinningOverlap(TH1* hNum, TH1* hDenom) {
  BinningOverlap result;
  result.hasOverlap = false;
  
  if (!hNum || !hDenom) return result;
  
  // Get bin edges
  TAxis* axisNum = hNum->GetXaxis();
  TAxis* axisDenom = hDenom->GetXaxis();
  
  double numXMin = axisNum->GetXmin();
  double numXMax = axisNum->GetXmax();
  double denomXMin = axisDenom->GetXmin();
  double denomXMax = axisDenom->GetXmax();
  
  // Find overlapping range
  double overlapMin = TMath::Max(numXMin, denomXMin);
  double overlapMax = TMath::Min(numXMax, denomXMax);
  
  if (overlapMin >= overlapMax) {
    std::cerr << "[Debug] No overlap: num range [" << numXMin << ", " << numXMax 
              << "], denom range [" << denomXMin << ", " << denomXMax << "]" << std::endl;
    return result;
  }
  
  result.hasOverlap = true;
  result.xMin = overlapMin;
  result.xMax = overlapMax;
  
  // Find bins in overlap region for numerator
  for (int i = 1; i <= hNum->GetNbinsX(); ++i) {
    double binLow = axisNum->GetBinLowEdge(i);
    double binUp = axisNum->GetBinUpEdge(i);
    if (binUp > overlapMin && binLow < overlapMax) {
      result.numBins.push_back(i);
    }
  }
  
  // Find bins in overlap region for denominator
  for (int i = 1; i <= hDenom->GetNbinsX(); ++i) {
    double binLow = axisDenom->GetBinLowEdge(i);
    double binUp = axisDenom->GetBinUpEdge(i);
    if (binUp > overlapMin && binLow < overlapMax) {
      result.denomBins.push_back(i);
    }
  }
  
  std::cerr << "[Debug] Overlap found: [" << overlapMin << ", " << overlapMax << "]" << std::endl;
  std::cerr << "[Debug] Numerator bins in overlap: " << result.numBins.size() << std::endl;
  std::cerr << "[Debug] Denominator bins in overlap: " << result.denomBins.size() << std::endl;
  
  return result;
}

// Calculate ratio histogram with binning matching
TH1* CalculateRatio(TH1* hNum, TH1* hDenom) {
  if (!hNum || !hDenom) return nullptr;
  
  // Check if binning is identical
  bool identicalBinning = true;
  if (hNum->GetNbinsX() != hDenom->GetNbinsX()) {
    identicalBinning = false;
  } else {
    // Check if bin edges match
    TAxis* axisNum = hNum->GetXaxis();
    TAxis* axisDenom = hDenom->GetXaxis();
    for (int i = 0; i <= hNum->GetNbinsX(); ++i) {
      if (TMath::Abs(axisNum->GetBinLowEdge(i) - axisDenom->GetBinLowEdge(i)) > 1e-6) {
        identicalBinning = false;
        break;
      }
    }
  }
  
  if (identicalBinning) {
    // Simple case: identical binning
    TH1* hRatio = (TH1*)hNum->Clone(Form("ratio_%s_%s", hNum->GetName(), hDenom->GetName()));
    hRatio->SetDirectory(nullptr);
    
    for (int i = 1; i <= hRatio->GetNbinsX(); ++i) {
      double valNum = hNum->GetBinContent(i);
      double errNum = hNum->GetBinError(i);
      double valDenom = hDenom->GetBinContent(i);
      double errDenom = hDenom->GetBinError(i);
      
      if (valDenom > 0) {
        double ratio = valNum / valDenom;
        double relErrNum = (valNum > 0) ? errNum / valNum : 0;
        double relErrDenom = errDenom / valDenom;
        double ratioErr = ratio * TMath::Sqrt(relErrNum * relErrNum + relErrDenom * relErrDenom);
        
        hRatio->SetBinContent(i, ratio);
        hRatio->SetBinError(i, ratioErr);
      } else {
        hRatio->SetBinContent(i, 0);
        hRatio->SetBinError(i, 0);
      }
    }
    
    return hRatio;
  } else {
    // Different binning: find overlap and calculate ratio in overlap region
    BinningOverlap overlap = FindBinningOverlap(hNum, hDenom);
    
    if (!overlap.hasOverlap || overlap.numBins.empty() || overlap.denomBins.empty()) {
      std::cerr << "[Error] No overlapping binning found between numerator and denominator" << std::endl;
      return nullptr;
    }
    
    // Create ratio histogram with numerator's binning in overlap region
    // Use numerator's binning as reference
    TH1* hRatio = (TH1*)hNum->Clone(Form("ratio_%s_%s", hNum->GetName(), hDenom->GetName()));
    hRatio->SetDirectory(nullptr);
    
    // Reset all bins
    for (int i = 0; i <= hRatio->GetNbinsX() + 1; ++i) {
      hRatio->SetBinContent(i, 0);
      hRatio->SetBinError(i, 0);
    }
    
    // Calculate ratio only in overlap region
    TAxis* axisNum = hNum->GetXaxis();
    TAxis* axisDenom = hDenom->GetXaxis();
    
    for (size_t idx = 0; idx < overlap.numBins.size(); ++idx) {
      int numBin = overlap.numBins[idx];
      double numBinCenter = axisNum->GetBinCenter(numBin);
      
      // Find corresponding bin in denominator
      int denomBin = axisDenom->FindBin(numBinCenter);
      if (denomBin < 1 || denomBin > hDenom->GetNbinsX()) continue;
      
      // Check if this denominator bin is in overlap region
      bool denomInOverlap = false;
      for (size_t j = 0; j < overlap.denomBins.size(); ++j) {
        if (overlap.denomBins[j] == denomBin) {
          denomInOverlap = true;
          break;
        }
      }
      if (!denomInOverlap) continue;
      
      double valNum = hNum->GetBinContent(numBin);
      double errNum = hNum->GetBinError(numBin);
      double valDenom = hDenom->GetBinContent(denomBin);
      double errDenom = hDenom->GetBinError(denomBin);
      
      if (valDenom > 0) {
        double ratio = valNum / valDenom;
        double relErrNum = (valNum > 0) ? errNum / valNum : 0;
        double relErrDenom = errDenom / valDenom;
        double ratioErr = ratio * TMath::Sqrt(relErrNum * relErrNum + relErrDenom * relErrDenom);
        
        hRatio->SetBinContent(numBin, ratio);
        hRatio->SetBinError(numBin, ratioErr);
      } else {
        hRatio->SetBinContent(numBin, 0);
        hRatio->SetBinError(numBin, 0);
      }
    }
    
    std::cerr << "[Info] Ratio calculated in overlap region [" << overlap.xMin 
              << ", " << overlap.xMax << "]" << std::endl;
    
    return hRatio;
  }
}

// ============================================================
// Main Function
// ============================================================

void DrawAnyRatioPlot() {
  std::cerr << "========================================" << std::endl;
  std::cerr << "Drawing Ratio Plot" << std::endl;
  std::cerr << "========================================" << std::endl;
  
  // 1. Load numerator histogram
  std::cerr << "[Info] Loading numerator from: " << kNumeratorFile.Data() << std::endl;
  std::cerr << "[Info] Histogram path: " << kNumeratorHistPath.Data() << std::endl;
  
  TFile* fNum = TFile::Open(kNumeratorFile.Data(), "READ");
  if (!fNum || fNum->IsZombie()) {
    std::cerr << "[Error] Cannot open numerator file: " << kNumeratorFile.Data() << std::endl;
    return;
  }
  
  TObject* objNum = fNum->Get(kNumeratorHistPath.Data());
  if (!objNum) {
    std::cerr << "[Error] Cannot find object '" << kNumeratorHistPath.Data() 
              << "' in numerator file: " << kNumeratorFile.Data() << std::endl;
    std::cerr << "[Info] Listing contents of numerator file:" << std::endl;
    fNum->ls();
    // fNum->Close();  // Keep file open for ROOT browser inspection
    return;
  }
  
  TH1* hNum = nullptr;
  TString numType = GetObjectType(objNum);
  
  // Try TH1 first
  hNum = dynamic_cast<TH1*>(objNum);
  if (!hNum) {
    // Try TGraph or TGraphErrors
    TGraph* grNum = dynamic_cast<TGraph*>(objNum);
    if (grNum) {
      std::cerr << "[Info] Converting " << numType << " to TH1D for numerator" << std::endl;
      hNum = GraphToHistogram(grNum);
      if (!hNum) {
        std::cerr << "[Error] Failed to convert " << numType << " to histogram" << std::endl;
        // fNum->Close();  // Keep file open for ROOT browser inspection
        return;
      }
    } else {
      std::cerr << "[Error] Object '" << kNumeratorHistPath.Data() 
                << "' is not a TH1, TGraph, or TGraphErrors (type: " << objNum->ClassName() << ")" << std::endl;
      // fNum->Close();  // Keep file open for ROOT browser inspection
      return;
    }
  }
  
  hNum->SetDirectory(nullptr);
  TH1* hNumClone = (TH1*)hNum->Clone("hNumerator_clone");
  hNumClone->SetDirectory(nullptr);
  std::cerr << "[Info] Loaded numerator: " << numType 
            << " -> " << GetObjectType(hNumClone) 
            << ", " << hNumClone->GetNbinsX() << " bins, range [" 
            << hNumClone->GetXaxis()->GetXmin() << ", " 
            << hNumClone->GetXaxis()->GetXmax() << "]" << std::endl;
  // fNum->Close();  // Keep file open for ROOT browser inspection
  
  // 2. Load denominator histogram
  std::cerr << "[Info] Loading denominator from: " << kDenominatorFile.Data() << std::endl;
  std::cerr << "[Info] Histogram path: " << kDenominatorHistPath.Data() << std::endl;
  
  TFile* fDenom = TFile::Open(kDenominatorFile.Data(), "READ");
  if (!fDenom || fDenom->IsZombie()) {
    std::cerr << "[Error] Cannot open denominator file: " << kDenominatorFile.Data() << std::endl;
    // delete hNumClone;  // Keep objects for ROOT browser inspection
    return;
  }
  
  TObject* objDenom = fDenom->Get(kDenominatorHistPath.Data());
  if (!objDenom) {
    std::cerr << "[Error] Cannot find object '" << kDenominatorHistPath.Data() 
              << "' in denominator file: " << kDenominatorFile.Data() << std::endl;
    std::cerr << "[Info] Listing contents of denominator file:" << std::endl;
    fDenom->ls();
    // fDenom->Close();  // Keep file open for ROOT browser inspection
    // delete hNumClone;  // Keep objects for ROOT browser inspection
    return;
  }
  
  TH1* hDenom = nullptr;
  TString denomType = GetObjectType(objDenom);
  
  // Try TH1 first
  hDenom = dynamic_cast<TH1*>(objDenom);
  if (!hDenom) {
    // Try TGraph or TGraphErrors
    TGraph* grDenom = dynamic_cast<TGraph*>(objDenom);
    if (grDenom) {
      std::cerr << "[Info] Converting " << denomType << " to TH1D for denominator" << std::endl;
      hDenom = GraphToHistogram(grDenom);
      if (!hDenom) {
        std::cerr << "[Error] Failed to convert " << denomType << " to histogram" << std::endl;
        // fDenom->Close();  // Keep file open for ROOT browser inspection
        // delete hNumClone;  // Keep objects for ROOT browser inspection
        return;
      }
    } else {
      std::cerr << "[Error] Object '" << kDenominatorHistPath.Data() 
                << "' is not a TH1, TGraph, or TGraphErrors (type: " << objDenom->ClassName() << ")" << std::endl;
      // fDenom->Close();  // Keep file open for ROOT browser inspection
      // delete hNumClone;  // Keep objects for ROOT browser inspection
      return;
    }
  }
  
  hDenom->SetDirectory(nullptr);
  TH1* hDenomClone = (TH1*)hDenom->Clone("hDenominator_clone");
  hDenomClone->SetDirectory(nullptr);
  std::cerr << "[Info] Loaded denominator: " << denomType 
            << " -> " << GetObjectType(hDenomClone) 
            << ", " << hDenomClone->GetNbinsX() << " bins, range [" 
            << hDenomClone->GetXaxis()->GetXmin() << ", " 
            << hDenomClone->GetXaxis()->GetXmax() << "]" << std::endl;
  // fDenom->Close();  // Keep file open for ROOT browser inspection
  
  // 3. Check binning and calculate ratio
  std::cerr << "[Info] Checking binning compatibility..." << std::endl;
  TH1* hRatio = CalculateRatio(hNumClone, hDenomClone);
  if (!hRatio) {
    std::cerr << "[Error] Failed to calculate ratio - no overlapping binning found" << std::endl;
    std::cerr << "[Debug] Numerator: " << hNumClone->GetNbinsX() << " bins, range [" 
              << hNumClone->GetXaxis()->GetXmin() << ", " << hNumClone->GetXaxis()->GetXmax() << "]" << std::endl;
    std::cerr << "[Debug] Denominator: " << hDenomClone->GetNbinsX() << " bins, range [" 
              << hDenomClone->GetXaxis()->GetXmin() << ", " << hDenomClone->GetXaxis()->GetXmax() << "]" << std::endl;
    // delete hNumClone;  // Keep objects for ROOT browser inspection
    // delete hDenomClone;  // Keep objects for ROOT browser inspection
    return;
  }
  std::cerr << "[Info] Ratio calculated successfully" << std::endl;
  
  // 5. Create canvas
  std::cerr << "[Info] Creating canvas with plot option " << kPlotOption << std::endl;
  static int nn = 0;
  Filipad2* pad = nullptr;
  TPad* yieldPad = nullptr;
  TPad* ratioPad = nullptr;
  
  if (kPlotOption == 1) {
    // Yield pad + ratio pad
    // Filipad2(const char* padName, int inID, float inRelSize, float inR, int inXOffset, int inYOffset, float inAspect, int ichop, int ichopPt)
    pad = new Filipad2("drawAnyRatioPlot", ++nn, 2.0f, 0.5f, 100, 50, 0.7f, 1, 1);
    if (!pad) {
      std::cerr << "[Error] Failed to create Filipad2" << std::endl;
      // delete hNumClone;  // Keep objects for ROOT browser inspection
      // delete hDenomClone;  // Keep objects for ROOT browser inspection
      // delete hRatio;  // Keep objects for ROOT browser inspection
      return;
    }
    pad->Draw();
    yieldPad = pad->GetPad(1);
    ratioPad = pad->GetPad(2);
    
    if (!yieldPad || !ratioPad) {
      std::cerr << "[Error] Failed to get pads from Filipad2" << std::endl;
      // delete pad;  // Keep objects for ROOT browser inspection
      // delete hNumClone;  // Keep objects for ROOT browser inspection
      // delete hDenomClone;  // Keep objects for ROOT browser inspection
      // delete hRatio;  // Keep objects for ROOT browser inspection
      return;
    }
    
    optFili(*yieldPad, 0, 0, kYieldLogX ? 1 : 0, kYieldLogY ? 1 : 0);
    optFili(*ratioPad, 0, 0, kRatioLogX ? 1 : 0, kRatioLogY ? 1 : 0);
  } else if (kPlotOption == 2) {
    // Ratio pad only - use single pad (ratio = 1.0 means full pad)
    pad = new Filipad2(kOutputName, ++nn, 2.0f, 1.0f, 100, 50, 0.7f, 1, 1);
    if (!pad) {
      std::cerr << "[Error] Failed to create Filipad2" << std::endl;
      // delete hNumClone;  // Keep objects for ROOT browser inspection
      // delete hDenomClone;  // Keep objects for ROOT browser inspection
      // delete hRatio;  // Keep objects for ROOT browser inspection
      return;
    }
    pad->Draw();
    ratioPad = pad->GetPad(1);
    
    if (!ratioPad) {
      std::cerr << "[Error] Failed to get pad from Filipad2" << std::endl;
      // delete pad;  // Keep objects for ROOT browser inspection
      // delete hNumClone;  // Keep objects for ROOT browser inspection
      // delete hDenomClone;  // Keep objects for ROOT browser inspection
      // delete hRatio;  // Keep objects for ROOT browser inspection
      return;
    }
    
    optFili(*ratioPad, 0, 0, kRatioLogX ? 1 : 0, kRatioLogY ? 1 : 0);
  } else {
    std::cerr << "[Error] Invalid plot option: " << kPlotOption << std::endl;
    // delete hNumClone;  // Keep objects for ROOT browser inspection
    // delete hDenomClone;  // Keep objects for ROOT browser inspection
    // delete hRatio;  // Keep objects for ROOT browser inspection
    return;
  }
  
  // 6. Draw yield pad (if option 1)
  if (kPlotOption == 1 && yieldPad) {
    yieldPad->cd();
    gPad->SetTicks(1, 1);
    
    // Draw denominator first (as reference)
    hoptset(*hDenomClone, 0, kDenominatorColor, kXMin, kXMax, kYieldYMin, kYieldYMax, 
            kMarkerSize, kLineWidth, 1, kDenominatorMarkerStyle);
    hset(*hDenomClone, kYieldPadXTitle, kYieldPadYTitle, 0.8, 1.0, 0.055, 0.055, 
         0.01, 0.005, 0.045, 0.040, 510, 1005);
    hDenomClone->GetYaxis()->SetTitleSize(0.065);
    hDenomClone->GetYaxis()->SetTitleOffset(1.0);
    hDenomClone->GetXaxis()->SetTitleSize(0.055);
    hDenomClone->GetXaxis()->SetTitleOffset(0.8);
    hDenomClone->GetXaxis()->SetLabelSize(0.055);
    hDenomClone->GetYaxis()->SetLabelSize(0.055);
    hDenomClone->Draw("lep");
    
    // Draw numerator
    hoptset(*hNumClone, 0, kNumeratorColor, kXMin, kXMax, kYieldYMin, kYieldYMax, 
            kMarkerSize, kLineWidth, 1, kNumeratorMarkerStyle);
    hNumClone->Draw("lep same");
    
    // Create legend for yield pad
    TLegend* legYield = new TLegend(0.50, 0.70, 0.85, 0.90, NULL, "brNDC");
    legYield->SetTextSize(0.040);
    legYield->SetBorderSize(0);
    legYield->SetFillColorAlpha(0, 0);
    legYield->SetTextFont(42);
    legYield->AddEntry(hNumClone, kNumeratorLegend.Data(), "lep");
    legYield->AddEntry(hDenomClone, kDenominatorLegend.Data(), "lep");
    legYield->Draw("same");
    
    // Additional legend on yield pad (if enabled)
    if (kUseYieldPadLegend && !kYieldPadLegendText.IsNull()) {
      TLegend* legYieldExtra = new TLegend(kYieldPadLegendX1, kYieldPadLegendY1, 
                                            kYieldPadLegendX2, kYieldPadLegendY2, NULL, "brNDC");
      legYieldExtra->SetTextSize(0.035);
      legYieldExtra->SetBorderSize(0);
      legYieldExtra->SetFillColorAlpha(0, 0);
      legYieldExtra->SetTextFont(42);
      legYieldExtra->AddEntry("", kYieldPadLegendText.Data(), "");
      legYieldExtra->Draw("same");
    }
  }
  
  // 7. Draw ratio pad
  if (ratioPad) {
    ratioPad->cd();
    gPad->SetTicks(1, 1);
    
    // Set up ratio histogram
    hoptset(*hRatio, 0, kNumeratorColor, kXMin, kXMax, kRatioYMin, kRatioYMax, 
            kMarkerSize, kLineWidth, 1, kNumeratorMarkerStyle);
    hset(*hRatio, kRatioPadXTitle, kRatioPadYTitle, 0.8, 1.0, 0.055, 0.055, 
         0.01, 0.005, 0.045, 0.040, 510, 1005);
    hRatio->GetXaxis()->SetTitleSize(0.065);
    hRatio->GetXaxis()->SetTitleOffset(1.0);
    hRatio->GetXaxis()->SetLabelSize(0.055);
    hRatio->GetYaxis()->SetTitleSize(0.065);
    hRatio->GetYaxis()->SetTitleOffset(1.0);
    hRatio->GetYaxis()->SetLabelSize(0.055);
    hRatio->Draw("lep");
    
    // Draw y=1 reference line
    TLine* line1 = new TLine(kXMin, 1.0, kXMax, 1.0);
    line1->SetLineColor(kBlack);
    line1->SetLineStyle(2);
    line1->SetLineWidth(2);
    line1->Draw("same");
    
    // Additional legend on ratio pad (if enabled)
    if (kUseRatioPadLegend && !kRatioPadLegendText.IsNull()) {
      TLegend* legRatioExtra = new TLegend(kRatioPadLegendX1, kRatioPadLegendY1, 
                                            kRatioPadLegendX2, kRatioPadLegendY2, NULL, "brNDC");
      legRatioExtra->SetTextSize(0.035);
      legRatioExtra->SetBorderSize(0);
      legRatioExtra->SetFillColorAlpha(0, 0);
      legRatioExtra->SetTextFont(42);
      legRatioExtra->AddEntry("", kRatioPadLegendText.Data(), "");
      legRatioExtra->Draw("same");
    }
  }
  
  // 8. Save
  pad->C->SaveAs(Form("%s.pdf", kOutputName.Data()));
  pad->C->SaveAs(Form("%s.root", kOutputName.Data()));
  
  std::cerr << "[Info] Plot saved: " << kOutputName.Data() << ".pdf" << std::endl;
  std::cerr << "========================================" << std::endl;
  std::cerr << "[Info] Objects kept in memory for ROOT browser inspection:" << std::endl;
  std::cerr << "  - hNumClone: numerator histogram" << std::endl;
  std::cerr << "  - hDenomClone: denominator histogram" << std::endl;
  std::cerr << "  - hRatio: ratio histogram" << std::endl;
  std::cerr << "  - pad: Filipad2 canvas" << std::endl;
  std::cerr << "  - fNum: numerator file (open)" << std::endl;
  std::cerr << "  - fDenom: denominator file (open)" << std::endl;
  
  // Cleanup - commented out to keep objects for ROOT browser inspection
  // delete hNumClone;
  // delete hDenomClone;
  // delete hRatio;
}
