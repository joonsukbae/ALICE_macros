///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Draw Systematic Uncertainty //////////
////////// author: Joonsuk Bae           ////////// 
////////// E-mail: jbae@cern.ch          //////////
////////// Last Modified: 2024           //////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

#include "../Filipad2.h"
#include <TROOT.h>
#include <TFile.h>
#include <TKey.h>
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
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <cmath>

// ============================================================
// Configuration Variables - Modify these to customize the plot
// ============================================================

// Default file path (using InvariantYieldResults instead of XsectionResults to avoid normalization errors)
TString kDefaultFile = "InvariantYieldResults/Run3_InvariantYield_data_497532_MC_596836.root"; // 2022 pass7 small, 25a2b tuner B, 1.5, temporary
// TString kDefaultFile = "InvariantYieldResults/Run3_InvariantYield_data_497532_MC_516969.root"; // 2022 pass7 small, 25a2b tuner A
// TString kDefaultFile = "InvariantYieldResults/Run3_InvariantYield_data_498133_MC_596836.root"; // 2022 pass7 full
TString kDefaultHistPath = "";  // Empty string - LoadHistogram will find first TH1 automatically (histogram name is "hist_XXX" from DrawJetsMCfJetQA.h)

// Source categories with variation file paths
// Each source can have one or more variations
// If constantUncertaintyPercent > 0, use constant value instead of files
struct SystematicSource {
  TString name;
  TString displayName;
  TString defaultLabel;  // Label for default configuration (e.g., "Default (* 1.5)")
  std::vector<TString> variationFiles;
  std::vector<TString> variationHistPaths;
  std::vector<TString> variationLabels;
  Color_t color;
  double constantUncertaintyPercent;  // If > 0, use this constant value (pT-independent)
};

std::vector<SystematicSource> kSystematicSources;

// Initialize systematic sources (placeholders - modify file paths as needed)
void InitializeSystematicSources() {
  kSystematicSources.clear();
  
  // 1. Track Efficiency
  SystematicSource trkEff;
  trkEff.name = "track_efficiency";
  trkEff.displayName = "Tracking efficiency";
  trkEff.defaultLabel = "Default (100%)";
  trkEff.variationFiles.push_back("InvariantYieldResults/Run3_InvariantYield_data_497532_MC_596832.root");
  trkEff.variationHistPaths.push_back("Run3_InvariantYield");
  trkEff.variationLabels.push_back("track efficiency - 99%");
  trkEff.color = kRed + 1;
  trkEff.constantUncertaintyPercent = 0.0;  // Use files, not constant
  kSystematicSources.push_back(trkEff);
  
  // 2. Track pT Resolution (2 variations)
  SystematicSource trkPtRes;
  trkPtRes.name = "track_pt_resolution";
  trkPtRes.displayName = "Track #it{p}_{T} resolution";
  trkPtRes.defaultLabel = "Default (* 1.5)";
  trkPtRes.variationFiles.push_back("InvariantYieldResults/Run3_InvariantYield_data_497532_MC_596838.root");
  trkPtRes.variationHistPaths.push_back("Run3_InvariantYield");
  trkPtRes.variationLabels.push_back("delta(Q/pT) * 1.2");
  trkPtRes.variationFiles.push_back("InvariantYieldResults/Run3_InvariantYield_data_497532_MC_596837.root");
  trkPtRes.variationHistPaths.push_back("Run3_InvariantYield");
  trkPtRes.variationLabels.push_back("delta(Q/pT) * 1.8");
  trkPtRes.color = kBlue + 1;
  trkPtRes.constantUncertaintyPercent = 0.0;  // Use files, not constant
  kSystematicSources.push_back(trkPtRes);
  
  // 3. Ambiguous Track (3 variations)
  SystematicSource ambTrk;
  ambTrk.name = "ambiguous_track";
  ambTrk.displayName = "Ambiguous track";
  ambTrk.defaultLabel = "Default (TimeMargin 500 ns)";
  ambTrk.variationFiles.push_back("InvariantYieldResults/Run3_InvariantYield_data_597125_MC_596833.root");
  ambTrk.variationHistPaths.push_back("Run3_InvariantYield");
  ambTrk.variationLabels.push_back("TimeMargin 125 ns (5 BC)");
  ambTrk.variationFiles.push_back("InvariantYieldResults/Run3_InvariantYield_data_597127_MC_596835.root");
  ambTrk.variationHistPaths.push_back("Run3_InvariantYield");
  ambTrk.variationLabels.push_back("TimeMargin 250 ns (10 BC)");
  ambTrk.variationFiles.push_back("InvariantYieldResults/Run3_InvariantYield_data_597126_MC_596834.root");
  ambTrk.variationHistPaths.push_back("Run3_InvariantYield");
  ambTrk.variationLabels.push_back("TimeMargin 1000 ns (40 BC)");
  ambTrk.color = kGreen + 2;
  ambTrk.constantUncertaintyPercent = 0.0;  // Use files, not constant
  kSystematicSources.push_back(ambTrk);
  
  // 4. Secondary Contamination (2 variations)
  SystematicSource secCont;
  secCont.name = "secondary_contamination";
  secCont.displayName = "Secondary contamination";
  secCont.defaultLabel = "Default (DCAz 2 cm)";
  // secCont.variationFiles.push_back("InvariantYieldResults/Run3_InvariantYield_data_619522_MC_619511.root");
  // secCont.variationHistPaths.push_back("Run3_InvariantYield");
  // secCont.variationLabels.push_back("DCAz 0.1 cm");
  secCont.variationFiles.push_back("InvariantYieldResults/Run3_InvariantYield_data_619713_MC_619512.root");
  secCont.variationHistPaths.push_back("Run3_InvariantYield");
  secCont.variationLabels.push_back("DCAz 0.2 cm");
  secCont.variationFiles.push_back("InvariantYieldResults/Run3_InvariantYield_data_619714_MC_619513.root");
  secCont.variationHistPaths.push_back("Run3_InvariantYield");
  secCont.variationLabels.push_back("DCAz 0.5 cm");
  secCont.color = kOrange + 7;
  secCont.constantUncertaintyPercent = 0.0;  // Use files, not constant
  kSystematicSources.push_back(secCont);
  
  // 5. Unfolding (loaded from DrawUnfoldingSystematicUncertainty.C output)
  SystematicSource unfolding;
  unfolding.name = "unfolding";
  unfolding.displayName = "Unfolding";
  unfolding.defaultLabel = "Default (SVD, K=12)";
  // Load from unfolding systematic uncertainty file
  TString unfoldingFile = "FinalSystematicUncertainty/UnfoldingSystematicUncertainty.root";
  TFile* unfoldingFileObj = TFile::Open(unfoldingFile.Data(), "READ");
  if (unfoldingFileObj && !unfoldingFileObj->IsZombie()) {
    TH1* hTotalUnfoldingUnc = (TH1*)unfoldingFileObj->Get("Total_unfolding_uncertainty");
    if (hTotalUnfoldingUnc) {
      // Use total unfolding uncertainty as a single variation
      // Save it temporarily and use constant uncertainty approach
      unfolding.constantUncertaintyPercent = 0.0;  // Will be loaded from file
      unfolding.variationFiles.push_back(unfoldingFile);  // Store file path
      unfolding.variationHistPaths.push_back("Total_unfolding_uncertainty");
      unfolding.variationLabels.push_back("Total unfolding uncertainty");
      std::cerr << "[Info] Unfolding uncertainty will be loaded from: " << unfoldingFile.Data() << std::endl;
    } else {
      std::cerr << "[Warning] Cannot find Total_unfolding_uncertainty in unfolding file, using placeholder" << std::endl;
      unfolding.constantUncertaintyPercent = 0.0;  // Placeholder
    }
    unfoldingFileObj->Close();
  } else {
    std::cerr << "[Warning] Cannot open unfolding uncertainty file, using placeholder" << std::endl;
    unfolding.constantUncertaintyPercent = 0.0;  // Placeholder
  }
  unfolding.color = kMagenta + 2;
  kSystematicSources.push_back(unfolding);
  
  // 6. Normalization (constant 10% uncertainty)
  SystematicSource norm;
  norm.name = "normalization";
  norm.displayName = "Normalization";
  norm.defaultLabel = "Default";
  norm.constantUncertaintyPercent = 5.0;  // pT-independent constant 10.3 -> 5%
  norm.color = kCyan + 2;
  kSystematicSources.push_back(norm);

  //  // 7. Event selection (1 variation)
  //  SystematicSource evtSel;
  //  evtSel.name = "event_selection";
  //  evtSel.displayName = "Event selection";
  //  evtSel.defaultLabel = "Default (sel8 / selMC)";
  //  evtSel.variationFiles.push_back("InvariantYieldResults/Run3_InvariantYield_data_619715_MC_619514.root");
  //  evtSel.variationHistPaths.push_back("Run3_InvariantYield");
  //  evtSel.variationLabels.push_back("sel8Full / selMCFull");
  //  evtSel.color = kOrange + 2;
  //  evtSel.constantUncertaintyPercent = 0.0;  // Use files, not constant
  //  kSystematicSources.push_back(evtSel);
}

// Plot ranges
double kXMin = 5.0;
double kXMax = 200.0;
double kYieldYMin = 1e-9;
double kYieldYMax = 1e-1;
double kRatioYMin = 0.8;
double kRatioYMax = 1.2;

// Axis titles
TString kYieldPadXTitle = "#it{p}_{T} (GeV/#it{c})";
TString kYieldPadYTitle = "1/N_{evt} d^{2}N/d#it{p}_{T}d#it{#eta}";  // Invariant Yield (no normalization error)
TString kRatioPadXTitle = "#it{p}_{T} (GeV/#it{c})";
TString kRatioPadYTitle = "Variation / Default";

// Output
TString kOutputDir = "FinalSystematicUncertainty";  // Output directory
TString kOutputName = "SystematicUncertainty";  // Output file name prefix

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
  
  double xMin = 1e10;
  double xMax = -1e10;
  for (int i = 0; i < graph->GetN(); ++i) {
    double x, y;
    graph->GetPoint(i, x, y);
    if (x < xMin) xMin = x;
    if (x > xMax) xMax = x;
  }
  
  if (xMin >= xMax) return nullptr;
  
  int nBins = TMath::Max(50, graph->GetN());
  TH1D* hist = new TH1D(Form("hist_from_graph_%p", graph), "", nBins, xMin, xMax);
  hist->SetDirectory(nullptr);
  
  for (int i = 0; i < graph->GetN(); ++i) {
    double x, y;
    graph->GetPoint(i, x, y);
    
    int bin = hist->FindBin(x);
    if (bin > 0 && bin <= hist->GetNbinsX()) {
      hist->SetBinContent(bin, y);
      
      TGraphErrors* grErr = dynamic_cast<TGraphErrors*>(graph);
      if (grErr) {
        double ey = grErr->GetErrorY(i);
        hist->SetBinError(bin, ey);
      } else {
        hist->SetBinError(bin, TMath::Sqrt(TMath::Abs(y)));
      }
    }
  }
  
  return hist;
}

// Load histogram from file
TH1* LoadHistogram(TString fileName, TString histPath) {
  TFile* file = TFile::Open(fileName.Data(), "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "[Error] Cannot open file: " << fileName.Data() << std::endl;
    return nullptr;
  }
  
  TObject* obj = nullptr;
  
  // If histPath is empty, find first TH1 automatically
  if (histPath.IsNull() || histPath == "") {
    std::cerr << "[Info] histPath is empty, searching for first TH1 in file..." << std::endl;
    TList* keys = file->GetListOfKeys();
    bool found = false;
    for (int i = 0; i < keys->GetSize(); ++i) {
      TKey* key = (TKey*)keys->At(i);
      TString className = key->GetClassName();
      if (className.Contains("TH1")) {
        obj = file->Get(key->GetName());
        if (obj) {
          std::cerr << "[Info] Found histogram: " << key->GetName() << " (" << className << ")" << std::endl;
          found = true;
          break;
        }
      }
    }
    if (!found) {
      std::cerr << "[Error] Cannot find any TH1 histogram in file: " << fileName.Data() << std::endl;
      // file->Close();  // Keep file open for ROOT browser inspection
      return nullptr;
    }
  } else {
    // Try to get histogram by name
    obj = file->Get(histPath.Data());
    if (!obj) {
      std::cerr << "[Warning] Cannot find object '" << histPath.Data() 
                << "' in file: " << fileName.Data() << std::endl;
      std::cerr << "[Info] Trying to find first TH1 in file..." << std::endl;
      // Try to find first TH1 in the file
      TList* keys = file->GetListOfKeys();
      bool found = false;
      for (int i = 0; i < keys->GetSize(); ++i) {
        TKey* key = (TKey*)keys->At(i);
        TString className = key->GetClassName();
        if (className.Contains("TH1")) {
          obj = file->Get(key->GetName());
          if (obj) {
            std::cerr << "[Info] Found histogram: " << key->GetName() << " (" << className << ")" << std::endl;
            found = true;
            break;
          }
        }
      }
      if (!found) {
        std::cerr << "[Error] Cannot find any TH1 histogram in file: " << fileName.Data() << std::endl;
        // file->Close();  // Keep file open for ROOT browser inspection
        return nullptr;
      }
    }
  }
  
  TH1* hist = nullptr;
  TString objType = GetObjectType(obj);
  
  hist = dynamic_cast<TH1*>(obj);
  if (!hist) {
    TGraph* gr = dynamic_cast<TGraph*>(obj);
    if (gr) {
      std::cerr << "[Info] Converting " << objType << " to TH1D" << std::endl;
      hist = GraphToHistogram(gr);
      if (!hist) {
        std::cerr << "[Error] Failed to convert " << objType << " to histogram" << std::endl;
        // file->Close();  // Keep file open for ROOT browser inspection
        return nullptr;
      }
    } else {
      std::cerr << "[Error] Object '" << histPath.Data() 
                << "' is not a TH1, TGraph, or TGraphErrors (type: " << obj->ClassName() << ")" << std::endl;
      // file->Close();  // Keep file open for ROOT browser inspection
      return nullptr;
    }
  }
  
  hist->SetDirectory(nullptr);
  TH1* histClone = (TH1*)hist->Clone(Form("hist_%s_%p", histPath.Data(), hist));
  histClone->SetDirectory(nullptr);
  
  // file->Close();  // Keep file open for ROOT browser inspection
  return histClone;
}

// Calculate ratio histogram
TH1* CalculateRatio(TH1* hNum, TH1* hDenom) {
  if (!hNum || !hDenom) return nullptr;
  
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
}

// Calculate relative uncertainty from ratio (|ratio - 1| * 100)
TH1* CalculateRelativeUncertainty(TH1* hRatio) {
  if (!hRatio) return nullptr;
  
  TH1* hRelUnc = (TH1*)hRatio->Clone(Form("relUnc_%s", hRatio->GetName()));
  hRelUnc->SetDirectory(nullptr);
  
  for (int i = 1; i <= hRelUnc->GetNbinsX(); ++i) {
    double ratio = hRatio->GetBinContent(i);
    double ratioErr = hRatio->GetBinError(i);
    
    double relUnc = TMath::Abs(ratio - 1.0) * 100.0;  // Convert to percentage
    double relUncErr = ratioErr * 100.0;  // Error in percentage
    
    hRelUnc->SetBinContent(i, relUnc);
    hRelUnc->SetBinError(i, relUncErr);
  }
  
  return hRelUnc;
}

// Calculate quadrature sum of systematic uncertainties
void CalculateQuadratureSum(const std::vector<TH1*>& hSysts, TH1* hSystResult) {
  if (hSysts.empty() || !hSystResult) return;
  
  int nBins = hSysts[0]->GetNbinsX();
  
  for (int i = 1; i <= nBins; ++i) {
    double sumOfSquares = 0.0;
    
    for (const auto& hist : hSysts) {
      double value = hist->GetBinContent(i);  // Already in percentage
      sumOfSquares += value * value;
    }
    
    double quadratureSum = TMath::Sqrt(sumOfSquares);
    hSystResult->SetBinContent(i, quadratureSum);
    hSystResult->SetBinError(i, 0);  // No error on total
  }
  
  // Smoothing using 3-bin moving average
  TH1* hSystSmoothed = (TH1*)hSysts[0]->Clone("hSystSmoothed");
  
  for (int i = 1; i <= nBins; ++i) {
    double smoothedValue = 0.0;
    
    if (i == 1) {
      smoothedValue = (hSystResult->GetBinContent(i) + hSystResult->GetBinContent(i + 1)) / 2.0;
    } else if (i == nBins) {
      smoothedValue = (hSystResult->GetBinContent(i) + hSystResult->GetBinContent(i - 1)) / 2.0;
    } else {
      smoothedValue = (hSystResult->GetBinContent(i - 1) +
                       hSystResult->GetBinContent(i) +
                       hSystResult->GetBinContent(i + 1)) / 3.0;
    }
    
    hSystSmoothed->SetBinContent(i, smoothedValue);
  }
  
  for (int i = 1; i <= nBins; ++i) {
    hSystResult->SetBinContent(i, hSystSmoothed->GetBinContent(i));
  }
  
  // delete hSystSmoothed;  // Keep for ROOT browser inspection
}

// ============================================================
// Main Function
// ============================================================

void DrawSystematicUncertainty() {
  std::cerr << "========================================" << std::endl;
  std::cerr << "Drawing Systematic Uncertainty" << std::endl;
  std::cerr << "========================================" << std::endl;
  
  // Create output directory at the beginning
  gSystem->MakeDirectory(kOutputDir.Data());
  
  // Initialize systematic sources
  InitializeSystematicSources();
  
  // Load default histogram
  std::cerr << "[Info] Loading default from: " << kDefaultFile.Data() << std::endl;
  TH1* hDefault = LoadHistogram(kDefaultFile, kDefaultHistPath);
  if (!hDefault) {
    std::cerr << "[Error] Failed to load default histogram" << std::endl;
    return;
  }
  
  static int nn = 0;
  std::vector<TH1*> hRelativeUncertainties;  // Store relative uncertainties for final plot
  std::vector<TString> sourceNames;
  std::vector<Color_t> sourceColors;
  
  // Process each systematic source
  for (size_t srcIdx = 0; srcIdx < kSystematicSources.size(); ++srcIdx) {
    const SystematicSource& source = kSystematicSources[srcIdx];
    std::cerr << "[Info] Processing source: " << source.displayName.Data() << std::endl;
    
    TH1* hSourceUnc = nullptr;
    
    // Special handling for unfolding: load from DrawUnfoldingSystematicUncertainty.C output
    if (source.name == "unfolding" && !source.variationFiles.empty()) {
      TString unfoldingFile = source.variationFiles[0];
      TString unfoldingHistPath = source.variationHistPaths[0];
      
      std::cerr << "[Info] Loading unfolding uncertainty from: " << unfoldingFile.Data() << std::endl;
      TH1* hUnfoldingUnc = LoadHistogram(unfoldingFile, unfoldingHistPath);
      
      if (hUnfoldingUnc) {
        // Rebin to match default binning if needed
        if (hUnfoldingUnc->GetNbinsX() != hDefault->GetNbinsX()) {
          std::cerr << "[Info] Rebinning unfolding uncertainty to match default binning" << std::endl;
          // Project to default binning
          TH1* hRebinned = (TH1*)hDefault->Clone(Form("hUnfoldingUnc_rebinned_%s", source.name.Data()));
          hRebinned->SetDirectory(nullptr);
          hRebinned->Reset();
          
          for (Int_t i = 1; i <= hRebinned->GetNbinsX(); ++i) {
            Double_t binCenter = hRebinned->GetBinCenter(i);
            Int_t srcBin = hUnfoldingUnc->GetXaxis()->FindBin(binCenter);
            if (srcBin >= 1 && srcBin <= hUnfoldingUnc->GetNbinsX()) {
              hRebinned->SetBinContent(i, hUnfoldingUnc->GetBinContent(srcBin));
              hRebinned->SetBinError(i, hUnfoldingUnc->GetBinError(srcBin));
            }
          }
          delete hUnfoldingUnc;
          hUnfoldingUnc = hRebinned;
        }
        hSourceUnc = hUnfoldingUnc;
        std::cerr << "[Info] Unfolding uncertainty loaded successfully" << std::endl;
      } else {
        std::cerr << "[Warning] Failed to load unfolding uncertainty, skipping..." << std::endl;
        continue;
      }
    }
    // Check if using constant uncertainty
    else if (source.constantUncertaintyPercent > 0) {
      // Create constant uncertainty histogram
      hSourceUnc = (TH1*)hDefault->Clone(Form("hSourceUnc_%s", source.name.Data()));
      hSourceUnc->SetDirectory(nullptr);
      hSourceUnc->Reset();
      
      int nBins = hSourceUnc->GetNbinsX();
      for (int i = 1; i <= nBins; ++i) {
        hSourceUnc->SetBinContent(i, source.constantUncertaintyPercent);
        hSourceUnc->SetBinError(i, 0);
      }
      
      std::cerr << "[Info] Using constant uncertainty: " << source.constantUncertaintyPercent << "%" << std::endl;
    } else {
      // Store all ratios for this source to find maximum deviation per bin
      std::vector<TH1*> sourceRatios;
      std::vector<TH1*> sourceVariations;  // Keep for plotting
      
      // Load all variations for this source
      for (size_t varIdx = 0; varIdx < source.variationFiles.size(); ++varIdx) {
        TString varFile = source.variationFiles[varIdx];
        TString varHistPath = source.variationHistPaths[varIdx];
        TString varLabel = source.variationLabels[varIdx];
        
        std::cerr << "[Info] Loading variation: " << varFile.Data() << std::endl;
        TH1* hVariation = LoadHistogram(varFile, varHistPath);
        if (!hVariation) {
          std::cerr << "[Warning] Failed to load variation, skipping..." << std::endl;
          continue;
        }
        
        sourceVariations.push_back(hVariation);
        
        // Calculate ratio
        TH1* hRatio = CalculateRatio(hVariation, hDefault);
        if (!hRatio) {
          std::cerr << "[Warning] Failed to calculate ratio, skipping..." << std::endl;
          continue;
        }
        
        sourceRatios.push_back(hRatio);
      }
      
      if (sourceRatios.empty()) {
        std::cerr << "[Warning] No valid ratios for source: " << source.displayName.Data() << ", skipping..." << std::endl;
        continue;
      }
      
      // Calculate uncertainty: for each bin, find maximum |ratio - 1| across all variations
      hSourceUnc = (TH1*)sourceRatios[0]->Clone(Form("hSourceUnc_%s", source.name.Data()));
      hSourceUnc->SetDirectory(nullptr);
      hSourceUnc->Reset();
      
      int nBins = hSourceUnc->GetNbinsX();
      for (int i = 1; i <= nBins; ++i) {
        double maxDeviation = 0.0;
        double maxDeviationErr = 0.0;
        
        // Find maximum |ratio - 1| across all variations for this bin
        for (size_t varIdx = 0; varIdx < sourceRatios.size(); ++varIdx) {
          TH1* hRatio = sourceRatios[varIdx];
          double ratio = hRatio->GetBinContent(i);
          double ratioErr = hRatio->GetBinError(i);
          
          // Calculate |ratio - 1| (deviation from unity)
          double deviation = TMath::Abs(ratio - 1.0);
          
          if (deviation > maxDeviation) {
            maxDeviation = deviation;
            maxDeviationErr = ratioErr;
          }
        }
        
        // Convert to percentage
        double relUnc = maxDeviation * 100.0;
        double relUncErr = maxDeviationErr * 100.0;
        
        hSourceUnc->SetBinContent(i, relUnc);
        hSourceUnc->SetBinError(i, relUncErr);
      }
      
      // Create one plot per source with all variations
      if (!sourceVariations.empty() && !sourceRatios.empty()) {
        TString padName = Form("Syst_%s", source.name.Data());
        padName.ReplaceAll(" ", "_");
        Filipad2* pad = new Filipad2(padName, ++nn, 2.0f, 0.5f, 100, 50, 0.7f, 1, 1);
        pad->Draw();
        
        TPad* yieldPad = pad->GetPad(1);
        TPad* ratioPad = pad->GetPad(2);
        
        if (yieldPad && ratioPad) {
          optFili(*yieldPad, 0, 0, 0, 1);  // logy = 1
          optFili(*ratioPad, 0, 0, 0, 0);  // logy = 0
          
          // Draw yield pad with all variations
          yieldPad->cd();
          gPad->SetTicks(1, 1);
          
          hoptset(*hDefault, 0, kBlack, kXMin, kXMax, kYieldYMin, kYieldYMax, 
                  0.8, 2, 1, 24);
          hset(*hDefault, kYieldPadXTitle, kYieldPadYTitle, 0.8, 1.0, 0.055, 0.055, 
               0.01, 0.005, 0.045, 0.040, 510, 1005);
          hDefault->Draw("lep");
          
          TLegend* legYield = new TLegend(0.50, 0.70, 0.85, 0.90, NULL, "brNDC");
          legYield->SetTextSize(0.040);
          legYield->SetBorderSize(0);
          legYield->SetFillColorAlpha(0, 0);
          legYield->SetTextFont(42);
          legYield->AddEntry(hDefault, source.defaultLabel.Data(), "lep");
          
          // Draw all variations
          // Use predefined color palette for variations to ensure valid colors
          // Colors are chosen to be distinct and within valid ROOT color range
          std::vector<Color_t> colorPalette = {
            kRed + 1,      // 0: Red
            kBlue + 1,     // 1: Blue
            kGreen + 2,    // 2: Green
            kMagenta + 2,  // 3: Magenta
            kOrange + 7,   // 4: Orange
            kCyan + 2,     // 5: Cyan
            kYellow + 1,   // 6: Yellow
            kPink + 1,     // 7: Pink
            kViolet + 2,   // 8: Violet
            kAzure + 2     // 9: Azure
          };
          
          // Start from source.color, then cycle through palette
          // If source.color is not in palette, start from index 0
          int baseColorIdx = 0;
          bool found = false;
          for (size_t i = 0; i < colorPalette.size(); ++i) {
            if (colorPalette[i] == source.color) {
              baseColorIdx = (i + 1) % colorPalette.size();
              found = true;
              break;
            }
          }
          // If source.color not found in palette, use interval-based approach
          // but ensure colors are within valid range (1-99 for ROOT)
          if (!found) {
            // Use source.color as base, but ensure valid range
            int startIdx = 0;
            for (size_t i = 0; i < colorPalette.size(); ++i) {
              if (TMath::Abs((int)colorPalette[i] - (int)source.color) < 3) {
                startIdx = (i + 1) % colorPalette.size();
                break;
              }
            }
            baseColorIdx = startIdx;
          }
          
          std::vector<Color_t> variationColors;
          for (size_t varIdx = 0; varIdx < sourceVariations.size(); ++varIdx) {
            int colorIdx = (baseColorIdx + varIdx) % colorPalette.size();
            variationColors.push_back(colorPalette[colorIdx]);
          }
          
          for (size_t varIdx = 0; varIdx < sourceVariations.size(); ++varIdx) {
            TH1* hVariation = sourceVariations[varIdx];
            Color_t varColor = variationColors[varIdx];
            
            hoptset(*hVariation, 0, varColor, kXMin, kXMax, kYieldYMin, kYieldYMax, 
                    0.8, 2, 1, 20 + varIdx);
            hVariation->Draw("lep same");
            legYield->AddEntry(hVariation, source.variationLabels[varIdx].Data(), "lep");
          }
          legYield->Draw("same");
          
          // Draw ratio pad with all variations
          ratioPad->cd();
          gPad->SetTicks(1, 1);
          
          TLegend* legRatio = new TLegend(0.50, 0.70, 0.85, 0.90, NULL, "brNDC");
          legRatio->SetTextSize(0.040);
          legRatio->SetBorderSize(0);
          legRatio->SetFillColorAlpha(0, 0);
          legRatio->SetTextFont(42);
          
          // Find actual min/max in all ratio histograms to adjust range if needed
          Double_t dataMin = 1e10;
          Double_t dataMax = -1e10;
          for (size_t varIdx = 0; varIdx < sourceRatios.size(); ++varIdx) {
            TH1* hRatio = sourceRatios[varIdx];
            for (Int_t i = 1; i <= hRatio->GetNbinsX(); ++i) {
              Double_t ratio = hRatio->GetBinContent(i);
              if (ratio > 0) {
                dataMin = TMath::Min(dataMin, ratio);
                dataMax = TMath::Max(dataMax, ratio);
              }
            }
          }
          // Adjust range if data exceeds default 0.8-1.2
          Double_t actualRatioMin = kRatioYMin;
          Double_t actualRatioMax = kRatioYMax;
          if (dataMin < kRatioYMin) actualRatioMin = TMath::Max(0.5, dataMin * 0.95);  // Add small margin
          if (dataMax > kRatioYMax) actualRatioMax = TMath::Min(2.0, dataMax * 1.05);  // Add small margin
          
          bool firstRatio = true;
          for (size_t varIdx = 0; varIdx < sourceRatios.size(); ++varIdx) {
            TH1* hRatio = sourceRatios[varIdx];
            Color_t varColor = (varIdx < variationColors.size()) ? variationColors[varIdx] : colorPalette[(baseColorIdx + varIdx) % colorPalette.size()];
            
            hoptset(*hRatio, 0, varColor, kXMin, kXMax, actualRatioMin, actualRatioMax, 
                    0.8, 2, 1, 20 + varIdx);
            hset(*hRatio, kRatioPadXTitle, kRatioPadYTitle, 0.8, 1.0, 0.055, 0.055, 
                 0.01, 0.005, 0.045, 0.040, 510, 1005);
            
            if (firstRatio) {
              hRatio->Draw("ep");  // Draw points only, no line connection
              firstRatio = false;
            } else {
              hRatio->Draw("ep same");  // Draw points only, no line connection
            }
            
            legRatio->AddEntry(hRatio, source.variationLabels[varIdx].Data(), "ep");
          }
          
          TLine* line1 = new TLine(kXMin, 1.0, kXMax, 1.0);
          line1->SetLineColor(kBlack);
          line1->SetLineStyle(2);
          line1->SetLineWidth(2);
          line1->Draw("same");
          
          legRatio->Draw("same");
          
          // Save source plot
          pad->C->SaveAs(Form("%s/%s_%s.pdf", kOutputDir.Data(), kOutputName.Data(), source.name.Data()));
          pad->C->SaveAs(Form("%s/%s_%s.root", kOutputDir.Data(), kOutputName.Data(), source.name.Data()));
        }
      }
    }
    
    // Store for final plot
    if (hSourceUnc) {
      hRelativeUncertainties.push_back(hSourceUnc);
      sourceNames.push_back(source.displayName);
      sourceColors.push_back(source.color);
    }
  }
  
  // Draw final systematic uncertainty plot
  if (hRelativeUncertainties.empty()) {
    std::cerr << "[Error] No relative uncertainties calculated" << std::endl;
    return;
  }
  
  std::cerr << "[Info] Drawing final systematic uncertainty plot" << std::endl;
  
  // Calculate total uncertainty
  TH1* hTotalUnc = (TH1*)hRelativeUncertainties[0]->Clone("hTotalSystematicUncertainty");
  hTotalUnc->Reset();
  CalculateQuadratureSum(hRelativeUncertainties, hTotalUnc);
  
  // Create canvas for final plot
  TCanvas* cFinal = new TCanvas("cSystematicUncertainty", "Systematic Uncertainty", 800, 600);
  cFinal->cd();
  gPad->SetTicks(1, 1);
  gPad->SetLogx(0);
  gPad->SetLogy(0);
  
  // Draw relative uncertainties
  const double kFinalYMin = 0.0;
  const double kFinalYMax = 30.0;

  bool firstHist = true;
  TLegend* legFinal = new TLegend(0.466165, 0.61913, 0.765664, 0.867826, NULL, "brNDC");
  legFinal->SetTextSize(0.040);
  legFinal->SetBorderSize(0);
  legFinal->SetFillColorAlpha(0, 0);
  legFinal->SetTextFont(42);
  
  for (size_t i = 0; i < hRelativeUncertainties.size(); ++i) {
    TH1* hist = hRelativeUncertainties[i];
    
    hist->SetLineWidth(5);
    hist->SetLineStyle(i + 1);
    hist->SetMarkerStyle(0);
    hist->SetFillStyle(0);
    hist->SetLineColor(sourceColors[i]);
    hist->GetXaxis()->SetRangeUser(kXMin, kXMax);
    hist->SetMinimum(kFinalYMin);
    hist->SetMaximum(kFinalYMax);
    hset(*hist, kYieldPadXTitle, "Rel. Uncertainties (%)", 1.3, 1.1, 0.05, 0.05, 0.01, 0.01,
         0.05, 0.06, 510, 505);
    
    if (firstHist) {
      hist->Draw("HIST");
      firstHist = false;
    } else {
      hist->Draw("HIST SAME");
    }
    
    legFinal->AddEntry(hist, sourceNames[i].Data(), "l");
  }
  
  // Draw total uncertainty
  hTotalUnc->SetLineWidth(6);
  hTotalUnc->SetLineStyle(1);
  hTotalUnc->SetLineColor(kBlack);
  hTotalUnc->SetMarkerStyle(0);
  hTotalUnc->SetFillStyle(0);
  hTotalUnc->GetXaxis()->SetRangeUser(kXMin, kXMax);
  hTotalUnc->SetMinimum(kFinalYMin);
  hTotalUnc->SetMaximum(kFinalYMax);
  hTotalUnc->Draw("HIST SAME");
  legFinal->AddEntry(hTotalUnc, "Total uncertainty", "l");
  
  legFinal->Draw();
  
  // Save final plot
  cFinal->SaveAs(Form("%s/%s_Final.pdf", kOutputDir.Data(), kOutputName.Data()));
  cFinal->SaveAs(Form("%s/%s_Final.root", kOutputDir.Data(), kOutputName.Data()));
  
  // Save systematic uncertainties to text file
  TString txtFileName = Form("%s/%s_Uncertainties.txt", kOutputDir.Data(), kOutputName.Data());
  std::ofstream txtFile(txtFileName.Data());
  if (!txtFile.is_open()) {
    std::cerr << "[Error] Cannot open file for writing: " << txtFileName.Data() << std::endl;
  } else {
    txtFile << "# Systematic Uncertainties (in %)" << std::endl;
    txtFile << "# Format: pT_low pT_high pT_center";
    for (size_t i = 0; i < sourceNames.size(); ++i) {
      txtFile << " " << sourceNames[i].Data();
    }
    txtFile << " Total" << std::endl;
    txtFile << std::fixed << std::setprecision(3);
    
    int nBins = hTotalUnc->GetNbinsX();
    for (int bin = 1; bin <= nBins; ++bin) {
      double pTLow = hTotalUnc->GetBinLowEdge(bin);
      double pTHigh = hTotalUnc->GetBinLowEdge(bin) + hTotalUnc->GetBinWidth(bin);
      double pTCenter = hTotalUnc->GetBinCenter(bin);
      
      // Only print bins within the plot range
      if (pTCenter < kXMin || pTCenter > kXMax) continue;
      
      txtFile << pTLow << " " << pTHigh << " " << pTCenter;
      
      // Write each source uncertainty
      for (size_t i = 0; i < hRelativeUncertainties.size(); ++i) {
        double unc = hRelativeUncertainties[i]->GetBinContent(bin);
        txtFile << " " << unc;
      }
      
      // Write total uncertainty
      double totalUnc = hTotalUnc->GetBinContent(bin);
      txtFile << " " << totalUnc << std::endl;
    }
    
    txtFile.close();
    std::cerr << "[Info] Systematic uncertainties saved to: " << txtFileName.Data() << std::endl;
  }
  
  std::cerr << "[Info] Systematic uncertainty plots saved" << std::endl;
  std::cerr << "========================================" << std::endl;
  std::cerr << "[Info] Objects kept in memory for ROOT browser inspection:" << std::endl;
  std::cerr << "  - hDefault: default histogram" << std::endl;
  std::cerr << "  - hRelativeUncertainties: vector of relative uncertainty histograms" << std::endl;
  std::cerr << "  - hTotalUnc: total systematic uncertainty" << std::endl;
}

// ============================================================
// Step-by-Step Workflow Function
// ============================================================
// This function runs the complete systematic uncertainty workflow:
// 1. MC Closure Test - determines optimal regularization parameters (d-vector for SVD, convergence for Bayes)
// 2. Unfolding Systematic Uncertainty - calculates unfolding-related systematics
// 3. DrawSystematicUncertainty - combines all sources into final uncertainty
//
// Usage: root -l -b -q 'DrawSystematicUncertainty.C("RunAll")'

void RunAllSystematicUncertainties() {
  std::cerr << "================================================================" << std::endl;
  std::cerr << "      COMPLETE SYSTEMATIC UNCERTAINTY WORKFLOW" << std::endl;
  std::cerr << "================================================================" << std::endl;
  std::cerr << std::endl;
  std::cerr << "This workflow runs the following steps in order:" << std::endl;
  std::cerr << "  Step 1: MC Closure Test (determines optimal SVD k and Bayes iter)" << std::endl;
  std::cerr << "  Step 2: Unfolding Systematic Uncertainty (method, iteration, prior dependence)" << std::endl;
  std::cerr << "  Step 3: Final Systematic Uncertainty (combines all sources)" << std::endl;
  std::cerr << std::endl;

  // Step 1: MC Closure Test
  std::cerr << "================================================================" << std::endl;
  std::cerr << "STEP 1: Running MC Closure Test" << std::endl;
  std::cerr << "================================================================" << std::endl;
  std::cerr << "Purpose: Determine optimal regularization parameters" << std::endl;
  std::cerr << "         SVD: d-vector method; Bayesian: convergence criterion" << std::endl;
  std::cerr << std::endl;

  // Load and execute DrawMcClosureTest
  gROOT->ProcessLine(".L DrawMcClosureTest.C");
  gROOT->ProcessLine("DrawMcClosureTest_JJMCOnly()");

  std::cerr << std::endl;
  std::cerr << "[Step 1 Complete] MC Closure Test finished" << std::endl;
  std::cerr << "Output: plots/MCClosureTest/MCClosure_*.root (contains optimalSVDk and optimalBayesIter)" << std::endl;
  std::cerr << std::endl;

  // Step 2: Unfolding Systematic Uncertainty
  std::cerr << "================================================================" << std::endl;
  std::cerr << "STEP 2: Running Unfolding Systematic Uncertainty" << std::endl;
  std::cerr << "================================================================" << std::endl;
  std::cerr << "Purpose: Calculate unfolding-related systematic uncertainties" << std::endl;
  std::cerr << "  - Method dependence (SVD vs Bayesian)" << std::endl;
  std::cerr << "  - Iteration dependence (k +/- 1)" << std::endl;
  std::cerr << "  - Prior dependence (POWHEG vs PYTHIA)" << std::endl;
  std::cerr << "  - MC closure non-closure" << std::endl;
  std::cerr << std::endl;

  // Load and execute DrawUnfoldingSystematicUncertainty
  gROOT->ProcessLine(".L DrawUnfoldingSystematicUncertainty.C");
  gROOT->ProcessLine("DrawUnfoldingSystematicUncertainty()");

  std::cerr << std::endl;
  std::cerr << "[Step 2 Complete] Unfolding Systematic Uncertainty finished" << std::endl;
  std::cerr << "Output: FinalSystematicUncertainty/UnfoldingSystematicUncertainty.root" << std::endl;
  std::cerr << std::endl;

  // Step 3: Final Systematic Uncertainty (all sources)
  std::cerr << "================================================================" << std::endl;
  std::cerr << "STEP 3: Running Final Systematic Uncertainty" << std::endl;
  std::cerr << "================================================================" << std::endl;
  std::cerr << "Purpose: Combine all systematic uncertainty sources" << std::endl;
  std::cerr << "  - Track efficiency" << std::endl;
  std::cerr << "  - Track pT resolution" << std::endl;
  std::cerr << "  - Ambiguous track" << std::endl;
  std::cerr << "  - Secondary contamination" << std::endl;
  std::cerr << "  - Unfolding (from Step 2)" << std::endl;
  std::cerr << "  - Normalization" << std::endl;
  std::cerr << std::endl;

  // Execute DrawSystematicUncertainty (already loaded)
  DrawSystematicUncertainty();

  std::cerr << std::endl;
  std::cerr << "[Step 3 Complete] Final Systematic Uncertainty finished" << std::endl;
  std::cerr << "Output: FinalSystematicUncertainty/SystematicUncertainty_Final.pdf" << std::endl;
  std::cerr << "        FinalSystematicUncertainty/SystematicUncertainty_Uncertainties.txt" << std::endl;
  std::cerr << std::endl;

  std::cerr << "================================================================" << std::endl;
  std::cerr << "      WORKFLOW COMPLETE" << std::endl;
  std::cerr << "================================================================" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Summary of outputs:" << std::endl;
  std::cerr << "  1. MC Closure Test:" << std::endl;
  std::cerr << "     - plots/MCClosureTest/MCClosure_*.root" << std::endl;
  std::cerr << "     - plots/MCClosureTest/<dataset>/Dvector.pdf" << std::endl;
  std::cerr << "     - plots/MCClosureTest/<dataset>/Bayes_deviation.pdf" << std::endl;
  std::cerr << std::endl;
  std::cerr << "  2. Unfolding Systematic Uncertainty:" << std::endl;
  std::cerr << "     - FinalSystematicUncertainty/UnfoldingSystematicUncertainty.root" << std::endl;
  std::cerr << "     - FinalSystematicUncertainty/UnfoldingSystematicUncertainty_*.pdf" << std::endl;
  std::cerr << "     - FinalSystematicUncertainty/UnfoldingSystematicUncertainty_*.txt" << std::endl;
  std::cerr << std::endl;
  std::cerr << "  3. Final Systematic Uncertainty:" << std::endl;
  std::cerr << "     - FinalSystematicUncertainty/SystematicUncertainty_Final.pdf" << std::endl;
  std::cerr << "     - FinalSystematicUncertainty/SystematicUncertainty_Uncertainties.txt" << std::endl;
  std::cerr << std::endl;
}
