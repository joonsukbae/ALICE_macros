#include "Filipad2.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TKey.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include <algorithm>
#include <cctype>
#include <iostream>
#include <vector>

using namespace std;

// Simple styling helpers (aligned with user's style)
void hset(TH1 &hid, TString xtit = "", TString ytit = "", double titoffx = 1.1,
          double titoffy = 1.2, double titsizex = 0.055, double titsizey = 0.055,
          double labeloffx = 0.01, double labeloffy = 0.005,
          double labelsizex = 0.05, double labelsizey = 0.05, int divx = 510,
          int divy = 510) {
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

void hoptset(TH1 &hid, Float_t N = 0, Color_t color = kBlack, Double_t minX = 0,
             Double_t maxX = 100, Double_t minY = 0, Double_t maxY = 1,
             Double_t MarkerSize = 1.0, Int_t markerStyle = 24, Int_t lineStyle = 1) {
  if (N != 0) {
    if (N == 1) {
      hid.Scale(1. / hid.Integral(), "width");
    } else if (N == 2) {
      hid.Scale(1. / hid.Integral(), "");
    } else {
      hid.Scale(1. / N, "width");
    }
  }
  hid.SetMarkerColor(color);
  hid.SetLineColor(color);
  hid.SetMarkerSize(MarkerSize);
  hid.SetMarkerStyle(markerStyle);
  hid.SetLineStyle(lineStyle);
}

void optFili(TPad &pid, Int_t gridx, Int_t gridy, Int_t logx, Int_t logy,
             Int_t logz = 0) {
  pid.SetGridy(gridx);
  pid.SetGridx(gridy);
  pid.SetLogx(logx);
  pid.SetLogy(logy);
  pid.SetLogz(logz);
}

static std::vector<int> ColorPalette = {
    kRed,         kBlue,       kGreen + 1, kMagenta,    kCyan + 1,
    kOrange + 1,  kYellow + 2, kAzure + 7, kViolet + 1, kSpring + 5,
    kPink + 1,    kTeal + 1,   kGray + 2,  kOrange + 3, kAzure + 2};

static inline TString LowerCopy(TString s) {
  std::string tmp = s.Data();
  std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c) { return std::tolower(c); });
  return TString(tmp.c_str());
}

// Check if histogram name looks like rho method related
static bool NameLooksLikeRhoMethod(const TString &name) {
  auto low = LowerCopy(name);
  return (low.Contains("h2_centrality_rhorandomconerandomtrackdirection") ||
          low.Contains("h2_centrality_rhorandomconewithoutleadingjet") ||
          low.Contains("h2_centrality_rhorandomconerandomtrackdirectionwithoutoneleadingjets") ||
          low.Contains("h2_centrality_rhorandomconerandomtrackdirectionwithouttwoleadingjets"));
}

// Get method name from histogram name
static TString GetMethodName(const TString &name) {
  auto low = LowerCopy(name);
  if (low.Contains("rhorandomconerandomtrackdirection") && !low.Contains("without")) 
    return "RC w/ Random Track Direction";
  if (low.Contains("rhorandomconewithoutleadingjet")) 
    return "RC w/o Leading Jet";
  if (low.Contains("rhorandomconerandomtrackdirectionwithoutoneleadingjets")) 
    return "RC w/ Random Track Direction w/o 1 Leading Jet";
  if (low.Contains("rhorandomconerandomtrackdirectionwithouttwoleadingjets")) 
    return "RC w/ Random Track Direction w/o 2 Leading Jets";
  return "Unknown";
}

// Project 2D histogram over centrality (X axis) to get MB
TH1 *ProjectRhoOverCentrality(TH2 *h2, const TString &outName) {
  if (!h2) return nullptr;
  
  // Project over full centrality range (X axis) to Y axis (rho)
  Int_t x1 = 1, x2 = h2->GetNbinsX();
  TH1 *h = h2->ProjectionY(outName, x1, x2, "e");
  return h;
}

void CollectRhoMethods(TDirectory *dir, std::vector<TH1 *> &outHists, 
                      std::vector<TString> &outLabels, std::vector<TString> &outMethods) {
  if (!dir) return;
  
  TIter next(dir->GetListOfKeys());
  while (TKey *key = (TKey *)next()) {
    TObject *obj = key->ReadObj();
    if (obj->InheritsFrom(TDirectory::Class())) {
      CollectRhoMethods((TDirectory *)obj, outHists, outLabels, outMethods);
      continue;
    }
    
    TString name = obj->GetName();
    
    // 디버깅: 모든 히스토그램 이름 출력
    if (obj->InheritsFrom(TH2::Class())) {
      std::cout << "Found TH2: " << name << std::endl;
    }
    
    if (!NameLooksLikeRhoMethod(name)) continue;

    if (obj->InheritsFrom(TH2::Class())) {
      TH2 *h2 = (TH2 *)obj;
      TH1 *proj = ProjectRhoOverCentrality(h2, Form("%s_MB", name.Data()));
      
      if (proj && proj->Integral() > 0) {
        outHists.push_back(proj);
        outLabels.push_back(name);
        outMethods.push_back(GetMethodName(name));
      }
    }
  }
}

void DrawRhoMethodsComparison(const char *filePath = \
  "/Users/js/cernbox/workspace/O2Physics/jets/data/AnalysisResults/LHC22o_apass7_minBias_small/sel8/jetSpectraCharged/AnalysisResults.root") {
  
  gStyle->SetOptStat(0); // Turn off statistics box
  
  auto file = TFile::Open(filePath, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Failed to open file: " << filePath << std::endl;
    return;
  }

  std::vector<TH1 *> hists;
  std::vector<TString> labels;
  std::vector<TString> methods;
  
  CollectRhoMethods(file, hists, labels, methods);

  if (hists.empty()) {
    std::cerr << "No rho method histograms found in file." << std::endl;
    file->Close();
    return;
  }

  std::cout << "Found " << hists.size() << " rho method histograms:" << std::endl;
  for (size_t i = 0; i < hists.size(); ++i) {
    std::cout << "  " << i+1 << ": " << labels[i] << " (" << methods[i] << ")" << std::endl;
  }

  // Find reference method (first one as reference)
  TH1 *refMethod = nullptr;
  int refIndex = 0; // Use first method as reference
  if (!hists.empty()) {
    refMethod = hists[0];
  }

  // Extract dataset info from file path
  TString datasetInfo = "";
  TString filePathStr(filePath);
  if (filePathStr.Contains("LHC22o")) {
    if (filePathStr.Contains("apass7")) {
      datasetInfo = "LHC22o_pass7";
    } else {
      datasetInfo = "LHC22o";
    }
  } else if (filePathStr.Contains("LHC24f3b")) {
    datasetInfo = "LHC24f3b";
  } else if (filePathStr.Contains("LHC24f3c")) {
    datasetInfo = "LHC24f3c";
  } else if (filePathStr.Contains("LHC25a2b")) {
    datasetInfo = "LHC25a2b";
  } else {
    datasetInfo = "Data";
  }

  // Create canvas
  TCanvas *c = new TCanvas("cRhoMethods", "Rho Methods Comparison", 1200, 800);
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.13);
  optFili(*c, 1, 1, 0, 1, 0);

  TLegend *leg = new TLegend(0.4, 0.65, 0.85, 0.95);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);

  // Find global min/max for Y axis
  double globalMinY = 1e9;
  double globalMaxY = -1e9;
  double minX = -10, maxX = 50; // Fixed X range

  for (size_t i = 0; i < hists.size(); ++i) {
    Color_t col = ColorPalette[i % ColorPalette.size()];
    hoptset(*hists[i], 1, col, 0, 0, 0, 0, 0.8, 20 + (i % 5), 1);
    
    double minY = hists[i]->GetMinimum(0);
    double maxY = hists[i]->GetMaximum();
    if (minY > 0) globalMinY = std::min(globalMinY, minY);
    globalMaxY = std::max(globalMaxY, maxY);
  }

  // Set axis titles and draw
  if (!hists.empty()) {
    hset(*hists[0], "#delta#it{#rho} (GeV/#it{c})", "1/N dN/d#delta#it{#rho}", 1.2, 1.0, 0.055, 0.055, 0.01, 0.005, 0.045, 0.045, 510, 505);
    hists[0]->GetXaxis()->SetRangeUser(minX, maxX);
    // hists[0]->GetYaxis()->SetRangeUser(std::max(1e-6, globalMinY * 0.8), globalMaxY * 1.5);
    hists[0]->GetYaxis()->SetRangeUser(1e-9, globalMaxY * 1.5);
    hists[0]->Draw("PE");

    for (size_t i = 1; i < hists.size(); ++i) {
      hists[i]->Draw("PESAME");
    }


    
    // Add legend entries
    for (size_t i = 0; i < hists.size(); ++i) {
      TString lbl = methods[i];
      if (i == refIndex) {
        lbl += " (Ref)";
      }
      leg->AddEntry(hists[i], lbl, "pe");
    }
    
    // Add dataset info to legend
    leg->AddEntry((TObject*)0, datasetInfo, "");
    leg->Draw();
    c->Update(); // Force canvas update to show in ROOT window
  }

    // Note: No need for second delta rho plot since the histograms already contain delta rho values

  c->SaveAs(Form("/Users/js/cernbox/workspace/O2Physics/macros/DrawJets/plots/AN_Charged-particle-jet-cross-section-in-pp-collisions-at-13.6-TeV/Figures/UEfluctuations/Rho_Methods_Comparison_%s.pdf", datasetInfo.Data()));
  c->SaveAs(Form("/Users/js/cernbox/workspace/O2Physics/macros/DrawJets/plots/AN_Charged-particle-jet-cross-section-in-pp-collisions-at-13.6-TeV/Figures/UEfluctuations/Rho_Methods_Comparison_%s.png", datasetInfo.Data()));

  // file->Close();
  
  std::cout << "Plots saved to plots/Rho_Methods_Comparison.pdf and .png" << std::endl;
}

// Function to draw multiple datasets separately
void DrawMultipleDatasets(const std::vector<const char*> &filePaths) {
  for (size_t i = 0; i < filePaths.size(); ++i) {
    std::cout << "Processing dataset " << i+1 << "/" << filePaths.size() << ": " << filePaths[i] << std::endl;
    DrawRhoMethodsComparison(filePaths[i]);
  }
}

// Example usage function
void DrawExampleDatasets() {
  std::vector<const char*> datasets = {
    "/Users/js/cernbox/workspace/O2Physics/jets/data/AnalysisResults/LHC22o_apass7_minBias_small/sel8/jetSpectraCharged/AnalysisResults.root",

    "/Users/js/cernbox/workspace/O2Physics/jets/mc/AnalysisResults/LHC24f3b/selMC/jetBackgroundAnalysis/AnalysisResults.root",
    "/Users/js/cernbox/workspace/O2Physics/jets/mc/AnalysisResults/LHC24f3c/selMC/jetBackgroundAnalysis/AnalysisResults.root",
    "/Users/js/cernbox/workspace/O2Physics/jets/mc/AnalysisResults/LHC25a2b/jetBackgroundAnalysis/AnalysisResults.root"
  };
  
  DrawMultipleDatasets(datasets);
}