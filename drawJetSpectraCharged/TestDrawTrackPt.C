// Test macro for drawing track pT distributions
// Independent from DrawJetsMCFilesTitles.h and other headers
// Author: Test
// Date: 2026-01-24

#if defined(__CLING__) || defined(__CINT__) || defined(__ROOTCLING__)

#include "Filipad2.h"
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>
#include <vector>

// ============================================
// Configuration
// ============================================
const TString mainDir = "~/cernbox/workspace/O2Physics/jets/AnalysisResults/";
const TString fRoot = "_AnalysisResults.root";

// Data file
// const TString dataFile = "498133" + fRoot; // 2022
const TString dataFile = "608785" + fRoot; // 2023
const TString dataDir = "jet-spectra-charged";
// const TString dataName = "2022 data (LHC22o pass7)";
const TString dataName = "2023 data (LHC23 pass4 Thin)";

// MC files
const std::vector<TString> mcFileNames = {
    // "515446" + fRoot, // MB MC, A tune
    // "555891" + fRoot, // MB MC, B tune
    // "594014" + fRoot, // MB MC, C tune
    // "604215" + fRoot, // MB MC, fix middle only DCA
    
    // "619659" + fRoot, // MB MC, my Tune
    // "619660" + fRoot, // MB MC, only DCA
    // "594275" + fRoot, // MB MC, fix middle C tune

    // "516969" + fRoot, // JJ MC, A tune
    // "596838" + fRoot, // JJ MC, B tune 1.2
    // "596836" + fRoot, // JJ MC, B tune 1.5
    // "596837" + fRoot, // JJ MC, B tune 1.8
    // "593757" + fRoot, // JJ MC, C tune
    // "502419" + fRoot, // JJ MC, D tune (file-based)

    "606088" + fRoot, // 2023 MB MC, A tune
    "594179" + fRoot, // 2023 MB MC, C tune
};

const std::vector<TString> mcDirs = {
    // "jet-spectra-charged_id34413", // MB MC, A tune
    // "jet-spectra-charged",          // MB MC, B tune
    // "jet-spectra-charged",          // MB MC, C tune
    "jet-spectra-charged",          // MB MC, fix middle C tune
    
    // "jet-spectra-charged",          // MB MC, my Tune
    // "jet-spectra-charged",          // MB MC, only DCA
    "jet-spectra-charged",          // MB MC, fix middle only DCA

    // "jet-spectra-charged",          // JJ MC, A tune
    // "jet-spectra-charged",          // JJ MC, B tune 1.2
    // "jet-spectra-charged",          // JJ MC, B tune 1.5
    // "jet-spectra-charged",          // JJ MC, B tune 1.8
    // "jet-spectra-charged",          // JJ MC, C tune
    // "jet-spectra-charged",          // JJ MC, D tune (file-based)
};

const std::vector<TString> mcNames = {
    // "MB MC, A tune",
    // "MB MC, B tune",
    // "MB MC, C tune",
    "2023 MB MC, w/o pT smearing",

    // "MB MC, my Tune",
    // "MB MC, only DCA",
    "2023 MB MC, #Delta(q/#it{p}_{T}) #times 1.5",
    
    // "JJ MC, A tune",
    // "JJ MC, B tune 1.2",
    // "JJ MC, B tune 1.5",
    // "JJ MC, B tune 1.8",
    // "JJ MC, C tune",
    // "JJ MC, D tune (file-based)",
};

// Color palette
const std::vector<Color_t> colors = {
    kBlack,      // Data
    kRed,        // MB MC, A tune
    kBlue,       // MB MC, B tune
    kGreen + 2,  // MB MC, C tune
    kMagenta,    // MB MC, my Tune
    kCyan,       // MB MC, only DCA
    kRed + 1,    // JJ MC, A tune
    kPink - 2,       // JJ MC, B tune 1.2
    kOrange + 7, // JJ MC, B tune 1.5
    kYellow + 2, // JJ MC, B tune 1.8
    kPink + 9,   // JJ MC, C tune
    kViolet - 6, // JJ MC, D tune
    kGray + 2,   // MB MC, my Tune
    kAzure + 2,   // MB MC, only DCA
    kPink - 7,   // JJ MC, A tune
    kSpring + 5, // JJ MC, B tune 1.2
    kGray + 2,   // JJ MC, B tune 1.5
    kAzure + 8,   // JJ MC, B tune 1.8
    kPink + 9,   // JJ MC, C tune
    kViolet - 6, // JJ MC, D tune
};

// Histogram name
const TString histName = "h_track_pt";

// Rebinning settings - edges must align with original h_track_pt binning:
// AxisSpec {200, -0.5, 199.5} => bin edges at -0.5, 0.5, 1.5, 2.5, ... 199.5
const Double_t Trackptbin[] = {0.5, 1.5, 3.5, 5.5, 7.5, 9.5, 14.5, 19.5, 24.5, 29.5, 39.5, 49.5, 59.5, 69.5, 79.5, 89.5, 99.5, 119.5, 139.5, 169.5, 199.5};
const Int_t nTrackptbin = sizeof(Trackptbin) / sizeof(Trackptbin[0]) - 1; // 20
const bool doRebin = true;

// Plot settings
const bool normalizeByEvents = true;
const double plotPtMin = 0.5;
const double plotPtMax = 199.5;

// ============================================
// Helper function: Get number of events
// ============================================
Double_t GetNevents(const char* fileName, const char* dirName) {
    TFile* file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        return 0;
    }
    
    // Try different event histogram names
    TH1* hEvents = nullptr;
    TString paths[] = {
        Form("%s/h_collisions_weighted", dirName),
        Form("%s/h_collisions", dirName),
        Form("%s/h_mcColl_counts_weight", dirName),
        Form("%s/h_mcColl_counts", dirName),
        Form("%s/h_mccollisions", dirName),
    };
    
    for (const auto& path : paths) {
        hEvents = (TH1*)file->Get(path);
        if (hEvents && hEvents->Integral() > 0) {
            break;
        }
    }
    
    if (!hEvents) {
        std::cerr << "Warning: No event histogram found in " << fileName << std::endl;
        file->Close();
        return 1.0; // Default to 1 if not found
    }
    
    // Get bin index 3 (typically the selected events bin)
    Double_t nevents = hEvents->GetBinContent(3);
    if (nevents <= 0) {
        nevents = hEvents->Integral();
    }
    
    file->Close();
    return nevents;
}

// ============================================
// Helper function: Load and normalize histogram
// ============================================
TH1* LoadTrackPt(const char* fileName, const char* dirName, const char* histLabel, 
                 Double_t nevents, Color_t color) {
    TFile* file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        return nullptr;
    }
    
    TString histPath = Form("%s/%s", dirName, histName.Data());
    TH1* h = (TH1*)file->Get(histPath);
    if (!h) {
        std::cerr << "Error: Cannot find histogram " << histPath << " in " << fileName << std::endl;
        file->Close();
        return nullptr;
    }
    
    // Clone and detach from file
    TH1* hClone = (TH1*)h->Clone(Form("TrackPt_%s", histLabel));
    hClone->SetDirectory(nullptr);
    
    // Rebin if requested (same as DrawJetsMCfTrackQA.h)
    if (doRebin) {
        std::cout << "[Info] Rebinning histogram: " << histLabel << std::endl;
        // Create non-const copy of bin edges array for Rebin
        std::vector<Double_t> binEdges(Trackptbin, Trackptbin + sizeof(Trackptbin) / sizeof(Trackptbin[0]));
        TH1* hRebinned = hClone->Rebin(nTrackptbin, Form("TrackPt_rebinned_%s", histLabel), binEdges.data());
        if (hRebinned) {
            hRebinned->SetDirectory(nullptr);
            delete hClone;  // Delete original
            hClone = hRebinned;  // Use rebinned version
        } else {
            std::cerr << "[Warning] Rebin failed for " << histLabel << ", using original" << std::endl;
        }
    }
    
    // Normalize by number of events
    if (normalizeByEvents && nevents > 0) {
        hClone->Scale(1.0 / nevents);
    }

    // Normalize by bin width for proper dN/dpT density
    for (int i = 1; i <= hClone->GetNbinsX(); ++i) {
        double w = hClone->GetBinWidth(i);
        if (w > 0) {
            hClone->SetBinContent(i, hClone->GetBinContent(i) / w);
            hClone->SetBinError(i, hClone->GetBinError(i) / w);
        }
    }
    
    // Set style
    hClone->SetLineColor(color);
    hClone->SetMarkerColor(color);
    hClone->SetMarkerStyle(20);
    hClone->SetLineWidth(2);
    
    // Set axis titles
    hClone->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hClone->GetYaxis()->SetTitle("1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    
    file->Close();
    return hClone;
}

// ============================================
// Main function
// ============================================
void TestDrawTrackPt() {
    std::cout << "========================================" << std::endl;
    std::cout << "TestDrawTrackPt: Starting..." << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Create canvas with Filipad2
    Filipad2* pad = new Filipad2("TestTrackPt", 1, 2, 0.4, 100, 50, 0.7, 1, 1);
    pad->Draw();
    
    TPad* mainPad = pad->GetPad(1);
    TPad* ratioPad = pad->GetPad(2);
    
    // Set up main pad
    mainPad->cd();
    mainPad->SetLogy();
    // x-axis is linear
    
    // Set up ratio pad
    ratioPad->cd();
    // x-axis is linear
    
    // Load data histogram
    std::cout << "\n[Info] Loading data histogram..." << std::endl;
    TString dataPath = mainDir + dataFile;
    Double_t dataNevents = GetNevents(dataPath.Data(), dataDir.Data());
    std::cout << "[Info] Data events: " << dataNevents << std::endl;
    
    TH1* hData = LoadTrackPt(dataPath.Data(), dataDir.Data(), "data", dataNevents, colors[0]);
    if (!hData) {
        std::cerr << "Error: Failed to load data histogram!" << std::endl;
        return;
    }
    
    // Set data style
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(1.2);
    hData->SetLineWidth(2);
    
    // Draw data on main pad
    mainPad->cd();
    hData->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hData->GetYaxis()->SetTitleSize(0.06);
    hData->GetYaxis()->SetTitleOffset(1.0);
    hData->Draw("PE");
    
    // Create legend
    TLegend* leg = new TLegend(0.4, 0.6, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(hData, dataName.Data(), "PE");
    
    // Load and draw MC histograms
    std::cout << "\n[Info] Loading MC histograms..." << std::endl;
    std::vector<TH1*> mcHists;
    std::vector<TH1*> ratioHists;
    
    for (size_t i = 0; i < mcFileNames.size(); ++i) {
        TString mcPath = mainDir + mcFileNames[i];
        std::cout << "[Info] Loading: " << mcNames[i] << " from " << mcFileNames[i] << std::endl;
        
        Double_t mcNevents = GetNevents(mcPath.Data(), mcDirs[i].Data());
        std::cout << "[Info]   Events: " << mcNevents << std::endl;
        
        TH1* hMC = LoadTrackPt(mcPath.Data(), mcDirs[i].Data(), 
                              Form("mc_%zu", i), mcNevents, colors[i + 1]);
        if (!hMC) {
            std::cerr << "[Warning] Failed to load " << mcNames[i] << std::endl;
            continue;
        }
        
        mcHists.push_back(hMC);
        
        // Draw on main pad
        mainPad->cd();
        hMC->Draw("PE SAME");
        leg->AddEntry(hMC, mcNames[i].Data(), "PE");
        
        // Create ratio histogram
        TH1* hRatio = (TH1*)hMC->Clone(Form("Ratio_%zu", i));
        hRatio->Divide(hData);
        hRatio->SetLineColor(colors[i + 1]);
        hRatio->SetMarkerColor(colors[i + 1]);
        hRatio->SetMarkerStyle(20);
        hRatio->SetLineWidth(2);
        hRatio->GetYaxis()->SetTitle("MC / Data");
        hRatio->GetYaxis()->SetTitleSize(0.10);
        hRatio->GetYaxis()->SetTitleOffset(0.5);
        hRatio->GetYaxis()->SetLabelSize(0.08);
        hRatio->GetYaxis()->SetRangeUser(0.0, 2.0);
        hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hRatio->GetXaxis()->SetTitleSize(0.08);
        hRatio->GetXaxis()->SetTitleOffset(1.0);
        hRatio->GetXaxis()->SetLabelSize(0.08);
        ratioHists.push_back(hRatio);
        
        // Draw on ratio pad
        ratioPad->cd();
        if (i == 0) {
            hRatio->Draw("PE");
        } else {
            hRatio->Draw("PE SAME");
        }
    }
    
    // Draw legend
    mainPad->cd();
    leg->Draw();
    
    // Draw reference line at y=1 on ratio pad
    ratioPad->cd();
    TLine* line = new TLine(plotPtMin, 1.0, plotPtMax, 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kBlack);
    line->Draw();
    
    // Update pads
    mainPad->Update();
    ratioPad->Update();

    // Save to PDF with run numbers in filename
    TString outDir = "plots/trackPtQA";
    gSystem->mkdir(outDir.Data(), kTRUE);
    TString runLabel = dataFile(0, dataFile.Index("_"));
    for (const auto& mc : mcFileNames) {
        runLabel += "_" + mc(0, mc.Index("_"));
    }
    TString outName = Form("%s/TrackPt_%s.pdf", outDir.Data(), runLabel.Data());
    pad->C->SaveAs(outName.Data());
    std::cout << "[Info] Saved to " << outName << std::endl;

    std::cout << "\n========================================" << std::endl;
    std::cout << "TestDrawTrackPt: Completed!" << std::endl;
    std::cout << "========================================" << std::endl;
}

#endif
