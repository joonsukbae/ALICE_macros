// Track eta and phi QA: Data vs MC comparison
// Uses h2_track_eta_track_phi from jet-spectra-charged
// Produces inclusive eta and phi distributions with ratio pads

#if defined(__CLING__) || defined(__CINT__) || defined(__ROOTCLING__)

#include "drawJetSpectraCharged/Filipad2.h"
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TString.h>
#include <TMath.h>
#include <TLine.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TLatex.h>
#include <iostream>
#include <vector>

// ============================================
// Configuration
// ============================================
const TString mainDir = "~/cernbox/workspace/O2Physics/jets/AnalysisResults/";
const TString fRoot = "_AnalysisResults.root";

// Data file
const TString dataFile = "498133" + fRoot; // 2022
// const TString dataFile = "608785" + fRoot; // 2023
const TString dataDir = "jet-spectra-charged";
const TString dataName = "2022 data (LHC22o pass7)";
// const TString dataName = "2023 data (LHC23 pass4 Thin)";

// MC files
const std::vector<TString> mcFileNames = {
    "604215" + fRoot, // MB MC, w/o pT smearing
    "594275" + fRoot, // MB MC, C tune
};

const std::vector<TString> mcDirs = {
    "jet-spectra-charged",
    "jet-spectra-charged",
};

const std::vector<TString> mcNames = {
    "MB MC, w/o pT smearing",
    "MB MC, #Delta(q/#it{p}_{T}) #times 1.5",
};

// Color palette
const std::vector<Color_t> colors = {
    kBlack,      // Data
    kRed,
    kBlue,
    kGreen + 2,
    kMagenta,
    kCyan,
    kOrange + 7,
};

// Histogram name
const TString histName = "h2_track_eta_track_phi";

// Rebin factors
const int kEtaRebin = 2;
const int kPhiRebin = 4;

// Output directory
const TString outDir = "drawJetSpectraCharged/plots/trackEtaPhiQA";

// ============================================
// Helper: Get number of events
// ============================================
Double_t GetNevents(const char* fileName, const char* dirName) {
    TFile* file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        return 0;
    }

    TH1* hEvents = nullptr;
    TString paths[] = {
        Form("%s/h_collisions_weighted", dirName),
        Form("%s/h_collisions", dirName),
    };

    for (const auto& path : paths) {
        hEvents = (TH1*)file->Get(path);
        if (hEvents && hEvents->Integral() > 0) break;
    }

    if (!hEvents) {
        std::cerr << "Warning: No event histogram found in " << fileName << std::endl;
        file->Close();
        return 1.0;
    }

    Double_t nevents = hEvents->GetBinContent(3);
    if (nevents <= 0) nevents = hEvents->Integral();

    file->Close();
    return nevents;
}

// ============================================
// Helper: Load TH2 and project eta or phi
// ============================================
TH1D* LoadProjection(const char* fileName, const char* dirName, const char* label,
                     Double_t nevents, Color_t color, bool projectEta) {
    TFile* file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        return nullptr;
    }

    TString histPath = Form("%s/%s", dirName, histName.Data());
    TH2* h2 = (TH2*)file->Get(histPath);
    if (!h2) {
        std::cerr << "Error: Cannot find " << histPath << " in " << fileName << std::endl;
        file->Close();
        return nullptr;
    }

    TH1D* h1 = nullptr;
    if (projectEta) {
        // h2 axes: X=eta, Y=phi
        h1 = h2->ProjectionX(Form("eta_%s", label));
    } else {
        // Project phi, restrict to [0, 2pi]
        int binYmin = h2->GetYaxis()->FindBin(0.0 + 1e-6);
        int binYmax = h2->GetYaxis()->FindBin(2 * TMath::Pi() - 1e-6);
        h1 = h2->ProjectionY(Form("phi_%s", label), 0, -1);
        h1->GetXaxis()->SetRangeUser(0, 2 * TMath::Pi());
    }

    if (!h1) {
        file->Close();
        return nullptr;
    }
    h1->SetDirectory(nullptr);

    // Rebin
    int rebin = projectEta ? kEtaRebin : kPhiRebin;
    if (rebin > 1) h1->Rebin(rebin);

    // Normalize by events
    if (nevents > 0) h1->Scale(1.0 / nevents);

    // Normalize by bin width
    for (int i = 1; i <= h1->GetNbinsX(); ++i) {
        double w = h1->GetBinWidth(i);
        if (w > 0) {
            h1->SetBinContent(i, h1->GetBinContent(i) / w);
            h1->SetBinError(i, h1->GetBinError(i) / w);
        }
    }

    // Style
    h1->SetLineColor(color);
    h1->SetMarkerColor(color);
    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(0.8);
    h1->SetLineWidth(2);
    h1->SetTitle("");

    file->Close();
    return h1;
}

// ============================================
// Draw one comparison (eta or phi)
// ============================================
void DrawComparison(bool projectEta) {
    TString varName = projectEta ? "eta" : "phi";
    TString xtitle = projectEta ? "#eta" : "#varphi (rad)";
    TString ytitle = projectEta ? "1/#it{N}_{evt} d#it{N}/d#eta"
                                : "1/#it{N}_{evt} d#it{N}/d#varphi";

    static int padID = 0;
    padID++;

    Filipad2* fpad = new Filipad2(Form("Track%s", varName.Data()), padID, 2, 0.4, 100, 50, 0.7, 1, 1);
    fpad->Draw();

    TPad* mainPad = fpad->GetPad(1);
    TPad* ratioPad = fpad->GetPad(2);

    // Load data
    TString dataPath = mainDir + dataFile;
    Double_t dataNevents = GetNevents(dataPath.Data(), dataDir.Data());
    std::cout << "[Info] Data events: " << dataNevents << std::endl;

    TH1D* hData = LoadProjection(dataPath.Data(), dataDir.Data(), "data",
                                  dataNevents, colors[0], projectEta);
    if (!hData) {
        std::cerr << "Error: Failed to load data!" << std::endl;
        return;
    }

    // Main pad
    mainPad->cd();
    hData->GetXaxis()->SetTitle(xtitle.Data());
    hData->GetYaxis()->SetTitle(ytitle.Data());
    hData->GetYaxis()->SetTitleSize(0.06);
    hData->GetYaxis()->SetTitleOffset(1.0);
    hData->Draw("PE");

    TLegend* leg = new TLegend(0.4, 0.65, 0.9, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(hData, dataName.Data(), "PE");

    // Load and draw MC
    std::vector<TH1D*> ratioHists;
    for (size_t i = 0; i < mcFileNames.size(); ++i) {
        TString mcPath = mainDir + mcFileNames[i];
        Double_t mcNevents = GetNevents(mcPath.Data(), mcDirs[i].Data());
        std::cout << "[Info] " << mcNames[i] << ": " << mcNevents << " events" << std::endl;

        TH1D* hMC = LoadProjection(mcPath.Data(), mcDirs[i].Data(),
                                    Form("mc_%zu", i), mcNevents, colors[i + 1], projectEta);
        if (!hMC) continue;

        mainPad->cd();
        hMC->Draw("PE SAME");
        leg->AddEntry(hMC, mcNames[i].Data(), "PE");

        // Ratio
        TH1D* hRatio = (TH1D*)hMC->Clone(Form("Ratio_%s_%zu", varName.Data(), i));
        hRatio->Divide(hData);
        hRatio->SetLineColor(colors[i + 1]);
        hRatio->SetMarkerColor(colors[i + 1]);
        hRatio->SetMarkerStyle(20);
        hRatio->SetMarkerSize(0.8);
        hRatio->SetLineWidth(2);
        hRatio->GetYaxis()->SetTitle("MC / Data");
        hRatio->GetYaxis()->SetTitleSize(0.10);
        hRatio->GetYaxis()->SetTitleOffset(0.5);
        hRatio->GetYaxis()->SetLabelSize(0.08);
        hRatio->GetYaxis()->SetRangeUser(0.9, 1.1);
        hRatio->GetXaxis()->SetTitle(xtitle.Data());
        hRatio->GetXaxis()->SetTitleSize(0.08);
        hRatio->GetXaxis()->SetTitleOffset(1.0);
        hRatio->GetXaxis()->SetLabelSize(0.08);
        ratioHists.push_back(hRatio);

        ratioPad->cd();
        if (i == 0) hRatio->Draw("PE");
        else        hRatio->Draw("PE SAME");
    }

    mainPad->cd();
    leg->Draw();

    // Reference line
    ratioPad->cd();
    double xmin = hData->GetXaxis()->GetXmin();
    double xmax = hData->GetXaxis()->GetXmax();
    if (!projectEta) { xmin = 0; xmax = 2 * TMath::Pi(); }
    TLine* line = new TLine(xmin, 1.0, xmax, 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kBlack);
    line->Draw();

    mainPad->Update();
    ratioPad->Update();

    // Save PDF with run numbers
    TString runLabel = dataFile(0, dataFile.Index("_"));
    for (const auto& mc : mcFileNames) {
        runLabel += "_" + mc(0, mc.Index("_"));
    }
    TString outName = Form("%s/Track%s_%s.pdf", outDir.Data(), varName.Data(), runLabel.Data());
    fpad->C->SaveAs(outName.Data());
    std::cout << "[Info] Saved to " << outName << std::endl;
}

// ============================================
// Main function
// ============================================
void DrawTrackEtaPhiQA() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    gSystem->mkdir(outDir.Data(), kTRUE);

    std::cout << "========================================" << std::endl;
    std::cout << "DrawTrackEtaPhiQA: Starting..." << std::endl;
    std::cout << "========================================" << std::endl;

    DrawComparison(true);   // eta
    DrawComparison(false);  // phi

    std::cout << "\n========================================" << std::endl;
    std::cout << "DrawTrackEtaPhiQA: Completed!" << std::endl;
    std::cout << "========================================" << std::endl;
}

#endif
