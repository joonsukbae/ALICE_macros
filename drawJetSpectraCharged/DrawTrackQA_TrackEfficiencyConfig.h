///////////////////////////////////////////////////
// DrawTrackQA_TrackEfficiencyConfig.h
// Configuration, constants, and shared utilities
// for the Track QA (track-efficiency) macro
///////////////////////////////////////////////////

#ifndef DRAWTRACKQACONFIG_H
#define DRAWTRACKQACONFIG_H

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
#include <TGraphErrors.h>
#include <TF1.h>
#include <iostream>
#include <vector>
#include <utility>

// ============================================================
// File Configuration (vector-based, DrawJetsMC pattern)
// ============================================================

// Base directory for all AnalysisResults files
TString trkqa_mainDir = "~/cernbox/workspace/O2Physics/jets/AnalysisResults/";

// Data reference
TString trkqa_refPath   = trkqa_mainDir + "496213_AnalysisResults.root";
TString trkqa_refDir    = "track-efficiency";
TString trkqa_DataLabel = "Data (2023 pass4)";

// MC files (vector)
std::vector<TString> trkqa_fileNames = {
    "593757_AnalysisResults.root",   // LHC25a2b, tuner C
    "628390_AnalysisResults.root",   // LHC23k4i, custom tuner
    // add more as needed
};
std::vector<TString> trkqa_histNames = {
    "2023 MB MC",
    "2023 JJ MC (pTHat 8.0x4.0)",
};
std::vector<TString> trkqa_McFileDirs = {
    "track-efficiency",
    "track-efficiency",
};

// Colors: [0]=Data, [1..]=MC files
std::vector<Color_t> trkqa_Colors = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+1};

// Marker styles: [0]=Data, [1..]=MC files
const int kDataMarker = 20;
const int kMCMarkers[] = {24, 25, 26, 27, 28, 30};
const double kMarkerSize = 0.6;

// ============================================================
// Process Switches (set to 0 to skip)
// ============================================================
const int TrackSelectionProcess  = 1;  // TPC rows, chi2, DCA 1D/2D/profile
const int EfficiencyProcess      = 1;  // efficiency, fake rate, secondary contamination
const int DCAProcess             = 1;  // DCA pT slices, width vs pT, tail fractions

// ============================================================
// Configuration Constants
// ============================================================
TString trkqa_outputDir = "plots/TrackQA_TrackEfficiency";

// pT range
const double kPtMin = 0.0;
const double kPtMax = 200.0;

// DCA pT slices
const int kNPtSlices = 8;
const double kPtSlicesLo[] = {0.15, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0};
const double kPtSlicesHi[] = {0.5,  1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0};

// DCA z cut values for tail fraction study
const int kNDCAzCuts = 3;
const double kDCAzCuts[] = {0.1, 0.2, 0.5};
const Color_t kDCAzCutColors[] = {(Color_t)(kRed+1), (Color_t)(kBlue+1), (Color_t)(kGreen+2)};
const int kDCAzCutStyles[] = {2, 1, 7};  // dashed, solid, long-dashed

// pT slice colors
const Color_t kSliceColors[] = {kBlack, (Color_t)(kRed+1), (Color_t)(kBlue+1), (Color_t)(kGreen+2),
                                (Color_t)(kMagenta+1), (Color_t)(kOrange+1), (Color_t)(kCyan+2), (Color_t)(kViolet+1)};

// ============================================================
// Histogram Names
// ============================================================

// Track selection 1D
const char* kHistTPCRows     = "h_trackselplot_tpccrossedrows";
const char* kHistTPCRowsFrac = "h_trackselplot_tpccrossedrowsoverfindable";
const char* kHistChi2TPC     = "h_trackselplot_chi2ncls_tpc";
const char* kHistChi2ITS     = "h_trackselplot_chi2ncls_its";
const char* kHistDCAxy       = "h_trackselplot_dcaxy";
const char* kHistDCAz        = "h_trackselplot_dcaz";

// Track selection 2D (pT vs property)
const char* kHist2DTPCRows     = "h2_trackselplot_pt_tpccrossedrows";
const char* kHist2DTPCRowsFrac = "h2_trackselplot_pt_tpccrossedrowsoverfindable";
const char* kHist2DChi2TPC     = "h2_trackselplot_pt_chi2ncls_tpc";
const char* kHist2DChi2ITS     = "h2_trackselplot_pt_chi2ncls_its";
const char* kHist2DDCAxy       = "h2_trackselplot_pt_dcaxy";
const char* kHist2DDCAz        = "h2_trackselplot_pt_dcaz";

// Efficiency 3D histograms
const char* kHistTruthLo = "h3_particle_pt_particle_eta_particle_phi_mcpartofinterest";
const char* kHistMatchLo = "h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary";
const char* kHistTruthHi = "h3_particle_pt_high_particle_eta_particle_phi_mcpartofinterest";
const char* kHistMatchHi = "h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_primary";
const char* kHistFake    = "h3_track_pt_track_eta_track_phi_nonassociatedtrack";
const char* kHistRecoPrim = "h3_track_pt_track_eta_track_phi_associatedtrack_primary";
const char* kHistRecoSec  = "h3_track_pt_track_eta_track_phi_associatedtrack_nonprimary";

// ============================================================
// Track Selection Variable Struct
// ============================================================
struct TrackQAVar {
    TString hist1D;
    TString hist2D;
    TString xTitle;
    TString shortName;
    double xMin, xMax;
    double ratioMin, ratioMax;
    bool logY1D;
};

std::vector<TrackQAVar> GetTrackQAVariables() {
    return {
        {kHistTPCRows, kHist2DTPCRows,
         "TPC crossed rows", "TPCCrossedRows",
         -0.5, 164.5, 0.5, 1.5, false},
        {kHistTPCRowsFrac, kHist2DTPCRowsFrac,
         "TPC crossed rows / findable", "TPCCrossedRowsOverFindable",
         0.0, 1.2, 0.5, 1.5, true},
        {kHistChi2TPC, kHist2DChi2TPC,
         "TPC #chi^{2}/cluster", "Chi2TPC",
         0.0, 10.0, 0.5, 1.5, false},
        {kHistChi2ITS, kHist2DChi2ITS,
         "ITS #chi^{2}/cluster", "Chi2ITS",
         0.0, 40.0, 0.5, 1.5, true},
        {kHistDCAxy, kHist2DDCAxy,
         "DCA_{xy} (cm)", "DCAxy",
         -1.0, 1.0, 0.0, 2.5, true},
        {kHistDCAz, kHist2DDCAz,
         "DCA_{z} (cm)", "DCAz",
         -4.0, 4.0, 0.0, 2.5, true},
    };
}

// ============================================================
// Shared Helper Functions
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
template <typename T>
T* GetHist(TDirectory* dir, const char* name) {
    T* h = dynamic_cast<T*>(dir->Get(name));
    if (!h) {
        std::cerr << "[Warning] Histogram not found: " << name << std::endl;
        return nullptr;
    }
    h->SetDirectory(0);
    return h;
}

// Open ROOT file safely
TFile* trkqa_OpenFile(const TString& path) {
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[Warning] Cannot open file: " << path << std::endl;
        if (f) delete f;
        return nullptr;
    }
    return f;
}

// Get directory from file
TDirectory* trkqa_GetDir(TFile* f, const TString& dirName) {
    TDirectory* dir = (TDirectory*)f->Get(dirName);
    if (!dir) {
        std::cerr << "[Warning] Directory not found: " << dirName << " in " << f->GetName() << std::endl;
    }
    return dir;
}

// Gaussian fit to core: iterative 2-sigma fit
// Returns {sigma, sigmaErr}
std::pair<double,double> FitCoreSigma(TH1* h, double fitRangeInit = 0) {
    double rms = h->GetStdDev();
    double mean = h->GetMean();
    double range = (fitRangeInit > 0) ? fitRangeInit : 2.0 * rms;

    TF1 gaus("gaus_tmp", "gaus", mean - range, mean + range);
    gaus.SetParameters(h->GetMaximum(), mean, rms * 0.5);
    h->Fit(&gaus, "QNR0");

    // Refit with refined range
    double sig1 = gaus.GetParameter(2);
    double mu1  = gaus.GetParameter(1);
    gaus.SetRange(mu1 - 2.0 * sig1, mu1 + 2.0 * sig1);
    h->Fit(&gaus, "QNR0");

    return {std::abs(gaus.GetParameter(2)), gaus.GetParError(2)};
}

// Build combined efficiency from low pT + high pT 3D histograms
// Returns merged efficiency histogram [0.15, 100] GeV with variable binning
TH1D* BuildCombinedEfficiency(TDirectory* dir, const char* truthLoName, const char* matchLoName,
                               const char* truthHiName, const char* matchHiName,
                               const char* suffix) {
    TH3* h3TruthLo = GetHist<TH3>(dir, truthLoName);
    TH3* h3MatchLo = GetHist<TH3>(dir, matchLoName);
    if (!h3TruthLo || !h3MatchLo) {
        if (h3TruthLo) delete h3TruthLo;
        if (h3MatchLo) delete h3MatchLo;
        return nullptr;
    }

    // |eta| < 0.9 cut
    int etaLo = h3TruthLo->GetYaxis()->FindBin(-0.899);
    int etaHi = h3TruthLo->GetYaxis()->FindBin(0.899);

    TH1* hTruthLo = h3TruthLo->ProjectionX(Form("hTruthLo_%s", suffix), etaLo, etaHi, 0, -1);
    TH1* hMatchLo = h3MatchLo->ProjectionX(Form("hMatchLo_%s", suffix), etaLo, etaHi, 0, -1);
    hTruthLo->SetDirectory(0);
    hMatchLo->SetDirectory(0);

    TH1* hEffLo = (TH1*)hMatchLo->Clone(Form("hEffLo_%s", suffix));
    hEffLo->SetDirectory(0);
    hEffLo->Divide(hMatchLo, hTruthLo, 1., 1., "B");

    // Try high pT
    TH3* h3TruthHi = GetHist<TH3>(dir, truthHiName);
    TH3* h3MatchHi = GetHist<TH3>(dir, matchHiName);
    TH1* hEffHi = nullptr;
    if (h3TruthHi && h3MatchHi) {
        int etaLoHi = h3TruthHi->GetYaxis()->FindBin(-0.899);
        int etaHiHi = h3TruthHi->GetYaxis()->FindBin(0.899);
        TH1* hTruthHi = h3TruthHi->ProjectionX(Form("hTruthHi_%s", suffix), etaLoHi, etaHiHi, 0, -1);
        TH1* hMatchHi = h3MatchHi->ProjectionX(Form("hMatchHi_%s", suffix), etaLoHi, etaHiHi, 0, -1);
        hTruthHi->SetDirectory(0);
        hMatchHi->SetDirectory(0);
        hEffHi = (TH1*)hMatchHi->Clone(Form("hEffHi_%s", suffix));
        hEffHi->SetDirectory(0);
        hEffHi->Divide(hMatchHi, hTruthHi, 1., 1., "B");
        delete hTruthHi; delete hMatchHi;
    }

    // Merge into variable-binned histogram
    const int nBinsLo = 200, nBinsHi = 18;
    const int nBinsCombined = nBinsLo + nBinsHi;
    double xbins[nBinsCombined + 1];
    for (int i = 0; i <= nBinsLo; i++) xbins[i] = 0.0 + i * 0.05;
    for (int i = 1; i <= nBinsHi; i++) xbins[nBinsLo + i] = 10.0 + i * 5.0;

    TH1D* hEffCombined = new TH1D(Form("hEffCombined_%s", suffix), "", nBinsCombined, xbins);
    hEffCombined->SetDirectory(0);

    for (int ib = 1; ib <= hEffLo->GetNbinsX(); ib++) {
        double pt = hEffLo->GetBinCenter(ib);
        int binC = hEffCombined->FindBin(pt);
        hEffCombined->SetBinContent(binC, hEffLo->GetBinContent(ib));
        hEffCombined->SetBinError(binC, hEffLo->GetBinError(ib));
    }
    if (hEffHi) {
        for (int ib = 1; ib <= hEffHi->GetNbinsX(); ib++) {
            double pt = hEffHi->GetBinCenter(ib);
            int binC = hEffCombined->FindBin(pt);
            hEffCombined->SetBinContent(binC, hEffHi->GetBinContent(ib));
            hEffCombined->SetBinError(binC, hEffHi->GetBinError(ib));
        }
    }

    // Cleanup
    delete h3TruthLo; delete h3MatchLo;
    delete hTruthLo; delete hMatchLo; delete hEffLo;
    if (h3TruthHi) delete h3TruthHi;
    if (h3MatchHi) delete h3MatchHi;
    if (hEffHi) delete hEffHi;

    return hEffCombined;
}

// Build fraction histogram from 3D: numerator / (numerator + denominator)
// Used for fake rate and secondary fraction
TH1* BuildFractionFromH3(TDirectory* dir, const char* numName, const char* denomName,
                          const char* suffix) {
    TH3* h3Num = GetHist<TH3>(dir, numName);
    TH3* h3Den = GetHist<TH3>(dir, denomName);
    if (!h3Num || !h3Den) {
        if (h3Num) delete h3Num;
        if (h3Den) delete h3Den;
        return nullptr;
    }

    int etaLo = h3Num->GetYaxis()->FindBin(-0.899);
    int etaHi = h3Num->GetYaxis()->FindBin(0.899);
    TH1* hNum = h3Num->ProjectionX(Form("hFracNum_%s", suffix), etaLo, etaHi, 0, -1);
    TH1* hDen = h3Den->ProjectionX(Form("hFracDen_%s", suffix), etaLo, etaHi, 0, -1);
    hNum->SetDirectory(0);
    hDen->SetDirectory(0);

    TH1* hDenom = (TH1*)hNum->Clone(Form("hFracDenom_%s", suffix));
    hDenom->Add(hDen);
    hDenom->SetDirectory(0);

    TH1* hFrac = (TH1*)hNum->Clone(Form("hFrac_%s", suffix));
    hFrac->SetDirectory(0);
    hFrac->Divide(hNum, hDenom, 1., 1., "B");

    delete h3Num; delete h3Den;
    delete hNum; delete hDen; delete hDenom;
    return hFrac;
}

#endif
