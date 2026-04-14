// Draw tracking efficiency from trackEfficiency task
// Efficiency = (particles with associated tracks) / (all primary particles)
// Uses true particle pT to avoid momentum smearing effects
//
// Note: This assumes acceptSplitCollisions = 0 (NonSplitOnly)

#if defined(__CLING__) || defined(__CINT__) || defined(__ROOTCLING__)

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TMath.h>
#include <TStyle.h>
#include <algorithm>
#include <vector>
#include <iostream>

// ============================================
// CONFIGURATION - Modify these as needed
// ============================================
const char* inputFileName = "../../../jets/AnalysisResults/594013_AnalysisResults.root";
const char* histPathPrefix = "track-efficiency";  // Histogram path prefix in ROOT file

// Eta and phi ranges for projection (use -999 to use full range)
const double etaMin = -0.9;
const double etaMax = 0.9;
const double phiMin = -999;  // -999 means use full range
const double phiMax = -999;  // -999 means use full range

// Plot settings
const double plotPtMin = 0.0;
const double plotPtMax = 100.0;
const bool drawLowPt = true;   // Draw efficiency for 0-10 GeV/c
const bool drawHighPt = true;  // Draw efficiency for 10-100 GeV/c

// Output settings
const char* outputDir = "./plots/TrackingEfficiency";
const char* outputFormat = "pdf";  // "pdf", "png", "root", etc.

// ============================================
// Histogram names
// ============================================
const char* hDenomLowPt = "h3_particle_pt_particle_eta_particle_phi_mcpartofinterest";
const char* hNumLowPt = "h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary";
const char* hDenomHighPt = "h3_particle_pt_high_particle_eta_particle_phi_mcpartofinterest";
const char* hNumHighPt = "h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_primary";

// ============================================
// Helper functions
// ============================================

TH1D* ProjectTH3ToPt(TH3F* h3, const char* name, const char* title, 
                     double etaMinProj = -999, double etaMaxProj = -999,
                     double phiMinProj = -999, double phiMaxProj = -999) {
  if (!h3) {
    std::cout << "Error: Input histogram is null" << std::endl;
    return nullptr;
  }

  // Get axis indices: 0=pt, 1=eta, 2=phi
  TAxis* ptAxis = h3->GetXaxis();
  TAxis* etaAxis = h3->GetYaxis();
  TAxis* phiAxis = h3->GetZaxis();

  // Find bin ranges for projection
  int etaBinMin = 1;
  int etaBinMax = etaAxis->GetNbins();
  int phiBinMin = 1;
  int phiBinMax = phiAxis->GetNbins();

  if (etaMinProj != -999 && etaMaxProj != -999) {
    etaBinMin = etaAxis->FindBin(etaMinProj);
    etaBinMax = etaAxis->FindBin(etaMaxProj);
  }
  if (phiMinProj != -999 && phiMaxProj != -999) {
    phiBinMin = phiAxis->FindBin(phiMinProj);
    phiBinMax = phiAxis->FindBin(phiMaxProj);
  }

  // Project to pT (axis 0) with eta and phi range
  TH1D* h1 = h3->ProjectionX(name, etaBinMin, etaBinMax, phiBinMin, phiBinMax);
  h1->SetTitle(title);
  h1->GetXaxis()->SetTitle(ptAxis->GetTitle());
  h1->GetYaxis()->SetTitle("Counts");

  return h1;
}

TH1D* CalculateEfficiency(TH1D* hNum, TH1D* hDenom, const char* name, const char* title) {
  if (!hNum || !hDenom) {
    std::cout << "Error: Numerator or denominator histogram is null" << std::endl;
    return nullptr;
  }

  TH1D* hEff = (TH1D*)hDenom->Clone(name);
  hEff->SetTitle(title);
  hEff->GetYaxis()->SetTitle("Tracking Efficiency");
  hEff->GetYaxis()->SetRangeUser(0.0, 1.1);
  hEff->Reset();

  int nBins = hEff->GetNbinsX();
  for (int i = 1; i <= nBins; i++) {
    double num = hNum->GetBinContent(i);
    double denom = hDenom->GetBinContent(i);
    
    if (denom > 0) {
      double eff = num / denom;
      double err = TMath::Sqrt(eff * (1 - eff) / denom);
      hEff->SetBinContent(i, eff);
      hEff->SetBinError(i, err);
    } else {
      hEff->SetBinContent(i, 0);
      hEff->SetBinError(i, 0);
    }
  }

  return hEff;
}

void DrawTrackingEfficiency() {
  // ============================================
  // Open input file
  // ============================================
  TFile* file = TFile::Open(inputFileName, "read");
  if (!file || file->IsZombie()) {
    std::cout << "Error: Cannot open file " << inputFileName << std::endl;
    return;
  }
  std::cout << "Opened file: " << inputFileName << std::endl;

  // ============================================
  // Get histograms
  // ============================================
  TString pathDenomLow = Form("%s/%s", histPathPrefix, hDenomLowPt);
  TString pathNumLow = Form("%s/%s", histPathPrefix, hNumLowPt);
  TString pathDenomHigh = Form("%s/%s", histPathPrefix, hDenomHighPt);
  TString pathNumHigh = Form("%s/%s", histPathPrefix, hNumHighPt);

  TH3F* h3DenomLow = (TH3F*)file->Get(pathDenomLow);
  TH3F* h3NumLow = (TH3F*)file->Get(pathNumLow);
  TH3F* h3DenomHigh = (TH3F*)file->Get(pathDenomHigh);
  TH3F* h3NumHigh = (TH3F*)file->Get(pathNumHigh);

  if (!h3DenomLow) std::cout << "Warning: Cannot find " << pathDenomLow << std::endl;
  if (!h3NumLow) std::cout << "Warning: Cannot find " << pathNumLow << std::endl;
  if (!h3DenomHigh) std::cout << "Warning: Cannot find " << pathDenomHigh << std::endl;
  if (!h3NumHigh) std::cout << "Warning: Cannot find " << pathNumHigh << std::endl;

  // ============================================
  // Project to pT
  // ============================================
  TH1D* hDenomLowProj = nullptr;
  TH1D* hNumLowProj = nullptr;
  TH1D* hDenomHighProj = nullptr;
  TH1D* hNumHighProj = nullptr;

  if (drawLowPt && h3DenomLow && h3NumLow) {
    hDenomLowProj = ProjectTH3ToPt(h3DenomLow, "hDenomLowProj", 
                                    "All Primary Particles (0-10 GeV/c);#it{p}_{T} (GeV/#it{c});Counts",
                                    etaMin, etaMax, phiMin, phiMax);
    hNumLowProj = ProjectTH3ToPt(h3NumLow, "hNumLowProj",
                                  "Primary Particles with Associated Tracks (0-10 GeV/c);#it{p}_{T} (GeV/#it{c});Counts",
                                  etaMin, etaMax, phiMin, phiMax);
  }

  if (drawHighPt && h3DenomHigh && h3NumHigh) {
    hDenomHighProj = ProjectTH3ToPt(h3DenomHigh, "hDenomHighProj",
                                     "All Primary Particles (10-100 GeV/c);#it{p}_{T} (GeV/#it{c});Counts",
                                     etaMin, etaMax, phiMin, phiMax);
    hNumHighProj = ProjectTH3ToPt(h3NumHigh, "hNumHighProj",
                                   "Primary Particles with Associated Tracks (10-100 GeV/c);#it{p}_{T} (GeV/#it{c});Counts",
                                   etaMin, etaMax, phiMin, phiMax);
  }

  // ============================================
  // Calculate efficiency
  // ============================================
  TH1D* hEffLow = nullptr;
  TH1D* hEffHigh = nullptr;
  TH1D* hEffMerged = nullptr;
  TH1D* hDenomMerged = nullptr;
  TH1D* hNumMerged = nullptr;

  if (hDenomLowProj && hNumLowProj) {
    hEffLow = CalculateEfficiency(hNumLowProj, hDenomLowProj, "hEffLow",
                                   "Tracking Efficiency (0-10 GeV/c);#it{p}_{T} (GeV/#it{c});Efficiency");
  }

  if (hDenomHighProj && hNumHighProj) {
    hEffHigh = CalculateEfficiency(hNumHighProj, hDenomHighProj, "hEffHigh",
                                    "Tracking Efficiency (10-100 GeV/c);#it{p}_{T} (GeV/#it{c});Efficiency");
  }

  // Merge low and high pT histograms
  if (hDenomLowProj && hNumLowProj && hDenomHighProj && hNumHighProj) {
    // Create merged denominator and numerator
    hDenomMerged = (TH1D*)hDenomLowProj->Clone("hDenomMerged");
    hDenomMerged->SetTitle("All Primary Particles;#it{p}_{T} (GeV/#it{c});Counts");
    hDenomMerged->Add(hDenomHighProj);
    
    hNumMerged = (TH1D*)hNumLowProj->Clone("hNumMerged");
    hNumMerged->SetTitle("Primary Particles with Associated Tracks;#it{p}_{T} (GeV/#it{c});Counts");
    hNumMerged->Add(hNumHighProj);
    
    // Calculate merged efficiency
    hEffMerged = CalculateEfficiency(hNumMerged, hDenomMerged, "hEffMerged",
                                     "Tracking Efficiency;#it{p}_{T} (GeV/#it{c});Efficiency");
  }

  // ============================================
  // Draw plots
  // ============================================
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas* c1 = new TCanvas("c1", "Tracking Efficiency", 1200, 800);
  c1->Divide(2, 2);

  // Panel 1: Merged Numerator and Denominator
  if (hDenomMerged && hNumMerged) {
    c1->cd(1);
    gPad->SetLogy();
    gPad->SetLogx();
    hDenomMerged->SetLineColor(kBlue);
    hDenomMerged->SetLineWidth(2);
    hNumMerged->SetLineColor(kRed);
    hNumMerged->SetLineWidth(2);
    hDenomMerged->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hDenomMerged->Draw("HIST");
    hNumMerged->Draw("HIST SAME");
    
    TLegend* leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg1->AddEntry(hDenomMerged, "All Primary Particles", "l");
    leg1->AddEntry(hNumMerged, "With Associated Tracks", "l");
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->Draw();
  }

  // Panel 2: Merged Efficiency (main plot)
  if (hEffMerged) {
    c1->cd(2);
    gPad->SetLogx();
    gPad->SetGridx();
    gPad->SetGridy();
    hEffMerged->SetLineColor(kBlack);
    hEffMerged->SetMarkerColor(kBlack);
    hEffMerged->SetMarkerStyle(20);
    hEffMerged->SetMarkerSize(0.8);
    hEffMerged->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hEffMerged->Draw("EP");
    
    TLine* line = new TLine(plotPtMin, 1.0, plotPtMax, 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kGray);
    line->Draw();
  }

  // Panel 3: Low pT Efficiency (for reference)
  if (hEffLow) {
    c1->cd(3);
    gPad->SetGridx();
    gPad->SetGridy();
    hEffLow->SetLineColor(kBlack);
    hEffLow->SetMarkerColor(kBlack);
    hEffLow->SetMarkerStyle(20);
    hEffLow->SetMarkerSize(0.8);
    hEffLow->GetXaxis()->SetRangeUser(plotPtMin, 10.0);
    hEffLow->Draw("EP");
    
    TLine* line1 = new TLine(plotPtMin, 1.0, 10.0, 1.0);
    line1->SetLineStyle(2);
    line1->SetLineColor(kGray);
    line1->Draw();
  }

  // Panel 4: High pT Efficiency (for reference)
  if (hEffHigh) {
    c1->cd(4);
    gPad->SetLogx();
    gPad->SetGridx();
    gPad->SetGridy();
    hEffHigh->SetLineColor(kBlack);
    hEffHigh->SetMarkerColor(kBlack);
    hEffHigh->SetMarkerStyle(20);
    hEffHigh->SetMarkerSize(0.8);
    hEffHigh->GetXaxis()->SetRangeUser(10.0, plotPtMax);
    hEffHigh->Draw("EP");
    
    TLine* line2 = new TLine(10.0, 1.0, plotPtMax, 1.0);
    line2->SetLineStyle(2);
    line2->SetLineColor(kGray);
    line2->Draw();
  }

  // Add title
  TLatex* title = new TLatex();
  title->SetNDC();
  title->SetTextSize(0.03);
  title->SetTextAlign(23);
  title->DrawLatex(0.5, 0.995, "Tracking Efficiency (acceptSplitCollisions = 0)");

  // ============================================
  // Save plots
  // ============================================
  TString outputPath = Form("%s/trackingEfficiency.%s", outputDir, outputFormat);
  c1->SaveAs(outputPath);
  std::cout << "Saved plot to: " << outputPath << std::endl;

  // Save merged efficiency plot (main plot)
  if (hEffMerged) {
    TCanvas* cMerged = new TCanvas("cMerged", "Tracking Efficiency Merged", 800, 600);
    gPad->SetLogx();
    gPad->SetGridx();
    gPad->SetGridy();
    hEffMerged->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hEffMerged->Draw("EP");
    TLine* line = new TLine(plotPtMin, 1.0, plotPtMax, 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kGray);
    line->Draw();
    title->DrawLatex(0.5, 0.995, "Tracking Efficiency (acceptSplitCollisions = 0)");
    TString outputPathMerged = Form("%s/trackingEfficiency_merged.%s", outputDir, outputFormat);
    cMerged->SaveAs(outputPathMerged);
    std::cout << "Saved merged efficiency plot to: " << outputPathMerged << std::endl;
  }

  // ============================================
  // Save to ROOT file for browser viewing
  // ============================================
  TString rootOutputPath = Form("%s/trackingEfficiency.root", outputDir);
  TFile* outputFile = new TFile(rootOutputPath, "RECREATE");
  if (outputFile && !outputFile->IsZombie()) {
    outputFile->cd();
    
    // Save canvases
    if (c1) {
      c1->Write("c_trackingEfficiency_4panel");
    }
    if (hEffMerged) {
      TCanvas* cMerged = new TCanvas("c_trackingEfficiency_merged", "Tracking Efficiency Merged", 800, 600);
      gPad->SetLogx();
      gPad->SetGridx();
      gPad->SetGridy();
      hEffMerged->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
      hEffMerged->Draw("EP");
      TLine* line = new TLine(plotPtMin, 1.0, plotPtMax, 1.0);
      line->SetLineStyle(2);
      line->SetLineColor(kGray);
      line->Draw();
      TLatex* title2 = new TLatex();
      title2->SetNDC();
      title2->SetTextSize(0.03);
      title2->SetTextAlign(23);
      title2->DrawLatex(0.5, 0.995, "Tracking Efficiency (acceptSplitCollisions = 0)");
      cMerged->Write();
    }
    
    // Save histograms
    if (hEffMerged) hEffMerged->Write("h_trackingEfficiency_merged");
    if (hEffLow) hEffLow->Write("h_trackingEfficiency_lowPt");
    if (hEffHigh) hEffHigh->Write("h_trackingEfficiency_highPt");
    if (hDenomMerged) hDenomMerged->Write("h_denominator_merged");
    if (hNumMerged) hNumMerged->Write("h_numerator_merged");
    if (hDenomLowProj) hDenomLowProj->Write("h_denominator_lowPt");
    if (hNumLowProj) hNumLowProj->Write("h_numerator_lowPt");
    if (hDenomHighProj) hDenomHighProj->Write("h_denominator_highPt");
    if (hNumHighProj) hNumHighProj->Write("h_numerator_highPt");
    
    outputFile->Close();
    std::cout << "Saved ROOT file to: " << rootOutputPath << std::endl;
    std::cout << "  You can now open this file in ROOT browser to view the histograms" << std::endl;
  } else {
    std::cout << "Warning: Could not create ROOT output file " << rootOutputPath << std::endl;
  }

  // ============================================
  // Print summary
  // ============================================
  std::cout << "\n============================================" << std::endl;
  std::cout << "Tracking Efficiency Summary" << std::endl;
  std::cout << "============================================" << std::endl;
  std::cout << "Input file: " << inputFileName << std::endl;
  std::cout << "Eta range: [" << etaMin << ", " << etaMax << "]" << std::endl;
  if (phiMin != -999 && phiMax != -999) {
    std::cout << "Phi range: [" << phiMin << ", " << phiMax << "]" << std::endl;
  } else {
    std::cout << "Phi range: Full range" << std::endl;
  }
  std::cout << "============================================" << std::endl;

  if (hEffLow) {
    std::cout << "\nLow pT (0-10 GeV/c) Efficiency:" << std::endl;
    int nBins = hEffLow->GetNbinsX();
    for (int i = 1; i <= nBins; i++) {
      double pt = hEffLow->GetBinCenter(i);
      if (pt >= plotPtMin && pt <= 10.0) {
        double eff = hEffLow->GetBinContent(i);
        double err = hEffLow->GetBinError(i);
        if (eff > 0) {
          std::cout << Form("  pT = %.2f GeV/c: %.4f +/- %.4f", pt, eff, err) << std::endl;
        }
      }
    }
  }

  if (hEffHigh) {
    std::cout << "\nHigh pT (10-100 GeV/c) Efficiency:" << std::endl;
    int nBins = hEffHigh->GetNbinsX();
    for (int i = 1; i <= nBins; i++) {
      double pt = hEffHigh->GetBinCenter(i);
      if (pt >= 10.0 && pt <= plotPtMax) {
        double eff = hEffHigh->GetBinContent(i);
        double err = hEffHigh->GetBinError(i);
        if (eff > 0) {
          std::cout << Form("  pT = %.2f GeV/c: %.4f +/- %.4f", pt, eff, err) << std::endl;
        }
      }
    }
  }

  file->Close();
  std::cout << "\nDone!" << std::endl;
}

#endif

