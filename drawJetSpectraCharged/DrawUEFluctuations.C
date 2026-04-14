///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// UE Fluctuation QA and Systematic //////////
////////// For MB-only pp analysis          //////////
////////// author: Joonsuk Bae              //////////
////////// Last Modified: 2025              //////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

#include "DrawJetsMCFilesTitles.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TLine.h>
#include <vector>
#include <string>
#include <iostream>
#include <cmath>

// ============================================================
// Configuration: R-dependent file list
// ============================================================

struct UEConfig {
  double R;
  TString dataFile;
  TString mcFile;
  int color;
  int markerStyle;
};

const int kNDatasets = 3;
const char* kPlotSaveNames[] = {
  "LHC22o-pass7-small_LHC24f3c",
  "LHC22o-pass7-small_LHC25a2b_tunerA",
  "LHC23-pass4-Thin_small_LHC23k4h"
};
const char* kDatasetLabels[] = {
  "2022 small + MB MC (LHC24f3c) tuner C",
  "2022 full + JJ MC (LHC25a2b) tuner A",
  "2023 small + MB MC (LHC23k4h)"
};

const char* kDirName = "jet-background-analysis-task";

std::vector<UEConfig> GetUEConfigs(int dataset) {
  std::vector<UEConfig> configs;
  const TString basePath = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults";

  const int colors[]  = {kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+1, kCyan+1, kYellow+1};
  const int markers[] = {20, 21, 22, 23, 24, 25, 26};
  const double rVals[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};

  // {dataID, mcID} for each R
  const char* ids[3][7][2] = {
    // 0: 2022 small data + MB MC (LHC24f3c) tuner C
    {{"592451","598076"}, {"592377","594013"}, {"592379","598075"},
     {"592378","594014"}, {"592380","598074"}, {"592381","598073"}, {"592382","594015"}},
    // 1: 2022 full data + JJ MC (LHC25a2b) tuner A
    {{"592451",""}, {"592377","599366"}, {"592379","599367"},
     {"592378","605663"}, {"592380","599368"}, {"592381","599369"}, {"592382","599370"}},
    // 2: 2023 small data + MB MC (LHC23k4h)
    {{"608780","606091"}, {"","606086"}, {"608779","606087"},
     {"608785","606088"}, {"608782","606089"}, {"608783","606090"}, {"608784","606092"}}
  };

  for (int i = 0; i < 7; ++i) {
    TString dataID = ids[dataset][i][0];
    TString mcID   = ids[dataset][i][1];
    if (dataID.Length() == 0 && mcID.Length() == 0) continue;
    configs.push_back(UEConfig{
        rVals[i],
        dataID.Length() > 0 ? Form("%s/%s_AnalysisResults.root", basePath.Data(), dataID.Data()) : "",
        mcID.Length() > 0   ? Form("%s/%s_AnalysisResults.root", basePath.Data(), mcID.Data())   : "",
        colors[i], markers[i]});
  }

  return configs;
}

// ============================================================
// Helper functions
// ============================================================

// Project 2D histogram to 1D by integrating over X axis (centrality)
TH1* ProjectY(TH2* h2, const char* name) {
  if (!h2) return nullptr;
  TH1* h1 = h2->ProjectionY(name);
  h1->SetDirectory(0);
  return h1;
}

// Get mean rho from h2_centrality_rho (integrate over centrality)
Double_t GetMeanRho(TH2* h2) {
  if (!h2) return -1.0;
  TH1* hProj = h2->ProjectionY("_temp_rho");
  Double_t mean = hProj->GetMean();
  delete hProj;
  return mean;
}

// Get RMS of delta pT distribution
Double_t GetDeltaPtRMS(TH2* h2, Double_t& rmsErr) {
  if (!h2) return -1.0;
  TH1* hProj = h2->ProjectionY("_temp_dpt");
  Double_t rms = hProj->GetRMS();
  rmsErr = hProj->GetRMSError();
  delete hProj;
  return rms;
}

// ============================================================
// Main plotting functions
// ============================================================

void DrawRhoDistribution(int dataset) {
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.05);

  auto configs = GetUEConfigs(dataset);
  const char* saveName = kPlotSaveNames[dataset];

  // Rho is a global event property (median pT density) — R-independent.
  // Show Data vs MC comparison using R=0.4 as representative.
  TCanvas *c = new TCanvas(Form("cRhoDist_%d", dataset), "UE #rho distribution: Data vs MC", 800, 600);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.05);

  TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);

  // Use R=0.4 config as representative (rho is R-independent)
  UEConfig *repCfg = nullptr;
  for (auto &cfg : configs) {
    if (fabs(cfg.R - 0.4) < 0.01) { repCfg = &cfg; break; }
  }
  if (!repCfg) return;

  // Data
  if (repCfg->dataFile.Length() > 0) {
    TFile *fData = TFile::Open(repCfg->dataFile, "READ");
    if (fData && !fData->IsZombie()) {
      TH2 *h2cent_rho = (TH2*)fData->Get(Form("%s/h2_centrality_rho", kDirName));
      if (h2cent_rho) {
        TH1 *hRho = ProjectY(h2cent_rho, Form("hRho_Data_%d", dataset));
        hRho->SetLineColor(kBlack);
        hRho->SetLineWidth(1);
        if (hRho->Integral() > 0) hRho->Scale(1.0 / hRho->Integral(), "width");
        hRho->GetXaxis()->SetRangeUser(0, 20);
        hRho->GetXaxis()->SetTitle("#rho_{UE} (GeV/area)");
        hRho->GetYaxis()->SetTitle("(1/N) dN/d#rho");
        hRho->GetYaxis()->SetTitleOffset(1.3);
        hRho->SetTitle("#rho distribution (R-independent)");
        hRho->Draw("HIST");
        leg->AddEntry(hRho, "Data", "l");
      }
      fData->Close();
    }
  }

  // MC
  if (repCfg->mcFile.Length() > 0) {
    TFile *fMC = TFile::Open(repCfg->mcFile, "READ");
    if (fMC && !fMC->IsZombie()) {
      TH2 *h2cent_rho = (TH2*)fMC->Get(Form("%s/h2_centrality_rho", kDirName));
      if (h2cent_rho) {
        TH1 *hRho = ProjectY(h2cent_rho, Form("hRho_MC_%d", dataset));
        hRho->SetLineColor(kRed+1);
        hRho->SetLineWidth(1);
        if (hRho->Integral() > 0) hRho->Scale(1.0 / hRho->Integral(), "width");
        hRho->Draw("HIST SAME");
        leg->AddEntry(hRho, "MB MC", "l");
      }
      fMC->Close();
    }
  }

  leg->Draw();
  c->SaveAs(Form("plots/UEFluctuations/%s/UE_RhoDistribution.pdf", saveName));
}

void DrawDeltaPtDistribution(int dataset) {
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.05);

  auto configs = GetUEConfigs(dataset);
  const char* saveName = kPlotSaveNames[dataset];

  TCanvas *c = new TCanvas(Form("cDeltaPt_%d", dataset), "#delta p_{T} distribution (Random Cone)", 1200, 600);
  c->Divide(2, 1);

  // Left pad: Data
  c->cd(1);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.05);

  TLegend *legData = new TLegend(0.55, 0.60, 0.88, 0.88);
  legData->SetBorderSize(0);
  legData->SetFillStyle(0);
  legData->SetTextSize(0.04);

  bool firstData = true;

  for (const auto &cfg : configs) {
    if (cfg.dataFile.Length() == 0) continue;
    TFile *f = TFile::Open(cfg.dataFile, "READ");
    if (!f || f->IsZombie()) continue;

    TH2 *h2cent_dpt = (TH2*)f->Get(Form("%s/h2_centrality_rhorandomcone", kDirName));
    if (!h2cent_dpt) { f->Close(); continue; }

    TH1 *hDpt = ProjectY(h2cent_dpt, Form("hDeltaPt_R%.1f_Data_%d", cfg.R, dataset));
    hDpt->SetLineColor(cfg.color);
    hDpt->SetLineWidth(1);

    if (hDpt->Integral() > 0) {
      hDpt->Scale(1.0 / hDpt->Integral(), "width");
    }

    if (firstData) {
      hDpt->GetXaxis()->SetRangeUser(-10, 50);
      hDpt->GetXaxis()->SetTitle("#delta p_{T} = p_{T,RC} - A_{RC} #times #rho (GeV/c)");
      hDpt->GetYaxis()->SetTitle("(1/N) dN/d#delta p_{T}");
      hDpt->GetYaxis()->SetTitleOffset(1.3);
      hDpt->SetTitle("Data (Random Cone)");
      hDpt->Draw("HIST");
      firstData = false;
    } else {
      hDpt->Draw("HIST SAME");
    }

    legData->AddEntry(hDpt, Form("R = %.1f", cfg.R), "l");
    f->Close();
  }
  legData->Draw();

  // Right pad: MC
  c->cd(2);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.05);

  TLegend *legMC = new TLegend(0.55, 0.60, 0.88, 0.88);
  legMC->SetBorderSize(0);
  legMC->SetFillStyle(0);
  legMC->SetTextSize(0.04);

  bool firstMC = true;

  for (const auto &cfg : configs) {
    if (cfg.mcFile.Length() == 0) continue;
    TFile *f = TFile::Open(cfg.mcFile, "READ");
    if (!f || f->IsZombie()) continue;

    TH2 *h2cent_dpt = (TH2*)f->Get(Form("%s/h2_centrality_rhorandomcone", kDirName));
    if (!h2cent_dpt) { f->Close(); continue; }

    TH1 *hDpt = ProjectY(h2cent_dpt, Form("hDeltaPt_R%.1f_MC_%d", cfg.R, dataset));
    hDpt->SetLineColor(cfg.color);
    hDpt->SetLineWidth(1);

    if (hDpt->Integral() > 0) {
      hDpt->Scale(1.0 / hDpt->Integral(), "width");
    }

    if (firstMC) {
      hDpt->GetXaxis()->SetRangeUser(-10, 50);
      hDpt->GetXaxis()->SetTitle("#delta p_{T} = p_{T,RC} - A_{RC} #times #rho (GeV/c)");
      hDpt->GetYaxis()->SetTitle("(1/N) dN/d#delta p_{T}");
      hDpt->GetYaxis()->SetTitleOffset(1.3);
      hDpt->SetTitle("MB MC (Random Cone)");
      hDpt->Draw("HIST");
      firstMC = false;
    } else {
      hDpt->Draw("HIST SAME");
    }

    legMC->AddEntry(hDpt, Form("R = %.1f", cfg.R), "l");
    f->Close();
  }
  legMC->Draw();

  c->SaveAs(Form("plots/UEFluctuations/%s/UE_DeltaPtDistribution.pdf", saveName));
}

void DrawDeltaPtRMSvsR(int dataset) {
  gStyle->SetOptStat(0);

  auto configs = GetUEConfigs(dataset);
  const char* saveName = kPlotSaveNames[dataset];

  TCanvas *c = new TCanvas(Form("cDeltaPtRMS_%d", dataset), "#sigma_{#delta p_{T}} vs R (for systematic)", 800, 600);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.05);

  std::vector<Double_t> rValsData, sigmaData, sigmaErrData;
  std::vector<Double_t> rValsMC, sigmaMC, sigmaErrMC;

  for (const auto &cfg : configs) {
    // Data
    if (cfg.dataFile.Length() > 0) {
      TFile *fData = TFile::Open(cfg.dataFile, "READ");
      if (fData && !fData->IsZombie()) {
        TH2 *h2cent_dpt = (TH2*)fData->Get(Form("%s/h2_centrality_rhorandomcone", kDirName));
        if (h2cent_dpt) {
          Double_t rmsErr = 0.0;
          Double_t rms = GetDeltaPtRMS(h2cent_dpt, rmsErr);
          rValsData.push_back(cfg.R);
          sigmaData.push_back(rms);
          sigmaErrData.push_back(rmsErr);
        }
        fData->Close();
      }
    }

    // MC
    if (cfg.mcFile.Length() > 0) {
      TFile *fMC = TFile::Open(cfg.mcFile, "READ");
      if (fMC && !fMC->IsZombie()) {
        TH2 *h2cent_dpt = (TH2*)fMC->Get(Form("%s/h2_centrality_rhorandomcone", kDirName));
        if (h2cent_dpt) {
          Double_t rmsErr = 0.0;
          Double_t rms = GetDeltaPtRMS(h2cent_dpt, rmsErr);
          rValsMC.push_back(cfg.R);
          sigmaMC.push_back(rms);
          sigmaErrMC.push_back(rmsErr);
        }
        fMC->Close();
      }
    }
  }

  TGraphErrors *gData = new TGraphErrors(rValsData.size());
  TGraphErrors *gMC = new TGraphErrors(rValsMC.size());

  for (size_t i = 0; i < rValsData.size(); ++i) {
    gData->SetPoint(i, rValsData[i], sigmaData[i]);
    gData->SetPointError(i, 0.0, sigmaErrData[i]);
  }
  for (size_t i = 0; i < rValsMC.size(); ++i) {
    gMC->SetPoint(i, rValsMC[i], sigmaMC[i]);
    gMC->SetPointError(i, 0.0, sigmaErrMC[i]);
  }

  gData->SetMarkerStyle(20);
  gData->SetMarkerColor(kBlack);
  gData->SetLineColor(kBlack);
  gData->SetLineWidth(1);
  gData->SetMarkerSize(1.3);
  gData->GetYaxis()->SetRangeUser(0.0, 1.5);

  gMC->SetMarkerStyle(24);
  gMC->SetMarkerColor(kRed);
  gMC->SetLineColor(kRed);
  gMC->SetLineWidth(1);
  gMC->SetMarkerSize(1.3);

  gData->GetXaxis()->SetTitle("R");
  gData->GetYaxis()->SetTitle("#sigma_{#delta p_{T}} (GeV/c)");
  gData->GetYaxis()->SetTitleOffset(1.3);
  gData->SetTitle("");
  gData->Draw("AP");
  gMC->Draw("P SAME");

  TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gData, "Data", "lep");
  leg->AddEntry(gMC, "MB MC", "lep");
  leg->Draw();

  // Save as ROOT file for systematic calculation
  TFile *fOut = TFile::Open(Form("plots/UEFluctuations/%s/UE_DeltaPtRMSvsR.root", saveName), "RECREATE");
  gData->Write("gSigmaDeltaPt_Data");
  gMC->Write("gSigmaDeltaPt_MC");
  fOut->Close();

  c->SaveAs(Form("plots/UEFluctuations/%s/UE_DeltaPtRMSvsR.pdf", saveName));

  std::cout << "========================================" << std::endl;
  std::cout << "UE Fluctuation RMS (for systematic)" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Data:" << std::endl;
  std::cout << "R\tRMS (GeV/c)" << std::endl;
  for (size_t i = 0; i < rValsData.size(); ++i) {
    std::cout << Form("%.1f\t%.4f", rValsData[i], sigmaData[i]) << std::endl;
  }
  std::cout << "MC:" << std::endl;
  std::cout << "R\tRMS (GeV/c)" << std::endl;
  for (size_t i = 0; i < rValsMC.size(); ++i) {
    std::cout << Form("%.1f\t%.4f", rValsMC[i], sigmaMC[i]) << std::endl;
  }
}

void DrawDeltaPtMethodComparison(int dataset) {
  gStyle->SetOptStat(0);

  auto configs = GetUEConfigs(dataset);
  const char* saveName = kPlotSaveNames[dataset];

  // Focus on R=0.4 for method comparison
  double targetR = 0.4;
  UEConfig *targetCfg = nullptr;
  for (auto &cfg : configs) {
    if (fabs(cfg.R - targetR) < 0.01) {
      targetCfg = &cfg;
      break;
    }
  }
  if (!targetCfg) return;

  TCanvas *c = new TCanvas(Form("cDeltaPtMethods_%d", dataset),
      Form("Method comparison: #delta p_{T} (R=%.1f)", targetR), 1200, 600);
  c->Divide(2, 1);

  // Data
  c->cd(1);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.05);

  TLegend *legData = new TLegend(0.55, 0.60, 0.88, 0.88);
  legData->SetBorderSize(0);
  legData->SetFillStyle(0);
  legData->SetTextSize(0.04);

  bool first = true;

  if (targetCfg->dataFile.Length() > 0) {
    TFile *fData = TFile::Open(targetCfg->dataFile, "READ");
    if (fData && !fData->IsZombie()) {
      // Standard Random Cone
      TH2 *h2_rc = (TH2*)fData->Get(Form("%s/h2_centrality_rhorandomcone", kDirName));
      if (h2_rc) {
        TH1 *h = ProjectY(h2_rc, Form("hRC_Data_%d", dataset));
        h->SetLineColor(kBlack);
        h->SetLineWidth(1);
        if (h->Integral() > 0) h->Scale(1.0 / h->Integral(), "width");
        h->GetXaxis()->SetRangeUser(-10, 50);
        h->GetXaxis()->SetTitle("#delta p_{T} (GeV/c)");
        h->GetYaxis()->SetTitle("(1/N) dN/d#delta p_{T}");
        h->GetYaxis()->SetTitleOffset(1.3);
        h->SetTitle("Data (R=0.4)");
        h->Draw("HIST");
        first = false;
        legData->AddEntry(h, "Random Cone", "l");
      }

      // Randomised track direction
      TH2 *h2_rand = (TH2*)fData->Get(Form("%s/h2_centrality_rhorandomconerandomtrackdirection", kDirName));
      if (h2_rand) {
        TH1 *h = ProjectY(h2_rand, Form("hRand_Data_%d", dataset));
        h->SetLineColor(kGreen+2);
        h->SetLineWidth(1);
        h->SetLineStyle(2);
        if (h->Integral() > 0) h->Scale(1.0 / h->Integral(), "width");
        h->Draw("HIST SAME");
        legData->AddEntry(h, "Randomised tracks", "l");
      }

      // Without leading jet
      TH2 *h2_noLead = (TH2*)fData->Get(Form("%s/h2_centrality_rhorandomconewithoutleadingjet", kDirName));
      if (h2_noLead) {
        TH1 *h = ProjectY(h2_noLead, Form("hNoLead_Data_%d", dataset));
        h->SetLineColor(kRed+1);
        h->SetLineWidth(1);
        h->SetLineStyle(3);
        if (h->Integral() > 0) h->Scale(1.0 / h->Integral(), "width");
        h->Draw("HIST SAME");
        legData->AddEntry(h, "RC w/o leading jet", "l");
      }

      legData->Draw();
      fData->Close();
    }
  }

  // MC
  c->cd(2);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.05);

  TLegend *legMC = new TLegend(0.55, 0.60, 0.88, 0.88);
  legMC->SetBorderSize(0);
  legMC->SetFillStyle(0);
  legMC->SetTextSize(0.04);

  first = true;

  if (targetCfg->mcFile.Length() > 0) {
    TFile *fMC = TFile::Open(targetCfg->mcFile, "READ");
    if (fMC && !fMC->IsZombie()) {
      TH2 *h2_rc = (TH2*)fMC->Get(Form("%s/h2_centrality_rhorandomcone", kDirName));
      if (h2_rc) {
        TH1 *h = ProjectY(h2_rc, Form("hRC_MC_%d", dataset));
        h->SetLineColor(kBlack);
        h->SetLineWidth(1);
        if (h->Integral() > 0) h->Scale(1.0 / h->Integral(), "width");
        h->GetXaxis()->SetRangeUser(-10, 50);
        h->GetXaxis()->SetTitle("#delta p_{T} (GeV/c)");
        h->GetYaxis()->SetTitle("(1/N) dN/d#delta p_{T}");
        h->GetYaxis()->SetTitleOffset(1.3);
        h->SetTitle("MB MC (R=0.4)");
        h->Draw("HIST");
        first = false;
        legMC->AddEntry(h, "Random Cone", "l");
      }

      TH2 *h2_rand = (TH2*)fMC->Get(Form("%s/h2_centrality_rhorandomconerandomtrackdirection", kDirName));
      if (h2_rand) {
        TH1 *h = ProjectY(h2_rand, Form("hRand_MC_%d", dataset));
        h->SetLineColor(kGreen+2);
        h->SetLineWidth(1);
        h->SetLineStyle(2);
        if (h->Integral() > 0) h->Scale(1.0 / h->Integral(), "width");
        h->Draw("HIST SAME");
        legMC->AddEntry(h, "Randomised tracks", "l");
      }

      TH2 *h2_noLead = (TH2*)fMC->Get(Form("%s/h2_centrality_rhorandomconewithoutleadingjet", kDirName));
      if (h2_noLead) {
        TH1 *h = ProjectY(h2_noLead, Form("hNoLead_MC_%d", dataset));
        h->SetLineColor(kRed+1);
        h->SetLineWidth(1);
        h->SetLineStyle(3);
        if (h->Integral() > 0) h->Scale(1.0 / h->Integral(), "width");
        h->Draw("HIST SAME");
        legMC->AddEntry(h, "RC w/o leading jet", "l");
      }

      legMC->Draw();
      fMC->Close();
    }
  }

  c->SaveAs(Form("plots/UEFluctuations/%s/UE_DeltaPtMethodComparison.pdf", saveName));
}

void DrawNtracksVsRho(int dataset) {
  gStyle->SetOptStat(0);

  auto configs = GetUEConfigs(dataset);
  const char* saveName = kPlotSaveNames[dataset];

  // h2_ntracks_rho is filled by processRho, which is R-independent.
  // Use R=0.4 as representative.
  UEConfig *repCfg = nullptr;
  for (auto &cfg : configs) {
    if (fabs(cfg.R - 0.4) < 0.01) { repCfg = &cfg; break; }
  }
  if (!repCfg) return;

  TCanvas *c = new TCanvas(Form("cNtracksRho_%d", dataset),
      "N_{tracks} vs #rho", 1200, 600);
  c->Divide(2, 1);

  // Data
  c->cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.15);
  gPad->SetLogz();

  if (repCfg->dataFile.Length() > 0) {
    TFile *fData = TFile::Open(repCfg->dataFile, "READ");
    if (fData && !fData->IsZombie()) {
      TH2 *h2 = (TH2*)fData->Get(Form("%s/h2_ntracks_rho", kDirName));
      if (h2) {
        TH2 *h2c = (TH2*)h2->Clone(Form("h2_ntracks_rho_Data_%d", dataset));
        h2c->SetDirectory(0);
        h2c->SetTitle("Data: N_{tracks} vs #rho");
        h2c->GetXaxis()->SetTitle("N_{tracks}");
        h2c->GetXaxis()->SetRangeUser(0, 200);
        h2c->GetYaxis()->SetTitle("#rho (GeV/area)");
        h2c->GetYaxis()->SetRangeUser(0, 20);
        h2c->GetYaxis()->SetTitleOffset(1.3);
        h2c->Draw("COLZ");
      }
      fData->Close();
    }
  }

  // MC
  c->cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.15);
  gPad->SetLogz();

  if (repCfg->mcFile.Length() > 0) {
    TFile *fMC = TFile::Open(repCfg->mcFile, "READ");
    if (fMC && !fMC->IsZombie()) {
      TH2 *h2 = (TH2*)fMC->Get(Form("%s/h2_ntracks_rho", kDirName));
      if (h2) {
        TH2 *h2c = (TH2*)h2->Clone(Form("h2_ntracks_rho_MC_%d", dataset));
        h2c->SetDirectory(0);
        h2c->SetTitle("MB MC: N_{tracks} vs #rho");
        h2c->GetXaxis()->SetTitle("N_{tracks}");
        h2c->GetXaxis()->SetRangeUser(0, 200);
        h2c->GetYaxis()->SetTitle("#rho (GeV/area)");
        h2c->GetYaxis()->SetRangeUser(0, 20);
        h2c->GetYaxis()->SetTitleOffset(1.3);
        h2c->Draw("COLZ");
      }
      fMC->Close();
    }
  }

  c->SaveAs(Form("plots/UEFluctuations/%s/UE_NtracksVsRho.pdf", saveName));

  // --- Profile: <rho> vs Ntracks ---
  TCanvas *cProf = new TCanvas(Form("cNtracksRhoProf_%d", dataset),
      "N_{tracks} vs <#rho_{UE}>", 800, 600);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.05);

  TLegend *leg = new TLegend(0.15, 0.70, 0.45, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);

  bool firstDrawn = false;

  // Data profile
  if (repCfg->dataFile.Length() > 0) {
    TFile *fData = TFile::Open(repCfg->dataFile, "READ");
    if (fData && !fData->IsZombie()) {
      TH2 *h2 = (TH2*)fData->Get(Form("%s/h2_ntracks_rho", kDirName));
      if (h2) {
        TProfile *prof = h2->ProfileX(Form("profRho_Data_%d", dataset));
        prof->SetDirectory(0);
        prof->SetLineColor(kBlack);
        prof->SetMarkerColor(kBlack);
        prof->SetMarkerStyle(20);
        prof->SetMarkerSize(0.8);
        prof->GetXaxis()->SetRangeUser(0, 200);
        prof->GetYaxis()->SetRangeUser(0, 10);
        prof->GetXaxis()->SetTitle("N_{tracks}");
        prof->GetYaxis()->SetTitle("<#rho_{UE}> (GeV/area)");
        prof->GetYaxis()->SetTitleOffset(1.3);
        prof->SetTitle("N_{tracks} vs <#rho_{UE}>");
        prof->Draw("E");
        leg->AddEntry(prof, "Data", "lep");
        firstDrawn = true;
      }
      fData->Close();
    }
  }

  // MC profile
  if (repCfg->mcFile.Length() > 0) {
    TFile *fMC = TFile::Open(repCfg->mcFile, "READ");
    if (fMC && !fMC->IsZombie()) {
      TH2 *h2 = (TH2*)fMC->Get(Form("%s/h2_ntracks_rho", kDirName));
      if (h2) {
        TProfile *prof = h2->ProfileX(Form("profRho_MC_%d", dataset));
        prof->SetDirectory(0);
        prof->SetLineColor(kRed+1);
        prof->SetMarkerColor(kRed+1);
        prof->SetMarkerStyle(24);
        prof->SetMarkerSize(0.8);
        if (!firstDrawn) {
          prof->GetXaxis()->SetRangeUser(0, 200);
          prof->GetYaxis()->SetRangeUser(0, 10);
          prof->GetXaxis()->SetTitle("N_{tracks}");
          prof->GetYaxis()->SetTitle("<#rho_{UE}> (GeV/area)");
          prof->GetYaxis()->SetTitleOffset(1.3);
          prof->SetTitle("N_{tracks} vs <#rho_{UE}>");
          prof->Draw("E");
        } else {
          prof->Draw("E SAME");
        }
        leg->AddEntry(prof, "MB MC", "lep");
      }
      fMC->Close();
    }
  }

  leg->Draw();
  cProf->SaveAs(Form("plots/UEFluctuations/%s/UE_NtracksVsMeanRho.pdf", saveName));
}

void DrawDeltaPtVsMultiplicity(int dataset) {
  gStyle->SetOptStat(0);

  auto configs = GetUEConfigs(dataset);
  const char* saveName = kPlotSaveNames[dataset];

  // Use R=0.4 as representative (delta pT uses fixed randomConeR, R-independent)
  UEConfig *repCfg = nullptr;
  for (auto &cfg : configs) {
    if (fabs(cfg.R - 0.4) < 0.01) { repCfg = &cfg; break; }
  }
  if (!repCfg) return;

  // Multiplicity (centrality percentile) bin edges
  const int nBins = 9;
  const double multBins[nBins + 1] = {0, 1, 5, 10, 20, 30, 50, 70, 100};

  // Three RC-based methods
  const int nMethods = 3;
  const char* histNames[nMethods] = {
    "h2_centrality_rhorandomcone",
    "h2_centrality_rhorandomconewithoutleadingjet",
    "h2_centrality_rhorandomconerandomtrackdirection"
  };
  const char* methodLabels[nMethods] = {
    "Random Cone",
    "RC w/o leading jet",
    "Randomised tracks"
  };
  const char* methodTags[nMethods] = {"RC", "RCnoLead", "RCrand"};

  for (int m = 0; m < nMethods; ++m) {
    // One canvas per method: 3x3 grid, each pad = one multiplicity class, Data vs MC
    TCanvas *c = new TCanvas(Form("cDptMult_%s_%d", methodTags[m], dataset),
        Form("#delta p_{T} vs mult: %s", methodLabels[m]), 1200, 1200);
    c->Divide(3, 3, 0.001, 0.001);

    // Open files once per method
    TFile *fData = nullptr;
    TFile *fMC = nullptr;
    TH2 *h2Data = nullptr;
    TH2 *h2MC = nullptr;

    if (repCfg->dataFile.Length() > 0) {
      fData = TFile::Open(repCfg->dataFile, "READ");
      if (fData && !fData->IsZombie())
        h2Data = (TH2*)fData->Get(Form("%s/%s", kDirName, histNames[m]));
    }
    if (repCfg->mcFile.Length() > 0) {
      fMC = TFile::Open(repCfg->mcFile, "READ");
      if (fMC && !fMC->IsZombie())
        h2MC = (TH2*)fMC->Get(Form("%s/%s", kDirName, histNames[m]));
    }

    for (int i = 0; i < nBins; ++i) {
      c->cd(i + 1);
      gPad->SetLogy();
      gPad->SetLeftMargin(0.14);
      gPad->SetRightMargin(0.03);
      gPad->SetTopMargin(0.08);
      gPad->SetBottomMargin(0.12);

      bool drawn = false;

      // Data projection
      if (h2Data) {
        int binLo = h2Data->GetXaxis()->FindBin(multBins[i] + 0.001);
        int binHi = h2Data->GetXaxis()->FindBin(multBins[i + 1] - 0.001);
        TH1 *hD = h2Data->ProjectionY(Form("hDptM_%s_D_%d_%d", methodTags[m], dataset, i), binLo, binHi);
        hD->SetDirectory(0);
        hD->SetLineColor(kBlack);
        hD->SetLineWidth(1);
        if (hD->Integral() > 0) hD->Scale(1.0 / hD->Integral(), "width");
        hD->GetXaxis()->SetRangeUser(-10, 50);
        hD->GetXaxis()->SetTitle("#delta p_{T} (GeV/c)");
        hD->GetXaxis()->SetTitleSize(0.05);
        hD->GetXaxis()->SetLabelSize(0.045);
        hD->GetYaxis()->SetTitle("(1/N) dN/d#delta p_{T}");
        hD->GetYaxis()->SetTitleSize(0.05);
        hD->GetYaxis()->SetTitleOffset(1.2);
        hD->GetYaxis()->SetLabelSize(0.045);
        hD->SetTitle(Form("%.0f-%.0f%%", multBins[i], multBins[i + 1]));
        hD->Draw("HIST");
        drawn = true;
      }

      // MC projection
      if (h2MC) {
        int binLo = h2MC->GetXaxis()->FindBin(multBins[i] + 0.001);
        int binHi = h2MC->GetXaxis()->FindBin(multBins[i + 1] - 0.001);
        TH1 *hM = h2MC->ProjectionY(Form("hDptM_%s_M_%d_%d", methodTags[m], dataset, i), binLo, binHi);
        hM->SetDirectory(0);
        hM->SetLineColor(kRed+1);
        hM->SetLineWidth(1);
        if (hM->Integral() > 0) hM->Scale(1.0 / hM->Integral(), "width");
        if (drawn) {
          hM->Draw("HIST SAME");
        } else {
          hM->GetXaxis()->SetRangeUser(-10, 50);
          hM->GetXaxis()->SetTitle("#delta p_{T} (GeV/c)");
          hM->GetXaxis()->SetTitleSize(0.05);
          hM->GetXaxis()->SetLabelSize(0.045);
          hM->GetYaxis()->SetTitle("(1/N) dN/d#delta p_{T}");
          hM->GetYaxis()->SetTitleSize(0.05);
          hM->GetYaxis()->SetTitleOffset(1.2);
          hM->GetYaxis()->SetLabelSize(0.045);
          hM->SetTitle(Form("%.0f-%.0f%%", multBins[i], multBins[i + 1]));
          hM->Draw("HIST");
          drawn = true;
        }
      }

      // Legend only in first pad
      if (i == 0) {
        TLegend *leg = new TLegend(0.50, 0.72, 0.93, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.045);
        leg->SetHeader(methodLabels[m]);
        // Use existing projected histograms for correct legend colors
        TH1 *hLegD = (TH1*)gDirectory->Get(Form("hDptM_%s_D_%d_0", methodTags[m], dataset));
        TH1 *hLegM = (TH1*)gDirectory->Get(Form("hDptM_%s_M_%d_0", methodTags[m], dataset));
        if (hLegD) leg->AddEntry(hLegD, "Data", "l");
        if (hLegM) leg->AddEntry(hLegM, "MC", "l");
        leg->Draw();
      }
    }

    if (fData) fData->Close();
    if (fMC) fMC->Close();

    c->SaveAs(Form("plots/UEFluctuations/%s/UE_DeltaPtVsMult_%s.pdf", saveName, methodTags[m]));
  }
}

// ============================================================
// Mean rho summary across all datasets
// ============================================================

void DrawMeanRhoSummary() {
  gStyle->SetOptStat(0);

  const TString basePath = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults";
  const char* dirNames[] = {"jet-background-analysis-task", "jet-background-analysis"};
  const int nDirs = 2;

  struct RhoFileInfo {
    const char* id;
    const char* label;
    int color;
  };
  const int nFiles = 5;
  RhoFileInfo files[nFiles] = {
    {"614374", "2022 Data",          kBlack},
    {"594014", "2022 MB MC",         kBlue+1},
    {"593757", "2022 JJ MC",         kRed+1},
    {"496213", "2023 Data",          kGreen+2},
    {"608785", "2023 MB MC",         kMagenta+1}
  };

  // Extract <rho> from each file
  double meanRho[nFiles], rmsRho[nFiles];
  bool found[nFiles];
  for (int i = 0; i < nFiles; ++i) {
    found[i] = false;
    meanRho[i] = 0;
    rmsRho[i] = 0;
    TString path = Form("%s/%s_AnalysisResults.root", basePath.Data(), files[i].id);
    TFile *f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) continue;

    TH2 *h2 = nullptr;
    for (int d = 0; d < nDirs; ++d) {
      h2 = (TH2*)f->Get(Form("%s/h2_centrality_rho", dirNames[d]));
      if (h2) break;
    }
    if (h2) {
      TH1 *hProj = h2->ProjectionY("_tmp_rho_summary");
      meanRho[i] = hProj->GetMean();
      rmsRho[i] = hProj->GetRMS();
      found[i] = true;
      delete hProj;
    }
    f->Close();
  }

  // Print summary
  printf("\n========== Mean Rho Summary ==========\n");
  printf("%-20s  %10s  %10s\n", "Dataset", "<rho>", "RMS(rho)");
  for (int i = 0; i < nFiles; ++i) {
    if (found[i])
      printf("%-20s  %10.4f  %10.4f\n", files[i].label, meanRho[i], rmsRho[i]);
    else
      printf("%-20s  NOT FOUND\n", files[i].label);
  }
  printf("======================================\n\n");

  // Bar chart
  TCanvas *c = new TCanvas("cMeanRho", "<#rho_{UE}> summary", 800, 600);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.18);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.08);
  gPad->SetGridy();

  TH1D *hBar = new TH1D("hMeanRhoBar", ";<#rho_{UE}> (GeV/area)", nFiles, 0, nFiles);
  hBar->SetFillColor(kAzure-4);
  hBar->SetLineColor(kBlack);
  hBar->SetBarWidth(0.6);
  hBar->SetBarOffset(0.2);

  double ymax = 0;
  for (int i = 0; i < nFiles; ++i) {
    hBar->GetXaxis()->SetBinLabel(i + 1, files[i].label);
    if (found[i]) {
      hBar->SetBinContent(i + 1, meanRho[i]);
      hBar->SetBinError(i + 1, rmsRho[i]);
      if (meanRho[i] + rmsRho[i] > ymax) ymax = meanRho[i] + rmsRho[i];
    }
  }
  hBar->GetXaxis()->SetLabelSize(0.05);
  hBar->GetYaxis()->SetRangeUser(0, ymax * 1.3);
  hBar->GetYaxis()->SetTitleOffset(1.2);
  hBar->Draw("BAR E");

  // Add value labels on top of bars
  TLatex latex;
  latex.SetTextSize(0.04);
  latex.SetTextAlign(21); // center-top
  for (int i = 0; i < nFiles; ++i) {
    if (found[i]) {
      latex.DrawLatex(i + 0.5, meanRho[i] + rmsRho[i] + ymax * 0.03,
                      Form("%.3f", meanRho[i]));
    }
  }

  TLatex title;
  title.SetNDC();
  title.SetTextSize(0.04);
  title.DrawLatex(0.15, 0.94, "Mean UE density #LT#rho_{UE}#GT (error bars = RMS)");

  gSystem->MakeDirectory("plots/UEFluctuations");
  c->SaveAs("plots/UEFluctuations/UE_MeanRhoSummary.pdf");
}

// ============================================================
// Main function
// ============================================================

void DrawUEFluctuations() {
  gSystem->MakeDirectory("plots");
  gSystem->MakeDirectory("plots/UEFluctuations");

  DrawMeanRhoSummary();

  for (int ds = 0; ds < kNDatasets; ++ds) {
    std::cout << "========================================" << std::endl;
    std::cout << "Dataset " << ds << ": " << kDatasetLabels[ds] << std::endl;
    std::cout << "========================================" << std::endl;

    gSystem->MakeDirectory(Form("plots/UEFluctuations/%s", kPlotSaveNames[ds]));

    DrawRhoDistribution(ds);
    DrawNtracksVsRho(ds);
    DrawDeltaPtDistribution(ds);
    DrawDeltaPtRMSvsR(ds);
    DrawDeltaPtMethodComparison(ds);
    DrawDeltaPtVsMultiplicity(ds);
  }

  std::cout << "All plots saved to plots/UEFluctuations/ directory" << std::endl;
}
