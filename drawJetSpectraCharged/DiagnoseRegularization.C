// DiagnoseRegularization.C — Scan k and Bayesian iterations to see full spread
// Shows whether current systematics (k±1, SVD vs Bayes) cover the real uncertainty
// Run: root -l -b -q DiagnoseRegularization.C
#include "Filipad2.h"
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include <RooUnfoldBayes.h>
#include <RooUnfoldResponse.h>
#include <RooUnfoldSvd.h>
#include <TSVDUnfold_local.h>

// === Same configuration as DrawUnfoldingSystematicUncertainty.C ===
const TString mainDir = "~/cernbox/workspace/O2Physics/jets/AnalysisResults/";
const TString fRoot = "_AnalysisResults.root";
TString dataFile = mainDir + "498133" + fRoot;
TString mcFile = mainDir + "516969" + fRoot;
const char* dataDir = "jet-spectra-charged";
const char* mcDir = "jet-spectra-charged_Nmax1p5";

const Double_t ptbin[21] = {5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20,
                             25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
const Double_t ptbinGen[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14,
                                16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
const Int_t nptBins = 20;
const Int_t nptBinsGen = 25;

// Kernel corrections (same as main macro)
struct KernelCorrections {
  TH1* hPurity;
  TH1* hEfficiency;
};

KernelCorrections ComputeKC(TH1* hMCReco, TH1* hMCTruth, TH2* h2Matched) {
  KernelCorrections out;
  TH1* hRecoMatched = h2Matched->ProjectionX("hRM_diag");
  TH1* hTruthMatched = h2Matched->ProjectionY("hTM_diag");
  hRecoMatched->SetDirectory(nullptr);
  hTruthMatched->SetDirectory(nullptr);

  out.hPurity = (TH1*)hMCReco->Clone("hPur_diag");
  out.hPurity->SetDirectory(nullptr);
  out.hPurity->Reset();
  for (Int_t i = 1; i <= out.hPurity->GetNbinsX(); ++i) {
    Double_t d = hMCReco->GetBinContent(i);
    Double_t n = hRecoMatched->GetBinContent(i);
    out.hPurity->SetBinContent(i, (d > 0) ? TMath::Min(1.0, TMath::Max(0.0, n / d)) : 0.0);
  }

  out.hEfficiency = (TH1*)hMCTruth->Clone("hEff_diag");
  out.hEfficiency->SetDirectory(nullptr);
  out.hEfficiency->Reset();
  for (Int_t i = 1; i <= out.hEfficiency->GetNbinsX(); ++i) {
    Double_t d = hMCTruth->GetBinContent(i);
    Double_t n = hTruthMatched->GetBinContent(i);
    out.hEfficiency->SetBinContent(i, (d > 0) ? TMath::Min(1.0, TMath::Max(0.0, n / d)) : 0.0);
  }
  return out;
}

TH1* ApplyPurity(const TH1* h, const TH1* p, const char* name) {
  TH1* out = (TH1*)h->Clone(name);
  out->SetDirectory(nullptr);
  for (Int_t i = 1; i <= out->GetNbinsX(); ++i) {
    out->SetBinContent(i, h->GetBinContent(i) * p->GetBinContent(i));
    out->SetBinError(i, h->GetBinError(i) * p->GetBinContent(i));
  }
  return out;
}

TH1* ApplyEffCorr(const TH1* h, const TH1* eff, const char* name) {
  TH1* out = (TH1*)h->Clone(name);
  out->SetDirectory(nullptr);
  for (Int_t i = 1; i <= out->GetNbinsX(); ++i) {
    Double_t e = eff->GetBinContent(i);
    if (e > 0) {
      out->SetBinContent(i, h->GetBinContent(i) / e);
      out->SetBinError(i, h->GetBinError(i) / e);
    } else {
      out->SetBinContent(i, 0);
      out->SetBinError(i, 0);
    }
  }
  return out;
}

void DiagnoseRegularization() {
  gSystem->Exec("mkdir -p DiagnoseRegularization");
  std::ofstream logFile("DiagnoseRegularization/log.txt");
  auto out = [&](const std::string& s) { std::cerr << s; logFile << s; };

  // Load files
  TFile* dF = TFile::Open(dataFile, "READ");
  TFile* mF = TFile::Open(mcFile, "READ");

  TH1* hDataRaw = (TH1*)dF->Get(Form("%s/h_jet_pt", dataDir));
  TH1* hMCRecoRaw = (TH1*)mF->Get(Form("%s/h_jet_pt", mcDir));
  TH1* hMCTruthRaw = (TH1*)mF->Get(Form("%s/h_jet_pt_part", mcDir));
  TH2* hRespRaw = (TH2*)mF->Get(Form("%s/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint", mcDir));

  hDataRaw->SetDirectory(nullptr);
  hMCRecoRaw->SetDirectory(nullptr);
  hMCTruthRaw->SetDirectory(nullptr);
  hRespRaw->SetDirectory(nullptr);

  // Rebin
  TH1* hData = hDataRaw->Rebin(nptBins, "hData_reg", ptbin);
  TH1* hMCReco = hMCRecoRaw->Rebin(nptBins, "hMCReco_reg", ptbin);
  TH1* hMCTruth = hMCTruthRaw->Rebin(nptBinsGen, "hMCTruth_reg", ptbinGen);
  hData->SetDirectory(nullptr);
  hMCReco->SetDirectory(nullptr);
  hMCTruth->SetDirectory(nullptr);

  TH2F* hResp = new TH2F("hResp_reg", "", nptBins, ptbin, nptBinsGen, ptbinGen);
  for (Int_t ix = 1; ix <= hRespRaw->GetNbinsX(); ++ix) {
    for (Int_t iy = 1; iy <= hRespRaw->GetNbinsY(); ++iy) {
      Double_t c = hRespRaw->GetBinContent(ix, iy);
      if (c <= 0) continue;
      Int_t xb = hResp->GetXaxis()->FindBin(hRespRaw->GetXaxis()->GetBinCenter(ix));
      Int_t yb = hResp->GetYaxis()->FindBin(hRespRaw->GetYaxis()->GetBinCenter(iy));
      if (xb < 1 || xb > nptBins || yb < 1 || yb > nptBinsGen) continue;
      hResp->SetBinContent(xb, yb, hResp->GetBinContent(xb, yb) + c);
    }
  }

  // Kernel corrections
  KernelCorrections kc = ComputeKC(hMCReco, hMCTruth, hResp);
  TH1* hDataMatched = ApplyPurity(hData, kc.hPurity, "hDataMatched_reg");

  // Build response
  TH1* hRM = hResp->ProjectionX("hRM_build");
  TH1* hTM = hResp->ProjectionY("hTM_build");
  hRM->SetDirectory(nullptr);
  hTM->SetDirectory(nullptr);
  // 3-arg constructor: directly sets _res=hResp, _mes=hRM, _tru=hTM
  // No Fill() loop — avoids doubling _mes/_tru that occurs with 2-arg + Fill pattern
  // Preserves original bin errors (Sumw2) from the 2D matrix
  RooUnfoldResponse* Response = new RooUnfoldResponse(hRM, hTM, hResp);

  out("=== Regularization Scan ===\n\n");

  // ============================================================
  // 1. SVD scan: k = 2 to 15
  // ============================================================
  out("--- SVD Unfolding: k scan ---\n");
  out("(particle-level yield per bin, after efficiency correction)\n\n");

  // Header: pT bins
  out(Form("%-6s", "k\\pT"));
  for (Int_t b = 6; b <= nptBinsGen; ++b) { // start from bin 6 = [5,6] GeV
    Double_t lo = hMCTruth->GetXaxis()->GetBinLowEdge(b);
    Double_t hi = hMCTruth->GetXaxis()->GetBinUpEdge(b);
    out(Form(" [%3.0f,%3.0f]", lo, hi));
  }
  out("\n");

  std::vector<Int_t> kValues = {2, 3, 4, 5, 6, 7, 8, 10, 12, 15};
  std::map<Int_t, TH1*> svdResults;

  for (Int_t k : kValues) {
    RooUnfoldSvd unfold(Response, hDataMatched, k);
    TH1* hUnf = (TH1*)unfold.Hreco()->Clone(Form("hSVD_k%d", k));
    hUnf->SetDirectory(nullptr);
    TH1* hCorr = ApplyEffCorr(hUnf, kc.hEfficiency, Form("hSVD_k%d_corr", k));
    svdResults[k] = hCorr;

    out(Form("k=%-3d ", k));
    for (Int_t b = 6; b <= nptBinsGen; ++b) {
      out(Form(" %9.3g", hCorr->GetBinContent(b)));
    }
    out("\n");
  }

  // Print ratio to k=5 (nominal)
  out("\n--- SVD ratio to nominal (k=5) ---\n");
  out(Form("%-6s", "k\\pT"));
  for (Int_t b = 6; b <= nptBinsGen; ++b) {
    Double_t lo = hMCTruth->GetXaxis()->GetBinLowEdge(b);
    Double_t hi = hMCTruth->GetXaxis()->GetBinUpEdge(b);
    out(Form("  [%3.0f,%3.0f]", lo, hi));
  }
  out("\n");

  TH1* hNomSVD = svdResults[5];
  for (Int_t k : kValues) {
    TH1* h = svdResults[k];
    out(Form("k=%-3d  ", k));
    for (Int_t b = 6; b <= nptBinsGen; ++b) {
      Double_t nom = hNomSVD->GetBinContent(b);
      Double_t var = h->GetBinContent(b);
      Double_t ratio = (nom > 0) ? var / nom : 0;
      out(Form("  %7.4f", ratio));
    }
    out("\n");
  }

  // ============================================================
  // 2. Bayesian scan: iter = 1 to 20
  // ============================================================
  out("\n--- Bayesian Unfolding: iteration scan ---\n");
  out("--- Ratio to SVD k=5 nominal ---\n");
  out(Form("%-8s", "iter\\pT"));
  for (Int_t b = 6; b <= nptBinsGen; ++b) {
    Double_t lo = hMCTruth->GetXaxis()->GetBinLowEdge(b);
    Double_t hi = hMCTruth->GetXaxis()->GetBinUpEdge(b);
    out(Form("  [%3.0f,%3.0f]", lo, hi));
  }
  out("\n");

  std::vector<Int_t> iterValues = {1, 2, 3, 4, 5, 6, 8, 10, 15, 20};
  std::map<Int_t, TH1*> bayesResults;

  for (Int_t it : iterValues) {
    RooUnfoldBayes unfold(Response, hDataMatched, it);
    TH1* hUnf = (TH1*)unfold.Hreco()->Clone(Form("hBayes_i%d", it));
    hUnf->SetDirectory(nullptr);
    TH1* hCorr = ApplyEffCorr(hUnf, kc.hEfficiency, Form("hBayes_i%d_corr", it));
    bayesResults[it] = hCorr;

    out(Form("iter=%-3d", it));
    for (Int_t b = 6; b <= nptBinsGen; ++b) {
      Double_t nom = hNomSVD->GetBinContent(b);
      Double_t var = hCorr->GetBinContent(b);
      Double_t ratio = (nom > 0) ? var / nom : 0;
      out(Form("  %7.4f", ratio));
    }
    out("\n");
  }

  // ============================================================
  // 3. Summary: max deviation from nominal per pT bin
  // ============================================================
  out("\n=== Summary: Max |deviation| from SVD k=5 per pT bin (%) ===\n");
  out(Form("%-30s", "Source"));
  for (Int_t b = 6; b <= nptBinsGen; ++b) {
    Double_t lo = hMCTruth->GetXaxis()->GetBinLowEdge(b);
    Double_t hi = hMCTruth->GetXaxis()->GetBinUpEdge(b);
    out(Form("  [%3.0f,%3.0f]", lo, hi));
  }
  out("\n");

  // Current systematics
  auto printMaxDev = [&](const char* label, std::vector<TH1*> vars) {
    out(Form("%-30s", label));
    for (Int_t b = 6; b <= nptBinsGen; ++b) {
      Double_t nom = hNomSVD->GetBinContent(b);
      Double_t maxDev = 0;
      for (auto* h : vars) {
        Double_t v = h->GetBinContent(b);
        if (nom > 0) maxDev = TMath::Max(maxDev, TMath::Abs(v / nom - 1.0));
      }
      out(Form("  %7.1f%%", maxDev * 100));
    }
    out("\n");
  };

  // Current: k±1
  printMaxDev("Current: SVD k±1 (k=4,6)", {svdResults[4], svdResults[6]});

  // Current: method (SVD k=5 vs Bayes iter=4)
  printMaxDev("Current: Bayes iter=4", {bayesResults[4]});

  // What SHOULD be covered: full k range
  printMaxDev("Full SVD range: k=2..15", {svdResults[2], svdResults[3], svdResults[4],
      svdResults[6], svdResults[7], svdResults[8], svdResults[10], svdResults[12], svdResults[15]});

  // Full Bayes range
  printMaxDev("Full Bayes range: iter=1..20", {bayesResults[1], bayesResults[2], bayesResults[3],
      bayesResults[4], bayesResults[5], bayesResults[6], bayesResults[8], bayesResults[10],
      bayesResults[15], bayesResults[20]});

  // Physically motivated: SVD k=5..10 (closure to halfway to data)
  printMaxDev("SVD k=5..10 (conservative)", {svdResults[6], svdResults[7], svdResults[8], svdResults[10]});

  // Bayesian iter=4..10
  printMaxDev("Bayes iter=4..10 (conservative)", {bayesResults[4], bayesResults[5], bayesResults[6],
      bayesResults[8], bayesResults[10]});

  // ============================================================
  // 4. Plot: ratio to nominal for all k and iter values
  // ============================================================
  Int_t colors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kOrange+1, kCyan+1, kViolet, kYellow+2, kGray+1};

  // SVD ratio plot
  TCanvas* cSVD = new TCanvas("cSVD", "SVD k scan", 1200, 600);
  cSVD->SetGridy();
  TH1* hFrame = (TH1*)hNomSVD->Clone("hFrame_svd");
  hFrame->Reset();
  hFrame->SetDirectory(nullptr);
  hFrame->GetXaxis()->SetRangeUser(5, 200);
  hFrame->GetYaxis()->SetRangeUser(0.5, 1.5);
  hFrame->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/#it{c})");
  hFrame->GetYaxis()->SetTitle("Ratio to SVD k=5");
  hFrame->SetTitle("SVD regularization scan");
  hFrame->Draw();

  TLine* line1 = new TLine(5, 1.0, 200, 1.0);
  line1->SetLineStyle(2);
  line1->Draw("same");

  TLegend* legSVD = new TLegend(0.12, 0.55, 0.35, 0.88);
  legSVD->SetBorderSize(0);
  legSVD->SetFillColorAlpha(0, 0);
  legSVD->SetTextSize(0.03);

  for (Int_t ik = 0; ik < (Int_t)kValues.size(); ++ik) {
    Int_t k = kValues[ik];
    TH1* hRatio = (TH1*)svdResults[k]->Clone(Form("hRatio_svd_k%d", k));
    hRatio->SetDirectory(nullptr);
    for (Int_t b = 1; b <= hRatio->GetNbinsX(); ++b) {
      Double_t nom = hNomSVD->GetBinContent(b);
      if (nom > 0) {
        hRatio->SetBinContent(b, svdResults[k]->GetBinContent(b) / nom);
        hRatio->SetBinError(b, 0);
      } else {
        hRatio->SetBinContent(b, 0);
      }
    }
    hRatio->SetLineColor(colors[ik % 10]);
    hRatio->SetLineWidth(k == 5 ? 3 : 2);
    hRatio->SetLineStyle(k == 5 ? 1 : 1);
    hRatio->SetMarkerColor(colors[ik % 10]);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.6);
    hRatio->GetXaxis()->SetRangeUser(5, 200);
    hRatio->Draw("hist same");
    legSVD->AddEntry(hRatio, Form("k=%d", k), "l");
  }
  legSVD->Draw("same");
  cSVD->SaveAs("DiagnoseRegularization/SVD_k_scan.pdf");

  // Bayesian ratio plot
  TCanvas* cBayes = new TCanvas("cBayes", "Bayes iter scan", 1200, 600);
  cBayes->SetGridy();
  TH1* hFrame2 = (TH1*)hFrame->Clone("hFrame_bayes");
  hFrame2->SetTitle("Bayesian iteration scan");
  hFrame2->Draw();
  line1->Draw("same");

  TLegend* legBayes = new TLegend(0.12, 0.55, 0.35, 0.88);
  legBayes->SetBorderSize(0);
  legBayes->SetFillColorAlpha(0, 0);
  legBayes->SetTextSize(0.03);

  for (Int_t ii = 0; ii < (Int_t)iterValues.size(); ++ii) {
    Int_t it = iterValues[ii];
    TH1* hRatio = (TH1*)bayesResults[it]->Clone(Form("hRatio_bayes_i%d", it));
    hRatio->SetDirectory(nullptr);
    for (Int_t b = 1; b <= hRatio->GetNbinsX(); ++b) {
      Double_t nom = hNomSVD->GetBinContent(b);
      if (nom > 0) {
        hRatio->SetBinContent(b, bayesResults[it]->GetBinContent(b) / nom);
        hRatio->SetBinError(b, 0);
      } else {
        hRatio->SetBinContent(b, 0);
      }
    }
    hRatio->SetLineColor(colors[ii % 10]);
    hRatio->SetLineWidth(it == 4 ? 3 : 2);
    hRatio->SetMarkerColor(colors[ii % 10]);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.6);
    hRatio->GetXaxis()->SetRangeUser(5, 200);
    hRatio->Draw("hist same");
    legBayes->AddEntry(hRatio, Form("iter=%d", it), "l");
  }
  legBayes->Draw("same");
  cBayes->SaveAs("DiagnoseRegularization/Bayes_iter_scan.pdf");

  logFile.close();
  std::cerr << "\n[Info] Log: DiagnoseRegularization/log.txt" << std::endl;
  std::cerr << "[Info] Plots: DiagnoseRegularization/SVD_k_scan.pdf, Bayes_iter_scan.pdf" << std::endl;
}
