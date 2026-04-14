// DiagnoseDataMCMismatch.C — Quick diagnostic to understand data-MC tension
// Run: root -l -b -q DiagnoseDataMCMismatch.C
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <iomanip>

void DiagnoseDataMCMismatch() {
  // === Same configuration as DrawUnfoldingSystematicUncertainty.C ===
  const TString mainDir = "~/cernbox/workspace/O2Physics/jets/AnalysisResults/";
  const TString fRoot = "_AnalysisResults.root";

  TString dataFile = mainDir + "498133" + fRoot;
  TString mcFile = mainDir + "516969" + fRoot;
  const char* dataDir = "jet-spectra-charged";
  const char* mcDir = "jet-spectra-charged_Nmax1p5";

  // Same binning
  const Double_t ptbin[21] = {5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20,
                               25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
  const Double_t ptbinGen[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14,
                                  16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
  const Int_t nptBins = 20;
  const Int_t nptBinsGen = 25;

  std::ofstream log("DiagnoseDataMCMismatch_log.txt");
  auto out = [&](const std::string& s) { std::cerr << s; log << s; };

  // ============================================================
  // 1. Open files, check basic histogram existence
  // ============================================================
  TFile* dF = TFile::Open(dataFile, "READ");
  TFile* mF = TFile::Open(mcFile, "READ");
  if (!dF || dF->IsZombie()) { out("[FATAL] Cannot open data file\n"); return; }
  if (!mF || mF->IsZombie()) { out("[FATAL] Cannot open MC file\n"); return; }

  out("=== Files opened successfully ===\n");
  out(Form("Data: %s / %s\n", dataFile.Data(), dataDir));
  out(Form("MC:   %s / %s\n", mcFile.Data(), mcDir));

  // ============================================================
  // 2. Load raw histograms (before rebinning)
  // ============================================================
  TH1* hDataRaw = (TH1*)dF->Get(Form("%s/h_jet_pt", dataDir));
  TH1* hMCRecoRaw = (TH1*)mF->Get(Form("%s/h_jet_pt", mcDir));
  TH1* hMCTruthRaw = (TH1*)mF->Get(Form("%s/h_jet_pt_part", mcDir));
  TH2* hRespRaw = (TH2*)mF->Get(Form("%s/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint", mcDir));

  // Check weighted versions
  TH1* hMCRecoWeighted = (TH1*)mF->Get(Form("%s/h_jet_pt_weighted", mcDir));
  TH1* hMCTruthWeighted = (TH1*)mF->Get(Form("%s/h_jet_pt_part_weighted", mcDir));
  TH2* hRespWeighted = (TH2*)mF->Get(Form("%s/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint_weighted", mcDir));

  out("\n=== Histogram availability ===\n");
  out(Form("  h_jet_pt (data):     %s  (entries=%.0f)\n", hDataRaw ? "YES" : "NO", hDataRaw ? hDataRaw->GetEntries() : 0));
  out(Form("  h_jet_pt (MC):       %s  (entries=%.0f)\n", hMCRecoRaw ? "YES" : "NO", hMCRecoRaw ? hMCRecoRaw->GetEntries() : 0));
  out(Form("  h_jet_pt_part (MC):  %s  (entries=%.0f)\n", hMCTruthRaw ? "YES" : "NO", hMCTruthRaw ? hMCTruthRaw->GetEntries() : 0));
  out(Form("  matched 2D (MC):     %s  (entries=%.0f)\n", hRespRaw ? "YES" : "NO", hRespRaw ? hRespRaw->GetEntries() : 0));
  out(Form("  h_jet_pt_weighted:           %s\n", hMCRecoWeighted ? "YES" : "NO"));
  out(Form("  h_jet_pt_part_weighted:      %s\n", hMCTruthWeighted ? "YES" : "NO"));
  out(Form("  matched 2D weighted:         %s\n", hRespWeighted ? "YES" : "NO"));

  // Check h_mccollisions
  TH1* hMCColl = (TH1*)mF->Get(Form("%s/h_mccollisions", mcDir));
  TH1* hMCCollW = (TH1*)mF->Get(Form("%s/h_mccollisions_weighted", mcDir));
  if (hMCColl) {
    out("\n=== h_mccollisions (unweighted) ===\n");
    for (Int_t i = 1; i <= hMCColl->GetNbinsX(); ++i)
      out(Form("  Bin %d (%.1f): %.0f\n", i, hMCColl->GetBinCenter(i), hMCColl->GetBinContent(i)));
  }
  if (hMCCollW) {
    out("\n=== h_mccollisions_weighted ===\n");
    for (Int_t i = 1; i <= hMCCollW->GetNbinsX(); ++i)
      out(Form("  Bin %d (%.1f): %.6g\n", i, hMCCollW->GetBinCenter(i), hMCCollW->GetBinContent(i)));
  }

  // Check data collisions
  TH1* hDataColl = (TH1*)dF->Get(Form("%s/h_collisions", dataDir));
  if (hDataColl) {
    out("\n=== h_collisions (data) ===\n");
    for (Int_t i = 1; i <= hDataColl->GetNbinsX(); ++i)
      out(Form("  Bin %d (%.1f): %.0f\n", i, hDataColl->GetBinCenter(i), hDataColl->GetBinContent(i)));
  }

  if (!hDataRaw || !hMCRecoRaw || !hMCTruthRaw || !hRespRaw) {
    out("[FATAL] Missing essential histograms\n");
    return;
  }

  // ============================================================
  // 3. Raw histogram ranges and binning
  // ============================================================
  out("\n=== Raw histogram binning ===\n");
  out(Form("  Data h_jet_pt:   nBins=%d, xMin=%.1f, xMax=%.1f, integral=%.6g\n",
      hDataRaw->GetNbinsX(), hDataRaw->GetXaxis()->GetXmin(), hDataRaw->GetXaxis()->GetXmax(), hDataRaw->Integral()));
  out(Form("  MC h_jet_pt:     nBins=%d, xMin=%.1f, xMax=%.1f, integral=%.6g\n",
      hMCRecoRaw->GetNbinsX(), hMCRecoRaw->GetXaxis()->GetXmin(), hMCRecoRaw->GetXaxis()->GetXmax(), hMCRecoRaw->Integral()));
  out(Form("  MC h_jet_pt_part: nBins=%d, xMin=%.1f, xMax=%.1f, integral=%.6g\n",
      hMCTruthRaw->GetNbinsX(), hMCTruthRaw->GetXaxis()->GetXmin(), hMCTruthRaw->GetXaxis()->GetXmax(), hMCTruthRaw->Integral()));
  out(Form("  MC matched 2D:   nBinsX=%d, nBinsY=%d, integral=%.6g\n",
      hRespRaw->GetNbinsX(), hRespRaw->GetNbinsY(), hRespRaw->Integral()));

  // ============================================================
  // 4. Rebin and compare detector-level shapes
  // ============================================================
  hDataRaw->SetDirectory(nullptr);
  hMCRecoRaw->SetDirectory(nullptr);
  hMCTruthRaw->SetDirectory(nullptr);
  hRespRaw->SetDirectory(nullptr);

  TH1* hData = hDataRaw->Rebin(nptBins, "hData_diag", ptbin);
  TH1* hMCReco = hMCRecoRaw->Rebin(nptBins, "hMCReco_diag", ptbin);
  TH1* hMCTruth = hMCTruthRaw->Rebin(nptBinsGen, "hMCTruth_diag", ptbinGen);

  // Rebin 2D
  TH2F* hResp = new TH2F("hResp_diag", "hResp_diag", nptBins, ptbin, nptBinsGen, ptbinGen);
  for (Int_t ix = 1; ix <= hRespRaw->GetNbinsX(); ++ix) {
    for (Int_t iy = 1; iy <= hRespRaw->GetNbinsY(); ++iy) {
      Double_t content = hRespRaw->GetBinContent(ix, iy);
      if (content <= 0) continue;
      Int_t xbin = hResp->GetXaxis()->FindBin(hRespRaw->GetXaxis()->GetBinCenter(ix));
      Int_t ybin = hResp->GetYaxis()->FindBin(hRespRaw->GetYaxis()->GetBinCenter(iy));
      if (xbin < 1 || xbin > nptBins || ybin < 1 || ybin > nptBinsGen) continue;
      hResp->SetBinContent(xbin, ybin, hResp->GetBinContent(xbin, ybin) + content);
    }
  }

  // Matched projections
  TH1* hRecoMatched = hResp->ProjectionX("hRecoMatched_diag");
  TH1* hTruthMatched = hResp->ProjectionY("hTruthMatched_diag");

  out("\n=== Detector-level comparison (rebinned) ===\n");
  out("  Data vs MC reco shape comparison (normalized to same integral in 5-200 GeV):\n");

  Double_t dataIntegral = hData->Integral();
  Double_t mcRecoIntegral = hMCReco->Integral();

  out(Form("  Data integral:    %.6g\n", dataIntegral));
  out(Form("  MC reco integral: %.6g\n", mcRecoIntegral));
  out(Form("  Ratio (Data/MC):  %.4f\n\n", dataIntegral / mcRecoIntegral));

  out(Form("  %-12s %-14s %-14s %-14s %-10s %-10s %-10s\n",
      "pT range", "Data", "MC reco", "MC matched", "Data/MC", "Purity", "Efficiency"));
  out("  " + std::string(90, '-') + "\n");

  for (Int_t i = 1; i <= nptBins; ++i) {
    Double_t lo = hData->GetXaxis()->GetBinLowEdge(i);
    Double_t hi = hData->GetXaxis()->GetBinUpEdge(i);
    Double_t d = hData->GetBinContent(i);
    Double_t mc = hMCReco->GetBinContent(i);
    Double_t matched = hRecoMatched->GetBinContent(i);
    Double_t ratio = (mc > 0) ? (d / dataIntegral) / (mc / mcRecoIntegral) : 0;
    Double_t purity = (mc > 0) ? matched / mc : 0;

    // Find corresponding truth bin for efficiency
    Int_t truthBin = hMCTruth->GetXaxis()->FindBin((lo + hi) / 2.0);
    Double_t truth = hMCTruth->GetBinContent(truthBin);
    Double_t truthMatched = hTruthMatched->GetBinContent(truthBin);
    Double_t efficiency = (truth > 0) ? truthMatched / truth : 0;

    out(Form("  [%3.0f,%3.0f]    %14.6g %14.6g %14.6g %10.4f %10.4f %10.4f\n",
        lo, hi, d, mc, matched, ratio, purity, efficiency));
  }

  // ============================================================
  // 5. Truth-level: MC truth and efficiency for low-pT bins (0-5 GeV)
  // ============================================================
  out("\n=== Truth-level bins (including sub-reco low-pT) ===\n");
  out(Form("  %-12s %-14s %-14s %-10s\n", "pT range", "MC truth", "MC matched", "Efficiency"));
  out("  " + std::string(50, '-') + "\n");
  for (Int_t i = 1; i <= nptBinsGen; ++i) {
    Double_t lo = hMCTruth->GetXaxis()->GetBinLowEdge(i);
    Double_t hi = hMCTruth->GetXaxis()->GetBinUpEdge(i);
    Double_t truth = hMCTruth->GetBinContent(i);
    Double_t matched = hTruthMatched->GetBinContent(i);
    Double_t eff = (truth > 0) ? matched / truth : 0;
    out(Form("  [%3.0f,%3.0f]    %14.6g %14.6g %10.6f\n", lo, hi, truth, matched, eff));
  }

  // ============================================================
  // 6. Response matrix diagonal strength
  // ============================================================
  out("\n=== Response matrix: diagonal fraction per reco bin ===\n");
  out("  (fraction of reco-bin entries on diagonal truth bin)\n");
  for (Int_t i = 1; i <= nptBins; ++i) {
    Double_t lo = hResp->GetXaxis()->GetBinLowEdge(i);
    Double_t hi = hResp->GetXaxis()->GetBinUpEdge(i);
    Double_t rowTotal = 0;
    Double_t diagContent = 0;
    for (Int_t j = 1; j <= nptBinsGen; ++j) {
      Double_t c = hResp->GetBinContent(i, j);
      rowTotal += c;
      // Check if truth bin overlaps with reco bin
      Double_t tlo = hResp->GetYaxis()->GetBinLowEdge(j);
      Double_t thi = hResp->GetYaxis()->GetBinUpEdge(j);
      if (tlo >= lo - 0.01 && thi <= hi + 0.01) {
        diagContent += c;
      }
    }
    Double_t diagFrac = (rowTotal > 0) ? diagContent / rowTotal : 0;
    out(Form("  [%3.0f,%3.0f]: total=%.6g, diagonal=%.6g, diagFrac=%.4f\n",
        lo, hi, rowTotal, diagContent, diagFrac));
  }

  // ============================================================
  // 7. Check if unweighted histograms are being used for JJ MC
  // ============================================================
  out("\n=== JJ MC weight check ===\n");
  if (hMCRecoWeighted && hMCRecoRaw) {
    Double_t unwInt = hMCRecoRaw->Integral();
    Double_t wInt = hMCRecoWeighted->Integral();
    out(Form("  Unweighted h_jet_pt integral: %.6g\n", unwInt));
    out(Form("  Weighted h_jet_pt integral:   %.6g\n", wInt));
    out(Form("  Ratio (weighted/unweighted):  %.6f\n", (unwInt > 0) ? wInt / unwInt : 0));
    if (TMath::Abs(wInt - unwInt) / TMath::Max(unwInt, 1.0) < 0.01) {
      out("  [WARNING] Weighted and unweighted are nearly identical — weighting may not be applied!\n");
    } else {
      out("  [OK] Weighted and unweighted differ significantly — JJ weights are applied.\n");
    }
  } else {
    out("  Weighted histograms not found. If this is JJ MC, the response matrix is WRONG.\n");
  }

  // ============================================================
  // 8. KEY CHECK: Are we loading UNWEIGHTED histograms for JJ MC?
  // ============================================================
  out("\n=== CRITICAL: Which histograms is LoadUnfoldingInputs loading? ===\n");
  out("  The macro loads:\n");
  out("    h_jet_pt      (unweighted) for MC reco\n");
  out("    h_jet_pt_part (unweighted) for MC truth\n");
  out("    h2_..._matchedgeo_mcdetaconstraint (unweighted) for response\n");
  out("  For JJ MC (LHC25a2b), these should be WEIGHTED versions!\n");

  if (hMCRecoWeighted) {
    Double_t ratio5 = 0, ratio50 = 0;
    Int_t bin5 = hMCRecoRaw->GetXaxis()->FindBin(5.5);
    Int_t bin50 = hMCRecoRaw->GetXaxis()->FindBin(50.5);
    if (hMCRecoRaw->GetBinContent(bin5) > 0)
      ratio5 = hMCRecoWeighted->GetBinContent(bin5) / hMCRecoRaw->GetBinContent(bin5);
    if (hMCRecoRaw->GetBinContent(bin50) > 0)
      ratio50 = hMCRecoWeighted->GetBinContent(bin50) / hMCRecoRaw->GetBinContent(bin50);
    out(Form("  weight/unweight ratio at pT~5.5 GeV:  %.6f\n", ratio5));
    out(Form("  weight/unweight ratio at pT~50.5 GeV: %.6f\n", ratio50));
    if (TMath::Abs(ratio5 - ratio50) > 0.01) {
      out("  [CRITICAL] pT-hat weighting changes shape significantly!\n");
      out("  Using unweighted histograms for JJ MC will give WRONG response matrix.\n");
    }
  }

  log.close();
  std::cerr << "\n[Info] Diagnostic log saved to DiagnoseDataMCMismatch_log.txt\n";
}
