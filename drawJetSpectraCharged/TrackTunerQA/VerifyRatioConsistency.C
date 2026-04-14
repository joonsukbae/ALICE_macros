// Quick verification: compare Data/MC ratios for sigma(1/pT) vs sigma(pT)/pT
// If filling is consistent, ratios should be identical

#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include <iostream>

void CheckHistogramPair(TH2* h2Data_1overpt, TH2* h2Data_ptoverpt,
                        TH2* h2MC_1overpt, TH2* h2MC_ptoverpt,
                        const char* label) {
  std::cout << "\n========================================" << std::endl;
  std::cout << "=== " << label << " ===" << std::endl;
  std::cout << "========================================" << std::endl;

  std::cout << "\n=== Histogram Info ===" << std::endl;
  std::cout << "sigma1overpt - Data entries: " << h2Data_1overpt->GetEntries()
            << ", MC entries: " << h2MC_1overpt->GetEntries() << std::endl;
  std::cout << "sigmapt      - Data entries: " << h2Data_ptoverpt->GetEntries()
            << ", MC entries: " << h2MC_ptoverpt->GetEntries() << std::endl;

  std::cout << "\n=== pT axis ranges ===" << std::endl;
  std::cout << "sigma1overpt: pT [" << h2Data_1overpt->GetXaxis()->GetXmin() << ", "
            << h2Data_1overpt->GetXaxis()->GetXmax() << "] with "
            << h2Data_1overpt->GetNbinsX() << " bins" << std::endl;
  std::cout << "sigmapt:      pT [" << h2Data_ptoverpt->GetXaxis()->GetXmin() << ", "
            << h2Data_ptoverpt->GetXaxis()->GetXmax() << "] with "
            << h2Data_ptoverpt->GetNbinsX() << " bins" << std::endl;

  std::cout << "\n=== Y-axis (sigma) ranges ===" << std::endl;
  std::cout << "sigma1overpt: [" << h2Data_1overpt->GetYaxis()->GetXmin() << ", "
            << h2Data_1overpt->GetYaxis()->GetXmax() << "] with "
            << h2Data_1overpt->GetNbinsY() << " bins" << std::endl;
  std::cout << "sigmapt:      [" << h2Data_ptoverpt->GetYaxis()->GetXmin() << ", "
            << h2Data_ptoverpt->GetYaxis()->GetXmax() << "] with "
            << h2Data_ptoverpt->GetNbinsY() << " bins" << std::endl;

  // Check for overflow in Y-axis
  std::cout << "\n=== Overflow check (Y-axis) ===" << std::endl;
  double overflow_1overpt_data = 0, overflow_ptoverpt_data = 0;
  double overflow_1overpt_mc = 0, overflow_ptoverpt_mc = 0;
  for (int i = 1; i <= h2Data_1overpt->GetNbinsX(); i++) {
    overflow_1overpt_data += h2Data_1overpt->GetBinContent(i, h2Data_1overpt->GetNbinsY()+1);
    overflow_ptoverpt_data += h2Data_ptoverpt->GetBinContent(i, h2Data_ptoverpt->GetNbinsY()+1);
    overflow_1overpt_mc += h2MC_1overpt->GetBinContent(i, h2MC_1overpt->GetNbinsY()+1);
    overflow_ptoverpt_mc += h2MC_ptoverpt->GetBinContent(i, h2MC_ptoverpt->GetNbinsY()+1);
  }
  std::cout << "sigma1overpt overflow - Data: " << overflow_1overpt_data << ", MC: " << overflow_1overpt_mc << std::endl;
  std::cout << "sigmapt overflow      - Data: " << overflow_ptoverpt_data << ", MC: " << overflow_ptoverpt_mc << std::endl;

  // Compare ratios at specific pT values
  std::cout << "\n=== Ratio comparison at specific pT bins ===" << std::endl;
  std::cout << "pT (GeV/c) | sigma(1/pT) Data/MC | sigma(pT)/pT Data/MC | Difference" << std::endl;
  std::cout << "-----------|---------------------|----------------------|------------" << std::endl;

  double maxDiff = 0;
  for (int i = 1; i <= h2Data_1overpt->GetNbinsX(); i++) {
    double pt = h2Data_1overpt->GetXaxis()->GetBinCenter(i);

    TH1* proj_1overpt_data = h2Data_1overpt->ProjectionY(Form("tmp1_%d", i), i, i);
    TH1* proj_1overpt_mc = h2MC_1overpt->ProjectionY(Form("tmp2_%d", i), i, i);
    TH1* proj_ptoverpt_data = h2Data_ptoverpt->ProjectionY(Form("tmp3_%d", i), i, i);
    TH1* proj_ptoverpt_mc = h2MC_ptoverpt->ProjectionY(Form("tmp4_%d", i), i, i);

    if (proj_1overpt_data->GetEntries() < 100 || proj_1overpt_mc->GetEntries() < 100) {
      delete proj_1overpt_data; delete proj_1overpt_mc;
      delete proj_ptoverpt_data; delete proj_ptoverpt_mc;
      continue;
    }

    double mean_1overpt_data = proj_1overpt_data->GetMean();
    double mean_1overpt_mc = proj_1overpt_mc->GetMean();
    double mean_ptoverpt_data = proj_ptoverpt_data->GetMean();
    double mean_ptoverpt_mc = proj_ptoverpt_mc->GetMean();

    double ratio_1overpt = (mean_1overpt_mc > 0) ? mean_1overpt_data / mean_1overpt_mc : 0;
    double ratio_ptoverpt = (mean_ptoverpt_mc > 0) ? mean_ptoverpt_data / mean_ptoverpt_mc : 0;
    double diff = (ratio_1overpt > 0) ? (ratio_ptoverpt - ratio_1overpt) / ratio_1overpt * 100 : 0;

    if (fabs(diff) > maxDiff) maxDiff = fabs(diff);

    std::cout << Form("%10.2f | %19.4f | %20.4f | %9.2f%%", pt, ratio_1overpt, ratio_ptoverpt, diff) << std::endl;

    // Also verify the mathematical relationship for first few bins
    if (i <= 3) {
      double expected_ptoverpt_data = mean_1overpt_data * pt;
      double expected_ptoverpt_mc = mean_1overpt_mc * pt;
      std::cout << Form("           | Expected sigmapt = sigma1overpt * pT:") << std::endl;
      std::cout << Form("           |   Data: expected=%.6f, actual=%.6f, diff=%.2f%%",
                        expected_ptoverpt_data, mean_ptoverpt_data,
                        (mean_ptoverpt_data - expected_ptoverpt_data) / expected_ptoverpt_data * 100) << std::endl;
      std::cout << Form("           |   MC:   expected=%.6f, actual=%.6f, diff=%.2f%%",
                        expected_ptoverpt_mc, mean_ptoverpt_mc,
                        (mean_ptoverpt_mc - expected_ptoverpt_mc) / expected_ptoverpt_mc * 100) << std::endl;
    }

    delete proj_1overpt_data; delete proj_1overpt_mc;
    delete proj_ptoverpt_data; delete proj_ptoverpt_mc;
  }

  std::cout << "\nMax ratio difference for " << label << ": " << maxDiff << "%" << std::endl;
}

void VerifyRatioConsistency(
  const char* dataFile = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/610951_AnalysisResults.root",
  const char* mcFile = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/611764_AnalysisResults.root",
  const char* directory = "track-efficiency"
) {
  TFile* fData = TFile::Open(dataFile, "READ");
  TFile* fMC = TFile::Open(mcFile, "READ");

  if (!fData || !fMC) {
    std::cerr << "Cannot open files" << std::endl;
    return;
  }

  TDirectory* dirData = (TDirectory*)fData->Get(directory);
  TDirectory* dirMC = (TDirectory*)fMC->Get(directory);

  // Get histograms - low pT range
  TH2* h2Data_1overpt = (TH2*)dirData->Get("h2_track_pt_track_sigma1overpt");
  TH2* h2Data_ptoverpt = (TH2*)dirData->Get("h2_track_pt_track_sigmapt");
  TH2* h2MC_1overpt = (TH2*)dirMC->Get("h2_track_pt_track_sigma1overpt");
  TH2* h2MC_ptoverpt = (TH2*)dirMC->Get("h2_track_pt_track_sigmapt");

  // Get histograms - high pT range
  TH2* h2Data_1overpt_high = (TH2*)dirData->Get("h2_track_pt_high_track_sigma1overpt");
  TH2* h2Data_ptoverpt_high = (TH2*)dirData->Get("h2_track_pt_high_track_sigmapt");
  TH2* h2MC_1overpt_high = (TH2*)dirMC->Get("h2_track_pt_high_track_sigma1overpt");
  TH2* h2MC_ptoverpt_high = (TH2*)dirMC->Get("h2_track_pt_high_track_sigmapt");

  // Check low pT histograms
  CheckHistogramPair(h2Data_1overpt, h2Data_ptoverpt, h2MC_1overpt, h2MC_ptoverpt, "LOW pT (0-10 GeV)");

  // Check high pT histograms
  CheckHistogramPair(h2Data_1overpt_high, h2Data_ptoverpt_high, h2MC_1overpt_high, h2MC_ptoverpt_high, "HIGH pT (10-100 GeV)");

  std::cout << "\n========================================" << std::endl;
  std::cout << "=== CONCLUSION ===" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "If ratios differ significantly, the mathematical relationship" << std::endl;
  std::cout << "sigma(pT)/pT = sigma(1/pT) * pT is NOT being honored." << std::endl;
  std::cout << "This means Data/MC ratios for the two quantities will differ." << std::endl;
  std::cout << "\nFor TrackTuner, use sigma(1/pT) as that's what the code expects." << std::endl;

  fData->Close();
  fMC->Close();
}
