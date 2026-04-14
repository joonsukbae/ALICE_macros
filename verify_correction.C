// Verify: Did TrackTuner actually change sigma(1/pT)?
void verify_correction() {
  // Baseline: DCA correction only, NO q/pT correction
  TFile* fBase = TFile::Open("../../jets/AnalysisResults/599366_AnalysisResults.root");
  // New test: Your custom CCDB q/pT correction
  TFile* fNew = TFile::Open("/Users/js/Downloads/TestMyFilePath_AnalysisResults.root");
  // Data for reference
  TFile* fData = TFile::Open("../../jets/AnalysisResults/493670_AnalysisResults.root");

  if (!fBase || !fNew || !fData) { printf("Cannot open files\n"); return; }

  TDirectory* dBase = (TDirectory*)fBase->Get("track-efficiency");
  TDirectory* dNew = (TDirectory*)fNew->Get("track-efficiency");
  TDirectory* dData = (TDirectory*)fData->Get("track-efficiency");

  // sigma(1/pT) histograms
  TH2* h2Base_s1pt = (TH2*)dBase->Get("h2_track_pt_track_sigma1overpt");
  TH2* h2New_s1pt = (TH2*)dNew->Get("h2_track_pt_track_sigma1overpt");
  TH2* h2Data_s1pt = (TH2*)dData->Get("h2_track_pt_track_sigma1overpt");

  // sigma(pT)/pT histograms
  TH2* h2Base_spt = (TH2*)dBase->Get("h2_track_pt_track_sigmapt");
  TH2* h2New_spt = (TH2*)dNew->Get("h2_track_pt_track_sigmapt");
  TH2* h2Data_spt = (TH2*)dData->Get("h2_track_pt_track_sigmapt");

  // High pT versions
  TH2* h2Base_s1pt_h = (TH2*)dBase->Get("h2_track_pt_high_track_sigma1overpt");
  TH2* h2New_s1pt_h = (TH2*)dNew->Get("h2_track_pt_high_track_sigma1overpt");
  TH2* h2Data_s1pt_h = (TH2*)dData->Get("h2_track_pt_high_track_sigma1overpt");

  TH2* h2Base_spt_h = (TH2*)dBase->Get("h2_track_pt_high_track_sigmapt");
  TH2* h2New_spt_h = (TH2*)dNew->Get("h2_track_pt_high_track_sigmapt");
  TH2* h2Data_spt_h = (TH2*)dData->Get("h2_track_pt_high_track_sigmapt");

  printf("=============================================================================\n");
  printf("  CRITICAL CHECK: Did TrackTuner actually modify the covariance matrix?\n");
  printf("=============================================================================\n\n");

  printf("File entries:\n");
  printf("  Baseline (no q/pT corr): %.0f entries\n", h2Base_s1pt ? h2Base_s1pt->GetEntries() : 0);
  printf("  NewTest (your CCDB):     %.0f entries\n", h2New_s1pt ? h2New_s1pt->GetEntries() : 0);
  printf("  Data:                    %.0f entries\n", h2Data_s1pt ? h2Data_s1pt->GetEntries() : 0);

  printf("\n=== sigma(1/pT) comparison (LOW pT) ===\n");
  printf("  pT | Baseline sig(1/pT) | NewTest sig(1/pT) | Data sig(1/pT) | New/Base | CHANGED?\n");
  printf("-----|---------------------|-------------------|----------------|----------|----------\n");

  vector<int> lowBins = {10, 30, 50, 70, 100};
  for (int bin : lowBins) {
    TH1* pBase = h2Base_s1pt->ProjectionY(Form("pb1_%d", bin), bin, bin);
    TH1* pNew = h2New_s1pt->ProjectionY(Form("pn1_%d", bin), bin, bin);
    TH1* pData = h2Data_s1pt->ProjectionY(Form("pd1_%d", bin), bin, bin);

    double pt = h2Base_s1pt->GetXaxis()->GetBinCenter(bin);
    double mBase = pBase->GetMean();
    double mNew = pNew->GetMean();
    double mData = pData->GetMean();
    double ratio = (mBase > 0) ? mNew / mBase : 0;
    const char* changed = (fabs(ratio - 1.0) > 0.01) ? "YES" : "NO";

    printf("%4.0f | %19.6f | %17.6f | %14.6f | %8.4f | %s\n",
           pt, mBase, mNew, mData, ratio, changed);

    delete pBase; delete pNew; delete pData;
  }

  printf("\n=== sigma(1/pT) comparison (HIGH pT) ===\n");
  printf("  pT | Baseline sig(1/pT) | NewTest sig(1/pT) | Data sig(1/pT) | New/Base | CHANGED?\n");
  printf("-----|---------------------|-------------------|----------------|----------|----------\n");

  vector<int> highBins = {10, 30, 50, 70, 90};
  for (int bin : highBins) {
    if (!h2Base_s1pt_h || !h2New_s1pt_h) continue;

    TH1* pBase = h2Base_s1pt_h->ProjectionY(Form("pbh1_%d", bin), bin, bin);
    TH1* pNew = h2New_s1pt_h->ProjectionY(Form("pnh1_%d", bin), bin, bin);
    TH1* pData = h2Data_s1pt_h->ProjectionY(Form("pdh1_%d", bin), bin, bin);

    double pt = h2Base_s1pt_h->GetXaxis()->GetBinCenter(bin);
    double mBase = pBase->GetMean();
    double mNew = pNew->GetMean();
    double mData = pData->GetMean();
    double ratio = (mBase > 0) ? mNew / mBase : 0;
    const char* changed = (fabs(ratio - 1.0) > 0.01) ? "YES" : "NO";

    printf("%4.0f | %19.6f | %17.6f | %14.6f | %8.4f | %s\n",
           pt, mBase, mNew, mData, ratio, changed);

    delete pBase; delete pNew; delete pData;
  }

  printf("\n=== sigma(pT)/pT comparison (LOW pT) ===\n");
  printf("  pT | Baseline sig(pT)/pT | NewTest sig(pT)/pT | Data sig(pT)/pT | New/Base | CHANGED?\n");
  printf("-----|----------------------|--------------------|-----------------|----------|----------\n");

  for (int bin : lowBins) {
    TH1* pBase = h2Base_spt->ProjectionY(Form("pbs_%d", bin), bin, bin);
    TH1* pNew = h2New_spt->ProjectionY(Form("pns_%d", bin), bin, bin);
    TH1* pData = h2Data_spt->ProjectionY(Form("pds_%d", bin), bin, bin);

    double pt = h2Base_spt->GetXaxis()->GetBinCenter(bin);
    double mBase = pBase->GetMean();
    double mNew = pNew->GetMean();
    double mData = pData->GetMean();
    double ratio = (mBase > 0) ? mNew / mBase : 0;
    const char* changed = (fabs(ratio - 1.0) > 0.01) ? "YES" : "NO";

    printf("%4.0f | %20.6f | %18.6f | %15.6f | %8.4f | %s\n",
           pt, mBase, mNew, mData, ratio, changed);

    delete pBase; delete pNew; delete pData;
  }

  printf("\n=== sigma(pT)/pT comparison (HIGH pT) ===\n");
  printf("  pT | Baseline sig(pT)/pT | NewTest sig(pT)/pT | Data sig(pT)/pT | New/Base | CHANGED?\n");
  printf("-----|----------------------|--------------------|-----------------|----------|----------\n");

  for (int bin : highBins) {
    if (!h2Base_spt_h || !h2New_spt_h) continue;

    TH1* pBase = h2Base_spt_h->ProjectionY(Form("pbsh_%d", bin), bin, bin);
    TH1* pNew = h2New_spt_h->ProjectionY(Form("pnsh_%d", bin), bin, bin);
    TH1* pData = h2Data_spt_h->ProjectionY(Form("pdsh_%d", bin), bin, bin);

    double pt = h2Base_spt_h->GetXaxis()->GetBinCenter(bin);
    double mBase = pBase->GetMean();
    double mNew = pNew->GetMean();
    double mData = pData->GetMean();
    double ratio = (mBase > 0) ? mNew / mBase : 0;
    const char* changed = (fabs(ratio - 1.0) > 0.01) ? "YES" : "NO";

    printf("%4.0f | %20.6f | %18.6f | %15.6f | %8.4f | %s\n",
           pt, mBase, mNew, mData, ratio, changed);

    delete pBase; delete pNew; delete pData;
  }

  printf("\n=============================================================================\n");
  printf("  VERDICT\n");
  printf("=============================================================================\n");
  printf("If New/Base ~ 1.0 for sigma(1/pT): TrackTuner is NOT applying correction!\n");
  printf("If New/Base > 1.0 for sigma(1/pT) but sigma(pT)/pT unchanged: Different issue!\n");

  fBase->Close();
  fNew->Close();
  fData->Close();
}
