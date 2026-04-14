// Direct inspection of CCDB graphs
void inspect_ccdb_direct() {
  const char* normalFile = "drawJetSpectraCharged/TrackTunerQA/Sigma1overPt_Data_LHC22o_pass7_globalTracks_MC_LHC25a2b.root";
  const char* extremeFile = "drawJetSpectraCharged/TrackTunerQA/Sigma1overPt_Data_LHC22o_pass7_globalTracks_MC_LHC25a2b_EXTREME_x5.root";

  TFile* fNormal = TFile::Open(normalFile);

  if (!fNormal) { printf("Cannot open normal file: %s\n", normalFile); return; }

  printf("=== Contents of NORMAL file ===\n");
  fNormal->ls();

  // Get the TList (CCDB format)
  TList* ccdbList = (TList*)fNormal->Get("ccdb_object");
  if (!ccdbList) {
    printf("ccdb_object TList not found!\n");
    return;
  }

  printf("\n=== Contents of ccdb_object TList ===\n");
  ccdbList->ls();

  TGraphErrors* grMC = (TGraphErrors*)ccdbList->FindObject("sigmaVsPtMc");
  TGraphErrors* grData = (TGraphErrors*)ccdbList->FindObject("sigmaVsPtData");

  if (!grMC) printf("sigmaVsPtMc NOT FOUND!\n");
  if (!grData) printf("sigmaVsPtData NOT FOUND!\n");

  if (!grMC || !grData) {
    printf("\nAll objects in TList:\n");
    TIter next(ccdbList);
    TObject* obj;
    while ((obj = next())) {
      printf("  %s (%s)\n", obj->GetName(), obj->ClassName());
    }
    return;
  }

  printf("\n=== NORMAL FILE: sigma(1/pT) values ===\n");
  printf("N points: MC=%d, Data=%d\n\n", grMC->GetN(), grData->GetN());
  printf("    pT (GeV/c) | sigma(1/pT)_MC | sigma(1/pT)_Data | Data/MC ratio\n");
  printf("   ------------|----------------|------------------|---------------\n");

  for (int i = 0; i < grMC->GetN(); i++) {
    double xMC, yMC, xData, yData;
    grMC->GetPoint(i, xMC, yMC);
    grData->GetPoint(i, xData, yData);
    double ratio = (yMC > 0) ? yData / yMC : 0;
    printf("   %11.2f | %14.6f | %16.6f | %13.4f\n", xMC, yMC, yData, ratio);
  }

  // Check extreme file
  TFile* fExtreme = TFile::Open(extremeFile);
  if (fExtreme) {
    TList* ccdbListExt = (TList*)fExtreme->Get("ccdb_object");
    if (ccdbListExt) {
      TGraphErrors* grMC_ext = (TGraphErrors*)ccdbListExt->FindObject("sigmaVsPtMc");
      TGraphErrors* grData_ext = (TGraphErrors*)ccdbListExt->FindObject("sigmaVsPtData");

      if (grMC_ext && grData_ext) {
        printf("\n=== EXTREME x5 FILE: sigma(1/pT) values ===\n");
        printf("N points: MC=%d, Data=%d\n\n", grMC_ext->GetN(), grData_ext->GetN());
        printf("    pT (GeV/c) | sigma(1/pT)_MC | sigma(1/pT)_Data | Data/MC ratio\n");
        printf("   ------------|----------------|------------------|---------------\n");

        for (int i = 0; i < grMC_ext->GetN(); i++) {
          double xMC, yMC, xData, yData;
          grMC_ext->GetPoint(i, xMC, yMC);
          grData_ext->GetPoint(i, xData, yData);
          double ratio = (yMC > 0) ? yData / yMC : 0;
          printf("   %11.2f | %14.6f | %16.6f | %13.4f\n", xMC, yMC, yData, ratio);
        }
      }
    }
    fExtreme->Close();
  }

  fNormal->Close();
}
