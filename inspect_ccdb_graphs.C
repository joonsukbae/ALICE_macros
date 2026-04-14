// Inspect the actual values in your CCDB q/pT graphs
void inspect_ccdb_graphs() {
  // Your CCDB file
  TFile* f = TFile::Open("../../jets/AnalysisResults/drawJetSpectraCharged/TrackTunerQA/qOverPtGraphs.root");

  if (!f || f->IsZombie()) {
    printf("Cannot open file!\n");
    return;
  }

  TGraphErrors* grMC = (TGraphErrors*)f->Get("sigmaVsPtMc");
  TGraphErrors* grData = (TGraphErrors*)f->Get("sigmaVsPtData");

  if (!grMC || !grData) {
    printf("Cannot find graphs!\n");
    printf("Available objects in file:\n");
    f->ls();
    return;
  }

  printf("=============================================================================\n");
  printf("  CCDB Graph Inspection: sigma(1/pT) vs pT\n");
  printf("=============================================================================\n\n");

  printf("Graph entries: MC=%d, Data=%d\n\n", grMC->GetN(), grData->GetN());

  printf("    pT (GeV/c) | sigma(1/pT)_MC | sigma(1/pT)_Data | Data/MC ratio\n");
  printf("   ------------|----------------|------------------|---------------\n");

  // Sample at specific pT values
  vector<double> ptValues = {0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 50.0, 70.0, 100.0};

  for (double pt : ptValues) {
    double valMC = grMC->Eval(pt);
    double valData = grData->Eval(pt);
    double ratio = (valMC > 0) ? valData / valMC : 0;

    printf("   %11.1f | %14.6f | %16.6f | %13.4f\n", pt, valMC, valData, ratio);
  }

  printf("\n=== Raw graph points (MC) ===\n");
  printf("  Point |    pT    | sigma(1/pT)_MC\n");
  printf("  ------|----------|---------------\n");
  for (int i = 0; i < grMC->GetN(); i++) {
    double x, y;
    grMC->GetPoint(i, x, y);
    printf("  %5d | %8.2f | %14.6f\n", i, x, y);
  }

  printf("\n=== Raw graph points (Data) ===\n");
  printf("  Point |    pT    | sigma(1/pT)_Data\n");
  printf("  ------|----------|------------------\n");
  for (int i = 0; i < grData->GetN(); i++) {
    double x, y;
    grData->GetPoint(i, x, y);
    printf("  %5d | %8.2f | %16.6f\n", i, x, y);
  }

  printf("\n=============================================================================\n");
  printf("  INTERPRETATION\n");
  printf("=============================================================================\n");
  printf("If Data/MC ratio is close to 1.0 everywhere, no correction will be applied.\n");
  printf("If ratio > 1, MC resolution will be DEGRADED (smeared) to match Data.\n");
  printf("If ratio < 1, something is wrong (Data should have worse resolution).\n");
  printf("\nCheck if the pT range of your graphs covers the analysis range!\n");

  f->Close();
}
