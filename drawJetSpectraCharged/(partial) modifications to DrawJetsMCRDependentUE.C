++ Update File: DrawJetsMCRDependentUE.C
@@
 inline TString GetOutputDirUE() {
   // Separate output directory for UE-subtracted R-dependent plots
   static const TString outputDir = "../plots/RDependentComparison/LHC24f3c_UEsub";
   return outputDir;
 }

// ======================================================================
// Run 2 UE-subtracted reference (HEPData Figure 3)
// ======================================================================

struct Run2UEXsecPoint {
  double pt;
  double sigma;
  double statErr;
  double sysErr;
};

// Hard-coded Run 2 UE-subtracted cross sections (Figure 3)
inline std::vector<Run2UEXsecPoint> GetRun2UEXsecData(double R) {
  std::vector<Run2UEXsecPoint> points;

  // R = 0.4
  if (TMath::Abs(R - 0.4) < 0.01) {
    points = {
      {5.5, 0.73271, 0.00014903, 0.062585}, {6.5, 0.39868, 9.8507e-05, 0.035592},
      {7.5, 0.2329, 6.8456e-05, 0.021378},  {8.5, 0.14445, 4.91e-05, 0.013835},
      {9.5, 0.094007, 3.4181e-05, 0.0094001},{11.0, 0.053479, 2.0988e-05, 0.0053987},
      {13.0, 0.026927, 1.2121e-05, 0.0027346},{15.0, 0.014953, 8.0655e-06, 0.0014948},
      {17.0, 0.0089322, 5.5281e-06, 0.00089048},{19.0, 0.0056478, 3.8943e-06, 0.00057549},
      {22.5, 0.0028768, 2.1461e-06, 0.00028901},{27.5, 0.0011753, 1.079e-06, 0.00011904},
      {35.0, 0.00041953, 4.7569e-07, 4.268e-05},{45.0, 0.00012699, 1.7316e-07, 1.3475e-05},
      {55.0, 4.8626e-05, 7.7014e-08, 5.2298e-06},{65.0, 2.1715e-05, 3.8654e-08, 2.3392e-06},
      {77.5, 9.3914e-06, 1.825e-08, 1.049e-06},{92.5, 3.8434e-06, 7.9487e-09, 4.2323e-07},
      {120.0, 1.1284e-06, 2.4285e-09, 1.2455e-07}
    };
  }
  // R = 0.7
  else if (TMath::Abs(R - 0.7) < 0.01) {
    points = {
      {5.5, 1.0815, 0.00035038, 0.082193}, {6.5, 0.62221, 0.00025633, 0.055258},
      {7.5, 0.39121, 0.00018964, 0.033611}, {8.5, 0.25658, 0.00014485, 0.022893},
      {9.5, 0.17185, 0.00010275, 0.01618},  {11.0, 0.099064, 6.12e-05, 0.0099809},
      {13.0, 0.048885, 3.5214e-05, 0.004959},{15.0, 0.026002, 2.3401e-05, 0.0027029},
      {17.0, 0.015078, 1.5457e-05, 0.0016587},{19.0, 0.0094971, 1.023e-05, 0.0010626},
      {22.5, 0.004699, 5.4377e-06, 0.00052095},{27.5, 0.0017396, 2.4159e-06, 0.00020834},
      {35.0, 0.00057631, 1.0265e-06, 6.9532e-05},{45.0, 0.00016294, 3.6086e-07, 2.0523e-05},
      {55.0, 6.1234e-05, 1.607e-07, 8.1226e-06},{65.0, 2.7226e-05, 8.117e-08, 3.7395e-06},
      {77.5, 1.1687e-05, 3.824e-08, 1.5777e-06},{92.5, 4.7354e-06, 1.6527e-08, 6.3909e-07},
      {120.0, 1.3679e-06, 4.9716e-09, 1.7564e-07}
    };
  }

  return points;
}

// Create Run2 UE-sub TGraphErrors (stat) and TGraphAsymmErrors (sys)
// with x-errors matched to Run 3 binning
inline std::pair<TGraphErrors*, TGraphAsymmErrors*>
CreateRun2UEGraphs(double R) {
  std::vector<Run2UEXsecPoint> points = GetRun2UEXsecData(R);
  if (points.empty()) return {nullptr, nullptr};

  const Double_t* ptbinGen = GetPtbinGen();
  int nBinsGen = 25; // ptbinGen has 26 edges

  int nPoints = points.size();
  TGraphErrors* grStat = new TGraphErrors(nPoints);
  TGraphAsymmErrors* grSys = new TGraphAsymmErrors(nPoints);

  grStat->SetName(Form("Run2UE_Stat_R%.1f", R));
  grSys->SetName(Form("Run2UE_Sys_R%.1f", R));

  grStat->SetLineColor(kBlack);
  grStat->SetMarkerColor(kBlack);
  grStat->SetMarkerStyle(20);
  grStat->SetLineWidth(2);

  grSys->SetFillColorAlpha(kGray + 1, 0.5);
  grSys->SetFillStyle(1001);
  grSys->SetLineWidth(0);
  grSys->SetLineColor(0);

  for (int i = 0; i < nPoints; ++i) {
    const auto& pt = points[i];

    // Match Run 3 bin center
    int binIdx = -1;
    for (int j = 0; j < nBinsGen; ++j) {
      double binCenter = 0.5 * (ptbinGen[j] + ptbinGen[j+1]);
      if (TMath::Abs(pt.pt - binCenter) < 0.11) {
        binIdx = j;
        break;
      }
    }

    double exLow = 0.0, exHigh = 0.0;
    if (binIdx >= 0 && binIdx < nBinsGen) {
      double binLow = ptbinGen[binIdx];
      double binHigh = ptbinGen[binIdx+1];
      exLow  = pt.pt - binLow;
      exHigh = binHigh - pt.pt;
    } else {
      // Fallback: symmetric errors from neighbor spacing
      if (i > 0 && i < nPoints - 1) {
        exLow  = 0.5 * (pt.pt - points[i-1].pt);
        exHigh = 0.5 * (points[i+1].pt - pt.pt);
      } else if (i == 0 && nPoints > 1) {
        exLow = exHigh = 0.5 * (points[i+1].pt - pt.pt);
      } else if (i == nPoints - 1 && nPoints > 1) {
        exLow = exHigh = 0.5 * (pt.pt - points[i-1].pt);
      }
    }

    grStat->SetPoint(i, pt.pt, pt.sigma);
    grStat->SetPointError(i, 0.5 * (exLow + exHigh), pt.statErr);

    grSys->SetPoint(i, pt.pt, pt.sigma);
    grSys->SetPointError(i, exLow, exHigh, pt.sysErr, pt.sysErr);
  }

  return {grStat, grSys};
}

// Helper: create ratio histogram (Run3 UE-sub / Run2 UE-sub)
inline TH1* CreateRatioTH1vsTGraphUE(TH1* hNum, TGraphErrors* grDenom,
                                     const char* name) {
  if (!hNum || !grDenom) return nullptr;

  TH1* hRatio = (TH1*)hNum->Clone(name);
  hRatio->Reset();
  hRatio->SetDirectory(0);
  hRatio->SetStats(0);

  int nPoints = grDenom->GetN();
  double* xDen = grDenom->GetX();
  double* yDen = grDenom->GetY();
  double* eDen = grDenom->GetEY();

  for (int i = 1; i <= hNum->GetNbinsX(); ++i) {
    double xLow = hNum->GetBinLowEdge(i);
    double xUp  = xLow + hNum->GetBinWidth(i);
    double yDenVal = 0.0;
    double eDenVal = 0.0;
    bool found = false;

    for (int j = 0; j < nPoints; ++j) {
      if (xDen[j] >= xLow && xDen[j] < xUp) {
        yDenVal = yDen[j];
        eDenVal = eDen[j];
        found = true;
        break;
      }
    }

    if (!found || yDenVal <= 0) continue;

    double yNum = hNum->GetBinContent(i);
    double eNum = hNum->GetBinError(i);
    if (yNum <= 0) continue;

    double ratio = yNum / yDenVal;
    double relNum = eNum / yNum;
    double relDen = eDenVal / yDenVal;
    double eRatio = ratio * TMath::Sqrt(relNum * relNum + relDen * relDen);

    hRatio->SetBinContent(i, ratio);
    hRatio->SetBinError(i, eRatio);
  }

  return hRatio;
}

// ======================================================================
// Draw R-dependent UE-subtracted cross sections (Run3 vs Run2 UE-sub)
// ======================================================================


