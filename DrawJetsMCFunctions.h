///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Draw macro for jet QA         //////////
////////// author: Joonsuk Bae           ////////// 
////////// E-mail: jbae@cern.ch          //////////
////////// Last Modified: 19 July 2024   //////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

////////////////
/// Set Bins ///
////////////////
#include <complex>
#include <vector>

void ALICEfigureLegend(const char *FigureLabel, double x11, double y11, double x12, double y12, double x21, double y21, double x22, double y22, double textsize=0.043);

void DrawJetsMC() {

  DrawHistos(fileNames, histNames, ColorPallete);

}

// define fns
Double_t Nevents(const char *fileName, const char *eventDir,
              const char *eventObj, Int_t ifMCP = 0) {
  auto file = TFile::Open(fileName, "open");
  auto Nevents = (TH1D *)file->Get(Form("%s/%s", eventDir, eventObj));

  Double_t nevents = 0;
  if (ifMCP == 0) {
    nevents = Nevents->GetBinContent(Nevents->FindBin(1.5));
  } else if (ifMCP == 1) {
    // auto NevtsMCP = (TH1D *) file->Get("bc-selection-task/hLumiTVX"); 
    // nevents = 1000 * NevtsMCP->Integral(1, NevtsMCP->GetNbinsX());
    // cout << "hLumiTVX: " << nevents << " / mb" << endl;
    nevents = Nevents->GetBinContent(Nevents->FindBin(0.5));
  }
  return nevents;
}
void PrintHistogramErrors(TH1 *hist) {
  int nBins = hist->GetNbinsX();
  for (int i = 1; i <= nBins; ++i) { // bin id starts from 1
    double binContent = hist->GetBinContent(i);
    double binError = hist->GetBinError(i);
    // std::cout << "Bin " << i << ": Content = " << binContent
    //           << ", Error = " << binError << std::endl;
  }
}

TH1 *DrawRatioTH1(TH1 *hNum, TH1 *hDenom) {
    TH1 *hRatio = (TH1 *)hNum->Clone("hRatio");
    hRatio->Reset(); 

    for (int i = 1; i <= hNum->GetNbinsX(); i++) {
        double xCenter = hNum->GetXaxis()->GetBinCenter(i); 
        double yNum = hNum->GetBinContent(i);
        double yNumErr = hNum->GetBinError(i); 

        int binDenom = hDenom->GetXaxis()->FindBin(xCenter);
        double yDenom = hDenom->GetBinContent(binDenom);

        if (yDenom != 0) {
            double ratio = yNum / yDenom;
            double ratioErr = ratio * (yNumErr / yNum); 

            hRatio->SetBinContent(i, ratio);
            hRatio->SetBinError(i, ratioErr);
        } else {
            hRatio->SetBinContent(i, 0);
            hRatio->SetBinError(i, 0); 
        }
    }

    return hRatio;
}

TH1 *DrawRatio(const char *ratioName, TH1 *refHist, TH1 *testHist,
               TString AxisTitleX, TString AxisTitleY, Color_t colorID,
               Double_t Ymin, Double_t Ymax, Double_t MarkerSize = 1) {
  TH1 *ratioHist = (TH1 *)refHist->Clone(ratioName);
  ratioHist = DrawRatioTH1(testHist, refHist);
  hset(*ratioHist, AxisTitleX, AxisTitleY, 1.2, 1.0, 0.07, 0.07, 0.01, 0.01,
       0.07, 0.07, 510, 505);
  // ratioHist->Divide(testHist, ratioHist, 1., 1., "B");
  ratioHist->GetYaxis()->SetRangeUser(Ymin, Ymax);
  // ratioHist->GetYaxis()->SetLimits(Ymin, Ymax);
  ratioHist->SetMarkerColor(colorID);
  ratioHist->SetLineColor(colorID);
  ratioHist->SetMarkerSize(MarkerSize);
  ratioHist->Draw("esame");

  return ratioHist;
}

// TH1 *DrawRatioTH1(TH1 *hNum, TH1 *hDenom, bool useSmallerBins = true) {
//     TH1 *hBase, *hCompare;

//     if (useSmallerBins) {
//         if (hNum->GetNbinsX() > hDenom->GetNbinsX()) {
//             hBase = hDenom;
//             hCompare = hNum;
//         } else {
//             hBase = hNum;
//             hCompare = hDenom;
//         }
//     } else {
//         hBase = hNum;
//         hCompare = hDenom;
//     }

//     TH1 *hRatio = (TH1 *)hBase->Clone("hRatio");
//     hRatio->Reset();

//     for (int i = 1; i <= hBase->GetNbinsX(); i++) {
//         double xCenter = hBase->GetXaxis()->GetBinCenter(i);
//         double yBase = hBase->GetBinContent(i);
//         double yBaseErr = hBase->GetBinError(i);

//         int binComp = hCompare->GetXaxis()->FindBin(xCenter);
//         double yComp = hCompare->GetBinContent(binComp);
//         double yCompErr = hCompare->GetBinError(binComp);

//         if (yComp != 0) {
//             double ratio = yBase / yComp;
//             double ratioErr = ratio * sqrt(pow(yBaseErr / yBase, 2) + pow(yCompErr / yComp, 2));

//             hRatio->SetBinContent(i, ratio);
//             hRatio->SetBinError(i, ratioErr);
//         } else {
//             hRatio->SetBinContent(i, 0);
//             hRatio->SetBinError(i, 0);
//         }
//     }

//     return hRatio;
// }
TH1 *MultiplyTH1(TH1 *hDef, TH1 *hComp) {
    TH1 *hBase, *hCompare;
    if (hDef->GetNbinsX() > hComp->GetNbinsX()) {
        hBase = hComp;
        hCompare = hDef;
    } else {
        hBase = hDef;
        hCompare = hComp;
    }

    TH1 *hMultiply = (TH1 *)hBase->Clone("hMultiply");
    hMultiply->Reset();
    
    for (int i = 1; i <= hBase->GetNbinsX(); i++) {
        double xCenter = hBase->GetXaxis()->GetBinCenter(i);
        double yBase = hBase->GetBinContent(i);
        double yBaseErr = hBase->GetBinError(i);
        
        int binComp = hCompare->GetXaxis()->FindBin(xCenter);
        double yComp = hCompare->GetBinContent(binComp);
        double yCompErr = hCompare->GetBinError(binComp);

        double product = yBase * yComp;
        double productErr = product * sqrt(pow(yBaseErr / yBase, 2) + pow(yCompErr / yComp, 2));

        hMultiply->SetBinContent(i, product);
        hMultiply->SetBinError(i, productErr);
    }

    return hMultiply;
}
// // Run 2
// TH1D* GraphToHistogram(TGraphErrors* graph) {
//   int nPoints = graph->GetN();
//   double* xValues = graph->GetX();
//   double* yValues = graph->GetY();
  
//   // Define histogram with same number of bins as the number of points in the graph
//   double xMin = graph->GetXaxis()->GetXmin();
//   double xMax = graph->GetXaxis()->GetXmax();
//   TH1D* hist = new TH1D("hist", "Histogram from Graph", nptBins, ptbin);

//   // Fill histogram with the values from the graph
//   for (int i = 0; i < nPoints; ++i) {
//     hist->SetBinContent(i+1, yValues[i]);
//   }
  
//   return hist;
// }
TH1D* GraphToHistogram(TGraphErrors* graph) {
    int nPoints = graph->GetN();
    double* xValues = graph->GetX();
    double* yValues = graph->GetY();
    double* yErrors = graph->GetEY();

    TH1D* hist = new TH1D("hist", "Histogram from Graph", nptBins, ptbin);

    for (int i = 0; i < nPoints; ++i) {
        int bin = hist->FindBin(xValues[i]);
        hist->SetBinContent(bin, yValues[i]);
        hist->SetBinError(bin, yErrors[i]);
    }

    return hist;
}
TGraphErrors *DrawRatioTGraph(TGraphErrors * gr, TH1 *hist) {
  int Npoint = gr->GetN();
  TGraphErrors *gr_ratio = new TGraphErrors(Npoint);
  double x, y1, ey1, y2, ex, ratio, ratio_err;
  for (auto i=0; i < Npoint; i++) {
    x = gr->GetX()[i];
    y1 = gr->GetY()[i];
    ey1 = gr->GetErrorY(i);
  
    auto histXbin = hist->GetXaxis()->FindBin(x);
    y2 = hist->GetBinContent(histXbin);
    ratio = y1 / y2;
    ratio_err = ey1 / y2;
    gr_ratio->SetPoint(i, x, ratio);
    gr_ratio->SetPointError(i, hist->GetXaxis()->GetBinWidth(histXbin)/2, ratio_err);
    
    // cout << "DrawRatioTgraph(npoint, x, y, y_err): (" << i << ", " << x << ", " << ratio << ", " << ratio_err << ")" << endl;
  }  

  return gr_ratio;
}
TGraphErrors *DrawTGraph(TGraphErrors * gr, TH1 *hist) {
  int Npoint = gr->GetN();
  TGraphErrors *gr_error = new TGraphErrors(Npoint);
  double x, ex1, y1, ey1;
  for (auto i=0; i < Npoint; i++) {
    x = gr->GetX()[i];
    y1 = gr->GetY()[i];
    ex1 = hist->GetXaxis()->GetBinWidth(i+5);
    cout << "(Npoint, x, y): " << i << ", " << x << ", " << y1 << endl;
    ey1 = gr->GetErrorY(i);

    gr_error->SetPoint(i, x, y1);
    gr_error->SetPointError(i, ex1, ey1);

  }

  return gr_error;
}

std::pair<double, double> getYAxisRange(TH1* histogram, Double_t Nevts, double Xmin, double Xmax, double UpMarginFactor = 5, double DownMarginFactor = 15) {
    if (!histogram) {
        std::cerr << "Invalid histogram pointer!" << std::endl;
        return {0, 1}; // 기본 범위 반환
    }

    double yMax = histogram->GetBinContent(histogram->FindBin(Xmin)) / Nevts;
    double yMin = histogram->GetBinContent(histogram->FindBin(Xmax)) / Nevts;
    yMin = yMin / DownMarginFactor;
    yMax = yMax * UpMarginFactor;

    return {yMin, yMax};
}

TGraphErrors *Run2Data(TH1 *hist) {
  // auto Run2Datafile = TFile::Open("~/cernbox/workspace/O2Physics/jets/Run2_YongzhenHou/HEPData-ins2026265-v1-Figure_3.root", "OPEN");// w/ UE sub.
  auto Run2Datafile = TFile::Open("~/cernbox/workspace/O2Physics/jets/Run2_YongzhenHou/HEPData-ins2026265-v1-Figure_A1.root", "OPEN");// w/o UE sub.
  // auto Run2DataGraph = (TGraphErrors *)Run2Datafile->Get("Figure 3/Graph1D_y3");
  auto Run2DataGraph = (TGraphErrors *)Run2Datafile->Get("Figure A1/Graph1D_y3");

  // TGraphErrors* Run2DataHist = DrawTGraph(Run2DataGraph, hist);
  int Npoint = Run2DataGraph->GetN();
  TGraphErrors *gr_error = new TGraphErrors(Npoint);
  double x, ex1, y1, ey1;
  for (auto i=0; i < Npoint; i++) {
    x = Run2DataGraph->GetX()[i];
    y1 = Run2DataGraph->GetY()[i];
    ex1 = hist->GetXaxis()->GetBinWidth(hist->GetXaxis()->FindBin(x))/2;
    ey1 = Run2DataGraph->GetErrorY(i);

    gr_error->SetPoint(i, x, y1);
    gr_error->SetPointError(i, ex1, ey1);

  }
  return gr_error;
}

#include <vector>
#include <TGraphErrors.h>

std::vector<TGraphErrors*> Run2Data_manual() {
    const int Npoint = 19;

    double pt[Npoint] = {5.5, 6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17, 19, 22.5, 27.5, 35, 45, 55, 65, 77.5, 92.5, 120};
    
#ifdef UESUB
    // New data for UESUB case
    double sigma[Npoint] = {0.73271, 0.39868, 0.2329, 0.14445, 0.094007, 0.053479, 0.026927, 0.014953, 0.0089322, 0.0056478,
                            0.0028768, 0.0011753, 0.00041953, 0.00012699, 4.8626e-05, 2.1715e-05, 9.3914e-06, 3.8434e-06, 1.1284e-06};

    double stat_err[Npoint] = {0.00014903, 9.8507e-05, 6.8456e-05, 4.91e-05, 3.4181e-05, 2.0988e-05, 1.2121e-05, 8.0655e-06, 5.5281e-06, 3.8943e-06,
                               2.1461e-06, 1.079e-06, 4.7569e-07, 1.7316e-07, 7.7014e-08, 3.8654e-08, 1.825e-08, 7.9487e-09, 2.4285e-09};

    double syst_err[Npoint] = {0.062585, 0.035592, 0.021378, 0.013835, 0.0094001, 0.0053987, 0.0027346, 0.0014948, 0.00089048, 0.00057549,
                               0.00028901, 0.00011904, 4.268e-05, 1.3475e-05, 5.2298e-06, 2.3392e-06, 1.049e-06, 4.2323e-07, 1.2455e-07};
#else
    // Original data for non-UESUB case
    double sigma[Npoint] = {1.0972, 0.58418, 0.33224, 0.20001, 0.12614, 0.069825, 0.034225, 0.018498, 0.010764, 0.0066421,
                            0.0033177, 0.0013313, 0.00046794, 0.00013985, 5.2994e-05, 2.3472e-05, 1.0087e-05, 4.1091e-06, 1.2028e-06};

    double stat_err[Npoint] = {0.00016396, 0.00010887, 7.5637e-05, 5.3822e-05, 3.7416e-05, 2.2965e-05, 1.3182e-05, 8.6563e-06, 5.8873e-06, 4.1364e-06,
                               2.2742e-06, 1.136e-06, 4.9698e-07, 1.7976e-07, 7.9569e-08, 3.9801e-08, 1.8747e-08, 8.1518e-09, 2.4881e-09};

    double syst_err[Npoint] = {0.074335, 0.042393, 0.024799, 0.015407, 0.010138, 0.0057313, 0.0029103, 0.001585, 0.00094719, 0.00061878,
                               0.00031195, 0.00012884, 4.6138e-05, 1.4546e-05, 5.6228e-06, 2.5033e-06, 1.1183e-06, 4.4943e-07, 1.3178e-07};
#endif

    TGraphErrors *gr_stat_error = new TGraphErrors(Npoint);
    TGraphErrors *gr_syst_error = new TGraphErrors(Npoint);

    for (int i = 0; i < Npoint; i++) {
        double ex1 = (pt[i+1] - pt[i]) / 2.0;

        gr_stat_error->SetPoint(i, pt[i], sigma[i]);
        gr_stat_error->SetPointError(i, ex1, stat_err[i]);

        gr_syst_error->SetPoint(i, pt[i], sigma[i]);
        gr_syst_error->SetPointError(i, ex1, syst_err[i]);
    }

    std::vector<TGraphErrors*> graphs = {gr_stat_error, gr_syst_error};
    return graphs;
}


TH1 *Run2MCgen() {
  auto Run2MCfile = TFile::Open("~/cernbox/workspace/O2Physics/jets/Run2_YongzhenHou/FinalsixGeneratorPP13TeV.root", "OPEN");
  // TH3 *h3Run2MCGraph = (TH3D *) Run2MCfile->Get("GenJetptetaphi04");
  TH3 *h3Run2MCGraph = (TH3D *) Run2MCfile->Get("GenJetCorrptetaphi04");

  TH1* Run2MCGraph = h3Run2MCGraph->Project3D("z");
  // TH1* Run2MCGraph = (TH1 *) Run2MCfile->Get("totJetpt04");
  Run2MCGraph->GetXaxis()->SetRangeUser(1,200);
  Run2MCGraph = Run2MCGraph->Rebin(nptBinsGen, "run2MCgraphRebined", ptbinGen);
  Run2MCGraph->Scale(1., "width");

  return Run2MCGraph;
}
TH1 *GetPYTHIA1360() {
  // auto File = TFile::Open("../../jets/mc/AnalysisResults/Systematics/13over136GenPythia/1360Charged/AnalysisResults.root", "OPEN");
  auto File = TFile::Open("../../jets/mc/AnalysisResults/Systematics/13over136GenPythia/1360Charged/AnalysisResults_1360new_highstat.root", "OPEN");
  TH1* hPythia = (TH1 *) File->Get("hJetPt");
  TH1* hNevt = (TH1 *) File->Get("hnevent");
  double nevt = hNevt->GetEntries();
  hPythia->GetXaxis()->SetRangeUser(1,200);
  // hPythia = hPythia->Rebin(nptBinsGen, "hPythia", ptbinGen);
  hPythia->Scale(77.77/nevt, "width");

  return hPythia;
}
TH1 *GetPYTHIA1300() {
  auto File = TFile::Open("../../jets/mc/AnalysisResults/Systematics/13over136GenPythia/1300Charged/AnalysisResults.root", "OPEN");
  TH1* hPythia = (TH1 *) File->Get("hJetPt");
  TH1* hNevt = (TH1 *) File->Get("hnevent");
  double nevt = hNevt->GetEntries();
  hPythia->GetXaxis()->SetRangeUser(1,200);
  // hPythia = hPythia->Rebin(nptBinsGen, "hPythia", ptbinGen);
  hPythia->Scale(1./nevt, "width");

  return hPythia;
}
TH1 *GetHERWIG1360() {
  // auto HerwigFile = TFile::Open("../../jets/mc/AnalysisResults/Systematics/AnalysisResults_herwig1360.root", "OPEN");
  // auto HerwigFile = TFile::Open("../../jets/mc/AnalysisResults/Systematics/Herwig/AnalysisResults1360.root", "OPEN"); // twice scale
  auto HerwigFile = TFile::Open("../../jets/mc/AnalysisResults/Systematics/Herwig/HERWIGMB_1360.root", "OPEN");
  // TH1* hHerwig = (TH1 *) HerwigFile->Get("hJetPt");
  TH1* hHerwig = (TH1 *) HerwigFile->Get("hJetPtX");
  TH1* hNevt = (TH1 *) HerwigFile->Get("hnevent");
  double nevt = hNevt->GetBinContent(1);
  hHerwig->GetXaxis()->SetRangeUser(1,200);
  // hHerwig = hHerwig->Rebin(nptBinsGen, "hHerwig", ptbinGen);
  // hHerwig->Scale(1./200000, "width");
  hHerwig->Scale(1./nevt, "width");

  return hHerwig;
}

void OutStatsTXT(TH1 *hist, const char *txtname) {
    TString dirPath = Form("plots/systematic/%s", histNames[0].Data());
    
    if (mkdir("plots", 0777) == -1 && errno != EEXIST) {
        std::cerr << "Error creating directory 'plots'\n";
    }
    if (mkdir("plots/systematic", 0777) == -1 && errno != EEXIST) {
        std::cerr << "Error creating directory 'plots/systematic'\n";
    }
    if (mkdir(dirPath.Data(), 0777) == -1 && errno != EEXIST) {
        std::cerr << "Error creating directory '" << dirPath << "'\n";
    }
    
    if (DRAWPLOTS) {
    std::cout << "OutStatsTXT: " << txtname << std::endl;
      std::ofstream outFile(Form("%s/SystErr_%s.txt", dirPath.Data(), txtname));

      for (int i = 1; i <= hist->GetNbinsX(); i++) {
          double binCenter = hist->GetXaxis()->GetBinCenter(i);
          double binContent = std::abs((hist->GetBinContent(i) - 1) * 100);
          TString outputLine = Form("Bin: %d, pT: %.2f, Error: %.5f %%\n", i, binCenter, binContent);

          outFile << outputLine;
      }

      outFile.close();
    }
}

void calculateQuadratureSum(const std::vector<TH1*>& hSysts, TH1* hSystResult) {
    int nBins = hSysts[0]->GetNbinsX();

    // Step 1: Calculate quadrature sum
    for (int i = 1; i <= nBins; ++i) {
        double sumOfSquares = 0.0;

        for (const auto& histo : hSysts) {
            double value = histo->GetBinContent(i) - 1;
            sumOfSquares += value * value;
        }

        double quadratureSum = std::sqrt(sumOfSquares);
        hSystResult->SetBinContent(i, quadratureSum + 1);  
    }

    // Step 2: Smoothing using 3-bin moving average
    TH1* hSystSmoothed = (TH1*) hSysts[0]->Clone("hSystSmoothed");  
    
    for (int i = 1; i <= nBins; ++i) {
        double smoothedValue = 0.0;
        
        if (i == 1) {
            smoothedValue = (hSystResult->GetBinContent(i) + hSystResult->GetBinContent(i + 1)) / 2.0;
        } else if (i == nBins) {
            smoothedValue = (hSystResult->GetBinContent(i) + hSystResult->GetBinContent(i - 1)) / 2.0;
        } else {
            smoothedValue = (hSystResult->GetBinContent(i - 1) +
                             hSystResult->GetBinContent(i) +
                             hSystResult->GetBinContent(i + 1)) / 3.0;
        }

        hSystSmoothed->SetBinContent(i, smoothedValue);
    }

    for (int i = 1; i <= nBins; ++i) {
        hSystResult->SetBinContent(i, hSystSmoothed->GetBinContent(i));
    }

    delete hSystSmoothed;
}
void DrawMultipleSources(std::vector<TH1*>& hSysts, TH1* hSystResult) {
    if (hSysts.size() < 7) {
        std::cerr << "Error: There should be at least 7 histograms!" << std::endl;
        return;
    }

    TCanvas* cSources = new TCanvas("cSources", "Systematic Uncertainty Sources", 800, 600);
    setpad(cSources, 0.02, 0.15, 0.12, 0.03);
    cSources->cd();

    TLegend* legend = new TLegend(0.622807,0.648696,0.799499,0.93913,NULL,"brNDC");
    legend->SetTextSize(0.043);
    legend->SetBorderSize(0);
    legend->SetFillColorAlpha(0, 0);

    bool firstHist = true;

    std::vector<TH1*> newHSysts;
    std::vector<TString> SourceNames = {"Tracking efficiency", "Track #it{p}_{T} resolution", "Unfolding", "Normalization", "Secondary particles", "Total uncertainty"};

    // Step 1: Create the Unfolding source by combining the histograms [2], [3], and [4] using Quadrature Sum
    TH1* hUnfolding = (TH1*)hSysts[2]->Clone("hUnfolding");
    hUnfolding->Reset();

    // Loop through each bin and calculate the quadrature sum of the unfolding sources
    for (int i = 1; i <= hUnfolding->GetNbinsX(); ++i) {
        double sumOfSquares = 0.0;

        // Loop through sources 2, 3, and 4 and subtract 1 from each value (since they are ratios)
        for (int j = 2; j <= 4; ++j) {
            double value = (hSysts[j]->GetBinContent(i) - 1.0) * 100;  // Subtract 1 and convert to percentage
            sumOfSquares += value * value;  // Add the square of the normalized value
        }

        double quadratureSum = std::sqrt(sumOfSquares);  // Calculate quadrature sum
        hUnfolding->SetBinContent(i, quadratureSum);  // Set the value for the unfolding source in %
    }

    // Step 2: Convert all sources to percentage by subtracting 1 and multiplying by 100
    for (int i = 0; i < hSysts.size(); ++i) {
        if (i == 2 || i == 3 || i == 4) continue;  // Skip the sources that were quadrature summed

        for (int bin = 1; bin <= hSysts[i]->GetNbinsX(); ++bin) {
            double newValue = (hSysts[i]->GetBinContent(bin) - 1.0) * 100;  // Subtract 1 and convert to percentage
            hSysts[i]->SetBinContent(bin, newValue);
        }
    }

    // Step 3: Convert hSystResult (Total Uncertainty) to percentage
    auto hSystTotal = (TH1 *) hSystResult->Clone("hSystTotal");
    for (int bin = 1; bin <= hSystTotal->GetNbinsX(); ++bin) {
        double newValue = (hSystTotal->GetBinContent(bin) - 1.0) * 100;  // Subtract 1 and convert to percentage
        hSystTotal->SetBinContent(bin, newValue);
    }

    // Step 4: Add the sources to the new histogram vector
    newHSysts.push_back(hSysts[0]);  // track efficiency
    newHSysts.push_back(hSysts[1]);  // track pT resolution
    newHSysts.push_back(hUnfolding); // Unfolding (quadrature sum of [2], [3], and [4])
    newHSysts.push_back(hSysts[5]);  // Normalization
    newHSysts.push_back(hSysts[6]);  // Secondary Particles
    newHSysts.push_back(hSystTotal); // Total Uncertainty (smoothing된 total error)

    // Step 5: Draw each histogram in the newHSysts vector
    for (size_t i = 0; i < newHSysts.size(); ++i) {
        TH1* hist = newHSysts[i];
        TString sourceName = SourceNames[i];

        hist->SetLineWidth(5);
        hist->SetLineStyle(i);
        hist->SetMarkerStyle(0);  // No markers
        hist->SetFillStyle(0);    // No fill
        hist->SetLineColor(ColorPallete[i]);  // Use the color palette

        // Custom histogram options
        hoptset(*hist, 0, ColorPallete[i], PlotPtMin, PlotPtMax, 0, 60.0, 0.6, 1, 2, 25);  // Adjust y-axis for percentage
        hset(*hist, JetPtDataFinalTitleX, "Rel. Uncertainties (%)", 1.3, 1.1, 0.05, 0.05, 0.01, 0.01,
             0.05, 0.06, 510, 505);

        if (firstHist) {
            hist->Draw("HIST");  // Draw first histogram normally
            firstHist = false;
        } else {
            hist->Draw("HIST SAME");  // Overlay the rest
        }

        // Add entry to the legend
        legend->AddEntry(hist, sourceName, "l");
    }

    // Draw the legend
    legend->Draw();
    ALICEfigureLegend("ALICE Preliminary", 0.110276,0.68,0.309524,0.96,0.109023,0.429565,0.309524,0.638261, 0.04);

    // Save the canvas if the DRAWPLOTS flag is set
    if (SYSTUNFOLD) {
        cSources->SaveAs(Form("%s/SystematicUncertaintySources.pdf", MakeDirName.Data()));
    }
}

TH1 *DrawUnfoldHerwig(TH1 * hPythiaInvY, TH1 *hPythiaTrue, TH1 *hPythiaReco, TH2 *h2Correlate, TH1 *hData, double Nevt) {
  std::cout << "Starting DrawUnfoldHerwig function." << std::endl;

  // Retrieve HERWIG/Pythia ratio and check its validity
  TH1 *hHerwig1360 = GetHERWIG1360();
  hHerwig1360->Scale(1./77.77);
  // TH1 *hPythia1360 = GetPYTHIA1360();
  TH1 *hHerwigPythiaRat = DrawRatioTH1(hPythiaInvY, hHerwig1360);

  std::cout << "Herwig/Pythia Ratio" << std::endl;
  for (int i = 1; i <= hHerwigPythiaRat->GetNbinsX(); i++) {
    std::cout << "(bin, pT, value): (" << i << ", " << hHerwigPythiaRat->GetBinCenter(i) << ", " << hHerwigPythiaRat->GetBinContent(i) << ")" << std::endl;
  }

  // Herwig의 pT 범위를 얻습니다.
  double minPT = hHerwigPythiaRat->GetXaxis()->GetXmin();
  double maxPT = hHerwigPythiaRat->GetXaxis()->GetXmax();

  // pT 범위에 맞춰 Pythia 히스토그램을 트림합니다.
  int binMin = hPythiaTrue->GetXaxis()->FindBin(minPT);
  int binMax = hPythiaTrue->GetXaxis()->FindBin(maxPT);

  TH1 *hPythiaTrueTrimmed = (TH1 *)hPythiaTrue->Clone("hPythiaTrueTrimmed");
  TH1 *hPythiaRecoTrimmed = (TH1 *)hPythiaReco->Clone("hPythiaRecoTrimmed");

  for (int i = 1; i < hPythiaTrueTrimmed->GetNbinsX(); i++) {
    if (i < binMin || i > binMax) {
      hPythiaTrueTrimmed->SetBinContent(i, 0);
      hPythiaTrueTrimmed->SetBinError(i, 0);
      hPythiaRecoTrimmed->SetBinContent(i, 0);
      hPythiaRecoTrimmed->SetBinError(i, 0);
    }
  }

  // Apply HERWIG/Pythia ratio
  TH1 *hTure = MultiplyTH1(hPythiaTrueTrimmed, hHerwigPythiaRat);
  // TH1 *hTure = (TH1 *) hHerwigPythiaRat->Clone(Form("hist_%i", ++n));
  TH1 *hReco = MultiplyTH1(hPythiaRecoTrimmed, hHerwigPythiaRat);
  // TH1 *hReco = (TH1 *) hPythiaRecoTrimmed->Clone(Form("hist_%i", ++n));
  if (!hTure || !hReco) {
    std::cerr << "Error: Cloning of hPythiaTrue or hPythiaReco failed." << std::endl;
    return nullptr;
  }
  std::cout << "hPythiaTrue and hPythiaReco successfully cloned and trimmed." << std::endl;

  // Clone the correlation histogram and scale it
  TH2* h2CorrelateHerwig = (TH2*)h2Correlate->Clone(Form("hist_%i", ++n));
  if (!h2CorrelateHerwig) {
    std::cerr << "Error: Cloning of h2Correlate failed." << std::endl;
    return nullptr;
  }
  std::cout << "h2Correlate successfully cloned." << std::endl;

  for (int y = 1; y <= h2CorrelateHerwig->GetNbinsY(); y++) {
    double fscale = hHerwigPythiaRat->GetBinContent(y);
    for (int x = 1; x <= h2CorrelateHerwig->GetNbinsX(); x++) {
      double content = h2CorrelateHerwig->GetBinContent(x, y);
      double error = h2CorrelateHerwig->GetBinError(x, y);
      h2CorrelateHerwig->SetBinContent(x, y, content * fscale);
      h2CorrelateHerwig->SetBinError(x, y, error * fscale);
    }
  }
  std::cout << "h2CorrelateHerwig successfully scaled." << std::endl;

  // Generate projections
  auto hMatchTrue = h2CorrelateHerwig->ProjectionY("hMatchTrue", binMin, binMax);
  auto hMatchReco = h2CorrelateHerwig->ProjectionX("hMatchReco", binMin, binMax);
  if (!hMatchTrue || !hMatchReco) {
    std::cerr << "Error: Projections hMatchTrue or hMatchReco failed." << std::endl;
    return nullptr;
  }
  std::cout << "Projections hMatchTrue and hMatchReco successfully created." << std::endl;

  // Create and manipulate response matrix
  TH2F *h2Res = (TH2F *)h2CorrelateHerwig->Clone(Form("hist_%i", ++n));
  TH1F *hFake = (TH1F *)hReco->Clone(Form("hist_%i", ++n));
  hFake->Add(hMatchReco, -1);
  TH1F *hMiss = (TH1F *)hTure->Clone(Form("hist_%i", ++n));
  hMiss->Add(hMatchTrue, -1);

  if (!h2Res || !hFake || !hMiss) {
    std::cerr << "Error: Cloning or manipulating histograms failed." << std::endl;
    return nullptr;
  }
  std::cout << "Response histograms successfully created and manipulated." << std::endl;

  // Build the response object
  RooUnfoldResponse *hResponse = new RooUnfoldResponse(hReco, hTure);
  if (!hResponse) {
    std::cerr << "Error: RooUnfoldResponse object creation failed." << std::endl;
    return nullptr;
  }
  std::cout << "RooUnfoldResponse object successfully created." << std::endl;

  for (auto i = binMin; i <= binMax; i++) {
    for (auto j = 1; j <= h2Res->GetNbinsY(); j++) {
      Double_t bincenx = h2Res->GetXaxis()->GetBinCenter(i);
      Double_t binceny = h2Res->GetYaxis()->GetBinCenter(j);
      Double_t bincont = h2Res->GetBinContent(i, j);
      for (auto k = 0;  k<bincont; k++){
      hResponse->Fill(bincenx, binceny);
      }
    }
  }
  std::cout << "RooUnfoldResponse filled with data from h2Res." << std::endl;


  for (auto i = binMin; i <= binMax; i++) {
    Double_t bincenx = hMiss->GetXaxis()->GetBinCenter(i);
    Double_t bincont = hMiss->GetBinContent(i);
    hResponse->Miss(bincenx, bincont);
  }
  std::cout << "RooUnfoldResponse filled with miss data." << std::endl;

  for (auto i = binMin; i <= binMax; i++) {
    Double_t bincenx = hFake->GetXaxis()->GetBinCenter(i);
    Double_t bincont = hFake->GetBinContent(i);
    hResponse->Fake(bincenx, bincont);
  }
  std::cout << "RooUnfoldResponse filled with fake data." << std::endl;

  // Perform the unfolding
  RooUnfoldBayes unfoldHerwig(hResponse, hData, 4);
  auto hUnfoldedDataHerwig = (TH1 *)unfoldHerwig.Hreco();
  if (!hUnfoldedDataHerwig) {
    std::cerr << "Error: Unfolding process failed." << std::endl;
    return nullptr;
  }
  std::cout << "Unfolding completed successfully." << std::endl;

  return hUnfoldedDataHerwig;
}
// TH1 *DrawUnfoldHerwig(TH1 *hPythiaTrue, TH1 *hPythiaReco, TH2 *h2Correlate, TH1 *hData, double Nevt) {
//   std::cout << "Starting DrawUnfoldHerwig function." << std::endl;

//   // Retrieve HERWIG/Pythia ratio and check its validity
//   TH1 *hHerwigPythiaRat = GetHERWIG1360();
//   TH1 *normPythia = (TH1 *)hPythiaTrue->Clone(Form("hist_%i", ++n));
//   normPythia->Scale(1./Nevt, "width");
//   DrawRatioTH1(hHerwigPythiaRat, normPythia);
//   // hHerwigPythiaRat->Divide(hHerwigPythiaRat, normPythia, 1, 1, "B");
//   cout << "Herwig/Pythia Ratio" <<endl;
//   for (int i=1; i<=hHerwigPythiaRat->GetNbinsX(); i++) {
//     cout << "(bin, pT, value): (" << i << ", " << hHerwigPythiaRat->GetBinCenter(i) << ", " << hHerwigPythiaRat->GetBinContent(i) << ")" << endl;
//   }

//   // Clone histograms and apply scaling
//   TH1 *hTure = MultiplyTH1(hTure, hHerwigPythiaRat);
//   TH1 *hReco = MultiplyTH1(hReco, hHerwigPythiaRat);
//   if (!hTure || !hReco) {
//     std::cerr << "Error: Cloning of hPythiaTrue or hPythiaReco failed." << std::endl;
//     return nullptr;
//   }
//   std::cout << "hPythiaTrue and hPythiaReco successfully cloned." << std::endl;

//   // hTure->Multiply(hTure, hHerwigPythiaRat, 1, 1, "B");
//   // hReco->Multiply(hReco, hHerwigPythiaRat, 1, 1, "B");

//   // Clone the correlation histogram and scale it
//   TH2* h2CorrelateHerwig = (TH2*)h2Correlate->Clone(Form("hist_%i", ++n));
//   if (!h2CorrelateHerwig) {
//     std::cerr << "Error: Cloning of h2Correlate failed." << std::endl;
//     return nullptr;
//   }
//   std::cout << "h2Correlate successfully cloned." << std::endl;

//   for (int y = 1; y <= h2CorrelateHerwig->GetNbinsY(); y++) {
//     double fscale = hHerwigPythiaRat->GetBinContent(y);
//     for (int x = 1; x <= h2CorrelateHerwig->GetNbinsX(); x++) {
//       double content = h2CorrelateHerwig->GetBinContent(x, y);
//       double error = h2CorrelateHerwig->GetBinError(x, y);
//       h2CorrelateHerwig->SetBinContent(x, y, content * fscale);
//       h2CorrelateHerwig->SetBinError(x, y, error * fscale);
//     }
//   }
//   std::cout << "h2CorrelateHerwig successfully scaled." << std::endl;

//   // Generate projections
//   auto hMatchTrue = h2CorrelateHerwig->ProjectionY("hMatchTrue", 1, h2CorrelateHerwig->GetNbinsX());
//   auto hMatchReco = h2CorrelateHerwig->ProjectionX("hMatchReco", 1, h2CorrelateHerwig->GetNbinsY());
//   if (!hMatchTrue || !hMatchReco) {
//     std::cerr << "Error: Projections hMatchTrue or hMatchReco failed." << std::endl;
//     return nullptr;
//   }
//   std::cout << "Projections hMatchTrue and hMatchReco successfully created." << std::endl;

//   // Create and manipulate response matrix
//   TH2F *h2Res = (TH2F *)h2CorrelateHerwig->Clone(Form("hist_%i", ++n));
//   TH1F *hFake = (TH1F *)hReco->Clone(Form("hist_%i", ++n));
//   hFake->Add(hMatchReco, -1);
//   TH1F *hMiss = (TH1F *)hTure->Clone(Form("hist_%i", ++n));
//   hMiss->Add(hMatchTrue, -1);

//   if (!h2Res || !hFake || !hMiss) {
//     std::cerr << "Error: Cloning or manipulating histograms failed." << std::endl;
//     return nullptr;
//   }
//   std::cout << "Response histograms successfully created and manipulated." << std::endl;

//   // Build the response object
//   RooUnfoldResponse *hResponse = new RooUnfoldResponse(hReco, hTure);
//   if (!hResponse) {
//     std::cerr << "Error: RooUnfoldResponse object creation failed." << std::endl;
//     return nullptr;
//   }
//   std::cout << "RooUnfoldResponse object successfully created." << std::endl;

//   for (auto i = 1; i <= h2Res->GetNbinsX(); i++) {
//     for (auto j = 1; j <= h2Res->GetNbinsY(); j++) {
//       Double_t bincenx = h2Res->GetXaxis()->GetBinCenter(i);
//       Double_t binceny = h2Res->GetYaxis()->GetBinCenter(j);
//       Double_t bincont = h2Res->GetBinContent(i, j);
//       hResponse->Fill(bincenx, binceny, bincont);
//     }
//   }
//   std::cout << "RooUnfoldResponse filled with data from h2Res." << std::endl;

//   for (auto i = 1; i <= hMiss->GetNbinsX(); i++) {
//     Double_t bincenx = hMiss->GetXaxis()->GetBinCenter(i);
//     Double_t bincont = hMiss->GetBinContent(i);
//     hResponse->Miss(bincenx, bincont);
//   }
//   std::cout << "RooUnfoldResponse filled with miss data." << std::endl;

//   for (auto i = 1; i <= hFake->GetNbinsX(); i++) {
//     Double_t bincenx = hFake->GetXaxis()->GetBinCenter(i);
//     Double_t bincont = hFake->GetBinContent(i);
//     hResponse->Fake(bincenx, bincont);
//   }
//   std::cout << "RooUnfoldResponse filled with fake data." << std::endl;

//   // Perform the unfolding
//   RooUnfoldBayes unfoldHerwig(hResponse, hData, 4);
//   auto hUnfoldedDataHerwig = (TH1 *)unfoldHerwig.Hreco();
//   if (!hUnfoldedDataHerwig) {
//     std::cerr << "Error: Unfolding process failed." << std::endl;
//     return nullptr;
//   }
//   std::cout << "Unfolding completed successfully." << std::endl;

//   return hUnfoldedDataHerwig;
// }
double CalculateLcurve(double residualNorm, double regularizationNorm) {
    return std::log(residualNorm) + std::log(regularizationNorm);
}

// double optimalK;
// void OptimizeRegularizationParameter(RooUnfoldResponse &responseMatrix, TH1 *dataJetPt, int maxK) {
//     std::vector<double> residualNorms;
//     std::vector<double> regularizationNorms;
//     std::vector<double> kValues;

//     // Loop over possible k-values to evaluate the L-curve
//     for (int k = 1; k <= maxK; ++k) {
//         RooUnfoldSvd unfoldSVD(&responseMatrix, dataJetPt, k);
//         TH1* unfolded = unfoldSVD.Hreco();

//         // Calculate residual norm (difference between measured and unfolded)
//         double residualNorm = 0.0;
//         double regularizationNorm = 0.0; // Initialize regularization norm

//         for (int i = 1; i <= dataJetPt->GetNbinsX(); ++i) {
//             double observed = dataJetPt->GetBinContent(i);
//             double expected = unfolded->GetBinContent(i);
//             double error = dataJetPt->GetBinError(i);

//             if (error != 0) {
//                 double residual = (observed - expected) / error;
//                 residualNorm += residual * residual;
//             }

//             // Assuming regularizationNorm is something like the sum of singular values or chi-squared approximation
//             regularizationNorm += std::abs(expected); // Simplified example, customize as needed
//         }

//         residualNorm = std::sqrt(residualNorm);

//         // Store the norms for L-curve analysis
//         residualNorms.push_back(residualNorm);
//         regularizationNorms.push_back(regularizationNorm);
//         kValues.push_back(k);

//         std::cout << "k = " << k << ": Residual Norm = " << residualNorm
//                   << ", Regularization Norm = " << regularizationNorm << std::endl;
//     }

//     // Analyze L-curve to find the optimal k
//     int optimalK = 1;
//     double minLcurveValue = CalculateLcurve(residualNorms[0], regularizationNorms[0]);
//     for (size_t i = 1; i < kValues.size(); ++i) {
//         double lcurveValue = CalculateLcurve(residualNorms[i], regularizationNorms[i]);
//         if (lcurveValue < minLcurveValue) {
//             minLcurveValue = lcurveValue;
//             optimalK = kValues[i];
//         }
//     }

//     std::cout << "Optimal k-value found: " << optimalK << std::endl;

//     // Plot L-curve
//     TCanvas *cLcurve = new TCanvas("cLcurve", "L-curve Analysis", 800, 600);
//     TGraph *lcurveGraph = new TGraph(kValues.size());
//     for (size_t i = 0; i < kValues.size(); ++i) {
//         lcurveGraph->SetPoint(i, std::log(residualNorms[i]), std::log(regularizationNorms[i]));
//     }
//     lcurveGraph->SetTitle("L-curve;log(Residual Norm);log(Regularization Norm)");
//     lcurveGraph->SetMarkerStyle(20);
//     lcurveGraph->Draw("ALP");

//     cLcurve->SaveAs("plots/AN_Charged-particle-jet-cross-section-in-pp-collisions-at-13.6-TeV/Figures/systematics/SVDUnfoldLcurve.pdf");

//     // Now perform the unfolding with the optimal k-value
//     RooUnfoldSvd unfoldSVDOptimal(&responseMatrix, dataJetPt, optimalK);
//     TH1* unfoldedOptimal = unfoldSVDOptimal.Hreco();

//     // Optionally, save the unfolded result
//     TCanvas *cUnfolded = new TCanvas("cUnfolded", "Unfolded Spectrum", 800, 600);
//     unfoldedOptimal->Draw();
//     cUnfolded->SaveAs("plots/AN_Charged-particle-jet-cross-section-in-pp-collisions-at-13.6-TeV/Figures/systematics/SVDUnfoldedOptimalK.pdf");
// }

TH1* DrawSecondaryContaimination(TH1* hCorrData) {
  std::vector<TH1*> hSystErrSecCon;
    Filipad2 *secPad = new Filipad2(++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
    secPad->Draw();
    TPad *secpad = secPad->GetPad(1);
    optFili(*secpad, 1, 1, 0, 1);
    TPad *secratpad = secPad->GetPad(2);
    optFili(*secratpad, 1, 1, 0, 0);

    TLegend *legSec = new TLegend(0.313397, 0.6, 0.578947, 0.843478, NULL, "brNDC");
    legSec->SetTextSize(0.05);
    legSec->SetBorderSize(0);
    legSec->SetFillColorAlpha(0, 0);

    secpad->cd();
    auto DefaultSec = (TH1 *) hCorrData->Clone(Form("hist_%i", nn));
    DefaultSec->Draw("pe");
    DefaultSec->GetXaxis()->SetRangeUser(PlotPtMin,PlotPtMax);
    legSec->AddEntry(DefaultSec, "Default", "lpe");

    TF1 *FitSec = new TF1("FitSec", "([2] + [3]*x) * pow(1 + x/([0]*[1]), -[1]) + ([6] + [7]*x) * pow(1 + x/([4]*[5]), -[5])", 5, PlotPtMax);
    FitSec->SetParameter(0, 0.7);
    FitSec->SetParameter(1, 7.1);
    FitSec->SetParameter(2, -210);
    FitSec->SetParameter(3, 63);
    FitSec->SetParameter(4, 0.7);
    FitSec->SetParameter(5, 7.1);
    FitSec->SetParameter(6, 210);
    FitSec->SetParameter(7, -62);
    DefaultSec->Fit(FitSec, "I");
    FitSec->SetLineColor(kBlack);
    FitSec->Draw("lsame");
    legSec->AddEntry(FitSec, "Fit Default", "lp");

    TF1 *FitSecUp = new TF1("FitSecUp", "([2] + [3]*1.01*x) * pow(1 + 1.01*x/([0]*[1]), -[1]) + ([6] + [7]*1.01*x) * pow(1 + 1.01*x/([4]*[5]), -[5])", 5, PlotPtMax);
    FitSecUp->SetParameters(FitSec->GetParameters());
    FitSecUp->SetLineColor(kBlue);
    FitSecUp->Draw("lsame");
    legSec->AddEntry(FitSecUp, "Fit +1%", "lp");

    TF1 *FitSecDown = new TF1("FitSecDown", "([2] + [3]*0.99*x) * pow(1 + 0.99*x/([0]*[1]), -[1]) + ([6] + [7]*0.99*x) * pow(1 + 0.99*x/([4]*[5]), -[5])", 5, PlotPtMax);
    FitSecDown->SetParameters(FitSec->GetParameters());
    FitSecDown->SetLineColor(kGreen);
    FitSecDown->Draw("lsame");
    legSec->AddEntry(FitSecDown, "Fit -1%", "lp");

    legSec->Draw();
    gStyle->SetOptFit(111);
    secpad->Update();

    TH1 *hist_default = new TH1D("hist_default", "Original Fit", nptBinsGen, ptbinGen);
    TH1 *hist_up = new TH1D("hist_up", "+1% Fit",  nptBinsGen, ptbinGen);
    TH1 *hist_down = new TH1D("hist_down", "-1% Fit",  nptBinsGen, ptbinGen);

    TH1 *hist_ratio_up = new TH1D("hist_ratio_up", "Ratio +1%",  nptBinsGen, ptbinGen);
    TH1 *hist_ratio_down = new TH1D("hist_ratio_down", "Ratio -1%",  nptBinsGen, ptbinGen);
    TH1 *hist_ratio_default = new TH1D("hist_ratio_default", "Ratio Default",  nptBinsGen, ptbinGen);

    for (int i = 1; i <= nptBinsGen; ++i) {
        double x = (ptbinGen[i-1] + ptbinGen[i]) / 2.0;

        double y_xsection = DefaultSec->GetBinContent(DefaultSec->FindBin(x));  
        double y_default = FitSec->Eval(x); 
        double y_up = FitSecUp->Eval(x);
        double y_down = FitSecDown->Eval(x);

        hist_default->SetBinContent(i, y_default);
        hist_up->SetBinContent(i, y_up);
        hist_down->SetBinContent(i, y_down);

        double ratio_up = y_up / y_default;
        double ratio_down = y_down / y_default;
        double ratio_default = y_default / y_xsection; 

        hist_ratio_up->SetBinContent(i, ratio_up);
        hist_ratio_down->SetBinContent(i, ratio_down);
        hist_ratio_default->SetBinContent(i, ratio_default);
    }

    secpad->cd();
    hist_default->SetMarkerStyle(24);
    hist_default->SetMarkerColor(kBlack);
    hset(*hist_default, JetPtDataFinalTitleX, "Comp. / Default", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
          0.1, 0.1, 510, 505);
    hist_default->Draw("P same");

    hist_up->SetMarkerStyle(24);
    hist_up->SetMarkerColor(kBlue);
    hist_up->Draw("P same");

    hist_down->SetMarkerStyle(24);
    hist_down->SetMarkerColor(kGreen);
    hist_down->Draw("P same");

    secratpad->cd();
    TLegend *legSecRat = new TLegend(0.586124,0.747826,0.861244,0.991304,NULL,"brNDC");
    legSecRat->SetTextSize(0.1);
    legSecRat->SetBorderSize(0);
    legSecRat->SetFillColorAlpha(0, 0);

    hist_ratio_default->SetMarkerStyle(25);
    hist_ratio_default->SetMarkerColor(kBlack);
    hset(*hist_ratio_default, JetPtDataFinalTitleX, "Comp. / Default", 1.2, 0.75, 0.1, 0.09, 0.01, 0.01, 0.1, 0.1, 510, 505);
    hoptset(*hist_ratio_default, 0, kBlack, PlotPtMin, PlotPtMax, 0.7, 1.3, 1, 1, 1, 25);
    legSecRat->AddEntry(hist_ratio_default, "Tsallis fit / Run 3", "p");
    hist_ratio_default->Draw("P same");

    hist_ratio_up->SetMarkerStyle(21);
    hist_ratio_up->SetMarkerColor(kBlue);
    hist_ratio_up->SetMarkerSize(0.5);
    hist_ratio_up->Draw("P same");

    hist_ratio_down->SetMarkerStyle(20);
    hist_ratio_down->SetMarkerColor(kGreen);
    hist_ratio_down->SetMarkerSize(0.5);
    hist_ratio_down->Draw("P same");

    legSecRat->Draw();
    // secpad->Update();

    if (SYSTUNFOLD) {
    secPad->C->Print(Form("%s/SystErrSecondaryContamination.pdf", MakeDirName.Data()));
    }

    return hist_ratio_up;
}

TH1 *GetTriggerEfficiency() {
  auto INELfile = TFile::Open("/Users/js/cernbox/workspace/O2Physics/jets/mc/AnalysisResults/LHC24f3/selMC/YieldCorrection/AnalysisResults_all.root", "read");
  auto hINELjet = (TH1 *)INELfile->Get("jet-finder-charged-qa/h_jet_pt_part");
  // auto hNevtINEL = (TH1 *)INELfile->Get("jet-finder-charged-qa/h_mccollisions");
  // double NevtINEL = hNevtINEL->GetBinContent(hNevtINEL->FindBin(1.5));
  // hINELjet->Scale(1./NevtINEL, "width");

  auto TRIGfile = TFile::Open("/Users/js/cernbox/workspace/O2Physics/jets/mc/AnalysisResults/LHC24f3/selMC/YieldCorrection/AnalysisResults_match.root", "read");
  auto hTRIGjet = (TH1 *)TRIGfile->Get("jet-finder-charged-qa/h_jet_pt_part");
  // auto hNevtTRIG = (TH1 *)TRIGfile->Get("jet-finder-charged-qa/h_mccollisions");
  // double NevtTRIG = hNevtTRIG->GetBinContent(hNevtTRIG->FindBin(2.5));
  // hTRIGjet->Scale(1./NevtTRIG, "width");

  TH1 *hTrigEff = (TH1 *) hTRIGjet->Clone(Form("hist_%i", ++n));
  hTrigEff->Divide(hTrigEff, hINELjet, 1, 1, "B");

  return hTrigEff;
}

TH1* ApplySystematicUncertainty(TH1* hist, TH1* errorHist) {
    TH1* histWithError = (TH1*)hist->Clone("histWithError");

    int nBins = hist->GetNbinsX();
    int nErrorBins = errorHist->GetNbinsX();

    for (int i = 1; i <= nBins; ++i) {
        double binContent = hist->GetBinContent(i);
        double statError = hist->GetBinError(i);
        
        if (i <= nErrorBins) {
            double sysUncertainty = errorHist->GetBinContent(i) - 1; 
            double sysError = binContent * sysUncertainty;

            double totalError = sqrt(pow(statError, 2) + pow(sysError, 2));
            histWithError->SetBinError(i, totalError);
        } else {
            histWithError->SetBinError(i, 0);
        }
    }

    return histWithError;
}
TH1 *Run3XSectionWoTrackTuner() {
  auto Run3file = TFile::Open("Run3_CrossSection_woTrackTuner.root", "read");
  auto hRun3Xsection = (TH1 *) Run3file->Get("Run3_CrossSection");
  return hRun3Xsection;
}
TH1 *Run3XSectionWTrackEff() {
  auto Run3file = TFile::Open("Run3_CrossSection_TrackingEfficiency.root", "read");
  auto hRun3Xsection = (TH1 *) Run3file->Get("Run3_CrossSection");
  return hRun3Xsection;
}
void ALICEfigureLegend(const char *FigureLabel, double x11, double y11, double x12, double y12, double x21, double y21, double x22, double y22, double textsize=0.043) {
      TLegend *ALICEleg2 =
          new TLegend(x21, y21, x22, y22, NULL,"brNDC");
      ALICEleg2->SetTextSize(textsize);
      ALICEleg2->SetBorderSize(0);
      ALICEleg2->SetTextAlign(12);
      ALICEleg2->SetFillColorAlpha(0,0); 
      ALICEleg2->AddEntry("", "|#it{#eta}_{jet}| < 0.5", "");
      ALICEleg2->AddEntry("", "Anti-#it{k}_{T}, #it{R} = 0.4", "");
      ALICEleg2->AddEntry("", "charged-particle jets", "");
      TLegend *ALICEleg1 =
          new TLegend(x11, y11, x12, y12, NULL,"brNDC");
      ALICEleg1->SetTextSize(textsize);
      ALICEleg1->SetBorderSize(0);
      ALICEleg1->SetTextAlign(12);
      ALICEleg1->SetFillColorAlpha(0,0); 
      ALICEleg1->AddEntry("", FigureLabel, "");
      ALICEleg1->AddEntry("", "pp #sqrt{#it{s}} = 13.6 TeV", "");
      ALICEleg1->AddEntry("", "#it{p}_{T, track} > 0.15 GeV/#it{c}", "");
      ALICEleg1->AddEntry("", "|#it{#eta}_{track}| < 0.9", "");

      ALICEleg1->Draw();
      ALICEleg2->Draw();
}