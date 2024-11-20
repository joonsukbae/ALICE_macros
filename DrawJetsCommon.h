///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Draw macro for jet QA         //////////
////////// author: Joonsuk Bae           ////////// 
////////// E-mail: jbae@cern.ch          //////////
////////// Last Modified: 08 Nov 2024    //////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

#include "BSHelper.cxx"
#include "Filipad2.h"
#include <__config>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <complex>
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
using namespace std;

std::vector<TString> fileNames;
std::vector<TString> histNames;

bool DRAWPLOTS;
double PlotPtMin;
double PlotPtMax;

// Double_t Trackptbin[17] = {0.15,  2,  4,  6,  8,  10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100};
Double_t Trackptbin[21] = {0.15,  2,  4,  6,  8,  10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 170, 200};
Double_t ptbin[21] = {5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
// Double_t ptbin[27] = {0, 1, 2, 3, 4, 5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200, 300};
Double_t ptbinGen[26] = {0, 1, 2, 3, 4, 5,  6,  7,  8,  9,  10, 12, 14, 16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
// Double_t ptbinGen[27] = {0, 1, 2, 3, 4, 5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200, 300};
Int_t nTrackptbin = sizeof(Trackptbin) / sizeof(Trackptbin[0]) - 1;
Int_t nptBins = sizeof(ptbin) / sizeof(ptbin[0]) - 1;
Int_t nptBinsGen = sizeof(ptbinGen) / sizeof(ptbinGen[0]) - 1;

std::vector<Color_t> ColorPallete = {
    kBlack,         kRed,         kBlue + 1, kGreen + 2, kOrange + 7,
    kMagenta + 2,   kTeal + 3,    kViolet + 2, kYellow + 3, kCyan - 6,
    kAzure + 2,     kPink - 7,    kSpring + 5, kGray + 2,  kAzure + 8
};

Double_t RBIN = 0.4;

Int_t n = 0;
Int_t nn = 0;
Int_t ii = 0;
void setpad(TVirtualPad *pad, Double_t tmargin = 0.02, Double_t bmargin = 0.15,
            Double_t lmargin = 0.13, Double_t rmargin = 0.05) {
  pad->SetTopMargin(tmargin);
  pad->SetLeftMargin(lmargin);
  pad->SetRightMargin(rmargin);
  pad->SetBottomMargin(bmargin);
  pad->SetName(Form("c%d", ++n));
}
template <typename T>
void hset(T &hid, TString xtit = "", TString ytit = "", double titoffx = 0.9,
          double titoffy = 1.2, double titsizex = 0.06, double titsizey = 0.06,
          double labeloffx = 0.01, double labeloffy = 0.001,
          double labelsizex = 0.05, double labelsizey = 0.05, int divx = 510,
          int divy = 510) {
  // hid.SetStats(0);

  hid.GetXaxis()->CenterTitle(1);
  hid.GetYaxis()->CenterTitle(1);

  hid.GetXaxis()->SetTitleOffset(titoffx);
  hid.GetYaxis()->SetTitleOffset(titoffy);

  hid.GetXaxis()->SetTitleSize(titsizex);
  hid.GetYaxis()->SetTitleSize(titsizey);

  hid.GetXaxis()->SetLabelOffset(labeloffx);
  hid.GetYaxis()->SetLabelOffset(labeloffy);

  hid.GetXaxis()->SetLabelSize(labelsizex);
  hid.GetYaxis()->SetLabelSize(labelsizey);

  hid.GetXaxis()->SetNdivisions(divx);
  hid.GetYaxis()->SetNdivisions(divy);

  hid.GetXaxis()->SetTitle(xtit);
  hid.GetYaxis()->SetTitle(ytit);
}
template <typename T>
void hoptset(T &hid, Double_t N = 1, Color_t color = kBlack, Double_t minX = 0,
             Double_t maxX = 100, Double_t minY = 0, Double_t maxY = 1,
             Double_t MarkerSize = .75, Int_t LineStyle = 1, Int_t LineWidth = 1, Int_t MarkerStyle = 20) {
  if (N != 0) {
    if (N == 1) {
      hid.Scale(1. / hid.Integral(), "width");
    } else if (N == 2) {
      hid.Scale(1. / hid.Integral(), "");
    } else {
      hid.Scale(1. / N, "width");
    }
  }

  hid.SetMarkerColor(color);
  hid.SetLineColor(color);
  hid.SetMarkerSize(MarkerSize);
  hid.SetMarkerStyle(MarkerStyle); 
  hid.SetLineStyle(LineStyle); 
  hid.SetLineWidth(LineWidth);

  hid.GetXaxis()->SetRangeUser(minX, maxX);
  hid.GetYaxis()->SetRangeUser(minY, maxY);

  hid.SetFillColorAlpha(color, 0.3);

}
void optFili(TPad &pid, Int_t gridx, Int_t gridy, Int_t logx, Int_t logy) {
  pid.SetGridy(gridx);
  pid.SetGridx(gridy);
  pid.SetLogx(logx);
  pid.SetLogy(logy);
}

std::vector<Double_t> NewBin(Int_t Nbins = 100, Double_t minBin = 0,
                             Double_t maxBin = 100) {
  std::vector<Double_t> Bin(Nbins + 1);
  for (Int_t i = 0; i <= Nbins; i++) {
    Bin[i] = minBin + i * (maxBin - minBin) / Nbins;
  }
  return Bin;
}

void ALICEfigureLegend(const char *FigureLabel, double x11, double y11, double x12, double y12, double x21, double y21, double x22, double y22, double textsize=0.043);

// define fns
Double_t Nevents(const char *fileName, const char *eventDir,
              const char *eventObj, Int_t ifMCP = 0) {
  auto file = TFile::Open(fileName, "open");
  auto Nevents = (TH1D *)file->Get(Form("%s/%s", eventDir, eventObj));

  Double_t nevents = 0;
  if (ifMCP == 0) {
    nevents = Nevents->GetBinContent(Nevents->FindBin(1.5));
  } else if (ifMCP == 1) {
    nevents = Nevents->GetBinContent(Nevents->FindBin(0.5));
  }
  return nevents;
}
void PrintHistogramErrors(TH1 *hist) {
  int nBins = hist->GetNbinsX();
  for (int i = 1; i <= nBins; ++i) { // bin id starts from 1
    double binContent = hist->GetBinContent(i);
    double binError = hist->GetBinError(i);
  }
}

TH1 *DrawRatio(const char *ratioName, TH1 *refHist, TH1 *testHist,
               TString AxisTitleX, TString AxisTitleY, Color_t colorID,
               Double_t Ymin, Double_t Ymax, Double_t MarkerSize = 1) {
  TH1 *ratioHist = (TH1 *)refHist->Clone(ratioName);
  hset(*ratioHist, AxisTitleX, AxisTitleY, 1.2, 1.0, 0.07, 0.07, 0.01, 0.01,
       0.07, 0.07, 510, 505);
  ratioHist->Divide(testHist, ratioHist, 1., 1., "B");
  ratioHist->GetYaxis()->SetRangeUser(Ymin, Ymax);
  // ratioHist->GetYaxis()->SetLimits(Ymin, Ymax);
  ratioHist->SetMarkerColor(colorID);
  ratioHist->SetLineColor(colorID);
  ratioHist->SetMarkerSize(MarkerSize);
  ratioHist->Draw("esame");

  return ratioHist;
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

TGraphErrors *Run2Data(TH1 *hist) {
  // auto Run2Datafile = TFile::Open("~/cernbox/workspace/O2Physics/jets/Run2_YongzhenHou/HEPData-ins2026265-v1-Figure_3.root", "OPEN");// w/ UE sub.
  auto Run2Datafile = TFile::Open("~/cernbox/workspace/O2Physics/jets/Run2_YongzhenHou/HEPData-ins2026265-v1-Figure_A1.root", "OPEN");// w/o UE sub.
  // auto Run2DataGraph = (TGraphErrors *)Run2Datafile->Get("Figure 3/Graph1D_y3");
  auto Run2DataGraph = (TGraphErrors *)Run2Datafile->Get("Figure A1/Graph1D_y3");

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

std::vector<TGraphErrors*> Run2Data_manual() {
    const int Npoint = 19; 
    
    double pt[Npoint] = {5.5, 6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17, 19, 22.5, 27.5, 35, 45, 55, 65, 77.5, 92.5, 120};
    
    double sigma[Npoint] = {1.0972, 0.58418, 0.33224, 0.20001, 0.12614, 0.069825, 0.034225, 0.018498, 0.010764, 0.0066421, 
                            0.0033177, 0.0013313, 0.00046794, 0.00013985, 5.2994e-05, 2.3472e-05, 1.0087e-05, 4.1091e-06, 1.2028e-06};

    double stat_err[Npoint] = {0.00016396, 0.00010887, 7.5637e-05, 5.3822e-05, 3.7416e-05, 2.2965e-05, 1.3182e-05, 8.6563e-06, 5.8873e-06, 4.1364e-06,
                               2.2742e-06, 1.136e-06, 4.9698e-07, 1.7976e-07, 7.9569e-08, 3.9801e-08, 1.8747e-08, 8.1518e-09, 2.4881e-09};

    double syst_err[Npoint] = {0.074335, 0.042393, 0.024799, 0.015407, 0.010138, 0.0057313, 0.0029103, 0.001585, 0.00094719, 0.00061878, 
                               0.00031195, 0.00012884, 4.6138e-05, 1.4546e-05, 5.6228e-06, 2.5033e-06, 1.1183e-06, 4.4943e-07, 1.3178e-07};

    TGraphErrors *gr_stat_error = new TGraphErrors(Npoint);
    TGraphErrors *gr_syst_error = new TGraphErrors(Npoint);

    for (int i = 0; i < Npoint; i++) {
        double ex1 = (ptbinGen[i+1] - ptbinGen[i]) / 2.0;

        gr_stat_error->SetPoint(i, pt[i], sigma[i]);
        gr_stat_error->SetPointError(i, ex1, stat_err[i]); 

        gr_syst_error->SetPoint(i, pt[i], sigma[i]);
        gr_syst_error->SetPointError(i, ex1, syst_err[i]); 
    }

    std::vector<TGraphErrors*> graphs = {gr_stat_error, gr_syst_error};
    return graphs;
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
// void DrawMultipleSources(std::vector<TH1*>& hSysts, TH1* hSystResult) {
//     if (hSysts.size() < 7) {
//         std::cerr << "Error: There should be at least 7 histograms!" << std::endl;
//         return;
//     }
//     TCanvas* cSources = new TCanvas("cSources", "Systematic Uncertainty Sources", 800, 600);
//     setpad(cSources, 0.02, 0.15, 0.12, 0.03);
//     cSources->cd();
//     TLegend* legend = new TLegend(0.622807,0.648696,0.799499,0.93913,NULL,"brNDC");
//     legend->SetTextSize(0.043);
//     legend->SetBorderSize(0);
//     legend->SetFillColorAlpha(0, 0);
//     bool firstHist = true;
//     std::vector<TH1*> newHSysts;
//     std::vector<TString> SourceNames = {"Tracking efficiency", "Track #it{p}_{T} resolution", "Unfolding", "Normalization", "Secondary particles", "Total uncertainty"};
//     // Step 1: Create the Unfolding source by combining the histograms [2], [3], and [4] using Quadrature Sum
//     TH1* hUnfolding = (TH1*)hSysts[2]->Clone("hUnfolding");
//     hUnfolding->Reset();
//     for (int i = 1; i <= hUnfolding->GetNbinsX(); ++i) {
//         double sumOfSquares = 0.0;
//         for (int j = 2; j <= 4; ++j) {
//             double value = (hSysts[j]->GetBinContent(i) - 1.0) * 100;  // Subtract 1 and convert to percentage
//             sumOfSquares += value * value;  // Add the square of the normalized value
//         }
//         double quadratureSum = std::sqrt(sumOfSquares);  // Calculate quadrature sum
//         hUnfolding->SetBinContent(i, quadratureSum);  // Set the value for the unfolding source in %
//     }
//     // Step 2: Convert all sources to percentage by subtracting 1 and multiplying by 100
//     for (int i = 0; i < hSysts.size(); ++i) {
//         if (i == 2 || i == 3 || i == 4) continue;  // Skip the sources that were quadrature summed

//         for (int bin = 1; bin <= hSysts[i]->GetNbinsX(); ++bin) {
//             double newValue = (hSysts[i]->GetBinContent(bin) - 1.0) * 100;  // Subtract 1 and convert to percentage
//             hSysts[i]->SetBinContent(bin, newValue);
//         }
//     }
//     // Step 3: Convert hSystResult (Total Uncertainty) to percentage
//     auto hSystTotal = (TH1 *) hSystResult->Clone("hSystTotal");
//     for (int bin = 1; bin <= hSystTotal->GetNbinsX(); ++bin) {
//         double newValue = (hSystTotal->GetBinContent(bin) - 1.0) * 100;  // Subtract 1 and convert to percentage
//         hSystTotal->SetBinContent(bin, newValue);
//     }
//     // Step 4: Add the sources to the new histogram vector
//     newHSysts.push_back(hSysts[0]);  // track efficiency
//     newHSysts.push_back(hSysts[1]);  // track pT resolution
//     newHSysts.push_back(hUnfolding); // Unfolding (quadrature sum of [2], [3], and [4])
//     newHSysts.push_back(hSysts[5]);  // Normalization
//     newHSysts.push_back(hSysts[6]);  // Secondary Particles
//     newHSysts.push_back(hSystTotal); // Total Uncertainty (smoothing된 total error)
//     // Step 5: Draw each histogram in the newHSysts vector
//     for (size_t i = 0; i < newHSysts.size(); ++i) {
//         TH1* hist = newHSysts[i];
//         TString sourceName = SourceNames[i];
//         hist->SetLineWidth(5);
//         hist->SetLineStyle(i);
//         hist->SetMarkerStyle(0);  // No markers
//         hist->SetFillStyle(0);    // No fill
//         hist->SetLineColor(ColorPallete[i]);  // Use the color palette
//         // Custom histogram options
//         hoptset(*hist, 0, ColorPallete[i], PlotPtMin, PlotPtMax, 0, 60.0, 0.6, 1, 2, 25);  // Adjust y-axis for percentage
//         hset(*hist, JetPtDataFinalTitleX, "Rel. Uncertainties (%)", 1.3, 1.1, 0.05, 0.05, 0.01, 0.01,
//              0.05, 0.06, 510, 505);
//         if (firstHist) {
//             hist->Draw("HIST");  // Draw first histogram normally
//             firstHist = false;
//         } else {
//             hist->Draw("HIST SAME");  // Overlay the rest
//         }
//         // Add entry to the legend
//         legend->AddEntry(hist, sourceName, "l");
//     }
//     // Draw the legend
//     legend->Draw();
//     ALICEfigureLegend("ALICE Preliminary", 0.110276,0.68,0.309524,0.96,0.109023,0.429565,0.309524,0.638261, 0.04);
//     // Save the canvas if the DRAWPLOTS flag is set
//     if (DRAWPLOTS) {
//         cSources->SaveAs(Form("%s/SystematicUncertaintySources.pdf", MakeDirName.Data()));
//     }
// }
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
  double minPT = hHerwigPythiaRat->GetXaxis()->GetXmin();
  double maxPT = hHerwigPythiaRat->GetXaxis()->GetXmax();

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