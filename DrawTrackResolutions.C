#include "BSHelper.cxx"
#include "Filipad2.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH3.h"
#include "TMath.h"
#include <__config>
#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;

// const char* MCfileNames =
// "../../jets/mc/jetfinderQA/AnalysisResults/LHC23d1k/AnalysisResults_hy_MatchColl_trkeff_trkres_150MeVtrkcut.root";
std::vector<TString> DatafileNames = {
    "../../jets/data/trackefficiency/AnalysisResults/LHC22o_apass6_MB_small/sel8Full_globalTracks/AnalysisResults.root",
    // "../../jets/data/trackefficiency/AnalysisResults/LHC22o_apass7_MB_sampling/sel8Full_globalTracks/AnalysisResults.root"
    // "../../jets/data/AnalysisResults/LHC22o_apass7_minBias_sampling/sel8/trackingEfficiency/AnalysisResults.root"
    "../../jets/data/AnalysisResults/LHC22o_apass7_minBias_small/sel8/Track100GeV/AnalysisResults.root"
};
std::vector<TString> MCfileNames = {
    // "/Users/jbae/cernbox/workspace/O2Physics/jets/mc/trackefficiency/AnalysisResults/LHC23d1k/AnalysisResults_HY_global_uniform_dcaZmax_2.root",
    // "/Users/jbae/cernbox/workspace/O2Physics/jets/mc/trackefficiency/AnalysisResults/LHC23d1k/AnalysisResults_HY_global_uniform_dcaZmax_2.root"
    // "/Users/jbae/cernbox/workspace/O2Physics/jets/mc/trackefficiency/AnalysisResults/LHC24d1_std/AnalysisResults_hy.root",
    // "/Users/jbae/cernbox/workspace/O2Physics/jets/mc/trackefficiency/AnalysisResults/LHC24d1_std/AnalysisResults_hy.root",
    // "/Users/jbae/cernbox/workspace/O2Physics/jets/mc/jetfinderQA/AnalysisResults/LHC24b1b/AnalysisResults_hy.root",
    "../../jets/mc/trackefficiency/AnalysisResults/LHC24b1b/selMC_globalTracks/AnalysisResults.root",
    // "../../jets/mc/trackefficiency/AnalysisResults/LHC24b1b/selMC_uniformTracks/AnalysisResults.root"
    // "../../jets/mc/trackefficiency/AnalysisResults/LHC24f3/selMC/AnalysisResults.root"
    "../../jets/mc/AnalysisResults/LHC24f3b/selMC/trackTuner/Track100GeV/AnalysisResults.root" // track tuned, track pT < 100 GeV, most recent pass7 anchored to MB, MCP: selMC w/o zvtx
    // "../../jets/mc/AnalysisResults/LHC24f3/selMC/AnalysisResults.root"
    };
const char* Dir = "track-efficiency";
// std::vector<TString> Dir = {"track-efficiency", "track-efficiency_uniform"};
// const char* HistName = "globalTracks"
std::vector<TString> DatahistNames = {"LHC22o_apass6_small", "Data(LHC22o)"};
std::vector<TString> MChistNames = {"LHC24b1b", "MC(LHC24f3b)"};
double Trackptbin[12] = {0, 2, 4, 6, 8, 10, 15, 20, 25, 30, 40, 50};
int nTrackptbin = sizeof(Trackptbin) / sizeof(Trackptbin[0]) - 1;

const char *TrackResolutionObj =
    "h2_particle_pt_track_pt_diif_associatedtrack_primary";
const char *TrackResolutionObj2 =
    "h2_particle_pt_track_pt_residual_associatedtrack_primary";
const char *TrackResolutionHighObj = "h2_particle_pt_high_track_pt_high_residual_associatedtrack_primary";
// const char *RatioTitleY = "MC / Data";
const char *TrackResolutionTitleX = "#it{p}_{T, track}^{Gen} - #it{p}_{T, "
                                    "track}^{Reco} / #it{p}_{T, track}^{Gen}";
const char *TrackResolutionTitleY = "1/N dN/d#it{p}_{T}";

std::vector<int> ColorPalette = {
    kRed,         kBlue,       kGreen,     kMagenta,   kCyan,
    kOrange,      kYellow + 2, kAzure + 7, kViolet,    kSpring,
    kPink,        kTeal,       kAzure,     kOrange + 3, kSpring + 8,
    kMagenta + 3, kYellow - 3, kRed - 4,   kGreen - 5, kBlue - 6
};
Int_t n = 0;
Int_t nn = 0;
void setpad(TVirtualPad *pad) {
  pad->SetTopMargin(0.02);
  pad->SetLeftMargin(0.13);
  pad->SetRightMargin(0.15);
  pad->SetBottomMargin(0.13);
  pad->SetName(Form("c%d", ++n));
}
void hset(TH1 &hid, TString xtit = "", TString ytit = "", double titoffx = 0.9,
          double titoffy = 1.2, double titsizex = 0.06, double titsizey = 0.06,
          double labeloffx = 0.01, double labeloffy = 0.001,
          double labelsizex = 0.05, double labelsizey = 0.05, int divx = 510,
          int divy = 510) { // hid.SetStats(0);

  // hid.GetXaxis()->CenterTitle(1);
  // hid.GetYaxis()->CenterTitle(1);

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
void hoptset(TH1 &hid, Float_t N = 1, Color_t color = kBlack, Double_t minX = 0,
             Double_t maxX = 100, Double_t minY = 0, Double_t maxY = 1,
             Double_t MarkerSize = 1) {
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
  hid.SetMarkerStyle(24);

  auto minXbin = hid.GetXaxis()->FindBin(minX);
  auto maxXbin = hid.GetXaxis()->FindBin(maxX);
  auto minYbin = hid.GetYaxis()->FindBin(minY);
  auto maxYbin = hid.GetYaxis()->FindBin(maxY);

  hid.GetXaxis()->SetRangeUser(minX, maxX);
  hid.GetYaxis()->SetRangeUser(minY, maxY);
  hid.SetFillColorAlpha(color, 0.3);
  hid.GetYaxis()->SetNdivisions(510);
}
void optFili(TPad &pid, Int_t gridx, Int_t gridy, Int_t logx, Int_t logy,
             Int_t logz) {
  pid.SetGridy(gridx);
  pid.SetGridx(gridy);
  pid.SetLogx(logx);
  pid.SetLogy(logy);
  pid.SetLogz(logz);
}

TH1 *H1TrkResCovMat(TString path, TString directory, TString histname, Color_t color, TLegend *legtrkres);
TH1 *DrawTrkResMC(const char *fileName, TString directory, const char *histName, Color_t color, TLegend *legtrkres);
// void H2TrackResolution(TString path, TString directory, TString histname,
//                        const double *trackptbin, int ntrackptbins);

void DrawTrackResolutions() {
  for (Int_t filei = 0;
       filei <
       std::min(DatafileNames.size(), std::min(MCfileNames.size(), std::min(DatahistNames.size(), MChistNames.size())));
       filei++) {
    
    auto canTrkRes = (TCanvas*) new TCanvas(Form("TrackMomentumResolution_%s", MCfileNames[filei].Data()), Form("TrackMomentumResolution_%s", MCfileNames[filei].Data()), 1000, 500);
    canTrkRes->cd();
    setpad(canTrkRes);
    TLegend *legtrkres = new TLegend(0.2,0.791579,0.44489,0.943158,NULL,"brNDC");
    legtrkres->SetTextSize(0.05);
    legtrkres->SetTextAlign(12);
    legtrkres->SetBorderSize(0);
    legtrkres->SetFillColorAlpha(0, 0);
    optFili(*canTrkRes, 1, 1, 0, 0, 0);

    TH1 *h1 = H1TrkResCovMat(DatafileNames[filei], Dir, DatahistNames[filei], kBlack, legtrkres);
    legtrkres->AddEntry(h1, Form("Data: covariance matrix (%s)", DatahistNames[filei].Data()), "pe");

    TH1 *h2 = H1TrkResCovMat(MCfileNames[filei], Dir, MChistNames[filei], kRed, legtrkres);
    legtrkres->AddEntry(h2, Form("MC: covariance matrix (%s)", MChistNames[filei].Data()), "pe");
    
    DrawTrkResMC(MCfileNames[filei], Dir, MChistNames[filei], kBlue, legtrkres);

    legtrkres->Draw();

    canTrkRes->SaveAs(Form("plots/AN_Charged-particle-jet-cross-section-in-pp-collisions-at-13.6-TeV/Figures/trackQA/%s_%s.pdf", DatahistNames[filei].Data(), MChistNames[filei].Data()));
  }
}

TH1 *H1TrkResCovMat(TString path, TString directory, TString histname, Color_t color, TLegend *legtrkres) {
    auto file = TFile::Open(path, "open");
    auto dir = (TDirectory *)file->Get(directory.Data());
    dir->cd();

    gStyle->SetOptStat(0);

    TH2 *h2trk_pt_trk_sigmapt = (TH2 *)dir->Get("h2_track_pt_track_sigmapt");
    TH1 *h1TrkResLow = new TH1F(Form("h1TrkResLow_%s", histname.Data()), Form("h1TrkResLow_%s", histname.Data()), 100, 0, 10);

    for (Int_t pti = 1; pti <= h2trk_pt_trk_sigmapt->GetNbinsX(); pti++) {
        TH1 *h1_trkpt_sigmapt = (TH1 *)h2trk_pt_trk_sigmapt->ProjectionY(Form("h1_trkpt_sigmapt_%s", histname.Data()), pti, pti, "e");
        // double binCenter = h2trk_pt_trk_sigmapt->GetXaxis()->GetBinCenter(pti);
        double binCenter = 1;
        double meansigma = h1_trkpt_sigmapt->GetMean() / binCenter;
        double errorsigma = h1_trkpt_sigmapt->GetMeanError() / binCenter;
        h1TrkResLow->SetBinContent(pti, meansigma);
        h1TrkResLow->SetBinError(pti, errorsigma);
    }

    // high pT 
    TH2 *h2trk_pt_high_trk_sigmapt = (TH2 *)dir->Get("h2_track_pt_high_track_sigmapt");
    TH1 *h1TrkResHigh = new TH1F(Form("h1TrkResHigh_%s", histname.Data()), Form("h1TrkResHigh_%s", histname.Data()), 90, 10, 100);

    for (Int_t pti = 1; pti <= h2trk_pt_high_trk_sigmapt->GetNbinsX(); pti++) {
        TH1 *h1_trkpthigh_sigmapt = (TH1 *)h2trk_pt_high_trk_sigmapt->ProjectionY(Form("h1_trkpthigh_sigmapt_%s", histname.Data()), pti, pti, "e");
        // double binCenter = h2trk_pt_high_trk_sigmapt->GetXaxis()->GetBinCenter(pti);
        double binCenter = 1;
        double meansigma = h1_trkpthigh_sigmapt->GetMean() / binCenter;
        double errorsigma = h1_trkpthigh_sigmapt->GetMeanError() / binCenter;
        h1TrkResHigh->SetBinContent(pti, meansigma);
        h1TrkResHigh->SetBinError(pti, errorsigma);
    }

    // new bin boundaries 
    std::vector<double> binEdges;
    for (double i = 0; i <= 10; i += 0.1) {
        binEdges.push_back(i);
    }
    for (double i = 11; i <= 100; i += 1) {
        binEdges.push_back(i);
    }

    TH1 *h1Combined = new TH1F(Form("h1Combined_%s", histname.Data()), Form("h1Combined_%s", histname.Data()), binEdges.size() - 1, &binEdges[0]);

    for (Int_t pti = 1; pti <= h1TrkResLow->GetNbinsX(); pti++) {
        double binCenter = h1TrkResLow->GetXaxis()->GetBinCenter(pti);
        h1Combined->Fill(binCenter, h1TrkResLow->GetBinContent(pti));
        int combinedBin = h1Combined->FindBin(binCenter);
        h1Combined->SetBinError(combinedBin, h1TrkResLow->GetBinError(pti));
    }

    for (Int_t pti = 1; pti <= h1TrkResHigh->GetNbinsX(); pti++) {
        double binCenter = h1TrkResHigh->GetXaxis()->GetBinCenter(pti);
        h1Combined->Fill(binCenter, h1TrkResHigh->GetBinContent(pti));
        int combinedBin = h1Combined->FindBin(binCenter);
        h1Combined->SetBinError(combinedBin, h1TrkResHigh->GetBinError(pti));
    }

    h1Combined->SetMarkerColor(color);
    h1Combined->SetLineColor(color);
    h1Combined->SetMarkerStyle(24);
    h1Combined->SetTitle("");
    hset(*h1Combined, "#it{p}_{T, track}^{true} (GeV/#it{c})", "#sigma(#it{p}_{T})/#it{p}_{T}", 1.0, 0.7, 0.06, 0.075, 0.01, 0.01, 0.05, 0.05, 510, 505);
    // h1Combined->GetXaxis()->SetTitle("#it{p}_{T, track} (GeV/#it{c})");
    // h1Combined->GetYaxis()->SetTitle("#sigma(#it{p}_{T})/#it{p}_{T}");

    hoptset(*h1Combined, 0, color, 0, 20, 0, 0.1, 1);
    h1Combined->Draw("PESAME");
    return h1Combined;
}

TH1 *DrawTrkResMC(const char *fileName, TString directory, const char *histName, Color_t color, TLegend *legtrkres) {
auto file = TFile::Open(fileName, "open");
auto dir = (TDirectory *)file->Get(directory.Data());
dir->cd();
TString histNameStr(histName);

TH2 *TrackResolution = nullptr;
TH2 *TrackResolutionTemp1 = (TH2 *)file->Get(Form("%s/%s", directory.Data(), TrackResolutionObj));
TH2 *TrackResolutionTemp2 = (TH2 *)file->Get(Form("%s/%s", directory.Data(), TrackResolutionObj2));
if (TrackResolutionTemp1) {
TrackResolution = TrackResolutionTemp1;
} else if (TrackResolutionTemp2) {
TrackResolution = TrackResolutionTemp2;
}

TH1* TrkMomRes = new TH1F(Form("TrkMomResMC_%s", histName), Form("TrkMomResMC_%s", histName), 500, 0., 10.);
TH1* TrkMomResFit = new TH1F(Form("TrkMomResFitMC_%s", histName), Form("TrkMomResFitMC_%s", histName), 500, 0., 10.);

for (Int_t pTi = 1; pTi <= TrackResolution->GetNbinsX(); pTi++) {
// for (Int_t pTi = 1; pTi <= 10; pTi++) {
  TH1D *trackResolution = (TH1D *)TrackResolution->ProjectionY(Form("TrackResY%i_%s", pTi + 1, histName), pTi, pTi);
  TF1 *fitfunc = new TF1("fit", "gaus", -2, 2);

  auto BinContent = (Double_t) trackResolution->GetRMS();
  auto BinError = (Double_t) trackResolution->GetMeanError();

  // new TCanvas;
  // trackResolution->Draw("pe");
  trackResolution->Fit(fitfunc, "QN");
  Double_t sigma = fitfunc->GetParameter(2);
  Double_t sigmaError = fitfunc->GetParError(2);

  TrkMomRes->SetBinContent(pTi, BinContent);
  TrkMomRes->SetBinError(pTi, BinError);

  TrkMomResFit->SetBinContent(pTi, sigma);
  TrkMomResFit->SetBinError(pTi, sigmaError);
}

hoptset(*TrkMomRes, 0, color, 0, 100, 0, 0.1, 1);
// TrkMomRes->Draw("pesame");
hoptset(*TrkMomResFit, 0, kBlue, 0, 20, 0, 0.1, 1);
TrkMomResFit->Draw("pesame");

TH2 *TrackResolutionHight = nullptr;
TrackResolutionHight = (TH2 *)file->Get(Form("%s/%s", directory.Data(), TrackResolutionHighObj));
if (TrackResolutionHight!=nullptr) {
  TH1* TrkMomResHigh = new TH1F(Form("TrkMomResHighMC_%s", histName), Form("TrkMomResHighMC_%s", histName), 18, 10., 100.);
  TH1* TrkMomResHighFit = new TH1F(Form("TrkMomResHighFitMC_%s", histName), Form("TrkMomResHighFitMC_%s", histName), 18, 10., 100.);

  for (Int_t pTj = 1; pTj <= TrackResolutionHight->GetNbinsX(); pTj++) {
  // for (Int_t pTj = 1; pTj <= 10; pTj++) {
    TH1D *trackResolution = (TH1D *)TrackResolutionHight->ProjectionY(Form("TrackResY%i_%s", pTj + 1, histName), pTj, pTj);
    TF1 *fitfunc = new TF1("fit", "gaus", -2, 2);

    auto BinContent = (Double_t) trackResolution->GetRMS();
    auto BinError = (Double_t) trackResolution->GetMeanError();

    // new TCanvas;
    // trackResolution->Draw("pe");
    trackResolution->Fit(fitfunc, "QN");
    Double_t sigma = fitfunc->GetParameter(2);
    Double_t sigmaError = fitfunc->GetParError(2);

    TrkMomResHigh->SetBinContent(pTj, BinContent);
    TrkMomResHigh->SetBinError(pTj, BinError);

    TrkMomResHighFit->SetBinContent(pTj, sigma);
    TrkMomResHighFit->SetBinError(pTj, sigmaError);
  }

  hoptset(*TrkMomResHigh, 0, color, 0, 100, 0, 0.1, 1);
  // TrkMomResHigh->Draw("pesame");
  hoptset(*TrkMomResHighFit, 0, kBlue, 0, 20, 0, 0.1, 1);
  TrkMomResHighFit->Draw("pesame");
}

// legtrkres->AddEntry(TrkMomRes, Form("MC: compare #it{p}_T^{gen} and #it{p}_T^{rec} (%s)", histName), "pe");
legtrkres->AddEntry(TrkMomResFit, Form("MC: #it{p}_{T} residual fit (%s)", histName), "pe");

return 0;
}