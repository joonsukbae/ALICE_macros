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

// const char* fileNames =
// "../../jets/mc/jetfinderQA/AnalysisResults/LHC23d1k/AnalysisResults_hy_MatchColl_trkeff_trkres_150MeVtrkcut.root";
std::vector<TString> fileNames = {
    // "/Users/jbae/cernbox/workspace/O2Physics/jets/mc/trackefficiency/AnalysisResults/LHC23d1k/AnalysisResults_HY_global_uniform_dcaZmax_2.root",
    // "/Users/jbae/cernbox/workspace/O2Physics/jets/mc/trackefficiency/AnalysisResults/LHC23d1k/AnalysisResults_HY_global_uniform_dcaZmax_2.root"
    // "/Users/jbae/cernbox/workspace/O2Physics/jets/mc/trackefficiency/AnalysisResults/LHC24d1_std/AnalysisResults_hy.root",
    // "/Users/jbae/cernbox/workspace/O2Physics/jets/mc/trackefficiency/AnalysisResults/LHC24d1_std/AnalysisResults_hy.root",
    // "/Users/jbae/cernbox/workspace/O2Physics/jets/mc/jetfinderQA/AnalysisResults/LHC24b1b/AnalysisResults_hy.root",
    "../../jets/mc/trackefficiency/AnalysisResults/LHC24b1b/selMC_uniformTracks/AnalysisResults.root",
    "../../jets/mc/trackefficiency/AnalysisResults/LHC24b1b/selMC_uniformTracks/AnalysisResults.root"
    };
// const char* Dir = "track-efficiency";
std::vector<TString> Dir = {"track-efficiency", "track-efficiency_uniform"};
// const char* HistName = "globalTracks"
std::vector<TString> histNames = {"LHC24b1b_selMC_globalTracks", "LHC24b1b_selMC_uniformTracks"};
double Trackptbin[12] = {0, 2, 4, 6, 8, 10, 15, 20, 25, 30, 40, 50};
int nTrackptbin = sizeof(Trackptbin) / sizeof(Trackptbin[0]) - 1;

const char *TrackResolutionObj =
    "h2_particle_pt_track_pt_diif_associatedtrack_primary";
const char *TrackResolutionObj2 =
    "h2_particle_pt_track_pt_residual_associatedtrack_primary";
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
  // cout<<Form("Nevets of %s: %f",hid.GetName(),N)<<endl;
  if (N == 1) {
    // cout<<Form("Nevets of %s: %f",hid.GetName(),N)<<endl;
    hid.Scale(1. / hid.Integral(), "width");
  } else if (N == 2) {
    hid.Scale(1. / hid.Integral(), "");
  } else {
    hid.Scale(1. / N, "width");
  }
  // cout<<Form("Entries of %s:",hid.GetName())<< hid.GetEntries() <<endl;

  hid.SetMarkerColor(color);
  hid.SetLineColor(color);
  hid.SetMarkerSize(MarkerSize);
  hid.SetMarkerStyle(22);

  auto minXbin = hid.GetXaxis()->FindBin(minX);
  auto maxXbin = hid.GetXaxis()->FindBin(maxX);
  auto minYbin = hid.GetYaxis()->FindBin(minY);
  auto maxYbin = hid.GetYaxis()->FindBin(maxY);

  hid.GetXaxis()->SetRangeUser(minX, maxX);
  hid.GetYaxis()->SetRangeUser(minY, maxY);
  hid.SetFillColorAlpha(color, 0.3);
  hid.GetYaxis()->SetNdivisions(505);
  // auto deltax = 0.1*(maxX-minX);
  // auto deltay = 0.1*(maxY-minY);
  // hid.GetXaxis()->SetLimits(minX - deltax, maxX + deltax);
  // hid.GetYaxis()->SetLimits(minY - deltay, maxY + deltay);
}
void optFili(TPad &pid, Int_t gridx, Int_t gridy, Int_t logx, Int_t logy,
             Int_t logz) {
  pid.SetGridy(gridx);
  pid.SetGridx(gridy);
  pid.SetLogx(logx);
  pid.SetLogy(logy);
  pid.SetLogz(logz);
}

void H1TrkResCovMat(TString path, TString directory, TString histname,
                       const double *trackptbin, int ntrackptbins);
void DrawHistos(TString path, TString directory, TString histname,
                int Color);
// void H2TrackResolution(TString path, TString directory, TString histname,
//                        const double *trackptbin, int ntrackptbins);

void DrawTrackResolutionPt() {
  for (Int_t filei = 0;
       filei <
       std::min(fileNames.size(), std::min(Dir.size(), histNames.size()));
       filei++) {

    H1TrkResCovMat(fileNames[filei], Dir[filei], histNames[filei],
                      Trackptbin, nTrackptbin);
    
    DrawHistos(fileNames[filei], Dir[filei], histNames[filei], ColorPalette[filei]);
    // H2TrackResolution(fileNames[filei], Dir[filei], histNames[filei],
    //                   Trackptbin, nTrackptbin);
  }
}

void H1TrkResCovMat(TString path, TString directory, TString histname,
                       const double *trackptbin, int ntrackptbins) {
  auto file = TFile::Open(path, "open");
  auto dir = (TDirectory *)file->Get(directory);
  dir->cd();

  gStyle->SetOptStat(0);

  auto h2trk_pt_trk_sigmapt =
      (TH2 *)dir->Get("h2_track_pt_track_sigmapt");

  TH1* h1TrkRes = new TH1F(Form("h1TrkRes_%s", histname.Data()), Form("h1TrkRes_%s", histname.Data()), 100, 0, 100);

  for (Int_t pti = 1; pti <= h2trk_pt_trk_sigmapt->GetNbinsX(); pti++) {
    // Int_t Lptbin = h2trk_pt_trk_sigmapt->GetYaxis()->FindBin(pti);
    auto h1_trkpt_sigmapt = (TH1*) h2trk_pt_trk_sigmapt->ProjectionY(Form("h1_trkpt_sigmapt_%s", histname.Data()), pti, pti, "e");
    Double_t entrysigma = h1_trkpt_sigmapt->GetEntries();
    Double_t meansigma = h1_trkpt_sigmapt->GetMean();
    // std::cout << "entrysigma (ptbin = " << pti <<"): " << entrysigma << std::endl; 
    Double_t errorsigma = h1_trkpt_sigmapt->GetMeanError();
    h1TrkRes->SetBinContent(pti, meansigma);
    h1TrkRes->SetBinError(pti, errorsigma);
  }

  auto cantrkres = new TCanvas(Form("trk_res_%s", histname.Data()),
                               Form("trk_res_%s", histname.Data()), 900, 800);
  auto legtrkres = new TLegend(0.55, 0.665, 0.85, 0.815, NULL, "brNDC");
  legtrkres->SetTextSize(0.03);
  legtrkres->SetTextAlign(33);
  legtrkres->SetBorderSize(0);
  cantrkres->cd();
  setpad(cantrkres);
  optFili(*cantrkres, 0, 0, 0, 0, 1);

  h1TrkRes->Draw("pe");
  // legtrkres->Draw();
  cantrkres->Print(
      Form("plots/TrackEnergyResolution_%s_linearX.pdf", histname.Data()));
}

TH1 *DrawTrkResMC(const char *fileName, TString directory, const char *histName,
                         TLegend *legend, TLegend *legend2, Color_t colorID) {
  auto file = TFile::Open(fileName, "open");
  auto dir = (TDirectory *)file->Get(directory);
  dir->cd();
  TString histNameStr(histName);

  TH2 *TrackResolution = nullptr;
  TH2 *TrackResolutionTemp1 = (TH2 *)file->Get(Form("%s/%s", directory, TrackResolutionObj));
  TH2 *TrackResolutionTemp2 = (TH2 *)file->Get(Form("%s/%s", directory, TrackResolutionObj2));

  if (TrackResolutionTemp1) {
    TrackResolution = TrackResolutionTemp1;
  } else if (TrackResolutionTemp2) {
    TrackResolution = TrackResolutionTemp2;
  }
  
  if (!TrackResolution) {
    std::cout << "Error: No valid TrackResolution object found." << std::endl;
  }

  std::vector<std::pair<double, double>> pTClasses = {
      {0.15, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 7}, {7, 10}};

  legend->AddEntry("", histName, "");
  for (Int_t i = 0; i < pTClasses.size(); i++) {
    TH1D *trackResolution = (TH1D *)TrackResolution->ProjectionY(
        Form("TrackResY%i_%s", i + 1, histName),
        TrackResolution->GetXaxis()->FindBin(pTClasses[i].first),
        TrackResolution->GetXaxis()->FindBin(pTClasses[i].second));
    trackResolution->GetXaxis()->SetTitle(TrackResolutionTitleX);
    trackResolution->GetYaxis()->SetTitle(TrackResolutionTitleY);
    legend->AddEntry(trackResolution,
                     Form("%.2f< #it{p}_{T, track}^{gen} <%.2f GeV/#it{c}, #mu: %.3f, #sigma: %.3f", pTClasses[i].first,
                          pTClasses[i].second, trackResolution->GetMean(), trackResolution->GetRMS()),
                     "pe2");
    // legend2->AddEntry(trackResolution, Form(“Mean: %.2f, Std Dev: %.2f”,
    // trackResolution->GetMean(), trackResolution->GetRMS()), “pe2”);
    // hset(*trackResolution, TrackResolutionTitleX, TrackPtTitleY, 0.7, 1.0,
    // 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*trackResolution, 2., ColorPalette[i], -1., 1., 1e-8, 1e0);
    trackResolution->Draw("esame");
  }

  return 0;
}

void DrawHistos(TString path, TString directory, TString histname,
                int Color) {
    TCanvas *canTrackResolution =
        new TCanvas(Form("TrackResolution_%s", histname),
                    Form("trackresolution%s", histname), 800, 600);
    canTrackResolution->SetLogy(1);
    canTrackResolution->SetGridx(1);
    canTrackResolution->SetGridy(1);
    canTrackResolution->Draw();
    TLegend *legtrackresolution =
        new TLegend(0.3, 0.12, 0.5, 0.4, NULL, "brNDC");
    legtrackresolution->SetTextSize(0.025);
    legtrackresolution->SetTextAlign(13);
    legtrackresolution->SetBorderSize(0);
    legtrackresolution->SetFillColorAlpha(0, 0);
    // TLegend *legtrackresolutionStat = new
    // TLegend(0.302005,0.113043,0.616541,0.330435,NULL,"brNDC");
    TLegend *legtrackresolutionStat =
        new TLegend(0.0, 0.0, 0.0, 0.0, NULL, "brNDC");
    legtrackresolutionStat->SetTextSize(0.03);
    legtrackresolutionStat->SetBorderSize(0);
    TString filePath = mainDir + fileNames[i];
    DrawTrkResMC(path, directory, histname,
                        legtrackresolution, legtrackresolutionStat,
                        Color);
    legtrackresolution->Draw();
    legtrackresolutionStat->Draw();
    canTrackResolution->Print(
        Form("plots/TrackResolution_R%.1f_%s.pdf", RBIN, histname));
}

// void H2TrackResolution(TString path, TString directory, TString histname,
//                        const double *trackptbin, int ntrackptbins) {
//   auto file = TFile::Open(path, "open");
//   auto dir = (TDirectory *)file->Get(directory);
//   dir->cd();

//   gStyle->SetOptStat("e");

//   auto h3trk_pt_trk_sigmapt =
//       (TH3 *)dir->Get("h3_centrality_track_pt_track_sigmapt");
//   auto htrk_sigma = (TH1 *)h3trk_pt_trk_sigmapt->Project3D("ZY");
//   // htrk_sigma->GetXaxis()->SetRangeUser(0, 50);
//   // htrk_sigma = htrk_sigma->Rebin(nTrackptbin, "htrk_mcp", Trackptbin);

//   auto cantrkeff = new TCanvas(Form("trk_eff_%s", histname.Data()),
//                                Form("trk_eff_%s", histname.Data()), 900, 800);
//   auto legtrkeff = new TLegend(0.55, 0.665, 0.85, 0.815, NULL, "brNDC");
//   legtrkeff->SetTextSize(0.03);
//   legtrkeff->SetTextAlign(33);
//   legtrkeff->SetBorderSize(0);
//   cantrkeff->cd();
//   setpad(cantrkeff);
//   //   cantrkeff->SetLogx(1);
//   optFili(*cantrkeff, 0, 0, 1, 0, 1);

//   htrk_sigma->SetTitle("");
//   legtrkeff->SetFillColorAlpha(0, 0);
//   legtrkeff->AddEntry("", "pp #sqrt{#it{s}}=13.6 TeV", "");
//   legtrkeff->AddEntry("", "2023 MB-MC", "");
//   legtrkeff->AddEntry("", Form("%s", histname.Data()), "");
//   hset(*htrk_sigma, "#it{p}_{T, track} (GeV)", "#sigma(#it{p}_{T})/#it{p}_{T}",
//        1.2, 1.2, 0.05, 0.05, 0.01, 0.01, 0.04, 0.04);
//   htrk_sigma->GetZaxis()->SetRangeUser(1, 1e10);
//   // htrk_sigma->GetXaxis()->SetTextAlign(33);
//   // htrk_sigma->GetYaxis()->SetTextAlign(33);
//   htrk_sigma->Draw("colz");

//   legtrkeff->Draw();
//   cantrkeff->Print(
//       Form("plots/TrackResolution_sigma_pt_%s_linearX.pdf", histname.Data()));
// }