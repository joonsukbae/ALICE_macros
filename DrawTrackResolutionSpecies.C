#include "BSHelper.cxx"
#include <__config>
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

enum SPECIES { PION = 0, KAON, PROTON, ELECTRON, NUM_SPECIES };
// enum SPECIES { POSITIVE = 0, NEGATIVE, NUM_SPECIES };
const char *speciesNames[NUM_SPECIES] = {"Pion", "Kaon", "Proton", "Electron"};
// const char *speciesNames[NUM_SPECIES] = {"POSITIVE", "NEGATIVE"};

Int_t n = 0;
Int_t nn = 0;
void setpad(TVirtualPad *pad) {
  pad->SetTopMargin(0.02);
  pad->SetLeftMargin(0.13);
  pad->SetRightMargin(0.05);
  pad->SetBottomMargin(0.15);
  pad->SetName(Form("c%d", ++n));
}
void hset(TH1 &hid, TString xtit = "", TString ytit = "", double titoffx = 0.9,
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
void hoptset(TH1 &hid, Float_t N = 1, Color_t color = kBlack, Double_t minX = 0,
             Double_t maxX = 100, Double_t minY = 0, Double_t maxY = 1,
             Double_t MarkerSize = 1.) {
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
void optFili(TPad &pid, Int_t gridx, Int_t gridy, Int_t logx, Int_t logy) {
  pid.SetGridy(gridx);
  pid.SetGridx(gridy);
  pid.SetLogx(logx);
  pid.SetLogy(logy);
}

Double_t RBIN = 0.4;
const char *Dir = "track-efficiency";
const char *EventObj = "h_collisions";
TString mainDir = "../../jets/mc/trackefficiency/";
std::vector<TString> fileNames = {
    // "AnalysisResults_HY_global_uniform.root"
    // "AnalysisResults_local_pvContributorTracks.root"
    // "../../jetfinderQA/AnalysisResults/LHC23d1k/AnalysisResults_hy_trkres_uniform.root",
    // "LHC23d1k/AnalysisResults_hy_MatchColl_trkeff_trkres.root"
    // "LHC23d1k/AnalysisResults_local_species.root"
    "AnalysisResults.root"
    };
std::vector<TString> histNames = {"LHC23d1k_globalTracks"};
std::vector<int> ColorPalette = {
    kRed,         kBlue,       kGreen,     kMagenta,    kCyan,
    kOrange,      kYellow + 2, kAzure + 7, kViolet,     kSpring,
    kPink,        kTeal,       kAzure,     kOrange + 3, kSpring + 8,
    kMagenta + 3, kYellow - 3, kRed - 4,   kGreen - 5,  kBlue - 6};
const char *TrackResolutionObj =
    // "h3_charged_particle_pt_track_pt_diif_associatedtrack_primary";
    "h3_species_particle_pt_track_pt_diif_associatedtrack_primary";
// const char *RatioTitleY = "MC / Data";
const char *TrackResolutionTitleX = "#it{p}_{T, track}^{Gen} - #it{p}_{T, "
                                    "track}^{Reco} / #it{p}_{T, track}^{Gen}";
const char *TrackResolutionTitleY = "1/N dN/d#it{p}_{T}";

// declare fns
void DrawHistos(const std::vector<TString> &fileNames,
                const std::vector<TString> &histNames,
                const std::vector<int> &ColorPalette);

void DrawTrackResolutionSpecies() {
  DrawHistos(fileNames, histNames, ColorPalette);
}

// define fns
float Nevents(const char *fileName, const char *eventDir,
              const char *eventObj) {
  auto file = TFile::Open(fileName, "open");
  auto dir = (TDirectory *)file->Get(eventDir);
  dir->cd();

  auto Nevents = (TH1D *)gROOT->FindObject(eventObj);
  float nevents = Nevents->GetBinContent(Nevents->FindBin(1.5));

  return nevents;
}

TH1 *DrawTrackResolution(const char *fileName, const char *histName,
                         Int_t particleId, TLegend *legend, TLegend *legend2,
                         Color_t colorID) {
  auto file = TFile::Open(fileName, "open");
  auto dir = (TDirectory *)file->Get(Dir);
  dir->cd();
  TString histNameStr(histName);

  TH3 *h3TrackResolution =
      (TH3 *)file->Get(Form("%s/%s", Dir, TrackResolutionObj));
  if (!h3TrackResolution) {
    std::cout << "Error: No valid TrackResolution object found." << std::endl;
  }

  std::vector<std::pair<double, double>> pTClasses = {
      {0.15, 1}, {1, 3}, {3, 10}};

  legend->AddEntry(
      "", (std::string(histName) + "_" + speciesNames[particleId]).c_str(), "");
  // for (Int_t species = PION; species <= ELECTRON; species++) {
  Int_t h3LBIN = h3TrackResolution->GetXaxis()->FindBin(particleId - 0.499999);
  Int_t h3RBIN = h3TrackResolution->GetXaxis()->FindBin(particleId + 0.499999);
  std::cout << "h3LBIN: " << h3LBIN << std::endl;
  std::cout << "h3RBIN: " << h3RBIN << std::endl;

  h3TrackResolution->GetXaxis()->SetRange(h3LBIN, h3RBIN);
  auto h2TrackResolution = (TH2 *)h3TrackResolution->Project3D("zy");

  TFile outFile(Form("plots/TrackResolutiont_%s_%s.root", histName,
                           speciesNames[particleId]), "RECREATE");
  h2TrackResolution->Write("h2TrackResolution");
  outFile.Close();


  for (Int_t i = 0; i < pTClasses.size(); i++) {
    TH1D *trackResolution = (TH1D *)h2TrackResolution->ProjectionY(
        Form("TrackResY%i_%s_%s", i + 1, histName, speciesNames[particleId]),
        h2TrackResolution->GetXaxis()->FindBin(pTClasses[i].first),
        h2TrackResolution->GetXaxis()->FindBin(pTClasses[i].second));
    trackResolution->GetXaxis()->SetTitle(TrackResolutionTitleX);
    trackResolution->GetYaxis()->SetTitle(TrackResolutionTitleY);
    legend->AddEntry(trackResolution,
                     Form("%.2f< #it{p}_{T, track}^{gen} <%.2f GeV/#it{c}",
                          pTClasses[i].first, pTClasses[i].second),
                     "pe2");
    // legend2->AddEntry(trackResolution, Form(“Mean: %.2f, Std Dev: %.2f”,
    // trackResolution->GetMean(), trackResolution->GetRMS()), “pe2”);
    // hset(*trackResolution, TrackResolutionTitleX, TrackPtTitleY, 0.7, 1.0,
    // 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*trackResolution, 2., ColorPalette[i], -1., 1., 1e-5, 1e0);
    trackResolution->Draw("esame");
  }
  // }

  return 0;
}

// operate fns
void DrawHistos(const std::vector<TString> &fileNames,
                const std::vector<TString> &histNames,
                const std::vector<int> &ColorPalette) {
  for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
    // gROOT->SetOptStat(0);
    for (Int_t species = PION; species <= ELECTRON; species++) {
    // for (Int_t species = POSITIVE; species <= NEGATIVE; species++) {

      TCanvas *canTrackResolution =
          new TCanvas(Form("TrackResolution_%s_%s", histNames[i].Data(),
                           speciesNames[species]),
                      Form("trackresolution%s_%s", histNames[i].Data(),
                           speciesNames[species]),
                      800, 600);
      canTrackResolution->SetLogy(1);
      canTrackResolution->SetGridx(1);
      canTrackResolution->SetGridy(1);
      canTrackResolution->Draw();
      TLegend *legtrackresolution =
          new TLegend(0.11, 0.5, 0.33, 0.89, NULL, "brNDC");
      legtrackresolution->SetTextSize(0.03);
      legtrackresolution->SetBorderSize(0);
      // TLegend *legtrackresolutionStat = new
      // TLegend(0.302005,0.113043,0.616541,0.330435,NULL,"brNDC");
      TLegend *legtrackresolutionStat =
          new TLegend(0.0, 0.0, 0.0, 0.0, NULL, "brNDC");
      legtrackresolutionStat->SetTextSize(0.03);
      legtrackresolutionStat->SetBorderSize(0);
      TString filePath = mainDir + fileNames[i];
      DrawTrackResolution(filePath.Data(), histNames[i].Data(), species,
                          legtrackresolution, legtrackresolutionStat,
                          ColorPalette[i]);
      legtrackresolution->Draw();
      legtrackresolutionStat->Draw();
      canTrackResolution->Print(Form("plots/TrackResolution_%s_%s.pdf",
                                     histNames[i].Data(),
                                     speciesNames[species]));
    }
  }
}