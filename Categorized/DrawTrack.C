#include "common.h"
#include "DrawJetsMC.h"

Int_t n = 0;
Int_t nn = 0;

TH1 *DrawTrackPt(const char *fileName, const char *histName, float Nevts,
                 TLegend *legend, Color_t colorID) {
  auto file = TFile::Open(fileName, "open");
  auto dir = (TDirectory *)file->Get(Dir);
  dir->cd();
  TH1 *TrackPt = (TH1 *)gROOT->FindObject(TrackPtObj);
  TrackPt =
      TrackPt->Rebin(nTrackptbin, Form("TrackPt_%s", histName), Trackptbin);
  legend->AddEntry(TrackPt, histName);
  hset(*TrackPt, TrackPtTitleX, TrackPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 510);
  hoptset(*TrackPt, Nevts, colorID, 0, 100, 1e-9, 1e0 + 0.05);
  TrackPt->Draw("esame");

  return TrackPt;
}

TH1 *DrawTrackEta(const char *fileName, const char *histName, float Nevts,
                  TLegend *legend, Color_t colorID) {
  auto file = TFile::Open(fileName, "open");
  auto dir = (TDirectory *)file->Get(Dir);
  dir->cd();
  TH1 *TrackEta = (TH1 *)gROOT->FindObject(TrackEtaObj);
  legend->AddEntry(TrackEta, histName);
  hset(*TrackEta, TrackEtaTitleX, TrackEtaTitleY, 0.9, 1.4, 0.05, 0.05, 0.01,
       0.01, 0.05, 0.05, 510, 510);
  hoptset(*TrackEta, Nevts, colorID, -0.9, 0.9, 0.46, 0.64);
  TrackEta->Draw("esame");

  return TrackEta;
}

TH1 *DrawTrackPhi(const char *fileName, const char *histName, float Nevts,
                  TLegend *legend, Color_t colorID) {
  auto file = TFile::Open(fileName, "open");
  auto dir = (TDirectory *)file->Get(Dir);
  dir->cd();
  TH1 *TrackPhi = (TH1 *)gROOT->FindObject(TrackPhiObj);
  std::vector<Double_t> newBins = NewBin(40, -1, 7);
  TrackPhi = TrackPhi->Rebin(40, Form("TrackPhi_%s", histName), newBins.data());
  legend->AddEntry(TrackPhi, histName);
  hset(*TrackPhi, TrackPhiTitleX, TrackPhiTitleY, 0.9, 1.4, 0.05, 0.05, 0.01,
       0.01, 0.05, 0.05, 510, 510);
  hoptset(*TrackPhi, Nevts, colorID, 0, 2 * TMath::Pi(), 0, 0.3);
  TrackPhi->Draw("esame");

  return TrackPhi;
}

void DrawTracks(const std::vector<TString> &fileNames,
                const std::vector<TString> &histNames,
                const std::vector<int> &ColorPallete) {
  if (TrackProcess == 1) {
    // Draw track pT
    Filipad2 *TrackPtPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    TrackPtPad->Draw();
    TPad *trackptpad = TrackPtPad->GetPad(1);
    optFili(*trackptpad, 1, 1, 0, 1);
    TPad *ratiotrackptpad = TrackPtPad->GetPad(2);
    optFili(*ratiotrackptpad, 1, 1, 0, 1);
    TLegend *legtrackpt =
        new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
    legtrackpt->SetTextSize(0.05);
    legtrackpt->SetBorderSize(0);
    trackptpad->cd();
    TH1 *TrackPtRatio =
        DrawTrackPt(refPath.Data(), refName, 1., legtrackpt, kBlack);
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      trackptpad->cd();
      TString filePath = mainDir + fileNames[i];
      TH1 *currentHist = DrawTrackPt(filePath.Data(), histNames[i].Data(), 1.,
                                     legtrackpt, ColorPallete[i]);
      TString ratioName =
          TString::Format("RatioHist_tpt_%s", fileNames[i].Data());
      ratiotrackptpad->cd();
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), TrackPtRatio, currentHist, TrackPtTitleX,
                    RatioTitleY, ColorPallete[i], 1e-1, 1e4);
    }
    trackptpad->cd();
    legtrackpt->Draw();
    TrackPtPad->C->Print(
        Form("plots/TrackPt_ITS%i.pdf", selITS));

    // Draw track eta
    Filipad2 *TrackEtaPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    TrackEtaPad->Draw();
    TPad *tracketapad = TrackEtaPad->GetPad(1);
    optFili(*tracketapad, 1, 1, 0, 0);
    TPad *ratiotracketapad = TrackEtaPad->GetPad(2);
    optFili(*ratiotracketapad, 1, 1, 0, 0);
    TLegend *legtracketa =
        new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
    legtracketa->SetTextSize(0.05);
    legtracketa->SetBorderSize(0);
    tracketapad->cd();
    TH1 *TrackEtaRatio =
        DrawTrackEta(refPath.Data(), refName, 1., legtracketa, kBlack);
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      tracketapad->cd();
      TString filePath = mainDir + fileNames[i];
      TH1 *currentHist = DrawTrackEta(filePath.Data(), histNames[i].Data(), 1.,
                                      legtracketa, ColorPallete[i]);
      TString ratioName =
          TString::Format("RatioHist_teta_%s", fileNames[i].Data());
      ratiotracketapad->cd();
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), TrackEtaRatio, currentHist,
                    TrackEtaTitleX, RatioTitleY, ColorPallete[i], 0.85, 1.15);
    }
    tracketapad->cd();
    legtracketa->Draw();
    TrackEtaPad->C->Print(
        Form("plots/TrackEta_ITS%i.pdf", selITS));

    // Draw track phi
    Filipad2 *TrackPhiPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    TrackPhiPad->Draw();
    TPad *trackphipad = TrackPhiPad->GetPad(1);
    optFili(*trackphipad, 1, 1, 0, 0);
    TPad *ratiotrackphipad = TrackPhiPad->GetPad(2);
    optFili(*ratiotrackphipad, 1, 1, 0, 0);
    TLegend *legtrackphi =
        new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
    legtrackphi->SetTextSize(0.05);
    legtrackphi->SetBorderSize(0);
    trackphipad->cd();
    TH1 *TrackPhiRatio =
        DrawTrackPhi(refPath.Data(), refName, 1, legtrackphi, kBlack);
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      trackphipad->cd();
      TString filePath = mainDir + fileNames[i];
      TH1 *currentHist = DrawTrackPhi(filePath.Data(), histNames[i].Data(), 1,
                                      legtrackphi, ColorPallete[i]);
      TString ratioName =
          TString::Format("RatioHist_tphi_%s", fileNames[i].Data());
      ratiotrackphipad->cd();
      float HistMean = currentHist->GetMean(1);
      cout << "(Dataset, Mean Y) : (" << histNames[i].Data() << ", " << HistMean
           << ")" << endl;
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), TrackPhiRatio, currentHist,
                    TrackPhiTitleX, RatioTitleY, ColorPallete[i], 0.5, 1.5);
    }
    trackphipad->cd();
    legtrackphi->Draw();
    TrackPhiPad->C->Print(
        Form("plots/TrackPhi_ITS%i.pdf", selITS));
  }
}

void DrawJetsMC() {
  DrawTracks(fileNames, histNames, ColorPallete);
}