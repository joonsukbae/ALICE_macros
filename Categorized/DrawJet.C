#include "common.h"
#include "DrawJetsMC.h"

Int_t n = 0;
Int_t nn = 0;

TH1 *DrawJetPt(const char *fileName, const char *histName, float Nevts,
               TLegend *legend, Color_t colorID, Int_t i = 0) {
  auto file = TFile::Open(fileName, "open");
  auto dir = (TDirectory *)file->Get(Dir);
  dir->cd();
  TH1 *JetPt = (TH1 *)gROOT->FindObject(JetPtObj);
  JetPt = JetPt->Rebin(nptBins, Form("JetPt_%s", histName), ptbin);
  legend->AddEntry(JetPt, histName);
  hset(*JetPt, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05,
       0.05, 510, 510);
  hoptset(*JetPt, Nevts, colorID, 0, 200, 1e-9, 1e0, 1.2 - 0.2 * i);
  JetPt->Draw("esame");

  return JetPt;
}

TH1 *DrawJetEta(const char *fileName, const char *histName, float Nevts,
                TLegend *legend, Color_t colorID) {
  auto file = TFile::Open(fileName, "open");
  auto dir = (TDirectory *)file->Get(Dir);
  dir->cd();

  TH3D *H3JetEta = (TH3D *)gROOT->FindObject(JetEtaObj);
  TH1 *JetEta = H3JetEta->ProjectionZ(
      Form("JetEta_%s", histName), H3JetEta->GetXaxis()->FindBin(RBIN + 1e-6),
      H3JetEta->GetXaxis()->FindBin(RBIN + 1e-5),
      H3JetEta->GetYaxis()->FindBin(40), H3JetEta->GetYaxis()->FindBin(80));

  legend->AddEntry(JetEta, histName);
  hset(*JetEta, JetEtaTitleX, JetEtaTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 510);
  hoptset(*JetEta, Nevts, colorID, -0.9 + RBIN, 0.9 - RBIN, 0., 2.);
  JetEta->Draw("esame");

  return JetEta;
}

TH1 *DrawJetPhi(const char *fileName, const char *histName, float Nevts,
                TLegend *legend, Color_t colorID) {
  auto file = TFile::Open(fileName, "open");
  auto dir = (TDirectory *)file->Get(Dir);
  dir->cd();

  TH3D *H3JetPhi = (TH3D *)gROOT->FindObject(JetPhiObj);
  TH1D *JetPhi = H3JetPhi->ProjectionZ(
      Form("JetPhi_%s", histName), H3JetPhi->GetXaxis()->FindBin(RBIN + 1e-6),
      H3JetPhi->GetXaxis()->FindBin(RBIN + 1e-5),
      H3JetPhi->GetYaxis()->FindBin(40), H3JetPhi->GetYaxis()->FindBin(80));
  if (JetPhi->GetBinWidth(0) != 0.1) {
    JetPhi->Rebin(4);
  } else {
    JetPhi->Rebin(2);
  }

  legend->AddEntry(JetPhi, histName);
  hset(*JetPhi, JetPhiTitleX, JetPhiTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 510);
  hoptset(*JetPhi, Nevts, colorID, 0, 2 * TMath::Pi(), 0., .4);
  JetPhi->Draw("esame");

  return JetPhi;
}

void DrawJets(const std::vector<TString> &fileNames,
              const std::vector<TString> &histNames,
              const std::vector<int> &ColorPallete) {
  if (JetProcess == 1) {
    // Draw jet pT MCD
    Filipad2 *JetPtPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    JetPtPad->Draw();
    TPad *jetptpad = JetPtPad->GetPad(1);
    optFili(*jetptpad, 1, 1, 0, 1);
    TPad *ratiojetptpad = JetPtPad->GetPad(2);
    optFili(*ratiojetptpad, 1, 1, 0, 0);
    TLegend *legjetpt =
        new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
    legjetpt->SetTextSize(0.05);
    legjetpt->SetBorderSize(0);
    jetptpad->cd();
    TH1 *JetPtRatio = DrawJetPt(refPath.Data(), refName, 1, legjetpt, kBlack);
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      jetptpad->cd();
      TString filePath = mainDir + fileNames[i];
      TH1 *currentHist = DrawJetPt(filePath.Data(), histNames[i].Data(), 1,
                                   legjetpt, ColorPallete[i], i);
      TString ratioName =
          TString::Format("RatioHist_jpt_%s", fileNames[i].Data());
      ratiojetptpad->cd();
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), JetPtRatio, currentHist, JetPtTitleX,
                    RatioTitleY, ColorPallete[i], 0., 1.5);
    }
    jetptpad->cd();
    legjetpt->Draw();
    JetPtPad->C->Print(
        Form("plots/JetPt_R%.1f_ITS%i.pdf", RBIN, selITS));

    // Draw jet eta
    Filipad2 *JetEtaPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    JetEtaPad->Draw();
    TPad *jetetapad = JetEtaPad->GetPad(1);
    optFili(*jetetapad, 1, 1, 0, 0);
    TPad *ratiojetetapad = JetEtaPad->GetPad(2);
    optFili(*ratiojetetapad, 1, 1, 0, 0);
    TLegend *legjeteta =
        new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
    TLegend *legjetetaPtRange =
        new TLegend(0.327751, 0.0521739, 0.598086, 0.121739, NULL, "brNDC");
    legjetetaPtRange->SetTextSize(0.05);
    legjetetaPtRange->SetBorderSize(0);
    legjetetaPtRange->AddEntry("", "40 < #it{p}_{T, jet} < 80 GeV", "");
    legjeteta->SetTextSize(0.05);
    legjeteta->SetBorderSize(0);
    jetetapad->cd();
    TH1 *JetEtaRatio =
        DrawJetEta(refPath.Data(), refName, 1., legjeteta, kBlack);
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      jetetapad->cd();
      TString filePath = mainDir + fileNames[i];
      TH1 *currentHist = DrawJetEta(filePath.Data(), histNames[i].Data(), 1.,
                                    legjeteta, ColorPallete[i]);
      TString ratioName =
          TString::Format("RatioHist_jeta_%s", fileNames[i].Data());
      ratiojetetapad->cd();
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), JetEtaRatio, currentHist, JetEtaTitleX,
                    RatioTitleY, ColorPallete[i], 0., 2.);
    }
    jetetapad->cd();
    legjeteta->Draw();
    legjetetaPtRange->Draw();
    JetEtaPad->C->Print(
        Form("plots/JetEta_R%.1f_ITS%i.pdf", RBIN, selITS));

    // Draw jet phi
    Filipad2 *JetPhiPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    JetPhiPad->Draw();
    TPad *jetphipad = JetPhiPad->GetPad(1);
    optFili(*jetphipad, 1, 1, 0, 0);
    TPad *ratiojetphipad = JetPhiPad->GetPad(2);
    optFili(*ratiojetphipad, 1, 1, 0, 0);
    TLegend *legjetphi =
        new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
    TLegend *legjetphiPtRange =
        new TLegend(0.327751, 0.0521739, 0.598086, 0.121739, NULL, "brNDC");
    legjetphiPtRange->SetTextSize(0.05);
    legjetphiPtRange->SetBorderSize(0);
    legjetphiPtRange->AddEntry("", "40 < #it{p}_{T, jet} < 80 GeV", "");
    legjetphi->SetTextSize(0.05);
    legjetphi->SetBorderSize(0);
    jetphipad->cd();
    TH1 *JetPhiRatio =
        DrawJetPhi(refPath.Data(), refName, 1., legjetphi, kBlack);
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      jetphipad->cd();
      TString filePath = mainDir + fileNames[i];
      TH1 *currentHist = DrawJetPhi(filePath.Data(), histNames[i].Data(), 1.,
                                    legjetphi, ColorPallete[i]);
      TString ratioName =
          TString::Format("RatioHist_jphi_%s", fileNames[i].Data());
      ratiojetphipad->cd();
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), JetPhiRatio, currentHist, JetPhiTitleX,
                    RatioTitleY, ColorPallete[i], 0., 2.);
    }
    jetphipad->cd();
    legjetphi->Draw();
    legjetphiPtRange->Draw();
    JetPhiPad->C->Print(
        Form("plots/JetPhi_R%.1f_ITS%i.pdf", RBIN, selITS));
  }
}

void DrawJetsMC() {
  DrawJets(fileNames, histNames, ColorPallete);
}