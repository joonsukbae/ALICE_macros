#include "common.h"
#include "DrawJetsMC.h"

Int_t n = 0;
Int_t nn = 0;

TH1 *DrawConstituentPt(const char *fileName, const char *histName, float Nevts,
                       TLegend *legend, Color_t colorID) {
  auto file = TFile::Open(fileName, "open");
  auto dir = (TDirectory *)file->Get(Dir);
  dir->cd();
  TH3D *H3ConstituentPt = (TH3D *)gROOT->FindObject(ConstPtObj);
  TH1 *ConstituentPt =
      H3ConstituentPt->ProjectionZ(Form("ConstituentPt_%s", histName),
                                   H3ConstituentPt->GetXaxis()->FindBin(RBIN),
                                   H3ConstituentPt->GetXaxis()->FindBin(RBIN),
                                   0, H3ConstituentPt->GetNbinsY());
  ConstituentPt =
      ConstituentPt->Rebin(nTrackptbin, "constituentpt", Trackptbin);
  legend->AddEntry(ConstituentPt, histName);
  hset(*ConstituentPt, ConstPtTitleX, ConstPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01,
       0.01, 0.05, 0.05, 510, 510);
  hoptset(*ConstituentPt, Nevts, colorID, 0, 80, 1e-8, 1);
  ConstituentPt->Draw("esame");

  return ConstituentPt;
}

TH1D *DrawConstituentEta(const char *fileName, const char *histName,
                         float Nevts, TLegend *legend, Color_t colorID) {
  auto file = TFile::Open(fileName, "open");
  auto dir = (TDirectory *)file->Get(Dir);
  dir->cd();
  TH3D *H3ConstituentEta = (TH3D *)gROOT->FindObject(ConstEtaObj);
  TH1D *ConstituentEta =
      H3ConstituentEta->ProjectionZ(Form("ConstituentEta_%s", histName),
                                    H3ConstituentEta->GetXaxis()->FindBin(RBIN),
                                    H3ConstituentEta->GetXaxis()->FindBin(RBIN),
                                    H3ConstituentEta->GetYaxis()->FindBin(40),
                                    H3ConstituentEta->GetYaxis()->FindBin(80));
  legend->AddEntry(ConstituentEta, histName);
  hset(*ConstituentEta, ConstEtaTitleX, ConstEtaTitleY, 0.9, 1.4, 0.05, 0.05,
       0.01, 0.01, 0.05, 0.05, 510, 510);
  hoptset(*ConstituentEta, Nevts, colorID, -0.5, 0.5, 0.4, 1.4);
  ConstituentEta->Draw("esame");

  return ConstituentEta;
}

TH1 *DrawConstituentPhi(const char *fileName, const char *histName, float Nevts,
                        TLegend *legend, Color_t colorID) {
  auto file = TFile::Open(fileName, "open");
  auto dir = (TDirectory *)file->Get(Dir);
  dir->cd();
  TH3D *H3ConstituentPhi = (TH3D *)gROOT->FindObject(ConstPhiObj);
  TH1 *ConstituentPhi =
      H3ConstituentPhi->ProjectionZ(Form("ConstituentPhi_%s", histName),
                                    H3ConstituentPhi->GetXaxis()->FindBin(RBIN),
                                    H3ConstituentPhi->GetXaxis()->FindBin(RBIN),
                                    0, H3ConstituentPhi->GetNbinsY());
  ConstituentPhi->Rebin(2);
  legend->AddEntry(ConstituentPhi, histName);
  hset(*ConstituentPhi, ConstPhiTitleX, ConstPhiTitleY, 0.9, 1.4, 0.05, 0.05,
       0.01, 0.01, 0.05, 0.05, 510, 510);
  hoptset(*ConstituentPhi, Nevts, colorID, 0, 2 * TMath::Pi(), 0., .4);
  ConstituentPhi->Draw("esame");

  return ConstituentPhi;
}

void DrawConstituents(const std::vector<TString> &fileNames,
                      const std::vector<TString> &histNames,
                      const std::vector<int> &ColorPallete) {
  if (ConstituentProcess == 1) {
    // Draw constituents pT
    Filipad2 *ConstPtPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    ConstPtPad->Draw();
    TPad *constptpad = ConstPtPad->GetPad(1);
    optFili(*constptpad, 1, 1, 0, 1);
    TPad *ratioconstptpad = ConstPtPad->GetPad(2);
    optFili(*ratioconstptpad, 1, 1, 0, 0);
    TLegend *legconstpt =
        new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
    legconstpt->SetTextSize(0.05);
    legconstpt->SetBorderSize(0);
    constptpad->cd();
    TH1 *ConstPtRatio =
        DrawConstituentPt(refPath.Data(), refName, 1, legconstpt, kBlack);
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      constptpad->cd();
      TString filePath = mainDir + fileNames[i];
      TH1 *currentHist = DrawConstituentPt(filePath.Data(), histNames[i].Data(),
                                           1, legconstpt, ColorPallete[i]);
      TString ratioName =
          TString::Format("RatioHist_cpt_%s", fileNames[i].Data());
      ratioconstptpad->cd();
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), ConstPtRatio, currentHist, ConstPtTitleX,
                    RatioTitleY, ColorPallete[i], 0., 2.1);
    }
    constptpad->cd();
    legconstpt->Draw();
    ConstPtPad->C->Print(
        Form("plots/ConstPt_R%.1f_ITS%i.pdf", RBIN, selITS));

    // Draw Constituents Eta
    Filipad2 *ConstEtaPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    ConstEtaPad->Draw();
    TPad *constetapad = ConstEtaPad->GetPad(1);
    optFili(*constetapad, 1, 1, 0, 0);
    TPad *ratioconstetapad = ConstEtaPad->GetPad(2);
    optFili(*ratioconstetapad, 1, 1, 0, 0);
    TLegend *legconsteta =
        new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
    TLegend *legconstetaPtRange =
        new TLegend(0.327751, 0.0521739, 0.598086, 0.121739, NULL, "brNDC");
    legconstetaPtRange->SetTextSize(0.05);
    legconstetaPtRange->SetBorderSize(0);
    legconstetaPtRange->AddEntry("", "40 < #it{p}_{T, jet} < 80 GeV", "");
    legconsteta->SetTextSize(0.05);
    legconsteta->SetBorderSize(0);
    constetapad->cd();
    TH1 *ConstEtaRatio =
        DrawConstituentEta(refPath.Data(), refName, 1, legconsteta, kBlack);
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      constetapad->cd();
      TString filePath = mainDir + fileNames[i];
      TH1 *currentHist =
          DrawConstituentEta(filePath.Data(), histNames[i].Data(), 1,
                             legconsteta, ColorPallete[i]);
      TString ratioName =
          TString::Format("RatioHist_ceta_%s", fileNames[i].Data());
      ratioconstetapad->cd();
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), ConstEtaRatio, currentHist,
                    ConstEtaTitleX, RatioTitleY, ColorPallete[i], 0.4, 1.6);
    }
    constetapad->cd();
    legconsteta->Draw();
    legconstetaPtRange->Draw();
    ConstEtaPad->C->Print(
        Form("plots/ConstEta_R%.1f_ITS%i.pdf", RBIN, selITS));

    // Draw Constituents Phi
    Filipad2 *ConstPhiPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    ConstPhiPad->Draw();
    TPad *constphipad = ConstPhiPad->GetPad(1);
    optFili(*constphipad, 1, 1, 0, 0);
    TPad *ratioconstphipad = ConstPhiPad->GetPad(2);
    optFili(*ratioconstphipad, 1, 1, 0, 0);
    TLegend *legconstphi =
        new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
    TLegend *legconstphiPtRange =
        new TLegend(0.327751, 0.0521739, 0.598086, 0.121739, NULL, "brNDC");
    legconstphiPtRange->SetTextSize(0.05);
    legconstphiPtRange->SetBorderSize(0);
    legconstphiPtRange->AddEntry("", "40 < #it{p}_{T, jet} < 80 GeV", "");
    legconstphi->SetTextSize(0.05);
    legconstphi->SetBorderSize(0);
    constphipad->cd();
    TH1 *ConstPhiRatio =
        DrawConstituentPhi(refPath.Data(), refName, 1, legconstphi, kBlack);
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      constphipad->cd();
      TString filePath = mainDir + fileNames[i];
      TH1 *currentHist =
          DrawConstituentPhi(filePath.Data(), histNames[i].Data(), 1,
                             legconstphi, ColorPallete[i]);
      TString ratioName =
          TString::Format("RatioHist_cphi_%s", fileNames[i].Data());
      ratioconstphipad->cd();
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), ConstPhiRatio, currentHist,
                    ConstPhiTitleX, RatioTitleY, ColorPallete[i], 0.2, 1.4);
    }
    constphipad->cd();
    legconstphi->Draw();
    legconstphiPtRange->Draw();
    ConstPhiPad->C->Print(
        Form("plots/ConstPhi_R%.1f_ITS%i.pdf", RBIN, selITS));
  }
}

void DrawJetsMC() {
  DrawConstituents(fileNames, histNames, ColorPallete);
}
 