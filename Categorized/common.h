#ifndef COMMON_H
#define COMMON_H

#include "BSHelper.cxx"
#include "Filipad2.h"
#include <__config>
#include <iostream>
#include <vector>
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
using namespace std;

enum {
  kJSbegin = 0,
  kMCD,
  kMCDMATCHED,
  kFAKE,
  kMCP,
  kMCPMATCHED,
  kMISS,
  kJSend
};

extern Int_t n;
extern Int_t nn;

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
             Double_t MarkerSize = 1) {
  if (N == 1) {
    hid.Scale(1. / hid.Integral(), "width");
  } else if (N == 2) {
    hid.Scale(1. / hid.Integral(), "");
  } else {
    hid.Scale(1. / N, "width");
  }

  hid.SetMarkerColor(color);
  hid.SetLineColor(color);
  hid.SetMarkerSize(MarkerSize);
  hid.SetMarkerStyle(22);

  hid.GetXaxis()->SetRangeUser(minX, maxX);
  hid.GetYaxis()->SetRangeUser(minY, maxY);
  hid.SetFillColorAlpha(color, 0.3);
  hid.GetYaxis()->SetNdivisions(505);
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

extern Double_t Trackptbin[15];
extern Double_t ptbin[21];
extern Double_t ptbinHEP[20];
extern Int_t nTrackptbin;
extern Int_t nptBins;
extern Int_t nptBinsHEP;
extern Double_t RBIN;
extern Int_t selITS;
extern Int_t Npthat;
extern const char *Dir;
extern const char *EventObj;
extern const char *EventWObj;
extern TString mainDir;
extern TString refDir;
extern TString refFile;
extern TString refPath;
extern const char *refName;
extern std::vector<TString> fileNames;
extern std::vector<TString> histNames;
extern std::vector<int> ColorPallete;
extern const char *TrackPtObj;
extern const char *TrackEtaObj;
extern const char *TrackPhiObj;
extern const char *ConstPtObj;
extern const char *ConstEtaObj;
extern const char *ConstPhiObj;
extern const char *JetPtObj;
extern const char *JetPtMCPObj;
extern const char *JetEtaObj;
extern const char *JetPhiObj;
extern const char *JetNtracksObj;
extern const char *JetAreaObj;
extern const char *JetResolutionObj;
extern const char *JetResolutionObjOld;
extern const char *RatioTitleY;
extern const char *TrackPtTitleX;
extern const char *TrackPtTitleY;
extern const char *TrackEtaTitleX;
extern const char *TrackEtaTitleY;
extern const char *TrackPhiTitleX;
extern const char *TrackPhiTitleY;
extern const char *ConstPtTitleX;
extern const char *ConstPtTitleY;
extern const char *ConstEtaTitleX;
extern const char *ConstEtaTitleY;
extern const char *ConstPhiTitleX;
extern const char *ConstPhiTitleY;
extern const char *JetPtTitleX;
extern const char *JetPtTitleY;
extern const char *JetPtGenTitleX;
extern const char *JetEtaTitleX;
extern const char *JetEtaTitleY;
extern const char *JetPhiTitleX;
extern const char *JetPhiTitleY;
extern const char *JetNtracksTitleX;
extern const char *JetNtracksTitleY;
extern const char *JRETitleX;
extern const char *JRETitleY;
extern const char *JRPTitleX;
extern const char *JRPTitleY;
extern const char *JetResolutionTitleX;
extern const char *JetResolutionTitleY;

#endif // COMMON_H