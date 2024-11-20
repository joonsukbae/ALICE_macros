#ifndef DRAWJETSMC_H
#define DRAWJETSMC_H

#include "BSHelper.cxx"
#include "Filipad2.h"
#include <__config>
#include <iostream>
#include <vector>
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

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

// Function declarations
void DrawJetMatching(const char *fileName, const char *DfileName,
                      const char *histName, float Nevts, TPad *jrep,
                      TPad *ratiojrep, TLegend *JRElegend, TPad *consistencyp,
                      TPad *ratioconsistp, TLegend *CONSISTlegend,
                      TPad *unfoldp, TPad *ratunfoldp, TLegend *unfoldlegend,
                      Color_t colorID, TPad *jrpp = nullptr,
                      TPad *ratiojrpp = nullptr, TLegend *JRPlegend = nullptr,
                      TFile *savefile = nullptr);
void DrawJetMatchings(const std::vector<TString> &fileNames,
                      const std::vector<TString> &histNames,
                      const std::vector<int> &ColorPallete);
void hset(TH1 &hid, TString xtit = "", TString ytit = "", double titoffx = 0.9,
          double titoffy = 1.2, double titsizex = 0.06, double titsizey = 0.06,
          double labeloffx = 0.01, double labeloffy = 0.001,
          double labelsizex = 0.05, double labelsizey = 0.05, int divx = 510,
          int divy = 510);
void hoptset(TH1 &hid, Float_t N = 1, Color_t color = kBlack, Double_t minX = 0,
             Double_t maxX = 100, Double_t minY = 0, Double_t maxY = 1,
             Double_t MarkerSize = 1);
void optFili(TPad &pid, Int_t gridx, Int_t gridy, Int_t logx, Int_t logy);
std::vector<Double_t> NewBin(Int_t Nbins = 100, Double_t minBin = 0,
                             Double_t maxBin = 100);
float Nevents(const char *fileName, const char *eventDir, const char *eventObj);
TH1 *DrawRatio(const char *ratioName, TH1 *refHist, TH1 *testHist, TString AxisTitleX,
               TString AxisTitleY, Color_t colorID, Double_t Ymin, Double_t Ymax,
               Double_t MarkerSize = 1);

#endif // DRAWJETSMC_H