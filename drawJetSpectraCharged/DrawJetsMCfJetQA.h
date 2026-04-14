///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Draw macro for jet QA         //////////
////////// author: Joonsuk Bae           ////////// 
////////// E-mail: jbae@cern.ch          //////////
////////// Last Modified: 19 July 2024   //////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

#include "TSVDUnfold_local.h"
#include "RooUnfoldSvd.h"

// d-vector k optimization (compact version for QA macro)
// Full implementation with bilateral stability: see DrawJetsMCRDependentHelpers.h
inline Int_t FindOptimalSvdK_Dvector(RooUnfoldResponse *Response, TH1 *hRecoData,
                                      Double_t /*stabilityThreshold*/ = 0.03) {
  if (!Response || !hRecoData) return 6;
  const Int_t nBins = hRecoData->GetNbinsX();
  const Int_t maxK = TMath::Min(nBins, 25);

  RooUnfoldSvd unfoldFull(Response, hRecoData, nBins);
  unfoldFull.Hreco();
  TSVDUnfold_local *svdImpl = unfoldFull.Impl();
  if (!svdImpl) return 6;
  TH1D *hD = svdImpl->GetD();
  if (!hD) return 6;

  const Int_t nD = hD->GetNbinsX();
  std::cout << "[Info] d-vector: ";
  for (Int_t i = 1; i <= nD; ++i) {
    Double_t absDi = TMath::Abs(hD->GetBinContent(i));
    std::cout << "|d_" << i << "|=" << Form("%.1f", absDi) << " ";
  }
  std::cout << std::endl;

  // Robust boundary: require 3+ out of 4 consecutive values below threshold
  // (avoids premature capping from isolated noise pockets, e.g. R=0.2 UE-sub)
  Int_t k = maxK;
  const Int_t windowSz = 4;
  const Int_t noiseReq = 3;
  for (Int_t i = 1; i <= nD - windowSz + 1; ++i) {
    Int_t nNoise = 0;
    for (Int_t j = 0; j < windowSz; ++j) {
      if (TMath::Abs(hD->GetBinContent(i + j)) < 1.0)
        nNoise++;
    }
    if (nNoise >= noiseReq) {
      k = TMath::Max(2, i - 1);
      break;
    }
  }
  if (k > 25) k = 25;
  std::cout << "[Info] Optimal SVD k (d-vector): " << k << std::endl;
  return k;
}

// Draw Jets
#include <cstdlib>
#include <cctype>
TH1 *DrawJetPt(const char *fileName, const char *histName, const char *Obj, Double_t Nevts,
               TLegend *legend, Color_t colorID, Int_t i = 0,
               const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  TH3 *h3JetPt = (TH3 *)file->Get(Form("%s/%s", Dir, Obj));
  TH1 *JetPt = h3JetPt->ProjectionX(Form("JetPt_%s", histName), 1, h3JetPt->GetNbinsY(), 1, h3JetPt->GetNbinsZ());
  if (REBINON) {
    JetPt = JetPt->Rebin(nptBins, Form("JetPt_%s", histName), ptbin);
  }
  legend->AddEntry(JetPt, histName);
  hset(*JetPt, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05,
       0.05, 510, 510);
  auto [yMin, yMax] = getYAxisRange(JetPt, Nevts, PlotPtMin, PlotPtMax); 
  hoptset(*JetPt, Nevts, colorID, PlotPtMin, PlotPtMax, yMin, yMax, colorID==kBlack? 1.2 : 1, 1, 2, colorID==kBlack? 21 : 20);
  JetPt->Draw("esame");
  // JetPt->SaveAs(Form("JetPt_%s.root", histName));

  return JetPt;
}
TH1 *DrawJetPtMCP(const char *fileName, const char *histName, const char *Obj, Double_t Nevts,
                  TLegend *legend, Color_t colorID, Int_t i = 0,
                  const char *Dir = nullptr) {
  TFile *file = TFile::Open(fileName, "open");
  if (!file || file->IsZombie()) {
    printf("[Error] Could not open file: %s\n", fileName);
    if (file) file->Close();
    return nullptr;
  }
  
  TH1 *JetPtMCP = nullptr;
  TString histPath = Form("%s/%s", Dir, Obj);
  
  // Try to load as TH1 first (1D histogram)
  TH1 *tempHist = (TH1 *)file->Get(histPath.Data());
  if (tempHist) {
    // Create independent histogram copy and detach from file
    JetPtMCP = (TH1 *)tempHist->Clone(Form("JetPtMCP_%s", histName));
    JetPtMCP->SetDirectory(0);
  } else {
    // Try to load as TH3 (3D histogram)
    TH3 *h3JetPtMCP = (TH3 *)file->Get(histPath.Data());
    if (h3JetPtMCP) {
      JetPtMCP = h3JetPtMCP->ProjectionX(Form("JetPtMCP_%s", histName), 1, h3JetPtMCP->GetNbinsY(), 1, h3JetPtMCP->GetNbinsZ());
    } else {
      printf("[Error] Could not find histogram: %s\n", histPath.Data());
      file->Close();
      return nullptr;
    }
  }
  
  file->Close();
  
  if (!JetPtMCP) {
    printf("[Error] Failed to load histogram\n");
    return nullptr;
  }
  
  if (REBINON) {
    TH1 *rebinnedHist = JetPtMCP->Rebin(nptBinsGen,Form("JetPtMCP_%s_rebinned",histName),ptbinGen);
    if (rebinnedHist) {
      delete JetPtMCP;
      JetPtMCP = rebinnedHist;
    }
  }
  
  legend->AddEntry(JetPtMCP, histName);
  hset(*JetPtMCP, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 510); 
  hoptset(*JetPtMCP, Nevts, colorID, 0, PlotPtMax, 1e-8, 1e-0, 1.2 - 0.2 * i);
  JetPtMCP->Draw("esame");

  return JetPtMCP;
}
TH1 *DrawJetEta(const char *fileName, const char *histName, const char *Obj, Double_t Nevts,
                TLegend *legend, Color_t colorID, const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  TH3D *H3JetEta = (TH3D *)file->Get(Form("%s/%s", Dir, Obj));
  TH1 *JetEta = H3JetEta->ProjectionY(Form("JetEta_%s", histName), 1, H3JetEta->GetNbinsY(), 1, H3JetEta->GetNbinsZ());

  // Rebin by 5: ensures bin edges at multiples of 0.1, aligning with fiducial |eta|<0.9-R
  TH1 *h1JetEta = JetEta->Rebin(5, Form("JetEta_rebinned_%s", histName));

  // TH1D *h1JetEta = new TH1D("ConstEta", "ConstEta", 90, -0.9, 0.9);
  // for (int ieta = 1; ieta <= 50; ieta++) {
  //   double etaContent = JetEta->GetBinContent(JetEta->GetXaxis()->FindBin(-0.51 + 0.02 * ieta));
  //   double etaError = JetEta->GetBinError(JetEta->GetXaxis()->FindBin(-0.51 + 0.02 * ieta));
  //   h1JetEta->SetBinContent(ieta+20, etaContent);
  //   h1JetEta->SetBinError(ieta+20, etaError);
  // }

  legend->AddEntry(h1JetEta, histName);
  hset(*h1JetEta, JetEtaTitleX, JetEtaTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 510);
  // NORMEVENTS? hoptset(*h1JetEta, Nevts, colorID, -0.9 + RBIN, 0.9 - RBIN, 0., 0.1) : hoptset(*h1JetEta, Nevts, colorID, -0.9 + RBIN, 0.9 - RBIN, 0., 1.5);
  NORMEVENTS? hoptset(*h1JetEta, Nevts, colorID, -1.0, 1.0, 0., 0.1) : hoptset(*h1JetEta, Nevts, colorID, -1.0, 1.0, 0., 1.5);
  h1JetEta->Draw("esame");

  return h1JetEta;
}
TH1 *DrawJetPhi(const char *fileName, const char *histName, const char *Obj, Double_t Nevts,
                TLegend *legend, Color_t colorID, const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  TH3D *H3JetPhi = (TH3D *)file->Get(Form("%s/%s", Dir, Obj));
  TH1D *JetPhi = H3JetPhi->ProjectionZ(Form("JetPhi_%s", histName), 1, H3JetPhi->GetNbinsY(), 1, H3JetPhi->GetNbinsZ());
  // Rebin by 4: finer bins (0.2 width) to reduce edge effects at phi~2pi boundary
  TH1 *h1JetPhi = JetPhi->Rebin(4, Form("JetPhi_rebinned_%s", histName));
  legend->AddEntry(h1JetPhi, histName, "p");
  hset(*h1JetPhi, JetPhiTitleX, JetPhiTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 510);
  // Display range clipped to avoid partially populated edge bins (phiAxis = {160, -1, 7})
  NORMEVENTS? hoptset(*h1JetPhi, Nevts, colorID, 0.01, 6.19, 0., 2e-2) : hoptset(*h1JetPhi, Nevts, colorID, 0.01, 6.19, 0., 0.3);
  h1JetPhi->SetMarkerStyle(20);
  h1JetPhi->SetMarkerColor(colorID);
  h1JetPhi->SetMarkerSize(0.8);
  h1JetPhi->Draw("esame");

  return h1JetPhi;
}
TH1 *DrawJetNtracks(const char *fileName, const char *histName, const char *Obj, Double_t Nevts,
                    TLegend *legend, Color_t colorID,
                    const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  TH2 *h2JetNtracks = (TH2 *)file->Get(Form("%s/%s", Dir, Obj));
  TH1 *JetNtracks = h2JetNtracks->ProjectionY(Form("JetNtracks_%s", histName), 1, h2JetNtracks->GetNbinsX());
  // if (REBINON) {
  //   JetNtracks->Rebin(2);
  // }
  legend->AddEntry(JetNtracks, histName);
  hset(*JetNtracks, JetNtracksTitleX, JetNtracksTitleY, 0.9, 1.4, 0.05, 0.05,
       0.01, 0.01, 0.05, 0.05, 510, 1005);
  NORMEVENTS? hoptset(*JetNtracks, Nevts, colorID, 0, 25, 5e-10, 1e-1) : hoptset(*JetNtracks, Nevts, colorID, 0, 40, 1e-8, 1e0);
  JetNtracks->Draw("esame");

  return JetNtracks;
}
TH2 *DrawJetArea(const char *fileName, const char *histName, const char *Obj, Double_t Nevts,
                 TLegend *legend, Color_t colorID, const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file: " << fileName << std::endl;
    return nullptr;
  }

  TH2 *JetArea = (TH2 *)file->Get(Form("%s/%s", Dir, Obj));
  if (!JetArea) {
    std::cerr << "Error projecting histogram: " << histName << std::endl;
    return nullptr;
  }

  legend->AddEntry(JetArea, histName, "l");
  hset(*JetArea, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 510);
  JetArea->Scale(1. / JetArea->Integral(), "width");
  JetArea->GetXaxis()->SetRangeUser(0., PlotPtMax);
  JetArea->GetYaxis()->SetRangeUser(0., 2.);
  JetArea->GetZaxis()->SetRangeUser(1e-7, 1e0);
  JetArea->Draw("colz");

  return JetArea;
}

// Draw Jet Matchings
// std::vector<TH1F*> UnfoldedHistos;
TH1 *DrawJetMatching(const char *MCfileName, const char *DfileName,
                     const char *histName, Double_t NevtsData, Double_t NevtsDataUnTrig, Double_t NevtsMCD, Double_t NevtsMCP, TPad *jrep,
                     TPad *ratiojrep, TLegend *JRElegend, TPad *consistencyp,
                     TPad *ratioconsistp, TLegend *CONSISTlegend, TPad *unfoldp,
                     TPad *ratunfoldp, TLegend *unfoldlegend, Color_t colorID,
                     TPad *jrpp = nullptr, TPad *ratiojrpp = nullptr,
                     TLegend *JRPlegend = nullptr, TFile *savefile = nullptr, const char *RefDir = nullptr,
                     const char *Dir = nullptr, const char *mccollcounterfile = nullptr,
                     int svdKReg = 4) {

  // Skim pT cut parsing removed to restore original behavior

  auto Dfile = TFile::Open(DfileName, "open");
  TH1 *DJetpt =
      (TH1 *)Dfile->Get(Form("%s/%s", RefDir, "h_jet_pt"));
  if (REBINON) {
    DJetpt = DJetpt->Rebin(nptBins, Form("DJetpt_%s", histName), ptbin);
  }
  double NjetsData = DJetpt->Integral(DJetpt->FindBin(5), DJetpt->GetNbinsX());
  cout << "The number of jets in measured Data: " << NjetsData << endl;

  auto MCfile = TFile::Open(MCfileName, "open");

  TH1 *JetMCPPtINEL = (TH1 *)MCfile->Get(Form("%s/%s", Dir, "h_jet_pt_part"));
  if (REBINON) {
    JetMCPPtINEL = JetMCPPtINEL->Rebin(nptBinsGen, Form("JetMCPPtINEL%s", histName), ptbinGen); 
  }
  double NjetsMCPINEL = JetMCPPtINEL->Integral(JetMCPPtINEL->FindBin(5), JetMCPPtINEL->GetNbinsX());
  cout << "The number of jets in true MC: " << NjetsMCPINEL << endl;

  TH1 *JetMCPPt = (TH1 *)MCfile->Get(Form("%s/%s", Dir, "h_jet_pt_part"));
  if (REBINON) {
    JetMCPPt = JetMCPPt->Rebin(nptBinsGen, Form("JetMCPPt_%s", histName), ptbinGen); 
  }
  double NjetsMCP = JetMCPPt->Integral(JetMCPPt->FindBin(5), JetMCPPt->GetNbinsX());
  cout << "The number of jets in true MC after TVX: " << NjetsMCP << endl;

  TH1 *JetMCDPt = (TH1 *)MCfile->Get(Form("%s/%s", Dir, "h_jet_pt"));
  if (REBINON) {
    JetMCDPt = JetMCDPt->Rebin(nptBins, Form("JetMCDPt_%s", histName), ptbin); 
  }
  double NjetsMCD = JetMCDPt->Integral(JetMCDPt->FindBin(5), JetMCDPt->GetNbinsX());
  cout << "The number of jets in reco MC: " << NjetsMCD << endl;

  TString histNameStr(histName);
  // Load 2D correlate histogram
  TH2 *HCorrelate2D = (TH2 *)MCfile->Get(
      Form("%s/%s", Dir, "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint"));
  if (!HCorrelate2D) {
    std::cerr << "[Error] No 2D correlate histogram found." << std::endl;
    return nullptr;
  }

  TH2F *h2HCorrelate = nullptr;
  if (REBINON) {
    // Custom rebin: 2D histogram (det pt on X, part pt on Y)
    h2HCorrelate = new TH2F(Form("hcorrelate_%s", histName),
                             Form("correlate_%s", histName),
                             nptBins, ptbin, nptBinsGen, ptbinGen);
    for (Int_t ix = 1; ix <= HCorrelate2D->GetNbinsX(); ++ix) {
      for (Int_t iy = 1; iy <= HCorrelate2D->GetNbinsY(); ++iy) {
        Double_t content = HCorrelate2D->GetBinContent(ix, iy);
        Double_t error = HCorrelate2D->GetBinError(ix, iy);
        Int_t xbin = h2HCorrelate->GetXaxis()->FindBin(HCorrelate2D->GetXaxis()->GetBinCenter(ix));
        Int_t ybin = h2HCorrelate->GetYaxis()->FindBin(HCorrelate2D->GetYaxis()->GetBinCenter(iy));
        Double_t currentContent = h2HCorrelate->GetBinContent(xbin, ybin);
        Double_t currentError = h2HCorrelate->GetBinError(xbin, ybin);
        h2HCorrelate->SetBinContent(xbin, ybin, currentContent + content);
        h2HCorrelate->SetBinError(xbin, ybin, std::sqrt(currentError * currentError + error * error));
      }
    }
  } else {
    // Use original binning
    h2HCorrelate = (TH2F *)HCorrelate2D->Clone(Form("hcorrelate_%s", histName));
  }

  // Project to get matched distributions
  TH1 *MCDMatchedpt = (TH1 *)h2HCorrelate->ProjectionX(Form("hMCDMatched_%s", histName), 1, h2HCorrelate->GetNbinsY(), "e");
  TH1 *MCPMatchedpt = (TH1 *)h2HCorrelate->ProjectionY(Form("hMCPMatched_%s", histName), 1, h2HCorrelate->GetNbinsX(), "e");

  // TH1 *MCPMatchedptSel = (TH1F *) HCorrelate->ProjectionY(Form("hMCPMatchedSel_%s", histName), RBIN, RBIN, HCorrelate->GetZaxis()->FindBin(5), HCorrelate->GetZaxis()->FindBin(100), "e");
  // MCPMatchedptSel = MCPMatchedptSel->Rebin(nptBins, Form("MCPptMatchedSel_%s", histName), ptbin); 
  // TH1 *MCPMatchedptSel = (TH1F *) h2HCorrelate->ProjectionY(Form("hMCPMatchedSel_%s", histName), h2HCorrelate->GetXaxis()->FindBin(5), h2HCorrelate->GetXaxis()->FindBin(100), "e"); 
  // MCPMatchedpt = MCPMatchedpt->Rebin(nptBins, Form("MCPptMatched_%s", histName), ptbin);

  TH2F *Respt = (TH2F *)h2HCorrelate->Clone(Form("hist_%i", ++n));
  TH1F *fake = (TH1F *)JetMCDPt->Clone(Form("hist_%i", ++n));
  fake->Add(MCDMatchedpt, -1);
  TH1F *miss = (TH1F *)JetMCPPt->Clone(Form("hist_%i", ++n));
  miss->Add(MCPMatchedpt, -1);

  // Kernel-style unfolding:
  // - Response is built from matched pairs only (no Fake/Miss filled into RooUnfoldResponse)
  // - Fake subtraction (purity) is applied to input BEFORE unfolding
  // - Efficiency correction is applied AFTER unfolding
  // 3-arg constructor: directly sets _res=Respt, _mes=MCDMatchedpt, _tru=MCPMatchedpt
  // No Fill() loop — avoids doubling _mes/_tru that occurs with 2-arg + Fill pattern
  // Preserves original bin errors (Sumw2) from the 2D matrix
  RooUnfoldResponse *Response = new RooUnfoldResponse(MCDMatchedpt, MCPMatchedpt, Respt);
  std::cout << "[Info] RooUnfoldResponse built with 3-arg constructor (matched pairs, raw counts)" << std::endl;

  TH2 *hResponseMatrix = Response->Hresponse();

  // TCanvas *canhResponseMatrix =
  //     new TCanvas(Form("hResponseMatrix_%s", histName),
  //                 Form("hResponseMatrix_%s", histName), 800, 800);
  // canhResponseMatrix->cd();
  // canhResponseMatrix->SetLogz(1);

  // auto hNormResponseMatrix = (TH2 *)hResponseMatrix->Clone(Form("hist_%i", ++n));
  // hNormResponseMatrix->Scale(1. / hNormResponseMatrix->Integral(), "width");
  // hNormResponseMatrix->GetZaxis()->SetRangeUser(1e-10, 1e0);
  // // hResponseMatrix->GetZaxis()->SetRangeUser(1e-3, 1e3);
  // hNormResponseMatrix->Draw("colz");

  std::cout << "[Info] Preparing kernel corrections for unfolding..." << std::endl;

  // Kernel-style corrections (computed from MC)
  TH1* hPurity = (TH1*)JetMCDPt->Clone(Form("hPurity_%s", histName));
  hPurity->SetDirectory(nullptr);
  hPurity->Reset();
  for (int i = 1; i <= hPurity->GetNbinsX(); ++i) {
    const double denom = JetMCDPt->GetBinContent(i);
    const double numer = MCDMatchedpt->GetBinContent(i);
    double p = (denom > 0) ? (numer / denom) : 0.0;
    if (p < 0) p = 0.0;
    if (p > 1) p = 1.0;
    hPurity->SetBinContent(i, p);
  }

  TH1* hEfficiency = (TH1*)JetMCPPt->Clone(Form("hEfficiency_%s", histName));
  hEfficiency->SetDirectory(nullptr);
  hEfficiency->Reset();
  for (int i = 1; i <= hEfficiency->GetNbinsX(); ++i) {
    const double denom = JetMCPPt->GetBinContent(i);
    const double numer = MCPMatchedpt->GetBinContent(i);
    double eff = (denom > 0) ? (numer / denom) : 0.0;
    if (eff < 0) eff = 0.0;
    if (eff > 1) eff = 1.0;
    hEfficiency->SetBinContent(i, eff);
  }

  // Purity correction (fake subtraction) BEFORE unfolding
  TH1* DJetptMatched = (TH1*)DJetpt->Clone(Form("DJetptMatched_%s", histName));
  DJetptMatched->SetDirectory(nullptr);
  for (int i = 1; i <= DJetptMatched->GetNbinsX(); ++i) {
    const double x = DJetpt->GetBinContent(i);
    const double ex = DJetpt->GetBinError(i);
    const double p = hPurity->GetBinContent(i);
    DJetptMatched->SetBinContent(i, x * p);
    DJetptMatched->SetBinError(i, ex * p);
  }

  // MC closure input (matched reco) for unfolding
  TH1* JetMCDPtMatched = (TH1*)JetMCDPt->Clone(Form("JetMCDPtMatched_%s", histName));
  JetMCDPtMatched->SetDirectory(nullptr);
  for (int i = 1; i <= JetMCDPtMatched->GetNbinsX(); ++i) {
    const double x = JetMCDPt->GetBinContent(i);
    const double ex = JetMCDPt->GetBinError(i);
    const double p = hPurity->GetBinContent(i);
    JetMCDPtMatched->SetBinContent(i, x * p);
    JetMCDPtMatched->SetBinError(i, ex * p);
  }

  // Data d-vector determines nominal k (standard ALICE practice)
  // MC closure svdKReg is logged for cross-check only
  Int_t dataK = FindOptimalSvdK_Dvector(Response, DJetptMatched);
  if (svdKReg > 0 && svdKReg != dataK) {
    std::cout << "[Info] MC closure k=" << svdKReg << " vs data d-vector k=" << dataK
              << " (using data, standard practice)" << std::endl;
  }
  std::cout << "[Info] Using data-optimized kreg=" << dataK << " for unfolding" << std::endl;

  RooUnfoldSvd unfoldCon(Response, JetMCDPtMatched, dataK);
  RooUnfoldSvd unfold(Response, DJetptMatched, dataK);
  RooUnfoldSvd unfold2(Response, DJetptMatched, TMath::Max(2, dataK - 1)); // systematic study: kreg - 1
  RooUnfoldSvd unfold3(Response, DJetptMatched, TMath::Min(20, dataK + 1)); // systematic study: kreg + 1

  auto hMCcorrected = (TH1F *)unfoldCon.Hreco();
  auto hDcorrected = (TH1F *)unfold.Hreco();

  cout << "[Kernel-style] Matched reco MC jet yield [10,140]: "
       << JetMCDPtMatched->Integral(JetMCDPtMatched->GetXaxis()->FindBin(PlotPtMin), JetMCDPtMatched->GetXaxis()->FindBin(PlotPtMax))
       << ", unfolded detected-truth yield [10,140]: "
       << hMCcorrected->Integral(hMCcorrected->GetXaxis()->FindBin(PlotPtMin), hMCcorrected->GetXaxis()->FindBin(PlotPtMax))
       << ", matched truth yield [10,140]: "
       << MCPMatchedpt->Integral(MCPMatchedpt->GetXaxis()->FindBin(PlotPtMin), MCPMatchedpt->GetXaxis()->FindBin(PlotPtMax))
       << endl;

  cout << "[Kernel-style] Raw data reco jet yield [10,140]: "
       << DJetpt->Integral(DJetpt->GetXaxis()->FindBin(PlotPtMin), DJetpt->GetXaxis()->FindBin(PlotPtMax))
       << ", purity-corrected (matched) yield [10,140]: "
       << DJetptMatched->Integral(DJetptMatched->GetXaxis()->FindBin(PlotPtMin), DJetptMatched->GetXaxis()->FindBin(PlotPtMax))
       << ", unfolded detected-truth yield [10,140]: "
       << hDcorrected->Integral(hDcorrected->GetXaxis()->FindBin(PlotPtMin), hDcorrected->GetXaxis()->FindBin(PlotPtMax))
       << endl;
  
  // Statistical errors from SVD covariance matrix diagonal
  // (aligned with DrawMcClosureTest.C — clean Ereco() extraction)
  TMatrixD covMatrix = unfold.Ereco();
  int nBins = hDcorrected->GetNbinsX();
  bool useCovMatrix = (covMatrix.GetNrows() >= nBins && covMatrix.GetNcols() >= nBins);
  for (int i = 1; i <= nBins; ++i) {
    int idx = i - 1;
    if (useCovMatrix && idx < covMatrix.GetNrows()) {
      double covDiag = covMatrix(idx, idx);
      if (covDiag > 0) {
        hDcorrected->SetBinError(i, TMath::Sqrt(covDiag));
      }
    }
  }

  // Efficiency correction (miss) AFTER unfolding
  for (int i = 1; i <= hDcorrected->GetNbinsX(); ++i) {
    const double eff = hEfficiency->GetBinContent(i);
    if (eff > 0) {
      hDcorrected->SetBinContent(i, hDcorrected->GetBinContent(i) / eff);
      hDcorrected->SetBinError(i, hDcorrected->GetBinError(i) / eff);
    } else {
      hDcorrected->SetBinContent(i, 0.0);
      hDcorrected->SetBinError(i, 0.0);
    }
  }
  if (hMCcorrected) {
    for (int i = 1; i <= hMCcorrected->GetNbinsX(); ++i) {
      const double eff = hEfficiency->GetBinContent(i);
      if (eff > 0) {
        hMCcorrected->SetBinContent(i, hMCcorrected->GetBinContent(i) / eff);
        hMCcorrected->SetBinError(i, hMCcorrected->GetBinError(i) / eff);
      } else {
        hMCcorrected->SetBinContent(i, 0.0);
        hMCcorrected->SetBinError(i, 0.0);
      }
    }
  }
  
  auto hDcorrected2 = (TH1F *)unfold2.Hreco(); // systematic study
  auto hDcorrected3 = (TH1F *)unfold3.Hreco(); // systematic study
  if (hDcorrected2) {
    for (int i = 1; i <= hDcorrected2->GetNbinsX(); ++i) {
      const double eff = hEfficiency->GetBinContent(i);
      if (eff > 0) {
        hDcorrected2->SetBinContent(i, hDcorrected2->GetBinContent(i) / eff);
        hDcorrected2->SetBinError(i, hDcorrected2->GetBinError(i) / eff);
      } else {
        hDcorrected2->SetBinContent(i, 0.0);
        hDcorrected2->SetBinError(i, 0.0);
      }
    }
  }
  if (hDcorrected3) {
    for (int i = 1; i <= hDcorrected3->GetNbinsX(); ++i) {
      const double eff = hEfficiency->GetBinContent(i);
      if (eff > 0) {
        hDcorrected3->SetBinContent(i, hDcorrected3->GetBinContent(i) / eff);
        hDcorrected3->SetBinError(i, hDcorrected3->GetBinError(i) / eff);
      } else {
        hDcorrected3->SetBinContent(i, 0.0);
        hDcorrected3->SetBinError(i, 0.0);
      }
    }
  }

  // int maxK = 20;
  // OptimizeRegularizationParameter(Response, DJetpt, maxK);
  
  RooUnfoldSvd unfoldSVD(Response, DJetptMatched, 16); // method comparison: fixed high kreg
  TH1* hDcorrectedSVD = unfoldSVD.Hreco();
  if (hDcorrectedSVD) {
    for (int i = 1; i <= hDcorrectedSVD->GetNbinsX(); ++i) {
      const double eff = hEfficiency->GetBinContent(i);
      if (eff > 0) {
        hDcorrectedSVD->SetBinContent(i, hDcorrectedSVD->GetBinContent(i) / eff);
        hDcorrectedSVD->SetBinError(i, hDcorrectedSVD->GetBinError(i) / eff);
      } else {
        hDcorrectedSVD->SetBinContent(i, 0.0);
        hDcorrectedSVD->SetBinError(i, 0.0);
      }
    }
  }

  //////////////////////////
  // Jet Momentum Resolution
  //////////////////////////
  cout << "Starting Jet Momentum Resolution..." << endl;
  
  TCanvas* canJMR = new TCanvas(Form("JMR_%s", histName), Form("Jet momentum resolution - %s", histName), 900, 800);
  setpad(canJMR, 0.02, 0.15, 0.15, 0.05);
  canJMR->cd();
  gPad->SetTicks(1,1);
  TLegend* legJMR = new TLegend(0.300111,0.450323,0.800111,0.553548, NULL, "brNDC");
    legJMR->SetTextSize(0.043);
    legJMR->SetBorderSize(0);
    legJMR->SetFillColorAlpha(0,0); 
  TH1* hjmr = (TH1*) JetMCPPt->Clone(Form("hist_%i", ++n));
  hjmr->Reset();
  Double_t jmrsum = 0;
  Double_t jmrsumweight = 0;
  for (Int_t GenPtBin = h2HCorrelate->GetYaxis()->FindBin(10); GenPtBin <= h2HCorrelate->GetNbinsY(); GenPtBin++) {
  // for (Int_t GenPtBin = 1; GenPtBin <= h2HCorrelate->GetNbinsY(); GenPtBin++) {
      TH1* hRecPt = (TH1*) h2HCorrelate->ProjectionX(Form("hRecPt_%s", histName), GenPtBin, GenPtBin, "e");
      Double_t RecPtRMS = hRecPt->GetRMS();
      Double_t RecPtRMSError = hRecPt->GetRMSError();
      Double_t GenPtCenter = h2HCorrelate->GetYaxis()->GetBinCenter(GenPtBin);
      Double_t GenPtBinWidth = h2HCorrelate->GetXaxis()->GetBinWidth(GenPtBin);

      Double_t JMRValue = RecPtRMS / GenPtCenter;
      Double_t JMRError = RecPtRMSError / GenPtCenter;

      hjmr->SetBinContent(GenPtBin, JMRValue);
      hjmr->SetBinError(GenPtBin, JMRError);

      jmrsum += JMRValue * GenPtBinWidth;
      jmrsumweight += GenPtBinWidth;
  }
//   for (Int_t GenPtBin = h2HCorrelate->GetYaxis()->FindBin(5); GenPtBin <= h2HCorrelate->GetNbinsY(); GenPtBin++) {
//     TH1* hRecPt = (TH1*) h2HCorrelate->ProjectionX(Form("hRecPt_%s", histName), GenPtBin, GenPtBin, "e");
    
//     TF1* crystalBall = new TF1("crystalBall", "[0]*ROOT::Math::crystalball_function(x, [1], [2], [3], [4])", hRecPt->GetXaxis()->GetXmin(), hRecPt->GetXaxis()->GetXmax());
//     // Set initial parameters for the fit: amplitude, alpha, n, mean, sigma
//     crystalBall->SetParameters(hRecPt->GetMaximum(), 1.5, 3.0, hRecPt->GetMean(), hRecPt->GetRMS());
//     hRecPt->Fit(crystalBall, "Q"); // Perform the fit, "Q" option for quiet mode

//     // Extract the sigma from the fit
//     Double_t RecPtSigma = crystalBall->GetParameter(4); // Sigma is usually the 5th parameter
//     Double_t RecPtSigmaError = crystalBall->GetParError(4); // Error on sigma

//     Double_t GenPtCenter = h2HCorrelate->GetYaxis()->GetBinCenter(GenPtBin);
//     Double_t GenPtBinWidth = h2HCorrelate->GetXaxis()->GetBinWidth(GenPtBin);

//     Double_t JMRValue = RecPtSigma / GenPtCenter;
//     Double_t JMRError = RecPtSigmaError / GenPtCenter;

//     hjmr->SetBinContent(GenPtBin, JMRValue);
//     hjmr->SetBinError(GenPtBin, JMRError);

//     jmrsum += JMRValue * GenPtBinWidth;
//     jmrsumweight += GenPtBinWidth;

//     delete crystalBall; // Clean up
// }
  hset(*hjmr, JetPtGenTitleX,"#sigma(#it{p}_{T, jet}^{reco})/#it{p}_{T, jet}^{true}", 1.2, 1.2, 0.05, 0.051, 0.01, 0.01, 0.05, 0.05, 510, 505);
  hjmr->GetXaxis()->SetRangeUser(PlotPtMin, PlotPtMax);
  hjmr->GetYaxis()->SetRangeUser(0., 1.8);
  hjmr->SetMarkerStyle(20);
  hjmr->SetMarkerSize(1.3);
  hjmr->SetLineWidth(2);
  hjmr->SetMarkerColor(kRed);
  hjmr->SetLineColor(kRed);
  hjmr->Draw("pE1");

  ALICEfigureLegend("ALICE WIP", 0.131403,0.698065,0.380846,0.948387, 0.482183,0.692903,0.732739,0.883871);
  legJMR->SetHeader("Jet Momentum Resolution", "C");
  legJMR->AddEntry("", Form("mean: %.3f", jmrsum / jmrsumweight), "");
  legJMR->Draw();
  
  if(DRAWPLOTS) {canJMR->Print(Form("%s/JetMomentumResolution_%s.pdf", MakeDirName.Data(), histName));}

  /////////////////////////////////////
  /// Jet Reconstruction Efficiency ///
  /////////////////////////////////////
  cout << "Starting Jet Reconstruction Efficiency..." << endl;
  
  jrep->cd();
  auto JREp = (TH1 *)JetMCPPt->Clone(Form("hist_%i", ++n));
  JREp->GetXaxis()->SetRangeUser(5,PlotPtMax);
  double entries_in_range = JREp->Integral();
  std::cout << "# entries of JetMCPPt > 0 GeV: " << JREp->GetEntries() << std::endl;
  std::cout << "# entries of JetMCPPt > 5 GeV: " << entries_in_range << std::endl;
  JRElegend->AddEntry("", histName, "");
  JRElegend->AddEntry(JREp, "Generated jets");
  hset(*JREp, JRETitleX, JRETitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05,
       0.05, 510, 505);
  // hoptset(*JREp, NevtsMCD / effTrigZvtx10, kBlack, 10, PlotPtMax, 3e-10, 1e-1);
  hoptset(*JREp, 0, kBlack, 10, PlotPtMax, 3e-10, 1e-1);
  auto mcpmatchedpt = (TH1 *)MCPMatchedpt->Clone(Form("hist_%i", ++n));
  mcpmatchedpt->GetXaxis()->SetRangeUser(0, PlotPtMax);
  std::cout << "# entries of MCPMatchedpt > 5 GeV: " << mcpmatchedpt->GetEntries() << std::endl;
  JRElegend->AddEntry(mcpmatchedpt, "Matched jets");
  hset(*mcpmatchedpt, JRETitleX, JRETitleY, 0.8, 1.4, 0.04, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  // hoptset(*mcpmatchedpt, NevtsMCD, kRed, 10, PlotPtMax, 3e-10, 1e-1);
  hoptset(*mcpmatchedpt, 0, kRed, 10, PlotPtMax, 3e-10, 1e-1);
  JREp->GetYaxis()->SetNdivisions(505);
  JREp->Draw("PE");
  mcpmatchedpt->Draw("PEsame");

  ratiojrep->cd();
  // enforce JetPtPad ratio X-axis (10-140) without touching styles of hist
  TH1F *__jetpt_ratio_frame = new TH1F(Form("frame_ratio_%s", histName), "", 100, PlotPtMin, PlotPtMax);
  __jetpt_ratio_frame->SetDirectory(nullptr);
  __jetpt_ratio_frame->SetMinimum(0.5);
  __jetpt_ratio_frame->SetMaximum(1.01);
  __jetpt_ratio_frame->GetXaxis()->SetTitle(JRETitleX);
  __jetpt_ratio_frame->GetYaxis()->SetTitle("#it{#varepsilon}_{reco}^{jet}");
  __jetpt_ratio_frame->SetLineColor(0);
  __jetpt_ratio_frame->SetMarkerSize(0);
  __jetpt_ratio_frame->Draw("AXIS");
  auto jre = (TH1 *)mcpmatchedpt->Clone(Form("hist_%i", ++n));
  jre->Divide(mcpmatchedpt, JREp, 1., 1., "B");
  hset(*jre, JRETitleX, "#it{#varepsilon}_{reco}^{jet}", 1.2, 1.0, 0.07, 0.07, 0.01,
       0.01, 0.07, 0.07, 510, 510);
  jre->SetMarkerColor(kRed);
  jre->SetLineColor(kRed);
  jre->SetMarkerSize(.7);
  jre->SetMarkerStyle(22);
  jre->GetXaxis()->SetRangeUser(PlotPtMin, PlotPtMax);
  jre->GetYaxis()->SetRangeUser(0.5, 1.01);
  jre->SetFillColorAlpha(kRed, 0.3);
  jre->Draw("pe");

  static TCanvas *sJREcanvas = nullptr; static TLegend *sJREcanLegend = nullptr; static bool sJREfirst = true;
  if (!sJREcanvas) {
    sJREcanvas = new TCanvas("JREcanvas_All", "Jet Reconstruction Efficiency - All datasets", 800, 600);
    sJREcanvas->cd();
    gPad->SetTicks(1, 1);
    setpad(sJREcanvas, 0.02, 0.15, 0.15, 0.05);
    sJREcanLegend = new TLegend(0.593985,0.393043,0.973684,0.733043,NULL,"brNDC");
    sJREcanLegend->SetBorderSize(0);
    sJREcanLegend->SetTextSize(0.04);
    sJREcanLegend->SetFillColorAlpha(0, 0);
  }
  sJREcanvas->cd();
  hoptset(*jre, 0, colorID, 5, PlotPtMax, 0.5, 1.01, 1, 1, 2, 20);
  hset(*jre, JRETitleX, "#it{#varepsilon}_{reco}^{jet}", 1.3, 1.0, 0.05, 0.07, 0.01,
       0.01, 0.05, 0.05, 510, 510);
  jre->SetMarkerStyle(20);
  jre->SetMarkerColor(colorID);
  jre->SetLineColor(colorID);
  jre->Draw(sJREfirst ? "pe" : "pe same");
  sJREcanLegend->AddEntry(jre, histName, "pe");
  sJREcanLegend->Draw();
  ALICEfigureLegend("ALICE WIP", 0.3, 0.2, 0.55, 0.52, 0.6, 0.2, 0.85, 0.44);
  sJREfirst = false;
  if (DRAWPLOTS) {
    sJREcanvas->SaveAs(Form("%s/JetReconstructionEfficiency_AllDatasets.pdf", MakeDirName.Data()));
  }


  ////////////////////////////////
  /// Jet Kinematic Efficiency ///
  ////////////////////////////////
  cout << "Starting Jet Kinematic Efficiency..." << endl;

  Double_t LcutKine = 7;
  Double_t RcutKine = 140;
    Filipad2 *JKEPad = new Filipad2(Form("Jet kinematic efficiency - %s", histName), ++nn, 2, 0.5, 100, 50, 0.7, 1, 1);
      JKEPad->Draw();
      TPad *jkepad = JKEPad->GetPad(1);
      optFili(*jkepad, 1, 1, 0, 1);
      TPad *ratiojkepads = JKEPad->GetPad(2);
      optFili(*ratiojkepads, 1, 1, 0, 0);
      TLegend *legjke =
          new TLegend(0.318182, 0.707826, 0.595694, 0.954783, NULL, "brNDC");
      legjke->SetTextSize(0.065);
      legjke->SetBorderSize(0);
      TLegend *legjke2 =
          new TLegend(0.488038, 0.481159, 0.535885, 0.701449, NULL, "brNDC");
      legjke2->SetTextSize(0.065);
      legjke2->SetBorderSize(0);
      legjke2->SetTextAlign(12);
      legjke2->AddEntry("", Form("%.f GeV #leq #it{p}_{T, jet}^{reco} #leq %.f GeV", LcutKine, RcutKine), "");
    jkepad->cd();
    TH1 *rec_total_window =
        (TH1 *)hResponseMatrix->ProjectionY("rec_total_window", 0, -1, "e");
        // (TH1 *)hResponseMatrix->ProjectionY("rec_total_window", 1, hResponseMatrix->GetNbinsX(), "e");
    //   TH1 *hmissed = (TH1 *)Response->Miss()->Clone("missed_events");
    TH1 *rec_selected_window =
        (TH1 *)hResponseMatrix->ProjectionY("rec_selected_window", hResponseMatrix->GetXaxis()->FindBin(LcutKine+1e-5), hResponseMatrix->GetXaxis()->FindBin(RcutKine-1e-5), "e");
    

    TH1 *JKEp = (TH1F *)rec_total_window->Clone(Form("hist_%i", ++n));
    legjke->AddEntry("", histName, "");
    legjke->AddEntry(JKEp, "Total #it{p}_{T, jet}^{reco} window");
    hset(*JKEp, JRETitleX, JRETitleY, 0.9, 1., 0.05, 0.07, 0.01, 0.01, 0.05,
    0.05,
         510, 510);
  hoptset(*JKEp, NevtsMCP, kBlack, 5, PlotPtMax, 3e-10, 1e-1);

    legjke->AddEntry(rec_selected_window, "Selected #it{p}_{T, jet}^{reco} window");
    hset(*rec_selected_window, JRETitleX, JRETitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
  hoptset(*rec_selected_window, NevtsMCP, kRed, 5, PlotPtMax, 1e-9, 1e-1);
    JKEp->GetYaxis()->SetNdivisions(505);
    JKEp->Draw("pe");
    rec_selected_window->Draw("pesame");

    ratiojkepads->cd();
    auto jke = (TH1F *)rec_selected_window->Clone(Form("hist_%i", ++n));
    jke->Divide(jke, JKEp, 1., 1., "B");
    hset(*jke, JRETitleX, "Kinematic efficiency", 1.2, 1.0, 0.06, 0.07, 0.01,
         0.01, 0.07, 0.07, 510, 510);
    jke->SetMarkerColor(kRed);
    jke->SetLineColor(kRed);
    jke->SetMarkerSize(.7);
    jke->SetMarkerStyle(22);
    jke->GetXaxis()->SetRangeUser(5, PlotPtMax);
    jke->GetYaxis()->SetRangeUser(0., 1.01);
    jke->SetFillColorAlpha(kRed, 0.3);
    jke->Draw("pe");

    jkepad->cd();
    legjke->Draw();
    legjke2->Draw();
    if(DRAWPLOTS) {JKEPad->C->Print(Form("%s/JKE_R%.1f_%s.pdf", MakeDirName.Data(), RBIN, histName));}

    static TCanvas *sJKEcanvas = nullptr; static TLegend *sJKEcanLegend = nullptr; static bool sJKEfirst = true;
    if (!sJKEcanvas) {
      sJKEcanvas = new TCanvas("JKEcanvas_All", "Jet Kinematic Efficiency - All datasets", 800, 600);
      sJKEcanvas->cd();
      gPad->SetTicks(1, 1);
      setpad(sJKEcanvas, 0.02, 0.15, 0.15, 0.05);
      sJKEcanLegend = new TLegend(0.593985,0.393043,0.973684,0.733043,NULL,"brNDC");
      sJKEcanLegend->SetBorderSize(0);
      sJKEcanLegend->SetTextSize(0.04);
      sJKEcanLegend->SetFillColorAlpha(0, 0);
    }
    sJKEcanvas->cd();
  hoptset(*jke, 0, colorID, 5, 200, 0.0, 1.01, 1, 1, 2, 20);
    hset(*jke, JRETitleX, "Kinematic efficiency", 1.3, 1.0, 0.05, 0.07, 0.01,
        0.01, 0.05, 0.05, 510, 510);
    jke->SetMarkerStyle(21);
    jke->SetMarkerColor(colorID);
    jke->SetLineColor(colorID);
    jke->Draw(sJKEfirst ? "pe" : "pe same");
    TLine *LineKineL = new TLine(PlotPtMin, 0.00, PlotPtMin, 1.01);
      LineKineL->SetLineColor(kGray + 3);
      LineKineL->SetLineStyle(2);
      LineKineL->SetLineWidth(2);
      LineKineL->Draw("lsame");
    TLine *LineKineR = new TLine(PlotPtMax, 0.00, PlotPtMax, 1.01);
      LineKineR->SetLineColor(kGray + 3);
      LineKineR->SetLineStyle(2);
      LineKineR->SetLineWidth(2);
      LineKineR->Draw("lsame");
    sJKEcanLegend->AddEntry(jke, histName, "pe");
    sJKEcanLegend->Draw();
    ALICEfigureLegend("ALICE WIP", 0.2, 0.46, 0.45, 0.70, 0.2, 0.26, 0.45, 0.44);
    sJKEfirst = false;
    if (DRAWPLOTS) {
      sJKEcanvas->SaveAs(Form("%s/JetKinematicEfficiency_AllDatasets.pdf", MakeDirName.Data()));
    }

  /////////////////////////////////
  /// Jet Reconstruction Purity ///
  /////////////////////////////////
  cout << "Starting Jet Reconstruction Purity..." << endl;

  jrpp->cd();
  auto JRPp = (TH1F *)JetMCDPt->Clone(Form("hist_%i", ++n));
  JRPlegend->AddEntry("", histName, "");
  JRPlegend->AddEntry(JRPp, "Detector level jets");
  hset(*JRPp, JRPTitleX, JRPTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05,
       0.05, 510, 505);
  // hoptset(*JRPp, NevtsMCD, kBlack, 5, PlotPtMax, 5e-9, 1e-3, 1, 1, 1, 20);
  hoptset(*JRPp, 0, kBlack, 5, PlotPtMax, 5e-9, 1e-3, 1, 1, 1, 20);
  auto mcdmatchedpt = (TH1F *)MCDMatchedpt->Clone(Form("hist_%i", ++n));
  JRPlegend->AddEntry(mcdmatchedpt, "Matched jets in Detector level");
  hset(*mcdmatchedpt, JRPTitleX, JRPTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  // hoptset(*mcdmatchedpt, NevtsMCD, kRed, 5, PlotPtMax, 3e-10, 1e-3, 1, 1, 1, 24);
  hoptset(*mcdmatchedpt, 0, kRed, 5, PlotPtMax, 3e-10, 1e-3, 1, 1, 1, 24);
  JRPp->GetYaxis()->SetNdivisions(505);
  JRPp->Draw("pe");
  mcdmatchedpt->Draw("pesame");
  ratiojrpp->cd();
  auto jrp = (TH1F *)MCDMatchedpt->Clone(Form("hist_%i", ++n));
  jrp->Divide(jrp, JetMCDPt, 1., 1., "B");
  hset(*jrp, JRPTitleX, "Jet purity", 1.2, 1.0, 0.07, 0.07, 0.01, 0.01, 0.07,
       0.07, 510, 505);
  jrp->SetMarkerColor(kRed);
  jrp->SetLineColor(kRed);
  jrp->SetMarkerSize(.7);
  jrp->SetMarkerStyle(24);
  jrp->GetXaxis()->SetRangeUser(0., PlotPtMax);
  jrp->GetYaxis()->SetRangeUser(0, 1.01);
  jrp->SetFillColorAlpha(kRed, 0.3);
  jrp->Draw("pe");

  static TCanvas *sJRPcanvas = nullptr; static TLegend *sJRPcanLegend = nullptr; static bool sJRPfirst = true;
  if (!sJRPcanvas) {
    sJRPcanvas = new TCanvas("JRPcanvas_All", "Jet Reconstruction Purity - All datasets", 800, 600);
    sJRPcanvas->cd();
    gPad->SetTicks(1, 1);
    setpad(sJRPcanvas, 0.02, 0.15, 0.15, 0.05);
    sJRPcanLegend = new TLegend(0.593985,0.393043,0.973684,0.733043,NULL,"brNDC");
    sJRPcanLegend->SetBorderSize(0);
    sJRPcanLegend->SetTextSize(0.04);
    sJRPcanLegend->SetFillColorAlpha(0, 0);
  }
  sJRPcanvas->cd();
  hoptset(*jrp, 0, colorID, 5, PlotPtMax, 0., 1.01, 1, 1, 2, 20);
  hset(*jrp, JRPTitleX, "Jet purity", 1.3, 1.0, 0.05, 0.07, 0.01,
      0.01, 0.05, 0.05, 510, 510);
  jrp->SetMarkerStyle(24);
  jrp->SetMarkerColor(colorID);
  jrp->SetLineColor(colorID);
  jrp->Draw(sJRPfirst ? "pe" : "pe same");
  sJRPcanLegend->AddEntry(jrp, histName, "pe");
  sJRPcanLegend->Draw();
  ALICEfigureLegend("ALICE WIP", 0.3, 0.2, 0.55, 0.52, 0.6, 0.2, 0.85, 0.44);
  sJRPfirst = false;
  if (DRAWPLOTS) {
    sJRPcanvas->SaveAs(Form("%s/JetReconstructionPurity_AllDatasets.pdf", MakeDirName.Data()));
  }

  ////////////////////////////////////////////////
  /// 2D pT correlation plot / Response Matrix ///
  ////////////////////////////////////////////////
  cout << "Starting 2D pT correlation plot / Response Matrix..." << endl;

  TCanvas *canResponseMatrix = new TCanvas(Form("Correlation_%s", histName),
                                        Form("Jet pT correlation - %s", histName), 1100, 1100);
  TLegend *legRM =
          new TLegend(0.0510018,0.735814,0.397086,0.995349,NULL,"brNDC");
  canResponseMatrix->cd();
  setpad(canResponseMatrix, 0.27, 0.12, 0.135, 0.2);
  // canResponseMatrix->SetLogx(1);
  // canResponseMatrix->SetLogy(1);
  canResponseMatrix->SetLogz(1);
  auto hRM = (TH2 *) hResponseMatrix->Clone(Form("hRM_%s", histName));
  hRM->GetXaxis()->SetTitleSize(0.05);
  hRM->GetYaxis()->SetTitleSize(0.05);
  hRM->GetZaxis()->SetTitleSize(0.05);
  hRM->GetXaxis()->SetTitleOffset(1.0);
  hRM->GetYaxis()->SetTitleOffset(1.3);
  hRM->GetZaxis()->SetTitleOffset(1.5);
  hRM->GetZaxis()->SetLabelOffset(0);
  hRM->GetXaxis()->SetLabelSize(0.045);
  hRM->GetYaxis()->SetLabelSize(0.05);
  hRM->GetZaxis()->SetLabelSize(0.045);
  hRM->GetXaxis()->SetNdivisions(505);
  hRM->GetYaxis()->SetNdivisions(505);
  hRM->GetZaxis()->SetNdivisions(505);
  // hRM->SetTitleOffset(0.1);
  hRM->SetTitle(Form("%s", histName));
  hRM->GetXaxis()->SetTitle("#it{p}_{T, jet}^{reco} (GeV/#it{c})");
  hRM->GetYaxis()->SetTitle("#it{p}_{T, jet}^{true} (GeV/#it{c})");
  hRM->GetZaxis()->SetTitle("Probability density");
  hRM->Scale(1. / h2HCorrelate->Integral(), "width");
  gStyle->SetPalette(kRainBow);
  hRM->Draw("colz");
  // Set the range after drawing the histogram
  // hRM->GetXaxis()->SetRangeUser(0., 300.);
  // hRM->GetYaxis()->SetRangeUser(0., 300.);
  hRM->GetZaxis()->SetRangeUser(1e-10, 1e0);

  // Compute diagonal performance metrics on matched (common-edge) bins only
  auto axX = hRM->GetXaxis();
  auto axY = hRM->GetYaxis();
  const int nxb = axX->GetNbins();
  const int nyb = axY->GetNbins();

  auto collectEdges = [](TAxis* ax){
    std::vector<double> e; e.reserve(ax->GetNbins()+1);
    for (int b = 1; b <= ax->GetNbins(); ++b) {
      e.push_back(ax->GetBinLowEdge(b));
    }
    e.push_back(ax->GetBinUpEdge(ax->GetNbins()));
    return e;
  };
  const auto ex = collectEdges(axX);
  const auto ey = collectEdges(axY);
  // intersection of edges with tolerance
  const double tol = 1e-9;
  std::vector<double> commonEdges;
  size_t ix = 0, iy = 0;
  while (ix < ex.size() && iy < ey.size()) {
    const double dx = ex[ix];
    const double dy = ey[iy];
    if (std::fabs(dx - dy) <= tol) { commonEdges.push_back(dx); ++ix; ++iy; }
    else if (dx < dy) { ++ix; }
    else { ++iy; }
  }
  // need at least two edges to form one common bin
  double totalSum = 0.0;
  double diagSum  = 0.0;
  double ddrSum   = 0.0; int ddrCnt = 0;
  if (commonEdges.size() >= 2) {
    // define the global common X-range across both axes
    const double xCommonMin = commonEdges.front();
    const double xCommonMax = commonEdges.back();
    const int ix1 = axX->FindFixBin(xCommonMin + 1e-9);
    const int ix2 = axX->FindFixBin(xCommonMax - 1e-9);

    for (size_t k = 0; k + 1 < commonEdges.size(); ++k) {
      const double low  = commonEdges[k];
      const double high = commonEdges[k+1];
      // find matching x-bin and y-bin whose exact edges match [low, high]
      auto findExactBin = [&](TAxis* ax)->int{
        const int nb = ax->GetNbins();
        for (int b = 1; b <= nb; ++b) {
          const double bl = ax->GetBinLowEdge(b);
          const double bu = ax->GetBinUpEdge(b);
          if (std::fabs(bl - low) <= tol && std::fabs(bu - high) <= tol) return b;
        }
        return 0;
      };
      const int bxDiag = findExactBin(axX);
      const int by     = findExactBin(axY);
      if (bxDiag == 0 || by == 0) continue; // skip if not exactly matched

      const double diagVal = hRM->GetBinContent(bxDiag, by);
      // Row sum across the entire common X-range for this matched Y row
      const double rowSum = hRM->Integral(ix1, ix2, by, by);

      if (rowSum <= 0.0) continue;
      totalSum += rowSum;
      diagSum  += diagVal;
      ddrSum   += diagVal / rowSum;
      ++ddrCnt;
    }
  }
  const double avgDDR = (ddrCnt > 0 ? ddrSum / ddrCnt : 0.0);
  const double offDiagFrac = (totalSum > 0.0 ? 1.0 - diagSum / totalSum : 1.0);

  auto XYline = new TF1("XYline", "x", 0, 300);
  XYline->SetLineColor(kBlack);
  XYline->SetLineStyle(2);
  XYline->Draw("same");

  legRM->SetTextSize(0.04);
  legRM->SetBorderSize(0);
  legRM->SetTextAlign(12);
  legRM->SetFillColorAlpha(0,0); 
  legRM->AddEntry("", "ALICE WIP", "");
  legRM->AddEntry("", "pp #sqrt{#it{s}} = 13.6 TeV, Response matrix", "");
  legRM->AddEntry("", "#it{p}_{T, track} > 0.15 GeV/#it{c}", "");
  legRM->AddEntry("", "|#it{#eta}_{track}| < 0.9, |#it{#eta}_{jet}| < 0.5", "");
  legRM->AddEntry("", "Anti-#it{k}_{T}, charged-particle jet, #it{R} = 0.4", "");
  legRM->Draw();

  // Annotate metrics on the canvas (NDC)
  TLatex latRM; latRM.SetNDC(true); latRM.SetTextSize(0.035);
  latRM.SetTextColor(kWhite); // 글씨 색상을 흰색으로 설정
  // latRM.DrawLatex(0.35, 0.7, Form("Avg DDR: %.2f", avgDDR));
  // latRM.DrawLatex(0.35, 0.6, Form("Off-diagonal fraction: %.1f%%", 100.0*offDiagFrac));

  canResponseMatrix->Update();
  if(DRAWRM) {canResponseMatrix->Print(Form("%s/JetPtCorrelation_R%.1f_%s.pdf",
                             MakeDirName.Data(), RBIN, histName));}

  /////////////////////////
  /// Consistency Check ///
  /////////////////////////
  cout << "Starting Consistency Check..." << endl;

  consistencyp->cd();
  auto rawMCP = (TH1F *)JetMCPPt->Clone(Form("hist_%i", ++n));
  CONSISTlegend->AddEntry(rawMCP, "MC Generated");
  hset(*rawMCP, JRETitleX, JRETitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05,
       0.05, 510, 505);
  hoptset(*rawMCP, NevtsMCD, kBlack, 10, PlotPtMax, 3e-10, 1e0);
  // rawMCP->Scale(1. / rawMCP->Integral(), "width");
  // rawMCP->SetMarkerSize(1.3);
  // rawMCP->SetMarkerStyle(26);
  rawMCP->Draw("pe");
  auto UnfoldMC = (TH1F *)hMCcorrected->Clone(Form("hist_%i", ++n));
  CONSISTlegend->AddEntry("", histName, "");
  CONSISTlegend->AddEntry(UnfoldMC, "Unfolded MC");

  hset(*UnfoldMC, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hoptset(*UnfoldMC, NevtsMCD, kRed, 10, PlotPtMax, 3e-10, 1e0, 0.75, 1, 1, 24);

  UnfoldMC->GetYaxis()->SetNdivisions(505);
  UnfoldMC->Draw("pesame");
  ratioconsistp->cd();
  auto ratconsist = (TH1F *)rawMCP->Clone(Form("hist_%i", ++n));
  ratconsist->Divide(ratconsist, UnfoldMC, 1., 1., "B");
  hset(*ratconsist, JetPtTitleX, "Data", 1.2, 1.0, 0.07, 0.07, 0.01, 0.01, 0.07,
       0.07, 510, 505);
  ratconsist->SetMarkerColor();
  ratconsist->SetLineColor(kRed);
  ratconsist->SetMarkerSize(.7);
  ratconsist->SetMarkerStyle(24);
  ratconsist->GetXaxis()->SetRangeUser(10., PlotPtMax);
  ratconsist->GetYaxis()->SetRangeUser(0.8, 1.2);
  ratconsist->Draw("pe");
  

  //////////////////////
  /// Data Unfolding ///
  //////////////////////
  cout << "Starting Data Unfolding..." << endl;

  // Data luminosity-calculator
  auto fMcLumiCounter = TFile::Open(McLumiCounterFile, "open"); 
  TH1* hMcLumiCounter = (TH1 *) fMcLumiCounter->Get("jet-luminosity-calculator/counter");
  auto NcollTVX = (double) hMcLumiCounter->GetBinContent(6);
  auto NcollTVXSelZvtx = (double) hMcLumiCounter->GetBinContent(7);
  // double NCollTVXSel = 2.52082e+09; // temporary value in 2024 Prel.
  float effTVXtoZvtxsel8 = NcollTVXSelZvtx / NcollTVX; // TVX -> (Zvtx+sel8) efficiency: 0.798 in 2022
  cout << "effTVXtoZvtxsel8 (MC): " << effTVXtoZvtxsel8 << endl;

  // Data luminosity histograms
  TH1 *hLumiTVX = (TH1 *) Dfile->Get(Form("%s", LumiTVXObj));
  Double_t lumiTVX_ub = hLumiTVX->Integral(1, hLumiTVX->GetNbinsX()); // μb⁻¹
  Double_t lumiTVX = lumiTVX_ub * 1e3; // mb⁻¹
  std::cout << "lumiTVX: " << lumiTVX << " /mb (" << lumiTVX_ub << " /ub)" << std::endl;

  TH1 *hNtvx = (TH1 *) Dfile->Get(Form("%s", NTVXObj));
  Double_t Ntvx = (Double_t) hNtvx->Integral(1, hNtvx->GetNbinsX());
  std::cout << "N_{TVX} (Data): " << Ntvx << std::endl;

  // BC-level luminosity after sel8 BC cuts (no collision reco bias, includes pileup correction)
  TH1 *hLumiTVXafterBC = (TH1 *) Dfile->Get(Form("%s", LumiTVXafBCcutsObj));
  Double_t lumiTVXafterBC_ub = 0;
  if (hLumiTVXafterBC) {
    lumiTVXafterBC_ub = hLumiTVXafterBC->Integral(1, hLumiTVXafterBC->GetNbinsX());
    std::cout << "lumiTVXafterBCcuts: " << lumiTVXafterBC_ub << " /ub" << std::endl;
  } else {
    std::cerr << "[Warning] hLumiTVXafterBCcuts not found, falling back to hLumiTVX" << std::endl;
    lumiTVXafterBC_ub = lumiTVX_ub;
  }

  // kNormSoft: rescale BC-level luminosity for correct σ_vis
  // σ_vis_CCDB_est ≈ N_TVX / hLumiTVX (includes ~2% pileup bias, cancels in ratio)
  // L_correct = hLumiTVXafterBCcuts × (σ_vis_CCDB_est / σ_vis_correct)
  Double_t sigmaVisCCDB_est = (lumiTVX_ub > 0) ? Ntvx / lumiTVX_ub / 1e3 : 0; // mb
  Double_t rescale = (kSigmaVis > 0) ? sigmaVisCCDB_est / kSigmaVis : 1.0;
  Double_t L_correct_ub = lumiTVXafterBC_ub * rescale;
  Double_t L_correct_mb = L_correct_ub * 1e3; // μb⁻¹ → mb⁻¹
  std::cout << "[kNormSoft] sigma_vis_CCDB (est): " << sigmaVisCCDB_est << " mb" << std::endl;
  std::cout << "[kNormSoft] sigma_vis_correct: " << kSigmaVis << " mb" << std::endl;
  std::cout << "[kNormSoft] rescale: " << rescale << std::endl;
  std::cout << "[kNormSoft] L_correct: " << L_correct_mb << " /mb" << std::endl;

  // Try both capitalizations (older tasks use Zvertex, newer use zvertex)
  TH1 *hZvtx10 = (TH1 *) Dfile->Get(Form("%s/h_collisions_Zvertex", RefDir));
  if (!hZvtx10) hZvtx10 = (TH1 *) Dfile->Get(Form("%s/h_collisions_zvertex", RefDir));
  // Also try the event-selection directory
  if (!hZvtx10) hZvtx10 = (TH1 *) Dfile->Get("event-selection-task/hColZaftSel");
  if (!hZvtx10) hZvtx10 = (TH1 *) Dfile->Get("eventselection-run3/eventselection/hColZaftSel");

  Double_t effZvtx10 = 1.0;
  if (hZvtx10) {
    TF1 *gausFitZvtx10 = new TF1("gausFitZvtx10", "gaus", -10, 10);
    gausFitZvtx10->SetParameter(0, hZvtx10->GetMaximum());
    gausFitZvtx10->SetParameter(1, hZvtx10->GetMean());
    gausFitZvtx10->SetParameter(2, hZvtx10->GetStdDev());
    hZvtx10->Fit(gausFitZvtx10, "Q0", "", -10, 10);
    Double_t Nzvtx10val = gausFitZvtx10->Integral(-10, 10);
    Double_t NzvtxTotal = gausFitZvtx10->Integral(-999, 990);
    if (NzvtxTotal > 0) effZvtx10 = Nzvtx10val / NzvtxTotal;
    cout << "effZvtx10 (Data, from Gaussian fit): " << effZvtx10 << endl;
  } else {
    std::cerr << "[Warning] h_collisions_Zvertex/zvertex not found in " << RefDir
              << " — using effZvtx10=1.0" << std::endl;
  }

  // Load histogram with priority: h_mcColl_counts_weighted > h_mcColl_counts > h_mccollisions
  // (Old format has priority over newer format)
  TH1 *hMcCollounts = nullptr;
  const char* histNames[] = {
    "h_mcColl_counts_weighted",  // Priority 1: old format, weighted
    "h_mcColl_counts",            // Priority 2: old format, non-weighted
    "h_mccollisions"              // Priority 3: newer format
  };
  const char* usedHistName = nullptr;
  
  for (int i = 0; i < 3; ++i) {
    hMcCollounts = (TH1 *) MCfile->Get(Form("%s/%s", Dir, histNames[i]));
    if (hMcCollounts) {
      usedHistName = histNames[i];
      std::cout << "Using " << usedHistName << " for MC collision counts" << std::endl;
      break;
    }
  }
  
  if (!hMcCollounts) {
    std::cerr << "Error: Cannot find any of the following histograms in " << Dir << ":" << std::endl;
    for (int i = 0; i < 3; ++i) {
      std::cerr << "  - " << histNames[i] << std::endl;
    }
    return nullptr;
  }
  Double_t Nmccollcount = hMcCollounts->GetBinContent(hMcCollounts->FindBin(5.5));
  std::cout << "Nmccollcount (selMC, bin 5.5): " << Nmccollcount << std::endl;
  Double_t NmccollAll = hMcCollounts->GetBinContent(hMcCollounts->FindBin(0.5));
  std::cout << "NmccollAll (INEL, bin 0.5): " << NmccollAll << std::endl;

  // Data track-efficiency
  auto fDatacollcount = TFile::Open(DataCollCounterFile, "open");
  TH1 *DataCollCounter = (TH1 *) fDatacollcount->Get(Form("%s/%s", "track-efficiency", EventObj));
  // Double_t NDatacollTVX = DataCollCounter->GetBinContent(DataCollCounter->FindBin(1.5));
  // Double_t NDatacollTVX = DataCollCounter->GetBinContent(DataCollCounter->FindBin(3.5)); // non-splited mcCollisions w/o EvSel
  Double_t NDatacollTVX = DataCollCounter->GetBinContent(DataCollCounter->FindBin(0.5));
  cout << "NDatacollTVX: " << NDatacollTVX << endl;
  // Double_t NDatacollsel8 = DataCollCounter->GetBinContent(DataCollCounter->FindBin(4.5));
  Double_t NDatacollsel8 = DataCollCounter->GetBinContent(DataCollCounter->FindBin(1.5));
  Double_t effTVXsel8 = NDatacollsel8 / NDatacollTVX;
  cout << "eff(sel8 / TVX): " << effTVXsel8 << endl;

  // MC track-efficiency
  auto fMccollcount = TFile::Open(mccollcounterfile, "open");
  TH1 *McCollCounter = (TH1 *) fMccollcount->Get("track-efficiency/hMcCollCutsCounts");
  Double_t NmccollINEL = McCollCounter->GetBinContent(McCollCounter->FindBin(1.5));
  Double_t NmccollTVX = McCollCounter->GetBinContent(McCollCounter->FindBin(2.5));
  cout << "NmccollINEL: " << NmccollINEL << endl;
  Double_t NmccollTVXnonSplit = McCollCounter->GetBinContent(McCollCounter->FindBin(3.5));
  Double_t effSplit = NmccollTVXnonSplit / NmccollTVX;
  cout << "effSplit: " << effSplit << endl;
  Double_t NmccollselMC = McCollCounter->GetBinContent(McCollCounter->FindBin(4.5));
  cout << "NmccollselMC: " << NmccollselMC << endl;
  Double_t effTVXselMC = NmccollselMC / NmccollTVX;
  cout << "eff(selMC / TVX): " << effTVXselMC << endl;

  // new TCanvas("hTrigEff", "hTrigEff", 800, 600);
  // auto hTriggerEfficiency = GetTriggerEfficiency();
  // hTriggerEfficiency = hTriggerEfficiency->Rebin(nptBinsGen, "hTriggerEfficiency", ptbinGen);
  // hTriggerEfficiency->Scale(1., "width");
  // hTriggerEfficiency->Draw("pe");

  ////////////////////////////////////////////////
  ///// Run 3 corrected Data and MC Invariant Yield /////
  ////////////////////////////////////////////////
  auto InvYData = (TH1F *)hDcorrected->Clone(Form("hist_%i", ++n));
  // hoptset(*InvYData, NCollTVXSel / effSplit / 1/*effTrigZvtx10*/, ColorPallete[1], PlotPtMin, PlotPtMax, 6e-9, 5e-3, 0.6, 1, 2, 24);
  // hoptset(*InvYData, NevtsData * 1.4341 /*effSplit*/ /1/*effTrigZvtx10*/, ColorPallete[1], PlotPtMin, PlotPtMax, 6e-9, 5e-3, 0.6, 1, 2, 24);
  hoptset(*InvYData, NevtsData, ColorPallete[1], PlotPtMin, PlotPtMax, 6e-9, 5e-3, 0.6, 1, 2, 24);
  
  // Store InvYData globally for comparison plot
  gInvYData = InvYData;

  auto InvYMC = (TH1*) JetMCPPt->Clone(Form("hist_%i", ++n));
  // hoptset(*InvYMC, NmccollselMC / 1/*effTrigZvtx10*/, ColorPallete[0], PlotPtMin, PlotPtMax, 3e-10, 1e-1, 0.75, 1, 2, 33);
  hoptset(*InvYMC, Nmccollcount * 1.4341 / 1/*effTrigZvtx10*/, ColorPallete[0], PlotPtMin, PlotPtMax, 3e-10, 1e-1, 0.75, 1, 2, 33);

  //////////////////////////////////////////////
  ///// Run 3 corrected Data Cross Section /////
  //////////////////////////////////////////////
  // kNormSoft: d²σ/dηdpT = N_jets / L_correct / ε_zvtx10 / Δη / ΔpT
  // No ε_INEL→TVX (TVX captures all jet-producing events)
  // No ε_sel8 (BC-level luminosity already includes sel8 BC cuts)
  double normFactor = 1.0 / L_correct_mb / effZvtx10 / deltaEta;
  cout << "[kNormSoft] normFactor = 1 / " << L_correct_mb << " / " << effZvtx10 << " / " << deltaEta << " = " << normFactor << endl;
  auto UnfoldData = (TH1 *)hDcorrected->Clone(Form("Run3_CrossSection_%s", histName));
  UnfoldData->Scale(normFactor, "width");
  hset(*UnfoldData, JetPtDataFinalTitleX, XSECTION ? XSectionTitleY : JetPtDataFinalTitleY, 0.9, 1.5, 0.05, 0.045, 0.01, 0.005, 0.045, 0.040, 510, 1005);
  hoptset(*UnfoldData, 0, ColorPallete[1], PlotPtMin, PlotPtMax, XSECTION? 3e-7 : 1e-9, XSECTION? 5e-1 : 1e-1, 0.6, 1, 2, 24);

  // Extract run numbers from file names
  // refFile is TString (e.g., "498133_AnalysisResults.root")
  TString refRunNumber = refFile;
  Ssiz_t refIndex = refRunNumber.Index("_AnalysisResults.root");
  if (refIndex != kNPOS) {
    refRunNumber.Remove(refIndex);
  }
  
  // MCfileName is const char* (e.g., "~/cernbox/.../515446_AnalysisResults.root")
  TString mcFileNameStr(MCfileName);
  TString mcBaseName = gSystem->BaseName(mcFileNameStr.Data()); // Get filename without path
  TString mcRunNumber = mcBaseName;
  Ssiz_t mcIndex = mcRunNumber.Index("_AnalysisResults.root");
  if (mcIndex != kNPOS) {
    mcRunNumber.Remove(mcIndex);
  }
  
  gSystem->MakeDirectory("InvariantYieldResults");
  InvYData->SaveAs(Form("InvariantYieldResults/Run3_InvariantYield_data_%s_MC_%s.root", refRunNumber.Data(), mcRunNumber.Data()));
  
  // Create output directory if it doesn't exist
  gSystem->MakeDirectory("XsectionResults");
  UnfoldData->SaveAs(Form("XsectionResults/Run3_CrossSection_data_%s_MC_%s.root", refRunNumber.Data(), mcRunNumber.Data()));


  // if (SYSTUNFOLD) {
  /////////////////////////////////////////////////////////
  ///// Systematic Study on Tracking Efficiency       /////
  /////////////////////////////////////////////////////////

  Filipad2 *SystTrkEffPad = new Filipad2(Form("Syst - Trk. Eff. - %s", histName), ++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
  SystTrkEffPad->Draw();
  TPad *systTrkEffpad = SystTrkEffPad->GetPad(1);
  optFili(*systTrkEffpad, 1, 1, 0, 1);
  TPad *systTrkEffRatpad = SystTrkEffPad->GetPad(2);
  optFili(*systTrkEffRatpad, 1, 1, 0, 0);

  // Legend for the plot
  TLegend *legsystTrkEff = new TLegend(0.313397, 0.6, 0.578947, 0.843478, NULL, "brNDC");
  legsystTrkEff->SetTextSize(0.05);
  legsystTrkEff->SetBorderSize(0);
  legsystTrkEff->SetFillColorAlpha(0, 0);

  // Draw the first histogram (without Tracking Efficiency - Default)
  systTrkEffpad->cd();
  auto defaultWoTrackEff = (TH1 *)UnfoldData->Clone(Form("hist_%i", ++n));  // Default case: without Tracking Efficiency
  legsystTrkEff->AddEntry(defaultWoTrackEff, "default", "lpe");
  hset(*defaultWoTrackEff, JetPtTitleX, XSectionTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
      0.05, 0.05, 510, 505);
  hoptset(*defaultWoTrackEff, 0, kBlack, PlotPtMin, PlotPtMax, 3e-7, 1e1, 0.7, 1, 2, 24);
  defaultWoTrackEff->Draw("pe");

  // Draw the second histogram (with Tracking Efficiency - Variation)
  auto UnfoldWTrackEff = (TH1 *)Run3XSectionWTrackEff();  // Modified function for with tracking efficiency
  legsystTrkEff->AddEntry(UnfoldWTrackEff, "Tracking efficiency -3%", "lpe");
  hoptset(*UnfoldWTrackEff, 0, kRed, PlotPtMin, PlotPtMax, 1e-9, 1e-1, 0.7, 1, 2, 25);
  UnfoldWTrackEff->Draw("pesame");

  legsystTrkEff->Draw();

  // Ratio plot for systematic uncertainty
  systTrkEffRatpad->cd();
  auto hRatWTrkEff = (TH1 *)UnfoldWTrackEff->Clone(Form("hist_%i", ++n));
  hRatWTrkEff->Divide(UnfoldWTrackEff, defaultWoTrackEff, 1, 1, "B");
  hset(*hRatWTrkEff, JetPtGenTitleX, "variation / default", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
      0.1, 0.1, 510, 505);
  hoptset(*hRatWTrkEff, 0, kRed, PlotPtMin, PlotPtMax, 0.98, 1.22, 0.7, 1, 2, 25);
  hRatWTrkEff->Draw("pe");

  // Save output for tracking efficiency systematic uncertainty
  OutStatsTXT(hRatWTrkEff, "TrackingEfficiency");

  if (SYSTUNFOLD) {
      SystTrkEffPad->C->SaveAs(Form("%s/SystErrTrackingEfficiency.pdf", MakeDirName.Data()));
  }

  /////////////////////////////////////////////////////////
  ///// Systematic Study on Track momentum resolution /////
  /////////////////////////////////////////////////////////

  Filipad2 *SystTrkResPad = new Filipad2(Form("Syst - Trk. Res. - %s", histName), ++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
      SystTrkResPad->Draw();
      TPad *systTrkRespad = SystTrkResPad->GetPad(1);
      optFili(*systTrkRespad, 1, 1, 0, 1);
      TPad *systTrkResRatpad = SystTrkResPad->GetPad(2);
      optFili(*systTrkResRatpad, 1, 1, 0, 0);
  TLegend *legsystTrkRes =
          new TLegend(0.313397,0.6,0.578947,0.843478,NULL,"brNDC");
          legsystTrkRes->SetTextSize(0.05);
          legsystTrkRes->SetBorderSize(0);
          legsystTrkRes->SetFillColorAlpha(0,0); 
  
  systTrkRespad->cd();
  auto defaultWtrackTuner = (TH1 *)UnfoldData->Clone(Form("hist_%i", ++n));
  legsystTrkRes->AddEntry(defaultWtrackTuner, "default", "lpe");
  hset(*defaultWtrackTuner, JetPtTitleX, XSectionTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hoptset(*defaultWtrackTuner, 0, kBlack, PlotPtMin, PlotPtMax, 3e-7, 1e0, 0.7, 1, 2, 24);
  defaultWtrackTuner->Draw("pe");

  auto UnfoldWoTrackTuner = (TH1 *) Run3XSectionWoTrackTuner(); // old track tuner comp.
  // auto UnfoldWoTrackTuner = (TH1 *) Run3XSectionWTrackPtSmear1p5(); // new pT smear comp.
  legsystTrkRes->AddEntry(UnfoldWoTrackTuner, "w/ track pT smearing +50%", "lpe");
  // hset(*UnfoldWoTrackTuner, JetPtTitleX, XSectionTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
  //      0.05, 0.05, 510, 505);
  hoptset(*UnfoldWoTrackTuner, 0, kRed, PlotPtMin, PlotPtMax, 1e-9, 1e0, 0.7, 1, 2, 25);
  UnfoldWoTrackTuner->Draw("pesame");

  legsystTrkRes->Draw();

  systTrkResRatpad->cd();
  auto hRatWoTrkTuner = (TH1 *) UnfoldWoTrackTuner->Clone(Form("hist_%i", ++n));
  hRatWoTrkTuner->Divide(hRatWoTrkTuner, defaultWtrackTuner, 1, 1, "B");
  hset(*hRatWoTrkTuner, JetPtGenTitleX, "variation / default", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
       0.1, 0.1, 510, 505);
  hoptset(*hRatWoTrkTuner, 0, kRed, PlotPtMin, PlotPtMax, 0.8, 1.2, 0.7, 1, 2, 25);
  hRatWoTrkTuner->Draw("pe");
  OutStatsTXT(hRatWoTrkTuner, "TrackPtResolution");

  // if(SYSTUNFOLD) {SystTrkResPad->C->SaveAs(Form("%s/SystErrTrackPtResolution.pdf", MakeDirName.Data()));} 
  // if(true) {SystTrkResPad->C->SaveAs("SystErrTrackPtResolution_temp.pdf");} 

  ///////////////////////////////////////////////////
  ///// Systematic Study on Unfolding Iteration /////
  ///////////////////////////////////////////////////

  Filipad2 *SystIterPad = new Filipad2(Form("Syst - Iteration - %s", histName), ++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
      SystIterPad->Draw();
      TPad *systiterpad = SystIterPad->GetPad(1);
      optFili(*systiterpad, 1, 1, 0, 1);
      TPad *systiterratpad = SystIterPad->GetPad(2);
      optFili(*systiterratpad, 1, 1, 0, 0);
  TLegend *legsystiter =
          new TLegend(0.313397,0.6,0.578947,0.843478,NULL,"brNDC");
          legsystiter->SetTextSize(0.05);
          legsystiter->SetBorderSize(0);
          legsystiter->SetFillColorAlpha(0,0); 
  
  systiterpad->cd();
  auto defaultUnfold = (TH1 *)hDcorrected->Clone(Form("hist_%i", ++n));
  legsystiter->AddEntry(defaultUnfold, Form("Data Unfolded (k_{reg}: %d)", dataK), "lpe");
  hset(*defaultUnfold, JetPtTitleX, JetPtDataFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hoptset(*defaultUnfold, NevtsData / effTVXselMC, kBlack, PlotPtMin, PlotPtMax, 1e-9, 1e-1, 0.7, 1, 2, 24);
  defaultUnfold->Draw("pe");

  auto Unfold2 = (TH1 *)hDcorrected2->Clone(Form("hist_%i", ++n));
  legsystiter->AddEntry(Unfold2, Form("Data Unfolded (k_{reg}: %d)", dataK - 1), "lpe");
  hset(*Unfold2, JetPtTitleX, JetPtDataFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hoptset(*Unfold2, NevtsData / effTVXselMC, kRed, PlotPtMin, PlotPtMax, 1e-9, 1e-1, 0.7, 1, 2, 27);
  Unfold2->Draw("pesame");

  auto Unfold3 = (TH1 *)hDcorrected3->Clone(Form("hist_%i", ++n));
  legsystiter->AddEntry(Unfold3, Form("Data Unfolded (k_{reg}: %d)", dataK + 1), "lpe");
  hset(*Unfold3, JetPtTitleX, JetPtDataFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hoptset(*Unfold3, NevtsData / effTVXselMC, kBlue, PlotPtMin, PlotPtMax, 1e-9, 1e-1, 0.7, 1, 2, 28);
  Unfold3->Draw("pesame");

  legsystiter->Draw();

  systiterratpad->cd();
  auto Unfoldrat2 = (TH1 *) Unfold2->Clone(Form("hist_%i", ++n));
  Unfoldrat2->Divide(Unfoldrat2, defaultUnfold, 1, 1, "B");
  hset(*Unfoldrat2, JetPtGenTitleX, "variation / default", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
       0.1, 0.1, 510, 505);
  hoptset(*Unfoldrat2, 0, kRed, PlotPtMin, PlotPtMax, 0.92, 1.08, 0.7, 1, 2, 27);
  Unfoldrat2->Draw("pe");
  OutStatsTXT(Unfoldrat2, Form("UnfoldKreg_%d", dataK - 1));

  auto Unfoldrat3 = (TH1 *) Unfold3->Clone(Form("hist_%i", ++n));
  Unfoldrat3->Divide(Unfoldrat3, defaultUnfold, 1, 1, "B");
  hset(*Unfoldrat3, JetPtGenTitleX, "variation / default", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
       0.1, 0.1, 510, 505);
  hoptset(*Unfoldrat3, 0, kBlue, PlotPtMin, PlotPtMax, 0.92, 1.08, 0.7, 1, 2, 28);
  Unfoldrat3->Draw("pesame");
  OutStatsTXT(Unfoldrat3, Form("UnfoldKreg_%d", dataK + 1));

  if(SYSTUNFOLD) {SystIterPad->C->SaveAs(Form("%s/SystErrUnfoldingIteration.pdf", MakeDirName.Data()));}

  ////////////////////////////////////////////////
  ///// Systematic Study on Unfolding Method /////
  ////////////////////////////////////////////////
  Filipad2 *SystMethodPad = new Filipad2(Form("Syst - Unfold method - %s", histName), ++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
  SystMethodPad->Draw();
  TPad *systmethodpad = SystMethodPad->GetPad(1);
  optFili(*systmethodpad, 1, 1, 0, 1);
  TPad *systmethodratpad = SystMethodPad->GetPad(2);
  optFili(*systmethodratpad, 1, 1, 0, 0);

  TLegend *legsystmethod =
      new TLegend(0.480861,0.617391,0.746411,0.86087,NULL,"brNDC");
  legsystmethod->SetTextSize(0.05);
  legsystmethod->SetBorderSize(0);
  legsystmethod->SetFillColorAlpha(0,0); 

  systmethodpad->cd();
  
  // Plot Bayesian Unfolding result (default)
  auto hBayesUnfold = (TH1F *)hDcorrected->Clone(Form("hist_%i", ++n));
  legsystmethod->AddEntry(hBayesUnfold, "Bayesian unfold");
  hset(*hBayesUnfold, JetPtTitleX, JetPtDataFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
        0.05, 0.05, 510, 505);
  hoptset(*hBayesUnfold, NevtsData / effTVXselMC, kBlack, PlotPtMin, PlotPtMax, 3e-10, 1e-1, 0.7, 1, 2, 24);
  hBayesUnfold->Draw("pe");

  // Plot SVD Unfolding result
  auto hSVDUnfold = (TH1F *)hDcorrectedSVD->Clone(Form("hist_%i", ++n));
  legsystmethod->AddEntry(hSVDUnfold, "SVD unfold");
  hset(*hSVDUnfold, JetPtTitleX, JetPtDataFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
        0.05, 0.05, 510, 505);
  hoptset(*hSVDUnfold, NevtsData / effTVXselMC, kBlue, PlotPtMin, PlotPtMax, 3e-10, 1e-1, 0.7, 1, 2, 27);
  hSVDUnfold->Draw("pesame");

  legsystmethod->Draw();

  systmethodratpad->cd();
  auto UnfoldratSVD = (TH1 *) hSVDUnfold->Clone(Form("hist_%i", ++n));
  UnfoldratSVD->Divide(UnfoldratSVD, hBayesUnfold, 1, 1, "B");
  hset(*UnfoldratSVD, JetPtGenTitleX, "SVD / Bayesian", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
        0.1, 0.1, 510, 505);
  hoptset(*UnfoldratSVD, 0, kBlue, PlotPtMin, PlotPtMax, 0.88, 1.12, 0.7, 1, 2, 27);
  UnfoldratSVD->Draw("pe");

  OutStatsTXT(UnfoldratSVD, "UnfoldRatSVD");
  if(SYSTUNFOLD) {SystMethodPad->C->SaveAs(Form("%s/SystErrUnfoldingMethod.pdf", MakeDirName.Data()));}


  ////////////////////////////////////////////////
  ///// Systematic Study on Unfolding Prior //////
  //////////////////////////////////////////////// 
  auto hUnfoldedDataHerwig = DrawUnfoldHerwig(InvYMC, JetMCPPt, JetMCDPt, h2HCorrelate, DJetpt, NevtsMCD / effTVXselMC);

  Filipad2 *SystPriorPad = new Filipad2(Form("Syst - Unfold prior dep. - %s", histName), ++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
      SystPriorPad->Draw();
      TPad *systpriorpad = SystPriorPad->GetPad(1);
      optFili(*systpriorpad, 1, 1, 0, 1);
      TPad *systpriorratpad = SystPriorPad->GetPad(2);
      optFili(*systpriorratpad, 1, 1, 0, 0);
  TLegend *legsystprior =
          new TLegend(0.480861,0.617391,0.746411,0.86087,NULL,"brNDC");
          legsystprior->SetTextSize(0.05);
          legsystprior->SetBorderSize(0);
          legsystprior->SetFillColorAlpha(0,0); 
  
  systpriorpad->cd();
  // auto defaultUnfold = (TH1 *)hDcorrected->Clone(Form("hist_%i", ++n));
  legsystprior->AddEntry(defaultUnfold, "PYTHIA8 prior");
  // hset(*defaultUnfold, JetPtTitleX, JetPtDataFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
  //      0.05, 0.05, 510, 505);
  hoptset(*defaultUnfold, 0, kBlack, PlotPtMin, PlotPtMax, 3e-10, 1e-1, 0.7, 1, 2, 24);
  defaultUnfold->Draw("pe");

  auto UnfoldHerwig = (TH1 *)hUnfoldedDataHerwig->Clone(Form("hist_%i", ++n));
  legsystprior->AddEntry(UnfoldHerwig, "Herwig prior");
  hset(*UnfoldHerwig, JetPtTitleX, JetPtDataFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hoptset(*UnfoldHerwig, NevtsData / effTVXselMC, ColorPallete[3], PlotPtMin, PlotPtMax, 3e-10, 1e-1, 0.7, 1, 2, 27);
  UnfoldHerwig->Draw("pesame");

  legsystprior->Draw();

  systpriorratpad->cd();
  // auto UnfoldratHerwig = (TH1 *) UnfoldHerwig->Clone(Form("hist_%i", ++n));
  // UnfoldratHerwig->Divide(UnfoldratHerwig, defaultUnfold, 1, 1, "B");
  auto UnfoldratHerwig = DrawRatioTH1(UnfoldHerwig, defaultUnfold); 
  hset(*UnfoldratHerwig, JetPtGenTitleX, "Herwig / PYTHIA8", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
       0.1, 0.1, 510, 505);
  hoptset(*UnfoldratHerwig, 0, ColorPallete[3], PlotPtMin, PlotPtMax, 0.5, 1.5, 0.7, 1, 2, 27);
  UnfoldratHerwig->Draw("pe");

  OutStatsTXT(UnfoldratHerwig, "UnfoldRatPriorHerwig");
  if (SYSTUNFOLD) {
  SystPriorPad->C->SaveAs(Form("%s/SystErrUnfoldingHerwig.pdf", MakeDirName.Data()));
  }

  /////////////////////////////////////////////////////////////////////
  ///// Systematic Uncertainty: Secondary Particles Contamination /////
  /////////////////////////////////////////////////////////////////////
  auto hSystErrSecCon = (TH1*) DrawSecondaryContaimination(UnfoldData);

  //////////////////////////////////////////
  ///// Total Systematic Uncertainty   /////
  ////////////////////////////////////////// 
  std::vector<TH1*> hSysts;
  hSysts.push_back(hRatWTrkEff);
  hSysts.push_back(hRatWoTrkTuner);
  // hSysts.push_back(Unfoldrat3); // negligeble
  // hSysts.push_back(Unfoldrat2); // negligeble
  hSysts.push_back(UnfoldratSVD);
  // hSysts.push_back(hSystErrSecCon); // constant
  hSysts.push_back(UnfoldratHerwig); // constant

  double SystMcClosure = 1.04;
  TH1 *hSystMcClosure = new TH1F(Form("hSystMcClosure_%s", histName), "hSystMcClosure", nptBinsGen, ptbinGen);
  for (int i=0; i<=hSystMcClosure->GetNbinsX(); i++) {
    hSystMcClosure->SetBinContent(i, SystMcClosure);
  }
  hSysts.push_back(hSystMcClosure);

  // double SystMcPrior = 1.05;
  // TH1 *hSystMcPrior = new TH1F("hSystMcPrior", "hSystMcPrior", nptBinsGen, ptbinGen);
  // for (int i=0; i<=hSystMcPrior->GetNbinsX(); i++) {
  //   hSystMcPrior->SetBinContent(i, SystMcPrior);
  // }
  // hSysts.push_back(hSystMcPrior);

  // double SystTrkPtRes = 1.02;
  // TH1 *hSystTrkPtRes = new TH1F("hSystTrkPtRes", "hSystTrkPtRes", nptBinsGen, ptbinGen);
  // for (int i=0; i<=hSystTrkPtRes->GetNbinsX(); i++) {
  //   hSystTrkPtRes->SetBinContent(i, SystTrkPtRes);
  // }
  // hSysts.push_back(hSystTrkPtRes);

  double SystNorm = 1.1;
  TH1 *hSystNorm = new TH1F(Form("hSystNorm_%s", histName), "hSystNorm", nptBinsGen, ptbinGen);
  for (int i=0; i<=hSystNorm->GetNbinsX(); i++) {
    hSystNorm->SetBinContent(i, SystNorm);
  }
  hSysts.push_back(hSystNorm);

  double SystSecCon = 1.05;
  TH1 *hSystSecCon = new TH1F(Form("hSystSecCon_%s", histName), "hSystSecCon", nptBinsGen, ptbinGen);
  for (int i=0; i<=hSystSecCon->GetNbinsX(); i++) {
    hSystSecCon->SetBinContent(i, SystSecCon);
  }
  hSysts.push_back(hSystSecCon);

  // auto SystTrkEffFile = TFile::Open("systematics/SystTrackingEfficiency.root", "read");
  // TIter nextStd(SystTrkEffFile->GetListOfKeys()); 
  // TKey *keyStd; 
  // TH1 *hSystTrkEff = nullptr; 
  // while ((keyStd = (TKey*)nextStd())) {
  //     TString name = keyStd->GetName();
  //     if (name.Contains("JetMCD")) {
  //         hSystTrkEff = (TH1*)SystTrkEffFile->Get(name);
  //         break;
  //     }
  // }
  // hSysts.push_back(hSystTrkEff);

  TH1* hSystResult = (TH1*)hSysts[0]->Clone(Form("quadratureSum_%s", histName));
    hSystResult->Reset();
  calculateQuadratureSum(hSysts, hSystResult);
  OutStatsTXT(hSystResult, "TotalSystematicUncertainty");

  DrawMultipleSources(hSysts, hSystResult);

  ///////////////////////////////////////////////////////////
  ///// Run 3 corrected Invariant Yield: Data vs gen MC /////
  ///////////////////////////////////////////////////////////
  Filipad2 *InvariantYieldPad = new Filipad2(Form("Invariant Yield - %s", histName), ++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
  InvariantYieldPad->Draw();
  TPad *yieldpad = InvariantYieldPad->GetPad(1);
  optFili(*yieldpad, 0, 0, 0, 1);
  TPad *yieldratpad = InvariantYieldPad->GetPad(2);
  optFili(*yieldratpad, 0, 0, 0, 0);

  TLegend *legendYield = new TLegend(0.543062,0.378882,0.80622,0.575155,NULL,"brNDC");
  legendYield->SetTextSize(0.04);
  legendYield->SetBorderSize(0);
  legendYield->SetFillColorAlpha(0, 0);
  TLegend *legendYield2 =
      new TLegend(0.160287,0.0434782,0.258373,0.234783,NULL,"brNDC");
  legendYield2->SetTextSize(0.043);
  legendYield2->SetBorderSize(0);
  legendYield2->SetTextAlign(12);
  legendYield2->SetFillColorAlpha(0,0); 
  legendYield2->AddEntry("", "|#it{#eta}_{jet}| < 0.5", "");
  legendYield2->AddEntry("", "Anti-#it{k}_{T}, #it{R} = 0.4", "");
  legendYield2->AddEntry("", "charged-particle jets", "");
  TLegend *legendYield3 =
      new TLegend(0.535885,0.664596,0.605263,0.915528,NULL,"brNDC");
  legendYield3->SetTextSize(0.043);
  legendYield3->SetBorderSize(0);
  legendYield3->SetTextAlign(12);
  legendYield3->SetFillColorAlpha(0,0); 
  legendYield3->AddEntry("", "ALICE WIP", "");
  legendYield3->AddEntry("", "pp #sqrt{#it{s}} = 13.6 TeV", "");
  legendYield3->AddEntry("", "#it{p}_{T, track} > 0.15 GeV/#it{c}", "");
  legendYield3->AddEntry("", "|#it{#eta}_{track}| < 0.9", "");

  yieldpad->cd();
  gPad->SetTicks(1, 1);

  hset(*InvYData, JetPtDataFinalTitleX, JetPtDataFinalTitleY, 0.9, 1.45, 0.05, 0.045, 0.01, 0.005,
       0.045, 0.040, 510, 505);
  InvYData->GetYaxis()->SetNdivisions(505);
  legendYield->AddEntry(InvYData, "Run 3,#kern[-0.7]{ }#sqrt{#it{s}}=13.6 TeV");
  InvYData->Draw("pe");

  TH1 *InvYdataWsyst = ApplySystematicUncertainty(InvYData, hSystResult);
  InvYdataWsyst->SetMarkerStyle(0);
  InvYdataWsyst->Draw("E2same");

  auto Ygraphs = Run2Data_manual();
  auto Ygrun2Syst = (TGraphErrors *) Ygraphs[1]->Clone(Form("syst_hist_%i", ++n));
  auto Yhrun2Syst = GraphToHistogram(Ygrun2Syst);
  auto Ygrun2Stat = (TGraphErrors *) Ygraphs[0]->Clone(Form("stat_hist_%i", ++n));
  auto Yhrun2Stat = GraphToHistogram(Ygrun2Stat);
  hoptset(*Yhrun2Syst, 0, ColorPallete[2], PlotPtMin, PlotPtMax, 1e-7, 1e0, 0.6, 1, 2, 25);
  Yhrun2Syst->Scale(1./ 58.1);
  Yhrun2Syst->Draw("pe2same"); 
  hoptset(*Yhrun2Stat, 0, ColorPallete[2], PlotPtMin, PlotPtMax, 1e-7, 1e0, 0.6, 1, 2, 25); 
  Yhrun2Stat->Scale(1./ 58.1);
  Yhrun2Stat->Draw("e same"); 
  legendYield->AddEntry(Yhrun2Syst, "Run 2,#kern[-0.7]{ }#sqrt{#it{s}}=13 TeV");

  InvYMC->Draw("pesame");
  legendYield->AddEntry(InvYMC, "PYTHIA8", "lep");

  legendYield->Draw();
  legendYield2->Draw();
  legendYield3->Draw();

  // auto hHerwig1360InvY = GetHERWIG1360();
  // hoptset(*hHerwig1360InvY, 0, ColorPallete[3], PlotPtMin, PlotPtMax, 1e-7, 1e1, 0.5, 1, 2, 24);
  // hHerwig1360InvY->Draw("pesame");

  yieldratpad->cd();
  gPad->SetTicks(1, 1);

  auto ratioYield = (TH1*)InvYMC->Clone(Form("hist_%i", ++n));
  ratioYield->Divide(InvYMC, InvYData, 1, 1, "B");
  hset(*ratioYield, JetPtDataFinalTitleX, "Comp. / Run 3 Data", 1.19, 0.76, 0.1, 0.085, 0.01, 0.01,
      0.1, 0.1, 510, 505);
  hoptset(*ratioYield, 0, ColorPallete[0], PlotPtMin, PlotPtMax, 0.6, 1.9, 0.65, 1, 2, 33);
  ratioYield->Draw("pe");

  auto InvYdataSelf = (TH1 *) InvYData->Clone(Form("hist_%i", ++n));
  InvYdataSelf->Divide(InvYdataSelf, InvYData, 1, 1, "B");
  TH1 *InvYDataSystRatio = ApplySystematicUncertainty(InvYdataSelf, hSystResult);
  InvYDataSystRatio->SetMarkerStyle(0);
  InvYDataSystRatio->Draw("E2same");

  TH1 *InvYDataStatRatio = (TH1 *)InvYData->Clone(Form("hist_stat_%i", n));
  InvYDataStatRatio->Reset(); 
  for (int i = 1; i <= InvYDataStatRatio->GetNbinsX(); ++i) {
      InvYDataStatRatio->SetBinContent(i, 1.0); 
      double error = InvYData->GetBinError(i) / InvYData->GetBinContent(i);
      InvYDataStatRatio->SetBinError(i, error); 
  }
  InvYDataStatRatio->SetMarkerStyle(0); 
  InvYDataStatRatio->Draw("E1same");

  auto ratioYieldRun2 = DrawRatioTH1(Yhrun2Syst, InvYData);
  hoptset(*ratioYieldRun2, 0, ColorPallete[2], PlotPtMin, PlotPtMax, 0.4, 2.7, 0.5, 1, 2, 25);
  ratioYieldRun2->Draw("pe2same");

  auto ratioYieldRun2Stat = DrawRatioTH1(Yhrun2Stat, InvYData);
  hoptset(*ratioYieldRun2Stat, 0, ColorPallete[2], PlotPtMin, PlotPtMax, 0.4, 2.7, 0.5, 1, 2, 25);
  ratioYieldRun2Stat->Draw("e same");

  // auto ratioYieldRun2 = DrawRatioTH1(Yrun2Data, InvYData);
  // hset(*ratioYieldRun2, JetPtGenTitleX, "MC / Data", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01, 0.1, 0.1, 510, 505);
  // hoptset(*ratioYieldRun2, 0, ColorPallete[2], PlotPtMin, PlotPtMax, 0.7, 1.9, 0.5, 1, 2, 25);
  // ratioYieldRun2->Draw("pe2same"); 

  OutStatsTXT(ratioYield, "Run3_vs_MC_InvariantYield_Ratio");

  if (DRAWPLOTS) {
    InvariantYieldPad->C->SaveAs(Form("%s/InvariantYield_R%.1f.pdf", MakeDirName.Data(), RBIN));
    InvariantYieldPad->C->SaveAs(Form("%s/InvariantYield_R%.1f.eps", MakeDirName.Data(), RBIN));
  }


  /////////////////////////////////////////
  ///// Cross Section: Data vs gen MC /////
  /////////////////////////////////////////
  unfoldp->cd();
  gPad->SetTicks(1, 1);
  unfoldlegend->AddEntry(UnfoldData, "Run 3,#kern[-0.7]{ }#sqrt{#it{s}}=13.6 TeV");
  UnfoldData->Draw("pe");
  // UnfoldData->SaveAs("ALICE_Run3_xsection(mb).root");

  TH1 *UnfoldDataWsyst = ApplySystematicUncertainty(UnfoldData, hSystResult);
  UnfoldDataWsyst->SetMarkerStyle(0);
  // UnfoldDataWsyst->Draw("E2same");
  // UnfoldDataWsyst->SaveAs("Run3_CrossSection_SystErr_wTrackTuner.root");

  auto graphs = Run2Data_manual();
  auto grun2Syst = (TGraphErrors *) graphs[1]->Clone(Form("syst_hist_%i", ++n));
  auto hrun2Syst = GraphToHistogram(grun2Syst);
  auto grun2Stat = (TGraphErrors *) graphs[0]->Clone(Form("stat_hist_%i", ++n));
  auto hrun2Stat = GraphToHistogram(grun2Stat);
  hoptset(*hrun2Syst, 0, ColorPallete[2], PlotPtMin, PlotPtMax, 1e-7, 1e0, 0.6, 1, 2, 25);
  hrun2Syst->Draw("pe2same"); 
  hoptset(*hrun2Stat, 0, ColorPallete[2], PlotPtMin, PlotPtMax, 1e-7, 1e0, 0.6, 1, 2, 25); 
  hrun2Stat->Draw("e same"); 
  unfoldlegend->AddEntry(hrun2Syst, "Run 2,#kern[-0.7]{ }#sqrt{#it{s}}=13 TeV");

  
  auto NormMCP = (TH1*) JetMCPPt->Clone(Form("hist_%i", ++n));
  // unfoldlegend->AddEntry(NormMCP, SEL8WINDOW? "PYTHIA8" : "MC true", "lep");
  // hoptset(*NormMCP, XSECTION ? ((NmccollINEL / SigmaINEL) * effTVXselMC) : (SEL8WINDOW? NevtsMCD / effTVXselMC : 1 /* TBC */), ColorPallete[0], PlotPtMin, PlotPtMax, XSECTION? 1e-7 : 1e-9, XSECTION? 1e1 : 1e-1, 0.5, 1, 2, 24);
  // MC truth cross-section: d²σ/dηdpT = N_jets_truth / NmccollAll(INEL) × σ_INEL / Δη / ΔpT
  // hoptset denominator N: h.Scale(1/N, "width"), so N = NmccollAll × deltaEta / SigmaINEL
  hoptset(*NormMCP, XSECTION ? (NmccollAll * deltaEta / SigmaINEL) : (SEL8WINDOW? Nmccollcount : 1 /* TBC */), ColorPallete[0], PlotPtMin, PlotPtMax, XSECTION? 1e-7 : 1e-9, XSECTION? 1e1 : 1e-1, 0.75, 1, 2, 33);
  // NormMCP->Draw("pesame");

  TH1 *PYTHIA1360 = GetPYTHIA1360();
  // unfoldlegend->AddEntry(PYTHIA1360, "PYTHIA8 13.6 TeV", "lep");
  hoptset(*PYTHIA1360, 0, ColorPallete[4], PlotPtMin, PlotPtMax,  1e-7, 1e1, 0.5, 1, 2, 24);
  hset(*PYTHIA1360, JetPtDataFinalTitleX, JetPtMCFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01);
  // PYTHIA1360->Draw("pesame");

  auto hPYTHIA13600 = DrawPythia13600();
  hoptset(*hPYTHIA13600, 0, kCyan, PlotPtMin, PlotPtMax,  1e-7, 1e1, 0.5, 1, 2, 24);
  unfoldlegend->AddEntry(hPYTHIA13600, "PYTHIA8 13.6 TeV", "lep");
  hPYTHIA13600->Draw("pesame");

  // TH1 *PYTHIA1300 = GetPYTHIA1300();
  // unfoldlegend->AddEntry(PYTHIA1300, "PYTHIA8 13.0 TeV", "lep");
  // hoptset(*PYTHIA1300, 0, kMagenta, PlotPtMin, PlotPtMax,  1e-7, 1e1, 0.5, 1, 2, 24);
  // hset(*PYTHIA1300, JetPtDataFinalTitleX, JetPtMCFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01);
  // if(XSECTION) {PYTHIA1300->Scale(77.23);}
  // PYTHIA1300->Draw("pesame");

  TH1 *HERWIG1360 = GetHERWIG1360();
  unfoldlegend->AddEntry(HERWIG1360, "Herwig", "lep");
  hoptset(*HERWIG1360, 0, ColorPallete[3], PlotPtMin, PlotPtMax,  1e-7, 1e1, 0.7, 1, 2, 34);
  hset(*HERWIG1360, JetPtDataFinalTitleX, JetPtMCFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01);
  // if(XSECTION) {HERWIG1360->Scale(SigmaINEL);}
  HERWIG1360->Draw("pesame");

  auto hChangwhan2022 = DrawChangwhan2022MB();
  unfoldlegend->AddEntry(hChangwhan2022, "Changwhan 2022", "lep");
  hoptset(*hChangwhan2022, 0, kMagenta, PlotPtMin, PlotPtMax, 0.4, 2.7, 0.4, 1, 2, 20);
  hChangwhan2022->Draw("pesame");

  auto hChangwhan2022JJMC = DrawChangwhan2022JJMC();
  unfoldlegend->AddEntry(hChangwhan2022JJMC, "Changwhan 2022 JJMC", "lep");
  hoptset(*hChangwhan2022JJMC, 0, kCyan+2, PlotPtMin, PlotPtMax, 0.4, 2.7, 0.4, 1, 2, 20);
  hChangwhan2022JJMC->Draw("pesame");

  auto hChangwhan2024 = DrawChangwhan2024MB();
  unfoldlegend->AddEntry(hChangwhan2024, "Changwhan 2024", "lep");
  hoptset(*hChangwhan2024, 0, kBrown, PlotPtMin, PlotPtMax, 0.4, 2.7, 0.4, 1, 2, 20);
  hChangwhan2024->Draw("pesame");


  ratunfoldp->cd();
  gPad->SetTicks(1, 1);
  
  // auto ratPYTHIA1300 = DrawRatioTH1(PYTHIA1300, UnfoldData);
  // hset(*ratPYTHIA1300, JetPtGenTitleX, "Comp. / Run 3", 1.2, 0.75, 0.1, 0.09, 0.01, 0.01,
  //      0.1, 0.1, 510, 505);
  // hoptset(*ratPYTHIA1300, 0, kMagenta, PlotPtMin, PlotPtMax, 0.4, 2.7, 0.4, 1, 2, 20);
  // ratPYTHIA1300->Draw("pe");

  // auto ratMCTVXData = (TH1 *)NormMCP->Clone(Form("hist_%i", ++n));
  // ratMCTVXData->Divide(ratMCTVXData, UnfoldData, 1., 1., "B");
  auto ratMCTVXData = DrawRatioTH1(NormMCP, UnfoldData);
  hset(*ratMCTVXData, JetPtDataFinalTitleX, "Comp. / Run 3 Data", 1.19, 0.76, 0.1, 0.085, 0.01, 0.01,
      0.1, 0.1, 510, 505);
  hoptset(*ratMCTVXData, 0, ColorPallete[0], PlotPtMin, PlotPtMax, 0.5, 1.5, 0.65, 1, 2, 33); // for jet-jet MC QA version
  // ratMCTVXData->Draw("pesame");

  // auto ratMCtrueData = (TH1F *)NormMCPINEL->Clone(Form("hist_%i", ++n));
  // ratMCtrueData->Divide(ratMCtrueData, UnfoldData, 1., 1., "B");
  // hset(*ratMCtrueData, JetPtDataFinalTitleX, "Comp. / Run 3", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
  //     0.1, 0.1, 510, 505);
  // hoptset(*ratMCtrueData, 0, ColorPallete[0], PlotPtMin, PlotPtMax, 0.4, 2.7, 0.4, 1, 2, 24);
  // ratMCtrueData->Draw("pesame");

  auto ratPYTHIA1360 = DrawRatioTH1(PYTHIA1360, UnfoldData);
  hset(*ratPYTHIA1360, JetPtDataFinalTitleX, "Comp. / Run 3", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
       0.1, 0.1, 510, 505);
  hoptset(*ratPYTHIA1360, 0, ColorPallete[4], PlotPtMin, PlotPtMax, 0.4, 2.7, 0.4, 1, 2, 20);
  // ratPYTHIA1360->Draw("pesame");

  auto ratPYTHIA13600 = DrawRatioTH1(hPYTHIA13600, UnfoldData);
  hset(*ratPYTHIA13600, JetPtDataFinalTitleX, "Comp. / Run 3", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
       0.1, 0.1, 510, 505);
  hoptset(*ratPYTHIA13600, 0, kCyan, PlotPtMin, PlotPtMax, 0.5, 1.5, 0.4, 1, 2, 20);
  ratPYTHIA13600->Draw("pesame");

  TH1 *ratHerwigRun3 = DrawRatioTH1(HERWIG1360, UnfoldData);
  hset(*ratHerwigRun3, JetPtDataFinalTitleX, "Herwig / Run 3", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
       0.1, 0.1, 510, 505);
  hoptset(*ratHerwigRun3, 0, ColorPallete[3], PlotPtMin, PlotPtMax, 0.4, 2.7, 0.6, 1, 2, 34);
  // ratHerwigRun3->Draw("pesame");

  auto ratChangwhan2022 = DrawRatioTH1(hChangwhan2022, UnfoldData);
  hoptset(*ratChangwhan2022, 0, kMagenta, PlotPtMin, PlotPtMax, 0.6, 1.25, 0.4, 1, 2, 20);
  hset(*ratChangwhan2022, JetPtDataFinalTitleX, "Comp. / Run 3 Data", 1.19, 0.76, 0.1, 0.085, 0.01, 0.01, 0.1, 0.1, 510, 505);
  ratChangwhan2022->Draw("pesame");

  auto ratChangwhan2022JJMC = DrawRatioTH1(hChangwhan2022JJMC, UnfoldData);
  hoptset(*ratChangwhan2022JJMC, 0, kCyan+2, PlotPtMin, PlotPtMax, 0.4, 2.7, 0.4, 1, 2, 20);
  ratChangwhan2022JJMC->Draw("pesame");

  auto ratChangwhan2024 = DrawRatioTH1(hChangwhan2024, UnfoldData);
  hoptset(*ratChangwhan2024, 0, kBrown, PlotPtMin, PlotPtMax, 0.4, 2.7, 0.4, 1, 2, 20);
  ratChangwhan2024->Draw("pesame");

  auto UnfoldDataSelf = (TH1 *) UnfoldData->Clone(Form("hist_%i", ++n));
  UnfoldDataSelf->Divide(UnfoldDataSelf, UnfoldData, 1, 1, "B");
  TH1 *UnfoldDataSystRatio = ApplySystematicUncertainty(UnfoldDataSelf, hSystResult);
  UnfoldDataSystRatio->SetMarkerStyle(0);
  // UnfoldDataSystRatio->Draw("E2same");

  TH1 *UnfoldDataStatRatio = (TH1 *)UnfoldData->Clone(Form("hist_stat_%i", n));
  UnfoldDataStatRatio->Reset(); 
  for (int i = 1; i <= UnfoldDataStatRatio->GetNbinsX(); ++i) {
      UnfoldDataStatRatio->SetBinContent(i, 1.0); 
      double error = UnfoldData->GetBinError(i) / UnfoldData->GetBinContent(i);
      UnfoldDataStatRatio->SetBinError(i, error); 
  }
  UnfoldDataStatRatio->SetMarkerStyle(0); 
  UnfoldDataStatRatio->Draw("E1same");

  auto ratRun2Run3data = DrawRatioTH1(hrun2Syst, UnfoldData);
  hoptset(*ratRun2Run3data, 0, ColorPallete[2], PlotPtMin, PlotPtMax, 0.4, 2.7, 0.5, 1, 2, 25);
  ratRun2Run3data->Draw("pe2same");

  auto ratRun2Run3dataStat = DrawRatioTH1(hrun2Stat, UnfoldData);
  hoptset(*ratRun2Run3dataStat, 0, ColorPallete[2], PlotPtMin, PlotPtMax, 0.4, 2.7, 0.5, 1, 2, 25);
  ratRun2Run3dataStat->Draw("e same");

  // // MC true / Data ratio
  // auto ratMCTrueData = (TH1F *)hTrueJetPt->Clone(Form("hist_%i", ++n));
  // ratMCTrueData->Divide(ratMCTrueData, UnfoldData, 1., 1., "B");
  // hset(*ratMCTrueData, JetPtGenTitleX, "MC True / Data", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
  //     0.1, 0.1, 510, 505);
  // hoptset(*ratMCTrueData, 0, ColorPallete[3], PlotPtMin, PlotPtMax, 0.5, 3.0);
  // ratMCTrueData->Draw("pesame");


  // UnfoldData->SaveAs("plots/systematic/UnfoldedData_standard.root");
  // UnfoldData->SaveAs("plots/systematic/UnfoldedData_trackingefficiency.root");

  // if (XSECTION) {
  // ////////////////////////////////////////////////////////
  // ///// Run 3 corrected data vs Run 2 corrected data /////
  // ////////////////////////////////////////////////////////
  // Filipad2 *Run2DataValPad = new Filipad2(++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
  //     Run2DataValPad->Draw();
  //     TPad *run2datavalpad = Run2DataValPad->GetPad(1);
  //     optFili(*run2datavalpad, 1, 1, 0, 1);
  //     TPad *run2datavalratpad = Run2DataValPad->GetPad(2);
  //     optFili(*run2datavalratpad, 1, 1, 0, 0);
  // TLegend *legRun2Data1 =
  //         new TLegend(0.476077,0.418634,0.741627,0.662112,NULL,"brNDC");
  //     legRun2Data1->SetTextSize(0.05);
  //     legRun2Data1->SetBorderSize(0);
  //     legRun2Data1->SetFillColorAlpha(0,0); 
  //     TLegend *legRun2Data2 =
  //         new TLegend(0.165072,0.0111801,0.210526,0.229814,NULL,"brNDC");
  //     legRun2Data2->SetTextSize(0.05);
  //     legRun2Data2->SetBorderSize(0);
  //     legRun2Data2->SetTextAlign(12);
  //     legRun2Data2->SetFillColorAlpha(0,0); 
  //     legRun2Data2->AddEntry("", "|#it{#eta}_{jet}| #leq 0.5", "");
  //     legRun2Data2->AddEntry("", "Anti-#it{k}_{T}, #it{R} = 0.4", "");
  //     TLegend *legRun2Data3 =
  //         new TLegend(0.543062,0.711801,0.626794,0.960248,NULL,"brNDC");
  //     legRun2Data3->SetTextSize(0.05);
  //     legRun2Data3->SetBorderSize(0);
  //     legRun2Data3->SetTextAlign(12);
  //     legRun2Data3->SetFillColorAlpha(0,0); 
  //     legRun2Data3->AddEntry("", "pp #sqrt{#it{s}} = 13.6 TeV", "");
  //     legRun2Data3->AddEntry("", "#it{p}_{T, track} #geq 0.15 GeV/#it{c}", "");
  //     legRun2Data3->AddEntry("", "|#it{#eta}_{track}| #leq 0.9", "");
  
  // run2datavalpad->cd();
  // hset(*UnfoldData, JetPtGenTitleX, XSECTION ? XSectionTitleY : JetPtDataFinalTitleY, 0.9, 1.1, 0.06, 0.06, 0.01, 0.01,
  //      0.05, 0.05, 510, 505);
  // UnfoldData->Draw("pe");

  // auto run2data = (TGraphErrors *) Run2Data(DJetpt)->Clone(Form("hist_%i", ++n));
  // hoptset(*run2data, 0, kBlue, PlotPtMin, PlotPtMax, 1e-7, 1e1);
  // run2data->Draw("pe2same");
  // legRun2Data1->AddEntry(UnfoldData, "Run 3 ALICE data");
  // legRun2Data1->AddEntry(run2data, "Run 2 ALICE data");
  // legRun2Data1->Draw();
  // legRun2Data2->Draw();
  // legRun2Data3->Draw();

  // run2datavalratpad->cd();
  // auto run2dataratio = (TGraphErrors *) DrawRatioTGraph(run2data, UnfoldData);
  // hoptset(*run2dataratio, 0, kBlue, PlotPtMin, PlotPtMax, 0, 2);
  // hset(*run2dataratio, JetPtGenTitleX, "Run 2 / Run 3", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
  //      0.1, 0.1, 510, 505);
  // run2dataratio->Draw("APE2");  

  // if (DRAWPLOTS) {Run2DataValPad->C->Print(Form("%s/UnfoldedJetPt_Run3_Run2_Data_R%.1f.pdf",
  //                              MakeDirName.Data(), RBIN));}

  
//   //////////////////////////////////////////
//   /////// Run 3 gen MC vs Run 2 gen MC /////
//   //////////////////////////////////////////
//   Filipad2 *Run2MCValPad = new Filipad2(++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
//       Run2MCValPad->Draw();
//       TPad *run2mcvalpad = Run2MCValPad->GetPad(1);
//       optFili(*run2mcvalpad, 1, 1, 0, 1);
//       TPad *run2mcvalratpad = Run2MCValPad->GetPad(2);
//       optFili(*run2mcvalratpad, 1, 1, 0, 0);
//   TLegend *legRun2MC1 =
//           new TLegend(0.476077,0.418634,0.741627,0.662112,NULL,"brNDC");
//       legRun2MC1->SetTextSize(0.05);
//       legRun2MC1->SetBorderSize(0);
//       legRun2MC1->SetFillColorAlpha(0,0);
  
//   run2mcvalpad->cd();
//   hset(*NormMCP, JetPtTitleX, XSECTION ? XSectionTitleY : JetPtMCFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
//        0.05, 0.05, 510, 505);
//   NormMCP->Draw("pe");

//   auto run2MCgen = Run2MCgen();
//   hoptset(*run2MCgen, 0, kBrown, PlotPtMin, PlotPtMax, 1e-7, 1e1);
//   run2MCgen->Draw("pesame");
//   legRun2MC1->AddEntry(NormMCP, "Run 3 ALICE MC true");
//   legRun2MC1->AddEntry(run2MCgen, "Run 2 ALICE MC true");
//   legRun2MC1->Draw();
//   legRun2Data2->Draw();
//   legRun2Data3->Draw();

//   run2mcvalratpad->cd();
//   auto run2MCratio = (TH1 *) run2MCgen->Clone(Form("hist_%i", ++n));
//   run2MCratio->Divide(run2MCgen, NormMCP, 1., 1., "B");
//   // run2MCratio->Scale(1.25, "");
//   hoptset(*run2MCratio, 0, kBrown, PlotPtMin, PlotPtMax, 0, 2);
//   hset(*run2MCratio, JetPtGenTitleX, "Run 2 / Run 3", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
//        0.1, 0.1, 510, 505);
//   run2MCratio->Draw("pe2same"); 

//   if (DRAWPLOTS) {Run2MCValPad->C->Print(Form("%s/UnfoldedJetPt_Run3_Run2_MC_R%.1f.pdf",
//                                MakeDirName.Data(), RBIN));}

//   // UnfoldData->SaveAs("run3UnfoldData.root");
//   // NormMCP->SaveAs("run3MCgen.root");
//   // run2data->SaveAs("run2data.root");
//   }


return jre;
}

// Response matrix slice comparison: project pT,det for specific pT,part ranges
// Overlays all MC files (area-normalized) to compare MB vs JJ response shapes
void DrawResponseSlices(const std::vector<TString> &fileNames,
                        const std::vector<TString> &histNames,
                        const std::vector<TString> &McFileDirectories,
                        const std::vector<Color_t> &colors,
                        const TString &outputDir) {

  const char *kResponseObj = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint";
  std::vector<std::pair<int, int>> ptPartSlices = {{16, 20}, {50, 60}};

  for (auto &slice : ptPartSlices) {
    TCanvas *can = new TCanvas(Form("ResponseSlice_%d_%d", slice.first, slice.second),
                               Form("ResponseSlice_%d_%d", slice.first, slice.second), 800, 600);
    can->Draw();
    setpad(can, 0.02, 0.12, 0.12, 0.02);

    TLegend *leg = new TLegend(0.50, 0.60, 0.92, 0.92, NULL, "brNDC");
    leg->SetTextSize(0.035);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0, 0);
    leg->SetHeader(Form("%d #leq #it{p}_{T,jet}^{true} < %d GeV/#it{c}", slice.first, slice.second));

    bool first = true;
    for (Int_t i = 0; i < (Int_t)std::min(fileNames.size(), histNames.size()); ++i) {
      TString filePath = mainDir + fileNames[i];
      TString fileDir = McFileDirectories[i];
      auto file = TFile::Open(filePath.Data(), "READ");
      if (!file || file->IsZombie()) { std::cerr << "Cannot open " << filePath << std::endl; continue; }

      TH2 *hResp = (TH2 *)file->Get(Form("%s/%s", fileDir.Data(), kResponseObj));
      if (!hResp) { std::cerr << "Missing " << kResponseObj << " in " << filePath << std::endl; file->Close(); continue; }

      Int_t binYlo = hResp->GetYaxis()->FindBin(slice.first + 1e-6);
      Int_t binYhi = hResp->GetYaxis()->FindBin(slice.second - 1e-6);
      TH1D *hSlice = hResp->ProjectionX(Form("slice_%d_%d_%s", slice.first, slice.second, histNames[i].Data()),
                                         binYlo, binYhi, "e");
      hSlice->SetDirectory(0);
      file->Close();

      if (hSlice->Integral() > 0) hSlice->Scale(1. / hSlice->Integral());

      hset(*hSlice, "#it{p}_{T,jet}^{reco} (GeV/#it{c})", "Probability", 1.10, 1.2, 0.047, 0.05,
           0.001, 0.001, 0.05, 0.05);
      hSlice->SetMarkerColor(colors[i + 1]);
      hSlice->SetLineColor(colors[i + 1]);
      hSlice->SetMarkerSize(0.8);
      hSlice->SetMarkerStyle(20 + i);
      hSlice->GetXaxis()->SetRangeUser(0, PlotPtMax);

      if (first) {
        hSlice->GetYaxis()->SetRangeUser(0., hSlice->GetMaximum() * 1.5);
        hSlice->Draw("e");
        first = false;
      } else {
        hSlice->Draw("esame");
      }
      leg->AddEntry(hSlice, histNames[i].Data(), "pe");
    }
    leg->Draw();
    if (DRAWPLOTS) { can->Print(Form("%s/ResponseSlice_ptPart_%d_%d.pdf", outputDir.Data(), slice.first, slice.second)); }
  }
}

// JES profile: mean (pT,true - pT,reco)/pT,true vs pT,true from h2_jet_pt_mcp_jet_pt_diff_matchedgeo
// Overlays all MC files with ratio panel (ratio to first file)
void DrawJESProfile(const std::vector<TString> &fileNames,
                    const std::vector<TString> &histNames,
                    const std::vector<TString> &McFileDirectories,
                    const std::vector<Color_t> &colors,
                    const TString &outputDir) {

  const char *kJESDiffObj = "h2_jet_pt_mcp_jet_pt_diff_matchedgeo";

  Filipad2 *JESPad = new Filipad2("JES profile", ++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
  JESPad->Draw();
  TPad *mainpad = JESPad->GetPad(1);
  optFili(*mainpad, 1, 1, 0, 0);
  TPad *ratiopad = JESPad->GetPad(2);
  optFili(*ratiopad, 1, 1, 0, 0);

  TLegend *leg = new TLegend(0.50, 0.68, 0.92, 0.95, NULL, "brNDC");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);

  TH1 *refProfile = nullptr;

  for (Int_t i = 0; i < (Int_t)std::min(fileNames.size(), histNames.size()); ++i) {
    TString filePath = mainDir + fileNames[i];
    TString fileDir = McFileDirectories[i];
    auto file = TFile::Open(filePath.Data(), "READ");
    if (!file || file->IsZombie()) { std::cerr << "Cannot open " << filePath << std::endl; continue; }

    TH2 *h2JES = (TH2 *)file->Get(Form("%s/%s", fileDir.Data(), kJESDiffObj));
    if (!h2JES) { std::cerr << "Missing " << kJESDiffObj << " in " << filePath << std::endl; file->Close(); continue; }

    TProfile *prof = h2JES->ProfileX(Form("JESprof_%s", histNames[i].Data()));
    prof->SetDirectory(0);
    file->Close();

    // Rebin to analysis binning for cleaner display
    TH1D *hProf = new TH1D(Form("JEShist_%s", histNames[i].Data()),
                            "", nptBinsGen, ptbinGen);
    for (int ib = 1; ib <= hProf->GetNbinsX(); ++ib) {
      double ptCenter = hProf->GetBinCenter(ib);
      int srcBin = prof->FindBin(ptCenter);
      if (prof->GetBinEntries(srcBin) > 0) {
        hProf->SetBinContent(ib, prof->GetBinContent(srcBin));
        hProf->SetBinError(ib, prof->GetBinError(srcBin));
      }
    }

    mainpad->cd();
    hset(*hProf, "#it{p}_{T,jet}^{true} (GeV/#it{c})",
         "#LT(#it{p}_{T}^{true} - #it{p}_{T}^{reco}) / #it{p}_{T}^{true}#GT",
         1.10, 1.2, 0.047, 0.05, 0.001, 0.001, 0.05, 0.05);
    hProf->SetMarkerColor(colors[i + 1]);
    hProf->SetLineColor(colors[i + 1]);
    hProf->SetMarkerSize(0.8);
    hProf->SetMarkerStyle(20 + i);
    hProf->GetXaxis()->SetRangeUser(5, PlotPtMax);
    hProf->GetYaxis()->SetRangeUser(-0.15, 0.15);

    if (i == 0) {
      hProf->Draw("e");
      refProfile = (TH1 *)hProf->Clone("JESref");
    } else {
      hProf->Draw("esame");
    }
    leg->AddEntry(hProf, histNames[i].Data(), "pe");

    if (refProfile && i > 0) {
      ratiopad->cd();
      DrawRatio(Form("JESratio_%s", histNames[i].Data()), refProfile, hProf,
                "#it{p}_{T,jet}^{true} (GeV/#it{c})", RatioTitleY,
                colors[i + 1], 0.5, 1.5);
    }
  }
  mainpad->cd();
  leg->Draw();
  if (DRAWPLOTS) { JESPad->C->Print(Form("%s/JESProfile.pdf", outputDir.Data())); }
}

TH1 *DrawJetResolution(const char *fileName, const char *histName, Double_t Nevt,
                       TLegend *legend, Color_t colorID,
                       const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  TString histNameStr(histName);

  TH2 *JetResolution = (TH2 *)file->Get(Form("%s/%s", Dir, JetResolutionObj));

  if (!JetResolution) {
    std::cout << "Error: No valid JetResolution object found." << std::endl;
    return nullptr;
  }

  // std::vector<std::pair<int, int>> pTClasses = {{10, 15}, {15, 20},
  // {20, 40}, {40, 60}};
  std::vector<std::pair<int, int>> pTClasses = {
      {5, 10}, {10, 15}, {20, 25}, {40, 50}, {85, 100}};

  for (Int_t i = 0; i < pTClasses.size(); i++) {
    TH1D *jetResolution = (TH1D *)JetResolution->ProjectionY(
        Form("JetResZ%i_%s", i + 1, histName),
        JetResolution->GetXaxis()->FindBin(pTClasses[i].first),
        JetResolution->GetXaxis()->FindBin(pTClasses[i].second));

    hset(*jetResolution, JetResolutionTitleX, "Probability density", 1.10, 1.2, 0.047, 0.05,
         0.001, 0.001, 0.05, 0.05);
    hoptset(*jetResolution, 2, ColorPallete[i], -0.5, 1.0, 1e-5, 10.0, 0.75, 1, 1, 20+i==24? 29 : 20+i);
    legend->AddEntry(jetResolution,
                     Form("#it{p}_{T, jet}^{true} [%i, %i] GeV/#it{c}, #mu: %.2f, #sigma: %.2f",
                          pTClasses[i].first, pTClasses[i].second, jetResolution->GetMean(), jetResolution->GetRMS()),
                     "pe");
    // jetResolution->GetMean(), jetResolution->GetRMS()), "pe2");
    // hset(*jetResolution, JetResolutionTitleX, JetPtTitleY, 0.7, 1.0, 0.05,
    // 0.05, 0.01, 0.01, 0.05, 0.05, 51);
    jetResolution->Draw("same");
  }
  
  // file->Close();
  return nullptr; // This function doesn't return a specific histogram
}
