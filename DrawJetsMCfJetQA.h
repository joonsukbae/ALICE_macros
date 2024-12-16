///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Draw macro for jet QA         //////////
////////// author: Joonsuk Bae           ////////// 
////////// E-mail: jbae@cern.ch          //////////
////////// Last Modified: 19 July 2024   //////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

// Draw Jets
#include <cstdlib>
TH1 *DrawJetPt(const char *fileName, const char *histName, const char *Obj, Double_t Nevts,
               TLegend *legend, Color_t colorID, Int_t i = 0,
               const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  TH1 *JetPt = (TH1 *)file->Get(Form("%s/%s", Dir, Obj));
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
  auto file = TFile::Open(fileName, "open");
  TH1 *JetPtMCP = (TH1 *)file->Get(Form("%s/%s", Dir, Obj));

  if (REBINON) {
    JetPtMCP = JetPtMCP->Rebin(nptBinsGen,Form("JetPtMCP_%s",histName),ptbinGen);
  }
  legend->AddEntry(JetPtMCP, histName);
  hset(*JetPtMCP, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 510); 
  hoptset(*JetPtMCP, Nevts, colorID, 0, PlotPtMax, 1e-12, 1e-1, 1.2 - 0.2 * i);
  JetPtMCP->Draw("esame");

  return JetPtMCP;
}
TH1 *DrawJetEta(const char *fileName, const char *histName, const char *Obj, Double_t Nevts,
                TLegend *legend, Color_t colorID, const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  TH3D *H3JetEta = (TH3D *)file->Get(Form("%s/%s", Dir, Obj));
  TH1 *JetEta = H3JetEta->ProjectionZ(
      Form("JetEta_%s", histName), H3JetEta->GetXaxis()->FindBin(RBIN + 1e-6),
      H3JetEta->GetXaxis()->FindBin(RBIN + 0.1 - 1e-5),
      1, H3JetEta->GetYaxis()->FindBin(PlotPtMax));

  TH1 *h1JetEta = JetEta->Rebin(JetEta->GetNbinsX() / 25, Form("JetEta_rebinned_%s", histName));

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
  TH1D *JetPhi = H3JetPhi->ProjectionZ(
      Form("JetPhi_%s", histName), H3JetPhi->GetXaxis()->FindBin(RBIN + 1e-6),
      H3JetPhi->GetXaxis()->FindBin(RBIN + 0.1 - 1e-5),
      1, H3JetPhi->GetNbinsY());
  TH1 *h1JetPhi = JetPhi->Rebin(JetPhi->GetNbinsX() / 20, Form("JetPhi_rebinned_%s", histName));
  legend->AddEntry(h1JetPhi, histName);
  hset(*h1JetPhi, JetPhiTitleX, JetPhiTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 510);
  NORMEVENTS? hoptset(*JetPhi, Nevts, colorID, 0, 2 * TMath::Pi(), 0., 2e-2) : hoptset(*JetPhi, Nevts, colorID, 0, 2 * TMath::Pi(), 0., 0.3);
  JetPhi->Draw("esame");

  return JetPhi;
}
TH1 *DrawJetNtracks(const char *fileName, const char *histName, const char *Obj, Double_t Nevts,
                    TLegend *legend, Color_t colorID,
                    const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  TH1 *JetNtracks = (TH1 *)file->Get(Form("%s/%s", Dir, Obj));
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

  TString path = Dir ? Form("%s/%s", Dir, Obj) : TString(Obj);
  TH3F *H3JetArea = (TH3F *)file->Get(path);
  if (!H3JetArea) {
    std::cerr << "Error getting histogram: " << path << std::endl;
    return nullptr;
  }

  auto JetArea = (TH2F *)H3JetArea->Project3D("zy");
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
                     TLegend *JRPlegend = nullptr, TFile *savefile = nullptr,
                     const char *Dir = nullptr) {

  auto Dfile = TFile::Open(DfileName, "open");
  TH1 *DJetpt =
      (TH1 *)Dfile->Get(Form("%s/%s", DataDirectory[0].Data(), "h_jet_pt"));
  if (REBINON) {
    DJetpt = DJetpt->Rebin(nptBins, Form("DJetpt_%s", histName), ptbin);
  }
  double NjetsData = DJetpt->Integral(DJetpt->FindBin(5), DJetpt->GetNbinsX());
  cout << "The number of jets in measured Data: " << NjetsData << endl;

  auto MCfile = TFile::Open(MCfileName, "open");

  TH1 *JetMCPPtINEL = (TH1 *)MCfile->Get(Form("%s/%s", "jet-finder-charged-qa", "h_jet_pt_part"));
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
  TH3 *HCorrelate = nullptr;
  TH3 *MatchTemp1 = (TH3 *)MCfile->Get(
      Form("%s/%s", Dir, "h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"));
  TH3 *MatchTemp2 =
      (TH3 *)MCfile->Get(Form("%s/%s", Dir, "h3_jet_r_jet_pt_part_jet_pt"));
  if (MatchTemp1) {
    HCorrelate = MatchTemp1;
  } else if (MatchTemp2) {
    HCorrelate = MatchTemp2;
  }
  // cout << "Nentries of 3DHcorrelate: " << HCorrelate->GetEntries() << endl;
  HCorrelate->GetXaxis()->SetRange(RBIN, RBIN);
  
  TH2F* h2HCorrelate = nullptr; 
  if (REBINON) {
    h2HCorrelate = new TH2F(Form("hcorrelate_%s", histName),
                                Form("R projected correlate_%s", histName),
                                nptBins, ptbin, nptBinsGen, ptbinGen);

    Int_t corrbin = HCorrelate->GetXaxis()->FindBin(RBIN + 1e-6);
    for (Int_t i = 1; i <= HCorrelate->GetNbinsY(); i++) {   // part.
      for (Int_t j = 1; j <= HCorrelate->GetNbinsZ(); j++) { // det.
        Double_t content = HCorrelate->GetBinContent(corrbin, i, j);
        Double_t error = HCorrelate->GetBinError(corrbin, i, j);

        Int_t binpart = h2HCorrelate->GetYaxis()->FindBin(
            HCorrelate->GetYaxis()->GetBinCenter(i));
        Int_t bin = h2HCorrelate->GetXaxis()->FindBin(
            HCorrelate->GetZaxis()->GetBinCenter(j));

        Double_t currentContent = h2HCorrelate->GetBinContent(bin, binpart);
        Double_t currentError = h2HCorrelate->GetBinError(bin, binpart);

        Double_t newContent = currentContent + content;
        Double_t newError = sqrt(pow(currentError, 2) +
                                pow(error, 2)); // Combine errors in quadrature

        h2HCorrelate->SetBinContent(bin, binpart, newContent);
        h2HCorrelate->SetBinError(bin, binpart, newError);
      }
    }
  } else {
    h2HCorrelate = (TH2F *)HCorrelate->Project3D(Form("%s_yze", histName));
  }

  TH1* MCDMatchedpt = nullptr;
  if (REBINON) {
    MCDMatchedpt = (TH1 *) h2HCorrelate->ProjectionX(Form("hMCDMatched_%s", histName), 1, h2HCorrelate->GetNbinsY(), "e");
    // MCDMatchedpt = MCDMatchedpt->Rebin(nptBins, Form("MCDptMatched_%s", histName), ptbin);
  } else {
    MCDMatchedpt = (TH1 *)HCorrelate->ProjectionZ(Form("hMCDMatched_%s", histName),
      HCorrelate->GetXaxis()->FindBin(RBIN + 1e-6),
      HCorrelate->GetXaxis()->FindBin(RBIN + 0.2 - 1e-6), 1,
      HCorrelate->GetNbinsY(), "e");
  }

  TH1 *MCPMatchedpt = (TH1 *) h2HCorrelate->ProjectionY(Form("hMCPMatched_%s", histName), 1, h2HCorrelate->GetNbinsX(), "e");

  // TH1 *MCPMatchedptSel = (TH1F *) HCorrelate->ProjectionY(Form("hMCPMatchedSel_%s", histName), RBIN, RBIN, HCorrelate->GetZaxis()->FindBin(5), HCorrelate->GetZaxis()->FindBin(100), "e");
  // MCPMatchedptSel = MCPMatchedptSel->Rebin(nptBins, Form("MCPptMatchedSel_%s", histName), ptbin); 
  // TH1 *MCPMatchedptSel = (TH1F *) h2HCorrelate->ProjectionY(Form("hMCPMatchedSel_%s", histName), h2HCorrelate->GetXaxis()->FindBin(5), h2HCorrelate->GetXaxis()->FindBin(100), "e"); 
  // MCPMatchedpt = MCPMatchedpt->Rebin(nptBins, Form("MCPptMatched_%s", histName), ptbin);

  TH2F *Respt = (TH2F *)h2HCorrelate->Clone(Form("hist_%i", ++n));
  TH1F *fake = (TH1F *)JetMCDPt->Clone(Form("hist_%i", ++n));
  fake->Add(MCDMatchedpt, -1);
  TH1F *miss = (TH1F *)JetMCPPt->Clone(Form("hist_%i", ++n));
  miss->Add(MCPMatchedpt, -1);

  RooUnfoldResponse *Response = new RooUnfoldResponse(JetMCDPt, JetMCPPt);
  for (auto i = 1; i <= Respt->GetNbinsX(); i++) {
    for (auto j = 1; j <= Respt->GetNbinsY(); j++) { // ptpair
      Double_t bincenx = Respt->GetXaxis()->GetBinCenter(i);
      Double_t binceny = Respt->GetYaxis()->GetBinCenter(j);
      Double_t bincont = Respt->GetBinContent(i, j);
      for (int k=1; k<=bincont; k++) {
        Response->Fill(bincenx, binceny);
      }
    }
  }
  for (auto i = 1; i <= miss->GetNbinsX(); i++) {
    Double_t bincenx = miss->GetXaxis()->GetBinCenter(i);
    Double_t bincont = miss->GetBinContent(i);
    Response->Miss(bincenx, bincont);
  }
  for (auto i = 1; i <= fake->GetNbinsX(); i++) {
    Double_t bincenx = fake->GetXaxis()->GetBinCenter(i);
    Double_t bincont = fake->GetBinContent(i);
    Response->Fake(bincenx, bincont);
  }

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

  RooUnfoldBayes unfoldCon(Response, JetMCDPt, 4);
  RooUnfoldBayes unfold(Response, DJetpt, 4);
  RooUnfoldBayes unfold2(Response, DJetpt, 5); // systematic study
  RooUnfoldBayes unfold3(Response, DJetpt, 6); // systematic study

  auto hMCcorrected = (TH1F *)unfoldCon.Hreco();
  auto hDcorrected = (TH1F *)unfold.Hreco();
  auto hDcorrected2 = (TH1F *)unfold2.Hreco(); // systematic study
  auto hDcorrected3 = (TH1F *)unfold3.Hreco(); // systematic study

  // int maxK = 20;
  // OptimizeRegularizationParameter(Response, DJetpt, maxK);
  
  RooUnfoldSvd unfoldSVD(Response, DJetpt, 16); // systematic study: SVD Unfolding
  TH1* hDcorrectedSVD = unfoldSVD.Hreco();

  //////////////////////////
  // Jet Momentum Resolution
  //////////////////////////
  TCanvas* canJMR = new TCanvas(Form("JMR_%s", histName), Form("JMR_%s", histName), 900, 800);
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

  ALICEfigureLegend("ALICE Simulation", 0.131403,0.698065,0.380846,0.948387, 0.482183,0.692903,0.732739,0.883871);
  legJMR->SetHeader("Jet Momentum Resolution", "C");
  legJMR->AddEntry("", Form("mean: %.3f", jmrsum / jmrsumweight), "");
  legJMR->Draw();
  
  if(DRAWPLOTS) {canJMR->Print(Form("%s/JetMomentumResolution_%s.pdf", MakeDirName.Data(), histName));}

  /////////////////////////////////////
  /// Jet Reconstruction Efficiency ///
  /////////////////////////////////////
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
  hoptset(*JREp, NevtsMCD / effTrigZvtx10, kBlack, 10, PlotPtMax, 3e-10, 1e-1);
  auto mcpmatchedpt = (TH1 *)MCPMatchedpt->Clone(Form("hist_%i", ++n));
  mcpmatchedpt->GetXaxis()->SetRangeUser(0, PlotPtMax);
  std::cout << "# entries of MCPMatchedpt > 5 GeV: " << mcpmatchedpt->GetEntries() << std::endl;
  JRElegend->AddEntry(mcpmatchedpt, "Matched jets");
  hset(*mcpmatchedpt, JRETitleX, JRETitleY, 0.8, 1.4, 0.04, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hoptset(*mcpmatchedpt, NevtsMCD, kRed, 10, PlotPtMax, 3e-10, 1e-1);
  JREp->GetYaxis()->SetNdivisions(505);
  JREp->Draw("PE");
  mcpmatchedpt->Draw("PEsame");

  ratiojrep->cd();
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

  auto JREcanvas = new TCanvas("JREcanvas", "JREcanvas", 800, 600);
  JREcanvas->cd();
  gPad->SetTicks(1, 1);
  setpad(JREcanvas, 0.02, 0.15, 0.15, 0.05);
  hoptset(*jre, 0, kRed, 5, PlotPtMax, 0.5, 1.01, 1, 1, 2, 20);
  hset(*jre, JRETitleX, "#it{#varepsilon}_{reco}^{jet}", 1.3, 1.0, 0.05, 0.07, 0.01,
       0.01, 0.05, 0.05, 510, 510);
  jre->Draw("pe");
  ALICEfigureLegend("ALICE Simulation", 0.3, 0.2, 0.55, 0.52, 0.6, 0.2, 0.85, 0.44);
  if (DRAWPLOTS) {
    JREcanvas->SaveAs(Form("%s/JetReconstructionEfficiency_%s.pdf", MakeDirName.Data(), histName));
  }


  ////////////////////////////////
  /// Jet Kinematic Efficiency ///
  ////////////////////////////////
  Double_t LcutKine = 7;
  Double_t RcutKine = 140;
    Filipad2 *JKEPad = new Filipad2(++nn, 2, 0.5, 100, 50, 0.7, 1, 1);
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

    auto JKEcanvas = new TCanvas("JKEcanvas", "JKEcanvas", 800, 600);
    JKEcanvas->cd();
    gPad->SetTicks(1, 1);
    setpad(JKEcanvas, 0.02, 0.15, 0.15, 0.05);
    hoptset(*jke, 0, kRed, 5, 200, 0.0, 1.01, 1, 1, 2, 20);
    hset(*jke, JRETitleX, "Kinematic efficiency", 1.3, 1.0, 0.05, 0.07, 0.01,
        0.01, 0.05, 0.05, 510, 510);
    jke->Draw("pe");
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
    
    ALICEfigureLegend("ALICE Simulation", 0.2, 0.46, 0.45, 0.70, 0.2, 0.26, 0.45, 0.44);
    if (DRAWPLOTS) {
      JKEcanvas->SaveAs(Form("%s/JetKinematicEfficiency_%s.pdf", MakeDirName.Data(), histName));
    }

  /////////////////////////////////
  /// Jet Reconstruction Purity ///
  /////////////////////////////////
  jrpp->cd();
  auto JRPp = (TH1F *)JetMCDPt->Clone(Form("hist_%i", ++n));
  JRPlegend->AddEntry("", histName, "");
  JRPlegend->AddEntry(JRPp, "Detector level jets");
  hset(*JRPp, JRPTitleX, JRPTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05,
       0.05, 510, 505);
  hoptset(*JRPp, NevtsMCD, kBlack, 5, PlotPtMax, 5e-9, 1e-3, 1, 1, 1, 20);
  // hoptset(*JRPp, 1, kBlack, 0, PlotPtMax, 1e-12, 1e-3);
  auto mcdmatchedpt = (TH1F *)MCDMatchedpt->Clone(Form("hist_%i", ++n));
  JRPlegend->AddEntry(mcdmatchedpt, "Matched jets in Detector level");
  hset(*mcdmatchedpt, JRPTitleX, JRPTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hoptset(*mcdmatchedpt, NevtsMCD, kRed, 5, PlotPtMax, 3e-10, 1e-3, 1, 1, 1, 24);
  // hoptset(*mcdmatchedpt, 1, kRed, 0, PlotPtMax, 1e-12, 1e-3);
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

  auto JRPcanvas = new TCanvas("JRPcanvas", "JRPcanvas", 800, 600);
  JRPcanvas->cd();
  gPad->SetTicks(1, 1);
  setpad(JRPcanvas, 0.02, 0.15, 0.15, 0.05);
  hoptset(*jrp, 0, kRed, 5, PlotPtMax, 0., 1.01, 1, 1, 2, 20);
  hset(*jrp, JRPTitleX, "Jet purity", 1.3, 1.0, 0.05, 0.07, 0.01,
      0.01, 0.05, 0.05, 510, 510);
  jrp->Draw("pe");
  ALICEfigureLegend("ALICE Simulation", 0.3, 0.2, 0.55, 0.52, 0.6, 0.2, 0.85, 0.44);
  if (DRAWPLOTS) {
    JRPcanvas->SaveAs(Form("%s/JetReconstructionPurity_%s.pdf", MakeDirName.Data(), histName));
  }

  ////////////////////////////////////////////////
  /// 2D pT correlation plot / Response Matrix ///
  ////////////////////////////////////////////////
  TCanvas *canResponseMatrix = new TCanvas(Form("Correlation_%s", histName),
                                        Form("c%s", histName), 1100, 1100);
  TLegend *legRM =
          new TLegend(0.0510018,0.735814,0.397086,0.995349,NULL,"brNDC");
  canResponseMatrix->cd();
  setpad(canResponseMatrix, 0.27, 0.12, 0.135, 0.2);
  // canResponseMatrix->SetLogx(1);
  // canResponseMatrix->SetLogy(1);
  canResponseMatrix->SetLogz(1);
  auto hRM = (TH2 *) hResponseMatrix->Clone("hRM");
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

  auto XYline = new TF1("XYline", "x", 0, 300);
  XYline->SetLineColor(kBlack);
  XYline->SetLineStyle(2);
  XYline->Draw("same");

  legRM->SetTextSize(0.04);
  legRM->SetBorderSize(0);
  legRM->SetTextAlign(12);
  legRM->SetFillColorAlpha(0,0); 
  legRM->AddEntry("", "ALICE Simulation", "");
  legRM->AddEntry("", "pp #sqrt{#it{s}} = 13.6 TeV, Response matrix", "");
  legRM->AddEntry("", "#it{p}_{T, track} > 0.15 GeV/#it{c}", "");
  legRM->AddEntry("", "|#it{#eta}_{track}| < 0.9, |#it{#eta}_{jet}| < 0.5", "");
  legRM->AddEntry("", "Anti-#it{k}_{T}, charged-particle jet, #it{R} = 0.4", "");
  legRM->Draw();

  canResponseMatrix->Update();
  if(DRAWRM) {canResponseMatrix->Print(Form("%s/JetPtCorrelation_R%.1f_%s.pdf",
                             MakeDirName.Data(), RBIN, histName));}

  /////////////////////////
  /// Consistency Check ///
  /////////////////////////
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

  auto lumiCounterFile = TFile::Open(DataCollCounterFile, "open");
  TH1* hLumiCounter = (TH1 *) lumiCounterFile->Get(LumiObj);
  auto NBCTVX = (double) hLumiCounter->GetBinContent(hLumiCounter->FindBin(1));
  auto NCollTVX = (double) hLumiCounter->GetBinContent(hLumiCounter->FindBin(5));
  auto NCollTVXSel = (double) hLumiCounter->GetBinContent(hLumiCounter->FindBin(6));
  double NBCsel = NBCTVX * (NCollTVXSel / NCollTVX);
  double luminositySel8 = NBCsel / SigmaTVX * P_mu;
  std::cout << "luminositySel8: " << luminositySel8 << " /mb" << std::endl;

  TH1 *hLumiTVX = (TH1 *) Dfile->Get(Form("%s", "bc-selection-task/hLumiTVX"));
  TH1 *hLumiTVXafBCcuts = (TH1 *) Dfile->Get(Form("%s", "bc-selection-task/hLumiTVXafterBCcuts"));
  auto lumiTVX = (Double_t) effTrigZvtx10 * 1000 * hLumiTVX->Integral(1, hLumiTVX->GetNbinsX());
  auto lumiTVXafBCcuts = (Double_t)  effTrigZvtx10 * 1000 * hLumiTVXafBCcuts->Integral(1, hLumiTVXafBCcuts->GetNbinsX()); // TVXafBCcuts doesn't have zvtx selection.
  std::cout << "lumiTVX: " << lumiTVX << " /mb" << std::endl;
  std::cout << "lumiTVXafBCcuts: " << lumiTVXafBCcuts << " /mb" << std::endl;

  TH1 *hMCLumiTVX = (TH1 *) MCfile->Get(Form("%s", "bc-selection-task/hLumiTVX"));
  auto MClumiTVX = (Double_t) 1000 * hMCLumiTVX->Integral(1, hMCLumiTVX->GetNbinsX());

  TH1 *hNtvx = (TH1 *) Dfile->Get(Form("%s", NTVXObj));
  Double_t Ntvx = (Double_t) hNtvx->Integral(1, hNtvx->GetNbinsX());
  std::cout << "N_{TVX}: " << Ntvx << std::endl;
  Double_t luminosity = effTrigZvtx10 * Ntvx / SigmaTVX * P_mu;
  std::cout << "luminosity_{int}: " << luminosity << " /mb" << std::endl;

  auto mccollCounterfileSel8 = TFile::Open(DataCollCounterFile, "open");
  // TH1 *DataCollCounter = (TH1 *) mccollCounterfileSel8->Get("track-efficiency/hMcCollCutsCounts");
  TH1 *DataCollCounter = (TH1 *) mccollCounterfileSel8->Get(Form("%s/%s", DataDirectory[0].Data(), EventObj));
  // Double_t NDatacollINEL = DataCollCounter->GetBinContent(DataCollCounter->FindBin(1.5));
  // Double_t NDatacollINEL = DataCollCounter->GetBinContent(DataCollCounter->FindBin(3.5)); // non-splited mcCollisions w/o EvSel
  Double_t NDatacollINEL = DataCollCounter->GetBinContent(DataCollCounter->FindBin(0.5));
  cout << "NDatacollINEL: " << NDatacollINEL << endl;
  // Double_t NDatacollsel8 = DataCollCounter->GetBinContent(DataCollCounter->FindBin(4.5));
  Double_t NDatacollsel8 = DataCollCounter->GetBinContent(DataCollCounter->FindBin(1.5));
  Double_t effINELsel8 = NDatacollsel8 / NDatacollINEL;
  cout << "eff(sel8 / INEL): " << effINELsel8 << endl;

  auto mccollCounterfile = TFile::Open(McCollCounterFile, "open");
  TH1 *McCollCounter = (TH1 *) mccollCounterfile->Get("track-efficiency/hMcCollCutsCounts");
  Double_t NmccollINEL = McCollCounter->GetBinContent(McCollCounter->FindBin(1.5));
  Double_t NmccollTVX = McCollCounter->GetBinContent(McCollCounter->FindBin(2.5));
  cout << "NmccollINEL: " << NmccollINEL << endl;
  Double_t NmccollTVXnonSplit = McCollCounter->GetBinContent(McCollCounter->FindBin(3.5));
  Double_t effSplit = NmccollTVXnonSplit / NmccollTVX;
  cout << "effSplit: " << effSplit << endl;
  Double_t NmccollselMC = McCollCounter->GetBinContent(McCollCounter->FindBin(4.5));
  cout << "NmccollselMC: " << NmccollselMC << endl;
  Double_t effINELselMC = NmccollselMC / NmccollTVX;
  cout << "eff(selMC / INEL): " << effINELselMC << endl;

  // new TCanvas("hTrigEff", "hTrigEff", 800, 600);
  // auto hTriggerEfficiency = GetTriggerEfficiency();
  // hTriggerEfficiency = hTriggerEfficiency->Rebin(nptBinsGen, "hTriggerEfficiency", ptbinGen);
  // hTriggerEfficiency->Scale(1., "width");
  // hTriggerEfficiency->Draw("pe");

  ////////////////////////////////////////////////
  ///// Run 3 corrected Data and MC Invariant Yield /////
  ////////////////////////////////////////////////
  auto InvYData = (TH1F *)hDcorrected->Clone(Form("hist_%i", ++n));
  hoptset(*InvYData, NCollTVXSel / effSplit / effTrigZvtx10, ColorPallete[1], PlotPtMin, PlotPtMax, 6e-9, 5e-3, 0.6, 1, 2, 24);

  auto InvYMC = (TH1*) JetMCPPt->Clone(Form("hist_%i", ++n));
  hoptset(*InvYMC, NmccollselMC / effTrigZvtx10, ColorPallete[0], PlotPtMin, PlotPtMax, 3e-10, 1e-1, 0.75, 1, 2, 33);

  //////////////////////////////////////////////
  ///// Run 3 corrected Data Cross Section /////
  //////////////////////////////////////////////
  auto UnfoldData = (TH1 *)hDcorrected->Clone("Run3_CrossSection");
  UnfoldData->Scale(1/luminositySel8 * effTrigZvtx10, "width");
  hset(*UnfoldData, JetPtDataFinalTitleX, XSECTION ? XSectionTitleY : JetPtDataFinalTitleY, 0.9, 1.5, 0.05, 0.045, 0.01, 0.005, 0.045, 0.040, 510, 1005);
  hoptset(*UnfoldData, 0, ColorPallete[1], PlotPtMin, PlotPtMax, XSECTION? 3e-7 : 1e-9, XSECTION? 5e-1 : 1e-1, 0.6, 1, 2, 24);
  // UnfoldData->SaveAs("Run3_CrossSection_wTrackTuner.root");


  // if (SYSTUNFOLD) {
  /////////////////////////////////////////////////////////
  ///// Systematic Study on Tracking Efficiency       /////
  /////////////////////////////////////////////////////////

  Filipad2 *SystTrkEffPad = new Filipad2(++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
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

  Filipad2 *SystTrkResPad = new Filipad2(++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
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
  legsystTrkRes->AddEntry(defaultWtrackTuner, "w/ track tuner (default)", "lpe");
  hset(*defaultWtrackTuner, JetPtTitleX, XSectionTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hoptset(*defaultWtrackTuner, 0, kBlack, PlotPtMin, PlotPtMax, 3e-7, 1e1, 0.7, 1, 2, 24);
  defaultWtrackTuner->Draw("pe");

  auto UnfoldWoTrackTuner = (TH1 *) Run3XSectionWoTrackTuner();
  legsystTrkRes->AddEntry(UnfoldWoTrackTuner, "w/o track tuner", "lpe");
  // hset(*UnfoldWoTrackTuner, JetPtTitleX, XSectionTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
  //      0.05, 0.05, 510, 505);
  hoptset(*UnfoldWoTrackTuner, 0, kRed, PlotPtMin, PlotPtMax, 1e-9, 1e-1, 0.7, 1, 2, 25);
  UnfoldWoTrackTuner->Draw("pesame");

  legsystTrkRes->Draw();

  systTrkResRatpad->cd();
  auto hRatWoTrkTuner = (TH1 *) UnfoldWoTrackTuner->Clone(Form("hist_%i", ++n));
  hRatWoTrkTuner->Divide(hRatWoTrkTuner, defaultWtrackTuner, 1, 1, "B");
  hset(*hRatWoTrkTuner, JetPtGenTitleX, "variation / default", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
       0.1, 0.1, 510, 505);
  hoptset(*hRatWoTrkTuner, 0, kRed, PlotPtMin, PlotPtMax, 0.95, 1.55, 0.7, 1, 2, 25);
  hRatWoTrkTuner->Draw("pe");
  OutStatsTXT(hRatWoTrkTuner, "TrackPtResolution");

  if(SYSTUNFOLD) {SystTrkResPad->C->SaveAs(Form("%s/SystErrTrackPtResolution.pdf", MakeDirName.Data()));} 

  ///////////////////////////////////////////////////
  ///// Systematic Study on Unfolding Iteration /////
  ///////////////////////////////////////////////////

  Filipad2 *SystIterPad = new Filipad2(++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
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
  legsystiter->AddEntry(defaultUnfold, "Data Unfolded (Iter: 4)", "lpe");
  hset(*defaultUnfold, JetPtTitleX, JetPtDataFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hoptset(*defaultUnfold, NevtsData / effINELselMC, kBlack, PlotPtMin, PlotPtMax, 1e-9, 1e-1, 0.7, 1, 2, 24);
  defaultUnfold->Draw("pe");

  auto Unfold2 = (TH1 *)hDcorrected2->Clone(Form("hist_%i", ++n));
  legsystiter->AddEntry(Unfold2, "Data Unfolded (Iter: 5)", "lpe");
  hset(*Unfold2, JetPtTitleX, JetPtDataFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hoptset(*Unfold2, NevtsData / effINELselMC, kRed, PlotPtMin, PlotPtMax, 1e-9, 1e-1, 0.7, 1, 2, 27);
  Unfold2->Draw("pesame");

  auto Unfold3 = (TH1 *)hDcorrected3->Clone(Form("hist_%i", ++n));
  legsystiter->AddEntry(Unfold3, "Data Unfolded (Iter: 6)", "lpe");
  hset(*Unfold3, JetPtTitleX, JetPtDataFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hoptset(*Unfold3, NevtsData / effINELselMC, kBlue, PlotPtMin, PlotPtMax, 1e-9, 1e-1, 0.7, 1, 2, 28);
  Unfold3->Draw("pesame");

  legsystiter->Draw();

  systiterratpad->cd();
  auto Unfoldrat2 = (TH1 *) Unfold2->Clone(Form("hist_%i", ++n));
  Unfoldrat2->Divide(Unfoldrat2, defaultUnfold, 1, 1, "B");
  hset(*Unfoldrat2, JetPtGenTitleX, "variation / default", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
       0.1, 0.1, 510, 505);
  hoptset(*Unfoldrat2, 0, kRed, PlotPtMin, PlotPtMax, 0.92, 1.08, 0.7, 1, 2, 27);
  Unfoldrat2->Draw("pe");
  OutStatsTXT(Unfoldrat2, "UnfoldIter_5");

  auto Unfoldrat3 = (TH1 *) Unfold3->Clone(Form("hist_%i", ++n));
  Unfoldrat3->Divide(Unfoldrat3, defaultUnfold, 1, 1, "B");
  hset(*Unfoldrat3, JetPtGenTitleX, "variation / default", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
       0.1, 0.1, 510, 505);
  hoptset(*Unfoldrat3, 0, kBlue, PlotPtMin, PlotPtMax, 0.92, 1.08, 0.7, 1, 2, 28);
  Unfoldrat3->Draw("pesame");
  OutStatsTXT(Unfoldrat3, "UnfoldIter_6");

  if(SYSTUNFOLD) {SystIterPad->C->SaveAs(Form("%s/SystErrUnfoldingIteration.pdf", MakeDirName.Data()));}

  ////////////////////////////////////////////////
  ///// Systematic Study on Unfolding Method /////
  ////////////////////////////////////////////////
  Filipad2 *SystMethodPad = new Filipad2(++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
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
  hoptset(*hBayesUnfold, NevtsData / effINELselMC, kBlack, PlotPtMin, PlotPtMax, 3e-10, 1e-1, 0.7, 1, 2, 24);
  hBayesUnfold->Draw("pe");

  // Plot SVD Unfolding result
  auto hSVDUnfold = (TH1F *)hDcorrectedSVD->Clone(Form("hist_%i", ++n));
  legsystmethod->AddEntry(hSVDUnfold, "SVD unfold");
  hset(*hSVDUnfold, JetPtTitleX, JetPtDataFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
        0.05, 0.05, 510, 505);
  hoptset(*hSVDUnfold, NevtsData / effINELselMC, kBlue, PlotPtMin, PlotPtMax, 3e-10, 1e-1, 0.7, 1, 2, 27);
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
  auto hUnfoldedDataHerwig = DrawUnfoldHerwig(InvYMC, JetMCPPt, JetMCDPt, h2HCorrelate, DJetpt, NevtsMCD / effINELselMC);

  Filipad2 *SystPriorPad = new Filipad2(++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
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
  hoptset(*UnfoldHerwig, NevtsData / effINELselMC, ColorPallete[3], PlotPtMin, PlotPtMax, 3e-10, 1e-1, 0.7, 1, 2, 27);
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
  TH1 *hSystMcClosure = new TH1F("hSystMcClosure", "hSystMcClosure", nptBinsGen, ptbinGen);
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
  TH1 *hSystNorm = new TH1F("hSystNorm", "hSystNorm", nptBinsGen, ptbinGen);
  for (int i=0; i<=hSystNorm->GetNbinsX(); i++) {
    hSystNorm->SetBinContent(i, SystNorm);
  }
  hSysts.push_back(hSystNorm);

  double SystSecCon = 1.05;
  TH1 *hSystSecCon = new TH1F("hSystSecCon", "hSystSecCon", nptBinsGen, ptbinGen);
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

  TH1* hSystResult = (TH1*)hSysts[0]->Clone("quadratureSum");
    hSystResult->Reset();
  calculateQuadratureSum(hSysts, hSystResult);
  OutStatsTXT(hSystResult, "TotalSystematicUncertainty");

  DrawMultipleSources(hSysts, hSystResult);

  ///////////////////////////////////////////////////////////
  ///// Run 3 corrected Invariant Yield: Data vs gen MC /////
  ///////////////////////////////////////////////////////////
  Filipad2 *InvariantYieldPad = new Filipad2(++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
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
  legendYield3->AddEntry("", "ALICE Preliminary", "");
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
  UnfoldDataWsyst->Draw("E2same");
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
  unfoldlegend->AddEntry(NormMCP, SEL8WINDOW? "PYTHIA8" : "MC true", "lep");
  // hoptset(*NormMCP, XSECTION ? ((NmccollINEL / SigmaINEL) * effINELselMC) : (SEL8WINDOW? NevtsMCD / effINELselMC : 1 /* TBC */), ColorPallete[0], PlotPtMin, PlotPtMax, XSECTION? 1e-7 : 1e-9, XSECTION? 1e1 : 1e-1, 0.5, 1, 2, 24);
  hoptset(*NormMCP, XSECTION ? (NmccollselMC / effTrigZvtx10 / SigmaTVX) : (SEL8WINDOW? NevtsMCD / effINELselMC : 1 /* TBC */), ColorPallete[0], PlotPtMin, PlotPtMax, XSECTION? 1e-7 : 1e-9, XSECTION? 1e1 : 1e-1, 0.75, 1, 2, 33);
  // hoptset(*NormMCP, 0, kGreen, PlotPtMin, PlotPtMax, XSECTION? 1e-7 : 1e-9, XSECTION? 1e1 : 1e-1);
  NormMCP->Draw("pesame");

  // TH1 *PYTHIA1360 = GetPYTHIA1360();
  // unfoldlegend->AddEntry(PYTHIA1360, "PYTHIA8 13.6 TeV", "lep");
  // hoptset(*PYTHIA1360, 0, ColorPallete[4], PlotPtMin, PlotPtMax,  1e-7, 1e1, 0.5, 1, 2, 24);
  // hset(*PYTHIA1360, JetPtDataFinalTitleX, JetPtMCFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01);
  // PYTHIA1360->Draw("pesame");

  // TH1 *PYTHIA1300 = GetPYTHIA1300();
  // unfoldlegend->AddEntry(PYTHIA1300, "PYTHIA8 13.0 TeV", "lep");
  // hoptset(*PYTHIA1300, 0, ColorPallete[5], PlotPtMin, PlotPtMax,  1e-7, 1e1, 0.5, 1, 2, 24);
  // hset(*PYTHIA1300, JetPtDataFinalTitleX, JetPtMCFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01);
  // if(XSECTION) {PYTHIA1300->Scale(77.23);}
  // PYTHIA1300->Draw("pesame");

  TH1 *HERWIG1360 = GetHERWIG1360();
  unfoldlegend->AddEntry(HERWIG1360, "Herwig", "lep");
  hoptset(*HERWIG1360, 0, ColorPallete[3], PlotPtMin, PlotPtMax,  1e-7, 1e1, 0.7, 1, 2, 34);
  hset(*HERWIG1360, JetPtDataFinalTitleX, JetPtMCFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01);
  // if(XSECTION) {HERWIG1360->Scale(SigmaINEL);}
  HERWIG1360->Draw("pesame");


  ratunfoldp->cd();
  gPad->SetTicks(1, 1);
  
  // auto ratPYTHIA1300 = DrawRatioTH1(PYTHIA1300, UnfoldData);
  // hset(*ratPYTHIA1300, JetPtGenTitleX, "Comp. / Run 3", 1.2, 0.75, 0.1, 0.09, 0.01, 0.01,
  //      0.1, 0.1, 510, 505);
  // hoptset(*ratPYTHIA1300, 0, ColorPallete[5], PlotPtMin, PlotPtMax, 0.4, 2.7, 0.4, 1, 2, 20);
  // ratPYTHIA1300->Draw("pe");

  // auto ratMCTVXData = (TH1 *)NormMCP->Clone(Form("hist_%i", ++n));
  // ratMCTVXData->Divide(ratMCTVXData, UnfoldData, 1., 1., "B");
  auto ratMCTVXData = DrawRatioTH1(NormMCP, UnfoldData);
  hset(*ratMCTVXData, JetPtDataFinalTitleX, "Comp. / Run 3 Data", 1.19, 0.76, 0.1, 0.085, 0.01, 0.01,
      0.1, 0.1, 510, 505);
  hoptset(*ratMCTVXData, 0, ColorPallete[0], PlotPtMin, PlotPtMax, 0.65, 2.8, 0.65, 1, 2, 33);
  ratMCTVXData->Draw("pesame");

  // auto ratMCtrueData = (TH1F *)NormMCPINEL->Clone(Form("hist_%i", ++n));
  // ratMCtrueData->Divide(ratMCtrueData, UnfoldData, 1., 1., "B");
  // hset(*ratMCtrueData, JetPtDataFinalTitleX, "Comp. / Run 3", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
  //     0.1, 0.1, 510, 505);
  // hoptset(*ratMCtrueData, 0, ColorPallete[0], PlotPtMin, PlotPtMax, 0.4, 2.7, 0.4, 1, 2, 24);
  // ratMCtrueData->Draw("pesame");

  // auto ratPYTHIA1360 = DrawRatioTH1(PYTHIA1360, UnfoldData);
  // hset(*ratPYTHIA1360, JetPtGenTitleX, "Comp. / Run 3", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
  //      0.1, 0.1, 510, 505);
  // hoptset(*ratPYTHIA1360, 0, ColorPallete[4], PlotPtMin, PlotPtMax, 0.4, 2.7, 0.4, 1, 2, 20);
  // ratPYTHIA1360->Draw("pesame");

  auto UnfoldDataSelf = (TH1 *) UnfoldData->Clone(Form("hist_%i", ++n));
  UnfoldDataSelf->Divide(UnfoldDataSelf, UnfoldData, 1, 1, "B");
  TH1 *UnfoldDataSystRatio = ApplySystematicUncertainty(UnfoldDataSelf, hSystResult);
  UnfoldDataSystRatio->SetMarkerStyle(0);
  UnfoldDataSystRatio->Draw("E2same");

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
  
  TH1 *ratHerwigRun3 = DrawRatioTH1(HERWIG1360, UnfoldData);
  // TH1 *ratHerwigRun3 = DrawRatioTH1(HERWIG1360, hrun2Syst);
  hset(*ratHerwigRun3, JetPtDataFinalTitleX, "Herwig / PYTHIA", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
       0.1, 0.1, 510, 505);
  hoptset(*ratHerwigRun3, 0, ColorPallete[3], PlotPtMin, PlotPtMax, 0.4, 2.7, 0.6, 1, 2, 34);
  ratHerwigRun3->Draw("pesame");

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
//   hoptset(*run2MCgen, 0, kViolet, PlotPtMin, PlotPtMax, 1e-7, 1e1);
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
//   hoptset(*run2MCratio, 0, kViolet, PlotPtMin, PlotPtMax, 0, 2);
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

TH1 *DrawJetResolution(const char *fileName, const char *histName, Double_t Nevt,
                       TLegend *legend, Color_t colorID,
                       const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  TString histNameStr(histName);

  TH3 *JetResolution = nullptr;
  TH3 *temp1 = (TH3 *)file->Get(Form("%s/%s", Dir, JetResolutionObjOld));
  TH3 *temp2 = (TH3 *)file->Get(Form("%s/%s", Dir, JetResolutionObj));

  if (temp1) {
    JetResolution = temp1;
  } else if (temp2) {
    JetResolution = temp2;
  }
  if (!JetResolution) {
    std::cout << "Error: No valid JetResolution object found." << std::endl;
  }

  std::vector<std::pair<int, int>> pTClasses = {{10, 15}, {15, 20},
  {20, 40}, {40, 60}};
  // std::vector<std::pair<int, int>> pTClasses = {
      // {5, 20}, {20, 40}, {40, 60}, {60, 85}};

  for (Int_t i = 0; i < pTClasses.size(); i++) {
    TH1D *jetResolution = (TH1D *)JetResolution->ProjectionZ(
        Form("JetResZ%i_%s", i + 1, histName),
        JetResolution->GetXaxis()->FindBin(RBIN - 1e-6),
        JetResolution->GetXaxis()->FindBin(RBIN + 1e-6),
        JetResolution->GetYaxis()->FindBin(pTClasses[i].first),
        JetResolution->GetYaxis()->FindBin(pTClasses[i].second));

    hset(*jetResolution, JetResolutionTitleX, "Probability density", 1.10, 1.2, 0.047, 0.05,
         0.001, 0.001, 0.05, 0.05);
    legend->AddEntry(jetResolution,
                     Form("#it{p}_{T, jet}^{true} [%i, %i] GeV/#it{c}, #mu: %.2f, #sigma: %.2f",
                          pTClasses[i].first, pTClasses[i].second, jetResolution->GetMean(), jetResolution->GetRMS()),
                     "pe");
    // jetResolution->GetMean(), jetResolution->GetRMS()), "pe2");
    // hset(*jetResolution, JetResolutionTitleX, JetPtTitleY, 0.7, 1.0, 0.05,
    // 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*jetResolution, 2, ColorPallete[i], -0.5, 1., 1e-5, 1e1);
    jetResolution->Draw("esame");
  }

  return 0;
}

