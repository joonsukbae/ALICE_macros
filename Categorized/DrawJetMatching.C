#include "common.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include "DrawJetsMC.h"

Int_t n = 0;
Int_t nn = 0;

TH1F *DrawJetMatching(const char *fileName, const char *DfileName,
                      const char *histName, float Nevts, TPad *jrep,
                      TPad *ratiojrep, TLegend *JRElegend, TPad *consistencyp,
                      TPad *ratioconsistp, TLegend *CONSISTlegend,
                      TPad *unfoldp, TPad *ratunfoldp, TLegend *unfoldlegend,
                      Color_t colorID, TPad *jrpp = nullptr,
                      TPad *ratiojrpp = nullptr, TLegend *JRPlegend = nullptr,
                      TFile *savefile = nullptr) {

  auto Dfile = TFile::Open(DfileName, "open");
  auto Ddir = (TDirectory *)Dfile->Get(Dir);
  Ddir->cd();
  TH1 *DJetpt = (TH1 *)gROOT->FindObject("h_jet_pt");
  cout << "Nbins DJetpt: " << DJetpt->GetNbinsX() << endl;

  auto file = TFile::Open(fileName, "open");
  auto dir = (TDirectory *)file->Get(Dir);
  dir->cd();
  TH1 *JetMCPPt = (TH1 *)file->Get(Form("%s/%s", Dir, "h_jet_pt_part"));
  TH1 *JetMCDPt = (TH1 *)file->Get(Form("%s/%s", Dir, "h_jet_pt"));
  TString histNameStr(histName);
  TH3 *HCorrelate = nullptr;
  TH3 *MatchTemp1 = (TH3 *)file->Get(
      Form("%s/%s", Dir, "h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"));
  TH3 *MatchTemp2 =
      (TH3 *)file->Get(Form("%s/%s", Dir, "h3_jet_r_jet_pt_part_jet_pt"));
  if (MatchTemp1) {
    HCorrelate = MatchTemp1;
  } else if (MatchTemp2) {
    HCorrelate = MatchTemp2;
  }
  HCorrelate->GetXaxis()->SetRange(RBIN, RBIN);
  auto h2HCorrelate = (TH2 *)HCorrelate->Project3D(Form("%s_yze", histName));

  TH2F *hcorrelate = new TH2F(Form("hcorrelate_%s", histName),
                              Form("R projected correlate_%s", histName),
                              nptBins, ptbin, nptBins, ptbin);
  Int_t corrbin = HCorrelate->GetXaxis()->FindBin(RBIN + 1e-6);
  for (Int_t i = 0; i <= HCorrelate->GetNbinsY(); i++) {
    for (Int_t j = 0; j <= HCorrelate->GetNbinsZ(); j++) {
      Double_t content = HCorrelate->GetBinContent(corrbin, i, j);
      Double_t error = HCorrelate->GetBinError(corrbin, i, j);

      Int_t binpart = hcorrelate->GetYaxis()->FindBin(
          HCorrelate->GetYaxis()->GetBinCenter(i));
      Int_t bin = hcorrelate->GetXaxis()->FindBin(
          HCorrelate->GetZaxis()->GetBinCenter(j));

      Double_t currentContent = hcorrelate->GetBinContent(bin, binpart);
      Double_t currentError = hcorrelate->GetBinError(bin, binpart);

      Double_t newContent = currentContent + content;
      Double_t newError = sqrt(pow(currentError, 2) +
                               pow(error, 2));

      hcorrelate->SetBinContent(bin, binpart, newContent);
      hcorrelate->SetBinError(bin, binpart, newError);
    }
  }

  TH1 *MCDMatchedpt = (TH1F *)HCorrelate->ProjectionZ(
      Form("hMCDMatched_%s", histName),
      HCorrelate->GetXaxis()->FindBin(RBIN + 1e-6),
      HCorrelate->GetXaxis()->FindBin(RBIN + 0.2 - 1e-6), 1,
      HCorrelate->GetNbinsY(), "e");
  TH1 *MCPMatchedpt =
      (TH1F *)HCorrelate->ProjectionY(Form("hMCPMatched_%s", histName), RBIN,
                                      RBIN, 1, HCorrelate->GetNbinsZ(), "e");

  TH2F *Respt = (TH2F *)h2HCorrelate->Clone();
  TH1F *fake = (TH1F *)JetMCDPt->Clone();
  fake->Add(MCDMatchedpt, -1);
  TH1F *miss = (TH1F *)JetMCPPt->Clone();
  miss->Add(MCPMatchedpt, -1);

  RooUnfoldResponse *Response = new RooUnfoldResponse(JetMCDPt, JetMCPPt);
  for (auto i = 1; i <= Respt->GetNbinsX(); i++) {
    for (auto j = 1; j <= Respt->GetNbinsY(); j++) {
      Double_t bincenx = Respt->GetXaxis()->GetBinCenter(i);
      Double_t binceny = Respt->GetYaxis()->GetBinCenter(j);
      Double_t bincont = Respt->GetBinContent(i, j);
      Response->Fill(bincenx, binceny, bincont);
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
  TH1 *MCPMatchedptSel = (TH1F *)hResponseMatrix->ProjectionY(
      Form("hMCPMatchedSel_%s", histName),
      hResponseMatrix->GetXaxis()->FindBin(5),
      hResponseMatrix->GetXaxis()->FindBin(100), "e");
  MCPMatchedptSel = MCPMatchedptSel->Rebin(
      nptBins, Form("MCPMatchedptSelRebin_%s", histName), ptbin);

  TCanvas *canhResponseMatrix =
      new TCanvas(Form("hResponseMatrix_%s", histName),
                  Form("hResponseMatrix_%s", histName), 800, 800);
  canhResponseMatrix->cd();
  canhResponseMatrix->SetLogz(1);

  auto hNormResponseMatrix = (TH2 *)hResponseMatrix->Clone();
  hNormResponseMatrix->Scale(1. / hNormResponseMatrix->Integral(), "width");
  hNormResponseMatrix->GetZaxis()->SetRangeUser(1e-10, 1e0);
  hNormResponseMatrix->Draw("colz");

  RooUnfoldBayes unfoldCon(Response, JetMCDPt, 4);
  RooUnfoldBayes unfold(Response, DJetpt, 4);
  auto hMCcorrected = (TH1F *)unfoldCon.Hreco();
  auto hDcorrected = (TH1F *)unfold.Hreco();

  jrep->cd();
  TH1 *rec_total_window =
      (TH1 *)hResponseMatrix->ProjectionY("rec_total_window", 0, -1);
  TH1 *hmissed = (TH1 *)Response->Miss()->Clone("missed_events");
  rec_total_window->Add(hmissed, 1);
  rec_total_window =
      rec_total_window->Rebin(nptBins, Form("JREpRebin_%s", histName), ptbin);
  TH1 *JREp = (TH1F *)rec_total_window->Clone();
  JRElegend->AddEntry("", histName, "");
  JRElegend->AddEntry(JREp, "Total #it{p}_{T, jet}^{reco} window");
  hset(*JREp, JRETitleX, JRETitleY, 0.9, 1., 0.05, 0.07, 0.01, 0.01, 0.05, 0.05,
       510, 510);
  hoptset(*JREp, Nevts, kBlack, 0 + 1e-6, 200, 5e-4, 1e-1);
  auto mcpmatchedpt = (TH1F *)MCPMatchedptSel->Clone();
  JRElegend->AddEntry(mcpmatchedpt, "Selected #it{p}_{T, jet}^{reco} window");
  hset(*mcpmatchedpt, JRETitleX, JRETitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 510);
  hoptset(*mcpmatchedpt, Nevts, kRed, 0 + 1e-6, 200, 1e-4, 1);
  JREp->GetYaxis()->SetNdivisions(505);
  JREp->Draw("pe");
  mcpmatchedpt->Draw("pesame");
  ratiojrep->cd();
  auto jre = (TH1F *)MCPMatchedptSel->Clone();
  auto KinJetMCPpt = (TH1F *)MCPMatchedpt->Clone();
  jre->Divide(jre, rec_total_window, 1., 1., "B");
  hset(*jre, JRETitleX, "Kinematic efficiency", 1.2, 1.0, 0.06, 0.07, 0.01,
       0.01, 0.07, 0.07, 510, 510);
  jre->SetMarkerColor(kRed);
  jre->SetLineColor(kRed);
  jre->SetMarkerSize(.7);
  jre->SetMarkerStyle(22);
  jre->GetXaxis()->SetRangeUser(0. + 1e-6, 200.);
  jre->GetYaxis()->SetRangeUser(0., 1.1);
  jre->SetFillColorAlpha(kRed, 0.3);
  jre->Draw("pe");

  jrpp->cd();
  auto JRPp = (TH1F *)JetMCDPt->Clone();
  JRPlegend->AddEntry("", histName, "");
  JRPlegend->AddEntry(JRPp, "Detector level jets");
  hset(*JRPp, JRPTitleX, JRPTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05,
       0.05, 510, 505);
  hoptset(*JRPp, Nevts, kBlack, 0, 200, 1e-12, 1e-3);
  auto mcdmatchedpt = (TH1F *)MCDMatchedpt->Clone();
  JRPlegend->AddEntry(mcdmatchedpt, "Matched jets in Detector level");
  hset(*mcdmatchedpt, JRPTitleX, JRPTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hoptset(*mcdmatchedpt, Nevts, kRed, 0, 200, 1e-12, 1e-3);
  JRPp->GetYaxis()->SetNdivisions(505);
  JRPp->Draw("pe");
  mcdmatchedpt->Draw("pesame");
  ratiojrpp->cd();
  auto jrp = (TH1F *)MCDMatchedpt->Clone();
  jrp->Divide(jrp, JetMCDPt, 1., 1., "B");
  hset(*jrp, JRPTitleX, "Jet Purity", 1.2, 1.0, 0.07, 0.07, 0.01, 0.01, 0.07,
       0.07, 510, 505);
  jrp->SetMarkerColor(kRed);
  jrp->SetLineColor(kRed);
  jrp->SetMarkerSize(.7);
  jrp->SetMarkerStyle(22);
  jrp->GetXaxis()->SetRangeUser(0., 200.);
  jrp->GetYaxis()->SetRangeUser(0., 1.);
  jrp->SetFillColorAlpha(kRed, 0.3);
  jrp->Draw("pe");

  TCanvas *canCorrelation = new TCanvas(Form("Correlation_%s", histName),
                                        Form("c%s", histName), 800, 800);
  canCorrelation->cd();
  canCorrelation->SetLogz(1);
  hcorrelate->GetXaxis()->SetTitleSize(0.033);
  hcorrelate->GetYaxis()->SetTitleSize(0.033);
  hcorrelate->GetXaxis()->SetLabelSize(0.025);
  hcorrelate->GetYaxis()->SetLabelSize(0.025);
  hcorrelate->GetZaxis()->SetLabelSize(0.025);
  hcorrelate->SetTitle(Form("%s", histName));
  hcorrelate->GetXaxis()->SetTitle("#it{p}_{T, jet}^{reco} (GeV/c)");
  hcorrelate->GetYaxis()->SetTitle("#it{p}_{T, jet}^{truth} (GeV/c)");
  hcorrelate->Scale(1. / hcorrelate->Integral(), "width");
  hcorrelate->GetZaxis()->SetRangeUser(1e-10, 1e0);
  hcorrelate->Draw("colz");
  canCorrelation->Print(Form(
      "plots/JetPtCorrelation_R%.1f_%s.pdf", RBIN, histName));

  consistencyp->cd();
  auto rawMCP = (TH1F *)JetMCPPt->Clone();
  CONSISTlegend->AddEntry(rawMCP, "MC Generated");
  hset(*rawMCP, JRETitleX, JRETitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05,
       0.05, 510, 505);
  rawMCP->Scale(1. / rawMCP->Integral(), "width");
  rawMCP->SetMarkerSize(1.3);
  rawMCP->SetMarkerStyle(26);
  rawMCP->Draw("pe");
  auto UnfoldMC = (TH1F *)hMCcorrected->Clone();
  CONSISTlegend->AddEntry("", histName, "");
  CONSISTlegend->AddEntry(UnfoldMC, "Unfolded MC");
  hset(*UnfoldMC, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hoptset(*UnfoldMC, 1., kRed, 0, 200, 1e-9, 1e0);
  UnfoldMC->GetYaxis()->SetNdivisions(505);
  UnfoldMC->Draw("pesame");
  ratioconsistp->cd();
  auto ratconsist = (TH1F *)rawMCP->Clone();
  ratconsist->Divide(ratconsist, UnfoldMC, 1., 1., "B");
  hset(*ratconsist, JetPtTitleX, "Data", 1.2, 1.0, 0.07, 0.07, 0.01, 0.01, 0.07,
       0.07, 510, 505);
  ratconsist->SetMarkerColor();
  ratconsist->SetLineColor(kRed);
  ratconsist->SetMarkerSize(.7);
  ratconsist->SetMarkerStyle(22);
  ratconsist->GetXaxis()->SetRangeUser(0., 200.);
  ratconsist->GetYaxis()->SetRangeUser(0.8, 1.2);
  ratconsist->Draw("pe");

  unfoldp->cd();
  auto UnfoldData = (TH1F *)hDcorrected->Clone();
  unfoldlegend->AddEntry("", histName, "");
  unfoldlegend->AddEntry(UnfoldData, "Unfolded Data");
  hset(*UnfoldData, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hoptset(*UnfoldData, 1., kRed, 0, 200, 1e-9, 1e0);
  UnfoldData->GetYaxis()->SetNdivisions(505);
  UnfoldData->Draw("pe");
  rawMCP->Draw("pesame");
  ratunfoldp->cd();
  auto ratunfold = (TH1F *)rawMCP->Clone();
  ratunfold->Divide(ratunfold, UnfoldData, 1., 1., "B");
  unfoldlegend->AddEntry(rawMCP, "MC Generated");
  hset(*ratunfold, JetPtTitleX, "MC / Data", 1.2, 1.0, 0.07, 0.07, 0.01, 0.01,
       0.07, 0.07, 510, 505);
  ratunfold->SetMarkerColor(kBlack);
  ratunfold->SetLineColor(kBlack);
  ratunfold->SetMarkerSize(.7);
  ratunfold->SetMarkerStyle(22);
  ratunfold->GetXaxis()->SetRangeUser(0., 200.);
  ratunfold->GetYaxis()->SetRangeUser(0., 1.6);
  ratunfold->Draw("pe");

  return jre;
}

void DrawJetMatchings(const std::vector<TString> &fileNames,
                      const std::vector<TString> &histNames,
                      const std::vector<int> &ColorPallete) {
  TFile *TSavefile = new TFile("TimeFrameEff.root", "RECREATE");
  for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
    Filipad2 *JREPad = new Filipad2(++nn, 2, 0.5, 100, 50, 0.7, 1, 1);
    JREPad->Draw();
    TPad *jrepad = JREPad->GetPad(1);
    optFili(*jrepad, 1, 1, 0, 1);
    TPad *ratiojrepads = JREPad->GetPad(2);
    optFili(*ratiojrepads, 1, 1, 0, 0);
    TLegend *legjre =
        new TLegend(0.526316, 0.697391, 0.949761, 0.94087, NULL, "brNDC");
    legjre->SetTextSize(0.065);
    legjre->SetBorderSize(0);
    TLegend *legjre2 =
        new TLegend(0.488038, 0.481159, 0.535885, 0.701449, NULL, "brNDC");
    legjre2->SetTextSize(0.065);
    legjre2->SetBorderSize(0);
    legjre2->SetTextAlign(12);
    legjre2->AddEntry("", "5 GeV #leq #it{p}_{T, jet}^{rec} #leq 100 GeV", "");

    Filipad2 *JRPPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    JRPPad->Draw();
    TPad *jrppad = JRPPad->GetPad(1);
    optFili(*jrppad, 1, 1, 0, 1);
    TPad *ratiojrppads = JRPPad->GetPad(2);
    optFili(*ratiojrppads, 1, 1, 0, 0);
    TLegend *legjrp =
        new TLegend(0.401914, 0.715942, 0.667464, 0.95942, NULL, "brNDC");
    legjrp->SetTextSize(0.05);
    legjrp->SetBorderSize(0);

    Filipad2 *ConsistencyPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    ConsistencyPad->Draw();
    TPad *consistpad = ConsistencyPad->GetPad(1);
    optFili(*consistpad, 1, 1, 0, 1);
    TPad *ratioconsistpad = ConsistencyPad->GetPad(2);
    optFili(*ratioconsistpad, 1, 1, 0, 0);
    TLegend *legconsist =
        new TLegend(0.401914, 0.715942, 0.667464, 0.95942, NULL, "brNDC");
    legconsist->SetTextSize(0.05);
    legconsist->SetBorderSize(0);

    Filipad2 *UnfoldPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    UnfoldPad->Draw();
    TPad *unfoldpad = UnfoldPad->GetPad(1);
    optFili(*unfoldpad, 1, 1, 0, 1);
    TPad *ratiounfoldpad = UnfoldPad->GetPad(2);
    optFili(*ratiounfoldpad, 1, 1, 0, 0);
    TLegend *legunfold =
        new TLegend(0.619617, 0.414493, 0.885167, 0.657971, NULL, "brNDC");
    legunfold->SetTextSize(0.05);
    legunfold->SetBorderSize(0);
    TLegend *legpt3 =
        new TLegend(0.184211, 0.04, 0.232057, 0.25913, NULL, "brNDC");
    legpt3->SetTextSize(0.062);
    legpt3->SetBorderSize(0);
    legpt3->SetTextAlign(12);
    legpt3->AddEntry("", "|#it{#eta}_{jet}| < 0.5", "");
    legpt3->AddEntry("", "Anti-#it{k}_{T}, #it{R} = 0.4", "");
    TLegend *legpt2 =
        new TLegend(0.287081, 0.698551, 0.368421, 0.947826, NULL, "brNDC");
    legpt2->SetTextSize(0.062);
    legpt2->SetBorderSize(0);
    legpt2->SetTextAlign(12);
    legpt2->AddEntry("", "pp #sqrt{#it{s}} = 13.6 TeV", "");
    legpt2->AddEntry("", "#it{p}_{T, track} > 0.15 GeV/#it{c}", "");
    legpt2->AddEntry("", "|#it{#eta}_{track}| < 0.9", "");

    TString filePath = mainDir + fileNames[i];
    DrawJetMatching(filePath.Data(), refPath.Data(), histNames[i].Data(),
                    Nevents(filePath.Data(), Dir, EventObj), jrepad,
                    ratiojrepads, legjre, consistpad, ratioconsistpad,
                    legconsist, unfoldpad, ratiounfoldpad, legunfold,
                    ColorPallete[i], jrppad, ratiojrppads, legjrp, TSavefile);

    jrepad->cd();
    legjre->Draw();
    legjre2->Draw();
    JREPad->C->Print(Form("plots/JRE_R%.1f_%s_ITS%i.pdf",
                          RBIN, histNames[i].Data(), selITS));

    jrppad->cd();
    legjrp->Draw();
    JRPPad->C->Print(Form("plots/JRP_R%.1f_%s_ITS%i.pdf",
                          RBIN, histNames[i].Data(), selITS));

    consistpad->cd();
    legconsist->Draw();
    ConsistencyPad->C->Print(
        Form("plots/ConsistencyCheck_R%.1f_%s_ITS%i.pdf",
             RBIN, histNames[i].Data(), selITS));

    unfoldpad->cd();
    legpt2->Draw();
    legpt3->Draw();
    legunfold->Draw();
    UnfoldPad->C->Print(
        Form("plots/UnfoldedJetPt_R%.1f_%s_ITS%i.pdf", RBIN,
             histNames[i].Data(), selITS));
  }

  for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
    TCanvas *canJetResolution =
        new TCanvas(Form("JetResolution_%s", histNames[i].Data()),
                    Form("jetresolution%s", histNames[i].Data()), 800, 600);
    canJetResolution->SetLogy(1);
    canJetResolution->SetGridx(1);
    canJetResolution->SetGridy(1);
    canJetResolution->Draw();
    TLegend *legjetresolution =
        new TLegend(0.3, 0.13, 0.616541, 0.40435, NULL, "brNDC");
    legjetresolution->SetTextSize(0.04);
    legjetresolution->SetBorderSize(0);
    TLegend *legjetresolutionStat =
        new TLegend(0.0, 0.0, 0.0, 0.0, NULL, "brNDC");
    legjetresolutionStat->SetTextSize(0.03);
    legjetresolutionStat->SetBorderSize(0);
    TString filePath = mainDir + fileNames[i];
    DrawJetResolution(filePath.Data(), histNames[i].Data(), legjetresolution,
                      legjetresolutionStat, ColorPallete[i]);
    legjetresolution->Draw();
    legjetresolutionStat->Draw();
    canJetResolution->Print(
        Form("plots/JetResolution_R%.1f_%s_ITS%i.pdf", RBIN,
             histNames[i].Data(), selITS));
  }

  TSavefile->Close();
}