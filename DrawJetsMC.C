#include "DrawJetsMCFilesTitles.h"
#include "DrawJetsMCFunctions.h"
#include "DrawJetsMCfTrackQA.h"
#include "DrawJetsMCfJetQA.h"

///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Draw macro for jet QA         //////////
////////// author: Joonsuk Bae           ////////// 
////////// E-mail: jbae@cern.ch          //////////
////////// Last Modified: 19 July 2024   //////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

// operate fns
void DrawHistos(const std::vector<TString> &fileNames,
                const std::vector<TString> &histNames,
                const std::vector<Color_t> &ColorPallete) {

  gSystem->MakeDirectory(MakeDirName.Data());
  Double_t NeventsDataUnTrig = 1.;
  Double_t NeventsData = 1.;
  Double_t NeventsMCD = 1.;
  Double_t NeventsMCP = 1.;
  if (NORMEVENTS) {
    NeventsDataUnTrig = Nevents(refPath.Data(),DataDirectory[0].Data(),EventObj, 1);
    NeventsData = Nevents(refPath.Data(),DataDirectory[0].Data(),EventObj, 0);
  }
  std::cout << "NeventsDataUnTrig: " << NeventsDataUnTrig << std::endl;
  std::cout << "NeventsData: " << NeventsData << std::endl;


  if (TrackProcess == 1) {
    // Draw track pT
    Filipad2 *TrackPtPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    TrackPtPad->Draw();
    TPad *trackptpad = TrackPtPad->GetPad(1);
    optFili(*trackptpad, 1, 1, 0, 1);
    TPad *ratiotrackptpad = TrackPtPad->GetPad(2);
    optFili(*ratiotrackptpad, 1, 1, 0, 0);
    TLegend *legtrackpt =
        new TLegend(0.330144,0.715942,0.598086,0.95942,NULL,"brNDC");
    legtrackpt->SetTextSize(0.05);
    legtrackpt->SetBorderSize(0);
    trackptpad->cd();
      TH1 *TrackPtRatio = DrawTrackPt(refPath.Data(), DataDatasetName, NeventsData, legtrackpt,
                                    ColorPallete[0], DataDirectory[0].Data());
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      trackptpad->cd();
      TString filePath = mainDir + fileNames[i];
      if (NORMEVENTS) {
        NeventsMCD = Nevents(filePath.Data(),Directory[0].Data(),EventObj);
      }
      std::cout << "NeventsMCD: " << NeventsMCD << std::endl;
      TH1 *currentHist =
          DrawTrackPt(filePath.Data(), histNames[i].Data(), NeventsMCD, legtrackpt,
                      ColorPallete[i+1], Directory[0].Data());
      TString ratioName =
          TString::Format("RatioHist_tpt_%s", fileNames[i].Data());
      ratiotrackptpad->cd();
      TH1 *ratioHist =
            DrawRatio(ratioName.Data(), TrackPtRatio, currentHist, TrackPtTitleX,
                      RatioTitleY, ColorPallete[i+1], 0, 1);
          // DrawRatioTH1(currentHist, TrackPtRatio);
      // hoptset(*ratioHist, 0, ColorPallete[i+1], 0, PlotPtMax, 0, 1, 1);
      ratioHist->Draw("esame");
    }
    trackptpad->cd();
    legtrackpt->Draw();
    if (DRAWPLOTS) {TrackPtPad->C->Print(Form("%s/TrackPt.pdf", MakeDirName.Data()));}

    // Draw track eta
    Filipad2 *TrackEtaPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    TrackEtaPad->Draw();
    TPad *tracketapad = TrackEtaPad->GetPad(1);
    optFili(*tracketapad, 1, 1, 0, 0);
    TPad *ratiotracketapad = TrackEtaPad->GetPad(2);
    optFili(*ratiotracketapad, 1, 1, 0, 0);
    TLegend *legtracketa =
        new TLegend(0.330144,0.715942,0.598086,0.95942,NULL,"brNDC");
    legtracketa->SetTextSize(0.05);
    legtracketa->SetBorderSize(0);
    tracketapad->cd();
    TH1 *TrackEtaRatio = DrawTrackEta(refPath.Data(), DataDatasetName, NeventsData, legtracketa,
                                      ColorPallete[0], DataDirectory[0].Data());
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      tracketapad->cd();
      TString filePath = mainDir + fileNames[i];
      if (NORMEVENTS) {
        NeventsMCD = Nevents(filePath.Data(),Directory[0].Data(),EventObj);
      }
      TH1 *currentHist =
          DrawTrackEta(filePath.Data(), histNames[i].Data(), NeventsMCD, legtracketa,
                       ColorPallete[i+1], Directory[0].Data());
      TString ratioName =
          TString::Format("RatioHist_teta_%s", fileNames[i].Data());
      ratiotracketapad->cd();
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), TrackEtaRatio, currentHist,
                    TrackEtaTitleX, RatioTitleY, ColorPallete[i+1], 0, 1);
    }
    tracketapad->cd();
    legtracketa->Draw();
    if (DRAWPLOTS) {TrackEtaPad->C->Print(Form("%s/TrackEta.pdf", MakeDirName.Data()));}

    // Draw track phi
    Filipad2 *TrackPhiPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    TrackPhiPad->Draw();
    TPad *trackphipad = TrackPhiPad->GetPad(1);
    optFili(*trackphipad, 1, 1, 0, 0);
    TPad *ratiotrackphipad = TrackPhiPad->GetPad(2);
    optFili(*ratiotrackphipad, 1, 1, 0, 0);
    TLegend *legtrackphi =
        new TLegend(0.330144,0.715942,0.598086,0.95942,NULL,"brNDC");
    legtrackphi->SetTextSize(0.05);
    legtrackphi->SetBorderSize(0);
    trackphipad->cd();
    // TH1* TrackPhiRatio = DrawTrackPhi(refPath.Data(), DataDatasetName,
    // Nevents(refPath.Data(),Dir,EventObj), legtrackphi, kRed);
    TH1 *TrackPhiRatio = DrawTrackPhi(refPath.Data(), DataDatasetName, NeventsData, legtrackphi,
                                      ColorPallete[0], DataDirectory[0].Data());
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      trackphipad->cd();
      TString filePath = mainDir + fileNames[i];
      if (NORMEVENTS) {
        NeventsMCD = Nevents(filePath.Data(),Directory[0].Data(),EventObj);
      }
      TH1 *currentHist =
          DrawTrackPhi(filePath.Data(), histNames[i].Data(), NeventsMCD, legtrackphi,
                       ColorPallete[i+1], Directory[0].Data());
      TString ratioName =
          TString::Format("RatioHist_tphi_%s", fileNames[i].Data());
      ratiotrackphipad->cd();
      float HistMean = currentHist->GetMean(1);
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), TrackPhiRatio, currentHist,
                    TrackPhiTitleX, RatioTitleY, ColorPallete[i+1], 0., 1.);
    }
    trackphipad->cd();
    legtrackphi->Draw();
    if (DRAWPLOTS) {TrackPhiPad->C->Print(Form("%s/TrackPhi.pdf", MakeDirName.Data()));}
  }

  if (ConstituentProcess == 1) {
    // Draw constituents pT
    Filipad2 *ConstPtPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    ConstPtPad->Draw();
    TPad *constptpad = ConstPtPad->GetPad(1);
    optFili(*constptpad, 1, 1, 0, 1);
    TPad *ratioconstptpad = ConstPtPad->GetPad(2);
    optFili(*ratioconstptpad, 1, 1, 0, 0);
    TLegend *legconstpt =
        new TLegend(0.330144,0.715942,0.598086,0.95942,NULL,"brNDC");
    legconstpt->SetTextSize(0.05);
    legconstpt->SetBorderSize(0);
    constptpad->cd();
    TH1 *ConstPtRatio = DrawConstituentPt(
        refPath.Data(), DataDatasetName, NeventsData, legconstpt, ColorPallete[0], DataDirectory[0].Data());
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      constptpad->cd();
      TString filePath = mainDir + fileNames[i];
      if (NORMEVENTS) {
        NeventsMCD = Nevents(filePath.Data(),Directory[0].Data(),EventObj);
      }
      TH1 *currentHist =
          DrawConstituentPt(filePath.Data(), histNames[i].Data(), NeventsMCD, legconstpt,
                            ColorPallete[i+1], Directory[0].Data());
      TString ratioName =
          TString::Format("RatioHist_cpt_%s", fileNames[i].Data());
      ratioconstptpad->cd();
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), ConstPtRatio, currentHist, ConstPtTitleX,
                    RatioTitleY, ColorPallete[i+1], 0., 5);
    }
    constptpad->cd();
    legconstpt->Draw();
    if (DRAWPLOTS) {ConstPtPad->C->Print(Form("%s/ConstPt_R%.1f.pdf",
                                    MakeDirName.Data(), RBIN));}

    // Draw Constituents Eta
    Filipad2 *ConstEtaPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    ConstEtaPad->Draw();
    TPad *constetapad = ConstEtaPad->GetPad(1);
    optFili(*constetapad, 1, 1, 0, 0);
    TPad *ratioconstetapad = ConstEtaPad->GetPad(2);
    optFili(*ratioconstetapad, 1, 1, 0, 0);
    TLegend *legconsteta =
        new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
    // TLegend *legconstetaPtRange =
    //     new TLegend(0.327751, 0.0521739, 0.598086, 0.121739, NULL, "brNDC");
    // legconstetaPtRange->SetTextSize(0.05);
    // legconstetaPtRange->SetBorderSize(0);
    // legconstetaPtRange->AddEntry("", "40 < #it{p}_{T, jet} < 80 GeV", "");
    legconsteta->SetTextSize(0.05);
    legconsteta->SetBorderSize(0);
    constetapad->cd();
    TH1 *ConstEtaRatio = DrawConstituentEta(
        refPath.Data(), DataDatasetName, NeventsData, legconsteta, ColorPallete[0], DataDirectory[0].Data());
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      constetapad->cd();
      TString filePath = mainDir + fileNames[i];
      if (NORMEVENTS) {
        NeventsMCD = Nevents(filePath.Data(),Directory[0].Data(),EventObj);
      }
      TH1 *currentHist =
          DrawConstituentEta(filePath.Data(), histNames[i].Data(), NeventsMCD,
                             legconsteta, ColorPallete[i+1], Directory[0].Data());
      TString ratioName =
          TString::Format("RatioHist_ceta_%s", fileNames[i].Data());
      ratioconstetapad->cd();
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), ConstEtaRatio, currentHist,
                    ConstEtaTitleX, RatioTitleY, ColorPallete[i+1], 0.5, 2.5);
    }
    constetapad->cd();
    legconsteta->Draw();
    // legconstetaPtRange->Draw();
    if (DRAWPLOTS) {ConstEtaPad->C->Print(Form("%s/ConstEta_R%.1f.pdf",
                                    MakeDirName.Data(), RBIN));}

    // Draw Constituents Phi
    Filipad2 *ConstPhiPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    ConstPhiPad->Draw();
    TPad *constphipad = ConstPhiPad->GetPad(1);
    optFili(*constphipad, 1, 1, 0, 0);
    TPad *ratioconstphipad = ConstPhiPad->GetPad(2);
    optFili(*ratioconstphipad, 1, 1, 0, 0);
    TLegend *legconstphi =
        new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
    // TLegend *legconstphiPtRange =
    //     new TLegend(0.327751, 0.0521739, 0.598086, 0.121739, NULL, "brNDC");
    // legconstphiPtRange->SetTextSize(0.05);
    // legconstphiPtRange->SetBorderSize(0);
    // legconstphiPtRange->AddEntry("", "40 < #it{p}_{T, jet} < 80 GeV", "");
    legconstphi->SetTextSize(0.05);
    legconstphi->SetBorderSize(0);
    constphipad->cd();
    TH1 *ConstPhiRatio = DrawConstituentPhi(
        refPath.Data(), DataDatasetName, NeventsData, legconstphi, ColorPallete[0], DataDirectory[0].Data());
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      constphipad->cd();
      TString filePath = mainDir + fileNames[i];
      if (NORMEVENTS) {
        NeventsMCD = Nevents(filePath.Data(),Directory[0].Data(),EventObj);
      }
      TH1 *currentHist =
          DrawConstituentPhi(filePath.Data(), histNames[i].Data(), NeventsMCD,
                             legconstphi, ColorPallete[i+1], Directory[0].Data());
      TString ratioName =
          TString::Format("RatioHist_cphi_%s", fileNames[i].Data());
      ratioconstphipad->cd();
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), ConstPhiRatio, currentHist,
                    ConstPhiTitleX, RatioTitleY, ColorPallete[i+1], 0.5, 2.5);
    }
    constphipad->cd();
    legconstphi->Draw();
    // legconstphiPtRange->Draw();
    if (DRAWPLOTS) {ConstPhiPad->C->Print(Form("%s/ConstPhi_R%.1f.pdf",
                                    MakeDirName.Data(), RBIN));}
  }

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

    TH1 *JetPtRatio = DrawJetPt(refPath.Data(), DataDatasetName, JetPtObj, NeventsData, legjetpt, ColorPallete[0], 0, DataDirectory[0].Data());

    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      jetptpad->cd();
      TString filePath = mainDir + fileNames[i];
      if (NORMEVENTS) {
        NeventsMCD = Nevents(filePath.Data(),Directory[0].Data(),EventObj);
        std::cout << "NeventsMCD_inJetProcess: " << NeventsMCD << std::endl;
      }
      TH1 *currentHist =
          DrawJetPt(filePath.Data(), histNames[i].Data(), JetPtWUEObj, NeventsMCD, legjetpt, ColorPallete[i+1], 0, Directory[0].Data());
      TString ratioName =
          TString::Format("RatioHist_jpt_%s", fileNames[i].Data());
      ratiojetptpad->cd();
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), JetPtRatio, currentHist, JetPtTitleX,
                    RatioTitleY, ColorPallete[i+1], 0., 5);
    }
    jetptpad->cd();
    legjetpt->Draw();
    if (DRAWPLOTS) {JetPtPad->C->Print(Form("%s/JetPt_R%.1f_.pdf",
                                    MakeDirName.Data(), RBIN));}

    // // Draw jet pT MCP (Monte Carlo Particle level)
    // // Create a Filipad2 object for plotting jet pT MCP with ratio panel
    // Filipad2 *JetPtMCPPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    // JetPtMCPPad->Draw();
    // TPad *jetptmcppad = JetPtMCPPad->GetPad(1); // Get upper pad for main plot
    // optFili(*jetptmcppad, 1, 1, 0, 1);
    // TPad *ratiojetptmcppad = JetPtMCPPad->GetPad(2); // Get lower pad for ratio plot
    // optFili(*ratiojetptmcppad, 1, 1, 0, 0);

    // // Create and configure legend
    // TLegend *legjetptmcp =
    //     new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
    // legjetptmcp->SetTextSize(0.05);
    // legjetptmcp->SetBorderSize(0);

    // // Draw reference histogram
    // jetptmcppad->cd();
    // TString JetPtMCPfilePath = mainDir + fileNames[0];
    // if (NORMEVENTS) {
    //   NeventsMCP = Nevents(JetPtMCPfilePath.Data(),Directory[0].Data(),EventObj, 1);
    // } 
    // std::cout << "NeventsMCP: " << NeventsMCP << std::endl;
    // TH1 *JetPtMCPRatio =
    //     DrawJetPtMCP(JetPtMCPfilePath.Data(), histNames[0].Data(), TrackPtObj, NeventsData,
    //                  legjetptmcp, kRed, 1, DataDirectory[0].Data());

    // // Loop over MC files and draw comparison histograms
    // for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
    //   jetptmcppad->cd();
    //   TString filePath = mainDir + fileNames[i];
    //   if (NORMEVENTS) {
    //     NeventsMCP = Nevents(filePath.Data(),Directory[0].Data(),EventObj, 1);
    //   }
    //   TH1 *currentHist =
    //       DrawJetPtMCP(filePath.Data(), histNames[i].Data(), JetPtMCPObj, NeventsMCP, legjetptmcp,
    //                    ColorPallete[i], i, Directory[0].Data());
      
    //   // Draw ratio plots starting from second file
    //   TString ratioName =
    //       TString::Format("RatioHist_jptmcp_%s", fileNames[i].Data());
    //   ratiojetptmcppad->cd();
    //   if (i >= 1) {
    //     TH1 *ratioHist = DrawRatio(ratioName.Data(), JetPtMCPRatio, currentHist,
    //                                JetPtGenTitleX, "MC / Anchored (19 Dec)",
    //                                ColorPallete[i], 0., 2., 1 - 0.3 * i);
    //   }
    // }
    // jetptmcppad->cd();
    // legjetptmcp->Draw();
    // if (DRAWPLOTS) {JetPtMCPPad->C->Print(Form("%s/JetPtMCP_R%.1f.pdf",
    //                                MakeDirName.Data(), RBIN));}

    // Draw jet eta
    Filipad2 *JetEtaPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    JetEtaPad->Draw();
    TPad *jetetapad = JetEtaPad->GetPad(1);
    optFili(*jetetapad, 1, 1, 0, 0);
    TPad *ratiojetetapad = JetEtaPad->GetPad(2);
    optFili(*ratiojetetapad, 1, 1, 0, 0);
    TLegend *legjeteta =
        new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
    // TLegend *legjetetaPtRange =
    //     new TLegend(0.327751, 0.0521739, 0.598086, 0.121739, NULL, "brNDC");
    // legjetetaPtRange->SetTextSize(0.05);
    // legjetetaPtRange->SetBorderSize(0);
    // legjetetaPtRange->AddEntry("", "40 < #it{p}_{T, jet} < 80 GeV", "");
    legjeteta->SetTextSize(0.05);
    legjeteta->SetBorderSize(0);
    jetetapad->cd();
    TH1 *JetEtaRatio = DrawJetEta(refPath.Data(), DataDatasetName, JetEtaObj, NeventsData, legjeteta,
                                  ColorPallete[0], DataDirectory[0].Data());
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      jetetapad->cd();
      TString filePath = mainDir + fileNames[i];
      if (NORMEVENTS) {
        NeventsMCD = Nevents(filePath.Data(),Directory[0].Data(),EventObj);
      }
      TH1 *currentHist =
          DrawJetEta(filePath.Data(), histNames[i].Data(), JetEtaWUEObj, NeventsMCD, legjeteta,
                     ColorPallete[i+1], Directory[0].Data());
      TString ratioName =
          TString::Format("RatioHist_jeta_%s", fileNames[i].Data());
      ratiojetetapad->cd();
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), JetEtaRatio, currentHist, JetEtaTitleX,
                    RatioTitleY, ColorPallete[i+1], 0., 2.2);
    }
    jetetapad->cd();
    legjeteta->Draw();
    // legjetetaPtRange->Draw();
    if (DRAWPLOTS) {JetEtaPad->C->Print(Form("%s/JetEta_R%.1f.pdf",
                                    MakeDirName.Data(), RBIN));}

    // Draw jet phi
    Filipad2 *JetPhiPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    JetPhiPad->Draw();
    TPad *jetphipad = JetPhiPad->GetPad(1);
    optFili(*jetphipad, 1, 1, 0, 0);
    TPad *ratiojetphipad = JetPhiPad->GetPad(2);
    optFili(*ratiojetphipad, 1, 1, 0, 0);
    TLegend *legjetphi =
        new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
    // TLegend *legjetphiPtRange =
    //     new TLegend(0.327751, 0.0521739, 0.598086, 0.121739, NULL, "brNDC");
    // legjetphiPtRange->SetTextSize(0.05);
    // legjetphiPtRange->SetBorderSize(0);
    // legjetphiPtRange->AddEntry("", "40 < #it{p}_{T, jet} < 80 GeV", "");
    legjetphi->SetTextSize(0.05);
    legjetphi->SetBorderSize(0);
    jetphipad->cd();
    TH1 *JetPhiRatio = DrawJetPhi(refPath.Data(), DataDatasetName, JetPhiObj, NeventsData, legjetphi,
                                  ColorPallete[0], DataDirectory[0].Data());
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      jetphipad->cd();
      TString filePath = mainDir + fileNames[i];
      if (NORMEVENTS) {
        NeventsMCD = Nevents(filePath.Data(),Directory[0].Data(),EventObj);
      }
      TH1 *currentHist =
          DrawJetPhi(filePath.Data(), histNames[i].Data(), JetPhiWUEObj, NeventsMCD, legjetphi,
                     ColorPallete[i+1], Directory[0].Data());
      TString ratioName =
          TString::Format("RatioHist_jphi_%s", fileNames[i].Data());
      ratiojetphipad->cd();
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), JetPhiRatio, currentHist, JetPhiTitleX,
                    RatioTitleY, ColorPallete[i+1], 0.5, 2.5);
    }
    jetphipad->cd();
    legjetphi->Draw();
    // legjetphiPtRange->Draw();
    if (DRAWPLOTS) {JetPhiPad->C->Print(Form("%s/JetPhi_R%.1f.pdf",
                                    MakeDirName.Data(), RBIN));}

    // Draw jet ntracks
    Filipad2 *JetNtracksPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    JetNtracksPad->Draw();
    TPad *jetntrackspad = JetNtracksPad->GetPad(1);
    optFili(*jetntrackspad, 1, 1, 0, 1);
    TPad *ratiojetntrackspad = JetNtracksPad->GetPad(2);
    optFili(*ratiojetntrackspad, 1, 1, 0, 0);
    TLegend *legjetntracks =
        new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
    legjetntracks->SetTextSize(0.05);
    legjetntracks->SetBorderSize(0);
    jetntrackspad->cd();
    TH1 *JetNtracksRatio =
        DrawJetNtracks(refPath.Data(), DataDatasetName, JetNtracksObj, NeventsData, legjetntracks, ColorPallete[0],
                       DataDirectory[0].Data());
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      jetntrackspad->cd();
      TString filePath = mainDir + fileNames[i];
      if (NORMEVENTS) {
        NeventsMCD = Nevents(filePath.Data(),Directory[0].Data(),EventObj);
      }
      TH1 *currentHist =
          DrawJetNtracks(filePath.Data(), histNames[i].Data(), JetNtracksWUEObj, NeventsMCD,
                         legjetntracks, ColorPallete[i+1], Directory[0].Data());
      TString ratioName =
          TString::Format("RatioHist_jntracks_%s", fileNames[i].Data());
      ratiojetntrackspad->cd();
      TH1 *ratioHist =
          DrawRatio(ratioName.Data(), JetNtracksRatio, currentHist,
                    JetNtracksTitleX, RatioTitleY, ColorPallete[i+1], 0., 3.1);
    }
    jetntrackspad->cd();
    legjetntracks->Draw();
    if (DRAWPLOTS) {JetNtracksPad->C->Print(Form("%s/JetNtracks_R%.1f.pdf",
                                    MakeDirName.Data(), RBIN));}

    // Draw jet area
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      TCanvas *canJetArea =
          new TCanvas(Form("JetArea_%s", histNames[i].Data()),
                      Form("jetarea%s", histNames[i].Data()), 600, 800);
      canJetArea->SetLogz(1);
      canJetArea->Draw();
      TLegend *legjetarea =
          new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
      legjetarea->SetTextSize(0.04);
      legjetarea->SetBorderSize(0);
      TString filePath = mainDir + fileNames[i];
      if (NORMEVENTS) {
        NeventsMCD = Nevents(filePath.Data(),Directory[0].Data(),EventObj);
      }
      TH2 *currentHist =
          DrawJetArea(filePath.Data(), histNames[i].Data(), JetAreaWUEObj, NeventsMCD, legjetarea,
                      ColorPallete[i+1], Directory[0].Data());
      legjetarea->Draw();
      if (DRAWPLOTS) {canJetArea->Print(Form("%s/JetArea_R%.1f_%s.pdf",
                                    MakeDirName.Data(), RBIN, histNames[i].Data()));}
    }
  }

  if (JetRhoProcess == 1) {
  // Draw <rho_UE> vs leadingjet pT
    TCanvas *canLeadingJetPtRho = new TCanvas("LeadingJetPtRho", "LeadingJetPtRho", 800, 700);
    gStyle->SetOptStat(0);
    canLeadingJetPtRho->Draw();
    setpad(canLeadingJetPtRho, 0.1, 0.15, 0.15);
    TLegend *legleadingjetptrho =
        new TLegend(0.521303,0.155556,0.858396,0.391111,NULL,"brNDC");
    legleadingjetptrho->SetTextSize(0.04);
    legleadingjetptrho->SetBorderSize(0);
    legleadingjetptrho->AddEntry("", DataDatasetName, "");

    // TH1 *LeadingJetPtRhoRatio = DrawLeadingJetPtRho(refPath.Data(), DataDatasetName, LeadingJetPtRhoObj,LeadingJetPtRhoMObj, NeventsData, legleadingjetptrho, ColorPallete[0], 0, DataDirectory[0].Data());

    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      TString filePath = mainDir + fileNames[i];
      if (NORMEVENTS) {
        NeventsMCD = Nevents(filePath.Data(),Directory[0].Data(),EventObj);
        std::cout << "NeventsMCD_inLeadingJetPtRhoProcess: " << NeventsMCD << std::endl;
      }
      DrawLeadingJetPtRho(filePath.Data(), histNames[i].Data(), NeventsMCD, legleadingjetptrho, ColorPallete[i+1], 0, Directory[0].Data());
    }
    legleadingjetptrho->Draw();
    if (DRAWPLOTS) {canLeadingJetPtRho->Print(Form("%s/LeadingJetPtRho_R%.1f.pdf",
                                    MakeDirName.Data(), RBIN));}
  }

  if (JetRandomConeProcess == 1) {
    TCanvas *canDeltaPtRandomCone = new TCanvas("DeltaPtRandomCone", "DeltaPtRandomCone", 800, 700);
    gStyle->SetOptStat(0);
    canDeltaPtRandomCone->Draw();
    optFili(*canDeltaPtRandomCone, 0, 0, 0, 1);
    setpad(canDeltaPtRandomCone, 0.1, 0.15, 0.15);
    TLegend *legdeltaptrandomcone =
        new TLegend(0.537594,0.604444,0.874687,0.884444,NULL,"brNDC");
    legdeltaptrandomcone->SetTextSize(0.04);
    legdeltaptrandomcone->SetBorderSize(0);
    legdeltaptrandomcone->AddEntry("", "LHC24f3 (local, 40 GB)", "");

    // TH1 *LeadingJetPtRhoRatio = DrawLeadingJetPtRho(refPath.Data(), DataDatasetName, LeadingJetPtRhoObj,LeadingJetPtRhoMObj, NeventsData, legleadingjetptrho, ColorPallete[0], 0, DataDirectory[0].Data());

    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      TString filePath = mainDir + fileNames[i];
      if (NORMEVENTS) {
        NeventsMCD = Nevents(filePath.Data(),Directory[0].Data(),EventObj);
        std::cout << "NeventsMCD_inDeltaPtRandomConeProcess: " << NeventsMCD << std::endl;
      }
      DrawDeltaRandomCone(filePath.Data(), histNames[i].Data(), RandomConeObj, RandomConeRandomTrackDirectionObj, RandomConeWoLeadingJetObj, RandomConeRandomTrackDirectionWoOneLeadingJetsObj, RandomConeRandomTrackDirectionWoTwoLeadingJetsObj,
      NeventsMCD, legdeltaptrandomcone, ColorPallete[i+1], 0, Directory[0].Data());
    }
    legdeltaptrandomcone->Draw();
    if (DRAWPLOTS) {canDeltaPtRandomCone->Print(Form("%s/DeltaPtRandomCone_%.1f_.pdf",
                                    MakeDirName.Data(), RBIN));}
  }

  if (JetMatchingProcess == 1) {
    TFile *TSavefile = new TFile("TimeFrameEff.root", "RECREATE");
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      // making JRE, Consistency Check, Unfolded Data results
      Filipad2 *JREPad = new Filipad2(++nn, 2, 0.5, 100, 50, 0.7, 1, 1);
      JREPad->Draw();
      TPad *jrepad = JREPad->GetPad(1);
      optFili(*jrepad, 1, 1, 0, 1);
      TPad *ratiojrepads = JREPad->GetPad(2);
      optFili(*ratiojrepads, 1, 1, 0, 0);
      TLegend *legjre =
          new TLegend(0.538278,0.624348,0.815789,0.874783,NULL,"brNDC");
      legjre->SetTextSize(0.065);
      legjre->SetBorderSize(0);
      TLegend *legjre2 =
          new TLegend(0.488038, 0.481159, 0.535885, 0.701449, NULL, "brNDC");
      legjre2->SetTextSize(0.065);
      legjre2->SetBorderSize(0);
      legjre2->SetTextAlign(12);
      // legjre2->AddEntry("", "5 GeV #leq #it{p}_{T, jet}^{rec} #leq 100 GeV",
                        // "");

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

      Filipad2 *UnfoldPad = new Filipad2(++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
      UnfoldPad->Draw();
      TPad *unfoldpad = UnfoldPad->GetPad(1);
      optFili(*unfoldpad, 0, 0, 0, 1);
      TPad *ratiounfoldpad = UnfoldPad->GetPad(2);
      optFili(*ratiounfoldpad, 0, 0, 0, 0);
      TLegend *legunfold =
          new TLegend(0.543062,0.361491,0.80622,0.575155,NULL,"brNDC");
      legunfold->SetTextSize(0.04);
      legunfold->SetBorderSize(0);
      legunfold->SetFillColorAlpha(0,0); 
      TLegend *legpt3 =
          new TLegend(0.160287,0.0434782,0.258373,0.234783,NULL,"brNDC");
      legpt3->SetTextSize(0.043);
      legpt3->SetBorderSize(0);
      legpt3->SetTextAlign(12);
      legpt3->SetFillColorAlpha(0,0); 
      legpt3->AddEntry("", "|#it{#eta}_{jet}| < 0.5", "");
      legpt3->AddEntry("", "Anti-#it{k}_{T}, #it{R} = 0.4", "");
      legpt3->AddEntry("", "charged-particle jets", "");
      TLegend *legpt2 =
          new TLegend(0.535885,0.664596,0.605263,0.915528,NULL,"brNDC");
      legpt2->SetTextSize(0.043);
      legpt2->SetBorderSize(0);
      legpt2->SetTextAlign(12);
      legpt2->SetFillColorAlpha(0,0); 
      legpt2->AddEntry("", "ALICE Preliminary", "");
      legpt2->AddEntry("", "pp #sqrt{#it{s}} = 13.6 TeV", "");
      legpt2->AddEntry("", "#it{p}_{T, track} > 0.15 GeV/#it{c}", "");
      legpt2->AddEntry("", "|#it{#eta}_{track}| < 0.9", "");

      TString filePath = mainDir + fileNames[i];
      if (NORMEVENTS) {
        NeventsMCD = Nevents(filePath.Data(),Directory[0].Data(),EventObj, 0);
        std::cout << "NeventsMCD_inJetMatchingProcess: " << NeventsMCD << std::endl;
        NeventsMCP = Nevents(filePath.Data(),Directory[0].Data(),EventObj, 1);
      }
      DrawJetMatching(filePath.Data(), refPath.Data(), histNames[i].Data(),
                      NeventsData, NeventsDataUnTrig, NeventsMCD, NeventsMCP,
                      jrepad, ratiojrepads, legjre, consistpad, ratioconsistpad,
                      legconsist, unfoldpad, ratiounfoldpad, legunfold,
                      ColorPallete[i], jrppad, ratiojrppads, legjrp, TSavefile,
                      Directory[0].Data());

      // DrawJetMatching(filePath.Data(), refPath.Data(), histNames[i].Data(),
      // 1, jrepad, ratiojrepads, legjre, consistpad, ratioconsistpad,
      // legconsist, unfoldpad, ratiounfoldpad, legunfold, ColorPallete[i],
      // jrppad, ratiojrppads, legjrp);

      jrepad->cd();
      legjre->Draw();
      legjre2->Draw();
      if (DRAWPLOTS) {JREPad->C->Print(Form("%s/JRE_R%.1f_%s.pdf", MakeDirName.Data(), RBIN, histNames[i].Data()));}
      // JREPad->C->Print(Form("JRE_R%.1f.root", MakeDirName.Data(), RBIN));

      jrppad->cd();
      legjrp->Draw();
      if (DRAWPLOTS) {JRPPad->C->Print(Form("%s/JRP_R%.1f_%s.pdf", MakeDirName.Data(), RBIN, histNames[i].Data()));}

      consistpad->cd();
      legconsist->Draw();
      if (DRAWPLOTS) {ConsistencyPad->C->Print(Form("%s/ConsistencyCheck_R%.1f_%s.pdf",
                                    MakeDirName.Data(), RBIN, histNames[i].Data()));}

      unfoldpad->cd();
      legpt2->Draw();
      legpt3->Draw();
      legunfold->Draw();
      if (DRAWPLOTS) {UnfoldPad->C->Print(Form("%s/CrossSection_Run3_R%.1f_%s.pdf",
                               MakeDirName.Data(), RBIN, histNames[i].Data()));}
    }

    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
      TCanvas *canJetResolution =
          new TCanvas(Form("JetResolution_%s", histNames[i].Data()),
                      Form("jetresolution%s", histNames[i].Data()), 800, 600);
      canJetResolution->SetLogy(1);
      setpad(canJetResolution, 0.02, 0.12, 0.12, 0.01);
      canJetResolution->Draw();
      TLegend *legjetresolution =
          new TLegend(0.50,0.668261,0.92,0.958261, NULL, "brNDC");
      legjetresolution->SetTextSize(0.035);
      legjetresolution->SetBorderSize(0);
      legjetresolution->SetTextAlign(12);
      legjetresolution->SetFillColorAlpha(0,0); 
      // legjetresolution->SetX1(legjetresolution->GetX1() - 0.1);
      legjetresolution->SetEntrySeparation(0.01);  // Decrease separation between entries
      legjetresolution->SetMargin(0.1);
      TString filePath = mainDir + fileNames[i];
      DrawJetResolution(filePath.Data(), histNames[i].Data(), 1., legjetresolution,
                        ColorPallete[i], Directory[0].Data());
      legjetresolution->Draw();
      ALICEfigureLegend("ALICE Simulation", 0.0889724,0.730435,0.338346,0.96, 0.31203,0.184348,0.612782,0.351304);
      if (DRAWPLOTS) {canJetResolution->Print(Form("%s/JetResolution_R%.1f_%s.pdf",
                                   MakeDirName.Data(), RBIN, histNames[i].Data()));}
    }

    TSavefile->Close();
  }
}
