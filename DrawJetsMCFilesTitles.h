///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Draw macro for jet QA         //////////
////////// author: Joonsuk Bae           ////////// 
////////// E-mail: jbae@cern.ch          //////////
////////// Last Modified: 08 Nov 2024    //////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

#include "BSHelper.cxx"
#include "Filipad2.h"
#include <__config>
#include <cstdlib>
#include <iostream>
#include <vector>
using namespace std;
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

///////////////////
/// plot switch ///
///////////////////
const auto TrackProcess = 0;
const auto ConstituentProcess = 1;
const auto JetProcess = 1;
const auto JetMatchingProcess = 1;
const bool REBINON = true;
const bool NORMEVENTS = true; // set normalization factor (true: Nevts, false: Nobs) 
#define NORMEVENT // set axis titles
const bool SEL8WINDOW = true; // true: final distribution constraint to sel8, false: final dist. in INEL
const bool XSECTION = false; // only set final corrected results (true: xsection normalization, false: Nevt normalization)
const bool DRAWPLOTS = true;
const bool DRAWRM = true;
const double PlotPtMin = 10;
const double PlotPtMax = 140;

const bool SYSTUNFOLD = false;

///////////////////////////////////////////////
/// pp 13.6 TeV normalization factors /////////
/// https://its.cern.ch/jira/browse/O2-3720 ///
/// https://alice-notes.web.cern.ch/system/files/notes/analysis/665/2018-09-23-INEL_norm_v2.pdf ///
///////////////////////////////////////////////
const Double_t SigmaINEL = 78.6; // mb
const Double_t SigmaTVX = 59.4; // mb
const Double_t InelNormFactor = SigmaINEL / SigmaTVX;
Double_t P_mu = (0.05 / (1 - TMath::Exp(-0.05))); // mu = average number of TVX-like collisions per BC
double effTrigZvtx10 = 0.971357; // estimated from local MC reco.

const std::vector<TString> DataDirectory = {"jet-finder-charged-qa"};
const std::vector<TString> Directory = {
    "jet-finder-charged-qa" // already CollMatched since 5 Aug 2024
    // "jet-finder-charged-qa_CollMatch"
    };
const TString mainDir = "../../jets/mc/jetfinderQA/AnalysisResults/";
const TString refDir =
    // "../../jets/data/jetfinderQA/AnalysisResults/LHC22o_apass6_MB_small/sel8_globalTracks/";
    // "../../jets/data/jetfinderQA/AnalysisResults/LHC22o_apass6_MB_small/sel8Full_globalTracks/";
    // "../../jets/data/jetfinderQA/AnalysisResults/LHC22o_apass6_MB_full/sel8Full_globalTracks/";
    // "../../jets/data/AnalysisResults/LHC22o_apass6_minBias_small/sel8/";
    // "../../jets/data/jetfinderQA/AnalysisResults/LHC22o_apass7_MB_small/sel8Full_globalTracks/";

    // "../../jets/data/AnalysisResults/LHC22o_apass7_minBias_small/sel8/"; // const < 200 GeV
    // "../../jets/data/AnalysisResults/LHC22o_apass7_minBias_small/sel8/Track100GeV/"; // [HP2024] const < 100 GeV
    // "../../jets/data/AnalysisResults/LHC22o_apass7_minBias_small/sel8/UEsub_TrackTuner_Const100GeV/AreaBasedDoSparse/"; // [QM2025] const < 100 GeV, UE-subtracted
    
    "../../jets/mc/AnalysisResults/LHC24f3b/selMC/trackTuner/Track100GeV/"; // track tuned, track pT < 100 GeV, MCP: selMC w/o zvtxGeV, MCP: selMC w/o zvtx
    // "../../jets/mc/AnalysisResults/LHC24g4/selMC/trackTuner/"; // Area-based method, do Sparse, track tuned, track pT < 100 GeV, MCP: selMC w/o zvtxGeV, MCP: selMC w/o zvtx
// TString refFile = "AnalysisResults.root";
TString refFile = "AnalysisResults.root";
//TString refFile = "CombinedResults.root";
TString refPath = refDir + refFile;
const char *refName = "LHC22o_apass7_MB_small";
std::vector<TString> fileNames = {
    // "../../AnalysisResults/LHC24f3b/selMC/AnalysisResults.root" // anchored to MB
    // "../../AnalysisResults/LHC24f3b/selMC/AnalysisResults.root" // MCP: selMC w/o zvtx
    // "../../AnalysisResults/LHC24f3b/selMC/trackTuner/Track100GeV/AnalysisResults.root" // [HP2024] track tuned, track pT < 100 GeV, MCP: selMC w/o zvtx
    // "../../AnalysisResults/LHC24f3b/selMC/TrackTuner_Const200GeV/AnalysisResults.root" // track tuned, track pT < 200 GeV, MCP: selMC w/o zvtx
    // "../../AnalysisResults/LHC24f3b/selMC/UEsub_TrackTuner_Const200GeV/AreaBased/AnalysisResults.root" // Area-based method, no Sparse, track tuned, track pT < 200 GeV, MCP: selMC w/o zvtx
    // "../../AnalysisResults/LHC24f3b/selMC/UEsub_TrackTuner_Const100GeV/AreaBasedDoSparse/AnalysisResults.root" // Area-based method, do Sparse, track tuned, track pT < 100 GeV, MCP: selMC w/o zvtxGeV, MCP: selMC w/o zvtx

    // Etc. for tests
    // "../../AnalysisResults/LHC24f3b/selMC/trackTuner/Track100GeV/DetPtLarger0GeV/AnalysisResults.root" // det jet pT > 0 GeV, track tuned, track pT < 100 GeV, MCP: selMC w/o zvtx
    // "../../../data/AnalysisResults/LHC22o_apass7_minBias_small/sel8/UEsub_TrackTuner_Const100GeV/AreaBasedDoSparse/CombinedResults.root" // [QM2025] Data, const < 100 GeV, UE-subtracted

    "../../AnalysisResults/LHC24g4/selMC/trackTuner/Track100GeV/AnalysisResults.root", // jet-jet MC
    "../AnalysisResults_MB-gap-2.root",
    "../AnalysisResults_MB-gap-3.root",
    "../AnalysisResults_MB-gap-4.root",
    "../AnalysisResults_MB-gap-5.root",
    "../AnalysisResults_MB-gap-6.root"

    // systematic uncertainties below:
    // "../../AnalysisResults/LHC24f3b/selMC/trackingEfficiency/AnalysisResults.root" // trackingEfficiency 97% uncertainty: track pT < 100 GeV, most recent pass7 anchored to MB, MCP: selMC w/o zvtx
    // "../../AnalysisResults/LHC24f3b/selMC/Track100GeV/AnalysisResults.root" // track pT resolution uncertainty: track pT < 100 GeV, most recent pass7 anchored to MB, MCP: selMC w/o zvtx
    // "../../AnalysisResults/LHC24f3/selMC_syst/trackingEfficiency/AnalysisResults.root"
};
const TString McCollCounterFile = "../../jets/mc/AnalysisResults/LHC24f3b/selMC/trackTuner/Track100GeV/AnalysisResults.root"; // [HP2024] track tuned, track pT < 100 GeV, MCP: selMC w/o zvtx
// const TString McCollCounterFile = "../../jets/mc/AnalysisResults/LHC24g4/selMC/trackEfficiency/AnalysisResults.root";
const TString DataCollCounterFile = "../../jets/data/AnalysisResults/LHC22o_apass7_minBias_small/sel8/AnalysisResults.root";
const TString DataDatasetName = "MB MC (LHC24f3b)";
// const TString MCDatasetName = "MC (LHC24f3b)";
const TString PlotSaveName = "MB-gaps-2to6";
const std::vector<TString> histNames = {
  // "LHC22-pass7-small_sel8_selMC_selectedWindow_woUEsub"
  // "LHC22-pass7-small_Track100GeV_sel8_selMC_selectedWindow_woUEsub"
  // "LHC22-pass7-small_Track100GeV_sel8_selMC_selectedWindow_wUEsub"
  // "HP-Approval_LHC22-pass7_LHC24f3b"
  // "LHC22-pass7-small_sel8_LHC24f3b_Track100GeV_selMC_selectedWindow_woUEsub_trackTuner"
  // "LHC22-pass7-small_sel8_LHC24g4_trackTuner_selMC_selectedWindow_woUEsub" // TrackTuner w/ track pT < 100 GeV
  // "LHC22-pass7-small_sel8_LHC24g4_selMC_selectedWindow_woUEsub"

  // etc.
  "Jet-Jet (LHC24g4)",
  "MB-gap-2",
  "MB-gap-3",
  "MB-gap-4",
  "MB-gap-5",
  "MB-gap-6"
};
const char* trackselection = "globalTracks";
const TString MakeDirName = "plots/AN_Charged-particle-jet-cross-section-in-pp-collisions-at-13.6-TeV/Figures/" + PlotSaveName;

const std::vector<Color_t> ColorPallete = {
    kBlack,         kRed,         kBlue + 1, kGreen + 2, kOrange + 7,
    kMagenta + 2,   kTeal + 3,    kViolet + 2, kYellow + 3, kCyan - 6,
    kAzure + 2,     kPink - 7,    kSpring + 5, kGray + 2,  kAzure + 8
};

// Double_t Trackptbin[17] = {0.15,  2,  4,  6,  8,  10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100};
const Double_t Trackptbin[21] = {0.15,  2,  4,  6,  8,  10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 170, 200};
const Double_t ptbin[21] = {5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
// Double_t ptbin[27] = {0, 1, 2, 3, 4, 5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200, 300};
const Double_t ptbinGen[26] = {0, 1, 2, 3, 4, 5,  6,  7,  8,  9,  10, 12, 14, 16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
// Double_t ptbinGen[27] = {0, 1, 2, 3, 4, 5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200, 300};
const Int_t nTrackptbin = sizeof(Trackptbin) / sizeof(Trackptbin[0]) - 1;
const Int_t nptBins = sizeof(ptbin) / sizeof(ptbin[0]) - 1;
const Int_t nptBinsGen = sizeof(ptbinGen) / sizeof(ptbinGen[0]) - 1;

const Double_t RBIN = 0.4;
const Int_t selITS = 1;
const Int_t Npthat = 2;

// no UE subtraction (Reference)
const char *TrackPtObj = "h3_centrality_track_pt_track_phi";
const char *TrackEtaObj = "h3_centrality_track_pt_track_eta";
const char *TrackPhiObj = "h3_centrality_track_pt_track_phi";
const char *ConstPtObj = "h3_jet_r_jet_pt_track_pt";
const char *ConstEtaObj = "h3_jet_r_jet_pt_track_eta";
const char *ConstPhiObj = "h3_jet_r_jet_pt_track_phi";
const char *JetPtObj = "h_jet_pt";
const char *JetPtMCDObj = "h_jet_pt";
const char *JetEtaObj = "h3_jet_r_jet_pt_jet_eta";
const char *JetPhiObj = "h3_jet_r_jet_pt_jet_phi";
const char *JetNtracksObj = "h_jet_ntracks";
const char *JetAreaObj = "h3_jet_r_jet_pt_jet_area";

// // UE subtraction (Comparison)
// const char *TrackPtWUEObj = "h3_centrality_track_pt_rhoareasubtracted";
// const char *TrackEtaWUEObj = "h2_centrality_track_eta_rhoareasubtracted";
// const char *TrackPhiWUEObj = "h2_centrality_track_phi_rhoareasubtracted";
// const char *ConstPtWUEObj = "h3_jet_r_jet_pt_track_pt_rhoareasubtracted";
// const char *ConstEtaWUEObj = "h3_jet_r_jet_pt_track_eta_rhoareasubtracted";
// const char *ConstPhiWUEObj = "h3_jet_r_jet_pt_track_phi_rhoareasubtracted";
// const char *JetPtWUEObj = "h_jet_pt_rhoareasubtracted";
// const char *JetPtMCDWUEObj = "h_jet_pt_rhoareasubtracted";
// const char *JetEtaWUEObj = "h3_jet_r_jet_pt_jet_eta_rhoareasubtracted";
// const char *JetPhiWUEObj = "h3_jet_r_jet_pt_jet_phi_rhoareasubtracted";
// const char *JetNtracksWUEObj = "h_jet_ntracks_rhoareasubtracted";
// const char *JetAreaWUEObj = "h3_jet_r_jet_pt_jet_area_rhoareasubtracted";

const char *TrackPtWUEObj = TrackPtObj;
const char *TrackEtaWUEObj = TrackEtaObj;
const char *TrackPhiWUEObj = TrackPhiObj;
const char *ConstPtWUEObj = ConstPtObj;
const char *ConstEtaWUEObj = ConstEtaObj;
const char *ConstPhiWUEObj = ConstPhiObj;
const char *JetPtWUEObj = JetPtObj;
const char *JetPtMCDWUEObj = JetPtMCDObj;
const char *JetEtaWUEObj = JetEtaObj;
const char *JetPhiWUEObj = JetPhiObj;
const char *JetNtracksWUEObj = JetNtracksObj;
const char *JetAreaWUEObj = JetAreaObj;

const char *LumiObj = "jet-luminosity-calculator/counter";
const char *EventObj = "h_collisions";
const char *EventWObj = "h_collisions_weighted";
const char *TrackPtObjOld = "h3_track_pt_track_eta_track_phi";
const char *TrackEtaObjOld = "h3_track_pt_track_eta_track_phi";
const char *TrackPhiObjOld = "h3_track_pt_track_eta_track_phi";
const char *JetPtMCPObj = "h_jet_pt_part";

const char *JetResolutionObj =
    "h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedgeo";
const char *JetResolutionObjOld = "h3_jet_r_jet_pt_part_jet_pt_diff";
const char *NTVXObj = "bc-selection-task/hCounterTVX";


///////////////////////
///// axis titles /////
///////////////////////
#ifdef NORMEVENT
  const char* NormDenomTrk = "1/#it{N}_{evt}";
  const char* NormDenomConst = "1/#it{N}_{evt}";
  const char* NormDenomJet = "1/#it{N}_{evt}";
  const char* NormDenomJetDataFinal = "1/#it{N}_{evt}";
  const char* NormDenomJetMCFinal = "1/#it{N}_{evt}";
#else
  const char* NormDenomTrk = "1/#it{N}_{trk}";
  const char* NormDenomConst = "1/#it{N}_{const}";
  const char* NormDenomJet = "1/#it{N}_{jet}";
  const char* NormDenomJetDataFinal = "1/#it{N}_{jet}";
  const char* NormDenomJetMCFinal = "1/#it{N}_{jet}";
#endif

TString RatioTitleY = "MC / Data";
TString TrackPtTitleX = "#it{p}_{T, track}^{reco} (GeV/#it{c})";
TString TrackPtTitleY = "1/#it{N}_{trk} d#it{N}/d#it{p}_{T}";
TString TrackEtaTitleX = "#it{#eta}_{track}";
TString TrackEtaTitleY = Form("%s d#it{N}/d#it{#eta}", NormDenomTrk);
TString TrackPhiTitleX = "#it{#varphi}_{track}";
TString TrackPhiTitleY = Form("%s d#it{N}/d#it{#varphi}", NormDenomTrk);
TString ConstPtTitleX = "#it{p}_{T, con}^{reco} (GeV/#it{c})";
TString ConstPtTitleY = "1/#it{N}_{const} d#it{N}/d#it{p}_{T}";
TString ConstEtaTitleX = "#it{#eta}_{con}";
TString ConstEtaTitleY = Form("%s d#it{N}/d#it{#eta}", NormDenomConst);
TString ConstPhiTitleX = "#it{#varphi}_{con}";
TString ConstPhiTitleY = Form("%s d#it{N}/d#it{#varphi}", NormDenomConst);
TString JetPtTitleX = "#it{p}_{T, jet}^{reco} (GeV/#it{c})";
TString JetPtTitleY = "1/#it{N}_{jet} d#it{N}/d#it{p}_{T}";
TString JetPtGenTitleX = "#it{p}_{T, jet}^{true} (GeV/#it{c})";
TString JetEtaTitleX = "#it{#eta}_{jet}";
TString JetEtaTitleY = Form("%s d#it{N}/d#it{#eta}", NormDenomJet);
TString JetPhiTitleX = "#it{#varphi}_{jet}";
TString JetPhiTitleY = Form("%s d#it{N}/d#it{#varphi}", NormDenomJet);
TString JetNtracksTitleX = "N_{jet tracks}";
TString JetNtracksTitleY = Form("%s d#it{N}/d#it{N}_{jet tracks}", NormDenomJet);
TString JRETitleX = "#it{p}_{T, jet}^{true} (GeV/#it{c})";
TString JRETitleY = Form("%s d#it{N}/d#it{p}_{T}", NormDenomJet);
TString JRPTitleX = "#it{p}_{T, jet}^{reco} (GeV/#it{c})";
TString JRPTitleY = Form("%s d#it{N}/d#it{p}_{T}", NormDenomJet);
TString XSectionTitleY = "d^{2}#sigma/d#it{p}_{T}d#it{#eta}  [mb (GeV/#it{c})^{-1}]";
TString JetPtDataFinalTitleX = "#it{p}_{T, jet}^{ch} (GeV/#it{c})";
TString JetPtDataFinalTitleY = Form("%s d#it{N}/d#it{p}_{T}", NormDenomJetDataFinal);
TString JetPtMCFinalTitleY = Form("%s d#it{N}/d#it{p}_{T}", NormDenomJetMCFinal);
TString JetResolutionTitleX = "(#it{p}_{T, jet}^{true} - #it{p}_{T, jet}^{reco}) / #it{p}_{T, jet}^{true}";
TString JetResolutionTitleY = Form("%s d#it{N}/d#it{p}_{T}", NormDenomJet);

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
Int_t n = 0;
Int_t nn = 0;
Int_t ii = 0;
void setpad(TVirtualPad *pad, Double_t tmargin = 0.02, Double_t bmargin = 0.15,
            Double_t lmargin = 0.13, Double_t rmargin = 0.05) {
  pad->SetTopMargin(tmargin);
  pad->SetLeftMargin(lmargin);
  pad->SetRightMargin(rmargin);
  pad->SetBottomMargin(bmargin);
  pad->SetName(Form("c%d", ++n));
}
template <typename T>
void hset(T &hid, TString xtit = "", TString ytit = "", double titoffx = 0.9,
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
template <typename T>
void hoptset(T &hid, Double_t N = 1, Color_t color = kBlack, Double_t minX = 0,
             Double_t maxX = 100, Double_t minY = 0, Double_t maxY = 1,
             Double_t MarkerSize = .75, Int_t LineStyle = 1, Int_t LineWidth = 1, Int_t MarkerStyle = 20) {
  if (N != 0) {
    if (N == 1) {
      hid.Scale(1. / hid.Integral(), "width");
    } else if (N == 2) {
      hid.Scale(1. / hid.Integral(), "");
    } else {
      hid.Scale(1. / N, "width");
    }
  }

  hid.SetMarkerColor(color);
  hid.SetLineColor(color);
  hid.SetMarkerSize(MarkerSize);
  hid.SetMarkerStyle(MarkerStyle); 
  hid.SetLineStyle(LineStyle); 
  hid.SetLineWidth(LineWidth);

  hid.GetXaxis()->SetRangeUser(minX, maxX);
  hid.GetYaxis()->SetRangeUser(minY, maxY);

  hid.SetFillColorAlpha(color, 0.3);

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


// declare fns
void DrawHistos(const std::vector<TString> &fileNames,
                const std::vector<TString> &histNames,
                const std::vector<Color_t> &ColorPallete);
// TLegend* legconstpt;

