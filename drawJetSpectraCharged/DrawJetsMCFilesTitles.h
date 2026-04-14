///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Draw macro for jet QA         //////////
////////// author: Joonsuk Bae           //////////
////////// E-mail: jbae@cern.ch          //////////
////////// Last Modified: 01 Sep 2025    //////////
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
const auto ConstituentProcess = 0;
const auto McpJetProcess = 0;
const auto JetProcess = 0;
const auto JetMatchingProcess = 1;
const auto ResponseMatrixQAProcess = 0;
const bool REBINON = true;
const bool NORMEVENTS =
    true;         // set normalization factor (true: Nevts, false: Nobs)
#define NORMEVENT // set axis titles
const bool SEL8WINDOW = true; // true: final distribution constraint to sel8,
                              // false: final dist. in INEL
const bool XSECTION = false; // only set final corrected results (true: xsection
                            // normalization, false: Nevt normalization)
const bool DRAWPLOTS = true;
const bool DRAWRM = true;
const double PlotPtMin = 5;
const double PlotPtMax = 140;

const bool SYSTUNFOLD = false;
const bool PrelimComparison = false;

const TString PrelimComparisonFile = "../22pass7PrelimValidation.root";

// Global variable for InvYData
TH1 *gInvYData = nullptr;

///////////////////////////////////////////////
/// pp 13.6 TeV normalization factors /////////
/// https://its.cern.ch/jira/browse/O2-3720 ///
/// https://alice-notes.web.cern.ch/system/files/notes/analysis/665/2018-09-23-INEL_norm_v2.pdf
/// ///
///////////////////////////////////////////////
const Double_t SigmaINEL = 78.6; // mb (PYTHIA8 Monash ND+SD+DD, 13.6 TeV)
// Correct visible TVX cross-sections (van der Meer scan)
const Double_t kSigmaVisTVX2022 = 53.0; // mb
const Double_t kSigmaVisTVX2023 = 53.4; // mb
const int kDatasetYear = 2023; // Current dataset year (LHC23_pass4_thin)
const Double_t kSigmaVis = (kDatasetYear >= 2023) ? kSigmaVisTVX2023 : kSigmaVisTVX2022;

const double effINELtoselMC =
    0.62539; // 603459/964920; estimated from MC - LHC24f3c.
const double effTrigZvtx10 =
    0.95608; // estimated from LHC22o-pass7-full Gaussian fit.
// const double effTVXtoZvtxsel8 = 0.789; // estimated from LHC22o-pass7.
const double effTVXtoZvtxsel8 = 0.79784948; // estimated from LHC22o-pass7-full.
const double effTVXtoselMC =
    0.83158946; // eff(ITSROBorder/TFBorder) estimated from LHC22o-pass7-full.
// const double effTVXtosel8 = 0.82160428; // estimated from LHC22o-pass7.
const double effTVXtosel8 =
    0.82160428; // eff(ITSROBorder) estimated from LHC22o-pass7-full.

const double effInvYselMCtoTVX =
    0.9226; // invariant yield efficiency from selMC to TVX estimated from MC
            // particle level in local - LHC24f3c.

const TString mainDir = "~/cernbox/workspace/O2Physics/jets/AnalysisResults/";
const TString fRoot = "_AnalysisResults.root";

// Data directory (ref)
// TString refFile = "473574" + fRoot; // small old. version?
// TString refFile = "497532" + fRoot; // 2022 pass7 small with same configuration as full 2022 pass7 (default)
// TString refFile = "619522" + fRoot; // 2022 pass7 small (systematic run) DCA 0.1
// TString refFile = "619713" + fRoot; // 2022 pass7 small (systematic run)
// TString refFile = "619714" + fRoot; // 2022 pass7 small (systematic run) DCA 0.5
// TString refFile = "619715" + fRoot; // 2022 pass7 small (systematic run) sel8Full
// TString refFile = "608781" + fRoot; // 2023 pass4 thin small (3 runs succeed out of 11)

// TString refFile = "498133" + fRoot; // full 2022 pass7 (default)
TString refFile = "496213" + fRoot; // full 2023 pass4 thin
// TString refFile = "484736" + fRoot; // full 2024 pass4 MB

// TString refFile = "628390" + fRoot; // 2023 MB MC (LHC23k4i)
TString refPath = mainDir + refFile;  // Commented out: refFile needs to be defined first
// const char *refName = "LHC22o_pass7";
const char *refName = "LHC23_pass4_thin";
// const char *refName = "2023 MB MC (noRCTsel, noAssoc, custom pTsmear)";
const std::vector<TString> RefFileDirectories = {"jet-spectra-charged"};

// MC directory (comp)
std::vector<TString> fileNames = {
    // "497532" + fRoot, // 2022 pass7 small

    // "594032" + fRoot, // 2023 pass4 thin small (8 runs succeed out of 11)

    // // "495943" + fRoot, // 2023 pass4 skimmed, fJetChLowPt>30 GeV
    // // "497151" + fRoot, // 2023 pass4 skimmed, fJetChLowPt>50 GeV
    // // "496213" + fRoot, // 2023 pass4 thin small
    // // "504362" + fRoot, // 2024 pass4 skimmed, fJetChLowPt>30 GeV
    // // "499055" + fRoot, // 2024 pass4 skimmed, fJetChLowPt>50 GeV
    // // "567621" + fRoot, // 2024 pass4 MB

    // // "515969" + fRoot, // LHC24f3c_fix, tuned w/o pT smearing

    // "515446" + fRoot, // LHC24f3c, tuner A
    // "555891" + fRoot, // LHC24f3c, tuner B
    // "594014" + fRoot, // LHC24f3c, tuner C
    // // // // "514593" + fRoot, // LHC24f3b, w/o tuner
    
    // "516969" + fRoot, // Provisional nominal configuration, LHC25a2b, tuner A
    // "605663" + fRoot, // Provisional nominal configuration, LHC25a2b, tuner A new
    // "596838" + fRoot, // LHC25a2b, tuner B, 1.2
    // "596836" + fRoot, // Default for systematic, LHC25a2b, tuner B, 1.5
    // "596837" + fRoot, // LHC25a2b, tuner B, 1.8
    // "593757" + fRoot, // LHC25a2b, tuner C
    // "502419" + fRoot, // LHC25a2b, tuner D (file-based)
    // "602733" + fRoot, // LHC25a2b, d(Q/pT)*reference, DCA corr auto-detect
    // "605664" + fRoot, // LHC25a2b, d(Q/pT)*2022based, DCA corr auto-detect

    // "515669" + fRoot, // LHC23k4h, tuned w/o pT smearing, tuner A
    // "594179" + fRoot, // LHC23k4h, tuned w/ pT smearing * 1.5, tuner C
    "628390" + fRoot, // LHC23k4i, custom tuner, woRCTselection

    // "LHC24f3c_local_DCAonly" + fRoot, // LHC24f3c local, DCAonly tuner
    // "LHC24f3c_local_pTsmearing1p5" + fRoot, // LHC24f3c local, constPtSmear 1.5
    // "LHC24f3c_local_myTuner" + fRoot, // LHC24f3c local, custom graph tuner

    // "LHC26b5_local_pTHat1p5_myTuner" + fRoot, // LHC26b5_local, CustomTuner, pTHat MCD/MCP (1.5/1.5)
    // "LHC26b5_local_pTHatScan" + fRoot, // LHC26b5_local, pTHat MCD/MCP (1.5/1.5)
    // "LHC26b5_local_pTHatScan" + fRoot, // LHC26b5_local, pTHat 2.0/2.0
    // "LHC26b5_local_pTHatScan" + fRoot, // LHC26b5_local, pTHat 4.0/4.0
    // "LHC26b5_local_pTHatScan" + fRoot, // LHC26b5_local, pTHat 8.0/4.0
    "LHC26b5_wopTHatCut" + fRoot, // LHC26b5_local, pTHat 8.0/4.0 without pTHat cut

    // // "532921test" + fRoot // LHC25a2b, tuned w/o pT smearing, exponentN=4,
    // // Nmax=1.5, process dummy off
    // // "570220" + fRoot, // LHC25a2b, tuned w/o pT smearing, pThat exp=4,Nmax=1.5

    // // Systematic Runs (LHC25a2b)
    // "596836" + fRoot, // Nominal configuration 
    // "596835" + fRoot, // Ambiguous track - TimeMargin 250 ns (10 BC) 
    // "596834" + fRoot, // Ambiguous track - TimeMargin 1000 ns (40 BC) 
    // "596833" + fRoot, // Ambiguous track - TimeMargin 125 ns (5 BC)
    // "596832" + fRoot, // Tracking efficiency - 99%
    // "619511" + fRoot, // Secondary contamination - DCAz 0.1 cm
    // "619512" + fRoot, // Secondary contamination - DCAz 0.2 cm
    // "619513" + fRoot, // Secondary contamination - DCAz 0.5 cm
    // "596838" + fRoot, // Track pT reesolution - delta(Q/pT) * 1.2
    // "596837" + fRoot, // Track pT reesolution - delta(Q/pT) * 1.8
    // "619514" + fRoot, // event selection - selMCFull 
    
};
const std::vector<TString> McFileDirectories = {
    // "jet-spectra-charged_id34410", // 2022 pass7 small

    // "jet-spectra-charged", // 2023 pass4 thin small

    // // "jet-spectra-charged", // 2023 pass4 skimmed, fJetChLowPt>30 GeV
    // // "jet-spectra-charged", // 2023 pass4 skimmed, fJetChHighPt>50 GeV
    // // "jet-spectra-charged", // 2023 pass4 thin small
    // // "jet-spectra-charged", // 2024 pass4 skimmed, fJetChLowPt>30 GeV
    // // "jet-spectra-charged", // 2024 pass4 skimmed, fJetChHighPt>50 GeV
    // // "jet-spectra-charged", // 2024 pass4 MB

    // "jet-spectra-charged_id34413", // LHC24f3c, tuner A
    // "jet-spectra-charged", // LHC24f3c, tuner B
    // "jet-spectra-charged", // LHC24f3c, tuner C

    // "jet-spectra-charged_Nmax1p5",          // Provisional nominal configuration, LHC25a2b, tuner A
    // "jet-spectra-charged",          // Provisional nominal configuration, LHC25a2b, tuner A new
    // "jet-spectra-charged",          // LHC25a2b, tuner B, 1.2
    // "jet-spectra-charged",          // LHC25a2b, tuner B, 1.5
    // "jet-spectra-charged",          // LHC25a2b, tuner B, 1.8
    // "jet-spectra-charged",          // LHC25a2b, tuner C
    // "jet-spectra-charged_Nmax1p5",          // LHC25a2b, tuner D (file-based)
    // "jet-spectra-charged",          // LHC25a2b, d(Q/pT)*reference, DCA corr auto-detect
    // "jet-spectra-charged",          // LHC25a2b, d(Q/pT)*2022based, DCA corr auto-detect

    // "jet-spectra-charged_id34413",          // LHC23k4h, tuned w/o pT smearing, tuner A
    // "jet-spectra-charged",          // LHC23k4h, tuned w/ pT smearing * 1.5, tuner C
    "jet-spectra-charged",          // LHC23k4i, custom tuner, woRCTselection

    // "jet-spectra-charged",          // LHC24f3c local, DCAonly tuner
    // "jet-spectra-charged",          // LHC24f3c local, constPtSmear 1.5
    // "jet-spectra-charged",          // LHC24f3c local, custom graph tuner

    // "jet-spectra-charged",          // LHC26b5_local, customTuner, pTHat MCD/MCP (1.5/1.5)
    // "jet-spectra-charged",          // LHC26b5_local, pTHat MCD/MCP (1.5/1.5)
    // "jet-spectra-charged_pTHat2",          // LHC26b5_local, pTHat MCD/MCP (2.0/2.0)
    // "jet-spectra-charged_pTHat4",          // LHC26b5_local, pTHat MCD/MCP (4.0/4.0)
    // "jet-spectra-charged",          // LHC26b5_local, pTHat MCD/MCP (8.0/4.0)
    "jet-spectra-charged",          // LHC26b5_local, pTHat MCD/MCP (8.0/4.0) without pTHat cut

    // // "jet-spectra-charged_Nmax2", // LHC25a2b, tuned w/o pT smearing, Nmax=2
    // // "jet-spectra-charged", // LHC25a2b, tuned w/o pT smearing, Nmax=4

    // // Systematic Runs (LHC25a2b)
    // "jet-spectra-charged", // default configuration
    // "jet-spectra-charged", // ambiguous track - TimeMargin 250 ns (10 BC)
    // "jet-spectra-charged", // ambiguous track - TimeMargin 1000 ns (40 BC)
    // "jet-spectra-charged", // ambiguous track - TimeMargin 125 ns (5 BC)
    // "jet-spectra-charged", // tracking efficiency - 99%
    // "jet-spectra-charged", // secondary contamination - DCAz 0.1 cm
    // "jet-spectra-charged", // secondary contamination - DCAz 0.2 cm
    // "jet-spectra-charged", // secondary contamination - DCAz 0.5 cm
    // "jet-spectra-charged", // track pT reesolution - delta(Q/pT) * 1.2
    // "jet-spectra-charged", // track pT reesolution - delta(Q/pT) * 1.8
    // "jet-spectra-charged", // event selection - selMCFull  

};
// --- SVD unfolding regularization parameter per MC run (from closure tests) ---
// Must match fileNames entries one-to-one
const std::vector<int> kSVDkReg = {
    4,  // PLACEHOLDER - 594179, LHC23k4h
    4,  // PLACEHOLDER - LHC24f3c local DCAonly
    4,  // PLACEHOLDER - LHC24f3c local pTsmearing1p5
    4,  // PLACEHOLDER - LHC24f3c local myTuner
};
const std::vector<TString> McCollCounterFile = {
    // mainDir + "484760" + fRoot, // 2022 pass7 small dummy run

    // mainDir + "485159" + fRoot, // 2023 pass4 thinn small

    // // trackEfficiency
    // // mainDir + "485159" + fRoot, // 2023 pass4 thin small
    // // mainDir + "567621" + fRoot, // 2024 pass4 MB
    // // mainDir + "485159" + fRoot, // LHC24f3c_fix
    
    // mainDir + "484760" + fRoot, // LHC24f3c, tuner A
    // mainDir + "484760" + fRoot, // LHC24f3c, tuner B
    // mainDir + "484760" + fRoot, // LHC24f3c, tuner C

    // mainDir + "484760" + fRoot, // LHC25a2b is unable at the moment
    // mainDir + "484760" + fRoot, // LHC25a2b is unable at the moment
    // mainDir + "484760" + fRoot, // LHC25a2b is unable at the moment
    // mainDir + "484760" + fRoot, // LHC25a2b is unable at the moment
    // mainDir + "484760" + fRoot, // LHC25a2b is unable at the moment
    // mainDir + "484760" + fRoot, // LHC25a2b is unable at the moment
    // mainDir + "484760" + fRoot, // LHC25a2b is unable at the moment
    // mainDir + "484760" + fRoot, // LHC25a2b, d(Q/pT)*reference, DCA corr auto-detect
    // mainDir + "484760" + fRoot, // LHC25a2b, d(Q/pT)*2022based, DCA corr auto-detect

    mainDir + "484760" + fRoot,  // LHC23k4h dummy
    // mainDir + "484760" + fRoot,  // LHC24f3c local DCAonly dummy
    // mainDir + "484760" + fRoot,  // LHC24f3c local pTsmearing1p5 dummy
    // mainDir + "484760" + fRoot,  // LHC24f3c local myTuner dummy

    // mainDir + "484760" + fRoot,  // LHC26b5 pTHat scan dummy
    // mainDir + "484760" + fRoot,  // LHC26b5 pTHat scan dummy
    // mainDir + "484760" + fRoot,  // LHC26b5 pTHat scan dummy
    // mainDir + "484760" + fRoot,  // LHC26b5 pTHat scan dummy
    // mainDir + "484760" + fRoot,  // LHC26b5 pTHat scan dummy
    mainDir + "484760" + fRoot,  // LHC26b5 pTHat scan dummy without pTHat cut

    // // Systematic Runs (LHC25a2b)
    // mainDir + "484760" + fRoot, // Nominal configuration
    // mainDir + "484760" + fRoot, // Ambiguous track - TimeMargin 250 ns (10 BC)
    // mainDir + "484760" + fRoot, // Ambiguous track - TimeMargin 1000 ns (40 BC)
    // mainDir + "484760" + fRoot, // Ambiguous track - TimeMargin 125 ns (5 BC)
    // mainDir + "484760" + fRoot, // Tracking efficiency - 99%
    // mainDir + "484760" + fRoot, // Secondary contamination - DCAz 0.1 cm
    // mainDir + "484760" + fRoot, // Secondary contamination - DCAz 0.2 cm
    // mainDir + "484760" + fRoot, // Secondary contamination - DCAz 0.5 cm
    // mainDir + "484760" + fRoot, // Track pT reesolution - delta(Q/pT) * 1.2
    // mainDir + "484760" + fRoot, // Track pT reesolution - delta(Q/pT) * 1.2
    // mainDir + "484760" + fRoot, // event selection - selMCFull 
};
// const TString DataLumiCounterFile = mainDir + "500187" + fRoot; // small -
// luminosity calculator const TString DataCollCounterFile = mainDir + "486304"
// + fRoot; // small - trackEfficiency
const TString McLumiCounterFile =
    mainDir + "593755" + fRoot; // LHC25a2b to full 2022o pass7 - luminosity calculator
    // mainDir + "594176" + fRoot; // LHC23k4h to full 2023 pass4 thin - luminosity calculator
const TString DataCollCounterFile =
    mainDir + "493670" + fRoot; // full 2022 pass7 - trackEfficiency
    // mainDir + "493670" + fRoot; // dummy run for full 2023 pass4 thin - trackEfficiency
    // mainDir + "484736" + fRoot; // full 2024 pass4 MB - trackEfficiency

// const TString DataDatasetName = "22 data (LHC22o-pass7 small)";
// const TString DataDatasetName = "22 data (LHC22o-pass7)";
const TString DataDatasetName = "2023 data (LHC23_pass4_Thin_small)";
// const TString DataDatasetName = "2023 MB MC (noRCTsel, noAssoc, custom pTsmear)";
// const TString DataDatasetName = "2024 data (pass1)";
// const TString PlotSaveName = "jetSectraCharged_LHC22o-pass7_LHC25a2b_tunerA_myTuner";
// const TString PlotSaveName = "jetSectraCharged_LHC22o-pass7_LHC24f3c_LHC25a2b";
// const TString PlotSaveName = "jetSectraCharged_LHC23_pass4_Thin_LHC23k4h";
const TString PlotSaveName = "jetSectraCharged_LHC23_pass4_Thin_LHC26b5_local";
// const TString PlotSaveName = "jetSectraCharged_LHC22o-pass7_LHC25a2b_SystematicRuns";
const std::vector<TString> histNames = {
    // "2022 pass7 small",
    
    // "2023 pass4 thin small 9 out of 11",

    // // "2023 pass4 skimmed, fJetChLowPt>30 GeV",
    // // "2023 pass4 skimmed, fJetChHighPt>55 GeV",
    // // "2023 pass4 thin small",
    // // "2024 pass4 skimmed, fJetChLowPt>30 GeV",
    // // "2024 pass4 skimmed, fJetChHighPt>55 GeV",

    // "MB MC, A tune",
    // "MB MC, B tune",
    // "MB MC, C tune",

    // "JJ MC, A tune",
    // "JJ MC, A tune new",
    // "JJ MC, B tune 1.2",
    // "2022 JJ MC (pTsmearing1p5)",
    // "JJ MC, B tune 1.8",
    // "2022 JJ MC",
    // "JJ MC, D tune (file-based)",
    // "JJ MC, d(qOverPt)-reference, DCA corr auto-detect",
    // "JJ MC, d(qOverPt)-2022based, DCA corr auto-detect",

    // "2023 MB MC (A tune)",
    // "2023 MB MC (RCTsel, Assoc, pTsmear 1.5)",
    "2023 MB MC (noRCTsel, noAssoc, custom pTsmear)",

    // "MB MC local (DCAonly)",
    // "MB MC local (pTsmearing1p5)",
    // "MB MC local (myTuner)",

    // "2023 JJ MC",
    // "2023 JJ MC (pTHat MCD: 1.5, MCP: 1.5)",
    // "2023 JJ MC (pTHat 2.0/2.0)",
    // "2023 JJ MC (pTHat 4.0/4.0)",
    // "2023 JJ MC (pTHat 8.0x4.0)",
    "2023 JJ MC (noRCTsel, noAssoc, custom pTsmear)",

    // // Systematic Runs (LHC25a2b)
    // "defaultConfiguration",
    // "ambiguousTrack-TimeMargin250ns",
    // "ambiguousTrack-TimeMargin1000ns",
    // "ambiguousTrack-TimeMargin125ns",
    // "trackingEfficiency-0.99",
    // "secondaryContamination-DCAz0.1cm",
    // "secondaryContamination-DCAz0.2cm",
    // "secondaryContamination-DCAz0.5cm",
    // "trackPtReesolution-dqOverPtX1.2",
    // "trackPtReesolution-dqOverPtX1.8",
    // "eventSelection-selMCFull",
};
const char *trackselection = "globalTracks";
const TString MakeDirName = "../plots/"
                            "AN_Charged-particle-jet-cross-section-in-pp-"
                            "collisions-at-13.6-TeV/Figures/" +
                            PlotSaveName;

const std::vector<Color_t> ColorPallete = {
    kBlack,       kRed,      kBlue + 1,   kGreen + 2,  kOrange + 7,
    kMagenta + 2, kTeal + 3, kViolet + 2, kYellow + 3, kCyan - 6,
    kAzure + 2,   kPink - 7, kSpring + 5, kGray + 2,   kAzure + 8};

// Double_t Trackptbin[17] = {0.15,  2,  4,  6,  8,  10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100}; 
const Double_t Trackptbin[21] = {0.15,  2,  4,  6,  8, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 170, 200};
// const Double_t Trackptbin[21] = {-0.5, 1.5,  3.5,  5.5,   7.5,   9.5,   14.5, 19.5, 24.5, 29.5, 39.5,  49.5,  59.5,  69.5, 79.5, 89.5, 99.5, 119.5, 139.5, 169.5, 199.5};
const Double_t ptbin[23] = {3,  4,  5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20,
                            25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
// Double_t ptbin[27] = {0, 1, 2, 3, 4, 5,  6,  7,  8,  9,  10, 12, 14,  16, 18,
// 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200, 300};
const Double_t ptbinGen[26] = {0,  1,  2,  3,  4,  5,   6,   7,  8,
                               9,  10, 12, 14, 16, 18,  20,  25, 30,
                               40, 50, 60, 70, 85, 100, 140, 200};
// Double_t ptbinGen[27] = {0, 1, 2, 3, 4, 5,  6,  7,  8,  9,  10, 12, 14,  16,
// 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200, 300};
const Int_t nTrackptbin = sizeof(Trackptbin) / sizeof(Trackptbin[0]) - 1;
const Int_t nptBins = sizeof(ptbin) / sizeof(ptbin[0]) - 1;
const Int_t nptBinsGen = sizeof(ptbinGen) / sizeof(ptbinGen[0]) - 1;

const Double_t RBIN = 0.4;
const Double_t deltaEta = 2.0 * (0.9 - RBIN); // jet eta acceptance: 2*etaMax
const Int_t selITS = 1;
const Int_t Npthat = 2;

// no UE subtraction (Reference)
const char *TrackPtObj = "h_track_pt";
const char *TrackEtaObj = "h2_track_eta_track_phi";
const char *TrackPhiObj = "h2_track_eta_track_phi";
const char *ConstPtObj = "h2_jet_pt_track_pt";
const char *ConstEtaObj = "h3_jet_r_jet_pt_track_eta"; // unnecessary?
const char *ConstPhiObj = "h3_jet_r_jet_pt_track_phi"; // unnecessary?
const char *JetPtObj = "h3_jet_pt_jet_eta_jet_phi";
const char *JetPtMCDObj = "h_jet_pt"; // unnecessary?
const char *JetEtaObj = "h3_jet_pt_jet_eta_jet_phi";
const char *JetPhiObj = "h3_jet_pt_jet_eta_jet_phi";
const char *JetNtracksObj = "h2_jet_pt_jet_ntracks;1";
const char *JetAreaObj = "h2_jet_pt_jet_area;1";

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

const char *EventObj = "h_collisions";
const char *EventWObj = "h_collisions_weighted";
const char *TrackPtObjOld = "h3_track_pt_track_eta_track_phi";
const char *TrackEtaObjOld = "h3_track_pt_track_eta_track_phi";
const char *TrackPhiObjOld = "h3_track_pt_track_eta_track_phi";
const char *JetPtMCPObj = "h_jet_pt_part";

const char *JetResolutionObj = "h2_jet_pt_mcp_jet_pt_diff_matchedgeo";
const char *LumiTVXObj = "eventselection-run3/luminosity/hLumiTVX";
const char *LumiTVXafBCcutsObj =
    "eventselection-run3/luminosity/hLumiTVXafterBCcuts";
const char *NTVXObj = "eventselection-run3/luminosity/hCounterTVX";

///////////////////////
///// axis titles /////
///////////////////////
#ifdef NORMEVENT
TString NormDenomTrk = "1/#it{N}_{evt}";
TString NormDenomConst = "1/#it{N}_{evt}";
TString NormDenomJet = "1/#it{N}_{evt}";
TString NormDenomJetDataFinal = "1/#it{N}_{evt}";
TString NormDenomJetMCFinal = "1/#it{N}_{evt}";
#else
TString NormDenomTrk = "1/#it{N}_{trk}";
TString NormDenomConst = "1/#it{N}_{const}";
TString NormDenomJet = "1/#it{N}_{jet}";
TString NormDenomJetDataFinal = "1/#it{N}_{jet}";
TString NormDenomJetMCFinal = "1/#it{N}_{jet}";
#endif

TString RatioTitleY = "MC / Data";
TString TrackPtTitleX = "#it{p}_{T, track}^{reco} (GeV/#it{c})";
TString TrackPtTitleY = Form("%s d#it{N}/d#it{p}_{T}", NormDenomTrk.Data());
TString TrackEtaTitleX = "#it{#eta}_{track}";
TString TrackEtaTitleY = Form("%s d#it{N}/d#it{#eta}", NormDenomTrk.Data());
TString TrackPhiTitleX = "#it{#varphi}_{track}";
TString TrackPhiTitleY = Form("%s d#it{N}/d#it{#varphi}", NormDenomTrk.Data());
TString ConstPtTitleX = "#it{p}_{T, con}^{reco} (GeV/#it{c})";
TString ConstPtTitleY = Form("%s d#it{N}/d#it{p}_{T}", NormDenomConst.Data());
TString ConstEtaTitleX = "#it{#eta}_{con}";
TString ConstEtaTitleY = Form("%s d#it{N}/d#it{#eta}", NormDenomConst.Data());
TString ConstPhiTitleX = "#it{#varphi}_{con}";
TString ConstPhiTitleY =
    Form("%s d#it{N}/d#it{#varphi}", NormDenomConst.Data());
TString JetPtTitleX = "#it{p}_{T, jet}^{reco} (GeV/#it{c})";
TString JetPtTitleY = Form("%s d#it{N}/d#it{p}_{T}", NormDenomJet.Data());
TString JetPtGenTitleX = "#it{p}_{T, jet}^{true} (GeV/#it{c})";
TString JetEtaTitleX = "#it{#eta}_{jet}";
TString JetEtaTitleY = Form("%s d#it{N}/d#it{#eta}", NormDenomJet.Data());
TString JetPhiTitleX = "#it{#varphi}_{jet}";
TString JetPhiTitleY = Form("%s d#it{N}/d#it{#varphi}", NormDenomJet.Data());
TString JetNtracksTitleX = "N_{jet tracks}";
TString JetNtracksTitleY =
    Form("%s d#it{N}/d#it{N}_{jet tracks}", NormDenomJet.Data());
TString JRETitleX = "#it{p}_{T, jet}^{true} (GeV/#it{c})";
TString JRETitleY = Form("%s d#it{N}/d#it{p}_{T}", NormDenomJet.Data());
TString JRPTitleX = "#it{p}_{T, jet}^{reco} (GeV/#it{c})";
TString JRPTitleY = Form("%s d#it{N}/d#it{p}_{T}", NormDenomJet.Data());
TString XSectionTitleY =
    "d^{2}#sigma/d#it{p}_{T}d#it{#eta}  [mb (GeV/#it{c})^{-1}]";
TString JetPtDataFinalTitleX = "#it{p}_{T, jet}^{ch} (GeV/#it{c})";
TString JetPtDataFinalTitleY =
    Form("%s d#it{N}/d#it{p}_{T}", NormDenomJetDataFinal.Data());
TString JetPtMCFinalTitleY =
    Form("%s d#it{N}/d#it{p}_{T}", NormDenomJetMCFinal.Data());
TString JetResolutionTitleX = "(#it{p}_{T, jet}^{true} - #it{p}_{T, "
                              "jet}^{reco}) / #it{p}_{T, jet}^{true}";
TString JetResolutionTitleY =
    Form("%s d#it{N}/d#it{p}_{T}", NormDenomJet.Data());

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
Int_t idata = 0;
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
             Double_t MarkerSize = .75, Int_t LineStyle = 1,
             Int_t LineWidth = 1, Int_t MarkerStyle = 20) {
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
