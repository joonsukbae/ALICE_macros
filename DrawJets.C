#include <__config>
#include <iostream>
#include <vector>
#include <dirent.h>
#include <regex>
#include "BSHelper.cxx"
#include "Filipad2.h"
using namespace std;
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

// plots switch
auto TrackProperties = 1;
auto ConstituentProperties = 1;
auto JetProperties = 1;

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
Int_t n=0;
Int_t nn=0;
void setpad(TVirtualPad *pad)
{
    pad->SetTopMargin(0.02);
    pad->SetLeftMargin(0.13);
    pad->SetRightMargin(0.05);
    pad->SetBottomMargin(0.15);
    pad->SetName(Form("c%d", ++n));
}
void hset(TH1& hid, TString xtit="", TString ytit="",
        double titoffx = 0.9, double titoffy = 1.2,
        double titsizex = 0.06, double titsizey = 0.06,
        double labeloffx = 0.01, double labeloffy = 0.001,
        double labelsizex = 0.05, double labelsizey = 0.05,
        int divx = 510, int divy=510)
{
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
void hoptset(TH1& hid, Float_t N = 1, Color_t color = kBlack, Double_t minX = 0, Double_t maxX = 100, Double_t minY = 0, Double_t maxY = 1) {
    // hid.Scale(1./N, "width");
    hid.Scale(1./hid.Integral(), "width");
    hid.SetMarkerColor(color);
    hid.SetLineColor(color);
    hid.SetMarkerSize(.5);
    hid.SetMarkerStyle(22);
    hid.GetXaxis()->SetRangeUser(minX,maxX);
    hid.GetYaxis()->SetRangeUser(minY,maxY);
    hid.SetFillColorAlpha(color, 0.3);
}
void optFili(TPad& pid, Int_t gridx, Int_t gridy, Int_t logx, Int_t logy) {
    pid.SetGridy(gridx);
    pid.SetGridx(gridy);
    pid.SetLogx(logx);
    pid.SetLogy(logy);
}

Double_t Trackptbin[15] = {0, 2, 4, 6, 8, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80};
Double_t ptbin[21] = {5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
Double_t ptbinHEP[20] = {5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140};
Int_t nTrackptbin = sizeof(Trackptbin) / sizeof(Trackptbin[0]) - 1;
Int_t nptBins = sizeof(ptbin) / sizeof(ptbin[0]) - 1;
Int_t nptBinsHEP = sizeof(ptbinHEP) / sizeof(ptbinHEP[0]) - 1;
// Double_t RBIN = 0.2;
Double_t RBIN = 0.4;
// Double_t RBIN = 0.6;
const char *Dir = "jet-finder-charged-qa";
const char *JEventDir = "jet-finder-charged-qa";
// const char *JEventDir = "jet-filter/spectra";
// const char *TEventDir = "jet-finder-charged-qa";
const char *JEventObj = "h_collisions";
// const char *JEventObj = "fProcessedEvents";
const char *TEventObj = "h_collisions";
TString mainDir = "../../jets/data/jetfinderQA/AnalysisResults/LHC22_pass4_highIR/";
// TString mainDir = "../../jets/data/jetfinderQA/AnalysisResults/LHC22_apass4/";
const char* DatasetName = "LHC22_pass4_highIR_vs_LHC22r_apass4_lowIR";
TString refDir = "../../jets/data/jetfinderQA/AnalysisResults/LHC22_apass4/LHC22r_apass4/";
// TString refDir = "../../jets/data/jetfinderQA/AnalysisResults/LHC22_apass4/LHC22r_apass4/";
TString refFile = "AnalysisResults.root";
// TString refFile = "AnalysisResults_LHC23d1k_HY_2.root";
TString refPath = refDir + refFile;
const char* refName = "LHC22r_apass4_lowIR (ref.)";

const char* TrackPtObjOld = "h_track_pt";
const char* TrackPtObj = "h2_centrality_track_pt";
const char* TrackEtaObjOld = "h_track_eta";
const char* TrackEtaObj = "h2_centrality_track_eta";
const char* TrackPhiObjOld = "h_track_phi";
const char* TrackPhiObj = "h2_centrality_track_phi";
const char* ConstPtObj = "h3_jet_r_jet_pt_track_pt";
const char* ConstEtaObj = "h3_jet_r_jet_pt_track_eta";
const char* ConstPhiObj = "h3_jet_r_jet_pt_track_phi";
const char* JetPtObj = "h_jet_pt";
const char* JetEtaObj = "h_jet_eta";
const char* JetPhiObj = "h_jet_phi";
const char* JetNtracksObj = "h_jet_ntracks";

const char* RatioTitleY = "Data / Ref.";
const char* TrackPtTitleX = "#it{p}_{T, track}^{reco} (GeV/c)";
const char* TrackPtTitleY = "1/N dN/d#it{p}_{T}";
const char* TrackEtaTitleX = "#it{#eta}_{track}";
const char* TrackEtaTitleY = "1/N dN/d#it{#eta}";
const char* TrackPhiTitleX = "#it{#varphi}_{track}";
const char* TrackPhiTitleY = "1/N dN/d#it{#varphi}";
const char* ConstPtTitleX = "#it{p}_{T, con}^{reco} (GeV/c)";
const char* ConstPtTitleY = "1/N dN/d#it{p}_{T}";
const char* ConstEtaTitleX = "#it{#eta}_{con}";
const char* ConstEtaTitleY = "1/N dN/d#it{#eta}";
const char* ConstPhiTitleX = "#it{#varphi}_{con}";
const char* ConstPhiTitleY = "1/N dN/d#it{#varphi}";
const char* JetPtTitleX = "#it{p}_{T, jet}^{reco} (GeV/c)";
const char* JetPtTitleY = "1/N dN/d#it{p}_{T}";
const char* JetEtaTitleX = "#it{#eta}_{jet}";
const char* JetEtaTitleY = "1/N dN/d#it{#eta}";
const char* JetPhiTitleX = "#it{#varphi}_{jet}";
const char* JetPhiTitleY = "1/N dN/d#it{#varphi}";
const char* JetNtracksTitleX = "N_{jet tracks}";
const char* JetNtracksTitleY = "1/N dN/dN_{jet tracks}";





// declare fns
void DrawHistos(const std::vector<TString>& fileNames, const std::vector<TString>& histNames, const std::vector<int>& ColorPallete);
// TLegend* legconstpt;

std::vector<TString> getFileList(const TString& directory) {
    std::vector<TString> fileNames;
    DIR* dirp = opendir(directory.Data());
    struct dirent* dp;
    std::regex pattern("AnalysisResults_([0-9]+)\\.root"); // 정규 표현식으로 run number가 있는 파일만 찾음

    while ((dp = readdir(dirp)) != nullptr) {
        TString fileName = dp->d_name;
        if (std::regex_match(fileName.Data(), pattern)) {
            fileNames.push_back(fileName);
        }
    }

    closedir(dirp);
    return fileNames;
}

std::vector<TString> getHistNames(const std::vector<TString>& fileNames) {
    std::vector<TString> histNames;
    for (const auto& fileName : fileNames) {
        std::smatch match;
        std::regex pattern("AnalysisResults_([0-9]+)\\.root");
        std::string fileNameStr = fileName.Data();
        if (std::regex_search(fileNameStr, match, pattern) && match.size() > 1) {
            histNames.push_back(Form("%s_", DatasetName) + match.str(1));
        }
    }
    return histNames;
}

void DrawJets() {
    std::vector<TString> fileNames = {"AnalysisResults.root", "LHC22m/AnalysisResults.root", "LHC22o/AnalysisResults.root", "LHC22p/AnalysisResults.root", "LHC22r/AnalysisResults.root", "LHC22t/AnalysisResults.root"};
    std::vector<TString> histNames = {"LHC22_pass4_highIR", "LHC22m_pass4_highIR","LHC22o_pass4_highIR", "LHC22p_pass4_highIR", "LHC22r_pass4_highIR", "LHC22t_pass4_highIR", "520259 (MC)","520471 (MC)", "520472 (MC)", "520473 (MC)"};
    // std::vector<TString> fileNames = {"AnalysisResults_520294.root", "../../../../../mc/jetfinderQA/AnalysisResults/LHC23d1k/AnalysisResults_520294.root"};
    // std::vector<TString> histNames = {"LHC22f_520294", "LHC23d1k_520294"};
    // auto fileNames = getFileList(mainDir);
    // auto histNames = getHistNames(fileNames);
    // for (const auto& fileName : fileNames) {
    //     std::cout << fileName << std::endl;
    // }
    // for (const auto& histName : histNames) {
    //     std::cout << histName << std::endl;
    // }
    // std::vector<TString> fileNames = {"LHC22m_apass4/AnalysisResults.root", "LHC22o_apass4/AnalysisResults.root", "LHC22q_apass4/AnalysisResults.root", "LHC22f_apass4/AnalysisResults.root"};
    // std::vector<TString> histNames = {"LHC22m", "LHC22o", "LHC22q", "LHC22f", "hist6"};
    std::vector<int> ColorPallete = {kRed, kBlue, kGreen+1, kMagenta+1, kCyan+1, kOrange+1, kYellow+2, kAzure+7,
    kViolet+1, kSpring-6, kPink+7, kTeal+3, kAzure+2, kOrange-3, kSpring+8,
    kMagenta-3, kYellow-3, kRed-4, kGreen-5, kBlue-6};

    DrawHistos(fileNames, histNames, ColorPallete);

    // return 0;
}

// define fns
float Nevents(const char* fileName, const char* eventDir, const char* eventObj) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(eventDir);
    dir->cd();

    auto Nevents = (TH1D *) gROOT->FindObject(eventObj);
	float nevents = Nevents->GetBinContent(Nevents->FindBin(0.5));
    
    return nevents;
}

TH1* DrawRatio(const char* ratioName, TH1* refHist, TH1* testHist, TString AxisTitleX, TString AxisTitleY, Color_t colorID, Double_t Ymin, Double_t Ymax) {
    TH1* ratioHist = (TH1*) refHist->Clone(ratioName);
    hset(*ratioHist, AxisTitleX, AxisTitleY, 1.2, 1.0, 0.07, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 510);
    ratioHist->Divide(testHist, ratioHist, 1., 1., "B");
    ratioHist->GetYaxis()->SetRangeUser(Ymin, Ymax);
    ratioHist->SetMarkerColor(colorID);
    ratioHist->SetLineColor(colorID);
    ratioHist->Draw("pesame");

    return ratioHist;
}
// Draw Tracks
TH1* DrawTrackPt(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();

    TH1 *TrackPt = nullptr;
    TH1 *temp1 = (TH1 *) gROOT->FindObject(TrackPtObjOld);
        std::cout << "temp1: "<< temp1 << std::endl;
    TH2 *temp2 = (TH2 *) gROOT->FindObject(TrackPtObj);
        std::cout << "temp2: "<< temp2 << std::endl;

    if (temp1) {
        TrackPt = temp1;
        std::cout << "Old ver. track exist" << std::endl;
    } else if (temp2) {
        TrackPt = temp2->ProjectionY(Form("trackpt_%s",histName), 0, temp2->GetNbinsX(), "e");
        std::cout << "New ver. track exist" << std::endl;
    }

    if (!TrackPt) {
    std::cout << "Error: No valid TrackPt object found." << std::endl;
    }

    // TH1 *TrackPt = (TH1 *) gROOT->FindObject(TrackPtObj);
  	TrackPt = TrackPt->Rebin(nptBins,Form("TrackPt_%s",histName),ptbin); 
    legend->AddEntry(TrackPt, histName);
    hset(*TrackPt, TrackPtTitleX, TrackPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*TrackPt, Nevts, colorID, 0, 100, 1e-7, 1e0);
    TrackPt->Draw("pesame");

    return TrackPt;
}
TH1* DrawTrackEta(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    
    TH1 *TrackEta = nullptr;
    TH1 *temptrketa1 = (TH1 *) file->Get(Form("%s/%s", Dir, TrackEtaObjOld));;
        std::cout << "temptrketa1: "<< temptrketa1 << std::endl;
    TH2 *temptrketa2 = (TH2 *) file->Get(Form("%s/%s", Dir, TrackEtaObj));
        std::cout << "temptrketa2: "<< temptrketa2 << std::endl;

    if (temptrketa1) {
        TrackEta = temptrketa1;
        std::cout << "Old ver. track eta exist" << std::endl;
    } else if (temptrketa2) {
        TrackEta = temptrketa2->ProjectionY(Form("tracketa_%s",histName), temptrketa2->GetXaxis()->FindBin(1e-6), temptrketa2->GetXaxis()->FindBin(100 - 1e-6), "e");
        std::cout << "New ver. track eta exist" << std::endl;
    }
    // temptrketa1 = nullptr;
    // temptrketa2 = nullptr;

    if (!TrackEta) {
    std::cout << "Error: No valid TrackEta object found." << std::endl;
    }
    // TH1 *TrackEta = (TH1 *) gROOT->FindObject(TrackEtaObj);
  	// TrackEta = TrackEta->Rebin(nptBins,Form("TrackEta_%s",histName),ptbin); 
    legend->AddEntry(TrackEta, histName);
    hset(*TrackEta, TrackEtaTitleX, TrackEtaTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*TrackEta, Nevts, colorID, -0.9, 0.9, .45, .7);
    TrackEta->Draw("pesame");

    std::cout << "End of DrawTrackEta loop" << std::endl;

    return TrackEta;
}
TH1* DrawTrackPhi(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();

    TH1 *TrackPhi = nullptr;
    TH1 *temptrkphi1 = (TH1 *) file->Get(Form("%s/%s", Dir, TrackPhiObjOld));;
        std::cout << "temptrkphi1: "<< temptrkphi1 << std::endl;
    TH2 *temptrkphi2 = (TH2 *) file->Get(Form("%s/%s", Dir, TrackPhiObj));;
        std::cout << "temptrkphi2: "<< temptrkphi2 << std::endl;

    if (temptrkphi1) {
        TrackPhi = temptrkphi1;
        std::cout << "Old ver. trackphi exist" << std::endl;
    } else if (temptrkphi2) {
        TrackPhi = temptrkphi2->ProjectionY(Form("trackphi_%s",histName), temptrkphi2->GetXaxis()->FindBin(1e-6), temptrkphi2->GetXaxis()->FindBin(100 - 1e-6), "e");
        std::cout << "New ver. track phi exist" << std::endl;
  	    TrackPhi->Rebin(2); 
    }
    temptrkphi1 = nullptr;
    temptrkphi2 = nullptr;

    if (!TrackPhi) {
    std::cout << "Error: No valid TrackEta object found." << std::endl;
    }

    // TH1 *TrackPhi = (TH1 *) gROOT->FindObject(TrackPhiObj);
    legend->AddEntry(TrackPhi, histName);
    hset(*TrackPhi, TrackPhiTitleX, TrackPhiTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*TrackPhi, Nevts, colorID, 0, 2*TMath::Pi(), 0.05, .3);
    TrackPhi->Draw("pesame");

    return TrackPhi;
}
// Draw Constituents
TH1* DrawConstituentPt(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH3D *H3ConstituentPt = (TH3D *) gROOT->FindObject(ConstPtObj);
    TH1 *ConstituentPt = H3ConstituentPt->ProjectionZ(Form("ConstituentPt_%s",histName),H3ConstituentPt->GetXaxis()->FindBin(RBIN), H3ConstituentPt->GetXaxis()->FindBin(RBIN), 0, H3ConstituentPt->GetNbinsY());
    ConstituentPt = ConstituentPt -> Rebin(nTrackptbin,"constituentpt",Trackptbin);
    legend->AddEntry(ConstituentPt, histName);
    hset(*ConstituentPt, ConstPtTitleX, ConstPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*ConstituentPt, Nevts, colorID, 0, 80, 1e-8, 1);
    ConstituentPt->Draw("pesame");

    return ConstituentPt;
}
TH1D* DrawConstituentEta(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH3D *H3ConstituentEta = (TH3D *) gROOT->FindObject(ConstEtaObj);
    TH1D *ConstituentEta = H3ConstituentEta->ProjectionZ(Form("ConstituentEta_%s",histName),H3ConstituentEta->GetXaxis()->FindBin(RBIN+1e-6), H3ConstituentEta->GetXaxis()->FindBin(RBIN+1e-5), 0, H3ConstituentEta->GetNbinsY());
    // ConstituentEta = ConstituentEta -> Rebin(nTrackptbin,"constituenteta",Trackptbin);
    legend->AddEntry(ConstituentEta, histName);
    hset(*ConstituentEta, ConstEtaTitleX, ConstEtaTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*ConstituentEta, Nevts, colorID, -0.9, 0.9, 0., 2.);
    ConstituentEta->Draw("pesame");

    return ConstituentEta;
}
TH1* DrawConstituentPhi(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID, Bool_t REBIN = true) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH3D *H3ConstituentPhi = (TH3D *) gROOT->FindObject(ConstPhiObj);
    TH1 *ConstituentPhi = H3ConstituentPhi->ProjectionZ(Form("ConstituentPhi_%s",histName),H3ConstituentPhi->GetXaxis()->FindBin(RBIN), H3ConstituentPhi->GetXaxis()->FindBin(RBIN), 0, H3ConstituentPhi->GetNbinsY());
    // ConstituentPhi = ConstituentPhi -> Rebin(nTrackptbin,"constituentphi",Trackptbin);
    if (REBIN) {
        ConstituentPhi->Rebin(2);
    }
    legend->AddEntry(ConstituentPhi, histName);
    hset(*ConstituentPhi, ConstPhiTitleX, ConstPhiTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*ConstituentPhi, Nevts, colorID, 0, 2*TMath::Pi(), 0, 0.3);
    ConstituentPhi->Draw("pesame");

    return ConstituentPhi;
}
// Draw Jets
TH1* DrawJetPt(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH1 *JetPt = (TH1 *) gROOT->FindObject(JetPtObj);
  	JetPt = JetPt->Rebin(nptBins,Form("JetPt_%s",histName),ptbin); 
    legend->AddEntry(JetPt, histName);
    hset(*JetPt, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*JetPt, Nevts, colorID, 0, 200, 1e-9, 1e-1);
    JetPt->Draw("pesame");

    return JetPt;
}
TH1* DrawJetEta(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH1 *JetEta = (TH1 *) gROOT->FindObject(JetEtaObj);
  	// JetEta = JetEta->Rebin(nptBins,Form("JetEta_%s",histName),ptbin); 
    legend->AddEntry(JetEta, histName);
    hset(*JetEta, JetEtaTitleX, JetEtaTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*JetEta, Nevts, colorID, -0.7, 0.7, 0.8, 1.5);
    JetEta->Draw("pesame");

    return JetEta;
}
TH1* DrawJetPhi(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID, Bool_t REBIN = true) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH1 *JetPhi = (TH1 *) gROOT->FindObject(JetPhiObj);
    if (REBIN) {
  	    JetPhi->Rebin(2); 
    }
    legend->AddEntry(JetPhi, histName);
    hset(*JetPhi, JetPhiTitleX, JetPhiTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*JetPhi, Nevts, colorID, 0, 2*TMath::Pi(), 0, .5);
    JetPhi->Draw("pesame");

    return JetPhi;
}
TH1* DrawJetNtracks(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH1 *JetNtracks = (TH1 *) gROOT->FindObject(JetNtracksObj);
    JetNtracks->Rebin(2); 
    legend->AddEntry(JetNtracks, histName);
    hset(*JetNtracks, JetNtracksTitleX, JetNtracksTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*JetNtracks, Nevts, colorID, 0, 40, 1e-10, 1e0);
    JetNtracks->Draw("pesame");

    return JetNtracks;
}

// TH1* DrawJetArea(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
//     auto file = TFile::Open(fileName, "open");
//     auto dir = (TDirectory*) file->Get(Dir);
//     dir->cd();
//     TH3D *H3JetArea = (TH3D *) gROOT->FindObject("h3_jet_r_jet_pt_jet_area");
//     auto JetArea = (TH1D *) H3JetArea->Project3D("zy");
//     legend->AddEntry(JetArea, histName);
//     hset(*JetPt, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
//     hoptset(*JetPt, Nevts, colorID, 0, 200, 1e-9, 1e-1);
//     JetPt->Draw("pesame");

//     return JetPt;
// }
// TH1* DrawJetResolution(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
//     auto file = TFile::Open(fileName, "open");
//     auto dir = (TDirectory*) file->Get(Dir);
//     dir->cd();
//     TH1 *JetPt = (TH1 *) gROOT->FindObject(JetPtObj);
//   	JetPt = JetPt->Rebin(nptBins,Form("JetPt_%s",histName),ptbin); 
//     legend->AddEntry(JetPt, histName);
//     hset(*JetPt, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
//     hoptset(*JetPt, Nevts, colorID, 0, 200, 1e-9, 1e-1);
//     JetPt->Draw("pesame");

//     return JetPt;
// }



// JetMatchingEfficiency
// JetMatchingPurity
// JetRM
// UnfoldedJetPt

// operate fns
void DrawHistos(const std::vector<TString>& fileNames, const std::vector<TString>& histNames, const std::vector<int>& ColorPallete) {
    std::string command = "mkdir -p plots/data/" + std::string(DatasetName);
    system(command.c_str());

    if(TrackProperties==1) {
        // Draw track pT
        Filipad2* TrackPtPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        TrackPtPad->Draw();
        TPad* trackptpad = TrackPtPad->GetPad(1); optFili(*trackptpad, 1,1,0,1);
        TPad* ratiotrackptpad = TrackPtPad->GetPad(2); optFili(*ratiotrackptpad, 1,1,0,0);
        TLegend *legtrackpt = new TLegend(0.538278,0.576812,0.811005,0.93913,NULL,"brNDC");
        legtrackpt -> SetTextSize(0.04);
        legtrackpt -> SetBorderSize(0);
        trackptpad->cd();
        TH1* TrackPtRatio = DrawTrackPt(refPath.Data(), refName, Nevents(refPath.Data(),Dir,TEventObj), legtrackpt, kBlack);
        for (size_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            trackptpad->cd(); 
            TString filePath = mainDir + fileNames[i];
            TH1* currentHist = DrawTrackPt(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,TEventObj), legtrackpt, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_tpt_%s", fileNames[i].Data());
            ratiotrackptpad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), TrackPtRatio, currentHist, TrackPtTitleX, RatioTitleY, ColorPallete[i], 0.5, 1.5);
        }
        trackptpad->cd();
        legtrackpt->Draw();
        // TrackPtPad->C->Print("TrackPt.pdf");
        TrackPtPad->C->Print(Form("plots/data/%s/TrackPt.pdf", DatasetName));
        
        // Draw track eta
        Filipad2* TrackEtaPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        TrackEtaPad->Draw();
        TPad* tracketapad = TrackEtaPad->GetPad(1); optFili(*tracketapad, 1,1,0,0);
        TPad* ratiotracketapad = TrackEtaPad->GetPad(2); optFili(*ratiotracketapad, 1,1,0,0);
        TLegend *legtracketa = new TLegend(0.538278,0.576812,0.811005,0.93913,NULL,"brNDC");
        legtracketa -> SetTextSize(0.04);
        legtracketa -> SetBorderSize(0);
        tracketapad->cd();
        TH1* TrackEtaRatio = DrawTrackEta(refPath.Data(), refName, Nevents(refPath.Data(),Dir,TEventObj), legtracketa, kBlack);
        for (size_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            tracketapad->cd(); 
            TString filePath = mainDir + fileNames[i];
            TH1* currentHist = DrawTrackEta(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,TEventObj), legtracketa, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_teta_%s", fileNames[i].Data());
            ratiotracketapad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), TrackEtaRatio, currentHist, TrackEtaTitleX, RatioTitleY, ColorPallete[i], 0.5, 1.5);
        }
        tracketapad->cd();
        legtracketa->Draw();
        TrackEtaPad->C->Print(Form("plots/data/%s/TrackEta.pdf", DatasetName));

        // Draw track phi
        Filipad2* TrackPhiPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        TrackPhiPad->Draw();
        TPad* trackphipad = TrackPhiPad->GetPad(1); optFili(*trackphipad, 1,1,0,0);
        TPad* ratiotrackphipad = TrackPhiPad->GetPad(2); optFili(*ratiotrackphipad, 1,1,0,0);
        TLegend *legtrackphi = new TLegend(0.538278,0.576812,0.811005,0.93913,NULL,"brNDC");
        legtrackphi -> SetTextSize(0.04);
        legtrackphi -> SetBorderSize(0);
        trackphipad->cd();
        TH1* TrackPhiRatio = DrawTrackPhi(refPath.Data(), refName, Nevents(refPath.Data(),Dir,TEventObj), legtrackphi, kBlack);
        for (size_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            trackphipad->cd(); 
            TString filePath = mainDir + fileNames[i];
            TH1* currentHist = DrawTrackPhi(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,TEventObj), legtrackphi, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_tphi_%s", fileNames[i].Data());
            ratiotrackphipad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), TrackPhiRatio, currentHist, TrackPhiTitleX, RatioTitleY, ColorPallete[i], 0.4, 1.2);
        }
        trackphipad->cd();
        legtrackphi->Draw();
        TrackPhiPad->C->Print(Form("plots/data/%s/TrackPhi.pdf", DatasetName));
    }

    if (ConstituentProperties==1) {
        // Draw constituents pT
        Filipad2* ConstPtPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        ConstPtPad->Draw();
        TPad* constptpad = ConstPtPad->GetPad(1); optFili(*constptpad, 1,1,0,1);
        TPad* ratioconstptpad = ConstPtPad->GetPad(2); optFili(*ratioconstptpad, 1,1,0,0);
        TLegend *legconstpt = new TLegend(0.538278,0.576812,0.811005,0.93913,NULL,"brNDC");
        legconstpt -> SetTextSize(0.04);
        legconstpt -> SetBorderSize(0);
        constptpad->cd();
        TH1* ConstPtRatio = DrawConstituentPt(refPath.Data(), refName, Nevents(refPath.Data(),JEventDir,JEventObj), legconstpt, kBlack);
        for (size_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            constptpad->cd(); 
            TString filePath = mainDir + fileNames[i];
            TH1* currentHist = DrawConstituentPt(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),JEventDir,JEventObj), legconstpt, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_cpt_%s", fileNames[i].Data());
            ratioconstptpad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), ConstPtRatio, currentHist, ConstPtTitleX, RatioTitleY, ColorPallete[i], 0.5, 1.5);
        }
        constptpad->cd();
        legconstpt->Draw();
        ConstPtPad->C->Print(Form("plots/data/%s/ConstPt_R%.1f.pdf", DatasetName, RBIN));
        
        // Draw Constituents Eta
        Filipad2* ConstEtaPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        ConstEtaPad->Draw();
        TPad* constetapad = ConstEtaPad->GetPad(1); optFili(*constetapad, 1, 1, 0 ,0);
        TPad* ratioconstetapad = ConstEtaPad->GetPad(2); optFili(*ratioconstetapad, 1, 1, 0 ,0);
        TLegend *legconsteta = new TLegend(0.538278,0.576812,0.811005,0.93913,NULL,"brNDC");
        legconsteta -> SetTextSize(0.04);
        legconsteta -> SetBorderSize(0);
        constetapad->cd();
        TH1* ConstEtaRatio = DrawConstituentEta(refPath.Data(), refName, Nevents(refPath.Data(),JEventDir,JEventObj), legconsteta, kBlack);
        for (size_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            constetapad->cd(); 
            TString filePath = mainDir + fileNames[i];
            TH1* currentHist = DrawConstituentEta(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),JEventDir,JEventObj), legconsteta, ColorPallete[i]);
            cout<<"(hist name, Nevents) = ("<< histNames[i].Data()<<", "<< Nevents(filePath.Data(),JEventDir,JEventObj)<<")"<<endl;
            TString ratioName = TString::Format("RatioHist_ceta_%s", fileNames[i].Data());
            ratioconstetapad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), ConstEtaRatio, currentHist, ConstEtaTitleX, RatioTitleY, ColorPallete[i], 0.4, 1.6);
        }
        constetapad->cd();
        legconsteta->Draw();
        ConstEtaPad->C->Print(Form("plots/data/%s/ConstEta_R%.1f.pdf", DatasetName, RBIN));

        // Draw Constituents Phi
        Filipad2* ConstPhiPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        ConstPhiPad->Draw();
        TPad* constphipad = ConstPhiPad->GetPad(1); optFili(*constphipad, 1, 1, 0 ,0);
        TPad* ratioconstphipad = ConstPhiPad->GetPad(2); optFili(*ratioconstphipad, 1, 1, 0 ,0);
        TLegend *legconstphi = new TLegend(0.538278,0.576812,0.811005,0.93913,NULL,"brNDC");
        legconstphi -> SetTextSize(0.04);
        legconstphi -> SetBorderSize(0);
        constphipad->cd();
        TH1* ConstPhiRatio = DrawConstituentPhi(refPath.Data(), refName, Nevents(refPath.Data(),JEventDir,JEventObj), legconstphi, kBlack, false);
        for (size_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            constphipad->cd(); 
            TString filePath = mainDir + fileNames[i];
            TH1* currentHist = DrawConstituentPhi(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),JEventDir,JEventObj), legconstphi, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_cphi_%s", fileNames[i].Data());
            ratioconstphipad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), ConstPhiRatio, currentHist, ConstPhiTitleX, RatioTitleY, ColorPallete[i], 0.2, 1.4);
        }
        constphipad->cd();
        legconstphi->Draw();
        ConstPhiPad->C->Print(Form("plots/data/%s/ConstPhi_R%.1f.pdf", DatasetName, RBIN));
    }
    
    if(JetProperties==1) {
        // Draw jet pT
        Filipad2* JetPtPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        JetPtPad->Draw();
        TPad* jetptpad = JetPtPad->GetPad(1); optFili(*jetptpad, 1,1,0,1);
        TPad* ratiojetptpad = JetPtPad->GetPad(2); optFili(*ratiojetptpad, 1,1,0,0);
        TLegend *legjetpt = new TLegend(0.538278,0.576812,0.811005,0.93913,NULL,"brNDC");
        legjetpt -> SetTextSize(0.04);
        legjetpt -> SetBorderSize(0);
        jetptpad->cd();
        TH1* JetPtRatio = DrawJetPt(refPath.Data(), refName, Nevents(refPath.Data(),JEventDir,JEventObj), legjetpt, kBlack);
        for (size_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            jetptpad->cd(); 
            TString filePath = mainDir + fileNames[i];
            TH1* currentHist = DrawJetPt(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),JEventDir,JEventObj), legjetpt, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_jpt_%s", fileNames[i].Data());
            ratiojetptpad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), JetPtRatio, currentHist, JetPtTitleX, RatioTitleY, ColorPallete[i], 0.5, 1.5);
        }
        jetptpad->cd();
        legjetpt->Draw();
        JetPtPad->C->Print(Form("plots/data/%s/JetPt_R%.1f.pdf", DatasetName, RBIN));
        
        // Draw jet eta
        Filipad2* JetEtaPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        JetEtaPad->Draw();
        TPad* jetetapad = JetEtaPad->GetPad(1); optFili(*jetetapad, 1,1,0,0);
        TPad* ratiojetetapad = JetEtaPad->GetPad(2); optFili(*ratiojetetapad, 1,1,0,0);
        TLegend *legjeteta = new TLegend(0.538278,0.576812,0.811005,0.93913,NULL,"brNDC");
        legjeteta -> SetTextSize(0.04);
        legjeteta -> SetBorderSize(0);
        jetetapad->cd();
        TH1* JetEtaRatio = DrawJetEta(refPath.Data(), refName, Nevents(refPath.Data(),JEventDir,JEventObj), legjeteta, kBlack);
        for (size_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            jetetapad->cd(); 
            TString filePath = mainDir + fileNames[i];
            TH1* currentHist = DrawJetEta(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),JEventDir,JEventObj), legjeteta, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_jeta_%s", fileNames[i].Data());
            ratiojetetapad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), JetEtaRatio, currentHist, JetEtaTitleX, RatioTitleY, ColorPallete[i], 0.5, 1.5);
        }
        jetetapad->cd();
        legjeteta->Draw();
        JetEtaPad->C->Print(Form("plots/data/%s/JetEta_R%.1f.pdf", DatasetName, RBIN));

        // Draw jet phi
        Filipad2* JetPhiPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        JetPhiPad->Draw();
        TPad* jetphipad = JetPhiPad->GetPad(1); optFili(*jetphipad, 1,1,0,0);
        TPad* ratiojetphipad = JetPhiPad->GetPad(2); optFili(*ratiojetphipad, 1,1,0,0);
        TLegend *legjetphi = new TLegend(0.538278,0.576812,0.811005,0.93913,NULL,"brNDC");
        legjetphi -> SetTextSize(0.04);
        legjetphi -> SetBorderSize(0);
        jetphipad->cd();
        TH1* JetPhiRatio = DrawJetPhi(refPath.Data(), refName, Nevents(refPath.Data(),JEventDir,JEventObj), legjetphi, kBlack, false);
        for (size_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            jetphipad->cd(); 
            TString filePath = mainDir + fileNames[i];
            TH1* currentHist = DrawJetPhi(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),JEventDir,JEventObj), legjetphi, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_jphi_%s", fileNames[i].Data());
            ratiojetphipad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), JetPhiRatio, currentHist, JetPhiTitleX, RatioTitleY, ColorPallete[i], 0.5, 1.5);
        }
        jetphipad->cd();
        legjetphi->Draw();
        JetPhiPad->C->Print(Form("plots/data/%s/JetPhi_R%.1f.pdf", DatasetName, RBIN));

        // Draw jet ntracks
        Filipad2* JetNtracksPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        JetNtracksPad->Draw();
        TPad* jetntrackspad = JetNtracksPad->GetPad(1); optFili(*jetntrackspad, 1,1,0,1);
        TPad* ratiojetntrackspad = JetNtracksPad->GetPad(2); optFili(*ratiojetntrackspad, 1,1,0,0);
        TLegend *legjetntracks = new TLegend(0.538278,0.576812,0.811005,0.93913,NULL,"brNDC");
        legjetntracks -> SetTextSize(0.04);
        legjetntracks -> SetBorderSize(0);
        jetntrackspad->cd();
        TH1* JetNtracksRatio = DrawJetNtracks(refPath.Data(), refName, Nevents(refPath.Data(),JEventDir,JEventObj), legjetntracks, kBlack);
        for (size_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            jetntrackspad->cd(); 
            TString filePath = mainDir + fileNames[i];
            TH1* currentHist = DrawJetNtracks(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),JEventDir,JEventObj), legjetntracks, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_jntracks_%s", fileNames[i].Data());
            ratiojetntrackspad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), JetNtracksRatio, currentHist, JetNtracksTitleX, RatioTitleY, ColorPallete[i], 0., 1.5);
        }
        jetntrackspad->cd();
        legjetntracks->Draw();
        JetNtracksPad->C->Print(Form("plots/data/%s/JetNtracks_R%.1f.pdf", DatasetName, RBIN));
    }


}
