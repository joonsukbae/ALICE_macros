#include <iostream>
#include <vector>
#include "BSHelper.cxx"
#include "Filipad2.h"
using namespace std;
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

// plots switch
auto TrackProcess = 1;
auto ConstituentProcess = 1;
auto JetProcess = 1;
auto JetMatchingProcess = 1;

enum {
  kCORRELATE = 0,
  kMCDMATCHED,
  kMCPMATCHED,
  kMATCHEDKIN,
  kFAKE,
  kMISS,
  kUNFOLDEDMC,
  kUNFOLDEDDATA,
  kEND
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
void hoptset(TH1& hid, Float_t N = 1, Color_t color = kBlack, Double_t minX = 0, Double_t maxX = 100, Double_t minY = 0, Double_t maxY = 1, Double_t MarkerSize = 1) {
    if (N == 1) {
        hid.Scale(1./hid.Integral(), "width");
    }
    else if (N==2) {
        hid.Scale(1./hid.Integral(), "");
    }
    else{
        hid.Scale(1./N, "width");
    }
 
    hid.SetMarkerColor(color);
    hid.SetLineColor(color);
    hid.SetMarkerSize(MarkerSize);
    hid.SetMarkerStyle(22);

    auto minXbin = hid.GetXaxis()->FindBin(minX);
    auto maxXbin = hid.GetXaxis()->FindBin(maxX);
    auto minYbin = hid.GetYaxis()->FindBin(minY);
    auto maxYbin = hid.GetYaxis()->FindBin(maxY);

    hid.GetXaxis()->SetRangeUser(minX,maxX);
    hid.GetYaxis()->SetRangeUser(minY,maxY);
    hid.SetFillColorAlpha(color, 0.3);
    hid.GetYaxis()->SetNdivisions(505);
}
void optFili(TPad& pid, Int_t gridx, Int_t gridy, Int_t logx, Int_t logy) {
    pid.SetGridy(gridx);
    pid.SetGridx(gridy);
    pid.SetLogx(logx);
    pid.SetLogy(logy);
}

std::vector<Double_t> NewBin(Int_t Nbins = 100, Double_t minBin = 0, Double_t maxBin = 100){
    std::vector<Double_t> Bin(Nbins + 1);
    for (Int_t i = 0; i <= Nbins; i++) {
        Bin[i] = minBin + i * (maxBin - minBin) / Nbins;
    }
    return Bin;
}
Double_t Trackptbin[15] = {0, 2, 4, 6, 8, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80};
Double_t ptbin[21] = {5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
Double_t ptbinHEP[20] = {5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140};
Int_t nTrackptbin = sizeof(Trackptbin) / sizeof(Trackptbin[0]) - 1;
Int_t nptBins = sizeof(ptbin) / sizeof(ptbin[0]) - 1;
Int_t nptBinsHEP = sizeof(ptbinHEP) / sizeof(ptbinHEP[0]) - 1;
Double_t RBIN = 0.4;
Int_t Npthat = 2;
const char *Dir = "jet-finder-charged-qa";
const char *EventObj = "h_collisions";
// TString mainDir = "../../jets/data/jetfinderQA/AnalysisResults/LHC22_apass4/R_0.2/";
TString mainDir = "../../jets/mc/jetfinderQA/AnalysisResults/";
// TString refDir = "../../jets/data/jetfinderQA/AnalysisResults/LHC22_apass4/R_0.2/LHC22f_apass4/";
TString refDir = "../../jets/data/jetfinderQA/AnalysisResults/LHC22_apass4/LHC22f_apass4/";
TString refFile = "AnalysisResults.root";
TString refPath = refDir + refFile;
const char* refName = "Raw Data (LHC22f)";
std::vector<TString> fileNames = {"AnalysisResults_LHC23d1k_HY_3.root", "AnalysisResults_LHC23d4_HY_4.root", "AnalysisResults_LHC23d4_HY_AfTrkWgt.root"};
std::vector<TString> histNames = {"Anchored MC", "Jet-MC (19 Dec)", "Jet-MC (17 Jan)"};
std::vector<int> ColorPallete = {kRed, kBlue, kGreen+1, kMagenta+1, kCyan+1, kOrange+1, kYellow+2, kAzure+7,
kViolet+1, kSpring-6, kPink+7, kTeal+3, kAzure+2, kOrange-3, kSpring+8,
kMagenta-3, kYellow-3, kRed-4, kGreen-5, kBlue-6};

const char* TrackPtObj = "h_track_pt";
const char* TrackEtaObj = "h_track_eta";
const char* TrackPhiObj = "h_track_phi";
const char* ConstPtObj = "h3_jet_r_jet_pt_track_pt";
const char* ConstEtaObj = "h3_jet_r_jet_pt_track_eta";
const char* ConstPhiObj = "h3_jet_r_jet_pt_track_phi";
const char* JetPtObj = "h_jet_pt";
const char* JetPtMCPObj = "h_jet_pt_part";
// const char* JetEtaObj = "h_jet_eta";
const char* JetEtaObj = "h3_jet_r_jet_pt_jet_eta";
// const char* JetPhiObj = "h_jet_phi";
const char* JetPhiObj = "h3_jet_r_jet_pt_jet_phi";
const char* JetNtracksObj = "h_jet_ntracks";
const char* JetAreaObj = "h3_jet_r_jet_pt_jet_area";
const char* JetResolutionObj = "h3_jet_r_jet_pt_part_jet_pt_diff";
const char* RatioTitleY = "MC / Data";
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
const char* JRETitleX = "#it{p}_{T, jet}^{truth} (GeV/c)";
const char* JRETitleY = "1/N dN/d#it{p}_{T}";
const char* JRPTitleX = "#it{p}_{T, jet}^{reco} (GeV/c)";
const char* JRPTitleY = "1/N dN/d#it{p}_{T}";
const char* JetResolutionTitleX = "#it{p}_{T, jet}^{Gen} - #it{p}_{T, jet}^{Reco} / #it{p}_{T, jet}^{Gen}";
const char* JetResolutionTitleY = "1/N dN/d#it{p}_{T}";


// declare fns
void DrawHistos(const std::vector<TString>& fileNames, const std::vector<TString>& histNames, const std::vector<int>& ColorPallete);
// TLegend* legconstpt;


void DrawJetsMC() {

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

TH1* DrawRatio(const char* ratioName, TH1* refHist, TH1* testHist, TString AxisTitleX, TString AxisTitleY, Color_t colorID, Double_t Ymin, Double_t Ymax, Double_t MarkerSize = 1) {
    TH1* ratioHist = (TH1*) refHist->Clone(ratioName);
    hset(*ratioHist, AxisTitleX, AxisTitleY, 1.2, 1.0, 0.07, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 505);
    ratioHist->Divide(testHist, ratioHist, 1., 1., "B");
    ratioHist->GetYaxis()->SetRangeUser(Ymin, Ymax);
    // ratioHist->GetYaxis()->SetLimits(Ymin, Ymax);
    ratioHist->SetMarkerColor(colorID);
    ratioHist->SetLineColor(colorID);
    ratioHist->SetMarkerSize(MarkerSize);
    ratioHist->Draw("esame");

    return ratioHist;
}
// Draw Tracks
TH1* DrawTrackPt(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH1 *TrackPt = (TH1 *) gROOT->FindObject(TrackPtObj);
  	TrackPt = TrackPt->Rebin(nTrackptbin,Form("TrackPt_%s",histName),Trackptbin); 
    legend->AddEntry(TrackPt, histName);
    hset(*TrackPt, TrackPtTitleX, TrackPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*TrackPt, Nevts, colorID, 0, 100, 1e-9, 1e0+0.05);
    TrackPt->Draw("esame");

    return TrackPt;
}
TH1* DrawTrackEta(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH1 *TrackEta = (TH1 *) gROOT->FindObject(TrackEtaObj);
  	// TrackEta = TrackEta->Rebin(nptBins,Form("TrackEta_%s",histName),ptbin); 
    legend->AddEntry(TrackEta, histName);
    hset(*TrackEta, TrackEtaTitleX, TrackEtaTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    // hoptset(*TrackEta, Nevts, colorID, -0.9, 0.9, 0, 25.);
    hoptset(*TrackEta, Nevts, colorID, -0.9, 0.9, 0.46, 0.64);
    TrackEta->Draw("esame");

    return TrackEta;
}
TH1 * DrawTrackPhi(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH1 *TrackPhi = (TH1 *) gROOT->FindObject(TrackPhiObj);
  	// TrackPhi->Rebin(10); 
  	std::vector<Double_t> newBins = NewBin(40, -1, 7);
    TrackPhi = TrackPhi->Rebin(40, Form("TrackPhi_%s", histName), newBins.data());
    legend->AddEntry(TrackPhi, histName);
    hset(*TrackPhi, TrackPhiTitleX, TrackPhiTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*TrackPhi, Nevts, colorID, 0, 2*TMath::Pi(), 0, 0.3);
    TrackPhi->Draw("esame");

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
    ConstituentPt->Draw("esame");

    return ConstituentPt;
}
TH1D* DrawConstituentEta(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH3D *H3ConstituentEta = (TH3D *) gROOT->FindObject(ConstEtaObj);
    TH1D *ConstituentEta = H3ConstituentEta->ProjectionZ(Form("ConstituentEta_%s",histName),H3ConstituentEta->GetXaxis()->FindBin(RBIN), H3ConstituentEta->GetXaxis()->FindBin(RBIN), H3ConstituentEta->GetYaxis()->FindBin(40), H3ConstituentEta->GetYaxis()->FindBin(80));
    cout<<"ConstituentEta :"<<ConstituentEta->GetEntries()<<endl;
    // TH1D *ConstituentEta = H3ConstituentEta->ProjectionZ(Form("ConstituentEta_%s",histName),H3ConstituentEta->GetXaxis()->FindBin(RBIN), H3ConstituentEta->GetXaxis()->FindBin(RBIN), 0, H3ConstituentEta->GetNbinsY());
    // ConstituentEta = ConstituentEta -> Rebin(nTrackptbin,"constituenteta",Trackptbin);
    legend->AddEntry(ConstituentEta, histName);
    hset(*ConstituentEta, ConstEtaTitleX, ConstEtaTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*ConstituentEta, Nevts, colorID, -0.5, 0.5, 0.4, 1.4);
    ConstituentEta->Draw("esame");

    return ConstituentEta;
}
TH1* DrawConstituentPhi(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH3D *H3ConstituentPhi = (TH3D *) gROOT->FindObject(ConstPhiObj);
    TH1 *ConstituentPhi = H3ConstituentPhi->ProjectionZ(Form("ConstituentPhi_%s",histName),H3ConstituentPhi->GetXaxis()->FindBin(RBIN), H3ConstituentPhi->GetXaxis()->FindBin(RBIN), 0, H3ConstituentPhi->GetNbinsY());
    // ConstituentPhi = ConstituentPhi -> Rebin(nTrackptbin,"constituentphi",Trackptbin);
    ConstituentPhi->Rebin(2);
    legend->AddEntry(ConstituentPhi, histName);
    hset(*ConstituentPhi, ConstPhiTitleX, ConstPhiTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*ConstituentPhi, Nevts, colorID, 0, 2*TMath::Pi(), 0., .4);
    ConstituentPhi->Draw("esame");

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
    hoptset(*JetPt, Nevts, colorID, 0, 200, 1e-9, 1e0);
    JetPt->Draw("esame");

    return JetPt;
}
TH1* DrawJetPtMCP(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID, Int_t i = 0) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH1 *JetPtMCP = (TH1 *) gROOT->FindObject(JetPtMCPObj);

  	// JetPtMCP = JetPtMCP->Rebin(nptBins,Form("JetPtMCP_%s",histName),ptbin); 
    legend->AddEntry(JetPtMCP, histName);
    hset(*JetPtMCP, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*JetPtMCP, Nevts, colorID, 0, 200, 1e-9, 1e0, 1 - 0.3*i);
    JetPtMCP->Draw("esame");

    return JetPtMCP;
}
TH1* DrawJetEta(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();

    TH3D *H3JetEta = (TH3D *) gROOT->FindObject(JetEtaObj); 
    TH1 *JetEta = H3JetEta->ProjectionZ(Form("JetEta_%s",histName),H3JetEta->GetXaxis()->FindBin(RBIN+1e-6), H3JetEta->GetXaxis()->FindBin(RBIN+1e-5), H3JetEta->GetYaxis()->FindBin(40), H3JetEta->GetYaxis()->FindBin(80));
    
    legend->AddEntry(JetEta, histName);
    hset(*JetEta, JetEtaTitleX, JetEtaTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*JetEta, Nevts, colorID, -0.9+RBIN, 0.9-RBIN, 0., 2.);
    JetEta->Draw("esame");

    return JetEta;
}
TH1* DrawJetPhi(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();

    TH3D *H3JetPhi = (TH3D *) gROOT->FindObject(JetPhiObj); 
    TH1D *JetPhi = H3JetPhi->ProjectionZ(Form("JetPhi_%s",histName),H3JetPhi->GetXaxis()->FindBin(RBIN+1e-6), H3JetPhi->GetXaxis()->FindBin(RBIN+1e-5), H3JetPhi->GetYaxis()->FindBin(40), H3JetPhi->GetYaxis()->FindBin(80));
    // TH1 *JetPhi = (TH1 *) gROOT->FindObject(JetPhiObj);
    if (JetPhi->GetBinWidth(0)!=0.1) {
  	    JetPhi->Rebin(4);
    }
    else {JetPhi->Rebin(2);}

    legend->AddEntry(JetPhi, histName);
    hset(*JetPhi, JetPhiTitleX, JetPhiTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*JetPhi, Nevts, colorID, 0, 2*TMath::Pi(), 0., .4);
    JetPhi->Draw("esame");

    return JetPhi;
}
TH1* DrawJetNtracks(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH1 *JetNtracks = (TH1 *) gROOT->FindObject(JetNtracksObj);
    // JetNtracks->Rebin(2); 
    legend->AddEntry(JetNtracks, histName);
    hset(*JetNtracks, JetNtracksTitleX, JetNtracksTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*JetNtracks, Nevts, colorID, 0, 40, 1e-12, 1e0);
    JetNtracks->Draw("esame");

    return JetNtracks;
}
TH1* DrawJetArea(const char* fileName, const char* histName, float Nevts, TLegend* legend, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH3F *H3JetArea = (TH3F *) gROOT->FindObject(JetAreaObj);
    auto JetArea = (TH1F *) H3JetArea->Project3D(Form("zy_%s",histName));
    legend->AddEntry("", histName, "");
    hset(*JetArea, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    JetArea->Scale(1./JetArea->Integral(), "width");
    JetArea->GetXaxis()->SetRangeUser(0., 200.);
    JetArea->GetYaxis()->SetRangeUser(0., 2.);
    JetArea->GetZaxis()->SetRangeUser(1e-7, 1e0);
    JetArea->Draw("colz");

    return JetArea;
}
TH1* DrawJetResolution(const char* fileName, const char* histName, TLegend* legend, TLegend* legend2, Color_t colorID) {
    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH3 *JetResolution = (TH3 *) gROOT->FindObject(JetResolutionObj);
    std::vector<std::pair<int, int>> pTClasses = {{5, 10}, {10, 15}, {15, 20}, {20, 40}, {40, 60}, {60, 80}, {80, 100}, {100, 140}, {140, 200}};

    legend->AddEntry("", histName, "");
    for (Int_t i = 0; i < pTClasses.size(); i++) {
        TH1D* jetResolution = (TH1D *) JetResolution->ProjectionZ(Form("JetResZ%i_%s", i+1, histName),
        JetResolution->GetXaxis()->FindBin(RBIN - 1e-6), JetResolution->GetXaxis()->FindBin(RBIN + 1e-6),
        JetResolution->GetYaxis()->FindBin(pTClasses[i].first), JetResolution->GetYaxis()->FindBin(pTClasses[i].second));
        
        jetResolution->GetXaxis()->SetTitle(JetResolutionTitleX);
        jetResolution->GetYaxis()->SetTitle(JetPtTitleY);
        legend->AddEntry(jetResolution, Form("%i< #it{p}_{T, jet}^{gen} <%i",pTClasses[i].first, pTClasses[i].second), "pe2");
        legend2->AddEntry(jetResolution, Form("Mean: %.2f, Std Dev: %.2f", jetResolution->GetMean(), jetResolution->GetRMS()), "pe2");
        // hset(*jetResolution, JetResolutionTitleX, JetPtTitleY, 0.7, 1.0, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
        hoptset(*jetResolution, 2., ColorPallete[i], -0.5, 1., 1e-5, 1e0);
        jetResolution->Draw("esame");
    }

    return 0;
}
// Draw Matched Result
std::vector<TH1*> MatchedMCPMCD;
std::vector<TH1*> JetMatching(const char* fileName, const char* DfileName, const char* histName) {
    
    auto Dfile = TFile::Open(DfileName, "open");
    auto Ddir = (TDirectory*) Dfile->Get(Dir);
    Ddir->cd();
    TH1 *DJetpt = (TH1 *) gROOT->FindObject("h_jet_pt");
    DJetpt = DJetpt->Rebin(nptBins, Form("DJetpt_%s", histName), ptbin);

    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    auto JetMCPPt = (TH1 *) gROOT->FindObject("h_jet_pt_part");
    JetMCPPt = JetMCPPt->Rebin(nptBins, Form("MCPpt_%s", histName), ptbin);
    auto JetMCDPt = (TH1 *) gROOT->FindObject("h_jet_pt");
    JetMCDPt = JetMCDPt->Rebin(nptBins, Form("MCDpt_%s", histName), ptbin);
    auto HCorrelate = (TH3 *) gROOT->FindObject("h3_jet_r_jet_pt_part_jet_pt");

    auto hcorrelate = new TH2(Form("hcorrelate_%s", histName), Form("R projected correlate_%s", histName), nptBins, ptbin, nptBins, ptbin);
    Int_t corrbin = HCorrelate->GetXaxis()->FindBin(RBIN + 1e-6);
    for (Int_t i=1; i<= HCorrelate->GetNbinsY(); i++) { // part.
        for (Int_t j=1; j<= HCorrelate->GetNbinsZ(); j++) { // det.
        Double_t content = HCorrelate->GetBinContent(corrbin, i, j);
        Int_t binpart = hcorrelate->GetYaxis()->FindBin(HCorrelate->GetYaxis()->GetBinCenter(i));
        Int_t bin = hcorrelate->GetXaxis()->FindBin(HCorrelate->GetZaxis()->GetBinCenter(j));
        Double_t currentContent = hcorrelate->GetBinContent(bin, binpart);
        Double_t newContent = currentContent + content;
        hcorrelate->SetBinContent(bin, binpart, newContent);
        }
    }

    auto *MCDMatchedpt = (TH1 *) hcorrelate->ProjectionX(Form("hMCDMatched_%s", histName), 1, hcorrelate->GetNbinsY(), "e");
    auto *MCPMatchedpt = (TH1 *) hcorrelate->ProjectionY(Form("hMCPMatched_%s", histName), 1, hcorrelate->GetNbinsX(), "e");   
    auto Respt = (TH2 *) hcorrelate->Clone();
    auto fake = (TH1*) JetMCDPt->Clone();
    fake->Add(MCDMatchedpt, -1);
    auto miss = (TH1*) JetMCPPt->Clone();
    miss->Add(MCPMatchedpt, -1);

    RooUnfoldResponse *Response = new RooUnfoldResponse(JetMCDPt, JetMCPPt);
    for (auto i = 1; i <= Respt->GetNbinsX(); i++)
    {
        for (auto j = 1; j <= Respt->GetNbinsY() ; j++)
        { //ptpair
        Double_t bincenx = Respt->GetXaxis()->GetBinCenter(i);
        Double_t binceny = Respt->GetYaxis()->GetBinCenter(j);
        Double_t bincont = Respt->GetBinContent(i, j);
        Response->Fill(bincenx, binceny, bincont);
        }
    }
    for (auto i = 1; i <= miss->GetNbinsX() ; i++)
    {
        Double_t bincenx = miss->GetXaxis()->GetBinCenter(i);
        Double_t bincont = miss->GetBinContent(i);
        Response->Miss(bincenx, bincont);
    }
    for (auto i = 1; i <= fake->GetNbinsX() ; i++)
    {
        Double_t bincenx = fake->GetXaxis()->GetBinCenter(i);
        Double_t bincont = fake->GetBinContent(i);
        Response->Fake(bincenx, bincont);
    }

    RooUnfoldBayes unfoldCon(Response, JetMCDPt, 4);
    RooUnfoldBayes unfold(Response, DJetpt, 4);
    auto hMCcorrected = (TH1 *) unfoldCon.Hreco();
    auto hDcorrected = (TH1 *) unfold.Hreco();

    
    MatchedMCPMCD.push_back(static_cast<TH1*>(JetMCDPt));
    MatchedMCPMCD.push_back(static_cast<TH1*>(hcorrelate));
    MatchedMCPMCD.push_back(static_cast<TH1*>(MCDMatchedpt));
    MatchedMCPMCD.push_back(static_cast<TH1*>(MCPMatchedpt));
    MatchedMCPMCD.push_back(static_cast<TH1*>(fake));
    MatchedMCPMCD.push_back(static_cast<TH1*>(miss));
    MatchedMCPMCD.push_back(static_cast<TH1*>(hMCcorrected));
    MatchedMCPMCD.push_back(static_cast<TH1*>(hDcorrected));

    return MatchedMCPMCD;
}
// Draw JRE
TH1F* DrawJRE(TH1F* JetMCPPt, TH1F* MCPMatchedpt, const char* histName, float Nevts, TPad* jrep, TPad* ratiojrep, TLegend* JRElegend, Color_t colorID = kBalck) {
    jrep->cd();
    auto JREp = (TH1F *) JetMCPPt->Clone();
    // JREp->Scale(1.,"width");
    JRElegend->AddEntry("", histName,"");
    JRElegend->AddEntry(JREp, "Particle level jets");
    hset(*JREp, JRETitleX, JRETitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    hoptset(*JREp, Nevts, kBlack, 0, 200, 1e-12, 1e-3);
    auto mcpmatchedpt = (TH1F *) MCPMatchedpt->Clone();
    JRElegend->AddEntry(mcpmatchedpt, "Matched jets in Particle level");
    hset(*mcpmatchedpt, JRETitleX, JRETitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    hoptset(*mcpmatchedpt, Nevts, kRed, 0, 200, 1e-12, 1e-3);
    JREp->GetYaxis()->SetNdivisions(505);
    JREp->Draw("pe");
    mcpmatchedpt->Draw("pesame");
    ratiojrep->cd();
    auto jre = (TH1F *) MCPMatchedpt -> Clone();
    jre -> Divide(jre, JetMCPPt, 1., 1., "B");
    // jre -> Divide(JREp);
    hset(*jre, JRETitleX, "JRE", 1.2, 1.0, 0.07, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 505);
    // hoptset(*jre, 1, kRed, 0, 200, 0, 1);
    // jre->Scale(1.,"width");
    jre->SetMarkerColor(kRed);
    jre->SetLineColor(kRed);
    jre->SetMarkerSize(.7);
    jre->SetMarkerStyle(22);
    jre->GetXaxis()->SetRangeUser(0., 200.);
    jre->GetYaxis()->SetRangeUser(0., 1.);
    jre->SetFillColorAlpha(kRed, 0.3);
    jre->Draw("pe");

    return 0;
}
// Draw RM
TH2F* DrawRM(TH2F* hRM,const char* histName, float Nevts, Color_t colorID = kBlack) {
    TCanvas* canRM = new TCanvas(Form("Correlation_%s", histName), Form("c%s", histName), 800, 800);
    canRM->cd();
    // canCorrelation->SetLogx(1);
    // canCorrelation->SetLogy(1);
    canRM->SetLogz(1);
    hRM->GetXaxis()->SetTitleSize(0.033); 
    hRM->GetYaxis()->SetTitleSize(0.033);
    hRM->GetXaxis()->SetLabelSize(0.025);
    hRM->GetYaxis()->SetLabelSize(0.025);
    hRM->GetZaxis()->SetLabelSize(0.025); 
    // hcorrelate->SetTitleOffset(0.1);
    hRM->SetTitle(Form("%s",histName));
    hRM->GetXaxis()->SetTitle("#it{p}_{T, jet}^{reco} (GeV/c)");
    hRM->GetYaxis()->SetTitle("#it{p}_{T, jet}^{truth} (GeV/c)");
    hRM->Scale(1./hRM->Integral(),"width");
    hRM->GetZaxis()->SetRangeUser(1e-10,1e0);
    hRM->Draw("colz");
    canRM->Print(Form("plots/2024Jan_KoALICE/JetPtCorrelation_R%.1f_%s.pdf", RBIN, histName));

    return 0;
}
// Draw Consistency Check
TH1F* DrawConsCheck(TH1F * hJetMCPPt, TH1F* hMCcorrected, const char* histName, float Nevts, TPad* consistencyp, TPad* ratioconsistp, TLegend* CONSISTlegend, Color_t colorID = kBlack) {
    consistencyp->cd();
    auto rawMCP = (TH1F *) hJetMCPPt->Clone();
    CONSISTlegend->AddEntry(rawMCP, "MC Generated");
    hset(*rawMCP, JRETitleX, JRETitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    rawMCP->Scale(1./rawMCP->Integral(), "width");
    rawMCP->SetMarkerSize(1.3);
    rawMCP->SetMarkerStyle(26);
    rawMCP->Draw("pe");
    auto UnfoldMC = (TH1F *) hMCcorrected->Clone();
    CONSISTlegend->AddEntry("", histName,"");
    CONSISTlegend->AddEntry(UnfoldMC, "Unfolded MC");
    hset(*UnfoldMC, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    hoptset(*UnfoldMC, 1., kRed, 0, 200, 1e-9, 1e0);
    UnfoldMC->GetYaxis()->SetNdivisions(505);
    UnfoldMC->Draw("pesame");
    ratioconsistp->cd();
    auto ratconsist = (TH1F *) hJetMCPPt -> Clone();
    ratconsist -> Divide(ratconsist, hMCcorrected, 1., 1., "B");
    hset(*ratconsist, JetPtTitleX, "Data", 1.2, 1.0, 0.07, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 505);
    ratconsist->SetMarkerColor();
    ratconsist->SetLineColor(kRed);
    ratconsist->SetMarkerSize(.7);
    ratconsist->SetMarkerStyle(22);
    ratconsist->GetXaxis()->SetRangeUser(0., 200.);
    ratconsist->GetYaxis()->SetRangeUser(0.8, 1.2);
    ratconsist->Draw("pe");

    return 0;
}
// Draw Unfold
TH1F* DrawUnfold(TH1F* hDcorrected, TH1F* rawMCP, const char* histName, float Nevts, TPad* unfoldp, TPad* ratunfoldp, TLegend* unfoldlegend, Color_t colorID = kBlack) {
    unfoldp->cd();
    auto UnfoldData = (TH1F *) hDcorrected->Clone();
    unfoldlegend->AddEntry("", histName,"");
    unfoldlegend->AddEntry(UnfoldData, "Unfolded Data");
    hset(*UnfoldData, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    hoptset(*UnfoldData, 1., kRed, 0, 200, 1e-9, 1e0);
    UnfoldData->GetYaxis()->SetNdivisions(505);
    UnfoldData->Draw("pe");
    rawMCP->Draw("pesame");
    ratunfoldp->cd();
    auto ratunfold = (TH1F *) rawMCP -> Clone();
    ratunfold -> Divide(ratunfold, UnfoldData, 1., 1., "B");
    unfoldlegend->AddEntry(rawMCP, "MC Generated");
    hset(*ratunfold, JetPtTitleX, "MC / Data", 1.2, 1.0, 0.07, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 505);
    ratunfold->SetMarkerColor(kBlack);
    ratunfold->SetLineColor(kBlack);
    ratunfold->SetMarkerSize(.7);
    ratunfold->SetMarkerStyle(22);
    ratunfold->GetXaxis()->SetRangeUser(0., 200.);
    ratunfold->GetYaxis()->SetRangeUser(0., 1.6);
    ratunfold->Draw("pe");
}

// Pad creation
std::map<std::string, Filipad2*> padMap;
std::map<std::string, TLegend*> legMap;
void CreatePadLegend(const std::string& padname, const std::string& legendname, double padRat = 0.4) {
    padMap[padname] = new Filipad2(++nn, 2, padRat, 100, 50, 0.7, 1, 1);
    legMap[legendname] = new Filipad2(++nn, 2, padRat, 100, 50, 0.7, 1, 1);

}
// Drawing Function
void DrawHistMac(const std::string& padName, const std::vector<TString>& histNames, TLegend* legend1, TLegend* legend2 = nullptr, TLegend* legend3 = nullptr) {
    Filipad2* pad = padMap[padName];
    pad->Draw();
    TPad* tpad1 = pad->GetPad(1); optFili(*tpad1, 1,1,0,1);
    TPad* tpadRat = pad->GetPad(2); optFili(*tpadRat, 1,1,0,1);
    TLegend *legendName1 = new TLegend(0.545455,0.715942,0.815789,0.95942,NULL,"brNDC");
    legendName1 -> SetTextSize(0.05);
    legendName1 -> SetBorderSize(0);
    tpad1->cd();
    // TH1* TrackPtRatio = DrawTrackPt(refPath.Data(), refName, Nevents(refPath.Data(),Dir,EventObj), legtrackpt, kBlack);
    TH1* TrackPtRatio = DrawTrackPt(refPath.Data(), refName, 1., legendName1, kBlack);
    for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
        tpad1->cd(); 
        TString filePath = mainDir + fileNames[i];
        // TH1* currentHist = DrawTrackPt(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,EventObj), legtrackpt, ColorPallete[i]);
        TH1* currentHist = DrawTrackPt(filePath.Data(), histNames[i].Data(), 1., legendName1, ColorPallete[i]);
        TString ratioName = TString::Format("RatioHist_tpt_%s", fileNames[i].Data());
        tpadRat->cd();
        TH1* ratioHist = DrawRatio(ratioName.Data(), TrackPtRatio, currentHist, TrackPtTitleX, RatioTitleY, ColorPallete[i], 1e-1, 1e4);
    }
    tpad1->cd();
    legendName1->Draw();
    padName->C->Print("plots/2024Jan_KoALICE/TrackPt.pdf");
}

// Operate fns
void DrawHistos(const std::vector<TString>& fileNames, const std::vector<TString>& histNames, const std::vector<int>& ColorPallete) {
    if(TrackProcess==1) {
        // Draw track pT
        Filipad2* TrackPtPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        TrackPtPad->Draw();
        TPad* trackptpad = TrackPtPad->GetPad(1); optFili(*trackptpad, 1,1,0,1);
        TPad* ratiotrackptpad = TrackPtPad->GetPad(2); optFili(*ratiotrackptpad, 1,1,0,1);
        TLegend *legtrackpt = new TLegend(0.545455,0.715942,0.815789,0.95942,NULL,"brNDC");
        legtrackpt -> SetTextSize(0.05);
        legtrackpt -> SetBorderSize(0);
        trackptpad->cd();
        // TH1* TrackPtRatio = DrawTrackPt(refPath.Data(), refName, Nevents(refPath.Data(),Dir,EventObj), legtrackpt, kBlack);
        TH1* TrackPtRatio = DrawTrackPt(refPath.Data(), refName, 1., legtrackpt, kBlack);
        for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            trackptpad->cd(); 
            TString filePath = mainDir + fileNames[i];
            // TH1* currentHist = DrawTrackPt(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,EventObj), legtrackpt, ColorPallete[i]);
            TH1* currentHist = DrawTrackPt(filePath.Data(), histNames[i].Data(), 1., legtrackpt, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_tpt_%s", fileNames[i].Data());
            ratiotrackptpad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), TrackPtRatio, currentHist, TrackPtTitleX, RatioTitleY, ColorPallete[i], 1e-1, 1e4);
        }
        trackptpad->cd();
        legtrackpt->Draw();
        TrackPtPad->C->Print("plots/2024Jan_KoALICE/TrackPt.pdf");
        
        // Draw track eta
        Filipad2* TrackEtaPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        TrackEtaPad->Draw();
        TPad* tracketapad = TrackEtaPad->GetPad(1); optFili(*tracketapad, 1,1,0,0);
        TPad* ratiotracketapad = TrackEtaPad->GetPad(2); optFili(*ratiotracketapad, 1,1,0,0);
        TLegend *legtracketa = new TLegend(0.545455,0.715942,0.815789,0.95942,NULL,"brNDC");
        legtracketa -> SetTextSize(0.05);
        legtracketa -> SetBorderSize(0);
        tracketapad->cd();
        // TH1* TrackEtaRatio = DrawTrackEta(refPath.Data(), refName, Nevents(refPath.Data(),Dir,EventObj), legtracketa, kBlack);
        TH1* TrackEtaRatio = DrawTrackEta(refPath.Data(), refName, 1., legtracketa, kBlack);
        for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            tracketapad->cd(); 
            TString filePath = mainDir + fileNames[i];
            // TH1* currentHist = DrawTrackEta(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,EventObj), legtracketa, ColorPallete[i]);
            TH1* currentHist = DrawTrackEta(filePath.Data(), histNames[i].Data(), 1., legtracketa, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_teta_%s", fileNames[i].Data());
            ratiotracketapad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), TrackEtaRatio, currentHist, TrackEtaTitleX, RatioTitleY, ColorPallete[i], 0.85, 1.15);
        }
        tracketapad->cd();
        legtracketa->Draw();
        TrackEtaPad->C->Print("plots/2024Jan_KoALICE/TrackEta.pdf");

        // Draw track phi
        Filipad2* TrackPhiPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        TrackPhiPad->Draw();
        TPad* trackphipad = TrackPhiPad->GetPad(1); optFili(*trackphipad, 1,1,0,0);
        TPad* ratiotrackphipad = TrackPhiPad->GetPad(2); optFili(*ratiotrackphipad, 1,1,0,0);
        TLegend *legtrackphi = new TLegend(0.545455,0.715942,0.815789,0.95942,NULL,"brNDC");
        legtrackphi -> SetTextSize(0.05);
        legtrackphi -> SetBorderSize(0);
        trackphipad->cd();
        // TH1* TrackPhiRatio = DrawTrackPhi(refPath.Data(), refName, Nevents(refPath.Data(),Dir,EventObj), legtrackphi, kBlack);
        TH1* TrackPhiRatio = DrawTrackPhi(refPath.Data(), refName, 1, legtrackphi, kBlack);
        for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            trackphipad->cd(); 
            TString filePath = mainDir + fileNames[i];
            // TH1* currentHist = DrawTrackPhi(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,EventObj), legtrackphi, ColorPallete[i]);
            TH1* currentHist = DrawTrackPhi(filePath.Data(), histNames[i].Data(), 1, legtrackphi, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_tphi_%s", fileNames[i].Data());
            ratiotrackphipad->cd();
            float HistMean = currentHist->GetMean(1);
            cout<< "(Dataset, Mean Y) : ("<< histNames[i].Data()<< ", " << HistMean <<")"<<endl;
            TH1* ratioHist = DrawRatio(ratioName.Data(), TrackPhiRatio, currentHist, TrackPhiTitleX, RatioTitleY, ColorPallete[i], 0.5, 1.5);
        }
        trackphipad->cd();
        legtrackphi->Draw();
        TrackPhiPad->C->Print("plots/2024Jan_KoALICE/TrackPhi.pdf");
    }

    if (ConstituentProcess==1) {
        // Draw constituents pT
        Filipad2* ConstPtPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        ConstPtPad->Draw();
        TPad* constptpad = ConstPtPad->GetPad(1); optFili(*constptpad, 1,1,0,1);
        TPad* ratioconstptpad = ConstPtPad->GetPad(2); optFili(*ratioconstptpad, 1,1,0,0);
        TLegend *legconstpt = new TLegend(0.545455,0.715942,0.815789,0.95942,NULL,"brNDC");
        legconstpt -> SetTextSize(0.05);
        legconstpt -> SetBorderSize(0);
        constptpad->cd();
        // TH1* ConstPtRatio = DrawConstituentPt(refPath.Data(), refName, Nevents(refPath.Data(),Dir,EventObj), legconstpt, kBlack);
        TH1* ConstPtRatio = DrawConstituentPt(refPath.Data(), refName, 1, legconstpt, kBlack);
        for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            constptpad->cd(); 
            TString filePath = mainDir + fileNames[i];
            // TH1* currentHist = DrawConstituentPt(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,EventObj), legconstpt, ColorPallete[i]);
            TH1* currentHist = DrawConstituentPt(filePath.Data(), histNames[i].Data(), 1, legconstpt, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_cpt_%s", fileNames[i].Data());
            ratioconstptpad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), ConstPtRatio, currentHist, ConstPtTitleX, RatioTitleY, ColorPallete[i], 0., 2.1);
        }
        constptpad->cd();
        legconstpt->Draw();
        ConstPtPad->C->Print(Form("plots/2024Jan_KoALICE/ConstPt_R%.1f.pdf",RBIN));
        
        // Draw Constituents Eta
        Filipad2* ConstEtaPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        ConstEtaPad->Draw();
        TPad* constetapad = ConstEtaPad->GetPad(1); optFili(*constetapad, 1, 1, 0 ,0);
        TPad* ratioconstetapad = ConstEtaPad->GetPad(2); optFili(*ratioconstetapad, 1, 1, 0 ,0);
        TLegend *legconsteta = new TLegend(0.545455,0.715942,0.815789,0.95942,NULL,"brNDC");
        TLegend *legconstetaPtRange = new TLegend(0.327751,0.0521739,0.598086,0.121739,NULL,"brNDC");
        legconstetaPtRange -> SetTextSize(0.05);
        legconstetaPtRange -> SetBorderSize(0);
        legconstetaPtRange->AddEntry("","40 < #it{p}_{T, jet} < 80 GeV","");
        legconsteta -> SetTextSize(0.05);
        legconsteta -> SetBorderSize(0);
        constetapad->cd();
        // TH1* ConstEtaRatio = DrawConstituentEta(refPath.Data(), refName, Nevents(refPath.Data(),Dir,EventObj), legconsteta, kBlack);
        TH1* ConstEtaRatio = DrawConstituentEta(refPath.Data(), refName, 1, legconsteta, kBlack);
        for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            constetapad->cd(); 
            TString filePath = mainDir + fileNames[i];
            // TH1* currentHist = DrawConstituentEta(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,EventObj), legconsteta, ColorPallete[i]);
            TH1* currentHist = DrawConstituentEta(filePath.Data(), histNames[i].Data(), 1, legconsteta, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_ceta_%s", fileNames[i].Data());
            ratioconstetapad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), ConstEtaRatio, currentHist, ConstEtaTitleX, RatioTitleY, ColorPallete[i], 0.4, 1.6);
        }
        constetapad->cd();
        legconsteta->Draw();
        legconstetaPtRange->Draw();
        ConstEtaPad->C->Print(Form("plots/2024Jan_KoALICE/ConstEta_R%.1f.pdf",RBIN));

        // Draw Constituents Phi
        Filipad2* ConstPhiPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        ConstPhiPad->Draw();
        TPad* constphipad = ConstPhiPad->GetPad(1); optFili(*constphipad, 1, 1, 0 ,0);
        TPad* ratioconstphipad = ConstPhiPad->GetPad(2); optFili(*ratioconstphipad, 1, 1, 0 ,0);
        TLegend *legconstphi = new TLegend(0.545455,0.715942,0.815789,0.95942,NULL,"brNDC");
        TLegend *legconstphiPtRange = new TLegend(0.327751,0.0521739,0.598086,0.121739,NULL,"brNDC");
        legconstphiPtRange -> SetTextSize(0.05);
        legconstphiPtRange -> SetBorderSize(0);
        legconstphiPtRange->AddEntry("","40 < #it{p}_{T, jet} < 80 GeV","");
        legconstphi -> SetTextSize(0.05);
        legconstphi -> SetBorderSize(0);
        constphipad->cd();
        // TH1* ConstPhiRatio = DrawConstituentPhi(refPath.Data(), refName, Nevents(refPath.Data(),Dir,EventObj), legconstphi, kBlack);
        TH1* ConstPhiRatio = DrawConstituentPhi(refPath.Data(), refName, 1, legconstphi, kBlack);
        for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            constphipad->cd(); 
            TString filePath = mainDir + fileNames[i];
            // TH1* currentHist = DrawConstituentPhi(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,EventObj), legconstphi, ColorPallete[i]);
            TH1* currentHist = DrawConstituentPhi(filePath.Data(), histNames[i].Data(), 1, legconstphi, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_cphi_%s", fileNames[i].Data());
            ratioconstphipad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), ConstPhiRatio, currentHist, ConstPhiTitleX, RatioTitleY, ColorPallete[i], 0.2, 1.4);
        }
        constphipad->cd();
        legconstphi->Draw();
        legconstphiPtRange->Draw();
        ConstPhiPad->C->Print(Form("plots/2024Jan_KoALICE/ConstPhi_R%.1f.pdf",RBIN));
    }
    
    if(JetProcess==1) {
        // Draw jet pT MCD
        Filipad2* JetPtPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        JetPtPad->Draw();
        TPad* jetptpad = JetPtPad->GetPad(1); optFili(*jetptpad, 1,1,0,1);
        TPad* ratiojetptpad = JetPtPad->GetPad(2); optFili(*ratiojetptpad, 1,1,0,0);
        TLegend *legjetpt = new TLegend(0.545455,0.715942,0.815789,0.95942,NULL,"brNDC");
        legjetpt -> SetTextSize(0.05);
        legjetpt -> SetBorderSize(0);
        jetptpad->cd();
        // TH1* JetPtRatio = DrawJetPt(refPath.Data(), refName, Nevents(refPath.Data(),Dir,EventObj), legjetpt, kBlack);
        TH1* JetPtRatio = DrawJetPt(refPath.Data(), refName, 1, legjetpt, kBlack);
        // TH1* JetPtRatio = DrawJetPt(refPath.Data(), refName, 1., legjetpt, kBlack);
        for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            jetptpad->cd(); 
            TString filePath = mainDir + fileNames[i];
            // TH1* currentHist = DrawJetPt(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,EventObj), legjetpt, ColorPallete[i]);
            TH1* currentHist = DrawJetPt(filePath.Data(), histNames[i].Data(), 1, legjetpt, ColorPallete[i]);
            // TH1* currentHist = DrawJetPt(filePath.Data(), histNames[i].Data(),1., legjetpt, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_jpt_%s", fileNames[i].Data());
            ratiojetptpad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), JetPtRatio, currentHist, JetPtTitleX, RatioTitleY, ColorPallete[i], 0., 1.5);
        }
        jetptpad->cd();
        legjetpt->Draw();
        JetPtPad->C->Print(Form("plots/2024Jan_KoALICE/JetPt_R%.1f.pdf",RBIN));

        //Draw jet pT MCP
        Filipad2* JetPtMCPPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        JetPtMCPPad->Draw();
        TPad* jetptmcppad = JetPtMCPPad->GetPad(1); optFili(*jetptmcppad, 1,1,0,1);
        TPad* ratiojetptmcppad = JetPtMCPPad->GetPad(2); optFili(*ratiojetptmcppad, 1,1,0,0);
        TLegend *legjetptmcp = new TLegend(0.545455,0.715942,0.815789,0.95942,NULL,"brNDC");
        legjetptmcp -> SetTextSize(0.05);
        legjetptmcp -> SetBorderSize(0);
        jetptmcppad->cd();
        TString JetPtMCPfilePath = mainDir + fileNames[0];
        TH1* JetPtMCPRatio = DrawJetPtMCP(JetPtMCPfilePath.Data(), histNames[2].Data(), Nevents(JetPtMCPfilePath.Data(),Dir,EventObj), legjetptmcp, kBlack);
        // TH1* JetPtMCPRatio = DrawJetPtMCP(JetPtMCPfilePath.Data(), histNames[2].Data(), 1, legjetptmcp, kBlack);
        // TH1* JetPtRatio = DrawJetPt(refPath.Data(), refName, 1., legjetpt, kBlack);
        for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            jetptmcppad->cd(); 
            TString filePath = mainDir + fileNames[i];
            TH1* currentHist = DrawJetPtMCP(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,EventObj), legjetptmcp, ColorPallete[i], i);
            // TH1* currentHist = DrawJetPtMCP(filePath.Data(), histNames[i].Data(), 1, legjetptmcp, ColorPallete[i], i);
            TString ratioName = TString::Format("RatioHist_jptmcp_%s", fileNames[i].Data());
            ratiojetptmcppad->cd();
            if (i>=1) {
                TH1* ratioHist = DrawRatio(ratioName.Data(), JetPtMCPRatio, currentHist, JetPtTitleX, "MC / Jet-MC (17 Jan)", ColorPallete[i], 0., 1.5, 1-0.3*i);
            }
        }
        jetptmcppad->cd();
        legjetptmcp->Draw();
        JetPtMCPPad->C->Print(Form("plots/2024Jan_KoALICE/JetPtMCP_R%.1f.pdf",RBIN));

    // cout<<Form("Debugger_%i",n++)<<endl;
        // Draw jet eta
        Filipad2* JetEtaPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        JetEtaPad->Draw();
        TPad* jetetapad = JetEtaPad->GetPad(1); optFili(*jetetapad, 1,1,0,0);
        TPad* ratiojetetapad = JetEtaPad->GetPad(2); optFili(*ratiojetetapad, 1,1,0,0);
        TLegend *legjeteta = new TLegend(0.545455,0.715942,0.815789,0.95942,NULL,"brNDC");
        TLegend *legjetetaPtRange = new TLegend(0.327751,0.0521739,0.598086,0.121739,NULL,"brNDC");
        legjetetaPtRange -> SetTextSize(0.05);
        legjetetaPtRange -> SetBorderSize(0);
        legjetetaPtRange->AddEntry("","40 < #it{p}_{T, jet} < 80 GeV","");
        legjeteta -> SetTextSize(0.05);
        legjeteta -> SetBorderSize(0);
        jetetapad->cd();
        // TH1* JetEtaRatio = DrawJetEta(refPath.Data(), refName, Nevents(refPath.Data(),Dir,EventObj), legjeteta, kBlack);
        TH1* JetEtaRatio = DrawJetEta(refPath.Data(), refName, 1., legjeteta, kBlack);
        for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            jetetapad->cd(); 
            TString filePath = mainDir + fileNames[i];
            // TH1* currentHist = DrawJetEta(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,EventObj), legjeteta, ColorPallete[i]);
            TH1* currentHist = DrawJetEta(filePath.Data(), histNames[i].Data(), 1., legjeteta, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_jeta_%s", fileNames[i].Data());
            ratiojetetapad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), JetEtaRatio, currentHist, JetEtaTitleX, RatioTitleY, ColorPallete[i], 0., 2.);
        }
        jetetapad->cd();
        legjeteta->Draw();
        legjetetaPtRange->Draw();
        JetEtaPad->C->Print(Form("plots/2024Jan_KoALICE/JetEta_R%.1f.pdf",RBIN));

        // Draw jet phi
        Filipad2* JetPhiPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        JetPhiPad->Draw();
        TPad* jetphipad = JetPhiPad->GetPad(1); optFili(*jetphipad, 1,1,0,0);
        TPad* ratiojetphipad = JetPhiPad->GetPad(2); optFili(*ratiojetphipad, 1,1,0,0);
        TLegend *legjetphi = new TLegend(0.545455,0.715942,0.815789,0.95942,NULL,"brNDC");
        TLegend *legjetphiPtRange = new TLegend(0.327751,0.0521739,0.598086,0.121739,NULL,"brNDC");
        legjetphiPtRange -> SetTextSize(0.05);
        legjetphiPtRange -> SetBorderSize(0);
        legjetphiPtRange->AddEntry("","40 < #it{p}_{T, jet} < 80 GeV","");
        legjetphi -> SetTextSize(0.05);
        legjetphi -> SetBorderSize(0);
        jetphipad->cd();
        // TH1* JetPhiRatio = DrawJetPhi(refPath.Data(), refName, Nevents(refPath.Data(),Dir,EventObj), legjetphi, kBlack);
        TH1* JetPhiRatio = DrawJetPhi(refPath.Data(), refName, 1., legjetphi, kBlack);
        for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            jetphipad->cd(); 
            TString filePath = mainDir + fileNames[i];
            // TH1* currentHist = DrawJetPhi(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,EventObj), legjetphi, ColorPallete[i]);
            TH1* currentHist = DrawJetPhi(filePath.Data(), histNames[i].Data(), 1., legjetphi, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_jphi_%s", fileNames[i].Data());
            ratiojetphipad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), JetPhiRatio, currentHist, JetPhiTitleX, RatioTitleY, ColorPallete[i], 0., 2.);
        }
        jetphipad->cd();
        legjetphi->Draw();
        legjetphiPtRange->Draw();
        JetPhiPad->C->Print(Form("plots/2024Jan_KoALICE/JetPhi_R%.1f.pdf",RBIN));

        // Draw jet ntracks
        Filipad2* JetNtracksPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1 );
        JetNtracksPad->Draw();
        TPad* jetntrackspad = JetNtracksPad->GetPad(1); optFili(*jetntrackspad, 1,1,0,1);
        TPad* ratiojetntrackspad = JetNtracksPad->GetPad(2); optFili(*ratiojetntrackspad, 1,1,0,0);
        TLegend *legjetntracks = new TLegend(0.545455,0.715942,0.815789,0.95942,NULL,"brNDC");
        legjetntracks -> SetTextSize(0.05);
        legjetntracks -> SetBorderSize(0);
        jetntrackspad->cd();
        TH1* JetNtracksRatio = DrawJetNtracks(refPath.Data(), refName, 1., legjetntracks, kBlack);
        for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            jetntrackspad->cd(); 
            TString filePath = mainDir + fileNames[i];
            TH1* currentHist = DrawJetNtracks(filePath.Data(), histNames[i].Data(), 1., legjetntracks, ColorPallete[i]);
            TString ratioName = TString::Format("RatioHist_jntracks_%s", fileNames[i].Data());
            ratiojetntrackspad->cd();
            TH1* ratioHist = DrawRatio(ratioName.Data(), JetNtracksRatio, currentHist, JetNtracksTitleX, RatioTitleY, ColorPallete[i], 0., 3.1);
        }
        jetntrackspad->cd();
        legjetntracks->Draw();
        JetNtracksPad->C->Print(Form("plots/2024Jan_KoALICE/JetNtracks_R%.1f.pdf",RBIN));

        // Draw jet area
        for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            TCanvas* canJetArea = new TCanvas(Form("JetArea_%s", histNames[i].Data()), Form("jetarea%s", histNames[i].Data()), 600, 800);
            canJetArea->SetLogz(1);
            canJetArea->Draw();
            TLegend *legjetarea = new TLegend(0.545455,0.715942,0.815789,0.95942,NULL,"brNDC");
            legjetarea -> SetTextSize(0.04);
            legjetarea -> SetBorderSize(0);
            TString filePath = mainDir + fileNames[i];
            // TH1* currentHist = DrawJetArea(filePath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,EventObj), legjetarea, ColorPallete[i]);
            TH1* currentHist = DrawJetArea(filePath.Data(), histNames[i].Data(), 1., legjetarea, ColorPallete[i]);
            legjetarea->Draw();
            canJetArea->Print(Form("plots/2024Jan_KoALICE/JetArea_R%.1f_%s.pdf",RBIN, histNames[i].Data()));
        }
    }

    if (JetMatchingProcess==1) {
        for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            TString filePath = mainDir + fileNames[i];

            // making JRE, Consistency Check, Unfolded Data results
            Filipad2* JREPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
            JREPad->Draw();
            TPad* jrepad = JREPad->GetPad(1); optFili(*jrepad, 1, 1, 0, 1);
            TPad* ratiojrepads = JREPad->GetPad(2); optFili(*ratiojrepads, 1, 1, 0, 0);
            TLegend* legjre = new TLegend(0.401914,0.715942,0.667464,0.95942,NULL,"brNDC");
            legjre -> SetTextSize(0.05);
            legjre -> SetBorderSize(0);
            DrawJetMatching(filePath.Data(), refPath.Data(), histNames[i].Data())[k];
            jrepad->cd();
            legjre->Draw();
            JREPad->C->Print(Form("plots/2024Jan_KoALICE/JRE_R%.1f_%s_ITS%i.pdf", RBIN, histNames[i].Data(), selITS));

            Filipad2* JRPPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
            JRPPad->Draw();
            TPad* jrppad = JRPPad->GetPad(1); optFili(*jrppad, 1, 1, 0, 1);
            TPad* ratiojrppads = JRPPad->GetPad(2); optFili(*ratiojrppads, 1, 1, 0, 0);
            TLegend* legjrp = new TLegend(0.401914,0.715942,0.667464,0.95942,NULL,"brNDC");
            legjrp -> SetTextSize(0.05);
            legjrp -> SetBorderSize(0);

            Filipad2* ConsistencyPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
            ConsistencyPad->Draw();
            TPad* consistpad = ConsistencyPad->GetPad(1); optFili(*consistpad, 1, 1, 0, 1);
            TPad* ratioconsistpad = ConsistencyPad->GetPad(2); optFili(*ratioconsistpad, 1, 1, 0, 0);
            TLegend* legconsist = new TLegend(0.401914,0.715942,0.667464,0.95942,NULL,"brNDC");
            legconsist -> SetTextSize(0.05);
            legconsist -> SetBorderSize(0);

            Filipad2* UnfoldPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
            UnfoldPad->Draw();
            TPad* unfoldpad = UnfoldPad->GetPad(1); optFili(*unfoldpad, 1, 1, 0, 1);
            TPad* ratiounfoldpad = UnfoldPad->GetPad(2); optFili(*ratiounfoldpad, 1, 1, 0, 0);
            TLegend* legunfold = new TLegend(0.619617,0.414493,0.885167,0.657971,NULL,"brNDC");
            legunfold -> SetTextSize(0.05);
            legunfold -> SetBorderSize(0);
            TLegend *legpt3 = new TLegend(0.184211,0.04,0.232057,0.25913,NULL,"brNDC");
            legpt3 -> SetTextSize(0.062);
            legpt3 -> SetBorderSize(0);
            legpt3->SetTextAlign(12);
            legpt3->AddEntry("","|#it{#eta}_{jet}| < 0.5","");
            legpt3->AddEntry("","Anti-#it{k}_{T}, #it{R} = 0.4","");
            TLegend *legpt2 = new TLegend(0.287081,0.698551,0.368421,0.947826,NULL,"brNDC");
            legpt2 -> SetTextSize(0.062);
            legpt2 -> SetBorderSize(0);
            legpt2->SetTextAlign(12);
            legpt2->AddEntry("","pp #sqrt{#it{s}} = 13.6 TeV","");
            legpt2->AddEntry("","#it{p}_{T, track} > 0.15 GeV/#it{c}","");
            legpt2->AddEntry("","|#it{#eta}_{track}| < 0.9","");

            TString filePath = mainDir + fileNames[i];
            DrawJetMatching(filePath.Data(), refPath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,EventObj), jrepad, ratiojrepads, legjre, consistpad, ratioconsistpad, legconsist, unfoldpad, ratiounfoldpad, legunfold, ColorPallete[i]);

            jrepad->cd();
            legjre->Draw();
            JREPad->C->Print(Form("plots/2024Jan_KoALICE/JRE_R%.1f_%s.pdf", RBIN, histNames[i].Data()));

            consistpad->cd();
            legconsist->Draw();
            ConsistencyPad->C->Print(Form("plots/2024Jan_KoALICE/ConsistencyCheck_R%.1f_%s.pdf", RBIN, histNames[i].Data())); 

            unfoldpad->cd();
            legpt2->Draw();
            legpt3->Draw();
            legunfold->Draw();
            UnfoldPad->C->Print(Form("plots/2024Jan_KoALICE/UnfoldedJetPt_R%.1f_%s.pdf", RBIN, histNames[i].Data())); 
        }

        for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {
            TCanvas* canJetResolution = new TCanvas(Form("JetResolution_%s", histNames[i].Data()), Form("jetresolution%s", histNames[i].Data()), 800, 600);
            canJetResolution->SetLogy(1);
            canJetResolution->SetGridx(1);
            canJetResolution->SetGridy(1);
            canJetResolution->Draw();
            TLegend *legjetresolution = new TLegend(0.578947,0.654783,0.879699,0.988696,NULL,"brNDC");
            legjetresolution -> SetTextSize(0.027);
            legjetresolution -> SetBorderSize(0);
            TLegend *legjetresolutionStat = new TLegend(0.302005,0.113043,0.616541,0.330435,NULL,"brNDC");
            legjetresolutionStat -> SetTextSize(0.03);
            legjetresolutionStat -> SetBorderSize(0);
            TString filePath = mainDir + fileNames[i];
            DrawJetResolution(filePath.Data(), histNames[i].Data(), legjetresolution, legjetresolutionStat, ColorPallete[i]);
            legjetresolution->Draw();
            legjetresolutionStat->Draw();
            canJetResolution->Print(Form("plots/2024Jan_KoALICE/JetResolution_R%.1f_%s.pdf",RBIN, histNames[i].Data()));
        }



    }

}
