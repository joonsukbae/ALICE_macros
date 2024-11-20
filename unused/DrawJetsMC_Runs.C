#include <iostream>
#include <vector>
#include "BSHelper.cxx"
#include "Filipad2.h"
using namespace std;
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

// plots switch
auto TrackProcess = 0;
auto ConstituentProcess = 0;
auto JetProcess = 0;
auto JetMatchingProcess = 1;

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
// Double_t RBIN = 0.2;
Double_t RBIN = 0.4;
Int_t selITS = 3;
// Double_t RBIN = 0.6;
Int_t Npthat = 2;
const char *Dir = "jet-finder-charged-qa";
const char *EventObj = "h_collisions";
// TString mainDir = "../../jets/data/jetfinderQA/AnalysisResults/LHC22_apass4/R_0.2/";
TString mainDir = "../../jets/mc/jetfinderQA/AnalysisResults/LHC23k2d/";
// TString refDir = "../../jets/data/jetfinderQA/AnalysisResults/LHC22_apass4/R_0.2/LHC22f_apass4/";
TString refDir = "../../jets/data/jetfinderQA/AnalysisResults/LHC22_apass4/LHC22f_apass4/";
TString refFile = "AnalysisResults.root";
TString refPath = refDir + refFile;
const char* refName = "Raw Data (LHC22f)";

std::vector<std::pair<TString, TString>> getFileAndHistNames(const TString& directory) {
    std::vector<std::pair<TString, TString>> fileAndHistNames;
    DIR* dir = opendir(directory.Data());
    if (!dir) {
        std::cerr << "Cannot open directory: " << directory << std::endl;
        return fileAndHistNames;
    }

    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        if (entry->d_type == DT_DIR && TString(entry->d_name).IsDigit()) {
            TString subDirName = entry->d_name;
            TString filePath = directory + "/" + subDirName + "/AnalysisResults.root";
            if (gSystem->AccessPathName(filePath, kFileExists) == 0) { // 파일 존재 확인
                fileAndHistNames.push_back(std::make_pair(filePath, subDirName)); // 파일 경로와 실행 번호 추가
            }
        }
    }
    closedir(dir);
    return fileAndHistNames;
}


std::vector<int> ColorPallete = {
    kRed, kBlue, kGreen+1, kMagenta+1, kCyan+1, kOrange+1, kYellow+2, kAzure+7,
    kViolet+1, kSpring-6, kPink+7, kTeal+3, kAzure+2, kOrange-3, kSpring+8,
    kMagenta-3, kYellow-3, kRed-4, kGreen-5, kBlue-6, kOrange+2, kSpring+9, 
    kTeal+1, kAzure+10, kViolet+2, kPink+10, kYellow-7, kRed+1, kGreen+4, 
    kBlue+2, kMagenta+4, kCyan+3, kOrange+7, kSpring-2, kTeal-7, kAzure-4, 
    kViolet-9, kPink+3, kYellow+5, kRed+3, kGreen-8, kBlue+10, kMagenta-6, 
    kCyan-9, kOrange+9, kSpring+4, kTeal+5, kAzure+1, kViolet+10, kPink+1
};

const char* JetPtObj = "h_jet_pt";
const char* JetPtMCPObj = "h_jet_pt_part";
const char* JetEtaObj = "h3_jet_r_jet_pt_jet_eta";
const char* JetPhiObj = "h3_jet_r_jet_pt_jet_phi";
const char* JetNtracksObj = "h_jet_ntracks";


const char* RatioTitleY = "MC / Data";
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


// declare fns
void DrawHistos(const std::vector<TString>& fileNames, const std::vector<TString>& histNames, const std::vector<int>& ColorPallete);
// TLegend* legconstpt;


void DrawJetsMC_Runs() {
    auto fileAndHistNames = getFileAndHistNames(mainDir);

    std::vector<TString> fileNames;
    std::vector<TString> histNames;
    for (const auto& pair : fileAndHistNames) {
        fileNames.push_back(pair.first);  // 파일 경로 저장
        histNames.push_back(pair.second); // 실행 번호를 히스토그램 이름으로 사용
    }

    // 테스트 출력 및 DrawHistos 함수 호출
    for (size_t i = 0; i < fileNames.size(); ++i) {
        std::cout << "File: " << fileNames[i] << ", Histogram: " << histNames[i] << std::endl;
    }
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
// Draw Jet Matchings
std::vector<TH1F*> UnfoldedHistos;
std::vector<TH1F*> DrawJetMatching(const char* fileName, const char* DfileName, const char* histName, float Nevts, TPad* jrep, TPad* ratiojrep, TLegend* JRElegend, TPad* consistencyp, TPad* ratioconsistp, TLegend* CONSISTlegend, TPad* unfoldp, TPad* ratunfoldp, TLegend* unfoldlegend, Color_t colorID, TPad* jrpp = nullptr, TPad* ratiojrpp = nullptr, TLegend* JRPlegend = nullptr) {
    
    auto Dfile = TFile::Open(DfileName, "open");
    auto Ddir = (TDirectory*) Dfile->Get(Dir);
    Ddir->cd();
    TH1 *DJetpt = (TH1 *) gROOT->FindObject("h_jet_pt");
    DJetpt = DJetpt->Rebin(nptBins, Form("DJetpt_%s", histName), ptbin);
    // cout<<"Nbins DJetpt: "<<DJetpt->GetNbinsX()<<endl;

    auto file = TFile::Open(fileName, "open");
    auto dir = (TDirectory*) file->Get(Dir);
    dir->cd();
    TH1 *JetMCPPt = (TH1 *) gROOT->FindObject("h_jet_pt_part");
    JetMCPPt = JetMCPPt->Rebin(nptBins, Form("MCPpt_%s", histName), ptbin);
    TH1 *JetMCDPt = (TH1 *) gROOT->FindObject("h_jet_pt");
    // cout<<"Nbins JetMCDPt: "<<JetMCDPt->GetNbinsX()<<endl;
    JetMCDPt = JetMCDPt->Rebin(nptBins, Form("MCDpt_%s", histName), ptbin);
    // cout<<"Nbins JetMCDPt: "<<JetMCDPt->GetNbinsX()<<endl;
    TString histNameStr(histName);
    TH3F* HCorrelate;
    // if (histNameStr.Contains("ITS")||histNameStr.Contains("hy")) {
      HCorrelate = (TH3F *) gROOT->FindObject("h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo");
    // } else {
    //   HCorrelate = (TH3F *) gROOT->FindObject("h3_jet_r_jet_pt_part_jet_pt");
    // }
    cout<<"Nentries of 3DHcorrelate: "<<HCorrelate->GetEntries()<<endl;

    TH2F *hcorrelate = new TH2F(Form("hcorrelate_%s", histName), Form("R projected correlate_%s", histName), nptBins, ptbin, nptBins, ptbin);
    Int_t corrbin = HCorrelate->GetXaxis()->FindBin(RBIN + 1e-6);
    for (Int_t i=1; i<= HCorrelate->GetNbinsY(); i++) { // part.
        for (Int_t j=1; j<= HCorrelate->GetNbinsZ(); j++) { // det.
        Double_t content = HCorrelate->GetBinContent(corrbin, i, j);
        Int_t binpart = hcorrelate->GetYaxis()->FindBin(HCorrelate->GetYaxis()->GetBinCenter(i));
        Int_t bin = hcorrelate->GetXaxis()->FindBin(HCorrelate->GetZaxis()->GetBinCenter(j));
        Double_t currentContent = hcorrelate->GetBinContent(bin, binpart);
        Double_t newContent = currentContent + content;
        // cout<< "Content in hcorrelate: ("<< bin <<","<< binpart<<","<<newContent<< ")"<<endl;
        hcorrelate->SetBinContent(bin, binpart, newContent);
        }
    }
    cout<<"Nentries of 2Dhcorrelate: "<<hcorrelate->GetBinContent(20,20)<<endl;

    TH1 *MCDMatchedpt = (TH1F *) hcorrelate->ProjectionX(Form("hMCDMatched_%s", histName), 1, hcorrelate->GetNbinsY(), "e");
    TH1 *MCPMatchedpt = (TH1F *) hcorrelate->ProjectionY(Form("hMCPMatched_%s", histName), 1, hcorrelate->GetNbinsX(), "e");   

    TH2F* Respt = (TH2F *) hcorrelate->Clone();
    TH1F* fake = (TH1F*) JetMCDPt->Clone();
    // cout<<"fake: "<<fake->GetBinContent(10)<<endl;
    fake->Add(MCDMatchedpt, -1);
    // cout<<"fake: "<<fake->GetBinContent(10)<<endl;
    TH1F* miss = (TH1F*) JetMCPPt->Clone();
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
    auto hMCcorrected = (TH1F *) unfoldCon.Hreco();
    auto hDcorrected = (TH1F *) unfold.Hreco();

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
    hoptset(*mcpmatchedpt, Nevts, colorID, 0, 200, 1e-12, 1e-3);
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
    jre->SetMarkerColor(colorID);
    jre->SetLineColor(colorID);
    jre->SetMarkerSize(.4);
    jre->SetMarkerStyle(22);
    jre->GetXaxis()->SetRangeUser(0., 200.);
    jre->GetYaxis()->SetRangeUser(0., 1.);
    jre->SetFillColorAlpha(colorID, 0.3);
    jre->Draw("pe");


    jrpp->cd();
    auto JRPp = (TH1F *) JetMCDPt->Clone();
    // JREp->Scale(1.,"width");
    JRPlegend->AddEntry("", histName,"");
    JRPlegend->AddEntry(JRPp, "Detector level jets");
    hset(*JRPp, JRPTitleX, JRPTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    hoptset(*JRPp, Nevts, kBlack, 0, 200, 1e-12, 1e-3);
    auto mcdmatchedpt = (TH1F *) MCDMatchedpt->Clone();
    JRPlegend->AddEntry(mcdmatchedpt, "Matched jets in Detector level");
    hset(*mcdmatchedpt, JRPTitleX, JRPTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
    hoptset(*mcdmatchedpt, Nevts, kRed, 0, 200, 1e-12, 1e-3);
    JRPp->GetYaxis()->SetNdivisions(505);
    JRPp->Draw("pe");
    mcdmatchedpt->Draw("pesame");
    ratiojrpp->cd();
    auto jrp = (TH1F *) MCDMatchedpt -> Clone();
    jrp -> Divide(jrp, JetMCDPt, 1., 1., "B");
    // jre -> Divide(JRPp);
    hset(*jrp, JRPTitleX, "JRP", 1.2, 1.0, 0.07, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 505);
    // hoptset(*jrp, 1, kRed, 0, 200, 0, 1);
    // jre->Scale(1.,"width");
    jrp->SetMarkerColor(kRed);
    jrp->SetLineColor(kRed);
    jrp->SetMarkerSize(.7);
    jrp->SetMarkerStyle(22);
    jrp->GetXaxis()->SetRangeUser(0., 200.);
    jrp->GetYaxis()->SetRangeUser(0., 1.);
    jrp->SetFillColorAlpha(kRed, 0.3);
    jrp->Draw("pe");

    return UnfoldedHistos;
}


// UnfoldedJetPt

// operate fns
void DrawHistos(const std::vector<TString>& fileNames, const std::vector<TString>& histNames, const std::vector<int>& ColorPallete) {
    if (JetMatchingProcess==1) {
            // making JRE, Consistency Check, Unfolded Data results
            Filipad2* JREPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
            JREPad->Draw();
            TPad* jrepad = JREPad->GetPad(1); optFili(*jrepad, 1, 1, 0, 1);
            TPad* ratiojrepads = JREPad->GetPad(2); optFili(*ratiojrepads, 1, 1, 0, 0);
            TLegend* legjre = new TLegend(0.401914,0.715942,0.667464,0.95942,NULL,"brNDC");
            legjre -> SetTextSize(0.05);
            legjre -> SetBorderSize(0);

            Filipad2* JRPPad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
            JRPPad->Draw();
            TPad* jrppad = JRPPad->GetPad(1); optFili(*jrppad, 1, 1, 0, 1);
            TPad* ratiojrppads = JRPPad->GetPad(2); optFili(*ratiojrppads, 1, 1, 0, 0);
            TLegend* legjrp = new TLegend(0.401914,0.715942,0.667464,0.95942,NULL,"brNDC");
            legjrp -> SetTextSize(0.05);
            legjrp -> SetBorderSize(0);
        
        for (Int_t i = 0; i < std::min(fileNames.size(), histNames.size()); ++i) {


            TString filePath =  fileNames[i];
            DrawJetMatching(filePath.Data(), refPath.Data(), histNames[i].Data(), Nevents(filePath.Data(),Dir,EventObj), jrepad, ratiojrepads, legjre, 0, 0, 0, 0, 0, 0, ColorPallete[i], jrppad, ratiojrppads, legjrp);

         }
            jrepad->cd();
            legjre->Draw();
            JREPad->C->Print(Form("plots/2024Jan_KoALICE/JRE_R%.1f.pdf", RBIN));

            jrppad->cd();
            legjrp->Draw();
            JRPPad->C->Print(Form("plots/2024Jan_KoALICE/JRP_R%.1f.pdf", RBIN));

    }

}
