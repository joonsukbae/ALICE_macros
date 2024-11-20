///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Draw macro for jet QA systematics //////
////////// author: Joonsuk Bae               ////// 
////////// E-mail: jbae@cern.ch              //////
////////// Last Modified: 19 July 2024       //////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

#include "../DrawJetsMCFilesTitles.h"
#include <__config>

TString HistDirectory = "jet-finder-charged-qa";

void GetTrackErr(TString defResult, TString compResult, TString outputName);
void GetConstErr(TString defResult, TString compResult, TString outputName);
void GetJetErr(TString defResult, TString compResult, TString outputName);
void GetSystematicError(TString stdResult, TString systResult);
TH1 *DrawSystematics(TH1 *hstd, TH1 *hsyst, TString FileName, const char *XaxisTitle, const char *YaxisTitle);

void GetTrackErr(TString defResult, TString compResult, TString outputName) {
    std::cout << "Opening default result file: " << defResult << std::endl;
    auto defFile = TFile::Open(defResult, "read");
    if (!defFile || defFile->IsZombie()) {
        std::cerr << "Error: Failed to open default result file." << std::endl;
        return;
    }

    std::cout << "Opening comparison result file: " << compResult << std::endl;
    auto compFile = TFile::Open(compResult, "read");
    if (!compFile || compFile->IsZombie()) {
        std::cerr << "Error: Failed to open comparison result file." << std::endl;
        return;
    }

    std::cout << "Accessing default track histogram..." << std::endl;
    auto h2DefTrack = (TH2 *)defFile->Get(Form("%s/%s", HistDirectory.Data(), TrackPtObj));
    if (!h2DefTrack) {
        std::cerr << "Error: Failed to retrieve default track histogram." << std::endl;
        return;
    }
    std::cout << "Successfully retrieved default track histogram." << std::endl;
    
    auto hDefTrack = (TH1*) h2DefTrack->ProjectionY("DefTrack", 1, h2DefTrack->GetNbinsX(), "e");
    if (!hDefTrack) {
        std::cerr << "Error: Failed to project default track histogram." << std::endl;
        return;
    }
    std::cout << "Successfully projected default track histogram." << std::endl;
    hDefTrack = hDefTrack->Rebin(nTrackptbin, "DefTrackRebin", Trackptbin);

    std::cout << "Accessing comparison track histogram..." << std::endl;
    auto h2CompTrack = (TH2 *)compFile->Get(Form("%s/%s", HistDirectory.Data(), TrackPtObj));
    if (!h2CompTrack) {
        std::cerr << "Error: Failed to retrieve comparison track histogram." << std::endl;
        return;
    }
    std::cout << "Successfully retrieved comparison track histogram." << std::endl;

    auto hCompTrack = (TH1*) h2CompTrack->ProjectionY("CompTrack", 1, h2CompTrack->GetNbinsX(), "e");
    if (!hCompTrack) {
        std::cerr << "Error: Failed to project comparison track histogram." << std::endl;
        return;
    }
    std::cout << "Successfully projected comparison track histogram." << std::endl;
    hCompTrack = hCompTrack->Rebin(nTrackptbin, "CompTrackRebin", Trackptbin);

    std::cout << "Drawing systematic comparison..." << std::endl;
    DrawSystematics(hDefTrack, hCompTrack, outputName.Data(), "#it{p}_{T} (GeV/c)", "dN/d#it{p}_{T}");

    std::cout << "Finished GetTrackErr." << std::endl;
}

void GetConstErr(TString defResult, TString compResult, TString outputName) {
    std::cout << "Opening default result file: " << defResult << std::endl;
    auto defFile = TFile::Open(defResult, "read");
    if (!defFile || defFile->IsZombie()) {
        std::cerr << "Error: Failed to open default result file." << std::endl;
        return;
    }

    std::cout << "Opening comparison result file: " << compResult << std::endl;
    auto compFile = TFile::Open(compResult, "read");
    if (!compFile || compFile->IsZombie()) {
        std::cerr << "Error: Failed to open comparison result file." << std::endl;
        return;
    }

    std::cout << "Accessing default track histogram..." << std::endl;
    auto h3DefConst = (TH3 *)defFile->Get(Form("%s/%s", HistDirectory.Data(), ConstPtObj));
    if (!h3DefConst) {
        std::cerr << "Error: Failed to retrieve default const histogram." << std::endl;
        return;
    }
    std::cout << "Successfully retrieved default const histogram." << std::endl;
    
    TH1 *hDefConst =  h3DefConst->ProjectionZ("DefConst",
                                   h3DefConst->GetXaxis()->FindBin(RBIN),
                                   h3DefConst->GetXaxis()->FindBin(RBIN),
                                   1, h3DefConst->GetNbinsY());
    if (!hDefConst) {
        std::cerr << "Error: Failed to project default const histogram." << std::endl;
        return;
    }
    std::cout << "Successfully projected default const histogram." << std::endl;
    hDefConst = hDefConst->Rebin(nTrackptbin, "DefConstRebin", Trackptbin);

    std::cout << "Accessing comparison const histogram..." << std::endl;
    auto h3CompConst = (TH3 *)compFile->Get(Form("%s/%s", HistDirectory.Data(), ConstPtObj));
    if (!h3CompConst) {
        std::cerr << "Error: Failed to retrieve comparison const histogram." << std::endl;
        return;
    }
    std::cout << "Successfully retrieved comparison const histogram." << std::endl;

    TH1 *hCompConst = h3CompConst->ProjectionZ("CompConst",
                                   h3CompConst->GetXaxis()->FindBin(RBIN),
                                   h3CompConst->GetXaxis()->FindBin(RBIN),
                                   1, h3CompConst->GetNbinsY());
    if (!hCompConst) {
        std::cerr << "Error: Failed to project comparison const histogram." << std::endl;
        return;
    }
    std::cout << "Successfully projected comparison const histogram." << std::endl;
    hCompConst = hCompConst->Rebin(nTrackptbin, "CompConstRebin", Trackptbin);

    std::cout << "Drawing systematic comparison..." << std::endl;
    DrawSystematics(hDefConst, hCompConst, outputName.Data(), "#it{p}_{T, constituent} (GeV/c)", "dN/d#it{p}_{T}");

    std::cout << "Finished GetConstErr." << std::endl;
}

void GetJetErr(TString defResult, TString compResult, TString outputName) {
    std::cout << "Opening default result file: " << defResult << std::endl;
    auto defFile = TFile::Open(defResult, "read");
    if (!defFile || defFile->IsZombie()) {
        std::cerr << "Error: Failed to open default result file." << std::endl;
        return;
    }

    std::cout << "Opening comparison result file: " << compResult << std::endl;
    auto compFile = TFile::Open(compResult, "read");
    if (!compFile || compFile->IsZombie()) {
        std::cerr << "Error: Failed to open comparison result file." << std::endl;
        return;
    }

    std::cout << "Accessing default track histogram..." << std::endl;
    auto hDefJet = (TH1 *)defFile->Get(Form("%s/%s", HistDirectory.Data(), JetPtObj));
    if (!hDefJet) {
        std::cerr << "Error: Failed to retrieve default Jet histogram." << std::endl;
        return;
    }
    std::cout << "Successfully retrieved default Jet histogram." << std::endl;
    hDefJet = hDefJet->Rebin(nptBins, "DefJetRebin", ptbin);

    std::cout << "Accessing comparison Jet histogram..." << std::endl;
    auto hCompJet = (TH1 *)compFile->Get(Form("%s/%s", HistDirectory.Data(), JetPtObj));
    if (!hCompJet) {
        std::cerr << "Error: Failed to retrieve comparison Jet histogram." << std::endl;
        return;
    }
    std::cout << "Successfully retrieved comparison Jet histogram." << std::endl;
    hCompJet = hCompJet->Rebin(nptBins, "CompJetRebin", ptbin);

    std::cout << "Drawing systematic comparison..." << std::endl;
    DrawSystematics(hDefJet, hCompJet, outputName.Data(), "#it{p}_{T, jet} (GeV/c)", "dN/d#it{p}_{T}");

    std::cout << "Finished GetJetErr." << std::endl;
}

void GetSystematicError(TString stdResult, TString systResult) {
    cout << "Opening standard result file: " << stdResult << endl;
    auto stdFile = TFile::Open(stdResult, "read");
    if (!stdFile || stdFile->IsZombie()) {
        cout << "Error: Failed to open standard result file." << endl;
        return;
    }

    TIter nextStd(stdFile->GetListOfKeys()); 
    TKey *keyStd; 
    TH1 *hstd = nullptr; 
    while ((keyStd = (TKey*)nextStd())) {
        TString name = keyStd->GetName();
        cout << "Checking object: " << name << endl;
        if (name.Contains("Run3_CrossSection")) {
            hstd = (TH1*)stdFile->Get(name);
            cout << "Found standard histogram: " << name << endl;
            break;
        }
    }
    if (!hstd) {
        cout << "Error: Run3_CrossSection histogram not found in standard result file." << endl;
    }
    // stdFile->Close();

    cout << "Opening systematic result file: " << systResult << endl;
    auto systFile = TFile::Open(systResult, "read");
    if (!systFile || systFile->IsZombie()) {
        cout << "Error: Failed to open systematic result file." << endl;
        return;
    }

    TIter nextSyst(systFile->GetListOfKeys());
    TKey *keySyst;
    TH1 *hsyst = nullptr; 
    while ((keySyst = (TKey*)nextSyst())) {
        TString name = keySyst->GetName();
        cout << "Checking object: " << name << endl;
        if (name.Contains("Run3_CrossSection")) {
            hsyst = (TH1*)systFile->Get(name);
            cout << "Found systematic histogram: " << name << endl;
            break;
        }
    }
    if (!hsyst) {
        cout << "Error: Run3_CrossSection histogram not found in systematic result file." << endl;
    }
    // systFile->Close();

    if (hstd && hsyst) {
        cout << "Drawing histograms..." << endl;
        DrawSystematics(hstd, hsyst, "TrackingEfficiency", "#it{p}_{T, jet} (GeV/c)", JetPtTitleY);
    } else {
        cout << "Error: Could not find histograms in one or both files." << endl;
    }
}

TH1 *DrawSystematics(TH1 *hstd, TH1 *hsyst, TString FileName, const char *XaxisTitle, const char *YaxisTitle) {
    Filipad2 *pad = new Filipad2(++nn, 2, 0.45, 100, 50, 0.7, 1, 1);
    pad->Draw();
    TPad *upperpad = pad->GetPad(1);
    TPad *lowerpad = pad->GetPad(2);
    TLegend *leg = new TLegend(0.449761,0.711801,0.715311,0.95528,NULL,"brNDC");
        leg->SetTextSize(0.05);
        leg->SetBorderSize(0);
        leg->SetFillColorAlpha(0,0);
    upperpad->cd();
    optFili(*upperpad, 1, 1, 0, 1);

    hstd->SetMarkerStyle(20);
    hstd->SetMarkerColor(kBlack);
    hstd->SetLineColor(kBlack);
    hstd->SetMarkerSize(1.3);
    hset(*hstd, XaxisTitle, YaxisTitle, 1.15, 1.1, 0.065, 0.065, 0.01, 0.01,
       0.08, 0.06, 510, 505);
    hstd->Draw("pe");
    hsyst->SetMarkerStyle(24);
    hsyst->SetMarkerColor(kRed);
    hsyst->SetLineColor(kRed);
    hsyst->SetMarkerSize(1.3);
    hsyst->Draw("pesame");

    leg->AddEntry(hstd, "Default", "pe");
    leg->AddEntry(hsyst, "Tracking efficiency -3%", "pe");
    leg->Draw();

    lowerpad->cd();
    optFili(*lowerpad, 1, 1, 0, 0);
    auto ratio = (TH1 *) hsyst->Clone();
    ratio->Divide(ratio, hstd, 1, 1, "B");
    hset(*ratio, XaxisTitle, "variation / default", 1.15, 1.05, 0.075, 0.075, 0.01, 0.01,
       0.08, 0.08, 510, 505);
    hoptset(*ratio, 0, kRed, 5, 300, 0.99, 1.11, 1, 1, 1, 24);
    ratio->Draw("pe");

    // fitting
    TF1 *fitFunc = new TF1("fitFunc", "pol0", 40, PlotPtMax);
    ratio->Fit(fitFunc, "", "", 40, PlotPtMax);
    fitFunc->Draw("same");

    std::ofstream outFile(Form("../plots/systematic/%s/SystErrTrmkEff_%s.txt", histNames[0].Data(), FileName.Data()));
    if (!outFile) {
        std::cerr << "Error: Could not open the file!" << std::endl;
        return nullptr;
    }

    for (int i = 1; i <= ratio->GetNbinsX(); i++) {
        double binCenter = ratio->GetXaxis()->GetBinCenter(i);
        double binContent = (ratio->GetBinContent(i) - 1) * 100;
        TString outputLine = Form("Bin: %d, pT: %.2f, Error: %.5f %%\n", i, binCenter, binContent);
        
        std::cout << outputLine;
        outFile << outputLine;
    }
    double avgContent = ((1 - hsyst->GetEntries()/hstd->GetEntries()) * 100);
    TString summarizeLline = Form("Total Error: %.5f %%\n", avgContent);
    std::cout << summarizeLline;
    outFile << summarizeLline;

    TFile* outputFile = TFile::Open(Form("Syst%s.root", FileName.Data()), "RECREATE");
    ratio->Write();
    outputFile->Close();

    outFile.close();

    pad->C->SaveAs(Form("%s/systematics/SystErr%s.pdf", MakeDirName.Data(), FileName.Data()));
}

void DrawJetsMCsyst_trackingEfficiency() {
    GetSystematicError("../plots/systematic/UnfoldedData_standard.root", "../plots/systematic/UnfoldedData_trackingefficiency.root");
    // GetTrackErr("../../../jets/mc/AnalysisResults/LHC24f3/selMC/AnalysisResults.root", "../../../jets/mc/AnalysisResults/LHC24f3/selMC_syst/trackingEfficiency/AnalysisResults.root", "TrackingEfficiency_Tracks");
    // GetConstErr("../../../jets/mc/AnalysisResults/LHC24f3/selMC/AnalysisResults.root", "../../../jets/mc/AnalysisResults/LHC24f3/selMC_syst/trackingEfficiency/AnalysisResults.root", "TrackingEfficiency_Constituents");
    // GetJetErr("../../../jets/mc/AnalysisResults/LHC24f3/selMC/AnalysisResults.root", "../../../jets/mc/AnalysisResults/LHC24f3/selMC_syst/trackingEfficiency/AnalysisResults.root", "TrackingEfficiency_Jets");
    // GetTrackErr("../../../jets/mc/jetfinderQA/AnalysisResults/LHC24f3/selMC_globalTracks/systematics/AnalysisResults_default.root", "../../../jets/mc/jetfinderQA/AnalysisResults/LHC24f3/selMC_globalTracks/systematics/AnalysisResults_randomUniform.root", "TrackingEfficiency_TracksUniform");
    // GetConstErr("../../../jets/mc/jetfinderQA/AnalysisResults/LHC24f3/selMC_globalTracks/systematics/AnalysisResults_default.root", "../../../jets/mc/jetfinderQA/AnalysisResults/LHC24f3/selMC_globalTracks/systematics/AnalysisResults_randomUniform.root", "TrackingEfficiency_ConstituentsUniform");
}