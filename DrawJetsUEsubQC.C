///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Draw macro for jet QA         //////////
////////// author: Joonsuk Bae           ////////// 
////////// E-mail: jbae@cern.ch          //////////
////////// Last Modified: 08 Nov 2024    //////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

#include "DrawJetsCommon.h"

TString JetPtTitleX = "#it{p}_{T, jet}^{reco} (GeV/#it{c})";
TString JetPtTitleY = "1/#it{N}_{jet} d#it{N}/d#it{p}_{T}";

const char *DataDirectory = "jet-finder-charged-qa";
const char *EventObj = "h_collisions";
const char *hJetPtWoUESubObj = "h_jet_pt";
const char *hJetPtWUESubObj = "h_jet_pt_rhoareasubtracted";
const char *h2JetPtRhoObj = "h2_centrality_rho";

void DrawPtComparison (const std::vector<TString> &fileNames, const std::vector<TString> &histNames);

void DrawJetsUEsubQC() {
    fileNames = {"~/cernbox/workspace/O2Physics/jets/data/AnalysisResults/LHC22o_apass7_minBias_small/sel8/UEsub_TrackTuner_Const100GeV/AreaBasedDoSparse/CombinedResults.root"};
    histNames = {"area-based-Sparse"};

    cout << "Debugger " << ++n << endl;

    DrawPtComparison(fileNames, histNames);
}

void DrawPtComparison(const std::vector<TString> &fileNames, const std::vector<TString> &histNames) {
    cout << "Debugger " << ++n << endl;
    Filipad2* ptcanvas = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    ptcanvas->Draw();
    TPad *ptpad = ptcanvas->GetPad(1);
    optFili(*ptpad, 1, 1, 0, 1);
    TPad *ratioptpad = ptcanvas->GetPad(2);
    optFili(*ratioptpad, 1, 1, 0, 0);
    cout << "Debugger " << ++n << endl;

    TH1 *hJetPtWoUESub = nullptr;

    for (Int_t i = 0; i < fileNames.size(); i++) {  // 수정: i < fileNames.size()
        TFile *filename = TFile::Open(fileNames[i].Data(), "open");
        if (!filename || filename->IsZombie()) {
            std::cerr << "Error: Could not open file " << fileNames[i].Data() << std::endl;
            continue;
        }
    cout << "Debugger " << ++n << endl;

        TH1 *Nevt = (TH1*) filename->Get(Form("%s/%s", DataDirectory, EventObj));
        if (!Nevt) {
            std::cerr << "Error: Could not find histogram " << EventObj << " in file " << fileNames[i].Data() << std::endl;
            filename->Close();
            continue;
        }
        Double_t nevt = Nevt->GetBinContent(Nevt->FindBin(1.5));
    cout << "Debugger " << ++n << endl;

        if (i==0) {
            hJetPtWoUESub = (TH1*) filename->Get(Form("%s/%s", DataDirectory, hJetPtWoUESubObj));
        }
        TH1* hJetPtWUESub = (TH1*) filename->Get(Form("%s/%s", DataDirectory, hJetPtWUESubObj));
        TH2 *h2JetPtRho = (TH2*) filename->Get(Form("%s/%s", DataDirectory, h2JetPtRhoObj));
            TH1 *hJetPtRho = (TH1*) h2JetPtRho->ProjectionY("hJetPtRho", 1, h2JetPtRho->GetNbinsX());
        
        if (!hJetPtWoUESub || !hJetPtWUESub || !hJetPtRho) {
            std::cerr << "Error: Could not find required histograms in file " << fileNames[i].Data() << std::endl;
            filename->Close();
            continue;
        }

        ptpad->cd();
        if (i==0) {
        hoptset(*hJetPtWoUESub, nevt, ColorPallete[3 * i], -10, 200, 1e-10, 1e0, 0.75, 1, 1, 20);
        hJetPtWoUESub->Draw("PESAME");
        }
        hoptset(*hJetPtWUESub, nevt, ColorPallete[3 * i + 1], -10, 200, 1e-10, 1e0, 0.75, 1, 1, 20);
        hJetPtWUESub->Draw("PESAME");
        hoptset(*hJetPtRho, 2.3312e9, ColorPallete[3 * i + 2], -10, 200, 1e-10, 1e0, 0.75, 1, 1, 20);
        hJetPtRho->Draw("PESAME");

        TLegend *legend = new TLegend(0.6, 0.6, 0.85, 0.85);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->AddEntry(hJetPtWoUESub, "w/o UE sub", "p");
        legend->AddEntry(hJetPtWUESub, "w/ UE sub", "p"); 
        legend->AddEntry(hJetPtRho, "#rho", "p");
        legend->Draw();

        ratioptpad->cd();
        TH1 *JetPtRat = DrawRatioTH1(hJetPtWUESub, hJetPtWoUESub);
        hoptset(*JetPtRat, 0, ColorPallete[3 * i], -10, 200, 0, 1, 0.75, 1, 1, 20);
        JetPtRat->Draw("PESAME");

        // filename->Close();  // 파일 닫기
    }
}