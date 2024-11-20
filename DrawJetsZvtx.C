///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Draw macro for jet QA         //////////
////////// author: Joonsuk Bae           ////////// 
////////// E-mail: jbae@cern.ch          //////////
////////// Last Modified: 19 Aug 2024    //////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////
#include "DrawJetsMCFilesTitles.h"
// #include "DrawJetsMCFunctions.h"

const char *objectZvtx = "h_collisions_zvtx";

TH1 *GetZvtx(const char *FileName, const char *objectName);
// void DrawComparison(TH1 *hDef, TH1 *hComp);

void DrawJetsZvtx() {
    TH1 *hzvtx = GetZvtx("/Users/js/cernbox/workspace/O2Physics/jets/mc/jetfinderQA/AnalysisResults_MCzvtx.root", objectZvtx);
    // TH1 *hmczvtx = GetZvtx("", objectZvtx);
    // DrawComparison(hzvtx, hmczvtx);

    double Ntrig = hzvtx->GetEntries();
    double Nzvtx10 = 0;
    int lBinZvtx = hzvtx->FindBin(-10); 
    int rBinZvtx = hzvtx->FindBin(10); 
    for (int i=lBinZvtx; i<=rBinZvtx; i++) {
        double Nevt = hzvtx->GetBinContent(i);
        Nzvtx10 += Nevt;
    }
    double effTrigZvtx10 = Nzvtx10 / Ntrig;
    cout << "(Ntrig, Nzvtx10, effTrigZvtx10): (" << Ntrig << ", " << Nzvtx10 << ", " << effTrigZvtx10 << ")" << endl;

    auto can1 = new TCanvas("c1", "c1", 1000, 800);
    setpad(can1);
    auto leg1 = new TLegend(0.151529,0.701613,0.40594,0.950484,NULL,"brNDC");
        leg1->SetTextSize(0.05);
        leg1->SetBorderSize(0);
        leg1->SetFillColorAlpha(0,0); 
    can1->cd();
    optFili(*can1, 1, 1, 0, 0);
    hset(*hzvtx, "#it{z}_{vtx} (cm)", "Normalized entries", 1, 1.3, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hoptset(*hzvtx, 1, kBlack, -30, 30, 0, 1e-1);
    hzvtx->Draw("pe");

    auto Fit = new TF1("fit", "gaus", -10, 10);
    hzvtx->Fit("fit");
    Fit->SetLineColor(kRed);
    Fit->Draw("l same");

    double Probability = Fit->Integral(-10, 10);
    cout << "Probability (-10 < zvtx < 10 cm): " << Probability << endl;
    leg1->AddEntry(Fit, Form("#int_{-10 cm}^{10 cm} z_{vtx} = %.3f", Probability), "l");
    leg1->Draw();

    can1->SaveAs(Form("%s/Zvtx.pdf", MakeDirName.Data()));

}

TH1 *GetZvtx(const char *FileName, const char *objectName) {
    auto File = TFile::Open(FileName, "read");
    auto hist = (TH1 *) File->Get(Form("jet-finder-charged-qa/%s", objectName));

    return hist;
}

// void DrawComparison(TH1 *hDef, TH1 *hComp) {
//     Filipad2 *Pad = new Filipad2(++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
//       Pad->Draw();
//       TPad *upPad = Pad->GetPad(1);
//       optFili(*upPad, 1, 1, 0, 1);
//       TPad *belowPad = Pad->GetPad(2);
//       optFili(*belowPad, 1, 1, 0, 0);
//   TLegend *leg1 =
//           new TLegend(0.313397,0.6,0.578947,0.843478,NULL,"brNDC");
//           leg1->SetTextSize(0.05);
//           leg1->SetBorderSize(0);
//           leg1->SetFillColorAlpha(0,0); 
  
//   upPad->cd();
//   leg1->AddEntry(hDef, "raw Data");
//   hoptset(*hDef, 0, kRed, 5, PlotPtMax, 3e-10, 1e-1, 1.4);
//   hDef->Draw("pe");

//   leg1->AddEntry(hComp, "reco MC");
//   hset(*hComp, JetPtTitleX, JetPtDataFinalTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
//        0.05, 0.05, 510, 505);
//   hoptset(*hComp, 0, kBlack, 5, PlotPtMax, 3e-10, 1e-1);
//   hComp->Draw("pesame");

//   leg1->Draw();

//   belowPad->cd();
//   auto hRatio = DrawRatioTH1(hDef, hComp); 
//   hset(*hRatio, JetPtGenTitleX, "Herwig / PYTHIA8", 1.2, 0.7, 0.1, 0.1, 0.01, 0.01,
//        0.1, 0.1, 510, 505);
//   hoptset(*hRatio, 0, ColorPallete[3], 5, PlotPtMax, 0.5, 1.5);
//   hRatio->Draw("pe");

// //   OutStatsTXT(hRatio, "UnfoldRatPriorHerwig");
// //   Pad->C->SaveAs("plots/AN_Charged-particle-jet-cross-section-in-pp-collisions-at-13.6-TeV/Figures/systematics/SystErrUnfoldingHerwig.pdf");
// }