#include <__config>
#include <iostream>
#include <vector>
#include "BSHelper.cxx"
#include "Filipad2.h"
using namespace std;



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
{		// hid.SetStats(0);

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

void optFili(TPad& pid, Int_t gridx, Int_t gridy, Int_t logx, Int_t logy) {
    pid.SetGridy(gridx);
    pid.SetGridx(gridy);
    pid.SetLogx(logx);
    pid.SetLogy(logy);
}

void DrawCollRecEff() {
    double Trackptbin[15] = {0, 2, 4, 6, 8, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80};
    int nTrackptbin = sizeof(Trackptbin) / sizeof(Trackptbin[0]) - 1;

    auto file = TFile::Open("CollRecEff.root", "open");
    auto Triggered = (TH1 *) file->Get("hMCPMatched_Jet-MC (local)");
    auto UnTriggered = (TH1 *) file->Get("hMCPMatched_Jet-MC (N=2)");
    
    auto hTimeFrameEff = (TH1 *) UnTriggered->Clone();
    hTimeFrameEff->Divide(hTimeFrameEff, Triggered, 1, 1, "B");
  	hTimeFrameEff = hTimeFrameEff->Rebin(nTrackptbin, "hTimeFrameEff", Trackptbin); 

    TCanvas* canTFE = new TCanvas("canTFE", "Canvas for Time Frame Efficiency", 800, 800);
    canTFE->SetLeftMargin(0.15);
    canTFE->cd();
    optFili(*canTFE, 1, 1, 0, 0);

    hTimeFrameEff->SetTitle("");
    hset(*hTimeFrameEff, "#it{p}_{T, jet}^{Gen} (GeV/c)", "Time frame efficiency = #varepsilon_{jet}^{reco, untrigger}/#varepsilon_{jet}^{reco, trigger}", 1.1, 1.5, 0.04, 0.045, 0.01, 0.0001, 0.035, 0.035, 510, 510);
    hTimeFrameEff->GetYaxis()->SetRangeUser(0.6, 1.0);
    hTimeFrameEff->Draw("pe");

}