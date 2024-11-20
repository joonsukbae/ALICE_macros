#include "DrawJetsMCFilesTitles.h"

void DrawJetPtRecComparison() {
    auto JetPtTuned = TFile::Open("JetPtMCDafTrkTuned.root", "open");
    auto jetpttuned = (TH1 *) JetPtTuned->Get("TrackPt_MB-pass7-sampling_sel8_selMC_selectedWindow");

    auto JetPtunTuned = TFile::Open("JetPtMCDbfTrkTuned.root", "open");
    auto jetptuntuned = (TH1 *) JetPtunTuned->Get("TrackPt_MB-pass7-sampling_sel8_selMC_selectedWindow");

    Filipad2 *TrackPtPad = new Filipad2(0, 2, 0.4, 100, 50, 0.7, 1, 1);
    TrackPtPad->Draw();
    TPad *trackptpad = TrackPtPad->GetPad(1);
    optFili(*trackptpad, 1, 1, 0, 1);
    trackptpad->cd();

    jetpttuned->SetMarkerSize(1.5);
    jetpttuned->SetMarkerColor(kRed);
    jetpttuned->SetLineColor(kRed);
    jetpttuned->Draw("pe");
    jetptuntuned->Draw("pesame");

    TPad *trackptpadrat = TrackPtPad->GetPad(2);
    trackptpadrat->cd();
    optFili(*trackptpadrat, 1, 1, 0, 0);

    auto ratio = (TH1 *) jetpttuned->Clone();
    ratio->Divide(ratio, jetptuntuned, 1, 1, "B");
    // ratio->SetMarkerColor(kRed);
    // ratio->SetLineColor(kRed);
    ratio->GetYaxis()->SetTitle("Tuned / pure");
    ratio->Draw("pe");


}