#include "Filipad2.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TString.h"
#include <iostream>

// Simple optFili function
void optFili(TPad &pid, Int_t gridx, Int_t gridy, Int_t logx, Int_t logy) {
  pid.SetGridy(gridx);
  pid.SetGridx(gridy);
  pid.SetLogx(logx);
  pid.SetLogy(logy);
}

// Simple hset function
void hset(TH1D &hid, TString xtit = "", TString ytit = "", double titoffx = 0.9,
          double titoffy = 1.0, double titsizex = 0.07, double titsizey = 0.07,
          double labeloffx = 0.01, double labeloffy = 0.01,
          double labelsizex = 0.07, double labelsizey = 0.07,
          int divx = 510, int divy = 505) {
  hid.GetXaxis()->CenterTitle(1);
  hid.GetYaxis()->CenterTitle(1);
  hid.GetXaxis()->SetTitleOffset(titoffx);
  hid.GetYaxis()->SetTitleOffset(titoffy);
  hid.GetXaxis()->SetTitleFont(43);
  hid.GetYaxis()->SetTitleFont(43);
  hid.GetXaxis()->SetTitleSize(titsizex);
  hid.GetYaxis()->SetTitleSize(titsizey);
  hid.GetXaxis()->SetLabelOffset(labeloffx);
  hid.GetYaxis()->SetLabelOffset(labeloffy);
  hid.GetXaxis()->SetLabelFont(43);
  hid.GetYaxis()->SetLabelFont(43);
  hid.GetXaxis()->SetLabelSize(labelsizex);
  hid.GetYaxis()->SetLabelSize(labelsizey);
  hid.GetXaxis()->SetNdivisions(divx);
  hid.GetYaxis()->SetNdivisions(divy);
  hid.GetXaxis()->SetTitle(xtit);
  hid.GetYaxis()->SetTitle(ytit);
}

const Double_t ptbinGen[26] = {0, 1, 2, 3, 4, 5,  6,  7,  8,  9,  10, 12, 14, 16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
const Int_t nptBinsGen = sizeof(ptbinGen) / sizeof(ptbinGen[0]) - 1;

///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Draw macro for PYTHIA comparison //////
////////// 13.6TeV vs 13TeV pp charged jets //////
////////// author: Assistant                 //////
////////// Last Modified: 2025               //////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

void DrawPYTHIAComparison() {
    
    // File paths
    TString file13TeV = "~/downloads/13000GeV_nevent_10_9.root";
    TString file136TeV = "~/downloads/13600GeV_nevent_10_9.root";
    
    // Open files
    TFile* f13 = TFile::Open(file13TeV.Data());
    TFile* f136 = TFile::Open(file136TeV.Data());
    
    if (!f13 || !f136) {
        cout << "Error: Cannot open files" << endl;
        return;
    }
    
    // Get histograms
    TH1D* hJetPt13 = (TH1D*)f13->Get("hJetPt");
    TH1D* hJetPt136 = (TH1D*)f136->Get("hJetPt");
    TH1D* hSigmaGen13 = (TH1D*)f13->Get("hSigmaGen");
    TH1D* hSigmaGen136 = (TH1D*)f136->Get("hSigmaGen");
    
    if (!hJetPt13 || !hJetPt136 || !hSigmaGen13 || !hSigmaGen136) {
        cout << "Error: Cannot find required histograms" << endl;
        f13->Close();
        f136->Close();
        return;
    }

    // Rebin histograms using the simple method
    hJetPt13 = (TH1D*)hJetPt13->Rebin(nptBinsGen, "hJetPt13Rebin", ptbinGen);
    hJetPt136 = (TH1D*)hJetPt136->Rebin(nptBinsGen, "hJetPt136Rebin", ptbinGen);
    
    // Get cross sections from first bin of hSigmaGen
    Double_t sigma13 = hSigmaGen13->GetBinContent(1) / 1000;
    Double_t sigma136 = hSigmaGen136->GetBinContent(1) / 1000;
    
    cout << "13TeV cross section: " << sigma13 << " mb" << endl;
    cout << "13.6TeV cross section: " << sigma136 << " mb" << endl;
    
    // Clone rebinned histograms for normalization
    TH1D* hJetPt13Norm = (TH1D*)hJetPt13->Clone("hJetPt13Norm");
    TH1D* hJetPt136Norm = (TH1D*)hJetPt136->Clone("hJetPt136Norm");
    
    // Normalize by cross section with width option
    hJetPt13Norm->Scale(sigma13/75.4, "width");
    hJetPt136Norm->Scale(sigma136/78.6, "width");
    
    // Create Filipad2 for plotting
    Filipad2 *JetPtPad = new Filipad2(1, 2, 0.4, 100, 50, 0.7, 1, 1);
    JetPtPad->Draw();
    
    TPad *jetptpad = JetPtPad->GetPad(1);
    optFili(*jetptpad, 1, 1, 0, 1);
    TPad *ratiojetptpad = JetPtPad->GetPad(2);
    optFili(*ratiojetptpad, 1, 1, 0, 0);
    
    // Create legend
    TLegend *legjetpt = new TLegend(0.545455, 0.715942, 0.815789, 0.95942, NULL, "brNDC");
    legjetpt->SetTextSize(0.05);
    legjetpt->SetBorderSize(0);
    
    // Draw main histograms
    jetptpad->cd();
    
    // Set histogram properties
    hJetPt13Norm->SetLineColor(kBlue);
    hJetPt13Norm->SetMarkerColor(kBlue);
    hJetPt13Norm->SetMarkerStyle(20);
    hJetPt13Norm->SetMarkerSize(0.6);
    
    hJetPt136Norm->SetLineColor(kRed);
    hJetPt136Norm->SetMarkerColor(kRed);
    hJetPt136Norm->SetMarkerStyle(21);
    hJetPt136Norm->SetMarkerSize(0.6);
    
    // Set axis titles and properties BEFORE drawing
    hJetPt13Norm->GetXaxis()->SetTitle("#it{p}_{T, jet}^{gen} (GeV/#it{c})");
    hJetPt13Norm->GetYaxis()->SetTitle("d^{2}#sigma/d#it{#eta}d#it{p}_{T, jet} (mb GeV/#it{c})^{-1}");
    hJetPt13Norm->GetXaxis()->SetTitleSize(0.07);
    hJetPt13Norm->GetYaxis()->SetTitleSize(0.07);
    hJetPt13Norm->GetXaxis()->SetLabelSize(0.07);
    hJetPt13Norm->GetYaxis()->SetLabelSize(0.07);
    hJetPt13Norm->GetXaxis()->SetTitleOffset(1.2);
    hJetPt13Norm->GetYaxis()->SetTitleOffset(1.0);
    
    // Draw histograms
    hJetPt13Norm->Draw("EP");
    hJetPt136Norm->Draw("EP same");
    
    // Add to legend
    legjetpt->AddEntry(hJetPt13Norm, "PYTIHA 13 TeV", "lep");
    legjetpt->AddEntry(hJetPt136Norm, "PYTHIA 13.6 TeV", "lep");
    legjetpt->Draw();
    
    // Create ratio histogram
    TH1D* hRatio = (TH1D*)hJetPt13Norm->Clone("hRatio");
    hRatio->Divide(hJetPt136Norm);
    
    // Draw ratio
    ratiojetptpad->cd();
    hRatio->SetLineColor(kBlue);  // Same color as 13TeV (numerator)
    hRatio->SetMarkerColor(kBlue);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.8);
    
    hRatio->GetXaxis()->SetTitle("#it{p}_{T, jet}^{gen} (GeV/#it{c})");
    hRatio->GetYaxis()->SetTitle("13 TeV / 13.6 TeV");
    hRatio->GetXaxis()->SetTitleSize(0.07);
    hRatio->GetYaxis()->SetTitleSize(0.07);
    hRatio->GetXaxis()->SetLabelSize(0.07);
    hRatio->GetYaxis()->SetLabelSize(0.07);
    hRatio->GetXaxis()->SetTitleOffset(1.2);
    hRatio->GetYaxis()->SetTitleOffset(1.0);
    hRatio->GetYaxis()->SetRangeUser(0.88, 1.02);
    
    // Draw horizontal line at 1
    TLine* line1 = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
    line1->SetLineStyle(2);
    line1->SetLineColor(kGray);
    
    hRatio->Draw("EP");
    line1->Draw("same");
    
    // Print plot
    JetPtPad->C->Print("PYTHIA_Comparison_13TeV_vs_136TeV.pdf");
    
    // // Clean up
    // f13->Close();
    // f136->Close();
    
    cout << "Plot saved as PYTHIA_Comparison_13TeV_vs_136TeV.pdf" << endl;
}
