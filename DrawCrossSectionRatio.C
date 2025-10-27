///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Draw macro for Cross Section Ratio  //////////
////////// author: Joonsuk Bae                 ////////// 
////////// E-mail: jbae@cern.ch                //////////
////////// Last Modified: 08 Nov 2024           //////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

#include "DrawJetsCommon.h"

void DrawCrossSectionRatio() {
    // 파일 경로 설정
    TString trackingEffFile = "Run3_CrossSection_TrackingEfficiency.root";
    TString referenceFile = "Run3_CrossSection.root";
    TString woTrackTunerFile = "Run3_CrossSection_WoTrackTuner.root";
    
    // histogram 이름 설정 (파일에서 확인된 이름 사용)
    TString histName = "Run3_CrossSection";
    
    // 출력 디렉토리 생성
    system("mkdir -p plots/ratio");
    
    // 파일 열기
    TFile* fTrackingEff = TFile::Open(trackingEffFile.Data(), "READ");
    TFile* fReference = TFile::Open(referenceFile.Data(), "READ");
    TFile* fWoTrackTuner = TFile::Open(woTrackTunerFile.Data(), "READ");
    
    if (!fTrackingEff || !fReference || !fWoTrackTuner) {
        std::cerr << "Error: Could not open one or more files!" << std::endl;
        return;
    }
    
    // histogram 가져오기
    TH1F* hTrackingEff = (TH1F*)fTrackingEff->Get(histName.Data());
    TH1F* hReference = (TH1F*)fReference->Get(histName.Data());
    TH1F* hWoTrackTuner = (TH1F*)fWoTrackTuner->Get(histName.Data());
    
    if (!hTrackingEff || !hReference || !hWoTrackTuner) {
        std::cerr << "Error: Could not find histogram '" << histName << "' in one or more files!" << std::endl;
        fTrackingEff->Close();
        fReference->Close();
        fWoTrackTuner->Close();
        return;
    }
    
    // Canvas 생성 (두 개의 pad로 나누기)
    TCanvas* c1 = new TCanvas("c1", "Cross Section Ratios", 1200, 600);
    c1->Divide(2, 1);
    
    // 첫 번째 ratio: Tracking Efficiency / Reference
    c1->cd(1);
    TPad* pad1 = (TPad*)c1->GetPad(1);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.05);
    pad1->SetTopMargin(0.05);
    pad1->SetBottomMargin(0.12);
    
    // Ratio histogram 생성 (Tracking Efficiency / Reference)
    TH1F* hRatio1 = (TH1F*)hTrackingEff->Clone("hRatio1");
    hRatio1->Reset();
    
    // Ratio 계산 (Tracking Efficiency / Reference)
    for (int i = 1; i <= hTrackingEff->GetNbinsX(); i++) {
        double xCenter = hTrackingEff->GetXaxis()->GetBinCenter(i);
        double yNum = hTrackingEff->GetBinContent(i);
        double yNumErr = hTrackingEff->GetBinError(i);
        
        int binDenom = hReference->GetXaxis()->FindBin(xCenter);
        double yDenom = hReference->GetBinContent(binDenom);
        double yDenomErr = hReference->GetBinError(binDenom);
        
        if (yDenom != 0) {
            double ratio = yNum / yDenom;
            // Error propagation: δ(ratio) = ratio * sqrt((δnum/num)² + (δdenom/denom)²)
            double ratioErr = ratio * sqrt(pow(yNumErr/yNum, 2) + pow(yDenomErr/yDenom, 2));
            
            hRatio1->SetBinContent(i, ratio);
            hRatio1->SetBinError(i, ratioErr);
        } else {
            hRatio1->SetBinContent(i, 0);
            hRatio1->SetBinError(i, 0);
        }
    }
    
    // Histogram 스타일링
    hset(*hRatio1, "Jet p_{T} (GeV/c)", "Ratio (Tracking Efficiency / Reference)", 
         1.2, 1.0, 0.07, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 505);
    
    hRatio1->SetMarkerColor(kRed);
    hRatio1->SetLineColor(kRed);
    hRatio1->SetMarkerSize(1.0);
    hRatio1->SetMarkerStyle(20);
    hRatio1->SetLineWidth(2);
    
    // Y축 범위 설정
    hRatio1->GetYaxis()->SetRangeUser(0.5, 1.5);
    
    // Grid 설정
    pad1->SetGridx();
    pad1->SetGridy();
    
    // Histogram 그리기
    hRatio1->Draw("E");
    
    // Reference line (y=1) 그리기
    TLine* refLine1 = new TLine(hRatio1->GetXaxis()->GetXmin(), 1.0, 
                                hRatio1->GetXaxis()->GetXmax(), 1.0);
    refLine1->SetLineColor(kBlack);
    refLine1->SetLineStyle(2);
    refLine1->SetLineWidth(2);
    refLine1->Draw("same");
    
    // Legend 생성
    TLegend* legend1 = new TLegend(0.65, 0.75, 0.85, 0.85);
    legend1->SetBorderSize(0);
    legend1->SetFillStyle(0);
    legend1->SetTextSize(0.04);
    legend1->AddEntry(hRatio1, "Tracking Efficiency / Reference", "p");
    legend1->AddEntry(refLine1, "Reference (y=1)", "l");
    legend1->Draw();
    
    // ALICE label 추가
    ALICEfigureLegend("ALICE Preliminary", 0.15, 0.85, 0.35, 0.95, 0.0, 0.0, 0.0, 0.0, 0.04);
    
    // 두 번째 ratio: WoTrackTuner / Reference
    c1->cd(2);
    TPad* pad2 = (TPad*)c1->GetPad(2);
    pad2->SetLeftMargin(0.12);
    pad2->SetRightMargin(0.05);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.12);
    
    // Ratio histogram 생성 (WoTrackTuner / Reference)
    TH1F* hRatio2 = (TH1F*)hWoTrackTuner->Clone("hRatio2");
    hRatio2->Reset();
    
    // Ratio 계산 (WoTrackTuner / Reference)
    for (int i = 1; i <= hWoTrackTuner->GetNbinsX(); i++) {
        double xCenter = hWoTrackTuner->GetXaxis()->GetBinCenter(i);
        double yNum = hWoTrackTuner->GetBinContent(i);
        double yNumErr = hWoTrackTuner->GetBinError(i);
        
        int binDenom = hReference->GetXaxis()->FindBin(xCenter);
        double yDenom = hReference->GetBinContent(binDenom);
        double yDenomErr = hReference->GetBinError(binDenom);
        
        if (yDenom != 0) {
            double ratio = yNum / yDenom;
            // Error propagation: δ(ratio) = ratio * sqrt((δnum/num)² + (δdenom/denom)²)
            double ratioErr = ratio * sqrt(pow(yNumErr/yNum, 2) + pow(yDenomErr/yDenom, 2));
            
            hRatio2->SetBinContent(i, ratio);
            hRatio2->SetBinError(i, ratioErr);
        } else {
            hRatio2->SetBinContent(i, 0);
            hRatio2->SetBinError(i, 0);
        }
    }
    
    // Histogram 스타일링
    hset(*hRatio2, "Jet p_{T} (GeV/c)", "Ratio (WoTrackTuner / Reference)", 
         1.2, 1.0, 0.07, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 505);
    
    hRatio2->SetMarkerColor(kBlue);
    hRatio2->SetLineColor(kBlue);
    hRatio2->SetMarkerSize(1.0);
    hRatio2->SetMarkerStyle(20);
    hRatio2->SetLineWidth(2);
    
    // Y축 범위 설정
    hRatio2->GetYaxis()->SetRangeUser(0.5, 1.5);
    
    // Grid 설정
    pad2->SetGridx();
    pad2->SetGridy();
    
    // Histogram 그리기
    hRatio2->Draw("E");
    
    // Reference line (y=1) 그리기
    TLine* refLine2 = new TLine(hRatio2->GetXaxis()->GetXmin(), 1.0, 
                                hRatio2->GetXaxis()->GetXmax(), 1.0);
    refLine2->SetLineColor(kBlack);
    refLine2->SetLineStyle(2);
    refLine2->SetLineWidth(2);
    refLine2->Draw("same");
    
    // Legend 생성
    TLegend* legend2 = new TLegend(0.65, 0.75, 0.85, 0.85);
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->SetTextSize(0.04);
    legend2->AddEntry(hRatio2, "WoTrackTuner / Reference", "p");
    legend2->AddEntry(refLine2, "Reference (y=1)", "l");
    legend2->Draw();
    
    // ALICE label 추가
    ALICEfigureLegend("ALICE Preliminary", 0.15, 0.85, 0.35, 0.95, 0.0, 0.0, 0.0, 0.0, 0.04);
    
    // Plot 저장
    c1->Print("plots/ratio/CrossSectionRatio.pdf");
    c1->Print("plots/ratio/CrossSectionRatio.png");
    
    // 통계 정보 출력
    std::cout << "=== Cross Section Ratio Analysis ===" << std::endl;
    std::cout << "Tracking Efficiency file: " << trackingEffFile.Data() << std::endl;
    std::cout << "Reference file: " << referenceFile.Data() << std::endl;
    std::cout << "WoTrackTuner file: " << woTrackTunerFile.Data() << std::endl;
    std::cout << "Histogram name: " << histName.Data() << std::endl;
    std::cout << "Number of bins: " << hRatio1->GetNbinsX() << std::endl;
    std::cout << "X-axis range: " << hRatio1->GetXaxis()->GetXmin() 
              << " - " << hRatio1->GetXaxis()->GetXmax() << " GeV/c" << std::endl;
    
    // 평균 ratio 계산 (첫 번째 ratio)
    double sumRatio1 = 0;
    int validBins1 = 0;
    for (int i = 1; i <= hRatio1->GetNbinsX(); i++) {
        if (hRatio1->GetBinContent(i) > 0) {
            sumRatio1 += hRatio1->GetBinContent(i);
            validBins1++;
        }
    }
    if (validBins1 > 0) {
        double avgRatio1 = sumRatio1 / validBins1;
        std::cout << "Average ratio (Tracking Efficiency / Reference): " << avgRatio1 
                  << " (from " << validBins1 << " valid bins)" << std::endl;
    }
    
    // 평균 ratio 계산 (두 번째 ratio)
    double sumRatio2 = 0;
    int validBins2 = 0;
    for (int i = 1; i <= hRatio2->GetNbinsX(); i++) {
        if (hRatio2->GetBinContent(i) > 0) {
            sumRatio2 += hRatio2->GetBinContent(i);
            validBins2++;
        }
    }
    if (validBins2 > 0) {
        double avgRatio2 = sumRatio2 / validBins2;
        std::cout << "Average ratio (WoTrackTuner / Reference): " << avgRatio2 
                  << " (from " << validBins2 << " valid bins)" << std::endl;
    }
    
    // 파일 닫기
    fTrackingEff->Close();
    fReference->Close();
    fWoTrackTuner->Close();
    
    std::cout << "Ratio plots saved to plots/ratio/CrossSectionRatio.pdf and .png" << std::endl;
}
