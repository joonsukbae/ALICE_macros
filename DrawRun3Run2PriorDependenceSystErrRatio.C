#include "Filipad2.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"
#include <iostream>
#include <vector>

using namespace std;

void DrawRun3Run2PriorDependenceSystErrRatio() {
  gStyle->SetOptStat(0);
  
  // 분모 데이터 (Default Dataset.csv)
  vector<double> x_denominator = {
    4.5, 5.5, 6.5, 7.5,
    8.5, 9.5, 11.0, 13.0,
    15.0, 17.0, 19.0, 22.5,
    27.5, 35.0, 45.0
  };
  
  vector<double> y_denominator = {
    0.9988486312399352, 0.9988486312399352, 0.9988486312399352, 1.0018679549114327,
    0.9988486312399352, 1.0018679549114327, 0.9887842190016098, 0.9948228663446053,
    0.9978421900161027, 1.013945249597423, 1.0219967793880833, 1.0240096618357484,
    1.027028985507246, 1.0300483091787436, 0.9535587761674715
  };
  
  // 분자 데이터 (Default Dataset 2.csv)
  vector<double> x_numerator = {
    5.5, 6.5, 7.5, 8.5,
    9.5, 11.0, 13.0, 15.0,
    17, 19, 22.5, 27.5,
    35., 45., 55, 65,
    77.5, 92.5, 120.0
  };
  
  vector<double> y_numerator = {
    0.9693950177935946, 0.9779359430604986, 0.9928825622775799, 1.0035587188612096,
    1.0163701067615651, 1.0355871886120989, 1.0633451957295375, 1.0825622775800712,
    1.084697508896797, 1.0825622775800712, 1.0654804270462632, 1.0483985765124562,
    1.0014234875444838, 0.9715302491103204, 0.9715302491103204, 0.9693950177935946,
    0.9202846975088974, 0.7836298932384338, 0.9032028469750895
  };
  
  // Ratio 계산을 위해 x 값이 겹치는 지점 찾기
  // 분자와 분모의 x 값이 정확히 일치하지 않으므로, 가장 가까운 점을 찾아서 보간하거나
  // 또는 각각의 x 범위에서 ratio를 계산
  
  // 방법 1: 분모의 x 값에 가장 가까운 분자 값을 찾아서 ratio 계산
  vector<double> x_ratio;
  vector<double> y_ratio;
  vector<double> y_ratio_err;
  
  for (size_t i = 0; i < x_denominator.size(); i++) {
    double x_den = x_denominator[i];
    double y_den = y_denominator[i];
    
    // 분자 데이터에서 가장 가까운 x 값 찾기
    double min_dist = 1e10;
    size_t closest_idx = 0;
    for (size_t j = 0; j < x_numerator.size(); j++) {
      double dist = TMath::Abs(x_numerator[j] - x_den);
      if (dist < min_dist) {
        min_dist = dist;
        closest_idx = j;
      }
    }
    
    // 거리가 너무 멀면 스킵 (예: 5 이상 차이나면)
    if (min_dist > 5.0) continue;
    
    double y_num = y_numerator[closest_idx];
    
    // 분모가 0이 아니면 ratio 계산
    if (y_den != 0) {
      double ratio = y_num / y_den;
      x_ratio.push_back(x_den);
      y_ratio.push_back(ratio);
      // 간단한 오차 전파 (더 정확한 오차 계산이 필요하면 수정)
      y_ratio_err.push_back(ratio * 0.01); // 임시로 1% 오차
    }
  }
  
  // 방법 2: 분자의 x 값에 가장 가까운 분모 값을 찾아서 ratio 계산 (더 많은 데이터 포인트)
  vector<double> x_ratio2;
  vector<double> y_ratio2;
  vector<double> y_ratio_err2;
  
  for (size_t i = 0; i < x_numerator.size(); i++) {
    double x_num = x_numerator[i];
    double y_num = y_numerator[i];
    
    // 분모 데이터에서 가장 가까운 x 값 찾기
    double min_dist = 1e10;
    size_t closest_idx = 0;
    for (size_t j = 0; j < x_denominator.size(); j++) {
      double dist = TMath::Abs(x_denominator[j] - x_num);
      if (dist < min_dist) {
        min_dist = dist;
        closest_idx = j;
      }
    }
    
    // 거리가 너무 멀면 스킵
    if (min_dist > 5.0) continue;
    
    double y_den = y_denominator[closest_idx];
    
    // 분모가 0이 아니면 ratio 계산
    if (y_den != 0) {
      double ratio = y_num / y_den;
      x_ratio2.push_back(x_num);
      y_ratio2.push_back(ratio);
      y_ratio_err2.push_back(ratio * 0.01); // 임시로 1% 오차
    }
  }
  
  // 더 많은 데이터 포인트를 가진 방법 2 사용
  int nPoints = x_ratio2.size();
  
  // 분자와 분모 그래프 생성
  int nNum = x_numerator.size();
  int nDen = x_denominator.size();
  
  TGraphErrors *gNumerator = new TGraphErrors(nNum);
  TGraphErrors *gDenominator = new TGraphErrors(nDen);
  
  for (int i = 0; i < nNum; i++) {
    gNumerator->SetPoint(i, x_numerator[i], y_numerator[i]);
    gNumerator->SetPointError(i, 0, 0);
  }
  
  for (int i = 0; i < nDen; i++) {
    gDenominator->SetPoint(i, x_denominator[i], y_denominator[i]);
    gDenominator->SetPointError(i, 0, 0);
  }
  
  // Ratio 그래프 생성
  TGraphErrors *gRatio = new TGraphErrors(nPoints);
  gRatio->SetName("gRatio");
  gRatio->SetTitle("");
  
  for (int i = 0; i < nPoints; i++) {
    gRatio->SetPoint(i, x_ratio2[i], y_ratio2[i]);
    gRatio->SetPointError(i, 0, y_ratio_err2[i]);
  }
  
  // Filipad2 생성
  Int_t nn = 0;
  Filipad2 *pad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
  pad->Draw();
  
  // 위 패드 (메인 플롯)
  TPad *toppad = pad->GetPad(1);
  toppad->cd();
  toppad->SetGridx(1);
  toppad->SetGridy(1);
  
  // 분모 그래프 스타일 설정
  gDenominator->SetMarkerColor(kBlue);
  gDenominator->SetLineColor(kBlue);
  gDenominator->SetMarkerStyle(20);
  gDenominator->SetMarkerSize(1.2);
  gDenominator->SetLineWidth(2);
  
  // 분자 그래프 스타일 설정
  gNumerator->SetMarkerColor(kRed);
  gNumerator->SetLineColor(kRed);
  gNumerator->SetMarkerStyle(20);
  gNumerator->SetMarkerSize(1.2);
  gNumerator->SetLineWidth(2);
  gNumerator->GetYaxis()->SetRangeUser(0.8, 1.2);
  
  // 분모 그래프 먼저 그리기
  gDenominator->GetXaxis()->SetTitle("x");
  gDenominator->GetYaxis()->SetTitle("Relative uncertainty (Prior dep.)");
  gDenominator->GetXaxis()->SetTitleOffset(1.2);
  gDenominator->GetYaxis()->SetTitleOffset(1.3);
  gDenominator->GetXaxis()->SetTitleSize(0.05);
  gDenominator->GetYaxis()->SetTitleSize(0.05);
  gDenominator->GetXaxis()->SetLabelSize(0.04);
  gDenominator->GetYaxis()->SetLabelSize(0.04);
  gDenominator->GetYaxis()->SetRangeUser(0.8, 1.2);
  gDenominator->Draw("AP");
  
  // 분자 그래프 그리기
  gNumerator->Draw("P SAME");
  
  // 위 패드 범례
  TLegend *legTop = new TLegend(0.3, 0.8, 0.6, 0.95);
  legTop->SetBorderSize(0);
  legTop->SetFillStyle(0);
  legTop->SetTextSize(0.045);
  legTop->AddEntry(gDenominator, "Run 2 (Epos / PYTHIA8)", "ple");
  legTop->AddEntry(gNumerator, "Run 3 Preliminary (Herwig / PYTHIA8)", "ple");
  legTop->Draw();
  
  // 아래 패드 (Ratio 플롯)
  TPad *ratiopad = pad->GetPad(2);
  ratiopad->cd();
  ratiopad->SetGridx(1);
  ratiopad->SetGridy(1);
  
  // Ratio 그래프 스타일 설정
  gRatio->SetMarkerColor(kRed);
  gRatio->SetLineColor(kRed);
  gRatio->SetMarkerStyle(20);
  gRatio->SetMarkerSize(1.2);
  gRatio->SetLineWidth(2);
  
  // Ratio 축 설정
  gRatio->GetXaxis()->SetTitle("#it{p}_{T, jet}^{truth} (GeV/c)");
  gRatio->GetYaxis()->SetTitle("Run 3 / Run 2");
  gRatio->GetXaxis()->SetTitleOffset(1.2);
  gRatio->GetYaxis()->SetTitleOffset(1.3);
  gRatio->GetXaxis()->SetTitleSize(0.05);
  gRatio->GetYaxis()->SetTitleSize(0.05);
  gRatio->GetXaxis()->SetLabelSize(0.04);
  gRatio->GetYaxis()->SetLabelSize(0.04);
  
  // Y 범위 설정 (1.0 기준선이 잘 보이도록)
  double ymin = 0.85;
  double ymax = 1.15;
  for (int i = 0; i < nPoints; i++) {
    double yval = y_ratio2[i];
    if (yval < ymin) ymin = yval;
    if (yval > ymax) ymax = yval;
  }
  ymin = ymin - 0.1 * (ymax - ymin);
  ymax = ymax + 0.1 * (ymax - ymin);
  if (ymin < 0) ymin = 0;
  gRatio->GetYaxis()->SetRangeUser(ymin, ymax);
  
  // Ratio 그리기
  gRatio->Draw("AP");
  
  // 1.0 기준선 추가
  TLine *line1 = new TLine(gRatio->GetXaxis()->GetXmin(), 1.0, 
                           gRatio->GetXaxis()->GetXmax(), 1.0);
  line1->SetLineColor(kBlack);
  line1->SetLineStyle(2);
  line1->SetLineWidth(1);
  line1->Draw("SAME");
  
  // 아래 패드 범례
  TLegend *legRatio = new TLegend(0.7, 0.85, 0.9, 0.95);
  legRatio->SetBorderSize(0);
  legRatio->SetFillStyle(0);
  legRatio->SetTextSize(0.06);
  legRatio->AddEntry(gRatio, "Run 3 / Run 2", "ple");
  legRatio->Draw();
  
  // 저장
  pad->C->SaveAs("ratio_plot.pdf");
  pad->C->SaveAs("ratio_plot.png");
  
  cout << "Ratio plot saved to ratio_plot.pdf and ratio_plot.png" << endl;
  cout << "Number of ratio points: " << nPoints << endl;
}
