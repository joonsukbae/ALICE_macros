#include "Filipad2.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <iostream>
#include <string>
#include <vector>

#include <RooUnfoldBayes.h>
#include <RooUnfoldResponse.h>

// Double_t ptbin[22] = {5,  6,  7,  8,  9,  10, 12, 14,  16,  18,  20,
//                       25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
// Double_t ptbinGen[27] = {0,  1,  2,  3,  4,  5,   6,   7,   8,
//                          9,  10, 12, 14, 16, 18,  20,  25,  30,
//                          40, 50, 60, 70, 85, 100, 140, 200};

Double_t ptbin[21] = {5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
Double_t ptbinGen[26] = {0, 1, 2, 3, 4, 5,  6,  7,  8,  9,  10, 12, 14, 16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
Int_t nptBins = sizeof(ptbin) / sizeof(ptbin[0]) - 1;
Int_t nptBinsGen = sizeof(ptbinGen) / sizeof(ptbinGen[0]) - 1;

TString JetPtTitleX = "#it{p}_{T, jet}^{reco} (GeV/c)";
TString JetPtTitleY = "1/N_{evt} dN/d#it{p}_{T}";

const char *fileName;
const char *Dataset;
const char *eventDir = "jet-finder-charged-qa";
const char *eventObj = "h_collisions";

Int_t n = 0;
Int_t nn = 0;

template <typename T>
void hset(T &hid, TString xtit = "", TString ytit = "", double titoffx = 0.9,
          double titoffy = 1.2, double titsizex = 0.06, double titsizey = 0.06,
          double labeloffx = 0.01, double labeloffy = 0.001,
          double labelsizex = 0.05, double labelsizey = 0.05, int divx = 510,
          int divy = 510) {
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

template <typename T>
void hoptset(T &hid, Double_t N = 1, Color_t color = kBlack, Double_t minX = 0,
             Double_t maxX = 100, Double_t minY = 0, Double_t maxY = 1,
             Double_t MarkerSize = 0.8) {
  if (N != 0) {
    if (N == 1) {
      hid.Scale(1. / hid.Integral(), "width");
    } else if (N == 2) {
      hid.Scale(1. / hid.Integral(), "");
    } else {
      hid.Scale(1. / N, "width");
    }
  }

  hid.SetMarkerColor(color);
  hid.SetLineColor(color);
  hid.SetLineWidth(2);
  hid.SetMarkerSize(MarkerSize);
  hid.SetMarkerStyle(20);

  hid.GetXaxis()->SetRangeUser(minX, maxX);
  hid.GetYaxis()->SetRangeUser(minY, maxY);
  hid.SetFillColorAlpha(color, 0.3);
  hid.GetYaxis()->SetNdivisions(505);
}

void optFili(TPad &pid, Int_t gridx, Int_t gridy, Int_t logx, Int_t logy) {
  pid.SetGridy(gridx);
  pid.SetGridx(gridy);
  pid.SetLogx(logx);
  pid.SetLogy(logy);
}

// Function prototypes
void LoadMCData(const std::string &fileName, TH1 *&hTrue, TH1 *&hReco,
                TH2 *&hMatch);
void DoCreateResponseMatrix(TH1 *hTrue1, TH1 *hReco1, TH2 *hMatch1,
                            RooUnfoldResponse *&Response);
void DoPerformUnfolding(RooUnfoldResponse *Response, TH1 *hReco, int &Niteration,
                        TH1 *&hUnfolded);
void DoPerformRefolding(RooUnfoldResponse *Response, TH1 *hTrue, TH1 *&hRefolded);
double GetNevent(const std::string &fileName, const char *eventDir,
              const char *eventObj, bool ifMCP = false);
void DoCompareHistograms(TH1 *hTrue, TH1 *hUnfolded, double Nevt, TString outputName, TString hist1name, TString hist2name, TString ratiotitle);

void DoClosureTest(const std::string &fileName1, const std::string &fileName2) {
  // Load data from files
  TH1 *hTrue1, *hReco1, *hTrue2, *hReco2;
  TH2 *hMatch1, *hMatch2;
  LoadMCData(fileName1, hTrue1, hReco1, hMatch1);
  LoadMCData(fileName2, hTrue2, hReco2, hMatch2);

  double Ntrue = GetNevent(fileName2, eventDir, eventObj, true);
  double Nreco = GetNevent(fileName2, eventDir, eventObj, false);
  std::cout << "(Ntrue, Nreco): " << Ntrue << ", " << Nreco << std::endl;

  // Create response matrix
  RooUnfoldResponse *Response;
  DoCreateResponseMatrix(hTrue1, hReco1, hMatch1, Response);

  // Perform unfolding on the second sample
  TH1 *hUnfolded;
  int bestIteration = 4; // default value
  DoPerformUnfolding(Response, hReco2, bestIteration, hUnfolded);
  DoCompareHistograms(hTrue2, hUnfolded, Ntrue, "MCunfolded.pdf", "True MC", "Unfolded", "Unfolded / True MC");

  TH1 *hRefolded;
  DoPerformRefolding(Response, hUnfolded, hRefolded);
  DoCompareHistograms(hReco2, hRefolded, Nreco, "MCrefolded.pdf", "Reco MC", "Refolded", "Refolded / Reco MC");

  // Clean up
  delete hTrue1;
  delete hReco1;
  delete hTrue2;
  delete hReco2;
  delete hMatch1;
  delete hMatch2;
  delete Response;
  delete hUnfolded;
  delete hRefolded;
}

// Function to load data from file
void LoadMCData(const std::string &fileName, TH1 *&hTrue, TH1 *&hReco,
                TH2 *&hMatch) {
  TFile *file = TFile::Open(fileName.c_str(), "open");
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file: " << fileName << std::endl;
    return;
  }

  hTrue = (TH1 *)file->Get("jet-finder-charged-qa/h_jet_pt_part");
  
  if (hTrue) {
    hTrue = hTrue->Rebin(nptBinsGen, "hTrueRebin", ptbinGen);
  } else {
    std::cerr << "Failed to load h_jet_pt_part" << std::endl;
    return;
  }

  hReco = (TH1 *)file->Get("jet-finder-charged-qa/h_jet_pt");
  if (hReco) {
    hReco = hReco->Rebin(nptBins, "hRecoRebin", ptbin);
  } else {
    std::cerr << "Failed to load h_jet_pt" << std::endl;
    return;
  }

  auto h3Match = (TH3 *)file->Get("jet-finder-charged-qa/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo");

  if (!h3Match) {
    std::cerr << "Failed to load h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo" << std::endl;
    return;
  }

    hMatch = new TH2D("hcorrelate", "R projected correlate", nptBins, ptbin, nptBinsGen, ptbinGen);

  Int_t corrbin = h3Match->GetXaxis()->FindBin(0.4 + 1e-6);
  for (Int_t i = 1; i <= h3Match->GetNbinsY(); i++) {   // part.
    for (Int_t j = 1; j <= h3Match->GetNbinsZ(); j++) { // det.
      Double_t content = h3Match->GetBinContent(corrbin, i, j);
      Double_t error = h3Match->GetBinError(corrbin, i, j);

      Int_t binpart = hMatch->GetYaxis()->FindBin(h3Match->GetYaxis()->GetBinCenter(i));
      Int_t bin = hMatch->GetXaxis()->FindBin(h3Match->GetZaxis()->GetBinCenter(j));

      Double_t currentContent = hMatch->GetBinContent(bin, binpart);
      Double_t currentError = hMatch->GetBinError(bin, binpart);

      Double_t newContent = currentContent + content;
      Double_t newError = sqrt(pow(currentError, 2) + pow(error, 2)); // Combine errors in quadrature

      hMatch->SetBinContent(bin, binpart, newContent);
      // for (auto k = 0;  k<newContent; k++){
      // hMatch->Fill(bin, binpart);
      // }
      hMatch->SetBinError(bin, binpart, newError);
    }
  }

  hTrue->SetDirectory(0);  // Detach from file
  hReco->SetDirectory(0);  // Detach from file
  hMatch->SetDirectory(0); // Detach from file
  file->Close();
  delete file;
}

// Function to create the response matrix
void DoCreateResponseMatrix(TH1 *hTrue1, TH1 *hReco1, TH2 *hMatch1,
                            RooUnfoldResponse *&Response) {
  auto hMatchReco =
      (TH1 *)hMatch1->ProjectionX("hRecoMatched", 1, hMatch1->GetNbinsY(), "e");
  auto hMatchTrue =
      (TH1 *)hMatch1->ProjectionY("hTrueMatched", 1, hMatch1->GetNbinsX(), "e");

  TH1F *fake = (TH1F *)hReco1->Clone(Form("%i", ++n));
  fake->Add(hMatchReco, -1);
  TH1F *miss = (TH1F *)hTrue1->Clone(Form("%i", ++n));
  miss->Add(hMatchTrue, -1);

  Response = new RooUnfoldResponse(hMatchReco, hMatchTrue);

  for (auto i = 1; i <= hMatch1->GetNbinsX(); i++) {
    for (auto j = 1; j <= hMatch1->GetNbinsY(); j++) { // ptpair
      Double_t bincenx = hMatch1->GetXaxis()->GetBinCenter(i);
      Double_t binceny = hMatch1->GetYaxis()->GetBinCenter(j);
      Double_t bincont = hMatch1->GetBinContent(i, j);
      // Response->Fill(bincenx, binceny, bincont);
      for (auto k = 0;  k<bincont; k++){
      Response->Fill(bincenx, binceny);
      }
    }
  }

  for (auto i = 1; i <= miss->GetNbinsX(); i++) {
    Double_t bincenx = miss->GetXaxis()->GetBinCenter(i);
    Double_t bincont = miss->GetBinContent(i);
    Response->Miss(bincenx, bincont);
  }

  for (auto i = 1; i <= fake->GetNbinsX(); i++) {
    Double_t bincenx = fake->GetXaxis()->GetBinCenter(i);
    Double_t bincont = fake->GetBinContent(i);
    Response->Fake(bincenx, bincont);
  }
}

// Function to perform unfolding
void DoPerformUnfolding(RooUnfoldResponse *Response, TH1 *hReco, int &Niteration,
                        TH1 *&hUnfolded) {
  double minChi2 = 1e10;
  int bestIteration = 1;

  for (int i = 1; i <= 10; ++i) {
    RooUnfoldBayes unfold(Response, hReco, i);
    TH1 *tempUnfolded = (TH1D *)unfold.Hreco();

    // Calculate chi2 between tempUnfolded and the true distribution
    double chi2 = tempUnfolded->Chi2Test(hReco, "UW CHI2");

    std::cout << "Iteration " << i << ": Chi2 = " << chi2 << std::endl;

    if (chi2 < minChi2) {
      minChi2 = chi2;
      bestIteration = i;
    }

    delete tempUnfolded;
  }

  std::cout << "Best iteration: " << bestIteration << std::endl;

  Niteration = bestIteration;
  RooUnfoldBayes unfold(Response, hReco, Niteration);
  hUnfolded = (TH1D *)unfold.Hreco();
}

// Function to perform refolding
void DoPerformRefolding(RooUnfoldResponse *Response, TH1 *hUnfolded, TH1 *&hRefolded) {
  hRefolded = (TH1 *)Response->ApplyToTruth(hUnfolded);
}

double GetNevent(const std::string &fileName, const char *eventDir,
                 const char *eventObj, bool ifMCP = false) {
  auto file = TFile::Open(fileName.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file: " << fileName << std::endl;
    return 0;
  }

  auto Nevents = (TH1D *)file->Get(Form("%s/%s", eventDir, eventObj));
  if (!Nevents) {
    std::cerr << "Error: Object not found in file: " << eventObj << std::endl;
    file->Close();
    delete file;
    return 0;
  }

  Double_t nevents = 0;
  if (ifMCP) {
    nevents = Nevents->GetBinContent(Nevents->FindBin(0.5));
  } else {
    nevents = Nevents->GetBinContent(Nevents->FindBin(1.5));
  }

  // 리소스 정리
  file->Close();
  delete file;

  return nevents;
}

// Function to compare histograms
void DoCompareHistograms(TH1 *hTrue, TH1 *hUnfolded, double Nevt, TString outputName, TString hist1name, TString hist2name, TString ratiotitle) {

  Filipad2 *ClosurePad = new Filipad2(++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
  ClosurePad->Draw();
  TPad *closurepad = ClosurePad->GetPad(1);
  optFili(*closurepad, 1, 1, 0, 1);
  TPad *ratioclosurepad = ClosurePad->GetPad(2);
  optFili(*ratioclosurepad, 1, 1, 0, 0);
  TLegend *legclosure =
      new TLegend(0.57177,0.530435,0.837321,0.773913,NULL,"brNDC");
  legclosure->SetTextSize(0.05);
  legclosure->SetBorderSize(0);

  closurepad->cd();
  auto hgen = (TH1 *) hTrue->Clone(Form("%i", ++n));
  hset(*hgen, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05,
       0.05, 510, 505);
  hoptset(*hgen, 0, kBlack, 5, 140, 1e1, 1e10, 1.2);
  auto hcorr = (TH1 *) hUnfolded->Clone(Form("%i", ++n));
  hset(*hcorr, JetPtTitleX, JetPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 505);
  hoptset(*hcorr, 0, kRed, 5, 140, 1e1, 1e10);
  
  hgen->Scale(1./Nevt, "width");
  hcorr->Scale(1./Nevt, "width");
  hgen->Draw("pe");
  hcorr->Draw("pesame");

  legclosure->AddEntry(hgen, hist1name, "pe");
  legclosure->AddEntry(hcorr, hist2name, "pe");
  legclosure->Draw();

  ratioclosurepad->cd();
  auto ratclosure = (TH1 *)hcorr->Clone(Form("%i", ++n));
  ratclosure->Divide(ratclosure, hgen, 1., 1., "B");
  hset(*ratclosure, JetPtTitleX, ratiotitle, 1.2, 1.0, 0.07, 0.07,
       0.01, 0.01, 0.07, 0.07, 510, 505);
  hoptset(*ratclosure, 0, kRed, 5, 140, 0.4, 1.6);
  ratclosure->Draw("pe");

  TString folderPath = Form("plots/AN_Charged-particle-jet-cross-section-in-pp-collisions-at-13.6-TeV/Figures/MCclosure/%s", Dataset);
  gSystem->Exec(Form("mkdir -p %s", folderPath.Data()));
  ClosurePad->C->SaveAs(Form("%s/%s", folderPath.Data(), outputName.Data()));
}

int DrawMCCloureTest() {
  // Example file names
  // Dataset 설정
  Dataset = "LHC24f3b_trackTuner_Track100GeV";
  std::string fileName1 =
      // "~/cernbox/workspace/O2Physics/jets/mc/AnalysisResults/"
      // // "LHC24f3/selMC/MCClosureTestSamples/"
      // "LHC24f3_trackTuner/selMC/MCClosureTestSamples/"
      // "merged_AnalysisResults1.root";
      "/Users/js/cernbox/workspace/O2Physics/jets/mc/AnalysisResults/LHC24f3b/selMC/TrackTuner/Track100GeV/MCClosureTestSamples/merged_AnalysisResults1.root";
  std::string fileName2 =
      // "~/cernbox/workspace/O2Physics/jets/mc/AnalysisResults/"
      // // "LHC24f3/selMC/MCClosureTestSamples/"
      // "LHC24f3_trackTuner/selMC/MCClosureTestSamples/"
      // "merged_AnalysisResults2.root";
      "/Users/js/cernbox/workspace/O2Physics/jets/mc/AnalysisResults/LHC24f3b/selMC/TrackTuner/Track100GeV/MCClosureTestSamples/merged_AnalysisResults2.root";

  // Run closure test
  DoClosureTest(fileName1, fileName2);

  return 0;
}