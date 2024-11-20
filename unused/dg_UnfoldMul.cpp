#include "BSHelper.cxx"
#include "Filipad2.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include "dndetaCommon.h"

void hset(TH1 &hid, TString xtit = "", TString ytit = "", double titoffx = 0.9,
          double titoffy = 1.2, double titsizex = 0.06, double titsizey = 0.06,
          double labeloffx = 0.01, double labeloffy = 0.001,
          double labelsizex = 0.05, double labelsizey = 0.05, int divx = 510,
          int divy = 510) {
  hid.SetStats(0);

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
Int_t n = 0;

void dg_UnfoldMul() {
  std::cout<< Form("Debugger_%i", n++) << std::endl;
  // TString dataname="MulLHC15f";
  // TString dataname="MulLHC16l";
  // TString mcname = "MulLHC15fMCPYTHIA8";
  // TString mcname = "MulLHC16lMCPYTHIA8";

  // auto clistmc = LoaddndetaResultList(mcname.Data(), "output");
  // auto clistdata = LoaddndetaResultList(dataname.Data(), "output");
  auto infile = TFile::Open(
      "AnalysisResults_MC.root",
      "open");
  auto infile2 = TFile::Open("AnalysisResults_Data.root",
                             "open");
  std::cout<< Form("Debugger_%i", n++) << std::endl;
  auto indir = (TDirectory *)infile->Get("multiplicity_dg");
  std::cout<< Form("Debugger_%i", n++) << std::endl;

  auto indir2 = (TDirectory *)infile2->Get("multiplicity_dg");
  if (!infile2) {
    std::cout << "no infile2" << std::endl;
  }
  /**********************************************************************************/

  auto recmc = (TH1D *)indir->Get("h_gen_particle_pt");
  auto truept = (TH1D *)indir->Get("h_recon_track_pt");
  auto res = (TH2D *)indir->Get("h2_matched_track_pt_particle_pt");
  std::cout<< Form("Debugger_%i", n++) << std::endl;
  /**********************************************************************************/

  auto matched_rec = (TH1D *) res->ProjectionX("matched_rec", 1, res->GetNbinsY());
  auto fake = (TH1D *)truept->Clone("fake");
  fake->Add(matched_rec, -1);
  auto matched_gen = (TH1D *) res->ProjectionY("matched_gen", 1, res->GetNbinsX());
  auto miss = (TH1D *)recmc->Clone("miss");
  miss->Add(matched_gen, -1);
  std::cout<< Form("Debugger_%i", n++) << std::endl;
  /**********************************************************************************/

  // auto miss = (TH1D*) clistmc->FindObject("hMultResponseMiss");
  // auto fake = (TH1D*) clistmc->FindObject("hMultResponseFake");
  auto recdata = (TH1D *)indir2->Get("h_data_track_pt");
  std::cout<< Form("Debugger_%i", n++) << std::endl;

  RooUnfoldResponse *Response = new RooUnfoldResponse(recmc, truept);

  for (auto i = 1; i <= res->GetNbinsX(); i++) {
    for (auto j = 1; j <= res->GetNbinsY(); j++) { // ptpair
      Double_t bincenx = res->GetXaxis()->GetBinCenter(i);
      Double_t binceny = res->GetYaxis()->GetBinCenter(j);
      Double_t bincont = res->GetBinContent(i, j);
      Response->Fill(bincenx, binceny, bincont);
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

  std::cout<< Form("Debugger_%i", n++) << std::endl;
  RooUnfoldBayes unfold(Response, recdata, 1);
  auto hcorrected = (TH1D *)unfold.Hreco();


  new TCanvas;

  Filipad2 *canvas = new Filipad2(1, 2, 0.5, 100, 50, 0.7, 1, 1);
  canvas->Draw();
  TPad *p = canvas->GetPad(1); // upper pad
  p->SetTickx();
  p->SetGridy(0);
  p->SetGridx(0);
  p->SetLogy(1);
  p->SetLogx(0);
  p->cd();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TLegend *leg =
      new TLegend(0.562201, 0.238261, 0.930622, 0.387826, "", "brNDC");
  leg->SetTextSize(0.062);
  leg->SetBorderSize(0);

  hcorrected->SetLineColor(1);
  hcorrected->Scale(1. / hcorrected->Integral(), "width");
  // hcorrected->SetMarkerStyle(24);
  hset(*hcorrected, "#it{N}_{ch}", "1/#it{N}_{evt} d#it{N}_{evt}/d#it{N}_{ch}",
       0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 510);
  hcorrected->GetXaxis()->SetRangeUser(0, 40);
  cout << hcorrected->GetMean() << endl;
  hcorrected->Draw("");
  
  recdata->SetLineColor(3);
  recdata->Scale(1. / recdata->Integral(), "width");
  // hcorrected->SetMarkerStyle(24);
  recdata->GetXaxis()->SetRangeUser(0, 40);
  cout << recdata->GetMean() << endl;
  recdata->Draw("same");

  truept->Scale(1. / truept->Integral(), "width");
  truept->SetLineColor(2);
  truept->SetMarkerColor(2);
  truept->SetLineWidth(2);
  truept->SetMarkerSize(3);
  truept->Draw("histpsame");

  leg->AddEntry(hcorrected, "Unfolded DATA", "l");
  leg->AddEntry(recdata, "raw DATA", "l");
  leg->AddEntry(truept, "PYTHIA8", "l");

  leg->Draw();

  p = canvas->GetPad(2); // lower pad
  p->SetTickx();
  p->SetGridy(1);
  p->SetGridx(0);
  p->SetLogy(0);
  ;
  p->SetLogx(0);
  p->cd();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  auto ratio = (TH1D *)hcorrected->Clone();
  ratio->Divide(truept);
  hset(*ratio, "#it{N}_{ch}", "DATA/PYTHIA8", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01,
       0.07, 0.07, 515, 505);
  ratio->GetXaxis()->SetRangeUser(0, 40);
  ratio->Draw("");

  TFile *f = new TFile("Mul.root", "recreate");
  truept->SetName("PYTHIA8");
  hcorrected->SetName("DATA");
  truept->Write();
  hcorrected->Write();
  f->Close();

  new TCanvas;
}