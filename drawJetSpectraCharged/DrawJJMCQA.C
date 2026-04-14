///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// JJ MC QA macro                //////////
////////// MB MC vs JJ MC comparison     //////////
////////// + unfolded spectra with data  //////////
////////// author: Joonsuk Bae           //////////
////////// E-mail: jbae@cern.ch          //////////
////////// Created: 06 Mar 2026          //////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

#include "BSHelper.cxx"
#include "Filipad2.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldSvd.h"
#include "TSVDUnfold_local.h"
#include <iostream>
#include <vector>
using namespace std;

///////////////////
/// plot switch ///
///////////////////
// Track-level
const bool kPlotTrackPt       = true;   // Track pT (from processTracksQC)
const bool kPlotTrackEtaPhi   = true;   // Track eta, phi
const bool kPlotTrackingEff   = true;   // Tracking efficiency (from track-efficiency task)
const bool kPlotTrackPtRes    = true;   // Track pT resolution (from track-efficiency task)
// Jet-level
const bool kPlotJetPtMCD      = true;   // Jet pT reco (MB vs JJ)
const bool kPlotJetPtMCP      = true;   // Jet pT particle (MB vs JJ)
const bool kPlotJetPtFineBin  = true;   // Fine-binned versions
const bool kPlotConstituentPt = true;   // Constituent pT (MCD and MCP)
const bool kPlotJetEtaPhi     = true;   // Jet eta, phi
const bool kPlotPtHat         = true;   // pTHat distribution after scaling
const bool kPlotResponseMatrix = true;  // Response matrix (2D colz)
const bool kPlotPurityEff     = true;   // Purity and efficiency
const bool kPlotKinEff        = true;   // Kinematic efficiency (response matrix window method)
const bool kPlotJetResolution = true;   // Jet pT resolution in pT classes
const bool kPlotJetNtracks    = true;   // Jet ntracks
// Unfolded
const bool kPlotUnfolded      = true;   // Unfolded raw pT spectrum comparison

const bool DRAWPLOTS          = true;
const bool NORMEVENTS         = true;
const bool REBINON            = true;

// Normalization point for shape comparison (GeV)
const Double_t kNormPt = 20.0;
// SVD regularization (data-driven d-vector in DrawJetMatching)
const Int_t kSVDkReg = 4;

///////////////////
/// Files setup ///
///////////////////
const TString mainDir = "~/cernbox/workspace/O2Physics/jets/AnalysisResults/";
const TString fRoot = "_AnalysisResults.root";

// MB MC
const TString kMBFile = "LHC23k4i_JJMCQA" + fRoot;
const TString kMBDir  = "jet-spectra-charged";
const TString kMBName = "2023 MB MC (DCA-only, noAssoc)";

// JJ MC
const TString kJJFile = "LHC26b5_JJMCQA" + fRoot;
const TString kJJDir  = "jet-spectra-charged";
const TString kJJName = "2023 JJ MC (DCA-only, noAssoc)";

// Data (for unfolded comparison)
const TString kDataFile = "496213" + fRoot;  // 2023 pass4 thin
const TString kDataDir  = "jet-spectra-charged";
const TString kDataName = "2023 data (LHC23_pass4_Thin)";

const TString PlotSaveName = "JJMC_QA_DCAonly";
const TString MakeDirName = "../plots/"
                            "AN_Charged-particle-jet-cross-section-in-pp-"
                            "collisions-at-13.6-TeV/Figures/" +
                            PlotSaveName;

///////////////////
/// Histogram objects from jetSpectraCharged
///////////////////
const char *EventObj     = "h_collisions";
const char *EventWObj    = "h_collisions_weighted";
const char *JetPt3DObj   = "h3_jet_pt_jet_eta_jet_phi";
const char *JetPt1DObj   = "h_jet_pt";
const char *JetPtMCPObj  = "h_jet_pt_part";
const char *TrackPtObj   = "h_track_pt";
const char *TrackEtaPhiObj = "h2_track_eta_track_phi";
const char *ConstPtMCDObj  = "h2_jet_pt_track_pt";
const char *ConstPtMCPObj  = "h2_jet_pt_part_track_pt_part";
const char *JetNtracksObj  = "h2_jet_pt_jet_ntracks";
const char *ResponseObj  = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint";
const char *ResolutionObj = "h2_jet_pt_mcp_jet_pt_diff_matchedgeo";
const char *PtHatCollObj = "h_coll_phat_weighted";
const char *PtHatMCCollObj = "h_mccoll_phat_weighted";

// track-efficiency task histograms
const char *kTrkEffDir = "track-efficiency";
const char *kHistTruthLo  = "h3_particle_pt_particle_eta_particle_phi_mcpartofinterest";
const char *kHistMatchLo  = "h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary";
const char *kHistTruthHi  = "h3_particle_pt_high_particle_eta_particle_phi_mcpartofinterest";
const char *kHistMatchHi  = "h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_primary";
const char *kHistTrackPtRes = "h2_particle_pt_track_pt_deltaptoverparticlept";

///////////////////
/// Bins        ///
///////////////////
const Double_t RBIN = 0.4;
const double PlotPtMin = 5;
const double PlotPtMax = 200;

const Double_t ptbin[23] = {3,  4,  5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20,
                            25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
const Int_t nptBins = sizeof(ptbin) / sizeof(ptbin[0]) - 1;

const Double_t ptbinGen[26] = {0,  1,  2,  3,  4,  5,   6,   7,  8,
                               9,  10, 12, 14, 16, 18,  20,  25, 30,
                               40, 50, 60, 70, 85, 100, 140, 200};
const Int_t nptBinsGen = sizeof(ptbinGen) / sizeof(ptbinGen[0]) - 1;

const Double_t Trackptbin[21] = {0.15, 2, 4, 6, 8, 10, 15, 20, 25, 30,
                                 40, 50, 60, 70, 80, 90, 100, 120, 140, 170, 200};
const Int_t nTrackptbin = sizeof(Trackptbin) / sizeof(Trackptbin[0]) - 1;

// pT classes for jet resolution
const double kResPtEdges[] = {10, 20, 40, 60, 100, 200};
const int kNResPtBins = sizeof(kResPtEdges) / sizeof(kResPtEdges[0]) - 1;
const Color_t kResPtColors[] = {kBlue + 1, kRed, kGreen + 2, kOrange + 7, kViolet + 2};

///////////////////
/// Style       ///
///////////////////
Int_t nn = 0;

template <typename T>
void hset(T &hid, TString xtit = "", TString ytit = "", double titoffx = 0.9,
          double titoffy = 1.2, double titsizex = 0.06, double titsizey = 0.06,
          double labeloffx = 0.01, double labeloffy = 0.001,
          double labelsizex = 0.05, double labelsizey = 0.05, int divx = 510,
          int divy = 510) {
  hid.GetXaxis()->CenterTitle(1);  hid.GetYaxis()->CenterTitle(1);
  hid.GetXaxis()->SetTitleOffset(titoffx);  hid.GetYaxis()->SetTitleOffset(titoffy);
  hid.GetXaxis()->SetTitleSize(titsizex);  hid.GetYaxis()->SetTitleSize(titsizey);
  hid.GetXaxis()->SetLabelOffset(labeloffx);  hid.GetYaxis()->SetLabelOffset(labeloffy);
  hid.GetXaxis()->SetLabelSize(labelsizex);  hid.GetYaxis()->SetLabelSize(labelsizey);
  hid.GetXaxis()->SetNdivisions(divx);  hid.GetYaxis()->SetNdivisions(divy);
  hid.GetXaxis()->SetTitle(xtit);  hid.GetYaxis()->SetTitle(ytit);
}

void optFili(TPad &pid, Int_t gridx, Int_t gridy, Int_t logx, Int_t logy) {
  pid.SetGridx(gridx);  pid.SetGridy(gridy);
  pid.SetLogx(logx);    pid.SetLogy(logy);
}

void StyleHist(TH1 *h, Color_t color, Int_t marker = 20, Double_t size = 0.8) {
  h->SetLineColor(color);  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker); h->SetMarkerSize(size); h->SetLineWidth(2);
}

//=============================================================================
// Get number of events
//=============================================================================
Double_t GetNevents(const char *fileName, const char *dir, Int_t ifMCP = 0) {
  TFile *file = TFile::Open(fileName, "read");
  if (!file || file->IsZombie()) return 1.0;
  Double_t nevents = 1.0;
  if (ifMCP == 0) {
    TH1D *hW = (TH1D *)file->Get(Form("%s/%s", dir, EventWObj));
    TH1D *hU = (TH1D *)file->Get(Form("%s/%s", dir, EventObj));
    TH1D *h = hW ? hW : hU;
    if (h) nevents = h->GetBinContent(h->FindBin(2.5));
  } else {
    TH1D *h3 = (TH1D *)file->Get(Form("%s/h_mcColl_counts_weight", dir));
    TH1D *h4 = (TH1D *)file->Get(Form("%s/h_mcColl_counts", dir));
    TH1D *h5 = (TH1D *)file->Get(Form("%s/h_mccollisions", dir));
    TH1D *h6 = (TH1D *)file->Get(Form("%s/h_mccollisions_weighted", dir));
    TH1D *h = h3 ? h3 : (h4 ? h4 : (h6 ? h6 : h5));
    if (h) nevents = h->GetBinContent(4);
  }
  file->Close();
  return nevents;
}

//=============================================================================
// Generic histogram loader (returns detached clone)
//=============================================================================
TH1 *LoadHist1D(const char *fileName, const char *dir, const char *obj, const char *tag) {
  TFile *file = TFile::Open(fileName, "read");
  if (!file || file->IsZombie()) return nullptr;
  TObject *raw = file->Get(Form("%s/%s", dir, obj));
  if (!raw) { file->Close(); return nullptr; }
  TH1 *h = nullptr;
  if (raw->InheritsFrom(TH3::Class())) {
    TH3 *h3 = (TH3 *)raw;
    h = h3->ProjectionX(Form("h_%s", tag), 1, h3->GetNbinsY(), 1, h3->GetNbinsZ());
  } else if (raw->InheritsFrom(TH1::Class())) {
    h = (TH1 *)raw->Clone(Form("h_%s", tag));
  }
  if (h) h->SetDirectory(0);
  file->Close();
  return h;
}

TH2 *LoadHist2D(const char *fileName, const char *dir, const char *obj, const char *tag) {
  TFile *file = TFile::Open(fileName, "read");
  if (!file || file->IsZombie()) return nullptr;
  TH2 *h = (TH2 *)file->Get(Form("%s/%s", dir, obj));
  if (!h) { file->Close(); return nullptr; }
  h = (TH2 *)h->Clone(Form("h2_%s", tag));
  h->SetDirectory(0);
  file->Close();
  return h;
}

//=============================================================================
// Load jet pT histogram (normalized, optionally rebinned)
//=============================================================================
TH1 *LoadJetPt(const char *fileName, const char *dir, const char *obj,
               Double_t Nevts, bool rebin, const char *tag) {
  TH1 *h = LoadHist1D(fileName, dir, obj, tag);
  if (!h) return nullptr;
  if (rebin) {
    bool isMCP = TString(obj).Contains("part");
    if (isMCP)
      h = h->Rebin(nptBinsGen, Form("%s_rb", tag), ptbinGen);
    else
      h = h->Rebin(nptBins, Form("%s_rb", tag), ptbin);
  }
  if (Nevts > 0) h->Scale(1.0 / Nevts, "width");
  return h;
}

//=============================================================================
// Normalize h2 to match h1 at a given pT point
//=============================================================================
Double_t NormalizeAtPt(TH1 *h1, TH1 *h2, Double_t pt) {
  Double_t v1 = h1->GetBinContent(h1->FindBin(pt));
  Double_t v2 = h2->GetBinContent(h2->FindBin(pt));
  if (v2 == 0) return 1.0;
  Double_t scale = v1 / v2;
  h2->Scale(scale);
  return scale;
}

//=============================================================================
// Make ratio histogram
//=============================================================================
TH1 *MakeRatio(TH1 *hNum, TH1 *hDenom, const char *name) {
  TH1 *r = (TH1 *)hNum->Clone(name);
  r->Reset();
  for (int i = 1; i <= hNum->GetNbinsX(); i++) {
    double yN = hNum->GetBinContent(i);
    double eN = hNum->GetBinError(i);
    int jD = hDenom->FindBin(hNum->GetBinCenter(i));
    double yD = hDenom->GetBinContent(jD);
    if (yD > 0 && yN > 0) {
      r->SetBinContent(i, yN / yD);
      r->SetBinError(i, (yN / yD) * (eN / yN));
    }
  }
  return r;
}

//=============================================================================
// Draw comparison panel: MB vs JJ with ratio
//=============================================================================
void DrawComparisonPanel(TH1 *hMB, TH1 *hJJ, const char *title,
                         const char *xtit, const char *ytit,
                         double xmin, double xmax,
                         double ratioMin, double ratioMax,
                         const char *saveName, bool logy = true,
                         bool normalizeAtRefPt = true) {
  TH1::AddDirectory(kFALSE);
  Filipad2 *pad = new Filipad2(title, ++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
  pad->Draw();

  TPad *upper = pad->GetPad(1);
  optFili(*upper, 1, 1, 0, logy ? 1 : 0);
  upper->cd();

  TH1 *hMBd = (TH1 *)hMB->Clone(Form("%s_MBd", title));
  TH1 *hJJd = (TH1 *)hJJ->Clone(Form("%s_JJd", title));
  if (normalizeAtRefPt) NormalizeAtPt(hMBd, hJJd, kNormPt);

  StyleHist(hMBd, kBlue + 1, 20);
  StyleHist(hJJd, kRed, 24);
  hset(*hMBd, xtit, ytit, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05);
  hMBd->GetXaxis()->SetRangeUser(xmin, xmax);
  hMBd->Draw("ep");
  hJJd->Draw("epsame");

  TLegend *leg = new TLegend(0.40, 0.70, 0.92, 0.94);
  leg->SetTextSize(0.04); leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->AddEntry(hMBd, kMBName, "pe");
  leg->AddEntry(hJJd, kJJName, "pe");
  if (normalizeAtRefPt)
    leg->AddEntry((TObject *)0, Form("Normalized at %.0f GeV", kNormPt), "");
  leg->Draw();

  TPad *lower = pad->GetPad(2);
  optFili(*lower, 1, 1, 0, 0);
  lower->cd();

  TH1 *hR = MakeRatio(hJJd, hMBd, Form("r_%s", title));
  StyleHist(hR, kRed, 24);
  hset(*hR, xtit, "JJ / MB", 0.9, 0.5, 0.1, 0.1, 0.01, 0.01, 0.08, 0.08);
  hR->GetXaxis()->SetRangeUser(xmin, xmax);
  hR->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  hR->Draw("ep");
  TLine *l = new TLine(xmin, 1.0, xmax, 1.0);
  l->SetLineStyle(2); l->Draw();

  if (DRAWPLOTS) pad->C->Print(Form("%s/%s", MakeDirName.Data(), saveName));
}

//=============================================================================
// Draw pTHat distribution
//=============================================================================
void DrawPtHat(const char *mbPath, const char *jjPath) {
  // Try weighted pTHat histograms
  TH1 *hMB = LoadHist1D(mbPath, kMBDir, PtHatCollObj, "pthat_MB");
  TH1 *hJJ = LoadHist1D(jjPath, kJJDir, PtHatCollObj, "pthat_JJ");
  // Fallback to MC collision pTHat
  if (!hMB) hMB = LoadHist1D(mbPath, kMBDir, PtHatMCCollObj, "pthat_mc_MB");
  if (!hJJ) hJJ = LoadHist1D(jjPath, kJJDir, PtHatMCCollObj, "pthat_mc_JJ");

  if (!hJJ) {
    cout << "[PtHat] No pTHat histogram found in JJ file" << endl;
    return;
  }

  TCanvas *c = new TCanvas("cPtHat", "pTHat after scaling", 800, 600);
  c->SetLogy(1); c->SetGridx(1); c->SetGridy(1);

  StyleHist(hJJ, kRed);
  hset(*hJJ, "#it{p}_{T,hard} (GeV/#it{c})", "Counts (weighted)", 0.9, 1.2);
  hJJ->Draw("hist");

  TLegend *leg = new TLegend(0.55, 0.72, 0.88, 0.88);
  leg->SetTextSize(0.04); leg->SetBorderSize(0);
  leg->AddEntry(hJJ, kJJName, "l");

  if (hMB) {
    StyleHist(hMB, kBlue + 1);
    hMB->Draw("histsame");
    leg->AddEntry(hMB, kMBName, "l");
  }
  leg->Draw();

  if (DRAWPLOTS) c->Print(Form("%s/PtHat_afterScaling.pdf", MakeDirName.Data()));
}

//=============================================================================
// Draw response matrix (2D colz)
//=============================================================================
void DrawResponseMatrix(const char *fileName, const char *dir, const char *name) {
  TH2 *h2 = LoadHist2D(fileName, dir, ResponseObj, Form("RM_%s", name));
  if (!h2) { cout << "[RM] Not found for " << name << endl; return; }

  TCanvas *c = new TCanvas(Form("cRM_%s", name), Form("Response %s", name), 800, 700);
  c->SetLogz(1);
  c->SetLeftMargin(0.12); c->SetRightMargin(0.15); c->SetBottomMargin(0.12);
  h2->GetXaxis()->SetRangeUser(0, 200);
  h2->GetYaxis()->SetRangeUser(0, 200);
  hset(*h2, "#it{p}_{T, jet}^{reco} (GeV/#it{c})", "#it{p}_{T, jet}^{true} (GeV/#it{c})");
  h2->SetTitle(Form("Response Matrix: %s", name));
  h2->Draw("colz");

  if (DRAWPLOTS) c->Print(Form("%s/ResponseMatrix_%s.pdf", MakeDirName.Data(), name));
}

void DrawPurityEffCombined(const char *mbPath, const char *jjPath) {
  // Load for both
  auto loadPE = [](const char *fName, const char *dir, const char *tag,
                   TH1 *&purity, TH1 *&eff) {
    TH2 *h2RM = LoadHist2D(fName, dir, ResponseObj, Form("PE2_RM_%s", tag));
    TH1 *hMCD = LoadHist1D(fName, dir, JetPt1DObj, Form("PE2_MCD_%s", tag));
    TH1 *hMCP = LoadHist1D(fName, dir, JetPtMCPObj, Form("PE2_MCP_%s", tag));
    if (!h2RM || !hMCD || !hMCP) { purity = nullptr; eff = nullptr; return; }
    if (REBINON) {
      hMCD = hMCD->Rebin(nptBins, Form("PE2_MCD_%s_rb", tag), ptbin);
      hMCP = hMCP->Rebin(nptBinsGen, Form("PE2_MCP_%s_rb", tag), ptbinGen);
    }
    TH2F *h2 = new TH2F(Form("PE2_RM2_%s", tag), "", nptBins, ptbin, nptBinsGen, ptbinGen);
    for (Int_t ix = 1; ix <= h2RM->GetNbinsX(); ++ix)
      for (Int_t iy = 1; iy <= h2RM->GetNbinsY(); ++iy) {
        Double_t c = h2RM->GetBinContent(ix, iy);
        Double_t e = h2RM->GetBinError(ix, iy);
        Int_t xb = h2->GetXaxis()->FindBin(h2RM->GetXaxis()->GetBinCenter(ix));
        Int_t yb = h2->GetYaxis()->FindBin(h2RM->GetYaxis()->GetBinCenter(iy));
        h2->SetBinContent(xb, yb, h2->GetBinContent(xb, yb) + c);
        h2->SetBinError(xb, yb, sqrt(h2->GetBinError(xb, yb) * h2->GetBinError(xb, yb) + e * e));
      }
    TH1 *hMCDm = h2->ProjectionX(Form("MCDm2_%s", tag));
    TH1 *hMCPm = h2->ProjectionY(Form("MCPm2_%s", tag));
    purity = (TH1 *)hMCD->Clone(Form("pur2_%s", tag)); purity->Reset();
    for (int i = 1; i <= hMCD->GetNbinsX(); i++) {
      double d = hMCD->GetBinContent(i), n = hMCDm->GetBinContent(i);
      if (d > 0) { double p = n / d; purity->SetBinContent(i, TMath::Max(0.0, TMath::Min(p, 1.0))); }
    }
    eff = (TH1 *)hMCP->Clone(Form("eff2_%s", tag)); eff->Reset();
    for (int i = 1; i <= hMCP->GetNbinsX(); i++) {
      double d = hMCP->GetBinContent(i), n = hMCPm->GetBinContent(i);
      if (d > 0) { double p = n / d; eff->SetBinContent(i, TMath::Max(0.0, TMath::Min(p, 1.0))); }
    }
  };

  TH1 *purMB, *effMB, *purJJ, *effJJ;
  loadPE(mbPath, kMBDir, "MB", purMB, effMB);
  loadPE(jjPath, kJJDir, "JJ", purJJ, effJJ);

  // Purity plot
  if (purMB && purJJ) {
    Filipad2 *p = new Filipad2("Purity", ++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    p->Draw();
    TPad *up = p->GetPad(1); optFili(*up, 1, 1, 0, 0); up->cd();
    StyleHist(purMB, kBlue + 1, 20); StyleHist(purJJ, kRed, 24);
    hset(*purMB, "#it{p}_{T, jet}^{reco} (GeV/#it{c})", "Purity");
    purMB->GetXaxis()->SetRangeUser(PlotPtMin, PlotPtMax);
    purMB->GetYaxis()->SetRangeUser(0, 1.1);
    purMB->Draw("ep"); purJJ->Draw("epsame");
    TLegend *lg = new TLegend(0.15, 0.15, 0.5, 0.35);
    lg->SetTextSize(0.04); lg->SetBorderSize(0);
    lg->AddEntry(purMB, kMBName, "pe"); lg->AddEntry(purJJ, kJJName, "pe");
    lg->Draw();
    TPad *lo = p->GetPad(2); optFili(*lo, 1, 1, 0, 0); lo->cd();
    TH1 *rP = MakeRatio(purJJ, purMB, "rPur");
    StyleHist(rP, kRed, 24);
    hset(*rP, "#it{p}_{T, jet}^{reco} (GeV/#it{c})", "JJ / MB");
    rP->GetXaxis()->SetRangeUser(PlotPtMin, PlotPtMax);
    rP->GetYaxis()->SetRangeUser(0.8, 1.2); rP->Draw("ep");
    TLine *l1 = new TLine(PlotPtMin, 1, PlotPtMax, 1); l1->SetLineStyle(2); l1->Draw();
    if (DRAWPLOTS) p->C->Print(Form("%s/Purity_R%.1f.pdf", MakeDirName.Data(), RBIN));
  }

  // Efficiency plot
  if (effMB && effJJ) {
    Filipad2 *p = new Filipad2("Efficiency", ++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
    p->Draw();
    TPad *up = p->GetPad(1); optFili(*up, 1, 1, 0, 0); up->cd();
    StyleHist(effMB, kBlue + 1, 20); StyleHist(effJJ, kRed, 24);
    hset(*effMB, "#it{p}_{T, jet}^{true} (GeV/#it{c})", "Efficiency");
    effMB->GetXaxis()->SetRangeUser(0, PlotPtMax);
    effMB->GetYaxis()->SetRangeUser(0, 1.1);
    effMB->Draw("ep"); effJJ->Draw("epsame");
    TLegend *lg = new TLegend(0.15, 0.15, 0.5, 0.35);
    lg->SetTextSize(0.04); lg->SetBorderSize(0);
    lg->AddEntry(effMB, kMBName, "pe"); lg->AddEntry(effJJ, kJJName, "pe");
    lg->Draw();
    TPad *lo = p->GetPad(2); optFili(*lo, 1, 1, 0, 0); lo->cd();
    TH1 *rE = MakeRatio(effJJ, effMB, "rEff");
    StyleHist(rE, kRed, 24);
    hset(*rE, "#it{p}_{T, jet}^{true} (GeV/#it{c})", "JJ / MB");
    rE->GetXaxis()->SetRangeUser(0, PlotPtMax);
    rE->GetYaxis()->SetRangeUser(0.8, 1.2); rE->Draw("ep");
    TLine *l1 = new TLine(0, 1, PlotPtMax, 1); l1->SetLineStyle(2); l1->Draw();
    if (DRAWPLOTS) p->C->Print(Form("%s/Efficiency_R%.1f.pdf", MakeDirName.Data(), RBIN));
  }
}

//=============================================================================
// Draw kinematic efficiency
// Definition: ProjectionY(selected reco pT window) / ProjectionY(all reco pT)
// from response matrix (matched jets only)
// i.e., fraction of true-level jets whose matched reco jet falls in analysis window
//=============================================================================
void DrawKinematicEfficiency(const char *fileName, const char *dir, const char *name,
                             TH1 *&outKE) {
  TH2 *h2RM = LoadHist2D(fileName, dir, ResponseObj, Form("KE_RM_%s", name));
  if (!h2RM) { cout << "[KinEff] Response matrix not found for " << name << endl; outKE = nullptr; return; }

  // Rebin response matrix
  TH2F *h2 = new TH2F(Form("KE_RM2_%s", name), "", nptBins, ptbin, nptBinsGen, ptbinGen);
  for (Int_t ix = 1; ix <= h2RM->GetNbinsX(); ++ix)
    for (Int_t iy = 1; iy <= h2RM->GetNbinsY(); ++iy) {
      Double_t c = h2RM->GetBinContent(ix, iy);
      Double_t e = h2RM->GetBinError(ix, iy);
      Int_t xb = h2->GetXaxis()->FindBin(h2RM->GetXaxis()->GetBinCenter(ix));
      Int_t yb = h2->GetYaxis()->FindBin(h2RM->GetYaxis()->GetBinCenter(iy));
      h2->SetBinContent(xb, yb, h2->GetBinContent(xb, yb) + c);
      h2->SetBinError(xb, yb, sqrt(h2->GetBinError(xb, yb) * h2->GetBinError(xb, yb) + e * e));
    }

  // Denominator: all matched true-level jets (project Y over full reco range)
  TH1 *hTotal = h2->ProjectionY(Form("KE_total_%s", name), 0, -1, "e");
  hTotal->SetDirectory(0);

  // Numerator: true-level jets matched to reco jets in analysis window [PlotPtMin, PlotPtMax]
  int binLow  = h2->GetXaxis()->FindBin(PlotPtMin + 1e-5);
  int binHigh = h2->GetXaxis()->FindBin(PlotPtMax - 1e-5);
  TH1 *hSelected = h2->ProjectionY(Form("KE_sel_%s", name), binLow, binHigh, "e");
  hSelected->SetDirectory(0);

  // Kinematic efficiency = selected / total
  outKE = (TH1 *)hSelected->Clone(Form("KE_%s", name));
  outKE->Divide(hSelected, hTotal, 1., 1., "B");
  outKE->SetDirectory(0);
}

void DrawKinEffCombined(const char *mbPath, const char *jjPath) {
  TH1 *keMB = nullptr, *keJJ = nullptr;
  DrawKinematicEfficiency(mbPath, kMBDir, "MB", keMB);
  DrawKinematicEfficiency(jjPath, kJJDir, "JJ", keJJ);

  if (!keMB || !keJJ) { cout << "[KinEff] Missing histograms" << endl; return; }

  Filipad2 *p = new Filipad2("KinEff", ++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
  p->Draw();

  TPad *up = p->GetPad(1); optFili(*up, 1, 1, 0, 0); up->cd();
  StyleHist(keMB, kBlue + 1, 20); StyleHist(keJJ, kRed, 24);
  hset(*keMB, "#it{p}_{T, jet}^{true} (GeV/#it{c})", "Kinematic efficiency");
  keMB->GetXaxis()->SetRangeUser(0, PlotPtMax);
  keMB->GetYaxis()->SetRangeUser(0, 1.1);
  keMB->Draw("ep"); keJJ->Draw("epsame");

  // Vertical lines at analysis pT window
  TLine *lL = new TLine(PlotPtMin, 0, PlotPtMin, 1.1);
  lL->SetLineStyle(2); lL->SetLineColor(kGray + 3); lL->SetLineWidth(2); lL->Draw();
  TLine *lR = new TLine(140, 0, 140, 1.1);
  lR->SetLineStyle(2); lR->SetLineColor(kGray + 3); lR->SetLineWidth(2); lR->Draw();

  TLegend *lg = new TLegend(0.15, 0.15, 0.55, 0.35);
  lg->SetTextSize(0.04); lg->SetBorderSize(0);
  lg->AddEntry(keMB, kMBName, "pe"); lg->AddEntry(keJJ, kJJName, "pe");
  lg->AddEntry(lL, Form("Analysis window [%.0f, 140] GeV", PlotPtMin), "l");
  lg->Draw();

  TPad *lo = p->GetPad(2); optFili(*lo, 1, 1, 0, 0); lo->cd();
  TH1 *rKE = MakeRatio(keJJ, keMB, "rKE");
  StyleHist(rKE, kRed, 24);
  hset(*rKE, "#it{p}_{T, jet}^{true} (GeV/#it{c})", "JJ / MB");
  rKE->GetXaxis()->SetRangeUser(0, PlotPtMax);
  rKE->GetYaxis()->SetRangeUser(0.8, 1.2); rKE->Draw("ep");
  TLine *l1 = new TLine(0, 1, PlotPtMax, 1); l1->SetLineStyle(2); l1->Draw();

  if (DRAWPLOTS) p->C->Print(Form("%s/KinematicEfficiency_R%.1f.pdf", MakeDirName.Data(), RBIN));
}

//=============================================================================
// Draw jet pT resolution in pT classes
//=============================================================================
void DrawJetResolution(const char *fileName, const char *dir, const char *name) {
  TH2 *h2 = LoadHist2D(fileName, dir, ResolutionObj, Form("JR_%s", name));
  if (!h2) { cout << "[JetRes] Not found for " << name << endl; return; }

  TCanvas *c = new TCanvas(Form("cJR_%s", name), Form("Resolution %s", name), 900, 600);
  c->SetLogy(1); c->SetGridx(1); c->SetGridy(1);

  TLegend *leg = new TLegend(0.55, 0.60, 0.88, 0.88);
  leg->SetTextSize(0.035); leg->SetBorderSize(0);
  leg->SetHeader(name);

  for (int ipt = 0; ipt < kNResPtBins; ipt++) {
    int bin1 = h2->GetXaxis()->FindBin(kResPtEdges[ipt]);
    int bin2 = h2->GetXaxis()->FindBin(kResPtEdges[ipt + 1]) - 1;
    TH1 *hSlice = h2->ProjectionY(Form("JR_%s_pt%d", name, ipt), bin1, bin2);
    if (hSlice->Integral() > 0) hSlice->Scale(1.0 / hSlice->Integral());
    StyleHist(hSlice, kResPtColors[ipt], 20, 0.6);
    hset(*hSlice, "(#it{p}_{T}^{true} - #it{p}_{T}^{reco}) / #it{p}_{T}^{true}",
         "Normalized counts");
    hSlice->GetXaxis()->SetRangeUser(-0.5, 0.5);
    if (ipt == 0) hSlice->Draw("hist");
    else hSlice->Draw("histsame");
    leg->AddEntry(hSlice, Form("%.0f < #it{p}_{T}^{true} < %.0f GeV",
                                kResPtEdges[ipt], kResPtEdges[ipt + 1]), "l");
  }
  leg->Draw();

  if (DRAWPLOTS) c->Print(Form("%s/JetResolution_%s.pdf", MakeDirName.Data(), name));
}

//=============================================================================
// Build tracking efficiency from track-efficiency task 3D histograms
// Efficiency = matched primary tracks / truth particles, |eta| < 0.9
//=============================================================================
TH1 *BuildTrackingEff(const char *fileName, const char *tag) {
  TFile *file = TFile::Open(fileName, "read");
  if (!file || file->IsZombie()) return nullptr;

  TH3 *h3Truth = (TH3 *)file->Get(Form("%s/%s", kTrkEffDir, kHistTruthLo));
  TH3 *h3Match = (TH3 *)file->Get(Form("%s/%s", kTrkEffDir, kHistMatchLo));
  if (!h3Truth || !h3Match) { file->Close(); return nullptr; }

  // |eta| < 0.9 cut
  int etaLo = h3Truth->GetYaxis()->FindBin(-0.899);
  int etaHi = h3Truth->GetYaxis()->FindBin(0.899);
  TH1 *hTruth = h3Truth->ProjectionX(Form("trkeff_truth_%s", tag), etaLo, etaHi, 0, -1);
  TH1 *hMatch = h3Match->ProjectionX(Form("trkeff_match_%s", tag), etaLo, etaHi, 0, -1);
  hTruth->SetDirectory(0); hMatch->SetDirectory(0);

  // Try high-pT histograms for extended range
  TH3 *h3TruthHi = (TH3 *)file->Get(Form("%s/%s", kTrkEffDir, kHistTruthHi));
  TH3 *h3MatchHi = (TH3 *)file->Get(Form("%s/%s", kTrkEffDir, kHistMatchHi));

  TH1 *hEff = (TH1 *)hMatch->Clone(Form("trkeff_%s", tag));
  hEff->Divide(hMatch, hTruth, 1., 1., "B");
  hEff->SetDirectory(0);

  // If high-pT available, build combined
  if (h3TruthHi && h3MatchHi) {
    int etaLoH = h3TruthHi->GetYaxis()->FindBin(-0.899);
    int etaHiH = h3TruthHi->GetYaxis()->FindBin(0.899);
    TH1 *hTruthHi = h3TruthHi->ProjectionX(Form("trkeff_truthhi_%s", tag), etaLoH, etaHiH, 0, -1);
    TH1 *hMatchHi = h3MatchHi->ProjectionX(Form("trkeff_matchhi_%s", tag), etaLoH, etaHiH, 0, -1);
    hTruthHi->SetDirectory(0); hMatchHi->SetDirectory(0);
    TH1 *hEffHi = (TH1 *)hMatchHi->Clone(Form("trkeffhi_%s", tag));
    hEffHi->Divide(hMatchHi, hTruthHi, 1., 1., "B");
    hEffHi->SetDirectory(0);

    // Merge: use low-pT bins up to 10 GeV, then high-pT bins
    const int nLo = 200, nHi = 18, nComb = nLo + nHi;
    double xbins[nComb + 1];
    for (int i = 0; i <= nLo; i++) xbins[i] = 0.0 + i * 0.05;
    for (int i = 1; i <= nHi; i++) xbins[nLo + i] = 10.0 + i * 5.0;
    TH1D *hComb = new TH1D(Form("trkeff_comb_%s", tag), "", nComb, xbins);
    hComb->SetDirectory(0);
    for (int ib = 1; ib <= hEff->GetNbinsX(); ib++) {
      int bc = hComb->FindBin(hEff->GetBinCenter(ib));
      hComb->SetBinContent(bc, hEff->GetBinContent(ib));
      hComb->SetBinError(bc, hEff->GetBinError(ib));
    }
    for (int ib = 1; ib <= hEffHi->GetNbinsX(); ib++) {
      int bc = hComb->FindBin(hEffHi->GetBinCenter(ib));
      hComb->SetBinContent(bc, hEffHi->GetBinContent(ib));
      hComb->SetBinError(bc, hEffHi->GetBinError(ib));
    }
    file->Close();
    return hComb;
  }

  file->Close();
  return hEff;
}

void DrawTrackingEffCombined(const char *mbPath, const char *jjPath) {
  TH1 *effMB = BuildTrackingEff(mbPath, "MB");
  TH1 *effJJ = BuildTrackingEff(jjPath, "JJ");
  if (!effMB && !effJJ) {
    cout << "[TrackingEff] track-efficiency histograms not found (task may not be included)" << endl;
    return;
  }

  Filipad2 *p = new Filipad2("TrackingEff", ++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
  p->Draw();
  TPad *up = p->GetPad(1); optFili(*up, 1, 1, 1, 0); up->cd();

  TH1 *hFirst = effMB ? effMB : effJJ;
  StyleHist(hFirst, effMB ? kBlue + 1 : kRed, effMB ? 20 : 24);
  hset(*hFirst, "#it{p}_{T}^{truth} (GeV/#it{c})", "Tracking efficiency");
  hFirst->GetXaxis()->SetRangeUser(0.15, 200);
  hFirst->GetYaxis()->SetRangeUser(0.0, 1.1);
  hFirst->Draw("ep");

  TLegend *lg = new TLegend(0.15, 0.15, 0.55, 0.35);
  lg->SetTextSize(0.04); lg->SetBorderSize(0);
  if (effMB) lg->AddEntry(effMB, kMBName, "pe");

  if (effJJ && effMB) {
    StyleHist(effJJ, kRed, 24);
    effJJ->Draw("epsame");
    lg->AddEntry(effJJ, kJJName, "pe");
  }
  lg->Draw();

  TPad *lo = p->GetPad(2); optFili(*lo, 1, 1, 1, 0); lo->cd();
  if (effMB && effJJ) {
    TH1 *rTE = MakeRatio(effJJ, effMB, "rTrkEff");
    StyleHist(rTE, kRed, 24);
    hset(*rTE, "#it{p}_{T}^{truth} (GeV/#it{c})", "JJ / MB");
    rTE->GetXaxis()->SetRangeUser(0.15, 200);
    rTE->GetYaxis()->SetRangeUser(0.9, 1.1); rTE->Draw("ep");
    TLine *l1 = new TLine(0.15, 1, 200, 1); l1->SetLineStyle(2); l1->Draw();
  }

  if (DRAWPLOTS) p->C->Print(Form("%s/TrackingEfficiency.pdf", MakeDirName.Data()));
}

//=============================================================================
// Draw track pT resolution from track-efficiency task
// h2_particle_pt_track_pt_deltaptoverparticlept: pT_truth vs (pT_reco-pT_truth)/pT_truth
//=============================================================================
void DrawTrackPtResolution(const char *fileName, const char *name) {
  TFile *file = TFile::Open(fileName, "read");
  if (!file || file->IsZombie()) return;
  TH2 *h2 = (TH2 *)file->Get(Form("%s/%s", kTrkEffDir, kHistTrackPtRes));
  if (!h2) {
    cout << "[TrackPtRes] h2_particle_pt_track_pt_deltaptoverparticlept not found for " << name << endl;
    file->Close(); return;
  }
  h2 = (TH2 *)h2->Clone(Form("trkres_%s", name));
  h2->SetDirectory(0);
  file->Close();

  // Profile: mean resolution vs pT
  TProfile *prof = h2->ProfileX(Form("trkres_prof_%s", name));
  prof->SetDirectory(0);

  TCanvas *c = new TCanvas(Form("cTrkRes_%s", name), Form("Track pT Resolution %s", name), 800, 600);
  c->SetLogx(); c->SetGridx(); c->SetGridy();
  c->SetLeftMargin(0.12); c->SetBottomMargin(0.12);
  hset(*prof, "#it{p}_{T}^{truth} (GeV/#it{c})",
       "#LT#Delta#it{p}_{T}/#it{p}_{T}^{truth}#GT");
  prof->GetXaxis()->SetRangeUser(0.15, 200.0);
  prof->SetMinimum(-0.05); prof->SetMaximum(0.05);
  prof->SetMarkerStyle(20); prof->SetMarkerSize(0.5);
  prof->SetMarkerColor(kBlue + 1); prof->SetLineColor(kBlue + 1);
  prof->SetStats(0);
  prof->DrawCopy("PE");
  TLine *l0 = new TLine(0.15, 0.0, 200.0, 0.0);
  l0->SetLineStyle(2); l0->SetLineColor(kGray + 1); l0->Draw();
  TLatex latex; latex.SetNDC(); latex.SetTextSize(0.035);
  latex.DrawLatex(0.15, 0.85, name);
  if (DRAWPLOTS) c->Print(Form("%s/TrackPtResolution_%s.pdf", MakeDirName.Data(), name));

  // 2D colz
  TCanvas *c2 = new TCanvas(Form("cTrkRes2D_%s", name), Form("Track pT Res 2D %s", name), 700, 600);
  c2->SetLogx(); c2->SetLogz();
  c2->SetLeftMargin(0.12); c2->SetRightMargin(0.15); c2->SetBottomMargin(0.12);
  hset(*h2, "#it{p}_{T}^{truth} (GeV/#it{c})", "#Delta#it{p}_{T}/#it{p}_{T}^{truth}");
  h2->GetXaxis()->SetRangeUser(0.15, 200.0);
  h2->GetYaxis()->SetRangeUser(-0.3, 0.3);
  h2->SetStats(0);
  h2->SetTitle(Form("Track #it{p}_{T} resolution: %s", name));
  h2->Draw("colz");
  if (DRAWPLOTS) c2->Print(Form("%s/TrackPtResolution2D_%s.pdf", MakeDirName.Data(), name));
}

//=============================================================================
// Draw unfolded spectrum comparison (MB MC vs JJ MC unfolding of data)
//=============================================================================
void DrawUnfoldedComparison(const char *mbPath, const char *jjPath,
                            const char *dataPath) {
  // Helper lambda for kernel-style unfolding
  auto doUnfold = [](const char *mcFile, const char *mcDir, const char *dataFile,
                     const char *dataDir, const char *tag, Double_t nEvtData) -> TH1 * {
    TH2 *h2RM = LoadHist2D(mcFile, mcDir, ResponseObj, Form("UF_RM_%s", tag));
    TH1 *hMCD = LoadHist1D(mcFile, mcDir, JetPt1DObj, Form("UF_MCD_%s", tag));
    TH1 *hMCP = LoadHist1D(mcFile, mcDir, JetPtMCPObj, Form("UF_MCP_%s", tag));
    TH1 *hData = LoadHist1D(dataFile, dataDir, JetPt1DObj, Form("UF_Data_%s", tag));
    if (!h2RM || !hMCD || !hMCP || !hData) return nullptr;

    if (REBINON) {
      hMCD = hMCD->Rebin(nptBins, Form("UF_MCD_%s_rb", tag), ptbin);
      hMCP = hMCP->Rebin(nptBinsGen, Form("UF_MCP_%s_rb", tag), ptbinGen);
      hData = hData->Rebin(nptBins, Form("UF_Data_%s_rb", tag), ptbin);
    }

    // Build rebinned response matrix
    TH2F *h2 = new TH2F(Form("UF_RM2_%s", tag), "", nptBins, ptbin, nptBinsGen, ptbinGen);
    for (Int_t ix = 1; ix <= h2RM->GetNbinsX(); ++ix)
      for (Int_t iy = 1; iy <= h2RM->GetNbinsY(); ++iy) {
        Double_t c = h2RM->GetBinContent(ix, iy);
        Double_t e = h2RM->GetBinError(ix, iy);
        Int_t xb = h2->GetXaxis()->FindBin(h2RM->GetXaxis()->GetBinCenter(ix));
        Int_t yb = h2->GetYaxis()->FindBin(h2RM->GetYaxis()->GetBinCenter(iy));
        h2->SetBinContent(xb, yb, h2->GetBinContent(xb, yb) + c);
        h2->SetBinError(xb, yb, sqrt(h2->GetBinError(xb, yb) * h2->GetBinError(xb, yb) + e * e));
      }

    TH1 *hMCDm = h2->ProjectionX(Form("UFm_MCD_%s", tag));
    TH1 *hMCPm = h2->ProjectionY(Form("UFm_MCP_%s", tag));

    // Purity
    TH1 *hPurity = (TH1 *)hMCD->Clone(Form("UFpur_%s", tag));
    hPurity->Reset();
    for (int i = 1; i <= hMCD->GetNbinsX(); i++) {
      double d = hMCD->GetBinContent(i), n = hMCDm->GetBinContent(i);
      double p = (d > 0) ? n / d : 0;
      hPurity->SetBinContent(i, TMath::Max(0.0, TMath::Min(p, 1.0)));
    }

    // Efficiency
    TH1 *hEff = (TH1 *)hMCP->Clone(Form("UFeff_%s", tag));
    hEff->Reset();
    for (int i = 1; i <= hMCP->GetNbinsX(); i++) {
      double d = hMCP->GetBinContent(i), n = hMCPm->GetBinContent(i);
      double p = (d > 0) ? n / d : 0;
      hEff->SetBinContent(i, TMath::Max(0.0, TMath::Min(p, 1.0)));
    }

    // Purity-correct data
    TH1 *hDataPur = (TH1 *)hData->Clone(Form("UFdp_%s", tag));
    for (int i = 1; i <= hDataPur->GetNbinsX(); i++) {
      hDataPur->SetBinContent(i, hData->GetBinContent(i) * hPurity->GetBinContent(i));
      hDataPur->SetBinError(i, hData->GetBinError(i) * hPurity->GetBinContent(i));
    }

    // SVD unfold
    RooUnfoldResponse *resp = new RooUnfoldResponse(hMCDm, hMCPm, h2);
    RooUnfoldSvd unfold(resp, hDataPur, kSVDkReg);
    TH1 *hUnfolded = (TH1 *)unfold.Hreco();
    hUnfolded->SetDirectory(0);

    // Efficiency correction
    for (int i = 1; i <= hUnfolded->GetNbinsX(); i++) {
      double eff = hEff->GetBinContent(i);
      if (eff > 0) {
        hUnfolded->SetBinContent(i, hUnfolded->GetBinContent(i) / eff);
        hUnfolded->SetBinError(i, hUnfolded->GetBinError(i) / eff);
      }
    }

    // Normalize: 1/Nevt dN/dpT/deta
    double deltaEta = 2.0 * (0.9 - RBIN);
    hUnfolded->Scale(1.0 / nEvtData / deltaEta, "width");

    return hUnfolded;
  };

  Double_t nEvtData = GetNevents(dataPath, kDataDir, 0);

  TH1 *hUF_MB = doUnfold(mbPath, kMBDir, dataPath, kDataDir, "MB", nEvtData);
  TH1 *hUF_JJ = doUnfold(jjPath, kJJDir, dataPath, kDataDir, "JJ", nEvtData);

  if (!hUF_MB || !hUF_JJ) {
    cout << "[Unfold] Failed to produce unfolded spectra" << endl;
    return;
  }

  // Draw comparison
  Filipad2 *pad = new Filipad2("Unfolded", ++nn, 2, 0.3, 100, 50, 0.7, 1, 1);
  pad->Draw();

  TPad *upper = pad->GetPad(1);
  optFili(*upper, 0, 0, 0, 1); upper->cd();

  StyleHist(hUF_MB, kBlue + 1, 20);
  StyleHist(hUF_JJ, kRed, 24);
  hset(*hUF_MB, "#it{p}_{T, jet}^{ch} (GeV/#it{c})",
       "1/#it{N}_{evt} d^{2}#it{N}/d#it{p}_{T}d#it{#eta}");
  hUF_MB->GetXaxis()->SetRangeUser(PlotPtMin, 140);
  hUF_MB->Draw("ep");
  hUF_JJ->Draw("epsame");

  TLegend *leg = new TLegend(0.45, 0.70, 0.92, 0.92);
  leg->SetTextSize(0.04); leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->AddEntry(hUF_MB, Form("Unfolded w/ %s", kMBName.Data()), "pe");
  leg->AddEntry(hUF_JJ, Form("Unfolded w/ %s", kJJName.Data()), "pe");
  leg->Draw();

  TLegend *leg2 = new TLegend(0.15, 0.04, 0.45, 0.22);
  leg2->SetTextSize(0.04); leg2->SetBorderSize(0); leg2->SetFillColorAlpha(0, 0);
  leg2->AddEntry("", "ALICE WIP", "");
  leg2->AddEntry("", "pp #sqrt{#it{s}} = 13.6 TeV", "");
  leg2->AddEntry("", Form("Anti-#it{k}_{T}, #it{R} = %.1f", RBIN), "");
  leg2->Draw();

  TPad *lower = pad->GetPad(2);
  optFili(*lower, 1, 1, 0, 0); lower->cd();

  TH1 *hR = MakeRatio(hUF_JJ, hUF_MB, "rUF");
  StyleHist(hR, kRed, 24);
  hset(*hR, "#it{p}_{T, jet}^{ch} (GeV/#it{c})", "JJ / MB unfolded", 0.9, 0.5, 0.1, 0.1);
  hR->GetXaxis()->SetRangeUser(PlotPtMin, 140);
  hR->GetYaxis()->SetRangeUser(0.8, 1.2);
  hR->Draw("ep");
  TLine *l = new TLine(PlotPtMin, 1, 140, 1); l->SetLineStyle(2); l->Draw();

  if (DRAWPLOTS) pad->C->Print(Form("%s/Unfolded_R%.1f_MBvsJJ.pdf", MakeDirName.Data(), RBIN));
}

//=============================================================================
// Main
//=============================================================================
void DrawJJMCQA() {
  TH1::AddDirectory(kFALSE);
  gSystem->MakeDirectory(MakeDirName.Data());

  TString mbPath = mainDir + kMBFile;
  TString jjPath = mainDir + kJJFile;
  TString dataPath = mainDir + kDataFile;

  Double_t nEvtMB_MCD = NORMEVENTS ? GetNevents(mbPath, kMBDir, 0) : 1.0;
  Double_t nEvtMB_MCP = NORMEVENTS ? GetNevents(mbPath, kMBDir, 1) : 1.0;
  Double_t nEvtJJ_MCD = NORMEVENTS ? GetNevents(jjPath, kJJDir, 0) : 1.0;
  Double_t nEvtJJ_MCP = NORMEVENTS ? GetNevents(jjPath, kJJDir, 1) : 1.0;

  cout << "=== Event counts ===" << endl;
  cout << "MB MCD: " << nEvtMB_MCD << ", MCP: " << nEvtMB_MCP << endl;
  cout << "JJ MCD: " << nEvtJJ_MCD << ", MCP: " << nEvtJJ_MCP << endl;

  // === TRACK LEVEL ===
  if (kPlotTrackPt) {
    TH1 *hMB = LoadHist1D(mbPath, kMBDir, TrackPtObj, "trkpt_MB");
    TH1 *hJJ = LoadHist1D(jjPath, kJJDir, TrackPtObj, "trkpt_JJ");
    if (hMB && hJJ) {
      if (REBINON) {
        hMB = hMB->Rebin(nTrackptbin, "trkpt_MB_rb", Trackptbin);
        hJJ = hJJ->Rebin(nTrackptbin, "trkpt_JJ_rb", Trackptbin);
      }
      if (NORMEVENTS) { hMB->Scale(1.0 / nEvtMB_MCD, "width"); hJJ->Scale(1.0 / nEvtJJ_MCD, "width"); }
      DrawComparisonPanel(hMB, hJJ, "TrackPt", "#it{p}_{T, track} (GeV/#it{c})",
                          "1/#it{N}_{evt} d#it{N}/d#it{p}_{T}", 0.15, 200, 0.5, 1.5, "TrackPt.pdf");
    } else cout << "[TrackPt] Histograms not found (processTracksQC may be removed)" << endl;
  }

  if (kPlotTrackEtaPhi) {
    TH2 *h2MB = LoadHist2D(mbPath, kMBDir, TrackEtaPhiObj, "trketaphi_MB");
    TH2 *h2JJ = LoadHist2D(jjPath, kJJDir, TrackEtaPhiObj, "trketaphi_JJ");
    if (h2MB && h2JJ) {
      // Eta
      TH1 *hEtaMB = h2MB->ProjectionX("trketa_MB"); TH1 *hEtaJJ = h2JJ->ProjectionX("trketa_JJ");
      if (NORMEVENTS) { hEtaMB->Scale(1.0 / nEvtMB_MCD, "width"); hEtaJJ->Scale(1.0 / nEvtJJ_MCD, "width"); }
      DrawComparisonPanel(hEtaMB, hEtaJJ, "TrackEta", "#it{#eta}_{track}",
                          "1/#it{N}_{evt} d#it{N}/d#it{#eta}", -1, 1, 0.8, 1.2, "TrackEta.pdf", false);
      // Phi
      TH1 *hPhiMB = h2MB->ProjectionY("trkphi_MB"); TH1 *hPhiJJ = h2JJ->ProjectionY("trkphi_JJ");
      if (NORMEVENTS) { hPhiMB->Scale(1.0 / nEvtMB_MCD, "width"); hPhiJJ->Scale(1.0 / nEvtJJ_MCD, "width"); }
      DrawComparisonPanel(hPhiMB, hPhiJJ, "TrackPhi", "#it{#varphi}_{track}",
                          "1/#it{N}_{evt} d#it{N}/d#it{#varphi}", 0, 6.3, 0.8, 1.2, "TrackPhi.pdf", false);
    } else cout << "[TrackEtaPhi] Histograms not found" << endl;
  }

  if (kPlotTrackingEff) DrawTrackingEffCombined(mbPath, jjPath);

  if (kPlotTrackPtRes) {
    DrawTrackPtResolution(mbPath, "MB_MC");
    DrawTrackPtResolution(jjPath, "JJ_MC");
  }

  // === JET LEVEL ===
  if (kPlotJetPtMCD) {
    TH1 *hMB = LoadJetPt(mbPath, kMBDir, JetPt3DObj, nEvtMB_MCD, true, "MB_MCD");
    TH1 *hJJ = LoadJetPt(jjPath, kJJDir, JetPt3DObj, nEvtJJ_MCD, true, "JJ_MCD");
    if (hMB && hJJ)
      DrawComparisonPanel(hMB, hJJ, "JetPtMCD", "#it{p}_{T, jet}^{reco} (GeV/#it{c})",
                          "1/#it{N}_{evt} d#it{N}/d#it{p}_{T}", PlotPtMin, PlotPtMax, 0.5, 1.5,
                          Form("JetPt_R%.1f_rebin.pdf", RBIN));
  }

  if (kPlotJetPtMCD && kPlotJetPtFineBin) {
    TH1 *hMB = LoadJetPt(mbPath, kMBDir, JetPt3DObj, nEvtMB_MCD, false, "MB_MCD_f");
    TH1 *hJJ = LoadJetPt(jjPath, kJJDir, JetPt3DObj, nEvtJJ_MCD, false, "JJ_MCD_f");
    if (hMB && hJJ)
      DrawComparisonPanel(hMB, hJJ, "JetPtMCD_fine", "#it{p}_{T, jet}^{reco} (GeV/#it{c})",
                          "1/#it{N}_{evt} d#it{N}/d#it{p}_{T}", PlotPtMin, PlotPtMax, 0.0, 2.0,
                          Form("JetPt_R%.1f_finebin.pdf", RBIN));
  }

  if (kPlotJetPtMCP) {
    TH1 *hMB = LoadJetPt(mbPath, kMBDir, JetPtMCPObj, nEvtMB_MCP, true, "MB_MCP");
    TH1 *hJJ = LoadJetPt(jjPath, kJJDir, JetPtMCPObj, nEvtJJ_MCP, true, "JJ_MCP");
    if (hMB && hJJ)
      DrawComparisonPanel(hMB, hJJ, "JetPtMCP", "#it{p}_{T, jet}^{true} (GeV/#it{c})",
                          "1/#it{N}_{evt} d#it{N}/d#it{p}_{T}", PlotPtMin, PlotPtMax, 0.5, 1.5,
                          Form("JetPtMCP_R%.1f.pdf", RBIN));
  }

  if (kPlotJetPtMCP && kPlotJetPtFineBin) {
    TH1 *hMB = LoadJetPt(mbPath, kMBDir, JetPtMCPObj, nEvtMB_MCP, false, "MB_MCP_f");
    TH1 *hJJ = LoadJetPt(jjPath, kJJDir, JetPtMCPObj, nEvtJJ_MCP, false, "JJ_MCP_f");
    if (hMB && hJJ)
      DrawComparisonPanel(hMB, hJJ, "JetPtMCP_fine", "#it{p}_{T, jet}^{true} (GeV/#it{c})",
                          "1/#it{N}_{evt} d#it{N}/d#it{p}_{T}", PlotPtMin, PlotPtMax, 0.0, 2.0,
                          Form("JetPtMCP_R%.1f_finebin.pdf", RBIN));
  }

  if (kPlotConstituentPt) {
    // MCD constituent pT
    TH2 *h2MB = LoadHist2D(mbPath, kMBDir, ConstPtMCDObj, "cpt_MCD_MB");
    TH2 *h2JJ = LoadHist2D(jjPath, kJJDir, ConstPtMCDObj, "cpt_MCD_JJ");
    if (h2MB && h2JJ) {
      TH1 *hMB = h2MB->ProjectionY("cpt_MCD_MB_py"); TH1 *hJJ = h2JJ->ProjectionY("cpt_MCD_JJ_py");
      if (NORMEVENTS) { hMB->Scale(1.0/nEvtMB_MCD, "width"); hJJ->Scale(1.0/nEvtJJ_MCD, "width"); }
      DrawComparisonPanel(hMB, hJJ, "ConstPtMCD", "#it{p}_{T, constituent}^{reco} (GeV/#it{c})",
                          "1/#it{N}_{evt} d#it{N}/d#it{p}_{T}", 0.15, 100, 0.5, 1.5, "ConstituentPt_MCD.pdf");
    }
    // MCP constituent pT
    h2MB = LoadHist2D(mbPath, kMBDir, ConstPtMCPObj, "cpt_MCP_MB");
    h2JJ = LoadHist2D(jjPath, kJJDir, ConstPtMCPObj, "cpt_MCP_JJ");
    if (h2MB && h2JJ) {
      TH1 *hMB = h2MB->ProjectionY("cpt_MCP_MB_py"); TH1 *hJJ = h2JJ->ProjectionY("cpt_MCP_JJ_py");
      if (NORMEVENTS) { hMB->Scale(1.0/nEvtMB_MCP, "width"); hJJ->Scale(1.0/nEvtJJ_MCP, "width"); }
      DrawComparisonPanel(hMB, hJJ, "ConstPtMCP", "#it{p}_{T, constituent}^{true} (GeV/#it{c})",
                          "1/#it{N}_{evt} d#it{N}/d#it{p}_{T}", 0.15, 100, 0.5, 1.5, "ConstituentPt_MCP.pdf");
    }
  }

  if (kPlotJetEtaPhi) {
    TFile *fMB = TFile::Open(mbPath, "read"); TFile *fJJ = TFile::Open(jjPath, "read");
    if (fMB && fJJ) {
      TH3 *h3MB = (TH3 *)fMB->Get(Form("%s/%s", kMBDir.Data(), JetPt3DObj));
      TH3 *h3JJ = (TH3 *)fJJ->Get(Form("%s/%s", kJJDir.Data(), JetPt3DObj));
      if (h3MB && h3JJ) {
        TH1 *hEtaMB = h3MB->ProjectionY("jeta_MB_py"); hEtaMB->SetDirectory(0);
        TH1 *hEtaJJ = h3JJ->ProjectionY("jeta_JJ_py"); hEtaJJ->SetDirectory(0);
        TH1 *hPhiMB = h3MB->ProjectionZ("jphi_MB_pz"); hPhiMB->SetDirectory(0);
        TH1 *hPhiJJ = h3JJ->ProjectionZ("jphi_JJ_pz"); hPhiJJ->SetDirectory(0);
        fMB->Close(); fJJ->Close();

        hEtaMB->Rebin(5); hEtaJJ->Rebin(5);
        if (NORMEVENTS) {
          hEtaMB->Scale(1.0/nEvtMB_MCD, "width"); hEtaJJ->Scale(1.0/nEvtJJ_MCD, "width");
          hPhiMB->Scale(1.0/nEvtMB_MCD, "width"); hPhiJJ->Scale(1.0/nEvtJJ_MCD, "width");
        }
        DrawComparisonPanel(hEtaMB, hEtaJJ, "JetEta", "#it{#eta}_{jet}",
                            "1/#it{N}_{evt} d#it{N}/d#it{#eta}", -1, 1, 0.8, 1.2,
                            Form("JetEta_R%.1f.pdf", RBIN), false);
        hPhiMB->Rebin(4); hPhiJJ->Rebin(4);
        DrawComparisonPanel(hPhiMB, hPhiJJ, "JetPhi", "#it{#varphi}_{jet}",
                            "1/#it{N}_{evt} d#it{N}/d#it{#varphi}", 0, 6.3, 0.8, 1.2,
                            Form("JetPhi_R%.1f.pdf", RBIN), false);
      } else { if (fMB) fMB->Close(); if (fJJ) fJJ->Close(); }
    }
  }

  if (kPlotJetNtracks) {
    TH2 *h2MB = LoadHist2D(mbPath, kMBDir, JetNtracksObj, "jnt_MB");
    TH2 *h2JJ = LoadHist2D(jjPath, kJJDir, JetNtracksObj, "jnt_JJ");
    if (h2MB && h2JJ) {
      TH1 *hMB = h2MB->ProjectionY("jnt_MB_py"); TH1 *hJJ = h2JJ->ProjectionY("jnt_JJ_py");
      if (NORMEVENTS) { hMB->Scale(1.0/nEvtMB_MCD, "width"); hJJ->Scale(1.0/nEvtJJ_MCD, "width"); }
      DrawComparisonPanel(hMB, hJJ, "JetNtracks", "N_{tracks}",
                          "1/#it{N}_{evt} d#it{N}/dN_{trk}", 0, 25, 0.5, 1.5, Form("JetNtracks_R%.1f.pdf", RBIN));
    }
  }

  if (kPlotPtHat) DrawPtHat(mbPath, jjPath);

  if (kPlotResponseMatrix) {
    DrawResponseMatrix(mbPath, kMBDir, "MB_MC");
    DrawResponseMatrix(jjPath, kJJDir, "JJ_MC");
  }

  if (kPlotPurityEff) DrawPurityEffCombined(mbPath, jjPath);

  if (kPlotKinEff) DrawKinEffCombined(mbPath, jjPath);

  if (kPlotJetResolution) {
    DrawJetResolution(mbPath, kMBDir, "MB_MC");
    DrawJetResolution(jjPath, kJJDir, "JJ_MC");
  }

  if (kPlotUnfolded) DrawUnfoldedComparison(mbPath, jjPath, dataPath);

  cout << endl << "=== QA plots saved to: " << MakeDirName << " ===" << endl;
}
