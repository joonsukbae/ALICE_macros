///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Draw Data Cross Section         //////////
////////// Using merged MB+JJ MC RM             //////////
////////// author: Joonsuk Bae           ////////// 
////////// E-mail: jbae@cern.ch          //////////
////////// Last Modified: 2025            //////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

#include "DrawJetsMCRDependentHelpers.h"
#include "BSHelper.cxx"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include <RooUnfoldBayes.h>
#include <RooUnfoldResponse.h>
#include <RooUnfoldSvd.h>
#include <TMatrixD.h>

// ============================================================
// Configuration
// ============================================================

// MB MC files (LHC24f3c) - used for low pT (< 20 GeV)
const char *kResponseFileMB =
    // "/Users/js/cernbox/workspace/O2Physics/AOD_copy_macro/merged_Response_LHC24f3c.root";
    // "~/cernbox/workspace/O2Physics/jets/AnalysisResults/515446_AnalysisResults.root"; // LHC24f3c, qOverPtData=1.0
    "~/cernbox/workspace/O2Physics/jets/AnalysisResults/555891_AnalysisResults.root"; // LHC24f3c, qOverPtData=1.5
    // "~/cernbox/workspace/O2Physics/jets/AnalysisResults/585837_AnalysisResults.root"; // LHC24f3c, file-based pT smearing

// JJ MC files (LHC25a2b) - used for high pT (>= 20 GeV)
const char *kResponseFileJJ =
    // "/Users/js/cernbox/workspace/O2Physics/AOD_copy_macro/merged_Response_LHC25a2b.root";
    "~/cernbox/workspace/O2Physics/jets/AnalysisResults/516969_AnalysisResults.root";

// Data file (Run 3, 2022 pass7)
const char *kDataFile = "~/cernbox/workspace/O2Physics/jets/AnalysisResults/498133_AnalysisResults.root";

// Merge boundary: 20 GeV (bin edge)
const Double_t kMergeBoundary = 20.0;

// Directory and object names
// const char *kDirJets = "jet-spectra-charged_id34413";  // MB MC directory
const char *kDirJets = "jet-spectra-charged";  // MB MC directory
const char *kDirJetsJJ = "jet-spectra-charged_Nmax1p5";  // JJ MC directory
const char *kDirJetsData = "jet-spectra-charged";  // Data directory
const char *kHistTrue = "h_jet_pt_part";
const char *kHistReco = "h_jet_pt";
const char *kHistResp2D = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint";

// SVD unfolding parameter: determined dynamically via d-vector + stability
// (computed per response matrix in UnfoldAndPlotCrossSection)

// Output directory
TString CrossSectionDirName = "plots/AN_Charged-particle-jet-cross-section-in-pp-collisions-at-13.6-TeV/Figures/CrossSection";

// ============================================================
// Helper functions
// ============================================================

// Event count helpers
Double_t GetEventCount(TFile *f, const char *dirName) {
  TDirectory *dir = dynamic_cast<TDirectory *>(f->Get(dirName));
  if (!dir) return 0.;
  TH1 *h = dynamic_cast<TH1 *>(dir->Get("h_collisions"));
  if (!h) return 0.;
  return h->GetBinContent(3);  // bin 3 = sel8
}

Double_t GetEventCountForTruth(TFile *f, const char *dirName, const char *histName = "h_mcColl_counts") {
  TDirectory *dir = dynamic_cast<TDirectory *>(f->Get(dirName));
  if (!dir) return 0.;
  TH1 *h = dynamic_cast<TH1 *>(dir->Get(histName));
  if (!h) return 0.;
  return h->GetBinContent(4);  // bin 4 = selMC
}

Double_t GetEventCountForUnfolded(TFile *f, const char *dirName, const char *histName = "h_collisions") {
  TDirectory *dir = dynamic_cast<TDirectory *>(f->Get(dirName));
  if (!dir) return 0.;
  TH1 *h = dynamic_cast<TH1 *>(dir->Get(histName));
  if (!h) return 0.;
  return h->GetBinContent(3);  // bin 3 = sel8
}

// ============================================================
// Load response sample functions
// ============================================================

// Load MB MC response sample only (from DrawMcClosureTest.C LoadResponseSample)
bool LoadResponseSampleMB(TH1D *&hTrue, TH1D *&hReco, TH2D *&hResp2D, Double_t &nEventsResp) {
  std::cerr << "========================================" << std::endl;
  std::cerr << "[Info] Loading MB MC response sample only" << std::endl;
  std::cerr << "  MB MC: " << kResponseFileMB << std::endl;
  std::cerr << "========================================" << std::endl;
  
  TFile *f = TFile::Open(kResponseFileMB, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "[Error] Cannot open MB MC response file: " << kResponseFileMB << std::endl;
    return false;
  }

  TDirectory *dir = dynamic_cast<TDirectory *>(f->Get(kDirJets));
  if (!dir) {
    std::cerr << "[Error] Directory not found: " << kDirJets << std::endl;
    f->Close();
    delete f;
    return false;
  }

  TObject *objTrue = dir->Get(kHistTrue);
  TObject *objReco = dir->Get(kHistReco);
  TObject *objResp = dir->Get(kHistResp2D);

  TH1D *hTrueIn = dynamic_cast<TH1D *>(objTrue);
  TH1D *hRecoIn = dynamic_cast<TH1D *>(objReco);
  TH2D *hRespIn = dynamic_cast<TH2D *>(objResp);

  // Convert TH1F/TH2F to TH1D/TH2D if needed (same as LoadMergedResponseSample)
  if (!hTrueIn) {
    TH1F *hTrueF = dynamic_cast<TH1F *>(objTrue);
    if (hTrueF) {
      if (hTrueF->GetXaxis()->GetXbins() && hTrueF->GetXaxis()->GetXbins()->GetSize() > 0) {
        hTrueIn = new TH1D("hTrueMB_converted", hTrueF->GetTitle(),
                          hTrueF->GetNbinsX(), hTrueF->GetXaxis()->GetXbins()->GetArray());
      } else {
        hTrueIn = new TH1D("hTrueMB_converted", hTrueF->GetTitle(),
                          hTrueF->GetNbinsX(), hTrueF->GetXaxis()->GetXmin(), 
                          hTrueF->GetXaxis()->GetXmax());
      }
      for (int i = 1; i <= hTrueF->GetNbinsX(); ++i) {
        hTrueIn->SetBinContent(i, hTrueF->GetBinContent(i));
        hTrueIn->SetBinError(i, hTrueF->GetBinError(i));
      }
    }
  }
  if (!hRecoIn) {
    TH1F *hRecoF = dynamic_cast<TH1F *>(objReco);
    if (hRecoF) {
      if (hRecoF->GetXaxis()->GetXbins() && hRecoF->GetXaxis()->GetXbins()->GetSize() > 0) {
        hRecoIn = new TH1D("hRecoMB_converted", hRecoF->GetTitle(),
                          hRecoF->GetNbinsX(), hRecoF->GetXaxis()->GetXbins()->GetArray());
      } else {
        hRecoIn = new TH1D("hRecoMB_converted", hRecoF->GetTitle(),
                          hRecoF->GetNbinsX(), hRecoF->GetXaxis()->GetXmin(), 
                          hRecoF->GetXaxis()->GetXmax());
      }
      for (int i = 1; i <= hRecoF->GetNbinsX(); ++i) {
        hRecoIn->SetBinContent(i, hRecoF->GetBinContent(i));
        hRecoIn->SetBinError(i, hRecoF->GetBinError(i));
      }
    }
  }
  if (!hRespIn) {
    TH2F *hRespF = dynamic_cast<TH2F *>(objResp);
    if (hRespF) {
      const Double_t *xbins = hRespF->GetXaxis()->GetXbins() && hRespF->GetXaxis()->GetXbins()->GetSize() > 0 
                             ? hRespF->GetXaxis()->GetXbins()->GetArray() : nullptr;
      const Double_t *ybins = hRespF->GetYaxis()->GetXbins() && hRespF->GetYaxis()->GetXbins()->GetSize() > 0 
                             ? hRespF->GetYaxis()->GetXbins()->GetArray() : nullptr;
      if (xbins && ybins) {
        hRespIn = new TH2D("hRespMB_converted", hRespF->GetTitle(),
                          hRespF->GetNbinsX(), xbins, hRespF->GetNbinsY(), ybins);
      } else if (xbins) {
        hRespIn = new TH2D("hRespMB_converted", hRespF->GetTitle(),
                          hRespF->GetNbinsX(), xbins, hRespF->GetNbinsY(), 
                          hRespF->GetYaxis()->GetXmin(), hRespF->GetYaxis()->GetXmax());
      } else if (ybins) {
        hRespIn = new TH2D("hRespMB_converted", hRespF->GetTitle(),
                          hRespF->GetNbinsX(), hRespF->GetXaxis()->GetXmin(), 
                          hRespF->GetXaxis()->GetXmax(), hRespF->GetNbinsY(), ybins);
      } else {
        hRespIn = new TH2D("hRespMB_converted", hRespF->GetTitle(),
                          hRespF->GetNbinsX(), hRespF->GetXaxis()->GetXmin(), 
                          hRespF->GetXaxis()->GetXmax(),
                          hRespF->GetNbinsY(), hRespF->GetYaxis()->GetXmin(), 
                          hRespF->GetYaxis()->GetXmax());
      }
      for (int i = 1; i <= hRespF->GetNbinsX(); ++i) {
        for (int j = 1; j <= hRespF->GetNbinsY(); ++j) {
          hRespIn->SetBinContent(i, j, hRespF->GetBinContent(i, j));
          hRespIn->SetBinError(i, j, hRespF->GetBinError(i, j));
        }
      }
    }
  }

  if (!hTrueIn || !hRecoIn || !hRespIn) {
    std::cerr << "[Error] Missing histograms in MB MC response file" << std::endl;
    f->Close();
    delete f;
    return false;
  }

  // Rebin to analysis binning (same as LoadMergedResponseSample)
  const Double_t *ptbin = GetPtbin();
  const Double_t *ptbinGen = GetPtbinGen();
  Int_t nptBins = GetNptBins();
  Int_t nptBinsGen = GetNptBinsGen();
  
  hTrue = dynamic_cast<TH1D *>(hTrueIn->Rebin(nptBinsGen, "hTrueMB_gen", ptbinGen));
  hReco = dynamic_cast<TH1D *>(hRecoIn->Rebin(nptBins, "hRecoMB_det", ptbin));

  // Build rebinned 2D matrix
  hResp2D = new TH2D("hResp2D_MB", "MB Response;#it{p}_{T}^{reco};#it{p}_{T}^{true}",
                     nptBins, ptbin, nptBinsGen, ptbinGen);

  for (int ix = 1; ix <= hRespIn->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hRespIn->GetNbinsY(); ++iy) {
      double content = hRespIn->GetBinContent(ix, iy);
      double error = hRespIn->GetBinError(ix, iy);
      if (content == 0.) continue;
      
      double xCenter = hRespIn->GetXaxis()->GetBinCenter(ix);
      double yCenter = hRespIn->GetYaxis()->GetBinCenter(iy);
      int xbin = hResp2D->GetXaxis()->FindBin(xCenter);
      int ybin = hResp2D->GetYaxis()->FindBin(yCenter);
      
      double currentContent = hResp2D->GetBinContent(xbin, ybin);
      double currentError = hResp2D->GetBinError(xbin, ybin);
      hResp2D->SetBinContent(xbin, ybin, currentContent + content);
      hResp2D->SetBinError(xbin, ybin, std::sqrt(currentError * currentError + error * error));
    }
  }

  hTrue->SetDirectory(nullptr);
  hReco->SetDirectory(nullptr);
  hResp2D->SetDirectory(nullptr);

  nEventsResp = GetEventCountForTruth(f, kDirJets);
  if (nEventsResp <= 0.) {
    nEventsResp = GetEventCount(f, kDirJets);
    if (nEventsResp <= 0.) {
      nEventsResp = hTrue->Integral();
    }
  }
  std::cerr << "[Info] MB MC event count: " << nEventsResp << std::endl;
  std::cerr << "[Info] MB MC histograms loaded (raw counts, no normalization)" << std::endl;
  std::cerr << "  hTrue integral: " << hTrue->Integral() << std::endl;
  std::cerr << "  hReco integral: " << hReco->Integral() << std::endl;
  std::cerr << "  hResp2D integral: " << hResp2D->Integral() << std::endl;
  std::cerr << "[Info] Event count normalization will be done AFTER unfolding for cross section calculation" << std::endl;

  f->Close();
  delete f;
  return true;
}

// Load JJ MC response sample only
bool LoadResponseSampleJJ(TH1D *&hTrue, TH1D *&hReco, TH2D *&hResp2D, Double_t &nEventsResp) {
  std::cerr << "========================================" << std::endl;
  std::cerr << "[Info] Loading JJ MC response sample only" << std::endl;
  std::cerr << "  JJ MC: " << kResponseFileJJ << std::endl;
  std::cerr << "========================================" << std::endl;
  
  TFile *f = TFile::Open(kResponseFileJJ, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "[Error] Cannot open JJ MC response file: " << kResponseFileJJ << std::endl;
    return false;
  }

  TDirectory *dir = dynamic_cast<TDirectory *>(f->Get(kDirJetsJJ));
  if (!dir) {
    std::cerr << "[Error] Directory not found: " << kDirJetsJJ << std::endl;
    f->Close();
    delete f;
    return false;
  }

  TObject *objTrue = dir->Get(kHistTrue);
  TObject *objReco = dir->Get(kHistReco);
  TObject *objResp = dir->Get(kHistResp2D);

  TH1D *hTrueIn = dynamic_cast<TH1D *>(objTrue);
  TH1D *hRecoIn = dynamic_cast<TH1D *>(objReco);
  TH2D *hRespIn = dynamic_cast<TH2D *>(objResp);

  // Convert TH1F/TH2F to TH1D/TH2D if needed (same as LoadResponseSampleMB)
  if (!hTrueIn) {
    TH1F *hTrueF = dynamic_cast<TH1F *>(objTrue);
    if (hTrueF) {
      if (hTrueF->GetXaxis()->GetXbins() && hTrueF->GetXaxis()->GetXbins()->GetSize() > 0) {
        hTrueIn = new TH1D("hTrueJJ_converted", hTrueF->GetTitle(),
                          hTrueF->GetNbinsX(), hTrueF->GetXaxis()->GetXbins()->GetArray());
      } else {
        hTrueIn = new TH1D("hTrueJJ_converted", hTrueF->GetTitle(),
                          hTrueF->GetNbinsX(), hTrueF->GetXaxis()->GetXmin(), 
                          hTrueF->GetXaxis()->GetXmax());
      }
      for (int i = 1; i <= hTrueF->GetNbinsX(); ++i) {
        hTrueIn->SetBinContent(i, hTrueF->GetBinContent(i));
        hTrueIn->SetBinError(i, hTrueF->GetBinError(i));
      }
    }
  }
  if (!hRecoIn) {
    TH1F *hRecoF = dynamic_cast<TH1F *>(objReco);
    if (hRecoF) {
      if (hRecoF->GetXaxis()->GetXbins() && hRecoF->GetXaxis()->GetXbins()->GetSize() > 0) {
        hRecoIn = new TH1D("hRecoJJ_converted", hRecoF->GetTitle(),
                          hRecoF->GetNbinsX(), hRecoF->GetXaxis()->GetXbins()->GetArray());
      } else {
        hRecoIn = new TH1D("hRecoJJ_converted", hRecoF->GetTitle(),
                          hRecoF->GetNbinsX(), hRecoF->GetXaxis()->GetXmin(), 
                          hRecoF->GetXaxis()->GetXmax());
      }
      for (int i = 1; i <= hRecoF->GetNbinsX(); ++i) {
        hRecoIn->SetBinContent(i, hRecoF->GetBinContent(i));
        hRecoIn->SetBinError(i, hRecoF->GetBinError(i));
      }
    }
  }
  if (!hRespIn) {
    TH2F *hRespF = dynamic_cast<TH2F *>(objResp);
    if (hRespF) {
      const Double_t *xbins = hRespF->GetXaxis()->GetXbins() && hRespF->GetXaxis()->GetXbins()->GetSize() > 0 
                             ? hRespF->GetXaxis()->GetXbins()->GetArray() : nullptr;
      const Double_t *ybins = hRespF->GetYaxis()->GetXbins() && hRespF->GetYaxis()->GetXbins()->GetSize() > 0 
                             ? hRespF->GetYaxis()->GetXbins()->GetArray() : nullptr;
      if (xbins && ybins) {
        hRespIn = new TH2D("hRespJJ_converted", hRespF->GetTitle(),
                          hRespF->GetNbinsX(), xbins, hRespF->GetNbinsY(), ybins);
      } else if (xbins) {
        hRespIn = new TH2D("hRespJJ_converted", hRespF->GetTitle(),
                          hRespF->GetNbinsX(), xbins, hRespF->GetNbinsY(), 
                          hRespF->GetYaxis()->GetXmin(), hRespF->GetYaxis()->GetXmax());
      } else if (ybins) {
        hRespIn = new TH2D("hRespJJ_converted", hRespF->GetTitle(),
                          hRespF->GetNbinsX(), hRespF->GetXaxis()->GetXmin(), 
                          hRespF->GetXaxis()->GetXmax(), hRespF->GetNbinsY(), ybins);
      } else {
        hRespIn = new TH2D("hRespJJ_converted", hRespF->GetTitle(),
                          hRespF->GetNbinsX(), hRespF->GetXaxis()->GetXmin(), 
                          hRespF->GetXaxis()->GetXmax(),
                          hRespF->GetNbinsY(), hRespF->GetYaxis()->GetXmin(), 
                          hRespF->GetYaxis()->GetXmax());
      }
      for (int i = 1; i <= hRespF->GetNbinsX(); ++i) {
        for (int j = 1; j <= hRespF->GetNbinsY(); ++j) {
          hRespIn->SetBinContent(i, j, hRespF->GetBinContent(i, j));
          hRespIn->SetBinError(i, j, hRespF->GetBinError(i, j));
        }
      }
    }
  }

  if (!hTrueIn || !hRecoIn || !hRespIn) {
    std::cerr << "[Error] Missing histograms in JJ MC response file" << std::endl;
    f->Close();
    delete f;
    return false;
  }

  // Rebin to analysis binning
  const Double_t *ptbin = GetPtbin();
  const Double_t *ptbinGen = GetPtbinGen();
  Int_t nptBins = GetNptBins();
  Int_t nptBinsGen = GetNptBinsGen();
  
  hTrue = dynamic_cast<TH1D *>(hTrueIn->Rebin(nptBinsGen, "hTrueJJ_gen", ptbinGen));
  hReco = dynamic_cast<TH1D *>(hRecoIn->Rebin(nptBins, "hRecoJJ_det", ptbin));

  // Build rebinned 2D matrix
  hResp2D = new TH2D("hResp2D_JJ", "JJ Response;#it{p}_{T}^{reco};#it{p}_{T}^{true}",
                     nptBins, ptbin, nptBinsGen, ptbinGen);

  for (int ix = 1; ix <= hRespIn->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hRespIn->GetNbinsY(); ++iy) {
      double content = hRespIn->GetBinContent(ix, iy);
      double error = hRespIn->GetBinError(ix, iy);
      if (content == 0.) continue;
      
      double xCenter = hRespIn->GetXaxis()->GetBinCenter(ix);
      double yCenter = hRespIn->GetYaxis()->GetBinCenter(iy);
      int xbin = hResp2D->GetXaxis()->FindBin(xCenter);
      int ybin = hResp2D->GetYaxis()->FindBin(yCenter);
      
      double currentContent = hResp2D->GetBinContent(xbin, ybin);
      double currentError = hResp2D->GetBinError(xbin, ybin);
      hResp2D->SetBinContent(xbin, ybin, currentContent + content);
      hResp2D->SetBinError(xbin, ybin, std::sqrt(currentError * currentError + error * error));
    }
  }

  hTrue->SetDirectory(nullptr);
  hReco->SetDirectory(nullptr);
  hResp2D->SetDirectory(nullptr);

  nEventsResp = GetEventCountForTruth(f, kDirJetsJJ, "h_mcColl_counts_weight");
  if (nEventsResp <= 0.) {
    nEventsResp = GetEventCount(f, kDirJetsJJ);
    if (nEventsResp <= 0.) {
      nEventsResp = hTrue->Integral();
    }
  }
  std::cerr << "[Info] JJ MC event count: " << nEventsResp << std::endl;
  std::cerr << "[Info] JJ MC histograms loaded (raw counts, no normalization)" << std::endl;
  std::cerr << "  hTrue integral: " << hTrue->Integral() << std::endl;
  std::cerr << "  hReco integral: " << hReco->Integral() << std::endl;
  std::cerr << "  hResp2D integral: " << hResp2D->Integral() << std::endl;
  std::cerr << "[Info] Event count normalization will be done AFTER unfolding for cross section calculation" << std::endl;

  f->Close();
  delete f;
  return true;
}

// ============================================================
// Load merged response sample (from DrawMcClosureTest.C)
// ============================================================

bool LoadMergedResponseSample(TH1D *&hTrue, TH1D *&hReco, TH2D *&hResp2D, 
                               Double_t &nEventsMB, Double_t &nEventsJJ) {
  std::cerr << "========================================" << std::endl;
  std::cerr << "[Info] Loading MERGED response sample" << std::endl;
  std::cerr << "  MB MC (low pT < 20 GeV): " << kResponseFileMB << std::endl;
  std::cerr << "  JJ MC (high pT >= 20 GeV): " << kResponseFileJJ << std::endl;
  std::cerr << "========================================" << std::endl;
  
  // Load MB MC
  TFile *fMB = TFile::Open(kResponseFileMB, "READ");
  if (!fMB || fMB->IsZombie()) {
    std::cerr << "[Error] Cannot open MB MC response file" << std::endl;
    return false;
  }
  
  TDirectory *dirMB = dynamic_cast<TDirectory *>(fMB->Get(kDirJets));
  if (!dirMB) {
    std::cerr << "[Error] Directory not found: " << kDirJets << std::endl;
    fMB->Close();
    delete fMB;
    return false;
  }
  
  TObject *objTrueMB = dirMB->Get(kHistTrue);
  TObject *objRecoMB = dirMB->Get(kHistReco);
  TObject *objRespMB = dirMB->Get(kHistResp2D);
  
  TH1D *hTrueInMB = dynamic_cast<TH1D *>(objTrueMB);
  TH1D *hRecoInMB = dynamic_cast<TH1D *>(objRecoMB);
  TH2D *hRespInMB = dynamic_cast<TH2D *>(objRespMB);
  
  // Convert TH1F/TH2F to TH1D/TH2D if needed
  if (!hTrueInMB) {
    TH1F *hTrueFMB = dynamic_cast<TH1F *>(objTrueMB);
    if (hTrueFMB) {
      if (hTrueFMB->GetXaxis()->GetXbins() && hTrueFMB->GetXaxis()->GetXbins()->GetSize() > 0) {
        hTrueInMB = new TH1D("hTrueMB_converted", hTrueFMB->GetTitle(),
                            hTrueFMB->GetNbinsX(), hTrueFMB->GetXaxis()->GetXbins()->GetArray());
      } else {
        hTrueInMB = new TH1D("hTrueMB_converted", hTrueFMB->GetTitle(),
                            hTrueFMB->GetNbinsX(), hTrueFMB->GetXaxis()->GetXmin(), 
                            hTrueFMB->GetXaxis()->GetXmax());
      }
      for (int i = 1; i <= hTrueFMB->GetNbinsX(); ++i) {
        hTrueInMB->SetBinContent(i, hTrueFMB->GetBinContent(i));
        hTrueInMB->SetBinError(i, hTrueFMB->GetBinError(i));
      }
    }
  }
  if (!hRecoInMB) {
    TH1F *hRecoFMB = dynamic_cast<TH1F *>(objRecoMB);
    if (hRecoFMB) {
      if (hRecoFMB->GetXaxis()->GetXbins() && hRecoFMB->GetXaxis()->GetXbins()->GetSize() > 0) {
        hRecoInMB = new TH1D("hRecoMB_converted", hRecoFMB->GetTitle(),
                             hRecoFMB->GetNbinsX(), hRecoFMB->GetXaxis()->GetXbins()->GetArray());
      } else {
        hRecoInMB = new TH1D("hRecoMB_converted", hRecoFMB->GetTitle(),
                             hRecoFMB->GetNbinsX(), hRecoFMB->GetXaxis()->GetXmin(), 
                             hRecoFMB->GetXaxis()->GetXmax());
      }
      for (int i = 1; i <= hRecoFMB->GetNbinsX(); ++i) {
        hRecoInMB->SetBinContent(i, hRecoFMB->GetBinContent(i));
        hRecoInMB->SetBinError(i, hRecoFMB->GetBinError(i));
      }
    }
  }
  if (!hRespInMB) {
    TH2F *hRespFMB = dynamic_cast<TH2F *>(objRespMB);
    if (hRespFMB) {
      const Double_t *xbinsMB = hRespFMB->GetXaxis()->GetXbins() && hRespFMB->GetXaxis()->GetXbins()->GetSize() > 0 
                                ? hRespFMB->GetXaxis()->GetXbins()->GetArray() : nullptr;
      const Double_t *ybinsMB = hRespFMB->GetYaxis()->GetXbins() && hRespFMB->GetYaxis()->GetXbins()->GetSize() > 0 
                                ? hRespFMB->GetYaxis()->GetXbins()->GetArray() : nullptr;
      if (xbinsMB && ybinsMB) {
        hRespInMB = new TH2D("hRespMB_converted", hRespFMB->GetTitle(),
                            hRespFMB->GetNbinsX(), xbinsMB, hRespFMB->GetNbinsY(), ybinsMB);
      } else if (xbinsMB) {
        hRespInMB = new TH2D("hRespMB_converted", hRespFMB->GetTitle(),
                            hRespFMB->GetNbinsX(), xbinsMB, hRespFMB->GetNbinsY(), 
                            hRespFMB->GetYaxis()->GetXmin(), hRespFMB->GetYaxis()->GetXmax());
      } else if (ybinsMB) {
        hRespInMB = new TH2D("hRespMB_converted", hRespFMB->GetTitle(),
                            hRespFMB->GetNbinsX(), hRespFMB->GetXaxis()->GetXmin(), 
                            hRespFMB->GetXaxis()->GetXmax(), hRespFMB->GetNbinsY(), ybinsMB);
      } else {
        hRespInMB = new TH2D("hRespMB_converted", hRespFMB->GetTitle(),
                            hRespFMB->GetNbinsX(), hRespFMB->GetXaxis()->GetXmin(), 
                            hRespFMB->GetXaxis()->GetXmax(),
                            hRespFMB->GetNbinsY(), hRespFMB->GetYaxis()->GetXmin(), 
                            hRespFMB->GetYaxis()->GetXmax());
      }
      for (int i = 1; i <= hRespFMB->GetNbinsX(); ++i) {
        for (int j = 1; j <= hRespFMB->GetNbinsY(); ++j) {
          hRespInMB->SetBinContent(i, j, hRespFMB->GetBinContent(i, j));
          hRespInMB->SetBinError(i, j, hRespFMB->GetBinError(i, j));
        }
      }
    }
  }
  
  if (!hTrueInMB || !hRecoInMB || !hRespInMB) {
    std::cerr << "[Error] Missing histograms in MB MC response file" << std::endl;
    fMB->Close();
    delete fMB;
    return false;
  }
  
  nEventsMB = GetEventCountForTruth(fMB, kDirJets);
  if (nEventsMB <= 0.) {
    nEventsMB = GetEventCount(fMB, kDirJets);
    if (nEventsMB <= 0.) {
      nEventsMB = hTrueInMB->Integral();
    }
  }
  
  // Load JJ MC
  TFile *fJJ = TFile::Open(kResponseFileJJ, "READ");
  if (!fJJ || fJJ->IsZombie()) {
    std::cerr << "[Error] Cannot open JJ MC response file" << std::endl;
    fMB->Close();
    delete fMB;
    return false;
  }
  
  TDirectory *dirJJ = dynamic_cast<TDirectory *>(fJJ->Get(kDirJetsJJ));
  if (!dirJJ) {
    std::cerr << "[Error] Directory not found: " << kDirJetsJJ << std::endl;
    fMB->Close();
    fJJ->Close();
    delete fMB;
    delete fJJ;
    return false;
  }
  
  TObject *objTrueJJ = dirJJ->Get(kHistTrue);
  TObject *objRecoJJ = dirJJ->Get(kHistReco);
  TObject *objRespJJ = dirJJ->Get(kHistResp2D);
  
  TH1D *hTrueInJJ = dynamic_cast<TH1D *>(objTrueJJ);
  TH1D *hRecoInJJ = dynamic_cast<TH1D *>(objRecoJJ);
  TH2D *hRespInJJ = dynamic_cast<TH2D *>(objRespJJ);
  
  // Convert TH1F/TH2F to TH1D/TH2D if needed (similar to MB MC)
  if (!hTrueInJJ) {
    TH1F *hTrueFJJ = dynamic_cast<TH1F *>(objTrueJJ);
    if (hTrueFJJ) {
      if (hTrueFJJ->GetXaxis()->GetXbins() && hTrueFJJ->GetXaxis()->GetXbins()->GetSize() > 0) {
        hTrueInJJ = new TH1D("hTrueJJ_converted", hTrueFJJ->GetTitle(),
                            hTrueFJJ->GetNbinsX(), hTrueFJJ->GetXaxis()->GetXbins()->GetArray());
      } else {
        hTrueInJJ = new TH1D("hTrueJJ_converted", hTrueFJJ->GetTitle(),
                            hTrueFJJ->GetNbinsX(), hTrueFJJ->GetXaxis()->GetXmin(), 
                            hTrueFJJ->GetXaxis()->GetXmax());
      }
      for (int i = 1; i <= hTrueFJJ->GetNbinsX(); ++i) {
        hTrueInJJ->SetBinContent(i, hTrueFJJ->GetBinContent(i));
        hTrueInJJ->SetBinError(i, hTrueFJJ->GetBinError(i));
      }
    }
  }
  if (!hRecoInJJ) {
    TH1F *hRecoFJJ = dynamic_cast<TH1F *>(objRecoJJ);
    if (hRecoFJJ) {
      if (hRecoFJJ->GetXaxis()->GetXbins() && hRecoFJJ->GetXaxis()->GetXbins()->GetSize() > 0) {
        hRecoInJJ = new TH1D("hRecoJJ_converted", hRecoFJJ->GetTitle(),
                             hRecoFJJ->GetNbinsX(), hRecoFJJ->GetXaxis()->GetXbins()->GetArray());
      } else {
        hRecoInJJ = new TH1D("hRecoJJ_converted", hRecoFJJ->GetTitle(),
                             hRecoFJJ->GetNbinsX(), hRecoFJJ->GetXaxis()->GetXmin(), 
                             hRecoFJJ->GetXaxis()->GetXmax());
      }
      for (int i = 1; i <= hRecoFJJ->GetNbinsX(); ++i) {
        hRecoInJJ->SetBinContent(i, hRecoFJJ->GetBinContent(i));
        hRecoInJJ->SetBinError(i, hRecoFJJ->GetBinError(i));
      }
    }
  }
  if (!hRespInJJ) {
    TH2F *hRespFJJ = dynamic_cast<TH2F *>(objRespJJ);
    if (hRespFJJ) {
      const Double_t *xbinsJJ = hRespFJJ->GetXaxis()->GetXbins() && hRespFJJ->GetXaxis()->GetXbins()->GetSize() > 0 
                                ? hRespFJJ->GetXaxis()->GetXbins()->GetArray() : nullptr;
      const Double_t *ybinsJJ = hRespFJJ->GetYaxis()->GetXbins() && hRespFJJ->GetYaxis()->GetXbins()->GetSize() > 0 
                                ? hRespFJJ->GetYaxis()->GetXbins()->GetArray() : nullptr;
      if (xbinsJJ && ybinsJJ) {
        hRespInJJ = new TH2D("hRespJJ_converted", hRespFJJ->GetTitle(),
                            hRespFJJ->GetNbinsX(), xbinsJJ, hRespFJJ->GetNbinsY(), ybinsJJ);
      } else if (xbinsJJ) {
        hRespInJJ = new TH2D("hRespJJ_converted", hRespFJJ->GetTitle(),
                            hRespFJJ->GetNbinsX(), xbinsJJ, hRespFJJ->GetNbinsY(), 
                            hRespFJJ->GetYaxis()->GetXmin(), hRespFJJ->GetYaxis()->GetXmax());
      } else if (ybinsJJ) {
        hRespInJJ = new TH2D("hRespJJ_converted", hRespFJJ->GetTitle(),
                            hRespFJJ->GetNbinsX(), hRespFJJ->GetXaxis()->GetXmin(), 
                            hRespFJJ->GetXaxis()->GetXmax(), hRespFJJ->GetNbinsY(), ybinsJJ);
      } else {
        hRespInJJ = new TH2D("hRespJJ_converted", hRespFJJ->GetTitle(),
                            hRespFJJ->GetNbinsX(), hRespFJJ->GetXaxis()->GetXmin(), 
                            hRespFJJ->GetXaxis()->GetXmax(),
                            hRespFJJ->GetNbinsY(), hRespFJJ->GetYaxis()->GetXmin(), 
                            hRespFJJ->GetYaxis()->GetXmax());
      }
      for (int i = 1; i <= hRespFJJ->GetNbinsX(); ++i) {
        for (int j = 1; j <= hRespFJJ->GetNbinsY(); ++j) {
          hRespInJJ->SetBinContent(i, j, hRespFJJ->GetBinContent(i, j));
          hRespInJJ->SetBinError(i, j, hRespFJJ->GetBinError(i, j));
        }
      }
    }
  }
  
  if (!hTrueInJJ || !hRecoInJJ || !hRespInJJ) {
    std::cerr << "[Error] Missing histograms in JJ MC response file" << std::endl;
    fMB->Close();
    fJJ->Close();
    delete fMB;
    delete fJJ;
    return false;
  }
  
  nEventsJJ = GetEventCountForTruth(fJJ, kDirJetsJJ, "h_mcColl_counts_weight");
  if (nEventsJJ <= 0.) {
    nEventsJJ = GetEventCount(fJJ, kDirJetsJJ);
    if (nEventsJJ <= 0.) {
      nEventsJJ = hTrueInJJ->Integral();
    }
  }
  
  // Rebin to analysis binning
  TH1D *hTrueMB = dynamic_cast<TH1D *>(hTrueInMB->Rebin(GetNptBinsGen(), "hTrueMB_gen", GetPtbinGen()));
  TH1D *hRecoMB = dynamic_cast<TH1D *>(hRecoInMB->Rebin(GetNptBins(), "hRecoMB_det", GetPtbin()));
  TH1D *hTrueJJ = dynamic_cast<TH1D *>(hTrueInJJ->Rebin(GetNptBinsGen(), "hTrueJJ_gen", GetPtbinGen()));
  TH1D *hRecoJJ = dynamic_cast<TH1D *>(hRecoInJJ->Rebin(GetNptBins(), "hRecoJJ_det", GetPtbin()));
  
  
  // CRITICAL: First rebin 2D response matrices to analysis binning
  // This must be done BEFORE normalization and merging
  const Double_t *ptbin = GetPtbin();
  const Double_t *ptbinGen = GetPtbinGen();
  Int_t nptBins = GetNptBins();
  Int_t nptBinsGen = GetNptBinsGen();
  
  TH2D *hResp2D_MB = new TH2D("hResp2D_MB_rebinned", "MB Response",
                              nptBins, ptbin, nptBinsGen, ptbinGen);
  TH2D *hResp2D_JJ = new TH2D("hResp2D_JJ_rebinned", "JJ Response",
                              nptBins, ptbin, nptBinsGen, ptbinGen);
  
  // Rebin MB MC response matrix
  for (int ix = 1; ix <= hRespInMB->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hRespInMB->GetNbinsY(); ++iy) {
      double content = hRespInMB->GetBinContent(ix, iy);
      double error = hRespInMB->GetBinError(ix, iy);
      if (content == 0.) continue;
      
      double xCenter = hRespInMB->GetXaxis()->GetBinCenter(ix);
      double yCenter = hRespInMB->GetYaxis()->GetBinCenter(iy);
      int xbin = hResp2D_MB->GetXaxis()->FindBin(xCenter);
      int ybin = hResp2D_MB->GetYaxis()->FindBin(yCenter);
      
      double currentContent = hResp2D_MB->GetBinContent(xbin, ybin);
      double currentError = hResp2D_MB->GetBinError(xbin, ybin);
      hResp2D_MB->SetBinContent(xbin, ybin, currentContent + content);
      hResp2D_MB->SetBinError(xbin, ybin, std::sqrt(currentError * currentError + error * error));
    }
  }
  
  // Rebin JJ MC response matrix
  for (int ix = 1; ix <= hRespInJJ->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hRespInJJ->GetNbinsY(); ++iy) {
      double content = hRespInJJ->GetBinContent(ix, iy);
      double error = hRespInJJ->GetBinError(ix, iy);
      if (content == 0.) continue;
      
      double xCenter = hRespInJJ->GetXaxis()->GetBinCenter(ix);
      double yCenter = hRespInJJ->GetYaxis()->GetBinCenter(iy);
      int xbin = hResp2D_JJ->GetXaxis()->FindBin(xCenter);
      int ybin = hResp2D_JJ->GetYaxis()->FindBin(yCenter);
      
      double currentContent = hResp2D_JJ->GetBinContent(xbin, ybin);
      double currentError = hResp2D_JJ->GetBinError(xbin, ybin);
      hResp2D_JJ->SetBinContent(xbin, ybin, currentContent + content);
      hResp2D_JJ->SetBinError(xbin, ybin, std::sqrt(currentError * currentError + error * error));
    }
  }
  
  
  // CRITICAL: Do NOT normalize before merging!
  // Merge first, then normalize with weighted average based on event counts
  // This ensures continuity at the boundary
  std::cerr << "[Info] MB MC event count: " << nEventsMB << std::endl;
  std::cerr << "[Info] JJ MC event count: " << nEventsJJ << std::endl;
  
  // Create merged histograms (unnormalized)
  hTrue = new TH1D("hTrue_merged", "Merged Truth", GetNptBinsGen(), GetPtbinGen());
  hReco = new TH1D("hReco_merged", "Merged Reco", GetNptBins(), GetPtbin());
  hResp2D = new TH2D("hResp2D_merged", "Merged Response;#it{p}_{T}^{reco};#it{p}_{T}^{true}",
                     GetNptBins(), GetPtbin(), GetNptBinsGen(), GetPtbinGen());
  
  hTrue->SetDirectory(nullptr);
  hReco->SetDirectory(nullptr);
  hResp2D->SetDirectory(nullptr);
  
  // CRITICAL FIX: Merge hResp2D FIRST, then derive hTrue and hReco from projections
  // This ensures consistency: hTrue and hReco must match hResp2D projections
  // Merge: Reco pT < 20 AND Truth pT < 20 → MB, otherwise → JJ
  // This approach matches DrawMcClosureTest.C and gives better closure test results
  // - Low pT reco measurements use MB MC response (optimized for low pT)
  // - High pT reco measurements use JJ MC response (optimized for high pT)
  int nMBBins = 0;
  int nJJBins = 0;
  double totalMBContent = 0.0;
  double totalJJContent = 0.0;
  
  // Debug: Check MB and JJ response matrix integrals before merging
  double hResp2D_MB_integral = hResp2D_MB->Integral();
  double hResp2D_JJ_integral = hResp2D_JJ->Integral();
  std::cerr << "[Debug] Before merging: hResp2D_MB integral=" << hResp2D_MB_integral 
            << ", hResp2D_JJ integral=" << hResp2D_JJ_integral << std::endl;
  
  // Count bins in each category for debugging
  int nMB_bins_total = 0, nJJ_bins_total = 0;
  double MB_content_all = 0.0, JJ_content_all = 0.0;
  
  // Also check where JJ MC content actually is
  double JJ_content_in_MB_region = 0.0;
  double JJ_content_in_JJ_region = 0.0;
  int nJJ_bins_in_MB_region = 0;
  int nJJ_bins_in_JJ_region = 0;
  
  for (int ix = 1; ix <= hResp2D->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hResp2D->GetNbinsY(); ++iy) {
      double recoHigh = hResp2D->GetXaxis()->GetBinUpEdge(ix);
      double truthHigh = hResp2D->GetYaxis()->GetBinUpEdge(iy);
      bool useMB = (recoHigh <= kMergeBoundary) && (truthHigh <= kMergeBoundary);
      
      double contentMB = hResp2D_MB->GetBinContent(ix, iy);
      double contentJJ = hResp2D_JJ->GetBinContent(ix, iy);
      
      if (useMB) {
        if (contentMB > 0) {
          nMB_bins_total++;
          MB_content_all += contentMB;
        }
        if (contentJJ > 0) {
          nJJ_bins_in_MB_region++;
          JJ_content_in_MB_region += contentJJ;
        }
      } else {
        if (contentJJ > 0) {
          nJJ_bins_total++;
          JJ_content_all += contentJJ;
          nJJ_bins_in_JJ_region++;
          JJ_content_in_JJ_region += contentJJ;
        }
      }
    }
  }
  std::cerr << "[Debug] Bin counting: nMB_bins=" << nMB_bins_total << " (content=" << MB_content_all 
            << "), nJJ_bins=" << nJJ_bins_total << " (content=" << JJ_content_all << ")" << std::endl;
  std::cerr << "[Debug] JJ MC distribution: In MB region (Reco<20 AND Truth<20): " 
            << nJJ_bins_in_MB_region << " bins, content=" << JJ_content_in_MB_region << std::endl;
  std::cerr << "[Debug] JJ MC distribution: In JJ region (otherwise): " 
            << nJJ_bins_in_JJ_region << " bins, content=" << JJ_content_in_JJ_region << std::endl;
  std::cerr << "[Debug] JJ MC total check: " << JJ_content_in_MB_region << " + " 
            << JJ_content_in_JJ_region << " = " << (JJ_content_in_MB_region + JJ_content_in_JJ_region) 
            << " (should equal " << hResp2D_JJ_integral << ")" << std::endl;
  
  for (int ix = 1; ix <= hResp2D->GetNbinsX(); ++ix) {
    double recoLow = hResp2D->GetXaxis()->GetBinLowEdge(ix);
    double recoHigh = hResp2D->GetXaxis()->GetBinUpEdge(ix);
    
    for (int iy = 1; iy <= hResp2D->GetNbinsY(); ++iy) {
      double truthLow = hResp2D->GetYaxis()->GetBinLowEdge(iy);
      double truthHigh = hResp2D->GetYaxis()->GetBinUpEdge(iy);
      
      // Check boundary: Reco pT < 20 AND Truth pT < 20 → MB MC
      bool useMB = (recoHigh <= kMergeBoundary) && (truthHigh <= kMergeBoundary);
      
      double contentMB = hResp2D_MB->GetBinContent(ix, iy);
      double errorMB = hResp2D_MB->GetBinError(ix, iy);
      double contentJJ = hResp2D_JJ->GetBinContent(ix, iy);
      double errorJJ = hResp2D_JJ->GetBinError(ix, iy);
      
      if (useMB) {
        // Use MB MC
        hResp2D->SetBinContent(ix, iy, contentMB);
        hResp2D->SetBinError(ix, iy, errorMB);
        if (contentMB > 0) {
          nMBBins++;
          totalMBContent += contentMB;
        }
      } else {
        // Use JJ MC
        hResp2D->SetBinContent(ix, iy, contentJJ);
        hResp2D->SetBinError(ix, iy, errorJJ);
        if (contentJJ > 0) {
          nJJBins++;
          totalJJContent += contentJJ;
        }
      }
    }
  }
  
  // Debug: Check merged response matrix integral
  double hResp2D_merged_integral = hResp2D->Integral();
  std::cerr << "[Debug] After merging: hResp2D merged integral=" << hResp2D_merged_integral 
            << ", MB content=" << totalMBContent << ", JJ content=" << totalJJContent << std::endl;
  std::cerr << "[Debug] Expected: MB content should be ~" << hResp2D_MB_integral 
            << " (only Reco<20 AND Truth<20 bins), JJ content should be ~" << hResp2D_JJ_integral 
            << " (all other bins)" << std::endl;
  std::cerr << "[Debug] Comparison: MB_content_all=" << MB_content_all 
            << " vs totalMBContent=" << totalMBContent << std::endl;
  std::cerr << "[Debug] Comparison: JJ_content_all=" << JJ_content_all 
            << " vs totalJJContent=" << totalJJContent << std::endl;
  
  // CRITICAL FIX: Derive hTrue and hReco from hResp2D projections to ensure consistency
  // This fixes the issue where hTrue and hResp2D True projection don't match
  // First, get projections from merged hResp2D
  std::cerr << "[Debug] Step 1: Getting projections from merged hResp2D..." << std::endl;
  TH1D *hTrueProj = dynamic_cast<TH1D *>(hResp2D->ProjectionY("hTrueProj_merged", 1, hResp2D->GetNbinsX(), "e"));
  TH1D *hRecoProj = dynamic_cast<TH1D *>(hResp2D->ProjectionX("hRecoProj_merged", 1, hResp2D->GetNbinsY(), "e"));
  hTrueProj->SetDirectory(nullptr);
  hRecoProj->SetDirectory(nullptr);
  
  double trueProjIntegral_before = hTrueProj->Integral();
  double recoProjIntegral_before = hRecoProj->Integral();
  std::cerr << "[Debug] Step 1 result: True projection integral=" << trueProjIntegral_before 
            << ", Reco projection integral=" << recoProjIntegral_before << std::endl;
  
  // Set hTrue and hReco to match projections (matched entries)
  std::cerr << "[Debug] Step 2: Setting hTrue and hReco to match projections..." << std::endl;
  for (int i = 1; i <= hTrue->GetNbinsX(); ++i) {
    int projBin = hTrueProj->GetXaxis()->FindBin(hTrue->GetXaxis()->GetBinCenter(i));
    if (projBin >= 1 && projBin <= hTrueProj->GetNbinsX()) {
      hTrue->SetBinContent(i, hTrueProj->GetBinContent(projBin));
      hTrue->SetBinError(i, hTrueProj->GetBinError(projBin));
    } else {
      hTrue->SetBinContent(i, 0.);
      hTrue->SetBinError(i, 0.);
    }
  }
  
  for (int i = 1; i <= hReco->GetNbinsX(); ++i) {
    int projBin = hRecoProj->GetXaxis()->FindBin(hReco->GetXaxis()->GetBinCenter(i));
    if (projBin >= 1 && projBin <= hRecoProj->GetNbinsX()) {
      hReco->SetBinContent(i, hRecoProj->GetBinContent(projBin));
      hReco->SetBinError(i, hRecoProj->GetBinError(projBin));
    } else {
      hReco->SetBinContent(i, 0.);
      hReco->SetBinError(i, 0.);
    }
  }
  
  double hTrueIntegral_after_proj = hTrue->Integral();
  double hRecoIntegral_after_proj = hReco->Integral();
  std::cerr << "[Debug] Step 2 result: hTrue integral=" << hTrueIntegral_after_proj 
            << ", hReco integral=" << hRecoIntegral_after_proj << std::endl;
  
  // Now add fake and miss entries from original 1D histograms
  // Merge reco histogram: Reco pT < 20 → MB, >= 20 → JJ (for fake calculation)
  TH1D *hRecoTotal = new TH1D("hRecoTotal", "Total Reco", GetNptBins(), GetPtbin());
  hRecoTotal->SetDirectory(nullptr);
  
  for (int i = 1; i <= hRecoTotal->GetNbinsX(); ++i) {
    double recoHigh = hRecoTotal->GetXaxis()->GetBinUpEdge(i);
    double recoCenter = hRecoTotal->GetXaxis()->GetBinCenter(i);
    
    if (recoHigh <= kMergeBoundary) {
      // Entire bin is below boundary → use MB MC
      int binMB = hRecoMB->GetXaxis()->FindBin(recoCenter);
      if (binMB >= 1 && binMB <= hRecoMB->GetNbinsX()) {
        hRecoTotal->SetBinContent(i, hRecoMB->GetBinContent(binMB));
        hRecoTotal->SetBinError(i, hRecoMB->GetBinError(binMB));
      }
    } else {
      // Bin is at or above boundary → use JJ MC
      int binJJ = hRecoJJ->GetXaxis()->FindBin(recoCenter);
      if (binJJ >= 1 && binJJ <= hRecoJJ->GetNbinsX()) {
        hRecoTotal->SetBinContent(i, hRecoJJ->GetBinContent(binJJ));
        hRecoTotal->SetBinError(i, hRecoJJ->GetBinError(binJJ));
      }
    }
  }
  
  // Merge truth histogram: Truth pT < 20 → MB, >= 20 → JJ (for miss calculation)
  // CRITICAL: hTrueTotal must include ALL truth entries, regardless of reco pT
  // This is different from hResp2D which only includes matched entries
  TH1D *hTrueTotal = new TH1D("hTrueTotal", "Total True", GetNptBinsGen(), GetPtbinGen());
  hTrueTotal->SetDirectory(nullptr);
  
  // Simply add MB and JJ truth histograms - they already contain all truth entries
  // This ensures we get ALL truth entries, not just matched ones
  for (int i = 1; i <= hTrueTotal->GetNbinsX(); ++i) {
    double truthHigh = hTrueTotal->GetXaxis()->GetBinUpEdge(i);
    double truthCenter = hTrueTotal->GetXaxis()->GetBinCenter(i);
    
    double content = 0.0;
    double error = 0.0;
    
    if (truthHigh <= kMergeBoundary) {
      // Entire bin is below boundary → use MB MC
      int binMB = hTrueMB->GetXaxis()->FindBin(truthCenter);
      if (binMB >= 1 && binMB <= hTrueMB->GetNbinsX()) {
        content = hTrueMB->GetBinContent(binMB);
        error = hTrueMB->GetBinError(binMB);
      }
    } else {
      // Bin is at or above boundary → use JJ MC
      int binJJ = hTrueJJ->GetXaxis()->FindBin(truthCenter);
      if (binJJ >= 1 && binJJ <= hTrueJJ->GetNbinsX()) {
        content = hTrueJJ->GetBinContent(binJJ);
        error = hTrueJJ->GetBinError(binJJ);
      }
    }
    
    hTrueTotal->SetBinContent(i, content);
    hTrueTotal->SetBinError(i, error);
  }
  
  std::cerr << "[Debug] Step 3a_detail: hTrueMB integral=" << hTrueMB->Integral() 
            << ", hTrueJJ integral=" << hTrueJJ->Integral() 
            << ", hTrueTotal integral=" << hTrueTotal->Integral() << std::endl;
  
  // Add fake and miss to hReco and hTrue
  // Fake = Reco total - Reco matched (from projection)
  std::cerr << "[Debug] Step 3: Calculating fake and miss..." << std::endl;
  double hRecoTotalIntegral = hRecoTotal->Integral();
  double hTrueTotalIntegral = hTrueTotal->Integral();
  std::cerr << "[Debug] Step 3a: hRecoTotal integral=" << hRecoTotalIntegral 
            << ", hTrueTotal integral=" << hTrueTotalIntegral << std::endl;
  
  TH1D *hFake = dynamic_cast<TH1D *>(hRecoTotal->Clone("hFake"));
  hFake->Add(hRecoProj, -1);
  double hFakeIntegral = hFake->Integral();
  std::cerr << "[Debug] Step 3b: hFake integral=" << hFakeIntegral << std::endl;
  
  TH1D *hMiss = dynamic_cast<TH1D *>(hTrueTotal->Clone("hMiss"));
  hMiss->Add(hTrueProj, -1);
  double hMissIntegral = hMiss->Integral();
  std::cerr << "[Debug] Step 3c: hMiss integral=" << hMissIntegral << std::endl;
  
  std::cerr << "[Debug] Step 4: Adding fake and miss to hReco and hTrue..." << std::endl;
  hReco->Add(hFake);
  hTrue->Add(hMiss);
  
  double hTrueIntegral_final = hTrue->Integral();
  double hRecoIntegral_final = hReco->Integral();
  std::cerr << "[Debug] Step 4 result: hTrue final integral=" << hTrueIntegral_final 
            << ", hReco final integral=" << hRecoIntegral_final << std::endl;
  std::cerr << "[Debug] Step 4 check: hTrue should equal hTrueTotal=" << hTrueTotalIntegral 
            << ", difference=" << (hTrueIntegral_final - hTrueTotalIntegral) << std::endl;
  std::cerr << "[Debug] Step 4 check: hReco should equal hRecoTotal=" << hRecoTotalIntegral 
            << ", difference=" << (hRecoIntegral_final - hRecoTotalIntegral) << std::endl;
  
  // Cleanup temporary histograms
  delete hRecoTotal;
  delete hTrueTotal;
  delete hFake;
  delete hMiss;
  
  // Debug: Check 1D merge consistency
  int nTrueMB = 0, nTrueJJ = 0;
  double trueMBContent = 0.0, trueJJContent = 0.0;
  int nRecoMB = 0, nRecoJJ = 0;
  double recoMBContent = 0.0, recoJJContent = 0.0;
  
  for (int i = 1; i <= hTrue->GetNbinsX(); ++i) {
    double truthHigh = hTrue->GetXaxis()->GetBinUpEdge(i);
    if (truthHigh <= kMergeBoundary) {
      double content = hTrue->GetBinContent(i);
      if (content > 0) {
        nTrueMB++;
        trueMBContent += content;
      }
    } else {
      double content = hTrue->GetBinContent(i);
      if (content > 0) {
        nTrueJJ++;
        trueJJContent += content;
      }
    }
  }
  
  for (int i = 1; i <= hReco->GetNbinsX(); ++i) {
    double recoHigh = hReco->GetXaxis()->GetBinUpEdge(i);
    if (recoHigh <= kMergeBoundary) {
      double content = hReco->GetBinContent(i);
      if (content > 0) {
        nRecoMB++;
        recoMBContent += content;
      }
    } else {
      double content = hReco->GetBinContent(i);
      if (content > 0) {
        nRecoJJ++;
        recoJJContent += content;
      }
    }
  }
  
  std::cerr << "[Debug] 1D Histogram Merge Check:" << std::endl;
  std::cerr << "  hTrue: MB bins=" << nTrueMB << " (content=" << trueMBContent << "), JJ bins=" << nTrueJJ << " (content=" << trueJJContent << ")" << std::endl;
  std::cerr << "  hReco: MB bins=" << nRecoMB << " (content=" << recoMBContent << "), JJ bins=" << nRecoJJ << " (content=" << recoJJContent << ")" << std::endl;
  
  // Cleanup temporary histograms
  delete hTrueProj;
  delete hRecoProj;
  
  
  // CRITICAL: Do NOT normalize by event count before merging!
  // Response matrix should be y-axis normalized (migration probability), not event-count normalized
  // Event count normalization will be done AFTER unfolding for cross section calculation
  // This matches DrawMcClosureTest.C behavior (no normalization in LoadMergedResponseSample)
  
  // Debug: Check 2D merge consistency with 1D
  // Verify that response matrix projection matches 1D histograms
  // Note: hTrueProj and hRecoProj were already created above and used, then deleted
  // Recreate them here for consistency check
  std::cerr << "[Debug] Step 5: Final consistency check..." << std::endl;
  TH1D *hRecoProjCheck = dynamic_cast<TH1D *>(hResp2D->ProjectionX("hRecoProj_check", 1, hResp2D->GetNbinsY(), "e"));
  TH1D *hTrueProjCheck = dynamic_cast<TH1D *>(hResp2D->ProjectionY("hTrueProj_check", 1, hResp2D->GetNbinsX(), "e"));
  
  double recoProjIntegral = hRecoProjCheck->Integral();
  double recoIntegral = hReco->Integral();
  double trueProjIntegral = hTrueProjCheck->Integral();
  double trueIntegral = hTrue->Integral();
  
  double recoRatio = (recoIntegral > 0) ? recoProjIntegral / recoIntegral : 0.0;
  double trueRatio = (trueIntegral > 0) ? trueProjIntegral / trueIntegral : 0.0;
  
  std::cerr << "[Debug] 2D Response Matrix Merge Check:" << std::endl;
  std::cerr << "  MB MC bins (Reco < " << kMergeBoundary << " AND Truth < " << kMergeBoundary << "): " << nMBBins << " (content=" << totalMBContent << ")" << std::endl;
  std::cerr << "  JJ MC bins (otherwise): " << nJJBins << " (content=" << totalJJContent << ")" << std::endl;
  std::cerr << "  Response matrix projection vs 1D histogram:" << std::endl;
  std::cerr << "    Reco: projection=" << recoProjIntegral << ", 1D=" << recoIntegral << ", ratio=" << recoRatio << std::endl;
  std::cerr << "    True: projection=" << trueProjIntegral << ", 1D=" << trueIntegral << ", ratio=" << trueRatio << std::endl;
  std::cerr << "    [NOTE] True ratio should be projection/(projection+miss) = " 
            << trueProjIntegral << "/" << trueIntegral << " = " << trueRatio << std::endl;
  std::cerr << "    [NOTE] Reco ratio should be projection/(projection+fake) = " 
            << recoProjIntegral << "/" << recoIntegral << " = " << recoRatio << std::endl;
  if (TMath::Abs(recoRatio - 1.0) > 0.1 || TMath::Abs(trueRatio - 1.0) > 0.1) {
    std::cerr << "  [WARNING] Large discrepancy between 2D projection and 1D histogram!" << std::endl;
    std::cerr << "  [WARNING] This is expected if there are significant fake/miss entries!" << std::endl;
  }
  
  delete hRecoProjCheck;
  delete hTrueProjCheck;
  
  std::cerr << "[Info] Response Matrix merge summary:" << std::endl;
  std::cerr << "  MB MC bins (Reco < " << kMergeBoundary << " AND Truth < " << kMergeBoundary << "): " << nMBBins << std::endl;
  std::cerr << "  JJ MC bins (otherwise): " << nJJBins << std::endl;
  std::cerr << "========================================" << std::endl;
  
  // Cleanup
  delete hTrueMB;
  delete hRecoMB;
  delete hTrueJJ;
  delete hRecoJJ;
  delete hResp2D_MB;
  delete hResp2D_JJ;
  fMB->Close();
  fJJ->Close();
  delete fMB;
  delete fJJ;
  
  std::cerr << "[Info] Merged response sample loaded successfully!" << std::endl;
  std::cerr << "  MB MC event count: " << nEventsMB << std::endl;
  std::cerr << "  JJ MC event count: " << nEventsJJ << std::endl;
  
  return true;
}

// ============================================================
// Build RooUnfoldResponse (kernel-style, same as DrawMcClosureTest.C)
// ============================================================

// Build RooUnfoldResponse from rebinned truth, reco and 2D matrix
// KERNEL-STYLE: Response is built from matched pairs only (no Fake/Miss filled)
// Purity (fake) and efficiency (miss) are handled externally
// CRITICAL: Normalize response matrix by y-axis (truth pT) to avoid scale issues when merging MB+JJ MC
RooUnfoldResponse *BuildResponse(TH1D *hTrue, TH1D *hReco, TH2D *hResp2D, double &scaleFactor,
                                  TH1D *&hPurity, TH1D *&hEfficiency) {

  // Set scaleFactor to 1.0 (no scaling, same as DrawMcClosureTest.C)
  scaleFactor = 1.0;
  hPurity = nullptr;
  hEfficiency = nullptr;
  
  // CRITICAL: Calculate matched projections from ORIGINAL (non-normalized) RM
  // These are needed to calculate fake and miss correctly
  TH1D *hRecoMatched =
      dynamic_cast<TH1D *>(hResp2D->ProjectionX("hRecoMatched", 1,
                                                hResp2D->GetNbinsY(), "e"));
  TH1D *hTrueMatched =
      dynamic_cast<TH1D *>(hResp2D->ProjectionY("hTrueMatched", 1,
                                                hResp2D->GetNbinsX(), "e"));

  // CRITICAL: Calculate fake and miss from ORIGINAL counts (before normalization)
  // Fake and miss must use original counts (not normalized probabilities)
  // because they represent absolute numbers of unmatched particles
  // Fake: reco particles without truth match (reco axis)
  // Miss: truth particles without reco match (truth axis)
  TH1D *fake = dynamic_cast<TH1D *>(hReco->Clone("hFake"));
  fake->Add(hRecoMatched, -1);
  TH1D *miss = dynamic_cast<TH1D *>(hTrue->Clone("hMiss"));
  miss->Add(hTrueMatched, -1);

  std::cerr << "[Info] Fake and miss calculated from original counts (before normalization)" << std::endl;

  // KERNEL-STYLE: Compute purity and efficiency for external application
  // Purity(pT_reco) = matchedReco / reco
  // Efficiency(pT_true) = matchedTruth / truth
  hPurity = dynamic_cast<TH1D *>(hReco->Clone("hPurity_response"));
  if (hPurity) {
    hPurity->SetDirectory(nullptr);
    hPurity->Reset();
    for (int i = 1; i <= hPurity->GetNbinsX(); ++i) {
      const double denom = hReco->GetBinContent(i);
      const double numer = hRecoMatched->GetBinContent(i);
      double p = (denom > 0) ? (numer / denom) : 0.0;
      if (p < 0) p = 0.0;
      if (p > 1) p = 1.0;
      hPurity->SetBinContent(i, p);
    }
  }
  hEfficiency = dynamic_cast<TH1D *>(hTrue->Clone("hEfficiency_response"));
  if (hEfficiency) {
    hEfficiency->SetDirectory(nullptr);
    hEfficiency->Reset();
    for (int i = 1; i <= hEfficiency->GetNbinsX(); ++i) {
      const double denom = hTrue->GetBinContent(i);
      const double numer = hTrueMatched->GetBinContent(i);
      double eff = (denom > 0) ? (numer / denom) : 0.0;
      if (eff < 0) eff = 0.0;
      if (eff > 1) eff = 1.0;
      hEfficiency->SetBinContent(i, eff);
    }
  }
  std::cerr << "[Info] Purity and efficiency computed for kernel-style corrections" << std::endl;
  
  // Clone the rebinned 2D matrix as Respt for normalization
  TH2D *Respt = dynamic_cast<TH2D *>(hResp2D->Clone("Respt"));
  
  // CRITICAL: Normalize response matrix by y-axis (truth pT) - each row sums to 1.0
  // This ensures MB MC and JJ MC can be merged without scale issues
  // Each truth bin (y-axis) row is normalized independently
  std::cerr << "[Info] Normalizing response matrix by y-axis (truth pT)..." << std::endl;
  for (int j = 1; j <= Respt->GetNbinsY(); ++j) {
    // Calculate row sum for this truth bin
    double rowSum = 0.0;
    double rowSumErrorSq = 0.0;
    for (int i = 1; i <= Respt->GetNbinsX(); ++i) {
      double content = Respt->GetBinContent(i, j);
      double error = Respt->GetBinError(i, j);
      rowSum += content;
      rowSumErrorSq += error * error;
    }
    
    // Normalize each bin in this row by row sum
    if (rowSum > 0.0) {
      for (int i = 1; i <= Respt->GetNbinsX(); ++i) {
        double content = Respt->GetBinContent(i, j);
        double error = Respt->GetBinError(i, j);
        // Normalize content
        Respt->SetBinContent(i, j, content / rowSum);
        // Normalize error (propagate uncertainty from normalization)
        double normalizedError = std::sqrt((error * error) / (rowSum * rowSum) + 
                                           (content * content * rowSumErrorSq) / (rowSum * rowSum * rowSum * rowSum));
        Respt->SetBinError(i, j, normalizedError);
      }
    }
  }
  std::cerr << "[Info] Response matrix normalized: each truth bin column sums to 1.0" << std::endl;

  // Create RooUnfoldResponse with MATCHED histograms only (no Fake/Miss)
  // 3-arg constructor: directly sets _res=Respt, _mes=hRecoMatched, _tru=hTrueMatched
  // No Fill() loop — avoids doubling _mes/_tru that occurs with 2-arg + Fill pattern
  // Preserves original bin errors (Sumw2) from the 2D matrix
  RooUnfoldResponse *resp = new RooUnfoldResponse(hRecoMatched, hTrueMatched, Respt);

  // KERNEL-STYLE: Do NOT add Miss/Fake to RooUnfoldResponse
  // Purity (fake) and efficiency (miss) are handled externally via hPurity and hEfficiency

  std::cerr << "[Info] RooUnfoldResponse built with 3-arg constructor (matched pairs, raw counts)" << std::endl;
  
  // Debug: Check fake/miss statistics
  double missTotal = 0.0, fakeTotal = 0.0;
  for (int i = 1; i <= miss->GetNbinsX(); ++i) {
    missTotal += miss->GetBinContent(i);
  }
  for (int i = 1; i <= fake->GetNbinsX(); ++i) {
    fakeTotal += fake->GetBinContent(i);
  }
  
  double recoMatchedIntegral = hRecoMatched->Integral();
  double trueMatchedIntegral = hTrueMatched->Integral();
  double recoIntegral = hReco->Integral();
  double trueIntegral = hTrue->Integral();
  
  std::cerr << "[Debug] Response Matrix Statistics:" << std::endl;
  std::cerr << "  Reco total=" << recoIntegral << ", matched=" << recoMatchedIntegral << ", fake=" << fakeTotal << std::endl;
  std::cerr << "  True total=" << trueIntegral << ", matched=" << trueMatchedIntegral << ", miss=" << missTotal << std::endl;
  if (recoIntegral > 0) {
    std::cerr << "  Fake fraction=" << (fakeTotal / recoIntegral) << std::endl;
  }
  if (trueIntegral > 0) {
    std::cerr << "  Miss fraction=" << (missTotal / trueIntegral) << std::endl;
  }
  if (missTotal > 0.1 * trueIntegral || fakeTotal > 0.1 * recoIntegral) {
    std::cerr << "  [WARNING] Large fake/miss fractions detected!" << std::endl;
  }
  
  delete Respt;
  // Note: hRecoMatched and hTrueMatched are kept for potential QA
  // delete hRecoMatched;
  // delete hTrueMatched;
  delete fake;
  delete miss;

  return resp;
}

// ============================================================
// Helper functions for reference data
// ============================================================

// Run 2 data points (hard-coded)
std::vector<TGraphErrors*> Run2Data_manual() {
  const int Npoint = 19;
  
  double pt[Npoint] = {5.5, 6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17, 19, 22.5, 27.5, 35, 45, 55, 65, 77.5, 92.5, 120};
  
  // Original data for non-UESUB case
  double sigma[Npoint] = {1.0972, 0.58418, 0.33224, 0.20001, 0.12614, 0.069825, 0.034225, 0.018498, 0.010764, 0.0066421,
                          0.0033177, 0.0013313, 0.00046794, 0.00013985, 5.2994e-05, 2.3472e-05, 1.0087e-05, 4.1091e-06, 1.2028e-06};
  
  double stat_err[Npoint] = {0.00016396, 0.00010887, 7.5637e-05, 5.3822e-05, 3.7416e-05, 2.2965e-05, 1.3182e-05, 8.6563e-06, 5.8873e-06, 4.1364e-06,
                             2.2742e-06, 1.136e-06, 4.9698e-07, 1.7976e-07, 7.9569e-08, 3.9801e-08, 1.8747e-08, 8.1518e-09, 2.4881e-09};
  
  double syst_err[Npoint] = {0.074335, 0.042393, 0.024799, 0.015407, 0.010138, 0.0057313, 0.0029103, 0.001585, 0.00094719, 0.00061878,
                             0.00031195, 0.00012884, 4.6138e-05, 1.4546e-05, 5.6228e-06, 2.5033e-06, 1.1183e-06, 4.4943e-07, 1.3178e-07};
  
  TGraphErrors *gr_stat_error = new TGraphErrors(Npoint);
  TGraphErrors *gr_syst_error = new TGraphErrors(Npoint);
  
  for (int i = 0; i < Npoint; i++) {
    double ex1 = (i < Npoint - 1) ? (pt[i+1] - pt[i]) / 2.0 : (pt[i] - pt[i-1]) / 2.0;
    
    gr_stat_error->SetPoint(i, pt[i], sigma[i]);
    gr_stat_error->SetPointError(i, ex1, stat_err[i]);
    
    gr_syst_error->SetPoint(i, pt[i], sigma[i]);
    gr_syst_error->SetPointError(i, ex1, syst_err[i]);
  }
  
  std::vector<TGraphErrors*> graphs = {gr_stat_error, gr_syst_error};
  return graphs;
}

// Convert TGraphErrors to TH1D
TH1D* GraphToHistogram(TGraphErrors* graph) {
  if (!graph) return nullptr;
  
  int nPoints = graph->GetN();
  const Double_t* ptbin = GetPtbin();
  int nptBins = GetNptBins();
  
  TH1D* hist = new TH1D(Form("hist_from_graph_%p", graph), "Histogram from Graph", nptBins, ptbin);
  hist->SetDirectory(nullptr);
  
  for (int i = 0; i < nPoints; ++i) {
    double x, y;
    graph->GetPoint(i, x, y);
    double ey = graph->GetErrorY(i);
    
    int bin = hist->FindBin(x);
    if (bin >= 1 && bin <= hist->GetNbinsX()) {
      hist->SetBinContent(bin, y);
      hist->SetBinError(bin, ey);
    }
  }
  
  return hist;
}

// Load PYTHIA 13600 data
TH1* DrawPythia13600() {
  std::cerr << "[Info] Loading PYTHIA 13600 data" << std::endl;
  
  TFile* fPythia13600 = TFile::Open("/Users/js/alice/pythiaGen/postprocess/results/pp_13600GeV_HardQCD_all_on_UE_ISR_FSR_on/PYTHIA_pp_13600_GeV.root", "read");
  
  if (!fPythia13600 || fPythia13600->IsZombie()) {
    std::cerr << "[Error] Cannot open PYTHIA file" << std::endl;
    return nullptr;
  }
  
  TH1* hPYTHIA = (TH1*)fPythia13600->Get("hJetPt");
  if (!hPYTHIA) {
    std::cerr << "[Error] hJetPt histogram not found in PYTHIA file" << std::endl;
    fPythia13600->Close();
    delete fPythia13600;
    return nullptr;
  }
  
  TH1* hNevents = (TH1*)fPythia13600->Get("hnevent");
  double Nevts = hNevents ? hNevents->GetBinContent(1) : 100000000.0;
  
  const Double_t ptbinPYTHIA[26] = {0,  1,  2,  3,  4,  5,   6,   7,  8,
    9,  10, 12, 14, 16, 18,  20,  25, 30,
    40, 50, 60, 70, 85, 100, 140, 200};
  
  TH1* hRebinned = hPYTHIA->Rebin(25, "hPYTHIA13600_rebinned", ptbinPYTHIA);
  hRebinned->SetDirectory(nullptr);
  
  hRebinned->Scale(1.0 / 100000000, "width");
  hRebinned->SetLineColor(kCyan);
  hRebinned->SetLineStyle(1);
  hRebinned->SetLineWidth(2);
  hRebinned->SetMarkerStyle(0);
  
  fPythia13600->Close();
  delete fPythia13600;
  
  return hRebinned;
}

// Load Changwhan 2022 JJ MC data
TH1* DrawChangwhan2022JJMC() {
  std::cerr << "[Info] Loading Changwhan 2022 JJ MC data" << std::endl;
  
  TFile* fChangwhan2022JJMC = TFile::Open("/Users/js/cernbox/workspace/O2Physics/macros/DrawJets/drawJetSpectraCharged/Changhwan/jetptincljet_unfolded_22o_25a2b_lumi.root", "read");
  
  if (!fChangwhan2022JJMC || fChangwhan2022JJMC->IsZombie()) {
    std::cerr << "[Error] Cannot open Changwhan JJ MC file" << std::endl;
    return nullptr;
  }
  
  TH1* hChangwhan2022JJMC = (TH1*)fChangwhan2022JJMC->Get("unfolded");
  if (!hChangwhan2022JJMC) {
    std::cerr << "[Error] unfolded histogram not found in Changwhan JJ MC file" << std::endl;
    fChangwhan2022JJMC->Close();
    delete fChangwhan2022JJMC;
    return nullptr;
  }
  
  // Clone and detach from file
  TH1* hClone = (TH1*)hChangwhan2022JJMC->Clone("hChangwhan2022JJMC_clone");
  hClone->SetDirectory(nullptr);
  
  fChangwhan2022JJMC->Close();
  delete fChangwhan2022JJMC;
  
  return hClone;
}

// ============================================================
// Load data file
// ============================================================

bool LoadData(TH1D *&hRecoData, Double_t &nEventsData, Double_t &Ntvx) {
  std::cerr << "[Info] Loading data file: " << kDataFile << std::endl;
  
  TFile *fData = TFile::Open(kDataFile, "READ");
  if (!fData || fData->IsZombie()) {
    std::cerr << "[Error] Cannot open data file" << std::endl;
    return false;
  }
  
  TDirectory *dirData = dynamic_cast<TDirectory *>(fData->Get(kDirJetsData));
  if (!dirData) {
    std::cerr << "[Error] Directory not found: " << kDirJetsData << std::endl;
    fData->Close();
    delete fData;
    return false;
  }
  
  TObject *objRecoData = dirData->Get(kHistReco);
  TH1D *hRecoInData = dynamic_cast<TH1D *>(objRecoData);
  
  if (!hRecoInData) {
    TH1F *hRecoFData = dynamic_cast<TH1F *>(objRecoData);
    if (hRecoFData) {
      if (hRecoFData->GetXaxis()->GetXbins() && hRecoFData->GetXaxis()->GetXbins()->GetSize() > 0) {
        hRecoInData = new TH1D("hRecoData_converted", hRecoFData->GetTitle(),
                               hRecoFData->GetNbinsX(), hRecoFData->GetXaxis()->GetXbins()->GetArray());
      } else {
        hRecoInData = new TH1D("hRecoData_converted", hRecoFData->GetTitle(),
                               hRecoFData->GetNbinsX(), hRecoFData->GetXaxis()->GetXmin(), 
                               hRecoFData->GetXaxis()->GetXmax());
      }
      for (int i = 1; i <= hRecoFData->GetNbinsX(); ++i) {
        hRecoInData->SetBinContent(i, hRecoFData->GetBinContent(i));
        hRecoInData->SetBinError(i, hRecoFData->GetBinError(i));
      }
    }
  }
  
  if (!hRecoInData) {
    std::cerr << "[Error] Cannot load reco histogram from data file" << std::endl;
    fData->Close();
    delete fData;
    return false;
  }
  
  // Rebin to analysis binning
  hRecoData = dynamic_cast<TH1D *>(hRecoInData->Rebin(GetNptBins(), "hRecoData_rebinned", GetPtbin()));
  hRecoData->SetDirectory(nullptr);
  
  // Get event count
  nEventsData = GetEventCountForUnfolded(fData, kDirJetsData);
  if (nEventsData <= 0.) {
    nEventsData = GetEventCount(fData, kDirJetsData);
    if (nEventsData <= 0.) {
      nEventsData = hRecoData->Integral();
    }
  }
  
  // Get Ntvx for cross section calculation
  // Use hCounterTVX instead of hLumiTVX
  TH1 *hNtvx = dynamic_cast<TH1 *>(fData->Get("eventselection-run3/luminosity/hCounterTVX"));
  if (hNtvx) {
    Ntvx = hNtvx->Integral(1, hNtvx->GetNbinsX());
    std::cerr << "[Info] Ntvx from hCounterTVX: " << Ntvx << std::endl;
  } else {
    std::cerr << "[Warning] Cannot find hCounterTVX histogram, trying fallback..." << std::endl;
    hNtvx = dynamic_cast<TH1 *>(fData->Get(GetLumiTVXObj()));
    if (hNtvx) {
      Ntvx = hNtvx->Integral(1, hNtvx->GetNbinsX());
      std::cerr << "[Warning] Using hLumiTVX as fallback: " << Ntvx << std::endl;
    } else {
      std::cerr << "[Error] Cannot find any NTVX histogram, using nEventsData as fallback" << std::endl;
      Ntvx = nEventsData;
    }
  }
  
  fData->Close();
  delete fData;
  
  std::cerr << "[Info] Data loaded successfully!" << std::endl;
  std::cerr << "  Event count: " << nEventsData << std::endl;
  std::cerr << "  Ntvx: " << Ntvx << std::endl;
  
  return true;
}

// ============================================================
// Helper function: Draw QA plots (purity, efficiency, response matrix)
// ============================================================

void DrawQAPlots(TH1D *hTrue, TH1D *hReco, TH2D *hResp2D, 
                 const char *plotTitle, const char *outputPrefix) {
  std::cerr << "========================================" << std::endl;
  std::cerr << "[Info] Drawing QA plots: " << plotTitle << std::endl;
  std::cerr << "========================================" << std::endl;
  
  // Create output directory
  TString outputDir = "plots/AN_Charged-particle-jet-cross-section-in-pp-collisions-at-13.6-TeV/Figures/QA";
  gSystem->MakeDirectory(outputDir.Data());
  
  // CRITICAL: Load original 2D correlation histogram and rebin it directly
  // This matches DrawJetsMCfJetQA.h approach where HCorrelate2D is loaded and rebinned directly
  // instead of using the already-rebinned hResp2D from LoadResponseSample
  TH2D *h2HCorrelate = nullptr;
  
  const Double_t *ptbin = GetPtbin();
  const Double_t *ptbinGen = GetPtbinGen();
  Int_t nptBins = GetNptBins();
  Int_t nptBinsGen = GetNptBinsGen();
  
  if (TString(outputPrefix) == "MergedMBJJ") {
    // For merged case, we need to merge MB and JJ 2D histograms
    // This is complex, so for now we'll use the passed hResp2D
    // TODO: Implement proper merging if needed
    h2HCorrelate = dynamic_cast<TH2D *>(hResp2D->Clone("h2HCorrelate_merged"));
    h2HCorrelate->SetDirectory(nullptr);
  } else if (TString(outputPrefix) == "MBonly") {
    TFile *f = TFile::Open(kResponseFileMB, "READ");
    if (f && !f->IsZombie()) {
      TDirectory *dir = dynamic_cast<TDirectory *>(f->Get(kDirJets));
      if (dir) {
        TH2 *HCorrelate2D = dynamic_cast<TH2 *>(dir->Get(kHistResp2D));
        if (HCorrelate2D) {
          if (GetREBINON()) {
            // Rebin directly from original (same as DrawJetsMCfJetQA.h)
            h2HCorrelate = new TH2D("h2HCorrelate_MB", "MB Correlation",
                                   nptBins, ptbin, nptBinsGen, ptbinGen);
            for (Int_t ix = 1; ix <= HCorrelate2D->GetNbinsX(); ++ix) {
              for (Int_t iy = 1; iy <= HCorrelate2D->GetNbinsY(); ++iy) {
                Double_t content = HCorrelate2D->GetBinContent(ix, iy);
                Double_t error = HCorrelate2D->GetBinError(ix, iy);
                // NOTE: Do NOT skip content == 0, same as DrawJetsMCfJetQA.h (line 232-242)
                Double_t xCenter = HCorrelate2D->GetXaxis()->GetBinCenter(ix);
                Double_t yCenter = HCorrelate2D->GetYaxis()->GetBinCenter(iy);
                Int_t xbin = h2HCorrelate->GetXaxis()->FindBin(xCenter);
                Int_t ybin = h2HCorrelate->GetYaxis()->FindBin(yCenter);
                Double_t currentContent = h2HCorrelate->GetBinContent(xbin, ybin);
                Double_t currentError = h2HCorrelate->GetBinError(xbin, ybin);
                h2HCorrelate->SetBinContent(xbin, ybin, currentContent + content);
                h2HCorrelate->SetBinError(xbin, ybin, std::sqrt(currentError * currentError + error * error));
              }
            }
          } else {
            h2HCorrelate = dynamic_cast<TH2D *>(HCorrelate2D->Clone("h2HCorrelate_MB"));
          }
          h2HCorrelate->SetDirectory(nullptr);
        }
      }
      f->Close();
      delete f;
    }
  } else if (TString(outputPrefix) == "JJonly") {
    TFile *f = TFile::Open(kResponseFileJJ, "READ");
    if (f && !f->IsZombie()) {
      TDirectory *dir = dynamic_cast<TDirectory *>(f->Get(kDirJetsJJ));
      if (dir) {
        TH2 *HCorrelate2D = dynamic_cast<TH2 *>(dir->Get(kHistResp2D));
        if (HCorrelate2D) {
          if (GetREBINON()) {
            // Rebin directly from original (same as DrawJetsMCfJetQA.h)
            h2HCorrelate = new TH2D("h2HCorrelate_JJ", "JJ Correlation",
                                   nptBins, ptbin, nptBinsGen, ptbinGen);
            for (Int_t ix = 1; ix <= HCorrelate2D->GetNbinsX(); ++ix) {
              for (Int_t iy = 1; iy <= HCorrelate2D->GetNbinsY(); ++iy) {
                Double_t content = HCorrelate2D->GetBinContent(ix, iy);
                Double_t error = HCorrelate2D->GetBinError(ix, iy);
                // NOTE: Do NOT skip content == 0, same as DrawJetsMCfJetQA.h (line 232-242)
                Double_t xCenter = HCorrelate2D->GetXaxis()->GetBinCenter(ix);
                Double_t yCenter = HCorrelate2D->GetYaxis()->GetBinCenter(iy);
                Int_t xbin = h2HCorrelate->GetXaxis()->FindBin(xCenter);
                Int_t ybin = h2HCorrelate->GetYaxis()->FindBin(yCenter);
                Double_t currentContent = h2HCorrelate->GetBinContent(xbin, ybin);
                Double_t currentError = h2HCorrelate->GetBinError(xbin, ybin);
                h2HCorrelate->SetBinContent(xbin, ybin, currentContent + content);
                h2HCorrelate->SetBinError(xbin, ybin, std::sqrt(currentError * currentError + error * error));
              }
            }
          } else {
            h2HCorrelate = dynamic_cast<TH2D *>(HCorrelate2D->Clone("h2HCorrelate_JJ"));
          }
          h2HCorrelate->SetDirectory(nullptr);
        }
      }
      f->Close();
      delete f;
    }
  }
  
  // Fallback to passed hResp2D if loading failed
  if (!h2HCorrelate) {
    std::cerr << "[Warning] Failed to load original 2D histogram, using passed hResp2D" << std::endl;
    h2HCorrelate = dynamic_cast<TH2D *>(hResp2D->Clone("h2HCorrelate_fallback"));
    h2HCorrelate->SetDirectory(nullptr);
  }
  
  // Get projections from rebinned correlation matrix (same as DrawJetsMCfJetQA.h)
  TH1D *hRecoMatched = dynamic_cast<TH1D *>(h2HCorrelate->ProjectionX("hRecoMatched", 1, h2HCorrelate->GetNbinsY(), "e"));
  TH1D *hTrueMatched = dynamic_cast<TH1D *>(h2HCorrelate->ProjectionY("hTrueMatched", 1, h2HCorrelate->GetNbinsX(), "e"));
  hRecoMatched->SetDirectory(nullptr);
  hTrueMatched->SetDirectory(nullptr);
  
  // CRITICAL: Load original h_jet_pt and h_jet_pt_part from MC files
  // This matches DrawJetsMCfJetQA.h approach where JetMCDPt and JetMCPPt are loaded directly
  TH1D *hRecoTotal = nullptr;
  TH1D *hTrueTotal = nullptr;
  
  if (TString(outputPrefix) == "MergedMBJJ") {
    // For merged case, need to load both MB and JJ MC files
    TFile *fMB = TFile::Open(kResponseFileMB, "READ");
    TFile *fJJ = TFile::Open(kResponseFileJJ, "READ");
    if (fMB && !fMB->IsZombie() && fJJ && !fJJ->IsZombie()) {
      TDirectory *dirMB = dynamic_cast<TDirectory *>(fMB->Get(kDirJets));
      TDirectory *dirJJ = dynamic_cast<TDirectory *>(fJJ->Get(kDirJetsJJ));
      if (dirMB && dirJJ) {
        TObject *objRecoMB = dirMB->Get(kHistReco);
        TObject *objRecoJJ = dirJJ->Get(kHistReco);
        TObject *objTrueMB = dirMB->Get(kHistTrue);
        TObject *objTrueJJ = dirJJ->Get(kHistTrue);
        
        TH1D *hRecoMB = dynamic_cast<TH1D *>(objRecoMB);
        TH1D *hRecoJJ = dynamic_cast<TH1D *>(objRecoJJ);
        TH1D *hTrueMB = dynamic_cast<TH1D *>(objTrueMB);
        TH1D *hTrueJJ = dynamic_cast<TH1D *>(objTrueJJ);
        
        // Convert TH1F to TH1D if needed
        if (!hRecoMB) {
          TH1F *hRecoFMB = dynamic_cast<TH1F *>(objRecoMB);
          if (hRecoFMB) {
            if (hRecoFMB->GetXaxis()->GetXbins() && hRecoFMB->GetXaxis()->GetXbins()->GetSize() > 0) {
              hRecoMB = new TH1D("hRecoMB_converted", hRecoFMB->GetTitle(),
                                hRecoFMB->GetNbinsX(), hRecoFMB->GetXaxis()->GetXbins()->GetArray());
            } else {
              hRecoMB = new TH1D("hRecoMB_converted", hRecoFMB->GetTitle(),
                                hRecoFMB->GetNbinsX(), hRecoFMB->GetXaxis()->GetXmin(), 
                                hRecoFMB->GetXaxis()->GetXmax());
            }
            for (int i = 1; i <= hRecoFMB->GetNbinsX(); ++i) {
              hRecoMB->SetBinContent(i, hRecoFMB->GetBinContent(i));
              hRecoMB->SetBinError(i, hRecoFMB->GetBinError(i));
            }
          }
        }
        if (!hRecoJJ) {
          TH1F *hRecoFJJ = dynamic_cast<TH1F *>(objRecoJJ);
          if (hRecoFJJ) {
            if (hRecoFJJ->GetXaxis()->GetXbins() && hRecoFJJ->GetXaxis()->GetXbins()->GetSize() > 0) {
              hRecoJJ = new TH1D("hRecoJJ_converted", hRecoFJJ->GetTitle(),
                                hRecoFJJ->GetNbinsX(), hRecoFJJ->GetXaxis()->GetXbins()->GetArray());
            } else {
              hRecoJJ = new TH1D("hRecoJJ_converted", hRecoFJJ->GetTitle(),
                                hRecoFJJ->GetNbinsX(), hRecoFJJ->GetXaxis()->GetXmin(), 
                                hRecoFJJ->GetXaxis()->GetXmax());
            }
            for (int i = 1; i <= hRecoFJJ->GetNbinsX(); ++i) {
              hRecoJJ->SetBinContent(i, hRecoFJJ->GetBinContent(i));
              hRecoJJ->SetBinError(i, hRecoFJJ->GetBinError(i));
            }
          }
        }
        if (!hTrueMB) {
          TH1F *hTrueFMB = dynamic_cast<TH1F *>(objTrueMB);
          if (hTrueFMB) {
            if (hTrueFMB->GetXaxis()->GetXbins() && hTrueFMB->GetXaxis()->GetXbins()->GetSize() > 0) {
              hTrueMB = new TH1D("hTrueMB_converted", hTrueFMB->GetTitle(),
                                hTrueFMB->GetNbinsX(), hTrueFMB->GetXaxis()->GetXbins()->GetArray());
            } else {
              hTrueMB = new TH1D("hTrueMB_converted", hTrueFMB->GetTitle(),
                                hTrueFMB->GetNbinsX(), hTrueFMB->GetXaxis()->GetXmin(), 
                                hTrueFMB->GetXaxis()->GetXmax());
            }
            for (int i = 1; i <= hTrueFMB->GetNbinsX(); ++i) {
              hTrueMB->SetBinContent(i, hTrueFMB->GetBinContent(i));
              hTrueMB->SetBinError(i, hTrueFMB->GetBinError(i));
            }
          }
        }
        if (!hTrueJJ) {
          TH1F *hTrueFJJ = dynamic_cast<TH1F *>(objTrueJJ);
          if (hTrueFJJ) {
            if (hTrueFJJ->GetXaxis()->GetXbins() && hTrueFJJ->GetXaxis()->GetXbins()->GetSize() > 0) {
              hTrueJJ = new TH1D("hTrueJJ_converted", hTrueFJJ->GetTitle(),
                                hTrueFJJ->GetNbinsX(), hTrueFJJ->GetXaxis()->GetXbins()->GetArray());
            } else {
              hTrueJJ = new TH1D("hTrueJJ_converted", hTrueFJJ->GetTitle(),
                                hTrueFJJ->GetNbinsX(), hTrueFJJ->GetXaxis()->GetXmin(), 
                                hTrueFJJ->GetXaxis()->GetXmax());
            }
            for (int i = 1; i <= hTrueFJJ->GetNbinsX(); ++i) {
              hTrueJJ->SetBinContent(i, hTrueFJJ->GetBinContent(i));
              hTrueJJ->SetBinError(i, hTrueFJJ->GetBinError(i));
            }
          }
        }
        
        if (hRecoMB && hRecoJJ && hTrueMB && hTrueJJ) {
          // Rebin if needed (same as LoadMergedResponseSample)
          const Double_t *ptbin = GetPtbin();
          const Double_t *ptbinGen = GetPtbinGen();
          Int_t nptBins = GetNptBins();
          Int_t nptBinsGen = GetNptBinsGen();
          
          if (GetREBINON()) {
            hRecoMB = dynamic_cast<TH1D *>(hRecoMB->Rebin(nptBins, "hRecoMB_rebinned", ptbin));
            hRecoJJ = dynamic_cast<TH1D *>(hRecoJJ->Rebin(nptBins, "hRecoJJ_rebinned", ptbin));
            hTrueMB = dynamic_cast<TH1D *>(hTrueMB->Rebin(nptBinsGen, "hTrueMB_rebinned", ptbinGen));
            hTrueJJ = dynamic_cast<TH1D *>(hTrueJJ->Rebin(nptBinsGen, "hTrueJJ_rebinned", ptbinGen));
          }
          
          // Sum MB and JJ for total
          hRecoTotal = dynamic_cast<TH1D *>(hRecoMB->Clone("hRecoTotal"));
          hRecoTotal->Add(hRecoJJ);
          hTrueTotal = dynamic_cast<TH1D *>(hTrueMB->Clone("hTrueTotal"));
          hTrueTotal->Add(hTrueJJ);
          hRecoTotal->SetDirectory(nullptr);
          hTrueTotal->SetDirectory(nullptr);
        }
      }
      fMB->Close();
      fJJ->Close();
      delete fMB;
      delete fJJ;
    }
  } else if (TString(outputPrefix) == "MBonly") {
    TFile *f = TFile::Open(kResponseFileMB, "READ");
    if (f && !f->IsZombie()) {
      TDirectory *dir = dynamic_cast<TDirectory *>(f->Get(kDirJets));
      if (dir) {
        TObject *objReco = dir->Get(kHistReco);
        TObject *objTrue = dir->Get(kHistTrue);
        TH1D *hRecoIn = dynamic_cast<TH1D *>(objReco);
        TH1D *hTrueIn = dynamic_cast<TH1D *>(objTrue);
        
        // Convert TH1F to TH1D if needed
        if (!hRecoIn) {
          TH1F *hRecoF = dynamic_cast<TH1F *>(objReco);
          if (hRecoF) {
            if (hRecoF->GetXaxis()->GetXbins() && hRecoF->GetXaxis()->GetXbins()->GetSize() > 0) {
              hRecoIn = new TH1D("hRecoMB_converted", hRecoF->GetTitle(),
                                hRecoF->GetNbinsX(), hRecoF->GetXaxis()->GetXbins()->GetArray());
            } else {
              hRecoIn = new TH1D("hRecoMB_converted", hRecoF->GetTitle(),
                                hRecoF->GetNbinsX(), hRecoF->GetXaxis()->GetXmin(), 
                                hRecoF->GetXaxis()->GetXmax());
            }
            for (int i = 1; i <= hRecoF->GetNbinsX(); ++i) {
              hRecoIn->SetBinContent(i, hRecoF->GetBinContent(i));
              hRecoIn->SetBinError(i, hRecoF->GetBinError(i));
            }
          }
        }
        if (!hTrueIn) {
          TH1F *hTrueF = dynamic_cast<TH1F *>(objTrue);
          if (hTrueF) {
            if (hTrueF->GetXaxis()->GetXbins() && hTrueF->GetXaxis()->GetXbins()->GetSize() > 0) {
              hTrueIn = new TH1D("hTrueMB_converted", hTrueF->GetTitle(),
                                hTrueF->GetNbinsX(), hTrueF->GetXaxis()->GetXbins()->GetArray());
            } else {
              hTrueIn = new TH1D("hTrueMB_converted", hTrueF->GetTitle(),
                                hTrueF->GetNbinsX(), hTrueF->GetXaxis()->GetXmin(), 
                                hTrueF->GetXaxis()->GetXmax());
            }
            for (int i = 1; i <= hTrueF->GetNbinsX(); ++i) {
              hTrueIn->SetBinContent(i, hTrueF->GetBinContent(i));
              hTrueIn->SetBinError(i, hTrueF->GetBinError(i));
            }
          }
        }
        
        if (hRecoIn && hTrueIn) {
          // Rebin if needed (same as LoadResponseSampleMB)
          const Double_t *ptbin = GetPtbin();
          const Double_t *ptbinGen = GetPtbinGen();
          Int_t nptBins = GetNptBins();
          Int_t nptBinsGen = GetNptBinsGen();
          
          if (GetREBINON()) {
            hRecoTotal = dynamic_cast<TH1D *>(hRecoIn->Rebin(nptBins, "hRecoTotal", ptbin));
            hTrueTotal = dynamic_cast<TH1D *>(hTrueIn->Rebin(nptBinsGen, "hTrueTotal", ptbinGen));
          } else {
            hRecoTotal = dynamic_cast<TH1D *>(hRecoIn->Clone("hRecoTotal"));
            hTrueTotal = dynamic_cast<TH1D *>(hTrueIn->Clone("hTrueTotal"));
          }
          hRecoTotal->SetDirectory(nullptr);
          hTrueTotal->SetDirectory(nullptr);
        }
      }
      f->Close();
      delete f;
    }
  } else if (TString(outputPrefix) == "JJonly") {
    TFile *f = TFile::Open(kResponseFileJJ, "READ");
    if (f && !f->IsZombie()) {
      TDirectory *dir = dynamic_cast<TDirectory *>(f->Get(kDirJetsJJ));
      if (dir) {
        // Load directly as TH1 (same as DrawJetsMCfJetQA.h line 210)
        // This allows TH1F or TH1D to be handled the same way
        TH1 *hRecoIn = (TH1 *)dir->Get(kHistReco);
        TH1 *hTrueIn = (TH1 *)dir->Get(kHistTrue);
        
        if (hRecoIn && hTrueIn) {
          // Rebin if needed (same as DrawJetsMCfJetQA.h line 211-213)
          const Double_t *ptbin = GetPtbin();
          const Double_t *ptbinGen = GetPtbinGen();
          Int_t nptBins = GetNptBins();
          Int_t nptBinsGen = GetNptBinsGen();
          
          if (GetREBINON()) {
            // Rebin directly (same as DrawJetsMCfJetQA.h)
            TH1 *hRecoRebinned = hRecoIn->Rebin(nptBins, "hRecoTotal", ptbin);
            TH1 *hTrueRebinned = hTrueIn->Rebin(nptBinsGen, "hTrueTotal", ptbinGen);
            // Convert to TH1D for consistency
            hRecoTotal = dynamic_cast<TH1D *>(hRecoRebinned);
            if (!hRecoTotal) {
              hRecoTotal = new TH1D("hRecoTotal", hRecoRebinned->GetTitle(),
                                   nptBins, ptbin);
              for (int i = 1; i <= hRecoRebinned->GetNbinsX(); ++i) {
                hRecoTotal->SetBinContent(i, hRecoRebinned->GetBinContent(i));
                hRecoTotal->SetBinError(i, hRecoRebinned->GetBinError(i));
              }
              delete hRecoRebinned;
            }
            hTrueTotal = dynamic_cast<TH1D *>(hTrueRebinned);
            if (!hTrueTotal) {
              hTrueTotal = new TH1D("hTrueTotal", hTrueRebinned->GetTitle(),
                                   nptBinsGen, ptbinGen);
              for (int i = 1; i <= hTrueRebinned->GetNbinsX(); ++i) {
                hTrueTotal->SetBinContent(i, hTrueRebinned->GetBinContent(i));
                hTrueTotal->SetBinError(i, hTrueRebinned->GetBinError(i));
              }
              delete hTrueRebinned;
            }
          } else {
            // Clone and convert to TH1D
            hRecoTotal = dynamic_cast<TH1D *>(hRecoIn->Clone("hRecoTotal"));
            if (!hRecoTotal) {
              hRecoTotal = new TH1D("hRecoTotal", hRecoIn->GetTitle(),
                                   hRecoIn->GetNbinsX(), hRecoIn->GetXaxis()->GetXmin(),
                                   hRecoIn->GetXaxis()->GetXmax());
              for (int i = 1; i <= hRecoIn->GetNbinsX(); ++i) {
                hRecoTotal->SetBinContent(i, hRecoIn->GetBinContent(i));
                hRecoTotal->SetBinError(i, hRecoIn->GetBinError(i));
              }
            }
            hTrueTotal = dynamic_cast<TH1D *>(hTrueIn->Clone("hTrueTotal"));
            if (!hTrueTotal) {
              hTrueTotal = new TH1D("hTrueTotal", hTrueIn->GetTitle(),
                                   hTrueIn->GetNbinsX(), hTrueIn->GetXaxis()->GetXmin(),
                                   hTrueIn->GetXaxis()->GetXmax());
              for (int i = 1; i <= hTrueIn->GetNbinsX(); ++i) {
                hTrueTotal->SetBinContent(i, hTrueIn->GetBinContent(i));
                hTrueTotal->SetBinError(i, hTrueIn->GetBinError(i));
              }
            }
          }
          hRecoTotal->SetDirectory(nullptr);
          hTrueTotal->SetDirectory(nullptr);
        }
      }
      f->Close();
      delete f;
    }
  }
  
  if (!hRecoTotal || !hTrueTotal) {
    std::cerr << "[Error] Failed to load original MC histograms for QA plots" << std::endl;
    delete hRecoMatched;
    delete hTrueMatched;
    return;
  }
  
  // Calculate Purity: matched reco / total reco (using original h_jet_pt, not fake/miss included)
  TH1D *hPurity = dynamic_cast<TH1D *>(hRecoMatched->Clone("hPurity"));
  hPurity->SetDirectory(nullptr);
  hPurity->Divide(hRecoMatched, hRecoTotal, 1., 1., "B");
  hPurity->SetStats(0);
  hPurity->SetTitle("");
  
  // Calculate Efficiency: matched truth / total truth (using original h_jet_pt_part, not fake/miss included)
  TH1D *hEfficiency = dynamic_cast<TH1D *>(hTrueMatched->Clone("hEfficiency"));
  hEfficiency->SetDirectory(nullptr);
  hEfficiency->Divide(hTrueMatched, hTrueTotal, 1., 1., "B");
  hEfficiency->SetStats(0);
  hEfficiency->SetTitle("");
  
  // Draw Purity
  TCanvas *canPurity = new TCanvas(Form("canPurity_%s", outputPrefix), "Jet Purity", 800, 600);
  setpad(canPurity, 0.02, 0.15, 0.15, 0.05);
  canPurity->Draw();
  canPurity->cd();
  gPad->SetGrid(0);
  hPurity->SetMarkerStyle(20);
  hPurity->SetMarkerColor(kRed);
  hPurity->SetLineColor(kRed);
  hPurity->SetLineWidth(2);
  hPurity->GetXaxis()->SetRangeUser(5, 140);
  hPurity->GetYaxis()->SetRangeUser(0.0, 1.1);
  hset(*hPurity, "#it{p}_{T, jet}^{reco} (GeV/#it{c})", "Jet purity", 1.3, 1.0, 0.05, 0.07, 0.01, 0.01, 0.05, 0.05, 510, 510);
  hPurity->Draw("pe");
  canPurity->Print(Form("%s/Purity_%s.pdf", outputDir.Data(), outputPrefix));
  delete canPurity;
  
  // Draw Efficiency
  TCanvas *canEfficiency = new TCanvas(Form("canEfficiency_%s", outputPrefix), "Jet Efficiency", 800, 600);
  setpad(canEfficiency, 0.02, 0.15, 0.15, 0.05);
  canEfficiency->Draw();
  canEfficiency->cd();
  gPad->SetGrid(0);
  hEfficiency->SetMarkerStyle(20);
  hEfficiency->SetMarkerColor(kBlue);
  hEfficiency->SetLineColor(kBlue);
  hEfficiency->SetLineWidth(2);
  hEfficiency->GetXaxis()->SetRangeUser(5, 140);
  hEfficiency->GetYaxis()->SetRangeUser(0.0, 1.1);
  hset(*hEfficiency, "#it{p}_{T, jet}^{true} (GeV/#it{c})", "#it{#varepsilon}_{reco}^{jet}", 1.3, 1.0, 0.05, 0.07, 0.01, 0.01, 0.05, 0.05, 510, 510);
  hEfficiency->Draw("pe");
  canEfficiency->Print(Form("%s/Efficiency_%s.pdf", outputDir.Data(), outputPrefix));
  delete canEfficiency;
  
  // Draw Response Matrix (use h2HCorrelate instead of hResp2D)
  TCanvas *canResp2D = new TCanvas(Form("canResp2D_%s", outputPrefix), "Response Matrix", 800, 600);
  canResp2D->SetLogz(1);
  setpad(canResp2D, 0.02, 0.15, 0.15, 0.15);
  canResp2D->Draw();
  canResp2D->cd();
  TH2D *hResp2DClone = dynamic_cast<TH2D *>(h2HCorrelate->Clone("hResp2DClone"));
  hResp2DClone->SetDirectory(nullptr);
  hResp2DClone->SetStats(0);
  hResp2DClone->SetTitle("");
  hset(*hResp2DClone, "#it{p}_{T, jet}^{reco} (GeV/#it{c})", "#it{p}_{T, jet}^{true} (GeV/#it{c})",
       1.2, 1.3, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 505);
  hResp2DClone->GetZaxis()->SetRangeUser(1e-0, 1e8);
  hResp2DClone->Draw("colz");
  canResp2D->Print(Form("%s/ResponseMatrix_%s.pdf", outputDir.Data(), outputPrefix));
  delete canResp2D;
  delete hResp2DClone;
  
  delete hRecoMatched;
  delete hTrueMatched;
  delete hPurity;
  delete hEfficiency;
  delete hRecoTotal;
  delete hTrueTotal;
  delete h2HCorrelate;
  
  std::cerr << "[Info] QA plots saved: " << outputPrefix << std::endl;
}

// ============================================================
// Helper function: Unfold and plot cross section for one response sample
// ============================================================

void UnfoldAndPlotCrossSection(TH1D *hTrue, TH1D *hReco, TH2D *hResp2D, 
                                TH1D *hRecoData, Double_t Ntvx,
                                const char *plotTitle, const char *outputFileName) {
  std::cerr << "========================================" << std::endl;
  std::cerr << "[Info] " << plotTitle << std::endl;
  std::cerr << "========================================" << std::endl;
  
  // Debug: Check input histograms before building response
  std::cerr << "[Debug] UnfoldAndPlotCrossSection input check:" << std::endl;
  std::cerr << "  hTrue integral: " << hTrue->Integral() << std::endl;
  std::cerr << "  hReco integral: " << hReco->Integral() << std::endl;
  std::cerr << "  hResp2D integral: " << hResp2D->Integral() << std::endl;
  std::cerr << "  hRecoData integral: " << hRecoData->Integral() << std::endl;
  
  TH1D *hTrueProjCheck = dynamic_cast<TH1D *>(hResp2D->ProjectionY("hTrueProj_input_check", 1, hResp2D->GetNbinsX(), "e"));
  TH1D *hRecoProjCheck = dynamic_cast<TH1D *>(hResp2D->ProjectionX("hRecoProj_input_check", 1, hResp2D->GetNbinsY(), "e"));
  std::cerr << "  hResp2D True projection: " << hTrueProjCheck->Integral() << std::endl;
  std::cerr << "  hResp2D Reco projection: " << hRecoProjCheck->Integral() << std::endl;
  std::cerr << "  hTrue / hTrueProj ratio: " << (hTrueProjCheck->Integral() > 0 ? hTrue->Integral() / hTrueProjCheck->Integral() : 0) << std::endl;
  std::cerr << "  hReco / hRecoProj ratio: " << (hRecoProjCheck->Integral() > 0 ? hReco->Integral() / hRecoProjCheck->Integral() : 0) << std::endl;
  delete hTrueProjCheck;
  delete hRecoProjCheck;
  
  // Build response matrix (kernel-style)
  double responseMatrixScale = 1.0;
  TH1D *hPurity = nullptr;
  TH1D *hEfficiency = nullptr;
  RooUnfoldResponse *Response = BuildResponse(hTrue, hReco, hResp2D, responseMatrixScale, hPurity, hEfficiency);

  if (!hPurity || !hEfficiency) {
    std::cerr << "[Error] Failed to compute purity/efficiency!" << std::endl;
    delete Response;
    return;
  }

  // KERNEL-STYLE: Apply purity correction to data before unfolding
  TH1D *hRecoDataMatched = dynamic_cast<TH1D *>(hRecoData->Clone("hRecoDataMatched"));
  hRecoDataMatched->SetDirectory(nullptr);
  for (int i = 1; i <= hRecoDataMatched->GetNbinsX(); ++i) {
    const double x = hRecoData->GetBinContent(i);
    const double ex = hRecoData->GetBinError(i);
    const double p = hPurity->GetBinContent(i);
    hRecoDataMatched->SetBinContent(i, x * p);
    hRecoDataMatched->SetBinError(i, ex * p);
  }
  std::cerr << "[Info] Applied purity correction to data (kernel-style)" << std::endl;

  // Determine optimal SVD k via d-vector + stability validation
  Int_t optimalSvdK = FindOptimalSvdK_Dvector(Response, hRecoDataMatched);
  std::cerr << "[Info] Performing SVD unfolding with k=" << optimalSvdK << std::endl;
  RooUnfoldSvd unfold(Response, hRecoDataMatched, optimalSvdK);
  
  // Check if unfolding produced a valid result
  TH1D *hUnfoldedCheck = dynamic_cast<TH1D *>(unfold.Hreco());
  if (!hUnfoldedCheck) {
    std::cerr << "[Error] SVD unfolding failed - result is null!" << std::endl;
    delete Response;
    return;
  }
  
  // Check for NaN/Inf values
  int nNaN = 0;
  int nInf = 0;
  for (int i = 1; i <= hUnfoldedCheck->GetNbinsX(); ++i) {
    double cont = hUnfoldedCheck->GetBinContent(i);
    if (TMath::IsNaN(cont)) nNaN++;
    if (!TMath::Finite(cont)) nInf++;
  }
  
  if (nNaN > 0 || nInf > 0) {
    std::cerr << "[Warning] Unfolded result contains invalid values: NaN=" << nNaN 
              << ", Inf=" << nInf << std::endl;
  } else {
    std::cerr << "[Info] SVD unfolding completed successfully" << std::endl;
  }
  
  TH1D *hUnfoldedDetected = dynamic_cast<TH1D *>(unfold.Hreco());
  if (!hUnfoldedDetected) {
    std::cerr << "[Error] Unfolding failed - hUnfoldedDetected is null!" << std::endl;
    delete Response;
    return;
  }
  hUnfoldedDetected->SetDirectory(nullptr);

  // KERNEL-STYLE: Apply efficiency correction after unfolding
  TH1D *hUnfolded = dynamic_cast<TH1D *>(hUnfoldedDetected->Clone("hUnfolded_effCorrected"));
  hUnfolded->SetDirectory(nullptr);
  for (int i = 1; i <= hUnfolded->GetNbinsX(); ++i) {
    const double x = hUnfoldedDetected->GetBinContent(i);
    const double ex = hUnfoldedDetected->GetBinError(i);
    const double pt = hUnfolded->GetXaxis()->GetBinCenter(i);
    const int effBin = hEfficiency->GetXaxis()->FindBin(pt);
    double eff = 0.0;
    if (effBin >= 1 && effBin <= hEfficiency->GetNbinsX()) {
      eff = hEfficiency->GetBinContent(effBin);
    }
    if (eff > 0) {
      hUnfolded->SetBinContent(i, x / eff);
      hUnfolded->SetBinError(i, ex / eff);
    } else {
      hUnfolded->SetBinContent(i, 0.0);
      hUnfolded->SetBinError(i, 0.0);
    }
  }
  std::cerr << "[Info] Applied efficiency correction to unfolded result (kernel-style)" << std::endl;
  
  // Check unfolded vs data ratio
  double dataIntegral = hRecoData->Integral();
  double unfoldedIntegral = hUnfolded->Integral();
  if (dataIntegral > 0) {
    double ratio = unfoldedIntegral / dataIntegral;
    std::cerr << "[Info] Unfolded/Data integral ratio: " << ratio << std::endl;
  
  // Debug: Detailed unfolding check
  std::cerr << "[Debug] Unfolding Quality Check:" << std::endl;
  std::cerr << "  Data integral: " << dataIntegral << std::endl;
  std::cerr << "  Unfolded integral: " << unfoldedIntegral << std::endl;
  if (dataIntegral > 0) {
    double ratio = unfoldedIntegral / dataIntegral;
    std::cerr << "  Ratio (Unfolded/Data): " << ratio << std::endl;
    if (ratio < 0.5 || ratio > 2.0) {
      std::cerr << "  [WARNING] Unusual ratio detected! Expected ~1.0" << std::endl;
    }
  }
  
  // Check bin-by-bin consistency
  int nZeroData = 0, nZeroUnfolded = 0;
  int nLargeDiff = 0;
  for (int ibin = 1; ibin <= hUnfolded->GetNbinsX(); ++ibin) {
    double dataCont = hRecoData->GetBinContent(ibin);
    double unfoldCont = hUnfolded->GetBinContent(ibin);
    if (dataCont == 0) nZeroData++;
    if (unfoldCont == 0) nZeroUnfolded++;
    if (dataCont > 0 && unfoldCont > 0) {
      double diff = TMath::Abs(unfoldCont - dataCont) / dataCont;
      if (diff > 0.5) nLargeDiff++;
    }
  }
  std::cerr << "  Zero bins: Data=" << nZeroData << ", Unfolded=" << nZeroUnfolded << std::endl;
  std::cerr << "  Bins with >50% difference: " << nLargeDiff << " / " << hUnfolded->GetNbinsX() << std::endl;
  
  }
  
  // Calculate errors from covariance matrix
  TMatrixD covMatrix = unfold.Ereco();
  int nBins = hUnfolded->GetNbinsX();
  bool useCovMatrix = (covMatrix.GetNrows() >= nBins && covMatrix.GetNcols() >= nBins);
  
  for (int i = 1; i <= hUnfolded->GetNbinsX(); ++i) {
    double err = hUnfolded->GetBinError(i);
    double cont = hUnfolded->GetBinContent(i);
    double absCont = TMath::Abs(cont);
    
    int matrixIdx = i - 1;
    double newErr = err;
    if (useCovMatrix && matrixIdx < covMatrix.GetNrows() && matrixIdx < covMatrix.GetNcols()) {
      double covDiag = covMatrix(matrixIdx, matrixIdx);
      if (covDiag > 0) {
        newErr = TMath::Sqrt(covDiag);
      }
    }
    
    // Check for invalid errors and fix
    if (TMath::IsNaN(newErr) || !TMath::Finite(newErr)) {
      if (absCont > 0) {
        newErr = TMath::Sqrt(absCont);
      } else {
        newErr = 0.;
      }
    } else if (newErr == 0. && absCont > 0) {
      newErr = TMath::Sqrt(absCont);
    } else if (newErr > 10. * absCont || newErr > 1e10) {
      double ptCenter = hUnfolded->GetXaxis()->GetBinCenter(i);
      int dataBin = hRecoData->GetXaxis()->FindBin(ptCenter);
      if (dataBin >= 1 && dataBin <= hRecoData->GetNbinsX()) {
        double dataContent = hRecoData->GetBinContent(dataBin);
        double dataError = hRecoData->GetBinError(dataBin);
        if (dataContent > 0 && dataError > 0 && absCont > 0) {
          newErr = dataError * TMath::Sqrt(absCont / dataContent);
          if (newErr > 0.5 * absCont) {
            newErr = 0.5 * absCont;
          }
        }
      }
      if (newErr <= 0 && absCont > 0) {
        newErr = TMath::Sqrt(absCont);
      }
      if (absCont > 0) {
        double minErr = 0.1 * TMath::Sqrt(absCont);
        if (newErr < minErr) {
          newErr = minErr;
        }
      }
    }
    
    hUnfolded->SetBinError(i, newErr);
  }
  
  // Calculate cross section
  TH1D *UnfoldData = dynamic_cast<TH1D *>(hUnfolded->Clone("Run3_CrossSection"));
  
  // Calculate integrated luminosity
  double L_int_mb = Ntvx / 53.0;  // mb^-1
  double L_int_pb = L_int_mb / 1e9;  // pb^-1 (1 mb^-1 = 10^9 pb^-1)
  
  // Normalization factor: (53 / Ntvx) / (efficiency factors)
  // Efficiency factors: 0.95608 (vertex), 0.835 (tracking), 0.96 (trigger)
//   double normFactor = 53. / Ntvx / 0.95608 / 0.835 / 0.96;
  double normFactor = 53. / Ntvx / 0.95608 / 0.835 / 0.96;
  
  // Convert unfolded counts to cross section: cross section = (counts * normFactor) / binWidth
  // Scale() with "width" option automatically divides by bin width for each bin
  UnfoldData->Scale(normFactor, "width");
  
  // Debug: Cross section calculation check
  std::cerr << "[Debug] Cross Section Calculation:" << std::endl;
  std::cerr << "  Ntvx: " << Ntvx << std::endl;
  std::cerr << "  L_int (pb^-1): " << L_int_pb << std::endl;
  std::cerr << "  normFactor: " << normFactor << std::endl;
  
  double crossSectionIntegral = UnfoldData->Integral("width");
  double unfoldedCountsIntegral = hUnfolded->Integral();
  std::cerr << "  Unfolded counts integral: " << unfoldedCountsIntegral << std::endl;
  std::cerr << "  Cross section integral (width-weighted): " << crossSectionIntegral << std::endl;
  
  // Check a few sample bins
  std::cerr << "  Sample bins (pT_low | pT_high | counts | cross_section):" << std::endl;
  for (int ibin = 1; ibin <= TMath::Min(5, UnfoldData->GetNbinsX()); ++ibin) {
    double pTLow = UnfoldData->GetXaxis()->GetBinLowEdge(ibin);
    double pTHigh = UnfoldData->GetXaxis()->GetBinUpEdge(ibin);
    double binWidth = UnfoldData->GetXaxis()->GetBinWidth(ibin);
    double counts = hUnfolded->GetBinContent(ibin);
    double crossSection = UnfoldData->GetBinContent(ibin);
    std::cerr << "    " << pTLow << " | " << pTHigh << " | " << counts << " | " << crossSection << std::endl;
  }
  
  
  // Plotting
  Filipad2 *CrossSectionPad = new Filipad2(Form("Cross Section - %s", plotTitle), ++GetNN(), 2, 0.3, 100, 50, 0.7, 1, 1);
  CrossSectionPad->Draw();
  TPad *unfoldp = CrossSectionPad->GetPad(1);
  optFili(*unfoldp, 0, 1, 0, 1);
  TPad *ratunfoldp = CrossSectionPad->GetPad(2);
  optFili(*ratunfoldp, 0, 0, 0, 0);
  
  TLegend *unfoldlegend = new TLegend(0.543062, 0.378882, 0.80622, 0.575155, NULL, "brNDC");
  unfoldlegend->SetTextSize(0.04);
  unfoldlegend->SetBorderSize(0);
  unfoldlegend->SetFillColorAlpha(0, 0);
  
  unfoldp->cd();
  gPad->SetTicks(1, 1);
  unfoldlegend->AddEntry(UnfoldData, "Run 3,#kern[-0.7]{ }#sqrt{#it{s}}=13.6 TeV");
  hset(*UnfoldData, GetJetPtGenTitleX(), "d^{2}#sigma/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]", 
       0.9, 1.5, 0.05, 0.045, 0.01, 0.005, 0.045, 0.040, 510, 1005);
  hoptset(*UnfoldData, 0, kRed, GetPlotPtMin(), GetPlotPtMax(), 
          3e-7, 5e-1, 0.6, 1, 2, 24);
  UnfoldData->Draw("pe");
  
  // Add Run 2 data
  auto graphs = Run2Data_manual();
  auto grun2Syst = (TGraphErrors *) graphs[1]->Clone("grun2Syst");
  auto hrun2Syst = GraphToHistogram(grun2Syst);
  auto grun2Stat = (TGraphErrors *) graphs[0]->Clone("grun2Stat");
  auto hrun2Stat = GraphToHistogram(grun2Stat);
  
  hoptset(*hrun2Syst, 0, kBlue, GetPlotPtMin(), GetPlotPtMax(), 3e-7, 5e-1, 0.6, 1, 2, 25);
  hrun2Syst->SetFillColorAlpha(kBlue, 0.3);
  hrun2Syst->SetMarkerStyle(0);
  hrun2Syst->Draw("pe2same");
  
  hoptset(*hrun2Stat, 0, kBlue, GetPlotPtMin(), GetPlotPtMax(), 3e-7, 5e-1, 0.6, 1, 2, 25);
  hrun2Stat->SetMarkerStyle(25);
  hrun2Stat->SetMarkerColor(kBlue);
  hrun2Stat->SetLineColor(kBlue);
  hrun2Stat->Draw("e same");
  
  unfoldlegend->AddEntry(hrun2Syst, "Run 2,#kern[-0.7]{ }#sqrt{#it{s}}=13 TeV");
  
  // Add PYTHIA 13600
  TH1 *hPYTHIA13600 = DrawPythia13600();
  if (hPYTHIA13600) {
    hoptset(*hPYTHIA13600, 0, kCyan, GetPlotPtMin(), GetPlotPtMax(), 3e-7, 5e-1, 0.5, 1, 2, 24);
    hset(*hPYTHIA13600, GetJetPtGenTitleX(), "d^{2}#sigma/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]", 
         0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hPYTHIA13600->Draw("pesame");
    unfoldlegend->AddEntry(hPYTHIA13600, "PYTHIA8 13.6 TeV", "lep");
  }
  
  // Add Changwhan JJ MC
  TH1 *hChangwhanJJMC = DrawChangwhan2022JJMC();
  if (hChangwhanJJMC) {
    hoptset(*hChangwhanJJMC, 0, kGreen+2, GetPlotPtMin(), GetPlotPtMax(), 3e-7, 5e-1, 0.5, 1, 2, 26);
    hset(*hChangwhanJJMC, GetJetPtGenTitleX(), "d^{2}#sigma/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]", 
         0.9, 1.4, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hChangwhanJJMC->Draw("pesame");
    unfoldlegend->AddEntry(hChangwhanJJMC, "Changwhan JJ MC", "lep");
  }
  
  unfoldlegend->Draw();
  ALICEfigureLegend("ALICE WIP", 0.131403, 0.698065, 0.380846, 0.948387, 
                    0.482183, 0.692903, 0.732739, 0.883871);
  
  // Ratio pad
  ratunfoldp->cd();
  gPad->SetTicks(1, 1);
  
  TH1F *frame = new TH1F(Form("hframe_ratio_%s", outputFileName), "", 100, GetPlotPtMin(), GetPlotPtMax());
  frame->SetDirectory(nullptr);
  frame->SetMinimum(0.8);
  frame->SetMaximum(1.2);
  hset(*frame, GetJetPtGenTitleX(), "Comp. / Run 3", 1.2, 0.75, 0.1, 0.09, 0.01, 0.01,
       0.1, 0.1, 510, 505);
  frame->SetLineColor(0);
  frame->SetMarkerSize(0);
  frame->Draw("AXIS");
  
  TLine *line = new TLine(GetPlotPtMin(), 1.0, GetPlotPtMax(), 1.0);
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->Draw("lsame");
  
  // Run 2 / Run 3 ratio
  auto ratRun2Run3data = DrawRatioTH1(hrun2Syst, UnfoldData);
  hoptset(*ratRun2Run3data, 0, kBlue, GetPlotPtMin(), GetPlotPtMax(), 0.4, 2.7, 0.5, 1, 2, 25);
  ratRun2Run3data->SetFillColorAlpha(kBlue, 0.3);
  ratRun2Run3data->SetMarkerStyle(0);
  ratRun2Run3data->Draw("pe2same");
  
  auto ratRun2Run3dataStat = DrawRatioTH1(hrun2Stat, UnfoldData);
  hoptset(*ratRun2Run3dataStat, 0, kBlue, GetPlotPtMin(), GetPlotPtMax(), 0.4, 2.7, 0.5, 1, 2, 25);
  ratRun2Run3dataStat->SetMarkerStyle(25);
  ratRun2Run3dataStat->SetMarkerColor(kBlue);
  ratRun2Run3dataStat->SetLineColor(kBlue);
  ratRun2Run3dataStat->Draw("e same");
  
  // PYTHIA 13600 / Run 3 ratio
  TH1 *ratPYTHIA13600 = nullptr;
  if (hPYTHIA13600) {
    ratPYTHIA13600 = DrawRatioTH1(hPYTHIA13600, UnfoldData);
    hoptset(*ratPYTHIA13600, 0, kCyan, GetPlotPtMin(), GetPlotPtMax(), 0.4, 2.7, 0.4, 1, 2, 24);
    ratPYTHIA13600->Draw("pesame");
  }
  
  // Changwhan JJ MC / Run 3 ratio
  TH1 *ratChangwhanJJMC = nullptr;
  if (hChangwhanJJMC) {
    ratChangwhanJJMC = DrawRatioTH1(hChangwhanJJMC, UnfoldData);
    hoptset(*ratChangwhanJJMC, 0, kGreen+2, GetPlotPtMin(), GetPlotPtMax(), 0.4, 2.7, 0.4, 1, 2, 26);
    ratChangwhanJJMC->Draw("pesame");
  }
  
  // Save plot
  if (GetDRAWPLOTS()) {
    CrossSectionPad->C->SaveAs(Form("%s/%s", CrossSectionDirName.Data(), outputFileName));
  }
  
  std::cerr << "[Info] Cross section plot saved: " << outputFileName << std::endl;
  
  // Cleanup
  delete Response;
  delete UnfoldData;
  delete hrun2Syst;
  delete hrun2Stat;
  delete grun2Syst;
  delete grun2Stat;
  delete frame;
  delete line;
  if (ratRun2Run3data) delete ratRun2Run3data;
  if (ratRun2Run3dataStat) delete ratRun2Run3dataStat;
  if (ratPYTHIA13600) delete ratPYTHIA13600;
  if (ratChangwhanJJMC) delete ratChangwhanJJMC;
}

// ============================================================
// Main function
// ============================================================

void DrawDataCrossSection() {
  std::cerr << "========================================" << std::endl;
  std::cerr << "Draw Data Cross Section" << std::endl;
  std::cerr << "Creating 3 plots:" << std::endl;
  std::cerr << "  1. Merged MB+JJ MC RM" << std::endl;
  std::cerr << "  2. MB MC only RM" << std::endl;
  std::cerr << "  3. JJ MC only RM" << std::endl;
  std::cerr << "SVD k: determined per response matrix (d-vector + stability)" << std::endl;
  std::cerr << "========================================" << std::endl;
  
  // Create output directory
  gSystem->Exec(Form("mkdir -p %s", CrossSectionDirName.Data()));
  
  // Load data (same for all 3 plots)
  TH1D *hRecoData = nullptr;
  Double_t nEventsData = 0.;
  Double_t Ntvx = 0.;
  
  if (!LoadData(hRecoData, nEventsData, Ntvx)) {
    std::cerr << "[Error] Failed to load data" << std::endl;
    return;
  }
  
  // Clone data histogram for each plot (since it might be modified)
  TH1D *hRecoDataMerged = dynamic_cast<TH1D *>(hRecoData->Clone("hRecoDataMerged"));
  TH1D *hRecoDataMB = dynamic_cast<TH1D *>(hRecoData->Clone("hRecoDataMB"));
  TH1D *hRecoDataJJ = dynamic_cast<TH1D *>(hRecoData->Clone("hRecoDataJJ"));
  hRecoDataMerged->SetDirectory(nullptr);
  hRecoDataMB->SetDirectory(nullptr);
  hRecoDataJJ->SetDirectory(nullptr);
  
  // ============================================================
  // Plot 1: Merged MB+JJ MC Response Matrix
  // ============================================================
  {
    TH1D *hTrue = nullptr;
    TH1D *hReco = nullptr;
    TH2D *hResp2D = nullptr;
    Double_t nEventsMB = 0.;
    Double_t nEventsJJ = 0.;
    
    if (!LoadMergedResponseSample(hTrue, hReco, hResp2D, nEventsMB, nEventsJJ)) {
      std::cerr << "[Error] Failed to load merged response sample" << std::endl;
    } else {
      // Draw QA plots
      DrawQAPlots(hTrue, hReco, hResp2D, "Merged MB+JJ MC RM", "MergedMBJJ");
      
      UnfoldAndPlotCrossSection(hTrue, hReco, hResp2D, hRecoDataMerged, Ntvx,
                                "Merged MB+JJ MC RM",
                                "CrossSection_MergedMBJJ.pdf");
      
      delete hTrue;
      delete hReco;
      delete hResp2D;
    }
  }
  
  // ============================================================
  // Plot 2: MB MC only Response Matrix
  // ============================================================
  {
    TH1D *hTrue = nullptr;
    TH1D *hReco = nullptr;
    TH2D *hResp2D = nullptr;
    Double_t nEventsResp = 0.;
    
    if (!LoadResponseSampleMB(hTrue, hReco, hResp2D, nEventsResp)) {
      std::cerr << "[Error] Failed to load MB MC response sample" << std::endl;
    } else {
      // Draw QA plots
      DrawQAPlots(hTrue, hReco, hResp2D, "MB MC only RM", "MBonly");
      
      UnfoldAndPlotCrossSection(hTrue, hReco, hResp2D, hRecoDataMB, Ntvx,
                                "MB MC only RM",
                                "CrossSection_MBonly.pdf");
      
      delete hTrue;
      delete hReco;
      delete hResp2D;
    }
  }
  
  // ============================================================
  // Plot 3: JJ MC only Response Matrix
  // ============================================================
  {
    TH1D *hTrue = nullptr;
    TH1D *hReco = nullptr;
    TH2D *hResp2D = nullptr;
    Double_t nEventsResp = 0.;
    
    if (!LoadResponseSampleJJ(hTrue, hReco, hResp2D, nEventsResp)) {
      std::cerr << "[Error] Failed to load JJ MC response sample" << std::endl;
    } else {
      // Draw QA plots
      DrawQAPlots(hTrue, hReco, hResp2D, "JJ MC only RM", "JJonly");
      
      UnfoldAndPlotCrossSection(hTrue, hReco, hResp2D, hRecoDataJJ, Ntvx,
                                "JJ MC only RM",
                                "CrossSection_JJonly.pdf");
      
      delete hTrue;
      delete hReco;
      delete hResp2D;
    }
  }
  
  // Cleanup
  delete hRecoData;
  delete hRecoDataMerged;
  delete hRecoDataMB;
  delete hRecoDataJJ;
  
  std::cerr << "========================================" << std::endl;
  std::cerr << "[Info] All 3 cross section plots completed!" << std::endl;
  std::cerr << "========================================" << std::endl;
}
