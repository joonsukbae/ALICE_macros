// Convert V0 mass resolution (K0s) to track pT resolution
// Based on approximation:
//   pT_track ≈ (1/2) × pT_K0s  (2-body symmetric decay)
//   σ_pT/pT ≈ C × σ_M/M_K0s    (C ≈ √2 ≈ 1.414)
//
// Usage:
//   .x ConvertV0MassToTrackPtResolution.C("plots/V0MassFits/V0MassFitResults.root", "plots/TrackPtResolution")

#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TMath.h>
#include <TString.h>

void ConvertV0MassToTrackPtResolutions(const char* inputFile = "plots/V0MassFits/V0MassFitResults.root",
                                      const char* outputDir = "plots/TrackPtResolution",
                                      double conversionFactor = TMath::Sqrt(2.0),  // C ≈ √2 ≈ 1.414
                                      double pTScaleFactor = 0.5,  // track pT ≈ 0.5 × K0s pT
                                      double k0sMass = 0.4976) {  // K0s mass in GeV/c^2
  
  // Open input file
  TFile* fIn = TFile::Open(inputFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    printf("[Error] Cannot open input file: %s\n", inputFile);
    return;
  }

  // Load graphs from file
  TGraphErrors* gRelSigma = dynamic_cast<TGraphErrors*>(fIn->Get("gRelSigmaVsPt"));
  TGraphErrors* gRelSigmaMC = dynamic_cast<TGraphErrors*>(fIn->Get("gRelSigmaVsPtMC"));
  
  if (!gRelSigma) {
    printf("[Error] Graph 'gRelSigmaVsPt' not found in file: %s\n", inputFile);
    fIn->Close();
    return;
  }

  printf("[Info] Loaded gRelSigmaVsPt with %d points\n", gRelSigma->GetN());
  if (gRelSigmaMC) {
    printf("[Info] Loaded gRelSigmaVsPtMC with %d points\n", gRelSigmaMC->GetN());
  }

  // Create output graphs for track pT resolution
  TGraphErrors* gTrackPtRes = new TGraphErrors();
  gTrackPtRes->SetName("gTrackPtResolution");
  gTrackPtRes->SetTitle("Track p_{T} resolution; p_{T}^{track} (GeV/c); #sigma_{p_{T}} / p_{T}");

  TGraphErrors* gTrackPtResMC = nullptr;
  if (gRelSigmaMC) {
    gTrackPtResMC = new TGraphErrors();
    gTrackPtResMC->SetName("gTrackPtResolutionMC");
    gTrackPtResMC->SetTitle("Track p_{T} resolution (MC); p_{T}^{track} (GeV/c); #sigma_{p_{T}} / p_{T}");
  }

  // Convert Data graph: K0s pT → track pT, σ_M/M → σ_pT/pT
  int nPoints = 0;
  for (int i = 0; i < gRelSigma->GetN(); ++i) {
    double xK0s, yMassRes;  // x: K0s pT, y: σ_M/M
    gRelSigma->GetPoint(i, xK0s, yMassRes);
    
    // Skip invalid points
    if (xK0s <= 0 || yMassRes <= 0) continue;
    
    // Convert pT: track pT ≈ 0.5 × K0s pT
    double xTrack = pTScaleFactor * xK0s;
    
    // Convert resolution: σ_pT/pT ≈ C × σ_M/M_K0s
    // Note: gRelSigma already contains σ_M/M, so we multiply by conversionFactor
    double yTrackRes = conversionFactor * yMassRes;
    
    // Error propagation for y: δ(σ_pT/pT) = C × δ(σ_M/M)
    double errYMassRes = gRelSigma->GetErrorY(i);
    double errYTrackRes = conversionFactor * errYMassRes;
    
    // Error for x: δ(pT_track) = 0.5 × δ(pT_K0s)
    double errXK0s = gRelSigma->GetErrorX(i);
    double errXTrack = pTScaleFactor * errXK0s;
    
    gTrackPtRes->SetPoint(nPoints, xTrack, yTrackRes);
    gTrackPtRes->SetPointError(nPoints, errXTrack, errYTrackRes);
    nPoints++;
  }

  // Convert MC graph if available
  int nPointsMC = 0;
  if (gRelSigmaMC) {
    for (int i = 0; i < gRelSigmaMC->GetN(); ++i) {
      double xK0s, yMassRes;
      gRelSigmaMC->GetPoint(i, xK0s, yMassRes);
      
      if (xK0s <= 0 || yMassRes <= 0) continue;
      
      double xTrack = pTScaleFactor * xK0s;
      double yTrackRes = conversionFactor * yMassRes;
      
      double errYMassRes = gRelSigmaMC->GetErrorY(i);
      double errYTrackRes = conversionFactor * errYMassRes;
      
      double errXK0s = gRelSigmaMC->GetErrorX(i);
      double errXTrack = pTScaleFactor * errXK0s;
      
      gTrackPtResMC->SetPoint(nPointsMC, xTrack, yTrackRes);
      gTrackPtResMC->SetPointError(nPointsMC, errXTrack, errYTrackRes);
      nPointsMC++;
    }
  }

  printf("[Info] Converted %d data points\n", nPoints);
  if (gTrackPtResMC) {
    printf("[Info] Converted %d MC points\n", nPointsMC);
  }

  // Create output directory
  gSystem->MakeDirectory(outputDir);

  // Draw and save plot
  TCanvas* c1 = new TCanvas("c_track_pt_res", "Track p_{T} Resolution", 800, 600);
  
  // Style data graph
  gTrackPtRes->SetMarkerStyle(21);
  gTrackPtRes->SetMarkerColor(kRed);
  gTrackPtRes->SetLineColor(kRed);
  gTrackPtRes->GetYaxis()->SetRangeUser(0.0, 0.15);
  gTrackPtRes->Draw("AP");

  // Add MC if available
  if (gTrackPtResMC) {
    gTrackPtResMC->SetMarkerStyle(20);
    gTrackPtResMC->SetMarkerColor(kBlack);
    gTrackPtResMC->SetLineColor(kBlack);
    gTrackPtResMC->Draw("P same");
  }

  // Legend
  TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->AddEntry(gTrackPtRes, "Data", "lpe");
  if (gTrackPtResMC) {
    leg->AddEntry(gTrackPtResMC, "MC", "lpe");
  }
  leg->Draw();

  // Add conversion formula annotation
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextSize(0.03);
  lat.SetTextColor(kGray+2);
  lat.DrawLatex(0.12, 0.25, Form("Conversion: p_{T}^{track} #approx %.2f #times p_{T}^{K^{0}_{S}}", pTScaleFactor));
  lat.DrawLatex(0.12, 0.20, Form("#sigma_{p_{T}}/p_{T} #approx %.3f #times #sigma_{M}/M_{K^{0}_{S}}", conversionFactor));
  lat.DrawLatex(0.12, 0.15, Form("M_{K^{0}_{S}} = %.4f GeV/c^{2}", k0sMass));

  c1->SaveAs(Form("%s/TrackPtResolution.pdf", outputDir));
  printf("[Info] Saved plot to: %s/TrackPtResolution.pdf\n", outputDir);
  delete c1;

  // Save graphs to ROOT file
  TFile* fOut = TFile::Open(Form("%s/TrackPtResolution.root", outputDir), "RECREATE");
  if (fOut && !fOut->IsZombie()) {
    gTrackPtRes->Write();
    if (gTrackPtResMC) {
      gTrackPtResMC->Write();
    }
    fOut->Close();
    printf("[Info] Saved graphs to: %s/TrackPtResolution.root\n", outputDir);
    delete fOut;
  }

  fIn->Close();
}

// Convenience function with default paths
void ConvertV0MassToTrackPtResolution() {
  ConvertV0MassToTrackPtResolutions("plots/V0MassFits/V0MassFitResults.root", 
                                    "plots/TrackPtResolution");
}

