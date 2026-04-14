// ValidateAllRRun.C
// Validates that the all-R run (614373, |eta|<0.2 for all R) is consistent with
// individual R runs (562003-562009, |eta|<0.9-R) at the overlapping |eta|<0.2 region.
//
// Category 1 (7 plots): Normalized per-event yield at |eta|<0.2 — AllR vs Individual
// Category 2 (7 plots): Per-unit-eta yield — fixed vs flexible eta from individual runs
// Bonus: Summary overlays, h_jet_pt vs h3 projection cross-check

#include "TFile.h"
#include "TH3.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"
#include "../Filipad2.h"

// R -> color mapping
Color_t GetColorForR_v(double R) {
  int idx = (int)(R * 10 + 0.5) - 1;
  static const Color_t map[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta+2, kCyan+2, kOrange+7};
  if (idx >= 0 && idx < 7) return map[idx];
  return kBlack;
}

// Jet pT bin edges
const int nPtBins = 20;
double ptBinEdges[nPtBins + 1] = {5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};

// ---------------------------------------------------------------------------
// Project TH3 (pT, eta, phi) -> 1D pT within |eta| < etaCut
// ---------------------------------------------------------------------------
TH1* ProjectH3toPt(TH3* h3, double etaCut, const char* name) {
  int binLo = h3->GetYaxis()->FindBin(-etaCut + 1e-6);
  int binHi = h3->GetYaxis()->FindBin( etaCut - 1e-6);
  TH1* hProj = h3->ProjectionX(name, binLo, binHi);
  hProj->SetDirectory(0);
  return hProj;
}

// ---------------------------------------------------------------------------
// Rebin to jet pT bins and divide by bin width -> dN/dpT
// ---------------------------------------------------------------------------
TH1* RebinToPtBins(TH1* h, const char* name) {
  TH1* hRebin = h->Rebin(nPtBins, name, ptBinEdges);
  hRebin->SetDirectory(0);
  // Divide by bin width for density
  for (int i = 1; i <= hRebin->GetNbinsX(); i++) {
    double w = hRebin->GetBinWidth(i);
    hRebin->SetBinContent(i, hRebin->GetBinContent(i) / w);
    hRebin->SetBinError(i, hRebin->GetBinError(i) / w);
  }
  return hRebin;
}

// ---------------------------------------------------------------------------
// Get number of selected events (last bin of h_collisions)
// ---------------------------------------------------------------------------
double GetNevents(TFile* f, const char* dir) {
  TH1* h = (TH1*)f->Get(Form("%s/h_collisions", dir));
  if (!h) {
    Printf("ERROR: h_collisions not found in directory '%s'!", dir);
    return 1.0;
  }
  double nev = h->GetBinContent(h->FindBin(3.5));  // after all cuts (occupancy)
  Printf("    Nevents from %s: %.0f", dir, nev);
  return nev;
}

// ---------------------------------------------------------------------------
void ValidateAllRRun() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TString baseDir = gSystem->ExpandPathName("~/cernbox/workspace/O2Physics/jets/AnalysisResults/");
  TString outDir  = "plots/ValidateAllR/";
  gSystem->mkdir(outDir, true);

  const int nR = 7;
  double Rvals[nR]       = {0.1,      0.2,      0.3,      0.4,      0.5,      0.6,      0.7};
  TString indivRunID[nR] = {"562004", "562003", "562005", "562006", "562007", "562008", "562009"};
  TString allRDir[nR]    = {"jet-spectra-charged_R01", "jet-spectra-charged_R02",
                            "jet-spectra-charged_R03", "jet-spectra-charged_R04",
                            "jet-spectra-charged_R05", "jet-spectra-charged_R06",
                            "jet-spectra-charged"};

  const double fixedEta = 0.2;
  const double deltaEtaFixed = 2.0 * fixedEta;  // = 0.4
  TString h3Name = "h3_jet_pt_jet_eta_jet_phi";

  // --- Open files ---
  TFile* fAllR = TFile::Open(baseDir + "614373_AnalysisResults.root", "READ");
  if (!fAllR || fAllR->IsZombie()) { Printf("ERROR: Cannot open 614373"); return; }

  TFile* fIndiv[nR];
  for (int i = 0; i < nR; i++) {
    fIndiv[i] = TFile::Open(baseDir + indivRunID[i] + "_AnalysisResults.root", "READ");
    if (!fIndiv[i] || fIndiv[i]->IsZombie()) {
      Printf("ERROR: Cannot open %s", indivRunID[i].Data()); return;
    }
  }

  // --- Sanity checks (Findings 2,3,4) ---
  Printf("=== Sanity checks ===");

  // Check h3 axis binning (Finding 3)
  TH3* h3check = (TH3*)fAllR->Get(Form("%s/%s", allRDir[0].Data(), h3Name.Data()));
  if (h3check) {
    Printf("  h3 axes: X(pT) %d bins [%.0f,%.0f], Y(eta) %d bins [%.1f,%.1f], Z(phi) %d bins [%.1f,%.1f]",
           h3check->GetXaxis()->GetNbins(), h3check->GetXaxis()->GetXmin(), h3check->GetXaxis()->GetXmax(),
           h3check->GetYaxis()->GetNbins(), h3check->GetYaxis()->GetXmin(), h3check->GetYaxis()->GetXmax(),
           h3check->GetZaxis()->GetNbins(), h3check->GetZaxis()->GetXmin(), h3check->GetZaxis()->GetXmax());
    if (TMath::Abs(h3check->GetYaxis()->GetXmin() + 1.0) > 1e-6 ||
        TMath::Abs(h3check->GetYaxis()->GetXmax() - 1.0) > 1e-6) {
      Printf("  WARNING: h3 eta axis range is NOT [-1,1]! FindBin logic may be wrong!");
    }
    // Note: FindBin + epsilon approach is robust to any nBins; no warning needed for nBins != 200
  }

  // Verify individual run eta ranges match expected |eta|<0.9-R (Finding 2)
  for (int i = 0; i < nR; i++) {
    TH3* h3 = (TH3*)fIndiv[i]->Get(Form("jet-spectra-charged/%s", h3Name.Data()));
    if (!h3) continue;
    TH1* hEta = h3->ProjectionY(Form("hEtaChk_%d", i));
    hEta->SetDirectory(0);
    double etaMaxPopulated = 0;
    for (int b = hEta->GetNbinsX(); b >= 1; b--) {
      if (hEta->GetBinContent(b) > 0) {
        etaMaxPopulated = hEta->GetXaxis()->GetBinUpEdge(b);
        break;
      }
    }
    double expectedEta = 0.9 - Rvals[i];
    double etaBinWidth = h3->GetYaxis()->GetBinWidth(1);
    Printf("  Indiv R=%.1f (%s): populated eta up to %.2f (expected %.1f + binwidth %.2f = %.2f)",
           Rvals[i], indivRunID[i].Data(), etaMaxPopulated,
           expectedEta, etaBinWidth, expectedEta + etaBinWidth);
    // Jets at eta = expectedEta fall into bin with upper edge = expectedEta + binWidth
    if (TMath::Abs(etaMaxPopulated - (expectedEta + etaBinWidth)) > etaBinWidth + 0.01) {
      Printf("  WARNING: R=%.1f eta range MISMATCH! Check R->RunID mapping!", Rvals[i]);
    }
    delete hEta;
  }

  // Verify AllR Nevents is consistent across all 7 directories (Finding 4)
  Printf("\n  AllR Nevents per R directory:");
  double nevAllR_check[nR];
  for (int i = 0; i < nR; i++) {
    TH1* hc = (TH1*)fAllR->Get(Form("%s/h_collisions", allRDir[i].Data()));
    if (hc) {
      nevAllR_check[i] = hc->GetBinContent(hc->FindBin(3.5));
      Printf("    %s: %.0f", allRDir[i].Data(), nevAllR_check[i]);
    } else {
      nevAllR_check[i] = -1;
      Printf("    %s: h_collisions NOT FOUND!", allRDir[i].Data());
    }
  }
  // Check consistency
  for (int i = 1; i < nR; i++) {
    if (nevAllR_check[i] > 0 && nevAllR_check[0] > 0 &&
        TMath::Abs(nevAllR_check[i] - nevAllR_check[0]) > 0.5) {
      Printf("  WARNING: Nevents differs between %s (%.0f) and %s (%.0f)!",
             allRDir[0].Data(), nevAllR_check[0], allRDir[i].Data(), nevAllR_check[i]);
    }
  }

  // --- Get event counts ---
  Printf("\n=== Event counts ===");
  double nevAllR = GetNevents(fAllR, "jet-spectra-charged");
  double nevIndiv[nR];
  for (int i = 0; i < nR; i++)
    nevIndiv[i] = GetNevents(fIndiv[i], "jet-spectra-charged");

  std::vector<Filipad2*> pads;

  // ==========================================================================
  // Category 1: Per-event yield at |eta|<0.2 — AllR vs Individual
  //   Normalize: (1/Nevt) * dN/dpT  [rebinned]
  // ==========================================================================
  Printf("\n=== Category 1: Per-event yield comparison ===");
  for (int i = 0; i < nR; i++) {
    double R = Rvals[i];
    Color_t col = GetColorForR_v(R);

    TH3* h3AllR  = (TH3*)fAllR->Get(Form("%s/%s", allRDir[i].Data(), h3Name.Data()));
    TH3* h3Indiv = (TH3*)fIndiv[i]->Get(Form("jet-spectra-charged/%s", h3Name.Data()));
    if (!h3AllR || !h3Indiv) { Printf("  SKIP R=%.1f: h3 missing", R); continue; }

    // Project to |eta|<0.2, rebin, normalize
    TH1* hAllR_proj  = ProjectH3toPt(h3AllR,  fixedEta, Form("hA1_%d", i));
    TH1* hIndiv_proj = ProjectH3toPt(h3Indiv, fixedEta, Form("hI1_%d", i));

    TH1* hAllR  = RebinToPtBins(hAllR_proj,  Form("hAllR_c1_R%02d",  (int)(R*10+0.5)));
    TH1* hIndiv = RebinToPtBins(hIndiv_proj, Form("hIndiv_c1_R%02d", (int)(R*10+0.5)));

    hAllR->Scale(1.0 / nevAllR);       // per-event yield
    hIndiv->Scale(1.0 / nevIndiv[i]);

    Printf("  R=%.1f: AllR yield(20-25)=%.6e, Indiv yield(20-25)=%.6e, ratio=%.4f",
           R, hAllR->GetBinContent(hAllR->FindBin(22.)),
           hIndiv->GetBinContent(hIndiv->FindBin(22.)),
           hIndiv->GetBinContent(hIndiv->FindBin(22.)) > 0 ?
           hAllR->GetBinContent(hAllR->FindBin(22.)) / hIndiv->GetBinContent(hIndiv->FindBin(22.)) : 0.);

    // --- Draw ---
    Filipad2* pad = new Filipad2(Form("Cat1_R%02d", (int)(R*10+0.5)), 1, 2.5, 0.35, 100+i*50, 100, 0.7);
    pads.push_back(pad);
    pad->Draw();

    TPad* p1 = pad->GetPad(1);
    p1->SetLogy(); p1->cd();

    double ymax = TMath::Max(hIndiv->GetMaximum(), hAllR->GetMaximum());
    TH1* hFrame = (TH1*)hIndiv->Clone(Form("hFr1_%d", i)); hFrame->Reset();
    hFrame->SetDirectory(0);
    hFrame->SetMaximum(ymax * 20); hFrame->SetMinimum(ymax * 1e-10);
    pad->Hset(hFrame, "", "(1/N_{evt}) dN/d#it{p}_{T} (GeV/#it{c})^{-1}",
              2.5, 2.0, 20, 20, 0.01, 0.001, 16, 16);
    hFrame->Draw();

    hIndiv->SetMarkerStyle(20); hIndiv->SetMarkerSize(1.0);
    hIndiv->SetMarkerColor(kBlack); hIndiv->SetLineColor(kBlack);
    hIndiv->Draw("same ep");

    hAllR->SetMarkerStyle(24); hAllR->SetMarkerSize(1.0);
    hAllR->SetMarkerColor(kRed); hAllR->SetLineColor(kRed);
    hAllR->Draw("same ep");

    TLegend* leg = new TLegend(0.35, 0.60, 0.93, 0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->SetTextFont(43); leg->SetTextSize(18);
    leg->SetHeader(Form("R = %.1f, |#eta_{jet}| < %.1f, per event", R, fixedEta));
    leg->AddEntry(hIndiv, Form("Indiv %s (N=%.2gM)", indivRunID[i].Data(), nevIndiv[i]/1e6), "lpe");
    leg->AddEntry(hAllR,  Form("AllR 614373 (N=%.2gM)", nevAllR/1e6), "lpe");
    leg->Draw();

    // Ratio
    TPad* p2 = pad->GetPad(2); p2->cd();
    TH1* hRat = (TH1*)hAllR->Clone(Form("hR1_%d", i));
    hRat->SetDirectory(0);
    hRat->Divide(hIndiv);
    hRat->SetMinimum(0.85); hRat->SetMaximum(1.15);
    hRat->SetMarkerStyle(20); hRat->SetMarkerSize(0.8);
    hRat->SetMarkerColor(col); hRat->SetLineColor(col);
    pad->Hset(hRat, "#it{p}_{T} (GeV/#it{c})", "AllR / Indiv",
              3.5, 2.0, 20, 20, 0.01, 0.001, 16, 16);
    hRat->Draw("ep");
    TLine* l = new TLine(ptBinEdges[0], 1.0, ptBinEdges[nPtBins], 1.0);
    l->SetLineStyle(2); l->Draw();

    pad->C->SaveAs(outDir + Form("Cat1_PerEvtYield_R%02d.pdf", (int)(R*10+0.5)));
    Printf("  -> Saved Cat1 R=%.1f", R);
  }

  // ==========================================================================
  // Category 2: Per-unit-eta yield — fixed vs flexible eta (individual runs)
  //   Normalize by delta_eta: spectrum / (2 * etaCut) → expected ratio = 1
  // ==========================================================================
  Printf("\n=== Category 2: Per-unit-eta yield ===");
  for (int i = 0; i < nR; i++) {
    double R = Rvals[i];
    double etaStd = 0.9 - R;
    double deltaEtaStd = 2.0 * etaStd;
    Color_t col = GetColorForR_v(R);

    TH3* h3Indiv = (TH3*)fIndiv[i]->Get(Form("jet-spectra-charged/%s", h3Name.Data()));
    if (!h3Indiv) { Printf("  SKIP R=%.1f", R); continue; }

    TH1* hFixed_proj = ProjectH3toPt(h3Indiv, fixedEta, Form("hF2_%d", i));
    TH1* hStd_proj   = ProjectH3toPt(h3Indiv, etaStd,   Form("hS2_%d", i));

    TH1* hFixed = RebinToPtBins(hFixed_proj, Form("hFixed_c2_R%02d", (int)(R*10+0.5)));
    TH1* hStd   = RebinToPtBins(hStd_proj,   Form("hStd_c2_R%02d",   (int)(R*10+0.5)));

    // Normalize by delta_eta (and Nevt for proper yield)
    hFixed->Scale(1.0 / (nevIndiv[i] * deltaEtaFixed));
    hStd->Scale(1.0 / (nevIndiv[i] * deltaEtaStd));

    Printf("  R=%.1f: dEta_fixed=%.1f, dEta_std=%.1f, yield ratio(20-25)=%.4f",
           R, deltaEtaFixed, deltaEtaStd,
           hStd->GetBinContent(hStd->FindBin(22.)) > 0 ?
           hFixed->GetBinContent(hFixed->FindBin(22.)) / hStd->GetBinContent(hStd->FindBin(22.)) : 0.);

    // --- Draw ---
    Filipad2* pad = new Filipad2(Form("Cat2_R%02d", (int)(R*10+0.5)), 1, 2.5, 0.35, 100+i*50, 200, 0.7);
    pads.push_back(pad);
    pad->Draw();

    TPad* p1 = pad->GetPad(1);
    p1->SetLogy(); p1->cd();

    double ymax = TMath::Max(hStd->GetMaximum(), hFixed->GetMaximum());
    TH1* hFrame = (TH1*)hStd->Clone(Form("hFr2_%d", i)); hFrame->Reset();
    hFrame->SetDirectory(0);
    hFrame->SetMaximum(ymax * 20); hFrame->SetMinimum(ymax * 1e-10);
    pad->Hset(hFrame, "", "(1/N_{evt}) d^{2}N/(d#it{p}_{T} d#eta) (GeV/#it{c})^{-1}",
              2.5, 2.0, 20, 20, 0.01, 0.001, 16, 16);
    hFrame->Draw();

    hStd->SetMarkerStyle(20); hStd->SetMarkerSize(1.0);
    hStd->SetMarkerColor(kBlack); hStd->SetLineColor(kBlack);
    hStd->Draw("same ep");

    hFixed->SetMarkerStyle(24); hFixed->SetMarkerSize(1.0);
    hFixed->SetMarkerColor(kRed); hFixed->SetLineColor(kRed);
    hFixed->Draw("same ep");

    TLegend* leg = new TLegend(0.35, 0.60, 0.93, 0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->SetTextFont(43); leg->SetTextSize(18);
    leg->SetHeader(Form("R = %.1f, run %s, per event per #eta", R, indivRunID[i].Data()));
    leg->AddEntry(hStd,   Form("|#eta| < %.1f / %.1f", etaStd, deltaEtaStd), "lpe");
    leg->AddEntry(hFixed, Form("|#eta| < %.1f / %.1f", fixedEta, deltaEtaFixed), "lpe");
    leg->Draw();

    // Ratio: binomial errors (Fixed is a subset of Std)
    TPad* p2 = pad->GetPad(2); p2->cd();
    // Use raw rebinned counts (no bin-width / normalization) for proper "B" errors
    TH1* hFixed_cnt = hFixed_proj->Rebin(nPtBins, Form("hFcnt_%d", i), ptBinEdges);
    hFixed_cnt->SetDirectory(0);
    TH1* hStd_cnt = hStd_proj->Rebin(nPtBins, Form("hScnt_%d", i), ptBinEdges);
    hStd_cnt->SetDirectory(0);
    TH1* hRat = (TH1*)hFixed_cnt->Clone(Form("hR2_%d", i));
    hRat->SetDirectory(0);
    hRat->Divide(hFixed_cnt, hStd_cnt, 1, 1, "B");  // binomial errors for subset
    hRat->Scale(deltaEtaStd / deltaEtaFixed);         // scale so expected ratio = 1.0
    hRat->SetMinimum(0.85); hRat->SetMaximum(1.15);
    hRat->SetMarkerStyle(20); hRat->SetMarkerSize(0.8);
    hRat->SetMarkerColor(col); hRat->SetLineColor(col);
    pad->Hset(hRat, "#it{p}_{T} (GeV/#it{c})", "Fixed / Std",
              3.5, 2.0, 20, 20, 0.01, 0.001, 16, 16);
    hRat->Draw("ep");
    TLine* l = new TLine(ptBinEdges[0], 1.0, ptBinEdges[nPtBins], 1.0);
    l->SetLineStyle(2); l->Draw();

    pad->C->SaveAs(outDir + Form("Cat2_PerUnitEta_R%02d.pdf", (int)(R*10+0.5)));
    Printf("  -> Saved Cat2 R=%.1f", R);
  }

  // ==========================================================================
  // Bonus 1: Summary overlay — all R from AllR run (per-event yield)
  // ==========================================================================
  Printf("\n=== Bonus: Summary overlays ===");
  {
    TCanvas* cSum = new TCanvas("cSumAllR", "AllR per-event yield", 800, 700);
    cSum->SetLogy();
    cSum->SetLeftMargin(0.14); cSum->SetBottomMargin(0.12); cSum->SetRightMargin(0.03);
    TLegend* leg = new TLegend(0.50, 0.45, 0.90, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->SetTextFont(43); leg->SetTextSize(20);
    leg->SetHeader("614373, |#eta| < 0.2, per event");

    bool first = true;
    for (int i = 0; i < nR; i++) {
      TH3* h3 = (TH3*)fAllR->Get(Form("%s/%s", allRDir[i].Data(), h3Name.Data()));
      if (!h3) continue;
      TH1* hp = ProjectH3toPt(h3, fixedEta, Form("hSA_%d", i));
      TH1* h = RebinToPtBins(hp, Form("hSumA_R%02d", (int)(Rvals[i]*10+0.5)));
      h->Scale(1.0 / nevAllR);
      h->SetMarkerStyle(20); h->SetMarkerSize(1.0);
      h->SetMarkerColor(GetColorForR_v(Rvals[i]));
      h->SetLineColor(GetColorForR_v(Rvals[i]));
      if (first) {
        h->SetMinimum(1e-12); h->SetMaximum(h->GetMaximum() * 100);
        h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        h->GetYaxis()->SetTitle("(1/N_{evt}) dN/d#it{p}_{T}");
        h->GetXaxis()->SetTitleFont(43); h->GetXaxis()->SetTitleSize(24);
        h->GetYaxis()->SetTitleFont(43); h->GetYaxis()->SetTitleSize(24);
        h->GetXaxis()->SetLabelFont(43); h->GetXaxis()->SetLabelSize(20);
        h->GetYaxis()->SetLabelFont(43); h->GetYaxis()->SetLabelSize(20);
        h->GetYaxis()->SetTitleOffset(1.6);
        h->Draw("ep"); first = false;
      } else { h->Draw("same ep"); }
      leg->AddEntry(h, Form("R = %.1f", Rvals[i]), "lpe");
    }
    leg->Draw();
    cSum->SaveAs(outDir + "Summary_AllR_PerEvtYield.pdf");
    Printf("  -> Saved Summary AllR");
  }

  // ==========================================================================
  // Bonus 2: Summary overlay — individual runs (per-event yield at |eta|<0.2)
  // ==========================================================================
  {
    TCanvas* cSum2 = new TCanvas("cSumInd", "Individual per-event yield", 800, 700);
    cSum2->SetLogy();
    cSum2->SetLeftMargin(0.14); cSum2->SetBottomMargin(0.12); cSum2->SetRightMargin(0.03);
    TLegend* leg = new TLegend(0.50, 0.40, 0.90, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->SetTextFont(43); leg->SetTextSize(20);
    leg->SetHeader("Individual runs, |#eta| < 0.2, per event");

    bool first = true;
    for (int i = 0; i < nR; i++) {
      TH3* h3 = (TH3*)fIndiv[i]->Get(Form("jet-spectra-charged/%s", h3Name.Data()));
      if (!h3) continue;
      TH1* hp = ProjectH3toPt(h3, fixedEta, Form("hSI_%d", i));
      TH1* h = RebinToPtBins(hp, Form("hSumI_R%02d", (int)(Rvals[i]*10+0.5)));
      h->Scale(1.0 / nevIndiv[i]);
      h->SetMarkerStyle(20); h->SetMarkerSize(1.0);
      h->SetMarkerColor(GetColorForR_v(Rvals[i]));
      h->SetLineColor(GetColorForR_v(Rvals[i]));
      if (first) {
        h->SetMinimum(1e-12); h->SetMaximum(h->GetMaximum() * 100);
        h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        h->GetYaxis()->SetTitle("(1/N_{evt}) dN/d#it{p}_{T}");
        h->GetXaxis()->SetTitleFont(43); h->GetXaxis()->SetTitleSize(24);
        h->GetYaxis()->SetTitleFont(43); h->GetYaxis()->SetTitleSize(24);
        h->GetXaxis()->SetLabelFont(43); h->GetXaxis()->SetLabelSize(20);
        h->GetYaxis()->SetLabelFont(43); h->GetYaxis()->SetLabelSize(20);
        h->GetYaxis()->SetTitleOffset(1.6);
        h->Draw("ep"); first = false;
      } else { h->Draw("same ep"); }
      leg->AddEntry(h, Form("R = %.1f (%s)", Rvals[i], indivRunID[i].Data()), "lpe");
    }
    leg->Draw();
    cSum2->SaveAs(outDir + "Summary_Individual_PerEvtYield.pdf");
    Printf("  -> Saved Summary Individual");
  }

  // ==========================================================================
  // Bonus 3: h_jet_pt vs h3 projection cross-check (individual runs)
  // Note: For AllR this is trivially 1:1 (both fill same |eta|<0.2).
  //       For individual runs it's also expected to match but is a better
  //       data integrity check since h3 has broader eta content.
  // ==========================================================================
  {
    TCanvas* cX = new TCanvas("cXchk", "h1D vs h3 proj (Indiv)", 1600, 800);
    cX->Divide(4, 2);
    for (int i = 0; i < nR; i++) {
      cX->cd(i + 1); gPad->SetLogy();
      gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.12);
      TH1* h1D = (TH1*)fIndiv[i]->Get("jet-spectra-charged/h_jet_pt");
      TH3* h3  = (TH3*)fIndiv[i]->Get(Form("jet-spectra-charged/%s", h3Name.Data()));
      if (!h1D || !h3) continue;
      h1D = (TH1*)h1D->Clone(Form("h1Dc_%d", i)); h1D->SetDirectory(0);
      TH1* hP = h3->ProjectionX(Form("hPF_%d", i)); hP->SetDirectory(0);
      h1D->SetMarkerStyle(20); h1D->SetMarkerSize(0.6);
      h1D->SetMarkerColor(kBlack); h1D->SetLineColor(kBlack);
      h1D->GetXaxis()->SetRangeUser(0, 150);
      h1D->GetXaxis()->SetTitleFont(43); h1D->GetXaxis()->SetTitleSize(18);
      h1D->GetYaxis()->SetTitleFont(43); h1D->GetYaxis()->SetTitleSize(18);
      h1D->GetXaxis()->SetLabelFont(43); h1D->GetXaxis()->SetLabelSize(14);
      h1D->GetYaxis()->SetLabelFont(43); h1D->GetYaxis()->SetLabelSize(14);
      h1D->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      h1D->GetYaxis()->SetTitle("Counts");
      gStyle->SetOptTitle(1);
      h1D->SetTitle(Form("R=%.1f (%s)", Rvals[i], indivRunID[i].Data()));
      h1D->Draw("ep");
      hP->SetMarkerStyle(24); hP->SetMarkerSize(0.6);
      hP->SetMarkerColor(kRed); hP->SetLineColor(kRed);
      hP->Draw("same ep");
      TLegend* lg = new TLegend(0.30, 0.72, 0.88, 0.88);
      lg->SetBorderSize(0); lg->SetFillStyle(0);
      lg->SetTextFont(43); lg->SetTextSize(14);
      lg->AddEntry(h1D, "h_jet_pt", "lpe");
      lg->AddEntry(hP, "h3 ProjX (full #eta)", "lpe");
      lg->Draw();
    }
    gStyle->SetOptTitle(0);
    cX->SaveAs(outDir + "Bonus_h1D_vs_h3Projection_Indiv.pdf");
    Printf("  -> Saved h1D vs h3 cross-check (Individual runs)");
  }

  // --- Cleanup ---
  fAllR->Close(); delete fAllR;
  for (int i = 0; i < nR; i++) { fIndiv[i]->Close(); delete fIndiv[i]; }

  Printf("\n=== Done. 17 plots saved to %s ===", outDir.Data());
}
