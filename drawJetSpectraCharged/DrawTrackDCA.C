// Draw Track DCA Distributions from trackEfficiency.cxx output
// Usage: root -l DrawTrackDCA.C

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include <iostream>
#include <vector>

// ============================================================
//  Input files — PUT YOUR ROOT FILE PATHS HERE
// ============================================================
const char* kBaseDir = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/";
// Data
const char* kData2022      = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/615987_AnalysisResults.root"; // 2022 data
const char* kData2023      = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/616559_AnalysisResults.root"; // 2023 data
// MC — processTrackSelectionHistograms was NOT enabled in these trains; re-run to get DCA plots
// const char* kJJMC2022      = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/617213_AnalysisResults.root"; // 2022 JJ MC — no h_trackselplot_*
// const char* kJJMC2022_10p  = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/617825_AnalysisResults.root"; // 2022 JJ MC 10% — no h_trackselplot_*
// const char* kMBMC2022      = "/Users/js/cernbox/workspace/O2Physics/jets/AnalysisResults/613497_AnalysisResults.root"; // 2022 MB MC — no h_trackselplot_*

// ============================================================
//  General configuration
// ============================================================
const char* kOutputDir = "plots/trackDCA";
const char* kTaskDir   = "track-efficiency"; // TDirectory inside AnalysisResults

// pT slices for DCA projections from 2D histograms
const double kPtSlicesLo[] = {0.15, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0};
const double kPtSlicesHi[] = {0.5,  1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0};
const int    kNPtSlices = 8;

const int kSliceColors[] = {kBlack, kRed+1, kBlue+1, kGreen+2,
                            kMagenta+1, kOrange+1, kCyan+2, kViolet+1};

// DCA_z cut values for secondary contamination study
const double kDCAzCuts[]    = {0.1, 0.2, 0.5};
const int    kDCAzCutColors[] = {kRed+1, kBlue+1, kGreen+2};
const int    kDCAzCutStyles[] = {2, 1, 7}; // dashed, solid, long-dashed
const int    kNDCAzCuts = 3;

// ============================================================
//  Input source definition
// ============================================================
struct Source {
    const char* path;
    const char* label;
    int   color;
    int   marker;
    int   lineStyle;
};

// --- Data sources ---
const Source kDataSources[] = {
    {kData2022, "Data 2022", kBlack,    20, 1},
    {kData2023, "Data 2023", kRed+1,    21, 1},
};
const int kNData = sizeof(kDataSources) / sizeof(Source);

// --- MC sources — uncomment when re-run with processTrackSelectionHistograms=true ---
// const Source kMCSources[] = {
//     {kJJMC2022,     "JJ MC 2022",      kBlue+1,    24, 1},
//     {kJJMC2022_10p, "JJ MC 2022 (10%)", kCyan+2,   25, 2},
//     {kMBMC2022,     "MB MC 2022",      kGreen+2,   26, 1},
// };
// const int kNMC = sizeof(kMCSources) / sizeof(Source);

// --- All sources combined — add MC back when available ---
const Source kAllSources[] = {
    {kData2022, "Data 2022", kBlack,    20, 1},
    {kData2023, "Data 2023", kRed+1,    21, 1},
};
const int kNAll = sizeof(kAllSources) / sizeof(Source);

// ============================================================
//  Helpers
// ============================================================
template <typename T>
T* GetHist(TFile* f, const char* name)
{
    TString path = Form("%s/%s", kTaskDir, name);
    T* h = dynamic_cast<T*>(f->Get(path));
    if (!h) std::cerr << "WARNING: " << path << " not found in " << f->GetName() << std::endl;
    return h;
}

TFile* OpenFile(const char* path)
{
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "WARNING: Cannot open " << path << std::endl;
        return nullptr;
    }
    return f;
}

void SetStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetTitleSize(0.045, "XYZ");
    gStyle->SetTitleOffset(1.1, "X");
    gStyle->SetTitleOffset(1.3, "Y");
}

void ApplyStyle(TH1* h, const Source& s)
{
    h->SetLineColor(s.color);
    h->SetMarkerColor(s.color);
    h->SetMarkerStyle(s.marker);
    h->SetMarkerSize(0.9);
    h->SetLineWidth(2);
    h->SetLineStyle(s.lineStyle);
}

void ApplyStyle(TGraphErrors* g, const Source& s)
{
    g->SetLineColor(s.color);
    g->SetMarkerColor(s.color);
    g->SetMarkerStyle(s.marker);
    g->SetMarkerSize(0.9);
    g->SetLineWidth(2);
    g->SetLineStyle(s.lineStyle);
}

// ============================================================
//  Gaussian fit to core: iterative 2-sigma fit
//  Returns {sigma, sigmaErr}. Fits +-2*RMS, then refines to +-2*sigma.
// ============================================================
std::pair<double,double> FitCoreSigma(TH1* h, double fitRangeInit = 0)
{
    double rms = h->GetStdDev();
    double mean = h->GetMean();
    double range = (fitRangeInit > 0) ? fitRangeInit : 2.0 * rms;

    TF1 gaus("gaus_tmp", "gaus", mean - range, mean + range);
    gaus.SetParameters(h->GetMaximum(), mean, rms * 0.5);
    h->Fit(&gaus, "QNR0");

    // Refit with refined range
    double sig1 = gaus.GetParameter(2);
    double mu1  = gaus.GetParameter(1);
    gaus.SetRange(mu1 - 2.0 * sig1, mu1 + 2.0 * sig1);
    h->Fit(&gaus, "QNR0");

    return {std::abs(gaus.GetParameter(2)), gaus.GetParError(2)};
}

// ============================================================
//  1) Draw 1D DCA_xy — all sources overlaid
// ============================================================
void DrawDCAxy1D(const Source* sources, int nSrc)
{
    auto* c = new TCanvas("cDCAxy1D", "", 700, 600);
    c->SetLogy();
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.12);

    auto* leg = new TLegend(0.60, 0.70, 0.92, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);

    TLatex sigTex; sigTex.SetNDC(); sigTex.SetTextSize(0.030);

    bool first = true;
    int nDrawn = 0;
    for (int i = 0; i < nSrc; ++i) {
        TFile* f = OpenFile(sources[i].path);
        if (!f) continue;
        auto* h = GetHist<TH1>(f, "h_trackselplot_dcaxy");
        if (!h) { delete f; continue; }
        h->SetDirectory(0);
        f->Close(); delete f;

        auto [sig, sigErr] = FitCoreSigma(h);

        if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
        ApplyStyle(h, sources[i]);
        h->GetXaxis()->SetTitle("DCA_{xy} (cm)");
        h->GetYaxis()->SetTitle("Self-normalized");
        h->GetXaxis()->SetRangeUser(-0.5, 0.5);

        h->Draw(first ? "hist" : "hist same");
        first = false;
        leg->AddEntry(h, Form("%s, #sigma = %.4f cm", sources[i].label, sig), "l");
        nDrawn++;
    }
    if (first) { delete c; return; }

    leg->Draw();
    TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035);
    tex.DrawLatex(0.16, 0.92, "DCA_{xy},  pp #sqrt{#it{s}} = 13.6 TeV");
    tex.SetTextSize(0.025);
    tex.DrawLatex(0.16, 0.87, "#sigma from Gaussian fit to core (#pm2#sigma)");

    c->SaveAs(Form("%s/DCAxy_1D.pdf", kOutputDir));
    delete c;
}

// ============================================================
//  2) Draw 1D DCA_z — all sources overlaid
// ============================================================
void DrawDCAz1D(const Source* sources, int nSrc)
{
    auto* c = new TCanvas("cDCAz1D", "", 700, 600);
    c->SetLogy();
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.12);

    auto* leg = new TLegend(0.60, 0.70, 0.92, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);

    bool first = true;
    for (int i = 0; i < nSrc; ++i) {
        TFile* f = OpenFile(sources[i].path);
        if (!f) continue;
        auto* h = GetHist<TH1>(f, "h_trackselplot_dcaz");
        if (!h) { delete f; continue; }
        h->SetDirectory(0);
        f->Close(); delete f;

        auto [sig, sigErr] = FitCoreSigma(h);

        if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
        ApplyStyle(h, sources[i]);
        h->GetXaxis()->SetTitle("DCA_{z} (cm)");
        h->GetYaxis()->SetTitle("Self-normalized");
        h->GetXaxis()->SetRangeUser(-2.0, 2.0);

        h->Draw(first ? "hist" : "hist same");
        first = false;
        leg->AddEntry(h, Form("%s, #sigma = %.4f cm", sources[i].label, sig), "l");
    }
    if (first) { delete c; return; }

    leg->Draw();
    TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035);
    tex.DrawLatex(0.16, 0.92, "DCA_{z},  pp #sqrt{#it{s}} = 13.6 TeV");
    tex.SetTextSize(0.025);
    tex.DrawLatex(0.16, 0.87, "#sigma from Gaussian fit to core (#pm2#sigma)");

    c->SaveAs(Form("%s/DCAz_1D.pdf", kOutputDir));
    delete c;
}

// ============================================================
//  3) Draw 2D pT vs DCA_xy (one canvas per source)
// ============================================================
void DrawDCAxy2D(const Source* sources, int nSrc)
{
    for (int i = 0; i < nSrc; ++i) {
        TFile* f = OpenFile(sources[i].path);
        if (!f) continue;
        auto* h2 = GetHist<TH2>(f, "h2_trackselplot_pt_dcaxy");
        if (!h2) { delete f; continue; }

        auto* c = new TCanvas(Form("cDCAxy2D_%d", i), "", 800, 600);
        c->SetLogz();
        c->SetLeftMargin(0.12);
        c->SetBottomMargin(0.12);
        c->SetRightMargin(0.14);

        h2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        h2->GetYaxis()->SetTitle("DCA_{xy} (cm)");
        h2->GetXaxis()->SetRangeUser(0, 50);
        h2->GetYaxis()->SetRangeUser(-0.5, 0.5);
        h2->Draw("colz");

        TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035);
        tex.DrawLatex(0.15, 0.92, Form("%s", sources[i].label));

        c->SaveAs(Form("%s/DCAxy_vs_pT_2D_%s.pdf", kOutputDir, sources[i].label));
        delete c;
        f->Close(); delete f;
    }
}

// ============================================================
//  4) Draw 2D pT vs DCA_z (one canvas per source)
// ============================================================
void DrawDCAz2D(const Source* sources, int nSrc)
{
    for (int i = 0; i < nSrc; ++i) {
        TFile* f = OpenFile(sources[i].path);
        if (!f) continue;
        auto* h2 = GetHist<TH2>(f, "h2_trackselplot_pt_dcaz");
        if (!h2) { delete f; continue; }

        auto* c = new TCanvas(Form("cDCAz2D_%d", i), "", 800, 600);
        c->SetLogz();
        c->SetLeftMargin(0.12);
        c->SetBottomMargin(0.12);
        c->SetRightMargin(0.14);

        h2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        h2->GetYaxis()->SetTitle("DCA_{z} (cm)");
        h2->GetXaxis()->SetRangeUser(0, 50);
        h2->GetYaxis()->SetRangeUser(-2.0, 2.0);
        h2->Draw("colz");

        TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035);
        tex.DrawLatex(0.15, 0.92, Form("%s", sources[i].label));

        c->SaveAs(Form("%s/DCAz_vs_pT_2D_%s.pdf", kOutputDir, sources[i].label));
        delete c;
        f->Close(); delete f;
    }
}

// ============================================================
//  5) Draw DCA_xy projections in pT slices (one source)
// ============================================================
void DrawDCAxyPtSlices(const Source& src)
{
    TFile* f = OpenFile(src.path);
    if (!f) return;
    auto* h2 = GetHist<TH2>(f, "h2_trackselplot_pt_dcaxy");
    if (!h2) { delete f; return; }

    auto* c = new TCanvas(Form("cDCAxySlices_%s", src.label), "", 800, 600);
    c->SetLogy();
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.12);

    auto* leg = new TLegend(0.62, 0.50, 0.92, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.028);
    leg->SetHeader(src.label);

    bool first = true;
    for (int i = 0; i < kNPtSlices; ++i) {
        int binLo = h2->GetXaxis()->FindBin(kPtSlicesLo[i] + 1e-6);
        int binHi = h2->GetXaxis()->FindBin(kPtSlicesHi[i] - 1e-6);

        TH1* hp = h2->ProjectionY(Form("hDCAxy_%s_%d", src.label, i), binLo, binHi);
        if (hp->GetEntries() < 1) { delete hp; continue; }
        hp->SetDirectory(0);
        if (hp->Integral() > 0) hp->Scale(1.0 / hp->Integral());

        hp->SetLineColor(kSliceColors[i % 8]);
        hp->SetLineWidth(2);
        hp->GetXaxis()->SetTitle("DCA_{xy} (cm)");
        hp->GetYaxis()->SetTitle("Normalized");
        hp->GetXaxis()->SetRangeUser(-0.3, 0.3);

        hp->Draw(first ? "hist" : "hist same");
        first = false;
        leg->AddEntry(hp, Form("%.1f < #it{p}_{T} < %.0f", kPtSlicesLo[i], kPtSlicesHi[i]), "l");
    }
    if (!first) {
        leg->Draw();
        c->SaveAs(Form("%s/DCAxy_pTslices_%s.pdf", kOutputDir, src.label));
    }
    delete c;
    f->Close(); delete f;
}

// ============================================================
//  6) Draw DCA_z projections in pT slices (one source)
// ============================================================
void DrawDCAzPtSlices(const Source& src)
{
    TFile* f = OpenFile(src.path);
    if (!f) return;
    auto* h2 = GetHist<TH2>(f, "h2_trackselplot_pt_dcaz");
    if (!h2) { delete f; return; }

    auto* c = new TCanvas(Form("cDCAzSlices_%s", src.label), "", 800, 600);
    c->SetLogy();
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.12);

    auto* leg = new TLegend(0.62, 0.50, 0.92, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.028);
    leg->SetHeader(src.label);

    bool first = true;
    for (int i = 0; i < kNPtSlices; ++i) {
        int binLo = h2->GetXaxis()->FindBin(kPtSlicesLo[i] + 1e-6);
        int binHi = h2->GetXaxis()->FindBin(kPtSlicesHi[i] - 1e-6);

        TH1* hp = h2->ProjectionY(Form("hDCAz_%s_%d", src.label, i), binLo, binHi);
        if (hp->GetEntries() < 1) { delete hp; continue; }
        hp->SetDirectory(0);
        if (hp->Integral() > 0) hp->Scale(1.0 / hp->Integral());

        hp->SetLineColor(kSliceColors[i % 8]);
        hp->SetLineWidth(2);
        hp->GetXaxis()->SetTitle("DCA_{z} (cm)");
        hp->GetYaxis()->SetTitle("Normalized");
        hp->GetXaxis()->SetRangeUser(-1.0, 1.0);

        hp->Draw(first ? "hist" : "hist same");
        first = false;
        leg->AddEntry(hp, Form("%.1f < #it{p}_{T} < %.0f", kPtSlicesLo[i], kPtSlicesHi[i]), "l");
    }
    if (!first) {
        leg->Draw();
        c->SaveAs(Form("%s/DCAz_pTslices_%s.pdf", kOutputDir, src.label));
    }
    delete c;
    f->Close(); delete f;
}

// ============================================================
//  7) Draw DCA_xy sigma vs pT — all sources overlaid
// ============================================================
void DrawDCAxyWidthVsPt(const Source* sources, int nSrc)
{
    auto* c = new TCanvas("cDCAxyWidth", "", 700, 600);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.12);
    c->SetLogx();

    auto* leg = new TLegend(0.55, 0.65, 0.92, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);

    bool first = true;
    for (int s = 0; s < nSrc; ++s) {
        TFile* f = OpenFile(sources[s].path);
        if (!f) continue;
        auto* h2 = GetHist<TH2>(f, "h2_trackselplot_pt_dcaxy");
        if (!h2) { delete f; continue; }

        std::vector<double> ptC(kNPtSlices), ptE(kNPtSlices);
        std::vector<double> sig(kNPtSlices), sigE(kNPtSlices);
        int n = 0;
        for (int i = 0; i < kNPtSlices; ++i) {
            int binLo = h2->GetXaxis()->FindBin(kPtSlicesLo[i] + 1e-6);
            int binHi = h2->GetXaxis()->FindBin(kPtSlicesHi[i] - 1e-6);
            TH1* hp = h2->ProjectionY(Form("hDCAxy_w_%d_%d", s, i), binLo, binHi);
            if (hp->GetEntries() < 50) { delete hp; continue; }
            ptC[n]  = 0.5 * (kPtSlicesLo[i] + kPtSlicesHi[i]);
            ptE[n]  = 0.5 * (kPtSlicesHi[i] - kPtSlicesLo[i]);
            auto [fitSig, fitSigErr] = FitCoreSigma(hp);
            sig[n]  = fitSig;
            sigE[n] = fitSigErr;
            n++;
            delete hp;
        }
        f->Close(); delete f;
        if (n < 2) continue;

        auto* gr = new TGraphErrors(n, ptC.data(), sig.data(), ptE.data(), sigE.data());
        ApplyStyle(gr, sources[s]);

        if (first) {
            gr->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            gr->GetYaxis()->SetTitle("#sigma_{Gauss}(DCA_{xy}) (cm)");
            gr->Draw("AP");
            first = false;
        } else {
            gr->Draw("P same");
        }
        leg->AddEntry(gr, sources[s].label, "lp");
    }
    if (first) { delete c; return; }

    leg->Draw();
    TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035);
    tex.DrawLatex(0.16, 0.92, "#sigma_{Gauss}(DCA_{xy}) vs #it{p}_{T}");

    c->SaveAs(Form("%s/DCAxy_width_vs_pT.pdf", kOutputDir));
    delete c;
}

// ============================================================
//  8) Draw DCA_z sigma vs pT — all sources overlaid
// ============================================================
void DrawDCAzWidthVsPt(const Source* sources, int nSrc)
{
    auto* c = new TCanvas("cDCAzWidth", "", 700, 600);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.12);
    c->SetLogx();

    auto* leg = new TLegend(0.55, 0.65, 0.92, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);

    bool first = true;
    for (int s = 0; s < nSrc; ++s) {
        TFile* f = OpenFile(sources[s].path);
        if (!f) continue;
        auto* h2 = GetHist<TH2>(f, "h2_trackselplot_pt_dcaz");
        if (!h2) { delete f; continue; }

        std::vector<double> ptC(kNPtSlices), ptE(kNPtSlices);
        std::vector<double> sig(kNPtSlices), sigE(kNPtSlices);
        int n = 0;
        for (int i = 0; i < kNPtSlices; ++i) {
            int binLo = h2->GetXaxis()->FindBin(kPtSlicesLo[i] + 1e-6);
            int binHi = h2->GetXaxis()->FindBin(kPtSlicesHi[i] - 1e-6);
            TH1* hp = h2->ProjectionY(Form("hDCAz_w_%d_%d", s, i), binLo, binHi);
            if (hp->GetEntries() < 50) { delete hp; continue; }
            ptC[n]  = 0.5 * (kPtSlicesLo[i] + kPtSlicesHi[i]);
            ptE[n]  = 0.5 * (kPtSlicesHi[i] - kPtSlicesLo[i]);
            auto [fitSig, fitSigErr] = FitCoreSigma(hp);
            sig[n]  = fitSig;
            sigE[n] = fitSigErr;
            n++;
            delete hp;
        }
        f->Close(); delete f;
        if (n < 2) continue;

        auto* gr = new TGraphErrors(n, ptC.data(), sig.data(), ptE.data(), sigE.data());
        ApplyStyle(gr, sources[s]);

        if (first) {
            gr->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            gr->GetYaxis()->SetTitle("#sigma_{Gauss}(DCA_{z}) (cm)");
            gr->Draw("AP");
            first = false;
        } else {
            gr->Draw("P same");
        }
        leg->AddEntry(gr, sources[s].label, "lp");
    }
    if (first) { delete c; return; }

    leg->Draw();
    TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035);
    tex.DrawLatex(0.16, 0.92, "#sigma_{Gauss}(DCA_{z}) vs #it{p}_{T}");

    c->SaveAs(Form("%s/DCAz_width_vs_pT.pdf", kOutputDir));
    delete c;
}

// ============================================================
//  9) DCA_z in pT slices with cut lines — per source
//     Shows where 0.1, 0.2, 0.5 cm cuts fall on the distribution
// ============================================================
void DrawDCAzCutLines(const Source& src)
{
    TFile* f = OpenFile(src.path);
    if (!f) return;
    auto* h2 = GetHist<TH2>(f, "h2_trackselplot_pt_dcaz");
    if (!h2) { delete f; return; }

    auto* c = new TCanvas(Form("cDCAzCuts_%s", src.label), "", 1200, 800);
    c->Divide(4, 2);

    for (int i = 0; i < kNPtSlices; ++i) {
        c->cd(i + 1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.14);
        gPad->SetBottomMargin(0.13);

        int binLo = h2->GetXaxis()->FindBin(kPtSlicesLo[i] + 1e-6);
        int binHi = h2->GetXaxis()->FindBin(kPtSlicesHi[i] - 1e-6);
        TH1* hp = h2->ProjectionY(Form("hDCAzCut_%s_%d", src.label, i), binLo, binHi);
        if (hp->GetEntries() < 1) continue;
        hp->SetDirectory(0);

        hp->SetLineColor(kBlack);
        hp->SetLineWidth(1);
        hp->GetXaxis()->SetTitle("DCA_{z} (cm)");
        hp->GetYaxis()->SetTitle("Counts");
        hp->GetXaxis()->SetRangeUser(-1.0, 1.0);
        hp->SetTitle(Form("%.1f < p_{T} < %.0f GeV/c", kPtSlicesLo[i], kPtSlicesHi[i]));
        gStyle->SetOptTitle(1);
        hp->Draw("hist");

        // Draw cut lines
        double ymin = hp->GetMinimum() > 0 ? hp->GetMinimum() : 0.5;
        double ymax = hp->GetMaximum() * 2;
        for (int j = 0; j < kNDCAzCuts; ++j) {
            auto* lineP = new TLine( kDCAzCuts[j], ymin,  kDCAzCuts[j], ymax);
            auto* lineN = new TLine(-kDCAzCuts[j], ymin, -kDCAzCuts[j], ymax);
            lineP->SetLineColor(kDCAzCutColors[j]);
            lineN->SetLineColor(kDCAzCutColors[j]);
            lineP->SetLineStyle(kDCAzCutStyles[j]);
            lineN->SetLineStyle(kDCAzCutStyles[j]);
            lineP->SetLineWidth(2);
            lineN->SetLineWidth(2);
            lineP->Draw();
            lineN->Draw();
        }

        // Fraction outside each cut
        double total = hp->Integral(1, hp->GetNbinsX());
        if (total > 0) {
            TLatex tex; tex.SetNDC(); tex.SetTextSize(0.045);
            for (int j = 0; j < kNDCAzCuts; ++j) {
                int bLo = hp->FindBin(-kDCAzCuts[j] + 1e-6);
                int bHi = hp->FindBin( kDCAzCuts[j] - 1e-6);
                double inside = hp->Integral(bLo, bHi);
                double fracOut = (total - inside) / total * 100;
                tex.SetTextColor(kDCAzCutColors[j]);
                tex.DrawLatex(0.50, 0.85 - j * 0.07,
                    Form("|cut|=%.1f: %.2f%% out", kDCAzCuts[j], fracOut));
            }
        }
    }
    gStyle->SetOptTitle(0);

    c->SaveAs(Form("%s/DCAz_CutLines_%s.pdf", kOutputDir, src.label));
    delete c;
    f->Close(); delete f;
}

// ============================================================
//  10) Fraction of tracks outside |DCA_z| cut vs pT
//      = secondary contamination estimate at each cut
// ============================================================
void DrawDCAzTailFractionVsPt(const Source* sources, int nSrc)
{
    auto* c = new TCanvas("cDCAzTail", "", 800, 600);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.12);
    c->SetLogx();
    c->SetLogy();
    c->SetGridy();

    auto* leg = new TLegend(0.50, 0.60, 0.92, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.025);

    bool first = true;
    int markerIdx = 0;
    const int markers[] = {20, 21, 22, 23, 24, 25, 26, 27};

    for (int s = 0; s < nSrc; ++s) {
        TFile* f = OpenFile(sources[s].path);
        if (!f) continue;
        auto* h2 = GetHist<TH2>(f, "h2_trackselplot_pt_dcaz");
        if (!h2) { delete f; continue; }

        for (int j = 0; j < kNDCAzCuts; ++j) {
            std::vector<double> ptC(kNPtSlices), ptE(kNPtSlices);
            std::vector<double> frac(kNPtSlices), fracE(kNPtSlices);
            int n = 0;
            for (int i = 0; i < kNPtSlices; ++i) {
                int binLo = h2->GetXaxis()->FindBin(kPtSlicesLo[i] + 1e-6);
                int binHi = h2->GetXaxis()->FindBin(kPtSlicesHi[i] - 1e-6);
                TH1* hp = h2->ProjectionY(Form("hDCAzTail_%d_%d_%d", s, j, i), binLo, binHi);
                double total = hp->Integral(1, hp->GetNbinsX());
                if (total < 100) { delete hp; continue; }

                int cutBinLo = hp->FindBin(-kDCAzCuts[j] + 1e-6);
                int cutBinHi = hp->FindBin( kDCAzCuts[j] - 1e-6);
                double inside = hp->Integral(cutBinLo, cutBinHi);
                double outside = total - inside;

                ptC[n]   = 0.5 * (kPtSlicesLo[i] + kPtSlicesHi[i]);
                ptE[n]   = 0.5 * (kPtSlicesHi[i] - kPtSlicesLo[i]);
                frac[n]  = outside / total * 100; // percent
                fracE[n] = sqrt(outside * (1 - outside/total)) / total * 100; // binomial error
                n++;
                delete hp;
            }
            if (n < 2) continue;

            auto* gr = new TGraphErrors(n, ptC.data(), frac.data(), ptE.data(), fracE.data());
            gr->SetLineColor(kDCAzCutColors[j]);
            gr->SetMarkerColor(kDCAzCutColors[j]);
            gr->SetMarkerStyle(markers[markerIdx % 8]);
            gr->SetMarkerSize(0.9);
            gr->SetLineWidth(2);
            gr->SetLineStyle(sources[s].color == kBlack ? 1 : 2); // solid for first source

            if (first) {
                gr->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                gr->GetYaxis()->SetTitle("Fraction outside |DCA_{z}| cut (%)");
                gr->GetYaxis()->SetRangeUser(0.01, 50);
                gr->Draw("AP");
                first = false;
            } else {
                gr->Draw("P same");
            }
            leg->AddEntry(gr, Form("%s, |DCA_{z}| < %.1f cm", sources[s].label, kDCAzCuts[j]), "lp");
            markerIdx++;
        }
        f->Close(); delete f;
    }
    if (first) { delete c; return; }

    leg->Draw();
    TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035);
    tex.DrawLatex(0.16, 0.92, "Track fraction outside |DCA_{z}| cut vs #it{p}_{T}");

    c->SaveAs(Form("%s/DCAz_TailFraction_vs_pT.pdf", kOutputDir));
    delete c;
}

// ============================================================
//  Main entry point
// ============================================================
void DrawTrackDCA()
{
    SetStyle();
    gSystem->mkdir(kOutputDir, true);

    // --- Switch each plot on/off ---
    bool doDCAxy1D            = true;  // 1D DCA_xy overlay (all sources)
    bool doDCAz1D             = true;  // 1D DCA_z overlay (all sources)
    bool doDCAxy2D            = true;  // 2D colz per source
    bool doDCAz2D             = true;  // 2D colz per source
    bool doDCAxyPtSlices      = true;  // pT-slice projections per source
    bool doDCAzPtSlices       = true;  // pT-slice projections per source
    bool doDCAxyWidthVsPt     = true;  // sigma(DCA_xy) vs pT overlay
    bool doDCAzWidthVsPt      = true;  // sigma(DCA_z)  vs pT overlay
    bool doDCAzCutLines       = true;  // DCA_z dist per pT slice with cut lines
    bool doDCAzTailFraction   = true;  // fraction outside DCA_z cut vs pT

    // 1D & sigma-vs-pT: overlay all sources for comparison
    if (doDCAxy1D)        DrawDCAxy1D(kAllSources, kNAll);
    if (doDCAz1D)         DrawDCAz1D(kAllSources, kNAll);
    if (doDCAxyWidthVsPt) DrawDCAxyWidthVsPt(kAllSources, kNAll);
    if (doDCAzWidthVsPt)  DrawDCAzWidthVsPt(kAllSources, kNAll);

    // 2D colz: one canvas per source (can't overlay colz)
    if (doDCAxy2D) DrawDCAxy2D(kAllSources, kNAll);
    if (doDCAz2D)  DrawDCAz2D(kAllSources, kNAll);

    // pT slices: one canvas per source
    if (doDCAxyPtSlices || doDCAzPtSlices) {
        for (int i = 0; i < kNAll; ++i) {
            if (doDCAxyPtSlices) DrawDCAxyPtSlices(kAllSources[i]);
            if (doDCAzPtSlices)  DrawDCAzPtSlices(kAllSources[i]);
        }
    }

    // Secondary contamination study: DCA_z cut variation
    if (doDCAzCutLines) {
        for (int i = 0; i < kNAll; ++i)
            DrawDCAzCutLines(kAllSources[i]);
    }
    if (doDCAzTailFraction) DrawDCAzTailFractionVsPt(kAllSources, kNAll);

    std::cout << "All plots saved to " << kOutputDir << "/" << std::endl;
}
