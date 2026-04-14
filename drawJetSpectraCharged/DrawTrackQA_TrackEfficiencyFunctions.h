///////////////////////////////////////////////////
// DrawTrackQA_TrackEfficiencyFunctions.h
// All drawing functions for Track QA (track-efficiency)
// 3 processes: TrackSelection, Efficiency, DCA
///////////////////////////////////////////////////

#ifndef DRAWTRACKQAFUNCTIONS_H
#define DRAWTRACKQAFUNCTIONS_H

// ============================================================
// Process 1: DrawTrackSelection
//   1D overlays, 2D COLZ side-by-side, profiles vs pT
// ============================================================

// --- 1a: 1D overlays (Data vs MC) with ratio ---
void DrawTrackSel1D(
    const std::vector<TString>& fileNames,
    const std::vector<TString>& histNames,
    const std::vector<TString>& McFileDirs,
    const std::vector<Color_t>& Colors,
    TString outputDir, TString refPath, TString refDir, TString DataLabel)
{
    const int nMC = fileNames.size();
    int nn = 200;  // offset for unique Filipad2 IDs
    std::vector<TrackQAVar> vars = GetTrackQAVariables();

    for (auto& v : vars) {
        // Load Data
        TFile* fData = trkqa_OpenFile(refPath);
        if (!fData) continue;
        TDirectory* dirData = trkqa_GetDir(fData, refDir);
        TH1* hData = dirData ? GetHist<TH1>(dirData, v.hist1D.Data()) : nullptr;
        fData->Close(); delete fData;
        if (!hData) continue;

        // Load MC
        std::vector<TH1*> hMC(nMC, nullptr);
        for (int i = 0; i < nMC; i++) {
            TFile* fMC = trkqa_OpenFile(trkqa_mainDir + fileNames[i]);
            if (!fMC) continue;
            TDirectory* dirMC = trkqa_GetDir(fMC, McFileDirs[i]);
            if (dirMC) hMC[i] = GetHist<TH1>(dirMC, v.hist1D.Data());
            fMC->Close(); delete fMC;
        }

        // Normalize to unit area
        if (hData->Integral() > 0) hData->Scale(1.0 / hData->Integral());
        for (int i = 0; i < nMC; i++)
            if (hMC[i] && hMC[i]->Integral() > 0) hMC[i]->Scale(1.0 / hMC[i]->Integral());

        // Create Filipad2 canvas
        Filipad2* fpad = new Filipad2(Form("cSel1D_%s", v.shortName.Data()), ++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
        fpad->Draw();

        // Upper pad
        TPad* p1 = fpad->GetPad(1);
        p1->cd();
        p1->SetTickx(); p1->SetTicky();
        if (v.logY1D) p1->SetLogy();

        hset(*hData, "", "Normalized",
             2.5, 1.5, 0.05, 0.05, 0.01, 0.001, 0.03, 0.03, 510, 510);
        hData->GetXaxis()->SetRangeUser(v.xMin, v.xMax);
        hData->SetMarkerStyle(kDataMarker);
        hData->SetMarkerSize(kMarkerSize);
        hData->SetMarkerColor(Colors[0]);
        hData->SetLineColor(Colors[0]);
        hData->SetStats(0);

        double ymax = hData->GetMaximum();
        for (int i = 0; i < nMC; i++)
            if (hMC[i]) ymax = TMath::Max(ymax, hMC[i]->GetMaximum());
        if (v.logY1D) { hData->SetMinimum(1e-6); hData->SetMaximum(ymax * 5); }
        else          { hData->SetMinimum(0);     hData->SetMaximum(ymax * 1.5); }

        hData->DrawCopy("PE");
        for (int i = 0; i < nMC; i++) {
            if (!hMC[i]) continue;
            hMC[i]->SetMarkerStyle(kMCMarkers[i % 6]);
            hMC[i]->SetMarkerSize(kMarkerSize);
            hMC[i]->SetMarkerColor(Colors[i+1]);
            hMC[i]->SetLineColor(Colors[i+1]);
            hMC[i]->DrawCopy("PE SAME");
        }

        TLegend* leg = new TLegend(0.50, 0.65, 0.92, 0.92);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
        leg->AddEntry(hData, DataLabel, "lp");
        for (int i = 0; i < nMC; i++)
            if (hMC[i]) leg->AddEntry(hMC[i], histNames[i], "lp");
        leg->Draw();

        TLatex latex; latex.SetNDC(); latex.SetTextSize(0.05);
        latex.DrawLatex(0.20, 0.85, v.xTitle);

        // Lower pad: ratio
        TPad* p2 = fpad->GetPad(2);
        p2->cd();
        p2->SetTickx(); p2->SetTicky();
        p2->SetGridy();

        bool firstR = true;
        for (int i = 0; i < nMC; i++) {
            if (!hMC[i]) continue;
            TH1* hRatio = (TH1*)hMC[i]->Clone(Form("hR1D_%s_%d", v.shortName.Data(), i));
            hRatio->Divide(hData);
            hRatio->SetDirectory(0);
            if (firstR) {
                hset(*hRatio, v.xTitle, "MC/Data",
                     2.5, 0.8, 0.10, 0.10, 0.01, 0.001, 0.05, 0.03, 510, 510);
                hRatio->GetXaxis()->SetRangeUser(v.xMin, v.xMax);
                hRatio->SetMinimum(v.ratioMin);
                hRatio->SetMaximum(v.ratioMax);
                hRatio->DrawCopy("PE");
                firstR = false;
            } else {
                hRatio->DrawCopy("PE SAME");
            }
            delete hRatio;
        }
        TLine* line = new TLine(v.xMin, 1.0, v.xMax, 1.0);
        line->SetLineStyle(2); line->SetLineColor(kGray+2);
        line->Draw();

        fpad->C->SaveAs(Form("%s/TrackSel_1D_%s.pdf", outputDir.Data(), v.shortName.Data()));

        delete hData;
        for (int i = 0; i < nMC; i++) if (hMC[i]) delete hMC[i];
    }
}

// --- 1b: 2D COLZ side-by-side (Data | MC1 | MC2 | ...) ---
void DrawTrackSel2D(
    const std::vector<TString>& fileNames,
    const std::vector<TString>& histNames,
    const std::vector<TString>& McFileDirs,
    const std::vector<Color_t>& Colors,
    TString outputDir, TString refPath, TString refDir, TString DataLabel)
{
    const int nMC = fileNames.size();
    std::vector<TrackQAVar> vars = GetTrackQAVariables();

    for (auto& v : vars) {
        int nPanels = 1 + nMC;  // Data + MC files
        int cWidth = 500 * nPanels;
        if (cWidth > 2400) cWidth = 2400;

        TCanvas* c = new TCanvas(Form("c2D_%s", v.shortName.Data()),
                                  Form("2D %s", v.shortName.Data()), cWidth, 500);
        c->Divide(nPanels, 1);

        // Data panel
        c->cd(1);
        gPad->SetLogx(); gPad->SetLogz();
        gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.15);

        TFile* fData = trkqa_OpenFile(refPath);
        if (fData) {
            TDirectory* dirData = trkqa_GetDir(fData, refDir);
            if (dirData) {
                TH2* h2 = GetHist<TH2>(dirData, v.hist2D.Data());
                if (h2) {
                    h2->SetTitle(Form("Data: %s", v.xTitle.Data()));
                    h2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                    h2->GetYaxis()->SetTitle(v.xTitle);
                    h2->GetXaxis()->SetRangeUser(0.15, kPtMax);
                    h2->GetYaxis()->SetRangeUser(v.xMin, v.xMax);
                    h2->SetStats(0);
                    h2->DrawCopy("COLZ");
                    delete h2;
                }
            }
            fData->Close(); delete fData;
        }

        // MC panels
        for (int i = 0; i < nMC; i++) {
            c->cd(i + 2);
            gPad->SetLogx(); gPad->SetLogz();
            gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.15);

            TFile* fMC = trkqa_OpenFile(trkqa_mainDir + fileNames[i]);
            if (!fMC) continue;
            TDirectory* dirMC = trkqa_GetDir(fMC, McFileDirs[i]);
            if (dirMC) {
                TH2* h2 = GetHist<TH2>(dirMC, v.hist2D.Data());
                if (h2) {
                    h2->SetTitle(Form("%s: %s", histNames[i].Data(), v.xTitle.Data()));
                    h2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                    h2->GetYaxis()->SetTitle(v.xTitle);
                    h2->GetXaxis()->SetRangeUser(0.15, kPtMax);
                    h2->GetYaxis()->SetRangeUser(v.xMin, v.xMax);
                    h2->SetStats(0);
                    h2->DrawCopy("COLZ");
                    delete h2;
                }
            }
            fMC->Close(); delete fMC;
        }

        c->SaveAs(Form("%s/TrackSel_2D_%s.pdf", outputDir.Data(), v.shortName.Data()));
        delete c;
    }
}

// --- 1c: Profile plots (mean property vs pT) with ratio ---
void DrawTrackSelProfiles(
    const std::vector<TString>& fileNames,
    const std::vector<TString>& histNames,
    const std::vector<TString>& McFileDirs,
    const std::vector<Color_t>& Colors,
    TString outputDir, TString refPath, TString refDir, TString DataLabel)
{
    const int nMC = fileNames.size();
    int nn = 400;
    std::vector<TrackQAVar> vars = GetTrackQAVariables();

    for (auto& v : vars) {
        // Data profile
        TFile* fData = trkqa_OpenFile(refPath);
        if (!fData) continue;
        TDirectory* dirData = trkqa_GetDir(fData, refDir);
        TProfile* profData = nullptr;
        if (dirData) {
            TH2* h2 = GetHist<TH2>(dirData, v.hist2D.Data());
            if (h2) {
                profData = h2->ProfileX(Form("profData_%s", v.shortName.Data()));
                profData->SetDirectory(0);
                delete h2;
            }
        }
        fData->Close(); delete fData;
        if (!profData) continue;

        // MC profiles
        std::vector<TProfile*> profMC(nMC, nullptr);
        for (int i = 0; i < nMC; i++) {
            TFile* fMC = trkqa_OpenFile(trkqa_mainDir + fileNames[i]);
            if (!fMC) continue;
            TDirectory* dirMC = trkqa_GetDir(fMC, McFileDirs[i]);
            if (dirMC) {
                TH2* h2 = GetHist<TH2>(dirMC, v.hist2D.Data());
                if (h2) {
                    profMC[i] = h2->ProfileX(Form("profMC%d_%s", i, v.shortName.Data()));
                    profMC[i]->SetDirectory(0);
                    delete h2;
                }
            }
            fMC->Close(); delete fMC;
        }

        // Canvas
        Filipad2* fpad = new Filipad2(Form("cProf_%s", v.shortName.Data()), ++nn, 2, 0.4, 100, 50, 0.7, 1, 1);
        fpad->Draw();

        TPad* p1 = fpad->GetPad(1);
        p1->cd();
        p1->SetTickx(); p1->SetTicky();
        p1->SetLogx();

        hset(*profData, "", Form("#LT%s#GT", v.xTitle.Data()),
             2.5, 1.5, 0.05, 0.05, 0.01, 0.001, 0.03, 0.03, 510, 510);
        profData->GetXaxis()->SetRangeUser(0.15, kPtMax);
        profData->SetMarkerStyle(kDataMarker);
        profData->SetMarkerSize(kMarkerSize);
        profData->SetMarkerColor(Colors[0]);
        profData->SetLineColor(Colors[0]);
        profData->SetStats(0);
        profData->DrawCopy("PE");

        for (int i = 0; i < nMC; i++) {
            if (!profMC[i]) continue;
            profMC[i]->SetMarkerStyle(kMCMarkers[i % 6]);
            profMC[i]->SetMarkerSize(kMarkerSize);
            profMC[i]->SetMarkerColor(Colors[i+1]);
            profMC[i]->SetLineColor(Colors[i+1]);
            profMC[i]->DrawCopy("PE SAME");
        }

        TLegend* leg = new TLegend(0.50, 0.65, 0.92, 0.92);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
        leg->AddEntry(profData, DataLabel, "lp");
        for (int i = 0; i < nMC; i++)
            if (profMC[i]) leg->AddEntry(profMC[i], histNames[i], "lp");
        leg->Draw();

        TLatex latex; latex.SetNDC(); latex.SetTextSize(0.05);
        latex.DrawLatex(0.20, 0.85, Form("#LT%s#GT vs #it{p}_{T}", v.xTitle.Data()));

        // Ratio pad
        TPad* p2 = fpad->GetPad(2);
        p2->cd();
        p2->SetTickx(); p2->SetTicky();
        p2->SetLogx(); p2->SetGridy();

        bool firstR = true;
        for (int i = 0; i < nMC; i++) {
            if (!profMC[i]) continue;
            TH1* hR = (TH1*)profMC[i]->Clone(Form("hRProf_%s_%d", v.shortName.Data(), i));
            hR->Divide(profData);
            hR->SetDirectory(0);
            if (firstR) {
                hset(*hR, "#it{p}_{T} (GeV/#it{c})", "MC/Data",
                     2.5, 0.8, 0.10, 0.10, 0.01, 0.001, 0.05, 0.03, 510, 510);
                hR->GetXaxis()->SetRangeUser(0.15, kPtMax);
                hR->SetMinimum(0.9);
                hR->SetMaximum(1.1);
                hR->DrawCopy("PE");
                firstR = false;
            } else {
                hR->DrawCopy("PE SAME");
            }
            delete hR;
        }
        TLine* line = new TLine(0.15, 1.0, kPtMax, 1.0);
        line->SetLineStyle(2); line->SetLineColor(kGray+2);
        line->Draw();

        fpad->C->SaveAs(Form("%s/TrackSel_Profile_%s.pdf", outputDir.Data(), v.shortName.Data()));

        delete profData;
        for (int i = 0; i < nMC; i++) if (profMC[i]) delete profMC[i];
    }
}

// Top-level dispatcher for Process 1
void DrawTrackSelection(
    const std::vector<TString>& fileNames,
    const std::vector<TString>& histNames,
    const std::vector<TString>& McFileDirs,
    const std::vector<Color_t>& Colors,
    TString outputDir, TString refPath, TString refDir, TString DataLabel)
{
    // Check if trackselplot histograms exist in Data
    TFile* fTest = trkqa_OpenFile(refPath);
    if (!fTest) return;
    TDirectory* dirTest = trkqa_GetDir(fTest, refDir);
    bool hasTrackSel = false;
    if (dirTest) {
        TH1* hTest = (TH1*)dirTest->Get(kHistTPCRows);
        if (hTest) { hasTrackSel = true; delete hTest; }
    }
    fTest->Close(); delete fTest;

    if (!hasTrackSel) {
        std::cerr << "[Warning] Track selection histograms (h_trackselplot_*) not available." << std::endl;
        std::cerr << "  processTrackSelectionHistograms was not enabled in these train configs." << std::endl;
        return;
    }

    std::cerr << "[Info] Drawing track selection 1D overlays..." << std::endl;
    DrawTrackSel1D(fileNames, histNames, McFileDirs, Colors, outputDir, refPath, refDir, DataLabel);

    std::cerr << "[Info] Drawing track selection 2D comparisons..." << std::endl;
    DrawTrackSel2D(fileNames, histNames, McFileDirs, Colors, outputDir, refPath, refDir, DataLabel);

    std::cerr << "[Info] Drawing track selection profiles..." << std::endl;
    DrawTrackSelProfiles(fileNames, histNames, McFileDirs, Colors, outputDir, refPath, refDir, DataLabel);
}

// ============================================================
// Process 2: DrawTrackEfficiency
//   Efficiency, fake rate, secondary contamination
//   Multi-file MC overlay
// ============================================================
void DrawTrackEfficiency(
    const std::vector<TString>& fileNames,
    const std::vector<TString>& histNames,
    const std::vector<TString>& McFileDirs,
    const std::vector<Color_t>& Colors,
    TString outputDir)
{
    const int nMC = fileNames.size();
    if (nMC == 0) return;

    // --- Tracking efficiency ---
    {
        TCanvas* cEff = new TCanvas("cEff", "Track Efficiency", 800, 600);
        cEff->SetTickx(); cEff->SetTicky();
        cEff->SetLogx();
        cEff->SetLeftMargin(0.12); cEff->SetBottomMargin(0.12);

        TLegend* leg = new TLegend(0.15, 0.15, 0.60, 0.15 + 0.05 * nMC);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.030);

        bool first = true;
        for (int i = 0; i < nMC; i++) {
            TFile* fMC = trkqa_OpenFile(trkqa_mainDir + fileNames[i]);
            if (!fMC) continue;
            TDirectory* dirMC = trkqa_GetDir(fMC, McFileDirs[i]);
            if (!dirMC) { fMC->Close(); delete fMC; continue; }

            TH1D* hEff = BuildCombinedEfficiency(dirMC, kHistTruthLo, kHistMatchLo,
                                                  kHistTruthHi, kHistMatchHi,
                                                  Form("eff_%d", i));
            fMC->Close(); delete fMC;
            if (!hEff) continue;

            hset(*hEff, "#it{p}_{T} (GeV/#it{c})", "Tracking efficiency",
                 1.2, 1.4, 0.04, 0.04, 0.01, 0.001, 0.03, 0.03, 510, 510);
            hEff->GetXaxis()->SetRangeUser(0.15, 100.0);
            hEff->SetMinimum(0.0);
            hEff->SetMaximum(1.1);
            hEff->SetMarkerStyle(kMCMarkers[i % 6]);
            hEff->SetMarkerSize(0.5);
            hEff->SetMarkerColor(Colors[i+1]);
            hEff->SetLineColor(Colors[i+1]);
            hEff->SetStats(0);
            hEff->DrawCopy(first ? "PE" : "PE SAME");
            first = false;
            leg->AddEntry(hEff, histNames[i], "lp");
            delete hEff;
        }
        if (!first) {
            TLine* l1 = new TLine(0.15, 1.0, 100.0, 1.0);
            l1->SetLineStyle(2); l1->SetLineColor(kGray+1); l1->Draw();
            leg->Draw();
            cEff->SaveAs(Form("%s/Efficiency_vs_pT.pdf", outputDir.Data()));
        }
        delete cEff;
    }

    // --- Fake rate ---
    {
        TCanvas* cFake = new TCanvas("cFake", "Fake Rate", 800, 600);
        cFake->SetLogx(); cFake->SetTickx(); cFake->SetTicky();
        cFake->SetLeftMargin(0.12); cFake->SetBottomMargin(0.12);

        TLegend* leg = new TLegend(0.50, 0.70, 0.92, 0.70 + 0.05 * nMC);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.030);

        bool first = true;
        for (int i = 0; i < nMC; i++) {
            TFile* fMC = trkqa_OpenFile(trkqa_mainDir + fileNames[i]);
            if (!fMC) continue;
            TDirectory* dirMC = trkqa_GetDir(fMC, McFileDirs[i]);
            if (!dirMC) { fMC->Close(); delete fMC; continue; }

            TH1* hFrac = BuildFractionFromH3(dirMC, kHistFake, kHistRecoPrim,
                                              Form("fake_%d", i));
            fMC->Close(); delete fMC;
            if (!hFrac) continue;

            hset(*hFrac, "#it{p}_{T} (GeV/#it{c})", "Fake track rate",
                 1.2, 1.4, 0.04, 0.04, 0.01, 0.001, 0.03, 0.03, 510, 510);
            hFrac->GetXaxis()->SetRangeUser(0.15, 10.0);
            hFrac->SetMinimum(0.0);
            hFrac->SetMaximum(0.15);
            hFrac->SetMarkerStyle(kMCMarkers[i % 6]);
            hFrac->SetMarkerSize(0.5);
            hFrac->SetMarkerColor(Colors[i+1]);
            hFrac->SetLineColor(Colors[i+1]);
            hFrac->SetStats(0);
            hFrac->DrawCopy(first ? "PE" : "PE SAME");
            first = false;
            leg->AddEntry(hFrac, histNames[i], "lp");
            delete hFrac;
        }
        if (!first) {
            leg->Draw();
            cFake->SaveAs(Form("%s/FakeRate_vs_pT.pdf", outputDir.Data()));
        }
        delete cFake;
    }

    // --- Secondary contamination ---
    {
        TCanvas* cSec = new TCanvas("cSec", "Secondary Fraction", 800, 600);
        cSec->SetLogx(); cSec->SetTickx(); cSec->SetTicky();
        cSec->SetLeftMargin(0.12); cSec->SetBottomMargin(0.12);

        TLegend* leg = new TLegend(0.50, 0.70, 0.92, 0.70 + 0.05 * nMC);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.030);

        bool first = true;
        for (int i = 0; i < nMC; i++) {
            TFile* fMC = trkqa_OpenFile(trkqa_mainDir + fileNames[i]);
            if (!fMC) continue;
            TDirectory* dirMC = trkqa_GetDir(fMC, McFileDirs[i]);
            if (!dirMC) { fMC->Close(); delete fMC; continue; }

            TH1* hFrac = BuildFractionFromH3(dirMC, kHistRecoSec, kHistRecoPrim,
                                              Form("sec_%d", i));
            fMC->Close(); delete fMC;
            if (!hFrac) continue;

            hset(*hFrac, "#it{p}_{T} (GeV/#it{c})", "Secondary fraction",
                 1.2, 1.4, 0.04, 0.04, 0.01, 0.001, 0.03, 0.03, 510, 510);
            hFrac->GetXaxis()->SetRangeUser(0.15, 10.0);
            hFrac->SetMinimum(0.0);
            hFrac->SetMaximum(0.15);
            hFrac->SetMarkerStyle(kMCMarkers[i % 6]);
            hFrac->SetMarkerSize(0.5);
            hFrac->SetMarkerColor(Colors[i+1]);
            hFrac->SetLineColor(Colors[i+1]);
            hFrac->SetStats(0);
            hFrac->DrawCopy(first ? "PE" : "PE SAME");
            first = false;
            leg->AddEntry(hFrac, histNames[i], "lp");
            delete hFrac;
        }
        if (!first) {
            leg->Draw();
            cSec->SaveAs(Form("%s/SecondaryFraction_vs_pT.pdf", outputDir.Data()));
        }
        delete cSec;
    }

    std::cerr << "[Info] DrawTrackEfficiency done" << std::endl;
}

// ============================================================
// Process 3: DrawDCAAnalysis
//   1D DCA overlays, pT slices, sigma vs pT, tail fractions
//   All files (Data + MC) treated as overlay sources
// ============================================================

// --- 3a: 1D DCA overlay (all sources) ---
void DrawDCA1D(const char* histName, const char* dcaLabel, const char* shortName,
               double xRange,
               const std::vector<TString>& fileNames,
               const std::vector<TString>& histNames,
               const std::vector<TString>& McFileDirs,
               const std::vector<Color_t>& Colors,
               TString outputDir, TString refPath, TString refDir, TString DataLabel)
{
    const int nMC = fileNames.size();
    auto* c = new TCanvas(Form("cDCA1D_%s", shortName), "", 700, 600);
    c->SetLogy();
    c->SetLeftMargin(0.13); c->SetBottomMargin(0.12);

    auto* leg = new TLegend(0.55, 0.65, 0.92, 0.92);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.028);

    bool first = true;

    // Data
    TFile* fData = trkqa_OpenFile(refPath);
    if (fData) {
        TDirectory* dirData = trkqa_GetDir(fData, refDir);
        if (dirData) {
            TH1* h = GetHist<TH1>(dirData, histName);
            if (h) {
                auto [sig, sigErr] = FitCoreSigma(h);
                if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
                h->SetLineColor(Colors[0]); h->SetMarkerColor(Colors[0]);
                h->SetMarkerStyle(kDataMarker); h->SetMarkerSize(0.9);
                h->SetLineWidth(2);
                h->GetXaxis()->SetTitle(Form("%s (cm)", dcaLabel));
                h->GetYaxis()->SetTitle("Self-normalized");
                h->GetXaxis()->SetRangeUser(-xRange, xRange);
                h->Draw("hist");
                first = false;
                leg->AddEntry(h, Form("%s, #sigma = %.4f", DataLabel.Data(), sig), "l");
            }
        }
        fData->Close(); delete fData;
    }

    // MC files
    for (int i = 0; i < nMC; i++) {
        TFile* fMC = trkqa_OpenFile(trkqa_mainDir + fileNames[i]);
        if (!fMC) continue;
        TDirectory* dirMC = trkqa_GetDir(fMC, McFileDirs[i]);
        if (dirMC) {
            TH1* h = GetHist<TH1>(dirMC, histName);
            if (h) {
                auto [sig, sigErr] = FitCoreSigma(h);
                if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
                h->SetLineColor(Colors[i+1]); h->SetMarkerColor(Colors[i+1]);
                h->SetMarkerStyle(kMCMarkers[i % 6]); h->SetMarkerSize(0.9);
                h->SetLineWidth(2);
                h->Draw(first ? "hist" : "hist same");
                first = false;
                leg->AddEntry(h, Form("%s, #sigma = %.4f", histNames[i].Data(), sig), "l");
            }
        }
        fMC->Close(); delete fMC;
    }

    if (!first) {
        leg->Draw();
        TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035);
        tex.DrawLatex(0.16, 0.92, Form("%s,  pp #sqrt{#it{s}} = 13.6 TeV", dcaLabel));
        tex.SetTextSize(0.025);
        tex.DrawLatex(0.16, 0.87, "#sigma from Gaussian fit to core (#pm2#sigma)");
        c->SaveAs(Form("%s/%s_1D.pdf", outputDir.Data(), shortName));
    }
    delete c;
}

// --- 3b: DCA pT slices per source ---
void DrawDCAPtSlices(const char* hist2DName, const char* dcaLabel, const char* shortName,
                     double dcaRange,
                     const std::vector<TString>& fileNames,
                     const std::vector<TString>& histNames,
                     const std::vector<TString>& McFileDirs,
                     TString outputDir, TString refPath, TString refDir, TString DataLabel)
{
    // Lambda to draw pT slices for one file
    auto drawSlicesForFile = [&](TFile* f, const TString& dirName, const TString& label) {
        TDirectory* dir = trkqa_GetDir(f, dirName);
        if (!dir) return;
        TH2* h2 = GetHist<TH2>(dir, hist2DName);
        if (!h2) return;

        auto* c = new TCanvas(Form("cSlice_%s_%s", shortName, label.Data()), "", 800, 600);
        c->SetLogy();
        c->SetLeftMargin(0.13); c->SetBottomMargin(0.12);

        auto* leg = new TLegend(0.62, 0.50, 0.92, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.028);
        leg->SetHeader(label);

        bool first = true;
        for (int i = 0; i < kNPtSlices; i++) {
            int binLo = h2->GetXaxis()->FindBin(kPtSlicesLo[i] + 1e-6);
            int binHi = h2->GetXaxis()->FindBin(kPtSlicesHi[i] - 1e-6);
            TH1* hp = h2->ProjectionY(Form("hSlice_%s_%s_%d", shortName, label.Data(), i), binLo, binHi);
            if (hp->GetEntries() < 1) { delete hp; continue; }
            hp->SetDirectory(0);
            if (hp->Integral() > 0) hp->Scale(1.0 / hp->Integral());
            hp->SetLineColor(kSliceColors[i % 8]);
            hp->SetLineWidth(2);
            hp->GetXaxis()->SetTitle(Form("%s (cm)", dcaLabel));
            hp->GetYaxis()->SetTitle("Normalized");
            hp->GetXaxis()->SetRangeUser(-dcaRange, dcaRange);
            hp->Draw(first ? "hist" : "hist same");
            first = false;
            leg->AddEntry(hp, Form("%.1f < #it{p}_{T} < %.0f", kPtSlicesLo[i], kPtSlicesHi[i]), "l");
        }
        if (!first) {
            leg->Draw();
            TString cleanLabel = label;
            cleanLabel.ReplaceAll(" ", "_");
            cleanLabel.ReplaceAll("(", "");
            cleanLabel.ReplaceAll(")", "");
            c->SaveAs(Form("%s/%s_pTslices_%s.pdf", outputDir.Data(), shortName, cleanLabel.Data()));
        }
        delete h2; delete c;
    };

    // Data
    TFile* fData = trkqa_OpenFile(refPath);
    if (fData) {
        drawSlicesForFile(fData, refDir, DataLabel);
        fData->Close(); delete fData;
    }

    // MC files
    const int nMC = fileNames.size();
    for (int i = 0; i < nMC; i++) {
        TFile* fMC = trkqa_OpenFile(trkqa_mainDir + fileNames[i]);
        if (!fMC) continue;
        drawSlicesForFile(fMC, McFileDirs[i], histNames[i]);
        fMC->Close(); delete fMC;
    }
}

// --- 3c: DCA sigma vs pT (all sources overlaid) ---
void DrawDCAWidthVsPt(const char* hist2DName, const char* dcaLabel, const char* shortName,
                      const std::vector<TString>& fileNames,
                      const std::vector<TString>& histNames,
                      const std::vector<TString>& McFileDirs,
                      const std::vector<Color_t>& Colors,
                      TString outputDir, TString refPath, TString refDir, TString DataLabel)
{
    const int nMC = fileNames.size();
    auto* c = new TCanvas(Form("cWidth_%s", shortName), "", 700, 600);
    c->SetLeftMargin(0.13); c->SetBottomMargin(0.12); c->SetLogx();

    auto* leg = new TLegend(0.55, 0.65, 0.92, 0.92);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.030);

    bool first = true;

    // Lambda for one source
    auto processSource = [&](TFile* f, const TString& dirName, Color_t col, int marker, const TString& label) {
        TDirectory* dir = trkqa_GetDir(f, dirName);
        if (!dir) return;
        TH2* h2 = GetHist<TH2>(dir, hist2DName);
        if (!h2) return;

        std::vector<double> ptC(kNPtSlices), ptE(kNPtSlices);
        std::vector<double> sig(kNPtSlices), sigE(kNPtSlices);
        int n = 0;
        for (int i = 0; i < kNPtSlices; i++) {
            int binLo = h2->GetXaxis()->FindBin(kPtSlicesLo[i] + 1e-6);
            int binHi = h2->GetXaxis()->FindBin(kPtSlicesHi[i] - 1e-6);
            TH1* hp = h2->ProjectionY(Form("hW_%s_%s_%d", shortName, label.Data(), i), binLo, binHi);
            if (hp->GetEntries() < 50) { delete hp; continue; }
            ptC[n] = 0.5 * (kPtSlicesLo[i] + kPtSlicesHi[i]);
            ptE[n] = 0.5 * (kPtSlicesHi[i] - kPtSlicesLo[i]);
            auto [fitSig, fitSigErr] = FitCoreSigma(hp);
            sig[n] = fitSig;
            sigE[n] = fitSigErr;
            n++;
            delete hp;
        }
        delete h2;
        if (n < 2) return;

        auto* gr = new TGraphErrors(n, ptC.data(), sig.data(), ptE.data(), sigE.data());
        gr->SetLineColor(col); gr->SetMarkerColor(col);
        gr->SetMarkerStyle(marker); gr->SetMarkerSize(0.9);
        gr->SetLineWidth(2);

        if (first) {
            gr->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            gr->GetYaxis()->SetTitle(Form("#sigma_{Gauss}(%s) (cm)", dcaLabel));
            gr->Draw("AP");
            first = false;
        } else {
            gr->Draw("P same");
        }
        leg->AddEntry(gr, label, "lp");
    };

    // Data
    TFile* fData = trkqa_OpenFile(refPath);
    if (fData) {
        processSource(fData, refDir, Colors[0], kDataMarker, DataLabel);
        fData->Close(); delete fData;
    }

    // MC
    for (int i = 0; i < nMC; i++) {
        TFile* fMC = trkqa_OpenFile(trkqa_mainDir + fileNames[i]);
        if (!fMC) continue;
        processSource(fMC, McFileDirs[i], Colors[i+1], kMCMarkers[i % 6], histNames[i]);
        fMC->Close(); delete fMC;
    }

    if (!first) {
        leg->Draw();
        TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035);
        tex.DrawLatex(0.16, 0.92, Form("#sigma_{Gauss}(%s) vs #it{p}_{T}", dcaLabel));
        c->SaveAs(Form("%s/%s_width_vs_pT.pdf", outputDir.Data(), shortName));
    }
    delete c;
}

// --- 3d: DCA_z tail fraction vs pT ---
void DrawDCAzTailFraction(
    const std::vector<TString>& fileNames,
    const std::vector<TString>& histNames,
    const std::vector<TString>& McFileDirs,
    const std::vector<Color_t>& Colors,
    TString outputDir, TString refPath, TString refDir, TString DataLabel)
{
    const int nMC = fileNames.size();
    auto* c = new TCanvas("cDCAzTail", "", 800, 600);
    c->SetLeftMargin(0.13); c->SetBottomMargin(0.12);
    c->SetLogx(); c->SetLogy(); c->SetGridy();

    auto* leg = new TLegend(0.45, 0.55, 0.92, 0.92);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.023);

    bool first = true;
    int markerIdx = 0;
    const int markers[] = {20, 21, 22, 23, 24, 25, 26, 27};

    // Lambda for one source x all DCA_z cuts
    auto processSource = [&](TFile* f, const TString& dirName, const TString& label, int srcColor) {
        TDirectory* dir = trkqa_GetDir(f, dirName);
        if (!dir) return;
        TH2* h2 = GetHist<TH2>(dir, kHist2DDCAz);
        if (!h2) return;

        for (int j = 0; j < kNDCAzCuts; j++) {
            std::vector<double> ptC(kNPtSlices), ptE(kNPtSlices);
            std::vector<double> frac(kNPtSlices), fracE(kNPtSlices);
            int n = 0;
            for (int i = 0; i < kNPtSlices; i++) {
                int binLo = h2->GetXaxis()->FindBin(kPtSlicesLo[i] + 1e-6);
                int binHi = h2->GetXaxis()->FindBin(kPtSlicesHi[i] - 1e-6);
                TH1* hp = h2->ProjectionY(Form("hTail_%s_%d_%d", label.Data(), j, i), binLo, binHi);
                double total = hp->Integral(1, hp->GetNbinsX());
                if (total < 100) { delete hp; continue; }
                int cutBinLo = hp->FindBin(-kDCAzCuts[j] + 1e-6);
                int cutBinHi = hp->FindBin( kDCAzCuts[j] - 1e-6);
                double inside = hp->Integral(cutBinLo, cutBinHi);
                double outside = total - inside;
                ptC[n]   = 0.5 * (kPtSlicesLo[i] + kPtSlicesHi[i]);
                ptE[n]   = 0.5 * (kPtSlicesHi[i] - kPtSlicesLo[i]);
                frac[n]  = outside / total * 100;
                fracE[n] = sqrt(outside * (1 - outside/total)) / total * 100;
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
            gr->SetLineStyle(srcColor == (int)Colors[0] ? 1 : 2);

            if (first) {
                gr->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                gr->GetYaxis()->SetTitle("Fraction outside |DCA_{z}| cut (%)");
                gr->GetYaxis()->SetRangeUser(0.01, 50);
                gr->Draw("AP");
                first = false;
            } else {
                gr->Draw("P same");
            }
            leg->AddEntry(gr, Form("%s, |DCA_{z}| < %.1f cm", label.Data(), kDCAzCuts[j]), "lp");
            markerIdx++;
        }
        delete h2;
    };

    // Data
    TFile* fData = trkqa_OpenFile(refPath);
    if (fData) {
        processSource(fData, refDir, DataLabel, Colors[0]);
        fData->Close(); delete fData;
    }

    // MC
    for (int i = 0; i < nMC; i++) {
        TFile* fMC = trkqa_OpenFile(trkqa_mainDir + fileNames[i]);
        if (!fMC) continue;
        processSource(fMC, McFileDirs[i], histNames[i], Colors[i+1]);
        fMC->Close(); delete fMC;
    }

    if (!first) {
        leg->Draw();
        TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035);
        tex.DrawLatex(0.16, 0.92, "Track fraction outside |DCA_{z}| cut vs #it{p}_{T}");
        c->SaveAs(Form("%s/DCAz_TailFraction_vs_pT.pdf", outputDir.Data()));
    }
    delete c;
}

// Top-level dispatcher for Process 3
void DrawDCAAnalysis(
    const std::vector<TString>& fileNames,
    const std::vector<TString>& histNames,
    const std::vector<TString>& McFileDirs,
    const std::vector<Color_t>& Colors,
    TString outputDir, TString refPath, TString refDir, TString DataLabel)
{
    std::cerr << "[Info] Drawing DCA_xy 1D overlay..." << std::endl;
    DrawDCA1D(kHistDCAxy, "DCA_{xy}", "DCA_xy", 0.5,
              fileNames, histNames, McFileDirs, Colors, outputDir, refPath, refDir, DataLabel);

    std::cerr << "[Info] Drawing DCA_z 1D overlay..." << std::endl;
    DrawDCA1D(kHistDCAz, "DCA_{z}", "DCA_z", 2.0,
              fileNames, histNames, McFileDirs, Colors, outputDir, refPath, refDir, DataLabel);

    std::cerr << "[Info] Drawing DCA_xy pT slices..." << std::endl;
    DrawDCAPtSlices(kHist2DDCAxy, "DCA_{xy}", "DCA_xy", 0.3,
                    fileNames, histNames, McFileDirs, outputDir, refPath, refDir, DataLabel);

    std::cerr << "[Info] Drawing DCA_z pT slices..." << std::endl;
    DrawDCAPtSlices(kHist2DDCAz, "DCA_{z}", "DCA_z", 1.0,
                    fileNames, histNames, McFileDirs, outputDir, refPath, refDir, DataLabel);

    std::cerr << "[Info] Drawing DCA_xy width vs pT..." << std::endl;
    DrawDCAWidthVsPt(kHist2DDCAxy, "DCA_{xy}", "DCA_xy",
                     fileNames, histNames, McFileDirs, Colors, outputDir, refPath, refDir, DataLabel);

    std::cerr << "[Info] Drawing DCA_z width vs pT..." << std::endl;
    DrawDCAWidthVsPt(kHist2DDCAz, "DCA_{z}", "DCA_z",
                     fileNames, histNames, McFileDirs, Colors, outputDir, refPath, refDir, DataLabel);

    std::cerr << "[Info] Drawing DCA_z tail fraction..." << std::endl;
    DrawDCAzTailFraction(fileNames, histNames, McFileDirs, Colors, outputDir, refPath, refDir, DataLabel);

    std::cerr << "[Info] DrawDCAAnalysis done" << std::endl;
}

#endif
