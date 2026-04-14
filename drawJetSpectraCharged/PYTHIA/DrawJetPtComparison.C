// Macro to draw hJetPt comparison between two PYTHIA files
// Upper pad: yield plot
// Lower pad: ratio plot

void DrawJetPtComparison13600GeVto13000GeV() {
    // File paths
    const char* file1_path = "/Users/js/alice/pythiaGen/postprocess/results/pp_13600GeV_HardQCD_all_on_UE_ISR_FSR_on/PYTHIA_pp_13600_GeV.root";
    const char* file2_path = "/Users/js/alice/pythiaGen/postprocess/results/pp_13000_GeV_HardQCD_all_on_UE_ISR_FSR_on/PYTHIA_pp_13000_GeV.root";
    
    // Open ROOT files
    TFile* file1 = TFile::Open(file1_path, "READ");
    TFile* file2 = TFile::Open(file2_path, "READ");
    
    if (!file1 || file1->IsZombie()) {
        std::cerr << "Error: Cannot open file " << file1_path << std::endl;
        return;
    }
    if (!file2 || file2->IsZombie()) {
        std::cerr << "Error: Cannot open file " << file2_path << std::endl;
        file1->Close();
        return;
    }
    
    // Get histograms
    TH1D* h1 = (TH1D*)file1->Get("hJetPt");
    TH1D* h2 = (TH1D*)file2->Get("hJetPt");
    
    if (!h1) {
        std::cerr << "Error: Cannot find hJetPt in file1" << std::endl;
        file1->Close();
        file2->Close();
        return;
    }
    if (!h2) {
        std::cerr << "Error: Cannot find hJetPt in file2" << std::endl;
        file1->Close();
        file2->Close();
        return;
    }
    
    // Get nevent histograms
    TH1D* hnevent1 = (TH1D*)file1->Get("hnevent");
    TH1D* hnevent2 = (TH1D*)file2->Get("hnevent");
    
    if (!hnevent1) {
        std::cerr << "Error: Cannot find hnevent in file1" << std::endl;
        file1->Close();
        file2->Close();
        return;
    }
    if (!hnevent2) {
        std::cerr << "Error: Cannot find hnevent in file2" << std::endl;
        file1->Close();
        file2->Close();
        return;
    }
    
    // Get nevent values
    Double_t nevent1 = hnevent1->GetBinContent(1);
    Double_t nevent2 = hnevent2->GetBinContent(1);
    
    if (nevent1 <= 0 || nevent2 <= 0) {
        std::cerr << "Error: Invalid nevent values" << std::endl;
        file1->Close();
        file2->Close();
        return;
    }
    
    // Define binning
    const Double_t ptbinGen[27] = {0,  1,  2,  3,  4,  5,   6,   7,  8,
                                   9,  10, 12, 14, 16, 18,  20,  25, 30,
                                   40, 50, 60, 70, 85, 100, 140, 200, 300};
    
    // Create rebinned histograms with custom binning
    TH1D* h1_rebinned = new TH1D("h1_rebinned", "", 25, ptbinGen);
    TH1D* h2_rebinned = new TH1D("h2_rebinned", "", 25, ptbinGen);
    
    // Fill rebinned histograms from original
    for (Int_t i = 1; i <= h1->GetNbinsX(); i++) {
        Double_t x = h1->GetBinCenter(i);
        Double_t content = h1->GetBinContent(i);
        Double_t error = h1->GetBinError(i);
        Int_t bin = h1_rebinned->FindBin(x);
        if (bin > 0 && bin <= h1_rebinned->GetNbinsX()) {
            h1_rebinned->SetBinContent(bin, h1_rebinned->GetBinContent(bin) + content);
            Double_t oldError = h1_rebinned->GetBinError(bin);
            h1_rebinned->SetBinError(bin, TMath::Sqrt(oldError*oldError + error*error));
        }
    }
    
    for (Int_t i = 1; i <= h2->GetNbinsX(); i++) {
        Double_t x = h2->GetBinCenter(i);
        Double_t content = h2->GetBinContent(i);
        Double_t error = h2->GetBinError(i);
        Int_t bin = h2_rebinned->FindBin(x);
        if (bin > 0 && bin <= h2_rebinned->GetNbinsX()) {
            h2_rebinned->SetBinContent(bin, h2_rebinned->GetBinContent(bin) + content);
            Double_t oldError = h2_rebinned->GetBinError(bin);
            h2_rebinned->SetBinError(bin, TMath::Sqrt(oldError*oldError + error*error));
        }
    }
    
    // Scale by bin width and hardcoded event count (100000000)
    h1_rebinned->Scale(1./100000000, "width");
    h2_rebinned->Scale(1./100000000, "width");
    
    // Remove titles
    h1_rebinned->SetTitle("");
    h2_rebinned->SetTitle("");
    
    // Clone for ratio calculation (13 TeV / 13.6 TeV)
    TH1D* h1_ratio = (TH1D*)h1_rebinned->Clone("h1_ratio");
    TH1D* h2_ratio = (TH1D*)h2_rebinned->Clone("h2_ratio");
    h1_ratio->SetTitle("");
    h2_ratio->SetTitle("");
    
    // Create Filipad canvas
    TCanvas* c = new TCanvas("c", "Jet p_{T} Comparison", 800, 800);
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.05);
    pad1->SetTopMargin(0.05);
    pad1->Draw();
    pad1->cd();
    
    // Upper pad: yield plot
    pad1->SetLogy();
    h1_rebinned->SetLineColor(kRed);
    h1_rebinned->SetMarkerColor(kRed);
    h1_rebinned->SetMarkerStyle(20);
    h1_rebinned->SetMarkerSize(1.0);
    h1_rebinned->SetLineWidth(2);
    
    h2_rebinned->SetLineColor(kBlue);
    h2_rebinned->SetMarkerColor(kBlue);
    h2_rebinned->SetMarkerStyle(21);
    h2_rebinned->SetMarkerSize(1.0);
    h2_rebinned->SetLineWidth(2);
    
    // Set axis labels and titles
    h1_rebinned->GetXaxis()->SetTitle("#it{p}_{T, jet}^{ch} (GeV/c)");
    h1_rebinned->GetYaxis()->SetTitle("d^{2}#sigma/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]");
    h1_rebinned->GetYaxis()->SetTitleSize(0.06);
    h1_rebinned->GetYaxis()->SetTitleOffset(1.15);
    h1_rebinned->GetXaxis()->SetLabelSize(0);
    h1_rebinned->GetXaxis()->SetTitleSize(0);
    
    // Draw
    h1_rebinned->Draw("PE");
    h2_rebinned->Draw("PE SAME");
    
    // Legend
    TLegend* leg = new TLegend(0.60, 0.5, 0.95, 0.7);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h1_rebinned, "pp 13.6 TeV", "lp");
    leg->AddEntry(h2_rebinned, "pp 13 TeV", "lp");
    leg->Draw();
    
    // PYTHIA description
    TLatex* pythia_text = new TLatex();
    pythia_text->SetNDC();
    pythia_text->SetTextFont(42);
    pythia_text->SetTextSize(0.05);
    pythia_text->SetTextAlign(12);
    pythia_text->DrawLatex(0.35, 0.90, "PYTHIA 8 Monash 2013");
    pythia_text->DrawLatex(0.35, 0.85, "HardQCD:all, #tau_{0} < 10 mm/c");
    
    // Lower pad: ratio plot
    c->cd();
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.05);
    pad2->Draw();
    pad2->cd();
    
    // Calculate ratio: 13 TeV / 13.6 TeV (h2/h1)
    h2_ratio->Divide(h1_ratio);
    h2_ratio->SetStats(0);
    h2_ratio->SetLineColor(kBlue);
    h2_ratio->SetMarkerColor(kBlue);
    h2_ratio->SetMarkerStyle(21);
    h2_ratio->SetMarkerSize(1.0);
    h2_ratio->SetLineWidth(2);
    
    // Set axis labels
    h2_ratio->GetXaxis()->SetTitle("#it{p}_{T, jet}^{ch} (GeV/c)");
    h2_ratio->GetXaxis()->SetTitleSize(0.12);
    h2_ratio->GetXaxis()->SetLabelSize(0.1);
    h2_ratio->GetYaxis()->SetTitle("13 / 13.6 TeV");
    h2_ratio->GetYaxis()->SetTitleSize(0.12);
    h2_ratio->GetYaxis()->SetTitleOffset(0.6);
    h2_ratio->GetYaxis()->SetLabelSize(0.1);
    h2_ratio->GetYaxis()->SetNdivisions(505);
    
    // Set range for ratio
    h2_ratio->SetMinimum(0.5);
    h2_ratio->SetMaximum(1.5);
    h2_ratio->GetYaxis()->SetRangeUser(0.89, 1.02);
    
    // Draw ratio
    h2_ratio->Draw("PE");
    
    // Draw reference line at y=1
    TLine* line = new TLine(h2_ratio->GetXaxis()->GetXmin(), 1.0, 
                           h2_ratio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kRed);
    line->Draw("SAME");
    
    // Save
    c->cd();
    c->SaveAs("JetPtComparison.pdf");
    c->SaveAs("JetPtComparison.png");
    
    // Save histograms to ROOT file
    TFile* outfile = new TFile("JetPtComparison.root", "RECREATE");
    c->Write("canvas");
    h1_rebinned->Write();
    h2_rebinned->Write();
    h2_ratio->Write();
    outfile->Close();
    
    std::cout << "Plot saved as JetPtComparison.pdf, JetPtComparison.png, and JetPtComparison.root" << std::endl;
    
    // Cleanup
    // file1->Close();
    // file2->Close();

    return;
}


// R-dependent jet pT comparison functions
void DrawJetPtComparisonR() {
    // File paths
    const char* file1_path = "~/alice/pythiaGen/postprocess/results/results_20251123_195353_results_pp_13600_GeV_HardQCD_all_on_UE_ISR_FSR_on/combined.root";
    const char* file2_path = "~/alice/pythiaGen/postprocess/results/results_20251123_195255_results_pp_13000_GeV_HardQCD_all_on_UE_ISR_FSR_on/combined.root";
    
    // Expand user path
    TString file1_expanded = gSystem->ExpandPathName(file1_path);
    TString file2_expanded = gSystem->ExpandPathName(file2_path);
    
    // Open ROOT files
    TFile* file1 = TFile::Open(file1_expanded.Data(), "READ");
    TFile* file2 = TFile::Open(file2_expanded.Data(), "READ");
    
    if (!file1 || file1->IsZombie()) {
        std::cerr << "Error: Cannot open file " << file1_expanded.Data() << std::endl;
        return;
    }
    if (!file2 || file2->IsZombie()) {
        std::cerr << "Error: Cannot open file " << file2_expanded.Data() << std::endl;
        file1->Close();
        return;
    }
    
    // Get nevent histograms
    TH1D* hnevent1 = (TH1D*)file1->Get("hnevent");
    TH1D* hnevent2 = (TH1D*)file2->Get("hnevent");
    
    if (!hnevent1 || !hnevent2) {
        std::cerr << "Error: Cannot find hnevent histograms" << std::endl;
        file1->Close();
        file2->Close();
        return;
    }
    
    Double_t nevent1 = hnevent1->GetBinContent(1);
    Double_t nevent2 = hnevent2->GetBinContent(1);
    
    if (nevent1 <= 0 || nevent2 <= 0) {
        std::cerr << "Error: Invalid nevent values" << std::endl;
        file1->Close();
        file2->Close();
        return;
    }
    
    // Define binning
    const Double_t ptbinGen[27] = {0,  1,  2,  3,  4,  5,   6,   7,  8,
                                   9,  10, 12, 14, 16, 18,  20,  25, 30,
                                   40, 50, 60, 70, 85, 100, 140, 200, 300};
    
    // R values
    const Int_t nR = 7;
    const Double_t R_values[nR] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};
    const char* R_names[nR] = {"R01", "R02", "R03", "R04", "R05", "R06", "R07"};
    
    // Color scheme for different R values
    const Int_t colors[nR] = {kRed, kBlue, kGreen+2, kMagenta, kCyan+1, kOrange+1, kViolet};
    
    // Arrays to store histograms
    TH1D* h1_R[nR];  // 13.6 TeV
    TH1D* h2_R[nR];  // 13 TeV
    TH1D* h1_R_rebinned[nR];
    TH1D* h2_R_rebinned[nR];
    TH1D* h1_R_ratio[nR];  // 13.6/13 ratio for each R
    TH1D* h2_R_ratio[nR];  // 13.6/13 ratio for each R
    
    // Load and process histograms for each R
    for (Int_t iR = 0; iR < nR; iR++) {
        TString hist_name = TString::Format("hJetPt_%s", R_names[iR]);
        
        h1_R[iR] = (TH1D*)file1->Get(hist_name.Data());
        h2_R[iR] = (TH1D*)file2->Get(hist_name.Data());
        
        if (!h1_R[iR] || !h2_R[iR]) {
            std::cerr << "Error: Cannot find " << hist_name.Data() << std::endl;
            file1->Close();
            file2->Close();
            return;
        }
        
        // Create rebinned histograms
        h1_R_rebinned[iR] = new TH1D(TString::Format("h1_R%02d_rebinned", iR), "", 25, ptbinGen);
        h2_R_rebinned[iR] = new TH1D(TString::Format("h2_R%02d_rebinned", iR), "", 25, ptbinGen);
        
        // Fill rebinned histograms
        for (Int_t i = 1; i <= h1_R[iR]->GetNbinsX(); i++) {
            Double_t x = h1_R[iR]->GetBinCenter(i);
            Double_t content = h1_R[iR]->GetBinContent(i);
            Double_t error = h1_R[iR]->GetBinError(i);
            Int_t bin = h1_R_rebinned[iR]->FindBin(x);
            if (bin > 0 && bin <= h1_R_rebinned[iR]->GetNbinsX()) {
                h1_R_rebinned[iR]->SetBinContent(bin, h1_R_rebinned[iR]->GetBinContent(bin) + content);
                Double_t oldError = h1_R_rebinned[iR]->GetBinError(bin);
                h1_R_rebinned[iR]->SetBinError(bin, TMath::Sqrt(oldError*oldError + error*error));
            }
        }
        
        for (Int_t i = 1; i <= h2_R[iR]->GetNbinsX(); i++) {
            Double_t x = h2_R[iR]->GetBinCenter(i);
            Double_t content = h2_R[iR]->GetBinContent(i);
            Double_t error = h2_R[iR]->GetBinError(i);
            Int_t bin = h2_R_rebinned[iR]->FindBin(x);
            if (bin > 0 && bin <= h2_R_rebinned[iR]->GetNbinsX()) {
                h2_R_rebinned[iR]->SetBinContent(bin, h2_R_rebinned[iR]->GetBinContent(bin) + content);
                Double_t oldError = h2_R_rebinned[iR]->GetBinError(bin);
                h2_R_rebinned[iR]->SetBinError(bin, TMath::Sqrt(oldError*oldError + error*error));
            }
        }
        
        // Scale by bin width and hardcoded event count (100000000)
        h1_R_rebinned[iR]->Scale(1./100000000, "width");
        h2_R_rebinned[iR]->Scale(1./100000000, "width");
        
        h1_R_rebinned[iR]->SetTitle("");
        h2_R_rebinned[iR]->SetTitle("");
        
        // Clone for ratio calculation (13/13.6)
        h1_R_ratio[iR] = (TH1D*)h1_R_rebinned[iR]->Clone(TString::Format("h1_R%02d_ratio", iR));
        h2_R_ratio[iR] = (TH1D*)h2_R_rebinned[iR]->Clone(TString::Format("h2_R%02d_ratio", iR));
        h2_R_ratio[iR]->Divide(h1_R_ratio[iR]);  // 13/13.6 (inverse)
        h2_R_ratio[iR]->SetTitle("");
    }
    
    // ========== Plot 1: 14 cross sections + 7 ratio plots ==========
    TCanvas* c1 = new TCanvas("c1", "R-dependent Jet pT Comparison", 1000, 800);
    TPad* pad1_upper = new TPad("pad1_upper", "pad1_upper", 0, 0.3, 1, 1.0);
    pad1_upper->SetBottomMargin(0.02);
    pad1_upper->SetLeftMargin(0.15);
    pad1_upper->SetRightMargin(0.05);
    pad1_upper->SetTopMargin(0.05);
    pad1_upper->Draw();
    pad1_upper->cd();
    pad1_upper->SetLogy();
    
    // Draw 14 histograms (7 R values × 2 energies)
    TLegend* leg1 = new TLegend(0.55, 0.5, 0.95, 0.88);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.025);
    leg1->SetNColumns(2);
    
    Bool_t first_drawn = kFALSE;
    for (Int_t iR = 0; iR < nR; iR++) {
        // 13.6 TeV
        h1_R_rebinned[iR]->SetLineColor(colors[iR]);
        h1_R_rebinned[iR]->SetMarkerColor(colors[iR]);
        h1_R_rebinned[iR]->SetMarkerStyle(20);  // filled circle
        h1_R_rebinned[iR]->SetMarkerSize(0.8);
        h1_R_rebinned[iR]->SetLineWidth(2);
        h1_R_rebinned[iR]->SetStats(0);  // Remove stat box
        h1_R_rebinned[iR]->GetXaxis()->SetTitle("#it{p}_{T, jet}^{ch} (GeV/c)");
        h1_R_rebinned[iR]->GetYaxis()->SetTitle("d^{2}#sigma/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]");
        h1_R_rebinned[iR]->GetYaxis()->SetTitleSize(0.06);
        h1_R_rebinned[iR]->GetYaxis()->SetTitleOffset(1.15);
        h1_R_rebinned[iR]->GetXaxis()->SetLabelSize(0);
        h1_R_rebinned[iR]->GetXaxis()->SetTitleSize(0);
        
        if (!first_drawn) {
            h1_R_rebinned[iR]->Draw("PE");
            first_drawn = kTRUE;
        } else {
            h1_R_rebinned[iR]->Draw("PE SAME");
        }
        leg1->AddEntry(h1_R_rebinned[iR], TString::Format("R=%.1f, 13.6 TeV", R_values[iR]), "lp");
        
        // 13 TeV
        h2_R_rebinned[iR]->SetLineColor(colors[iR]);
        h2_R_rebinned[iR]->SetMarkerColor(colors[iR]);
        h2_R_rebinned[iR]->SetMarkerStyle(24);  // open circle
        h2_R_rebinned[iR]->SetMarkerSize(0.8);
        h2_R_rebinned[iR]->SetLineWidth(2);
        h2_R_rebinned[iR]->SetLineStyle(2);
        h2_R_rebinned[iR]->SetStats(0);  // Remove stat box
        h2_R_rebinned[iR]->Draw("PE SAME");
        leg1->AddEntry(h2_R_rebinned[iR], TString::Format("R=%.1f, 13 TeV", R_values[iR]), "lp");
    }
    
    leg1->Draw();
    
    // PYTHIA description
    TLatex* pythia_text1 = new TLatex();
    pythia_text1->SetNDC();
    pythia_text1->SetTextFont(42);
    pythia_text1->SetTextSize(0.04);
    pythia_text1->SetTextAlign(12);
    pythia_text1->DrawLatex(0.30, 0.92, "PYTHIA 8 Monash 2013");
    pythia_text1->DrawLatex(0.30, 0.87, "HardQCD:all, #tau_{0} < 10 mm/c");
    
    // Lower pad: ratio plots
    c1->cd();
    TPad* pad1_lower = new TPad("pad1_lower", "pad1_lower", 0, 0.0, 1, 0.3);
    pad1_lower->SetTopMargin(0.02);
    pad1_lower->SetBottomMargin(0.25);
    pad1_lower->SetLeftMargin(0.15);
    pad1_lower->SetRightMargin(0.05);
    pad1_lower->Draw();
    pad1_lower->cd();
    
    first_drawn = kFALSE;
    for (Int_t iR = 0; iR < nR; iR++) {
        h2_R_ratio[iR]->SetLineColor(colors[iR]);
        h2_R_ratio[iR]->SetMarkerColor(colors[iR]);
        h2_R_ratio[iR]->SetMarkerStyle(20);
        h2_R_ratio[iR]->SetMarkerSize(0.8);
        h2_R_ratio[iR]->SetLineWidth(2);
        h2_R_ratio[iR]->SetStats(0);
        h2_R_ratio[iR]->GetXaxis()->SetTitle("#it{p}_{T, jet}^{ch} (GeV/c)");
        h2_R_ratio[iR]->GetXaxis()->SetTitleSize(0.12);
        h2_R_ratio[iR]->GetXaxis()->SetLabelSize(0.1);
        h2_R_ratio[iR]->GetYaxis()->SetTitle("13 / 13.6 TeV");
        h2_R_ratio[iR]->GetYaxis()->SetTitleSize(0.12);
        h2_R_ratio[iR]->GetYaxis()->SetTitleOffset(0.6);
        h2_R_ratio[iR]->GetYaxis()->SetLabelSize(0.1);
        h2_R_ratio[iR]->GetYaxis()->SetNdivisions(505);
        h2_R_ratio[iR]->GetYaxis()->SetRangeUser(0.95, 1.01);
        
        if (!first_drawn) {
            h2_R_ratio[iR]->Draw("PE");
            first_drawn = kTRUE;
        } else {
            h2_R_ratio[iR]->Draw("PE SAME");
        }
    }
    
    // Reference line
    TLine* line1 = new TLine(h2_R_ratio[0]->GetXaxis()->GetXmin(), 1.0, 
                            h2_R_ratio[0]->GetXaxis()->GetXmax(), 1.0);
    line1->SetLineStyle(2);
    line1->SetLineColor(kBlack);
    line1->Draw("SAME");
    
    c1->cd();
    c1->SaveAs("JetPtComparison_R_dependent.pdf");
    c1->SaveAs("JetPtComparison_R_dependent.png");
    
    // Save histograms to ROOT file (use copies to avoid affecting original histograms)
    {
        TFile* outfile1 = TFile::Open("JetPtComparison_R_dependent.root", "RECREATE");
        if (!outfile1 || outfile1->IsZombie()) {
            std::cerr << "Error: Cannot create ROOT file JetPtComparison_R_dependent.root" << std::endl;
        } else {
            outfile1->cd();
            c1->Write("canvas", TObject::kOverwrite);
            for (Int_t iR = 0; iR < nR; iR++) {
                // Explicit histogram names: 13600GeV = 13.6 TeV, 13000GeV = 13 TeV
                TString h13600_name = TString::Format("hPYTHIA_13600GeV_R%02d_rebinned", iR);
                TH1D* h13600_copy = (TH1D*)h1_R_rebinned[iR]->Clone(h13600_name);
                h13600_copy->SetDirectory(0);  // Keep in memory, not attached to file
                h13600_copy->Write(h13600_name, TObject::kOverwrite);
                
                TString h13000_name = TString::Format("hPYTHIA_13000GeV_R%02d_rebinned", iR);
                TH1D* h13000_copy = (TH1D*)h2_R_rebinned[iR]->Clone(h13000_name);
                h13000_copy->SetDirectory(0);
                h13000_copy->Write(h13000_name, TObject::kOverwrite);
                
                TString hRatio_name = TString::Format("hPYTHIA_Ratio_13000GeV_over_13600GeV_R%02d", iR);
                TH1D* hRatio_copy = (TH1D*)h2_R_ratio[iR]->Clone(hRatio_name);
                hRatio_copy->SetDirectory(0);
                hRatio_copy->Write(hRatio_name, TObject::kOverwrite);
            }
            outfile1->Write();
            outfile1->Flush();
            outfile1->Close();
            delete outfile1;
        }
    }
    
    // ========== Plot 2: sigma(R=0.1) / sigma(R=X) ==========
    TCanvas* c2 = new TCanvas("c2", "Ratio to R=0.1", 800, 600);
    c2->SetLogy(0);
    c2->SetBottomMargin(0.15);  // Increase bottom margin for x-axis title
    
    TH1D* h1_R01_over_R[nR];
    TH1D* h2_R01_over_R[nR];
    
    for (Int_t iR = 0; iR < nR; iR++) {
        h1_R01_over_R[iR] = (TH1D*)h1_R_rebinned[0]->Clone(TString::Format("h1_R01_over_R%02d", iR));
        h1_R01_over_R[iR]->Divide(h1_R_rebinned[iR]);
        h1_R01_over_R[iR]->SetLineColor(colors[iR]);
        h1_R01_over_R[iR]->SetMarkerColor(colors[iR]);
        h1_R01_over_R[iR]->SetMarkerStyle(20);  // filled circle for 13.6 TeV
        h1_R01_over_R[iR]->SetMarkerSize(1.0);
        h1_R01_over_R[iR]->SetLineWidth(2);
        h1_R01_over_R[iR]->SetTitle("");
        h1_R01_over_R[iR]->SetStats(0);
        
        h2_R01_over_R[iR] = (TH1D*)h2_R_rebinned[0]->Clone(TString::Format("h2_R01_over_R%02d", iR));
        h2_R01_over_R[iR]->Divide(h2_R_rebinned[iR]);
        h2_R01_over_R[iR]->SetLineColor(colors[iR]);
        h2_R01_over_R[iR]->SetMarkerColor(colors[iR]);
        h2_R01_over_R[iR]->SetMarkerStyle(24);  // open circle for 13 TeV
        h2_R01_over_R[iR]->SetMarkerSize(1.0);
        h2_R01_over_R[iR]->SetLineWidth(2);
        h2_R01_over_R[iR]->SetLineStyle(2);
        h2_R01_over_R[iR]->SetTitle("");
        h2_R01_over_R[iR]->SetStats(0);
    }
    
    TLegend* leg2_136 = new TLegend(0.4, 0.6, 0.65, 0.9);  // 13.6 TeV legend (left by 0.2)
    leg2_136->SetBorderSize(0);
    leg2_136->SetFillStyle(0);
    leg2_136->SetTextSize(0.03);
    
    TLegend* leg2_13 = new TLegend(0.65, 0.6, 0.9, 0.9);  // 13 TeV legend (right by 0.1)
    leg2_13->SetBorderSize(0);
    leg2_13->SetFillStyle(0);
    leg2_13->SetTextSize(0.03);
    
    first_drawn = kFALSE;
    for (Int_t iR = 0; iR < nR; iR++) {
        if (iR <= 0) continue;  // Skip R=0.1 / R=0.1 (only X > 0.1)
        
        h1_R01_over_R[iR]->GetXaxis()->SetTitle("#it{p}_{T, jet}^{ch} (GeV/c)");
        h1_R01_over_R[iR]->GetXaxis()->SetTitleSize(0.05);
        h1_R01_over_R[iR]->GetXaxis()->SetTitleOffset(1.0);
        h1_R01_over_R[iR]->GetXaxis()->SetLabelSize(0.04);
        h1_R01_over_R[iR]->GetYaxis()->SetTitle("#sigma(R=0.1) / #sigma(R=X)");
        h1_R01_over_R[iR]->GetYaxis()->SetTitleSize(0.05);
        h1_R01_over_R[iR]->GetYaxis()->SetTitleOffset(0.8);
        h1_R01_over_R[iR]->GetYaxis()->SetLabelSize(0.04);
        h1_R01_over_R[iR]->GetYaxis()->SetRangeUser(0.0, 1.5);
        
        if (!first_drawn) {
            h1_R01_over_R[iR]->Draw("PE");
            first_drawn = kTRUE;
        } else {
            h1_R01_over_R[iR]->Draw("PE SAME");
        }
        leg2_136->AddEntry(h1_R01_over_R[iR], TString::Format("R=%.1f, 13.6 TeV", R_values[iR]), "lp");
        
        h2_R01_over_R[iR]->Draw("PE SAME");
        leg2_13->AddEntry(h2_R01_over_R[iR], TString::Format("R=%.1f, 13 TeV", R_values[iR]), "lp");
    }
    
    leg2_136->Draw();
    leg2_13->Draw();
    
    // PYTHIA description
    TLatex* pythia_text2 = new TLatex();
    pythia_text2->SetNDC();
    pythia_text2->SetTextFont(42);
    pythia_text2->SetTextSize(0.03);
    pythia_text2->SetTextAlign(12);
    pythia_text2->DrawLatex(0.15, 0.87, "PYTHIA 8 Monash 2013");
    pythia_text2->DrawLatex(0.15, 0.83, "HardQCD:all, #tau_{0} < 10 mm/c");
    
    c2->SaveAs("JetPtComparison_R01_over_R.pdf");
    c2->SaveAs("JetPtComparison_R01_over_R.png");
    
    // Save histograms to ROOT file (use copies to avoid affecting original histograms)
    {
        TFile* outfile2 = TFile::Open("JetPtComparison_R01_over_R.root", "RECREATE");
        if (!outfile2 || outfile2->IsZombie()) {
            std::cerr << "Error: Cannot create ROOT file JetPtComparison_R01_over_R.root" << std::endl;
        } else {
            outfile2->cd();
            c2->Write("canvas", TObject::kOverwrite);
            for (Int_t iR = 0; iR < nR; iR++) {
                if (iR <= 0) continue;
                TString h1_name = TString::Format("h1_R01_over_R%02d", iR);
                TH1D* h1_copy = (TH1D*)h1_R01_over_R[iR]->Clone(h1_name);
                h1_copy->SetDirectory(0);  // Keep in memory, not attached to file
                h1_copy->Write(h1_name, TObject::kOverwrite);
                
                TString h2_name = TString::Format("h2_R01_over_R%02d", iR);
                TH1D* h2_copy = (TH1D*)h2_R01_over_R[iR]->Clone(h2_name);
                h2_copy->SetDirectory(0);
                h2_copy->Write(h2_name, TObject::kOverwrite);
            }
            outfile2->Write();
            outfile2->Flush();
            outfile2->Close();
            delete outfile2;
        }
    }
    
    // ========== Plot 3: sigma(R=0.2) / sigma(R=X) ==========
    TCanvas* c3 = new TCanvas("c3", "Ratio to R=0.2", 800, 600);
    c3->SetLogy(0);
    c3->SetBottomMargin(0.15);  // Increase bottom margin for x-axis title
    
    TH1D* h1_R02_over_R[nR];
    TH1D* h2_R02_over_R[nR];
    
    for (Int_t iR = 0; iR < nR; iR++) {
        h1_R02_over_R[iR] = (TH1D*)h1_R_rebinned[1]->Clone(TString::Format("h1_R02_over_R%02d", iR));
        h1_R02_over_R[iR]->Divide(h1_R_rebinned[iR]);
        h1_R02_over_R[iR]->SetLineColor(colors[iR]);
        h1_R02_over_R[iR]->SetMarkerColor(colors[iR]);
        h1_R02_over_R[iR]->SetMarkerStyle(20);  // filled circle for 13.6 TeV
        h1_R02_over_R[iR]->SetMarkerSize(1.0);
        h1_R02_over_R[iR]->SetLineWidth(2);
        h1_R02_over_R[iR]->SetTitle("");
        h1_R02_over_R[iR]->SetStats(0);
        
        h2_R02_over_R[iR] = (TH1D*)h2_R_rebinned[1]->Clone(TString::Format("h2_R02_over_R%02d", iR));
        h2_R02_over_R[iR]->Divide(h2_R_rebinned[iR]);
        h2_R02_over_R[iR]->SetLineColor(colors[iR]);
        h2_R02_over_R[iR]->SetMarkerColor(colors[iR]);
        h2_R02_over_R[iR]->SetMarkerStyle(24);  // open circle for 13 TeV
        h2_R02_over_R[iR]->SetMarkerSize(1.0);
        h2_R02_over_R[iR]->SetLineWidth(2);
        h2_R02_over_R[iR]->SetLineStyle(2);
        h2_R02_over_R[iR]->SetTitle("");
        h2_R02_over_R[iR]->SetStats(0);
    }
    
    TLegend* leg3_136 = new TLegend(0.4, 0.6, 0.65, 0.9);  // 13.6 TeV legend (left by 0.2)
    leg3_136->SetBorderSize(0);
    leg3_136->SetFillStyle(0);
    leg3_136->SetTextSize(0.03);
    
    TLegend* leg3_13 = new TLegend(0.65, 0.6, 0.9, 0.9);  // 13 TeV legend (right by 0.1)
    leg3_13->SetBorderSize(0);
    leg3_13->SetFillStyle(0);
    leg3_13->SetTextSize(0.03);
    
    first_drawn = kFALSE;
    for (Int_t iR = 0; iR < nR; iR++) {
        if (iR <= 1) continue;  // Skip R=0.2 / R=0.1 and R=0.2 / R=0.2 (only X > 0.2)
        
        h1_R02_over_R[iR]->GetXaxis()->SetTitle("#it{p}_{T, jet}^{ch} (GeV/c)");
        h1_R02_over_R[iR]->GetXaxis()->SetTitleSize(0.05);
        h1_R02_over_R[iR]->GetXaxis()->SetTitleOffset(1.0);
        h1_R02_over_R[iR]->GetXaxis()->SetLabelSize(0.04);
        h1_R02_over_R[iR]->GetYaxis()->SetTitle("#sigma(R=0.2) / #sigma(R=X)");
        h1_R02_over_R[iR]->GetYaxis()->SetTitleSize(0.05);
        h1_R02_over_R[iR]->GetYaxis()->SetTitleOffset(0.8);
        h1_R02_over_R[iR]->GetYaxis()->SetLabelSize(0.04);
        h1_R02_over_R[iR]->GetYaxis()->SetRangeUser(0.0, 1.5);
        
        if (!first_drawn) {
            h1_R02_over_R[iR]->Draw("PE");
            first_drawn = kTRUE;
        } else {
            h1_R02_over_R[iR]->Draw("PE SAME");
        }
        leg3_136->AddEntry(h1_R02_over_R[iR], TString::Format("R=%.1f, 13.6 TeV", R_values[iR]), "lp");
        
        h2_R02_over_R[iR]->Draw("PE SAME");
        leg3_13->AddEntry(h2_R02_over_R[iR], TString::Format("R=%.1f, 13 TeV", R_values[iR]), "lp");
    }
    
    leg3_136->Draw();
    leg3_13->Draw();
    
    // PYTHIA description
    TLatex* pythia_text3 = new TLatex();
    pythia_text3->SetNDC();
    pythia_text3->SetTextFont(42);
    pythia_text3->SetTextSize(0.03);
    pythia_text3->SetTextAlign(12);
    pythia_text3->DrawLatex(0.15, 0.87, "PYTHIA 8 Monash 2013");
    pythia_text3->DrawLatex(0.15, 0.83, "HardQCD:all, #tau_{0} < 10 mm/c");
    
    c3->SaveAs("JetPtComparison_R02_over_R.pdf");
    c3->SaveAs("JetPtComparison_R02_over_R.png");
    
    // Save histograms to ROOT file (use copies to avoid affecting original histograms)
    {
        TFile* outfile3 = TFile::Open("JetPtComparison_R02_over_R.root", "RECREATE");
        if (!outfile3 || outfile3->IsZombie()) {
            std::cerr << "Error: Cannot create ROOT file JetPtComparison_R02_over_R.root" << std::endl;
        } else {
            outfile3->cd();
            c3->Write("canvas", TObject::kOverwrite);
            for (Int_t iR = 0; iR < nR; iR++) {
                if (iR <= 1) continue;
                TString h1_name = TString::Format("h1_R02_over_R%02d", iR);
                TH1D* h1_copy = (TH1D*)h1_R02_over_R[iR]->Clone(h1_name);
                h1_copy->SetDirectory(0);  // Keep in memory, not attached to file
                h1_copy->Write(h1_name, TObject::kOverwrite);
                
                TString h2_name = TString::Format("h2_R02_over_R%02d", iR);
                TH1D* h2_copy = (TH1D*)h2_R02_over_R[iR]->Clone(h2_name);
                h2_copy->SetDirectory(0);
                h2_copy->Write(h2_name, TObject::kOverwrite);
            }
            outfile3->Write();
            outfile3->Flush();
            outfile3->Close();
            delete outfile3;
        }
    }
    
    std::cout << "All R-dependent plots saved." << std::endl;
    
    // Cleanup
    file1->Close();
    file2->Close();
}


void DrawJetPtComparison() {
    // DrawJetPtComparison13600GeVto13000GeV();
    DrawJetPtComparisonR();
}