#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TError.h>

// Draw and compute Int(-10,10)/Int(-Inf,Inf) using a fit performed only in [-10, 10]
// Usage inside ROOT:  .x DrawZVertexExtrapolation.C("/path/to/AnalysisResults.root")
void DrawZVertexExtrapolation(const char* inputFilePath,
                              const char* histogramPath = "jet-spectra-charged/h_collisions_Zvertex",
                              double fitMin = -10.0,
                              double fitMax = 10.0,
                              const char* canvasName = "cZVtxExtrapolation") {
	// Quiet ROOT warnings for missing objects
	gErrorIgnoreLevel = kWarning;

	if (!inputFilePath || !*inputFilePath) {
		::Error("DrawZVertexExtrapolation", "inputFilePath is empty.");
		return;
	}

	// Open file
	TFile* inputFile = TFile::Open(inputFilePath, "READ");
	if (!inputFile || inputFile->IsZombie()) {
		::Error("DrawZVertexExtrapolation", "Failed to open file: %s", inputFilePath);
		return;
	}

	// Fetch histogram
	TH1* hist = dynamic_cast<TH1*>(inputFile->Get(histogramPath));
	if (!hist) {
		::Error("DrawZVertexExtrapolation", "Histogram not found at path: %s", histogramPath);
		inputFile->Close();
		return;
	}

	// Prepare canvas
	TCanvas* canvas = new TCanvas(canvasName, canvasName, 900, 700);
	canvas->SetMargin(0.12, 0.04, 0.12, 0.08);

	// Style histogram
	hist->SetLineColor(kBlack);
	hist->SetMarkerStyle(20);
	hist->SetMarkerSize(1.0);
	hist->SetMarkerColor(kBlack);
	if (hist->GetXaxis()) hist->GetXaxis()->SetTitle("z_{vtx} (cm)");
	if (hist->GetYaxis()) hist->GetYaxis()->SetTitle("Counts");

	// Draw histogram first
	hist->Draw("E");

	// Get raw histogram mean to fix during fit
	double histMean = hist->GetMean();

	// Fit only in [-10, 10] with mean fixed to raw histogram mean
	TF1* fitFunc = new TF1("fitFunc", "gaus", fitMin, fitMax);
	// Set initial mean to histogram mean and fix it
	fitFunc->SetParameter(1, histMean);  // Set mean parameter
	fitFunc->FixParameter(1, histMean);  // Fix mean parameter
	// Fit with mean fixed (only amplitude and sigma will be fitted)
	hist->Fit(fitFunc, "RQ0");
	// Expand function drawing range to full axis for visualization
	double axisMin = hist->GetXaxis() ? hist->GetXaxis()->GetXmin() : -50.0;
	double axisMax = hist->GetXaxis() ? hist->GetXaxis()->GetXmax() : 50.0;
	fitFunc->SetRange(axisMin, axisMax);
	fitFunc->SetLineColor(kRed+1);
	fitFunc->SetLineWidth(2);
	fitFunc->Draw("SAME");

	// Compute integrals using the fitted function
	// Numerator: Integral in [-10, 10]
	double numerator = fitFunc->Integral(fitMin, fitMax);

	// Denominator: Integral in (-Inf, +Inf) — approximate with wide numeric bounds
	// For a Gaussian, [-1e6, 1e6] is practically equivalent to (-Inf, +Inf)
	double denomLower = -1.0e6;
	double denomUpper =  1.0e6;
	double denominator = fitFunc->Integral(denomLower, denomUpper);

	double ratio = (denominator > 0.0) ? (numerator / denominator) : 0.0;

	// Get fitted Gaussian parameters
	double constant = fitFunc->GetParameter(0);
	double mean = fitFunc->GetParameter(1);
	double sigma = fitFunc->GetParameter(2);

	// Annotate results on canvas
	TLatex latex;
	latex.SetNDC(true);
	latex.SetTextSize(0.03);
	// latex.DrawLatex(0.15, 0.85, Form("Histogram: %s", histogramPath));
	// latex.DrawLatex(0.15, 0.80, Form("Fit range: [%.1f, %.1f] cm (gaus)", fitMin, fitMax));
	
	// Display Gaussian function
	latex.DrawLatex(0.15, 0.92, Form("f(x) = %.3e #times exp(-0.5#times((x-%.3f)/%.3f)^{2})", constant, mean, sigma));
	latex.DrawLatex(0.15, 0.87, Form("Mean = %.3f cm, #sigma = %.3f cm", mean, sigma));
	
	latex.DrawLatex(0.15, 0.78, Form("#int_{%.0f}^{%.0f} f(x) dx = %.3e", fitMin, fitMax, numerator));
	latex.DrawLatex(0.15, 0.72, Form("#int_{-#infty}^{+#infty} f(x) dx #approx %.3e", denominator));
	latex.SetTextSize(0.035);
	latex.SetTextColor(kBlue+2);
	latex.DrawLatex(0.15, 0.62, Form("#int_{-10}^{10}f/#int_{-#infty}^{+#infty}f = %.5f", ratio));

	// Optional legend
	TLegend* legend = new TLegend(0.75, 0.55, 0.95, 0.70);
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);
	legend->SetTextSize(0.03);
	legend->AddEntry(hist, "Data", "lep");
	legend->AddEntry(fitFunc, "Fit (gaus) [#pm10 cm]", "l");
	legend->Draw();

	canvas->Modified();
	canvas->Update();

	// Keep file open while canvas is alive so that histogram stays valid in interactive sessions.
	// If running in batch and you want to save, uncomment below lines and provide an output path.
	canvas->SaveAs("plots/ZVertexExtrapolation.pdf");

	// Note: Not closing inputFile here on purpose to avoid deleting objects if user keeps the canvas.
	// If you prefer explicit ownership, you can Clone() the histogram and then close the file safely.
}

void DrawZVertexEfficiency() {
    DrawZVertexExtrapolation("../../../jets/AnalysisResults/498133_AnalysisResults.root");
}