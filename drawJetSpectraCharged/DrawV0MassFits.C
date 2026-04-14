// Draw V0 invariant mass per pT bin and fit Gaussian core to extract mean and sigma
// author: Joonsuk Bae (request), implemented by assistant

#if defined(__CLING__) || defined(__CINT__) || defined(__ROOTCLING__)

#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TString.h>
#include <THnSparse.h>
#include <TArrayD.h>
#include <vector>
#include <string>

// Fit models from Gijs' slides
static const char *kPol1ExpGausExp = "[0]*TMath::Gaus(x,[1],[2]) * (x > ([1]-[2]*[3])) * (x < ([1]+[2]*[4])) + [0]*TMath::Exp((x - [1] + [2]*[3]/2.)/([2]/[3])) * (x <= [1]-[2]*[3]) + [0]*TMath::Exp(-(x - [1] - [2]*[4]/2.)/([2]/[4])) * (x >= [1]+[2]*[4]) + [5]+[6]*x";
static const char *kPol1GausGaus   = "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[1],[4]) + [5]+[6]*x";
static const char *kPol1GausGausExp= "[0]*TMath::Gaus(x,[1],[2]) * (x < ([1]+[2]*[3])) + [0]*TMath::Exp(-(x - [1] - [2]*[3]/2.)/([2]/[3])) * (x >= [1]+[2]*[3]) + [4]*TMath::Gaus(x,[1],[5]) + [6]+[7]*x";

struct FitResult { double pt; double ptLow; double ptHigh; double mean; double sigma; double meanErr; double sigmaErr; double meanConv; double sigmaConv; double chi2ndf; int status; };

static std::vector<double> MakePtBins() {
  std::vector<double> edges;
  for (int i=0;i<=10;i++) edges.push_back(i);        // 0-10: 1 GeV intervals
  for (int p=12;p<=20;p+=2) edges.push_back(p);      // 10-20: 2 GeV intervals
  for (int p=25;p<=30;p+=5) edges.push_back(p);      // 20-30: 5 GeV intervals
  for (int p=40;p<=50;p+=10) edges.push_back(p);     // 30-50: 10 GeV intervals
  return edges;
}

static TF1* BuildFunction(const char *modelName, double xMin, double xMax) {
  TString m(modelName);
  if (m.EqualTo("pol1ExpGausExp")) {
    auto f = new TF1("pol1ExpGausExp", kPol1ExpGausExp, xMin, xMax);
    return f;
  } else if (m.EqualTo("pol1GausGaus")) {
    auto f = new TF1("pol1GausGaus", kPol1GausGaus, xMin, xMax);
    return f;
  } else if (m.EqualTo("pol1GausGausExp")) {
    auto f = new TF1("pol1GausGausExp", kPol1GausGausExp, xMin, xMax);
    return f;
  }
  // default
  auto f = new TF1("pol1ExpGausExp", kPol1ExpGausExp, xMin, xMax);
  return f;
}

static void SeedParameters(TF1 *f, const char *modelName, TH1 *h, const char *v0Type = "K0S") {
  double amp = h->GetMaximum();
  double mass0, sigma0;
  TString v0(v0Type);
  if (v0.EqualTo("K0S")) {
    mass0 = 0.4976; // K0S mass (GeV/c^2)
    sigma0 = 0.005;
  } else if (v0.EqualTo("Lambda")) {
    mass0 = 1.1157; // Lambda mass (GeV/c^2)
    sigma0 = 0.002;
  } else if (v0.EqualTo("AntiLambda")) {
    mass0 = 1.1157; // AntiLambda mass (GeV/c^2)
    sigma0 = 0.002;
  } else {
    mass0 = 0.4976; // default to K0S
    sigma0 = 0.005;
  }
  TString m(modelName);
  if (m.EqualTo("pol1ExpGausExp")) {
    f->SetParameter(0, amp);
    f->SetParameter(1, mass0);
    f->SetParameter(2, sigma0);
    f->SetParameter(3, 1.0);
    f->SetParameter(4, 1.0);
    f->SetParameter(5, h->GetBinContent(1));
    f->SetParameter(6, 0.0);
  } else if (m.EqualTo("pol1GausGaus")) {
    f->SetParameter(0, 0.7*amp);
    f->SetParameter(1, mass0);
    f->SetParameter(2, sigma0);
    f->SetParameter(3, 0.3*amp);
    f->SetParameter(4, 1.6*sigma0);
    f->SetParameter(5, h->GetBinContent(1));
    f->SetParameter(6, 0.0);
  } else if (m.EqualTo("pol1GausGausExp")) {
    f->SetParameter(0, 0.7*amp);
    f->SetParameter(1, mass0);
    f->SetParameter(2, sigma0);
    f->SetParameter(3, 1.0);
    f->SetParameter(4, 0.3*amp);
    f->SetParameter(5, 1.6*sigma0);
    f->SetParameter(6, h->GetBinContent(1));
    f->SetParameter(7, 0.0);
  }
}

static FitResult FitOne(TH1 *hMass, double ptLow, double ptHigh, const char *modelName, const char *outDir, int idx, const char *v0Type = "K0S") {
  FitResult r{}; r.ptLow=ptLow; r.ptHigh=ptHigh; r.pt=0.5*(ptLow+ptHigh);
  
  // Set mass window based on V0 type
  double xMin, xMax;
  TString v0(v0Type);
  if (v0.EqualTo("K0S")) {
    xMin = 0.43; xMax = 0.57; // K0S mass window
  } else if (v0.EqualTo("Lambda") || v0.EqualTo("AntiLambda")) {
    xMin = 1.05; xMax = 1.18; // Lambda mass window
  } else {
    xMin = 0.43; xMax = 0.57; // default to K0S
  }
  
  auto f = BuildFunction(modelName, xMin, xMax);
  SeedParameters(f, modelName, hMass, v0Type);
  hMass->SetDirectory(nullptr);

  TCanvas *c = new TCanvas(Form("c_v0_%d", idx), Form("V0 mass %.1f-%.1f GeV/c", ptLow, ptHigh), 800, 600);
  
  // Set axis labels based on V0 type
  if (v0.EqualTo("K0S")) {
    hMass->GetXaxis()->SetTitle("m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2})");
  } else if (v0.EqualTo("Lambda")) {
    hMass->GetXaxis()->SetTitle("m_{p#pi^{-}} (GeV/#it{c}^{2})");
  } else if (v0.EqualTo("AntiLambda")) {
    hMass->GetXaxis()->SetTitle("m_{#bar{p}#pi^{+}} (GeV/#it{c}^{2})");
  } else {
    hMass->GetXaxis()->SetTitle("m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2})");
  }
  hMass->GetYaxis()->SetTitle("Counts");
  hMass->Draw("E");

  auto fitRes = hMass->Fit(f, "SR");
  int status = (fitRes ? fitRes->Status() : -1);
  r.status = status;

  // Extract mean/sigma from the Gaussian core: parameters [1] and [2]
  r.mean = f->GetParameter(1);
  r.sigma = f->GetParameter(2);
  r.meanErr = f->GetParError(1);
  r.sigmaErr = f->GetParError(2);
  r.chi2ndf = (f->GetNDF() > 0 ? f->GetChisquare()/f->GetNDF() : -1);

  f->SetLineColor(kRed+1);
  f->SetLineWidth(2);
  f->Draw("same");

  // Overlay Gaussian component only (signal core)
  TF1 *fGaus = new TF1(Form("gaus_core_%d", idx), "[0]*TMath::Gaus(x,[1],[2])", xMin, xMax);
  fGaus->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
  fGaus->SetLineColor(kBlue+1);
  fGaus->SetLineStyle(2);
  fGaus->SetLineWidth(2);
  fGaus->Draw("same");

  // Legend on the left
  TLegend *leg = new TLegend(0.14, 0.75, 0.47, 0.90);
  leg->SetBorderSize(0);
  leg->AddEntry(f, "Full fit (signal+background)", "l");
  leg->AddEntry(fGaus, "Core Gaussian", "l");
  leg->Draw();

  // Compute effective mean/sigma from full convolution shape (signal = total minus pol1 background)
  // Identify background parameters index depending on model
  double b0 = 0.0, b1 = 0.0;
  TString m(modelName);
  if (m.EqualTo("pol1ExpGausExp") || m.EqualTo("pol1GausGaus")) {
    b0 = f->GetParameter(5);
    b1 = f->GetParameter(6);
  } else if (m.EqualTo("pol1GausGausExp")) {
    b0 = f->GetParameter(6);
    b1 = f->GetParameter(7);
  }
  double mu = r.mean;
  double sig = TMath::Abs(r.sigma);
  double a = TMath::Max(xMin, mu - 5.0*sig);
  double b = TMath::Min(xMax, mu + 5.0*sig);
  const int N = 2000;
  double sumW = 0.0, sumWX = 0.0, sumWX2 = 0.0;
  for (int i = 0; i < N; ++i) {
    double x = a + (b - a) * (i + 0.5) / N;
    double y = f->Eval(x) - (b0 + b1 * x);
    if (y < 0) y = 0; // avoid negative after background subtraction
    sumW += y;
    sumWX += y * x;
    sumWX2 += y * x * x;
  }
  if (sumW > 0) {
    double meanConv = sumWX / sumW;
    double varConv = sumWX2 / sumW - meanConv * meanConv;
    r.meanConv = meanConv;
    r.sigmaConv = (varConv > 0 ? TMath::Sqrt(varConv) : 0.0);
  } else {
    r.meanConv = 0.0;
    r.sigmaConv = 0.0;
  }

  TLatex lat; lat.SetNDC(true); lat.SetTextSize(0.04);
  lat.DrawLatex(0.14, 0.70, Form("%.1f < p_{T} < %.1f GeV/c", ptLow, ptHigh));
  lat.DrawLatex(0.14, 0.65, Form("Gauss: #mu = %.5f, #sigma = %.5f", r.mean, r.sigma));
  lat.DrawLatex(0.14, 0.60, Form("#chi^{2}/ndf = %.2f", r.chi2ndf));

  gSystem->MakeDirectory(outDir);
  // Use idx >= 1000 to identify MC plots and add "MC" suffix to filename
  const char* suffix = (idx >= 1000) ? "_MC" : "";
  c->SaveAs(Form("%s/V0MassFit_%.1f_%.1f%s.pdf", outDir, ptLow, ptHigh, suffix));

  delete c;
  return r;
}

static TObject* FindObjectRecursive(TDirectory* dir, const char* nameSubstr) {
  if (!dir) return nullptr;
  TIter next(dir->GetListOfKeys());
  while (TObject* keyObj = next()) {
    auto key = dynamic_cast<TKey*>(keyObj);
    if (!key) continue;
    TObject* obj = key->ReadObj();
    if (!obj) continue;
    // Match by substring on full path or name
    TString objName = obj->GetName();
    TString objClass = obj->ClassName();
    if (objName.Contains(nameSubstr) && (objClass.Contains("TH2") || objClass.Contains("TH2F") || objClass.Contains("TH2D"))) {
      return obj;
    }
    auto subdir = dynamic_cast<TDirectory*>(obj);
    if (subdir) {
      TObject* found = FindObjectRecursive(subdir, nameSubstr);
      if (found) return found;
    }
  }
  return nullptr;
}

static void ListCandidatesRecursive(TDirectory* dir, const char* nameFilter = "V0") {
  if (!dir) return;
  TString dirPath = dir->GetPath();
  TIter next(dir->GetListOfKeys());
  while (TObject* keyObj = next()) {
    auto key = dynamic_cast<TKey*>(keyObj);
    if (!key) continue;
    TObject* obj = key->ReadObj();
    if (!obj) continue;
    TString name = obj->GetName();
    TString cls = obj->ClassName();
    if (name.Contains(nameFilter)) {
      printf("[Hint] %s/%s (%s)\n", dirPath.Data(), name.Data(), cls.Data());
    }
    auto subdir = dynamic_cast<TDirectory*>(obj);
    if (subdir) ListCandidatesRecursive(subdir, nameFilter);
  }
}

// Main entry (internal function)
void DrawV0MassFitsInternal(const char *RunNumber,
                           const char *histPath,
                           const char *outDir,
                           const char *fitModel,
                           const char *v0Type,
                           const char *mcRunNumber,
                           const char *mcHistPath) {
  if (strlen(RunNumber)==0) {
    printf("[Warn] RunNumber not provided. Please pass the run number (e.g. 426828)\n");
    return;
  }

  // Load data file
  auto f = TFile::Open(Form("../../../jets/AnalysisResults/%s_AnalysisResults.root", RunNumber), "read");
  if (!f || f->IsZombie()) {
    printf("[Error] Cannot open data file: %s\n", Form("../../../jets/AnalysisResults/%s_AnalysisResults.root", RunNumber));
    return;
  }

  // Load MC file if provided
  TFile *fMC = nullptr;
  if (strlen(mcRunNumber) > 0) {
    fMC = TFile::Open(Form("../../../jets/AnalysisResults/%s_AnalysisResults.root", mcRunNumber), "read");
    if (!fMC || fMC->IsZombie()) {
      printf("[Warning] Cannot open MC file: %s\n", Form("../../../jets/AnalysisResults/%s_AnalysisResults.root", mcRunNumber));
      fMC = nullptr;
    } else {
      printf("[Info] MC file loaded: %s\n", mcRunNumber);
    }
  }

  // Load data histogram
  std::vector<TString> candidates;
  candidates.push_back(histPath);
  candidates.push_back("jet-fragmentation_id28293/data/V0/V0PtMass");

  TH2 *h2raw = nullptr;
  THnSparseT<TArrayD> *hsparse = nullptr;
  for (const auto &cand : candidates) {
    TObject *obj = f->Get(cand);
    if (!obj) continue;
    h2raw = dynamic_cast<TH2*>(obj);
    if (!h2raw) h2raw = dynamic_cast<TH2F*>(obj);
    if (!h2raw) h2raw = dynamic_cast<TH2D*>(obj);
    if (!h2raw) hsparse = dynamic_cast<THnSparseT<TArrayD>*>(obj);
    if (h2raw) {
      printf("[Info] Loaded data histogram: %s\n", cand.Data());
      break;
    } else if (hsparse) {
      printf("[Info] Loaded data THnSparse: %s\n", cand.Data());
      break;
    }
  }

  // Load MC histogram if MC file is available
  TH2 *h2rawMC = nullptr;
  THnSparseT<TArrayD> *hsparseMC = nullptr;
  if (fMC) {
    std::vector<TString> mcCandidates;
    mcCandidates.push_back(mcHistPath);
    mcCandidates.push_back("jet-fragmentation/mcd/V0/V0PtMass");
    
    for (const auto &cand : mcCandidates) {
      TObject *obj = fMC->Get(cand);
      if (!obj) continue;
      h2rawMC = dynamic_cast<TH2*>(obj);
      if (!h2rawMC) h2rawMC = dynamic_cast<TH2F*>(obj);
      if (!h2rawMC) h2rawMC = dynamic_cast<TH2D*>(obj);
      if (!h2rawMC) hsparseMC = dynamic_cast<THnSparseT<TArrayD>*>(obj);
      if (h2rawMC) {
        printf("[Info] Loaded MC histogram: %s\n", cand.Data());
        break;
      } else if (hsparseMC) {
        printf("[Info] Loaded MC THnSparse: %s\n", cand.Data());
        break;
      }
    }
  }
  // Fallback: recursive search by name substring (only if nothing found yet)
  if (!h2raw && !hsparse) {
    TObject* found = FindObjectRecursive(f, "V0PtMass");
    h2raw = dynamic_cast<TH2*>(found);
    if (!h2raw) h2raw = dynamic_cast<TH2F*>(found);
    if (!h2raw) h2raw = dynamic_cast<TH2D*>(found);
    if (!h2raw && found) hsparse = dynamic_cast<THnSparseT<TArrayD>*>(found);
    if (h2raw) {
      printf("[Info] Found histogram by recursive search: %s (%s)\n", h2raw->GetName(), h2raw->ClassName());
    } else if (hsparse) {
      printf("[Info] Found THnSparse by recursive search: %s (%s)\n", hsparse->GetName(), hsparse->ClassName());
    }
  }
  if (!h2raw && !hsparse) {
    printf("[Error] Data histogram not found. Tried paths including: %s\n", histPath);
    printf("[Error] File: %s\n", Form("../../../jets/AnalysisResults/%s_AnalysisResults.root", RunNumber));
    printf("[Info] Listing possible candidates containing 'V0' to help locate the histogram...\n");
    ListCandidatesRecursive(f, "V0");
    f->Close();
    if (fMC) fMC->Close();
    return;
  }

  // Check MC histogram if MC file was provided
  if (fMC && !h2rawMC && !hsparseMC) {
    printf("[Warning] MC histogram not found. Tried paths including: %s\n", mcHistPath);
    printf("[Warning] MC File: %s\n", Form("../../../jets/AnalysisResults/%s_AnalysisResults.root", mcRunNumber));
    printf("[Info] Listing possible MC candidates containing 'V0' to help locate the histogram...\n");
    ListCandidatesRecursive(fMC, "V0");
  }

  std::vector<double> ptEdges = MakePtBins();
  const int nBins = (int)ptEdges.size()-1;

  std::vector<FitResult> results;
  results.reserve(nBins);

  // Summary graphs for data
  TGraphErrors *gSigma = new TGraphErrors(nBins);
  TGraphErrors *gMean  = new TGraphErrors(nBins);
  TGraphErrors *gSigmaConv = new TGraphErrors(nBins);
  TGraphErrors *gRelSigma = new TGraphErrors(nBins);
  TGraphErrors *gRelSigmaConv = new TGraphErrors(nBins);
  gSigma->SetName("gSigmaVsPt");
  gMean->SetName("gMeanVsPt");
  gSigmaConv->SetName("gSigmaConvVsPt");
  gRelSigma->SetName("gRelSigmaVsPt");
  gRelSigmaConv->SetName("gRelSigmaConvVsPt");

  // Summary graphs for MC (if available)
  TGraphErrors *gSigmaMC = nullptr;
  TGraphErrors *gMeanMC  = nullptr;
  TGraphErrors *gSigmaConvMC = nullptr;
  TGraphErrors *gRelSigmaMC = nullptr;
  TGraphErrors *gRelSigmaConvMC = nullptr;
  
  if (h2rawMC || hsparseMC) {
    gSigmaMC = new TGraphErrors(nBins);
    gMeanMC = new TGraphErrors(nBins);
    gSigmaConvMC = new TGraphErrors(nBins);
    gRelSigmaMC = new TGraphErrors(nBins);
    gRelSigmaConvMC = new TGraphErrors(nBins);
    gSigmaMC->SetName("gSigmaVsPtMC");
    gMeanMC->SetName("gMeanVsPtMC");
    gSigmaConvMC->SetName("gSigmaConvVsPtMC");
    gRelSigmaMC->SetName("gRelSigmaVsPtMC");
    gRelSigmaConvMC->SetName("gRelSigmaConvVsPtMC");
  }

  for (int i=0;i<nBins;i++) {
    double lo = ptEdges[i];
    double hi = ptEdges[i+1];

    // Process data (no pT limit)
    TH1 *hProj = nullptr;
    if (h2raw) {
      int binLo = h2raw->GetXaxis()->FindBin(lo+1e-6);
      int binHi = h2raw->GetXaxis()->FindBin(hi-1e-6);
      if (binHi < binLo) binHi = binLo;
      hProj = h2raw->ProjectionY(Form("hMass_pt_%.1f_%.1f", lo, hi), binLo, binHi);
    } else if (hsparse) {
      // Axis 0: pT, Axis 1: mass (per Gijs)
      auto axPt = hsparse->GetAxis(0);
      auto axM  = hsparse->GetAxis(1);
      int binLo = axPt->FindBin(lo+1e-6);
      int binHi = axPt->FindBin(hi-1e-6);
      if (binHi < binLo) binHi = binLo;
      axPt->SetRange(binLo, binHi);
      // Project mass axis
      hProj = hsparse->Projection(1);
      hProj->SetName(Form("hMass_pt_%.1f_%.1f", lo, hi));
      // Clear range for next loop
      axPt->SetRange(0, 0);
    }

    FitResult r{};
    r.ptLow=lo; r.ptHigh=hi; r.pt=0.5*(lo+hi); r.status=-2;
    
    // Draw individual mass peak plots for QA up to 50 GeV/c
    // For QA purposes, create plots even with fewer entries
    if (hProj && hi <= 50.0) {
      if (hProj->GetEntries() >= 1) {
        r = FitOne(hProj, lo, hi, fitModel, outDir, i, v0Type);
      }
    }
    results.push_back(r);
    
    // Fill data summary graphs only for pT <= 20 GeV/c (for resolution plots)
    if (hi <= 20.0) {
      gSigma->SetPoint(i, r.pt, r.sigma);
      gSigma->SetPointError(i, 0.5*(hi-lo), r.sigmaErr);
      gMean->SetPoint(i, r.pt, r.mean);
      gMean->SetPointError(i, 0.5*(hi-lo), r.meanErr);
      gSigmaConv->SetPoint(i, r.pt, r.sigmaConv);
      gSigmaConv->SetPointError(i, 0.5*(hi-lo), 0.0);
      // relative resolution sigma/mean
      double relG = (r.mean != 0.0) ? r.sigma / r.mean : 0.0;
      double relGErr = 0.0;
      if (r.mean != 0.0) {
        double dfdSigma = 1.0 / r.mean;                  // ∂(σ/μ)/∂σ
        double dfdMu    = -r.sigma / (r.mean * r.mean);  // ∂(σ/μ)/∂μ
        double var = dfdSigma*dfdSigma*(r.sigmaErr*r.sigmaErr)
                   + dfdMu*dfdMu*(r.meanErr*r.meanErr); // covariance term omitted
        relGErr = (var > 0.0 ? TMath::Sqrt(var) : 0.0);
      }
      double relC = (r.meanConv != 0.0) ? r.sigmaConv / r.meanConv : 0.0;
      gRelSigma->SetPoint(i, r.pt, relG);
      gRelSigma->SetPointError(i, 0.5*(hi-lo), relGErr);
      gRelSigmaConv->SetPoint(i, r.pt, relC);
      gRelSigmaConv->SetPointError(i, 0.5*(hi-lo), 0.0);
    }

    // Process MC if available
    // For individual plots: process up to 50 GeV/c (30-40, 40-50 bins)
    // For summary graphs: only fill points up to 20 GeV/c
    if (gSigmaMC && hi <= 50.0) {
      TH1 *hProjMC = nullptr;
      if (h2rawMC) {
        int binLo = h2rawMC->GetXaxis()->FindBin(lo+1e-6);
        int binHi = h2rawMC->GetXaxis()->FindBin(hi-1e-6);
        if (binHi < binLo) binHi = binLo;
        hProjMC = h2rawMC->ProjectionY(Form("hMassMC_pt_%.1f_%.1f", lo, hi), binLo, binHi);
      } else if (hsparseMC) {
        // Axis 0: pT, Axis 1: mass (per Gijs)
        auto axPt = hsparseMC->GetAxis(0);
        auto axM  = hsparseMC->GetAxis(1);
        int binLo = axPt->FindBin(lo+1e-6);
        int binHi = axPt->FindBin(hi-1e-6);
        if (binHi < binLo) binHi = binLo;
        axPt->SetRange(binLo, binHi);
        // Project mass axis
        hProjMC = hsparseMC->Projection(1);
        hProjMC->SetName(Form("hMassMC_pt_%.1f_%.1f", lo, hi));
        // Clear range for next loop
        axPt->SetRange(0, 0);
      }

      FitResult rMC{};
      rMC.ptLow=lo; rMC.ptHigh=hi; rMC.pt=0.5*(lo+hi); rMC.status=-2;
      
      // Draw individual mass peak plots for MC up to 50 GeV/c
      // For QA purposes, create plots even with fewer entries
      if (hProjMC && hProjMC->GetEntries() >= 1) {
        rMC = FitOne(hProjMC, lo, hi, fitModel, outDir, i+1000, v0Type); // Use different index for MC
      }
      
      // Fill MC summary graphs only for pT <= 20 GeV/c
      if (hi <= 20.0) {
        gSigmaMC->SetPoint(i, rMC.pt, rMC.sigma);
        gSigmaMC->SetPointError(i, 0.5*(hi-lo), rMC.sigmaErr);
        gMeanMC->SetPoint(i, rMC.pt, rMC.mean);
        gMeanMC->SetPointError(i, 0.5*(hi-lo), rMC.meanErr);
        gSigmaConvMC->SetPoint(i, rMC.pt, rMC.sigmaConv);
        gSigmaConvMC->SetPointError(i, 0.5*(hi-lo), 0.0);
        // relative resolution sigma/mean
        double relGMC = (rMC.mean != 0.0) ? rMC.sigma / rMC.mean : 0.0;
        double relGErrMC = 0.0;
        if (rMC.mean != 0.0) {
          double dfdSigma = 1.0 / rMC.mean;
          double dfdMu    = -rMC.sigma / (rMC.mean * rMC.mean);
          double var = dfdSigma*dfdSigma*(rMC.sigmaErr*rMC.sigmaErr)
                     + dfdMu*dfdMu*(rMC.meanErr*rMC.meanErr);
          relGErrMC = (var > 0.0 ? TMath::Sqrt(var) : 0.0);
        }
        double relCMC = (rMC.meanConv != 0.0) ? rMC.sigmaConv / rMC.meanConv : 0.0;
        gRelSigmaMC->SetPoint(i, rMC.pt, relGMC);
        gRelSigmaMC->SetPointError(i, 0.5*(hi-lo), relGErrMC);
        gRelSigmaConvMC->SetPoint(i, rMC.pt, relCMC);
        gRelSigmaConvMC->SetPointError(i, 0.5*(hi-lo), 0.0);
      }
    }
  }

  // Draw summary graphs
  gSystem->MakeDirectory(outDir);

  {
    // Resolution vs pT: overlay relative resolutions (sigma/mean)
    TCanvas *c1 = new TCanvas("c_sigma", "Resolution vs pT", 800, 600);
    TString v0(v0Type);
    TString title;
    if (v0.EqualTo("K0S")) {
      title = "K^{0}_{S} mass resolution; p_{T} (GeV/c); #sigma / #mu";
    } else if (v0.EqualTo("Lambda")) {
      title = "#Lambda mass resolution; p_{T} (GeV/c); #sigma / #mu";
    } else if (v0.EqualTo("AntiLambda")) {
      title = "#bar{#Lambda} mass resolution; p_{T} (GeV/c); #sigma / #mu";
    } else {
      title = "V0 mass resolution; p_{T} (GeV/c); #sigma / #mu";
    }
    gRelSigma->SetTitle(title.Data());
    gRelSigma->SetMarkerStyle(21);  // Data Core: 빨간 찬 박스
    gRelSigma->SetMarkerColor(kRed);
    gRelSigma->SetLineColor(kRed);
    gRelSigma->Draw("AP");
    
    // Set y-axis range to 0-0.1
    gRelSigma->GetHistogram()->GetYaxis()->SetRangeUser(0.0, 0.1);

    gRelSigmaConv->SetMarkerStyle(25);  // Data RMS: 빨간 빈 박스
    gRelSigmaConv->SetMarkerColor(kRed);
    gRelSigmaConv->SetLineColor(kRed);
    gRelSigmaConv->Draw("P same");

    // Add MC if available
    if (gRelSigmaMC) {
      gRelSigmaMC->SetMarkerStyle(20);  // MC Core: 검정 찬 박스
      gRelSigmaMC->SetMarkerColor(kBlack);
      gRelSigmaMC->SetLineColor(kBlack);
      gRelSigmaMC->Draw("P same");

      gRelSigmaConvMC->SetMarkerStyle(24);  // MC RMS: 검정 빈 박스
      gRelSigmaConvMC->SetMarkerColor(kBlack);
      gRelSigmaConvMC->SetLineColor(kBlack);
      gRelSigmaConvMC->Draw("P same");
    }

    auto leg = new TLegend(0.6, 0.67, 0.88, 0.87);
    leg->SetBorderSize(0);
    leg->AddEntry(gRelSigma, "Data Core #sigma/#mu", "lpe");
    leg->AddEntry(gRelSigmaConv, "Data RMS #sigma/#mu", "lpe");
    if (gRelSigmaMC) {
      leg->AddEntry(gRelSigmaMC, "MC Core #sigma/#mu", "lpe");
      leg->AddEntry(gRelSigmaConvMC, "MC RMS #sigma/#mu", "lpe");
    }
    leg->Draw();
    
    // Add convolution formula annotation
    TLatex convFormula;
    convFormula.SetNDC(true);
    convFormula.SetTextSize(0.03);
    convFormula.SetTextColor(kGray+2);
    convFormula.DrawLatex(0.12, 0.85, "RMS = #sqrt{<x^{2}> - <x>^{2}} from signal-only distribution");

    c1->SaveAs(Form("%s/V0MassFitResolutionVsPt.pdf", outDir));
    delete c1;
  }

  // Ratio plots: Data/MC for Core (sigma/mu) and RMS (sigmaConv/muConv)
  if (gRelSigmaMC) {
    // Core ratio (Data/MC)
    // Match points by x-value (pT) to ensure correct ratio calculation
    int nRatioPoints = 0;
    TGraphErrors *gRatioCore = new TGraphErrors();
    gRatioCore->SetName("gRelSigmaRatioCoreDataOverMC");
    
    for (int ipD = 0; ipD < gRelSigma->GetN(); ++ipD) {
      double xD, yD; gRelSigma->GetPoint(ipD, xD, yD);
      
      // Find matching MC point with same x-value (pT)
      int ipM = -1;
      for (int j = 0; j < gRelSigmaMC->GetN(); ++j) {
        double xM, yM; gRelSigmaMC->GetPoint(j, xM, yM);
        if (TMath::Abs(xD - xM) < 0.01) { // Allow small tolerance for floating point comparison
          ipM = j;
          break;
        }
      }
      
      if (ipM >= 0) {
        double xM, yM; gRelSigmaMC->GetPoint(ipM, xM, yM);
        double eDy = gRelSigma->GetErrorY(ipD);
        double eMy = gRelSigmaMC->GetErrorY(ipM);
        double ratio = (yM != 0.0 ? yD / yM : 0.0);
        double eratio = 0.0;
        if (yD > 0.0 && yM > 0.0) {
          double relD = (eDy / yD);
          double relM = (eMy / yM);
          double var = relD*relD + relM*relM; // covariance omitted
          eratio = ratio * (var > 0.0 ? TMath::Sqrt(var) : 0.0);
        }
        
        // Calculate x error from pt bin width
        double xErr = 0.0;
        for (int k = 0; k < (int)ptEdges.size()-1; ++k) {
          double lo = ptEdges[k];
          double hi = ptEdges[k+1];
          if (TMath::Abs(xD - 0.5*(lo+hi)) < 0.01) {
            xErr = 0.5*(hi-lo);
            break;
          }
        }
        
        gRatioCore->SetPoint(nRatioPoints, xD, ratio);
        gRatioCore->SetPointError(nRatioPoints, xErr, eratio);
        nRatioPoints++;
      } else {
        printf("[Warning] Core ratio: No matching MC point found for Data pT = %.2f GeV/c\n", xD);
      }
    }
    TCanvas *cRC = new TCanvas("c_ratio_core", "Data/MC Ratio (Core)", 800, 600);
    gRatioCore->SetTitle("Data/MC ratio (Core); p_{T} (GeV/c); (Data/MC) #sigma/#mu");
    gRatioCore->SetMarkerStyle(21);
    gRatioCore->SetMarkerColor(kBlue+2);
    gRatioCore->SetLineColor(kBlue+2);
    gRatioCore->GetYaxis()->SetRangeUser(0.0, 2.0);
    gRatioCore->Draw("AP");
    cRC->SaveAs(Form("%s/V0MassFitResolutionVsPt_RatioCore.pdf", outDir));
    delete cRC;

    // RMS ratio (Data/MC)
    // Match points by x-value (pT) to ensure correct ratio calculation
    int nRatioPointsRMS = 0;
    TGraphErrors *gRatioRMS = new TGraphErrors();
    gRatioRMS->SetName("gRelSigmaRatioRMSDataOverMC");
    
    for (int ipD = 0; ipD < gRelSigmaConv->GetN(); ++ipD) {
      double xD, yD; gRelSigmaConv->GetPoint(ipD, xD, yD);
      
      // Find matching MC point with same x-value (pT)
      int ipM = -1;
      for (int j = 0; j < gRelSigmaConvMC->GetN(); ++j) {
        double xM, yM; gRelSigmaConvMC->GetPoint(j, xM, yM);
        if (TMath::Abs(xD - xM) < 0.01) { // Allow small tolerance for floating point comparison
          ipM = j;
          break;
        }
      }
      
      if (ipM >= 0) {
        double xM, yM; gRelSigmaConvMC->GetPoint(ipM, xM, yM);
        double eDy = gRelSigmaConv->GetErrorY(ipD);
        double eMy = gRelSigmaConvMC->GetErrorY(ipM);
        double ratio = (yM != 0.0 ? yD / yM : 0.0);
        double eratio = 0.0;
        if (yD > 0.0 && yM > 0.0) {
          double relD = (eDy / yD);
          double relM = (eMy / yM);
          double var = relD*relD + relM*relM;
          eratio = ratio * (var > 0.0 ? TMath::Sqrt(var) : 0.0);
        }
        
        gRatioRMS->SetPoint(nRatioPointsRMS, xD, ratio);
        gRatioRMS->SetPointError(nRatioPointsRMS, 0.0, eratio);
        nRatioPointsRMS++;
      } else {
        printf("[Warning] RMS ratio: No matching MC point found for Data pT = %.2f GeV/c\n", xD);
      }
    }
    TCanvas *cRR = new TCanvas("c_ratio_rms", "Data/MC Ratio (RMS)", 800, 600);
    gRatioRMS->SetTitle("Data/MC ratio (RMS); p_{T} (GeV/c); (Data/MC) RMS/mean");
    gRatioRMS->SetMarkerStyle(25);
    gRatioRMS->SetMarkerColor(kRed+1);
    gRatioRMS->SetLineColor(kRed+1);
    gRatioRMS->GetYaxis()->SetRangeUser(0.0, 2.0);
    gRatioRMS->Draw("AP");
    cRR->SaveAs(Form("%s/V0MassFitResolutionVsPt_RatioRMS.pdf", outDir));
    delete cRR;
  }

  {
    TCanvas *c2 = new TCanvas("c_mean", "Mean vs pT", 800, 600);
    TString v0(v0Type);
    TString meanTitle;
    if (v0.EqualTo("K0S")) {
      meanTitle = "K^{0}_{S} mass mean; p_{T} (GeV/c); mean (GeV/c^{2})";
    } else if (v0.EqualTo("Lambda")) {
      meanTitle = "#Lambda mass mean; p_{T} (GeV/c); mean (GeV/c^{2})";
    } else if (v0.EqualTo("AntiLambda")) {
      meanTitle = "#bar{#Lambda} mass mean; p_{T} (GeV/c); mean (GeV/c^{2})";
    } else {
      meanTitle = "V0 mass mean; p_{T} (GeV/c); mean (GeV/c^{2})";
    }
    gMean->SetTitle(meanTitle.Data());
    gMean->SetMarkerStyle(21);
    gMean->Draw("AP");
    c2->SaveAs(Form("%s/V0MeanVsPt.pdf", outDir));
    delete c2;
  }

  // Store to a ROOT file
  TFile *fout = TFile::Open(Form("%s/V0MassFitResults.root", outDir), "RECREATE");
  if (fout && !fout->IsZombie()) {
    gSigma->Write();
    gMean->Write();
    gSigmaConv->Write();
    gRelSigma->Write();
    gRelSigmaConv->Write();
    
    // Write MC graphs if available
    if (gSigmaMC) {
      gSigmaMC->Write();
      gMeanMC->Write();
      gSigmaConvMC->Write();
      gRelSigmaMC->Write();
      gRelSigmaConvMC->Write();
    }
    
    fout->Close();
    delete fout;
  }

  f->Close();
  if (fMC) fMC->Close();
}

// Public interface functions
void DrawV0MassFits(const char *RunNumberData, const char *RunNumberMC) {
  // Defaults: fixed histogram paths per user's layout
  DrawV0MassFitsInternal(RunNumberData,
                        // "jet-fragmentation_id28293/data/V0/V0PtMass",
                        "jet-fragmentation_id28293/data/V0/V0PtMass",
                        "plots/V0MassFits",
                        "pol1ExpGausExp",
                        "K0S",
                        RunNumberMC,
                        "jet-fragmentation/mcd/V0/V0PtMass");
}

void DrawV0MassFits(const char *RunNumberData) {
  DrawV0MassFitsInternal(RunNumberData,
                         "jet-fragmentation_id28293/data/V0/V0PtMass",
                         "plots/V0MassFits",
                         "pol1ExpGausExp",
                         "K0S",
                         "",
                         "jet-fragmentation/mcd/V0/V0PtMass");
}

// Full control interface
void DrawV0MassFits(const char *RunNumber,
                    const char *histPath,
                    const char *outDir,
                    const char *fitModel,
                    const char *v0Type,
                    const char *mcRunNumber,
                    const char *mcHistPath) {
  DrawV0MassFitsInternal(RunNumber, histPath, outDir, fitModel, v0Type, mcRunNumber, mcHistPath);
}

#else

// This file is intended to be executed as a ROOT macro (cling). To run:
// root -l 'DrawV0MassFits.C("426828")'  // Data only (K0S default)
// root -l 'DrawV0MassFits.C("426828", "520161")'  // Data + MB MC (K0S default)
// root -l 'DrawV0MassFits.C("426828", "523547")'  // Data + JJ MC (K0S default)
// root -l 'DrawV0MassFits.C("426828", "jet-fragmentation_id28293/data/V0/V0PtMass", "plots/V0MassFits", "pol1ExpGausExp", "Lambda", "520161", "jet-fragmentation/mcd/V0/V0PtMass")'  // Full control

#endif


