// DrawNormalizationSchematic.C
// Clean LaTeX-style normalization efficiency tables
// Usage: root -l DrawNormalizationSchematic.C

#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>

const char* kOutputDir = "plots/ModelComparisons";
const int kF = 42;
const int kFB = 62;

void T(double x, double y, const char* t, float s, int a, Color_t c, int f=42) {
  TLatex* l = new TLatex(x, y, t);
  l->SetNDC(false); l->SetTextSize(s); l->SetTextAlign(a);
  l->SetTextColor(c); l->SetTextFont(f); l->Draw();
}
void HL(double x1, double x2, double y, Color_t c=kBlack, int w=1) {
  TLine* l = new TLine(x1, y, x2, y);
  l->SetLineColor(c); l->SetLineWidth(w); l->Draw();
}

void DrawNormalizationSchematic() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gSystem->mkdir(kOutputDir, true);

  // ============================================================
  // Main canvas: Two tables side by side
  // ============================================================
  TCanvas* c = new TCanvas("eff", "Efficiency Tables", 1100, 750);
  c->SetFillColor(kWhite);
  c->Range(0, 0, 110, 75);

  // Title equation
  T(55, 73, "#bf{Normalization:}   #sigma = #frac{1}{L_{int}} #times #frac{N_{jet}^{sel}}{#Delta#it{p}_{T} #Delta#it{#eta}} #times #frac{1}{#varepsilon_{total}}", 0.025, 22, kBlack, kF);

  // ============================================================
  // LEFT TABLE: Data-Driven (608785, R=0.4)
  // ============================================================
  double LX = 2, LW = 50;
  double LY = 68;
  float rh = 3.2;  // row height
  float fs = 0.019; // font size
  float fsh = 0.017; // header font

  // Columns
  double lc1 = LX + 1.5;   // step name
  double lc2 = LX + 33;    // ε_step
  double lc3 = LX + 43;    // ε_cumul

  T(LX + LW/2, LY + 1.5, "Data-Driven (BC-level)", 0.024, 22, kBlue+2, kFB);
  T(LX + LW/2, LY - 0.5, "608785, R = 0.4", 0.015, 22, kGray+1, kF);

  double y = LY - 2;
  HL(LX, LX+LW, y, kBlack, 2);
  T(lc1, y-rh/2, "Selection step", fsh, 12, kBlack, kFB);
  T(lc2, y-rh/2, "#varepsilon_{step}", fsh, 22, kBlack, kFB);
  T(lc3, y-rh/2, "#varepsilon_{cumul}", fsh, 22, kBlack, kFB);
  y -= rh; HL(LX, LX+LW, y, kBlack, 2);

  // BC level
  T(lc1-0.5, y-rh/2, "#scale[0.8]{BC level}", fsh, 12, kGray+1, kFB);
  y -= rh*0.7;

  struct R { const char* n; const char* es; const char* ec; Color_t c; };
  R dRows[] = {
    {"TVX trigger",                    "in #sigma_{vis}", "in #sigma_{vis}", kGray+1},
    {"sel8 BC (TFB + ITSROF)",         "0.873",           "0.873",           kBlue+2},
    {"RCT quality (CBT hadronPID)",    "0.375",           "0.327",           kBlue+2},
  };
  for (int i = 0; i < 3; i++) {
    T(lc1, y-rh/2, dRows[i].n, fs, 12, kBlack, kF);
    T(lc2, y-rh/2, dRows[i].es, fs, 22, dRows[i].c, kF);
    T(lc3, y-rh/2, dRows[i].ec, fs, 22, dRows[i].c, kF);
    y -= rh; HL(LX+1, LX+LW-1, y, kGray, 1);
  }

  // Transition
  HL(LX, LX+LW, y, kGray+1, 1);
  T(lc1-0.5, y-rh/2, "#scale[0.8]{BC #rightarrow Collision}", fsh, 12, kGray+1, kFB);
  y -= rh*0.7;

  T(lc1, y-rh/2, "BC #rightarrow PV reco", fs, 12, kBlack, kF);
  T(lc2, y-rh/2, "0.951", fs, 22, kRed+1, kF);
  T(lc3, y-rh/2, "#color[2]{not applied}", 0.016, 22, kRed+1, kF);
  y -= rh; HL(LX+1, LX+LW-1, y, kGray, 1);

  // Collision level
  HL(LX, LX+LW, y, kGray+1, 1);
  T(lc1-0.5, y-rh/2, "#scale[0.8]{Collision level}", fsh, 12, kGray+1, kFB);
  y -= rh*0.7;

  T(lc1, y-rh/2, "z-vertex |z| < 10 cm", fs, 12, kBlack, kF);
  T(lc2, y-rh/2, "0.949", fs, 22, kBlue+2, kF);
  T(lc3, y-rh/2, "0.311", fs, 22, kBlue+2, kF);
  y -= rh; HL(LX+1, LX+LW-1, y, kGray, 1);

  // Total
  HL(LX, LX+LW, y, kBlack, 2);
  T(lc1, y-rh/2-0.3, "#bf{Total}", 0.022, 12, kBlack, kFB);
  T(lc3, y-rh/2-0.3, "#bf{0.311}", 0.024, 22, kBlue+2, kFB);
  y -= rh+0.5; HL(LX, LX+LW, y, kBlack, 2);

  T(LX+1, y-1.5, "Sources: eventselection-run3/luminosity/", 0.013, 12, kGray+1, kF);
  T(LX+1, y-3.2, "#sigma_{vis} = 53.4 mb (vdM 2023)", 0.013, 12, kGray+1, kF);

  // ============================================================
  // RIGHT TABLE: MC-Driven (605134, sel8)
  // ============================================================
  double RX = 56, RW = 52;
  double RY = 68;

  double rc1 = RX + 1.5;
  double rc2 = RX + 21;  // ε_step^evt
  double rc3 = RX + 29;  // ε_cum^evt
  double rc4 = RX + 38;  // ε_step^jet
  double rc5 = RX + 47;  // ε_cum^jet

  T(RX + RW/2, RY + 1.5, "MC-Driven (event & jet level)", 0.024, 22, kRed+1, kFB);
  T(RX + RW/2, RY - 0.5, "605134, sel8, R = 0.4", 0.015, 22, kGray+1, kF);

  y = RY - 2;
  HL(RX, RX+RW, y, kBlack, 2);

  // Two-row header
  T(rc1, y-rh*0.35, "Selection", 0.016, 12, kBlack, kFB);
  T(rc1, y-rh*0.75, "step", 0.016, 12, kBlack, kFB);
  T((rc2+rc3)/2, y-rh*0.3, "Event", 0.015, 22, kRed+1, kFB);
  T(rc2, y-rh*0.7, "#varepsilon_{s}", 0.015, 22, kRed+1, kF);
  T(rc3, y-rh*0.7, "#varepsilon_{c}", 0.015, 22, kRed+1, kF);
  T((rc4+rc5)/2, y-rh*0.3, "Jet [20-30]", 0.015, 22, kBlue+2, kFB);
  T(rc4, y-rh*0.7, "#varepsilon_{s}", 0.015, 22, kBlue+2, kF);
  T(rc5, y-rh*0.7, "#varepsilon_{c}", 0.015, 22, kBlue+2, kF);
  y -= rh; HL(RX, RX+RW, y, kBlack, 2);

  struct MR { const char* n; const char* se; const char* ce; const char* sj; const char* cj; };
  MR mRows[] = {
    {"INEL (all)",       "1.000", "1.000", "1.000", "1.000"},
    {"Has reco coll.",   "0.748", "0.748", "0.951", "0.951"},
    {"Split coll.",      "0.994", "0.744", "0.982", "0.935"},
    {"TVX trigger",      "0.898", "0.668", "0.987", "0.923"},
    {"TF border",        "0.975", "0.651", "0.973", "0.897"},
    {"ITS ROF border",   "0.845", "0.550", "0.844", "0.757"},
    {"z-vtx |z|<10",     "0.806", "0.443", "1.000", "0.757"},
  };
  for (int i = 0; i < 7; i++) {
    T(rc1, y-rh/2, mRows[i].n, fs, 12, kBlack, kF);
    T(rc2, y-rh/2, mRows[i].se, 0.017, 22, kRed+1, kF);
    T(rc3, y-rh/2, mRows[i].ce, 0.017, 22, kRed+1, kF);
    T(rc4, y-rh/2, mRows[i].sj, 0.017, 22, kBlue+2, kF);
    T(rc5, y-rh/2, mRows[i].cj, 0.017, 22, kBlue+2, kF);
    y -= rh; HL(RX+1, RX+RW-1, y, kGray, 1);
  }

  // Total
  HL(RX, RX+RW, y, kBlack, 2);
  T(rc1, y-rh/2-0.3, "#bf{Total}", 0.020, 12, kBlack, kFB);
  T(rc3, y-rh/2-0.3, "#bf{0.443}", 0.022, 22, kRed+1, kFB);
  T(rc5, y-rh/2-0.3, "#bf{0.757}", 0.022, 22, kBlue+2, kFB);
  y -= rh+0.5; HL(RX, RX+RW, y, kBlack, 2);

  T(RX+1, y-1.5, "Source: jet-cross-section-efficiency (605134)", 0.013, 12, kGray+1, kF);
  T(RX+1, y-3.2, "#varepsilon_{zvtx}^{jet} = 1: jet finder already applied", 0.013, 12, kGray+1, kF);
  T(RX+1, y-4.9, "#sigma_{INEL} = 78.58 mb (PYTHIA8 Monash)", 0.013, 12, kGray+1, kF);

  // ============================================================
  // Bottom: Comparison
  // ============================================================
  HL(5, 105, 6, kBlack, 2);
  T(55, 4, "#varepsilon_{total}^{data} = #bf{0.311}          #varepsilon_{total}^{MC,jet} = #bf{0.757}          #varepsilon_{total}^{MC,evt} = #bf{0.443}", 0.022, 22, kBlack, kF);
  T(55, 1.5, "#varepsilon^{data} captures sel8+RCT+zvtx (BC-level).  #varepsilon^{MC} captures reco+TVX+sel8+zvtx (event/jet-level).", 0.015, 22, kGray+2, kF);
  HL(5, 105, 0.5, kBlack, 2);

  c->Print(Form("%s/NormalizationSchematic_EffTable.pdf", kOutputDir));
  std::cout << "Saved: " << kOutputDir << "/NormalizationSchematic_EffTable.pdf" << std::endl;
}
