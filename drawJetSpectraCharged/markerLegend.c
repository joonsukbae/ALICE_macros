{
//========= Macro generated from object: TPave/Legend of markers/lines/boxes to represent obj's
//========= by ROOT version6.36.04
   
   TLegend *leg = new TLegend(0.488294, 0.414815, 0.924749, 0.896296, nullptr, "brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.0466667);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *legentry = leg->AddEntry("dataR04_5","ALICE data","pe");
   legentry->SetMarkerStyle(20);
   legentry->SetMarkerSize(0.8);
   legentry->SetTextFont(42);
   legentry = leg->AddEntry("dR04_MC truth (PYTHIA8 Monash)_5","MC truth (PYTHIA8 Monash)","l");
   legentry->SetLineColor(TColor::GetColor("#666666"));
   legentry->SetLineStyle(7);
   legentry->SetLineWidth(2);
   legentry->SetTextFont(42);
   legentry = leg->AddEntry("dR04_PYTHIA8 Monash (HY)_5","PYTHIA8 Monash (HY)","l");
   legentry->SetLineColor(TColor::GetColor("#009999"));
   legentry->SetLineWidth(2);
   legentry->SetTextFont(42);
   legentry = leg->AddEntry("dR04_PYTHIA8 Rope (HY)_5","PYTHIA8 Rope (HY)","l");
   legentry->SetLineColor(TColor::GetColor("#ff3333"));
   legentry->SetLineWidth(2);
   legentry->SetTextFont(42);
   legentry = leg->AddEntry("dR04_PYTHIA8 Shoving (HY)_5","PYTHIA8 Shoving (HY)","l");
   legentry->SetLineColor(TColor::GetColor("#006600"));
   legentry->SetLineWidth(2);
   legentry->SetTextFont(42);
   legentry = leg->AddEntry("dR04_PYTHIA8 Monash (MB, KIAF)_5","PYTHIA8 Monash (MB, KIAF)","l");
   legentry->SetLineColor(TColor::GetColor("#0000ff"));
   legentry->SetLineWidth(2);
   legentry->SetTextFont(42);
   legentry = leg->AddEntry("dR04_PYTHIA8 Monash (JJ, KIAF)_5","PYTHIA8 Monash (JJ, KIAF)","l");
   legentry->SetLineColor(TColor::GetColor("#000066"));
   legentry->SetLineStyle(2);
   legentry->SetLineWidth(2);
   legentry->SetTextFont(42);
   legentry = leg->AddEntry("dR04_POWHEG NLO (Hadi)_5","POWHEG NLO (Hadi)","l");
   legentry->SetLineColor(TColor::GetColor("#9933ff"));
   legentry->SetLineStyle(3);
   legentry->SetLineWidth(3);
   legentry->SetTextFont(42);
   leg->Draw();
}
