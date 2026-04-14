{
//========= Macro generated from object: TPave/Legend of markers/lines/boxes to represent obj's
//========= by ROOT version6.36.04
   
   TLegend *leg = new TLegend(0.20976, 0.712003, 0.792212, 0.88709, nullptr, "brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *legentry = leg->AddEntry("dataRR_0_3","ALICE data","lpe");
   legentry->SetLineWidth(2);
   legentry->SetMarkerStyle(20);
   legentry->SetMarkerSize(0.6);
   legentry->SetTextFont(42);
   legentry = leg->AddEntry("dRR_MonashMB_3","PYTHIA8 Monash (MB, KIAF)","l");
   legentry->SetLineColor(TColor::GetColor("#00cccc"));
   legentry->SetLineStyle(3);
   legentry->SetLineWidth(2);
   legentry->SetTextFont(42);
   legentry = leg->AddEntry("dRR_Monash_3","PYTHIA8 Monash (JJ, KIAF)","l");
   legentry->SetLineColor(TColor::GetColor("#ff0000"));
   legentry->SetLineStyle(9);
   legentry->SetLineWidth(2);
   legentry->SetTextFont(42);
   legentry = leg->AddEntry("dRRPW_3","POWHEG NLO","l");
   legentry->SetLineColor(TColor::GetColor("#cc00cc"));
   legentry->SetLineStyle(7);
   legentry->SetLineWidth(2);
   legentry->SetTextFont(42);
   leg->Draw();
}
