{
//========= Macro generated from object: TPave/Legend of markers/lines/boxes to represent obj's
//========= by ROOT version6.32.06
   
   TLegend *leg = new TLegend(0.576961,0.16,0.847071,0.403478,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("LeadingJetPtRhoProj","Rhosparse Median","lpf");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = 1179;
   color = new TColor(ci, 1, 0, 0, " ", 0.3);
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("LeadingJetPtRhoMProj","Perpendicular Cone","lpf");

   ci = 1180;
   color = new TColor(ci, 0.8, 0, 0, " ", 0.3);
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#cc0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#cc0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
}
