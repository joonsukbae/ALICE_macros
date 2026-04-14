// Draw Tracks
TH1 *DrawTrackPt(const char *fileName, const char *histName, Double_t Nevts,
                 TLegend *legend, Color_t colorID, const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "OPEN");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open file " << fileName << std::endl;
    return nullptr;
  }

  TH1* TrackPt = (TH1*) file->Get(Form("%s/%s", Dir, TrackPtObj));
  if (!TrackPt) {
    std::cerr << "[Warning] " << TrackPtObj << " not found in " << Dir << " of " << fileName << std::endl;
    return nullptr;
  }
  if (REBINON) {
    std::cout << "Rebinning histogram: " << histName << std::endl;
    TrackPt = TrackPt->Rebin(nTrackptbin, Form("TrackPt_%s", histName), Trackptbin);
  }

  legend->AddEntry(TrackPt, histName);
  hset(*TrackPt, TrackPtTitleX, TrackPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01, 0.01,
       0.05, 0.05, 510, 510);
  NORMEVENTS? hoptset(*TrackPt, Nevts, colorID, 0, 200, 2e-10, 1e0 + 0.05, 1, 1, 2, colorID==kBlack? 21 : 20) : hoptset(*TrackPt, Nevts, colorID, 0, 200, 2e-10, 1e0 + 0.05, 0.6, 1, 2, 24);

  std::cout << "Drawing histogram: " << histName << std::endl;
  TrackPt->Draw("pesame");

  return TrackPt;
}
TH1 *DrawTrackEta(const char *fileName, const char *histName, Double_t Nevts,
                  TLegend *legend, Color_t colorID, const char *Dir = nullptr) {
  // 파일 열기 체크
  auto file = TFile::Open(fileName, "open");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open file " << fileName << std::endl;
    return nullptr;
  }

  // 히스토그램 가져오기 체크
  TH2 *h2TrackEta = (TH2 *)file->Get(Form("%s/%s", Dir, TrackEtaObj));
  if (!h2TrackEta) {
    std::cerr << "Error: " << TrackEtaObj << " could not be found in directory " << Dir << " of file " << fileName << std::endl;
    return nullptr;
  }

  // Projection 체크
  TH1* TrackEta = (TH1*) h2TrackEta->ProjectionX(Form("TrackEta_%s", histName), 1, h2TrackEta->GetNbinsY(), "e");
  if (!TrackEta) {
    std::cerr << "Error: Failed to create projection for " << histName << std::endl;
    return nullptr;
  }


  // Rebin by 5: ensures bin edges at multiples of 0.1, aligning with track |eta|<0.9 boundary
  TH1 *h1TrackEta = TrackEta->Rebin(5, Form("TrackEta_rebinned_%s", histName));

  // 레전드와 스타일 설정
  legend->AddEntry(h1TrackEta, histName);
  hset(*h1TrackEta, TrackEtaTitleX, TrackEtaTitleY, 0.9, 1.4, 0.05, 0.05, 0.01,
       0.01, 0.05, 0.05, 510, 510);

  // 정규화 및 그리기 옵션 설정
  if (NORMEVENTS) {
    std::cout << "Normalizing with " << Nevts << " events" << std::endl;
    hoptset(*h1TrackEta, Nevts, colorID, -0.9, 0.9, 0, 16);
  } else {
    hoptset(*h1TrackEta, Nevts, colorID, -0.9, 0.9, 0., 1.5);
  }

  h1TrackEta->Draw("esame");
  std::cout << "Successfully drew track eta histogram for " << histName << std::endl;

  return h1TrackEta;
}
TH1 *DrawTrackPhi(const char *fileName, const char *histName, Double_t Nevts,
                  TLegend *legend, Color_t colorID, const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  TH2 *h2TrackPhi = (TH2 *)file->Get(Form("%s/%s", Dir, TrackPhiObj));
  if (!h2TrackPhi) {
    std::cerr << "[Warning] " << TrackPhiObj << " not found in " << Dir << " of " << fileName << std::endl;
    return nullptr;
  }
  TH1* TrackPhi = (TH1*) h2TrackPhi->ProjectionY(Form("TrackPhi_%s", histName), 1, h2TrackPhi->GetNbinsX(), "e");
  if (REBINON) {
    // TrackPhi->Rebin(10);
    std::vector<Double_t> newBins = NewBin(40, -1, 7);
    TrackPhi = TrackPhi->Rebin(40, Form("TrackPhi_%s", histName), newBins.data());
  }
  legend->AddEntry(TrackPhi, histName);
  hset(*TrackPhi, TrackPhiTitleX, TrackPhiTitleY, 0.9, 1.4, 0.05, 0.05, 0.01,
       0.01, 0.05, 0.05, 510, 510);
  NORMEVENTS? hoptset(*TrackPhi, Nevts, colorID, 0, 2 * TMath::Pi(), 0, 5) : hoptset(*TrackPhi, Nevts, colorID, 0, 2 * TMath::Pi(), 0, 0.3);
  TrackPhi->Draw("esame");

  return TrackPhi;
}
// Draw Constituents
TH1 *DrawConstituentPt(const char *fileName, const char *histName, Double_t Nevts,
                       TLegend *legend, Color_t colorID,
                       const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open file " << fileName << std::endl;
    return nullptr;
  }
  TH2 *H2ConstituentPt = (TH2 *)file->Get(Form("%s/%s", Dir, ConstPtObj));
  cout << "sucessfully loaded H2ConstituentPt from in " << Dir << endl;
  if (!H2ConstituentPt) {
    std::cerr << "Error: " << ConstPtObj << " could not be found in directory " << Dir << " of file " << fileName << std::endl;
    return nullptr;
  }
  TH1 *ConstituentPt = H2ConstituentPt->ProjectionY(Form("ConstituentPt_%s", histName), 1, H2ConstituentPt->GetNbinsX());
  cout << "sucessfully loaded ConstituentPt in " << Dir << endl;
  if (!ConstituentPt) {
    std::cerr << "Error: Failed to create projection for " << histName << std::endl;
    return nullptr;
  }
  if (REBINON) {
    ConstituentPt =
        ConstituentPt->Rebin(nTrackptbin, Form("ConstPt_%s", histName), Trackptbin);
  }
  legend->AddEntry(ConstituentPt, histName);
  hset(*ConstituentPt, ConstPtTitleX, ConstPtTitleY, 0.9, 1.4, 0.05, 0.05, 0.01,
       0.01, 0.05, 0.05, 510, 510);
  hoptset(*ConstituentPt, Nevts, colorID, 0, 100, 2e-10, 1);
  ConstituentPt->Draw("esame");

  return ConstituentPt;
}
TH1D *DrawConstituentEta(const char *fileName, const char *histName,
                         Double_t Nevts, TLegend *legend, Color_t colorID,
                         const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  TH3D *H3ConstituentEta = (TH3D *)file->Get(Form("%s/%s", Dir, ConstEtaObj));
  TH1D *ConstituentEta =
      H3ConstituentEta->ProjectionZ(Form("ConstituentEta_%s", histName),
                                    H3ConstituentEta->GetXaxis()->FindBin(RBIN),
                                    H3ConstituentEta->GetXaxis()->FindBin(RBIN),
                                    1,
                                    H3ConstituentEta->GetNbinsY());
  float nZconEta = H3ConstituentEta->GetNbinsZ();
  TH1D *h1ConstEta = new TH1D("ConstEta", "ConstEta", nZconEta, -1., 1.);
  for (int ieta = 1; ieta <= nZconEta; ieta++) {
    double etaContent = ConstituentEta->GetBinContent(ConstituentEta->FindBin(-1 + (2/nZconEta) * ieta));
    double etaError = ConstituentEta->GetBinError(ConstituentEta->FindBin(-1 + (2/nZconEta) * ieta));
    h1ConstEta->SetBinContent(ieta+1, etaContent);
    h1ConstEta->SetBinError(ieta+1, etaError);
  }

  legend->AddEntry(h1ConstEta, histName);
  hset(*h1ConstEta, ConstEtaTitleX, ConstEtaTitleY, 0.9, 1.4, 0.05, 0.05,
       0.01, 0.01, 0.05, 0.05, 510, 510);
  NORMEVENTS? hoptset(*h1ConstEta, Nevts, colorID, -0.5, 0.5, 0., 0.3) : hoptset(*h1ConstEta, Nevts, colorID, -0.5, 0.5, 0., 1.5); 
  h1ConstEta->Draw("esame");

  return h1ConstEta;
}
TH1 *DrawConstituentPhi(const char *fileName, const char *histName, Double_t Nevts,
                        TLegend *legend, Color_t colorID,
                        const char *Dir = nullptr) {
  auto file = TFile::Open(fileName, "open");
  TH3D *H3ConstituentPhi = (TH3D *)file->Get(Form("%s/%s", Dir, ConstPhiObj));
  TH1 *ConstituentPhi =
      H3ConstituentPhi->ProjectionZ(Form("ConstituentPhi_%s", histName),
                                    H3ConstituentPhi->GetXaxis()->FindBin(RBIN),
                                    H3ConstituentPhi->GetXaxis()->FindBin(RBIN),
                                    1, H3ConstituentPhi->GetNbinsY());
  if (REBINON) {
    ConstituentPhi = ConstituentPhi -> Rebin(2);
  }
  ConstituentPhi->Rebin(2);
  legend->AddEntry(ConstituentPhi, histName);
  hset(*ConstituentPhi, ConstPhiTitleX, ConstPhiTitleY, 0.9, 1.4, 0.05, 0.05,
       0.01, 0.01, 0.05, 0.05, 510, 510);
  NORMEVENTS? hoptset(*ConstituentPhi, Nevts, colorID, 0, 2 * TMath::Pi(), 0., 0.06) : hoptset(*ConstituentPhi, Nevts, colorID, 0, 2 * TMath::Pi(), 0., 0.3);
  ConstituentPhi->Draw("esame");

  return ConstituentPhi;
}
