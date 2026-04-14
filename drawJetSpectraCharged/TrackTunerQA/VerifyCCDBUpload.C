// Verify CCDB Upload
// Purpose: Verify that the file downloaded from CCDB contains the correct structure
// Author: Auto-generated
// Date: 2026-01-29

#include "TFile.h"
#include "TList.h"
#include "TGraphErrors.h"
#include <iostream>
#include <string>

void VerifyCCDBUpload(const char* downloadedFile, 
                      const char* expectedPath = "",
                      bool verbose = true)
{
  std::cout << "========================================" << std::endl;
  std::cout << "Verifying CCDB Upload" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Downloaded file: " << downloadedFile << std::endl;
  if (strlen(expectedPath) > 0) {
    std::cout << "Expected CCDB path: " << expectedPath << std::endl;
  }
  std::cout << std::endl;
  
  // Open downloaded file
  TFile* file = TFile::Open(downloadedFile, "READ");
  if (!file || file->IsZombie() || !file->IsOpen()) {
    std::cerr << "❌ Error: Cannot open file " << downloadedFile << std::endl;
    return;
  }
  
  std::cout << "✅ File opened successfully" << std::endl;
  std::cout << "   File size: " << file->GetSize() << " bytes" << std::endl;
  
  // List all objects in the file
  if (verbose) {
    std::cout << "\n=== File Contents ===" << std::endl;
    file->ls();
  }
  
  // Get ccdb_object TList
  TList* ccdb_object = (TList*)file->Get("ccdb_object");
  if (!ccdb_object) {
    std::cerr << "❌ Error: Cannot find 'ccdb_object' TList in file" << std::endl;
    file->Close();
    return;
  }
  
  std::cout << "\n✅ Found 'ccdb_object' TList" << std::endl;
  
  if (verbose) {
    std::cout << "\n=== ccdb_object Contents ===" << std::endl;
    ccdb_object->ls();
  }
  
  // Verify required graphs
  TGraphErrors* gMC = (TGraphErrors*)ccdb_object->FindObject("sigmaVsPtMc");
  TGraphErrors* gData = (TGraphErrors*)ccdb_object->FindObject("sigmaVsPtData");
  
  bool allOK = true;
  
  if (!gMC) {
    std::cerr << "❌ Error: sigmaVsPtMc not found!" << std::endl;
    allOK = false;
  } else {
    std::cout << "\n✅ sigmaVsPtMc found:" << std::endl;
    std::cout << "   Points: " << gMC->GetN() << std::endl;
    if (gMC->GetN() > 0) {
      Double_t x, y;
      gMC->GetPoint(0, x, y);
      std::cout << "   First point: pT = " << x << " GeV/c, sigma = " << y << std::endl;
      gMC->GetPoint(gMC->GetN()-1, x, y);
      std::cout << "   Last point: pT = " << x << " GeV/c, sigma = " << y << std::endl;
    }
  }
  
  if (!gData) {
    std::cerr << "❌ Error: sigmaVsPtData not found!" << std::endl;
    allOK = false;
  } else {
    std::cout << "\n✅ sigmaVsPtData found:" << std::endl;
    std::cout << "   Points: " << gData->GetN() << std::endl;
    if (gData->GetN() > 0) {
      Double_t x, y;
      gData->GetPoint(0, x, y);
      std::cout << "   First point: pT = " << x << " GeV/c, sigma = " << y << std::endl;
      gData->GetPoint(gData->GetN()-1, x, y);
      std::cout << "   Last point: pT = " << x << " GeV/c, sigma = " << y << std::endl;
    }
  }
  
  // Check pT ranges
  if (gMC && gData && gMC->GetN() > 0 && gData->GetN() > 0) {
    Double_t xMC_min, yMC, xMC_max;
    Double_t xData_min, yData, xData_max;
    gMC->GetPoint(0, xMC_min, yMC);
    gMC->GetPoint(gMC->GetN()-1, xMC_max, yMC);
    gData->GetPoint(0, xData_min, yData);
    gData->GetPoint(gData->GetN()-1, xData_max, yData);
    
    std::cout << "\n=== pT Range Comparison ===" << std::endl;
    std::cout << "MC pT range: " << xMC_min << " - " << xMC_max << " GeV/c" << std::endl;
    std::cout << "Data pT range: " << xData_min << " - " << xData_max << " GeV/c" << std::endl;
    
    if (xMC_max < 50.0 || xData_max < 50.0) {
      std::cout << "⚠️  Warning: pT range may be too small for jet analysis (need up to ~200 GeV/c)" << std::endl;
    }
  }
  
  // Summary
  std::cout << "\n========================================" << std::endl;
  if (allOK) {
    std::cout << "✅ VERIFICATION PASSED" << std::endl;
    std::cout << "   File structure is correct" << std::endl;
    std::cout << "   Ready for production CCDB upload" << std::endl;
  } else {
    std::cout << "❌ VERIFICATION FAILED" << std::endl;
    std::cout << "   File structure is incorrect" << std::endl;
    std::cout << "   Do NOT upload to production CCDB" << std::endl;
  }
  std::cout << "========================================" << std::endl;
  
  file->Close();
}

// Convenience function to verify both files
void VerifyBothCCDBUploads()
{
  std::cout << "========================================" << std::endl;
  std::cout << "Verifying Both CCDB Uploads" << std::endl;
  std::cout << "========================================" << std::endl;
  
  // You need to download the files from CCDB first and provide the paths
  std::cout << "\n⚠️  Please download files from CCDB and provide paths:" << std::endl;
  std::cout << "   Example:" << std::endl;
  // std::cout << "   VerifyCCDBUpload(\"~/Downloads/[ID1]\", \"Users/j/jbae/qOverPtGraphs/Sigma1overPt_Data_LHC22o-pass7_globalTracks_MC_LHC24f3c\");" << std::endl;
  std::cout << "   VerifyCCDBUpload(\"~/Downloads/[ID2]\", \"Users/j/jbae/qOverPtGraphs/Sigma1overPt_Data_LHC22o-pass7_globalTracks_MC_LHC25a2b\");" << std::endl;
}
