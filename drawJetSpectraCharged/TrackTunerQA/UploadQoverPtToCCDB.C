// Upload Q/pT Correction Graphs to CCDB
// Purpose: Upload TList (ccdb_object) containing sigmaVsPtMc and sigmaVsPtData to CCDB
// Uses CcdbApi::storeAsTFileAny to properly set File name field
// Author: Auto-generated
// Date: 2026-01-29

#include "TFile.h"
#include "TList.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "CCDB/CcdbApi.h"
#include <iostream>
#include <map>
#include <string>

void UploadQoverPtToCCDB(const char* rootFilePath, 
                         const char* ccdbPath,
                         bool testCCDB = true,
                         const char* uploadedBy = "Joonsuk Bae",
                         const char* description = "Q/pT correction graphs from covariance matrix")
{
  std::cout << "=== Uploading Q/pT Correction File to CCDB ===" << std::endl;
  std::cout << "ROOT file: " << rootFilePath << std::endl;
  std::cout << "CCDB path: " << ccdbPath << std::endl;
  std::cout << "Test CCDB: " << (testCCDB ? "Yes" : "No") << std::endl;
  
  // Verify file exists
  if (gSystem->AccessPathName(rootFilePath)) {
    std::cerr << "[Error] File does not exist: " << rootFilePath << std::endl;
    return;
  }
  
  // Open ROOT file to verify structure
  TFile* inputFile = TFile::Open(rootFilePath, "READ");
  if (!inputFile || inputFile->IsZombie() || !inputFile->IsOpen()) {
    std::cerr << "[Error] Could not open ROOT file: " << rootFilePath << std::endl;
    return;
  }
  
  std::cout << "\n=== Verifying file structure ===" << std::endl;
  
  // Get ccdb_object TList
  TList* ccdb_object = (TList*)inputFile->Get("ccdb_object");
  if (!ccdb_object) {
    std::cerr << "[Error] Could not find 'ccdb_object' in ROOT file" << std::endl;
    inputFile->Close();
    return;
  }
  
  std::cout << "[OK] Found ccdb_object TList" << std::endl;
  
  // Verify contents
  TGraphErrors* gMC = (TGraphErrors*)ccdb_object->FindObject("sigmaVsPtMc");
  TGraphErrors* gData = (TGraphErrors*)ccdb_object->FindObject("sigmaVsPtData");
  
  if (!gMC || !gData) {
    std::cerr << "[Error] Required graphs not found in ccdb_object!" << std::endl;
    if (!gMC) std::cerr << "  - sigmaVsPtMc not found" << std::endl;
    if (!gData) std::cerr << "  - sigmaVsPtData not found" << std::endl;
    inputFile->Close();
    return;
  }
  
  std::cout << "[OK] sigmaVsPtMc: " << gMC->GetN() << " points" << std::endl;
  std::cout << "[OK] sigmaVsPtData: " << gData->GetN() << " points" << std::endl;
  
  inputFile->Close();
  
  // Extract file name from path
  TString rootFilePathStr(rootFilePath);
  TString fileName = gSystem->BaseName(rootFilePathStr.Data());
  
  // Prepare metadata for CcdbApi
  std::map<std::string, std::string> metadata;
  metadata["UploadedBy"] = uploadedBy;
  metadata["Email"] = "jbae@cern.ch";
  metadata["Description"] = description;
  metadata["ObjectType"] = "TList";
  metadata["Contains"] = "sigmaVsPtMc,sigmaVsPtData";
  metadata["FileName"] = fileName.Data();  // Add file name to metadata
  
  // Determine host URL
  std::string hostUrl = testCCDB ? "http://ccdb-test.cern.ch:8080" : "http://alice-ccdb.cern.ch";

  // Validity range
  // Use a wide validity range so that the object is found regardless of the CCDB timestamp
  // used in grid/hyperloop jobs (often tied to run time rather than upload time).
  // This mirrors the common convention used by existing TrackTuner reference objects.
  // NOTE: CCDB treats Valid-from=0 specially in some setups; use 1ms like many reference objects.
  Long64_t startTimestampMs = 1; // 01 Jan 1970 00:00:00 UTC + 1ms
  Long64_t endTimestampMs = 2524604400000LL; // 01 Jan 2050 00:00:00 CET (approx), as used in CCDB listings
  
  std::cout << "\n=== Uploading to CCDB using CcdbApi ===" << std::endl;
  std::cout << "Host: " << hostUrl << std::endl;
  std::cout << "CCDB Path: " << ccdbPath << std::endl;
  std::cout << "File Name: " << fileName << std::endl;
  std::cout << "Start timestamp: " << startTimestampMs << " (" << TTimeStamp(startTimestampMs/1000).AsString() << ")" << std::endl;
  std::cout << "End timestamp: " << endTimestampMs << " (" << TTimeStamp(endTimestampMs/1000).AsString() << ")" << std::endl;
  
  // Initialize CCDB API
  o2::ccdb::CcdbApi ccdb;
  ccdb.init(hostUrl);
  
  // Upload TList using storeAsTFileAny
  // Note: storeAsTFileAny may not populate File name field, but we include it in metadata
  int res = ccdb.storeAsTFileAny(ccdb_object, ccdbPath, metadata, startTimestampMs, endTimestampMs);
  
  if (res != 0) {
    std::cerr << "[Error] Upload to CCDB failed. Return code: " << res << std::endl;
  } else {
    std::cout << "[OK] Successfully uploaded to CCDB." << std::endl;
    std::cout << "   Path: " << ccdbPath << std::endl;
    if (testCCDB) {
      std::cout << "   Test URL: http://ccdb-test.cern.ch:8080/browse/" << ccdbPath << std::endl;
    } else {
      std::cout << "   Production URL: http://alice-ccdb.cern.ch/browse/" << ccdbPath << std::endl;
    }
    std::cout << "   Note: File name field may still be empty in web interface," << std::endl;
    std::cout << "         but TrackTuner can find files using path and nameFileQoverPt." << std::endl;
  }
}

// Main function to upload both files
void UploadQoverPtFiles(bool testCCDB = true)
{
  std::cout << "========================================" << std::endl;
  std::cout << "Upload Q/pT Correction Files to CCDB" << std::endl;
  std::cout << "========================================" << std::endl;
  
  // File paths: Files are now in the same directory as the macro
  const char* mbMCFile = "Sigma1overPt_Data_LHC22o_pass7_globalTracks_MC_LHC24f3c.root";
  const char* jjMCFile = "Sigma1overPt_Data_LHC22o_pass7_globalTracks_MC_LHC25a2b.root";
  
  // CCDB paths: Each file goes into its own subdirectory to avoid conflicts
  // This ensures TrackTuner can find the correct file even if File name field is empty
  // Format: Users/j/jbae/qOverPtGraphs/<subdirectory>/
  // TrackTuner will use: pathFileQoverPt = "Users/j/jbae/qOverPtGraphs/<subdirectory>"
  //                      nameFileQoverPt = "Sigma1overPt_Data_LHC22o_pass7_globalTracks_MC_<MC_TAG>.root"
  const char* mbMCCcdbPath = "Users/j/jbae/qOverPtGraphs/Data_LHC22o_pass7_MC_LHC24f3c";
  const char* jjMCCcdbPath = "Users/j/jbae/qOverPtGraphs/Data_LHC22o_pass7_MC_LHC25a2b";
  
  if (testCCDB) {
    std::cout << "\n[Info] Uploading to TEST CCDB first." << std::endl;
    std::cout << "   Please verify the upload before using production CCDB." << std::endl;
  } else {
    std::cout << "\n[Warning] Uploading to PRODUCTION CCDB." << std::endl;
    std::cout << "   Make sure you have tested in test CCDB first!" << std::endl;
  }
  
  // // Upload MB MC file
  // std::cout << "\n========================================" << std::endl;
  // std::cout << "Uploading MB MC file" << std::endl;
  // std::cout << "========================================" << std::endl;
  // UploadQoverPtToCCDB(mbMCFile, mbMCCcdbPath, testCCDB, 
  //                     "Joonsuk Bae", 
  //                     "Q/pT correction graphs: Data (LHC22o-pass7) vs MB MC (LHC24f3c)");
  
  // Upload JJ MC file
  std::cout << "\n========================================" << std::endl;
  std::cout << "Uploading JJ MC file" << std::endl;
  std::cout << "========================================" << std::endl;
  UploadQoverPtToCCDB(jjMCFile, jjMCCcdbPath, testCCDB,
                      "Joonsuk Bae",
                      "Q/pT correction graphs: Data (LHC22o-pass7) vs JJ MC (LHC25a2b)");
  
  std::cout << "\n========================================" << std::endl;
  std::cout << "Upload complete!" << std::endl;
  std::cout << "========================================" << std::endl;
}

// Default execution (for direct ROOT macro call)
// Usage: root -l -b -q 'UploadQoverPtToCCDB.C(true)'  // test CCDB
//        root -l -b -q 'UploadQoverPtToCCDB.C(false)' // production CCDB
void UploadQoverPtToCCDB(bool testCCDB = true) {
  UploadQoverPtFiles(testCCDB);
}
