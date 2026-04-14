// Download Q/pT Correction Graphs from CCDB
// Purpose: Download TList (ccdb_object) from CCDB and save to ROOT file for verification
// Author: Auto-generated
// Date: 2026-01-29

#include "TFile.h"
#include "TList.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "CCDB/CcdbApi.h"
#include <iostream>
#include <map>
#include <string>

void DownloadFromCCDB(const char* ccdbPath,
                      const char* outputFile,
                      bool testCCDB = true,
                      long timestamp = -1)  // -1 means use current time
{
  std::cout << "=== Downloading from CCDB ===" << std::endl;
  std::cout << "CCDB path: " << ccdbPath << std::endl;
  std::cout << "Output file: " << outputFile << std::endl;
  std::cout << "Test CCDB: " << (testCCDB ? "Yes" : "No") << std::endl;
  
  // Determine host URL
  std::string hostUrl = testCCDB ? "http://ccdb-test.cern.ch:8080" : "http://alice-ccdb.cern.ch";
  
  // Get timestamp (use current time if not specified)
  if (timestamp < 0) {
    timestamp = TTimeStamp().GetSec() * 1000LL; // Current time in milliseconds
  }
  
  std::cout << "Host: " << hostUrl << std::endl;
  std::cout << "Timestamp: " << timestamp << " (" << TTimeStamp(timestamp/1000).AsString() << ")" << std::endl;
  
  // Initialize CCDB API
  o2::ccdb::CcdbApi ccdb;
  ccdb.init(hostUrl);
  
  // Prepare metadata (empty for retrieval)
  std::map<std::string, std::string> metadata;
  
  std::cout << "\n=== Retrieving TList from CCDB ===" << std::endl;
  
  // Retrieve TList from CCDB
  TList* ccdb_object = ccdb.retrieveFromTFileAny<TList>(ccdbPath, metadata, timestamp);
  
  if (!ccdb_object) {
    std::cerr << "❌ Error: Could not retrieve TList from CCDB!" << std::endl;
    std::cerr << "   Path: " << ccdbPath << std::endl;
    std::cerr << "   Timestamp: " << timestamp << std::endl;
    return;
  }
  
  std::cout << "✅ Successfully retrieved TList from CCDB!" << std::endl;
  
  // List contents
  std::cout << "\n=== TList Contents ===" << std::endl;
  ccdb_object->ls();
  
  // Verify contents
  TGraphErrors* gMC = (TGraphErrors*)ccdb_object->FindObject("sigmaVsPtMc");
  TGraphErrors* gData = (TGraphErrors*)ccdb_object->FindObject("sigmaVsPtData");
  
  if (!gMC || !gData) {
    std::cerr << "❌ Error: Required graphs not found in TList!" << std::endl;
    if (!gMC) std::cerr << "  - sigmaVsPtMc not found" << std::endl;
    if (!gData) std::cerr << "  - sigmaVsPtData not found" << std::endl;
    return;
  }
  
  std::cout << "\n✅ sigmaVsPtMc: " << gMC->GetN() << " points" << std::endl;
  std::cout << "✅ sigmaVsPtData: " << gData->GetN() << " points" << std::endl;
  
  // Save to ROOT file
  std::cout << "\n=== Saving to ROOT file ===" << std::endl;
  TFile* outFile = TFile::Open(outputFile, "RECREATE");
  if (!outFile || outFile->IsZombie()) {
    std::cerr << "❌ Error: Could not create output file " << outputFile << std::endl;
    return;
  }
  
  outFile->WriteObject(ccdb_object, "ccdb_object");
  outFile->Close();
  
  std::cout << "✅ Successfully saved to " << outputFile << std::endl;
  std::cout << "\n=== Verification ===" << std::endl;
  std::cout << "You can now verify the downloaded file using:" << std::endl;
  std::cout << "  root -l -b -q 'VerifyCCDBUpload.C(\"" << outputFile << "\")'" << std::endl;
}

// Convenience function to download both files
void DownloadQoverPtFiles(bool testCCDB = true)
{
  std::cout << "========================================" << std::endl;
  std::cout << "Download Q/pT Correction Files from CCDB" << std::endl;
  std::cout << "========================================" << std::endl;
  
  // CCDB paths
  const char* mbMCCcdbPath = "Users/j/jbae/qOverPtGraphs/Data_LHC22o_pass7_MC_LHC24f3c";
  const char* jjMCCcdbPath = "Users/j/jbae/qOverPtGraphs/Data_LHC22o_pass7_MC_LHC25a2b";
  
  // Output file names
  const char* mbMCOutputFile = "downloaded_MB_MC.root";
  const char* jjMCOutputFile = "downloaded_JJ_MC.root";
  
  if (testCCDB) {
    std::cout << "\n⚠️  Downloading from TEST CCDB!" << std::endl;
  } else {
    std::cout << "\n⚠️  Downloading from PRODUCTION CCDB!" << std::endl;
  }
  
  // // Download MB MC file
  // std::cout << "\n========================================" << std::endl;
  // std::cout << "Downloading MB MC file" << std::endl;
  // std::cout << "========================================" << std::endl;
  // DownloadFromCCDB(mbMCCcdbPath, mbMCOutputFile, testCCDB);
  
  // Download JJ MC file
  std::cout << "\n========================================" << std::endl;
  std::cout << "Downloading JJ MC file" << std::endl;
  std::cout << "========================================" << std::endl;
  DownloadFromCCDB(jjMCCcdbPath, jjMCOutputFile, testCCDB);
  
  std::cout << "\n========================================" << std::endl;
  std::cout << "Download complete!" << std::endl;
  std::cout << "========================================" << std::endl;
}

// Default execution (for direct ROOT macro call)
// Usage: root -l -b -q 'DownloadFromCCDB.C(true)'  // test CCDB
//        root -l -b -q 'DownloadFromCCDB.C(false)' // production CCDB
void DownloadFromCCDB(bool testCCDB = true) {
  DownloadQoverPtFiles(testCCDB);
}
