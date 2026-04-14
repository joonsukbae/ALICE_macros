// Upload an EXTREME Q/pT correction graphs file to CCDB (for TrackTuner debugging).
//
// This is a thin wrapper around UploadQoverPtToCCDB.C.
// It uploads the file produced by CreateExtremeQoverPtGraphsFile.C to a NEW CCDB path:
//   Users/j/jbae/qOverPtGraphs/Data_LHC22o_pass7_MC_LHC25a2b_extreme_<scaleTag>
//
// Usage examples:
//   root -l -b -q 'UploadExtremeQoverPtToCCDB.C(true,  5.0)'  # test CCDB
//   root -l -b -q 'UploadExtremeQoverPtToCCDB.C(false, 5.0)'  # production CCDB

#include "TString.h"
#include <cmath>
#include <iostream>

// Reuse the existing uploader implementation.
#include "UploadQoverPtToCCDB.C"

static TString MakeScaleTagExt(double scale)
{
  TString tag;
  if (std::abs(scale - std::round(scale)) < 1e-9) {
    tag = Form("x%.0f", std::round(scale));
  } else {
    tag = Form("x%.3f", scale);
    while (tag.EndsWith("0")) {
      tag.Chop();
    }
    if (tag.EndsWith(".")) {
      tag.Chop();
    }
  }
  tag.ReplaceAll(".", "p");
  return tag;
}

void UploadExtremeQoverPtToCCDB(bool testCCDB = true, double scale = 5.0)
{
  if (!(scale > 0.0)) {
    std::cerr << "[Error] scale must be > 0, got " << scale << std::endl;
    return;
  }

  const TString scaleTag = MakeScaleTagExt(scale);
  const TString rootFile = Form("Sigma1overPt_Data_LHC22o_pass7_globalTracks_MC_LHC25a2b_EXTREME_%s.root", scaleTag.Data());
  const TString ccdbPath = Form("Users/j/jbae/qOverPtGraphs/Data_LHC22o_pass7_MC_LHC25a2b_extreme_%s", scaleTag.Data());
  const TString desc = Form("EXTREME DEBUG qOverPtGraphs: sigmaVsPtData = %.6g * sigmaVsPtMc (constant)", scale);

  std::cout << "=== UploadExtremeQoverPtToCCDB ===" << std::endl;
  std::cout << "File:     " << rootFile << std::endl;
  std::cout << "CCDB path: " << ccdbPath << std::endl;
  std::cout << "Target:   " << (testCCDB ? "TEST" : "PRODUCTION") << std::endl;
  std::cout << "Desc:     " << desc << std::endl;

  UploadQoverPtToCCDB(rootFile.Data(), ccdbPath.Data(), testCCDB, "Joonsuk Bae", desc.Data());
}
