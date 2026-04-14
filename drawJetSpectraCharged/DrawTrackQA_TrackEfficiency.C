///////////////////////////////////////////////////
///////////////////////////////////////////////////
////////// Track QA (track-efficiency)     //////////
////////// Data vs MC comparison           //////////
////////// author: Joonsuk Bae             //////////
////////// E-mail: jbae@cern.ch            //////////
////////// Last Modified: 2026-03-05       //////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////
//
// Usage: root -l DrawTrackQA_TrackEfficiency.C
//
// Consolidates DrawTrackQA.C, DrawTrackingEfficiency.C,
// and DrawTrackDCA.C into one modular macro.
//
// Configure files in DrawTrackQA_TrackEfficiencyConfig.h
// Toggle processes with TrackSelectionProcess, EfficiencyProcess, DCAProcess

#include "DrawTrackQA_TrackEfficiencyConfig.h"
#include "DrawTrackQA_TrackEfficiencyFunctions.h"

void DrawTrackQA_TrackEfficiency() {
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // Create output directory
    gSystem->mkdir(trkqa_outputDir, true);

    std::cerr << "========================================" << std::endl;
    std::cerr << " Track QA (track-efficiency)" << std::endl;
    std::cerr << " Data: " << trkqa_refPath << std::endl;
    std::cerr << " MC files: " << trkqa_fileNames.size() << std::endl;
    std::cerr << " Output: " << trkqa_outputDir << std::endl;
    std::cerr << "========================================" << std::endl;

    // Validate vector sizes
    int nMC = trkqa_fileNames.size();
    if ((int)trkqa_histNames.size() != nMC || (int)trkqa_McFileDirs.size() != nMC) {
        std::cerr << "[ERROR] Vector size mismatch: fileNames=" << nMC
                  << " histNames=" << trkqa_histNames.size()
                  << " McFileDirs=" << trkqa_McFileDirs.size() << std::endl;
        return;
    }
    if ((int)trkqa_Colors.size() < nMC + 1) {
        std::cerr << "[ERROR] Colors vector too small: need " << nMC + 1
                  << " (Data + " << nMC << " MC), have " << trkqa_Colors.size() << std::endl;
        return;
    }

    // Process 1: Track selection variables (1D, 2D, profiles)
    if (TrackSelectionProcess) {
        std::cerr << "\n--- Process 1: Track Selection ---" << std::endl;
        DrawTrackSelection(trkqa_fileNames, trkqa_histNames, trkqa_McFileDirs,
                            trkqa_Colors, trkqa_outputDir,
                            trkqa_refPath, trkqa_refDir, trkqa_DataLabel);
    }

    // Process 2: Tracking efficiency, fake rate, secondary contamination (MC only)
    if (EfficiencyProcess) {
        std::cerr << "\n--- Process 2: Efficiency ---" << std::endl;
        DrawTrackEfficiency(trkqa_fileNames, trkqa_histNames, trkqa_McFileDirs,
                             trkqa_Colors, trkqa_outputDir);
    }

    // Process 3: DCA detailed analysis
    if (DCAProcess) {
        std::cerr << "\n--- Process 3: DCA Analysis ---" << std::endl;
        DrawDCAAnalysis(trkqa_fileNames, trkqa_histNames, trkqa_McFileDirs,
                         trkqa_Colors, trkqa_outputDir,
                         trkqa_refPath, trkqa_refDir, trkqa_DataLabel);
    }

    std::cerr << "\n========================================" << std::endl;
    std::cerr << " All done! Output in: " << trkqa_outputDir << std::endl;
    std::cerr << "========================================" << std::endl;
}
