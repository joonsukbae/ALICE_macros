///////////////////////////////////////////////////
////////// Test Macro for New Analysis       /////
///////////////////////////////////////////////////

#include "DrawJetsMC_New.C"

void TestNewAnalysis() {
    std::cout << "=== Testing New Analysis Framework ===" << std::endl;
    
    // 출력 디렉토리 설정
    SetOutputDirectory("plots/NewAnalysis_Test/");
    
    // 모든 프로세스 활성화
    EnableAllProcesses();
    
    // 단계별 테스트
    std::cout << "\n1. Testing Track pT comparison..." << std::endl;
    DrawTrackPtOnly();
    
    std::cout << "\n2. Testing Track Resolution comparison..." << std::endl;
    DrawTrackResolutionOnly();
    
    std::cout << "\n3. Testing Jet pT comparison..." << std::endl;
    DrawJetPtOnly();
    
    std::cout << "\n4. Testing Jet Eta comparison..." << std::endl;
    DrawJetEtaOnly();
    
    std::cout << "\n5. Testing Jet Phi comparison..." << std::endl;
    DrawJetPhiOnly();
    
    std::cout << "\n=== Test Completed ===" << std::endl;
}

void QuickTest() {
    std::cout << "=== Quick Test: Track pT Only ===" << std::endl;
    
    SetOutputDirectory("plots/QuickTest/");
    DisableAllProcesses();
    gPlotConfig.trackProcess = true;
    
    DrawTrackPtOnly();
    
    std::cout << "=== Quick Test Completed ===" << std::endl;
}