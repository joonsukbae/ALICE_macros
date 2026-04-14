// ============================================
// UE-subtracted Analysis Entry Point (wrapper)
// ============================================
//
// This file exists so you can run UE-subtracted mode directly with:
//   root -l DrawJetsMCRDependentUE.C
//
// The actual implementation lives in DrawJetsMCRDependent.C.

{
  gROOT->LoadMacro("DrawJetsMCRDependent.C+");
  DrawJetsMCRDependentUE();
}
