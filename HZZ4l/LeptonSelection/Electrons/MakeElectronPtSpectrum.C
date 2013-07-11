//================================================================================================
//
// HWW selection macro
//
//  * plots distributions associated with selected events
//  * prints list of selected events from data
//  * outputs ROOT files of events passing selection for each sample, 
//    which can be processed by plotSelect.C
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 2D histograms
#include <TGraphAsymmErrors.h>     
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

// // define structures to read in ntuple
// #include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
// #include "EWKAna/Ntupler/interface/TEventInfo.hh"
// #include "EWKAna/Ntupler/interface/TElectron.hh"
// #include "EWKAna/Ntupler/interface/TPhoton.hh"
// #include "EWKAna/Ntupler/interface/TMuon.hh"
// #include "EWKAna/Ntupler/interface/TJet.hh"

// // lumi section selection with JSON files
// #include "MitCommon/DataFormats/interface/Types.h"
// #include "MitAna/DataCont/interface/RunLumiRangeMap.h"
// #include "MitCommon/MathTools/interface/MathUtils.h"
// #include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
// #include "MitHiggs/Utils/interface/EfficiencyUtils.h"
// #include "MitHiggs/Utils/interface/PlotUtils.h"

#include "HiggsAna/CommonData/interface/ElectronTree.h"

#endif

//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}

//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void FillLeptonPtSpectrum()
{  


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *histRealElePtSource = new TH1F("RealElePtSource", "; p_{T} [GeV/c] ; Number of Events ",  50, 0 , 100);
  TH1F *histRealElePtTarget = new TH1F("RealElePtTarget", "; p_{T} [GeV/c] ; Number of Events ",  50, 0 , 100);
  TH1F *histFakeElePtSource = new TH1F("FakeElePtSource", "; p_{T} [GeV/c] ; Number of Events ",  50, 0 , 100);


  //*****************************************************************************************
  //Setup
  //*****************************************************************************************

  //Variables
  Float_t                 fWeight;
  UInt_t                  fRunNumber;
  UInt_t                  fLumiSectionNumber;
  UInt_t                  fEventNumber;
  Float_t                 fElePt; 
  Float_t                 fEleEta; 
  Float_t                 fElePhi; 
  Float_t                 fEleSCEt; 
  Float_t                 fEleSCEta; 
  Float_t                 fEleSCPhi; 
  Float_t                 fElePFIso; 
  
  //CutBased Variables
  Float_t                 fEleSigmaIEtaIEta; 
  Float_t                 fEleDEtaIn; 
  Float_t                 fEleDPhiIn; 
  Float_t                 fEleHoverE; 
  Float_t                 fEleD0; 
  Float_t                 fEleDZ; 
  Float_t                 fEleFBrem; 
  Float_t                 fEleEOverP; 

  //Additional Vars used in Likelihood
  Float_t                 fEleESeedClusterOverPout; 
  Float_t                 fEleSigmaIPhiIPhi; 
  Float_t                 fEleNBrem; 
  Float_t                 fEleOneOverEMinusOneOverP; 
  Float_t                 fEleESeedClusterOverPIn; 
  Float_t                 fEleIP3d; 
  Float_t                 fEleIP3dSig; 

  Float_t                 fEleHcalDepth1OverEcal;
  Float_t                 fEleHcalDepth2OverEcal;
  Float_t                 fEledEtaCalo;
  Float_t                 fEledPhiCalo;
  Float_t                 fElePreShowerOverRaw;
  Float_t                 fEleCovIEtaIPhi;
  Float_t                 fEleSCEtaWidth;
  Float_t                 fEleSCPhiWidth;
  Float_t                 fEleGsfTrackChi2OverNdof;
  Float_t                 fEleR9;

  Float_t                 fEleSeedEMaxOverE;
  Float_t                 fEleSeedETopOverE;
  Float_t                 fEleSeedEBottomOverE;
  Float_t                 fEleSeedELeftOverE;
  Float_t                 fEleSeedERightOverE;
  Float_t                 fEleSeedE2ndOverE;
  Float_t                 fEleSeedE2x5RightOverE;
  Float_t                 fEleSeedE2x5LeftOverE;
  Float_t                 fEleSeedE2x5TopOverE;
  Float_t                 fEleSeedE2x5BottomOverE;
  Float_t                 fEleSeedE2x5MaxOverE;
  Float_t                 fEleSeedE1x3OverE;
  Float_t                 fEleSeedE3x1OverE;
  Float_t                 fEleSeedE1x5OverE;
  Float_t                 fEleSeedE2x2OverE;
  Float_t                 fEleSeedE3x2OverE;
  Float_t                 fEleSeedE3x3OverE;
  Float_t                 fEleSeedE4x4OverE;
  Float_t                 fEleSeedE5x5OverE;

  //Isolation Variables
  Float_t                 fEleChargedIso03; 
  Float_t                 fEleNeutralHadronIso03; 
  Float_t                 fEleGammaIso03; 
  Float_t                 fEleChargedIso04; 
  Float_t                 fEleNeutralHadronIso04; 
  Float_t                 fEleGammaIso04; 
  Float_t                 fEleChargedIso04FromOtherVertices; 
  Float_t                 fEleNeutralHadronIso04_10Threshold; 
  Float_t                 fEleGammaIso04_10Threshold; 
  Float_t                 fEleTrkIso03; 
  Float_t                 fEleEMIso03; 
  Float_t                 fEleHadIso03; 
  Float_t                 fEleTrkIso04; 
  Float_t                 fEleEMIso04; 
  Float_t                 fEleHadIso04; 
  Float_t                 fRho; 
  Float_t                 fNVertices; 


  
  //*****************************************************************************************
  //RealEleSourceTree
  //*****************************************************************************************
  TFile *RealEleSourceFile = new TFile("ElectronSelectionTraining.Real.root", "READ");
  TTree *RealEleSourceTree = (TTree*)RealEleSourceFile->Get("Electrons");
  RealEleSourceTree->SetBranchAddress( "weight", &fWeight);
  RealEleSourceTree->SetBranchAddress( "run", &fRunNumber);
  RealEleSourceTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  RealEleSourceTree->SetBranchAddress( "event", &fEventNumber);
  RealEleSourceTree->SetBranchAddress( "pt", &fElePt); 
  RealEleSourceTree->SetBranchAddress( "eta", &fEleEta); 
  RealEleSourceTree->SetBranchAddress( "phi", &fElePhi); 
  RealEleSourceTree->SetBranchAddress( "scet", &fEleSCEt); 
  RealEleSourceTree->SetBranchAddress( "sceta", &fEleSCEta); 
  RealEleSourceTree->SetBranchAddress( "scphi", &fEleSCPhi); 
  RealEleSourceTree->SetBranchAddress( "pfiso", &fElePFIso); 
  RealEleSourceTree->SetBranchAddress( "SigmaIEtaIEta", &fEleSigmaIEtaIEta); 
  RealEleSourceTree->SetBranchAddress( "DEtaIn", &fEleDEtaIn); 
  RealEleSourceTree->SetBranchAddress( "DPhiIn", &fEleDPhiIn); 
  RealEleSourceTree->SetBranchAddress( "HoverE", &fEleHoverE); 
  RealEleSourceTree->SetBranchAddress( "D0", &fEleD0); 
  RealEleSourceTree->SetBranchAddress( "DZ", &fEleDZ); 
  RealEleSourceTree->SetBranchAddress( "FBrem", &fEleFBrem); 
  RealEleSourceTree->SetBranchAddress( "EOverP", &fEleEOverP); 
  RealEleSourceTree->SetBranchAddress( "ESeedClusterOverPout", &fEleESeedClusterOverPout); 
  RealEleSourceTree->SetBranchAddress( "SigmaIPhiIPhi", &fEleSigmaIPhiIPhi); 
  RealEleSourceTree->SetBranchAddress( "NBrem", &fEleNBrem); 
  RealEleSourceTree->SetBranchAddress( "OneOverEMinusOneOverP", &fEleOneOverEMinusOneOverP); 
  RealEleSourceTree->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
  RealEleSourceTree->SetBranchAddress( "IP3d", &fEleIP3d); 
  RealEleSourceTree->SetBranchAddress( "IP3dSig", &fEleIP3dSig); 
  RealEleSourceTree->SetBranchAddress( "HcalDepth1OverEcal", &fEleHcalDepth1OverEcal); 
  RealEleSourceTree->SetBranchAddress( "HcalDepth2OverEcal", &fEleHcalDepth2OverEcal); 
  RealEleSourceTree->SetBranchAddress( "dEtaCalo", &fEledEtaCalo); 
  RealEleSourceTree->SetBranchAddress( "dPhiCalo", &fEledPhiCalo); 
  RealEleSourceTree->SetBranchAddress( "PreShowerOverRaw", &fElePreShowerOverRaw); 
  RealEleSourceTree->SetBranchAddress( "CovIEtaIPhi", &fEleCovIEtaIPhi); 
  RealEleSourceTree->SetBranchAddress( "SCEtaWidth", &fEleSCEtaWidth); 
  RealEleSourceTree->SetBranchAddress( "SCPhiWidth", &fEleSCPhiWidth); 
  RealEleSourceTree->SetBranchAddress( "GsfTrackChi2OverNdof", &fEleGsfTrackChi2OverNdof); 
  RealEleSourceTree->SetBranchAddress( "R9", &fEleR9); 
  RealEleSourceTree->SetBranchAddress( "SeedEMaxOverE", &fEleSeedEMaxOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedETopOverE", &fEleSeedETopOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedEBottomOverE", &fEleSeedEBottomOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedELeftOverE", &fEleSeedELeftOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedERightOverE", &fEleSeedERightOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE2ndOverE", &fEleSeedE2ndOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE2x5RightOverE", &fEleSeedE2x5RightOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE2x5LeftOverE", &fEleSeedE2x5LeftOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE2x5TopOverE", &fEleSeedE2x5TopOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE2x5BottomOverE", &fEleSeedE2x5BottomOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE2x5MaxOverE", &fEleSeedE2x5MaxOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE1x3OverE", &fEleSeedE1x3OverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE3x1OverE", &fEleSeedE3x1OverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE1x5OverE", &fEleSeedE1x5OverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE2x2OverE", &fEleSeedE2x2OverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE3x2OverE", &fEleSeedE3x2OverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE3x3OverE", &fEleSeedE3x3OverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE4x4OverE", &fEleSeedE4x4OverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE5x5OverE", &fEleSeedE5x5OverE); 
  RealEleSourceTree->SetBranchAddress( "ChargedIso03", &fEleChargedIso03); 
  RealEleSourceTree->SetBranchAddress( "NeutralHadronIso03", &fEleNeutralHadronIso03); 
  RealEleSourceTree->SetBranchAddress( "GammaIso03", &fEleGammaIso03); 
  RealEleSourceTree->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
  RealEleSourceTree->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
  RealEleSourceTree->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
  RealEleSourceTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  RealEleSourceTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  RealEleSourceTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  RealEleSourceTree->SetBranchAddress( "TrkIso03", &fEleTrkIso03); 
  RealEleSourceTree->SetBranchAddress( "EMIso03", &fEleEMIso03); 
  RealEleSourceTree->SetBranchAddress( "HadIso03", &fEleHadIso03); 
  RealEleSourceTree->SetBranchAddress( "TrkIso04", &fEleTrkIso04); 
  RealEleSourceTree->SetBranchAddress( "EMIso04", &fEleEMIso04); 
  RealEleSourceTree->SetBranchAddress( "HadIso04", &fEleHadIso04); 
  RealEleSourceTree->SetBranchAddress( "Rho", &fRho); 
  RealEleSourceTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < RealEleSourceTree->GetEntries(); ientry++) {       	
    RealEleSourceTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
    
    histRealElePtSource->Fill(fElePt);
  } 
  


  //*****************************************************************************************
  //RealEleTargetTree
  //*****************************************************************************************
  TFile *RealEleTargetFile = new TFile("ElectronSelectionTraining.HWW115.root", "READ");
  TTree *RealEleTargetTree = (TTree*)RealEleTargetFile->Get("Electrons");
  RealEleTargetTree->SetBranchAddress( "weight", &fWeight);
  RealEleTargetTree->SetBranchAddress( "run", &fRunNumber);
  RealEleTargetTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  RealEleTargetTree->SetBranchAddress( "event", &fEventNumber);
  RealEleTargetTree->SetBranchAddress( "pt", &fElePt); 
  RealEleTargetTree->SetBranchAddress( "eta", &fEleEta); 
  RealEleTargetTree->SetBranchAddress( "phi", &fElePhi); 
  RealEleTargetTree->SetBranchAddress( "scet", &fEleSCEt); 
  RealEleTargetTree->SetBranchAddress( "sceta", &fEleSCEta); 
  RealEleTargetTree->SetBranchAddress( "scphi", &fEleSCPhi); 
  RealEleTargetTree->SetBranchAddress( "pfiso", &fElePFIso); 
  RealEleTargetTree->SetBranchAddress( "SigmaIEtaIEta", &fEleSigmaIEtaIEta); 
  RealEleTargetTree->SetBranchAddress( "DEtaIn", &fEleDEtaIn); 
  RealEleTargetTree->SetBranchAddress( "DPhiIn", &fEleDPhiIn); 
  RealEleTargetTree->SetBranchAddress( "HoverE", &fEleHoverE); 
  RealEleTargetTree->SetBranchAddress( "D0", &fEleD0); 
  RealEleTargetTree->SetBranchAddress( "DZ", &fEleDZ); 
  RealEleTargetTree->SetBranchAddress( "FBrem", &fEleFBrem); 
  RealEleTargetTree->SetBranchAddress( "EOverP", &fEleEOverP); 
  RealEleTargetTree->SetBranchAddress( "ESeedClusterOverPout", &fEleESeedClusterOverPout); 
  RealEleTargetTree->SetBranchAddress( "SigmaIPhiIPhi", &fEleSigmaIPhiIPhi); 
  RealEleTargetTree->SetBranchAddress( "NBrem", &fEleNBrem); 
  RealEleTargetTree->SetBranchAddress( "OneOverEMinusOneOverP", &fEleOneOverEMinusOneOverP); 
  RealEleTargetTree->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
  RealEleTargetTree->SetBranchAddress( "IP3d", &fEleIP3d); 
  RealEleTargetTree->SetBranchAddress( "IP3dSig", &fEleIP3dSig); 
  RealEleTargetTree->SetBranchAddress( "HcalDepth1OverEcal", &fEleHcalDepth1OverEcal); 
  RealEleTargetTree->SetBranchAddress( "HcalDepth2OverEcal", &fEleHcalDepth2OverEcal); 
  RealEleTargetTree->SetBranchAddress( "dEtaCalo", &fEledEtaCalo); 
  RealEleTargetTree->SetBranchAddress( "dPhiCalo", &fEledPhiCalo); 
  RealEleTargetTree->SetBranchAddress( "PreShowerOverRaw", &fElePreShowerOverRaw); 
  RealEleTargetTree->SetBranchAddress( "CovIEtaIPhi", &fEleCovIEtaIPhi); 
  RealEleTargetTree->SetBranchAddress( "SCEtaWidth", &fEleSCEtaWidth); 
  RealEleTargetTree->SetBranchAddress( "SCPhiWidth", &fEleSCPhiWidth); 
  RealEleTargetTree->SetBranchAddress( "GsfTrackChi2OverNdof", &fEleGsfTrackChi2OverNdof); 
  RealEleTargetTree->SetBranchAddress( "R9", &fEleR9); 
  RealEleTargetTree->SetBranchAddress( "SeedEMaxOverE", &fEleSeedEMaxOverE); 
  RealEleTargetTree->SetBranchAddress( "SeedETopOverE", &fEleSeedETopOverE); 
  RealEleTargetTree->SetBranchAddress( "SeedEBottomOverE", &fEleSeedEBottomOverE); 
  RealEleTargetTree->SetBranchAddress( "SeedELeftOverE", &fEleSeedELeftOverE); 
  RealEleTargetTree->SetBranchAddress( "SeedERightOverE", &fEleSeedERightOverE); 
  RealEleTargetTree->SetBranchAddress( "SeedE2ndOverE", &fEleSeedE2ndOverE); 
  RealEleTargetTree->SetBranchAddress( "SeedE2x5RightOverE", &fEleSeedE2x5RightOverE); 
  RealEleTargetTree->SetBranchAddress( "SeedE2x5LeftOverE", &fEleSeedE2x5LeftOverE); 
  RealEleTargetTree->SetBranchAddress( "SeedE2x5TopOverE", &fEleSeedE2x5TopOverE); 
  RealEleTargetTree->SetBranchAddress( "SeedE2x5BottomOverE", &fEleSeedE2x5BottomOverE); 
  RealEleTargetTree->SetBranchAddress( "SeedE2x5MaxOverE", &fEleSeedE2x5MaxOverE); 
  RealEleTargetTree->SetBranchAddress( "SeedE1x3OverE", &fEleSeedE1x3OverE); 
  RealEleTargetTree->SetBranchAddress( "SeedE3x1OverE", &fEleSeedE3x1OverE); 
  RealEleTargetTree->SetBranchAddress( "SeedE1x5OverE", &fEleSeedE1x5OverE); 
  RealEleTargetTree->SetBranchAddress( "SeedE2x2OverE", &fEleSeedE2x2OverE); 
  RealEleTargetTree->SetBranchAddress( "SeedE3x2OverE", &fEleSeedE3x2OverE); 
  RealEleTargetTree->SetBranchAddress( "SeedE3x3OverE", &fEleSeedE3x3OverE); 
  RealEleTargetTree->SetBranchAddress( "SeedE4x4OverE", &fEleSeedE4x4OverE); 
  RealEleTargetTree->SetBranchAddress( "SeedE5x5OverE", &fEleSeedE5x5OverE); 
  RealEleTargetTree->SetBranchAddress( "ChargedIso03", &fEleChargedIso03); 
  RealEleTargetTree->SetBranchAddress( "NeutralHadronIso03", &fEleNeutralHadronIso03); 
  RealEleTargetTree->SetBranchAddress( "GammaIso03", &fEleGammaIso03); 
  RealEleTargetTree->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
  RealEleTargetTree->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
  RealEleTargetTree->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
  RealEleTargetTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  RealEleTargetTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  RealEleTargetTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  RealEleTargetTree->SetBranchAddress( "TrkIso03", &fEleTrkIso03); 
  RealEleTargetTree->SetBranchAddress( "EMIso03", &fEleEMIso03); 
  RealEleTargetTree->SetBranchAddress( "HadIso03", &fEleHadIso03); 
  RealEleTargetTree->SetBranchAddress( "TrkIso04", &fEleTrkIso04); 
  RealEleTargetTree->SetBranchAddress( "EMIso04", &fEleEMIso04); 
  RealEleTargetTree->SetBranchAddress( "HadIso04", &fEleHadIso04); 
  RealEleTargetTree->SetBranchAddress( "Rho", &fRho); 
  RealEleTargetTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < RealEleTargetTree->GetEntries(); ientry++) {       	
    RealEleTargetTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
    
    histRealElePtTarget->Fill(fElePt);
  } //loop over electrons
  


  //*****************************************************************************************
  //FakeEleTargetTree
  //*****************************************************************************************
  TFile *FakeEleSourceFile = new TFile("ElectronSelectionTraining.Fakes.root", "READ");
  TTree *FakeEleSourceTree = (TTree*)FakeEleSourceFile->Get("Electrons");
  FakeEleSourceTree->SetBranchAddress( "weight", &fWeight);
  FakeEleSourceTree->SetBranchAddress( "run", &fRunNumber);
  FakeEleSourceTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  FakeEleSourceTree->SetBranchAddress( "event", &fEventNumber);
  FakeEleSourceTree->SetBranchAddress( "pt", &fElePt); 
  FakeEleSourceTree->SetBranchAddress( "eta", &fEleEta); 
  FakeEleSourceTree->SetBranchAddress( "phi", &fElePhi); 
  FakeEleSourceTree->SetBranchAddress( "scet", &fEleSCEt); 
  FakeEleSourceTree->SetBranchAddress( "sceta", &fEleSCEta); 
  FakeEleSourceTree->SetBranchAddress( "scphi", &fEleSCPhi); 
  FakeEleSourceTree->SetBranchAddress( "pfiso", &fElePFIso); 
  FakeEleSourceTree->SetBranchAddress( "SigmaIEtaIEta", &fEleSigmaIEtaIEta); 
  FakeEleSourceTree->SetBranchAddress( "DEtaIn", &fEleDEtaIn); 
  FakeEleSourceTree->SetBranchAddress( "DPhiIn", &fEleDPhiIn); 
  FakeEleSourceTree->SetBranchAddress( "HoverE", &fEleHoverE); 
  FakeEleSourceTree->SetBranchAddress( "D0", &fEleD0); 
  FakeEleSourceTree->SetBranchAddress( "DZ", &fEleDZ); 
  FakeEleSourceTree->SetBranchAddress( "FBrem", &fEleFBrem); 
  FakeEleSourceTree->SetBranchAddress( "EOverP", &fEleEOverP); 
  FakeEleSourceTree->SetBranchAddress( "ESeedClusterOverPout", &fEleESeedClusterOverPout); 
  FakeEleSourceTree->SetBranchAddress( "SigmaIPhiIPhi", &fEleSigmaIPhiIPhi); 
  FakeEleSourceTree->SetBranchAddress( "NBrem", &fEleNBrem); 
  FakeEleSourceTree->SetBranchAddress( "OneOverEMinusOneOverP", &fEleOneOverEMinusOneOverP); 
  FakeEleSourceTree->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
  FakeEleSourceTree->SetBranchAddress( "IP3d", &fEleIP3d); 
  FakeEleSourceTree->SetBranchAddress( "IP3dSig", &fEleIP3dSig); 
  FakeEleSourceTree->SetBranchAddress( "HcalDepth1OverEcal", &fEleHcalDepth1OverEcal); 
  FakeEleSourceTree->SetBranchAddress( "HcalDepth2OverEcal", &fEleHcalDepth2OverEcal); 
  FakeEleSourceTree->SetBranchAddress( "dEtaCalo", &fEledEtaCalo); 
  FakeEleSourceTree->SetBranchAddress( "dPhiCalo", &fEledPhiCalo); 
  FakeEleSourceTree->SetBranchAddress( "PreShowerOverRaw", &fElePreShowerOverRaw); 
  FakeEleSourceTree->SetBranchAddress( "CovIEtaIPhi", &fEleCovIEtaIPhi); 
  FakeEleSourceTree->SetBranchAddress( "SCEtaWidth", &fEleSCEtaWidth); 
  FakeEleSourceTree->SetBranchAddress( "SCPhiWidth", &fEleSCPhiWidth); 
  FakeEleSourceTree->SetBranchAddress( "GsfTrackChi2OverNdof", &fEleGsfTrackChi2OverNdof); 
  FakeEleSourceTree->SetBranchAddress( "R9", &fEleR9); 
  FakeEleSourceTree->SetBranchAddress( "SeedEMaxOverE", &fEleSeedEMaxOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedETopOverE", &fEleSeedETopOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedEBottomOverE", &fEleSeedEBottomOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedELeftOverE", &fEleSeedELeftOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedERightOverE", &fEleSeedERightOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2ndOverE", &fEleSeedE2ndOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x5RightOverE", &fEleSeedE2x5RightOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x5LeftOverE", &fEleSeedE2x5LeftOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x5TopOverE", &fEleSeedE2x5TopOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x5BottomOverE", &fEleSeedE2x5BottomOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x5MaxOverE", &fEleSeedE2x5MaxOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE1x3OverE", &fEleSeedE1x3OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE3x1OverE", &fEleSeedE3x1OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE1x5OverE", &fEleSeedE1x5OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x2OverE", &fEleSeedE2x2OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE3x2OverE", &fEleSeedE3x2OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE3x3OverE", &fEleSeedE3x3OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE4x4OverE", &fEleSeedE4x4OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE5x5OverE", &fEleSeedE5x5OverE); 
  FakeEleSourceTree->SetBranchAddress( "ChargedIso03", &fEleChargedIso03); 
  FakeEleSourceTree->SetBranchAddress( "NeutralHadronIso03", &fEleNeutralHadronIso03); 
  FakeEleSourceTree->SetBranchAddress( "GammaIso03", &fEleGammaIso03); 
  FakeEleSourceTree->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
  FakeEleSourceTree->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
  FakeEleSourceTree->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
  FakeEleSourceTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  FakeEleSourceTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  FakeEleSourceTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  FakeEleSourceTree->SetBranchAddress( "TrkIso03", &fEleTrkIso03); 
  FakeEleSourceTree->SetBranchAddress( "EMIso03", &fEleEMIso03); 
  FakeEleSourceTree->SetBranchAddress( "HadIso03", &fEleHadIso03); 
  FakeEleSourceTree->SetBranchAddress( "TrkIso04", &fEleTrkIso04); 
  FakeEleSourceTree->SetBranchAddress( "EMIso04", &fEleEMIso04); 
  FakeEleSourceTree->SetBranchAddress( "HadIso04", &fEleHadIso04); 
  FakeEleSourceTree->SetBranchAddress( "Rho", &fRho); 
  FakeEleSourceTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < FakeEleSourceTree->GetEntries(); ientry++) {       	
    FakeEleSourceTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
    
    histFakeElePtSource->Fill(fElePt);
  } //loop over electrons
  

  TFile *file = new TFile("ElectronPtSpectrum.root", "UPDATE");
  file->cd();
  file->WriteTObject(histRealElePtSource, histRealElePtSource->GetName(), "WriteDelete");  
  file->WriteTObject(histRealElePtTarget, histRealElePtTarget->GetName(), "WriteDelete");  
  file->WriteTObject(histFakeElePtSource, histFakeElePtSource->GetName(), "WriteDelete");  
  file->Close();
  delete file;



  gBenchmark->Show("WWTemplate");       
} 

//*************************************************************************************************
//Compute Reweight Factors
//*************************************************************************************************
void  MakeReweightFactors() {
  
  TFile* file = new TFile("ElectronPtSpectrum.root","READ");
  TH1F *histRealElePtSource = (TH1F*)file->Get("RealElePtSource");
  TH1F *histRealElePtTarget = (TH1F*)file->Get("RealElePtTarget");
  TH1F *histFakeElePtSource = (TH1F*)file->Get("FakeElePtSource");
  TH1F *histFakeElePtTarget = (TH1F*)file->Get("FakeElePtTarget");

  NormalizeHist(histRealElePtSource);
  NormalizeHist(histRealElePtTarget);
  NormalizeHist(histFakeElePtSource);
  NormalizeHist(histFakeElePtTarget);


  TH1F *RealElectronPtReweightFactor = (TH1F*)histRealElePtSource->Clone("RealElectronPtReweightFactor");
  RealElectronPtReweightFactor->SetBinContent(0,1.0);
  for(UInt_t a=1; a < RealElectronPtReweightFactor->GetXaxis()->GetNbins()+2; ++a) {
    if (histRealElePtSource->GetBinContent(a)>0) {
      RealElectronPtReweightFactor->SetBinContent(a,histRealElePtTarget->GetBinContent(a) / histRealElePtSource->GetBinContent(a));
      cout << a << " " << histRealElePtTarget->GetBinContent(a) << " / " << histRealElePtSource->GetBinContent(a) << " = " << histRealElePtTarget->GetBinContent(a) / histRealElePtSource->GetBinContent(a) << endl;
    } else {
      RealElectronPtReweightFactor->SetBinContent(a,1.0);
    }
  }

  TH1F *FakeElectronPtReweightFactor = (TH1F*)histFakeElePtSource->Clone("FakeElectronPtReweightFactor");
  FakeElectronPtReweightFactor->SetBinContent(0,1.0);
  for(UInt_t a=1; a < FakeElectronPtReweightFactor->GetXaxis()->GetNbins()+2; ++a) {
    if (histFakeElePtSource->GetBinContent(a)>0) {
      FakeElectronPtReweightFactor->SetBinContent(a,histFakeElePtTarget->GetBinContent(a) / histFakeElePtSource->GetBinContent(a));
    } else {
      FakeElectronPtReweightFactor->SetBinContent(a,1.0);
    }
  }

  TFile *outfile = new TFile("ElectronPtSpectrum.root", "UPDATE");
  outfile->cd();
  outfile->WriteTObject(RealElectronPtReweightFactor, RealElectronPtReweightFactor->GetName(), "WriteDelete");  
  outfile->WriteTObject(FakeElectronPtReweightFactor, FakeElectronPtReweightFactor->GetName(), "WriteDelete");  
  outfile->Close();
  delete outfile;

  TCanvas *cv = new TCanvas("cv","cv", 800, 600);
  histRealElePtSource->GetXaxis()->SetTitleOffset(1.05);
  histRealElePtSource->GetYaxis()->SetTitleOffset(1.4);
  histRealElePtSource->Draw("hist");
  cv->SaveAs("SignalPt_Source.png");

  histRealElePtTarget->GetXaxis()->SetTitleOffset(1.05);
  histRealElePtTarget->GetYaxis()->SetTitleOffset(1.4);
  histRealElePtTarget->Draw("hist");
  cv->SaveAs("SignalPt_Target.png");

  histFakeElePtSource->GetXaxis()->SetTitleOffset(1.05);
  histFakeElePtSource->GetYaxis()->SetTitleOffset(1.4);
  histFakeElePtSource->Draw("hist");
  cv->SaveAs("BkgPt_Source.png");

  histFakeElePtTarget->GetXaxis()->SetTitleOffset(1.05);
  histFakeElePtTarget->GetYaxis()->SetTitleOffset(1.4);
  histFakeElePtTarget->Draw("hist");
  cv->SaveAs("BkgPt_Target.png");


}


//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void UpdateNtupleWeights()
{  


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TFile* inputfile = new TFile("ElectronPtSpectrum.root","READ");
  TH1F *RealElectronPtReweightFactor = (TH1F*)inputfile->Get("RealElectronPtReweightFactor");
  RealElectronPtReweightFactor->SetDirectory(0);
  TH1F *FakeElectronPtReweightFactor = (TH1F*)inputfile->Get("FakeElectronPtReweightFactor");
  FakeElectronPtReweightFactor->SetDirectory(0);
  inputfile->Close();
  delete inputfile;


  //*****************************************************************************************
  //Setup
  //*****************************************************************************************

  //Variables
  Float_t                 fWeight;
  UInt_t                  fRunNumber;
  UInt_t                  fLumiSectionNumber;
  UInt_t                  fEventNumber;
  Float_t                 fElePt; 
  Float_t                 fEleEta; 
  Float_t                 fElePhi; 
  Float_t                 fEleSCEt; 
  Float_t                 fEleSCEta; 
  Float_t                 fEleSCPhi; 
  Float_t                 fElePFIso; 
  
  //CutBased Variables
  Float_t                 fEleSigmaIEtaIEta; 
  Float_t                 fEleDEtaIn; 
  Float_t                 fEleDPhiIn; 
  Float_t                 fEleHoverE; 
  Float_t                 fEleD0; 
  Float_t                 fEleDZ; 
  Float_t                 fEleFBrem; 
  Float_t                 fEleEOverP; 

  //Additional Vars used in Likelihood
  Float_t                 fEleESeedClusterOverPout; 
  Float_t                 fEleSigmaIPhiIPhi; 
  Float_t                 fEleNBrem; 
  Float_t                 fEleOneOverEMinusOneOverP; 
  Float_t                 fEleESeedClusterOverPIn; 
  Float_t                 fEleIP3d; 
  Float_t                 fEleIP3dSig; 

  Float_t                 fEleHcalDepth1OverEcal;
  Float_t                 fEleHcalDepth2OverEcal;
  Float_t                 fEledEtaCalo;
  Float_t                 fEledPhiCalo;
  Float_t                 fElePreShowerOverRaw;
  Float_t                 fEleCovIEtaIPhi;
  Float_t                 fEleSCEtaWidth;
  Float_t                 fEleSCPhiWidth;
  Float_t                 fEleGsfTrackChi2OverNdof;
  Float_t                 fEleR9;

  Float_t                 fEleSeedEMaxOverE;
  Float_t                 fEleSeedETopOverE;
  Float_t                 fEleSeedEBottomOverE;
  Float_t                 fEleSeedELeftOverE;
  Float_t                 fEleSeedERightOverE;
  Float_t                 fEleSeedE2ndOverE;
  Float_t                 fEleSeedE2x5RightOverE;
  Float_t                 fEleSeedE2x5LeftOverE;
  Float_t                 fEleSeedE2x5TopOverE;
  Float_t                 fEleSeedE2x5BottomOverE;
  Float_t                 fEleSeedE2x5MaxOverE;
  Float_t                 fEleSeedE1x3OverE;
  Float_t                 fEleSeedE3x1OverE;
  Float_t                 fEleSeedE1x5OverE;
  Float_t                 fEleSeedE2x2OverE;
  Float_t                 fEleSeedE3x2OverE;
  Float_t                 fEleSeedE3x3OverE;
  Float_t                 fEleSeedE4x4OverE;
  Float_t                 fEleSeedE5x5OverE;


  //Isolation Variables
  Float_t                 fEleStandardLikelihood; 
  Float_t                 fElePFMVA; 
  Float_t                 fEleChargedIso03; 
  Float_t                 fEleNeutralHadronIso03; 
  Float_t                 fEleGammaIso03; 
  Float_t                 fEleChargedIso04; 
  Float_t                 fEleNeutralHadronIso04; 
  Float_t                 fEleGammaIso04; 
  Float_t                 fEleChargedIso04FromOtherVertices; 
  Float_t                 fEleNeutralHadronIso04_10Threshold; 
  Float_t                 fEleGammaIso04_10Threshold; 
  Float_t                 fEleTrkIso03; 
  Float_t                 fEleEMIso03; 
  Float_t                 fEleHadIso03; 
  Float_t                 fEleTrkIso04; 
  Float_t                 fEleEMIso04; 
  Float_t                 fEleHadIso04; 
  Float_t                 fRho; 
  Float_t                 fNVertices; 


  //*****************************************************************************************
  //RealEleSourceTree Output
  //*****************************************************************************************
  TFile *RealEleSourceOutputFile = new TFile("ElectronSelectionTraining.Real.weighted.root", "RECREATE");
  TTree *RealEleSourceOutputTree = new TTree("Electrons","Electrons");
  RealEleSourceOutputTree->SetAutoFlush(0);

  RealEleSourceOutputTree->Branch("weight",&fWeight,"weight/F");
  RealEleSourceOutputTree->Branch("run",&fRunNumber,"run/i");
  RealEleSourceOutputTree->Branch("lumi",&fLumiSectionNumber,"lumi/i");
  RealEleSourceOutputTree->Branch("event",&fEventNumber,"event/i");
  RealEleSourceOutputTree->Branch("pt",&fElePt,"pt/F"); 
  RealEleSourceOutputTree->Branch("eta",&fEleEta,"eta/F"); 
  RealEleSourceOutputTree->Branch("phi",&fElePhi,"phi/F"); 
  RealEleSourceOutputTree->Branch("scet",&fEleSCEt,"scet/F"); 
  RealEleSourceOutputTree->Branch("sceta",&fEleSCEta,"sceta/F"); 
  RealEleSourceOutputTree->Branch("scphi",&fEleSCPhi,"scphi/F"); 
  RealEleSourceOutputTree->Branch("pfiso",&fElePFIso,"pfiso/F"); 
  
  //CutBased Variables
  RealEleSourceOutputTree->Branch("SigmaIEtaIEta",&fEleSigmaIEtaIEta,"SigmaIEtaIEta/F"); 
  RealEleSourceOutputTree->Branch("DEtaIn",&fEleDEtaIn,"DEtaIn/F"); 
  RealEleSourceOutputTree->Branch("DPhiIn",&fEleDPhiIn,"DPhiIn/F"); 
  RealEleSourceOutputTree->Branch("HoverE",&fEleHoverE,"HoverE/F"); 
  RealEleSourceOutputTree->Branch("D0",&fEleD0,"D0/F"); 
  RealEleSourceOutputTree->Branch("DZ",&fEleDZ,"DZ/F"); 
  RealEleSourceOutputTree->Branch("FBrem",&fEleFBrem,"FBrem/F"); 
  RealEleSourceOutputTree->Branch("EOverP",&fEleEOverP,"EOverP/F"); 

  //Additional Vars used in Likelihood
  RealEleSourceOutputTree->Branch("ESeedClusterOverPout",&fEleESeedClusterOverPout,"ESeedClusterOverPout/F"); 
  RealEleSourceOutputTree->Branch("SigmaIPhiIPhi",&fEleSigmaIPhiIPhi,"SigmaIPhiIPhi/F"); 
  RealEleSourceOutputTree->Branch("NBrem",&fEleNBrem,"NBrem/F"); 
  RealEleSourceOutputTree->Branch("OneOverEMinusOneOverP",&fEleOneOverEMinusOneOverP,"OneOverEMinusOneOverP/F"); 
  RealEleSourceOutputTree->Branch("ESeedClusterOverPIn",&fEleESeedClusterOverPIn,"ESeedClusterOverPIn/F"); 
  RealEleSourceOutputTree->Branch("IP3d",&fEleIP3d,"IP3d/F"); 
  RealEleSourceOutputTree->Branch("IP3dSig",&fEleIP3dSig,"IP3dSig/F"); 

  RealEleSourceOutputTree->Branch("HcalDepth1OverEcal",&fEleHcalDepth1OverEcal,"HcalDepth1OverEcal/F"); 
  RealEleSourceOutputTree->Branch("HcalDepth2OverEcal",&fEleHcalDepth2OverEcal,"HcalDepth2OverEcal/F"); 
  RealEleSourceOutputTree->Branch("dEtaCalo",&fEledEtaCalo,"dEtaCalo/F"); 
  RealEleSourceOutputTree->Branch("dPhiCalo",&fEledPhiCalo,"dPhiCalo/F"); 
  RealEleSourceOutputTree->Branch("PreShowerOverRaw",&fElePreShowerOverRaw,"PreShowerOverRaw/F"); 
  RealEleSourceOutputTree->Branch("CovIEtaIPhi",&fEleCovIEtaIPhi,"CovIEtaIPhi/F"); 
  RealEleSourceOutputTree->Branch("SCEtaWidth",&fEleSCEtaWidth,"SCEtaWidth/F"); 
  RealEleSourceOutputTree->Branch("SCPhiWidth",&fEleSCPhiWidth,"SCPhiWidth/F"); 
  RealEleSourceOutputTree->Branch("GsfTrackChi2OverNdof",&fEleGsfTrackChi2OverNdof,"GsfTrackChi2OverNdof/F"); 
  RealEleSourceOutputTree->Branch("R9",&fEleR9,"R9/F"); 

  RealEleSourceOutputTree->Branch("SeedEMaxOverE",&fEleSeedEMaxOverE,"SeedEMaxOverE/F"); 
  RealEleSourceOutputTree->Branch("SeedETopOverE",&fEleSeedETopOverE,"SeedETopOverE/F"); 
  RealEleSourceOutputTree->Branch("SeedEBottomOverE",&fEleSeedEBottomOverE,"SeedEBottomOverE/F"); 
  RealEleSourceOutputTree->Branch("SeedELeftOverE",&fEleSeedELeftOverE,"SeedELeftOverE/F"); 
  RealEleSourceOutputTree->Branch("SeedERightOverE",&fEleSeedERightOverE,"SeedERightOverE/F"); 
  RealEleSourceOutputTree->Branch("SeedE2ndOverE",&fEleSeedE2ndOverE,"SeedE2ndOverE/F"); 
  RealEleSourceOutputTree->Branch("SeedE2x5RightOverE",&fEleSeedE2x5RightOverE,"SeedE2x5RightOverE/F"); 
  RealEleSourceOutputTree->Branch("SeedE2x5LeftOverE",&fEleSeedE2x5LeftOverE,"SeedE2x5LeftOverE/F"); 
  RealEleSourceOutputTree->Branch("SeedE2x5TopOverE",&fEleSeedE2x5TopOverE,"SeedE2x5TopOverE/F"); 
  RealEleSourceOutputTree->Branch("SeedE2x5BottomOverE",&fEleSeedE2x5BottomOverE,"SeedE2x5BottomOverE/F"); 
  RealEleSourceOutputTree->Branch("SeedE2x5MaxOverE",&fEleSeedE2x5MaxOverE,"SeedE2x5MaxOverE/F"); 
  RealEleSourceOutputTree->Branch("SeedE1x3OverE",&fEleSeedE1x3OverE,"SeedE1x3OverE/F"); 
  RealEleSourceOutputTree->Branch("SeedE3x1OverE",&fEleSeedE3x1OverE,"SeedE3x1OverE/F"); 
  RealEleSourceOutputTree->Branch("SeedE1x5OverE",&fEleSeedE1x5OverE,"SeedE1x5OverE/F"); 
  RealEleSourceOutputTree->Branch("SeedE2x2OverE",&fEleSeedE2x2OverE,"SeedE2x2OverE/F"); 
  RealEleSourceOutputTree->Branch("SeedE3x2OverE",&fEleSeedE3x2OverE,"SeedE3x2OverE/F"); 
  RealEleSourceOutputTree->Branch("SeedE3x3OverE",&fEleSeedE3x3OverE,"SeedE3x3OverE/F"); 
  RealEleSourceOutputTree->Branch("SeedE4x4OverE",&fEleSeedE4x4OverE,"SeedE4x4OverE/F"); 
  RealEleSourceOutputTree->Branch("SeedE5x5OverE",&fEleSeedE5x5OverE,"SeedE5x5OverE/F"); 


  //Isolation Variables
  RealEleSourceOutputTree->Branch("StandardLikelihood",&fEleStandardLikelihood,"StandardLikelihood/F"); 
  RealEleSourceOutputTree->Branch("PFMVA",&fElePFMVA,"PFMVA/F"); 
  RealEleSourceOutputTree->Branch("ChargedIso03",&fEleChargedIso03,"ChargedIso03/F"); 
  RealEleSourceOutputTree->Branch("NeutralHadronIso03",&fEleNeutralHadronIso03,"NeutralHadronIso03/F"); 
  RealEleSourceOutputTree->Branch("GammaIso03",&fEleGammaIso03,"GammaIso03/F"); 
  RealEleSourceOutputTree->Branch("ChargedIso04",&fEleChargedIso04,"ChargedIso04/F"); 
  RealEleSourceOutputTree->Branch("NeutralHadronIso04",&fEleNeutralHadronIso04,"NeutralHadronIso04/F"); 
  RealEleSourceOutputTree->Branch("GammaIso04",&fEleGammaIso04,"GammaIso04/F"); 
  RealEleSourceOutputTree->Branch("ChargedIso04FromOtherVertices",&fEleChargedIso04FromOtherVertices,"ChargedIso04FromOtherVertices/F"); 
  RealEleSourceOutputTree->Branch("NeutralHadronIso04_10Threshold",&fEleNeutralHadronIso04_10Threshold,"NeutralHadronIso04_10Threshold/F"); 
  RealEleSourceOutputTree->Branch("GammaIso04_10Threshold",&fEleGammaIso04_10Threshold,"GammaIso04_10Threshold/F"); 
  RealEleSourceOutputTree->Branch("TrkIso03",&fEleTrkIso03,"TrkIso03/F"); 
  RealEleSourceOutputTree->Branch("EMIso03",&fEleEMIso03,"EMIso03/F"); 
  RealEleSourceOutputTree->Branch("HadIso03",&fEleHadIso03,"HadIso03/F"); 
  RealEleSourceOutputTree->Branch("TrkIso04",&fEleTrkIso04,"TrkIso04/F"); 
  RealEleSourceOutputTree->Branch("EMIso04",&fEleEMIso04,"EMIso04/F"); 
  RealEleSourceOutputTree->Branch("HadIso04",&fEleHadIso04,"HadIso04/F"); 
  RealEleSourceOutputTree->Branch("Rho",&fRho,"Rho/F"); 
  RealEleSourceOutputTree->Branch("NVertices",&fNVertices,"NVertices/F"); 

  
  //*****************************************************************************************
  //RealEleSourceTree
  //*****************************************************************************************
  TFile *RealEleSourceFile = new TFile("ElectronSelectionTraining.Real.root", "READ");
  TTree *RealEleSourceTree = (TTree*)RealEleSourceFile->Get("Electrons");
  RealEleSourceTree->SetBranchAddress( "weight", &fWeight);
  RealEleSourceTree->SetBranchAddress( "run", &fRunNumber);
  RealEleSourceTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  RealEleSourceTree->SetBranchAddress( "event", &fEventNumber);
  RealEleSourceTree->SetBranchAddress( "pt", &fElePt); 
  RealEleSourceTree->SetBranchAddress( "eta", &fEleEta); 
  RealEleSourceTree->SetBranchAddress( "phi", &fElePhi); 
  RealEleSourceTree->SetBranchAddress( "scet", &fEleSCEt); 
  RealEleSourceTree->SetBranchAddress( "sceta", &fEleSCEta); 
  RealEleSourceTree->SetBranchAddress( "scphi", &fEleSCPhi); 
  RealEleSourceTree->SetBranchAddress( "pfiso", &fElePFIso); 
  RealEleSourceTree->SetBranchAddress( "SigmaIEtaIEta", &fEleSigmaIEtaIEta); 
  RealEleSourceTree->SetBranchAddress( "DEtaIn", &fEleDEtaIn); 
  RealEleSourceTree->SetBranchAddress( "DPhiIn", &fEleDPhiIn); 
  RealEleSourceTree->SetBranchAddress( "HoverE", &fEleHoverE); 
  RealEleSourceTree->SetBranchAddress( "D0", &fEleD0); 
  RealEleSourceTree->SetBranchAddress( "DZ", &fEleDZ); 
  RealEleSourceTree->SetBranchAddress( "FBrem", &fEleFBrem); 
  RealEleSourceTree->SetBranchAddress( "EOverP", &fEleEOverP); 
  RealEleSourceTree->SetBranchAddress( "ESeedClusterOverPout", &fEleESeedClusterOverPout); 
  RealEleSourceTree->SetBranchAddress( "SigmaIPhiIPhi", &fEleSigmaIPhiIPhi); 
  RealEleSourceTree->SetBranchAddress( "NBrem", &fEleNBrem); 
  RealEleSourceTree->SetBranchAddress( "OneOverEMinusOneOverP", &fEleOneOverEMinusOneOverP); 
  RealEleSourceTree->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
  RealEleSourceTree->SetBranchAddress( "IP3d", &fEleIP3d); 
  RealEleSourceTree->SetBranchAddress( "IP3dSig", &fEleIP3dSig); 
  RealEleSourceTree->SetBranchAddress( "HcalDepth1OverEcal", &fEleHcalDepth1OverEcal); 
  RealEleSourceTree->SetBranchAddress( "HcalDepth2OverEcal", &fEleHcalDepth2OverEcal); 
  RealEleSourceTree->SetBranchAddress( "dEtaCalo", &fEledEtaCalo); 
  RealEleSourceTree->SetBranchAddress( "dPhiCalo", &fEledPhiCalo); 
  RealEleSourceTree->SetBranchAddress( "PreShowerOverRaw", &fElePreShowerOverRaw); 
  RealEleSourceTree->SetBranchAddress( "CovIEtaIPhi", &fEleCovIEtaIPhi); 
  RealEleSourceTree->SetBranchAddress( "SCEtaWidth", &fEleSCEtaWidth); 
  RealEleSourceTree->SetBranchAddress( "SCPhiWidth", &fEleSCPhiWidth); 
  RealEleSourceTree->SetBranchAddress( "GsfTrackChi2OverNdof", &fEleGsfTrackChi2OverNdof); 
  RealEleSourceTree->SetBranchAddress( "R9", &fEleR9); 
  RealEleSourceTree->SetBranchAddress( "SeedEMaxOverE", &fEleSeedEMaxOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedETopOverE", &fEleSeedETopOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedEBottomOverE", &fEleSeedEBottomOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedELeftOverE", &fEleSeedELeftOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedERightOverE", &fEleSeedERightOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE2ndOverE", &fEleSeedE2ndOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE2x5RightOverE", &fEleSeedE2x5RightOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE2x5LeftOverE", &fEleSeedE2x5LeftOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE2x5TopOverE", &fEleSeedE2x5TopOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE2x5BottomOverE", &fEleSeedE2x5BottomOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE2x5MaxOverE", &fEleSeedE2x5MaxOverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE1x3OverE", &fEleSeedE1x3OverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE3x1OverE", &fEleSeedE3x1OverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE1x5OverE", &fEleSeedE1x5OverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE2x2OverE", &fEleSeedE2x2OverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE3x2OverE", &fEleSeedE3x2OverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE3x3OverE", &fEleSeedE3x3OverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE4x4OverE", &fEleSeedE4x4OverE); 
  RealEleSourceTree->SetBranchAddress( "SeedE5x5OverE", &fEleSeedE5x5OverE); 
  RealEleSourceTree->SetBranchAddress( "StandardLikelihood", &fEleStandardLikelihood); 
  RealEleSourceTree->SetBranchAddress( "PFMVA", &fElePFMVA); 
  RealEleSourceTree->SetBranchAddress( "ChargedIso03", &fEleChargedIso03); 
  RealEleSourceTree->SetBranchAddress( "NeutralHadronIso03", &fEleNeutralHadronIso03); 
  RealEleSourceTree->SetBranchAddress( "GammaIso03", &fEleGammaIso03); 
  RealEleSourceTree->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
  RealEleSourceTree->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
  RealEleSourceTree->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
  RealEleSourceTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  RealEleSourceTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  RealEleSourceTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  RealEleSourceTree->SetBranchAddress( "TrkIso03", &fEleTrkIso03); 
  RealEleSourceTree->SetBranchAddress( "EMIso03", &fEleEMIso03); 
  RealEleSourceTree->SetBranchAddress( "HadIso03", &fEleHadIso03); 
  RealEleSourceTree->SetBranchAddress( "TrkIso04", &fEleTrkIso04); 
  RealEleSourceTree->SetBranchAddress( "EMIso04", &fEleEMIso04); 
  RealEleSourceTree->SetBranchAddress( "HadIso04", &fEleHadIso04); 
  RealEleSourceTree->SetBranchAddress( "Rho", &fRho); 
  RealEleSourceTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < RealEleSourceTree->GetEntries(); ientry++) {       	
    RealEleSourceTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Eelectron " << ientry << endl;
    fWeight = RealElectronPtReweightFactor->GetBinContent(RealElectronPtReweightFactor->GetXaxis()->FindFixBin(fElePt));
    RealEleSourceOutputTree->Fill();
  } 
  RealEleSourceOutputFile->Write();
  RealEleSourceOutputFile->Close();
  delete RealEleSourceFile;

  //*****************************************************************************************
  //FakeEleSourceTree Output
  //*****************************************************************************************
  TFile *FakeEleSourceOutputFile = new TFile("ElectronSelectionTraining.Fake.weighted.root", "RECREATE");
  TTree *FakeEleSourceOutputTree = new TTree("Electrons","Electrons");
  FakeEleSourceOutputTree->SetAutoFlush(0);

  FakeEleSourceOutputTree->Branch("weight",&fWeight,"weight/F");
  FakeEleSourceOutputTree->Branch("run",&fRunNumber,"run/i");
  FakeEleSourceOutputTree->Branch("lumi",&fLumiSectionNumber,"lumi/i");
  FakeEleSourceOutputTree->Branch("event",&fEventNumber,"event/i");
  FakeEleSourceOutputTree->Branch("pt",&fElePt,"pt/F"); 
  FakeEleSourceOutputTree->Branch("eta",&fEleEta,"eta/F"); 
  FakeEleSourceOutputTree->Branch("phi",&fElePhi,"phi/F"); 
  FakeEleSourceOutputTree->Branch("scet",&fEleSCEt,"scet/F"); 
  FakeEleSourceOutputTree->Branch("sceta",&fEleSCEta,"sceta/F"); 
  FakeEleSourceOutputTree->Branch("scphi",&fEleSCPhi,"scphi/F"); 
  FakeEleSourceOutputTree->Branch("pfiso",&fElePFIso,"pfiso/F"); 
  
  //CutBased Variables
  FakeEleSourceOutputTree->Branch("SigmaIEtaIEta",&fEleSigmaIEtaIEta,"SigmaIEtaIEta/F"); 
  FakeEleSourceOutputTree->Branch("DEtaIn",&fEleDEtaIn,"DEtaIn/F"); 
  FakeEleSourceOutputTree->Branch("DPhiIn",&fEleDPhiIn,"DPhiIn/F"); 
  FakeEleSourceOutputTree->Branch("HoverE",&fEleHoverE,"HoverE/F"); 
  FakeEleSourceOutputTree->Branch("D0",&fEleD0,"D0/F"); 
  FakeEleSourceOutputTree->Branch("DZ",&fEleDZ,"DZ/F"); 
  FakeEleSourceOutputTree->Branch("FBrem",&fEleFBrem,"FBrem/F"); 
  FakeEleSourceOutputTree->Branch("EOverP",&fEleEOverP,"EOverP/F"); 

  //Additional Vars used in Likelihood
  FakeEleSourceOutputTree->Branch("ESeedClusterOverPout",&fEleESeedClusterOverPout,"ESeedClusterOverPout/F"); 
  FakeEleSourceOutputTree->Branch("SigmaIPhiIPhi",&fEleSigmaIPhiIPhi,"SigmaIPhiIPhi/F"); 
  FakeEleSourceOutputTree->Branch("NBrem",&fEleNBrem,"NBrem/F"); 
  FakeEleSourceOutputTree->Branch("OneOverEMinusOneOverP",&fEleOneOverEMinusOneOverP,"OneOverEMinusOneOverP/F"); 
  FakeEleSourceOutputTree->Branch("ESeedClusterOverPIn",&fEleESeedClusterOverPIn,"ESeedClusterOverPIn/F"); 
  FakeEleSourceOutputTree->Branch("IP3d",&fEleIP3d,"IP3d/F"); 
  FakeEleSourceOutputTree->Branch("IP3dSig",&fEleIP3dSig,"IP3dSig/F"); 

  FakeEleSourceOutputTree->Branch("HcalDepth1OverEcal",&fEleHcalDepth1OverEcal,"HcalDepth1OverEcal/F"); 
  FakeEleSourceOutputTree->Branch("HcalDepth2OverEcal",&fEleHcalDepth2OverEcal,"HcalDepth2OverEcal/F"); 
  FakeEleSourceOutputTree->Branch("dEtaCalo",&fEledEtaCalo,"dEtaCalo/F"); 
  FakeEleSourceOutputTree->Branch("dPhiCalo",&fEledPhiCalo,"dPhiCalo/F"); 
  FakeEleSourceOutputTree->Branch("PreShowerOverRaw",&fElePreShowerOverRaw,"PreShowerOverRaw/F"); 
  FakeEleSourceOutputTree->Branch("CovIEtaIPhi",&fEleCovIEtaIPhi,"CovIEtaIPhi/F"); 
  FakeEleSourceOutputTree->Branch("SCEtaWidth",&fEleSCEtaWidth,"SCEtaWidth/F"); 
  FakeEleSourceOutputTree->Branch("SCPhiWidth",&fEleSCPhiWidth,"SCPhiWidth/F"); 
  FakeEleSourceOutputTree->Branch("GsfTrackChi2OverNdof",&fEleGsfTrackChi2OverNdof,"GsfTrackChi2OverNdof/F"); 
  FakeEleSourceOutputTree->Branch("R9",&fEleR9,"R9/F"); 

  FakeEleSourceOutputTree->Branch("SeedEMaxOverE",&fEleSeedEMaxOverE,"SeedEMaxOverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedETopOverE",&fEleSeedETopOverE,"SeedETopOverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedEBottomOverE",&fEleSeedEBottomOverE,"SeedEBottomOverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedELeftOverE",&fEleSeedELeftOverE,"SeedELeftOverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedERightOverE",&fEleSeedERightOverE,"SeedERightOverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedE2ndOverE",&fEleSeedE2ndOverE,"SeedE2ndOverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedE2x5RightOverE",&fEleSeedE2x5RightOverE,"SeedE2x5RightOverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedE2x5LeftOverE",&fEleSeedE2x5LeftOverE,"SeedE2x5LeftOverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedE2x5TopOverE",&fEleSeedE2x5TopOverE,"SeedE2x5TopOverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedE2x5BottomOverE",&fEleSeedE2x5BottomOverE,"SeedE2x5BottomOverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedE2x5MaxOverE",&fEleSeedE2x5MaxOverE,"SeedE2x5MaxOverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedE1x3OverE",&fEleSeedE1x3OverE,"SeedE1x3OverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedE3x1OverE",&fEleSeedE3x1OverE,"SeedE3x1OverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedE1x5OverE",&fEleSeedE1x5OverE,"SeedE1x5OverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedE2x2OverE",&fEleSeedE2x2OverE,"SeedE2x2OverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedE3x2OverE",&fEleSeedE3x2OverE,"SeedE3x2OverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedE3x3OverE",&fEleSeedE3x3OverE,"SeedE3x3OverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedE4x4OverE",&fEleSeedE4x4OverE,"SeedE4x4OverE/F"); 
  FakeEleSourceOutputTree->Branch("SeedE5x5OverE",&fEleSeedE5x5OverE,"SeedE5x5OverE/F"); 

  //Isolation Variables
  FakeEleSourceOutputTree->Branch("StandardLikelihood",&fEleStandardLikelihood,"StandardLikelihood/F"); 
  FakeEleSourceOutputTree->Branch("PFMVA",&fElePFMVA,"PFMVA/F"); 
  FakeEleSourceOutputTree->Branch("ChargedIso03",&fEleChargedIso03,"ChargedIso03/F"); 
  FakeEleSourceOutputTree->Branch("NeutralHadronIso03",&fEleNeutralHadronIso03,"NeutralHadronIso03/F"); 
  FakeEleSourceOutputTree->Branch("GammaIso03",&fEleGammaIso03,"GammaIso03/F"); 
  FakeEleSourceOutputTree->Branch("ChargedIso04",&fEleChargedIso04,"ChargedIso04/F"); 
  FakeEleSourceOutputTree->Branch("NeutralHadronIso04",&fEleNeutralHadronIso04,"NeutralHadronIso04/F"); 
  FakeEleSourceOutputTree->Branch("GammaIso04",&fEleGammaIso04,"GammaIso04/F"); 
  FakeEleSourceOutputTree->Branch("ChargedIso04FromOtherVertices",&fEleChargedIso04FromOtherVertices,"ChargedIso04FromOtherVertices/F"); 
  FakeEleSourceOutputTree->Branch("NeutralHadronIso04_10Threshold",&fEleNeutralHadronIso04_10Threshold,"NeutralHadronIso04_10Threshold/F"); 
  FakeEleSourceOutputTree->Branch("GammaIso04_10Threshold",&fEleGammaIso04_10Threshold,"GammaIso04_10Threshold/F"); 
  FakeEleSourceOutputTree->Branch("TrkIso03",&fEleTrkIso03,"TrkIso03/F"); 
  FakeEleSourceOutputTree->Branch("EMIso03",&fEleEMIso03,"EMIso03/F"); 
  FakeEleSourceOutputTree->Branch("HadIso03",&fEleHadIso03,"HadIso03/F"); 
  FakeEleSourceOutputTree->Branch("TrkIso04",&fEleTrkIso04,"TrkIso04/F"); 
  FakeEleSourceOutputTree->Branch("EMIso04",&fEleEMIso04,"EMIso04/F"); 
  FakeEleSourceOutputTree->Branch("HadIso04",&fEleHadIso04,"HadIso04/F"); 
  FakeEleSourceOutputTree->Branch("Rho",&fRho,"Rho/F"); 
  FakeEleSourceOutputTree->Branch("NVertices",&fNVertices,"NVertices/F"); 


  //*****************************************************************************************
  //FakeEleTargetTree
  //*****************************************************************************************
  TFile *FakeEleSourceFile = new TFile("ElectronSelectionTraining.Fakes.root", "READ");
  TTree *FakeEleSourceTree = (TTree*)FakeEleSourceFile->Get("Electrons");
  FakeEleSourceTree->SetBranchAddress( "weight", &fWeight);
  FakeEleSourceTree->SetBranchAddress( "run", &fRunNumber);
  FakeEleSourceTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  FakeEleSourceTree->SetBranchAddress( "event", &fEventNumber);
  FakeEleSourceTree->SetBranchAddress( "pt", &fElePt); 
  FakeEleSourceTree->SetBranchAddress( "eta", &fEleEta); 
  FakeEleSourceTree->SetBranchAddress( "phi", &fElePhi); 
  FakeEleSourceTree->SetBranchAddress( "scet", &fEleSCEt); 
  FakeEleSourceTree->SetBranchAddress( "sceta", &fEleSCEta); 
  FakeEleSourceTree->SetBranchAddress( "scphi", &fEleSCPhi); 
  FakeEleSourceTree->SetBranchAddress( "pfiso", &fElePFIso); 
  FakeEleSourceTree->SetBranchAddress( "SigmaIEtaIEta", &fEleSigmaIEtaIEta); 
  FakeEleSourceTree->SetBranchAddress( "DEtaIn", &fEleDEtaIn); 
  FakeEleSourceTree->SetBranchAddress( "DPhiIn", &fEleDPhiIn); 
  FakeEleSourceTree->SetBranchAddress( "HoverE", &fEleHoverE); 
  FakeEleSourceTree->SetBranchAddress( "D0", &fEleD0); 
  FakeEleSourceTree->SetBranchAddress( "DZ", &fEleDZ); 
  FakeEleSourceTree->SetBranchAddress( "FBrem", &fEleFBrem); 
  FakeEleSourceTree->SetBranchAddress( "EOverP", &fEleEOverP); 
  FakeEleSourceTree->SetBranchAddress( "ESeedClusterOverPout", &fEleESeedClusterOverPout); 
  FakeEleSourceTree->SetBranchAddress( "SigmaIPhiIPhi", &fEleSigmaIPhiIPhi); 
  FakeEleSourceTree->SetBranchAddress( "NBrem", &fEleNBrem); 
  FakeEleSourceTree->SetBranchAddress( "OneOverEMinusOneOverP", &fEleOneOverEMinusOneOverP); 
  FakeEleSourceTree->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
  FakeEleSourceTree->SetBranchAddress( "IP3d", &fEleIP3d); 
  FakeEleSourceTree->SetBranchAddress( "IP3dSig", &fEleIP3dSig); 
  FakeEleSourceTree->SetBranchAddress( "HcalDepth1OverEcal", &fEleHcalDepth1OverEcal); 
  FakeEleSourceTree->SetBranchAddress( "HcalDepth2OverEcal", &fEleHcalDepth2OverEcal); 
  FakeEleSourceTree->SetBranchAddress( "dEtaCalo", &fEledEtaCalo); 
  FakeEleSourceTree->SetBranchAddress( "dPhiCalo", &fEledPhiCalo); 
  FakeEleSourceTree->SetBranchAddress( "PreShowerOverRaw", &fElePreShowerOverRaw); 
  FakeEleSourceTree->SetBranchAddress( "CovIEtaIPhi", &fEleCovIEtaIPhi); 
  FakeEleSourceTree->SetBranchAddress( "SCEtaWidth", &fEleSCEtaWidth); 
  FakeEleSourceTree->SetBranchAddress( "SCPhiWidth", &fEleSCPhiWidth); 
  FakeEleSourceTree->SetBranchAddress( "GsfTrackChi2OverNdof", &fEleGsfTrackChi2OverNdof); 
  FakeEleSourceTree->SetBranchAddress( "R9", &fEleR9); 
  FakeEleSourceTree->SetBranchAddress( "SeedEMaxOverE", &fEleSeedEMaxOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedETopOverE", &fEleSeedETopOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedEBottomOverE", &fEleSeedEBottomOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedELeftOverE", &fEleSeedELeftOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedERightOverE", &fEleSeedERightOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2ndOverE", &fEleSeedE2ndOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x5RightOverE", &fEleSeedE2x5RightOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x5LeftOverE", &fEleSeedE2x5LeftOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x5TopOverE", &fEleSeedE2x5TopOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x5BottomOverE", &fEleSeedE2x5BottomOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x5MaxOverE", &fEleSeedE2x5MaxOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE1x3OverE", &fEleSeedE1x3OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE3x1OverE", &fEleSeedE3x1OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE1x5OverE", &fEleSeedE1x5OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x2OverE", &fEleSeedE2x2OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE3x2OverE", &fEleSeedE3x2OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE3x3OverE", &fEleSeedE3x3OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE4x4OverE", &fEleSeedE4x4OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE5x5OverE", &fEleSeedE5x5OverE); 
  FakeEleSourceTree->SetBranchAddress( "StandardLikelihood", &fEleStandardLikelihood); 
  FakeEleSourceTree->SetBranchAddress( "PFMVA", &fElePFMVA); 
  FakeEleSourceTree->SetBranchAddress( "ChargedIso03", &fEleChargedIso03); 
  FakeEleSourceTree->SetBranchAddress( "NeutralHadronIso03", &fEleNeutralHadronIso03); 
  FakeEleSourceTree->SetBranchAddress( "GammaIso03", &fEleGammaIso03); 
  FakeEleSourceTree->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
  FakeEleSourceTree->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
  FakeEleSourceTree->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
  FakeEleSourceTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  FakeEleSourceTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  FakeEleSourceTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  FakeEleSourceTree->SetBranchAddress( "TrkIso03", &fEleTrkIso03); 
  FakeEleSourceTree->SetBranchAddress( "EMIso03", &fEleEMIso03); 
  FakeEleSourceTree->SetBranchAddress( "HadIso03", &fEleHadIso03); 
  FakeEleSourceTree->SetBranchAddress( "TrkIso04", &fEleTrkIso04); 
  FakeEleSourceTree->SetBranchAddress( "EMIso04", &fEleEMIso04); 
  FakeEleSourceTree->SetBranchAddress( "HadIso04", &fEleHadIso04); 
  FakeEleSourceTree->SetBranchAddress( "Rho", &fRho); 
  FakeEleSourceTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < FakeEleSourceTree->GetEntries(); ientry++) {       	
    FakeEleSourceTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Eelectron " << ientry << endl;
    
    fWeight = FakeElectronPtReweightFactor->GetBinContent(FakeElectronPtReweightFactor->GetXaxis()->FindFixBin(fElePt));
    FakeEleSourceOutputTree->Fill();
  } //loop over electrons
  FakeEleSourceOutputFile->Write();
  FakeEleSourceOutputFile->Close();
  delete FakeEleSourceFile;




  gBenchmark->Show("WWTemplate");       
} 



//*************************************************************************************************
//Main Function
//*************************************************************************************************
void MakeElectronPtSpectrum() {
  FillLeptonPtSpectrum();
  MakeReweightFactors();
  UpdateNtupleWeights();

}
