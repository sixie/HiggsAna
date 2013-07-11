//root -l HiggsAna/HZZ4l/HLL/CreateEfficiencyMap.C+'("HZZEfficiencyMap_HZZ125.root","HZZ125"0)'


//================================================================================================
//
// Create Efficiency Map
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
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

// data structs
#include "HiggsAna/HZZ4l/HLL/LeptonResolutionData/ElectronEfficiencyMap.h"
#include "HiggsAna/HZZ4l/HLL/LeptonResolutionData/MuonEfficiencyMap.h"
#include "LeptonResponseMap.hh"

#endif



//=== MAIN MACRO =================================================================================================

void LoadLeptonPtModel(Int_t Option = 0) 
{  

//   TFile *LeptonResponseFile = new TFile("HiggsAna/HZZ4l/HLL/LeptonResolutionData/PtResolutionModel_ZZ.root","READ");
  TFile *LeptonResponseFile = 0;

  if (Option == 0) {
    LeptonResponseFile = new TFile("PtResolutionModel_ZZ.root","READ");
  } else if (Option == 1) {
    LeptonResponseFile = new TFile("PtResolutionModel_ZZ.root","UPDATE");
  }
  assert(LeptonResponseFile);

  const UInt_t NPtBins = 14;
  const UInt_t NEtaBins = 16;
  Double_t ptBins[15] = {5,7,8,9,10,12,14,16,18,20,25,30,35,40,50};
  Double_t etaBins[17] = {0,0.2,0.4,0.6,0.8,1,1.2,1.4442,1.566,1.8,2,2.1,2.2,2.3,2.4,2.5,2.6};

  TH2F* DoubleSidedCBShapeParamArray_Electrons_mean = 0;
  TH2F* DoubleSidedCBShapeParamArray_Electrons_sigma = 0;
  TH2F* DoubleSidedCBShapeParamArray_Electrons_alphaL = 0;
  TH2F* DoubleSidedCBShapeParamArray_Electrons_nL = 0;
  TH2F* DoubleSidedCBShapeParamArray_Electrons_alphaR = 0;
  TH2F* DoubleSidedCBShapeParamArray_Electrons_nR = 0;

  if (Option == 0) {
    DoubleSidedCBShapeParamArray_Electrons_mean = (TH2F*)LeptonResponseFile->Get("DoubleSidedCBShapeParamArray_Electrons_mean");
    DoubleSidedCBShapeParamArray_Electrons_sigma = (TH2F*)LeptonResponseFile->Get("DoubleSidedCBShapeParamArray_Electrons_sigma");
    DoubleSidedCBShapeParamArray_Electrons_alphaL = (TH2F*)LeptonResponseFile->Get("DoubleSidedCBShapeParamArray_Electrons_alphaL");
    DoubleSidedCBShapeParamArray_Electrons_nL = (TH2F*)LeptonResponseFile->Get("DoubleSidedCBShapeParamArray_Electrons_nL");
    DoubleSidedCBShapeParamArray_Electrons_alphaR = (TH2F*)LeptonResponseFile->Get("DoubleSidedCBShapeParamArray_Electrons_alphaR");
    DoubleSidedCBShapeParamArray_Electrons_nR = (TH2F*)LeptonResponseFile->Get("DoubleSidedCBShapeParamArray_Electrons_nR");
  }

  if (Option == 1) {
    DoubleSidedCBShapeParamArray_Electrons_mean = new TH2F( "DoubleSidedCBShapeParamArray_Electrons_mean", "", NPtBins, 0, NPtBins, NEtaBins, 0, NEtaBins);
    DoubleSidedCBShapeParamArray_Electrons_sigma = new TH2F( "DoubleSidedCBShapeParamArray_Electrons_sigma", "", NPtBins, 0, NPtBins, NEtaBins, 0, NEtaBins);
    DoubleSidedCBShapeParamArray_Electrons_alphaL = new TH2F( "DoubleSidedCBShapeParamArray_Electrons_alphaL", "", NPtBins, 0, NPtBins, NEtaBins, 0, NEtaBins);
    DoubleSidedCBShapeParamArray_Electrons_nL = new TH2F( "DoubleSidedCBShapeParamArray_Electrons_nL", "", NPtBins, 0, NPtBins, NEtaBins, 0, NEtaBins);
    DoubleSidedCBShapeParamArray_Electrons_alphaR = new TH2F( "DoubleSidedCBShapeParamArray_Electrons_alphaR", "", NPtBins, 0, NPtBins, NEtaBins, 0, NEtaBins);
    DoubleSidedCBShapeParamArray_Electrons_nR = new TH2F( "DoubleSidedCBShapeParamArray_Electrons_nR", "", NPtBins, 0, NPtBins, NEtaBins, 0, NEtaBins);
  }

  for (uint i=0; i < NPtBins+2; ++i) {
    for (uint j=0; j < NEtaBins+2; ++j) {
      RooWorkspace *w = (RooWorkspace*)LeptonResponseFile->Get(Form("LeptonPtResolutionModel_%s_PtBin%d_EtaBin%d","Electrons",i, j));
      if (!w) { 
        cout << "cannot load workspace: " << Form("LeptonPtResolutionModel_%s_PtBin%d_EtaBin%d","Electrons",i, j) << endl;
        assert(w);
      }
      RooRealVar     *mean  = w->var("mean");
      RooRealVar     *sigma = w->var("sigma");
      RooRealVar     *alphaL= w->var("alphaL");
      RooRealVar     *nL    = w->var("nL");
      RooRealVar     *alphaR= w->var("alphaR");
      RooRealVar     *nR    = w->var("nR");
      
      DoubleSidedCBShapeParamArray_Electrons_mean->SetBinContent(i,j,mean->getVal());
      DoubleSidedCBShapeParamArray_Electrons_sigma->SetBinContent(i,j,sigma->getVal());
      DoubleSidedCBShapeParamArray_Electrons_alphaL->SetBinContent(i,j,alphaL->getVal());
      DoubleSidedCBShapeParamArray_Electrons_nL->SetBinContent(i,j,nL->getVal());
      DoubleSidedCBShapeParamArray_Electrons_alphaR->SetBinContent(i,j,alphaR->getVal());
      DoubleSidedCBShapeParamArray_Electrons_nR->SetBinContent(i,j,nR->getVal());
      //cout << "Model Bin " << i << " " << j << " : " << mean->getVal() << " : " << DoubleSidedCBShapeParamArray_Electrons_mean->GetBinContent(i,j) << endl;
    }
  }


  for (uint i=0; i < NPtBins+2; ++i) {
    for (uint j=0; j < NEtaBins+2; ++j) {
      cout << "Bin " << i << " " << j << " : " << DoubleSidedCBShapeParamArray_Electrons_mean->GetBinContent(i,j) << endl;
    }
  }

  if (Option == 1) {
    LeptonResponseFile->WriteTObject(DoubleSidedCBShapeParamArray_Electrons_mean, "DoubleSidedCBShapeParamArray_Electrons_mean", "WriteDelete");
    LeptonResponseFile->WriteTObject(DoubleSidedCBShapeParamArray_Electrons_sigma, "DoubleSidedCBShapeParamArray_Electrons_sigma", "WriteDelete");
    LeptonResponseFile->WriteTObject(DoubleSidedCBShapeParamArray_Electrons_alphaL, "DoubleSidedCBShapeParamArray_Electrons_alphaL", "WriteDelete");
    LeptonResponseFile->WriteTObject(DoubleSidedCBShapeParamArray_Electrons_nL, "DoubleSidedCBShapeParamArray_Electrons_nL", "WriteDelete");
    LeptonResponseFile->WriteTObject(DoubleSidedCBShapeParamArray_Electrons_alphaR, "DoubleSidedCBShapeParamArray_Electrons_alphaR", "WriteDelete");
    LeptonResponseFile->WriteTObject(DoubleSidedCBShapeParamArray_Electrons_nR, "DoubleSidedCBShapeParamArray_Electrons_nR", "WriteDelete");
  } 



} 


