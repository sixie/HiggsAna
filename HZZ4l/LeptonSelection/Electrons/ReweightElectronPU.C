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
void DoFakeElectronPUReweighting(string TargetFileName, string SourceFilename, string OutputFilename)
{  


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *histRhoRealElectronTarget = new TH1F("RhoRealElectronTarget", "; #rho [GeV] ; Number of Events ",  50, 0 , 50);
  TH1F *histRhoFakeElectronSource = new TH1F("RhoFakeElectronSource", "; #rho [GeV] ; Number of Events ",  50, 0 , 50);


  //*****************************************************************************************
  //Generate Target PU distribution from Real Electron Sample
  //*****************************************************************************************
  ElectronTree RealEleTargetTree;
  RealEleTargetTree.LoadTree(TargetFileName.c_str());
  RealEleTargetTree.InitTree();

  for(UInt_t ientry=0; ientry < RealEleTargetTree.tree_->GetEntries(); ientry++) {       	
    RealEleTargetTree.tree_->GetEntry(ientry);    
    if (ientry % 100000 == 0) cout << "Real Electron Sample " << ientry << endl;
    histRhoRealElectronTarget->Fill(RealEleTargetTree.fRho, RealEleTargetTree.fWeight);
  }


  //*****************************************************************************************
  //Generate Source PU distribution from Fake Electron Sample
  //*****************************************************************************************
  ElectronTree FakeEleSourceTree;
  FakeEleSourceTree.LoadTree(SourceFilename.c_str());
  FakeEleSourceTree.InitTree();

  for(UInt_t ientry=0; ientry < FakeEleSourceTree.tree_->GetEntries(); ientry++) {       	
    FakeEleSourceTree.tree_->GetEntry(ientry);    
    if (ientry % 100000 == 0) cout << "Fake Electron Sample " << ientry << endl;
    histRhoFakeElectronSource->Fill(FakeEleSourceTree.fRho, FakeEleSourceTree.fWeight);
  }


  //*****************************************************************************************
  //Create PU reweight histogram
  //*****************************************************************************************
  NormalizeHist(histRhoRealElectronTarget);
  NormalizeHist(histRhoFakeElectronSource);

  TH1F *FakeElectronPUReweightFactor = (TH1F*)histRhoFakeElectronSource->Clone("FakeElectronPUReweightFactor");
  FakeElectronPUReweightFactor->SetBinContent(0,1.0);
  for(UInt_t a=1; a < FakeElectronPUReweightFactor->GetXaxis()->GetNbins()+2; ++a) {
    if (histRhoFakeElectronSource->GetBinContent(a)>0) {
      FakeElectronPUReweightFactor->SetBinContent(a,histRhoRealElectronTarget->GetBinContent(a) / histRhoFakeElectronSource->GetBinContent(a));
    } else {
      FakeElectronPUReweightFactor->SetBinContent(a,1.0);
    }
  }



  //*************************************************************************************************
  //Output tree
  //*************************************************************************************************
  TFile *FakeEleOutputFile = new TFile(OutputFilename.c_str(), "RECREATE");
  TTree *FakeEleOutputTree = FakeEleSourceTree.tree_->CloneTree(0);
  FakeEleOutputTree->SetAutoFlush(0);

  //*****************************************************************************************
  //FakeEleSourceTree
  //*****************************************************************************************
  for(UInt_t ientry=0; ientry < FakeEleSourceTree.tree_->GetEntries(); ientry++) {       	
    FakeEleSourceTree.tree_->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Electron " << ientry << endl;
    FakeEleSourceTree.fWeight = FakeEleSourceTree.fWeight * FakeElectronPUReweightFactor->GetBinContent(FakeElectronPUReweightFactor->GetXaxis()->FindFixBin(FakeEleSourceTree.fRho));
    FakeEleOutputTree->Fill();
  }
  FakeEleOutputFile->Write();
  FakeEleOutputFile->Close();

  gBenchmark->Show("WWTemplate");       
}



//*************************************************************************************************
//Main Function
//*************************************************************************************************
void ReweightElectronPU() {

  DoFakeElectronPUReweighting("ElectronSelectionTraining.Real.weighted.root", "ElectronSelectionTraining.Fake.weighted.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.root" );


}

void ReweightElectronPU(string TargetFileName, string SourceFilename, string OutputFilename) {

  DoFakeElectronPUReweighting(TargetFileName, SourceFilename, OutputFilename);

}
