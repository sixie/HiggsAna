//================================================================================================
//
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <iostream>                 // standard I/O

#include "HiggsAna/CommonData/interface/ElectronTree.h"

#endif

//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void DoSplit(string inputFileName, string outputFileName, Int_t Option)
{  


  //*****************************************************************************************
  // Load Input
  //*****************************************************************************************
  ElectronTree eleTree;
  eleTree.LoadTree(inputFileName.c_str());
  eleTree.InitTree();
  cout << "Events in the ntuple: " << eleTree.tree_->GetEntries() << endl;


  //*************************************************************************************************
  //Output tree
  //*************************************************************************************************
  cout << "Output File : " << outputFileName << endl;
  TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
  outputFile->cd();
  TTree *outputTree = eleTree.tree_->CloneTree(0);

  for(UInt_t ientry=0; ientry < eleTree.tree_->GetEntries(); ientry++) {       	
    eleTree.tree_->GetEntry(ientry);
    if (ientry % 100000 == 0) cout << "Electron " << ientry << endl;

    if (Option == 0) {
      if (eleTree.fEventNumber % 2 == 0) {
        outputTree->Fill();
      }
    } 
    if (Option == 1) {
      if (eleTree.fEventNumber % 2 != 0) {
        outputTree->Fill();
      }
    }
  } 
  cout << "Events in output ntuple: " << outputTree->GetEntries() << endl;
  outputFile->Write();
  outputFile->Close();
  delete outputFile;

} 



//*************************************************************************************************
//Main Function
//*************************************************************************************************
void SplitElectronNtuples() {

  DoSplit("ElectronSelectionTraining.Fake_ZPlusJet_2012.root","ElectronSelectionTraining.Fake_ZPlusJet_2012.Training.root",0);
  DoSplit("ElectronSelectionTraining.Fake_ZPlusJet_2012.root","ElectronSelectionTraining.Fake_ZPlusJet_2012.Testing.root",1);
  
}

void SplitElectronNtuples(string inputFileName, string outputFileName, Int_t Option) {

  DoSplit(inputFileName,outputFileName,Option);
  
}


