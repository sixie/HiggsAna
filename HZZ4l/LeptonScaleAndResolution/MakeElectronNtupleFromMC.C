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

// define structures to read in ntuple
#include "HiggsAna/DataTree/interface/HiggsAnaDefs.hh"
#include "HiggsAna/DataTree/interface/TEventInfo.hh"
#include "HiggsAna/DataTree/interface/TElectron.hh"
#include "HiggsAna/DataTree/interface/TPhoton.hh"
#include "HiggsAna/DataTree/interface/TMuon.hh"
#include "HiggsAna/DataTree/interface/TJet.hh"
#include "HiggsAna/DataTree/interface/TPFCandidate.hh"
#include "HiggsAna/DataTree/interface/TGenParticle.hh"
#include "HiggsAna/Utils/LeptonIDCuts.hh"
#include "CITCommon/CommonData/interface/ElectronTree.h"
#include "HiggsAna/Utils/LeptonTools.hh"
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "HiggsAna/HZZ4l/Utils/LeptonSelection.hh"


#endif

//*************************************************************************************************
//Convert int to string
//*************************************************************************************************
string IntToString(int i) {
  char temp[100];
  sprintf(temp, "%d", i);
  string str = temp;
  return str;
}


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
void MakeNtuple(const string inputFilename,  const string outputFilename, Int_t Option = 0);

//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = TFile::Open(infname,"read");
  assert(inf);

  TTree* t = (TTree*)inf->Get(tname);
  
  if (!t) {
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("HwwNtuplerMod");
    if (!dir) {
      cout << "Cannot get Directory HwwNtuplerMod from file " << infname << endl;
      assert(dir);
    }
    t = (TTree*)dir->Get(tname);
  }

  if (!t) {
    cout << "Cannot get Tree with name " << tname << " from file " << infname << endl;
  }
  assert(t);


  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}


//=== MAIN MACRO =================================================================================================
void MakeElectronNtupleFromMC(Int_t Sample = 0) {
 
  //HZZ 
  if (Sample == 1) {
    MakeNtuple("/data/smurf/sixie/hist/HZZ4lNtuples_old/mc/AllNtuple_HZZ4lNtuple_s12-h125zz4l-gf-v9_noskim_0004.root","ElectronSelectionTraining.s11HZZ125.Training.root", 0);
  }

  //HWW
  if (Sample == 2) {
    
  }

} 


void MakeElectronNtupleFromMC(const string inputFilename, const string outputFilename, Int_t Option = 0) {
  MakeNtuple(inputFilename,outputFilename, Option);
}



void MakeNtuple(const string inputFilename, const string outputFilename, Int_t Option)
{  
  gBenchmark->Start("WWTemplate");

  //*****************************************************************************************
  //Define Data Era
  //*****************************************************************************************
  UInt_t DataEra = kDataEra_NONE;
  DataEra = kDataEra_2012_MC;

  if (Option == 1) {
    DataEra = kDataEra_NONE;    
  }

  cout << "using DataEra = " << DataEra << endl;

  //*****************************************************************************************
  //Setup MVA
  //*****************************************************************************************
  EGammaMvaEleEstimator *eleIDMVA = new EGammaMvaEleEstimator();
  vector<string> weightFiles;

  weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat1.weights.xml")).c_str());
  weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat2.weights.xml")).c_str());
  weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat3.weights.xml")).c_str());
  weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat4.weights.xml")).c_str());
  weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat5.weights.xml")).c_str());
  weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat6.weights.xml")).c_str());
  eleIDMVA->initialize( "BDT", EGammaMvaEleEstimator::kNonTrig,  kTRUE, weightFiles);

  //*****************************************************************************************
  //Setup Output Tree
  //*****************************************************************************************
  TFile *outputFile = new TFile(outputFilename.c_str(), "RECREATE");
  citana::ElectronTree *eleTree = new citana::ElectronTree;
  eleTree->CreateTree();
  eleTree->tree_->SetAutoFlush(0);

  UInt_t NElectronsFilled = 0;
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  higgsana::TEventInfo *info    = new higgsana::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("higgsana::TElectron");
  TClonesArray *pfcandidateArr = new TClonesArray("higgsana::TPFCandidate");
  TClonesArray *genparticleArr = new TClonesArray("higgsana::TGenParticle");

  Int_t NEvents = 0;



    //********************************************************
    // Get Tree
    //********************************************************
    eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
    TBranch *infoBr;
    TBranch *electronBr;
    TBranch *muonBr;
    TBranch *pfcandidateBr;
    TBranch *genparticleBr;


    //*****************************************************************************************
    //Loop over Data Tree
    //*****************************************************************************************
    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr); pfcandidateBr = eventTree->GetBranch("PFCandidate");
    eventTree->SetBranchAddress("GenParticle", &genparticleArr); genparticleBr = eventTree->GetBranch("GenParticle");

    cout << "InputFile " << inputFilename << " --- Total Events : " << eventTree->GetEntries() << endl;
//     for(UInt_t ientry=0; ientry < 10000; ientry++) {       	
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 10000 == 0) cout << "Event " << ientry << endl;

      Double_t eventweight = info->eventweight;

      //********************************************************
      // Load the branches
      //********************************************************
      electronArr->Clear(); 
      pfcandidateArr->Clear(); 
      genparticleArr->Clear(); 
      electronBr->GetEntry(ientry);
      pfcandidateBr->GetEntry(ientry);
      genparticleBr->GetEntry(ientry);


      //********************************************************
      // Pileup Energy Density
      //********************************************************
      Double_t rhoEleIso = 0;
      UInt_t EleEAEra = 0;

      if (DataEra == kDataEra_2011_MC) {     
        if (!(isnan(info->RhoKt6PFJetsForIso25) || 
              isinf(info->RhoKt6PFJetsForIso25))) {
          rhoEleIso = info->RhoKt6PFJetsForIso25;
        }
        EleEAEra = kDataEra_2011_Data;
      } else if (DataEra == kDataEra_2012_MC) {
        if (!(isnan(info->RhoKt6PFJets) || 
              isinf(info->RhoKt6PFJets))) {
          rhoEleIso = info->RhoKt6PFJets;
        }
        EleEAEra = kDataEra_2012_Data;
      }


      //********************************************************
      // Loop Over Electrons
      //********************************************************
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const higgsana::TElectron *ele = (higgsana::TElectron*)((*electronArr)[i]);

        //********************************************************
        // Select MC Truth Electrons
        //********************************************************              
        if (!((UInt_t(abs(max(0,ele->isMCReal))) & 1) == 1)) continue;
        
        //********************************************************
        // Do HZZ Selection
        //********************************************************              
        Bool_t passID = (PassEleHZZ4lPreselection(ele) && PassEleHZZ4lICHEP2012ID(ele, eleIDMVA));
        vector<const higgsana::TPFCandidate*> photonsToVeto;
        Bool_t passIsolation = PassEleHZZ4lICHEP2012Iso(ele,i,pfcandidateArr,rhoEleIso,EleEAEra,photonsToVeto);
        if (Option == 1) {
          if (!(passID && passIsolation)) continue;
        }

        //********************************************************************************************
        //do not fill electrons that decayed from a tau
        //********************************************************************************************
        double minDR = 9999;
        Bool_t IsPromptGenElectron = kFALSE;
        const higgsana::TGenParticle *matchedGenElectron =0;
        for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
          const higgsana::TGenParticle *gen = (higgsana::TGenParticle*)((*genparticleArr)[k]);

          //status 1 match
          if (abs(gen->pdgid) == 11 && (gen->status == 1)) {
            double tmpDR = higgsana::deltaR( ele->eta, ele->phi, gen->eta , gen->phi);
            if (tmpDR < minDR && tmpDR < 0.3) {
              minDR = tmpDR;
              matchedGenElectron = gen;
            }
          }
        }
        if (!matchedGenElectron || abs(matchedGenElectron->motherPdgID) == 15) continue;

        //Fill These Electrons
        NElectronsFilled++;
        FillElectronTree( eleTree, ele, i, pfcandidateArr, eleIDMVA, genparticleArr, rhoEleIso, EleEAEra, 
                          info->nPV0, info->runNum, info->lumiSec, info->evtNum, 1.0);
        
      } //loop over electrons

    }

    cout << "Total Electrons: " << NElectronsFilled << endl;


    delete info;
    delete electronArr;
    delete pfcandidateArr;
    delete genparticleArr;

    cout << "Total Electrons: " << NElectronsFilled << endl;

    outputFile->Write();
    outputFile->Close();
  
    delete outputFile;
    if (eleTree) delete eleTree;

    gBenchmark->Show("WWTemplate");       
} 


