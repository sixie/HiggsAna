//root -l -b -q EWKAna/Hww/FakeRate/ComputeElectronFakeRate_Data.C+\(\)
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
#include "HiggsAna/Utils/LeptonIDCuts.hh"
#include "HiggsAna/DataTree/interface/Types.h"
#include "HiggsAna/Utils/LeptonTools.hh"
#include "HiggsAna/CommonData/interface/ElectronTree.h"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

#endif


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
void MakeFakeElectronTrainingNtuple(Int_t Option = 0) {

  if (Option == 0) {
    MakeNtuple("2012Data","/data/blue/sixie/LeptonSelection/Electrons/ElectronSelectionTraining.QCDFakes_2012.root", Option);
  }
  if (Option == 1) {
    MakeNtuple("2012Data","/data/blue/sixie/LeptonSelection/ElectronsNoRhoCorr/ElectronSelectionTraining.QCDFakes_2012.root", Option);
  }


}


void MakeNtuple(const string inputFilename, const string outputFilename, Int_t Option)
{  
  gBenchmark->Start("WWTemplate");

  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  TFile *outputFile = new TFile(outputFilename.c_str(), "RECREATE");
  ElectronTree *eleTree = new ElectronTree;
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
  TClonesArray *muonArr = new TClonesArray("higgsana::TMuon");
  TClonesArray *jetArr = new TClonesArray("higgsana::TJet");
  TClonesArray *photonArr = new TClonesArray("higgsana::TPhoton");
  TClonesArray *pfcandidateArr = new TClonesArray("higgsana::TPFCandidate");
  
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("/data/smurf/data/Winter11_4700ipb/auxiliar/hww.Full2011.json"); 
  rlrm.AddJSONFile("/data/blue/sixie/HZZ4l/auxiliar/2012/Cert_190456-196531_8TeV_PromptReco_Collisions12_JSON.txt");

  Int_t NEvents = 0;

  UInt_t DataEra = kDataEra_NONE;

  vector<string> inputfiles;
  if (inputFilename == "LIST") {
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-del-m10-v1.FakeTriggerSkim.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-del-pr-v4.FakeTriggerSkim.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-del-a05-v1.FakeTriggerSkim.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-del-o03-v1.FakeTriggerSkim.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11b-del-pr-v1.FakeTriggerSkim.root");    
  } 
  else if (inputFilename == "2012Data") {
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r12a-del-pr-v1_FakeRateTriggerSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r12b-del-pr-v1_FakeRateTriggerSkimmed.root");
    DataEra = kDataEra_2012_MC;
  }
  else {
    inputfiles.push_back(inputFilename);
  }

  for (UInt_t f = 0; f < inputfiles.size(); ++f) {

    //********************************************************
    // Get Tree
    //********************************************************
    eventTree = getTreeFromFile(inputfiles[f].c_str(),"Events"); 
    TBranch *infoBr;
    TBranch *electronBr;
    TBranch *muonBr;
    TBranch *jetBr;
    TBranch *photonBr;
    TBranch *pfcandidateBr;


    //*****************************************************************************************
    //Loop over Data Tree
    //*****************************************************************************************
    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("Photon", &photonArr);     photonBr = eventTree->GetBranch("Photon");
    eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");
    eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr); pfcandidateBr = eventTree->GetBranch("PFCandidate");

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      Double_t eventweight = info->eventweight;

      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...

      NEvents++;

      //********************************************************
      // Load the branches
      //********************************************************
      electronArr->Clear(); 
      muonArr->Clear(); 
      photonArr->Clear(); 
      jetArr->Clear(); 
      pfcandidateArr->Clear(); 
      electronBr->GetEntry(ientry);
      muonBr->GetEntry(ientry);
      photonBr->GetEntry(ientry);
      jetBr->GetEntry(ientry);
      pfcandidateBr->GetEntry(ientry);


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
      // TcMet
      //********************************************************
      TVector3 pfMet;        
      if(info->pfMEx!=0 || info->pfMEy!=0) {       
        pfMet.SetXYZ(info->pfMEx, info->pfMEy, 0);
      }
      Double_t met = pfMet.Pt();

      Int_t NElectrons = electronArr->GetEntries();
 

      //********************************************************
      // Event Selection Cuts
      //********************************************************
      //veto events with more than 1 reco electron
      if (NElectrons > 1) continue;
      //met cut removed W events
      if (met > 20) continue;


      //******************************************************************************
      //loop over electrons 
      //******************************************************************************
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const higgsana::TElectron *ele = (higgsana::TElectron*)((*electronArr)[i]);
 

        //make cut on dz
        if (fabs(ele->dz) > 0.1) continue;

        //protect against pathologies
        if (TMath::IsNaN(ele->sigiPhiiPhi)) {
          cout << "Pathological SigmaIPhiIPhi : " 
               << info->runNum << " " << info->lumiSec << " " << info->evtNum << endl;
          continue;
        }
        
        //********************************************************
        //find leading jet in the event
        //********************************************************
        Double_t leadingJetPt = -1;
        //pass event selection     
        for(Int_t j=0; j<jetArr->GetEntries(); j++) {
          const higgsana::TJet *jet = (higgsana::TJet*)((*jetArr)[j]);        
          if (jet->pt > leadingJetPt &&
              higgsana::deltaR(jet->eta, jet->phi, ele->eta, ele->phi) > 1.0) {
            leadingJetPt = jet->pt;          
          }
        }
      
        //Fill These Electrons
        NElectronsFilled++;
        
        if (Option == 0) {
          FillElectronTree( eleTree, ele, pfcandidateArr, rhoEleIso, EleEAEra, 
                            info->nPV0, info->runNum, info->lumiSec, info->evtNum);
        } else if (Option == 1) {
          FillElectronTree( eleTree, ele, pfcandidateArr, rhoEleIso, kDataEra_NONE, 
                            info->nPV0, info->runNum, info->lumiSec, info->evtNum);
        }

      } //loop over electrons

    } //end loop over data  

    cout << "Total Electrons: " << NElectronsFilled << endl;

  } //end loop over files

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;


  cout << "Total Electrons: " << NElectronsFilled << endl;
  outputFile->Write();
  outputFile->Close();

  gBenchmark->Show("WWTemplate");       
} 


