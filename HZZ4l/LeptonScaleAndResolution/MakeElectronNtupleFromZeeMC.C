//root -l -b -q EWKAna/Hww/FakeRate/ComputeElectronFakeRate_Data.C+\(\)
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
#include "HiggsAna/Utils/LeptonIDCuts.hh"
#include "HiggsAna/DataTree/interface/Types.h"
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"

#include "CITCommon/CommonData/interface/ElectronTree.h"
#include "HiggsAna/Utils/LeptonTools.hh"

#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

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
void MakeElectronNtupleFromZeeMC() {

  MakeNtuple("LIST","ElectronSelectionTraining.Real.Training.root");
  MakeNtuple("LIST","ElectronSelectionTraining.Real.Testing.root");

}


void MakeElectronNtupleFromZeeMC(const string inputFilename, const string outputFilename, Int_t Option) {

  MakeNtuple(inputFilename, outputFilename, Option);

}


void MakeNtuple(const string inputFilename, const string outputFilename, Int_t Option)
{  
  gBenchmark->Start("WWTemplate");

  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  TFile *outputFile = new TFile(outputFilename.c_str(), "RECREATE");
  citana::ElectronTree *eleTree = new citana::ElectronTree;
  eleTree->CreateTree();
  eleTree->tree_->SetAutoFlush(0);

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
  TClonesArray *genparticleArr = new TClonesArray("higgsana::TGenParticle");

  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("/afs/cern.ch/work/s/sixie/public/HZZ4l/auxiliar/2011/2011Combined.json"); 
  rlrm.AddJSONFile("/afs/cern.ch/work/s/sixie/public/HZZ4l/auxiliar/2012/Cert_Full2012_53X_JSON.txt");

  Int_t NEvents = 0;

  UInt_t DataEra = kDataEra_NONE;

  vector<string> inputfiles;
  if (inputFilename == "LIST") {
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-del-a05-v1.TightPlusReco.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-del-m10-v1.TightPlusReco.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-del-o03-v1.TightPlusReco.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-del-pr-v4.TightPlusReco.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11b-del-pr-v1.TightPlusReco.root");
  } else {
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
    TBranch *genparticleBr;


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
    eventTree->SetBranchAddress("GenParticle", &genparticleArr); genparticleBr = eventTree->GetBranch("GenParticle");

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      NEvents++;

      //********************************************************
      // Load the branches
      //********************************************************
      electronArr->Clear(); 
      muonArr->Clear(); 
      photonArr->Clear(); 
      jetArr->Clear(); 
      pfcandidateArr->Clear(); 
      genparticleArr->Clear(); 
      electronBr->GetEntry(ientry);
      muonBr->GetEntry(ientry);
      photonBr->GetEntry(ientry);
      jetBr->GetEntry(ientry);
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


      Int_t NElectrons = electronArr->GetEntries();

      //dilepton preselection
      if (NElectrons < 2) continue;


      //******************************************************************************
      //loop over electron pairs
      //******************************************************************************

      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const higgsana::TElectron *tag = (higgsana::TElectron*)((*electronArr)[i]);

        //********************************************************
        // Select MC Truth Electrons
        //********************************************************              
        if (!((UInt_t(abs(max(0,tag->isMCReal))) & 1) == 1)) continue;

        //Tighter cuts on the tag to reduce bkg
	if(tag->pt          < 20)  continue;
	if(fabs(tag->eta) > 2.5) continue;
        if (!passCutBasedEleID(tag,ComputeElePFIso04(tag,pfcandidateArr,rhoEleIso,EleEAEra))) continue;

        for(Int_t j=0; j<electronArr->GetEntries(); j++) {
	  if(i==j) continue;

          const higgsana::TElectron *probe = (higgsana::TElectron*)((*electronArr)[j]);

          //********************************************************
          // Select MC Truth Electrons
          //********************************************************              
          if (!((UInt_t(abs(max(0,probe->isMCReal))) & 1) == 1)) continue;

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
              double tmpDR = higgsana::deltaR( probe->eta, probe->phi, gen->eta , gen->phi);
              if (tmpDR < minDR && tmpDR < 0.3) {
                minDR = tmpDR;
                matchedGenElectron = gen;
              }
            }
          }
          if (!matchedGenElectron || abs(matchedGenElectron->motherPdgID) == 15) continue;



          if(probe->q == tag->q) continue;
	  if(fabs(probe->eta) > 2.5) continue;	  

          higgsana::FourVectorM lepton1;
          higgsana::FourVectorM lepton2;
          lepton1.SetCoordinates(tag->pt, tag->eta, tag->phi, 0.51099892e-3 );
          lepton2.SetCoordinates(probe->pt, probe->eta, probe->phi, 0.51099892e-3 );
          higgsana::FourVectorM dilepton = lepton1+lepton2;

          //select Z peak
          if (dilepton.M() > 75.0 && dilepton.M() < 105.0) {

            if (fabs(probe->dz) >= 0.1) continue;

            //Fill Probe Electron
            NElectronsFilled++;
            FillElectronTree( eleTree, probe, j, pfcandidateArr, eleIDMVA, genparticleArr, rhoEleIso, EleEAEra, 
                              info->nPV0, info->runNum, info->lumiSec, info->evtNum, dilepton.M());
            
          }
        }
      }

    }

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


