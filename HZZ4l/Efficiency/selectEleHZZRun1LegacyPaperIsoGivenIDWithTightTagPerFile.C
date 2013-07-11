//================================================================================================
//
// Select probes for electron working point (ID + iso) efficiency with Tag&Probe method
//
//  * outputs ROOT file with a TTree of probes
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for 4-vector calculations
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
#include "HiggsAna/DataTree/interface/TMuon.hh"
#include "HiggsAna/DataTree/interface/TJet.hh"
#include "HiggsAna/DataTree/interface/TGenInfo.hh"
#include "HiggsAna/DataTree/interface/TGenParticle.hh"
#include "HiggsAna/DataTree/interface/TPFCandidate.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// Helper functions for Electron ID selection
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"
#include "HiggsAna/Utils/CommonDefs.hh"
#include "HiggsAna/Utils/LeptonIDCuts.hh"
#include "HiggsAna/HZZ4l/Utils/LeptonSelection.hh"
#include "HiggsAna/HZZ4l/Utils/HZZDefs.hh"
#include "HiggsAna/HZZ4l/Utils/FSRRecovery.hh"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

// structure for output ntuple
#include "CITCommon/CommonData/interface/EffData.hh" 
#endif


//=== MAIN MACRO ================================================================================================= 

void selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTagPerFile(const string inputfile,          // input file
                              const string outputfile,         // output directory
                              const Bool_t  matchGen = kFALSE, // match to generator muons
                              Int_t dataType = 0,              // del = 0, sel = 1, mc = -1
                              Int_t DataEraInput = 2
  ) {
  gBenchmark->Start("selectEleMVAIsoTightGivenID");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  UInt_t DataEra = kDataEra_NONE;
  if (DataEraInput == 1) DataEra = kDataEra_2011_MC;
  if (DataEraInput == 2) DataEra = kDataEra_2012_MC;
  if (DataEraInput == 11) DataEra = kDataEra_2011_Data;
  if (DataEraInput == 12) DataEra = kDataEra_2012_Data;

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
  

  // mass region
  Double_t massLo = 40;
  Double_t massHi = 200;

  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  Double_t nProbes = 0;
  
  //
  // Set up output ntuple
  //
  TFile *outFile = new TFile(outputfile.c_str(),"RECREATE"); 
  TTree *outTree = new TTree("Events","Events");
  EffData data;
  outTree->Branch("Events",&data.mass,"mass/F:pt:eta:phi:weight:q/I:npv/i:npu:pass:runNum:lumiSec:evtNum:rho/F");

  TFile *infile=0;
  TTree *eventTree=0;
  
  // Data structures to store info from TTrees
  higgsana::TEventInfo *info  = new higgsana::TEventInfo();
  higgsana::TGenInfo *gen     = new higgsana::TGenInfo();
  TClonesArray *electronArr = new TClonesArray("higgsana::TElectron");
  TClonesArray *pfcandidateArr = new TClonesArray("higgsana::TPFCandidate");
  TClonesArray *muonArr = new TClonesArray("higgsana::TMuon");

  // Read input file and get the TTrees
  cout << "Processing " << inputfile << "..." << endl;
  infile = TFile::Open(inputfile.c_str(),"read");
  assert(infile);



  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kFALSE;
  mithep::RunLumiRangeMap rlrm;
  if (!matchGen) {
    hasJSON = kTRUE;
    rlrm.AddJSONFile("/afs/cern.ch/work/s/sixie/public/HZZ4l/auxiliar/2012/Cert_Full2012_53X_JSON.txt"); 
    rlrm.AddJSONFile("/afs/cern.ch/work/s/sixie/public/HZZ4l/auxiliar/2011/2011Combined.json");
  }
    
  eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
  eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("Muon", &muonArr);         TBranch *muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr);  TBranch *pfcandidateBr = eventTree->GetBranch("PFCandidate");
  cout << "NEvents = " << eventTree->GetEntries() << endl;

  TBranch *genBr = 0;
  if(matchGen) {
    eventTree->SetBranchAddress("Gen", &gen);
    genBr = eventTree->GetBranch("Gen");
  }
    
  Double_t weight = 1;

  // loop over events
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    if (ientry % 100000 == 0) cout << "Processed Event " << ientry << endl;
    infoBr->GetEntry(ientry);

    // check for certified runs
    mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
    if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  
 
    //***********************************************************
    // Definition of Pileup Energy density
    //***********************************************************
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

    //use only odd numbered events to evaluate efficiency for data. even numbered events were used for training
    //if (info->evtNum % 2 == 0 && !matchGen) continue;


    // trigger requirement               
    Bool_t passTrigger = kFALSE;
    if(dataType == 0) {
      if ((info->triggerBits & kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50) == kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50) passTrigger = kTRUE;
      if ((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) == kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) passTrigger = kTRUE;
      if ((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) == kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) passTrigger = kTRUE;
    } else if(dataType == 1) {
      if(DataEraInput == 2) {
        if(info->triggerBits & kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50)  continue;
        if(info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30)  continue;
        if(info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17)                         continue;
      }
        
      if ((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) == kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) passTrigger = kTRUE;
    }
    if(dataType != -1 && !passTrigger) continue;     
      
    // good vertex requirement
    if(!(info->hasGoodPV)) continue;

    if(matchGen) genBr->GetEntry(ientry);
      
    electronArr->Clear();
    muonArr->Clear(); 
    pfcandidateArr->Clear(); 
    electronBr->GetEntry(ientry);
    muonBr->GetEntry(ientry);
    pfcandidateBr->GetEntry(ientry);

    //********************************************************
    //Low Met Requirement
    //********************************************************
    TVector3 met;        
    if(info->pfMEx!=0 || info->pfMEy!=0) {       
      met.SetXYZ(info->pfMEx, info->pfMEy, 0);
    }
    if (met.Pt() > 25) continue;
      

    //********************************************************
    //Loop over TAG electrons
    //********************************************************
    for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {

      const higgsana::TElectron *tag = (higgsana::TElectron*)((*electronArr)[i]);
	
      if(matchGen) {
        Bool_t match1 = (higgsana::deltaR(tag->eta, tag->phi, gen->eta_1, gen->phi_1) < 0.5);
        Bool_t match2 = (higgsana::deltaR(tag->eta, tag->phi, gen->eta_2, gen->phi_2) < 0.5);
        if(!match1 && !match2)
          continue;
      }

      if(tag->pt          < 20)  continue;
      if(fabs(tag->scEta) > 2.5) continue;

      if (!PassEleSimpleCutsVeryTight(tag,pfcandidateArr,rhoEleIso,EleEAEra)) continue;

      if(dataType == 0 &&
         !((info->triggerBits & kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50) && (tag->hltMatchBits & kHLTObject_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT)) &&
         !((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (tag->hltMatchBits & kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT)) &&
         !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (tag->hltMatchBits & kHLTObject_Ele32_CaloIdL_CaloIsoVL)) &&
         !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (tag->hltMatchBits & kHLTObject_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT)) &&
         !((info->triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL) && (tag->hltMatchBits & kHLTObject_Ele17_CaloIdL_CaloIsoVL)) ) 
        continue;
      
      if (dataType == 1) {
        if(dataType == 1 &&
           !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLTObject_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && tag->pt > 30) &&
           !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLTObject_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && tag->pt > 35) &&
           !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLTObject_Ele52_CaloIdVT_TrkIdT) && tag->pt > 60) &&
           !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLTObject_Ele27_WP80) && tag->pt > 30) && 
           !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLTObject_Ele32_WP70) && tag->pt > 35)  
          )
          continue;       
      }

      const Double_t m = 0.000511;
      TLorentzVector vtag;
      vtag.SetPtEtaPhiM(tag->pt, tag->eta, tag->phi, m);
        
      for(Int_t j=0; j<electronArr->GetEntriesFast(); j++) {
        if(i==j) continue;
	  
        const higgsana::TElectron *probe = (higgsana::TElectron*)((*electronArr)[j]);
        if(probe->q == tag->q) continue;
	  
// 	  if(typev[ifile]==eDiEl &&
// 	     !((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (probe->hltMatchBits & kHLTObject_SC8)) &&
// 	     !((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (probe->hltMatchBits & kHLTObject_Ele8)) &&
// 	     !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (probe->hltMatchBits & kHLTObject_SC17)))
// 	    continue;
	  
        if(matchGen) {
          Bool_t match1 = (higgsana::deltaR(probe->eta, probe->phi, gen->eta_1, gen->phi_1) < 0.5);
          Bool_t match2 = (higgsana::deltaR(probe->eta, probe->phi, gen->eta_2, gen->phi_2) < 0.5);
          if(!match1 && !match2)
            continue;
        }
	  
        if(fabs(probe->scEta) > 2.5) continue;

        TLorentzVector vprobe;
        vprobe.SetPtEtaPhiM(probe->pt, probe->eta, probe->phi, m);
	  
        TLorentzVector vdielectron = vtag + vprobe;
        if((vdielectron.M()<massLo) || (vdielectron.M()>massHi)) continue;	  	  

        //for probes with pT < 10, require the Ele20_SC4 trigger
        if (probe->pt < 10) {
          if(dataType == 0) {

            if (!((info->triggerBits & kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50) == kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50)) continue;
            if (!((info->triggerBits & kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50) && (tag->hltMatchBits & kHLTObject_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT))) continue;
          }
        }

        nProbes++;

        Bool_t passID = (PassEleHZZ4lPreselection(probe) && PassEleHZZ4lRun1LegacyPaperID(probe, eleIDMVA));
        
        vector<const higgsana::TPFCandidate*> photonsToVeto;
        Bool_t passIsolation = PassEleHZZ4lICHEP2012Iso(probe,j,pfcandidateArr,rhoEleIso,EleEAEra,photonsToVeto);
        if (!passID) continue;

        //******************
        //PASS
        //******************
        Bool_t pass = passID && passIsolation;

        // Fill tree
        data.mass   = vdielectron.M();
        data.pt     = probe->pt;
        data.eta    = probe->scEta;
        data.phi    = probe->phi;
        data.weight = weight;
        data.q      = probe->q;
        data.npv    = info->nPV0;
        data.npu    = info->nPUEvents;
        data.pass   = (pass) ? 1 : 0;
        data.rho    = rhoEleIso;
        outTree->Fill();	  
      }
    }
  }

  delete infile;
  infile=0, eventTree=0;      
  delete info;
  delete gen;
  delete electronArr;
  
     
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  cout << " Number of probes selected: " << nProbes << endl;
  
  outFile->Write();
  outFile->Close();
  delete outFile;
  
  cout << endl;
  cout << "  <> Output saved in " << outputfile << endl;    
  cout << endl;  
      
  gBenchmark->Show("selectEleLHEffTP"); 
}
