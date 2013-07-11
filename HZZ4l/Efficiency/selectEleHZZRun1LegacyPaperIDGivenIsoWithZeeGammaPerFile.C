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


Bool_t MatchedToStatus1Ele( const higgsana::TElectron *ele , TClonesArray *genparticleArr) {
  Bool_t matched = kFALSE;
  for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
    const higgsana::TGenParticle *gen = (higgsana::TGenParticle*)((*genparticleArr)[k]);

    //match to status 1 or status 3
    if (abs(gen->pdgid) == 11 && (gen->status == 1)) {
      if ( higgsana::deltaR( ele->eta, ele->phi, gen->eta , gen->phi) < 0.3) {
        matched = kTRUE;
        break;
      }
    }
  }
  return matched;
}


Bool_t MatchedToStatus1Photon( const higgsana::TPhoton *pho , TClonesArray *genparticleArr) {
  Bool_t matched = kFALSE;
  for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
    const higgsana::TGenParticle *gen = (higgsana::TGenParticle*)((*genparticleArr)[k]);

    //match to status 1 or status 3
    if (abs(gen->pdgid) == 22) {
//       cout << "gen pho " << gen->status << " " << gen->pt << " " << gen->eta << " " << gen->phi << " : " << pho->et << " " << pho->eta << " " << pho->phi << endl;
    }

    if (abs(gen->pdgid) == 22 && (gen->status == 1)) {
      if ( higgsana::deltaR( pho->eta, pho->phi, gen->eta , gen->phi) < 0.3) {
        matched = kTRUE;
        break;
      }
    }
  }
//   cout << "matched " << matched << endl;
  return matched;
}




//=== MAIN MACRO ================================================================================================= 

void selectEleHZZRun1LegacyPaperIDGivenIsoWithZeeGammaPerFile(const string inputfile,          // input file
                                                              const string outputfile,         // output directory
                                                              const Bool_t  matchGen = kFALSE, // match to generator muons
                                                              Int_t dataType = 0,              // del = 0, sel = 1, mc sig = -1, mc bkg = -2
                                                              Int_t DataEraInput = 2
  ) {

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
  TClonesArray *genparticleArr = new TClonesArray("higgsana::TGenParticle");
  TClonesArray *electronArr = new TClonesArray("higgsana::TElectron");
  TClonesArray *photonArr = new TClonesArray("higgsana::TPhoton");
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
  }
    
  eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
  eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("Photon", &photonArr); TBranch *photonBr = eventTree->GetBranch("Photon");
  eventTree->SetBranchAddress("Muon", &muonArr);         TBranch *muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr);  TBranch *pfcandidateBr = eventTree->GetBranch("PFCandidate");
  cout << "NEvents = " << eventTree->GetEntries() << endl;

  TBranch *genparticleBr;
  if(matchGen) {
    eventTree->SetBranchAddress("GenParticle", &genparticleArr);
    genparticleBr = eventTree->GetBranch("GenParticle");
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
    //Don't need this for zee gamma
    //if (info->evtNum % 2 == 0 && !matchGen) continue;

    // trigger requirement               
    Bool_t passTrigger = kFALSE;
    if(dataType == 0 ) {
      if ((info->triggerBits & kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50) == kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50) passTrigger = kTRUE;
      if ((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) == kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) passTrigger = kTRUE;
      if ((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) == kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) passTrigger = kTRUE;
    } else if(dataType == 1 ) {
      if(info->triggerBits & kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50)  continue;
      if(info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30)  continue;
      if(info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17)                         continue;
 
      if ((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) == kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) passTrigger = kTRUE;
    } else if (dataType < 0) {
      if ((info->triggerBits & kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50) == kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50) passTrigger = kTRUE;
      if ((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) == kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) passTrigger = kTRUE;
      if ((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) == kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) passTrigger = kTRUE;
      if ((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) == kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) passTrigger = kTRUE;
    }
    if(!passTrigger) continue;     



    // good vertex requirement
    if(!(info->hasGoodPV)) continue;

    electronArr->Clear();
    muonArr->Clear(); 
    pfcandidateArr->Clear(); 
    genparticleArr->Clear(); 
    electronBr->GetEntry(ientry);
    photonBr->GetEntry(ientry);
    muonBr->GetEntry(ientry);
    pfcandidateBr->GetEntry(ientry);
    if(matchGen) {
      genparticleBr->GetEntry(ientry);
    }  


    //********************************************************
    //Loop over TAG electrons
    //********************************************************
    vector<Int_t> probeAlreadyUsed;
    for(Int_t i=0; i<electronArr->GetEntriesFast(); ++i) {
      probeAlreadyUsed.push_back(kFALSE);
    }

    for(Int_t i=0; i<electronArr->GetEntriesFast(); ++i) {

      const higgsana::TElectron *tag = (higgsana::TElectron*)((*electronArr)[i]);
	
      Bool_t TagIsEle = MatchedToStatus1Ele(tag, genparticleArr);

      if(tag->pt          < 20)  continue;
      if(fabs(tag->scEta) > 2.5) continue;
//       if(dataType == 1) {
//         if(tag->pt          < 30)  continue;
//       }

      if (!passCutBasedEleID(tag, ComputeElePFIso04(tag,pfcandidateArr,rhoEleIso,EleEAEra))) continue;

      if(dataType == 0 &&
         !((info->triggerBits & kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50) && (tag->hltMatchBits & kHLTObject_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT)) &&
         !((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (tag->hltMatchBits & kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT)) &&
         !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (tag->hltMatchBits & kHLTObject_Ele32_CaloIdL_CaloIsoVL)) &&
         !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (tag->hltMatchBits & kHLTObject_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT))         
        )
        continue;
      
      if(dataType == 1 &&
         !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLTObject_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT)) &&
         !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLTObject_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT)) &&
         !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLTObject_Ele52_CaloIdVT_TrkIdT)) &&
         !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLTObject_Ele65_CaloIdVT_TrkIdT)) &&
         !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLTObject_Ele80_CaloIdVT_TrkIdT)) &&
         !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLTObject_Ele27_WP80)) &&
         !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLTObject_Ele32_WP70))
        ) 
        continue;
        

      const Double_t m = 0.000511;
      TLorentzVector vtag;
      vtag.SetPtEtaPhiM(tag->pt, tag->eta, tag->phi, m);
        
      //Find Photon
      for(Int_t j=0; j<photonArr->GetEntriesFast(); ++j) {
      
        const higgsana::TPhoton *pho = (higgsana::TPhoton*)((*photonArr)[j]);
        TLorentzVector vpho;
        vpho.SetPtEtaPhiM(pho->et, pho->eta, pho->phi, 0);
        
        if (pho->et < 10) continue;
        if (fabs(pho->eta) > 2.5) continue;
    
        if (!passPhotonSimpleCuts(pho)) continue;        

        Bool_t PhotonIsReal = MatchedToStatus1Photon(pho, genparticleArr);

        if (fabs(tag->eta - pho->eta) < 0.15 && higgsana::deltaR( tag->eta, tag->phi, pho->eta, pho->phi) < 0.7) continue;                

        //Find Probes
        for(Int_t k=0; k<electronArr->GetEntriesFast(); ++k) {

          //no duplicates
          if(k==i) continue;    
          if (probeAlreadyUsed[k]) continue;

          const higgsana::TElectron *probe = (higgsana::TElectron*)((*electronArr)[k]);
          if(probe->q == tag->q) continue;
          TLorentzVector vprobe;
          vprobe.SetPtEtaPhiM(probe->pt, probe->eta, probe->phi, m);

          Bool_t ProbeIsEle = MatchedToStatus1Ele(probe, genparticleArr);


          //High Efficiency Pixel Veto
          double minDEta = fmin( fabs(tag->eta - pho->eta) , fabs(probe->eta - pho->eta ));
          double minDPhi = fmin( higgsana::deltaPhi(tag->phi, pho->phi) , higgsana::deltaPhi( probe->phi, pho->phi));
          double minDR = fmin( higgsana::deltaR(tag->eta,tag->phi, pho->eta, pho->phi) , higgsana::deltaR( probe->eta, probe->phi, pho->eta, pho->phi));
          Bool_t passPixelVeto = kTRUE;
          if (fabs(pho->eta) < 1.479) {
            if (minDEta > 0.04 || minDPhi > 0.3) {
              if (pho->hasPixelSeed) passPixelVeto = kFALSE;
            }
          } else {
            if (minDEta > 0.08 || minDPhi > 0.3) {
              if (pho->hasPixelSeed) passPixelVeto = kFALSE;
            }
          }

          //more eegamma selection cuts          
          if (fabs( probe->eta - pho->eta) < 0.15 && higgsana::deltaR( probe->eta, probe->phi, pho->eta, pho->phi) < 0.7) continue;

          //Selection cuts
          if (!((vtag + vpho + vprobe).M() + (vtag+vprobe).M() < 180)) continue;
          if (!((vtag+vprobe).M() > 40)) continue;
          if (!(higgsana::deltaR(probe->eta, probe->phi, pho->eta, pho->phi) < 1.5)) continue;
          if (!passPixelVeto) continue;

          //cut to emulate trigger.
          if (!((vtag+vpho).M() > 50)) continue;


          //optional cuts for tighter selection
          if (dataType == 0 || dataType == 1) {
            if (!(higgsana::deltaR(probe->eta, probe->phi, pho->eta, pho->phi) < 0.7)) continue;
          }
            
          //gen matching for MC
          if (matchGen) {
            if (dataType == -1) {
              if (!(TagIsEle && PhotonIsReal && ProbeIsEle)) continue;
              if (!(higgsana::deltaR(probe->eta, probe->phi, pho->eta, pho->phi) < 0.7)) continue;
            }
            if (dataType == -2) {
              if ((TagIsEle && PhotonIsReal && ProbeIsEle)) continue;
            }
          }

          TLorentzVector vdielectrongamma = vtag + vpho + vprobe;
          if((vdielectrongamma.M()<massLo) || (vdielectrongamma.M()>massHi)) continue;	  	  
 

          //find pfphotons to remove the photon footprint
          vector<const higgsana::TPFCandidate*> photonsToVeto;

          Double_t pfphotonMinDR = 9999;
          const higgsana::TPFCandidate *photonPFCandidate = 0;
          for( uint p=0; p< uint(pfcandidateArr->GetEntries()); ++p ) {
            const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*pfcandidateArr)[p]);
            if (higgsana::deltaR(pho->eta,pho->phi,pf->eta,pf->phi) < pfphotonMinDR) {
              photonPFCandidate = pf;
              pfphotonMinDR = higgsana::deltaR(pho->eta,pho->phi,pf->eta,pf->phi);
            }
          }
          if (pfphotonMinDR < 0.1) photonsToVeto.push_back(photonPFCandidate);

          for( uint p=0; p< uint(pfcandidateArr->GetEntries()); ++p ) {
            const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*pfcandidateArr)[p]);
            if (pf->pfType != eGamma) continue;
            if (fabs(pho->eta) < 1.479) {
              if (fabs(pho->eta - pf->eta) < 0.015) {
                photonsToVeto.push_back(pf);
              }
            } else {
              if (higgsana::deltaR(pho->eta,pho->phi,pf->eta,pf->phi) < 0.07) {
                photonsToVeto.push_back(pf);
              }
            }
          }


          //Fill probe
          nProbes++;
          Bool_t passID = (PassEleHZZ4lPreselection(probe) && PassEleHZZ4lRun1LegacyPaperID(probe, eleIDMVA));
          Bool_t passIsolation = PassEleHZZ4lICHEP2012Iso(probe,j,pfcandidateArr,rhoEleIso,EleEAEra,photonsToVeto);
          if (!passIsolation) continue;

          //debug printout
          cout << "RunLumiEvent: " << info->runNum << " " <<  info->lumiSec << " " << info->evtNum << " : " << probe->pt << " " << probe->eta << " " << probe->phi << endl;
          Bool_t test = PassEleHZZ4lICHEP2012ID(probe, eleIDMVA, true);
          cout << "pass : " << PassEleHZZ4lPreselection(probe) << " , " << test << endl;

          //******************
          //PASS
          //******************
          Bool_t pass = passID && passIsolation;

          // Fill tree
          data.mass   = vdielectrongamma.M();
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
          probeAlreadyUsed[k] = kTRUE;
          
        }


      }


    }
  }

  delete infile;
  infile=0, eventTree=0;      
  delete info;
  delete genparticleArr;
  delete electronArr;
  delete photonArr;
  delete muonArr;
  
     
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
