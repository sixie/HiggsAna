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



// print event dump
void MakeNtuple(const string inputFilename,  const string outputFilename, Int_t dataType = 0);

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
void MakeElectronNtupleFromZeeGammaMC() {

  MakeNtuple("LIST","ElectronSelectionTraining.Real.Training.root");
  MakeNtuple("LIST","ElectronSelectionTraining.Real.Testing.root");

}


void MakeElectronNtupleFromZeeGammaMC(const string inputFilename, const string outputFilename, Int_t dataType) {

  MakeNtuple(inputFilename, outputFilename, dataType );

}


void MakeNtuple(const string inputFilename, const string outputFilename, Int_t dataType)
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
  TClonesArray *genparticleArr = new TClonesArray("higgsana::TGenParticle");
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
      pfcandidateArr->Clear(); 
      genparticleArr->Clear(); 
      electronBr->GetEntry(ientry);
      muonBr->GetEntry(ientry);
      photonBr->GetEntry(ientry);
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
      // trigger requirement               
      //********************************************************
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



      Int_t NElectrons = electronArr->GetEntries();

      //dilepton preselection
      if (NElectrons < 2) continue;


      //********************************************************
      //Loop over TAG electrons
      //********************************************************
      vector<Int_t> probeAlreadyUsed;
      for(Int_t i=0; i<electronArr->GetEntriesFast(); ++i) {
        probeAlreadyUsed.push_back(kFALSE);
      }

      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const higgsana::TElectron *tag = (higgsana::TElectron*)((*electronArr)[i]);

        Bool_t TagIsEle = MatchedToStatus1Ele(tag, genparticleArr);
        if (!TagIsEle) continue;

        //Tighter cuts on the tag to reduce bkg
	if(tag->pt          < 20)  continue;
	if(fabs(tag->eta) > 2.5) continue;
//       if(dataType == 1) {
//         if(tag->pt          < 30)  continue;
//       }

        if (!passCutBasedTightEleID(tag,ComputeElePFIso04(tag,pfcandidateArr,rhoEleIso,EleEAEra))) continue;

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
          if (!PhotonIsReal) continue;

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
            if (!ProbeIsEle) continue;

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
            
            TLorentzVector vdielectrongamma = vtag + vpho + vprobe;
            if((vdielectrongamma.M()< 60) || (vdielectrongamma.M()>120)) continue;

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

            NElectronsFilled++;
            FillElectronTree( eleTree, probe, j, pfcandidateArr, eleIDMVA, genparticleArr, rhoEleIso, EleEAEra, 
                              info->nPV0, info->runNum, info->lumiSec, info->evtNum, vdielectrongamma.M());
  
          } //loop over probes
        } //loop over photons
      } //loop over tags
    } //end loop over events

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


