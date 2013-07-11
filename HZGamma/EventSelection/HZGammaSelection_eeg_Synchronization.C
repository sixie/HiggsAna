//root -l EWKAna/HZZ4l/Selection/HZZ4lSelection.C+\(\"\"\)

//================================================================================================
//
// HZZ4l selection macro
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
#include "HiggsAna/DataTree/interface/TMuon.hh"
#include "HiggsAna/DataTree/interface/TJet.hh"
#include "HiggsAna/DataTree/interface/TGenParticle.hh"
#include "HiggsAna/DataTree/interface/TPFCandidate.hh"

// #include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"
#include "HiggsAna/Utils/CommonDefs.hh"
#include "HiggsAna/Utils/LeptonIDCuts.hh"
#include "HiggsAna/HZGamma/Utils/LeptonSelection.hh"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// output data structs
// #include "HiggsAna/HZZ4l/interface/HZZKinematics.hh"
// #include "HiggsAna/HZZ4l/interface/HZZGenInfo.hh"
// #include "HiggsAna/HZZ4l/interface/HZZ4lDefs.hh"

// #include "TMVAGui.C"
// #include "TMVA/Tools.h"
// #include "TMVA/Reader.h"
// #include "TMVA/MethodCuts.h"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================

//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
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

void HZGammaSelection_eeg_Synchronization(Int_t Option = 0) 
{  
  gBenchmark->Start("HZZTemplate");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  Bool_t printDebug = kFALSE;
  //Double_t lumi = 1092;              // luminosity (pb^-1)
  Bool_t doFSR = kTRUE;

  Bool_t isData = kFALSE;
  UInt_t DataEra = kDataEra_2012_MC;

//   if (Option == 0) {
//     DataEra = kDataEra_2011_MC;
//   } else if (Option == 1) {
//     DataEra = kDataEra_2012_MC;
//   }

  vector<vector<string> > inputFiles;

  inputFiles.push_back(vector<string>());
  if (Option == 0) {
    inputFiles.back().push_back("/afs/cern.ch/work/s/sixie/public/HZGamma/Synchronization/BACON/eeg/s12-zllm50-2-v9_002C5B35-519B-E111-862D-001E67398025.root");
  } 
  if (Option == 1) {
    inputFiles.back().push_back("/afs/cern.ch/work/s/sixie/public/HZGamma/Synchronization/BACON/mumug/BACON_s12-zllm50-2-v9_EAF43999-8D9B-E111-A418-003048D4610E.root");
  }

  vector<string> processNames;
  processNames.push_back("HZZ120");
//   processNames.push_back("Data");
//   processNames.push_back("ZZ");
//   processNames.push_back("WZ");
//   processNames.push_back("ttbar");
//   processNames.push_back("WW"); 

  assert(processNames.size() == inputFiles.size());


  //--------------------------------------------------------------------------------------------------------------
  // Yields
  //==============================================================================================================  
  const UInt_t NSelectionStages = 12;
  Int_t NEventsTotal = 0;
  Int_t NEventsPassSelectionStages[NSelectionStages];
  for (UInt_t n=0; n < NSelectionStages; ++n) {
    NEventsPassSelectionStages[n] = 0;
  }


  
  //--------------------------------------------------------------------------------------------------------------
  // output ntuple structure
  //==============================================================================================================  
//   HZZKinematics kinematics;
//   HZZGenInfo    geninfo;


  //--------------------------------------------------------------------------------------------------------------
  // Pileup Reweighting
  //==============================================================================================================  
//   TFile *fPUFile = TFile::Open("/data/smurf/sixie/Pileup/weights/PileupReweighting.Fall11_To_Full2011.root");
//   TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
//   assert(fhDPU);
//   fhDPU->SetDirectory(0);
//   delete fPUFile;


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  ofstream eventList("ReferenceSelection_sixie.eventList.txt");
  ofstream eleList("ElectronDump_sixie.txt");
  ofstream phoList("PhotonDump_sixie.txt");
  //
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  higgsana::TEventInfo *info    = new higgsana::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("higgsana::TElectron");
  TClonesArray *muonArr = new TClonesArray("higgsana::TMuon");
  TClonesArray *jetArr = new TClonesArray("higgsana::TJet");
  TClonesArray *photonArr = new TClonesArray("higgsana::TPhoton");
  TClonesArray *genparticleArr = new TClonesArray("higgsana::TGenParticle");
  TClonesArray *pfcandidateArr = new TClonesArray("higgsana::TPFCandidate");
  
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  mithep::RunLumiRangeMap rlrm;
  if (isData) {
     rlrm.AddJSONFile("/afs/cern.ch/work/s/sixie/public/HZZ4l/auxiliar/2012/Cert_190456-196531_8TeV_PromptReco_Collisions12_JSON.txt"); 
  }

  for (UInt_t qq = 0; qq<inputFiles.size() ; ++qq) { 
    for (UInt_t f = 0; f < inputFiles[qq].size() ; ++f) {
      //********************************************************
      // Get Tree
      //********************************************************
      cout << "Reading File " << inputFiles[qq][f] << endl;
      eventTree = getTreeFromFile(inputFiles[qq][f].c_str(),"Events"); 
      TBranch *infoBr;
      TBranch *electronBr;
      TBranch *muonBr;
      TBranch *jetBr;
      TBranch *photonBr;
      //TBranch *genparticleBr;
      TBranch *pfcandidateBr;


      //*****************************************************************************************
      //Loop over muon Data Tree
      //*****************************************************************************************
      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");
      eventTree->SetBranchAddress("Photon", &photonArr);     photonBr = eventTree->GetBranch("Photon");
//       eventTree->SetBranchAddress("GenParticle", &genparticleArr);         genparticleBr = eventTree->GetBranch("GenParticle");
      eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr);         pfcandidateBr = eventTree->GetBranch("PFCandidate");
  
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
        printDebug = kFALSE;

        infoBr->GetEntry(ientry);
        if (ientry % 10000 == 0) cout << "Event " << ientry << endl;
	
        //Use Only Testing events
        //if (info->evtNum % 2 == 0 ) continue;
        //Use Only Training events
        //if (info->evtNum % 2 != 0 ) continue;
       

        //********************************************************
        //double mynpu = TMath::Min((double)info->nPUEvents,34.999);
        //Int_t npuxbin = fhDPU->GetXaxis()->FindBin(mynpu);
        //double npuWeight = fhDPU->GetBinContent(npuxbin);
        //********************************************************


        //***********************************************************
        // Definition of Pileup Energy density
        //***********************************************************

        
        Double_t rhoMuonIso = 0;
        Double_t rhoEleIso = 0;
        Double_t rhoPhoIso = 0;
        UInt_t MuonEAEra = 0;
        UInt_t EleEAEra = 0;
        UInt_t PhoEAEra = 0;

        if (DataEra == kDataEra_2011_MC) {
     
          if (!(isnan(info->RhoKt6PFJetsForIso25) || 
                isinf(info->RhoKt6PFJetsForIso25))) {
            rhoMuonIso = info->RhoKt6PFJetsForIso25;
            rhoEleIso = info->RhoKt6PFJetsForIso25;
            rhoPhoIso = info->RhoKt6PFJetsForIso25;
          }
          MuonEAEra = kDataEra_2011_Data;
          EleEAEra = kDataEra_2011_Data;
          PhoEAEra = kDataEra_2011_Data;
        } else if (DataEra == kDataEra_2012_MC) {

          if (!(isnan(info->RhoKt6PFJetsCentralNeutral) || 
                isinf(info->RhoKt6PFJetsCentralNeutral))) {
            rhoMuonIso = info->RhoKt6PFJetsCentralNeutral;
          }

          if (!(isnan(info->RhoKt6PFJets) || 
                isinf(info->RhoKt6PFJets))) {
            rhoEleIso = info->RhoKt6PFJets;
            rhoPhoIso = info->RhoKt6PFJets;
          }

          MuonEAEra = kDataEra_2012_Data;
          EleEAEra = kDataEra_2012_Data;
          PhoEAEra = kDataEra_2012_Data;
        }

        if(printDebug) cout << "Muon Isolation Rho: " << rhoMuonIso << endl;
        if(printDebug) cout << "Muon EA Era: " << MuonEAEra << endl;
        if(printDebug) cout << "Ele Isolation Rho: " << rhoEleIso << endl;
        if(printDebug) cout << "Ele EA Era: " << EleEAEra << endl;
        if(printDebug) cout << "Pho Isolation Rho: " << rhoPhoIso << endl;
        if(printDebug) cout << "Pho EA Era: " << PhoEAEra << endl;

   

        //********************************************************
        // Printdebug
        //********************************************************
        if ((0 == 1) 
             || (info->evtNum == 64036294)   
          ) printDebug = kTRUE;
        


        //********************************************************
        // Load the branches
        //********************************************************
        electronArr->Clear(); 
        muonArr->Clear(); 
        jetArr->Clear();  
       
        genparticleArr->Clear(); 
        pfcandidateArr->Clear(); 
        electronBr->GetEntry(ientry);
        muonBr->GetEntry(ientry);
        jetBr->GetEntry(ientry);
        photonBr->GetEntry(ientry);
//         genparticleBr->GetEntry(ientry);
        pfcandidateBr->GetEntry(ientry);


        if (printDebug) {
          cout << endl;
          cout << "****************************************************************************\n";
          cout << "Event " << info->runNum << " " << info->lumiSec << " " << info->evtNum << " ";
        }


        //********************************************************
        // Met
        //********************************************************
        TVector3 pfMet;        
        if(info->pfMEx!=0 || info->pfMEy!=0) {       
          pfMet.SetXYZ(info->pfMEx, info->pfMEy, 0);
        }



        //********************************************************
        // Selection Steps
        //********************************************************
        Bool_t PassSelectionStages[NSelectionStages];
        for (UInt_t n=0; n < NSelectionStages; ++n) {
          PassSelectionStages[n] = kFALSE;
        }

        //***************************************************************
        // For data apply JSON
        //***************************************************************
        // check for certified runs
//         mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
//         if(isData && !rlrm.HasRunLumi(rl)) continue;  
        

        //***************************************************************
        // Selection Stage 0 : 2 reco electrons, pt > 20, 10, pass JSON
        //***************************************************************
        int nlep_above_20=0;
        int nlep_above_10=0;
//         for(uint i=0; i<uint(muonArr->GetEntries()); i++) 
//         {
//           const higgsana::TMuon *mu = (higgsana::TMuon*)((*muonArr)[i]);        
//           if( !( (mu->typeBits & kGlobal) == kGlobal || 
//                  (mu->typeBits & kTracker) == kTracker) 
//             ) continue;
//           if( fabs(mu->eta) > 2.4 ) continue; 
//           if( mu->pt > 10 )  nlep_above_10++;
//           if( mu->pt > 20 )  nlep_above_20++;
//         }
        for(uint i=0; i<uint(electronArr->GetEntries()); i++) 
        {
          const higgsana::TElectron *ele = (higgsana::TElectron*)((*electronArr)[i]);        
//           if( fabs(ele->eta) > 2.5 ) continue; 
          if( ele->pt > 10 )  nlep_above_10++;
          if( ele->pt > 20 )  nlep_above_20++;
        }
        if( nlep_above_10 >= 2 && nlep_above_20 >= 1) {
          PassSelectionStages[0] = kTRUE;
        } else {
          PassSelectionStages[0] = kFALSE;                   
        }

        //***************************************************************
        // Selection Stage 1 : Vertex found
        //***************************************************************
        //I think my ntuples always require vertex
        PassSelectionStages[1] = kTRUE;
        

        //***************************************************************
        // Selection Stage 2 : HLT selection
        //***************************************************************
        if (DataEra == kDataEra_2011_MC) {
          if (
            (info->triggerBits & kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL) == kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
            ||
            (info->triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL) == kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL
            ||
            (info->triggerBits & kHLT_DoubleMu7) == kHLT_DoubleMu7
            ||
            (info->triggerBits & kHLT_Mu13_Mu8) == kHLT_Mu13_Mu8
            ||
            (info->triggerBits & kHLT_Mu17_Mu8) == kHLT_Mu17_Mu8
            ||
            (info->triggerBits & kHLT_IsoMu24) == kHLT_IsoMu24
            ) {
            PassSelectionStages[2] = kTRUE;
          } else {
            PassSelectionStages[2] = kFALSE;
          }

          if (printDebug) {
            cout << "DEBUG Trigger: "
                 << Bool_t((info->triggerBits & kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL) == kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL) << " "
                 << Bool_t((info->triggerBits & kHLT_Mu17_Mu8) == kHLT_Mu17_Mu8) << " "
                 << endl;
            cout << PassSelectionStages[2] << endl;
          }   


        } else if (DataEra == kDataEra_2012_MC) {
          if (
            (info->triggerBits & kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL) == kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
            ||
            (info->triggerBits & kHLT_Mu17_Mu8) == kHLT_Mu17_Mu8
            ||
            (info->triggerBits & kHLT_Mu17_TkMu8) == kHLT_Mu17_TkMu8
            ) {
            PassSelectionStages[2] = kTRUE;
          } else {
            PassSelectionStages[2] = kFALSE;
          }

          if (printDebug) {
            cout << "DEBUG Trigger: "
                 << Bool_t((info->triggerBits & kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL) == kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL) << " "
                 << Bool_t((info->triggerBits & kHLT_Mu17_Mu8) == kHLT_Mu17_Mu8) << " "
                 << Bool_t((info->triggerBits & kHLT_Mu17_TkMu8) == kHLT_Mu17_TkMu8) << " "
                 << Bool_t((info->triggerBits & kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL) == kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL) << " "
                  << Bool_t((info->triggerBits & kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL) == kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL) << " "
                << endl;
            cout << PassSelectionStages[2] << endl;
          }
        }
        if (printDebug) cout << "PassSelectionStages[2] = " << PassSelectionStages[2] << endl;


        //***********************************************************
        // Lepton Selection
        //***********************************************************
        vector<Int_t> leptonType;
        vector<UInt_t> leptonIndex;
        Int_t NLepPassID = 0;
        Int_t NLepPassIDIso = 0;
//         //    
//         if( printDebug ) cout << "\tnMuons: " << muonArr->GetEntries() << endl;
//         //----------------------------------------------------
//         for(uint i=0; i<uint(muonArr->GetEntries()); i++) 
//         {
//           const higgsana::TMuon *mu = (higgsana::TMuon*)((*muonArr)[i]);      

//           if (!((mu->typeBits & kGlobal) == kGlobal)) continue;

//           if( printDebug ) { 
//             cout << "muon:: " << mu->pt << " " << mu->eta << " " << mu->phi << " : " ;
//           }

//           if (! ( mu->pt > 7 && fabs(mu->eta) < 2.4)) {
//             if( printDebug )cout << endl;
//             continue;
//           }
//           if( printDebug )cout << " pass preselection : ";
          
//           if (!PassMuonHZGammaID(mu, (DataEra == kDataEra_2011_MC) )) {
//             if( printDebug )cout << endl;
//             continue;
//           } 
//           if( printDebug )cout << " pass ID : ";
//           NLepPassID++;      

//           if (! PassMuonHZGammaIso(mu, pfcandidateArr, rhoMuonIso,MuonEAEra, printDebug)) {
//             continue; 
//             if( printDebug )cout << " fail Iso : \n";
//           }
//           NLepPassIDIso++;      
//           if( printDebug )cout << " pass Iso : \n";


//           leptonType.push_back(13);
//           leptonIndex.push_back(i);
//         }
    
        //
        if( printDebug ) { cout << "\tnElectron: " << electronArr->GetEntries() << endl; }
        // --------------------------------------------------------------------------------
        for(uint i=0; i<uint(electronArr->GetEntries()); i++) 
        {
          const higgsana::TElectron *ele = (higgsana::TElectron*)((*electronArr)[i]);

          if( printDebug ) { 
            cout << "ele:: " << ele->pt << " " << ele->eta << " " << ele->phi << " : " ;
          }

          //Clean ID'ed and global muons from electrons if they fall within dR < 0.05
          Bool_t overlapsWithMuon = kFALSE;
          for(uint j=0; j<uint(muonArr->GetEntries()); j++) {
            const higgsana::TMuon *tmpmu = (higgsana::TMuon*)((*muonArr)[j]);
            if(!((tmpmu->typeBits & kGlobal) == kGlobal)) continue;
            if(!(fabs(tmpmu->eta) < 2.4)) continue;

            if (higgsana::deltaR(tmpmu->eta,tmpmu->phi, ele->eta,ele->phi) < 0.05) {
              overlapsWithMuon = kTRUE;
              break;
            }
          }
          if (overlapsWithMuon) {
            if( printDebug ) cout << " : Overlaps with Muon : ";
            if( printDebug ) cout << endl;              
            continue;
          }

          //preselection
          if (!( ele->pt > 7 && fabs(ele->eta) < 2.5)) {
            continue;
          }

          if (PassSelectionStages[0] && PassSelectionStages[1] && PassSelectionStages[2]) {
            eleList << info->evtNum << " " << ele->pt << " " << ele->eta << " " << PassEleHZGammaID(ele) << " " << PassEleHZGammaIso(ele, i, pfcandidateArr, rhoEleIso,EleEAEra, printDebug) << "\n";
          }


          if (!PassEleHZGammaID(ele, printDebug)) {
            if( printDebug )cout << " Fail ID " << endl;
            continue;
          }
          if( printDebug )cout << " pass ID : ";
          NLepPassID++;      

          if (! PassEleHZGammaIso(ele, i ,pfcandidateArr, rhoEleIso,EleEAEra, printDebug)) {
          if( printDebug )cout << " fail Iso : \n ";
            continue; 
          }
          NLepPassIDIso++;      
          if( printDebug )cout << " pass Iso : \n";
   
          leptonType.push_back(11);
          leptonIndex.push_back(i);
        }

        if( printDebug ) cout << "\n NLepPassID " << NLepPassID << " : NLepPassIDIso " << NLepPassIDIso << endl;

        //***********************************************************
        // Selection Stage 3: Lepton Preselection + ID
        //***********************************************************
        if (NLepPassID >= 2) {
          PassSelectionStages[3] = kTRUE;
        } else {
          PassSelectionStages[3] = kFALSE;
        }

        //***********************************************************
        // Selection Stage 4: Lepton Preselection + ID + Iso
        //***********************************************************
        if (NLepPassIDIso >= 2) {
          PassSelectionStages[4] = kTRUE;
        } else {
          PassSelectionStages[4] = kFALSE;
        }


        //***********************************************************
        // Photon Selection
        //***********************************************************        
        vector<UInt_t> photonIndex;
        Int_t NPhoPassSelection = 0;
        Int_t NPhoPassDRCut = 0;

        //    
        if( printDebug ) cout << "\n\tnPhotons: " << photonArr->GetEntries() << endl;
        //----------------------------------------------------
        for(uint i=0; i<uint(photonArr->GetEntries()); i++) 
        {
          const higgsana::TPhoton *pho = (higgsana::TPhoton*)((*photonArr)[i]);      
          
          if( printDebug ) cout << "\npho:: " << pho->et << " " << pho->eta << " " << pho->phi << " : " ;

          if (!(pho->et > 15)) continue;
          if (!(fabs(pho->scEta) < 2.5)) continue;


          if (PassSelectionStages[0] && PassSelectionStages[1] && PassSelectionStages[2]) {
            phoList << info->evtNum << " " << pho->et << " " << pho->eta << " " << PassPhotonHZGammaID(pho, (DataEra == kDataEra_2011_MC)) << " " << PassPhotonHZGammaIso( pho, pfcandidateArr, rhoEleIso,EleEAEra) << "\n";
          }

          if( printDebug ) cout << " pass presel : " ;

          if (!PassPhotonHZGammaID(pho, (DataEra == kDataEra_2011_MC), printDebug)) {
            if( printDebug ) cout << "fail ID \n";
            continue;
          }
          if( printDebug ) cout << " pass ID : " ;
          if (!PassPhotonHZGammaIso( pho, pfcandidateArr, rhoEleIso,EleEAEra, printDebug)) {
            if( printDebug ) cout << "fail iso \n";
            continue;
          } 
          if( printDebug ) cout << " pass iso : " ;
          
          NPhoPassSelection++;

          //Look for lepton overlaps
          Bool_t IsSeparatedFromLeptons = kTRUE;
          for (uint l=0; l<leptonType.size(); l++) {          
            TLorentzVector v;
            if (leptonType[l] == 11) {
              const higgsana::TElectron *tmpele = (higgsana::TElectron*)((*electronArr)[leptonIndex[l]]);  
              v.SetPtEtaPhiM( tmpele->pt, tmpele->eta, tmpele->phi, ELECTRONMASS);
            } else {
              const higgsana::TMuon *tmpmu = (higgsana::TMuon*)((*muonArr)[leptonIndex[l]]);      
              v.SetPtEtaPhiM( tmpmu->pt, tmpmu->eta, tmpmu->phi, MUONMASS);
            }
           
            if (higgsana::deltaR(pho->eta, pho->phi, v.Eta(), v.Phi()) <= 0.7) {
              IsSeparatedFromLeptons = kFALSE;
              break;
            }
          }


          if (IsSeparatedFromLeptons) {
            if( printDebug ) cout << " pass DR cut : \n" ;
            NPhoPassDRCut++;
            photonIndex.push_back(i);
          }
        }

        if( printDebug ) cout << "\n NPhoPassSelection " << NPhoPassSelection << " : NPhoPassDRCut " << NPhoPassDRCut << endl;

        //***********************************************************
        // Selection Stage 5: Photon Selection 
        //***********************************************************
        if (NPhoPassSelection >= 1) {
          PassSelectionStages[5] = kTRUE;
        } else {
          PassSelectionStages[5] = kFALSE;
        }

        //***********************************************************
        // Selection Stage 6:  Lepton - Photon dR cut
        //***********************************************************
        if (NPhoPassDRCut >= 1) {
          PassSelectionStages[6] = kTRUE;
        } else {
          PassSelectionStages[6] = kFALSE;
        }
 

        //***********************************************************
        // Selection Stage 7:  Find Best Z Candidate
        //***********************************************************
        Int_t ZLeptonType = 0;
        Int_t ZLepton1Index = -1;
        Int_t ZLepton2Index = -1;
        TLorentzVector ZLepton1FourVector;
        TLorentzVector ZLepton2FourVector;
        Double_t BestZCandMass = 0;

        for (uint l=0; l<leptonType.size(); l++) {          
          TLorentzVector v1;
          Int_t Lep1Charge = 0;
          if (leptonType[l] == 11) {
            const higgsana::TElectron *tmpele1 = (higgsana::TElectron*)((*electronArr)[leptonIndex[l]]);  
            v1.SetPtEtaPhiM( tmpele1->pt, tmpele1->eta, tmpele1->phi, ELECTRONMASS);
            Lep1Charge = tmpele1->q;
          } else {
            const higgsana::TMuon *tmpmu1 = (higgsana::TMuon*)((*muonArr)[leptonIndex[l]]);      
            v1.SetPtEtaPhiM( tmpmu1->pt, tmpmu1->eta, tmpmu1->phi, MUONMASS);
            Lep1Charge = tmpmu1->q;
          }

          for (uint k=l+1; k<leptonType.size(); k++) {          
            Int_t Lep2Charge = 0;
            TLorentzVector v2;
            if (leptonType[k] == 11) {
              const higgsana::TElectron *tmpele2= (higgsana::TElectron*)((*electronArr)[leptonIndex[k]]);  
              v2.SetPtEtaPhiM( tmpele2->pt, tmpele2->eta, tmpele2->phi, ELECTRONMASS);
              Lep2Charge = tmpele2->q;
            } else {
              const higgsana::TMuon *tmpmu2 = (higgsana::TMuon*)((*muonArr)[leptonIndex[k]]);      
              v2.SetPtEtaPhiM( tmpmu2->pt, tmpmu2->eta, tmpmu2->phi, MUONMASS);
              Lep2Charge = tmpmu2->q;
            }

            if ( PassSelectionStages[0] && PassSelectionStages[1] && PassSelectionStages[2] && PassSelectionStages[3] && PassSelectionStages[4] && PassSelectionStages[5] && PassSelectionStages[6]) {
              cout << k << " " << l << " : " << leptonType[k] << " " << leptonType[l] << " : " << Lep1Charge << " " << Lep2Charge << endl;
            }

            //same flavor, opposite charge
            if( !(leptonType[k] == leptonType[l])) continue;
            if( !(Lep1Charge != Lep2Charge)) continue;

            //best Z candidate
            if ( fabs((v1+v2).M() - 91.2) < fabs (BestZCandMass - 91.2) ) {
              ZLeptonType = leptonType[l];
              ZLepton1Index = l;
              ZLepton2Index = k;
              BestZCandMass = (v1+v2).M();
              Double_t lepMass = (leptonType[l] == 1) ? ELECTRONMASS : MUONMASS;
              ZLepton1FourVector.SetPtEtaPhiM( v1.Pt(), v1.Eta(), v1.Phi(), lepMass);
              ZLepton2FourVector.SetPtEtaPhiM( v2.Pt(), v2.Eta(), v2.Phi(), lepMass);
            }

          }
        }

        TLorentzVector  ZFourVector = ZLepton1FourVector+ZLepton2FourVector;
        if ((ZLepton1Index >= 0 && ZLepton2Index >= 0)) {
          PassSelectionStages[7] = kTRUE;
        } else {
          PassSelectionStages[8] = kFALSE;
        }

        //***********************************************************
        // Selection Stage 8:  Z mass cut
        //***********************************************************
        if (ZFourVector.M() > 50) {
          PassSelectionStages[8] = kTRUE;
        } else {
          PassSelectionStages[8] = kFALSE;        
        }

        //***********************************************************
        // Selection Stage 9  Zg Candidate
        //***********************************************************
        Int_t PhotonIndex = -1;
        TLorentzVector PhotonFourVector;
        Double_t BestCandidatePhotonPt = 0;

        for (uint p=0; p<photonIndex.size(); p++) {          
          const higgsana::TPhoton *tmppho = (higgsana::TPhoton*)((*photonArr)[photonIndex[p]]);
          
          if (tmppho->et > BestCandidatePhotonPt) {
            PhotonIndex= p;
            PhotonFourVector.SetPtEtaPhiM( tmppho->et, tmppho->eta, tmppho->phi, 0.0);
            BestCandidatePhotonPt = tmppho->et;
          }
        }

        TLorentzVector HiggsFourVector = ZFourVector+PhotonFourVector;
        if (PhotonIndex >= 0) {
          PassSelectionStages[9] = kTRUE;
        } else {
          PassSelectionStages[9] = kFALSE;
        }

    

        //***********************************************************
        // Selection Stage 10:  Z mass cut
        //***********************************************************
        if (HiggsFourVector.M() > 110) {
          PassSelectionStages[10] = kTRUE;
        } else {
          PassSelectionStages[10] = kFALSE;        
        }

        //***********************************************************
        // Selection Stage 11:  Z mass cut
        //***********************************************************
        if (HiggsFourVector.M() < 180) {
          PassSelectionStages[11] = kTRUE;
        } else {
          PassSelectionStages[11] = kFALSE;        
        }



        //********************************************************
        // Dump Stuff
        //********************************************************



        //*********************************************************
        //Tabulate Yields
        //*********************************************************
        NEventsTotal++;
        for (UInt_t n=0; n < NSelectionStages; ++n) {

          Bool_t passUpToCurrentStage = kTRUE;

          for (UInt_t q=0; q <= n; ++q) {
            passUpToCurrentStage = passUpToCurrentStage && PassSelectionStages[q];
          }

          if (passUpToCurrentStage) {
            eventList << "1";
            NEventsPassSelectionStages[n]++;            
          } else {
            eventList << "0";
          }
          eventList << " ";
        }

        eventList << endl;
      } //end loop over data     
    }
  }

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;


  //********************************************
  //Show Yields
  //********************************************
  cout << "Total Events: " << NEventsTotal << endl;
  for (UInt_t n=0; n < NSelectionStages; ++n) {

    cout << "Selection Stage " << n << " : " 
         << NEventsPassSelectionStages[n]  
         << endl;    
  }

   
  gBenchmark->Show("WWTemplate");       
} 



