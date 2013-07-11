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

#include "Common/MyTools.hh"        // miscellaneous helper functions
 
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
#include "EffData.hh" 
#endif

//=== Functions  ================================================================================================= 

Bool_t MatchedToGen( const higgsana::TElectron *ele , TClonesArray *genparticleArr) {
  Bool_t matched = kFALSE;
  for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
    const higgsana::TGenParticle *gen = (higgsana::TGenParticle*)((*genparticleArr)[k]);

    //match to status 1 or status 3
    if (abs(gen->pdgid) == 11 && (gen->status == 1 || gen->status == 3)) {
      if ( higgsana::deltaR( ele->eta, ele->phi, gen->eta , gen->phi) < 0.3) {
        matched = kTRUE;
        break;
      }
    }
  }
  return matched;
}

Bool_t MatchedToGen( const higgsana::TMuon *mu , TClonesArray *genparticleArr) {
  Bool_t matched = kFALSE;
  for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
    const higgsana::TGenParticle *gen = (higgsana::TGenParticle*)((*genparticleArr)[k]);

    //match to status 1 or status 3
    if (abs(gen->pdgid) == 13 && (gen->status == 1 || gen->status == 3)) {
      if ( higgsana::deltaR( ele->eta, ele->phi, gen->eta , gen->phi) < 0.3) {
        matched = kTRUE;
        break;
      }
    }
  }
  return matched;
}

//=== MAIN MACRO ================================================================================================= 

void selectEleHZZICHEP2012WP_Zto4l(const TString conf,              // input file
                                   const TString outputDir,         // output directory
                                   const Bool_t  matchGen = kFALSE, // match to generator muons
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
  eleIDMVA->initialize( "BDT", EGammaMvaEleEstimator::kNonTrig, kTRUE, weightFiles);
  



  // mass region
  Double_t massLo;
  Double_t massHi;

  Double_t lumi;              // luminosity (pb^-1)
  
  vector<TString>  fnamev;    // sample files
  vector<Int_t>    typev;     // dataset type 
  vector<Double_t> xsecv;     // per file cross section
  vector<TString>  jsonv;     // per file JSON file

  //
  // parse .conf file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') { 
      state++; 
      continue; 
    }
    
    if(state==0) {  // general settings
      stringstream ss1(line); ss1 >> lumi;
      getline(ifs,line);
      stringstream ss2(line); ss2 >> massLo >> massHi; 
      
    } else if(state==1) {  // define data sample
      string fname;
      Int_t type;
      Double_t xsec;
      string json;
      stringstream ss(line);
      ss >> fname >> type >> xsec >> json;
      fnamev.push_back(fname);
      typev.push_back(type);
      xsecv.push_back(xsec);
      jsonv.push_back(json);        
    }
  }
  ifs.close();
  
  enum { eMC, eMuEl, eDiMu, eMu, eDiEl, eEl };  // dataset type
  
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  Double_t nProbes = 0;
  
  //
  // Set up output ntuple
  //
  gSystem->mkdir(outputDir,kTRUE);
  TFile *outFile = new TFile(outputDir+TString("/probes.root"),"RECREATE"); 
  TTree *outTree = new TTree("Events","Events");
  EffData data;
  outTree->Branch("Events",&data.mass,"mass/F:pt:eta:phi:weight:q/I:npv/i:npu:pass:runNum:lumiSec:evtNum:rho/F");

  TFile *infile=0;
  TTree *eventTree=0;
  
  // Data structures to store info from TTrees
  higgsana::TEventInfo *info  = new higgsana::TEventInfo();
  TClonesArray *genparticleArr = new TClonesArray("higgsana::TGenParticle");
  TClonesArray *electronArr = new TClonesArray("higgsana::TElectron");
  TClonesArray *pfcandidateArr = new TClonesArray("higgsana::TPFCandidate");
  TClonesArray *muonArr = new TClonesArray("higgsana::TMuon");

  // loop over files  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = TFile::Open(fnamev[ifile],"read");
    assert(infile);

    Bool_t hasJSON = kFALSE;
    mithep::RunLumiRangeMap rlrm;
    if(jsonv[ifile].CompareTo("NONE")!=0) { 
      hasJSON = kTRUE;
      rlrm.AddJSONFile(jsonv[ifile].Data()); 
    }
    
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
    eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("Muon", &muonArr);         TBranch *muonBr = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr);  TBranch *pfcandidateBr = eventTree->GetBranch("PFCandidate");
    cout << "NEvents = " << eventTree->GetEntries() << endl;

    //TBranch *genBr = 0;
    TBranch *genparticleBr;
    if(matchGen) {
      //eventTree->SetBranchAddress("Gen", &gen);
      //genBr = eventTree->GetBranch("Gen");
      eventTree->SetBranchAddress("GenParticle", &genparticleArr);
      genparticleBr = eventTree->GetBranch("GenParticle");
    }
    
    // Determine maximum number of events to consider
    // *** CASES ***
    // <> lumi < 0 => use all events in the sample
    // <> xsec = 0 => for data (use all events)
    const Double_t xsec = xsecv[ifile];
    Double_t weight = 1;
    if(lumi>0) { 
      if(xsec>0) { weight = lumi*xsec/(Double_t)eventTree->GetEntries(); }      
    }

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

      // trigger requirement 
      // already done at skim step
      
      // good vertex requirement
      if(!(info->hasGoodPV)) continue;

      genparticleArr->Clear(); 
      electronArr->Clear();
      muonArr->Clear(); 
      pfcandidateArr->Clear(); 
      
      electronBr->GetEntry(ientry);
      muonBr->GetEntry(ientry);
      pfcandidateBr->GetEntry(ientry);
      if(matchGen) {
        //genBr->GetEntry(ientry);
        genparticleBr->GetEntry(ientry);
      }  


      //********************************************************
      //Make Lepton Vector
      //********************************************************
      vector<Int_t> leptonType;
      vector<UInt_t> leptonIndex;
      for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
        const higgsana::TElectron *ele = (higgsana::TElectron*)((*electronArr)[i]);
        if (ele->pt > 7 && fabs(ele->eta) < 2.5) {
          leptonType.push_back(11 * ele->q);
          leptonIndex.push_back(i);
        }
      }
      for(uint i=0; i<uint(muonArr->GetEntries()); i++) 
      {
        const higgsana::TMuon *mu = (higgsana::TMuon*)((*muonArr)[i]);      
        if (mu->pt > 5 && fabs(mu->eta) < 2.4) {
          leptonType.push_back(13 * mu->q);
          leptonIndex.push_back(i);
        }
      }
      

      //********************************************************
      //Loop over all 4-lepton combinations
      //********************************************************
      vector<const higgsana::TPFCandidate*> photonsToVeto; //don't do FSR for now
      
      for(Int_t lepIndex1=0; lepIndex1<leptonType.size(); ++lepIndex1) {

        TLorentzVector Lep1Vector;
          //Tag1 Lepton Requirements
        Bool_t tag1Pass = kFALSE;
        if (abs(leptonType[lepIndex1]) == 11) {
          const higgsana::TElectron *ele1 = (higgsana::TElectron*)((*electronArr)[leptonIndex[lepIndex1]]);
          if (matchGen && !MatchedToGen(ele1, genparticleArr ) ) continue; 
          tag1Pass = PassEleHZZ4lPreselection(ele1) && 
            PassEleHZZ4lICHEP2012ID(ele1, eleIDMVA, printDebug) && 
            PassEleHZZ4lICHEP2012Iso(ele1,pfcandidateArr,rhoEleIso,EleEAEra,photonsToVeto);
          Lep1Vector.SetPtEtaPhiM( ele1->pt, ele1->eta, ele1->phi, ELECTRONMASS);
        } else if (abs(leptonType[lepIndex1]) == 13) {
          const higgsana::TMuon *mu1 = (higgsana::TMuon*)((*muonArr)[leptonIndex[lepIndex1]]);
          if (matchGen && !MatchedToGen(mu1, genparticleArr ) ) continue; 
          tag1Pass = PassMuonHZZ4lPreselection(mu1) && 
            PassMuonHZZ4lICHEP2012ID(mu1,i,pfcandidateArr) &&
            PassMuonHZZ4lICHEP2012Iso(mu1, i, pfcandidateArr,rhoMuonIso,MuonEAEra,photonsToVeto);
          Lep1Vector.SetPtEtaPhiM( mu1->pt, mu1->eta, mu1->phi, MUONMASS);
        }
        if (!tag1Pass) continue;

        for(Int_t lepIndex2=0; lepIndex2<leptonType.size(); ++lepIndex2) {

          TLorentzVector Lep2Vector;
          
          //don't use duplicate leptons
          if (lepIndex2 == lepIndex1) continue;

          //Tag2 Lepton Requirements
          Bool_t tag2Pass = kFALSE;          
          if (abs(leptonType[lepIndex2]) == 11) {
            const higgsana::TElectron *ele2 = (higgsana::TElectron*)((*electronArr)[leptonIndex[lepIndex2]]);
            if (matchGen && !MatchedToGen(ele2, genparticleArr ) ) continue; 
            tag2Pass = PassEleHZZ4lPreselection(ele2) && 
              PassEleHZZ4lICHEP2012ID(ele2, eleIDMVA, printDebug) && 
              PassEleHZZ4lICHEP2012Iso(ele2,pfcandidateArr,rhoEleIso,EleEAEra,photonsToVeto);
            Lep2Vector.SetPtEtaPhiM( ele2->pt, ele2->eta, ele2->phi, ELECTRONMASS);
          } else if (abs(leptonType[lepIndex3]) == 13) {
            const higgsana::TMuon *mu2 = (higgsana::TMuon*)((*muonArr)[leptonIndex[lepIndex2]]);
            if (matchGen && !MatchedToGen(mu2, genparticleArr ) ) continue; 
            tag2Pass = PassMuonHZZ4lPreselection(mu2) && 
              PassMuonHZZ4lICHEP2012ID(mu2,i,pfcandidateArr) &&
              PassMuonHZZ4lICHEP2012Iso(mu2, i, pfcandidateArr,rhoMuonIso,MuonEAEra,photonsToVeto);
            Lep2Vector.SetPtEtaPhiM( mu2->pt, mu2->eta, mu2->phi, MUONMASS);
          }
          if (!tag2Pass) continue;

          for(Int_t lepIndex3=0; lepIndex3<leptonType.size(); ++lepIndex3) {

            TLorentzVector Lep3Vector;

            //don't use duplicate leptons
            if (lepIndex3 == lepIndex1 || lepIndex3 == lepIndex2) continue;

            //Tag2 Lepton Requirements
            Bool_t tag3Pass = kFALSE;
            if (abs(leptonType[lepIndex3]) == 11) {
              const higgsana::TElectron *ele3 = (higgsana::TElectron*)((*electronArr)[leptonIndex[lepIndex3]]);
              if (matchGen && !MatchedToGen(ele3, genparticleArr ) ) continue; 
              tag3Pass = PassEleHZZ4lPreselection(ele3) && 
                PassEleHZZ4lICHEP2012ID(ele3, eleIDMVA, printDebug) && 
                PassEleHZZ4lICHEP2012Iso(ele3,pfcandidateArr,rhoEleIso,EleEAEra,photonsToVeto);
              Lep3Vector.SetPtEtaPhiM( ele3->pt, ele3->eta, ele3->phi, ELECTRONMASS);
            } else if (abs(leptonType[lepIndex3]) == 13) {
              const higgsana::TMuon *mu3 = (higgsana::TMuon*)((*muonArr)[leptonIndex[lepIndex3]]);
              if (matchGen && !MatchedToGen(mu3, genparticleArr ) ) continue; 
              tag3Pass = PassMuonHZZ4lPreselection(mu3) && 
                PassMuonHZZ4lICHEP2012ID(mu3,i,pfcandidateArr) &&
                PassMuonHZZ4lICHEP2012Iso(mu3, i, pfcandidateArr,rhoMuonIso,MuonEAEra,photonsToVeto);
              Lep3Vector.SetPtEtaPhiM( mu3->pt, mu3->eta, mu3->phi, MUONMASS);
           }
            if (!tag3Pass) continue;

            for(Int_t lepIndex4=0; lepIndex4<leptonType.size(); ++lepIndex4) {

              TLorentzVector Lep4Vector;
          
              //don't use duplicate leptons
              if (lepIndex4 == lepIndex1 || lepIndex4 == lepIndex2 || lepIndex4 == lepIndex3) continue;

              //probe only electrons
              if (!(abs(leptonType[lepIndex4]) == 11)) continue;

              const higgsana::TElectron *ele4 = (higgsana::TElectron*)((*electronArr)[leptonIndex[lepIndex4]]);
              if (matchGen && !MatchedToGen(ele4, genparticleArr ) ) continue; 
              Lep4Vector.SetPtEtaPhiM( ele4->pt, ele4->eta, ele4->phi, ELECTRONMASS);

              //***********************************************************
              // Some More Selection Cuts
              //***********************************************************
              
              //charge-flavor requirement
              if (!(leptonType[lepIndex1] + leptonType[lepIndex2] + 
                    leptonType[lepIndex3] + leptonType[lepIndex4] == 0)
                ) continue;
              
              //pt1,pt2 cuts
              Int_t NLeptonsAbove20 = 0;
              Int_t NLeptonsAbove10 = 0;
              if (Lep1Vector.Pt() > 20) NLeptonsAbove20++;
              if (Lep1Vector.Pt() > 10) NLeptonsAbove10++;
              if (Lep2Vector.Pt() > 20) NLeptonsAbove20++;
              if (Lep2Vector.Pt() > 10) NLeptonsAbove10++;
              if (Lep3Vector.Pt() > 20) NLeptonsAbove20++;
              if (Lep3Vector.Pt() > 10) NLeptonsAbove10++;
              if (Lep4Vector.Pt() > 20) NLeptonsAbove20++;
              if (Lep5Vector.Pt() > 10) NLeptonsAbove10++;
              if (!(NLeptonsAbove20>=1 && NLeptonsAbove10 >= 2)) continue;

              //resonance killing
              if ( (Lep1Vector+Lep2Vector).M() < 4 ) continue;
              if ( (Lep1Vector+Lep3Vector).M() < 4 ) continue;
              if ( (Lep1Vector+Lep4Vector).M() < 4 ) continue;
              if ( (Lep2Vector+Lep3Vector).M() < 4 ) continue;
              if ( (Lep2Vector+Lep4Vector).M() < 4 ) continue;
              if ( (Lep3Vector+Lep4Vector).M() < 4 ) continue;

              //mass window
              TLorentzVector FourLeptonVector = Lep1Vector + Lep2Vector + Lep3Vector + Lep4Vector;
              if((FourLeptonVector.M()<40) || (FourLeptonVector.M()>200)) continue;

              nProbes++;
              
              const higgsana::TElectron *probe = ele4;
              Bool_t passID = (PassEleHZZ4lPreselection(probe) && PassEleHZZ4lICHEP2012ID(probe, eleIDMVA));
              Bool_t passIsolation = PassEleHZZ4lICHEP2012Iso(probe,pfcandidateArr,rhoEleIso,EleEAEra,photonsToVeto);
              
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
      }

    }
    delete infile;
    infile=0, eventTree=0;    
  }
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
  if(lumi>0) {
    cout << " L_int = " << lumi << "/pb" << endl;
    cout << endl;
  }
  cout << " Number of probes selected: " << nProbes << endl;
  
  outFile->Write();
  outFile->Close();
  delete outFile;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("selectEleLHEffTP"); 
}
