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
#include "HiggsAna/Ntupler/interface/HiggsAnaDefs.hh"
#include "HiggsAna/Ntupler/interface/TEventInfo.hh"
#include "HiggsAna/Ntupler/interface/TElectron.hh"
#include "HiggsAna/Ntupler/interface/TPhoton.hh"
#include "HiggsAna/Ntupler/interface/TMuon.hh"
#include "HiggsAna/Ntupler/interface/TJet.hh"
#include "HiggsAna/Ntupler/interface/TGenInfo.hh"
#include "HiggsAna/Ntupler/interface/TGenParticle.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// Helper functions for Electron ID selection
#include "HiggsAna/Utils/LeptonIDCuts.hh"
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodSwitches.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodMeasurements.h"

// structure for output ntuple
#include "EffData.hh" 
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



//=== MAIN MACRO ================================================================================================= 

void selectEleDenominatorEffTP(const TString conf,              // input file
                      const TString outputDir,         // output directory
		      const Bool_t  matchGen = kFALSE  // match to generator muons
) {
  gBenchmark->Start("selectEleLHEffTP");

  TH1F *DRStatus1ToStatus3_Status1Matched = new TH1F("DRStatus1ToStatus3_Status1Matched", ";DR;Number of Events", 100, 0, 0.2);
  TH1F *DRStatus1ToStatus3_Status3Matched = new TH1F("DRStatus1ToStatus3_Status3Matched", ";DR;Number of Events", 100, 0, 0.2);
  TH1F *Mass_Status1Matched = new TH1F("Mass_Status1Matched", ";Mass;Number of Events", 35, 50, 120);
  TH1F *Mass_Status3Matched = new TH1F("Mass_Status3Matched", ";Mass;Number of Events", 35, 50, 120);

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  //*****************************************************************************************
  //Setup LH
  //*****************************************************************************************
  TFile *fileLH = TFile::Open("/home/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronLikelihoodPdfs_MC.root");
  TDirectory *EB0lt15dir = fileLH->GetDirectory("/");
  TDirectory *EB1lt15dir = fileLH->GetDirectory("/");
  TDirectory *EElt15dir = fileLH->GetDirectory("/");
  TDirectory *EB0gt15dir = fileLH->GetDirectory("/");
  TDirectory *EB1gt15dir = fileLH->GetDirectory("/"); 
  TDirectory *EEgt15dir = fileLH->GetDirectory("/");

  LikelihoodSwitches defaultSwitches;
  defaultSwitches.m_useEoverP = false;
  defaultSwitches.m_useDeltaEta = true;
  defaultSwitches.m_useDeltaPhi = true;
  defaultSwitches.m_useHoverE = false;        
  defaultSwitches.m_useSigmaEtaEta = true;
  defaultSwitches.m_useSigmaPhiPhi = true;
  defaultSwitches.m_useFBrem = true;
  defaultSwitches.m_useOneOverEMinusOneOverP = true;
 
  ElectronLikelihood *LH = new ElectronLikelihood(&(*EB0lt15dir),&(*EB1lt15dir), &(*EElt15dir), 
                                                  &(*EB0gt15dir), &(*EB1gt15dir), &(*EEgt15dir),
                                                  defaultSwitches,
                                                  std::string("class"),std::string("class"),true,true);



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
  
  Double_t NTotalProbes = 0;
  Double_t NGenMatchedStatus1Probes = 0;
  Double_t NGenMatchedStatus3Probes = 0;
  Double_t NNotMatched = 0;
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
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  mithep::TGenInfo *gen     = new mithep::TGenInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *genparticleArr = new TClonesArray("mithep::TGenParticle");
  
  // loop over files  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]); 
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
    eventTree->SetBranchAddress("GenParticle", &genparticleArr); TBranch *genparticleBr = eventTree->GetBranch("GenParticle");      
    cout << "NEvents = " << eventTree->GetEntries() << endl;

    TBranch *genBr = 0;
    if(matchGen) {
      eventTree->SetBranchAddress("Gen", &gen);
      genBr = eventTree->GetBranch("Gen");    
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
      if(matchGen) {
        genparticleBr->GetEntry(ientry);
      }

      // check for certified runs
      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  

      // trigger requirement               
      ULong_t  trigger = 0;
      if(typev[ifile]==eDiEl) {
        trigger = kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 |
		  kHLT_Ele32_CaloIdL_CaloIsoVL_SC17 |
	          kHLT_Ele17_CaloIdL_CaloIsoVL;
        
      } else if(typev[ifile]==eEl) {
        if(info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30)  continue;
	if(info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17)                         continue;
	if(info->triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL)                              continue;
        
	trigger = kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT;
      }
      if(typev[ifile]!=eMC && !(info->triggerBits & trigger)) continue;     
      
      vector<Double_t> GenPstatus;
      vector<Double_t> GenPPt;
      vector<Double_t> GenPEta;
      vector<Double_t> GenPPhi;
      if (typev[ifile] == eMC) {
        for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
          const mithep::TGenParticle *gen = (mithep::TGenParticle*)((*genparticleArr)[k]);
          if (abs(gen->pdgid) == 11 ) {
            GenPstatus.push_back(gen->status);
            GenPPt.push_back(gen->pt); 
            GenPEta.push_back(gen->eta); 
            GenPPhi.push_back(gen->phi); 
          }
        }
      }


      // good vertex requirement
      if(!(info->hasGoodPV)) continue;
      
      if(matchGen) {
        genBr->GetEntry(ientry);
      }

      electronArr->Clear();
      electronBr->GetEntry(ientry);
      for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
        const mithep::TElectron *tag = (mithep::TElectron*)((*electronArr)[i]);
	
	if(matchGen) {

// 	  Bool_t match1 = (toolbox::deltaR(tag->eta, tag->phi, gen->eta_1, gen->phi_1) < 0.5);
// 	  Bool_t match2 = (toolbox::deltaR(tag->eta, tag->phi, gen->eta_2, gen->phi_2) < 0.5);
// 	  if(!match1 && !match2)
// 	    continue;

          Bool_t match = kFALSE;
          
//           for(Int_t q=0; q< GenPPt.size() ; ++q) {
//             if (GenPstatus[q] == 1 && toolbox::deltaR(tag->eta, tag->phi, GenPEta[q], GenPPhi[q]) < 0.5) match = kTRUE;
//           }
//           if (!match) continue;
	}
	
	if(tag->pt          < 20)  continue;
	if(fabs(tag->scEta) > 2.5) continue;
 	if(!passCutBasedTightEleID(tag))        continue;
	
	if(typev[ifile]==eDiEl &&
	   !((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (tag->hltMatchBits & kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT)) &&
	   !((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (tag->hltMatchBits & kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT)) &&
	   !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (tag->hltMatchBits & kHLTObject_Ele32_CaloIdL_CaloIsoVL)) &&
	   !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (tag->hltMatchBits & kHLTObject_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT)) &&
	   !((info->triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL) && (tag->hltMatchBits & kHLTObject_Ele17_CaloIdL_CaloIsoVL)) ) 
	  continue;
      
        if(typev[ifile]==eEl &&
	   !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLTObject_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT)) &&
	   !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLTObject_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT)) &&
	   !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLTObject_Ele52_CaloIdVT_TrkIdT)) ) 
	  continue;
        
	const Double_t m = 0.000511;
	TLorentzVector vtag;
	vtag.SetPtEtaPhiM(tag->pt, tag->eta, tag->phi, m);
        
	for(Int_t j=0; j<electronArr->GetEntriesFast(); j++) {
	  if(i==j) continue;
	  
	  const mithep::TElectron *probe = (mithep::TElectron*)((*electronArr)[j]);
	  if(probe->q == tag->q) continue;
	  
	  if(typev[ifile]==eDiEl &&
	     !((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (probe->hltMatchBits & kHLTObject_SC8)) &&
	     !((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (probe->hltMatchBits & kHLTObject_Ele8)) &&
	     !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (probe->hltMatchBits & kHLTObject_SC17)))
	    continue;
	  

	  if(fabs(probe->scEta) > 2.5) continue;	  
	  TLorentzVector vprobe;
	  vprobe.SetPtEtaPhiM(probe->pt, probe->eta, probe->phi, m);
	  
	  TLorentzVector vdielectron = vtag + vprobe;
	  if((vdielectron.M()<massLo) || (vdielectron.M()>massHi)) continue;	  	  


          Double_t status1MinDR = 9999;
          Double_t status1Eta;
          Double_t status1Phi;
          Double_t status3MinDR = 9999;
          Double_t status3Eta;
          Double_t status3Phi;
          if (typev[ifile] == eMC) {
            for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
              const mithep::TGenParticle *gen = (mithep::TGenParticle*)((*genparticleArr)[k]);
              if (abs(gen->pdgid) == 11 ) {
                if (gen->status == 1 &&  toolbox::deltaR(probe->eta, probe->phi, gen->eta, gen->phi) < status1MinDR) {
                  status1MinDR = toolbox::deltaR(probe->eta, probe->phi, gen->eta, gen->phi);
                  status1Eta = gen->eta;
                  status1Phi = gen->phi;
                }
                if (gen->status == 3 &&  toolbox::deltaR(probe->eta, probe->phi, gen->eta, gen->phi) < status3MinDR) {
                  status3MinDR = toolbox::deltaR(probe->eta, probe->phi, gen->eta, gen->phi);
                  status3Eta = gen->eta;
                  status3Phi = gen->phi;
                }
              }
            }
          }

          if (status1MinDR > 0.5 || status3MinDR > 0.5) {
//            cout << "UNMATCHED ELE: " << probe->pt << " " << probe->eta << " " << probe->phi << " : " << probe->eta << " " << probe->scPhi << endl;          
            for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
              const mithep::TGenParticle *gen = (mithep::TGenParticle*)((*genparticleArr)[k]);
              if (abs(gen->pdgid) == 11) {
//                cout << "GEN: " << gen->pdgid << " : " << gen->status << " " << gen->pt << " " << gen->eta << " " << gen->phi << " : " << toolbox::deltaR(probe->eta, probe->scPhi, gen->eta, gen->phi) << endl;
              }
              if (abs(gen->pdgid) == 11 ) {
                if (gen->status == 1 &&  toolbox::deltaR(probe->eta, probe->phi, gen->eta, gen->phi) < status1MinDR) {
                  status1MinDR = toolbox::deltaR(probe->eta, probe->phi, gen->eta, gen->phi);
                  status1Eta = gen->eta;
                  status1Phi = gen->phi;
                }
                if (gen->status == 3 &&  toolbox::deltaR(probe->eta, probe->scPhi, gen->eta, gen->phi) < status3MinDR) {
                  status3MinDR = toolbox::deltaR(probe->eta, probe->scPhi, gen->eta, gen->phi);
                  status3Eta = gen->eta;
                  status3Phi = gen->phi;
                }
              }
            }
          }


          NTotalProbes++;
	  if(matchGen) {
// 	    Bool_t match1 = (toolbox::deltaR(probe->eta, probe->phi, gen->eta_1, gen->phi_1) < 0.5);
// 	    Bool_t match2 = (toolbox::deltaR(probe->eta, probe->phi, gen->eta_2, gen->phi_2) < 0.5);
// 	    if(!match1 && !match2)
// 	      continue;

            Bool_t matchStatus1 = kFALSE;          
            Bool_t matchStatus3 = kFALSE;          
            for(Int_t q=0; q< GenPPt.size() ; ++q) {
              if (GenPstatus[q] == 1 && toolbox::deltaR(probe->scEta, probe->scPhi, GenPEta[q], GenPPhi[q]) < 0.1) matchStatus1 = kTRUE;
              if (GenPstatus[q] == 3 && toolbox::deltaR(probe->scEta, probe->scPhi, GenPEta[q], GenPPhi[q]) < 0.1) matchStatus3 = kTRUE;
            }
            if (matchStatus1) {
              NGenMatchedStatus1Probes++;
              if (probe->pt > 10 && probe->pt < 20 && fabs(probe->scEta) < 1.479) {
                DRStatus1ToStatus3_Status1Matched->Fill( toolbox::deltaR(status1Eta, status1Phi,status3Eta, status3Phi));
                Mass_Status1Matched->Fill(vdielectron.M());
              }
            }
            if (matchStatus3) {
              NGenMatchedStatus3Probes++;
              if (probe->pt > 10 && probe->pt < 20 && fabs(probe->scEta) < 1.479) {
                DRStatus1ToStatus3_Status3Matched->Fill( toolbox::deltaR(status1Eta, status1Phi,status3Eta, status3Phi));
                Mass_Status3Matched->Fill(vdielectron.M());
              }
            }
            if (!matchStatus1 && !matchStatus3) {
              NNotMatched++;
            }

            if (!matchStatus1) continue;
	  }

	  nProbes++;
	  
	  Bool_t pass = kTRUE; //passElectronDenominatorV4(probe);

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
          data.rho    = info->PileupEnergyDensity;
	  outTree->Fill();	  
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
      
  cout << "Total Probes:" << NTotalProbes << endl;
  cout << "Status1 Matched Probes:" << NGenMatchedStatus1Probes << endl;
  cout << "Status3 Matched Probes:" << NGenMatchedStatus3Probes << endl;
  cout << "NotMatched: " << NNotMatched << endl;

  NormalizeHist(  DRStatus1ToStatus3_Status1Matched);
  NormalizeHist(  DRStatus1ToStatus3_Status3Matched);
  NormalizeHist(  Mass_Status1Matched);
  NormalizeHist(  Mass_Status3Matched);

  TFile *outputFile = new TFile("GenMatching.root", "UPDATE");
  outputFile->WriteTObject(DRStatus1ToStatus3_Status1Matched, DRStatus1ToStatus3_Status1Matched->GetName(), "WriteDelete");
  outputFile->WriteTObject(DRStatus1ToStatus3_Status3Matched, DRStatus1ToStatus3_Status3Matched->GetName(), "WriteDelete");
  outputFile->WriteTObject(Mass_Status1Matched, Mass_Status1Matched->GetName(), "WriteDelete");
  outputFile->WriteTObject(Mass_Status3Matched, Mass_Status3Matched->GetName(), "WriteDelete");
  outputFile->Close();

  gBenchmark->Show("selectEleLHEffTP"); 
}
