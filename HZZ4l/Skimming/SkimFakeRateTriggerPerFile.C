// *************************************
// From EOS
//
// cmsLs /store/user/sixie/hist/AllNtuple/cern/filefi/028/ >! fileListA
// cmsPfn `cat fileListA | awk '{print $5}'` | sed 's/?svcClass=default//' | sed 's/root:\/\/eoscms\/\/eos\/cms\/store\/user\/sixie\/hist\/AllNtuple\/cern\/filefi\/028\///' > ! fileListA2
//
// cmsLs /store/group/phys_higgs/cmshzz4l/sixie/hist/AllNtuple/cern/filefi/028/ >! fileListB
// cmsPfn `cat fileListB | awk '{print $5}'` | sed 's/?svcClass=default//' | sed 's/root:\/\/eoscms\/\/eos\/cms\/store\/group\/phys_higgs\/cmshzz4l\/sixie\/hist\/AllNtuple\/cern\/filefi\/028\///' > ! fileListB2
//
// foreach f(`cat fileListB2 | grep HZZ4lNtuple | grep r12b-dmu-pr-v1`)
// root -l -b -q HiggsAna/HZZ4l/Skimming/SkimFakeRateTriggerPerFile.C+\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshzz4l/sixie/hist/AllNtuple/cern/filefi/028/$f\",\"/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/$f.FakeRateTriggerSkimmed.root\",0\)
// end
//
// cmsLs /store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/028/ >! fileListEGamma
// cmsPfn `cat fileListEGamma | awk '{print $5}'` | sed 's/?svcClass=default//' | sed 's/root:\/\/eoscms\/\/eos\/cms\/store\/group\/phys_egamma\/sixie\/hist\/AllNtuple\/cern\/filefi\/028\///' > ! fileListEGamma2
//
// foreach f(`cat fileListEGamma2 | grep HZZ4lNtuple | grep r12b-del-pr-v1`)
// root -l -b -q HiggsAna/HZZ4l/Skimming/SkimFakeRateTriggerPerFile.C+\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/028/$f\",\"/tmp/sixie/$f.FakeRateTriggerSkimmed.root\",1\)
// end
//
// cmsLs /store/group/phys_muon/sixie/hist/AllNtuple/cern/filefi/028/ >! fileListMuon
// cmsPfn `cat fileListMuon | awk '{print $5}'` | sed 's/?svcClass=default//' | sed 's/root:\/\/eoscms\/\/eos\/cms\/store\/group\/phys_muon\/sixie\/hist\/AllNtuple\/cern\/filefi\/028\///' > ! fileListMuon2
//
// foreach f(`cat fileListMuon2 | grep HZZ4lNtuple | grep r12b-dmu-pr-v1`)
// root -l -b -q HiggsAna/HZZ4l/Skimming/SkimFakeRateTriggerPerFile.C+\(\"root://eoscms//eos/cms/store/group/phys_muon/sixie/hist/AllNtuple/cern/filefi/028/$f\",\"/tmp/sixie/$f.FakeRateTriggerSkimmed.root\",0\)
// end
// *************************************
// From BLUE
//
//foreach f ( /data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-del-pr-v1_noskim_????.root )
//root -l -b -q HiggsAna/HZZ4l/Skimming/SkimFakeRateTriggerPerFile.C+\(\"$f\",\"$f.FakeRateTriggerSkimmed.root\",1\)
//end
//
// *************************************


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TVector3.h>               
#include <TChain.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TBenchmark.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <TMath.h>
#endif

// define structures to read in ntuple
#include "HiggsAna/DataTree/interface/HiggsAnaDefs.hh"
#include "HiggsAna/DataTree/interface/TEventInfo.hh"
#include "HiggsAna/DataTree/interface/TMuon.hh"
#include "HiggsAna/DataTree/interface/TElectron.hh"
#include "HiggsAna/DataTree/interface/TJet.hh"
#include "HiggsAna/DataTree/interface/TPhoton.hh"
#include "HiggsAna/DataTree/interface/TPFCandidate.hh"
#include "HiggsAna/DataTree/interface/TGenParticle.hh"
// #include "HiggsAna/Utils/LeptonIDCuts.hh"

//=== FUNCTION DECLARATIONS ======================================================================================
Bool_t passHLT(ULong_t triggerBits, UInt_t runNum, Int_t TriggerSelection);
Bool_t passElectronDenominatorCuts(const higgsana::TElectron *ele);
Bool_t passMuonDenominatorCuts(const higgsana::TMuon *mu);

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

Double_t getNEvents(string filename) {

  //Get Number of Events in the Sample
  TFile *file = new TFile(filename.c_str(),"READ");
  if (!file) {
    cout << "Could not open file " << filename << endl;
    return 0;
  }

  TDirectory *dir = (TDirectory*)file->FindObjectAny("AnaFwkMod");
  if (!dir) {
    cout << "Could not find directory AnaFwkMod"
         << " in file " << filename << endl;
    delete file;
    return 0;
  }

  TH1F *hist = (TH1F*)dir->Get("hDAllEvents");
  if (!hist) {
    cout << "Could not find histogram hDEvents in directory AnaFwkMod"
         << " in file " << filename << endl;
    delete dir;
    delete file;
    return 0;
  }  
  return hist->Integral();
}

// Main macro function
//--------------------------------------------------------------------------------------------------
void SkimFakeRateTriggerPerFile(string inputFilename, string outputFilename, Int_t TriggerSelection) 
{
  gBenchmark->Start("SkimNtuples");
    
  TTree::SetMaxTreeSize(kMaxLong64);
  
  // Don't write TObject part of the objects
  higgsana::TEventInfo::Class()->IgnoreTObjectStreamer();
  higgsana::TMuon::Class()->IgnoreTObjectStreamer();
  higgsana::TElectron::Class()->IgnoreTObjectStreamer();
  higgsana::TJet::Class()->IgnoreTObjectStreamer();
  higgsana::TPhoton::Class()->IgnoreTObjectStreamer();

  // Data structures to store info from TTrees
  higgsana::TEventInfo *info  = new higgsana::TEventInfo();
  TClonesArray *muonArr     = new TClonesArray("higgsana::TMuon");
  TClonesArray *electronArr     = new TClonesArray("higgsana::TElectron");
  TClonesArray *jetArr    = new TClonesArray("higgsana::TJet");
  TClonesArray *photonArr     = new TClonesArray("higgsana::TPhoton");
  TClonesArray *pfcandidateArr     = new TClonesArray("higgsana::TPFCandidate");
 
  UInt_t nInputEvts = 0;
  UInt_t nPassEvts  = 0;
  UInt_t nEventsTotal = 0;

  TFile* outfile = new TFile(outputFilename.c_str(), "RECREATE");
  
  //
  // Initialize data trees and structs
  // 
  TTree *outEventTree = new TTree("Events","Events"); 
  outEventTree->Branch("Info",     &info);
  outEventTree->Branch("Muon",     &muonArr);
  outEventTree->Branch("Electron", &electronArr);
  outEventTree->Branch("PFJet",    &jetArr);
  outEventTree->Branch("Photon",   &photonArr);
  outEventTree->Branch("PFCandidate",   &pfcandidateArr);


  // list input ntuple files to be skimmed
  
  cout << "Skimming " << inputFilename << endl;
  TTree *eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
  nEventsTotal += getNEvents(inputFilename.c_str()); 
  assert(eventTree);
    
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Muon",     &muonArr);     TBranch *muonBr     = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("PFJet",    &jetArr);      TBranch *jetBr      = eventTree->GetBranch("PFJet");
  eventTree->SetBranchAddress("Photon",    &photonArr);  TBranch *photonBr   = eventTree->GetBranch("Photon");
  eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr);  TBranch *pfcandidateBr   = eventTree->GetBranch("PFCandidate");
     
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) { 
    infoBr->GetEntry(ientry);
    muonArr->Clear();     muonBr->GetEntry(ientry);
    electronArr->Clear(); electronBr->GetEntry(ientry);      
    jetArr->Clear(); jetBr->GetEntry(ientry);
    photonArr->Clear(); photonBr->GetEntry(ientry);
    pfcandidateArr->Clear(); pfcandidateBr->GetEntry(ientry);

    if (ientry % 100000 == 0) cout << "Events: " << ientry << endl;

    nInputEvts++;
      
    Bool_t keep = kTRUE;
    if (!passHLT(info->triggerBits, info->runNum, TriggerSelection)) keep = kFALSE;

    TVector3 met;        
    if(info->pfMEx!=0 || info->pfMEy!=0) {       
      met.SetXYZ(info->pfMEx, info->pfMEy, 0);
    }
//     if (met.Pt() > 20) keep = kFALSE;

//     //events with only 1 reco lepton
//     if(TriggerSelection == 0) {
//       Int_t NDenominatorMuons = 0;
//       Int_t NRecoMuons = 0;
//       for (UInt_t i=0; i< UInt_t(muonArr->GetEntries()) ; ++i) {
//         const higgsana::TMuon *mu = (higgsana::TMuon*)((*muonArr)[i]);              
//         if (mu->pt > 10) NRecoMuons++;
//         if (passMuonDenominatorCuts(mu)) NDenominatorMuons++;
//       }
        
//       if (NRecoMuons > 1) keep = kFALSE;
//       if (!(NDenominatorMuons == 1 )) keep = kFALSE;
//     }
//     if(TriggerSelection == 1) {
//       Int_t NDenominatorElectrons = 0;
//       Int_t NRecoElectrons = 0;
//       for (UInt_t i=0; i< UInt_t(electronArr->GetEntries()) ; ++i) {
//         const higgsana::TElectron *ele = (higgsana::TElectron*)((*electronArr)[i]);  
//         if (ele->pt > 10) NRecoElectrons++;
//         if (passElectronDenominatorCuts(ele)) NDenominatorElectrons++;
//       }
        
//       if (NRecoElectrons > 1) keep = kFALSE;
//       if (!(NDenominatorElectrons == 1 )) keep = kFALSE;
//     }


    if(keep) {
      outEventTree->Fill();
      nPassEvts++;
    }
  }
  outfile->Write();
  outfile->Close();
  
  delete info;
  delete muonArr;
  delete electronArr;
  delete jetArr;
    
  std::cout << outputFilename << " created!" << std::endl;
  std::cout << " >>> Total Number of Events: " << nEventsTotal << std::endl;
  std::cout << " >>> Events processed: " << nInputEvts << std::endl;
  std::cout << " >>>   Events passing: " << nPassEvts << std::endl;

  gBenchmark->Show("SkimNtuples");
}  

Bool_t passHLT(ULong_t triggerBits, UInt_t runNum, Int_t TriggerSelection) {


  Bool_t pass = kFALSE;
  if (TriggerSelection == -1 || TriggerSelection == 0) {
    if ( (triggerBits & kHLT_Mu8) == kHLT_Mu8 ) pass = kTRUE;
    if ( (triggerBits & kHLT_Mu15) == kHLT_Mu15 ) pass = kTRUE;    
    if ( (triggerBits & kHLT_Mu17) == kHLT_Mu17 ) pass = kTRUE;    
  }
  if (TriggerSelection == -1 || TriggerSelection == 1) {
    //it's electron data      
    if ( (triggerBits & kHLT_Ele8) == kHLT_Ele8) pass = kTRUE;
    if ( (triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL) == kHLT_Ele8_CaloIdL_CaloIsoVL) pass = kTRUE;
    if ( (triggerBits & kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL) == kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL) pass = kTRUE;
    if ( (triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40) == kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40) pass = kTRUE;
    if ( (triggerBits & kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30) == kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30) pass = kTRUE;
    if ( (triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL) == kHLT_Ele17_CaloIdL_CaloIsoVL) pass = kTRUE;
    if ( (triggerBits & kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL) == kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL) pass = kTRUE;
    if ( (triggerBits & kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30) == kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30) pass = kTRUE;

  }
  return pass;
}


Bool_t passElectronDenominatorCuts( const higgsana::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  // Use Looser Denominator

//   //Barrel 
//   if (fabs(ele->eta) < 1.5) {
//     if (! ( (0==0)
//             && fabs(ele->dz) < 0.1
//           )
//       ) {
//       pass = kFALSE;
//     }      
//   }
//   //Endcap
//   else if (fabs(ele->eta) > 1.5) {
//     if (! (  (0==0)
//              && fabs(ele->dz) < 0.1
//           )
//       ) {
//       pass = kFALSE;
//     }
//   } else {
//     pass = kFALSE;
//     return pass;
//   }
  
  return pass;
}



Bool_t passMuonDenominatorCuts(const higgsana::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  if (!(
        fabs(mu->dz) < 0.1
        )
    ) {
    pass = kFALSE;
  }
 
  return pass;
}
