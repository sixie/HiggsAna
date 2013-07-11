// *************************************
// From EOS
//
// cmsLs /store/user/sixie/hist/AllNtuple/cern/filefi/028/ >! fileList
// cmsPfn `cat fileList | awk '{print $5}'` | sed 's/?svcClass=default//' | sed 's/root:\/\/eoscms\/\/eos\/cms\/store\/user\/sixie\/hist\/AllNtuple\/cern\/filefi\/028\///' > ! fileList2
//
// cmsLs /store/group/phys_higgs/cmshzz4l/sixie/hist/AllNtuple/cern/filefi/028/ >! fileListB
// cmsPfn `cat fileListB | awk '{print $5}'` | sed 's/?svcClass=default//' | sed 's/root:\/\/eoscms\/\/eos\/cms\/store\/group\/phys_higgs\/cmshzz4l\/sixie\/hist\/AllNtuple\/cern\/filefi\/028\///' > ! fileListB2
//
// foreach f(`cat fileListA2 | grep HZZ4lNtuple | grep r12a-del-pr-v1`)
// root -l -b -q HiggsAna/HZZ4l/Skimming/SkimTightPlusRecoPerFile.C+\(\"root://eoscms//eos/cms/store/user/sixie/hist/AllNtuple/cern/filefi/028/$f\",\"/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/$f.TightPlusRecoSkimmed.root\"\)
// end
// foreach f(`cat fileListEGamma2 | grep HZZ4lNtuple | grep r12b-del-pr-v1`)
// root -l -b -q HiggsAna/HZZ4l/Skimming/SkimTightPlusRecoPerFile.C+\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/028/$f\",\"/tmp/sixie/$f.TightPlusRecoSkimmed.root\"\)
// end
// foreach f(`cat fileListMuon2 | grep HZZ4lNtuple | grep r12b-smu-pr-v1`)
// root -l -b -q HiggsAna/HZZ4l/Skimming/SkimTightPlusRecoPerFile.C+\(\"root://eoscms//eos/cms/store/group/phys_muon/sixie/hist/AllNtuple/cern/filefi/028/$f\",\"/tmp/sixie/$f.TightPlusRecoSkimmed.root\"\)
// end
// foreach f(`cat fileListEGamma2 | grep HZZ4lNtuple | grep r12b-sel-pr-v1`)
// root -l -b -q HiggsAna/HZZ4l/Skimming/SkimTightPlusRecoPerFile.C+\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/028/$f\",\"/tmp/sixie/$f.TightPlusRecoSkimmed.root\"\)
// end
//
// *************************************
// From BLUE
//
//foreach f ( /data/blue/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r11a-del-m10-v1_noskim_????.root )
//root -l -b -q HiggsAna/HZZ4l/Skimming/SkimTightPlusRecoPerFile.C+\(\"$f\",\"$f.TightPlusRecoSkimmed.root\"\)
//end
//
// *************************************



#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
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
#endif

// define structures to read in ntuple
#include "HiggsAna/DataTree/interface/HiggsAnaDefs.hh"
#include "HiggsAna/DataTree/interface/TEventInfo.hh"
#include "HiggsAna/DataTree/interface/TMuon.hh"
#include "HiggsAna/DataTree/interface/TElectron.hh"
#include "HiggsAna/DataTree/interface/TJet.hh"
#include "HiggsAna/DataTree/interface/TPhoton.hh"
#include "HiggsAna/DataTree/interface/TPFCandidate.hh"
#include "HiggsAna/Utils/LeptonIDCuts.hh"


//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  cout << "--- Open file " << infname << endl;
  
  TFile* inf = TFile::Open(infname,"read");
  assert(inf);

  TTree* t = (TTree*)inf->Get(tname);
  
//   if (!t) {
//     TDirectory *dir = (TDirectory*)inf->FindObjectAny("HwwNtuplerMod");
//     if (!dir) {
//       cout << "Cannot get Directory ZeeAnalysisMod from file " << infname << endl;
//       assert(dir);
//     }
//     t = (TTree*)dir->Get(tname);
//   }

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
void SkimTightPlusRecoPerFile(string inputFilename, string outputFilename) 
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

  cout << "Skimming " << inputFilename << "..." << endl;
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
    
    //****************************************************
    //Rho for isolation pileup correction
    //****************************************************
    Double_t rhoMuonIso = 0;
    Double_t rhoEleIso = 0;
    UInt_t MuonEAEra = kDataEra_2012_Data;
    UInt_t EleEAEra = kDataEra_2012_Data;
    if (!(isnan(info->RhoKt6PFJetsCentralNeutral) || 
          isinf(info->RhoKt6PFJetsCentralNeutral))) {
      rhoMuonIso = info->RhoKt6PFJetsCentralNeutral;
    }    
    if (!(isnan(info->RhoKt6PFJets) || 
          isinf(info->RhoKt6PFJets))) {
      rhoEleIso = info->RhoKt6PFJets;
    }  


    Int_t NTightMuons = 0;
    Int_t NRecoMuons = 0;
    for(Int_t i=0; i<muonArr->GetEntries(); i++) {
      const higgsana::TMuon *mu = (higgsana::TMuon*)((*muonArr)[i]);      
      if ( fabs(mu->eta) < 2.4
           && 
           mu->pt > 5.0
        ) {
        NRecoMuons++;

        if ( passMuonID(mu, ComputeMuonPFIso04(mu,pfcandidateArr,rhoMuonIso,MuonEAEra))) NTightMuons++;
      }
    }
    Int_t NTightElectrons = 0;
    Int_t NRecoElectrons = 0;
    for(Int_t i=0; i<electronArr->GetEntries(); i++) {
      const higgsana::TElectron *ele = (higgsana::TElectron*)((*electronArr)[i]);
      if ( fabs(ele->eta) < 2.5
           && 
           ele->pt > 5.0
        ) {
        NRecoElectrons++;
        if (passCutBasedEleID(ele,ComputeElePFIso04(ele,pfcandidateArr,rhoEleIso,EleEAEra))) NTightElectrons++;
      }
    }

    //One Tight AND Two Loose
    if (!( NTightElectrons + NTightMuons >= 1  
           && NRecoElectrons + NRecoMuons >= 2
          ) ) {
      keep = kFALSE;
    }
      
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

