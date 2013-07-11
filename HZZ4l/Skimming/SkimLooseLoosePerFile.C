//foreach f ( /home/sixie/hist/AllNtuple/HWWNtuple/data/unprocessed/AllNtuple_HWWNtuple_r11a-del-pr-v4_noskim_????.root )
//root -l -b -q EWKAna/Hww/Skim/SkimLooseLoosePerFile.C+\(\"$f\",\"$f.LooseLooseSkimmed.root\"\)
//end



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
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TJet.hh"
#include "EWKAna/Ntupler/interface/TPhoton.hh"
#include "EWKAna/Utils/LeptonIDCuts.hh"

//=== FUNCTION DECLARATIONS ======================================================================================
Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu, Double_t fRho);
Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele);

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
      cout << "Cannot get Directory ZeeAnalysisMod from file " << infname << endl;
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
void SkimLooseLoosePerFile(string inputFilename, string outputFilename) 
{
  gBenchmark->Start("SkimNtuples");
    
  TTree::SetMaxTreeSize(kMaxLong64);
  
  // Don't write TObject part of the objects
  mithep::TEventInfo::Class()->IgnoreTObjectStreamer();
  mithep::TMuon::Class()->IgnoreTObjectStreamer();
  mithep::TElectron::Class()->IgnoreTObjectStreamer();
  mithep::TJet::Class()->IgnoreTObjectStreamer();
  mithep::TPhoton::Class()->IgnoreTObjectStreamer();

  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  TClonesArray *muonArr     = new TClonesArray("mithep::TMuon");
  TClonesArray *electronArr     = new TClonesArray("mithep::TElectron");
  TClonesArray *jetArr    = new TClonesArray("mithep::TJet");
  TClonesArray *photonArr     = new TClonesArray("mithep::TPhoton");
    
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
     
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) { 
    infoBr->GetEntry(ientry);
    muonArr->Clear();     muonBr->GetEntry(ientry);
    electronArr->Clear(); electronBr->GetEntry(ientry);      
    jetArr->Clear(); jetBr->GetEntry(ientry);
    photonArr->Clear(); photonBr->GetEntry(ientry);
     
    if (ientry % 100000 == 0) cout << "Events: " << ientry << endl;

    nInputEvts++;
      
    Bool_t keep = kTRUE;
    
    Double_t rho = 0;
    if (!(TMath::IsNaN(info->PileupEnergyDensity) || isinf(info->PileupEnergyDensity))) rho = info->PileupEnergyDensity;


    Int_t NDenominatorMuons = 0;
    for(Int_t i=0; i<muonArr->GetEntries(); i++) {
      const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);      
      if ( passMuonDenominatorCuts(mu, rho)
           &&
           fabs(mu->eta) < 2.4
           && 
           mu->pt > 10.0
        ) {
        NDenominatorMuons++;
      }
    }
    Int_t NDenominatorElectrons = 0;
    for(Int_t i=0; i<electronArr->GetEntries(); i++) {
      const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
      if ( passElectronDenominatorCuts(ele)
           &&
           fabs(ele->eta) < 2.5
           && 
           ele->pt > 10.0
        ) {
        NDenominatorElectrons++;
      }
    }

    //One Tight AND Two Loose
    if (!(NDenominatorElectrons + NDenominatorMuons >= 2 )) {
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





Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(ele->pt > 10 && fabs(ele->eta) < 2.5)) pass = kFALSE;


  //Barrel 
  if (fabs(ele->scEta) < 1.479) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.007
            && fabs(ele->deltaPhiIn) < 0.15
            && ele->HoverE < 0.12
            && ele->nExpHitsInner <= 0
            && passConversionVeto(ele->isConv)
            && fabs(ele->dz) < 0.1
            && (ele->trkIso03) / ele->pt < 0.20
            && (ele->emIso03) / ele->pt < 0.20
            && (ele->hadIso03) / ele->pt < 0.20
              
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else {
    if (! (  (0==0)
             && ele->sigiEtaiEta < 0.03
             && fabs(ele->deltaEtaIn) < 0.009
             && fabs(ele->deltaPhiIn) < 0.10
            && ele->HoverE < 0.10
             && ele->nExpHitsInner <= 0
             && passConversionVeto(ele->isConv)
             && fabs(ele->dz) < 0.1
             && (ele->trkIso03) / ele->pt < 0.20
             && (ele->emIso03) / ele->pt < 0.20
             && (ele->hadIso03) / ele->pt < 0.20
          )
      ) {
      pass = kFALSE;
    }
  } 

  return pass;
}




Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu, Double_t fRho) {
  
  Bool_t pass = kTRUE;

  if (! 
      ( (
          (Bool_t(mu->typeBits & kGlobal) 
           && mu->muNchi2 < 10.0
           && (mu->nValidHits > 0)
           && (mu->nMatch > 1 )
         )
           || 
          ( mu->typeBits & kTracker            
          && Bool_t(mu->qualityBits & kTMLastStationTight) 
            )
        )
        && mu->typeBits & kTracker
        && mu->nTkHits > 10
        && (mu->nPixHits > 0)
        && fabs(mu->d0) < 0.2
        && fabs(mu->dz) < 0.1
        && (
          ((mu->ChargedIso03 + mu->NeutralIso03_05Threshold - fRho*MuonEffectiveArea(kMuNeutralIso03,mu->eta)) / mu->pt < 1.0)
          ||
          ((mu->trkIso03 + mu->emIso03 + mu->hadIso03 
            - fRho*MuonEffectiveArea(kMuEMIso03,mu->eta) 
            - fRho*MuonEffectiveArea(kMuHadIso03,mu->eta)) / mu->pt < 1.0
            )
          )
        && (mu->pterr / mu->pt < 0.1)
        && mu->TrkKink  < 20
        )
    ) pass = kFALSE;

 
  return pass;
}
