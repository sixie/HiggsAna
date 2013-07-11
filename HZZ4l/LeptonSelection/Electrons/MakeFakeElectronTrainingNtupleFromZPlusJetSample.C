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
#include "HiggsAna/DataTree/interface/Types.h"
#include "HiggsAna/Utils/LeptonIDCuts.hh"
#include "HiggsAna/Utils/LeptonTools.hh"
#include "HiggsAna/Utils/CommonTools.hh"
#include "HiggsAna/CommonData/interface/ElectronTree.h"
#include "HiggsAna/HZZ4l/Utils/LeptonSelection.hh"

#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

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
void MakeFakeElectronTrainingNtupleFromZPlusJetSample(Int_t Option = 0) {

//   MakeNtuple("LIST","ElectronSelectionTraining.Fake_ZPlusJet.root",-1);
  if (Option == 0) {
    MakeNtuple("2012Data","/data/blue/sixie/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_ZPlusJet_2012.root", Option);
  }
  if (Option == 1) {
    MakeNtuple("2012Data","/data/blue/sixie/LeptonSelection/ElectronsNoRhoCorr/ElectronSelectionTraining.Fake_ZPlusJet_2012.root", Option);
  }

}


void MakeNtuple(const string inputFilename, const string outputFilename)
{  
  gBenchmark->Start("WWTemplate");

  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  TFile *outputFile = new TFile(outputFilename.c_str(), "RECREATE");
  ElectronTree *eleTree = new ElectronTree;
  eleTree->CreateTree();
  eleTree->tree_->SetAutoFlush(0);

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

  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("/data/smurf/data/Winter11_4700ipb/auxiliar/hww.Full2011.json"); 
  rlrm.AddJSONFile("/data/blue/sixie/HZZ4l/auxiliar/2012/Cert_190456-196531_8TeV_PromptReco_Collisions12_JSON.txt");

  Int_t NEvents = 0;

  UInt_t DataEra = kDataEra_NONE;

  vector<string> inputfiles;
  if (inputFilename == "LIST") {
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-del-m10-v1_noskim.ZPlusFakeSkimmed.root");  
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-del-pr-v4_noskim.ZPlusFakeSkimmed.root");  
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-del-a05-v1_noskim.ZPlusFakeSkimmed.root");  
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-del-o03-v1_noskim.ZPlusFakeSkimmed.root");  
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11b-del-pr-v1_noskim.ZPlusFakeSkimmed.root");  
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-smu-m10-v1_noskim.ZPlusFakeSkimmed.root");  
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-smu-pr-v4_noskim.ZPlusFakeSkimmed.root");  
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-smu-a05-v1_noskim.ZPlusFakeSkimmed.root");  
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-smu-o03-v1_noskim.ZPlusFakeSkimmed.root");  
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11b-smu-pr-v1_noskim.ZPlusFakeSkimmed.root");
    DataEra = kDataEra_2011_MC;
  } else if (inputFilename == "2012Data") {
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r12a-del-pr-v1_ZPlusFakeSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r12a-dmu-pr-v1_ZPlusFakeSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r12b-del-pr-v1_ZPlusFakeSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r12b-dmu-pr-v1_ZPlusFakeSkimmed.root");
    DataEra = kDataEra_2012_MC;
  } else {
    inputfiles.push_back(inputFilename);
  }

  //DataEra
  if (Option == 1) {
    DataEra = kDataEra_NONE;    
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

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      Double_t eventweight = info->eventweight;

      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...

      NEvents++;

      //********************************************************
      // Load the branches
      //********************************************************
      electronArr->Clear(); 
      muonArr->Clear(); 
      photonArr->Clear(); 
      jetArr->Clear(); 
      pfcandidateArr->Clear(); 
      electronBr->GetEntry(ientry);
      muonBr->GetEntry(ientry);
      photonBr->GetEntry(ientry);
      jetBr->GetEntry(ientry);
      pfcandidateBr->GetEntry(ientry);

      //********************************************************
      // Pileup Energy Density
      //********************************************************
      Double_t rhoMuonIso = 0;
      Double_t rhoEleIso = 0;
      UInt_t MuonEAEra = 0;
      UInt_t EleEAEra = 0;

      if (DataEra == kDataEra_2011_MC) {     
        if (!(isnan(info->RhoKt6PFJetsForIso25) || 
              isinf(info->RhoKt6PFJetsForIso25))) {
          rhoMuonIso = info->RhoKt6PFJetsForIso25;
          rhoEleIso = info->RhoKt6PFJetsForIso25;
        }
        EleEAEra = kDataEra_2011_Data;
      } else if (DataEra == kDataEra_2012_MC) {
        if (!(isnan(info->RhoKt6PFJetsCentralNeutral) || 
              isinf(info->RhoKt6PFJetsCentralNeutral))) {
          rhoMuonIso = info->RhoKt6PFJetsCentralNeutral;
        }
        if (!(isnan(info->RhoKt6PFJets) || 
              isinf(info->RhoKt6PFJets))) {
          rhoEleIso = info->RhoKt6PFJets;
        }
        MuonEAEra = kDataEra_2012_Data;
        EleEAEra = kDataEra_2012_Data;
      }



      //********************************************************
      // TcMet
      //********************************************************
      TVector3 pfMet;        
      if(info->pfMEx!=0 || info->pfMEy!=0) {       
        pfMet.SetXYZ(info->pfMEx, info->pfMEy, 0);
      }
      Double_t met = pfMet.Pt();

      //********************************************************
      // Z-Selection
      //********************************************************
      Int_t NElectrons = 0;
      vector<const higgsana::TMuon*> goodMuons;
      vector<const higgsana::TElectron*> goodElectrons;
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const higgsana::TMuon *mu = (higgsana::TMuon*)((*muonArr)[i]);
        if( mu->pt > 20.0 && 
          passMuonID(mu, ComputeMuonPFIso04(mu,pfcandidateArr,rhoMuonIso,MuonEAEra))
          ) {
          goodMuons.push_back(mu);
        }
      }
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const higgsana::TElectron *ele = (higgsana::TElectron*)((*electronArr)[i]);
        if (ele->pt > 5) NElectrons++;
        if( ele->pt > 20.0 && 
            passCutBasedEleID(ele,ComputeElePFIso04(ele,pfcandidateArr,rhoEleIso,EleEAEra))
          ) {
          goodElectrons.push_back(ele);
        }
      }

      Int_t NZCandidates = 0;
      const higgsana::TMuon* ZMuon1 = 0;
      const higgsana::TMuon* ZMuon2 = 0;
      const higgsana::TElectron *ZEle1 = 0;
      const higgsana::TElectron *ZEle2 = 0;
      Double_t ZPt = 0;
      Double_t ZMass = 0;
      for(Int_t i=0; i<goodMuons.size(); i++) {
        higgsana::FourVectorM mu1;
        mu1.SetCoordinates(goodMuons[i]->pt, goodMuons[i]->eta, goodMuons[i]->phi, 105.658369e-3 );
        for(Int_t j=i+1; j<goodMuons.size(); j++) {
          higgsana::FourVectorM mu2;
          mu2.SetCoordinates(goodMuons[j]->pt, goodMuons[j]->eta, goodMuons[j]->phi, 105.658369e-3 );
          higgsana::FourVectorM dilepton = mu1+mu2;
          if (dilepton.M() > 80 && dilepton.M() < 105) {
            ZMuon1 = goodMuons[i];
            ZMuon2 = goodMuons[j];
            ZPt = dilepton.Pt();
            ZMass = dilepton.M();
            NZCandidates++;
          }
        }
      }
      for(Int_t i=0; i<goodElectrons.size(); i++) {
        higgsana::FourVectorM ele1;
        ele1.SetCoordinates(goodElectrons[i]->pt, goodElectrons[i]->eta, goodElectrons[i]->phi, 0.51099892e-3 );
        for(Int_t j=i+1; j<goodElectrons.size(); j++) {
          higgsana::FourVectorM ele2;
          ele2.SetCoordinates(goodElectrons[j]->pt, goodElectrons[j]->eta, goodElectrons[j]->phi, 0.51099892e-3 );
          higgsana::FourVectorM dilepton = ele1+ele2;
          if (dilepton.M() > 80 && dilepton.M() < 105) {
            ZEle1 = goodElectrons[i];
            ZEle2 = goodElectrons[j];
            ZPt = dilepton.Pt();
            ZMass = dilepton.M();
            NZCandidates++;
          }
        }
      }

      //********************************************************
      // Define Z Lepton Type
      //********************************************************
      if (NZCandidates != 1) continue;
      Int_t ZDecayType = 0;
      if (ZMuon1 && ZMuon2) ZDecayType = 13;
      else if (ZEle1 && ZEle2) ZDecayType = 11;
      else cout << "Error. Z Leptons not properly found.\n";            

      //********************************************************
      // Veto ZZ events
      //********************************************************
      if ((ZDecayType == 11 &&  NElectrons > 3) 
          ||
          (ZDecayType == 13 &&  NElectrons > 1) 
        ) continue;
      
      //MET Cut
      if (met > 20) continue;

//       cout << ZDecayType << " : ";
//       if ( ZEle1 ) cout << ZEle1->pt << " " << ZEle1->eta << " " << ZEle1->phi << " :";
//       if ( ZEle2 ) cout << ZEle2->pt << " " << ZEle2->eta << " " << ZEle2->phi << " :";
//       cout << endl;

      //******************************************************************************
      //loop over electrons 
      //******************************************************************************
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const higgsana::TElectron *ele = (higgsana::TElectron*)((*electronArr)[i]);

//         cout << ele->pt << " " << ele->eta << " " << ele->phi << " : " ;
//         if (ZEle1)
//           cout << (ele == ZEle1) << " " << deltaR(ele->eta, ele->phi, ZEle1->eta, ZEle1->phi) << " : ";
//         if (ZEle2)
//           cout << (ele == ZEle2) << " " << deltaR(ele->eta, ele->phi, ZEle2->eta, ZEle2->phi) << " : ";
//         cout << endl;

//         //Remove Z Electrons
//         if (ZDecayType == 11 && (ele == ZEle1 || ele == ZEle2)) continue;

        //make cut on dz
        if (fabs(ele->dz) > 0.1) continue;
        
        //veto against muons
        if (ZEle1 && deltaR(ele->eta, ele->phi, ZEle1->eta, ZEle1->phi) < 0.3) continue;
        if (ZEle2 && deltaR(ele->eta, ele->phi, ZEle2->eta, ZEle2->phi) < 0.3) continue;
        if (ZMuon1 && deltaR(ele->eta, ele->phi, ZMuon1->eta, ZMuon1->phi) < 0.3) continue;
        if (ZMuon2 && deltaR(ele->eta, ele->phi, ZMuon2->eta, ZMuon2->phi) < 0.3) continue;

        //***********************
        //Fill Electron
        //***********************
        NElectronsFilled++;
        FillElectronTree( eleTree, ele, pfcandidateArr, rhoEleIso, EleEAEra, 
                          info->nPV0, info->runNum, info->lumiSec, info->evtNum);

        
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


