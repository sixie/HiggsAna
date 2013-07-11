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
#include "HiggsAna/Ntupler/interface/HiggsAnaDefs.hh"
#include "HiggsAna/Ntupler/interface/TEventInfo.hh"
#include "HiggsAna/Ntupler/interface/TElectron.hh"
#include "HiggsAna/Ntupler/interface/TPhoton.hh"
#include "HiggsAna/Ntupler/interface/TMuon.hh"
#include "HiggsAna/Ntupler/interface/TJet.hh"
#include "HiggsAna/Ntupler/interface/TPFCandidate.hh"
#include "HiggsAna/Utils/LeptonIDCuts.hh"

#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "HiggsAna/HZZ4l/Utils/LeptonSelection.hh"

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
Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu);
void MakeNtuple(const string inputFilename,  const string outputFilename, Int_t TreeType);

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
void MakeFakeMuonTrainingNtupleFromZPlusJetSample(Int_t type = -1) {

  if (type == -1) {
    MakeNtuple("LIST","MuonSelection.Fake_ZPlusJet.root", -1);
  }

}


void MakeNtuple(const string inputFilename, const string outputFilename, Int_t TreeType)
{  
  gBenchmark->Start("WWTemplate");

  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  TFile *outputFile = new TFile(outputFilename.c_str(), "RECREATE");
  TTree *muTree = new TTree("Muons","Muons");
  muTree->SetAutoFlush(0);

  //Variables
  Float_t                 fWeight;
  UInt_t                  fRunNumber;
  UInt_t                  fLumiSectionNumber;
  UInt_t                  fEventNumber;
  Float_t                 fMuPt; 
  Float_t                 fMuEta; 
  Float_t                 fMuPhi; 
  Float_t                 fMuPFIso; 
  
  //CutBased Variables
  Float_t                 fMuTkNchi2; 
  Float_t                 fMuGlobalNchi2; 
  Float_t                 fMuNValidHits; 
  Float_t                 fMuNTrackerHits; 
  Float_t                 fMuNPixelHits; 
  Float_t                 fMuNMatches; 
  Float_t                 fMuD0; 

  //Additional Vars used in Likelihood
  Float_t                 fMuIP3d; 
  Float_t                 fMuIP3dSig; 
  Float_t                 fMuTrkKink; 
  Float_t                 fMuGlobalKink; 
  Float_t                 fMuSegmentCompatibility; 
  Float_t                 fMuCaloCompatibility; 
  Float_t                 fMuHadEnergy; 
  Float_t                 fMuHoEnergy; 
  Float_t                 fMuEmEnergy; 
  Float_t                 fMuHadS9Energy; 
  Float_t                 fMuHoS9Energy; 
  Float_t                 fMuEmS9Energy; 

  //Isolation Variables
  Float_t                 fMuChargedIso03; 
  Float_t                 fMuChargedIso03FromOtherVertices; 
  Float_t                 fMuNeutralIso03_05Threshold; 
  Float_t                 fMuNeutralIso03_10Threshold; 
  Float_t                 fMuChargedIso04; 
  Float_t                 fMuChargedIso04FromOtherVertices; 
  Float_t                 fMuNeutralIso04_05Threshold; 
  Float_t                 fMuNeutralIso04_10Threshold; 
  Float_t                 fMuTrkIso03; 
  Float_t                 fMuEMIso03; 
  Float_t                 fMuHadIso03; 
  Float_t                 fMuTrkIso05; 
  Float_t                 fMuEMIso05; 
  Float_t                 fMuHadIso05; 
  Float_t                 fRho; 
  Float_t                 fNVertices; 

  //More Isolation Variables
  Float_t fPFIso_MeanEta;
  Float_t fPFIso_MeanPhi;
  Float_t fPFIso_SigmaEtaEta;
  Float_t fPFIso_SigmaPhiPhi;
  Float_t fPFIso_SigmaEtaPhi;
  Float_t fChargedIso_DR0p0To0p1;
  Float_t fChargedIso_DR0p1To0p2;
  Float_t fChargedIso_DR0p2To0p3;
  Float_t fChargedIso_DR0p3To0p4;
  Float_t fChargedIso_DR0p4To0p5;
  Float_t fChargedIso_DR0p5To0p7;
  Float_t fChargedIso_DR0p7To1p0;
  Float_t fChargedIso_MeanEta;
  Float_t fChargedIso_MeanPhi;
  Float_t fChargedIso_SigmaEtaEta;
  Float_t fChargedIso_SigmaPhiPhi;
  Float_t fChargedIso_SigmaEtaPhi;
  Float_t fGammaIso_DR0p0To0p1;
  Float_t fGammaIso_DR0p1To0p2;
  Float_t fGammaIso_DR0p2To0p3;
  Float_t fGammaIso_DR0p3To0p4;
  Float_t fGammaIso_DR0p4To0p5;
  Float_t fGammaIso_DR0p5To0p7;
  Float_t fGammaIso_DR0p7To1p0;
  Float_t fGammaIso_MeanEta;
  Float_t fGammaIso_MeanPhi;
  Float_t fGammaIso_SigmaEtaEta;
  Float_t fGammaIso_SigmaPhiPhi;
  Float_t fGammaIso_SigmaEtaPhi;
  Float_t fNeutralHadronIso_DR0p0To0p1;
  Float_t fNeutralHadronIso_DR0p1To0p2;
  Float_t fNeutralHadronIso_DR0p2To0p3;
  Float_t fNeutralHadronIso_DR0p3To0p4;
  Float_t fNeutralHadronIso_DR0p4To0p5;
  Float_t fNeutralHadronIso_DR0p5To0p7;
  Float_t fNeutralHadronIso_DR0p7To1p0;
  Float_t fNeutralHadronIso_MeanEta;
  Float_t fNeutralHadronIso_MeanPhi;
  Float_t fNeutralHadronIso_SigmaEtaEta;
  Float_t fNeutralHadronIso_SigmaPhiPhi;
  Float_t fNeutralHadronIso_SigmaEtaPhi;
  Float_t fDirectionalPFIso;
  Float_t fDirectionalChargedIso;
  Float_t fDirectionalGammaIso;
  Float_t fDirectionalNeutralHadronIso;


  muTree->Branch("weight",&fWeight,"weight/F");
  muTree->Branch("run",&fRunNumber,"run/i");
  muTree->Branch("lumi",&fLumiSectionNumber,"lumi/i");
  muTree->Branch("event",&fEventNumber,"event/i");
  muTree->Branch("pt",&fMuPt,"pt/F"); 
  muTree->Branch("eta",&fMuEta,"eta/F"); 
  muTree->Branch("phi",&fMuPhi,"phi/F"); 
  muTree->Branch("pfiso",&fMuPFIso,"pfiso/F"); 
  
  //CutBased Variables
  muTree->Branch("TkNchi2",&fMuTkNchi2,"TkNchi2/F"); 
  muTree->Branch("GlobalNchi2",&fMuGlobalNchi2,"GlobalNchi2/F"); 
  muTree->Branch("NValidHits",&fMuNValidHits,"NValidHits/F"); 
  muTree->Branch("NTrackerHits",&fMuNTrackerHits,"NTrackerHits/F"); 
  muTree->Branch("NPixelHits",&fMuNPixelHits,"NPixelHits/F"); 
  muTree->Branch("NMatches",&fMuNMatches,"NMatches/F"); 
  muTree->Branch("D0",&fMuD0,"D0/F"); 

  //Additional Vars used in Likelihood
  muTree->Branch("IP3d",&fMuIP3d,"IP3d/F"); 
  muTree->Branch("IP3dSig",&fMuIP3dSig,"IP3dSig/F"); 
  muTree->Branch("TrkKink",&fMuTrkKink,"TrkKink/F"); 
  muTree->Branch("GlobalKink",&fMuGlobalKink,"GlobalKink/F"); 
  muTree->Branch("SegmentCompatibility",&fMuSegmentCompatibility,"SegmentCompatibility/F"); 
  muTree->Branch("CaloCompatibility",&fMuCaloCompatibility,"CaloCompatibility/F"); 
  muTree->Branch("HadEnergy",&fMuHadEnergy,"HadEnergy/F"); 
  muTree->Branch("HoEnergy",&fMuHoEnergy,"HoEnergy/F"); 
  muTree->Branch("EmEnergy",&fMuEmEnergy,"EmEnergy/F"); 
  muTree->Branch("HadS9Energy",&fMuHadS9Energy,"HadS9Energy/F"); 
  muTree->Branch("HoS9Energy",&fMuHoS9Energy,"HoS9Energy/F"); 
  muTree->Branch("EmS9Energy",&fMuEmS9Energy,"EmS9Energy/F"); 

  //Isolation Variables
  muTree->Branch("ChargedIso03",&fMuChargedIso03,"ChargedIso03/F"); 
  muTree->Branch("ChargedIso03FromOtherVertices",&fMuChargedIso03FromOtherVertices,"ChargedIso03FromOtherVertices/F"); 
  muTree->Branch("NeutralIso03_05Threshold",&fMuNeutralIso03_05Threshold,"NeutralIso03_05Threshold/F"); 
  muTree->Branch("NeutralIso03_10Threshold",&fMuNeutralIso03_10Threshold,"NeutralIso03_10Threshold/F"); 
  muTree->Branch("ChargedIso04",&fMuChargedIso04,"ChargedIso04/F"); 
  muTree->Branch("ChargedIso04FromOtherVertices",&fMuChargedIso04FromOtherVertices,"ChargedIso04FromOtherVertices/F"); 
  muTree->Branch("NeutralIso04_05Threshold",&fMuNeutralIso04_05Threshold,"NeutralIso04_05Threshold/F"); 
  muTree->Branch("NeutralIso04_10Threshold",&fMuNeutralIso04_10Threshold,"NeutralIso04_10Threshold/F"); 
  muTree->Branch("TrkIso03",&fMuTrkIso03,"TrkIso03/F"); 
  muTree->Branch("EMIso03",&fMuEMIso03,"EMIso03/F"); 
  muTree->Branch("HadIso03",&fMuHadIso03,"HadIso03/F"); 
  muTree->Branch("TrkIso05",&fMuTrkIso05,"TrkIso05/F"); 
  muTree->Branch("EMIso05",&fMuEMIso05,"EMIso05/F"); 
  muTree->Branch("HadIso05",&fMuHadIso05,"HadIso05/F"); 
  muTree->Branch("Rho",&fRho,"Rho/F"); 
  muTree->Branch("NVertices",&fNVertices,"NVertices/F"); 

  //More Isolation Variables
  muTree->Branch("PFIso_MeanEta",&fPFIso_MeanEta,"PFIso_MeanEta/F");
  muTree->Branch("PFIso_MeanPhi",&fPFIso_MeanPhi,"fPFIso_MeanPhi/F");
  muTree->Branch("PFIso_SigmaEtaEta",&fPFIso_SigmaEtaEta,"PFIso_SigmaEtaEta/F");
  muTree->Branch("PFIso_SigmaPhiPhi",&fPFIso_SigmaPhiPhi,"PFIso_SigmaPhiPhi/F");
  muTree->Branch("PFIso_SigmaEtaPhi",&fPFIso_SigmaEtaPhi,"PFIso_SigmaEtaPhi/F");
  muTree->Branch("ChargedIso_DR0p0To0p1",&fChargedIso_DR0p0To0p1,"ChargedIso_DR0p0To0p1/F");
  muTree->Branch("ChargedIso_DR0p1To0p2",&fChargedIso_DR0p1To0p2,"ChargedIso_DR0p1To0p2/F");
  muTree->Branch("ChargedIso_DR0p2To0p3",&fChargedIso_DR0p2To0p3,"ChargedIso_DR0p2To0p3/F");
  muTree->Branch("ChargedIso_DR0p3To0p4",&fChargedIso_DR0p3To0p4,"ChargedIso_DR0p3To0p4/F");
  muTree->Branch("ChargedIso_DR0p4To0p5",&fChargedIso_DR0p4To0p5,"ChargedIso_DR0p4To0p5/F");
  muTree->Branch("ChargedIso_DR0p5To0p7",&fChargedIso_DR0p5To0p7,"ChargedIso_DR0p5To0p7/F");
  muTree->Branch("ChargedIso_DR0p7To1p0",&fChargedIso_DR0p7To1p0,"ChargedIso_DR0p7To1p0/F");
  muTree->Branch("ChargedIso_MeanEta",&fChargedIso_MeanEta,"ChargedIso_MeanEta/F");
  muTree->Branch("ChargedIso_MeanPhi",&fChargedIso_MeanPhi,"fChargedIso_MeanPhi/F");
  muTree->Branch("ChargedIso_SigmaEtaEta",&fChargedIso_SigmaEtaEta,"ChargedIso_SigmaEtaEta/F");
  muTree->Branch("ChargedIso_SigmaPhiPhi",&fChargedIso_SigmaPhiPhi,"ChargedIso_SigmaPhiPhi/F");
  muTree->Branch("ChargedIso_SigmaEtaPhi",&fChargedIso_SigmaEtaPhi,"ChargedIso_SigmaEtaPhi/F");
  muTree->Branch("GammaIso_DR0p0To0p1",&fGammaIso_DR0p0To0p1,"GammaIso_DR0p0To0p1/F");
  muTree->Branch("GammaIso_DR0p1To0p2",&fGammaIso_DR0p1To0p2,"GammaIso_DR0p1To0p2/F");
  muTree->Branch("GammaIso_DR0p2To0p3",&fGammaIso_DR0p2To0p3,"GammaIso_DR0p2To0p3/F");
  muTree->Branch("GammaIso_DR0p3To0p4",&fGammaIso_DR0p3To0p4,"GammaIso_DR0p3To0p4/F");
  muTree->Branch("GammaIso_DR0p4To0p5",&fGammaIso_DR0p4To0p5,"GammaIso_DR0p4To0p5/F");
  muTree->Branch("GammaIso_DR0p5To0p7",&fGammaIso_DR0p5To0p7,"GammaIso_DR0p5To0p7/F");
  muTree->Branch("GammaIso_DR0p7To1p0",&fGammaIso_DR0p7To1p0,"GammaIso_DR0p7To1p0/F");
  muTree->Branch("GammaIso_MeanEta",&fGammaIso_MeanEta,"GammaIso_MeanEta/F");
  muTree->Branch("GammaIso_MeanPhi",&fGammaIso_MeanPhi,"fGammaIso_MeanPhi/F");
  muTree->Branch("GammaIso_SigmaEtaEta",&fGammaIso_SigmaEtaEta,"GammaIso_SigmaEtaEta/F");
  muTree->Branch("GammaIso_SigmaPhiPhi",&fGammaIso_SigmaPhiPhi,"GammaIso_SigmaPhiPhi/F");
  muTree->Branch("GammaIso_SigmaEtaPhi",&fGammaIso_SigmaEtaPhi,"GammaIso_SigmaEtaPhi/F");
  muTree->Branch("NeutralHadronIso_DR0p0To0p1",&fNeutralHadronIso_DR0p0To0p1,"NeutralHadronIso_DR0p0To0p1/F");
  muTree->Branch("NeutralHadronIso_DR0p1To0p2",&fNeutralHadronIso_DR0p1To0p2,"NeutralHadronIso_DR0p1To0p2/F");
  muTree->Branch("NeutralHadronIso_DR0p2To0p3",&fNeutralHadronIso_DR0p2To0p3,"NeutralHadronIso_DR0p2To0p3/F");
  muTree->Branch("NeutralHadronIso_DR0p3To0p4",&fNeutralHadronIso_DR0p3To0p4,"NeutralHadronIso_DR0p3To0p4/F");
  muTree->Branch("NeutralHadronIso_DR0p4To0p5",&fNeutralHadronIso_DR0p4To0p5,"NeutralHadronIso_DR0p4To0p5/F");
  muTree->Branch("NeutralHadronIso_DR0p5To0p7",&fNeutralHadronIso_DR0p5To0p7,"NeutralHadronIso_DR0p5To0p7/F");
  muTree->Branch("NeutralHadronIso_DR0p7To1p0",&fNeutralHadronIso_DR0p7To1p0,"NeutralHadronIso_DR0p7To1p0/F");
  muTree->Branch("NeutralHadronIso_MeanEta",&fNeutralHadronIso_MeanEta,"NeutralHadronIso_MeanEta/F");
  muTree->Branch("NeutralHadronIso_MeanPhi",&fNeutralHadronIso_MeanPhi,"fNeutralHadronIso_MeanPhi/F");
  muTree->Branch("NeutralHadronIso_SigmaEtaEta",&fNeutralHadronIso_SigmaEtaEta,"NeutralHadronIso_SigmaEtaEta/F");
  muTree->Branch("NeutralHadronIso_SigmaPhiPhi",&fNeutralHadronIso_SigmaPhiPhi,"NeutralHadronIso_SigmaPhiPhi/F");
  muTree->Branch("NeutralHadronIso_SigmaEtaPhi",&fNeutralHadronIso_SigmaEtaPhi,"NeutralHadronIso_SigmaEtaPhi/F");
  muTree->Branch("DirectionalPFIso",&fDirectionalPFIso,"DirectionalPFIso/F");
  muTree->Branch("DirectionalChargedIso",&fDirectionalChargedIso,"DirectionalChargedIso/F");
  muTree->Branch("DirectionalGammaIso",&fDirectionalGammaIso,"DirectionalGammaIso/F");
  muTree->Branch("DirectionalNeutralHadronIso",&fDirectionalNeutralHadronIso,"DirectionalNeutralHadronIso/F");


  UInt_t NMuonsFilled = 0;
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  TClonesArray *photonArr = new TClonesArray("mithep::TPhoton");
  TClonesArray *pfcandidateArr = new TClonesArray("mithep::TPFCandidate");

  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("/data/smurf/data/Winter11_4700ipb/auxiliar/hww.Full2011.json"); 

  Int_t NEvents = 0;

  vector<string> inputfiles;
  if (inputFilename == "LIST") {
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-del-a05-v1.ZPlusFakeSkim.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-del-m10-v1.ZPlusFakeSkim.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-del-o03-v1.ZPlusFakeSkim.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-del-pr-v4.ZPlusFakeSkim.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11b-del-pr-v1.ZPlusFakeSkim.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-dmu-a05-v1.ZPlusFakeSkim.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-dmu-m10-v1.ZPlusFakeSkim.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-dmu-o03-v1.ZPlusFakeSkim.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-dmu-pr-v4.ZPlusFakeSkim.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11b-dmu-pr-v1.ZPlusFakeSkim.root");
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

      //Split into Training and Test Trees
      if (TreeType == 0 && info->evtNum % 2 != 0 ) continue;
      if (TreeType == 1 && info->evtNum % 2 == 0 ) continue;

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
      Double_t rho = 0;
      if (!(TMath::IsNaN(info->PileupEnergyDensity) || isinf(info->PileupEnergyDensity))) rho = info->PileupEnergyDensity;


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
      vector<const mithep::TMuon*> goodMuons;
      vector<const mithep::TElectron*> goodElectrons;
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
        if(passMuonID(mu)) goodMuons.push_back(mu);
      }
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
        if (ele->pt > 5) NElectrons++;
        if(passCutBasedEleID(ele)) goodElectrons.push_back(ele);
      }

      Int_t NZCandidates = 0;
      const mithep::TMuon* ZMuon1 = 0;
      const mithep::TMuon* ZMuon2 = 0;
      const mithep::TElectron *ZEle1 = 0;
      const mithep::TElectron *ZEle2 = 0;
      Double_t ZPt = 0;
      Double_t ZMass = 0;
      for(Int_t i=0; i<goodMuons.size(); i++) {
        mithep::FourVectorM mu1;
        mu1.SetCoordinates(goodMuons[i]->pt, goodMuons[i]->eta, goodMuons[i]->phi, 105.658369e-3 );
        for(Int_t j=i+1; j<goodMuons.size(); j++) {
          mithep::FourVectorM mu2;
          mu2.SetCoordinates(goodMuons[j]->pt, goodMuons[j]->eta, goodMuons[j]->phi, 105.658369e-3 );
          mithep::FourVectorM dilepton = mu1+mu2;
          if (dilepton.M() > 75 && dilepton.M() < 105) {
            ZMuon1 = goodMuons[i];
            ZMuon2 = goodMuons[j];
            ZPt = dilepton.Pt();
            ZMass = dilepton.M();
            NZCandidates++;
          }
        }
      }
      for(Int_t i=0; i<goodElectrons.size(); i++) {
        mithep::FourVectorM ele1;
        ele1.SetCoordinates(goodElectrons[i]->pt, goodElectrons[i]->eta, goodElectrons[i]->phi, 0.51099892e-3 );
        for(Int_t j=i+1; j<goodElectrons.size(); j++) {
          mithep::FourVectorM ele2;
          ele2.SetCoordinates(goodElectrons[j]->pt, goodElectrons[j]->eta, goodElectrons[j]->phi, 0.51099892e-3 );
          mithep::FourVectorM dilepton = ele1+ele2;
          if (dilepton.M() > 75 && dilepton.M() < 105) {
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

      //********************************************************
      // Find Electrons and Muons passing ID
      //********************************************************
      vector<Int_t> GenLeptonType;
      vector<Int_t> GenLeptonIndex;
      vector<Double_t> GenLeptonEta;
      vector<Double_t> GenLeptonPhi;
      vector<Int_t> IdentifiedLeptonType;
      vector<Int_t> IdentifiedLeptonIndex;

      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);      
        
        //match to gen muon
        if ((UInt_t(abs(max(0,mu->isMCReal))) & 2) == 2) {
          GenLeptonType.push_back(13);
          GenLeptonIndex.push_back(i);
          GenLeptonEta.push_back(mu->eta);
          GenLeptonPhi.push_back(mu->phi);
        }

        if ( mu->pt > 5.0
             &&
             fabs(mu->eta) < 2.4
             && 
             ( (mu->typeBits & kGlobal) == kGlobal || (mu->typeBits & kTracker) == kTracker)
             &&
             mu->nTkHits > 10
             &&
             passMuonIP_HZZ2011(mu)  //Only Include for Veto Option 1
          ) {
          IdentifiedLeptonType.push_back(13);
          IdentifiedLeptonIndex.push_back(i);
        }
      }
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
     
        //match to gen muon
        if ((UInt_t(abs(max(0,ele->isMCReal))) & 1) == 1) {
          GenLeptonType.push_back(11);
          GenLeptonIndex.push_back(i);
          GenLeptonEta.push_back(ele->eta);
          GenLeptonPhi.push_back(ele->phi);
        }

        if ( (0==0)
             &&
             ele->pt > 5.0
             && 
             fabs(ele->eta) < 2.5
             &&     
             PassCiCID(ele,1) //Veto Option1
//              &&     
//              passEleWP95ID(ele) //Veto Option2
          ) {
          IdentifiedLeptonType.push_back(11);
          IdentifiedLeptonIndex.push_back(i);
        }
      }



      //******************************************************************************
      //loop over muons
      //******************************************************************************
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);

        //make cut on dz
        if (fabs(mu->dz) > 0.2) continue;
        
        //veto against muons
        if (ZEle1 && mithep::MathUtils::DeltaR(mu->eta, mu->phi, ZEle1->eta, ZEle1->phi) < 0.3) continue;
        if (ZEle2 && mithep::MathUtils::DeltaR(mu->eta, mu->phi, ZEle2->eta, ZEle2->phi) < 0.3) continue;
        if (ZMuon1 && mithep::MathUtils::DeltaR(mu->eta, mu->phi, ZMuon1->eta, ZMuon1->phi) < 0.3) continue;
        if (ZMuon2 && mithep::MathUtils::DeltaR(mu->eta, mu->phi, ZMuon2->eta, ZMuon2->phi) < 0.3) continue;

        //pass denominator cuts
        if (!passMuonDenominatorCuts(mu)) continue;

        //Fill These Muons
        fWeight = 1.0;
        fRunNumber = info->runNum;
        fLumiSectionNumber = info->lumiSec;
        fEventNumber = info->evtNum;
        fMuPt = mu->pt; 
        fMuEta = mu->eta; 
        fMuPhi = mu->phi; 
        fMuPFIso = mu->ChargedIso03 + mu->NeutralIso03_05Threshold - info->PileupEnergyDensity*MuonEffectiveArea(kMuNeutralIso03,mu->eta);
   
        //CutBased Variables
        fMuTkNchi2 = mu->tkNchi2 ; 
        fMuGlobalNchi2 = mu->muNchi2 ; 
        fMuNValidHits = mu->nValidHits; 
        fMuNTrackerHits = mu->nTkHits; 
        fMuNPixelHits = mu->nPixHits; 
        fMuNMatches = mu->nMatch ; 
        fMuD0 = mu->d0 ; 

        //Additional Vars 
        fMuIP3d = mu->ip3d ; 
        fMuIP3dSig = mu->ip3dSig ; 
        fMuTrkKink = mu->TrkKink ; 
        fMuGlobalKink = mu->GlobalKink ; 
        fMuSegmentCompatibility = mu->SegmentCompatibility ; 
        fMuCaloCompatibility = mu->CaloCompatilibity ; 
        fMuHadEnergy = mu->HadEnergy; 
        fMuHoEnergy = mu->HoEnergy; 
        fMuEmEnergy = mu->EmEnergy; 
        fMuHadS9Energy = mu->HadS9Energy; 
        fMuHoS9Energy = mu->HoS9Energy; 
        fMuEmS9Energy = mu->EmS9Energy; 



        //Isolation Variables
        fMuChargedIso03 = mu->ChargedIso03 ; 
        fMuChargedIso03FromOtherVertices = mu->ChargedIso03FromOtherVertices ; 
        fMuNeutralIso03_05Threshold = mu->NeutralIso03_05Threshold ; 
        fMuNeutralIso03_10Threshold = mu->NeutralIso03_10Threshold ; 
        fMuChargedIso04 = mu->ChargedIso04 ; 
        fMuChargedIso04FromOtherVertices = mu->ChargedIso04FromOtherVertices ; 
        fMuNeutralIso04_05Threshold = mu->NeutralIso04_05Threshold ; 
        fMuNeutralIso04_10Threshold = mu->NeutralIso04_10Threshold ; 
        fMuTrkIso03 = mu->trkIso03; 
        fMuEMIso03 = mu->emIso03; 
        fMuHadIso03 = mu->hadIso03; 
        fMuTrkIso05 = mu->trkIso05; 
        fMuEMIso05 = mu->emIso05; 
        fMuHadIso05 = mu->hadIso05; 
        fRho = info->PileupEnergyDensity; 
        fNVertices = info->nPV0; 

        //*************************************************
        //Isolation STUFF
        //*************************************************
        Double_t tmpPFIso_NormalizedMeanEta  = 0;
        Double_t tmpPFIso_NormalizedMeanPhi  = 0;

        Double_t tmpPFIso_MeanEta_numerator  = 0;
        Double_t tmpPFIso_MeanPhi_numerator  = 0;
        Double_t tmpPFIso_SigmaEtaEta_numerator  = 0;
        Double_t tmpPFIso_SigmaPhiPhi_numerator  = 0;
        Double_t tmpPFIso_SigmaEtaPhi_numerator  = 0;
        Double_t tmpPFIso_Covariance_denominator  = 0;

        Double_t tmpChargedIso_DR0p0To0p05  = 0;
        Double_t tmpChargedIso_DR0p05To0p1  = 0;
        Double_t tmpChargedIso_DR0p1To0p15  = 0;
        Double_t tmpChargedIso_DR0p15To0p2  = 0;
        Double_t tmpChargedIso_DR0p2To0p25  = 0;
        Double_t tmpChargedIso_DR0p25To0p3  = 0;
        Double_t tmpChargedIso_DR0p3To0p4  = 0;
        Double_t tmpChargedIso_DR0p4To0p5  = 0;
        Double_t tmpChargedIso_DR0p5To0p7  = 0;
        Double_t tmpChargedIso_DR0p7To1p0  = 0;
        Double_t tmpChargedIso_MeanEta_numerator  = 0;
        Double_t tmpChargedIso_MeanPhi_numerator  = 0;
        Double_t tmpChargedIso_SigmaEtaEta_numerator  = 0;
        Double_t tmpChargedIso_SigmaPhiPhi_numerator  = 0;
        Double_t tmpChargedIso_SigmaEtaPhi_numerator  = 0;
        Double_t tmpChargedIso_Covariance_denominator  = 0;
        Double_t tmpGammaIso_DR0p0To0p05  = 0;
        Double_t tmpGammaIso_DR0p05To0p1  = 0;
        Double_t tmpGammaIso_DR0p1To0p15  = 0;
        Double_t tmpGammaIso_DR0p15To0p2  = 0;
        Double_t tmpGammaIso_DR0p2To0p25  = 0;
        Double_t tmpGammaIso_DR0p25To0p3  = 0;
        Double_t tmpGammaIso_DR0p3To0p4  = 0;
        Double_t tmpGammaIso_DR0p4To0p5  = 0;
        Double_t tmpGammaIso_DR0p5To0p7  = 0;
        Double_t tmpGammaIso_DR0p7To1p0  = 0;
        Double_t tmpGammaIso_MeanEta_numerator  = 0;
        Double_t tmpGammaIso_MeanPhi_numerator  = 0;
        Double_t tmpGammaIso_SigmaEtaEta_numerator  = 0;
        Double_t tmpGammaIso_SigmaPhiPhi_numerator  = 0;
        Double_t tmpGammaIso_SigmaEtaPhi_numerator  = 0;
        Double_t tmpGammaIso_Covariance_denominator  = 0;
        Double_t tmpNeutralHadronIso_DR0p0To0p05  = 0;
        Double_t tmpNeutralHadronIso_DR0p05To0p1  = 0;
        Double_t tmpNeutralHadronIso_DR0p1To0p15  = 0;
        Double_t tmpNeutralHadronIso_DR0p15To0p2  = 0;
        Double_t tmpNeutralHadronIso_DR0p2To0p25  = 0;
        Double_t tmpNeutralHadronIso_DR0p25To0p3  = 0;
        Double_t tmpNeutralHadronIso_DR0p3To0p4  = 0;
        Double_t tmpNeutralHadronIso_DR0p4To0p5  = 0;
        Double_t tmpNeutralHadronIso_DR0p5To0p7  = 0;
        Double_t tmpNeutralHadronIso_DR0p7To1p0  = 0;
        Double_t tmpNeutralHadronIso_MeanEta_numerator  = 0;
        Double_t tmpNeutralHadronIso_MeanPhi_numerator  = 0;
        Double_t tmpNeutralHadronIso_SigmaEtaEta_numerator  = 0;
        Double_t tmpNeutralHadronIso_SigmaPhiPhi_numerator  = 0;
        Double_t tmpNeutralHadronIso_SigmaEtaPhi_numerator  = 0;
        Double_t tmpNeutralHadronIso_Covariance_denominator  = 0;

 
        for(Int_t k=0; k<pfcandidateArr->GetEntries(); ++k) {
          const mithep::TPFCandidate *pf = (mithep::TPFCandidate*)((*pfcandidateArr)[k]);
          if (pf->matchedObjectType == 13 && pf->matchedObjectIndex == i) continue;
          Double_t deta = (mu->eta - pf->eta);
          Double_t dphi = mithep::MathUtils::DeltaPhi(Double_t(mu->phi),Double_t(pf->phi));
          Double_t dr = fabs(mithep::MathUtils::DeltaR(mu->phi,mu->eta, pf->phi, pf->eta));
          if (dr > 1.0) continue;

          Bool_t IsLeptonFootprint = kFALSE;
          //************************************************************
          // Lepton Footprint Removal
          //************************************************************
          if (dr < 1.0) {
            //Check for id'ed leptons
            for (Int_t q=0; q < IdentifiedLeptonType.size() ; ++q) {
              if (IdentifiedLeptonType[q] == 11) {
                const mithep::TElectron *tmpele = (mithep::TElectron*)((*electronArr)[IdentifiedLeptonIndex[q]]);
                if (IdentifiedLeptonType[q] == 11 && pf->matchedObjectType == 11 && pf->matchedObjectIndex == IdentifiedLeptonIndex[q]) {
                  IsLeptonFootprint = kTRUE;
                }
                if (pf->q != 0 && fabs(tmpele->scEta) > 1.479 && mithep::MathUtils::DeltaR(tmpele->phi,tmpele->eta, pf->phi, pf->eta) < 0.015) IsLeptonFootprint = kTRUE;
                if (pf->pfType == eGamma) {
                  if (fabs(tmpele->scEta) > 1.479) {
                    if (mithep::MathUtils::DeltaR(tmpele->phi,tmpele->eta, pf->phi, pf->eta) < 0.08) IsLeptonFootprint = kTRUE;
                  }
                }
              }
              if (IdentifiedLeptonType[q] == 13) {
                const mithep::TMuon *tmpmu = (mithep::TMuon*)((*muonArr)[IdentifiedLeptonIndex[q]]);
                if (IdentifiedLeptonType[q] == 13 && pf->matchedObjectType == 13 && pf->matchedObjectIndex == IdentifiedLeptonIndex[q]) {
                  IsLeptonFootprint = kTRUE;
                 } 
                if (pf->q != 0 && mithep::MathUtils::DeltaR(tmpmu->phi,tmpmu->eta, pf->phi, pf->eta) < 0.01) IsLeptonFootprint = kTRUE;
              }
            }
          }

          if (IsLeptonFootprint) {
            continue;
          }

          //charged
          if (!IsLeptonFootprint) {
            if (pf->q != 0) {
              if (abs(pf->dz) > 0.2) continue;
              //************************************************************
              // Veto any PFmuon, or PFEle
              if (pf->pfType == eElectron || pf->pfType == eMuon) continue;
              //************************************************************
              //************************************************************
              // Footprint Veto
              if (dr < 0.01) continue;
              //************************************************************
                
              if (dr < 0.05) tmpChargedIso_DR0p0To0p05 += pf->pt;
              if (dr >= 0.05 && dr < 0.10) tmpChargedIso_DR0p05To0p1 += pf->pt;
              if (dr >= 0.10 && dr < 0.15) tmpChargedIso_DR0p1To0p15 += pf->pt;
              if (dr >= 0.15 && dr < 0.20) tmpChargedIso_DR0p15To0p2 += pf->pt;
              if (dr >= 0.20 && dr < 0.25) tmpChargedIso_DR0p2To0p25 += pf->pt;
              if (dr >= 0.25 && dr < 0.3) tmpChargedIso_DR0p25To0p3 += pf->pt;
              if (dr >= 0.3 && dr < 0.4) tmpChargedIso_DR0p3To0p4 += pf->pt;
              if (dr >= 0.4 && dr < 0.5) tmpChargedIso_DR0p4To0p5 += pf->pt;
              if (dr >= 0.5 && dr < 0.7) tmpChargedIso_DR0p5To0p7 += pf->pt;
              if (dr >= 0.7 && dr < 1.0) tmpChargedIso_DR0p7To1p0 += pf->pt;
              tmpChargedIso_Covariance_denominator += pf->pt;
              tmpChargedIso_MeanEta_numerator += pf->pt * deta;
              tmpChargedIso_MeanPhi_numerator += pf->pt * dphi;
              tmpChargedIso_SigmaEtaEta_numerator += pf->pt * deta*deta;
              tmpChargedIso_SigmaPhiPhi_numerator += pf->pt * dphi*dphi;
              tmpChargedIso_SigmaEtaPhi_numerator += pf->pt * deta*dphi;
            }
            else if (pf->pfType == eGamma) {
              //************************************************************
              // Footprint Veto
              // if (dr < 0.05) continue;
              //************************************************************
              if (dr < 0.05) tmpGammaIso_DR0p0To0p05 += pf->pt;
              if (dr >= 0.05 && dr < 0.10) tmpGammaIso_DR0p05To0p1 += pf->pt;
              if (dr >= 0.10 && dr < 0.15) tmpGammaIso_DR0p1To0p15 += pf->pt;
              if (dr >= 0.15 && dr < 0.20) tmpGammaIso_DR0p15To0p2 += pf->pt;
              if (dr >= 0.20 && dr < 0.25) tmpGammaIso_DR0p2To0p25 += pf->pt;
              if (dr >= 0.25 && dr < 0.3) tmpGammaIso_DR0p25To0p3 += pf->pt;
              if (dr >= 0.3 && dr < 0.4) tmpGammaIso_DR0p3To0p4 += pf->pt;
              if (dr >= 0.4 && dr < 0.5) tmpGammaIso_DR0p4To0p5 += pf->pt;
              if (dr >= 0.5 && dr < 0.7) tmpGammaIso_DR0p5To0p7 += pf->pt;
              if (dr >= 0.7 && dr < 1.0) tmpGammaIso_DR0p7To1p0 += pf->pt;
              tmpGammaIso_Covariance_denominator += pf->pt;
              tmpGammaIso_MeanEta_numerator += pf->pt * deta;
              tmpGammaIso_MeanPhi_numerator += pf->pt * dphi;
              tmpGammaIso_SigmaEtaEta_numerator += pf->pt * deta*deta;
              tmpGammaIso_SigmaPhiPhi_numerator += pf->pt * dphi*dphi;
              tmpGammaIso_SigmaEtaPhi_numerator += pf->pt * deta*dphi;
            }
            else {
              if (dr < 0.05) tmpNeutralHadronIso_DR0p0To0p05 += pf->pt;
              if (dr >= 0.05 && dr < 0.10) tmpNeutralHadronIso_DR0p05To0p1 += pf->pt;
              if (dr >= 0.10 && dr < 0.15) tmpNeutralHadronIso_DR0p1To0p15 += pf->pt;
              if (dr >= 0.15 && dr < 0.20) tmpNeutralHadronIso_DR0p15To0p2 += pf->pt;
              if (dr >= 0.20 && dr < 0.25) tmpNeutralHadronIso_DR0p2To0p25 += pf->pt;
              if (dr >= 0.25 && dr < 0.3) tmpNeutralHadronIso_DR0p25To0p3 += pf->pt;
              if (dr >= 0.3 && dr < 0.4) tmpNeutralHadronIso_DR0p3To0p4 += pf->pt;
              if (dr >= 0.4 && dr < 0.5) tmpNeutralHadronIso_DR0p4To0p5 += pf->pt;
              if (dr >= 0.5 && dr < 0.7) tmpNeutralHadronIso_DR0p5To0p7 += pf->pt;
              if (dr >= 0.7 && dr < 1.0) tmpNeutralHadronIso_DR0p7To1p0 += pf->pt;
              tmpNeutralHadronIso_Covariance_denominator += pf->pt;
              tmpNeutralHadronIso_MeanEta_numerator += pf->pt * deta;
              tmpNeutralHadronIso_MeanPhi_numerator += pf->pt * dphi;
              tmpNeutralHadronIso_SigmaEtaEta_numerator += pf->pt * deta*deta;
              tmpNeutralHadronIso_SigmaPhiPhi_numerator += pf->pt * dphi*dphi;
              tmpNeutralHadronIso_SigmaEtaPhi_numerator += pf->pt * deta*dphi;
            }

            //For combined PFIso Shape
            tmpPFIso_Covariance_denominator += pf->pt;
            tmpPFIso_MeanEta_numerator += pf->pt * deta;
            tmpPFIso_MeanPhi_numerator += pf->pt * dphi;
            tmpPFIso_SigmaEtaEta_numerator += pf->pt * deta*deta;
            tmpPFIso_SigmaPhiPhi_numerator += pf->pt * dphi*dphi;
            tmpPFIso_SigmaEtaPhi_numerator += pf->pt * deta*dphi;

            //For directional iso calculation  
            if (dr > 0) {
              tmpPFIso_NormalizedMeanEta += pf->pt * deta / dr;
              tmpPFIso_NormalizedMeanPhi += pf->pt * dphi / dr;
            }
          }

        }

        //***********************
        //Fill isolation rings
        //***********************
        fChargedIso_DR0p0To0p1   = TMath::Min((tmpChargedIso_DR0p0To0p05 + tmpChargedIso_DR0p05To0p1)/mu->pt, 2.5);
        fChargedIso_DR0p1To0p2   = TMath::Min((tmpChargedIso_DR0p1To0p15 + tmpChargedIso_DR0p15To0p2)/mu->pt, 2.5);
        fChargedIso_DR0p2To0p3 = TMath::Min((tmpChargedIso_DR0p2To0p25 + tmpChargedIso_DR0p25To0p3)/mu->pt, 2.5);
        fChargedIso_DR0p3To0p4 = TMath::Min((tmpChargedIso_DR0p3To0p4)/mu->pt, 2.5);
        fChargedIso_DR0p4To0p5 = TMath::Min((tmpChargedIso_DR0p4To0p5)/mu->pt, 2.5);
        fChargedIso_DR0p5To0p7 = TMath::Min((tmpChargedIso_DR0p5To0p7)/mu->pt, 2.5);
        fChargedIso_DR0p7To1p0 = TMath::Min((tmpChargedIso_DR0p7To1p0)/mu->pt, 2.5);
        fGammaIso_DR0p0To0p1   = TMath::Max(TMath::Min((tmpGammaIso_DR0p0To0p05 + tmpGammaIso_DR0p05To0p1 + tmpGammaIso_DR0p1To0p15 - rho*MuonEffectiveArea(kMuGammaIsoDR0p0To0p1, mu->eta))/mu->pt, 2.5), 0.0);
        fGammaIso_DR0p1To0p2   = TMath::Max(TMath::Min((tmpGammaIso_DR0p1To0p15 + tmpGammaIso_DR0p15To0p2 - rho*MuonEffectiveArea(kMuGammaIsoDR0p1To0p2, mu->eta))/mu->pt, 2.5), 0.0);
        fGammaIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpGammaIso_DR0p2To0p25 + tmpGammaIso_DR0p25To0p3 - rho*MuonEffectiveArea(kMuGammaIsoDR0p2To0p3, mu->eta))/mu->pt, 2.5), 0.0);
        fGammaIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpGammaIso_DR0p3To0p4 - rho*MuonEffectiveArea(kMuGammaIsoDR0p3To0p4, mu->eta))/mu->pt, 2.5), 0.0);
        fGammaIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpGammaIso_DR0p4To0p5 - rho*MuonEffectiveArea(kMuGammaIsoDR0p4To0p5, mu->eta))/mu->pt, 2.5), 0.0);
        fGammaIso_DR0p5To0p7 = TMath::Max(TMath::Min((tmpGammaIso_DR0p5To0p7 - rho*MuonEffectiveArea(kMuGammaIsoDR0p5To0p7, mu->eta))/mu->pt, 2.5), 0.0);
        fGammaIso_DR0p7To1p0 = TMath::Max(TMath::Min((tmpGammaIso_DR0p7To1p0 - rho*MuonEffectiveArea(kMuGammaIsoDR0p7To1p0, mu->eta))/mu->pt, 2.5), 0.0);
        fNeutralHadronIso_DR0p0To0p1   = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p0To0p05 + tmpNeutralHadronIso_DR0p05To0p1 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p0To0p1, mu->eta))/mu->pt, 2.5), 0.0);
        fNeutralHadronIso_DR0p1To0p2   = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p1To0p15 + tmpNeutralHadronIso_DR0p15To0p2 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p1To0p2, mu->eta))/mu->pt, 2.5), 0.0);
        fNeutralHadronIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p2To0p25 + tmpNeutralHadronIso_DR0p25To0p3 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p2To0p3, mu->eta))/mu->pt, 2.5), 0.0);
        fNeutralHadronIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p3To0p4 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p3To0p4, mu->eta))/mu->pt, 2.5), 0.0);
        fNeutralHadronIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p4To0p5 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p4To0p5, mu->eta))/mu->pt, 2.5), 0.0);
        fNeutralHadronIso_DR0p5To0p7 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p5To0p7 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p5To0p7, mu->eta))/mu->pt, 2.5), 0.0);
        fNeutralHadronIso_DR0p7To1p0 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p7To1p0 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p7To1p0, mu->eta))/mu->pt, 2.5), 0.0);
        
        Double_t tmpPFIso_MeanEta  = (tmpPFIso_Covariance_denominator>0) ? tmpPFIso_MeanEta_numerator/tmpPFIso_Covariance_denominator : 0;
        Double_t tmpPFIso_MeanPhi  = (tmpPFIso_Covariance_denominator>0) ? tmpPFIso_MeanPhi_numerator/tmpPFIso_Covariance_denominator : 0;
        Double_t tmpPFIso_SigmaEtaEta  = (tmpPFIso_Covariance_denominator>0) ? sqrt(tmpPFIso_SigmaEtaEta_numerator/tmpPFIso_Covariance_denominator) : 0;
        Double_t tmpPFIso_SigmaPhiPhi  = (tmpPFIso_Covariance_denominator>0) ? sqrt(tmpPFIso_SigmaPhiPhi_numerator/tmpPFIso_Covariance_denominator) : 0;
        Double_t tmpPFIso_CovEtaPhi  = (tmpPFIso_Covariance_denominator>0) ? (tmpPFIso_SigmaEtaPhi_numerator/tmpPFIso_Covariance_denominator) : 0;
        Double_t tmpPFIso_SigmaEtaPhi = 0;
        if (tmpPFIso_SigmaEtaEta*tmpPFIso_SigmaPhiPhi != 0) tmpPFIso_SigmaEtaPhi = tmpPFIso_CovEtaPhi / (tmpPFIso_SigmaEtaEta*tmpPFIso_SigmaPhiPhi);
        else tmpPFIso_SigmaEtaPhi = (tmpPFIso_CovEtaPhi < 0) ? -1 : (tmpPFIso_CovEtaPhi > 0);

        Double_t tmpChargedIso_MeanEta  = (tmpChargedIso_Covariance_denominator>0) ? tmpChargedIso_MeanEta_numerator/tmpChargedIso_Covariance_denominator : 0;
        Double_t tmpChargedIso_MeanPhi  = (tmpChargedIso_Covariance_denominator>0) ? tmpChargedIso_MeanPhi_numerator/tmpChargedIso_Covariance_denominator : 0;
        Double_t tmpChargedIso_SigmaEtaEta  = (tmpChargedIso_Covariance_denominator>0) ? sqrt(tmpChargedIso_SigmaEtaEta_numerator/tmpChargedIso_Covariance_denominator) : 0;
        Double_t tmpChargedIso_SigmaPhiPhi  = (tmpChargedIso_Covariance_denominator>0) ? sqrt(tmpChargedIso_SigmaPhiPhi_numerator/tmpChargedIso_Covariance_denominator) : 0;
        Double_t tmpChargedIso_CovEtaPhi  = (tmpChargedIso_Covariance_denominator>0) ? (tmpChargedIso_SigmaEtaPhi_numerator/tmpChargedIso_Covariance_denominator) : 0;
        Double_t tmpChargedIso_SigmaEtaPhi = 0;
        if (tmpChargedIso_SigmaEtaEta*tmpChargedIso_SigmaPhiPhi != 0) tmpChargedIso_SigmaEtaPhi = tmpChargedIso_CovEtaPhi / (tmpChargedIso_SigmaEtaEta*tmpChargedIso_SigmaPhiPhi);
        else tmpChargedIso_SigmaEtaPhi = (tmpChargedIso_CovEtaPhi < 0) ? -1 : (tmpChargedIso_CovEtaPhi > 0);

        Double_t tmpGammaIso_MeanEta  = (tmpGammaIso_Covariance_denominator>0) ? tmpGammaIso_MeanEta_numerator/tmpGammaIso_Covariance_denominator : 0;
        Double_t tmpGammaIso_MeanPhi  = (tmpGammaIso_Covariance_denominator>0) ? tmpGammaIso_MeanPhi_numerator/tmpGammaIso_Covariance_denominator : 0;
        Double_t tmpGammaIso_SigmaEtaEta  = (tmpGammaIso_Covariance_denominator>0) ? sqrt(tmpGammaIso_SigmaEtaEta_numerator/tmpGammaIso_Covariance_denominator) : 0;
        Double_t tmpGammaIso_SigmaPhiPhi  = (tmpGammaIso_Covariance_denominator>0) ? sqrt(tmpGammaIso_SigmaPhiPhi_numerator/tmpGammaIso_Covariance_denominator) : 0;
        Double_t tmpGammaIso_CovEtaPhi  = (tmpGammaIso_Covariance_denominator>0) ? (tmpGammaIso_SigmaEtaPhi_numerator/tmpGammaIso_Covariance_denominator) : 0;
        Double_t tmpGammaIso_SigmaEtaPhi = 0;
        if (tmpGammaIso_SigmaEtaEta*tmpGammaIso_SigmaPhiPhi != 0) tmpGammaIso_SigmaEtaPhi = tmpGammaIso_CovEtaPhi / (tmpGammaIso_SigmaEtaEta*tmpGammaIso_SigmaPhiPhi);
        else tmpGammaIso_SigmaEtaPhi = (tmpGammaIso_CovEtaPhi < 0) ? -1 : (tmpGammaIso_CovEtaPhi > 0);

        Double_t tmpNeutralHadronIso_MeanEta  = (tmpNeutralHadronIso_Covariance_denominator>0) ? tmpNeutralHadronIso_MeanEta_numerator/tmpNeutralHadronIso_Covariance_denominator : 0;
        Double_t tmpNeutralHadronIso_MeanPhi  = (tmpNeutralHadronIso_Covariance_denominator>0) ? tmpNeutralHadronIso_MeanPhi_numerator/tmpNeutralHadronIso_Covariance_denominator : 0;
        Double_t tmpNeutralHadronIso_SigmaEtaEta  = (tmpNeutralHadronIso_Covariance_denominator>0) ? sqrt(tmpNeutralHadronIso_SigmaEtaEta_numerator/tmpNeutralHadronIso_Covariance_denominator) : 0;
        Double_t tmpNeutralHadronIso_SigmaPhiPhi  = (tmpNeutralHadronIso_Covariance_denominator>0) ? sqrt(tmpNeutralHadronIso_SigmaPhiPhi_numerator/tmpNeutralHadronIso_Covariance_denominator) : 0;
        Double_t tmpNeutralHadronIso_CovEtaPhi  = (tmpNeutralHadronIso_Covariance_denominator>0) ? (tmpNeutralHadronIso_SigmaEtaPhi_numerator/tmpNeutralHadronIso_Covariance_denominator) : 0;
        Double_t tmpNeutralHadronIso_SigmaEtaPhi = 0;
        if (tmpNeutralHadronIso_SigmaEtaEta*tmpNeutralHadronIso_SigmaPhiPhi != 0) tmpNeutralHadronIso_SigmaEtaPhi = tmpNeutralHadronIso_CovEtaPhi / (tmpNeutralHadronIso_SigmaEtaEta*tmpNeutralHadronIso_SigmaPhiPhi);
        else tmpNeutralHadronIso_SigmaEtaPhi = (tmpNeutralHadronIso_CovEtaPhi < 0) ? -1 : (tmpNeutralHadronIso_CovEtaPhi > 0);

        //***********************
        //Fill Parton Shower Shapes
        //***********************
        fPFIso_MeanEta = tmpPFIso_MeanEta;
        fPFIso_MeanPhi = tmpPFIso_MeanPhi;
        fPFIso_SigmaEtaEta = tmpPFIso_SigmaEtaEta;
        fPFIso_SigmaPhiPhi = tmpPFIso_SigmaPhiPhi;
        fPFIso_SigmaEtaPhi = tmpPFIso_SigmaEtaPhi;
        fChargedIso_MeanEta = tmpChargedIso_MeanEta;
        fChargedIso_MeanPhi = tmpChargedIso_MeanPhi;
        fChargedIso_SigmaEtaEta = tmpChargedIso_SigmaEtaEta;
        fChargedIso_SigmaPhiPhi = tmpChargedIso_SigmaPhiPhi;
        fChargedIso_SigmaEtaPhi = tmpChargedIso_SigmaEtaPhi;
        fGammaIso_MeanEta = tmpGammaIso_MeanEta;
        fGammaIso_MeanPhi = tmpGammaIso_MeanPhi;
        fGammaIso_SigmaEtaEta = tmpGammaIso_SigmaEtaEta;
        fGammaIso_SigmaPhiPhi = tmpGammaIso_SigmaPhiPhi;
        fGammaIso_SigmaEtaPhi = tmpGammaIso_SigmaEtaPhi;
        fNeutralHadronIso_MeanEta = tmpNeutralHadronIso_MeanEta;
        fNeutralHadronIso_MeanPhi = tmpNeutralHadronIso_MeanPhi;
        fNeutralHadronIso_SigmaEtaEta = tmpNeutralHadronIso_SigmaEtaEta;
        fNeutralHadronIso_SigmaPhiPhi = tmpNeutralHadronIso_SigmaPhiPhi;
        fNeutralHadronIso_SigmaEtaPhi = tmpNeutralHadronIso_SigmaEtaPhi;
               

        //Compute directional isolation
        TVector2 PFIsoCentroid(tmpPFIso_NormalizedMeanEta,tmpPFIso_NormalizedMeanPhi);
        Double_t tmpPFIsoAngleSqrSum = 0;
        Double_t tmpChargedIsoAngleSqrSum = 0;
        Double_t tmpGammaIsoAngleSqrSum = 0;
        Double_t tmpNeutralHadronIsoAngleSqrSum = 0;
        for(Int_t k=0; k<pfcandidateArr->GetEntries(); ++k) {
          const mithep::TPFCandidate *pf = (mithep::TPFCandidate*)((*pfcandidateArr)[k]);
          if (pf->matchedObjectType == 13 && pf->matchedObjectIndex == i) continue;
          Double_t deta = (mu->eta - pf->eta);
          Double_t dphi = mithep::MathUtils::DeltaPhi(Double_t(mu->phi),Double_t(pf->phi));
          Double_t dr = fabs(mithep::MathUtils::DeltaR(mu->phi,mu->eta, pf->phi, pf->eta));
          if (dr > 1.0) continue;
                            
          TVector2 *tmpPFAngle;
          if (dr > 0) {
            tmpPFAngle = new TVector2(pf->pt*deta/dr, pf->pt*dphi/dr);
          } else {
            continue;
          }

          //if one of the vectors is 0, then assign 0 to the angle -> don't do anything
          if (tmpPFAngle->Mod()*PFIsoCentroid.Mod() == 0) continue;
          //if the angle is 0, don't do anything. to avoid acos(1) = nan;
          if (fabs((tmpPFAngle->Px()*PFIsoCentroid.Px() + tmpPFAngle->Py()*PFIsoCentroid.Py())/(tmpPFAngle->Mod()*PFIsoCentroid.Mod()) - 1) < 0.001) continue;
          
          if (pf->q != 0) {
            if (abs(pf->dz) > 0.2) continue;
            if (pf->pfType == eElectron || pf->pfType == eMuon) continue;
            if (dr < 0.010) continue;
            tmpChargedIsoAngleSqrSum += pow(acos((tmpPFAngle->Px()*PFIsoCentroid.Px() + tmpPFAngle->Py()*PFIsoCentroid.Py())/(tmpPFAngle->Mod()*PFIsoCentroid.Mod())),2);                
          }
          else if (pf->pfType == eGamma) {
            tmpGammaIsoAngleSqrSum += pow(acos((tmpPFAngle->Px()*PFIsoCentroid.Px() + tmpPFAngle->Py()*PFIsoCentroid.Py())/(tmpPFAngle->Mod()*PFIsoCentroid.Mod())),2);                
          }
          else {
            tmpNeutralHadronIsoAngleSqrSum += pow(acos((tmpPFAngle->Px()*PFIsoCentroid.Px() + tmpPFAngle->Py()*PFIsoCentroid.Py())/(tmpPFAngle->Mod()*PFIsoCentroid.Mod())),2);                
          }
          tmpPFIsoAngleSqrSum += pow(acos((tmpPFAngle->Px()*PFIsoCentroid.Px() + tmpPFAngle->Py()*PFIsoCentroid.Py())/(tmpPFAngle->Mod()*PFIsoCentroid.Mod())),2);                              
          delete tmpPFAngle;
        }

        //***********************
        //fill directional iso hists
        //***********************
        fDirectionalPFIso = TMath::Min(tmpPFIsoAngleSqrSum,300.0);
        fDirectionalChargedIso = TMath::Min(tmpChargedIsoAngleSqrSum,200.0);
        fDirectionalGammaIso = TMath::Min(tmpGammaIsoAngleSqrSum,200.0);
        fDirectionalNeutralHadronIso = TMath::Min(tmpNeutralHadronIsoAngleSqrSum,200.0);

        NMuonsFilled++;
        muTree->Fill();

      } //loop over muons

    }

    cout << "Total Muons: " << NMuonsFilled << endl;

  } //end loop over files

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;


  cout << "Total Muons: " << NMuonsFilled << endl;
  outputFile->Write();
  outputFile->Close();

  gBenchmark->Show("WWTemplate");       
} 


//--------------------------------------------------------------------------------------------------
Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(mu->pt > 5 && fabs(mu->eta) < 2.4)) pass = kFALSE;

  if (! 
      ( ((mu->typeBits & kGlobal == kGlobal) || (mu->typeBits & kTracker == kTracker))
        && fabs(mu->dz) < 0.2
        && (mu->ChargedIso03) / mu->pt < 0.7
        )
    )
    pass = kFALSE;    

  return pass;
}


