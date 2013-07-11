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
#include <MitStyle.h>

// define structures to read in ntuple
#include "HiggsAna/Ntupler/interface/HiggsAnaDefs.hh"
#include "HiggsAna/Ntupler/interface/TEventInfo.hh"
#include "HiggsAna/Ntupler/interface/TElectron.hh"
#include "HiggsAna/Ntupler/interface/TPhoton.hh"
#include "HiggsAna/Ntupler/interface/TMuon.hh"
#include "HiggsAna/Ntupler/interface/TJet.hh"
#include "HiggsAna/Ntupler/interface/TPFCandidate.hh"
#include "HiggsAna/Utils/LeptonIDCuts.hh"
#include "HiggsAna/HZZ4l/Utils/LeptonSelection.hh"


// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include "MitHiggs/Utils/interface/EfficiencyUtils.h"
#include "MitHiggs/Utils/interface/PlotUtils.h"
#include "MitPhysics/Utils/interface/MuonIDMVA.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

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
Bool_t passMuonDenominatorCuts(const mithep::TMuon *ele, Int_t DenominatorType);

void DoComputeMuonFakeRate(const string inputFilename,
                               const string label, 
                               const string outputFilename,  
                               const string minioutputFilename, Int_t Option = 0);

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
void ComputeMuonFakeRate_Data_ZSample(Int_t Option=0) {


//************************************
//MVA 
//************************************

//------------------
//Fully Inclusive
//------------------
  if (Option==0) {
     DoComputeMuonFakeRate("LIST","MuonFakeRate","MuonFakeRate.BaselineDetIso.ZSample.root", "FakeRates_Muon_BaselineDetIso_ZSample", 1);
  } 
  if (Option==1) {
     DoComputeMuonFakeRate("LIST","MuonFakeRate","MuonFakeRate.HZZIDMVAAndIsoMVALooseWP.ZSample.root", "FakeRates_Muon_HZZIDMVAAndIsoMVALooseWP_ZSample", 10);     
  } 
  if (Option==2) {
     DoComputeMuonFakeRate("LIST","MuonFakeRate","MuonFakeRate.HZZIDMVAAndIsoMVATightWP.ZSample.root", "FakeRates_Muon_HZZIDMVAAndIsoMVATightWP_ZSample", 11);     
  } 

}



void DoComputeMuonFakeRate(const string inputFilename,
                               const string label, 
                               const string outputFilename, 
                               const string smurfOutputFilename,
                               Int_t Option)
{  
  gBenchmark->Start("WWTemplate");

  Double_t LUMI = 26.5;

  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  mithep::MuonIDMVA *muonIDMVA = new mithep::MuonIDMVA();
  vector<string> muonIDMVAWeightFiles;
  muonIDMVAWeightFiles.push_back("MitPhysics/data/MuonMVAWeights/MuonIDMVA_BDTG_V0_barrel_lowpt.weights.xml");
  muonIDMVAWeightFiles.push_back("MitPhysics/data/MuonMVAWeights/MuonIDMVA_BDTG_V0_barrel_highpt.weights.xml");
  muonIDMVAWeightFiles.push_back("MitPhysics/data/MuonMVAWeights/MuonIDMVA_BDTG_V0_endcap_lowpt.weights.xml");
  muonIDMVAWeightFiles.push_back("MitPhysics/data/MuonMVAWeights/MuonIDMVA_BDTG_V0_endcap_highpt.weights.xml");
  muonIDMVAWeightFiles.push_back("MitPhysics/data/MuonMVAWeights/MuonIDMVA_BDTG_V0_tracker.weights.xml");
  muonIDMVAWeightFiles.push_back("MitPhysics/data/MuonMVAWeights/MuonIDMVA_BDTG_V0_global.weights.xml");
  muonIDMVA->Initialize( "MuonIDMVA", mithep::MuonIDMVA::kIDV0,
                         kTRUE, muonIDMVAWeightFiles);

  mithep::MuonIDMVA *muonIsoMVA = new mithep::MuonIDMVA();
  vector<string> muonIsoMVAWeightFiles;
  muonIsoMVAWeightFiles.push_back("MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_barrel_lowpt.weights.xml");
  muonIsoMVAWeightFiles.push_back("MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_barrel_highpt.weights.xml");
  muonIsoMVAWeightFiles.push_back("MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_endcap_lowpt.weights.xml");
  muonIsoMVAWeightFiles.push_back("MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_endcap_highpt.weights.xml");
  muonIsoMVAWeightFiles.push_back("MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_tracker.weights.xml");
  muonIsoMVAWeightFiles.push_back("MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_global.weights.xml");
  muonIsoMVA->Initialize( "MuonIsoMVA", mithep::MuonIDMVA::kIsoRingsV0,
                          kTRUE, muonIsoMVAWeightFiles);

  Float_t                 varChargedIso_DR0p0To0p1;
    Float_t                 varChargedIso_DR0p1To0p2;
    Float_t                 varChargedIso_DR0p2To0p3;
    Float_t                 varChargedIso_DR0p3To0p4;
    Float_t                 varChargedIso_DR0p4To0p5;
    Float_t                 varGammaIso_DR0p0To0p1;
    Float_t                 varGammaIso_DR0p1To0p2;
    Float_t                 varGammaIso_DR0p2To0p3;
    Float_t                 varGammaIso_DR0p3To0p4;
    Float_t                 varGammaIso_DR0p4To0p5;
    Float_t                 varNeutralHadronIso_DR0p0To0p1;
    Float_t                 varNeutralHadronIso_DR0p1To0p2;
    Float_t                 varNeutralHadronIso_DR0p2To0p3;
    Float_t                 varNeutralHadronIso_DR0p3To0p4;
    Float_t                 varNeutralHadronIso_DR0p4To0p5;
    Float_t                 varMuTkNchi2;
    Float_t                 varMuGlobalNchi2;
    Float_t                 varMuNValidHits;
    Float_t                 varMuNTrackerHits;
    Float_t                 varMuNPixelHits;
    Float_t                 varMuNMatches;
    Float_t                 varMuTrkKink;
    Float_t                 varMuSegmentCompatibility;
    Float_t                 varMuCaloCompatibility;
    Float_t                 varMuHadEnergy;
    Float_t                 varMuEmEnergy;
    Float_t                 varMuHadS9Energy;
    Float_t                 varMuEmS9Energy;

    Float_t varMuPt, varMuEta;

  TMVA::Reader  *MuonIDReader[6];
  for(UInt_t j=0; j<6; ++j) {
    MuonIDReader[j] = new TMVA::Reader( "!Color:!Silent:Error" );  
  }
  TMVA::Reader  *MuonIsoReader[6];
  for(UInt_t j=0; j<6; ++j) {
    MuonIsoReader[j] = new TMVA::Reader( "!Color:!Silent:Error" );  
  }

  for(UInt_t j=0; j<6; ++j) {
    MuonIsoReader[j]->AddVariable( "ChargedIso_DR0p0To0p1",         &varChargedIso_DR0p0To0p1           );
    MuonIsoReader[j]->AddVariable( "ChargedIso_DR0p1To0p2",         &varChargedIso_DR0p1To0p2           );
    MuonIsoReader[j]->AddVariable( "ChargedIso_DR0p2To0p3",       &varChargedIso_DR0p2To0p3         );
    MuonIsoReader[j]->AddVariable( "ChargedIso_DR0p3To0p4",        &varChargedIso_DR0p3To0p4          );
    MuonIsoReader[j]->AddVariable( "ChargedIso_DR0p4To0p5",        &varChargedIso_DR0p4To0p5          );
    MuonIsoReader[j]->AddVariable( "GammaIso_DR0p0To0p1",         &varGammaIso_DR0p0To0p1           );
    MuonIsoReader[j]->AddVariable( "GammaIso_DR0p1To0p2",         &varGammaIso_DR0p1To0p2           );
    MuonIsoReader[j]->AddVariable( "GammaIso_DR0p2To0p3",       &varGammaIso_DR0p2To0p3         );
    MuonIsoReader[j]->AddVariable( "GammaIso_DR0p3To0p4",        &varGammaIso_DR0p3To0p4          );
    MuonIsoReader[j]->AddVariable( "GammaIso_DR0p4To0p5",          &varGammaIso_DR0p4To0p5          );
    MuonIsoReader[j]->AddVariable( "NeutralHadronIso_DR0p0To0p1",         &varNeutralHadronIso_DR0p0To0p1           );
    MuonIsoReader[j]->AddVariable( "NeutralHadronIso_DR0p1To0p2",         &varNeutralHadronIso_DR0p1To0p2           );
    MuonIsoReader[j]->AddVariable( "NeutralHadronIso_DR0p2To0p3",       &varNeutralHadronIso_DR0p2To0p3         );
    MuonIsoReader[j]->AddVariable( "NeutralHadronIso_DR0p3To0p4",        &varNeutralHadronIso_DR0p3To0p4          );
    MuonIsoReader[j]->AddVariable( "NeutralHadronIso_DR0p4To0p5",        &varNeutralHadronIso_DR0p4To0p5          );
  }
  for(UInt_t j=0; j<6; ++j) {
    if (j==0) MuonIsoReader[j]->BookMVA( "MuonIsoMVA_V0_BDTG method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/MuonIsoMVA/MuonIsoMVA_BDTG_barrel_lowpt_V0.weights.xml");  
    if (j==1) MuonIsoReader[j]->BookMVA( "MuonIsoMVA_V0_BDTG method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/MuonIsoMVA/MuonIsoMVA_BDTG_barrel_highpt_V0.weights.xml"); 
    if (j==2) MuonIsoReader[j]->BookMVA( "MuonIsoMVA_V0_BDTG method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/MuonIsoMVA/MuonIsoMVA_BDTG_endcap_lowpt_V0.weights.xml");  
    if (j==3) MuonIsoReader[j]->BookMVA( "MuonIsoMVA_V0_BDTG method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/MuonIsoMVA/MuonIsoMVA_BDTG_endcap_highpt_V0.weights.xml"); 
    if (j==4) MuonIsoReader[j]->BookMVA( "MuonIsoMVA_V0_BDTG method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/MuonIsoMVA/MuonIsoMVA_BDTG_tracker_V0.weights.xml");  
    if (j==5) MuonIsoReader[j]->BookMVA( "MuonIsoMVA_V0_BDTG method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/MuonIsoMVA/MuonIsoMVA_BDTG_global_V0.weights.xml"); 
  }

  for(UInt_t j=0; j<6; ++j) {
    MuonIDReader[j]->AddVariable( "TkNchi2",               &varMuTkNchi2      );
    if(j!=4) MuonIDReader[j]->AddVariable( "GlobalNchi2",  &varMuGlobalNchi2  );
    if(j!=4) MuonIDReader[j]->AddVariable( "NValidHits",   &varMuNValidHits   );
    MuonIDReader[j]->AddVariable( "NTrackerHits",          &varMuNTrackerHits );
    MuonIDReader[j]->AddVariable( "NPixelHits",            &varMuNPixelHits   );
    if(j!=5) MuonIDReader[j]->AddVariable( "NMatches",     &varMuNMatches     );
    MuonIDReader[j]->AddVariable( "TrkKink",               &varMuTrkKink      );
    MuonIDReader[j]->AddVariable( "SegmentCompatibility",  &varMuSegmentCompatibility );
    MuonIDReader[j]->AddVariable( "CaloCompatibility",     &varMuCaloCompatibility    );
    MuonIDReader[j]->AddVariable( "HadEnergy",             &varMuHadEnergy    );
    MuonIDReader[j]->AddVariable( "EmEnergy",              &varMuEmEnergy     );
    MuonIDReader[j]->AddVariable( "HadS9Energy",           &varMuHadS9Energy  );
    MuonIDReader[j]->AddVariable( "EmS9Energy",            &varMuEmS9Energy   );
  }

  for(UInt_t j=0; j<6; ++j) {
    if (j==0) MuonIDReader[j]->BookMVA( "MuonIDMVA_V0_BDTG method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/MuonIDMVA/MuonIDMVA_BDTG_barrel_lowpt_V2.weights.xml");                
    if (j==1) MuonIDReader[j]->BookMVA( "MuonIDMVA_V0_BDTG method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/MuonIDMVA/MuonIDMVA_BDTG_barrel_highpt_V2.weights.xml");                
    if (j==2) MuonIDReader[j]->BookMVA( "MuonIDMVA_V0_BDTG method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/MuonIDMVA/MuonIDMVA_BDTG_endcap_lowpt_V2.weights.xml");                
    if (j==3) MuonIDReader[j]->BookMVA( "MuonIDMVA_V0_BDTG method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/MuonIDMVA/MuonIDMVA_BDTG_endcap_highpt_V2.weights.xml");                
    if (j==4) MuonIDReader[j]->BookMVA( "MuonIDMVA_V0_BDTG method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/MuonIDMVA/MuonIDMVA_BDTG_tracker_V2.weights.xml");                
    if (j==5) MuonIDReader[j]->BookMVA( "MuonIDMVA_V0_BDTG method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/MuonIDMVA/MuonIDMVA_BDTG_global_V2.weights.xml");                
  }

  //*****************************************************************************************
  //Define Pt bins
  //*****************************************************************************************
  vector<double> ptbins;
//   ptbins.push_back(10);  
//   ptbins.push_back(80);  


  ptbins.push_back(0);  
  ptbins.push_back(5);  
  ptbins.push_back(7);  
  ptbins.push_back(10);  
  ptbins.push_back(15);  
  ptbins.push_back(20);  
  ptbins.push_back(25);  
  ptbins.push_back(30);  
  ptbins.push_back(35);  
  ptbins.push_back(100);  


  vector<double> etabins;
  etabins.push_back(0.0);
  etabins.push_back(0.25);
  etabins.push_back(0.5);
  etabins.push_back(0.75);
  etabins.push_back(1.0);
  etabins.push_back(1.25);
  etabins.push_back(1.479);
  etabins.push_back(1.75);
  etabins.push_back(2.0);
  etabins.push_back(2.25);
  etabins.push_back(2.5);
  etabins.push_back(3.0);


  vector<double> phibins;
  phibins.push_back(-3.25);
  phibins.push_back(-2.75);
  phibins.push_back(-2.25);
  phibins.push_back(-1.75);
  phibins.push_back(-1.25);
  phibins.push_back(-0.75);
  phibins.push_back(-0.25);
  phibins.push_back(0.25);
  phibins.push_back(0.75);
  phibins.push_back(1.25);
  phibins.push_back(1.75);
  phibins.push_back(2.25);
  phibins.push_back(2.75);
  phibins.push_back(3.25);

  vector<double> nvtxbins;
  nvtxbins.push_back(0);
  nvtxbins.push_back(2);
  nvtxbins.push_back(4);
  nvtxbins.push_back(6);
  nvtxbins.push_back(8);
  nvtxbins.push_back(10);
  nvtxbins.push_back(12);
  nvtxbins.push_back(14);
  nvtxbins.push_back(16);
  nvtxbins.push_back(18);
  nvtxbins.push_back(20);
  nvtxbins.push_back(22);
  nvtxbins.push_back(24);
  nvtxbins.push_back(26);



  vector<double> rhobins;
  rhobins.push_back(0);
  rhobins.push_back(2);
  rhobins.push_back(4);
  rhobins.push_back(6);
  rhobins.push_back(8);
  rhobins.push_back(10);
  rhobins.push_back(12);
  rhobins.push_back(14);
  rhobins.push_back(16);
  rhobins.push_back(18);
  rhobins.push_back(20);
  rhobins.push_back(22);
  rhobins.push_back(24);
  rhobins.push_back(26);



  vector<double> ptbins2D;
  ptbins2D.push_back(5);  
  ptbins2D.push_back(10);  
  ptbins2D.push_back(35);  
  vector<double> etabins2D;
  etabins2D.push_back(0.0);
  etabins2D.push_back(1.5);
  etabins2D.push_back(2.5);


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *histLeadingJetPt = new TH1F("leadingJetPt" , "; p_{T} [GeV/c] ; Number of Events ",  200, 0 , 200);

  //3D array, indices give: [denominatorType][SampleType][ptThreshold]

  vector<vector<vector<TH1F*> > > DenominatorVector_Pt;
  vector<vector<vector<TH1F*> > > DenominatorVector_Eta;
  vector<vector<vector<TH1F*> > > DenominatorVector_Phi;
  vector<vector<vector<TH1F*> > > DenominatorVector_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Rho;
  vector<vector<vector<TH2F*> > > DenominatorVector_EtaPt;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt;
  vector<vector<vector<TH1F*> > > NumeratorVector_Eta;
  vector<vector<vector<TH1F*> > > NumeratorVector_Phi;
  vector<vector<vector<TH1F*> > > NumeratorVector_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Rho;
  vector<vector<vector<TH2F*> > > NumeratorVector_EtaPt;
  vector<vector<vector<TH1F*> > > LeptonJetPt;
  vector<vector<vector<TH1F*> > > DenominatorIsolation;

  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Barrel_Rho;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Endcap_Rho;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Barrel_Rho;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Endcap_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Barrel_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Endcap_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Barrel_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Endcap_Rho;


  vector<Double_t> denominatorType;
//   denominatorType.push_back(1);
//   denominatorType.push_back(2);
  denominatorType.push_back(3);
//   denominatorType.push_back(4);
//   denominatorType.push_back(5);
  vector<string> sampleLabel;
  sampleLabel.push_back("ZSample");
  vector<Double_t> ptThreshold;
  ptThreshold.push_back(0);
  //ptThreshold.push_back(5);
  //ptThreshold.push_back(10);
  //ptThreshold.push_back(15);
  
  for (UInt_t denominatorTypeIndex = 0; denominatorTypeIndex < denominatorType.size(); ++denominatorTypeIndex) {
    vector<vector<TH1F*> > tmpDenominatorVector_Pt;
    vector<vector<TH1F*> > tmpDenominatorVector_Eta;
    vector<vector<TH1F*> > tmpDenominatorVector_Phi;
    vector<vector<TH1F*> > tmpDenominatorVector_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Rho;
    vector<vector<TH2F*> > tmpDenominatorVector_EtaPt;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt;
    vector<vector<TH1F*> > tmpNumeratorVector_Eta;
    vector<vector<TH1F*> > tmpNumeratorVector_Phi;
    vector<vector<TH1F*> > tmpNumeratorVector_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Rho;
    vector<vector<TH2F*> > tmpNumeratorVector_EtaPt;
    vector<vector<TH1F*> > tmpLeptonJetPt;
    vector<vector<TH1F*> > tmpDenominatorIsolation;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Barrel_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Barrel_Rho;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Endcap_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Endcap_Rho;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Barrel_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Barrel_Rho;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Endcap_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Endcap_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Barrel_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Barrel_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Endcap_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Endcap_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Barrel_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Barrel_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Endcap_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Endcap_Rho;


    for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
      vector<TH1F*>  tmptmpDenominatorVector_Pt;
      vector<TH1F*>  tmptmpDenominatorVector_Eta;
      vector<TH1F*>  tmptmpDenominatorVector_Phi;
      vector<TH1F*>  tmptmpDenominatorVector_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Rho;
      vector<TH2F*>  tmptmpDenominatorVector_EtaPt;
      vector<TH1F*>  tmptmpNumeratorVector_Pt;
      vector<TH1F*>  tmptmpNumeratorVector_Eta;
      vector<TH1F*>  tmptmpNumeratorVector_Phi;
      vector<TH1F*>  tmptmpNumeratorVector_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Rho;
      vector<TH2F*>  tmptmpNumeratorVector_EtaPt;
      vector<TH1F*>  tmptmpLeptonJetPt;
      vector<TH1F*>  tmptmpDenominatorIsolation;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Barrel_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Barrel_Rho;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Endcap_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Endcap_Rho;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Barrel_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Barrel_Rho;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Endcap_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Endcap_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Barrel_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Barrel_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Endcap_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Endcap_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Barrel_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Barrel_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Endcap_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Endcap_Rho;
      for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {
        TH1F *histDenominator_Pt = new TH1F(("histDenominator_Pt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str(), "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histDenominator_Eta = new TH1F(("histDenominator_Eta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #eta ; Number of Events ",  100, 0 , 3.0);
        TH1F *histDenominator_Phi = new TH1F(("histDenominator_Phi_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #phi ; Number of Events ",  100, -3.2 , 3.2);
        TH1F *histDenominator_NVtx = new TH1F(("histDenominator_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Rho = new TH1F(("histDenominator_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH2F *histDenominator_EtaPt = new TH2F(("histDenominator_EtaPt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #eta ; p_{T} [GeV/c] ; Number of Events ", 100, 0, 3.0,  100, 0 , 100);
        TH1F *histNumerator_Pt = new TH1F(("histNumerator_Pt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histNumerator_Eta = new TH1F(("histNumerator_Eta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #eta ; Number of Events ",  100, 0.0 , 3.0);
        TH1F *histNumerator_Phi = new TH1F(("histNumerator_Phi_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #phi ; Number of Events ",  100, -3.2 , 3.2);
        TH1F *histNumerator_NVtx = new TH1F(("histNumerator_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Rho = new TH1F(("histNumerator_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH2F *histNumerator_EtaPt = new TH2F(("histNumerator_EtaPt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #eta ; p_{T} [GeV/c] ; Number of Events ", 100, 0, 3.0,  100, 0 , 100);        
        TH1F *histLeptonJetPt = new TH1F(("histLeptonJetPt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histDenominatorIsolation = new TH1F(("histDenominatorIsolation_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; PF RelIso ; Number of Events ",  100, 0 , 1.0);

        TH1F *histDenominator_Pt10To20_Barrel_NVtx = new TH1F(("histDenominator_Pt10To20_Barrel_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt10To20_Barrel_Rho = new TH1F(("histDenominator_Pt10To20_Barrel_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histDenominator_Pt10To20_Endcap_NVtx = new TH1F(("histDenominator_Pt10To20_Endcap_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt10To20_Endcap_Rho = new TH1F(("histDenominator_Pt10To20_Endcap_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histDenominator_Pt20ToInf_Barrel_NVtx = new TH1F(("histDenominator_Pt20ToInf_Barrel_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt20ToInf_Barrel_Rho = new TH1F(("histDenominator_Pt20ToInf_Barrel_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histDenominator_Pt20ToInf_Endcap_NVtx = new TH1F(("histDenominator_Pt20ToInf_Endcap_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt20ToInf_Endcap_Rho = new TH1F(("histDenominator_Pt20ToInf_Endcap_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt10To20_Barrel_NVtx = new TH1F(("histNumerator_Pt10To20_Barrel_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt10To20_Barrel_Rho = new TH1F(("histNumerator_Pt10To20_Barrel_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt10To20_Endcap_NVtx = new TH1F(("histNumerator_Pt10To20_Endcap_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt10To20_Endcap_Rho = new TH1F(("histNumerator_Pt10To20_Endcap_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt20ToInf_Barrel_NVtx = new TH1F(("histNumerator_Pt20ToInf_Barrel_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt20ToInf_Barrel_Rho = new TH1F(("histNumerator_Pt20ToInf_Barrel_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt20ToInf_Endcap_NVtx = new TH1F(("histNumerator_Pt20ToInf_Endcap_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt20ToInf_Endcap_Rho = new TH1F(("histNumerator_Pt20ToInf_Endcap_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);



        tmptmpDenominatorVector_Pt.push_back(histDenominator_Pt);
        tmptmpDenominatorVector_Eta.push_back(histDenominator_Eta);
        tmptmpDenominatorVector_Phi.push_back(histDenominator_Phi);
        tmptmpDenominatorVector_NVtx.push_back(histDenominator_NVtx);
        tmptmpDenominatorVector_Rho.push_back(histDenominator_Rho);
        tmptmpDenominatorVector_EtaPt.push_back(histDenominator_EtaPt);
        tmptmpNumeratorVector_Pt.push_back(histNumerator_Pt);
        tmptmpNumeratorVector_Eta.push_back(histNumerator_Eta);
        tmptmpNumeratorVector_Phi.push_back(histNumerator_Phi);
        tmptmpNumeratorVector_NVtx.push_back(histNumerator_NVtx);
        tmptmpNumeratorVector_Rho.push_back(histNumerator_Rho);
        tmptmpNumeratorVector_EtaPt.push_back(histNumerator_EtaPt);
        tmptmpLeptonJetPt.push_back(histLeptonJetPt);
        tmptmpDenominatorIsolation.push_back(histDenominatorIsolation);
        tmptmpDenominatorVector_Pt10To20_Barrel_NVtx.push_back(histDenominator_Pt10To20_Barrel_NVtx);
        tmptmpDenominatorVector_Pt10To20_Barrel_Rho.push_back(histDenominator_Pt10To20_Barrel_Rho);
        tmptmpDenominatorVector_Pt10To20_Endcap_NVtx.push_back(histDenominator_Pt10To20_Endcap_NVtx);
        tmptmpDenominatorVector_Pt10To20_Endcap_Rho.push_back(histDenominator_Pt10To20_Endcap_Rho);
        tmptmpDenominatorVector_Pt20ToInf_Barrel_NVtx.push_back(histDenominator_Pt20ToInf_Barrel_NVtx);
        tmptmpDenominatorVector_Pt20ToInf_Barrel_Rho.push_back(histDenominator_Pt20ToInf_Barrel_Rho);
        tmptmpDenominatorVector_Pt20ToInf_Endcap_NVtx.push_back(histDenominator_Pt20ToInf_Endcap_NVtx);
        tmptmpDenominatorVector_Pt20ToInf_Endcap_Rho.push_back(histDenominator_Pt20ToInf_Endcap_Rho);
        tmptmpNumeratorVector_Pt10To20_Barrel_NVtx.push_back(histNumerator_Pt10To20_Barrel_NVtx);
        tmptmpNumeratorVector_Pt10To20_Barrel_Rho.push_back(histNumerator_Pt10To20_Barrel_Rho);
        tmptmpNumeratorVector_Pt10To20_Endcap_NVtx.push_back(histNumerator_Pt10To20_Endcap_NVtx);
        tmptmpNumeratorVector_Pt10To20_Endcap_Rho.push_back(histNumerator_Pt10To20_Endcap_Rho);
        tmptmpNumeratorVector_Pt20ToInf_Barrel_NVtx.push_back(histNumerator_Pt20ToInf_Barrel_NVtx);
        tmptmpNumeratorVector_Pt20ToInf_Barrel_Rho.push_back(histNumerator_Pt20ToInf_Barrel_Rho);
        tmptmpNumeratorVector_Pt20ToInf_Endcap_NVtx.push_back(histNumerator_Pt20ToInf_Endcap_NVtx);
        tmptmpNumeratorVector_Pt20ToInf_Endcap_Rho.push_back(histNumerator_Pt20ToInf_Endcap_Rho);

      }
      tmpDenominatorVector_Pt.push_back(tmptmpDenominatorVector_Pt);
      tmpDenominatorVector_Eta.push_back(tmptmpDenominatorVector_Eta);
      tmpDenominatorVector_Phi.push_back(tmptmpDenominatorVector_Phi);
      tmpDenominatorVector_NVtx.push_back(tmptmpDenominatorVector_NVtx);
      tmpDenominatorVector_Rho.push_back(tmptmpDenominatorVector_Rho);
      tmpDenominatorVector_EtaPt.push_back(tmptmpDenominatorVector_EtaPt);
      tmpNumeratorVector_Pt.push_back(tmptmpNumeratorVector_Pt);
      tmpNumeratorVector_Eta.push_back(tmptmpNumeratorVector_Eta);
      tmpNumeratorVector_Phi.push_back(tmptmpNumeratorVector_Phi);
      tmpNumeratorVector_NVtx.push_back(tmptmpNumeratorVector_NVtx);
      tmpNumeratorVector_Rho.push_back(tmptmpNumeratorVector_Rho);
      tmpNumeratorVector_EtaPt.push_back(tmptmpNumeratorVector_EtaPt);
      tmpLeptonJetPt.push_back(tmptmpLeptonJetPt);
      tmpDenominatorIsolation.push_back(tmptmpDenominatorIsolation);
      tmpDenominatorVector_Pt10To20_Barrel_NVtx.push_back(tmptmpDenominatorVector_Pt10To20_Barrel_NVtx);
      tmpDenominatorVector_Pt10To20_Barrel_Rho.push_back(tmptmpDenominatorVector_Pt10To20_Barrel_Rho);
      tmpDenominatorVector_Pt10To20_Endcap_NVtx.push_back(tmptmpDenominatorVector_Pt10To20_Endcap_NVtx);
      tmpDenominatorVector_Pt10To20_Endcap_Rho.push_back(tmptmpDenominatorVector_Pt10To20_Endcap_Rho);
      tmpDenominatorVector_Pt20ToInf_Barrel_NVtx.push_back(tmptmpDenominatorVector_Pt20ToInf_Barrel_NVtx);
      tmpDenominatorVector_Pt20ToInf_Barrel_Rho.push_back(tmptmpDenominatorVector_Pt20ToInf_Barrel_Rho);
      tmpDenominatorVector_Pt20ToInf_Endcap_NVtx.push_back(tmptmpDenominatorVector_Pt20ToInf_Endcap_NVtx);
      tmpDenominatorVector_Pt20ToInf_Endcap_Rho.push_back(tmptmpDenominatorVector_Pt20ToInf_Endcap_Rho);
      tmpNumeratorVector_Pt10To20_Barrel_NVtx.push_back(tmptmpNumeratorVector_Pt10To20_Barrel_NVtx);
      tmpNumeratorVector_Pt10To20_Barrel_Rho.push_back(tmptmpNumeratorVector_Pt10To20_Barrel_Rho);
      tmpNumeratorVector_Pt10To20_Endcap_NVtx.push_back(tmptmpNumeratorVector_Pt10To20_Endcap_NVtx);
      tmpNumeratorVector_Pt10To20_Endcap_Rho.push_back(tmptmpNumeratorVector_Pt10To20_Endcap_Rho);
      tmpNumeratorVector_Pt20ToInf_Barrel_NVtx.push_back(tmptmpNumeratorVector_Pt20ToInf_Barrel_NVtx);
      tmpNumeratorVector_Pt20ToInf_Barrel_Rho.push_back(tmptmpNumeratorVector_Pt20ToInf_Barrel_Rho);
      tmpNumeratorVector_Pt20ToInf_Endcap_NVtx.push_back(tmptmpNumeratorVector_Pt20ToInf_Endcap_NVtx);
      tmpNumeratorVector_Pt20ToInf_Endcap_Rho.push_back(tmptmpNumeratorVector_Pt20ToInf_Endcap_Rho);
    }
    DenominatorVector_Pt.push_back(tmpDenominatorVector_Pt);
    DenominatorVector_Eta.push_back(tmpDenominatorVector_Eta);
    DenominatorVector_Phi.push_back(tmpDenominatorVector_Phi);
    DenominatorVector_NVtx.push_back(tmpDenominatorVector_NVtx);
    DenominatorVector_Rho.push_back(tmpDenominatorVector_Rho);
    DenominatorVector_EtaPt.push_back(tmpDenominatorVector_EtaPt);
    NumeratorVector_Pt.push_back(tmpNumeratorVector_Pt);
    NumeratorVector_Eta.push_back(tmpNumeratorVector_Eta);
    NumeratorVector_Phi.push_back(tmpNumeratorVector_Phi);
    NumeratorVector_NVtx.push_back(tmpNumeratorVector_NVtx);
    NumeratorVector_Rho.push_back(tmpNumeratorVector_Rho);
    NumeratorVector_EtaPt.push_back(tmpNumeratorVector_EtaPt);
    LeptonJetPt.push_back(tmpLeptonJetPt);
    DenominatorIsolation.push_back(tmpDenominatorIsolation);
    DenominatorVector_Pt10To20_Barrel_NVtx.push_back(tmpDenominatorVector_Pt10To20_Barrel_NVtx);
    DenominatorVector_Pt10To20_Barrel_Rho.push_back(tmpDenominatorVector_Pt10To20_Barrel_Rho);
    DenominatorVector_Pt10To20_Endcap_NVtx.push_back(tmpDenominatorVector_Pt10To20_Endcap_NVtx);
    DenominatorVector_Pt10To20_Endcap_Rho.push_back(tmpDenominatorVector_Pt10To20_Endcap_Rho);
    DenominatorVector_Pt20ToInf_Barrel_NVtx.push_back(tmpDenominatorVector_Pt20ToInf_Barrel_NVtx);
    DenominatorVector_Pt20ToInf_Barrel_Rho.push_back(tmpDenominatorVector_Pt20ToInf_Barrel_Rho);
    DenominatorVector_Pt20ToInf_Endcap_NVtx.push_back(tmpDenominatorVector_Pt20ToInf_Endcap_NVtx);
    DenominatorVector_Pt20ToInf_Endcap_Rho.push_back(tmpDenominatorVector_Pt20ToInf_Endcap_Rho);
    NumeratorVector_Pt10To20_Barrel_NVtx.push_back(tmpNumeratorVector_Pt10To20_Barrel_NVtx);
    NumeratorVector_Pt10To20_Barrel_Rho.push_back(tmpNumeratorVector_Pt10To20_Barrel_Rho);
    NumeratorVector_Pt10To20_Endcap_NVtx.push_back(tmpNumeratorVector_Pt10To20_Endcap_NVtx);
    NumeratorVector_Pt10To20_Endcap_Rho.push_back(tmpNumeratorVector_Pt10To20_Endcap_Rho);
    NumeratorVector_Pt20ToInf_Barrel_NVtx.push_back(tmpNumeratorVector_Pt20ToInf_Barrel_NVtx);
    NumeratorVector_Pt20ToInf_Barrel_Rho.push_back(tmpNumeratorVector_Pt20ToInf_Barrel_Rho);
    NumeratorVector_Pt20ToInf_Endcap_NVtx.push_back(tmpNumeratorVector_Pt20ToInf_Endcap_NVtx);
    NumeratorVector_Pt20ToInf_Endcap_Rho.push_back(tmpNumeratorVector_Pt20ToInf_Endcap_Rho);
  }

  ofstream eventListFile("eventList.txt");

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
  } else if (inputFilename == "RUN2011A") {
  }
  else {
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
    eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr);         pfcandidateBr = eventTree->GetBranch("PFCandidate");
 
    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
//       cout << "start event " << ientry << " : " << info->runNum << " " << info->lumiSec << " " << info->evtNum << " \n";
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
      
      //For MVA, only use odd event numbers because even event numbers were used for training
      if (Option >= 20 && Option < 30) {
//         if (info->evtNum % 2 == 0) continue;

//         if (Option == 21) if (!(info->nPV0 >= 0 && info->nPV0 <= 2)) continue;
//         if (Option == 22) if (!(info->nPV0 >= 3 && info->nPV0 <= 5)) continue;
//         if (Option == 23) if (!(info->nPV0 >= 6 && info->nPV0 <= 10)) continue;

      }

      Double_t eventweight = info->eventweight;

      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...

      Double_t rho = 0;
      if (!(TMath::IsNaN(info->PileupEnergyDensity) || isinf(info->PileupEnergyDensity))) rho = info->PileupEnergyDensity;

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
      Int_t NLeptons = 0;
      Int_t NMuons = 0;
      vector<const mithep::TMuon*> goodMuons;
      vector<const mithep::TElectron*> goodElectrons;
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
        if (mu->pt > 5 && mu->typeBits > 0) { NLeptons++; NMuons++; }
        if(passMuonID(mu)) goodMuons.push_back(mu);
      }
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
        if (ele->pt > 5) { NLeptons++; }
        if(passCutBasedEleID(ele)) goodElectrons.push_back(ele);
      }

      Int_t NZCandidates = 0;
      const mithep::TMuon* ZMuon1 = 0;
      const mithep::TMuon* ZMuon2 = 0;
      const mithep::TElectron *ZEle1 = 0;
      const mithep::TElectron *ZEle2 = 0;
      Double_t ZPt = 0;
      Double_t ZMass = 0;
      for(UInt_t i=0; i<goodMuons.size(); i++) {
        mithep::FourVectorM mu1;
        mu1.SetCoordinates(goodMuons[i]->pt, goodMuons[i]->eta, goodMuons[i]->phi, 105.658369e-3 );
        for(UInt_t j=i+1; j<goodMuons.size(); j++) {
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
      for(UInt_t i=0; i<goodElectrons.size(); i++) {
        mithep::FourVectorM ele1;
        ele1.SetCoordinates(goodElectrons[i]->pt, goodElectrons[i]->eta, goodElectrons[i]->phi, 0.51099892e-3 );
        for(UInt_t j=i+1; j<goodElectrons.size(); j++) {
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
      // Veto Events with more than 3 lepton candidates
      //********************************************************
      if ((ZDecayType == 11 &&  NMuons > 1) 
          ||
          (ZDecayType == 13 &&  NMuons > 3) 
        ) continue;
      //if (NLeptons > 3) continue;




      //*****************************************************************************
      // Find Electrons and Muons passing minimum ID, for isolation veto
      //*****************************************************************************
      vector<UInt_t> IdentifiedLeptonType;
      vector<Int_t> IdentifiedLeptonIndex;
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);      
        
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
     
        if ( (0==0)
             &&
             ele->pt > 5.0
             && 
             fabs(ele->eta) < 2.5
             &&     
             PassCiCID(ele,1) 
          ) {
          IdentifiedLeptonType.push_back(11);
          IdentifiedLeptonIndex.push_back(i);
        }
      }

      
      //********************************************************
      // Loop over electrons
      //********************************************************
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
        
        //Remove Z Electrons
        if (ZDecayType == 13 && (mu == ZMuon1 || mu == ZMuon2)) continue;
        if (ZEle1 && mithep::MathUtils::DeltaR(mu->eta, mu->phi, ZEle1->eta, ZEle1->phi) < 0.3) continue;
        if (ZEle2 && mithep::MathUtils::DeltaR(mu->eta, mu->phi, ZEle2->eta, ZEle2->phi) < 0.3) continue;
        if (ZMuon1 && mithep::MathUtils::DeltaR(mu->eta, mu->phi, ZMuon1->eta, ZMuon1->phi) < 0.3) continue;
        if (ZMuon2 && mithep::MathUtils::DeltaR(mu->eta, mu->phi, ZMuon2->eta, ZMuon2->phi) < 0.3) continue;

        Double_t leadingJetPt = -1;
        Double_t leptonJetPt = -1;
        //pass event selection     
        for(Int_t j=0; j<jetArr->GetEntries(); j++) {
          const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[j]);        
          if (jet->pt > leadingJetPt &&
              mithep::MathUtils::DeltaR(jet->phi, jet->eta, mu->phi, mu->eta) > 1.0) {
            leadingJetPt = jet->pt;          
          }
          if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, mu->phi, mu->eta) < 0.5) {
            leptonJetPt = jet->pt;
          }
        }
      
        //if there's no jet ( == isolated lepton?) then take pt of lepton
        if (leptonJetPt < 0) {
          leptonJetPt = mu->pt;
        }

        for ( UInt_t denominatorTypeIndex = 0 ; denominatorTypeIndex < denominatorType.size() ; ++denominatorTypeIndex ) {
          for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
            for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {
          
              //********************************************************
              // Event Selection Cuts
              //********************************************************

              if (met > 20) continue;              
              Bool_t PassZPtSelection = kFALSE;
              if (ZPt > ptThreshold[ptThresholdIndex]) {
                PassZPtSelection = kTRUE;               
              }
              if (!PassZPtSelection) continue;      
              if ((mu->pt < 5)) continue;
              if ((mu->pt > 35)) continue;

              if (passMuonDenominatorCuts(mu, denominatorType[denominatorTypeIndex])) {


                DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(mu->pt,eventweight);
                DenominatorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(mu->eta),eventweight);
                DenominatorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(mu->phi,eventweight);
                DenominatorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                DenominatorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                DenominatorVector_EtaPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(mu->eta),mu->pt, eventweight);
              
                if (mu->pt < 20 && fabs(mu->eta) < 1.479) {
                  DenominatorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                } 
                if (mu->pt < 20 && fabs(mu->eta) >= 1.479) {
                  DenominatorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                } 
                if (mu->pt >= 20 && fabs(mu->eta) < 1.479) {
                  DenominatorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                } 
                if (mu->pt >= 20 && fabs(mu->eta) >= 1.479) {
                  DenominatorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                } 
              

                //*************************************************
                //Compute Isolation
                //*************************************************
                Bool_t passIsoMVALooseWP = kFALSE;
                Bool_t passIsoMVATightWP = kFALSE;


                Double_t tmpChargedIso_DR0p0To0p1  = 0;
                Double_t tmpChargedIso_DR0p1To0p2  = 0;
                Double_t tmpChargedIso_DR0p2To0p3  = 0;
                Double_t tmpChargedIso_DR0p3To0p4  = 0;
                Double_t tmpChargedIso_DR0p4To0p5  = 0;
                Double_t tmpChargedIso_DR0p5To0p7  = 0;
                Double_t tmpGammaIso_DR0p0To0p1  = 0;
                Double_t tmpGammaIso_DR0p1To0p2  = 0;
                Double_t tmpGammaIso_DR0p2To0p3  = 0;
                Double_t tmpGammaIso_DR0p3To0p4  = 0;
                Double_t tmpGammaIso_DR0p4To0p5  = 0;
                Double_t tmpGammaIso_DR0p5To0p7  = 0;
                Double_t tmpNeutralHadronIso_DR0p0To0p1  = 0;
                Double_t tmpNeutralHadronIso_DR0p1To0p2  = 0;
                Double_t tmpNeutralHadronIso_DR0p2To0p3  = 0;
                Double_t tmpNeutralHadronIso_DR0p3To0p4  = 0;
                Double_t tmpNeutralHadronIso_DR0p4To0p5  = 0;
                Double_t tmpNeutralHadronIso_DR0p5To0p7  = 0;

                //Loop over PF Candidates
                for(Int_t k=0; k<pfcandidateArr->GetEntries(); ++k) {
                  const mithep::TPFCandidate *pf = (mithep::TPFCandidate*)((*pfcandidateArr)[k]);
                  if (pf->matchedObjectType == 13 && pf->matchedObjectIndex == i) continue;

                  Double_t deta = (mu->eta - pf->eta);
                  Double_t dphi = mithep::MathUtils::DeltaPhi(Double_t(mu->phi),Double_t(pf->phi));
                  Double_t dr = mithep::MathUtils::DeltaR(mu->phi,mu->eta, pf->phi, pf->eta);
                  if (dr > 1.0) continue;

                  Bool_t IsLeptonFootprint = kFALSE;
                  //************************************************************
                  // Lepton Footprint Removal
                  //************************************************************
                  if (dr < 1.0) {
                    //Check for id'ed leptons
                    for (UInt_t q=0; q < IdentifiedLeptonType.size() ; ++q) {
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

                      if (dr < 0.1) tmpChargedIso_DR0p0To0p1 += pf->pt;                  
                      if (dr >= 0.1 && dr < 0.2) tmpChargedIso_DR0p1To0p2 += pf->pt;
                      if (dr >= 0.2 && dr < 0.3) tmpChargedIso_DR0p2To0p3 += pf->pt;
                      if (dr >= 0.3 && dr < 0.4) tmpChargedIso_DR0p3To0p4 += pf->pt;
                      if (dr >= 0.4 && dr < 0.5) tmpChargedIso_DR0p4To0p5 += pf->pt;
                      if (dr >= 0.5 && dr < 0.7) tmpChargedIso_DR0p5To0p7 += pf->pt;
                    }
                    else if (pf->pfType == eGamma) {
                      if (info->runNum == 165993 && info->lumiSec == 1032 && info->evtNum == 1132718343) {cout << "PFGamma: " << dr << " " << pf->pt << endl;}

                      if (dr < 0.1) tmpGammaIso_DR0p0To0p1 += pf->pt;                  
                      if (dr >= 0.1 && dr < 0.2) tmpGammaIso_DR0p1To0p2 += pf->pt;
                      if (dr >= 0.2 && dr < 0.3) tmpGammaIso_DR0p2To0p3 += pf->pt;
                      if (dr >= 0.3 && dr < 0.4) tmpGammaIso_DR0p3To0p4 += pf->pt;
                      if (dr >= 0.4 && dr < 0.5) tmpGammaIso_DR0p4To0p5 += pf->pt;
                      if (dr >= 0.5 && dr < 0.7) tmpGammaIso_DR0p5To0p7 += pf->pt; 
                    }
                    else {
                      if (info->runNum == 165993 && info->lumiSec == 1032 && info->evtNum == 1132718343) {cout << "PFNeutralHad: " << dr << " " << pf->pt << endl;}
                      if (dr < 0.1) tmpNeutralHadronIso_DR0p0To0p1 += pf->pt;                  
                      if (dr >= 0.1 && dr < 0.2) tmpNeutralHadronIso_DR0p1To0p2 += pf->pt;
                      if (dr >= 0.2 && dr < 0.3) tmpNeutralHadronIso_DR0p2To0p3 += pf->pt;
                      if (dr >= 0.3 && dr < 0.4) tmpNeutralHadronIso_DR0p3To0p4 += pf->pt;
                      if (dr >= 0.4 && dr < 0.5) tmpNeutralHadronIso_DR0p4To0p5 += pf->pt;
                      if (dr >= 0.5 && dr < 0.7) tmpNeutralHadronIso_DR0p5To0p7 += pf->pt; 
                    }
                  } //no lepton footprint          
                }

                Int_t EffectiveAreaVersion = -1;  
                EffectiveAreaVersion = 0;  //2011 Data

                varMuPt = mu->pt;
                varMuEta = mu->eta;
                varChargedIso_DR0p0To0p1   = TMath::Min((tmpChargedIso_DR0p0To0p1)/mu->pt, 2.5);
                varChargedIso_DR0p1To0p2   = TMath::Min((tmpChargedIso_DR0p1To0p2)/mu->pt, 2.5);
                varChargedIso_DR0p2To0p3 = TMath::Min((tmpChargedIso_DR0p2To0p3)/mu->pt, 2.5);
                varChargedIso_DR0p3To0p4 = TMath::Min((tmpChargedIso_DR0p3To0p4)/mu->pt, 2.5);
                varChargedIso_DR0p4To0p5 = TMath::Min((tmpChargedIso_DR0p4To0p5)/mu->pt, 2.5);
                varGammaIso_DR0p0To0p1   = TMath::Max(TMath::Min((tmpGammaIso_DR0p0To0p1 - rho*MuonEffectiveArea(kMuGammaIsoDR0p0To0p1, mu->eta, EffectiveAreaVersion))/mu->pt, 2.5), 0.0);
                varGammaIso_DR0p1To0p2   = TMath::Max(TMath::Min((tmpGammaIso_DR0p1To0p2 - rho*MuonEffectiveArea(kMuGammaIsoDR0p1To0p2, mu->eta, EffectiveAreaVersion))/mu->pt, 2.5), 0.0);
                varGammaIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpGammaIso_DR0p2To0p3 - rho*MuonEffectiveArea(kMuGammaIsoDR0p2To0p3, mu->eta, EffectiveAreaVersion))/mu->pt, 2.5), 0.0);
                varGammaIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpGammaIso_DR0p3To0p4 - rho*MuonEffectiveArea(kMuGammaIsoDR0p3To0p4, mu->eta, EffectiveAreaVersion))/mu->pt, 2.5), 0.0);
                varGammaIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpGammaIso_DR0p4To0p5 - rho*MuonEffectiveArea(kMuGammaIsoDR0p4To0p5, mu->eta, EffectiveAreaVersion))/mu->pt, 2.5), 0.0);
                varNeutralHadronIso_DR0p0To0p1   = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p0To0p1 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p0To0p1, mu->eta, EffectiveAreaVersion))/mu->pt, 2.5), 0.0);
                varNeutralHadronIso_DR0p1To0p2   = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p1To0p2 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p1To0p2, mu->eta, EffectiveAreaVersion))/mu->pt, 2.5), 0.0);
                varNeutralHadronIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p2To0p3 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p2To0p3, mu->eta, EffectiveAreaVersion))/mu->pt, 2.5), 0.0);
                varNeutralHadronIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p3To0p4 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p3To0p4, mu->eta, EffectiveAreaVersion))/mu->pt, 2.5), 0.0);
                varNeutralHadronIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p4To0p5 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p4To0p5, mu->eta, EffectiveAreaVersion))/mu->pt, 2.5), 0.0);

                if (info->runNum == 165993 && info->lumiSec == 1032 && info->evtNum == 1132718343) {
                  cout << varChargedIso_DR0p0To0p1 << " " 
                       << varChargedIso_DR0p1To0p2 << " " 
                       << varChargedIso_DR0p2To0p3 << " " 
                       << varChargedIso_DR0p3To0p4 << " " 
                       << varChargedIso_DR0p4To0p5 << " " 
                       << varGammaIso_DR0p0To0p1   << " " 
                       << varGammaIso_DR0p1To0p2   << " " 
                       << varGammaIso_DR0p2To0p3   << " " 
                       << varGammaIso_DR0p3To0p4   << " " 
                       << varGammaIso_DR0p4To0p5   << " " 
                       << varNeutralHadronIso_DR0p0To0p1   << " " 
                       << varNeutralHadronIso_DR0p1To0p2   << " " 
                       << varNeutralHadronIso_DR0p2To0p3 << " " 
                       << varNeutralHadronIso_DR0p3To0p4 << " " 
                       << varNeutralHadronIso_DR0p4To0p5 << " " 
                       << endl;
                }


                Bool_t varIsGlobal = ((mu->typeBits & 1) == 1);
                Bool_t varIsTracker = ((mu->typeBits & 2) == 2);

                //evaluate MVA
                Int_t isomvabin = -1;
                if (varIsGlobal && varIsTracker && fabs(mu->eta) < 1.5 && mu->pt < 10) isomvabin = 0;
                else if (varIsGlobal && varIsTracker && fabs(mu->eta) < 1.5 && mu->pt >= 10) isomvabin = 1;
                else if (varIsGlobal && varIsTracker && fabs(mu->eta) >= 1.5 && mu->pt < 10) isomvabin = 2;
                else if (varIsGlobal && varIsTracker && fabs(mu->eta) >= 1.5 && mu->pt >= 10) isomvabin = 3;
                else if (!varIsGlobal && varIsTracker) isomvabin = 4;
                else if (varIsGlobal && !varIsTracker) isomvabin = 5;
                else { cout << "error\n"; assert(0);}
                Double_t MuIsoMVA = muonIsoMVA->MVAValue_IsoRings(
                  varMuPt,
                  varMuEta,
                  varChargedIso_DR0p0To0p1, 
                  varChargedIso_DR0p1To0p2, 
                  varChargedIso_DR0p2To0p3, 
                  varChargedIso_DR0p3To0p4, 
                  varChargedIso_DR0p4To0p5, 
                  varGammaIso_DR0p0To0p1, 
                  varGammaIso_DR0p1To0p2, 
                  varGammaIso_DR0p2To0p3, 
                  varGammaIso_DR0p3To0p4, 
                  varGammaIso_DR0p4To0p5, 
                  varNeutralHadronIso_DR0p0To0p1, 
                  varNeutralHadronIso_DR0p1To0p2, 
                  varNeutralHadronIso_DR0p2To0p3, 
                  varNeutralHadronIso_DR0p3To0p4, 
                  varNeutralHadronIso_DR0p4To0p5,
                  kFALSE
                  );                
                //Double_t MuIsoMVA = MuonIsoReader[isomvabin]->EvaluateMVA( "MuonIsoMVA_V0_BDTG method" );            
                if (passMuIsoMVA_LooseWP(mu->pt, mu->eta, varIsGlobal , varIsTracker, MuIsoMVA) 
                    && (varChargedIso_DR0p0To0p1 + varChargedIso_DR0p1To0p2 + varChargedIso_DR0p2To0p3 < 0.7)
                  ) {
                  passIsoMVALooseWP = kTRUE;
                }
                if (passMuIsoMVA_TightWP(mu->pt, mu->eta, varIsGlobal , varIsTracker, MuIsoMVA)
                    && (varChargedIso_DR0p0To0p1 + varChargedIso_DR0p1To0p2 + varChargedIso_DR0p2To0p3 < 0.7)
                  ) {
                  passIsoMVATightWP = kTRUE;
                }
              
        
                //*************************************************
                //Compute ID MVA
                //*************************************************
                Bool_t passIDMVALooseWP = kFALSE;
                Bool_t passIDMVATightWP = kFALSE;

                varMuTkNchi2   = mu->tkNchi2;      
                varMuGlobalNchi2  = mu->muNchi2;      
                varMuNValidHits   =  mu->nValidHits;   
                varMuNTrackerHits   = mu->nTkHits;   
                varMuNPixelHits = mu->nPixHits;     
                varMuNMatches   = mu->nMatch;   
                varMuTrkKink     = mu->TrkKink; 
                varMuSegmentCompatibility   = mu->SegmentCompatibility;     
                varMuCaloCompatibility    =    mu->CaloCompatilibity;
                varMuHadEnergy = mu->HadEnergy;
                varMuEmEnergy     = mu->EmEnergy;
                varMuHadS9Energy   = mu->HadS9Energy;       
                varMuEmS9Energy    = mu->EmS9Energy;   		
                
                varMuPt = mu->pt;
                varMuEta = mu->eta;

 
                //evaluate MVA
                Int_t idmvabin = -1;
                if (varIsGlobal && varIsTracker && fabs(mu->eta) < 1.5 && mu->pt < 10) idmvabin = 0;
                else if (varIsGlobal && varIsTracker && fabs(mu->eta) < 1.5 && mu->pt >= 10) idmvabin = 1;
                else if (varIsGlobal && varIsTracker && fabs(mu->eta) >= 1.5 && mu->pt < 10) idmvabin = 2;
                else if (varIsGlobal && varIsTracker && fabs(mu->eta) >= 1.5 && mu->pt >= 10) idmvabin = 3;
                else if (!varIsGlobal && varIsTracker) idmvabin = 4;
                else if (varIsGlobal && !varIsTracker) idmvabin = 5;
                else { cout << "error\n"; assert(0);}
                Double_t MuIDMVA = muonIDMVA->MVAValue_ID( varMuPt,
                                                           varMuEta,
                                                           varIsGlobal,
                                                           varIsTracker,
                                                           varMuTkNchi2,
                                                           varMuGlobalNchi2,
                                                           varMuNValidHits,
                                                           varMuNTrackerHits,
                                                           varMuNPixelHits,
                                                           varMuNMatches,
                                                           varMuTrkKink,
                                                           varMuSegmentCompatibility,
                                                           varMuCaloCompatibility,
                                                           varMuHadEnergy,
                                                           varMuEmEnergy,
                                                           varMuHadS9Energy,
                                                           varMuEmS9Energy,
                                                           kFALSE);
                //Double_t MuIDMVA = MuonIDReader[idmvabin]->EvaluateMVA( "MuonIDMVA_V0_BDTG method" );            
                if (passMuIDMVA_LooseWP(mu->pt, mu->eta, varIsGlobal , varIsTracker, MuIDMVA)) {
                  passIDMVALooseWP = kTRUE;
                }
                if (passMuIDMVA_TightWP(mu->pt, mu->eta, varIsGlobal , varIsTracker, MuIDMVA)) {
                  passIDMVATightWP = kTRUE;
                }
                                             
                if (mu->pt >= 10 && fabs(mu->eta) < 1.5 && ((mu->typeBits & 1) == 1) && ((mu->typeBits & 2) == 2)) {
                  eventListFile << info->runNum << " " 
                                << info->lumiSec << " " 
                                << info->evtNum<< " "
                                << " : "
                                << mu->pt << " " << mu->eta << " " << mu->phi << " : "
                                << passMuonID_HZZ2011(mu) << " " 
                                << passMuonIP_HZZ2011(mu) << " " 
                                << PassMuDetIso025(mu,rho) << " " 
                                << MuIDMVA << "," << MuIsoMVA << " : " 

                                << varChargedIso_DR0p0To0p1   << " "
                                << varChargedIso_DR0p1To0p2   << " "
                                << varChargedIso_DR0p2To0p3 << " "
                                << varChargedIso_DR0p3To0p4 << " "
                                << varChargedIso_DR0p4To0p5 << " "
                                << varGammaIso_DR0p0To0p1   << " "
                                << varGammaIso_DR0p1To0p2   << " "
                                << varGammaIso_DR0p2To0p3 << " "
                                << varGammaIso_DR0p3To0p4 << " "
                                << varGammaIso_DR0p4To0p5 << " "
                                << varNeutralHadronIso_DR0p0To0p1   << " "
                                << varNeutralHadronIso_DR0p1To0p2   << " "
                                << varNeutralHadronIso_DR0p2To0p3 << " "
                                << varNeutralHadronIso_DR0p3To0p4 << " "
                                << varNeutralHadronIso_DR0p4To0p5 << " : "
                    
                                << Bool_t(passMuonID_HZZ2011(mu) && passMuonIP_HZZ2011(mu) && PassMuDetIso025(mu,rho)) << " "
                                << Bool_t(passIDMVALooseWP && passIsoMVALooseWP)
                                << endl;
                }



                if (info->evtNum == 122200967) {
                  cout << "DEBUG\n";
                  cout << "DONE DEBUG\n";
                }



                if (
                  (Option == 1 && (passMuonID_HZZ2011(mu) &&  passMuonIP_HZZ2011(mu) 
                                   && PassMuDetIso025(mu,rho)
                    ))
                  ||
                  (Option == 10 && passIDMVALooseWP && passIsoMVALooseWP && fabs(mu->ip3dSig) < 4.0)
                  ||
                  (Option == 11 && passIDMVATightWP && passIsoMVATightWP && fabs(mu->ip3dSig) < 4.0)
                  ) {
                  NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(mu->pt, eventweight);
                  NumeratorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(mu->eta), eventweight);
                  NumeratorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(mu->phi, eventweight);
                  NumeratorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  NumeratorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  NumeratorVector_EtaPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(mu->eta), mu->pt , eventweight);   
        
                  if (mu->pt < 20 && fabs(mu->eta) < 1.479) {
                    NumeratorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  }
                  if (mu->pt < 20 && fabs(mu->eta) >= 1.479) {
                    NumeratorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  } 
                  if (mu->pt >= 20 && fabs(mu->eta) < 1.479) {
                    NumeratorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  } 
                  if (mu->pt >= 20 && fabs(mu->eta) >= 1.479) {
                    NumeratorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  }
                }
              }

            } //loop over ptThresholds
          } //loop over sample types
        } //loop over denominator types

      } //loop over electrons

//      cout << "dddata " << " : " << info->runNum << " " << info->lumiSec << " " << info->evtNum << " ";

    } //end loop over data    
    cout << "done " << inputfiles[f] << endl;

  } //end loop over files

  eventListFile.close();

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;
  
  
  //*****************************************************************************************
  //Make Efficiency Plots
  //*****************************************************************************************
  
  for (UInt_t denominatorTypeIndex = 0; denominatorTypeIndex < denominatorType.size(); ++denominatorTypeIndex) {
    for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
      for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {

        Bool_t printDebug = kFALSE;
        if (denominatorType[denominatorTypeIndex] == 4 
            // && sampleLabel[sampleTypeIndex] == "Ele8CaloIdLCaloIsoVLSample" && ptThreshold[ptThresholdIndex] == 30
          ) printDebug = kTRUE;
        cout << label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_EtaPt" << endl;
        
        Int_t ErrorType = 2; //Clopper Pearson errors
        TGraphAsymmErrors *efficiency_pt = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt", ptbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_eta = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Eta", etabins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_phi = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Phi", phibins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_nvtx = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_rho = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        mithep::TH2DAsymErr *efficiency_EtaPt = 
          mithep::EfficiencyUtils::createEfficiencyHist2D(NumeratorVector_EtaPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
                                                          DenominatorVector_EtaPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_EtaPt",etabins2D, ptbins2D,  ErrorType, printDebug);
        

        TGraphAsymmErrors *efficiency_Pt10To20_Barrel_nvtx = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Barrel_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt10To20_Barrel_rho = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Barrel_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt10To20_Endcap_nvtx = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Endcap_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt10To20_Endcap_rho = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Endcap_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Barrel_nvtx = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Barrel_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Barrel_rho = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Barrel_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Endcap_nvtx = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Endcap_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Endcap_rho = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Endcap_Rho", rhobins, ErrorType, -99, -99, 0, 1);

        //rebin hists
        NumeratorVector_EtaPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex] = mithep::PlotUtils::rebin(NumeratorVector_EtaPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex],etabins2D, ptbins2D);
        DenominatorVector_EtaPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex] = mithep::PlotUtils::rebin(DenominatorVector_EtaPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex],etabins2D, ptbins2D);

        TFile *file = new TFile(outputFilename.c_str(), "UPDATE");
        file->cd();
        
        file->WriteTObject(efficiency_pt, efficiency_pt->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_eta, efficiency_eta->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_phi, efficiency_phi->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_nvtx, efficiency_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_rho, efficiency_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_EtaPt, efficiency_EtaPt->GetName(), "WriteDelete");
        file->WriteTObject(NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(NumeratorVector_EtaPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], NumeratorVector_EtaPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(DenominatorVector_EtaPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_EtaPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(LeptonJetPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], LeptonJetPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(DenominatorIsolation[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorIsolation[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        
        file->WriteTObject(efficiency_Pt10To20_Barrel_nvtx, efficiency_Pt10To20_Barrel_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt10To20_Barrel_rho, efficiency_Pt10To20_Barrel_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt10To20_Endcap_nvtx, efficiency_Pt10To20_Endcap_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt10To20_Endcap_rho, efficiency_Pt10To20_Endcap_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Barrel_nvtx, efficiency_Pt20ToInf_Barrel_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Barrel_rho, efficiency_Pt20ToInf_Barrel_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Endcap_nvtx, efficiency_Pt20ToInf_Endcap_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Endcap_rho, efficiency_Pt20ToInf_Endcap_rho->GetName(), "WriteDelete");

        file->Close();
        
        //*****************************************************
        // Combined Fake Rate
        //*****************************************************
        if (denominatorType[denominatorTypeIndex] == 4 && 
            (sampleLabel[sampleTypeIndex] == "Ele8CaloIdLCaloIsoVLCombinedSample"
              ) 
          ) {
          
          TH2F *eff = 0;
          TH2F *effErrorLow = 0;
          TH2F *effErrorHigh = 0;
          
          file = new TFile((smurfOutputFilename + ".root").c_str(), "UPDATE");

          mithep::EfficiencyUtils::createEfficiencyHist2D(NumeratorVector_EtaPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
                                                          DenominatorVector_EtaPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
                                                          "ElectronFakeRate_V"+IntToString(denominatorType[denominatorTypeIndex])+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_EtaPt", 
                                                          ptbins2D, etabins2D, ErrorType, file);
        
          file->Close();
          delete file;
        }

      }
    }
  }



  gBenchmark->Show("WWTemplate");       
} 




Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu, Int_t DenominatorType) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(fabs(mu->eta) < 2.5)) pass = kFALSE;


  if (DenominatorType == 1) {
    //Barrel 
    if(!( fabs(mu->dz) < 0.1
          && ( ((mu->typeBits & kGlobal) == kGlobal) 
               || ( ((mu->typeBits & kTracker) == kTracker)
                    && ((mu->qualityBits & kAllArbitrated) == kAllArbitrated)
                 )
            )
          && mu->trkIso03/mu->pt < 0.7
         )
      ) {
      pass = kFALSE;
    }
  }

  if (DenominatorType == 2) {
    //Barrel 
    if(!( fabs(mu->dz) < 0.1   
          && ( ((mu->typeBits & kGlobal) == kGlobal) 
               || ( ((mu->typeBits & kTracker) == kTracker)
                    && ((mu->qualityBits & kAllArbitrated) == kAllArbitrated)
                 )
            )               
          && fabs(mu->ip3dSig) < 4.0
          && mu->trkIso03/mu->pt < 0.7
         )
      ) {
      pass = kFALSE;
    }
  }

  if (DenominatorType == 3) {
    //Barrel 
    if(!( fabs(mu->dz) < 0.1   
          && ( ((mu->typeBits & kGlobal) == kGlobal) 
               && ((mu->typeBits & kTracker) == kTracker)
            )               
          && fabs(mu->ip3dSig) < 4.0
          && mu->trkIso03/mu->pt < 0.7
         )
      ) {
      pass = kFALSE;
    }
  }

  if (DenominatorType == 4) {
    //Barrel 
    if(!( fabs(mu->dz) < 0.1   
          && ( !((mu->typeBits & kGlobal) == kGlobal) 
               && ((mu->typeBits & kTracker) == kTracker)
            )               
          && fabs(mu->ip3dSig) < 4.0
          && mu->trkIso03/mu->pt < 0.7
         )
      ) {
      pass = kFALSE;
    }
  }
  if (DenominatorType == 5) {
    //Barrel 
    if(!( fabs(mu->dz) < 0.1   
          && ( ((mu->typeBits & kGlobal) == kGlobal) 
               && !((mu->typeBits & kTracker) == kTracker)
            )               
          && fabs(mu->ip3dSig) < 4.0
          && mu->trkIso03/mu->pt < 0.7
         )
      ) {
      pass = kFALSE;
    }
  }
  
  return pass;
}


