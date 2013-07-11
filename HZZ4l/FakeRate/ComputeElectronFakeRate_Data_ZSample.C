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
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"

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
Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele, Int_t DenominatorType);

void DoComputeElectronFakeRate(const string inputFilename,
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
void ComputeElectronFakeRate_Data_ZSample(Int_t Option=0) {



//------------------
//Fully Inclusive
//------------------
  if (Option==0) {
     DoComputeElectronFakeRate("LIST","ElectronFakeRate","ElectronFakeRate.CiCTightDetIso.ZSample.root", "FakeRates_Electron_CiCTightDetIso_ZSample", 1);
  } 
  if (Option==1) {
     DoComputeElectronFakeRate("LIST","ElectronFakeRate","ElectronFakeRate.HZZIDMVAAndIsoMVALoose.ZSample.root", "FakeRates_Electron_HZZIDMVAAndIsoMVALoose_ZSample", 10);     
  } 
  if (Option==2) {
     DoComputeElectronFakeRate("LIST","ElectronFakeRate","ElectronFakeRate.HZZIDMVAAndIsoMVATight.ZSample.root", "FakeRates_Electron_HZZIDMVAAndIsoMVATight_ZSample", 11);     
  } 

  if (Option==3) {
     DoComputeElectronFakeRate("LIST","ElectronFakeRate","ElectronFakeRate.MVAIDDetIso.ZSample.root", "FakeRates_Electron_MVAIDDetIso_ZSample", 3);     
  }
  if (Option==4) {
     DoComputeElectronFakeRate("LIST","ElectronFakeRate","ElectronFakeRate.CiCTightPFIso.ZSample.root", "FakeRates_Electron_CiCTightPFIso_ZSample", 4);     
  }
  if (Option==5) {
     DoComputeElectronFakeRate("LIST","ElectronFakeRate","ElectronFakeRate.CiCTightMVAIso.ZSample.root", "FakeRates_Electron_CiCTightMVAIso_ZSample", 5);
  }



  if (Option==10) {
     DoComputeElectronFakeRate("2012Data","ElectronFakeRate","ElectronFakeRate.CiCTightDetIso.ZSample.2012Data.root", "FakeRates_Electron_CiCTightDetIso_ZSample_2012Data", 1);
  } 
  if (Option==11) {
     DoComputeElectronFakeRate("2012Data","ElectronFakeRate","ElectronFakeRate.HZZIDMVAAndIsoMVALoose.2012Data.ZSample.root", "FakeRates_Electron_HZZIDMVAAndIsoMVALoose_ZSample_2012Data", 10);     
  } 
  if (Option==12) {
     DoComputeElectronFakeRate("2012Data","ElectronFakeRate","ElectronFakeRate.HZZIDMVAAndIsoMVATight.2012Data.ZSample.root", "FakeRates_Electron_HZZIDMVAAndIsoMVATight_ZSample_2012Data", 11);     
  }


}



void DoComputeElectronFakeRate(const string inputFilename,
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
  mithep::ElectronIDMVA *electronIDMVA = new mithep::ElectronIDMVA();
  vector<string> electronIDMVAWeightFiles;
  electronIDMVAWeightFiles.push_back("MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml");
  electronIDMVAWeightFiles.push_back("MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml");
  electronIDMVAWeightFiles.push_back("MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml");
  electronIDMVAWeightFiles.push_back("MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml");
  electronIDMVAWeightFiles.push_back("MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml");
  electronIDMVAWeightFiles.push_back("MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml");
  electronIDMVA->Initialize( "ElectronIDMVA", mithep::ElectronIDMVA::kIDEGamma2012NonTrigV0,
                             kTRUE, electronIDMVAWeightFiles);

  mithep::ElectronIDMVA *electronIsoMVA = new mithep::ElectronIDMVA();
  vector<string> electronIsoMVAWeightFiles;
  electronIsoMVAWeightFiles.push_back("MitPhysics/data/ElectronMVAWeights/ElectronIso_BDTG_V0_BarrelPt5To10.weights.xml");
  electronIsoMVAWeightFiles.push_back("MitPhysics/data/ElectronMVAWeights/ElectronIso_BDTG_V0_EndcapPt5To10.weights.xml");
  electronIsoMVAWeightFiles.push_back("MitPhysics/data/ElectronMVAWeights/ElectronIso_BDTG_V0_BarrelPt10ToInf.weights.xml");
  electronIsoMVAWeightFiles.push_back("MitPhysics/data/ElectronMVAWeights/ElectronIso_BDTG_V0_EndcapPt10ToInf.weights.xml");
  electronIsoMVA->Initialize( "ElectronIsoMVA", mithep::ElectronIDMVA::kIsoRingsV0,
                              kTRUE, electronIsoMVAWeightFiles);

 

  Float_t varfbrem;
  Float_t vardeta;
  Float_t vardphi;
  Float_t varsee;
  Float_t varetawidth;
  Float_t varphiwidth;
  Float_t varHoE;
  Float_t varEoP;
  Float_t vare1x5e5x5;
  Float_t varEoPout;
  Float_t vardetacalo;
  Float_t varkfchi2;
  Float_t varkfhits;
  Float_t varspp;
  Float_t varIoEmIoP;
  Float_t varR9;
  Float_t vargsfchi2;
  Float_t varPreShowerOverRaw;
  Float_t varkfhitsall;
  Float_t vareleEoPout;
  Float_t varOneMinuse1x5e5x5;
  Float_t varIoEmIoP_Emanuele;
  Float_t vardetaabs;
  Float_t vardphiabs;
  Float_t vardetacaloabs;
  Float_t varChargedIso_DR0p0To0p1;
  Float_t varChargedIso_DR0p1To0p2;
  Float_t varChargedIso_DR0p2To0p3;
  Float_t varChargedIso_DR0p3To0p4;
  Float_t varChargedIso_DR0p4To0p5;
  Float_t varChargedIso_DR0p5To0p7;
  Float_t varGammaIso_DR0p0To0p1;
  Float_t varGammaIso_DR0p1To0p2;
  Float_t varGammaIso_DR0p2To0p3;
  Float_t varGammaIso_DR0p3To0p4;
  Float_t varGammaIso_DR0p4To0p5;
  Float_t varGammaIso_DR0p5To0p7;
  Float_t varNeutralHadronIso_DR0p0To0p1;
  Float_t varNeutralHadronIso_DR0p1To0p2;
  Float_t varNeutralHadronIso_DR0p2To0p3;
  Float_t varNeutralHadronIso_DR0p3To0p4;
  Float_t varNeutralHadronIso_DR0p4To0p5;
  Float_t varNeutralHadronIso_DR0p5To0p7;
  Float_t varPt, varEta;

  TMVA::Reader  *ElectronIDReader[6];
  for(UInt_t j=0; j<6; ++j) {
    ElectronIDReader[j] = new TMVA::Reader( "!Color:!Silent:Error" );  
  }
  TMVA::Reader  *ElectronIsoReader[4];
  for(UInt_t j=0; j<4; ++j) {
    ElectronIsoReader[j] = new TMVA::Reader( "!Color:!Silent:Error" );  
  }

  for(UInt_t j=0; j<4; ++j) {
    ElectronIsoReader[j]->AddVariable( "ChargedIso_DR0p0To0p1",         &varChargedIso_DR0p0To0p1           );
    ElectronIsoReader[j]->AddVariable( "ChargedIso_DR0p1To0p2",         &varChargedIso_DR0p1To0p2           );
    ElectronIsoReader[j]->AddVariable( "ChargedIso_DR0p2To0p3",       &varChargedIso_DR0p2To0p3         );
    ElectronIsoReader[j]->AddVariable( "ChargedIso_DR0p3To0p4",        &varChargedIso_DR0p3To0p4          );
    ElectronIsoReader[j]->AddVariable( "ChargedIso_DR0p4To0p5",        &varChargedIso_DR0p4To0p5          );
    ElectronIsoReader[j]->AddVariable( "GammaIso_DR0p0To0p1",         &varGammaIso_DR0p0To0p1           );
    ElectronIsoReader[j]->AddVariable( "GammaIso_DR0p1To0p2",         &varGammaIso_DR0p1To0p2           );
    ElectronIsoReader[j]->AddVariable( "GammaIso_DR0p2To0p3",       &varGammaIso_DR0p2To0p3         );
    ElectronIsoReader[j]->AddVariable( "GammaIso_DR0p3To0p4",        &varGammaIso_DR0p3To0p4          );
    ElectronIsoReader[j]->AddVariable( "GammaIso_DR0p4To0p5",          &varGammaIso_DR0p4To0p5          );
    ElectronIsoReader[j]->AddVariable( "NeutralHadronIso_DR0p0To0p1",         &varNeutralHadronIso_DR0p0To0p1           );
    ElectronIsoReader[j]->AddVariable( "NeutralHadronIso_DR0p1To0p2",         &varNeutralHadronIso_DR0p1To0p2           );
    ElectronIsoReader[j]->AddVariable( "NeutralHadronIso_DR0p2To0p3",       &varNeutralHadronIso_DR0p2To0p3         );
    ElectronIsoReader[j]->AddVariable( "NeutralHadronIso_DR0p3To0p4",        &varNeutralHadronIso_DR0p3To0p4          );
    ElectronIsoReader[j]->AddVariable( "NeutralHadronIso_DR0p4To0p5",        &varNeutralHadronIso_DR0p4To0p5          );
    ElectronIsoReader[j]->AddSpectator( "eta",   &varEta );
    ElectronIsoReader[j]->AddSpectator( "pt" ,   &varPt  );
  }
  for(UInt_t j=0; j<4; ++j) {
    if (j==0) ElectronIsoReader[j]->BookMVA( "ElectronIsoMVA_BDTG_V0 method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/ElectronIsoMVAWeights_TrainedWithWJets/ElectronIsoMVA_V0_BarrelPt5To10_V0_BDTG.weights.xml");   
    if (j==1) ElectronIsoReader[j]->BookMVA( "ElectronIsoMVA_BDTG_V0 method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/ElectronIsoMVAWeights_TrainedWithWJets/ElectronIsoMVA_V0_EndcapPt5To10_V0_BDTG.weights.xml");   
    if (j==2) ElectronIsoReader[j]->BookMVA( "ElectronIsoMVA_BDTG_V0 method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/ElectronIsoMVAWeights_TrainedWithWJets/ElectronIsoMVA_V0_BarrelPt10ToInf_V0_BDTG.weights.xml");   
    if (j==3) ElectronIsoReader[j]->BookMVA( "ElectronIsoMVA_BDTG_V0 method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/ElectronIsoMVAWeights_TrainedWithWJets/ElectronIsoMVA_V0_EndcapPt10ToInf_V0_BDTG.weights.xml");   
  }

  for(UInt_t j=0; j<6; ++j) {
    //****************
    //V1
    //****************
    ElectronIDReader[j]->AddVariable( "fbrem", &varfbrem );
    ElectronIDReader[j]->AddVariable( "kfchi2", &varkfchi2 );
    ElectronIDReader[j]->AddVariable( "kfhitsall", &varkfhitsall );
    ElectronIDReader[j]->AddVariable( "gsfchi2", &vargsfchi2 );
    ElectronIDReader[j]->AddVariable( "deta", &vardetaabs );
    ElectronIDReader[j]->AddVariable( "dphi", &vardphiabs );
    ElectronIDReader[j]->AddVariable( "detacalo", &vardetacaloabs );
    ElectronIDReader[j]->AddVariable( "see", &varsee );
    ElectronIDReader[j]->AddVariable( "spp", &varspp );
    ElectronIDReader[j]->AddVariable( "etawidth", &varetawidth );
    ElectronIDReader[j]->AddVariable( "phiwidth", &varphiwidth );
    ElectronIDReader[j]->AddVariable( "e1x5e5x5", &varOneMinuse1x5e5x5 );
    ElectronIDReader[j]->AddVariable( "R9", &varR9 );
    ElectronIDReader[j]->AddVariable( "HoE", &varHoE );
    ElectronIDReader[j]->AddVariable( "EoP", &varEoP );
    ElectronIDReader[j]->AddVariable( "IoEmIoP", &varIoEmIoP_Emanuele );
    ElectronIDReader[j]->AddVariable( "EoPout", &varEoPout );
    if (j==2 || j==5) {
      ElectronIDReader[j]->AddVariable( "PreShowerOverRaw", &varPreShowerOverRaw );
    }
    ElectronIDReader[j]->AddSpectator( "eta",   &varEta );
    ElectronIDReader[j]->AddSpectator( "pt" ,   &varPt  );
  }

  for(UInt_t j=0; j<6; ++j) {
    if (j==0) ElectronIDReader[j]->BookMVA( "ElectronIDMVA_BDTG_V1 method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/SavedWeights20120418/ElectronIDNonTrigMVAAlternativeWeights/Electrons_BDTG_NonTrigV1_Cat1.weights.xml");                
    if (j==1) ElectronIDReader[j]->BookMVA( "ElectronIDMVA_BDTG_V1 method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/SavedWeights20120418/ElectronIDNonTrigMVAAlternativeWeights/Electrons_BDTG_NonTrigV1_Cat2.weights.xml");                
    if (j==2) ElectronIDReader[j]->BookMVA( "ElectronIDMVA_BDTG_V1 method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/SavedWeights20120418/ElectronIDNonTrigMVAAlternativeWeights/Electrons_BDTG_NonTrigV1_Cat3.weights.xml");                
    if (j==3) ElectronIDReader[j]->BookMVA( "ElectronIDMVA_BDTG_V1 method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/SavedWeights20120418/ElectronIDNonTrigMVAAlternativeWeights/Electrons_BDTG_NonTrigV1_Cat4.weights.xml");                
    if (j==4) ElectronIDReader[j]->BookMVA( "ElectronIDMVA_BDTG_V1 method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/SavedWeights20120418/ElectronIDNonTrigMVAAlternativeWeights/Electrons_BDTG_NonTrigV1_Cat5.weights.xml");
    if (j==5) ElectronIDReader[j]->BookMVA( "ElectronIDMVA_BDTG_V1 method", "/data/blue/sixie/HZZ4l/LeptonSelection/data/SavedWeights20120418/ElectronIDNonTrigMVAAlternativeWeights/Electrons_BDTG_NonTrigV1_Cat6.weights.xml");
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
  etabins2D.push_back(0.8);
  etabins2D.push_back(1.479);
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
  //denominatorType.push_back(1);
  denominatorType.push_back(2);
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
  rlrm.AddJSONFile("/data/blue/sixie/HZZ4l/auxiliar/2012/Cert_190456-191276_8TeV_PromptReco_Collisions12_JSON.txt"); 

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
  } else if (inputFilename == "2012Data") {
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r12a-del-pr-v1_noskim.ZPlusFakeSkimmed.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim.ZPlusFakeSkimmed.root");
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
      Int_t NElectrons = 0;
      vector<const mithep::TMuon*> goodMuons;
      vector<const mithep::TElectron*> goodElectrons;
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
        if (mu->pt > 5) NLeptons++;
        if(passMuonID(mu)) goodMuons.push_back(mu);
      }
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
        if (ele->pt > 5) { NLeptons++; NElectrons++; }
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
      if ((ZDecayType == 11 &&  NElectrons > 3) 
          ||
          (ZDecayType == 13 &&  NElectrons > 1) 
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
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
        
        //Remove Z Electrons
        if (ZDecayType == 11 && (ele == ZEle1 || ele == ZEle2)) continue;
        if (ZEle1 && mithep::MathUtils::DeltaR(ele->eta, ele->phi, ZEle1->eta, ZEle1->phi) < 0.3) continue;
        if (ZEle2 && mithep::MathUtils::DeltaR(ele->eta, ele->phi, ZEle2->eta, ZEle2->phi) < 0.3) continue;
        if (ZMuon1 && mithep::MathUtils::DeltaR(ele->eta, ele->phi, ZMuon1->eta, ZMuon1->phi) < 0.3) continue;
        if (ZMuon2 && mithep::MathUtils::DeltaR(ele->eta, ele->phi, ZMuon2->eta, ZMuon2->phi) < 0.3) continue;

        Double_t leadingJetPt = -1;
        Double_t leptonJetPt = -1;
        //pass event selection     
        for(Int_t j=0; j<jetArr->GetEntries(); j++) {
          const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[j]);        
          if (jet->pt > leadingJetPt &&
              mithep::MathUtils::DeltaR(jet->phi, jet->eta, ele->phi, ele->eta) > 1.0) {
            leadingJetPt = jet->pt;          
          }
          if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, ele->phi, ele->eta) < 0.5) {
            leptonJetPt = jet->pt;
          }
        }
      
        //if there's no jet ( == isolated lepton?) then take pt of lepton
        if (leptonJetPt < 0) {
          leptonJetPt = ele->pt;
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
              if ((ele->pt < 5)) continue;
              if ((ele->pt > 35)) continue;

              if (passElectronDenominatorCuts(ele, denominatorType[denominatorTypeIndex])) {

                if (ele->pt < 10 && fabs(ele->scEta) < 0.8) {
                  eventListFile << info->runNum << " " 
                                << info->lumiSec << " " 
                                << info->evtNum<< " "
                                << " : "
                                << ele->pt << " " << ele->eta << " " << ele->phi << " "
                                << PassCiCID(ele,1)
                                << endl;
                }

                DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt,eventweight);
                DenominatorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(ele->eta),eventweight);
                DenominatorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->phi,eventweight);
                DenominatorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                DenominatorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                DenominatorVector_EtaPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(ele->eta),ele->pt, eventweight);
              
                if (ele->pt < 20 && fabs(ele->eta) < 1.479) {
                  DenominatorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                } 
                if (ele->pt < 20 && fabs(ele->eta) >= 1.479) {
                  DenominatorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                } 
                if (ele->pt >= 20 && fabs(ele->eta) < 1.479) {
                  DenominatorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                } 
                if (ele->pt >= 20 && fabs(ele->eta) >= 1.479) {
                  DenominatorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                } 
              

                //*************************************************
                //Compute Isolation
                //*************************************************
                Bool_t passIsoMVALoose = kFALSE;
                Bool_t passIsoMVATight = kFALSE;
                Bool_t passIsoMVASameBkgAsDetIso025 = kFALSE;
                Bool_t passPFIsoSameBkgAsDetIso025 = kFALSE;

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
                  if (pf->matchedObjectType == 11 && pf->matchedObjectIndex == i) continue;

                  Double_t deta = (ele->eta - pf->eta);
                  Double_t dphi = mithep::MathUtils::DeltaPhi(Double_t(ele->phi),Double_t(pf->phi));
                  Double_t dr = mithep::MathUtils::DeltaR(ele->phi,ele->eta, pf->phi, pf->eta);
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
                      if (fabs(ele->scEta) > 1.479 && dr < 0.015) continue;
                      //************************************************************

                      if (dr < 0.1) tmpChargedIso_DR0p0To0p1 += pf->pt;                  
                      if (dr >= 0.1 && dr < 0.2) tmpChargedIso_DR0p1To0p2 += pf->pt;
                      if (dr >= 0.2 && dr < 0.3) tmpChargedIso_DR0p2To0p3 += pf->pt;
                      if (dr >= 0.3 && dr < 0.4) tmpChargedIso_DR0p3To0p4 += pf->pt;
                      if (dr >= 0.4 && dr < 0.5) tmpChargedIso_DR0p4To0p5 += pf->pt;
                      if (dr >= 0.5 && dr < 0.7) tmpChargedIso_DR0p5To0p7 += pf->pt;
                    }
                    else if (pf->pfType == eGamma) {
                      //************************************************************
                      // Footprint Veto
                      if (fabs(ele->scEta) > 1.479) {
                        if (mithep::MathUtils::DeltaR(ele->phi,ele->eta, pf->phi, pf->eta) < 0.08) continue;
                      }
                      //************************************************************
                      if (dr < 0.1) tmpGammaIso_DR0p0To0p1 += pf->pt;                  
                      if (dr >= 0.1 && dr < 0.2) tmpGammaIso_DR0p1To0p2 += pf->pt;
                      if (dr >= 0.2 && dr < 0.3) tmpGammaIso_DR0p2To0p3 += pf->pt;
                      if (dr >= 0.3 && dr < 0.4) tmpGammaIso_DR0p3To0p4 += pf->pt;
                      if (dr >= 0.4 && dr < 0.5) tmpGammaIso_DR0p4To0p5 += pf->pt;
                      if (dr >= 0.5 && dr < 0.7) tmpGammaIso_DR0p5To0p7 += pf->pt; 
                    }
                    else {
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

                varPt = ele->pt;
                varEta = ele->scEta;
                varChargedIso_DR0p0To0p1   = TMath::Min((tmpChargedIso_DR0p0To0p1)/ele->pt, 2.5);
                varChargedIso_DR0p1To0p2   = TMath::Min((tmpChargedIso_DR0p1To0p2)/ele->pt, 2.5);
                varChargedIso_DR0p2To0p3 = TMath::Min((tmpChargedIso_DR0p2To0p3)/ele->pt, 2.5);
                varChargedIso_DR0p3To0p4 = TMath::Min((tmpChargedIso_DR0p3To0p4)/ele->pt, 2.5);
                varChargedIso_DR0p4To0p5 = TMath::Min((tmpChargedIso_DR0p4To0p5)/ele->pt, 2.5);
                varChargedIso_DR0p5To0p7 = TMath::Min((tmpChargedIso_DR0p5To0p7)/ele->pt, 2.5);
                varGammaIso_DR0p0To0p1   = TMath::Max(TMath::Min((tmpGammaIso_DR0p0To0p1 - rho*ElectronEffectiveArea(kEleGammaIsoDR0p0To0p1, ele->scEta, EffectiveAreaVersion))/ele->pt, 2.5), 0.0);
                varGammaIso_DR0p1To0p2   = TMath::Max(TMath::Min((tmpGammaIso_DR0p1To0p2 - rho*ElectronEffectiveArea(kEleGammaIsoDR0p1To0p2, ele->scEta, EffectiveAreaVersion))/ele->pt, 2.5), 0.0);
                varGammaIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpGammaIso_DR0p2To0p3 - rho*ElectronEffectiveArea(kEleGammaIsoDR0p2To0p3, ele->scEta, EffectiveAreaVersion))/ele->pt, 2.5), 0.0);
                varGammaIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpGammaIso_DR0p3To0p4 - rho*ElectronEffectiveArea(kEleGammaIsoDR0p3To0p4, ele->scEta, EffectiveAreaVersion))/ele->pt, 2.5), 0.0);
                varGammaIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpGammaIso_DR0p4To0p5 - rho*ElectronEffectiveArea(kEleGammaIsoDR0p4To0p5, ele->scEta, EffectiveAreaVersion))/ele->pt, 2.5), 0.0);
                varGammaIso_DR0p5To0p7 = TMath::Max(TMath::Min((tmpGammaIso_DR0p5To0p7 - rho*ElectronEffectiveArea(kEleGammaIsoDR0p5To0p7, ele->scEta, EffectiveAreaVersion))/ele->pt, 2.5), 0.0);
                varNeutralHadronIso_DR0p0To0p1   = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p0To0p1 - rho*ElectronEffectiveArea(kEleNeutralHadronIsoDR0p0To0p1, ele->scEta, EffectiveAreaVersion))/ele->pt, 2.5), 0.0);
                varNeutralHadronIso_DR0p1To0p2   = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p1To0p2 - rho*ElectronEffectiveArea(kEleNeutralHadronIsoDR0p1To0p2, ele->scEta, EffectiveAreaVersion))/ele->pt, 2.5), 0.0);
                varNeutralHadronIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p2To0p3 - rho*ElectronEffectiveArea(kEleNeutralHadronIsoDR0p2To0p3, ele->scEta, EffectiveAreaVersion))/ele->pt, 2.5), 0.0);
                varNeutralHadronIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p3To0p4 - rho*ElectronEffectiveArea(kEleNeutralHadronIsoDR0p3To0p4, ele->scEta, EffectiveAreaVersion))/ele->pt, 2.5), 0.0);
                varNeutralHadronIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p4To0p5 - rho*ElectronEffectiveArea(kEleNeutralHadronIsoDR0p4To0p5, ele->scEta, EffectiveAreaVersion))/ele->pt, 2.5), 0.0);
                varNeutralHadronIso_DR0p5To0p7 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p5To0p7 - rho*ElectronEffectiveArea(kEleNeutralHadronIsoDR0p5To0p7, ele->scEta, EffectiveAreaVersion))/ele->pt, 2.5), 0.0);

                //evaluate MVA
                Int_t mvabin = -1;
                if (varPt <= 10 && fabs(varEta) < 1.479) mvabin = 0;
                else if (varPt <= 10 && fabs(varEta) >= 1.479) mvabin = 1;
                else if (varPt > 10 && fabs(varEta) < 1.479) mvabin = 2;
                else if (varPt > 10 && fabs(varEta) >= 1.479) mvabin = 3;
                else { cout << "error\n"; assert(0);}
                Double_t EleIsoMVA = electronIsoMVA->MVAValue_IsoRings(
                  varPt,
                  varEta,
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
                //Double_t EleIsoMVA = ElectronIsoReader[mvabin]->EvaluateMVA( "ElectronIsoMVA_BDTG_V0 method" );
                if (passEleIsoMVA_LooseWP(ele->pt, ele->scEta, EleIsoMVA) 
                    && (varChargedIso_DR0p0To0p1 + varChargedIso_DR0p1To0p2 + varChargedIso_DR0p2To0p3 < 0.7)
                  ) {
                  passIsoMVALoose = kTRUE;
                }
                if (passEleIsoMVA_TightWP(ele->pt, ele->scEta, EleIsoMVA)
                    && (varChargedIso_DR0p0To0p1 + varChargedIso_DR0p1To0p2 + varChargedIso_DR0p2To0p3 < 0.7)
                  ) {
                  passIsoMVATight = kTRUE;
                }
                passIsoMVASameBkgAsDetIso025 = passEleIsoMVASameBkgAsDetIso025(ele->pt, ele->scEta, EleIsoMVA);
                passPFIsoSameBkgAsDetIso025 = passElePFIsoSameBkgAsDetIso025(ele->pt, ele->scEta, varChargedIso_DR0p0To0p1 + varChargedIso_DR0p1To0p2 + varChargedIso_DR0p2To0p3 + varGammaIso_DR0p0To0p1 + varGammaIso_DR0p1To0p2 + varGammaIso_DR0p2To0p3 + varNeutralHadronIso_DR0p0To0p1 + varNeutralHadronIso_DR0p1To0p2+ varNeutralHadronIso_DR0p2To0p3);
                

                //*************************************************
                //Compute ID MVA
                //*************************************************
                Bool_t passIDMVALoose = kFALSE;
                Bool_t passIDMVATight = kFALSE;
                Bool_t passIDMVASameBkgAsCiCTight = kFALSE;

                varfbrem = max(double(ele->fBrem),-1.0);
                vardeta = ele->deltaEtaIn;
                vardphi = ele->deltaPhiIn;
                varsee = ele->sigiEtaiEta;
                varetawidth = ele->SCEtaWidth;
                varphiwidth = ele->SCPhiWidth;
                varHoE = ele->HoverE;
                varEoP =  min(double(ele->EOverP), 20.0);
                vare1x5e5x5 = ele->SeedE1x5OverE/ele->SeedE5x5OverE;
                varEoPout = min(double(ele->ESeedClusterOverPout),20.0);
                vardetacalo = ele->dEtaCalo ;
                varkfchi2 = min(double(ele->KFTrackChi2OverNdof),10.0);
                varkfhits = ele->KFTrackNHits ;
                varspp = TMath::Sqrt(ele->sigiPhiiPhi) ;
                varIoEmIoP = (1.0/(ele->scEt * TMath::CosH(ele->scEta)) - 1/ele->p);
                varR9 = min(double(ele->R9), 5.0) ;
                vargsfchi2 = min(double(ele->GsfTrackChi2OverNdof),200.0);
                varPreShowerOverRaw = ele->PreShowerOverRaw;
                varkfhitsall = ele->KFTrackNHits;
                vareleEoPout = ele->ESeedClusterOverPout ;    
                varOneMinuse1x5e5x5 = min(max(1 - double(ele->SeedE1x5OverE/ele->SeedE5x5OverE) , -1.0),2.0);
                varIoEmIoP_Emanuele = (1 - ele->EOverP)/(ele->scEt * TMath::CosH(ele->scEta));
                vardetaabs = min(fabs(double(ele->deltaEtaIn)),0.06);
                vardphiabs = min(fabs(double(ele->deltaPhiIn)),0.6);
                vardetacaloabs = min(fabs(double(ele->dEtaCalo)),0.2);

                //evaluate MVA
                Int_t idmvabin = -1;
                if (varPt <= 10 && fabs(varEta) < 0.8) idmvabin = 0;
                else if (varPt <= 10 && fabs(varEta) >= 0.8 && fabs(varEta) < 1.479) idmvabin = 1;
                else if (varPt <= 10 && fabs(varEta) >= 1.479) idmvabin = 2;
                else if (varPt > 10 && fabs(varEta) < 0.8) idmvabin = 3;
                else if (varPt > 10 && fabs(varEta) >= 0.8 && fabs(varEta) < 1.479) idmvabin = 4;
                else if (varPt > 10 && fabs(varEta) >= 1.479) idmvabin = 5;
                Double_t EleIDMVA = electronIDMVA->MVAValue_IDNonTrig(
                  varPt,
                  varEta,
                  varfbrem,
                  varkfchi2,
                  varkfhitsall,
                  vargsfchi2,
                  vardetaabs, 
                  vardphiabs,
                  vardetacaloabs,
                  varsee,
                  varspp,
                  varetawidth,
                  varphiwidth,
                  varOneMinuse1x5e5x5,
                  varR9,
                  varHoE,
                  varEoP,
                  varIoEmIoP_Emanuele,
                  varEoPout,
                  varPreShowerOverRaw,
                  kFALSE);                
                //Double_t EleIDMVA = ElectronIDReader[idmvabin]->EvaluateMVA( "ElectronIDMVA_BDTG_V1 method" );
                if (passEleIDMVA_LooseWP(ele->pt, ele->scEta, EleIDMVA)) {
                  passIDMVALoose = kTRUE;
                }              
                if (passEleIDMVA_TightWP(ele->pt, ele->scEta, EleIDMVA)) {
                  passIDMVATight = kTRUE;
                }
                passIDMVASameBkgAsCiCTight = passEleIDMVASameBkgAsCiCTight(ele->pt, ele->scEta, EleIDMVA);
         
   
                if (
                  (Option == 1 && (PassCiCID(ele,1) && PassEleDetIso025(ele,rho) 
                                   && fabs(ele->ip3dSig) < 4.0 && ele->nExpHitsInner <= 1
                    ))
                  ||
                  (Option == 3 && passIDMVASameBkgAsCiCTight && PassEleDetIso025(ele,rho) && fabs(ele->ip3dSig) < 4.0 && ele->nExpHitsInner <= 1)
                  ||
                  (Option == 4 && PassCiCID(ele,1) && passPFIsoSameBkgAsDetIso025 && fabs(ele->ip3dSig) < 4.0 && ele->nExpHitsInner <= 1)
                  ||
                  (Option == 5 && PassCiCID(ele,1) && passIsoMVASameBkgAsDetIso025 && fabs(ele->ip3dSig) < 4.0 && ele->nExpHitsInner <= 1)
                  ||
                  (Option == 10 && passIDMVALoose && passIsoMVALoose && fabs(ele->ip3dSig) < 4.0 && ele->nExpHitsInner <= 1)
                  ||
                  (Option == 11 && passIDMVATight && passIsoMVATight && fabs(ele->ip3dSig) < 4.0 && ele->nExpHitsInner <= 1)
                  ) {
                  NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt, eventweight);
                  NumeratorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(ele->eta), eventweight);
                  NumeratorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->phi, eventweight);
                  NumeratorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  NumeratorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  NumeratorVector_EtaPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(ele->eta), ele->pt , eventweight);   
        
                  if (ele->pt < 20 && fabs(ele->eta) < 1.479) {
                    NumeratorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  }
                  if (ele->pt < 20 && fabs(ele->eta) >= 1.479) {
                    NumeratorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  } 
                  if (ele->pt >= 20 && fabs(ele->eta) < 1.479) {
                    NumeratorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  } 
                  if (ele->pt >= 20 && fabs(ele->eta) >= 1.479) {
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




Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele, Int_t DenominatorType) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(fabs(ele->eta) < 2.5)) pass = kFALSE;


  if (DenominatorType == 1) {
    //Barrel 
    if(!( fabs(ele->dz) < 0.1   
         )
      ) {
      pass = kFALSE;
    }
  }

  if (DenominatorType == 2) {
    //Barrel 
    if(!( ele->ip3dSig < 4.0
          && ele->trkIso03/ele->pt < 0.7
          && fabs(ele->dz) < 0.1   
          && ele->nExpHitsInner <= 1
         )
      ) {
      pass = kFALSE;
    }
  }
  
  return pass;
}


