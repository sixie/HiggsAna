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
#include "HiggsAna/Ntupler/interface/HiggsAnaDefs.hh"
#include "HiggsAna/Ntupler/interface/TEventInfo.hh"
#include "HiggsAna/Ntupler/interface/TElectron.hh"
#include "HiggsAna/Ntupler/interface/TMuon.hh"
#include "HiggsAna/Ntupler/interface/TJet.hh"
#include "HiggsAna/Ntupler/interface/TGenParticle.hh"
#include "HiggsAna/Ntupler/interface/TPFCandidate.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodSwitches.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodMeasurements.h"
#include "HiggsAna/HZZ4l/Utils/LeptonSelection.hh"
#include "HiggsAna/Utils/LeptonIDCuts.hh"

// output data structs
#include "HiggsAna/HZZ4l/interface/HZZKinematics.hh"
#include "HiggsAna/HZZ4l/interface/HZZGenInfo.hh"
#include "HiggsAna/HZZ4l/interface/HZZ4lDefs.hh"

#include "TMVAGui.C"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

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

void HZZ4lSelection_MVAIsoStudy(const string Label = "", Int_t Option = 0) 
{  
  gBenchmark->Start("HZZTemplate");
  string label = Label;
  if (Label != "") label = "_" + Label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  Bool_t printDebug = kFALSE;
  Double_t lumi = 1092;              // luminosity (pb^-1)


  vector<vector<string> > inputFiles;
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/data/blue/sixie/ntuples/HZZ4l/data/HZZNtuple_r11a-data.1092ipb.FourRecoLeptonSkim.root");

  inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-h115zz4l-gf-v14b-pu_noskim_0000.root");
   inputFiles.back().push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-h120zz4l-gf-v14b-pu_noskim_0000.root");
//  inputFiles.back().push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-h130zz4l-gf-v14b-pu_noskim_0000.root");
  //weight = 16.63*0.000169 / (301479) = 9.32227452e-9 / pb^-1
  //weight(4.6fb^-1) = 4.29e-5

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
  Double_t NGenAccepted_4E = 0;
  Double_t NGenAccepted_4M = 0;
  Double_t NGenAccepted_2E2M = 0;

  Double_t NBaselineID_4E = 0;
  Double_t NBaselineID_4M = 0;
  Double_t NBaselineID_2E2M = 0;

  Double_t NBaselinePairwise_4E = 0;
  Double_t NBaselinePairwise_4M = 0;
  Double_t NBaselinePairwise_2E2M = 0;

  Double_t NDetIso025_4E = 0;
  Double_t NDetIso025_4M = 0;
  Double_t NDetIso025_2E2M = 0;

  Double_t NMVAIsoLoose_4E = 0;
  Double_t NMVAIsoLoose_4M = 0;
  Double_t NMVAIsoLoose_2E2M = 0;

  Double_t NMVAIsoTight_4E = 0;
  Double_t NMVAIsoTight_4M = 0;
  Double_t NMVAIsoTight_2E2M = 0;

  //Selected
  Double_t NSelectedBaselineID_4E = 0;
  Double_t NSelectedBaselineID_4M = 0;
  Double_t NSelectedBaselineID_2E2M = 0;

  Double_t NSelectedBaselinePairwise_4E = 0;
  Double_t NSelectedBaselinePairwise_4M = 0;
  Double_t NSelectedBaselinePairwise_2E2M = 0;

  Double_t NSelectedDetIso025_4E = 0;
  Double_t NSelectedDetIso025_4M = 0;
  Double_t NSelectedDetIso025_2E2M = 0;

  Double_t NSelectedMVAIsoLoose_4E = 0;
  Double_t NSelectedMVAIsoLoose_4M = 0;
  Double_t NSelectedMVAIsoLoose_2E2M = 0;

  Double_t NSelectedMVAIsoTight_4E = 0;
  Double_t NSelectedMVAIsoTight_4M = 0;
  Double_t NSelectedMVAIsoTight_2E2M = 0;


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  vector <TH1F*>  fHZZ4lSelection; 
  vector <TH1F*>  fHZZTo4ESelection; 
  vector <TH1F*>  fHZZTo4MuSelection; 
  vector <TH1F*>  fHZZTo2E2MuSelection; 
  vector <TH1F*>  fZ1Mass; 
  vector <TH1F*>  fZ2Mass; 
  vector <TH1F*>  fZZMass; 
  for (int q=0; q<processNames.size() ; ++q) {
    TH1F *tmpHZZ4lSelection = new TH1F(("hHZZ4lSelection"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHZZTo4ESelection = new TH1F(("hHZZTo4ESelection"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHZZTo4MuSelection = new TH1F(("hHZZTo4MuSelection"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHZZTo2E2MuSelection = new TH1F(("hHZZTo2E2MuSelection"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpZ1Mass = new TH1F(("hZ1Mass"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 40, 0, 200);
    TH1F *tmpZ2Mass = new TH1F(("hZ2Mass"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 40, 0, 200);
    TH1F *tmpZZMass = new TH1F(("hZZMass"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 40, 0, 600);
    tmpHZZ4lSelection->Sumw2();
    tmpHZZTo4ESelection->Sumw2();
    tmpHZZTo4MuSelection->Sumw2();
    tmpHZZTo2E2MuSelection->Sumw2();   
    tmpZ1Mass->Sumw2();   
    tmpZ2Mass->Sumw2();   
    tmpZZMass->Sumw2();   
    
    fHZZ4lSelection.push_back(tmpHZZ4lSelection);
    fHZZTo4ESelection.push_back(tmpHZZTo4ESelection);
    fHZZTo4MuSelection.push_back(tmpHZZTo4MuSelection);
    fHZZTo2E2MuSelection.push_back(tmpHZZTo2E2MuSelection);
    fZ1Mass.push_back(tmpZ1Mass);
    fZ2Mass.push_back(tmpZ2Mass);
    fZZMass.push_back(tmpZZMass);
  }


  vector<string> CutLabel;
  CutLabel.push_back("Dilepton20/10");

  
  //--------------------------------------------------------------------------------------------------------------
  // output ntuple structure
  //==============================================================================================================  
  HZZKinematics kinematics;
  HZZGenInfo    geninfo;

  TFile *fOutputFile = new TFile(("HZZOutput"+label+".root").c_str(), "RECREATE");
  TTree *fOutputTree = new TTree("hzz4l","hzz4l");
  fOutputTree->Branch("geninfo",&geninfo);
  fOutputTree->Branch("kinematics",&kinematics);



  //--------------------------------------------------------------------------------------------------------------
  // Pileup Reweighting
  //==============================================================================================================  
  TFile *fPUFile = TFile::Open("/data/smurf/sixie/Pileup/weights/PileupReweighting.Fall11_To_Full2011.root");
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;
  Int_t EffectiveAreaVersion = -1;  
  EffectiveAreaVersion = 2;  //Fall11MC


  //--------------------------------------------------------------------------------------------------------------
  // Load IsoMVAs
  //==============================================================================================================  
  TMVA::Reader *ElectronIsoReader = new TMVA::Reader( "!Color:!Silent" );    
  TMVA::Reader *MuonIsoReader = new TMVA::Reader( "!Color:!Silent" );    

  Float_t                 varChargedIso_DR0p0To0p2;
  Float_t                 varGammaIso_DR0p0To0p2;
  Float_t                 varNeutralHadronIso_DR0p0To0p2;


  Float_t                 varChargedIso_DR0p0To0p1;
  Float_t                 varChargedIso_DR0p1To0p2;
  Float_t                 varChargedIso_DR0p2To0p3;
  Float_t                 varChargedIso_DR0p3To0p4;
  Float_t                 varChargedIso_DR0p4To0p5;
  Float_t                 varChargedIso_DR0p5To0p7;
  Float_t                 varGammaIso_DR0p0To0p1;
  Float_t                 varGammaIso_DR0p1To0p2;
  Float_t                 varGammaIso_DR0p2To0p3;
  Float_t                 varGammaIso_DR0p3To0p4;
  Float_t                 varGammaIso_DR0p4To0p5;
  Float_t                 varGammaIso_DR0p5To0p7;
  Float_t                 varNeutralHadronIso_DR0p0To0p1;
  Float_t                 varNeutralHadronIso_DR0p1To0p2;
  Float_t                 varNeutralHadronIso_DR0p2To0p3;
  Float_t                 varNeutralHadronIso_DR0p3To0p4;
  Float_t                 varNeutralHadronIso_DR0p4To0p5;
  Float_t                 varNeutralHadronIso_DR0p5To0p7;

  Float_t                 varPFIso_MeanEta;
  Float_t                 varPFIso_MeanPhi;
  Float_t                 varPFIso_SigmaEtaEta;
  Float_t                 varPFIso_SigmaPhiPhi;
  Float_t                 varPFIso_SigmaEtaPhi;
  Float_t                 varChargedIso_MeanEta;
  Float_t                 varChargedIso_MeanPhi;
  Float_t                 varChargedIso_SigmaEtaEta;
  Float_t                 varChargedIso_SigmaPhiPhi;
  Float_t                 varChargedIso_SigmaEtaPhi;
  Float_t                 varGammaIso_MeanEta;
  Float_t                 varGammaIso_MeanPhi;
  Float_t                 varGammaIso_SigmaEtaEta;
  Float_t                 varGammaIso_SigmaPhiPhi;
  Float_t                 varGammaIso_SigmaEtaPhi;
  Float_t                 varNeutralHadronIso_MeanEta;
  Float_t                 varNeutralHadronIso_MeanPhi;
  Float_t                 varNeutralHadronIso_SigmaEtaEta;
  Float_t                 varNeutralHadronIso_SigmaPhiPhi;
  Float_t                 varNeutralHadronIso_SigmaEtaPhi;
  Float_t                 varDirectionalPFIso;
  Float_t                 varDirectionalChargedIso;
  Float_t                 varDirectionalGammaIso;
  Float_t                 varDirectionalNeutralHadronIso;
  

  ElectronIsoReader->AddVariable( "ChargedIso_DR0p0To0p1",         &varChargedIso_DR0p0To0p1           );
  ElectronIsoReader->AddVariable( "ChargedIso_DR0p1To0p2",         &varChargedIso_DR0p1To0p2           );
  ElectronIsoReader->AddVariable( "ChargedIso_DR0p2To0p3",       &varChargedIso_DR0p2To0p3         );
  ElectronIsoReader->AddVariable( "ChargedIso_DR0p3To0p4",        &varChargedIso_DR0p3To0p4          );
  ElectronIsoReader->AddVariable( "ChargedIso_DR0p4To0p5",        &varChargedIso_DR0p4To0p5          );
  ElectronIsoReader->AddVariable( "ChargedIso_DR0p5To0p7",        &varChargedIso_DR0p5To0p7          );
  ElectronIsoReader->AddVariable( "GammaIso_DR0p0To0p1",         &varGammaIso_DR0p0To0p1           );
  ElectronIsoReader->AddVariable( "GammaIso_DR0p1To0p2",         &varGammaIso_DR0p1To0p2           );
  ElectronIsoReader->AddVariable( "GammaIso_DR0p2To0p3",       &varGammaIso_DR0p2To0p3         );
  ElectronIsoReader->AddVariable( "GammaIso_DR0p3To0p4",        &varGammaIso_DR0p3To0p4          );
  ElectronIsoReader->AddVariable( "GammaIso_DR0p4To0p5",          &varGammaIso_DR0p4To0p5          );
  ElectronIsoReader->AddVariable( "GammaIso_DR0p5To0p7",          &varGammaIso_DR0p5To0p7          );
  ElectronIsoReader->AddVariable( "NeutralHadronIso_DR0p0To0p1",         &varNeutralHadronIso_DR0p0To0p1           );
  ElectronIsoReader->AddVariable( "NeutralHadronIso_DR0p1To0p2",         &varNeutralHadronIso_DR0p1To0p2           );
  ElectronIsoReader->AddVariable( "NeutralHadronIso_DR0p2To0p3",       &varNeutralHadronIso_DR0p2To0p3         );
  ElectronIsoReader->AddVariable( "NeutralHadronIso_DR0p3To0p4",        &varNeutralHadronIso_DR0p3To0p4          );
  ElectronIsoReader->AddVariable( "NeutralHadronIso_DR0p4To0p5",        &varNeutralHadronIso_DR0p4To0p5          );
  ElectronIsoReader->AddVariable( "NeutralHadronIso_DR0p5To0p7",        &varNeutralHadronIso_DR0p5To0p7          );

//   MuonIsoReader->AddVariable( "ChargedIso_DR0To0p2",         &varChargedIso_DR0p0To0p2           );
//   MuonIsoReader->AddVariable( "ChargedIso_DR0p2To0p3",       &varChargedIso_DR0p2To0p3         );
//   MuonIsoReader->AddVariable( "ChargedIso_DR0p3To0p4",        &varChargedIso_DR0p3To0p4          );
//   MuonIsoReader->AddVariable( "ChargedIso_DR0p4To0p5",        &varChargedIso_DR0p4To0p5          );
//   MuonIsoReader->AddVariable( "ChargedIso_DR0p5To0p7",        &varChargedIso_DR0p5To0p7          );
//   MuonIsoReader->AddVariable( "GammaIso_DR0To0p2",         &varGammaIso_DR0p0To0p2           );
//   MuonIsoReader->AddVariable( "GammaIso_DR0p2To0p3",       &varGammaIso_DR0p2To0p3         );
//   MuonIsoReader->AddVariable( "GammaIso_DR0p3To0p4",        &varGammaIso_DR0p3To0p4          );
//   MuonIsoReader->AddVariable( "GammaIso_DR0p4To0p5",          &varGammaIso_DR0p4To0p5          );
//   MuonIsoReader->AddVariable( "GammaIso_DR0p5To0p7",          &varGammaIso_DR0p5To0p7          );
//   MuonIsoReader->AddVariable( "NeutralHadronIso_DR0To0p2",         &varNeutralHadronIso_DR0p0To0p2           );
//   MuonIsoReader->AddVariable( "NeutralHadronIso_DR0p2To0p3",       &varNeutralHadronIso_DR0p2To0p3         );
//   MuonIsoReader->AddVariable( "NeutralHadronIso_DR0p3To0p4",        &varNeutralHadronIso_DR0p3To0p4          );
//   MuonIsoReader->AddVariable( "NeutralHadronIso_DR0p4To0p5",        &varNeutralHadronIso_DR0p4To0p5          );
//   MuonIsoReader->AddVariable( "NeutralHadronIso_DR0p5To0p7",        &varNeutralHadronIso_DR0p5To0p7          );

  //We have to use all of these vars because the weight file had it
  MuonIsoReader->AddVariable( "PFIso_MeanEta",                      &varPFIso_MeanEta                        );
  MuonIsoReader->AddVariable( "PFIso_MeanPhi",                      &varPFIso_MeanPhi                        );
  MuonIsoReader->AddVariable( "PFIso_SigmaEtaEta",                  &varPFIso_SigmaEtaEta                    );
  MuonIsoReader->AddVariable( "PFIso_SigmaPhiPhi",                  &varPFIso_SigmaPhiPhi                    );
  MuonIsoReader->AddVariable( "PFIso_SigmaEtaPhi",                  &varPFIso_SigmaEtaPhi                    );
  MuonIsoReader->AddVariable( "ChargedIso_DR0To0p2",         &varChargedIso_DR0p0To0p2           );
  MuonIsoReader->AddVariable( "ChargedIso_DR0p2To0p3",       &varChargedIso_DR0p2To0p3         );
  MuonIsoReader->AddVariable( "ChargedIso_DR0p3To0p4",        &varChargedIso_DR0p3To0p4          );
  MuonIsoReader->AddVariable( "ChargedIso_DR0p4To0p5",        &varChargedIso_DR0p4To0p5          );
  MuonIsoReader->AddVariable( "ChargedIso_DR0p5To0p7",        &varChargedIso_DR0p5To0p7          );
  MuonIsoReader->AddVariable( "ChargedIso_MeanEta",                 &varChargedIso_MeanEta                   );
  MuonIsoReader->AddVariable( "ChargedIso_MeanPhi",                 &varChargedIso_MeanPhi                   );
  MuonIsoReader->AddVariable( "ChargedIso_SigmaEtaEta",             &varChargedIso_SigmaEtaEta               );
  MuonIsoReader->AddVariable( "ChargedIso_SigmaPhiPhi",             &varChargedIso_SigmaPhiPhi               );
  MuonIsoReader->AddVariable( "ChargedIso_SigmaEtaPhi",             &varChargedIso_SigmaEtaPhi               );
  MuonIsoReader->AddVariable( "GammaIso_DR0To0p2",         &varGammaIso_DR0p0To0p2           );
  MuonIsoReader->AddVariable( "GammaIso_DR0p2To0p3",       &varGammaIso_DR0p2To0p3         );
  MuonIsoReader->AddVariable( "GammaIso_DR0p3To0p4",        &varGammaIso_DR0p3To0p4          );
  MuonIsoReader->AddVariable( "GammaIso_DR0p4To0p5",        &varGammaIso_DR0p4To0p5          );
  MuonIsoReader->AddVariable( "GammaIso_DR0p5To0p7",        &varGammaIso_DR0p5To0p7          );
  MuonIsoReader->AddVariable( "GammaIso_MeanEta",                 &varGammaIso_MeanEta                   );
  MuonIsoReader->AddVariable( "GammaIso_MeanPhi",                 &varGammaIso_MeanPhi                   );
  MuonIsoReader->AddVariable( "GammaIso_SigmaEtaEta",             &varGammaIso_SigmaEtaEta               );
  MuonIsoReader->AddVariable( "GammaIso_SigmaPhiPhi",             &varGammaIso_SigmaPhiPhi               );
  MuonIsoReader->AddVariable( "GammaIso_SigmaEtaPhi",             &varGammaIso_SigmaEtaPhi               );
  MuonIsoReader->AddVariable( "NeutralHadronIso_DR0To0p2",         &varNeutralHadronIso_DR0p0To0p2           );
  MuonIsoReader->AddVariable( "NeutralHadronIso_DR0p2To0p3",       &varNeutralHadronIso_DR0p2To0p3         );
  MuonIsoReader->AddVariable( "NeutralHadronIso_DR0p3To0p4",        &varNeutralHadronIso_DR0p3To0p4          );
  MuonIsoReader->AddVariable( "NeutralHadronIso_DR0p4To0p5",        &varNeutralHadronIso_DR0p4To0p5          );
  MuonIsoReader->AddVariable( "NeutralHadronIso_DR0p5To0p7",        &varNeutralHadronIso_DR0p5To0p7          );
  MuonIsoReader->AddVariable( "NeutralHadronIso_MeanEta",                 &varNeutralHadronIso_MeanEta                   );
  MuonIsoReader->AddVariable( "NeutralHadronIso_MeanPhi",                 &varNeutralHadronIso_MeanPhi                   );
  MuonIsoReader->AddVariable( "NeutralHadronIso_SigmaEtaEta",             &varNeutralHadronIso_SigmaEtaEta               );
  MuonIsoReader->AddVariable( "NeutralHadronIso_SigmaPhiPhi",             &varNeutralHadronIso_SigmaPhiPhi               );
  MuonIsoReader->AddVariable( "NeutralHadronIso_SigmaEtaPhi",             &varNeutralHadronIso_SigmaEtaPhi               );
  MuonIsoReader->AddVariable( "DirectionalPFIso",                   &varDirectionalPFIso                     );
  MuonIsoReader->AddVariable( "DirectionalChargedIso",              &varDirectionalChargedIso                );
  MuonIsoReader->AddVariable( "DirectionalGammaIso",                &varDirectionalGammaIso                  );
  MuonIsoReader->AddVariable( "DirectionalNeutralHadronIso",        &varDirectionalNeutralHadronIso          );

  Float_t varPt, varEta;
  ElectronIsoReader->AddSpectator( "eta",   &varEta );
  ElectronIsoReader->AddSpectator( "pt" ,   &varPt  );
  MuonIsoReader->AddSpectator( "eta",   &varEta );
  MuonIsoReader->AddSpectator( "pt" ,   &varPt  );

  ElectronIsoReader->BookMVA( "ElectronIsoMVA_BDTG_V0 method", "/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/weights/ElectronIsoMVA_BDTCat_BDTG_V0.weights.xml");   
  MuonIsoReader->BookMVA( "MuIsoMVA_BDTG_V1 method", "/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IsoMVA/weights_Save/MuonIsoMVA_BDTCat_BDTG_V1.weights.xml");   


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TFile *inputFile=0;
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  TClonesArray *genparticleArr = new TClonesArray("mithep::TGenParticle");
  TClonesArray *pfcandidateArr = new TClonesArray("mithep::TPFCandidate");
  
  for (int q = 0; q<inputFiles.size() ; ++q) { 
    for (int f = 0; f < inputFiles[q].size() ; ++f) {
      //********************************************************
      // Get Tree
      //********************************************************
      cout << "Reading File " << inputFiles[q][f] << endl;
      eventTree = getTreeFromFile(inputFiles[q][f].c_str(),"Events"); 
      TBranch *infoBr;
      TBranch *electronBr;
      TBranch *muonBr;
      TBranch *jetBr;
      TBranch *genparticleBr;
      TBranch *pfcandidateBr;


      //*****************************************************************************************
      //Loop over muon Data Tree
      //*****************************************************************************************
      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");
      eventTree->SetBranchAddress("GenParticle", &genparticleArr);         genparticleBr = eventTree->GetBranch("GenParticle");
      eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr);         pfcandidateBr = eventTree->GetBranch("PFCandidate");
  
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
        infoBr->GetEntry(ientry);
        if (ientry % 10000 == 0) cout << "Event " << ientry << endl;
	
        //Use Only Testing events
        if (info->evtNum % 2 == 0 ) continue;
        //Use Only Training events
        //if (info->evtNum % 2 != 0 ) continue;
       

        //********************************************************
        double mynpu = TMath::Min((double)info->nPUEvents,34.999);
        Int_t npuxbin = fhDPU->GetXaxis()->FindBin(mynpu);
        double npuWeight = fhDPU->GetBinContent(npuxbin);
        //********************************************************
        //********************************************************
        // Pileup Energy Density
        //********************************************************
        Double_t rho = 0;
        if (!(TMath::IsNaN(info->PileupEnergyDensity) || isinf(info->PileupEnergyDensity))) rho = info->PileupEnergyDensity;
        

        //********************************************************
        // Printdebug
        //********************************************************
        printDebug = kFALSE;
        if ((0 == 1) 
            || (info->evtNum == 9316)           
          ) printDebug = kTRUE;
        



        //********************************************************
        // Load the branches
        //********************************************************
        electronArr->Clear(); 
        muonArr->Clear(); 
        jetArr->Clear();  //if (info->evtNum % 2 != 0 ) continue;
       
        genparticleArr->Clear(); 
        pfcandidateArr->Clear(); 
        electronBr->GetEntry(ientry);
        muonBr->GetEntry(ientry);
        jetBr->GetEntry(ientry);
        genparticleBr->GetEntry(ientry);
        pfcandidateBr->GetEntry(ientry);

//         Double_t eventweight = info->eventweight * lumi * npuWeight;
        Double_t eventweight = npuWeight;


        //********************************************************
        // GenInfo
        //********************************************************
        vector<Double_t> pdgid;
        vector<Double_t> pt;
        vector<Double_t> eta;
        vector<Double_t> phi;
        vector<Double_t> pt_ele;
        vector<Double_t> eta_ele;
        vector<Double_t> phi_ele;
        vector<Double_t> pt_mu;
        vector<Double_t> eta_mu;
        vector<Double_t> phi_mu;

        Double_t minGenLeptonDR = 9999;
        for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
          const mithep::TGenParticle *gen = (mithep::TGenParticle*)((*genparticleArr)[k]);


          if (gen->status == 1 && (abs(gen->pdgid) == 11 || abs(gen->pdgid) == 13 )) {
            if (abs(gen->pdgid) == 11) {
              if (!(fabs(gen->eta) < 2.5 && gen->pt > 7.0)) continue;
            }
            if (abs(gen->pdgid) == 13) {
              if (!(fabs(gen->eta) < 2.4 && gen->pt > 5.0)) continue;
            }
            pdgid.push_back(gen->pdgid); 
            pt.push_back(gen->pt); 
            eta.push_back(gen->eta); 
            phi.push_back(gen->phi); 
            if (abs(gen->pdgid) == 11) {
              pt_ele.push_back(gen->pt);
              eta_ele.push_back(gen->eta);
              phi_ele.push_back(gen->phi);
            }
            if (abs(gen->pdgid) == 13) {
              pt_mu.push_back(gen->pt);
              eta_mu.push_back(gen->eta);
              phi_mu.push_back(gen->phi);
            }
          }

          //find pair of gen leptons that are closest in DR
          if (!(gen->status == 1 && (abs(gen->pdgid) == 11 || abs(gen->pdgid) == 13))) continue;
          for(Int_t l=k+1; l<genparticleArr->GetEntries(); l++) {
            const mithep::TGenParticle *p2 = (mithep::TGenParticle*)((*genparticleArr)[l]);
            if (!(p2->status == 1 && (abs(p2->pdgid) == 11 || abs(p2->pdgid) == 13))) continue;
            if (minGenLeptonDR > mithep::MathUtils::DeltaR(gen->phi,gen->eta,p2->phi,p2->eta)) {
              minGenLeptonDR = mithep::MathUtils::DeltaR(gen->phi,gen->eta,p2->phi,p2->eta);
            }
          }
        }
        //******************************************************************
        //Remove events with two leptons that are very close together
        //******************************************************************
        if (minGenLeptonDR < 0.4) continue;


        geninfo.id_1_a = 0;
        geninfo.id_1_b = 0;
        geninfo.id_2_a = 0;
        geninfo.id_2_b = 0;
        geninfo.pt_1_a = 0;
        geninfo.pt_1_b = 0;
        geninfo.pt_2_a = 0;
        geninfo.pt_2_b = 0;
        geninfo.eta_1_a = 0;
        geninfo.eta_1_b = 0;
        geninfo.eta_2_a = 0;
        geninfo.eta_2_b = 0;

        if (pdgid.size() == 4) {
          if (pt_ele.size() == 4) {
            geninfo.id_1_a = 11;
            geninfo.id_1_b = 11;
            geninfo.id_2_a = 11;
            geninfo.id_2_b = 11;
            geninfo.pt_1_a = pt_ele[0];
            geninfo.pt_1_b = pt_ele[1];
            geninfo.pt_2_a = pt_ele[2];
            geninfo.pt_2_b = pt_ele[3];
            geninfo.eta_1_a = eta_ele[0];
            geninfo.eta_1_b = eta_ele[1];
            geninfo.eta_2_a = eta_ele[2];
            geninfo.eta_2_b = eta_ele[3];
          }
          if (pt_mu.size() == 4) {
            geninfo.id_1_a = 13;
            geninfo.id_1_b = 13;
            geninfo.id_2_a = 13;
            geninfo.id_2_b = 13;
            geninfo.pt_1_a = pt_mu[0];
            geninfo.pt_1_b = pt_mu[1];
            geninfo.pt_2_a = pt_mu[2];
            geninfo.pt_2_b = pt_mu[3];
            geninfo.eta_1_a = eta_mu[0];
            geninfo.eta_1_b = eta_mu[1];
            geninfo.eta_2_a = eta_mu[2];
            geninfo.eta_2_b = eta_mu[3];
          }
          if (pt_ele.size() == 2 && pt_mu.size() == 2) {
            geninfo.id_1_a = 11;
            geninfo.id_1_b = 11;
            geninfo.id_2_a = 13;
            geninfo.id_2_b = 13;
            geninfo.pt_1_a = pt_ele[0];
            geninfo.pt_1_b = pt_ele[1];
            geninfo.pt_2_a = pt_mu[0];
            geninfo.pt_2_b = pt_mu[1];
            geninfo.eta_1_a = eta_ele[0];
            geninfo.eta_1_b = eta_ele[1];
            geninfo.eta_2_a = eta_mu[0];
            geninfo.eta_2_b = eta_mu[1];
          }
        }

        if (pt_ele.size() == 4) {
          NGenAccepted_4E += eventweight;
        } else if (pt_mu.size() == 4) {
          NGenAccepted_4M += eventweight;
        } else if (pt_ele.size() == 2 && pt_mu.size() == 2) {
          NGenAccepted_2E2M += eventweight;
        }

        //********************************************************
        // Met
        //********************************************************
        TVector3 pfMet;        
        if(info->pfMEx!=0 || info->pfMEy!=0) {       
          pfMet.SetXYZ(info->pfMEx, info->pfMEy, 0);
        }

        //********************************************************
        // Lepton Selection
        //********************************************************
        Int_t NLeptons = 0;
        vector<Int_t> leptonType;
        vector<Int_t> leptonIndex;
        vector<Double_t> leptonPt;
        vector<Double_t> leptonEta;
        vector<Double_t> leptonPhi;
        vector<Int_t> leptonCharge;

        Int_t NJets = 0;
        const mithep::TJet *leadingJet = 0;
    
        for(Int_t i=0; i<muonArr->GetEntries(); i++) {
          const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);      
          
          if ( (0==0)
               &&
               mu->pt > 5.0
               &&
               fabs(mu->eta) < 2.4
               && 
               passMuonID_HZZ2011(mu)
               &&
               passMuonIP_HZZ2011(mu)                              
            ) {
            leptonPt.push_back(mu->pt);
            leptonEta.push_back(mu->eta);
            leptonPhi.push_back(mu->phi);
            leptonType.push_back(13);
            leptonIndex.push_back(i);  
            leptonCharge.push_back(mu->q);
          }
        }

        for(Int_t i=0; i<electronArr->GetEntries(); i++) {
          const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);

          //Do we do this?
          Bool_t isMuonOverlap = kFALSE;
          for (int k=0; k<leptonPt.size(); ++k) {
            if ( leptonType[k] == 13 
                 && mithep::MathUtils::DeltaR(ele->phi, ele->eta, leptonPhi[k],leptonEta[k]) < 0.1
              ) {
              isMuonOverlap = kTRUE; 
              break;
            }
          }

          if ( (0==0)
               &&
               ele->pt > 7.0
               && 
               fabs(ele->eta) < 2.5
               &&                                   
               !isMuonOverlap   
               &&     
               PassCiCID(ele,1)
            ) {
            leptonPt.push_back(ele->pt);
            leptonEta.push_back(ele->eta);
            leptonPhi.push_back(ele->phi);
            leptonType.push_back(11);
            leptonIndex.push_back(i);
            leptonCharge.push_back(ele->q);
          }
        }

        //sort leptons
        Int_t tempType;
        Int_t tempIndex;
        Double_t tempPt;
        Double_t tempEta;
        Double_t tempPhi;
        Int_t tempCharge;
        for (int l=0; l<leptonIndex.size(); l++) {
          for (int k=0; k < leptonIndex.size() - 1; k++) {
            if (leptonPt[k+1] > leptonPt[k]) {
              tempType = leptonType[k];
              tempIndex = leptonIndex[k];
              tempPt = leptonPt[k];
              tempEta = leptonEta[k];
              tempPhi = leptonPhi[k];
              tempCharge = leptonCharge[k];
          
              leptonType[k] = leptonType[k+1];
              leptonIndex[k] = leptonIndex[k+1];
              leptonPt[k] = leptonPt[k+1];
              leptonEta[k] = leptonEta[k+1];
              leptonPhi[k] = leptonPhi[k+1];
              leptonCharge[k] = leptonCharge[k+1];

              leptonType[k+1] = tempType;
              leptonIndex[k+1] = tempIndex;
              leptonPt[k+1] = tempPt;
              leptonEta[k+1] = tempEta;
              leptonPhi[k+1] = tempPhi;
              leptonCharge[k+1] = tempCharge;
          
            }
          }
        }

        //remove 5lepton events just to be not confused for the pairwise isolation
        if (leptonPt.size() > 4) continue; 
        if (leptonPt.size() < 4) continue; 

        //******************************************************************************
        //How many events have 4 ID'ed leptons?
        //******************************************************************************
        if (leptonPt.size() >= 4) {
          if (pt_ele.size() == 4) {
            NBaselineID_4E += eventweight;
          } else if (pt_mu.size() == 4) {
            NBaselineID_4M += eventweight;
          } else if (pt_ele.size() == 2 && pt_mu.size() == 2) {
            NBaselineID_2E2M += eventweight;
          }
        }

        //******************************************************************************
        //DetIso/pt < 0.25 Selection
        //******************************************************************************
        Int_t NLeptonsPassDetIso025 = 0;
        for(int i = 0; i < leptonPt.size(); ++i) {
          if (leptonType[i] == 11) {
            const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[leptonIndex[i]]);
            Double_t DetIso03PUCorrection = 0;
            if (fabs(ele->scEta) < 1.5) DetIso03PUCorrection = rho * (0.078 + 0.026);
            if ( (ele->trkIso03 + ele->emIso03 + ele->hadIso03 - DetIso03PUCorrection) / ele->pt < 0.25) NLeptonsPassDetIso025 += eventweight;
          } else if (leptonType[i] == 13) {
            const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[leptonIndex[i]]);   
            Double_t DetIso03PUCorrection = 0;
            if (fabs(mu->eta) < 1.5) DetIso03PUCorrection = rho * (0.087 + 0.042);
            if ( (mu->trkIso03 + mu->emIso03 + mu->hadIso03 - DetIso03PUCorrection) / mu->pt < 0.25) NLeptonsPassDetIso025 += eventweight;
          } else {
            cout << "weird lepton type : " << leptonType[i] << endl;
          }          
        }
        if (NLeptonsPassDetIso025 >= 4) {
          if (pt_ele.size() == 4) {
            NDetIso025_4E += eventweight;
          } else if (pt_mu.size() == 4) {
            NDetIso025_4M += eventweight;
          } else if (pt_ele.size() == 2 && pt_mu.size() == 2) {
            NDetIso025_2E2M += eventweight;
          }
        }

        //******************************************************************************
        //Pairwise Iso Selection
        //******************************************************************************
        Bool_t PassPairwiseIso = kTRUE;
        for(int i = 0; i < leptonPt.size(); ++i) {
          Double_t iso1 = 0;
          if (leptonType[i] == 11) {
            const mithep::TElectron *ele1 = (mithep::TElectron*)((*electronArr)[leptonIndex[i]]);
            Double_t DetIso03PUCorrection = 0;
            if (fabs(ele1->scEta) < 1.5) DetIso03PUCorrection = rho * (0.078 + 0.026);
            iso1 = (ele1->trkIso03 + ele1->emIso03 + ele1->hadIso03 - DetIso03PUCorrection) / ele1->pt;
          } else if (leptonType[i] == 13) {
            const mithep::TMuon *mu1 = (mithep::TMuon*)((*muonArr)[leptonIndex[i]]);   
            Double_t DetIso03PUCorrection = 0;
            if (fabs(mu1->eta) < 1.5) DetIso03PUCorrection = rho * (0.087 + 0.042);
            iso1 = (mu1->trkIso03 + mu1->emIso03 + mu1->hadIso03 - DetIso03PUCorrection) / mu1->pt;
          } else {
            cout << "weird lepton type : " << leptonType[i] << endl;
          }

          for(int j = 0; j < leptonPt.size(); ++j) {
            Double_t iso2 = 0;
            if (leptonType[j] == 11) {
              const mithep::TElectron *ele2 = (mithep::TElectron*)((*electronArr)[leptonIndex[j]]);
              Double_t DetIso03PUCorrection = 0;
              if (fabs(ele2->scEta) < 1.5) DetIso03PUCorrection = rho * (0.078 + 0.026);
              iso2 = (ele2->trkIso03 + ele2->emIso03 + ele2->hadIso03 - DetIso03PUCorrection) / ele2->pt;
            } else if (leptonType[j] == 13) {
              const mithep::TMuon *mu2 = (mithep::TMuon*)((*muonArr)[leptonIndex[j]]);   
              Double_t DetIso03PUCorrection = 0;
              if (fabs(mu2->eta) < 1.5) DetIso03PUCorrection = rho * (0.087 + 0.042);
              iso2 = (mu2->trkIso03 + mu2->emIso03 + mu2->hadIso03 - DetIso03PUCorrection) / mu2->pt;
            } else {
              cout << "weird lepton type : " << leptonType[j] << endl;
            }

            if (iso1 + iso2 > 0.35) PassPairwiseIso = kFALSE;
          }
        }
        if (PassPairwiseIso) {
          if (pt_ele.size() == 4) {
            NBaselinePairwise_4E += eventweight;
          } else if (pt_mu.size() == 4) {
            NBaselinePairwise_4M += eventweight;
          } else if (pt_ele.size() == 2 && pt_mu.size() == 2) {
            NBaselinePairwise_2E2M += eventweight;
          }
        }


        //******************************************************************************
        //MVA Iso Selection
        //******************************************************************************
        Int_t NLeptonsPassMVAIsoLoose = 0;
        Int_t NLeptonsPassMVAIsoTight = 0;
        for(int i = 0; i < leptonPt.size(); ++i) {
          if (leptonType[i] == 11) {
            const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[leptonIndex[i]]);

            
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
              if (pf->matchedObjectType == 11 && pf->matchedObjectIndex == leptonIndex[i]) continue;

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
                for (Int_t q=0; q < leptonType.size() ; ++q) {
                  if (leptonType[q] == 11) {
                    const mithep::TElectron *tmpele = (mithep::TElectron*)((*electronArr)[leptonIndex[q]]);
                    if (leptonType[q] == 11 && pf->matchedObjectType == 11 && pf->matchedObjectIndex == leptonIndex[q]) {
                      IsLeptonFootprint = kTRUE;
                    }
                    if (pf->q != 0 && fabs(tmpele->scEta) > 1.479 && mithep::MathUtils::DeltaR(tmpele->phi,tmpele->eta, pf->phi, pf->eta) < 0.015) IsLeptonFootprint = kTRUE;
                    if (pf->pfType == eGamma) {
                      if (fabs(tmpele->scEta) > 1.479) {
                        if (mithep::MathUtils::DeltaR(tmpele->phi,tmpele->eta, pf->phi, pf->eta) < 0.08) IsLeptonFootprint = kTRUE;
                      }
                    }
                  }
                  if (leptonType[q] == 13) {
                    const mithep::TMuon *tmpmu = (mithep::TMuon*)((*muonArr)[leptonIndex[q]]);
                    if (leptonType[q] == 13 && pf->matchedObjectType == 13 && pf->matchedObjectIndex == leptonIndex[q]) {
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

            //***********************
            //Fill isolation rings
            //***********************
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
            Double_t EleIsoMVA = ElectronIsoReader->EvaluateMVA( "ElectronIsoMVA_BDTG_V0 method" );
            if (passEleIsoMVASameAsDetIso025(ele->pt, ele->scEta, EleIsoMVA)) {
              NLeptonsPassMVAIsoLoose++;
            }
            if (passEleIsoMVASameAsDetIso015(ele->pt, ele->scEta, EleIsoMVA)) {
              NLeptonsPassMVAIsoTight++;
            }

          } else if (leptonType[i] == 13) {
            const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[leptonIndex[i]]);   

            Double_t tmpChargedIso_DR0p0To0p2  = 0;   
            Double_t tmpChargedIso_DR0p2To0p3  = 0;
            Double_t tmpChargedIso_DR0p3To0p4  = 0;
            Double_t tmpChargedIso_DR0p4To0p5  = 0;
            Double_t tmpChargedIso_DR0p5To0p7  = 0;
            Double_t tmpGammaIso_DR0p0To0p2  = 0;   
            Double_t tmpGammaIso_DR0p2To0p3  = 0;
            Double_t tmpGammaIso_DR0p3To0p4  = 0;
            Double_t tmpGammaIso_DR0p4To0p5  = 0;
            Double_t tmpGammaIso_DR0p5To0p7  = 0;
            Double_t tmpNeutralHadronIso_DR0p0To0p2  = 0;   
            Double_t tmpNeutralHadronIso_DR0p2To0p3  = 0;
            Double_t tmpNeutralHadronIso_DR0p3To0p4  = 0;
            Double_t tmpNeutralHadronIso_DR0p4To0p5  = 0;
            Double_t tmpNeutralHadronIso_DR0p5To0p7  = 0;

            //Loop over PF Candidates
            for(Int_t k=0; k<pfcandidateArr->GetEntries(); ++k) {
              const mithep::TPFCandidate *pf = (mithep::TPFCandidate*)((*pfcandidateArr)[k]);
              if (pf->matchedObjectType == 13 && pf->matchedObjectIndex == leptonIndex[i]) continue;

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
                for (Int_t q=0; q < leptonType.size() ; ++q) {
                  if (leptonType[q] == 11) {
                    const mithep::TElectron *tmpele = (mithep::TElectron*)((*electronArr)[leptonIndex[q]]);
                    if (leptonType[q] == 11 && pf->matchedObjectType == 11 && pf->matchedObjectIndex == leptonIndex[q]) {
                      IsLeptonFootprint = kTRUE;
                    }
                    if (pf->q != 0 && fabs(tmpele->scEta) > 1.479 && mithep::MathUtils::DeltaR(tmpele->phi,tmpele->eta, pf->phi, pf->eta) < 0.015) IsLeptonFootprint = kTRUE;
                    if (pf->pfType == eGamma) {
                      if (fabs(tmpele->scEta) > 1.479) {
                        if (mithep::MathUtils::DeltaR(tmpele->phi,tmpele->eta, pf->phi, pf->eta) < 0.08) IsLeptonFootprint = kTRUE;
                      }
                    }
                  }
                  if (leptonType[q] == 13) {
                    const mithep::TMuon *tmpmu = (mithep::TMuon*)((*muonArr)[leptonIndex[q]]);
                    if (leptonType[q] == 13 && pf->matchedObjectType == 13 && pf->matchedObjectIndex == leptonIndex[q]) {
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

                  //cout << "PFCharged: " << pf->pt << " " << abs(pf->dz) << " " << dr << endl;

                  if (dr < 0.2) tmpChargedIso_DR0p0To0p2 += pf->pt;              
                  if (dr >= 0.2 && dr < 0.3) tmpChargedIso_DR0p2To0p3 += pf->pt;
                  if (dr >= 0.3 && dr < 0.4) tmpChargedIso_DR0p3To0p4 += pf->pt;
                  if (dr >= 0.4 && dr < 0.5) tmpChargedIso_DR0p4To0p5 += pf->pt;
                  if (dr >= 0.5 && dr < 0.7) tmpChargedIso_DR0p5To0p7 += pf->pt;
                }
                else if (pf->pfType == eGamma) {
                  //************************************************************
                  // Footprint Veto
//               if (dr < 0.07) continue;
                  //************************************************************
                  //cout << "PFGamma: " << pf->pt << " " << abs(pf->dz) << " " << dr << endl;
                  if (dr < 0.2) tmpGammaIso_DR0p0To0p2 += pf->pt;              
                  if (dr >= 0.2 && dr < 0.3) tmpGammaIso_DR0p2To0p3 += pf->pt;
                  if (dr >= 0.3 && dr < 0.4) tmpGammaIso_DR0p3To0p4 += pf->pt;
                  if (dr >= 0.4 && dr < 0.5) tmpGammaIso_DR0p4To0p5 += pf->pt;
                  if (dr >= 0.5 && dr < 0.7) tmpGammaIso_DR0p5To0p7 += pf->pt;

                }
                else {
                  //************************************************************
                  // Footprint Veto
//               if (dr < 0.10) continue;
                  //************************************************************
                  //cout << "PFNeutralHad: " << pf->pt << " " << abs(pf->dz) << " " << dr << endl;
                  if (dr < 0.2) tmpNeutralHadronIso_DR0p0To0p2 += pf->pt;              
                  if (dr >= 0.2 && dr < 0.3) tmpNeutralHadronIso_DR0p2To0p3 += pf->pt;
                  if (dr >= 0.3 && dr < 0.4) tmpNeutralHadronIso_DR0p3To0p4 += pf->pt;
                  if (dr >= 0.4 && dr < 0.5) tmpNeutralHadronIso_DR0p4To0p5 += pf->pt;
                  if (dr >= 0.5 && dr < 0.7) tmpNeutralHadronIso_DR0p5To0p7 += pf->pt;
                }      
              } //no lepton footprint
            }

            //***********************
            //Fill isolation rings
            //***********************

            varPt = mu->pt;
            varEta = mu->eta;

 //            varChargedIso_DR0p0To0p2   = TMath::Min((tmpChargedIso_DR0p0To0p2)/mu->pt, 2.5);
//             varChargedIso_DR0p2To0p3 = TMath::Min((tmpChargedIso_DR0p2To0p3)/mu->pt, 2.5);
//             varChargedIso_DR0p3To0p4 = TMath::Min((tmpChargedIso_DR0p3To0p4)/mu->pt, 2.5);
//             varChargedIso_DR0p4To0p5 = TMath::Min((tmpChargedIso_DR0p4To0p5)/mu->pt, 2.5);
//             varChargedIso_DR0p5To0p7 = TMath::Min((tmpChargedIso_DR0p5To0p7)/mu->pt, 2.5);
//             varGammaIso_DR0p0To0p2   = TMath::Max(TMath::Min((tmpGammaIso_DR0p0To0p2 - rho*MuonEffectiveArea(kMuGammaIsoDR0p0To0p2, mu->eta, 0))/mu->pt, 2.5), 0.0);
//             varGammaIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpGammaIso_DR0p2To0p3 - rho*MuonEffectiveArea(kMuGammaIsoDR0p2To0p3, mu->eta, 0))/mu->pt, 2.5), 0.0);
//             varGammaIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpGammaIso_DR0p3To0p4 - rho*MuonEffectiveArea(kMuGammaIsoDR0p3To0p4, mu->eta, 0))/mu->pt, 2.5), 0.0);
//             varGammaIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpGammaIso_DR0p4To0p5 - rho*MuonEffectiveArea(kMuGammaIsoDR0p4To0p5, mu->eta, 0))/mu->pt, 2.5), 0.0);
//             varGammaIso_DR0p5To0p7 = TMath::Max(TMath::Min((tmpGammaIso_DR0p5To0p7 - rho*MuonEffectiveArea(kMuGammaIsoDR0p5To0p7, mu->eta, 0))/mu->pt, 2.5), 0.0);
//             varNeutralHadronIso_DR0p0To0p2   = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p0To0p2 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p0To0p2, mu->eta, 0))/mu->pt, 2.5), 0.0);
//             varNeutralHadronIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p2To0p3 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p2To0p3, mu->eta, 0))/mu->pt, 2.5), 0.0);
//             varNeutralHadronIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p3To0p4 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p3To0p4, mu->eta, 0))/mu->pt, 2.5), 0.0);
//             varNeutralHadronIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p4To0p5 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p4To0p5, mu->eta, 0))/mu->pt, 2.5), 0.0);
//             varNeutralHadronIso_DR0p5To0p7 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p5To0p7 - rho*MuonEffectiveArea(kMuNeutralHadronIsoDR0p5To0p7, mu->eta, 0))/mu->pt, 2.5), 0.0);

            
            //Evaluate MVA
            Double_t MuIsoMVA = MuonIsoReader->EvaluateMVA( "MuIsoMVA_BDTG_V1 method" );
            if (passMuonIsoMVASameAsDetIso025(mu->pt, mu->eta, MuIsoMVA)) {
              NLeptonsPassMVAIsoLoose++;
            }
            if (passMuonIsoMVASameAsDetIso015(mu->pt, mu->eta, MuIsoMVA)) {
              NLeptonsPassMVAIsoTight++;
            }

          } else {
            cout << "weird lepton type : " << leptonType[i] << endl;
          }          
        }

        //Count Passing Events
        if (NLeptonsPassMVAIsoLoose >= 4) {
          if (pt_ele.size() == 4) {
            NMVAIsoLoose_4E += eventweight;
          } else if (pt_mu.size() == 4) {
            NMVAIsoLoose_4M += eventweight;
          } else if (pt_ele.size() == 2 && pt_mu.size() == 2) {
            NMVAIsoLoose_2E2M += eventweight;
          }
        }
        if (NLeptonsPassMVAIsoTight >= 4) {
          if (pt_ele.size() == 4) {
            NMVAIsoTight_4E += eventweight;
          } else if (pt_mu.size() == 4) {
            NMVAIsoTight_4M += eventweight;
          } else if (pt_ele.size() == 2 && pt_mu.size() == 2) {
            NMVAIsoTight_2E2M += eventweight;
          }
        }



        //******************************************************************************
        //Event Selection
        //Initialize Event Tree
        //******************************************************************************
        kinematics.l1type = 0;
        kinematics.l2type = 0;
        kinematics.l3type = 0;
        kinematics.l4type = 0;
        kinematics.l1pt = 0;
        kinematics.l2pt = 0;
        kinematics.l3pt = 0;
        kinematics.l4pt = 0;
        kinematics.l1eta = 0;
        kinematics.l2eta = 0;
        kinematics.l3eta = 0;
        kinematics.l4eta = 0;
        kinematics.l1phi = 0;
        kinematics.l2phi = 0;
        kinematics.l3phi = 0;
        kinematics.l4phi = 0;
        kinematics.Z1pt = 0;
        kinematics.Z2pt = 0;
        kinematics.ZZpt = 0;
        kinematics.Z1eta = 0;
        kinematics.Z2eta = 0;
        kinematics.ZZeta = 0;
        kinematics.mZ1 = 0; 
        kinematics.mZ2 = 0; 
        kinematics.m4l = 0;
        kinematics.channel = 0;


        Bool_t passSelection = kTRUE;

        //trick to break execution if event doesn't satisfy certain requirements
        for (Int_t q=0; q<=0;q++) {

          //******************************************************************************
          //Z1 Selection
          //******************************************************************************
          if (leptonPt.size() < 2) break;

          Int_t Z1LeptonPlusIndex = -1;
          Int_t Z1LeptonMinusIndex = -1;
          Double_t BestZ1Mass = -1;

          for(int i = 0; i < leptonPt.size(); ++i) {
            for(int j = i+1; j < leptonPt.size(); ++j) {
              if (!(leptonPt[i] > 20.0 || leptonPt[j] > 20.0)) continue; //one lepton must be 20 GeV
              if (!(leptonPt[i] > 10.0 && leptonPt[j] > 10.0)) continue; //both leptons must be 10 GeV
              if (leptonCharge[i] == leptonCharge[j]) continue;          //require opp sign
              if (fabs(leptonType[i]) != fabs(leptonType[j])) continue;  //require same flavor

              //Make Z1 hypothesis
              mithep::FourVectorM leptonPlus;
              mithep::FourVectorM leptonMinus;

              Double_t Lepton1Mass = 0.51099892e-3;
              if (leptonType[i] != 11) Lepton1Mass = 105.658369e-3;
              Double_t Lepton2Mass = 0.51099892e-3;
              if (leptonType[j] != 11) Lepton2Mass = 105.658369e-3;

              if (leptonCharge[i] > 0 ) {
                leptonPlus.SetCoordinates(leptonPt[i], leptonEta[i], leptonPhi[i], Lepton1Mass );
                leptonMinus.SetCoordinates(leptonPt[j], leptonEta[j], leptonPhi[j], Lepton2Mass );
              } else {
                leptonPlus.SetCoordinates(leptonPt[j], leptonEta[j], leptonPhi[j], Lepton2Mass );
                leptonMinus.SetCoordinates(leptonPt[i], leptonEta[i], leptonPhi[i], Lepton1Mass );
              }
              mithep::FourVectorM dilepton = leptonPlus+leptonMinus;
            
              if (BestZ1Mass < 0 || fabs(dilepton.M() - 91.1876) < fabs(BestZ1Mass - 91.1876)) {
                if (dilepton.M() > 50) {
                  BestZ1Mass = dilepton.M();
                  if (leptonCharge[i] > 0) {
                    Z1LeptonPlusIndex = i;
                    Z1LeptonMinusIndex = j;
                  } else {
                    Z1LeptonPlusIndex = j;
                    Z1LeptonMinusIndex = i;
                  }
                }
              }
            }
          }
        
          if (Z1LeptonPlusIndex == -1) break;

          //******************************************************************************
          //Make Z1 Candidate FourVector
          //******************************************************************************
          mithep::FourVectorM Z1LeptonPlus;
          mithep::FourVectorM Z1LeptonMinus;
          if (leptonType[Z1LeptonPlusIndex] == 11) {
            Z1LeptonPlus.SetCoordinates(leptonPt[Z1LeptonPlusIndex], leptonEta[Z1LeptonPlusIndex], leptonPhi[Z1LeptonPlusIndex], 0.51099892e-3 );
          } else {
            Z1LeptonPlus.SetCoordinates(leptonPt[Z1LeptonPlusIndex], leptonEta[Z1LeptonPlusIndex], leptonPhi[Z1LeptonPlusIndex], 105.658369e-3 );
          }
          if (leptonType[Z1LeptonMinusIndex] == 11) {
            Z1LeptonMinus.SetCoordinates(leptonPt[Z1LeptonMinusIndex], leptonEta[Z1LeptonMinusIndex], leptonPhi[Z1LeptonMinusIndex], 0.51099892e-3 );
          } else {
            Z1LeptonMinus.SetCoordinates(leptonPt[Z1LeptonMinusIndex], leptonEta[Z1LeptonMinusIndex], leptonPhi[Z1LeptonMinusIndex], 105.658369e-3 );
          }
          mithep::FourVectorM Z1Candidate = Z1LeptonPlus+Z1LeptonMinus;
  

          //******************************************************************************
          //Z2 Selection
          //******************************************************************************
          Int_t Z2LeptonPlusIndex = -1;
          Int_t Z2LeptonMinusIndex = -1;
          Double_t BestZ2Mass = -1;
          for(int i = 0; i < leptonPt.size(); ++i) {
            for(int j = i+1; j < leptonPt.size(); ++j) {
              if (i == Z1LeptonPlusIndex || i == Z1LeptonMinusIndex) continue; //skip Z1 leptons
              if (j == Z1LeptonPlusIndex || j == Z1LeptonMinusIndex) continue; //skip Z1 leptons
              if (leptonCharge[i] == leptonCharge[j]) continue;         //require opp sign
              if (fabs(leptonType[i]) != fabs(leptonType[j])) continue; //require same flavor
            
              //Make Z2 hypothesis
              mithep::FourVectorM leptonPlus;
              mithep::FourVectorM leptonMinus;

              Double_t Lepton1Mass = 0.51099892e-3;
              if (leptonType[i] != 11) Lepton1Mass = 105.658369e-3;
              Double_t Lepton2Mass = 0.51099892e-3;
              if (leptonType[j] != 11) Lepton2Mass = 105.658369e-3;

              if (leptonCharge[i] > 0 ) {
                leptonPlus.SetCoordinates(leptonPt[i], leptonEta[i], leptonPhi[i], Lepton1Mass );
                leptonMinus.SetCoordinates(leptonPt[j], leptonEta[j], leptonPhi[j], Lepton2Mass );
              } else {
                leptonPlus.SetCoordinates(leptonPt[j], leptonEta[j], leptonPhi[j], Lepton2Mass );
                leptonMinus.SetCoordinates(leptonPt[i], leptonEta[i], leptonPhi[i], Lepton1Mass );
              }

              mithep::FourVectorM dilepton = leptonPlus+leptonMinus;
              mithep::FourVectorM fourLepton = Z1Candidate + dilepton;

              if (!(dilepton.M() > 12.0)) continue;
              if (!(fourLepton.M() > 100.0)) continue;
            
              //for 4e and 4mu, require at least 1 of the other opp sign lepton pairs have mass > 12
              if (fabs(leptonType[i]) == fabs(leptonType[Z1LeptonPlusIndex])) {
                mithep::FourVectorM pair1 = Z1LeptonPlus+leptonMinus;
                mithep::FourVectorM pair2 = Z1LeptonMinus+leptonPlus;
                if (!(pair1.M() > 12 || pair2.M() > 12)) continue;
              }
            
              //Disambiguiation is done by choosing the pair with the largest ptMax, and largest ptMin
              if (Z2LeptonPlusIndex < 0) {
                if (leptonCharge[i] > 0) {
                  Z2LeptonPlusIndex = i;
                  Z2LeptonMinusIndex = j;
                } else {
                  Z2LeptonPlusIndex = j;
                  Z2LeptonMinusIndex = i;
                }
              } else {
                Double_t BestPairPtMax = leptonPt[Z2LeptonPlusIndex];               
                Double_t BestPairPtMin = leptonPt[Z2LeptonMinusIndex]; 
                if (leptonPt[Z2LeptonMinusIndex] > BestPairPtMax) {
                  BestPairPtMax = leptonPt[Z2LeptonMinusIndex];
                  BestPairPtMin = leptonPt[Z2LeptonPlusIndex];
                }

                Double_t CurrentPairPtMax = leptonPt[i];               
                Double_t CurrentPairPtMin = leptonPt[j]; 
                if (leptonPt[j] > CurrentPairPtMax) {
                  CurrentPairPtMax = leptonPt[j];
                  CurrentPairPtMin = leptonPt[i];
                }

                if (CurrentPairPtMax > BestPairPtMax) {
                  if (leptonCharge[i] > 0) {
                    Z2LeptonPlusIndex = i;
                    Z2LeptonMinusIndex = j;
                  } else {
                    Z2LeptonPlusIndex = j;
                    Z2LeptonMinusIndex = i;
                  }
                } else if (CurrentPairPtMax  == BestPairPtMax) {
                  if (CurrentPairPtMin > BestPairPtMin) {
                    if (leptonCharge[i] > 0) {
                      Z2LeptonPlusIndex = i;
                      Z2LeptonMinusIndex = j;
                    } else {
                      Z2LeptonPlusIndex = j;
                      Z2LeptonMinusIndex = i;
                    }                  
                  }
                }
              }            
            }
          }
 
          if (Z2LeptonPlusIndex == -1) break;

          //******************************************************************************
          //Make Z2 Candidate FourVector
          //******************************************************************************
          mithep::FourVectorM Z2LeptonPlus;
          mithep::FourVectorM Z2LeptonMinus;
          if (leptonType[Z2LeptonPlusIndex] == 11) {
            Z2LeptonPlus.SetCoordinates(leptonPt[Z2LeptonPlusIndex], leptonEta[Z2LeptonPlusIndex], leptonPhi[Z2LeptonPlusIndex], 0.51099892e-3 );
          } else {
            Z2LeptonPlus.SetCoordinates(leptonPt[Z2LeptonPlusIndex], leptonEta[Z2LeptonPlusIndex], leptonPhi[Z2LeptonPlusIndex], 105.658369e-3 );
          }
          if (leptonType[Z2LeptonMinusIndex] == 11) {
            Z2LeptonMinus.SetCoordinates(leptonPt[Z2LeptonMinusIndex], leptonEta[Z2LeptonMinusIndex], leptonPhi[Z2LeptonMinusIndex], 0.51099892e-3 );
          } else {
            Z2LeptonMinus.SetCoordinates(leptonPt[Z2LeptonMinusIndex], leptonEta[Z2LeptonMinusIndex], leptonPhi[Z2LeptonMinusIndex], 105.658369e-3 );
          }
          mithep::FourVectorM Z2Candidate = Z2LeptonPlus+Z2LeptonMinus;
          mithep::FourVectorM ZZSystem = Z1Candidate + Z2Candidate;


          //***************************************************************
          // remaining kinematic cuts 
          //***************************************************************
          if (Z1Candidate.M() > 120) break;
          if (Z2Candidate.M() > 120) break;
          if (Z2Candidate.M() < 12) break;
                   
          //***************************************************************
          // Define Channel Type
          //***************************************************************             
          Int_t Channel = -1;
          if (leptonType[Z1LeptonPlusIndex] == 11 && leptonType[Z2LeptonPlusIndex] == 11 ) {
            Channel = kFourEle;
          } else if (leptonType[Z1LeptonPlusIndex] == 13 && leptonType[Z2LeptonPlusIndex] == 13 ) {
            Channel = kFourMu;
          } else if (leptonType[Z1LeptonPlusIndex] == 11 && leptonType[Z2LeptonPlusIndex] == 13 ) {
            Channel = kTwoEleTwoMu;
          } else if (leptonType[Z1LeptonPlusIndex] == 13 && leptonType[Z2LeptonPlusIndex] == 11 ) {
            Channel = kTwoMuTwoEle;
          }
        

          //***************************************************************
          // Fill Kinematics Ntuple
          //***************************************************************             
          kinematics.l1type = leptonType[Z1LeptonPlusIndex];
          kinematics.l2type = leptonType[Z1LeptonMinusIndex];
          kinematics.l3type = leptonType[Z2LeptonPlusIndex];
          kinematics.l4type = leptonType[Z2LeptonMinusIndex];
          kinematics.l1pt = leptonPt[Z1LeptonPlusIndex];
          kinematics.l2pt = leptonPt[Z1LeptonMinusIndex];
          kinematics.l3pt = leptonPt[Z2LeptonPlusIndex];
          kinematics.l4pt = leptonPt[Z2LeptonMinusIndex];
          kinematics.l1eta = leptonEta[Z1LeptonPlusIndex];
          kinematics.l2eta = leptonEta[Z1LeptonPlusIndex];
          kinematics.l3eta = leptonEta[Z2LeptonMinusIndex];
          kinematics.l4eta = leptonEta[Z2LeptonMinusIndex];
          kinematics.l1phi = leptonEta[Z1LeptonPlusIndex];
          kinematics.l2phi = leptonEta[Z1LeptonMinusIndex];
          kinematics.l3phi = leptonEta[Z2LeptonPlusIndex];
          kinematics.l4phi = leptonEta[Z2LeptonMinusIndex];
          kinematics.Z1pt = Z1Candidate.Pt();
          kinematics.Z2pt = Z2Candidate.Pt();
          kinematics.ZZpt = ZZSystem.Pt();
          kinematics.Z1eta = Z1Candidate.Eta();
          kinematics.Z2eta = Z2Candidate.Eta();
          kinematics.ZZeta = ZZSystem.Eta();
          kinematics.mZ1 = Z1Candidate.M(); 
          kinematics.mZ2 = Z2Candidate.M(); 
          kinematics.m4l = ZZSystem.M();
          kinematics.channel = Channel;

          //***************************************************************
          // Count Events
          //***************************************************************             
          //ID Only
          if (pt_ele.size() == 4) {
            NSelectedBaselineID_4E += eventweight;
          } else if (pt_mu.size() == 4) {
            NSelectedBaselineID_4M += eventweight;
          } else if (pt_ele.size() == 2 && pt_mu.size() == 2) {
            NSelectedBaselineID_2E2M += eventweight;
          }

          //DetIso025
          if (NLeptonsPassDetIso025 >= 4) {
            if (pt_ele.size() == 4) {
              NSelectedDetIso025_4E += eventweight;
            } else if (pt_mu.size() == 4) {
              NSelectedDetIso025_4M += eventweight;
            } else if (pt_ele.size() == 2 && pt_mu.size() == 2) {
              NSelectedDetIso025_2E2M += eventweight;
            }
          }
          //Pairwise
          if (PassPairwiseIso) {
            if (pt_ele.size() == 4) {
              NSelectedBaselinePairwise_4E += eventweight;
            } else if (pt_mu.size() == 4) {
              NSelectedBaselinePairwise_4M += eventweight;
            } else if (pt_ele.size() == 2 && pt_mu.size() == 2) {
              NSelectedBaselinePairwise_2E2M += eventweight;
            }
          }

          //MVA Iso
          if (NLeptonsPassMVAIsoLoose >= 4) {
            if (pt_ele.size() == 4) {
              NSelectedMVAIsoLoose_4E += eventweight;
            } else if (pt_mu.size() == 4) {
              NSelectedMVAIsoLoose_4M += eventweight;
            } else if (pt_ele.size() == 2 && pt_mu.size() == 2) {
              NSelectedMVAIsoLoose_2E2M += eventweight;
            }
          }
          if (NLeptonsPassMVAIsoTight >= 4) {
            if (pt_ele.size() == 4) {
              NSelectedMVAIsoTight_4E += eventweight;
            } else if (pt_mu.size() == 4) {
              NSelectedMVAIsoTight_4M += eventweight;
            } else if (pt_ele.size() == 2 && pt_mu.size() == 2) {
              NSelectedMVAIsoTight_2E2M += eventweight;
            }
          }



        } // dummy loop over current one event

        //Fill Tree
        if (Option == 0) {
          fOutputTree->Fill();
        } else if (Option == 1) {
          if (PassPairwiseIso) {
            fOutputTree->Fill();
          }
        } else if (Option == 2) {
          if (NLeptonsPassMVAIsoLoose >= 4) {
            fOutputTree->Fill();
          }
        } else if (Option == 3) {
          if (NLeptonsPassMVAIsoTight >= 4) {
            fOutputTree->Fill();
          }
        } else if (Option == 100) {
          if (NLeptonsPassMVAIsoLoose >= 4 && !PassPairwiseIso) {
            fOutputTree->Fill();
          }
        }


      } //end loop over data     
    }
  }

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;


  cout << "Number of GenAccepted Events: " << NGenAccepted_4E << " " << NGenAccepted_4M << " " << NGenAccepted_2E2M << endl;
  cout << "Number of BaselineID'd Events: " << NBaselineID_4E << " " << NBaselineID_4M << " " << NBaselineID_2E2M << endl;
  cout << "Number of BaselinePairwise'd Events: " << NBaselinePairwise_4E << " " << NBaselinePairwise_4M << " " << NBaselinePairwise_2E2M << endl;
  cout << "Number of DetIso025'd Events: " << NDetIso025_4E << " " << NDetIso025_4M << " " << NDetIso025_2E2M << endl;
  cout << "Number of MVAIsoLoose'd Events: " << NMVAIsoLoose_4E << " " << NMVAIsoLoose_4M << " " << NMVAIsoLoose_2E2M << endl;
  cout << "Number of MVAIsoTight'd Events: " << NMVAIsoTight_4E << " " << NMVAIsoTight_4M << " " << NMVAIsoTight_2E2M << endl;
  cout << "Selected Number of BaselineID'd Events: " << NSelectedBaselineID_4E << " " << NSelectedBaselineID_4M << " " << NSelectedBaselineID_2E2M << endl;
  cout << "Selected Number of BaselinePairwise'd Events: " << NSelectedBaselinePairwise_4E << " " << NSelectedBaselinePairwise_4M << " " << NSelectedBaselinePairwise_2E2M << endl;
  cout << "Selected Number of DetIso025'd Events: " << NSelectedDetIso025_4E << " " << NSelectedDetIso025_4M << " " << NSelectedDetIso025_2E2M << endl;
  cout << "Selected Number of MVAIsoLoose'd Events: " << NSelectedMVAIsoLoose_4E << " " << NSelectedMVAIsoLoose_4M << " " << NSelectedMVAIsoLoose_2E2M << endl;
  cout << "Selected Number of MVAIsoTight'd Events: " << NSelectedMVAIsoTight_4E << " " << NSelectedMVAIsoTight_4M << " " << NSelectedMVAIsoTight_2E2M << endl;
  cout << endl;
  cout << "Efficiency of BaselineID: " << NBaselineID_4E/NGenAccepted_4E << " " << NBaselineID_4M/NGenAccepted_4M << " " << NBaselineID_2E2M/NGenAccepted_2E2M << endl;
  cout << "Efficiency of BaselinePairwise: " << NBaselinePairwise_4E/NGenAccepted_4E << " " << NBaselinePairwise_4M/NGenAccepted_4M << " " << NBaselinePairwise_2E2M/NGenAccepted_2E2M << endl;
  cout << "Efficiency of DetIso025: " << NDetIso025_4E/NGenAccepted_4E << " " << NDetIso025_4M/NGenAccepted_4M << " " << NDetIso025_2E2M/NGenAccepted_2E2M << endl;
  cout << "Efficiency of MVAIsoLoose: " << NMVAIsoLoose_4E/NGenAccepted_4E << " " << NMVAIsoLoose_4M/NGenAccepted_4M << " " << NMVAIsoLoose_2E2M/NGenAccepted_2E2M << endl;
  cout << "Efficiency of MVAIsoTight: " << NMVAIsoTight_4E/NGenAccepted_4E << " " << NMVAIsoTight_4M/NGenAccepted_4M << " " << NMVAIsoTight_2E2M/NGenAccepted_2E2M << endl;

  cout << "Efficiency of Selected BaselineID: " << NSelectedBaselineID_4E/NGenAccepted_4E << " " << NSelectedBaselineID_4M/NGenAccepted_4M << " " << NSelectedBaselineID_2E2M/NGenAccepted_2E2M << endl;
  cout << "Efficiency of Selected BaselinePairwise: " << NSelectedBaselinePairwise_4E/NGenAccepted_4E << " " << NSelectedBaselinePairwise_4M/NGenAccepted_4M << " " << NSelectedBaselinePairwise_2E2M/NGenAccepted_2E2M << endl;
  cout << "Efficiency of Selected DetIso025: " << NSelectedDetIso025_4E/NGenAccepted_4E << " " << NSelectedDetIso025_4M/NGenAccepted_4M << " " << NSelectedDetIso025_2E2M/NGenAccepted_2E2M << endl;
  cout << "Efficiency of Selected MVAIsoLoose: " << NSelectedMVAIsoLoose_4E/NGenAccepted_4E << " " << NSelectedMVAIsoLoose_4M/NGenAccepted_4M << " " << NSelectedMVAIsoLoose_2E2M/NGenAccepted_2E2M << endl;
  cout << "Efficiency of Selected MVAIsoTight: " << NSelectedMVAIsoTight_4E/NGenAccepted_4E << " " << NSelectedMVAIsoTight_4M/NGenAccepted_4M << " " << NSelectedMVAIsoTight_2E2M/NGenAccepted_2E2M << endl;

  fOutputFile->cd();
  fOutputFile->Write();
  fOutputFile->Close();
  

  
  gBenchmark->Show("WWTemplate");       
} 



