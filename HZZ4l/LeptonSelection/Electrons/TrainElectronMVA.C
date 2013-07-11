#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodCategory.h"
#endif

void TrainElectronMVA( 
  TString SignalTrainingFile   = "ElectronSelectionTraining.Real.Training.root",
  TString SignalTestingFile    = "ElectronSelectionTraining.Real.Testing.root",
  TString BkgTrainingFile      = "ElectronSelectionTraining.Fake.Training.root",
  TString BkgTestingFile       = "ElectronSelectionTraining.Fake.Testing.root",
  TString label                = "default",
  string  versionLabel         = "",
  Int_t   Option               = -1,
  TString myMethodList         = "BDTG"
) {
  
  TMVA::Tools::Instance();
  TString outfileName = Label + ".root";
  TFile* outputFile = TFile::Open(outfileName, "RECREATE");
  TMVA::Factory *factory = new TMVA::Factory("TMVA", outputFile, "!V:!Silent");
  
  factory->AddVariable("see", 'F');
  factory->AddVariable("deta", 'F');
  factory->AddVariable("dphi", 'F');  
  factory->AddVariable("d0", 'F');
  factory->AddVariable("fbrem", 'F');
  factory->AddVariable("EoP", 'F');
  factory->AddVariable("EoPout", 'F');
  factory->AddVariable("spp", 'F');
  factory->AddVariable("IoEmIoP", 'F');
  factory->AddVariable("EoPin", 'F');
  factory->AddVariable("ip3d", 'F');
  factory->AddVariable("ip3ds", 'F');

  factory->AddSpectator( "eta");  
  factory->AddSpectator( "pt");
  factory->AddSpectator( "EventNumberParity");
  

  //Split using parity of event number
  TFile* inputSignalTraining = TFile::Open("ElectronSelectionTraining.Real.Training.root");   
  TFile* inputBkgTraining = TFile::Open("ElectronSelectionTraining.Fakes.Training.root");
  TFile* inputSignalTesting = TFile::Open("ElectronSelectionTraining.Real.Testing.root");   
  TFile* inputBkgTesting = TFile::Open("ElectronSelectionTraining.Fakes.Testing.root");
  TTree *signalTraining     = (TTree*)inputSignalTraining->Get("Electrons");
  TTree *backgroundTraining = (TTree*)inputBkgTraining->Get("Electrons");
  TTree *signalTesting     = (TTree*)inputSignalTesting->Get("Electrons");
  TTree *backgroundTesting = (TTree*)inputBkgTesting->Get("Electrons");  
  factory->AddSignalTree    (signalTraining,1.0,TMVA::Types::kTraining);
  factory->AddBackgroundTree(backgroundTraining,1.0,TMVA::Types::kTraining);
  factory->AddSignalTree    (signalTesting,1.0,TMVA::Types::kTesting);
  factory->AddBackgroundTree(backgroundTesting,1.0,TMVA::Types::kTesting);
  factory->PrepareTrainingAndTestTree("pt > 10 && pt <= 20 && abs(dz) < 0.1 && abs(d0) < 0.02 && matchConv == 0 && missHits == 0 && combPFIsoHWW < 0.4 && DenomFakeSmurf == 1 ",
                                      "pt > 10 && pt <= 20 && abs(dz) < 0.1 && abs(d0) < 0.02 && matchConv == 0 && missHits == 0 && combPFIsoHWW < 0.4 && DenomFakeSmurf == 1 ",
				      "nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:!V" );


 
  //*****************************************************************
  //Adaptive
  //*****************************************************************

  TMVA::MethodBase* bdtCat_Adaptive = factory->BookMethod( TMVA::Types::kCategory, "BDTCat_Adaptive","" );
  TMVA::MethodCategory* category_BDTAdaptive = dynamic_cast<TMVA::MethodCategory*>(bdtCat_Adaptive);
  category_BDTAdaptive->AddMethod( "abs(eta) <= 0.8",
			 "see:deta:dphi:d0:fbrem:EoP:EoPout:spp:IoEmIoP:EoPin:ip3d:ip3ds",
			 TMVA::Types::kBDT, 
			 "Category_BDTAdaptive_1",
			 "!H:!V:NTrees=800:NNodesMax=1000:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:nCuts=2000:SeparationType=GiniIndex:PruneMethod=CostComplexity:PruneStrength=5" );
  


  category_BDTAdaptive->AddMethod( "abs(eta) > 0.8 && abs(eta)<=1.485",
			 "see:deta:dphi:d0:fbrem:EoP:EoPout:spp:IoEmIoP:EoPin:ip3d:ip3ds",
			 TMVA::Types::kBDT, 
			 "Category_BDTAdaptive_2",
			 "!H:!V:NTrees=800:NNodesMax=1000:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:nCuts=2000:SeparationType=GiniIndex:PruneMethod=CostComplexity:PruneStrength=5" );
  
 
 
  category_BDTAdaptive->AddMethod( "abs(eta)> 1.485",
			 "see:deta:dphi:d0:fbrem:EoP:EoPout:spp:IoEmIoP:EoPin:ip3d:ip3ds",
			 TMVA::Types::kBDT, 
			 "Category_BDTAdaptive_3",
			 "!H:!V:NTrees=800:NNodesMax=1000:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:nCuts=2000:SeparationType=GiniIndex:PruneMethod=CostComplexity:PruneStrength=5" );
  



  //*****************************************************************
  //BDTG
  //*****************************************************************

  TMVA::MethodBase* bdtCat_BDTG = factory->BookMethod( TMVA::Types::kCategory, "BDTCat_BDTG","" );
  TMVA::MethodCategory* category_BDTG = dynamic_cast<TMVA::MethodCategory*>(bdtCat_BDTG);
  category_BDTG->AddMethod( "abs(eta) <= 0.8",
			 "see:deta:dphi:d0:fbrem:EoP:EoPout:spp:IoEmIoP:EoPin:ip3d:ip3ds",
			 TMVA::Types::kBDT, 
			 "Category_BDTG_1",
			 "!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:NNodesMax=5:UseNvars=4:PruneStrength=0.0:PruneMethod=nopruning:MaxDepth=6" );
  


  category_BDTG->AddMethod( "abs(eta) > 0.8 && abs(eta)<=1.485",
			 "see:deta:dphi:d0:fbrem:EoP:EoPout:spp:IoEmIoP:EoPin:ip3d:ip3ds",
			 TMVA::Types::kBDT, 
			 "Category_BDTG_2",
			 "!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:NNodesMax=5:UseNvars=4:PruneStrength=0.0:PruneMethod=nopruning:MaxDepth=6" );
  
 
 
  category_BDTG->AddMethod( "abs(eta)> 1.485",
			 "see:deta:dphi:d0:fbrem:EoP:EoPout:spp:IoEmIoP:EoPin:ip3d:ip3ds",
			 TMVA::Types::kBDT, 
			 "Category_BDTG_3",
			 "!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:NNodesMax=5:UseNvars=4:PruneStrength=0.0:PruneMethod=nopruning:MaxDepth=6" );
  
  

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  
  outputFile->Close();

  delete factory;
}
