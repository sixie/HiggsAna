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

void TrainElectronIsoMVA( string Label,
                          Int_t Option = -1,
                          string versionLabel = "") {
  
  string label = versionLabel;
  if (Label != "") label = Label + "_";

  TMVA::Tools::Instance();
  TFile* outputFile = TFile::Open( (label+versionLabel+".root").c_str(), "RECREATE");
  TMVA::Factory *factory = new TMVA::Factory((label+versionLabel).c_str(), outputFile, "!V:!Silent");
  
  factory->AddSpectator( "eta");  
  factory->AddSpectator( "pt");
  factory->AddVariable("ChargedIso_DR0p0To0p1", 'F');
  factory->AddVariable("ChargedIso_DR0p1To0p2", 'F');
  factory->AddVariable("ChargedIso_DR0p2To0p3", 'F');
  factory->AddVariable("ChargedIso_DR0p3To0p4", 'F');
  factory->AddVariable("ChargedIso_DR0p4To0p5", 'F');

  factory->AddVariable("GammaIso_DR0p0To0p1", 'F');
  factory->AddVariable("GammaIso_DR0p1To0p2", 'F');
  factory->AddVariable("GammaIso_DR0p2To0p3", 'F');
  factory->AddVariable("GammaIso_DR0p3To0p4", 'F');
  factory->AddVariable("GammaIso_DR0p4To0p5", 'F');

  factory->AddVariable("NeutralHadronIso_DR0p0To0p1", 'F');
  factory->AddVariable("NeutralHadronIso_DR0p1To0p2", 'F');
  factory->AddVariable("NeutralHadronIso_DR0p2To0p3", 'F');
  factory->AddVariable("NeutralHadronIso_DR0p3To0p4", 'F');
  factory->AddVariable("NeutralHadronIso_DR0p4To0p5", 'F');



  //Split using parity of event number
//   TFile* inputSignalTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/ElectronSelectionTraining.s11HZZ.Training.root");   
//   TFile* inputBkgTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/ElectronSelectionTraining.Fake_ZPlusJet.Training.root");
//   TFile* inputSignalTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/ElectronSelectionTraining.s11HZZ.Training.root");   
//   TFile* inputBkgTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/ElectronSelectionTraining.Fake_ZPlusJet.Testing.root");

// //   //WJets vs ZJets
//   TFile* inputSignalTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/ElectronSelectionTraining.Fake_WPlusJet.root");   
//   TFile* inputBkgTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/ElectronSelectionTraining.Fake_ZPlusJet.root");
//   TFile* inputSignalTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/ElectronSelectionTraining.Fake_WPlusJet.root");   
//   TFile* inputBkgTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/ElectronSelectionTraining.Fake_ZPlusJet.root");

  //WENu vs WMuNu
 //  TFile* inputSignalTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/ElectronSelectionTraining.Fake_WToENuPlusJet.root");   
//   TFile* inputBkgTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/ElectronSelectionTraining.Fake_WToMuNuPlusJet.root");
//   TFile* inputSignalTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/ElectronSelectionTraining.Fake_WToENuPlusJet.root");   
//   TFile* inputBkgTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/ElectronSelectionTraining.Fake_WToMuNuPlusJet.root");

  //HZZ vs WJets 
  TFile* inputSignalTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/ElectronSelectionTraining.s11HZZ.Training.root");   
  TFile* inputBkgTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/ElectronSelectionTraining.Fake_WPlusJet.root");
  TFile* inputSignalTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/ElectronSelectionTraining.s11HZZ.Training.root");   
  TFile* inputBkgTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/ElectronSelectionTraining.Fake_WPlusJet.root");

  TTree *signalTraining     = (TTree*)inputSignalTraining->Get("Electrons");
  TTree *backgroundTraining = (TTree*)inputBkgTraining->Get("Electrons");
  TTree *signalTesting     = (TTree*)inputSignalTesting->Get("Electrons");
  TTree *backgroundTesting = (TTree*)inputBkgTesting->Get("Electrons");  

  assert(signalTraining);
  assert(signalTesting);
  assert(backgroundTraining);
  assert(backgroundTesting);
  factory->AddSignalTree    (signalTraining,1.0,TMVA::Types::kTraining);
  factory->AddBackgroundTree(backgroundTraining,1.0,TMVA::Types::kTraining);
  factory->AddSignalTree    (signalTesting,1.0,TMVA::Types::kTesting);
  factory->AddBackgroundTree(backgroundTesting,1.0,TMVA::Types::kTesting);

  if (Option == -1) {
    factory->PrepareTrainingAndTestTree("pt > 5 && pt <= 35",
                                        "pt > 5 && pt <= 35",
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:!V" );
  }
  if (Option == 0) {
    factory->PrepareTrainingAndTestTree("pt > 5 && pt <= 10 && abs(eta) <= 1.479",
                                        "pt > 5 && pt <= 10 && abs(eta) <= 1.479",
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:!V" );
  }
  if (Option == 1) {
    factory->PrepareTrainingAndTestTree("pt > 5 && pt <= 10 && abs(eta) > 1.479",
                                        "pt > 5 && pt <= 10 && abs(eta) > 1.479",
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:!V" );
  }
  if (Option == 2) {
    factory->PrepareTrainingAndTestTree("pt > 10 && pt <= 35 && abs(eta) <= 1.479",
                                        "pt > 10 && pt <= 35 && abs(eta) <= 1.479",
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:!V" );
  }
  if (Option == 3) {
    factory->PrepareTrainingAndTestTree("pt > 10 && pt <= 35 && abs(eta) > 1.479",
                                        "pt > 10 && pt <= 35 && abs(eta) > 1.479",
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:!V" );
  }

  //*****************************************************************
  //V0: Only RelPt Rings with size 0.10
  //*****************************************************************
  if (versionLabel == "" || versionLabel == "V0" || versionLabel == "V1" || versionLabel == "V2" ) {

    if (Option >= 0) {
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=5:UseNvars=4:PruneStrength=0.0:PruneMethod=nopruning:MaxDepth=6" );
//      factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
//                           "H:!V:CreateMVAPdfs:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=30:NSmoothBkg[0]=30:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=40:Nbins=200" );
    }
    else if (Option == -1) {
      TMVA::MethodBase* bdtCat_BDTG_V0 = factory->BookMethod( TMVA::Types::kCategory, (label+"BDTCat_BDTG_V0").c_str(),"" );
      TMVA::MethodCategory* category_BDTG_V0 = dynamic_cast<TMVA::MethodCategory*>(bdtCat_BDTG_V0);
      category_BDTG_V0->AddMethod( "pt <= 10 && abs(eta) <= 1.5",
                                   "ChargedIso_DR0p0To0p1:ChargedIso_DR0p1To0p2:ChargedIso_DR0p2To0p3:ChargedIso_DR0p3To0p4:ChargedIso_DR0p4To0p5:GammaIso_DR0p0To0p1:GammaIso_DR0p1To0p2:GammaIso_DR0p2To0p3:GammaIso_DR0p3To0p4:GammaIso_DR0p4To0p5:NeutralHadronIso_DR0p0To0p1:NeutralHadronIso_DR0p1To0p2:NeutralHadronIso_DR0p2To0p3:NeutralHadronIso_DR0p3To0p4:NeutralHadronIso_DR0p4To0p5",
                                   TMVA::Types::kBDT, 
                                   "Category_BDTG_V0_BarrelPt5To10",
                                   "!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=5:UseNvars=4:PruneStrength=0.0:PruneMethod=nopruning:MaxDepth=6" );
      
      category_BDTG_V0->AddMethod( "pt <= 10 && abs(eta) > 1.5",
                                   "ChargedIso_DR0p0To0p1:ChargedIso_DR0p1To0p2:ChargedIso_DR0p2To0p3:ChargedIso_DR0p3To0p4:ChargedIso_DR0p4To0p5:GammaIso_DR0p0To0p1:GammaIso_DR0p1To0p2:GammaIso_DR0p2To0p3:GammaIso_DR0p3To0p4:GammaIso_DR0p4To0p5:NeutralHadronIso_DR0p0To0p1:NeutralHadronIso_DR0p1To0p2:NeutralHadronIso_DR0p2To0p3:NeutralHadronIso_DR0p3To0p4:NeutralHadronIso_DR0p4To0p5",
                                   TMVA::Types::kBDT, 
                                   "Category_BDTG_V0_EndcapPt5To10",
                                   "!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=5:UseNvars=4:PruneStrength=0.0:PruneMethod=nopruning:MaxDepth=6" );
      
      category_BDTG_V0->AddMethod( "pt > 10 && pt <= 20 && abs(eta) <= 1.5",
                                   "ChargedIso_DR0p0To0p1:ChargedIso_DR0p1To0p2:ChargedIso_DR0p2To0p3:ChargedIso_DR0p3To0p4:ChargedIso_DR0p4To0p5:GammaIso_DR0p0To0p1:GammaIso_DR0p1To0p2:GammaIso_DR0p2To0p3:GammaIso_DR0p3To0p4:GammaIso_DR0p4To0p5:NeutralHadronIso_DR0p0To0p1:NeutralHadronIso_DR0p1To0p2:NeutralHadronIso_DR0p2To0p3:NeutralHadronIso_DR0p3To0p4:NeutralHadronIso_DR0p4To0p5",
                                   TMVA::Types::kBDT, 
                                   "Category_BDTG_V0_BarrelPt10To20",
                                   "!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=5:UseNvars=4:PruneStrength=0.0:PruneMethod=nopruning:MaxDepth=6" );
      
      category_BDTG_V0->AddMethod( "pt > 10 && pt <= 20 && abs(eta) > 1.5",
                                   "ChargedIso_DR0p0To0p1:ChargedIso_DR0p1To0p2:ChargedIso_DR0p2To0p3:ChargedIso_DR0p3To0p4:ChargedIso_DR0p4To0p5:GammaIso_DR0p0To0p1:GammaIso_DR0p1To0p2:GammaIso_DR0p2To0p3:GammaIso_DR0p3To0p4:GammaIso_DR0p4To0p5:NeutralHadronIso_DR0p0To0p1:NeutralHadronIso_DR0p1To0p2:NeutralHadronIso_DR0p2To0p3:NeutralHadronIso_DR0p3To0p4:NeutralHadronIso_DR0p4To0p5",
                                   TMVA::Types::kBDT, 
                                   "Category_BDTG_V0_EndcapPt10To20",
                                   "!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=5:UseNvars=4:PruneStrength=0.0:PruneMethod=nopruning:MaxDepth=6" );
      
      category_BDTG_V0->AddMethod( "pt > 20 && abs(eta) <= 1.5",
                                   "ChargedIso_DR0p0To0p1:ChargedIso_DR0p1To0p2:ChargedIso_DR0p2To0p3:ChargedIso_DR0p3To0p4:ChargedIso_DR0p4To0p5:GammaIso_DR0p0To0p1:GammaIso_DR0p1To0p2:GammaIso_DR0p2To0p3:GammaIso_DR0p3To0p4:GammaIso_DR0p4To0p5:NeutralHadronIso_DR0p0To0p1:NeutralHadronIso_DR0p1To0p2:NeutralHadronIso_DR0p2To0p3:NeutralHadronIso_DR0p3To0p4:NeutralHadronIso_DR0p4To0p5",
                                   TMVA::Types::kBDT, 
                                   "Category_BDTG_V0_BarrelPt20ToInf",
                                   "!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=5:UseNvars=4:PruneStrength=0.0:PruneMethod=nopruning:MaxDepth=6" );
      
      category_BDTG_V0->AddMethod( "pt > 20 && abs(eta) > 1.5",
                                 "ChargedIso_DR0p0To0p1:ChargedIso_DR0p1To0p2:ChargedIso_DR0p2To0p3:ChargedIso_DR0p3To0p4:ChargedIso_DR0p4To0p5:GammaIso_DR0p0To0p1:GammaIso_DR0p1To0p2:GammaIso_DR0p2To0p3:GammaIso_DR0p3To0p4:GammaIso_DR0p4To0p5:NeutralHadronIso_DR0p0To0p1:NeutralHadronIso_DR0p1To0p2:NeutralHadronIso_DR0p2To0p3:NeutralHadronIso_DR0p3To0p4:NeutralHadronIso_DR0p4To0p5",
                                 TMVA::Types::kBDT, 
                                 "Category_BDTG_V0_EndcapPt20ToInf",
                                 "!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=5:UseNvars=4:PruneStrength=0.0:PruneMethod=nopruning:MaxDepth=6" );
    }
  }




  
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  
  outputFile->Close();

  delete factory;
}
