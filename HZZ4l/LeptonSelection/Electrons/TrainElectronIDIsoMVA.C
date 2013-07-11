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

void TrainElectronIDIsoMVA( string Label,
                            Int_t Option = -1,
                            string versionLabel = "") {
  
  string label = versionLabel;
  if (Label != "") label = Label + "_";
  
  TMVA::Tools::Instance();
  TFile* outputFile = TFile::Open( (label+versionLabel+".root").c_str(), "RECREATE");
  TMVA::Factory *factory = new TMVA::Factory((label+versionLabel).c_str(), outputFile, "!V:!Silent");
  

  factory->AddVariable("fbrem", 'F');
  factory->AddVariable("kfchi2", 'F');
  factory->AddVariable("kfhits", 'F');
  factory->AddVariable("gsfchi2", 'F');
  factory->AddVariable("deta", 'F');
  factory->AddVariable("dphi", 'F');
  factory->AddVariable("detacalo", 'F');
  factory->AddVariable("see", 'F');
  factory->AddVariable("spp", 'F');
  factory->AddVariable("etawidth", 'F');
  factory->AddVariable("phiwidth", 'F');
  factory->AddVariable("e1x5e5x5", 'F');
  factory->AddVariable("R9", 'F');
  factory->AddVariable("HoE", 'F');
  factory->AddVariable("EoP", 'F');
  factory->AddVariable("IoEmIoP", 'F');
  factory->AddVariable("EoPout", 'F');
  if (Option == 2 || Option == 5) {
    factory->AddVariable("PreShowerOverRaw", 'F');
  }
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


  factory->AddSpectator( "eta");  
  factory->AddSpectator( "pt");



  //Split using parity of event number
//   TFile* inputSignalTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.s11HZZ.Training.root");   
//   TFile* inputBkgTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.Fake_ZPlusJet.Training.root");
//   TFile* inputSignalTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.s11HZZ.Training.root");   
//   TFile* inputBkgTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.Fake_ZPlusJet.Testing.root");


  //WJets vs ZJets
//   TFile* inputSignalTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.Fake_WPlusJet.root");   
//   TFile* inputBkgTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.Fake_ZPlusJet.root");
//   TFile* inputSignalTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.Fake_WPlusJet.root");   
//   TFile* inputBkgTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.Fake_ZPlusJet.root");

  //WENu vs WMuNu
//   TFile* inputSignalTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.Fake_WToMuNuPlusJet.root");   
//   TFile* inputBkgTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.Fake_WToENuPlusJet.root");
//   TFile* inputSignalTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.Fake_WToMuNuPlusJet.root");   
//   TFile* inputBkgTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.Fake_WToENuPlusJet.root");


  //HZZ vs WJets 
//   TFile* inputSignalTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.s11HZZ.Training.root");   
//   TFile* inputBkgTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.Fake_WPlusJet.root");
//   TFile* inputSignalTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.s11HZZ.Training.root");   
//   TFile* inputBkgTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.Fake_WPlusJet.root");

  //HZZ vs ZJets MC
  TFile* inputSignalTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.s11HZZ.Training.root");   
  TFile* inputBkgTraining = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.Fall11ZJetsBkg.root");
  TFile* inputSignalTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.s11HZZ.Training.root");   
  TFile* inputBkgTesting = TFile::Open("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/ElectronSelectionTraining.Fake_ZPlusJet.root");

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
    factory->PrepareTrainingAndTestTree("pt > 5 && pt <= 10 && abs(eta) <= 0.8",
                                        "pt > 5 && pt <= 10 && abs(eta) <= 0.8",
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:!V" );
  }
  if (Option == 1) {
    factory->PrepareTrainingAndTestTree("pt > 5 && pt <= 10 && abs(eta) > 0.8 && abs(eta) <= 1.479",
                                        "pt > 5 && pt <= 10 && abs(eta) > 0.8 && abs(eta) <= 1.479",
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:!V" );
  }
  if (Option == 2) {
    factory->PrepareTrainingAndTestTree("pt > 5 && pt <= 10 && abs(eta) > 1.479",
                                        "pt > 5 && pt <= 10 && abs(eta) > 1.479",
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:!V" );
  }
  if (Option == 3) {
    factory->PrepareTrainingAndTestTree("pt > 10 && pt <= 35 && abs(eta) <= 0.8",
                                        "pt > 10 && pt <= 35 && abs(eta) <= 0.8",
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:!V" );
  }
  if (Option == 4) {
    factory->PrepareTrainingAndTestTree("pt > 10 && pt <= 35 && abs(eta) > 0.8 && abs(eta) <= 1.479",
                                        "pt > 10 && pt <= 35 && abs(eta) > 0.8 && abs(eta) <= 1.479",
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:!V" );
  }
  if (Option == 5) {
    factory->PrepareTrainingAndTestTree("pt > 10 && pt <= 35 && abs(eta) > 1.479",
                                        "pt > 10 && pt <= 35 && abs(eta) > 1.479",
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:!V" );
  }


  //*****************************************************************
  //V0: Only RelPt Rings with size 0.10
  //*****************************************************************
  if (versionLabel == "" || versionLabel == "V0" || versionLabel == "V1" ) {

    if (Option >= 0) {
//       factory->BookMethod( TMVA::Types::kBDT, "BDTG",
//                            "!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=5:UseNvars=4:PruneStrength=0.0:PruneMethod=nopruning:MaxDepth=6" );
        factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
                             "H:!V:CreateMVAPdfs:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=30:NSmoothBkg[0]=30:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=40:Nbins=200" );
    }
  }

  
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  
  outputFile->Close();

  delete factory;
}
