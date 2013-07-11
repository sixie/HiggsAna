//root -l -q -b $CMSSW_BASE/src/HiggsAna/BAMBUNtupler/macros/runAllNtupler2012.C+\(\"0000\",\"noskim\",\"r12a-del-pr-v1\",\"cern/filefi/028\",\"/afs/cern.ch/user/s/sixie/catalog\",\"AllNtupler\",10000,-1\) 
//root -l -q -b $CMSSW_BASE/src/HiggsAna/BAMBUNtupler/macros/runAllNtupler2012.C+\(\"0000\",\"noskim\",\"s12-zjets-m50-v5\",\"t2mit/filefi/026\",\"/net/hisrv0001/home/sixie/catalog\",\"AllNtupler\",10000,6\) 

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/PDFProducerMod.h"
#include "MitPhysics/Mods/interface/HKFactorProducer.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/CaloMetCorrectionMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PFTauIDMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/PFTauCleaningMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/PartonFlavorHistoryMod.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/MetCol.h" 
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/CaloMetCol.h"
#include "MitPhysics/SelMods/interface/GenericSelMod.h"
#include "HiggsAna/BAMBUNtupler/interface/HwwNtuplerMod.hh"

#endif

using namespace mithep;

//==================================================================================================
/*
 * Triggers of interest
 *

  kHLT_Mu9                                   = 0x0000001,
  kHLT_Mu11                                  = 0x0000002,
  kHLT_Mu13_v1                               = 0x0000004,
  kHLT_Mu15_v1                               = 0x0000008,
  kHLT_DoubleMu5_v1                          = 0x0000010,
  kHLT_Jet15U                                = 0x0000020,
  kHLT_Jet30U                                = 0x0000040,
  kHLT_Jet50U                                = 0x0000080,
  kHLT_Photon10_L1R                          = 0x0000100,
  kHLT_Photon10_Cleaned_L1R                  = 0x0000200,
  kHLT_Photon15_L1R                          = 0x0000400,
  kHLT_Photon15_Cleaned_L1R                  = 0x0000800,
  kHLT_Photon20_L1R                          = 0x0001000,
  kHLT_Photon20_Cleaned_L1R                  = 0x0002000,
  kHLT_Photon30_Cleaned_L1R                  = 0x0004000,
  kHLT_Ele15_SW_L1R                          = 0x0008000,
  kHLT_Ele15_LW_L1R                          = 0x0010000,
  kHLT_Ele15_SW_CaloEleId_L1R                = 0x0020000,
  kHLT_Ele17_SW_L1R                          = 0x0040000,
  kHLT_Ele17_SW_CaloEleId_L1R                = 0x0080000,       
  kHLT_Ele17_SW_TightEleId_L1R               = 0x0100000,
  kHLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1 = 0x0200000,
  kHLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1  = 0x0400000,
  kHLT_Ele17_SW_TighterEleIdIsol_L1R_v2      = 0x0800000,
  kHLT_DoubleEle10_SW_L1R                    = 0x1000000,
  kHLT_DoubleEle15_SW_L1R_v1                 = 0x2000000

  */
    
//==================================================================================================
/*
 * Run on a BAMBU fileset
 *
 * Example usage:
 *   root -l -q -b runZeeNtupler.C+\(\"0000\",\"p10-zee-v26\",\"cern/filler/014a\",\"/home/ceballos/catalog\",1,0,1,-1,0,1\)
 *
 * Output file name has standard format: <dataset>_<fileset>_ntuple.root
 *
 */
void runAllNtupler2012(
  const char *fileset  = "",
  const char *skim         = "noskim",
  const char *dataset    = "s8-ttbar-id9",
  const char *book       = "mit/filler/006",
  const char *catalogDir = "/home/mitprod/catalog",
  const char *outputName = "HwwHiggsNtupleMaker",
  int   nEvents          = -1,
  int   sampleID         = -1
  )
{
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 3;

  //******************************************************************
  //Set up Options
  //******************************************************************
  bool doHWWNtupler       = false;
  bool doSmurfNtupler     = false;
  bool doFourLeptonNtupler = true;
  bool usePDFProducer     = false;
  string pdfSetName = "";

  bool isData             = false;
  bool isDataMuonElectron = false;
  bool isDataDMuon        = false;
  bool isDataSMuon        = false;
  bool isDataDElectron    = false;
  bool isDataSElectron    = false;
  bool isDataSPhoton      = false;
  Bool_t isEmbeddedSample = kFALSE;
  bool applyMllGenCut     = false;
  bool applyZVetoGenCut   = false;

  double fIntRadius = 0.0;
  double ptJetCut = 30.0;
  double etaJetCut = 5.0;

  int processid = -999999999;
  TString MCType = "kMCTypeUndef";
  Bool_t applyPartonFlavorFilter = kFALSE;
  Bool_t applyISRFilter          = kFALSE;
  Bool_t applyVVFilter           = kFALSE;
  Int_t fakeRatePredictionType = 0;
  Bool_t useSelectGenLeptons   = kTRUE;
  Bool_t isPhotonControlSample = kFALSE;

  int fDecay = sampleID; 
  if (sampleID < 0) fDecay = (-1) * (abs(sampleID) % 10);
  if (sampleID == 10333) fDecay = 3;
  
  Int_t runSkim = 0; if (sampleID <= -1 && sampleID >= -9) runSkim = 1; 
  if (sampleID <= -11 && sampleID >= -19) runSkim = 2; 
  if (sampleID <= -21 && sampleID >= -29) runSkim = 3; 
  if (sampleID <= -31 && sampleID >= -39) runSkim = 4; 

  if (fDecay < 0) {
    isData = true;
    if (fDecay == -1) isDataDMuon = true;
    if (fDecay == -2) isDataDElectron = true;
    if (fDecay == -3) isDataMuonElectron = true;
    if (fDecay == -4) isDataSMuon = true;
    if (fDecay == -5) isDataSElectron = true;
    if (fDecay == -6) isDataSPhoton = true;
    if (fDecay == -7) isEmbeddedSample = true;
  }

  if (fDecay == 14) applyVVFilter = kTRUE;
  if (fDecay > 20000) {
    fakeRatePredictionType = 1;
    useSelectGenLeptons    = kFALSE;
  }
  if (isDataSPhoton || fDecay == 10666 || fDecay == 11666) {
    isPhotonControlSample = kTRUE;
    doHWWNtupler = kFALSE;
  }

  if (isEmbeddedSample) {
    fIntRadius          = 0.05;
    doHWWNtupler        = false;
    doFourLeptonNtupler = false;
  }

  cout << "Summarize Run Options: " << "fDecay == " << fDecay << " "
       << "runSkim == " << runSkim << " "
       << endl;



  //******************************************************************
  //Modules
  //******************************************************************

  // Generator info
  GeneratorMod *GeneratorMod1 = new GeneratorMod;
  GeneratorMod1->SetPrintDebug(kFALSE);
  GeneratorMod1->SetPtLeptonMin(0.0);
  GeneratorMod1->SetEtaLeptonMax(2.7);
  GeneratorMod1->SetPtPhotonMin(0.0);
  GeneratorMod1->SetEtaPhotonMax(2.7);
  GeneratorMod1->SetPtRadPhotonMin(5.0);
  GeneratorMod1->SetEtaRadPhotonMax(2.7);
  GeneratorMod1->SetIsData(isData);
  GeneratorMod1->SetFillHist(!isData);
  if(applyMllGenCut == kTRUE){
    GeneratorMod1->SetPdgIdCut(23);
    GeneratorMod1->SetMassMaxCut(50.);
  }
  else if(applyZVetoGenCut == kTRUE){
    GeneratorMod1->SetPdgIdCut(23);
    GeneratorMod1->SetMassMinCut(20000.);
    GeneratorMod1->SetMassMaxCut(20000.);
  }
  GeneratorMod1->SetApplyISRFilter(applyISRFilter);
  GeneratorMod1->SetApplyVVFilter(applyVVFilter);
  GeneratorMod1->SetAllowWWEvents(kTRUE);
  GeneratorMod1->SetAllowWZEvents(kFALSE);
  GeneratorMod1->SetAllowZZEvents(kFALSE);

  PartonFlavorHistoryMod *PartonFlavorHistoryMod1 = new PartonFlavorHistoryMod;
  PartonFlavorHistoryMod1->SetMCSampleType(MCType);
  PartonFlavorHistoryMod1->SetApplyPartonFlavorFilter(applyPartonFlavorFilter);

  // HLT info
  HLTMod *hltmod = new HLTMod;
  hltmod->AddTrigger("HLT_Mu15_v2");
  hltmod->AddTrigger("HLT_IsoMu17_v5");
  hltmod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2");
  hltmod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2");
  hltmod->AddTrigger("HLT_DoubleMu7_v1");
  hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v2");
  hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v2"); 
  hltmod->AddTrigger("!HLT_Mu15_v2");
  hltmod->SetAbortIfNotAccepted(kFALSE);
  hltmod->SetTrigObjsName("myhltobjs");


  //------------------------------------------------------------------------------------------------
  // Run RunLumiSelectionMod
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSelection = new RunLumiSelectionMod;      
  runLumiSelection->SetAcceptMC(!isData);
  runLumiSelection->SetAcceptAll(kTRUE);
  runLumiSelection->SetAbortIfNotAccepted(kFALSE);

  //------------------------------------------------------------------------------------------------
  // PV filter selection
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof(4);
  goodPVFilterMod->SetMaxAbsZ(24.0);
  goodPVFilterMod->SetMaxRho(2.0);
  if (sampleID == 333 || sampleID > 10000) { goodPVFilterMod->SetVertexesName("DAPrimaryVertexes"); cout << "use DAPrimaryVertexes\n";}
  else goodPVFilterMod->SetVertexesName("PrimaryVertexes");

  Bool_t isFastSim = kFALSE;


  //***********************************************************************************8
  //Lepton Selection
  //***********************************************************************************8

  // Object ID and Cleaning Sequence
  MuonIDMod *muonID1 = new MuonIDMod;
  muonID1->SetPrintMVADebugInfo(kFALSE);
  muonID1->SetClassType("GlobalTracker");
  muonID1->SetIDType("WWMuIdV3");
  muonID1->SetIsoType("PFIso");
  muonID1->SetApplyD0Cut(kFALSE);
  muonID1->SetApplyDZCut(kTRUE);
  muonID1->SetWhichVertex(0);
  muonID1->SetIntRadius(fIntRadius);

  ElectronIDMod *electronID1 = new ElectronIDMod;
  electronID1->SetPrintMVADebugInfo(kFALSE);
  electronID1->SetIDType("VBTFWorkingPoint80Id");
  electronID1->SetIsoType("PFIso");
  electronID1->SetChargeFilter(kFALSE);
  electronID1->SetApplyD0Cut(kTRUE);
  electronID1->SetApplyDZCut(kTRUE);
  electronID1->SetWhichVertex(0);
  electronID1->SetApplyConversionFilterType2(kFALSE);
  electronID1->SetApplyConversionFilterType1(kTRUE);
  electronID1->SetNExpectedHitsInnerCut(0);
  electronID1->SetIntRadius(fIntRadius);


  PhotonIDMod *photonIDMod1 = new PhotonIDMod;
  photonIDMod1->SetIsoType("MITPUCorrected");
  photonIDMod1->SetHadOverEmMax(0.05);
  photonIDMod1->SetApplyPixelSeed(kFALSE);
  photonIDMod1->SetApplyElectronVetoConvRecovery(kTRUE);
  photonIDMod1->SetApplyConversionId(kTRUE);

  PFTauIDMod *pftauIDMod1 = new PFTauIDMod;
  pftauIDMod1->SetPFTausName("HPSTaus");
  pftauIDMod1->SetIsHPSSel(kFALSE);


  //***********************************************************************************8
  //Jet Selection
  //***********************************************************************************8

  const char *jetInput1 = "AKt5PFJets";
  PublisherMod<PFJet,Jet> *pubJet1 = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet1->SetInputName(jetInput1);
  pubJet1->SetOutputName(Form("Pub%s",jetInput1));

  JetCorrectionMod *jetCorr1_ntuple = new JetCorrectionMod;
  jetCorr1_ntuple->AddCorrectionFromFile("/afs/cern.ch/user/s/sixie/CMSSW_analysis/src/MitPhysics/data/START52_V9_L1FastJet_AK5PF.txt.txt"); 
  jetCorr1_ntuple->AddCorrectionFromFile("/afs/cern.ch/user/s/sixie/CMSSW_analysis/src/MitPhysics/data/START52_V9_L2Relative_AK5PF.txt"); 
  jetCorr1_ntuple->AddCorrectionFromFile("/afs/cern.ch/user/s/sixie/CMSSW_analysis/src/MitPhysics/data/START52_V9_L3Absolute_AK5PF.txt");
  if(isData == true){
    jetCorr1_ntuple->AddCorrectionFromFile("/afs/cern.ch/user/s/sixie/CMSSW_analysis/src/MitPhysics/data/START52_V9_L2L3Residual_AK5PF.txt");  
  }
  jetCorr1_ntuple->SetInputName(pubJet1->GetOutputName());
  jetCorr1_ntuple->SetCorrectedName("CorrectedJets_ntuple");

  JetIDMod *theJetID1_ntuple = new JetIDMod;
  theJetID1_ntuple->SetInputName(jetCorr1_ntuple->GetOutputName());
  theJetID1_ntuple->SetPtCut(20.0);
  theJetID1_ntuple->SetEtaMaxCut(etaJetCut);
  theJetID1_ntuple->SetJetEEMFractionMinCut(0.00);
  theJetID1_ntuple->SetOutputName("GoodJets_ntuple");
  theJetID1_ntuple->SetApplyBetaCut(kFALSE);



  //***********************************************************************************8
  //Cleaning
  //***********************************************************************************8
  ElectronCleaningMod *electronCleaning1 = new ElectronCleaningMod;
  PhotonCleaningMod *photonCleaningMod1 = new PhotonCleaningMod;
  PFTauCleaningMod *pftauCleaningMod1 = new PFTauCleaningMod;

  JetCleaningMod *theJetCleaning1_ntuple = new JetCleaningMod;
  theJetCleaning1_ntuple->SetGoodJetsName("CorrectedJets_ntuple");
  theJetCleaning1_ntuple->SetCleanJetsName("CleanJets_ntuple");
  if (isPhotonControlSample) {
    theJetCleaning1_ntuple->SetApplyPhotonRemoval(kTRUE);
  }

  MergeLeptonsMod *merger1 = new MergeLeptonsMod;
  merger1->SetMuonsName(muonID1->GetOutputName());
  merger1->SetElectronsName(electronCleaning1->GetOutputName());

  JetIDMod *theJetID2_ntuple = new JetIDMod;
  theJetID2_ntuple->SetInputName(jetCorr1_ntuple->GetOutputName());
  theJetID2_ntuple->SetPtCut(0.0);
  theJetID2_ntuple->SetEtaMaxCut(etaJetCut);
  theJetID2_ntuple->SetJetEEMFractionMinCut(0.00);
  theJetID2_ntuple->SetOutputName("GoodJetsNoPtCut_ntuple");
  theJetID2_ntuple->SetApplyBetaCut(kFALSE);

  JetCleaningMod *theJetCleaning2_ntuple = new JetCleaningMod;
  theJetCleaning2_ntuple->SetGoodJetsName(jetCorr1_ntuple->GetOutputName());
  theJetCleaning2_ntuple->SetCleanJetsName("CleanJetsNoPtCut_ntuple");
  if (isPhotonControlSample) {
    theJetCleaning2_ntuple->SetApplyPhotonRemoval(kTRUE);
  }

  //***********************************************************************************8
  //MET
  //***********************************************************************************8
  const char *metInput = "TCMet";
  PublisherMod<Met,Met> *pubMet = new PublisherMod<Met,Met>("MetPub");
  pubMet->SetInputName(metInput);
  pubMet->SetOutputName(Form("Pub%s",metInput));

  const char *metPFInput = "PFMet";
  PublisherMod<PFMet,Met> *pubPFMet = new PublisherMod<PFMet,Met>("MetPFPub");
  pubPFMet->SetInputName(metPFInput);
  pubPFMet->SetOutputName(Form("Pub%s",metPFInput));

  //***********************************************************************************8
  //Fakeable Objects Definition
  //***********************************************************************************8

  // Lepton ID with loose requirements
  MuonIDMod *muonIDFakeable = new MuonIDMod;
  muonIDFakeable->SetClassType("GlobalTracker");
  muonIDFakeable->SetIDType("WWMuIdV3");
  muonIDFakeable->SetIsoType("PFIso");
  muonIDFakeable->SetApplyD0Cut(kTRUE);
  muonIDFakeable->SetApplyDZCut(kTRUE);
  muonIDFakeable->SetD0Cut(0.20);
  muonIDFakeable->SetPFIsoCut(0.40);
  muonIDFakeable->SetWhichVertex(0);
  muonIDFakeable->SetIntRadius(fIntRadius);
  muonIDFakeable->SetCleanMuonsName("CleanMuonsFakeable");

  ElectronIDMod *electronIDFakeable = new ElectronIDMod;
  electronIDFakeable->SetIDType("VBTFWorkingPointFakeableId");
  electronIDFakeable->SetIsoType("TrackJura");
  electronIDFakeable->SetTrackIsoCut(0.2);
  electronIDFakeable->SetEcalJurIsoCut(0.2);
  electronIDFakeable->SetHcalIsoCut(0.2);
  electronIDFakeable->SetApplyConversionFilterType1(kTRUE);
  electronIDFakeable->SetApplyConversionFilterType2(kFALSE);
  electronIDFakeable->SetChargeFilter(kFALSE);
  electronIDFakeable->SetApplyD0Cut(kTRUE);
  electronIDFakeable->SetApplyDZCut(kTRUE);
  electronIDFakeable->SetNExpectedHitsInnerCut(0);
  electronIDFakeable->SetD0Cut(0.02);
  electronIDFakeable->SetWhichVertex(0);
  electronIDFakeable->SetIntRadius(fIntRadius);
  electronIDFakeable->SetGoodElectronsName("GoodElectronsFakeable");

  ElectronCleaningMod *electronCleaningFakeable = new ElectronCleaningMod;
  electronCleaningFakeable->SetCleanMuonsName(muonIDFakeable->GetOutputName());
  electronCleaningFakeable->SetGoodElectronsName(electronIDFakeable->GetOutputName());
  electronCleaningFakeable->SetCleanElectronsName("CleanElectronsFakeable");

  MergeLeptonsMod *mergerFakeable = new MergeLeptonsMod;
  mergerFakeable->SetMuonsName(muonIDFakeable->GetOutputName());
  mergerFakeable->SetElectronsName(electronCleaningFakeable->GetOutputName());
  mergerFakeable->SetMergedName("MergedLeptonsFakeable");

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile =  TString("./") + TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  printf("\nRoot output: %s\n\n",rootFile.Data());  



  //------------------------------------------------------------------------------------------------
  //
  // HZZ4l Ntupler
  //
  //------------------------------------------------------------------------------------------------
  TString rootFileHZZ4lNtuple = TString("./");
  rootFileHZZ4lNtuple += TString(outputName);
  rootFileHZZ4lNtuple += TString("_HZZ4lNtuple_") + TString(dataset) + TString("_") + TString(skim); 
  rootFileHZZ4lNtuple += TString("_") + TString(fileset);
  rootFileHZZ4lNtuple += TString(".root");

  HwwNtuplerMod *hzz4lNtuplerMod = new HwwNtuplerMod;
  hzz4lNtuplerMod->SetOutputName(rootFileHZZ4lNtuple.Data());          // output ntuple file name
  hzz4lNtuplerMod->SetUseGen(!isData);              // look at generator information (must set to kFALSE if MCParticle collection do not exist)
  hzz4lNtuplerMod->SetSkipIfHLTFail(kFALSE);  // skip to next event if no HLT accept
  hzz4lNtuplerMod->SetFSRMode(kFALSE);
  hzz4lNtuplerMod->SetMuonPtMin(3);
  hzz4lNtuplerMod->SetMuonPtMax(99999999);
  hzz4lNtuplerMod->SetMuonEtaMin(-3);
  hzz4lNtuplerMod->SetMuonEtaMax(3);
  hzz4lNtuplerMod->SetElePtMin(5);
  hzz4lNtuplerMod->SetElePtMax(99999999);
  hzz4lNtuplerMod->SetEleEtaMin(-3);
  hzz4lNtuplerMod->SetEleEtaMax(3);
  hzz4lNtuplerMod->SetCleanJetsName(theJetCleaning1_ntuple->GetOutputName());
  hzz4lNtuplerMod->SetCleanJetsNoPtCutName(jetCorr1_ntuple->GetOutputName());
  hzz4lNtuplerMod->SetJetPtMin(7);
  hzz4lNtuplerMod->SetComputePDFWeights(usePDFProducer);
  hzz4lNtuplerMod->SetPDFName(pdfSetName.c_str());


  //***********************************************************************************************************
  //Single Muons
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu24_eta2p1_v11",                        kHLT_IsoMu24, kHLTObject_IsoMu24 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu24_eta2p1_v12",                        kHLT_IsoMu24, kHLTObject_IsoMu24 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu24_eta2p1_v13",                        kHLT_IsoMu24, kHLTObject_IsoMu24 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu24_eta2p1_v14",                        kHLT_IsoMu24, kHLTObject_IsoMu24 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu24_eta2p1_v15",                        kHLT_IsoMu24, kHLTObject_IsoMu24 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu30_eta2p1_v11",                        kHLT_IsoMu24, kHLTObject_IsoMu30 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f30QL3crIsoFiltered10" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu30_eta2p1_v12",                        kHLT_IsoMu24, kHLTObject_IsoMu30 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f30QL3crIsoFiltered10" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu30_eta2p1_v13",                        kHLT_IsoMu24, kHLTObject_IsoMu30 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f30QL3crIsoRhoFiltered0p15" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu30_eta2p1_v14",                        kHLT_IsoMu24, kHLTObject_IsoMu30 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f30QL3crIsoRhoFiltered0p15" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu30_eta2p1_v15",                        kHLT_IsoMu24, kHLTObject_IsoMu30 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f30QL3crIsoRhoFiltered0p15" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu34_eta2p1_v9",                         kHLT_IsoMu24, kHLTObject_IsoMu34 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f34QL3crIsoFiltered10" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu34_eta2p1_v10",                        kHLT_IsoMu24, kHLTObject_IsoMu34 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f34QL3crIsoFiltered10" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu34_eta2p1_v11",                        kHLT_IsoMu24, kHLTObject_IsoMu34 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f34QL3crIsoRhoFiltered0p15" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu34_eta2p1_v12",                        kHLT_IsoMu24, kHLTObject_IsoMu34 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f34QL3crIsoRhoFiltered0p15" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu34_eta2p1_v13",                        kHLT_IsoMu24, kHLTObject_IsoMu34 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f34QL3crIsoRhoFiltered0p15" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu40_eta2p1_v6",                         kHLT_IsoMu24, kHLTObject_IsoMu40 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f40QL3crIsoFiltered10" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu40_eta2p1_v7",                         kHLT_IsoMu24, kHLTObject_IsoMu40 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f40QL3crIsoFiltered10" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu40_eta2p1_v8",                         kHLT_IsoMu24, kHLTObject_IsoMu40 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f40QL3crIsoRhoFiltered0p15" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu40_eta2p1_v9",                         kHLT_IsoMu24, kHLTObject_IsoMu40 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f40QL3crIsoRhoFiltered0p15" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu40_eta2p1_v10",                         kHLT_IsoMu24, kHLTObject_IsoMu40 , "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f40QL3crIsoRhoFiltered0p15" );

  //these were added after Run2011B - i guess after the issues with endcap trigger muons were fixed
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu24_v15",                        kHLT_IsoMu24, kHLTObject_IsoMu24 , "hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu24_v16",                        kHLT_IsoMu24, kHLTObject_IsoMu24 , "hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu24_v17",                        kHLT_IsoMu24, kHLTObject_IsoMu24 , "hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu30_v9",                         kHLT_IsoMu24, kHLTObject_IsoMu24 , "hltL3crIsoL1sMu16L1f0L2f16QL3f30QL3crIsoRhoFiltered0p15" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu30_v10",                        kHLT_IsoMu24, kHLTObject_IsoMu24 , "hltL3crIsoL1sMu16L1f0L2f16QL3f30QL3crIsoRhoFiltered0p15" );
  hzz4lNtuplerMod->AddTrigger("HLT_IsoMu30_v11",                        kHLT_IsoMu24, kHLTObject_IsoMu24 , "hltL3crIsoL1sMu16L1f0L2f16QL3f30QL3crIsoRhoFiltered0p15" );

  //
  //***********************************************************************************************************


  //***********************************************************************************************************
 //Single Electron Triggers
  hzz4lNtuplerMod->AddTrigger("HLT_Ele27_WP80_v8",                             kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, kHLTObject_Ele27_WP80, "hltEle27WP80TrackIsoFilter");
  hzz4lNtuplerMod->AddTrigger("HLT_Ele27_WP80_v9",                             kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, kHLTObject_Ele27_WP80, "hltEle27WP80TrackIsoFilter");
  hzz4lNtuplerMod->AddTrigger("HLT_Ele27_WP80_v10",                            kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, kHLTObject_Ele27_WP80, "hltEle27WP80TrackIsoFilter");
  hzz4lNtuplerMod->AddTrigger("HLT_Ele27_WP80_v11",                            kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, kHLTObject_Ele27_WP80, "hltEle27WP80TrackIsoFilter");

  //W-specific trigger 
  hzz4lNtuplerMod->AddTrigger("HLT_Ele27_WP80_PFMET_MT50_v2",                  kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, kHLTObject_Ele27_WP80, "hltEle27WP80TrackIsoFilter");
  hzz4lNtuplerMod->AddTrigger("HLT_Ele27_WP80_PFMET_MT50_v3",                  kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, kHLTObject_Ele27_WP80, "hltEle27WP80TrackIsoFilter");
  hzz4lNtuplerMod->AddTrigger("HLT_Ele27_WP80_PFMET_MT50_v4",                  kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, kHLTObject_Ele27_WP80, "hltEle27WP80TrackIsoFilter");
  hzz4lNtuplerMod->AddTrigger("HLT_Ele27_WP80_PFMET_MT50_v5",                  kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, kHLTObject_Ele27_WP80, "hltEle27WP80TrackIsoFilter");
  hzz4lNtuplerMod->AddTrigger("HLT_Ele27_WP80_PFMET_MT50_v6",                  kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, kHLTObject_Ele27_WP80, "hltEle27WP80TrackIsoFilter");
  //
  //***********************************************************************************************************



  //***********************************************************************************************************
  //Main Dielectron Triggers
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15", kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16", kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17", kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18", kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter" );


  //
  //***********************************************************************************************************

  //***********************************************************************************************************
  //Main Dimuon Triggers
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_Mu8_v16",                          kHLT_Mu17_Mu8,   kHLTObject_Mu17, "hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17", kHLTObject_Mu8, "hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8");
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_Mu8_v17",                          kHLT_Mu17_Mu8,   kHLTObject_Mu17, "hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17", kHLTObject_Mu8, "hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8");
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_Mu8_v18",                          kHLT_Mu17_Mu8,   kHLTObject_Mu17, "hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17", kHLTObject_Mu8, "hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8");
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_Mu8_v19",                          kHLT_Mu17_Mu8,   kHLTObject_Mu17, "hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17", kHLTObject_Mu8, "hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8");
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_Mu8_v20",                          kHLT_Mu17_Mu8,   kHLTObject_Mu17, "hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17", kHLTObject_Mu8, "hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8");
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_Mu8_v21",                          kHLT_Mu17_Mu8,   kHLTObject_Mu17, "hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17", kHLTObject_Mu8, "hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8");
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_TkMu8_v9",                         kHLT_Mu17_TkMu8, kHLTObject_Mu17, "hltL3fL1sMu10MuOpenL1f0L2f10L3Filtered17", kHLTObject_TkMu8, "hltDiMuonGlbFiltered17TrkFiltered8");
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_TkMu8_v10",                         kHLT_Mu17_TkMu8, kHLTObject_Mu17, "hltL3fL1sMu10MuOpenL1f0L2f10L3Filtered17", kHLTObject_TkMu8, "hltDiMuonGlbFiltered17TrkFiltered8");
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_TkMu8_v11",                         kHLT_Mu17_TkMu8, kHLTObject_Mu17, "hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17", kHLTObject_TkMu8, "hltDiMuonGlbFiltered17TrkFiltered8");
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_TkMu8_v12",                         kHLT_Mu17_TkMu8, kHLTObject_Mu17, "hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17", kHLTObject_TkMu8, "hltDiMuonGlbFiltered17TrkFiltered8");
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_TkMu8_v13",                         kHLT_Mu17_TkMu8, kHLTObject_Mu17, "hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17", kHLTObject_TkMu8, "hltDiMuonGlbFiltered17TrkFiltered8");
  //
  //***********************************************************************************************************


  //***********************************************************************************************************
  //Main EMu Triggers
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4",       kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Mu17, "hltL1Mu12EG7L3MuFiltered17", kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter"  );
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5",       kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Mu17, "hltL1Mu12EG7L3MuFiltered17", kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter"  );
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",       kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Mu17, "hltL1Mu12EG7L3MuFiltered17", kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter"  );
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",       kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Mu17, "hltL1Mu12EG7L3MuFiltered17", kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter"  );
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",       kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Mu17, "hltL1Mu12EG7L3MuFiltered17", kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter"  );


  hzz4lNtuplerMod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4",       kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Mu8, "hltL1MuOpenEG12L3Filtered8", kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter"  );
  hzz4lNtuplerMod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5",       kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Mu8, "hltL1MuOpenEG12L3Filtered8", kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter"  );
  hzz4lNtuplerMod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",       kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Mu8, "hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8", kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter"  );
  hzz4lNtuplerMod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",       kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Mu8, "hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8", kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter"  );
  hzz4lNtuplerMod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",       kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Mu8, "hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8", kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter"  );
  //
  //***********************************************************************************************************










  //***********************************************************************************************************
  //T&P triggers

  hzz4lNtuplerMod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v3",                   kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50, kHLTObject_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT, "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v4",                   kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50, kHLTObject_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT, "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v5",                   kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50, kHLTObject_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT, "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v6",                   kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50, kHLTObject_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT, "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter" );

  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v3",                   kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter"  );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v4",                   kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter"  );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v5",                   kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter"  );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v6",                   kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter"  );

  hzz4lNtuplerMod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v3",                   kHLT_Ele32_CaloIdL_CaloIsoVL_SC17, kHLTObject_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter"  );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v4",                   kHLT_Ele32_CaloIdL_CaloIsoVL_SC17, kHLTObject_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter"  );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v5",                   kHLT_Ele32_CaloIdL_CaloIsoVL_SC17, kHLTObject_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter"  );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v6",                   kHLT_Ele32_CaloIdL_CaloIsoVL_SC17, kHLTObject_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter"  );

  //
  //***********************************************************************************************************







  //***********************************************************************************************************
  //Fake Rate triggers
  hzz4lNtuplerMod->AddTrigger("HLT_Mu8_v16",                                kHLT_Mu8 ,  kHLTObject_Mu8 , "hltL3fL1sMu3L3Filtered8"  );
  hzz4lNtuplerMod->AddTrigger("HLT_Mu8_v17",                                kHLT_Mu8 ,  kHLTObject_Mu8 , "hltL3fL1sMu3L3Filtered8"  );
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_v3",                                kHLT_Mu17,  kHLTObject_Mu17, "hltL3fL1sMu12L3Filtered17");
  hzz4lNtuplerMod->AddTrigger("HLT_Mu17_v4",                                kHLT_Mu17,  kHLTObject_Mu17, "hltL3fL1sMu12L3Filtered17");


  hzz4lNtuplerMod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v14",                            kHLT_Ele8_CaloIdL_CaloIsoVL,  kHLTObject_Ele8_CaloIdL_CaloIsoVL , "hltEle8CaloIdLCaloIsoVLPixelMatchFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v15",                            kHLT_Ele8_CaloIdL_CaloIsoVL,  kHLTObject_Ele8_CaloIdL_CaloIsoVL , "hltEle8CaloIdLCaloIsoVLPixelMatchFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v16",                            kHLT_Ele8_CaloIdL_CaloIsoVL,  kHLTObject_Ele8_CaloIdL_CaloIsoVL , "hltEle8CaloIdLCaloIsoVLPixelMatchFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v17",                            kHLT_Ele8_CaloIdL_CaloIsoVL,  kHLTObject_Ele8_CaloIdL_CaloIsoVL , "hltEle8CaloIdLCaloIsoVLPixelMatchFilter" );



  hzz4lNtuplerMod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v3",      kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle8TightIdLooseIsoTrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v4",      kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle8TightIdLooseIsoTrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v5",      kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle8TightIdLooseIsoTrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v6",      kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle8TightIdLooseIsoTrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v7",      kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle8TightIdLooseIsoTrackIsoFilter" );

  hzz4lNtuplerMod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12",           kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle8TightIdLooseIsoTrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13",           kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle8TightIdLooseIsoTrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14",           kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle8TightIdLooseIsoTrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15",           kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle8TightIdLooseIsoTrackIsoFilter" );

  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v14",                           kHLT_Ele17_CaloIdL_CaloIsoVL, kHLTObject_Ele17_CaloIdL_CaloIsoVL , "hltEle17CaloIdLCaloIsoVLPixelMatchFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v15",                           kHLT_Ele17_CaloIdL_CaloIsoVL, kHLTObject_Ele17_CaloIdL_CaloIsoVL , "hltEle17CaloIdLCaloIsoVLPixelMatchFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v16",                           kHLT_Ele17_CaloIdL_CaloIsoVL, kHLTObject_Ele17_CaloIdL_CaloIsoVL , "hltEle17CaloIdLCaloIsoVLPixelMatchFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v17",                           kHLT_Ele17_CaloIdL_CaloIsoVL, kHLTObject_Ele17_CaloIdL_CaloIsoVL , "hltEle17CaloIdLCaloIsoVLPixelMatchFilter" );

  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v3",     kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v4",     kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v5",     kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v6",     kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v7",     kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter" );

  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3",           kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4",           kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5",           kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter" );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",           kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter" );






  //
  //***********************************************************************************************************

  //Add L1Triggers
  hzz4lNtuplerMod->AddL1Trigger("L1_SingleEG5",  kL1_SingleEG5 );
  hzz4lNtuplerMod->AddL1Trigger("L1_SingleEG12", kL1_SingleEG12  );                  
  hzz4lNtuplerMod->AddL1Trigger("L1_SingleEG20", kL1_SingleEG20  );                  
  hzz4lNtuplerMod->AddL1Trigger("L1_SingleEG30", kL1_SingleEG30  );                  
  hzz4lNtuplerMod->AddL1Trigger("L1_SingleMu3",  kL1_SingleMu3 );                     
  hzz4lNtuplerMod->AddL1Trigger("L1_SingleMu7 ", kL1_SingleMu7 );                    
  hzz4lNtuplerMod->AddL1Trigger("L1_SingleMu10", kL1_SingleMu10 );                    
  hzz4lNtuplerMod->AddL1Trigger("L1_SingleMu12", kL1_SingleMu12 );                   
  hzz4lNtuplerMod->AddL1Trigger("L1_SingleMu20", kL1_SingleMu20 );
  hzz4lNtuplerMod->AddL1Trigger("L1_SingleMu25", kL1_SingleMu25 );
  hzz4lNtuplerMod->AddL1Trigger("L1_MuOpen_EG12",kL1_MuOpen_EG12 );                 
  hzz4lNtuplerMod->AddL1Trigger("L1_MuOpen_EG5", kL1_MuOpen_EG5  );                   
  hzz4lNtuplerMod->AddL1Trigger("L1_Mu3_EG5",    kL1_Mu3_EG5     );                    
  hzz4lNtuplerMod->AddL1Trigger("L1_Mu5_EG12",   kL1_Mu5_EG12    );                    
  hzz4lNtuplerMod->AddL1Trigger("L1_Mu7_EG5",    kL1_Mu7_EG5 );                    
  hzz4lNtuplerMod->AddL1Trigger("L1_Mu12_EG5",   kL1_Mu12_EG5 );                      
  hzz4lNtuplerMod->AddL1Trigger("L1_DoubleMu0",  kL1_DoubleMu0 );                     
  hzz4lNtuplerMod->AddL1Trigger("L1_DoubleMu5",  kL1_DoubleMu5 );                    
  hzz4lNtuplerMod->AddL1Trigger("L1_DoubleEG3",  kL1_DoubleEG3  );                    
  hzz4lNtuplerMod->AddL1Trigger("L1_DoubleEG5",  kL1_DoubleEG5  );                    
  hzz4lNtuplerMod->AddL1Trigger("L1_DoubleEG_12_5", kL1_DoubleEG_12_5  );            
  hzz4lNtuplerMod->AddL1Trigger("L1_Mu3_Jet20_Central", kL1_Mu3_Jet20_Central  );    
  hzz4lNtuplerMod->AddL1Trigger("L1_EG5_Jet36_deltaPhi1",  kL1_EG5_Jet36_deltaPhi1 );

  hzz4lNtuplerMod->SetPrintHLT(kTRUE); // print HLT table at start of analysis?
  

  //Add 2010 Triggers
  hzz4lNtuplerMod->AddTrigger("HLT_DoubleMu3", kHLT_DoubleMu3 );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele10_LW_EleId_L1R", kHLT_Ele10_LW );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele15_SW_CaloEleId_L1R", kHLT_Ele15_SW );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele15_SW_EleId_L1R"    , kHLT_Ele15_SW );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele15_SW_L1R"          , kHLT_Ele15_SW );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2", kHLT_Ele17_SW );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1", kHLT_Ele17_SW );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1", kHLT_Ele17_SW );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3", kHLT_Ele17_SW );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2", kHLT_Ele17_SW );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v1", kHLT_Ele17_SW );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_SW_TighterEleId_L1R_v1", kHLT_Ele17_SW );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_SW_TightEleIdIsol_L1R_v1", kHLT_Ele17_SW );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_SW_TightEleId_L1R", kHLT_Ele17_SW );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_SW_CaloEleId_L1R", kHLT_Ele17_SW );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_SW_EleId_L1R", kHLT_Ele17_SW );
  hzz4lNtuplerMod->AddTrigger("HLT_Ele17_SW_LooseEleId_L1R", kHLT_Ele17_SW );



  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // Chain modules together
  if(isData == false){
    GeneratorMod1->Add(PartonFlavorHistoryMod1);
    PartonFlavorHistoryMod1->Add(runLumiSelection);
    runLumiSelection->Add(goodPVFilterMod);
  }
  else {
    GeneratorMod1->Add(runLumiSelection);
    runLumiSelection->Add(goodPVFilterMod);
  }

  //Standard Sequence
  goodPVFilterMod->Add(muonID1);
  muonID1->Add(electronID1);
  electronID1->Add(photonIDMod1);
  photonIDMod1->Add(pftauIDMod1);
  pftauIDMod1->Add(pubJet1);
  pubJet1->Add(jetCorr1_ntuple);
  jetCorr1_ntuple->Add(electronCleaning1);
  electronCleaning1->Add(photonCleaningMod1);
  photonCleaningMod1->Add(pftauCleaningMod1);
  pftauCleaningMod1->Add(theJetCleaning1_ntuple);
  theJetCleaning1_ntuple->Add(theJetCleaning2_ntuple);
  theJetCleaning2_ntuple->Add(merger1);
  merger1->Add(pubMet);
  pubMet->Add(pubPFMet);

  //------------------------------------------------------------------------------------------------
  // FourLepton Skim + Ntupler
  //------------------------------------------------------------------------------------------------
  if (doFourLeptonNtupler) {
    pubPFMet->Add(hzz4lNtuplerMod);
  }

  //------------------------------------------------------------------------------------------------
  //
  // setup analysis object
  //
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kFALSE);
  if(nEvents >= 0) 
    ana->SetProcessNEvents(nEvents);

  ana->AddSuperModule(GeneratorMod1);
  ana->SetPrintScale(100);

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  printf("\nRely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n\n",book,dataset,skim,fileset);
  Catalog *c = new Catalog(catalogDir);
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(book,dataset,fileset);
  else 
    d = c->FindDataset(book,skimdataset.Data(),fileset);
  ana->AddDataset(d);

  //sync file
  //ana->AddFile("/data/blue/khahn/bambu/backport/s12-h120zz4l-gf-v9/0CAA68E2-3491-E111-9F03-003048FFD760.root");



  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  ana->SetOutputName(rootFile.Data());
  ana->SetCacheSize(0);


  //
  // run analysis after successful initialisation
  //
  ana->Run(!gROOT->IsBatch());
}

