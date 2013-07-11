//================================================================================================
//
// Select Zee events for scale and resolution studies
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for 4-vector calculations
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "Common/MyTools.hh"        // miscellaneous helper functions

// define structures to read in ntuple
#include "HiggsAna/DataTree/interface/HiggsAnaDefs.hh"
#include "HiggsAna/DataTree/interface/TEventInfo.hh"
#include "HiggsAna/DataTree/interface/TElectron.hh"
#include "HiggsAna/DataTree/interface/TMuon.hh"
#include "HiggsAna/DataTree/interface/TJet.hh"
#include "HiggsAna/DataTree/interface/TGenInfo.hh"
#include "HiggsAna/DataTree/interface/TGenParticle.hh"
#include "HiggsAna/DataTree/interface/TPFCandidate.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// Helper functions for Electron ID selection
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"
#include "HiggsAna/Utils/CommonDefs.hh"
#include "HiggsAna/Utils/LeptonIDCuts.hh"
#include "HiggsAna/HZZ4l/Utils/LeptonSelection.hh"
#include "HiggsAna/HZZ4l/Utils/HZZDefs.hh"
#include "HiggsAna/HZZ4l/Utils/FSRRecovery.hh"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

// Evaluating function for regression
#include "EGamma/EGammaAnalysisTools/interface/ElectronEnergyRegressionEvaluate.h"

// structure for output ntuple
#include "CITCommon/CommonData/interface/ZeeEventTree.h"
#endif


//=== MAIN MACRO ================================================================================================= 

void MakeZeeEventNtuples(const string inputfile,    // input file
		const string outputfile,   // output file
		const Bool_t  matchGen = kFALSE, // match to generator
		const string PUReweightFile = "",
		Int_t PDType = 0,
		Int_t DataEraInput = 2
		) {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  UInt_t DataEra = kDataEra_NONE;
  if (DataEraInput == 1) DataEra = kDataEra_2011_MC;
  if (DataEraInput == 2) DataEra = kDataEra_2012_MC;
  if (DataEraInput == 11) DataEra = kDataEra_2011_Data;
  if (DataEraInput == 12) DataEra = kDataEra_2012_Data;

  //*****************************************************************************************
  //Setup MVA
  //*****************************************************************************************
  EGammaMvaEleEstimator *eleIDMVA = new EGammaMvaEleEstimator();
  vector<string> weightFiles;

  weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat1.weights.xml")).c_str());
  weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat2.weights.xml")).c_str());
  weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat3.weights.xml")).c_str());
  weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat4.weights.xml")).c_str());
  weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat5.weights.xml")).c_str());
  weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat6.weights.xml")).c_str());

  eleIDMVA->initialize( "BDT", EGammaMvaEleEstimator::kNonTrig,
		  kTRUE, weightFiles);


  //********************************************************
  // Pileup Reweighting
  //********************************************************
  TFile *fPUFile = 0;
  if (PUReweightFile != "") {
    fPUFile = TFile::Open(PUReweightFile.c_str());
  }
  TH1D *fhDPU = 0;
  if (fPUFile) fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  if(fhDPU) fhDPU->SetDirectory(0);
  if (fPUFile) delete fPUFile;


  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t applyGoodLumi = kTRUE;
  if (matchGen) applyGoodLumi = kFALSE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("/afs/cern.ch/work/s/sixie/public/HZZ4l/auxiliar/2011/2011Combined.json"); 
  rlrm.AddJSONFile("/afs/cern.ch/work/s/sixie/public/HZZ4l/auxiliar/2012/Cert_190456-196531_8TeV_PromptReco_Collisions12_JSON.txt");


  //********************************************************
  // mass region
  //********************************************************
  Double_t massLo = 40;
  Double_t massHi = 200;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  Double_t nEvents = 0;

  //*****************************************************************************************
  // Set up output ntuple
  //*****************************************************************************************
  TFile *outFile = new TFile(outputfile.c_str(),"RECREATE"); 

  citana::ZeeEventTree *zeeEventTree = new citana::ZeeEventTree;
  zeeEventTree->CreateTree();
  //   zeeEventTree->tree_->SetAutoFlush(0);


  //*****************************************************************************************
  // Read Input File
  //*****************************************************************************************
  TFile *infile=0;
  TTree *eventTree=0;

  // Data structures to store info from TTrees
  higgsana::TEventInfo *info  = new higgsana::TEventInfo();
  higgsana::TGenInfo *gen     = new higgsana::TGenInfo();
  TClonesArray *electronArr = new TClonesArray("higgsana::TElectron");
  TClonesArray *pfcandidateArr = new TClonesArray("higgsana::TPFCandidate");
  TClonesArray *muonArr = new TClonesArray("higgsana::TMuon");

  // Read input file and get the TTrees
  cout << "Processing " << inputfile << "..." << endl;
  infile = TFile::Open(inputfile.c_str(),"read");
  assert(infile);

  eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
  eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("Muon", &muonArr);         TBranch *muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr);  TBranch *pfcandidateBr = eventTree->GetBranch("PFCandidate");
  cout << "NEvents = " << eventTree->GetEntries() << endl;

  TBranch *genBr = 0;
  if(matchGen) {
    eventTree->SetBranchAddress("Gen", &gen);
    genBr = eventTree->GetBranch("Gen");
  }

  // Setting up environments for evaluating the regression
  ElectronEnergyRegressionEvaluate *evaluator_NoTrkVar = new ElectronEnergyRegressionEvaluate();
  ElectronEnergyRegressionEvaluate *evaluator_WithTrkVar = new ElectronEnergyRegressionEvaluate();
  ElectronEnergyRegressionEvaluate *evaluator_NoTrkVarTwoPtBins = new ElectronEnergyRegressionEvaluate();
  ElectronEnergyRegressionEvaluate *evaluator_WithTrkVarTwoPtBins = new ElectronEnergyRegressionEvaluate();

  evaluator_NoTrkVar->initialize ("/afs/cern.ch/work/a/atakeda/public/release/CMSSW_5_2_4_patch1/src/weights_Summer12HZZ4l/weightFile_V00.root",
		  ElectronEnergyRegressionEvaluate::kNoTrkVar);
  evaluator_WithTrkVar->initialize   ("/afs/cern.ch/work/a/atakeda/public/release/CMSSW_5_2_4_patch1/src/weights_Summer12HZZ4l/weightFile_V01.root",
		  ElectronEnergyRegressionEvaluate::kWithTrkVar);
  evaluator_NoTrkVarTwoPtBins->initialize   ("/afs/cern.ch/work/a/atakeda/public/release/CMSSW_5_2_4_patch1/src/weights_Summer12HZZ4l/weightFile_V10.root",
		  ElectronEnergyRegressionEvaluate::kNoTrkVarTwoPtBins);
  evaluator_WithTrkVarTwoPtBins->initialize     ("/afs/cern.ch/work/a/atakeda/public/release/CMSSW_5_2_4_patch1/src/weights_Summer12HZZ4l/weightFile_V11.root",
		  ElectronEnergyRegressionEvaluate::kWithTrkVarTwoPtBins);

  // Checking initialization
  assert(evaluator_NoTrkVar->isInitialized());
  assert(evaluator_WithTrkVar->isInitialized());
  assert(evaluator_NoTrkVarTwoPtBins->isInitialized());
  assert(evaluator_WithTrkVarTwoPtBins->isInitialized());


  // loop over events
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    if (ientry % 100000 == 0) cout << "Processed Event " << ientry << endl;
    infoBr->GetEntry(ientry);

    // check for certified runs
    mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
    if(applyGoodLumi && !rlrm.HasRunLumi(rl)) continue;  

    //***********************************************************
    // Definition of Pileup Energy density
    //***********************************************************
    Double_t rhoEleIso = 0;
    UInt_t EleEAEra = 0;
    if (DataEra == kDataEra_2011_MC) {

      if (!(isnan(info->RhoKt6PFJetsForIso25) || 
			      isinf(info->RhoKt6PFJetsForIso25))) {
	rhoEleIso = info->RhoKt6PFJetsForIso25;
      }
      EleEAEra = kDataEra_2011_Data;

    } else if (DataEra == kDataEra_2012_MC) {

      if (!(isnan(info->RhoKt6PFJets) || 
			      isinf(info->RhoKt6PFJets))) {
	rhoEleIso = info->RhoKt6PFJets;
      }
      EleEAEra = kDataEra_2012_Data;
    }

    //***********************************************************
    // Pileup Weight
    //***********************************************************
    double npuWeight = 1;
    if (fhDPU) {
      double mynpu = TMath::Min((double)info->nPUEvents,34.999);
      Int_t npuxbin = fhDPU->GetXaxis()->FindBin(mynpu);
      npuWeight = fhDPU->GetBinContent(npuxbin);
    }


    // trigger requirement               
    Bool_t passTrigger = kFALSE;
    UInt_t triggerBits = 0;
    if(PDType == 0) {
      if ((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) == kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) {
	passTrigger = kTRUE;
	triggerBits |= citana::ZeeEventTree::kEle17SC8Trigger;
      }
      if ((info->triggerBits & kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50) == kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50) {
	passTrigger = kTRUE;
	triggerBits |= citana::ZeeEventTree::kEle20SC4Trigger;
      }
      if ((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) == kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) {
	passTrigger = kTRUE;
	triggerBits |= citana::ZeeEventTree::kEle32SC17Trigger;
      }
      if ((info->triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL) == kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL) {
	passTrigger = kTRUE;
	triggerBits |= citana::ZeeEventTree::kEle17Ele8Loose;
      }
      if ((info->triggerBits & kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL) == kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL) {
	passTrigger = kTRUE;
	triggerBits |= citana::ZeeEventTree::kEle17Ele8Tight;
      }
    } else if(PDType == 1) {
      //if it's from single ele PD, then remove overlap by vetoing double ele triggers
      if ((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) == kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30)  continue;
      if ((info->triggerBits & kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50) == kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50)  continue;
      if ((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) == kHLT_Ele32_CaloIdL_CaloIsoVL_SC17)  continue;
      if ((info->triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL) == kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL)  continue;
      if ((info->triggerBits & kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL) == kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL)  continue;

      if ((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) == kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) {
	passTrigger = kTRUE;
	triggerBits |= citana::ZeeEventTree::kSingleEleTight;
      }
    }
    if(PDType != -1 && !passTrigger) continue;     

    // good vertex requirement
    if(!(info->hasGoodPV)) continue;

    if(matchGen) genBr->GetEntry(ientry);

    electronArr->Clear();
    muonArr->Clear(); 
    pfcandidateArr->Clear(); 
    electronBr->GetEntry(ientry);
    muonBr->GetEntry(ientry);
    pfcandidateBr->GetEntry(ientry);

    //********************************************************
    //Loop over both electrons
    //********************************************************
    Int_t Electron1Index = -1;
    Int_t Electron2Index = -1;
    Double_t Electron1Pt = -1;
    Double_t Electron2Pt = -1;
    Double_t Mass = 0;

    for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
      const higgsana::TElectron *ele1 = (higgsana::TElectron*)((*electronArr)[i]);

      if(matchGen) {
	Bool_t match1 = (higgsana::deltaR(ele1->eta, ele1->phi, gen->eta_1, gen->phi_1) < 0.5);
	Bool_t match2 = (higgsana::deltaR(ele1->eta, ele1->phi, gen->eta_2, gen->phi_2) < 0.5);
	if(!match1 && !match2)
	  continue;
      }
      if(fabs(ele1->scEta) > 2.5) continue;

      //make four vector for electron1
      const Double_t m = 0.000511;
      TLorentzVector vele1;
      vele1.SetPtEtaPhiM(ele1->pt, ele1->eta, ele1->phi, m);

      for(Int_t j=i+1; j<electronArr->GetEntriesFast(); j++) {

	const higgsana::TElectron *ele2 = (higgsana::TElectron*)((*electronArr)[j]);

	if(matchGen) {
	  Bool_t match1 = (higgsana::deltaR(ele2->eta, ele2->phi, gen->eta_1, gen->phi_1) < 0.5);
	  Bool_t match2 = (higgsana::deltaR(ele2->eta, ele2->phi, gen->eta_2, gen->phi_2) < 0.5);
	  if(!match1 && !match2)
	    continue;
	}

	//eta cuts 
	if(fabs(ele2->scEta) > 2.5) continue;

	//charge requirement
	if(ele2->q == ele1->q) continue;

	//pt cuts
	if (!(ele1->pt > 10 || ele2->pt > 10)) continue;
	if (!(ele1->pt > 10 && ele2->pt > 10)) continue;


	TLorentzVector vele2;
	vele2.SetPtEtaPhiM(ele2->pt, ele2->eta, ele2->phi, m);

	TLorentzVector vdielectron = vele1 + vele2;
	if((vdielectron.M()<massLo) || (vdielectron.M()>massHi)) continue;	  	  

	//choose pair with highest pt 
	double eleptmax = ele1->pt; 
	double eleptmin = ele2->pt; 
	Int_t eleptmaxIndex = i;
	Int_t eleptminIndex = j;
	if (ele2->pt > ele1->pt) {
	  eleptmax = ele2->pt;
	  eleptmin = ele1->pt;
	  eleptmaxIndex = j;
	  eleptminIndex = i;
	}

	if (eleptmax > Electron1Pt || eleptmin > Electron2Pt) {            
	  Electron1Index = eleptmaxIndex;
	  Electron2Index = eleptminIndex;
	  Electron1Pt = eleptmax;
	  Electron2Pt = eleptmin;
	  Mass = vdielectron.M();
	}
      }
    }

    if (!(Electron1Index >= 0 && Electron2Index >= 0)) continue;

    const higgsana::TElectron *Ele1 = (higgsana::TElectron*)((*electronArr)[Electron1Index]);
    const higgsana::TElectron *Ele2 = (higgsana::TElectron*)((*electronArr)[Electron2Index]);

    //find matching gen leptons
    Double_t GenPt1 = -1;
    Double_t GenPt2 = -1;
    if(matchGen) {
      if (higgsana::deltaR(Ele1->eta, Ele1->phi, gen->eta_1, gen->phi_1) < higgsana::deltaR(Ele1->eta, Ele1->phi, gen->eta_2, gen->phi_2)) {
	if (higgsana::deltaR(Ele1->eta, Ele1->phi, gen->eta_1, gen->phi_1) < 0.5) {
	  GenPt1 = gen->pt_1;
	  GenPt2 = gen->pt_2;
	}
      } else {
	if (higgsana::deltaR(Ele1->eta, Ele1->phi, gen->eta_2, gen->phi_2) < 0.5) {
	  GenPt1 = gen->pt_2;
	  GenPt2 = gen->pt_1;
	}
      }
    }


    //empty vector of photons 
    vector<const higgsana::TPFCandidate*> photonsToVeto;

    //******************
    //Fill Zee Event
    //******************
    zeeEventTree->fWeight = npuWeight;
    zeeEventTree->fRunNumber = info->runNum;
    zeeEventTree->fLumiSectionNumber = info->lumiSec;
    zeeEventTree->fEventNumber = info->evtNum;
    zeeEventTree->fRho = rhoEleIso;
    zeeEventTree->fNVertices = info->nPV0;
    zeeEventTree->fEventTriggerBits = triggerBits;
    zeeEventTree->fMass = Mass;

    zeeEventTree->fEle1Pt = Ele1->pt; 
    zeeEventTree->fEle1Eta = Ele1->eta; 
    zeeEventTree->fEle1Phi = Ele1->phi; 
    zeeEventTree->fEle1SCEt = Ele1->scEt; 
    zeeEventTree->fEle1SCEta = Ele1->scEta; 
    zeeEventTree->fEle1SCPhi = Ele1->scPhi; 
    zeeEventTree->fEle1GenPt = GenPt1; 
    zeeEventTree->fEle1EnergyCorrAndSmeared = Ele1->EcalEnergy; 
    zeeEventTree->fEle1Energy = Ele1->EcalEnergy; 
    zeeEventTree->fEle1EnergyRegression = 0; 
    zeeEventTree->fEle1Charge = Ele1->q; 
    zeeEventTree->fEle1HZZICHEP2012IDMVA = EvaluateEleHZZ4lICHEP2012IDMVA(Ele1, eleIDMVA); 
    zeeEventTree->fEle1PFIso04 = ComputeElePFIso04(Ele1,pfcandidateArr,rhoEleIso,EleEAEra);
    zeeEventTree->fEle1R9 = Ele1->R9;
    zeeEventTree->fEle1PassLooseSimpleCuts = PassEleSimpleCutsLoose( Ele1,pfcandidateArr,rhoEleIso,EleEAEra); 
    zeeEventTree->fEle1PassMediumSimpleCuts = PassEleSimpleCutsMedium( Ele1,pfcandidateArr,rhoEleIso,EleEAEra); 
    zeeEventTree->fEle1PassTightSimpleCuts = PassEleSimpleCutsTight( Ele1,pfcandidateArr,rhoEleIso,EleEAEra);
    zeeEventTree->fEle1PassHZZICHEP2012 = (PassEleHZZ4lICHEP2012ID(Ele1,eleIDMVA) && PassEleHZZ4lICHEP2012Iso(Ele1,pfcandidateArr,rhoEleIso,EleEAEra,photonsToVeto)); 

    zeeEventTree->fEle2Pt = Ele2->pt; 
    zeeEventTree->fEle2Eta = Ele2->eta; 
    zeeEventTree->fEle2Phi = Ele2->phi; 
    zeeEventTree->fEle2SCEt = Ele2->scEt; 
    zeeEventTree->fEle2SCEta = Ele2->scEta; 
    zeeEventTree->fEle2SCPhi = Ele2->scPhi; 
    zeeEventTree->fEle2GenPt = GenPt2; 
    zeeEventTree->fEle2EnergyCorrAndSmeared = Ele2->EcalEnergy; 
    zeeEventTree->fEle2Energy = Ele2->EcalEnergy; 
    zeeEventTree->fEle2EnergyRegression = 0; 
    zeeEventTree->fEle2Charge = Ele2->q; 
    zeeEventTree->fEle2HZZICHEP2012IDMVA = EvaluateEleHZZ4lICHEP2012IDMVA(Ele2, eleIDMVA); 
    zeeEventTree->fEle2PFIso04 = ComputeElePFIso04(Ele2,pfcandidateArr,rhoEleIso,EleEAEra);
    zeeEventTree->fEle2R9 = Ele2->R9;
    zeeEventTree->fEle2PassLooseSimpleCuts = PassEleSimpleCutsLoose( Ele2,pfcandidateArr,rhoEleIso,EleEAEra); 
    zeeEventTree->fEle2PassMediumSimpleCuts = PassEleSimpleCutsMedium( Ele2,pfcandidateArr,rhoEleIso,EleEAEra); 
    zeeEventTree->fEle2PassTightSimpleCuts = PassEleSimpleCutsTight( Ele2,pfcandidateArr,rhoEleIso,EleEAEra);
    zeeEventTree->fEle2PassHZZICHEP2012 = (PassEleHZZ4lICHEP2012ID(Ele2,eleIDMVA) && PassEleHZZ4lICHEP2012Iso(Ele2,pfcandidateArr,rhoEleIso,EleEAEra,photonsToVeto)); 

    /* =================== Using regression evaluators to compute regression energies ===================== */

    double Ele1sep;
    if (Ele1->sigiEtaiEta*Ele1->sigiPhiiPhi > 0) {
      Ele1sep = Ele1->CovIEtaIPhi/(Ele1->sigiEtaiEta*Ele1->sigiPhiiPhi);
    } else if (Ele1->CovIEtaIPhi>0) {
      Ele1sep = 1.0; 
    } else {
      Ele1sep = -1.0; 
    }
    double Ele2sep;
    if (Ele2->sigiEtaiEta*Ele2->sigiPhiiPhi > 0) {
      Ele2sep = Ele2->CovIEtaIPhi/(Ele2->sigiEtaiEta*Ele2->sigiPhiiPhi);
    } else if (Ele2->CovIEtaIPhi>0) {
      Ele2sep = 1.0; 
    } else {
      Ele2sep = -1.0; 
    }

    zeeEventTree->fEle1EnergyRegressionNoTrkVar = evaluator_NoTrkVar->regressionValueNoTrkVar(
		    Ele1->SCRawEnergy,
		    Ele1->scEta,
		    Ele1->scPhi,
		    Ele1->R9,
		    Ele1->E5x5,
		    Ele1->SCEtaWidth,
		    Ele1->SCPhiWidth,
		    Ele1->nBrem + 1,             //nclusters = nbrem + 1
		    Ele1->HoverE,
		    rhoEleIso,        // as calculated above,
		    info->nPV0,     // see above
		    Ele1->EtaSeed,
		    Ele1->PhiSeed,
		    Ele1->ESeed,
		    Ele1->E3x3Seed,
		    Ele1->E5x5Seed,
		    Ele1->sigiEtaiEta,
		    Ele1->sigiPhiiPhi,
		    Ele1sep,
		    Ele1->EMaxSeed,
		    Ele1->E2ndSeed,
		    Ele1->ETopSeed,
		    Ele1->EBottomSeed,
		    Ele1->ELeftSeed,
		    Ele1->ERightSeed,
		    Ele1->E2x5MaxSeed,
		    Ele1->E2x5TopSeed,
		    Ele1->E2x5BottomSeed,
		    Ele1->E2x5LeftSeed,
		    Ele1->E2x5RightSeed,
		    Ele1->pt,
		    Ele1->IEtaSeed,
		    Ele1->IPhiSeed,
		    Ele1->EtaCrySeed,
		    Ele1->PhiCrySeed,
		    Ele1->PreShowerOverRaw );

    zeeEventTree->fEle1EnergyRegressionNoTrkVarTwoPtBins = evaluator_NoTrkVarTwoPtBins->regressionValueNoTrkVar(
		    Ele1->SCRawEnergy,
		    Ele1->scEta,
		    Ele1->scPhi,
		    Ele1->R9,
		    Ele1->E5x5,
		    Ele1->SCEtaWidth,
		    Ele1->SCPhiWidth,
		    Ele1->nBrem + 1,             //nclusters = nbrem + 1
		    Ele1->HoverE,
		    rhoEleIso,        // as calculated above,
		    info->nPV0,     // see above
		    Ele1->EtaSeed,
		    Ele1->PhiSeed,
		    Ele1->ESeed,
		    Ele1->E3x3Seed,
		    Ele1->E5x5Seed,
		    Ele1->sigiEtaiEta,
		    Ele1->sigiPhiiPhi,
		    Ele1sep,
		    Ele1->EMaxSeed,
		    Ele1->E2ndSeed,
		    Ele1->ETopSeed,
		    Ele1->EBottomSeed,
		    Ele1->ELeftSeed,
		    Ele1->ERightSeed,
		    Ele1->E2x5MaxSeed,
		    Ele1->E2x5TopSeed,
		    Ele1->E2x5BottomSeed,
		    Ele1->E2x5LeftSeed,
		    Ele1->E2x5RightSeed,
		    Ele1->pt,
		    Ele1->IEtaSeed,
		    Ele1->IPhiSeed,
		    Ele1->EtaCrySeed,
		    Ele1->PhiCrySeed,
		    Ele1->PreShowerOverRaw );

    zeeEventTree->fEle1EnergyRegressionWithTrkVar = evaluator_WithTrkVar->regressionValueWithTrkVar(
		    Ele1->SCRawEnergy,
		    Ele1->scEta,
		    Ele1->scPhi,
		    Ele1->R9,
		    Ele1->E5x5,
		    Ele1->SCEtaWidth,
		    Ele1->SCPhiWidth,
		    Ele1->nBrem + 1,             //nclusters = nbrem + 1
		    Ele1->HoverE,
		    rhoEleIso,        // as calculated above,
		    info->nPV0,     // see above
		    Ele1->EtaSeed,
		    Ele1->PhiSeed,
		    Ele1->ESeed,
		    Ele1->E3x3Seed,
		    Ele1->E5x5Seed,
		    Ele1->sigiEtaiEta,
		    Ele1->sigiPhiiPhi,
		    Ele1sep,
		    Ele1->EMaxSeed,
		    Ele1->E2ndSeed,
		    Ele1->ETopSeed,
		    Ele1->EBottomSeed,
		    Ele1->ELeftSeed,
		    Ele1->ERightSeed,
		    Ele1->E2x5MaxSeed,
		    Ele1->E2x5TopSeed,
		    Ele1->E2x5BottomSeed,
		    Ele1->E2x5LeftSeed,
		    Ele1->E2x5RightSeed,
		    Ele1->pt,
		    Ele1->pIn,
		    Ele1->fBrem,
		    Ele1->q,
		    Ele1->EOverP,
		    Ele1->IEtaSeed,
		    Ele1->IPhiSeed,
		    Ele1->EtaCrySeed,
		    Ele1->PhiCrySeed,
		    Ele1->PreShowerOverRaw );

    zeeEventTree->fEle1EnergyRegressionWithTrkVarTwoPtBins = evaluator_WithTrkVarTwoPtBins->regressionValueWithTrkVar(
		    Ele1->SCRawEnergy,
		    Ele1->scEta,
		    Ele1->scPhi,
		    Ele1->R9,
		    Ele1->E5x5,
		    Ele1->SCEtaWidth,
		    Ele1->SCPhiWidth,
		    Ele1->nBrem + 1,             //nclusters = nbrem + 1
		    Ele1->HoverE,
		    rhoEleIso,        // as calculated above,
		    info->nPV0,     // see above
		    Ele1->EtaSeed,
		    Ele1->PhiSeed,
		    Ele1->ESeed,
		    Ele1->E3x3Seed,
		    Ele1->E5x5Seed,
		    Ele1->sigiEtaiEta,
		    Ele1->sigiPhiiPhi,
		    Ele1sep,
		    Ele1->EMaxSeed,
		    Ele1->E2ndSeed,
		    Ele1->ETopSeed,
		    Ele1->EBottomSeed,
		    Ele1->ELeftSeed,
		    Ele1->ERightSeed,
		    Ele1->E2x5MaxSeed,
		    Ele1->E2x5TopSeed,
		    Ele1->E2x5BottomSeed,
		    Ele1->E2x5LeftSeed,
		    Ele1->E2x5RightSeed,
		    Ele1->pt,
		    Ele1->pIn,
		    Ele1->fBrem,
		    Ele1->q,
		    Ele1->EOverP,
		    Ele1->IEtaSeed,
		    Ele1->IPhiSeed,
		    Ele1->EtaCrySeed,
		    Ele1->PhiCrySeed,
		    Ele1->PreShowerOverRaw );

    // Now for electron 2
    zeeEventTree->fEle2EnergyRegressionNoTrkVar = evaluator_NoTrkVar->regressionValueNoTrkVar(
		    Ele2->SCRawEnergy,
		    Ele2->scEta,
		    Ele2->scPhi,
		    Ele2->R9,
		    Ele2->E5x5,
		    Ele2->SCEtaWidth,
		    Ele2->SCPhiWidth,
		    Ele2->nBrem + 1,             //nclusters = nbrem + 1
		    Ele2->HoverE,
		    rhoEleIso,        // as calculated above,
		    info->nPV0,     // see above
		    Ele2->EtaSeed,
		    Ele2->PhiSeed,
		    Ele2->ESeed,
		    Ele2->E3x3Seed,
		    Ele2->E5x5Seed,
		    Ele2->sigiEtaiEta,
		    Ele2->sigiPhiiPhi,
                    Ele2sep,
		    Ele2->EMaxSeed,
		    Ele2->E2ndSeed,
		    Ele2->ETopSeed,
		    Ele2->EBottomSeed,
		    Ele2->ELeftSeed,
		    Ele2->ERightSeed,
		    Ele2->E2x5MaxSeed,
		    Ele2->E2x5TopSeed,
		    Ele2->E2x5BottomSeed,
		    Ele2->E2x5LeftSeed,
		    Ele2->E2x5RightSeed,
		    Ele2->pt,
		    Ele2->IEtaSeed,
		    Ele2->IPhiSeed,
		    Ele2->EtaCrySeed,
		    Ele2->PhiCrySeed,
		    Ele2->PreShowerOverRaw );

    zeeEventTree->fEle2EnergyRegressionNoTrkVarTwoPtBins = evaluator_NoTrkVarTwoPtBins->regressionValueNoTrkVar(
		    Ele2->SCRawEnergy,
		    Ele2->scEta,
		    Ele2->scPhi,
		    Ele2->R9,
		    Ele2->E5x5,
		    Ele2->SCEtaWidth,
		    Ele2->SCPhiWidth,
		    Ele2->nBrem + 1,             //nclusters = nbrem + 1
		    Ele2->HoverE,
		    rhoEleIso,        // as calculated above,
		    info->nPV0,     // see above
		    Ele2->EtaSeed,
		    Ele2->PhiSeed,
		    Ele2->ESeed,
		    Ele2->E3x3Seed,
		    Ele2->E5x5Seed,
		    Ele2->sigiEtaiEta,
		    Ele2->sigiPhiiPhi,
                    Ele2sep,
		    Ele2->EMaxSeed,
		    Ele2->E2ndSeed,
		    Ele2->ETopSeed,
		    Ele2->EBottomSeed,
		    Ele2->ELeftSeed,
		    Ele2->ERightSeed,
		    Ele2->E2x5MaxSeed,
		    Ele2->E2x5TopSeed,
		    Ele2->E2x5BottomSeed,
		    Ele2->E2x5LeftSeed,
		    Ele2->E2x5RightSeed,
		    Ele2->pt,
		    Ele2->IEtaSeed,
		    Ele2->IPhiSeed,
		    Ele2->EtaCrySeed,
		    Ele2->PhiCrySeed,
		    Ele2->PreShowerOverRaw );

    zeeEventTree->fEle2EnergyRegressionWithTrkVar = evaluator_WithTrkVar->regressionValueWithTrkVar(
		    Ele2->SCRawEnergy,
		    Ele2->scEta,
		    Ele2->scPhi,
		    Ele2->R9,
		    Ele2->E5x5,
		    Ele2->SCEtaWidth,
		    Ele2->SCPhiWidth,
		    Ele2->nBrem + 1,             //nclusters = nbrem + 1
		    Ele2->HoverE,
		    rhoEleIso,        // as calculated above,
		    info->nPV0,     // see above
		    Ele2->EtaSeed,
		    Ele2->PhiSeed,
		    Ele2->ESeed,
		    Ele2->E3x3Seed,
		    Ele2->E5x5Seed,
		    Ele2->sigiEtaiEta,
		    Ele2->sigiPhiiPhi,
                    Ele2sep,
		    Ele2->EMaxSeed,
		    Ele2->E2ndSeed,
		    Ele2->ETopSeed,
		    Ele2->EBottomSeed,
		    Ele2->ELeftSeed,
		    Ele2->ERightSeed,
		    Ele2->E2x5MaxSeed,
		    Ele2->E2x5TopSeed,
		    Ele2->E2x5BottomSeed,
		    Ele2->E2x5LeftSeed,
		    Ele2->E2x5RightSeed,
		    Ele2->pt,
		    Ele2->pIn,
		    Ele2->fBrem,
		    Ele2->q,
		    Ele2->EOverP,
		    Ele2->IEtaSeed,
		    Ele2->IPhiSeed,
		    Ele2->EtaCrySeed,
		    Ele2->PhiCrySeed,
		    Ele2->PreShowerOverRaw );

    zeeEventTree->fEle2EnergyRegressionWithTrkVarTwoPtBins = evaluator_WithTrkVarTwoPtBins->regressionValueWithTrkVar(
		    Ele2->SCRawEnergy,
		    Ele2->scEta,
		    Ele2->scPhi,
		    Ele2->R9,
		    Ele2->E5x5,
		    Ele2->SCEtaWidth,
		    Ele2->SCPhiWidth,
		    Ele2->nBrem + 1,             //nclusters = nbrem + 1
		    Ele2->HoverE,
		    rhoEleIso,        // as calculated above,
		    info->nPV0,     // see above
		    Ele2->EtaSeed,
		    Ele2->PhiSeed,
		    Ele2->ESeed,
		    Ele2->E3x3Seed,
		    Ele2->E5x5Seed,
		    Ele2->sigiEtaiEta,
		    Ele2->sigiPhiiPhi,
                    Ele2sep,
		    Ele2->EMaxSeed,
		    Ele2->E2ndSeed,
		    Ele2->ETopSeed,
		    Ele2->EBottomSeed,
		    Ele2->ELeftSeed,
		    Ele2->ERightSeed,
		    Ele2->E2x5MaxSeed,
		    Ele2->E2x5TopSeed,
		    Ele2->E2x5BottomSeed,
		    Ele2->E2x5LeftSeed,
		    Ele2->E2x5RightSeed,
		    Ele2->pt,
		    Ele2->pIn,
		    Ele2->fBrem,
		    Ele2->q,
		    Ele2->EOverP,
		    Ele2->IEtaSeed,
		    Ele2->IPhiSeed,
		    Ele2->EtaCrySeed,
		    Ele2->PhiCrySeed,
		    Ele2->PreShowerOverRaw );

    zeeEventTree->tree_->Fill();
    nEvents++;

  }


  // Cleaning up the evaluators
  delete evaluator_NoTrkVar;
  delete evaluator_WithTrkVar;
  delete evaluator_NoTrkVarTwoPtBins;
  delete evaluator_WithTrkVarTwoPtBins;
  delete infile;
  infile=0, eventTree=0;    

  delete info;
  delete gen;
  delete electronArr;


  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  cout << " Number of events selected: " << nEvents << endl;

  outFile->Write();
  outFile->Close();
  delete outFile;

  cout << endl;
  cout << "  <> Output saved in " << outputfile << endl;    
  cout << endl;  

}
