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

// #include "Common/MyTools.hh"        // miscellaneous helper functions

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
#include "HiggsAna/HZGamma/Utils/LeptonSelection.hh"

// Evaluating function for regression
//#include "EGamma/EGammaAnalysisTools/interface/ElectronEnergyRegressionEvaluate.h"

// structure for output ntuple
#include "CITCommon/CommonData/interface/ZmumuGammaEventTree.h"
#endif


//=== MAIN MACRO ================================================================================================= 

void MakeZmumuGammaEventNtuples(const string inputfile,    // input file
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


  double eventWeight = 1;
//53X MC
//   if (PDType == -1) eventWeight = 10000 * 1177.6*3 / 30391580.0; //mll50
//   if (PDType == -2) eventWeight = 10000 * 3969.61*3.0*0.0498327*1.45 / 1700850.0 ; //mll1050
//52X MC
  if (PDType == -1) eventWeight = 10000 * 1177.6*3 / 26821030.0; //mll50
  if (PDType == -2) eventWeight = 10000 * 3969.61*3.0*0.0498327*1.45 / 6046518.0 ; //mll1050

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
  rlrm.AddJSONFile("/afs/cern.ch/work/s/sixie/public/HZGamma/auxiliar/Cert_190456-196531_8TeV_29Jun2012ReReco_Collisions12_JSON.txt"); 
  rlrm.AddJSONFile("/afs/cern.ch/work/s/sixie/public/HZGamma/auxiliar/Cert_190782-190949_8TeV_29Jun2012ReReco-recover_Collisions12_JSON.txt");
  rlrm.AddJSONFile("/afs/cern.ch/work/s/sixie/public/HZGamma/auxiliar/Cert_Full2012_53X_JSON.txt");


  //********************************************************
  // mass region
  //********************************************************
//   Double_t massLo = 40;
//   Double_t massHi = 200;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  Double_t nEvents = 0;

  //*****************************************************************************************
  // Set up output ntuple
  //*****************************************************************************************
  TFile *outFile = new TFile(outputfile.c_str(),"RECREATE"); 

  citana::ZmumuGammaEventTree *zmumuGammaEventTree = new citana::ZmumuGammaEventTree;
  zmumuGammaEventTree->CreateTree();
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
  TClonesArray *genparticleArr = new TClonesArray("higgsana::TGenParticle");
  TClonesArray *muonArr = new TClonesArray("higgsana::TMuon");
  TClonesArray *photonArr = new TClonesArray("higgsana::TPhoton");

  // Read input file and get the TTrees
  cout << "Processing " << inputfile << "..." << endl;
  infile = TFile::Open(inputfile.c_str(),"read");
  assert(infile);

  eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
  eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("Muon", &muonArr);         TBranch *muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("Photon", &photonArr);     TBranch *photonBr = eventTree->GetBranch("Photon");
  eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr);  TBranch *pfcandidateBr = eventTree->GetBranch("PFCandidate");
  cout << "NEvents = " << eventTree->GetEntries() << endl;

  TBranch *genparticleBr = 0;
  if(matchGen) {
    eventTree->SetBranchAddress("GenParticle", &genparticleArr);         
    genparticleBr = eventTree->GetBranch("GenParticle");
  }


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
    Double_t rhoMuonIso = 0;
    Double_t rhoEleIso = 0;
    Double_t rhoPhoIso = 0;
    UInt_t MuonEAEra = 0;
    UInt_t EleEAEra = 0;
    UInt_t PhoEAEra = 0;

    if (DataEra == kDataEra_2011_MC) {
     
      if (!(isnan(info->RhoKt6PFJetsForIso25) || 
            isinf(info->RhoKt6PFJetsForIso25))) {
        rhoMuonIso = info->RhoKt6PFJetsForIso25;
        rhoEleIso = info->RhoKt6PFJetsForIso25;
        rhoPhoIso = info->RhoKt6PFJetsForIso25;
      }
      MuonEAEra = kDataEra_2011_Data;
      EleEAEra = kDataEra_2011_Data;
      PhoEAEra = kDataEra_2011_Data;
    } else if (DataEra == kDataEra_2012_MC) {

      if (!(isnan(info->RhoKt6PFJetsCentralNeutral) || 
            isinf(info->RhoKt6PFJetsCentralNeutral))) {
        rhoMuonIso = info->RhoKt6PFJetsCentralNeutral;
      }

      if (!(isnan(info->RhoKt6PFJets) || 
            isinf(info->RhoKt6PFJets))) {
        rhoEleIso = info->RhoKt6PFJets;
        rhoPhoIso = info->RhoKt6PFJets;
      }

      MuonEAEra = kDataEra_2012_Data;
      EleEAEra = kDataEra_2012_Data;
      PhoEAEra = kDataEra_2012_Data;
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
      if ( (info->triggerBits & kHLT_DoubleMu7) == kHLT_DoubleMu7 ) {
	passTrigger = kTRUE;
      }
      if ( (info->triggerBits & kHLT_Mu13_Mu8) == kHLT_Mu13_Mu8 ) {
	passTrigger = kTRUE;
      }
      if ( (info->triggerBits & kHLT_Mu17_Mu8) == kHLT_Mu17_Mu8 ) {
	passTrigger = kTRUE;
      }
    } else if(PDType == 1) {

      //if it's from single ele PD, then remove overlap by vetoing double ele triggers
      if ((info->triggerBits & kHLT_DoubleMu7) == kHLT_DoubleMu7)  continue;
      if ((info->triggerBits & kHLT_Mu13_Mu8) == kHLT_Mu13_Mu8)  continue;
      if ((info->triggerBits & kHLT_Mu17_Mu8) == kHLT_Mu17_Mu8)  continue;

      if ( (info->triggerBits & kHLT_IsoMu24) == kHLT_IsoMu24 ) {
	passTrigger = kTRUE;
      }

    }

//     if( PDType != -1 && PDType != -2 &&  !passTrigger) continue;     
     
    // good vertex requirement
    if(!(info->hasGoodPV)) continue;

    if(matchGen) genparticleBr->GetEntry(ientry);
      

    electronArr->Clear();
    muonArr->Clear(); 
    photonArr->Clear(); 
    pfcandidateArr->Clear(); 
    electronBr->GetEntry(ientry);
    muonBr->GetEntry(ientry);
    photonBr->GetEntry(ientry);
    pfcandidateBr->GetEntry(ientry);

    //********************************************************
    //Loop over muons
    //********************************************************
//     Int_t Muon1Index = -1;
//     Int_t Muon2Index = -1;
//     Double_t Muon1Pt = -1;
//     Double_t Muon2Pt = -1;
//     Double_t Mass = 0;

    for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
      const higgsana::TMuon *mu1 = (higgsana::TMuon*)((*muonArr)[i]);

      if(fabs(mu1->eta) > 2.4) continue;
      if (!PassMuonHZGammaID(mu1, (DataEra == kDataEra_2011_MC) )) continue;
      if (! (ComputeMuonPFIsoRings(mu1, pfcandidateArr, rhoMuonIso,MuonEAEra, kPFChargedIso, 0.0, 0.4) / mu1->pt < 0.2)) continue;

//       cout << "mu: " << mu1->pt << " " << mu1->eta << " " << mu1->phi << " : " << Bool_t((mu1->typeBits & kGlobal) == kGlobal) << " " << mu1->PassPFId << " " << mu1->muNchi2 << " " << mu1->nValidHits << " " << mu1->nMatch << " " << mu1->d0 << " " << mu1->dz << " " << mu1->nPixHits << " " << mu1->trkLayers << " " << PassMuonHZGammaID(mu1, (DataEra == kDataEra_2011_MC) )  << endl;

      //make four vector for muon1
      TLorentzVector vmu1;
      vmu1.SetPtEtaPhiM(mu1->pt, mu1->eta, mu1->phi, MUONMASS);

      for(Int_t j=i+1; j<muonArr->GetEntriesFast(); j++) {

	const higgsana::TMuon *mu2 = (higgsana::TMuon*)((*muonArr)[j]);

	if(fabs(mu2->eta) > 2.4) continue;
        if (!PassMuonHZGammaID(mu2, (DataEra == kDataEra_2011_MC) )) continue;
        if (! (ComputeMuonPFIsoRings(mu2, pfcandidateArr, rhoMuonIso,MuonEAEra, kPFChargedIso, 0.0, 0.4) / mu2->pt < 0.2)) continue;

	//charge requirement
	if(mu2->q == mu1->q) continue;

	//pt cuts
	if (!(mu1->pt > 10 || mu2->pt > 10)) continue;
	if (!(mu1->pt > 10 && mu2->pt > 10)) continue;


	TLorentzVector vmu2;
	vmu2.SetPtEtaPhiM(mu2->pt, mu2->eta, mu2->phi, MUONMASS);

 	TLorentzVector vdimuon = vmu1 + vmu2;

        for(Int_t p=0; p<photonArr->GetEntriesFast(); ++p) {
          
          const higgsana::TPhoton *pho = (higgsana::TPhoton*)((*photonArr)[p]);
          
//           cout << "Photon: " << pho->energyRegressionHgg2012/cosh(pho->eta) << " " << pho->eta << " " << pho->phi << " : " 
//                << PassPhotonHZGammaID(pho, (DataEra == kDataEra_2011_MC)) << endl;

          if (!(pho->energyRegressionHgg2012/cosh(pho->eta) > 10)) continue;
          if (!(fabs(pho->scEta) < 2.5)) continue;
          if ( fabs(pho->scEta) > 1.4442 && fabs(pho->scEta) < 1.566) continue;
          //noisy photon removal
          if (pho->scEta > -1.78 && pho->scEta < -1.75 && pho->scPhi > 1.36 && pho->scPhi < 1.39) continue;
          
          if (!PassPhotonHZGammaID(pho, (DataEra == kDataEra_2011_MC))) {
            continue;
          }
          
          //ignore this for now. no iso applied.
//           if (!PassPhotonHZGammaIso( pho, pfcandidateArr, rhoPhoIso,PhoEAEra)) {
//             continue;
//           }

          TLorentzVector vpho;
          //Use SC Energy
          //vpho.SetPtEtaPhiM(pho->et, pho->eta, pho->phi, 0);
          //Use regression energy
          vpho.SetPtEtaPhiM(pho->energyRegressionHgg2012/cosh(pho->eta) , pho->eta, pho->phi, 0);

         //******************************************************************
          //find matching gen leptons
          //******************************************************************
    
          const higgsana::TGenParticle *genMu1 = 0;
          const higgsana::TGenParticle *genMu2 = 0;
          const higgsana::TGenParticle *genPhoton = 0;
          for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
            const higgsana::TGenParticle *gen = (higgsana::TGenParticle*)((*genparticleArr)[k]);
            if (gen->status == 1 && (abs(gen->pdgid) == 13 ) ) {

              //find gen mu1
              if (higgsana::deltaR(gen->eta, gen->phi, mu1->eta, mu1->phi) < 0.1) {
                if (!genMu1 || higgsana::deltaR(gen->eta, gen->phi, mu1->eta, mu1->phi) < higgsana::deltaR(genMu1->eta, genMu1->phi, mu1->eta, mu1->phi)) {
                  genMu1 = gen;
                }
              }

              //find gen mu2
              if (higgsana::deltaR(gen->eta, gen->phi, mu2->eta, mu2->phi) < 0.1) {
                if (!genMu2 || higgsana::deltaR(gen->eta, gen->phi, mu2->eta, mu2->phi) < higgsana::deltaR(genMu2->eta, genMu2->phi, mu2->eta, mu2->phi)) {
                  genMu2 = gen;
                }
              }
            } //end if muon

            //find photon
            if (gen->status == 1 && (abs(gen->pdgid) == 22 ) ) {
              if (higgsana::deltaR(gen->eta, gen->phi, pho->eta, pho->phi) < 0.1) {
                if (!genPhoton || higgsana::deltaR(gen->eta, gen->phi, pho->eta, pho->phi) < higgsana::deltaR(genPhoton->eta, genPhoton->phi, pho->eta, pho->phi)) {
                  genPhoton = gen;
                }
              }
            }
          }

          //******************
          //Fill ZmumuGamma Event
          //******************
//           cout << "Trigger: " 
//                << Bool_t ( (info->triggerBits & kHLT_DoubleMu7) == kHLT_DoubleMu7 ) << " " 
//                << Bool_t( (info->triggerBits & kHLT_Mu13_Mu8) == kHLT_Mu13_Mu8 ) << " "
//                << Bool_t( (info->triggerBits & kHLT_Mu17_Mu8) == kHLT_Mu17_Mu8 ) << " "
//                << Bool_t( (info->triggerBits &  kHLT_Mu17_TkMu8) ==  kHLT_Mu17_TkMu8 ) << " "
//                << " : " << passTrigger
//                << endl;


          zmumuGammaEventTree->fWeight = npuWeight * eventWeight;
          zmumuGammaEventTree->fRunNumber = info->runNum;
          zmumuGammaEventTree->fLumiSectionNumber = info->lumiSec;
          zmumuGammaEventTree->fEventNumber = info->evtNum;
          zmumuGammaEventTree->fNPU = info->nPUEvents;
          zmumuGammaEventTree->fRho = rhoEleIso; 
          zmumuGammaEventTree->fNVertices = info->nPV0; 
          zmumuGammaEventTree->fMass = (vmu1+vmu2+vpho).M();
          zmumuGammaEventTree->fDileptonMass = (vmu1+vmu2).M();
          zmumuGammaEventTree->fPhotonPt = vpho.Pt(); 
          zmumuGammaEventTree->fPhotonEta = vpho.Eta(); 
          zmumuGammaEventTree->fPhotonPhi = vpho.Phi(); 
          zmumuGammaEventTree->fPhotonR9 = pho->R9; 
          zmumuGammaEventTree->fPhotonIsEB = (fabs(pho->scEta) < 1.566); 
          zmumuGammaEventTree->fPhotonHoE = pho->HoESingleTower; 
          zmumuGammaEventTree->fSCEt = pho->scEt; 
          zmumuGammaEventTree->fSCE = pho->scEt*cosh(pho->scEta); 
          zmumuGammaEventTree->fSCRawE = pho->scEt*cosh(pho->scEta);  
          zmumuGammaEventTree->fSCEta = pho->scEta; 
          zmumuGammaEventTree->fSCEtaWidth = 0; //for cluster correction calculations only
          zmumuGammaEventTree->fSCPhiWidth = 0; //for cluster correction calculations only
          zmumuGammaEventTree->fPhoToTrackDeltaR = 0; // do we need this? //old CiC selection ele veto
          zmumuGammaEventTree->fPhoPassElectronVeto = pho->passConversionSafeEleVeto; 
          zmumuGammaEventTree->fPhoHasMatchedConversion = 0;  //do we need this? no.
          zmumuGammaEventTree->fPhoHasPixelMatch = pho->hasPixelSeed; 
          zmumuGammaEventTree->fPhoSigmaIEtaIEta = pho->sigiEtaiEta; 
          zmumuGammaEventTree->fPhoTrackIso = pho->trkIso03Hollow; 
          zmumuGammaEventTree->fPhoEcalIso = pho->emIso03; 
          zmumuGammaEventTree->fPhoHcalIso = pho->hadIso03; 
          zmumuGammaEventTree->fPhoPdgId = 0;  //what's this used for? use deltaR circle. then use closest pt particle within circle. 0.1. mostly used to identify if it's a real photon or not.
          zmumuGammaEventTree->fPhoCrackCorr = 0; //part of cluster corrections.
          zmumuGammaEventTree->fMu1Pt = mu1->pt; 
          zmumuGammaEventTree->fMu1Eta = mu1->eta; 
          zmumuGammaEventTree->fMu1Phi = mu1->phi; 
          zmumuGammaEventTree->fMu1Charge = mu1->q; 
          zmumuGammaEventTree->fMu1TrackChi2 = 0; //do we need this? useful to see kinks in muons, but not needed for now
          zmumuGammaEventTree->fMu1TrackNormalizedChi2 = mu1->tkNchi2; 
          zmumuGammaEventTree->fMu1DeltaR = higgsana::deltaR(pho->eta, pho->phi, mu1->eta, mu1->phi); 
          zmumuGammaEventTree->fMu1CalEnergyEm = mu1->EmEnergy; //not used
          zmumuGammaEventTree->fMu1CalEnergyEmMax = 0; //what's this? not used
          zmumuGammaEventTree->fMu1CalEnergyEmHad = 0; //what's this? not used
          zmumuGammaEventTree->fMu2Pt = mu2->pt; 
          zmumuGammaEventTree->fMu2Eta= mu2->eta; 
          zmumuGammaEventTree->fMu2Phi = mu2->phi; 
          zmumuGammaEventTree->fMu2Charge = mu2->q; 
          zmumuGammaEventTree->fMu2TrackChi2 = 0;  //do we need this?
          zmumuGammaEventTree->fMu2TrackNormalizedChi2 = mu2->tkNchi2; 
          zmumuGammaEventTree->fMu2DeltaR = higgsana::deltaR(pho->eta, pho->phi, mu2->eta, mu2->phi);  
          zmumuGammaEventTree->fMu2CalEnergyEm = mu2->EmEnergy; //not used
          zmumuGammaEventTree->fMu2CalEnergyEmMax = 0; //what's this?
          zmumuGammaEventTree->fMu2CalEnergyEmHad = 0; //what's this?

          zmumuGammaEventTree->fMinDeltaEta = fmin( fabs(mu1->eta - pho->eta), fabs(mu2->eta - pho->eta) );
          zmumuGammaEventTree->fMinDeltaR = fmin(higgsana::deltaR(mu1->eta, mu1->phi, pho->eta, pho->phi),higgsana::deltaR(mu2->eta, mu2->phi, pho->eta, pho->phi));
          zmumuGammaEventTree->fMinDeltaPhi = fmin(higgsana::deltaPhi(mu1->phi, pho->phi),higgsana::deltaPhi(mu2->phi, pho->phi));
          zmumuGammaEventTree->fkRatio = 0; //what's that? photon energy calculated from on-shell z assumption to reconstructed photon energy
          zmumuGammaEventTree->fPreshowerE = 0; //also for cluster corrections

          if (genMu1) {
            zmumuGammaEventTree->fGenMu1Pt = genMu1->pt; 
            zmumuGammaEventTree->fGenMu1Eta = genMu1->eta; 
            zmumuGammaEventTree->fGenMu1Phi = genMu1->phi; 
          }
          if (genMu2) {
            zmumuGammaEventTree->fGenMu2Pt = genMu2->pt; 
            zmumuGammaEventTree->fGenMu2Eta = genMu2->eta; 
            zmumuGammaEventTree->fGenMu2Phi = genMu2->phi;
          }

          if (genPhoton) {
            zmumuGammaEventTree->fGenPhoE = genPhoton->pt * cosh(genPhoton->eta);  
            zmumuGammaEventTree->fGenPhoEt = genPhoton->pt; 
            zmumuGammaEventTree->fGenPhoEta = genPhoton->eta; 
            zmumuGammaEventTree->fGenPhoPhi = genPhoton->phi; 
            zmumuGammaEventTree->fGenPhoMotherPdgId = genPhoton->motherPdgID; 

            if (fabs(genPhoton->motherPdgID) == 13) {
              zmumuGammaEventTree->fIsFSR = kTRUE;
            } else {
              zmumuGammaEventTree->fIsFSR = kFALSE;
            }
            if ( (fabs(genPhoton->motherPdgID) >= 1 &&  fabs(genPhoton->motherPdgID) <= 6) ||  fabs(genPhoton->motherPdgID) == 21)  {
              zmumuGammaEventTree->fIsISR = kTRUE;
            } else {
              zmumuGammaEventTree->fIsISR = kFALSE;
            }
          }

          //Do we need these? Don't need them for now
          zmumuGammaEventTree->fPhoIEtaX = 0;       
          zmumuGammaEventTree->fPhoIPhiY = 0; 
          zmumuGammaEventTree->fMuNearIEtaX = 0; 
          zmumuGammaEventTree->fMuNearIPhiY = 0; 
          zmumuGammaEventTree->fMuNearIsEB = 0; 
          zmumuGammaEventTree->fMuNearIndex = 0; 

          zmumuGammaEventTree->tree_->Fill();
          if ((vmu1+vmu2).M() > 50 && (vmu1+vmu2+vpho).M() > 75 && (vmu1+vmu2+vpho).M() < 105 && (vmu1+vmu2).M()+(vmu1+vmu2+vpho).M() < 180) {
            nEvents++;
          }


        } //loop over photon

      } //loop over 2nd muon
    } //loop over first muon

  } //loop over events


    // Cleaning up the evaluators
  delete infile;
  infile=0, eventTree=0;    

  delete info;
  delete gen;
  delete muonArr;


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
