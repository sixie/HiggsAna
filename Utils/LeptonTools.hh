#ifndef HIGGSANA_UTILS_LEPTONTOOLS_HH
#define HIGGSANA_UTILS_LEPTONTOOLS_HH

#include <TMath.h>
#include <TClonesArray.h>          
#include "HiggsAna/DataTree/interface/TElectron.hh"
#include "HiggsAna/DataTree/interface/TMuon.hh"
#include "HiggsAna/DataTree/interface/TPFCandidate.hh"
#include "HiggsAna/DataTree/interface/TGenParticle.hh"
#include "CITCommon/CommonData/interface/ElectronTree.h"
#include "CITCommon/CommonData/interface/MuonTree.h"
#include "HiggsAna/Utils/LeptonIDCuts.hh"
#include "HiggsAna/DataTree/interface/HiggsAnaDefs.hh"
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"

//**************************************************************************************************************
//Function Definitions
//**************************************************************************************************************
Double_t EvaluateEGammaNonTriggeringMVA( const higgsana::TElectron *ele, EGammaMvaEleEstimator* eleMVAEstimator,
                                         Bool_t printDebug = kFALSE);  


//**************************************************************************************************************
//Function Implementations
//**************************************************************************************************************
Double_t EvaluateEGammaNonTriggeringMVA( const higgsana::TElectron *ele, EGammaMvaEleEstimator* eleMVAEstimator,
                                Bool_t printDebug) {


  double _fbrem = (ele->fBrem < -1. ) ? -1. : ele->fBrem;

  double _kftrk_chisq = ele->KFTrackChi2OverNdof;
  if( _kftrk_chisq > 10. ) _kftrk_chisq = 10.;

  double _kftrk_nhits = ele->KFTrackNLayersWithMeasurement; 

  double _gsftrk_chisq = ele->GsfTrackChi2OverNdof;
  if( _gsftrk_chisq > 200. ) _gsftrk_chisq = 200.;
 
  double _deta = fabs(ele->deltaEtaIn);
  if( _deta > 0.06 ) _deta = 0.06;

  double _dphi = fabs(ele->deltaPhiIn);
  if( _dphi > 0.6 ) _dphi = 0.6;

  double _detacalo = fabs(ele->dEtaCalo);
  if( _detacalo > 0.2 ) _detacalo = 0.2;

  double _sigieie = ele->sigiEtaiEta;

  double _sigiphiiphi = ele->sigiPhiiPhi;
  if( isnan(_sigiphiiphi ) ) _sigiphiiphi = 0.;

  double _etawidth = ele->SCEtaWidth; 

  double _phiwidth = ele->SCPhiWidth;

  double _e1x5e5x5 = 1 - ele->SeedE1x5OverE5x5;
  if( _e1x5e5x5 < -1. ) _e1x5e5x5 = -1.;
  if( _e1x5e5x5 >  2. ) _e1x5e5x5 =  2.;

  double _r9 = (ele->R9 > 5 ) ? 5 : ele->R9;

  double _h_o_e = ele->HoverE;

  double _e_o_p = (ele->EOverP > 20. ) ? 20 : ele->EOverP;

  double _eeleclu_o_pout = (ele->EEleClusterOverPout > 20. ) ? 20. : ele->EEleClusterOverPout;


  double electronP = (ele->pt*TMath::CosH(ele->eta)); //I think this is right
  double _IoEmIoP =  (double)(1.0/ele->EcalEnergy) - (double)(1.0/electronP);

  double _epreoraw = ele->PreShowerOverRaw;
  
  double mvaval = eleMVAEstimator->mvaValue( _fbrem,
					 _kftrk_chisq,
					 _kftrk_nhits,
					 _gsftrk_chisq,
					 _deta,
					 _dphi,
					 _detacalo,
					 _sigieie,
					 _sigiphiiphi,
					 _etawidth,
					 _phiwidth,
					 _e1x5e5x5,
					 _r9,
					 _h_o_e,
					 _e_o_p,
					 _IoEmIoP,
					 _eeleclu_o_pout,
					 _epreoraw,
					 ele->scEta,	
					 ele->pt,
					 printDebug );
  
  return mvaval;

}



void FillElectronTree(citana::ElectronTree *eleTree,
                      const higgsana::TElectron *ele, 
                      Int_t eleIndex,
                      TClonesArray *pfCandidates, 
                      Double_t rho, UInt_t DataEra,
                      UInt_t NPV,
                      UInt_t runNum,
                      UInt_t lumiSec,
                      UInt_t evtNum,
                      Float_t weight = 1.0
                      
  ) {

  eleTree->fWeight = weight;
  eleTree->fRunNumber = runNum;
  eleTree->fLumiSectionNumber = lumiSec;
  eleTree->fEventNumber = evtNum;
  eleTree->fEleEventNumberParity = (evtNum % 2 == 0);
  eleTree->fElePt = ele->pt; 
  eleTree->fEleEta = ele->eta; 
  eleTree->fElePhi = ele->phi; 
  eleTree->fEleSCEt = ele->scEt; 
  eleTree->fEleSCEta = ele->scEta; 
  eleTree->fEleSCPhi = ele->scPhi; 
  eleTree->fEleEcalEnergy = ele->EcalEnergy; 
  eleTree->fEleIsEcalDriven = ele->isEcalDriven;
  eleTree->fEleTriggerBit = 0;
  eleTree->fRho = rho; 
  eleTree->fNVertices = NPV; 
  eleTree->fEleD0 = ele->d0; 
  eleTree->fEleDZ = ele->dz; 
  eleTree->fEleIP3d = ele->ip3d; 
  eleTree->fEleIP3dSig = ele->ip3dSig; 
  eleTree->fEleMatchedConversion = !passConversionVeto(ele->isConv);
  eleTree->fEleConvDCot = ele->partnerDeltaCot;
  eleTree->fEleConvDist = ele->partnerDist;
  eleTree->fEleNMissHits = ele->nExpHitsInner;
  eleTree->fEleNBrem = ele->nBrem; 

  if (ele->fBrem < -1) {
    eleTree->fEleFBrem = -1;
  } else {
    eleTree->fEleFBrem = ele->fBrem; 
  }

  if (ele->EOverP > 20) {
    eleTree->fEleEOverP = 20.0;
  } else {
    eleTree->fEleEOverP = ele->EOverP; 
  }

  eleTree->fEleESeedClusterOverPIn = ele->ESeedClusterOverPIn; 
  eleTree->fEleESeedClusterOverPout = ele->ESeedClusterOverPout; 

  if (ele->EEleClusterOverPout > 20) {
    eleTree->fEleEEleClusterOverPout = 20.0;
  } else {
    eleTree->fEleEEleClusterOverPout = ele->EEleClusterOverPout; 
  }

  eleTree->fEleOneOverEMinusOneOverP = (1.0/(ele->scEt*cosh(ele->scEta)) - 1.0 / ele->pIn);

  if ( fabs(ele->deltaEtaIn) > 0.6) {
    eleTree->fEleDEtaIn = 0.6;
  } else {
    eleTree->fEleDEtaIn = fabs(ele->deltaEtaIn);
  }

  eleTree->fEleDPhiIn = ele->deltaPhiIn; 

  if (fabs(ele->dEtaCalo) > 0.2) {
    eleTree->fEledEtaCalo = 0.2;
  } else {
    eleTree->fEledEtaCalo = fabs(ele->dEtaCalo);
  }

  eleTree->fEledPhiCalo = ele->dPhiCalo;
  eleTree->fEleSigmaIEtaIEta = ele->sigiEtaiEta; 
  if (!std::isnan(ele->sigiPhiiPhi)) {
    eleTree->fEleSigmaIPhiIPhi = ele->sigiPhiiPhi; 
  } else {
    eleTree->fEleSigmaIPhiIPhi = 0.0;
  }
  if (eleTree->fEleSigmaIEtaIEta*eleTree->fEleSigmaIPhiIPhi > 0) {
    eleTree->fEleSigmaIEtaIPhi = ele->CovIEtaIPhi/(eleTree->fEleSigmaIEtaIEta*eleTree->fEleSigmaIPhiIPhi);
  } else if (ele->CovIEtaIPhi>0) {
    eleTree->fEleSigmaIEtaIPhi = 1.0; 
  } else {
    eleTree->fEleSigmaIEtaIPhi = -1.0; 
  }
  eleTree->fEleSCEtaWidth = ele->SCEtaWidth;
  eleTree->fEleSCPhiWidth = ele->SCPhiWidth;

  if (ele->R9 > 5) {
    eleTree->fEleR9 = 5.0;
  } else {
    eleTree->fEleR9 = ele->R9;
  }

  eleTree->fElePreShowerOverRaw = ele->PreShowerOverRaw;
  eleTree->fEleHoverE = ele->HoverE; 

  if (ele->GsfTrackChi2OverNdof > 200) {
    eleTree->fEleGsfTrackChi2OverNdof = 200;
  } else {
    eleTree->fEleGsfTrackChi2OverNdof = ele->GsfTrackChi2OverNdof;
  }
  if (ele->KFTrackChi2OverNdof > 10) {
    eleTree->fEleKFTrackChi2OverNDoF = 10;
  } else {
    eleTree->fEleKFTrackChi2OverNDoF = ele->KFTrackChi2OverNdof;
  }
  eleTree->fEleKFTrackNHits = ele->KFTrackNHits;
  eleTree->fEleKFTrackNLayersWithMeasurement = ele->KFTrackNLayersWithMeasurement;
  if (ele->SeedE1x5OverE5x5 >= 0) {
    eleTree->fEleOneMinusSeedE1x5OverE5x5 = 1.0 - ele->SeedE1x5OverE5x5;
    if (1.0 - ele->SeedE1x5OverE5x5 > 2.0) {
      eleTree->fEleOneMinusSeedE1x5OverE5x5 = 2.0;
    }
  } else {
    eleTree->fEleOneMinusSeedE1x5OverE5x5 = -1.0;
  }
  eleTree->fElePFMVA = ele->mva;
  eleTree->fEleTrkIso03 = ele->trkIso03; 
  eleTree->fEleEMIso03 = ele->emIso03; 
  eleTree->fEleHadIso03 = ele->hadIso03; 
  eleTree->fEleTrkIso04 = ele->trkIso04; 
  eleTree->fEleEMIso04 = ele->emIso04; 
  eleTree->fEleHadIso04 = ele->hadIso04; 
  eleTree->fElePFIso04 = ComputeElePFIso04( ele, eleIndex, pfCandidates, rho, DataEra); 
  eleTree->fChargedIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.0, 0.1)/ele->pt;
  eleTree->fChargedIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.1, 0.2)/ele->pt;
  eleTree->fChargedIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.2, 0.3)/ele->pt;
  eleTree->fChargedIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.3, 0.4)/ele->pt;
  eleTree->fChargedIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.4, 0.5)/ele->pt;
  eleTree->fGammaIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.0, 0.1)/ele->pt;
  eleTree->fGammaIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.1, 0.2)/ele->pt;
  eleTree->fGammaIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.2, 0.3)/ele->pt;
  eleTree->fGammaIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.3, 0.4)/ele->pt;
  eleTree->fGammaIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.4, 0.5)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.0, 0.1)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.1, 0.2)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.2, 0.3)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.3, 0.4)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.4, 0.5)/ele->pt;
  eleTree->fElePassTriggerDenominator = passElectronTriggerDenominator(ele);

  //**********************************************
  //Fill Variables for Regression
  //**********************************************
  eleTree->fIsEB = ele->isEB;
  eleTree->fIsEE = ele->isEE;  //temporary fix for my bug in filling ele->isEE.
  
  eleTree->fSCRawEnergy = ele->SCRawEnergy ;
  eleTree->fNClusters = ele->nBrem+1 ;
  eleTree->fEtaSeed = ele->EtaSeed ;
  eleTree->fPhiSeed = ele->PhiSeed ;
  eleTree->fESeed = ele->ESeed ;
  eleTree->fE3x3Seed = ele->E3x3Seed ;
  eleTree->fE5x5Seed = ele->E5x5Seed ;
  eleTree->fEMaxSeed = ele->EMaxSeed ;
  eleTree->fE2ndSeed = ele->E2ndSeed ;
  eleTree->fETopSeed = ele->ETopSeed ;
  eleTree->fEBottomSeed = ele->EBottomSeed ;
  eleTree->fELeftSeed = ele->ELeftSeed ;
  eleTree->fERightSeed = ele->ERightSeed ;
  eleTree->fE2x5MaxSeed = ele->E2x5MaxSeed ;
  eleTree->fE2x5TopSeed = ele->E2x5TopSeed ;
  eleTree->fE2x5BottomSeed = ele->E2x5BottomSeed ;
  eleTree->fE2x5LeftSeed = ele->E2x5LeftSeed ;
  eleTree->fE2x5RightSeed = ele->E2x5RightSeed ;
  eleTree->fIEtaSeed = ele->IEtaSeed ;
  eleTree->fIPhiSeed = ele->IPhiSeed ;
  eleTree->fEtaCrySeed = ele->EtaCrySeed ;
  eleTree->fPhiCrySeed = ele->PhiCrySeed ;
  eleTree->fEcalEnergyError = ele->EcalEnergyError;
  eleTree->fGsfTrackPIn = ele->pIn ;
  eleTree->fCharge = ele->q ;
  eleTree->fTrackMomentumError = fmin(ele->TrackMomentumError,500.0);
  eleTree->fEleClassification = ele->Classification;

  //**********************************************
  //Make final consistency checks before filling
  //**********************************************
  if (TMath::IsNaN(ele->sigiPhiiPhi)) {
    cout << "Problem with sigiPhiiPhi: NaN. exit.\n";
    assert(0);
  }

  //***********************
  //Fill Electron
  //***********************
  eleTree->tree_->Fill();
}




void FillElectronTree(citana::ElectronTree *eleTree,
                      const higgsana::TElectron *ele, 
                      Int_t eleIndex,
                      TClonesArray *pfCandidates, 
                      vector<const higgsana::TPFCandidate*> particlesToVeto,
                      Double_t rho, UInt_t DataEra,
                      UInt_t NPV,
                      UInt_t runNum,
                      UInt_t lumiSec,
                      UInt_t evtNum,
                      Float_t weight = 1.0
                      
  ) {

  eleTree->fWeight = weight;
  eleTree->fRunNumber = runNum;
  eleTree->fLumiSectionNumber = lumiSec;
  eleTree->fEventNumber = evtNum;
  eleTree->fEleEventNumberParity = (evtNum % 2 == 0);
  eleTree->fElePt = ele->pt; 
  eleTree->fEleEta = ele->eta; 
  eleTree->fElePhi = ele->phi; 
  eleTree->fEleSCEt = ele->scEt; 
  eleTree->fEleSCEta = ele->scEta; 
  eleTree->fEleSCPhi = ele->scPhi; 
  eleTree->fEleEcalEnergy = ele->EcalEnergy; 
  eleTree->fEleIsEcalDriven = ele->isEcalDriven;
  eleTree->fEleTriggerBit = 0;
  eleTree->fRho = rho; 
  eleTree->fNVertices = NPV; 
  eleTree->fEleD0 = ele->d0; 
  eleTree->fEleDZ = ele->dz; 
  eleTree->fEleIP3d = ele->ip3d; 
  eleTree->fEleIP3dSig = ele->ip3dSig; 
  eleTree->fEleMatchedConversion = !passConversionVeto(ele->isConv);
  eleTree->fEleConvDCot = ele->partnerDeltaCot;
  eleTree->fEleConvDist = ele->partnerDist;
  eleTree->fEleNMissHits = ele->nExpHitsInner;
  eleTree->fEleNBrem = ele->nBrem; 

  if (ele->fBrem < -1) {
    eleTree->fEleFBrem = -1;
  } else {
    eleTree->fEleFBrem = ele->fBrem; 
  }

  if (ele->EOverP > 20) {
    eleTree->fEleEOverP = 20.0;
  } else {
    eleTree->fEleEOverP = ele->EOverP; 
  }

  eleTree->fEleESeedClusterOverPIn = ele->ESeedClusterOverPIn; 
  eleTree->fEleESeedClusterOverPout = ele->ESeedClusterOverPout; 

  if (ele->EEleClusterOverPout > 20) {
    eleTree->fEleEEleClusterOverPout = 20.0;
  } else {
    eleTree->fEleEEleClusterOverPout = ele->EEleClusterOverPout; 
  }

  eleTree->fEleOneOverEMinusOneOverP = (1.0/ele->EcalEnergy - 1.0 / ele->pIn);

  if ( fabs(ele->deltaEtaIn) > 0.6) {
    eleTree->fEleDEtaIn = 0.6;
  } else {
    eleTree->fEleDEtaIn = fabs(ele->deltaEtaIn);
  }

  eleTree->fEleDPhiIn = ele->deltaPhiIn; 

  if (fabs(ele->dEtaCalo) > 0.2) {
    eleTree->fEledEtaCalo = 0.2;
  } else {
    eleTree->fEledEtaCalo = ele->dEtaCalo;
  }

  eleTree->fEledPhiCalo = ele->dPhiCalo;
  eleTree->fEleSigmaIEtaIEta = ele->sigiEtaiEta; 
  if (!std::isnan(ele->sigiPhiiPhi)) {
    eleTree->fEleSigmaIPhiIPhi = ele->sigiPhiiPhi; 
  } else {
    eleTree->fEleSigmaIPhiIPhi = 0.0;
  }
  if (eleTree->fEleSigmaIEtaIEta*eleTree->fEleSigmaIPhiIPhi > 0) {
    eleTree->fEleSigmaIEtaIPhi = ele->CovIEtaIPhi/(eleTree->fEleSigmaIEtaIEta*eleTree->fEleSigmaIPhiIPhi);
  } else if (ele->CovIEtaIPhi>0) {
    eleTree->fEleSigmaIEtaIPhi = 1.0; 
  } else {
    eleTree->fEleSigmaIEtaIPhi = -1.0; 
  }
  eleTree->fEleSCEtaWidth = ele->SCEtaWidth;
  eleTree->fEleSCPhiWidth = ele->SCPhiWidth;

  if (ele->R9 > 5) {
    eleTree->fEleR9 = 5.0;
  } else {
    eleTree->fEleR9 = ele->R9;
  }

  eleTree->fElePreShowerOverRaw = ele->PreShowerOverRaw;
  eleTree->fEleHoverE = ele->HoverE; 

  if (ele->GsfTrackChi2OverNdof > 200) {
    eleTree->fEleGsfTrackChi2OverNdof = 200;
  } else {
    eleTree->fEleGsfTrackChi2OverNdof = ele->GsfTrackChi2OverNdof;
  }
  if (ele->KFTrackChi2OverNdof > 10) {
    eleTree->fEleKFTrackChi2OverNDoF = 10;
  } else {
    eleTree->fEleKFTrackChi2OverNDoF = ele->KFTrackChi2OverNdof;
  }
  eleTree->fEleKFTrackNHits = ele->KFTrackNHits;
  eleTree->fEleKFTrackNLayersWithMeasurement = ele->KFTrackNLayersWithMeasurement;
  if (ele->SeedE1x5OverE5x5 >= 0) {
    eleTree->fEleOneMinusSeedE1x5OverE5x5 = 1.0 - ele->SeedE1x5OverE5x5;
    if (1.0 - ele->SeedE1x5OverE5x5 > 2.0) {
      eleTree->fEleOneMinusSeedE1x5OverE5x5 = 2.0;
    }
  } else {
    eleTree->fEleOneMinusSeedE1x5OverE5x5 = -1.0;
  }
  eleTree->fElePFMVA = ele->mva;
  eleTree->fEleTrkIso03 = ele->trkIso03; 
  eleTree->fEleEMIso03 = ele->emIso03; 
  eleTree->fEleHadIso03 = ele->hadIso03; 
  eleTree->fEleTrkIso04 = ele->trkIso04; 
  eleTree->fEleEMIso04 = ele->emIso04; 
  eleTree->fEleHadIso04 = ele->hadIso04; 
  eleTree->fElePFIso04 = ComputeElePFIso04( ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra); 
  eleTree->fChargedIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.0, 0.1)/ele->pt;
  eleTree->fChargedIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.1, 0.2)/ele->pt;
  eleTree->fChargedIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.2, 0.3)/ele->pt;
  eleTree->fChargedIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.3, 0.4)/ele->pt;
  eleTree->fChargedIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.4, 0.5)/ele->pt;
  eleTree->fGammaIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.0, 0.1)/ele->pt;
  eleTree->fGammaIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.1, 0.2)/ele->pt;
  eleTree->fGammaIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.2, 0.3)/ele->pt;
  eleTree->fGammaIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.3, 0.4)/ele->pt;
  eleTree->fGammaIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.4, 0.5)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.0, 0.1)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.1, 0.2)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.2, 0.3)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.3, 0.4)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.4, 0.5)/ele->pt;
  eleTree->fElePassTriggerDenominator = passElectronTriggerDenominator(ele);

  //**********************************************
  //Fill Variables for Regression
  //**********************************************
  eleTree->fIsEB = ele->isEB;
  eleTree->fIsEE = Bool_t(fabs(ele->scEta) > 1.479);  //temporary fix for my bug in filling ele->isEE.
  
  eleTree->fSCRawEnergy = ele->SCRawEnergy ;
  eleTree->fNClusters = ele->nBrem+1 ;
  eleTree->fEtaSeed = ele->EtaSeed ;
  eleTree->fPhiSeed = ele->PhiSeed ;
  eleTree->fESeed = ele->ESeed ;
  eleTree->fE3x3Seed = ele->E3x3Seed ;
  eleTree->fE5x5Seed = ele->E5x5Seed ;
  eleTree->fEMaxSeed = ele->EMaxSeed ;
  eleTree->fE2ndSeed = ele->E2ndSeed ;
  eleTree->fETopSeed = ele->ETopSeed ;
  eleTree->fEBottomSeed = ele->EBottomSeed ;
  eleTree->fELeftSeed = ele->ELeftSeed ;
  eleTree->fERightSeed = ele->ERightSeed ;
  eleTree->fE2x5MaxSeed = ele->E2x5MaxSeed ;
  eleTree->fE2x5TopSeed = ele->E2x5TopSeed ;
  eleTree->fE2x5BottomSeed = ele->E2x5BottomSeed ;
  eleTree->fE2x5LeftSeed = ele->E2x5LeftSeed ;
  eleTree->fE2x5RightSeed = ele->E2x5RightSeed ;
  eleTree->fIEtaSeed = ele->IEtaSeed ;
  eleTree->fIPhiSeed = ele->IPhiSeed ;
  eleTree->fEtaCrySeed = ele->EtaCrySeed ;
  eleTree->fPhiCrySeed = ele->PhiCrySeed ;
  eleTree->fEcalEnergyError = ele->EcalEnergyError;
  eleTree->fGsfTrackPIn = ele->pIn ;
  eleTree->fCharge = ele->q ;
  eleTree->fTrackMomentumError = fmin(ele->TrackMomentumError,500.0);
  eleTree->fEleClassification = ele->Classification;


  //**********************************************
  //Make final consistency checks before filling
  //**********************************************
  if (TMath::IsNaN(ele->sigiPhiiPhi)) {
    cout << "Problem with sigiPhiiPhi: NaN. exit.\n";
    assert(0);
  }

  //***********************
  //Fill Electron
  //***********************
  eleTree->tree_->Fill();
}


void FillElectronTree(citana::ElectronTree *eleTree,
                      const higgsana::TElectron *ele, 
                      Int_t eleIndex,
                      TClonesArray *pfCandidates, 
                      TClonesArray *genParticles, 
                      Double_t rho, UInt_t DataEra,
                      UInt_t NPV,
                      UInt_t runNum,
                      UInt_t lumiSec,
                      UInt_t evtNum,
                      Float_t weight = 1.0
                      
  ) {

  eleTree->fWeight = weight;
  eleTree->fRunNumber = runNum;
  eleTree->fLumiSectionNumber = lumiSec;
  eleTree->fEventNumber = evtNum;
  eleTree->fEleEventNumberParity = (evtNum % 2 == 0);
  eleTree->fElePt = ele->pt; 
  eleTree->fEleEta = ele->eta; 
  eleTree->fElePhi = ele->phi; 
  eleTree->fEleSCEt = ele->scEt; 
  eleTree->fEleSCEta = ele->scEta; 
  eleTree->fEleSCPhi = ele->scPhi; 
  eleTree->fEleEcalEnergy = ele->EcalEnergy; 
  eleTree->fEleIsEcalDriven = ele->isEcalDriven;
  eleTree->fEleTriggerBit = 0;
  eleTree->fRho = rho; 
  eleTree->fNVertices = NPV; 
  eleTree->fEleD0 = ele->d0; 
  eleTree->fEleDZ = ele->dz; 
  eleTree->fEleIP3d = ele->ip3d; 
  eleTree->fEleIP3dSig = ele->ip3dSig; 
  eleTree->fEleMatchedConversion = !passConversionVeto(ele->isConv);
  eleTree->fEleConvDCot = ele->partnerDeltaCot;
  eleTree->fEleConvDist = ele->partnerDist;
  eleTree->fEleNMissHits = ele->nExpHitsInner;
  eleTree->fEleNBrem = ele->nBrem; 

  if (ele->fBrem < -1) {
    eleTree->fEleFBrem = -1;
  } else {
    eleTree->fEleFBrem = ele->fBrem; 
  }

  if (ele->EOverP > 20) {
    eleTree->fEleEOverP = 20.0;
  } else {
    eleTree->fEleEOverP = ele->EOverP; 
  }

  eleTree->fEleESeedClusterOverPIn = ele->ESeedClusterOverPIn; 
  eleTree->fEleESeedClusterOverPout = ele->ESeedClusterOverPout; 

  if (ele->EEleClusterOverPout > 20) {
    eleTree->fEleEEleClusterOverPout = 20.0;
  } else {
    eleTree->fEleEEleClusterOverPout = ele->EEleClusterOverPout; 
  }

  eleTree->fEleOneOverEMinusOneOverP = (1.0/ele->EcalEnergy - 1.0 / ele->pIn);

  if ( fabs(ele->deltaEtaIn) > 0.6) {
    eleTree->fEleDEtaIn = 0.6;
  } else {
    eleTree->fEleDEtaIn = fabs(ele->deltaEtaIn);
  }

  eleTree->fEleDPhiIn = ele->deltaPhiIn; 

  if (fabs(ele->dEtaCalo) > 0.2) {
    eleTree->fEledEtaCalo = 0.2;
  } else {
    eleTree->fEledEtaCalo = ele->dEtaCalo;
  }

  eleTree->fEledPhiCalo = ele->dPhiCalo;
  eleTree->fEleSigmaIEtaIEta = ele->sigiEtaiEta; 
  if (!std::isnan(ele->sigiPhiiPhi)) {
    eleTree->fEleSigmaIPhiIPhi = ele->sigiPhiiPhi; 
  } else {
    eleTree->fEleSigmaIPhiIPhi = 0.0;
  }
  if (eleTree->fEleSigmaIEtaIEta*eleTree->fEleSigmaIPhiIPhi > 0) {
    eleTree->fEleSigmaIEtaIPhi = ele->CovIEtaIPhi/(eleTree->fEleSigmaIEtaIEta*eleTree->fEleSigmaIPhiIPhi);
  } else if (ele->CovIEtaIPhi>0) {
    eleTree->fEleSigmaIEtaIPhi = 1.0; 
  } else {
    eleTree->fEleSigmaIEtaIPhi = -1.0; 
  }
  eleTree->fEleSCEtaWidth = ele->SCEtaWidth;
  eleTree->fEleSCPhiWidth = ele->SCPhiWidth;

  if (ele->R9 > 5) {
    eleTree->fEleR9 = 5.0;
  } else {
    eleTree->fEleR9 = ele->R9;
  }

  eleTree->fElePreShowerOverRaw = ele->PreShowerOverRaw;
  eleTree->fEleHoverE = ele->HoverE; 

  if (ele->GsfTrackChi2OverNdof > 200) {
    eleTree->fEleGsfTrackChi2OverNdof = 200;
  } else {
    eleTree->fEleGsfTrackChi2OverNdof = ele->GsfTrackChi2OverNdof;
  }
  if (ele->KFTrackChi2OverNdof > 10) {
    eleTree->fEleKFTrackChi2OverNDoF = 10;
  } else {
    eleTree->fEleKFTrackChi2OverNDoF = ele->KFTrackChi2OverNdof;
  }
  eleTree->fEleKFTrackNHits = ele->KFTrackNHits;
  eleTree->fEleKFTrackNLayersWithMeasurement = ele->KFTrackNLayersWithMeasurement;
  if (ele->SeedE1x5OverE5x5 >= 0) {
    eleTree->fEleOneMinusSeedE1x5OverE5x5 = 1.0 - ele->SeedE1x5OverE5x5;
    if (1.0 - ele->SeedE1x5OverE5x5 > 2.0) {
      eleTree->fEleOneMinusSeedE1x5OverE5x5 = 2.0;
    }
  } else {
    eleTree->fEleOneMinusSeedE1x5OverE5x5 = -1.0;
  }
  eleTree->fElePFMVA = ele->mva;
  eleTree->fEleTrkIso03 = ele->trkIso03; 
  eleTree->fEleEMIso03 = ele->emIso03; 
  eleTree->fEleHadIso03 = ele->hadIso03; 
  eleTree->fEleTrkIso04 = ele->trkIso04; 
  eleTree->fEleEMIso04 = ele->emIso04; 
  eleTree->fEleHadIso04 = ele->hadIso04; 
  eleTree->fElePFIso04 = ComputeElePFIso04( ele, eleIndex, pfCandidates, rho, DataEra); 
  eleTree->fChargedIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.0, 0.1)/ele->pt;
  eleTree->fChargedIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.1, 0.2)/ele->pt;
  eleTree->fChargedIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.2, 0.3)/ele->pt;
  eleTree->fChargedIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.3, 0.4)/ele->pt;
  eleTree->fChargedIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.4, 0.5)/ele->pt;
  eleTree->fGammaIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.0, 0.1)/ele->pt;
  eleTree->fGammaIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.1, 0.2)/ele->pt;
  eleTree->fGammaIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.2, 0.3)/ele->pt;
  eleTree->fGammaIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.3, 0.4)/ele->pt;
  eleTree->fGammaIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.4, 0.5)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.0, 0.1)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.1, 0.2)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.2, 0.3)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.3, 0.4)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.4, 0.5)/ele->pt;
  eleTree->fElePassTriggerDenominator = passElectronTriggerDenominator(ele);


  //**********************************************
  //Fill Variables for Regression
  //**********************************************
  eleTree->fIsEB = ele->isEB;
  eleTree->fIsEE = Bool_t(fabs(ele->scEta) > 1.479);  //temporary fix for my bug in filling ele->isEE.
  
  eleTree->fSCRawEnergy = ele->SCRawEnergy ;
  eleTree->fNClusters = ele->nBrem+1 ;
  eleTree->fEtaSeed = ele->EtaSeed ;
  eleTree->fPhiSeed = ele->PhiSeed ;
  eleTree->fESeed = ele->ESeed ;
  eleTree->fE3x3Seed = ele->E3x3Seed ;
  eleTree->fE5x5Seed = ele->E5x5Seed ;
  eleTree->fEMaxSeed = ele->EMaxSeed ;
  eleTree->fE2ndSeed = ele->E2ndSeed ;
  eleTree->fETopSeed = ele->ETopSeed ;
  eleTree->fEBottomSeed = ele->EBottomSeed ;
  eleTree->fELeftSeed = ele->ELeftSeed ;
  eleTree->fERightSeed = ele->ERightSeed ;
  eleTree->fE2x5MaxSeed = ele->E2x5MaxSeed ;
  eleTree->fE2x5TopSeed = ele->E2x5TopSeed ;
  eleTree->fE2x5BottomSeed = ele->E2x5BottomSeed ;
  eleTree->fE2x5LeftSeed = ele->E2x5LeftSeed ;
  eleTree->fE2x5RightSeed = ele->E2x5RightSeed ;
  eleTree->fIEtaSeed = ele->IEtaSeed ;
  eleTree->fIPhiSeed = ele->IPhiSeed ;
  eleTree->fEtaCrySeed = ele->EtaCrySeed ;
  eleTree->fPhiCrySeed = ele->PhiCrySeed ;
  eleTree->fEcalEnergyError = ele->EcalEnergyError;
  eleTree->fGsfTrackPIn = ele->pIn ;
  eleTree->fCharge = ele->q ;
  eleTree->fTrackMomentumError = fmin(ele->TrackMomentumError,500.0);
  eleTree->fEleClassification = ele->Classification;

  //**********************************************
  //Find the gen level electrons
  //**********************************************
  double_t tmpGenEnergyStatus1 = -1;
  double_t tmpGenEnergyStatus3 = -1;
  double_t tmpGenEnergySupercluster = -1;

  double minDRStatus1 = 9999;
  Int_t matchedStatus1ElectronIndex = -1;
  double minDRStatus3 = 9999;  
  Int_t matchedStatus3ElectronIndex = -1;
  for(Int_t k=0; k<genParticles->GetEntries(); k++) {
    const higgsana::TGenParticle *gen = (higgsana::TGenParticle*)((*genParticles)[k]);

    //status 1 match
    if (abs(gen->pdgid) == 11 && (gen->status == 1)) {
      double tmpDR = higgsana::deltaR( ele->eta, ele->phi, gen->eta , gen->phi);
      if (tmpDR < minDRStatus1 && tmpDR < 0.3) {
        minDRStatus1 = tmpDR;
        TLorentzVector v1; 
        v1.SetPtEtaPhiM( gen->pt, gen->eta, gen->phi, ELECTRONMASS);
        tmpGenEnergyStatus1 = v1.E();
        matchedStatus1ElectronIndex = k;
//         cout << "status1 - Check energy calculation: " << gen->pt << " " << gen->eta << " : " << v1.E() << " ---> " << gen->pt / TMath::CosH(gen->eta) << endl;
      }
    }


    //status 3 match
    if (abs(gen->pdgid) == 11 && (gen->status == 3)) {
      double tmpDR = higgsana::deltaR( ele->eta, ele->phi, gen->eta , gen->phi);
      if (tmpDR < minDRStatus3 && tmpDR < 0.3) {
        minDRStatus3 = tmpDR;
        TLorentzVector v1; 
        v1.SetPtEtaPhiM( gen->pt, gen->eta, gen->phi, ELECTRONMASS);
        tmpGenEnergyStatus3 = v1.E();
        matchedStatus3ElectronIndex = k;
//         cout << "status 3 - Check energy calculation: " << gen->pt << " " << gen->eta << " : " << v1.E() << endl;
      }
    }
  }  

  

  TLorentzVector genElectronSupercluster;
  const higgsana::TGenParticle *genEleStatus1 = (higgsana::TGenParticle*)((*genParticles)[matchedStatus1ElectronIndex]);
  genElectronSupercluster.SetPtEtaPhiM( genEleStatus1->pt, genEleStatus1->eta, genEleStatus1->phi, ELECTRONMASS);

//   cout << "Status1 electron : " << genEleStatus1->pt << " " << genEleStatus1->eta << " " << genEleStatus1->phi << " , energy = " << genElectronSupercluster.E() << endl;

  for(Int_t k=0; k<genParticles->GetEntries(); k++) {
    const higgsana::TGenParticle *gen = (higgsana::TGenParticle*)((*genParticles)[k]);
    if (gen->pdgid == 22 && (gen->status == 1)) {
      //within deltaPhi road search
      if( fabs(gen->eta - genEleStatus1->eta) < 0.03 && higgsana::deltaPhi(gen->phi, genEleStatus1->phi) < 0.3) {
        TLorentzVector tmpPhoton;
        tmpPhoton.SetPtEtaPhiM( gen->pt, gen->eta, gen->phi, 0);
        TLorentzVector vtmp = genElectronSupercluster + tmpPhoton;
        genElectronSupercluster.SetPtEtaPhiM( vtmp.Pt(), vtmp.Eta(), vtmp.Phi(), vtmp.M());        
//         cout << "add photon : " << gen->pt << " " << gen->eta << " " << gen->phi << " -> " << genElectronSupercluster.Pt() << " " << genElectronSupercluster.Eta() << " " << genElectronSupercluster.Phi() << endl;
      }
    }
  }
//   cout << "Final GenLevel superclustered electron : " << genElectronSupercluster.Pt() << " " << genElectronSupercluster.Eta() << " " << genElectronSupercluster.Phi() << " , energy = " << genElectronSupercluster.E() << endl;
//   cout << endl;

     
  eleTree->fGeneratedEnergy = genElectronSupercluster.E() ;
  eleTree->fGeneratedEnergyStatus1 = tmpGenEnergyStatus1 ;
  eleTree->fGeneratedEnergyStatus3 = tmpGenEnergyStatus3 ;

  //debug
//   if (tmpGenEnergyStatus3 < 2.0) {
//     cout << "Electron : " << ele->pt << " " << ele->eta << " " << ele->phi << endl;
//     cout << "Matched Status1 : " << genEleStatus1->pt << " " << genEleStatus1->eta << " " << genEleStatus1->phi << endl;
//     for(Int_t k=0; k<genParticles->GetEntries(); k++) {
//       const higgsana::TGenParticle *gen = (higgsana::TGenParticle*)((*genParticles)[k]);
//       cout << "Particle " << k << " : " << gen->pdgid << " " << gen->status << " : " << gen->motherPdgID << " : " << gen->pt << " " << gen->eta << " " << gen->phi << endl;
//     }
//   }

  //**********************************************
  //Make final consistency checks before filling
  //**********************************************
  if (TMath::IsNaN(ele->sigiPhiiPhi)) {
    cout << "Problem with sigiPhiiPhi: NaN. exit.\n";
    assert(0);
  }

  //***********************
  //Fill Electron
  //***********************
  eleTree->tree_->Fill();
}








void FillElectronTree(citana::ElectronTree *eleTree,
                      const higgsana::TElectron *ele, 
                      Int_t eleIndex,
                      TClonesArray *pfCandidates, 
                      EGammaMvaEleEstimator *eleIDMVA,
                      Double_t rho, UInt_t DataEra,
                      UInt_t NPV,
                      UInt_t runNum,
                      UInt_t lumiSec,
                      UInt_t evtNum,
                      Float_t weight = 1.0
                      
  ) {

  eleTree->fWeight = weight;
  eleTree->fRunNumber = runNum;
  eleTree->fLumiSectionNumber = lumiSec;
  eleTree->fEventNumber = evtNum;
  eleTree->fEleEventNumberParity = (evtNum % 2 == 0);
  eleTree->fElePt = ele->pt; 
  eleTree->fEleEta = ele->eta; 
  eleTree->fElePhi = ele->phi; 
  eleTree->fEleSCEt = ele->scEt; 
  eleTree->fEleSCEta = ele->scEta; 
  eleTree->fEleSCPhi = ele->scPhi; 
  eleTree->fEleEcalEnergy = ele->EcalEnergy; 
  eleTree->fEleIsEcalDriven = ele->isEcalDriven;
  eleTree->fEleTriggerBit = 0;
  eleTree->fRho = rho; 
  eleTree->fNVertices = NPV; 
  eleTree->fEleD0 = ele->d0; 
  eleTree->fEleDZ = ele->dz; 
  eleTree->fEleIP3d = ele->ip3d; 
  eleTree->fEleIP3dSig = ele->ip3dSig; 
  eleTree->fEleMatchedConversion = !passConversionVeto(ele->isConv);
  eleTree->fEleConvDCot = ele->partnerDeltaCot;
  eleTree->fEleConvDist = ele->partnerDist;
  eleTree->fEleNMissHits = ele->nExpHitsInner;
  eleTree->fEleNBrem = ele->nBrem; 

  if (ele->fBrem < -1) {
    eleTree->fEleFBrem = -1;
  } else {
    eleTree->fEleFBrem = ele->fBrem; 
  }

  if (ele->EOverP > 20) {
    eleTree->fEleEOverP = 20.0;
  } else {
    eleTree->fEleEOverP = ele->EOverP; 
  }

  eleTree->fEleESeedClusterOverPIn = ele->ESeedClusterOverPIn; 
  eleTree->fEleESeedClusterOverPout = ele->ESeedClusterOverPout; 

  if (ele->EEleClusterOverPout > 20) {
    eleTree->fEleEEleClusterOverPout = 20.0;
  } else {
    eleTree->fEleEEleClusterOverPout = ele->EEleClusterOverPout; 
  }

  eleTree->fEleOneOverEMinusOneOverP = (1.0/ele->EcalEnergy - 1.0 / ele->pIn);

  if ( fabs(ele->deltaEtaIn) > 0.6) {
    eleTree->fEleDEtaIn = 0.6;
  } else {
    eleTree->fEleDEtaIn = fabs(ele->deltaEtaIn);
  }

  eleTree->fEleDPhiIn = ele->deltaPhiIn; 

  if (fabs(ele->dEtaCalo) > 0.2) {
    eleTree->fEledEtaCalo = 0.2;
  } else {
    eleTree->fEledEtaCalo = ele->dEtaCalo;
  }

  eleTree->fEledPhiCalo = ele->dPhiCalo;
  eleTree->fEleSigmaIEtaIEta = ele->sigiEtaiEta; 
  if (!std::isnan(ele->sigiPhiiPhi)) {
    eleTree->fEleSigmaIPhiIPhi = ele->sigiPhiiPhi; 
  } else {
    eleTree->fEleSigmaIPhiIPhi = 0.0;
  }
  if (eleTree->fEleSigmaIEtaIEta*eleTree->fEleSigmaIPhiIPhi > 0) {
    eleTree->fEleSigmaIEtaIPhi = ele->CovIEtaIPhi/(eleTree->fEleSigmaIEtaIEta*eleTree->fEleSigmaIPhiIPhi);
  } else if (ele->CovIEtaIPhi>0) {
    eleTree->fEleSigmaIEtaIPhi = 1.0; 
  } else {
    eleTree->fEleSigmaIEtaIPhi = -1.0; 
  }
  eleTree->fEleSCEtaWidth = ele->SCEtaWidth;
  eleTree->fEleSCPhiWidth = ele->SCPhiWidth;

  if (ele->R9 > 5) {
    eleTree->fEleR9 = 5.0;
  } else {
    eleTree->fEleR9 = ele->R9;
  }

  eleTree->fElePreShowerOverRaw = ele->PreShowerOverRaw;
  eleTree->fEleHoverE = ele->HoverE; 

  if (ele->GsfTrackChi2OverNdof > 200) {
    eleTree->fEleGsfTrackChi2OverNdof = 200;
  } else {
    eleTree->fEleGsfTrackChi2OverNdof = ele->GsfTrackChi2OverNdof;
  }
  if (ele->KFTrackChi2OverNdof > 10) {
    eleTree->fEleKFTrackChi2OverNDoF = 10;
  } else {
    eleTree->fEleKFTrackChi2OverNDoF = ele->KFTrackChi2OverNdof;
  }
  eleTree->fEleKFTrackNHits = ele->KFTrackNHits;
  eleTree->fEleKFTrackNLayersWithMeasurement = ele->KFTrackNLayersWithMeasurement;
  if (ele->SeedE1x5OverE5x5 >= 0) {
    eleTree->fEleOneMinusSeedE1x5OverE5x5 = 1.0 - ele->SeedE1x5OverE5x5;
    if (1.0 - ele->SeedE1x5OverE5x5 > 2.0) {
      eleTree->fEleOneMinusSeedE1x5OverE5x5 = 2.0;
    }
  } else {
    eleTree->fEleOneMinusSeedE1x5OverE5x5 = -1.0;
  }
  eleTree->fElePFMVA = ele->mva;
  eleTree->fEleTrkIso03 = ele->trkIso03; 
  eleTree->fEleEMIso03 = ele->emIso03; 
  eleTree->fEleHadIso03 = ele->hadIso03; 
  eleTree->fEleTrkIso04 = ele->trkIso04; 
  eleTree->fEleEMIso04 = ele->emIso04; 
  eleTree->fEleHadIso04 = ele->hadIso04; 
  eleTree->fElePFIso04 = ComputeElePFIso04( ele, eleIndex, pfCandidates, rho, DataEra); 
  eleTree->fChargedIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.0, 0.1)/ele->pt;
  eleTree->fChargedIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.1, 0.2)/ele->pt;
  eleTree->fChargedIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.2, 0.3)/ele->pt;
  eleTree->fChargedIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.3, 0.4)/ele->pt;
  eleTree->fChargedIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.4, 0.5)/ele->pt;
  eleTree->fGammaIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.0, 0.1)/ele->pt;
  eleTree->fGammaIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.1, 0.2)/ele->pt;
  eleTree->fGammaIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.2, 0.3)/ele->pt;
  eleTree->fGammaIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.3, 0.4)/ele->pt;
  eleTree->fGammaIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.4, 0.5)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.0, 0.1)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.1, 0.2)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.2, 0.3)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.3, 0.4)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.4, 0.5)/ele->pt;
  eleTree->fElePassTriggerDenominator = passElectronTriggerDenominator(ele);

  //**********************************************
  //Fill Variables for Regression
  //**********************************************
  eleTree->fIsEB = ele->isEB;
  eleTree->fIsEE = Bool_t(fabs(ele->scEta) > 1.479);  //temporary fix for my bug in filling ele->isEE.
  
  eleTree->fSCRawEnergy = ele->SCRawEnergy ;
  eleTree->fNClusters = ele->nBrem+1 ;
  eleTree->fEtaSeed = ele->EtaSeed ;
  eleTree->fPhiSeed = ele->PhiSeed ;
  eleTree->fESeed = ele->ESeed ;
  eleTree->fE3x3Seed = ele->E3x3Seed ;
  eleTree->fE5x5Seed = ele->E5x5Seed ;
  eleTree->fEMaxSeed = ele->EMaxSeed ;
  eleTree->fE2ndSeed = ele->E2ndSeed ;
  eleTree->fETopSeed = ele->ETopSeed ;
  eleTree->fEBottomSeed = ele->EBottomSeed ;
  eleTree->fELeftSeed = ele->ELeftSeed ;
  eleTree->fERightSeed = ele->ERightSeed ;
  eleTree->fE2x5MaxSeed = ele->E2x5MaxSeed ;
  eleTree->fE2x5TopSeed = ele->E2x5TopSeed ;
  eleTree->fE2x5BottomSeed = ele->E2x5BottomSeed ;
  eleTree->fE2x5LeftSeed = ele->E2x5LeftSeed ;
  eleTree->fE2x5RightSeed = ele->E2x5RightSeed ;
  eleTree->fIEtaSeed = ele->IEtaSeed ;
  eleTree->fIPhiSeed = ele->IPhiSeed ;
  eleTree->fEtaCrySeed = ele->EtaCrySeed ;
  eleTree->fPhiCrySeed = ele->PhiCrySeed ;
  eleTree->fEcalEnergyError = ele->EcalEnergyError;
  eleTree->fGsfTrackPIn = ele->pIn ;
  eleTree->fCharge = ele->q ;
  eleTree->fTrackMomentumError = fmin(ele->TrackMomentumError,500.0);
  eleTree->fEleClassification = ele->Classification;

  //**********************************************
  //MVA values
  //**********************************************
  eleTree->fEGammaNonTriggeringMVA = EvaluateEGammaNonTriggeringMVA(ele,eleIDMVA, false);

  //**********************************************
  //Make final consistency checks before filling
  //**********************************************
  if (TMath::IsNaN(ele->sigiPhiiPhi)) {
    cout << "Problem with sigiPhiiPhi: NaN. exit.\n";
    assert(0);
  }

  //***********************
  //Fill Electron
  //***********************
  eleTree->tree_->Fill();
}




void FillElectronTree(citana::ElectronTree *eleTree,
                      const higgsana::TElectron *ele, 
                      Int_t eleIndex,
                      TClonesArray *pfCandidates, 
                      vector<const higgsana::TPFCandidate*> particlesToVeto,
                      EGammaMvaEleEstimator *eleIDMVA,
                      Double_t rho, UInt_t DataEra,
                      UInt_t NPV,
                      UInt_t runNum,
                      UInt_t lumiSec,
                      UInt_t evtNum,
                      Float_t weight = 1.0
                      
  ) {

  eleTree->fWeight = weight;
  eleTree->fRunNumber = runNum;
  eleTree->fLumiSectionNumber = lumiSec;
  eleTree->fEventNumber = evtNum;
  eleTree->fEleEventNumberParity = (evtNum % 2 == 0);
  eleTree->fElePt = ele->pt; 
  eleTree->fEleEta = ele->eta; 
  eleTree->fElePhi = ele->phi; 
  eleTree->fEleSCEt = ele->scEt; 
  eleTree->fEleSCEta = ele->scEta; 
  eleTree->fEleSCPhi = ele->scPhi; 
  eleTree->fEleEcalEnergy = ele->EcalEnergy; 
  eleTree->fEleIsEcalDriven = ele->isEcalDriven;
  eleTree->fEleTriggerBit = 0;
  eleTree->fRho = rho; 
  eleTree->fNVertices = NPV; 
  eleTree->fEleD0 = ele->d0; 
  eleTree->fEleDZ = ele->dz; 
  eleTree->fEleIP3d = ele->ip3d; 
  eleTree->fEleIP3dSig = ele->ip3dSig; 
  eleTree->fEleMatchedConversion = !passConversionVeto(ele->isConv);
  eleTree->fEleConvDCot = ele->partnerDeltaCot;
  eleTree->fEleConvDist = ele->partnerDist;
  eleTree->fEleNMissHits = ele->nExpHitsInner;
  eleTree->fEleNBrem = ele->nBrem; 

  if (ele->fBrem < -1) {
    eleTree->fEleFBrem = -1;
  } else {
    eleTree->fEleFBrem = ele->fBrem; 
  }

  if (ele->EOverP > 20) {
    eleTree->fEleEOverP = 20.0;
  } else {
    eleTree->fEleEOverP = ele->EOverP; 
  }

  eleTree->fEleESeedClusterOverPIn = ele->ESeedClusterOverPIn; 
  eleTree->fEleESeedClusterOverPout = ele->ESeedClusterOverPout; 

  if (ele->EEleClusterOverPout > 20) {
    eleTree->fEleEEleClusterOverPout = 20.0;
  } else {
    eleTree->fEleEEleClusterOverPout = ele->EEleClusterOverPout; 
  }

  eleTree->fEleOneOverEMinusOneOverP = (1.0/ele->EcalEnergy - 1.0 / ele->pIn);

  if ( fabs(ele->deltaEtaIn) > 0.6) {
    eleTree->fEleDEtaIn = 0.6;
  } else {
    eleTree->fEleDEtaIn = fabs(ele->deltaEtaIn);
  }

  eleTree->fEleDPhiIn = ele->deltaPhiIn; 

  if (fabs(ele->dEtaCalo) > 0.2) {
    eleTree->fEledEtaCalo = 0.2;
  } else {
    eleTree->fEledEtaCalo = ele->dEtaCalo;
  }

  eleTree->fEledPhiCalo = ele->dPhiCalo;
  eleTree->fEleSigmaIEtaIEta = ele->sigiEtaiEta; 
  if (!std::isnan(ele->sigiPhiiPhi)) {
    eleTree->fEleSigmaIPhiIPhi = ele->sigiPhiiPhi; 
  } else {
    eleTree->fEleSigmaIPhiIPhi = 0.0;
  }
  if (eleTree->fEleSigmaIEtaIEta*eleTree->fEleSigmaIPhiIPhi > 0) {
    eleTree->fEleSigmaIEtaIPhi = ele->CovIEtaIPhi/(eleTree->fEleSigmaIEtaIEta*eleTree->fEleSigmaIPhiIPhi);
  } else if (ele->CovIEtaIPhi>0) {
    eleTree->fEleSigmaIEtaIPhi = 1.0; 
  } else {
    eleTree->fEleSigmaIEtaIPhi = -1.0; 
  }
  eleTree->fEleSCEtaWidth = ele->SCEtaWidth;
  eleTree->fEleSCPhiWidth = ele->SCPhiWidth;

  if (ele->R9 > 5) {
    eleTree->fEleR9 = 5.0;
  } else {
    eleTree->fEleR9 = ele->R9;
  }

  eleTree->fElePreShowerOverRaw = ele->PreShowerOverRaw;
  eleTree->fEleHoverE = ele->HoverE; 

  if (ele->GsfTrackChi2OverNdof > 200) {
    eleTree->fEleGsfTrackChi2OverNdof = 200;
  } else {
    eleTree->fEleGsfTrackChi2OverNdof = ele->GsfTrackChi2OverNdof;
  }
  if (ele->KFTrackChi2OverNdof > 10) {
    eleTree->fEleKFTrackChi2OverNDoF = 10;
  } else {
    eleTree->fEleKFTrackChi2OverNDoF = ele->KFTrackChi2OverNdof;
  }
  eleTree->fEleKFTrackNHits = ele->KFTrackNHits;
  eleTree->fEleKFTrackNLayersWithMeasurement = ele->KFTrackNLayersWithMeasurement;
  if (ele->SeedE1x5OverE5x5 >= 0) {
    eleTree->fEleOneMinusSeedE1x5OverE5x5 = 1.0 - ele->SeedE1x5OverE5x5;
    if (1.0 - ele->SeedE1x5OverE5x5 > 2.0) {
      eleTree->fEleOneMinusSeedE1x5OverE5x5 = 2.0;
    }
  } else {
    eleTree->fEleOneMinusSeedE1x5OverE5x5 = -1.0;
  }
  eleTree->fElePFMVA = ele->mva;
  eleTree->fEleTrkIso03 = ele->trkIso03; 
  eleTree->fEleEMIso03 = ele->emIso03; 
  eleTree->fEleHadIso03 = ele->hadIso03; 
  eleTree->fEleTrkIso04 = ele->trkIso04; 
  eleTree->fEleEMIso04 = ele->emIso04; 
  eleTree->fEleHadIso04 = ele->hadIso04; 
  eleTree->fElePFIso04 = ComputeElePFIso04( ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra); 
  eleTree->fChargedIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.0, 0.1)/ele->pt;
  eleTree->fChargedIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.1, 0.2)/ele->pt;
  eleTree->fChargedIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.2, 0.3)/ele->pt;
  eleTree->fChargedIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.3, 0.4)/ele->pt;
  eleTree->fChargedIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.4, 0.5)/ele->pt;
  eleTree->fGammaIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.0, 0.1)/ele->pt;
  eleTree->fGammaIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.1, 0.2)/ele->pt;
  eleTree->fGammaIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.2, 0.3)/ele->pt;
  eleTree->fGammaIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.3, 0.4)/ele->pt;
  eleTree->fGammaIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.4, 0.5)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.0, 0.1)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.1, 0.2)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.2, 0.3)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.3, 0.4)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.4, 0.5)/ele->pt;
  eleTree->fElePassTriggerDenominator = passElectronTriggerDenominator(ele);

  //**********************************************
  //Fill Variables for Regression
  //**********************************************
  eleTree->fIsEB = ele->isEB;
  eleTree->fIsEE = Bool_t(fabs(ele->scEta) > 1.479);  //temporary fix for my bug in filling ele->isEE.
  
  eleTree->fSCRawEnergy = ele->SCRawEnergy ;
  eleTree->fNClusters = ele->nBrem+1 ;
  eleTree->fEtaSeed = ele->EtaSeed ;
  eleTree->fPhiSeed = ele->PhiSeed ;
  eleTree->fESeed = ele->ESeed ;
  eleTree->fE3x3Seed = ele->E3x3Seed ;
  eleTree->fE5x5Seed = ele->E5x5Seed ;
  eleTree->fEMaxSeed = ele->EMaxSeed ;
  eleTree->fE2ndSeed = ele->E2ndSeed ;
  eleTree->fETopSeed = ele->ETopSeed ;
  eleTree->fEBottomSeed = ele->EBottomSeed ;
  eleTree->fELeftSeed = ele->ELeftSeed ;
  eleTree->fERightSeed = ele->ERightSeed ;
  eleTree->fE2x5MaxSeed = ele->E2x5MaxSeed ;
  eleTree->fE2x5TopSeed = ele->E2x5TopSeed ;
  eleTree->fE2x5BottomSeed = ele->E2x5BottomSeed ;
  eleTree->fE2x5LeftSeed = ele->E2x5LeftSeed ;
  eleTree->fE2x5RightSeed = ele->E2x5RightSeed ;
  eleTree->fIEtaSeed = ele->IEtaSeed ;
  eleTree->fIPhiSeed = ele->IPhiSeed ;
  eleTree->fEtaCrySeed = ele->EtaCrySeed ;
  eleTree->fPhiCrySeed = ele->PhiCrySeed ;
  eleTree->fEcalEnergyError = ele->EcalEnergyError;
  eleTree->fGsfTrackPIn = ele->pIn ;
  eleTree->fCharge = ele->q ;
  eleTree->fTrackMomentumError = fmin(ele->TrackMomentumError,500.0);
  eleTree->fEleClassification = ele->Classification;

  //**********************************************
  //MVA values
  //**********************************************
  eleTree->fEGammaNonTriggeringMVA = EvaluateEGammaNonTriggeringMVA(ele,eleIDMVA, false);


  //**********************************************
  //Make final consistency checks before filling
  //**********************************************
  if (TMath::IsNaN(ele->sigiPhiiPhi)) {
    cout << "Problem with sigiPhiiPhi: NaN. exit.\n";
    assert(0);
  }

  //***********************
  //Fill Electron
  //***********************
  eleTree->tree_->Fill();
}


void FillElectronTree(citana::ElectronTree *eleTree,
                      const higgsana::TElectron *ele, 
                      Int_t eleIndex,
                      TClonesArray *pfCandidates, 
                      EGammaMvaEleEstimator *eleIDMVA,
                      TClonesArray *genParticles, 
                      Double_t rho, UInt_t DataEra,
                      UInt_t NPV,
                      UInt_t runNum,
                      UInt_t lumiSec,
                      UInt_t evtNum,
                      Float_t weight = 1.0
                      
  ) {

  eleTree->fWeight = weight;
  eleTree->fRunNumber = runNum;
  eleTree->fLumiSectionNumber = lumiSec;
  eleTree->fEventNumber = evtNum;
  eleTree->fEleEventNumberParity = (evtNum % 2 == 0);
  eleTree->fElePt = ele->pt; 
  eleTree->fEleEta = ele->eta; 
  eleTree->fElePhi = ele->phi; 
  eleTree->fEleSCEt = ele->scEt; 
  eleTree->fEleSCEta = ele->scEta; 
  eleTree->fEleSCPhi = ele->scPhi; 
  eleTree->fEleEcalEnergy = ele->EcalEnergy; 
  eleTree->fEleIsEcalDriven = ele->isEcalDriven;
  eleTree->fEleTriggerBit = 0;
  eleTree->fRho = rho; 
  eleTree->fNVertices = NPV; 
  eleTree->fEleD0 = ele->d0; 
  eleTree->fEleDZ = ele->dz; 
  eleTree->fEleIP3d = ele->ip3d; 
  eleTree->fEleIP3dSig = ele->ip3dSig; 
  eleTree->fEleMatchedConversion = !passConversionVeto(ele->isConv);
  eleTree->fEleConvDCot = ele->partnerDeltaCot;
  eleTree->fEleConvDist = ele->partnerDist;
  eleTree->fEleNMissHits = ele->nExpHitsInner;
  eleTree->fEleNBrem = ele->nBrem; 

  if (ele->fBrem < -1) {
    eleTree->fEleFBrem = -1;
  } else {
    eleTree->fEleFBrem = ele->fBrem; 
  }

  if (ele->EOverP > 20) {
    eleTree->fEleEOverP = 20.0;
  } else {
    eleTree->fEleEOverP = ele->EOverP; 
  }

  eleTree->fEleESeedClusterOverPIn = ele->ESeedClusterOverPIn; 
  eleTree->fEleESeedClusterOverPout = ele->ESeedClusterOverPout; 

  if (ele->EEleClusterOverPout > 20) {
    eleTree->fEleEEleClusterOverPout = 20.0;
  } else {
    eleTree->fEleEEleClusterOverPout = ele->EEleClusterOverPout; 
  }

  eleTree->fEleOneOverEMinusOneOverP = (1.0/ele->EcalEnergy - 1.0 / ele->pIn);

  if ( fabs(ele->deltaEtaIn) > 0.6) {
    eleTree->fEleDEtaIn = 0.6;
  } else {
    eleTree->fEleDEtaIn = fabs(ele->deltaEtaIn);
  }

  eleTree->fEleDPhiIn = ele->deltaPhiIn; 

  if (fabs(ele->dEtaCalo) > 0.2) {
    eleTree->fEledEtaCalo = 0.2;
  } else {
    eleTree->fEledEtaCalo = ele->dEtaCalo;
  }

  eleTree->fEledPhiCalo = ele->dPhiCalo;
  eleTree->fEleSigmaIEtaIEta = ele->sigiEtaiEta; 
  if (!std::isnan(ele->sigiPhiiPhi)) {
    eleTree->fEleSigmaIPhiIPhi = ele->sigiPhiiPhi; 
  } else {
    eleTree->fEleSigmaIPhiIPhi = 0.0;
  }
  if (eleTree->fEleSigmaIEtaIEta*eleTree->fEleSigmaIPhiIPhi > 0) {
    eleTree->fEleSigmaIEtaIPhi = ele->CovIEtaIPhi/(eleTree->fEleSigmaIEtaIEta*eleTree->fEleSigmaIPhiIPhi);
  } else if (ele->CovIEtaIPhi>0) {
    eleTree->fEleSigmaIEtaIPhi = 1.0; 
  } else {
    eleTree->fEleSigmaIEtaIPhi = -1.0; 
  }
  eleTree->fEleSCEtaWidth = ele->SCEtaWidth;
  eleTree->fEleSCPhiWidth = ele->SCPhiWidth;

  if (ele->R9 > 5) {
    eleTree->fEleR9 = 5.0;
  } else {
    eleTree->fEleR9 = ele->R9;
  }

  eleTree->fElePreShowerOverRaw = ele->PreShowerOverRaw;
  eleTree->fEleHoverE = ele->HoverE; 

  if (ele->GsfTrackChi2OverNdof > 200) {
    eleTree->fEleGsfTrackChi2OverNdof = 200;
  } else {
    eleTree->fEleGsfTrackChi2OverNdof = ele->GsfTrackChi2OverNdof;
  }
  if (ele->KFTrackChi2OverNdof > 10) {
    eleTree->fEleKFTrackChi2OverNDoF = 10;
  } else {
    eleTree->fEleKFTrackChi2OverNDoF = ele->KFTrackChi2OverNdof;
  }
  eleTree->fEleKFTrackNHits = ele->KFTrackNHits;
  eleTree->fEleKFTrackNLayersWithMeasurement = ele->KFTrackNLayersWithMeasurement;
  if (ele->SeedE1x5OverE5x5 >= 0) {
    eleTree->fEleOneMinusSeedE1x5OverE5x5 = 1.0 - ele->SeedE1x5OverE5x5;
    if (1.0 - ele->SeedE1x5OverE5x5 > 2.0) {
      eleTree->fEleOneMinusSeedE1x5OverE5x5 = 2.0;
    }
  } else {
    eleTree->fEleOneMinusSeedE1x5OverE5x5 = -1.0;
  }
  eleTree->fElePFMVA = ele->mva;
  eleTree->fEleTrkIso03 = ele->trkIso03; 
  eleTree->fEleEMIso03 = ele->emIso03; 
  eleTree->fEleHadIso03 = ele->hadIso03; 
  eleTree->fEleTrkIso04 = ele->trkIso04; 
  eleTree->fEleEMIso04 = ele->emIso04; 
  eleTree->fEleHadIso04 = ele->hadIso04; 
  eleTree->fElePFIso04 = ComputeElePFIso04( ele, eleIndex, pfCandidates, rho, DataEra); 
  eleTree->fChargedIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.0, 0.1)/ele->pt;
  eleTree->fChargedIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.1, 0.2)/ele->pt;
  eleTree->fChargedIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.2, 0.3)/ele->pt;
  eleTree->fChargedIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.3, 0.4)/ele->pt;
  eleTree->fChargedIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFChargedIso, 0.4, 0.5)/ele->pt;
  eleTree->fGammaIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.0, 0.1)/ele->pt;
  eleTree->fGammaIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.1, 0.2)/ele->pt;
  eleTree->fGammaIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.2, 0.3)/ele->pt;
  eleTree->fGammaIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.3, 0.4)/ele->pt;
  eleTree->fGammaIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFGammaIso, 0.4, 0.5)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.0, 0.1)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.1, 0.2)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.2, 0.3)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.3, 0.4)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, eleIndex, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.4, 0.5)/ele->pt;
  eleTree->fElePassTriggerDenominator = passElectronTriggerDenominator(ele);


  //**********************************************
  //Fill Variables for Regression
  //**********************************************
  eleTree->fIsEB = ele->isEB;
  eleTree->fIsEE = Bool_t(fabs(ele->scEta) > 1.479);  //temporary fix for my bug in filling ele->isEE.
  
  eleTree->fSCRawEnergy = ele->SCRawEnergy ;
  eleTree->fNClusters = ele->nBrem+1 ;
  eleTree->fEtaSeed = ele->EtaSeed ;
  eleTree->fPhiSeed = ele->PhiSeed ;
  eleTree->fESeed = ele->ESeed ;
  eleTree->fE3x3Seed = ele->E3x3Seed ;
  eleTree->fE5x5Seed = ele->E5x5Seed ;
  eleTree->fEMaxSeed = ele->EMaxSeed ;
  eleTree->fE2ndSeed = ele->E2ndSeed ;
  eleTree->fETopSeed = ele->ETopSeed ;
  eleTree->fEBottomSeed = ele->EBottomSeed ;
  eleTree->fELeftSeed = ele->ELeftSeed ;
  eleTree->fERightSeed = ele->ERightSeed ;
  eleTree->fE2x5MaxSeed = ele->E2x5MaxSeed ;
  eleTree->fE2x5TopSeed = ele->E2x5TopSeed ;
  eleTree->fE2x5BottomSeed = ele->E2x5BottomSeed ;
  eleTree->fE2x5LeftSeed = ele->E2x5LeftSeed ;
  eleTree->fE2x5RightSeed = ele->E2x5RightSeed ;
  eleTree->fIEtaSeed = ele->IEtaSeed ;
  eleTree->fIPhiSeed = ele->IPhiSeed ;
  eleTree->fEtaCrySeed = ele->EtaCrySeed ;
  eleTree->fPhiCrySeed = ele->PhiCrySeed ;
  eleTree->fEcalEnergyError = ele->EcalEnergyError;
  eleTree->fGsfTrackPIn = ele->pIn ;
  eleTree->fCharge = ele->q ;
  eleTree->fTrackMomentumError = fmin(ele->TrackMomentumError,500.0);
  eleTree->fEleClassification = ele->Classification;

  //**********************************************
  //Find the gen level electrons
  //**********************************************
  double_t tmpGenEnergyStatus1 = -1;
  double_t tmpGenEnergyStatus3 = -1;
  double_t tmpGenEnergySupercluster = -1;

  double minDRStatus1 = 9999;
  Int_t matchedStatus1ElectronIndex = -1;
  double minDRStatus3 = 9999;  
  Int_t matchedStatus3ElectronIndex = -1;
  for(Int_t k=0; k<genParticles->GetEntries(); k++) {
    const higgsana::TGenParticle *gen = (higgsana::TGenParticle*)((*genParticles)[k]);

    //status 1 match
    if (abs(gen->pdgid) == 11 && (gen->status == 1)) {
      double tmpDR = higgsana::deltaR( ele->eta, ele->phi, gen->eta , gen->phi);
      if (tmpDR < minDRStatus1 && tmpDR < 0.3) {
        minDRStatus1 = tmpDR;
        TLorentzVector v1; 
        v1.SetPtEtaPhiM( gen->pt, gen->eta, gen->phi, ELECTRONMASS);
        tmpGenEnergyStatus1 = v1.E();
        matchedStatus1ElectronIndex = k;
//         cout << "status1 - Check energy calculation: " << gen->pt << " " << gen->eta << " : " << v1.E() << " ---> " << gen->pt / TMath::CosH(gen->eta) << endl;
      }
    }


    //status 3 match
    if (abs(gen->pdgid) == 11 && (gen->status == 3)) {
      double tmpDR = higgsana::deltaR( ele->eta, ele->phi, gen->eta , gen->phi);
      if (tmpDR < minDRStatus3 && tmpDR < 0.3) {
        minDRStatus3 = tmpDR;
        TLorentzVector v1; 
        v1.SetPtEtaPhiM( gen->pt, gen->eta, gen->phi, ELECTRONMASS);
        tmpGenEnergyStatus3 = v1.E();
        matchedStatus3ElectronIndex = k;
//         cout << "status 3 - Check energy calculation: " << gen->pt << " " << gen->eta << " : " << v1.E() << endl;
      }
    }
  }  

  

  TLorentzVector genElectronSupercluster;
  const higgsana::TGenParticle *genEleStatus1 = (higgsana::TGenParticle*)((*genParticles)[matchedStatus1ElectronIndex]);
  genElectronSupercluster.SetPtEtaPhiM( genEleStatus1->pt, genEleStatus1->eta, genEleStatus1->phi, ELECTRONMASS);

//   cout << "Status1 electron : " << genEleStatus1->pt << " " << genEleStatus1->eta << " " << genEleStatus1->phi << " , energy = " << genElectronSupercluster.E() << endl;

  for(Int_t k=0; k<genParticles->GetEntries(); k++) {
    const higgsana::TGenParticle *gen = (higgsana::TGenParticle*)((*genParticles)[k]);
    if (gen->pdgid == 22 && (gen->status == 1)) {
      //within deltaPhi road search
      if( fabs(gen->eta - genEleStatus1->eta) < 0.03 && higgsana::deltaPhi(gen->phi, genEleStatus1->phi) < 0.3) {
        TLorentzVector tmpPhoton;
        tmpPhoton.SetPtEtaPhiM( gen->pt, gen->eta, gen->phi, 0);
        TLorentzVector vtmp = genElectronSupercluster + tmpPhoton;
        genElectronSupercluster.SetPtEtaPhiM( vtmp.Pt(), vtmp.Eta(), vtmp.Phi(), vtmp.M());        
//         cout << "add photon : " << gen->pt << " " << gen->eta << " " << gen->phi << " -> " << genElectronSupercluster.Pt() << " " << genElectronSupercluster.Eta() << " " << genElectronSupercluster.Phi() << endl;
      }
    }
  }
//   cout << "Final GenLevel superclustered electron : " << genElectronSupercluster.Pt() << " " << genElectronSupercluster.Eta() << " " << genElectronSupercluster.Phi() << " , energy = " << genElectronSupercluster.E() << endl;
//   cout << endl;

     
  eleTree->fGeneratedEnergy = genElectronSupercluster.E() ;
  eleTree->fGeneratedEnergyStatus1 = tmpGenEnergyStatus1 ;
  eleTree->fGeneratedEnergyStatus3 = tmpGenEnergyStatus3 ;

  //debug
//   if (tmpGenEnergyStatus3 < 2.0) {
//     cout << "Electron : " << ele->pt << " " << ele->eta << " " << ele->phi << endl;
//     cout << "Matched Status1 : " << genEleStatus1->pt << " " << genEleStatus1->eta << " " << genEleStatus1->phi << endl;
//     for(Int_t k=0; k<genParticles->GetEntries(); k++) {
//       const higgsana::TGenParticle *gen = (higgsana::TGenParticle*)((*genParticles)[k]);
//       cout << "Particle " << k << " : " << gen->pdgid << " " << gen->status << " : " << gen->motherPdgID << " : " << gen->pt << " " << gen->eta << " " << gen->phi << endl;
//     }
//   }

  //**********************************************
  //MVA values
  //**********************************************
  eleTree->fEGammaNonTriggeringMVA = EvaluateEGammaNonTriggeringMVA(ele,eleIDMVA, false);

  //**********************************************
  //Make final consistency checks before filling
  //**********************************************
  if (TMath::IsNaN(ele->sigiPhiiPhi)) {
    cout << "Problem with sigiPhiiPhi: NaN. exit.\n";
    assert(0);
  }

  //***********************
  //Fill Electron
  //***********************
  eleTree->tree_->Fill();
}




void FillMuonTree(citana::MuonTree *muTree,
                      const higgsana::TMuon *mu, 
                      TClonesArray *pfCandidates, 
                      Double_t rho, UInt_t DataEra,
                      UInt_t NPV,
                      UInt_t runNum,
                      UInt_t lumiSec,
                      UInt_t evtNum,
                      Float_t weight = 1.0
                      
  ) {


  //Fill These Muons

  muTree->fWeight = weight;
  muTree->fRunNumber = runNum;
  muTree->fLumiSectionNumber = lumiSec;
  muTree->fEventNumber = evtNum;
  muTree->fMuEventNumberParity = (evtNum % 2 == 0);
  muTree->fRho = rho; 
  muTree->fNVertices = NPV; 

  muTree->fMuPt = mu->pt; 
  muTree->fMuEta = mu->eta; 
  muTree->fMuPhi = mu->phi; 

  muTree->fMuTypeBits = mu->typeBits;
  muTree->fIsAllArbitrated = ((mu->qualityBits & kAllArbitrated) == kAllArbitrated);
  muTree->fMuTkNchi2 = mu->tkNchi2; 
  muTree->fMuGlobalNchi2 = mu->muNchi2; 
  muTree->fMuNValidHits = mu->nValidHits; 
  muTree->fMuNTrackerHits = mu->nTkHits; 
  muTree->fMuNPixelHits = mu->nPixHits; 
  muTree->fMuNMatches = mu->nMatch ; 
  muTree->fMuD0 = mu->d0; 

  //Additional Vars 
  muTree->fMuIP3d = mu->ip3d ; 
  muTree->fMuIP3dSig = mu->ip3dSig ; 
  muTree->fMuTrkKink = mu->TrkKink ; 
  muTree->fMuGlobalKink = mu->GlobalKink ; 
  muTree->fMuSegmentCompatibility = mu->SegmentCompatibility ; 
  muTree->fMuCaloCompatibility = mu->CaloCompatilibity ; 
  muTree->fMuHadEnergy = mu->HadEnergy; 
  muTree->fMuHoEnergy = mu->HoEnergy; 
  muTree->fMuEmEnergy = mu->EmEnergy; 
  muTree->fMuHadS9Energy = mu->HadS9Energy; 
  muTree->fMuHoS9Energy = mu->HoS9Energy; 
  muTree->fMuEmS9Energy = mu->EmS9Energy; 

  //Isolation Variables
  muTree->fMuTrkIso03 = mu->trkIso03; 
  muTree->fMuEMIso03 = mu->emIso03; 
  muTree->fMuHadIso03 = mu->hadIso03; 
  muTree->fMuTrkIso05 = mu->trkIso05; 
  muTree->fMuEMIso05 = mu->emIso05; 
  muTree->fMuHadIso05 = mu->hadIso05; 
  muTree->fMuPFIso04 = ComputeMuonPFIso04( mu, pfCandidates, rho, DataEra); 

  muTree->fChargedIso_DR0p0To0p1 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFChargedIso, 0.0, 0.1);
  muTree->fChargedIso_DR0p1To0p2 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFChargedIso, 0.1, 0.2);
  muTree->fChargedIso_DR0p2To0p3 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFChargedIso, 0.2, 0.3);
  muTree->fChargedIso_DR0p3To0p4 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFChargedIso, 0.3, 0.4);
  muTree->fChargedIso_DR0p4To0p5 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFChargedIso, 0.4, 0.5);
  muTree->fGammaIso_DR0p0To0p1 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFGammaIso, 0.0, 0.1);
  muTree->fGammaIso_DR0p1To0p2 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFGammaIso, 0.1, 0.2);
  muTree->fGammaIso_DR0p2To0p3 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFGammaIso, 0.2, 0.3);
  muTree->fGammaIso_DR0p3To0p4 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFGammaIso, 0.3, 0.4);
  muTree->fGammaIso_DR0p4To0p5 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFGammaIso, 0.4, 0.5);
  muTree->fNeutralHadronIso_DR0p0To0p1 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.0, 0.1);
  muTree->fNeutralHadronIso_DR0p1To0p2 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.1, 0.2);
  muTree->fNeutralHadronIso_DR0p2To0p3 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.2, 0.3);
  muTree->fNeutralHadronIso_DR0p3To0p4 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.3, 0.4);
  muTree->fNeutralHadronIso_DR0p4To0p5 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.4, 0.5);
  muTree->fMuPassTriggerDenominator = passMuonTriggerDenominator(mu, pfCandidates, rho, DataEra);

  //***********************
  //Fill Muon
  //***********************
  muTree->tree_->Fill();



}



void FillMuonTree(citana::MuonTree *muTree,
                  const higgsana::TMuon *mu, 
                  TClonesArray *pfCandidates, 
                  vector<const higgsana::TPFCandidate*> particlesToVeto,
                  Double_t rho, UInt_t DataEra,
                  UInt_t NPV,
                  UInt_t runNum,
                  UInt_t lumiSec,
                  UInt_t evtNum,
                  Float_t weight = 1.0
                  
  ) {


  //Fill These Muons

  muTree->fWeight = weight;
  muTree->fRunNumber = runNum;
  muTree->fLumiSectionNumber = lumiSec;
  muTree->fEventNumber = evtNum;
  muTree->fMuEventNumberParity = (evtNum % 2 == 0);
  muTree->fRho = rho; 
  muTree->fNVertices = NPV; 

  muTree->fMuPt = mu->pt; 
  muTree->fMuEta = mu->eta; 
  muTree->fMuPhi = mu->phi; 

  muTree->fMuTypeBits = mu->typeBits;
  muTree->fIsAllArbitrated = ((mu->qualityBits & kAllArbitrated) == kAllArbitrated);
  muTree->fMuTkNchi2 = mu->tkNchi2; 
  muTree->fMuGlobalNchi2 = mu->muNchi2; 
  muTree->fMuNValidHits = mu->nValidHits; 
  muTree->fMuNTrackerHits = mu->nTkHits; 
  muTree->fMuNPixelHits = mu->nPixHits; 
  muTree->fMuNMatches = mu->nMatch ; 
  muTree->fMuD0 = mu->d0; 

  //Additional Vars 
  muTree->fMuIP3d = mu->ip3d ; 
  muTree->fMuIP3dSig = mu->ip3dSig ; 
  muTree->fMuTrkKink = mu->TrkKink ; 
  muTree->fMuGlobalKink = mu->GlobalKink ; 
  muTree->fMuSegmentCompatibility = mu->SegmentCompatibility ; 
  muTree->fMuCaloCompatibility = mu->CaloCompatilibity ; 
  muTree->fMuHadEnergy = mu->HadEnergy; 
  muTree->fMuHoEnergy = mu->HoEnergy; 
  muTree->fMuEmEnergy = mu->EmEnergy; 
  muTree->fMuHadS9Energy = mu->HadS9Energy; 
  muTree->fMuHoS9Energy = mu->HoS9Energy; 
  muTree->fMuEmS9Energy = mu->EmS9Energy; 

  //Isolation Variables
  muTree->fMuTrkIso03 = mu->trkIso03; 
  muTree->fMuEMIso03 = mu->emIso03; 
  muTree->fMuHadIso03 = mu->hadIso03; 
  muTree->fMuTrkIso05 = mu->trkIso05; 
  muTree->fMuEMIso05 = mu->emIso05; 
  muTree->fMuHadIso05 = mu->hadIso05; 
  muTree->fMuPFIso04 = ComputeMuonPFIso04( mu, pfCandidates, rho, DataEra); 

  muTree->fChargedIso_DR0p0To0p1 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.0, 0.1);
  muTree->fChargedIso_DR0p1To0p2 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.1, 0.2);
  muTree->fChargedIso_DR0p2To0p3 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.2, 0.3);
  muTree->fChargedIso_DR0p3To0p4 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.3, 0.4);
  muTree->fChargedIso_DR0p4To0p5 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.4, 0.5);
  muTree->fGammaIso_DR0p0To0p1 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.0, 0.1);
  muTree->fGammaIso_DR0p1To0p2 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.1, 0.2);
  muTree->fGammaIso_DR0p2To0p3 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.2, 0.3);
  muTree->fGammaIso_DR0p3To0p4 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.3, 0.4);
  muTree->fGammaIso_DR0p4To0p5 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.4, 0.5);
  muTree->fNeutralHadronIso_DR0p0To0p1 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.0, 0.1);
  muTree->fNeutralHadronIso_DR0p1To0p2 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.1, 0.2);
  muTree->fNeutralHadronIso_DR0p2To0p3 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.2, 0.3);
  muTree->fNeutralHadronIso_DR0p3To0p4 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.3, 0.4);
  muTree->fNeutralHadronIso_DR0p4To0p5 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.4, 0.5);
  muTree->fMuPassTriggerDenominator = passMuonTriggerDenominator(mu, pfCandidates, rho, DataEra);

  //***********************
  //Fill Muon
  //***********************
  muTree->tree_->Fill();



}



#endif
