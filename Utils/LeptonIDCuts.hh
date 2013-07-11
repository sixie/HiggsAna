#ifndef LEPTONIDCUTS_HH
#define LEPTONIDCUTS_HH

#include "HiggsAna/Utils/CommonDefs.hh"
#include "HiggsAna/Utils/CommonTools.hh"
#include "HiggsAna/Utils/IsolationPileupCorrections.hh"
#include "HiggsAna/DataTree/interface/TPhoton.hh"
#include "HiggsAna/DataTree/interface/TElectron.hh"
#include "HiggsAna/DataTree/interface/TMuon.hh"
#include "HiggsAna/DataTree/interface/HiggsAnaDefs.hh"
#include <cassert>
#include "TMath.h"

Double_t ComputeMuonPFIso04( const higgsana::TMuon *muon,
                             TClonesArray *pfCandidates, 
                             Double_t rho, UInt_t DataEra,
                             Bool_t printDebug = kFALSE);
Double_t ComputeMuonPFIso04( const higgsana::TMuon *muon,
                             TClonesArray *pfCandidates,
                             vector<const higgsana::TPFCandidate*> particlesToVeto,
                             Double_t rho, UInt_t DataEra,
                             Bool_t printDebug = kFALSE);

Double_t ComputeElePFIso04( const higgsana::TElectron *ele, 
                            TClonesArray *pfCandidates, 
                            Double_t rho, UInt_t DataEra, 
                            Bool_t printDebug = kFALSE);
Double_t ComputeElePFIso04( const higgsana::TElectron *ele, 
                            Int_t eleIndex,
                            TClonesArray *pfCandidates, 
                            Double_t rho, UInt_t DataEra, 
                            Bool_t printDebug = kFALSE);
Double_t ComputeElePFIso04( const higgsana::TElectron *ele, 
                            Int_t eleIndex,
                            TClonesArray *pfCandidates, 
                            vector<const higgsana::TPFCandidate*> particlesToVeto,
                            Double_t rho, UInt_t DataEra, 
                            Bool_t printDebug = kFALSE);
Double_t ComputePhotonPFIso03( const higgsana::TPhoton *pho,
                             TClonesArray *pfCandidates, 
                             Double_t rho, UInt_t DataEra,
                             Bool_t printDebug = kFALSE);

Double_t ComputeMuonPFIsoRings( const higgsana::TMuon *muon,
                                TClonesArray *pfCandidates, 
                                Double_t rho, UInt_t DataEra,
                                UInt_t pfType,
                                Double_t minDR,
                                Double_t maxDR,
                                Bool_t printDebug = kFALSE);
Double_t ComputeMuonPFIsoRings( const higgsana::TMuon *muon,
                                TClonesArray *pfCandidates, 
                                vector<const higgsana::TPFCandidate*> particlesToVeto,
                                Double_t rho, UInt_t DataEra,
                                UInt_t pfType,
                                Double_t minDR,
                                Double_t maxDR,
                                Bool_t printDebug = kFALSE);

Double_t ComputeElePFIsoRings( const higgsana::TElectron *ele, 
                               Int_t eleIndex,
                               TClonesArray *pfCandidates, 
                               vector<const higgsana::TPFCandidate*> particlesToVeto,
                               Double_t rho, UInt_t DataEra,
                               UInt_t pfIsoType,
                               Double_t minDR,
                               Double_t maxDR,
                               Bool_t printDebug = kFALSE);
Double_t ComputeElePFIsoRings( const higgsana::TElectron *ele, 
                               Int_t eleIndex,
                               TClonesArray *pfCandidates, 
                               Double_t rho, UInt_t DataEra,
                               UInt_t pfIsoType,
                               Double_t minDR,
                               Double_t maxDR,
                               Bool_t printDebug = kFALSE);

Double_t ComputePhotonNonCorrectedIsoDr03(const higgsana::TPFCandidate *photon, 
                             UInt_t leptonType,
                             UInt_t leptonIndex,
                             TClonesArray *fPFCandidates,
                             Bool_t printDebug = kFALSE);

Double_t ComputePhotonDBetaCorrectedRelativeIsoDr03(const higgsana::TPFCandidate *photon, 
                                       UInt_t leptonType,
                                       UInt_t leptonIndex,
                                       TClonesArray *fPFCandidates,
                                       Bool_t printDebug = kFALSE);


Bool_t passConversionVeto(Int_t isConv);
Bool_t passMuonID(const higgsana::TMuon *muon, Double_t pfIso03);
Bool_t passCutBasedEleID(const higgsana::TElectron *electron, Double_t pfIso04);
Bool_t passCutBasedTightEleID(const higgsana::TElectron *electron);
Bool_t passCutBasedTightEleIDNoIso(const higgsana::TElectron *electron);
Bool_t passElectronTriggerDenominator(const higgsana::TElectron *electron);


Bool_t PassEleSimpleCutsLoose( const higgsana::TElectron *ele,   
                               TClonesArray *pfCandidates, 
                               Double_t rho, UInt_t DataEra, 
                               Bool_t printDebug = kFALSE);  
Bool_t PassEleSimpleCutsMedium( const higgsana::TElectron *ele,   
                                TClonesArray *pfCandidates, 
                                Double_t rho, UInt_t DataEra, 
                                Bool_t printDebug = kFALSE);  
Bool_t PassEleSimpleCutsTight( const higgsana::TElectron *ele,   
                                TClonesArray *pfCandidates, 
                                Double_t rho, UInt_t DataEra, 
                                Bool_t printDebug = kFALSE);  
Bool_t PassEleSimpleCutsVeryTight( const higgsana::TElectron *ele,   
                                   TClonesArray *pfCandidates, 
                                   Double_t rho, UInt_t DataEra, 
                                   Bool_t printDebug = kFALSE);  

Bool_t passPhotonSimpleCuts(const higgsana::TPhoton *pho);


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
Double_t ComputeMuonPFIso04( const higgsana::TMuon *muon,
                           TClonesArray *pfCandidates, 
                           Double_t rho, UInt_t DataEra,
                           Bool_t printDebug) {

  Double_t fChargedIso  = 0.0;
  Double_t fGammaIso  = 0.0;
  Double_t fNeutralHadronIso  = 0.0;

  //
  //Loop over PF Candidates
  //
  for(UInt_t k=0; int(k) < pfCandidates->GetEntries(); ++k) {

    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*pfCandidates)[k]);

    //*****************************
    // veto Pileup PF Candidates
    //*****************************
    if (!pf->IsPFNoPU) continue;
    
    Double_t dr = higgsana::deltaR(muon->eta, muon->phi, pf->eta, pf->phi);
    if (dr > 0.4) continue;

    //***********************************************************************
    // Charged Iso 
    //***********************************************************************
    if (pf->q != 0) {      
      //don't include pf electrons or pf muons
      if (pf->pfType == eElectron || pf->pfType == eMuon) continue;
      if (dr > 0.0001) {
        if( printDebug ) cout << "charged:: " << pf->pt << " : " << pf->dz << " : "
                              << dr << endl;
        fChargedIso += pf->pt;
      }
    }
    
    //***********************************************************************
    // Gamma Iso 
    //***********************************************************************
    else if (pf->pfType == eGamma) {
      if( pf->pt > 0.5 && dr > 0.01) {
        if( printDebug ) cout << "gamma:: " << pf->pt << " " 
                              << dr << endl;
        
        fGammaIso += pf->pt;
      }
    }
    
    //***********************************************************************
    // Other Neutrals
    //***********************************************************************
    else {
      if( pf->pt > 0.5 && dr > 0.01) {
        if( printDebug ) cout << "neutral hadron:: " << pf->pt << " " 
                              << dr << endl;        
	fNeutralHadronIso += pf->pt;
      }
    }
  }

  if(printDebug) cout << "rho: " << rho << endl;

  double pfIso = fChargedIso + fmax(0.0,(fGammaIso + fNeutralHadronIso 
                                         -rho*MuonEffectiveArea(kMuGammaAndNeutralHadronIso04,
                                                                muon->eta,DataEra)));
  
  if( printDebug ) { 
    cout << "PFiso: " << pfIso
	 << "\tfChargedIso: " << fChargedIso
	 << "\tfGammaIso: " << fGammaIso
	 << "\tfNeutralHadronIso: " << fNeutralHadronIso
	 << endl;
  }

  return pfIso;

}


//--------------------------------------------------------------------------------------------------
Double_t ComputeMuonPFIso04( const higgsana::TMuon *muon,
                             TClonesArray *pfCandidates, 
                             vector<const higgsana::TPFCandidate*> particlesToVeto,
                             Double_t rho, UInt_t DataEra,
                             Bool_t printDebug) {

  Double_t fChargedIso  = 0.0;
  Double_t fGammaIso  = 0.0;
  Double_t fNeutralHadronIso  = 0.0;

  //
  //Loop over PF Candidates
  //
  for(UInt_t k=0; int(k) < pfCandidates->GetEntries(); ++k) {

    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*pfCandidates)[k]);

    //*****************************
    // veto Pileup PF Candidates
    //*****************************
    if (!pf->IsPFNoPU) continue;
    
    //*******************************************
    // don't include particles on veto list
    //*******************************************
    bool isVeto = false;
    for( UInt_t p=0; p<particlesToVeto.size(); p++ ) {
      if( pf == particlesToVeto[p] ) { 
	isVeto = true;
	break;
      }
    }
    if( isVeto ) continue;

    Double_t dr = higgsana::deltaR(muon->eta, muon->phi, pf->eta, pf->phi);
    if (dr > 0.4) continue;

    //***********************************************************************
    // Charged Iso 
    //***********************************************************************
    if (pf->q != 0) {      
      //don't include pf electrons or pf muons
      if (pf->pfType == eElectron || pf->pfType == eMuon) continue;
      if( printDebug ) cout << "charged:: " << pf->pt << " " 
                           << dr << endl;
      fChargedIso += pf->pt;
    }
    
    //***********************************************************************
    // Gamma Iso  
    //***********************************************************************
    else if (pf->pfType == eGamma) {
      if( pf->pt > 0.5 ) {
      if( printDebug ) cout << "gamma:: " << pf->pt << " " 
                           << dr << endl;
        fGammaIso += pf->pt;
      }
    }
    
    //***********************************************************************
    // Other Neutrals
    //***********************************************************************
    else {
      if( pf->pt > 0.5 ) {
        if( printDebug ) cout << "neutral:: " << pf->pt << " " 
                              << dr << endl;
        
	fNeutralHadronIso += pf->pt;
      }
    }
  }

  if(printDebug) cout << "rho: " << rho << endl;

  double pfIso = fChargedIso + fmax(0.0,(fGammaIso + fNeutralHadronIso 
                                         -rho*MuonEffectiveArea(kMuGammaAndNeutralHadronIso04,
                                                                muon->eta,DataEra)));
  
  if( printDebug ) { 
    cout << "PFiso: " << pfIso
	 << "\tfChargedIso: " << fChargedIso
	 << "\tfGammaIso: " << fGammaIso
	 << "\tfNeutralHadronIso: " << fNeutralHadronIso
	 << endl;
  }

  return pfIso;

}

//--------------------------------------------------------------------------------------------------
Double_t ComputeElePFIso04( const higgsana::TElectron *ele, 
                            TClonesArray *pfCandidates, 
                            Double_t rho, UInt_t DataEra, 
                            Bool_t printDebug) {

  return ComputeElePFIso04(ele, -1, pfCandidates, rho, DataEra, printDebug);
}


//--------------------------------------------------------------------------------------------------
Double_t ComputeElePFIso04( const higgsana::TElectron *ele, 
                            Int_t eleIndex,
                            TClonesArray *pfCandidates, 
                            Double_t rho, UInt_t DataEra, 
                            Bool_t printDebug) {
 
  
  Double_t fChargedIso  = 0.0;
  Double_t fGammaIso  = 0.0;
  Double_t fNeutralHadronIso  = 0.0;
  
  //
  //Loop over PF Candidates
  //
  for(UInt_t k=0; int(k) < pfCandidates->GetEntries(); ++k) {

    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*pfCandidates)[k]);

    //*****************************
    // veto Pileup PF Candidates
    //*****************************
    if (!pf->IsPFNoPU) continue;

    Double_t dr = higgsana::deltaR(ele->eta, ele->phi, pf->eta, pf->phi);
    if (dr > 0.4) continue;

    if(printDebug) { 
      cout << "pf :: type: " << pf->pfType << "\tpt: " << pf->pt << "\tdR: " << dr;
      cout << "\tdZ: " << pf->dz;
      cout << endl;
    }

    //***********************************************************************
    // Charged Iso 
    //***********************************************************************
    if (pf->q != 0) {

      // Veto any PFmuon, or PFEle
      if ( pf->pfType == eElectron || pf->pfType == eMuon) {
        continue;
      }

      // Footprint Veto
      if (fabs(ele->scEta) > 1.479 && dr < 0.015) continue;

      if( printDebug) cout << "charged:: pt: " << pf->pt 
                           << "\ttype: " << pf->pfType 
                           << endl;

      fChargedIso += pf->pt;
    }

    //***********************************************************************
    // Gamma Iso 
    //***********************************************************************
    else if (pf->pfType == eGamma) {
      if( printDebug) cout << "gamma:: " << pf->pt << " : " << pf->matchedObjectType << " " <<  pf->matchedObjectIndex << " ? " << eleIndex << " : "
                           << dr << endl;
      if (eleIndex >= 0) {
        //proper veto of electron footprint
        if (pf->matchedObjectType == 11 && Int_t(pf->matchedObjectIndex) == eleIndex) continue;
        //this is a temporary fix
        if (fabs(ele->scEta) < 1.479) {
          if (fabs(ele->eta - pf->eta) < 0.015) continue;
        }
      }

      if (fabs(ele->scEta) > 1.479) {
        if (higgsana::deltaR(ele->eta,ele->phi, pf->eta, pf->phi) < 0.08) continue;
      }
      fGammaIso += pf->pt;
    }

    //***********************************************************************
    // Other Neutrals
    //***********************************************************************
    else {
      if( printDebug ) cout << "neutral:: " << pf->pt << " " 
                           << dr << endl;
      fNeutralHadronIso += pf->pt;
    }
  }

  if( printDebug ) cout << "rho: " << rho << endl;

  double pfIso = fChargedIso + fmax(0.0,(fGammaIso + fNeutralHadronIso 
					-rho*ElectronEffectiveArea(kEleGammaAndNeutralHadronIso04,
								   ele->eta,DataEra)));


  if( printDebug ) { 
    cout << "PFiso: " << pfIso
	 << "\tfChargedIso: " << fChargedIso
	 << "\tfGammaIso: " << fGammaIso
	 << "\tfNeutralHadronIso: " << fNeutralHadronIso
	 << endl;
  }

  return pfIso;

}


//--------------------------------------------------------------------------------------------------
Double_t ComputeElePFIso04( const higgsana::TElectron *ele, 
                            Int_t eleIndex,
                            TClonesArray *pfCandidates, 
                            vector<const higgsana::TPFCandidate*> particlesToVeto,
                            Double_t rho, UInt_t DataEra, 
                            Bool_t printDebug) {
 
  
  Double_t fChargedIso  = 0.0;
  Double_t fGammaIso  = 0.0;
  Double_t fNeutralHadronIso  = 0.0;
  
  //
  //Loop over PF Candidates
  //
  for(UInt_t k=0; int(k) < pfCandidates->GetEntries(); ++k) {

    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*pfCandidates)[k]);

    //*****************************
    // veto Pileup PF Candidates
    //*****************************
    if (!pf->IsPFNoPU) continue;

    //*******************************************
    // don't include particles on veto list
    //*******************************************
    bool isVeto = false;
    for( UInt_t p=0; p<particlesToVeto.size(); p++ ) {
      if( pf == particlesToVeto[p] ) { 
	isVeto = true;
	break;
      }
    }
    if( isVeto ) continue;

    Double_t dr = higgsana::deltaR(ele->eta, ele->phi, pf->eta, pf->phi);
    if (dr > 0.4) continue;

    if(printDebug) { 
      cout << "pf :: type: " << pf->pfType << "\tpt: " << pf->pt << "\tdR: " << dr;
      cout << "\tdZ: " << pf->dz;
      cout << endl;
    }

    //***********************************************************************
    // Charged Iso 
    //***********************************************************************
    if (pf->q != 0) {

      // Veto any PFmuon, or PFEle
      if ( pf->pfType == eElectron || pf->pfType == eMuon) {
        continue;
      }

      // Footprint Veto
      if (fabs(ele->scEta) > 1.479 && dr < 0.015) continue;

      if( printDebug) cout << "charged:: pt: " << pf->pt 
                           << "\ttype: " << pf->pfType
                           << endl;

      fChargedIso += pf->pt;
    }

    //***********************************************************************
    // Gamma Iso 
    //***********************************************************************
    else if (pf->pfType == eGamma) {

      if (eleIndex >= 0) {
        //proper veto of electron footprint
        if (pf->matchedObjectType == 11 && Int_t(pf->matchedObjectIndex) == eleIndex) continue;
        //this is a temporary fix
        if (fabs(ele->scEta) < 1.479) {
          if (fabs(ele->eta - pf->eta) < 0.015) continue;
        }
      }

      if (fabs(ele->scEta) > 1.479) {
        if (higgsana::deltaR(ele->eta,ele->phi, pf->eta, pf->phi) < 0.08) continue;
      }
      if( printDebug) cout << "gamma:: " << pf->pt << " " 
                           << dr << endl;
      fGammaIso += pf->pt;
    }

    //***********************************************************************
    // Other Neutrals
    //***********************************************************************
    else {
      if( printDebug ) cout << "neutral:: " << pf->pt << " " 
                           << dr << endl;
      fNeutralHadronIso += pf->pt;
    }
  }

  if( printDebug ) cout << "rho: " << rho << endl;

  double pfIso = fChargedIso + fmax(0.0,(fGammaIso + fNeutralHadronIso 
					-rho*ElectronEffectiveArea(kEleGammaAndNeutralHadronIso04,
								   ele->eta,DataEra)));


  if( printDebug ) { 
    cout << "PFiso: " << pfIso
	 << "\tfChargedIso: " << fChargedIso
	 << "\tfGammaIso: " << fGammaIso
	 << "\tfNeutralHadronIso: " << fNeutralHadronIso
	 << endl;
  }

  return pfIso;

}


//--------------------------------------------------------------------------------------------------
Double_t ComputePhotonPFIso03( const higgsana::TPhoton *pho,
                               TClonesArray *pfCandidates, 
                               UInt_t pfIsoType,
                               Double_t rho, UInt_t DataEra,
                               Bool_t printDebug) {

  Double_t fIso = 0.0;

  //
  //Loop over PF Candidates
  //
  for(UInt_t k=0; int(k) < pfCandidates->GetEntries(); ++k) {

    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*pfCandidates)[k]);

    //*****************************
    // veto Pileup PF Candidates
    //*****************************
    if (!pf->IsPFNoPU) continue;
    
    Double_t dr = higgsana::deltaR(pho->eta, pho->phi, pf->eta, pf->phi);
    if (dr > 0.3) continue;

    //***********************************************************************
    // Charged Iso 
    //***********************************************************************
    if (pf->q != 0) {   
      //don't include pf electrons or pf muons
      if (pf->pfType == eElectron || pf->pfType == eMuon) continue;
      if (fabs(pf->dz) > 0.2) continue;
      if (dr <= 0.02) continue;
      if (pfIsoType == kPFChargedIso) {
        if( printDebug) cout << "charged:: " << pf->pt << " : "
                             << pf->dz << " : " 
                             << dr << endl;
        fIso += pf->pt;
      }
    }
    
    //***********************************************************************
    // Gamma Iso 
    //***********************************************************************
    else if (pf->pfType == eGamma) {
      if ( fabs(pho->scEta) < 1.479) {
        if ( fabs(pho->eta - pf->eta) <= 0.015) continue;
        if ( fabs(pho->scEta - pf->eta) <= 0.015) continue;
      } else {
        if ( dr <= 0.00864*fabs(sinh(pho->scEta))*4) continue;
      }
      if (pfIsoType == kPFGammaIso) {
        if( printDebug) cout << "gamma:: " << pf->pt << " : " << pho->eta << " " << pho->scEta << " " << pf->eta  << " : "
                             << dr << endl;
        fIso += pf->pt;
      }
    }
    
    //***********************************************************************
    // Other Neutrals
    //***********************************************************************
    else {
      if (pfIsoType == kPFNeutralHadronIso) {
      if( printDebug) cout << "neutral hadron:: " << pf->pt << " " 
                           << dr << endl;
        fIso += pf->pt;
      }
    }
  }
  
  if(printDebug) cout << "rho: " << rho << endl;
  

  UInt_t EAType = 0;
  if (pfIsoType == kPFChargedIso) {
    EAType = kPhotonChargedIso03;
  } else if (pfIsoType == kPFGammaIso) {
    EAType = kPhotonGammaIso03;
  } else if (pfIsoType == kPFNeutralHadronIso) {
    EAType = kPhotonNeutralHadronIso03;
  } else {
    EAType = kPhotonNoCorrection;
  }

  double iso = fmax(0.0, fIso - rho*PhotonEffectiveArea(EAType, pho->scEta,DataEra));
  
  if( printDebug ) { 
    cout << "uncorrected iso: " << fIso
	 << "\tcorrected iso: " << iso
	 << endl;
  }

  return iso;


}



//--------------------------------------------------------------------------------------------------
Double_t ComputeMuonPFIsoRings( const higgsana::TMuon *muon,
                                TClonesArray *pfCandidates, 
                                Double_t rho, UInt_t DataEra,
                                UInt_t pfIsoType,
                                Double_t minDR,
                                Double_t maxDR,
                                Bool_t printDebug) {

  Double_t fIso = 0.0;

  //
  //Loop over PF Candidates
  //
  for(UInt_t k=0; int(k) < pfCandidates->GetEntries(); ++k) {

    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*pfCandidates)[k]);

    //*****************************
    // veto Pileup PF Candidates
    //*****************************
    if (!pf->IsPFNoPU) continue;
    
    Double_t dr = higgsana::deltaR(muon->eta, muon->phi, pf->eta, pf->phi);
    if (dr > 0.5) continue;
    if (dr < minDR) continue;
    if (dr >= maxDR) continue;

    //***********************************************************************
    // Charged Iso 
    //***********************************************************************
    if (pf->q != 0) {      
      //don't include pf electrons or pf muons
      if (pf->pfType == eElectron || pf->pfType == eMuon) continue;

      if (pfIsoType == kPFChargedIso) fIso += pf->pt;
    }
    
    //***********************************************************************
    // Gamma Iso 
    //***********************************************************************
    else if (pf->pfType == eGamma) {
      if( pf->pt > 0.5 ) {
        if (pfIsoType == kPFGammaIso) fIso += pf->pt;
      }
    }
    
    //***********************************************************************
    // Other Neutrals
    //***********************************************************************
    else {
      if( pf->pt > 0.5 ) {
	if (pfIsoType == kPFNeutralHadronIso) fIso += pf->pt;
      }
    }
  }

  if(printDebug) cout << "rho: " << rho << endl;

  UInt_t EAType = 0;
  if (pfIsoType == kPFGammaIso) {
    if (minDR == 0.0 && maxDR == 0.1) EAType = kMuGammaIsoDR0p0To0p1;
    else if (minDR == 0.1 && maxDR == 0.2) EAType = kMuGammaIsoDR0p1To0p2;
    else if (minDR == 0.2 && maxDR == 0.3) EAType = kMuGammaIsoDR0p2To0p3;
    else if (minDR == 0.3 && maxDR == 0.4) EAType = kMuGammaIsoDR0p3To0p4;
    else if (minDR == 0.4 && maxDR == 0.5) EAType = kMuGammaIsoDR0p4To0p5;
    else EAType = kMuNoCorrection;
  } else if (pfIsoType == kPFNeutralHadronIso) {
    if (minDR == 0.0 && maxDR == 0.1) EAType = kMuNeutralHadronIsoDR0p0To0p1;
    else if (minDR == 0.1 && maxDR == 0.2) EAType = kMuNeutralHadronIsoDR0p1To0p2;
    else if (minDR == 0.2 && maxDR == 0.3) EAType = kMuNeutralHadronIsoDR0p2To0p3;
    else if (minDR == 0.3 && maxDR == 0.4) EAType = kMuNeutralHadronIsoDR0p3To0p4;
    else if (minDR == 0.4 && maxDR == 0.5) EAType = kMuNeutralHadronIsoDR0p4To0p5;
    else EAType = kMuNoCorrection;
  } else {
    EAType = kMuNoCorrection;
  }

  double iso = fmax(0.0, fIso - rho*MuonEffectiveArea(EAType, muon->eta,DataEra));
  
  if( printDebug ) { 
    cout << "uncorrected iso: " << fIso
	 << "\tcorrected iso: " << iso
	 << endl;
  }

  return iso;

}


//--------------------------------------------------------------------------------------------------
Double_t ComputeMuonPFIsoRings( const higgsana::TMuon *muon,
                                TClonesArray *pfCandidates, 
                                vector<const higgsana::TPFCandidate*> particlesToVeto,
                                Double_t rho, UInt_t DataEra,
                                UInt_t pfIsoType,
                                Double_t minDR,
                                Double_t maxDR,
                                Bool_t printDebug) {

  Double_t fIso = 0.0;

  //
  //Loop over PF Candidates
  //
  for(UInt_t k=0; int(k) < pfCandidates->GetEntries(); ++k) {

    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*pfCandidates)[k]);

    //*****************************
    // veto Pileup PF Candidates
    //*****************************
    if (!pf->IsPFNoPU) continue;
    
    //*******************************************
    // don't include particles on veto list
    //*******************************************
    bool isVeto = false;
    for( UInt_t p=0; p<particlesToVeto.size(); p++ ) {
      if( pf == particlesToVeto[p] ) { 
	isVeto = true;
	break;
      }
    }
    if( isVeto ) continue;

    Double_t dr = higgsana::deltaR(muon->eta, muon->phi, pf->eta, pf->phi);
    if (dr > 0.5) continue;
    if (dr < minDR) continue;
    if (dr >= maxDR) continue;

    //***********************************************************************
    // Charged Iso 
    //***********************************************************************
    if (pf->q != 0) {      
      //don't include pf electrons or pf muons
      if (pf->pfType == eElectron || pf->pfType == eMuon) continue;

      if (pfIsoType == kPFChargedIso) fIso += pf->pt;
    }
    
    //***********************************************************************
    // Gamma Iso 
    //***********************************************************************
    else if (pf->pfType == eGamma) {
      if( pf->pt > 0.5 ) {
        if (pfIsoType == kPFGammaIso) fIso += pf->pt;
      }
    }
    
    //***********************************************************************
    // Other Neutrals
    //***********************************************************************
    else {
      if( pf->pt > 0.5 ) {
	if (pfIsoType == kPFNeutralHadronIso) fIso += pf->pt;
      }
    }
  }

  if(printDebug) cout << "rho: " << rho << endl;

  UInt_t EAType = 0;
  if (pfIsoType == kPFGammaIso) {
    if (minDR == 0.0 && maxDR == 0.1) EAType = kMuGammaIsoDR0p0To0p1;
    else if (minDR == 0.1 && maxDR == 0.2) EAType = kMuGammaIsoDR0p1To0p2;
    else if (minDR == 0.2 && maxDR == 0.3) EAType = kMuGammaIsoDR0p2To0p3;
    else if (minDR == 0.3 && maxDR == 0.4) EAType = kMuGammaIsoDR0p3To0p4;
    else if (minDR == 0.4 && maxDR == 0.5) EAType = kMuGammaIsoDR0p4To0p5;
    else EAType = kMuNoCorrection;
  } else if (pfIsoType == kPFNeutralHadronIso) {
    if (minDR == 0.0 && maxDR == 0.1) EAType = kMuNeutralHadronIsoDR0p0To0p1;
    else if (minDR == 0.1 && maxDR == 0.2) EAType = kMuNeutralHadronIsoDR0p1To0p2;
    else if (minDR == 0.2 && maxDR == 0.3) EAType = kMuNeutralHadronIsoDR0p2To0p3;
    else if (minDR == 0.3 && maxDR == 0.4) EAType = kMuNeutralHadronIsoDR0p3To0p4;
    else if (minDR == 0.4 && maxDR == 0.5) EAType = kMuNeutralHadronIsoDR0p4To0p5;
    else EAType = kMuNoCorrection;
  } else {
    EAType = kMuNoCorrection;
  }

  double iso = fmax(0.0, fIso - rho*MuonEffectiveArea(EAType, muon->eta,DataEra));
  
  if( printDebug ) { 
    cout << "uncorrected iso: " << fIso
	 << "\tcorrected iso: " << iso
	 << endl;
  }

  return iso;

}



//--------------------------------------------------------------------------------------------------
Double_t ComputeElePFIsoRings( const higgsana::TElectron *ele, 
                               Int_t eleIndex,
                               TClonesArray *pfCandidates, 
                               Double_t rho, UInt_t DataEra,
                               UInt_t pfIsoType,
                               Double_t minDR,
                               Double_t maxDR,
                               Bool_t printDebug) {
 
  Double_t fIso = 0.0;
  
  //
  //Loop over PF Candidates
  //
  for(UInt_t k=0; int(k) < pfCandidates->GetEntries(); ++k) {

    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*pfCandidates)[k]);

    //*****************************
    // veto Pileup PF Candidates
    //*****************************
    if (!pf->IsPFNoPU) continue;

    Double_t dr = higgsana::deltaR(ele->eta, ele->phi, pf->eta, pf->phi);
    if (dr > 0.5) continue;
    if (dr < minDR) continue;
    if (dr >= maxDR) continue;

    if(printDebug) { 
      cout << "pf :: type: " << pf->pfType << "\tpt: " << pf->pt << "\tdR: " << dr;
      cout << "\tdZ: " << pf->dz;
      cout << endl;
    }

    //***********************************************************************
    // Charged Iso 
    //***********************************************************************
    if (pf->q != 0) {

      // Veto any PFmuon, or PFEle
      if ( pf->pfType == eElectron || pf->pfType == eMuon) {
        continue;
      }

      // Footprint Veto
      if (fabs(ele->scEta) > 1.479 && dr < 0.015) continue;

      if( printDebug) cout << "charged:: pt: " << pf->pt 
                           << "\ttype: " << pf->pfType
                           << endl;

      if (pfIsoType == kPFChargedIso) fIso += pf->pt;
    }

    //***********************************************************************
    // Gamma Iso 
    //***********************************************************************
    else if (pf->pfType == eGamma) {

      if (eleIndex >= 0) {
        //proper veto of electron footprint
        if (pf->matchedObjectType == 11 && Int_t(pf->matchedObjectIndex) == eleIndex) continue;
        //this is a temporary fix
        if (fabs(ele->scEta) < 1.479) {
          if (fabs(ele->eta - pf->eta) < 0.015) continue;
        }
      }

      if (fabs(ele->scEta) > 1.479) {
        if (higgsana::deltaR(ele->eta,ele->phi, pf->eta, pf->phi) < 0.08) continue;
      }
      if( printDebug) cout << "gamma:: " << pf->pt << " " 
                           << dr << endl;
      if (pfIsoType == kPFGammaIso) fIso += pf->pt;
    }

    //***********************************************************************
    // Other Neutrals
    //***********************************************************************
    else {
      if( printDebug ) cout << "neutral:: " << pf->pt << " " 
                           << dr << endl;
      if (pfIsoType == kPFNeutralHadronIso) fIso += pf->pt;
    }
  }

  if( printDebug ) cout << "rho: " << rho << endl;

  UInt_t EAType = 0;
  if (pfIsoType == kPFGammaIso) {
    if (minDR == 0.0 && maxDR == 0.1) EAType = kEleGammaIsoDR0p0To0p1;
    else if (minDR == 0.1 && maxDR == 0.2) EAType = kEleGammaIsoDR0p1To0p2;
    else if (minDR == 0.2 && maxDR == 0.3) EAType = kEleGammaIsoDR0p2To0p3;
    else if (minDR == 0.3 && maxDR == 0.4) EAType = kEleGammaIsoDR0p3To0p4;
    else if (minDR == 0.4 && maxDR == 0.5) EAType = kEleGammaIsoDR0p4To0p5;
    else EAType = kEleNoCorrection;
  } else if (pfIsoType == kPFNeutralHadronIso) {
    if (minDR == 0.0 && maxDR == 0.1) EAType = kEleNeutralHadronIsoDR0p0To0p1;
    else if (minDR == 0.1 && maxDR == 0.2) EAType = kEleNeutralHadronIsoDR0p1To0p2;
    else if (minDR == 0.2 && maxDR == 0.3) EAType = kEleNeutralHadronIsoDR0p2To0p3;
    else if (minDR == 0.3 && maxDR == 0.4) EAType = kEleNeutralHadronIsoDR0p3To0p4;
    else if (minDR == 0.4 && maxDR == 0.5) EAType = kEleNeutralHadronIsoDR0p4To0p5;
    else EAType = kEleNoCorrection;
  } else {
    EAType = kEleNoCorrection;
  }

  double iso = fmax(0.0, fIso - rho*ElectronEffectiveArea(EAType, ele->scEta,DataEra));


  if( printDebug ) { 
    cout << "uncorrected iso: " << fIso
	 << "\tcorrected iso: " << iso
	 << endl;
  }

  return iso;

}


//--------------------------------------------------------------------------------------------------
Double_t ComputeElePFIsoRings( const higgsana::TElectron *ele, 
                               Int_t eleIndex,
                               TClonesArray *pfCandidates, 
                               vector<const higgsana::TPFCandidate*> particlesToVeto,
                               Double_t rho, UInt_t DataEra,
                               UInt_t pfIsoType,
                               Double_t minDR,
                               Double_t maxDR,
                               Bool_t printDebug) {
 
  Double_t fIso = 0.0;
  
  //
  //Loop over PF Candidates
  //
  for(UInt_t k=0; int(k) < pfCandidates->GetEntries(); ++k) {

    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*pfCandidates)[k]);

    //*****************************
    // veto Pileup PF Candidates
    //*****************************
    if (!pf->IsPFNoPU) continue;

    //*******************************************
    // don't include particles on veto list
    //*******************************************
    bool isVeto = false;
    for( UInt_t p=0; p<particlesToVeto.size(); p++ ) {
      if( pf == particlesToVeto[p] ) { 
	isVeto = true;
	break;
      }
    }
    if( isVeto ) continue;

    Double_t dr = higgsana::deltaR(ele->eta, ele->phi, pf->eta, pf->phi);
    if (dr > 0.5) continue;
    if (dr < minDR) continue;
    if (dr >= maxDR) continue;

    if(printDebug) { 
      cout << "pf :: type: " << pf->pfType << "\tpt: " << pf->pt << "\tdR: " << dr;
      cout << "\tdZ: " << pf->dz;
      cout << endl;
    }

    //***********************************************************************
    // Charged Iso 
    //***********************************************************************
    if (pf->q != 0) {

      // Veto any PFmuon, or PFEle
      if ( pf->pfType == eElectron || pf->pfType == eMuon) {
        continue;
      }

      // Footprint Veto
      if (fabs(ele->scEta) > 1.479 && dr < 0.015) continue;

      if( printDebug) cout << "charged:: pt: " << pf->pt 
                           << "\ttype: " << pf->pfType
                           << endl;

      if (pfIsoType == kPFChargedIso) fIso += pf->pt;
    }

    //***********************************************************************
    // Gamma Iso 
    //***********************************************************************
    else if (pf->pfType == eGamma) {

      if (eleIndex >= 0) {
        //proper veto of electron footprint
        if (pf->matchedObjectType == 11 && Int_t(pf->matchedObjectIndex) == eleIndex) continue;
        //this is a temporary fix
        if (fabs(ele->scEta) < 1.479) {
          if (fabs(ele->eta - pf->eta) < 0.015) continue;
        }
      }

      if (fabs(ele->scEta) > 1.479) {
        if (higgsana::deltaR(ele->eta,ele->phi, pf->eta, pf->phi) < 0.08) continue;
      }
      if( printDebug) cout << "gamma:: " << pf->pt << " " 
                           << dr << endl;
      if (pfIsoType == kPFGammaIso) fIso += pf->pt;
    }

    //***********************************************************************
    // Other Neutrals
    //***********************************************************************
    else {
      if( printDebug ) cout << "neutral:: " << pf->pt << " " 
                           << dr << endl;
      if (pfIsoType == kPFNeutralHadronIso) fIso += pf->pt;
    }
  }

  if( printDebug ) cout << "rho: " << rho << endl;

  UInt_t EAType = 0;
  if (pfIsoType == kPFGammaIso) {
    if (minDR == 0.0 && maxDR == 0.1) EAType = kEleGammaIsoDR0p0To0p1;
    else if (minDR == 0.1 && maxDR == 0.2) EAType = kEleGammaIsoDR0p1To0p2;
    else if (minDR == 0.2 && maxDR == 0.3) EAType = kEleGammaIsoDR0p2To0p3;
    else if (minDR == 0.3 && maxDR == 0.4) EAType = kEleGammaIsoDR0p3To0p4;
    else if (minDR == 0.4 && maxDR == 0.5) EAType = kEleGammaIsoDR0p4To0p5;
    else EAType = kEleNoCorrection;
  } else if (pfIsoType == kPFNeutralHadronIso) {
    if (minDR == 0.0 && maxDR == 0.1) EAType = kEleNeutralHadronIsoDR0p0To0p1;
    else if (minDR == 0.1 && maxDR == 0.2) EAType = kEleNeutralHadronIsoDR0p1To0p2;
    else if (minDR == 0.2 && maxDR == 0.3) EAType = kEleNeutralHadronIsoDR0p2To0p3;
    else if (minDR == 0.3 && maxDR == 0.4) EAType = kEleNeutralHadronIsoDR0p3To0p4;
    else if (minDR == 0.4 && maxDR == 0.5) EAType = kEleNeutralHadronIsoDR0p4To0p5;
    else EAType = kEleNoCorrection;
  } else {
    EAType = kEleNoCorrection;
  }

  double iso = fmax(0.0, fIso - rho*ElectronEffectiveArea(EAType, ele->scEta,DataEra));


  if( printDebug ) { 
    cout << "uncorrected iso: " << fIso
	 << "\tcorrected iso: " << iso
	 << endl;
  }

  return iso;

}


Double_t ComputePhotonNonCorrectedIsoDr03(const higgsana::TPFCandidate *photon, 
                                     UInt_t leptonType,
                                     UInt_t leptonIndex,
                                     TClonesArray *fPFCandidates,
                                     Bool_t printDebug) {
  

  //
  // final iso 
  //
  Double_t fChargedIso  = 0.0;
  Double_t fGammaIso  = 0.0;
  Double_t fNeutralHadronIso  = 0.0;
  Double_t fpfPU  = 0.0;

  //
  // Loop over PF Candidates
  //
  for(int k=0; k<fPFCandidates->GetEntries(); ++k) {

    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*fPFCandidates)[k]);

    Double_t dr = higgsana::deltaR(photon->eta,photon->phi, pf->eta, pf->phi);
    if (dr > 0.3) continue;


    if (printDebug) cout << "photon pfiso: pfcandidate " << k << " : " << pf->pt << " " << pf->eta << " " << pf->phi 
                         << " : " << pf->matchedObjectType << " ==? " << leptonType << " : " 
                         << pf->matchedObjectIndex << " ==? " << leptonIndex << " : " 
                         << dr 
                         << " : " << pf->IsPFNoPU << " , " <<  pf->q << " "
                         << endl;
    
    if( !(pf->IsPFNoPU) && pf->q != 0 ) {
      if( pf->pt >= 0.2 && dr > 0.01 )
	fpfPU += pf->pt;
      continue;
    }
    
    //
    // skip this photon
    //
    if( pf == photon ) continue;
      
    //
    // Charged Iso 
    //
    if (pf->q != 0 ) {

      // if it matches the nearby lepton, don't include it in the iso
      if (pf->matchedObjectType == leptonType && pf->matchedObjectIndex == leptonIndex) continue;

      if( dr > 0.01 && pf->pt >= 0.2 )
	fChargedIso += pf->pt;
    }
    
    //
    // Gamma Iso 
    //
    else if (pf->pfType == eGamma) {
      if( pf->pt > 0.5 && dr > 0.01) 
	fGammaIso += pf->pt;
    }
    
    //
    // Other Neutrals
    //
    else {
      if( pf->pt > 0.5 && dr > 0.01) 
	fNeutralHadronIso += pf->pt;
    }
    
  }
  
  if( printDebug ) { 
    cout << "photon dbetaIso :: " << endl;
    cout << "\tfChargedIso: " << fChargedIso
	 << "\tfGammaIso: " << fGammaIso
	 << "\tfNeutralHadronIso: " << fNeutralHadronIso
	 << "\tfpfPU: " << fpfPU 
	 << endl;
  }
  double pfIso = fChargedIso + fGammaIso + fNeutralHadronIso + fpfPU;
  return pfIso/photon->pt;


}



Double_t ComputePhotonDBetaCorrectedRelativeIsoDr03(const higgsana::TPFCandidate *photon, 
                                       UInt_t leptonType,
                                       UInt_t leptonIndex,
                                       TClonesArray *fPFCandidates,
                                       Bool_t printDebug) {
  

  //
  // final iso 
  //
  Double_t fChargedIso  = 0.0;
  Double_t fGammaIso  = 0.0;
  Double_t fNeutralHadronIso  = 0.0;
  Double_t fpfPU  = 0.0;

  //
  // Loop over PF Candidates
  //
  for(int k=0; k<fPFCandidates->GetEntries(); ++k) {

    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*fPFCandidates)[k]);

    Double_t dr = higgsana::deltaR(photon->eta,photon->phi, pf->eta, pf->phi);
    if (dr > 0.3) continue;

    if( !(pf->IsPFNoPU) && pf->q != 0 ) {
      if( pf->pt >= 0.2 && dr > 0.01 )
	fpfPU += pf->pt;
      continue;
    }
    
    //
    // skip this photon
    //
    if( pf == photon ) continue;
      
    //
    // Charged Iso 
    //
    if (pf->q != 0 ) {

      // if it matches the nearby lepton, don't include it in the iso
      if (pf->matchedObjectType == leptonType && pf->matchedObjectIndex == leptonIndex) continue;

      if( dr > 0.01 && pf->pt >= 0.2 )
	fChargedIso += pf->pt;
    }
    
    //
    // Gamma Iso 
    //
    else if (pf->pfType == eGamma) {
      if( pf->pt > 0.5 && dr > 0.01) 
	fGammaIso += pf->pt;
    }
    
    //
    // Other Neutrals
    //
    else {
      if( pf->pt > 0.5 && dr > 0.01) 
	fNeutralHadronIso += pf->pt;
    }
    
  }
  
  if( printDebug ) { 
    cout << "photon dbetaIso :: " << endl;
    cout << "\tfChargedIso: " << fChargedIso
	 << "\tfGammaIso: " << fGammaIso
	 << "\tfNeutralHadronIso: " << fNeutralHadronIso
	 << "\tfpfPU: " << fpfPU 
	 << endl;
  }
  double pfIso = fChargedIso + fGammaIso + fNeutralHadronIso - 0.5*fpfPU;
  return pfIso/photon->pt;

}



//--------------------------------------------------------------------------------------------------
Bool_t passConversionVeto(Int_t isConv) {
 
  Bool_t pass = kFALSE;
  pass = ( (UInt_t(isConv) & 512) == 512);

  return pass; 
}




//--------------------------------------------------------------------------------------------------
Bool_t passMuonID(const higgsana::TMuon *muon, Double_t pfIso04)
{

  if(muon->nTkHits	  < 11)    return kFALSE;
  if(muon->nPixHits	  < 1)     return kFALSE;
  if(muon->pterr/muon->pt > 0.1)   return kFALSE;
  if(fabs(muon->dz)       > 0.1)   return kFALSE;
  if(muon->TrkKink        >= 20)   return kFALSE;

  Bool_t isGlobal  = ((muon->typeBits & kGlobal) == kGlobal) && (muon->muNchi2 < 10) && (muon->nMatch > 1) && (muon->nValidHits > 0);
  Bool_t isTracker = ((muon->typeBits & kTracker) == kTracker) && (muon->qualityBits & kTMLastStationTight);
  if(!isGlobal && !isTracker) return kFALSE;

  if(fabs(muon->d0)>0.02)   return kFALSE;
  if(fabs(muon->eta)<1.479) return (pfIso04<0.25*(muon->pt));
  else                      return (pfIso04<0.25*(muon->pt));

}



//--------------------------------------------------------------------------------------------------
Bool_t passCutBasedEleID(const higgsana::TElectron *electron, Double_t pfIso04)
{

  if(fabs(electron->d0) > 0.02) return kFALSE;
  if(fabs(electron->dz) > 0.1)  return kFALSE;
  
  // conversion rejection
  if(electron->nExpHitsInner > 0)           return kFALSE;
  if(!passConversionVeto(electron->isConv)) return kFALSE;
     
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {
    // barrel
    if(pfIso04 > 0.15*(electron->pt)) return kFALSE;
     
    if(electron->pt>20) {
      if(electron->sigiEtaiEta	    > 0.01)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.004) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.06)  return kFALSE;
    
    } else {
      if(electron->sigiEtaiEta	    > 0.01)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.004) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.03)  return kFALSE;
    }
  
  } else {
    // endcap
    if(pfIso04 > 0.15*(electron->pt)) return kFALSE;
     
    if(electron->pt>20) {
      if(electron->sigiEtaiEta	    > 0.03)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.007) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.03)  return kFALSE;
    
    } else {
      if(electron->sigiEtaiEta	    > 0.03)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.005) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.02)  return kFALSE;
    }
  }

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t passCutBasedTightEleID(const higgsana::TElectron *electron, Double_t pfIso04)
{

  if(fabs(electron->d0) > 0.02) return kFALSE;
  if(fabs(electron->dz) > 0.1)  return kFALSE;
  
  // conversion rejection
  if(electron->nExpHitsInner > 0)           return kFALSE;
  if(!passConversionVeto(electron->isConv)) return kFALSE;
     
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {
    // barrel
    if(pfIso04 > 0.10*(electron->pt)) return kFALSE;
     
    if(electron->pt>20) {
      if(electron->sigiEtaiEta	    > 0.01)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.004) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.06)  return kFALSE;
      if(electron->HoverE	    > 0.04)  return kFALSE;
    
    } else {
      if(electron->sigiEtaiEta	    > 0.01)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.004) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.03)  return kFALSE;
      if(electron->HoverE	    > 0.025) return kFALSE;    
    }
  
  } else {
    // endcap
    if(pfIso04 > 0.10*(electron->pt)) return kFALSE;
     
    if(electron->pt>20) {
      if(electron->sigiEtaiEta	    > 0.03)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.007) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.03)  return kFALSE;
      if(electron->HoverE	    > 0.10)  return kFALSE;
    
    } else {
      if(electron->sigiEtaiEta	    > 0.03)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.005) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.02)  return kFALSE;
      if(electron->HoverE	    > 0.10)  return kFALSE;      
    }
  }

  return kTRUE;
}


//--------------------------------------------------------------------------------------------------
Bool_t passCutBasedTightEleIDNoIso(const higgsana::TElectron *electron)
{

  if(fabs(electron->d0) > 0.02) return kFALSE;
  if(fabs(electron->dz) > 0.1)  return kFALSE;
  
  // conversion rejection
  if(electron->nExpHitsInner > 0)           return kFALSE;
  if(!passConversionVeto(electron->isConv)) return kFALSE;
     
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {
     
    if(electron->pt>20) {
      if(electron->sigiEtaiEta	    > 0.01)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.004) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.06)  return kFALSE;
      if(electron->HoverE	    > 0.04)  return kFALSE;
    
    } else {
      if(electron->sigiEtaiEta	    > 0.01)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.004) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.03)  return kFALSE;
      if(electron->HoverE	    > 0.025) return kFALSE;    
    }
  
  } else {
     
    if(electron->pt>20) {
      if(electron->sigiEtaiEta	    > 0.03)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.007) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.03)  return kFALSE;
      if(electron->HoverE	    > 0.10)  return kFALSE;
    
    } else {
      if(electron->sigiEtaiEta	    > 0.03)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.005) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.02)  return kFALSE;
      if(electron->HoverE	    > 0.10)  return kFALSE;      
    }
  }

  return kTRUE;
}


Bool_t passElectronTriggerDenominator(const higgsana::TElectron *ele) {

  Bool_t pass = kTRUE;
  if (fabs(ele->eta) >= 2.5) pass = kFALSE;
  
  //Barrel 
  if (fabs(ele->scEta) < 1.479) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01
            && fabs(ele->deltaEtaIn) < 0.007
            && fabs(ele->deltaPhiIn) < 0.15
            && ele->HoverE < 0.12
            && ele->nExpHitsInner <= 0
            && passConversionVeto(ele->isConv)
            && fabs(ele->d0) < 0.02
            && fabs(ele->dz) < 0.1
            && ( ele->trkIso03 ) / ele->pt < 0.2
            && ( fmax(ele->emIso03 - 1.0,0.0) ) / ele->pt < 0.20
            && (ele->hadIso03) / ele->pt < 0.20
          )
      ) {
      pass = kFALSE;
    }
  }
  //Endcap
  else {
    if (! (  (0==0)
             && ele->sigiEtaiEta < 0.03
             && fabs(ele->deltaEtaIn) < 0.009
             && fabs(ele->deltaPhiIn) < 0.10
             && ele->HoverE < 0.10
             && ele->nExpHitsInner <= 0
             && passConversionVeto(ele->isConv)
             && fabs(ele->d0) < 0.02
             && fabs(ele->dz) < 0.1
             && (ele->trkIso03 ) / ele->pt < 0.2
             && (ele->emIso03 ) / ele->pt < 0.20
             && (ele->hadIso03) / ele->pt < 0.20

          )
      ) {
      pass = kFALSE;
    }
  } 

  return pass;
}

Bool_t passMuonTriggerDenominator(const higgsana::TMuon *mu, TClonesArray *pfCandidates, 
                                  Double_t rho, UInt_t DataEra ) {

  Bool_t pass = kTRUE;
  if (fabs(mu->eta) >= 2.4) pass = kFALSE;
  
  if (! 
      ( (
          (Bool_t(mu->typeBits & kGlobal) 
           && mu->muNchi2 < 10.0
           && (mu->nValidHits > 0)
           && (mu->nMatch > 1 )
            )
          || 
          ( mu->typeBits & kTracker            
            && Bool_t(mu->qualityBits & kTMLastStationTight) 
            )
        )
        && mu->typeBits & kTracker
        && mu->nTkHits > 10
        && ( mu->nPixHits > 0)          
        && fabs(mu->d0) < 0.2
        && fabs(mu->dz) < 0.1          
        && ComputeMuonPFIso04( mu, pfCandidates,  rho, DataEra, kFALSE) / mu->pt < 0.7
        && ( mu->pterr / mu->pt < 0.1)
        && (mu->TrkKink < 20)
        )
    ) pass = kFALSE;    

  return pass;
}

Bool_t PassEleSimpleCutsLoose( const higgsana::TElectron *ele,   
                               TClonesArray *pfCandidates, 
                               Double_t rho, UInt_t DataEra, 
                               Bool_t printDebug) {

  Bool_t pass = kTRUE;
  if (ele->isEB) {
    if (!(
          fabs(ele->deltaEtaIn) < 0.007
          && fabs(ele->deltaPhiIn) < 0.15
          && ele->sigiEtaiEta < 0.01
          && ele->HoverE < 0.12
          && fabs(ele->d0) < 0.02
          && fabs(ele->dz) < 0.2
          && fabs(1.0/(ele->scEt*TMath::CosH(ele->scEta)) - ele->EOverP/ele->EcalEnergy) < 0.05
          && ele->nExpHitsInner <= 1
          && ((ele->isConv & 512) == 512)
          && ComputeElePFIso04( ele, pfCandidates, rho, DataEra, printDebug)/ele->pt < 0.20
          )
      ) pass = kFALSE;
  } else {
    if (!(
          fabs(ele->deltaEtaIn) < 0.009
          && fabs(ele->deltaPhiIn) < 0.10
          && ele->sigiEtaiEta < 0.03
          && ele->HoverE < 0.10
          && fabs(ele->d0) < 0.02
          && fabs(ele->dz) < 0.2
          && fabs(1.0/(ele->scEt*TMath::CosH(ele->scEta)) - ele->EOverP/ele->EcalEnergy) < 0.05
          && ele->nExpHitsInner <= 1
          && ((ele->isConv & 512) == 512)
          && ComputeElePFIso04( ele, pfCandidates, rho, DataEra, printDebug)/ele->pt < 0.15
          )
      ) pass = kFALSE;
  }

  return pass;

}

Bool_t PassEleSimpleCutsMedium( const higgsana::TElectron *ele,   
                               TClonesArray *pfCandidates, 
                               Double_t rho, UInt_t DataEra, 
                               Bool_t printDebug) {

  Bool_t pass = kTRUE;
  if (ele->isEB) {
    if (!(
          fabs(ele->deltaEtaIn) < 0.004
          && fabs(ele->deltaPhiIn) < 0.06
          && ele->sigiEtaiEta < 0.01
          && ele->HoverE < 0.12
          && fabs(ele->d0) < 0.02
          && fabs(ele->dz) < 0.1
          && fabs(1.0/(ele->scEt*TMath::CosH(ele->scEta)) - ele->EOverP/ele->EcalEnergy) < 0.05
          && ele->nExpHitsInner <= 1
          && ((ele->isConv & 512) == 512)
          && ComputeElePFIso04( ele, pfCandidates, rho, DataEra, printDebug)/ele->pt < 0.20
          )
      ) pass = kFALSE;
  } else {
    if (!(
          fabs(ele->deltaEtaIn) < 0.007
          && fabs(ele->deltaPhiIn) < 0.03
          && ele->sigiEtaiEta < 0.03
          && ele->HoverE < 0.10
          && fabs(ele->d0) < 0.02
          && fabs(ele->dz) < 0.1
          && fabs(1.0/(ele->scEt*TMath::CosH(ele->scEta)) - ele->EOverP/ele->EcalEnergy) < 0.05
          && ele->nExpHitsInner <= 1
          && ((ele->isConv & 512) == 512)
          && ComputeElePFIso04( ele, pfCandidates, rho, DataEra, printDebug)/ele->pt < 0.15
          )
      ) pass = kFALSE;
  }

  return pass;

}

Bool_t PassEleSimpleCutsTight( const higgsana::TElectron *ele,   
                               TClonesArray *pfCandidates, 
                               Double_t rho, UInt_t DataEra, 
                               Bool_t printDebug) {

  Bool_t pass = kTRUE;
  if (ele->isEB) {
    if (!(
          fabs(ele->deltaEtaIn) < 0.004
          && fabs(ele->deltaPhiIn) < 0.03
          && ele->sigiEtaiEta < 0.01
          && ele->HoverE < 0.12
          && fabs(ele->d0) < 0.02
          && fabs(ele->dz) < 0.1
          && fabs(1.0/(ele->scEt*TMath::CosH(ele->scEta)) - ele->EOverP/ele->EcalEnergy) < 0.05
          && ele->nExpHitsInner <= 0
          && ((ele->isConv & 512) == 512)
          && ComputeElePFIso04( ele, pfCandidates, rho, DataEra, printDebug)/ele->pt < 0.15
          )
      ) pass = kFALSE;
  } else {
    if (!(
          fabs(ele->deltaEtaIn) < 0.005
          && fabs(ele->deltaPhiIn) < 0.02
          && ele->sigiEtaiEta < 0.03
          && ele->HoverE < 0.10
          && fabs(ele->d0) < 0.02
          && fabs(ele->dz) < 0.1
          && fabs(1.0/(ele->scEt*TMath::CosH(ele->scEta)) - ele->EOverP/ele->EcalEnergy) < 0.05
          && ele->nExpHitsInner <= 0
          && ((ele->isConv & 512) == 512)
          && ComputeElePFIso04( ele, pfCandidates, rho, DataEra, printDebug)/ele->pt < 0.10
          )
      ) pass = kFALSE;
  }

  return pass;

}


Bool_t PassEleSimpleCutsVeryTight( const higgsana::TElectron *ele,   
                               TClonesArray *pfCandidates, 
                               Double_t rho, UInt_t DataEra, 
                               Bool_t printDebug) {

  Bool_t pass = kTRUE;
  if (ele->isEB) {
    if (!(
          fabs(ele->deltaEtaIn) < 0.004
          && fabs(ele->deltaPhiIn) < 0.03
          && ele->sigiEtaiEta < 0.01
          && ele->HoverE < 0.12
          && fabs(ele->d0) < 0.02
          && fabs(ele->dz) < 0.1
          && fabs(1.0/(ele->scEt*TMath::CosH(ele->scEta)) - ele->EOverP/ele->EcalEnergy) < 0.05
          && ele->nExpHitsInner <= 0
          && ((ele->isConv & 512) == 512)
          && ComputeElePFIso04( ele, pfCandidates, rho, DataEra, printDebug)/ele->pt < 0.075
          )
      ) pass = kFALSE;
  } else {
    if (!(
          fabs(ele->deltaEtaIn) < 0.005
          && fabs(ele->deltaPhiIn) < 0.02
          && ele->sigiEtaiEta < 0.03
          && ele->HoverE < 0.10
          && fabs(ele->d0) < 0.02
          && fabs(ele->dz) < 0.1
          && fabs(1.0/(ele->scEt*TMath::CosH(ele->scEta)) - ele->EOverP/ele->EcalEnergy) < 0.05
          && ele->nExpHitsInner <= 0
          && ((ele->isConv & 512) == 512)
          && ComputeElePFIso04( ele, pfCandidates, rho, DataEra, printDebug)/ele->pt < 0.075
          )
      ) pass = kFALSE;
  }

  return pass;

}

Bool_t passPhotonSimpleCuts(const higgsana::TPhoton *pho) {

  Bool_t pass = kTRUE;

  //Barrel 
  if (fabs(pho->eta) < 1.479) {
    if (! ( pho->sigiEtaiEta < 0.0105
            && pho->R9 > 0.9
          )
      ) {
      pass = kFALSE;
    }
  }
  //Endcap
  else {
    if (! (  pho->sigiEtaiEta < 0.028
             && pho->R9 > 0.9
          )
      ) {
      pass = kFALSE;
    }
  }

  return pass;
}

#endif
