#ifndef HIGGSANA_HZZ4L_UTILS_LEPTONSELECTION_HH
#define HIGGSANA_HZZ4L_UTILS_LEPTONSELECTION_HH

#include "TLorentzVector.h"
#include "HiggsAna/Utils/CommonDefs.hh"
#include "HiggsAna/Utils/CommonTools.hh"
#include "HiggsAna/Utils/IsolationPileupCorrections.hh"
#include "HiggsAna/DataTree/interface/TElectron.hh"
#include "HiggsAna/DataTree/interface/TMuon.hh"
#include "HiggsAna/DataTree/interface/HiggsAnaDefs.hh"
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"

Int_t Classify(const higgsana::TElectron *ele);
Int_t Classify(Double_t eta, Double_t eOverP, Double_t fBrem, Bool_t isEB, Bool_t isEcalDriven);

bool compute_cut(double x, double et, double cut_min, double cut_max, bool gtn = false);

Bool_t passMuonID_HZZ2011(const higgsana::TMuon *muon);
Bool_t passMuonIP_HZZ2011(const higgsana::TMuon *muon);
Int_t  PassCiCID(const higgsana::TElectron *ele, const Int_t typeCuts = 1);
Bool_t passEleWP95ID(const higgsana::TElectron *electron);
Bool_t PassMuonHZZ4lPreselectionWithoutPtCut( const higgsana::TMuon *muon);  
Bool_t PassMuonHZZ4lPreselection( const higgsana::TMuon *muon);  
Bool_t PassMuonHZZ4lICHEP2012ID( const higgsana::TMuon *muon, UInt_t muonIndex, TClonesArray *pfCandidates);  
Bool_t PassMuonHZZ4lICHEP2012Iso( const higgsana::TMuon *muon, UInt_t muonIndex, TClonesArray *pfCandidates, 
                                  Double_t rho, UInt_t DataEra, 
                                  vector<const higgsana::TPFCandidate*> photonsToVeto, Bool_t printDebug = kFALSE);
Bool_t PassEleHZZ4lPreselectionWithoutPtCut( const higgsana::TElectron *ele);  
Bool_t PassEleHZZ4lPreselection( const higgsana::TElectron *ele);  
Double_t EvaluateEleHZZ4lICHEP2012IDMVA( const higgsana::TElectron *ele, EGammaMvaEleEstimator* eleMVAEstimator,
                                         Bool_t printDebug = kFALSE);  
Bool_t PassEleHZZ4lICHEP2012ID( const higgsana::TElectron *ele, EGammaMvaEleEstimator* eleMVAEstimator,
                                Bool_t printDebug = kFALSE);  
Bool_t PassEleHZZ4lRun1LegacyPaperID( const higgsana::TElectron *ele, EGammaMvaEleEstimator* eleMVAEstimator,
                                Bool_t printDebug = kFALSE);  
Bool_t PassEleHZZ4lICHEP2012Iso( const higgsana::TElectron *ele, UInt_t eleIndex,
                                 TClonesArray *pfCandidates, 
                                 Double_t rho, UInt_t DataEra, 
                                 vector<const higgsana::TPFCandidate*> photonsToVeto, Bool_t printDebug = kFALSE);  
Bool_t PassEleHZZ4lICHEP2012TighterIso( const higgsana::TElectron *ele, UInt_t eleIndex,
                                        TClonesArray *pfCandidates, 
                                        Double_t rho, UInt_t DataEra, 
                                        vector<const higgsana::TPFCandidate*> photonsToVeto, Bool_t printDebug = kFALSE);  

//=== FUNCTION DEFINITIONS ======================================================================================

Bool_t PassMuonHZZ4lPreselectionWithoutPtCut( const higgsana::TMuon *muon) {

  Bool_t pass = kTRUE;
  if ( fabs(muon->ip3dSig) >= 100) pass = kFALSE;
  if ( fabs(muon->eta) > 2.4) pass = kFALSE;
  if ( ! ( (muon->typeBits & kGlobal) == kGlobal 
           || 
           (muon->typeBits & kTracker) == kTracker
         )
    ) pass = kFALSE;
  
  if ( fabs(muon->d0) >= 0.5) pass = kFALSE;
  if ( fabs(muon->dz) >= 1.0) pass = kFALSE;

  return pass;
}

Bool_t PassMuonHZZ4lPreselection( const higgsana::TMuon *muon) {

  Bool_t pass = kTRUE;
  if ( fabs(muon->ip3dSig) >= 100) pass = kFALSE;
  if ( muon->pt < 5) pass = kFALSE;
  if ( fabs(muon->eta) > 2.4) pass = kFALSE;
  if ( ! ( (muon->typeBits & kGlobal) == kGlobal 
           || 
           (muon->typeBits & kTracker) == kTracker
         )
    ) pass = kFALSE;
  
  if ( fabs(muon->d0) >= 0.5) pass = kFALSE;
  if ( fabs(muon->dz) >= 1.0) pass = kFALSE;

  return pass;
}

Bool_t PassMuonHZZ4lICHEP2012ID( const higgsana::TMuon *muon, UInt_t muonIndex, TClonesArray *pfCandidates) {

  Bool_t pass = kFALSE;

  // check that it matches to a PF muon
  for( UInt_t i=0; i < UInt_t(pfCandidates->GetEntries()); i++ ) {

    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*pfCandidates)[i]);        

    if (!pf->IsPFNoPU) continue;
    
    if( pf->matchedObjectType == 13 && pf->matchedObjectIndex == muonIndex 
        && pf->pfType == eMuon ) {
      pass = kTRUE;
      break;
    }  
  }

  //apply IP3D at this stage too
  if (!(fabs(muon->ip3dSig) < 4)) pass = kFALSE;
  
  return pass;
}


Bool_t PassMuonHZZ4lICHEP2012Iso( const higgsana::TMuon *muon, UInt_t muonIndex,
                                  TClonesArray *pfCandidates, 
                                  Double_t rho, UInt_t DataEra, 
                                  vector<const higgsana::TPFCandidate*> photonsToVeto,
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
    
    //*****************************
    // veto FSR recovered photons
    //*****************************
    bool vetoPhoton = false;
    for( UInt_t p=0; p<photonsToVeto.size(); p++ ) {
      if( pf == photonsToVeto[p] ) { 
	vetoPhoton = true;
	break;
      }
    }
    if( vetoPhoton ) continue;

//     Double_t deta = (muon->eta - pf->eta);
//     Double_t dphi = fabs(deltaPhi(muon->phi, pf->phi));
    Double_t dr = higgsana::deltaR(muon->eta, muon->phi, pf->eta, pf->phi);
    if (dr > 0.4) continue;

    //remove pf candidate matching to reco muon
    if (pf->matchedObjectType == 13 && pf->matchedObjectIndex == muonIndex) continue;


    //***********************************************************************
    // Charged Iso 
    //***********************************************************************
    if (pf->q != 0) {      
      //don't include pf electrons or pf muons
      if (pf->pfType == eElectron || pf->pfType == eMuon) continue;

      fChargedIso += pf->pt;
    }
    
    //***********************************************************************
    // Gamma Iso 
    //***********************************************************************
    else if (pf->pfType == eGamma) {
      if( pf->pt > 0.5 )
      fGammaIso += pf->pt;
    }
    
    //***********************************************************************
    // Other Neutrals
    //***********************************************************************
    else {
      if( pf->pt > 0.5 ) 
	fNeutralHadronIso += pf->pt;
    }
  }

  if(printDebug) cout << "rho: " << rho << endl;

  //modify
  TLorentzVector  tmpvec;
  tmpvec.SetPtEtaPhiM(muon->pt,muon->eta,muon->phi,MUONMASS);
  for( UInt_t p=0; p<photonsToVeto.size(); p++ ) {
    const higgsana::TPFCandidate *pfphoton  = photonsToVeto[p];
    TLorentzVector pfphotonvec;
    pfphotonvec.SetPtEtaPhiM(pfphoton->pt,pfphoton->eta,pfphoton->phi,0.);
    tmpvec += pfphotonvec;
  }

  double pfIso = fChargedIso + fmax(0.0,(fGammaIso + fNeutralHadronIso 
                                         -rho*MuonEffectiveArea(kMuGammaAndNeutralHadronIso04,
                                                                tmpvec.Eta(),DataEra)));
								//muon->Eta(),DataEra)));
  
  if( printDebug ) { 
    cout << "PFiso: " << pfIso
	 << "\tfChargedIso: " << fChargedIso
	 << "\tfGammaIso: " << fGammaIso
	 << "\tfNeutralHadronIso: " << fNeutralHadronIso
	 << endl;
  }


  Bool_t pass = ( (pfIso / muon->pt) < 0.40 );
  return pass;

}



Bool_t PassEleHZZ4lPreselectionWithoutPtCut( const higgsana::TElectron *ele) {

  Bool_t pass = kTRUE;
  if ( fabs(ele->ip3dSig) >= 100) pass = kFALSE;
  if ( fabs(ele->eta) >= 2.5) pass = kFALSE;
  if ( ele->nExpHitsInner > 1) pass = kFALSE;
  
  if ( fabs(ele->d0) >= 0.5) pass = kFALSE;
  if ( fabs(ele->dz) >= 1.0) pass = kFALSE;

  return pass;


}


Bool_t PassEleHZZ4lPreselection( const higgsana::TElectron *ele) {

  Bool_t pass = kTRUE;
  if ( fabs(ele->ip3dSig) >= 100) pass = kFALSE;
  if ( ele->pt < 7) pass = kFALSE;
  if ( fabs(ele->eta) >= 2.5) pass = kFALSE;
  if ( ele->nExpHitsInner > 1) pass = kFALSE;
  
  if ( fabs(ele->d0) >= 0.5) pass = kFALSE;
  if ( fabs(ele->dz) >= 1.0) pass = kFALSE;

  return pass;


}

  

Double_t EvaluateEleHZZ4lICHEP2012IDMVA( const higgsana::TElectron *ele, EGammaMvaEleEstimator* eleMVAEstimator,
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



Bool_t PassEleHZZ4lICHEP2012ID( const higgsana::TElectron *ele, EGammaMvaEleEstimator* eleMVAEstimator,
                                Bool_t printDebug) {

  double mvaval = EvaluateEleHZZ4lICHEP2012IDMVA(ele,eleMVAEstimator, printDebug);
  
  bool pass = false;
  
  Int_t subdet = 0;
  if (fabs(ele->scEta) < 0.8) subdet = 0;
  else if (fabs(ele->scEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ele->pt > 10) ptBin = 1;
	
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0; // eta<0.8, pt<10
  if (subdet == 1 && ptBin == 0) MVABin = 1; // 0.8<eta<1.479, pt<10
  if (subdet == 2 && ptBin == 0) MVABin = 2; // eta>1.478, pt<10
  if (subdet == 0 && ptBin == 1) MVABin = 3; // eta<0.8, pt>10
  if (subdet == 1 && ptBin == 1) MVABin = 4; // 0.8<eta<1.479, pt>10
  if (subdet == 2 && ptBin == 1) MVABin = 5; // eta>1.478, pt>10    

  if( MVABin == 0 && mvaval > 0.470 ) pass = true;
  if( MVABin == 1 && mvaval > 0.004 ) pass = true;
  if( MVABin == 2 && mvaval > 0.295 ) pass = true;
  if( MVABin == 3 && mvaval > 0.500 ) pass = true;
  if( MVABin == 4 && mvaval > 0.120 ) pass = true;
  if( MVABin == 5 && mvaval > 0.600 ) pass = true;


  //apply IP3D at this stage too
  if (!(fabs(ele->ip3dSig) < 4)) pass = kFALSE;

  return pass;

}

Bool_t PassEleHZZ4lRun1LegacyPaperID( const higgsana::TElectron *ele, EGammaMvaEleEstimator* eleMVAEstimator,
                                      Bool_t printDebug) {
  
  double mvaval = EvaluateEleHZZ4lICHEP2012IDMVA(ele,eleMVAEstimator, printDebug);
  
  bool pass = false;
  
  Int_t subdet = 0;
  if (fabs(ele->scEta) < 0.8) subdet = 0;
  else if (fabs(ele->scEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ele->pt > 10) ptBin = 1;
	
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0; // eta<0.8, pt<10
  if (subdet == 1 && ptBin == 0) MVABin = 1; // 0.8<eta<1.479, pt<10
  if (subdet == 2 && ptBin == 0) MVABin = 2; // eta>1.478, pt<10
  if (subdet == 0 && ptBin == 1) MVABin = 3; // eta<0.8, pt>10
  if (subdet == 1 && ptBin == 1) MVABin = 4; // 0.8<eta<1.479, pt>10
  if (subdet == 2 && ptBin == 1) MVABin = 5; // eta>1.478, pt>10    

  if( MVABin == 0 && mvaval > 0.470 ) pass = true;
  if( MVABin == 1 && mvaval > 0.004 ) pass = true;
  if( MVABin == 2 && mvaval > 0.295 ) pass = true;
  if( MVABin == 3 && mvaval > -0.34 ) pass = true;
  if( MVABin == 4 && mvaval > -0.65 ) pass = true;
  if( MVABin == 5 && mvaval > 0.600 ) pass = true;


  //apply IP3D at this stage too
  if (!(fabs(ele->ip3dSig) < 4)) pass = kFALSE;

  return pass;

}



Bool_t PassEleHZZ4lICHEP2012Iso( const higgsana::TElectron *ele, 
                                 UInt_t eleIndex,
                                 TClonesArray *pfCandidates, 
                                 Double_t rho, UInt_t DataEra, 
                                 vector<const higgsana::TPFCandidate*> photonsToVeto, 
                                 Bool_t printDebug) {
 
  
  Double_t fChargedIso  = 0.0;
  Double_t fGammaIso  = 0.0;
  Double_t fNeutralHadronIso  = 0.0;

  
  //
  //Loop over PF Candidates
  //
  for(UInt_t k=0; int(k) < pfCandidates->GetEntries(); ++k) {

    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*pfCandidates)[k]);

//     Double_t deta = (ele->eta - pf->eta);
//     Double_t dphi = fabs(deltaPhi(ele->phi, pf->phi));
    Double_t dr = higgsana::deltaR(ele->eta, ele->phi, pf->eta, pf->phi);

    if(printDebug) { 
      cout << "pf :: type: " << pf->pfType << "\tpt: " << pf->pt << "\tdR: " << dr;
      cout << "\tdZ: " << pf->dz
           << "\tPFnoPU: " << pf->IsPFNoPU ;
      cout << endl;
    }

    //*****************************
    // veto Pileup PF Candidates
    //*****************************
    if (!pf->IsPFNoPU) continue;

    //*****************************
    // veto FSR recovered photons
    //*****************************
    bool vetoPhoton = false;
    for( UInt_t p=0; p<photonsToVeto.size(); p++ ) {
      if( pf == photonsToVeto[p] ) { 
	vetoPhoton = true;
	break;
      }
    }

    if( vetoPhoton ) continue;
    if (dr > 0.4) continue;


    //***********************************************************************************************************
    // veto PF candidates that match with the electron supercluster if electron has more than 1 missing hit
    // there's also some cut on pf.mva_nothing_gamma() > 0.99, but let me ignore this for now. 
    // emulate this by requiring it's a pfgamma
    //***********************************************************************************************************
    if (ele->nExpHitsInner > 0) {
      if (pf->pfType == eGamma && pf->matchedObjectType == 11 && pf->matchedObjectIndex == eleIndex) {
        continue;
      }
    }

    //
    // sync : I don't think theyre doing this ...
    //
    //     if ( (pf->HasTrackerTrk() && (pf->TrackerTrk() == ele->TrackerTrk())) ||
    // 	 (pf->HasGsfTrk() && (pf->GsfTrk() == ele->GsfTrk()))) { 
    //       if( ctrl.debug ) cout << "\tskipping, matches to the electron ..."  << endl;
    //       continue;
    //     }

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

  Bool_t pass = ( (pfIso / ele->pt) < 0.40 );
  return pass;


}




Bool_t PassEleHZZ4lICHEP2012TighterIso( const higgsana::TElectron *ele, 
                                 UInt_t eleIndex,
                                 TClonesArray *pfCandidates, 
                                 Double_t rho, UInt_t DataEra, 
                                 vector<const higgsana::TPFCandidate*> photonsToVeto, 
                                 Bool_t printDebug) {
 
  
  Double_t fChargedIso  = 0.0;
  Double_t fGammaIso  = 0.0;
  Double_t fNeutralHadronIso  = 0.0;

  
  //
  //Loop over PF Candidates
  //
  for(UInt_t k=0; int(k) < pfCandidates->GetEntries(); ++k) {

    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*pfCandidates)[k]);

//     Double_t deta = (ele->eta - pf->eta);
//     Double_t dphi = fabs(deltaPhi(ele->phi, pf->phi));
    Double_t dr = higgsana::deltaR(ele->eta, ele->phi, pf->eta, pf->phi);

    if(printDebug) { 
      cout << "pf :: type: " << pf->pfType << "\tpt: " << pf->pt << "\tdR: " << dr;
      cout << "\tdZ: " << pf->dz
           << "\tPFnoPU: " << pf->IsPFNoPU ;
      cout << endl;
    }

    //*****************************
    // veto Pileup PF Candidates
    //*****************************
    if (!pf->IsPFNoPU) continue;

    //*****************************
    // veto FSR recovered photons
    //*****************************
    bool vetoPhoton = false;
    for( UInt_t p=0; p<photonsToVeto.size(); p++ ) {
      if( pf == photonsToVeto[p] ) { 
	vetoPhoton = true;
	break;
      }
    }

    if( vetoPhoton ) continue;
    if (dr > 0.4) continue;


    //***********************************************************************************************************
    // veto PF candidates that match with the electron supercluster if electron has more than 1 missing hit
    // there's also some cut on pf.mva_nothing_gamma() > 0.99, but let me ignore this for now. 
    // emulate this by requiring it's a pfgamma
    //***********************************************************************************************************
    if (ele->nExpHitsInner > 0) {
      if (pf->pfType == eGamma && pf->matchedObjectType == 11 && pf->matchedObjectIndex == eleIndex) {
        continue;
      }
    }

    //
    // sync : I don't think theyre doing this ...
    //
    //     if ( (pf->HasTrackerTrk() && (pf->TrackerTrk() == ele->TrackerTrk())) ||
    // 	 (pf->HasGsfTrk() && (pf->GsfTrk() == ele->GsfTrk()))) { 
    //       if( ctrl.debug ) cout << "\tskipping, matches to the electron ..."  << endl;
    //       continue;
    //     }

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

  Bool_t pass = ( (pfIso / ele->pt) < 0.15 );
  return pass;


}





Int_t Classify(Double_t eta, Double_t eOverP, Double_t fBrem, Bool_t isEB, Bool_t isEcalDriven) {
  
  int cat = -1;
  if (isEB == kTRUE) {
    if ((fBrem >= 0.12) && (eOverP > 0.9) && (eOverP < 1.2))
      cat = 0;
    else if (((eta >  .445   && eta <  .45  ) ||
  	      (eta >  .79    && eta <  .81  ) ||
  	      (eta > 1.137   && eta < 1.157 ) ||
  	      (eta > 1.47285 && eta < 1.4744)))
      cat = 6;
    else if (!isEcalDriven)
      cat = 8;
    else if (fBrem < 0.12)
      cat = 1;
    else
      cat = 2;
  } else {
    if ((fBrem >= 0.2) && (eOverP > 0.82) && (eOverP < 1.22))
      cat = 3;
    else if (eta > 1.5 && eta <  1.58)
      cat = 7;
    else if (!isEcalDriven)
      cat = 8;
    else if (fBrem < 0.2)
      cat = 4;
    else
      cat = 5;
  }

  return cat;
}

Int_t Classify(const higgsana::TElectron *ele) {
  return Classify(ele->scEta, ele->EOverP, ele->fBrem, ele->isEB, ele->isEcalDriven);
}


//--------------------------------------------------------------------------------------------------
bool compute_cut(double x, double et, double cut_min, double cut_max, bool gtn) {

  float et_min = 10;
  float et_max = 40;

  bool accept = false;
  float cut = cut_max; //  the cut at et=40 GeV

  if(et < et_max) {
    cut = cut_min + (1/et_min - 1/et)*(cut_max - cut_min)/(1/et_min - 1/et_max);
  } 
  
  if(et < et_min) {
    cut = cut_min;
  } 

  if(gtn) {   // useful for e/p cut which is gt
    accept = (x >= cut);
  } 
  else {
    accept = (x <= cut);
  }

  return accept;
}


//--------------------------------------------------------------------------------------------------
Int_t PassCiCID( Double_t pt, Double_t scEt, Double_t scEta, Bool_t isEcalDriven, Bool_t isEB, 
                 Double_t fBrem, Double_t EOverP, Double_t hOverE, Double_t sigmaee, Double_t deltaPhiIn, Double_t deltaEtaIn, 
                 Double_t eSeedOverPin, Double_t mishits, 
                 Double_t tkIso, Double_t ecalIso, Double_t hcalIso, 
                 Double_t partnerDist, Double_t partnerDeltaCot, 
                 Double_t d0, Double_t ip3dSig,
                 const Int_t typeCuts, Bool_t debug = kFALSE) {


  if (debug) {
    cout <<  pt << " "  <<  scEt << " "  <<  scEta << " " << isEcalDriven << " " << isEB << " " 
         <<  fBrem << " "  <<  EOverP << " "  <<  hOverE << " "  <<  sigmaee << " "  <<  deltaPhiIn << " "  <<  deltaEtaIn << " " 
         <<  eSeedOverPin << " "  <<  mishits << " " 
         <<  tkIso << " "  <<  ecalIso << " "  <<  hcalIso << " " 
         <<  partnerDist << " "  <<  partnerDeltaCot << " " 
         <<  d0 << " "  <<  ip3dSig << " " << typeCuts << endl;
  }

  int cat = Classify(scEta, EOverP, fBrem, isEB, isEcalDriven);
//   int eb;

//   if (isEB) 
//     eb = 0;
//   else 
//     eb = 1; 

  // Medium cuts
  Double_t cutdcotdistMedium[9] = {
  3.32e-02, 2.92e-02, 2.49e-02, 3.92e-02, 3.41e-02, 3.96e-02, 2.91e-02, 3.95e-02, 7.71e-03};
  Double_t cutdetainMedium[9] = {
  1.33e-02, 4.48e-03, 9.22e-03, 1.54e-02, 7.26e-03, 1.24e-02, 1.29e-02, 3.84e-02, 1.88e-02};
  Double_t cutdetainlMedium[9] = {
  1.21e-02, 4.22e-03, 9.18e-03, 1.61e-02, 6.45e-03, 1.16e-02, 1.23e-02, 6.20e-02, 2.43e-02};
  Double_t cutdphiinMedium[9] = {
  7.09e-02, 2.43e-01, 2.96e-01, 7.98e-02, 2.35e-01, 2.76e-01, 3.42e-01, 4.04e-01, 2.99e-01};
  Double_t cutdphiinlMedium[9] = {
  7.42e-02, 2.43e-01, 2.97e-01, 9.12e-02, 2.26e-01, 2.76e-01, 3.34e-01, 5.58e-01, 2.91e-01};
  Double_t cuteseedopcorMedium[9] = {
  6.42e-01, 9.44e-01, 4.53e-01, 7.62e-01, 3.67e-01, 5.57e-01, 1.98e-01, 9.15e-01, 6.28e-02};
  Double_t cutfmishitsMedium[9] = {
  4.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01};
  Double_t cuthoeMedium[9] = {
  1.96e-01, 6.30e-02, 1.48e-01, 3.66e-01, 5.66e-02, 1.45e-01, 4.29e-01, 4.28e-01, 3.99e-01};
  Double_t cuthoelMedium[9] = {
  2.19e-01, 6.19e-02, 1.47e-01, 3.58e-01, 4.61e-02, 1.46e-01, 3.26e-01, 3.81e-01, 3.89e-01};
  Double_t cutip_gsfMedium[9] = {
  2.45e-02, 9.74e-02, 1.48e-01, 5.49e-02, 5.65e-01, 3.33e-01, 2.04e-01, 5.41e-01, 1.21e-01};
  Double_t cutip_gsflMedium[9] = {
  1.92e-02, 9.81e-02, 1.33e-01, 4.34e-02, 5.65e-01, 3.24e-01, 2.33e-01, 4.30e-01, 6.44e-02};
  Double_t cutiso_sumMedium[9] = {
  1.44e+01, 1.12e+01, 1.09e+01, 1.08e+01, 6.35e+00, 9.78e+00, 1.30e+01, 1.62e+01, 1.96e+00};
  Double_t cutiso_sumoetMedium[9] = {
  1.01e+01, 6.41e+00, 6.00e+00, 8.14e+00, 3.90e+00, 4.76e+00, 6.86e+00, 6.48e+00, 1.74e+01};
  Double_t cutiso_sumoetlMedium[9] = {
  9.44e+00, 7.67e+00, 7.15e+00, 7.34e+00, 3.35e+00, 4.70e+00, 8.32e+00, 7.55e+00, 6.25e+00};
  Double_t cutseeMedium[9] = {
  1.30e-02, 1.09e-02, 1.18e-02, 3.94e-02, 3.04e-02, 3.28e-02, 1.00e-02, 3.73e-02, 6.69e-02};
  Double_t cutseelMedium[9] = {
  1.42e-02, 1.11e-02, 1.29e-02, 4.32e-02, 2.96e-02, 3.82e-02, 1.01e-02, 4.45e-02, 1.19e-01};

  // Tight cuts
  Double_t cutdcotdistTight[9] = {
  2.68e-02, 2.36e-02, 2.21e-02, 3.72e-02, 3.17e-02, 3.61e-02, 2.55e-02, 3.75e-02, 2.16e-04};
  Double_t cutdetainTight[9] = {
  8.92e-03, 3.96e-03, 8.50e-03, 1.34e-02, 6.27e-03, 1.05e-02, 1.12e-02, 3.09e-02, 1.88e-02};
  Double_t cutdetainlTight[9] = {
  9.23e-03, 3.77e-03, 8.70e-03, 1.39e-02, 5.60e-03, 9.40e-03, 1.07e-02, 6.20e-02, 4.10e-03};
  Double_t cutdphiinTight[9] = {
  6.37e-02, 1.53e-01, 2.90e-01, 7.69e-02, 1.81e-01, 2.34e-01, 3.42e-01, 3.93e-01, 2.84e-01};
  Double_t cutdphiinlTight[9] = {
  6.92e-02, 2.33e-01, 2.96e-01, 8.65e-02, 1.85e-01, 2.76e-01, 3.34e-01, 3.53e-01, 2.90e-01};
  Double_t cuteseedopcorTight[9] = {
  6.52e-01, 9.69e-01, 9.12e-01, 7.79e-01, 3.67e-01, 6.99e-01, 3.28e-01, 9.67e-01, 5.89e-01};
  Double_t cutfmishitsTight[9] = {
  4.50e+00, 1.50e+00, 5.00e-01, 1.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
  Double_t cuthoeTight[9] = {
  1.74e-01, 4.88e-02, 1.46e-01, 3.64e-01, 4.93e-02, 1.45e-01, 4.29e-01, 4.20e-01, 3.99e-01};
  Double_t cuthoelTight[9] = {
  2.19e-01, 5.25e-02, 1.47e-01, 3.57e-01, 4.25e-02, 1.45e-01, 3.26e-01, 3.80e-01, 1.32e-01};
  Double_t cutip_gsfTight[9] = {
  1.58e-02, 8.25e-02, 1.15e-01, 4.05e-02, 5.40e-01, 1.51e-01, 7.74e-02, 4.17e-01, 7.80e-02};
  Double_t cutip_gsflTight[9] = {
  1.27e-02, 6.26e-02, 9.68e-02, 3.02e-02, 5.65e-01, 1.46e-01, 7.90e-02, 4.10e-01, 4.79e-02};
  Double_t cutiso_sumTight[9] = {
  1.23e+01, 9.77e+00, 1.01e+01, 9.77e+00, 6.13e+00, 7.55e+00, 1.30e+01, 1.62e+01, 1.78e+00};
  Double_t cutiso_sumoetTight[9] = {
  7.75e+00, 5.45e+00, 5.67e+00, 5.97e+00, 3.17e+00, 3.86e+00, 6.06e+00, 5.31e+00, 1.05e+01};
  Double_t cutiso_sumoetlTight[9] = {
  7.56e+00, 5.08e+00, 5.77e+00, 5.74e+00, 2.37e+00, 3.32e+00, 4.97e+00, 5.46e+00, 3.82e+00};
  Double_t cutseeTight[9] = {
  1.16e-02, 1.07e-02, 1.08e-02, 3.49e-02, 2.89e-02, 3.08e-02, 9.87e-03, 3.37e-02, 4.40e-02};
  Double_t cutseelTight[9] = {
  1.27e-02, 1.08e-02, 1.13e-02, 4.19e-02, 2.81e-02, 3.02e-02, 9.76e-03, 4.28e-02, 2.98e-02};
  
  // SuperTight cuts
  Double_t cutdcotdistSuperTight[9] = {
  2.11e-02, 1.86e-02, 1.55e-02, 3.40e-02, 2.85e-02, 3.32e-02, 1.64e-02, 3.75e-02, 1.30e-04};
  Double_t cutdetainSuperTight[9] = {
  7.84e-03, 3.67e-03, 7.00e-03, 1.28e-02, 5.65e-03, 9.53e-03, 1.08e-02, 2.97e-02, 7.24e-03};
  Double_t cutdetainlSuperTight[9] = {
  7.61e-03, 3.28e-03, 6.57e-03, 1.03e-02, 5.05e-03, 8.55e-03, 1.07e-02, 2.94e-02, 4.10e-03};
  Double_t cutdphiinSuperTight[9] = {
  4.83e-02, 7.39e-02, 2.38e-01, 5.74e-02, 1.29e-01, 2.13e-01, 3.31e-01, 3.93e-01, 2.84e-01};
  Double_t cutdphiinlSuperTight[9] = {
  5.79e-02, 7.21e-02, 2.18e-01, 7.70e-02, 1.41e-01, 2.11e-01, 2.43e-01, 3.53e-01, 2.89e-01};
  Double_t cuteseedopcorSuperTight[9] = {
  7.32e-01, 9.77e-01, 9.83e-01, 8.55e-01, 4.31e-01, 7.35e-01, 4.18e-01, 9.99e-01, 5.89e-01};
  Double_t cutfmishitsSuperTight[9] = {
  3.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
  Double_t cuthoeSuperTight[9] = {
  9.19e-02, 4.11e-02, 1.42e-01, 3.35e-01, 3.82e-02, 1.41e-01, 4.29e-01, 4.01e-01, 3.99e-01};
  Double_t cuthoelSuperTight[9] = {
  7.51e-02, 3.81e-02, 1.41e-01, 3.32e-01, 3.10e-02, 1.43e-01, 2.35e-01, 3.80e-01, 1.32e-01};
  Double_t cutip_gsfSuperTight[9] = {
  1.42e-02, 2.66e-02, 1.06e-01, 3.38e-02, 3.23e-01, 1.07e-01, 7.74e-02, 2.32e-01, 7.80e-02};
  Double_t cutip_gsflSuperTight[9] = {
  1.15e-02, 2.72e-02, 8.41e-02, 2.49e-02, 4.17e-01, 1.02e-01, 7.90e-02, 1.69e-01, 4.79e-02};
  Double_t cutiso_sumSuperTight[9] = {
  8.95e+00, 8.18e+00, 8.75e+00, 7.47e+00, 5.43e+00, 5.87e+00, 8.16e+00, 1.02e+01, 1.78e+00};
  Double_t cutiso_sumoetSuperTight[9] = {
  6.45e+00, 5.14e+00, 4.99e+00, 5.21e+00, 2.65e+00, 3.12e+00, 4.52e+00, 4.72e+00, 3.68e+00};
  Double_t cutiso_sumoetlSuperTight[9] = {
  6.02e+00, 3.96e+00, 4.23e+00, 4.73e+00, 1.99e+00, 2.64e+00, 3.72e+00, 3.81e+00, 1.44e+00};
  Double_t cutseeSuperTight[9] = {
  1.09e-02, 1.05e-02, 1.05e-02, 3.24e-02, 2.81e-02, 2.95e-02, 9.77e-03, 2.75e-02, 2.95e-02};
  Double_t cutseelSuperTight[9] = {
  1.12e-02, 1.05e-02, 1.07e-02, 3.51e-02, 2.75e-02, 2.87e-02, 9.59e-03, 2.67e-02, 2.98e-02};

  // HyperTight1 cuts
  Double_t cutdcotdistHyperTight1[9] = {
  1.48e-02, 1.50e-02, 8.25e-03, 3.16e-02, 2.85e-02, 3.15e-02, 6.62e-03, 3.48e-02, 3.63e-06};
  Double_t cutdetainHyperTight1[9] = {
  6.51e-03, 3.51e-03, 5.53e-03, 9.16e-03, 5.30e-03, 8.28e-03, 1.08e-02, 2.97e-02, 7.24e-03};
  Double_t cutdetainlHyperTight1[9] = {
  6.05e-03, 3.23e-03, 4.93e-03, 8.01e-03, 4.93e-03, 7.91e-03, 1.03e-02, 2.94e-02, 4.10e-03};
  Double_t cutdphiinHyperTight1[9] = {
  4.83e-02, 4.91e-02, 2.30e-01, 3.48e-02, 7.44e-02, 2.04e-01, 9.95e-02, 3.93e-01, 2.84e-01};
  Double_t cutdphiinlHyperTight1[9] = {
  4.74e-02, 4.51e-02, 2.18e-01, 2.99e-02, 7.37e-02, 2.11e-01, 9.99e-02, 3.53e-01, 2.89e-01};
  Double_t cuteseedopcorHyperTight1[9] = {
  7.72e-01, 9.90e-01, 1.01e+00, 8.55e-01, 9.11e-01, 7.72e-01, 9.17e-01, 1.06e+00, 7.63e-01};
  Double_t cutfmishitsHyperTight1[9] = {
  3.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
  Double_t cuthoeHyperTight1[9] = {
  6.17e-02, 3.70e-02, 1.41e-01, 2.91e-01, 3.82e-02, 1.34e-01, 4.19e-01, 3.87e-01, 3.93e-01};
  Double_t cuthoelHyperTight1[9] = {
  4.43e-02, 3.57e-02, 1.41e-01, 2.81e-01, 3.07e-02, 1.28e-01, 2.27e-01, 3.80e-01, 1.32e-01};
  Double_t cutip_gsfHyperTight1[9] = {
  1.21e-02, 1.76e-02, 6.01e-02, 2.96e-02, 1.74e-01, 9.70e-02, 7.74e-02, 1.33e-01, 7.80e-02};
  Double_t cutip_gsflHyperTight1[9] = {
  1.01e-02, 1.56e-02, 6.87e-02, 2.13e-02, 1.25e-01, 8.16e-02, 7.90e-02, 1.30e-01, 4.79e-02};
  Double_t cutiso_sumHyperTight1[9] = {
  7.92e+00, 6.85e+00, 7.87e+00, 6.77e+00, 4.47e+00, 5.28e+00, 6.57e+00, 1.02e+01, 1.78e+00};
  Double_t cutiso_sumoetHyperTight1[9] = {
  5.20e+00, 3.93e+00, 3.88e+00, 4.10e+00, 2.40e+00, 2.43e+00, 3.49e+00, 3.94e+00, 3.01e+00};
  Double_t cutiso_sumoetlHyperTight1[9] = {
  4.18e+00, 3.12e+00, 3.44e+00, 3.25e+00, 1.77e+00, 2.06e+00, 2.83e+00, 3.12e+00, 1.43e+00};
  Double_t cutseeHyperTight1[9] = {
  1.05e-02, 1.04e-02, 1.01e-02, 3.24e-02, 2.80e-02, 2.85e-02, 9.67e-03, 2.61e-02, 2.95e-02};
  Double_t cutseelHyperTight1[9] = {
  1.04e-02, 1.03e-02, 1.01e-02, 3.04e-02, 2.74e-02, 2.78e-02, 9.58e-03, 2.54e-02, 2.83e-02};

  // HyperTight2 cuts
  Double_t cutdcotdistHyperTight2[9] = {
  1.15e-02, 1.07e-02, 4.01e-03, 2.97e-02, 2.85e-02, 3.10e-02, 9.34e-04, 3.40e-02, 2.82e-07};
  Double_t cutdetainHyperTight2[9] = {
  5.29e-03, 2.56e-03, 4.89e-03, 7.89e-03, 5.30e-03, 7.37e-03, 8.91e-03, 9.36e-03, 5.94e-03};
  Double_t cutdetainlHyperTight2[9] = {
  4.48e-03, 2.59e-03, 4.42e-03, 6.54e-03, 4.93e-03, 6.98e-03, 8.49e-03, 9.06e-03, -4.81e-03};
  Double_t cutdphiinHyperTight2[9] = {
  2.41e-02, 3.83e-02, 1.48e-01, 2.91e-02, 3.15e-02, 1.57e-01, 8.90e-02, 1.02e-01, 2.81e-01};
  Double_t cutdphiinlHyperTight2[9] = {
  2.13e-02, 3.79e-02, 1.25e-01, 2.24e-02, 3.69e-02, 1.64e-01, 9.99e-02, 9.23e-02, 2.37e-01};
  Double_t cuteseedopcorHyperTight2[9] = {
  1.03e+00, 9.95e-01, 1.03e+00, 1.01e+00, 9.46e-01, 9.03e-01, 9.97e-01, 1.14e+00, 8.00e-01};
  Double_t cutfmishitsHyperTight2[9] = {
  1.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01};
  Double_t cuthoeHyperTight2[9] = {
  4.94e-02, 3.45e-02, 1.40e-01, 2.02e-01, 3.82e-02, 1.19e-01, 1.23e-01, 3.82e-01, 2.50e-01};
  Double_t cuthoelHyperTight2[9] = {
  4.04e-02, 3.42e-02, 1.31e-01, 1.85e-01, 3.01e-02, 1.27e-01, 2.27e-01, 3.80e-01, 1.32e-01};
  Double_t cutip_gsfHyperTight2[9] = {
  1.14e-02, 1.38e-02, 5.29e-02, 1.87e-02, 1.31e-01, 8.63e-02, 7.74e-02, 1.04e-01, 2.42e-02};
  Double_t cutip_gsflHyperTight2[9] = {
  9.83e-03, 1.35e-02, 4.27e-02, 1.72e-02, 1.25e-01, 7.92e-02, 7.90e-02, 1.30e-01, 3.40e-02};
  Double_t cutiso_sumHyperTight2[9] = {
  6.40e+00, 5.77e+00, 6.54e+00, 5.22e+00, 3.86e+00, 4.63e+00, 6.31e+00, 1.02e+01, 1.78e+00};
  Double_t cutiso_sumoetHyperTight2[9] = {
  4.03e+00, 3.03e+00, 3.24e+00, 3.13e+00, 2.05e+00, 2.01e+00, 2.99e+00, 3.44e+00, 2.76e+00};
  Double_t cutiso_sumoetlHyperTight2[9] = {
  3.08e+00, 2.31e+00, 2.84e+00, 2.53e+00, 1.65e+00, 1.72e+00, 2.34e+00, 3.11e+00, 1.35e+00};
  Double_t cutseeHyperTight2[9] = {
  1.03e-02, 1.03e-02, 9.88e-03, 3.03e-02, 2.79e-02, 2.79e-02, 9.67e-03, 2.52e-02, 2.58e-02};
  Double_t cutseelHyperTight2[9] = {
  1.02e-02, 1.02e-02, 9.80e-03, 2.90e-02, 2.74e-02, 2.75e-02, 9.58e-03, 2.49e-02, 2.50e-02};

  // HyperTight3 cuts
  Double_t cutdcotdistHyperTight3[9] = {
  9.63e-03, 5.11e-03, 1.95e-04, 2.97e-02, 2.85e-02, 2.18e-02, 2.61e-05, 2.57e-02, 2.82e-07};
  Double_t cutdetainHyperTight3[9] = {
  4.86e-03, 2.29e-03, 4.40e-03, 7.79e-03, 4.07e-03, 6.33e-03, 7.70e-03, 7.93e-03, 5.94e-03};
  Double_t cutdetainlHyperTight3[9] = {
  4.48e-03, 2.30e-03, 4.14e-03, 6.04e-03, 3.87e-03, 6.09e-03, 7.97e-03, 8.04e-03, -4.81e-03};
  Double_t cutdphiinHyperTight3[9] = {
  2.41e-02, 2.88e-02, 7.39e-02, 2.91e-02, 1.91e-02, 1.14e-01, 3.61e-02, 8.92e-02, 2.81e-01};
  Double_t cutdphiinlHyperTight3[9] = {
  1.95e-02, 3.42e-02, 8.06e-02, 2.22e-02, 2.26e-02, 9.73e-02, 4.51e-02, 9.23e-02, 2.37e-01};
  Double_t cuteseedopcorHyperTight3[9] = {
  1.07e+00, 1.01e+00, 1.08e+00, 1.01e+00, 9.69e-01, 9.10e-01, 1.04e+00, 1.20e+00, 8.00e-01};
  Double_t cutfmishitsHyperTight3[9] = {
  5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01};
  Double_t cuthoeHyperTight3[9] = {
  3.52e-02, 3.45e-02, 1.33e-01, 1.88e-01, 2.72e-02, 1.19e-01, 9.28e-02, 2.46e-01, 2.50e-01};
  Double_t cuthoelHyperTight3[9] = {
  4.04e-02, 3.40e-02, 1.31e-01, 1.84e-01, 2.64e-02, 1.18e-01, 9.76e-02, 2.53e-01, 1.32e-01};
  Double_t cutip_gsfHyperTight3[9] = {
  1.14e-02, 1.26e-02, 3.79e-02, 1.68e-02, 1.21e-01, 5.29e-02, 7.74e-02, 3.35e-02, 2.42e-02};
  Double_t cutip_gsflHyperTight3[9] = {
  9.83e-03, 1.18e-02, 3.59e-02, 1.56e-02, 1.20e-01, 5.36e-02, 7.90e-02, 2.88e-02, 3.40e-02};
  Double_t cutiso_sumHyperTight3[9] = {
  5.40e+00, 5.41e+00, 5.88e+00, 4.32e+00, 3.86e+00, 4.33e+00, 5.87e+00, 9.05e+00, 1.78e+00};
  Double_t cutiso_sumoetHyperTight3[9] = {
  3.03e+00, 2.50e+00, 2.58e+00, 2.44e+00, 1.91e+00, 1.76e+00, 2.92e+00, 3.13e+00, 2.76e+00};
  Double_t cutiso_sumoetlHyperTight3[9] = {
  2.36e+00, 2.02e+00, 2.29e+00, 1.89e+00, 1.65e+00, 1.69e+00, 2.03e+00, 2.79e+00, 1.35e+00};
  Double_t cutseeHyperTight3[9] = {
  1.03e-02, 1.01e-02, 9.84e-03, 2.89e-02, 2.74e-02, 2.73e-02, 9.47e-03, 2.44e-02, 2.58e-02};
  Double_t cutseelHyperTight3[9] = {
  1.02e-02, 1.00e-02, 9.73e-03, 2.79e-02, 2.73e-02, 2.69e-02, 9.40e-03, 2.46e-02, 2.50e-02};

  // HyperTight4 cuts
  Double_t cutdcotdistHyperTight4[9] = {
  2.70e-04, 1.43e-04, 1.95e-04, 2.64e-03, 2.82e-02, 1.64e-02, 2.61e-05, 2.57e-02, 2.82e-07};
  Double_t cutdetainHyperTight4[9] = {
  2.44e-03, 1.67e-03, 2.26e-03, 3.43e-03, 3.51e-03, 3.52e-03, 2.98e-03, 4.79e-03, 5.94e-03};
  Double_t cutdetainlHyperTight4[9] = {
  2.34e-03, 1.29e-03, 2.30e-03, 3.30e-03, 3.61e-03, 3.84e-03, 2.53e-03, 3.66e-03, -4.81e-03};
  Double_t cutdphiinHyperTight4[9] = {
  8.44e-03, 5.21e-03, 2.18e-02, 1.39e-02, 7.82e-03, 1.52e-02, 2.59e-02, 3.87e-02, 2.81e-01};
  Double_t cutdphiinlHyperTight4[9] = {
  5.77e-03, 3.20e-03, 2.85e-02, 2.22e-02, 7.00e-03, 1.84e-02, 2.91e-02, 4.40e-02, 2.37e-01};
  Double_t cuteseedopcorHyperTight4[9] = {
  1.15e+00, 1.01e+00, 1.21e+00, 1.07e+00, 9.69e-01, 9.10e-01, 1.08e+00, 1.36e+00, 8.00e-01};
  Double_t cutfmishitsHyperTight4[9] = {
  5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01};
  Double_t cuthoeHyperTight4[9] = {
  2.39e-02, 2.68e-02, 2.12e-02, 1.03e-01, 9.92e-03, 7.07e-02, 7.12e-02, 1.48e-01, 2.50e-01};
  Double_t cuthoelHyperTight4[9] = {
  2.87e-02, 1.94e-02, 2.16e-02, 5.68e-02, 1.35e-02, 4.04e-02, 7.98e-02, 1.50e-01, 1.32e-01};
  Double_t cutip_gsfHyperTight4[9] = {
  7.61e-03, 5.22e-03, 3.79e-02, 1.02e-02, 4.62e-02, 1.82e-02, 7.74e-02, 3.35e-02, 2.42e-02};
  Double_t cutip_gsflHyperTight4[9] = {
  7.81e-03, 4.25e-03, 3.08e-02, 1.04e-02, 2.35e-02, 2.45e-02, 7.90e-02, 2.88e-02, 3.40e-02};
  Double_t cutiso_sumHyperTight4[9] = {
  5.40e+00, 5.41e+00, 5.88e+00, 4.32e+00, 3.86e+00, 4.33e+00, 5.86e+00, 9.05e+00, 1.78e+00};
  Double_t cutiso_sumoetHyperTight4[9] = {
  2.53e+00, 2.10e+00, 1.87e+00, 1.84e+00, 1.79e+00, 1.61e+00, 2.53e+00, 1.98e+00, 2.76e+00};
  Double_t cutiso_sumoetlHyperTight4[9] = {
  2.28e+00, 2.02e+00, 2.04e+00, 1.69e+00, 1.65e+00, 1.61e+00, 2.03e+00, 1.82e+00, 1.35e+00};
  Double_t cutseeHyperTight4[9] = {
  9.99e-03, 9.61e-03, 9.65e-03, 2.75e-02, 2.61e-02, 2.64e-02, 9.18e-03, 2.44e-02, 2.58e-02};
  Double_t cutseelHyperTight4[9] = {
  9.66e-03, 9.69e-03, 9.58e-03, 2.73e-02, 2.66e-02, 2.66e-02, 8.64e-03, 2.46e-02, 2.50e-02};

  Double_t cutdcotdist[9];
  Double_t cutdetain[9];
  Double_t cutdetainl[9];
  Double_t cutdphiin[9];
  Double_t cutdphiinl[9];
  Double_t cuteseedopcor[9];
  Double_t cutfmishits[9];
  Double_t cuthoe[9];
  Double_t cuthoel[9];
  Double_t cutip_gsf[9];
  Double_t cutip_gsfl[9];
  Double_t cutiso_sum[9];
  Double_t cutiso_sumoet[9];
  Double_t cutiso_sumoetl[9];
  Double_t cutsee[9];
  Double_t cutseel[9];
  if	 (typeCuts == 0) {
    memcpy(cutdcotdist   ,cutdcotdistMedium   ,sizeof(cutdcotdistMedium));
    memcpy(cutdetain     ,cutdetainMedium     ,sizeof(cutdetainMedium));
    memcpy(cutdetainl    ,cutdetainlMedium    ,sizeof(cutdetainlMedium));
    memcpy(cutdphiin     ,cutdphiinMedium     ,sizeof(cutdphiinMedium));
    memcpy(cutdphiinl    ,cutdphiinlMedium    ,sizeof(cutdphiinlMedium));
    memcpy(cuteseedopcor ,cuteseedopcorMedium ,sizeof(cuteseedopcorMedium));
    memcpy(cutfmishits   ,cutfmishitsMedium   ,sizeof(cutfmishitsMedium));
    memcpy(cuthoe        ,cuthoeMedium	      ,sizeof(cuthoeMedium));
    memcpy(cuthoel       ,cuthoelMedium	      ,sizeof(cuthoelMedium));
    memcpy(cutip_gsf     ,cutip_gsfMedium     ,sizeof(cutip_gsfMedium));
    memcpy(cutip_gsfl    ,cutip_gsflMedium    ,sizeof(cutip_gsflMedium));
    memcpy(cutiso_sum    ,cutiso_sumMedium    ,sizeof(cutiso_sumMedium));
    memcpy(cutiso_sumoet ,cutiso_sumoetMedium ,sizeof(cutiso_sumoetMedium));
    memcpy(cutiso_sumoetl,cutiso_sumoetlMedium,sizeof(cutiso_sumoetlMedium));
    memcpy(cutsee        ,cutseeMedium	      ,sizeof(cutseeMedium));
    memcpy(cutseel       ,cutseelMedium	      ,sizeof(cutseelMedium));
  }
  else if(typeCuts == 1) {
    memcpy(cutdcotdist   ,cutdcotdistTight   ,sizeof(cutdcotdistTight));
    memcpy(cutdetain     ,cutdetainTight     ,sizeof(cutdetainTight));
    memcpy(cutdetainl    ,cutdetainlTight    ,sizeof(cutdetainlTight));
    memcpy(cutdphiin     ,cutdphiinTight     ,sizeof(cutdphiinTight));
    memcpy(cutdphiinl    ,cutdphiinlTight    ,sizeof(cutdphiinlTight));
    memcpy(cuteseedopcor ,cuteseedopcorTight ,sizeof(cuteseedopcorTight));
    memcpy(cutfmishits   ,cutfmishitsTight   ,sizeof(cutfmishitsTight));
    memcpy(cuthoe        ,cuthoeTight	     ,sizeof(cuthoeTight));
    memcpy(cuthoel       ,cuthoelTight	     ,sizeof(cuthoelTight));
    memcpy(cutip_gsf     ,cutip_gsfTight     ,sizeof(cutip_gsfTight));
    memcpy(cutip_gsfl    ,cutip_gsflTight    ,sizeof(cutip_gsflTight));
    memcpy(cutiso_sum    ,cutiso_sumTight    ,sizeof(cutiso_sumTight));
    memcpy(cutiso_sumoet ,cutiso_sumoetTight ,sizeof(cutiso_sumoetTight));
    memcpy(cutiso_sumoetl,cutiso_sumoetlTight,sizeof(cutiso_sumoetlTight));
    memcpy(cutsee        ,cutseeTight	     ,sizeof(cutseeTight));
    memcpy(cutseel       ,cutseelTight	     ,sizeof(cutseelTight));
  }
  else if(typeCuts == 2) {
    memcpy(cutdcotdist   ,cutdcotdistSuperTight   ,sizeof(cutdcotdistSuperTight));
    memcpy(cutdetain     ,cutdetainSuperTight     ,sizeof(cutdetainSuperTight));
    memcpy(cutdetainl    ,cutdetainlSuperTight    ,sizeof(cutdetainlSuperTight));
    memcpy(cutdphiin     ,cutdphiinSuperTight     ,sizeof(cutdphiinSuperTight));
    memcpy(cutdphiinl    ,cutdphiinlSuperTight    ,sizeof(cutdphiinlSuperTight));
    memcpy(cuteseedopcor ,cuteseedopcorSuperTight ,sizeof(cuteseedopcorSuperTight));
    memcpy(cutfmishits   ,cutfmishitsSuperTight   ,sizeof(cutfmishitsSuperTight));
    memcpy(cuthoe        ,cuthoeSuperTight	  ,sizeof(cuthoeSuperTight));
    memcpy(cuthoel       ,cuthoelSuperTight	  ,sizeof(cuthoelSuperTight));
    memcpy(cutip_gsf     ,cutip_gsfSuperTight     ,sizeof(cutip_gsfSuperTight));
    memcpy(cutip_gsfl    ,cutip_gsflSuperTight    ,sizeof(cutip_gsflSuperTight));
    memcpy(cutiso_sum    ,cutiso_sumSuperTight    ,sizeof(cutiso_sumSuperTight));
    memcpy(cutiso_sumoet ,cutiso_sumoetSuperTight ,sizeof(cutiso_sumoetSuperTight));
    memcpy(cutiso_sumoetl,cutiso_sumoetlSuperTight,sizeof(cutiso_sumoetlSuperTight));
    memcpy(cutsee        ,cutseeSuperTight	  ,sizeof(cutseeSuperTight));
    memcpy(cutseel       ,cutseelSuperTight	  ,sizeof(cutseelSuperTight));
  }
  else if(typeCuts == 3) {
    memcpy(cutdcotdist   ,cutdcotdistHyperTight1   ,sizeof(cutdcotdistHyperTight1));
    memcpy(cutdetain     ,cutdetainHyperTight1     ,sizeof(cutdetainHyperTight1));
    memcpy(cutdetainl    ,cutdetainlHyperTight1    ,sizeof(cutdetainlHyperTight1));
    memcpy(cutdphiin     ,cutdphiinHyperTight1     ,sizeof(cutdphiinHyperTight1));
    memcpy(cutdphiinl    ,cutdphiinlHyperTight1    ,sizeof(cutdphiinlHyperTight1));
    memcpy(cuteseedopcor ,cuteseedopcorHyperTight1 ,sizeof(cuteseedopcorHyperTight1));
    memcpy(cutfmishits   ,cutfmishitsHyperTight1   ,sizeof(cutfmishitsHyperTight1));
    memcpy(cuthoe        ,cuthoeHyperTight1	  ,sizeof(cuthoeHyperTight1));
    memcpy(cuthoel       ,cuthoelHyperTight1	  ,sizeof(cuthoelHyperTight1));
    memcpy(cutip_gsf     ,cutip_gsfHyperTight1     ,sizeof(cutip_gsfHyperTight1));
    memcpy(cutip_gsfl    ,cutip_gsflHyperTight1    ,sizeof(cutip_gsflHyperTight1));
    memcpy(cutiso_sum    ,cutiso_sumHyperTight1    ,sizeof(cutiso_sumHyperTight1));
    memcpy(cutiso_sumoet ,cutiso_sumoetHyperTight1 ,sizeof(cutiso_sumoetHyperTight1));
    memcpy(cutiso_sumoetl,cutiso_sumoetlHyperTight1,sizeof(cutiso_sumoetlHyperTight1));
    memcpy(cutsee        ,cutseeHyperTight1	  ,sizeof(cutseeHyperTight1));
    memcpy(cutseel       ,cutseelHyperTight1	  ,sizeof(cutseelHyperTight1));
  }
  else if(typeCuts == 4) {
    memcpy(cutdcotdist   ,cutdcotdistHyperTight2   ,sizeof(cutdcotdistHyperTight2));
    memcpy(cutdetain     ,cutdetainHyperTight2     ,sizeof(cutdetainHyperTight2));
    memcpy(cutdetainl    ,cutdetainlHyperTight2    ,sizeof(cutdetainlHyperTight2));
    memcpy(cutdphiin     ,cutdphiinHyperTight2     ,sizeof(cutdphiinHyperTight2));
    memcpy(cutdphiinl    ,cutdphiinlHyperTight2    ,sizeof(cutdphiinlHyperTight2));
    memcpy(cuteseedopcor ,cuteseedopcorHyperTight2 ,sizeof(cuteseedopcorHyperTight2));
    memcpy(cutfmishits   ,cutfmishitsHyperTight2   ,sizeof(cutfmishitsHyperTight2));
    memcpy(cuthoe        ,cuthoeHyperTight2	  ,sizeof(cuthoeHyperTight2));
    memcpy(cuthoel       ,cuthoelHyperTight2	  ,sizeof(cuthoelHyperTight2));
    memcpy(cutip_gsf     ,cutip_gsfHyperTight2     ,sizeof(cutip_gsfHyperTight2));
    memcpy(cutip_gsfl    ,cutip_gsflHyperTight2    ,sizeof(cutip_gsflHyperTight2));
    memcpy(cutiso_sum    ,cutiso_sumHyperTight2    ,sizeof(cutiso_sumHyperTight2));
    memcpy(cutiso_sumoet ,cutiso_sumoetHyperTight2 ,sizeof(cutiso_sumoetHyperTight2));
    memcpy(cutiso_sumoetl,cutiso_sumoetlHyperTight2,sizeof(cutiso_sumoetlHyperTight2));
    memcpy(cutsee        ,cutseeHyperTight2	  ,sizeof(cutseeHyperTight2));
    memcpy(cutseel       ,cutseelHyperTight2	  ,sizeof(cutseelHyperTight2));
  }
  else if(typeCuts == 5) {
    memcpy(cutdcotdist   ,cutdcotdistHyperTight3   ,sizeof(cutdcotdistHyperTight3));
    memcpy(cutdetain     ,cutdetainHyperTight3     ,sizeof(cutdetainHyperTight3));
    memcpy(cutdetainl    ,cutdetainlHyperTight3    ,sizeof(cutdetainlHyperTight3));
    memcpy(cutdphiin     ,cutdphiinHyperTight3     ,sizeof(cutdphiinHyperTight3));
    memcpy(cutdphiinl    ,cutdphiinlHyperTight3    ,sizeof(cutdphiinlHyperTight3));
    memcpy(cuteseedopcor ,cuteseedopcorHyperTight3 ,sizeof(cuteseedopcorHyperTight3));
    memcpy(cutfmishits   ,cutfmishitsHyperTight3   ,sizeof(cutfmishitsHyperTight3));
    memcpy(cuthoe        ,cuthoeHyperTight3	  ,sizeof(cuthoeHyperTight3));
    memcpy(cuthoel       ,cuthoelHyperTight3	  ,sizeof(cuthoelHyperTight3));
    memcpy(cutip_gsf     ,cutip_gsfHyperTight3     ,sizeof(cutip_gsfHyperTight3));
    memcpy(cutip_gsfl    ,cutip_gsflHyperTight3    ,sizeof(cutip_gsflHyperTight3));
    memcpy(cutiso_sum    ,cutiso_sumHyperTight3    ,sizeof(cutiso_sumHyperTight3));
    memcpy(cutiso_sumoet ,cutiso_sumoetHyperTight3 ,sizeof(cutiso_sumoetHyperTight3));
    memcpy(cutiso_sumoetl,cutiso_sumoetlHyperTight3,sizeof(cutiso_sumoetlHyperTight3));
    memcpy(cutsee        ,cutseeHyperTight3	  ,sizeof(cutseeHyperTight3));
    memcpy(cutseel       ,cutseelHyperTight3	  ,sizeof(cutseelHyperTight3));
  }
  else if(typeCuts == 6) {
    memcpy(cutdcotdist   ,cutdcotdistHyperTight4   ,sizeof(cutdcotdistHyperTight4));
    memcpy(cutdetain     ,cutdetainHyperTight4     ,sizeof(cutdetainHyperTight4));
    memcpy(cutdetainl    ,cutdetainlHyperTight4    ,sizeof(cutdetainlHyperTight4));
    memcpy(cutdphiin     ,cutdphiinHyperTight4     ,sizeof(cutdphiinHyperTight4));
    memcpy(cutdphiinl    ,cutdphiinlHyperTight4    ,sizeof(cutdphiinlHyperTight4));
    memcpy(cuteseedopcor ,cuteseedopcorHyperTight4 ,sizeof(cuteseedopcorHyperTight4));
    memcpy(cutfmishits   ,cutfmishitsHyperTight4   ,sizeof(cutfmishitsHyperTight4));
    memcpy(cuthoe        ,cuthoeHyperTight4	  ,sizeof(cuthoeHyperTight4));
    memcpy(cuthoel       ,cuthoelHyperTight4	  ,sizeof(cuthoelHyperTight4));
    memcpy(cutip_gsf     ,cutip_gsfHyperTight4     ,sizeof(cutip_gsfHyperTight4));
    memcpy(cutip_gsfl    ,cutip_gsflHyperTight4    ,sizeof(cutip_gsflHyperTight4));
    memcpy(cutiso_sum    ,cutiso_sumHyperTight4    ,sizeof(cutiso_sumHyperTight4));
    memcpy(cutiso_sumoet ,cutiso_sumoetHyperTight4 ,sizeof(cutiso_sumoetHyperTight4));
    memcpy(cutiso_sumoetl,cutiso_sumoetlHyperTight4,sizeof(cutiso_sumoetlHyperTight4));
    memcpy(cutsee        ,cutseeHyperTight4	  ,sizeof(cutseeHyperTight4));
    memcpy(cutseel       ,cutseelHyperTight4	  ,sizeof(cutseelHyperTight4));
  }
  else {
    return 0;
  }
  
  const int ncuts = 10;
  std::vector<bool> cut_results(ncuts, false);
  
  float iso_sum = tkIso + ecalIso + hcalIso;
  if(fabs(scEta)>1.5) 
    iso_sum += (fabs(scEta)-1.5)*1.09;
  
//   float iso_sumoet = iso_sum*(40./scEt);
  
  float eseedopincor = eSeedOverPin + fBrem;
  if(fBrem < 0)
    eseedopincor = eSeedOverPin;

  //float dist = (TMath::Abs(partnerDist)      == -9999.? 9999:TMath::Abs(partnerDist));
  //float dcot = (TMath::Abs(partnerDeltaCot) == -9999.? 9999:TMath::Abs(partnerDeltaCot));

//   float dcotdistcomb = ((0.04 - std::max(dist, dcot)) > 0?(0.04 - std::max(dist, dcot)):0);

//   Double_t ip = 99999;
//   ip = d0;

  if (debug) {
    cout << "eseedopincor = " << eseedopincor << " >? " <<  cuteseedopcor[cat] << endl;
  }

  for (int cut=0; cut<ncuts; cut++) {
    switch (cut) {
    case 0:
      cut_results[cut] = compute_cut(fabs(deltaEtaIn), scEt, cutdetainl[cat], cutdetain[cat]);
      break;
    case 1:
      cut_results[cut] = compute_cut(fabs(deltaPhiIn), scEt, cutdphiinl[cat], cutdphiin[cat]);
      break;
    case 2:
      cut_results[cut] = (eseedopincor > cuteseedopcor[cat]);
      break;
    case 3:
      cut_results[cut] = compute_cut(hOverE, scEt, cuthoel[cat], cuthoe[cat]);
      break;
    case 4:
      cut_results[cut] = compute_cut(sigmaee, scEt, cutseel[cat], cutsee[cat]);
      break;
    case 5:
      // cut_results[cut] = compute_cut(iso_sumoet, scEt, cutiso_sumoetl[cat], cutiso_sumoet[cat]);
      cut_results[cut] = kTRUE;
      break;
    case 6:
      //cut_results[cut] = (iso_sum < cutiso_sum[cat]);
      cut_results[cut] = (tkIso/pt < 0.7);
      break;
    case 7:
      //cut_results[cut] = compute_cut(fabs(ip), scEt, cutip_gsfl[cat], cutip_gsf[cat]);
      //cut_results[cut] = (ip3dSig < 4);
      cut_results[cut] = kTRUE;
      break;
    case 8:
      cut_results[cut] = (mishits < cutfmishits[cat]);
      break;
    case 9:
//       cut_results[cut] = (dcotdistcomb < cutdcotdist[cat]);
      cut_results[cut] = kTRUE;
      break;
    }
  }
    
  if (debug) {
    cout << cut_results[0]  << " " 
         <<  cut_results[1]  << " " 
         <<  cut_results[2]  << " " 
         <<  cut_results[3]  << " " 
         <<  cut_results[4] << " "
         <<  cut_results[6]<< " " 
         <<  cut_results[7]<< " "   
         <<  cut_results[8]
         << endl;
  }

  Bool_t result = kFALSE;
  if (cut_results[0] && cut_results[1] && cut_results[2] && cut_results[3] && cut_results[4]
      &&
      cut_results[6]
      && 
      cut_results[7]
      && 
      cut_results[8]
    ) {
    result = kTRUE;
  }

  if (debug) {
    cout << "result = " << result << endl;
  }

  return result;
}


//--------------------------------------------------------------------------------------------------
Int_t PassCiCID(const higgsana::TElectron *ele,
                const Int_t typeCuts, Bool_t debug ) {
  return PassCiCID(ele->pt, ele->scEt, ele->scEta, ele->isEcalDriven, ele->isEB, 
                   ele->fBrem, ele->EOverP, ele->HoverE, ele->sigiEtaiEta, ele->deltaPhiIn, ele->deltaEtaIn, 
                   ele->ESeedClusterOverPIn, ele->nExpHitsInner, 
                   ele->trkIso03, ele->emIso04, ele->hadIso04, 
                   ele->partnerDist, ele->partnerDeltaCot, 
                   ele->d0, ele->ip3dSig, 
                   typeCuts, debug);

//   return PassCiCID(ele->pt, ele->scEt, ele->scEta, 1,  Bool_t(fabs(ele->eta) < 1.479), 
//                    ele->fBrem, ele->EOverP, ele->HoverE, ele->sigiEtaiEta, ele->deltaPhiIn, ele->deltaEtaIn, 
//                    ele->ESeedClusterOverPIn, 0, 
//                    0,0,0,
//                    0,0,
//                    ele->d0, ele->ip3dSig, 
//                    typeCuts, debug);
}

//--------------------------------------------------------------------------------------------------
Int_t PassCiCID(const higgsana::TElectron *ele,
                const Int_t typeCuts ) {
  return PassCiCID(ele->pt, ele->scEt, ele->scEta, ele->isEcalDriven, ele->isEB, 
                   ele->fBrem, ele->EOverP, ele->HoverE, ele->sigiEtaiEta, ele->deltaPhiIn, ele->deltaEtaIn, 
                   ele->ESeedClusterOverPIn, ele->nExpHitsInner, 
                   ele->trkIso03, ele->emIso04, ele->hadIso04, 
                   ele->partnerDist, ele->partnerDeltaCot, 
                   ele->d0, ele->ip3dSig, 
                   typeCuts, kFALSE);
//   return PassCiCID(ele->pt, ele->scEt, ele->scEta, 1,  Bool_t(fabs(ele->eta) < 1.479), 
//                    ele->fBrem, ele->EOverP, ele->HoverE, ele->sigiEtaiEta, ele->deltaPhiIn, ele->deltaEtaIn, 
//                    ele->ESeedClusterOverPIn, 0, 
//                    0,0,0,
//                    0,0,
//                    ele->d0, ele->ip3dSig, 
//                    typeCuts, kFALSE);
}

//--------------------------------------------------------------------------------------------------
Int_t PassEleDetIso025(const higgsana::TElectron *ele, Double_t rho) {

  Bool_t pass = kFALSE;
  Double_t DetIso03PUCorrection = 0;
  if (fabs(ele->scEta) < 1.5) DetIso03PUCorrection = rho * (0.078 + 0.026);  
  else DetIso03PUCorrection = rho * (0.046 + 0.072);
  Double_t iso = (ele->trkIso03 + ele->emIso03 + ele->hadIso03 - DetIso03PUCorrection) / ele->pt;
  if (iso < 0.25 && ele->trkIso03/ele->pt < 0.7) pass = kTRUE;
  
  return pass;

}



//--------------------------------------------------------------------------------------------------
Bool_t passMuonID_HZZ2011(const higgsana::TMuon *muon)
{
  
  if(!(muon->typeBits & kGlobal))  return kFALSE;
  if(muon->nTkHits	      < 11 )    return kFALSE;

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t passMuonIP_HZZ2011(const higgsana::TMuon *muon)
{
  if(fabs(muon->ip3dSig)  > 4) return kFALSE;    

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Int_t PassMuDetIso025(const higgsana::TMuon *mu, Double_t rho) {

  Bool_t pass = kFALSE;
  Double_t DetIso03PUCorrection = 0;
  if (fabs(mu->eta) < 1.5) DetIso03PUCorrection = rho * (0.087 + 0.042);  
  else DetIso03PUCorrection = rho * (0.049 + 0.059);
  Double_t iso = (mu->trkIso03 + mu->emIso03 + mu->hadIso03 - DetIso03PUCorrection) / mu->pt;
  if (iso < 0.25 && mu->trkIso03/mu->pt < 0.7) pass = kTRUE;
  
  return pass;

}

//--------------------------------------------------------------------------------------------------
Bool_t passEleWP95ID(const higgsana::TElectron *electron)
{

  if(fabs(electron->d0) > 0.04) return kFALSE;
  if(fabs(electron->dz) > 0.2)  return kFALSE;
  
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {
     
    if(electron->pt>20) {
      if(electron->sigiEtaiEta	    > 0.01)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.007) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.8)  return kFALSE;
      if(electron->HoverE	    > 0.15)  return kFALSE;
    
    } else {
      if(electron->sigiEtaiEta	    > 0.012)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.007) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.8)  return kFALSE;
      if(electron->HoverE	    > 0.15) return kFALSE;    
    }
  
  } else {
     
    if(electron->pt>20) {
      if(electron->sigiEtaiEta	    > 0.03)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.010) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.8)  return kFALSE;
    } else {
      if(electron->sigiEtaiEta	    > 0.032)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.010) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.8)  return kFALSE;
    }
  }

  return kTRUE;
}

Bool_t passMuonIsoMVASameAsDetIso025( Double_t fPt, Double_t fEta, Double_t MVAValue) {

  Int_t subdet = 0;
  if (fabs(fEta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (fPt > 10.0) ptBin = 1;
  if (fPt > 20.0) ptBin = 2;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;
  if (subdet == 0 && ptBin == 2) MVABin = 4;
  if (subdet == 1 && ptBin == 2) MVABin = 5;

  Double_t MVACut = -999;
  if (MVABin == 0) MVACut = -0.4050;
  if (MVABin == 1) MVACut = -0.6022;
  if (MVABin == 2) MVACut = 0.5402;
  if (MVABin == 3) MVACut = 0.4338;
  if (MVABin == 4) MVACut = 0.7962;
  if (MVABin == 5) MVACut = 0.7422;

  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}

Bool_t passMuonIsoMVASameAsDetIso015( Double_t fPt, Double_t fEta, Double_t MVAValue) {

  Int_t subdet = 0;
  if (fabs(fEta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (fPt > 10.0) ptBin = 1;
  if (fPt > 20.0) ptBin = 2;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;
  if (subdet == 0 && ptBin == 2) MVABin = 4;
  if (subdet == 1 && ptBin == 2) MVABin = 5;

  Double_t MVACut = -999;
  if (MVABin == 0) MVACut = -0.0854;
  if (MVABin == 1) MVACut = -0.4158;
  if (MVABin == 2) MVACut = 0.7958;
  if (MVABin == 3) MVACut = 0.7058;
  if (MVABin == 4) MVACut = 0.9602;
  if (MVABin == 5) MVACut = 0.9482;

  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}


Bool_t passEleIsoMVASameBkgAsDetIso025( Double_t fPt, Double_t fEta, Double_t MVAValue) {

  Int_t subdet = 0;
  if (fabs(fEta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (fPt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;

  Double_t MVACut = -999;
  if (MVABin == 0) MVACut = -0.623;
  if (MVABin == 1) MVACut = -0.740;
  if (MVABin == 2) MVACut = -0.055;
  if (MVABin == 3) MVACut = 0.197;

  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}


Bool_t passEleIsoMVASameAsBkgDetIso015( Double_t fPt, Double_t fEta, Double_t MVAValue) {

  Int_t subdet = 0;
  if (fabs(fEta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (fPt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;

  Double_t MVACut = -999;
  if (MVABin == 0) MVACut = -0.277;
  if (MVABin == 1) MVACut = -0.566;
  if (MVABin == 2) MVACut = 0.525;
  if (MVABin == 3) MVACut = 0.600;

  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}


Bool_t passEleIsoMVA_LooseWP( Double_t fElePt, Double_t fEleEta, Double_t IsoMVAValue) {

  Int_t subdet = 0;
  if (fabs(fEleEta) < 0.8) subdet = 0;
  else if (fabs(fEleEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (fElePt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t IsoMVACut = -999;
  //Cut Values for MVA-V10 (Detector Based Iso)
  if (MVABin == 0) {
    IsoMVACut = 0.385;
  }
  if (MVABin == 1) {
    IsoMVACut = -0.083;
  }
  if (MVABin == 2) {
    IsoMVACut = -0.573;
  }
  if (MVABin == 3) {
    IsoMVACut = 0.413;
  }
  if (MVABin == 4) {
    IsoMVACut = 0.271;
  }
  if (MVABin == 5) {
    IsoMVACut = 0.135;
  }

  if (IsoMVAValue > IsoMVACut) return kTRUE;
  return kFALSE;
}



Bool_t passEleIDMVA_LooseWP( Double_t fElePt, Double_t fEleEta, Double_t IDMVAValue) {

  Int_t subdet = 0;
  if (fabs(fEleEta) < 0.8) subdet = 0;
  else if (fabs(fEleEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (fElePt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t IDMVACut = -999;
  if (MVABin == 0) {
    IDMVACut =  0.369;
  }
  if (MVABin == 1) {
    IDMVACut = -0.025;
  }
  if (MVABin == 2) {
    IDMVACut = 0.531;
  }
  if (MVABin == 3) {
    IDMVACut = 0.735;
  }
  if (MVABin == 4) {
    IDMVACut = 0.467;
  }
  if (MVABin == 5) {
    IDMVACut = 0.795;
  }

  if (IDMVAValue > IDMVACut) return kTRUE;
  return kFALSE;
}


Bool_t passEleIsoMVA_TightWP( Double_t fElePt, Double_t fEleEta, Double_t IsoMVAValue) {

  Int_t subdet = 0;
  if (fabs(fEleEta) < 0.8) subdet = 0;
  else if (fabs(fEleEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (fElePt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t IsoMVACut = -999;
  //Cut Values for MVA-V10 (Detector Based Iso)
  if (MVABin == 0) {
    IsoMVACut = 0.553;
  }
  if (MVABin == 1) {
    IsoMVACut = -0.237;
  }
  if (MVABin == 2) {
    IsoMVACut = -0.573;
  }
  if (MVABin == 3) {
    IsoMVACut = 0.521;
  }
  if (MVABin == 4) {
    IsoMVACut = 0.531;
  }
  if (MVABin == 5) {
    IsoMVACut = 0.493;
  }

  if (IsoMVAValue > IsoMVACut) return kTRUE;
  return kFALSE;
}



Bool_t passEleIDMVA_TightWP( Double_t fElePt, Double_t fEleEta, Double_t IDMVAValue) {

  Int_t subdet = 0;
  if (fabs(fEleEta) < 0.8) subdet = 0;
  else if (fabs(fEleEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (fElePt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t IDMVACut = -999;
  if (MVABin == 0) {
    IDMVACut =  0.093;
  }
  if (MVABin == 1) {
    IDMVACut = 0.451;
  }
  if (MVABin == 2) {
    IDMVACut = 0.595;
  }
  if (MVABin == 3) {
    IDMVACut = 0.881;
  }
  if (MVABin == 4) {
    IDMVACut = 0.731;
  }
  if (MVABin == 5) {
    IDMVACut = 0.819;
  }

  if (IDMVAValue > IDMVACut) return kTRUE;
  return kFALSE;
}




Bool_t passMuIsoMVA_LooseWP( Double_t fMuPt, Double_t fMuEta, 
                     Bool_t varIsGlobal, Bool_t varIsTracker,
                     Double_t IsoMVAValue
  ) {

  Int_t subdet = 0;
  if (fabs(fMuEta) < 1.5) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (fMuPt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (varIsGlobal && varIsTracker && subdet == 0 && ptBin == 0) MVABin = 0;
  else if (varIsGlobal && varIsTracker && subdet == 0 && ptBin == 1) MVABin = 1;
  else if (varIsGlobal && varIsTracker && subdet == 1 && ptBin == 0) MVABin = 2;
  else if (varIsGlobal && varIsTracker && subdet == 1 && ptBin == 1) MVABin = 3;
  else if (!varIsGlobal && varIsTracker) MVABin = 4;
  else if (varIsGlobal && !varIsTracker) MVABin = 5;
  assert(MVABin >= 0 );

  Double_t IsoMVACut = -999;
  if (MVABin == 0) {
    IsoMVACut = -0.615;
  }
  if (MVABin == 1) {
    IsoMVACut = 0.235;
  }
  if (MVABin == 2) {
    IsoMVACut = -0.825;
  }
  if (MVABin == 3) {
    IsoMVACut = 0.155;
  }
  if (MVABin == 4) {
    IsoMVACut = -0.965;
  }
  if (MVABin == 5) {
    IsoMVACut = -0.989;
  }

  if (IsoMVAValue > IsoMVACut) return kTRUE;
  return kFALSE;
}


Bool_t passMuIDMVA_LooseWP( Double_t fMuPt, Double_t fMuEta, 
                    Bool_t varIsGlobal, Bool_t varIsTracker,
                    Double_t IDMVAValue) {

  Int_t subdet = 0;
  if (fabs(fMuEta) < 1.5) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (fMuPt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (varIsGlobal && varIsTracker && subdet == 0 && ptBin == 0) MVABin = 0;
  else if (varIsGlobal && varIsTracker && subdet == 0 && ptBin == 1) MVABin = 1;
  else if (varIsGlobal && varIsTracker && subdet == 1 && ptBin == 0) MVABin = 2;
  else if (varIsGlobal && varIsTracker && subdet == 1 && ptBin == 1) MVABin = 3;
  else if (!varIsGlobal && varIsTracker) MVABin = 4;
  else if (varIsGlobal && !varIsTracker) MVABin = 5;
  assert(MVABin >= 0 );


  Double_t IDMVACut = -999;
  if (MVABin == 0) {
    IDMVACut =  -0.825;
  }
  if (MVABin == 1) {
    IDMVACut = -0.745;
  }
  if (MVABin == 2) {
    IDMVACut = -0.895;
  }
  if (MVABin == 3) {
    IDMVACut = -0.595;
  }
  if (MVABin == 4) {
    IDMVACut = -0.865;
  }
  if (MVABin == 5) {
    IDMVACut = -0.979;
  }

  if (IDMVAValue > IDMVACut ) return kTRUE;
  return kFALSE;
}


Bool_t passMuIsoMVA_TightWP( Double_t fMuPt, Double_t fMuEta, 
                     Bool_t varIsGlobal, Bool_t varIsTracker,
                     Double_t IsoMVAValue
  ) {

  Int_t subdet = 0;
  if (fabs(fMuEta) < 1.5) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (fMuPt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (varIsGlobal && varIsTracker && subdet == 0 && ptBin == 0) MVABin = 0;
  else if (varIsGlobal && varIsTracker && subdet == 0 && ptBin == 1) MVABin = 1;
  else if (varIsGlobal && varIsTracker && subdet == 1 && ptBin == 0) MVABin = 2;
  else if (varIsGlobal && varIsTracker && subdet == 1 && ptBin == 1) MVABin = 3;
  else if (!varIsGlobal && varIsTracker) MVABin = 4;
  else if (varIsGlobal && !varIsTracker) MVABin = 5;
  assert(MVABin >= 0 );

  Double_t IsoMVACut = -999;
  if (MVABin == 0) {
    IsoMVACut = -0.455;
  }
  if (MVABin == 1) {
    IsoMVACut = 0.695;
  }
  if (MVABin == 2) {
    IsoMVACut = -0.765;
  }
  if (MVABin == 3) {
    IsoMVACut = 0.565;
  }
  if (MVABin == 4) {
    IsoMVACut = -0.965;
  }
  if (MVABin == 5) {
    IsoMVACut = -0.989;
  }

  if (IsoMVAValue > IsoMVACut) return kTRUE;
  return kFALSE;
}


Bool_t passMuIDMVA_TightWP( Double_t fMuPt, Double_t fMuEta, 
                    Bool_t varIsGlobal, Bool_t varIsTracker,
                    Double_t IDMVAValue) {

  Int_t subdet = 0;
  if (fabs(fMuEta) < 1.5) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (fMuPt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (varIsGlobal && varIsTracker && subdet == 0 && ptBin == 0) MVABin = 0;
  else if (varIsGlobal && varIsTracker && subdet == 0 && ptBin == 1) MVABin = 1;
  else if (varIsGlobal && varIsTracker && subdet == 1 && ptBin == 0) MVABin = 2;
  else if (varIsGlobal && varIsTracker && subdet == 1 && ptBin == 1) MVABin = 3;
  else if (!varIsGlobal && varIsTracker) MVABin = 4;
  else if (varIsGlobal && !varIsTracker) MVABin = 5;
  assert(MVABin >= 0 );


  Double_t IDMVACut = -999;
  if (MVABin == 0) {
    IDMVACut =  -0.685;
  }
  if (MVABin == 1) {
    IDMVACut = -0.485;
  }
  if (MVABin == 2) {
    IDMVACut = -0.765;
  }
  if (MVABin == 3) {
    IDMVACut = -0.895;
  }
  if (MVABin == 4) {
    IDMVACut = -0.865;
  }
  if (MVABin == 5) {
    IDMVACut = -0.979;
  }

  if (IDMVAValue > IDMVACut ) return kTRUE;
  return kFALSE;
}

Bool_t passEleIDMVASameBkgAsCiCTight( Double_t fElePt, Double_t fEleSCEta, Double_t MVAValue) {

  Int_t subdet = 0;
  if (fabs(fEleSCEta) < 0.8) subdet = 0;
  else if (fabs(fEleSCEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (fElePt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t MVACut = -999;
  if (MVABin == 0) MVACut = 0.624;
  if (MVABin == 1) MVACut = 0.205;
  if (MVABin == 2) MVACut = 0.539; 
  if (MVABin == 3) MVACut = 0.789;
  if (MVABin == 4) MVACut = 0.413;
  if (MVABin == 5) MVACut = 0.755;  

  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}


Bool_t passElePFIsoSameBkgAsDetIso025( Double_t fPt, Double_t fEta, Double_t PFIso) {

  Int_t subdet = 0;
  if (fabs(fEta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (fPt > 10.0) ptBin = 1;

  Int_t Bin = -1;
  if (subdet == 0 && ptBin == 0) Bin = 0;
  if (subdet == 1 && ptBin == 0) Bin = 1;
  if (subdet == 0 && ptBin == 1) Bin = 2;
  if (subdet == 1 && ptBin == 1) Bin = 3;

  Double_t CutValue = -999;
  //Cut Values for MVA-V10 (Detector Based Iso)
//   if (Bin == 0) CutValue = 0.41;
//   if (Bin == 1) CutValue = 0.37;
//   if (Bin == 2) CutValue = 0.41;
//   if (Bin == 3) CutValue = 0.29;
  if (Bin == 0) CutValue = 0.3;
  if (Bin == 1) CutValue = 0.3;
  if (Bin == 2) CutValue = 0.3;
  if (Bin == 3) CutValue = 0.3;

  if (PFIso < CutValue) return kTRUE;
  return kFALSE;
}

#endif
