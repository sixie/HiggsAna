#ifndef HIGGSANA_HZZ4L_UTILS_LEPTONSELECTION_HH
#define HIGGSANA_HZZ4L_UTILS_LEPTONSELECTION_HH

#include "TLorentzVector.h"
#include "HiggsAna/Utils/CommonDefs.hh"
#include "HiggsAna/Utils/CommonTools.hh"
#include "HiggsAna/Utils/IsolationPileupCorrections.hh"
#include "HiggsAna/DataTree/interface/TElectron.hh"
#include "HiggsAna/DataTree/interface/TMuon.hh"
#include "HiggsAna/DataTree/interface/HiggsAnaDefs.hh"

Bool_t PassMuonHZGammaID( const higgsana::TMuon *muon, Bool_t is2011 = kFALSE);  
Bool_t PassMuonHZGammaIso( const higgsana::TMuon *muon, TClonesArray *pfCandidates, 
                           Double_t rho, UInt_t DataEra, 
                           Bool_t printDebug = kFALSE);
Bool_t PassEleHZGammaID( const higgsana::TElectron *ele, 
                         Bool_t printDebug = kFALSE);  
Bool_t PassEleHZGammaIso( const higgsana::TElectron *ele, 
                          Int_t eleIndex,
                          TClonesArray *pfCandidates, 
                          Double_t rho, UInt_t DataEra, 
                          vector<const higgsana::TPFCandidate*> photonsToVeto, 
                          Bool_t printDebug = kFALSE);  
Bool_t PassPhotonHZGammaID( const higgsana::TPhoton *pho, Bool_t is2011 = kFALSE, Bool_t printDebug = kFALSE);  
Bool_t PassPhotonHZGammaIso( const higgsana::TPhoton *pho,     
                             TClonesArray *pfCandidates, 
                             Double_t rho, UInt_t DataEra, 
                             Bool_t printDebug = kFALSE);

//=== FUNCTION DEFINITIONS ======================================================================================


Bool_t PassMuonHZGammaID( const higgsana::TMuon *muon, Bool_t is2011 ) {

  Bool_t pass = kTRUE;

  if(!((muon->typeBits & kGlobal) == kGlobal)) pass = kFALSE;
  if (!muon->PassPFId)                         pass = kFALSE;
  if(muon->muNchi2 >= 10)                      pass = kFALSE;
  if(muon->nValidHits <= 0)                    pass = kFALSE;
  if(muon->nMatch <= 1)                        pass = kFALSE;
  if(fabs(muon->d0)>0.2)                       pass = kFALSE;
  if(fabs(muon->dz)       > 0.5)               pass = kFALSE;
  if(muon->nPixHits	  < 1)                 pass = kFALSE;

  if (is2011) {
    if(muon->trkLayers <= 8)                   pass = kFALSE; 
  } else {
    if(muon->trkLayers <= 5)                   pass = kFALSE; 
  }

  return pass;
}


Bool_t PassMuonHZGammaIso( const higgsana::TMuon *muon,
                           TClonesArray *pfCandidates, 
                           Double_t rho, UInt_t DataEra, 
                           Bool_t printDebug) {

  Double_t pfIso = ComputeMuonPFIso04(muon, pfCandidates, rho, DataEra, printDebug);
  Bool_t pass = ( (pfIso / muon->pt) < 0.12 );
  return pass;

}





Bool_t PassEleHZGammaID( const higgsana::TElectron *ele, Bool_t printDebug) {

  Bool_t pass = kTRUE;
  if (fabs(ele->eta) < 1.566) {
    if (!(
          fabs(ele->deltaEtaIn) < 0.007
          && fabs(ele->deltaPhiIn) < 0.15
          && ele->sigiEtaiEta < 0.01
          && ele->HoverE < 0.12
          && fabs(ele->d0) < 0.02
          && fabs(ele->dz) < 0.2
          && fabs(1.0/ele->EcalEnergy - 1.0/ele->pIn) < 0.05
          && ele->nExpHitsInner <= 1
          && ((ele->isConv & 1024) == 1024)
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
          && fabs(1.0/ele->EcalEnergy - 1.0/ele->pIn) < 0.05
          && ele->nExpHitsInner <= 1
          && ((ele->isConv & 1024) == 1024)
          )
      ) pass = kFALSE;
  }

  if (printDebug) {
    cout << "\nEle: " << ele->pt << " " << ele->eta << " " << ele->phi << " : " << ele->deltaEtaIn << " " << ele->deltaPhiIn << " " << ele->sigiEtaiEta << " " << ele->HoverE << " " << ele->d0 << " " << ele->dz << " " << fabs(1.0/ele->EcalEnergy - 1.0/ele->pIn) << " " << ele->nExpHitsInner << " " << Bool_t((ele->isConv & 1024) == 1024) << " : " << Bool_t(fabs(ele->eta) < 1.566) << " : " << pass << endl;
//     cout << fabs(1.0/(ele->scEt*TMath::CosH(ele->scEta)) - ele->EOverP/ele->EcalEnergy) << " "
//          << fabs(1.0/ele->EcalEnergy - 1.0/ele->pIn) << " " 
//          << fabs(1.0/ele->EcalEnergy - ele->EOverP/ele->EcalEnergy) << " "
//          << fabs(1.0/ele->EcalEnergy - 1.0/(ele->pt * cosh(ele->eta)) ) << " "
//          << endl;
  }

  return pass;  

}
 

Bool_t PassEleHZGammaIso( const higgsana::TElectron *ele, 
                          Int_t eleIndex,
                          TClonesArray *pfCandidates, 
                          Double_t rho, UInt_t DataEra, 
                          Bool_t printDebug) {
    
  Double_t pfIso = ComputeElePFIso04(ele, eleIndex, pfCandidates, rho, DataEra, printDebug);
  Bool_t pass = ( (pfIso / ele->pt) < 0.40 );
  return pass;

}

Bool_t PassPhotonHZGammaID( const higgsana::TPhoton *pho, Bool_t is2011, Bool_t printDebug) {

  Bool_t pass = kTRUE;

  if (fabs(pho->scEta) < 1.479) {
    if (!(pho->passConversionSafeEleVeto)) pass = kFALSE;
    if (!(pho->HoESingleTower < 0.05)) pass = kFALSE;
    if (!(pho->sigiEtaiEta < 0.012)) pass = kFALSE;
  } else {
    if (!(pho->passConversionSafeEleVeto)) pass = kFALSE;
    if (!(pho->HoESingleTower < 0.05)) pass = kFALSE;
    if (!(pho->sigiEtaiEta < 0.034)) pass = kFALSE;
  }

  //some noise cleaning for 2012 data
  if (!is2011) {
    if (pho->scEta > -1.78 && pho->scEta < -1.75 && pho->scEta > 1.36 && pho->scEta < 1.39) pass = kFALSE;    
  }

  if (printDebug) {
    cout << "\nEle: " << pho->et << " " << pho->eta << " " << pho->phi << " : " << pho->scEta << " " << pho->passConversionSafeEleVeto << " " << pho->HoESingleTower << " " << pho->sigiEtaiEta << " : " << pass << endl;
  }

  return pass;

}

Bool_t PassPhotonHZGammaIso( const higgsana::TPhoton *pho,     
                             TClonesArray *pfCandidates, 
                             Double_t rho, UInt_t DataEra, 
                             Bool_t printDebug) {

  Bool_t pass = kTRUE;


  // isolation cuts
  Double_t pfChargedIso = ComputePhotonPFIso03( pho, pfCandidates,kPFChargedIso, rho, DataEra, printDebug);
  Double_t pfGammaIso = ComputePhotonPFIso03( pho, pfCandidates,kPFGammaIso, rho, DataEra, printDebug);
  Double_t pfNeutralHadronIso = ComputePhotonPFIso03( pho, pfCandidates,kPFNeutralHadronIso, rho, DataEra, printDebug);
  
  if (fabs(pho->scEta) < 1.479) {
    if (!(pfChargedIso       < (2.6                  ))) pass = kFALSE;
    if (!(pfGammaIso         < (1.3 + 0.005 * pho->et))) pass = kFALSE;
    if (!(pfNeutralHadronIso < (3.5 + 0.04  * pho->et))) pass = kFALSE;
  } else {
    if (!(pfChargedIso       < (2.3                  ))) pass = kFALSE;
    if (!(pfNeutralHadronIso < (2.9 + 0.04  * pho->et))) pass = kFALSE;
  }

  if (printDebug) {
    cout << "pfiso: " <<  pfChargedIso << " " << pfGammaIso << " " << pfNeutralHadronIso << endl;
  }


  return pass;

}


#endif
