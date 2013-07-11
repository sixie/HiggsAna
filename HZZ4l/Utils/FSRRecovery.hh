#ifndef HIGGSANA_HZZ4L_UTILS_FSRRECOVERY_HH
#define HIGGSANA_HZZ4L_UTILS_FSRRECOVERY_HH

#include <TLorentzVector.h>        
#include <TClonesArray.h>        
#include "HiggsAna/Utils/CommonDefs.hh"
#include "HiggsAna/Utils/CommonTools.hh"
#include "HiggsAna/Utils/LeptonIDCuts.hh"
#include "HiggsAna/Utils/IsolationPileupCorrections.hh"
#include "HiggsAna/DataTree/interface/TElectron.hh"
#include "HiggsAna/DataTree/interface/TMuon.hh"
#include "HiggsAna/DataTree/interface/TPFCandidate.hh"
#include "HiggsAna/DataTree/interface/HiggsAnaDefs.hh"

Bool_t FSRRecovery_TypeI ( higgsana::TElectron * el, 
			   const uint electronOriginalIndex,
			   const uint electronIndex,
			   vector<Bool_t> &FSRRecoveryAttempted,
			   vector<Int_t> &lepType,
			   vector<TLorentzVector> &lepvec,
			   TClonesArray *pfArr,
			   TLorentzVector * Zvec,
			   vector<const higgsana::TPFCandidate*> &photonsToVeto,
                           Bool_t printDebug = kFALSE ) ;

Bool_t FSRRecovery_TypeI ( higgsana::TMuon * mu, 
			   const uint muonOriginalIndex,
			   const uint muonIndex,
			   vector<Bool_t> &FSRRecoveryAttempted,
			   vector<Int_t> &lepType,
			   vector<TLorentzVector> &lepvec,
			   TClonesArray *pfArr,
			   TLorentzVector * Zvec,
			   vector<const higgsana::TPFCandidate*> &photonsToVeto,
                           Bool_t printDebug = kFALSE );

bool recover_typeII_Photon( higgsana::TMuon * mu, 
			    const uint muonIndex,
                            vector<Bool_t> &FSRRecoveryAttempted,
                            vector<TLorentzVector> &lepvec,
                            Bool_t printDebug = kFALSE );


//=== FUNCTION DEFINITIONS ======================================================================================

Bool_t FSRRecovery_TypeI ( higgsana::TElectron * el, 
			   const uint electronOriginalIndex,
			   const uint electronIndex,
			   vector<Bool_t> &FSRRecoveryAttempted,
			   vector<Int_t> &lepType,
			   vector<TLorentzVector> &lepvec,
			   TClonesArray *pfArr,
			   TLorentzVector * Zvec,
			   vector<const higgsana::TPFCandidate*> &photonsToVeto,
                           Bool_t printDebug ) {

  if( FSRRecoveryAttempted[electronIndex] ) return false;

  vector<uint> photonIndices; 
  for( uint i=0; i< uint(pfArr->GetEntries()); i++ ) { 

    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*pfArr)[i]);
    if (!pf->IsPFNoPU) continue;

    if( pf->pfType == eGamma &&
	pf->pt > 2.0 && fabs(pf->eta) < 2.4 ) {

      if( printDebug ) std::cerr << "FSR :: pass preselection ... pt: "<< pf->pt << std::endl;
      float dR = higgsana::deltaR(pf->eta,pf->phi,lepvec[electronIndex].Eta(), lepvec[electronIndex].Phi());
      if( printDebug ) std::cerr << "FSR :: dR = " 
                                 << dR
                                 << std::endl;
      
      //
      // veto if close to an electron SC 
      //
      bool flagEleSC = false;
      for( uint j=0; j<lepvec.size(); j++ ) { 
	if( !lepType[j] == 11 )      continue;
      
        //not necessary because all of them are ID'd and preselected at this stage
	//if( !(lepvec[j].status.looseIDAndPre()) ) continue;
        
	double eeta=lepvec[j].Eta(); double ephi=lepvec[j].Phi(); 
	float dPhi = higgsana::deltaPhi(pf->phi,ephi);
	float dEta = fabs(pf->eta-eeta);
	float tmpDR = higgsana::deltaR(pf->eta,pf->phi, eeta, ephi);
	if(printDebug) cout << "FSR :: comparing to ele, dPhi: " << dPhi 
			    << "\tdEta: " << dEta 
			    << "\tetaPH: " << pf->eta
			    << "\tetaELH: " << eeta
			    << "\ttmpDR:" << tmpDR << endl;
    	if( (dPhi<2.&& dEta<0.05) || tmpDR<0.15 ) { 
          flagEleSC = true;
          break;
	}
	if( flagEleSC ) break;
      }
      if( flagEleSC ) continue;
      if( printDebug ) std::cerr << "FSR :: not matched to an ele SC ... " << std::endl;

      
      //
      // check that input electron is the closest lepton to this photon
      //
      bool FoundCloserLepton=kFALSE;
      for( uint j=0; j<lepvec.size(); j++ ) { 
	if( j == electronIndex ) continue;

        //unncessary
	//if( !(lepvec[j].status.looseIDAndPre()) ) continue;

	float tmp_dR =  higgsana::deltaR(pf->eta,pf->phi, lepvec[j].Eta(), lepvec[j].Phi());
	if( tmp_dR < dR ) {
	  if(printDebug) cout << "FSR :: found closer lepton (j="<<j<<" : "
			      <<tmp_dR<<" vs "<<dR<<") skipping..." << endl;   
          FoundCloserLepton=kTRUE;
          break;
	}
      }
      if( FoundCloserLepton ) continue;


      //
      // Z mass OK?
      //
      TLorentzVector pvec;
      pvec.SetPtEtaPhiM( pf->pt, pf->eta, pf->phi, 0.0);
      float newMass = (pvec + *Zvec).M(); 
      if( !( newMass > 4.   && 
	     newMass < 100. && 
	     (fabs(newMass-ZMASS) < fabs(Zvec->M()-ZMASS))
	     ) ) continue;
      if( printDebug ) std::cerr << "FSR :: improved Zmass  ... " <<
	Zvec->M() << " -> " << newMass << std::endl;

      
      //
      // "keep all photons close to one of the 4L electrons ..."
      //
      if( dR < 0.07 ) { 
	if( printDebug ) std::cerr << "FSR :: dR < 0.07, pushing  ... " << std::endl;
	photonIndices.push_back(i);
      }

      //      
      // "need tighter cuts for other photons ..."
      //
      //if( dR < 0.5 && pf->pt > 4. && ComputePhotonDBetaCorrectedRelativeIsoDr03(pf, 11, electronOriginalIndex, pfArr, printDebug) < 1.0) {
      if( dR < 0.5 && pf->pt > 4. && ComputePhotonNonCorrectedIsoDr03(pf, 11, electronOriginalIndex, pfArr, printDebug) < 1.0) {
	if( printDebug ) std::cerr << "FSR :: tighter cuts, pushing  ... " << std::endl;
	photonIndices.push_back(i);
      }
    }
  }

  float highest_pt  = -1;   int highest_pt_index=-1;
  float smallest_dR = 999.; int smallest_dR_index=-1;
  for( uint i=0; i<photonIndices.size(); i++ ) {
    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)(pfArr->At(photonIndices[i]));
    float dR = higgsana::deltaR(pf->eta,pf->phi, el->eta, el->phi);
    if( pf->pt > highest_pt ) { 
      highest_pt_index = photonIndices[i];
      highest_pt       = pf->pt;
    }
    if( dR < smallest_dR ) { 
      smallest_dR_index = photonIndices[i];
      smallest_dR       = dR;
    }
  }

  const higgsana::TPFCandidate *thepf;
  if( highest_pt > 4. ) { 
    thepf  = (const higgsana::TPFCandidate*)(pfArr->At(highest_pt_index));
    // "... remove it from lepton isolation ..."
    //    PFnoPUflag[highest_pt_index] = 0;
    // TMP, commented flip above for FSR study 
    // gammaMatches[highest_pt_index].push_back(lepvec[electronIndex].index);
    photonsToVeto.push_back(thepf);
  } else if( smallest_dR != 999. ) {
    thepf  = (const higgsana::TPFCandidate*)(pfArr->At(smallest_dR_index));
    // "... remove it from lepton isolation ..."
    // PFnoPUflag[smallest_dR_index] = 0;
    // TMP, commented flip above for FSR study 
    // gammaMatches[smallest_dR_index].push_back(lepvec[electronIndex].index);
    photonsToVeto.push_back(thepf);
  } else { 
    return false;
  }
  
  if( thepf != NULL ) { 
    // add to the electron
    TLorentzVector elvec,phvec,newelvec;
    elvec.SetPtEtaPhiM( el->pt, el->eta, el->phi, ELECTRONMASS);
    phvec.SetPtEtaPhiM( thepf->pt, thepf->eta, thepf->phi, 0.0);
    newelvec = elvec+phvec;
    // don't update the electron object, just simplelepton
    //     el->SetPtEtaPhi        (newelvec.Pt(),
    // 			    newelvec.Eta(),
    // 			    newelvec.Phi());
    lepvec[electronIndex] += phvec;
    FSRRecoveryAttempted[electronIndex] = kTRUE;
    return true;      
  }
  return false;

}


//--------------------------------------------------------------------------------------------------
// typeI = PF IDed photons.  NB : repurpose PFnoPUflag, flip for recovered photons  
// so that they are skipped in the isolation calculation
//--------------------------------------------------------------------------------------------------
Bool_t FSRRecovery_TypeI ( higgsana::TMuon * mu, 
			   const uint muonOriginalIndex,
			   const uint muonIndex,
			   vector<Bool_t> &FSRRecoveryAttempted,
			   vector<Int_t> &lepType,
			   vector<TLorentzVector> &lepvec,
			   TClonesArray *pfArr,
			   TLorentzVector * Zvec,
			   vector<const higgsana::TPFCandidate*> &photonsToVeto,
                           Bool_t printDebug ) 
//--------------------------------------------------------------------------------------------------
{
  if( FSRRecoveryAttempted[muonIndex] ) return false;

  vector<uint> photonIndices; 
  for( uint i=0; i<uint(pfArr->GetEntries()); i++ ) { 
    
    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)((*pfArr)[i]);

    if (!pf->IsPFNoPU) continue;

    if( pf->pfType == eGamma &&
	pf->pt > 2.0 && fabs(pf->eta) < 2.4 ) {

      if( printDebug ) std::cout << "FSR :: pass preselection ... pt: "<< pf->pt << std::endl;

      float dR = higgsana::deltaR(pf->eta,pf->phi,lepvec[muonIndex].Eta(), lepvec[muonIndex].Phi());
      if( printDebug ) std::cout << "FSR :: dR = " 
                                 << dR
                                 << std::endl;

      //
      // veto if close to an electron SC 
      //
      bool flagEleSC = false;
      for( uint j=0; j<lepvec.size(); j++ ) { 

	if( !(lepType[j] == 11 ))      continue;

        //not necessary because all of them are ID'd and preselected at this stage
	//if( !(lepvec[j].status.looseIDAndPre()) ) continue;

	double eeta=lepvec[j].Eta(); double ephi=lepvec[j].Phi(); 
	float dPhi = higgsana::deltaPhi(pf->phi,ephi);
	float dEta = fabs(pf->eta-eeta);
	float tmpDR = higgsana::deltaR(pf->eta,pf->phi, eeta, ephi);
	if(printDebug) cout << "FSR :: comparing to ele, dPhi: " << dPhi 
			    << "\tdEta: " << dEta 
			    << "\tetaPH: " << pf->eta
			    << "\tetaELH: " << eeta
			    << "\ttmpDR:" << tmpDR << endl;

	if( (dPhi<2.&& dEta<0.05) || tmpDR<0.15 ) { 
	    flagEleSC = true;
	    break;
	}
	if( flagEleSC ) break;
      }
      if( flagEleSC ) {
        if(printDebug) cout << "matches to electron. fails \n";
        continue;
      }
      if( printDebug ) std::cout << "FSR :: not matched to an ele SC ... " << std::endl;


      //
      // check that input muon is the closest lepton to this photon
      //
      bool FoundCloserLepton=kFALSE;
      for( uint j=0; j<lepvec.size(); j++ ) { 
	if( j == muonIndex ) continue;

        //unncessary
	//if( !(lepvec[j].status.looseIDAndPre()) ) continue;

	float tmp_dR =  higgsana::deltaR(pf->eta,pf->phi, lepvec[j].Eta(), lepvec[j].Phi());
	if( tmp_dR < dR ) {
	  if(printDebug) cout << "FSR :: found closer lepton (j="<<j<<" : "
			      <<tmp_dR<<" vs "<<dR<<") skipping..." << endl;   
          FoundCloserLepton=kTRUE;
	  break;
	}
      }
      if( FoundCloserLepton ) continue;

      if(printDebug) cout << "no closer lepton\n";
      
      //
      // Z mass OK?
      //
      TLorentzVector pvec;
      pvec.SetPtEtaPhiM( pf->pt, pf->eta, pf->phi, 0.0);
      float newMass = (pvec + *Zvec).M(); 

      if(printDebug) cout << "Zvec: " << Zvec->Pt() << " " << Zvec->Eta() << " " << Zvec->Phi() << " " << endl;
      if(printDebug) cout <<  pf->pt << " " << pf->eta << " " <<  pf->phi << endl;
      if(printDebug) cout << newMass << " " << fabs(newMass-ZMASS)  << " <? " << fabs(Zvec->M()-ZMASS) << endl;

      if( !( newMass > 4.   && 
	     newMass < 100. && 
	     (fabs(newMass-ZMASS) < fabs(Zvec->M()-ZMASS))
	     ) ) continue;
      if( printDebug ) std::cout << "FSR :: improved Zmass  ... " <<
	Zvec->M() << " -> " << newMass << std::endl;

      //
      // "keep all photons close to one of the 4L muons ..."
      //
      if( dR < 0.07 ) { 
	if( printDebug ) std::cout << "FSR :: dR < 0.07, pushing  ... " << std::endl;
	photonIndices.push_back(i);
      }

      //      
      // "need tighter cuts for other photons ..."
      //
      if( printDebug ) std::cout << "FSR :: pass tighter?, pT: " << pf->pt << std::endl;
      //if( dR < 0.5 && pf->pt > 4. && ComputePhotonDBetaCorrectedIsoDr03(pf, 13, muonOriginalIndex, pfArr, printDebug) < 1.0) {
      if( dR < 0.5 && pf->pt > 4. && ComputePhotonNonCorrectedIsoDr03(pf, 13, muonOriginalIndex, pfArr, printDebug) < 1.0) {
	if( printDebug ) std::cout << "FSR :: tighter cuts, pushing  ... " << std::endl;
	photonIndices.push_back(i);
      }
    }
  }

  float highest_pt  = -1;   int highest_pt_index=-1;
  float smallest_dR = 999.; int smallest_dR_index=-1;
  for( uint i=0; i<photonIndices.size(); i++ ) {
    const higgsana::TPFCandidate *pf = (higgsana::TPFCandidate*)(pfArr->At(photonIndices[i]));
    float dR = higgsana::deltaR(pf->eta,pf->phi, mu->eta, mu->phi);
    if( pf->pt > highest_pt ) { 
      highest_pt_index = photonIndices[i];
      highest_pt       = pf->pt;
    }
    if( dR < smallest_dR ) { 
      smallest_dR_index = photonIndices[i];
      smallest_dR       = dR;
    }
  }

  const higgsana::TPFCandidate *thepf;
  if( highest_pt > 4. ) { 
    if(printDebug) std::cout << "FSR :: taking highest pt gamma, index = " <<  highest_pt_index << endl;
    thepf  = (const higgsana::TPFCandidate*)(pfArr->At(highest_pt_index));
    // "... remove it from lepton isolation ..."
    //    PFnoPUflag[highest_pt_index] = 0;
    // TMP, commented flip above for FSR study 
    // gammaMatches[highest_pt_index].push_back(lepvec[muonIndex].index);
    photonsToVeto.push_back(thepf);
  } else if( smallest_dR != 999. ) {
    if(printDebug) std::cout << "FSR :: taking smallest dR gamma, index = " <<  highest_pt_index << endl;
    thepf  = (const higgsana::TPFCandidate*)(pfArr->At(smallest_dR_index));
    // "... remove it from lepton isolation ..."
    //    PFnoPUflag[smallest_dR_index] = 0;
    // TMP, commented flip above for FSR study 
    //gammaMatches[smallest_dR_index].push_back(lepvec[muonIndex].index);
    photonsToVeto.push_back(thepf);
  } else { 
    return false;
  }

  TLorentzVector pvec;
  if( thepf != NULL ) { 
    // add to the muon
    if( printDebug ) cout << "FSR :: before return, oldpT=" << mu->pt << endl;
    TLorentzVector muvec,phvec,newmuvec;
    muvec.SetPtEtaPhiM( mu->pt, mu->eta, mu->phi, MUONMASS);
    phvec.SetPtEtaPhiM( thepf->pt, thepf->eta, thepf->phi, 0.0);
    pvec = phvec;
    newmuvec = muvec+phvec;
    // don't update the muon object, just simplelepton
    //     mu->SetPtEtaPhi        (newmuvec.Pt(),
    // 			    newmuvec.Eta(),
    // 			    newmuvec.Phi());
    lepvec[muonIndex] += phvec;
    FSRRecoveryAttempted[muonIndex] = kTRUE;
    return true;      
  }
  return false;
}


//--------------------------------------------------------------------------------------------------
// typeII = "PFClusters linked to muons"
//--------------------------------------------------------------------------------------------------
bool recover_typeII_Photon( higgsana::TMuon * mu, 
			    const uint muonIndex,
                            vector<Bool_t> &FSRRecoveryAttempted,
                            vector<TLorentzVector> &lepvec,
                            Bool_t printDebug ) 
//--------------------------------------------------------------------------------------------------
{
  if( FSRRecoveryAttempted[muonIndex] ) return false;

  if ( mu->PFMuonEEcal >= 2.0 && mu->PFMuonEtEcal >= 2.0 ) {
    // don't update the muon object, just simplelepton
    //       mu->SetPtEtaPhi        (mu->pt+phpt,
    // 			      mu->Eta(),
    // 			      mu->Phi());

    TLorentzVector pvec;
    pvec.SetPtEtaPhiM( mu->PFMuonEtEcal , mu->eta, mu->phi, 0.0);
    lepvec[muonIndex] += pvec;
    if(printDebug) cout << "FSR :: t2, new pt " <<  lepvec[muonIndex].Pt() << endl;
    FSRRecoveryAttempted[muonIndex] = true;  
    return true;
  }
  
  return false;
}



#endif
