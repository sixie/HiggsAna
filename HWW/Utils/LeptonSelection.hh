#ifndef HIGGSANA_HZZ4L_UTILS_LEPTONSELECTION_HH
#define HIGGSANA_HZZ4L_UTILS_LEPTONSELECTION_HH

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

Double_t IDIsoMVAV4( const higgsana::TElectron *ele,
                     TClonesArray *pfCandidates, 
                     Double_t rho, UInt_t DataEra,  
                     EGammaMvaEleEstimator* eleMVAEstimator, 
                     Bool_t printDebug = kFALSE);

Bool_t PassEleHWWICHEP2012IDMVA( const higgsana::TElectron *ele,                                  
                                 EGammaMvaEleEstimator* eleMVAEstimator,
                                 Bool_t printDebug = kFALSE);  

Bool_t PassEleHWWIDIsoMVAV3( const higgsana::TElectron *ele,                                  
                             TClonesArray *pfCandidates, 
                             Double_t rho, UInt_t DataEra, 
                             EGammaMvaEleEstimator* eleMVAEstimator,
                             Bool_t printDebug = kFALSE);
Bool_t PassEleHWWIDIsoMVAV4( const higgsana::TElectron *ele,                                  
                             TClonesArray *pfCandidates, 
                             Double_t rho, UInt_t DataEra, 
                             EGammaMvaEleEstimator* eleMVAEstimator,
                             Bool_t printDebug = kFALSE );

Bool_t PassEleHWWIDIsoMVAV5( const higgsana::TElectron *ele,                                  
                             TClonesArray *pfCandidates, 
                             Double_t rho, UInt_t DataEra, 
                             EGammaMvaEleEstimator* eleMVAEstimator,
                             Bool_t printDebug = kFALSE );


//=== FUNCTION DEFINITIONS ======================================================================================
Bool_t PassEleHWWICHEP2012IDMVA( const higgsana::TElectron *ele,                                      
                                 EGammaMvaEleEstimator* eleMVAEstimator, 
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

  double _OneMinusE1x5E5x5 = 1 - ele->SeedE1x5OverE5x5;
  if( _OneMinusE1x5E5x5 < -1. ) _OneMinusE1x5E5x5 = -1.;
  if( _OneMinusE1x5E5x5 >  2. ) _OneMinusE1x5E5x5 =  2.;

  double _r9 = (ele->R9 > 5 ) ? 5 : ele->R9;

  double _h_o_e = ele->HoverE;

  double _e_o_p = (ele->EOverP > 20. ) ? 20 : ele->EOverP;

  double _IoEmIoP =  1.0/(ele->EcalEnergy) - 1/(ele->pt*TMath::CosH(ele->eta)); 

  double _eeleclu_o_pout = (ele->EEleClusterOverPout > 20. ) ? 20. : ele->EEleClusterOverPout;

  double _epreoraw = ele->PreShowerOverRaw;
  
  double _d0 = ele->d0;
  double _ip3d = ele->ip3d;
  
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
                                             _OneMinusE1x5E5x5,
                                             _r9,
                                             _h_o_e,
                                             _e_o_p,
                                             _IoEmIoP,
                                             _eeleclu_o_pout,
                                             _epreoraw,
                                             _d0,
                                             _ip3d,
                                             ele->scEta,	
                                             ele->pt,
                                             printDebug );
  
  bool pass = false;
  
  Int_t subdet = 0;
  if (fabs(ele->scEta) < 0.8) subdet = 0;
  else if (fabs(ele->scEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ele->pt > 20) ptBin = 1;
	
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0; // eta<0.8, pt<20
  if (subdet == 1 && ptBin == 0) MVABin = 1; // 0.8<eta<1.479, pt<20
  if (subdet == 2 && ptBin == 0) MVABin = 2; // eta>1.478, pt<20
  if (subdet == 0 && ptBin == 1) MVABin = 3; // eta<0.8, pt>20
  if (subdet == 1 && ptBin == 1) MVABin = 4; // 0.8<eta<1.479, pt>20
  if (subdet == 2 && ptBin == 1) MVABin = 5; // eta>1.478, pt>20    

  if( MVABin == 0 && mvaval > 0.000 ) pass = true;
  if( MVABin == 1 && mvaval > 0.100 ) pass = true;
  if( MVABin == 2 && mvaval > 0.620 ) pass = true;
  if( MVABin == 3 && mvaval > 0.940 ) pass = true;
  if( MVABin == 4 && mvaval > 0.850 ) pass = true;
  if( MVABin == 5 && mvaval > 0.920 ) pass = true;

  return pass;

}






Bool_t PassEleHWWIDIsoMVAV3( const higgsana::TElectron *ele,
                                        TClonesArray *pfCandidates, 
                                        Double_t rho, UInt_t DataEra,  
                                        EGammaMvaEleEstimator* eleMVAEstimator, 
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

  double _OneMinusE1x5E5x5 = 1 - ele->SeedE1x5OverE5x5;
  if( _OneMinusE1x5E5x5 < -1. ) _OneMinusE1x5E5x5 = -1.;
  if( _OneMinusE1x5E5x5 >  2. ) _OneMinusE1x5E5x5 =  2.;

  double _r9 = (ele->R9 > 5 ) ? 5 : ele->R9;

  double _h_o_e = ele->HoverE;

  double _e_o_p = (ele->EOverP > 20. ) ? 20 : ele->EOverP;

  double _IoEmIoP =  (1.0/ele->EcalEnergy - 1.0 / ele->pIn);

  double _eeleclu_o_pout = (ele->EEleClusterOverPout > 20. ) ? 20. : ele->EEleClusterOverPout;

  double _epreoraw = ele->PreShowerOverRaw;
  
  double _d0 = ele->d0;
  double _ip3d = ele->ip3d;

  double _ChargedIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.0, 0.1);
  double _ChargedIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.1, 0.2);
  double _ChargedIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.2, 0.3);
  double _ChargedIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.3, 0.4);
  double _ChargedIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.4, 0.5);
  double _GammaIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.0, 0.1);
  double _GammaIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.1, 0.2);
  double _GammaIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.2, 0.3);
  double _GammaIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.3, 0.4);
  double _GammaIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.4, 0.5);
  double _NeutralHadronIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.0, 0.1);
  double _NeutralHadronIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.1, 0.2);
  double _NeutralHadronIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.2, 0.3);
  double _NeutralHadronIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.3, 0.4);
  double _NeutralHadronIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.4, 0.5);

  double _rho = rho;

  double mvaval = eleMVAEstimator->IDIsoCombinedMvaValue( _fbrem,
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
                                             _OneMinusE1x5E5x5,
                                             _r9,
                                             _h_o_e,
                                             _e_o_p,
                                             _IoEmIoP,
                                             _eeleclu_o_pout,
                                             _epreoraw,
                                             _d0,
                                             _ip3d,
                                             _ChargedIso_DR0p0To0p1,
                                             _ChargedIso_DR0p1To0p2,
                                             _ChargedIso_DR0p2To0p3,
                                             _ChargedIso_DR0p3To0p4,
                                             _ChargedIso_DR0p4To0p5,
                                             _GammaIso_DR0p0To0p1,
                                             _GammaIso_DR0p1To0p2,
                                             _GammaIso_DR0p2To0p3,
                                             _GammaIso_DR0p3To0p4,
                                             _GammaIso_DR0p4To0p5,
                                             _NeutralHadronIso_DR0p0To0p1,
                                             _NeutralHadronIso_DR0p1To0p2,
                                             _NeutralHadronIso_DR0p2To0p3,
                                             _NeutralHadronIso_DR0p3To0p4,
                                             _NeutralHadronIso_DR0p4To0p5,
                                             _rho,
                                             ele->scEta,	
                                             ele->pt,
                                             printDebug );
  
  bool pass = false;
  
  Int_t subdet = 0;
  if (fabs(ele->scEta) < 0.8) subdet = 0;
  else if (fabs(ele->scEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ele->pt > 20) ptBin = 1;
	
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0; // eta<0.8, pt<20
  if (subdet == 1 && ptBin == 0) MVABin = 1; // 0.8<eta<1.479, pt<20
  if (subdet == 2 && ptBin == 0) MVABin = 2; // eta>1.478, pt<20
  if (subdet == 0 && ptBin == 1) MVABin = 3; // eta<0.8, pt>20
  if (subdet == 1 && ptBin == 1) MVABin = 4; // 0.8<eta<1.479, pt>20
  if (subdet == 2 && ptBin == 1) MVABin = 5; // eta>1.478, pt>20    

  if( MVABin == 0 && mvaval > 0.118 ) pass = true;
  if( MVABin == 1 && mvaval > 0.214 ) pass = true;
  if( MVABin == 2 && mvaval > 0.516 ) pass = true;
  if( MVABin == 3 && mvaval > 0.935 ) pass = true;
  if( MVABin == 4 && mvaval > 0.889 ) pass = true;
  if( MVABin == 5 && mvaval > 0.872 ) pass = true;

  return pass;

}



Bool_t PassEleHWWIDIsoMVAV4( const higgsana::TElectron *ele,
                                        TClonesArray *pfCandidates, 
                                        Double_t rho, UInt_t DataEra,  
                                        EGammaMvaEleEstimator* eleMVAEstimator, 
                                        Bool_t printDebug) {
  
  double mvaval = IDIsoMVAV4( ele, pfCandidates, 
                              rho, DataEra,  
                              eleMVAEstimator, 
                              printDebug);

  bool pass = false;
  
  Int_t subdet = 0;
  if (fabs(ele->scEta) < 0.8) subdet = 0;
  else if (fabs(ele->scEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ele->pt > 20) ptBin = 1;
	
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0; // eta<0.8, pt<20
  if (subdet == 1 && ptBin == 0) MVABin = 1; // 0.8<eta<1.479, pt<20
  if (subdet == 2 && ptBin == 0) MVABin = 2; // eta>1.478, pt<20
  if (subdet == 0 && ptBin == 1) MVABin = 3; // eta<0.8, pt>20
  if (subdet == 1 && ptBin == 1) MVABin = 4; // 0.8<eta<1.479, pt>20
  if (subdet == 2 && ptBin == 1) MVABin = 5; // eta>1.478, pt>20    

  if( MVABin == 0 && mvaval > 0.106 ) pass = true;
  if( MVABin == 1 && mvaval > 0.143 ) pass = true;
  if( MVABin == 2 && mvaval > 0.414 ) pass = true;
  if( MVABin == 3 && mvaval > 0.935 ) pass = true;
  if( MVABin == 4 && mvaval > 0.889 ) pass = true;
  if( MVABin == 5 && mvaval > 0.855 ) pass = true;

  return pass;

}


Double_t IDIsoMVAV4( const higgsana::TElectron *ele,
                     TClonesArray *pfCandidates, 
                     Double_t rho, UInt_t DataEra,  
                     EGammaMvaEleEstimator* eleMVAEstimator, 
                     Bool_t printDebug) {
  

  double _fbrem = (ele->fBrem < -1. ) ? -1. : ele->fBrem;

  double _kftrk_chisq = ele->KFTrackChi2OverNdof;
  if( _kftrk_chisq > 10. ) _kftrk_chisq = 10.;

  double _kftrk_nhits = ele->KFTrackNLayersWithMeasurement; 

  double _gsftrk_chisq = ele->GsfTrackChi2OverNdof;
  if( _gsftrk_chisq > 200. ) _gsftrk_chisq = 200.;
 
  double _deta = fabs(ele->deltaEtaIn);
  if( _deta > 0.06 ) _deta = 0.06;

  double _dphi = ele->deltaPhiIn;

  double _detacalo = ele->dEtaCalo;  

  double _sigieie = ele->sigiEtaiEta;

  double _sigiphiiphi = ele->sigiPhiiPhi;
  if( isnan(_sigiphiiphi ) ) _sigiphiiphi = 0.;

  double _etawidth = ele->SCEtaWidth; 

  double _phiwidth = ele->SCPhiWidth;

  double _OneMinusE1x5E5x5 = 1 - ele->SeedE1x5OverE5x5;
  if( _OneMinusE1x5E5x5 < -1. ) _OneMinusE1x5E5x5 = -1.;
  if( _OneMinusE1x5E5x5 >  2. ) _OneMinusE1x5E5x5 =  2.;

  double _r9 = (ele->R9 > 5 ) ? 5 : ele->R9;

  double _h_o_e = ele->HoverE;

  double _e_o_p = (ele->EOverP > 20. ) ? 20 : ele->EOverP;

  double _IoEmIoP =  (1.0/ele->EcalEnergy - 1.0 / ele->pIn);

  double _eeleclu_o_pout = (ele->EEleClusterOverPout > 20. ) ? 20. : ele->EEleClusterOverPout;

  double _epreoraw = ele->PreShowerOverRaw;
  
  double _d0 = ele->d0;
  double _ip3d = ele->ip3d;

  double _ChargedIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.0, 0.1)/ele->pt;
  double _ChargedIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.1, 0.2)/ele->pt;
  double _ChargedIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.2, 0.3)/ele->pt;
  double _ChargedIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.3, 0.4)/ele->pt;
  double _ChargedIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.4, 0.5)/ele->pt;
  double _GammaIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.0, 0.1)/ele->pt;
  double _GammaIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.1, 0.2)/ele->pt;
  double _GammaIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.2, 0.3)/ele->pt;
  double _GammaIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.3, 0.4)/ele->pt;
  double _GammaIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.4, 0.5)/ele->pt;
  double _NeutralHadronIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.0, 0.1)/ele->pt;
  double _NeutralHadronIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.1, 0.2)/ele->pt;
  double _NeutralHadronIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.2, 0.3)/ele->pt;
  double _NeutralHadronIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.3, 0.4)/ele->pt;
  double _NeutralHadronIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.4, 0.5)/ele->pt;

  double _rho = rho;

  double mvaval = eleMVAEstimator->IDIsoCombinedMvaValue( _fbrem,
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
                                             _OneMinusE1x5E5x5,
                                             _r9,
                                             _h_o_e,
                                             _e_o_p,
                                             _IoEmIoP,
                                             _eeleclu_o_pout,
                                             _epreoraw,
                                             _d0,
                                             _ip3d,
                                             _ChargedIso_DR0p0To0p1,
                                             _ChargedIso_DR0p1To0p2,
                                             _ChargedIso_DR0p2To0p3,
                                             _ChargedIso_DR0p3To0p4,
                                             _ChargedIso_DR0p4To0p5,
                                             _GammaIso_DR0p0To0p1,
                                             _GammaIso_DR0p1To0p2,
                                             _GammaIso_DR0p2To0p3,
                                             _GammaIso_DR0p3To0p4,
                                             _GammaIso_DR0p4To0p5,
                                             _NeutralHadronIso_DR0p0To0p1,
                                             _NeutralHadronIso_DR0p1To0p2,
                                             _NeutralHadronIso_DR0p2To0p3,
                                             _NeutralHadronIso_DR0p3To0p4,
                                             _NeutralHadronIso_DR0p4To0p5,
                                             _rho,
                                             ele->scEta,	
                                             ele->pt,
                                             printDebug );
  
  return mvaval;
}



Bool_t PassEleHWWIDIsoMVAV5( const higgsana::TElectron *ele,
                                        TClonesArray *pfCandidates, 
                                        Double_t rho, UInt_t DataEra,  
                                        EGammaMvaEleEstimator* eleMVAEstimator, 
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

  double _OneMinusE1x5E5x5 = 1 - ele->SeedE1x5OverE5x5;
  if( _OneMinusE1x5E5x5 < -1. ) _OneMinusE1x5E5x5 = -1.;
  if( _OneMinusE1x5E5x5 >  2. ) _OneMinusE1x5E5x5 =  2.;

  double _r9 = (ele->R9 > 5 ) ? 5 : ele->R9;

  double _h_o_e = ele->HoverE;

  double _e_o_p = (ele->EOverP > 20. ) ? 20 : ele->EOverP;

  double _IoEmIoP =  (1.0/ele->EcalEnergy - 1.0 / ele->pIn);

  double _eeleclu_o_pout = (ele->EEleClusterOverPout > 20. ) ? 20. : ele->EEleClusterOverPout;

  double _epreoraw = ele->PreShowerOverRaw;
  
  double _d0 = ele->d0;
  double _ip3d = ele->ip3d;

  double _ChargedIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.0, 0.1);
  double _ChargedIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.1, 0.2);
  double _ChargedIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.2, 0.3);
  double _ChargedIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.3, 0.4);
  double _ChargedIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.4, 0.5);
  double _GammaIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.0, 0.1);
  double _GammaIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.1, 0.2);
  double _GammaIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.2, 0.3);
  double _GammaIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.3, 0.4);
  double _GammaIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.4, 0.5);
  double _NeutralHadronIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.0, 0.1);
  double _NeutralHadronIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.1, 0.2);
  double _NeutralHadronIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.2, 0.3);
  double _NeutralHadronIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.3, 0.4);
  double _NeutralHadronIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.4, 0.5);

  double _rho = rho;

  double mvaval = eleMVAEstimator->IDIsoCombinedMvaValue( _fbrem,
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
                                             _OneMinusE1x5E5x5,
                                             _r9,
                                             _h_o_e,
                                             _e_o_p,
                                             _IoEmIoP,
                                             _eeleclu_o_pout,
                                             _epreoraw,
                                             _d0,
                                             _ip3d,
                                             _ChargedIso_DR0p0To0p1,
                                             _ChargedIso_DR0p1To0p2,
                                             _ChargedIso_DR0p2To0p3,
                                             _ChargedIso_DR0p3To0p4,
                                             _ChargedIso_DR0p4To0p5,
                                             _GammaIso_DR0p0To0p1,
                                             _GammaIso_DR0p1To0p2,
                                             _GammaIso_DR0p2To0p3,
                                             _GammaIso_DR0p3To0p4,
                                             _GammaIso_DR0p4To0p5,
                                             _NeutralHadronIso_DR0p0To0p1,
                                             _NeutralHadronIso_DR0p1To0p2,
                                             _NeutralHadronIso_DR0p2To0p3,
                                             _NeutralHadronIso_DR0p3To0p4,
                                             _NeutralHadronIso_DR0p4To0p5,
                                             _rho,
                                             ele->scEta,	
                                             ele->pt,
                                             printDebug );
  
  bool pass = false;
  
  Int_t subdet = 0;
  if (fabs(ele->scEta) < 0.8) subdet = 0;
  else if (fabs(ele->scEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ele->pt > 20) ptBin = 1;
	
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0; // eta<0.8, pt<20
  if (subdet == 1 && ptBin == 0) MVABin = 1; // 0.8<eta<1.479, pt<20
  if (subdet == 2 && ptBin == 0) MVABin = 2; // eta>1.478, pt<20
  if (subdet == 0 && ptBin == 1) MVABin = 3; // eta<0.8, pt>20
  if (subdet == 1 && ptBin == 1) MVABin = 4; // 0.8<eta<1.479, pt>20
  if (subdet == 2 && ptBin == 1) MVABin = 5; // eta>1.478, pt>20    

  if( MVABin == 0 && mvaval > 0.116 ) pass = true;
  if( MVABin == 1 && mvaval > 0.211 ) pass = true;
  if( MVABin == 2 && mvaval > 0.511 ) pass = true;
  if( MVABin == 3 && mvaval > 0.935 ) pass = true;
  if( MVABin == 4 && mvaval > 0.888 ) pass = true;
  if( MVABin == 5 && mvaval > 0.871 ) pass = true;

  return pass;

}

#endif
