#ifndef ISOLATIONPILEUPCORRECTIONS_HH
#define ISOLATIONPILEUPCORRECTIONS_HH

#include <cassert>
#include "TMath.h"

//*************************************************
//=== Muon Effective Area Pileup Corrections  ====
//*************************************************
enum {        
        kMuNoCorrection,
        kMuTrkIso03, 
	kMuEcalIso03, 
	kMuHcalIso03, 
	kMuTrkIso05, 
	kMuEcalIso05, 
	kMuHcalIso05, 
	kMuChargedIso03, 
	kMuGammaIso03, 
	kMuNeutralHadronIso03, 
	kMuGammaAndNeutralHadronIso03,
	kMuGammaIso03Tight, 
	kMuNeutralHadronIso03Tight, 
	kMuGammaAndNeutralHadronIso03Tight,
	kMuChargedIso04, 
	kMuGammaIso04, 
	kMuNeutralHadronIso04, 
	kMuGammaAndNeutralHadronIso04,
	kMuGammaIso04Tight, 
	kMuNeutralHadronIso04Tight, 
	kMuGammaAndNeutralHadronIso04Tight,
	kMuGammaIsoDR0p0To0p1,
	kMuGammaIsoDR0p1To0p2,
	kMuGammaIsoDR0p2To0p3,
	kMuGammaIsoDR0p3To0p4,
	kMuGammaIsoDR0p4To0p5,
	kMuNeutralHadronIsoDR0p0To0p1,
	kMuNeutralHadronIsoDR0p1To0p2,
	kMuNeutralHadronIsoDR0p2To0p3,
	kMuNeutralHadronIsoDR0p3To0p4,
	kMuNeutralHadronIsoDR0p4To0p5,
        kMuHadEnergy, 
        kMuHoEnergy, 
        kMuEmEnergy, 
        kMuHadS9Energy, 
        kMuHoS9Energy, 
        kMuEmS9Energy,
      };



Double_t MuonEffectiveArea(UInt_t type, Double_t Eta, 
                           UInt_t TargetDataEra) {

  Double_t EffectiveArea = 0;

  if (type == kMuNoCorrection) {
    return 0.0;
  }
  if (TargetDataEra == kDataEra_NONE) {
    return 0.0;
  }
  
  //2012 Data Effective Areas
  else if (TargetDataEra == kDataEra_2012_Data) {
    if (type == kMuGammaIso04){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.50419;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.30582;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.19765;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.28723;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.52529;
      if (fabs(Eta) >= 2.3 )  		      EffectiveArea = 0.48818;
    }
    if (type == kMuNeutralHadronIso04){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.16580;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.25904;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.24695;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.22021;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.34045;
      if (fabs(Eta) >= 2.3 )  		      EffectiveArea = 0.21592;
    }
    if (type == kMuGammaAndNeutralHadronIso04){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.674;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.565;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.442;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.515;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.821;
      if (fabs(Eta) >= 2.3 )  		      EffectiveArea = 0.660;
    }
    if (type == kMuGammaAndNeutralHadronIso03){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.382;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.317;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.242;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.326;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.462;
      if (fabs(Eta) >= 2.3 )  		      EffectiveArea = 0.372;
    }
    if (type == kMuGammaAndNeutralHadronIso04Tight){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.340;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.310;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.315;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.415;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.658;
      if (fabs(Eta) >= 2.3 )  		      EffectiveArea = 0.405;
    }
    if (type == kMuGammaAndNeutralHadronIso03Tight){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.207;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.183;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.177;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.271;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.348;
      if (fabs(Eta) >= 2.3 )  		      EffectiveArea = 0.246;
    }
  }

  //2011 Data Effective Areas
  else if (TargetDataEra == kDataEra_2011_Data) {
    
    if (type == kMuGammaIsoDR0p0To0p1) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.004;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.002;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.002;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.000;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.000;
      if (fabs(Eta) >= 2.3 ) EffectiveArea = 0.005;
    }
    if (type == kMuGammaIsoDR0p1To0p2) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.011;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.008;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.005;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.008;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.008;
      if (fabs(Eta) >= 2.3 ) EffectiveArea = 0.011;
    }
    if (type == kMuGammaIsoDR0p2To0p3) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.023;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.016;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.010;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.014;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.017;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.021;
    }
    if (type == kMuGammaIsoDR0p3To0p4) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.036;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.026;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.017;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.023;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.028;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.032;
    }
    if (type == kMuGammaIsoDR0p4To0p5) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.051;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.037;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.028;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.033;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.042;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.052;
    }
    if (type == kMuNeutralHadronIsoDR0p0To0p1) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.002;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.001;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.001;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.001;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.005;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.007;
    }
    if (type == kMuNeutralHadronIsoDR0p1To0p2) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.005;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.008;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.009;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.009;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.010;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.014;
    }
    if (type == kMuNeutralHadronIsoDR0p2To0p3) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.010;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.015;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.017;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.017;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.019;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.024;
    }
    if (type == kMuNeutralHadronIsoDR0p3To0p4) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.015;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.021;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.024;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.032;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.038;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.038;
    }
    if (type == kMuNeutralHadronIsoDR0p4To0p5) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.020;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.026;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.033;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.045;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.051;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.114;
    }
    /// BEGIN FROM SLIDE 11 OF  https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=188494
    /// NOTE: to be used with the rho from ALL pf candidates within |eta|<2.5
    if (type == kMuGammaIso03){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.049;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.5 ) EffectiveArea = 0.030;
      if (fabs(Eta) >= 1.5 && fabs(Eta) < 2.0 ) EffectiveArea = 0.022;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.034;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.041;
      if (fabs(Eta) >= 2.3 )  		    EffectiveArea = 0.048;
    }
    if (type == kMuGammaIso04){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.085;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.5 ) EffectiveArea = 0.052;
      if (fabs(Eta) >= 1.5 && fabs(Eta) < 2.0 ) EffectiveArea = 0.038;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.055;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.070;
      if (fabs(Eta) >= 2.3 )  		    EffectiveArea = 0.081;
    }
    if (type == kMuNeutralHadronIso03){
  	if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.027;
  	if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.5 ) EffectiveArea = 0.039;
  	if (fabs(Eta) >= 1.5 && fabs(Eta) < 2.0 ) EffectiveArea = 0.044;
  	if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.047;
  	if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.055;
  	if (fabs(Eta) >= 2.3 )		      EffectiveArea = 0.065;
    }
    if (type == kMuNeutralHadronIso04){
  	if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.046;
  	if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.5 ) EffectiveArea = 0.067;
  	if (fabs(Eta) >= 1.5 && fabs(Eta) < 2.0 ) EffectiveArea = 0.074;
  	if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.083;
  	if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.095;
  	if (fabs(Eta) >= 2.3 )		      EffectiveArea = 0.105;
    }
    if (type == kMuGammaAndNeutralHadronIso03){
  	if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.076;
  	if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.5 ) EffectiveArea = 0.070;
  	if (fabs(Eta) >= 1.5 && fabs(Eta) < 2.0 ) EffectiveArea = 0.067;
  	if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.082;
  	if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.097;
  	if (fabs(Eta) >= 2.3 )		      EffectiveArea = 0.115;
    }
    if (type == kMuGammaAndNeutralHadronIso04){
  	if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.132;
  	if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.5 ) EffectiveArea = 0.120;
  	if (fabs(Eta) >= 1.5 && fabs(Eta) < 2.0 ) EffectiveArea = 0.114;
  	if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.139;
  	if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.168;
  	if (fabs(Eta) >= 2.3 )		      EffectiveArea = 0.189;
    }
    /// END FROM SLIDE 11 OF  https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=188494

  } 
  
  //Fall11 MC Effective Areas
  else if (TargetDataEra == kDataEra_2011_MC) {
    if (type == kMuGammaIsoDR0p0To0p1) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.004;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.002;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.003;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.009;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.003;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.011;
    }
    if (type == kMuGammaIsoDR0p1To0p2) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.012;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.008;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.006;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.012;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.019;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.024;
    }
    if (type == kMuGammaIsoDR0p2To0p3) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.026;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.020;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.012;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.022;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.027;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.034;
    }
    if (type == kMuGammaIsoDR0p3To0p4) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.042;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.033;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.022;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.036;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.059;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.068;
    }
    if (type == kMuGammaIsoDR0p4To0p5) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.060;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.043;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.036;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.055;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.092;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.115;
    }
    if (type == kMuNeutralHadronIsoDR0p0To0p1) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.002;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.004;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.004;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.004;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.010;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.014;
    }
    if (type == kMuNeutralHadronIsoDR0p1To0p2) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.005;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.007;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.009;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.009;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.015;
      if (fabs(Eta) >= 2.3  ) 		      EffectiveArea = 0.017;
    }
    if (type == kMuNeutralHadronIsoDR0p2To0p3) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.009;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.015;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.016;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.018;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.022;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.026;
    }
    if (type == kMuNeutralHadronIsoDR0p3To0p4) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.013;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.021;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.026;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.032;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.037;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.042;
    }
    if (type == kMuNeutralHadronIsoDR0p4To0p5) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.017;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.026;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.035;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.046;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.063;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.135;
    }
  }

  else if (TargetDataEra == kDataEra_2012_MC) {
    
    //don't have these yet
    EffectiveArea = 0;
  }


  return EffectiveArea;  

}



enum {
  kEleNoCorrection,
  kEleChargedIso03, 
  kEleNeutralHadronIso03, 
  kEleGammaIso03, 
  kEleGammaIsoVetoEtaStrip03, 
  kEleChargedIso04, 
  kEleNeutralHadronIso04, 
  kEleGammaAndNeutralHadronIso04,
  kEleGammaIso04, 
  kEleGammaIsoVetoEtaStrip04, 
  kEleNeutralHadronIso007, 
  kEleNeutralIso04, 
  kEleHoverE, 
  kEleHcalDepth1OverEcal, 
  kEleHcalDepth2OverEcal,
  kEleGammaIsoDR0p0To0p1,
  kEleGammaIsoDR0p1To0p2,
  kEleGammaIsoDR0p2To0p3,
  kEleGammaIsoDR0p3To0p4,
  kEleGammaIsoDR0p4To0p5,
  kEleNeutralHadronIsoDR0p0To0p1,
  kEleNeutralHadronIsoDR0p1To0p2,
  kEleNeutralHadronIsoDR0p2To0p3,
  kEleNeutralHadronIsoDR0p3To0p4,
  kEleNeutralHadronIsoDR0p4To0p5
};





Double_t ElectronEffectiveArea(UInt_t type, Double_t SCEta, UInt_t TargetDataEra) {

  Double_t EffectiveArea = 0;

  //NoCorrections
  if (type == kEleNoCorrection) {
    return 0.0;
  }    
  if (TargetDataEra == kDataEra_NONE) {
    return 0.0;
  }

  //2012 Data Effective Areas
  else if (TargetDataEra == kDataEra_2012_Data) {
    if (type == kEleGammaAndNeutralHadronIso04) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.19;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.25;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.12;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.21;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.27;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.44;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.52;
    }				
    if (type == kEleGammaIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.051;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.032;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.006;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.007;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.024;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.013;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.013;
    }
    if (type == kEleGammaIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.013;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.013;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.021;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.052;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.066;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.043;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.102;
    }
    if (type == kEleGammaIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.026;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.028;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.041;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.034;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.042;
    }
    if (type == kEleGammaIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.039;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.032;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.024;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.053;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.059;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.069;
    }
    if (type == kEleGammaIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.059;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.045;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.033;
     if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.043;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.056;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.065;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.074;
    }
    if (type == kEleNeutralHadronIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.001;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.006;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.001;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.001;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.003;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.008;
    }
    if (type == kEleNeutralHadronIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.002;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.001;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.004;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.003;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.005;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.006;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.010;
    }
    if (type == kEleNeutralHadronIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.007;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.009;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.009;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.007;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.018;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.028;
    }
    if (type == kEleNeutralHadronIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.011;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.008;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.018;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.026;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.063;
    }
    if (type == kEleNeutralHadronIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.011;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.023;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.038;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.051;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.143;
    }
  } 

  //2011 Data Effective Areas
  else if (TargetDataEra == kDataEra_2011_Data) {
    if (type == kEleGammaAndNeutralHadronIso04) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.180;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.200;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.150;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.190;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.210;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.220;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.290;
    }
    if (type == kEleGammaIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.033;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.005;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.007;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.004;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.000;
    }
    if (type == kEleGammaIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.019;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.042;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.041;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.035;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.041;
    }
    if (type == kEleGammaIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.014;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.029;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.039;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.042;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.048;
    }
    if (type == kEleGammaIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.036;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.029;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.029;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.042;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.047;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.054;
    }
    if (type == kEleGammaIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.051;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.038;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.028;
     if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.036;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.047;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.057;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.059;
    }
    if (type == kEleNeutralHadronIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.001;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.002;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.002;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.000;
    }
    if (type == kEleNeutralHadronIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.005;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.008;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.008;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.006;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.003;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.001;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.003;
    }
    if (type == kEleNeutralHadronIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.014;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.019;
    }
    if (type == kEleNeutralHadronIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.015;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.021;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.025;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.030;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.036;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.038;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.084;
    }
    if (type == kEleNeutralHadronIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.027;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.035;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.045;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.051;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.107;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.228;
    }
  } 
  
  //Fall11 MC Effective Areas
  else if (TargetDataEra == kDataEra_2011_MC) {
    if (type == kEleGammaIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.014;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.004;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.021;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.012;
    }
    if (type == kEleGammaIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.011;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.015;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.042;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.055;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.068;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.067;
    }
    if (type == kEleGammaIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.024;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.038;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.051;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.066;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.080;
    }
    if (type == kEleGammaIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.040;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.032;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.021;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.047;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.066;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.083;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.123;
    }
    if (type == kEleGammaIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.059;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.041;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.037;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.057;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.095;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.123;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.133;
    }
    if (type == kEleNeutralHadronIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.002;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.003;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.000;
    }
    if (type == kEleNeutralHadronIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.006;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.008;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.006;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.005;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.002;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.007;
    }
    if (type == kEleNeutralHadronIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.009;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.014;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.018;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.021;
    }
    if (type == kEleNeutralHadronIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.013;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.019;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.027;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.035;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.037;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.043;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.110;
    }
    if (type == kEleNeutralHadronIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.027;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.036;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.045;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.057;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.123;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.220;
    }
  }

  //Summer12 MC Effective Areas
  else if (TargetDataEra == kDataEra_2012_MC) {
    
    //have not computed these yet.
    EffectiveArea = 0.0;

  }

  return EffectiveArea;  
}




//*************************************************
//=== Photon Effective Area Pileup Corrections  ====
//*************************************************
enum {        
        kPhotonNoCorrection,
        kPhotonChargedIso03, 
	kPhotonGammaIso03, 
	kPhotonNeutralHadronIso03
      };



Double_t PhotonEffectiveArea(UInt_t type, Double_t Eta, 
                           UInt_t TargetDataEra) {

  Double_t EffectiveArea = 0;

  if (type == kMuNoCorrection) {
    return 0.0;
  }
  if (TargetDataEra == kDataEra_NONE) {
    return 0.0;
  }
  
  //Data Effective Areas
  else {
    if (type == kPhotonChargedIso03){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.012;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.010;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.014;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.012;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.016;
      if (fabs(Eta) >= 2.3 && fabs(Eta) < 2.4 )   EffectiveArea = 0.020;
      if (fabs(Eta) >= 2.4 )  		          EffectiveArea = 0.012;
    }
    if (type == kPhotonNeutralHadronIso03){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.030;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.057;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.039;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.015;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.024;
      if (fabs(Eta) >= 2.3 && fabs(Eta) < 2.4 )   EffectiveArea = 0.039;
      if (fabs(Eta) >= 2.4 )  		          EffectiveArea = 0.072;
    }
    if (type == kPhotonGammaIso03){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.148;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.130;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.112;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.216;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.262;
      if (fabs(Eta) >= 2.3 && fabs(Eta) < 2.4 )   EffectiveArea = 0.260;
      if (fabs(Eta) >= 2.4 )  		          EffectiveArea = 0.266;
    }
  }

  return EffectiveArea;  

}



#endif
