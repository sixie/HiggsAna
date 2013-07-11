#ifndef HIGGSANA_UTILS_COMMONTOOLS_HH
#define HIGGSANA_UTILS_COMMONTOOLS_HH

#include <TMath.h>

namespace higgsana
{

  Double_t deltaR(const Double_t eta1, const Double_t phi1, const Double_t eta2, const Double_t phi2);
  Double_t deltaPhi(const Double_t phi1, const Double_t phi2);
  Double_t Eta2Theta(Double_t eta);
  Double_t Theta2Eta(Double_t theta);

//-------------------------------------------------------------------------------------------------------
  Double_t deltaR(const Double_t eta1, const Double_t phi1, const Double_t eta2, const Double_t phi2) { 
    Double_t dphi = deltaPhi(phi1,phi2);
    Double_t deta = eta1 - eta2;
    return sqrt( dphi*dphi + deta*deta);
  }

//-------------------------------------------------------------------------------------------------------
  Double_t deltaPhi(const Double_t phi1, const Double_t phi2) {
    Double_t dphi = phi1-phi2;
    while (dphi > TMath::Pi())
      dphi -= TMath::TwoPi();
    while (dphi <= -TMath::Pi())
      dphi += TMath::TwoPi();

    return dphi;
  }

//-------------------------------------------------------------------------------------------------------
  Double_t Eta2Theta(Double_t eta) {
    return 2.*TMath::ATan(exp(-eta)); 
  }

//-------------------------------------------------------------------------------------------------------
  Double_t Theta2Eta(Double_t theta) 
  { 
    return -TMath::Log(TMath::Tan(theta/2.)); 
  }

}

#endif
