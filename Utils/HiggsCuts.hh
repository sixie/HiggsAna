#ifndef HIGGSCUTS_HH
#define HIGGSCUTS_HH

#include <iostream>
#include <cassert>

Double_t _mH[18]          = { 115, 120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500, 550, 600 };

Double_t _Pt1_0jet[18]    = {  20,  20,  25,  25,  27,  30,  34,  36,  38,  40,  55,  70,  80,  90, 110, 120, 130, 140 };
Double_t _Pt2_0jet[18]    = {  10,  10,  10,  15,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25 };
Double_t _Mll_0jet[18]    = {  40,  40,  45,  45,  50,  50,  50,  60,  80,  90, 150, 200, 250, 300, 350, 400, 450, 500 };
Double_t _DPhi_0jet[18]   = { 115, 115,  90,  90,  90,  60,  60,  70,  90, 100, 140, 175, 175, 175, 175, 175, 175, 175 };
Double_t _MtLow_0jet[18]  = {  70,  70,  75,  80,  80,  90, 110, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120 };
Double_t _MtHigh_0jet[18] = { 110, 120, 125, 130, 150, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500, 550, 600 };

Double_t _Pt1_1jet[18]    = {  20,  20,  25,  25,  27,  30,  34,  36,  38,  40,  55,  70,  80,  90, 110, 120, 130, 140 };
Double_t _Pt2_1jet[18]    = {  10,  10,  10,  15,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25 };
Double_t _Mll_1jet[18]    = {  40,  40,  45,  45,  50,  50,  50,  60,  80,  90, 150, 200, 250, 300, 350, 400, 450, 500 };
Double_t _DPhi_1jet[18]   = { 115, 115,  90,  90,  90,  60,  60,  70,  90, 100, 140, 175, 175, 175, 175, 175, 175, 175 };
Double_t _MtLow_1jet[18]  = {  70,  70,  75,  80,  80,  90, 110, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120 };
Double_t _MtHigh_1jet[18] = { 110, 120, 125, 130, 150, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500, 550, 600 };

Bool_t passHCuts(const Double_t mH,
                 const Int_t    njet,
		 const Double_t pt1,
		 const Double_t pt2,
		 const Double_t mll,
		 const Double_t dphi,
		 const Double_t mt
) {

  Int_t im=-1;
  for(Int_t i=0; i<18; i++) {
    if(mH == _mH[i]) {
      im = i;
      break;
    }
  }
  assert(im>=0);
  assert(njet<3);
  
  if(njet==0) {
    if(pt1  < _Pt1_0jet[im])  return kFALSE;
    if(pt2  < _Pt2_0jet[im])  return kFALSE;
    if(mll  > _Mll_0jet[im])  return kFALSE;
    if(dphi > _DPhi_0jet[im]) return kFALSE;
    
    if(mt < _MtLow_0jet[im] || mt > _MtHigh_0jet[im]) return kFALSE; 
    
  } else if(njet==1) {
    if(pt1  < _Pt1_1jet[im])  return kFALSE;
    if(pt2  < _Pt2_1jet[im])  return kFALSE;
    if(mll  > _Mll_1jet[im])  return kFALSE;
    if(dphi > _DPhi_1jet[im]) return kFALSE;
    
    if(mt < _MtLow_1jet[im] || mt > _MtHigh_1jet[im]) return kFALSE;   
  
  } else if(njet==2) {
    if(mH<=200 && mll>100) return kFALSE;
  }
  
  return kTRUE;
}

void printHCuts(const Double_t mH, const Int_t njet) {
  Int_t im=-1;
  for(Int_t i=0; i<18; i++) {
    if(mH == _mH[i]) {
      im = i;
      break;
    }
  }
  assert(im>=0);
  assert(njet<3);
  
  cout << "   mH = " << mH << ", njet = " << njet << endl;
  if(njet==0) {
    cout << " *  leading lepton pT > " << _Pt1_0jet[im] << endl;
    cout << " * trailing lepton pT > " << _Pt2_0jet[im] << endl;
    cout << " *      dilepton mass < " << _Mll_0jet[im] << endl;
    cout << " *      dilepton dphi < " << _DPhi_0jet[im] << endl;
    cout << " *  transverse mass in [" << _MtLow_0jet[im] << ", " << _MtHigh_0jet[im] << "]" << endl; 
    
  } else if(njet==1) {
    cout << " *  leading lepton pT > " << _Pt1_1jet[im] << endl;
    cout << " * trailing lepton pT > " << _Pt2_1jet[im] << endl;
    cout << " *      dilepton mass < " << _Mll_1jet[im] << endl;
    cout << " *      dilepton dphi < " << _DPhi_1jet[im] << endl;
    cout << " *  transverse mass in [" << _MtLow_1jet[im] << ", " << _MtHigh_1jet[im] << "]" << endl;   
  
  } else if(njet==2) {
    if(mH<=200) cout << " * dilepton mass < 100" << endl;
  }
}
#endif
