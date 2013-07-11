#ifndef HZZCUTS_HH
#define HZZCUTS_HH

#include <iostream>
#include <cassert>

Double_t _mH[4]          = { 200, 250, 300, 400 };

Double_t _Met_0jet[4]    = {  50,  60,  70,  70 };
Double_t _MtLow_0jet[4]  = { 180, 220, 260, 320 };
Double_t _MtHigh_0jet[4] = { 220, 260, 320, 450 };

Double_t _Met_1jet[4]    = {  60,  60, 100, 100};
Double_t _MtLow_1jet[4]  = { 180, 180, 260, 300 };
Double_t _MtHigh_1jet[4] = { 200, 260, 320, 450 };

Double_t _Met_2jet[4]    = {  80,  80, 100, 100 };
Double_t _MtLow_2jet[4]  = { 180, 180, 240, 300 };
Double_t _MtHigh_2jet[4] = { 200, 260, 320, 450 };


Bool_t passHZZCuts(const Double_t mH,
                   const Int_t    njet,
		   const Double_t met,
		   const Double_t mt
) {

  Int_t im=-1;
  for(Int_t i=0; i<4; i++) {
    if(mH == _mH[i]) {
      im = i;
      break;
    }
  }
  assert(im>=0);
  
  if(njet==0) {
    if(met < _Met_0jet[im])    return kFALSE;
    if(mt  < _MtLow_0jet[im])  return kFALSE;
    if(mt  > _MtHigh_0jet[im]) return kFALSE;
  
  } else if(njet==1) {
    if(met < _Met_1jet[im])    return kFALSE;
    if(mt  < _MtLow_1jet[im])  return kFALSE;
    if(mt  > _MtHigh_1jet[im]) return kFALSE;
  
  } else if(njet>=2) {
    if(met < _Met_2jet[im])    return kFALSE;
    if(mt  < _MtLow_2jet[im])  return kFALSE;
    if(mt  > _MtHigh_2jet[im]) return kFALSE;
  }
    
  return kTRUE;
}
#endif
