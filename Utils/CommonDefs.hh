#ifndef HIGGSANA_UTILS_COMMONDEFS_HH
#define HIGGSANA_UTILS_COMMONDEFS_HH

//*******************************************
//=== Analysis Eras  ====
//*******************************************
enum { 
  kDataEra_NONE,
  kDataEra_2011_MC,
  kDataEra_2012_MC,
  kDataEra_2011_Data,
  kDataEra_2012_Data,
};

//*******************************************
//=== PFIso Types  ====
//*******************************************
enum { 
  kPFChargedIso,
  kPFGammaIso,
  kPFNeutralHadronIso
};

//*******************************************
// Some Constants
//*******************************************
const Double_t MUONMASS     = 105.658369e-3;
const Double_t ELECTRONMASS = 0.51099892e-3;
const Double_t ZMASS       = 91.1876;

#endif
