#ifndef HIGGSANA_NTUPLER_HIGGSANADEFS_HH 
#define HIGGSANA_NTUPLER_HIGGSANADEFS_HH

enum EMuType 
{ 
  kGlobal     = 1, 
  kTracker    = 2, 
  kStandalone = 4
};

enum EQualityBit
{ 
  // descriptions from DataFormats/MuonReco/interface/MuonSelectors.h
  kAll  			    = 0x000001,  // dummy options - always true
  kAllGlobalMuons		    = 0x000002,  // checks isGlobalMuon flag
  kAllStandAloneMuons		    = 0x000004,  // checks isStandAloneMuon flag
  kAllTrackerMuons		    = 0x000008,  // checks isTrackerMuon flag
  kTrackerMuonArbitrated	    = 0x000010,  // resolve ambiguity of sharing segments
  kAllArbitrated		    = 0x000020,  // all muons with the tracker muon arbitrated
  kGlobalMuonPromptTight	    = 0x000040,  // global muons with tighter fit requirements
  kTMLastStationLoose		    = 0x000080,  // penetration depth loose selector
  kTMLastStationTight		    = 0x000100,  // penetration depth tight selector
  kTM2DCompatibilityLoose	    = 0x000200,  // likelihood based loose selector
  kTM2DCompatibilityTight	    = 0x000400,  // likelihood based tight selector
  kTMOneStationLoose		    = 0x000800,  // require one well matched segment
  kTMOneStationTight		    = 0x001000,  // require one well matched segment
  kTMLastStationOptimizedLowPtLoose = 0x002000,  // combination of TMLastStation and TMOneStation
  kTMLastStationOptimizedLowPtTight = 0x004000,  // combination of TMLastStation and TMOneStation
  kGMTkChiCompatibility 	    = 0x008000,  // require tk stub have good chi2 relative to glb track
  kGMStaChiCompatibility	    = 0x010000,  // require sta stub have good chi2 compatibility relative to glb track
  kGMTkKinkTight		    = 0x020000,  // require a small kink value in the tracker stub
  kTMLastStationAngLoose	    = 0x040000,  // TMLastStationLoose with additional angular cuts
  kTMLastStationAngTight	    = 0x080000,  // TMLastStationTight with additional angular cuts
  kTMOneStationAngLoose 	    = 0x100000,  // TMOneStationLoose with additional angular cuts
  kTMOneStationAngTight 	    = 0x200000,  // TMOneStationTight with additional angular cuts
  //The two algorithms that follow are identical to what were known as
  //TMLastStationOptimizedLowPt* (sans the Barrel) as late as revision
  //1.7 of this file. The names were changed because indeed the low pt
  //optimization applies only to the barrel region, whereas the sel-
  //ectors above are more efficient at low pt in the endcaps, which is
  //what we feel is more suggestive of the algorithm name. This will be
  //less confusing for future generations of CMS members, I hope...
  kTMLastStationOptimizedBarrelLowPtLoose = 0x400000,  // combination of TMLastStation and TMOneStation but with low pT optimization in barrel only
  kTMLastStationOptimizedBarrelLowPtTight = 0x800000   // combination of TMLastStation and TMOneStation but with low pT optimization in barrel only
}; 

enum ETriggerBit
{  
  kHLT_Mu24                                                = 1UL<<0,
  kHLT_Mu30                                                = 1UL<<1,
  kHLT_Mu8                                                 = 1UL<<2,
  kHLT_Mu15                                                = 1UL<<3,

  kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT              = 1UL<<4,
  kHLT_Ele8                                                = 1UL<<5,
  kHLT_Ele8_CaloIdL_CaloIsoVL                              = 1UL<<6,
  kHLT_Ele17_CaloIdL_CaloIsoVL                             = 1UL<<7,

  kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40                        = 1UL<<8,
  kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL       = 1UL<<9,
  kHLT_Mu8_Jet40                                           = 1UL<<10,
  kHLT_Mu8_Photon20_CaloIdVT_IsoT                          = 1UL<<11,
 
  kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL      = 1UL<<12,
  kHLT_Ele32_CaloIdL_CaloIsoVL_SC17                        = 1UL<<13,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 = 1UL<<14,
  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL  = 1UL<<15,

  kHLT_DoubleMu7                                           = 1UL<<16,
  kHLT_DoubleMu6                                           = 1UL<<17,
  kHLT_Mu17_Ele8_CaloIdL                                   = 1UL<<18,
  kHLT_Mu8_Ele17_CaloIdL                                   = 1UL<<19,

  kHLT_Photon26_IsoVL_Photon18                             = 1UL<<20,
  kHLT_Photon26_CaloIdL_IsoVL_Photon18                     = 1UL<<21,	
  kHLT_Photon26_IsoVL_Photon18_IsoVL                       = 1UL<<22,
  kHLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL       = 1UL<<23,

  kHLT_IsoMu12                                             = 1UL<<24,
  kHLT_IsoMu17                                             = 1UL<<25,
  kHLT_IsoMu24                                             = 1UL<<26,  
  kHLT_IsoMu30                                             = 1UL<<27,

  kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL             = 1UL<<28,
  kHLT_Mu17_Mu8                                            = 1UL<<29,
  kHLT_Mu3                                                 = 1UL<<30,
  kHLT_Mu5                                                 = 1UL<<31,

  kHLT_Jet30                                               = 1UL<<32,
  kHLT_Jet60                                               = 1UL<<33,
  kHLT_Jet80                                               = 1UL<<34,
  kHLT_Jet110                                              = 1UL<<35,

  kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50 = 1UL<<36,
  kHLT_Jet370                                              = 1UL<<37,

  kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL             = 1UL<<40,
  kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30       = 1UL<<41,
  kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL            = 1UL<<42,
  kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30      = 1UL<<43,

  kHLT_Mu13_Mu8                                            = 1UL<<44,
  kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL        = 1UL<<45,
  kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL        = 1UL<<46,

  kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL                         = 1UL<<48,
  kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL                         = 1UL<<49,
  kHLT_Mu17_TkMu8                                          = 1UL<<50,
  kHLT_Mu17                                                = 1UL<<51,

  kHLT_DoubleMu3                                           = 1UL<<52,
  kHLT_Ele10_LW                                            = 1UL<<53,
  kHLT_Ele15_SW                                            = 1UL<<54,
  kHLT_Ele17_SW                                            = 1UL<<55


};


enum ETriggerObjectBit
{  
  kHLTObject_Mu24                                                = 1UL<<0, 
  kHLTObject_Mu30                                                = 1UL<<1, 
  kHLTObject_Mu8                                                 = 1UL<<2, 
  kHLTObject_Mu15                                                = 1UL<<3, 

  kHLTObject_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT              = 1UL<<4, 
  kHLTObject_Ele8                                                = 1UL<<5, 
  kHLTObject_Ele8_CaloIdL_CaloIsoVL                              = 1UL<<6, 
  kHLTObject_Ele17_CaloIdL_CaloIsoVL                             = 1UL<<7, 

  kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL            = 1UL<<8, 
  kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL             = 1UL<<9, 
  kHLTObject_Ele32_CaloIdL_CaloIsoVL                             = 1UL<<10,
  kHLTObject_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT               = 1UL<<11,
 
  kHLTObject_Photon20_CaloIdVT_IsoT                              = 1UL<<12,
  kHLTObject_SC17                                                = 1UL<<13,
  kHLTObject_SC8                                                 = 1UL<<14,
  kHLTObject_Ele17                                               = 1UL<<15,

  kHLTObject_Ele8_CaloIdL                                        = 1UL<<16,
  kHLTObject_Ele17_CaloIdL                                       = 1UL<<17,
  kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT            = 1UL<<18,
  kHLTObject_Ele52_CaloIdVT_TrkIdT                               = 1UL<<19,

  kHLTObject_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT              = 1UL<<20,
  kHLTObject_Photon26_CaloIdL_IsoVL                              = 1UL<<21,  
  kHLTObject_Ele32_WP70                                          = 1UL<<22,
  kHLTObject_Ele32_WP80                                          = 1UL<<23,

  kHLTObject_IsoMu12                                             = 1UL<<24,
  kHLTObject_IsoMu17                                             = 1UL<<25,
  kHLTObject_IsoMu24                                             = 1UL<<26,  
  kHLTObject_IsoMu30                                             = 1UL<<27,

  kHLTObject_Mu7                                                 = 1UL<<28,
  kHLTObject_Mu17                                                = 1UL<<29,
  kHLTObject_Mu13                                                = 1UL<<30,
  kHLTObject_Mu3                                                 = 1UL<<31,

  kHLTObject_Jet30                                               = 1UL<<32,
  kHLTObject_Jet60                                               = 1UL<<33,
  kHLTObject_Jet80                                               = 1UL<<34,
  kHLTObject_Jet110                                              = 1UL<<35,

  kHLTObject_IsoMu34                                             = 1UL<<36,
  kHLTObject_IsoMu40                                             = 1UL<<37,

  kHLTObject_Mu5                                                 = 1UL<<40,
  kHLTObject_Mu40                                                = 1UL<<41,
  kHLTObject_TkMu8                                               = 1UL<<42, 

  kHLTObject_Ele65_CaloIdVT_TrkIdT                               = 1UL<<44,
  kHLTObject_Ele80_CaloIdVT_TrkIdT                               = 1UL<<45,
  kHLTObject_Ele27_WP80                                          = 1UL<<46,
  kHLTObject_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT            = 1UL<<47,

  kHLTObject_Ele8_CaloIdT_CaloIsoVL                              = 1UL<<48,
  kHLTObject_Ele17_CaloIdT_CaloIsoVL                             = 1UL<<49

};






enum L1TriggerBit
{  
  kL1_SingleEG5                                            = 1UL<<0,
  kL1_SingleEG12                                           = 1UL<<1,
  kL1_SingleEG20                                           = 1UL<<2,
  kL1_SingleEG30                                           = 1UL<<3,
  
  kL1_SingleMu3                                            = 1UL<<4,
  kL1_SingleMu7                                            = 1UL<<5,
  kL1_SingleMu10                                           = 1UL<<6,
  kL1_SingleMu12                                           = 1UL<<7,

  kL1_MuOpen_EG12                                          = 1UL<<8,
  kL1_MuOpen_EG5                                           = 1UL<<9,
  kL1_Mu5_EG12                                             = 1UL<<10,
  kL1_Mu7_EG5                                              = 1UL<<11,
 
  kL1_Mu12_EG5                                             = 1UL<<12,
  kL1_DoubleMu0                                            = 1UL<<13,
  kL1_DoubleMu5                                            = 1UL<<14,
  kL1_Mu3_EG5                                              = 1UL<<15,

  kL1_DoubleEG3                                            = 1UL<<16,
  kL1_DoubleEG5                                            = 1UL<<17,
  kL1_DoubleEG_12_5                                        = 1UL<<18,


  kL1_Mu3_Jet20_Central                                    = 1UL<<20,
  kL1_EG5_Jet36_deltaPhi1                                  = 1UL<<21,	     


  kL1_SingleMu20                                           = 1UL<<24,
  kL1_SingleMu25                                           = 1UL<<25


};

enum EPFType { 
  eX = 0,          //unidentified
  eHadron,         //charged hadron
  eElectron,       //electron
  eMuon,           //muon
  eGamma,          //photon
  eNeutralHadron,  //neutral hadron
  eHadronHF,       //hadron in HF
  eEGammaHF        //EM object in HF
};

#endif

