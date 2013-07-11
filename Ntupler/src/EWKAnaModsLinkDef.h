#ifndef HIGGSANA_NTUPLER_LINKDEF_H
#define HIGGSANA_NTUPLER_LINKDEF_H
#include "HiggsAna/Ntupler/interface/HwwNtuplerMod.hh"
#include "HiggsAna/Ntupler/interface/BambuGenDumperMod.hh"
#include "HiggsAna/Ntupler/interface/TEventInfo.hh"
#include "HiggsAna/Ntupler/interface/TGenInfo.hh"
#include "HiggsAna/Ntupler/interface/TGenParticle.hh"
#include "HiggsAna/Ntupler/interface/TMuon.hh"
#include "HiggsAna/Ntupler/interface/TElectron.hh"
#include "HiggsAna/Ntupler/interface/TJet.hh"
#include "HiggsAna/Ntupler/interface/TPhoton.hh"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::HwwNtuplerMod+;
#pragma link C++ class mithep::BambuGenDumperMod+;
#pragma link C++ class mithep::TEventInfo+;
#pragma link C++ class mithep::TGenInfo+;
#pragma link C++ class mithep::TMuon+;
#pragma link C++ class mithep::TElectron+;
#pragma link C++ class mithep::TJet+;
#pragma link C++ class mithep::TPhoton+;
#pragma link C++ class mithep::TGenParticle+;
#pragma link C++ class mithep::TPFCandidate+;
#endif
