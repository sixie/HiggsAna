#ifndef HIGGSANA_NTUPLER_LINKDEF_H
#define HIGGSANA_NTUPLER_LINKDEF_H
#include "HiggsAna/DataTree/interface/TEventInfo.hh"
#include "HiggsAna/DataTree/interface/TGenInfo.hh"
#include "HiggsAna/DataTree/interface/TGenParticle.hh"
#include "HiggsAna/DataTree/interface/TMuon.hh"
#include "HiggsAna/DataTree/interface/TElectron.hh"
#include "HiggsAna/DataTree/interface/TJet.hh"
#include "HiggsAna/DataTree/interface/TCaloJet.hh"
#include "HiggsAna/DataTree/interface/TPhoton.hh"
#include "HiggsAna/DataTree/interface/TPFCandidate.hh"
#include "HiggsAna/DataTree/interface/Types.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace higgsana;

#pragma link C++ class higgsana::TEventInfo+;
#pragma link C++ class higgsana::TGenInfo+;
#pragma link C++ class higgsana::TMuon+;
#pragma link C++ class higgsana::TElectron+;
#pragma link C++ class higgsana::TJet+;
#pragma link C++ class higgsana::TCaloJet+;
#pragma link C++ class higgsana::TPhoton+;
#pragma link C++ class higgsana::TGenParticle+;
#pragma link C++ class higgsana::TPFCandidate+;
#pragma link C++ typedef higgsana::FourVector;
#pragma link C++ typedef higgsana::FourVectorM;
#pragma link C++ typedef higgsana::FourVectorE;
#endif
