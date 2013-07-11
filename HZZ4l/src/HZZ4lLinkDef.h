#ifndef HIGGSANA_HZZ4L_LINKDEF_H
#define HIGGSANA_HZZ4L_LINKDEF_H
#include "HiggsAna/HZZ4l/interface/HZZKinematics.hh"
#include "HiggsAna/HZZ4l/interface/HZZGenInfo.hh"
#include "HiggsAna/HZZ4l/interface/HZZEfficiencyMap.hh"
#include "HiggsAna/HZZ4l/interface/RooDoubleSidedCBShape.hh"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ class HZZKinematics+;
#pragma link C++ class HZZGenInfo+;
#pragma link C++ class HZZEfficiencyMap+;
#pragma link C++ class RooDoubleSidedCBShape+;
#endif
