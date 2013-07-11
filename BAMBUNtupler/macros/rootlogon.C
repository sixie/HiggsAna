// $Id: rootlogon.C,v 1.2 2013/04/18 07:51:57 sixie Exp $

{

  if (gSystem->Getenv("CMSSW_VERSION")) {

    TString rfitpath("/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms5/include");
    TString path = gSystem->GetIncludePath();
    path += "-I. -I$ROOTSYS/src -I";
    path += rfitpath;
    gSystem->SetIncludePath(path.Data());
 
    TString str = gSystem->GetMakeSharedLib();
    if (str.Contains("-m32")==0 && str.Contains("-m64")==0) {
      str.ReplaceAll("g++", "g++ -m32");
      gSystem->SetMakeSharedLib(str);
    }
  }

  gROOT->Macro("$CMSSW_BASE/src/MitAna/macros/setRootEnv.C+");
//   gROOT->Macro("MitStyle.C+");
//   gROOT->Macro("Smurf/Core/LeptonScaleLookup.cc+");
//  gROOT->Macro("Smurf/HZZ/goodrun.cc+");

  loadmylib("MitCommon","DataFormats");
  loadmylib("MitPhysics","Mods");
  loadmylib("MitPhysics","Utils");
  loadmylib("MitHiggs","Utils");
  loadmylib("MitHiggs","HwwMods");
  loadmylib("HiggsAna","BAMBUNtupler");
  loadmylib("HiggsAna","HZZ4l");
  loadmylib("EGamma","EGammaAnalysisTools");
  loadmylib("MitHZGamma","Mods");

  gSystem->Load("libFWCoreFWLite.so"); AutoLibraryLoader::enable();
  gSystem->Load("libDataFormatsFWLite.so");

  //gROOT->Macro("HiggsAna/HZZ4l/HLL/RooDoubleSidedCBShape.cc+");


}
