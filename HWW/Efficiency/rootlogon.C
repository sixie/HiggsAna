{    
  if(gSystem->Getenv("CMSSW_VERSION")) {
    TString rfitpath("$ROOFITSYS/include");
    TString path = gSystem->GetIncludePath();
    path += "-I. -I$ROOTSYS/src -I";
    path += rfitpath;
    gSystem->SetIncludePath(path.Data());
    
    TString str = gSystem->GetMakeSharedLib();
    if (str.Contains("-m32")==0 && str.Contains("-m64")==0) {
      str.ReplaceAll("g++", "g++ -m32");
      gSystem->SetMakeSharedLib(str);
    }      
    
    gROOT->Macro("$CMSSW_BASE/src/MitAna/macros/setRootEnv.C+");
    
    gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc462/libMitCommonDataFormats.so");
    gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc462/libHiggsAnaDataTree.so");
    gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc462/libEGammaEGammaAnalysisTools.so");

    gROOT->Macro("$CMSSW_BASE/src/CITCommon/FitModels/RooVoigtianShape.cc+");
    gROOT->Macro("$CMSSW_BASE/src//CITCommon/FitModels/RooErf.cc+");
    gROOT->Macro("$CMSSW_BASE/src//CITCommon/FitModels/RooCMSShape.cc+");
    gROOT->Macro("$CMSSW_BASE/src/HiggsAna/Utils/CPlot.cc+");
    gROOT->Macro("$CMSSW_BASE/src/HiggsAna/Utils/PlotStyle.cc+"); 
    gROOT->Macro("$CMSSW_BASE/src/HiggsAna/Utils/CEffUser1D.cc+");
    gROOT->Macro("$CMSSW_BASE/src/HiggsAna/Utils/CEffUser2D.cc+"); 
    
  } 

// else {

//     TString mypath("/home/sixie/Include");
//     TString path = gSystem->GetIncludePath();
//     path += " -I";
//     path += mypath;
//     gSystem->SetIncludePath(path.Data());      
//     gROOT->Macro(mypath+TString("/CITCommon/FitModels/RooVoigtianShape.cc+"));
//     gROOT->Macro(mypath+TString("/CITCommon/FitModels/RooErf.cc+"));
//     gROOT->Macro(mypath+TString("/CITCommon/FitModels/RooCMSShape.cc+"));
// //    gROOT->Macro(mypath+TString("/CITCommon/FitModels/RooVoigtianShape.cc+"));
// //    gROOT->Macro(mypath+TString("/CITCommon/FitModels/RooCMSShape.cc+"));
    
//     gROOT->Macro(mypath+TString("/HiggsAna/Utils/CPlot.cc+"));
//     gROOT->Macro(mypath+TString("/HiggsAna/Utils/MitStyleRemix.cc+")); 
//     gROOT->Macro(mypath+TString("/HiggsAna/Utils/CEffUser1D.cc+"));
//     gROOT->Macro(mypath+TString("/HiggsAna/Utils/CEffUser2D.cc+")); 
//   }
               
  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
