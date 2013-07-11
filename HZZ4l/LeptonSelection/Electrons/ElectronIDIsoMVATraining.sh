#####################################################################################
#Make Muon Ntuples
#####################################################################################
root -l -b -q HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeRealElectronTrainingNtuple.C+
root -l -b -q HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeFakeElectronTrainingNtuple.C+
root -l -b -q HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeMCElectronTrainingNtuple.C+

#####################################################################################
#Do pt reweighting
#####################################################################################
root -l -b -q EWKAna/Hww/LeptonSelection/MakeElectronPtSpectrum.C+
root -l -b -q EWKAna/Hww/LeptonSelection/ReweightElectronPU.C+



#####################################################################################
#MVA Training 
#####################################################################################

#Add Variables
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA","V0")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA","V1")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA","V2")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA","V3")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA","V4")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA","V5")'

#wjets trained
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V0_CentralPt5To10",0,"V0")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V0_TransitionPt5To10",1,"V0")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V0_EndcapPt5To10",2,"V0")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V0_CentralPt10ToInf",3,"V0")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V0_TransitionPt10ToInf",4,"V0")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V0_EndcapPt10ToInf",5,"V0")'

#Z+jets trained
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V1_CentralPt5To10",0,"V1")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V1_TransitionPt5To10",1,"V1")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V1_EndcapPt5To10",2,"V1")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V1_CentralPt10ToInf",3,"V1")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V1_TransitionPt10ToInf",4,"V1")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V1_EndcapPt10ToInf",5,"V1")'

#z+jets MC trained
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V2_CentralPt5To10",0,"V2")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V2_TransitionPt5To10",1,"V2")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V2_EndcapPt5To10",2,"V2")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V2_CentralPt10ToInf",3,"V2")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V2_TransitionPt10ToInf",4,"V2")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_V2_EndcapPt10ToInf",5,"V2")'


root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDIsoMVA.C+'("ElectronIDIsoMVA_Test_CentralPt5To10",0,"V0")'


#####################################################################################
#MVA Evaluate
#####################################################################################

#FAKES
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDIsoMVA.C+'("ElectronSelectionTraining.Real.Testing.root","output/ElectronSelectionTraining.Real.Testing.root","ElectronIDIsoMVA")'

#HZZ Signal
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDIsoMVA.C+'("ElectronSelectionTraining.f11HZZ.Testing.root","output/ElectronSelectionTraining.f11HZZ.Testing.root","ElectronIDIsoMVA")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDIsoMVA.C+'("ElectronSelectionTraining.s11HZZ.Testing.root","output/ElectronSelectionTraining.s11HZZ.Testing.root","ElectronIDIsoMVA")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDIsoMVA.C+'("ElectronSelectionTraining.HZZ.Testing.NoOverlapLeptons.root","output/ElectronSelectionTraining.HZZ.Testing.NoOverlapLeptons.root","ElectronIDIsoMVA")'

#Z+Jets Data
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDMVA.C+'("ElectronSelectionTraining.Fake_ZPlusJet_2012.root","output/ElectronSelectionTraining.Fake_ZPlusJet_2012.root","ElectronIDMVA")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIsoMVA.C+'("ElectronSelectionTraining.Fake_ZPlusJet_2012.root","output/ElectronSelectionTraining.Fake_ZPlusJet_2012.root","ElectronIsoMVA")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDIsoMVA.C+'("ElectronSelectionTraining.Fake_ZPlusJet_2012.root","output/ElectronSelectionTraining.Fake_ZPlusJet_2012.root","ElectronIDIsoMVA")'


#Z+Jets MC
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDMVA.C+'("ElectronSelectionTraining.s12ZJets52X.root","output/ElectronSelectionTraining.s12ZJets52X.root","ElectronIDMVA")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIsoMVA.C+'("ElectronSelectionTraining.s12ZJets52X.root","output/ElectronSelectionTraining.s12ZJets52X.root","ElectronIsoMVA")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDIsoMVA.C+'("ElectronSelectionTraining.s12ZJets52X.root","output/ElectronSelectionTraining.s12ZJets52X.root","ElectronIDIsoMVA")'

root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDMVA.C+'("ElectronSelectionTraining.s12ZJets51X.root","output/ElectronSelectionTraining.s12ZJets51X.root","ElectronIDMVA")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIsoMVA.C+'("ElectronSelectionTraining.s12ZJets51X.root","output/ElectronSelectionTraining.s12ZJets51X.root","ElectronIsoMVA")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDIsoMVA.C+'("ElectronSelectionTraining.s12ZJets51X.root","output/ElectronSelectionTraining.s12ZJets51X.root","ElectronIDIsoMVA")'



#####################################################################################
#Performance Plots
#####################################################################################
root -l  /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.Fakes.Testing.root","BarrelPtBin0",0)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.Fakes.Testing.root","EndcapPtBin0",1)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.Fakes.Testing.root","BarrelPtBin1",2)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.Fakes.Testing.root","EndcapPtBin1",3)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.Fakes.Testing.root","BarrelPtBin2",4)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.Fakes.Testing.root","EndcapPtBin2",5)'



#####################################################################################
#Performance Plots HZZ Vs ZJets MC
#####################################################################################
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDIsoMVAPerformancePlots.C+'("ElectronSelectionTraining.HZZ.Testing.root","ElectronSelectionTraining.Fake_ZPlusJet.Testing.root","HZZVsZJetsData_BarrelPtBin0",0)'

#####################################################################################
#Performance Plots HZZ Vs ZJets Data
#####################################################################################
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.f11HZZ.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.Fake_ZPlusJet.Testing.root","HZZVsZJetsData_BarrelPtBin0",0)'



#####################################################################################
#Performance Plots ZeeMC Vs Z+Fake Data
#####################################################################################
root -l  /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.s12ZJets52X.root","/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.Fake_ZPlusJet_2012.root","HZZVsZJetsData_BarrelPtBin0",0)'

root -l  /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.s12ZJets51X.root","/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDIsoMVA/output/ElectronSelectionTraining.Fake_ZPlusJet_2012.root","HZZVsZJetsData_BarrelPtBin0",0)'


#####################################################################################
#Performance Plots with HWW115
#####################################################################################
root -l  /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIDIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDIsoMVA/output/MuonSelection.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDIsoMVA/output/MuonSelection.Fake.Testing.root","BarrelPtBin0",0)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIDIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDIsoMVA/output/MuonSelection.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDIsoMVA/output/MuonSelection.Fake.Testing.root","EndcapPtBin0",1)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIDIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDIsoMVA/output/MuonSelection.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDIsoMVA/output/MuonSelection.Fake.Testing.root","BarrelPtBin1",2)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIDIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDIsoMVA/output/MuonSelection.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDIsoMVA/output/MuonSelection.Fake.Testing.root","EndcapPtBin1",3)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIDIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDIsoMVA/output/MuonSelection.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDIsoMVA/output/MuonSelection.Fake.Testing.root","BarrelPtBin2",4)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIDIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDIsoMVA/output/MuonSelection.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDIsoMVA/output/MuonSelection.Fake.Testing.root","EndcapPtBin2",5)'





