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
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA","V0")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA","V1")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA","V2")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA","V3")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA","V4")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA","V5")'

#Train in bins
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA_V0_BarrelPt5To10",0,"V0")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA_V0_EndcapPt5To10",1,"V0")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA_V0_BarrelPt10ToInf",2,"V0")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA_V0_EndcapPt10ToInf",3,"V0")'


root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA_V1_BarrelPt5To10",0,"V1")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA_V1_EndcapPt5To10",1,"V1")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA_V1_BarrelPt10ToInf",2,"V1")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA_V1_EndcapPt10ToInf",3,"V1")'


root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA_V2_BarrelPt5To10",0,"V2")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA_V2_EndcapPt5To10",1,"V2")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA_V2_BarrelPt10ToInf",2,"V2")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIsoMVA.C+'("ElectronIsoMVA_V2_EndcapPt10ToInf",3,"V2")'


#####################################################################################
#MVA Evaluate
#####################################################################################

#FAKES
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIsoMVA.C+'("ElectronSelectionTraining.Real.Testing.root","output/ElectronSelectionTraining.Real.Testing.root","ElectronIsoMVA")'

#HZZ Signal
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIsoMVA.C+'("ElectronSelectionTraining.f11HZZ.Testing.root","output/ElectronSelectionTraining.f11HZZ.Testing.root","ElectronIsoMVA")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIsoMVA.C+'("ElectronSelectionTraining.s11HZZ.Testing.root","output/ElectronSelectionTraining.s11HZZ.Testing.root","ElectronIsoMVA")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIsoMVA.C+'("ElectronSelectionTraining.HZZ.Testing.NoOverlapLeptons.root","output/ElectronSelectionTraining.HZZ.Testing.NoOverlapLeptons.root","ElectronIsoMVA")'

#Z+Jets Data
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIsoMVA.C+'("ElectronSelectionTraining.Fake_ZPlusJet.Testing.root","output/ElectronSelectionTraining.Fake_ZPlusJet.Testing.root","ElectronIsoMVA")'
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIsoMVA.C+'("ElectronSelectionTraining.Fake_ZPlusJet.root","output/ElectronSelectionTraining.Fake_ZPlusJet.root","ElectronIsoMVA")'

#Z+Jets MC
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIsoMVA.C+'("ElectronSelectionTraining.Fall11ZJetsBkg.root","output/ElectronSelectionTraining.Fall11ZJetsBkg.root","ElectronIsoMVA")'



#####################################################################################
#Performance Plots
#####################################################################################
root -l  /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/output/ElectronSelectionTraining.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/output/ElectronSelectionTraining.Fakes.Testing.root","BarrelPtBin0",0)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/output/ElectronSelectionTraining.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/output/ElectronSelectionTraining.Fakes.Testing.root","EndcapPtBin0",1)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/output/ElectronSelectionTraining.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/output/ElectronSelectionTraining.Fakes.Testing.root","BarrelPtBin1",2)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/output/ElectronSelectionTraining.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/output/ElectronSelectionTraining.Fakes.Testing.root","EndcapPtBin1",3)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/output/ElectronSelectionTraining.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/output/ElectronSelectionTraining.Fakes.Testing.root","BarrelPtBin2",4)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/output/ElectronSelectionTraining.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/output/ElectronSelectionTraining.Fakes.Testing.root","EndcapPtBin2",5)'



#####################################################################################
#Performance Plots HZZ Vs ZJets MC
#####################################################################################
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIsoMVAPerformancePlots.C+'("ElectronSelectionTraining.HZZ.Testing.root","ElectronSelectionTraining.Fake_ZPlusJet.Testing.root","HZZVsZJetsData_BarrelPtBin0",0)'

#####################################################################################
#Performance Plots HZZ Vs ZJets Data
#####################################################################################
root -l -q -b /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/output/ElectronSelectionTraining.f11HZZ.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IsoMVA/output/ElectronSelectionTraining.Fake_ZPlusJet.Testing.root","HZZVsZJetsData_BarrelPtBin0",0)'





#####################################################################################
#Performance Plots with HWW115
#####################################################################################
root -l  /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IsoMVA/output/MuonSelection.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IsoMVA/output/MuonSelection.Fake.Testing.root","BarrelPtBin0",0)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IsoMVA/output/MuonSelection.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IsoMVA/output/MuonSelection.Fake.Testing.root","EndcapPtBin0",1)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IsoMVA/output/MuonSelection.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IsoMVA/output/MuonSelection.Fake.Testing.root","BarrelPtBin1",2)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IsoMVA/output/MuonSelection.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IsoMVA/output/MuonSelection.Fake.Testing.root","EndcapPtBin1",3)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IsoMVA/output/MuonSelection.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IsoMVA/output/MuonSelection.Fake.Testing.root","BarrelPtBin2",4)'
root -l  /data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIsoMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IsoMVA/output/MuonSelection.Real.Testing.root","/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IsoMVA/output/MuonSelection.Fake.Testing.root","EndcapPtBin2",5)'







root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("/data/blue/sixie/releases/analysis/CMSSW_4_4_2/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IsoMVA/output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.Fake.EndcapPtBin0_V10.root","EndcapPtBin0_HWW115",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.Fake.BarrelPtBin1_V10.root","BarrelPtBin1_HWW115",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.Fake.EndcapPtBin1_V10.root","EndcapPtBin1_HWW115",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.Fake.BarrelPtBin2_V10.root","BarrelPtBin2_HWW115",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.Fake.EndcapPtBin2_V10.root","EndcapPtBin2_HWW115",5)'



root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW130_V10.root","output/MuonNtuple.Fake.BarrelPtBin0_V10.root","BarrelPtBin0_HWW130",0)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW130_V10.root","output/MuonNtuple.Fake.EndcapPtBin0_V10.root","EndcapPtBin0_HWW130",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW130_V10.root","output/MuonNtuple.Fake.BarrelPtBin1_V10.root","BarrelPtBin1_HWW130",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW130_V10.root","output/MuonNtuple.Fake.EndcapPtBin1_V10.root","EndcapPtBin1_HWW130",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW130_V10.root","output/MuonNtuple.Fake.BarrelPtBin2_V10.root","BarrelPtBin2_HWW130",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW130_V10.root","output/MuonNtuple.Fake.EndcapPtBin2_V10.root","EndcapPtBin2_HWW130",5)'


#####################################################################################
#Performance Plots with Fall11Zmm
#####################################################################################
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Fall11Zmm_V10.root","output/MuonNtuple.Fake.BarrelPtBin0_V10.root","BarrelPtBin0_Fall11Zmm",0)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Fall11Zmm_V10.root","output/MuonNtuple.Fake.EndcapPtBin0_V10.root","EndcapPtBin0_Fall11Zmm",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Fall11Zmm_V10.root","output/MuonNtuple.Fake.BarrelPtBin1_V10.root","BarrelPtBin1_Fall11Zmm",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Fall11Zmm_V10.root","output/MuonNtuple.Fake.EndcapPtBin1_V10.root","EndcapPtBin1_Fall11Zmm",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Fall11Zmm_V10.root","output/MuonNtuple.Fake.BarrelPtBin2_V10.root","BarrelPtBin2_Fall11Zmm",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Fall11Zmm_V10.root","output/MuonNtuple.Fake.EndcapPtBin2_V10.root","EndcapPtBin2_Fall11Zmm",5)'



#####################################################################################
#Performance Plots HWW115 Vs WJets MC
#####################################################################################
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Summer11Skimmed_V10.root","BarrelPtBin0_HWW115Summer11WJets",0)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Summer11Skimmed_V10.root","EndcapPtBin0_HWW115Summer11WJets",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Summer11Skimmed_V10.root","BarrelPtBin1_HWW115Summer11WJets",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Summer11Skimmed_V10.root","EndcapPtBin1_HWW115Summer11WJets",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Summer11Skimmed_V10.root","BarrelPtBin2_HWW115Summer11WJets",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Summer11Skimmed_V10.root","EndcapPtBin2_HWW115Summer11WJets",5)'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Winter10Skimmed_V10.root","BarrelPtBin0_HWW115Winter10WJets",0)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Winter10Skimmed_V10.root","EndcapPtBin0_HWW115Winter10WJets",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Winter10Skimmed_V10.root","BarrelPtBin1_HWW115Winter10WJets",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Winter10Skimmed_V10.root","EndcapPtBin1_HWW115Winter10WJets",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Winter10Skimmed_V10.root","BarrelPtBin2_HWW115Winter10WJets",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Winter10Skimmed_V10.root","EndcapPtBin2_HWW115Winter10WJets",5)'


#####################################################################################
#Comparison Plots
#####################################################################################


root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet0LowPt.root","output/ElectronNtuple.HWW115.Subdet0LowPt.root","Data","HWW115","Subdet0LowPt_Real",0)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet1LowPt.root","output/ElectronNtuple.HWW115.Subdet1LowPt.root","Data","HWW115","Subdet1LowPt_Real",1)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet2LowPt.root","output/ElectronNtuple.HWW115.Subdet2LowPt.root","Data","HWW115","Subdet2LowPt_Real",2)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet0HighPt.root","output/ElectronNtuple.HWW115.Subdet0HighPt.root","Data","HWW115","Subdet0HighPt_Real",3)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet1HighPt.root","output/ElectronNtuple.HWW115.Subdet1HighPt.root","Data","HWW115","Subdet1HighPt_Real",4)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet2HighPt.root","output/ElectronNtuple.HWW115.Subdet2HighPt.root","Data","HWW115","Subdet2HighPt_Real",5)'



root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtupleC.Real.Subdet0LowPt.root","output/ElectronNtupleC.HWW115.Subdet0LowPt.root","Data","HWW115","Subdet0LowPt_C_Real",0)'
