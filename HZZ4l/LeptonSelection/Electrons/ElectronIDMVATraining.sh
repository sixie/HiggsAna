#####################################################################################
#Do pt reweighting
#####################################################################################
root -l -b -q EWKAna/Hww/LeptonSelection/MakeElectronPtSpectrum.C+


#####################################################################################
#Do PU reweighting
#####################################################################################
root -l -b -q $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/ReweightElectronPU.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Summer12HZZ.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WPlusJet_2012.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WPlusJet_2012.PtAndPUWeighted.root")'



#####################################################################################
#Split into training and testing ntuples
#####################################################################################

root -l -b -q $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/SplitElectronNtuples.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Summer12HZZ.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Summer12HZZ.Training.root",0)'
root -l -b -q $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/SplitElectronNtuples.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Summer12HZZ.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Summer12HZZ.Testing.root",1)'
root -l -b -q $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/SplitElectronNtuples.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_MC.AllNtuple_HZZ4lNtuple_s12-h125zz4l-gf-v9_noskim_0000.root.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Summer12HZZ125.Training.root",0)'
root -l -b -q $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/SplitElectronNtuples.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_MC.AllNtuple_HZZ4lNtuple_s12-h125zz4l-gf-v9_noskim_0000.root.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Summer12HZZ125.Testing.root",1)'
root -l -b -q $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/SplitElectronNtuples.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WPlusJet_2012.PtAndPUWeighted.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WPlusJet_2012.PtAndPUWeighted.Training.root",0)'
root -l -b -q $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/SplitElectronNtuples.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WPlusJet_2012.PtAndPUWeighted.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WPlusJet_2012.PtAndPUWeighted.Testing.root",1)'
root -l -b -q $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/SplitElectronNtuples.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_ZPlusJet_2012.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_ZPlusJet_2012.Training.root",0)'
root -l -b -q $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/SplitElectronNtuples.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_ZPlusJet_2012.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_ZPlusJet_2012.Testing.root",1)'


#####################################################################################
#MVA Training 
#####################################################################################

#Train in bins
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDMVA.C+'("ElectronIDMVA_V1_CentralPt5To10",0,"V1")'
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDMVA.C+'("ElectronIDMVA_V1_TransitionPt5To10",1,"V1")'
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDMVA.C+'("ElectronIDMVA_V1_EndcapPt5To10",2,"V1")'
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDMVA.C+'("ElectronIDMVA_V1_CentralPt10ToInf",3,"V1")'
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDMVA.C+'("ElectronIDMVA_V1_TransitionPt10ToInf",4,"V1")'
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDMVA.C+'("ElectronIDMVA_V1_EndcapPt10ToInf",5,"V1")'

root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDMVA.C+'("ElectronIDMVA_V2_CentralPt5To10",0,"V2")'
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDMVA.C+'("ElectronIDMVA_V2_TransitionPt5To10",1,"V2")'
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDMVA.C+'("ElectronIDMVA_V2_EndcapPt5To10",2,"V2")'
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDMVA.C+'("ElectronIDMVA_V2_CentralPt10ToInf",3,"V2")'
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDMVA.C+'("ElectronIDMVA_V2_TransitionPt10ToInf",4,"V2")'
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/TrainElectronIDMVA.C+'("ElectronIDMVA_V2_EndcapPt10ToInf",5,"V2")'



#####################################################################################
#MVA Evaluate
#####################################################################################

#FAKES
#root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDMVA.C+'("ElectronSelectionTraining.Real.Testing.root","output/ElectronSelectionTraining.Real.Testing.root","ElectronIDMVA")'

#HZZ Signal
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDMVA.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Summer12HZZ.Training.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Summer12HZZ.Training.root","ElectronIDMVA")'
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDMVA.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Summer12HZZ.Testing.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Summer12HZZ.Testing.root","ElectronIDMVA")'
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDMVA.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Summer12HZZ125.Training.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Summer12HZZ125.Training.root","ElectronIDMVA")'
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDMVA.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Summer12HZZ125.Testing.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Summer12HZZ125.Testing.root","ElectronIDMVA")'

#Z+Jets Data
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDMVA.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_ZPlusJet_2012.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Fake_ZPlusJet_2012.root","ElectronIDMVA")'
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDMVA.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_ZPlusJet_2012.Testing.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Fake_ZPlusJet_2012.Testing.root","ElectronIDMVA")'

#Z+Jets MC
#root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/EvaluateElectronIDMVA.C+'("ElectronSelectionTraining.Fall11ZJetsBkg.root","output/ElectronSelectionTraining.Fall11ZJetsBkg.root","ElectronIDMVA")'



#####################################################################################
#Performance Plots
#####################################################################################
root -l  $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDMVAPerformancePlots.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Summer12HZZ.Testing.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Fake_ZPlusJet_2012.root","EtaBin0PtBin0",0)'
root -l  $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDMVAPerformancePlots.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Summer12HZZ.Testing.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Fake_ZPlusJet_2012.Testing.root","EtaBin1PtBin0",1)'
root -l  $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDMVAPerformancePlots.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Summer12HZZ.Testing.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Fake_ZPlusJet_2012.Testing.root","EtaBin2PtBin0",2)'
root -l  $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDMVAPerformancePlots.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Summer12HZZ.Testing.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Fake_ZPlusJet_2012.Testing.root","EtaBin0PtBin1",3)'
root -l  $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDMVAPerformancePlots.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Summer12HZZ.Testing.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Fake_ZPlusJet_2012.Testing.root","EtaBin1PtBin1",4)'
root -l  $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDMVAPerformancePlots.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Summer12HZZ.Testing.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Fake_ZPlusJet_2012.Testing.root","EtaBin2PtBin1",5)'


root -l  $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDMVAPerformancePlots.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Summer12HZZ125.Testing.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Fake_ZPlusJet_2012.Testing.root","EtaBin0PtBin0",0)'
root -l  $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDMVAPerformancePlots.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Summer12HZZ125.Testing.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Fake_ZPlusJet_2012.Testing.root","EtaBin0PtBin1",3)'



#####################################################################################
#Performance Plots HZZ Vs ZJets MC
#####################################################################################
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDMVAPerformancePlots.C+'("ElectronSelectionTraining.HZZ.Testing.root","ElectronSelectionTraining.Fake_ZPlusJet.Testing.root","HZZVsZJetsData_BarrelPtBin0",0)'

#####################################################################################
#Performance Plots HZZ Vs ZJets Data
#####################################################################################
root -l -q -b $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeElectronIDMVAPerformancePlots.C+'("$CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.f11HZZ.Testing.root","$CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/output/ElectronSelectionTraining.Fake_ZPlusJet.root","HZZVsZJetsData_BarrelPtBin0",0)'





#####################################################################################
#Performance Plots with HWW115
#####################################################################################
root -l  $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIDMVAPerformancePlots.C+'("$CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDMVA/output/MuonSelection.Real.Testing.root","$CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDMVA/output/MuonSelection.Fake.Testing.root","BarrelPtBin0",0)'
root -l  $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIDMVAPerformancePlots.C+'("$CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDMVA/output/MuonSelection.Real.Testing.root","$CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDMVA/output/MuonSelection.Fake.Testing.root","EndcapPtBin0",1)'
root -l  $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIDMVAPerformancePlots.C+'("$CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDMVA/output/MuonSelection.Real.Testing.root","$CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDMVA/output/MuonSelection.Fake.Testing.root","BarrelPtBin1",2)'
root -l  $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIDMVAPerformancePlots.C+'("$CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDMVA/output/MuonSelection.Real.Testing.root","$CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDMVA/output/MuonSelection.Fake.Testing.root","EndcapPtBin1",3)'
root -l  $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIDMVAPerformancePlots.C+'("$CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDMVA/output/MuonSelection.Real.Testing.root","$CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDMVA/output/MuonSelection.Fake.Testing.root","BarrelPtBin2",4)'
root -l  $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/MakeMuonIDMVAPerformancePlots.C+'("$CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDMVA/output/MuonSelection.Real.Testing.root","$CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDMVA/output/MuonSelection.Fake.Testing.root","EndcapPtBin2",5)'







root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("$CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Muons/IDMVA/output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.Fake.EndcapPtBin0_V10.root","EndcapPtBin0_HWW115",1)'
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

root -l  $CMSSW_BASE/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/CompareElectronDistributions.C+'("/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WPlusJet_2012.PtAndPUWeighted.Training.root","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_ZPlusJet_2012.Testing.root","WPlusFake","ZPlusFake","EtaBin0PtBin1",3)'

