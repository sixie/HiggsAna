#--------------------------------------------------------------
#  Electron Efficiencies
#==============================================================


#--------------------------------------------------------------
#  Electron Triggers
#==============================================================

####################################
# 
# Double Electron Triggers
#
# HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL : 150000 -> 170053
# HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL : 150000 -> 170053
# HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL  170054 -> 999999
# 
#
# root -l -q selectDoubleElectronLeadingLegTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_Full2011\",0\)
# root -l -q selectDoubleElectronLeadingLegTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_150000-170053\",1\)
# root -l -q selectDoubleElectronLeadingLegTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_170054-999999\",2\)
#
# root -l -q selectDoubleElectronTrailingLegTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_Full2011\",0\)
# root -l -q selectDoubleElectronTrailingLegTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_150000-170053\",1\)
# root -l -q selectDoubleElectronTrailingLegTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_170054-999999\",2\)
#
#root -l -b -q plotEff.C+\(\"elDoubleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_Full2011/probes.root\",\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_Full2011/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elDoubleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_150000-170053/probes.root\",\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_150000-170053/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elDoubleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_170054-999999/probes.root\",\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_170054-999999/basic_76_106\",\"all\",1,0,\"\"\)
#
#root -l -b -q plotEff.C+\(\"elDoubleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_Full2011/probes.root\",\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_Full2011/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elDoubleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_150000-170053/probes.root\",\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_150000-170053/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elDoubleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_170054-999999/probes.root\",\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_170054-999999/basic_76_106\",\"all\",1,0,\"\"\)
#
# root -l -q makeTriggerEfficiency.C+'("Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_Full2011/basic_76_106/eff.root","Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_Full2011/basic_76_106/","h2_results_electron_double_leadingleg","ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_Full2011")'
# root -l -q makeTriggerEfficiency.C+'("Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_Full2011/basic_76_106/eff.root","Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_Full2011/basic_76_106/","h2_results_electron_double_trailingleg","ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_Full2011")'
#
#####################################
####################################
# 
# MuEG Triggers : Electron Leg
#
# HLT_Mu17_Ele8_CaloIdL : 150000 -> 173198
# HLT_Mu8_Ele17_CaloIdL : 150000 -> 170053
# HLT_Mu8_Ele17_CaloIdT_CaloIsoVL : 170054 -> 999999
# HLT_Mu17_Ele8_CaloIdT_CaloIsoVL : 173199 -> 999999
# 
#
# root -l -q selectMuEGElectronLeadingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_Full2011\",0\)
# root -l -q selectMuEGElectronLeadingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_150000-170053\",1\)
# root -l -q selectMuEGElectronLeadingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_170054-173199\",2\)
# root -l -q selectMuEGElectronLeadingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_173200-999999\",3\)
#
# root -l -q selectMuEGElectronTrailingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_Full2011\",0\)
# root -l -q selectMuEGElectronTrailingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_150000-170053\",1\)
# root -l -q selectMuEGElectronTrailingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_170054-173199\",2\)
# root -l -q selectMuEGElectronTrailingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_173200-999999\",3\)
#
# root -l -q selectMuEGElectronTrailingLegTrigEffTP.C+\(\"data_sel_Test.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_Test\",3\)
#
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_Full2011/probes.root\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_Full2011/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_150000-170053/probes.root\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_150000-170053/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_170054-173199/probes.root\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_170054-173199/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_173200-999999/probes.root\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_173200-999999/basic_76_106\",\"all\",1,0,\"\"\)
#
#
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_Full2011/probes.root\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_Full2011/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_150000-170053/probes.root\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_150000-170053/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_170054-173199/probes.root\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_170054-173199/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_173200-999999/probes.root\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_173200-999999/basic_76_106\",\"all\",1,0,\"\"\)
#
#
#
# root -l -q makeTriggerEfficiency.C+'("Data_ElectronMVAIDIsoCombined_DoubleEleSeeded_Full2011/basic_76_106/eff.root","Data_ElectronMVAIDIsoCombined_DoubleEleSeeded_Full2011/basic_76_106/","h2_results_electron_double","ElectronMVAIDIsoCombined_DoubleEle_Full2011")'
#
#####################################
#
# Single Electron Triggers
# 
# HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT : 150000 -> 164237
# HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT : 165085 -> 166967
# HLT_Ele52_CaloIdVT_TrkIdT                  : 166968 -> 170053
# HLT_Ele65_CaloIdVT_TrkIdT                  : 170054 -> 178380
# HLT_Ele80_CaloIdVT_TrkIdT                  : 178381 -> 999999
#
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_Full2011\",0\)
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-164237\",1\)
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_165085-166967\",2\)
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_166968-170053\",3\)
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_170054-178380\",4\)
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_178381-999999\",5\)
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-166967\",10\)
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-170053\",11\)
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-178380\",12\)
#
#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_Full2011/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_Full2011/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-164237/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-164237/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_165085-166967/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_165085-166967/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_166968-170053/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_166968-170053/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_170054-178380/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_170054-178380/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_178381-999999/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_178381-999999/basic_76_106\",\"all\",1,0,\"\"\)

#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-166967/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-166967/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-170053/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-170053/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-178380/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-178380/basic_76_106\",\"all\",1,0,\"\"\)


#
# root -l -q makeTriggerEfficiency.C+'("Data_ElectronMVAIDIsoCombined_SingleEleSeeded_Full2011/basic_76_106/eff.root","Data_ElectronMVAIDIsoCombined_SingleEleSeeded_Full2011/basic_76_106/","h2_results_electron_single","ElectronMVAIDIsoCombined_SingleEle_Full2011")'
#
####################################


#--------------------------------------------------------------
#  Electron HWW ICHEP2012 WP
#==============================================================

####################################
# 
# Select Monte Carlo
#
# root -l -q selectEleHWWICHEP2012WP.C+\(\"s12-zee.conf\",\"results/Summer12_EleHWWICHEP2012WP\",1,2\)
#
####################################

####################################
# 
# Select 2012 data
#
#root -l -q selectEleHWWICHEP2012WP.C+\(\"data_el.2012.conf\",\"results/Data2012_EleHWWICHEP2012WP\"\)
#
####################################

####################################
# 
# Measure Efficiency in MC and Data
#
#
##########
#2012
#########
#root -l -b -q plotEff.C+\(\"el0.fine.bins\",0,0,0,0,\"results/Summer12_EleHWWICHEP2012WP/probes.root\",\"results/Summer12_EleHWWICHEP2012WP/basic_76_106_ReweightedToFull2011_FineBinning\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root\"\)
#root -l -b -q plotEff.C+\(\"el0.fine.bins\",2,1,2,3,\"results/Data2012_EleHWWICHEP2012WP/probes.root\",\"results/Data2012_EleHWWICHEP2012WP/basic2_76_106_FineBinning\",\"all\",1,0,\"results/Summer12_EleHWWICHEP2012WP/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root\"\)
#
#root -l -b -q plotEff.C+\(\"el0.bins\",0,0,0,0,\"results/Summer12_EleHWWICHEP2012WP/probes.root\",\"results/Summer12_EleHWWICHEP2012WP/basic_76_106_ReweightedToFull2011\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root\"\)
#root -l -b -q plotEff.C+\(\"el0.bins\",2,1,2,3,\"results/Data2012_EleHWWICHEP2012WP/probes.root\",\"results/Data2012_EleHWWICHEP2012WP/basic2_76_106\",\"all\",1,0,\"results/Summer12_EleHWWICHEP2012WP/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root\"\)
#
####################################

####################################
# 
# Compute Efficiency Scale Factors
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2012_EleHWWICHEP2012WP/basic2_76_106/eff.root","results/Summer12_EleHWWICHEP2012WP/basic_76_106_ReweightedToFull2011/eff.root","results/Data2012_EleHWWICHEP2012WP/basic2_76_106/","h2_results_electron_selection","EleHWWICHEP2012WP_Full2012")'

#
####################################







#--------------------------------------------------------------
#  Electron HWW IDIsoCombined MVA V3 WP
#==============================================================

####################################
# 
# Select Monte Carlo
#
# root -l -q selectEleHWWIDIsoMVAV3WP.C+\(\"s12-zee.conf\",\"results/Summer12_EleHWWIDIsoMVAV3WP\",1,2\)
# root -l -q /afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HWW/Efficiency/selectEleHWWIDIsoMVAV3WP.C+\(\"root://eoscms//eos/cms/store/user/sixie/hist/AllNtuple/cern/filefi/028/AllNtuple_HZZ4lNtuple_s12-zllm50-2-v9_noskim_0000.root\",\"/afs/cern.ch/work/s/sixie/public/HWW/Efficiency/Electrons/Summer12_EleHWWIDIsoMVAV3WP/probes.AllNtuple_HZZ4lNtuple_s12-zllm50-2-v9_noskim_0000.root\",1,-1,2\)
#
####################################

####################################
# 
# Select 2012 data
#
#root -l -q selectEleHWWIDIsoMVAV3WP.C+\(\"data_el.2012.conf\",\"results/Data2012_EleHWWIDIsoMVAV3WP\"\)
#
####################################

####################################
# 
# Measure Efficiency in MC and Data
#
#
##########
#2012
#########
#root -l -b -q plotEff.C+\(\"el0.fine.bins\",0,0,0,0,\"results/Summer12_EleHWWIDIsoMVAV3WP/probes.root\",\"results/Summer12_EleHWWIDIsoMVAV3WP/basic_76_106_ReweightedToFull2011_FineBinning\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root\"\)
#root -l -b -q plotEff.C+\(\"el0.fine.bins\",2,1,2,3,\"results/Data2012_EleHWWIDIsoMVAV3WP/probes.root\",\"results/Data2012_EleHWWIDIsoMVAV3WP/basic2_76_106_FineBinning\",\"all\",1,0,\"results/Summer12_EleHWWIDIsoMVAV3WP/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root\"\)
#
#root -l -b -q plotEff.C+\(\"el0.bins\",0,0,0,0,\"results/Summer12_EleHWWIDIsoMVAV3WP/probes.root\",\"results/Summer12_EleHWWIDIsoMVAV3WP/basic_76_106_ReweightedToFull2011\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root\"\)
#root -l -b -q plotEff.C+\(\"el0.bins\",2,1,2,3,\"results/Data2012_EleHWWIDIsoMVAV3WP/probes.root\",\"results/Data2012_EleHWWIDIsoMVAV3WP/basic2_76_106\",\"all\",1,0,\"results/Summer12_EleHWWIDIsoMVAV3WP/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root\"\)
#
####################################

####################################
# 
# Compute Efficiency Scale Factors
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2012_EleHWWIDIsoMVAV3WP/basic2_76_106/eff.root","results/Summer12_EleHWWIDIsoMVAV3WP/basic_76_106_ReweightedToFull2011/eff.root","results/Data2012_EleHWWIDIsoMVAV3WP/basic2_76_106/","h2_results_electron_selection","EleHWWIDIsoMVAV3WP_Full2012")'
#
####################################



#--------------------------------------------------------------
#  Electron HWW IDIsoCombined MVA V4 WP
#==============================================================

####################################
# 
# Select Monte Carlo
#
# root -l -q selectEleHWWIDIsoMVAV4WP.C+\(\"s12-zee.conf\",\"results/Summer12_EleHWWIDIsoMVAV4WP\",1,2\)
#
####################################

####################################
# 
# Select 2012 data
#
#root -l -q selectEleHWWIDIsoMVAV4WP.C+\(\"data_el.2012.conf\",\"results/Data2012_EleHWWIDIsoMVAV4WP\"\)
#
####################################

####################################
# 
# Measure Efficiency in MC and Data
#
##########
#2012
#########
#root -l -b -q plotEff.C+\(\"el0.fine.bins\",0,0,0,0,\"results/Summer12_EleHWWIDIsoMVAV4WP/probes.root\",\"results/Summer12_EleHWWIDIsoMVAV4WP/basic_76_106_ReweightedToFull2011_FineBinning\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root\"\)
#root -l -b -q plotEff.C+\(\"el0.fine.bins\",2,1,2,3,\"results/Data2012_EleHWWIDIsoMVAV4WP/probes.root\",\"results/Data2012_EleHWWIDIsoMVAV4WP/basic2_76_106_FineBinning\",\"all\",1,0,\"results/Summer12_EleHWWIDIsoMVAV4WP/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root\"\)
#
#root -l -b -q plotEff.C+\(\"el0.bins\",0,0,0,0,\"results/Summer12_EleHWWIDIsoMVAV4WP/probes.root\",\"results/Summer12_EleHWWIDIsoMVAV4WP/basic_76_106_ReweightedToFull2011\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root\"\)
#root -l -b -q plotEff.C+\(\"el0.bins\",2,1,2,3,\"results/Data2012_EleHWWIDIsoMVAV4WP/probes.root\",\"results/Data2012_EleHWWIDIsoMVAV4WP/basic2_76_106\",\"all\",1,0,\"results/Summer12_EleHWWIDIsoMVAV4WP/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root\"\)
#
####################################

####################################
# 
# Compute Efficiency Scale Factors
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2012_EleHWWIDIsoMVAV4WP/basic2_76_106/eff.root","results/Summer12_EleHWWIDIsoMVAV4WP/basic_76_106_ReweightedToFull2011/eff.root","results/Data2012_EleHWWIDIsoMVAV4WP/basic2_76_106/","h2_results_electron_selection","EleHWWIDIsoMVAV4WP_Full2012")'
#
####################################


