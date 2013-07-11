#--------------------------------------------------------------
#  Electron Efficiencies
#==============================================================

#--------------------------------------------------------------
#  Electron HZZ Run1 Legacy Paper Tight Tag + Met cut
#==============================================================

####################################
# 
# Measure Efficiency in MC and Data
#
#
###################################


#######################################################
#Run1 Legacy Paper  2011 Data
#######################################################
##############
#High Pt bins
##############
#root -l -b -q plotEff.C+\(\"el0.bins\",0,0,0,0,\"results/Fall11_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"results/Fall11_EleHZZRun1LegacyPaperWPWithTightTag/basic_76_106_ReweightedToRun2011\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2011/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"el0.bins\",2,1,2,3,\"results/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"results/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106\",\"all\",1,0,\"results/Fall11_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2011/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)
##############
#Low Pt bins
##############
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Fall11_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"results/Fall11_EleHZZRun1LegacyPaperWPWithTightTag/basic_76_106_LowPt_ReweightedToRun2011\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2011/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",2,1,2,3,\"results/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"results/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_LowPt\",\"all\",1,0,\"results/Fall11_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2011/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",2,1,2,3,\"results/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/probes.sel.root\",\"results/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_7To10GeV\",\"all\",1,0,\"results/Fall11_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2011/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)




##############################################
#With Turn on effect modeled for background
##############################################
#root -l -b -q plotEff.C+\(\"el0.bins\",2,1,2,2,\"results/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"results/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_WithTurnOn\",\"all\",1,0,\"results/Fall11_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2011/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)





#######################################################
#Run1 Legacy Paper  2012 Data
######################################################
##############
#High Pt bins
##############
#root -l -b -q plotEff.C+\(\"el0.bins\",0,0,0,0,\"results/Summer12_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"results/Summer12_EleHZZRun1LegacyPaperWPWithTightTag/basic_76_106_ReweightedToRun1LegacyPaper\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#root -l -b -q plotEff.C+\(\"el0.bins\",2,1,2,3,\"results/Data2012_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"results/Data2012_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106\",\"all\",1,0,\"results/Summer12_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
##############
#Low Pt bins
##############
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Summer12_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"results/Summer12_EleHZZRun1LegacyPaperWPWithTightTag/basic_76_106_LowPt_ReweightedToRun1LegacyPaper\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",2,1,2,3,\"results/Data2012_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"results/Data2012_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_LowPt\",\"all\",1,0,\"results/Summer12_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)


##############################################
#With Turn on effect modeled for background
##############################################
#root -l -b -q plotEff.C+\(\"el0.bins\",2,1,2,2,\"results/Data2012_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"results/Data2012_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_WithTurnOn\",\"all\",1,0,\"results/Summer12_EleHZZRun1LegacyPaperWPWithTightTag/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)

#######################################################
#######################################################


####################################
# 
# Compute Efficiency Scale Factors
#
#####
#2011
#####
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106/eff.root","results/Fall11_EleHZZRun1LegacyPaperWPWithTightTag/basic_76_106_ReweightedToRun2011/eff.root","results/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106/","h2_results_electron_selection","EleHZZRun1LegacyPaperWPWithTightTag_Run1LegacyPaper_2011Data")'
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_LowPt/eff.root","results/Fall11_EleHZZRun1LegacyPaperWPWithTightTag/basic_76_106_LowPt_ReweightedToRun2011/eff.root","results/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_LowPt/","h2_results_electron_selection","EleHZZRun1LegacyPaperWPWithTightTag_LowPt_Run1LegacyPaper_2011Data")'
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_7To10GeV/eff.root","results/Fall11_EleHZZRun1LegacyPaperWPWithTightTag/basic_76_106_LowPt_ReweightedToRun2011/eff.root","results/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_7To10GeV/","h2_results_electron_selection","EleHZZRun1LegacyPaperWPWithTightTag_7To10GeV_Run1LegacyPaper_2011Data")'


#####
#2012
#####
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2012_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_WithTurnOn/eff.root","results/Summer12_EleHZZRun1LegacyPaperWPWithTightTag/basic_76_106_ReweightedToRun1LegacyPaper/eff.root","results/Data2012_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_WithTurnOn/","h2_results_electron_selection","EleHZZRun1LegacyPaperWPWithTightTag_Run1LegacyPaper_2012Data")'
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2012_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_LowPt_WithTurnOn/eff.root","results/Summer12_EleHZZRun1LegacyPaperWPWithTightTag/basic_76_106_LowPt_ReweightedToRun1LegacyPaper/eff.root","results/Data2012_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_LowPt_WithTurnOn/","h2_results_electron_selection","EleHZZRun1LegacyPaperWPWithTightTag_LowPt_Run1LegacyPaper_2012Data")'
####################################











#--------------------------------------------------------------
#  Electron HZZ Run1 Legacy Paper IDGivenIso Tight Tag + Met cut
#==============================================================
#######################################################
#Run1 Legacy Paper 2011 Data
######################################################
##############
#Low Pt bins
##############


#####
# Use BWxCB signal model
#####
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Fall11_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/probes.root\",\"results/Fall11_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic_76_106_LowPt_ReweightedToRun1LegacyPaper\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2011/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",2,1,1,3,\"results/Data2011_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/probes.sel.root\",\"results/Data2011_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic2_76_106_LowPt\",\"all\",1,0,\"results/Fall11_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2011/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)
#
# Compute Efficiency Scale Factors
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2011_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic2_76_106_LowPt/eff.root","results/Fall11_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic_76_106_LowPt_ReweightedToRun1LegacyPaper/eff.root","results/Data2011_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic2_76_106_LowPt/","h2_results_electron_selection","EleHZZRun1LegacyPaperIDGivenIsoWithTightTag_Run1LegacyPaper")'


#####
# Use MC template signal model
#####
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",2,1,2,3,\"results/Data2011_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/probes.sel.root\",\"results/Data2011_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic2_76_106_LowPt_MCTemplate\",\"all\",1,0,\"results/Fall11_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2011/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)
#
# Compute Efficiency Scale Factors
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2011_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic2_76_106_LowPt_MCTemplate/eff.root","results/Fall11_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic_76_106_LowPt_ReweightedToRun1LegacyPaper/eff.root","results/Data2011_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic2_76_106_LowPt_MCTemplate/","h2_results_electron_selection","EleHZZRun1LegacyPaperIDGivenIsoWithTightTag_Run1LegacyPaper_2011Data")'



#######################################################
#Run1 Legacy Paper 2012 Data
######################################################
##############
#Low Pt bins
##############


#####
# Use BWxCB signal model
#####
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Summer12_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/probes.root\",\"results/Summer12_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic_76_106_LowPt_ReweightedToRun1LegacyPaper\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",2,1,1,3,\"results/Data2012_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/probes.root\",\"results/Data2012_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic2_76_106_LowPt\",\"all\",1,0,\"results/Summer12_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#
# Compute Efficiency Scale Factors
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2012_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic2_76_106_LowPt/eff.root","results/Summer12_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic_76_106_LowPt_ReweightedToRun1LegacyPaper/eff.root","results/Data2012_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic2_76_106_LowPt/","h2_results_electron_selection","EleHZZRun1LegacyPaperIDGivenIsoWithTightTag_Run1LegacyPaper")'


#####
# Use MC template signal model
#####
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",2,1,2,3,\"results/Data2012_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/probes.root\",\"results/Data2012_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic2_76_106_LowPt_MCTemplate\",\"all\",1,0,\"results/Summer12_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#
# Compute Efficiency Scale Factors
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2012_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic2_76_106_LowPt_MCTemplate/eff.root","results/Summer12_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic_76_106_LowPt_ReweightedToRun1LegacyPaper/eff.root","results/Data2012_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/basic2_76_106_LowPt_MCTemplate/","h2_results_electron_selection","EleHZZRun1LegacyPaperIDGivenIsoWithTightTag_Run1LegacyPaper")'







#--------------------------------------------------------------
#  Electron HZZ HCP2012 IsoGivenID Tight Tag + Met cut
#==============================================================
#######################################################
#Run1 Legacy Paper 2011 Data
######################################################
##############
#Low Pt bins
##############
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Fall11_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/probes.root\",\"results/Fall11_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/basic_76_106_LowPt_ReweightedToRun1LegacyPaper\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2011/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)
####################
#Use BWxCB
####################
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",2,1,1,3,\"results/Data2011_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/probes.sel.root\",\"results/Data2012_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/basic2_76_106_LowPt\",\"all\",1,0,\"results/Fall11_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2011/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)
#
####################
#Use MC Template
####################
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",2,1,2,3,\"results/Data2011_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/probes.sel.root\",\"results/Data2011_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/basic2_76_106_LowPt_MCTemplate\",\"all\",1,0,\"results/Fall11_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2011/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)
#


# 
# Compute Efficiency Scale Factors
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2011_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/basic2_76_106_LowPt_MCTemplate/eff.root","results/Fall11_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/basic_76_106_LowPt_ReweightedToRun1LegacyPaper/eff.root","results/Data2011_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/basic2_76_106_LowPt_MCTemplate/","h2_results_electron_selection","EleHZZRun1LegacyPaperIsoGivenIDWithTightTag_Run1LegacyPaper_2011Data")'
#



#######################################################
#Run1 Legacy Paper 2012 Data
######################################################
##############
#Low Pt bins
##############
#
#root -l -b -q plotEff.C+\(\"el0.bins\",0,0,0,0,\"results/Summer12_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/probes.root\",\"results/Summer12_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/basic_76_106_ReweightedToRun1LegacyPaper\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#root -l -b -q plotEff.C+\(\"el0.bins\",2,1,2,3,\"results/Data2012_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/probes.root\",\"results/Data2012_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/basic2_76_106\",\"all\",1,0,\"results/Summer12_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Summer12_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/probes.root\",\"results/Summer12_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/basic_76_106_LowPt_ReweightedToRun1LegacyPaper\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
####################
#Use BWxCB
####################
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",2,1,1,3,\"results/Data2012_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/probes.root\",\"results/Data2012_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/basic2_76_106_LowPt\",\"all\",1,0,\"results/Summer12_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#
####################
#Use MC Template
####################
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",2,1,2,3,\"results/Data2012_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/probes.root\",\"results/Data2012_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/basic2_76_106_LowPt_MCTemplate\",\"all\",1,0,\"results/Summer12_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#


# 
# Compute Efficiency Scale Factors
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2012_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/basic2_76_106_LowPt_MCTemplate/eff.root","results/Summer12_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/basic_76_106_LowPt_ReweightedToRun1LegacyPaper/eff.root","results/Data2012_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/basic2_76_106_LowPt_MCTemplate/","h2_results_electron_selection","EleHZZRun1LegacyPaperIsoGivenIDWithTightTag_Run1LegacyPaper")'
#




#--------------------------------------------------------------
#  Electron HZZ Moriond2013 IDGivenTighterIso 
#==============================================================

#######################################################
#2012 HCP
######################################################
#
#root -l -b -q plotEff.C+\(\"el0.bins\",0,0,0,0,\"results/Summer12_EleHZZMoriond2013IDGivenTighterIso/probes.root\",\"results/Summer12_EleHZZMoriond2013IDGivenTighterIso/basic_76_106_ReweightedToMoriond2013\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#root -l -b -q plotEff.C+\(\"el0.bins\",2,1,2,3,\"results/Data2012_EleHZZMoriond2013IDGivenTighterIso/probes.root\",\"results/Data2012_EleHZZMoriond2013IDGivenTighterIso/basic2_76_106\",\"all\",1,0,\"results/Summer12_EleHZZMoriond2013IDGivenTighterIso/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Summer12_EleHZZMoriond2013IDGivenTighterIso/probes.root\",\"results/Summer12_EleHZZMoriond2013IDGivenTighterIso/basic_76_106_LowPt_ReweightedToMoriond2013\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",2,1,2,3,\"results/Data2012_EleHZZMoriond2013IDGivenTighterIso/probes.root\",\"results/Data2012_EleHZZMoriond2013IDGivenTighterIso/basic2_76_106_LowPt\",\"all\",1,0,\"results/Summer12_EleHZZMoriond2013IDGivenTighterIso/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#

# 
# Compute Efficiency Scale Factors
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2012_EleHZZMoriond2013IDGivenTighterIso/basic2_76_106_LowPt/eff.root","results/Summer12_EleHZZMoriond2013IDGivenTighterIso/basic_76_106_LowPt_ReweightedToMoriond2013/eff.root","results/Data2012_EleHZZMoriond2013IDGivenTighterIso/basic2_76_106_LowPt/","h2_results_electron_selection","EleHZZMoriond2013IDGivenTighterIso_Moriond2013")'
#


#######################################################
#Run1 Legacy Paper 2012 Data
######################################################
##############
#Low Pt bins
##############


#####
# Use BWxCB signal model
#####
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Summer12_EleHZZRun1LegacyPaperIDGivenTighterIso/probes.root\",\"results/Summer12_EleHZZRun1LegacyPaperIDGivenTighterIso/basic_76_106_LowPt_ReweightedToRun1LegacyPaper\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",2,1,1,3,\"results/Data2012_EleHZZRun1LegacyPaperIDGivenTighterIso/probes.root\",\"results/Data2012_EleHZZRun1LegacyPaperIDGivenTighterIso/basic2_76_106_LowPt\",\"all\",1,0,\"results/Summer12_EleHZZRun1LegacyPaperIDGivenTighterIso/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#
# Compute Efficiency Scale Factors
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2012_EleHZZRun1LegacyPaperIDGivenTighterIso/basic2_76_106_LowPt/eff.root","results/Summer12_EleHZZRun1LegacyPaperIDGivenTighterIso/basic_76_106_LowPt_ReweightedToRun1LegacyPaper/eff.root","results/Data2012_EleHZZRun1LegacyPaperIDGivenTighterIso/basic2_76_106_LowPt/","h2_results_electron_selection","EleHZZRun1LegacyPaperIDGivenTighterIso_Run1LegacyPaper")'


#####
# Use MC template signal model
#####
#root -l -b -q plotEff.C+\(\"el0.lowpt.bins\",2,1,2,3,\"results/Data2012_EleHZZRun1LegacyPaperIDGivenTighterIso/probes.root\",\"results/Data2012_EleHZZRun1LegacyPaperIDGivenTighterIso/basic2_76_106_LowPt_MCTemplate\",\"all\",1,0,\"results/Summer12_EleHZZRun1LegacyPaperIDGivenTighterIso/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#
# Compute Efficiency Scale Factors
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2012_EleHZZRun1LegacyPaperIDGivenTighterIso/basic2_76_106_LowPt_MCTemplate/eff.root","results/Summer12_EleHZZRun1LegacyPaperIDGivenTighterIso/basic_76_106_LowPt_ReweightedToRun1LegacyPaper/eff.root","results/Data2012_EleHZZRun1LegacyPaperIDGivenTighterIso/basic2_76_106_LowPt_MCTemplate/","h2_results_electron_selection","EleHZZRun1LegacyPaperIDGivenTighterIso_Run1LegacyPaper")'





#--------------------------------------------------------------
#  Electron HZZ Moriond 2013 With ZeeGamma
#==============================================================

####################################
# 
# Measure Efficiency in MC and Data
#
#
##########
#2011
#########
#
#root -l -b -q plotEff_ZeeGamma.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Fall11_EleHZZICHEP2012WPWithZeeGamma/probes.root\",\"results/Fall11_EleHZZICHEP2012WPWithZeeGamma/basic_76_106_LowPt_ReweightedToFull2011\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2011/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff_ZeeGamma.C+\(\"el0.lowpt.bins\",2,0,2,7,\"results/Data2011_EleHZZICHEP2012WPWithZeeGamma/probes.root\",\"results/Data2011_EleHZZICHEP2012WPWithZeeGamma/basic2_76_106_LowPt\",\"all\",1,0,\"results/Fall11_EleHZZICHEP2012WPWithZeeGamma/probes.root\",\"results/Fall11_EleHZZICHEP2012WPWithZeeGammaBkg/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2011/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)
#
##########
#2012
#########
#
#root -l -b -q plotEff_ZeeGamma.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Summer12_EleHZZICHEP2012WPWithZeeGamma/probes.root\",\"results/Summer12_EleHZZICHEP2012WPWithZeeGamma/basic_76_106_ReweightedTo2012C_LowPt\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY_To_Run2012C.root\"\)
#root -l -b -q plotEff_ZeeGamma.C+\(\"el0.lowpt.bins\",2,0,2,7,\"results/Data2012_EleHZZICHEP2012WPWithZeeGamma/probes.root\",\"results/Data2012_EleHZZICHEP2012WPWithZeeGamma/basic2_76_106_LowPt\",\"all\",1,0,\"results/Summer12_EleHZZICHEP2012WPWithZeeGamma/probes.root\",\"results/Summer12_EleHZZICHEP2012WPWithZeeGammaBkg/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY_To_Run2012.root\"\)
#
#
##########
#2012C
#########
#
#root -l -b -q plotEff_ZeeGamma.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Summer12_EleHZZICHEP2012WPWithZeeGamma/probes.root\",\"results/Summer12_EleHZZICHEP2012WPWithZeeGamma/basic_76_106_ReweightedTo2012C_LowPt\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY_To_Run2012C.root\"\)
#root -l -b -q plotEff_ZeeGamma.C+\(\"el0.lowpt.bins\",2,0,2,7,\"results/Data2012_EleHZZICHEP2012WPWithZeeGamma/probes.root\",\"results/Data2012_EleHZZICHEP2012WPWithZeeGamma/basic2_76_106_LowPt\",\"all\",1,0,\"results/Summer12_EleHZZICHEP2012WPWithZeeGamma/probes.root\",\"results/Summer12_EleHZZICHEP2012WPWithZeeGammaBkg/probes.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY_To_Run2012C.root\"\)
#
#
##########
#2013 Moriond
#########
#
#root -l -b -q plotEff_ZeeGamma.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Summer12_EleHZZMoriond2013WPWithZeeGamma/probes.root\",\"results/Summer12_EleHZZMoriond2013WPWithZeeGamma/basic_76_106_LowPt_ReweightedToMoriond2013\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#root -l -b -q plotEff_ZeeGamma.C+\(\"el0.lowpt.bins\",2,0,2,7,\"results/Data2012_EleHZZMoriond2013WPWithZeeGamma/probes.root\",\"results/Data2012_EleHZZMoriond2013WPWithZeeGamma/basic2_76_106_LowPt\",\"all\",1,0,\"results/Summer12_EleHZZMoriond2013WPWithZeeGamma/probes.root\",\"results/Summer12_EleHZZMoriond2013WPWithZeeGamma/probes.bkg.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#

####################################

####################################
# 
# Compute Efficiency Scale Factors
#
#######
2011
#######
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2011_EleHZZICHEP2012WPWithZeeGamma/basic2_76_106_LowPt/eff.root","results/Fall11_EleHZZICHEP2012WPWithZeeGamma/basic_76_106_LowPt_ReweightedToFull2011/eff.root","results/Data2011_EleHZZICHEP2012WPWithZeeGamma/basic2_76_106_LowPt/","h2_results_electron_selection","EleHZZICHEP2012WPWithZeeGamma_Full2011")'
#
#
#######
2012
#######
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2012_EleHZZMoriond2013WPWithZeeGamma/basic2_76_106_LowPt/eff.root","results/Summer12_EleHZZMoriond2013WPWithZeeGamma/basic_76_106_LowPt_ReweightedToMoriond2013/eff.root","results/Data2012_EleHZZMoriond2013WPWithZeeGamma/basic2_76_106_LowPt/","h2_results_electron_selection","EleHZZMoriond2013WPWithZeeGamma_Moriond2013")'
####################################




#--------------------------------------------------------------
#  Electron HZZ Moriond2013 Iso Only With ZeeGamma 
#==============================================================

####################################
# 
# Measure Efficiency in MC and Data
#
#
##########
#2012
#########
#
#root -l -b -q plotEff_ZeeGamma.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Summer12_EleHZZMoriond2013IsoWithZeeGamma/probes.root\",\"results/Summer12_EleHZZMoriond2013IsoWithZeeGamma/basic_76_106_ReweightedToMoriond2013_LowPt\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#root -l -b -q plotEff_ZeeGamma.C+\(\"el0.lowpt.bins\",1,0,1,7,\"results/Data2012_EleHZZMoriond2013IsoWithZeeGamma/probes.root\",\"results/Data2012_EleHZZMoriond2013IsoWithZeeGamma/basic2_76_106_LowPt\",\"all\",1,0,\"results/Summer12_EleHZZMoriond2013IsoWithZeeGamma/probes.root\",\"results/Summer12_EleHZZMoriond2013IsoWithZeeGamma/probes.bkg.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#
####################################

####################################
# 
# Compute Efficiency Scale Factors
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2011_EleHZZMoriond2013IsoWithZeeGamma/basic2_76_106/eff.root","results/Fall11_EleHZZMoriond2013IsoWithZeeGamma/basic_76_106_ReweightedToFull2011/eff.root","results/Data2011_EleHZZMoriond2013IsoWithZeeGamma/basic2_76_106/","h2_results_electron_selection","EleHZZMoriond2013IsoWithZeeGamma_Full2011")'
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2012_EleHZZMoriond2013IsoWithZeeGamma/basic2_76_106_LowPt/eff.root","results/Summer12_EleHZZMoriond2013IsoWithZeeGamma/basic_76_106_ReweightedToMoriond2013_LowPt/eff.root","results/Data2012_EleHZZMoriond2013IsoWithZeeGamma/basic2_76_106_LowPt/","h2_results_electron_selection","EleHZZMoriond2013IsoWithZeeGamma_Full2012")'
####################################




#--------------------------------------------------------------
#  Electron HZZ Moriond2013 ID Given Iso With ZeeGamma
#==============================================================

####################################
# 
# Measure Efficiency in MC and Data
#
#
##########
#2012
#########
#
#root -l -b -q plotEff_ZeeGamma.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Summer12_EleHZZMoriond2013IDGivenIsoWithZeeGamma/probes.root\",\"results/Summer12_EleHZZMoriond2013IDGivenIsoWithZeeGamma/basic_76_106_ReweightedToMoriond2013_LowPt\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
###############
# Do Bkg Fit
###############
#root -l -b -q plotEff_ZeeGamma.C+\(\"el0.lowpt.bins\",1,0,2,7,\"results/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/probes.root\",\"results/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/basic2_76_106_LowPt\",\"all\",1,0,\"results/Summer12_EleHZZMoriond2013IDGivenIsoWithZeeGamma/probes.root\",\"results/Summer12_EleHZZMoriond2013IDGivenIsoWithZeeGamma/probes.bkg.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#
#
###############
# No Bkg Assumption
###############
#root -l -b -q plotEff_ZeeGamma.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/probes.root\",\"results/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/basic2_76_106_LowPt\",\"all\",1,0,\"results/Summer12_EleHZZMoriond2013IDGivenIsoWithZeeGamma/probes.root\",\"results/Summer12_EleHZZMoriond2013IDGivenIsoWithZeeGamma/probes.bkg.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#
####################################

####################################
# 
# Compute Efficiency Scale Factors
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2011_EleHZZMoriond2013IDGivenIsoWithZeeGamma/basic2_76_106/eff.root","results/Fall11_EleHZZMoriond2013IDGivenIsoWithZeeGamma/basic_76_106_ReweightedToFull2011/eff.root","results/Data2011_EleHZZMoriond2013IDGivenIsoWithZeeGamma/basic2_76_106/","h2_results_electron_selection","EleHZZMoriond2013IDGivenIsoWithZeeGamma_Full2011")'
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/basic2_76_106_LowPt/eff.root","results/Summer12_EleHZZMoriond2013IDGivenIsoWithZeeGamma/basic_76_106_ReweightedToMoriond2013_LowPt/eff.root","results/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/basic2_76_106_LowPt/","h2_results_electron_selection","EleHZZMoriond2013IDGivenIsoWithZeeGamma_Full2012")'
####################################




#--------------------------------------------------------------
#  Electron HZZ Moriond2013 Iso Given ID With ZeeGamma
#==============================================================

####################################
# 
# Measure Efficiency in MC and Data
#
#
##########
#2012
#########
#
#root -l -b -q plotEff_ZeeGamma.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Summer12_EleHZZMoriond2013IsoGivenIDWithZeeGamma/probes.root\",\"results/Summer12_EleHZZMoriond2013IsoGivenIDWithZeeGamma/basic_76_106_ReweightedToMoriond2013_LowPt\",\"all\",1,0,\"\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#
###############
# Do Bkg Fit
###############
#root -l -b -q plotEff_ZeeGamma.C+\(\"el0.lowpt.bins\",2,0,2,7,\"results/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/probes.root\",\"results/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/basic2_76_106_LowPt\",\"all\",1,0,\"results/Summer12_EleHZZMoriond2013IsoGivenIDWithZeeGamma/probes.root\",\"results/Summer12_EleHZZMoriond2013IsoGivenIDWithZeeGamma/probes.bkg.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#
###############
# No Bkg Assumption
###############
#root -l -b -q plotEff_ZeeGamma.C+\(\"el0.lowpt.bins\",0,0,0,0,\"results/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/probes.root\",\"results/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/basic2_76_106_LowPt\",\"all\",1,0,\"results/Summer12_EleHZZMoriond2013IsoGivenIDWithZeeGamma/probes.root\",\"results/Summer12_EleHZZMoriond2013IsoGivenIDWithZeeGamma/probes.bkg.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\"\)
#
####################################

####################################
# 
# Compute Efficiency Scale Factors
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2011_EleHZZMoriond2013IsoGivenIDWithZeeGamma/basic2_76_106/eff.root","results/Fall11_EleHZZMoriond2013IsoGivenIDWithZeeGamma/basic_76_106_ReweightedToFull2011/eff.root","results/Data2011_EleHZZMoriond2013IsoGivenIDWithZeeGamma/basic2_76_106/","h2_results_electron_selection","EleHZZMoriond2013IsoGivenIDWithZeeGamma_Full2011")'
#
# root -l -q makeEfficiencyScaleFactors.C+'("results/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/basic2_76_106_LowPt/eff.root","results/Summer12_EleHZZMoriond2013IsoGivenIDWithZeeGamma/basic_76_106_ReweightedToMoriond2013_LowPt/eff.root","results/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/basic2_76_106_LowPt/","h2_results_electron_selection","EleHZZMoriond2013IsoGivenIDWithZeeGamma_Full2012")'
####################################





#---------------------------------------------------------------------------------------------------
#  Compute Scale Factor Combining Zee, Zee Low Pt Bins, and ZeeGamma Low Pt Bins
#==============================================================-------------------------------------


#########
#2011
#########
# root -l -q makeEfficiencyScaleFactors.C+'("StandardEffBins.root","results/Data2011_EleHZZICHEP2012WPWithTightTag/basic2_76_106/eff.root","results/Data2011_EleHZZICHEP2012WPWithTightTag/basic2_76_106_LowPt/eff.root","results/Data2011_EleHZZICHEP2012WPWithZeeGamma/basic2_76_106_LowPt/eff.root","results/Fall11_EleHZZICHEP2012WPWithTightTag/basic_76_106_ReweightedToRun2011/eff.root","results/Fall11_EleHZZICHEP2012WPWithTightTag/basic_76_106_LowPt_ReweightedToRun2011/eff.root","results/Fall11_EleHZZICHEP2012WPWithZeeGamma/basic_76_106_LowPt_ReweightedToFull2011/eff.root","results/Data2011_EleHZZICHEP2012WPWithZeeGamma/basic2_76_106_LowPt/","heff_electron_selection","EleHZZICHEP2012WPMixed_Full2011")'

##################
#2013 Moriond
##################
# root -l -q makeEfficiencyScaleFactors.C+'("StandardEffBins.root","results/Data2012_EleHZZHCP2012WPWithTightTag/basic2_76_106/eff.root","results/Data2012_EleHZZHCP2012WPWithTightTag/basic2_76_106_LowPt/eff.root","results/Data2012_EleHZZHCP2012WPWithZeeGamma/basic2_76_106_LowPt/eff.root","results/Summer12_EleHZZHCP2012WPWithTightTag/basic_76_106_ReweightedToHCP2012/eff.root","results/Summer12_EleHZZHCP2012WPWithTightTag/basic_76_106_LowPt_ReweightedToHCP2012/eff.root","results/Summer12_EleHZZHCP2012WPWithZeeGamma/basic_76_106_LowPt_ReweightedToHCP2012/eff.root","results/Data2012_EleHZZHCP2012WPWithCombinedMethods/","heff_electron_selection","EleHZZHCP2012WPMixed_HCP2012")'

#---------------------------------------------------------------------------------------------------
#  Compute Scale Factor Combining Zee, Zee Low Pt Bins, and Adhoc Factorized Zee for 7-10 bins
#==============================================================-------------------------------------

########################################################################
#Run1 Legacy Paper
########################################################################
###############
#2011 Data
###############
# mkdir -p results/Data2011_EleHZZRun1LegacyPaperWPWithCombinedMethods/
# root -l -q makeEfficiencyScaleFactors.C+'("StandardEffBins.root","results/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_WithTurnOn/eff.root","results/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_LowPt/eff.root","results/Fall11_EleHZZRun1LegacyPaperWPWithTightTag/basic_76_106_ReweightedToRun2011/eff.root","results/Fall11_EleHZZRun1LegacyPaperWPWithTightTag/basic_76_106_LowPt_ReweightedToRun2011/eff.root",2011,"results/Data2011_EleHZZRun1LegacyPaperWPWithCombinedMethods/","heff_electron_selection","EleHZZRun1LegacyPaperWPFromZeeWithLowestPtBinFactorized_Run1LegacyPaper_2011Data")'


###############
#2012 Data
###############
# mkdir -p results/Data2012_EleHZZRun1LegacyPaperWPWithCombinedMethods/
# root -l -q makeEfficiencyScaleFactors.C+'("StandardEffBins.root","results/Data2012_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_WithTurnOn/eff.root","results/Data2012_EleHZZRun1LegacyPaperWPWithTightTag/basic2_76_106_LowPt/eff.root","results/Summer12_EleHZZRun1LegacyPaperWPWithTightTag/basic_76_106_ReweightedToRun1LegacyPaper/eff.root","results/Summer12_EleHZZRun1LegacyPaperWPWithTightTag/basic_76_106_LowPt_ReweightedToRun1LegacyPaper/eff.root",2012,"results/Data2012_EleHZZRun1LegacyPaperWPWithCombinedMethods/","heff_electron_selection","EleHZZRun1LegacyPaperWPFromZeeWithLowestPtBinFactorized_Run1LegacyPaper_2012Data")'


#---------------------------------------------------------------------------------------------------
#  Plot efficiency comparisons
#==============================================================-------------------------------------

#root -l plotEfficiency.C+

