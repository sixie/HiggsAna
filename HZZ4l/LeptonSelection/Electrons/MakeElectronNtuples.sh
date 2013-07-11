#######################################################################################
# QCDFake Sample
#######################################################################################
root -l -b -q HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeFakeElectronTrainingNtuple.C+


#######################################################################################
# Z+Fake Sample
#######################################################################################
root -l -b -q HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeFakeElectronTrainingNtupleFromZPlusJetSample.C+


#######################################################################################
# W+Fake Sample
#######################################################################################



######################
# Run2012A SingleMu
######################
foreach f ( `ls /data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu*Tight* | sed 's/\/data\/smurf\/sixie\/hist\/HZZ4lNtuples\/data\/unprocessed\///' ` )
root -l -b -q HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeFakeElectronTrainingNtupleFromWPlusJetSample.C+\(\"/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/$f\",\"/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToMuNuPlusJet.${f}.root\"\)
end

######################
# Run2012B SingleMu
######################
foreach f ( `ls /data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu*Tight* | sed 's/\/data\/smurf\/sixie\/hist\/HZZ4lNtuples\/data\/unprocessed\///' ` )
root -l -b -q HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeFakeElectronTrainingNtupleFromWPlusJetSample.C+\(\"/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/$f\",\"/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToMuNuPlusJet.${f}.root\"\)
end


######################
# Run2012A SingleEle
######################
foreach f ( `ls /data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel*Tight* | sed 's/\/data\/smurf\/sixie\/hist\/HZZ4lNtuples\/data\/unprocessed\///' ` )
root -l -b -q HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeFakeElectronTrainingNtupleFromWPlusJetSample.C+\(\"/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/$f\",\"/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToENuPlusJet.${f}.root\"\)
end

#Something got screwed up here. I'm getting signal electrons in this sample...
#don't use this for now.

######################
# Run2012B SingleEle
######################
foreach f ( `ls /data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel*Tight* | sed 's/\/data\/smurf\/sixie\/hist\/HZZ4lNtuples\/data\/unprocessed\///' ` )
root -l -b -q HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeFakeElectronTrainingNtupleFromWPlusJetSample.C+\(\"/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/$f\",\"/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToENuPlusJet.${f}.root\"\)
end


# Merge
hadd -f /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToENuPlusJet_2012A.root  /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToENuPlusJet.AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_*Tight*.root
hadd -f /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToENuPlusJet_2012B.root  /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToENuPlusJet.AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_*Tight*.root
hadd -f /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToMuNuPlusJet_2012A.root  /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToMuNuPlusJet.AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_*Tight*.root
hadd -f /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToMuNuPlusJet_2012B.root  /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToMuNuPlusJet.AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_*Tight*.root
hadd -f /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WPlusJet_2012.root /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToENuPlusJet_2012A.root /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToENuPlusJet_2012B.root /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToMuNuPlusJet_2012A.root /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToMuNuPlusJet_2012B.root

# Cleanup
mv /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToMuNuPlusJet.* /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToENuPlusJet.* /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/done/


#######################################################################################
# MC Fake Samples
#######################################################################################

#Z+jets
#ttbar
#W+jets




#######################################################################################
# Data Z Tag and Probe Sample
#######################################################################################

#Run2012A
foreach f ( `ls /data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-del*Tight* | sed 's/\/data\/smurf\/sixie\/hist\/HZZ4lNtuples\/data\/unprocessed\///' ` )
root -l -b -q HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeRealElectronTrainingNtuple.C+\(\"/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/$f\",\"/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_ZeeTagAndProbe.${f}.root\"\)
end

#Run2012B
foreach f ( `ls /data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-del*Tight* | sed 's/\/data\/smurf\/sixie\/hist\/HZZ4lNtuples\/data\/unprocessed\///' ` )
root -l -b -q HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeRealElectronTrainingNtuple.C+\(\"/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/$f\",\"/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_ZeeTagAndProbe.${f}.root\"\)
end

# Merge
hadd -f /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_ZeeTagAndProbe.root /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_ZeeTagAndProbe.AllNtuple_HZZ4lNtuple_*root
# Cleanup
mv /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_ZeeTagAndProbe.AllNtuple_HZZ4lNtuple_*root /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/done/



#######################################################################################
# MC Real Samples
#######################################################################################

#DY->ee

#HZZ
foreach f ( `cat fileListB | grep HZZ | awk '{print $5}' | sed 's/\/store\/group\/phys_higgs\/cmshzz4l\/sixie\/hist\/AllNtuple\/cern\/filefi\/028\///' | grep "h...zz4l" ` )
root -l -b -q HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeMCElectronTrainingNtuple.C+\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshzz4l/sixie/hist/AllNtuple/cern/filefi/028/$f\",\"/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_MC.${f}.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root\",kTRUE\)
end

#Merge them together : don't use 125 sample
hadd -f /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Summer12HZZ.root /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_MC.AllNtuple_HZZ4lNtuple_s12-h120zz4l* /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_MC.AllNtuple_HZZ4lNtuple_s12-h120zz4l* /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_MC.AllNtuple_HZZ4lNtuple_s12-h121zz4l* /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_MC.AllNtuple_HZZ4lNtuple_s12-h123zz4l* /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_MC.AllNtuple_HZZ4lNtuple_s12-h124zz4l* /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_MC.AllNtuple_HZZ4lNtuple_s12-h126zz4l* /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_MC.AllNtuple_HZZ4lNtuple_s12-h127zz4l* /data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_MC.AllNtuple_HZZ4lNtuple_s12-h130zz4l* 



#HWW
foreach f ( `cat fileListHiggsHWW2| grep "HZZ4lNtuple" | grep ww `)
root -l -b -q HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeMCElectronTrainingNtuple.C+\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshww/sixie/hist/AllNtuple/cern/filefi/028/$f\",\"/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Real_MC.${f}.root\",\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root\",kTRUE\)
end

#Merge them together: don't use 125 sample
hadd -f /data/blue/sixie/LeptonSelection/Electrons/ElectronSelectionTraining.Summer12HWW.root /data/blue/sixie/LeptonSelection/Electrons/ElectronSelectionTraining.Real_MC.AllNtuple_HZZ4lNtuple_s12-h115ww* /data/blue/sixie/LeptonSelection/Electrons/ElectronSelectionTraining.Real_MC.AllNtuple_HZZ4lNtuple_s12-h120ww* /data/blue/sixie/LeptonSelection/Electrons/ElectronSelectionTraining.Real_MC.AllNtuple_HZZ4lNtuple_s12-h130ww*
hadd -f /data/blue/sixie/LeptonSelection/Electrons/ElectronSelectionTraining.Summer12HWW125.root /data/blue/sixie/LeptonSelection/Electrons/ElectronSelectionTraining.Real_MC.AllNtuple_HZZ4lNtuple_s12-h125ww* 
