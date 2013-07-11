#!/bin/tcsh
#===================================================================================================
# Submit a set of jobs to run over a given dataset, splitting the jobs according to the filesets.
#
# Version 1.0                                                                      November 14, 2008
#===================================================================================================




####################################################
#
# HZZ Run1LegacyPaper Selection + Run1LegacyPaper Data
#
####################################################


#############
#42X MC
#############

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep f11-zeem20-powheg `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaper_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaper_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWPPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",1,-1,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Fall11_EleHZZRun1LegacyPaper/
end

#############
#53X MC
#############

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep s12-zllm50-v7a `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaper_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaper_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWPPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",1,-1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZRun1LegacyPaper/
end


#############
#Data
#############
##Run2011##
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r11`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaper_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaper_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWPPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",0,0,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2011_EleHZZRun1LegacyPaper/
    sleep 1
end

##Run2012A##
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12a-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaper_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaper_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWPPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaper/
    sleep 1
end

##Run2012B##
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12b-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaper_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaper_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWPPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaper/
    sleep 1
end

##Run2012C##
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12c-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaper_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaper_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWPPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaper/
    sleep 1
end

##Run2012D##
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12d-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaper_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaper_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWPPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaper/
    sleep 1
end



#############
#
# HZZ Run1LegacyPaper ID Given Tighter Iso
#
#############


#############
#MC
#############
#42X MC
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep f11-zeem20-powheg`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIDGivenTighterIso_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIDGivenTighterIso_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIDGivenTighterIsoPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",1,-1,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Fall11_EleHZZRun1LegacyPaperIDGivenTighterIso/
  sleep 1
end

#53X MC
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep s12-zllm50-v7a`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIDGivenTighterIso_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIDGivenTighterIso_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIDGivenTighterIsoPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",1,-1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZRun1LegacyPaperIDGivenTighterIso/
  sleep 1
end



#############
#Data
#############
##Run2011A#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r11`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIDGivenTighterIso_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIDGivenTighterIso_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIDGivenTighterIsoPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",0,0,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2011_EleHZZRun1LegacyPaperIDGivenTighterIso/
  sleep 1
end

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep sel | grep HZZ | grep r11`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIDGivenTighterIso_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIDGivenTighterIso_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIDGivenTighterIsoPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",0,1,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2011_EleHZZRun1LegacyPaperIDGivenTighterIso/
  sleep 1
end


##Run2012A#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12a-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIDGivenTighterIso_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIDGivenTighterIso_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIDGivenTighterIsoPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaperIDGivenTighterIso/
  sleep 1
end

##Run2012B#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12b-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIDGivenTighterIso_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIDGivenTighterIso_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIDGivenTighterIsoPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaperIDGivenTighterIso/
  sleep 1
end

##Run2012C#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12c-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIDGivenTighterIso_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIDGivenTighterIso_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIDGivenTighterIsoPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaperIDGivenTighterIso/
  sleep 1
end

##Run2012D#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12d-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIDGivenTighterIso_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIDGivenTighterIso_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIDGivenTighterIsoPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaperIDGivenTighterIso/
  sleep 1
end






#############
#
# HZZ Run1LegacyPaper WP with Tight Tag Cuts & MET Cut
#
#############

#############
#MC
#############
#42X MC
foreach file(` cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ | grep f11-zeem20-powheg `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperWPWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperWPWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWPWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",1,-1,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Fall11_EleHZZRun1LegacyPaperWPWithTightTag/
  sleep 1
end


#53X MC
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep s12-zllm50-v7a`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperWPWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperWPWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWPWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",1,-1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZRun1LegacyPaperWPWithTightTag/
  sleep 1
end



#############
#Data
#############
##Run2011#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r11 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperWPWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperWPWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWPWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",0,0,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/
  sleep 1
end

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep sel | grep HZZ | grep r11 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperWPWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperWPWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWPWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",0,1,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2011_EleHZZRun1LegacyPaperWPWithTightTag/
  sleep 1
end


##Run2012A#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12a-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperWPWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperWPWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWPWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaperWPWithTightTag/
  sleep 1
end

##Run2012B#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12b-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperWPWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperWPWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWPWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaperWPWithTightTag/
  sleep 1
end

##Run2012C#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12c-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperWPWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperWPWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWPWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaperWPWithTightTag/
  sleep 1
end

##Run2012D#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12d-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperWPWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperWPWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWPWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaperWPWithTightTag/
  sleep 1
end




#############
#
# HZZ Run1LegacyPaper ID Given Iso WP with Tight Tag Cuts & MET Cut
#
#############


#############
#MC
#############

#42X MC
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep f11-zeem20-powheg`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",1,-1,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Fall11_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/
  sleep 1
end

#53X MC
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep s12-zllm50-v7a`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",1,-1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/
  sleep 1
end



#############
#Data
#############
##Run2011#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r11`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",0,0,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2011_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep sel | grep HZZ | grep r11`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",0,1,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2011_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/
  sleep 1
end



##Run2012A#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12a-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/
  sleep 1
end

##Run2012B#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12b-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/
  sleep 1
end

##Run2012C#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12c-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/
  sleep 1
end

##Run2012D#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12d-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIDGivenIsoWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaperIDGivenIsoWithTightTag/
  sleep 1
end





#############
#
# HZZ Run1LegacyPaper Iso Given ID WP with Tight Tag Cuts & MET Cut
#
#############


#############
#MC
#############
#42X MC
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep f11-zeem20-powheg`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",1,-1,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Fall11_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/
  sleep 1
end

#53X MC
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep s12-zllm50-v7a`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",1,-1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/
  sleep 1
end



#############
#Data
#############
##Run2011#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r11`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",0,0,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2011_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep sel | grep HZZ | grep r11`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",0,1,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2011_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/
  sleep 1
end


##Run2012A#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12a-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/
  sleep 1
end

##Run2012B#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12b-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/
  sleep 1
end

##Run2012C#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12c-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/
  sleep 1
end

##Run2012D#
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep del | grep HZZ | grep r12d-del-j22-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTag_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTag_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperIsoGivenIDWithTightTagPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZRun1LegacyPaperIsoGivenIDWithTightTag/
  sleep 1
end




#############
#
# HZZ Run1LegacyPaper Z->4l Electron Tag and Probe
#
#############

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListB2 | grep HZZ | grep powheg | grep 4e `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaper_Zto4l_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaper_Zto4l_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWP_Zto4lPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshzz4l/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"probes.${file}\",1,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZRun1LegacyPaperWP_Zto4l/
end


foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListB2 | grep HZZ | grep powheg | grep 2e2m `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaper_Zto4l_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaper_Zto4l_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWP_Zto4lPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshzz4l/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"probes.${file}\",1,3,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZRun1LegacyPaperWP_Zto4l/
end


foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListStoreUser2 | grep HZZ | grep zll`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZRun1LegacyPaper_Zto4l_${file}.out -J Efficiency_selectEleHZZRun1LegacyPaper_Zto4l_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZRun1LegacyPaperWP_Zto4lPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshzz4l/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZRun1LegacyPaperWP_Zto4l/
end



