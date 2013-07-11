#!/bin/tcsh
#===================================================================================================
# Submit a set of jobs to run over a given dataset, splitting the jobs according to the filesets.
#
# Version 1.0                                                                      November 14, 2008
#===================================================================================================




#############
#
# HZZ HCP2012 WP with Z->ee+gamma
#
#############



#############
#MC - 42X
#############

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ4lNtuple | grep f11-zeem20-powheg `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",1,-1,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Fall11_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ4lNtuple | grep f11-zeem20-powheg `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGammaBkg_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGammaBkg_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",1,-2,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Fall11_EleHZZMoriond2013WPWithZeeGammaBkg/
  sleep 1
end



#############
#MC - 52X
#############

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListHiggsHZZ4l2 | grep HZZ4lNtuple | grep s12-zllm50-2-v9 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshzz4l/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"probes.${file}\",1,-1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListHiggsHZZ4l2 | grep HZZ4lNtuple | grep s12-zllm50-2-v9 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGammaBkg_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGammaBkg_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshzz4l/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"probes.${file}\",1,-2,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013WPWithZeeGammaBkg/
  sleep 1
end



#############
#MC - 53X
#############

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep s12-zllm50-v7a `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",1,-1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep s12-zllm50-v7a `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGammaBkg_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGammaBkg_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",1,-2,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013WPWithZeeGammaBkg/
  sleep 1
end



#############
#Data - Run2011 
#############
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep del | grep HZZ | grep r11 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",0,0,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2011_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end


foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep sel | grep HZZ | grep r11`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"probes.${file}\",0,1,1\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2011_SingleEle_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end



#############
#Data - RUn2012 A,B
#############
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ | grep r12a-del-j13-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ | grep r12a-del-a06-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ | grep r12b-del-j13-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end


foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12a-sel-j13-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12a-sel-a06-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12b-sel-j13-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end



#############
#Data - RUn2012 C
#############
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12c-del-a24-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12c-del-pr-v2 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end


foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12c-sel-a24-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12c-sel-pr-v2`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end


#############
#Data - RUn2012 D
#############

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12d-del-pr-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12d-sel-pr-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013WPWithZeeGamma/
  sleep 1
end




#############
#
# HZZ Moriond2013 IDGivenIso with Z->ee+gamma
#
#############



#############
#MC - 53X
#############

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep s12-zllm50-v7a `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IDGivenIsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",1,-1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IDGivenIsoWithZeeGamma/
  sleep 1
end

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep s12-zllm50-v7a `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IDGivenIsoWithZeeGammaBkg_${file}.out -J Efficiency_selectEleHZZMoriond2013IDGivenIsoWithZeeGammaBkg_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IDGivenIsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",1,-2,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IDGivenIsoWithZeeGammaBkg/
  sleep 1
end



#############
#Data - RUn2012 A,B
#############
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ | grep r12a-del-j13-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IDGivenIsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ | grep r12a-del-a06-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IDGivenIsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ | grep r12b-del-j13-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IDGivenIsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/
  sleep 1
end


foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12a-sel-j13-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IDGivenIsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IDGivenIsoWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12a-sel-a06-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IDGivenIsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IDGivenIsoWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12b-sel-j13-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IDGivenIsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IDGivenIsoWithZeeGamma/
  sleep 1
end


#############
#Data - RUn2012 C
#############
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12c-del-a24-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IDGivenIsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12c-del-pr-v2 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IDGivenIsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/
  sleep 1
end


foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12c-sel-a24-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IDGivenIsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IDGivenIsoWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12c-sel-pr-v2`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IDGivenIsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IDGivenIsoWithZeeGamma/
  sleep 1
end

#############
#Data - RUn2012 D
#############
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12d-del-pr-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IDGivenIsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/
  sleep 1
end

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12d-sel-pr-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IDGivenIsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IDGivenIsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IDGivenIsoWithZeeGamma/
  sleep 1
end





#Merge
hadd -f /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/probes.del.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/*.root
hadd -f /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/probes.sel.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IDGivenIsoWithZeeGamma/*.root
hadd -f /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/probes.root  /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/probes.del.root  /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IDGivenIsoWithZeeGamma/probes.sel.root
hadd -f /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IDGivenIsoWithZeeGamma/probes.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IDGivenIsoWithZeeGamma/*.root
hadd -f /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IDGivenIsoWithZeeGammaBkg/probes.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IDGivenIsoWithZeeGammaBkg/*.root



#############
#
# HZZ Moriond2013 IsoGivenID with Z->ee+gamma
#
#############


#############
#MC - 53X
#############

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep s12-zllm50-v7a `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoGivenIDWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",1,-1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IsoGivenIDWithZeeGamma/
  sleep 1
end

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep s12-zllm50-v7a `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoGivenIDWithZeeGammaBkg_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoGivenIDWithZeeGammaBkg_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoGivenIDWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",1,-2,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IsoGivenIDWithZeeGammaBkg/
  sleep 1
end



#############
#Data - RUn2012 A,B
#############
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ | grep r12a-del-j13-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoGivenIDWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ | grep r12a-del-a06-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoGivenIDWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ | grep r12b-del-j13-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoGivenIDWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/
  sleep 1
end


foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12a-sel-j13-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoGivenIDWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IsoGivenIDWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12a-sel-a06-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoGivenIDWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IsoGivenIDWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12b-sel-j13-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoGivenIDWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IsoGivenIDWithZeeGamma/
  sleep 1
end


#############
#Data - RUn2012 C
#############
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12c-del-a24-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoGivenIDWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12c-del-pr-v2 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoGivenIDWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/
  sleep 1
end


foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12c-sel-a24-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoGivenIDWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IsoGivenIDWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12c-sel-pr-v2`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoGivenIDWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IsoGivenIDWithZeeGamma/
  sleep 1
end

#############
#Data - RUn2012 D
#############
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12d-del-pr-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoGivenIDWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/
  sleep 1
end

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12d-sel-pr-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoGivenIDWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoGivenIDWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IsoGivenIDWithZeeGamma/
  sleep 1
end





#Merge
hadd -f /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/probes.del.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/*.root
hadd -f /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/probes.sel.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IsoGivenIDWithZeeGamma/*.root
hadd -f /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/probes.root  /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/probes.del.root  /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoGivenIDWithZeeGamma/probes.sel.root
hadd -f /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IsoGivenIDWithZeeGamma/probes.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IsoGivenIDWithZeeGamma/*.root
hadd -f /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IsoGivenIDWithZeeGammaBkg/probes.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IsoGivenIDWithZeeGammaBkg/*.root



#############
#
# HZZ Moriond2013 Iso with Z->ee+gamma
#
#############



#############
#MC - 53X
#############

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep s12-zllm50-v7a `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",1,-1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IsoWithZeeGamma/
  sleep 1
end

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple | grep s12-zllm50-v7a `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoWithZeeGammaBkg_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoWithZeeGammaBkg_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",1,-2,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IsoWithZeeGammaBkg/
  sleep 1
end



#############
#Data - RUn2012 A,B
#############
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ | grep r12a-del-j13-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ | grep r12a-del-a06-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ | grep r12b-del-j13-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoWithZeeGamma/
  sleep 1
end


foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12a-sel-j13-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IsoWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12a-sel-a06-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IsoWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12b-sel-j13-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IsoWithZeeGamma/
  sleep 1
end


#############
#Data - RUn2012 C
#############
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12c-del-a24-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12c-del-pr-v2 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoWithZeeGamma/
  sleep 1
end


foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12c-sel-a24-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IsoWithZeeGamma/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12c-sel-pr-v2`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IsoWithZeeGamma/
  sleep 1
end


#############
#Data - RUn2012 D
#############
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12d-del-pr-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoWithZeeGamma/
  sleep 1
end

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep HZZ | grep r12d-sel-pr-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013IsoWithZeeGamma_${file}.out -J Efficiency_selectEleHZZMoriond2013IsoWithZeeGamma_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012IsoWithZeeGammaPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IsoWithZeeGamma/
  sleep 1
end



#Merge
hadd -f /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoWithZeeGamma/probes.del.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoWithZeeGamma/*.root
hadd -f /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoWithZeeGamma/probes.sel.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013IsoWithZeeGamma/*.root
hadd -f /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoWithZeeGamma/probes.root  /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoWithZeeGamma/probes.del.root  /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013IsoWithZeeGamma/probes.sel.root
hadd -f /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IsoWithZeeGamma/probes.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IsoWithZeeGamma/*.root
hadd -f /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IsoWithZeeGammaBkg/probes.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013IsoWithZeeGammaBkg/*.root







#############
#
# HZZ Moriond2013 WP with Z->ee+gamma - Selection Studies
#
#############




#############
#MC
#############

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/s12-zee.conf | grep AllNtuple | awk '{print $1}' | sed 's/root:\/\/eoscms\/\/eos\/cms\/store\/group\/phys_higgs\/cmshzz4l\/sixie\/hist\/AllNtuple\/cern\/filefi\/028\///' `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGammaStudy_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGammaStudy_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaStudyPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshzz4l/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"probes.${file}\",1,-1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZMoriond2013WPWithZeeGammaStudy/
  sleep 1
end


#############
#Data
#############
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/data_el.2012.EOS.conf | grep AllNtuple | awk '{print $1}' | sed 's/root:\/\/eoscms\/\/eos\/cms\/store\/group\/phys_egamma\/sixie\/hist\/AllNtuple\/cern\/filefi\/028\///' `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGammaStudy_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGammaStudy_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaStudyPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"probes.${file}\",0,0,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_EleHZZMoriond2013WPWithZeeGammaStudy/
  sleep 1
end


foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis//src/fileListEGamma2 | grep sel | grep HZZ `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/selectEleHZZMoriond2013WPWithZeeGammaStudy_${file}.out -J Efficiency_selectEleHZZMoriond2013WPWithZeeGammaStudy_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/Efficiency/ selectEleHZZHCP2012WPWithZeeGammaStudyPerFile.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"probes.${file}\",0,1,2\) probes.${file} /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Data2012_SingleEle_EleHZZMoriond2013WPWithZeeGammaStudy/
  sleep 1
end
