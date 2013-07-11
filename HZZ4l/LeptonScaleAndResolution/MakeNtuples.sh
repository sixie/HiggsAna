############################################################################################
# Make Ntuples used for energy scale correction, and resolution smearing
############################################################################################



#MC
root -l HiggsAna/HZZ4l/LeptonScaleAndResolution/MakeZeeEventNtuples.C+\(\"root://eoscms//eos/cms/store/user/sixie/hist/AllNtuple/cern/filefi/028/AllNtuple_HZZ4lNtuple_s12-zllm50-2-v9_noskim_0000.root\",\"/data/blue/sixie/LeptonScaleAndResolution/Electrons/ZeeNtuple.AllNtuple_HZZ4lNtuple_s12-zllm50-2-v9_noskim_0000.root.root\",kTRUE,\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root\",0,2\)

foreach f ( `cat fileListHiggsHZZ4l2 | grep "HZZ4lNtuple" | grep s12-zllm50-2-v9`)
root -l -b -q HiggsAna/HZZ4l/LeptonScaleAndResolution/MakeZeeEventNtuples.C+\(\"root://eoscms//eos/cms/store/user/sixie/hist/AllNtuple/cern/filefi/028/$f\",\"/data/blue/sixie/LeptonScaleAndResolution/Electrons/ZeeNtuple.${f}.root\",kTRUE,\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root\",0,2\)
end



#2012 Data : Tight Plus Reco Skimmed
root -l HiggsAna/HZZ4l/LeptonScaleAndResolution/MakeZeeEventNtuples.C+\(\"/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-del-pr-v1_noskim_0000.root.TightPlusRecoSkimmed.root\",\"/data/blue/sixie/LeptonScaleAndResolution/Electrons/ZeeNtuple.AllNtuple_HZZ4lNtuple_r12a-del-pr-v1_noskim_0000.root.TightPlusRecoSkimmed.root.root\",kFALSE,\"\",0,2\)

foreach f ( `ls /data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-del*Tight* | sed 's/\/data\/smurf\/sixie\/hist\/HZZ4lNtuples\/data\/unprocessed\///' ` )
root -l -b -q HiggsAna/HZZ4l/LeptonScaleAndResolution/MakeZeeEventNtuples.C+\(\"/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/$f\",\"/data/blue/sixie/LeptonScaleAndResolution/Electrons/ZeeNtuple.${f}.root\",kFALSE,\"\",0,2\)
end
foreach f ( `ls /data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-del*Tight* | sed 's/\/data\/smurf\/sixie\/hist\/HZZ4lNtuples\/data\/unprocessed\///' ` )
root -l -b -q HiggsAna/HZZ4l/LeptonScaleAndResolution/MakeZeeEventNtuples.C+\(\"/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/$f\",\"/data/blue/sixie/LeptonScaleAndResolution/Electrons/ZeeNtuple.${f}.root\",kFALSE,\"\",0,2\)
end



#2012 Data : Non-skimmed
foreach f ( `cat fileListEGamma2 | grep "HZZ4lNtuple" | grep r12a-del` )
root -l -b -q HiggsAna/HZZ4l/LeptonScaleAndResolution/MakeZeeEventNtuples.C+\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/028/$f\",\"/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple.${f}.root\",kFALSE,\"\",0,2\)
end
foreach f ( `cat fileListEGamma2 | grep "HZZ4lNtuple" | grep r12b-del` )
root -l -b -q HiggsAna/HZZ4l/LeptonScaleAndResolution/MakeZeeEventNtuples.C+\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/028/$f\",\"/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple.${f}.root\",kFALSE,\"\",0,2\)
end

##############################################################
#52X 2012 Data: on lxbatch
##############################################################

#52X MC
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/fileListHiggsHZZ4l2 | grep "HZZ4lNtuple" | grep s12-zllm50-2-v9`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZeeEventNtuples/MakeZeeEventNtuples_${file}.out -J MakeZeeEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeZeeEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshzz4l/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"ZeeNtuple.${file}\",kTRUE,\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY_To_Run2012HCP.root\",0,1\) ZeeNtuple.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/
  sleep 1
end

#Run2012A
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12a-del-pr-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZeeEventNtuples/MakeZeeEventNtuples_${file}.out -J MakeZeeEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeZeeEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"ZeeNtuple.${file}\",kFALSE,\"\",0,1\) ZeeNtuple.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/
  sleep 1
end

#Run2012B
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12b-del-pr-v1`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZeeEventNtuples/MakeZeeEventNtuples_${file}.out -J MakeZeeEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeZeeEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"ZeeNtuple.${file}\",kFALSE,\"\",0,1\) ZeeNtuple.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/
  sleep 1
end


##############################################################
#2011 Data: on lxbatch
##############################################################

#Fall11 MC
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep f11-zeem20-powheg-v14b-bp `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZeeEventNtuples/MakeZeeEventNtuples_${file}.out -J MakeZeeEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeZeeEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"ZeeNtuple.${file}\",kTRUE,\"\",0,1\) ZeeNtuple.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/
  sleep 1
end

#Run2011A
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r11a-del-j16-v1-bp `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZeeEventNtuples/MakeZeeEventNtuples_${file}.out -J MakeZeeEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeZeeEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"ZeeNtuple.${file}\",kFALSE,\"\",0,1\) ZeeNtuple.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/
  sleep 1
end


#Run2011B
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r11b-del-j16-v1-bp `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZeeEventNtuples/MakeZeeEventNtuples_${file}.out -J MakeZeeEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeZeeEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"ZeeNtuple.${file}\",kFALSE,\"\",0,1\) ZeeNtuple.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/
  sleep 1
end


##############################################################
#53X 2012 Data: on lxbatch
##############################################################


#53X MC
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep s12-zllm50-v7a`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZeeEventNtuples/MakeZeeEventNtuples_${file}.out -J MakeZeeEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeZeeEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZeeNtuple.${file}\",kTRUE,\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_HCP2012.root\",0,2\) ZeeNtuple.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/
  sleep 1
end

#Run2012A
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12a-del-j13-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZeeEventNtuples/MakeZeeEventNtuples_${file}.out -J MakeZeeEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeZeeEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZeeNtuple.${file}\",kFALSE,\"\",0,2\) ZeeNtuple.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12a-del-a06-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZeeEventNtuples/MakeZeeEventNtuples_${file}.out -J MakeZeeEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeZeeEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZeeNtuple.${file}\",kFALSE,\"\",0,2\) ZeeNtuple.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/
  sleep 1
end

#Run2012B
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12b-del-j13-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZeeEventNtuples/MakeZeeEventNtuples_${file}.out -J MakeZeeEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeZeeEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZeeNtuple.${file}\",kFALSE,\"\",0,2\) ZeeNtuple.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/
  sleep 1
end

#Run2012C
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12c-del-a24-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZeeEventNtuples/MakeZeeEventNtuples_${file}.out -J MakeZeeEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeZeeEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZeeNtuple.${file}\",kFALSE,\"\",0,2\) ZeeNtuple.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12c-del-pr-v2 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZeeEventNtuples/MakeZeeEventNtuples_${file}.out -J MakeZeeEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeZeeEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZeeNtuple.${file}\",kFALSE,\"\",0,2\) ZeeNtuple.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/
  sleep 1
end

#Run2012D
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12d-del-pr-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZeeEventNtuples/MakeZeeEventNtuples_${file}.out -J MakeZeeEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeZeeEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZeeNtuple.${file}\",kFALSE,\"\",0,2\) ZeeNtuple.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/
  sleep 1
end

##############################################################
#52X 2012 Data Jun29 Rereco: on lxbatch
##############################################################


#Run2012A
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12a-del-j13-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZeeEventNtuples/MakeZeeEventNtuples_${file}.out -J MakeZeeEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeZeeEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZeeNtuple.${file}\",kFALSE,\"\",0,2\) ZeeNtuple.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12a-del-a06-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZeeEventNtuples/MakeZeeEventNtuples_${file}.out -J MakeZeeEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeZeeEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZeeNtuple.${file}\",kFALSE,\"\",0,2\) ZeeNtuple.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/
  sleep 1
end



############################################################################################
# Make Ntuples used for training energy regression
#
# Submit to lxbatch
#
############################################################################################

##Higgs -> ZZ signal

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListHiggsHZZ4l2 | grep HZZ4lNtuple |  grep "h...zz" `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromMC_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromMC_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromMC.C +\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshzz4l/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.${file}\",0\) ElectronNtuples.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end

## 52X: Z->ee
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListHiggsHZZ4l2 | grep HZZ4lNtuple |  grep s12-zllm50-2-v9 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromMC_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromMC_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromMC.C +\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshzz4l/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"ElectronNtuples.${file}\",0\) ElectronNtuples.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end

## 53X: Z->ee
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple |  grep s12-zllm50-v7a `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromMC_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromMC_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromMC.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.${file}\",0\) ElectronNtuples.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end


## ZZ
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListHiggsHZZ4l2 | grep HZZ4lNtuple | grep powheg | grep -v 4m `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromMC_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromMC_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromMC.C +\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshzz4l/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"ElectronNtuples.${file}\",0\) ElectronNtuples.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end



## Zee MC
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple |  grep s12-zllm50-v7a`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromMC_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromMC_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeMC.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.ZeeTP.${file}\",0\) ElectronNtuples.ZeeTP.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end


## Zee Data
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12a-del-j13-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromData_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromData_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeData.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.${file}\",0\) ElectronNtuples.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12a-del-a06-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromData_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromData_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeData.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.${file}\",0\) ElectronNtuples.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12b-del-j13-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromData_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromData_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeData.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.${file}\",0\) ElectronNtuples.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12c-del-a24-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromData_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromData_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeData.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.${file}\",0\) ElectronNtuples.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12c-del-pr-v2 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromData_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromData_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeData.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.${file}\",0\) ElectronNtuples.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12d-del-pr-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromData_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromData_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeData.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.${file}\",0\) ElectronNtuples.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end


## ZeeGamma MC
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple |  grep s12-zllm50-v7a `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromMC_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromMC_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeGammaMC.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.ZeeGamma.${file}\",-1\) ElectronNtuples.ZeeGamma.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end


## ZeeGamma Data
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12a-del `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromData_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromData_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeGammaData.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.ZeeGamma.${file}\",0\) ElectronNtuples.ZeeGamma.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12b-del `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromData_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromData_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeGammaData.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.ZeeGamma.${file}\",0\) ElectronNtuples.ZeeGamma.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12c-del-a24-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromData_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromData_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeGammaData.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.ZeeGamma.${file}\",0\) ElectronNtuples.ZeeGamma.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12c-del-pr-v2 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromData_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromData_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeGammaData.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.ZeeGamma.${file}\",0\) ElectronNtuples.ZeeGamma.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12d-del-pr-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromData_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromData_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeGammaData.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.ZeeGamma.${file}\",0\) ElectronNtuples.ZeeGamma.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end



foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12a-sel `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromData_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromData_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeGammaData.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.ZeeGamma.${file}\",1\) ElectronNtuples.ZeeGamma.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12b-sel `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromData_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromData_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeGammaData.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.ZeeGamma.${file}\",1\) ElectronNtuples.ZeeGamma.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12c-sel-a24-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromData_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromData_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeGammaData.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.ZeeGamma.${file}\",1\) ElectronNtuples.ZeeGamma.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12c-sel-pr-v2 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromData_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromData_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeGammaData.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.ZeeGamma.${file}\",1\) ElectronNtuples.ZeeGamma.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep r12d-sel-pr-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromData_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromData_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromZeeGammaData.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ElectronNtuples.ZeeGamma.${file}\",1\) ElectronNtuples.ZeeGamma.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end




#Merge
hadd -f /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.Summer12DY.root /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.*s12-zllm50-2-v9*.root

hadd -f /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.Summer12HZZ1254l.root /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.*s12-h125zz4l-gf-v9*.root  

hadd -f /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.Summer12HZZ4l.root /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.*s12-h120zz4l-gf-v9*.root  /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.*s12-h121zz4l-gf-v9*.root  /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.*s12-h123zz4l-gf-v9*.root  /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.*s12-h124zz4l-gf-v9*.root  /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.*s12-h127zz4l-gf-v9*.root  /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.*s12-h130zz4l-gf-v9*.root 

hadd -f /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.Summer12ZZ.root /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.*s12-zz2e2m-powheg-v9*.root /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.*s12-zz4e-powheg-v9*.root


############################################################################################
# Make Ntuples used for training energy regression
#
# 2011 Data and MC
#
############################################################################################

#Zee MC

foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/fileListEGamma2 | grep HZZ4lNtuple |  grep f11-zeem20-powheg-v14b-bp`)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/LeptonScaleAndResolution/MakeElectronNtupleFromMC_${file}.out -J LeptonScaleAndResolution_MakeElectronNtupleFromMC_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/user/s/sixie/CMSSW_analysis/src/HiggsAna/HZZ4l/LeptonScaleAndResolution/ MakeElectronNtupleFromMC.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/025/${file}\",\"ElectronNtuples.${file}\",0\) ElectronNtuples.${file} /afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/
  sleep 1
end
