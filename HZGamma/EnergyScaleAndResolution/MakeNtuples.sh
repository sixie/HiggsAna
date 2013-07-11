############################################################################################
# Make Ntuples used for energy scale correction, and resolution smearing
############################################################################################


##############################################################
#53X MC
##############################################################
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep s12-zllm50-v7a `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZmumuGammaEventNtuples/MakeZmumuGammaEventNtuples_${file}.out -J MakeZmumuGammaEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/HiggsAna/HZGamma/EnergyScaleAndResolution/ MakeZmumuGammaEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZmumuGammaNtuple.${file}\",kTRUE,\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\",-1,2\) ZmumuGammaNtuple.${file} /afs/cern.ch/work/s/sixie/public/EnergyScaleAndResolution/ZmumuGammaEvents/
  sleep 1
end

foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/fileListEGamma2 | grep "HZZ4lNtuple" | grep s12-zllm1050-v7a `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZmumuGammaEventNtuples/MakeZmumuGammaEventNtuples_${file}.out -J MakeZmumuGammaEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/HiggsAna/HZGamma/EnergyScaleAndResolution/ MakeZmumuGammaEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_egamma/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZmumuGammaNtuple.${file}\",kTRUE,\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY53X_To_Full2012.root\",-2,2\) ZmumuGammaNtuple.${file} /afs/cern.ch/work/s/sixie/public/EnergyScaleAndResolution/ZmumuGammaEvents/
  sleep 1
end


##############################################################
#52X MC
##############################################################
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/fileListHiggsHZZ4l2 | grep "HZZ4lNtuple" | grep s12-zllm50-2-v9 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZmumuGammaEventNtuples/MakeZmumuGammaEventNtuples_${file}.out -J MakeZmumuGammaEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/HiggsAna/HZGamma/EnergyScaleAndResolution/ MakeZmumuGammaEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshzz4l/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"ZmumuGammaNtuple.${file}\",kTRUE,\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY_To_Run2012AB.root\",-1,2\) ZmumuGammaNtuple.${file} /afs/cern.ch/work/s/sixie/public/EnergyScaleAndResolution/ZmumuGammaEvents/
  sleep 1
end

foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/fileListHiggsHZZ4l2 | grep "HZZ4lNtuple" | grep s12-zllm1050-v9 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZmumuGammaEventNtuples/MakeZmumuGammaEventNtuples_${file}.out -J MakeZmumuGammaEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/HiggsAna/HZGamma/EnergyScaleAndResolution/ MakeZmumuGammaEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshzz4l/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"ZmumuGammaNtuple.${file}\",kTRUE,\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/PileupReweighting.Summer12DY_To_Run2012AB.root\",-2,2\) ZmumuGammaNtuple.${file} /afs/cern.ch/work/s/sixie/public/EnergyScaleAndResolution/ZmumuGammaEvents/
  sleep 1
end


##############################################################
# 52X Data
##############################################################

#Run2012A
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/fileListMuon2 | grep "HZZ4lNtuple" | grep r12a-dmu-j29 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZmumuGammaEventNtuples/MakeZmumuGammaEventNtuples_${file}.out -J MakeZmumuGammaEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/HiggsAna/HZGamma/EnergyScaleAndResolution/ MakeZmumuGammaEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_muon/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"ZmumuGammaNtuple.${file}\",kFALSE,\"\",0,2\) ZmumuGammaNtuple.${file} /afs/cern.ch/work/s/sixie/public/EnergyScaleAndResolution/ZmumuGammaEvents/
  sleep 1
end

#Run2012B
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/fileListMuon2 | grep "HZZ4lNtuple" | grep r12b-dmu-j29-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZmumuGammaEventNtuples/MakeZmumuGammaEventNtuples_${file}.out -J MakeZmumuGammaEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/HiggsAna/HZGamma/EnergyScaleAndResolution/ MakeZmumuGammaEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_muon/sixie/hist/AllNtuple/cern/filefi/028/${file}\",\"ZmumuGammaNtuple.${file}\",kFALSE,\"\",0,2\) ZmumuGammaNtuple.${file} /afs/cern.ch/work/s/sixie/public/EnergyScaleAndResolution/ZmumuGammaEvents/
  sleep 1
end


##############################################################
# 53X Data
##############################################################

#Run2012A
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/fileListMuon2 | grep "HZZ4lNtuple" | grep r12a-dmu-j13-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZmumuGammaEventNtuples/MakeZmumuGammaEventNtuples_${file}.out -J MakeZmumuGammaEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/HiggsAna/HZGamma/EnergyScaleAndResolution/ MakeZmumuGammaEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_muon/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZmumuGammaNtuple.${file}\",kFALSE,\"\",0,2\) ZmumuGammaNtuple.${file} /afs/cern.ch/work/s/sixie/public/EnergyScaleAndResolution/ZmumuGammaEvents/
  sleep 1
end
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/fileListMuon2 | grep "HZZ4lNtuple" | grep r12a-dmu-a06-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZmumuGammaEventNtuples/MakeZmumuGammaEventNtuples_${file}.out -J MakeZmumuGammaEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/HiggsAna/HZGamma/EnergyScaleAndResolution/ MakeZmumuGammaEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_muon/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZmumuGammaNtuple.${file}\",kFALSE,\"\",0,2\) ZmumuGammaNtuple.${file} /afs/cern.ch/work/s/sixie/public/EnergyScaleAndResolution/ZmumuGammaEvents/
  sleep 1
end

#Run2012B
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/fileListMuon2 | grep "HZZ4lNtuple" | grep r12b-dmu-j13-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZmumuGammaEventNtuples/MakeZmumuGammaEventNtuples_${file}.out -J MakeZmumuGammaEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/HiggsAna/HZGamma/EnergyScaleAndResolution/ MakeZmumuGammaEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_muon/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZmumuGammaNtuple.${file}\",kFALSE,\"\",0,2\) ZmumuGammaNtuple.${file} /afs/cern.ch/work/s/sixie/public/EnergyScaleAndResolution/ZmumuGammaEvents/
  sleep 1
end

#Run2012C
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/fileListMuon2 | grep "HZZ4lNtuple" | grep r12c-dmu-a24-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZmumuGammaEventNtuples/MakeZmumuGammaEventNtuples_${file}.out -J MakeZmumuGammaEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/HiggsAna/HZGamma/EnergyScaleAndResolution/ MakeZmumuGammaEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_muon/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZmumuGammaNtuple.${file}\",kFALSE,\"\",0,2\) ZmumuGammaNtuple.${file} /afs/cern.ch/work/s/sixie/public/EnergyScaleAndResolution/ZmumuGammaEvents/
  sleep 1
end
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/fileListMuon2 | grep "HZZ4lNtuple" | grep r12c-dmu-pr-v2 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZmumuGammaEventNtuples/MakeZmumuGammaEventNtuples_${file}.out -J MakeZmumuGammaEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/HiggsAna/HZGamma/EnergyScaleAndResolution/ MakeZmumuGammaEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_muon/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZmumuGammaNtuple.${file}\",kFALSE,\"\",0,2\) ZmumuGammaNtuple.${file} /afs/cern.ch/work/s/sixie/public/EnergyScaleAndResolution/ZmumuGammaEvents/
  sleep 1
end

#Run2012D
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/fileListMuon2 | grep "HZZ4lNtuple" | grep r12d-dmu-pr-v1 `)
  echo $file
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/MakeZmumuGammaEventNtuples/MakeZmumuGammaEventNtuples_${file}.out -J MakeZmumuGammaEventNtuples_${file} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_3_patch3/src/HiggsAna/HZGamma/EnergyScaleAndResolution/ MakeZmumuGammaEventNtuples.C +\(\"root://eoscms//eos/cms/store/group/phys_muon/sixie/hist/AllNtuple/cern/filefi/029/${file}\",\"ZmumuGammaNtuple.${file}\",kFALSE,\"\",0,2\) ZmumuGammaNtuple.${file} /afs/cern.ch/work/s/sixie/public/EnergyScaleAndResolution/ZmumuGammaEvents/
  sleep 1
end




############################################################################################
# Make Ntuples used for training energy regression
#
# Submit to lxbatch
#
############################################################################################

