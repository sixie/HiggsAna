################################################################
#Option 0 : stats toys
################################################################
foreach seed(`seq 1 1000`)
  echo $seed
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/ZeeGammaMassFitSystematicStudyToys_Option0_${seed}.out -J ZeeGammaMassFitSystematicStudyToys_Option0_${seed} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/ ZeeGammaMassFitSystematicStudy.C +\(\"/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/results/Data2012_EleHZZICHEP2012WPWithZeeGamma/basic2_76_106_LowPt/plots/FitWorkspaceFile_etapt_0.root\",${seed},0,1\) EffToyResults_Option0_Seed${seed}.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/ZeeGammaMassFitSystematicStudyToys/ /afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/rootlogon.C
end

################################################################
#Option 1 : mean -1
################################################################
foreach seed(`seq 1 1000`)
  echo $seed
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/ZeeGammaMassFitSystematicStudyToys_Option1_${seed}.out -J ZeeGammaMassFitSystematicStudyToys_Option1_${seed} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/ ZeeGammaMassFitSystematicStudy.C +\(\"/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/results/Data2012_EleHZZICHEP2012WPWithZeeGamma/basic2_76_106_LowPt/plots/FitWorkspaceFile_etapt_0.root\",${seed},1,1\) EffToyResults_Option1_Seed${seed}.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/ZeeGammaMassFitSystematicStudyToys/ /afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/rootlogon.C
end

################################################################
#Option 2 : mean +1
################################################################
foreach seed(`seq 1 1000`)
  echo $seed
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/ZeeGammaMassFitSystematicStudyToys_Option2_${seed}.out -J ZeeGammaMassFitSystematicStudyToys_Option2_${seed} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/ ZeeGammaMassFitSystematicStudy.C +\(\"/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/results/Data2012_EleHZZICHEP2012WPWithZeeGamma/basic2_76_106_LowPt/plots/FitWorkspaceFile_etapt_0.root\",${seed},2,1\) EffToyResults_Option2_Seed${seed}.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/ZeeGammaMassFitSystematicStudyToys/ /afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/rootlogon.C
end

################################################################
#Option 3 : sigma -1
################################################################
foreach seed(`seq 1 1000`)
  echo $seed
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/ZeeGammaMassFitSystematicStudyToys_Option3_${seed}.out -J ZeeGammaMassFitSystematicStudyToys_Option3_${seed} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/ ZeeGammaMassFitSystematicStudy.C +\(\"/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/results/Data2012_EleHZZICHEP2012WPWithZeeGamma/basic2_76_106_LowPt/plots/FitWorkspaceFile_etapt_0.root\",${seed},3,1\) EffToyResults_Option3_Seed${seed}.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/ZeeGammaMassFitSystematicStudyToys/ /afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/rootlogon.C
end

################################################################
#Option 4 : sigma +1
################################################################
foreach seed(`seq 1 1000`)
  echo $seed
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/ZeeGammaMassFitSystematicStudyToys_Option4_${seed}.out -J ZeeGammaMassFitSystematicStudyToys_Option4_${seed} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/ ZeeGammaMassFitSystematicStudy.C +\(\"/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/results/Data2012_EleHZZICHEP2012WPWithZeeGamma/basic2_76_106_LowPt/plots/FitWorkspaceFile_etapt_0.root\",${seed},4,1\) EffToyResults_Option4_Seed${seed}.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/ZeeGammaMassFitSystematicStudyToys/ /afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/rootlogon.C
end

################################################################
#Option 10 : stat toys from workspace
################################################################
foreach seed(`seq 1 1000`)
  echo $seed
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Efficiency/ZeeGammaMassFitSystematicStudyToys_Option10_${seed}.out -J ZeeGammaMassFitSystematicStudyToys_Option10_${seed} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob.csh /afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/ ZeeGammaMassFitSystematicStudy.C +\(\"/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/results/Data2012_EleHZZICHEP2012WPWithZeeGamma/basic2_76_106_LowPt/plots/FitWorkspaceFile_etapt_0.root\",${seed},10,1\) EffToyResults_Option10_Seed${seed}.root /afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/ZeeGammaMassFitSystematicStudyToys/ /afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/rootlogon.C
end

