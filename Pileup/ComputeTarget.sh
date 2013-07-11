##########################
#Check out proper code
##########################
cvs co -r V03-03-16 RecoLuminosity/LumiDB
scram b -k
cd RecoLuminosity/LumiDB/scripts



##########################
#Compute Target pileup true
##########################

#full 2012
pileupCalc.py -i /afs/cern.ch/work/s/sixie/public/HZZ4l/auxiliar/2012/Cert_190456-200601_8TeV_PromptReco_Collisions12_JSON_v2.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-201678_corr.txt  --calcMode true  --minBiasXsec 69400 --maxPileupBin 70 PileupTarget_190456To200601.true.root 

#Run2012A/B
pileupCalc.py -i /afs/cern.ch/work/s/sixie/public/HZZ4l/auxiliar/2012/Cert_190456-196531_8TeV_PromptReco_Collisions12_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-201678_corr.txt  --calcMode true  --minBiasXsec 69400 --maxPileupBin 70 PileupTarget_190456To196531.true.root 

#only Run2012C
pileupCalc.py -i /afs/cern.ch/work/s/sixie/public/HZZ4l/auxiliar/2012/Cert_190456-200601_8TeV_PromptReco_Collisions12_JSON_v2.198049-200601.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-201678_corr.txt  --calcMode true  --minBiasXsec 69400 --maxPileupBin 70 PileupTarget_198049To200601.true.root 


##########################
#Compute Target pileup observed
##########################

#Run2012A/B
root -l -b -q HiggsAna/Pileup/MakePileupObserved.C+'("PileupTarget_190456To196531.true.69400.root","PileupTarget_190456To196531.obs.69400.root")'

#Run2012C
root -l -b -q HiggsAna/Pileup/MakePileupObserved.C+'("PileupTarget_198049To200601.true.69400.root","PileupTarget_198049To200601.obs.69400.root")'

#Run2012ABC
root -l -b -q HiggsAna/Pileup/MakePileupObserved.C+'("PileupTarget_190456To200601.true.69400.root","PileupTarget_190456To200601.obs.69400.root")'


##########################
# Compute Target pileup observed
##########################
root -l HiggsAna/Pileup/PileupReweighting.C+'(-1)'
