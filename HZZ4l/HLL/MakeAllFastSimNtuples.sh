#*************************************************************************************************
######################################################################
## 
## Make EfficiencyMap Ntuples
##
######################################################################
#*************************************************************************************************

######################################################################
## ZZ Bkg
######################################################################
root -l HiggsAna/HZZ4l/HLL/MakeEfficiencyMapNtuple.C+'("ZZTest",-1)'

root -l HiggsAna/HZZ4l/HLL/MakeEfficiencyMapNtuple.C+'("ZZ",0)'


######################################################################
## Signal
######################################################################
root -l HiggsAna/HZZ4l/HLL/MakeEfficiencyMapNtuple.C+'("HZZ120",120)'
root -l HiggsAna/HZZ4l/HLL/MakeEfficiencyMapNtuple.C+'("HZZ121",121)'
root -l HiggsAna/HZZ4l/HLL/MakeEfficiencyMapNtuple.C+'("HZZ122",122)'
root -l HiggsAna/HZZ4l/HLL/MakeEfficiencyMapNtuple.C+'("HZZ123",123)'
root -l HiggsAna/HZZ4l/HLL/MakeEfficiencyMapNtuple.C+'("HZZ124",124)'
root -l HiggsAna/HZZ4l/HLL/MakeEfficiencyMapNtuple.C+'("HZZ125",125)'
root -l HiggsAna/HZZ4l/HLL/MakeEfficiencyMapNtuple.C+'("HZZ126",126)'
root -l HiggsAna/HZZ4l/HLL/MakeEfficiencyMapNtuple.C+'("HZZ127",127)'
root -l HiggsAna/HZZ4l/HLL/MakeEfficiencyMapNtuple.C+'("HZZ130",130)'



#*************************************************************************************************
######################################################################
## 
## Make Fast Sim Validation Ntuples
##
######################################################################
#*************************************************************************************************




######################################################################
## ZZ Bkg
######################################################################
root -l HiggsAna/HZZ4l/HLL/MakeFastSimValidationNtuple.C+'("ZZ",0)'


######################################################################
## Signal
######################################################################
root -l HiggsAna/HZZ4l/HLL/MakeFastSimValidationNtuple.C+'("HZZ120",120)'
root -l HiggsAna/HZZ4l/HLL/MakeFastSimValidationNtuple.C+'("HZZ121",121)'
root -l HiggsAna/HZZ4l/HLL/MakeFastSimValidationNtuple.C+'("HZZ122",122)'
root -l HiggsAna/HZZ4l/HLL/MakeFastSimValidationNtuple.C+'("HZZ123",123)'
root -l HiggsAna/HZZ4l/HLL/MakeFastSimValidationNtuple.C+'("HZZ124",124)'
root -l HiggsAna/HZZ4l/HLL/MakeFastSimValidationNtuple.C+'("HZZ125",125)'
root -l HiggsAna/HZZ4l/HLL/MakeFastSimValidationNtuple.C+'("HZZ126",126)'
root -l HiggsAna/HZZ4l/HLL/MakeFastSimValidationNtuple.C+'("HZZ127",127)'
root -l HiggsAna/HZZ4l/HLL/MakeFastSimValidationNtuple.C+'("HZZ130",130)'

