#!/bin/sh
root -l -q 'findFailedJobs.C("WJetsToLNu_HT_MC2016")'
root -l -q 'findFailedJobs.C("WJetsToLNu_HT_MC2017")'
root -l -q 'findFailedJobs.C("WJetsToLNu_HT_MC2018")'
root -l -q 'findFailedJobs.C("ZJetsToNuNu_HT_MC2016")'
root -l -q 'findFailedJobs.C("ZJetsToNuNu_HT_MC2017")'
root -l -q 'findFailedJobs.C("ZJetsToNuNu_HT_MC2018")'

root -l -q 'findFailedJobs.C("TTJets_DiLept_MC2016")'
root -l -q 'findFailedJobs.C("TTJets_DiLept_MC2017")'
root -l -q 'findFailedJobs.C("TTJets_DiLept_MC2018")'
#root -l -q 'findFailedJobs.C("TTJets_HT_MC2016")'
#root -l -q 'findFailedJobs.C("TTJets_HT_MC2017")'
#root -l -q 'findFailedJobs.C("TTJets_HT_MC2018")'
root -l -q 'findFailedJobs.C("TTJets_SingleLeptFromT_MC2016")'
root -l -q 'findFailedJobs.C("TTJets_SingleLeptFromT_MC2017")'
root -l -q 'findFailedJobs.C("TTJets_SingleLeptFromT_MC2018")'

hadd -f TTJets_MC2016.root TTJets_DiLept_MC2016.root TTJets_SingleLeptFromT_MC2016.root 
hadd -f TTJets_MC2017.root TTJets_DiLept_MC2017.root TTJets_SingleLeptFromT_MC2017.root
hadd -f TTJets_MC2018.root TTJets_DiLept_MC2018.root TTJets_SingleLeptFromT_MC2018.root

root -l -q 'findFailedJobs.C("ST__MC2016")'
root -l -q 'findFailedJobs.C("ST__MC2017")'
root -l -q 'findFailedJobs.C("ST__MC2018")'

root -l -q 'findFailedJobs.C("MET_Run2016")'
root -l -q 'findFailedJobs.C("MET_Run2017")'
root -l -q 'findFailedJobs.C("MET_Run2018")'

#root -l -q 'findFailedJobs.C("WGJets_MonoPhoton_PtG_MC2016")'
#root -l -q 'findFailedJobs.C("WGJets_MonoPhoton_PtG_MC2017")'
#root -l -q 'findFailedJobs.C("WGJets_MonoPhoton_PtG_MC2018")'
#root -l -q 'findFailedJobs.C("DYJetsToLL_M-50_HT_MC2016")'
#root -l -q 'findFailedJobs.C("DYJetsToLL_M-50_HT_MC2017")'
#root -l -q 'findFailedJobs.C("DYJetsToLL_M-50_HT_MC2018")'
#root -l -q 'findFailedJobs.C("GJets_DR-0p4_HT_MC2016")'
#root -l -q 'findFailedJobs.C("GJets_DR-0p4_HT_MC2017")'
#root -l -q 'findFailedJobs.C("GJets_DR-0p4_HT_MC2018")'
#root -l -q 'findFailedJobs.C("QCD_HT_MC2016")'
#root -l -q 'findFailedJobs.C("QCD_HT_MC2017")'
#root -l -q 'findFailedJobs.C("QCD_HT_MC2018")'
#root -l -q 'findFailedJobs.C("TTGJets__MC2016")'
#root -l -q 'findFailedJobs.C("TTGJets__MC2017")'
#root -l -q 'findFailedJobs.C("TTGJets__MC2018")'
