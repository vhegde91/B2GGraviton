#!/bin/sh    

root -l -q 'splitRunList.C("WJetsToLNu_HT_MC2016.txt",1)'
root -l -q 'splitRunList.C("WJetsToLNu_HT_MC2017.txt",1)'
root -l -q 'splitRunList.C("WJetsToLNu_HT_MC2018.txt",1)'

root -l -q 'splitRunList.C("ZJetsToNuNu_HT_MC2016.txt",1)'
root -l -q 'splitRunList.C("ZJetsToNuNu_HT_MC2017.txt",1)'
root -l -q 'splitRunList.C("ZJetsToNuNu_HT_MC2018.txt",1)'

root -l -q 'splitRunList.C("TTJets_DiLept_MC2016.txt",1)'
root -l -q 'splitRunList.C("TTJets_DiLept_MC2017.txt",1)'
root -l -q 'splitRunList.C("TTJets_DiLept_MC2018.txt",1)'
root -l -q 'splitRunList.C("TTJets_HT_MC2016.txt",1)'
root -l -q 'splitRunList.C("TTJets_HT_MC2017.txt",1)'
root -l -q 'splitRunList.C("TTJets_HT_MC2018.txt",1)'
root -l -q 'splitRunList.C("TTJets_SingleLeptFromT_MC2016.txt",1)'
root -l -q 'splitRunList.C("TTJets_SingleLeptFromT_MC2017.txt",1)'
root -l -q 'splitRunList.C("TTJets_SingleLeptFromT_MC2018.txt",1)'

root -l -q 'splitRunList.C("ST__MC2016.txt",1)'
root -l -q 'splitRunList.C("ST__MC2017.txt",1)'
root -l -q 'splitRunList.C("ST__MC2018.txt",1)'

root -l -q 'splitRunList.C("RareBG_MC2016.txt",1)'
root -l -q 'splitRunList.C("RareBG_MC2017.txt",1)'
root -l -q 'splitRunList.C("RareBG_MC2018.txt",1)'

root -l -q 'splitRunList.C("vbfGrav_1200_MC2016.txt",1)'
root -l -q 'splitRunList.C("vbfGrav_800_MC2017.txt",1)'
root -l -q 'splitRunList.C("vbfGrav_800_MC2018.txt",1)'


#root -l -q 'splitRunList.C("MET_Run2016.txt",1)'
#root -l -q 'splitRunList.C("MET_Run2017.txt",1)'
#root -l -q 'splitRunList.C("MET_Run2018.txt",1)'

#root -l -q 'splitRunList.C("WGJets_MonoPhoton_PtG_MC2017.txt",1)'
#root -l -q 'splitRunList.C("WGJets_MonoPhoton_PtG_MC2017.txt",1)'
#root -l -q 'splitRunList.C("WGJets_MonoPhoton_PtG_MC2018.txt",1)'
#root -l -q 'splitRunList.C("DYJetsToLL_M-50_HT_MC2017.txt",1)'
#root -l -q 'splitRunList.C("DYJetsToLL_M-50_HT_MC2017.txt",1)'
#root -l -q 'splitRunList.C("DYJetsToLL_M-50_HT_MC2018.txt",1)'
#root -l -q 'splitRunList.C("GJets_DR-0p4_HT_MC2017.txt",1)'
#root -l -q 'splitRunList.C("GJets_DR-0p4_HT_MC2017.txt",1)'
#root -l -q 'splitRunList.C("GJets_DR-0p4_HT_MC2018.txt",1)'
#root -l -q 'splitRunList.C("QCD_HT_MC2017.txt",1)'
#root -l -q 'splitRunList.C("QCD_HT_MC2017.txt",1)'
#root -l -q 'splitRunList.C("QCD_HT_MC2018.txt",1)'
#root -l -q 'splitRunList.C("TTGJets__MC2017.txt",1)'
#root -l -q 'splitRunList.C("TTGJets__MC2017.txt",1)'
#root -l -q 'splitRunList.C("TTGJets__MC2018.txt",1)'
