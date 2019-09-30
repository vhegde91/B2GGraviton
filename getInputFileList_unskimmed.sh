#!/bin/sh
#for sample in DYJetsToLL_M-50_HT GJets_DR-0p4_HT QCD_HT TTGJets_ TTJets_DiLept TTJets_HT TTJets_SingleLeptFromT WGJets_MonoPhoton_PtG WJetsToLNu_HT ZJetsToNuNu_HT
for sample in ZJetsToNuNu_HT WJetsToLNu_HT
do
    grep $sample /uscms/home/vhegde/inputFileList_V17.txt | grep Summer16 > unskimmed_${sample}_MC2016.txt
    echo "root -l -q 'splitRunList.C(\"unskimmed_${sample}_MC2016.txt\",30)'"

    grep $sample /uscms/home/vhegde/inputFileList_V17.txt | grep Fall17 > unskimmed_${sample}_MC2017.txt
    echo "root -l -q 'splitRunList.C(\"unskimmed_${sample}_MC2017.txt\",30)'"

    grep $sample /uscms/home/vhegde/inputFileList_V17.txt | grep Autumn18 > unskimmed_${sample}_MC2018.txt
    echo "root -l -q 'splitRunList.C(\"unskimmed_${sample}_MC2018.txt\",30)'"
done