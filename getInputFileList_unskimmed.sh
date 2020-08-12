#!/bin/sh
#for sample in DYJetsToLL_M-50_HT GJets_DR-0p4_HT QCD_HT TTGJets_ TTJets_DiLept TTJets_HT TTJets_SingleLeptFromT WGJets_MonoPhoton_PtG WJetsToLNu_HT ZJetsToNuNu_HT
#for sample in ZJetsToNuNu_HT WJetsToLNu_HT
for sample in TTJets_DiLept TTJets_HT TTJets_SingleLeptFromT
do
    grep $sample /uscms/home/vhegde/inputFileList_V17.txt | grep Summer16 > unskimmed_${sample}_MC2016.txt
    echo "root -l -q 'splitRunList.C(\"unskimmed_${sample}_MC2016.txt\",30)'"

    grep $sample /uscms/home/vhegde/inputFileList_V17.txt | grep Fall17 > unskimmed_${sample}_MC2017.txt
    echo "root -l -q 'splitRunList.C(\"unskimmed_${sample}_MC2017.txt\",30)'"

    grep $sample /uscms/home/vhegde/inputFileList_V17.txt | grep Autumn18 > unskimmed_${sample}_MC2018.txt
    echo "root -l -q 'splitRunList.C(\"unskimmed_${sample}_MC2018.txt\",30)'"
done

rm unskimmed_Rare_MC2016.txt unskimmed_Rare_MC2017.txt unskimmed_Rare_MC2018.txt

for sample in Summer16v3.WWTo1L1Nu2Q_ Summer16v3.WWZ_ Summer16v3.WZTo1L1Nu2Q_ Summer16v3.WZTo1L3Nu_ Summer16v3.WZZ_ Summer16v3.ZZTo2L2Q_ Summer16v3.ZZTo2Q2Nu_ Summer16v3.ZZZ_ Summer16v3.TTTT_ Summer16v3.TTWJetsToLNu_ Summer16v3.TTWJetsToQQ_ Summer16v3.TTGJets_ Summer16v3.TTZToLLNuNu_ Summer16v3.TTZToQQ
do
    grep $sample /uscms/home/vhegde/inputFileList_V17.txt | grep Summer16 >> unskimmed_Rare_MC2016.txt
    echo "root -l -q 'splitRunList.C(\"unskimmed_Rare_MC2016.txt\",30)'"
done

for sample in Fall17.WWTo1L1Nu2Q_ Summer16v3.WWZ_ Fall17.WZTo1L1Nu2Q_ Fall17.WZTo1L3Nu_ Fall17.WZZ_ Fall17.ZZTo2L2Q_ Summer16v3.ZZTo2Q2Nu_ Fall17.ZZZ_ Fall17.TTTT_ Fall17.TTWJetsToLNu_ Fall17.TTWJetsToQQ_ Fall17.TTGJets_ Fall17.TTZToLLNuNu_ Fall17.TTZToQQ_
do
    grep $sample /uscms/home/vhegde/inputFileList_V17.txt >> unskimmed_Rare_MC2017.txt
    echo "root -l -q 'splitRunList.C(\"unskimmed_Rare_MC2017.txt\",30)'"
done

for sample in Autumn18.WWTo1L1Nu2Q_ Summer16v3.WWZ_ Fall17.WZTo1L1Nu2Q_ Autumn18.WZTo1L3Nu_ Fall17.WZZ_ Autumn18.ZZTo2L2Q_ Summer16v3.ZZTo2Q2Nu_ Fall17.ZZZ_ Fall17.TTTT_ Autumn18.TTWJetsToLNu_ Autumn18.TTWJetsToQQ_ Autumn18.TTGJets_ Autumn18.TTZToLLNuNu_ Autumn18.TTZToQQ_
do
    grep $sample /uscms/home/vhegde/inputFileList_V17.txt >> unskimmed_Rare_MC2018.txt
    echo "root -l -q 'splitRunList.C(\"unskimmed_Rare_MC2018.txt\",30)'"
done
