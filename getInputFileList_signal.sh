#!/bin/sh
# for sample in mG1000 mG3000 mRad1000 mRad3000 mWprime1000 mWprime3000
# do
#     grep $sample a.txt | grep 2016 | grep PrivateSamples.VBF > unskimmed_vbf_${sample}_MC2016.txt
#     grep $sample a.txt | grep 2017 | grep PrivateSamples.VBF > unskimmed_vbf_${sample}_MC2017.txt
#     grep $sample a.txt | grep 2018 | grep PrivateSamples.VBF > unskimmed_vbf_${sample}_MC2018.txt
#     # echo "./signalRegGraviton unskimmed_vbf_${sample}_MC2016.txt unskimmed_WtFix_vbf_${sample}_MC2016.root Sig_MC_2016"
#     # echo "./signalRegGraviton unskimmed_vbf_${sample}_MC2017.txt unskimmed_WtFix_vbf_${sample}_MC2017.root Sig_MC_2017"
#     # echo "./signalRegGraviton unskimmed_vbf_${sample}_MC2018.txt unskimmed_WtFix_vbf_${sample}_MC2018.root Sig_MC_2018"
#     echo "root -l -q 'splitRunList.C(\"unskimmed_vbf_${sample}_MC2016.txt\",1)'"
#     echo "root -l -q 'splitRunList.C(\"unskimmed_vbf_${sample}_MC2017.txt\",1)'"
#     echo "root -l -q 'splitRunList.C(\"unskimmed_vbf_${sample}_MC2018.txt\",1)'"
# done
for sample in mG1000 mG3000 mRad1000 mRad3000 mWprime1000 mWprime3000
do
    grep $sample a.txt | grep 2016 | grep PrivateSamples.ggF > unskimmed_ggf_${sample}_MC2016.txt
    grep $sample a.txt | grep 2017 | grep PrivateSamples.ggF > unskimmed_ggf_${sample}_MC2017.txt
    grep $sample a.txt | grep 2018 | grep PrivateSamples.ggF > unskimmed_ggf_${sample}_MC2018.txt
    echo "root -l -q 'splitRunList.C(\"unskimmed_ggf_${sample}_MC2016.txt\",100)'"
    echo "root -l -q 'splitRunList.C(\"unskimmed_ggf_${sample}_MC2017.txt\",100)'"
    echo "root -l -q 'splitRunList.C(\"unskimmed_ggf_${sample}_MC2018.txt\",100)'"
    # echo "./signalRegGraviton unskimmed_ggf_${sample}_MC2016.txt unskimmed_WtFix_ggf_${sample}_MC2016.root Sig_MC_2016"
    # echo "./signalRegGraviton unskimmed_ggf_${sample}_MC2017.txt unskimmed_WtFix_ggf_${sample}_MC2017.root Sig_MC_2017"
    # echo "./signalRegGraviton unskimmed_ggf_${sample}_MC2018.txt unskimmed_WtFix_ggf_${sample}_MC2018.root Sig_MC_2018"
done
