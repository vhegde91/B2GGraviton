#!/bin/sh    
root -l -q 'splitRunList.C("unskimmed_ZJetsToNuNu_HT_MC2016.txt",30)'
root -l -q 'splitRunList.C("unskimmed_ZJetsToNuNu_HT_MC2017.txt",30)'
root -l -q 'splitRunList.C("unskimmed_ZJetsToNuNu_HT_MC2018.txt",30)'
root -l -q 'splitRunList.C("unskimmed_WJetsToLNu_HT_MC2016.txt",30)'
root -l -q 'splitRunList.C("unskimmed_WJetsToLNu_HT_MC2017.txt",30)'
root -l -q 'splitRunList.C("unskimmed_WJetsToLNu_HT_MC2018.txt",30)'