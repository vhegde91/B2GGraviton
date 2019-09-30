#!/bin/sh    
root -l -q 'findFailedJobs.C("unskimmed_ZJetsToNuNu_HT_MC2016")'
root -l -q 'findFailedJobs.C("unskimmed_ZJetsToNuNu_HT_MC2017")'
root -l -q 'findFailedJobs.C("unskimmed_ZJetsToNuNu_HT_MC2018")'
root -l -q 'findFailedJobs.C("unskimmed_WJetsToLNu_HT_MC2016")'
root -l -q 'findFailedJobs.C("unskimmed_WJetsToLNu_HT_MC2017")'
root -l -q 'findFailedJobs.C("unskimmed_WJetsToLNu_HT_MC2018")'