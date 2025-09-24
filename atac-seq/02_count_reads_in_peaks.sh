#! /bin/bash
#$ -l h_vmem=5g
#$ -cwd
#$ -q broad
#$ -l h_rt=1:00:00

#ABC Model Pipeline
#Mohamad Najia


source /broad/software/scripts/useuse 

use BEDTools


bamFile=$1
peakFile=$2
outputFile=$3

bedtools multicov -bams $bamFile -bed $peakFile > $outputFile
