#! /bin/bash
#Count ATAC-seq reads in consensus ATAC peaks
#Mohamad Najia


#declare variables 
project_dir=/broad/blainey_lab_storage/MN/20250122_AA_Tcell_in_vivo_SS2_ATAC_nextseq/CD25_CAR-T_Manuscript/atac-seq
countsDir=$project_dir/counts
peakFile=$project_dir/AA_PB_T_in_vivo_ATAC_hg38_consensusPeaks_minSamples2_minScore1.narrowPeak
samplesheet=$project_dir/samplesheet.txt

cd $countsDir


while IFS=$'\t' read	SAMPLE	BAM
do
    echo "Submitting job: $SAMPLE"

    outputFile=$countsDir/${SAMPLE}_counts.bed

    qsub -N $SAMPLE -j y -o $countsDir/$SAMPLE.counts.log $project_dir/02_count_reads_in_peaks.sh $BAM $peakFile $outputFile

done < $samplesheet






