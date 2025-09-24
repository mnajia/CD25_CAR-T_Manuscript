#! /bin/bash
#Build a custom kallisto index
#Mohamad Najia


output_dir=/broad/blainey_lab_storage/MN/20250122_AA_Tcell_in_vivo_SS2_ATAC_nextseq/CD25_CAR-T_Manuscript/rna-seq/kallisto_custom_index
index_name=kallisto_CAR_custom_index.idx
cdna=$output_dir/hg19.annot.cdna.CAR.custom
k=31

cd $output_dir
cat CAR.constructs.cdna.fa /broad/blainey_lab/Mo/genomes/hg19.annot.cdna > $cdna


qsub -N kallisto_index -o $output_dir/out.log -e $output_dir/err.log -m ea -M mnajia@broadinstitute.org /broad/blainey_lab/Mo/scripts/build_kallisto_index.sh $index_name $cdna $k
