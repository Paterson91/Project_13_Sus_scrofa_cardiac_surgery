#!/bin/bash

#Setting up some variables
#Location of STAR index file

Working_dir=/mnt/Scratch/Sus_scrofa_melanie/.
STAR_genome=/mnt/Scratch/Swine_flu/Sus_scrofa/STAR/
GTF_in=/mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic_3.gtf
threads_in=100


cd $Working_dir/fastq

############################################################################
############################### Untrimmed fastqc run #######################


mkdir -p fastqc_output/
fastqc *_sorted.fastq.gz -outdir=fastqc_output/ -t $threads_in

############################################################################
############################### Trim Script ################################

for i in `ls *_R1_sorted.fastq.gz | sed 's/_R1_sorted.fastq.gz//'`
	do
		echo ""
		echo ">>>>>>>>>>>>>>>>>>>> $i Trimming begun <<<<<<<<<<<<<<<<<<<<"
		echo ""
		/usr/bin/bbmap/bbduk.sh \
		in1=$i\_R1_sorted.fastq.gz \
		in2=$i\_R2_sorted.fastq.gz \
		out1=$i\_R1_trimmed.fastq.gz \
		out2=$i\_R2_trimmed.fastq.gz \
		ref=/usr/bin/bbmap/resources/adapters.fa \
		tpe tbo
		echo ">>>>>>>>>>>>>>>>>>> $i Trimming Complete <<<<<<<<<<<<<<<<<<<<"
		echo ""
	done
command; echo "QC Complete" | mail -s "Trimming Complete" a.paterson@bristol.ac.uk

############################################################################
############################# Trimmed fastqc run ###########################

mkdir -p fastqc_output_trimmed/
fastqc *trimmed.fastq.gz -outdir=fastqc_output_trimmed/ -t $threads_in

############################################################################
################ MultiQC Collation of both Trimmed and Untrimmed ###########

multiqc .
command; echo "QC Complete" | mail -s "QC Complete" a.paterson@bristol.ac.uk
cd ..

#INDEXING SUS SCROFA GENOME (NCBI)

if test -f "$STAR_genome"; then
  echo "STAR Indexed Genome exists. At file path $STAR_genome"
else
  echo "STAR Indexed Genome does not exist. Creating new..."
  STAR --runThreadN 50 \
  --runMode genomeGenerate \
  --genomeDir $STAR_genome \
  --genomeFastaFiles /mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic.fna \
  --sjdbGTFfile /mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic.gff \
  --sjdbOverhang 74 \
  --limitGenomeGenerateRAM=500000000000
fi

if [ $? -ne 0 ]
then
    command; echo "Indexing Error" | mail -s "Indexing Error" a.paterson@bristol.ac.uk
else
    command; echo "Indexing Complete" | mail -s "Indexing Complete" a.paterson@bristol.ac.uk
fi

#RUNNING ALIGNMENTS

mkdir -p STAR_output

for i in `find . -name "*_R1_trimmed.fastq.gz" | sed 's/_R1_trimmed.fastq.gz//'`
	do
		echo ""
		echo ">>>>>>>>>>>>>>>>>>>> $i 1st Pass Alignment begun <<<<<<<<<<<<<<<<<<<<"
		echo ""
		STAR --genomeDir /mnt/Scratch/Swine_flu/Sus_scrofa/STAR \
    --runMode alignReads \
		--readFilesIn $i\_R1_trimmed.fastq.gz $i\_R2_trimmed.fastq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix STAR_output/`basename $i` \
    --outSJfilterReads Unique \
		--outSAMtype BAM SortedByCoordinate --runThreadN $threads_in
		echo ""
		echo ">>>>>>>>>>>>>>>>>>>> $i 1st Pass Alignment complete <<<<<<<<<<<<<<<<<<<<"
		echo ""
	done

#Parsing GFF3 to GTF for featureCounts
#gffread GCF_000003025.6_Sscrofa11.1_genomic.gff -T -F -o GCF_000003025.6_Sscrofa11.1_genomic_2.gtf
#perl -ne 'chomp; @a=split/\t/; %h=split(/ /,$a[8]); $a[8]=join(" ",("gene_name",$h{"gene_name"},"transcript_id",$h{"transcript_id"})); \
#print join("\t",@a),"\n";' GCF_000003025.6_Sscrofa11.1_genomic_2.gtf > GCF_000003025.6_Sscrofa11.1_genomic_3.gtf

#COUNTING FEATURES

cd STAR_output

featureCounts \
-p -g gene_name -t exon \
-a $GTF_in \
-o Counts_stranded_gtf.gene -F GTF -T 64 \
-s 0 \
*.bam

#Tidy Output - Remove header, negative controls, strand information
tail -n +2 Counts_stranded_gtf.gene | cut -f1,7- | sed 's/Aligned.sortedByCoord.out.bam//g'> Counts_stranded_gtf.gene_tidied

#Vim to remove BAM extensions...
#CANT SED REMOVE BAM EXTENSION?
#tail -n +2 Counts_stranded_gtf.gene | cut -f1,7- | sed 's/-.*[^ ]//'> Counts_stranded_gtf.gene_tidied
#Differential Expression Analysis - DESeq2
#NOTE; Edit config.R to change options

#Rscript DESeq2.R config.R
