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
		--outSAMtype BAM SortedByCoordinate --runThreadN 100
		echo ""
		echo ">>>>>>>>>>>>>>>>>>>> $i 1st Pass Alignment complete <<<<<<<<<<<<<<<<<<<<"
		echo ""
	done
