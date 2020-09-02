featureCounts \
-p -g gene_name -t exon \
-a /mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic_3.gtf \
-o Counts_stranded_gtf.gene -F GTF -T 64 \
-s 0 \
*.bam
