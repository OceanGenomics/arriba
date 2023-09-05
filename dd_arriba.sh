sample=$1


/data008/users/yshen/arriba_debug/arriba -x /data008/users/yshen/benchmark_results/CTAT_FUSIONTRANS/STARfusion_results/$sample/Aligned.out.bam \
       -o /data008/users/yshen/arriba_debug/output3/$sample/fusions.tsv \
       -a /data008/users/yshen/arriba_v2.4.0/hg19/hs37d5viral.fa \
       -g /data008/users/yshen/arriba_v2.4.0/hg19/GENCODE19.gtf \
       -b /data008/users/yshen/arriba_v2.4.0/database/blacklist_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz \
       -k /data008/users/yshen/arriba_v2.4.0/database/known_fusions_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz \
       -t /data008/users/yshen/arriba_v2.4.0/database/known_fusions_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz \
       -p /data008/users/yshen/arriba_v2.4.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3 \
       #-f homologs
       #-X 
       #-O /data008/users/yshen/arriba_debug/test_output/fusions.discarded.tsv \
