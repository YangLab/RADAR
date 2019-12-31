
## start gfClient server
user_bin=/picb/rnomics4/rotation/fuzhican/bin;
ref_genome_path_blat=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/hg38_all.2bit
#Port=$(random 49152 65000)
Port=49253
${user_bin}/gfServer start 127.0.0.1 $Port -repMatch=2253 -stepSize=5 -log=./gfServer_untrans.log $ref_genome_path_blat



