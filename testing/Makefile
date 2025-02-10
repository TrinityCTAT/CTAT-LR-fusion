
simple: test_excl_shortreads

both:  test_incl_shortreads test_excl_shortreads

all_long_read_tests: test_excl_shortreads test_excl_shortreads_set_num_total_reads  test_excl_shortreads_bam_input test_oncoreads test_oncoreads_chim_candidates_only


all: all_long_read_tests test_incl_shortreads  test_longOnly_fqgz test_longOnly_fq test_incl_coding_effect


test_excl_shortreads:
	../ctat-LR-fusion -T transcripts.fa --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short --min_J 1  --min_sumJS 1 --vis --extract_fusion_LR_fasta fusion_ev_LR.reads.fa --min_trans_overlap_length 75

test_excl_shortreads_reg_mm2:
	../ctat-LR-fusion -T transcripts.fa --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short.reg_mm2 --min_J 1  --min_sumJS 1 --vis --extract_fusion_LR_fasta fusion_ev_LR.reads.fa --min_trans_overlap_length 75 --no_ctat_mm2


test_excl_shortreads_bam_input:
	../ctat-LR-fusion --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short.LR_bam --min_J 1  --min_sumJS 1 --vis --extract_fusion_LR_fasta fusion_ev_LR.reads.fa --min_trans_overlap_length 75 --LR_bam transcripts.mm2.bam

test_excl_shortreads_set_num_total_reads:
	../ctat-LR-fusion -T transcripts.fa --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short.setNumReads --min_J 1  --min_sumJS 1 --vis --extract_fusion_LR_fasta fusion_ev_LR.reads.fa --min_trans_overlap_length 75 --num_total_reads 500000



test_incl_shortreads:
	../ctat-LR-fusion -T transcripts.fa --left_fq reads_1.fq.gz --right_fq reads_2.fq.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.incl_short --vis

test_incl_shortreads_SEonly:
	../ctat-LR-fusion -T transcripts.fa --left_fq reads_1.fq.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.incl_short_SEonly --vis


test_incl_coding_effect: test_excl_shortreads
	../ctat-LR-fusion -T transcripts.fa --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short --min_J 1  --min_sumJS 1 --vis --extract_fusion_LR_fasta fusion_ev_LR.reads.fa --examine_coding_effect


test_incl_FI_params:
	../ctat-LR-fusion -T transcripts.fa --left_fq reads_1.fq.gz --right_fq reads_2.fq.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.incl_short_custFI --vis  --FI_extra_params "--STAR_xtra_params '--limitBAMsortRAM 51419850732' " 


test_longOnly_fqgz:
	../ctat-LR-fusion -T transcripts.fq.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short.fqgz --min_J 1  --min_sumJS 1 --min_novel_junction_support 1

test_longOnly_fq: test_excl_shortreads


test_oncoreads:
	../ctat-LR-fusion -T oncofusionreads.fasta.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.oncoreads --vis

test_oncoreads_chim_candidates_only:
	../ctat-LR-fusion -T oncofusionreads.fasta.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.oncoreads.chim_candidates_only --chim_candidates_only


test_docker_LR_only:
	../Docker/run_docker.sh ctat-LR-fusion -T /data/transcripts.fa --genome_lib_dir /ctat_genome_lib -o ctat_LR_fusion_outdir.LR-only.docker --vis --min_trans_overlap_length 75


test_docker_LR_and_SR:
	../Docker/run_docker.sh ctat-LR-fusion -T /data/transcripts.fa --genome_lib_dir /ctat_genome_lib --min_J 1  --min_sumJS 1 --left_fq reads_1.fq.gz --right_fq reads_2.fq.gz -o ctat_LR_fusion_outdir.LR_and_SR.docker --vis --min_trans_overlap_length 75


test_singularity:
	singularity exec -e -B `pwd`:/data -B ${CTAT_GENOME_LIB}:/ctat_genome_lib ../Docker/ctat_lr_fusion.simg ctat-LR-fusion -T /data/transcripts.fa --genome_lib_dir /ctat_genome_lib --min_J 1  --min_sumJS 1 --left_fq reads_1.fq.gz --right_fq reads_2.fq.gz -o ctat_LR_fusion_outdir.singularity --vis

clean:
	rm -rf ctat_LR_fusion_outdir.* 
	sudo rm -rf ./ctat_LR_fusion_outdir.*.docker





