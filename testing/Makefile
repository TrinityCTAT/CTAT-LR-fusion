
simple: test_excl_shortreads

both:  test_incl_shortreads test_excl_shortreads

basic_long_read_tests: test_longOnly_fa test_longOnly_fq test_longOnly_fqgz test_longOnly_fagz test_longOnly_fa_nophase2

all_long_read_tests: test_excl_shortreads test_excl_shortreads_set_num_total_reads  test_excl_shortreads_bam_input test_oncoreads test_oncoreads_chim_candidates_only test_longOnly_fa test_longOnly_fq test_longOnly_fqgz test_longOnly_fagz test_longOnly_fa_nophase2 test_longOnly_fa_max_rigor test_longOnly_fa_incl_targets test_longOnly_fa_only_targets test_single_cell


all: all_long_read_tests test_incl_shortreads  test_longOnly_fqgz test_longOnly_fq test_incl_coding_effect

#########################################################
####################### No short reads ##################

test_excl_shortreads:
	../ctat-LR-fusion -T transcripts.fa --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short  --vis --extract_fusion_LR_fasta fusion_ev_LR.reads.fa --min_trans_overlap_length 75

test_longOnly_fq:
	../ctat-LR-fusion -T transcripts.fq --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short.fq  --min_novel_junction_support 1 --min_trans_overlap_length 75

test_longOnly_fqgz:
	../ctat-LR-fusion -T transcripts.fq.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short.fqgz --min_novel_junction_support 1 --min_trans_overlap_length 75

test_longOnly_fa: test_excl_shortreads

test_longOnly_fagz: 
	../ctat-LR-fusion -T transcripts.fa.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short.fagz --min_novel_junction_support 1 --min_trans_overlap_length 75


test_longOnly_fa_max_rigor:
	../ctat-LR-fusion -T transcripts.fa --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.LRonly_max_rigor  --vis --extract_fusion_LR_fasta fusion_ev_LR.reads.fa --min_trans_overlap_length 75 --max_rigor


test_longOnly_fa_incl_targets:
	../ctat-LR-fusion -T transcripts.fa --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.LRonly_incl_targets  --vis --extract_fusion_LR_fasta fusion_ev_LR.reads.fa --min_trans_overlap_length 75 --max_rigor --incl_fusion_targets fusion_supplement_list.txt

test_longOnly_fa_only_targets:
	../ctat-LR-fusion -T transcripts.fa --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.LRonly_only_targets  --vis --extract_fusion_LR_fasta fusion_ev_LR.reads.fa --min_trans_overlap_length 75 --max_rigor --only_fusion_targets fusion_supplement_list.txt


test_longOnly_fa_nophase2:
	../ctat-LR-fusion -T transcripts.fa --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short.nophase2  --min_num_LR 2  --vis --extract_fusion_LR_fasta fusion_ev_LR.reads.fa --min_trans_overlap_length 75

test_excl_shortreads_reg_mm2:
	../ctat-LR-fusion -T transcripts.fa --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short.reg_mm2 --vis --extract_fusion_LR_fasta fusion_ev_LR.reads.fa --min_trans_overlap_length 75 --no_ctat_mm2


test_excl_shortreads_bam_input:
	../ctat-LR-fusion --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short.LR_bam  --vis --extract_fusion_LR_fasta fusion_ev_LR.reads.fa --min_trans_overlap_length 75 --LR_bam transcripts.mm2.bam

test_excl_shortreads_set_num_total_reads:
	../ctat-LR-fusion -T transcripts.fa --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short.setNumReads --vis --extract_fusion_LR_fasta fusion_ev_LR.reads.fa --min_trans_overlap_length 75 --num_total_reads 500000

###########################################################
#################### Single cell ##########################
# synthetic single-cell inputs: the bulk test reads with fake CB:Z: and UB:Z: tags added,
# which ctat-LR-fusion should detect and encode into the read names as CB^UMI^readname

test_single_cell: test_sc_fq test_sc_LR_bam

transcripts.sc.fq.gz: transcripts.fq.gz
	./sc_test_data.py make_fastq transcripts.fq.gz transcripts.sc.fq.gz

transcripts.sc.mm2.bam: transcripts.mm2.bam
	./sc_test_data.py make_bam transcripts.mm2.bam transcripts.sc.mm2.bam

test_sc_fq: transcripts.sc.fq.gz
	../ctat-LR-fusion -T transcripts.sc.fq.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.sc_fq --min_novel_junction_support 1 --min_trans_overlap_length 75
	./sc_test_data.py check ctat_LR_fusion_outdir.sc_fq/ctat-LR-fusion.fusion_predictions.tsv
	../util/sc/fusion_cell_umi_summary.py --fusions ctat_LR_fusion_outdir.sc_fq/ctat-LR-fusion.fusion_predictions.tsv --output_prefix ctat_LR_fusion_outdir.sc_fq/ctat-LR-fusion.sc

test_sc_LR_bam: transcripts.sc.mm2.bam
	../ctat-LR-fusion --LR_bam transcripts.sc.mm2.bam --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.sc_LR_bam --min_trans_overlap_length 75
	./sc_test_data.py check ctat_LR_fusion_outdir.sc_LR_bam/ctat-LR-fusion.fusion_predictions.tsv


###########################################################
#################### With short reads too #################

test_incl_shortreads:
	../ctat-LR-fusion -T transcripts.fa --left_fq reads_1.fq.gz --right_fq reads_2.fq.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.incl_short --vis

test_incl_shortreads_SEonly:
	../ctat-LR-fusion -T transcripts.fa --left_fq reads_1.fq.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.incl_short_SEonly --vis



########### other ###############


test_incl_coding_effect: test_excl_shortreads
	../ctat-LR-fusion -T transcripts.fa --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short --vis --extract_fusion_LR_fasta fusion_ev_LR.reads.fa --examine_coding_effect


test_incl_FI_params:
	../ctat-LR-fusion -T transcripts.fa --left_fq reads_1.fq.gz --right_fq reads_2.fq.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.incl_short_custFI --vis  --FI_extra_params "--STAR_xtra_params '--limitBAMsortRAM 51419850732' " 


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
	rm -f transcripts.sc.fq.gz transcripts.sc.mm2.bam
	sudo rm -rf ./ctat_LR_fusion_outdir.*.docker





