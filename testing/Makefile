
simple: test_excl_shortreads

both:  test_incl_shortreads test_excl_shortreads

all: test_incl_shortreads test_excl_shortreads test_longOnly_fqgz test_longOnly_fq


test_excl_shortreads:
	../ctat-LR-fusion -T transcripts.fa --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short --min_J 1  --min_sumJS 1 --vis --extract_fusion_LR_fasta fusion_ev_LR.reads.fa



test_incl_shortreads:
	../ctat-LR-fusion -T transcripts.fa --left_fq reads_1.fq.gz --right_fq reads_2.fq.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.incl_short --vis

test_incl_shortreads_SEonly:
	../ctat-LR-fusion -T transcripts.fa --left_fq reads_1.fq.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.incl_short_SEonly --vis



test_longOnly_fqgz:
	../ctat-LR-fusion -T transcripts.fq.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short.fqgz --min_J 1  --min_sumJS 1 --min_novel_junction_support 1

test_longOnly_fq:
	../ctat-LR-fusion -T transcripts.fq --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short.fq --min_J 1  --min_sumJS 1 --min_novel_junction_support 1


test_oncoreads:
	../ctat-LR-fusion -T oncofusionreads.fasta.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.oncoreads --vis

test_oncoreads_chim_candidates_only:
	../ctat-LR-fusion -T oncofusionreads.fasta.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.oncoreads.chim_candidates_only --chim_candidates_only


test_docker:
	../Docker/run_docker.sh ctat-LR-fusion -T /data/transcripts.fa --genome_lib_dir /ctat_genome_lib --min_J 1  --min_sumJS 1 --left_fq reads_1.fq.gz --right_fq reads_2.fq.gz -o ctat_LR_fusion_outdir.docker --vis


test_singularity:
	singularity exec -e -B `pwd`:/data -B ${CTAT_GENOME_LIB}:/ctat_genome_lib ../Docker/ctat_lr_fusion.simg ctat-LR-fusion -T /data/transcripts.fa --genome_lib_dir /ctat_genome_lib --min_J 1  --min_sumJS 1 --left_fq reads_1.fq.gz --right_fq reads_2.fq.gz -o ctat_LR_fusion_outdir.singularity --vis

clean:
	rm -rf ./ctat_LR_fusion_outdir.incl_short ./ctat_LR_fusion_outdir.excl_short* ./ctat_LR_fusion_outdir.singularity ./ctat_LR_fusion_outdir.oncoreads
	sudo rm -rf ./ctat_LR_fusion_outdir.docker





