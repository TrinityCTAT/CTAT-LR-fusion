
simple: test_excl_shortreads

both:  test_incl_shortreads test_excl_shortreads

all: test_incl_shortreads test_excl_shortreads test_longOnly_fqgz


test_incl_shortreads:
	../ctat-LR-fusion -T transcripts.fa --left_fq reads_1.fq.gz --right_fq reads_2.fq.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.incl_short --vis


test_excl_shortreads:
	../ctat-LR-fusion -T transcripts.fa --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short --min_J 1  --min_sumJS 1 --vis


test_longOnly_fqgz:
	../ctat-LR-fusion -T transcripts.fq.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o ctat_LR_fusion_outdir.excl_short.fqgz --min_J 1  --min_sumJS 1 --min_novel_junction_support 1

test_docker:
	../Docker/run_docker.sh ctat-LR-fusion -T /data/transcripts.fa --genome_lib_dir /ctat_genome_lib --min_J 1  --min_sumJS 1 --left_fq reads_1.fq.gz --right_fq reads_2.fq.gz -o ctat_LR_fusion_outdir.docker --vis


test_singularity:
	singularity exec -e -B `pwd`:/data -B ${CTAT_GENOME_LIB}:/ctat_genome_lib ../Docker/ctat_lr_fusion.simg ctat-LR-fusion -T /data/transcripts.fa --genome_lib_dir /ctat_genome_lib --min_J 1  --min_sumJS 1 --left_fq reads_1.fq.gz --right_fq reads_2.fq.gz -o ctat_LR_fusion_outdir.singularity --vis

clean:
	rm -rf ./ctat_LR_fusion_outdir.incl_short ./ctat_LR_fusion_outdir.excl_short*
	sudo rm -rf ./ctat_LR_fusion_outdir.docker





