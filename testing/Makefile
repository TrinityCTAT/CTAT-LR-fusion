all: test_incl_shortreads test_excl_shortreads test_longOnly_fqgz


test_incl_shortreads:
	../GMAP-fusion -T transcripts.fa --left_fq reads_1.fq.gz --right_fq reads_2.fq.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o gmap_fusion_outdir.incl_short


test_excl_shortreads:
	../GMAP-fusion -T transcripts.fa --genome_lib_dir ${CTAT_GENOME_LIB} -o gmap_fusion_outdir.excl_short --min_J 1  --min_sumJS 1


test_longOnly_fqgz:
	../GMAP-fusion -T transcripts.fq.gz --genome_lib_dir ${CTAT_GENOME_LIB} -o gmap_fusion_outdir.excl_short.fqgz --min_J 1  --min_sumJS 1

clean:
	rm -rf ./gmap_fusion_outdir.incl_short ./gmap_fusion_outdir.excl_short*




