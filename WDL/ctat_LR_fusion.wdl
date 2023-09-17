version 1.0


workflow ctat_LR_fusion_wf {

    input {
       String sample_name
       File transcripts
       File genome_lib_tar_mm2_only
       File genome_lib_tar_with_STAR_idx
       Int min_per_id=90
       Int min_J = 1
       Int min_sumJS = 1    
       Int min_novel_junction_support = 1
       File? illumina_left_fq
       File? illumina_right_fq
      
       String docker="trinityctat/ctat_lr_fusion:latest"
       Int cpu = 10
       String memory="50G"
       Int preemptible = 0
       Int maxRetries = 0
       Float disk_space_multiplier = 3.0

      
     }
    
     call CTAT_LR_FUSION_TASK {
        input:
          sample_name=sample_name,
          transcripts=transcripts,
          genome_lib_tar= if defined(illumina_left_fq) then genome_lib_tar_with_STAR_idx else genome_lib_tar_mm2_only,
          min_per_id=min_per_id,
          min_J=min_J,
          min_sumJS=min_sumJS,    
          min_novel_junction_support=min_novel_junction_support,
          illumina_left_fq=illumina_left_fq,
	      illumina_right_fq=illumina_right_fq,

          docker=docker,
          cpu=cpu,
          memory=memory,
          preemptible=preemptible,
          maxRetries=maxRetries,
          disk_space_multiplier=disk_space_multiplier    
     }
}


task CTAT_LR_FUSION_TASK {

    input {
       String sample_name
       File transcripts
       File genome_lib_tar
       Int min_per_id
       Int min_J
       Int min_sumJS    
       Int min_novel_junction_support
       File? illumina_left_fq
       File? illumina_right_fq

       String docker
       Int cpu
       String memory
       Int preemptible
       Int maxRetries
       Float disk_space_multiplier
  }

  Int disk_space = ceil(size(genome_lib_tar, "GB") * disk_space_multiplier)
  
  command <<<

    set -ex

    # untar the genome lib
    tar xvf ~{genome_lib_tar}
    rm ~{genome_lib_tar}
    
    # ctat-LR-fusion

    ctat-LR-fusion --version

    ctat-LR-fusion -T ~{transcripts} \
                --genome_lib_dir ctat_genome_lib_build_dir \
                --min_J ~{min_J}  --min_sumJS ~{min_sumJS} --min_novel_junction_support ~{min_novel_junction_support} \
                --min_per_id ~{min_per_id} \
                --CPU ~{cpu} \
                --vis \
                ~{"--left_fq " + illumina_left_fq} ~{"--right_fq " + illumina_right_fq } \
                -o ctat_LR_fusion_outdir


    mv ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.preliminary.tsv ~{sample_name}.ctat-LR-fusion.fusion_predictions.preliminary.tsv
    mv ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.preliminary.abridged.tsv ~{sample_name}.ctat-LR-fusion.fusion_predictions.preliminary.abridged.tsv
  
    mv ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.tsv ~{sample_name}.ctat-LR-fusion.fusion_predictions.tsv
    mv ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.abridged.tsv ~{sample_name}.ctat-LR-fusion.fusion_predictions.abridged.tsv 

    mv ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_inspector_web.html ~{sample_name}.ctat-LR-fusion.fusion_inspector_web.html

    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.genome.fa ~{sample_name}.ctat-LR-fusion.igv.genome.fa
    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.genome.fa.fai ~{sample_name}.ctat-LR-fusion.igv.genome.fa.fai
    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.annot.gtf ~{sample_name}.ctat-LR-fusion.igv.annot.gtf
    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.annot.bed ~{sample_name}.ctat-LR-fusion.igv.annot.bed
    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.LR.sorted.bam ~{sample_name}.ctat-LR-fusion.igv.LR.sorted.bam
    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.LR.sorted.bam.bai ~{sample_name}.ctat-LR-fusion.igv.LR.sorted.bam.bai
    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.pfam.bed ~{sample_name}.ctat-LR-fusion.igv.pfam.bed
    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.seqsimilar.bed ~{sample_name}.ctat-LR-fusion.igv.seqsimilar.bed
    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.LR.breakoint.roi.bed ~{sample_name}.ctat-LR-fusion.igv.LR.breakoint.roi.bed

    tar -zcvf ~{sample_name}.ctat-LR-fusion.igv.tar.gz ~{sample_name}.ctat-LR-fusion.igv.*
    

    
    >>>
    
    output {
      File fusion_report="~{sample_name}.ctat-LR-fusion.fusion_predictions.tsv"
      File fusion_report_abridged="~{sample_name}.ctat-LR-fusion.fusion_predictions.abridged.tsv"

      File prelim_fusion_report="~{sample_name}.ctat-LR-fusion.fusion_predictions.preliminary.tsv"
      File prelim_fusion_report_abridged="~{sample_name}.ctat-LR-fusion.fusion_predictions.preliminary.abridged.tsv"

      File fusion_report_html="~{sample_name}.ctat-LR-fusion.fusion_inspector_web.html"
      File igv_tar="~{sample_name}.ctat-LR-fusion.igv.tar.gz"
    }
    

    runtime {
            docker: "~{docker}"
            disks: "local-disk " + disk_space + " HDD"
            memory: "~{memory}"
            cpu: cpu
            preemptible: preemptible
            maxRetries: maxRetries
    }
}

