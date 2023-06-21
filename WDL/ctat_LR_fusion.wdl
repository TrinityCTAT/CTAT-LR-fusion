version 1.0


workflow ctat_LR_fusion_wf {

    input {
       String sample_name
       File transcripts
       File genome_lib_tar
       Int min_per_id=90
       Int min_J = 1
       Int min_sumJS = 1    
       Int min_novel_junction_support = 1
       
      
       String docker="trinityctat/minimap2fusion:latest"
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
          genome_lib_tar=genome_lib_tar,
          min_per_id=min_per_id,
          min_J=min_J,
          min_sumJS=min_sumJS,    
          min_novel_junction_support=min_novel_junction_support,

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
                -o ctat_LR_fusion_outdir


    mv ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.preliminary.tsv ~{sample_name}.ctat-LR-fusion.fusion_predictions.preliminary.tsv
  
    mv ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.tsv ~{sample_name}.ctat-LR-fusion.fusion_predictions.tsv 
    


    >>>
    
    output {
      File fusion_report="~{sample_name}.ctat-LR-fusion.fusion_predictions.tsv"
      File prelim_fusion_report="~{sample_name}.ctat-LR-fusion.fusion_predictions.preliminary.tsv"
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

