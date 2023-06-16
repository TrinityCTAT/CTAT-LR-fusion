version 1.0


workflow minimap2_fusion_wf {

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
    
     call MINIMAP2_FUSION_TASK {
        input:
          sample_name=sample_name,
          transcripts=transcripts,
          genome_lib_tar=genome_lib_tar,
          min_per_id=min_per_id,
          docker=docker,
          cpu=cpu,
          memory=memory,
          preemptible=preemptible,
	      maxRetries=maxRetries,
          disk_space_multiplier=disk_space_multiplier    
     }
}


task MINIMAP2_FUSION_TASK {

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
    
    # minimap2-fusion

    minimap2-fusion --version

    minimap2-fusion -T ~{transcripts} \
                --genome_lib_dir ctat_genome_lib_build_dir \
                --min_J ~{min_J}  --min_sumJS ~{min_sumJS} --min_novel_junction_support ~{min_novel_junction_support} \
                --min_per_id ~{min_per_id} \
                -o minimap2_fusion_outdir


    mv minimap2_fusion_outdir/minimap2-fusion.fusion_predictions.tsv ~{sample_name}.minimap2-fusion.fusion_predictions.tsv 

    >>>
    
    output {
      File fusion_report="~{sample_name}.minimap2-fusion.fusion_predictions.tsv"
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

