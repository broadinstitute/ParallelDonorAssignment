version 1.0

workflow donor_assign {
    call generate_regions
    output {
        Array[String] regions_list = generate_regions.regions_list
    }
}

task generate_regions {
    input {
        File BAI
        String BAM_PATH
        Int num_splits
        String docker_image = 'http://us.gcr.io/landerlab-atacseq-200218/donor_assign:latest'
    }
    command {
        gsutil cat ${BAM_PATH} | samtools view -H > header.sam
        python3 split_regions.py header.sam ${BAI} ${num_splits} # for now, not going to pipe because we have other prints in there ## > list_of_regions.txt

    }
    output {
        Array[String] regions_list = read_lines("list_of_regions.txt")
    }

    runtime {
        docker: docker_image
        cpu: 4
        memory: "32GB"
    }

    # scatter(r in regions){
    #     call extract_bam_region{
    #         input{
    #             BAM
    #             BAI
    #             r
    #         }
    #     }
    # }
}