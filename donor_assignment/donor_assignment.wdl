version 1.0

workflow donor_assign {
    input {
        File BAI
        File BAM
        Int num_splits
    }
    call generate_regions {
        input:
            BAI = BAI, 
            BAM_PATH = "~{BAM}",
            num_splits = num_splits,
    }

    scatter (region in generate_regions.regions_list){
        call count_region {
            input:
                BAI = BAI,
                BAM_PATH = "~{BAM}",
                region = region
        }
    }

    output {
        Array[String] regions_list = generate_regions.regions_list
    }
}

task generate_regions {
    input {
        File BAI
        String BAM_PATH
        Int num_splits
        String docker_image = 'us.gcr.io/landerlab-atacseq-200218/donor_assign:0.3'
    }
    command {
        gsutil cat ${BAM_PATH} | samtools view -H > header.sam
        python3 /app/donor_assignment/split_regions.py header.sam ${BAI} ${num_splits} # for now, not going to pipe because we have other prints in there ## > list_of_regions.txt
        head -5 list_of_regions.txt > only_five_regions.txt
    }
    output {
        Array[String] regions_list = read_lines("only_five_regions.txt")
    }

    runtime {
        docker: docker_image
        cpu: 4
        memory: "32GB"
    }
}

task count_region {
    input {
        File BAI
        String BAM_PATH
        String region
        String docker_image = 'us.gcr.io/landerlab-atacseq-200218/donor_assign:0.3'
    }

    command {
        ## from https://support.terra.bio/hc/en-us/community/posts/16214505476507-How-to-run-samtools-on-gs-object-directly-to-get-a-BAM-slice-fast-for-example-
        ## write the GCP token in a file
        gcloud auth print-access-token > token.txt
        ## point the HTS ENV variable to that file
        export HTS_AUTH_LOCATION="token.txt"
        samtools view -X ${BAI} -b ${BAM_PATH} -o region.bam region
        ls -l region.bam
        echo $GOOGLE_APPLICATION_CREDENTIALS
    }

    runtime {
        docker: docker_image
        cpu: 4
        memory: "32GB"
    }
}