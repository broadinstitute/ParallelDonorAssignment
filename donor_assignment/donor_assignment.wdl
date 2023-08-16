version 1.0

workflow donor_assign {
    input {
        File BAI
        File BAM
        Int num_splits
        File VCF
        File donor_list_file
        String docker_image = 'us.gcr.io/landerlab-atacseq-200218/donor_assign:0.19'
        String git_branch = "bgzblocks"
    }
    call generate_regions {
        input:
            BAI = BAI, 
            BAM_PATH = "~{BAM}",
            num_splits = num_splits,
            docker_image = docker_image,
            git_branch = git_branch
    }

    scatter (region in generate_regions.regions_list){
        call region_donor_log_likelihoods {
            input:
                BAI = BAI,
                BAM_PATH = "~{BAM}",
                region = region,
                VCF_PATH = "~{VCF}",
                donor_list_file = donor_list_file,
                docker_image = docker_image,
                git_branch = git_branch
        }
    }

    call gather_region_donor_log_likelihoods{
        input:
            barcode_log_likelihood = region_donor_log_likelihoods.barcode_log_likelihood,
            docker_image = docker_image,
            git_branch = git_branch
    }

    output {
        Array[String] regions_list = generate_regions.regions_list
        Array[File] region_barcode_log_likelihood_list = region_donor_log_likelihoods.barcode_log_likelihood 
        File total_barcode_donor_likelihoods = gather_region_donor_log_likelihoods.total_barcode_donor_likelihoods
    }
}

task generate_regions {
    input {
        File BAI
        String BAM_PATH
        Int num_splits
        String docker_image
        String git_branch
    }
    command {
        set -ex
        (git clone https://github.com/broadinstitute/ParallelDonorAssignment.git /app ; cd /app ; git checkout ${git_branch})
        gsutil cat ${BAM_PATH} | samtools view -H > header.sam
        python3 /app/donor_assignment/split_regions.py header.sam ${BAI} ${num_splits}
        head -5 list_of_regions.txt > only_five_regions.txt
    }
    output {
        Array[String] regions_list = read_lines("only_five_regions.txt")
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "16GB"
        preemptible: 1
    }
}

task region_donor_log_likelihoods {
    input {
        File BAI
        String BAM_PATH
        String region
        String VCF_PATH
        File donor_list_file
        String docker_image
        String git_branch
        String chrom_region = sub(region, " .*", "")
        String file_region = sub(region, ".* bytes:", "")
    }

    command {
        set -ex
        (git clone https://github.com/broadinstitute/ParallelDonorAssignment.git /app ; cd /app ; git checkout ${git_branch})
        sh /app/monitor_script.sh &

        # use gsutil instead of htslib for stability
        gsutil cat ${BAM_PATH} | samtools view -H -O bam > region.bam
        gsutil cat -r ${file_region} ${BAM_PATH} >> region.bam

        gsutil cp ${VCF_PATH} full.vcf.gz
        gsutil cp ${VCF_PATH}.tbi full.vcf.gz.tbi
        bcftools view -O z -o region.vcf.gz full.vcf.gz ${chrom_region}

        ls -l region.bam region.vcf.gz
        python3 /app/donor_assignment/count_reads_on_variants.py region.bam region.vcf.gz
        python3 /app/donor_assignment/likelihood_per_region.py results.tsv.gz ${donor_list_file} region.vcf.gz ${chrom_region}
    }

    output {
      File barcode_log_likelihood = glob("barcode_log_likelihood_*.txt.gz")[0]
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "32GB"
        preemptible: 1
    }
}


task gather_region_donor_log_likelihoods {
    input {
        Array[String] barcode_log_likelihood
        String docker_image
        String git_branch
    }

    command {
        (git clone https://github.com/broadinstitute/ParallelDonorAssignment.git /app ; cd /app ; git checkout ${git_branch})
        sh /app/monitor_script.sh &
        python3 /app/donor_assignment/gather_barcode_likelihoods.py ~{write_lines(barcode_log_likelihood)} total_barcode_donor_likelihoods.txt.gz
    }

    output {
        File total_barcode_donor_likelihoods = 'total_barcode_donor_likelihoods.txt.gz'
    }

    runtime {
        docker: docker_image
        cpu: 4
        memory: "32GB"
        preemptible: 1
    }
}
