version 1.0

workflow donor_assign {
    input {
        File BAI
        File BAM
        Int num_splits
        File VCF
        File donor_list_file
        String likelihood_method
        String docker_image = 'us.gcr.io/landerlab-atacseq-200218/donor_assign:0.20'
        String git_branch = "dropulation_likelihoods"
    }

    Int bam_split_size = floor(size(BAM, "GB") / num_splits)
    Int vcf_size = floor(size(BAM, "GB"))

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
                bam_size = bam_split_size,
                vcf_size = vcf_size,
                likelihood_method = likelihood_method,
                docker_image = docker_image,
                git_branch = git_branch
        }
    }

    call gather_region_donor_log_likelihoods{
        input:
            barcode_log_likelihood = region_donor_log_likelihoods.barcode_log_likelihood,
            files_size = size(region_donor_log_likelihoods.barcode_log_likelihood, "GB"),
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
        python3 -u /app/donor_assignment/split_regions.py header.sam ${BAI} ${num_splits}
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
        Int bam_size
        Int vcf_size
        String likelihood_method
        String docker_image
        String git_branch
    }

    String chrom_region = sub(region, " .*", "")
    String file_region = sub(region, ".* bytes:", "")
    Int disk_size = bam_size * 2 + vcf_size + 100

    command {
        set -ex
        (git clone https://github.com/broadinstitute/ParallelDonorAssignment.git /app ; cd /app ; git checkout ${git_branch})
        bash /app/monitor_script.sh &

        # use gsutil instead of htslib for stability
        gsutil cat ${BAM_PATH} | samtools view -H -O bam > region.bam
        gsutil cat -r ${file_region} ${BAM_PATH} >> region.bam

        gsutil -q cp ${VCF_PATH} full.vcf.gz
        gsutil -q cp ${VCF_PATH}.tbi full.vcf.gz.tbi
        bcftools view -O z -o region.vcf.gz full.vcf.gz ${chrom_region}

        ls -l region.bam region.vcf.gz
        python3 -u /app/donor_assignment/count_reads_on_variants.py region.bam region.vcf.gz
        python3 -u /app/donor_assignment/likelihood_per_region.py results.tsv.gz ${donor_list_file} region.vcf.gz ${chrom_region} ${likelihood_method}
    }

    output {
      File barcode_log_likelihood = glob("barcode_log_likelihood_*.txt.gz")[0]
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "32GB"
        preemptible: 2
        disks: "local-disk ~{disk_size} HDD"
    }
}


task gather_region_donor_log_likelihoods {
    input {
        Array[String] barcode_log_likelihood
        Float files_size
        String docker_image
        String git_branch
    }

    Int disk_size = ceil(files_size)

    command {
        (git clone https://github.com/broadinstitute/ParallelDonorAssignment.git /app ; cd /app ; git checkout ${git_branch})
        bash /app/monitor_script.sh &
        python3 -u /app/donor_assignment/gather_barcode_likelihoods.py ~{write_lines(barcode_log_likelihood)} total_barcode_donor_likelihoods.txt.gz
    }

    output {
        File total_barcode_donor_likelihoods = 'total_barcode_donor_likelihoods.txt.gz'
    }

    runtime {
        docker: docker_image
        cpu: 4
        memory: "32GB"
        preemptible: 1
        disks: "local-disk ~{disk_size} HDD"
    }
}
