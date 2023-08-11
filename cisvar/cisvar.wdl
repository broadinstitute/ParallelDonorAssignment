version 1.0

workflow cisvar {
    call est_cisvar
    output {
        File donor_weights = est_cisvar.donor_weights
        File snp_estimates = est_cisvar.snp_estimates
    }
}

task est_cisvar {
    input {
        File BAM
        File VCF
        File donor_list
        Int min_coverage = 20
    }
    command {
        samtools sort ${BAM} -o sorted.bam
        python3 /app/cisvar/cisvar.py ${donor_list} ${VCF} ${BAM} ${min_coverage}
    }
    output {
        File donor_weights = "donor_weights.txt"
        File snp_estimates = "cisvar_estimates.txt.gz"
    }
    runtime {
        docker: "us.gcr.io/landerlab-atacseq-200218/donor_assign:0.6"
        cpu: 4
        memory: "32GB"
    }
}
