version 1.0

workflow cisvar {
    input {
        File BAM
        File VCF
    }

    call est_cisvar {
        input:
        BAM="~{BAM}",
        VCF="~{VCF}",
        bam_size=size(BAM, "GB"),
        vcf_size=size(VCF, "GB")
    }
    output {
        File donor_weights = est_cisvar.donor_weights
        File snp_estimates = est_cisvar.snp_estimates
    }
}

task est_cisvar {
    input {
        String BAM
        String VCF
        File donor_list
        Int min_coverage = 20
        Float bam_size
        Float vcf_size
    }

    Int disk_size = ceil(bam_size * 2.5)

    command {
        set -ex
        ## Work from google buckets directly for faster startup
        ## from https://support.terra.bio/hc/en-us/community/posts/16214505476507-How-to-run-samtools-on-gs-object-directly-to-get-a-BAM-slice-fast-for-example-
        gcloud auth print-access-token > token.txt
        export HTS_AUTH_LOCATION="token.txt"
        samtools view -d vW:1 ${BAM} -O BAM | samtools sort - -o sorted.bam
        # check file size
        SIZE=$(du -sk sorted.bam | sed "s/[^0-9].*//")
        if ((SIZE<1000)) ; then
            >&2 echo "Bam file shorter than 1meg - did you run STAR in WASP mode?"
            exit 1;
        fi
        pip3 install scipy --break-system-packages
        python3 /app/cisvar/cisvar.py ${donor_list} ${VCF} sorted.bam ${min_coverage}
    }
    output {
        File donor_weights = "donor_weights.txt"
        File snp_estimates = "cisvar_estimates.txt.gz"
    }
    runtime {
        docker: "us.gcr.io/landerlab-atacseq-200218/donor_assign:0.7"
        cpu: 4
        memory: "64GB"
        disks: "local_disk ~{disk_size}GB HDD"
        preemptible: 1
    }
}
