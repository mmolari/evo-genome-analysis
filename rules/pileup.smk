pileup_config = config["pileup"]


rule map_reads:
    input:
        fastq=in_fld + "/reads/{sample_id}.fastq.gz",
        ref=in_fld + "/references/{ref_id}.fa",
    output:
        sam=out_fld + "/mapped_reads/{ref_id}/{sample_id}.sam",
        bam=out_fld + "/mapped_reads/{ref_id}/{sample_id}.bam",
        bai=out_fld + "/mapped_reads/{ref_id}/{sample_id}.bam.bai",
    conda:
        "../conda_envs/pileup.yml"
    shell:
        """
        minimap2 -a -x map-ont {input.ref} {input.fastq} > {output.sam}
        samtools sort {output.sam} > {output.bam}
        samtools index {output.bam}
        """


rule map_summary:
    input:
        bam=rules.map_reads.output.bam,
    output:
        txt=out_fld + "/mapped_reads/{ref_id}/{sample_id}.info.txt",
    conda:
        "../conda_envs/pileup.yml"
    shell:
        """
        python3 scripts/pileup/map_summary.py \
            --bam {input.bam} > {output.txt}
        """


rule extract_nonprimary:
    input:
        bam=rules.map_reads.output.bam,
    output:
        csv=out_fld + "/pileup/{ref_id}/non_primary/{sample_id}/non_primary.csv",
    conda:
        "../conda_envs/pileup.yml"
    shell:
        """
        python3 scripts/pileup/extract_nonprimary.py \
            --bam {input.bam} --csv {output.csv}
        """


rule extract_unmapped:
    input:
        bam=rules.map_reads.output.bam,
    output:
        csv=out_fld + "/pileup/{ref_id}/non_primary/{sample_id}/unmapped.csv",
        fastq=out_fld + "/pileup/{ref_id}/non_primary/{sample_id}/unmapped.fastq.gz",
    conda:
        "../conda_envs/pileup.yml"
    shell:
        """
        python3 scripts/pileup/extract_unmapped.py \
            --bam {input.bam} --csv {output.csv} --fastq {output.fastq}
        """


rule build_pileup:
    input:
        bam=rules.map_reads.output.bam,
    output:
        pileup=out_fld + "/pileup/{ref_id}/{rec_id}/{sample_id}/allele_counts.npz",
        insertions=out_fld + "/pileup/{ref_id}/{rec_id}/{sample_id}/insertions.pkl.gz",
        clips=out_fld + "/pileup/{ref_id}/{rec_id}/{sample_id}/clips.pkl.gz",
    params:
        min_q=pileup_config["qual_min"],
        min_L=pileup_config["clip_minL"],
        out_dir=lambda w: out_fld + f"/pileup/{w.ref_id}/{w.rec_id}/{w.sample_id}",
    conda:
        "../conda_envs/pileup.yml"
    shell:
        """
        python3 scripts/pileup/build_pileup.py \
            --bam_file {input.bam} \
            --out_dir {params.out_dir} \
            --ref_record_name {wildcards.rec_id} \
            --qual_min {params.min_q} \
            --clip_minL {params.min_L} \
        """


rule extract_cons_gap_freqs:
    input:
        pileup=rules.build_pileup.output.pileup,
        ref=rules.map_reads.input.ref,
    output:
        npz=out_fld + "/pileup/{ref_id}/{rec_id}/{sample_id}/cons_gap_freqs.npz",
    conda:
        "../conda_envs/pileup.yml"
    shell:
        """
        python3 scripts/pileup/extract_freqs.py \
            --ref_genome {input.ref} --pileup {input.pileup} \
            --record {wildcards.rec_id} --npz {output.npz} \
            --verbose
        """
