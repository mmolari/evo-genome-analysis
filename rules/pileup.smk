pileup_config = config["pileup"]


rule subsample:
    message:
        "Subsampling {input.reads} to {params.target_bp} using filtlong (or skipping subsampling)."
    input:
        reads=in_fld + "/reads/{sample_id}.fastq.gz",
    output:
        subsampled=out_fld + "/subsampled_reads/{sample_id}.fastq.gz",
    params:
        target_bp=pileup_config.get("subsample_target_bp", None),
        min_length=2000,
    conda:
        "../conda_envs/pileup.yml"
    shell:
        """
        if [ -z "{params.target_bp}" ]; then
            # If no target_bp is set, just copy the original reads
            cp {input.reads} {output.subsampled}
        else
            # Run filtlong to subsample the reads
            filtlong \
                --min_length {params.min_length} \
                --target_bases {params.target_bp} \
                {input.reads} \
                | gzip > {output.subsampled}
        fi
        """


rule map_reads:
    input:
        fastq=lambda wildcards: (
            out_fld + "/subsampled_reads/{sample_id}.fastq.gz"
            if pileup_config.get("subsample_target_bp")
            else in_fld + "/reads/{sample_id}.fastq.gz"
        ),
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
        out_dir=lambda w: (
            out_fld + f"/pileup/{w.ref_id}/{w.rec_id}/{w.sample_id}"
        ).replace(" ", ""),
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


rule extract_counts:
    input:
        pileups=lambda w: expand(
            rules.build_pileup.output.pileup,
            sample_id=pileups[w.ref_id],
            allow_missing=True,
        ),
        ref=rules.map_reads.input.ref,
    output:
        cov=out_fld + "/pileup/{ref_id}/{rec_id}/coverage.npz",
        cons=out_fld + "/pileup/{ref_id}/{rec_id}/consensus.npz",
        gap=out_fld + "/pileup/{ref_id}/{rec_id}/gaps.npz",
    params:
        rec_id=lambda w: w.rec_id,
    conda:
        "../conda_envs/pileup.yml"
    shell:
        """
        python3 scripts/pileup/extract_counts.py \
            --pileups {input.pileups} \
            --ref_seq {input.ref} \
            --ref_record {params.rec_id} \
            --cov_ct {output.cov} \
            --cons_ct {output.cons} \
            --gap_ct {output.gap} \
        """


rule extract_insertion_counts:
    input:
        insertions=lambda w: expand(
            rules.build_pileup.output.insertions,
            sample_id=pileups[w.ref_id],
            allow_missing=True,
        ),
        ref=rules.map_reads.input.ref,
    output:
        insertions=out_fld + "/pileup/{ref_id}/{rec_id}/insertions.npz",
    params:
        rec_id=lambda w: w.rec_id,
    conda:
        "../conda_envs/pileup.yml"
    shell:
        """
        python3 scripts/pileup/extract_insertion_counts.py \
            --ref_fasta {input.ref} \
            --ref_record {params.rec_id} \
            --ins_dicts {input.insertions} \
            --ins_ct {output.insertions} \
        """


rule extract_clip_counts:
    input:
        clips=lambda w: expand(
            rules.build_pileup.output.clips,
            sample_id=pileups[w.ref_id],
            allow_missing=True,
        ),
        ref=rules.map_reads.input.ref,
    output:
        clips=out_fld + "/pileup/{ref_id}/{rec_id}/clips.npz",
    params:
        rec_id=lambda w: w.rec_id,
    conda:
        "../conda_envs/pileup.yml"
    shell:
        """
        python3 scripts/pileup/extract_clip_counts.py \
            --ref_fasta {input.ref} \
            --ref_record {params.rec_id} \
            --clip_dicts {input.clips} \
            --clip_ct {output.clips} \
        """


rule pileup_all:
    input:
        [
            expand(
                rules.extract_counts.output,
                ref_id=ref,
                rec_id=ref_records[ref],
            )
            for ref, reads in pileups.items()
        ],
        [
            expand(
                rules.extract_insertion_counts.output,
                ref_id=ref,
                rec_id=ref_records[ref],
            )
            for ref, reads in pileups.items()
        ],
        [
            expand(
                rules.extract_clip_counts.output,
                ref_id=ref,
                rec_id=ref_records[ref],
            )
            for ref, reads in pileups.items()
        ],
        [
            expand(rules.extract_unmapped.output, ref_id=ref, sample_id=reads)
            for ref, reads in pileups.items()
        ],
        [
            expand(rules.extract_nonprimary.output, ref_id=ref, sample_id=reads)
            for ref, reads in pileups.items()
        ],
        [
            expand(rules.map_summary.output, ref_id=ref, sample_id=reads)
            for ref, reads in pileups.items()
        ],
