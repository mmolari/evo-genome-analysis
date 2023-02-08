rule map_reads:
    input:
        fa=in_fld + "/reads/{read_id}.fastq.gz",
        ref=in_fld + "/references/{ref_id}.fa",
    output:
        sam=out_fld + "/mapped_reads/{ref_id}/{read_id}.sam",
        bam=out_fld + "/mapped_reads/{ref_id}/{read_id}.bam",
        bai=out_fld + "/mapped_reads/{ref_id}/{read_id}.bam.bai",
    conda:
        "../conda_envs/pileup.yml"
    shell:
        """
        minimap2 -a -x map-ont {input.ref} {input.fa} > {output.sam}
        samtools sort {output.sam} > {output.bam}
        samtools index {output.bam}
        """


# rule demultiplex_mappings_stats:
# txt/json file with n. / len of primary mappings stratified by reference
# and
#     input:
#     output:
#     conda:
#     shell:

# rule extract_nonprimary_mappings:
# csv with secondary/supplementary mappings info, with all the corresponding references.
#     input:
#     output:
#     conda:
#     shell:


rule build_pileup:
    input:
        bam=rules.map_reads.output.bam,
    output:
        pileup=out_fld + "/pileup/{ref_id}/{rec_id}/{read_id}/allele_counts.npz",
        insertions=out_fld + "/pileup/{ref_id}/{rec_id}/{read_id}/insertions.pkl.gz",
        clips=out_fld + "/pileup/{ref_id}/{rec_id}/{read_id}/clips.pkl.gz",
    params:
        min_q=config["qual_min"],
        min_L=config["clip_minL"],
        out_dir=lambda w: out_fld + f"/pileup/{w.ref_id}/{w.rec_id}/{w.read_id}",
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
