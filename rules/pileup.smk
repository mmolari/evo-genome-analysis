rule map_reads:
    input:
        fa=f"{in_fld}/reads/{read_id}.fastq.gz",
        ref=f"{inf_fld}/references/{ref_id}.fa",
    output:
        sam=f"{out_fld}/mapped_reads/{ref_id}/{read_id}.sam",
        bam=f"{out_fld}/mapped_reads/{ref_id}/{read_id}.bam",
        bai=f"{out_fld}/mapped_reads/{ref_id}/{read_id}.bam.bai",
    conda:
        "conda_envs/pileup.yml"
    shell:
        """
        minimap2 -a -x map-ont {input.ref} {input.fa} > {output.sam}
        samtools sort {output.sam} > {output.bam}
        samtools index {output.bam}
        """
