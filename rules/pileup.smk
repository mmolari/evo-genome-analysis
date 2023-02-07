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
