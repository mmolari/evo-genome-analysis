rule concat_positions:
    input:
        gaps=rules.plot_gaps.output,
        cons=rules.plot_consensus.output,
        ins=rules.plot_insertions.output,
        clips=rules.plot_clips.output,
    output:
        pos=out_fld + "/annotations/{ref_id}/{rec_id}/positions.csv",
    conda:
        "../conda_envs/plots.yml"
    shell:
        """
        python scripts/annotations/concat_positions.py \
            --gaps {input.gaps} \
            --cons {input.cons} \
            --ins {input.ins} \
            --clips {input.clips} \
            --out {output.pos}
        """


rule assing_annotations:
    input:
        ref=in_fld + "/references/{ref_id}.gbk",
        pos=rules.concat_positions.output,
    output:
        hit=out_fld + "/annotations/{ref_id}/{rec_id}/annotations_hit.csv",
        miss=out_fld + "/annotations/{ref_id}/{rec_id}/annotations_miss.csv",
    conda:
        "../conda_envs/plots.yml"
    shell:
        """
        python scripts/annotations/assign_annotations.py \
            --ref {input.ref} \
            --record {wildcards.rec_id} \
            --pos {input.pos} \
            --out_hit {output.hit} \
            --out_miss {output.miss}
        """


rule annotate_all:
    input:
        [
            expand(
                rules.assing_annotations.output,
                ref_id=ref,
                rec_id=ref_records[ref],
            )
            for ref, reads in pileups.items()
        ],
