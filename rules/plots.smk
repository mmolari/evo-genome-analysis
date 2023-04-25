rule coverage:
    input:
        pileups=lambda w: expand(
            rules.build_pileup.output.pileup,
            sample_id=pileups[w.ref_id],
            allow_missing=True,
        ),
    output:
        png=out_fld + "/figs/{ref_id}/{rec_id}/coverage.png",
        html=out_fld + "/figs/{ref_id}/{rec_id}/coverage.html",
    params:
        maxbins=config["plots"]["coverage-maxbins"],
    conda:
        "../conda_envs/plots.yml"
    shell:
        """
        python3 scripts/plots/coverage.py \
            --samples {input.pileups} \
            --maxbins {params.maxbins} \
            --out_png {output.png} \
            --out_html {output.html}
        """


rule plot_all:
    input:
        [
            expand(
                rules.coverage.output,
                ref_id=ref,
                rec_id=ref_records[ref],
            )
            for ref, reads in pileups.items()
        ],
