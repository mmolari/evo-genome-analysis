rule coverage:
    input:
        cov=rules.extract_counts.output.cov,
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
            --coverage {input.cov} \
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
