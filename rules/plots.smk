plot_config = config["plots"]


rule plot_coverage:
    input:
        cov=rules.extract_counts.output.cov,
    output:
        png=out_fld + "/figs/{ref_id}/{rec_id}/coverage.png",
        html=out_fld + "/figs/{ref_id}/{rec_id}/coverage.html",
    params:
        maxbins=plot_config["coverage-maxbins"],
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


rule plot_gaps:
    input:
        gap=rules.extract_counts.output.gap,
        cov=rules.extract_counts.output.cov,
    output:
        fld=directory(out_fld + "/figs/{ref_id}/{rec_id}/gaps"),
    params:
        freq_thr=plot_config["gaps"]["freq-threshold"],
        cov_thr=plot_config["gaps"]["coverage-threshold"],
        n_top_trajs=plot_config["gaps"]["n-top-trajs"],
    conda:
        "../conda_envs/plots.yml"
    shell:
        """
        python3 scripts/plots/gaps.py \
            --freq_thr {params.freq_thr} \
            --cov_thr {params.cov_thr} \
            --n_top_trajs {params.n_top_trajs} \
            --gap_npz {input.gap} \
            --cov_npz {input.cov} \
            --plot_fld {output.fld} \
        """


rule plot_all:
    input:
        [
            expand(
                rules.plot_coverage.output,
                ref_id=ref,
                rec_id=ref_records[ref],
            )
            for ref, reads in pileups.items()
        ],
        [
            expand(
                rules.plot_gaps.output,
                ref_id=ref,
                rec_id=ref_records[ref],
            )
            for ref, reads in pileups.items()
        ],
