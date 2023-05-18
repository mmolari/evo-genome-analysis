plot_config = config["plots"]


rule plot_coverage:
    input:
        cov=rules.extract_counts.output.cov,
    output:
        fld=directory(out_fld + "/figs/{ref_id}/{rec_id}/coverage"),
    params:
        maxbins=plot_config["coverage-maxbins"],
    conda:
        "../conda_envs/plots.yml"
    shell:
        """
        python3 scripts/plots/coverage.py \
            --coverage {input.cov} \
            --maxbins {params.maxbins} \
            --out_fld {output.fld} \
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


rule plot_consensus:
    input:
        cons=rules.extract_counts.output.cons,
        cov=rules.extract_counts.output.cov,
    output:
        fld=directory(out_fld + "/figs/{ref_id}/{rec_id}/consensus"),
    params:
        freq_thr=plot_config["consensus"]["freq-threshold"],
        cov_thr=plot_config["consensus"]["coverage-threshold"],
        n_top_trajs=plot_config["consensus"]["n-top-trajs"],
    conda:
        "../conda_envs/plots.yml"
    shell:
        """
        python3 scripts/plots/consensus.py \
            --freq_thr {params.freq_thr} \
            --cov_thr {params.cov_thr} \
            --n_top_trajs {params.n_top_trajs} \
            --cons_npz {input.cons} \
            --cov_npz {input.cov} \
            --plot_fld {output.fld} \
        """


rule plot_insertions:
    input:
        ins=rules.extract_insertion_counts.output.insertions,
        cov=rules.extract_counts.output.cov,
    output:
        fld=directory(out_fld + "/figs/{ref_id}/{rec_id}/insertions"),
    params:
        freq_thr=plot_config["insertions"]["freq-threshold"],
        cov_thr=plot_config["insertions"]["coverage-threshold"],
        n_top_trajs=plot_config["insertions"]["n-top-trajs"],
        cov_window=plot_config["insertions"]["coverage-window"],
        cov_fraction=plot_config["insertions"]["coverage-fraction"],
    conda:
        "../conda_envs/plots.yml"
    shell:
        """
        python3 scripts/plots/insertions.py \
            --freq_thr {params.freq_thr} \
            --cov_thr {params.cov_thr} \
            --n_top_trajs {params.n_top_trajs} \
            --ins_npz {input.ins} \
            --cov_npz {input.cov} \
            --plot_fld {output.fld} \
            --cov_window {params.cov_window} \
            --cov_fraction {params.cov_fraction} \
        """


rule plot_clips:
    input:
        clips=rules.extract_clip_counts.output.clips,
        cov=rules.extract_counts.output.cov,
    output:
        fld=directory(out_fld + "/figs/{ref_id}/{rec_id}/clips"),
    params:
        freq_thr=plot_config["clips"]["freq-threshold"],
        cov_thr=plot_config["clips"]["coverage-threshold"],
        n_top_trajs=plot_config["clips"]["n-top-trajs"],
        cov_window=plot_config["clips"]["coverage-window"],
        cov_fraction=plot_config["clips"]["coverage-fraction"],
    conda:
        "../conda_envs/plots.yml"
    shell:
        """
        python3 scripts/plots/clips.py \
            --freq_thr {params.freq_thr} \
            --cov_thr {params.cov_thr} \
            --n_top_trajs {params.n_top_trajs} \
            --clips_npz {input.clips} \
            --cov_npz {input.cov} \
            --plot_fld {output.fld} \
            --cov_window {params.cov_window} \
            --cov_fraction {params.cov_fraction} \
        """


rule plot_supplementary:
    input:
        csv=rules.extract_nonprimary.output.csv,
    output:
        fld=directory(out_fld + "/figs/{ref_id}/non_primary/{sample_id}/supplementary"),
    params:
        L_thr=plot_config["supplementary"]["length-threshold"],
    conda:
        "../conda_envs/plots.yml"
    shell:
        """
        python3 scripts/plots/supplementary.py \
            --L_thr {params.L_thr} \
            --in_csv {input.csv} \
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
        [
            expand(
                rules.plot_consensus.output,
                ref_id=ref,
                rec_id=ref_records[ref],
            )
            for ref, reads in pileups.items()
        ],
        [
            expand(
                rules.plot_insertions.output,
                ref_id=ref,
                rec_id=ref_records[ref],
            )
            for ref, reads in pileups.items()
        ],
        [
            expand(
                rules.plot_clips.output,
                ref_id=ref,
                rec_id=ref_records[ref],
            )
            for ref, reads in pileups.items()
        ],
        [
            expand(rules.plot_supplementary.output, ref_id=ref, sample_id=reads)
            for ref, reads in pileups.items()
        ],
