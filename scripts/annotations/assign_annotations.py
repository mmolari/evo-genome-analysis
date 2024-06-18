import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO, SeqFeature
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", help="genbank reference file")
    parser.add_argument("--record", help="genabnk record name")
    parser.add_argument("--pos", help="position csv table")
    parser.add_argument("--out_hit", help="output csv annotation table for hits")
    parser.add_argument("--out_miss", help="output csv annotation table for missed")
    return parser.parse_args()


def parse_genbank(fname, record_name, pos_list):
    # turn to 0-based positions
    positions = [SeqFeature.ExactPosition(p - 1) for p in pos_list]
    results = defaultdict(list)
    seq = None
    record_found = False
    found_records = []
    with open(fname) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            found_records.append(record.name)
            if record.name != record_name:
                continue
            record_found = True
            for feature in record.features:
                if feature.type == "source":
                    seq = feature.extract(record.seq)
                    continue
                for pos in positions:
                    if int(pos) in feature:
                        results[pos].append(feature)

    assert record_found, f"Record {record_name} not found in {fname} -> {found_records}"
    return results, seq


def parse_cds(feature, seq, pos):
    """
    extract:
    - aminoacid
    - aminoacid position
    - codon
    - position in codon
    """
    w = "none"
    # 0-based
    assert type(pos) == SeqFeature.ExactPosition
    # only simple location implemented
    if type(feature.location) != SeqFeature.SimpleLocation:
        print(f"# not simple location!")
        print(f"{pos=}")
        print(f"{feature=}")
        w = "not simple location!"
    b, e = feature.location.start, feature.location.end

    # nucleotide position in the cds
    if feature.location.strand == 1:
        cds_pos = pos - b  # 0-based position in cds (from start of translation)
        codon_pos = cds_pos % 3  # 0-based position in codon
        x = pos - codon_pos  # 0-based position of codon start
        codon = seq[x : x + 3]  # codon sequence
        aa = codon.translate()  # aminoacid
        transl_nucl = seq[pos]  # nucleotide in translated direction
    else:
        cds_pos = e - 1 - pos  # 0-based position in cds (from start of translation)
        codon_pos = cds_pos % 3  # 0-based position in codon
        x = pos + codon_pos  # 0-based position of codon start
        codon = seq[x - 2 : x + 1].reverse_complement()  # codon sequence
        aa = codon.translate()  # aminoacid
        transl_nucl = seq[pos : pos + 1].reverse_complement()[
            0
        ]  # nucleotide in translated direction

    codon = "".join(codon)
    n_codon = cds_pos // 3

    res = {
        "cds_pos": int(cds_pos) + 1,  # position of nucleotide in cds (1-based)
        "codon_pos": int(codon_pos) + 1,  # position of nucleotide in codon (1-based)
        "codon": codon,  # codon (string) in translated direction
        "aa": aa[0],  # aminoacid (string)
        "n_aa": int(n_codon) + 1,  # position of aminoacid in protein (1-based)
        "cds_len": len(feature.location),  # length of cds
        "protein_len": len(feature.location) // 3,  # length of protein
        "transl_nucl": transl_nucl,  # nucleotide in translated direction
        "warning": w, # warning if not simple location
    }

    return res


def parse_feature(feature):
    res = {
        "type": feature.type,  # feature type
        "loc_type": type(feature.location).__name__,  # location type
        "start": int(feature.location.start) + 1,  # 1-based start
        "end": int(feature.location.end),  # 1-based end
        "strand": feature.location.strand,  # strand
    }
    match res["type"]:
        case "CDS":
            for key in ["gene", "product", "locus_tag"]:
                if key in feature.qualifiers:
                    res[key] = feature.qualifiers[key][0]
        case "gene":
            for key in ["gene", "locus_tag"]:
                if key in feature.qualifiers:
                    res[key] = feature.qualifiers[key][0]
        case "misc_RNA":
            for key in ["product", "locus_tag"]:
                if key in feature.qualifiers:
                    res[key] = feature.qualifiers[key][0]
        case "tRNA":
            for key in ["product", "locus_tag"]:
                if key in feature.qualifiers:
                    res[key] = feature.qualifiers[key][0]
        case "tmRNA":
            for key in ["product", "locus_tag"]:
                if key in feature.qualifiers:
                    res[key] = feature.qualifiers[key][0]
        case "rRNA":
            for key in ["product", "locus_tag"]:
                if key in feature.qualifiers:
                    res[key] = feature.qualifiers[key][0]
        case _:
            raise ValueError(
                f"Processing of feature type {res['type']} not implemented"
            )
    return res


def parse_features(feature_dict, seq):
    results = []
    for pos, features in feature_dict.items():
        for feature in features:
            res = {}
            res["pos"] = int(pos) + 1  # 1-based position in genome
            res |= parse_feature(feature)
            assert type(pos) == SeqFeature.ExactPosition
            res["nucl_fwd"] = seq[pos]
            if res["type"] == "CDS":
                res |= parse_cds(feature, seq, pos)

            results.append(res)
    df = pd.DataFrame(results)
    return df


def find_prev_next_CDS(fname, record_name, pos):
    # turn to 0-based positions
    pos = SeqFeature.ExactPosition(pos - 1)
    f_right = None
    d_right = np.inf
    f_left = None
    d_left = np.inf
    with open(fname) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            if record.name != record_name:
                continue
            for feature in record.features:
                if feature.type == "source":
                    L = len(feature)
                if feature.type != "CDS":
                    continue
                s = feature.location.start
                e = feature.location.end
                assert s < e
                dr = (s - pos) % L
                dl = (pos - e) % L
                if dr < d_right:
                    f_right = feature
                    d_right = dr
                if dl < d_left:
                    f_left = feature
                    d_left = dl

    return {
        "f_fwd": f_right,
        "d_fwd": d_right,
        "f_rev": f_left,
        "d_rev": d_left,
    }


def parse_off_feature_position(fname, record, positions):
    df = []
    N = len(positions)
    for n, pos in enumerate(positions):

        res = find_prev_next_CDS(fname, record, pos)
        df.append(
            {
                "pos": pos,
                "direction": "fwd",
                "distance": res["d_fwd"],
            }
            | parse_feature(res["f_fwd"])
        )
        df.append(
            {
                "pos": pos,
                "direction": "rev",
                "distance": res["d_rev"],
            }
            | parse_feature(res["f_rev"])
        )

        # progress bar
        if (n % 10 == 0) or (n == N - 1):
            end = "\n" if n == N - 1 else "\r"
            frac = int(round((n + 1) / N, 2) * 100)
            print(f"[{'#'*(frac//2):50s}] {frac}% ({n+1}/{N})", end=end)
    df = pd.DataFrame(df)
    return df


if __name__ == "__main__":
    args = parse_args()
    pos_df = pd.read_csv(args.pos)
    positions = pos_df["position"].unique()

    feature_dict, seq = parse_genbank(args.ref, args.record, positions)
    df_hit = parse_features(feature_dict, seq)

    pos_hit = df_hit["pos"].unique()
    pos_miss = set(positions) - set(pos_hit)
    pos_miss = np.sort(list(pos_miss))

    print(f"positions on feature: {len(pos_hit)}")
    print(f"positions off feature: {len(pos_miss)}")

    df_miss = parse_off_feature_position(args.ref, args.record, pos_miss)

    df_hit.to_csv(args.out_hit, index=False)
    df_miss.to_csv(args.out_miss, index=False)
