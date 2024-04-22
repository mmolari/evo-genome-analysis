# %%
import pandas as pd
import numpy as np
from Bio import SeqIO, SeqFeature
from collections import defaultdict


def parse_genbank(fname, record_name, pos_list):
    # turn to 0-based positions
    positions = [SeqFeature.ExactPosition(p - 1) for p in pos_list]
    results = defaultdict(list)
    seq = None
    record_found = False
    with open(fname) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            if record.id != record_name:
                continue
            record_found = True
            for feature in record.features:
                if feature.type == "source":
                    seq = feature.extract(record.seq)
                    continue
                for pos in positions:
                    if int(pos) in feature:
                        results[pos].append(feature)

    assert record_found, f"Record {record_name} not found in {fname}"
    return results, seq


def parse_cds(feature, seq, pos):
    """
    extract:
    - aminoacid
    - aminoacid position
    - codon
    - position in codon
    """
    # 0-based
    assert type(pos) == SeqFeature.ExactPosition
    assert (
        type(feature.location) == SeqFeature.SimpleLocation
    )  # only simple location implemented
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
    }

    return res


def parse_features(feature_dict, seq):
    results = []
    for pos, features in feature_dict.items():
        for feature in features:
            res = {
                "pos": int(pos) + 1,  # 1-based position in genome
                "type": feature.type,  # feature type
                "loc_type": type(feature.location).__name__,  # location type
                "start": int(feature.location.start) + 1,  # 1-based start
                "end": int(feature.location.end),  # 1-based end
                "strand": feature.location.strand,  # strand
            }
            assert type(pos) == SeqFeature.ExactPosition
            res["nucl_fwd"] = seq[pos]

            match res["type"]:
                case "CDS":
                    for key in ["gene", "product", "locus_tag"]:
                        if key in feature.qualifiers:
                            res[key] = feature.qualifiers[key][0]
                    res |= parse_cds(feature, seq, pos)
                case "gene":
                    for key in ["gene", "locus_tag"]:
                        if key in feature.qualifiers:
                            res[key] = feature.qualifiers[key][0]
                case _:
                    raise ValueError(
                        f"Processing of feature type {res['type']} not implemented"
                    )

            results.append(res)
    df = pd.DataFrame(results)
    return df


# positions in 1-based numbering
pos = [4043149, 2832558, 2953318]

fname = "../test_data/input-test-leo/references/ST131-I.gbk"
rec_name = "R1_4963_kbp"
res, seq = parse_genbank(
    fname,
    rec_name,
    pos,
)

fdf = parse_features(res, seq)
fdf

# %%
