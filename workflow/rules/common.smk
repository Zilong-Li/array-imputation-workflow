import os
import pandas as pd
from collections import defaultdict

# dict : {"chr1": {"vcf": path, "phased": "no", "trios":"path","geneticmap":"path"}, ...}
REFPANEL = (
    pd.read_csv(config["phasing"]["refpanel1"], sep="\t", dtype=str)
    .set_index("chr")
    .to_dict(orient="index")
)

if os.path.exists(config["phasing"]["refpanel2"]):
    REFPANEL2 = (
        pd.read_csv(config["phasing"]["refpanel2"], sep="\t", dtype=str)
        .set_index("chr")
        .to_dict(orient="index")
    )


CHUNKSIZE = config["imputation"]["chunksize"]

# programs
BCFTOOLS = config["bcftools"]
PLINK = config["plink"]
SHAPEIT2 = config["shapeit2"]
IMPUTE2 = config["impute2"]

dpos = defaultdict(list)
with open(config["bim"], "r") as f:
    for line in f:
        tmp = line.rstrip().split()
        dpos[tmp[0]].append(int(tmp[3]))


chroms = REFPANEL.keys()
chroms = ["chr21"]

RUN = config["scenario"]


wildcard_constraints:
    chrom="|".join(REFPANEL.keys()),


def get_all_results():
    if RUN == "all":
        return get_phasing_results(), get_imputation_results(), get_merging_results()
    elif RUN == "phasing":
        return get_phasing_results()
    elif RUN == "imputation":
        return get_imputation_results()
    elif RUN == "merge":
        return get_merging_results()
    else:
        raise (ValueError("illegal scenario provided.\n"))


def get_imputation_results():
    return expand(rules.convert_impute2_formats.output, chrom=chroms)


def get_phasing_results():
    res = expand(rules.prepare_ref1.output, chrom=chroms)
    if os.path.exists(config["phasing"]["refpanel2"]):
        res.append(expand(rules.prepare_ref2.output, chrom=chroms))
    return res


def get_merging_results():
    res = expand(rules.merge_unphasedvcf_ref1.output, chrom=chroms)
    res.append(expand(rules.merge_phasedvcf_ref1.output, chrom=chroms))
    if os.path.exists(config["phasing"]["refpanel2"]):
        res.append(expand(rules.merge_unphasedvcf_ref1_ref2.output, chrom=chroms))
        res.append(expand(rules.merge_phasedvcf_ref1_ref2.output, chrom=chroms))
    return res


def get_regions_list_per_chrom(chrom):
    """split chr into chunks given chunksize; return starts, ends' pairs"""
    pos = sorted(dpos.get(chrom))
    s, e = pos[0], pos[-1]
    n = int((e - s) / CHUNKSIZE) + 1
    if (n - 1) * CHUNKSIZE == e - s:
        n = n - 1
    starts, ends = [], []
    for i in range(n):
        ps = CHUNKSIZE * i + s
        pe = CHUNKSIZE * (i + 1) + s - 1
        pe = e if pe > e else pe
        starts.append(ps)
        ends.append(pe)
    return starts, ends
