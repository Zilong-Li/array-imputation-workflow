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
BCFTOOLS = config['bcftools']
PLINK = config['plink']
SHAPEIT2 = config['shapeit2']
IMPUTE2 = config['impute2']


wildcard_constraints:
    chrom="|".join(REFPANEL.keys()),


def get_all_results():
    return expand(rules.run_prephasing.output, chrom=REFPANEL.keys())


def get_regions_list_per_chrom(wildcards):
    """split chr into chunks given chunksize; return starts, ends' pairs"""
    d = defaultdict(list)
    with open(config["bim"], "r") as f:
        for line in f:
            tmp = line.rstrip().split()
            tmp[0] = tmp[0].replace("chr", "")
            d[tmp[0]].append(int(tmp[3]))
    chrom = f"{wildcards.chrom}"
    pos = sorted(d.get(chrom))
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
