
OUTDIR = "results/subset_refpanel"

def get_samples_list_comma():
    return ",".join(SAMPLES.keys())


def get_ref_vcf(wildcards):
    return REFPANEL[wildcards.chrom]["vcf"]


rule subset_refpanel:
    input:
        get_ref_vcf,
    output:
        [
            os.path.join(OUTDIR, "{chrom}." + ext)
            for ext in ["bcf", "hap.gz", "legend.gz", "sites.vcf.gz", "sites.tsv.gz"]
        ],
    params:
        samples=get_samples_list_comma(),
    threads: 2
    shell:
        """
        ./workflow/scripts/prep-refs.sh {wildcards.chrom} {input} {OUTDIR} ^{params.samples} {threads}
        """
