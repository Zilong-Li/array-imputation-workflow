
MERGING = "results/merge"

rule merge_unphasedvcf_ref1:
    input:
        vcf=rules.convert_impute2_formats.output.unphased,
        ref1=rules.phasing_ref1.output.vcf,
    output:
        vcf=os.path.join(MERGING, "merge_refpanel1_unphased_impute2_{chrom}.vcf.gz"),
    log:
        os.path.join(MERGING, "merge_refpanel1_unphased_impute2_{chrom}.llog"),
    threads: 4
    shell:
        """
        {BCFTOOLS} merge {input.vcf} {input.ref1} | {BCFTOOLS} annotate -x INFO --set-id '%CHROM:%POS:%REF:%FIRST_ALT' --threads {threads} -Oz -o {output.vcf} && {BCFTOOLS} index -f {output.vcf}
        """


rule merge_phasedvcf_ref1:
    input:
        vcf=rules.convert_impute2_formats.output.phased,
        ref1=rules.phasing_ref1.output.vcf,
    output:
        vcf=os.path.join(MERGING, "merge_refpanel1_phased_impute2_{chrom}.vcf.gz"),
    log:
        os.path.join(MERGING, "merge_refpanel1_phased_impute2_{chrom}.llog"),
    threads: 4
    shell:
        """
        {BCFTOOLS} merge {input.vcf} {input.ref1} | {BCFTOOLS} annotate -x INFO --set-id '%CHROM:%POS:%REF:%FIRST_ALT' --threads {threads} -Oz -o {output.vcf} && {BCFTOOLS} index -f {output.vcf}
        """

rule merge_unphasedvcf_ref1_ref2:
    input:
        vcf=rules.convert_impute2_formats.output.unphased,
        ref1=rules.phasing_ref1.output.vcf,
        ref2=rules.phasing_ref2.output.vcf,
    output:
        vcf=os.path.join(MERGING, "merge_refpanel12_unphased_impute2_{chrom}.vcf.gz"),
    log:
        os.path.join(MERGING, "merge_refpanel12_unphased_impute2_{chrom}.llog"),
    threads: 4
    shell:
        """
        {BCFTOOLS} merge {input.vcf} {input.ref1} {input.ref2} | {BCFTOOLS} annotate -x INFO --set-id '%CHROM:%POS:%REF:%FIRST_ALT' --threads {threads} -Oz -o {output.vcf} && {BCFTOOLS} index -f {output.vcf}
        """

rule merge_phasedvcf_ref1_ref2:
    input:
        vcf=rules.convert_impute2_formats.output.phased,
        ref1=rules.phasing_ref1.output.vcf,
        ref2=rules.phasing_ref2.output.vcf,
    output:
        vcf=os.path.join(MERGING, "merge_refpanel12_phased_impute2_{chrom}.vcf.gz"),
    log:
        os.path.join(MERGING, "merge_refpanel12_phased_impute2_{chrom}.llog"),
    threads: 4
    shell:
        """
        {BCFTOOLS} merge {input.vcf} {input.ref1} {input.ref2} | {BCFTOOLS} annotate -x INFO --set-id '%CHROM:%POS:%REF:%FIRST_ALT' --threads {threads} -Oz -o {output.vcf} && {BCFTOOLS} index -f {output.vcf}
        """
