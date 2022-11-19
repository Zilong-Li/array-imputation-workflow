

rule merge_unphased_vcfs:
    input:
        vcf=rules.convert_genhaps2vcf.output.v1,
        ref=,
    output:
        fix=temp(),
        vcf=,
    threads: 4
    shell:
        """
        {FixVCF} {input.ref} {input.vcf} {output.fix} && {BCFTOOLS} index -f {output.fix}
        {BCFTOOLS} merge {input.ref} {output.fix} | {BCFTOOLS} annotate -x INFO --set-id %CHROM:%POS:%REF:%FIRST_ALT --threads {threads} -Oz -o {output.vcf} && {BCFTOOLS} index -f {output.vcf}
        """

rule merge_phased_vcfs:
    input:
        vcf=rules.convert_genhaps2vcf.output.v2,
        ref=,
    output:
        fix=temp(),
        vcf=,
    threads: 4
    shell:
        """
        {FixVCF} {input.ref} {input.vcf} {output.fix} && {BCFTOOLS} index -f {output.fix}
        {BCFTOOLS} merge {input.ref} {output.fix} | {BCFTOOLS} annotate -x INFO --set-id %CHROM:%POS:%REF:%FIRST_ALT --threads {threads} -Oz -o {output.vcf} && {BCFTOOLS} index -f {output.vcf}
        """
