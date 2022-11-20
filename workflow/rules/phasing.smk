PHASING = "results/phasing"


rule phasing_ref1:
    output:
        haps=temp(os.path.join(PHASING, "refpanel1_{chrom}.haps")),
        sample=temp(os.path.join(PHASING, "refpanel1_{chrom}.sample")),
    log:
        os.path.join(PHASING, "refpanel1_{chrom}.llog"),
    params:
        bfile=lambda wildcards: PHASING + "/refpanel1_" + wildcards.chrom,
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
        a1=temp(lambda wildcards: REFPANEL[wildcards.chrom]["vcf"] + ".txt"),
        maps=lambda wildcards: REFPANEL[wildcards.chrom]["geneticmap"],
        phased=lambda wildcards: REFPANEL[wildcards.chrom]["phased"],
        fam=lambda wildcards: REFPANEL[wildcards.chrom]["trios"],
        intype=config["phasing"]["input"],
        shapeit2=config["phasing"]["shapeit2"],
        out=lambda wildcards, output: output[0][:-5],
    threads: 20
    shell:
        """
        (
        if [ {params.phased} == "yes" ];then \
            bcftools convert --hapsample {output.haps},{output.sample} {params.vcf} \
        ; else \
            {BCFTOOLS} query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n' {params.vcf} > {params.a1} && \
            {PLINK} --{params.intype} {params.vcf} --a1-allele {params.a1} 4 3 \# --make-bed --out {params.bfile} && \
            if [ -s {params.fam} ];then cp {params.fam} {params.bfile}.fam;fi && \
            {SHAPEIT2} {params.shapeit2} -B {params.bfile} -M {params.maps} -O {params.out} --thread {threads} \
        ; fi
        ) &> {log}
        """


rule prepare_ref1:
    input:
        haps=rules.phasing_ref1.output.haps,
        sample=rules.phasing_ref1.output.sample,
    output:
        hap=os.path.join(PHASING, "refpanel1_{chrom}.hap"),
        leg=os.path.join(PHASING, "refpanel1_{chrom}.leg"),
        sam=os.path.join(PHASING, "refpanel1_{chrom}.sam"),
    log:
        os.path.join(PHASING, "refpanel1_{chrom}.llog"),
    params:
        haps=lambda wildcards, input: input[0][:-5],
    shell:
        """
        {SHAPEIT2} -convert --input-haps {params.haps} --output-ref {output.hap} {output.leg} {output.sam} &> {log}
        """


rule phasing_ref2:
    output:
        haps=temp(os.path.join(PHASING, "refpanel2_{chrom}.haps")),
        sample=temp(os.path.join(PHASING, "refpanel2_{chrom}.sample")),
    params:
        bfile=lambda wildcards: PHASING + "/refpanel2_" + wildcards.chrom,
        vcf=lambda wildcards: REFPANEL2[wildcards.chrom]["vcf"],
        a1=temp(lambda wildcards: REFPANEL2[wildcards.chrom]["vcf"] + ".txt"),
        maps=lambda wildcards: REFPANEL2[wildcards.chrom]["geneticmap"],
        phased=lambda wildcards: REFPANEL2[wildcards.chrom]["phased"],
        fam=lambda wildcards: REFPANEL2[wildcards.chrom]["trios"],
        intype=config["phasing"]["input"],
        shapeit2=config["phasing"]["shapeit2"],
        out=lambda wildcards, output: output[0][:-5],
    threads: 20
    shell:
        """
        if [ {params.phased} == "yes" ];then \
            bcftools convert --hapsample {output.haps},{output.sample} {params.vcf} \
        ; else \
            {BCFTOOLS} query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n' {params.vcf} > {params.a1} && \
            {PLINK} --{params.intype} {params.vcf} --a1-allele {params.a1} 4 3 \# --make-bed --out {params.bfile} && \
            if [ -s {params.fam} ];then cp {params.fam} {params.bfile}.fam;fi && \
            {SHAPEIT2} {params.shapeit2} -B {params.bfile} -M {params.maps} -O {params.out} --thread {threads} \
        ; fi
        """

rule prepare_ref2:
    input:
        haps=rules.phasing_ref2.output.haps,
        sample=rules.phasing_ref2.output.sample,
    output:
        hap=os.path.join(PHASING, "refpanel2_{chrom}.hap"),
        leg=os.path.join(PHASING, "refpanel2_{chrom}.leg"),
        sam=os.path.join(PHASING, "refpanel2_{chrom}.sam"),
    log:
        os.path.join(PHASING, "refpanel2_{chrom}.llog"),
    params:
        haps=lambda wildcards, input: input[0][:-5],
    shell:
        """
        {SHAPEIT2} -convert --input-haps {params.haps} --output-ref {output.hap} {output.leg} {output.sam} &> {log}
        """
