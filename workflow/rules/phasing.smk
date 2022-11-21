PHASING = "results/phasing"


rule phasing_ref1:
    output:
        haps=temp(os.path.join(PHASING, "refpanel1_{chrom}.haps")),
        sample=temp(os.path.join(PHASING, "refpanel1_{chrom}.sample")),
        vcf=os.path.join(PHASING, "refpanel1_{chrom}.vcf.gz"),
    log:
        os.path.join(PHASING, "phasing_ref1_{chrom}.llog"),
    params:
        bfile=lambda wildcards: PHASING + "/refpanel1_" + wildcards.chrom,
        vcf=lambda wildcards: os.path.abspath(REFPANEL[wildcards.chrom]["vcf"]),
        a1=temp(lambda wildcards: REFPANEL[wildcards.chrom]["vcf"] + ".txt"),
        maps=lambda wildcards: REFPANEL[wildcards.chrom]["geneticmap"],
        phased=lambda wildcards: REFPANEL[wildcards.chrom]["phased"],
        fam=lambda wildcards: REFPANEL[wildcards.chrom]["trios"],
        intype=config["phasing"]["input"],
        shapeit2=config["phasing"]["shapeit2"],
        out=lambda wildcards, output: output[0][:-5],
    threads: 40
    shell:
        """
        (
        if [ {params.phased} == "yes" ];then \
            ln -sf {params.vcf} {output.vcf} && {BCFTOOLS} index -f {output.vcf} && {BCFTOOLS} convert --hapsample {output.haps},{output.sample} {params.vcf} \
        ; else \
            {BCFTOOLS} query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n' {params.vcf} > {params.a1} && \
            {PLINK} --{params.intype} {params.vcf} --a1-allele {params.a1} 4 3 \# --make-bed --out {params.bfile} --allow-extra-chr --output-chr chr26 && \
            if [ -s {params.fam} ];then cp {params.fam} {params.bfile}.fam;fi && \
            {SHAPEIT2} {params.shapeit2} -B {params.bfile} -M {params.maps} -O {params.out} --thread {threads} && \
            awk 'NR>2 {{$1=$2; $4=0; $5=0}};1' {params.out}.sample > {params.out}.samples && \
            awk '$1="{wildcards.chrom}:"$3"_"$4"_"$5, $2="{wildcards.chrom}:"$3"_"$4"_"$5' {params.out}.haps > {params.out}.haps.tmp && mv {params.out}.haps.tmp {params.out}.haps &&  \
            {BCFTOOLS} convert --vcf-ids --hapsample2vcf {params.out}.haps,{params.out}.samples | {BCFTOOLS} annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz -o {output.vcf} --threads 4 && {BCFTOOLS} index -f {output.vcf} \
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
        os.path.join(PHASING, "prepare_ref1_{chrom}.llog"),
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
        vcf=os.path.join(PHASING, "refpanel2_{chrom}.vcf.gz"),
    log:
        os.path.join(PHASING, "phasing_ref2_{chrom}.llog"),
    params:
        bfile=lambda wildcards: PHASING + "/refpanel2_" + wildcards.chrom,
        vcf=lambda wildcards: os.path.abspath(REFPANEL2[wildcards.chrom]["vcf"]),
        a1=temp(lambda wildcards: REFPANEL2[wildcards.chrom]["vcf"] + ".txt"),
        maps=lambda wildcards: REFPANEL2[wildcards.chrom]["geneticmap"],
        phased=lambda wildcards: REFPANEL2[wildcards.chrom]["phased"],
        fam=lambda wildcards: REFPANEL2[wildcards.chrom]["trios"],
        intype=config["phasing"]["input"],
        shapeit2=config["phasing"]["shapeit2"],
        out=lambda wildcards, output: output[0][:-5],
    threads: 40
    shell:
        """
        (
        if [ {params.phased} == "yes" ];then \
            ln -sf {params.vcf} {output.vcf} && {BCFTOOLS} index -f {output.vcf} && {BCFTOOLS} convert --hapsample {output.haps},{output.sample} {params.vcf} \
        ; else \
            {BCFTOOLS} query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n' {params.vcf} > {params.a1} && \
            {PLINK} --{params.intype} {params.vcf} --a1-allele {params.a1} 4 3 \# --make-bed --out {params.bfile} --allow-extra-chr --output-chr chr26 && \
            if [ -s {params.fam} ];then cp {params.fam} {params.bfile}.fam;fi && \
            {SHAPEIT2} {params.shapeit2} -B {params.bfile} -M {params.maps} -O {params.out} --thread {threads} && \
            awk 'NR>2 {{$1=$2; $4=0; $5=0}};1' {params.out}.sample > {params.out}.samples && \
            awk '$1="{wildcards.chrom}:"$3"_"$4"_"$5 && $2="{wildcards.chrom}:"$3"_"$4"_"$5' {params.out}.haps > {params.out}.haps.tmp && mv {params.out}.haps.tmp {params.out}.haps &&  \
            {BCFTOOLS} convert --vcf-ids --hapsample2vcf {params.out}.haps,{params.out}.samples | {BCFTOOLS} annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz -o {output.vcf} --threads 4 && {BCFTOOLS} index -f {output.vcf} \
        ; fi
        ) &> {log}
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
        os.path.join(PHASING, "prepare_ref2_{chrom}.llog"),
    params:
        haps=lambda wildcards, input: input[0][:-5],
    shell:
        """
        {SHAPEIT2} -convert --input-haps {params.haps} --output-ref {output.hap} {output.leg} {output.sam} &> {log}
        """
