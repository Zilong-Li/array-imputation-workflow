

rule phasing_ref1:
    output:
        haps=os.path.join(OUTDIR, "refpanel1_{chrom}.haps"),
        samples=os.path.join(OUTDIR, "refpanel1_{chrom}.samples"),
    params:
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
        a1=temp(lambda wildcards: REFPANEL[wildcards.chrom]["vcf"] + ".txt"),
        maps=lambda wildcards: REFPANEL[wildcards.chrom]["geneticmap"],
        phased=lambda wildcards: REFPANEL[wildcards.chrom]["phased"],
        intype=config["phasing"]["input"],
        shapeit2=config["phasing"]["shapeit2"],
        bfile=OUTDIR + "/refpanel1",
        out=lambda wildcards, output: output[0][:-4],
    threads: 40
    shell:
        """
        if [ {params.phased} == "yes" ];then bcftools convert --hapsample {input} {output}; \
        else \
            {BCFTOOLS} {params.vcf} query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > {params.a1}
            {PLINK} --{params.intype} {params.vcf} --a1-allele {params.a1} 4 3 # --make-bed --out {params.bfile}; \
            {SHAPEIT2} {params.shapeit2} -B {params.bfile} -M {params.maps} -O {params.out} --thread {threads}; \
        fi
        """
