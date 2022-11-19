
rule split_chrs:
    input:
        bed=config["bed"],
        bim=config["bim"],
        fam=config["fam"],
    output:
        bed=os.path.join(OUTDIR, DATANAME + "_{chrom}.bed"),
        bim=os.path.join(OUTDIR, DATANAME + "_{chrom}.bim"),
        fam=os.path.join(OUTDIR, DATANAME + "_{chrom}.fam"),
    params:
        bfile=lambda wildcards, input: input[0][:-4],
        out=lambda wildcards, output: output[0][:-4],
    shell:
        """
        {PLINK} --bfile {params.bfile} --allow-extra-chr --chr {wildcards.chrom} --keep-allele-order --make-bed --out {params.out}
        """


rule check_alignment:
    input:
        bed=rules.split_chrs.output.bed,
        bim=rules.split_chrs.output.bim,
        fam=rules.split_chrs.output.fam,
        hap=rules.phasing_ref1.output.hap,
        leg=rules.phasing_ref1.output.leg,
        sam=rules.phasing_ref1.output.sam,
        maps=os.path.join(MapDir, "genetic_map_hg38_chr{chrom}.txt"),
    output:
        os.path.join(OutByChr, OutPrefix + "_{chrom}.alignment.snp.strand.exclude"),
    params:
        log=os.path.join(OutByChr, OutPrefix + "_{chrom}.alignment"),
        bfile=lambda wildcards, input: input[0][:-4],
    shell:
        """
        {SHAPEIT2} -check -B {params.bfile} -R {input.hap} {input.leg} {input.sam} -M {input.maps} --output-log {params.log} || true
        """


rule run_prephasing:
    input:
        bed=rules.split_chrs.output.bed,
        bim=rules.split_chrs.output.bim,
        fam=rules.split_chrs.output.fam,
        hap=rules.phasing_ref1.output.hap,
        leg=rules.phasing_ref1.output.leg,
        sam=rules.phasing_ref1.output.sam,
        maps=os.path.join(MapDir, "genetic_map_hg38_chr{chrom}.txt"),
        exclude=os.path.join(
            OutByChr, OutPrefix + "_chr{chrom}.alignment.snp.strand.exclude"
        ),
    output:
        os.path.join(OutPrephasing, OutPrefix + "_chr{chrom}.prephased.sample"),
        os.path.join(OutPrephasing, OutPrefix + "_chr{chrom}.prephased.haps"),
    params:
        refs=os.path.join(InPhasing, WgsPrefix + "_{chrom}"),
        bed=os.path.join(OutByChr, OutPrefix + "_chr{chrom}"),
        maps=os.path.join(MapDir, "genetic_map_hg38_chr{chrom}.txt"),
        out=os.path.join(OutPrephasing, OutPrefix + "_chr{chrom}.prephased"),
    threads: 40
    shell:
        """
        {SHAPEIT2} --duohmm -W 5 -B {params.bed} -R {params.refs}.hap {params.refs}.leg {params.refs}.sam -M {params.maps} --exclude-snp {input} -O {params.out} --thread {threads}
        """


rule convert_format:
    input:
        os.path.join(OutPrephasing, OutPrefix + "_chr{chrom}.prephased.haps"),
    output:
        os.path.join(OutPrephasing, OutPrefix + "_chr{chrom}.prephased.vcf.gz"),
    params:
        out=os.path.join(OutPrephasing, OutPrefix + "_chr{chrom}.prephased"),
    shell:
        """
        {SHAPEIT2} -convert --input-haps {params.out} --output-vcf {params.out}.vcf.gz && \
        gzip -dc {params.out}.vcf.gz | bgzip -c >{params.out}.t.vcf.gz && mv {params.out}.t.vcf.gz {params.out}.vcf.gz && \
        tabix -f {params.out}.vcf.gz && echo done
        """


def get_imputed_chunks(wildcards):
    """get all imputed results for chunks of 5M regions"""
    chrom = f"{wildcards.chrom}"
    gens = expand(
        str(rules.run_imputation_byregion.output.gen),
        region=get_chunks_list(chrom),
        allow_missing=True,
    )
    info = [f"{gen}_info" for gen in gens]
    haps = [f"{gen}_haps" for gen in gens]
    allele_probs = [f"{gen}_allele_probs" for gen in gens]
    return {"gen": gens, "info": info, "haps": haps, "allele_probs": allele_probs}



rule run_imputation_byregion:
    input:
        haps=os.path.join(OutPrephasing, OutPrefix + "_chr{chrom}.prephased.haps"),
        hap1=os.path.join(In1KGP, "chr{chrom}.forMega.hap.gz"),
        leg1=os.path.join(In1KGP, "chr{chrom}.forMega.legend.gz"),
        hap2=os.path.join(InPhasing, WgsPrefix + "_{chrom}.hap"),
        leg2=os.path.join(InPhasing, WgsPrefix + "_{chrom}.leg"),
        maps=os.path.join(MapDir, "genetic_map_hg38_chr{chrom}.txt"),
    output:
        gen=temp(os.path.join(OutWithBoth, OutPrefix + ".impute2.chr{chrom}-{region}")),
        info=temp(
            os.path.join(OutWithBoth, OutPrefix + ".impute2.chr{chrom}-{region}_info")
        ),
        haps=temp(
            os.path.join(OutWithBoth, OutPrefix + ".impute2.chr{chrom}-{region}_haps")
        ),
        allele_probs=temp(
            os.path.join(
                OutWithBoth, OutPrefix + ".impute2.chr{chrom}-{region}_allele_probs"
            )
        ),
        summary=temp(
            os.path.join(
                OutWithBoth, OutPrefix + ".impute2.chr{chrom}-{region}_summary"
            )
        ),
        warning=temp(
            os.path.join(
                OutWithBoth, OutPrefix + ".impute2.chr{chrom}-{region}_warnings"
            )
        ),
        infobysample=temp(
            os.path.join(
                OutWithBoth, OutPrefix + ".impute2.chr{chrom}-{region}_info_by_sample"
            )
        ),
    threads: 1
    benchmark:
        os.path.join(
            OutWithBoth, OutPrefix + ".impute2.chr{chrom}-{region}.benchmark.txt"
        )
    run:
        rg = f"{wildcards.region}"
        ps = rg.split("-")[0]
        pe = rg.split("-")[1]
        shell(
            """
        impute2 -phase -use_prephased_g -known_haps_g {input.haps} -h {input.hap2} -l {input.leg2} -m {input.maps} -int {ps} {pe} -Ne 20000 -o {output}
        grep "no SNPs in the imputation interval" {output}_summary && touch {output} || true
        """
        )



rule concat_imputed_chunks:
    input:
        unpack(get_imputed_chunks),
    output:
        gen=protected(os.path.join(OutWithBoth, OutPrefix + ".impute2.chr{chrom}.gz")),
        info=protected(
            os.path.join(OutWithBoth, OutPrefix + ".impute2.chr{chrom}.info.gz")
        ),
        haps=protected(
            os.path.join(OutWithBoth, OutPrefix + ".impute2.chr{chrom}.haps.gz")
        ),
        allele_probs=protected(
            os.path.join(
                OutWithBoth, OutPrefix + ".impute2.chr{chrom}.allele_probs.gz"
            )
        ),
        lst=os.path.join(OutWithBoth, OutPrefix + ".impute2.chr{chrom}.chunks.list"),
    run:
        gather = zip(input.gen, input.info, input.haps, input.allele_probs)
        with open(output.lst, "w") as filelist:
            for gen, info, haps, allele_probs in gather:
                print("\t".join([gen, info, haps, allele_probs]), file=filelist)
        shell(
            """
        awk '{{print $1;}}' {output.lst} | xargs cat |gzip -c > {output.gen}
        awk '{{print $2;}}' {output.lst} | xargs cat |gzip -c > {output.info}
        awk '{{print $3;}}' {output.lst} | xargs cat |gzip -c > {output.haps}
        awk '{{print $4;}}' {output.lst} | xargs cat |gzip -c > {output.allele_probs}
        """
        )


rule convert_genhaps2vcf:
    input:
        haps=os.path.join(OutWithBoth, OutPrefix + ".impute2.chr{chrom}.haps.gz"),
        gen=os.path.join(OutWithBoth, OutPrefix + ".impute2.chr{chrom}.gz"),
        sample=os.path.join(OutPrephasing, OutPrefix + "_chr{chrom}.prephased.sample"),
    output:
        v1=os.path.join(OutWithBoth, OutPrefix + ".impute2.chr{chrom}.unphased.vcf.gz"),
        v2=os.path.join(OutWithBoth, OutPrefix + ".impute2.chr{chrom}.phased.vcf.gz"),
        v3=temp(
            os.path.join(OutWithBoth, OutPrefix + ".impute2.chr{chrom}.phased.tmp.gz")
        ),
    params:
        os.path.join(OutWithBoth, OutPrefix + ".impute2.chr{chrom}"),
    shell:
        """
        {PLINK} --gen {input.gen} --sample {input.sample} --hard-call-threshold 0.1 --oxford-single-chr {wildcards.chrom} --recode vcf-iid bgz --out {params}.unphased && {BCFTOOLS} index -f {output.v1}
        cp -f {input.sample} {params}.sample
        gzip -dc {input.haps} >{params}.haps && \
        {Shapeit} -convert --input-haps {params} --output-vcf {output.v2} && rm -f {params}.haps && \
        zcat {output.v2} | awk 'OFS="\t" {{if($0!~/^#/)$1={wildcards.chrom};print}}' | bgzip -c >{output.v3} && mv {output.v3} {output.v2} && \
        {BCFTOOLS} index -f {output.v2} && echo done
        """


rule calc_info_af:
    input:
        gen=os.path.join(OutWithBoth, OutPrefix + ".impute2.chr{chrom}.gz"),
        info=os.path.join(OutWithBoth, OutPrefix + ".impute2.chr{chrom}.info.gz"),
    output:
        os.path.join(OutWithBoth, OutPrefix + ".impute2.chr{chrom}.af.gz"),
    shell:
        """
        {CalcuAF} {input.gen} {input.info} {wildcards.chrom} |gzip -c > {output}
        """
