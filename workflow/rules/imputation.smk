IMPUTATION = "results/imputation"
DATANAME = os.path.basename(config["bed"])[:-4]


def get_impute2_output_chunks(wildcards):
    starts, ends = get_regions_list_per_chrom(wildcards.chrom)
    gen = expand(
        rules.run_impute2_byregion.output.gen,
        zip,
        start=starts,
        end=ends,
        allow_missing=True,
    )
    info = expand(
        rules.run_impute2_byregion.output.info,
        zip,
        start=starts,
        end=ends,
        allow_missing=True,
    )
    info2 = expand(
        rules.run_impute2_byregion.output.info2,
        zip,
        start=starts,
        end=ends,
        allow_missing=True,
    )
    haps = expand(
        rules.run_impute2_byregion.output.haps,
        zip,
        start=starts,
        end=ends,
        allow_missing=True,
    )
    probs = expand(
        rules.run_impute2_byregion.output.probs,
        zip,
        start=starts,
        end=ends,
        allow_missing=True,
    )
    return {"gen": gen, "info": info, "info2": info2, "haps": haps, "probs": probs}


rule split_chrs:
    input:
        bed=config["bed"],
        bim=config["bim"],
        fam=config["fam"],
    output:
        bed=os.path.join(IMPUTATION, DATANAME + "_{chrom}.bed"),
        bim=os.path.join(IMPUTATION, DATANAME + "_{chrom}.bim"),
        fam=os.path.join(IMPUTATION, DATANAME + "_{chrom}.fam"),
    log:
        os.path.join(IMPUTATION, DATANAME + "_{chrom}.llog"),
    params:
        bfile=lambda wildcards, input: input[0][:-4],
        out=lambda wildcards, output: output[0][:-4],
    shell:
        """
        {PLINK} --bfile {params.bfile} --allow-extra-chr --chr {wildcards.chrom} --keep-allele-order --make-bed --out {params.out} &> {log}
        """


rule check_alignment:
    input:
        bed=rules.split_chrs.output.bed,
        bim=rules.split_chrs.output.bim,
        fam=rules.split_chrs.output.fam,
        hap=rules.prepare_ref1.output.hap,
        leg=rules.prepare_ref1.output.leg,
        sam=rules.prepare_ref1.output.sam,
    output:
        os.path.join(IMPUTATION, "checks", "{chrom}.snp.strand.exclude"),
    log:
        os.path.join(IMPUTATION, "checks", "{chrom}.check_alignment.llog"),
    params:
        maps=lambda wildcards: REFPANEL[wildcards.chrom]["geneticmap"],
        bfile=lambda wildcards, input: input[0][:-4],
        out=lambda wildcards, output: output[0][:-19],
    shell:
        """
        {SHAPEIT2} -check -B {params.bfile} -R {input.hap} {input.leg} {input.sam} -M {params.maps} --output-log {params.out} || true &> {log}
        """


rule run_prephasing:
    input:
        bed=rules.split_chrs.output.bed,
        bim=rules.split_chrs.output.bim,
        fam=rules.split_chrs.output.fam,
        hap=rules.prepare_ref1.output.hap,
        leg=rules.prepare_ref1.output.leg,
        sam=rules.prepare_ref1.output.sam,
        excl=rules.check_alignment.output,
    output:
        haps=os.path.join(IMPUTATION, "prephasing", "{chrom}.prephasing.haps"),
        sample=os.path.join(IMPUTATION, "prephasing", "{chrom}.prephasing.sample"),
    log:
        os.path.join(IMPUTATION, "prephasing", "{chrom}.prephasing.llog"),
    params:
        maps=lambda wildcards: REFPANEL[wildcards.chrom]["geneticmap"],
        bfile=lambda wildcards, input: input[0][:-4],
        out=lambda wildcards, output: output[0][:-5],
        shapeit2=config["imputation"]["shapeit2"],
    threads: 20
    shell:
        """
        {SHAPEIT2} {params.shapeit2} -B {params.bfile} -R {input.hap} {input.leg} {input.sam} -M {params.maps} --exclude-snp {input.excl} -O {params.out} --thread {threads} &> {log}
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


rule run_impute2_byregion:
    input:
        haps=rules.run_prephasing.output.haps,
        hap1=rules.prepare_ref1.output.hap,
        leg1=rules.prepare_ref1.output.leg,
        hap2=rules.prepare_ref2.output.hap,
        leg2=rules.prepare_ref2.output.leg,
    output:
        gen=temp(os.path.join(IMPUTATION, "impute2", "{chrom}-{start}-{end}")),
        info=temp(os.path.join(IMPUTATION, "impute2", "{chrom}-{start}-{end}_info")),
        haps=temp(os.path.join(IMPUTATION, "impute2", "{chrom}-{start}-{end}_haps")),
        summary=temp(
            os.path.join(IMPUTATION, "impute2", "{chrom}-{start}-{end}_summary")
        ),
        probs=temp(
            os.path.join(IMPUTATION, "impute2", "{chrom}-{start}-{end}_allele_probs")
        ),
        info2=temp(
            os.path.join(IMPUTATION, "impute2", "{chrom}-{start}-{end}_info_by_sample")
        ),
    log:
        os.path.join(IMPUTATION, "impute2", "{chrom}-{start}-{end}.llog"),
    params:
        maps=lambda wildcards: REFPANEL[wildcards.chrom]["geneticmap"],
        impute2=config["imputation"]["impute2"],
    benchmark:
        os.path.join(IMPUTATION, "impute2", "{chrom}-{start}-{end}.benchmark.txt")
    threads: 1
    shell:
        """
        (
        {IMPUTE2} {params.impute2} -phase -use_prephased_g -known_haps_g {input.haps} -merge_ref_panels -merge_ref_panels_output_ref {output.gen}_merged_ref -h {input.hap1} {input.hap2} -l {input.leg1} {input.leg2} -m {params.maps} -int {wildcards.start} {wildcards.end} -o {output.gen}
        grep "no SNPs in the imputation interval" {output.summary} && touch {output} || true
        ) &> {log}
        """


rule ligate_impute2_chunks:
    input:
        unpack(get_impute2_output_chunks),
    output:
        gen=os.path.join(IMPUTATION, "impute2", "impute2.{chrom}.gen.gz"),
        haps=os.path.join(IMPUTATION, "impute2", "impute2.{chrom}.haps.gz"),
        info=os.path.join(IMPUTATION, "impute2", "impute2.{chrom}.info.gz"),
        info2=os.path.join(IMPUTATION, "impute2", "impute2.{chrom}.info_by_sample.gz"),
        probs=os.path.join(IMPUTATION, "impute2", "impute2.{chrom}.allele_probs.gz"),
    log:
        os.path.join(IMPUTATION, "impute2", "impute2.{chrom}.ligate.llog"),
    shell:
        """
        ( \
        echo {input.gen} | tr ' ' '\n'  | xargs cat | awk '$2={wildcards.chrom}":"$3"_"$4"_"$5' | gzip -c > {output.gen} && \
        echo {input.info} | tr ' ' '\n'  | xargs cat | gzip -c > {output.info} && \
        echo {input.info2} | tr ' ' '\n'  | xargs cat | gzip -c > {output.info2} && \
        echo {input.haps} | tr ' ' '\n'  | xargs cat | gzip -c > {output.haps} && \
        echo {input.probs} | tr ' ' '\n'  | xargs cat | gzip -c > {output.probs} \
        ) &> {log}
        """


rule convert_formats:
    input:
        gen=rules.ligate_impute2_chunks.output.gen,
        haps=rules.ligate_impute2_chunks.output.haps,
        sample=rules.run_prephasing.output.sample,
    output:
        vcf1=os.path.join(IMPUTATION, "impute2", "impute2.{chrom}.unphased.vcf.gz"),
        vcf2=os.path.join(IMPUTATION, "impute2", "impute2.{chrom}.phased.vcf.gz"),
    log:
        os.path.join(IMPUTATION, "impute2", "impute2.{chrom}.convert_formats.llog"),
    params:
        samples=os.path.join(IMPUTATION, "impute2", "impute2.{chrom}.samples"),
    shell:
        """
        (
        awk 'NR>2 {{$1=$2; $4=0; $5=0}};1' {input.sample} > {params.samples} && \
        {BCFTOOLS} convert -G {input.gen},{params.samples} --threads 4 -Oz -o {output.vcf1} && {BCFTOOLS} index -f {output.vcf1} && \
        {SHAPEIT2} -convert --input-haps {input.haps} --output-vcf  {output.vcf2} && {BCFTOOLS} index -f {output.vcf2} \
        ) &> {log}
        """
