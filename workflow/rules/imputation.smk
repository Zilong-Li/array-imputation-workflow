IMPUTATION = "results/imputation"
DATANAME = os.path.basename(config["bed"])[:-4]


def get_impute2_output_chunks(wildcards):
    starts, ends = get_regions_list_per_chrom(wildcards.chrom)
    if config["phasing"].get("refpanel2"):
        gen = expand(
            rules.run_impute2_byregion_refpanel12.output.gen,
            zip,
            start=starts,
            end=ends,
            allow_missing=True,
        )
        info = expand(
            rules.run_impute2_byregion_refpanel12.output.info,
            zip,
            start=starts,
            end=ends,
            allow_missing=True,
        )
        info2 = expand(
            rules.run_impute2_byregion_refpanel12.output.info2,
            zip,
            start=starts,
            end=ends,
            allow_missing=True,
        )
        haps = expand(
            rules.run_impute2_byregion_refpanel12.output.haps,
            zip,
            start=starts,
            end=ends,
            allow_missing=True,
        )
        probs = expand(
            rules.run_impute2_byregion_refpanel12.output.probs,
            zip,
            start=starts,
            end=ends,
            allow_missing=True,
        )
    else:
        gen = expand(
            rules.run_impute2_byregion_refpanel1.output.gen,
            zip,
            start=starts,
            end=ends,
            allow_missing=True,
        )
        info = expand(
            rules.run_impute2_byregion_refpanel1.output.info,
            zip,
            start=starts,
            end=ends,
            allow_missing=True,
        )
        info2 = expand(
            rules.run_impute2_byregion_refpanel1.output.info2,
            zip,
            start=starts,
            end=ends,
            allow_missing=True,
        )
        haps = expand(
            rules.run_impute2_byregion_refpanel1.output.haps,
            zip,
            start=starts,
            end=ends,
            allow_missing=True,
        )
        probs = expand(
            rules.run_impute2_byregion_refpanel1.output.probs,
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


rule run_impute2_byregion_refpanel1:
    input:
        haps=rules.run_prephasing.output.haps,
        hap1=rules.prepare_ref1.output.hap,
        leg1=rules.prepare_ref1.output.leg,
    output:
        gen=temp(os.path.join(IMPUTATION, "impute2", "refpanel1", "{chrom}-{start}-{end}")),
        info=temp(os.path.join(IMPUTATION, "impute2", "refpanel1", "{chrom}-{start}-{end}_info")),
        haps=temp(os.path.join(IMPUTATION, "impute2", "refpanel1", "{chrom}-{start}-{end}_haps")),
        summary=temp(
            os.path.join(IMPUTATION, "impute2", "refpanel1", "{chrom}-{start}-{end}_summary")
        ),
        probs=temp(
            os.path.join(IMPUTATION, "impute2", "refpanel1", "{chrom}-{start}-{end}_allele_probs")
        ),
        info2=temp(
            os.path.join(IMPUTATION, "impute2", "refpanel1", "{chrom}-{start}-{end}_info_by_sample")
        ),
    log:
        os.path.join(IMPUTATION, "impute2", "refpanel1", "{chrom}-{start}-{end}.llog"),
    params:
        maps=lambda wildcards: REFPANEL[wildcards.chrom]["geneticmap"],
        impute2=config["imputation"]["impute2"],
    benchmark:
        os.path.join(IMPUTATION, "impute2", "{chrom}-{start}-{end}.benchmark.txt")
    threads: 1
    shell:
        """
        (
        {IMPUTE2} {params.impute2} -phase -use_prephased_g -known_haps_g {input.haps} -h {input.hap1} -l {input.leg1} -m {params.maps} -int {wildcards.start} {wildcards.end} -o {output.gen}
        grep "no SNPs in the imputation interval" {output.summary} && touch {output} || true
        ) &> {log}
        """


rule run_impute2_byregion_refpanel12:
    input:
        haps=rules.run_prephasing.output.haps,
        hap1=rules.prepare_ref1.output.hap,
        leg1=rules.prepare_ref1.output.leg,
        hap2=rules.prepare_ref2.output.hap,
        leg2=rules.prepare_ref2.output.leg,
    output:
        gen=temp(os.path.join(IMPUTATION, "impute2", "refpanel12", "{chrom}-{start}-{end}")),
        info=temp(os.path.join(IMPUTATION, "impute2", "refpanel12", "{chrom}-{start}-{end}_info")),
        haps=temp(os.path.join(IMPUTATION, "impute2", "refpanel12", "{chrom}-{start}-{end}_haps")),
        summary=temp(
            os.path.join(IMPUTATION, "impute2", "refpanel12", "{chrom}-{start}-{end}_summary")
        ),
        probs=temp(
            os.path.join(IMPUTATION, "impute2", "refpanel12", "{chrom}-{start}-{end}_allele_probs")
        ),
        info2=temp(
            os.path.join(IMPUTATION, "impute2", "refpanel12", "{chrom}-{start}-{end}_info_by_sample")
        ),
    log:
        os.path.join(IMPUTATION, "impute2", "refpanel12", "{chrom}-{start}-{end}.llog"),
    params:
        maps=lambda wildcards: REFPANEL[wildcards.chrom]["geneticmap"],
        impute2=config["imputation"]["impute2"],
    benchmark:
        os.path.join(IMPUTATION, "impute2", "{chrom}-{start}-{end}.benchmark.txt")
    threads: 1
    shell:
        """
        (
        {IMPUTE2} {params.impute2} -phase -use_prephased_g -known_haps_g {input.haps} -merge_ref_panels -merge_ref_panels_output_ref {output.gen}_merged_ref -h {input.hap2} {input.hap1} -l {input.leg2} {input.leg1} -m {params.maps} -int {wildcards.start} {wildcards.end} -o {output.gen}
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
        echo {input.gen} | tr ' ' '\n'  | xargs cat | awk '$2="{wildcards.chrom}:"$3"_"$4"_"$5' | gzip -c > {output.gen} && \
        echo {input.info} | tr ' ' '\n'  | xargs cat | awk 'NR>1 && $1="{wildcards.chrom}:"$3"_"$4"_"$5;1' |gzip -c > {output.info} && \
        echo {input.info2} | tr ' ' '\n'  | xargs cat | gzip -c > {output.info2} && \
        echo {input.haps} | tr ' ' '\n'  | xargs cat | awk '$1="{wildcards.chrom}:"$3"_"$4"_"$5, $2="{wildcards.chrom}:"$3"_"$4"_"$5' |gzip -c > {output.haps} && \
        echo {input.probs} | tr ' ' '\n'  | xargs cat |awk '$1="{wildcards.chrom}:"$3"_"$4"_"$5' | gzip -c > {output.probs} \
        ) &> {log}
        """


rule convert_impute2_formats:
    input:
        gen=rules.ligate_impute2_chunks.output.gen,
        haps=rules.ligate_impute2_chunks.output.haps,
        sample=rules.run_prephasing.output.sample,
    output:
        unphased=os.path.join(IMPUTATION, "impute2", "impute2.{chrom}.unphased.vcf.gz"),
        unphasedcsi=os.path.join(
            IMPUTATION, "impute2", "impute2.{chrom}.unphased.vcf.gz.csi"
        ),
        phased=os.path.join(IMPUTATION, "impute2", "impute2.{chrom}.phased.vcf.gz"),
        phasedcsi=os.path.join(
            IMPUTATION, "impute2", "impute2.{chrom}.phased.vcf.gz.csi"
        ),
    log:
        os.path.join(IMPUTATION, "impute2", "impute2.{chrom}.convert_formats.llog"),
    params:
        samples=os.path.join(IMPUTATION, "impute2", "impute2.{chrom}.samples"),
    shell:
        """
        (
        awk 'NR>2 {{$1=$2; $4=0; $5=0}};1' {input.sample} > {params.samples} && \
        {BCFTOOLS} convert -G {input.gen},{params.samples} | {dosage} -i - | {BCFTOOLS} annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz -o {output.unphased} && {BCFTOOLS} index -f {output.unphased} && \
        {BCFTOOLS} convert --vcf-ids --hapsample2vcf {input.haps},{params.samples} | {BCFTOOLS} annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz -o {output.phased} && {BCFTOOLS} index -f {output.phased} \
        ) &> {log}
        """
