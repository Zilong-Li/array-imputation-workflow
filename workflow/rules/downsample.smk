
OUTDIR = "results/downsample"


rule downsample_bam:
    output:
        bam=os.path.join(OUTDIR, "{sample}.{depth}.bam"),
        bai=os.path.join(OUTDIR, "{sample}.{depth}.bam.bai"),
    run:
        inbam = SAMPLE[wildcards.sample]["bam"]
        depbam = SAMPLE[wildcards.sample]["depth"]
        frac = wildcards.depth / depbam
        shell(
            """ samtools view -s {frac} -o {output.bam} {inbam} && samtools index {output.bam} """
        )
