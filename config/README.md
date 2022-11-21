# General Settings

To configure this workflow, modify [config/config.yaml](config.yaml) according to your needs, following the explanations provided in the file.

# Reference Panels Sheet

There are 5 mandatory columns `chr`,`vcf`,`phased`,`trios` and`geneticmap` needed to be defined in the reference panel sheet. The `phased` column only validates `yes` and `no` string. The `trios` column defines a path of plink `fam` file or `na`. The `geneticmap`defines a path of genetic map used by `shapeite2` and `impute2`.

- `refpanel1.tsv` is mandatory. The `vcf` can be phased or unphased. This is primarily used as panel for prephasing the chip array data.
- `refpanel2.tsv` is optional. The `vcf` can be phased or unphased. If this is defined in [config.yaml](config.yaml), then `impute2` will use both `refpanel1` and `refpanel2` as reference panels.

# Input Files

The chip array genotype data must be in `plink` formats required by `shapeite2`. Therefore you must define `bed`, `bim` and `fam` in [config.yaml](config.yaml).

# Scenarios

Several scenarios are defined in [config.yaml](config.yaml). One can choose to run each scenario individually or `all`.


``` sh
snakemake -j60 --config scenario=phasing -n
```

