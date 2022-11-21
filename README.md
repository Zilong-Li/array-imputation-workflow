# Snakemake Imputation Workflow For Chip Array Data

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)


## Usage

Make sure you have all dependencies solved. Also, make sure all data are using the same reference coordinate, especially paying attention to the chromosome names. All inputs are configured by [config/config.yaml](config/config.yaml). Check [config/README](config/README.md) for details.

## Dependencies

- [shapeit2](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)
- [impute2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html)
- [bcftools](http://www.htslib.org/download/)
- [plink](https://www.cog-genomics.org/plink/1.9/)(1.9)
- [dosage](https://github.com/Zilong-Li/vcfpp/blob/main/tools/dosage.cpp)
- pandas (python)
- awk (linux)
