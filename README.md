# Snakemake Imputation Workflow For Chip Array Data

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)


## Usage

Make sure you have all dependencies solved. Also, make sure all data are using the same reference coordinate, especially paying attention to the chromosome names. All inputs are configured by [config/config.yaml](config/config.yaml). Check [config/README.md](config/README.md) for details.

## Dependencies

- shapeit2
- impute2
- bcftools
- plink (1.9)
- pandas (python)
- awk (linux)
- [dosage](https://github.com/Zilong-Li/vcfpp/blob/main/tools/dosage.cpp)
