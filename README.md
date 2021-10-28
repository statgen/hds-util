# hds-util
This is a post-processing tool for [Minimac4](https://github.com/statgen/Minimac4) and [Michigan Imputation Server](https://imputationserver.readthedocs.io/en/latest/) (MIS). It can generate FORMAT fields from HDS, convert from the SAV file format to BCF or VCF, and paste together sample groups that were split due to MIS sample size limit.  

## Installation

```bash
# cget can be installed with `pip3 install --user cget`
cget install -f ./requirements.txt
mkdir build; cd build
cmake -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake -DCMAKE_BUILD_TYPE=Release ..
make
make install
```

## Usage

```
# Generate GT and DS format fields and convert to BCF file format.
hds-util in.sav -f GT,GP -O bcf -o out.bcf

# Paste samples together and recompute estimatad r-square across all samples.
hds-util in1.sav in2.sav in3.sav > merged.sav

# Paste samples, genearte GT and DS while keeping HDS, and filter variants with R2<0.1.  
hds-util -f GT,DS,HDS -m 0.1 in1.sav in2.sav in3.sav > merged.sav
```

## Pasting Samples
To impute datasets that exceed the MIS maximum sample size, array VCFs must be split into sample group files. The imputed sample groups can be paste together using hds-util. The site list for each sample group file must match, so the mininum r-square threshold must be disabled when submitting the imputation job. An r-square filter can be applied in hds-util with `--min-r2 <threshold>`.

## Field Generation Formulas

Below are the formulas for calculating other FORMAT fields from HDS values where `x` is the first haplotype dosage and `y` is the second.

### Diploid
```
DS = x+y
GT = round(x), round(y)
GP = (1-x)(1-y), x(1-y)+y(1-x), xy
SD = x(1-x)+y(1-y)
```

### Haploid
```
DS = x
GT = round(x)
GP = 1-x, x
SD = x(1-x)
```
