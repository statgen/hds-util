# hds-util

## Installation

```bash
# cget can be installed with `pip3 install --user cget`
cget install -f ./requirements.txt
mkdir build; cd build
cmake -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake -DCMAKE_BUILD_TYPE=Release ..
make
make install
```

# Usage

```
# Generate GT and DS format fields and convert to BCF file format.
hds-util in.sav -f GT,GP -O bcf -o out.bcf

# Paste samples together and recompute estimatad r-square across all samples.
hds-util in1.sav in2.sav in3.sav > merged.sav

# Paste samples, genearte GT and DS while keeping HDS, and filter variants with R2<0.1.  
hds-util -f GT,DS,HDS -m 0.1 in1.sav in2.sav in3.sav > merged.sav
```

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
