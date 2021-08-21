# hds-util

## Installation

```bash
cget install -f ./requirements.txt
mkdir build; cd build
cmake -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake -DCMAKE_BUILD_TYPE=Release ..
make
make install
```

## R-squared Formulas

TODO

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
