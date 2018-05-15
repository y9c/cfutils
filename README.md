# cfutils

**Chromatogram File Utils**

for sanger sequencing data visualizing, alignment, mutation calling, and trimming etc.

## Demo

![plot chromatogram with mutation](https://raw.githubusercontent.com/yech1990/cfutils/master/data/plot.png)

code to generate demo above

```python
import matplotlib.pyplot as plt
from cfutils.align import align
from cfutils.parser import parse_abi
from cfutils.show import highlight_base, plot_chromatograph
seq = parse_abi('./data/B5-M13R_B07.ab1', trim=True)
subject_fasta = './data/3kref.fa'
mutations = align(seq, subject_fasta, ignore_ambig=True)
# select the first one
selected_mutation = mutations[0][2]
fig, ax = plt.subplots(1, 1, figsize=(15, 6))
plot_chromatograph(
    seq, ax, xlim=[selected_mutation - 10, selected_mutation + 10])
highlight_base(selected_mutation, seq, ax)
plt.savefig('./test.pdf')
```

## How to install?

### form pypi

*(use this way ONLY, if you don't know what't going on)*

```bash
pip install cfutils
```

### manipulate the soource code

clone from github

```bash
git clone git@github.com:yech1990/cfutils.git 
```

install dependance

```bash
make init
```

do unittest

```bash
make test
```

## How to use?
 
```python
import cfutils as cf

```

## ChangeLog

- build as python package for pypi
- fix bug that highlihgting wrong base

## TODO

- [ ] call mutation by alignment and plot Chromatogram graphic
- [ ] add a doc
- [x] change xaxis by peak location
