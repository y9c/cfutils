# cfutils

**Chromatogram File Utils**

for sanger sequencing data visualizing, alignment, mutation calling, and trimming etc.

## Demo

![plot chromatogram with mutation](https://raw.githubusercontent.com/yech1990/cfutils/master/data/plot.png)

code to generate demo above

```python
import matplotlib.pyplot as plt
from cfutils.parser import parse_abi, parse_fasta
from cfutils.align import align
from cfutils.show import highlight_base, plot_chromatograph

query_record = parse_abi('./data/B5-M13R_B07.ab1', trim=False)
subject_record = parse_fasta('./data/3kref.fa')
mutations = align(query_record, subject_record, ignore_ambig=True)
selected_mutation = mutations[5][2]
print(selected_mutation)
fig, ax = plt.subplots(1, 1, figsize=(15, 6))
plot_chromatograph(
    query_record, ax, xlim=[selected_mutation - 10, selected_mutation + 10])
highlight_base(selected_mutation, query_record, ax)
plt.savefig('./test.pdf')
```

## How to install?

### form pypi

*(use this way ONLY, if you don't know what't going on)*

```bash
pip install cfutils
```

### manipulate the source code

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

in command line

```bash
cfutils mut --query ./data/B5-M13R_B07.ab1 --subject data/3kref.fa
```

as python modual

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
- [ ] fix bug that chromatogram switch pos after trim
- [x] wrap as a cli app
- [ ] return quality score in output
- [ ] fix issue that selected base is not in the middle
