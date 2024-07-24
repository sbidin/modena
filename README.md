# modena

Modena is a nanopore-based computational method for detecting a wide spectrum
of epigenetic and epitranscriptomic modifications.

It uses an unsupervised learning approach, namely resampling of nanopore
signals followed by the Kuiper test. Unlike other unsupervised tools,
classification is performed by 1D clustering of scores into two groups.

### setup
```shell
$ git clone https://github.com/sbidin/modena.git
$ cd modena
$ poetry install
$ poetry run python -m modena --help # See options.
```

### inputs
Both datasets need to be supplied in `blow5` or `slow5` format, alongside their
`f5c resquiggle` output `tsv` files. If your dataset is in single/multi `fast5`
format, or `pod5` format, you can apply conversions using one of the following
tools:

* single-`fast5` to multi-`fast5`: [ont_fast5_api](https://github.com/nanoporetech/ont_fast5_api?tab=readme-ov-file#single_to_multi_fast5)
* multi-`fast5` to `blow5`/`slow5`: [slow5tools](https://github.com/hasindu2008/slow5tools?tab=readme-ov-file#usage)
* `pod5` to `blow5`/`slow5`: [blue-crab](https://github.com/Psy-Fer/blue-crab?tab=readme-ov-file#usage)

To resquiggle your data with `f5c`, [install f5c](https://hasindu2008.github.io/f5c/docs/quick-start) and run [the resquiggle command](https://hasindu2008.github.io/f5c/docs/commands#resquiggle):
```shell
$ f5c resquiggle data.fastq data.blow5 > resquiggled.tsv
```

### example modena usage
Both datasets need a `blow5` or `slow5` file and a corresponding
`f5c`-resquiggled `tsv` file.
```shell
$ poetry run python -m modena -1 a.blow5 a.tsv -2 b.blow5 b.tsv -o out.tsv
$ poetry run python -m modena --help # See here for more options.
```

### output format
Modena outputs a simple `tsv` file with four columns:
* position, `int`, 1-based
* coverage, `int`, a count of all reads that contributed to the signal
* distance, `float`, a two-sample Kuiper-test measure
* label, `str`, `"pos"` or `"neg"`, separating positions into two clusters

### notes
Modena v1 supported only `fast5` files resquiggled by
[tombo](https://nanoporetech.github.io/tombo/index.html). If you have an older
dataset processed by tombo and wish to use it, visit [tag v1.0.0](https://github.com/sbidin/modena/tree/v1.0.0).
