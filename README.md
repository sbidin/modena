# nodclust

A nanopore-based computational method for detecting a wide spectrum of
epigenetic/epitranscriptomic modifications.

This application is the implementation of an unsupervised approach for detecting
a diverse range of epigenetic and epitranscriptomic modifications. Using a
combination of the Kuiper test and Jenks 1D clustering, this method exhibits an
average F1-score of 0.746, significantly outperforming rival methods whose
average F1-scores fluctuate between 0.244 and 0.512. We depart from the common
paradigm by using one-dimensional clustering to establish the classification
threshold.

### Install and run locally
```shell
$ python -m pip install -e . # Install.
$ python -m nodclust path/dataset1 path/dataset2 output.bed # Run.
$ python -m nodclust --help # See options.
```

### Install and run via docker
```shell
$ docker build . -t nodclust # Build docker image.
$ docker run -v `pwd`:`pwd` -w `pwd` --rm -it nodclust ./dataset1 ./dataset2 output.bed # Run.
$ docker run --rm -it nodclust --help # See options.
```

### Supported input formats
Only single-FAST5 files processed by
[tombo](https://nanoporetech.github.io/tombo/index.html) are supported. If your
input is in multi-FAST5 format, convert it to single-FAST5 by running the
`multi_to_single_fast5` utility [available
here](https://github.com/nanoporetech/ont_fast5_api). If your inputs don't
contain basecalls, annotate them by running tombo's `preprocess
annotate_raws_with_fastqs` subcommand [as shown
here](https://nanoporetech.github.io/tombo/examples.html?highlight=annotate_raw_with_fastqs).
Finally, make sure your inputs have been run through tombo's `resquiggle`
command [as explained
here](https://nanoporetech.github.io/tombo/examples.html?highlight=resquiggle).

As an example of the above, depending on your inputs, you might run the following:
```shell
$ multi_to_single_fast5 --input_path path/multis --save_path path/singles --recursive
$ tombo preprocess annotate_raws_with_fastqs --fast5-basedir path/singles --fastq-filenames reads.fastq
$ tombo resquiggle path/singles genome.fasta --processes 4 --num-most-common-errors 5
```

### Options
A quick overview of the available options can be seen via `--help`.

```text
Usage: nodclust [OPTIONS] DATASET1 DATASET2 OUTPUT_BED

  Compare two datasets & output an annotated bedMethyl file.

Options:
  -a, --acid TEXT              Filter by acid, dna or rna
  -c, --chromosome TEXT        Filter by chromosome regex
  --force-acid                 Force read files as specified by --acid
  -f, --from-position INTEGER  Filter by minimum position (inclusive)
  -m, --min-coverage INTEGER   Filter by minimum coverage (default 5)
  --no-distance-sum            Don't sum neighbour position distances
  --random-seed INTEGER        Force a random seed, for reproducibility
  -r, --resample-size INTEGER  Signal resample size; 0 to disable (default 15)
  -s, --strand TEXT            Filter by strand, '+' or '-'
  -t, --to-position INTEGER    Filter by maximum position (inclusive)
  --help                       Show this message and exit.
```
