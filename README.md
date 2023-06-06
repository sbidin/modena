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
$ python -m nodclust compare path/to/dataset/one and/path/to/dataset/two # Run.
$ python -m nodclust --help # See options.
```

### Install and run via docker
```shell
$ docker build . -t nodclust # Build docker image.
$ docker run -v `pwd`:`pwd` -w `pwd` --rm -it nodclust compare ./dataset/one ./and/two # Run.
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
