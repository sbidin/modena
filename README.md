# signals

TODO a more detailed description. Ivan neka ovdje uskoči.

### Install and run via docker
```shell
$ docker build . -t signals # Build docker image.
$ docker run -v `pwd`:`pwd` -w `pwd` --rm -it signals compare ./dataset/one ./and/two # Run.
$ docker run --rm -it signals --help # See options.
```

### Install and run locally
```shell
$ python -m pip install -e . # Install.
$ python -m signals compare path/to/dataset/one and/path/to/dataset/two # Run.
$ python -m signals --help # See options.
```

### Supported input formats
Only single-FAST5 files processed by
[tombo](https://nanoporetech.github.io/tombo/index.html) are supported. If your
datasets haven't been processed by tombo already, follow these steps:

* If your input is in multi-FAST55 format, convert it to single-FAST5 by running
  the `multi_to_single_fast5` utility [available
  here](https://github.com/nanoporetech/ont_fast5_api).
* If your inputs don't contain basecalls, annotate them by running tombo's
  `preprocess annotate_raws_with_fastqs` subcommand [as shown
  here](https://nanoporetech.github.io/tombo/examples.html?highlight=annotate_raw_with_fastqs).
* Run your inputs through tombo's `resquiggle` command [as explained
  here](https://nanoporetech.github.io/tombo/examples.html?highlight=resquiggle).

As an example of the above, after installing tombo, run the following:
```shell
$ multi_to_single_fast5 --input_path path/multis --save_path path/singles --recursive
$ tombo preprocess annotate_raws_with_fastqs --fast5-basedir path/singles --fastq-filenames reads.fastq
$ tombo resquiggle path/singles genome.fasta --processes 4 --num-most-common-errors 5
```
