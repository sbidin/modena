# signals

TODO a more detailed description.


### install and run via docker
```shell
$ docker build . -t signals # Build docker image.
$ docker run --rm -it signals --help # See options.
$ docker run -v `pwd`:`pwd` -w `pwd` --rm -it signals ./dataset/one ./and/two # Run.
```

### install and run locally
```shell
$ python -m pip install -e . # Install.
$ python -m signals path/to/dataset/one and/path/to/dataset/two # Run.
$ python -m signals --help # See more options here.
```

### inputs must be re-squiggled
Only re-squiggled FAST5 datasets are supported. If your dataset hasn't been
resquiggled, install
[tombo](https://github.com/nanoporetech/tombo#getting-started) and apply [the
tombo re-squiggle
algorithm](https://nanoporetech.github.io/tombo/resquiggle.html) on it first.
```shell
$ tombo resquiggle path/to/fast5s/ genome.fasta --processes 4 --num-most-common-errors 5
```
