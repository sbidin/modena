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
$ python -m signals path/to/dataset/one path/to/dataset/two # Run.
$ python -m signals --help # See more options here.
```
