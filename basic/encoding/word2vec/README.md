# Text encoding using wordvec

## Quick-run text encoding

1. Compile and quickly generate vector binaries on the small `canadianLegal.txt` dataset. Feel free to swap with your own textual dataset. The bigger the better results!

```console
$ make
$ ./generateVectorbin.sh
```

2. Interactive console to view vector distances of generated `vectors.bin` file

```console
$ ./distance.out vectors.bin
```

## References

[Word2Vec original archive](https://code.google.com/archive/p/word2vec/)

