CC = gcc
#Using -Ofast instead of -O3 might result in faster code, but is supported only by newer GCC versions
CFLAGS = -lm -pthread -O3 -march=native -Wall -funroll-loops -Wno-unused-result

all: word2vec distance

word2vec : word2vec.c
	$(CC) word2vec.c -o word2vec.out $(CFLAGS)
distance : distance.c
	$(CC) distance.c -o distance.out $(CFLAGS)

clean:
	rm -rf word2vec.out distance.out