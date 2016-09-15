CC=gcc
GPP=g++
MAKE=make
CFLAGS= -O3
LIBS= -lm

all:    pfilt seq2mtx seqfeat

pfilt:
	$(CC) $(CFLAGS) -o pfilt pfilt.c -lm

seq2mtx:
	$(CC) $(CFLAGS) -o seq2mtx seq2mtx.c -lm

seqfeat:
	$(GPP) features.cpp -o features

clean:
	-rm -f ./pfilt
	-rm -f ./seq2mtx
	-rm -f ./features
	-rm -f ../bin/pfilt
	-rm -f ../bin/seq2mtx
	-rm -f ../bin/features

install:
	mv pfilt ../bin
	mv seq2mtx ../bin
	mv features ../bin
