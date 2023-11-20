#FC=ifort
#FFLAGS= -O3 -fopenmp -no-wrap-margin -assume buffered_io

FC=gfortran
FFLAGS=-O3 -fopenmp

#FFLAGS=-O3 -g -fopenmp -Wall -fwhole-file -fcheck=all -pedantic -fbacktrace -fall-intrinsics -ffpe-trap=overflow -fbounds-check

#FFLAGS= -O0 -g -fopenmp -fimplicit-none -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wunused-parameter -fwhole-file -fcheck=all -std=f2008 -pedantic -fbacktrace -fall-intrinsics -ffpe-trap=invalid,zero,overflow -fbounds-check -Wuninitialized -Wfatal-errors

# Define the version of spglib to be used
export SPG_VERSION=1.14.1

#LDFLAGS=$(FFLAGS) -mkl
#LDFLAGS=$(FFLAGS) -L$(PREFIX)/external/spglib/lib -mkl -lsymspg
#LDFLAGS=$(FFLAGS) -llapack -lopenblas
LDFLAGS=$(FFLAGS) -L$(PREFIX)/external/spglib/lib -llapack -lopenblas -lsymspg
LD=$(FC)
PREFIX=$(PWD)

export

all: sheap external

external: spglib cabal

sheap:
	(cd src/; make)

spglib:
	(cd external/spglib; make)

cabal:
	(cd external/cabal; make)

install:
	(cp src/sheap bin/)
	(cp external/cabal/cabal-sheap bin/)

neat:
	(cd src; make clean)
	(cd external/spglib; make clean)
	(cd external/cabal; make clean)

clean: neat
	(rm -f bin/sheap)

dist: clean
	tar -czf ../sheap-`date "+%d%m%Y"`.tgz  --exclude=".*" -C .. sheap
