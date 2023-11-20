.DEFAULT:
	-touch $@
all: cabal
cabal.o: cabal.f90 spglib_f08.o niggli_reduce.o constants.o
	$(FC) $(FFLAGS) -c cabal.f90
spglib_f08.o: spglib_f08.f90
	$(FC) $(FFLAGS) -c ../spglib/spglib-$(SPG_VERSION)/example/spglib_f08.f90
niggli_reduce.o: niggli_reduce.f90 constants.o
	$(FC) $(FFLAGS) -c niggli_reduce.f90
constants.o: $(PREFIX)/src/constants.f90
constants.o: $(PREFIX)/src/constants.f90
	$(FC) $(FFLAGS) -c $(PREFIX)/src/constants.f90

SRC = cabal.f90 ../spglib/spglib-master/example/spglib_f08.f90 niggli_reduce.f90 $(PREFIX)/src/constants.f90
OBJ = cabal.o spglib_f08.o niggli_reduce.o constants.o
clean: neat
	-rm -f cabal-sheap
neat:
	-rm -f $(OBJ) *.mod
cabal: $(OBJ) 
	$(LD) $(OBJ) -o cabal-sheap $(LDFLAGS)
