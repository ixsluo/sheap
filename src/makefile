.DEFAULT:
	-touch $@
all: sheap
sheap.o: sheap.f90 optim.o model.o control.o rng.o constants.o
	$(FC) $(FFLAGS) -c sheap.f90
constants.o: constants.f90
	$(FC) $(FFLAGS) -c constants.f90
rng.o: rng.f90 constants.o
	$(FC) $(FFLAGS) -c rng.f90
random.o: random.f90 rng.o constants.o
	$(FC) $(FFLAGS) -c random.f90
control.o: control.f90 constants.o
	$(FC) $(FFLAGS) -c control.f90
model.o: model.f90 control.o random.o rng.o constants.o
	$(FC) $(FFLAGS) -c model.f90
loss.o: loss.f90 model.o control.o constants.o
	$(FC) $(FFLAGS) -c loss.f90
optim.o: optim.f90 loss.o model.o random.o control.o constants.o
	$(FC) $(FFLAGS) -c optim.f90
SRC = sheap.f90 optim.o loss.f90 model.f90 control.f90 random.f90 rng.f90 constants.f90
OBJ = sheap.o optim.o loss.o model.o control.o random.o rng.o constants.o
clean: neat
	-rm -f sheap
neat:
	-rm -f $(OBJ) *.mod
sheap: $(OBJ) 
	$(LD) $(OBJ) -o sheap $(LDFLAGS)
