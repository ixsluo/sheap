.DEFAULT:

all: get libsymspg.a

SPGLIB=spglib-$(SPG_VERSION)

# The naming of the archive is like v1.14.1.zip or master.zip
ifeq ($(SPG_VERSION), master)
  SPG_ARCHIVE=$(SPG_VERSION).zip
else
  SPG_ARCHIVE=v$(SPG_VERSION).zip
endif

OBJ = $(SPGLIB)/src/arithmetic.o $(SPGLIB)/src/cell.o $(SPGLIB)/src/debug.o $(SPGLIB)/src/delaunay.o $(SPGLIB)/src/hall_symbol.o $(SPGLIB)/src/kgrid.o \
    $(SPGLIB)/src/kpoint.o $(SPGLIB)/src/mathfunc.o $(SPGLIB)/src/niggli.o $(SPGLIB)/src/overlap.o $(SPGLIB)/src/pointgroup.o $(SPGLIB)/src/primitive.o \
    $(SPGLIB)/src/refinement.o $(SPGLIB)/src/site_symmetry.o $(SPGLIB)/src/sitesym_database.o $(SPGLIB)/src/spacegroup.o \
    $(SPGLIB)/src/spg_database.o $(SPGLIB)/src/spglib.o $(SPGLIB)/src/spin.o $(SPGLIB)/src/symmetry.o $(SPGLIB)/src/determination.o 

get:
	-wget -c --timeout=10 --tries=1 -nv https://github.com/atztogo/spglib/archive/v$(SPG_VERSION).zip
	-unzip -qq -u $(SPG_ARCHIVE) '$(SPGLIB)/src/*' '$(SPGLIB)/example/*'

libsymspg.a: $(OBJ)
	-ar -cr lib/libsymspg.a $(OBJ)
	-mkdir -p include
	-cp $(SPGLIB)/src/spglib.h include/spglib.h

# Explicitly clean master.zip
clean:
	-rm -fr $(SPGLIB)
	-rm -f master.zip
