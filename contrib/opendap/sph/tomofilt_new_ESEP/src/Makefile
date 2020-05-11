
# Added -w to quiet warnings about older fortran constructs. jhrg 5/6/20
FF        = gfortran -w -ffixed-line-length-none -fno-automatic -O3
CC        = gcc -O3

BIN       = $(TOMOFILT)/bin
LIB       = $(TOMOFILT)/lib

all: mkexpmatxy invexpandxy sphexp sphadd

# The original Makefile built all of the binaries, but opendap's script
# only needs some of them. jhrg 5/7/20
original: mkexpmatxy invexpandxy mk3d_res mk3d_res_ESEP sphadd spt2sph \
	sphexp linegrep depmap raw2xyz-eqdist sph2v_input pwrspecsph correlatorsph

mkexpmatxy:
	${FF} -o $(BIN)/mkexpmatxy mkexpmatxy.f \
	$(LIB)/libHJVH.a \
	$(LIB)/libSphFunc.a \
	$(LIB)/libgeo.a

invexpandxy:
	${FF} -o $(BIN)/invexpandxy invexpandxy.f \
	$(LIB)/libHJVH.a \
	$(LIB)/libSphFunc.a \
	$(LIB)/libgeo.a

sphadd:
	${FF} -o $(BIN)/sphadd sphadd.f \
	$(LIB)/libHJVH.a \
	$(LIB)/libSphFunc.a \
	$(LIB)/libgeo.a

sphexp:
	${FF} -o $(BIN)/sphexp sphexp.f \
	$(LIB)/libHJVH.a \
	$(LIB)/libSphFunc.a \
	$(LIB)/libgeo.a

# Below this are the 'original' executables. jhrg 5/7/20
mk3d_res:
	${FF} -o $(BIN)/mk3d_res mk3d_res.f \
	$(LIB)/libHJVH.a \
	$(LIB)/libSphFunc.a \
	$(LIB)/libgeo.a

mk3d_res_ESEP:
	${FF} -o $(BIN)/mk3d_res_ESEP mk3d_res_ESEP.f \
	$(LIB)/libHJVH.a \
	$(LIB)/libSphFunc.a \
	$(LIB)/libgeo.a

spt2sph:
	${FF} -o $(BIN)/spt2sph spt2sph.f \
	$(LIB)/libHJVH.a \
	$(LIB)/libSphFunc.a \
	$(LIB)/libgeo.a

linegrep:
	${CC} -o $(BIN)/linegrep linegrep.c

depmap:
	${FF} -o $(BIN)/depmap depmap.f ${LIB}/libS20.a

raw2xyz-eqdist: 
	${FF} -o $(BIN)/raw2xyz-eqdist raw2xyz-eqdist.f

sph2v_input:
	${FF} -o $(BIN)/sph2v_input sph2v_input.f

pwrspecsph:
	${FF} -o $(BIN)/pwrspecsph pwrspecsph.f

correlatorsph:
	${FF} -o $(BIN)/correlatorsph correlatorsph.f

