
# Added -w to quiet warnings about older fortran constructs. jhrg 5/6/20
FF        = gfortran -w -ffixed-line-length-none -fno-automatic -O3
CC        = gcc -O3

BIN       = $(TOMOFILT)/bin
LIB       = $(TOMOFILT)/lib

all: check_tomofilt opendap_convert mkexpmatxy invexpandxy sphexp sphadd

check_tomofilt:
	@if test -n "$$TOMOFILT"; then true; else echo "Define TOMOFILT"; false; fi

# The original Makefile built all of the binaries, but opendap's script
# only needs some of them. jhrg 5/7/20
original: mkexpmatxy invexpandxy mk3d_res mk3d_res_ESEP sphadd spt2sph \
	sphexp linegrep depmap raw2xyz-eqdist sph2v_input pwrspecsph correlatorsph

opendap_convert:
	${FF} -o $(BIN)/opendap_convert opendap_convert.f $(LIB)/libSphFunc_opendap.a

mkexpmatxy:
	${FF} -o $(BIN)/mkexpmatxy mkexpmatxy.f $(LIB)/libSphFunc_opendap.a

invexpandxy:
	${FF} -o $(BIN)/invexpandxy invexpandxy.f $(LIB)/libSphFunc_opendap.a

sphadd:
	${FF} -o $(BIN)/sphadd sphadd.f $(LIB)/libSphFunc_opendap.a

sphexp:
	${FF} -o $(BIN)/sphexp_opendap sphexp_opendap.f \
	$(LIB)/libSphFunc_opendap.a

