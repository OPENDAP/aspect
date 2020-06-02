#-- Gfortran

# Added -w to quiet warnings about older fortran constructs. jhrg 5/6/20
FCOMPILE = gfortran -w -c -ffixed-line-length-none -fno-automatic -O3

MAKE = make
# Added -rv for OSX. jhrg 5/6/20
# AR_FLAGS = rvU
AR_FLAGS = -rv

.f.a:
	${FCOMPILE}  $<
	ar $(ARFLAGS) $@ $*.o
	rm -f $*.o
.F.a:
	${FCOMPILE}  $<
	ar $(ARFLAGS) $@ $*.o
	rm -f $*.o

LIB = libSphFunc_opendap
LIBDIR = ..

# Added jhrg 5/6/20
all: $(LIBDIR)/$(LIB).a

$(LIBDIR)/$(LIB).a: \
	$(LIBDIR)/$(LIB).a(ahouse2.o) \
	$(LIBDIR)/$(LIB).a(dot.o) \
	$(LIBDIR)/$(LIB).a(istlen.o) \
	$(LIBDIR)/$(LIB).a(rmod2.o) \
	$(LIBDIR)/$(LIB).a(wsphhead.o) \
	$(LIBDIR)/$(LIB).a(legndr.o) \
	$(LIBDIR)/$(LIB).a(normylm.o) \
	$(LIBDIR)/$(LIB).a(ylm.o)
	ranlib $(LIBDIR)/$(LIB).a

# Added jhrg 5/6/20
clean:
	-rm *.o
	-rm $(LIBDIR)/$(LIB).a
