#-- Gfortran

# Added -w to quiet warnings about older fortran constructs. jhrg 5/6/20
FCOMPILE = gfortran -w -c -ffixed-line-length-none -static -fno-automatic -mcmodel=medium -m64 -O3
# FFCOMPILE = gfortran -c  -ffixed-line-length-none         -fno-automatic -O3

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

LIB = libHJVH
LIBDIR = ..

# Added jhrg 5/6/20
all: $(LIBDIR)/$(LIB).a

$(LIBDIR)/$(LIB).a: \
	$(LIBDIR)/$(LIB).a(addata.o) \
	$(LIBDIR)/$(LIB).a(ahouse.o) \
	$(LIBDIR)/$(LIB).a(ahouse2.o) \
	$(LIBDIR)/$(LIB).a(convpack.o) \
	$(LIBDIR)/$(LIB).a(dcopy.o) \
	$(LIBDIR)/$(LIB).a(scopy.o) \
	$(LIBDIR)/$(LIB).a(dot.o) \
	$(LIBDIR)/$(LIB).a(getfdp.o) \
	$(LIBDIR)/$(LIB).a(integ.o) \
	$(LIBDIR)/$(LIB).a(istlen.o) \
	$(LIBDIR)/$(LIB).a(inittrf.o) \
	$(LIBDIR)/$(LIB).a(modtrf.o) \
	$(LIBDIR)/$(LIB).a(rmod.o) \
	$(LIBDIR)/$(LIB).a(rmod2.o) \
	$(LIBDIR)/$(LIB).a(saxpy.o) \
	$(LIBDIR)/$(LIB).a(sdot.o) \
	$(LIBDIR)/$(LIB).a(splh.o) \
	$(LIBDIR)/$(LIB).a(strgrep.o) \
	$(LIBDIR)/$(LIB).a(rsphhead.o) \
	$(LIBDIR)/$(LIB).a(wsphhead.o) \
	$(LIBDIR)/$(LIB).a(wspthead.o) \
	$(LIBDIR)/$(LIB).a(wint2ch.o)

	ranlib $(LIBDIR)/$(LIB).a

# Added jhrg 5/6/20
clean:
	-rm *.o
	-rm $(LIBDIR)/$(LIB).a
