
# New Makefile to build/install the TOMOFILT ESEP code.

all: libs srcs

libs:
	(cd lib && $(MAKE) $(MFLAGS))

srcs:
	(cd src && $(MAKE) $(MFLAGS))

clean:
	-rm *.sph
	-rm *.spt
	(cd geodyn && $(MAKE) $(MFLAGS) clean)

libclean: clean
	(cd lib && $(MAKE) $(MFLAGS) clean)

execs = correlatorsph depmap invexpandxy linegrep mk3d_res mk3d_res_ESEP \
	mkexpmatxy pwrspecsph raw2xyz-eqdist sph2v_input sphadd sphexp spt2sph

reallyclean: libclean
	-(cd bin && rm $(execs))
