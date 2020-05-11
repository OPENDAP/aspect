
# These are the libraries we need for the 'opendap' version of the
# dofilt _ESEP script. Below the 'oritignal' target wlll build the
# original set of libraries needed by the original binaries
SUBDIRS = libgeo libHJVH libSphFunc libS20

all:
	@for m in $(SUBDIRS); do \
		(cd $$m && $(MAKE) $(MFLAGS) all); \
	done

original:
	@for m in $(SUBDIRS); do \
		(cd $$m && $(MAKE) $(MFLAGS) all); \
	done

clean:
	@for m in $(SUBDIRS); do \
		(echo $$m && cd $$m && $(MAKE) $(MFLAGS) clean); \
	done
