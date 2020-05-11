
# Clean these directories to see what is generated and to keep
# generated data out of github. jhrg 5/7/20

clean: clean-dat clean-raw
	-rm *.sph
	-rm *.spt

DATA_DIRS = outputfiles comparefiles slices pwrfiles

clean-dat:
	-for d in $(DATA_DIRS); do (cd $$d && rm *.dat); done

clean-raw:
	-(cd rawfiles && rm *.raw)
