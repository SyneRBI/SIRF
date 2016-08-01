all:
	@if [ ! -r testdata.h5 ]; then \
	ismrmrd_generate_cartesian_shepp_logan -r 2; \
	fi;
	@if [ ! -r testdata_a2.h5 ]; then \
	ismrmrd_generate_cartesian_shepp_logan -a 2; \
	fi;

