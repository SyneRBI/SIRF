# Copyright Nov 2016 Casper da Costa-Luis
# Apache License 2.0

# This file builds to ./build and installs to ./build/dist[bin|lib|include]
.PHONY:
	all
	clean
	distclean

all: build
	$(MAKE) -C build install

build:
	mkdir -p $@ && cd $@ && cmake -DCMAKE_INSTALL_PREFIX=dist ..

clean: build
	$(MAKE) -C build clean

distclean:
	rm -rf build

DC=docker
XA=xargs -r

dbuild:
	$(DC) build --build-arg mainUser=sirf -t casperdcl/sirf .

drun:
	$(DC) run -it --rm --name sirf_container casperdcl/sirf

dclean:
	$(DC) ps -a -q | $(XA) $(DC) rm
	$(DC) volume ls -q -f "dangling=true" | $(XA) $(DC) volume rm
	$(DC) images -q -f "dangling=true" | $(XA) $(DC) rmi
