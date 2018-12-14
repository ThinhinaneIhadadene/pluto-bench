#include directories of benchmarks
DIRS = cvtColor \
	   edgeDetect \
	   gaussian \
		 convolution\
		 heat3d\
#		 warp_affine

all: orig par dist distopt distopt_fop distopt_foifi

orig:
	@-for d in $(DIRS); do \
		make -C $$d $@; \
		done

par:
	@-for d in $(DIRS); do \
		make -C $$d par; \
		done

dist:
	@-for d in $(DIRS); do \
		make -C $$d dist; \
		done

distopt:
	@-for d in $(DIRS); do \
	make -C $$d distopt; \
	done

distopt_fop:
	@-for d in $(DIRS); do \
	make -C $$d distopt_fop; \
	done

distopt_foifi:
	@-for d in $(DIRS); do \
	make -C $$d distopt_foifi; \
	done

clean:
	@-for d in $(DIRS); do \
		make -C $$d clean; \
		done
