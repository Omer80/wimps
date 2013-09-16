.PHONY: clean init

ifeq ($(MAKECMDGOALS),init)
init:
#	$(MAKE) -C sources/micromegas
	$(MAKE) -C sources
#	echo \#define LPATH \"`pwd`\" > microRun/lib/lpath.h
#	$(MAKE) -C microRun  main=main.c
else 

ifeq ($(MAKECMDGOALS),clean)
clean:
#	$(MAKE) -C sources/micromegas clean
	$(MAKE) -C sources clean
	$(MAKE) -C main clean
#	rm  -f microRun/lib/lpath.h	
#	$(MAKE) -C microRun clean 
else

micrO = sources/micromegas

all:
	$(MAKE) -C main  nmhdecay nmspec
endif

endif
