FORT    =  mpif90
FORT77  =  mpif90
NICE    = -fdefault-real-8
LIB     =
OPTS    =  -O3 -fbacktrace -ffree-line-length-none
#PROF   = -pg
#DEBUG   =  -g

OBJECTS =  temp.o prms.o data.o main.o si.o out.o timeint.o stats.o nlist.o ran1.o tools.o par.o parallel.o knock.o

mdrun2:  $(OBJECTS) Makefile
	$(FORT) $(LIB) -o mdrun2 $(PROF) $(OBJECTS) $(LIB)

main.o:   main.f90 si.o prms.o data.o timeint.o par.o parallel.o
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  main.f90

stats.o:   stats.f90 si.o prms.o data.o parallel.o
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  stats.f90

temp.o:   temp.f90 prms.o stats.o data.o parallel.o par.o
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  temp.f90

timeint.o:   timeint.f90 si.o prms.o data.o knock.o out.o par.o stats.o temp.o
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  timeint.f90

knock.o:   knock.f90 prms.o data.o
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  knock.f90

si.o:  si.f90 prms.o nlist.o par.o
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  si.f90

nlist.o:   nlist.f90 prms.o data.o
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  nlist.f90

out.o:  out.f90 prms.o stats.o temp.o parallel.o par.o
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  out.f90

data.o:  data.f90 prms.o tools.o
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  data.f90

par.o:  par.f90 prms.o parallel.o
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  par.f90

tools.o:  tools.f90 prms.o parallel.o par.o
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  tools.f90

prms.o:  prms.f90 parallel.o
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  prms.f90

parallel.o:  parallel.f90
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  parallel.f90

ran1.o:  ran1.f
	$(FORT77) -c $(OPTS) $(NICE) ran1.f

installlocal:
	make mdrun2
	mv mdrun2 ../bin

install:
	. /etc/environ.sh; use -e -r openmpi-1.6.3-gnu-4.7.2; make mdrun2
	mv mdrun2 ../bin

ra:
	rm D/*out* D/wav*.plt?????

clean:
	rm -f *.o *.mod

distclean: clean
	rm -f ../bin/mdrun2

