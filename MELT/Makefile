LFLAGS= -L/home/chekkim/LibGlace/versions/netcdf/intelmpi/netcdf-4.7.2_intel21/lib
IFLAGS= -I/home/chekkim/LibGlace/versions/netcdf/intelmpi/netcdf-4.7.2_intel21/include
#LFLAGS=-L/usr/lib/aarch64-linux-gnu
#IFLAGS = -I/usr/include

all: compile


SRC=./MELT/
BLD='.'

# MELT (including the PICO model)
melt: $(SRC)/melt.F90 $(SRC)/DistanceModule.F90
	elmerf90 -check all -debug $(NCINC) $^ -o $(BLD)/$@ $(FFLAGS) $(NCLIB) -f

compile:
	elmerf90 DistanceModule.F90 PICO.F90 Melt.F90 -o Melt $(LFLAGS) $(IFLAGS) 

