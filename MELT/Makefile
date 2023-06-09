#LFLAGS= -L/home/chekkim/LibGlace/versions/netcdf/intelmpi/netcdf-4.7.2_intel21/lib
#IFLAGS= -I/home/chekkim/LibGlace/versions/netcdf/intelmpi/netcdf-4.7.2_intel21/include
LFLAGS=-L/usr/lib/aarch64-linux-gnu
IFLAGS = -I/usr/include

all: compile

SRC=.
BLD=.

# MELT (including the PICO model)
melt: $(SRC)/PICO.F90 $(SRC)/DistanceModule.F90 $(SRC)/Melt.F90
	elmerf90 $(NCINC) $^ -o $(BLD)/$@ $(FFLAGS) $(NCLIB)


melt2:
	elmerf90 DistanceModule.F90 PICO.F90 Melt.F90 -o Melt $(LFLAGS) $(IFLAGS)
