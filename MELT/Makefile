LFLAGS= -L/home/chekkim/LibGlace/versions/netcdf/intelmpi/netcdf-4.7.2_intel21/lib
IFLAGS= -I/home/chekkim/LibGlace/versions/netcdf/intelmpi/netcdf-4.7.2_intel21/include
#LFLAGS=-L/usr/lib/aarch64-linux-gnu
#IFLAGS = -I/usr/include

all: compile

compile:
	elmerf90 DistanceModule.F90 PICO.F90 Melt.F90 -o Melt $(LFLAGS) $(IFLAGS) 
	#elmerf90 IceFrontLocation.F90 -o IceFrontLocation
	#elmerf90 FrontThickness_mask.F90 -o FrontThickness_mask
	#elmerf90 USF_CondFront.F90  -o  USF_CondFront 
	#elmerf90 USF_CondDistance.F90 -o USF_CondDistance
