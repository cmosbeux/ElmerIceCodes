LFLAGS= -L/home/chekkim/LibGlace/versions/netcdf/intelmpi/netcdf-4.7.2_intel18/lib
IFLAGS= -I/home/chekkim/LibGlace/versions/netcdf/intelmpi/netcdf-4.7.2_intel18/include

all: compile

compile:  
        elmerf90 boxmodel_solver.F90 -o boxmodel_solver $(LFLAGS) $(IFLAGS)
        elmerf90 FrontThickness_mask.F90 -o FrontThickness_mask
        elmerf90 USF_CondFront.F90  -o  USF_CondFront 
        elmerf90 USF_CondDistance.F90 -o USF_CondDistance

