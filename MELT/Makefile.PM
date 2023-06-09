
SRC=./

BLD=../MY_BLD/

NCINC=-I$(NETCDFFORTRAN_INCDIR) -I$(NETCDFC_INCDIR)

NCLIB=-L$(NETCDFFORTRAN_LIBDIR) -L$(NETCDFC_LIBDIR) -lnetcdf -lnetcdff

VPATH=$(BLD)

FFLAGS=-I$(XIOS_INC) -O2 -fpe0 -ftrapuv -g -traceback

LIBS= USF_ExitCondition boxmodel_solver_t FrontThickness_mask USF_CondFront USF_CondDistance SlipCoef nearestpoint GetConnectedAreas passive_cond ExcludeAreas UGridDataReader

all: $(LIBS)

#---------------------------------------------------
# USR function
USF_ExitCondition:${SRC}/USF_ExitCondition.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

UGridDataReader:$(SRC)/UGridDataReader.F90
	elmerf90 $(NCINC) $^ -o $(BLD)/$@ $(FFLAGS) $(NCLIB)

#---------------------------------------------------
# Post Porcessing
# PP define basin
nearestpoint:$(SRC)/nearestpoint.F90
	elmerf90 $(NCINC) $^ -o $(BLD)/$@ $(FFLAGS) $(NCLIB)
#
#---------------------------------------------------

# BMB 
# BMB: PICO
boxmodel_solver_t:$(SRC)/boxmodel_solver_t.F90
	elmerf90 $(NCINC) $^ -o $(BLD)/$@ $(FFLAGS) $(NCLIB)

FrontThickness_mask:$(SRC)/FrontThickness_mask.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

USF_CondFront:$(SRC)/USF_CondFront.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

USF_CondDistance:$(SRC)/USF_CondDistance.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

#---------------------------------------------------
# Slip Coefficient correction Coulomb
SlipCoef:$(SRC)/SlipCoef.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

clean:
	cd $(BLD) && rm -f $(LIBS)

#---------------------------------------------------
# Processing icebergs
GetConnectedAreas: $(SRC)/GetConnectedAreas.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

passive_cond: $(SRC)/passive_cond.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

ExcludeAreas: $(SRC)/ExcludeAreas.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

clean:
	cd $(BLD) && rm -f $(LIBS)