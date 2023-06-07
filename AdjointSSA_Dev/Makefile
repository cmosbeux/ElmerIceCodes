
EXECUTABLES := AdjointSSASolvers

SRC=src
OBJ=objects

_SOURCES := AdjointSolverSSA.F90 \
	DJDp_Adjoint_SSA.F90 \
	CostSolver_SSA_Adjoint.F90 \
	CostSolver_Disc.F90 \
	CostSolver_Disc_Reg.F90 \
	CostSolver_SSA_Regularisation.F90 \
	CostSolver_SSA_FluxDiv_Adjoint.F90 \
	SSASolver.F90

SOURCES = $(patsubst $(SRC)/%, % ,$(_SOURCES))


_SOLVERS := $(patsubst %.F90, %.o ,$(SOURCES))
SOLVERS = $(patsubst %, $(OBJ)/% ,$(_SOLVERS))

compile: $(EXECUTABLES)

AdjointSSASolvers: $(SOLVERS)
	elmerf90 $^ -o $@

$(SOLVERS): $(OBJ)/%.o: $(SRC)/%.F90
	@if [ ! -d "./$(OBJ)" ]; then \
	     mkdir $(OBJ) ;\
	fi
	elmerf90 -c $< -o $@


clean: 
	rm -f $(EXECUTABLES)
	rm -f $(SOLVERS)
