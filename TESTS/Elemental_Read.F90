SUBROUTINE ElemRead( Model,Solver,dt,Transient )
  !------------------------------------------------------------------------------
  !USE CoordinateSystems
  !USE MeshUtils
  USE Netcdf
  USE DefUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt

  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  TYPE(Mesh_t),POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Solver_t),POINTER :: PSolver
  TYPE(Nodes_t) :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  TYPE(Element_t),POINTER ::  Element
 
  TYPE(Variable_t),POINTER, SAVE :: BasinVar=>NULL()
  INTEGER, POINTER, SAVE :: BasinPerm(:)
  REAL(KIND=dp), POINTER, SAVE :: Basin(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='Test'


  !------------------------------------------------------------------------------
  ! 1- Read constants and parameters of the simulation :
  !------------------------------------------------------------------------------

  !- Simulation parameters (idealized ocean forcing) :

  Params => GetSolverParams()

  Mesh => Model % Mesh

  Nmax = Solver % Mesh % NumberOfNodes

  ! get the number of the basin to associate values
  BasinVar => VariableGet( Model % Mesh % Variables, 'basins')
  IF (.NOT.ASSOCIATED(BasinVar)) &
       &    CALL FATAL(SolverName,'basins not found')

  IF ( BasinVar % TYPE .NE. Variable_on_element) &
       &   CALL FATAL(SolverName,'basins not a vairable on element')

  BasinPerm => BasinVar % Perm
  Basin => BasinVar % Values

CALL INFO(Trim(SolverName),'READ VARIABLE', Level =5)

DO e=1,Solver % NumberOfActiveElements
    Element => GetActiveElement(e)
    CALL GetElementNodes( ElementNodes, Element, Solver)
    n = GetElementNOFNodes()
    NodeIndexes => Element % NodeIndexes
    Indexx => Element % ElementIndex

    b = NINT(MAXVAL(Basin(BasinPerm(Indexx))))
    write(*,*) b

END SUBROUTINE ElemRead




