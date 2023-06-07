!***********************************************************************!
                                                                        !   
SUBROUTINE VarIntegrator( Model,Solver,dt,TransientSimulation )         !
                                                                        !
!***********************************************************************!

!Declaration
!-----------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!-----------------------------------------------------------------------

  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  TYPE(ParEnv_t), pointer :: MyPE, PEs  

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!  
  character(len=max_name_len), parameter :: defaultVarfile = 'Var_int.dat'
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,VarFile
  CHARACTER(LEN=MAX_NAME_LEN) :: BCName,HSolName,VarSolName
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: TimeVar,HSol
  TYPE(ValueList_t), POINTER :: BC,SolverParams
  TYPE(Nodes_t) :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  REAL(KIND=dp), POINTER :: H(:)
  INTEGER, POINTER :: NodeIndexes(:), HPerm(:)
  Logical :: Firsttime=.TRUE.,Found,Parallel,stat,Gotit
  INTEGER :: i,j,k,l,t,n,NMAX,DIM,ierr
  REAL(kind=dp) :: Vol, Vol_S, time       
  REAL(kind=dp) :: Bu,Bv,u,v,w,s,coeff,SqrtElementMetric,x
  REAL(KIND=dp) :: NodeH(Model % MaxElementNodes),Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)

  save Firsttime,Parallel
  save SolverName,HSolName,VarSolName,VarFile
  save ElementNodes     

!FirsTime initialisation : parallel run and usefull Variables  
!----------------------------------------------------------------------
   If (Firsttime) THEN
    
     ! Check for parallel run 
     Parallel = .FALSE.
     
     IF (ParEnv % PEs > 1)  THEN
       Parallel = .TRUE.
     END IF
  
    SolverParams => GetSolverParams()
    VarSolName = GetString( SolverParams,'Integrated Variable Name', Found)
    
    If(.NOT.Found) THEN
      CALL FATAL(SolverName,'Keyword >Integrated Variable Name< not found in section >Solver<')
    END IF

!End of First visit
  End IF
  Firsttime=.FALSE.


!Get Variables for each time step
!-----------------------------------------------------------------------
  HSol => VariableGet( Solver % Mesh % Variables, VarSolName  )
  IF ( ASSOCIATED( HSol ) ) THEN
        H => HSol % Values
        HPerm => HSol % Perm
  ELSE
        WRITE(Message,'(A,A,A)') 'No variable >',VarSolName,' < found'
        CALL FATAL(SolverName,Message)
  END IF
  
  
! Loop on boundary elements
!-----------------------------------------------------------------------
  
  DO t=1,Solver % Mesh % NumberOfBoundaryElements
    
    Element => GetBoundaryElement(t)
    BC => GetBC()
    IF ( .NOT. ASSOCIATED(BC) ) CYCLE

    !Boundary where to apply integral is to be determined by bed or surface 
    BCName = GetString( BC, 'Name', Found)
    IF ((BCNAME /= 'bed')) CYCLE 
    
    CALL GetElementNodes( ElementNodes )
    n = GetElementNOFNodes()
    NodeIndexes => Element % NodeIndexes
    NodeH(1:n) = H(HPerm(NodeIndexes(1:n)))
    
    ! Numerical integration loop (n Gauss Point by Element)
    IntegStuff = GaussPoints( Element )
    DO i=1,IntegStuff % n
       U = IntegStuff % u(i)
       V = IntegStuff % v(i)
       W = IntegStuff % w(i)

       ! Basis function values
       stat = ElementInfo( Element, ElementNodes, U, V, W, SqrtElementMetric, Basis, dBasisdx )

       x = SUM( ElementNodes % x(1:n) * Basis(1:n) )
       s = 1.0d0

       s = s * SqrtElementMetric * IntegStuff % s(i)

       coeff = SUM(NodeH(1:n) * Basis(1:n))
       Vol=Vol+coeff*s
    
    END DO
    ! End of numerical integration loop
  END DO
  ! End of loop on boundary elements 

!Writing step
!----------------------------------------

  TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )
  Time = TimeVar % Values(1)
  
  Varfile = listgetstring( solver % values,'Integrated Filename',found )
  IF (.NOT. Found) VarFile = DefaultVarFile 
  
  IF (Parallel) THEN
    CALL MPI_ALLREDUCE(Vol, Vol_S, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF (ParEnv % MyPE == 0) THEN
      OPEN (12, file=VarFile, position='append')
      WRITE (12,'(e13.5,2x,e15.8,2x,e15.8,2x,e15.8)') Time, Vol_S
      CLOSE (12)
    END IF
  
  ELSE
    OPEN (12, file=VarFile, position='append')
    WRITE (12,'(e13.5,2x,e15.8,2x,e15.8,2x,e15.8)') Time, Vol
    CLOSE (12)
  END IF

RETURN

!*****************************************************************************!
END SUBROUTINE VarIntegrator                                                  !
!*****************************************************************************!

