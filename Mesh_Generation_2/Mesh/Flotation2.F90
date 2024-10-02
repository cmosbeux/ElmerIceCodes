!!! Apply the flotation criterion to update Zb and Zs from the ice thickness H
!
!   OUTPUT Variables:
!     groundedmask
!     Zb
!     Zs
!   
!   INPUT Variable:
!     H
!     bedrock 
!
! PARAMETERS:
!   Constants: 
!     zsea
!     rhow
!   Material:
!      SSA Mean Density
!     
!------------------------------------------------------------------------------
SUBROUTINE Flotation2( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE CoordinateSystems
  USE MeshUtils
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
  TYPE(Variable_t),POINTER :: Var
  TYPE(Variable_t),POINTER :: ZbVar,ZsVar
  TYPE(Variable_t),POINTER :: HVar,BedVar
  TYPE(Variable_t),POINTER :: GMVar
  TYPE(Element_t), POINTER :: Element,Parent
  TYPE(ValueList_t),POINTER :: BodyForce,Material, Params

  REAL(KIND=dp),DIMENSION(:),ALLOCATABLE,SAVE :: Density
  REAL(KIND=dp) :: zsea,rhow,rhoi
  REAL(KIND=dp) :: H,zb,zs,bedrock
  REAL(KIND=dp),PARAMETER :: EPS=0.1_dp


  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER :: Unground,UngroundTot
  INTEGER :: t,i,j,n,np
  INTEGER :: ierr

  LOGICAL,SAVE :: Initialized = .FALSE.
  LOGICAL :: Found
  LOGICAL :: GotIt

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='Flotation'
  CHARACTER(LEN=MAX_NAME_LEN) :: ZbName,ZsName,HName

#ifndef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: RealTime,CPUTime
#endif
  REAL(KIND=dp) :: at,st,at0,st0
  INTEGER :: ndow
  LOGICAL :: SolverTiming=.TRUE.

!----------------------------------------------------------------------------
  IF (SolverTiming) THEN
     at0=RealTime()
     st0=CPUTime()
  ENDIF

!------------------------------------------------------------------------------
  Mesh => Model % Mesh

  Params => Solver % Values

!!! get required variables Zb,Zs,H
  ZbName = GetString(Params, 'Bottom Surface Name', GotIt)
  IF (GotIt) THEN
    CALL INFO(SolverName, 'Bottom Surface Name found', level=4)
  ELSE
    CALL INFO(SolverName, 'Bottom Surface Name not found - using default Zb', level=1)
    WRITE(ZbName,'(A)') 'Zb'
  END IF
  zbVar => VariableGet( Model % Mesh % Variables, ZbName,UnFoundFatal=.TRUE.)
  
  ZsName = GetString(Params, 'Top Surface Name', GotIt)
  IF (GotIt) THEN
    CALL INFO(SolverName, 'Top Surface Name found', level=4)
  ELSE
    CALL INFO(SolverName, 'Top Surface Name not found - using default Zs', level=1)
    WRITE(ZsName,'(A)') 'Zs'
  END IF
  zsVar => VariableGet( Model % Mesh % Variables, ZsName,UnFoundFatal=.TRUE.)

  HName = GetString(Params, 'Thickness Variable Name', GotIt)
  IF (GotIt) THEN
    CALL INFO(SolverName, 'Thickness Variable Name found', level=4)
  ELSE
    CALL INFO(SolverName, 'Thickness Variable  Name not found - using default H', level=1)
    WRITE(HName,'(A)') 'H'
  END IF
  HVar => VariableGet( Model % Mesh % Variables, HName, UnFoundFatal=.TRUE.)

!!
!! get variables GLMAsk,bedrock
  GMVar => VariableGet( Model % Mesh % Variables, 'GroundedMask',UnFoundFatal=.TRUE.)
  IF ((ParEnv % PEs>1).AND.(.NOT.ASSOCIATED(Solver%Matrix))) &
       CALL FATAL(SolverName,'Solver%Matrix should be associated to update GLMask')

  BedVar => VariableGet( Model % Mesh % Variables, 'bedrock',UnFoundFatal=.TRUE.)

!!! Do some initialisation/allocation
  IF ((.NOT.Initialized).OR.Mesh%Changed) THEN

    IF (Initialized) deallocate(Density)
    
    N=Model % MaxElementNodes
    allocate(Density(N))

    IF (ListGetLogical(Params,'Grounded Mask Re-init')) &
     GMVar % Values = +1.0_dp

    Initialized = .TRUE.
  END IF
!!

 zsea = GetCReal( Model % Constants, 'Sea Level', Found )
 If (.NOT.Found) THEN
    WRITE(Message,'(A)') 'Constant >Sea Level< not found. &
       &Setting to 0.0'
    CALL INFO(SolverName, Message, level=3)
    zsea = 0._dp
 End if
 rhow = GetCReal( Model % Constants, 'water density', Found )
 If (.NOT.Found) THEN
    WRITE(Message,'(A)') 'Constant Water Density not found. &
       &Setting to 1.03225e-18'
    CALL INFO(SolverName, Message, level=3)
    rhow = 1.03225d-18
 END IF

!! GL is grounded
! WHERE(abs(GMVar % Values).LT.EPS) GMVar % Values = +1.0_dp
  !$omp parallel do
  Do t=1,Solver % Mesh % NumberOfNodes
     IF (GMVar % Perm (t) == 0) CYCLE
     IF (GMVar % Values(GMVar % Perm (t)).GE.0._dp) THEN
        GMVar % Values(GMVar % Perm (t))=+1._dp
     ELSE
        GMVar % Values(GMVar % Perm (t))=-1._dp
     END IF
  end do
  !$omp end parallel do


!$omp parallel do private(Element,n,NodeIndexes,Parent,np,Material,Density,i,bedrock,j,H,rhoi,zb)
 !! Boundary nodes are grounded if bedrock above sea level or (GL or floating)
 Do t=1,Model % NumberOfBoundaryElements

   Element => GetBoundaryElement(t)
   n = GetElementNOFNodes(Element)
   NodeIndexes => Element % NodeIndexes

   Parent => Element % BoundaryInfo % Left
   IF (.NOT.ASSOCIATED(Parent)) Parent => Element % BoundaryInfo % Right
   IF (.NOT.ASSOCIATED(Parent)) CALL FATAL(SolverName,'Unable to found parent element')
   np = Parent % TYPE % NumberOfNodes

   Material => GetMaterial(Parent)
   Density(1:np) = ListGetReal( Material, 'SSA Mean Density',np,Parent% NodeIndexes, UnfoundFatal=.TRUE.)

   Do i=1,n
     bedrock=BedVar%Values(BedVar%Perm(NodeIndexes(i)))
     IF (bedrock.GT.zsea) THEN
       GMVar % Values(GMVar % Perm(NodeIndexes(i)))=1._dp
     ELSE
       Do j=1,np
          IF (NodeIndexes(i).EQ.Parent% NodeIndexes(j)) EXIT
       End Do
       IF (NodeIndexes(i).NE.Parent% NodeIndexes(j)) &
          CALL FATAL(SolverName,'Node Index not found in Parent')
         
       H=HVar%Values(HVar%Perm(NodeIndexes(i)))
       rhoi=Density(j)
       zb=zsea-H*rhoi/rhow
       IF (zb.LE.bedrock) THEN
         GMVar % Values(GMVar % Perm(NodeIndexes(i)))=0._dp
       ELSE
         GMVar % Values(GMVar % Perm(NodeIndexes(i)))=-1._dp
       END IF
     ENDIF
   End DO
 End do
 !$omp end parallel do

 
!$omp parallel do private(Element,n,NodeIndexes,Material,Density,i,H,bedrock,rhoi,zb,zs)
 ! check if floating element needs to be grounded (GL can advance)
 Do t=1,Solver % Mesh % NumberOfBulkElements
   Element => Solver % Mesh % Elements(t)
   n = GetElementNOFNodes(Element)
   NodeIndexes => Element % NodeIndexes

   Material => GetMaterial(Element)
   Density(1:n) = &
    ListGetReal( Material, 'SSA Mean Density',n, NodeIndexes,UnFoundFatal=.True.)

   Do i=1,n
     H=HVar%Values(HVar%Perm(NodeIndexes(i)))
     bedrock=BedVar%Values(BedVar%Perm(NodeIndexes(i)))
     ! Grounded node zb=bedrock
     IF (GMVar % Values(GMVar % Perm(NodeIndexes(i))).GT.-EPS) THEN
       zb = bedrock
     ELSE !evaluate floatation to see if the node has to be grounded
      rhoi=Density(i)
      zb=zsea-H*rhoi/rhow
      IF (zb.LE.bedrock) THEN
          GMVar % Values(GMVar % Perm(NodeIndexes(i)))=+1._dp
          zb=bedrock
       END IF
     END IF
     ! Update zb and zs variables
     zs=zb+H
     ZbVar%Values(ZbVar%Perm(NodeIndexes(i)))=zb
     ZsVar%Values(ZsVar%Perm(NodeIndexes(i)))=zs
   End Do
 End DO
!$omp end parallel do
 
 ndow = 0
 DO WHILE (.TRUE.) 
  ndow = ndow + 1

  !$omp parallel do private(Element,n,NodeIndexes,i)
   ! GL is the boundary between floating and grounded Elements
   Do t=1,Solver % Mesh % NumberOfBulkElements
      Element => Solver % Mesh % Elements(t)
      n = GetElementNOFNodes(Element)
      NodeIndexes => Element % NodeIndexes

     IF (ANY(GMVar % Values(GMVar % Perm(NodeIndexes(1:n))).LT.-EPS)) THEN
        Do i=1,n
          IF (GMVar % Values(GMVar % Perm(NodeIndexes(i))).GT.EPS) &
             GMVar % Values(GMVar % Perm(NodeIndexes(i)))=0._dp
        End Do
     END IF
   End do
  !$omp end parallel do

   ! Get Parallel Min
   IF ( ParEnv % PEs>1 ) CALL ParallelSumVector( Solver % Matrix, GMVar%Values ,1 )

   ! GL retreat
   Unground=0

  !$omp parallel do &
  !$omp& private(Element,n,NodeIndexes,Material,Density,i,H,rhoi,zb,bedrock,zs) &
  !$omp& reduction(+:Unground) 
   Do t=1,Solver % Mesh % NumberOfBulkElements
      Element => Solver % Mesh % Elements(t)
      n = GetElementNOFNodes(Element)
      NodeIndexes => Element % NodeIndexes
      
      ! Is there a GL node
      IF (.NOT.ANY(abs(GMVar % Values(GMVar % Perm(NodeIndexes(1:n)))).LT.EPS)) CYCLE

      Material => GetMaterial(Element)
      Density(1:n) = &
        ListGetReal( Material, 'SSA Mean Density',n, NodeIndexes,UnFoundFatal=.True.)

      Do i=1,n
       IF (abs(GMVar % Values(GMVar % Perm(NodeIndexes(i)))).GT.EPS) CYCLE
       H=HVar%Values(HVar%Perm(NodeIndexes(i)))
       rhoi=Density(i)
       zb=zsea-H*rhoi/rhow
       bedrock=BedVar%Values(BedVar%Perm(NodeIndexes(i)))
       IF (zb.GT.bedrock) THEN
        GMVar % Values(GMVar % Perm(NodeIndexes(i)))=-1._dp
        Unground=Unground+1

        zs=zb+H
        ZbVar%Values(ZbVar%Perm(NodeIndexes(i)))=zb
        ZsVar%Values(ZsVar%Perm(NodeIndexes(i)))=zs
       END IF
      End Do
      
   End do
   !$omp end parallel do
 
   UngroundTot=Unground
   IF ( ParEnv % PEs>1 ) THEN
     CALL ParallelSumVector( Solver % Matrix, GMVar%Values ,1 )
     CALL MPI_ALLREDUCE( Unground, UngroundTot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
   END IF
   IF (UngroundTot.EQ.0) EXIT

  End do

  IF (SolverTiming) THEN
     at=RealTime()
     st=CPUTime()
     PRINT *,'REAL TIME: ',at-at0
     PRINT *,'CPU  TIME: ',st-st0
     PRINT *,'Number of WHILE LOOPS: ',ndow
  ENDIF


!------------------------------------------------------------------------------
END SUBROUTINE Flotation2
!------------------------------------------------------------------------------
