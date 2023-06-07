!------------------------------------------------------------------------------
SUBROUTINE NodalNudging( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
  USE ParticleUtils
  USE GeneralUtils
  USE DefUtils
  USE Interpolation
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!  
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,CostFile,UsedDataFile
  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName,VariableName,ObsFileName, NudgedVarName
  CHARACTER(LEN=MAX_NAME_LEN) :: ObsDensityName, TopSurfName, BottomSurfName
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: TimeVar,CostVar
  TYPE(Variable_t), POINTER :: PointerToVariable, Variable, NudgedVar, ObsDensity, &
                               TopSurfVar, BottomSurfVar
  TYPE(ValueList_t), POINTER :: SolverParams, Material 
  TYPE(Nodes_t) :: ElementNodes

  REAL(KIND=dp), POINTER :: VariableValues(:), NVar(:), Density(:) ,&
                            TopSurf(:), BottomSurf(:)
  REAL(Kind=dp) :: LowerLimit
  !REAL(KIND=dp), ALLOCATABLE :: LowerLimit(:)  
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER, POINTER :: NPerm(:),Permutation(:), DensityPerm(:), TopSurfPerm(:),& 
        BottomSurfPerm(:)
  Logical :: Firsttime=.true.,FirstRound=.true., Found,Parallel,ParallelFile,stat,Gotit
  Logical :: FoundLowerLimit
  INTEGER :: i,j,k,l,s,t,n,NMAX,DIM,ierr,c,ok
  REAL(kind=dp) :: Cost,Cost_S
  REAL(kind=dp) :: Coord(3),UVW(3),coeff,SqrtElementMetric,x
  REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model %MaxElementNodes,3)
  INTEGER,SAVE :: VDOFs
  INTEGER :: STDOFs
  REAL(KIND=dp) :: Distance, sigma 

  REAL(KIND=dp),ALLOCATABLE,SAVE :: V(:)
  REAL(KIND=dp),ALLOCATABLE,SAVE :: xobs(:,:),Vobs(:,:)  
  INTEGER,ALLOCATABLE,SAVE :: InElement(:)
  INTEGER :: ElementIndex
  INTEGER,SAVE :: NTOT=-1,NTOT_S
  integer,SAVE :: nobs, M
  LOGICAL,SAVE :: SAVE_USED_DATA=.False.
  CHARACTER*10 :: date,temps

  INTEGER,PARAMETER :: IO=12

  save Firsttime,FirstRound,Parallel
  save SolverName,CostSolName,CostFile,VariableName
  save ElementNodes, LowerLimit

!----------------------------------------------------------------------------- 

  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  STDOFs = PointerToVariable % DOFs
  VDOFs = 1 !a changer si obs vectorielle (v(x,y) ou autre!!!! 
  SolverParams => GetSolverParams()

  ! Get usfull information
  DIM=GetInteger(SolverParams ,'Problem Dimension',Found)
  If (.NOT.Found) then
     CALL WARN(SolverName,'Keyword >Problem Dimension< not found, assume DIM = CoordinateSystemDimension()')
     DIM = CoordinateSystemDimension()
  Endif
  
  ! Variance for gausian callback to obs
  sigma =  GetConstReal(SolverParams ,'Variance Sigma',Found)
  IF (.NOT.Found) THEN
     sigma = 5000.0_dp
     CALL INFO(SolverName, &
                    '>Sigma< Not Found. Setting to 1000.0', Level=5)
  END IF


  If (Firsttime) then
!!!!!!! Check for parallel run 
   Parallel = .FALSE.
   IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
           IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
                   Parallel = .TRUE.
           END IF
   END IF

   WRITE(SolverName, '(A)') 'NodalNudging'

 !! Read file and get the observations
   ObsFileName =  GetString( SolverParams,'Observation File Name', Found)
   IF(.NOT.Found) THEN
       CALL FATAL(SolverName,'Keyword >Observation File Name< not found in section >Solver<')
   END IF
   ParallelFile = .False.
   ParallelFile = GetLogical(SolverParams,'Parallel Observation Files', Found)
   if (Parallel.AND.ParallelFile) &
    write(ObsFileName,'(A,A,I0)') trim(ObsFileName),'.',ParEnv % MyPE


   SAVE_USED_DATA = GetLogical(SolverParams,'Save used data', Found)

   open(IO,file=trim(ObsFileName),status = 'old',iostat = ok)
   if(ok /= 0) then
       write(message,'(A,A)') 'Unable to open file ',TRIM(ObsFileName)
       CALL Fatal(Trim(SolverName),Trim(message))
   end if
   nobs=0
   do while(ok == 0)
     read(io,*,iostat = ok)
     if (ok == 0) nobs = nobs + 1
   end do
   close(IO)

   N = model % MaxElementNodes
   M = Solver % Mesh % NumberOfNodes

   !Allocation of node postion vector, obs postion, obs value, mod value at obs 
   allocate(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N))   
   allocate(xobs(nobs,3),Vobs(nobs,VDOFs),V(VDOFs),InElement(nobs))
   !Allocatation of array for ice thickness limit
   !allocate(LowerLimit(n))

   InElement(:)=-1
   Vobs=0.0_dp
   xobs=0.0_dp
   open(IO,file=trim(ObsFileName))
   do i=1,nobs
     read(IO,*) (xobs(i,j),j=1,DIM),(Vobs(i,j),j=1,VDOFs) 
   end do
   close(IO)

 Firsttime=.False.
 END IF !End of firsttime
   
  ! Take information about the Nudged Variable : need value for each step
  NudgedVarName = GetString( SolverParams,'Nudged Variable Name', Found)
  If(.NOT.Found) THEN
       CALL FATAL(SolverName,'Keyword >Nudged Variable Name< not found in &
                section >Solver<')
  END IF

  NudgedVar => VariableGet( Solver % Mesh % Variables, NudgedVarName  )
  IF ( ASSOCIATED( NudgedVar ) ) THEN
        NVar => NudgedVar % Values
        NPerm => NudgedVar % Perm
  ELSE
        WRITE(Message,'(A,A,A)') 'No nudged variable >',NudgedVarName,' < found'
        CALL FATAL(SolverName,Message)
  END IF

  ! Take information about the top surface 
  TopSurfName = GetString( SolverParams,'Top Surface Name', Found)
  If(.NOT.Found) THEN
       CALL INFO (SolverName,'Keyword >Top Surface  Name< not found in &
                section >Solver<')
  ELSE  
        TopSurfVar => VariableGet( Solver % Mesh % Variables, TopSurfName  )
        IF ( ASSOCIATED( TopSurfVar ) ) THEN
                TopSurf => TopSurfVar % Values
                TopSurfPerm => TopSurfVar % Perm
        ELSE
                WRITE(Message,'(A,A,A)') 'No top surface variable >',TopSurfName,' < found'
                CALL FATAL (SolverName,Message)
        END IF 
  END IF
  
  ! Take information about the bottom surface
  BottomSurfName = GetString( SolverParams,'Bottom Surface Name', Found)
  If(.NOT.Found) THEN
       CALL FATAL(SolverName,'Keyword >Bottom Surface Name< not found in &
                section >Solver<')
  ELSE
        BottomSurfVar => VariableGet( Solver % Mesh % Variables, BottomSurfName  )
        IF ( ASSOCIATED( BottomSurfVar ) ) THEN
                BottomSurf => BottomSurfVar % Values
                BottomSurfPerm => BottomSurfVar % Perm
        ELSE
                WRITE(Message,'(A,A,A)') 'No variable >',BottomSurfName,' < found'
                CALL FATAL (SolverName,Message)
        END IF
  END IF

  !Variable giving the observation Density
  ObsDensityName = GetString( SolverParams,'Observation Density Variable Name', Found)
  If(.NOT.Found) THEN
      CALL INFO (SolverName,'Keyword >Observation Density Variable Name< & 
                not found in section >Solver<')
  ELSE
       ObsDensity => VariableGet( Solver % Mesh % Variables, ObsDensityName  )
       IF ( ASSOCIATED( ObsDensity ) ) THEN
                Density => ObsDensity % Values
                DensityPerm => ObsDensity % Perm
       ELSE
                WRITE(Message,'(A,A,A)') 'No variable >',ObsDensityName,' <found'
                CALL FATAL(SolverName,Message)
       END IF
  END IF
  
  CALL StartAdvanceOutput(SolverName,'Nodal Nudging')

  VariableValues(:) = 0.0_dp
  !Loop on observations
  DO s=1,nobs
        CALL AdvanceOutput(s,nobs)
        
        !Need to find in which Element the data point resides
        IF (FirstRound) THEN
                ElementIndex=-1  ! don't know why but if i don't reset ElmentIndex it fails
                Coord=0._dp
                Coord(1:DIM)=xobs(s,1:DIM)
                CALL LocateParticleInMeshOctree(ElementIndex, Coord) !Find particle in mesh       
                If (ElementIndex.NE.0) InElement(s)=ElementIndex  ! Fill array
        ENDIF 
        
        !now we know in which element is each observation by reading the array 
        ! InElement(:) of the size of the number of obs. 

        !Data Point has been found in one element
        IF (InElement(s)>0) THEN ! If < 0 the particle is not in the mesh
                Element => GetActiveElement(InElement(s))
                n = GetElementNOFNodes()
                NodeIndexes => Element % NodeIndexes

                ! set coords of highest occuring dim. to zero (to get correct path element)
                !--------------------------------------------------------------------------
                ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes(1:n))
                IF (DIM == 1) THEN !1D 
                  ElementNodes %  y(1:n) = 0.0_dp
                  ElementNodes % z(1:n) = 0.0_dp
                ELSE IF (DIM == 2) THEN !2D
                  ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes(1:n))
                  ElementNodes % z(1:n) = 0.0_dp
                ELSE
                  WRITE(Message,'(a,i1,a)')&
                  'It is not possible to apply nudging with DOFs=',&
                              DIM, ' . Aborting'
                  CALL Fatal( SolverName, Message)
                END IF
                ! Check if obs is really in Element
                IF (.NOT.PointInElement(Element,ElementNodes,xobs(s,1:3),UVW))  THEN
                  CALL FATAL(SolverName,'Point was supposed to be found in this element')
  
                ELSE
                  stat = ElementInfo(Element,ElementNodes,UVW(1),UVW(2),UVW(3),SqrtElementMetric, &
                             Basis,dBasisdx )
                  DO k=1,n !Loop on node of local element
                    DO j=1,VDOFs
                                V(j)=SUM(NVar(VDOFs*(NPerm(NodeIndexes(1:n))-1)+j)*basis(1:n))
                                !compute distance to each obs
                                IF (DIM==1) THEN
                                        Distance = abs(xobs(s,1) - ElementNodes % x(k))
                                ELSE 
                                        Distance = ((xobs(s,1) - ElementNodes % x(k))**2  &
                                                   +(xobs(s,2) - ElementNodes % y(k))**2)**0.5 
                                END IF 
                                coeff = exp(-0.5*(Distance/sigma)**2)   
                                VariableValues(Permutation(Nodeindexes(k))) = VariableValues(Permutation(Nodeindexes(k)))+ coeff*(V(j)-Vobs(s,j))                  
                                !Obs density only computed first round
                                IF (FirstRound) THEN
                                  Density(DensityPerm(nodeindexes(k))) = Density(DensityPerm(nodeindexes(k))) +1  
                                END IF
                    END DO
                  END DO
                END IF
        END IF
  END DO
  FirstRound=.False.

  !Get Info about minimum ice thickness
  !Material => GetMaterial()
  !LowerLimit=0.0_dp
  !LowerLimit(1:n) = GetReal(Material,'Minimum Ice Thickness', FoundLowerLimit)
   !        IF (.Not.FoundLowerLimit)  &
   !           CALL FATAL( SolverName, '> Minimum Ice Thickness < not given in Material' )

  Material => GetMaterial() 
  LowerLimit = GetConstReal(Material, 'Minimum Ice Thickness', FoundLowerLimit)
           IF (.Not.FoundLowerLimit)  &
              CALL FATAL( SolverName, '> Minimum Ice Thickness < not given in Material' )
 
  !Allocation of Callback value for each node
  DO t=1,Solver % Mesh % NumberOfNodes
        IF (Density(DensityPerm(t)) == 0) THEN
           CYCLE
        ELSE
           !Avoid ice thickness < lowerlimit
           VariableValues(Permutation(t)) = VariableValues(Permutation(t))/Density(DensityPerm(t)) 
           IF ( (BottomSurfName .eq. NudgedVarName) .AND. ((TopSurf(TopSurfPerm(t)) - VariableValues(Permutation(t))) < LowerLimit) ) THEN
                VariableValues(Permutation(t)) = TopSurf(TopSurfPerm(t)) - LowerLimit
                WRITE(Message,'(A)') &
                'Bottom Surf. going greater than Top Surf. >> corrected' 
                CALL INFO(SolverName, Message)
           ELSE IF ( (TopSurfName .eq. NudgedVarName) .AND. (( VariableValues(Permutation(t)) - BottomSurf(BottomSurfPerm(t))) < LowerLimit) ) THEN
                VariableValues(Permutation(t)) = BottomSurf(BottomSurfPerm(t)) + LowerLimit 
                WRITE(Message,'(A)') &
                'Top Surf. going lower than Bottom Surf. >> corrected'            
                CALL INFO(SolverName, Message)
           ELSE IF (((TopSurfName .ne. NudgedVarName) .AND. (BottomSurfName .ne. NudgedVarName) &
                .AND. VariableValues(Permutation(t)) < LowerLimit)) THEN
                !Assume NudgedVariable is the ice thickness
                VariableValues(Permutation(t)) = LowerLimit 
                WRITE(Message,'(A)') &
                'Ice thickness going under Lower limit >> corrected'                      
                CALL INFO(SolverName, Message) 
           END IF 
                
        END IF
  END DO

RETURN
!------------------------------------------------------------------------------
END SUBROUTINE NodalNudging 
!------------------------------------------------------------------------------
!******************************************************************************
