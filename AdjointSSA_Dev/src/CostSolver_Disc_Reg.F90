!
!  CALL to quadrant tree will be problematiq if we work on a lower dimension
!  than the mesh.....
! * 
! *****************************************************************************
!
! *****************************************************************************
SUBROUTINE CostSolver_Disc_Regularisation( Model,Solver,dt,TransientSimulation )
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
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: DefaultCostFile = 'CostOfT.dat'
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,CostFile,UsedDataFile
  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName,VariableName,ObsFileName, GradSolName
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: TimeVar,CostVar
  TYPE(Variable_t), POINTER :: VelocitybSol,Variable, DJDVariable
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Nodes_t) :: ElementNodes
  REAL(KIND=dp), POINTER :: Vb(:),Values(:), DJDValues(:)
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER, POINTER :: VbPerm(:),Perm(:), DJDPerm(:)
  Logical :: Firsttime=.true.,Found,Parallel,ParallelFile,stat,Gotit, Reset
  integer :: i,j,k,l,s,t,n,NMAX,DIM,ierr,c,ok
  real(kind=dp) :: Cost,Cost_S, Lambda
  real(kind=dp) :: Coord(3),UVW(3),coeff,SqrtElementMetric,x
  REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
  INTEGER,SAVE :: VDOFs
  REAL(KIND=dp),ALLOCATABLE,SAVE :: V(:)
  REAL(KIND=dp),ALLOCATABLE,SAVE :: xobs(:,:),Vobs(:,:)
  INTEGER,ALLOCATABLE,SAVE :: InElement(:)
  INTEGER :: ElementIndex
  INTEGER,SAVE :: NTOT=-1,NTOT_S
  integer,SAVE :: nobs, nobs_found
  LOGICAL,SAVE :: FirstRound=.True.,SAVE_USED_DATA=.False.
  CHARACTER*10 :: date,temps

  INTEGER,PARAMETER :: IO=12

  save Firsttime,Parallel 
  save SolverName,CostSolName,CostFile,VariableName, GradSolName, Lambda
  save ElementNodes


  !PRINT*, firsttime
  WRITE(SolverName, '(A)') 'CostSolver_DataObs'

  SolverParams => GetSolverParams()
  DIM=GetInteger(SolverParams ,'Problem Dimension',Found)
  If (.NOT.Found) then
     CALL WARN(SolverName,'Keyword >Problem Dimension< not found, assume DIM = CoordinateSystemDimension()')
     DIM = CoordinateSystemDimension()
  Endif

  If (Firsttime) then
    N = model % MaxElementNodes
    allocate(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N))
  
!!!!!!! Check for parallel run 
    Parallel = .FALSE.
    IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
            IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
                    Parallel = .TRUE.
            END IF
    END IF

    
!!!!!!!!!!! get Solver Variables
   
  !Print*, SolverName 

  CostFile = ListGetString(Solver % Values,'Cost Filename',Found )
    IF (.NOT. Found) CostFile = DefaultCostFile
    CALL DATE_AND_TIME(date,temps)
    If (Parallel) then
        if (ParEnv % MyPe.EQ.0) then
           OPEN (IO, FILE=CostFile)
                    write(IO,'(a1,a2,a1,a2,a1,a4,5x,a2,a1,a2,a1,a2)') '#',date(5:6),'/',date(7:8),'/',date(1:4), &
                                 temps(1:2),':',temps(3:4),':',temps(5:6)
           CLOSE(IO)
         End if
    Else
           OPEN (IO, FILE=CostFile)
                    write(IO,'(a1,a2,a1,a2,a1,a4,5x,a2,a1,a2,a1,a2)') '#',date(5:6),'/',date(7:8),'/',date(1:4), &
                                 temps(1:2),':',temps(3:4),':',temps(5:6)
           CLOSE(IO)
    End if

   !Name of the Cost Function
   CostSolName =  GetString( SolverParams,'Cost Variable Name', Found)
          IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Cost Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >CostValue<')
                    WRITE(CostSolName,'(A)') 'CostValue'
          END IF
   !Name of the variable
   VariableName =  GetString( SolverParams,'Observed Variable Name', Found)
   IF(.NOT.Found) THEN
       CALL FATAL(SolverName,'Keyword >Observed Variable Name< not found in section >Solver<')
   END IF

   Variable => VariableGet( Solver % Mesh % Variables, Trim(VariableName)  )
   IF ( .NOT.ASSOCIATED( Variable ) ) THEN
       WRITE(Message,'(A,A,A)') & 
                       'No variable >',Trim(VariableName),' < found'
       CALL FATAL(SolverName,Message)
   END IF

   !Get Variable DOFs
   VDOFs=Variable%DOFs

   !Name of the gradient of variable 
   GradSolName =  GetString( SolverParams,'Gradient Variable Name', Found)
   IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Gradient Variable Name< not found in section >Solver<')
              CALL WARN(SolverName,'Taking default value >DJDZb<')
              WRITE(GradSolName,'(A)') 'DJDZb'
   END IF

   !Optional weighting term
   Lambda =  GetConstReal( SolverParams,'Lambda', Found)
   IF(.NOT.Found) THEN
           CALL WARN(SolverName,'Keyword >Lambda< not found  in section >Solver<')
           CALL WARN(SolverName,'Taking default value Lambda=1.0')
           Lambda = 1.0
   End if


!! Do we need to reset cost and DJDVar to 0? Default YES
   Reset =  GetLogical( SolverParams,'Reset Cost Value', Found)
            IF(.NOT.Found) Reset=.True.

 !! Get the obs
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

   allocate(xobs(nobs,3),Vobs(nobs,VDOFs),V(VDOFs),InElement(nobs))
   InElement(:)=-1

   Vobs=0.0_dp
   xobs=0.0_dp
   open(IO,file=trim(ObsFileName))
   do i=1,nobs
     read(IO,*) (xobs(i,j),j=1,DIM),(Vobs(i,j),j=1,VDOFs)
   end do
   close(IO)
   
  !!! End of First visit
    Firsttime=.false.
  Endif

    
!Initialization of variable and derivative (names are saved at firsttime)
    Variable => VariableGet( Solver % Mesh % Variables, TRIM(VariableName)  )
    IF ( ASSOCIATED( Variable ) ) THEN
            Values => Variable % Values
            Perm => Variable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable >',VariableName,' < found'
            CALL FATAL(SolverName,Message)
    END IF
    
    DJDVariable => VariableGet( Solver % Mesh % Variables, TRIM(GradSolName))
    IF ( ASSOCIATED( DJDVariable ) ) THEN
            DJDValues => DJDVariable % Values
            DJDPerm => DJDVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable >',GradSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF
    IF (Reset) DJDValues=0.0_dp

   CALL StartAdvanceOutput(SolverName,'Compute cost')

   Cost=0.0_dp

   !Start Loop on observations to compute Cost and Der. Cost
   nobs_found =0  

   IF (FirstRound) then
     Do s=1,nobs
        CALL AdvanceOutput(s,nobs)
        !Need to find in which Element the data point resides
        ElementIndex=-1  ! don't know why but if i don't reset ElmentIndex it fails
        Coord=0._dp
        Coord(1:DIM)=xobs(s,1:DIM)
        CALL LocateParticleInMeshOctree(ElementIndex,Coord)
        If (ElementIndex.NE.0) InElement(s)=ElementIndex
        IF ((s==nobs).AND.(NTOT < 0)) NTOT=COUNT(InElement(:)>0)     
  
     END Do
   ENDIF !End if FirstRound
   
   Do s=1,nobs 
     CALL AdvanceOutput(s,nobs)
     
     ! Data Point has been found in one element
      IF (InElement(s)>0) THEN
          nobs_found = nobs_found+1
         Element => GetActiveElement(InElement(s))
         n = GetElementNOFNodes()
         NodeIndexes => Element % NodeIndexes
    ! set coords of highest occuring dimension to zero (to get correct path element)
          !-------------------------------------------------------------------------------
         ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
         IF (DIM == 1) THEN !1D SSA
            ElementNodes % y(1:n) = 0.0_dp
            ElementNodes % z(1:n) = 0.0_dp
         ELSE IF (DIM == 2) THEN !2D SSA
            ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
            ElementNodes % z(1:n) = 0.0_dp
         ELSE
            WRITE(Message,'(a,i1,a)')&
                'It is not possible to compute SSA problems with DOFs=',&
                DIM, ' . Aborting'
            CALL Fatal( SolverName, Message)
         END IF
         !Print*, xobs(s,1:2), UVW
         !IF (.NOT.PointInElement(Element,ElementNodes,xobs(s,1:3),UVW))  THEN
         IF (.NOT.PointInElement(Element,ElementNodes,xobs(s,1:DIM),UVW))  THEN
              CALL FATAL(SolverName,'Point was supposed to be found in this element')
         ELSE
            stat = ElementInfo( Element,ElementNodes,UVW(1),UVW(2),UVW(3),SqrtElementMetric, &
                              Basis,dBasisdx )
           ! Variable at obs point
            Do i=1,VDOFs
              V(i)=SUM(Values(VDOFs*(Perm(NodeIndexes(1:n))-1)+i)*basis(1:n))
            End do
            ! Update cost
            Do i=1,VDOFs
              Cost=Cost+0.5*Lambda*(V(i)-Vobs(s,i))*(V(i)-Vobs(s,i))/NTOT
            End do

            !Derivative of Cost at nodal Point
            Do j=1,n
              Do i=1,VDOFs
               k=VDOFS*(DJDPerm(NodeIndexes(j))-1)+i
               DJDValues(k)=DJDValues(k)+Lambda*(V(i)-Vobs(s,i))*basis(j)/NTOT
              End do
            End Do

          END IF

         ELSE

            WRITE(Message,'(a,I0,a)')&
                'Data Point',s,'found in no element'
            CALL Info( SolverName, Message,level=15)
         END IF

    END DO !Do on s
    WRITE(Message, '(A,A,A)')&
                'Number of observations in mesh:', nobs_found , 'points'
    CALL Info( SolverName, Message,level=1) 

     
!    if (NTOT < 0) NTOT=COUNT(InElement(:)>0)
    !Normalization of Cost and Gradient
!    Cost = Cost/NTOT
 !   DJDValues(:) = DJDValues(:)/NTOT
     !print*, 'nobs of zb: ', NTOT

!! Save Cost as a function of time
    TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )

    IF (Parallel) THEN
           CALL MPI_ALLREDUCE(NTOT,NTOT_S,1,&
                  MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)  
   CALL MPI_ALLREDUCE(Cost,Cost_S,1,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  
          CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
          IF (ASSOCIATED(CostVar)) THEN
               IF (Reset) then
                 CostVar % Values(1)=Cost_S
               Else
                 CostVar % Values(1)=CostVar % Values(1)+Cost_S
               End IF 
         END IF
         IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) then
                 OPEN (IO, FILE=CostFile,POSITION='APPEND')
                 write(IO,'(e13.5,2x,e15.8,2x,e15.8)') TimeVar % Values(1),Cost_S,sqrt(2.0*Cost_S/(lambda))
                 CLOSE(IO)
                 write(Message,'(A,A,I0)') trim(SolverName),'total number of data points:',NTOT_S
                 call INFO(SolverName,Message,level=3)
         End if
   ELSE
            CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
            IF (ASSOCIATED(CostVar)) THEN
                 IF (Reset) then
                    CostVar % Values(1)=Cost
                 Else
                    CostVar % Values(1)=CostVar % Values(1)+Cost
                 Endif
            END IF
            OPEN (IO, FILE=CostFile,POSITION='APPEND')
            write(IO,'(e13.5,2x,e15.8,2x,e15.8)') TimeVar % Values(1),Cost,sqrt(2.0*Cost/(lambda))
            close(IO)
            write(Message,'(A,A,I0)') trim(SolverName),'total number of data points:',NTOT
            call INFO(SolverName,Message,level=3)
   END IF
   
   If (FirstRound) THEN
     IF (SAVE_USED_DATA) THEN
      If (Parallel) then
       write(UsedDataFile,'(A,A,I0)') trim(SolverName),'.useddata.',ParEnv % MyPE
      Else
       write(UsedDataFile,'(A,A)') trim(SolverName),'.useddata'
      End if

       open(IO,File=UsedDataFile)
       Do s=1,nobs
         If (InElement(s)>0) then
            if ((DIM.eq.1).AND.(VDOFS.EQ.1)) Then
               write(IO,'(e13.5,2x,e15.8)') (xobs(s,i),i=1,DIM),(Vobs(s,i),i=1,VDOFs)
            Else if ((DIM.eq.2).AND.(VDOFS.EQ.2)) Then
               write(IO,'(e15.8,2x,e15.8,2x,e15.8,2x,e15.8)') (xobs(s,i),i=1,DIM),(Vobs(s,i),i=1,VDOFs)
            End if
         End if
       End do
       close(IO)
     END IF
   End if
   FirstRound=.False.
   Return
!------------------------------------------------------------------------------
END SUBROUTINE CostSolver_Disc_Regularisation
!------------------------------------------------------------------------------
! *****************************************************************************
