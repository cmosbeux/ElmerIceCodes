SUBROUTINE boxmodel_solver( Model,Solver,dt,Transient )
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
  TYPE(Variable_t),POINTER :: MeltVar=>NULL(), GMVar=>NULL(), DepthVar=>NULL(), SVar=>NULL(), TVar=>NULL(), BoxVar=>NULL() 
  TYPE(Variable_t),POINTER :: isfslopeVar=>NULL(), distGLVar=>NULL(), distIFVar=>NULL()
  TYPE(Variable_t),POINTER :: TimeVar=>NULL()
  TYPE(Nodes_t) :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  TYPE(Element_t),POINTER ::  Element


  REAL(kind=dp),ALLOCATABLE :: VisitedNode(:),db(:),Basis(:),dBasisdx(:,:) ,Depth(:)
  REAL(kind=dp) :: u,v,w,SqrtElementMetric,s

  INTEGER , POINTER :: MeltPerm(:), GMPerm(:), DepthPerm(:),NodeIndexes(:), SPerm(:), TPerm(:), BPerm(:), loc(:), isfslopePerm(:), distGLPerm(:), distIFPerm(:), Indexx
  INTEGER , DIMENSION(:), ALLOCATABLE :: boxes(:)
  REAL(KIND=dp) , POINTER :: Melt(:),GM(:),isfslope(:), Boxnumber(:), DATAPointer(:,:), distGL(:), distIF(:),DepthVal(:)

  LOGICAL ::  stat, Found,UnFoundFatal=.TRUE.

  CHARACTER(len=MAX_NAME_LEN) ::  variabletype, VariableName
  CHARACTER(len = 200) :: meltValue
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='PICO', FName, DataF, DepthName

  INTEGER :: node,  e, t, n, i, j,  knd, kk, ii,  b,  ierr,  NetcdfStatus,varid,ncid
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER :: maxbastmp, nD,DIM

  REAL(KIND=dp) :: localInteg, Integ, Integ_Reduced, zzz, tmp1, xbox, Tstar,   &
       &                 Area_Reduced, g1, g2,nn,   &
       &                sn, distmax, max_Reduced, TT, SS, qqq_Reduced, totalmelt
  REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: basinmax,basin_Reduced , S0, T0,qqq, localunity,rr
  REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE ::  Mbox, Tbox, Sbox, Zbox, Abox

!!!! SAVE
  TYPE(Variable_t),POINTER, SAVE :: BasinVar=>NULL()
  LOGICAL,SAVE :: Initialized = .FALSE.,  ExtrudedMesh=.False.
  LOGICAL, SAVE :: Firsttime=.TRUE.,llGL, Parallel 
  INTEGER, SAVE :: boxmax, Nmax, MaxBas
  INTEGER :: tmeanid, nlen,nnn
  INTEGER, POINTER, SAVE :: BasinPerm(:)
  REAL(KIND=dp), POINTER, SAVE :: Basin(:)
  REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, SAVE :: S_mean, T_mean
  REAL(KIND=dp), SAVE :: sealevel, lbd1, lbd2, lbd3, meltfac, K, gT,  rhostar, CC,beta, alpha, mskcrit


  !------------------------------------------------------------------------------
  ! 1- Read constants and parameters of the simulation :
  !------------------------------------------------------------------------------

  !- Simulation parameters (idealized ocean forcing) :

  Params => GetSolverParams()

  Mesh => Model % Mesh

  Nmax = Solver % Mesh % NumberOfNodes

  !- get the number of the basin to associate values
  BasinVar => VariableGet( Model % Mesh % Variables, 'basins')
  IF (.NOT.ASSOCIATED(BasinVar)) &
       &    CALL FATAL(SolverName,'basins not found')

  !IF ( BasinVar % TYPE .NE. Variable_on_element) &
  !     &   CALL FATAL(SolverName,'basins not a vairable on element')

  BasinPerm => BasinVar % Perm
  Basin => BasinVar % Values

  IF (Firsttime) THEN
     Firsttime=.False.

     llGL=ListGetLogical( Model % Simulation, 'Grounding Line Melt', UnFoundFatal=UnFoundFatal )

     !- General :

     !- collect constant values from SIF
     sealevel = ListGetCReal( Model % Constants, 'Sea Level',UnFoundFatal=UnFoundFatal)
     meltfac=ListGetCReal( Model % Constants, 'Melt factor',UnFoundFatal=UnFoundFatal)

     lbd1 = ListGetCReal( Model % Constants, 'Liquidus slope',UnFoundFatal=UnFoundFatal )
     lbd2 = ListGetCReal( Model % Constants, 'Liquidus intercept',UnFoundFatal=UnFoundFatal)
     lbd3 = ListGetCReal( Model % Constants, 'Liquidus pressure coeff',UnFoundFatal=UnFoundFatal)

     Parallel = .FALSE. !(ParEnv %PEs > 1)
     maxbastmp=MAXVAL(NINT(Basin))

     IF (Parallel) THEN
        CALL MPI_ALLREDUCE(maxbastmp,MaxBas,1,MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
     ELSE
        MaxBas=maxbastmp
     END IF

     !ALLOCATE(S_mean(MaxBas), T_mean(MaxBas))
     
     !- get values for each basin from NetCDF file
     DataF = ListGetString( Params, 'Data File', Found, UnFoundFatal )

     NetCDFstatus = NF90_OPEN(DataF,NF90_NOWRITE,ncid)
     NetCDFstatus = nf90_inq_dimid(ncid,'number_of_basins' , tmeanid)
     NetCDFstatus = nf90_inquire_dimension(ncid, tmeanid , len=nlen)

     ALLOCATE(S_mean(nlen), T_mean(nlen))

     NetCDFstatus = nf90_inq_varid(ncid,'T_mean',varid)
     NetCDFstatus = nf90_get_var(ncid, varid,T_mean)
     IF ( NetCDFstatus /= NF90_NOERR ) THEN
        CALL Fatal(Trim(SolverName), &
             'Unable to get netcdf variable T_mean')
     END IF

     NetCDFstatus = nf90_inq_varid(ncid,'S_mean',varid)
     NetCDFstatus = nf90_get_var(ncid, varid,S_mean)
     IF ( NetCDFstatus /= NF90_NOERR ) THEN
        CALL Fatal(Trim(SolverName), &
             'Unable to get netcdf variable S_mean')
     END IF

     NetCDFstatus = nf90_inq_varid(ncid,'max box',varid)
     NetCDFstatus = nf90_get_var(ncid, varid,boxmax)
     IF ( NetCDFstatus /= NF90_NOERR ) THEN
        CALL Fatal(Trim(SolverName), &
             'Unable to get netcdf max box')
     END IF

     NetCDFstatus = nf90_inq_varid(ncid,'Circulation_Parameter',varid)
     NetCDFstatus = nf90_get_var(ncid, varid,CC)
     IF ( NetCDFstatus /= NF90_NOERR ) THEN
        CALL Fatal(Trim(SolverName), &
             'Unable to get netcdf Circulation_Parameter')
     END IF

     NetCDFstatus = nf90_inq_varid(ncid,'Effective Exchange Velocity',varid)
     NetCDFstatus = nf90_get_var(ncid, varid,gT)
     IF ( NetCDFstatus /= NF90_NOERR ) THEN
        CALL Fatal(Trim(SolverName), &
             'Unable to get netcdf Effective Exchange Velocity')
     END IF

     NetCDFstatus = nf90_inq_varid(ncid,'Thermal Expansion coeff',varid)
     NetCDFstatus = nf90_get_var(ncid, varid,alpha)
     IF ( NetCDFstatus /= NF90_NOERR ) THEN
        CALL Fatal(Trim(SolverName), &
             'Unable to get netcdf Thermal Expansion coeff')
     END IF

     NetCDFstatus = nf90_inq_varid(ncid,'Salinity contraction coeff',varid)
     NetCDFstatus = nf90_get_var(ncid, varid,beta)
     IF ( NetCDFstatus /= NF90_NOERR ) THEN
        CALL Fatal(Trim(SolverName), &
             'Unable to get netcdf Salinity contraction coeff')
     END IF

     NetCDFstatus = nf90_inq_varid(ncid,'EOS ref Density',varid)
     NetCDFstatus = nf90_get_var(ncid, varid,rhostar)
     IF ( NetCDFstatus /= NF90_NOERR ) THEN
        CALL Fatal(Trim(SolverName), &
             'Unable to get netcdf EOS ref Density')
     END IF
     NetCDFstatus = nf90_close(ncid)

     IF ( llGL ) THEN
        mskcrit =  0.5 ! Melt is at the Grounding Line and floating points
     else
        mskcrit = -0.5 ! No melt at the Grounding Line, only floating points
     ENDif
     CALL INFO(Trim(SolverName),'END FIRST TIME', Level =5)
  END IF
  !-- End of first time check

  CALL INFO(Trim(SolverName),'START', Level =5)

  BoxVar => VariableGet( Model % Mesh % Variables, 'Boxes')
  IF (.NOT.ASSOCIATED(BoxVar)) &
       &    CALL FATAL(SolverName,'boxes not found')

  MeltVar => VariableGet( Model % Mesh % Variables, 'Melt')
  IF (.NOT.ASSOCIATED(MeltVar)) &
       &    CALL FATAL(SolverName,'Melt not found')

  GMVar => VariableGet( Model % Mesh % Variables, 'GroundedMask')
  IF (.NOT.ASSOCIATED(GMVar)) &
       &    CALL FATAL(SolverName,'GroundedMask not found')

  !---- MODIFY FOR STOKES
  DepthName = ListGetString(Params, 'Bottom Surface Variable Name', UnFoundFatal=.TRUE.)
  DepthVar => VariableGet( Model % Mesh % Variables, DepthName,UnFoundFatal=.TRUE.)
  IF (.NOT.ASSOCIATED(DepthVar)) &
       &    CALL FATAL(SolverName,'>Bottom Surface Variable Name< not found')
  !----

  !- get distance to GL
  distGLVar => VariableGet( Model % Mesh % Variables, 'distGL')
  IF (.NOT.ASSOCIATED(distGLVar)) &
       &    CALL FATAL(SolverName,'distGL not found')

  !-- get distance to Ice Front
  distIFVar => VariableGet( Model % Mesh % Variables, 'distIF')
  IF (.NOT.ASSOCIATED(distIFVar)) &
       &    CALL FATAL(SolverName,'distIF not found')

  !Get Permutations and values
  BPerm => BoxVar % Perm
  Boxnumber  => BoxVar % Values

  MeltPerm => MeltVar % Perm
  Melt => MeltVar % Values

  GMPerm => GMVar % Perm
  GM => GMVar % Values

  DepthPerm => DepthVar % Perm
  DepthVal => DepthVar % Values
  Allocate (Depth(SIZE(DepthVal)))
  Depth = sealevel - DepthVal     ! Depth < 0 under sea level

  ALLOCATE( Zbox(boxmax,MaxBas), Abox(boxmax,MaxBas), Tbox(boxmax,MaxBas), Sbox(boxmax,MaxBas), Mbox(boxmax,MaxBas),qqq(MaxBas), T0(MaxBas), S0(MaxBas))
  ALLOCATE(basin_Reduced(MaxBas),basinmax(MaxBas),boxes(MaxBas))
  ALLOCATE(VisitedNode(Nmax),rr(Solver % NumberOfActiveElements),localunity(Solver % NumberOfActiveElements),Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3))

  CALL INFO(Trim(SolverName),'LOAD Variables', Level =5)

  T0 = T_mean 
  S0 = S_mean 
  distGLPerm => distGLVar % Perm
  distGL => distGLVar % Values

  distIFPerm => distIFVar % Perm
  distIF => distIFVar % Values

  VisitedNode=0.0_dp
  Boxnumber(:)=0.0_dp
  totalmelt=0.0_dp
  Abox(:,:)=0.0_dp
  Zbox(:,:)=0.0_dp
  basinmax(:)=0.0_dp
  distmax = 0.0_dp
  basin_Reduced(:)=0.0_dp
  boxes(:)=0.0_dp
  melt(:)=0.0_dp
  localunity(:) = 0.0_dp
  rr(:)=0.0

  ! first loop on element to determine the number of boxes per basins
  CALL INFO(Trim(SolverName),'START BOXES', Level =5)
  nnn = 0
  b=1
  DO e=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(e)
     CALL GetElementNodes( ElementNodes )
     n = GetElementNOFNodes()
     NodeIndexes => Element % NodeIndexes
     Indexx => Element % ElementIndex
     
     ! -check if floating or melting (look at groundedmask)
     IF ( ANY( GM(GMPerm(NodeIndexes(:))) .GE. mskcrit ) ) CYCLE
     
     !b = NINT(Basin(BasinPerm(Indexx)))
     !IF (b.LT.1 .OR. b.GT.18) CYCLE !basin number > 0 & <19
     !check maximal distance to GL at current element (for each basin)
     nnn = nnn+1

     IF (basinmax(b) < MAXVAL(distGL(distGLPerm(NodeIndexes(1:n))))) THEN
        basinmax(b) = MAXVAL(distGL(distGLPerm(NodeIndexes(1:n))))
     END IF

     !check maximal distance to GL for all the basins
     IF (distmax < MAXVAL(distGL(distGLPerm(NodeIndexes(1:n))))) THEN
        distmax = MAXVAL(distGL(distGLPerm(NodeIndexes(1:n))))
     END IF
  END DO
  !write(*,*) 'b, distmax, basinmax, BasMax = ', b , distmax, basinmax(1), MaxBas
  
  ! deal with parallelization
  IF (Parallel) THEN
     CALL MPI_ALLREDUCE(distmax,max_Reduced,1,MPI_DOUBLE_PRECISION,MPI_MAX,ELMER_COMM_WORLD,ierr)
     distmax = max_Reduced
  END IF

  !DO b=1,MaxBas
     IF (Parallel) THEN
        CALL MPI_ALLREDUCE(basinmax(b),basin_Reduced(b),1,MPI_DOUBLE_PRECISION,MPI_MAX,ELMER_COMM_WORLD,ierr)
        basinmax(b) = basin_Reduced(b)
     END IF
  !END DO
  !write(*,*) 'dist max to GL = ', distmax, basinmax(b), b
  !compute the number of boxes per basin (Eq. (9) in Reese et al., 2018)
  boxes = 1+NINT(SQRT(basinmax/distmax)*(boxmax-1))
  !write(*,*) 'basinmax, boxes = ', basinmax, boxes
  CALL INFO(TRIM(SolverName),'Boxes DONE', Level =5)
  
  !call sleep(10)
  !-----------------------------------------------------
  CALL INFO(TRIM(SolverName),'START Area Computation', Level = 5)
  !- Calculate total area of each box (Ak in Reese et al., 2018):
  ! second loop on element
  b=1
  DO e=1,Solver % NumberOfActiveElements

     Element => GetActiveElement(e)
     CALL GetElementNodes( ElementNodes )
     n = GetElementNOFNodes()
     NodeIndexes => Element % NodeIndexes
     Indexx => Element % ElementIndex

     VisitedNode(NodeIndexes(1:n))=VisitedNode(NodeIndexes(1:n))+1.0_dp
     IF ( ANY( GM(GMPerm(NodeIndexes(:))) .GE. mskcrit ) ) CYCLE   ! leave the loop if grounded node in the element

     b = 1 !NINT(Basin(BasinPerm(Indexx)))
     !IF (b.LT.1 .OR. b.GT.18) CYCLE !basin number > 0
     nD = boxes(b)
     !write(*,*) 'num boxes = ', b, nD, boxes
     !non dimensional relative distance to the GL (Eq. (10) in Reese et al., 2018)
     rr(Indexx) = (SUM(distGL(distGLPerm(NodeIndexes(:))))/MAX(1,SIZE(distGL(distGLPerm(NodeIndexes(:))))))   &
          &             / (SUM( distGL(distGLPerm(NodeIndexes(:)))) / MAX(1,SIZE(distGL(distGLPerm(NodeIndexes(:))))) &
          & + SUM( distIF(distIFPerm(NodeIndexes(:)))) / MAX(1,SIZE(distGL(distGLPerm(NodeIndexes(:))))))
     
     IntegStuff = GaussPoints( Element )
     DO t=1,IntegStuff % n
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)
        stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
             Basis,dBasisdx )
        s = SqrtElementMetric * IntegStuff % s(t)
        localunity(Indexx) = localunity(Indexx) + s * SUM(Basis(1:n))  !(surface of the element)
     END DO
    IF (isnan(localunity(Indexx))) THEN
        write(*,*) 'localunity NAN found and changed to 0'
        localunity(Indexx) = 1.0_dp
    ENDIF
     ! check condition for interaction of the box with grid cell coordinate to teh ice cell (Eq. (11))
     DO kk=1,nD
        IF ( rr(Indexx) .GT. 1.0-SQRT(1.0*(nD-kk+1)/nD) .AND. rr(Indexx) .LE. 1.0-SQRT(1.0*(nD-kk)/nD) ) THEN
           Abox(kk,b) = Abox(kk,b) + localunity(Indexx)  !air of box kk in basin b
           Boxnumber(BPerm(Indexx))=kk
        ENDIF
     ENDDO
  END DO !end of loop on elements
  
  !Deal again with parallelization (I deleted a loop on the number of basin to treat 1 basin case)
  ! TO DO: a clean job about the 1 basin case!!!!
  !DO b=1,MaxBas
   b=1
   nD=boxes(b)
     DO kk=1,nD
        IF (Parallel) THEN
           !write(*,*) kk
           CALL MPI_ALLREDUCE(Abox(kk,b),Area_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
           Abox(kk,b) = Area_Reduced
        END IF
     ENDDO
  !END DO
  !write(*,*) 'Abox= ', Abox
  CALL INFO(TRIM(SolverName),'Area Computation DONE', Level =5)

  !write(*,*) 'Air boxes = ', Abox
!-----------------------------------------------------
  ! Compute Tbox, Sbox and qqq and melt for each element of the first box (B1)
  ! Third loop on element
  ! We solve for x = -g1.(Tstar + x - ay) and   (Eqs. A6 and A7)
  CALL INFO(TRIM(SolverName),'START First Box', Level =5)
  Tbox(:,:)=0.d0 ; Sbox(:,:)=0.d0 ; qqq(:)=0.d0
  DO e=1,Solver % NumberOfActiveElements
     b=1
     Element => GetActiveElement(e)
     CALL GetElementNodes( ElementNodes )
     n = GetElementNOFNodes()
     NodeIndexes => Element % NodeIndexes
     Indexx => Element % ElementIndex
     !Only for first box
     IF (  Boxnumber(BPerm(Indexx))==1 ) THEN
        b = 1 !NINT(Basin(BasinPerm(Indexx)))
        zzz = SUM(Depth(DepthPerm(NodeIndexes(:))))/SIZE(NodeIndexes(:)) !mean depth of an element
        Tstar = lbd1*S0(b) + lbd2 + lbd3*zzz - T0(b)  !NB: Tstar should be < 0  (Temperature at the ice-ocean interface; Eq. (5))
        g1 = gT * Abox(1,b)   !exchange velocity
        tmp1 = g1 / (CC*rhostar*(beta*S0(b)*meltfac-alpha))
        sn = (0.5*tmp1)**2 - tmp1*Tstar
        !write(*,*) '(Zb, T*)', zzz, Tstar
        ! to avoid negative discriminent (no solution for x otherwise) :
        IF( sn .lt. 0.d0 ) THEN
           xbox = 0.d0
        else
           xbox = - 0.5*tmp1 + SQRT(sn) ! standard solution (Reese et al)
        ENDif
        TT = T0(b) - xbox
        SS = S0(b) - xbox*S0(b)*meltfac
        Tbox(1,b) = Tbox(1,b) + TT * localunity(Indexx)
        Sbox(1,b) = Sbox(1,b) + SS * localunity(Indexx)
        !write(*,*) 'Tbox, Sbox', Tbox(1,b), Sbox(1,b)
        qqq(b) = qqq(b) + CC*rhostar*(beta*(S0(b)-SS)-alpha*(T0(b)-TT)) * localunity(Indexx) !flux (per basin)
        Melt(MeltPerm(Indexx)) = - gT * meltfac * ( lbd1*SS + lbd2 + lbd3*zzz - TT )
        totalmelt=totalmelt+Melt(MeltPerm(Indexx))* localunity(Indexx)
     END IF

  END DO
  ! deal with parallelisation of loop 3
  !DO b=1,MaxBas
     b = 1
     nD=boxes(b)
     IF (Parallel) THEN
        CALL MPI_ALLREDUCE(Tbox(1,b),Integ_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(Sbox(1,b),Area_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(qqq(b),qqq_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
        Tbox(1,b) = Integ_Reduced
        Sbox(1,b) = Area_Reduced
        qqq(b) = qqq_Reduced
     END IF
  !ENDDO
    
  Tbox(1,1:MaxBas) = Tbox(1,1:MaxBas) / Abox(1,1:MaxBas)
  Sbox(1,1:MaxBas) = Sbox(1,1:MaxBas) / Abox(1,1:MaxBas)
  qqq(1:MaxBas) = qqq(1:MaxBas) / Abox(1,1:MaxBas)
  CALL INFO(TRIM(SolverName),'First Box DONE', Level =5)
  !write(*,*) 'Tbox, Sbox, Melt =', Tbox(1,b), Sbox(1,b), totalmelt

  !-----------------------------------------------------
  ! Temperature, salinity and melt in possible other boxes (B2,...,Bn)
  ! 4th loops on the element
  CALL INFO(TRIM(SolverName),'Start Other Boxes', Level =5)
  DO kk=2,boxmax
     DO e=1,Solver % NumberOfActiveElements
        Element => GetActiveElement(e)
        CALL GetElementNodes( ElementNodes )
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
             Indexx => Element % ElementIndex

        IF (  Boxnumber(BPerm(Indexx))==kk ) THEN
           b= 1!NINT(Basin(BasinPerm(Indexx)))
           zzz = SUM(Depth(DepthPerm(NodeIndexes(:))))/SIZE(NodeIndexes(:)) !mean depth of an element
           Tstar = lbd1*Sbox(kk-1,b) + lbd2 + lbd3*zzz - Tbox(kk-1,b)
           g1  = gT * Abox(kk,b)
           g2  = g1 * meltfac
           xbox = - g1 * Tstar / ( qqq(b) + g1 - g2*lbd1*Sbox(kk-1,b) )
           TT = Tbox(kk-1,b) - xbox
           SS = Sbox(kk-1,b) - xbox*Sbox(kk-1,b)*meltfac
           Tbox(kk,b) =  Tbox(kk,b) + TT * localunity(Indexx)
           Sbox(kk,b) =  Sbox(kk,b) + SS * localunity(Indexx)
           Melt(MeltPerm(Indexx)) = - gT * meltfac * ( lbd1*SS + lbd2 + lbd3*zzz - TT )
           totalmelt=totalmelt+Melt(MeltPerm(Indexx))* localunity(Indexx)
        END IF
     END DO
        b = 1
        !DO b=1,MaxBas
        nD=boxes(b)
        IF (Parallel) THEN
           CALL MPI_ALLREDUCE(Tbox(kk,b),Integ_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
           CALL MPI_ALLREDUCE(Sbox(kk,b),Area_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
           Tbox(kk,b) = Integ_Reduced
           Sbox(kk,b) = Area_Reduced
        END IF
     !ENDDO
     Tbox(kk,1:MaxBas) = Tbox(kk,1:MaxBas) / Abox(kk,1:MaxBas)
     Sbox(kk,1:MaxBas) = Sbox(kk,1:MaxBas) / Abox(kk,1:MaxBas)
  END DO

  IF (Parallel) THEN
     CALL MPI_ALLREDUCE(TotalMelt,Integ_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) THEN
        WRITE(meltValue,'(F20.2)') Integ_Reduced 
        Message='TOTAL_MELT_RATE: '//meltValue
        CALL INFO(SolverName,Message,Level=1)
     END IF
  ELSE
     Integ_Reduced = Integ
  ENDIF
 
  write(*,*) 'Abox, Tbox, Sbox, Melt =', Abox(:,b), Tbox(:,b), Sbox(:,b), totalmelt
  CALL INFO(TRIM(SolverName),'Other Boxes DONE', Level =5)

  DEALLOCATE( Zbox, Abox, Tbox, Sbox, Mbox, T0, S0,rr,localunity)

  CALL INFO(TRIM(SolverName),'Dealocation1', Level =5)
  DEALLOCATE(basin_Reduced,basinmax,boxes)
  CALL INFO(TRIM(SolverName),'Dealocation2', Level =5)
  DEALLOCATE(VisitedNode, Basis, dBasisdx)
  CALL INFO(TRIM(SolverName),'Dealocation3', Level =5)
  CALL INFO(TRIM(SolverName),'PICO MELT CALCULATION DONE', Level =5)
END SUBROUTINE boxmodel_solver


