SUBROUTINE quadratic_solver( Model,Solver,dt,Transient )
  !------------------------------------------------------------------------------
  USE DefUtils
  USE Netcdf

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: Transient
  REAL(KIND=dp) :: dt

  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='QUADRATIC', ParamName

  TYPE(Mesh_t),POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params

  !! variables associated to elements and IPs
  TYPE(Element_t),POINTER ::  Element
  TYPE(Nodes_t),SAVE :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  REAL(kind=dp) :: u, v, w, s, SqrtElementMetric
  REAL(kind=dp),ALLOCATABLE,SAVE :: Basis(:), dBasisdx(:,:)
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER :: Nmax
  INTEGER :: n

  !! Mandatory variables and associated poiners
  TYPE(Variable_t),POINTER :: MeltVar,GMVar,DepthVar,SalVar,TempVar,TFVar,BasinVar
  REAL(KIND=dp), POINTER ::  Melt(:),GM(:),DepthVal(:),Sal(:),Temp(:),TF(:),Basin(:)
  INTEGER , POINTER :: MeltPerm(:),GMPerm(:),DepthPerm(:),SalPerm(:),TempPerm(:),TFPerm(:),BasinPerm(:)
  REAL(kind=dp),ALLOCATABLE :: Depth(:)

  !! Variable relative to netcdf
  INTEGER :: NetcdfStatus,varid,ncid,dimid
  INTEGER, SAVE :: nlen
  CHARACTER(LEN=MAX_NAME_LEN) :: Data_Temp
  REAL(KIND=dp), DIMENSION(:), ALLOCATABLE,SAVE :: delta_T
  
  !! Physical Parameters
  REAL(KIND=dp), SAVE :: sealevel, lambda1, lambda2, lambda3, yearinday, K, beta_S, rhoi, Lf, rhow, cpw, gravity, coriolis, mskcrit
  LOGICAL :: llGL
  LOGICAL,SAVE :: HaveCorrection
  LOGICAL :: DEBUG=.TRUE.

  REAL(KIND=dp) ::  Integ_Reduced, totalmelt, Sloc, Tloc, Surf, ze, slope, Tfr, me

  INTEGER ::  e, t, Indexx, ierr

  LOGICAL ::  stat, Found,UnFoundFatal=.TRUE.
  LOGICAL, SAVE :: Firsttime=.TRUE.
  LOGICAL, SAVE :: Parallel

  CHARACTER(len = 200) :: meltValue

  !------------------------------------------------------------------------------
  ! 
  !------------------------------------------------------------------------------
  Params => GetSolverParams()
  ParamName = ListGetString(Params, 'Slope method', UnFoundFatal=.TRUE.)

  Mesh => Model % Mesh

  Nmax = Solver % Mesh % NumberOfNodes

  
  !------------------------------------------------------------------------------
  ! Get mandatory variables
  !------------------------------------------------------------------------------
  MeltVar => VariableGet(Model % Mesh % Variables, 'Melt',UnFoundFatal=.TRUE.)
  MeltPerm => MeltVar % Perm
  Melt => MeltVar % Values

  ! Check Melt on Elements
  IF (MeltVar % TYPE /= Variable_on_elements) THEN
     CALL FATAL(SolverName, 'Melt is not a variable on elements')
  END IF

  GMVar => VariableGet(Model % Mesh % Variables, 'GroundedMask',UnFoundFatal=.TRUE.)
  GMPerm => GMVar % Perm
  GM => GMVar % Values

  DepthVar => VariableGet(Model % Mesh % Variables, 'Zb', UnFoundFatal=.TRUE.)
  DepthPerm => DepthVar % Perm
  DepthVal => DepthVar % Values

  SalVar => VariableGet(Model % Mesh % Variables, 'salinity', UnFoundFatal=.TRUE.)
  SalPerm => SalVar % Perm
  Sal => SalVar % Values

  TempVar => VariableGet(Model % Mesh % Variables, 'temperature', UnFoundFatal=.TRUE.)
  TempPerm => TempVar % Perm
  Temp => TempVar % Values

  BasinVar => VariableGet(Model % Mesh % Variables, 'basins',UnFoundFatal=.FALSE.)
  BasinPerm => BasinVar % Perm
  Basin => BasinVar % Values

  TFVar => VariableGet(Model % Mesh % Variables, 'TF', UnFoundFatal=.FALSE.)
  TFPerm => TFVar % Perm
  TF => TFVar % Values

  ! Sanity check
  IF (Solver % MeshChanged) &
     CALL FATAL(SolverName,'Mesh has changed not supported ...')

  !------------------------------------------------------------------------------
  ! 1- Initialisation, Read constants and parameters of the simulation :
  !------------------------------------------------------------------------------
  IF (Firsttime) THEN
     CALL INFO(Trim(SolverName),'BEGIN FIRST TIME', Level = 5)
     Firsttime=.False.

     ! - Grounding line :
     llGL          = ListGetLogical( Params, 'Grounding Line Melt', UnFoundFatal = UnFoundFatal )
     IF ( llGL ) THEN
        mskcrit =  0.5 ! Melt is at the Grounding Line and floating points
     ELSE
        mskcrit = -0.5 ! No melt at the Grounding Line, only floating points
     ENDIF

     !- General :
     sealevel      = ListGetCReal( Model % Constants, 'Sea Level', UnFoundFatal = UnFoundFatal )
     yearinday     = ListGetCReal( Model % Constants, 'Calendar', UnFoundFatal = UnFoundFatal )
     rhoi          = ListGetCReal( Model % Constants, 'Ice density', UnFoundFatal = UnFoundFatal )
     Lf            = ListGetCReal( Model % Constants, 'Ice fusion latent heat', UnFoundFatal = UnFoundFatal )
     rhow          = ListGetCReal( Model % Constants, 'Water Density', UnFoundFatal = UnFoundFatal )
     cpw           = ListGetCReal( Model % Constants, 'Sea Water Specific heat', UnFoundFatal = UnFoundFatal )
     gravity       = ListGetCReal( Model % Constants, 'Gravity', UnFoundFatal = UnFoundFatal )
     lambda1       = ListGetCReal( Model % Constants, 'Liquidus slope', UnFoundFatal = UnFoundFatal )
     lambda2       = ListGetCReal( Model % Constants, 'Liquidus intercept', UnFoundFatal = UnFoundFatal )
     lambda3       = ListGetCReal( Model % Constants, 'Liquidus pressure coeff', UnFoundFatal = UnFoundFatal )
     beta_S        = ListGetCReal( Model % Constants, 'Salinity coefficient', UnFoundFatal = UnFoundFatal )
     coriolis      = ListGetCReal( Model % Constants, 'Coriolis', UnFoundFatal = UnFoundFatal )
     
     Parallel = (ParEnv %PEs > 1)
     
     ! - Quadratic
     K  = ListGetCReal( Model % Constants, 'Coeff_K', UnFoundFatal = UnFoundFatal )
     
     !  Correction of temperature
     HaveCorrection = ListGetLogical( Params, 'Correction', Found )
     IF (.NOT.Found) HaveCorrection=.FALSE.
     IF (DEBUG) THEN
         IF (ParEnv % MyPE == 0) THEN
             PRINT *,"Correction found val",Found,HaveCorrection
         END IF
     ENDIF
     IF (HaveCorrection) THEN
         DATA_Temp = ListGetString( Params, 'Data file', Found )
         NetCDFstatus = NF90_OPEN( Trim(Data_Temp), NF90_NOWRITE, ncid )
         ! basin
         NetCDFstatus = nf90_inq_dimid( ncid, 'basins', dimid)
         IF (NetCDFstatus /= NF90_NOERR ) THEN
             CALL Fatal(Trim(SolverName), "dim <basins> not found")
         END IF
         NetCDFstatus = nf90_inquire_dimension( ncid, dimid , len = nlen )
         
         ! Allocate T correction
         ALLOCATE( delta_T(nlen) )

         ! Correction of temperature for each basin
         NetCDFstatus = nf90_inq_varid( ncid, 'delta_T', varid)
         NetCDFstatus = nf90_get_var( ncid, varid, delta_T )
         IF ( NetCDFstatus /= NF90_NOERR ) THEN
             CALL Fatal(Trim(SolverName), 'Unable to get netcdf variable delta_T')
         END IF

         ! close file
         NetCDFstatus=nf90_close(ncid)
     END IF

     ! ALLOCATE arrays with mesh dimensions
     ALLOCATE(Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3))

     CALL INFO(Trim(SolverName),'END FIRST TIME', Level = 5)
  END IF


  !------------------------------------------------------------------------------
  ! 2- Calculation of melt on floating elements :
  !------------------------------------------------------------------------------

  CALL INFO(Trim(SolverName),'MELT CALCULATION', Level = 5)
  
  ! Initialisation
  totalmelt = 0.0_dp
  Melt(:) = 0.0_dp
  TF(:) = 0.0_dp

  ! expression of depth (depth < 0 under sealevel)
  ALLOCATE(Depth(SIZE(DepthVal)))
  Depth = sealevel - DepthVal
  
  ! Expression of thermal forcing at nodes for comparison
  TF = Temp - (lambda1 * Sal + lambda2 + lambda3 * Depth)

  ! Loop on elements
  DO e = 1,Solver % NumberOfActiveElements
     Element => GetActiveElement(e)
     CALL GetElementNodes(ElementNodes, Element, Solver)
     n = GetElementNOFNodes()
     NodeIndexes => Element % NodeIndexes
     Indexx = Element % ElementIndex

     ! we keep only floating element
     IF ( ANY( GM(GMPerm(NodeIndexes(1:n))) .GE. mskcrit ) ) CYCLE
     
     ! conversion nodal value to elemental values
     IntegStuff = GaussPoints( Element )
     ! initialisation
     surf = 0.0_dp
     Tloc = 0.0_dp
     Sloc = 0.0_dp
     slope = 0.0_dp
     Tfr= 0.0_dp
     ze = 0.0_dp
     me = 0.0_dp
     ! Loop on the integration points PI
     DO t = 1,IntegStuff % n
        U = IntegStuff % u(t) ! coordinate x of PI
        V = IntegStuff % v(t) ! coordinate y of PI
        W = IntegStuff % w(t) ! coordinate z of PI
        S = IntegStuff % s(t) ! weigth of PI
        stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric,Basis,dBasisdx)
        ! Surface of the element (sum on integration points)
        surf = surf + S * SqrtElementMetric
        ! T,S & z evaluated on integration points
        ze   = SUM(Depth(DepthPerm(NodeIndexes(1:n))) * Basis(1:n))
        Tloc = SUM(Temp(TempPerm(NodeIndexes(1:n))) * Basis(1:n))
        Sloc = SUM(Sal(SalPerm(NodeIndexes(1:n))) * Basis(1:n))
        ! choice of the slope : local slope (calculation Elmer) ou Antarctic slope (constant value)
        SELECT CASE (TRIM(ParamName))
           CASE ('global')
             slope = ListGetCReal( Model % Constants, 'slope', UnFoundFatal = UnFoundFatal ) ! slope here corresponds to sin(atan(theta))
             slope = tan(asin(slope)) 
           CASE ('local')
             slope  = sqrt( SUM(Depth(DepthPerm(NodeIndexes(1:n))) * dBasisdx(1:n,1))**2 + SUM(Depth(DepthPerm(NodeIndexes(1:n))) * dBasisdx(1:n,2))**2 )
           CASE DEFAULT
             CALL FATAL(SolverName,'slope for parameterisation does not exist.')
        END SELECT
        ! Calculation of freezing point
        Tfr = lambda1 * Sloc + lambda2 + lambda3 * ze
        ! Calculation of melt
        IF (HaveCorrection) THEN
             IF ( nint(Basin(BasinPerm(Indexx))) .GT. nlen ) &
                 CALL FATAL(SolverName, 'number of basins does not match.')

             me = K * (rhow / rhoi) * (cpw / Lf)**2 * beta_S * Sloc * (abs(gravity) / (2 * abs(coriolis))) * sin(atan(slope)) * (abs(Tloc - Tfr + delta_T(Basin(BasinPerm(Indexx))))) * (Tloc - Tfr + delta_T(Basin(BasinPerm(Indexx))))
        ELSE
             me = K * (rhow / rhoi) * (cpw / Lf)**2 * beta_S * Sloc * (abs(gravity) / (2 * abs(coriolis))) * sin(atan(slope)) * (abs(Tloc - Tfr)) * (Tloc - Tfr)
        END IF
        Melt(MeltPerm(Indexx)) = Melt(MeltPerm(Indexx)) + me * S * SqrtElementMetric
      END DO
    ! Calculation of total melt
    totalMelt = totalMelt + Melt(MeltPerm(Indexx))
    ! Melt per surface units
    Melt(MeltPerm(Indexx)) = Melt(MeltPerm(Indexx)) / surf
    END DO
     
    CALL INFO(Trim(SolverName),'END MELT CALCULATION', Level = 5)
    
    ! parallel 
    IF (Parallel) THEN
      CALL MPI_ALLREDUCE(totalMelt,Integ_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
    ELSE
      Integ_Reduced = totalMelt
    ENDIF

    CALL INFO(SolverName,"----------------------------------------", Level=1)
    WRITE(meltValue,'(F20.2)') Integ_Reduced * 917 * 365 / 1.0e12 ! convert m/d in Gt/y
    Message = 'QUADRATIC INTEGRATED BASAL MELT [Gt/y]: '//meltValue
    CALL INFO(SolverName,Message,Level=1)
    CALL INFO(SolverName,"----------------------------------------", Level=1)
  
    ! reverse signe for Elmer (loss of mass (i.e. melt) is negative)
    Melt = -Melt


END SUBROUTINE quadratic_solver

