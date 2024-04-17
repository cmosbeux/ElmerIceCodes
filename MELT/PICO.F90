!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! *
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! *
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
! *  Boston, MA 02110-1301, USA.
! *
!/*****************************************************************************/
!  
!    Authors: Urruty, Caillet, Gillet-Chaulet, Mosbeux
!    Email:   cmosbeux@univ-gremoble-alpes.fr
!    Web:     http://www.csc.fi/elmer
!    Address: CSC - IT Center for Science Ltd.
!             Keilaranta 14
!             02101 Espoo, Finland
!  
!    Original Date: 2021
!    Modification:  31 May 2023 (Works both in 2D (SSA) and 3D (Stokes))
!/*****************************************************************************/

!------------------------------------------------------------------------------
MODULE PICO
!------------------------------------------------------------------------------
  !USE CoordinateSystems
  !USE MeshUtils
  USE Netcdf
  USE DefUtils

  CONTAINS 

  !------------------------------------------------------------------------------
  SUBROUTINE boxmodel_solver( Model,Solver,dt,Transient )
  !------------------------------------------------------------------------------
  
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(Model_t)  :: Model
    TYPE(Solver_t), TARGET :: Solver
    LOGICAL ::  Transient
    REAL(KIND=dp) :: dt

    !------------------------------------------------------------------------------
    !    Local variables
    !------------------------------------------------------------------------------
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='PICO', DepthName

    TYPE(Mesh_t),POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Params

    !! variables associated to elements and IPs
    TYPE(Element_t),POINTER ::  Element
    TYPE(Nodes_t),SAVE :: ElementNodes
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    REAL(kind=dp) :: u, v, w, SqrtElementMetric, s, surf
    REAL(kind=dp),ALLOCATABLE,SAVE :: Basis(:), dBasisdx(:,:)
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: Nmax
    INTEGER :: n

    !! Mandatory variables and associated pointers
    TYPE(Variable_t),POINTER :: BasinVar, BoxVar, MeltVar, GMVar, DepthVar, distGLVar, distIFVar
    REAL(KIND=dp), POINTER ::  Basin(:), Boxnumber(:), Melt(:), GM(:), DepthVal(:), distGL(:), distIF(:)
    INTEGER , POINTER :: BasinPerm(:), BPerm(:), MeltPerm(:), GMPerm(:), DepthPerm(:), distGLPerm(:), distIFPerm(:)
    REAL(KIND=dp),ALLOCATABLE :: Depth(:)

    !! Variables related to netcd (containing input variables)
    INTEGER :: NetcdfStatus,varid,ncid
    CHARACTER(LEN=MAX_NAME_LEN) :: DataFT, DataFS
    INTEGER :: tmeanid, nlen
    INTEGER,SAVE :: nTime
    INTEGER :: nTime2

    INTEGER,SAVE :: VisitedTimes=0, TimeOffset
    INTEGER :: TimePoint
    REAL(KIND=dp) :: Time

    !! Physical Parameters
    REAL(KIND=dp), SAVE :: sealevel, lbd1, lbd2, lbd3, yearinday, meltfac, K, gT,  rhostar, CC,beta, alpha, mskcrit, scaling
    INTEGER, SAVE :: boxmax,MaxBas
    LOGICAL :: llGL, PanAntarctic, MeltNodal

    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: S_mean, T_mean
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: Zbox, Abox, Tbox, Sbox
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, SAVE :: qqq, T0, S0, basin_Reduced, basinmax, delta_T
    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: boxes
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, SAVE :: localunity,rr

    REAL(KIND=dp) ::  Integ_Reduced, zzz, tmp1, xbox, Tstar,   &
        &                Area_Reduced, g1, g2, &
        &                sn, max_Reduced, TT, SS, qqq_Reduced, totalmelt
    REAL(KIND=dp) :: distmax,dmax

    INTEGER ::  e, t, kk, b
    INTEGER :: nD
    INTEGER :: maxbastmp
    INTEGER :: ierr
    INTEGER ::  Indexx

    LOGICAL ::  stat, Found,UnFoundFatal=.TRUE.
    LOGICAL, SAVE :: Firsttime=.TRUE.
    LOGICAL, SAVE :: Parallel

    CHARACTER(len = 200) :: meltValue

    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    Params => GetSolverParams()

    Mesh => Model % Mesh  

    Nmax = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
    
    !------------------------------------------------------------------------------
    ! Get mandatory variables
    !------------------------------------------------------------------------------

    BasinVar => VariableGet(Model % Mesh % Variables, 'Basins', UnFoundFatal=.TRUE.)
    BasinPerm => BasinVar % Perm
    Basin => BasinVar % Values

    BoxVar => VariableGet(Model % Mesh % Variables, 'Boxes', UnFoundFatal=.TRUE.)
    BPerm => BoxVar % Perm
    Boxnumber => BoxVar % Values

    MeltVar => VariableGet(Model % Mesh % Variables, 'Melt', UnFoundFatal=.TRUE.)
    MeltPerm => MeltVar % Perm
    Melt => MeltVar % Values

    ! cm: for now, the 3d version works only with nodal variable, this should be fixed/cleaned in the future
    MeltNodal = ListGetLogical( Params,'Nodal Melt',Found) 
    IF (MeltNodal) THEN
      IF (BasinVar % TYPE == Variable_on_elements .OR. BoxVar % TYPE == Variable_on_elements &
        & .OR. MeltVar % TYPE == Variable_on_elements) THEN
        CALL FATAL(SolverName, 'Basins, Boxes, or Melt is not a variable on nodes')
      END IF
    ELSE
      IF (BasinVar % TYPE /= Variable_on_elements .OR. BoxVar % TYPE /= Variable_on_elements &
        & .OR. MeltVar % TYPE /= Variable_on_elements) THEN
        CALL FATAL(SolverName, 'Basins, Boxes, or Melt is not a variable on elements')
      END IF
    END IF

    GMVar => VariableGet( Model % Mesh % Variables, 'GroundedMask',UnFoundFatal=.TRUE.)
    GMPerm => GMVar % Perm
    GM => GMVar % Values

    ! cm: for modulablility, we ask the user to give the name of the bottom surface variable 
    DepthName = ListGetString(Params, 'Bottom Surface Variable Name', UnFoundFatal=.TRUE.)
    DepthVar => VariableGet( Model % Mesh % Variables, DepthName,UnFoundFatal=.TRUE.)
    IF (.NOT.ASSOCIATED(DepthVar)) &
        &    CALL FATAL(SolverName,'>Bottom Surface Variable Name< not found')
    DepthPerm => DepthVar % Perm
    DepthVal => DepthVar % Values
    
    distGLVar => VariableGet( Model % Mesh % Variables, 'distGL',UnFoundFatal=.TRUE.)
    distGLPerm => distGLVar % Perm
    distGL => distGLVar % Values

    distIFVar => VariableGet( Model % Mesh % Variables, 'distIF',UnFoundFatal=.TRUE.)
    distIFPerm => distIFVar % Perm
    distIF => distIFVar % Values

    ! Sanity check
    IF (Solver % MeshChanged) &
      CALL FATAL(SolverName,'Mesh has changed not supported ...')

    !------------------------------------------------------------------------------
    ! 1 - Initialisation, Read constants and parameters of the simulation :
    !------------------------------------------------------------------------------
    IF (Firsttime) THEN
      Firsttime=.False.

      ! - Grounding line :
      llGL      = ListGetLogical( Params, 'Grounding Line Melt', UnFoundFatal = UnFoundFatal )
      ! - Offset for reading data :
      TimeOffset= ListGetInteger( Params, 'Time Counter start', UnFoundFatal = UnFoundFatal )

      ! - PanAntarctic or Regional (important for the number of basins)
      PanAntarctic = ListGetLogical( Params, 'PanAntarctic', UnFoundFatal=UnFoundFatal )

      !- General :
      sealevel = ListGetCReal( Model % Constants, 'Sea Level', UnFoundFatal = UnFoundFatal )
      lbd1     = ListGetCReal( Model % Constants, 'Liquidus slope', UnFoundFatal = UnFoundFatal )
      lbd2     = ListGetCReal( Model % Constants, 'Liquidus intercept', UnFoundFatal = UnFoundFatal )
      lbd3     = ListGetCReal( Model % Constants, 'Liquidus pressure coeff', UnFoundFatal = UnFoundFatal )
      yearinday = ListGetCReal( Model % Constants, 'Calendar', UnFoundFatal = UnFoundFatal )

      ! - PICO : 
      boxmax   = ListGetInteger( Model % Constants, 'Nb Boxes', UnFoundFatal = UnFoundFatal )
      CC       = ListGetCReal( Model % Constants, 'Overturning Coefficient', UnFoundFatal = UnFoundFatal )
      gT       = ListGetCReal( Model % Constants, 'Temperature Exchange Velocity', UnFoundFatal = UnFoundFatal )
      alpha    = ListGetCReal( Model % Constants, 'Thermal Expansion Coefficient EOS', UnFoundFatal = UnFoundFatal )
      beta     = ListGetCReal( Model % Constants, 'Salinity Contraction Coefficient EOS', UnFoundFatal = UnFoundFatal )
      rhostar  = ListGetCReal( Model % Constants, 'In Situ Density EOS', UnFoundFatal = UnFoundFatal )
      meltfac  = ListGetCReal( Model % Constants, 'Melt Factor', UnFoundFatal = UnFoundFatal )

      scaling  = ListGetCReal( Model % Constants, 'Scaling', Found )
      IF(.NOT.Found) THEN
          CALL WARN(SolverName,'Keyword >Scaling< not found  in section >Constant<')
          CALL WARN(SolverName,'Taking default value Scaling = 1')
          Scaling = 1.0
      END IF

      
      !cy: orignal version looks at the maximal basin number and loop over the basins. this works for PanAntarctic
      ! simulations but not for region simulations where we do not simulate all the basins
      maxbastmp = MAXVAL(NINT(Basin))

      Parallel = (ParEnv %PEs > 1)
      IF (Parallel) THEN
          CALL MPI_ALLREDUCE(maxbastmp,MaxBas,1,MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
      ELSE
          MaxBas=maxbastmp
      END IF

      !! Allocate arrays related to PICO
      ALLOCATE(Zbox(boxmax,MaxBas), Abox(boxmax,MaxBas), Tbox(boxmax,MaxBas), Sbox(boxmax,MaxBas))
      ALLOCATE(qqq(MaxBas), T0(MaxBas), S0(MaxBas))
      ALLOCATE(basin_Reduced(MaxBas),basinmax(MaxBas),boxes(MaxBas))

      !! ALLOCATE arrays with mesh dimensions
      !ALLOCATE(rr(Solver % NumberOfActiveElements), localunity(Solver % NumberOfActiveElements), Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3))
      !cy: changed allocation for nodal variables
      ALLOCATE(rr(Nmax),localunity(Nmax),Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3))
      ALLOCATE(Depth(SIZE(DepthVal)))


      !!------------------------------------------------------------------------------
      ! GET FORCINGS
      !!------------------------------------------------------------------------------     

      DataFT = ListGetString( Params, 'data file', Found, UnFoundFatal )

      NetCDFstatus = NF90_OPEN( Trim(DataFT), NF90_NOWRITE, ncid )
      NetCDFstatus = nf90_inq_dimid( ncid, 'number_of_basins' , tmeanid)
      IF (NetCDFstatus /= NF90_NOERR ) THEN
          CALL Fatal(Trim(SolverName), &
            "dim <number_of_basins> not found")
      ENDIF

      NetCDFstatus = nf90_inquire_dimension( ncid, tmeanid , len = nlen )
    
      IF (nlen.NE.MaxBas .AND. PanAntarctic) THEN
        CALL Fatal(Trim(SolverName),"Number of basins does not agree")
      ELSE IF (nlen.EQ.MaxBas .AND. PanAntarctic) THEN
        CALL INFO(Trim(SolverName), 'PanAntarctic Simulation', Level = 5)
      ELSE
        CALL INFO(Trim(SolverName), 'Regional Simulation', Level = 5)
      ENDIF

      !!! check if we have a time dimension
      NetCDFstatus = nf90_inq_dimid( ncid, 'time' , varid)
      IF (NetCDFstatus.EQ.NF90_NOERR ) THEN
          NetCDFstatus = nf90_inquire_dimension( ncid, varid , len = nTime )
      ELSE
          nTime = 1
      ENDIF
  
      ALLOCATE( T_mean(nlen,nTime) )
      ALLOCATE( S_mean(nlen,nTime) )
      ALLOCATE( delta_T(nlen) ) 

      !initialize delta T to zero to make sure it does not come at a weird value when not in the netcdf
      delta_T(:) = 0.0_dp  
      NetCDFstatus = nf90_inq_varid( ncid, 'delta_T', varid)
      IF ( NetCDFstatus.EQ.NF90_NOERR ) THEN
          NetCDFstatus = nf90_get_var( ncid, varid, delta_T )
      ELSE
          CALL INFO(Trim(SolverName), 'Unable to get netcdf variable delta_T. Set to 0.0', Level = 3)
          delta_T(:) = 0.0_dp   
      END IF

      NetCDFstatus = nf90_inq_varid( ncid, 'T_mean', varid)
      NetCDFstatus = nf90_get_var( ncid, varid, T_mean )
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
          CALL Fatal(Trim(SolverName), &
              'Unable to get netcdf variable T_mean')
      END IF

      NetCDFstatus = nf90_inq_varid( ncid, 'S_mean', varid)
      NetCDFstatus = nf90_get_var( ncid, varid, S_mean )
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
          CALL Fatal(Trim(SolverName), &
              'Unable to get netcdf variable S_mean')
      END IF

      ! close file
      NetCDFstatus=nf90_close(ncid)
      
      IF ( llGL ) THEN
          mskcrit =  0.5 ! Melt is at the Grounding Line and floating points
      else
          mskcrit = -0.5 ! No melt at the Grounding Line, only floating points
      ENDif
      CALL INFO(Trim(SolverName),'END FIRST TIME', Level =5)
    END IF

    CALL INFO(Trim(SolverName),'START', Level =5)

    Depth = sealevel - DepthVal     ! Depth < 0 under sea level

    IF (nTime.GT.1) THEN
    ! get time index
        VisitedTimes = VisitedTimes + 1
        IF( ListGetLogical( Params, "Is Time Counter", Found ) ) THEN
          TimePoint = VisitedTimes + TimeOffset
        ELSE
          TimePoint = ListGetInteger( Params, "Time Index", Found )
          IF (.NOT.Found) THEN
            Time = GetTime()
            dt = GetTimeStepSize()
            !everything is computed for daily units but we can scale it depending on the target unit
            TimePoint = floor((time/(yearinday/scaling))-(dt/(yearinday/scaling))/2) + 1 + TimeOffset
          END IF
        END IF
        !TimePoint = max(1,min(TimePoint,nTime))
        CALL INFO(Trim(SolverName),"Use Time Index: "//I2S(TimePoint), Level=3)
    ELSE
      TimePoint=1
    ENDIF

    IF (TimePoint.GT.nTime) THEN
         CALL Fatal(Trim(SolverName), &
              'Unable to read the line of T and S in netCDF')
    END IF

    T0(1:MaxBas) = T_mean(1:MaxBas,TimePoint) + delta_T(1:MaxBas)
    IF (VisitedTimes == 1) THEN
      write(*,*) T0(1)
     END IF  

     S0(1:MaxBas) = S_mean(1:MaxBas,TimePoint)

    Boxnumber(:) = 0.0_dp
    totalmelt = 0.0_dp
    Abox(:,:) = 0.0_dp
    Zbox(:,:) = 0.0_dp
    basinmax(:) = 0.0_dp
    distmax = 0.0_dp
    basin_Reduced(:) = 0.0_dp
    boxes(:) = 0.0_dp
    Melt(:) = 0.0_dp
    localunity(:) = 0.0_dp
    rr(:) = 0.0_dp

    !!------------------------------------------------------------------------------
    ! 2 - DEFINE BOXES FOR EACH BASIN
    !!------------------------------------------------------------------------------
    
    CALL INFO(Trim(SolverName),'START BOXES', Level =5)

    ! first loop on element to determine the number of boxes per basins
    DO e=1,Solver % NumberOfActiveElements

      Element => GetActiveElement(e)
      CALL GetElementNodes( ElementNodes, Element, Solver)
      n = GetElementNOFNodes()
      NodeIndexes => Element % NodeIndexes
      Indexx = Element % ElementIndex
      
      ! check if floating or melting (look at groundedmask)
      IF ( ANY( GM(GMPerm(NodeIndexes(:))) .GE. mskcrit ) ) CYCLE

      IF (MeltNodal) THEN
        b = NINT(MAXVAL(Basin(BasinPerm(NodeIndexes(1:n)))))
      ELSE
        b = NINT(Basin(BasinPerm(Indexx)))
      ENDIF
  
      ! check maximal distance to GL at current element (for each basin)
      dmax = MAXVAL(distGL(distGLPerm(NodeIndexes(1:n))))
      IF (basinmax(b) < dmax ) THEN
          basinmax(b) = dmax
      END IF

      !check maximal distance to GL for all the basins
      dmax = MAXVAL(distGL(distGLPerm(NodeIndexes(1:n))))
      IF (distmax < dmax) THEN
          distmax = dmax
      END IF
    END DO
    
    IF (Parallel) THEN
      CALL MPI_ALLREDUCE(distmax,max_Reduced,1,MPI_DOUBLE_PRECISION,MPI_MAX,ELMER_COMM_WORLD,ierr)
      distmax= max_Reduced

      CALL MPI_ALLREDUCE(basinmax,basin_Reduced,MaxBas,MPI_DOUBLE_PRECISION,MPI_MAX,ELMER_COMM_WORLD,ierr)
      basinmax(1:MaxBas) = basin_Reduced(1:MaxBas)
    END IF
    
    !compute the number of boxes per basin (Eq. (9) in Reese et al., 2018)
    boxes = 1+NINT(SQRT(basinmax/distmax)*(boxmax-1))

    CALL INFO(TRIM(SolverName),'BOXES DONE', Level = 5)

    !!------------------------------------------------------------------------------
    ! 3 - DEFINE AREA OF EACH BOX
    !!------------------------------------------------------------------------------
    !- Calculate total area of each box (Ak in Reese et al., 2018):
    ! second loop on elements
    CALL INFO(TRIM(SolverName),'START Area Computation', Level = 5)

    DO e=1,Solver % NumberOfActiveElements

      Element => GetActiveElement(e)
      CALL GetElementNodes( ElementNodes, Element, Solver )
      n = GetElementNOFNodes()
      NodeIndexes => Element % NodeIndexes
      Indexx = Element % ElementIndex

      IF ( ANY( GM(GMPerm(NodeIndexes(:))) .GE. mskcrit ) ) CYCLE    ! leave the loop if grounded point in the element

      IF (MeltNodal) THEN
        b = NINT(MAXVAL(Basin(BasinPerm(NodeIndexes(1:n)))))
        rr(NodeIndexes(:)) = (SUM(distGL(distGLPerm(NodeIndexes(:)))) /MAX(1,SIZE(distGL(distGLPerm(NodeIndexes(:))))))   &
        &      / (SUM(distGL(distGLPerm(NodeIndexes(:)))) / MAX(1,SIZE(distGL(distGLPerm(NodeIndexes(:))))) &
        & + SUM( distIF(distIFPerm(NodeIndexes(:)))) / MAX(1,SIZE(distGL(distGLPerm(NodeIndexes(:))))))
      ELSE
        b = NINT(Basin(BasinPerm(Indexx)))
        !non dimensional relative distance to the GL (Eq. (10) in Reese et al., 2018)
        rr(Indexx) = (SUM(distGL(distGLPerm(NodeIndexes(:))))/MAX(1,SIZE(distGL(distGLPerm(NodeIndexes(:))))))   &
        &             / (SUM( distGL(distGLPerm(NodeIndexes(:)))) / MAX(1,SIZE(distGL(distGLPerm(NodeIndexes(:))))) &
        & + SUM( distIF(distIFPerm(NodeIndexes(:)))) / MAX(1,SIZE(distGL(distGLPerm(NodeIndexes(:))))))
      ENDIF
      
      nD=boxes(b)
      
      IntegStuff = GaussPoints( Element )
      DO t=1,IntegStuff % n
          U = IntegStuff % u(t)
          V = IntegStuff % v(t)
          W = IntegStuff % w(t)
          stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
              Basis,dBasisdx )
          s = SqrtElementMetric * IntegStuff % s(t)
        ! Surface of the element
          IF (MeltNodal) THEN
            localunity(NodeIndexes(1:n)) = localunity(NodeIndexes(1:n)) +  s * Basis(1:n)
          ELSE
            localunity(Indexx) = localunity(Indexx) + s * SUM(Basis(1:n))
          ENDIF
      END DO
    
      ! check cond for box interation with grid cell coordinate (Eq. (11))
      IF (MeltNodal) THEN
        DO kk=1,nD
          IF (MAXVAL(rr(NodeIndexes(1:n))) .GT. (1.0 - SQRT(1.0*(nD-kk+1)/nD)) &
            & .AND. MAXVAL(rr(NodeIndexes(1:n))) .LE. (1.0 - SQRT(1.0*(nD-kk)/nD))) THEN
            Abox(kk,b) = Abox(kk,b) + SUM(localunity(NodeIndexes(1:n))) / SIZE(NodeIndexes(:))   ! air of box kk in basin b
            Boxnumber(BPerm(NodeIndexes(1:n))) = kk
          ENDIF
        
        ENDDO
      ELSE
        DO kk=1,nD
          IF ( rr(Indexx) .GT. 1.0-SQRT(1.0*(nD-kk+1)/nD) .AND. rr(Indexx) .LE. 1.0-SQRT(1.0*(nD-kk)/nD) ) THEN
            Abox(kk,b) = Abox(kk,b) + localunity(Indexx)  !air of box kk in basin b
            Boxnumber(BPerm(Indexx)) = kk
          ENDIF
        ENDDO
      ENDIF
    END DO  !end of loop on elements
    
    ! cm: Here we do a loop over the number of basins. Works for (1,...,Nmax) but not for, e.g., (3,7,12)
    IF (Parallel) THEN
      DO b=1,MaxBas
        nD=boxes(b)
        DO kk=1,nD
          CALL MPI_ALLREDUCE(Abox(kk,b),Area_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
          Abox(kk,b) = Area_Reduced
          WRITE(message, '(A,I0,A,I0,A,F7.1,A)') 'Area(', kk, ',', b, ') = ', Abox(kk, b)/1e6 , ' km2'
          CALL INFO(TRIM(SolverName), message, Level = 5)
        ENDDO
      END DO
    END IF
    
    CALL INFO(TRIM(SolverName),'Area Computation DONE', Level = 5)
    
    !!------------------------------------------------------------------------------
    ! 4 - COMPUTE THE MELT IN BOXE 1 (based on input parameters and depth)
    !!------------------------------------------------------------------------------
    ! Compute Tbox, Sbox and qqq and melt for each element of the first box (B1)
    ! We solve for x = -g1.(Tstar + x - ay) and   (Eqs. A6 and A7)
    ! Third loop on elements
    CALL INFO(TRIM(SolverName),'STARTING Melt Boxe 1', Level = 5)
    Tbox(:,:)=0.d0 ; Sbox(:,:)=0.d0 ; qqq(:)=0.d0
    DO e=1,Solver % NumberOfActiveElements
      Element => GetActiveElement(e)
      CALL GetElementNodes( ElementNodes )
      n = GetElementNOFNodes()
      NodeIndexes => Element % NodeIndexes
      Indexx = Element % ElementIndex

      !Check that we are we work in 2D or 3D
      IF (MeltNodal .AND. MAXVAL(Boxnumber(BPerm(NodeIndexes(1:n))))==1 ) THEN
        b = NINT(MAXVAL(Basin(BasinPerm(NodeIndexes(1:n)))))
        surf = SUM(localunity(NodeIndexes(1:n)))/n 
      ELSEIF (.NOT. MeltNodal .AND. Boxnumber(BPerm(Indexx))==1 ) THEN
        b = NINT(Basin(BasinPerm(Indexx)))
        surf = localunity(Indexx)
      ELSE
        CYCLE
      ENDIF
      
      zzz = SUM(Depth(DepthPerm(NodeIndexes(1:n))))/n !mean depth of an element
      Tstar = lbd1*S0(b) + lbd2 + lbd3*zzz - T0(b)  !NB: Tstar should be < 0 (Temperature at the ice-ocean interface; Eq. (5))
      g1 = gT * Abox(1,b) !exchange velocity
      tmp1 = g1 / (CC*rhostar*(beta*S0(b)*meltfac-alpha))
      sn = (0.5*tmp1)**2 - tmp1*Tstar
      ! to avoid negative discriminent (no solution for x otherwise) :
      IF( sn .LT. 0.d0 ) THEN
        xbox = 0.d0
      ELSE
        xbox = - 0.5*tmp1 + SQRT(sn) ! standard solution (Reese et al)
      END IF
      TT = T0(b) - xbox 
      SS = S0(b) - xbox*S0(b)*meltfac 
      Tbox(1,b) = Tbox(1,b) + TT *  surf
      Sbox(1,b) = Sbox(1,b) + SS *  surf
      qqq(b) = qqq(b) + CC*rhostar*(beta*(S0(b)-SS)-alpha*(T0(b)-TT)) * surf !flux (per basin)

      IF (MeltNodal) THEN
        Melt(MeltPerm(NodeIndexes(1:n))) = - gT * meltfac * ( lbd1*SS + lbd2 + lbd3*zzz - TT ) 
        totalmelt = totalmelt &
          & + SUM(Melt(MeltPerm(NodeIndexes(1:n))))/n &
          & * SUM(localunity(NodeIndexes(1:n)))/n
      ELSE 
        Melt(MeltPerm(Indexx)) = - gT * meltfac * ( lbd1*SS + lbd2 + lbd3*zzz - TT )
        totalmelt = totalmelt + Melt(MeltPerm(Indexx)) * localunity(Indexx)
      ENDIF 
    END DO
    
    DO b=1,MaxBas
      nD=boxes(b)
      IF (Parallel) THEN
        CALL MPI_ALLREDUCE(Tbox(1,b),Integ_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(Sbox(1,b),Area_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(qqq(b),qqq_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
        Tbox(1,b) = Integ_Reduced
        Sbox(1,b) = Area_Reduced
        qqq(b) = qqq_Reduced
      END IF
    ENDDO

    Tbox(1,1:MaxBas) = Tbox(1,1:MaxBas) / Abox(1,1:MaxBas)
    Sbox(1,1:MaxBas) = Sbox(1,1:MaxBas) / Abox(1,1:MaxBas)
    qqq(1:MaxBas) = qqq(1:MaxBas) / Abox(1,1:MaxBas)

    CALL INFO(TRIM(SolverName),'Melt Boxe 1 DONE', Level = 5)

    !!------------------------------------------------------------------------------
    ! 5 - COMPUTE THE MELT IN OTHER BOXES (based on input parameters and depth)
    !!------------------------------------------------------------------------------
    ! Temperature, salinity and melt in possible other boxes (B2,...,Bn)
    ! 4 loops on the elements
    CALL INFO(TRIM(SolverName),'START Melt Other Boxes', Level = 5)
    DO kk=2,boxmax
      DO e=1,Solver % NumberOfActiveElements
        Element => GetActiveElement(e)
        CALL GetElementNodes( ElementNodes )
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
        Indexx = Element % ElementIndex

        !Check that that melt is nodal or elemental
        IF (MeltNodal .AND. MAXVAL(Boxnumber(BPerm(NodeIndexes(1:n))))==kk ) THEN
          b = NINT(MAXVAL(Basin(BasinPerm(NodeIndexes(1:n)))))
          surf = SUM(localunity(NodeIndexes(1:n)))/n 
        ELSEIF (.NOT. MeltNodal .AND. Boxnumber(BPerm(Indexx))==kk ) THEN
          b = NINT(Basin(BasinPerm(Indexx)))
          surf = localunity(Indexx)
        ELSE
          CYCLE
        ENDIF

        zzz = SUM(Depth(DepthPerm(NodeIndexes(:))))/n !mean depth of an element
        Tstar = lbd1*Sbox(kk-1,b) + lbd2 + lbd3*zzz - Tbox(kk-1,b)
        g1  = gT * Abox(kk,b)
        g2  = g1 * meltfac
        xbox = - g1 * Tstar / ( qqq(b) + g1 - g2*lbd1*Sbox(kk-1,b) )
        TT = Tbox(kk-1,b) - xbox
        SS = Sbox(kk-1,b) - xbox*Sbox(kk-1,b)*meltfac
        Tbox(kk,b) =  Tbox(kk,b) + TT * surf
        Sbox(kk,b) =  Sbox(kk,b) + SS * surf

        IF (MeltNodal) THEN
          Melt(MeltPerm(NodeIndexes(1:n))) = - gT * meltfac * ( lbd1*SS + lbd2 + lbd3*zzz - TT )  
          totalmelt = totalmelt &
            & + SUM(Melt(MeltPerm(NodeIndexes(1:n))))/n &
            & * SUM(localunity(NodeIndexes(1:n)))/n
        ELSE 
          Melt(MeltPerm(Indexx)) = - gT * meltfac * ( lbd1*SS + lbd2 + lbd3*zzz - TT )
          totalmelt = totalmelt + Melt(MeltPerm(Indexx)) * localunity(Indexx)
        ENDIF 
       
      END DO

      DO b=1,MaxBas
          nD=boxes(b)
          IF (Parallel) THEN
            CALL MPI_ALLREDUCE(Tbox(kk,b),Integ_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
            CALL MPI_ALLREDUCE(Sbox(kk,b),Area_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
            Tbox(kk,b) = Integ_Reduced
            Sbox(kk,b) = Area_Reduced
          END IF
      ENDDO

      Tbox(kk,1:MaxBas) = Tbox(kk,1:MaxBas) / Abox(kk,1:MaxBas)
      Sbox(kk,1:MaxBas) = Sbox(kk,1:MaxBas) / Abox(kk,1:MaxBas)

    END DO

    IF (Parallel) THEN
      CALL MPI_ALLREDUCE(TotalMelt,Integ_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
    ELSE
      Integ_Reduced = TotalMelt
    ENDIF
    
    CALL INFO(TRIM(SolverName),'Melt Other Boxes DONE', Level = 5)

    CALL INFO(SolverName,"----------------------------------------", Level=1)
    WRITE(meltValue,'(F20.2)') Integ_Reduced*0.917/1.0e9 * Scaling
    Message = 'PICO INTEGRATED BASAL MELT [Gt/a]: '//meltValue   ! 0.917/1.0e9 to convert m3/a in Gt/a
    CALL INFO(SolverName, Message, Level=1)
    WRITE(meltValue,'(F20.2)') Integ_Reduced*0.917/1.0e9
    Message = 'PICO INTEGRATED BASAL MELT [Gt/j] (rho=917): '//meltValue 
    CALL INFO(SolverName, Message, Level=1)
    CALL INFO(SolverName,"----------------------------------------", Level=1)
    
    ! reverse signe for Elmer (loss of mass (i.e. melt) is negative)
    Melt = -Melt

    ! scaling if indicated
    Melt = Scaling * Melt

  END SUBROUTINE boxmodel_solver

END MODULE PICO

