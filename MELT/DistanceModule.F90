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
!    Authors: Juha Ruokolainen
!    Email:   Juha.Ruokolainen@csc.fi
!    Web:     http://www.csc.fi/elmer
!    Address: CSC - IT Center for Science Ltd.
!             Keilaranta 14
!             02101 Espoo, Finland
!  
!    Original Date: 09 Nov 2007
!    Modification:  31 May 2023 (reuse the solver as a module for melt param.)
!/*****************************************************************************/


!------------------------------------------------------------------------------
MODULE Distance
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE

    CONTAINS
    
    !--------------------------------------------------------------------
    !> Initialization for the primary solver
    !--------------------------------------------------------------------
    SUBROUTINE DistanceSolver_init( Model,Solver,dt,TransientSimulation, DistName)
      USE DefUtils
      IMPLICIT NONE
    !------------------------------------------------------------------------------
      TYPE(Solver_t) :: Solver
      TYPE(Model_t) :: Model
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation
    !------------------------------------------------------------------------------
      TYPE(ValueList_t), POINTER :: SolverParams
      LOGICAL :: Found
      CHARACTER(LEN=MAX_NAME_LEN) :: DistName
    
      SolverParams => GetSolverParams()
    
      IF( .NOT. ListCheckPresent(SolverParams,'Nonlinear System Convergence Tolerance') ) THEN
        CALL ListAddConstReal(SolverParams,'Nonlinear System Convergence Tolerance',1.0d-8) 
      END IF	
    
    END SUBROUTINE DistanceSolver_init	
    

    !------------------------------------------------------------------------------
    !> Geometric variant of the DistanceSolver.
    !> \ingroup Solvers
    !------------------------------------------------------------------------------
    SUBROUTINE DistanceSolver2( Model,Solver,dt,TransientSimulation, DistName )
    !------------------------------------------------------------------------------
      TYPE(Solver_t) :: Solver
      TYPE(Model_t) :: Model
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation
    !------------------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------------------
      LOGICAL ::  Found, Found1, Found2
      TYPE(Element_t),POINTER :: Element
    
      TYPE(ValueList_t), POINTER :: SolverParams, BC, BodyForce
    
      INTEGER :: n, m, nb, nnb, nd, t, i,j,k,istat, active, MaxIter, maxnode,comm
      REAL(KIND=dp) :: Pnorm,Norm,RelaxDT,TOL,x0,y0,z0,dist
      TYPE(Mesh_t), POINTER :: Mesh

      TYPE(Variable_t),POINTER :: distVar
      REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), buf(:) ,xp(:),yp(:),zp(:),bdist(:),bd(:),&
          condition(:), work(:)
      REAL(KIND=dp), POINTER :: DistVal(:)
      INTEGER, POINTER :: gPerm(:), ibuf(:), aperm(:),bperm(:),cperm(:), DistPerm(:)
      LOGICAL :: DummyDistance 

      CHARACTER(LEN=MAX_NAME_LEN) :: DistName
    
      SAVE STIFF, FORCE
    !------------------------------------------------------------------------------
    
      CALL Info('DistanceSolver2','-------------------------------------------',Level=5)
      CALL Info('DistanceSolver2','Using geometric distance solver',Level=5)
      CALL Info('DistanceSolver2','-------------------------------------------',Level=5)
      CALL ResetTimer('DistanceSolver1')
    
      Mesh => GetMesh()
    
      CALL Info('DistanceSolver2','Working with mesh: '&
          //TRIM(Mesh % Name),Level=5)
      CALL Info('DistanceSolver2','Solving for variable: '&
          //TRIM(DistName),Level=5)
    
      DistVar => VariableGet(Model % Mesh % Variables, DistName, UnFoundFatal=.TRUE.)
      DistPerm => DistVar % Perm
      DistVal => DistVar % Values

      n = Mesh % NumberOfNodes
      ALLOCATE( aperm(n), bperm(n), bdist(n) )
      aperm = 0; bperm = 0; bdist=0
    
      m = Mesh % MaxElementNodes
      ALLOCATE( condition(m), work(m) )
      condition = -1.0_dp; work = 0.0_dp
    
      SolverParams => GetSolverParams()
    
      DummyDistance = GetLogical( SolverParams,'Dummy Distance Computation',Found ) 
    
      nb = 0
      DO t=1,Solver % NumberOfActiveElements
         Element => GetActiveElement(t)
         IF ( .NOT.ASSOCIATED(Element) ) CYCLE
    
         BodyForce => GetBodyForce()
         IF(ASSOCIATED(BodyForce)) THEN
           nd = GetElementNOFNodes(Element)
           work(1:nd) = ListGetReal(BodyForce, TRIM(DistName), nd, Element % NodeIndexes, Found )
           IF ( Found) THEN
              condition(1:nd) = GetReal(BodyForce, TRIM(DistName) // " Condition", Found )        
              DO i=1,nd             
                 j = Element % NodeIndexes(i)
    
                 IF (Found .AND. condition(i) < 0.0) CYCLE

                 IF ( bperm(j) == 0 ) THEN
                    nb = nb + 1
                    aperm(nb) = j
                    bperm(j)  = nb
                    bdist(j) = work(i)
                 END IF
              END DO
           END IF
         END IF
      END DO
    
      !WRITE(*,*) 'Conditions Found, nb of nodes', MAXVAL(bdist(:)), nd
    
      ! Test of Boundary Elements
      DO t=1,Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(t)
    
        IF ( .NOT. ActiveBoundaryElement(Element) ) CYCLE
        BC => GetBC(Element)
        nd = GetElementNOFNodes(Element)
    
        work(1:nd) = GetReal(BC, DistName, Found )
        IF ( Found .OR. GetLogical( BC, 'Noslip Wall BC', Found1 ) ) THEN
          DO i=1,nd
            j = Element % NodeIndexes(i)
            IF ( bperm(j) == 0 ) THEN
              nb = nb + 1
              aperm(nb) = j
              bperm(j)  = nb
              bdist(j) = work(i)
            END IF
          END DO
        END IF
      END DO
    

      IF( ParEnv % PEs == 1 ) THEN
        IF( nb == 0 ) THEN
          CALL Warn('DistanceSolver2','No known distances given for the distance solver!')
          RETURN
        ELSE 
          WRITE( Message,'(A,I0)') 'Number of fixed nodes on bulk: ',nb
          CALL Info('DistanceSolver2',Message,Level=5)
        END IF
      END IF
    
    
      IF ( ParEnv % PEs > 1 ) THEN
        comm = Solver % Matrix % Comm
        gPerm => Mesh % ParallelInfo % GlobalDOFs
    
        maxnode = ParallelReduction(MAXVAL(gPerm(aperm(1:nb))),2)
        ALLOCATE( ibuf(maxnode),cperm(maxnode) );
    
        cperm=0
        cperm(gperm(aperm(1:nb))) = 1
        CALL MPI_ALLREDUCE(cperm, ibuf, maxnode, MPI_INTEGER, MPI_SUM, comm, i)
    
        nnb = 0; cperm=0
        DO i=1,maxnode
         IF ( ibuf(i)/=0 ) THEN
           nnb=nnb+1
           cperm(i)=nnb
           ibuf(nnb)=ibuf(i)
         END IF
        END DO
    
        ALLOCATE( xp(nnb), yp(nnb), zp(nnb), bd(nnb), buf(nnb) )
    
        buf=0._dp
        buf(cperm(gperm(aperm(1:nb)))) = Mesh % Nodes % x(aperm(1:nb))
        CALL MPI_ALLREDUCE(buf, xp, nnb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, i)
        xp(1:nnb) = xp(1:nnb)/ibuf(1:nnb)
    
        buf=0._dp
        buf(cperm(gperm(aperm(1:nb)))) = Mesh % Nodes % y(aperm(1:nb))
        CALL MPI_ALLREDUCE(buf, yp, nnb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, i)
        yp(1:nnb) = yp(1:nnb)/ibuf(1:nnb)
    
        buf=0._dp
        buf(cperm(gperm(aperm(1:nb)))) = Mesh % Nodes % z(aperm(1:nb))
        CALL MPI_ALLREDUCE(buf, zp, nnb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, i)
        zp(1:nnb) = zp(1:nnb)/ibuf(1:nnb)
    
        buf=0._dp
        buf(cperm(gperm(aperm(1:nb)))) = bdist(aperm(1:nb))
        CALL MPI_ALLREDUCE(buf, bd, nnb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, i)
        bd(1:nnb) = bd(1:nnb)/ibuf(1:nnb)
    
        nb = nnb
        DEALLOCATE(ibuf,buf,cperm)
      ELSE
        ALLOCATE( xp(nb),yp(nb),zp(nb),bd(nb) )
    
        DO i=1,nb
          xp(i) = Solver % Mesh % Nodes % x(aperm(i))
          yp(i) = Solver % Mesh % Nodes % y(aperm(i))
          zp(i) = Solver % Mesh % Nodes % z(aperm(i))
          bd(i) = bdist(aperm(i))
        END DO
      END IF
        
      IF( DummyDistance ) THEN
        CALL distcomp0()
      ELSE
        CALL distcomp()
      END IF
    
      IF( ParEnv % PEs == 1 ) THEN
        WRITE( Message,'(A,ES12.3)') 'Maximum distance from given nodes: ',MAXVAL(DistVal)
        CALL Info('DistanceSolver2',Message,Level=8)
      END IF
    
    
      DEALLOCATE(xp,yp,zp,aperm,bperm,condition)
      DistVar % Norm = SQRT(SUM(DistVal**2))
    
      CALL InvalidateVariable( CurrentModel % Meshes, Solver % Mesh, DistName )	
    
      CALL CheckTimer('DistanceSolver2',Delete=.TRUE.)
      CALL Info('DistanceSolver2','All done')
    
    CONTAINS
    
      SUBROUTINE distcomp0
    
        INTEGER :: i,j,k,cnt
        REAL(KIND=dp) :: xl,yl,zl
    
        DistVal = HUGE(1._dp)
        DO i=1,n
          k = DistPerm(i)
          IF ( k <= 0 ) CYCLE
     
          IF ( bperm(i) /= 0 ) THEN
              DistVal(k) = 0._dp
          ELSE
            xl = Mesh % Nodes % x(i)
            yl = Mesh % Nodes % y(i)
            zl = Mesh % Nodes % z(i)
            DO j=1,nb
              dist = (xp(j)-xl)**2 + (yp(j)-yl)**2 + (zp(j)-zl)**2
              DistVal(k) = MIN(DistVal(k), dist)
            END DO
          END IF
        END DO
      END SUBROUTINE distcomp0
    
      SUBROUTINE distcomp
    
        REAL(KIND=dp), ALLOCATABLE :: xxp(:),yyp(:),zzp(:),d(:),dd(:), bbd(:)
        REAL(KIND=dp) :: xl,yl,zl,dl,c(3),ldist
        INTEGER :: i,j,k,l,nnb,cnt,minl
        INTEGER, ALLOCATABLE :: dperm(:), near(:)
    
        nnb = n+nb
        ALLOCATE( xxp(nnb), yyp(nnb), zzp(nnb), d(nnb), dperm(nnb), &
                  near(nnb), dd(nb), bbd(nnb) )
    
        CALL RANDOM_NUMBER(c)
    
        xxp(1:n) = Mesh % Nodes % x
        yyp(1:n) = Mesh % Nodes % y
        zzp(1:n) = Mesh % Nodes % z
        bbd(1:n) = bdist(1:n)
    
        xxp(n+1:nnb) = xp
        yyp(n+1:nnb) = yp
        zzp(n+1:nnb) = zp
        bbd(n+1:nnb) = bd
     
        xl = MAXVAL(xxp)-MINVAL(xxp)
        yl = MAXVAL(yyp)-MINVAL(yyp)
        zl = MAXVAL(zzp)-MINVAL(zzp)
        xxp = xxp + xl*(2*c(1)-1)
        yyp = yyp + yl*(2*c(2)-1)
        zzp = zzp + zl*(2*c(3)-1)
    
        d = SQRT(xxp**2 + yyp**2 + zzp**2)
        dperm = (/ (i,i=1,nnb) /)
    
        CALL SortR(nnb,dperm,d)
    
        j = 0
        DO i=1,nnb
          near(i)=j
          IF ( dperm(i)>n ) THEN
            j=j+1
            xp(j)=xxp(dperm(i))
            yp(j)=yyp(dperm(i))
            zp(j)=zzp(dperm(i))
            dd(j)=d(i)
            bd(j)=bbd(dperm(i))
          END IF
        END DO
    
        DO i=1,nnb
          j=dperm(i)
          IF ( j>n ) CYCLE
    
          k = DistPerm(j)
          IF ( k <= 0 ) CYCLE
    
          IF ( bperm(j)/=0 ) THEN
            DistVal(k) = bdist(j)
          ELSE
            dl = d(i)
            xl = xxp(j)
            yl = yyp(j)
            zl = zzp(j)
            dist = HUGE(1._dp)
            DO l=near(i)+1,nb
              IF( (dd(l)-dl)**2>dist ) EXIT
              ldist = (xp(l)-xl)**2+(yp(l)-yl)**2+(zp(l)-zl)**2
              IF ( dist>ldist) THEN
                minl = l; dist=ldist
              END IF
            END DO
    
            DO l=near(i),1,-1
              IF( (dd(l)-dl)**2>dist ) EXIT
              ldist = (xp(l)-xl)**2+(yp(l)-yl)**2+(zp(l)-zl)**2
              IF ( dist>ldist) THEN
                minl = l; dist=ldist
              END IF
            END DO
    
            DistVal(k) = SQRT(dist)+bd(minl)
    
          END  IF
        END DO
        DEALLOCATE( xxp, yyp, zzp, d, dperm, near, dd, bbd )
      END SUBROUTINE distcomp
    !------------------------------------------------------------------------------
    END SUBROUTINE DistanceSolver2
    !------------------------------------------------------------------------------

    
    !------------------------------------------------------------------------------
    SUBROUTINE IceFrontMask( Model,Solver,dt,TransientSimulation )
    !------------------------------------------------------------------------------
        USE DefUtils
    
        IMPLICIT NONE
        !------------------------------------------------------------------------------
        TYPE(Solver_t) :: Solver
        TYPE(Model_t) :: Model
    
        REAL(KIND=dp) :: dt
        LOGICAL :: TransientSimulation
    
        !------------------------------------------------------------------------------
        ! Local variables
        !------------------------------------------------------------------------------
    
        TYPE(Mesh_t),POINTER :: Mesh
        TYPE(Element_t),POINTER :: Element
        TYPE(ValueList_t), POINTER :: Material, SolverParams, BC
        TYPE(Variable_t), POINTER :: PointerToVariable=>NULL(), HVar=>NULL()
        TYPE(Nodes_t), SAVE :: Nodes
    
        TYPE(ValueList_t), POINTER :: ParentMaterial
        TYPE(Element_t), POINTER :: BoundaryElement, ParentElement
        INTEGER :: other_body_id, nboundary, nparent,BoundaryElementNode, ParentElementNode, istat, k, kk
    
        LOGICAL :: AllocationsDone = .FALSE., GotIt, stat,UnFoundFatal=.TRUE.
        LOGICAL :: FirstTime = .True., NewTime, IsFront
        INTEGER :: DIM,  Nmax,n,t,i
        INTEGER, POINTER :: Permutation(:), HPerm(:), NodeIndexes(:)
        REAL(KIND=dp), POINTER :: VariableValues(:)
    
        LOGICAL :: Parallel
    
        CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'Ice Front Location'
    
        SAVE AllocationsDone, DIM, SolverName
        !------------------------------------------------------------------------------
    
        Mesh => Solver % Mesh
    
        PointerToVariable => VariableGet(Model % Mesh % Variables, 'FrontMask', UnFoundFatal=.TRUE.)
        Permutation  => PointerToVariable % Perm
        VariableValues => PointerToVariable % Values
    
        CALL INFO(SolverName, 'Computing Real Calving front mask', level=3)
    
    
        Parallel = .FALSE.
        IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
        IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
            Parallel = .TRUE.
        END IF
        END IF
    
        !--------------------------------------------------------------
        ! Allocate some permanent storage:
        !--------------------------------------------------------------
    
        k = 0
        kk = 0
        DO t=1,Model % NumberOfBoundaryElements
    
        Element => GetBoundaryElement(t)
        n = GetElementNOFNodes(Element)
        NodeIndexes => Element % NodeIndexes
        !grounded node where Calving Front is false
        BC => GetBC(Element)
        
        DIM = CoordinateSystemDimension()
    
        IsFront = ListGetLogical(BC,'Ice Front',GotIt)
        IF (Isfront) THEN
            VariableValues(Permutation(NodeIndexes(1:n))) = 1.0_dp
        ELSE
            CYCLE
        END IF
    
        END DO
    
        IF ( ParEnv % PEs>1 ) CALL ParallelSumVector( Solver % Matrix, VariableValues, 1 )
    
        CALL INFO( SolverName , 'Done')
        
    END SUBROUTINE IceFrontMask
    !------------------------------------------------------------------------------       

END MODULE Distance
    