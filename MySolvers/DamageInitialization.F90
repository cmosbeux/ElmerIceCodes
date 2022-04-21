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
! *****************************************************************************/
!
!/******************************************************************************
! ******************************************************************************
! *
! *  Authors: Cyrille Mosbeux
! *  Email:   cyrille.mosbeux@univ-grenoble-alpes.fr
! *  Web:     http://www.csc.fi/elmer
! *
! *  Original Date: 28 Mars 2022
! *
! *****************************************************************************/


SUBROUTINE DamageInitializationSolver( Model,Solver,dt,TransientSimulation )

USE DefUtils
IMPLICIT NONE !----------------------------------------------------------------

TYPE(Solver_t) :: Solver
TYPE(Model_t) :: Model
REAL(KIND=dp) :: dt
LOGICAL :: TransientSimulation

!----------------------------------------------------------------
! Local variables
!----------------------------------------------------------------

TYPE(Element_t),POINTER :: Element
LOGICAL :: AllocationsDone = .FALSE., Found, Converged, &
    LimitSolution, FoundLowerLimit, FoundUpperLimit, UnFoundFatal=.TRUE.

LOGICAL, ALLOCATABLE ::  LimitedSolution(:,:), ActiveNode(:,:)

INTEGER :: i, n, t, k, l, istat, other_body_id, iter, NonlinearIter, Active, NMAX, MMAX
INTEGER :: dummyInt, Indexes(128)
INTEGER :: CorrectedLowerLimit,CorrectedUpperLimit
INTEGER , POINTER :: NodeIndexes(:), DPerm(:)
REAL(KIND=dp), POINTER :: Damage(:), ForceVector(:), PointerToResidualVector(:)

REAL(KIND=dp) :: Norm, PrevNorm, LinearTol, NonlinearTol, RelativeChange, OriginalValue
TYPE(ValueList_t), POINTER :: BodyForce, Material, BC, SolverParams
TYPE(Variable_t), POINTER :: Dresidual

REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), LOAD(:), FORCE(:), &
    UpperLimit(:), LowerLimit(:), ResidualVector(:), OldValues(:), OldRHS(:), &
    StiffVector(:)

TYPE(Matrix_t), POINTER :: Systemmatrix

CHARACTER(LEN=MAX_NAME_LEN) :: BoundaryType
CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName = 'DamageInitialization'
CHARACTER(LEN=MAX_NAME_LEN) :: VariableName

SAVE MASS, STIFF, LOAD, FORCE,&
    AllocationsDone, PrevNorm, UpperLimit, LowerLimit, VariableName, LimitedSolution, ActiveNode, ResidualVector, &
    OldValues, OldRHS, StiffVector


! Get variable name
!-----------------------
VariableName = TRIM(Solver % Variable % Name)
Damage => Solver % Variable % Values     ! Nodal values
IF (.NOT.ASSOCIATED(Damage)) CALL Fatal(SolverName, 'Variable values not associated')
DPerm => Solver % Variable % Perm       ! Permutations index

SystemMatrix => Solver % Matrix
ForceVector => Solver % Matrix % RHS

! Allocate some permanent storage, this is done first time
!--------------------------------------------------------------
IF ( .NOT. AllocationsDone ) THEN
    NMAX = Solver % Mesh % MaxElementNodes !big enough for elemental arrays
    MMAX = Model % Mesh % NumberOfNodes
    K = SIZE( SystemMatrix % Values )
    L = SIZE( SystemMatrix % RHS )

    WRITE(*,*) 'NMAX, MMAX:', NMAX, MMAX
    ALLOCATE( FORCE(NMAX), LOAD(NMAX), MASS(NMAX,NMAX), STIFF(NMAX,NMAX), &
        LimitedSolution(MMAX,2), &
        Lowerlimit(MMAX), &
        UpperLimit(MMAX), &
        ActiveNode(MMAX,2), &
        ResidualVector(L), &
        OldValues(K),&
        OldRHS(L),&
        StiffVector(L),&
        STAT=istat )

    IF ( istat /= 0 ) THEN
        CALL Fatal( 'Damage initialization','Memory allocation error for matrix/vectors.' )
    END IF

    CALL Info(SolverName,'Memory allocations done' )
    AllocationsDone = .TRUE.
    ActiveNode = .FALSE.
    ResidualVector = 0.0_dp
    Damage = 0.0_dp
END IF
LimitedSolution=.FALSE.


!    Get variables for the residual
!------------------------------------
DResidual => VariableGet( Model % Mesh % Variables, TRIM(VariableName) // ' Residual',UnFoundFatal=UnFoundFatal)
PointerToResidualVector => DResidual % Values


!Read in solver parameters
!-------------------------
SolverParams => GetSolverParams()
IF (.NOT. ASSOCIATED(SolverParams)) &
    CALL FATAL('Damage initialization','No Solver section found')

LinearTol = GetConstReal( SolverParams, &
       'Linear System Convergence Tolerance',    Found )
  IF ( .NOT.Found ) THEN
     CALL Fatal(SolverName, 'No >Linear System Convergence Tolerance< found')
  END IF


NonlinearIter = GetInteger(SolverParams, & 'Nonlinear System Max Iterations', Found)
IF ( .NOT.Found ) NonlinearIter = 1
NonlinearTol = GetConstReal( SolverParams, 'Nonlinear System Convergence Tolerance', Found )
IF ( .NOT.Found ) NonlinearTol = 1.0D-03

!Get limiter on the solution (hard limit are Damage > 0 and Damage<1)
LimitSolution = GetLogical( SolverParams, 'Limit Solution', Found )
IF ( .NOT. Found ) LimitSolution = .FALSE.
IF (LimitSolution) THEN
      CALL Info(SolverName, 'Keyword > Limit Solution < found. Solution will be limited',Level=6)
ELSE
      CALL Info(SolverName, 'No keyword > Limit Solution < found. Solution will not be limited',Level=6)
END IF


!----------------------------------------------------------------
! Nonlinear iteration loop
!----------------------------------------------------------------
DO iter=1,NonlinearIter
    Converged = .FALSE.
    WRITE(Message,'(A,I5,A,I5)') 'Nonlinear iteration no.',iter,&
          ' of max. ', NonlinearIter
     CALL INFO('DamageInitialization', Message,level=1)
    !Initialize the system and do the assembly:
    !------------------------------------------
    CALL DefaultInitialize()

!----------------------------------------------------------------
! Assembly for the domain
!----------------------------------------------------------------
    CorrectedUpperLimit = 0
    CorrectedLowerLimit = 0
    DO t=1,Solver % NumberOfActiveElements

        !Get element info
        !-----------------
        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
    
        !Get material parameters
        !------------------------
        Material => GetMaterial()
        IF (.NOT. ASSOCIATED(Material)) THEN
        WRITE(Message,'(A,I5,A)') 'No material for bulk element no. ',t,' found.'
            CALL FATAL('DamageInitialization',Message)
        END IF
    
        !Limited Solution
        !-------------------
        IF (LimitSolution) THEN
            UpperLimit(Element % Nodeindexes(1:n)) = ListGetReal(Material,TRIM(VariableName) // ' Upper Limit', n, Element % NodeIndexes, FoundUpperLimit)
            LimitedSolution(Element % Nodeindexes(1:n), 1) = FoundUpperLimit

            LowerLimit(Element % Nodeindexes(1:n)) = ListGetReal(Material,TRIM(VariableName) // ' Lower Limit', n, Element % NodeIndexes, FoundLowerLimit)
            LimitedSolution(Element % Nodeindexes(1:n), 2) = FoundLowerLimit
        END IF

        !Get load for force vector (source)
        !----------------------------------
        LOAD = 0.0d0
        BodyForce => GetBodyForce()
        IF ( ASSOCIATED(BodyForce) ) LOAD(1:n) = GetReal( BodyForce, 'Damage Source', Found )
        
!        DO k=1,n
!            IF (LOAD(k)<0.0) WRITE(*,*) LOAD(k)
!        END DO


        !Get element local matrix and rhs vector:
        !----------------------------------------
        CALL LocalMatrix( MASS, STIFF, FORCE, LOAD, Element, n,TransientSimulation)

        !Update global matrix and rhs vector from local matrix and vector
        !---------------------------------------------------------------
        IF ( TransientSimulation ) CALL Default1stOrderTime( MASS,STIFF,FORCE )
        !------------------------------------------------------------------------------
        !      Update global matrix and rhs vector from local matrix & vector
        !------------------------------------------------------------------------------
        CALL DefaultUpdateEquations( STIFF, FORCE )

    !----------------------------------------------------------------
    END DO ! end Assembly for the domain
    !----------------------------------------------------------------

    CALL DefaultFinishBulkAssembly()

!----------------------------------------------------------------
! Assembly for the Neumann boundary conditions (Not used right now)
!----------------------------------------------------------------
!    DO t=1, Solver % Mesh % NumberOfBoundaryElements
!        Element => GetBoundaryElement(t)
!        IF ( .NOT. ActiveBoundaryElement() ) CYCLE
!        n = GetElementNOFNodes()
!        ! no evaluation of Neumann BC’s on points
!        IF ( GetElementFamily() == 1 ) CYCLE
!        BC => GetBC()
!
!        FORCE = 0.0d00
!        MASS = 0.0d00
!        STIFF = 0.0d00
!
!        ! check type of boundary and set BC accordingly
!        !----------------------------------------------
!        BoundaryType = GetString(BC,'Boundary Type',Found)
!        IF (.NOT. Found) CYCLE
!    END DO ! end Assembly for Neumann boundary conditions
!----------------------------------------------------------------

    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    !------------------------------------------------------------------------------
    ! limit solution
    !------------------------------------------------------------------------------

    IF (LimitSolution) THEN
        OldValues = SystemMatrix % Values
        OldRHS = ForceVector

        ! manipulation of the matrix
        !---------------------------
        DO i=1,Model % Mesh % NumberOfNodes
            k = DPerm(i)
            IF ((k > 0) .AND. (ActiveNode(i,1) .OR. ActiveNode(i,2))) THEN
                CALL ZeroRow( SystemMatrix, k )
                CALL SetMatrixElement( SystemMatrix, k, k, 1.0_dp )
                IF (ActiveNode(i,1)) THEN
                    SystemMatrix % RHS(k) = LowerLimit(i) +1e-5
                ELSE IF (ActiveNode(i,2)) THEN
                    SystemMatrix % RHS(k) = UpperLimit(i)
                END IF
            END IF
        END DO
    END IF

    ! Solve the system
    ! -----------------------------------------------------------------------------

    PrevNorm = Solver % Variable % Norm
    Norm = DefaultSolve()

    ! compute relative change of norm
    !------------------------------------------------------------------------------
    IF ( PrevNorm + Norm /= 0.0d0 ) THEN
        RelativeChange = 2.0d0 * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
    ELSE
        RelativeChange = 0.0d0
    END IF
    
    WRITE( Message, * ) 'Result Norm : ',Norm
    CALL Info( 'DamageInitialization', Message, Level=4 )
    WRITE( Message, * ) 'Relative Change : ', RelativeChange
    CALL Info( 'DamageInitialization', Message, Level=4 )


    ! compute residual
    !------------------------------------------------------------------------------
    !-----------------------------
    ! determine "active" nodes set
    !-----------------------------
    IF (LimitSolution) THEN

        SystemMatrix % Values = OldValues
        ForceVector = OldRHS

        !to modify for parallel
        CALL CRS_MatrixVectorMultiply( SystemMatrix, Damage, StiffVector)
        ResidualVector =  StiffVector - ForceVector

        DO i=1,Model % NumberOfNodes
            l= DPerm(i)
             IF (l<1) CYCLE
             !---------------------------------------------------------
             ! if upper limit is exceeded, manipulate matrix in any case
             !----------------------------------------------------------
             IF ((LimitedSolution(i,1)).AND.(Damage(l)-LowerLimit(i)<0.0_dp)) THEN
                ActiveNode(i,1) = .TRUE.
                WRITE(*,*) Damage(l)
             END IF
             IF ((LimitedSolution(i,2)).AND.(Damage(l)-UpperLimit(i)>0.0_dp)) THEN
                ActiveNode(i,2) = .TRUE.
             END IF

             IF ( LimitedSolution(i,1) .AND. ResidualVector(l) < -LinearTol &
                      .AND. iter>1 ) ActiveNode(i,1) = .FALSE.
             IF ( LimitedSolution(i,2) .AND. ResidualVector(l) >  LinearTol &
                      .AND. iter>1 ) ActiveNode(i,2) = .FALSE.

             IF( .NOT.ActiveNode(i,1) .AND. .NOT.ActiveNode(i,2) ) THEN
                PointerToResidualVector(DResidual % Perm(i)) = 0.0_dp
             ELSE
                PointerToResidualVector(DResidual % Perm(i)) = ResidualVector(l)
            ENDIF
        END DO
    
        ! write info about corrected nodes
        WRITE(Message,'(a,e13.6,a,e13.6)') &
       'Max/min Damage:', MAXVAL(Damage(:)),'/', MINVAL( Damage(:))
        CALL Info( 'DamageInitialization', Message, Level=4 )

    END IF

    ! check for convergence
    ! ----------------------
    IF ( RelativeChange < NonlinearTol ) THEN ! NO
        Converged = .TRUE.
        EXIT
    ELSE ! YES
        PrevNorm = Norm
    END IF

!----------------------------------------------------------------
END DO ! of the nonlinear iteration loop
!----------------------------------------------------------------

! write info about corrected nodes and convergence
!---------------------------------------------------
WRITE(Message,'(a,e13.6,a,e13.6)') &
'Max/min Damage:', MAXVAL(Damage(:)),'/', MINVAL( Damage(:))
CALL Info( 'DamageInitialization', Message, Level=1 )
WRITE ( Message, * ) 'Number of constrained points (lower/upper): ', COUNT(ActiveNode(:,1)), COUNT(ActiveNode(:,2))

CALL Info( 'DamageInitialization', Message, Level=1 )

IF ((.NOT.Converged) .AND. (NonlinearIter > 1)) THEN
    WRITE( Message, * ) 'Nonlinear solution has not converged','Relative Change=',RelativeChange,'>',NonlinearTol
    CALL Warn('DamageInitialization', Message)
ELSE
    WRITE( Message, * ) 'Nonlinear solution has converged after ', iter,' steps.'
    CALL Info('DamageInitialization',Message,Level=1)
END IF

!----------------------------------------------------------------
!internal subroutines of DamageInitialization Solver
!----------------------------------------------------------------
CONTAINS

!----------------------------------------------------------------
SUBROUTINE LocalMatrix(MASS, STIFF, FORCE, LOAD, Element, n, TransientSimulation)

IMPLICIT NONE
!----------------------------------------------------------------
REAL(KIND=dp), DIMENSION(:,:) :: MASS, STIFF
REAL(KIND=dp), DIMENSION(:) :: FORCE, LOAD
INTEGER :: n
TYPE(Element_t), POINTER :: Element
LOGICAL :: TransientSimulation
!----------------------------------------------------------------
REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3)
REAL(KIND=dp) :: detJ, LoadAtIP
LOGICAL :: Stat, getSecondDerivatives
INTEGER :: t,i,j,DIM
TYPE(GaussIntegrationPoints_t) :: IP
TYPE(Nodes_t) :: Nodes
SAVE Nodes


DIM = CoordinateSystemDimension()
CALL GetElementNodes( Nodes )

STIFF = 0.0d0
FORCE = 0.0d0
MASS = 0.0d0

!Numerical integration:
!----------------------
IP = GaussPoints( Element )

!----------------------------------------------------------------
!  Loop over Gauss-points (element Integration)
!----------------------------------------------------------------
DO t=1,IP % n
    !Basis function values & derivatives at the integration point:
    !-------------------------------------------------------------
    getSecondDerivatives = .FALSE.
    stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
        IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, getSecondDerivatives)

    !The source term at the integration point:
    !-----------------------------------------
    LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )


    DO j=1,n
        FORCE(j) = FORCE(j) + IP % s(t) * DetJ * LoadAtIP * Basis(j)
        DO i=1,n

            IF (TransientSimulation) THEN
                MASS(i,j) = MASS(i,j)+ IP % s(t) * DetJ * Basis(i)*Basis(j)
            END IF
            STIFF(i,j) = STIFF(i,j) + IP % s(t) * DetJ * SUM(dBasisdx(i,1:DIM) * dBasisdx(j,1:DIM))
        END DO
    END DO

!----------------------------------------------------------------
END DO ! end Loop over Gauss-points (element Integration)
!----------------------------------------------------------------
END SUBROUTINE LocalMatrix
!----------------------------------------------------------------



!----------------------------------------------------------------
SUBROUTINE BoundaryCondition(LOAD, FORCE, Element, n)

IMPLICIT NONE
!----------------------------------------------------------------
REAL(KIND=dp), DIMENSION(:) :: FORCE, LOAD
INTEGER :: n
TYPE(Element_t), POINTER :: Element
!----------------------------------------------------------------
REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3)
REAL(KIND=dp) :: detJ, LoadAtIP
LOGICAL :: stat, getSecondDerivatives
INTEGER :: t,j
TYPE(GaussIntegrationPoints_t) :: IP
TYPE(Nodes_t) :: Nodes
SAVE Nodes
!----------------------------------------------------------------

CALL GetElementNodes( Nodes )
FORCE = 0.0d0
!Numerical integration:
!----------------------
IP = GaussPoints( Element )

!-----------------------------------------------------------------
! Loop over Gauss-points (boundary element Integration)
!-----------------------------------------------------------------
DO t=1,IP % n
    !Basis function values & derivatives at the integration point:
    !-------------------------------------------------------------
    getSecondDerivatives = .FALSE.
    stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
    IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, getSecondDerivatives)
    !The source term at the integration point:
    !----------------------
    LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )
    DO j=1,n
        FORCE(j) = FORCE(j) + IP % s(t) * DetJ * LoadAtIP * Basis(j)
    END DO
END DO
!----------------------------------------------------------------
END SUBROUTINE BoundaryCondition
!----------------------------------------------------------------

!-----------------------------------------------------------------
END SUBROUTINE DamageInitializationSolver
!-----------------------------------------------------------------
