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
LOGICAL :: AllocationsDone = .FALSE., Found, Converged
INTEGER :: n, t, istat, other_body_id, iter, NonlinearIter
REAL(KIND=dp) :: Norm, PrevNorm=0.0d00, NonlinearTol, RelativeChange
TYPE(ValueList_t), POINTER :: BodyForce, Material, BC, SolverParams
REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), LOAD(:), FORCE(:)

CHARACTER(LEN=MAX_NAME_LEN) :: BoundaryType
SAVE MASS, STIFF, LOAD, FORCE,&
    AllocationsDone, PrevNorm

! Allocate some permanent storage, this is done first time
!--------------------------------------------------------------
IF ( .NOT. AllocationsDone ) THEN
    N = Solver % Mesh % MaxElementNodes !big enough for elemental arrays
    ALLOCATE( FORCE(N), LOAD(N), MASS(N,N), STIFF(N,N), STAT=istat )
    IF ( istat /= 0 ) THEN
        CALL Fatal( 'MyHeatSolve','Memory allocation error for matrix/vectors.' )
    END IF
    AllocationsDone = .TRUE.
END IF

!Read in solver parameters
!-------------------------
SolverParams => GetSolverParams()
IF (.NOT. ASSOCIATED(SolverParams)) &
    CALL FATAL('Damage initialization','No Solver section found')

NonlinearIter = GetInteger(SolverParams, & 'Nonlinear System Max Iterations', Found)
IF ( .NOT.Found ) NonlinearIter = 1
NonlinearTol = GetConstReal( SolverParams, 'Nonlinear System Convergence Tolerance', Found )
IF ( .NOT.Found ) NonlinearTol = 1.0D-03


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
     DO t=1,Solver % NumberOfActiveElements
        ! get element info
        !-----------------
        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        ! get material parameters
        !----------------------
        Material => GetMaterial()
        IF (.NOT. ASSOCIATED(Material)) THEN
        WRITE(Message,'(A,I5,A)') 'No material for bulk element no. ',t,' found.'
           CALL FATAL('DamageInitialization',Message)
        END IF

        !Get load for force vector !------------------------- LOAD = 0.0d0
        BodyForce => GetBodyForce()
        IF ( ASSOCIATED(BodyForce) ) LOAD(1:n) = GetReal( BodyForce, 'Damage Source', Found )

        !Get element local matrix and rhs vector:
        !----------------------------------------
        CALL LocalMatrix( MASS, STIFF, FORCE, LOAD, Element, n,TransientSimulation)

        !Update global matrix and rhs vector from local matrix and vector
        !---------------------------------------------------------------
        IF ( TransientSimulation ) THEN
            CALL Default1stOrderTime( MASS,STIFF,FORCE )
        END IF

        CALL DefaultUpdateEquations( STIFF, FORCE )
!----------------------------------------------------------------
    END DO ! end Assembly for the domain
!----------------------------------------------------------------

!----------------------------------------------------------------
! Assembly for the Neumann boundaru conditions (Not used right now)
!----------------------------------------------------------------

    DO t=1, Solver % Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(t)
        IF ( .NOT.ActiveBoundaryElement() ) CYCLE
        n = GetElementNOFNodes()
        ! no evaluation of Neumann BC’s on points
        IF ( GetElementFamily() == 1 ) CYCLE
        BC => GetBC()

        FORCE = 0.0d00
        MASS = 0.0d00
        STIFF = 0.0d00
    
        ! check type of boundary and set BC accordingly
        !----------------------------------------------
        BoundaryType = GetString(BC,'Boundary Type',Found)
        IF (.NOT. Found) CYCLE

!----------------------------------------------------------------
    END DO ! end Assembly for Neumann boundary conditions
!----------------------------------------------------------------

    CALL DefaultFinishAssembly()
    
    ! call Elmer Solver routine for Dirichlet BCs
    !-------------------------------------------------
    CALL DefaultDirichletBCs()

    
    ! Solve the system
    ! ----------------
    Norm = DefaultSolve()

    ! compute relative change of norm
    ! -------------------------------
    IF ( PrevNorm + Norm /= 0.0d0 ) THEN
        RelativeChange = 2.0d0 * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
    ELSE
            RelativeChange = 0.0d0
    END IF
    
    WRITE( Message, * ) 'Result Norm : ',Norm
    CALL Info( 'DamageInitialization', Message, Level=4 )
    WRITE( Message, * ) 'Relative Change : ', RelativeChange
    CALL Info( 'DamageInitialization', Message, Level=4 )
    
    ! do we have to do another round?
    ! -------------------------------
    IF ( RelativeChange < NonlinearTol ) THEN ! NO
        Converged = .TRUE.
        EXIT
    ELSE ! YES
        PrevNorm = Norm
    END IF
!----------------------------------------------------------------
END DO ! of the nonlinear iteration
!----------------------------------------------------------------

! has non-linear solution converged?
! ----------------------------------
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
                MASS(i,j) = MASS(i,j)+ IP % s(t) * DetJ * & Basis(i)*Basis(j)
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
