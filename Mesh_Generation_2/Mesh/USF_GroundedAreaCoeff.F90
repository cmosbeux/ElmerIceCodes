!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  A USF that can be used with SaveScalars to get the grounded area;
!   tested in 2D plane view applications with SSA
!
!  SaveScalars:
!      Variable 5 = "groundedmask"
!      Operator 5 = "volume"
!      Coefficient 5 = "GroundedAreaCoeff"
!  Material:
!       GroundedAreaCoeff = Variable groundedmask
!         real procedure "USF_GroundedAreaCoeff" "GroundedAreaCoeff" 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION GroundedAreaCoeff(Model,nodenumber,VarIn) RESULT(VarOut)
       USE DefUtils
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn,VarOut

       TYPE(Element_t),POINTER :: Element
       REAL(KIND=dp), ALLOCATABLE,SAVE :: LocalGL(:)
       LOGICAL,SAVE:: AllocationDone=.FALSE.
       INTEGER :: n

       IF (.NOT.AllocationDone) THEN 
         allocate(LocalGL(Model % Mesh % MaxElementNodes))
         AllocationDone=.TRUE.
       END IF

       Element => GetCurrentElement()
       n=Element%TYPE%NumberOfNodes
       CALL GetLocalSolution(LocalGL,'groundedmask',Element)
       IF (ALL(LocalGL(1:n).GT.-0.5)) THEN
          VarOut=1._dp
       ELSE
          VarOut=0._dp
       END IF

       End FUNCTION GroundedAreaCoeff
