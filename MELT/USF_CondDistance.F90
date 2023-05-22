FUNCTION CondDistance(Model,nodenumber,VarIn) RESULT(VarOut)

USE types
USE CoordinateSystems
USE SolverUtils
USE ElementDescription
USE DefUtils

IMPLICIT NONE

TYPE(Model_t) :: Model
INTEGER :: nodenumber
REAL(kind=dp) :: VarIn !GroundedMask
REAL(kind=dp) :: VarOut

REAL(kind=dp),parameter :: cond=0.1_dp

if (abs(VarIn).GT.cond) then
  VarOut=-1.0 ! is not GL
else
  VarOut=1.0  ! is GL
endif

END FUNCTION CondDistance
