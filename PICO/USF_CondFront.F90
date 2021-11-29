FUNCTION CondFront(Model,nodenumber,VarIn) RESULT(VarOut)

USE types
USE CoordinateSystems
USE SolverUtils
USE ElementDescription
USE DefUtils

IMPLICIT NONE

TYPE(Model_t) :: Model
INTEGER :: nodenumber
REAL(kind=dp) :: VarIn
REAL(kind=dp) :: VarOut

REAL(kind=dp),parameter :: cond=0.5_dp

IF (VarIn .GE. cond) THEN
  VarOut = 1.0_dp !is front
ELSE
  VarOut = -1.0_dp !is not front
ENDIF

END FUNCTION CondFront
