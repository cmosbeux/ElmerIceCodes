      FUNCTION DistanceCondGL(Model,nodenumber,VarIn) RESULT(VarOut)
      USE DefUtils
      implicit none
      !-----------------
      TYPE(Model_t) :: Model
      INTEGER :: nodenumber
      REAL(kind=dp) :: VarIn,VarOut

      IF ((VarIn.LT.0.5).AND.(VarIn.GT.-0.5)) THEN
        VarOut = +1.0
      ELSE
        VarOut = -1.0
      END IF
      End FUNCTION DistanceCondGL

      FUNCTION dfunction(Model,nodenumber,VarIn) RESULT(VarOut)
      USE DefUtils
      implicit none
      !-----------------
      TYPE(Model_t) :: Model
      INTEGER :: nodenumber
      REAL(kind=dp):: VarIn,VarOut

      REAL(kind=dp),parameter :: c=2._dp/9._dp
      REAL(kind=dp),SAVE :: a,b
      REAL(kind=dp) :: a2,d

      LOGICAL,SAVE :: FirstVisit=.TRUE.

      IF (FirstVisit) THEN
        a=ListGetConstReal( Model % Constants,'GL mesh grad',UnFoundFatal=.TRUE.)
        b=ListGetConstReal( Model % Constants,'GL mesh size',UnFoundFatal=.TRUE.)
      END IF

       a2=a**2
       d=a*VarIn+b

       VarOut=-log(d)/(c*a2)

      END FUNCTION dfunction

      FUNCTION Hmin_d(Model,nodenumber,VarIn) RESULT(VarOut)
      USE DefUtils
      implicit none
      !-----------------
      TYPE(Model_t) :: Model
      INTEGER :: nodenumber
      REAL(kind=dp):: VarIn,VarOut
      REAL(kind=dp),SAVE :: Extent,HMinIn,HMinOut
      
       LOGICAL,SAVE :: FirstVisit=.TRUE.

      IF (FirstVisit) THEN
        Extent=ListGetConstReal( Model % Constants,'Hmin Extent',UnFoundFatal=.TRUE.)
        HMinIn=ListGetConstReal( Model % Constants,'Hmin in',UnFoundFatal=.TRUE.)
        HMinOut=ListGetConstReal( Model % Constants,'Hmin out',UnFoundFatal=.TRUE.)
      END IF

      IF (VarIn.LT.Extent) THEN
         VarOut = HMinIn
      Else
         Varout = HMinOut
      End if

      END 
