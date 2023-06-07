!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author  : Cyrille Mosbeux
! This user function computes a parametrized melt rate as a function of the ice
! base depth and the location (as a function of X in (X,Y) stereopolar coords )
! in the Ross Ice Shelf embayment. The X dependence takes advantage of the
! differences between the West and East bathymetry and the correlated High
! Salinity Surface Water circulation


       FUNCTION Melt2Depth(Model,nodenumber,VarIn) RESULT(melt)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(3) !x, z, gmask
       REAL(kind=dp) :: x, z, gmask
       REAL(kind=dp) :: K, ri, rw, Lf, Cw, Tz1, Tz2, z1, z2
       REAL(kind=dp) :: T, slope,  melt, we_boundary

       x = VarIn(1)
       z = VarIn(2)
       gmask = VarIn(3)
       K = 15.77
       Cw = 4.218
       Lf = 335.0
       rw = 1028.
       ri = 917.
       z1 = -400.0
       z2 = -1100.0
       Tz1= 0.25
       Tz2 = 3.5

       ! groundedmask test
       melt = 0.
       IF ( gmask > 0. ) RETURN

       ! account for the boundary
       we_boundary = 75000.
       IF ( x < we_boundary ) THEN
           Tz2 = 1.0
       ENDIF
       ! account for depth
       IF ( z > z1 ) THEN
           T = Tz1
       ELSEIF ( z < z2 ) THEN
           T = Tz2
       ELSE
           ! linear interpolation
           slope = (Tz2 - Tz1) / (z2 - z1)
           T = Tz1 + z * slope
       ENDIF

       melt = 8*K*rw*Cw*T*abs(T)/(ri*Lf)


       END FUNCTION Melt2Depth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
