SUBROUTINE nearestpoint( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
!USE CoordinateSystems
!USE MeshUtils
USE DefUtils
USE Netcdf


IMPLICIT NONE
!------------------------------------------------------------------------------
TYPE(Solver_t), TARGET :: Solver
TYPE(Model_t) :: Model
REAL(KIND=dp) :: dt
LOGICAL :: Transient

TYPE(ValueList_t), POINTER :: Params
TYPE(Variable_t), POINTER :: Var
TYPE(Nodes_t) :: ElementNodes
TYPE(Element_t),POINTER ::  Element
REAL(KIND=dp), POINTER :: Values(:)
INTEGER, POINTER :: Perm(:), NodeIndexes(:), Indexxx


CHARACTER(LEN=MAX_NAME_LEN) :: TargetVariableName,VariableName,DataF
CHARACTER(LEN=MAX_NAME_LEN) :: Xdim,dimName,FillName,Name,variabletype
CHARACTER(LEN=MAX_NAME_LEN) :: FName
CHARACTER(LEN=MAX_NAME_LEN),parameter :: &
                         SolverName='nearestpoint'
LOGICAL :: GotVar,GotTVar,Found,UnFoundFatal=.TRUE.
LOGICAL :: CheckBBox, NETCDFFormat
LOGICAL :: HaveFillv


INTEGER :: NetcdfStatus,varid,ncid
INTEGER :: NoVar,k,e
INTEGER :: nx,ny,n,node
INTEGER :: XMinIndex,XMaxIndex
INTEGER :: YMinIndex,YMaxIndex, Indexx, Indexy


REAL(KIND=dp) ::  Xb, Yb, dx, dy
REAL(KIND=DP),allocatable :: xx(:),yy(:),DEM(:,:)
REAL(KIND=DP) :: fillv
REAL(KIND=DP) :: Val
REAL(KIND=DP) :: xmin,xmax,ymin,ymax

CALL INFO(Trim(SolverName),'start solver', Level=5)

Params => GetSolverParams()

NoVar=0
GotVar=.True.

DO WHILE(GotVar)
   NoVar = NoVar + 1
   ! netcdf reading
   
   WRITE(Name,'(A,I0)') 'Variable ',NoVar

   
   VariableName = ListGetString( Params, TRIM(Name), GotVar )
   IF (.NOT.GotVar) exit

   WRITE(Name,'(A,I0)') 'Target Variable ',NoVar
   TargetVariableName=ListGetString( Params, TRIM(Name), GotTVar)
   IF (.NOT.GotTVar) TargetVariableName=VariableName

   Var => VariableGet(Model %  Mesh % Variables, TargetVariableName )
   IF(.NOT.ASSOCIATED(Var)) Then
      CALL VariableAddVector(Model % Mesh % Variables,Model % Mesh,Solver,TargetVariableName,1)
      Var => VariableGet(Model %  Mesh % Variables, TargetVariableName )
   ENDIF
   Values => Var % Values
   Perm => Var % Perm

   WRITE (FName,'(A,I0,A)') 'Variable ',NoVar,' Data File'
      
   DataF = ListGetString( Params, TRIM(FName), Found, UnFoundFatal )
   k = INDEX( DataF,'.nc' )
   NETCDFFormat = ( k /= 0 )

   IF (NETCDFFormat) then
      CALL INFO(Trim(SolverName),'Data File is in netcdf format', Level=5)
   Else
      CALL INFO(Trim(SolverName),'Data File is in ascii', Level=5)
   Endif

   NetCDFstatus = NF90_OPEN(trim(DataF),NF90_NOWRITE,ncid)
   IF ( NetCDFstatus /= NF90_NOERR ) THEN
      CALL Fatal(Trim(SolverName), &
           'Unable to open NETCDF File')
   END IF
   WRITE (dimName,'(A,I0,A)') 'Variable ',&
        NoVar,' x-dim Name'
   Xdim=ListGetString( Params, TRIM(dimName), Found )
   if (.NOT.Found) Xdim='x'
   NetCDFstatus = nf90_inq_dimid(ncid, trim(Xdim) , varid)
   IF ( NetCDFstatus /= NF90_NOERR ) THEN
      CALL Fatal(Trim(SolverName), &
           'Unable to  get netcdf x-dim Id')
   ENDIF
   NetCDFstatus = nf90_inquire_dimension(ncid,varid,len=nx)
   IF ( NetCDFstatus /= NF90_NOERR ) THEN
      CALL Fatal(Trim(SolverName), &
           'Unable to  get netcdf nx')
   ENDIF
   WRITE (dimName,'(A,I0,A)') 'Variable ',&
        NoVar,' y-dim Name'
   Xdim=ListGetString( Params, TRIM(dimName), Found )
   if (.NOT.Found) Xdim='y'
   NetCDFstatus = nf90_inq_dimid(ncid, trim(Xdim) , varid)
   IF ( NetCDFstatus /= NF90_NOERR ) THEN
      CALL Fatal(Trim(SolverName), &
           'Unable to  get netcdf y-dim Id')
   ENDIF
   NetCDFstatus = nf90_inquire_dimension(ncid,varid,len=ny)
   IF ( NetCDFstatus /= NF90_NOERR ) THEN
      CALL Fatal(Trim(SolverName), &
           'Unable to  get netcdf ny')
   ENDIF

   !! allocate good size
   allocate(xx(nx),yy(ny))

   !! Get X variable
   WRITE (dimName,'(A,I0,A)') 'Variable ',&
        NoVar,'x-Var Name'
   Xdim=ListGetString( Params, TRIM(dimName), Found )
   if (.NOT.Found) Xdim='x'
   NetCDFstatus = nf90_inq_varid(ncid,trim(Xdim),varid)
   IF ( NetCDFstatus /= NF90_NOERR ) THEN
      CALL Fatal(Trim(SolverName), &
           'Unable to get netcdf x-variable id')
   ENDIF
   NetCDFstatus = nf90_get_var(ncid, varid,xx)
   IF ( NetCDFstatus /= NF90_NOERR ) THEN
      CALL Fatal(Trim(SolverName), &
           'Unable to get netcdf x-variable ')
   ENDIF
   !! Get Y variable
   WRITE (dimName,'(A,I0,A)') 'Variable ',&
        NoVar,'y-Var Name'
   Xdim=ListGetString( Params, TRIM(dimName), Found )
   if (.NOT.Found) Xdim='y'
   NetCDFstatus = nf90_inq_varid(ncid,trim(Xdim),varid)
   IF ( NetCDFstatus /= NF90_NOERR ) THEN
      CALL Fatal(Trim(SolverName), &
           'Unable to get netcdf y-variable id')
   ENDIF
   NetCDFstatus = nf90_get_var(ncid, varid,yy)
   IF ( NetCDFstatus /= NF90_NOERR ) THEN
      CALL Fatal(Trim(SolverName), &
           'Unable to get netcdf y-variable')
   ENDIF

   !! Check that there is data within the domain
   IF ((MAXVAL(xx).LT.xmin).OR.(MINVAL(xx).GT.xmax)&
        .OR.(MAXVAL(yy).LT.ymin).OR.(MINVAL(yy).GT.ymax)) &
        CALL Fatal(Trim(SolverName), &
        'No data within model domain')

   XMinIndex = 1
   XMaxIndex = nx
   YMinIndex = 1
   YMaxIndex = ny

   write(message,'(A,I0,A,I0,A)') 'NETCDF: reading nx=',nx,&
        ' and ny=',ny,' data points'
   CALL INFO(Trim(SolverName),Trim(message),Level=5)
   write(message,*) 'X Indexes: ',&
        XMinIndex,XMaxIndex,xx(XMinIndex),xx(XMaxIndex)
   CALL INFO(Trim(SolverName),Trim(message),Level=10)
   write(message,*) 'Y Indexes: ', &
        YMinIndex,YMaxIndex,yy(YMinIndex),yy(YMaxIndex)
   CALL INFO(Trim(SolverName),Trim(message),Level=10)

   dx=(xx(XMaxIndex)-xx(XMinIndex))/(XMaxIndex-1)
   dy=(yy(YMaxIndex)-yy(YMinIndex))/(YMaxIndex-1)
   xx=xx-dx/2
   yy=yy-dy/2
   
   allocate(DEM(nx,ny))

   !! Get the variable
   NetCDFstatus = nf90_inq_varid(ncid,TRIM(VariableName),varid)
   IF ( NetCDFstatus /= NF90_NOERR ) THEN
      CALL Fatal(Trim(SolverName), &
           'Unable to get netcdf variable id')
   ENDIF
   NetCDFstatus = nf90_get_var(ncid, varid,DEM(:,:),&
        start = (/ XMinIndex, YMinIndex /),     &
        count = (/ nx,ny/))
   IF ( NetCDFstatus /= NF90_NOERR ) THEN
      CALL Fatal(Trim(SolverName), &
           'Unable to get netcdf variable')
   ENDIF
   HaveFillV=.True.
   NetCDFstatus = nf90_get_att(ncid, varid,"_FillValue",fillv)
   IF ( NetCDFstatus /= NF90_NOERR ) THEN
      HaveFillV=.False.
   ENDIF
   WRITE (FillName,'(A,I0,A)') 'Variable ',NoVar,' Fill Value'
   Val=ListGetConstReal(Params, TRIM(FillName) , Found)
   IF (Found) THEN
      HaveFillV=.TRUE.
      fillv=Val
   ENDIF

   !! Close NETCDF
   NetCDFstatus = nf90_close(ncid)   
 
   ! loop on the elements
   DO e=1,Solver % NumberOfActiveElements
      Element => GetActiveElement(e)
      CALL GetElementNodes( ElementNodes )
      n = GetElementNOFNodes(Element)
      SELECT CASE ((Var % TYPE))

      CASE(Variable_on_elements)! case element       
       !print *,'coord', Xcoord,Ycoord
         ! barycentre computing
       Indexxx => Element % ElementIndex
       Xb = SUM(ElementNodes % x(1:n))/n
       Yb = SUM(ElementNodes % y(1:n))/n
      
       Indexx=CEILING((Xb-xx(XMinIndex))/dx+1)
       Indexy=CEILING((Yb-yy(YMinIndex))/dy+1)
       Values(Perm(Indexxx))=NINT(DEM(Indexx,Indexy)) ! valeur que je donne à l'élément

      CASE (Variable_on_nodes) 
        NodeIndexes => Element % NodeIndexes
        Do node=1,n
          Xb =  ElementNodes % x(node)
          Yb = ElementNodes % y(node)

          Indexx=CEILING((Xb-xx(XMinIndex))/dx+1)
          Indexy=CEILING((Yb-yy(YMinIndex))/dy+1)
          !print *, Indexx 

          Values(Perm(NodeIndexes(node)))=INT(DEM(Indexx,Indexy)) ! valeur que je donne au noeud
        end Do
      END SELECT
   end DO
   
end DO

CALL INFO(Trim(SolverName), &
     '-----ALL DONE----------',Level=5)

end SUBROUTINE nearestpoint
