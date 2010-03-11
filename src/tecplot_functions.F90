subroutine i_ordered(name,data,n,ndim)
#IFDEF USE_TECIO
  implicit none
  integer         , intent(in) :: n,ndim
  character*(*)   , intent(in) :: name
  double precision, intent(in) :: data(n,ndim)

  ! Working 
  character(80)                :: var_names 

  INTEGER ZoneType /0/
  integer I,idim

  ! These Never Change
  integer VIsDouble                / 1/
  INTEGER ICellMax                 / 0/
  INTEGER JCellMax                 / 0/
  INTEGER KCellMax                 / 0/
  INTEGER DIsDouble                / 1/
  double precision  SolTime      / 0.0/
  INTEGER StrandID                 / 0/      !/* StaticZone */
  INTEGER ParentZn                 / 0/
  INTEGER IsBlock                  / 1/      !/* Block */
  INTEGER NFConns                  / 0/
  INTEGER FNMode                   / 0/
  INTEGER TotalNumFaceNodes        / 0/
  INTEGER TotalNumBndryFaces       / 0/
  INTEGER TotalNumBndryConnections / 0/
  INTEGER ShrConn                  / 0/
  INTEGER PassiveVarList(ndim)
  INTEGER ValueLocation(ndim)
  INTEGER ShareVarFromZone(ndim)
  INTEGER TECDAT112,TECZNE112
  PassiveVarList(:) = 0
  ValueLocation(:)  = 1
  ShareVarFromZone(:) = 0

  I = TECZNE112(name,ZoneType,n,1,1,ICellMax,JCellMax,KCellMax,&
       SolTime,StrandID,ParentZn,IsBlock,NFConns,FNMode,&
       TotalNumFaceNodes,TotalNumBndryFaces,TotalNumBndryConnections,&
       PassiveVarList, ValueLocation, ShareVarFromZone,ShrConn)
  do idim=1,ndim
     I   = TECDAT112(n, data(:,idim), DIsDouble)
  end do
#ENDIF
end subroutine i_ordered

subroutine ij_ordered(name,data,n,m,ndim)
#IFDEF USE_TECIO
  implicit none
  integer         , intent(in) :: n,m,ndim
  character*(*)   , intent(in) :: name
  double precision, intent(in) :: data(n,m,ndim)

  ! Working 
  character(80)                :: var_names 

  INTEGER ZoneType /0/
  integer I,idim

  ! These Never Change
  integer VIsDouble                / 1/
  INTEGER ICellMax                 / 0/
  INTEGER JCellMax                 / 0/
  INTEGER KCellMax                 / 0/
  INTEGER DIsDouble                / 1/
  double precision  SolTime      / 0.0/
  INTEGER StrandID                 / 0/      !/* StaticZone */
  INTEGER ParentZn                 / 0/
  INTEGER IsBlock                  / 1/      !/* Block */
  INTEGER NFConns                  / 0/
  INTEGER FNMode                   / 0/
  INTEGER TotalNumFaceNodes        / 0/
  INTEGER TotalNumBndryFaces       / 0/
  INTEGER TotalNumBndryConnections / 0/
  INTEGER ShrConn                  / 0/
  INTEGER PassiveVarList(ndim)
  INTEGER ValueLocation(ndim)
  INTEGER ShareVarFromZone(ndim)
  INTEGER TECDAT112,TECZNE112

  PassiveVarList(:) = 0
  ValueLocation(:)  = 1
  ShareVarFromZone(:) = 0

  I = TECZNE112(name,ZoneType,n,m,1,ICellMax,JCellMax,KCellMax,&
       SolTime,StrandID,ParentZn,IsBlock,NFConns,FNMode,&
       TotalNumFaceNodes,TotalNumBndryFaces,TotalNumBndryConnections,&
       PassiveVarList, ValueLocation, ShareVarFromZone,ShrConn)
  do idim=1,ndim
     I   = TECDAT112(n*m, data(:,:,idim), DIsDouble)
  end do
#ENDIF
end subroutine ij_ordered

subroutine ijk_ordered(name,data,n,m,l,ndim)
#IFDEF USE_TECIO
  implicit none
  integer         , intent(in) :: n,m,l,ndim
  character*(*)   , intent(in) :: name
  double precision, intent(in) :: data(n,m,l,ndim)

  ! Working 
  character(80)                :: var_names 

  INTEGER ZoneType /0/
  integer I,idim

  ! These Never Change
  integer VIsDouble                / 1/
  INTEGER ICellMax                 / 0/
  INTEGER JCellMax                 / 0/
  INTEGER KCellMax                 / 0/
  INTEGER DIsDouble                / 1/
  double precision  SolTime      / 0.0/
  INTEGER StrandID                 / 0/      !/* StaticZone */
  INTEGER ParentZn                 / 0/
  INTEGER IsBlock                  / 1/      !/* Block */
  INTEGER NFConns                  / 0/
  INTEGER FNMode                   / 0/
  INTEGER TotalNumFaceNodes        / 0/
  INTEGER TotalNumBndryFaces       / 0/
  INTEGER TotalNumBndryConnections / 0/
  INTEGER ShrConn                  / 0/
  INTEGER PassiveVarList(ndim)
  INTEGER ValueLocation(ndim)
  INTEGER ShareVarFromZone(ndim)
  INTEGER TECDAT112,TECZNE112

  PassiveVarList(:) = 0
  ValueLocation(:)  = 1
  ShareVarFromZone(:) = 0

  I = TECZNE112(name,ZoneType,n,m,l,ICellMax,JCellMax,KCellMax,&
       SolTime,StrandID,ParentZn,IsBlock,NFConns,FNMode,&
       TotalNumFaceNodes,TotalNumBndryFaces,TotalNumBndryConnections,&
       PassiveVarList, ValueLocation, ShareVarFromZone,ShrConn)
  do idim=1,ndim
     I   = TECDAT112(n*m, data(:,:,:,idim), DIsDouble)
  end do
#ENDIF
end subroutine ijk_ordered

subroutine open_tecplot(fname,ndim)
 
#IFDEF USE_TECIO

  character*(*)   , intent(in) :: fname
  integer         , intent(in) :: ndim
  ! Working 
  character(80)                :: var_names 
  ! These Never Change
  integer VIsDouble                /1/
  integer debug                    /1/
  integer FileType                 /0/ 
  INTEGER ZoneType                 /0/
  integer I
  INTEGER TECINI112
  if (ndim == 1) then
     var_names = "X"//char(0)
  else if (ndim == 2) then
     var_names = "X,Y"//char(0)
  else
     var_names = "X,Y,Z"//char(0)
  end if

  I = TECINI112("pyPSG Data"//char(0),var_names,fname,"."//char(0),&
       FileType,Debug,VIsDouble)
#ENDIF
end subroutine open_tecplot

subroutine close_tecplot()
  integer I,TECEND112
#IFDEF USE_TECIO
  I = TECEND112()
#ENDIF
end subroutine close_tecplot

subroutine tecplot_test(available)
  
  integer, intent(out) :: available
  
#IFDEF USE_TECIO
  available = 1
#ELSE
  available = 0
#ENDIF

end subroutine tecplot_test
