program test_test

  integer*4 debug,VIsDouble,FileType,I,III
 
  double precision X1(4),Y1(4),P1(4),X2(4),Y2(4),P2(4)

  INTEGER*4 ICellMax                 / 0/
  INTEGER*4 JCellMax                 / 0/
  INTEGER*4 KCellMax                 / 0/
  INTEGER*4 DIsDouble                / 1/
  double precision  SolTime    / 360.0/
  INTEGER*4 StrandID                 / 0/      !/* StaticZone */
  INTEGER*4 ParentZn                 / 0/
  INTEGER*4 IsBlock                  / 1/      !/* Block */
  INTEGER*4 NFConns                  / 0/
  INTEGER*4 FNMode                   / 0/
  INTEGER*4 TotalNumFaceNodes        / 1/
  INTEGER*4 TotalNumBndryFaces       / 1/
  INTEGER*4 TotalNumBndryConnections / 1/
  INTEGER*4 ShrConn                  / 0/
  INTEGER*4 PassiveVarList(3) /0,0,0/
  INTEGER*4 ValueLocation(3) /1,1,1/
  INTEGER*4 ShareVarFromZone(3)/ 0,0,0/
  !*Ordered Zone Parameters*/
  INTEGER*4 IMax /2/
  INTEGER*4 JMax /2/
  INTEGER*4 KMax /1/
!  Ordered Zone */
  INTEGER*4 ZoneType /0/

  Debug = 0
  FileType = 0
  I = 0

  I = TECINI112("IJ Ordered Zones"//char(0),"X Y P"//char(0), &
       "ij_ordered.plt"//char(0),"."//char(0),FileType,Debug,VIsDouble)

  X1(1) = .125
  Y1(1) = .5
  P1(1) = 5

  X1(2) = .625
  Y1(2) = .5
  P1(2) = 7.5

  X1(3) = .125
  Y1(3) = .875
  P1(3) = 10

  X1(4) = .625
  Y1(4) = .875
  P1(4) = 7.5

  X2(1) = .375
  Y2(1) = .125
  P2(1) = 5

  X2(2) = .875
  Y2(2) = .125
  P2(2) = 7.5

  X2(3) = .375
  Y2(3) = .5
  P2(3) = 10

  X2(4) = .875
  Y2(4) = .5
  P2(4) = 7.5

  I = TECZNE112("Ordered Zone"//char(0),ZoneType,IMax,JMax,KMax,ICellMax,&
       JCellMax,KCellMax,SolTime,StrandID,ParentZn,IsBlock,NFConns,FNMode,&
       TotalNumFaceNodes,TotalNumBndryFaces,TotalNumBndryConnections,&
       PassiveVarList, ValueLocation, ShareVarFromZone,ShrConn)
 
  III = IMax*JMax*KMax
  I   = TECDAT112(III, X1, DIsDouble)
  I   = TECDAT112(III, Y1, DIsDouble)
  I   = TECDAT112(III, P1, DIsDouble)

  I = TECZNE112("Ordered Zone2"//char(0),ZoneType,IMax,JMax,KMax,ICellMax,&
       JCellMax,KCellMax,SolTime,StrandID,ParentZn,IsBlock,NFConns,FNMode,&
       TotalNumFaceNodes,TotalNumBndryFaces,TotalNumBndryConnections,&
       PassiveVarList, ValueLocation, ShareVarFromZone,ShrConn)

  I   = TECDAT112(III, X2, DIsDouble)
  I   = TECDAT112(III, Y2, DIsDouble)
  I   = TECDAT112(III, P2, DIsDouble)

  I = TECEND112()
end program test_test
