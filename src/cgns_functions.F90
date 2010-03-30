subroutine cgns_test()

#ifdef USE_CGNS 
  print *,'CGNS enabled in this build'
#else
  print *,'CGNS not enabled in this build'
#endif

end subroutine cgns_test

subroutine open_cgns(filename,cg,nzones)

#ifdef USE_CGNS

  implicit none
  include 'cgnslib_f.h'

  character*32,intent(in) :: filename
  integer , intent(out) ::cg,nzones

  character*32 :: name
  integer base,j,ier,nbases
  character*32 basename
  integer  CellDim, PhysDim

  call cg_open_f(filename, MODE_READ, cg, ier)
  if (ier .eq. ERROR) call cg_error_exit_f

  call cg_nbases_f(cg, nbases, ier)
  if (ier .eq. ERROR) call cg_error_exit_f

  if (nbases .gt. 1) then
     print *, 'Warning: This program only reads the first base in a cgns file'
  end if

  base = 1

  call cg_base_read_f(cg, base, basename, CellDim, PhysDim, ier)
  if (ier .eq. ERROR) call cg_error_exit_f
  if (cellDim .ne. 3 .or. PhysDim .ne. 3) then
     print *,'The Cells must be hexahedreal in 3 dimensions'
     stop
  end if
  !       *** base attribute:  GOTO base node
  call cg_goto_f(cg, base, ier, 'end')
  if (ier .eq. ERROR) call cg_error_exit_f
  call cg_nzones_f(cg, base, nzones, ier)
  if (ier .eq. ERROR) call cg_error_exit_f

#endif
end subroutine open_cgns

subroutine read_cgns_zone_shape(cg,zone,zoneshape)
#ifdef USE_CGNS
  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg,zone
  integer, intent(out) :: zoneshape(3)

  ! Working
  integer base,zonesize(6),ier
  character*32, zonename

  base = 1
  call cg_zone_read_f(cg, base, zone, zonename, zonesize, ier)
  zoneshape(1:3) = zonesize(1:3)

#endif
end subroutine read_cgns_zone_shape


subroutine read_cgns_zone(cg,zone,X,nx,ny,nz,faceBCs)
#ifdef USE_CGNS
  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg,zone,nx,ny,nz
  double precision, intent(out):: X(nx,ny,nz,3)
  integer, intent(out) :: faceBCs(6)
  ! Working
  integer base,zonesize(6),ier,zonetype,start(3),upper(3),nbocos,boco
  character*32, zonename,boconame
  integer ptset_type,normalIndex(3),NormalListFlag,datatype,ndataset
  integer faceID,points(6),bocotype,npts
  double precision data_double(6)
  ! Read the coordinates
  base = 1
  start(1) = 1
  start(2) = 1
  start(3) = 1
  upper(1) = nx
  upper(2) = ny
  upper(3) = nz
  faceBcs(:) = 0
  call cg_coord_read_f(cg,base,zone,'CoordinateX',RealDouble,start,upper,X(:,:,:,1),ier)
  call cg_coord_read_f(cg,base,zone,'CoordinateY',RealDouble,start,upper,X(:,:,:,2),ier)
  call cg_coord_read_f(cg,base,zone,'CoordinateZ',RealDouble,start,upper,X(:,:,:,3),ier)

  ! Get the boundary conditions

  call cg_nbocos_f(cg, base, zone, nbocos, ier)
  if (ier .eq. ERROR) call cg_error_exit_f

  do boco=1, nbocos
     call cg_boco_info_f(cg, base, zone, boco, boconame,bocotype,&
          ptset_type,npts,NormalIndex,NormalListFlag,datatype,ndataset,ier)
     if (ier .eq. ERROR) call cg_error_exit_f

     call cg_boco_read_f(cg, base, zone, boco, points,data_double, ier)
     if (BCTypeName(bocotype) == 'BCWallViscous' .or. &
          BCTypeName(bocotype) == 'BCInflow' .or. &
          BCTypeName(bocotype) == 'BCOutflow' .or. &
          BCTypeName(bocotype) == 'BCFarfield') then
        
        ! Make sure its a SURFACE BC
        if ( (points(4)-points(1) > 0 .and. points(5)-points(2) > 0 ) .or. &
             (points(4)-points(1) > 0 .and. points(6)-points(3) > 0 ) .or. &
             (points(5)-points(2) > 0 .and. points(6)-points(3) > 0 )) then

           if (points(6) == points(3) .and. points(3) == 1) then ! Face 1 - Zmin
              faceID = 1
           else if (points(6) == points(3) .and. points(3) == nz) then ! face 2 Zmax
              faceID = 2
           else if (points(4) == points(1) .and. points(4) == 1) then ! Face 3 Imin
              faceID = 3
           else if (points(4) == points(1) .and. points(4) == nx) then ! Face 4 I max
              faceID = 4
           else if (points(5) == points(2) .and. points(5) == 1) then ! Face 5 J min
              faceID = 5
           else if (points(5) == points(2) .and. points(5) == ny) then ! Face 6 J max
              faceID = 6
           end if
           
           if (BCTypeName(bocotype) == 'BCWallViscous' ) then
              faceBCs(faceID) = 1
           end if
           if (BCTypeName(bocotype) == 'BCInflow' ) then
              faceBCs(faceID) = 2
           end if
           if (BCTypeName(bocotype) == 'BCOutflow' ) then
              faceBCs(faceID) = 2
           end if
           if (BCTypeName(bocotype) == 'BCFarfield' ) then
              faceBCs(faceID) = 2
           end if
           
        end if
     end if
  end do
#endif
   end subroutine read_cgns_zone

   subroutine close_cgns(cg)
#ifdef USE_CGNS
     implicit none
     include 'cgnslib_f.h'
     integer, intent(in) :: cg
     integer ier
     call cg_close_f(cg, ier)
     if (ier .eq. ERROR) call cg_error_exit_f
#endif
   end subroutine close_cgns
