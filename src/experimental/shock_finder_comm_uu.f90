
module SurfaceData

  implicit none

  logical :: xextent=.false.
  logical :: yextent=.false.
  logical :: zextent=.false.
  integer :: x_ghost=0
  integer :: y_ghost=0
  integer :: z_ghost=0
  integer :: npoints=16
  integer :: ndimensions=2
  integer :: problem_dimensions=2
  integer :: nfaces=8
  integer :: nface_types=2
  integer, allocatable, dimension(:,:) :: points
  integer, allocatable, dimension(:,:) :: normals
  integer, allocatable, dimension(:) :: faces
  integer, allocatable, dimension(:) :: face_types
  character(len=60), allocatable, dimension(:) :: area_elements
endmodule
!***********************************************************************
program shock_finder3D
!
  use Cparam
  use SurfaceData
!
  implicit none
!
  integer :: unitno=6
  integer :: rotation=0
!
  call read_surfaceinfo

  if ( ndimensions /= problem_dimensions ) then
     print "(a,i1,a,i1,a)", "The compiled in shock profile integral is ",ndimensions, &
                     "D but the problem is ", problem_dimensions,"D. STOPPING."
     STOP
  endif

  xextent=(nxgrid/=1)
  yextent=(nygrid/=1)
  zextent=(nzgrid/=1)

  if (xextent) x_ghost=3
  if (yextent) y_ghost=3
  if (zextent) z_ghost=3
!
! Rotations of a 2D profile:
!
! 0 = xy
! 1 = xz
! 2 = yz
!
!
! Fluid body calculation
!
  call make_calc_body(unitno)
!
! Local edge calculation (awaiting contributions from neighbours)
!
  call make_calc_internalboundary(unitno)
!
! Neighbour edge calculation (contributions for neighbours)
!
  call make_calc_externalboundary(unitno)
!
! Tidy up
!
  deallocate(points)
  deallocate(normals)
  deallocate(faces)
  deallocate(face_types)
  deallocate(area_elements)
endprogram shock_finder3D
!***********************************************************************
subroutine make_calc_body(unitno)
  use Cparam
  integer :: unitno
  write(unitno,"(a)") "!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)"
  write(unitno,"(a)") "!"
  write(unitno,"(a)") "!***********************************************************************"
  write(unitno,"(a)") "  subroutine shock_calc_body(f)"
  write(unitno,"(a)") "    use Cdata"
  write(unitno,"(a)") "    use Mpicomm, only: stop_it"
  write(unitno,"(a)") "    real, dimension (mx,my,mz,mfarray) :: f"
  call declare_facefactors(unitno)
  write(unitno,"(a)") "    if (.not.lgauss_integral_comm_uu) call fatal_error( &"
  write(unitno,"(a)") "      'shock_calc_body', &"
  write(unitno,"(a)") "      'Attempted to use communicated uu shock profile with non-communicated uu code!')"
  if ((nxgrid/=1).and.(nygrid/=1).and.(nzgrid/=1)) then
    call evaluate_facefactors(unitno,0)
    call evaluate_integral(unitno,0,-3,+3,-3,+3,-3,+3,'f(l1:l2,m,n,ishock)','l1:l2','m','n')
  elseif ((nxgrid/=1).and.(nygrid/=1)) then
    call evaluate_facefactors(unitno,0)
    call evaluate_integral(unitno,0,-3,+3,-3,+3,-3,+3,'f(l1:l2,m,n1,ishock)','l1:l2','m','n1')
  elseif ((nxgrid/=1).and.(nzgrid/=1)) then
    call evaluate_facefactors(unitno,1)
    call evaluate_integral(unitno,1,-3,+3,-3,+3,-3,+3,'f(l1:l2,m1,n,ishock)','l1:l2','m1','n')
  elseif ((nzgrid/=1).and.(nygrid/=1)) then
    call evaluate_facefactors(unitno,2)
    call evaluate_integral(unitno,2,-3,+3,-3,+3,-3,+3,'f(l1,m,n,ishock)','l1','m','n')
  else
    write(unitno,"(a)") "      call fatal_error('shock_calc_body','Case not implemented')"
  endif
  write(unitno,"(a)") "    endsubroutine shock_calc_body"
endsubroutine make_calc_body
!***********************************************************************
subroutine make_calc_internalboundary(unitno)
  use Cparam
  integer :: unitno,i,j,k
  logical :: ibound, jbound, kbound
  character (len=4) :: istr='',jstr='',kstr=''
  character (len=5) :: pistr='',nistr=''
  character (len=5) :: pjstr='',njstr=''
  character (len=5) :: pkstr='',nkstr=''

  write(unitno,"(a)") "!***********************************************************************"
  write(unitno,"(a)") "  subroutine shock_calc_internalboundary(f)"
  write(unitno,"(a)") "    use Cdata"
  write(unitno,"(a)") "    use Mpicomm, only: stop_it"
  write(unitno,"(a)") "    real, dimension (mx,my,mz,mfarray) :: f"

  write(unitno,"(a)") "      call stop_it('shock_calc_internalboundary: Should not used in communicated uu shock')"
  write(unitno,"(a)") "      if (NO_WARN) print*,f"
  write(unitno,"(a)") "  endsubroutine shock_calc_internalboundary"

endsubroutine make_calc_internalboundary
!***********************************************************************
subroutine make_calc_externalboundary(unitno)
  use Cparam
  integer :: unitno,i,j,k
  integer :: imin, imax
  integer :: jmin, jmax
  integer :: kmin, kmax
  logical :: ibound, jbound, kbound
  character (len=4) :: istr='',jstr='',kstr=''
  character (len=5) :: pistr='',nistr=''
  character (len=5) :: pjstr='',njstr=''
  character (len=5) :: pkstr='',nkstr=''

  write(unitno,"(a)") "!  called: make_calc_externalboundary"
  write(unitno,"(a)") "  subroutine shock_calc_externalboundary(f)"
  write(unitno,"(a)") "    use Cdata"
  write(unitno,"(a)") "    use Mpicomm, only: stop_it"
  write(unitno,"(a)") "    real, dimension (mx,my,mz,mfarray) :: f"
  write(unitno,"(a)") "      call stop_it('shock_calc_externalboundary: Should not used in communicated uu shock')"
  write(unitno,"(a)") "      if (NO_WARN) print*,f"

  write(unitno,"(a)") "  endsubroutine shock_calc_externalboundary"
endsubroutine make_calc_externalboundary
!***********************************************************************
subroutine evaluate_facefactors(unitno,rotation)
   use SurfaceData

   integer :: rotation, unitno
   integer :: i,face_type,dir1,dir2,dir3
   character(len=2), dimension(3) :: dmesh = (/ 'dx','dy','dz' /)

   dir1=mod(3-rotation,3)+1
   dir2=mod(4-rotation,3)+1
   dir3=mod(5-rotation,3)+1
   print "(a,a)","      dA = ",dmesh(dir1)
   print "(a,a)","      dB = ",dmesh(dir2)
   print "(a,a)","      dC = ",dmesh(dir3)
   do face_type=1,nface_types
     print "(a17,i1,a)","      face_factor",face_type," = -"//trim(area_elements(face_type))//"* dxmin**2"
   enddo  !Face_types
endsubroutine evaluate_facefactors
!***********************************************************************
subroutine declare_facefactors(unitno)
!
   use SurfaceData
!
   integer :: unitno
   integer :: face_type
!
   write(unitno,"(a)") "     real :: dA,dB,dC"
   do face_type=1,nface_types
     write(unitno,"(a,a11,i1)") "     real :: ","face_factor",face_type
   enddo
!
endsubroutine declare_facefactors
!***********************************************************************
subroutine offset_name(vname,offset,offname)
   character(len=*), intent(inout) :: vname,offname
   integer, intent(in) :: offset
   integer :: colon
   character(len=1) :: offset_sgn = '+'
   character(len=20) :: vname1, vname2
   if (offset==0)  then
     offname=trim(vname)
     return
   endif

   if (offset>0) offset_sgn='+'
   if (offset<0) offset_sgn='-'

   colon=SCAN(vname, ':')
   if (colon==0) then
     write (offname,"(a,a1,i1)") &
        trim(vname),offset_sgn,abs(offset)
     return
   else
     vname1=vname(1:colon-1)
     vname2=vname(colon+1:len_trim(vname))
     write (offname,"(a,a1,i1,a1,a,a1,i1)") &
        trim(vname1),offset_sgn,abs(offset),":", &
        trim(vname2),offset_sgn,abs(offset)
   endif



endsubroutine offset_name
!***********************************************************************
subroutine evaluate_integral(unitno,rotation,imin,imax,jmin,jmax, &
    kmin,kmax,intname,iname,jname,kname)
!
   use SurfaceData
!
!   character(len=*), optional, intent(in) :: intname, iname, jname, kname
   integer :: rotation, unitno
   integer :: ncontrib
   integer :: i,face_type,comp_it,comp,pnt,mesh,dir1,dir2,dir3
   logical :: lfirst=.true., lfirst_term=.true.
   character :: sgn = '+'
   character(len=3), dimension(3) :: vel = (/ 'iux','iuy','iuz' /)
   character(len=*) :: intname !='integral'
   character(len=*) :: iname, jname, kname
   character(len=50) :: offname
   logical :: lskip_open_bracket = .true.
!
!   if (present(intname)) intname_=intname
!   if (present(iname)  ) iname_=iname
!   if (present(jname)  ) jname_=jname
!   if (present(kname)  ) kname_=kname
!
   dir1=mod(0+rotation,3)+1
   dir2=mod(1+rotation,3)+1
   dir3=mod(2+rotation,3)+1
   lfirst=.true.
   do face_type=1,nface_types
     ncontrib=0.
     lfirst_term=.true.
     if (lfirst) then
       write(unitno,*) "       ",trim(intname)," = ( &"
     else
       if (.not. lskip_open_bracket) write(unitno,*) "                + ( &"
     endif

     lskip_open_bracket=.false.
     do pnt=1,npoints
       if ( face_types(faces(pnt)) /= face_type ) cycle
       If (      (points(dir1,pnt)<imin).or.(points(dir1,pnt)>imax)   &
           .or. (points(dir2,pnt)<jmin).or.(points(dir2,pnt)>jmax)   &
           .or. (points(dir3,pnt)<kmin).or.(points(dir3,pnt)>kmax) ) cycle
       do comp_it=0,2
         comp=mod(comp_it+rotation,3)+1
         do mesh=1,3
           !mesh=mod(mesh_it+rotation,3)+1
         enddo
         if (normals(comp,pnt)/=0) then
           ncontrib=ncontrib+1
           if (normals(comp,pnt)>0) then
             sgn='+'
           else
             sgn='-'
           endif
           write (unitno,"(a,a1,a4)",ADVANCE='no') "                    ", sgn, " f( "
           call offset_name(trim(iname),points(dir1,pnt),offname)
           write (unitno,"(a)",ADVANCE='no') trim(offname)//"   , "
           call offset_name(trim(jname),points(dir2,pnt),offname)
           write (unitno,"(a)",ADVANCE='no') trim(offname)//"   , "
           call offset_name(trim(kname),points(dir3,pnt),offname)
           write (unitno,"(a)",ADVANCE='no') trim(offname)//"   , "

           write (unitno,"(a3,a3)") vel(comp_it+1),") &"
           lfirst_term=.false.
         endif
       enddo
     enddo
!
     if (ncontrib==0) then
       if (face_type==nface_types) then
         write (unitno,"(a)") "                    0. ) "
       else
         lskip_open_bracket=.true.
       endif
     else
       if (face_type==nface_types) then
         write (unitno,"(a,a11,i1)") &
             "                   ) * ", "face_factor",face_type
       else
         write (unitno,"(a,a11,i1,a)") &
             "                   ) * ", "face_factor",face_type," &"
       endif
     endif
     lfirst=.false.
   enddo  !Face_types
!
endsubroutine evaluate_integral
!***********************************************************************
subroutine read_surfaceinfo
!
   use SurfaceData
!
   integer :: pnt, face_type
   character(len=20) :: header
!
! Read surface information
!
   read*,header
   read*,npoints
   read*,header
   read*,ndimensions
   read*,header
   read*,nfaces
   read*,header
   read*,nface_types

   problem_dimensions=0
   if (nxgrid/=1) problem_dimensions=problem_dimensions+1
   if (nygrid/=1) problem_dimensions=problem_dimensions+1
   if (nzgrid/=1) problem_dimensions=problem_dimensions+1
!
   allocate( points(3,npoints) )
   allocate( normals(3,npoints) )
   allocate( faces(npoints) )
   allocate( face_types(nfaces) )
   allocate( area_elements(nface_types) )
!
   read*,header
   do pnt=1,npoints
     read*,points(:,pnt)
   enddo
   read*,header
   do pnt=1,npoints
     read*,normals(:,pnt)
   enddo
   read*,header
   read*,faces(:)
   read*,header
   read*,face_types(:)
   read*,header
   do face_type=1,nface_types
     read*,area_elements(face_type)
   enddo
endsubroutine read_surfaceinfo
!***********************************************************************
