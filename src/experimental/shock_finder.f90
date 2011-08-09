
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
program shock_finder2D
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
endprogram shock_finder2D
!***********************************************************************
subroutine make_calc_body(unitno)
!
  use Cparam
  use SurfaceData
!
  implicit none
!
  integer :: unitno
  write(unitno,"(a)") "!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)"
  write(unitno,"(a)") "!"
  write(unitno,"(a)") "!***********************************************************************"
  write(unitno,"(a)") "  subroutine shock_calc_body(f)"
  write(unitno,"(a)") "    use Cdata"
  write(unitno,"(a)") "    use Mpicomm, only: stop_it"
  write(unitno,"(a)") "    real, dimension (mx,my,mz,mfarray) :: f"
  call declare_facefactors(unitno)
  write(unitno,"(a)") "    integer :: i,j,k"

  if ( ndimensions /= problem_dimensions ) then
     write (unitno,"(a)") "call fatal_error('shock_calc_body',& "
     write (unitno,"(a,i1,a,i1,a)") "'The compiled in shock profile integral is ", &
                      ndimensions, &
                     "D but the problem is ", problem_dimensions,"D. STOPPING.')"
  endif
  if ((nxgrid/=1).and.(nygrid/=1).and.(nzgrid/=1)) then
!  write(unitno,"(a)") "      if (lroot) print*,'shock_calc_body: The compiled in shock profile integral is 2D'"
!  write(unitno,"(a)") "      call stop_it('shock_calc_body: This problem is 3D - INCONSITENCY')"
    call evaluate_facefactors(unitno,0)
    write(unitno,"(a)") "      do k=n1i,n2i"
    write(unitno,"(a)") "      do j=m1i,m2i"
    write(unitno,"(a)") "      do i=l1i,l2i"
    call evaluate_integral(unitno,0,-3,+3,-3,+3,-3,+3,'f(i,j,k,ishock)','i','j','k')
    write(unitno,"(a)") "      enddo ; enddo ; enddo"
  elseif ((nxgrid/=1).and.(nygrid/=1)) then
    call evaluate_facefactors(unitno,0)
    write(unitno,"(a)") "      do j=m1i,m2i"
    write(unitno,"(a)") "      do i=l1i,l2i"
    call evaluate_integral(unitno,0,-3,+3,-3,+3,-3,+3,'f(l1,j,k,ishock)','i','j','n1')
    write(unitno,"(a)") "      enddo ; enddo"
  elseif ((nxgrid/=1).and.(nzgrid/=1)) then
    call evaluate_facefactors(unitno,1)
    write(unitno,"(a)") "      do k=n1i,n2i"
    write(unitno,"(a)") "      do i=l1i,l2i"
    call evaluate_integral(unitno,1,-3,+3,-3,+3,-3,+3,'f(l1,j,k,ishock)','i','m1','k')
    write(unitno,"(a)") "      enddo ; enddo"
  elseif ((nzgrid/=1).and.(nygrid/=1)) then
    call evaluate_facefactors(unitno,2)
    write(unitno,"(a)") "      do k=n1i,n2i"
    write(unitno,"(a)") "      do j=m1i,m2i"
    call evaluate_integral(unitno,2,-3,+3,-3,+3,-3,+3,'f(l1,j,k,ishock)','l1','j','k')
    write(unitno,"(a)") "      enddo ; enddo"
  else
    write(unitno,"(a)") "      call stop_it('shock_calc_body: Case not implemented')"
  endif
  write(unitno,"(a)") "    endsubroutine shock_calc_body"
endsubroutine make_calc_body
!***********************************************************************
subroutine make_calc_3d_int_j_ji(unitno)
  use Cparam
!
  implicit none
!
  integer :: unitno,i,j,k
  character (len=4) :: istr='',jstr='',kstr=''
  character (len=5) :: pistr='',nistr=''
  character (len=5) :: pjstr='',njstr=''
  character (len=5) :: pkstr='',nkstr=''
  write(unitno,"(a)") "!***********************************************************************"
  write(unitno,"(a)") "  subroutine shock_calc_int_j_ji(f)"
  write(unitno,"(a)") "    use Cdata"
  write(unitno,"(a)") "    use Mpicomm, only: stop_it"
  write(unitno,"(a)") "    real, dimension (mx,my,mz,mfarray) :: f"
  call declare_facefactors(unitno)
  write(unitno,"(a)") "    integer :: k"

  call evaluate_facefactors(unitno,0)
! Top/bottom
  pjstr=""
  njstr=""
  do j=0,2
    jstr=itoa(j)
    if (j/=0) njstr='-'//trim(jstr)
    if (j/=0) pjstr='+'//trim(jstr)
    write(unitno,"(a)") "      do k=n1i,n2i"
    call evaluate_integral(unitno,2,-3,3,-j,+3,-3,+3, &
        'f(l1i:l2i,m1'//trim(pjstr)//',k,ishock)','l1i:l2i','m1'//trim(pjstr),'k')
    call evaluate_integral(unitno,2,-3,3,-3,+j,-3,+3, &
        'f(l1i:l2i,m2'//trim(njstr)//',k,ishock)','l1i:l2i','m2'//trim(njstr),'k')
    pistr=""
    nistr=""
    do i=0,2
      istr=itoa(i)
      if (i/=0) nistr='-'//trim(istr)
      if (i/=0) pistr='+'//trim(istr)
      call evaluate_integral(unitno,2,-i,3,-j,+3,-3,+3, &
        'f(l1'//trim(pistr)//',m1'//trim(pjstr)//',k,ishock)',&
          'l1'//trim(pistr),'m1'//trim(pjstr),'k')
      call evaluate_integral(unitno,2,-3,i,-3,+j,-3,+3, &
        'f(l2'//trim(nistr)//',m2'//trim(njstr)//',k,ishock)', &
          'l2'//trim(nistr),'m2'//trim(njstr),'k')
      call evaluate_integral(unitno,2,-i,3,-3,+j,-3,+3, &
        'f(l1'//trim(pistr)//',m2'//trim(njstr)//',k,ishock)',&
          'l1'//trim(pistr),'m2'//trim(njstr),'k')
      call evaluate_integral(unitno,2,-3,i,-j,+3,-3,+3, &
        'f(l2'//trim(nistr)//',m1'//trim(pjstr)//',k,ishock)', &
          'l2'//trim(nistr),'m1'//trim(pjstr),'k')
    enddo
    write(unitno,"(a)") "      enddo"
  enddo
  write(unitno,"(a)") "    endsubroutine shock_calc_int_j_ji"
endsubroutine make_calc_3d_int_j_ji
!***********************************************************************
subroutine make_calc_3d_int_k_ki(unitno)
  use Cparam
!
  implicit none
!
  integer :: unitno,i,j,k
  character (len=4) :: istr='',jstr='',kstr=''
  character (len=5) :: pistr='',nistr=''
  character (len=5) :: pjstr='',njstr=''
  character (len=5) :: pkstr='',nkstr=''
  write(unitno,"(a)") "!***********************************************************************"
  write(unitno,"(a)") "  subroutine shock_calc_int_k_ki(f)"
  write(unitno,"(a)") "    use Cdata"
  write(unitno,"(a)") "    use Mpicomm, only: stop_it"
  write(unitno,"(a)") "    real, dimension (mx,my,mz,mfarray) :: f"
  call declare_facefactors(unitno)
  write(unitno,"(a)") "    integer :: j"

  call evaluate_facefactors(unitno,0)
! Left/right
  pkstr=""
  nkstr=""
  do k=0,2
    kstr=itoa(k)
    if (k/=0) nkstr='-'//trim(kstr)
    if (k/=0) pkstr='+'//trim(kstr)
    write(unitno,"(a)") "      do j=m1i,m2i"
    call evaluate_integral(unitno,2,-3,3,-3,+3,-k,+3, &
        'f(l1i:l2i,j,n1'//trim(pkstr)//',ishock)','l1i:l2i','j','n1'//trim(pkstr))
    call evaluate_integral(unitno,2,-3,3,-3,+3,-3,+k, &
        'f(l1i:l2i,j,n2'//trim(nkstr)//',ishock)','l1i:l2i','j','n2'//trim(nkstr))
    pistr=""
    nistr=""
    do i=0,2
      istr=itoa(i)
      if (i/=0) nistr='-'//trim(istr)
      if (i/=0) pistr='+'//trim(istr)
      call evaluate_integral(unitno,2,-i,3,-3,+3,-k,+3, &
        'f(l1'//trim(pistr)//',j,n1'//trim(pkstr)//',ishock)',&
          'l1'//trim(pistr),'j','n1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,i,-3,+3,-3,+k, &
        'f(l2'//trim(nistr)//',j,n2'//trim(nkstr)//',ishock)',&
          'l2'//trim(nistr),'j','n2'//trim(nkstr))
      call evaluate_integral(unitno,2,-i,3,-3,+3,-3,+k, &
        'f(l1'//trim(pistr)//',j,n2'//trim(nkstr)//',ishock)',&
          'l1'//trim(pistr),'j','n2'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,i,-3,+3,-k,+3, &
        'f(l2'//trim(nistr)//',j,n1'//trim(pkstr)//',ishock)',&
          'l2'//trim(nistr),'j','n1'//trim(pkstr))
    enddo
    write(unitno,"(a)") "      enddo"
  enddo

  write(unitno,"(a)") "    endsubroutine shock_calc_int_k_ki"
endsubroutine make_calc_3d_int_k_ki
!***********************************************************************
subroutine make_calc_3d_int_kj_kji(unitno,k)
  use Cparam
!
  implicit none
!
  integer :: unitno,i,j,k
  character (len=4) :: istr='',jstr='',kstr=''
  character (len=5) :: pistr='',nistr=''
  character (len=5) :: pjstr='',njstr=''
  character (len=5) :: pkstr='',nkstr=''
  write(unitno,"(a)") "!***********************************************************************"
  write(unitno,"(a,i1,a)") "  subroutine shock_calc_int_kj_kji_",k,"(f)"
  write(unitno,"(a)") "    use Cdata"
  write(unitno,"(a)") "    use Mpicomm, only: stop_it"
  write(unitno,"(a)") "    real, dimension (mx,my,mz,mfarray) :: f"
  call declare_facefactors(unitno)

  call evaluate_facefactors(unitno,0)
! Corners
  pkstr=""
  nkstr=""
  pjstr=""
  njstr=""
!  do k=0,2
  do j=0,2
  kstr=itoa(k)
  jstr=itoa(j)
  if (k/=0) nkstr='-'//trim(kstr)
  if (k/=0) pkstr='+'//trim(kstr)
  if (j/=0) njstr='-'//trim(jstr)
  if (j/=0) pjstr='+'//trim(jstr)
  call evaluate_integral(unitno,2,-3,3,-j,+3,-k,+3, &
      'f(l1i:l2i,m1'//trim(pjstr)//',n1'//trim(pkstr)//',ishock)', &
      'l1i:l2i','m1'//trim(pjstr),'n1'//trim(pkstr))
  call evaluate_integral(unitno,2,-3,3,-3,+j,-k,+3, &
      'f(l1i:l2i,m2'//trim(njstr)//',n1'//trim(pkstr)//',ishock)', &
      'l1i:l2i','m2'//trim(njstr),'n1'//trim(pkstr))
  call evaluate_integral(unitno,2,-3,3,-3,+j,-3,+k, &
      'f(l1i:l2i,m2'//trim(njstr)//',n2'//trim(nkstr)//',ishock)', &
      'l1i:l2i','m2'//trim(njstr),'n2'//trim(nkstr))
  call evaluate_integral(unitno,2,-3,3,-j,+3,-3,+k, &
      'f(l1i:l2i,m1'//trim(pjstr)//',n2'//trim(nkstr)//',ishock)', &
      'l1i:l2i','m1'//trim(pjstr),'n2'//trim(nkstr))
    pistr=""
    nistr=""
    do i=0,2
      istr=itoa(i)
      if (i/=0) nistr='-'//trim(istr)
      if (i/=0) pistr='+'//trim(istr)
      call evaluate_integral(unitno,2,-i,3,-j,+3,-k,+3, &
          'f(l1'//trim(pistr)//',m1'//trim(pjstr)//',n1'//trim(pkstr)//',ishock)', &
          'l1'//trim(pistr),'m1'//trim(pjstr),'n1'//trim(pkstr))
      call evaluate_integral(unitno,2,-i,3,-3,+j,-k,+3, &
          'f(l1'//trim(pistr)//',m2'//trim(njstr)//',n1'//trim(pkstr)//',ishock)', &
          'l1'//trim(pistr),'m2'//trim(njstr),'n1'//trim(pkstr))
      call evaluate_integral(unitno,2,-i,3,-3,+j,-3,+k, &
          'f(l1'//trim(pistr)//',m2'//trim(njstr)//',n2'//trim(nkstr)//',ishock)', &
          'l1'//trim(pistr),'m2'//trim(njstr),'n2'//trim(nkstr))
      call evaluate_integral(unitno,2,-i,3,-j,+3,-3,+k, &
          'f(l1'//trim(pistr)//',m1'//trim(pjstr)//',n2'//trim(nkstr)//',ishock)', &
          'l1'//trim(pistr),'m1'//trim(pjstr),'n2'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,i,-j,+3,-k,+3, &
          'f(l2'//trim(nistr)//',m1'//trim(pjstr)//',n1'//trim(pkstr)//',ishock)', &
          'l2'//trim(nistr),'m1'//trim(pjstr),'n1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,i,-3,+j,-k,+3, &
          'f(l2'//trim(nistr)//',m2'//trim(njstr)//',n1'//trim(pkstr)//',ishock)', &
          'l2'//trim(nistr),'m2'//trim(njstr),'n1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,i,-3,+j,-3,+k, &
          'f(l2'//trim(nistr)//',m2'//trim(njstr)//',n2'//trim(nkstr)//',ishock)', &
          'l2'//trim(nistr),'m2'//trim(njstr),'n2'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,i,-j,+3,-3,+k, &
          'f(l2'//trim(nistr)//',m1'//trim(pjstr)//',n2'//trim(nkstr)//',ishock)', &
          'l2'//trim(nistr),'m1'//trim(pjstr),'n2'//trim(nkstr))
    enddo
  enddo
!  enddo

  write(unitno,"(a,i1)") "    endsubroutine shock_calc_int_kj_kji_",k
endsubroutine make_calc_3d_int_kj_kji
!***********************************************************************
subroutine make_calc_internalboundary(unitno)
  use Cparam
!
  implicit none
!
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
  call declare_facefactors(unitno)
  write(unitno,"(a)") "    integer :: i,j,k"
  if ((nxgrid/=1).and.(nygrid/=1).and.(nzgrid/=1)) then
!    write(unitno,"(a)") "    real, dimension(l2i-l1i+1) :: integral1,integral2"
    write(unitno,"(a)") "!  make_calc_internalboundary: 3D CASE"
    call evaluate_facefactors(unitno,0)

    write(unitno,"(a)") " call shock_calc_int_j_ji(f)"
    write(unitno,"(a)") " call shock_calc_int_k_ki(f)"
    write(unitno,"(a)") " call shock_calc_int_kj_kji_0(f)"
    write(unitno,"(a)") " call shock_calc_int_kj_kji_1(f)"
    write(unitno,"(a)") " call shock_calc_int_kj_kji_2(f)"

!  write(unitno,"(a)") "      if (lroot) print*,'shock_calc_body: The compiled in shock profile integral is 2D'"
!  write(unitno,"(a)") "      call stop_it('shock_calc_body: This problem is 3D - INCONSITENCY')"
  elseif ((nxgrid/=1).and.(nygrid/=1)) then
  write(unitno,"(a)") "    real :: integral1,integral2"
  call evaluate_facefactors(unitno,0)
! Top/bottom
  do j=0,2
  jstr=itoa(j)
  if (j/=0) njstr='-'//trim(jstr)
  if (j/=0) pjstr='+'//trim(jstr)
  write(unitno,"(a)") "      do i=l1i+1,l2i-1"
  call evaluate_integral(unitno,0,-3,+3,-j,+3,0,0, &
      'integral1','i','m1'//trim(pjstr),'n1')
  call evaluate_integral(unitno,0,-3,+3,-3,+j,0,0, &
      'integral2','i','m2'//trim(njstr),'n1')
  write(unitno,"(a)") "        f(i,m1"//trim(pjstr)//",n1,ishock)=integral1"
  write(unitno,"(a)") "        f(i,m2"//trim(njstr)//",n1,ishock)=integral2"
  write(unitno,"(a)") "      enddo"
  enddo
! Left/right
  do i=0,2
  istr=itoa(i)
  if (i/=0) nistr='-'//trim(istr)
  if (i/=0) pistr='+'//trim(istr)
  write(unitno,"(a)") "      do j=m1i+1,m2i-1"
  call evaluate_integral(unitno,0,-i,+3,-3,+3,0,0, &
      'integral1','l1'//trim(pistr),'j','n1')
  call evaluate_integral(unitno,0,-3,+i,-3,+3,0,0, &
      'integral2','l2'//trim(nistr),'j','n1')
  write(unitno,"(a)") "        f(l1"//trim(pistr)//",j,n1,ishock)=integral1"
  write(unitno,"(a)") "        f(l2"//trim(nistr)//",j,n1,ishock)=integral2"
  write(unitno,"(a)") "      enddo"
  enddo
! Corners
  do j=0,2
  do i=0,2
  istr=itoa(i)
  jstr=itoa(j)
  if (i/=0) nistr='-'//trim(istr)
  if (i/=0) pistr='+'//trim(istr)
  if (j/=0) njstr='-'//trim(jstr)
  if (j/=0) pjstr='+'//trim(jstr)
  call evaluate_integral(unitno,0,-i,+3,-j,+3,0,0, &
      'f(l1'//trim(pistr)//',m1'//trim(pjstr)//',n1,ishock)', &
      'l1'//trim(pistr),'m1'//trim(pjstr),'n1')
  call evaluate_integral(unitno,0,-3,+i,-j,+3,0,0, &
      'f(l2'//trim(nistr)//',m1'//trim(pjstr)//',n1,ishock)', &
      'l2'//trim(nistr),'m1'//trim(pjstr),'n1')
  call evaluate_integral(unitno,0,-3,+i,-3,+j,0,0, &
      'f(l2'//trim(nistr)//',m2'//trim(njstr)//',n1,ishock)', &
      'l2'//trim(nistr),'m2'//trim(njstr),'n1')
  call evaluate_integral(unitno,0,-i,+3,-3,+j,0,0, &
      'f(l1'//trim(pistr)//',m2'//trim(njstr)//',n1,ishock)', &
      'l1'//trim(pistr),'m2'//trim(njstr),'n1')
  enddo
  enddo
!
  elseif ((nxgrid/=1).and.(nzgrid/=1)) then
  write(unitno,"(a)") "    real :: integral1,integral2"
  call evaluate_facefactors(unitno,1)
! Top/bottom
  do k=0,2
  kstr=itoa(k)
  if (k/=0) nkstr='-'//trim(kstr)
  if (k/=0) pkstr='+'//trim(kstr)
  write(unitno,"(a)") "      do i=l1i+1,l2i-1"
  call evaluate_integral(unitno,1,-3,+3,0,0,-k,+3, &
      'integral1','i','m1','n1'//trim(pkstr))
  call evaluate_integral(unitno,1,-3,+3,0,0,-3,+k, &
      'integral2','i','m1','n2'//trim(nkstr))
  write(unitno,"(a)") "        f(i,m1,n1"//trim(pkstr)//",ishock)=integral1"
  write(unitno,"(a)") "        f(i,m1,n2"//trim(nkstr)//",ishock)=integral2"
  write(unitno,"(a)") "      enddo"
  enddo
! Left/right
  do i=0,2
  istr=itoa(i)
  if (i/=0) nistr='-'//trim(istr)
  if (i/=0) pistr='+'//trim(istr)
  write(unitno,"(a)") "      do k=n1i+1,n2i-1"
  call evaluate_integral(unitno,1,-i,+3,0,0,-3,+3, &
      'integral1','l1'//trim(pistr),'m1','k')
  call evaluate_integral(unitno,1,-3,+i,0,0,-3,+3, &
      'integral2','l2'//trim(nistr),'m1','k')
  write(unitno,"(a)") "        f(l1"//trim(pistr)//",m1,k,ishock)=integral1"
  write(unitno,"(a)") "        f(l2"//trim(nistr)//",m1,k,ishock)=integral2"
  write(unitno,"(a)") "      enddo"
  enddo
! Corners
  do k=0,2
  do i=0,2
  istr=itoa(i)
  kstr=itoa(k)
  if (i/=0) nistr='-'//trim(istr)
  if (i/=0) pistr='+'//trim(istr)
  if (k/=0) nkstr='-'//trim(kstr)
  if (k/=0) pkstr='+'//trim(kstr)
  call evaluate_integral(unitno,1,-i,+3,0,0,-k,+3, &
      'f(l1'//trim(pistr)//',m1,n1'//trim(pkstr)//',ishock)', &
      'l1'//trim(pistr),'m1','n1'//trim(pkstr))
  call evaluate_integral(unitno,1,-3,+i,0,0,-k,+3, &
      'f(l2'//trim(nistr)//',m1,n1'//trim(pkstr)//',ishock)', &
      'l2'//trim(nistr),'m1','n1'//trim(pkstr))
  call evaluate_integral(unitno,1,-3,+i,0,0,-3,+k, &
      'f(l2'//trim(nistr)//',m1,n2'//trim(nkstr)//',ishock)', &
      'l2'//trim(nistr),'m1','n2'//trim(nkstr))
  call evaluate_integral(unitno,1,-i,+3,0,0,-3,+k, &
      'f(l1'//trim(pistr)//',m1,n2'//trim(nkstr)//',ishock)', &
      'l1'//trim(pistr),'m1','n2'//trim(nkstr))
  enddo
  enddo
!
  elseif ((nzgrid/=1).and.(nygrid/=1)) then
  write(unitno,"(a)") "    real :: integral1,integral2"
  call evaluate_facefactors(unitno,2)
! Top/bottom
  do j=0,2
  jstr=itoa(j)
  if (j/=0) njstr='-'//trim(jstr)
  if (j/=0) pjstr='+'//trim(jstr)
  write(unitno,"(a)") "      do k=n1i+1,n2i-1"
  call evaluate_integral(unitno,2,0,0,-j,+3,-3,+3, &
      'integral1','l1','m1'//trim(pjstr),'k')
  call evaluate_integral(unitno,2,0,0,-3,+j,-3,+3, &
      'integral2','l1','m2'//trim(njstr),'k')
  write(unitno,"(a)") "        f(l1,m1"//trim(pjstr)//",k,ishock)=integral1"
  write(unitno,"(a)") "        f(l1,m2"//trim(njstr)//",k,ishock)=integral2"
  write(unitno,"(a)") "      enddo"
  enddo
! Left/right
  do k=0,2
  kstr=itoa(k)
  if (k/=0) nkstr='-'//trim(kstr)
  if (k/=0) pkstr='+'//trim(kstr)
  write(unitno,"(a)") "      do j=m1i+1,m2i-1"
  call evaluate_integral(unitno,2,0,0,-3,+3,-k,+3, &
      'integral1','l1','j','n1'//trim(pkstr))
  call evaluate_integral(unitno,2,0,0,-3,+3,-3,+k, &
      'integral2','l1','j','n2'//trim(nkstr))
  write(unitno,"(a)") "        f(l1,j,n1"//trim(pkstr)//",ishock)=integral1"
  write(unitno,"(a)") "        f(l1,j,n2"//trim(nkstr)//",ishock)=integral2"
  write(unitno,"(a)") "      enddo"
  enddo
! Corners
  do k=0,2
  do j=0,2
  kstr=itoa(k)
  jstr=itoa(j)
  if (k/=0) nkstr='-'//trim(kstr)
  if (k/=0) pkstr='+'//trim(kstr)
  if (j/=0) njstr='-'//trim(jstr)
  if (j/=0) pjstr='+'//trim(jstr)
  call evaluate_integral(unitno,2,0,0,-j,+3,-k,+3, &
      'f(l1,m1'//trim(pjstr)//',n1'//trim(pkstr)//',ishock)', &
      'l1','m1'//trim(pjstr),'n1'//trim(pkstr))
  call evaluate_integral(unitno,2,0,0,-3,+j,-k,+3, &
      'f(l1,m2'//trim(njstr)//',n1'//trim(pkstr)//',ishock)', &
      'l1','m2'//trim(njstr),'n1'//trim(pkstr))
  call evaluate_integral(unitno,2,0,0,-3,+j,-3,+k, &
      'f(l1,m2'//trim(njstr)//',n2'//trim(nkstr)//',ishock)', &
      'l1','m2'//trim(njstr),'n2'//trim(nkstr))
  call evaluate_integral(unitno,2,0,0,-j,+3,-3,+k, &
      'f(l1,m1'//trim(pjstr)//',n2'//trim(nkstr)//',ishock)', &
      'l1','m1'//trim(pjstr),'n2'//trim(nkstr))
  enddo
  enddo
!
  else
  write(unitno,"(a)") "      call stop_it('shock_calc_body: Case not implemented')"
  endif
  write(unitno,"(a)") "  endsubroutine shock_calc_internalboundary"

  if ((nxgrid/=1).and.(nygrid/=1).and.(nzgrid/=1)) then
    call make_calc_3d_int_k_ki(unitno)
    call make_calc_3d_int_j_ji(unitno)
    call make_calc_3d_int_kj_kji(unitno,0)
    call make_calc_3d_int_kj_kji(unitno,1)
    call make_calc_3d_int_kj_kji(unitno,2)
  endif
endsubroutine make_calc_internalboundary
!***********************************************************************
subroutine make_calc_3d_ext_j_ji(unitno)
  use Cparam
!
  implicit none
!
  integer :: unitno,i,j,k
  character (len=4) :: istr='',jstr='',kstr=''
  character (len=5) :: pistr='',nistr=''
  character (len=5) :: pjstr='',njstr=''
  character (len=5) :: pkstr='',nkstr=''
  write(unitno,"(a)") "!***********************************************************************"
  write(unitno,"(a)") "  subroutine shock_calc_ext_j_ji(f)"
  write(unitno,"(a)") "    use Cdata"
  write(unitno,"(a)") "    use Mpicomm, only: stop_it"
  write(unitno,"(a)") "    real, dimension (mx,my,mz,mfarray) :: f"
  call declare_facefactors(unitno)
  write(unitno,"(a)") "    integer :: k"

  call evaluate_facefactors(unitno,0)
! Top/bottom
  pjstr=""
  njstr=""
  do j=0,2
    jstr=itoa(j)
    if (j/=0) njstr='-'//trim(jstr)
    if (j/=0) pjstr='+'//trim(jstr)
    write(unitno,"(a)") "      do k=n1i,n2i"
    call evaluate_integral(unitno,2,-3,3,3-j,+3,-3,+3, &
        'f(l1i:l2i,1'//trim(pjstr)//',k,ishock)','l1i:l2i','1'//trim(pjstr),'k')
    call evaluate_integral(unitno,2,-3,3,-3,-3+j,-3,+3, &
        'f(l1i:l2i,my'//trim(njstr)//',k,ishock)','l1i:l2i','my'//trim(njstr),'k')
    pistr=""
    nistr=""
    do i=0,2
      istr=itoa(i)
      if (i/=0) nistr='-'//trim(istr)
      if (i/=0) pistr='+'//trim(istr)
      call evaluate_integral(unitno,2,-i,3,3-j,+3,-3,+3, &
        'f(l1'//trim(pistr)//',1'//trim(pjstr)//',k,ishock)',&
          'l1'//trim(pistr),'1'//trim(pjstr),'k')
      call evaluate_integral(unitno,2,-3,+i,-3,-3+j,-3,+3, &
        'f(l2'//trim(nistr)//',my'//trim(njstr)//',k,ishock)', &
          'l2'//trim(nistr),'my'//trim(njstr),'k')
      call evaluate_integral(unitno,2,-i,3,-3,-3+j,-3,+3, &
        'f(l1'//trim(pistr)//',my'//trim(njstr)//',k,ishock)',&
          'l1'//trim(pistr),'my'//trim(njstr),'k')
      call evaluate_integral(unitno,2,-3,+i,3-j,+3,-3,+3, &
        'f(l2'//trim(nistr)//',1'//trim(pjstr)//',k,ishock)', &
          'l2'//trim(nistr),'1'//trim(pjstr),'k')


      call evaluate_integral(unitno,2,3-i,3,3-j,+3,-3,+3, &
        'f(1'//trim(pistr)//',1'//trim(pjstr)//',k,ishock)',&
          '1'//trim(pistr),'1'//trim(pjstr),'k')
      call evaluate_integral(unitno,2,-3,-3+i,-3,-3+j,-3,+3, &
        'f(mx'//trim(nistr)//',my'//trim(njstr)//',k,ishock)', &
          'mx'//trim(nistr),'my'//trim(njstr),'k')
      call evaluate_integral(unitno,2,3-i,3,-3,-3+j,-3,+3, &
        'f(1'//trim(pistr)//',my'//trim(njstr)//',k,ishock)',&
          '1'//trim(pistr),'my'//trim(njstr),'k')
      call evaluate_integral(unitno,2,-3,-3+i,3-j,+3,-3,+3, &
        'f(mx'//trim(nistr)//',1'//trim(pjstr)//',k,ishock)', &
          'mx'//trim(nistr),'1'//trim(pjstr),'k')
    enddo
    write(unitno,"(a)") "      enddo"
  enddo
  write(unitno,"(a)") "    endsubroutine shock_calc_ext_j_ji"
endsubroutine make_calc_3d_ext_j_ji
!***********************************************************************
subroutine make_calc_3d_ext_k_ki(unitno)
  use Cparam
!
  implicit none
!
  integer :: unitno,i,j,k
  character (len=4) :: istr='',jstr='',kstr=''
  character (len=5) :: pistr='',nistr=''
  character (len=5) :: pjstr='',njstr=''
  character (len=5) :: pkstr='',nkstr=''
  write(unitno,"(a)") "!***********************************************************************"
  write(unitno,"(a)") "  subroutine shock_calc_ext_k_ki(f)"
  write(unitno,"(a)") "    use Cdata"
  write(unitno,"(a)") "    use Mpicomm, only: stop_it"
  write(unitno,"(a)") "    real, dimension (mx,my,mz,mfarray) :: f"
  call declare_facefactors(unitno)
  write(unitno,"(a)") "    integer :: j"

  call evaluate_facefactors(unitno,0)
! Left/right
  pkstr=""
  nkstr=""
  do k=0,2
    kstr=itoa(k)
    if (k/=0) nkstr='-'//trim(kstr)
    if (k/=0) pkstr='+'//trim(kstr)
    write(unitno,"(a)") "      do j=m1i,m2i"
    call evaluate_integral(unitno,2,-3,3,-3,+3,3-k,+3, &
        'f(l1i:l2i,j,1'//trim(pkstr)//',ishock)','l1i:l2i','j','1'//trim(pkstr))
    call evaluate_integral(unitno,2,-3,3,-3,+3,-3,-3+k, &
        'f(l1i:l2i,j,mz'//trim(nkstr)//',ishock)','l1i:l2i','j','mz'//trim(nkstr))
    pistr=""
    nistr=""
    do i=0,2
      istr=itoa(i)
      if (i/=0) nistr='-'//trim(istr)
      if (i/=0) pistr='+'//trim(istr)
      call evaluate_integral(unitno,2,-i,3,-3,+3,3-k,+3, &
        'f(l1'//trim(pistr)//',j,1'//trim(pkstr)//',ishock)',&
          'l1'//trim(pistr),'j','1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,+i,-3,+3,-3,-3+k, &
        'f(l2'//trim(nistr)//',j,mz'//trim(nkstr)//',ishock)',&
          'l2'//trim(nistr),'j','mz'//trim(nkstr))
      call evaluate_integral(unitno,2,-i,3,-3,+3,-3,-3+k, &
        'f(l1'//trim(pistr)//',j,mz'//trim(nkstr)//',ishock)',&
          'l1'//trim(pistr),'j','mz'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,+i,-3,+3,3-k,+3, &
        'f(l2'//trim(nistr)//',j,1'//trim(pkstr)//',ishock)',&
          'l2'//trim(nistr),'j','1'//trim(pkstr))

      call evaluate_integral(unitno,2,3-i,3,-3,+3,3-k,+3, &
        'f(1'//trim(pistr)//',j,1'//trim(pkstr)//',ishock)',&
          '1'//trim(pistr),'j','1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,-3+i,-3,+3,-3,-3+k, &
        'f(mx'//trim(nistr)//',j,mz'//trim(nkstr)//',ishock)',&
          'mx'//trim(nistr),'j','mz'//trim(nkstr))
      call evaluate_integral(unitno,2,3-i,3,-3,+3,-3,-3+k, &
        'f(1'//trim(pistr)//',j,mz'//trim(nkstr)//',ishock)',&
          '1'//trim(pistr),'j','mz'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,-3+i,-3,+3,3-k,+3, &
        'f(mx'//trim(nistr)//',j,1'//trim(pkstr)//',ishock)',&
          'mx'//trim(nistr),'j','1'//trim(pkstr))
    enddo
    write(unitno,"(a)") "      enddo"
  enddo

  write(unitno,"(a)") "    endsubroutine shock_calc_ext_k_ki"
endsubroutine make_calc_3d_ext_k_ki
!***********************************************************************
subroutine make_calc_3d_ext_kj_kji(unitno,k)
  use Cparam
!
  implicit none
!
  integer :: unitno,i,j,k
  character (len=4) :: istr='',jstr='',kstr=''
  character (len=5) :: pistr='',nistr=''
  character (len=5) :: pjstr='',njstr=''
  character (len=5) :: pkstr='',nkstr=''
  write(unitno,"(a)") "!***********************************************************************"
  write(unitno,"(a,i1,a)") "  subroutine shock_calc_ext_kj_kji_",k,"(f)"
  write(unitno,"(a)") "    use Cdata"
  write(unitno,"(a)") "    use Mpicomm, only: stop_it"
  write(unitno,"(a)") "    real, dimension (mx,my,mz,mfarray) :: f"
  call declare_facefactors(unitno)

  call evaluate_facefactors(unitno,0)
! Corners
  pkstr=""
  nkstr=""
  pjstr=""
  njstr=""
!  do k=0,2
  do j=0,2
  kstr=itoa(k)
  jstr=itoa(j)
  if (k/=0) nkstr='-'//trim(kstr)
  if (k/=0) pkstr='+'//trim(kstr)
  if (j/=0) njstr='-'//trim(jstr)
  if (j/=0) pjstr='+'//trim(jstr)
  call evaluate_integral(unitno,2,-3,3,-j,+3,3-k,+3, &
      'f(l1i:l2i,m1'//trim(pjstr)//',1'//trim(pkstr)//',ishock)', &
      'l1i:l2i','m1'//trim(pjstr),'1'//trim(pkstr))
  call evaluate_integral(unitno,2,-3,3,-3,+j,3-k,+3, &
      'f(l1i:l2i,m2'//trim(njstr)//',1'//trim(pkstr)//',ishock)', &
      'l1i:l2i','m2'//trim(njstr),'1'//trim(pkstr))
  call evaluate_integral(unitno,2,-3,3,-3,+j,-3,-3+k, &
      'f(l1i:l2i,m2'//trim(njstr)//',mz'//trim(nkstr)//',ishock)', &
      'l1i:l2i','m2'//trim(njstr),'mz'//trim(nkstr))
  call evaluate_integral(unitno,2,-3,3,-j,+3,-3,-3+k, &
      'f(l1i:l2i,m1'//trim(pjstr)//',mz'//trim(nkstr)//',ishock)', &
      'l1i:l2i','m1'//trim(pjstr),'mz'//trim(nkstr))

  call evaluate_integral(unitno,2,-3,3,3-j,+3,-k,+3, &
      'f(l1i:l2i,1'//trim(pjstr)//',n1'//trim(pkstr)//',ishock)', &
      'l1i:l2i','1'//trim(pjstr),'n1'//trim(pkstr))
  call evaluate_integral(unitno,2,-3,3,-3,-3+j,-k,+3, &
      'f(l1i:l2i,my'//trim(njstr)//',n1'//trim(pkstr)//',ishock)', &
      'l1i:l2i','my'//trim(njstr),'n1'//trim(pkstr))
  call evaluate_integral(unitno,2,-3,3,-3,-3+j,-3,+k, &
      'f(l1i:l2i,my'//trim(njstr)//',n2'//trim(nkstr)//',ishock)', &
      'l1i:l2i','my'//trim(njstr),'n2'//trim(nkstr))
  call evaluate_integral(unitno,2,-3,3,3-j,+3,-3,+k, &
      'f(l1i:l2i,1'//trim(pjstr)//',n2'//trim(nkstr)//',ishock)', &
      'l1i:l2i','1'//trim(pjstr),'n2'//trim(nkstr))

  call evaluate_integral(unitno,2,-3,3,3-j,+3,3-k,+3, &
      'f(l1i:l2i,1'//trim(pjstr)//',1'//trim(pkstr)//',ishock)', &
      'l1i:l2i','1'//trim(pjstr),'1'//trim(pkstr))
  call evaluate_integral(unitno,2,-3,3,-3,-3+j,3-k,+3, &
      'f(l1i:l2i,my'//trim(njstr)//',1'//trim(pkstr)//',ishock)', &
      'l1i:l2i','my'//trim(njstr),'1'//trim(pkstr))
  call evaluate_integral(unitno,2,-3,3,-3,-3+j,-3,-3+k, &
      'f(l1i:l2i,my'//trim(njstr)//',mz'//trim(nkstr)//',ishock)', &
      'l1i:l2i','my'//trim(njstr),'mz'//trim(nkstr))
  call evaluate_integral(unitno,2,-3,3,3-j,+3,-3,-3+k, &
      'f(l1i:l2i,1'//trim(pjstr)//',mz'//trim(nkstr)//',ishock)', &
      'l1i:l2i','1'//trim(pjstr),'mz'//trim(nkstr))
    pistr=""
    nistr=""
    do i=0,2
      istr=itoa(i)
      if (i/=0) nistr='-'//trim(istr)
      if (i/=0) pistr='+'//trim(istr)
      call evaluate_integral(unitno,2,-i,3,-j,+3,3-k,+3, &
          'f(l1'//trim(pistr)//',m1'//trim(pjstr)//',1'//trim(pkstr)//',ishock)', &
          'l1'//trim(pistr),'m1'//trim(pjstr),'1'//trim(pkstr))
      call evaluate_integral(unitno,2,-i,3,-3,+j,3-k,+3, &
          'f(l1'//trim(pistr)//',m2'//trim(njstr)//',1'//trim(pkstr)//',ishock)', &
          'l1'//trim(pistr),'m2'//trim(njstr),'1'//trim(pkstr))
      call evaluate_integral(unitno,2,-i,3,-3,+j,-3,-3+k, &
          'f(l1'//trim(pistr)//',m2'//trim(njstr)//',mz'//trim(nkstr)//',ishock)', &
          'l1'//trim(pistr),'m2'//trim(njstr),'mz'//trim(nkstr))
      call evaluate_integral(unitno,2,-i,3,-j,+3,-3,-3+k, &
          'f(l1'//trim(pistr)//',m1'//trim(pjstr)//',mz'//trim(nkstr)//',ishock)', &
          'l1'//trim(pistr),'m1'//trim(pjstr),'mz'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,+i,-j,+3,3-k,+3, &
          'f(l2'//trim(nistr)//',m1'//trim(pjstr)//',1'//trim(pkstr)//',ishock)', &
          'l2'//trim(nistr),'m1'//trim(pjstr),'1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,+i,-3,+j,3-k,+3, &
          'f(l2'//trim(nistr)//',m2'//trim(njstr)//',1'//trim(pkstr)//',ishock)', &
          'l2'//trim(nistr),'m2'//trim(njstr),'1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,+i,-3,+j,-3,-3+k, &
          'f(l2'//trim(nistr)//',m2'//trim(njstr)//',mz'//trim(nkstr)//',ishock)', &
          'l2'//trim(nistr),'m2'//trim(njstr),'mz'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,+i,-j,+3,-3,-3+k, &
          'f(l2'//trim(nistr)//',m1'//trim(pjstr)//',mz'//trim(nkstr)//',ishock)', &
          'l2'//trim(nistr),'m1'//trim(pjstr),'mz'//trim(nkstr))

      call evaluate_integral(unitno,2,-i,3,3-j,+3,-k,+3, &
          'f(l1'//trim(pistr)//',1'//trim(pjstr)//',n1'//trim(pkstr)//',ishock)', &
          'l1'//trim(pistr),'1'//trim(pjstr),'n1'//trim(pkstr))
      call evaluate_integral(unitno,2,-i,3,-3,-3+j,-k,+3, &
          'f(l1'//trim(pistr)//',my'//trim(njstr)//',n1'//trim(pkstr)//',ishock)', &
          'l1'//trim(pistr),'my'//trim(njstr),'n1'//trim(pkstr))
      call evaluate_integral(unitno,2,-i,3,-3,-3+j,-3,+k, &
          'f(l1'//trim(pistr)//',my'//trim(njstr)//',n2'//trim(nkstr)//',ishock)', &
          'l1'//trim(pistr),'my'//trim(njstr),'n2'//trim(nkstr))
      call evaluate_integral(unitno,2,-i,3,3-j,+3,-3,+k, &
          'f(l1'//trim(pistr)//',1'//trim(pjstr)//',n2'//trim(nkstr)//',ishock)', &
          'l1'//trim(pistr),'1'//trim(pjstr),'n2'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,+i,3-j,+3,-k,+3, &
          'f(l2'//trim(nistr)//',1'//trim(pjstr)//',n1'//trim(pkstr)//',ishock)', &
          'l2'//trim(nistr),'1'//trim(pjstr),'n1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,+i,-3,-3+j,-k,+3, &
          'f(l2'//trim(nistr)//',my'//trim(njstr)//',n1'//trim(pkstr)//',ishock)', &
          'l2'//trim(nistr),'my'//trim(njstr),'n1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,+i,-3,-3+j,-3,+k, &
          'f(l2'//trim(nistr)//',my'//trim(njstr)//',n2'//trim(nkstr)//',ishock)', &
          'l2'//trim(nistr),'my'//trim(njstr),'n2'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,+i,3-j,+3,-3,+k, &
          'f(l2'//trim(nistr)//',1'//trim(pjstr)//',n2'//trim(nkstr)//',ishock)', &
          'l2'//trim(nistr),'1'//trim(pjstr),'n2'//trim(nkstr))

      call evaluate_integral(unitno,2,-i,3,3-j,+3,3-k,+3, &
          'f(l1'//trim(pistr)//',1'//trim(pjstr)//',1'//trim(pkstr)//',ishock)', &
          'l1'//trim(pistr),'1'//trim(pjstr),'1'//trim(pkstr))
      call evaluate_integral(unitno,2,-i,3,-3,-3+j,3-k,+3, &
          'f(l1'//trim(pistr)//',my'//trim(njstr)//',1'//trim(pkstr)//',ishock)', &
          'l1'//trim(pistr),'my'//trim(njstr),'1'//trim(pkstr))
      call evaluate_integral(unitno,2,-i,3,-3,-3+j,-3,-3+k, &
          'f(l1'//trim(pistr)//',my'//trim(njstr)//',mz'//trim(nkstr)//',ishock)', &
          'l1'//trim(pistr),'my'//trim(njstr),'mz'//trim(nkstr))
      call evaluate_integral(unitno,2,-i,3,3-j,+3,-3,-3+k, &
          'f(l1'//trim(pistr)//',1'//trim(pjstr)//',mz'//trim(nkstr)//',ishock)', &
          'l1'//trim(pistr),'1'//trim(pjstr),'mz'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,+i,3-j,+3,3-k,+3, &
          'f(l2'//trim(nistr)//',1'//trim(pjstr)//',1'//trim(pkstr)//',ishock)', &
          'l2'//trim(nistr),'1'//trim(pjstr),'1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,+i,-3,-3+j,3-k,+3, &
          'f(l2'//trim(nistr)//',my'//trim(njstr)//',1'//trim(pkstr)//',ishock)', &
          'l2'//trim(nistr),'my'//trim(njstr),'1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,+i,-3,-3+j,-3,-3+k, &
          'f(l2'//trim(nistr)//',my'//trim(njstr)//',mz'//trim(nkstr)//',ishock)', &
          'l2'//trim(nistr),'my'//trim(njstr),'mz'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,+i,3-j,+3,-3,-3+k, &
          'f(l2'//trim(nistr)//',1'//trim(pjstr)//',mz'//trim(nkstr)//',ishock)', &
          'l2'//trim(nistr),'1'//trim(pjstr),'mz'//trim(nkstr))




      call evaluate_integral(unitno,2,3-i,3,-j,+3,3-k,+3, &
          'f(1'//trim(pistr)//',m1'//trim(pjstr)//',1'//trim(pkstr)//',ishock)', &
          '1'//trim(pistr),'m1'//trim(pjstr),'1'//trim(pkstr))
      call evaluate_integral(unitno,2,3-i,3,-3,+j,3-k,+3, &
          'f(1'//trim(pistr)//',m2'//trim(njstr)//',1'//trim(pkstr)//',ishock)', &
          '1'//trim(pistr),'m2'//trim(njstr),'1'//trim(pkstr))
      call evaluate_integral(unitno,2,3-i,3,-3,+j,-3,-3+k, &
          'f(1'//trim(pistr)//',m2'//trim(njstr)//',mz'//trim(nkstr)//',ishock)', &
          '1'//trim(pistr),'m2'//trim(njstr),'mz'//trim(nkstr))
      call evaluate_integral(unitno,2,3-i,3,-j,+3,-3,-3+k, &
          'f(1'//trim(pistr)//',m1'//trim(pjstr)//',mz'//trim(nkstr)//',ishock)', &
          '1'//trim(pistr),'m1'//trim(pjstr),'mz'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,-3+i,-j,+3,3-k,+3, &
          'f(mx'//trim(nistr)//',m1'//trim(pjstr)//',1'//trim(pkstr)//',ishock)', &
          'mx'//trim(nistr),'m1'//trim(pjstr),'1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,-3+i,-3,+j,3-k,+3, &
          'f(mx'//trim(nistr)//',m2'//trim(njstr)//',1'//trim(pkstr)//',ishock)', &
          'mx'//trim(nistr),'m2'//trim(njstr),'1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,-3+i,-3,+j,-3,-3+k, &
          'f(mx'//trim(nistr)//',m2'//trim(njstr)//',mz'//trim(nkstr)//',ishock)', &
          'mx'//trim(nistr),'m2'//trim(njstr),'mz'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,-3+i,-j,+3,-3,-3+k, &
          'f(mx'//trim(nistr)//',m1'//trim(pjstr)//',mz'//trim(nkstr)//',ishock)', &
          'mx'//trim(nistr),'m1'//trim(pjstr),'mz'//trim(nkstr))

      call evaluate_integral(unitno,2,3-i,3,3-j,+3,-k,+3, &
          'f(1'//trim(pistr)//',1'//trim(pjstr)//',n1'//trim(pkstr)//',ishock)', &
          '1'//trim(pistr),'1'//trim(pjstr),'n1'//trim(pkstr))
      call evaluate_integral(unitno,2,3-i,3,-3,-3+j,-k,+3, &
          'f(1'//trim(pistr)//',my'//trim(njstr)//',n1'//trim(pkstr)//',ishock)', &
          '1'//trim(pistr),'my'//trim(njstr),'n1'//trim(pkstr))
      call evaluate_integral(unitno,2,3-i,3,-3,-3+j,-3,+k, &
          'f(1'//trim(pistr)//',my'//trim(njstr)//',n2'//trim(nkstr)//',ishock)', &
          '1'//trim(pistr),'my'//trim(njstr),'n2'//trim(nkstr))
      call evaluate_integral(unitno,2,3-i,3,3-j,+3,-3,+k, &
          'f(1'//trim(pistr)//',1'//trim(pjstr)//',n2'//trim(nkstr)//',ishock)', &
          '1'//trim(pistr),'1'//trim(pjstr),'n2'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,-3+i,3-j,+3,-k,+3, &
          'f(mx'//trim(nistr)//',1'//trim(pjstr)//',n1'//trim(pkstr)//',ishock)', &
          'mx'//trim(nistr),'1'//trim(pjstr),'n1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,-3+i,-3,-3+j,-k,+3, &
          'f(mx'//trim(nistr)//',my'//trim(njstr)//',n1'//trim(pkstr)//',ishock)', &
          'mx'//trim(nistr),'my'//trim(njstr),'n1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,-3+i,-3,-3+j,-3,+k, &
          'f(mx'//trim(nistr)//',my'//trim(njstr)//',n2'//trim(nkstr)//',ishock)', &
          'mx'//trim(nistr),'my'//trim(njstr),'n2'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,-3+i,3-j,+3,-3,+k, &
          'f(mx'//trim(nistr)//',1'//trim(pjstr)//',n2'//trim(nkstr)//',ishock)', &
          'mx'//trim(nistr),'1'//trim(pjstr),'n2'//trim(nkstr))

      call evaluate_integral(unitno,2,3-i,3,3-j,+3,3-k,+3, &
          'f(1'//trim(pistr)//',1'//trim(pjstr)//',1'//trim(pkstr)//',ishock)', &
          '1'//trim(pistr),'1'//trim(pjstr),'1'//trim(pkstr))
      call evaluate_integral(unitno,2,3-i,3,-3,-3+j,3-k,+3, &
          'f(1'//trim(pistr)//',my'//trim(njstr)//',1'//trim(pkstr)//',ishock)', &
          '1'//trim(pistr),'my'//trim(njstr),'1'//trim(pkstr))
      call evaluate_integral(unitno,2,3-i,3,-3,-3+j,-3,-3+k, &
          'f(1'//trim(pistr)//',my'//trim(njstr)//',mz'//trim(nkstr)//',ishock)', &
          '1'//trim(pistr),'my'//trim(njstr),'mz'//trim(nkstr))
      call evaluate_integral(unitno,2,3-i,3,3-j,+3,-3,-3+k, &
          'f(1'//trim(pistr)//',1'//trim(pjstr)//',mz'//trim(nkstr)//',ishock)', &
          '1'//trim(pistr),'1'//trim(pjstr),'mz'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,-3+i,3-j,+3,3-k,+3, &
          'f(mx'//trim(nistr)//',1'//trim(pjstr)//',1'//trim(pkstr)//',ishock)', &
          'mx'//trim(nistr),'1'//trim(pjstr),'1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,-3+i,-3,-3+j,3-k,+3, &
          'f(mx'//trim(nistr)//',my'//trim(njstr)//',1'//trim(pkstr)//',ishock)', &
          'mx'//trim(nistr),'my'//trim(njstr),'1'//trim(pkstr))
      call evaluate_integral(unitno,2,-3,-3+i,-3,-3+j,-3,-3+k, &
          'f(mx'//trim(nistr)//',my'//trim(njstr)//',mz'//trim(nkstr)//',ishock)', &
          'mx'//trim(nistr),'my'//trim(njstr),'mz'//trim(nkstr))
      call evaluate_integral(unitno,2,-3,-3+i,3-j,+3,-3,-3+k, &
          'f(mx'//trim(nistr)//',1'//trim(pjstr)//',mz'//trim(nkstr)//',ishock)', &
          'mx'//trim(nistr),'1'//trim(pjstr),'mz'//trim(nkstr))

    enddo
  enddo
!  enddo

  write(unitno,"(a,i1)") "    endsubroutine shock_calc_ext_kj_kji_",k
endsubroutine make_calc_3d_ext_kj_kji
!***********************************************************************
subroutine make_calc_externalboundary(unitno)
  use Cparam
!
  implicit none
!
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
  call declare_facefactors(unitno)
  write(unitno,"(a)") "    integer :: i,j,k"

  if ((nxgrid/=1).and.(nygrid/=1).and.(nzgrid/=1)) then
!    write(unitno,"(a)") "    real, dimension(l2i-l1i+1) :: integral1,integral2"
    write(unitno,"(a)") "!  make_calc_externalboundary: 3D CASE"
    call evaluate_facefactors(unitno,0)

    write(unitno,"(a)") " call shock_calc_ext_j_ji(f)"
    write(unitno,"(a)") " call shock_calc_ext_k_ki(f)"
    write(unitno,"(a)") " call shock_calc_ext_kj_kji_0(f)"
    write(unitno,"(a)") " call shock_calc_ext_kj_kji_1(f)"
    write(unitno,"(a)") " call shock_calc_ext_kj_kji_2(f)"

!    write(unitno,"(a)") "!  make_calc_externalboundary: 3D CASE"
!    write(unitno,"(a)") "      if (lroot) print*,'shock_calc_body: The compiled in shock profile integral is 2D'"
!    write(unitno,"(a)") "      call stop_it('shock_calc_body: This problem is 3D - INCONSITENCY')"
  elseif ((nxgrid/=1).and.(nygrid/=1)) then
    write(unitno,"(a)") "    real :: integral1,integral2"
    write(unitno,"(a)") "!  make_calc_externalboundary: 2D CASE (xy)"
    call evaluate_facefactors(unitno,0)
! Top/bottom
  write(unitno,"(a)") "!  make_calc_externalboundary: ! Top/bottom"
    do j=0,2
  write(unitno,"(a,i2)") "!  make_calc_externalboundary: ! Top/bottom: j = ",j
      jstr=itoa(j)
      if (j/=0) njstr='-'//trim(jstr)
      if (j/=0) pjstr='+'//trim(jstr)
      write(unitno,"(a)") "      do i=l1i+1,l2i-1"
      call evaluate_integral(unitno,0,-3,+3,3-j, 3,0,0, &
          'integral1','i','1 '//trim(pjstr),'n1')
      call evaluate_integral(unitno,0,-3,+3,-3,j-3,0,0, &
          'integral2','i','my'//trim(njstr),'n1')
      write(unitno,"(a)") "        f(i,1 "//trim(pjstr)//",n1,ishock)=integral1"
      write(unitno,"(a)") "        f(i,my"//trim(njstr)//",n1,ishock)=integral2"
      write(unitno,"(a)") "      enddo"
    enddo
! Left/right
  write(unitno,"(a)") "!  make_calc_externalboundary: ! Left/right"
    do i=0,2
  write(unitno,"(a,i2)") "!  make_calc_externalboundary: ! Left/right: i = ",i
      istr=itoa(i)
      if (i/=0) nistr='-'//trim(istr)
      if (i/=0) pistr='+'//trim(istr)
      write(unitno,"(a)") "      do j=m1i+1,m2i-1"
      call evaluate_integral(unitno,0,3-i, 3,-3,+3,0,0, &
          'integral1','1 '//trim(pistr),'j','n1')
      call evaluate_integral(unitno,0,-3,i-3,-3,+3,0,0, &
          'integral2','mx'//trim(nistr),'j','n1')
      write(unitno,"(a)") "        f(1 "//trim(pistr)//",j,n1,ishock)=integral1"
      write(unitno,"(a)") "        f(mx"//trim(nistr)//",j,n1,ishock)=integral2"
      write(unitno,"(a)") "      enddo"
    enddo
! Near Corners
  write(unitno,"(a)") "!  make_calc_externalboundary: ! Near Corners"
  do j=0,2
  do i=0,2
  write(unitno,"(a,i2,a,i2)") "!  make_calc_externalboundary: ! Near Corners: j = ",j,", i = ",i
  jstr=itoa(j)
  istr=itoa(i)
  if (j/=0) njstr='-'//trim(jstr)
  if (j/=0) pjstr='+'//trim(jstr)
  if (i/=0) nistr='-'//trim(istr)
  if (i/=0) pistr='+'//trim(istr)
  call evaluate_integral(unitno,2,-i,+3,3-j,+3,0,0, &
      'f(l1'//trim(pistr)//',1 '//trim(pjstr)//',n1,ishock)', &
      'l1'//trim(pistr),'1 '//trim(pjstr),'n1')
  call evaluate_integral(unitno,2,-3,+i,3-j,+3,0,0, &
      'f(l2'//trim(nistr)//',m1'//trim(pjstr)//',n1,ishock)', &
      'l2'//trim(nistr),'1 '//trim(pjstr),'n1')
  call evaluate_integral(unitno,2,-3,+i,-3,j-3,0,0, &
      'f(l2'//trim(nistr)//',m2'//trim(njstr)//',n1,ishock)', &
      'l2'//trim(nistr),'my'//trim(njstr),'n1')
  call evaluate_integral(unitno,2,-i,+3,-3,j-3,0,0, &
      'f(l1'//trim(pistr)//',m2'//trim(njstr)//',n1,ishock)', &
      'l1'//trim(pistr),'my'//trim(njstr),'n1')

  call evaluate_integral(unitno,2,3-i,+3,-j,+3,0,0, &
      'f(l1'//trim(pistr)//',m2'//trim(pjstr)//',n1,ishock)', &
      '1 '//trim(pistr),'m1'//trim(pjstr),'n1')
  call evaluate_integral(unitno,2,-3,i-3,-j,+3,0,0, &
      'f(l2'//trim(nistr)//',m1'//trim(pjstr)//',n1,ishock)', &
      'mx'//trim(nistr),'m1'//trim(pjstr),'n1')
  call evaluate_integral(unitno,2,-3,i-3,-3,+j,0,0, &
      'f(l2'//trim(nistr)//',m2'//trim(njstr)//',n1,ishock)', &
      'mx'//trim(nistr),'m2'//trim(njstr),'n1')
  call evaluate_integral(unitno,2,3-i,+3,-3,+j,0,0, &
      'f(l1'//trim(pistr)//',m2'//trim(njstr)//',n1,ishock)', &
      '1 '//trim(pistr),'m2'//trim(njstr),'n1')
  enddo
  enddo
! Corners
  write(unitno,"(a)") "!  make_calc_externalboundary: ! Corners"
    do j=0,2
    do i=0,2
  write(unitno,"(a,i2,a,i2)") "!  make_calc_externalboundary: ! Corners: j = ",j,", i = ",i
      istr=itoa(i)
      jstr=itoa(j)
      if (i/=0) nistr='-'//trim(istr)
      if (i/=0) pistr='+'//trim(istr)
      if (j/=0) njstr='-'//trim(jstr)
      if (j/=0) pjstr='+'//trim(jstr)
      call evaluate_integral(unitno,0,3-i, 3,3-j, 3,0,0, &
          'f(1 '//trim(pistr)//',1 '//trim(pjstr)//',n1,ishock)', &
          '1 '//trim(pistr),'1 '//trim(pjstr),'n1')
      call evaluate_integral(unitno,0,-3,i-3,3-j, 3,0,0, &
          'f(mx'//trim(nistr)//',1 '//trim(pjstr)//',n1,ishock)', &
          'mx'//trim(nistr),'1 '//trim(pjstr),'n1')
      call evaluate_integral(unitno,0,-3,i-3,-3,j-3,0,0, &
          'f(mx'//trim(nistr)//',my'//trim(njstr)//',n1,ishock)', &
          'mx'//trim(nistr),'my'//trim(njstr),'n1')
      call evaluate_integral(unitno,0,3-i, 3,-3,j-3,0,0, &
          'f(1 '//trim(pistr)//',my'//trim(njstr)//',n1,ishock)', &
          '1 '//trim(pistr),'my'//trim(njstr),'n1')
    enddo
    enddo
!
  elseif ((nxgrid/=1).and.(nzgrid/=1)) then
    write(unitno,"(a)") "    real :: integral1,integral2"
    write(unitno,"(a)") "!  make_calc_externalboundary: 2D CASE (xz)"
    call evaluate_facefactors(unitno,1)
! Top/bottom
    do k=0,2
      kstr=itoa(k)
      if (k/=0) nkstr='-'//trim(kstr)
      if (k/=0) pkstr='+'//trim(kstr)
      write(unitno,"(a)") "      do i=l1i+1,l2i-1"
      call evaluate_integral(unitno,1,-3,+3,0,0,3-k, 3, &
          'integral1','i','m1','1 '//trim(pkstr))
      call evaluate_integral(unitno,1,-3,+3,0,0,-3,k-3, &
          'integral2','i','m1','mz'//trim(nkstr))
      write(unitno,"(a)") "        f(i,m1,1 "//trim(pkstr)//",ishock)=integral1"
      write(unitno,"(a)") "        f(i,m1,mz"//trim(nkstr)//",ishock)=integral2"
      write(unitno,"(a)") "      enddo"
    enddo
! Left/right
    do i=0,2
      istr=itoa(i)
      if (i/=0) nistr='-'//trim(istr)
      if (i/=0) pistr='+'//trim(istr)
      write(unitno,"(a)") "      do k=n1i+1,n2i-1"
      call evaluate_integral(unitno,1,3-i, 3,0,0,-3,+3, &
          'integral1','1 '//trim(pistr),'m1','k')
      call evaluate_integral(unitno,1,-3,i-3,0,0,-3,+3, &
          'integral2','mx'//trim(nistr),'m1','k')
      write(unitno,"(a)") "        f(1 "//trim(pistr)//",m1,k,ishock)=integral1"
      write(unitno,"(a)") "        f(mx"//trim(nistr)//",m1,k,ishock)=integral2"
      write(unitno,"(a)") "      enddo"
      enddo
! Near Corners
!  do k=0,2
!  do i=0,2
!  kstr=itoa(k)
!  istr=itoa(i)
!  if (k/=0) nkstr='-'//trim(kstr)
!  if (k/=0) pkstr='+'//trim(kstr)
!  if (i/=0) nistr='-'//trim(istr)
!  if (i/=0) pistr='+'//trim(istr)
!  call evaluate_integral(unitno,2,-i,+3,0,0,3-k,+3,'f(l1'//trim(pistr)//',m1,1 '//trim(pkstr)//',ishock)','l1'//trim(pistr),'m1','1 '//trim(pkstr))
!  call evaluate_integral(unitno,2,-3,+i,0,0,3-k,+3,'f(l2'//trim(nistr)//',m1,1 '//trim(pkstr)//',ishock)','l2'//trim(nistr),'m1','1 '//trim(pkstr))
!  call evaluate_integral(unitno,2,-3,+i,0,0,-3,k-3,'f(l2'//trim(nistr)//',m1,mz'//trim(nkstr)//',ishock)','l2'//trim(nistr),'m1','mz'//trim(nkstr))
!  call evaluate_integral(unitno,2,-i,+3,0,0,-3,k-3,'f(l1'//trim(pistr)//',m1,mz'//trim(nkstr)//',ishock)','l1'//trim(pistr),'m1','mz'//trim(nkstr))
!
!  call evaluate_integral(unitno,2,3-i,+3,0,0,-k,+3,'f(1 '//trim(pistr)//',m1,n1'//trim(pkstr)//',ishock)','1 '//trim(pistr),'m1','n1'//trim(pkstr))
!  call evaluate_integral(unitno,2,-3,i-3,0,0,-k,+3,'f(mx'//trim(nistr)//',m1,n1'//trim(pkstr)//',ishock)','mx'//trim(nistr),'m1','n1'//trim(pkstr))
!  call evaluate_integral(unitno,2,-3,i-3,0,0,-3,+k,'f(mx'//trim(nistr)//',m1,n2'//trim(nkstr)//',ishock)','mx'//trim(nistr),'m1','n2'//trim(nkstr))
!  call evaluate_integral(unitno,2,3-i,+3,0,0,-3,+k,'f(1 '//trim(pistr)//',m1,n2'//trim(nkstr)//',ishock)','1 '//trim(pistr),'m1','n2'//trim(nkstr))
!  enddo
!  enddo
! Corners
    do k=0,2
    do i=0,2
      istr=itoa(i)
      kstr=itoa(k)
      if (i/=0) nistr='-'//trim(istr)
      if (i/=0) pistr='+'//trim(istr)
      if (k/=0) nkstr='-'//trim(kstr)
      if (k/=0) pkstr='+'//trim(kstr)
      call evaluate_integral(unitno,1,3-i, 3,0,0,3-k, 3, &
          'f(1 '//trim(pistr)//',m1,1 '//trim(pkstr)//',ishock)', &
          '1 '//trim(pistr),'m1','1 '//trim(pkstr))
      call evaluate_integral(unitno,1,-3,i-3,0,0,3-k, 3, &
          'f(mx'//trim(nistr)//',m1,1 '//trim(pkstr)//',ishock)', &
          'mx'//trim(nistr),'m1','1 '//trim(pkstr))
      call evaluate_integral(unitno,1,-3,i-3,0,0,-3,k-3, &
          'f(mx'//trim(nistr)//',m1,mz'//trim(nkstr)//',ishock)', &
          'mx'//trim(nistr),'m1','mz'//trim(nkstr))
      call evaluate_integral(unitno,1,3-i, 3,0,0,-3,k-3, &
          'f(1 '//trim(pistr)//',m1,mz'//trim(nkstr)//',ishock)', &
          '1 '//trim(pistr),'m1','mz'//trim(nkstr))
    enddo
    enddo
!
  elseif ((nzgrid/=1).and.(nygrid/=1)) then
    write(unitno,"(a)") "    real :: integral1,integral2"
    write(unitno,"(a)") "!  make_calc_externalboundary: 2D CASE (yz)"
    call evaluate_facefactors(unitno,2)
! Top/bottom
    do k=0,2
      kstr=itoa(k)
      if (k/=0) nkstr='-'//trim(kstr)
      if (k/=0) pkstr='+'//trim(kstr)
      write(unitno,"(a)") "      do j=m1i+1,m2i-1"
      call evaluate_integral(unitno,2,0,0,-3,+3,3-k,+3, &
          'integral1','l1','j','1 '//trim(pkstr))
      call evaluate_integral(unitno,2,0,0,-3,+3,-3,k-3, &
          'integral2','l1','j','mz'//trim(nkstr))
      write(unitno,"(a)") "        f(l1,j,1 "//trim(pkstr)//",ishock)=integral1"
      write(unitno,"(a)") "        f(l1,j,mz"//trim(nkstr)//",ishock)=integral2"
      write(unitno,"(a)") "      enddo"
    enddo
! Left/right
    do j=0,2
      jstr=itoa(j)
      if (j/=0) njstr='-'//trim(jstr)
      if (j/=0) pjstr='+'//trim(jstr)
      write(unitno,"(a)") "      do k=n1i+1,n2i-1"
      call evaluate_integral(unitno,2,0,0,3-j,+3,-3,+3, &
          'integral1','l1','1 '//trim(pjstr),'k')
      call evaluate_integral(unitno,2,0,0,-3,j-3,-3,+3, &
          'integral2','l1','my'//trim(njstr),'k')
      write(unitno,"(a)") "        f(l1,1 "//trim(pjstr)//",k,ishock)=integral1"
      write(unitno,"(a)") "        f(l1,my"//trim(njstr)//",k,ishock)=integral2"
      write(unitno,"(a)") "      enddo"
    enddo
! Near Corners
  do k=0,2
  do j=0,2
  kstr=itoa(k)
  jstr=itoa(j)
  if (k/=0) nkstr='-'//trim(kstr)
  if (k/=0) pkstr='+'//trim(kstr)
  if (j/=0) njstr='-'//trim(jstr)
  if (j/=0) pjstr='+'//trim(jstr)
  call evaluate_integral(unitno,2,0,0,-j,+3,3-k,+3, &
      'f(l1,m1'//trim(pjstr)//',1 '//trim(pkstr)//',ishock)', &
      'l1','m1'//trim(pjstr),'1 '//trim(pkstr))
  call evaluate_integral(unitno,2,0,0,-3,+j,3-k,+3, &
      'f(l1,m2'//trim(njstr)//',1 '//trim(pkstr)//',ishock)', &
      'l1','m2'//trim(njstr),'1 '//trim(pkstr))
  call evaluate_integral(unitno,2,0,0,-3,+j,-3,k-3, &
      'f(l1,m2'//trim(njstr)//',mz'//trim(nkstr)//',ishock)', &
      'l1','m2'//trim(njstr),'mz'//trim(nkstr))
  call evaluate_integral(unitno,2,0,0,-j,+3,-3,k-3, &
      'f(l1,m1'//trim(pjstr)//',mz'//trim(nkstr)//',ishock)', &
      'l1','m1'//trim(pjstr),'mz'//trim(nkstr))

!  call evaluate_integral(unitno,2,0,0,3-j,+3,-k,+3,'f(l1,1 '//trim(pjstr)//',n1'//trim(pkstr)//',ishock)','l1','1 '//trim(pjstr),'n1'//trim(pkstr))
!  call evaluate_integral(unitno,2,0,0,-3,j-3,-k,+3,'f(l1,my'//trim(njstr)//',n1'//trim(pkstr)//',ishock)','l1','my'//trim(njstr),'n1'//trim(pkstr))
!  call evaluate_integral(unitno,2,0,0,-3,j-3,-3,+k,'f(l1,my'//trim(njstr)//',n2'//trim(nkstr)//',ishock)','l1','my'//trim(njstr),'n2'//trim(nkstr))
!  call evaluate_integral(unitno,2,0,0,3-j,+3,-3,+k,'f(l1,1 '//trim(pjstr)//',n2'//trim(nkstr)//',ishock)','l1','1 '//trim(pjstr),'n2'//trim(nkstr))
  enddo
  enddo
! Corners
    do k=0,2
    do j=0,2
      kstr=itoa(k)
      jstr=itoa(j)
      if (k/=0) nkstr='-'//trim(kstr)
      if (k/=0) pkstr='+'//trim(kstr)
      if (j/=0) njstr='-'//trim(jstr)
      if (j/=0) pjstr='+'//trim(jstr)
      call evaluate_integral(unitno,2,0,0,3-j,+3,3-k,+3, &
          'f(l1,1 '//trim(pjstr)//',1 '//trim(pkstr)//',ishock)', &
          'l1','1 '//trim(pjstr),'1 '//trim(pkstr))
      call evaluate_integral(unitno,2,0,0,-3,j-3,3-k,+3, &
          'f(l1,my'//trim(njstr)//',1 '//trim(pkstr)//',ishock)', &
          'l1','my'//trim(njstr),'1 '//trim(pkstr))
      call evaluate_integral(unitno,2,0,0,-3,j-3,-3,k-3, &
          'f(l1,my'//trim(njstr)//',mz'//trim(nkstr)//',ishock)', &
          'l1','my'//trim(njstr),'mz'//trim(nkstr))
      call evaluate_integral(unitno,2,0,0,3-j,+3,-3,k-3, &
          'f(l1,1 '//trim(pjstr)//',mz'//trim(nkstr)//',ishock)', &
          'l1','1 '//trim(pjstr),'mz'//trim(nkstr))
    enddo
    enddo
!
  else
    write(unitno,"(a)") "      call stop_it('shock_calc_body: Case not implemented')"
  endif
  write(unitno,"(a)") "  endsubroutine shock_calc_externalboundary"

  if ((nxgrid/=1).and.(nygrid/=1).and.(nzgrid/=1)) then
    call make_calc_3d_ext_k_ki(unitno)
    call make_calc_3d_ext_j_ji(unitno)
    call make_calc_3d_ext_kj_kji(unitno,0)
    call make_calc_3d_ext_kj_kji(unitno,1)
    call make_calc_3d_ext_kj_kji(unitno,2)
  endif
endsubroutine make_calc_externalboundary
!***********************************************************************
subroutine evaluate_facefactors(unitno,rotation)
   use SurfaceData
!
  implicit none
!
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
  implicit none
!
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
!
  implicit none
!
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
!
  implicit none
!
!   character(len=*), optional, intent(in) :: intname, iname, jname, kname
   integer :: imin,imax,jmin,jmax,kmin,kmax
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
       If (     (points(dir1,pnt)<imin).or.(points(dir1,pnt)>imax)   &
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
   use Cparam
   use SurfaceData
!
  implicit none
!
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
