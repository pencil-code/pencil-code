! $Id: testfield_xz.f90 20934 2013-08-22 19:27:54Z mreinhardt@nordita.org $
!
!  This modules deals with all aspects of testfield fields; if no
!  testfield fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testfield relevant subroutines listed in here.
!
!  NOTE: since the fall of 2007 we have been using the routine
!  testfield_z.f90, not this one! For more information, please
!  contact Axel Brandenburg <brandenb@nordita.org>
!
!  Alex Hubbard and Matthias Rheinhardt have then developed the
!  present module from inactive/testfield.f90 rather than the
!  inactive/testfield_xz.f90 that also exists.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ltestfield = .true.
! CPARAM integer, parameter :: njtest=9
!
! MVAR CONTRIBUTION 27
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Testfield
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Testfield_general
!
  implicit none
!
  include '../testfield.h'
!
! run parameters
!
  logical :: lflucts_with_xyaver=.false.
  real    :: ktestfield_x=1., ktestfield_y=1., xx0=0., yy0=0.
!
  namelist /testfield_run_pars/       &
           B_ext,                     &
           reinitialize_aatest,       &
           lsoca,                     &
           etatest,                   &
           itestfield,                &
           ktestfield_x,              &
           ktestfield_y,              &
           xx0,                       &
           yy0,                       &
           daainit,                   &
           linit_aatest,              &
           rescale_aatest,            &
           lflucts_with_xyaver
!
! diagnostic variables
!
  integer, parameter :: n_cdiags = 44, len_cdiags = 8
  character(LEN=len_cdiags), dimension(n_cdiags) :: cdiags = &
  (/ 'alp11   ','alp21   ','alp31   ','alp12   ','alp22   ','alp32   ','alp13   ','alp23   ','alp33   ',&   ! DIAG_DOC: $\alpha_{ij}$       
     'eta111  ','eta211  ','eta311  ','eta121  ','eta221  ','eta321  ','eta131  ','eta231  ','eta331  ',&   ! DIAG_DOC: $\eta_{ijk}$
     'eta113  ','eta213  ','eta313  ','eta123  ','eta223  ','eta323  ','eta133  ','eta233  ','eta333  ',&
     'alp11cc ','alp11cs ','alp11sc ','alp11ss ','eta123cc','eta123cs','eta123sc','eta123ss',           &   ! DIAG_DOC: $\alpha_{11,\rm hh},$ 
                                                                                                            ! Diag_DOC: $\eta_{11,\rm hh}, {\rm h}={\rm c,s}$
     'E11     ','E21     ','E31     ','E12     ','E22     ','E32     ','E13     ','E23     ','E33     '  /) ! DIAG_DOC: ${\cal E}^i_j$
!
  integer, dimension(n_cdiags):: idiags=0, idiags_y=0, idiags_xy=0
  integer, parameter :: idiag_Eij_start=36, idiag_Eij_stop=idiag_Eij_start+9-1,idiag_base_end=27
!
  integer, dimension(4) :: idiag_alp11f, idiag_eta123f            
  equivalence(idiags(idiag_base_end+1),idiag_alp11f), (idiags(idiag_base_end+5),idiag_eta123f)      ! alias names for selected diagnostics
!
!  work variables
!
  real, dimension (nx)              :: cx,sx
  real, dimension (ny)              :: cy,sy
  real, dimension (nx,ny,3,3)       :: Minv
  real, dimension (nx,ny,3,njtest)  :: uxbtestm
  logical, dimension(idiag_base_end):: twod_need_1d, twod_need_2d
  logical                           :: needed2d_1d, needed2d_2d
!
  contains
!
!***********************************************************************
    subroutine initialize_testfield(f,lstarting)
!
!  Perform any post-parameter-read initialization
!
!   2-jun-05/axel: adapted from magnetic
!  27-jun-13/MR  : set itestfield='1' as it is only implemented case
!                  set lcalc_uumeanxy=.true., itestfield now string     
      use Cdata
      use Hydro, only: lcalc_uumeanxy
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, intent(in) :: lstarting
!
       real, dimension (nx)        :: c2x, cx1
       real, dimension (ny)        :: c2y, cy1
       integer :: i,j
!
!  set to zero and then rescale the testfield
!  (in future, could call something like init_aa_simple)
!
      if (reinitialize_aatest) f(:,:,:,iaatest:iaatest+ntestfield-1)=0.
!
      if (linit_aatest) lrescaling_testfield=.true.
!
      itestfield='1'
!
!  xx and yy for calculating diffusive part of emf
!
      cx=cos(ktestfield_x*(x(l1:l2)+xx0))
      sx=sin(ktestfield_x*(x(l1:l2)+xx0))
!
      cy=cos(ktestfield_y*(y(m1:m2)+yy0))
      sy=sin(ktestfield_y*(y(m1:m2)+yy0))
!
      !c2x=cos(2*ktestfield_x*(x(l1:l2)+xx0))
      !c2z=cos(2*ktestfield_z*(z(n1:n2)+zz0))
!
      cx1=1./cx
      cy1=1./cy
!
      do i=1,nx
      do j=1,ny
!
!        Minv(i,j,1,:) = (/ 0.5*(c2x(i)+c2z(j))*cx1(i)*cz1(j),              sx(i)*cz1(j),              sz(j)*cx1(i) /)
        Minv(i,j,1,:) = (/ (1.- sx(i)**2 - sy(j)**2)*cx1(i)*cy1(j),              sx(i)*cy1(j),              sy(j)*cx1(i) /)
        Minv(i,j,2,:) = (/              -sx(i)*cy1(j)/ktestfield_x, cx(i)*cy1(j)/ktestfield_x,                        0. /)
        Minv(i,j,3,:) = (/              -sy(j)*cx1(i)/ktestfield_y,                        0., cy(j)*cx1(i)/ktestfield_y /)
!
!        Minv(i,j,:,:) = RESHAPE((/  &
!                                  (/ (1.- sx(i)**2 - sz(j)**2)*cx1(i)*cz1(j),        sx(i)*cz1(j),              sz(j)*cx1(i) /),&
!                                  (/              -sx(i)*cz1(j)/ktestfield_x, cx(i)*cz1(j)/ktestfield_x,                  0. /),&
!                                  (/              -sz(j)*cx1(i)/ktestfield_z,                  0., cz(j)*cx1(i)/ktestfield_z /) &
!                                /), (/3,3/), ORDER=(/ 2, 1 /))
      enddo
      enddo
!
      lcalc_uumeanxy = .true.
!
!  write testfield information to a file (for convenient post-processing)
!
      if (lroot) then
        open(1,file=trim(datadir)//'/testfield_info.dat',STATUS='unknown')
        write(1,'(a,i1)') 'lsoca='  ,merge(1,0,lsoca)
        write(1,'(a,a)') 'itestfield=',itestfield
        write(1,'(a,2(f3.0))') 'ktestfield_x,y=', ktestfield_x, ktestfield_y
        close(1)
      endif
!
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_testfield
!***********************************************************************
    subroutine read_testfield_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=testfield_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=testfield_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_testfield_run_pars
!***********************************************************************
    subroutine write_testfield_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=testfield_run_pars)
!
    endsubroutine write_testfield_run_pars
!***********************************************************************
    subroutine daatest_dt(f,df,p)
!
!  testfield evolution
!  calculate da^(q)/dt=uxB^(q)+eta*del2A^(q), where q=1,...,9
!
!   3-jun-05/axel: coded
!  27-jun-13/MR  : correct calculation of uufluct intro'd
!                  moved calculation of xy-averaged quantities to 
!                  calc_coefficients, completed
!
      use Cdata
      use Diagnostics
      use Mpicomm, only: stop_it
      use Sub
      use Hydro, only: lcalc_uumeanxy, uumxy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: uxB,bbtest,btest,uxbtest,uufluct
      real, dimension (nx,3) :: del2Atest
      integer :: jtest,j,i
!
      intent(in)     :: f,p
      intent(inout)  :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'daatest_dt: SOLVE'
!
      if (headtt) then
        if (iaxtest /= 0) call identify_bcs('Axtest',iaxtest)
        if (iaytest /= 0) call identify_bcs('Aytest',iaytest)
        if (iaztest /= 0) call identify_bcs('Aztest',iaztest)
      endif
!
!  do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further down in the file.
!
      do jtest=1,njtest
!
        iaxtest=iaatest+3*(jtest-1)
        iaztest=iaxtest+2
!
        call del2v(f,iaxtest,del2Atest)
!
        select case (itestfield)
          case ('1'); call set_bbtest (bbtest,jtest)
          case ('2'); call set_bbtest2(bbtest,jtest)
          case ('3'); call set_bbtest3(bbtest,jtest)
          case ('4'); call set_bbtest4(bbtest,jtest)
          case default; call fatal_error('daatest_dt','undefined itestfield')
        endselect
!
!  add an external field, if present
!
        if (B_ext(1)/=0.) bbtest(:,1)=bbtest(:,1)+B_ext(1)
        if (B_ext(2)/=0.) bbtest(:,2)=bbtest(:,2)+B_ext(2)
        if (B_ext(3)/=0.) bbtest(:,3)=bbtest(:,3)+B_ext(3)
!
        uufluct = p%uu-uumxy(l1:l2,n,:)

        call cross_mn(uufluct,bbtest,uxB)
        df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest)+etatest*del2Atest+uxB 
!
        if (.not.lsoca) then

          call curl(f,iaxtest,btest)
          call cross_mn(uufluct,btest,uxbtest)
!
!  subtract average emf
!
          df(l1:l2,m,n,iaxtest:iaztest)= df(l1:l2,m,n,iaxtest:iaztest) &
                                        +uxbtest-uxbtestm(:,n-n1+1,:,jtest)
        endif
!
!  diffusive time step, just take the max of diffus_eta (if existent)
!  and whatever is calculated here
!
        if (lfirst.and.ldt) diffus_eta=max(diffus_eta,etatest*dxyz_2)
      
      enddo
!
    endsubroutine daatest_dt
!***********************************************************************
    subroutine get_slices_testfield(f,slices)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_testfield
!***********************************************************************
    subroutine testfield_after_boundary(f,p)
!
!  calculate <uxb>, which is needed when lsoca=.false.
!
!  21-jan-06/axel: coded
!
      use Cdata
      use Sub
      use Hydro, only: calc_pencils_hydro
      use Mpicomm, only: stop_it,mpiallreduce_sum     !,mpiallreduce_sum_arr
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: btest,uxbtest
      integer :: jtest,j,nxy=nxgrid*nygrid, nl
      logical :: headtt_save
      real :: fac
      real, dimension (nx,ny,3) :: uxbtestm1
!
      intent(in) :: f
!
      logical :: need_output
!
!  In this routine we will reset headtt after the first pencil,
!  so we need to reset it afterwards.
!
      headtt_save=headtt
      fac=1./ny
      if ((ldiagnos .and. needed2d_1d) .or. &
          (l2davgfirst .and. needed2d_2d)) then
        need_output=.true. ; else ; need_output=.false. ; endif
!
!  do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further up in the file.
!
      do jtest=1,njtest
!
        iaxtest=iaatest+3*(jtest-1)
        iaztest=iaxtest+2

        if (lsoca .and.(.not. need_output)) then
          uxbtestm(:,:,:,jtest)=0.
        else
!
          do n=n1,n2
!
            nl=n-n1+1
!
            uxbtestm(:,nl,:,jtest)=0.
!
            do m=m1,m2
!
              call curl(f,iaxtest,btest)
              call cross_mn(f(l1:l2,m,n,iux:iuz),btest,uxbtest)
!
!  without SOCA, the alpha tensor is anisotropic even for the standard
!  Roberts flow. To check that this is because of the averaging that
!  enters, we allow ourselves to turn it off with lflucts_with_xyaver=.true.
!  It is off by default.
!
              do j=1,3
                if (lflucts_with_xyaver) then
                  uxbtestm(:,nl,j,jtest)=spread(sum( &
                    uxbtestm(:,nl,j,jtest)+fac*uxbtest(:,j),1),1,nx)/nx
                else
                  uxbtestm(:,nl,j,jtest)= &
                    uxbtestm(:,nl,j,jtest)+fac*uxbtest(:,j)
                endif
              enddo
              headtt=.false.
            enddo
          enddo
!
!  do communication along z
!
          call mpiallreduce_sum(uxbtestm(:,:,:,jtest),uxbtestm1,(/nx,ny,3/),idir=3)
!         or
!         call mpiallreduce_sum_arr(uxbtestm(1,1,1,jtest),uxbtestm1,nx*nz*3,idir=2)       !avoids copy
          uxbtestm(:,:,:,jtest)=uxbtestm1/nprocz
!
        endif
      enddo
!
!  reset headtt
!
      headtt=headtt_save
!
      if (need_output) call calc_coefficients
!
    endsubroutine testfield_after_boundary
!***********************************************************************
    subroutine calc_coefficients
!
!  26-feb-13/MR: determination of y-averaged components of alpha completed
!   6-mar-13/MR: internal order of etas changed; calls to save_name, *sum_mn_name 
!                simplified
!   7-mar-13/MR: further shortened by introduction of do-loops in calculating
!                temp_array
!  27-jun-13/MR: avoided calculation of pencil-case, introduced calculation of mean EMF
!   
    Use Diagnostics
    Use Sub, only: fourier_single_mode
!
    integer :: i, j, ij, count, nl, jtest
    real, dimension(2,nx) :: temp_fft_y
    real, dimension(2,2) :: temp_fft
    real, dimension (:,:,:), allocatable :: temp_array
    logical, dimension(idiag_base_end) :: need_temp
    integer, dimension (idiag_base_end) :: twod_address
!
! Mean EMF
!
    do n=n1,n2 
!
      nl = n-n1+1
      jtest = 1
!
      do i=idiag_Eij_start,idiag_Eij_stop,3                            ! over all testfields
        do j=1,3                                                       ! over vector components
!          call ysum_mn_name_xz(uxbtestm(:,nl,j,jtest),idiags_xz(i+j-1))
          call zsum_mn_name_xy(uxbtestm(:,nl,j,jtest),idiags_xy(i+j-1))
        enddo
        jtest = jtest+1
      enddo
!
    enddo

    if (l2davgfirst .and. needed2d_2d) need_temp=twod_need_2d
    if (ldiagnos .and. needed2d_1d) need_temp=need_temp .or. twod_need_1d
!
    count=0
    twod_address=1     !to ensure valid indices even when a slot is unused (flagged by need_temp)

    do j=1,idiag_base_end
      if (need_temp(j)) then
        count=count+1
        twod_address(j)=count
      endif
    enddo
!
    if (count==0) return

    allocate(temp_array(nx,ny,count))
!
    select case (itestfield)
    case ('1')
!
      do n=1,ny
        do i=1,3 
!
          if (need_temp(i))   & !(idiag_alpi1*/=0) &
            temp_array(:,n,twod_address(i))= &
                Minv(:,n,1,1)*uxbtestm(:,n,i,1)+Minv(:,n,1,2)*uxbtestm(:,n,i,2)+Minv(:,n,1,3)*uxbtestm(:,n,i,3)
!
          if (need_temp(3+i)) & !(idiag_alpi2*/=0) &
            temp_array(:,n,twod_address(3+i))= &
                Minv(:,n,1,1)*uxbtestm(:,n,i,4)+Minv(:,n,1,2)*uxbtestm(:,n,i,5)+Minv(:,n,1,3)*uxbtestm(:,n,i,6)
!
          if (need_temp(6+i)) & !(idiag_alpi3*/=0) &
            temp_array(:,n,twod_address(6+i))= &
                Minv(:,n,1,1)*uxbtestm(:,n,i,7)+Minv(:,n,1,2)*uxbtestm(:,n,i,8)+Minv(:,n,1,3)*uxbtestm(:,n,i,9)
!
          if (need_temp(9+i)) & !(idiag_etai11*/=0) &
            temp_array(:,n,twod_address(9+i))= &
                Minv(:,n,2,1)*uxbtestm(:,n,i,1)+Minv(:,n,2,2)*uxbtestm(:,n,i,2)+Minv(:,n,2,3)*uxbtestm(:,n,i,3)
!
          if (need_temp(12+i)) & !(idiag_etai21*/=0) &
            temp_array(:,n,twod_address(12+i))= &
                Minv(:,n,2,1)*uxbtestm(:,n,i,4)+Minv(:,n,2,2)*uxbtestm(:,n,i,5)+Minv(:,n,2,3)*uxbtestm(:,n,i,6)
!
          if (need_temp(15+i)) & !(idiag_etai31*/=0) &
            temp_array(:,n,twod_address(15+i))= &
                Minv(:,n,2,1)*uxbtestm(:,n,i,7)+Minv(:,n,2,2)*uxbtestm(:,n,i,8)+Minv(:,n,2,3)*uxbtestm(:,n,i,9)
!
          if (need_temp(18+i)) & !(idiag_etai13*/=0) &
            temp_array(:,n,twod_address(18+i))= &
                Minv(:,n,3,1)*uxbtestm(:,n,i,1)+Minv(:,n,3,2)*uxbtestm(:,n,i,2)+Minv(:,n,3,3)*uxbtestm(:,n,i,3)
!
          if (need_temp(21+i)) & !(idiag_etai23*/=0) &
            temp_array(:,n,twod_address(21+i))= &
                Minv(:,n,3,1)*uxbtestm(:,n,i,4)+Minv(:,n,3,2)*uxbtestm(:,n,i,5)+Minv(:,n,3,3)*uxbtestm(:,n,i,6)
!
          if (need_temp(24+i)) & !(idiag_etai33*/=0) &
            temp_array(:,n,twod_address(24+i))= &
                Minv(:,n,3,1)*uxbtestm(:,n,i,7)+Minv(:,n,3,2)*uxbtestm(:,n,i,8)+Minv(:,n,3,3)*uxbtestm(:,n,i,9)
        enddo
      enddo
!
    case default
      call fatal_error('calc_coefficients','Calculation of coefficients not implemented for itestfield /= 1')
      temp_array=0.
!
    end select
!
    if (ldiagnos .and. needed2d_1d) then
!
      if (any(idiag_alp11f/=0)) then
        call fourier_single_mode(temp_array(:,:,twod_address(1)), &
            (/nx,ny/), 1., 3, temp_fft_y, l2nd=.true.)
        if (lroot) then
          call fourier_single_mode(temp_fft_y, (/2,nx/), 1., 1, temp_fft, l2nd=.true.)
          ij=1
          do i=1,2
            do j=1,2
              call save_name(temp_fft(i,j), idiag_alp11f(ij))
              ij=ij+1
            enddo
          enddo 
        endif
      endif
!
      if (any(idiag_eta123f/=0)) then
        call fourier_single_mode(temp_array(:,:,twod_address(22)), &
            (/nx,ny/), 1., 3, temp_fft_y, l2nd=.true.)
        if (lroot) then
          call fourier_single_mode(temp_fft_y, (/2,nx/), 1., 1, temp_fft, l2nd=.true.)
          ij=1
          do i=1,2
            do j=1,2
              call save_name(temp_fft(i,j), idiag_eta123f(ij))
              ij=ij+1
            enddo
          enddo 
        endif
      endif
!
    endif
!
    temp_array = ny*temp_array    ! ny multiplied because we are in the following only in an n loop
!
    if (ldiagnos .and. needed2d_1d) then
      do n=n1,n2
        
        nl = n-n1+1
        do i=1,size(twod_address)
          call sum_mn_name(temp_array(:,nl,twod_address(i)), idiags(i))
        enddo

      enddo
    endif
!
    if (l1davgfirst .and. needed2d_1d) then  !!!TBC
      do n=n1,n2
!
        nl = n-n1+1
        do i=1,size(twod_address)
          call yzsum_mn_name_x(temp_array(:,nl,twod_address(i)), idiags_y(i))
        enddo
!
      enddo
    endif
!
    if (l2davgfirst .and. needed2d_2d) then
      do n=n1,n2 
!
        nl = n-n1+1
        do i=1,size(twod_address)
          call zsum_mn_name_xy(temp_array(:,nl,twod_address(i)),idiags_xy(i))
        enddo
!
      enddo
    endif
!
    deallocate(temp_array)
!
    endsubroutine calc_coefficients
!***********************************************************************
    subroutine set_bbtest(bbtest,jtest)
!
!  set testfield
!
!   3-jun-05/axel: coded
!
      use Cdata
      use Sub
!
      real, dimension (nx,3) :: bbtest
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: bbtest
!
      integer :: nl
!
!  set bbtest for each of the 9 cases
!
      nl = n-n1+1
!
      select case (jtest)
!
      case (1); bbtest(:,1)=cx*cy(nl); bbtest(:,2)=0.; bbtest(:,3)=0.
      case (2); bbtest(:,1)=sx*cy(nl); bbtest(:,2)=0.; bbtest(:,3)=0.
      case (3); bbtest(:,1)=cx*sy(nl); bbtest(:,2)=0.; bbtest(:,3)=0.
!
      case (4); bbtest(:,1)=0.; bbtest(:,2)=cx*cy(nl); bbtest(:,3)=0.
      case (5); bbtest(:,1)=0.; bbtest(:,2)=sx*cy(nl); bbtest(:,3)=0.
      case (6); bbtest(:,1)=0.; bbtest(:,2)=cx*sy(nl); bbtest(:,3)=0.
!
      case (7); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cx*cy(nl)
      case (8); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sx*cy(nl)
      case (9); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cx*sy(nl)
!
      case default; bbtest(:,:)=0.
!
      endselect
!
    endsubroutine set_bbtest
!***********************************************************************
    subroutine set_bbtest2(bbtest,jtest)
!
!  set alternative testfield
!
!  10-jun-05/axel: adapted from set_bbtest
!
      use Cdata
      use Sub
!
      real, dimension (nx,3) :: bbtest
      real, dimension (nx) :: cx,sx,cy,sy,xy
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: bbtest
!
!  xx and yy for calculating diffusive part of emf
!
      cx=cos(x(l1:l2))
      sx=sin(x(l1:l2))
      cy=cos(y(m))
      sy=sin(y(m))
      xy=cx*cy
!
!  set bbtest for each of the 9 cases
!
      select case (jtest)
      case (1); bbtest(:,1)=xy; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (2); bbtest(:,1)=0.; bbtest(:,2)=xy; bbtest(:,3)=0.
      case (3); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=xy
      case (4); bbtest(:,1)=sy; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (5); bbtest(:,1)=0.; bbtest(:,2)=sy; bbtest(:,3)=0.
      case (6); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sy
      case (7); bbtest(:,1)=sx; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (8); bbtest(:,1)=0.; bbtest(:,2)=sx; bbtest(:,3)=0.
      case (9); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sx
      case default; bbtest(:,:)=0.
      endselect
!
    endsubroutine set_bbtest2
!***********************************************************************
    subroutine set_bbtest4(bbtest,jtest)
!
!  set testfield using constant and linear functions
!
!  15-jun-05/axel: adapted from set_bbtest3
!
      use Cdata
      use Sub
!
      real, dimension (nx,3) :: bbtest
      real, dimension (nx) :: xx,yy
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: bbtest
!
!  xx and yy for calculating diffusive part of emf
!
      xx=x(l1:l2)
      yy=y(m)
!
!  set bbtest for each of the 9 cases
!
      select case (jtest)
      case (1); bbtest(:,1)=1.; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (2); bbtest(:,1)=0.; bbtest(:,2)=1.; bbtest(:,3)=0.
      case (3); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=1.
      case (4); bbtest(:,1)=yy; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (5); bbtest(:,1)=0.; bbtest(:,2)=yy; bbtest(:,3)=0.
      case (6); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=yy
      case (7); bbtest(:,1)=xx; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (8); bbtest(:,1)=0.; bbtest(:,2)=xx; bbtest(:,3)=0.
      case (9); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=xx
      case default; bbtest(:,:)=0.
      endselect
!
    endsubroutine set_bbtest4
!***********************************************************************
    subroutine set_bbtest3(bbtest,jtest)
!
!  set alternative testfield
!
!  10-jun-05/axel: adapted from set_bbtest
!
      use Cdata
      use Sub
!
      real, dimension (nx,3) :: bbtest
      real, dimension (nx) :: cx,sx,cy,sy
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: bbtest
!
!  xx and yy for calculating diffusive part of emf
!
      cx=cos(x(l1:l2))
      sx=sin(x(l1:l2))
      cy=cos(y(m))
      sy=sin(y(m))
!
!  set bbtest for each of the 9 cases
!
      select case (jtest)
      case (1); bbtest(:,1)=1.; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (2); bbtest(:,1)=0.; bbtest(:,2)=1.; bbtest(:,3)=0.
      case (3); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=1.
      case (4); bbtest(:,1)=sy; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (5); bbtest(:,1)=0.; bbtest(:,2)=sy; bbtest(:,3)=0.
      case (6); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cy
      case (7); bbtest(:,1)=cx; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (8); bbtest(:,1)=0.; bbtest(:,2)=sx; bbtest(:,3)=0.
      case (9); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sx
      case default; bbtest(:,:)=0.
      endselect
!
    endsubroutine set_bbtest3
!***********************************************************************
    subroutine rprint_testfield(lreset,lwrite)
!
!  reads and registers print parameters relevant for testfield fields
!
!   3-jun-05/axel: adapted from rprint_magnetic
!  26-feb-13/MR  : output of ntestfield in index.pro added
!   6-mar-13/MR  : alternative parse_name used
!
      use Cdata
      use Diagnostics
      use Sub, only: loptest
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname,inamey,inamexy,i
!
!  reset everything in case of RELOAD
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
!
        idiags=0; idiags_y=0; idiags_xy=0
!
      endif
!
!  check for those quantities that we want to evaluate online
!
!  volume averages
!
      do iname=1,nname
        do i=1,size(idiags)
          call parse_name(iname,cname,cform,cdiags(i),idiags(i))
        enddo
      enddo
!
!  xy-averages
!
      do inamey=1,nnamey
        do i=1,size(idiags_y)
          call parse_name(inamey,cnamey,cformy,trim(cdiags(i))//'y',idiags_y(i))
        enddo
      enddo
!
!  y-averages
!
      do inamexy=1,nnamexy
        do i=1,size(idiags_xy)
          call parse_name(inamexy,cnamexy,cformxy,trim(cdiags(i))//'xy',idiags_xy(i))
        enddo
      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!
      if (loptest(lwrite)) then
!
        write(3,*) 'iaatest=',iaatest
        write(3,*) 'ntestfield=',ntestfield
        write(3,*) 'nnamey=',nnamey
        write(3,*) 'nnamexy=',nnamexy
        write(3,*) 'nnamexy=',nnamexy
!
      endif
!
      call diagnos_interdep
!
    endsubroutine rprint_testfield
!***********************************************************************
    subroutine diagnos_interdep
!
! 2d-dependences
!
    twod_need_2d = idiags_xy(1:idiag_base_end)/=0
!
    needed2d_2d = any(twod_need_2d)
!
!  2d dependencies of 0 or 1-d averages
!
    twod_need_1d = idiags(1:idiag_base_end)/=0 .or. idiags_y(1:idiag_base_end)/=0
!
    needed2d_1d=any(twod_need_1d)
!
  endsubroutine diagnos_interdep
!***********************************************************************
endmodule Testfield
