! $Id$
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

module Testfield

  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../testfield.h'

  character (len=labellen) :: initaatest='zero'

  ! input parameters
  real, dimension(3) :: B_ext=(/0.,0.,0./)
  real, dimension (nx,3) :: bbb
  real :: amplaa=0., kx_aa=1.,ky_aa=1.,kz_aa=1.
  real :: taainit=0.,daainit=0.
  logical :: reinitialize_aatest=.false.
  logical :: xextent=.true.,zextent=.true.,lsoca=.false.,lset_bbtest2=.false.
  logical :: linit_aatest=.false., lflucts_with_xyaver=.false.
  integer :: itestfield=1
  real :: ktestfield_x=1., ktestfield_z=1., xx0=0., zz0=0., kanalyze_x=-1.
  real :: kanalyze_z=-1.
  integer, parameter :: mtestfield=3*njtest
  integer :: naainit

  namelist /testfield_init_pars/ &
       B_ext,xextent,zextent,initaatest

  ! run parameters
  real :: etatest=0.
  real, dimension(njtest) :: rescale_aatest=0.
  namelist /testfield_run_pars/ &
       B_ext,reinitialize_aatest,xextent,zextent,lsoca, &
       lset_bbtest2,etatest,itestfield,ktestfield_x, &
       ktestfield_z,xx0,zz0,kanalyze_x,kanalyze_z,daainit, &
       linit_aatest, lflucts_with_xyaver, &
       rescale_aatest

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_alp11=0,idiag_alp21=0,idiag_alp31=0
  integer :: idiag_alp12=0,idiag_alp22=0,idiag_alp32=0
  integer :: idiag_alp13=0,idiag_alp23=0,idiag_alp33=0
  integer :: idiag_eta111=0,idiag_eta211=0,idiag_eta311=0
  integer :: idiag_eta112=0,idiag_eta212=0,idiag_eta312=0
  integer :: idiag_eta121=0,idiag_eta221=0,idiag_eta321=0
  integer :: idiag_eta122=0,idiag_eta222=0,idiag_eta322=0
  integer :: idiag_eta131=0,idiag_eta231=0,idiag_eta331=0
  integer :: idiag_eta132=0,idiag_eta232=0,idiag_eta332=0
  integer :: idiag_alp11cc=0,idiag_alp11cs=0,idiag_alp11sc=0,idiag_alp11ss=0
  integer :: idiag_eta122cc=0,idiag_eta122cs=0,idiag_eta122sc=0,idiag_eta122ss=0
  integer :: idiag_alp11z=0,idiag_alp21z=0,idiag_alp31z=0
  integer :: idiag_alp12z=0,idiag_alp22z=0,idiag_alp32z=0
  integer :: idiag_alp13z=0,idiag_alp23z=0,idiag_alp33z=0
  integer :: idiag_eta111z=0,idiag_eta211z=0,idiag_eta311z=0
  integer :: idiag_eta121z=0,idiag_eta221z=0,idiag_eta321z=0
  integer :: idiag_eta131z=0,idiag_eta231z=0,idiag_eta331z=0
  integer :: idiag_eta113z=0,idiag_eta213z=0,idiag_eta313z=0
  integer :: idiag_eta123z=0,idiag_eta223z=0,idiag_eta323z=0
  integer :: idiag_eta133z=0,idiag_eta233z=0,idiag_eta333z=0
  integer :: idiag_alp11xz=0,idiag_alp21xz=0,idiag_alp31xz=0
  integer :: idiag_alp12xz=0,idiag_alp22xz=0,idiag_alp32xz=0
  integer :: idiag_alp13xz=0,idiag_alp23xz=0,idiag_alp33xz=0
  integer :: idiag_eta111xz=0,idiag_eta211xz=0,idiag_eta311xz=0
  integer :: idiag_eta112xz=0,idiag_eta212xz=0,idiag_eta312xz=0
  integer :: idiag_eta121xz=0,idiag_eta221xz=0,idiag_eta321xz=0
  integer :: idiag_eta122xz=0,idiag_eta222xz=0,idiag_eta322xz=0
  integer :: idiag_eta131xz=0,idiag_eta231xz=0,idiag_eta331xz=0
  integer :: idiag_eta132xz=0,idiag_eta232xz=0,idiag_eta332xz=0
  integer :: idiag_eta113xz=0,idiag_eta213xz=0,idiag_eta313xz=0
  integer :: idiag_eta123xz=0,idiag_eta223xz=0,idiag_eta323xz=0
  integer :: idiag_eta133xz=0,idiag_eta233xz=0,idiag_eta333xz=0
  integer :: idiag_alp11exz=0,idiag_alp21exz=0,idiag_alp31exz=0
  integer :: idiag_alp12exz=0,idiag_alp22exz=0,idiag_alp32exz=0
  integer :: idiag_alp13exz=0,idiag_alp23exz=0,idiag_alp33exz=0

  real, dimension (nx,nz,3,njtest) :: uxbtestm
  real, dimension (nx)        :: cx,sx
  real, dimension (nz)        :: cz,sz
  real, dimension (nx,nz,3,3) :: Minv
  logical :: needed2d_1d, needed2d_2d
  logical, dimension(27) :: twod_need_1d, twod_need_2d

  contains

!***********************************************************************
    subroutine register_testfield()
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaatest, etc; increase nvar accordingly
!
!   3-jun-05/axel: adapted from register_magnetic
!
      use Cdata
      use Mpicomm
      use Sub
!
      integer :: j
!
!  Set first and last index of text field
!  Note: iaxtest, iaytest, and iaztest are initialized to the first test field.
!  These values are used in this form in start, but later overwritten.
!
      iaatest=nvar+1
      iaxtest=iaatest
      iaytest=iaatest+1
      iaztest=iaatest+2
      ntestfield=mtestfield
      nvar=nvar+ntestfield
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_testfield: nvar = ', nvar
        print*, 'register_testfield: iaatest = ', iaatest
      endif
!
!  Put variable names in array
!
      do j=iaatest,nvar
        varname(j) = 'aatest'
      enddo
      iaztestpq=nvar
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id$")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_testfield: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',aatest $'
          if (nvar == mvar) write(4,*) ',aatest'
        else
          write(4,*) ',aatest $'
        endif
        write(15,*) 'aatest = fltarr(mx,my,mz,ntestfield)*one'
      endif
!
    endsubroutine register_testfield
!***********************************************************************
    subroutine initialize_testfield(f,lstarting)
!
!  Perform any post-parameter-read initialization
!
!   2-jun-05/axel: adapted from magnetic
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, intent(in) :: lstarting
!
       real, dimension (nx)        :: c2x, cx1
       real, dimension (nz)        :: c2z, cz1
       integer :: i,j
!
!  set to zero and then rescale the testfield
!  (in future, could call something like init_aa_simple)
!
      if (reinitialize_aatest) then
        f(:,:,:,iaatest:iaatest+ntestfield-1)=0.
      endif
!
      if (linit_aatest) then
        lrescaling_testfield=.true.
      endif
!
!  xx and zz for calculating diffusive part of emf
!
      cx=cos(ktestfield_x*(x(l1:l2)+xx0))
      sx=sin(ktestfield_x*(x(l1:l2)+xx0))

      cz=cos(ktestfield_z*(z(n1:n2)+zz0))
      sz=sin(ktestfield_z*(z(n1:n2)+zz0))

      !c2x=cos(2*ktestfield_x*(x(l1:l2)+xx0))
      !c2z=cos(2*ktestfield_z*(z(n1:n2)+zz0))

      cx1=1./cx
      cz1=1./cz

      do i=1,nx
      do j=1,nz

!        Minv(i,j,1,:) = (/ 0.5*(c2x(i)+c2z(j))*cx1(i)*cz1(j),              sx(i)*cz1(j),              sz(j)*cx1(i) /)
        Minv(i,j,1,:) = (/ (1.- sx(i)**2 - sz(j)**2)*cx1(i)*cz1(j),              sx(i)*cz1(j),              sz(j)*cx1(i) /)
        Minv(i,j,2,:) = (/              -sx(i)*cz1(j)/ktestfield_x, cx(i)*cz1(j)/ktestfield_x,                        0. /)
        Minv(i,j,3,:) = (/              -sz(j)*cx1(i)/ktestfield_z,                        0., cz(j)*cx1(i)/ktestfield_z /)
!
!        Minv(i,j,:,:) = RESHAPE((/  &
!                                  (/ (1.- sx(i)**2 - sz(j)**2)*cx1(i)*cz1(j),        sx(i)*cz1(j),              sz(j)*cx1(i) /),&
!                                  (/              -sx(i)*cz1(j)/ktestfield_x, cx(i)*cz1(j)/ktestfield_x,                  0. /),&
!                                  (/              -sz(j)*cx1(i)/ktestfield_z,                  0., cz(j)*cx1(i)/ktestfield_z /) &
!                                /), (/3,3/), ORDER=(/ 2, 1 /))
      enddo
      enddo
!
!  write testfield information to a file (for convenient post-processing)
!
      if (lroot) then
        open(1,file=trim(datadir)//'/testfield_info.dat',STATUS='unknown')
        write(1,'(a,i1)') 'xextent=',merge(1,0,xextent)
        write(1,'(a,i1)') 'zextent=',merge(1,0,zextent)
        write(1,'(a,i1)') 'lsoca='  ,merge(1,0,lsoca)
        write(1,'(a,i2)') 'itestfield=',itestfield
        write(1,'(a,2(f3.0))') 'ktestfield_x,z=', ktestfield_x, ktestfield_z
        close(1)
      endif
!
      call keep_compiler_quiet(lstarting)
    endsubroutine initialize_testfield
!***********************************************************************
    subroutine init_aatest(f)
!
!  initialise testfield; called from start.f90
!
!   2-jun-05/axel: adapted from magnetic
!
      use Cdata
      use Mpicomm
      use Gravity, only: gravz
      use Sub
      use Initcond
      use InitialCondition, only: initial_condition_aatest
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      select case (initaatest)

      case ('zero', '0'); f(:,:,:,iaatest:iaatest+ntestfield-1)=0.

      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'init_aatest: check initaatest: ', trim(initaatest)
        call stop_it("")

      endselect
!
!  Interface for user's own subroutine
!
      if (linitial_condition) call initial_condition_aatest(f)
!
    endsubroutine init_aatest
!***********************************************************************
    subroutine pencil_criteria_testfield()
!
!   All pencils that the Testfield module depends on are specified here.
!
!  26-jun-05/anders: adapted from magnetic
!
      use Cdata
!
      lpenc_requested(i_uu)=.true.
!
    endsubroutine pencil_criteria_testfield
!***********************************************************************
    subroutine pencil_interdep_testfield(lpencil_in)
!
!  Interdependency among pencils from the Testfield module is specified here.
!
!  26-jun-05/anders: adapted from magnetic
!
      use Cdata
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_testfield
!***********************************************************************
    subroutine read_testfield_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=testfield_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=testfield_init_pars,ERR=99)
      endif

99    return
    endsubroutine read_testfield_init_pars
!***********************************************************************
    subroutine write_testfield_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=testfield_init_pars)

    endsubroutine write_testfield_init_pars
!***********************************************************************
    subroutine read_testfield_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=testfield_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=testfield_run_pars,ERR=99)
      endif

99    return
    endsubroutine read_testfield_run_pars
!***********************************************************************
    subroutine write_testfield_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=testfield_run_pars)

    endsubroutine write_testfield_run_pars
!***********************************************************************
    subroutine daatest_dt(f,df,p)
!
!  testfield evolution
!  calculate da^(q)/dt=uxB^(q)+eta*del2A^(q), where q=1,...,9
!
!   3-jun-05/axel: coded
!
      use Cdata
      use Diagnostics
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      real, dimension (nx,3) :: uxB,bbtest,btest,uxbtest,duxbtest
      real, dimension (nx,3) :: del2Atest
      integer :: jtest,j
!
      intent(in)     :: f,p
      intent(inout)  :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'daatest_dt: SOLVE'
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

          iaxtest=iaatest+3*(jtest-1)
          iaztest=iaxtest+2

          call del2v(f,iaxtest,del2Atest)

          select case (itestfield)
          case (1); call set_bbtest(bbtest,jtest)
          case (2); call set_bbtest2(bbtest,jtest)
          case (3); call set_bbtest3(bbtest,jtest)
          case (4); call set_bbtest4(bbtest,jtest)
          case default; call fatal_error('daatest_dt','undefined itestfield')

          endselect
!
!  add an external field, if present
!
          if (B_ext(1)/=0.) bbtest(:,1)=bbtest(:,1)+B_ext(1)
          if (B_ext(2)/=0.) bbtest(:,2)=bbtest(:,2)+B_ext(2)
          if (B_ext(3)/=0.) bbtest(:,3)=bbtest(:,3)+B_ext(3)
!
          call cross_mn(p%uu,bbtest,uxB)

          if (lsoca) then
            df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest) &
              +uxB+etatest*del2Atest
          else
            call curl(f,iaxtest,btest)
            call cross_mn(p%uu,btest,uxbtest)
!
!  subtract average emf
!
            do j=1,3
              duxbtest(:,j)=uxbtest(:,j)-uxbtestm(:,n-n1+1,j,jtest)
            enddo
!
!  advance test field equation
!
            df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest) &
              +uxB+etatest*del2Atest+duxbtest
          endif
!
!  diffusive time step, just take the max of diffus_eta (if existent)
!  and whatever is calculated here
!
          if (lfirst.and.ldt) then
            diffus_eta=max(diffus_eta,etatest*dxyz_2)
          endif
!
!  calculate alpha, begin by calculating uxbtest (if not already done above because of SOCA)
!
!          if ((ldiagnos.or.l1davgfirst).and.lsoca) then
!              call curl(f,iaxtest,btest)
!              call cross_mn(p%uu,btest,uxbtest)
!           endif
!
!  in the following block, we have already swapped the 4-6 entries with 7-9
!
     enddo
!
!  in the following block, we have already swapped the 4-6 entries with 7-9
!
     if (.false.) then           !l1davgfirst) then
       select case (jtest)
       case (1)
         call xysum_mn_name_z(uxbtest(:,1),idiag_alp11z)
         call xysum_mn_name_z(uxbtest(:,2),idiag_alp21z)
         call xysum_mn_name_z(uxbtest(:,3),idiag_alp31z)
       case (2)
         call xysum_mn_name_z(uxbtest(:,1),idiag_alp12z)
         call xysum_mn_name_z(uxbtest(:,2),idiag_alp22z)
         call xysum_mn_name_z(uxbtest(:,3),idiag_alp32z)
       case (3)
         call xysum_mn_name_z(uxbtest(:,1),idiag_alp13z)
         call xysum_mn_name_z(uxbtest(:,2),idiag_alp23z)
         call xysum_mn_name_z(uxbtest(:,3),idiag_alp33z)
       case (4)
         call xysum_mn_name_z(uxbtest(:,1),idiag_eta113z)
         call xysum_mn_name_z(uxbtest(:,2),idiag_eta213z)
         call xysum_mn_name_z(uxbtest(:,3),idiag_eta313z)
       case (5)
         call xysum_mn_name_z(uxbtest(:,1),idiag_eta123z)
         call xysum_mn_name_z(uxbtest(:,2),idiag_eta223z)
         call xysum_mn_name_z(uxbtest(:,3),idiag_eta323z)
       case (6)
         call xysum_mn_name_z(uxbtest(:,1),idiag_eta133z)
         call xysum_mn_name_z(uxbtest(:,2),idiag_eta233z)
         call xysum_mn_name_z(uxbtest(:,3),idiag_eta333z)
       case (7)
         call xysum_mn_name_z(uxbtest(:,1),idiag_eta111z)
         call xysum_mn_name_z(uxbtest(:,2),idiag_eta211z)
         call xysum_mn_name_z(uxbtest(:,3),idiag_eta311z)
       case (8)
         call xysum_mn_name_z(uxbtest(:,1),idiag_eta121z)
         call xysum_mn_name_z(uxbtest(:,2),idiag_eta221z)
         call xysum_mn_name_z(uxbtest(:,3),idiag_eta321z)
       case (9)
         call xysum_mn_name_z(uxbtest(:,1),idiag_eta131z)
         call xysum_mn_name_z(uxbtest(:,2),idiag_eta231z)
         call xysum_mn_name_z(uxbtest(:,3),idiag_eta331z)
       end select
     endif
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
      integer :: jtest,j,nxy=nxgrid*nygrid, nscan
      logical :: headtt_save
      real :: fac
      real, dimension (nx,nz,3) :: uxbtestm1
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

        iaxtest=iaatest+3*(jtest-1)
        iaztest=iaxtest+2
        if (lsoca .and.(.not. need_output)) then
          uxbtestm(:,:,:,jtest)=0.
        else

          do n=n1,n2

            nscan=n-n1+1

            uxbtestm(:,nscan,:,jtest)=0.

            do m=m1,m2

              call calc_pencils_hydro(f,p)
              call curl(f,iaxtest,btest)
              call cross_mn(p%uu,btest,uxbtest)
!
!  without SOCA, the alpha tensor is anisotropic even for the standard
!  Roberts flow. To check that this is because of the averaging that
!  enters, we allow outselves to turn it off with lflucts_with_xyaver=.true.
!  It is off by default.
!
              do j=1,3
                if (lflucts_with_xyaver) then
                  uxbtestm(:,nscan,j,jtest)=spread(sum( &
                    uxbtestm(:,nscan,j,jtest)+fac*uxbtest(:,j),1),1,nx)/nx
                else
                  uxbtestm(:,nscan,j,jtest)= &
                    uxbtestm(:,nscan,j,jtest)+fac*uxbtest(:,j)
                endif
              enddo
              headtt=.false.
            enddo
          enddo
!
!  do communication along y
!
          call mpiallreduce_sum(uxbtestm(:,:,:,jtest),uxbtestm1,(/nx,nz,3/),idir=2)
!         or
!         call mpiallreduce_sum_arr(uxbtestm(1,1,1,jtest),uxbtestm1,nx*nz*3,idir=2)       !avoids copy
          uxbtestm(:,:,:,jtest)=uxbtestm1/nprocy
!
        endif
      enddo
!
!  reset headtt
!
      headtt=headtt_save
!
      if (need_output) call calc_coefficients

    endsubroutine testfield_after_boundary
!***********************************************************************
    subroutine calc_coefficients
!
    Use Diagnostics
    Use Sub, only: fourier_single_mode
!
    integer :: j, count
    real, dimension(2,nx) :: temp_fft_z
    real, dimension(2,2) :: temp_fft
    real, dimension (:,:,:), allocatable :: temp_array
    logical, dimension(27) :: need_temp
    integer, dimension (27) :: twod_address
!
    if (l2davgfirst .and. needed2d_2d) need_temp=twod_need_2d
    if (ldiagnos .and. needed2d_1d) need_temp=need_temp .or. twod_need_1d
!
    count=0
    twod_address=0
    do j=1,27
      if (need_temp(j)) then
        count=count+1
        twod_address(j)=count
      endif
    enddo
!
    allocate(temp_array(nx,nz,count))
!
     select case (itestfield)
    case (1)

      do n=1,nz

      if (need_temp(1)) & !(idiag_alp11xz/=0) &
        temp_array(:,n,twod_address(1))= &
            Minv(:,n,1,1)*uxbtestm(:,n,1,1)+Minv(:,n,1,2)*uxbtestm(:,n,1,2)+Minv(:,n,1,3)*uxbtestm(:,n,1,3)

      if (need_temp(2)) & !(idiag_alp21xz/=0) &
        temp_array(:,n,twod_address(2))= &
            Minv(:,n,1,1)*uxbtestm(:,n,2,1)+Minv(:,n,1,2)*uxbtestm(:,n,2,2)+Minv(:,n,1,3)*uxbtestm(:,n,2,3)

      if (need_temp(3)) & !(idiag_alp31xz/=0) &
        temp_array(:,n,twod_address(3))= &
            Minv(:,n,1,1)*uxbtestm(:,n,3,1)+Minv(:,n,1,2)*uxbtestm(:,n,3,2)+Minv(:,n,1,3)*uxbtestm(:,n,3,3)

      if (need_temp(4)) & !(idiag_alp12xz/=0) &
        temp_array(:,n,twod_address(4))= &
            Minv(:,n,1,1)*uxbtestm(:,n,1,4)+Minv(:,n,1,2)*uxbtestm(:,n,1,5)+Minv(:,n,1,3)*uxbtestm(:,n,1,6)

      if (need_temp(5)) & !(idiag_alp22xz/=0) &
        temp_array(:,n,twod_address(5))= &
            Minv(:,n,1,1)*uxbtestm(:,n,2,4)+Minv(:,n,1,2)*uxbtestm(:,n,2,5)+Minv(:,n,1,3)*uxbtestm(:,n,2,6)

      if (need_temp(6)) & !(idiag_alp32xz/=0) &
        temp_array(:,n,twod_address(6))= &
            Minv(:,n,1,1)*uxbtestm(:,n,3,4)+Minv(:,n,1,2)*uxbtestm(:,n,3,5)+Minv(:,n,1,3)*uxbtestm(:,n,3,6)

      if (need_temp(7)) & !(idiag_alp13xz/=0) &
        temp_array(:,n,twod_address(7))= &
            Minv(:,n,1,1)*uxbtestm(:,n,1,7)+Minv(:,n,1,2)*uxbtestm(:,n,1,8)+Minv(:,n,1,3)*uxbtestm(:,n,1,9)

      if (need_temp(8)) & !(idiag_alp23xz/=0) &
        temp_array(:,n,twod_address(8))= &
            Minv(:,n,1,1)*uxbtestm(:,n,2,7)+Minv(:,n,1,2)*uxbtestm(:,n,2,8)+Minv(:,n,1,3)*uxbtestm(:,n,2,9)

      if (need_temp(9)) & !(idiag_alp33xz/=0) &
        temp_array(:,n,twod_address(9))= &
            Minv(:,n,1,1)*uxbtestm(:,n,3,7)+Minv(:,n,1,2)*uxbtestm(:,n,3,8)+Minv(:,n,1,3)*uxbtestm(:,n,3,9)
!

      if (need_temp(10)) & !(idiag_eta111xz/=0) &
        temp_array(:,n,twod_address(10))= &
            Minv(:,n,2,1)*uxbtestm(:,n,1,1)+Minv(:,n,2,2)*uxbtestm(:,n,1,2)+Minv(:,n,2,3)*uxbtestm(:,n,1,3)

      if (need_temp(11)) & !(idiag_eta112xz/=0) &
        temp_array(:,n,twod_address(11))= &
            Minv(:,n,3,1)*uxbtestm(:,n,1,1)+Minv(:,n,3,2)*uxbtestm(:,n,1,2)+Minv(:,n,3,3)*uxbtestm(:,n,1,3)

      if (need_temp(12)) & !(idiag_eta211xz/=0) &
        temp_array(:,n,twod_address(12))= &
            Minv(:,n,2,1)*uxbtestm(:,n,2,1)+Minv(:,n,2,2)*uxbtestm(:,n,2,2)+Minv(:,n,2,3)*uxbtestm(:,n,2,3)

      if (need_temp(13)) & !(idiag_eta212xz/=0) &
        temp_array(:,n,twod_address(13))= &
            Minv(:,n,3,1)*uxbtestm(:,n,2,1)+Minv(:,n,3,2)*uxbtestm(:,n,2,2)+Minv(:,n,3,3)*uxbtestm(:,n,2,3)

      if (need_temp(14)) & !(idiag_eta311xz/=0) &
        temp_array(:,n,twod_address(14))= &
            Minv(:,n,2,1)*uxbtestm(:,n,3,1)+Minv(:,n,2,2)*uxbtestm(:,n,3,2)+Minv(:,n,2,3)*uxbtestm(:,n,3,3)

      if (need_temp(15)) & !(idiag_eta312xz/=0) &
        temp_array(:,n,twod_address(15))= &
            Minv(:,n,3,1)*uxbtestm(:,n,3,1)+Minv(:,n,3,2)*uxbtestm(:,n,3,2)+Minv(:,n,3,3)*uxbtestm(:,n,3,3)

!

      if (need_temp(16)) & !(idiag_eta121xz/=0) &
        temp_array(:,n,twod_address(16))= &
            Minv(:,n,2,1)*uxbtestm(:,n,1,4)+Minv(:,n,2,2)*uxbtestm(:,n,1,5)+Minv(:,n,2,3)*uxbtestm(:,n,1,6)

      if (need_temp(17)) & !(idiag_eta122xz/=0) &
        temp_array(:,n,twod_address(17))= &
            Minv(:,n,3,1)*uxbtestm(:,n,1,4)+Minv(:,n,3,2)*uxbtestm(:,n,1,5)+Minv(:,n,3,3)*uxbtestm(:,n,1,6)

      if (need_temp(18)) & !(idiag_eta221xz/=0) &
        temp_array(:,n,twod_address(18))= &
            Minv(:,n,2,1)*uxbtestm(:,n,2,4)+Minv(:,n,2,2)*uxbtestm(:,n,2,5)+Minv(:,n,2,3)*uxbtestm(:,n,2,6)

      if (need_temp(19)) & !(idiag_eta222xz/=0) &
        temp_array(:,n,twod_address(19))= &
            Minv(:,n,3,1)*uxbtestm(:,n,2,4)+Minv(:,n,3,2)*uxbtestm(:,n,2,5)+Minv(:,n,3,3)*uxbtestm(:,n,2,6)

      if (need_temp(20)) & !(idiag_eta321xz/=0) &
        temp_array(:,n,twod_address(20))= &
            Minv(:,n,2,1)*uxbtestm(:,n,3,4)+Minv(:,n,2,2)*uxbtestm(:,n,3,5)+Minv(:,n,2,3)*uxbtestm(:,n,3,6)

      if (need_temp(21)) & !(idiag_eta322xz/=0) &
        temp_array(:,n,twod_address(21))= &
            Minv(:,n,3,1)*uxbtestm(:,n,3,4)+Minv(:,n,3,2)*uxbtestm(:,n,3,5)+Minv(:,n,3,3)*uxbtestm(:,n,3,6)


      if (need_temp(22)) & !(idiag_eta131xz/=0) &
        temp_array(:,n,twod_address(22))= &
            Minv(:,n,2,1)*uxbtestm(:,n,1,7)+Minv(:,n,2,2)*uxbtestm(:,n,1,8)+Minv(:,n,2,3)*uxbtestm(:,n,1,9)

      if (need_temp(23)) & !(idiag_eta132xz/=0) &
        temp_array(:,n,twod_address(23))= &
            Minv(:,n,3,1)*uxbtestm(:,n,1,7)+Minv(:,n,3,2)*uxbtestm(:,n,1,8)+Minv(:,n,3,3)*uxbtestm(:,n,1,9)

      if (need_temp(24)) & !(idiag_eta231xz/=0) &
        temp_array(:,n,twod_address(24))= &
            Minv(:,n,2,1)*uxbtestm(:,n,2,7)+Minv(:,n,2,2)*uxbtestm(:,n,2,8)+Minv(:,n,2,3)*uxbtestm(:,n,2,9)

      if (need_temp(25)) & !(idiag_eta232xz/=0) &
        temp_array(:,n,twod_address(25))= &
            Minv(:,n,3,1)*uxbtestm(:,n,2,7)+Minv(:,n,3,2)*uxbtestm(:,n,2,8)+Minv(:,n,3,3)*uxbtestm(:,n,2,9)

      if (need_temp(26)) & !(idiag_eta331xz/=0) &
        temp_array(:,n,twod_address(26))= &
            Minv(:,n,2,1)*uxbtestm(:,n,3,7)+Minv(:,n,2,2)*uxbtestm(:,n,3,8)+Minv(:,n,2,3)*uxbtestm(:,n,3,9)

      if (need_temp(27)) & !(idiag_eta332xz/=0) &
        temp_array(:,n,twod_address(27))= &
            Minv(:,n,3,1)*uxbtestm(:,n,3,7)+Minv(:,n,3,2)*uxbtestm(:,n,3,8)+Minv(:,n,3,3)*uxbtestm(:,n,3,9)

    enddo
    case default

    end select


    if (l2davgfirst .and. needed2d_2d) then
!
    do n=n1,n2 ! ny multiplied in because we are not in the mn loop
!
      if (idiag_alp11xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(1)),idiag_alp11xz)

      if (idiag_alp21xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(2)),idiag_alp21xz)

      if (idiag_alp31xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(3)),idiag_alp31xz)

      if (idiag_alp12xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(4)),idiag_alp12xz)

      if (idiag_alp22xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(5)),idiag_alp22xz)

      if (idiag_alp32xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(6)),idiag_alp32xz)

      if (idiag_alp13xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(7)),idiag_alp13xz)

      if (idiag_alp23xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(8)),idiag_alp23xz)

      if (idiag_alp33xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(9)),idiag_alp33xz)

!

      if (idiag_eta111xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(10)),idiag_eta111xz)

      if (idiag_eta112xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(11)),idiag_eta112xz)

      if (idiag_eta211xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(12)),idiag_eta211xz)

      if (idiag_eta212xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(13)),idiag_eta212xz)

      if (idiag_eta311xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(14)),idiag_eta311xz)

      if (idiag_eta312xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(15)),idiag_eta312xz)



      if (idiag_eta121xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(16)),idiag_eta121xz)

      if (idiag_eta122xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(17)),idiag_eta122xz)

      if (idiag_eta221xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(18)),idiag_eta221xz)

      if (idiag_eta222xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(19)),idiag_eta222xz)

      if (idiag_eta321xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(20)),idiag_eta321xz)

      if (idiag_eta322xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(21)),idiag_eta322xz)


      if (idiag_eta131xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(22)),idiag_eta131xz)

      if (idiag_eta132xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(23)),idiag_eta132xz)

      if (idiag_eta231xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(24)),idiag_eta231xz)

      if (idiag_eta232xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(25)),idiag_eta232xz)

      if (idiag_eta331xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(26)),idiag_eta331xz)

      if (idiag_eta332xz/=0) &
        call ysum_mn_name_xz(ny*temp_array(:,n-n1+1,twod_address(27)),idiag_eta332xz)
    enddo
    endif
!
    if (ldiagnos .and. needed2d_1d) then

      do n=n1,n2

        if (idiag_alp11/=0) call sum_mn_name(ny*temp_array(:,n-n1+1,twod_address(1)), idiag_alp11)
        if (idiag_alp22/=0) call sum_mn_name(ny*temp_array(:,n-n1+1,twod_address(5)), idiag_alp22)
        if (idiag_alp33/=0) call sum_mn_name(ny*temp_array(:,n-n1+1,twod_address(9)), idiag_alp33)


        if (idiag_eta122/=0) call sum_mn_name(ny*temp_array(:,n-n1+1, twod_address(17)), idiag_eta122)

      enddo

      if (idiag_alp11cc/=0 .or. idiag_alp11cs/=0 .or. idiag_alp11sc/=0 .or. idiag_alp11ss/=0) then
        call fourier_single_mode(temp_array(:,:,twod_address(1)), &
            (/nx,nz/), 1., 3, temp_fft_z, l2nd=.true.)
        if (lroot) then
          call fourier_single_mode(temp_fft_z, (/2,nx/), 1., 1, temp_fft, l2nd=.true.)
          if (idiag_alp11cc/=0) call save_name(temp_fft(1,1), idiag_alp11cc)
          if (idiag_alp11cs/=0) call save_name(temp_fft(1,2), idiag_alp11cs)
          if (idiag_alp11sc/=0) call save_name(temp_fft(2,1), idiag_alp11sc)
          if (idiag_alp11ss/=0) call save_name(temp_fft(2,2), idiag_alp11ss)
        endif
      endif

      if (idiag_eta122cc/=0 .or. idiag_eta122cs/=0 .or. idiag_eta122sc/=0 .or. idiag_eta122ss/=0) then
        call fourier_single_mode(temp_array(:,:,twod_address(17)), &
            (/nx,nz/), 1., 3, temp_fft_z, l2nd=.true.)
        if (lroot) then
          call fourier_single_mode(temp_fft_z, (/2,nx/), 1., 1, temp_fft, l2nd=.true.)
          if (idiag_eta122cc/=0) call save_name(temp_fft(1,1), idiag_eta122cc)
          if (idiag_eta122cs/=0) call save_name(temp_fft(1,2), idiag_eta122cs)
          if (idiag_eta122sc/=0) call save_name(temp_fft(2,1), idiag_eta122sc)
          if (idiag_eta122ss/=0) call save_name(temp_fft(2,2), idiag_eta122ss)
        endif
      endif

    endif
!
    deallocate(temp_array)
!
    endsubroutine calc_coefficients
!***********************************************************************
    subroutine rescaling_testfield(f)
!
!  Rescale testfield by factor rescale_aatest(jtest),
!  which could be different for different testfields
!
!  18-may-08/axel: rewrite from rescaling as used in magnetic
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=fnlen) :: file
      logical :: ltestfield_out
      logical, save :: lfirst_call=.true.
      integer :: j,jtest
!
      intent(inout) :: f
!
! reinitialize aatest periodically if requested
!
      if (linit_aatest) then
        file=trim(datadir)//'/tinit_aatest.dat'
        if (lfirst_call) then
          call read_snaptime(trim(file),taainit,naainit,daainit,t)
          if (taainit==0 .or. taainit < t-daainit) then
            taainit=t+daainit
          endif
          lfirst_call=.false.
        endif
!
!  Do only one xy plane at a time (for cache efficiency)
!
        if (t >= taainit) then
          do jtest=1,njtest
            iaxtest=iaatest+3*(jtest-1)
            iaztest=iaxtest+2
            do j=iaxtest,iaztest
              do n=n1,n2
                f(l1:l2,m1:m2,n,j)=rescale_aatest(jtest)*f(l1:l2,m1:m2,n,j)
              enddo
            enddo
          enddo
          call update_snaptime(file,taainit,naainit,daainit,t,ltestfield_out)
        endif
      endif
!
    endsubroutine rescaling_testfield
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
!  set bbtest for each of the 9 cases
!
      select case (jtest)

      case (1); bbtest(:,1)=cx*cz(n-n1+1); bbtest(:,2)=0.; bbtest(:,3)=0.
      case (2); bbtest(:,1)=sx*cz(n-n1+1); bbtest(:,2)=0.; bbtest(:,3)=0.
      case (3); bbtest(:,1)=cx*sz(n-n1+1); bbtest(:,2)=0.; bbtest(:,3)=0.

      case (4); bbtest(:,1)=0.; bbtest(:,2)=cx*cz(n-n1+1); bbtest(:,3)=0.
      case (5); bbtest(:,1)=0.; bbtest(:,2)=sx*cz(n-n1+1); bbtest(:,3)=0.
      case (6); bbtest(:,1)=0.; bbtest(:,2)=cx*sz(n-n1+1); bbtest(:,3)=0.

      case (7); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cx*cz(n-n1+1)
      case (8); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sx*cz(n-n1+1)
      case (9); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cx*sz(n-n1+1)

      case default; bbtest(:,:)=0.

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
      real, dimension (nx) :: cx,sx,cz,sz,xz
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: bbtest
!
!  xx and zz for calculating diffusive part of emf
!
      cx=cos(x(l1:l2))
      sx=sin(x(l1:l2))
      cz=cos(z(n))
      sz=sin(z(n))
      xz=cx*cz
!
!  set bbtest for each of the 9 cases
!
      select case (jtest)
      case (1); bbtest(:,1)=xz; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (2); bbtest(:,1)=0.; bbtest(:,2)=xz; bbtest(:,3)=0.
      case (3); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=xz
      case (4); bbtest(:,1)=sz; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (5); bbtest(:,1)=0.; bbtest(:,2)=sz; bbtest(:,3)=0.
      case (6); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sz
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
      real, dimension (nx) :: xx,zz
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: bbtest
!
!  xx and zz for calculating diffusive part of emf
!
      xx=x(l1:l2)
      zz=z(n)
!
!  set bbtest for each of the 9 cases
!
      select case (jtest)
      case (1); bbtest(:,1)=1.; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (2); bbtest(:,1)=0.; bbtest(:,2)=1.; bbtest(:,3)=0.
      case (3); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=1.
      case (4); bbtest(:,1)=zz; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (5); bbtest(:,1)=0.; bbtest(:,2)=zz; bbtest(:,3)=0.
      case (6); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=zz
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
      real, dimension (nx) :: cx,sx,cz,sz
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: bbtest
!
!  xx and zz for calculating diffusive part of emf
!
      cx=cos(x(l1:l2))
      sx=sin(x(l1:l2))
      cz=cos(z(n))
      sz=sin(z(n))
!
!  set bbtest for each of the 9 cases
!
      select case (jtest)
      case (1); bbtest(:,1)=1.; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (2); bbtest(:,1)=0.; bbtest(:,2)=1.; bbtest(:,3)=0.
      case (3); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=1.
      case (4); bbtest(:,1)=sz; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (5); bbtest(:,1)=0.; bbtest(:,2)=sz; bbtest(:,3)=0.
      case (6); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cz
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
!
      use Cdata
      use Diagnostics
!
      integer :: iname,inamez,inamexz
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of RELOAD
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_alp11=0; idiag_alp21=0; idiag_alp31=0
        idiag_alp12=0; idiag_alp22=0; idiag_alp32=0
        idiag_alp13=0; idiag_alp23=0; idiag_alp33=0
        idiag_eta111=0; idiag_eta211=0; idiag_eta311=0
        idiag_eta112=0; idiag_eta212=0; idiag_eta312=0
        idiag_eta121=0; idiag_eta221=0; idiag_eta321=0
        idiag_eta122=0; idiag_eta222=0; idiag_eta322=0
        idiag_eta131=0; idiag_eta231=0; idiag_eta331=0
        idiag_eta132=0; idiag_eta232=0; idiag_eta332=0
        idiag_alp11cc=0;idiag_alp11cs=0;idiag_alp11sc=0;idiag_alp11ss=0
        idiag_eta122cc=0;idiag_eta122cs=0;idiag_eta122sc=0;idiag_eta122ss=0
        idiag_alp11z=0; idiag_alp21z=0; idiag_alp31z=0
        idiag_alp12z=0; idiag_alp22z=0; idiag_alp32z=0
        idiag_alp13z=0; idiag_alp23z=0; idiag_alp33z=0
        idiag_eta111z=0; idiag_eta211z=0; idiag_eta311z=0
        idiag_eta121z=0; idiag_eta221z=0; idiag_eta321z=0
        idiag_eta131z=0; idiag_eta231z=0; idiag_eta331z=0
        idiag_eta113z=0; idiag_eta213z=0; idiag_eta313z=0
        idiag_eta123z=0; idiag_eta223z=0; idiag_eta323z=0
        idiag_eta133z=0; idiag_eta233z=0; idiag_eta333z=0
        idiag_alp11xz=0; idiag_alp21xz=0; idiag_alp31xz=0
        idiag_alp12xz=0; idiag_alp22xz=0; idiag_alp32xz=0
        idiag_alp13xz=0; idiag_alp23xz=0; idiag_alp33xz=0
        idiag_eta111xz=0; idiag_eta211xz=0; idiag_eta311xz=0
        idiag_eta112xz=0; idiag_eta212xz=0; idiag_eta312xz=0
        idiag_eta121xz=0; idiag_eta221xz=0; idiag_eta321xz=0
        idiag_eta122xz=0; idiag_eta222xz=0; idiag_eta322xz=0
        idiag_eta131xz=0; idiag_eta231xz=0; idiag_eta331xz=0
        idiag_eta132xz=0; idiag_eta232xz=0; idiag_eta332xz=0
        idiag_eta113xz=0; idiag_eta213xz=0; idiag_eta313xz=0
        idiag_eta123xz=0; idiag_eta223xz=0; idiag_eta323xz=0
        idiag_eta133xz=0; idiag_eta233xz=0; idiag_eta333xz=0
        idiag_alp11exz=0; idiag_alp21exz=0; idiag_alp31exz=0
        idiag_alp12exz=0; idiag_alp22exz=0; idiag_alp32exz=0
        idiag_alp13exz=0; idiag_alp23exz=0; idiag_alp33exz=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'alp11',idiag_alp11)
        call parse_name(iname,cname(iname),cform(iname),'alp21',idiag_alp21)
        call parse_name(iname,cname(iname),cform(iname),'alp31',idiag_alp31)
        call parse_name(iname,cname(iname),cform(iname),'alp12',idiag_alp12)
        call parse_name(iname,cname(iname),cform(iname),'alp22',idiag_alp22)
        call parse_name(iname,cname(iname),cform(iname),'alp32',idiag_alp32)
        call parse_name(iname,cname(iname),cform(iname),'alp13',idiag_alp13)
        call parse_name(iname,cname(iname),cform(iname),'alp23',idiag_alp23)
        call parse_name(iname,cname(iname),cform(iname),'alp33',idiag_alp33)
        call parse_name(iname,cname(iname),cform(iname),'eta111',idiag_eta111)
        call parse_name(iname,cname(iname),cform(iname),'eta112',idiag_eta112)
        call parse_name(iname,cname(iname),cform(iname),'eta211',idiag_eta211)
        call parse_name(iname,cname(iname),cform(iname),'eta212',idiag_eta212)
        call parse_name(iname,cname(iname),cform(iname),'eta311',idiag_eta311)
        call parse_name(iname,cname(iname),cform(iname),'eta312',idiag_eta312)
        call parse_name(iname,cname(iname),cform(iname),'eta121',idiag_eta121)
        call parse_name(iname,cname(iname),cform(iname),'eta122',idiag_eta122)
        call parse_name(iname,cname(iname),cform(iname),'eta221',idiag_eta221)
        call parse_name(iname,cname(iname),cform(iname),'eta222',idiag_eta222)
        call parse_name(iname,cname(iname),cform(iname),'eta321',idiag_eta321)
        call parse_name(iname,cname(iname),cform(iname),'eta322',idiag_eta322)
        call parse_name(iname,cname(iname),cform(iname),'eta131',idiag_eta131)
        call parse_name(iname,cname(iname),cform(iname),'eta132',idiag_eta132)
        call parse_name(iname,cname(iname),cform(iname),'eta231',idiag_eta231)
        call parse_name(iname,cname(iname),cform(iname),'eta232',idiag_eta232)
        call parse_name(iname,cname(iname),cform(iname),'eta331',idiag_eta331)
        call parse_name(iname,cname(iname),cform(iname),'eta332',idiag_eta332)
        call parse_name(iname,cname(iname),cform(iname),'alp11cc',idiag_alp11cc)
        call parse_name(iname,cname(iname),cform(iname),'alp11cs',idiag_alp11cs)
        call parse_name(iname,cname(iname),cform(iname),'alp11sc',idiag_alp11sc)
        call parse_name(iname,cname(iname),cform(iname),'alp11ss',idiag_alp11ss)
        call parse_name(iname,cname(iname),cform(iname),'eta122cc',idiag_eta122cc)
        call parse_name(iname,cname(iname),cform(iname),'eta122cs',idiag_eta122cs)
        call parse_name(iname,cname(iname),cform(iname),'eta122sc',idiag_eta122sc)
        call parse_name(iname,cname(iname),cform(iname),'eta122ss',idiag_eta122ss)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp11z',idiag_alp11z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp21z',idiag_alp21z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp31z',idiag_alp31z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp12z',idiag_alp12z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp22z',idiag_alp22z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp32z',idiag_alp32z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp13z',idiag_alp13z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp23z',idiag_alp23z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp33z',idiag_alp33z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta111z',idiag_eta111z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta211z',idiag_eta211z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta311z',idiag_eta311z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta121z',idiag_eta121z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta221z',idiag_eta221z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta321z',idiag_eta321z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta131z',idiag_eta131z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta231z',idiag_eta231z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta331z',idiag_eta331z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta113z',idiag_eta113z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta213z',idiag_eta213z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta313z',idiag_eta313z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta123z',idiag_eta123z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta223z',idiag_eta223z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta323z',idiag_eta323z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta133z',idiag_eta133z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta233z',idiag_eta233z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta333z',idiag_eta333z)
      enddo
!
!  check for those quantities for which we want y-averages
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp11xz',idiag_alp11xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp21xz',idiag_alp21xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp31xz',idiag_alp31xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp12xz',idiag_alp12xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp22xz',idiag_alp22xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp32xz',idiag_alp32xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp13xz',idiag_alp13xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp23xz',idiag_alp23xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp33xz',idiag_alp33xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta111xz',idiag_eta111xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta112xz',idiag_eta112xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta211xz',idiag_eta211xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta212xz',idiag_eta212xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta311xz',idiag_eta311xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta312xz',idiag_eta312xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta121xz',idiag_eta121xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta122xz',idiag_eta122xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta221xz',idiag_eta221xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta222xz',idiag_eta222xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta321xz',idiag_eta321xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta322xz',idiag_eta322xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta131xz',idiag_eta131xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta132xz',idiag_eta132xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta231xz',idiag_eta231xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta232xz',idiag_eta232xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta331xz',idiag_eta331xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta332xz',idiag_eta332xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta113xz',idiag_eta113xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta213xz',idiag_eta213xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta313xz',idiag_eta313xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta123xz',idiag_eta123xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta223xz',idiag_eta223xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta323xz',idiag_eta323xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta133xz',idiag_eta133xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta233xz',idiag_eta233xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta333xz',idiag_eta333xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp11exz',idiag_alp11exz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp21exz',idiag_alp21exz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp31exz',idiag_alp31exz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp12exz',idiag_alp12exz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp22exz',idiag_alp22exz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp32exz',idiag_alp32exz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp13exz',idiag_alp13exz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp23exz',idiag_alp23exz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp33exz',idiag_alp33exz)
      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!
      if (lwr) then
        write(3,*) 'idiag_alp11=',idiag_alp11
        write(3,*) 'idiag_alp21=',idiag_alp21
        write(3,*) 'idiag_alp31=',idiag_alp31
        write(3,*) 'idiag_alp12=',idiag_alp12
        write(3,*) 'idiag_alp22=',idiag_alp22
        write(3,*) 'idiag_alp32=',idiag_alp32
        write(3,*) 'idiag_alp13=',idiag_alp13
        write(3,*) 'idiag_alp23=',idiag_alp23
        write(3,*) 'idiag_alp33=',idiag_alp33
        write(3,*) 'idiag_eta111=',idiag_eta111
        write(3,*) 'idiag_eta112=',idiag_eta112
        write(3,*) 'idiag_eta211=',idiag_eta211
        write(3,*) 'idiag_eta212=',idiag_eta212
        write(3,*) 'idiag_eta311=',idiag_eta311
        write(3,*) 'idiag_eta312=',idiag_eta312
        write(3,*) 'idiag_eta121=',idiag_eta121
        write(3,*) 'idiag_eta122=',idiag_eta122
        write(3,*) 'idiag_eta221=',idiag_eta221
        write(3,*) 'idiag_eta222=',idiag_eta222
        write(3,*) 'idiag_eta321=',idiag_eta321
        write(3,*) 'idiag_eta322=',idiag_eta322
        write(3,*) 'idiag_eta131=',idiag_eta131
        write(3,*) 'idiag_eta132=',idiag_eta132
        write(3,*) 'idiag_eta231=',idiag_eta231
        write(3,*) 'idiag_eta232=',idiag_eta232
        write(3,*) 'idiag_eta331=',idiag_eta331
        write(3,*) 'idiag_eta332=',idiag_eta332
        write(3,*) 'idiag_alp11cc=',idiag_alp11cc
        write(3,*) 'idiag_alp11cs=',idiag_alp11cs
        write(3,*) 'idiag_alp11sc=',idiag_alp11sc
        write(3,*) 'idiag_alp11ss=',idiag_alp11ss
        write(3,*) 'idiag_eta122cc=',idiag_eta122cc
        write(3,*) 'idiag_eta122cs=',idiag_eta122cs
        write(3,*) 'idiag_eta122sc=',idiag_eta122sc
        write(3,*) 'idiag_eta122ss=',idiag_eta122ss
        write(3,*) 'idiag_alp11z=',idiag_alp11z
        write(3,*) 'idiag_alp21z=',idiag_alp21z
        write(3,*) 'idiag_alp31z=',idiag_alp31z
        write(3,*) 'idiag_alp12z=',idiag_alp12z
        write(3,*) 'idiag_alp22z=',idiag_alp22z
        write(3,*) 'idiag_alp32z=',idiag_alp32z
        write(3,*) 'idiag_alp13z=',idiag_alp13z
        write(3,*) 'idiag_alp23z=',idiag_alp23z
        write(3,*) 'idiag_alp33z=',idiag_alp33z
        write(3,*) 'idiag_eta111z=',idiag_eta111z
        write(3,*) 'idiag_eta211z=',idiag_eta211z
        write(3,*) 'idiag_eta311z=',idiag_eta311z
        write(3,*) 'idiag_eta121z=',idiag_eta121z
        write(3,*) 'idiag_eta221z=',idiag_eta221z
        write(3,*) 'idiag_eta321z=',idiag_eta321z
        write(3,*) 'idiag_eta131z=',idiag_eta131z
        write(3,*) 'idiag_eta231z=',idiag_eta231z
        write(3,*) 'idiag_eta331z=',idiag_eta331z
        write(3,*) 'idiag_eta113z=',idiag_eta113z
        write(3,*) 'idiag_eta213z=',idiag_eta213z
        write(3,*) 'idiag_eta313z=',idiag_eta313z
        write(3,*) 'idiag_eta123z=',idiag_eta123z
        write(3,*) 'idiag_eta223z=',idiag_eta223z
        write(3,*) 'idiag_eta323z=',idiag_eta323z
        write(3,*) 'idiag_eta133z=',idiag_eta133z
        write(3,*) 'idiag_eta233z=',idiag_eta233z
        write(3,*) 'idiag_eta333z=',idiag_eta333z
        write(3,*) 'idiag_alp11xz=',idiag_alp11xz
        write(3,*) 'idiag_alp21xz=',idiag_alp21xz
        write(3,*) 'idiag_alp31xz=',idiag_alp31xz
        write(3,*) 'idiag_alp12xz=',idiag_alp12xz
        write(3,*) 'idiag_alp22xz=',idiag_alp22xz
        write(3,*) 'idiag_alp32xz=',idiag_alp32xz
        write(3,*) 'idiag_alp13xz=',idiag_alp13xz
        write(3,*) 'idiag_alp23xz=',idiag_alp23xz
        write(3,*) 'idiag_alp33xz=',idiag_alp33xz
        write(3,*) 'idiag_eta111xz=',idiag_eta111xz
        write(3,*) 'idiag_eta211xz=',idiag_eta211xz
        write(3,*) 'idiag_eta311xz=',idiag_eta311xz
        write(3,*) 'idiag_eta112xz=',idiag_eta112xz
        write(3,*) 'idiag_eta212xz=',idiag_eta212xz
        write(3,*) 'idiag_eta312xz=',idiag_eta312xz
        write(3,*) 'idiag_eta121xz=',idiag_eta121xz
        write(3,*) 'idiag_eta221xz=',idiag_eta221xz
        write(3,*) 'idiag_eta321xz=',idiag_eta321xz
        write(3,*) 'idiag_eta122xz=',idiag_eta122xz
        write(3,*) 'idiag_eta222xz=',idiag_eta222xz
        write(3,*) 'idiag_eta322xz=',idiag_eta322xz
        write(3,*) 'idiag_eta131xz=',idiag_eta131xz
        write(3,*) 'idiag_eta231xz=',idiag_eta231xz
        write(3,*) 'idiag_eta331xz=',idiag_eta331xz
        write(3,*) 'idiag_eta132xz=',idiag_eta132xz
        write(3,*) 'idiag_eta232xz=',idiag_eta232xz
        write(3,*) 'idiag_eta332xz=',idiag_eta332xz
        write(3,*) 'idiag_eta113xz=',idiag_eta113xz
        write(3,*) 'idiag_eta213xz=',idiag_eta213xz
        write(3,*) 'idiag_eta313xz=',idiag_eta313xz
        write(3,*) 'idiag_eta123xz=',idiag_eta123xz
        write(3,*) 'idiag_eta223xz=',idiag_eta223xz
        write(3,*) 'idiag_eta323xz=',idiag_eta323xz
        write(3,*) 'idiag_eta133xz=',idiag_eta133xz
        write(3,*) 'idiag_eta233xz=',idiag_eta233xz
        write(3,*) 'idiag_eta333xz=',idiag_eta333xz
        write(3,*) 'idiag_alp11exz=',idiag_alp11exz
        write(3,*) 'idiag_alp21exz=',idiag_alp21exz
        write(3,*) 'idiag_alp31exz=',idiag_alp31exz
        write(3,*) 'idiag_alp12exz=',idiag_alp12exz
        write(3,*) 'idiag_alp22exz=',idiag_alp22exz
        write(3,*) 'idiag_alp32exz=',idiag_alp32exz
        write(3,*) 'idiag_alp13exz=',idiag_alp13exz
        write(3,*) 'idiag_alp23exz=',idiag_alp23exz
        write(3,*) 'idiag_alp33exz=',idiag_alp33exz
        write(3,*) 'iaatest=',iaatest
        write(3,*) 'nnamez=',nnamez
        write(3,*) 'nnamexy=',nnamexy
        write(3,*) 'nnamexz=',nnamexz
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
    twod_need_2d=.false.
!
    if (idiag_alp11xz/=0) &
        twod_need_2d(1)=.true.
    if (idiag_alp21xz/=0) &
        twod_need_2d(2)=.true.
    if (idiag_alp31xz/=0) &
        twod_need_2d(3)=.true.
    if (idiag_alp12xz/=0) &
        twod_need_2d(4)=.true.
    if (idiag_alp22xz/=0) &
        twod_need_2d(5)=.true.
    if (idiag_alp32xz/=0) &
        twod_need_2d(6)=.true.
    if (idiag_alp13xz/=0) &
        twod_need_2d(7)=.true.
    if (idiag_alp23xz/=0) &
        twod_need_2d(8)=.true.
    if (idiag_alp33xz/=0) &
        twod_need_2d(9)=.true.
!
    if (idiag_eta111xz/=0) &
        twod_need_2d(10)=.true.
    if (idiag_eta112xz/=0) &
        twod_need_2d(11)=.true.
    if (idiag_eta211xz/=0) &
        twod_need_2d(12)=.true.
    if (idiag_eta212xz/=0) &
        twod_need_2d(13)=.true.
    if (idiag_eta311xz/=0) &
        twod_need_2d(14)=.true.
    if (idiag_eta312xz/=0) &
        twod_need_2d(15)=.true.
!
    if (idiag_eta121xz/=0) &
        twod_need_2d(16)=.true.
    if (idiag_eta122xz/=0) &
        twod_need_2d(17)=.true.
    if (idiag_eta221xz/=0) &
        twod_need_2d(18)=.true.
    if (idiag_eta222xz/=0) &
        twod_need_2d(19)=.true.
    if (idiag_eta321xz/=0) &
        twod_need_2d(20)=.true.
    if (idiag_eta322xz/=0) &
        twod_need_2d(21)=.true.
!
    if (idiag_eta131xz/=0) &
        twod_need_2d(22)=.true.
    if (idiag_eta132xz/=0) &
        twod_need_2d(23)=.true.
    if (idiag_eta231xz/=0) &
        twod_need_2d(24)=.true.
    if (idiag_eta232xz/=0) &
        twod_need_2d(25)=.true.
    if (idiag_eta331xz/=0) &
        twod_need_2d(26)=.true.
    if (idiag_eta332xz/=0) &
        twod_need_2d(27)=.true.
!
    needed2d_2d = any(twod_need_2d)
!
!  2d dependencies of 0 or 1-d averages
!
    twod_need_1d=.false.
!
    if (idiag_alp11/=0) &
        twod_need_1d(1)=.true.
    if (idiag_alp21/=0) &
        twod_need_1d(2)=.true.
    if (idiag_alp31/=0) &
        twod_need_1d(3)=.true.
    if (idiag_alp12/=0) &
        twod_need_1d(4)=.true.
    if (idiag_alp22/=0) &
        twod_need_1d(5)=.true.
    if (idiag_alp32/=0) &
        twod_need_1d(6)=.true.
    if (idiag_alp13/=0) &
        twod_need_1d(7)=.true.
    if (idiag_alp23/=0) &
        twod_need_1d(8)=.true.
    if (idiag_alp33/=0) &
        twod_need_1d(9)=.true.
!
    if (idiag_eta111/=0) &
        twod_need_1d(10)=.true.
    if (idiag_eta112/=0) &
        twod_need_1d(11)=.true.
    if (idiag_eta211/=0) &
        twod_need_1d(12)=.true.
    if (idiag_eta212/=0) &
        twod_need_1d(13)=.true.
    if (idiag_eta311/=0) &
        twod_need_1d(14)=.true.
    if (idiag_eta312/=0) &
        twod_need_1d(15)=.true.
!
    if (idiag_eta121/=0) &
        twod_need_1d(16)=.true.
    if (idiag_eta122/=0) &
        twod_need_1d(17)=.true.
    if (idiag_eta221/=0) &
        twod_need_1d(18)=.true.
    if (idiag_eta222/=0) &
        twod_need_1d(19)=.true.
    if (idiag_eta321/=0) &
        twod_need_1d(20)=.true.
    if (idiag_eta322/=0) &
        twod_need_1d(21)=.true.
!
    if (idiag_eta131/=0) &
        twod_need_1d(22)=.true.
    if (idiag_eta132/=0) &
        twod_need_1d(23)=.true.
    if (idiag_eta231/=0) &
        twod_need_1d(24)=.true.
    if (idiag_eta232/=0) &
        twod_need_1d(25)=.true.
    if (idiag_eta331/=0) &
        twod_need_1d(26)=.true.
    if (idiag_eta332/=0) &
        twod_need_1d(27)=.true.
!
    needed2d_1d=any(twod_need_1d)
!
  endsubroutine diagnos_interdep
!***********************************************************************
endmodule Testfield
