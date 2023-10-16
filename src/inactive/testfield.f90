! $Id$

!  This modules deals with all aspects of testfield fields; if no
!  testfield fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testfield relevant subroutines listed in here.
!
!  NOTE: since the fall of 2007 we have been using the routine
!  testfield_z.f90, not this one! For more information, please
!  contact Axel Brandenburg <brandenb@nordita.org>
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ltestfield = .true.
! CPARAM logical, parameter :: ltestfield_z  = .true.
! CPARAM logical, parameter :: ltestfield_xy = .false.
! CPARAM logical, parameter :: ltestfield_xz = .false.
!
! MVAR CONTRIBUTION 36
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Testfield

  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  implicit none

  include '../testfield.h'

  character (len=labellen) :: initaatest='zero'

  ! input parameters
  real, dimension(3) :: B_ext=(/0.,0.,0./)
  real, dimension (nx,3) :: bbb
  real :: amplaa=0., kx_aa=1.,ky_aa=1.,kz_aa=1.
  real :: taainit=0.,daainit=0.
  logical :: reinitialize_aatest=.false.
  logical :: xextent=.true.,zextent=.true.,lsoca=.true.,lset_bbtest2=.false.
  logical :: linit_aatest=.false.
  integer :: itestfield=1
  real :: ktestfield=1.
  integer :: naainit

  namelist /testfield_init_pars/ &
       B_ext,xextent,zextent,initaatest

  ! run parameters
  real :: etatest=0.
  real, dimension(njtest) :: rescale_aatest=0.
  namelist /testfield_run_pars/ &
       B_ext,reinitialize_aatest,xextent,zextent,lsoca, &
       lset_bbtest2,etatest,itestfield,ktestfield,daainit, &
       linit_aatest, &
       rescale_aatest

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_alp11=0,idiag_alp21=0,idiag_alp31=0
  integer :: idiag_alp12=0,idiag_alp22=0,idiag_alp32=0
  integer :: idiag_alp13=0,idiag_alp23=0,idiag_alp33=0
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
  integer :: idiag_eta121xz=0,idiag_eta221xz=0,idiag_eta321xz=0
  integer :: idiag_eta131xz=0,idiag_eta231xz=0,idiag_eta331xz=0
  integer :: idiag_eta113xz=0,idiag_eta213xz=0,idiag_eta313xz=0
  integer :: idiag_eta123xz=0,idiag_eta223xz=0,idiag_eta323xz=0
  integer :: idiag_eta133xz=0,idiag_eta233xz=0,idiag_eta333xz=0
  integer :: idiag_alp11exz=0,idiag_alp21exz=0,idiag_alp31exz=0
  integer :: idiag_alp12exz=0,idiag_alp22exz=0,idiag_alp32exz=0
  integer :: idiag_alp13exz=0,idiag_alp23exz=0,idiag_alp33exz=0

  real, dimension (mz,3,njtest) :: uxbtestm

  contains
!
!***********************************************************************
    subroutine register_testfield
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaatest, etc; increase nvar accordingly
!
!   3-jun-05/axel: adapted from register_magnetic
!
      use FArrayManager
      use Mpicomm
      use Sub
!
!  Register test field.
!
      ntestfield = 3*njtest
      call farray_register_pde('aatest',iaatest,array=ntestfield)
      call farray_index_append('ntestfield',ntestfield)
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Writing files for use with IDL.
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
    subroutine initialize_testfield(f)
!
!  Perform any post-parameter-read initialization
!
!   2-jun-05/axel: adapted from magnetic
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  set to zero and then rescale the testfield
!  (in future, could call something like init_aa_simple)
!
      if (reinitialize_aatest) then
        f(:,:,:,iaatest:iaatest+3*njtest-1)=0.
      endif
!
!  write testfield information to a file (for convenient post-processing)
!
      if (lroot) then
        open(1,file=trim(datadir)//'/testfield_info.dat',STATUS='unknown')
        write(1,'(a,i1)') 'xextent=',merge(1,0,xextent)
        write(1,'(a,i1)') 'zextent=',merge(1,0,zextent)
        write(1,'(a,i1)') 'lsoca='  ,merge(1,0,lsoca)
        write(1,'(a,i2)') 'itestfield=',itestfield
        write(1,'(a,f3.0)') 'ktestfield=',ktestfield
        close(1)
      endif
!
    endsubroutine initialize_testfield
!***********************************************************************
    subroutine init_aatest(f)
!
!  initialise testfield; called from start.f90
!
!   2-jun-05/axel: adapted from magnetic
!
      use InitialCondition, only: initial_condition_aatest
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      select case (initaatest)

      case ('zero', '0'); f(:,:,:,iaatest:iaatest+3*njtest-1)=0.

      case default
        call fatal_error("init_aatest",'no such initaatest: '//trim(initaatest))
      endselect
!
!  Interface for user's own subroutine
!
      if (linitial_condition) call initial_condition_aatest(f)
!
    endsubroutine init_aatest
!***********************************************************************
    subroutine pencil_criteria_testfield
!
!   All pencils that the Testfield module depends on are specified here.
!
!  26-jun-05/anders: adapted from magnetic
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
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_testfield
!***********************************************************************
    subroutine read_testfield_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=testfield_init_pars, IOSTAT=iostat)
!
    endsubroutine read_testfield_init_pars
!***********************************************************************
    subroutine write_testfield_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=testfield_init_pars)
!
    endsubroutine write_testfield_init_pars
!***********************************************************************
    subroutine read_testfield_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=testfield_run_pars, IOSTAT=iostat)
!
    endsubroutine read_testfield_run_pars
!***********************************************************************
    subroutine write_testfield_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=testfield_run_pars)
!
    endsubroutine write_testfield_run_pars
!***********************************************************************
    subroutine daatest_dt(f,df,p)
!
!  testfield evolution
!  calculate da^(q)/dt=uxB^(q)+eta*del2A^(q), where q=1,...,9
!
!   3-jun-05/axel: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      real, dimension (nx,3) :: uxB,bbtest,btest,uxbtest,duxbtest
      real, dimension (nx,3) :: del2Atest
      real, dimension (nx) :: diffus_eta
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
        if ((jtest>= 1.and.jtest<= 3)&
        .or.(jtest>= 4.and.jtest<= 6.and.zextent)&
        .or.(jtest>=10.and.jtest<=12.and.xextent.and.(.not.lset_bbtest2))&
        .or.(jtest>= 7.and.jtest<= 9.and.xextent)) then
!       if ((jtest>= 1.and.jtest<= 3)&
!       .or.(jtest>= 4.and.jtest<= 6.and.xextent)&
!       .or.(jtest>=10.and.jtest<=12.and.xextent.and.(.not.lset_bbtest2))&
!       .or.(jtest>= 7.and.jtest<= 9.and.zextent)) then
          iaxtest=iaatest+3*(jtest-1)
          iaztest=iaxtest+2
          call del2v(f,iaxtest,del2Atest)
          select case (itestfield)
          case (1); call set_bbtest(bbtest,jtest,ktestfield)
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
              duxbtest(:,j)=uxbtest(:,j)-uxbtestm(n,j,jtest)
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
            diffus_eta=etatest*dxyz_2
            maxdiffus = max(maxdiffus,diffus_eta)
          endif
        endif
      enddo

      call calc_diagnostics_testfield(f,p)
!
    endsubroutine daatest_dt
!***********************************************************************
    subroutine calc_diagnostics_testfield(f,p)
!
      use Diagnostics
      use Sub, only: cross_mn, curl

      type (pencil_case) :: p
      real, dimension (mx,my,mz,mfarray) :: f

      real, dimension (nx,3) :: btest,uxbtest
      integer :: jtest

      do jtest=1,njtest

        if ((jtest>= 1.and.jtest<= 3) &
        .or.(jtest>= 4.and.jtest<= 6.and.zextent) &
        .or.(jtest>=10.and.jtest<=12.and.xextent.and.(.not.lset_bbtest2)) &
        .or.(jtest>= 7.and.jtest<= 9.and.xextent)) then
!       if ((jtest>= 1.and.jtest<= 3)&
!       .or.(jtest>= 4.and.jtest<= 6.and.xextent)&
!       .or.(jtest>=10.and.jtest<=12.and.xextent.and.(.not.lset_bbtest2))&
!       .or.(jtest>= 7.and.jtest<= 9.and.zextent)) then
          iaxtest=iaatest+3*(jtest-1)
          iaztest=iaxtest+2
!
!  calculate alpha, begin by calculating uxbtest (if not already done above)
!
          if ((ldiagnos.or.l1davgfirst).and.lsoca) then
            call curl(f,iaxtest,btest)
            call cross_mn(p%uu,btest,uxbtest)
           endif
!
!  in the following block, we have already swapped the 4-6 entries with 7-9
!
          if (ldiagnos) then  
            select case (jtest)
            case (1)
              if (idiag_alp11/=0) call sum_mn_name(uxbtest(:,1),idiag_alp11)
              if (idiag_alp21/=0) call sum_mn_name(uxbtest(:,2),idiag_alp21)
              if (idiag_alp31/=0) call sum_mn_name(uxbtest(:,3),idiag_alp31)
              if (idiag_alp11xz/=0) call ysum_mn_name_xz(uxbtest(:,1),idiag_alp11xz)
              if (idiag_alp21xz/=0) call ysum_mn_name_xz(uxbtest(:,2),idiag_alp21xz)
              if (idiag_alp31xz/=0) call ysum_mn_name_xz(uxbtest(:,3),idiag_alp31xz)
            case (2)
              if (idiag_alp12/=0) call sum_mn_name(uxbtest(:,1),idiag_alp12)
              if (idiag_alp22/=0) call sum_mn_name(uxbtest(:,2),idiag_alp22)
              if (idiag_alp32/=0) call sum_mn_name(uxbtest(:,3),idiag_alp32)
              if (idiag_alp12xz/=0) call ysum_mn_name_xz(uxbtest(:,1),idiag_alp12xz)
              if (idiag_alp22xz/=0) call ysum_mn_name_xz(uxbtest(:,2),idiag_alp22xz)
              if (idiag_alp32xz/=0) call ysum_mn_name_xz(uxbtest(:,3),idiag_alp32xz)
            case (3)
              if (idiag_alp13/=0) call sum_mn_name(uxbtest(:,1),idiag_alp13)
              if (idiag_alp23/=0) call sum_mn_name(uxbtest(:,2),idiag_alp23)
              if (idiag_alp33/=0) call sum_mn_name(uxbtest(:,3),idiag_alp33)
              if (idiag_alp13xz/=0) call ysum_mn_name_xz(uxbtest(:,1),idiag_alp13xz)
              if (idiag_alp23xz/=0) call ysum_mn_name_xz(uxbtest(:,2),idiag_alp23xz)
              if (idiag_alp33xz/=0) call ysum_mn_name_xz(uxbtest(:,3),idiag_alp33xz)
            case (4)
              if (idiag_eta113xz/=0) call ysum_mn_name_xz(uxbtest(:,1),idiag_eta113xz)
              if (idiag_eta213xz/=0) call ysum_mn_name_xz(uxbtest(:,2),idiag_eta213xz)
              if (idiag_eta313xz/=0) call ysum_mn_name_xz(uxbtest(:,3),idiag_eta313xz)
            case (5)
              if (idiag_eta123xz/=0) call ysum_mn_name_xz(uxbtest(:,1),idiag_eta123xz)
              if (idiag_eta223xz/=0) call ysum_mn_name_xz(uxbtest(:,2),idiag_eta223xz)
              if (idiag_eta323xz/=0) call ysum_mn_name_xz(uxbtest(:,3),idiag_eta323xz)
            case (6)
              if (idiag_eta133xz/=0) call ysum_mn_name_xz(uxbtest(:,1),idiag_eta133xz)
              if (idiag_eta233xz/=0) call ysum_mn_name_xz(uxbtest(:,2),idiag_eta233xz)
              if (idiag_eta333xz/=0) call ysum_mn_name_xz(uxbtest(:,3),idiag_eta333xz)
            case (7)
              if (idiag_eta111xz/=0) call ysum_mn_name_xz(uxbtest(:,1),idiag_eta111xz)
              if (idiag_eta211xz/=0) call ysum_mn_name_xz(uxbtest(:,2),idiag_eta211xz)
              if (idiag_eta311xz/=0) call ysum_mn_name_xz(uxbtest(:,3),idiag_eta311xz)
            case (8)              
              if (idiag_eta121xz/=0) call ysum_mn_name_xz(uxbtest(:,1),idiag_eta121xz)
              if (idiag_eta221xz/=0) call ysum_mn_name_xz(uxbtest(:,2),idiag_eta221xz)
              if (idiag_eta321xz/=0) call ysum_mn_name_xz(uxbtest(:,3),idiag_eta321xz)
            case (9)
              if (idiag_eta131xz/=0) call ysum_mn_name_xz(uxbtest(:,1),idiag_eta131xz)
              if (idiag_eta231xz/=0) call ysum_mn_name_xz(uxbtest(:,2),idiag_eta231xz)
              if (idiag_eta331xz/=0) call ysum_mn_name_xz(uxbtest(:,3),idiag_eta331xz)
            case (10)
              if (idiag_alp11exz/=0) call ysum_mn_name_xz(uxbtest(:,1),idiag_alp11exz)
              if (idiag_alp21exz/=0) call ysum_mn_name_xz(uxbtest(:,2),idiag_alp21exz)
              if (idiag_alp31exz/=0) call ysum_mn_name_xz(uxbtest(:,3),idiag_alp31exz)
            case (11)
              if (idiag_alp12exz/=0) call ysum_mn_name_xz(uxbtest(:,1),idiag_alp12exz)
              if (idiag_alp22exz/=0) call ysum_mn_name_xz(uxbtest(:,2),idiag_alp22exz)
              if (idiag_alp32exz/=0) call ysum_mn_name_xz(uxbtest(:,3),idiag_alp32exz)
            case (12)
              if (idiag_alp13exz/=0) call ysum_mn_name_xz(uxbtest(:,1),idiag_alp13exz)
              if (idiag_alp23exz/=0) call ysum_mn_name_xz(uxbtest(:,2),idiag_alp23exz)
              if (idiag_alp33exz/=0) call ysum_mn_name_xz(uxbtest(:,3),idiag_alp33exz)
            end select
          endif
!
!  in the following block, we have already swapped the 4-6 entries with 7-9
!
          if (l1davgfirst) then
            select case (jtest)
            case (1)
              if (idiag_alp11z/=0) call xysum_mn_name_z(uxbtest(:,1),idiag_alp11z)
              if (idiag_alp21z/=0) call xysum_mn_name_z(uxbtest(:,2),idiag_alp21z)
              if (idiag_alp31z/=0) call xysum_mn_name_z(uxbtest(:,3),idiag_alp31z)
            case (2)
              if (idiag_alp12z/=0) call xysum_mn_name_z(uxbtest(:,1),idiag_alp12z)
              if (idiag_alp22z/=0) call xysum_mn_name_z(uxbtest(:,2),idiag_alp22z)
              if (idiag_alp32z/=0) call xysum_mn_name_z(uxbtest(:,3),idiag_alp32z)
            case (3)
              if (idiag_alp13z/=0) call xysum_mn_name_z(uxbtest(:,1),idiag_alp13z)
              if (idiag_alp23z/=0) call xysum_mn_name_z(uxbtest(:,2),idiag_alp23z)
              if (idiag_alp33z/=0) call xysum_mn_name_z(uxbtest(:,3),idiag_alp33z)
            case (4)
              if (idiag_eta113z/=0) call xysum_mn_name_z(uxbtest(:,1),idiag_eta113z)
              if (idiag_eta213z/=0) call xysum_mn_name_z(uxbtest(:,2),idiag_eta213z)
              if (idiag_eta313z/=0) call xysum_mn_name_z(uxbtest(:,3),idiag_eta313z)
            case (5)
              if (idiag_eta123z/=0) call xysum_mn_name_z(uxbtest(:,1),idiag_eta123z)
              if (idiag_eta223z/=0) call xysum_mn_name_z(uxbtest(:,2),idiag_eta223z)
              if (idiag_eta323z/=0) call xysum_mn_name_z(uxbtest(:,3),idiag_eta323z)
            case (6)
              if (idiag_eta133z/=0) call xysum_mn_name_z(uxbtest(:,1),idiag_eta133z)
              if (idiag_eta233z/=0) call xysum_mn_name_z(uxbtest(:,2),idiag_eta233z)
              if (idiag_eta333z/=0) call xysum_mn_name_z(uxbtest(:,3),idiag_eta333z)
            case (7)
              if (idiag_eta111z/=0) call xysum_mn_name_z(uxbtest(:,1),idiag_eta111z)
              if (idiag_eta211z/=0) call xysum_mn_name_z(uxbtest(:,2),idiag_eta211z)
              if (idiag_eta311z/=0) call xysum_mn_name_z(uxbtest(:,3),idiag_eta311z)
            case (8)
              if (idiag_eta121z/=0) call xysum_mn_name_z(uxbtest(:,1),idiag_eta121z)
              if (idiag_eta221z/=0) call xysum_mn_name_z(uxbtest(:,2),idiag_eta221z)
              if (idiag_eta321z/=0) call xysum_mn_name_z(uxbtest(:,3),idiag_eta321z)
            case (9)
              if (idiag_eta131z/=0) call xysum_mn_name_z(uxbtest(:,1),idiag_eta131z)
              if (idiag_eta231z/=0) call xysum_mn_name_z(uxbtest(:,2),idiag_eta231z)
              if (idiag_eta331z/=0) call xysum_mn_name_z(uxbtest(:,3),idiag_eta331z)
            end select
          endif
        endif
      enddo
!
      call keep_compiler_quiet(p)
!
    endsubroutine calc_diagnostics_testfield
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
    subroutine testfield_before_boundary(f)
!
!  Actions to take before boundary conditions are set.
!
!    4-oct-18/axel+nishant: adapted from testflow
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine testfield_before_boundary
!***********************************************************************
    subroutine testfield_after_boundary(f)
!
!  calculate <uxb>, which is needed when lsoca=.false.
!
!  21-jan-06/axel: coded
!
      use Sub
      use Hydro, only: calc_pencils_hydro
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: btest,uxbtest
      integer :: jtest,j,nxy=nxgrid*nygrid
      logical :: headtt_save
      real :: fac
!
!  In this routine we will reset headtt after the first pencil,
!  so we need to reset it afterwards.
!
      headtt_save=headtt
      fac=1./nxy
!
!  do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further up in the file.
!
      do jtest=1,njtest
        if ((jtest>= 1.and.jtest<= 3)&
        .or.(jtest>= 4.and.jtest<= 6.and.zextent)&
        .or.(jtest>=10.and.jtest<=12.and.xextent.and.(.not.lset_bbtest2))&
        .or.(jtest>= 7.and.jtest<= 9.and.xextent)) then
!       if ((jtest>= 1.and.jtest<= 3)&
!       .or.(jtest>= 4.and.jtest<= 6.and.xextent)&
!       .or.(jtest>=10.and.jtest<=12.and.xextent.and.(.not.lset_bbtest2))&
!       .or.(jtest>= 7.and.jtest<= 9.and.zextent)) then
          iaxtest=iaatest+3*(jtest-1)
          iaztest=iaxtest+2
          if (lsoca) then
            uxbtestm(:,:,jtest)=0.
          else
            do n=n1,n2
              uxbtestm(n,:,jtest)=0.
              do m=m1,m2
                call calc_pencils_hydro(f,p)
                call curl(f,iaxtest,btest)
                call cross_mn(p%uu,btest,uxbtest)
                do j=1,3
                  uxbtestm(n,j,jtest)=uxbtestm(n,j,jtest)+fac*sum(uxbtest(:,j))
                enddo
                headtt=.false.
              enddo
            enddo
!
!  do communication for array of size mz,3,ntestfield/3=mz*ntestfield
!  (Could do the same in momentum removal procedure.)
!
!  real, dimension (mz,3,ntestfield/3) :: fsum_tmp,fsum
!           fsum_tmp=uxbtestm
!           call mpiallreduce_sum(fsum_tmp,fsum)
!           uxbtestm=fsum
!
          endif
        endif
      enddo
!
!  reset headtt
!
      headtt=headtt_save
!
    endsubroutine testfield_after_boundary
!***********************************************************************
    subroutine rescaling_testfield(f)
!
!  Rescale testfield by factor rescale_aatest(jtest),
!  which could be different for different testfields
!
!  18-may-08/axel: rewrite from rescaling as used in magnetic
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=fnlen) :: file
      logical :: ltestfield_out
      integer,save :: ifirst=0
      integer :: j,jtest
!
      intent(inout) :: f
!
! reinitialize aatest periodically if requested
!
      if (linit_aatest) then
        file=trim(datadir)//'/tinit_aatest.dat'
        if (ifirst==0) then
          call read_snaptime(trim(file),taainit,naainit,daainit,t)
          if (taainit==0 .or. taainit < t-daainit) then
            taainit=t+daainit
          endif
          ifirst=1
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
    subroutine set_bbtest(bbtest,jtest,ktestfield)
!
!  set testfield
!
!   3-jun-05/axel: coded
!
      use Sub
!
      real, dimension (nx,3) :: bbtest
      real, dimension (nx) :: cx,sx,cz,sz
      integer :: jtest
      real :: ktestfield
!
      intent(in)  :: jtest,ktestfield
      intent(out) :: bbtest
!
!  xx and zz for calculating diffusive part of emf
!
      cx=cos(ktestfield*x(l1:l2))
      sx=sin(ktestfield*x(l1:l2))
      cz=cos(ktestfield*z(n))
      sz=sin(ktestfield*z(n))
!
!  set bbtest for each of the 9 cases
!
      select case (jtest)
      case (1); bbtest(:,1)=sz; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (2); bbtest(:,1)=0.; bbtest(:,2)=sz; bbtest(:,3)=0.
      case (3); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sz
      case (4); bbtest(:,1)=cz; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (5); bbtest(:,1)=0.; bbtest(:,2)=cz; bbtest(:,3)=0.
      case (6); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cz
      case (7); bbtest(:,1)=cx; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (8); bbtest(:,1)=0.; bbtest(:,2)=cx; bbtest(:,3)=0.
      case (9); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cx
      case (10); bbtest(:,1)=sx; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (11); bbtest(:,1)=0.; bbtest(:,2)=sx; bbtest(:,3)=0.
      case (12); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sx
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
      use Diagnostics
!
      integer :: iname,inamez,inamexz
      logical :: lreset
      logical, optional :: lwrite
!
!  reset everything in case of RELOAD
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_alp11=0; idiag_alp21=0; idiag_alp31=0
        idiag_alp12=0; idiag_alp22=0; idiag_alp32=0
        idiag_alp13=0; idiag_alp23=0; idiag_alp33=0
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
        idiag_eta121xz=0; idiag_eta221xz=0; idiag_eta321xz=0
        idiag_eta131xz=0; idiag_eta231xz=0; idiag_eta331xz=0
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
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta211xz',idiag_eta211xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta311xz',idiag_eta311xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta121xz',idiag_eta121xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta221xz',idiag_eta221xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta321xz',idiag_eta321xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta131xz',idiag_eta131xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta231xz',idiag_eta231xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta331xz',idiag_eta331xz)
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
    endsubroutine rprint_testfield

endmodule Testfield
