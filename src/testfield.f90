! $Id: testfield.f90,v 1.2 2005-06-06 20:25:09 brandenb Exp $

!  This modules deals with all aspects of testfield fields; if no
!  testfield fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testfield relevant subroutines listed in here.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 27
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Testfield

  use Cparam

  implicit none

  character (len=labellen) :: initaatest='zero'
  ! input parameters
  real, dimension (nx,3) :: bbb
  real :: amplaa=0., kx_aa=1.,ky_aa=1.,kz_aa=1.
  logical :: reinitalize_aatest=.false.

  namelist /testfield_init_pars/ &
       initaatest

  ! run parameters
  real :: etatest=0.
  namelist /testfield_run_pars/ &
       etatest,reinitalize_aatest

  ! other variables (needs to be consistent with reset list below)
  integer :: i_alp11=0,i_alp21=0,i_alp31=0
  integer :: i_alp12=0,i_alp22=0,i_alp32=0
  integer :: i_alp13=0,i_alp23=0,i_alp33=0
  integer :: i_alp11z=0,i_alp21z=0,i_alp31z=0
  integer :: i_alp12z=0,i_alp22z=0,i_alp32z=0
  integer :: i_alp13z=0,i_alp23z=0,i_alp33z=0

  contains

!***********************************************************************
    subroutine register_testfield()
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaatest, etc; increase nvar accordingly
!
!  1-may-02/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_aa called twice')
      first = .false.
!
      ltestfield = .true.
      iaatest = nvar+1              ! indices to access aa
      nvar = nvar+27            ! added 27 variables
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_testfield: nvar = ', nvar
        print*, 'register_testfield: iaatest = ', iaatest
      endif
!
!  Put variable names in array
!
      varname(iax) = 'ax'
      varname(iay) = 'ay'
      varname(iaz) = 'az'
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: testfield.f90,v 1.2 2005-06-06 20:25:09 brandenb Exp $")
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
          if (nvar < mvar) write(4,*) ',aa $'
          if (nvar == mvar) write(4,*) ',aa'
        else
          write(4,*) ',aa $'
        endif
        write(15,*) 'aa = fltarr(mx,my,mz,3)*one'
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
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
!  set to zero and then rescale the testfield
!  (in future, could call something like init_aa_simple)
!
      if (reinitalize_aatest) then
        f(:,:,:,iaatest:iaatest+26)=0.
      endif
!
    endsubroutine initialize_testfield
!***********************************************************************
    subroutine init_aatest(f,xx,yy,zz)
!
!  initialise testfield; called from start.f90
!
!   2-jun-05/axel: adapted from magnetic
!
      use Cdata
      use Mpicomm
      use Density
      use Gravity, only: gravz
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz,tmp,prof
      real, dimension (nx,3) :: bb
      real, dimension (nx) :: b2,fact
      real :: beq2
!
      select case(initaatest)

      case('zero', '0'); f(:,:,:,iaatest:iaatest+26)=0.

      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'init_aatest: check initaatest: ', trim(initaatest)
        call stop_it("")

      endselect
!
    endsubroutine init_aatest
!***********************************************************************
    subroutine daatest_dt(f,df,uu)
!
!  testfield evolution
!  calculate dA/dt=uxBtest+eta*del2Atest
!
!   3-jun-05/axel: coded
!
      use Cdata
      use Sub
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: bb,aa,uxB,uu,bbtest,btest,uxbtest
      real, dimension (nx,3) :: del2Atest
      integer :: jtest
!
      intent(in)     :: f,uu
      intent(inout)  :: df     
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'daatest_dt: SOLVE'
      if (headtt) then
        call identify_bcs('Axtest',iaxtest)
        call identify_bcs('Aytest',iaytest)
        call identify_bcs('Aztest',iaztest)
      endif
!
!  do each test field at a time
!
      do jtest=1,9
        iaxtest=iaatest+3*(jtest-1)
        iaztest=iaxtest+2
        call del2v(f,iaxtest,del2Atest)
        call set_bbtest(bbtest,jtest)
        call cross_mn(uu,bbtest,uxB)
        df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest) &
          +uxB+etatest*del2Atest
!
!  calculate uxbtest
!
        call curl(f,iaxtest,btest)
        call cross_mn(uu,btest,uxbtest)
!
!  calculate alpha
!
        if (ldiagnos) then
          select case(jtest)
          case(1)
            if (i_alp11/=0) call sum_mn_name(uxbtest(:,1),i_alp11)
            if (i_alp21/=0) call sum_mn_name(uxbtest(:,2),i_alp21)
            if (i_alp31/=0) call sum_mn_name(uxbtest(:,3),i_alp31)
            if (i_alp11z/=0) call xysum_mn_name_z(uxbtest(:,1),i_alp11z)
          case(2)
            if (i_alp12/=0) call sum_mn_name(uxbtest(:,1),i_alp12)
            if (i_alp22/=0) call sum_mn_name(uxbtest(:,2),i_alp22)
            if (i_alp32/=0) call sum_mn_name(uxbtest(:,3),i_alp32)
          case(3)
            if (i_alp13/=0) call sum_mn_name(uxbtest(:,1),i_alp13)
            if (i_alp23/=0) call sum_mn_name(uxbtest(:,2),i_alp23)
            if (i_alp33/=0) call sum_mn_name(uxbtest(:,3),i_alp33)
          end select
        endif
      enddo
!
    endsubroutine daatest_dt
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
      select case(jtest)
      case(1); bbtest(:,1)=1.; bbtest(:,2)=0.; bbtest(:,3)=0.
      case(2); bbtest(:,1)=0.; bbtest(:,2)=1.; bbtest(:,3)=0.
      case(3); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=1.
      case default; bbtest(:,:)=0.
      endselect
!
    endsubroutine set_bbtest
!***********************************************************************
    subroutine rprint_testfield(lreset,lwrite)
!
!  reads and registers print parameters relevant for testfield fields
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Cdata
      use Sub
!
      integer :: iname,inamez,ixy,irz
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
        i_alp11=0; i_alp21=0; i_alp31=0
        i_alp12=0; i_alp22=0; i_alp32=0
        i_alp13=0; i_alp23=0; i_alp33=0
        i_alp11z=0; i_alp21z=0; i_alp31z=0
        i_alp12z=0; i_alp22z=0; i_alp32z=0
        i_alp13z=0; i_alp23z=0; i_alp33z=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'alp11',i_alp11)
        call parse_name(iname,cname(iname),cform(iname),'alp21',i_alp21)
        call parse_name(iname,cname(iname),cform(iname),'alp31',i_alp31)
        call parse_name(iname,cname(iname),cform(iname),'alp12',i_alp12)
        call parse_name(iname,cname(iname),cform(iname),'alp22',i_alp22)
        call parse_name(iname,cname(iname),cform(iname),'alp32',i_alp32)
        call parse_name(iname,cname(iname),cform(iname),'alp13',i_alp13)
        call parse_name(iname,cname(iname),cform(iname),'alp23',i_alp23)
        call parse_name(iname,cname(iname),cform(iname),'alp33',i_alp33)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp11z',i_alp11z)
      enddo
!
!  check for those quantities for which we want z-averages
!
!     do ixy=1,nnamexy
!       call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bxmxy',i_bxmxy)
!       call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bymxy',i_bymxy)
!       call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bzmxy',i_bzmxy)
!     enddo
!
!  write column, i_XYZ, where our variable XYZ is stored
!
      if (lwr) then
        write(3,*) 'i_alp11=',i_alp11
        write(3,*) 'i_alp21=',i_alp21
        write(3,*) 'i_alp31=',i_alp31
        write(3,*) 'i_alp12=',i_alp12
        write(3,*) 'i_alp22=',i_alp22
        write(3,*) 'i_alp32=',i_alp32
        write(3,*) 'i_alp13=',i_alp13
        write(3,*) 'i_alp23=',i_alp23
        write(3,*) 'i_alp33=',i_alp33
        write(3,*) 'i_alp11z=',i_alp11z
        write(3,*) 'iaatest=',iaatest
        write(3,*) 'nnamez=',nnamez
        write(3,*) 'nnamexy=',nnamexy
      endif
!
    endsubroutine rprint_testfield

endmodule Testfield
