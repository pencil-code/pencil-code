! $Id: testfield.f90,v 1.10 2005-06-13 19:25:55 brandenb Exp $

!  This modules deals with all aspects of testfield fields; if no
!  testfield fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testfield relevant subroutines listed in here.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 36
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
  logical :: xextent=.true.,zextent=.true.,lsoca=.true.,lset_bbtest2=.false.

  namelist /testfield_init_pars/ &
       xextent,zextent,initaatest

  ! run parameters
  real :: etatest=0.
  namelist /testfield_run_pars/ &
       reinitalize_aatest,xextent,zextent,lsoca, &
       lset_bbtest2,etatest

  ! other variables (needs to be consistent with reset list below)
  integer :: i_alp11=0,i_alp21=0,i_alp31=0
  integer :: i_alp12=0,i_alp22=0,i_alp32=0
  integer :: i_alp13=0,i_alp23=0,i_alp33=0
  integer :: i_alp11z=0,i_alp21z=0,i_alp31z=0
  integer :: i_alp12z=0,i_alp22z=0,i_alp32z=0
  integer :: i_alp13z=0,i_alp23z=0,i_alp33z=0
  integer :: i_eta111z=0,i_eta211z=0,i_eta311z=0
  integer :: i_eta121z=0,i_eta221z=0,i_eta321z=0
  integer :: i_eta131z=0,i_eta231z=0,i_eta331z=0
  integer :: i_eta113z=0,i_eta213z=0,i_eta313z=0
  integer :: i_eta123z=0,i_eta223z=0,i_eta323z=0
  integer :: i_eta133z=0,i_eta233z=0,i_eta333z=0
  integer :: i_alp11xz=0,i_alp21xz=0,i_alp31xz=0
  integer :: i_alp12xz=0,i_alp22xz=0,i_alp32xz=0
  integer :: i_alp13xz=0,i_alp23xz=0,i_alp33xz=0
  integer :: i_eta111xz=0,i_eta211xz=0,i_eta311xz=0
  integer :: i_eta121xz=0,i_eta221xz=0,i_eta321xz=0
  integer :: i_eta131xz=0,i_eta231xz=0,i_eta331xz=0
  integer :: i_eta113xz=0,i_eta213xz=0,i_eta313xz=0
  integer :: i_eta123xz=0,i_eta223xz=0,i_eta323xz=0
  integer :: i_eta133xz=0,i_eta233xz=0,i_eta333xz=0
  integer :: i_alp11exz=0,i_alp21exz=0,i_alp31exz=0
  integer :: i_alp12exz=0,i_alp22exz=0,i_alp32exz=0
  integer :: i_alp13exz=0,i_alp23exz=0,i_alp33exz=0

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
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_aa called twice')
      first = .false.
!
      ltestfield = .true.
      iaatest = nvar+1              ! indices to access aa
!---  nvar = nvar+27            ! added 27 variables
      nvar = nvar+36            ! added 36 variables
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
           "$Id: testfield.f90,v 1.10 2005-06-13 19:25:55 brandenb Exp $")
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
!  calculate da^(q)/dt=uxB^(q)+eta*del2A^(q), where q=1,...,9
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
!  do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!
      do jtest=1,12
        if ((jtest>= 1.and.jtest<= 3)&
        .or.(jtest>= 4.and.jtest<= 6.and.xextent)&
        .or.(jtest>=10.and.jtest<=12.and.xextent)&
        .or.(jtest>= 7.and.jtest<= 9.and.zextent)) then
          iaxtest=iaatest+3*(jtest-1)
          iaztest=iaxtest+2
          call del2v(f,iaxtest,del2Atest)
          if (lset_bbtest2) then
            call set_bbtest2(bbtest,jtest)
          else
            call set_bbtest(bbtest,jtest)
          endif
          call cross_mn(uu,bbtest,uxB)
          if (lsoca) then
            df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest) &
              +uxB+etatest*del2Atest
          else
            call curl(f,iaxtest,btest)
            call cross_mn(uu,btest,uxbtest)
            df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest) &
              +uxB+etatest*del2Atest+uxbtest
          endif
!
!  calculate alpha, begin by calculating uxbtest (if not already done above)
!
          if (ldiagnos) then
            if (lsoca) then
              call curl(f,iaxtest,btest)
              call cross_mn(uu,btest,uxbtest)
            endif
            select case(jtest)
            case(1)
              if (i_alp11/=0) call sum_mn_name(uxbtest(:,1),i_alp11)
              if (i_alp21/=0) call sum_mn_name(uxbtest(:,2),i_alp21)
              if (i_alp31/=0) call sum_mn_name(uxbtest(:,3),i_alp31)
              if (i_alp11z/=0) call xysum_mn_name_z(uxbtest(:,1),i_alp11z)
              if (i_alp21z/=0) call xysum_mn_name_z(uxbtest(:,2),i_alp21z)
              if (i_alp31z/=0) call xysum_mn_name_z(uxbtest(:,3),i_alp31z)
              if (i_alp11xz/=0) call ysum_mn_name_xz(uxbtest(:,1),i_alp11xz)
              if (i_alp21xz/=0) call ysum_mn_name_xz(uxbtest(:,2),i_alp21xz)
              if (i_alp31xz/=0) call ysum_mn_name_xz(uxbtest(:,3),i_alp31xz)
            case(2)
              if (i_alp12/=0) call sum_mn_name(uxbtest(:,1),i_alp12)
              if (i_alp22/=0) call sum_mn_name(uxbtest(:,2),i_alp22)
              if (i_alp32/=0) call sum_mn_name(uxbtest(:,3),i_alp32)
              if (i_alp12z/=0) call xysum_mn_name_z(uxbtest(:,1),i_alp12z)
              if (i_alp22z/=0) call xysum_mn_name_z(uxbtest(:,2),i_alp22z)
              if (i_alp32z/=0) call xysum_mn_name_z(uxbtest(:,3),i_alp32z)
              if (i_alp12xz/=0) call ysum_mn_name_xz(uxbtest(:,1),i_alp12xz)
              if (i_alp22xz/=0) call ysum_mn_name_xz(uxbtest(:,2),i_alp22xz)
              if (i_alp32xz/=0) call ysum_mn_name_xz(uxbtest(:,3),i_alp32xz)
            case(3)
              if (i_alp13/=0) call sum_mn_name(uxbtest(:,1),i_alp13)
              if (i_alp23/=0) call sum_mn_name(uxbtest(:,2),i_alp23)
              if (i_alp33/=0) call sum_mn_name(uxbtest(:,3),i_alp33)
              if (i_alp13z/=0) call xysum_mn_name_z(uxbtest(:,1),i_alp13z)
              if (i_alp23z/=0) call xysum_mn_name_z(uxbtest(:,2),i_alp23z)
              if (i_alp33z/=0) call xysum_mn_name_z(uxbtest(:,3),i_alp33z)
              if (i_alp13xz/=0) call ysum_mn_name_xz(uxbtest(:,1),i_alp13xz)
              if (i_alp23xz/=0) call ysum_mn_name_xz(uxbtest(:,2),i_alp23xz)
              if (i_alp33xz/=0) call ysum_mn_name_xz(uxbtest(:,3),i_alp33xz)
            case(4)
              if (i_eta111z/=0) call xysum_mn_name_z(uxbtest(:,1),i_eta111z)
              if (i_eta211z/=0) call xysum_mn_name_z(uxbtest(:,2),i_eta211z)
              if (i_eta311z/=0) call xysum_mn_name_z(uxbtest(:,3),i_eta311z)
              if (i_eta111xz/=0) call ysum_mn_name_xz(uxbtest(:,1),i_eta111xz)
              if (i_eta211xz/=0) call ysum_mn_name_xz(uxbtest(:,2),i_eta211xz)
              if (i_eta311xz/=0) call ysum_mn_name_xz(uxbtest(:,3),i_eta311xz)
            case(5)
              if (i_eta121z/=0) call xysum_mn_name_z(uxbtest(:,1),i_eta121z)
              if (i_eta221z/=0) call xysum_mn_name_z(uxbtest(:,2),i_eta221z)
              if (i_eta321z/=0) call xysum_mn_name_z(uxbtest(:,3),i_eta321z)
              if (i_eta121xz/=0) call ysum_mn_name_xz(uxbtest(:,1),i_eta121xz)
              if (i_eta221xz/=0) call ysum_mn_name_xz(uxbtest(:,2),i_eta221xz)
              if (i_eta321xz/=0) call ysum_mn_name_xz(uxbtest(:,3),i_eta321xz)
            case(6)
              if (i_eta131z/=0) call xysum_mn_name_z(uxbtest(:,1),i_eta131z)
              if (i_eta231z/=0) call xysum_mn_name_z(uxbtest(:,2),i_eta231z)
              if (i_eta331z/=0) call xysum_mn_name_z(uxbtest(:,3),i_eta331z)
              if (i_eta131xz/=0) call ysum_mn_name_xz(uxbtest(:,1),i_eta131xz)
              if (i_eta231xz/=0) call ysum_mn_name_xz(uxbtest(:,2),i_eta231xz)
              if (i_eta331xz/=0) call ysum_mn_name_xz(uxbtest(:,3),i_eta331xz)
            case(7)
              if (i_eta113z/=0) call xysum_mn_name_z(uxbtest(:,1),i_eta113z)
              if (i_eta213z/=0) call xysum_mn_name_z(uxbtest(:,2),i_eta213z)
              if (i_eta313z/=0) call xysum_mn_name_z(uxbtest(:,3),i_eta313z)
              if (i_eta113xz/=0) call ysum_mn_name_xz(uxbtest(:,1),i_eta113xz)
              if (i_eta213xz/=0) call ysum_mn_name_xz(uxbtest(:,2),i_eta213xz)
              if (i_eta313xz/=0) call ysum_mn_name_xz(uxbtest(:,3),i_eta313xz)
            case(8)
              if (i_eta123z/=0) call xysum_mn_name_z(uxbtest(:,1),i_eta123z)
              if (i_eta223z/=0) call xysum_mn_name_z(uxbtest(:,2),i_eta223z)
              if (i_eta323z/=0) call xysum_mn_name_z(uxbtest(:,3),i_eta323z)
              if (i_eta123xz/=0) call ysum_mn_name_xz(uxbtest(:,1),i_eta123xz)
              if (i_eta223xz/=0) call ysum_mn_name_xz(uxbtest(:,2),i_eta223xz)
              if (i_eta323xz/=0) call ysum_mn_name_xz(uxbtest(:,3),i_eta323xz)
            case(9)
              if (i_eta133z/=0) call xysum_mn_name_z(uxbtest(:,1),i_eta133z)
              if (i_eta233z/=0) call xysum_mn_name_z(uxbtest(:,2),i_eta233z)
              if (i_eta333z/=0) call xysum_mn_name_z(uxbtest(:,3),i_eta333z)
              if (i_eta133xz/=0) call ysum_mn_name_xz(uxbtest(:,1),i_eta133xz)
              if (i_eta233xz/=0) call ysum_mn_name_xz(uxbtest(:,2),i_eta233xz)
              if (i_eta333xz/=0) call ysum_mn_name_xz(uxbtest(:,3),i_eta333xz)
            case(10)
              if (i_alp11exz/=0) call ysum_mn_name_xz(uxbtest(:,1),i_alp11exz)
              if (i_alp21exz/=0) call ysum_mn_name_xz(uxbtest(:,2),i_alp21exz)
              if (i_alp31exz/=0) call ysum_mn_name_xz(uxbtest(:,3),i_alp31exz)
            case(11)
              if (i_alp12exz/=0) call ysum_mn_name_xz(uxbtest(:,1),i_alp12exz)
              if (i_alp22exz/=0) call ysum_mn_name_xz(uxbtest(:,2),i_alp22exz)
              if (i_alp32exz/=0) call ysum_mn_name_xz(uxbtest(:,3),i_alp32exz)
            case(12)
              if (i_alp13exz/=0) call ysum_mn_name_xz(uxbtest(:,1),i_alp13exz)
              if (i_alp23exz/=0) call ysum_mn_name_xz(uxbtest(:,2),i_alp23exz)
              if (i_alp33exz/=0) call ysum_mn_name_xz(uxbtest(:,3),i_alp33exz)
            end select
          endif
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
      select case(jtest)
      case(1); bbtest(:,1)=sz; bbtest(:,2)=0.; bbtest(:,3)=0.
      case(2); bbtest(:,1)=0.; bbtest(:,2)=sz; bbtest(:,3)=0.
      case(3); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sz
      case(4); bbtest(:,1)=cx; bbtest(:,2)=0.; bbtest(:,3)=0.
      case(5); bbtest(:,1)=0.; bbtest(:,2)=cx; bbtest(:,3)=0.
      case(6); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cx
      case(7); bbtest(:,1)=cz; bbtest(:,2)=0.; bbtest(:,3)=0.
      case(8); bbtest(:,1)=0.; bbtest(:,2)=cz; bbtest(:,3)=0.
      case(9); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cz
      case(10); bbtest(:,1)=sx; bbtest(:,2)=0.; bbtest(:,3)=0.
      case(11); bbtest(:,1)=0.; bbtest(:,2)=sx; bbtest(:,3)=0.
      case(12); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sx
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
      select case(jtest)
      case(1); bbtest(:,1)=xz; bbtest(:,2)=0.; bbtest(:,3)=0.
      case(2); bbtest(:,1)=0.; bbtest(:,2)=xz; bbtest(:,3)=0.
      case(3); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=xz
      case(4); bbtest(:,1)=sx; bbtest(:,2)=0.; bbtest(:,3)=0.
      case(5); bbtest(:,1)=0.; bbtest(:,2)=sx; bbtest(:,3)=0.
      case(6); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sx
      case(7); bbtest(:,1)=sz; bbtest(:,2)=0.; bbtest(:,3)=0.
      case(8); bbtest(:,1)=0.; bbtest(:,2)=sz; bbtest(:,3)=0.
      case(9); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sz
      case default; bbtest(:,:)=0.
      endselect
!
    endsubroutine set_bbtest2
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
      select case(jtest)
      case(1); bbtest(:,1)=1.; bbtest(:,2)=0.; bbtest(:,3)=0.
      case(2); bbtest(:,1)=0.; bbtest(:,2)=1.; bbtest(:,3)=0.
      case(3); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=1.
      case(4); bbtest(:,1)=cx; bbtest(:,2)=0.; bbtest(:,3)=0.
      case(5); bbtest(:,1)=0.; bbtest(:,2)=sx; bbtest(:,3)=0.
      case(6); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sx
      case(7); bbtest(:,1)=sz; bbtest(:,2)=0.; bbtest(:,3)=0.
      case(8); bbtest(:,1)=0.; bbtest(:,2)=sz; bbtest(:,3)=0.
      case(9); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cz
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
      use Sub
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
        i_alp11=0; i_alp21=0; i_alp31=0
        i_alp12=0; i_alp22=0; i_alp32=0
        i_alp13=0; i_alp23=0; i_alp33=0
        i_alp11z=0; i_alp21z=0; i_alp31z=0
        i_alp12z=0; i_alp22z=0; i_alp32z=0
        i_alp13z=0; i_alp23z=0; i_alp33z=0
        i_eta111z=0; i_eta211z=0; i_eta311z=0
        i_eta121z=0; i_eta221z=0; i_eta321z=0
        i_eta131z=0; i_eta231z=0; i_eta331z=0
        i_eta113z=0; i_eta213z=0; i_eta313z=0
        i_eta123z=0; i_eta223z=0; i_eta323z=0
        i_eta133z=0; i_eta233z=0; i_eta333z=0
        i_alp11xz=0; i_alp21xz=0; i_alp31xz=0
        i_alp12xz=0; i_alp22xz=0; i_alp32xz=0
        i_alp13xz=0; i_alp23xz=0; i_alp33xz=0
        i_eta111xz=0; i_eta211xz=0; i_eta311xz=0
        i_eta121xz=0; i_eta221xz=0; i_eta321xz=0
        i_eta131xz=0; i_eta231xz=0; i_eta331xz=0
        i_eta113xz=0; i_eta213xz=0; i_eta313xz=0
        i_eta123xz=0; i_eta223xz=0; i_eta323xz=0
        i_eta133xz=0; i_eta233xz=0; i_eta333xz=0
        i_alp11exz=0; i_alp21exz=0; i_alp31exz=0
        i_alp12exz=0; i_alp22exz=0; i_alp32exz=0
        i_alp13exz=0; i_alp23exz=0; i_alp33exz=0
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
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp21z',i_alp21z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp31z',i_alp31z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp12z',i_alp12z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp22z',i_alp22z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp32z',i_alp32z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp13z',i_alp13z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp23z',i_alp23z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp33z',i_alp33z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta111z',i_eta111z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta211z',i_eta211z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta311z',i_eta311z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta121z',i_eta121z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta221z',i_eta221z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta321z',i_eta321z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta131z',i_eta131z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta231z',i_eta231z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta331z',i_eta331z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta113z',i_eta113z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta213z',i_eta213z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta313z',i_eta313z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta123z',i_eta123z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta223z',i_eta223z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta323z',i_eta323z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta133z',i_eta133z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta233z',i_eta233z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta333z',i_eta333z)
      enddo
!
!  check for those quantities for which we want y-averages
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp11xz',i_alp11xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp21xz',i_alp21xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp31xz',i_alp31xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp12xz',i_alp12xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp22xz',i_alp22xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp32xz',i_alp32xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp13xz',i_alp13xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp23xz',i_alp23xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp33xz',i_alp33xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta111xz',i_eta111xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta211xz',i_eta211xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta311xz',i_eta311xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta121xz',i_eta121xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta221xz',i_eta221xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta321xz',i_eta321xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta131xz',i_eta131xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta231xz',i_eta231xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta331xz',i_eta331xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta113xz',i_eta113xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta213xz',i_eta213xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta313xz',i_eta313xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta123xz',i_eta123xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta223xz',i_eta223xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta323xz',i_eta323xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta133xz',i_eta133xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta233xz',i_eta233xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'eta333xz',i_eta333xz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp11exz',i_alp11exz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp21exz',i_alp21exz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp31exz',i_alp31exz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp12exz',i_alp12exz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp22exz',i_alp22exz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp32exz',i_alp32exz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp13exz',i_alp13exz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp23exz',i_alp23exz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'alp33exz',i_alp33exz)
      enddo
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
        write(3,*) 'i_alp21z=',i_alp21z
        write(3,*) 'i_alp31z=',i_alp31z
        write(3,*) 'i_alp12z=',i_alp12z
        write(3,*) 'i_alp22z=',i_alp22z
        write(3,*) 'i_alp32z=',i_alp32z
        write(3,*) 'i_alp13z=',i_alp13z
        write(3,*) 'i_alp23z=',i_alp23z
        write(3,*) 'i_alp33z=',i_alp33z
        write(3,*) 'i_eta111z=',i_eta111z
        write(3,*) 'i_eta211z=',i_eta211z
        write(3,*) 'i_eta311z=',i_eta311z
        write(3,*) 'i_eta121z=',i_eta121z
        write(3,*) 'i_eta221z=',i_eta221z
        write(3,*) 'i_eta321z=',i_eta321z
        write(3,*) 'i_eta131z=',i_eta131z
        write(3,*) 'i_eta231z=',i_eta231z
        write(3,*) 'i_eta331z=',i_eta331z
        write(3,*) 'i_eta113z=',i_eta113z
        write(3,*) 'i_eta213z=',i_eta213z
        write(3,*) 'i_eta313z=',i_eta313z
        write(3,*) 'i_eta123z=',i_eta123z
        write(3,*) 'i_eta223z=',i_eta223z
        write(3,*) 'i_eta323z=',i_eta323z
        write(3,*) 'i_eta133z=',i_eta133z
        write(3,*) 'i_eta233z=',i_eta233z
        write(3,*) 'i_eta333z=',i_eta333z
        write(3,*) 'i_alp11xz=',i_alp11xz
        write(3,*) 'i_alp21xz=',i_alp21xz
        write(3,*) 'i_alp31xz=',i_alp31xz
        write(3,*) 'i_alp12xz=',i_alp12xz
        write(3,*) 'i_alp22xz=',i_alp22xz
        write(3,*) 'i_alp32xz=',i_alp32xz
        write(3,*) 'i_alp13xz=',i_alp13xz
        write(3,*) 'i_alp23xz=',i_alp23xz
        write(3,*) 'i_alp33xz=',i_alp33xz
        write(3,*) 'i_eta111xz=',i_eta111xz
        write(3,*) 'i_eta211xz=',i_eta211xz
        write(3,*) 'i_eta311xz=',i_eta311xz
        write(3,*) 'i_eta121xz=',i_eta121xz
        write(3,*) 'i_eta221xz=',i_eta221xz
        write(3,*) 'i_eta321xz=',i_eta321xz
        write(3,*) 'i_eta131xz=',i_eta131xz
        write(3,*) 'i_eta231xz=',i_eta231xz
        write(3,*) 'i_eta331xz=',i_eta331xz
        write(3,*) 'i_eta113xz=',i_eta113xz
        write(3,*) 'i_eta213xz=',i_eta213xz
        write(3,*) 'i_eta313xz=',i_eta313xz
        write(3,*) 'i_eta123xz=',i_eta123xz
        write(3,*) 'i_eta223xz=',i_eta223xz
        write(3,*) 'i_eta323xz=',i_eta323xz
        write(3,*) 'i_eta133xz=',i_eta133xz
        write(3,*) 'i_eta233xz=',i_eta233xz
        write(3,*) 'i_eta333xz=',i_eta333xz
        write(3,*) 'i_alp11exz=',i_alp11exz
        write(3,*) 'i_alp21exz=',i_alp21exz
        write(3,*) 'i_alp31exz=',i_alp31exz
        write(3,*) 'i_alp12exz=',i_alp12exz
        write(3,*) 'i_alp22exz=',i_alp22exz
        write(3,*) 'i_alp32exz=',i_alp32exz
        write(3,*) 'i_alp13exz=',i_alp13exz
        write(3,*) 'i_alp23exz=',i_alp23exz
        write(3,*) 'i_alp33exz=',i_alp33exz
        write(3,*) 'iaatest=',iaatest
        write(3,*) 'nnamez=',nnamez
        write(3,*) 'nnamexy=',nnamexy
        write(3,*) 'nnamexz=',nnamexz
      endif
!
    endsubroutine rprint_testfield

endmodule Testfield
