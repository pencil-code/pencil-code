! $Id: noentropy.f90,v 1.62 2004-07-03 02:13:14 theine Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Entropy

  !
  ! isothermal case; almost nothing to do
  !

  use Cparam
  use Cdata
  use Hydro

  implicit none

  integer :: dummyss          ! We cannot define empty namelists
  namelist /entropy_init_pars/ dummyss
  namelist /entropy_run_pars/  dummyss

  ! run parameters
  real :: hcond0=0.,hcond1=impossible,chi=impossible
  real :: Fbot=impossible,FbotKbot=impossible,Kbot=impossible
  real :: Ftop=impossible,FtopKtop=impossible
  logical :: lmultilayer=.true.
  logical :: lcalc_heatcond_constchi=.false.

  ! other variables (needs to be consistent with reset list below)
  integer :: i_dtc=0,i_ssm=0,i_ugradpm=0

  contains

!***********************************************************************
    subroutine register_entropy()
!
!  no energy equation is being solved; use isothermal equation of state
!  28-mar-02/axel: dummy routine, adapted from entropy.f of 6-nov-01.
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_ent called twice')
      first = .false.
!
      lentropy = .false.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: noentropy.f90,v 1.62 2004-07-03 02:13:14 theine Exp $")
!
    endsubroutine register_entropy
!***********************************************************************
    subroutine initialize_entropy(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded 
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      logical :: lstarting
!
      if (ip == 0) print*,f,lstarting ! keep compiler quiet
!
    endsubroutine initialize_entropy
!***********************************************************************
    subroutine init_ss(f,xx,yy,zz)
!
!  initialise entropy; called from start.f90
!  28-mar-02/axel: dummy routine, adapted from entropy.f of 6-nov-01.
!  24-nov-02/tony: renamed for consistancy (i.e. init_[varaible name])
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      if(ip==1) print*,f,xx,yy,zz  !(to remove compiler warnings)
    endsubroutine init_ss
!***********************************************************************
    subroutine dss_dt(f,df,uu,glnrho,divu,rho1,lnrho,cs2,TT1,shock,gshock,bb,bij)
!
!  28-mar-02/axel: dummy routine, adapted from entropy.f of 6-nov-01.
!  19-may-02/axel: added isothermal pressure gradient
!   9-jun-02/axel: pressure gradient added to du/dt already here
!
      use Density
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: bij
      real, dimension (nx,3) :: uu,glnrho,gshock,bb
      real, dimension (nx) :: lnrho,rho1,divu,cs2,TT1,uglnrho,shock,rho
      integer :: j,ju
!
      intent(in) :: f,uu,glnrho,rho1,shock,gshock
      intent(out) :: cs2,TT1  !(df is dummy)
!
!  sound speed squared and inverse temperature
!  note: this is also correct for gamma=1
!
      cs2=cs20*exp(gamma1*(lnrho-lnrho0))
      TT1=gamma1/cs2

!
!  ``cs2/dx^2'' for timestep
!
      if (lfirst.and.ldt) advec_cs2=cs2*dxyz_2
      if (headtt.or.ldebug) print*,'dss_dt: max(advec_cs2) =',maxval(advec_cs2)
!
!  subtract isothermal/polytropic pressure gradient term in momentum equation
!
      if (lhydro) then
        do j=1,3
          ju=j+iuu-1
          df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)-cs2*glnrho(:,j)
        enddo
      endif
!
!  Calculate entropy related diagnostics
!
      if (ldiagnos) then
        if (i_dtc/=0) call max_mn_name(sqrt(advec_cs2)/cdt,i_dtc,l_dt=.true.)
        if (i_ugradpm/=0) then
          rho=1./rho1
          call dot_mn(uu,glnrho,uglnrho)
          call sum_mn_name(rho*cs2*uglnrho,i_ugradpm)
        endif
      endif
!
      if(ip==1) print*,f,df,uu,divu,rho1,shock,gshock,bb,bij  !(compiler)
    endsubroutine dss_dt
!***********************************************************************
    subroutine rprint_entropy(lreset,lwrite)
!
!  reads and registers print parameters relevant to entropy
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Cdata
      use Sub
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_dtc=0; i_ssm=0; i_ugradpm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtc',i_dtc)
        call parse_name(iname,cname(iname),cform(iname),'ugradpm',i_ugradpm)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_dtc=',i_dtc
        write(3,*) 'i_ssm=',i_ssm
        write(3,*) 'i_ugradpm=',i_ugradpm
        write(3,*) 'nname=',nname
        write(3,*) 'iss=',iss
      endif
!
    endsubroutine rprint_entropy
!***********************************************************************
    subroutine heatcond(x,y,z,hcond)
!
!  calculate the heat conductivity lambda
!  NB: if you modify this profile, you *must* adapt gradloghcond below.
!
!  23-jan-02/wolf: coded
!  28-mar-02/axel: dummy routine, adapted from entropy.f of 6-nov-01.
!
      use Cdata, only: ip
!
      real, dimension (nx) :: x,y,z
      real, dimension (nx) :: hcond
      if(ip==1) print*,x,y,z,hcond  !(to remove compiler warnings)
!
    endsubroutine heatcond
!***********************************************************************
    subroutine gradloghcond(x,y,z,glhc)
!
!  calculate grad(log lambda), where lambda is the heat conductivity
!  NB: *Must* be in sync with heatcond() above.
!
!  23-jan-02/wolf: coded
!  28-mar-02/axel: dummy routine, adapted from entropy.f of 6-nov-01.
!
      use Cdata, only: ip
!
      real, dimension (nx) :: x,y,z
      real, dimension (nx,3) :: glhc
      if(ip==1) print*,x,y,z,glhc  !(to remove compiler warnings)
!
    endsubroutine gradloghcond
!***********************************************************************
endmodule Entropy
