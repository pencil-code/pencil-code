! $Id: entropy_const.f90,v 1.3 2004-06-22 04:00:02 brandenb Exp $

!  This module is for systems with spatially fixed entropy
!  distribution. This implies Ds/Dt=u.grads only, which is used
!  in Ds/Dt=(1/gamma)*Dlnp/Dt-Dlnrho/Dt, for example.
!  This procedure has been used in the context of accretion discs
!  (see von Rekowski et al., 2003, A&A 398, 825). The shock jump
!  relations are modified (see Sect 9.3.6 of Brandenburg 2003, in
!  "Computational aspects...", ed. Ferriz-Mas & Nunez, Taylor & Francis,
!  or astro-ph/0109497.

!  This current implementation has temporarly been used to immitate
!  a corona in solar context.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Entropy

  use Cparam
  use Cdata
  use Hydro

  implicit none

  real :: dummy,xss_corona,wss_corona,ss_corona
  real, dimension(nx) :: profxss,dprofxss
  character (len=labellen) :: iss_profile='nothing'

  ! parameters necessary for consistency with other modules,
  ! but not directly used.
  real :: hcond0=0.,hcond1=impossible,chi=impossible
  real :: Fbot=impossible,FbotKbot=impossible,Kbot=impossible
  real :: Ftop=impossible,FtopKtop=impossible
  real :: ss_const=1.
  logical :: lmultilayer=.true.
  logical :: lcalc_heatcond_constchi=.false.
  character (len=labellen), dimension(ninit) :: initss='nothing'

  ! input parameters
  namelist /entropy_init_pars/ &
      initss, ss_const

  ! run parameters
  namelist /entropy_run_pars/ &
      dummy

  ! other variables (needs to be consistent with reset list below)
  integer :: i_dtc=0,i_eth=0,i_ethdivum=0,i_ssm=0,i_ugradpm=0, i_ethtot=0
  integer :: i_dtchi=0
  integer :: i_ssmphi=0

  contains

!***********************************************************************
    subroutine register_entropy()
!
!  initialise variables which should know that we solve an entropy
!  equation: iss, etc; increase nvar accordingly
!
!  6-nov-01/wolf: coded
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
      lentropy = .true.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: entropy_const.f90,v 1.3 2004-06-22 04:00:02 brandenb Exp $")
!
    endsubroutine register_entropy
!***********************************************************************
    subroutine initialize_entropy()
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  21-jul-2002/wolf: coded
!
      use Cdata
!
      integer :: ierr

  ! run parameters
  namelist /entropy_initialize_pars/ &
       iss_profile,xss_corona,wss_corona,ss_corona
!
!  read initialization parameters
!
      if(lroot) print*,'read entropy_initialize_pars'
      open(1,FILE='run.in',FORM='formatted',STATUS='old')
      read(1,NML=entropy_initialize_pars,ERR=99,IOSTAT=ierr)
      close(1)
!
!  set profiles
!
      if (iss_profile=='nothing') then
        profxss=1.
      elseif (iss_profile=='corona') then
        profxss=ss_corona*.5*(1.+tanh((x(l1:l2)-xss_corona)/wss_corona))
        dprofxss=ss_corona*.5/(wss_corona*cosh((x(l1:l2)-xss_corona)/wss_corona)**2)
      endif
!
      if(lroot) then
        print*,'initialize_entropy: iss_profile=',iss_profile
        print*,'initialize_entropy: wss_corona,ss_corona=',xss_corona,wss_corona,ss_corona
      endif

      return
99    if(lroot) print*,'i/o error'
!
    endsubroutine initialize_entropy
!***********************************************************************
    subroutine init_ss(f,xx,yy,zz)
!
!  initialise entropy; called from start.f90
!  07-nov-2001/wolf: coded
!  24-nov-2002/tony: renamed for consistancy (i.e. init_[variable name]) 
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      intent(in) :: xx,yy,zz
!
      do iinit=1,ninit
        select case(initss(iinit))

          case('zero', '0'); f(:,:,:,iss) = 0.
          case('const_ss'); f(:,:,:,iss)=f(:,:,:,iss)+ss_const

        endselect
      enddo

      if(ip==0) print*,f,xx,yy  !(to keep compiler quiet)
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
      Use Ionization, only: pressure_gradient
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: bij
      real, dimension (nx,3) :: uu,glnrho,gshock,bb
      real, dimension (nx) :: lnrho,rho1,divu,cs2,TT1,uglnrho,shock,rho
      real, dimension (nx) :: ss,cp1tilde
      integer :: j,ju
!
      intent(in) :: f,uu,glnrho,rho1,shock,gshock
      intent(out) :: cs2,TT1  !(df is dummy)
!
!  sound speed squared and inverse temperature
!
      TT1=0.
      ss=profxss
      call pressure_gradient(f,cs2,cp1tilde)
      print*, it, m, n, maxval(cs2), maxval(ss)
      if (lfirst.and.ldt) call max_for_dt(cs2,maxadvec2)
!
!  subtract isothermal/polytropic pressure gradient term in momentum equation
!
      if (lhydro) then
        df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-cs2*(glnrho(:,1)+dprofxss)
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-cs2*(glnrho(:,2))
        df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-cs2*(glnrho(:,3))
      endif
!
!  Calculate entropy related diagnostics
!
      if (ldiagnos) then
        if (i_dtc/=0) call max_mn_name(sqrt(cs2)/dxmin/cdt,i_dtc,l_dt=.true.)
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
      integer :: iname,irz
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
        i_dtc=0; i_eth=0; i_ethdivum=0; i_ssm=0; i_ugradpm=0; i_ethtot=0
        i_dtchi=0
        i_ssmphi=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtc',i_dtc)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',i_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'ethtot',i_ethtot)
        call parse_name(iname,cname(iname),cform(iname),'ethdivum',i_ethdivum)
        call parse_name(iname,cname(iname),cform(iname),'eth',i_eth)
        call parse_name(iname,cname(iname),cform(iname),'ssm',i_ssm)
        call parse_name(iname,cname(iname),cform(iname),'ugradpm',i_ugradpm)
      enddo
!
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'ssmphi',i_ssmphi)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_dtc=',i_dtc
        write(3,*) 'i_dtchi=',i_dtchi
        write(3,*) 'i_ethtot=',i_ethtot
        write(3,*) 'i_ethdivum=',i_ethdivum
        write(3,*) 'i_eth=',i_eth
        write(3,*) 'i_ssm=',i_ssm
        write(3,*) 'i_ugradpm=',i_ugradpm
        write(3,*) 'nname=',nname
        write(3,*) 'iss=',iss
        write(3,*) 'i_ssmphi=',i_ssmphi
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

