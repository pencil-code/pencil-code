! $Id: noentropy.f90,v 1.33 2002-11-12 16:07:59 dobler Exp $

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
  real, dimension (nx) :: cs2,TT1 ! Can't make this scalar, as daa_dt uses it

  ! other variables (needs to be consistent with reset list below)
  integer :: i_ssm=0,i_ugradpm=0

  contains

!***********************************************************************
    subroutine register_ent()
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
           "$Id: noentropy.f90,v 1.33 2002-11-12 16:07:59 dobler Exp $")
!
    endsubroutine register_ent
!***********************************************************************
    subroutine init_ent(f,xx,yy,zz)
!
!  initialise entropy; called from start.f90
!  28-mar-02/axel: dummy routine, adapted from entropy.f of 6-nov-01.
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      if(ip==1) print*,f,xx,yy,zz  !(to remove compiler warnings)
    endsubroutine init_ent
!***********************************************************************
    subroutine dss_dt(f,df,uu,glnrho,rho1,lnrho,cs2,TT1)
!
!  28-mar-02/axel: dummy routine, adapted from entropy.f of 6-nov-01.
!  19-may-02/axel: added isothermal pressure gradient
!   9-jun-02/axel: pressure gradient added to du/dt already here
!
      use Density
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: uu,glnrho
      real, dimension (nx) :: lnrho,rho1,cs2,TT1,uglnrho
      integer :: j,ju
!
      intent(in) :: f,uu,glnrho,rho1
      intent(out) :: cs2,TT1  !(df is dummy)
!
!  sound speed squared and inverse temperature
!
      TT1=0.
      if (gamma==1.) then
        cs2=cs20
      else
        cs2=cs20*exp(gamma1*lnrho)
      endif
      if (lfirst.and.ldt) maxadvec2=amax1(maxadvec2,cs2)
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
        if (i_ugradpm/=0) then
          call dot_mn(uu,glnrho,uglnrho)
          call sum_mn_name(cs2*uglnrho,i_ugradpm)
        endif
      endif
!
      if(ip==1) print*,f,df,uu,rho1  !(compiler)
    endsubroutine dss_dt
!***********************************************************************
    subroutine ss_run_hook()
!
!  called by run.f90 after reading parameters, but before the time loop
!
!
    endsubroutine ss_run_hook
!***********************************************************************
    subroutine rprint_entropy(lreset)
!
!  reads and registers print parameters relevant to entropy
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Cdata
      use Sub
!
      integer :: iname
      logical :: lreset
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_ssm=0; i_ugradpm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ugradpm',i_ugradpm)
      enddo
!
!  write column where which magnetic variable is stored
!
      write(3,*) 'i_ssm=',i_ssm
      write(3,*) 'i_ugradpm=',i_ugradpm
      write(3,*) 'nname=',nname
      write(3,*) 'ient=',ient
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
    subroutine bc_ss_flux(f,topbot)
!
!  dummy routine for entropy boundary condition
!
!  12-jun-2002/axel: coded
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
    endsubroutine bc_ss_flux
!***********************************************************************
    subroutine bc_ss_temp_old(f,topbot)
!
!  dummy routine for entropy boundary condition
!
!  12-jun-2002/axel: coded
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
    endsubroutine bc_ss_temp_old
!***********************************************************************
    subroutine bc_ss_temp_x(f,topbot)
!
!  dummy routine for entropy boundary condition
!
!  12-jun-2002/axel: coded
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
    endsubroutine bc_ss_temp_x
!***********************************************************************
    subroutine bc_ss_temp_y(f,topbot)
!
!  dummy routine for entropy boundary condition
!
!  12-jun-2002/axel: coded
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
    endsubroutine bc_ss_temp_y
!***********************************************************************
    subroutine bc_ss_temp_z(f,topbot)
!
!  dummy routine for entropy boundary condition
!
!  12-jun-2002/axel: coded
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
    endsubroutine bc_ss_temp_z
!***********************************************************************
    subroutine bc_ss_stemp_x(f,topbot)
!
!  dummy routine for entropy boundary condition
!
!  12-jun-2002/axel: coded
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
    endsubroutine bc_ss_stemp_x
!***********************************************************************
    subroutine bc_ss_stemp_y(f,topbot)
!
!  dummy routine for entropy boundary condition
!
!  12-jun-2002/axel: coded
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
    endsubroutine bc_ss_stemp_y
!***********************************************************************
    subroutine bc_ss_stemp_z(f,topbot)
!
!  dummy routine for entropy boundary condition
!
!  12-jun-2002/axel: coded
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
    endsubroutine bc_ss_stemp_z
!***********************************************************************
    subroutine bc_ss_energy(f,topbot)
!
!  dummy routine for entropy boundary condition
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
    endsubroutine bc_ss_energy
!***********************************************************************

endmodule Entropy
