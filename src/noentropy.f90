module Entropy

  !
  ! isothermal case; almost nothing to do
  !

  use Cparam

  implicit none

  real, dimension (nx) :: cs2,TT1 ! Can't make this scalar, as daa_dt uses it

  ! input parameters
  namelist /entropy_init_pars/ &
       cs0,gamma,rho0

  ! run parameters
  namelist /entropy_run_pars/ &
       cs0

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
           "$RCSfile: noentropy.f90,v $", &
           "$Revision: 1.13 $", &
           "$Date: 2002-05-31 20:43:45 $")
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
      cs2 = cs20                ! (Really needed?)
      if(ip==1) print*,f,xx,yy,zz  !(to remove compiler warnings)
!
    endsubroutine init_ent
!***********************************************************************
    subroutine dss_dt(f,df,uu,uij,divu,rho1,glnrho,gpprho,cs2,TT1,chi)
!
!  28-mar-02/axel: dummy routine, adapted from entropy.f of 6-nov-01.
!  19-may-02/axel: added isothermal pressure gradient
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3,3) :: uij
      real, dimension (nx,3) :: uu,glnrho,gpprho
      real, dimension (nx) :: divu,rho1,chi
      real, dimension (nx) :: cs2,TT1
      logical, save :: first=.true.
!
      intent(in) :: f,uu,uij,divu,rho1,glnrho
      intent(out) :: gpprho,cs2,TT1,chi !(df is dummy)
!
      if (first) then
        TT1 = 0.
        chi = 0.
        cs2 = cs20
        if (gamma /= 1) then
          if (lroot) print*, 'Noentropy, thus resetting gamma to 1'
          gamma = 1
        endif
        first=.false.
      endif
      gpprho=cs20*glnrho
      if(ip==1) print*,f,df,uu,uij,divu,rho1,glnrho,gpprho  !(to remove compiler warnings)
!
    endsubroutine dss_dt
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
