module Entropy

  !
  ! isothermal case; almost nothing to do
  !

  use Cparam

  implicit none

  real, dimension (nx) :: cs2,TT1 ! Can't make this scalar, as daa_dt uses it
  integer :: ient

  integer :: dummy              ! We cannot define empty namelists
  namelist /entropy_init_pars/ dummy
  namelist /entropy_run_pars/  dummy

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
           "$Revision: 1.6 $", &
           "$Date: 2002-05-11 12:18:48 $")
!
    endsubroutine register_ent
!***********************************************************************
    subroutine init_ent(f,init,ampl,xx,yy,zz)
!
!  initialise entropy; called from start.f90
!  28-mar-02/axel: dummy routine, adapted from entropy.f of 6-nov-01.
!
      use Cdata
      use sub, only: step
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,r,p,xx,yy,zz
      real, dimension (mz) :: stp
      real :: ampl,beta1,cs2int,ssint
      integer :: init
!
      cs2 = cs20                ! (Really needed?)
!
    endsubroutine init_ent
!***********************************************************************
    subroutine dss_dt(f,df,uu,uij,divu,rho1,glnrho,gpprho,cs2,TT1,chi)
!
!  28-mar-02/axel: dummy routine, adapted from entropy.f of 6-nov-01.
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
      if (first) then
        cs2 = cs20
        TT1 = 1./cs20           ! Assumes R/mu = 1, ehich is probably OK,
                                ! as the usual c_p = 1 does not work in the
                                ! isothermal case
        if (gamma /= 1) print*, 'Noentropy, thus resetting gamma to 1'
        gamma = 1
        first=.false.
      else
        print*,"Noentropy dss_dt: This can't happen"
      endif
!
    endsubroutine dss_dt
!***********************************************************************
    subroutine heatcond(x,y,z,hcond)
!
!  calculate the heat conductivity lambda
!  NB: if you modify this profile, you *must* adapt gradloghcond below.
!  23-jan-2002/wolf: coded
!  28-mar-02/axel: dummy routine, adapted from entropy.f of 6-nov-01.
!
      use Cdata, only: nx,lgravz,lgravr,z0,z1,z2,ztop, &
           hcond0,hcond1,hcond2,whcond
      use Sub, only: step
!
      real, dimension (nx) :: x,y,z
      real, dimension (nx) :: hcond
!
    endsubroutine heatcond
!***********************************************************************
    subroutine gradloghcond(x,y,z,glhc)
!
!  calculate grad(log lambda), where lambda is the heat conductivity
!  NB: *Must* be in sync with heatcond() above.
!  23-jan-2002/wolf: coded
!  28-mar-02/axel: dummy routine, adapted from entropy.f of 6-nov-01.
!
      use Cdata, only: nx,lgravz,lgravr,z0,z1,z2,ztop, &
           hcond0,hcond1,hcond2,whcond
      use Sub, only: der_step
!
      real, dimension (nx) :: x,y,z
      real, dimension (nx,3) :: glhc
!
    endsubroutine gradloghcond
!***********************************************************************

endmodule Entropy
