module Entropy

  use Cparam

  implicit none

  integer :: ient

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
           "$Revision: 1.1 $", &
           "$Date: 2002-03-28 18:51:49 $")
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
use Mpicomm
use IO
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,r,p,xx,yy,zz
      real, dimension (mz) :: stp
      real :: ampl,beta1,cs2int,ssint
      integer :: init
!
    endsubroutine init_ent
!***********************************************************************
    subroutine dss_dt(f,df,uu,uij,divu,rho1,glnrho,gpprho,cs2,chi)
!
!  28-mar-02/axel: dummy routine, adapted from entropy.f of 6-nov-01.
!
      use Cdata
      use Mpicomm
      use Sub
      use Global
      use Slices
      use IO
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3,3) :: uij,sij
      real, dimension (nx,3) :: uu,glnrho,gpprho,gss,glnT,glnTlambda,glhc
      real, dimension (nx) :: divu,rho1,cs2,chi
      real, dimension (nx) :: x_mn,y_mn,z_mn,r_mn
      real, dimension (nx) :: ugss,thdiff,del2ss,del2lnrho,sij2,g2
      real, dimension (nx) :: ss,lnrho,TT1,lambda
      real, dimension (nx) :: heat,prof
      real :: ssref,z_prev=-1.23e20
      integer :: i,j
!
      save :: z_prev,lambda,glhc
!
      intent(in) :: f,uu,uij,divu,rho1,glnrho
      intent(out) :: df,gpprho,cs2,chi
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
