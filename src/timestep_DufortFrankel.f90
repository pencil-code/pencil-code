! $Id$
!
!  DuFort-Frankel time stepping routine. As is well known, this proceedure
!  has disadvantages, but it is well suited to deal with the time step constrain
!  near the center in spherical polar coordinates.
!
module Timestep
!
  use Cparam
  use Cdata
!
  implicit none
!
  private
!
  public :: time_step
!
  contains
!***********************************************************************
    subroutine time_step(f,df,p)
!
!  Runge Kutta advance, accurate to order itorder.
!  At the moment, itorder can be 1, 2, or 3.
!
!   2-apr-01/axel: coded
!  14-sep-01/axel: moved itorder to cdata
!
      use Mpicomm
      use Cdata
      use Equ
      use Particles_main
      use BorderProfiles, only: border_quenching
      use Interstellar, only: calc_snr_damp_int
      use Shear, only: advance_shear
      use Special, only: special_after_timestep
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df,fold
      real, dimension (nx) :: diffarr,d2o,d2n,d2s,d2e,d2w
      real, dimension (nx) :: cdt1,cdt2,fold_tmp
      type (pencil_case) :: p
      real :: ds
      real :: dt1, dt1_local, dt1_last=0.0
      integer :: j
!
!  Set up df and ds for each time sub.
!
      lfirst=.true.
      df=0.
      ds=1.
!
!  make sure fold exists at start-up
!
      if (it==1) fold=f
!
!  Change df according to the chosen physics modules.
!
      call pde(f,df,p)
!
!  If we are in the first time substep we need to calculate timestep dt.
!  Only do it on the root processor, then broadcast dt to all others.
!
      if (ldt) then
        dt1_local=maxval(dt1_max(1:nx))
        !Timestep growth limiter
        if (real(ddt) > 0.) dt1_local=max(dt1_local,dt1_last)
        call mpireduce_max(dt1_local,dt1)
        if (lroot) dt=1.0/dt1
        !Timestep growth limiter
        if (ddt/=0.) dt1_last=dt1_local/ddt
        call mpibcast_real(dt,1)
      endif
      if (ip<=6) print*,'TIMESTEP: iproc,dt=',iproc,dt  !(all have same dt?)
!
!  Time evolution of grid variables.
!
diffarr=1.
!
      do m=m1,m2
        if (coord_system=='cartesian') then
          d2o=-2.*(dx_1(l1:l2)**2+dy_1(m)**2)
          d2n=dx_1(l1:l2)**2
          d2s=dx_1(l1:l2)**2
          d2e=dy_1(m)**2
          d2w=dy_1(m)**2
        elseif (coord_system=='spherical') then
          d2o=-2.*dx_1(l1:l2)**2-2.*r2_mn*dy_1(m)**2-r2_mn*sin2th(m)
          d2n=dx_1(l1:l2)**2+r1_mn*dx_1(l1:l2)
          d2s=dx_1(l1:l2)**2-r1_mn*dx_1(l1:l2)
          d2e=dy_1(m)**2+.5*tanth(m)*dy_1(m)*r2_mn
          d2w=dy_1(m)**2-.5*tanth(m)*dy_1(m)*r2_mn
        endif
!
!  copy current f to fold_tmp, and set it to fold after it has been used
!
        do j=1,mvar; do n=n1,n2;
          if (lborder_profiles) call border_quenching(f,df,j,dt_beta_ts(itsub))
          cdt1=.5*(1./dt+diffarr*d2o)
          cdt2=2./(1./dt-diffarr*d2o)
          fold_tmp=f(l1:l2,m,n,j)
          f(l1:l2,m,n,j)=(df(l1:l2,m,n,j)+diffarr*( &
             d2n*f(l1+1:l2+1,m,n,j)+d2s*f(l1-1:l2-1,m,n,j) &
            +d2e*f(l1:l2,m+1,n,j)+d2w*f(l1:l2,m-1,n,j) &
            )+cdt1*fold(l1:l2,m,n,j))*cdt2
          fold(l1:l2,m,n,j)=fold_tmp
        enddo; enddo
      enddo
!
      if (lspecial) &
          call special_after_timestep(f,df,dt_beta_ts(itsub)*ds)
!
!  Increase time.
!
      t=t+dt*ds
!
    endsubroutine time_step
!***********************************************************************
endmodule Timestep
