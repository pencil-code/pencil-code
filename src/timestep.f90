! $Id: timestep.f90,v 1.26 2004-07-05 20:50:10 theine Exp $

module Timestep

  use Cparam
  use Cdata

  implicit none

!
!  border_prof_[x-z] could be of size n[x-z], but having the same
!  length as f() (in the given dimension) gives somehow more natural code.
!
  real, dimension(mx) :: border_prof_x=1.0
  real, dimension(my) :: border_prof_y=1.0
  real, dimension(mz) :: border_prof_z=1.0

  contains

!***********************************************************************
    subroutine rk_2n(f,df)
!
!  Runge Kutta advance, accurate to order itorder
!  At the moment, itorder can be 1, 2, or 3.
!
!   2-apr-01/axel: coded
!  14-sep-01/axel: moved itorder to cdata
!
      use Mpicomm
      use Cdata
      use Equ
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: ds
      integer :: j
!
!  coefficients for up to order 3
!
      if (itorder==1) then
        alpha=(/ 0., 0., 0. /)
        beta =(/ 1., 0., 0. /)
      elseif (itorder==2) then
        alpha=(/ 0., -1./2., 0. /)
        beta=(/ 1./2.,  1.,  0. /)
      elseif (itorder==3) then
        !alpha=(/0., -2./3., -1./)
        !beta=(/1./3., 1., 1./2./)
        !  use coefficients of Williamson (1980)
        alpha=(/  0. ,  -5./9., -153./128. /)
        beta=(/ 1./3., 15./16.,    8./15.  /)
      else
        if (lroot) print*,'Not implemented: itorder=',itorder
        call mpifinalize
      endif
!
!  dt_beta may be needed in other modules (like Dustdensity) for fixed dt
!
      if (.not. ldt) dt_beta=dt*beta
!
      do itsub=1,itorder
        if (itsub==1) then
          lfirst=.true.
          df=0.
          ds=0.
        else
          lfirst=.false.
          df=alpha(itsub)*df  !(could be subsumed into pde, but could be dangerous!)
          ds=alpha(itsub)*ds
        endif

        call pde(f,df)
        ds=ds+1.
!
!  If we are in the first step we need to calculate timestep dt.
!  This is done here because it uses UUmax which was calculated in pde.
!  Only do it on the root processor, then broadcast dt to all others.
!
        if (lroot.and.lfirst.and.ldt) dt=1.0/maxval(dt1_max)

        if (lfirst) call mpibcast_real(dt,1)
        if (ldt) dt_beta=dt*beta
        if (ip<=6) print*,'TIMESTEP: iproc,dt=',iproc,dt  !(all have same dt?)
!
!  do this loop in pencils, for cache efficiency
!
        do j=1,mvar
        do n=n1,n2
        do m=m1,m2
          f(l1:l2,m,n,j)=f(l1:l2,m,n,j)+dt_beta(itsub)*df(l1:l2,m,n,j) &
                        *border_prof_x(l1:l2)*border_prof_y(m)*border_prof_z(n)
        enddo
        enddo
        enddo
        t=t+dt_beta(itsub)*ds
      enddo
!
    endsubroutine rk_2n
!***********************************************************************
    subroutine border_profiles()
!
!  Position-dependent quenching factor that multiplies rhs of pde
!  by a factor that goes gradually to zero near the boundaries.
!  border_frac is a 3-D array, separately for all three directions.
!  border_frac=1 would affect everything between center and border.
!
      use Cdata

      real, dimension(nx) :: xi
      real, dimension(ny) :: eta
      real, dimension(nz) :: zeta
      real :: border_width,lborder,uborder
!
!  x-direction
!
      border_prof_x(l1:l2)=1

      if ((border_frac_x(1)>0) .and. (.not. lperi(1))) then
        border_width=border_frac_x(1)*Lxyz(1)/2
        lborder=xyz0(1)+border_width
        xi=1-max(lborder-x(l1:l2),0.0)/border_width
        border_prof_x(l1:l2)=min(border_prof_x(l1:l2),xi**2*(3-2*xi))
      endif

      if ((border_frac_x(2)>0) .and. (.not. lperi(1))) then
        border_width=border_frac_x(2)*Lxyz(1)/2
        uborder=xyz1(1)-border_width
        xi=1-max(x(l1:l2)-uborder,0.0)/border_width
        border_prof_x(l1:l2)=min(border_prof_x(l1:l2),xi**2*(3-2*xi))
      endif
!
!  y-direction
!
      border_prof_y(m1:m2)=1

      if ((border_frac_y(1)>0) .and. (.not. lperi(2))) then
        border_width=border_frac_y(1)*Lxyz(2)/2
        lborder=xyz0(2)+border_width
        eta=1-max(lborder-y(m1:m2),0.0)/border_width
        border_prof_y(m1:m2)=min(border_prof_y(m1:m2),eta**2*(3-2*eta))
      endif

      if ((border_frac_y(2)>0) .and. (.not. lperi(2))) then
        border_width=border_frac_y(2)*Lxyz(2)/2
        uborder=xyz1(2)-border_width
        eta=1-max(y(m1:m2)-uborder,0.0)/border_width
        border_prof_y(m1:m2)=min(border_prof_y(m1:m2),eta**2*(3-2*eta))
      endif
!
!  z-direction
!
      border_prof_z(n1:n2)=1

      if ((border_frac_z(1)>0) .and. (.not. lperi(3))) then
        border_width=border_frac_z(1)*Lxyz(3)/2
        lborder=xyz0(3)+border_width
        zeta=1-max(lborder-z(n1:n2),0.0)/border_width
        border_prof_z(n1:n2)=min(border_prof_z(n1:n2),zeta**2*(3-2*zeta))
      endif

      if ((border_frac_z(2)>0) .and. (.not. lperi(3))) then
        border_width=border_frac_z(2)*Lxyz(3)/2
        uborder=xyz1(3)-border_width
        zeta=1-max(z(n1:n2)-uborder,0.0)/border_width
        border_prof_z(n1:n2)=min(border_prof_z(n1:n2),zeta**2*(3-2*zeta))
      endif
!
    endsubroutine border_profiles
!***********************************************************************

endmodule Timestep
