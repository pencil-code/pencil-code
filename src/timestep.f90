module Timestep

  implicit none

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
      use Forcing
      use Cdata
      use Equ
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (3) :: alpha,beta,dt_beta
      real :: ds
      integer :: i,j
!
!HPF$ ALIGN (:,:,:,*) WITH tmpl :: f,df
!HPF$ SHADOW(0,3,3,0) :: f,df
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
      do i=1,itorder
        if (i==1) then
          lfirst=.true.
          df=0.
          ds=0.
        else
          lfirst=.false.
          df=alpha(i)*df  !(this could be subsumed into pde, but could be dangerous!)
          ds=alpha(i)*ds
        endif
        call pde(f,df)
        ds=ds+1.
!
!  if we are in the first step we need to calculate timestep
!  This can only be done using UUmax which was calculated in pde.
!  This only be done if we are on the root processor.
!  Then need to broadcast dt to all processors.
!
        if (lfirst.and.lroot) then
          if (ldt) dt = cdt*dxmin/UUmax
          if (ip<7) print*,'dt,cdt,dx,dy,dz,UUmax=',dt,cdt,dx,dy,dz,UUmax
        endif
        call mpibcast_real(dt,1)
        dt_beta=dt*beta
        if (ip<=6) print*,'TIMESTEP: iproc,dt=',iproc,dt  !(all have same dt?)
!
!  do this loop in slice, for cache efficiency
!
        do j=1,mvar
        do n=n1,n2
        do m=m1,m2
          f(l1:l2,m,n,j)=f(l1:l2,m,n,j)+dt_beta(i)*df(l1:l2,m,n,j)
        enddo
        enddo
        enddo
        t=t+dt_beta(i)*ds
      enddo
!
    endsubroutine rk_2n
!***********************************************************************

endmodule Timestep

