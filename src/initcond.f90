! $Id: initcond.f90,v 1.8 2002-08-15 19:09:06 nilshau Exp $ 

module Initcond 
 
!  This module contains code used by the corresponding physics
!  modules to set up various initial conditions (stratitication,
!  perturbations, and other structures). This module is not used
!  during run time (although it is used by the physics modules that
!  are used both during run time and for the initial condition).

  use Cdata

  implicit none

  interface gaunoise            ! Overload the `set_random' function
    module procedure gaunoise_vect
    module procedure gaunoise_scal
  endinterface

  contains

!***********************************************************************
    subroutine sinxsinz(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave
!
!  26-jul-02/axel: coded
!
      integer :: i,j
      real, dimension (mx,my,mz,mvar) :: f
      real,optional :: kx,ky,kz
      real :: ampl,kx1=pi/2.,ky1=0.,kz1=pi/2.
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  sinx(kx*x)*sin(kz*z)
!
      if (present(kx)) kx1=kx
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
      if (ampl==0) then
        if (lroot) print*,'ampl=0 in sinx*sinz wave; kx,kz=',kx1,kz1
      else
        if (lroot) print*,'sinx*sinz wave; ampl,kx,kz=',ampl,kx1,kz1
        j=i+1
        f(:,:,:,j)=f(:,:,:,j)+ampl*(spread(spread(sin(kx1*x),2,my),3,mz)&
                                   *spread(spread(sin(kz1*z),1,mx),2,my))
      endif
!
    endsubroutine sinxsinz
!***********************************************************************
    subroutine wave(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave
!
!   6-jul-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  set x-wave
!
      if (present(kx)) then
        k=kx
        if (ampl==0) then
          if (lroot) print*,'ampl=0 in wave; kx=',k
        else
          if (lroot) print*,'wave: kx,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(sin(k*x),2,my),3,mz)
        endif
      endif
!
!  set y-wave
!
      if (present(ky)) then
        k=ky
        if (ampl==0) then
          if (lroot) print*,'ampl=0 in wave; ky=',k
        else
          if (lroot) print*,'wave: ky,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(sin(k*y),1,mx),3,mz)
        endif
      endif
!
!  set z-wave
!
      if (present(kz)) then
        k=kz
        if (ampl==0) then
          if (lroot) print*,'ampl=0 in wave; kz=',k
        else
          if (lroot) print*,'wave: kz,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(sin(k*z),1,mx),2,my)
        endif
      endif
!
    endsubroutine wave
!***********************************************************************
    subroutine beltrami(ampl,f,i,kx,ky,kz)
!
!  Beltrami field (as initial condition)
!
!  19-jun-02/axel: coded
!   5-jul-02/axel: made additive (if called twice), kx,ky,kz are optional
!
      integer :: i,j
      real, dimension (mx,my,mz,mvar) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.,fac
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  set x-dependent Beltrami field
!
      if (present(kx)) then
        k=kx; if(k==0) print*,'k must not be zero!'; fac=sqrt(abs(ampl/k))
        if (ampl==0) then
          if (lroot) print*,'ampl=0 in Beltrami field; kx=',k
        elseif (ampl>0) then
          if (lroot) print*,'Beltrami field (pos-hel): kx,i=',k,i
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*x),2,my),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*x),2,my),3,mz)
        elseif (ampl<0) then
          if (lroot) print*,'Beltrami field (neg-hel): kx,i=',k,i
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*x),2,my),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*x),2,my),3,mz)
        endif
      endif
!
!  set y-dependent Beltrami field
!
      if (present(ky)) then
        k=ky; if(k==0) print*,'k must not be zero!'; fac=sqrt(abs(ampl/k))
        if (ampl==0) then
          if (lroot) print*,'ampl=0 in Beltrami field; ky=',k
        elseif (ampl>0) then
          if (lroot) print*,'Beltrami field (pos-hel): ky,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*y),1,mx),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*y),1,mx),3,mz)
        elseif (ampl<0) then
          if (lroot) print*,'Beltrami field (neg-hel): ky,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*y),1,mx),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*y),1,mx),3,mz)
        endif
      endif
!
!  set z-dependent Beltrami field
!
      if (present(kz)) then
        k=kz; if(k==0) print*,'k must not be zero!'; fac=sqrt(abs(ampl/k))
        if (ampl==0) then
          if (lroot) print*,'ampl=0 in Beltrami field; kz=',k
        elseif (ampl>0) then
          if (lroot) print*,'Beltrami field (pos-hel): kz,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*z),1,mx),2,my)
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*z),1,mx),2,my)
        elseif (ampl<0) then
          if (lroot) print*,'Beltrami field (neg-hel): kz,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*z),1,mx),2,my)
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*z),1,mx),2,my)
        endif
      endif
!
    endsubroutine beltrami
!***********************************************************************
    subroutine planet(ampl,f,xx,yy,zz,eps,radius,gamma)
!
!  Ellipsoidal planet solution (Goldreich, Narayan, Goodman 1987)
!
!   6-jul-02/axel: coded
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz,rr2,psi,hh,h0
      real :: ampl,sigma2,sigma,delta2,delta,eps,radius,hmax,hmin
      real :: gamma,eps2,rad2,radius2
      integer :: i,j
!
!  calculate sigma
!
      print*,'planet: qshear,eps=',qshear,eps
      eps2=eps**2
      radius2=radius**2
      sigma2=2*qshear/(1.-eps2)
      if (sigma2<0.) then
        print*,'sigma2<0 not allowed; choose another value of eps_planet'
      else
        sigma=sqrt(sigma2)
      endif
!
!  calculate delta
!
      delta2=(2.-sigma)*sigma
      print*,'planet: sigma,delta2=',sigma,delta2
      if (delta2<=0.) then
        print*,'delta2<=0 not allowed'
      else
        delta=sqrt(delta2)
      endif
!
!  calculate psi, hh, and h0
!
!      h0=ampl**(gamma-1.)-zz**2/delta2
!      h0=0.
      rr2=xx**2+eps2*yy**2  !+zz**2/Rz**2
      hh=+.5*delta2*Omega**2*(radius2-rr2)
!
!  limit dynamical range to 1:100
!
      hmax=maxval(hh)
      hmin=ampl**(gamma-1)*hmax
      hh=amax1(hh,hmin)
!
      if (gamma<=1.) print*,'must have gamma>1 for planet solution'
!
      do i=l1,l2
        do j=m1,m2
          rad2=x(i)**2+eps2*y(j)**2
!          if (rad2<radius2) then
          if (hh(i,j,4)>hmin) then
            f(i,j,:,iux)=   eps2*sigma *Omega*yy(i,j,:)
            f(i,j,:,iuy)=(qshear-sigma)*Omega*xx(i,j,:)
          else
            f(i,j,:,iux)=0.
            f(i,j,:,iuy)=0.
          endif
        enddo
      enddo
      f(:,:,:,ilnrho)=alog(hh)/(gamma-1.)
!
    endsubroutine planet
!***********************************************************************
    subroutine crazy(ampl,f,i)
!
!  A crazy initial condition
!  (was useful to initialize all points with finite values)
!
!  19-may-02/axel: coded
!
      integer :: i,j
      real, dimension (mx,my,mz,mvar) :: f
      real :: ampl
!
      if (lroot) print*, 'sinusoidal magnetic field: for debugging purposes'
      j=j; f(:,:,:,j)=f(:,:,:,j)+ampl*&
        spread(spread(sin(2*x),2,my),3,mz)*&
        spread(spread(sin(3*y),1,mx),3,mz)*&
        spread(spread(cos(1*z),1,mx),2,my)
      j=i+1; f(:,:,:,j)=f(:,:,:,j)+ampl*&
        spread(spread(sin(5*x),2,my),3,mz)*&
        spread(spread(sin(1*y),1,mx),3,mz)*&
        spread(spread(cos(2*z),1,mx),2,my)
      j=j+2; f(:,:,:,j)=f(:,:,:,j)+ampl*&
        spread(spread(sin(3*x),2,my),3,mz)*&
        spread(spread(sin(4*y),1,mx),3,mz)*&
        spread(spread(cos(2*z),1,mx),2,my)
!
    endsubroutine crazy
!***********************************************************************
    subroutine htube(ampl,f,i,xx,yy,zz,radius,epsilon_nonaxi)
!
!  Horizontal flux tube (for vector potential)
!
!   7-jun-02/axel+vladimir: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,xx,yy,zz,modulate
      real :: ampl,radius,epsilon_nonaxi,ky
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'set variable to zero; i=',i
      else
        ky=2*pi/Ly
        print*,'implement y-dependent flux tube in xz-plane; i=',i
        print*,'radius,epsilon_nonaxi=',radius,epsilon_nonaxi
        modulate=1.+epsilon_nonaxi*sin(ky*yy)
! completely quenched "gaussian"
        tmp=.5*ampl/modulate*exp(-(xx**2+zz**2)/(max((radius*modulate)**2-xx**2-zz**2,1e-6)))
        if ((ip<=8).and.lroot) print*,'horizontal flux tube: i=',i
        f(:,:,:,i  )=+zz*tmp
        f(:,:,:,i+1)=0.
        f(:,:,:,i+2)=-xx*tmp
      endif
!
    endsubroutine htube
!***********************************************************************
    subroutine hlayer(ampl,f,i,xx,yy,zz,zflayer,width)
!
!  Horizontal flux layer (for vector potential)
!
!  19-jun-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: ampl,zflayer,width
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'HLAYER: set variable to zero; i=',i
      else
        if (lroot) print*,'horizontal flux layer; i=',i
        if ((ip<=16).and.lroot) print*,'ampl,width=',ampl,width
        f(:,:,:,i  )=0.
        f(:,:,:,i+1)=ampl*tanh((zz-zflayer)/width)
        f(:,:,:,i+2)=0.
      endif
!
      if (ip==1) print*,xx,yy
    endsubroutine hlayer
!***********************************************************************
    subroutine uniform_x(ampl,f,i,xx,yy,zz)
!
!  Uniform B_x field (for vector potential)
!
!  19-jun-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: ampl
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'set variable to zero; i=',i
      else
        print*,'uniform x-field ; i=',i
        if ((ip<=16).and.lroot) print*,'ampl=',ampl
        f(:,:,:,i  )=0.
        f(:,:,:,i+1)=-ampl*zz
        f(:,:,:,i+2)=0.
      endif
!
      if (ip==1) print*,xx,yy
    endsubroutine uniform_x
!***********************************************************************
    subroutine uniform_y(ampl,f,i,xx,yy,zz)
!
!  Uniform B_y field (for vector potential)
!
!  27-jul-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: ampl
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'set variable to zero; i=',i
      else
        print*,'uniform x-field ; i=',i
        if ((ip<=16).and.lroot) print*,'ampl=',ampl
        f(:,:,:,i  )=ampl*zz
        f(:,:,:,i+1)=0.
        f(:,:,:,i+2)=0.
      endif
!
      if (ip==1) print*,xx,yy
    endsubroutine uniform_y
!***********************************************************************
    subroutine vfield(ampl,f,i,xx)
!
!  Vertical field, for potential field test
!
!  14-jun-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx
      real :: ampl,kx
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'set variable to zero; i=',i
      else
        kx=2*pi/Lx
        print*,'implement x-dependent vertical field'
        if ((ip<=8).and.lroot) print*,'x-dependent vertical field'
        f(:,:,:,i  )=0.
        f(:,:,:,i+1)=ampl*sin(kx*xx)
        f(:,:,:,i+2)=0.
      endif
!
    endsubroutine vfield
!***********************************************************************
    subroutine gaunoise_vect(ampl,f,i1,i2)
!
!  Write snapshot file of penciled vector data (for debugging).
!
!  23-may-02/axel: coded
!
      integer :: i,i1,i2
      real, dimension (mx,my,mz) :: r,p,tmp
      real, dimension (mx,my,mz,mvar) :: f
      real :: ampl
!
!  set gaussian random noise vector
!
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'set variable to zero; i1,i2=',i1,i2
      else
        if ((ip<=8).and.lroot) print*,'set_random_vect: i1,i2=',i1,i2
        do i=i1,i2
          if (modulo(i-i1,2)==0) then
            call random_number(r)
            call random_number(p)
            tmp=sqrt(-2*alog(r))*sin(2*pi*p)
          else
            tmp=sqrt(-2*alog(r))*cos(2*pi*p)
          endif
          !call smooth_3d(tmp,ismo)  !(may want to smooth)
          f(:,:,:,i)=ampl*tmp
          if (lroot) print*,'set gaussian noise: variable i=',i
        enddo
      endif
!
    endsubroutine gaunoise_vect
!***********************************************************************
    subroutine gaunoise_scal(ampl,f,i)
!
!  Write snapshot file of penciled vector data (for debugging).
!
!  23-may-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz) :: r,p,tmp
      real, dimension (mx,my,mz,mvar) :: f
      real :: ampl
!
!  set gaussian random noise vector
!
      if ((ip<=8).and.lroot) print*,'set_random_scal: i=',i
      call random_number(r)
      call random_number(p)
      tmp=sqrt(-2*alog(r))*sin(2*pi*p)
      !call smooth_3d(tmp,ismo)  !(may want to smooth)
      f(:,:,:,i)=ampl*tmp
      print*,'set gaussian noise: variable i=',i
!
    endsubroutine gaunoise_scal

endmodule Initcond
