! $Id: initcond.f90,v 1.25 2003-03-02 17:54:40 brandenb Exp $ 

module Initcond 
 
!  This module contains code used by the corresponding physics
!  modules to set up various initial conditions (stratitication,
!  perturbations, and other structures). This module is not used
!  during run time (although it is used by the physics modules that
!  are used both during run time and for the initial condition).

  use Cdata
  use General

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
        f(:,:,:,j)=f(:,:,:,j)+ampl*(spread(spread(cos(kx1*x),2,my),3,mz)&
                                   *spread(spread(cos(ky1*y),1,mx),3,mz)&
                                   *spread(spread(cos(kz1*z),1,mx),2,my))
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
    subroutine jump(f,i,fleft,fright,width,dir)
!
!  jump
!
!  19-sep-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx) :: prof
      real :: fleft,fright,width
      character(len=*) :: dir
!
!  jump; check direction
!
      select case(dir)
!
      case('x')
        prof=fleft+(fright-fleft)*.5*(1.+tanh(x/width))
        f(:,:,:,i)=f(:,:,:,i)+spread(spread(prof,2,my),3,mz)
!
      case default
        print*,'jump: no default value'
!
      endselect
!
    endsubroutine jump
!***********************************************************************
    subroutine bjump(f,i,fleft,fright,width,dir)
!
!  jump in B-field (in terms of magnetic vector potential)
!
!   9-oct-02/wolf+axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx) :: prof,alog_cosh_xwidth
      real :: fleft,fright,width
      character(len=*) :: dir
!
!  jump; check direction
!
      select case(dir)
!
!  use correct signs when calling this routine
!  Ay=+int Bz dx
!  Az=-int By dx
!
!  alog(cosh(x/width)) = 
!
      case('x')
        alog_cosh_xwidth=abs(x/width)+alog(.5*(1.+exp(-2*abs(x/width))))
        prof=.5*(fright+fleft)*x &
            +.5*(fright-fleft)*width*alog_cosh_xwidth
        f(:,:,:,i)=f(:,:,:,i)-spread(spread(prof,2,my),3,mz)
      case default
        print*,'jump: no default value'
!
      endselect
!
    endsubroutine bjump
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
    subroutine planet(ampl,f,xx,yy,zz,eps,radius,gamma,ztop)
!
!  Ellipsoidal planet solution (Goldreich, Narayan, Goodman 1987)
!
!   6-jul-02/axel: coded
!  22-feb032/axel: fixed 3-D background solution for enthalpy
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz,rr2,hh
      real :: ampl,sigma2,sigma,delta2,delta,eps,radius,hmax
      real :: gamma,eps2,radius2
      real :: cs20=1.,gamma1,ztop
      integer :: l,m,n
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
!  calculate hh
!  add continuous vertical stratification to horizontal planet solution
!
      rr2=xx**2+eps2*yy**2
      hh=+.5*delta2*Omega**2*amax1(radius2-rr2,0.)+.5*Omega**2*(ztop**2-zz**2)
      print*,'planet: ztop=',ztop,minval(hh)
!
!  limit dynamical to 1% of maximum value (can always add a constant)
!
      hh=amax1(hh,0.)
      hmax=maxval(hh)
      hh=hh+ampl*hmax
!
      rr2=xx**2+eps2*yy**2+zz**2/delta2
      do n=n1,n2
        do m=m1,m2
          do l=l1,l2
            if (rr2(l,m,n)<radius2) then
              f(l,m,n,iux)=   eps2*sigma *Omega*yy(l,m,n)
              f(l,m,n,iuy)=(qshear-sigma)*Omega*xx(l,m,n)
            else
              f(l,m,n,iux)=0.
              f(l,m,n,iuy)=0.
            endif
          enddo
        enddo
      enddo
!
!  calculate density, depending on what gamma is
!
      print*,'planet: hmin=',minval(hh)
      if(gamma<=1.) print*,'must have gamma>1 for planet solution'
      if(gamma==1.) then
        f(:,:,:,ilnrho)=hh/cs20
        print*,'planet solution for gamma=1'
      else
        gamma1=gamma-1.
        f(:,:,:,ilnrho)=alog(gamma1*hh/cs20)/gamma1
        print*,'planet solution for gamma=',gamma
      endif
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
    subroutine htube(ampl,f,i1,i2,xx,yy,zz,radius,epsilon_nonaxi)
!
!  Horizontal flux tube (for vector potential, or passive scalar)
!
!   7-jun-02/axel+vladimir: coded
!  11-sep-02/axel: allowed for scalar field (if i1=i2)
!
      integer :: i1,i2
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,xx,yy,zz,modulate
      real :: ampl,radius,epsilon_nonaxi,ky
!
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'set variable to zero; i1,i2=',i1,i2
      else
        ky=2*pi/Ly
        if(lroot) then
          print*,'implement y-dependent flux tube in xz-plane; i1,i2=',i1,i2
          print*,'radius,epsilon_nonaxi=',radius,epsilon_nonaxi
        endif
        modulate=1.+epsilon_nonaxi*sin(ky*yy)
!
! completely quenched "gaussian"
!
        tmp=.5*ampl/modulate*exp(-(xx**2+zz**2)/(max((radius*modulate)**2-xx**2-zz**2,1e-6)))
!
!  check whether vector or scalar
!
        if(i1==i2) then
          if(lroot) print*,'htube: set scalar'
          f(:,:,:,i1)=tmp
        elseif(i1+2==i2) then
          if(lroot) print*,'htube: set vector'
          f(:,:,:,i1 )=+zz*tmp
          f(:,:,:,i1+1)=0.
          f(:,:,:,i1+2)=-xx*tmp
        else
          if(lroot) print*,'htube: bad value of i2=',i2
        endif
      endif
!
    endsubroutine htube
!***********************************************************************
    subroutine htube2(ampl,f,i1,i2,xx,yy,zz,radius,epsilon_nonaxi)
!
!  Horizontal flux tube (for vector potential, or passive scalar)
!
!   7-jun-02/axel+vladimir: coded
!  11-sep-02/axel: allowed for scalar field (if i1=i2)
!
      integer :: i1,i2
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,xx,yy,zz,modulate
      real :: ampl,radius,epsilon_nonaxi,ky
!
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'set variable to zero; i1,i2=',i1,i2
      else
        ky=2*pi/Ly
        if(lroot) then
          print*,'implement y-dependent flux tube in xz-plane; i1,i2=',i1,i2
          print*,'radius,epsilon_nonaxi=',radius,epsilon_nonaxi
        endif
!
!  constant, when epsilon_nonaxi; otherwise modulation about zero
!
        if(epsilon_nonaxi==0) then
          modulate=1.
        else
          modulate=epsilon_nonaxi*sin(ky*yy)
        endif
!
! completely quenched "gaussian"
!
        tmp=.5*ampl*modulate*exp(-(xx**2+zz**2)/radius**2)
!
!  check whether vector or scalar
!
        if(i1==i2) then
          if(lroot) print*,'htube: set scalar'
          f(:,:,:,i1)=tmp
        elseif(i1+2==i2) then
          if(lroot) print*,'htube: set vector'
          f(:,:,:,i1 )=+zz*tmp
          f(:,:,:,i1+1)=0.
          f(:,:,:,i1+2)=-xx*tmp
        else
          if(lroot) print*,'htube: bad value of i2=',i2
        endif
      endif
!
    endsubroutine htube2
!***********************************************************************
    subroutine magsupport(ampl,f,zz,gravz,cs0,rho0)
!
!  magnetically supported horizontal flux layer
!  (for aa):  By^2 = By0^2 * exp(-z/H),
!  where H=2*cs20/abs(gravz) and ampl=cs0*sqrt(2*rho0)
!  should be used when invoking this routine.
!  Here, ampl=pmag/pgas.
!
!   7-dec-02/axel: coded
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: zz
      real :: ampl,H,A0,gravz,cs0,rho0,lnrho0
!
      if (ampl==0) then
        if (lroot) print*,'magsupport: do nothing'
      else
        lnrho0=alog(rho0)
        H=(1+ampl)*cs0**2/abs(gravz)
        A0=-2*H*ampl*cs0*sqrt(2*rho0)
        if (lroot) print*,'magsupport; H,A0=',H,A0
        f(:,:,:,iaa)=A0*exp(-.5*zz/H)
        f(:,:,:,ilnrho)=lnrho0-zz/H
      endif
!
    endsubroutine magsupport
!***********************************************************************
    subroutine hfluxlayer(ampl,f,i,xx,yy,zz,zflayer,width)
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
    endsubroutine hfluxlayer
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
            call random_number_wrapper(r)
            call random_number_wrapper(p)
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
      call random_number_wrapper(r)
      call random_number_wrapper(p)
      tmp=sqrt(-2*alog(r))*sin(2*pi*p)
      f(:,:,:,i)=ampl*tmp
      print*,'set gaussian noise: variable i=',i
!
    endsubroutine gaunoise_scal
!***********************************************************************
    subroutine trilinear(ampl,f,ivar,xx,yy,zz)
!
!  Produce a profile that is linear in any non-periodic direction, but
!  periodic in periodic ones (for testing purposes).
!
!  5-nov-02/wolf: coded
! 23-nov-02/axel: included scaling factor ampl, corrected lperi argument
!
      integer :: ivar
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,xx,yy,zz
      real :: ampl
!
      if (lroot) print*, 'uu: trilinear in ', ivar
!
!  x direction
!
      if (lperi(1)) then
        tmp = sin(2*pi/Lx*(xx-xyz0(1)-0.25*Lxyz(1)))
      else
        tmp = xx
      endif
!
!  y direction
!
      if (lperi(2)) then
        tmp = tmp + 10*sin(2*pi/Ly*(yy-xyz0(2)-0.25*Lxyz(2)))
      else
        tmp = tmp + 10*yy
      endif
!
!  z direction
!
      if (lperi(3)) then
        tmp = tmp + 100*sin(2*pi/Lz*(zz-xyz0(3)-0.25*Lxyz(3)))
      else
        tmp = tmp + 100*zz
      endif
!
      f(:,:,:,ivar) = ampl*tmp
!
    endsubroutine trilinear
!***********************************************************************
    subroutine cos_cos_sin(ampl,f,ivar,xx,yy,zz)
!
!  Produce a profile that is linear in any non-periodic direction, but
!  periodic in periodic ones (for testing purposes).
!
!  7-dec-02/axel: coded
!
      integer :: ivar
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: ampl,kx,ky,kz
!
      if (lroot) print*, 'uu: trilinear in ', ivar
!
      kx=2*pi/Lx*3
      ky=2*pi/Ly*3
      kz=pi/Lz
      f(:,:,:,ivar) = ampl*cos(kx*xx)*cos(ky*yy)*sin(kz*zz)
!
    endsubroutine cos_cos_sin
!***********************************************************************
    subroutine tor_pert(ampl,f,ivar,xx,yy,zz)
!
!  Produce a profile that is periodic in the y- and z-directions.
!  For testing the Balbus-Hawley instability of a toroidal magnetic field
!
!  12-feb-03/ulf: coded
!
      integer :: ivar
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: ampl,ky,kz
!
      if (lroot) print*, 'uu: trilinear in ', ivar
!
      ky=2*pi/Ly
      kz=2.*pi/Lz
      f(:,:,:,ivar) = ampl*cos(ky*yy)*cos(kz*zz)
!
      print*,'xx(1,1,1)=',xx(1,1,1) !(to keep compiler quiet)
    endsubroutine tor_pert
!***********************************************************************

endmodule Initcond
