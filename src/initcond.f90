! $Id: initcond.f90,v 1.90 2003-11-13 14:40:27 ajohan Exp $ 

module Initcond 
 
!  This module contains code used by the corresponding physics
!  modules to set up various initial conditions (stratitication,
!  perturbations, and other structures). This module is not used
!  during run time (although it is used by the physics modules that
!  are used both during run time and for the initial condition).

  use Cdata
  use General
  use Mpicomm

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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
        if (lroot) print*,'sinxsinz: ampl=0 in sinx*sinz wave; kx,kz=',kx1,kz1
      else
        if (lroot) print*,'sinxsinz: sinx*sinz wave; ampl,kx,kz=',ampl,kx1,kz1
        j=i+1
        f(:,:,:,j)=f(:,:,:,j)+ampl*(spread(spread(cos(kx1*x),2,my),3,mz)&
                                   *spread(spread(cos(ky1*y),1,mx),3,mz)&
                                   *spread(spread(cos(kz1*z),1,mx),2,my))
      endif
!
    endsubroutine sinxsinz
!***********************************************************************
    subroutine hat(ampl,f,i,width,kx,ky,kz)
!
!  hat bump
!
!   2-may-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar+maux) :: f
      real,optional :: kx,ky,kz
      real :: ampl,width,k=1.,width2,k2
!
!  prepare
!
      width2=width**2
!
!  set x-hat
!
      if (present(kx)) then
        k=kx
        k2=k**2
        if (ampl==0) then
          if (lroot) print*,'hat: ampl=0; kx=',k
        else
          if (lroot) print*,'hat: kx,i,ampl=',k,i,ampl
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(.5+.5*tanh(k2*(width2-x**2)),2,my),3,mz)
        endif
      endif
!
!  set y-hat
!
      if (present(ky)) then
        k=ky
        k2=k**2
        if (ampl==0) then
          if (lroot) print*,'hat: ampl=0; ky=',k
        else
          if (lroot) print*,'hat: ky,i,ampl=',k,i,ampl
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(.5+.5*tanh(k2*(width2-y**2)),1,mx),3,mz)
        endif
      endif
!
!  set z-hat
!
      if (present(kz)) then
        k=kz
        k2=k**2
        if (ampl==0) then
          if (lroot) print*,'hat: ampl=0; kz=',k
        else
          if (lroot) print*,'hat: kz,i,ampl=',k,i,ampl
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(.5+.5*tanh(k2*(width2-z**2)),1,mx),2,my)
        endif
      endif
!
    endsubroutine hat
!***********************************************************************
    subroutine gaussian(ampl,f,i,kx,ky,kz)
!
!  gaussian bump
!
!   2-may-03/axel: coded
!  20-sep-03/axel: added 1/2 factor in defn; hopefully ok with everyone?
!
      integer :: i
      real, dimension (mx,my,mz,mvar+maux) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.
!
!  wavenumber k
!
!  set x-wave
!
      if (present(kx)) then
        k=kx
        if (ampl==0) then
          if (lroot) print*,'gaussian: ampl=0; kx=',k
        else
          if (lroot) print*,'gaussian: kx,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(exp(-.5*(k*x)**2),2,my),3,mz)
        endif
      endif
!
!  set y-wave
!
      if (present(ky)) then
        k=ky
        if (ampl==0) then
          if (lroot) print*,'gaussian: ampl=0; ky=',k
        else
          if (lroot) print*,'gaussian: ky,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(exp(-.5*(k*y)**2),1,mx),3,mz)
        endif
      endif
!
!  set z-wave
!
      if (present(kz)) then
        k=kz
        if (ampl==0) then
          if (lroot) print*,'gaussian: ampl=0; kz=',k
        else
          if (lroot) print*,'gaussian: kz,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(exp(-.5*(k*z)**2),1,mx),2,my)
        endif
      endif
!
    endsubroutine gaussian
!***********************************************************************
    subroutine gaussian3d(ampl,f,i,xx,yy,zz,radius)
!
!  gaussian 3-D bump
!
!  28-may-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz) :: xx,yy,zz
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: ampl,radius,radius21
!
      radius21=1./radius**2
      f(:,:,:,i)=f(:,:,:,i)+ampl*exp(-(xx**2+yy**2+zz**2)*radius21)
!
    endsubroutine gaussian3d
!***********************************************************************
    subroutine parabola(ampl,f,i,kx,ky,kz)
!
!  gaussian bump
!
!   2-may-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar+maux) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.
!
!  wavenumber k
!
!  set x-wave
!
      if (present(kx)) then
        k=kx
        if (ampl==0) then
          if (lroot) print*,'parabola: ampl=0; kx=',k
        else
          if (lroot) print*,'parabola: kx,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread((-(k*x)**2),2,my),3,mz)
        endif
      endif
!
!  set y-wave
!
      if (present(ky)) then
        k=ky
        if (ampl==0) then
          if (lroot) print*,'parabola: ampl=0; ky=',k
        else
          if (lroot) print*,'parabola: ky,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread((-(k*y)**2),1,mx),3,mz)
        endif
      endif
!
!  set z-wave
!
      if (present(kz)) then
        k=kz
        if (ampl==0) then
          if (lroot) print*,'parabola: ampl=0; kz=',k
        else
          if (lroot) print*,'parabola: kz,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread((-(k*z)**2),1,mx),2,my)
        endif
      endif
!
    endsubroutine parabola
!***********************************************************************
    subroutine wave(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave
!
!   6-jul-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar+maux) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.
!
!  wavenumber k
!
!  set x-wave
!
      if (present(kx)) then
        k=kx
        if (ampl==0) then
          if (lroot) print*,'wave: ampl=0; kx=',k
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
          if (lroot) print*,'wave: ampl=0; ky=',k
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
          if (lroot) print*,'wave: ampl=0; kz=',k
        else
          if (lroot) print*,'wave: kz,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(sin(k*z),1,mx),2,my)
        endif
      endif
!
    endsubroutine wave
!***********************************************************************
    subroutine wave_uu(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave
!
!  14-apr-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar+maux) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.
!
!  wavenumber k
!
!  set x-wave
!
      if (present(kx)) then
        k=kx
        if (ampl==0) then
          if (lroot) print*,'wave_uu: ampl=0; kx=',k
        else
          if (lroot) print*,'wave_uu: kx,i=',k,i
          f(:,:,:,i)=alog(1.+ampl*spread(spread(sin(k*x),2,my),3,mz)*f(:,:,:,iux))
        endif
      endif
!
!  set y-wave
!
      if (present(ky)) then
        k=ky
        if (ampl==0) then
          if (lroot) print*,'wave_uu: ampl=0; ky=',k
        else
          if (lroot) print*,'wave_uu: ky,i=',k,i
          f(:,:,:,i)=alog(1.+ampl*spread(spread(sin(k*y),1,mx),3,mz)*f(:,:,:,iuy))
        endif
      endif
!
!  set z-wave
!
      if (present(kz)) then
        k=kz
        if (ampl==0) then
          if (lroot) print*,'wave_uu: ampl=0; kz=',k
        else
          if (lroot) print*,'wave_uu: kz,i=',k,i,iuz
          f(:,:,:,i)=alog(1.+ampl*spread(spread(sin(k*z),1,mx),2,my)*f(:,:,:,iuz))
        endif
      endif
!
    endsubroutine wave_uu
!***********************************************************************
    subroutine modes(ampl,coef,f,i,kx,ky,kz,xx,yy,zz)
!
!  mode
!
!  30-oct-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: tmp,xx,yy,zz
      complex :: coef
      complex :: ii=(0.,1.)
      real :: ampl,kx,ky,kz
!
      f(:,:,:,i)=ampl*real(coef*exp(ii*(kx*xx+ky*yy+kz*zz)))
!
    endsubroutine modes
!***********************************************************************
    subroutine modev(ampl,coef,f,i,kx,ky,kz,xx,yy,zz)
!
!  mode
!
!  30-oct-03/axel: coded
!
      integer :: i,ivv
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: tmp,xx,yy,zz
      complex, dimension (3) :: coef
      complex :: ii=(0.,1.)
      real :: ampl,kx,ky,kz
!
      do ivv=0,2
        f(:,:,:,ivv+i)=ampl*real(coef(ivv+1)*exp(ii*(kx*xx+ky*yy+kz*zz)))
      enddo
!
    endsubroutine modev
!***********************************************************************
    subroutine modeb(ampl,coefb,f,i,kx,ky,kz,xx,yy,zz)
!
!  mode
!
!  30-oct-03/axel: coded
!
      integer :: i,ivv
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: tmp,xx,yy,zz
      complex, dimension (3) :: coef,coefb
      complex :: ii=(0.,1.)
      real :: ampl,kx,ky,kz,k2
!
      print*,'print Ak coefficients from Bk coefficients'
      k2=kx**2+ky**2+kz**2
      coef(1)=(ky*coefb(3)-kz*coefb(2))/k2
      coef(2)=(kz*coefb(1)-kx*coefb(3))/k2
      coef(3)=(kx*coefb(2)-ky*coefb(1))/k2
      print*,'coef=',coef
!
      do ivv=0,2
        f(:,:,:,ivv+i)=ampl*real(coef(ivv+1)*exp(ii*(kx*xx+ky*yy+kz*zz)))
      enddo
!
    endsubroutine modeb
!***********************************************************************
    subroutine jump(f,i,fleft,fright,width,dir)
!
!  jump
!
!  19-sep-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx) :: profx
      real, dimension (my) :: profy
      real, dimension (mz) :: profz
      real :: fleft,fright,width
      character(len=*) :: dir
!
!  jump; check direction
!
      select case(dir)
!
      case('x')
        profx=fleft+(fright-fleft)*.5*(1.+tanh(x/width))
        f(:,:,:,i)=f(:,:,:,i)+spread(spread(profx,2,my),3,mz)
!
      case('y')
        profy=fleft+(fright-fleft)*.5*(1.+tanh(y/width))
        f(:,:,:,i)=f(:,:,:,i)+spread(spread(profy,1,mx),3,mz)
!
      case('z')
        profz=fleft+(fright-fleft)*.5*(1.+tanh(z/width))
        f(:,:,:,i)=f(:,:,:,i)+spread(spread(profz,1,mx),2,my)
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
        print*,'bjump: no default value'
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.,fac
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  set x-dependent Beltrami field
!
      if (present(kx)) then
        k=kx; if(k==0) print*,'beltrami: k must not be zero!'; fac=sqrt(abs(ampl/k))
        if (ampl==0) then
          if (lroot) print*,'beltrami: ampl=0; kx=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami: Beltrami field (pos-hel): kx,i=',k,i
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*x),2,my),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*x),2,my),3,mz)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami: Beltrami field (neg-hel): kx,i=',k,i
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*x),2,my),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*x),2,my),3,mz)
        endif
      endif
!
!  set y-dependent Beltrami field
!
      if (present(ky)) then
        k=ky; if(k==0) print*,'beltrami: k must not be zero!'; fac=sqrt(abs(ampl/k))
        if (ampl==0) then
          if (lroot) print*,'beltrami: ampl=0; ky=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami: Beltrami field (pos-hel): ky,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*y),1,mx),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*y),1,mx),3,mz)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami: Beltrami field (neg-hel): ky,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*y),1,mx),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*y),1,mx),3,mz)
        endif
      endif
!
!  set z-dependent Beltrami field
!
      if (present(kz)) then
        k=kz; if(k==0) print*,'beltrami: k must not be zero!'; fac=sqrt(abs(ampl/k))
        if (ampl==0) then
          if (lroot) print*,'beltrami: ampl=0; kz=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami: Beltrami field (pos-hel): kz,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*z),1,mx),2,my)
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*z),1,mx),2,my)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami: Beltrami field (neg-hel): kz,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*z),1,mx),2,my)
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*z),1,mx),2,my)
        endif
      endif
!
    endsubroutine beltrami
!***********************************************************************
    subroutine soundwave(ampl,f,i,kx,ky,kz)
!
!  sound wave (as initial condition)
!
!   2-aug-02/axel: adapted from Beltrami
!
      integer :: i
      real, dimension (mx,my,mz,mvar+maux) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.,fac
!
!  wavenumber k
!
!  set x-dependent sin wave
!
      if (present(kx)) then
        k=kx; if(k==0) print*,'soundwave: k must not be zero!'; fac=sqrt(abs(ampl/k))
        if (ampl==0) then
          if (lroot) print*,'soundwave: ampl=0; kx=',k
        else
          if (lroot) print*,'soundwave: kx,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(sin(k*x),2,my),3,mz)
        endif
      endif
!
!  set y-dependent sin wave field
!
      if (present(ky)) then
        k=ky; if(k==0) print*,'soundwave: k must not be zero!'; fac=sqrt(abs(ampl/k))
        if (ampl==0) then
          if (lroot) print*,'soundwave: ampl=0; ky=',k
        else
          if (lroot) print*,'soundwave: ky,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(sin(k*y),1,mx),3,mz)
        endif
      endif
!
!  set z-dependent sin wave field
!
      if (present(kz)) then
        k=kz; if(k==0) print*,'soundwave: k must not be zero!'; fac=sqrt(abs(ampl/k))
        if (ampl==0) then
          if (lroot) print*,'soundwave: ampl=0; kz=',k
        else
          if (lroot) print*,'soundwave: kz,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(sin(k*z),1,mx),2,my)
        endif
      endif
!
    endsubroutine soundwave
!***********************************************************************
    subroutine stratification(ampl,f,xx,yy,zz)
!
!  read mean stratification from "stratification.dat"
!
!   8-apr-03/axel: coded
!
      use Mpicomm
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real, dimension (mz) :: lnrho0,SS0
      real :: ztmp,ampl
      logical :: exist
!
!  read mean stratification and write into array
!  if file is not found in run directory, search under trim(directory)
!
      inquire(file='stratification.dat',exist=exist)
      if(exist) then
        open(19,file='stratification.dat')
      else
        inquire(file=trim(directory)//'/stratification.ascii',exist=exist)
        if(exist) then
          open(19,file=trim(directory)//'/stratification.ascii')
        else
          call stop_it('stratification: *** error *** - no input file')
        endif
      endif
!
!  read data
!
      do n=1,nz*ipz
        read(19,*) ztmp,lnrho0(n),SS0(n)
      enddo
!
      do n=1,mz
        read(19,*) ztmp,lnrho0(n),SS0(n)
        if(ip<5) print*,"stratification: ",ztmp,lnrho0(n),SS0(n)
        f(:,:,n,ilnrho)=lnrho0(n)
        f(:,:,n,iss)=SS0(n)
      enddo
      close(19)
!
      if(ip==0) print*,ampl,xx,yy,zz !(to keep compiler quiet)
    endsubroutine stratification
!***********************************************************************
    subroutine planet_hc(ampl,f,xx,yy,zz,eps,radius,gamma,cs20,rho0,width)
!
!  Ellipsoidal planet solution (Goldreich, Narayan, Goodman 1987)
!
!   6-jul-02/axel: coded
!  22-feb-03/axel: fixed 3-D background solution for enthalpy
!  26-Jul-03/anders: Revived from June 1 version
!
      use Mpicomm, only: mpireduce_sum, mpibcast_real_scl
      
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz,hh,xi
      real, dimension (mx,my) :: delS
      real :: ampl,sigma2,sigma,delta2,delta,eps,radius,a_ell,b_ell,c_ell
      real :: gamma,cs20,gamma1,eps2,radius2,width
      real :: lnrhosum_thisbox,rho0
      real, dimension(1) :: lnrhosum_thisbox_tmp,lnrhosum_wholebox
!
!  calculate sigma
!
      if (lroot) print*,'planet_hc: qshear,eps=',qshear,eps
      eps2=eps**2
      radius2=radius**2
      sigma2=2*qshear/(1.-eps2)
      if (sigma2<0.) then
        if (lroot) print*, &
          'planet_hc: sigma2<0 not allowed; choose another value of eps_planet'
      else
        sigma=sqrt(sigma2)
      endif
!
!  calculate delta
!
      delta2=(2.-sigma)*sigma
      if (lroot) print*,'planet_hc: sigma,delta2,radius=',sigma,delta2,radius
      if (delta2<=0.) then
        if (lroot) print*,'planet_hc: delta2<=0 not allowed'
      else
        delta=sqrt(delta2)
      endif
!
!  calculate gamma1
!
      gamma1=gamma-1.
      if (lroot) print*,'planet_hc: gamma=',gamma
!
!  calculate hh
!  xi=1 inside vortex, and 0 outside
!
      hh=+.5*delta2*Omega**2*(radius2-xx**2-eps2*yy**2)-.5*Omega**2*zz**2
      xi=.5+.5*tanh(hh/width)
!
!  ellipse parameters
!
      b_ell = radius
      a_ell = radius/eps
      c_ell = radius*delta
      if(lroot) print*,'planet_hc: Ellipse axes (b_ell,a_ell,c_ell)=', &
          b_ell,a_ell,c_ell
!
!  add continuous vertical stratification to horizontal planet solution
!  NOTE: if width is too small, the vertical integration below may fail.
!
      if (lroot) print*,"planet_hc: integrate hot corona"
      hh(:,:,n2)=1.  !(initial condition)
      f(:,:,:,iss)=-alog(ampl)*xi
      do n=n2-1,n1,-1
        delS=f(:,:,n+1,iss)-f(:,:,n,iss)
        hh(:,:,n)=(hh(:,:,n+1)*(1.-.5*delS)+ &
             Omega**2*.5*(z(n)+z(n+1))*dz)/(1.+.5*delS)
      enddo
!
!  Calculate velocities (Kepler speed subtracted)
!
      f(:,:,:,iux)=   eps2*sigma *Omega*yy*xi
      f(:,:,:,iuy)=(qshear-sigma)*Omega*xx*xi
!
      if (lroot) print*,'planet_hc: hmin=',minval(hh(l1:l2,m1:m2,n1:n2))

      if(gamma1<0. .and. lroot) &
          print*,'planet_hc: must have gamma>1 for planet solution'
!
!  calculate density, depending on what gamma is
!
      if(lentropy) then
        f(l1:l2,m1:m2,n1:n2,ilnrho)= &
             (alog(gamma1*hh(l1:l2,m1:m2,n1:n2)/cs20) &
             -gamma*f(l1:l2,m1:m2,n1:n2,iss))/gamma1
        if (lroot) &
          print*,'planet_hc: planet solution with entropy for gamma=',gamma
      else
        if(gamma==1.) then
          f(l1:l2,m1:m2,n1:n2,ilnrho)=hh(l1:l2,m1:m2,n1:n2)/cs20
          if (lroot) print*,'planet_hc: planet solution for gamma=1'
        else
          f(l1:l2,m1:m2,n1:n2,ilnrho)=&
               alog(gamma1*hh(l1:l2,m1:m2,n1:n2)/cs20)/gamma1
          if (lroot) print*,'planet_hc: planet solution for gamma=',gamma
        endif
      endif
!
!  Use average density of box as unit density
!
      lnrhosum_thisbox = sum(f(l1:l2,m1:m2,n1:n2,ilnrho))
      if (ip<14) &
        print*,'planet_hc: iproc,lnrhosum_thisbox=',iproc,lnrhosum_thisbox
!
!  Must put sum_thisbox in 1-dimensional array 
!
      lnrhosum_thisbox_tmp = (/ lnrhosum_thisbox /)
!
!  Add sum_thisbox up for all processors, deliver to root
!
      call mpireduce_sum(lnrhosum_thisbox_tmp,lnrhosum_wholebox,1)
      if (lroot .and. ip<14) &
        print*,'planet_hc: lnrhosum_wholebox=',lnrhosum_wholebox
!
!  Calculate <rho> and send to all processors
!      
      if (lroot) rho0 = exp(-lnrhosum_wholebox(1)/(nxgrid*nygrid*nzgrid))
      call mpibcast_real_scl(rho0,1)
      if (ip<14) print*,'planet_hc: iproc,rho0=',iproc,rho0
!
!  Multiply density by rho0 (divide by <rho>)
!      
      f(l1:l2,m1:m2,n1:n2,ilnrho) = f(l1:l2,m1:m2,n1:n2,ilnrho) + alog(rho0)
!
    endsubroutine planet_hc
!***********************************************************************
    subroutine planet(rbound,f,xx,yy,zz,eps,radius,gamma,cs20,rho0,width,hh0)
!
!  Cylindrical planet solution (Goldreich, Narayan, Goodman 1987)
!
!   jun-03/anders: coded (adapted from old 'planet', now 'planet_hc')
!
      use Mpicomm, only: mpireduce_sum, mpibcast_real_scl

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz,hh,xi,r_ell
      real :: rbound,sigma2,sigma,delta2,delta,eps,radius
      real :: gamma,eps2,radius2,width,a_ell,b_ell,c_ell
      real :: gamma1,ztop,cs20,hh0
      real :: lnrhosum_thisbox,rho0
      real, dimension(1) :: lnrhosum_thisbox_tmp,lnrhosum_wholebox
!
!  calculate sigma
!
      if (lroot) print*,'planet: qshear,eps=',qshear,eps
      eps2=eps**2
      radius2=radius**2
      sigma2=2*qshear/(1.-eps2)
      if (sigma2<0. .and. lroot) then
        print*, &
            'planet: sigma2<0 not allowed; choose another value of eps_planet'
      else
        sigma=sqrt(sigma2)
      endif

      gamma1=gamma-1.
!
!  calculate delta
!
      delta2=(2.-sigma)*sigma
      if (lroot) print*,'planet: sigma,delta2,radius=',sigma,delta2,radius
      if (delta2<0. .and. lroot) then
        print*,'planet: delta2<0 not allowed'
      else
        delta=sqrt(delta2)
      endif

      ztop=z0+lz
      if (lroot) print*,'planet: ztop=', ztop
!
!  Cylinder vortex 3-D solution (b_ell along x, a_ell along y)
!
      b_ell = radius
      a_ell = radius/eps
      c_ell = radius*delta
      r_ell = sqrt(xx**2/b_ell**2+yy**2/a_ell**2)
      if(lroot) print*,'planet: Ellipse axes (b_ell,a_ell,c_ell)=', &
          b_ell,a_ell,c_ell
!
!  xi is 1 inside vortex and 0 outside
!
      xi = 1./(exp((1/width)*(r_ell-rbound))+1.)
      if(lroot) print*,'planet: width,rbound',width,rbound
!
!  Calculate enthalpy inside vortex
!
      hh = 0.5*delta2*Omega**2*(radius2-xx**2-eps2*yy**2) &
           -0.5*Omega**2*zz**2 + 0.5*Omega**2*ztop**2 + hh0
!
!  Calculate enthalpy outside vortex
!
      where (r_ell .gt. 1) hh=-0.5*Omega**2*zz**2 + 0.5*Omega**2*ztop**2 + hh0 
!
!  Calculate velocities (Kepler speed subtracted)
!
      f(:,:,:,iux)=   eps2*sigma *Omega*yy*xi
      f(:,:,:,iuy)=(qshear-sigma)*Omega*xx*xi
!
!  calculate density, depending on what gamma is
!
      if (lroot) print*,'planet: hh0,hmin',&
           hh0,minval(hh(l1:l2,m1:m2,n1:n2))
      if (lentropy .and. lroot) print*,'planet: smin,smax', &
           minval(f(:,:,:,iss)), maxval(f(:,:,:,iss))
      if (gamma1<0. .and. lroot) & 
           print*,'planet: must have gamma>1 for planet solution'
!
!  have to use explicit indices here, because ghostzones are not set
!
      if(lentropy) then
        f(l1:l2,m1:m2,n1:n2,ilnrho) = &
            (alog(gamma1*hh(l1:l2,m1:m2,n1:n2)/cs20) &
            - gamma*f(l1:l2,m1:m2,n1:n2,iss))/gamma1
        if (lroot) print*,'planet: planet solution for gamma=',gamma
      else
      if(gamma==1.) then
        f(l1:l2,m1:m2,n1:n2,ilnrho) = hh(l1:l2,m1:m2,n1:n2)/cs20
          if (lroot) print*,'planet: planet solution for gamma=1'
        else
          f(l1:l2,m1:m2,n1:n2,ilnrho) = &
               alog(gamma1*hh(l1:l2,m1:m2,n1:n2)/cs20)/gamma1
          if (lroot) print*,'planet: planet solution for gamma=',gamma
        endif
      endif
!
!  Use average density of box as unit density
!
      lnrhosum_thisbox = sum(f(l1:l2,m1:m2,n1:n2,ilnrho))
      if (ip<14) &
        print*,'planet_hc: iproc,lnrhosum_thisbox=',iproc,lnrhosum_thisbox
!
!  Must put sum_thisbox in 1-dimensional array 
!
      lnrhosum_thisbox_tmp = (/ lnrhosum_thisbox /)
!
!  Add sum_thisbox up for all processors, deliver to root
!
      call mpireduce_sum(lnrhosum_thisbox_tmp,lnrhosum_wholebox,1)
      if (lroot .and. ip<14) &
        print*,'planet_hc: lnrhosum_wholebox=',lnrhosum_wholebox
!
!  Calculate <rho> and send to all processors
!      
      if (lroot) rho0 = exp(-lnrhosum_wholebox(1)/(nxgrid*nygrid*nzgrid))
      call mpibcast_real_scl(rho0,1)
      if (ip<14) print*,'planet_hc: iproc,rho0=',iproc,rho0
!
!  Multiply density by rho0 (divide by <rho>)
!      
      f(l1:l2,m1:m2,n1:n2,ilnrho) = f(l1:l2,m1:m2,n1:n2,ilnrho) + alog(rho0)
!      
    endsubroutine planet
!***********************************************************************
    subroutine baroclinic(f,xx,yy,zz,gamma,t_fct_type,rhomid_fct_type, &
                          TT_til0,alpha_TT,rho0,beta_rho)
!
!  Baroclinic shearing sheet initial condition
!  11-nov-03/anders: coded
!
      integer :: i,j,k,izmid,hr_q,mz_hr,izmid_hr
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz,I_int,rho,TT_til,dlnrhodx
      real, dimension (1000*mz) :: z_hr,I_int_z,TT_til_hr,dPhidz_hr
      real :: gamma,TT_til0,alpha_TT,rho0,beta_rho,dz_hr
      character (len=labellen):: t_fct_type,rhomid_fct_type
!
!  Temperature
!
      select case(t_fct_type)
        case('linear_z')
          TT_til = TT_til0 + alpha_TT * abs(zz)
          TT_til_hr = TT_til0 + alpha_TT * abs(z_hr)
!
!  Need high resolution in z to calculate integral
!
          hr_q     = 1000
          mz_hr    = hr_q*mz
!
!  Index number of mid-plane
!
          izmid    = mz/2
          izmid_hr = mz_hr/2
!
!  Calculate high resolution z array
!      
          dz_hr = (zz(1,1,mz)-zz(1,1,1))/mz_hr
          z_hr(1)=zz(1,1,1)
          do i=2,mz_hr; z_hr(i)=z_hr(i-1)+dz_hr; enddo
!
!  The potential is Omega**2*( -r0^2 + r0*xx - xx^2 + 0.5*zz^2 )
!
          dPhidz_hr = Omega**2*z_hr
!
!  Integral for calculating hydrostatic equlibrium density
!
          I_int_z(izmid_hr:mz_hr) = spline_integral(z_hr(izmid_hr:mz_hr), &
              1/TT_til_hr(izmid_hr:mz_hr)*dPhidz_hr(izmid_hr:mz_hr))
          I_int_z(izmid_hr:1:-1) = spline_integral(z_hr(izmid_hr:1:-1), &
              1/TT_til_hr(izmid_hr:1:-1)*dPhidz_hr(izmid_hr:1:-1))
          print*,I_int_z(1),I_int_z(mz_hr)
!
!  Put high resolution integral back into normal resolution array
!
          do i=1,mx
            do j=1,my
              do k=1,mz
                I_int(i,j,k) = I_int_z(hr_q*(k-1)+1)
              enddo
            enddo
          enddo
      endselect
!
!  Density in the mid-plane
!
      select case(rhomid_fct_type)
        case ('linear_x')
          f(:,:,izmid,ilnrho) = alog(rho0 + beta_rho * xx(:,:,izmid))
      endselect
!
!  Solution to hydrostatic equlibrium in the z-direction
!
      do k=1,mz
        f(:,:,k,ilnrho) = f(:,:,izmid,ilnrho) - &
            alog(TT_til(:,:,k)/TT_til(:,:,izmid)) - I_int(:,:,k)
      enddo

      rho = exp(f(:,:,:,ilnrho))
     
      do j=0,my
        do k=0,mz
          dlnrhodx(:,j,k) = spline_derivative(xx(:,j,k),f(:,j,k,ilnrho))
        enddo
      enddo
!
!  Toroidal velocity comes from hyd. stat. eq. equ. in the x-direction
!  (potential term vansihes when velocities measured relative to -1.5 O x.
!
      f(:,:,:,iuy) = 1/(2*Omega)*TT_til*dlnrhodx(:,:,:)
!
!  Entropy is calculated from temperature and density ( P = T~*rho )
!
      f(:,:,:,iss) = 1/gamma*alog(TT_til(:,:,:)*rho(:,:,:)) - f(:,:,:,ilnrho)
!
    endsubroutine baroclinic
!***********************************************************************
    subroutine crazy(ampl,f,i)
!
!  A crazy initial condition
!  (was useful to initialize all points with finite values)
!
!  19-may-02/axel: coded
!
      integer :: i,j
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: ampl
!
      if (lroot) print*, 'crazy: sinusoidal magnetic field: for debugging purposes'
      j=i; f(:,:,:,j)=f(:,:,:,j)+ampl*&
        spread(spread(sin(2*x),2,my),3,mz)*&
        spread(spread(sin(3*y),1,mx),3,mz)*&
        spread(spread(cos(1*z),1,mx),2,my)
      j=i+1; f(:,:,:,j)=f(:,:,:,j)+ampl*&
        spread(spread(sin(5*x),2,my),3,mz)*&
        spread(spread(sin(1*y),1,mx),3,mz)*&
        spread(spread(cos(2*z),1,mx),2,my)
      j=i+2; f(:,:,:,j)=f(:,:,:,j)+ampl*&
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: tmp,xx,yy,zz,modulate
      real :: ampl,radius,epsilon_nonaxi,ky
!
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'htube: set variable to zero; i1,i2=',i1,i2
      else
        ky=2*pi/Ly
        if(lroot) then
          print*,'htube: implement y-dependent flux tube in xz-plane; i1,i2=',i1,i2
          print*,'htube: radius,epsilon_nonaxi=',radius,epsilon_nonaxi
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: tmp,xx,yy,zz,modulate
      real :: ampl,radius,epsilon_nonaxi,ky
!
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'htube2: set variable to zero; i1,i2=',i1,i2
      else
        ky=2*pi/Ly
        if(lroot) then
          print*,'htube2: implement y-dependent flux tube in xz-plane; i1,i2=',i1,i2
          print*,'htube2: radius,epsilon_nonaxi=',radius,epsilon_nonaxi
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
          if(lroot) print*,'htube2: set scalar'
          f(:,:,:,i1)=tmp
        elseif(i1+2==i2) then
          if(lroot) print*,'htube2: set vector'
          f(:,:,:,i1 )=+zz*tmp
          f(:,:,:,i1+1)=0.
          f(:,:,:,i1+2)=-xx*tmp
        else
          if(lroot) print*,'htube2: bad value of i2=',i2
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: zz
      real :: ampl,H,A0,gravz,cs0,rho0,lnrho0
!
      if (ampl==0) then
        if (lroot) print*,'magsupport: do nothing'
      else
        lnrho0=alog(rho0)
        H=(1+ampl)*cs0**2/abs(gravz)
        A0=-2*H*ampl*cs0*sqrt(2*rho0)
        if (lroot) print*,'magsupport: H,A0=',H,A0
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: ampl,zflayer,width
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'hfluxlayer: set variable to zero; i=',i
      else
        if (lroot) print*,'hfluxlayer: horizontal flux layer; i=',i
        if ((ip<=16).and.lroot) print*,'hfluxlayer: ampl,width=',ampl,width
        f(:,:,:,i  )=0.
        f(:,:,:,i+1)=ampl*tanh((zz-zflayer)/width)
        f(:,:,:,i+2)=0.
      endif
!
      if (ip==1) print*,xx,yy
    endsubroutine hfluxlayer
!***********************************************************************
    subroutine halfcos_x(ampl,f,i,xx,yy,zz)
!
!  Uniform B_x field (for vector potential)
!
!  19-jun-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: ampl,kz,zbot
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'halfcos_x: set variable to zero; i=',i
      else
        print*,'halscos_x: half cosine x-field ; i=',i
        kz=0.5*pi/Lz
        zbot=xyz0(3)
        ! ztop=xyz0(3)+Lxyz(3)
        if ((ip<=16).and.lroot) print*,'halfcos_x: ampl,kz=',ampl,kz
        f(:,:,:,i  )=0.
        f(:,:,:,i+1)=-ampl*sin(kz*(zz-zbot))
        f(:,:,:,i+2)=0.
      endif
!
      if (ip==1) print*,xx,yy
    endsubroutine halfcos_x
!***********************************************************************
    subroutine uniform_x(ampl,f,i,xx,yy,zz)
!
!  Uniform B_x field (for vector potential)
!
!  19-jun-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: ampl
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'uniform_x: set variable to zero; i=',i
      else
        print*,'uniform_x: uniform x-field ; i=',i
        if ((ip<=16).and.lroot) print*,'uniform_x: ampl=',ampl
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: ampl
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'uniform_y: set variable to zero; i=',i
      else
        print*,'uniform_y: uniform y-field ; i=',i
        if ((ip<=16).and.lroot) print*,'uniform_y: ampl=',ampl
        f(:,:,:,i  )=ampl*zz
        f(:,:,:,i+1)=0.
        f(:,:,:,i+2)=0.
      endif
!
      if (ip==1) print*,xx,yy
    endsubroutine uniform_y
!***********************************************************************
    subroutine uniform_z(ampl,f,i,xx,yy,zz)
!
!  Uniform B_x field (for vector potential)
!
!  24-jul-03/axel: adapted from uniform_x
!
      integer :: i
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: ampl
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'uniform_z: set variable to zero; i=',i
      else
        print*,'uniform_z: uniform z-field ; i=',i
        if ((ip<=16).and.lroot) print*,'uniform_z: ampl=',ampl
        f(:,:,:,i  )=0.
        f(:,:,:,i+1)=+ampl*xx
        f(:,:,:,i+2)=0.
      endif
!
      if (ip==1) print*,yy,zz
    endsubroutine uniform_z
!***********************************************************************
    subroutine vfield(ampl,f,i,xx)
!
!  Vertical field, for potential field test
!
!  14-jun-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx
      real :: ampl,kx
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'vfield: set variable to zero; i=',i
      else
        kx=2*pi/Lx
        print*,'vfield: implement x-dependent vertical field'
        if ((ip<=8).and.lroot) print*,'vfield: x-dependent vertical field'
        f(:,:,:,i  )=0.
        f(:,:,:,i+1)=ampl*sin(kx*xx)
        f(:,:,:,i+2)=0.
      endif
!
    endsubroutine vfield
!***********************************************************************
    subroutine gaunoise_vect(ampl,f,i1,i2)
!
!  Add Gaussian (= normally distributed) white noise for variables i1:i2
!
!  23-may-02/axel: coded
!  10-sep-03/axel: result only *added* to whatever f array had before
!
      integer :: i,i1,i2
      real, dimension (mx,my,mz) :: r,p,tmp
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: ampl
!
!  set gaussian random noise vector
!
      if (ampl==0) then
        if (lroot) print*,'gaunoise_vect: ampl=0 for i1,i2=',i1,i2
      else
        if ((ip<=8).and.lroot) print*,'gaunoise_vect: i1,i2=',i1,i2
        do i=i1,i2
          if (modulo(i-i1,2)==0) then
            call random_number_wrapper(r)
            call random_number_wrapper(p)
            tmp=sqrt(-2*alog(r))*sin(2*pi*p)
          else
            tmp=sqrt(-2*alog(r))*cos(2*pi*p)
          endif
          !call smooth_3d(tmp,ismo)  !(may want to smooth)
          f(:,:,:,i)=f(:,:,:,i)+ampl*tmp
          if (lroot) print*,'gaunoise_vect: variable i=',i
        enddo
      endif
!
    endsubroutine gaunoise_vect
!***********************************************************************
    subroutine gaunoise_scal(ampl,f,i)
!
!  Ad Gaussian (= normally distributed) white noise for variable i
!
!  23-may-02/axel: coded
!  10-sep-03/axel: result only *added* to whatever f array had before
!
      integer :: i
      real, dimension (mx,my,mz) :: r,p,tmp
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: ampl
!
!  set gaussian random noise vector
!
      if (ampl==0) then
        if (lroot) print*,'gaunoise_scal: ampl=0 for i=',i
      else
        if ((ip<=8).and.lroot) print*,'gaunoise_scal: i=',i
        call random_number_wrapper(r)
        call random_number_wrapper(p)
        tmp=sqrt(-2*alog(r))*sin(2*pi*p)
        f(:,:,:,i)=f(:,:,:,i)+ampl*tmp
        if (lroot) print*,'gaunoise_scal: variable i=',i
      endif
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: tmp,xx,yy,zz
      real :: ampl
!
      if (lroot) print*, 'trilinear: ivar = ', ivar
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: ampl,kx,ky,kz
!
      if (lroot) print*, 'cos_cos_sin: ivar = ', ivar
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: ampl,ky,kz
!
      if (lroot) print*, 'tor_pert: sinusoidal modulation of ivar = ', ivar
!
      ky=2*pi/Ly
      kz=2.*pi/Lz
      f(:,:,:,ivar) = ampl*cos(ky*yy)*cos(kz*zz)
!
      if (ip==0) print*,'xx(1,1,1)=',xx(1,1,1) !(to keep compiler quiet)
    endsubroutine tor_pert
!***********************************************************************
    subroutine diffrot(ampl,f,ivar,xx,yy,zz)
!
!  Set up profile for differential rotation
!
!  16-jul-03/axel: coded
!
      integer :: ivar
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: ampl
!
      if (lroot) print*, 'diffrot: sinusoidal modulation of ivar = ', ivar
!
      f(:,:,:,ivar) = ampl*cos(xx)*cos(zz)
!
      if (ip==0) print*,'yy(1,1,1)=',yy(1,1,1) !(to keep compiler quiet)
    endsubroutine diffrot
!***********************************************************************
    subroutine olddiffrot(ampl,f,ivar,xx,yy,zz)
!
!  Set up profile for differential rotation
!
!  16-jul-03/axel: coded
!
      integer :: ivar
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: ampl,kx,kz
!
      if (lroot) print*, 'olddiffrot: sinusoidal modulation of ivar = ', ivar
!
      kx=.5*pi/Lx
      kz=.5*pi/Lz
      f(:,:,:,ivar) = ampl*sin(kx*xx)*cos(kz*zz)
!
      if(ip==0) print*,yy !(to keep compiler quiet)
    endsubroutine olddiffrot
!***********************************************************************
    subroutine powern(ampl,initpower,cutoff,f,i1,i2)
!
!   Produces k^initpower*exp(-k**2/cutoff**2)  spectrum.
!   Still just one processor (but can be remeshed afterwards).
!
!   07-may-03/tarek: coded
!
      integer :: i,i1,i2
      real, dimension (nx,ny,nz) :: k2
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx,ny,nz) :: u_re,u_im
      real :: ampl,initpower,cutoff
 
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'powern: set variable to zero; i1,i2=',i1,i2
      else
        call gaunoise_vect(ampl,f,i1,i2) ! which has a k^2. spectrum

        if ((initpower.ne.2.).or.(cutoff.ne.0.)) then

          k2=      (spread(spread(cshift &
            ((/(i-(nx+1)/2,i=0,nx-1)/),+(nx+1)/2)*2*pi/Lx,2,ny),3,nz))**2.
          k2= k2 + (spread(spread(cshift &
            ((/(i-(ny+1)/2,i=0,ny-1)/),+(ny+1)/2)*2*pi/Ly,1,nx),3,nz))**2.
          k2= k2 + (spread(spread(cshift &
            ((/(i-(nz+1)/2,i=0,nz-1)/),+(nz+1)/2)*2*pi/Lz,1,nx),2,ny))**2.

          k2(1,1,1) = 1.  ! Avoid division by zero 

          do i=i1,i2
            u_re=f(l1:l2,m1:m2,n1:n2,i)
            u_im=0. 
            !  fft of gausian noise w/ k^2 spectrum
            call transform_fftpack(u_re,u_im,1)
            ! change to k^n spectrum
            u_re =(k2)**(.25*initpower-.5)*u_re 
            u_im =(k2)**(.25*initpower-.5)*u_im 
            ! cutoff (changed to hyperviscous cutoff filter)
            if (cutoff .ne. 0.) then
              u_re = u_re*exp(-(k2/cutoff**2.)**2) 
              u_im = u_im*exp(-(k2/cutoff**2.)**2) 
            endif   
            ! back to real space 
            call transform_fftpack(u_re,u_im,-1)
            f(l1:l2,m1:m2,n1:n2,i)=u_re
            
            if (lroot .and. (cutoff.eq.0)) then 
              print*,'powern: k^',initpower,' spectrum : var  i=',i
            else
              print*,'powern: with cutoff : k^n*exp(-k^4/k0^4) w/ n=', &
                     initpower,', k0 =',cutoff,' : var  i=',i
            endif 
          enddo !i            
      endif !(initpower.ne.2.).or.(cutoff.ne.0.)

    endif !(ampl.eq.0)
!
    endsubroutine powern

!***********************************************************************
    subroutine power_randomphase(ampl,initpower,cutoff,f,i1,i2)
!
!   Produces k^initpower*exp(-k**2/cutoff**2)  spectrum.
!   Still just one processor (but can be remeshed afterwards).
!
!   07-may-03/tarek: coded
!
      integer :: i,i1,i2
      real, dimension (nx,ny,nz) :: k2
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx,ny,nz) :: u_re,u_im,r
      real :: ampl,initpower,mhalf,cutoff
 
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'powern: set variable to zero; i1,i2=',i1,i2
      else
!
!  calculate k^2
!
        k2=      (spread(spread(cshift &
          ((/(i-(nx+1)/2,i=0,nx-1)/),+(nx+1)/2)*2*pi/Lx,2,ny),3,nz))**2.
        k2= k2 + (spread(spread(cshift &
          ((/(i-(ny+1)/2,i=0,ny-1)/),+(ny+1)/2)*2*pi/Ly,1,nx),3,nz))**2.
        k2= k2 + (spread(spread(cshift &
          ((/(i-(nz+1)/2,i=0,nz-1)/),+(nz+1)/2)*2*pi/Lz,1,nx),2,ny))**2.

        k2(1,1,1) = 1.  ! Avoid division by zero 
!
!  To get shell integrated power spectrum E ~ k^n, we need u ~ k^m
!  and since E(k) ~ u^2 k^2 we have n=2m+2, so m=n/2-1.
!  Further, since we operate on k^2, we need m/2 (called mhalf below)
!
        mhalf=.5*(.5*initpower-1)
!
!  generate all 3 velocity components separately
!
        do i=i1,i2
          ! generate k^n spectrum with random phase (between -pi and pi)
          call random_number_wrapper(r); u_re=ampl*k2**mhalf*cos(pi*(2*r-1))
          call random_number_wrapper(r); u_im=ampl*k2**mhalf*sin(pi*(2*r-1))
          ! cutoff (changed to hyperviscous cutoff filter)
          if (cutoff .ne. 0.) then
            u_re = u_re*exp(-(k2/cutoff**2.)**2) 
            u_im = u_im*exp(-(k2/cutoff**2.)**2) 
          endif   
          ! back to real space 
          call transform_fftpack(u_re,u_im,-1)
          f(l1:l2,m1:m2,n1:n2,i)=u_re
          
          if (lroot .and. (cutoff.eq.0)) then 
            print*,'powern: k^',initpower,' spectrum : var  i=',i
          else
            print*,'powern: with cutoff : k^n*exp(-k^4/k0^4) w/ n=', &
                   initpower,', k0 =',cutoff,' : var  i=',i
          endif 
        enddo !i            

      endif !(ampl.eq.0)
!
    endsubroutine power_randomphase

!***********************************************************************
endmodule Initcond
