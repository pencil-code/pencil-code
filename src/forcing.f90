! $Id: forcing.f90,v 1.11 2002-06-02 21:01:59 brandenb Exp $

module Forcing

!! Module for forcing in the Navier-Stokes equation.

  use Cdata

  implicit none

  integer :: iforce=2,iforce2=0,kfountain=5
  real :: force=0.,relhel=1.,height=pi,fountain=1.

  integer :: dummy              ! We cannot define empty namelists
  namelist /forcing_init_pars/ dummy

  namelist /forcing_run_pars/ &
       iforce,force,relhel,height, &
       iforce2,kfountain,fountain

  contains

!***********************************************************************
    subroutine register_forcing()
!
!  add forcing in timestep()
!  11-may-2002/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_forcing called twice')
      first = .false.
!
      lforcing = .true.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$RCSfile: forcing.f90,v $", &
           "$Revision: 1.11 $", &
           "$Date: 2002-06-02 21:01:59 $")
!
    endsubroutine register_forcing
!***********************************************************************
    subroutine addforce(f)
!
!  add forcing at the end of each time step
!  Since forcing is constant during one time step,
!  this can be added as an Euler 1st order step
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
!
!  calculate and add forcing function
!
      if(iforce==1) call forcing_irro(f)
      if(iforce==2) call forcing_hel(f)
      if(iforce==3) call forcing_fountain(f)
!
!  add *additional* forcing function
!
      if(iforce2==1) call forcing_fountain(f)
!
    endsubroutine addforce
!***********************************************************************
    subroutine forcing_irro(f)
!
!  add acoustic forcing function, using a set of precomputed wavevectors
!  This forcing drives pressure waves
!
!  10-sep-01/axel: coded
!
      use Mpicomm
      use Cdata
!
      real :: phase,ffnorm
      real, save :: kav
      real, dimension (2) :: fran
      real, dimension (mx,my,mz,mvar) :: f
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      complex, dimension (3) :: ikk
      integer, dimension(mk), save :: kkx,kky,kkz
      integer, save :: ifirst,nk
      integer :: ik,j,jf
!
      if (ifirst==0) then
        if (lroot) print*,'opening k.dat'
        open(9,file='k.dat')
        read(9,*) nk,kav
        if (lroot) print*,'average k=',kav
        if(nk.gt.mk) then
          if (lroot) print*,'dimension mk in forcing1 is insufficient'
          print*,'nk=',nk,'mk=',mk
          call mpifinalize
        end if
        read(9,*) (kkx(ik),ik=1,nk)
        read(9,*) (kky(ik),ik=1,nk)
        read(9,*) (kkz(ik),ik=1,nk)
        close(9)
        ifirst=1
      endif
!
      call random_number(fran)
      phase=pi*(2*fran(1)-1.)
      ik=nk*.9999*fran(2)+1
      if (ip<=6) print*,'ik,phase,kk=',ik,phase,kkx(ik),kky(ik),kkz(ik),dt,ifirst
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      ffnorm=force*sqrt(kav/dt)*dt
      fx=exp(cmplx(0.,kkx(ik)*x+phase))*ffnorm
      fy=exp(cmplx(0.,kky(ik)*y))
      fz=exp(cmplx(0.,kkz(ik)*z))
!
      ikk(1)=cmplx(0.,kkx(ik))
      ikk(2)=cmplx(0.,kky(ik))
      ikk(3)=cmplx(0.,kkz(ik))
!
      do j=1,3
        jf=j+iux-1
        do n=n1,n2
        do m=m1,m2
          f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+real(ikk(j)*fx(l1:l2)*fy(m)*fz(n))
        enddo
        enddo
      enddo
!
    endsubroutine forcing_irro
!***********************************************************************
    subroutine forcing_hel(f)
!
!  add helical forcing function, using a set of precomputed wavevectors
!
!  10-apr-00/axel: coded
!
      use Mpicomm
      use Cdata
      use Sub
      use Hydro
!
      real :: phase,ffnorm
      real, save :: kav
      real, dimension (2) :: fran
      real, dimension (mz) :: tmpz
      real, dimension (mx,my,mz,mvar) :: f
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      complex, dimension (3) :: coef
      integer, dimension(mk), save :: kkx,kky,kkz
      integer, save :: ifirst,nk
      integer :: ik,j,jf,kx,ky,kz,kex,key,kez,kkex,kkey,kkez
      real :: k2,k,ex,ey,ez,kde,sig=1.,fact
!
      if (ifirst==0) then
        if (lroot) print*,'opening k.dat'
        open(9,file='k.dat')
        read(9,*) nk,kav
        if (lroot) print*,'average k=',kav
        if(nk.gt.mk) then
          if (lroot) print*,'dimension mk in forcing1 is insufficient'
          print*,'nk=',nk,'mk=',mk
          call mpifinalize
        end if
        read(9,*) (kkx(ik),ik=1,nk)
        read(9,*) (kky(ik),ik=1,nk)
        read(9,*) (kkz(ik),ik=1,nk)
        close(9)
        ifirst=1
      endif
!
!  generate random coefficients -1 < fran < 1
!  ff=force*Re(exp(i(kx+phase)))
!  |k_i| < akmax
!
      call random_number(fran)
      phase=pi*(2*fran(1)-1.)
      ik=nk*.9999*fran(2)+1
      if (ip<=6) print*,'ik,phase,kk=',ik,phase,kkx(ik),kky(ik),kkz(ik),dt,ifirst
!
      kx=kkx(ik)
      ky=kky(ik)
      kz=kkz(ik)
      if(ip.le.4) print*, 'kx,ky,kz=',kx,ky,kz
!
      k2=float(kx**2+ky**2+kz**2)
      k=sqrt(k2)
!
!  pick e1 if kk not parallel to ee1. ee2 else.
!
      if((ky.eq.0).and.(kz.eq.0)) then
        ex=0; ey=1; ez=0
      else
        ex=1; ey=0; ez=0
      endif
!
!  k.e
!
      kde=kx*ex+ky*ey+kz*ez
!
!  k x e
!
      kex=ky*ez-kz*ey
      key=kz*ex-kx*ez
      kez=kx*ey-ky*ex
!
!  k x (k x e)
!
      kkex=ky*kez-kz*key
      kkey=kz*kex-kx*kez
      kkez=kx*key-ky*kex
!
!  ik x (k x e) + i*phase
!
!  Normalise ff; since we don't know dt yet, we finalize this
!  within timestep where dt is determined and broadcast.
!
!  This does already include the new sqrt(2) factor (missing in B01).
!  So, in order to reproduce the 0.1 factor mentioned in B01
!  we have to set force=0.07.
!
      ffnorm=sqrt(2.)*k*sqrt(k2-kde**2)/sqrt(kav*cs0**3)
      if (ip.le.12) print*,'k,kde,ffnorm,kav,dt,cs0=',k,kde,ffnorm,kav,dt,cs0
      if (ip.le.12) print*,'k*sqrt(k2-kde**2)=',k*sqrt(k2-kde**2)
      write(21,'(f10.4,3i3,f7.3)') t,kx,ky,kz,phase
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      fact=force/ffnorm*sqrt(dt)
!
!  The wavevector is for the case where Lx=Ly=Lz=2pi. If that is not the
!  case one needs to scale by 2pi/Lx, etc.
!
      fx=exp(cmplx(0.,kx*x+phase))*fact
      fy=exp(cmplx(0.,ky*y))
      fz=exp(cmplx(0.,kz*z))
!
!  add possibility of z-profile in the forcing
!
      if (height/=0.) then
        tmpz=(z/height)**2
        fz=fz*exp(-tmpz**5/amax1(1.-tmpz,1e-5))
      endif
!
      if (ip.le.5) print*,'fx=',fx
      if (ip.le.5) print*,'fy=',fy
      if (ip.le.5) print*,'fz=',fz
!
!  prefactor
!
      sig=relhel
      coef(1)=cmplx(k*float(kex),sig*float(kkex))
      coef(2)=cmplx(k*float(key),sig*float(kkey))
      coef(3)=cmplx(k*float(kez),sig*float(kkez))
      if (ip.le.5) print*,'coef=',coef
!
      do j=1,3
        jf=j+iux-1
        do n=n1,n2
        do m=m1,m2
          f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+real(coef(j)*fx(l1:l2)*fy(m)*fz(n))
        enddo
        enddo
      enddo
!
      if (ip.le.12) print*,'forcing OK'
!
    endsubroutine forcing_hel
!***********************************************************************
    subroutine forcing_roberts(f)
!
!  add some artificial fountain flow
!  (to check for example small scale magnetic helicity loss)
!
!  30-may-02/axel: coded
!
      use Mpicomm
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx) :: sxx,cxx
      real, dimension (mx) :: sx,cx
      real, dimension (my) :: sy,cy
      real, dimension (mz) :: sz,cz,tmpz,gz,gg,ss=1.,gz1
      real :: kx,ky,kz,ffnorm,fac
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      kx=kfountain
      ky=kfountain
      kz=1.
!
      sx=sin(kx*x); cx=cos(kx*x)
      sy=sin(ky*y); cy=cos(ky*y)
      sz=sin(kz*z); cz=cos(kz*z)
!
!  abbreviation
!
      sxx=sx(l1:l2)
      cxx=cx(l1:l2)
!
!  g(z) and g'(z)
!  use z-profile to cut off
!
      if (height/=0.) then
        tmpz=(z/height)**2
        gz=sz*exp(-tmpz**5/amax1(1.-tmpz,1e-5))
      endif
!
      fac=1./(60.*dz)
      gg(1:3)=0.; gg(mz-2:mz)=0. !!(border pts to zero)
      gg(4:mz-3)=fac*(45.*(gz(5:mz-2)-gz(3:mz-4)) &
                      -9.*(gz(6:mz-1)-gz(2:mz-5)) &
                         +(gz(7:mz)  -gz(1:mz-6)))
!
!  make sign antisymmetric
!
      where(z<0) ss=-1.
      gz1=-ss*gz !!(negative for z>0)
      ffnorm=fountain*nu*dt
!
!  set forcing function
!
      do n=n1,n2
      do m=m1,m2
        f(l1:l2,m,n,iux)=f(l1:l2,m,n,iux)+ffnorm*(+sxx*cy(m)*gz1(n)+cxx*sy(m)*gg(n))
        f(l1:l2,m,n,iuy)=f(l1:l2,m,n,iuy)+ffnorm*(-cxx*sy(m)*gz1(n)+sxx*cy(m)*gg(n))
        f(l1:l2,m,n,iuz)=f(l1:l2,m,n,iuz)+ffnorm*sxx*sy(m)*gz(n)*2.
      enddo
      enddo
!
    endsubroutine forcing_roberts
!***********************************************************************
    subroutine forcing_fountain(f)
!
!  add some artificial fountain flow
!  (to check for example small scale magnetic helicity loss)
!
!  30-may-02/axel: coded
!
      use Mpicomm
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx) :: sxx,cxx
      real, dimension (mx) :: sx,cx
      real, dimension (my) :: sy,cy
      real, dimension (mz) :: sz,cz,tmpz,gz,gg,ss=1.,gz1
      real :: kx,ky,kz,ffnorm,fac
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      kx=kfountain
      ky=kfountain
      kz=1.
!
      sx=sin(kx*x); cx=cos(kx*x)
      sy=sin(ky*y); cy=cos(ky*y)
      sz=sin(kz*z); cz=cos(kz*z)
!
!  abbreviation
!
      sxx=sx(l1:l2)
      cxx=cx(l1:l2)
!
!  g(z) and g'(z)
!  use z-profile to cut off
!
      if (height/=0.) then
        tmpz=(z/height)**2
        gz=sz*exp(-tmpz**5/amax1(1.-tmpz,1e-5))
      endif
!
      fac=1./(60.*dz)
      gg(1:3)=0.; gg(mz-2:mz)=0. !!(border pts to zero)
      gg(4:mz-3)=fac*(45.*(gz(5:mz-2)-gz(3:mz-4)) &
                      -9.*(gz(6:mz-1)-gz(2:mz-5)) &
                         +(gz(7:mz)  -gz(1:mz-6)))
!
!  make sign antisymmetric
!
      where(z<0) ss=-1.
      gz1=-ss*gz !!(negative for z>0)
      ffnorm=fountain*nu*kfountain**2*dt
!
!  set forcing function
!
      do n=n1,n2
      do m=m1,m2
        f(l1:l2,m,n,iux)=f(l1:l2,m,n,iux)+ffnorm*(cxx*sy(m)*gg(n))
        f(l1:l2,m,n,iuy)=f(l1:l2,m,n,iuy)+ffnorm*(sxx*cy(m)*gg(n))
        f(l1:l2,m,n,iuz)=f(l1:l2,m,n,iuz)+ffnorm*sxx*sy(m)*gz(n)*2.
      enddo
      enddo
!
    endsubroutine forcing_fountain
!***********************************************************************

endmodule Forcing
