module Forcing

!! Module for forcing in the Navier-Stokes equation.

  use Cdata

  implicit none

  namelist /forcing_run_pars/ &
       iforce,force,relhel

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
           "$Revision: 1.6 $", &
           "$Date: 2002-05-26 16:42:58 $")
!
    endsubroutine register_forcing
!***********************************************************************
    subroutine addforce(f)
!
      use Cdata
!
!  add forcing at the end of each time step
!  Since forcing is constant during one time step,
!  this can be added as an Euler 1st order step
!
      real, dimension (mx,my,mz,mvar) :: f
!
!  calculate and add forcing function
!
      if(iforce==1) call forcing1(f)
      if(iforce==2) call forcing2(f)
!
    endsubroutine addforce
!***********************************************************************
    subroutine forcing1(f)
!
!  add acoustic forcing function, using a set of precomputed wavevectors
!  This forcing drives pressure waves
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
    endsubroutine forcing1
!***********************************************************************
    subroutine forcing2(f)
!
!  add helical forcing function, using a set of precomputed wavevectors
!
      use Mpicomm
      use Cdata
      use Sub
!
      real :: phase,ffnorm
      real, save :: kav
      real, dimension (2) :: fran
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
    endsubroutine forcing2
!***********************************************************************

endmodule Forcing
