! $Id: forcing.f90,v 1.24 2002-08-18 12:04:40 brandenb Exp $

module Forcing

!  Module for forcing in the Navier-Stokes equation
!  (or, in special cases, in the entropy equation).

  use Cdata

  implicit none

  real :: force=0.,relhel=1.,height_ff=0.,r_ff=0.,fountain=1.,width_ff=.5
  real :: dforce=0.,radius_ff
  integer :: kfountain=5
  character (len=labellen) :: iforce='zero', iforce2='zero'

  integer :: dummy              ! We cannot define empty namelists
  namelist /forcing_init_pars/ dummy

  namelist /forcing_run_pars/ &
       iforce,force,relhel,height_ff,r_ff,width_ff, &
       iforce2,kfountain,fountain, &
       dforce,radius_ff

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
           "$Id: forcing.f90,v 1.24 2002-08-18 12:04:40 brandenb Exp $")
!
    endsubroutine register_forcing
!***********************************************************************
    subroutine forcing_run_hook()
!
!  read seed field parameters
!
      use Cdata
      use Sub, only: inpui
!
      if (lroot.and.ip<14) print*, 'reading seed file'
      call inpui(trim(directory)//'/seed.dat',seed,nseed)
      call random_seed(put=seed(1:nseed))
!
    endsubroutine forcing_run_hook
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
      if (headtt.or.ldebug) print*,'FORCING: addforce started'
!
!  calculate and add forcing function
!
      select case(iforce)
      case ('zero'); if (headt) print*,'No forcing'
      case ('irrotational');  call forcing_irro(f)
      case ('helical', '2');  call forcing_hel(f)
      case ('fountain', '3'); call forcing_fountain(f)
      case ('horiz-shear');   call forcing_hshear(f)
      case ('twist');         call forcing_twist(f)
      case ('diffrot');       call forcing_diffrot(f)
      case ('blobs');         call forcing_blobs(f)
      case default; if(lroot) print*,'No such forcing iforce=',trim(iforce)
      endselect
!
!  add *additional* forcing function
!
      select case(iforce2)
      case ('zero'); if(headtt .and. lroot) print*,'No additional forcing'
      case ('irrotational'); call forcing_irro(f)
      case ('helical');      call forcing_hel(f)
      case ('fountain');     call forcing_fountain(f)
      case ('horiz-shear');  call forcing_hshear(f)
      case default; if(lroot) print*,'No such forcing iforce2=',trim(iforce2)
      endselect
!
      if (headtt.or.ldebug) print*,'FORCING: done addforce'
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
      integer, parameter :: mk=3000
      integer, dimension(mk), save :: kkx,kky,kkz
      integer, save :: ifirst,nk
      integer :: ik,j,jf
!
      if (ifirst==0) then
        if (lroot) print*,'irrotational forcing; opening k.dat'
        open(9,file='k.dat')
        read(9,*) nk,kav
        if (lroot) print*,'average k=',kav
        if(nk.gt.mk) then
          if (lroot) print*,'dimension mk in forcing_irro is insufficient'
          print*,'nk=',nk,'mk=',mk
          call mpifinalize
        end if
        read(9,*) (kkx(ik),ik=1,nk)
        read(9,*) (kky(ik),ik=1,nk)
        read(9,*) (kkz(ik),ik=1,nk)
        close(9)
      endif
      ifirst=ifirst+1
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
      real, dimension (nx) :: radius,tmpx
      real, dimension (mz) :: tmpz
      real, dimension (mx,my,mz,mvar) :: f
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      complex, dimension (3) :: coef
      integer, parameter :: mk=3000
      integer, dimension(mk), save :: kkx,kky,kkz
      integer, save :: ifirst,nk
      integer :: ik,j,jf
      real :: kx0,kx,ky,kz,k2,k
      real :: ex,ey,ez,kde,sig=1.,fact,kex,key,kez,kkex,kkey,kkez
!
      if (ifirst==0) then
        if (lroot) print*,'helical forcing; opening k.dat'
        open(9,file='k.dat')
        read(9,*) nk,kav
        if (lroot) print*,'average k=',kav
        if(nk.gt.mk) then
          if (lroot) print*,'dimension mk in forcing_hel is insufficient'
          print*,'nk=',nk,'mk=',mk
          call mpifinalize
        end if
        read(9,*) (kkx(ik),ik=1,nk)
        read(9,*) (kky(ik),ik=1,nk)
        read(9,*) (kkz(ik),ik=1,nk)
        close(9)
      endif
      ifirst=ifirst+1
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
      kx0=kkx(ik)
      ky=kky(ik)
      kz=kkz(ik)
!
!  in the shearing sheet approximation, kx = kx0 - St*k_y.
!  Here, St=-deltay/Lx
!
      if (Sshear==0.) then
        kx=kx0
      else
        kx=kx0+ky*deltay/Lx
      endif
!
      if(headt.or.ip<5) print*, 'kx0,kx,ky,kz=',kx0,kx,ky,kz
      k2=kx**2+ky**2+kz**2
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
      write(21,'(f10.4,5f8.2)') t,kx0,kx,ky,kz,phase
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
      fx=exp(cmplx(0.,2*pi/Lx*kx*x+phase))*fact
      fy=exp(cmplx(0.,2*pi/Ly*ky*y))
      fz=exp(cmplx(0.,2*pi/Lz*kz*z))
!
!  possibly multiply forcing by z-profile
!
      if (height_ff/=0.) then
        if (lroot .and. ifirst==1) print*,'forcing_hel: include z-profile'
        tmpz=(z/height_ff)**2
        fz=fz*exp(-tmpz**5/amax1(1.-tmpz,1e-5))
      endif
!
!  possibly multiply forcing by sgn(z) and radial profile
!
      if (r_ff/=0.) then
        if (lroot .and. ifirst==1) &
             print*,'forcing_hel: applying sgn(z)*xi(r) profile'
        !
        ! only z-dependent part can be done here; radial stuff needs to go
        ! into the loop
        !
        tmpz = tanh(z/width_ff)
        fz = fz*tmpz
      endif
!
      if (ip.le.5) print*,'fx=',fx
      if (ip.le.5) print*,'fy=',fy
      if (ip.le.5) print*,'fz=',fz
!
!  prefactor
!
      sig=relhel
      coef(1)=cmplx(k*kex,sig*kkex)
      coef(2)=cmplx(k*key,sig*kkey)
      coef(3)=cmplx(k*kez,sig*kkez)
      if (ip.le.5) print*,'coef=',coef
!
! loop the two cases separately, so we don't check for r_ff during
! each loop cycle which could inhibit (pseudo-)vectorisation
!
      if (r_ff == 0) then       ! no radial profile
        do j=1,3
          jf=j+iux-1
          do n=n1,n2
            do m=m1,m2
              f(l1:l2,m,n,jf) = &
                   f(l1:l2,m,n,jf)+real(coef(j)*fx(l1:l2)*fy(m)*fz(n))
            enddo
          enddo
        enddo
      else                      ! with radial profile
        do j=1,3
          jf=j+iux-1
          do n=n1,n2
            sig = relhel*tmpz(n)
            coef(1)=cmplx(k*kex,sig*kkex)
            coef(2)=cmplx(k*key,sig*kkey)
            coef(3)=cmplx(k*kez,sig*kkez)
            do m=m1,m2
              radius = sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
              tmpx = 0.5*(1.-tanh((radius-r_ff)/width_ff))
              f(l1:l2,m,n,jf) = &
                   f(l1:l2,m,n,jf) + real(coef(j)*tmpx*fx(l1:l2)*fy(m)*fz(n))
            enddo
          enddo
        enddo
      endif
!
      if (ip.le.12) print*,'forcing OK'
!
    endsubroutine forcing_hel
!***********************************************************************
    subroutine forcing_hel_noshear(f)
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
      real, dimension (nx) :: radius,tmpx
      real, dimension (mz) :: tmpz
      real, dimension (mx,my,mz,mvar) :: f
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      complex, dimension (3) :: coef
      integer, parameter :: mk=3000
      integer, dimension(mk), save :: kkx,kky,kkz
      integer, save :: ifirst,nk
      integer :: ik,j,jf,kx,ky,kz,kex,key,kez,kkex,kkey,kkez
      real :: k2,k,ex,ey,ez,kde,sig=1.,fact
!
      if (ifirst==0) then
        if (lroot) print*,'helical forcing; opening k.dat'
        open(9,file='k.dat')
        read(9,*) nk,kav
        if (lroot) print*,'average k=',kav
        if(nk.gt.mk) then
          if (lroot) print*,'dimension mk in forcing_hel is insufficient'
          print*,'nk=',nk,'mk=',mk
          call mpifinalize
        end if
        read(9,*) (kkx(ik),ik=1,nk)
        read(9,*) (kky(ik),ik=1,nk)
        read(9,*) (kkz(ik),ik=1,nk)
        close(9)
      endif
      ifirst=ifirst+1
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
      fx=exp(cmplx(0.,2*pi/Lx*kx*x+phase))*fact
      fy=exp(cmplx(0.,2*pi/Ly*ky*y))
      fz=exp(cmplx(0.,2*pi/Lz*kz*z))
!
!  possibly multiply forcing by z-profile
!
      if (height_ff/=0.) then
        if (lroot .and. ifirst==1) print*,'forcing_hel: include z-profile'
        tmpz=(z/height_ff)**2
        fz=fz*exp(-tmpz**5/amax1(1.-tmpz,1e-5))
      endif
!
!  possibly multiply forcing by sgn(z) and radial profile
!
      if (r_ff/=0.) then
        if (lroot .and. ifirst==1) &
             print*,'forcing_hel: applying sgn(z)*xi(r) profile'
        !
        ! only z-dependent part can be done here; radial stuff needs to go
        ! into the loop
        !
        tmpz = tanh(z/width_ff)
        fz = fz*tmpz
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
! loop the two cases separately, so we don't check for r_ff during
! each loop cycle which could inhibit (pseudo-)vectorisation
!
      if (r_ff == 0) then       ! no radial profile
        do j=1,3
          jf=j+iux-1
          do n=n1,n2
            do m=m1,m2
              f(l1:l2,m,n,jf) = &
                   f(l1:l2,m,n,jf)+real(coef(j)*fx(l1:l2)*fy(m)*fz(n))
            enddo
          enddo
        enddo
      else                      ! with radial profile
        do j=1,3
          jf=j+iux-1
          do n=n1,n2
            sig = relhel*tmpz(n)
            coef(1)=cmplx(k*float(kex),sig*float(kkex))
            coef(2)=cmplx(k*float(key),sig*float(kkey))
            coef(3)=cmplx(k*float(kez),sig*float(kkez))
            do m=m1,m2
              radius = sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
              tmpx = 0.5*(1.-tanh((radius-r_ff)/width_ff))
              f(l1:l2,m,n,jf) = &
                   f(l1:l2,m,n,jf) + real(coef(j)*tmpx*fx(l1:l2)*fy(m)*fz(n))
            enddo
          enddo
        enddo
      endif
!
      if (ip.le.12) print*,'forcing OK'
!
    endsubroutine forcing_hel_noshear
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
!  identify ourselves
!
      if (headtt.or.ldebug) print*,'Roberts forcing'
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
      if (height_ff/=0.) then
        tmpz=(z/height_ff)**2
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
!  identify ourselves
!
      if (headtt.or.ldebug) print*,'fountain forcing'
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
      if (height_ff/=0.) then
        tmpz=(z/height_ff)**2
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
    subroutine forcing_hshear(f)
!
!  add horizontal shear
!
!  19-jun-02/axel+bertil: coded
!
      use Mpicomm
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx) :: fx
      real, dimension (mz) :: fz
      real :: kx,ffnorm
!
!  need to multiply by dt (for Euler step).
!  Define fz with ghost zones, so fz(n) is the correct position
!
      kx=2*pi/Lx
      fx=cos(kx*x(l1:l2))
      fz=1./cosh(z/width_ff)**2
      ffnorm=force*dt  !(dt for the timestep)
!
!  add to velocity (here only y-component)
!
      do n=n1,n2
      do m=m1,m2
        f(l1:l2,m,n,iuy)=f(l1:l2,m,n,iuy)+ffnorm*fx*fz(n)
      enddo
      enddo
!
    endsubroutine forcing_hshear
!***********************************************************************
    subroutine forcing_twist(f)
!
!  add circular twisting motion, (ux, 0, uz)
!
!  19-jul-02/axel: coded
!
      use Mpicomm
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx,nz) :: xx,zz,r2,tmp,fx,fz
      real :: ffnorm,ry2,fy,ytwist1,ytwist2
!
!  identifier
!
      if(headt) print*,'forcing_twist: r_ff,width_ff=',r_ff,width_ff
!
!  need to multiply by dt (for Euler step).
!
      ffnorm=force*dt  !(dt for the timestep)
!
!  add to velocity
!  calculate r2=(x^2+z^2)/r^2
!
      xx=spread(x(l1:l2),2,nz)
      zz=spread(z(n1:n2),1,nx)
      if (r_ff==0.) then
        if(lroot) print*,'forcing_twist: division by r_ff=0!!'
      endif
      r2=(xx**2+zz**2)/r_ff**2
      tmp=exp(-r2/amax1(1.-r2,1e-5))*ffnorm
      fx=-zz*tmp
      fz=+xx*tmp
!
!  have opposite twists at
!
      y0=xyz0(2)
      ytwist1=y0+0.25*Ly
      ytwist2=y0+0.75*Ly
!
      do m=m1,m2
        !
        ! first twister
        !
        ry2=((y(m)-ytwist1)/width_ff)**2
        fy=exp(-ry2/amax1(1.-ry2,1e-5))
        f(l1:l2,m,n1:n2,iux)=f(l1:l2,m,n1:n2,iux)+fy*fx
        f(l1:l2,m,n1:n2,iuz)=f(l1:l2,m,n1:n2,iuz)+fy*fz
        !
        ! second twister
        !
        ry2=((y(m)-ytwist2)/width_ff)**2
        fy=exp(-ry2/amax1(1.-ry2,1e-5))
        f(l1:l2,m,n1:n2,iux)=f(l1:l2,m,n1:n2,iux)-fy*fx
        f(l1:l2,m,n1:n2,iuz)=f(l1:l2,m,n1:n2,iuz)-fy*fz
      enddo
!
    endsubroutine forcing_twist
!***********************************************************************
    subroutine forcing_diffrot(f)
!
!  add differential rotation
!
!  26-jul-02/axel: coded
!
      use Mpicomm
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx,nz) :: fx,fz,tmp
      real :: ffnorm,ffnorm2,kx,kz
!
!  identifier
!
      if(headt) print*,'forcing_diffrot'
!
!  need to multiply by dt (for Euler step).
!
      ffnorm=force*dt  !(dt for the timestep)
!
!  prepare velocity, Uy=sinkx*sinkz
!
      kx=.5*pi/Lx
      kz=.5*pi/Lz
      fx=spread(sin(kx*x(l1:l2)),2,nz)
      fz=spread(sin(kz*z(n1:n2)),1,nx)
!
!  this forcing term is balanced by diffusion operator;
!  need to multiply by nu*k^2
!
      ffnorm2=ffnorm*nu*(kx**2+kz**2)
      tmp=ffnorm2*fx*fz
!
!  add
!
      do m=m1,m2
        f(l1:l2,m,n1:n2,iuy)=f(l1:l2,m,n1:n2,iuy)+tmp
      enddo
!
    endsubroutine forcing_diffrot
!***********************************************************************
    subroutine forcing_blobs(f)
!
!  add blobs in entropy every dforce time units
!
!  28-jul-02/axel: coded
!
      !use Mpicomm
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f
      real, save :: tforce=0.
      integer, save :: ifirst=0
      integer, save :: nforce=0
      logical :: lforce
      character (len=4) :: ch
      character (len=80) :: file
!
!  identifier
!
      if(headt) print*,'forcing_blobs'
!
!  the last forcing time is recorded in tforce.dat
!
      file='tmp/tforce.dat'
      if (ifirst==0) then
        call out1 (trim(file),tforce,nforce,dforce,t)
        ifirst=1
      endif
!
!  Check whether we want to do forcing at this time.
!
      call out2 (trim(file),tforce,nforce,dforce,t,lforce,ch,.true.)
      if (lforce) then
        call blob(force,f,ient,radius_ff,0.,0.,.5)
      endif
!
    endsubroutine forcing_blobs
!***********************************************************************

endmodule Forcing
