! $Id: forcing.f90,v 1.50 2003-08-13 15:30:07 mee Exp $

module Forcing

!  Module for forcing in the Navier-Stokes equation
!  (or, in special cases, in the entropy equation).

  use Cdata
  use General

  implicit none

  real :: force=0.,force2=0.
  real :: relhel=1.,height_ff=0.,r_ff=0.,fountain=1.,width_ff=.5
  real :: dforce=0.,radius_ff,k1_ff=1.,slope_ff=0.,work_ff=0.
  real :: tforce_stop=impossible
  integer :: kfountain=5
  logical :: lwork_ff=.false.,lmomentum_ff=.false.
  character (len=labellen) :: iforce='zero', iforce2='zero'

  integer :: dummy              ! We cannot define empty namelists
  namelist /forcing_init_pars/ dummy

  namelist /forcing_run_pars/ &
       iforce,force,relhel,height_ff,r_ff,width_ff, &
       iforce2,force2,kfountain,fountain,tforce_stop, &
       dforce,radius_ff,k1_ff,slope_ff,work_ff,lmomentum_ff

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
      logical :: lstart
!
      if (.not. first) call stop_it('register_forcing: called twice')
      first = .false.
!
      lforcing = .true.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: forcing.f90,v 1.50 2003-08-13 15:30:07 mee Exp $")
!
    endsubroutine register_forcing
!***********************************************************************
    subroutine initialize_forcing(lstart)
!
!  read seed field parameters
!  nothing done from start.f90 (lstart=.true.)
!
      use Cdata
      use Sub, only: inpui
!
      logical, save :: first=.true.
      logical :: lstart
!
      if (lstart) then
        if(ip<4) print*,'initialize_forcing: not needed in start'
      else
        if (first) then
           if (lroot.and.ip<14) print*, 'initialize_forcing: reading seed file'
           call inpui(trim(directory)//'/seed.dat',seed,nseed)
           call random_seed_wrapper(put=seed(1:nseed))
        endif
!
!  check whether we want constant forcing at each timestep,
!  in which case lwork_ff is set to true.
!
        if(work_ff/=0.) then
          force=1.
          lwork_ff=.true.
          if(lroot) print*,'initialize_forcing: reset force=1., because work_ff is set'
        endif
      endif
!
    endsubroutine initialize_forcing
!***********************************************************************
    subroutine addforce(f)
!
!  add forcing at the end of each time step
!  Since forcing is constant during one time step,
!  this can be added as an Euler 1st order step
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
!  if t>tforce_stop, then turn off forcing
!  This can be useful for producing good initial conditions
!  for turbulent decay experiments.
!
      if (t>tforce_stop) then
        if (headtt.or.ldebug) print*,'addforce: t>tforce_stop; no forcing'
      else
        if (headtt.or.ldebug) print*,'addforce: addforce started'
!
!  calculate and add forcing function
!
        select case(iforce)
        case ('zero'); if (headt) print*,'addforce: No forcing'
        case ('irrotational');  call forcing_irro(f)
        case ('helical', '2');  call forcing_hel(f)
        case ('fountain', '3'); call forcing_fountain(f)
        case ('horiz-shear');   call forcing_hshear(f)
        case ('twist');         call forcing_twist(f)
        case ('diffrot');       call forcing_diffrot(f,force)
        case ('blobs');         call forcing_blobs(f)
        case ('hel_smooth');    call forcing_hel_smooth(f)
        case default; if(lroot) print*,'addforce: No such forcing iforce=',trim(iforce)
        endselect
!
!  add *additional* forcing function
!
        select case(iforce2)
        case ('zero'); if(headtt .and. lroot) print*,'addforce: No additional forcing'
        case ('irrotational'); call forcing_irro(f)
        case ('helical');      call forcing_hel(f)
        case ('fountain');     call forcing_fountain(f)
        case ('horiz-shear');  call forcing_hshear(f)
        case ('diffrot');      call forcing_diffrot(f,force2)
        case default; if(lroot) print*,'addforce: No such forcing iforce2=',trim(iforce2)
        endselect
!
        if (headtt.or.ldebug) print*,'addforce: done addforce'
      endif
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
        if (lroot) print*,'forcing_irro: opening k.dat'
        open(9,file='k.dat')
        read(9,*) nk,kav
        if (lroot) print*,'forcing_irro: average k=',kav
        if(nk.gt.mk) then
          if (lroot) print*,'forcing_irro: dimension mk in forcing_irro is insufficient'
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
      call random_number_wrapper(fran)
      phase=pi*(2*fran(1)-1.)
      ik=nk*.9999*fran(2)+1
      if (ip<=6) print*,'forcing_irro: ik,phase,kk=',ik,phase,kkx(ik),kky(ik),kkz(ik),dt,ifirst
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
!  Add helical forcing function, using a set of precomputed wavevectors.
!  The relative helicity of the forcing function is determined by the factor
!  sigma, called here also relhel. If it is +1 or -1, the forcing is a fully
!  helical Beltrami wave of positive or negative helicity. For |relhel| < 1
!  the helicity less than maximum. For relhel=0 the forcing is nonhelical.
!  The forcing function is now normalized to unity (also for |relhel| < 1).
!
!  10-apr-00/axel: coded
!   3-sep-02/axel: introduced k1_ff, to rescale forcing function if k1/=1.
!  25-sep-02/axel: preset force_ampl to unity (in case slope is not controlled)
!   9-nov-02/axel: corrected normalization factor for the case |relhel| < 1.
!
      use Mpicomm
      use Cdata
      use General
      use Sub
      use Hydro
!
      real :: phase,ffnorm
      real, save :: kav
      real, dimension (2) :: fran
      real, dimension (nx) :: radius,tmpx,rho1
      real, dimension (mz) :: tmpz
      real, dimension (mx,my,mz,mvar+maux) :: f
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      complex, dimension (3) :: coef
      logical, dimension (3), save :: extent
      integer, parameter :: mk=3000
      integer, dimension(mk), save :: kkx,kky,kkz
      integer, save :: ifirst,nk
      integer :: ik,j,jf
      real :: kx0,kx,ky,kz,k2,k,force_ampl=1.
      real :: ex,ey,ez,kde,sig=1.,fact,kex,key,kez,kkex,kkey,kkez
!
      if (ifirst==0) then
        if (lroot) print*,'forcing_hel: opening k.dat'
        open(9,file='k.dat')
        read(9,*) nk,kav
        if (lroot) print*,'forcing_hel: average k=',kav
        if(nk.gt.mk) then
          if (lroot) print*,'forcing_hel: dimension mk in forcing_hel is insufficient'
          print*,'nk=',nk,'mk=',mk
          call mpifinalize
        end if
        read(9,*) (kkx(ik),ik=1,nk)
        read(9,*) (kky(ik),ik=1,nk)
        read(9,*) (kkz(ik),ik=1,nk)
        close(9)
        extent(1)=nx.ne.1
        extent(2)=ny.ne.1
        extent(3)=nz.ne.1
      endif
      ifirst=ifirst+1
!
!  generate random coefficients -1 < fran < 1
!  ff=force*Re(exp(i(kx+phase)))
!  |k_i| < akmax
!
      call random_number_wrapper(fran)
      phase=pi*(2*fran(1)-1.)
      ik=nk*.9999*fran(2)+1
      if(ip<=6) print*,'forcing_hel: ik,phase,kk=',ik,phase,kkx(ik),kky(ik),kkz(ik),dt,ifirst
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
      if(headt.or.ip<5) print*, 'forcing_hel: kx0,kx,ky,kz=',kx0,kx,ky,kz
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
!  Normalize ff; since we don't know dt yet, we finalize this
!  within timestep where dt is determined and broadcast.
!
!  This does already include the new sqrt(2) factor (missing in B01).
!  So, in order to reproduce the 0.1 factor mentioned in B01
!  we have to set force=0.07.
!
!  Furthermore, for |relhel| < 1, sqrt(2) should be replaced by
!  sqrt(1.+relhel**2). This is done now (9-nov-02).
!  This means that the previous value of force=0.07 (for relhel=0)
!  should now be replaced by 0.05.
!
!  Note: kav is not to be scaled with k1_ff (forcing should remain
!  unaffected when changing k1_ff).
!
      ffnorm=sqrt(1.+relhel**2) &
        *k*sqrt(k2-kde**2)/sqrt(kav*cs0**3)*(k/kav)**slope_ff
      if (ip.le.9) print*,'forcing_hel: k,kde,ffnorm,kav,dt,cs0=',k,kde,ffnorm,kav,dt,cs0
      if (ip.le.9) print*,'forcing_hel: k*sqrt(k2-kde**2)=',k*sqrt(k2-kde**2)
      !!(debug...) write(21,'(f10.4,5f8.2)') t,kx0,kx,ky,kz,phase
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
      fx=exp(cmplx(0.,kx*k1_ff*x+phase))*fact
      fy=exp(cmplx(0.,ky*k1_ff*y))
      fz=exp(cmplx(0.,kz*k1_ff*z))
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
      if (ip.le.5) print*,'forcing_hel: fx=',fx
      if (ip.le.5) print*,'forcing_hel: fy=',fy
      if (ip.le.5) print*,'forcing_hel: fz=',fz
!
!  prefactor
!
      coef(1)=cmplx(k*kex,relhel*kkex)
      coef(2)=cmplx(k*key,relhel*kkey)
      coef(3)=cmplx(k*kez,relhel*kkez)
      if (ip.le.5) print*,'forcing_hel: coef=',coef
!
!  In the past we always forced the du/dt, but in some cases
!  it may be better to force rho*du/dt (if lmomentum_ff=.true.)
!  For compatibility with earlier results, lmomentum_ff=.false. by default.
!
      if(lmomentum_ff) then
        rho1=exp(-f(l1:l2,m,n,ilnrho))
      else
        rho1=1.
      endif
!
!  loop the two cases separately, so we don't check for r_ff during
!  each loop cycle which could inhibit (pseudo-)vectorisation
!
      if (r_ff == 0) then       ! no radial profile
        if (lwork_ff) call calc_force_ampl(f,fx,fy,fz,coef,force_ampl)
        do j=1,3
          if(extent(j)) then
            jf=j+iux-1
            do n=n1,n2
              do m=m1,m2
                f(l1:l2,m,n,jf) = f(l1:l2,m,n,jf) &
                  +force_ampl*rho1*real(coef(j)*fx(l1:l2)*fy(m)*fz(n))
              enddo
            enddo
          endif
        enddo
      else                      ! with radial profile
        do j=1,3
          if(extent(j)) then
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
                     f(l1:l2,m,n,jf) + rho1*real(coef(j)*tmpx*fx(l1:l2)*fy(m)*fz(n))
              enddo
            enddo
          endif
        enddo
      endif
!
      if (ip.le.9) print*,'forcing_hel: forcing OK'
!
    endsubroutine forcing_hel
!***********************************************************************
    subroutine calc_force_ampl(f,fx,fy,fz,coef,force_ampl)
!
!  calculates the coefficient for a forcing that satisfies
!  <rho*u*f> = constant.
!
!   7-sep-02/axel: coded
!
      use Cdata
      use Sub
      use Hydro
      use Mpicomm
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx,3) :: uu
      real, dimension (nx) :: rho,udotf
      real, dimension (1) :: fsum_tmp,fsum
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      complex, dimension (3) :: coef
      real :: rho_uu_ff,force_ampl
      integer :: j
!
      rho_uu_ff=0.
      do n=n1,n2
        do m=m1,m2
          rho=exp(f(l1:l2,m,n,ilnrho))
          uu=f(l1:l2,m,n,iux:iuz)
          udotf=0.
          do j=1,3
            udotf=udotf+uu(:,j)*real(coef(j)*fx(l1:l2)*fy(m)*fz(n))
          enddo
          rho_uu_ff=rho_uu_ff+sum(rho*udotf)
        enddo
      enddo
!
!  on different processors, this result needs to be communicated
!  to other processors
!
      fsum_tmp(1)=rho_uu_ff
      call mpireduce_sum(fsum_tmp,fsum,1)
      if(lroot) rho_uu_ff=rho_uu_ff/(nw*ncpus)
      call mpibcast_real(rho_uu_ff,1)
!
!  scale forcing function
!  but do this only when rho_uu_ff>0.; never allow it to change sign
!
      if(headt) print*,'calc_force_ampl: divide forcing function by rho_uu_ff=',rho_uu_ff
      force_ampl=work_ff/(.1+amax1(0.,rho_uu_ff))
!
    endsubroutine calc_force_ampl
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
        if (lroot) print*,'force_hel_noshear: opening k.dat'
        open(9,file='k.dat')
        read(9,*) nk,kav
        if (lroot) print*,'force_hel_noshear: average k=',kav
        if(nk.gt.mk) then
          if (lroot) print*,'force_hel_noshear: dimension mk in forcing_hel is insufficient'
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
      call random_number_wrapper(fran)
      phase=pi*(2*fran(1)-1.)
      ik=nk*.9999*fran(2)+1
      if (ip<=6) print*,'force_hel_noshear: ik,phase,kk=',ik,phase,kkx(ik),kky(ik),kkz(ik),dt,ifirst
!
      kx=kkx(ik)
      ky=kky(ik)
      kz=kkz(ik)
      if(ip.le.4) print*, 'force_hel_noshear: kx,ky,kz=',kx,ky,kz
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
      if (ip.le.12) print*,'force_hel_noshear: k,kde,ffnorm,kav,dt,cs0=',k,kde,ffnorm,kav,dt,cs0
      if (ip.le.12) print*,'force_hel_noshear: k*sqrt(k2-kde**2)=',k*sqrt(k2-kde**2)
      !!(debug...) write(21,'(f10.4,3i3,f7.3)') t,kx,ky,kz,phase
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
        if (lroot .and. ifirst==1) print*,'forcing_hel_noshear: include z-profile'
        tmpz=(z/height_ff)**2
        fz=fz*exp(-tmpz**5/amax1(1.-tmpz,1e-5))
      endif
!
!  possibly multiply forcing by sgn(z) and radial profile
!
      if (r_ff/=0.) then
        if (lroot .and. ifirst==1) &
             print*,'forcing_hel_noshear: applying sgn(z)*xi(r) profile'
        !
        ! only z-dependent part can be done here; radial stuff needs to go
        ! into the loop
        !
        tmpz = tanh(z/width_ff)
        fz = fz*tmpz
      endif
!
      if (ip.le.5) print*,'force_hel_noshear: fx=',fx
      if (ip.le.5) print*,'force_hel_noshear: fy=',fy
      if (ip.le.5) print*,'force_hel_noshear: fz=',fz
!
!  prefactor
!
      sig=relhel
      coef(1)=cmplx(k*float(kex),sig*float(kkex))
      coef(2)=cmplx(k*float(key),sig*float(kkey))
      coef(3)=cmplx(k*float(kez),sig*float(kkez))
      if (ip.le.5) print*,'force_hel_noshear: coef=',coef
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
      if (ip.le.12) print*,'force_hel_noshear: forcing OK'
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: sxx,cxx
      real, dimension (mx) :: sx,cx
      real, dimension (my) :: sy,cy
      real, dimension (mz) :: sz,cz,tmpz,gz,gg,ss=1.,gz1
      real :: kx,ky,kz,ffnorm,fac
!
!  identify ourselves
!
      if (headtt.or.ldebug) print*,'forcing_roberts: ENTER'
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: sxx,cxx
      real, dimension (mx) :: sx,cx
      real, dimension (my) :: sy,cy
      real, dimension (mz) :: sz,cz,tmpz,gz,gg,ss=1.,gz1
      real :: kx,ky,kz,ffnorm,fac
!
!  identify ourselves
!
      if (headtt.or.ldebug) print*,'forcing_fountain: ENTER'
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: fx
      real, dimension (mz) :: fz
      real :: kx,ffnorm
!
!  need to multiply by dt (for Euler step).
!  Define fz with ghost zones, so fz(n) is the correct position
!  Comment: Brummell et al have a polynomial instead!
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
    subroutine forcing_diffrot(f,force_ampl)
!
!  add differential rotation
!
!  26-jul-02/axel: coded
!
      use Mpicomm
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx,nz) :: fx,fz,tmp
      real :: force_ampl,ffnorm,ffnorm2,kx,kz
!
!  identifier
!
      if(headt) print*,'forcing_diffrot: ENTER'
!
!  need to multiply by dt (for Euler step).
!
      ffnorm=force_ampl*dt  !(dt for the timestep)
!
!  prepare velocity, Uy=cosx*cosz
!
      fx=spread(cos(x(l1:l2)),2,nz)
      fz=spread(cos(z(n1:n2)),1,nx)
!
!  this forcing term is balanced by diffusion operator;
!  need to multiply by nu*k^2, but k=sqrt(1+1) for the forcing
!
      ffnorm2=ffnorm*nu*2
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, save :: tforce=0.
      integer, save :: ifirst=0
      integer, save :: nforce=0
      logical :: lforce
      character (len=4) :: ch
      character (len=135) :: file
!
!  identifier
!
      if(headt) print*,'forcing_blobs: ENTER'
!
!  the last forcing time is recorded in tforce.dat
!
      file=trim(datadir)//'/tforce.dat'
      if (ifirst==0) then
        call read_snaptime(trim(file),tforce,nforce,dforce,t)
        ifirst=1
      endif
!
!  Check whether we want to do forcing at this time.
!
      call update_snaptime(file,tforce,nforce,dforce,t,lforce,ch,ENUM=.true.)
      if (lforce) then
        call blob(force,f,iss,radius_ff,0.,0.,.5)
      endif
!
    endsubroutine forcing_blobs
!***********************************************************************
    subroutine forcing_hel_smooth(f)
!
      use Mpicomm
      use Cdata
      use Hydro
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,3) :: force1,force2,force_vec
      real :: phase1,phase2,p_weight
      real :: kx01,ky1,kz1,kx02,ky2,kz2
      integer, parameter :: mk=3000
      integer, dimension(mk), save :: kkx,kky,kkz
      integer, save :: ifirst,nk
      integer :: ik1,ik2,ik
      real, save :: kav
!
      if (ifirst==0) then
         if (lroot) print*,'forcing_hel_smooth: opening k.dat'
         open(9,file='k.dat')
         read(9,*) nk,kav
         if (lroot) print*,'forcing_hel_smooth: average k=',kav
         if(nk.gt.mk) then
            if (lroot) print*, &
                 'forcing_hel_smooth: dimension mk in forcing_hel_smooth is insufficient'
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
!  Re-calculate forcing wave numbers if necessary
!
      !tsforce is set to -10 in cdata.f90. It should also be saved in a file
      !so that it will be read again on restarts.
      if (t .gt. tsforce) then  
         if (tsforce .lt. 0) then
            call random_number_wrapper(fran1)
         else
            fran1=fran2
         endif
         call random_number_wrapper(fran2)
         tsforce=t+dtforce
      endif
      phase1=pi*(2*fran1(1)-1.)
      ik1=nk*.9999*fran1(2)+1
      kx01=kkx(ik1)
      ky1=kky(ik1)
      kz1=kkz(ik1)
      phase2=pi*(2*fran2(1)-1.)
      ik2=nk*.9999*fran2(2)+1
      kx02=kkx(ik2)
      ky2=kky(ik2)
      kz2=kkz(ik2)

!
!  Calculate forcing function
!
      call hel_vec(f,kx01,ky1,kz1,phase1,kav,ifirst,force1)
      call hel_vec(f,kx02,ky2,kz2,phase2,kav,ifirst,force2)
!
!  Determine weight parameter
!
      p_weight=(tsforce-t)/dtforce
!
!  Add forcing
!      
      force_vec=p_weight*force1+(1-p_weight)*force2
      f(l1:l2,m1:m2,n1:n2,1:3)= &
           f(l1:l2,m1:m2,n1:n2,1:3)+force_vec(l1:l2,m1:m2,n1:n2,:)
!
    end subroutine forcing_hel_smooth
!***********************************************************************
    subroutine hel_vec(f,kx0,ky,kz,phase,kav,ifirst,force1)
!
!  Add helical forcing function, using a set of precomputed wavevectors.
!  The relative helicity of the forcing function is determined by the factor
!  sigma, called here also relhel. If it is +1 or -1, the forcing is a fully
!  helical Beltrami wave of positive or negative helicity. For |relhel| < 1
!  the helicity less than maximum. For relhel=0 the forcing is nonhelical.
!  The forcing function is now normalized to unity (also for |relhel| < 1).
!
!  10-apr-00/axel: coded
!   3-sep-02/axel: introduced k1_ff, to rescale forcing function if k1/=1.
!  25-sep-02/axel: preset force_ampl to unity (in case slope is not controlled)
!   9-nov-02/axel: corrected normalization factor for the case |relhel| < 1.
!  17-jan-03/nils: adapted from forcing_hel
!
      use Mpicomm
      use Cdata
      use General
      use Sub
      use Hydro
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: phase,ffnorm
      real :: kav
      real, dimension (nx) :: radius,tmpx
      real, dimension (mz) :: tmpz
      real, dimension (mx,my,mz,3) :: force1
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      complex, dimension (3) :: coef
      integer :: j,jf
      integer :: ifirst
      real :: kx0,kx,ky,kz,k2,k,force_ampl=1.
      real :: ex,ey,ez,kde,sig=1.,fact,kex,key,kez,kkex,kkey,kkez
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
      if(headt.or.ip<5) print*, 'hel_vec: kx0,kx,ky,kz=',kx0,kx,ky,kz
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
!  Normalize ff; since we don't know dt yet, we finalize this
!  within timestep where dt is determined and broadcast.
!
!  This does already include the new sqrt(2) factor (missing in B01).
!  So, in order to reproduce the 0.1 factor mentioned in B01
!  we have to set force=0.07.
!
!  Furthermore, for |relhel| < 1, sqrt(2) should be replaced by
!  sqrt(1.+relhel**2). This is done now (9-nov-02).
!  This means that the previous value of force=0.07 (for relhel=0)
!  should now be replaced by 0.05.
!
!  Note: kav is not to be scaled with k1_ff (forcing should remain
!  unaffected when changing k1_ff).
!
      ffnorm=sqrt(1.+relhel**2) &
        *k*sqrt(k2-kde**2)/sqrt(kav*cs0**3)*(k/kav)**slope_ff
      if (ip.le.9) print*,'hel_vec: k,kde,ffnorm,kav,dt,cs0=',k,kde, &
                                                        ffnorm,kav,dt,cs0
      if (ip.le.9) print*,'hel_vec: k*sqrt(k2-kde**2)=',k*sqrt(k2-kde**2)
      !!(debug...) write(21,'(f10.4,5f8.2)') t,kx0,kx,ky,kz,phase
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
      fx=exp(cmplx(0.,kx*k1_ff*x+phase))*fact
      fy=exp(cmplx(0.,ky*k1_ff*y))
      fz=exp(cmplx(0.,kz*k1_ff*z))
!
!  possibly multiply forcing by z-profile
!
      if (height_ff/=0.) then
        if (lroot .and. ifirst==1) print*,'hel_vec: include z-profile'
        tmpz=(z/height_ff)**2
        fz=fz*exp(-tmpz**5/amax1(1.-tmpz,1e-5))
      endif
!
!  possibly multiply forcing by sgn(z) and radial profile
!
      if (r_ff/=0.) then
        if (lroot .and. ifirst==1) &
             print*,'hel_vec: applying sgn(z)*xi(r) profile'
        !
        ! only z-dependent part can be done here; radial stuff needs to go
        ! into the loop
        !
        tmpz = tanh(z/width_ff)
        fz = fz*tmpz
      endif
!
      if (ip.le.5) print*,'hel_vec: fx=',fx
      if (ip.le.5) print*,'hel_vec: fy=',fy
      if (ip.le.5) print*,'hel_vec: fz=',fz
!
!  prefactor
!
      coef(1)=cmplx(k*kex,relhel*kkex)
      coef(2)=cmplx(k*key,relhel*kkey)
      coef(3)=cmplx(k*kez,relhel*kkez)
      if (ip.le.5) print*,'hel_vec: coef=',coef
!
! loop the two cases separately, so we don't check for r_ff during
! each loop cycle which could inhibit (pseudo-)vectorisation
!
      if (r_ff == 0) then       ! no radial profile
        if (lwork_ff) call calc_force_ampl(f,fx,fy,fz,coef,force_ampl)
        do j=1,3
          jf=j+iux-1
          do n=n1,n2
            do m=m1,m2
               force1(l1:l2,m,n,jf) = &
                +force_ampl*real(coef(j)*fx(l1:l2)*fy(m)*fz(n))
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
              force1(l1:l2,m,n,jf) =  real(coef(j)*tmpx*fx(l1:l2)*fy(m)*fz(n))
            enddo
          enddo
        enddo
      endif
!
      if (ip.le.9) print*,'hel_vec: forcing OK'
!
    end subroutine hel_vec
!***********************************************************************

endmodule Forcing
