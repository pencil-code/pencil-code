! $Id: interstellar.f90,v 1.46 2003-08-19 21:39:55 mee Exp $

!  This modules contains the routines for SNe-driven ISM simulations.
!  Still in development. 

module Interstellar

  use Cparam
  use Cdata
  use Density

  implicit none

  real :: x_SN,y_SN,z_SN,rho_SN,lnrho_SN,yH_SN,TT_SN,ss_SN,ee_SN,ampl_SN=1.0
  integer :: l_SN,m_SN,n_SN
  real, dimension(nx) :: dr2_SN     ! Pencil storing radius to SN
  real :: t_next_SNI=0.0, t_interval_SNI=3.64e-3,h_SNI=0.325,h_SNII=0.09
  real :: tau_cloud=2e-2, r_SNI=3.e+4, r_SNII=4.e+3
  integer, parameter :: ninterstellarsave=1
  real, dimension(ninterstellarsave) :: interstellarsave
  real, parameter :: rho_crit=1.,TT_crit=4000.
  real, parameter :: frac_converted=0.02,frac_heavy=0.10,mass_SN=10.
  real, parameter :: rho_min=1.e-6

  ! normalisation factors for 1-d, 2-d, and 3-d profiles like exp(-r^6)
  real, parameter, dimension(3) :: &
                        cnorm_SN = (/ 1.85544 , 2.80538 , 3.71213666 /) 
  ! ( 1d: 2    int_0^infty exp(-(r/a)^6)     dr) / a
  !   2d: 2 pi int_0^infty exp(-(r/a)^6) r   dr) / a^2
  !   3d: 4 pi int_0^infty exp(-(r/a)^6) r^2 dr) / a^3 )
  ! ( cf. 3.128289613 -- from where ?!? )
  ! NB: 1d and 2d results just from numerical integration -- calculate
  !      exact integrals at some point...

!  cp1=1/cp used to convert TT (and ss) into interstellar code units
!  (useful, as many conditions conveniently expressed in terms of TT)
!  code units based on:
!    [length]  = 1kpc  = 3.09 10^21 cm
!    [time]    = 1Gyr  = 3.15 10^16 s             !no, on [u]=1km/s...
!    [rho]     =       = 1.00 10^-24 g/cm^3
!  Lambdaunits converts coolH into interstellar code units.
!   (this should really just be incorporated into coolH coefficients)
!  NB: will start using thermodynamics, and unit_length, etc., imminently...

!  real, parameter :: cp1=27.8   !=R * gamma / (mu * (gamma-1))  27.8 
  real, parameter :: TTunits=46.6
  double precision :: unit_Lambda
  real, parameter :: tosolarMkpc3=1.483e7

  real, parameter :: rhoUV_cgs=0.1
  real, parameter :: TUV_cgs=7000.,T0UV_cgs=12000.,cUV_cgs=5.e-4
  double precision, parameter, dimension(6) ::  &
  coolT_cgs=(/ 300.D0,     2000.D0,    8000.D0,    1.D5,    4.D7,     1.D9 /),  &
  coolH_cgs=(/ 2.2380D-32, 1.0012D-30, 4.6240D-36, 1.7800D-18, 3.2217D-27, 0.D0   /)
  
  real :: rhoUV,TUV,T0UV,cUV
  real, dimension(6) :: coolT, &
    coolB=(/ 2.,       1.5,      2.867,    -.65,    0.5,      0.   /)
  double precision, dimension(6) :: coolH 
  integer :: iproc_SN,ipy_SN,ipz_SN
  logical :: ltestSN = .false.  ! If set .true. SN are only exploded at the
                              ! origin and ONLY the type I scheme is used
                              ! Used in kompaneets test etc.


  ! input parameters
  real :: TT_SN_min=1.e7
  logical :: lnever_move_mass
!ajwm: disable SNII for debugging 
  logical :: lSNI=.true., lSNII=.false.
  logical :: laverageSNheating = .false.
 ! real, parameter :: TT_SN_min=1.e8    ! vary for debug, tests
 
  integer :: dummy 
  namelist /interstellar_init_pars/ dummy

  ! run parameters
  logical:: uniform_zdist_SNI = .false.
  namelist /interstellar_run_pars/ &
      t_next_SNI,t_interval_SNI,h_SNI,ampl_SN,tau_cloud, &
      uniform_zdist_SNI, ltestSN, TT_SN_min, lnever_move_mass, &
      lSNI, lSNII, laverageSNheating

  contains

!***********************************************************************
    subroutine register_interstellar()
!
!  19-nov-02/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_interstellar called twice')
      first = .false.
!
      linterstellar = .true.
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_interstellar: ENTER'
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: interstellar.f90,v 1.46 2003-08-19 21:39:55 mee Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_interstellar: nvar > mvar')
      endif
!
    endsubroutine register_interstellar
!***********************************************************************
    subroutine initialize_interstellar(lstart)
!
!  Perform any post-parameter-read initialization eg. set derived 
!  parameters
!
!  24-nov-02/tony: coded
!
!  read parameters from seed.dat and interstellar.dat
!
      use Cdata
      use General
      use Sub, only: inpui,inpup
      use Mpicomm, only: mpibcast_real
      use Ionization, only: getmu
!
      logical, save :: first=.true.
      logical :: lstart
      logical :: exist
      real :: mu

      if (first) then
         if (.not. lstart) then
            call inpui(trim(directory)//'/seed.dat',seed,nseed)
            if (lroot.and.ip<14) then
               print*, 'initialize_interstellar: reading seed file'
               print*, 'initialize_interstellar: nseed,seed',nseed,seed(1:nseed)
            endif
            call random_seed_wrapper(put=seed(1:nseed))
         endif
!
!AB: comment please why interstellar should be read in and what it does
!
         if (lroot) then
            inquire(file=trim(datadir)//'/interstellar.dat',exist=exist)
            if (exist) then 
               if (ip<=14) print*, 'initialize_interstellar: read interstellar.dat'
                call inpup(trim(datadir)//'/interstellar.dat',  &
                    interstellarsave,ninterstellarsave)
               if (ip<=14) print*, 'initialize_interstellar: t_next_SNI', &
                    interstellarsave(1)
            else
               interstellarsave(1)=t_next_SNI
            endif

         endif
         call mpibcast_real(interstellarsave,1)
         t_next_SNI=interstellarsave(1)
      endif

      if (lroot.and.uniform_zdist_SNI) then
         print*,'initialize_interstellar: using UNIFORM z-distribution of SNI'
      endif
      
      call getmu(mu) 
      if (unit_system=='cgs') then
        unit_Lambda = unit_energy * unit_velocity**3 * unit_time**2 * &
                       mu**2 * m_H**2
!unit_energy / (unit_density*unit_time*unit_mass)
      elseif (unit_system=='SI') then
        unit_Lambda = unit_energy * unit_velocity**3 * unit_time**2 * &
                       mu**2 * m_H**2 * 1D-2
!ajwm Constant factor 1D-2 IS NOT CORRECT... NEED TO RECALCULATE
      endif
      print*,'initialize_interstellar: unit_Lambda',unit_Lambda
      coolH = coolH_cgs / unit_Lambda 
      coolT = coolT_cgs / unit_temperature

      if (lroot.and.ip<14) then
        print*,'initialize_interstellar: nseed,seed',nseed,seed(1:nseed)
        print*,'initialize_interstellar: finished'
      endif
!
    endsubroutine initialize_interstellar
!***********************************************************************
    subroutine calc_heat_cool_interstellar(df,rho1,TT,TT1,yH)
!
!  This routine calculates and applies the optically thin cooling function
!  together with UV heating.
!
!  We may want to move it to the entropy module for good, because its use
!  is not restricted to interstellar runs (could be used for solar corona).
!  Also, it doesn't pose an extra load on memory usage or compile time.
!  (We should allow that UV heating can be turned off; so rhoUV should
!  be made an input parameter.)
!
!  19-nov-02/graeme: adapted from calc_heat_cool
!  10-aug-03/axel: TT is used as input
!
      use Cdata
      use Mpicomm
      use Density, only : rho0
      use Sub
      use Ionization
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      real, dimension (nx), intent(in) :: rho1,TT,TT1,yH
      real, dimension (nx) :: heat,cool,rho
      real :: norm
      integer :: i
!
!  identifier
!
      if(headtt) print*,'calc_heat_cool_interstellar: ENTER'

      rho=1./rho1
!
!  define T in K, for calculation of both UV heating and radiative cooling
!
!  add T-dept radiative cooling, from Rosen et al., ApJ, 413, 137, 1993
!  cooling is Lambda*rho^2, with (eq 7)
!     Lambda=coolH(i)*TT*coolB(i),   for coolT(i) <= T < coolT(i+1)
!  nb: our coefficients coolH(i) differ from those in Rosen et al. by
!   factor (mu mp)^2, with mu=1.2, since Rosen works in number density, n.
!   (their cooling = Lambda*n^2,  rho=mu mp n.)
!  The factor Lambdaunits converts from cgs units to code units.
!  We've increased coolT(1), to avoid creating gas too cold to resolve.
!    
      cool=0.0
      do i=1,5
        where (coolT(i) <= TT(:) .and. TT(:) < coolT(i+1))            &
           cool(:)=cool(:) + yH*coolH(i)*rho**2*(TT(:)*unit_temperature)**coolB(i)
      enddo
!
!  add UV heating, cf. Wolfire et al., ApJ, 443, 152, 1995
!  with the values above, this gives about 0.012 erg/g/s (T < ~1.e4 K)
!  nb: need rho0 from density_[init/run]_pars, if i want to implement
!      the the arm/interarm scaling.
!
    heat=0.0
!ajwm: DISABLE UV HEATING -- Requires reformulation
!   heat(:)=rhoUV*(rho0/1.38)**1.4*unit_Lambda*coolH(3)*TUV**coolB(3)*   &
!                               0.5*(1.0+tanh(cUV*(T0UV-TT(:))))
!

    if (laverageSNheating) then
!ajwm: need to do unit_system stuff with scale heights etc.
!  Average SN heating 
!   (due to SNI)
       norm = sqrt(2*pi)*ampl_SN
       heat(:)=heat(:)+ r_SNI*rho1*norm/h_SNI*exp(-(z(l1:l2)/h_SNI)**2)
!   (due to SNII)
       heat(:)=heat(:)+ r_SNII*rho1*norm/h_SNII*exp(-(z(l1:l2)/h_SNII)**2)
    endif
!  heat and cool were calculated in terms of de/dt [~ erg/g/s], 
!  so just multiply by TT1 to get ds/dt [~ erg/g/s/K]:
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + &
                      TT1(:)*(heat(:) - cool(:))
!
    endsubroutine calc_heat_cool_interstellar
!***********************************************************************
    subroutine check_SN(f)
!
!  Checks for SNe, and implements appropriately:
!   relevant subroutines in entropy.f90
!
    use Cdata
!
    real, dimension(mx,my,mz,mvar+maux) :: f
    logical :: l_SNI=.false.   !only allow SNII if no SNI this step
                               !(may not be worth keeping)
!
    intent(inout) :: f
!
!  identifier  
!
      if(headtt) print*,'check_SN: ENTER'
!
!  Do separately for SNI (simple scheme) and SNII (Boris' scheme)
!
    if (lSNI) call check_SNI (f,l_SNI)
    if (lSNII) call check_SNII(f,l_SNI)
!
    endsubroutine check_SN
!***********************************************************************
    subroutine check_SNI(f,l_SNI)
!
!  If time for next SNI, then implement, and calculate time of subsequent SNI
!
    use Cdata
    use Mpicomm
    use General
!
    real, dimension(mx,my,mz,mvar+maux) :: f
    real, dimension(1) :: franSN
    logical :: l_SNI
!
    intent(inout) :: f,l_SNI
!
!  identifier
!
    if(headtt) print*,'check_SNI: ENTER'
!
    l_SNI=.false.
    if (t >= t_next_SNI) then
       call position_SNI(f)
       call explode_SN(f,1)
       !  pre-determine time for next SNI
       if (lroot) then
          call random_number_wrapper(franSN)   
          t_next_SNI=t_next_SNI + (1.0 + 0.4*(franSN(1)-0.5)) * t_interval_SNI
          if (ip<14) print*,'check_SNI: Next SNI at time = ',t_next_SNI
          interstellarsave(1)=t_next_SNI
       endif
       call mpibcast_real(t_next_SNI,1)
       l_SNI=.true.
    endif
!
    endsubroutine check_SNI
!***********************************************************************
    subroutine check_SNII(f,l_SNI)
!
!  Check for SNII, via self-regulating scheme.
!
    use Cdata
    use General
    use Mpicomm
! 
    real, dimension(mx,my,mz,mvar+maux) :: f
    real, dimension(nx) :: lnrho,rho,rho_cloud,ss,TT
!    real :: lnrho,rho,rho_cloud,ss,TT
    real :: mass_cloud,mass_cloud_dim,freq_SNII,prob_SNII,rate_SNII,dv
    real, dimension(1) :: franSN,fsum1,fsum1_tmp,fmpi1
    real, dimension(ncpus) :: mass_cloud_byproc
    integer :: icpu
    logical :: l_SNI
!
    intent(in) :: l_SNI
    intent(inout) :: f
!
!  identifier
!
    if(headtt) print*,'check_SNII: ENTER'
!
    if (.not. l_SNI) then         ! only do if no SNI this step
!
!  NB: currently no 'nzskip' mechanism to prevent SNII occurring near
!   top or bottom of box.  Non-trivial to implement with nprocz > 1 -- and
!   shouldn't be a real problem, as mass near boundaries should be low.
       mass_cloud=0.0
       do n=n1,n2
           do m=m1,m2
             lnrho(:)=f(l1:l2,m,n,ilnrho)
             rho(:)=exp(lnrho(:))
             ss(:)=f(l1:l2,m,n,iss)
!ajwm             TT(:)=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss(:))/gamma1*cp1
             rho_cloud(:)=0.0
             where (rho(:) >= rho_crit .and. TT(:) <= TT_crit)   &
                  rho_cloud(:)=rho(:)
             mass_cloud=mass_cloud+sum(rho_cloud(:))
          enddo
       enddo
       fsum1_tmp=(/ mass_cloud /)
       !print*,'check_SNII, iproc,fsum1_tmp:',iproc,fsum1_tmp(1)
       call mpireduce_sum(fsum1_tmp,fsum1,1) 
       call mpibcast_real(fsum1,1)
       dv=1.
       if (nxgrid/=1) dv=dv*dx
       if (nygrid/=1) dv=dv*dy
       if (nzgrid/=1) dv=dv*dz
       mass_cloud_dim=fsum1(1)*dv*tosolarMkpc3
       !print*,'check_SNII: iproc,fsum1:',iproc,fsum1(1)
       ! need convert to dimensional units, for rate/probability calculation only. 
       ! don't overwrite mass_cloud (on individual processors), as it's re-used.
       !if (lroot .and. ip < 14) &
       !     print*,'check_SNII: mass_cloud_dim:',mass_cloud_dim
       !
       freq_SNII=frac_heavy*frac_converted*mass_cloud_dim/mass_SN/tau_cloud
       prob_SNII=freq_SNII*dt
       rate_SNII=freq_SNII*1e-3
       if (Lxyz(1)/=0.) rate_SNII=rate_SNII/Lxyz(1)
       if (Lxyz(2)/=0.) rate_SNII=rate_SNII/Lxyz(2)
       if (lroot) call random_number_wrapper(franSN)   
       call mpibcast_real(franSN,1)
       ! if (lroot .and. ip < 14) &
       !      print*,'check_SNII: rate,prob,rnd:',rate_SNII,prob_SNII,franSN(1)
       if (franSN(1) <= prob_SNII) then
          !  position_SNII needs the mass_clouds for each processor;  
          !   communicate and store them here, to avoid recalculation.
          mass_cloud_byproc(:)=0.0
          ! use non-root broadcasts for the communication...
          do icpu=1,ncpus
             fmpi1=mass_cloud
             call mpibcast_real_nonroot(fmpi1,1,icpu-1)
             mass_cloud_byproc(icpu)=fmpi1(1)
          enddo
          ! if (lroot.and.ip<14) print*,'check_SNII: mass_cloud_byproc:',mass_cloud_byproc
          call position_SNII(f,mass_cloud_byproc)
          call explode_SN(f,2)
       endif
    endif
    !
  endsubroutine check_SNII
!***********************************************************************
    subroutine position_SNI(f)
!
!   determine position for next SNI (w/ fixed scale-height)
!
    use Cdata
    use Mpicomm
    use General
    use Ionization
!
    real, dimension(mx,my,mz,mvar+maux) :: f
!
    real, dimension(nzgrid) :: cum_prob_SNI
    real :: zn
    real, dimension(3) :: fran3, fmpi3
    real, dimension(4) :: fmpi4
    real, dimension(2) :: fmpi2
    integer :: i, l, nzskip=10   !prevent SNI from being too close to boundaries
    integer, dimension(4) :: impi4
!
    !intent(in) :: f
    intent(inout) :: f
!
!  identifier
!
    if(headtt) print*,'position_SNI: ENTER'
    
    !
    !  NB: this routine not usefully parallelised -- could all be done on root.
    !
    !  Lx,Ly,Lz,x0,y0,z0 not necessarily (correctly) set in cdata module. 
    !  (dx,dy,dz OK.  Lx,Ly,Lz assume periodic domains, not true for z.
    !   x0,y0,z0 not set (& note that their value also depends on periodicity.) )
    !  The following assumes x, y periodic, z non-periodic (see start.f90)
    Lx=Lxyz(1);       Ly=Lxyz(2);       Lz=Lxyz(3)
    x0=xyz0(1)+.5*dx; y0=xyz0(2)+.5*dy; z0=xyz0(3)
    
    if (lperi(1)) then; x0=xyz0(1)+.5*dx; else; x0=xyz0(1); endif
    if (lperi(2)) then; y0=xyz0(2)+.5*dy; else; y0=xyz0(2); endif
    if (lperi(3)) then; z0=xyz0(3)+.5*dz; else; z0=xyz0(3); endif
    
    !if (lroot) print*, 'position_SNI: dx,dy,dz=',dx,dy,dz
    !if (lroot) print*, 'position_SNI: x0,y0,z0,Lx,Ly,Lz=',x0,y0,z0,Lx,Ly,Lz
    
    !
    !  Pick SN position (x_SN,y_SN,z_SN)
    !
    call random_number(fran3)    ! get 3 random numbers
                                 ! on all processors to keep
                                 ! rnd. generators in sync
    if (lroot) then
       if (ltestSN) then
          i=int(nxgrid/2)+1
          x_SN=x0+(i-1)*dx
          l_SN=(i-1)
       else
          i=int(fran3(1)*nxgrid)+1
          x_SN=x0+(i-1)*dx
          l_SN=(i-1)
       endif
       !if (lroot) print*, 'position_SNI: fran3(1),x = ',fran3(1),x_SN

       if (ltestSN) then
          i=int(nygrid/2)+1
          y_SN=y0+(i-1)*dy
          ipy_SN=(i-1)/ny  ! uses integer division
          m_SN=(i-1)-(ipy_SN*ny)+nghost
          !if (lroot) print*, 'position_SNI: fran3(2),y_SN,ipy_SN=',fran3(2),y_SN,ipy_SN
       else
          i=int(fran3(2)*nygrid)+1
          y_SN=y0+(i-1)*dy
          ipy_SN=(i-1)/ny  ! uses integer division
          m_SN=(i-1)-(ipy_SN*ny)+nghost
          !if (lroot) print*, 'position_SNI: fran3(2),i,y_SN,ipy_SN=',fran3(2),i,y_SN,ipy_SN
       endif

       if (ltestSN) then
          i=int(nzgrid/2)+1
          z_SN=z0+(i-1)*dz
          ipz_SN=(i-1)/nz   ! uses integer division
          n_SN=(i-1)-(ipz_SN*nz)+nghost
       elseif (uniform_zdist_SNI) then
          i=int(fran3(3)*nzgrid)+1
          z_SN=z0+(i-1)*dz
          ipz_SN=(i-1)/nz   ! uses integer division
          n_SN=(i-1)-(ipz_SN*nz)+nghost
       else
          !
          !  Cumulative probability funcion in z currently calculated each time.
          !  It's constant, and could be stored (and calculated in init)
          cum_prob_SNI(1:nzskip)=0.0
          do n=nzskip+1,nzgrid-nzskip
             zn=z0+(n-1)*dz
             cum_prob_SNI(n)=cum_prob_SNI(n-1)+exp(-(zn/h_SNI)**2)
          enddo
          cum_prob_SNI=cum_prob_SNI/cum_prob_SNI(nzgrid-nzskip)
          !  The following should never be needed, but just in case floating point 
          !  errors ever lead to cum_prob_SNI(nzgrid-nzskip) < rnd < 1.
          cum_prob_SNI(nzgrid-nzskip+1:nzgrid)=1.0   
          !if (lroot) print*, 'position_SNI: cum_prob_SNI=',cum_prob_SNI
          
          do i=nzskip+1,nzgrid-nzskip
             if (cum_prob_SNI(i-1) <= fran3(3) .and. fran3(3) < cum_prob_SNI(i)) &
                  then
                z_SN=z0+(i-1)*dz
                ipz_SN=(i-1)/nz  ! uses integer division
                n_SN=(i-1)-(ipz_SN*nz)+nghost
                exit
             endif
          enddo
       endif
       iproc_SN=ipz_SN*nprocy + ipy_SN
       !if (lroot) print*, 'position_SNI: fran3(3),z_SN,ipz_SN,iproc_SN=', &
       !                fran3(3),z_SN,ipz_SN,iproc_SN
    endif

!
!  Broadcast position to all processors from root;
!  also broadcast iproc_SN, needed for later broadcast of rho_SN.
!
    fmpi3=(/ x_SN, y_SN, z_SN /)
    call mpibcast_real(fmpi3,3)
    x_SN=fmpi3(1); y_SN=fmpi3(2); z_SN=fmpi3(3)
!
    impi4=(/ iproc_SN, l_SN, m_SN, n_SN /)
    call mpibcast_int(impi4,4)
    iproc_SN=impi4(1)
    l_SN=impi4(2)
    m_SN=impi4(3)
    n_SN=impi4(4)
!
!  With current SN scheme, we need rho at the SN location.
!
    rho_SN=0.0
    if (iproc==iproc_SN) then
      !print*, 'position_SNI:',l_SN,m_SN,n_SN,ilnrho,f(l_SN,m_SN,n_SN,ilnrho)
      lnrho_SN=f(l_SN,m_SN,n_SN,ilnrho)
!      rho_SN=exp(lnrho_SN)
! calculate TT_SN here, for later use in explode_SN
      ss_SN=f(l_SN,m_SN,n_SN,iss)

! NEED TO USE IONISATION CALCS 

!      call ioncalc(lnrho,ss,cs2,TT1,cp1tilde, &
!        Temperature=TT_SN_new,InternalEnergy=ee,IonizationFrac=yH)
!      TT_SN=cs20*exp(gamma1*(lnrho_SN-lnrho0)+gamma*ss_SN)/gamma1*cp1
      call ionget(lnrho_SN,ss_SN,yH_SN,TT_SN)
      call thermodynamics(lnrho_SN,ss_SN,yH_SN,TT_SN,ee=ee_SN)

      if (lroot.and.ip<14) &
           print*, 'position_SNI:',l_SN,m_SN,n_SN,x_SN,y_SN,z_SN,rho_SN,TT_SN
    endif
!
!  Broadcast lnrho_SN, TT_SN, etc. to all processors.
!
    fmpi4=(/ lnrho_SN, ee_SN, ss_SN, TT_SN /)
    call mpibcast_real_nonroot(fmpi4,4,iproc_SN)
    lnrho_SN=fmpi4(1); ee_SN=fmpi4(2); ss_SN=fmpi4(3); TT_SN=fmpi4(4)
    rho_SN=exp(lnrho_SN);

    if (lroot.and.ip<=14) &
         print*, 'position_SNI:',iproc_SN,x_SN,y_SN,z_SN,rho_SN,TT_SN,ee_SN,ss_SN
!
    endsubroutine position_SNI
!***********************************************************************
    subroutine position_SNII(f,mass_cloud_byproc)
!
!  Determine position for next SNII (using Boris' scheme)
!  It seems impractical to sort all high density points across all processors;
!  instead, we just construct cumulative pdfs that allow us to pick a processor,
!  and then a point on that processor, with probability proportional to rho.
!  As a result, the SN position is *not* independent of ncpus (or of nprocy 
!  and nprocz).  (It is repeatable given fixed nprocy/z though.)
!
    use Cdata
    use General
    use Mpicomm
!
    real, dimension(mx,my,mz,mvar+maux) :: f
    real, dimension(ncpus) :: mass_cloud_byproc
    real, dimension(0:ncpus) :: cum_prob_byproc
    real, dimension(1) :: franSN
    real, dimension(5) :: fmpi5
    real :: mass_cloud,cum_mass,cum_prob_onproc
    real :: lnrho,rho,ss,TT
    integer :: icpu, l
!
    !intent(in) :: f,mass_cloud_byproc
    intent(in) :: mass_cloud_byproc
    intent(inout) :: f
!
!  identifier
!
      if(headtt) print*,'position_SNII: ENTER'
!
!  Construct cumulative distribution function, using mass_cloud_byproc.
!  NB: icpu=iproc+1 (iproc in [0,ncpus-1], icpu in [1,ncpus] )
!
    cum_prob_byproc=0.0
    do icpu=1,ncpus
      mass_cloud=mass_cloud_byproc(icpu)
      cum_prob_byproc(icpu)=cum_prob_byproc(icpu-1)+mass_cloud_byproc(icpu)
    enddo
    cum_prob_byproc(:)=cum_prob_byproc(:)/cum_prob_byproc(ncpus)
    if (lroot) then
      print*,'position_SNII: mass_cloud_byproc=',mass_cloud_byproc
      print*,'position_SNII: cum_prob_byproc=',cum_prob_byproc
      print*,'position_SNII: mass_cloud=',mass_cloud
    endif
!     
!  Use random number to detemine which processor SN is on.
!  (Use root processor for rand, to ensure repeatability.)
!
    if (lroot) call random_number_wrapper(franSN)   
    call mpibcast_real(franSN,1)
    do icpu=1,ncpus
      if (cum_prob_byproc(icpu-1) <= franSN(1) .and.                      &
           franSN(1) < cum_prob_byproc(icpu)) then
        iproc_SN=icpu-1 
        exit
      endif
    enddo
    if (lroot) print*, 'position_SNII: franSN(1),iproc_SN=',franSN(1),iproc_SN
!
!  Use random number to pick SNII location on the right processor.
!  (No obvious reason to re-use the original random number for this.)
!    franSN(1)=(franSN(1)-cum_prob_byproc(iproc_SN)) /                      &
!              (cum_prob_byproc(iproc_SN+1)-cum_prob_byproc(iproc_SN))
!
    if (lroot) call random_number_wrapper(franSN)   
    call mpibcast_real(franSN,1)
    if (iproc == iproc_SN) then
      cum_mass=0.0
find_SN: do n=n1,n2
           do m=m1,m2
             do l=l1,l2
               lnrho=f(l,m,n,ilnrho)
               rho=exp(lnrho)
               ss=f(l,m,n,iss)
!ajwm: should use thermodynamics subroutine but need to resolve temperature unit issue first
!               TT=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss)/gamma1*cp1
               if (rho >= rho_crit .and. TT <= TT_crit) then
                 cum_mass=cum_mass+rho
                 cum_prob_onproc=cum_mass/mass_cloud
                 if (franSN(1) <= cum_prob_onproc) then
                    x_SN=x(l); y_SN=y(m); z_SN=z(n); rho_SN=rho; TT_SN=TT
                    print*,'position_SNII: cum_mass,cum_prob_onproc,franSN(1)=', &
                                          cum_mass,cum_prob_onproc,franSN(1)
                    exit find_SN
                 endif
               endif
             enddo
           enddo
         enddo find_SN
         print*,'position_SNII: iproc_SN,x_SN,y_SN,z_SN,rho_SN=',iproc,x_SN,y_SN,z_SN,rho_SN
     endif
    fmpi5=(/ x_SN, y_SN, z_SN, rho_SN, TT_SN /)
    call mpibcast_real_nonroot(fmpi5,5,iproc_SN)
    x_SN=fmpi5(1); y_SN=fmpi5(2); z_SN=fmpi5(3)
    rho_SN=fmpi5(4); TT_SN=fmpi5(5)
!
    if (lroot) &
         print*, 'position_SNII: iproc,x_SN,y_SN,z_SN,rho_SN,TT_SN=', &
                                                             iproc,x_SN,y_SN,z_SN,rho_SN,TT_SN
!
    endsubroutine position_SNII
!***********************************************************************
    subroutine explode_SN(f,itype_SN)
      !
      !  Implement SN (of either type), at pre-calculated position
      !  (This can all be made more efficient, after debugging.)
      !
      !  ??-nov-02/grs : coded from GalaxyCode                        
      !  20-may-03/tony: pencil formulation and broken into subroutines
      !
      use Cdata
      use Mpicomm
      use Ionization
      !
      real, intent(inout), dimension(mx,my,mz,mvar+maux) :: f
      integer, intent(in) :: itype_SN

      real :: width_SN,width_shell_outer,width_shell_inner,c_SN
      real :: profile_integral, mass_shell, mass_gain
      real :: EE_SN=0.,EE2_SN=0.,rho_SN_new,lnrho_SN_new,ss_SN_new,yH_SN_new,TT_SN_new,dv
      
      integer, parameter :: point_width=4
     ! integer :: point_width=8            ! vary for debug, tests
      real, dimension(nx) :: deltarho, deltaEE, ee_old
      real, dimension(nx) :: rho_old, TT_old, TT_new, ss_new,lnrho_new
      real, dimension(1) :: fmpi1, fmpi1_tmp
      real, dimension(2) :: fmpi2, fmpi2_tmp
     ! the following not really wanted, but are required if we want
     ! to use the thermodynamics function in the [No]Ionisation module
      real, dimension(nx) ::  lnrho_old, ss_old, TT1_old, cs2_old, yH_old

      logical :: lmove_mass=.false.
      integer :: idim
          
      !
      !  identifier
      !
      if(lroot.and.ip<14) print*,'explode_SN: itype_SN=',itype_SN
      !
      width_SN=point_width*dxmin      
      idim=0                         !allow for 1-d, 2-d and 3-d profiles...
      if (nxgrid /=1) idim=idim+1
      if (nygrid /=1) idim=idim+1
      if (nzgrid /=1) idim=idim+1
      c_SN=ampl_SN/cnorm_SN(idim)/width_SN**idim
      dv=1.
      if (nxgrid/=1) dv=dv*dx
      if (nygrid/=1) dv=dv*dy
      if (nzgrid/=1) dv=dv*dz

      if (lroot.and.ip<=14) print*,'explode_SN: width_SN,c_SN=', width_SN,c_SN
        
      !
      !  Now deal with (if nec.) mass relocation
      !

      call perturb_energy(lnrho_SN,(ee_SN*rho_SN)+c_SN,ss_SN_new,TT_SN_new)

      if(lroot) print*, &
         'explode_SN: TT_SN, TT_SN_new, TT_SN_min =', &
                                TT_SN,TT_SN_new,TT_SN_min

      if (TT_SN_new < TT_SN_min) then
         lmove_mass=.not.lnever_move_mass
         ! lmove_mass=.false.  ! use to switch off for debug...

         ! The bit that BREAKS the pencil formulation...
         ! must know the total moved mass BEFORE attempting mass relocation 

         rho_SN_new=c_SN/TT_SN_min*TTunits

         ! ASSUME: SN will fully ionize the gas at its centre
         call getdensity((ee_SN*rho_SN)+c_SN,TT_SN_min,1.,rho_SN_new)
         lnrho_SN_new=alog(rho_SN_new)

         call perturb_energy(lnrho_SN_new,(ee_SN*rho_SN)+c_SN,ss_SN_new,TT_SN_new)

         if(lroot) print*, &
            'explode_SN: Relocate mass... TT_SN_new, rho_SN_new=', &
                                                     TT_SN_new,rho_SN_new

         call calcmassprofileintegral_SN(f,width_SN,profile_integral)
         fmpi1_tmp=(/ profile_integral /)
         call mpireduce_sum(fmpi1_tmp,fmpi1,1) 
         call mpibcast_real(fmpi1,1)
         profile_integral=fmpi1(1)*dv
         mass_shell=-(rho_SN_new-rho_SN)*profile_integral
         if (lroot.and.ip<14) &
           print*, 'explode_SN: mass_shell=',mass_shell
         mass_gain=0.
      endif
      

      EE_SN=0. !; EE_SN2=0.
      do n=n1,n2
         do m=m1,m2
            call proximity_SN()
            deltarho(:)=0.
            if (lmove_mass) then
               call makecavity_SN(deltarho,width_SN,rho_SN_new-rho_SN, &
                      mass_shell,cnorm_SN(idim),idim,mass_gain)
            endif

            call injectenergy_SN(deltaEE,width_SN,c_SN,EE_SN)
            ! Apply perturbations
            lnrho_old=f(l1:l2,m,n,ilnrho)
            rho_old=exp(lnrho_old)
            ss_old=f(l1:l2,m,n,iss)

            call ionget(f,yH_old,TT_old)
            call thermodynamics(lnrho_old,ss_old,yH_old,TT_old,ee=ee_old)

            ! use amax1 with rho_min to ensure rho doesn't go negative
            if (lmove_mass) then
              lnrho_new(:)=alog(amax1(rho_old(:)+deltarho(:),rho_min))
            else
              lnrho_new(:)=lnrho_old(:)
            endif
  
            call perturb_energy(lnrho_new, &
                                   (ee_old*rho_old)+deltaEE,ss_new,TT_new)

            if (lmove_mass) f(l1:l2,m,n,ilnrho)=lnrho_new
            f(l1:l2,m,n,iss)=ss_new

       enddo
      enddo
      

      ! Sum and share diagnostics etc. amongst processors
      fmpi2_tmp=(/ mass_gain, EE_SN /)
      call mpireduce_sum(fmpi2_tmp,fmpi2,4) 
      call mpibcast_real(fmpi2,4)
      mass_gain=fmpi2(1)*dv
      EE_SN=fmpi2(2)*dv
! Extra debug - no longer calculated 
!      EE2_SN=fmpi3(3)*dv; 

      if (lroot.and.ip<14) print*, &
           'explode_SN: mass_gain=',mass_gain
     
      if (lroot) then
         open(1,file=trim(datadir)//'/time_series.dat',position='append')
         write(1,'(a,1e11.3," ",i1," ",i2," ",1e11.3," ",1e11.3," ",1e11.3," ",1e11.3," ",1e11.3,a)')  &
                     '#ExplodeSN: (t,type,iproc,x,y,z,rho,energy)=(', &
                     t,itype_SN,iproc_SN,x_SN,y_SN,z_SN,rho_SN,EE_SN, ')'
         write(6,'(a,1e11.3," ",i1," ",i2," ",1e11.3," ",1e11.3," ",1e11.3," ",1e11.3," ",1e11.3,a)')  &
                     '#ExplodeSN: (t,type,iproc,x,y,z,rho,energy)=(', &
                     t,itype_SN,iproc_SN,x_SN,y_SN,z_SN,rho_SN,EE_SN, ')'
         close(1)
      endif
      
    endsubroutine explode_SN

!***********************************************************************
    subroutine calcmassprofileintegral_SN(f,width_SN,profile_integral)
!
!  Calculate integral of mass cavity profile  
!
!  22-may-03/ajwm: coded
!
      use CData

      real, intent(in), dimension(mx,my,mz,mvar+maux) :: f
      real, intent(in) :: width_SN
      real, intent(out) :: profile_integral
      real :: dx_SN_in,dx_SN_out_x0,dx_SN_out_x1,dy_SN_in,dy_SN_out_y
      real :: dy_SN_out_x0a,dy_SN_out_x0b,dy_SN_out_x1a,dy_SN_out_x1b
      real :: dz_SN,yshift
      integer :: l,mshift,il,im,in
     
      !
      !  Obtain distance to SN
      !

      profile_integral=0.
      do n=n1,n2
         in=n-nghost
         dz_SN=abs(z(n)-z_SN)
         do m=m1,m2
            im=m-nghost
            !  consider all possible positions in xy plane, to get shortest
            !  can be made more efficient later
            !  y-separations with no boundaries crossed in x
            dy_SN_in=abs(y(m)-y_SN)                      !dyi
            dy_SN_out_y=Ly-dy_SN_in                      !dyii
            !  y-separations across sliding periodic boundary at x=x0
            yshift=y(m)-deltay
            mshift=0
            if (yshift < y0) then 
               mshift=int((y0-yshift)/Ly)+1
               yshift=yshift+mshift*Ly
            elseif (yshift > y0+Ly) then      
               mshift=int((yshift-(y0+Ly))/Ly)+1
               yshift=yshift-mshift*Ly
            endif
            dy_SN_out_x0a=abs(yshift-y_SN)               !dyiii
            dy_SN_out_x0b=Ly-dy_SN_out_x0a               !dyiv
            !  y-separations across sliding periodic boundary at x=x1
            yshift=y(m)+deltay
            mshift=0
            if (yshift < y0) then 
               mshift=int((y0-yshift)/Ly)+1
               yshift=yshift+mshift*Ly
            elseif (yshift > y0+Ly) then
               mshift=int((yshift-(y0+Ly))/Ly)+1
               yshift=yshift-mshift*Ly
            endif
            dy_SN_out_x1a=abs(yshift-y_SN)               !dyv
            dy_SN_out_x1b=Ly-dy_SN_out_x1a               !dyvi
            do l=l1,l2
               il=l-nghost
               !  x-separations associated with each of the above
               dx_SN_in=abs(x(l)-x_SN)                    !dxi=dxii
               dx_SN_out_x0=Lx+(x(l)-x_SN)                !dxiii=dxiv
               dx_SN_out_x1=Lx-(x(l)-x_SN)                !dxv=dxvi
               dr2_SN(il)=min( dx_SN_in**2 + dy_SN_in**2,             &
                    dx_SN_in**2 + dy_SN_out_y**2,          &
                    dx_SN_out_x0**2 + dy_SN_out_x0a**2,    &
                    dx_SN_out_x0**2 + dy_SN_out_x0b**2,    &
                    dx_SN_out_x1**2 + dy_SN_out_x1a**2,    &
                    dx_SN_out_x1**2 + dy_SN_out_x1b**2 )   &
                    + dz_SN**2
            enddo
            profile_integral = profile_integral + sum(exp(-(dr2_SN(:)/width_SN**2)**3))
         enddo
      enddo

    endsubroutine calcmassprofileintegral_SN
!***********************************************************************
    subroutine proximity_SN()
!
!  Calculate pencil of distance to SN explosion site
!
!  20-may-03/ajwm: extracted from explode_SN code written by grs
!  22-may-03/ajwm: pencil formulation
!
!
      use CData

      real :: dx_SN_in,dx_SN_out_x0,dx_SN_out_x1,dy_SN_in,dy_SN_out_y
      real :: dy_SN_out_x0a,dy_SN_out_x0b,dy_SN_out_x1a,dy_SN_out_x1b
      real :: dz_SN,yshift
      integer :: l,mshift,il,im,in
     
      !
      !  Obtain distance to SN
      !
         im=m-nghost
         in=n-nghost
         dz_SN=abs(z(n)-z_SN)
         !  consider all possible positions in xy plane, to get shortest
         !  can be made more efficient later
         !  y-separations with no boundaries crossed in x
         dy_SN_in=abs(y(m)-y_SN)                      !dyi
         dy_SN_out_y=Ly-dy_SN_in                      !dyii
         !  y-separations across sliding periodic boundary at x=x0
         yshift=y(m)-deltay
         mshift=0
         if (yshift < y0) then 
            mshift=int((y0-yshift)/Ly)+1
            yshift=yshift+mshift*Ly
         elseif (yshift > y0+Ly) then      
            mshift=int((yshift-(y0+Ly))/Ly)+1
            yshift=yshift-mshift*Ly
         endif
         dy_SN_out_x0a=abs(yshift-y_SN)               !dyiii
         dy_SN_out_x0b=Ly-dy_SN_out_x0a               !dyiv
         !  y-separations across sliding periodic boundary at x=x1
         yshift=y(m)+deltay
         mshift=0
         if (yshift < y0) then 
            mshift=int((y0-yshift)/Ly)+1
            yshift=yshift+mshift*Ly
         elseif (yshift > y0+Ly) then
            mshift=int((yshift-(y0+Ly))/Ly)+1
            yshift=yshift-mshift*Ly
         endif
         dy_SN_out_x1a=abs(yshift-y_SN)               !dyv
         dy_SN_out_x1b=Ly-dy_SN_out_x1a               !dyvi
         do l=l1,l2
            il=l-nghost
            !  x-separations associated with each of the above
            dx_SN_in=abs(x(l)-x_SN)                    !dxi=dxii
            dx_SN_out_x0=Lx+(x(l)-x_SN)                !dxiii=dxiv
            dx_SN_out_x1=Lx-(x(l)-x_SN)                !dxv=dxvi
            dr2_SN(il)=min( dx_SN_in**2 + dy_SN_in**2,             &
                 dx_SN_in**2 + dy_SN_out_y**2,          &
                 dx_SN_out_x0**2 + dy_SN_out_x0a**2,    &
                 dx_SN_out_x0**2 + dy_SN_out_x0b**2,    &
                 dx_SN_out_x1**2 + dy_SN_out_x1a**2,    &
                 dx_SN_out_x1**2 + dy_SN_out_x1b**2 )   &
                 + dz_SN**2
         enddo
       
    endsubroutine proximity_SN
!***********************************************************************
    subroutine makecavity_SN(deltarho,width_SN,depth,mass_shell, &
                             cnorm_dim,idim,mass_gain)   
      use CData
      !
      real, intent(in) :: width_SN, depth, mass_shell, cnorm_dim
      real, intent(inout) :: mass_gain
      real, intent(out), dimension(nx) :: deltarho
      integer, intent(in) :: idim
      !
      real, dimension(nx) :: profile_shell_outer,profile_shell_inner
      real, dimension(1) :: fsum1,fsum1_tmp
      real :: width_shell_outer, width_shell_inner, c_shell
      ! real, dimension(1) :: fsum1,fsum1_tmp

      width_shell_outer=2.0*width_SN
      width_shell_inner=1.5*width_SN

      deltarho(:) =  depth*exp(-(dr2_SN(:)/width_SN**2)**3)
      
      c_shell=mass_shell /                                  &
           (cnorm_dim*((width_shell_outer**idim)-(width_shell_inner**idim)))

!      if (lroot) print*, &
!           'explode_SN, c_shell:',c_shell
      !  add missing mass back into shell
      
      profile_shell_outer(:)=                              &
           exp(-(dr2_SN(:)/width_shell_outer**2)**3)
      profile_shell_inner(:)=                              &
           exp(-(dr2_SN(:)/width_shell_inner**2)**3)
      
      deltarho(:)=deltarho(:) + c_shell *      &
           (profile_shell_outer(:) - profile_shell_inner(:))
      mass_gain=mass_gain + sum(deltarho(:))  
      
    endsubroutine makecavity_SN

!***********************************************************************
    subroutine injectenergy_SN(deltaEE,width_SN,c_SN,EE_SN)
      use CData
      !
      real, intent(in) :: width_SN,c_SN
      real, intent(inout) :: EE_SN
      real, intent(out), dimension(nx) :: deltaEE
      !
      real, dimension(nx) :: profile_SN
      
      ! Whether mass moved or not, inject energy.
      !

      profile_SN=exp(-(dr2_SN(:)/width_SN**2)**3)

      deltaEE(:)=c_SN*profile_SN(:) ! spatial energy density 
      EE_SN=EE_SN+sum(deltaEE(:))   

    endsubroutine injectenergy_SN

endmodule interstellar
