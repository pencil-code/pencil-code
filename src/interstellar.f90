! $Id: interstellar.f90,v 1.21 2003-05-29 11:19:30 ngrs Exp $

!  This modules contains the routines for SNe-driven ISM simulations.
!  Still in development. 

module Interstellar

  use Cparam
  use Cdata
  use Density

  implicit none

  real :: x_SN,y_SN,z_SN,rho_SN,lnrho_SN,TT_SN,ss_SN,ampl_SN=1.0
  integer :: l_SN,m_SN,n_SN
  real, dimension(nx) :: dr2_SN     ! Pencil storing radius to SN
  real :: t_next_SNI=0.0,t_interval_SNI=3.64e-3,h_SNI=0.325
  real :: tau_cloud=2e-2
  integer, parameter :: ninterstellarsave=1
  real, dimension(ninterstellarsave) :: interstellarsave
  real, parameter :: rho_crit=1.,TT_crit=4000.
  real, parameter :: frac_converted=0.02,frac_heavy=0.10,mass_SN=10.
  real, parameter :: cnorm_SN = 3.71213666 !  3.128289613 
                         ! (int exp(-(r/a)^6) 4\pi r^2 dr) / (a^3)
                         ! for gaussian of width a

!  cp1=1/cp used to convert TT (and ss) into interstellar code units
!  (useful, as many conditions conveniently expressed in terms of TT)
!  code units based on:
!    [length]  = 1kpc  = 3.09 10^21 cm
!    [time]    = 1Gyr  = 3.15 10^16 s             !no, on [u]=1km/s...
!    [rho]     =       = 1.00 10^-24 g/cm^3
!  Lambdaunits converts coolH into interstellar code units.
!   (this should really just be incorporated into coolH coefficients)

  real, parameter :: cp1=27.8   !=R * gamma / (mu * (gamma-1))  27.8 
  real, parameter :: TTunits=46.6,tosolarMkpc3=1.483e7
  real, parameter :: Lambdaunits=3.29e-18
  real, parameter :: rhoUV=0.1,TUV=7000.,T0UV=12000.,cUV=5.e-4
  real, parameter, dimension(6) ::                                          &
        coolT=(/ 500.,     2000.,    8000.,    1.e5,    4.e7,     1.e9 /),  &
        coolH=(/ 5.595e15, 2.503e17, 1.156e12, 4.45e29, 8.054e20, 0.   /),  &
        coolB=(/ 2.,       1.5,      2.867,    -.65,    0.5,      0.   /)
  integer :: iproc_SN,ipy_SN,ipz_SN

  ! input parameters
 
  integer :: dummy 
  namelist /interstellar_init_pars/ dummy

  ! run parameters
  namelist /interstellar_run_pars/ &
      t_next_SNI,t_interval_SNI,h_SNI,ampl_SN,tau_cloud

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
        print*, 'Register_interstellar'
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: interstellar.f90,v 1.21 2003-05-29 11:19:30 ngrs Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_interstellar: nvar > mvar')
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
!
      logical, save :: first=.true.
      logical :: lstart
      logical :: exist

      if (first) then
         if (lroot.and.ip<14) then
            print*, 'initialize_interstellar: reading seed file'
            print*, 'initialize_interstellar: nseed,seed',nseed,seed(1:nseed)
         endif
         if (.not. lstart) then
            call inpui(trim(directory)//'/seed.dat',seed,nseed)
            call random_seed_wrapper(put=seed(1:nseed))
         endif
!
         if (lroot) then
            inquire(file=trim(datadir)//'/interstellar.dat',exist=exist)
            if (exist) then 
               if (ip<14) print*, 'initialize_interstellar: read interstellar.dat'
                call inpup(trim(datadir)//'/interstellar.dat',  &
                    interstellarsave,ninterstellarsave)
               if (ip<14) print*, 'initialize_interstellar: t_next_SNI',t_next_SNI
            else
               interstellarsave(1)=t_next_SNI
            endif
         endif
         call mpibcast_real(interstellarsave,1)
         t_next_SNI=interstellarsave(1)
      endif

      if (lroot.and.ip<14) then
         print*,'initialize_interstellar: nseed,seed',nseed,seed(1:nseed)
         print*,'initialize_interstellar: finished'
      endif
!
    endsubroutine initialize_interstellar
!***********************************************************************
    subroutine calc_heat_cool_interstellar(df,rho1,TT1)
!
!  adapted from calc_heat_cool
!
      use Cdata
      use Mpicomm
      use Density, only : rho0
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: rho1,TT1
      real, dimension (nx) :: TT,heat,cool 
      integer :: i
!
      intent(in) :: rho1,TT1
      intent(inout) :: df
!
!  identifier
!
      if(headtt) print*,'calc_heat_cool_interstellar'
!
!  define T in K, for calculation of both UV heating and radiative cooling
!
      TT=cp1/TT1      ! leave for now (used several times)
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
        where (coolT(i) <= TT(:) .and. TT(:) < coolT(i+1))                 &
           cool(:)=cool(:) + Lambdaunits/rho1(:)*coolH(i)*TT(:)**coolB(i)
      enddo
!! the following may be easier to optimise? (or use do independent?)
!!      where (coolT(1) <= TT(:) .and. TT(:) < coolT(2)) 
!!        cool(:)=cool(:) + Lambdaunits/rho1(:)*coolH(1)*TT(:)**coolB(1)
!!      elsewhere (coolT(2) <= TT(:) .and. TT(:) < coolT(3)) 
!!        cool(:)=cool(:) + Lambdaunits/rho1(:)*coolH(2)*TT(:)**coolB(2)
!!      elsewhere (coolT(3) <= TT(:) .and. TT(:) < coolT(4)) 
!!        cool(:)=cool(:) + Lambdaunits/rho1(:)*coolH(3)*TT(:)**coolB(3)
!!      elsewhere (coolT(4) <= TT(:) .and. TT(:) < coolT(5)) 
!!        cool(:)=cool(:) + Lambdaunits/rho1(:)*coolH(4)*TT(:)**coolB(4)
!!      elsewhere (coolT(5) <= TT(:)) 
!!        cool(:)=cool(:) + Lambdaunits/rho1(:)*coolH(5)*TT(:)**coolB(5)
!!      endwhere
!
!  add UV heating, cf. Wolfire et al., ApJ, 443, 152, 1995
!  with the values above, this gives about 0.012 erg/g/s (T < ~1.e4 K)
!  nb: need rho0 from density_[init/run]_pars, if i want to implement
!      the the arm/interarm scaling.
!
      heat(:)=rhoUV*(rho0/1.38)**1.4*Lambdaunits*coolH(3)*TUV**coolB(3)*   &
                               0.5*(1.0+tanh(cUV*(T0UV-TT(:))))
!
!  heat and cool were calculated in terms of de/dt [~ erg/g/s], 
!  so just multiply by TT1 to get ds/dt [~ erg/g/s/K]:
!
      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + TT1(:)*(heat(:) - cool(:))
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
    real, dimension(mx,my,mz,mvar) :: f
    logical :: l_SNI=.false.   !only allow SNII if no SNI this step
                               !(may not be worth keeping)
!
    intent(inout) :: f
!
!  identifier  
!
      if(headtt) print*,'check_SN'
!
!  Do separately for SNI (simple scheme) and SNII (Boris' scheme)
!
    call check_SNI (f,l_SNI)
    ! call check_SNII(f,l_SNI)  !ngrs: leave out for faster debug
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
    real, dimension(mx,my,mz,mvar) :: f
    real, dimension(1) :: franSN
    logical :: l_SNI
!
    intent(inout) :: f,l_SNI
!
!  identifier
!
    if(headtt) print*,'check_SNI'
!
    l_SNI=.false.
    if (t >= t_next_SNI) then
      call position_SNI(f)
      call explode_SN(f,1)
!  pre-determine time for next SNI
      if (lroot) then
        call random_number_wrapper(franSN)   
        t_next_SNI=t_next_SNI + (1.0 + 0.4*(franSN(1)-0.5)) * t_interval_SNI
        if (ip<14) print*,'Next SNI at time: ',t_next_SNI
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
    real, dimension(mx,my,mz,mvar) :: f
    real, dimension(nx) :: lnrho,rho,rho_cloud,ss,TT
!    real :: lnrho,rho,rho_cloud,ss,TT
    real :: mass_cloud,mass_cloud_dim,freq_SNII,prob_SNII,rate_SNII
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
    if(headtt) print*,'check_SNII'
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
             ss(:)=f(l1:l2,m,n,ient)
             TT(:)=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss(:))/gamma1*cp1
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
       mass_cloud_dim=fsum1(1)*(dx*dy*dz)*tosolarMkpc3
       !print*,'check_SNII, iproc,fsum1:',iproc,fsum1(1)
       ! need convert to dimensional units, for rate/probability calculation only. 
       ! don't overwrite mass_cloud (on individual processors), as it's re-used.
       !if (lroot .and. ip < 14) &
       !     print*,'check_SNII, mass_cloud_dim:',mass_cloud_dim
       !
       freq_SNII=frac_heavy*frac_converted*mass_cloud_dim/mass_SN/tau_cloud
       prob_SNII=freq_SNII*dt
       rate_SNII=freq_SNII*1e-3
       if (Lxyz(1)/=0.) rate_SNII=rate_SNII/Lxyz(1)
       if (Lxyz(2)/=0.) rate_SNII=rate_SNII/Lxyz(2)
       if (lroot) call random_number_wrapper(franSN)   
       call mpibcast_real(franSN,1)
       ! if (lroot .and. ip < 14) &
       !      print*,'check_SNII, rate,prob,rnd:',rate_SNII,prob_SNII,franSN(1)
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
          ! if (lroot.and.ip<14) print*,'check_SNII, mass_cloud_byproc:',mass_cloud_byproc
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
!
    real, dimension(mx,my,mz,mvar) :: f
!
    real, dimension(nzgrid) :: cum_prob_SNI
    real :: zn
    real, dimension(3) :: fran3,fmpi3
    real, dimension(2) :: fmpi2
    integer :: i, l, nzskip=10   !prevent SNI from being too close to boundaries
    integer, dimension(1) :: impi
    logical :: lfound,lfoundx,lfoundy,lfoundz
!
    !intent(in) :: f
    intent(inout) :: f
!
!  identifier
!
      if(headtt) print*,'position_SNI'
!
!  NB: this routine not usefully parallelised -- could all be done on root.
!
!  Lx,Ly,Lz,x0,y0,z0 not necessarily (correctly) set in cdata module. 
!  (dx,dy,dz OK.  Lx,Ly,Lz assume periodic domains, not true for z.
!   x0,y0,z0 not set (& note that their value also depends on periodicity.) )
!  The following assumes x, y periodic, z non-periodic (see start.f90)
    Lx=Lxyz(1);       Ly=Lxyz(2);       Lz=Lxyz(3)
    !dx=Lx/nxgrid;     dy=Ly/nygrid;     dz=Lz/(nzgrid-1)    !already OK
    x0=xyz0(1)+.5*dx; y0=xyz0(2)+.5*dy; z0=xyz0(3)
    !if (lroot) print*, 'dx,dy,dz',dx,dy,dz
    !if (lroot) print*, 'x0,y0,z0,Lx,Ly,Lz',x0,y0,z0,Lx,Ly,Lz
!
!  Cumulative probability funcion in z currently calculated each time.
!  It's constant, and could be stored (and calculated in init)
!
    cum_prob_SNI(1:nzskip)=0.0
    do n=nzskip+1,nzgrid-nzskip
      zn=z0+(n-1)*dz
      cum_prob_SNI(n)=cum_prob_SNI(n-1)+exp(-(zn/h_SNI)**2)
    enddo
    cum_prob_SNI=cum_prob_SNI/cum_prob_SNI(nzgrid-nzskip)
!  The following should never be needed, but just in case floating point 
!  errors ever lead to cum_prob_SNI(nzgrid-nzskip) < rnd < 1.
    cum_prob_SNI(nzgrid-nzskip+1:nzgrid)=1.0   
    !if (lroot) print*, 'cum_prob_SNI',cum_prob_SNI
!
!  Pick SN position (x_SN,y_SN,z_SN)
!
    call random_number(fran3)    ! get 3 random numbers
    i=int(fran3(1)*nxgrid)+1
    x_SN=x0+(i-1)*dx
    !if (lroot) print*, 'x',fran3(1),x_SN
!
    i=int(fran3(2)*nygrid)+1
    y_SN=y0+(i-1)*dy
    ipy_SN=(i-1)/ny
    !if (lroot) print*, 'y',fran3(2),i,y_SN,ipy_SN
!
    do n=nzskip+1,nzgrid-nzskip
      if (cum_prob_SNI(n-1) <= fran3(3) .and. fran3(3) < cum_prob_SNI(n)) &
        z_SN=z0+(n-1)*dz
    enddo
    ipz_SN=(z_SN-z0)*nprocz/Lz
    iproc_SN=ipz_SN*nprocy + ipy_SN
    !if (lroot) print*, 'z',fran3(3),z_SN,ipz_SN,iproc_SN
!
!  Broadcast position to all processors from root;
!  also broadcast iproc_SN, needed for later broadcast of rho_SN.
!
    fmpi3=(/ x_SN, y_SN, z_SN /)
    call mpibcast_real(fmpi3,3)
    x_SN=fmpi3(1); y_SN=fmpi3(2); z_SN=fmpi3(3)
!
    impi=iproc_SN
    call mpibcast_int(impi,1)
    iproc_SN=impi(1)
!
!  With current SN scheme, we need rho at the SN location.
!
    rho_SN=0.0
    lfound=.false.; lfoundx=.false.; lfoundy=.false.; lfoundz=.false.
    if (iproc==iproc_SN) then
      do l=l1,l2; if (abs(x(l)-x_SN) < 1e-6) then; l_SN=l; lfoundx=.true.; endif; enddo
      do m=m1,m2; if (abs(y(m)-y_SN) < 1e-6) then; m_SN=m; lfoundy=.true.; endif; enddo
      do n=n1,n2; if (abs(z(n)-z_SN) < 1e-6) then; n_SN=n; lfoundz=.true.; endif; enddo
      !print*,'l,m,n_SN',l_SN,m_SN,n_SN
      lfound=(lfoundx .and. lfoundy .and. lfoundz)
      if (.not. lfound) print*,'position_SNI: SN not found!'
      !print*, 'position_SNI:',l_SN,m_SN,n_SN,ilnrho,f(l_SN,m_SN,n_SN,ilnrho)
      lnrho_SN=f(l_SN,m_SN,n_SN,ilnrho)
      rho_SN=exp(lnrho_SN)
! calculate TT_SN here, for later use in explode_SN
      ss_SN=f(l_SN,m_SN,n_SN,ient)
! NEED TO USE IONISATION CALCS 

!      call ioncalc(lnrho,ss,cs2,TT1,cp1tilde, &
!        Temperature=TT_SN_new,InternalEnergy=ee,IonizationFrac=yH)
      TT_SN=cs20*exp(gamma1*(lnrho_SN-lnrho0)+gamma*ss_SN)/gamma1*cp1
      if (lroot.and.ip<14) &
           print*, 'position_SNI:',l_SN,m_SN,n_SN,x_SN,y_SN,z_SN,rho_SN,TT_SN
    endif
!
!  Broadcast rho_SN, TT_SN to all processors.
!
    fmpi2=(/ rho_SN, TT_SN /)
    call mpibcast_real_nonroot(fmpi2,2,iproc_SN)
    rho_SN=fmpi2(1); TT_SN=fmpi2(2)
    if (lroot.and.ip<14) &
         print*, 'position_SNI:',iproc_SN,x_SN,y_SN,z_SN,rho_SN,TT_SN
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
    real, dimension(mx,my,mz,mvar) :: f
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
      if(headtt) print*,'position_SNI'
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
      print*,'position_SNII, mass_cloud_byproc:',mass_cloud_byproc
      print*,'position_SNII, cum_prob_byproc:',cum_prob_byproc
      print*,'position_SNII, mass_cloud:',mass_cloud
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
    if (lroot) print*, 'position_SNII, franSN(1),iproc_SN:',franSN(1),iproc_SN
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
               ss=f(l,m,n,ient)
               TT=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss)/gamma1*cp1
               if (rho >= rho_crit .and. TT <= TT_crit) then
                 cum_mass=cum_mass+rho
                 cum_prob_onproc=cum_mass/mass_cloud
                 if (franSN(1) <= cum_prob_onproc) then
                    x_SN=x(l); y_SN=y(m); z_SN=z(n); rho_SN=rho; TT_SN=TT
                    print*,'position_SNII,cum_mass,cum_prob_onproc,franSN(1):', &
                                          cum_mass,cum_prob_onproc,franSN(1)
                    exit find_SN
                 endif
               endif
             enddo
           enddo
         enddo find_SN
         print*,'position_SNII (iproc_SN):',iproc,x_SN,y_SN,z_SN,rho_SN
     endif
    fmpi5=(/ x_SN, y_SN, z_SN, rho_SN, TT_SN /)
    call mpibcast_real_nonroot(fmpi5,5,iproc_SN)
    x_SN=fmpi5(1); y_SN=fmpi5(2); z_SN=fmpi5(3)
    rho_SN=fmpi5(4); TT_SN=fmpi5(5)
!
    if (lroot) &
         print*, 'position_SNII (root):',iproc,x_SN,y_SN,z_SN,rho_SN,TT_SN
!
    endsubroutine position_SNII
!***********************************************************************
    subroutine explode_SN(f,itype_SN)
      !
      !  Implement SN (of either type), at pre-calculated position
      !  (This can all be made more efficient, after debugging.)
      !
      !  ??-nov-02/grs : coded from GalaxyCode                        
      !  20-may-03/ajwm: pencil formulation and broken into subroutines
      !
      use Cdata
      use Mpicomm
      use Ionization
      !
      real, intent(inout), dimension(mx,my,mz,mvar) :: f
      integer, intent(in) :: itype_SN

      real :: width_SN,width_shell_outer,width_shell_inner,c_SN
      real :: mass_check=0., mass_shell=0., profile_integral
      real :: EE_SN=0.,EE2_SN=0.,rho_SN_new,TT_SN_new,dv

      real, parameter :: TT_limit=1.e7
     ! real, parameter :: TT_limit=1.e4    ! make weaker, for debug
      
      integer, parameter :: point_width=4
     ! integer :: point_width=8            ! make larger, for debug
      real, dimension(nx) :: deltarho, deltaEE, e_old, rho_old, e_new
      real, dimension(4) :: fmpi4, fmpi4_tmp
      real, dimension(nx) ::  lnrho_old, ss_old,TT1, cs2, cp1tilde ! not really wanted but sideeffect of
                                                                   ! using the (No)Ionisation module
                                                                   ! thermodynamics function
      integer:: move_mass=0      
          
      !
      !  identifier
      !
      if(lroot.and.ip<14) print*,'explode_SN, itype_SN:',itype_SN
      !
      width_SN=point_width*dxmin      
      c_SN=ampl_SN/(cnorm_SN*(width_SN**3))      !normalision for SN profile
      if (lroot.and.ip<14) print*,'width_SN,c_SN', width_SN,c_SN
        
      !
      !  Now deal with (if nec.) mass relocation
      !

      TT_SN_new=TT_SN+(c_SN*rho_SN*gamma*cp1)
      if(lroot) print*,'explode_SN - TT_SN, TT_SN_new :',TT_SN,TT_SN_new

      if (TT_SN_new < TT_limit) then
         ! move_mass=1  
         move_mass=0  !ngrs: switch off for debug...

         ! The bit that BREAKS the pencil formulation...
         ! must know the total moved mass BEFORE attempting mass relocation 
         ! could be put into a conditional
         call calcmassprofileintegral_SN(profile_integral,width_SN)      
         if (lroot.and.ip<14) print*, 'moving mass...'
      endif
!     if (lroot) print*, &
!           'explode_SN, TT_SN_new,TT_limit:',TT_SN_new,TT_limit
      

      do n=n1,n2
         do m=m1,m2
            call proximity_SN()
            if (move_mass.eq.1) then
               rho_SN_new=c_SN/TT_limit*TTunits
               mass_shell=(rho_SN_new-rho_SN)*profile_integral
               call makecavity_SN(deltarho,width_SN,rho_SN_new-rho_SN,mass_shell,mass_check)
            endif

            call injectenergy_SN(deltaEE,width_SN,ampl_SN,EE_SN)
            ! Apply perturbations
            ! (\delta s)/cp = (1/gamma)*ln(1+(\delta e / e)) - (gamma1/gamma)*ln(1+(\delta rho / rho))

            lnrho_old=f(l1:l2,m,n,ilnrho)
            rho_old=exp(lnrho_old)
            ss_old=f(l1:l2,m,n,ient)

            call thermodynamics(lnrho_old,ss_old,cs2,TT1,cp1tilde, InternalEnergy=e_old)

!            cs2=cs20*exp(gamma1*(lnrho_old-lnrho0)+gamma*ss_old)
!            e_new = cs20 * exp(gamma*f(l1:l2,m,n,ient))*((rho_old/rho0)**(gamma-2.))/(gamma*gamma1)

!            if (abs(sum(e_old-e_new)).gt.1.) then
!               print*,e_old-e_new
!            endif

!ajwm - 4\pi ??? What's that all about?...
            if (move_mass.eq.1) f(l1:l2,m,n,ilnrho)=alog(rho_old+deltarho)
            f(l1:l2,m,n,ient)=ss_old + &
                 ( alog(1.+ (deltaEE  &                           ! / 12.56637061
                            *(rho_old+deltarho) / e_old)) ) / gamma  
            if (move_mass.eq.1) f(l1:l2,m,n,ient) = f(l1:l2,m,n,ient) &
                                          - (gamma1*alog(1 + deltarho / rho_old) )/ gamma

            call thermodynamics(f(l1:l2,m,n,ilnrho),f(l1:l2,m,n,ient),cs2,TT1,cp1tilde,InternalEnergy=e_new)

            EE2_SN=EE2_SN+sum((e_new*exp(f(l1:l2,m,n,ilnrho)))-(e_old*rho_old))

            if (sum(deltaEE).gt.1.) then
               print*, sum(e_new*exp(f(l1:l2,m,n,ilnrho))),sum(e_old*rho_old),sum(e_new*exp(f(l1:l2,m,n,ilnrho)))-sum(e_old*rho_old),sum(deltaEE)
            endif
       enddo
      enddo
      

      ! Sum and share diagnostics etc. amongst processors
      dv=1.
      if (nxgrid/=1) dv=dv*dx
      if (nygrid/=1) dv=dv*dy
      if (nzgrid/=1) dv=dv*dz

      fmpi4_tmp=(/ mass_check, EE_SN, mass_shell, EE2_SN /)
      call mpireduce_sum(fmpi4_tmp,fmpi4,4) 
      call mpibcast_real(fmpi4,4)
      mass_check=fmpi4(1)*dv
      EE_SN=fmpi4(2)*dv; 
      EE2_SN=fmpi4(4)*dv; 
      mass_shell=fmpi4(3)*dv; 

      if (lroot.and.ip<14) print*, &
           'explode_SN, masses:',mass_shell,mass_check

      ! if (lroot.and.ip<14) print*,'profile_check', &
      !      (profile_check)**(1./3.),cnorm_SN
     
      if (lroot) then
         open(1,file=trim(datadir)//'/time_series.dat',position='append')
         write(1,'(a,1e11.3," ",i1," ",i2," ",1e11.3," ",1e11.3," ",1e11.3," ",1e11.3," ",1e11.3," ",1e11.3,a)')  &
              '#ExplodeSN: (t,type,iproc,x,y,z,rho,energy)=(', &
              t,itype_SN,iproc_SN,x_SN,y_SN,z_SN,rho_SN,EE_SN,EE2_SN, ')'
         write(6,'(a,1e11.3," ",i1," ",i2," ",1e11.3," ",1e11.3," ",1e11.3," ",1e11.3," ",1e11.3," ",1e11.3,a)')  &
              '#ExplodeSN: (t,type,iproc,x,y,z,rho,energy)=(', &
              t,itype_SN,iproc_SN,x_SN,y_SN,z_SN,rho_SN,EE_SN,EE2_SN, ')'
         close(1)
      endif
      
    endsubroutine explode_SN

!***********************************************************************
    subroutine calcmassprofileintegral_SN(profile_integral,width_SN)
!
!  Calculate integral of mass cavity profile  
!
!  22-may-03/ajwm: coded
!
!
      use CData

      real, intent(in) :: width_SN
      real :: dx_SN_in,dx_SN_out_x0,dx_SN_out_x1,dy_SN_in,dy_SN_out_y
      real :: dy_SN_out_x0a,dy_SN_out_x0b,dy_SN_out_x1a,dy_SN_out_x1b
      real :: dz_SN,yshift
      real, intent(out) :: profile_integral
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
            profile_integral = profile_integral + sum(exp(-min((dr2_SN(:)/width_SN**2)**3, 75.)))
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
    subroutine makecavity_SN(deltarho,width_SN,depth,mass_shell,mass_check)
      use CData
      !  (Can't do all of this within a single pencil-sweep, as I need a global
      !   sum of the relocated mass, for the current scheme.)
      !  (Would be much simpler if mass relocation was abandoned -- could do
      !   the following on pencils, and wouldn't need to communicate rho_SN, TT_SN)
      real, intent(in) :: width_SN, depth, mass_shell
      real, intent(inout) :: mass_check
      !
      real, dimension(nx) :: profile_shell_outer,profile_shell_inner, profile_cavity
      real, dimension(1) :: fsum1,fsum1_tmp
      real :: width_shell_outer, width_shell_inner, c_shell
      ! real, dimension(1) :: fsum1,fsum1_tmp
      real :: cnorm_SN=3.128289613    ! (int exp(-(r/a)^6) 4\pi r^2 dr) / (a^3)
                                      ! for gaussian of width a

      real, intent(out), dimension(nx) :: deltarho

      width_shell_outer=2.0*width_SN
      width_shell_inner=1.5*width_SN

      deltarho(:) =  depth*exp(-min((dr2_SN(:)/width_SN**2)**3, 75.))
      
      c_shell=mass_shell /                                  &
           (cnorm_SN*((width_shell_outer**3)-(width_shell_inner**3)))

!      if (lroot) print*, &
!           'explode_SN, c_shell:',c_shell
      !  add missing mass back into shell
      
      profile_shell_outer(:)=                              &
           exp(-min((dr2_SN(:)/width_shell_outer**2)**3, 75.))
      profile_shell_inner(:)=                              &
           exp(-min((dr2_SN(:)/width_shell_inner**2)**3, 75.))
      
      deltarho(:)=deltarho(:) + c_shell *      &
           (profile_shell_outer(:) - profile_shell_inner(:))
      
      mass_check=mass_check + c_shell *                    &
           sum(profile_shell_outer(:) - profile_shell_inner(:))
      
    endsubroutine makecavity_SN

!***********************************************************************
    subroutine injectenergy_SN(deltaEE,width_SN,energy_SN,EE_SN)
      use CData
      !  (Can't do all of this within a single pencil-sweep, as I need a global
      !   sum of the relocated mass, for the current scheme.)
      !  (Would be much simpler if mass relocation was abandoned -- could do
      !   the following on pencils, and wouldn't need to communicate rho_SN, TT_SN)
      real, intent(in) :: width_SN, energy_SN
      !
      real, dimension(nx) :: profile_SN
      real, intent(inout) :: EE_SN
      real, intent(out), dimension(nx) :: deltaEE

      real :: c_SN
      !
      ! Whether mass moved or not, inject energy.
      !

      c_SN=energy_SN/(cnorm_SN*(width_SN**3))      !normalision for SN profile
!      if (lroot) print*,'width_SN,c_SN', width_SN,c_SN

      profile_SN=exp(-min((dr2_SN(:)/width_SN**2)**3, 75.))

      deltaEE(:)=c_SN*profile_SN(:) ! spatial energy density 
      EE_SN=EE_SN+sum(c_SN*profile_SN(:))   ! EE in (code) erg, not erg/g!

    endsubroutine injectenergy_SN

endmodule interstellar
