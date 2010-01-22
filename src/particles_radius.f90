! $Id$
!
!  This module takes care of everything related to particle radius.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 1
! CPARAM logical, parameter :: lparticles_radius=.true.
!
!***************************************************************
module Particles_radius
!
  use Cdata
  use Messages
  use Particles_cdata
  use Particles_sub
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_radius.h'
!
  real :: vthresh_sweepup=-1.0, deltavp12_floor=0.0
  real, dimension (ninit) :: ap0=0.0
  real :: tstart_sweepup_par=0.0, cdtps=0.2, cdtpc=0.2
  real :: tstart_condensation_par=0.0
  real :: apmin=0.0, latent_heat_SI=2.257e6, alpha_cond=1.0, alpha_cond1=1.0
  real :: diffusion_coefficient=1.0, diffusion_coefficient1=1.0
  real :: tau_damp_evap=0.0, tau_damp_evap1=0.0
  logical :: lsweepup_par=.true., lcondensation_par=.false.
  logical :: llatent_heat=.true.
  character (len=labellen), dimension(ninit) :: initap='nothing'
  character (len=labellen) :: condensation_coefficient_type='constant'
!
  namelist /particles_radius_init_pars/ &
      initap, ap0, rhops, vthresh_sweepup, deltavp12_floor, &
      lsweepup_par, lcondensation_par, tstart_sweepup_par, cdtps, apmin, &
      condensation_coefficient_type, alpha_cond, diffusion_coefficient, &
      tau_damp_evap, llatent_heat, cdtpc
!
  namelist /particles_radius_run_pars/ &
      rhops, vthresh_sweepup, deltavp12_floor, &
      lsweepup_par, lcondensation_par, tstart_sweepup_par, cdtps, apmin, &
      condensation_coefficient_type, alpha_cond, diffusion_coefficient, &
      tau_damp_evap, llatent_heat, cdtpc
!
  integer :: idiag_apm=0, idiag_ap2m=0, idiag_apmin=0, idiag_apmax=0
  integer :: idiag_dvp12m=0, idiag_dtsweepp=0
!
  contains
!***********************************************************************
    subroutine register_particles_radius()
!
!  Set up indices for access to the fp and dfp arrays.
!
!  22-aug-05/anders: coded
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Index for particle radius.
!
      iap=npvar+1
!
!  Increase npvar accordingly.
!
      npvar=npvar+1
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_particles: npvar > mpvar','')
      endif
!
    endsubroutine register_particles_radius
!***********************************************************************
    subroutine initialize_particles_radius(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  22-aug-05/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
!  Calculate the number density of bodies within a superparticle.
!
      if (npart_radii>1 .and. &
          (.not. lcartesian_coords .or. lparticles_nbody .or. &
          lparticles_number .or. lparticles_spin)) then 
        call fatal_error('initialize_particles_radius: npart_radii > 1','')
      else
        mp_swarm=4/3.0*pi*rhops*ap0(1)**3
        if (lroot) print*, 'initialize_particles_radius: '// &
            'mass per dust grain mp_swarm=', mp_swarm
      endif
!
!  Short hand for spherical particle prefactor.
!
      four_pi_rhops_over_three=4*pi*rhops/3.0
!
!  Inverse coefficients.
!
      alpha_cond1=1/alpha_cond
      diffusion_coefficient1=1/diffusion_coefficient
      if (tau_damp_evap/=0.0) tau_damp_evap1=1/tau_damp_evap
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_radius
!***********************************************************************
    subroutine set_particle_radius(f,fp,npar_low,npar_high,init)
!
!  Set radius of new particles.
!
!  18-sep-09/nils: adapted from init_particles_radius
!
      use General, only: random_number_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: npar_low,npar_high
      logical, optional :: init
      logical :: initial
      real :: radius_fraction
!
      integer :: j,ind
!
      initial=.false.
      if (present(init)) then
        if (init) initial=.true.
      endif
!
      do j=1,ninit
        select case (initap(j))
!
        case ('nothing')
          if (initial.and.lroot.and.j==1) &
              print*, 'set_particles_radius: nothing'
!
        case ('constant')
          if (initial.and.lroot) print*, 'set_particles_radius: constant radius'
          call random_number_wrapper(radius_fraction)
          ind=ceiling(npart_radii*radius_fraction)
          fp(npar_low:npar_high,iap)=ap0(ind)
!
        endselect
      enddo
!
      call keep_compiler_quiet(f)
!
    endsubroutine set_particle_radius
!***********************************************************************
    subroutine pencil_criteria_par_radius()
!
!  All pencils that the Particles_radius module depends on are specified here.
!
!  21-nov-06/anders: coded
!
      if (lsweepup_par) then
        lpenc_requested(i_uu)=.true.
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_cc)=.true.
      endif
      if (lcondensation_par) then
        lpenc_requested(i_cs2)=.true.
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_ppvap)=.true.
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_cc)=.true.
        lpenc_requested(i_cc1)=.true.
        lpenc_requested(i_np)=.true.
        lpenc_requested(i_rhop)=.true.
        if (ltemperature) then
          lpenc_requested(i_cv1)=.true.
          lpenc_requested(i_TT1)=.true.
        endif
      endif
!
    endsubroutine pencil_criteria_par_radius
!***********************************************************************
    subroutine dap_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle radius.
!
!  22-aug-05/anders: coded
!
      use Particles_number
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
      logical :: lfirstcall=.true., lheader
!
      intent (in) :: f, fp
      intent (out) :: dfp
!
!  Print out header information in first time step.
!
      lheader=lfirstcall.and.lroot.and.(.not.lpencil_check_at_work)
!
!  Identify module.
!
      if (lheader) print*,'dap_dt_pencil: Calculate dap/dt'
!
      if (lsweepup_par) &
          call dap_dt_sweepup_pencil(f,df,fp,dfp,p,ineargrid)
      if (lcondensation_par) &
          call dap_dt_condensation_pencil(f,df,fp,dfp,p,ineargrid)
!
      lfirstcall=.false.
!
    endsubroutine dap_dt_pencil
!***********************************************************************
    subroutine dap_dt_sweepup_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Increase in particle radius due to sweep-up of small grains in the gas.
!
!  22-aug-05/anders: coded
!
      use Diagnostics, only: max_mn_name
      use Particles_number
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (nx) :: dt1_sweepup
      real :: deltavp
      integer :: k, ix0, ix
!
      intent (in) :: f, fp
      intent (out) :: dfp
!
!  Increase in particle radius due to sweep-up of small grains in the gas.
!
      if (t>=tstart_sweepup_par) then
!
        if (lfirst.and.ldt) dt1_sweepup=0.0
!
        if (npar_imn(imn)/=0) then
          do k=k1_imn(imn),k2_imn(imn)
            ix0=ineargrid(k,1)
            ix=ix0-nghost
!  No interpolation needed here.
!  Relative speed.
            deltavp=sqrt( &
                (fp(k,ivpx)-p%uu(ix,1))**2 + &
                (fp(k,ivpy)-p%uu(ix,2))**2 + &
                (fp(k,ivpz)-p%uu(ix,3))**2 )
            if (deltavp12_floor/=0.0) &
                deltavp=sqrt(deltavp**2+deltavp12_floor**2)
!  Allow boulders to sweep up small grains if relative velocity not too high.
            if (deltavp<=vthresh_sweepup .or. vthresh_sweepup<0.0) then
              if (.not. lpscalar) then
                call fatal_error('dap_dt', &
                    'must have passive scalar module for sweep-up')
              endif
!  Radius increase due to sweep-up.
              dfp(k,iap) = dfp(k,iap) + 0.25*deltavp*p%cc(ix)*p%rho(ix)*rhops1
!
!  Deplete gas of small grains.
!
              if (lparticles_number) np_swarm=fp(k,inpswarm)
              if (lpscalar_nolog) then
                df(ix0,m,n,icc) = df(ix0,m,n,icc) - &
                    np_swarm*pi*fp(k,iap)**2*deltavp*p%cc(ix)
              else
                df(ix0,m,n,ilncc) = df(ix0,m,n,ilncc) - &
                    np_swarm*pi*fp(k,iap)**2*deltavp
              endif
!
!  Time-step contribution of sweep-up.
!
              if (lfirst.and.ldt) then
                dt1_sweepup(ix) = dt1_sweepup(ix) + &
                    np_swarm*pi*fp(k,iap)**2*deltavp
              endif
!
            endif
!
            if (ldiagnos) then
              if (idiag_dvp12m/=0) call sum_par_name((/deltavp/),idiag_dvp12m)
            endif
          enddo
        endif
!
!  Time-step contribution of sweep-up.
!
          if (lfirst.and.ldt) then
            dt1_sweepup=dt1_sweepup/cdtps
            dt1_max=max(dt1_max,dt1_sweepup)
            if (ldiagnos.and.idiag_dtsweepp/=0) &
                call max_mn_name(dt1_sweepup,idiag_dtsweepp,l_dt=.true.)
          endif
!
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine dap_dt_sweepup_pencil
!***********************************************************************
    subroutine dap_dt_condensation_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Change in particle radius due to condensation of monomers from the gas and
!  evaporation of solid/liquid droplets.
!
!  15-jan-10/anders: coded
!
      use Particles_number
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (nx) :: ap_equi, vth, dt1_condensation, rhovap
      real, dimension (nx) :: total_surface_area, ppsat, supsatratio1
      real, dimension (nx) :: rhocond_tot, rhosat, np_total
      real :: dapdt, drhocdt, alpha_cond_par
      integer :: k, ix, ix0
!
      intent (in) :: f, fp
      intent (out) :: dfp
!
!  Change in particle radius due to condensation and evaporation.
!
      if (t>=tstart_condensation_par) then
!
        if (npar_imn(imn)/=0) then
          rhovap=p%cc*p%rho
          vth=sqrt(p%cs2)   ! Should actually be thermal speed of vapor
          ppsat=6.035e11*exp(-5938*p%TT1)  ! Valid for water
          supsatratio1=ppsat/p%ppvap
          rhosat=supsatratio1*rhovap
          rhocond_tot=p%rhop+rhovap
          if (lfirst.and.ldt) then
            np_total=0.0
            total_surface_area=0.0
            dt1_condensation=0.0
          endif
          do k=k1_imn(imn),k2_imn(imn)
            ix0=ineargrid(k,1)
            ix=ix0-nghost
            if (.not. lpscalar) then
              call fatal_error('dap_dt', &
                  'must have passive scalar module for condensation')
            endif
!
!  The condensation/evaporation mass flux is
!
!    F = vth*rhovap*alpha
!
!  where vth is the thermal speed, rhopvap is the mass density of vapor,
!  and alpha \in [0,1] is the condensation coefficient. Various authors argue
!  for different choices of alpha. Following Barkstrom 1978 we take
!
!    alpha = 1/(vth*ap/D + 1/alpha0)
!
!  where D is the diffusion coefficient of vapor and ap the particle radius.
!  Small particles, with ap<D/(alpha0*vth), have alpha=alpha0, while the
!  condenation coefficient decreases linearly with particle size for larger
!  particles. Barkstrom 1978 use alpha0=0.04.
!
            select case (condensation_coefficient_type)
              case ('constant')
                alpha_cond_par=alpha_cond
              case ('size-dependent')
                alpha_cond_par=1/(vth(ix)*fp(k,iap)*diffusion_coefficient1+ &
                    alpha_cond1)
              case default
                if (lroot) print*, 'dap_dt_condensation_pencil: '// &
                    'invalid condensation coefficient type'
                call fatal_error('dap_dt_condensation_pencil','')
            endselect
!
!  Radius increase by condensation or decrease by evaporation.
!
            if (fp(k,iap)<apmin) then
!
!  Do not allow particles to become smaller than a minimum radius.
!
              dapdt=-tau_damp_evap1*(fp(k,iap)-apmin)
              if (lfirst .and. ldt) then
                dt1_condensation(ix)=max(dt1_condensation(ix),tau_damp_evap1)
              endif
            else
              dapdt=0.25*vth(ix)*rhovap(ix)*rhops1*(1.0-supsatratio1(ix))* &
                  alpha_cond_par
!
!  Damp approach to minimum size. The radius decreases linearly with time in
!  the limit of small particles; therefore we need to damp the evaporation to
!  avoid small time-steps.
!
              if (dapdt<0.0) then
                dapdt=dapdt*min(1.0,(fp(k,iap)/apmin-1.0)**2)
                if (lfirst .and. ldt) then
                  dt1_condensation(ix)=max(dt1_condensation(ix), &
                      abs(dapdt/(fp(k,iap)-apmin)))
                endif
              endif
            endif
!
            dfp(k,iap)=dfp(k,iap)+dapdt
!
!  Vapor monomers are added to the gas or removed from the gas.
!
            if (lparticles_number) np_swarm=fp(k,inpswarm)
            if (lfirst.and.ldt) np_total(ix)=np_total(ix)+np_swarm
            drhocdt=-dapdt*4*pi*fp(k,iap)**2*rhops*np_swarm
            if (lpscalar_nolog) then
              df(ix0,m,n,icc)  =df(ix0,m,n,icc)   + p%rho1(ix)*drhocdt
            else
              df(ix0,m,n,ilncc)=df(ix0,m,n,ilncc) + p%rho1(ix)*drhocdt*p%cc1(ix)
            endif
!
!  Release latent heat to gas / remove heat from gas.
!
            if (ltemperature.and.llatent_heat) then
              df(ix0,m,n,ilnTT)=df(ix0,m,n,ilnTT) - &
                  latent_heat_SI*p%rho1(ix)*p%TT1(ix)*p%cv1(ix)*drhocdt
            endif
!
            if (lfirst.and.ldt) total_surface_area(ix)=total_surface_area(ix)+ &
                4*pi*fp(k,iap)**2*np_swarm*alpha_cond_par
          enddo
!
!  Time-step contribution of condensation.
!
          if (lfirst.and.ldt) then
            ap_equi=((p%rhop+(1.0-supsatratio1)*rhovap)/ &
                (4.0/3.0*pi*rhops*np_swarm*p%np))**(1.0/3.0)
            do ix=1,nx
              if (rhocond_tot(ix)>rhosat(ix)) then
                dt1_condensation(ix) = max(total_surface_area(ix)*vth(ix), &
                    pi*vth(ix)*np_total(ix)*ap_equi(ix)**2*alpha_cond)
              endif
            enddo
          endif
        endif
!
!  Time-step contribution of condensation.
!
        if (lfirst.and.ldt) then
          dt1_condensation=dt1_condensation/cdtpc
          dt1_max=max(dt1_max,dt1_condensation)
        endif
!
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine dap_dt_condensation_pencil
!***********************************************************************
    subroutine dap_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle radius.
!
!  21-nov-06/anders: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
!  Diagnostic output.
!
      if (ldiagnos) then
        if (idiag_apm/=0) call sum_par_name(fp(1:npar_loc,iap),idiag_apm)
        if (idiag_ap2m/=0) call sum_par_name(fp(1:npar_loc,iap)**2,idiag_ap2m,lsqrt=.true.)
        if (idiag_apmin/=0) &
            call max_par_name(-fp(1:npar_loc,iap),idiag_apmin,lneg=.true.)
        if (idiag_apmax/=0) &
            call max_par_name(fp(1:npar_loc,iap),idiag_apmax)
      endif
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dap_dt
!***********************************************************************
    subroutine read_particles_rad_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
      integer :: i
!
      if (present(iostat)) then
        read(unit,NML=particles_radius_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_radius_init_pars,ERR=99)
      endif
!
!  Find how many different particle radii we are using. This must be done
!  because not all parts of the code are adapted to work with more than one
!  particle radius.
!
      do i=1,ninit
        if (ap0(i)/=0) then
          npart_radii=npart_radii+1
        endif
      enddo
!
99    return
!
    endsubroutine read_particles_rad_init_pars
!***********************************************************************
    subroutine write_particles_rad_init_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_radius_init_pars)
!
    endsubroutine write_particles_rad_init_pars
!***********************************************************************
    subroutine read_particles_rad_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_radius_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_radius_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_rad_run_pars
!***********************************************************************
    subroutine write_particles_rad_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_radius_run_pars)
!
    endsubroutine write_particles_rad_run_pars
!***********************************************************************
    subroutine rprint_particles_radius(lreset,lwrite)
!
!  Read and register print parameters relevant for particles radius.
!
!  22-aug-05/anders: coded
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
      logical :: lwr
!
!  Write information to index.pro.
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      if (lwr) write(3,*) 'iap=', iap
!
!  Reset everything in case of reset.
!
      if (lreset) then
        idiag_apm=0; idiag_ap2m=0; idiag_apmin=0; idiag_apmax=0
        idiag_dvp12m=0; idiag_dtsweepp=0
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot.and.ip<14) &
          print*, 'rprint_particles_radius: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'apm',idiag_apm)
        call parse_name(iname,cname(iname),cform(iname),'ap2m',idiag_ap2m)
        call parse_name(iname,cname(iname),cform(iname),'apmin',idiag_apmin)
        call parse_name(iname,cname(iname),cform(iname),'apmax',idiag_apmax)
        call parse_name(iname,cname(iname),cform(iname),'dvp12m',idiag_dvp12m)
        call parse_name(iname,cname(iname),cform(iname),'dtsweepp',idiag_dtsweepp)
      enddo
!
    endsubroutine rprint_particles_radius
!***********************************************************************
endmodule Particles_radius
