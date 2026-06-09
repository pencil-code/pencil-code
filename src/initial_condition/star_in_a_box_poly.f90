! $Id$
!
!  Piecewise polytopic initial condition for star-in-a-box models
!
!  11-apr-23/pete: coded
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!** AUTOMATIC REFERENCE-LINK.TEX GENERATION ********************
! Declare relevant citations from pencil-code/doc/citations/ref.bib for this module.
! The entries are taken from pencil-code/doc/citations/notes.tex
!
! 2006ApJ...638..336D,% Dobler, Stix & Brandenburg "Magnetic Field Generation in Fully Convective Rotating Spheres"
!
!***************************************************************
!
module InitialCondition
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use EquationOfState
  use Sub, only: step, der_step
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: star_luminosity=1.0, Rstar=1.0
  real :: r1_poly, width_npoly, max_rheat = 1.0
  real :: npoly1=1.5, npoly2=3., npoly_jump=1.0, nad=1.5, K0=1.0
  real :: npolyK1=1.5, npolyK2=3., wheat=0.1, r1_K=1.0, width_K=1.0, K_floor=0.0
  character (len=labellen) :: Kproftype='none'
  character (len=labellen) :: strat_type='polytropic'
!
  namelist /initial_condition_pars/ &
      star_luminosity, Rstar, nad, npoly1, npoly2, npolyK1, npolyK2, &
      r1_poly, width_npoly, npoly_jump, strat_type, wheat, r1_K, width_K, &
      K_floor, Kproftype, max_rheat
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  07-may-09/wlad: coded
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_all(f,profiles)
!
!  Initializes all the f arrays in one call. This subroutine is called last.
!
!  03-mar-23/pete: coded
!
      use SharedVariables, only: get_shared_variable
      use EquationOfState, only: rho0, cs20, get_gamma_etc
      use General, only: safe_character_assign
      use FArrayManager
      use Sub, only: poly, register_report_aux
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray), intent(inout):: f
      real, dimension (nx,*),             optional, intent(out)  :: profiles
!
      real, dimension (2*nxgrid) :: rr_sph, npoly, gnpoly, gg_r, dTdr_global, dlnTdr_global, prof
      real, dimension (2*nxgrid) :: rho_global, TT_global, heat, lumi
      real, dimension (2*nxgrid) :: hcond_global, ghcond_global
      real, dimension (2*nxgrid) :: FF, GG, FF_prime, GG_prime, ggg_r
      real, dimension (2*nxgrid) :: lnrho_global, dlnrhodr_global, cs2_global, ss_global
      real, dimension (2*nxgrid) :: tmp1, tmp2, trans, K_env
      real :: rmax, drr, rr_tmp, delr=0.01, gamma, cv
      real :: hcond_tmp, ghcond_tmp, q, rr1
      real, target :: cs2cool=0.0
      real, dimension (:,:), pointer :: cpot
      real, dimension (:,:), pointer :: cpot2
      real, dimension (:), pointer :: g_r
      real, pointer :: r1_pot1, r0_pot, n_pot, g0
      integer :: ir, i, j, n, m, ix, ierr, nr, jr
      integer, parameter :: unit=1
      logical, pointer :: lss_running_aver, lhcond_global, lcool_prof_as_var
!      character (len=labellen), dimension(:), pointer :: ipotential
!
!  Retrieve stuff needed from gravity and entropy modules
!
      call get_shared_variable('cpot', cpot)
      call get_shared_variable('cpot2', cpot2)
      call get_shared_variable('g_r', g_r)
      call get_shared_variable('r1_pot1',r1_pot1)
      call get_shared_variable('r0_pot',r0_pot)
      call get_shared_variable('n_pot',n_pot)
      call get_shared_variable('g0',g0)
      call get_shared_variable('lss_running_aver',lss_running_aver)
      call get_shared_variable('lhcond_global',lhcond_global)
      call get_shared_variable('lcool_prof_as_var',lcool_prof_as_var)
      call get_gamma_etc(gamma,cv=cv)
!
!  Compute rr_sph, i.e., half-diagonal with radius going from zero to maximum
!  value at the corners. Add a little bit of safety margin with the last term.
!
      rmax = sqrt((3./4.)*Lx**2) + 0.1*Lx ! Assuming the domain is a cube!
      nr = 2*nxgrid
      do ir=1,nr
         rr_sph(ir) = rmax*real(ir-1)/real(nr-1)
      enddo
      drr = rr_sph(2) - rr_sph(1)
!
!  Decide wheter cpot or cpot2 is used to calculate gg_r
!
      if (any(cpot2 /= 0.)) then
         gg_r = poly( (/   cpot2(2,1),  2*cpot2(3,1),  3*cpot2(4,1),  4*cpot2(5,1),  &
                         5*cpot2(6,1),  6*cpot2(7,1),  7*cpot2(8,1),  8*cpot2(9,1),  &
                         9*cpot2(10,1),10*cpot2(11,1),11*cpot2(12,1),12*cpot2(13,1), &
                        13*cpot2(14,1),14*cpot2(15,1),15*cpot2(16,1) /), rr_sph)
      else
!
!  Compute gravitational acceleration (does not yet work for flattened potentials)
!      
         gg_r = - rr_sph * poly( (/ 2*(cpot(1,1)*cpot(4,1)-cpot(2,1)), &
                     3*(cpot(1,1)*cpot(5,1)-cpot(3,1)), &
                     4*cpot(1,1)*cpot(3,1), &
                     cpot(5,1)*cpot(2,1)-cpot(3,1)*cpot(4,1), &
                     2*cpot(2,1)*cpot(3,1), &
                     cpot(3,1)**2  /), rr_sph) &
                     / poly( (/ 1., 0., cpot(4,1), cpot(5,1), &
                     cpot(3,1) /), rr_sph)**2
!
!  Compute derivative of gg_r for later use
!
!  cpot = a0, a2, a3, b2, b3 
!
         FF = poly( (/ 2*(cpot(1,1)*cpot(4,1)-cpot(2,1)), &
                     3*(cpot(1,1)*cpot(5,1)-cpot(3,1)), &
                     4*cpot(1,1)*cpot(3,1), &
                     cpot(5,1)*cpot(2,1)-cpot(3,1)*cpot(4,1), &
                     2*cpot(2,1)*cpot(3,1), &
                     cpot(3,1)**2  /), rr_sph)
!
         FF_prime = poly( (/ 3*(cpot(1,1)*cpot(5,1)-3*cpot(3,1)), &
                     8*cpot(1,1)*cpot(3,1), &
                     3*(cpot(2,1)*cpot(5,1) - cpot(3,1)*cpot(4,1)), &
                     8*cpot(2,1)*cpot(3,1), &
                     5*cpot(3,1)**2 /), rr_sph)
!
         GG = poly( (/ 1., 0., cpot(4,1), cpot(5,1), &
                     cpot(3,1) /), rr_sph)
         GG_prime = poly( (/ 0., 2*cpot(4,1), 3*cpot(5,1), &
                     4*cpot(3,1) /), rr_sph)

         ggg_r = (GG**2*(FF + rr_sph*FF_prime) - 2*rr_sph*FF*GG*GG_prime)/(GG**4)
!
      endif
!
      if (g0 .ne. 0.) then
         gg_r= -g0*r1_pot1**(3*n_pot)*rr_sph**(n_pot-1)*(r0_pot**n_pot*(r1_pot1**(2*n_pot) &
                         + rr_sph**(2*n_pot))**(.5) + (r1_pot1*rr_sph)**n_pot)**((-1-n_pot)/n_pot) &
                         / (r1_pot1**(2*n_pot) + rr_sph**(2*n_pot))**((-1+2*n_pot)/(2*n_pot))
      endif
!
!  Compute profile of "polytropic index" npoly and its gradient for later use
!
      npoly = npoly1 + (npoly2-npoly1)*step(rr_sph,r1_poly,width_npoly)
!      gnpoly = (npoly2-npoly1)*der_step(rr_sph,r1_poly,width_npoly)
!
!  Set up temperature gradient
!
      dTdr_global = gg_r/(cv*(gamma-1)*(npoly+1))
      where (rr_sph >= Rstar) dTdr_global=0.
!
!  Integrate temperature from the known value at the centre to the surface
!
     TT_global(1)=cs20/(cv*gamma*(gamma-1.))
     do j=2,nr
        TT_global(j)=TT_global(j-1)+0.5*(dTdr_global(j)+dTdr_global(j-1))*drr
     enddo
!
!  Compute gradient of logarithmic temperature
!
     dlnTdr_global=dTdr_global/TT_global
!
!  Gradient of logarithmic density assuming hydrostatic equilibrium
!
     dlnrhodr_global=-dlnTdr_global+gg_r/(cv*(gamma-1)*TT_global)
!
!  Integrate temperature from the known value at the centre to the surface
!    
     lnrho_global(1)=log(rho0)
     do j=2,nr
        lnrho_global(j)=lnrho_global(j-1)+0.5*(dlnrhodr_global(j)+dlnrhodr_global(j-1))*drr
     enddo
!
!  Compute speed of sound and specific entropy
!
     cs2_global=TT_global*cv*gamma*(gamma-1.)
     ss_global=log(cs2_global/cs20)/gamma - & 
             (gamma-1.)/(gamma)*(lnrho_global-log(rho0))
!
!  Now obsolete Gaussian-kernel mapping for lnrho and ss for reference
!
!     do ix=1,mx
!        do iy=1,my
!           do iz=1,mz
!              rr_tmp = sqrt(x(ix)**2+y(iy)**2+z(iz)**2)
!              prof = exp(-(rr_sph-rr_tmp)**2/(2*delr**2))/(delr*sqrt(2*pi))
!              f(ix,iy,iz,ilnrho) = rmax*sum(prof*lnrho_global)/real(nr)
!              f(ix,iy,iz,iss) = rmax*sum(prof*ss_global)/real(nr)
!           enddo
!        enddo
!     enddo
!
     cs2cool=TT_global(nr)*cv*gamma*(gamma-1.)
     if (lroot) print*,'cs2cool = ',cs2cool
!
!  Integrate luminosity from the heating function
!
!PJK: inactive for the time being
     ! heat = 4*pi*rr_sph**2*(star_luminosity/(2*pi*wheat**2)**1.5)*exp(-rr_sph**2/(2.*wheat**2))
     ! lumi(1)=0.
     ! do j=2,nr
     !    lumi(j)=lumi(j-1)+0.5*(heat(j)+heat(j-1))*drr
     ! enddo
!
!  Compute heat conductivity and its gradient (re-use npoly here)
!
     if (lhcond_global) then
          select case (Kproftype)
!
          case ('picewise-poly') ! piecewise polytropic
            if (lroot) print*, 'initial_condition_all: piecewise polytropic K-profile'
            npoly = npolyK1 + (npolyK2-npolyK1)*step(rr_sph,r1_poly,width_npoly)
!            gnpoly = (npolyK2-npolyK1)*der_step(rr_sph,r1_poly,width_npoly)
            hcond_global = -star_luminosity/(4.*pi*(rr_sph+tini)**2*(gg_r+tini))*cv*(gamma-1)*(npoly+1)
            hcond_global(1) = hcond_global(2)
!
          case ('const+gaussian') ! K constant in the core, Gaussian decay outside
            if (lroot) print*, 'initial_condition_all: const+Gaussian K-profile'
            hcond_global = (1. - K_floor)*exp(-((rr_sph-r1_K)/width_K)**2) + K_floor
            where (rr_sph < r1_K) hcond_global = 1.
!
          case default
!
!  Catch unknown values
!
            if (lroot) print*, 'initial_condition_all: '//&
                 'No such value for Kproftype: ', trim(Kproftype)
            call stop_it("")
!
          endselect
!
!  Luminosity corresponding to the radiative flux from hcond_global profile
!
        lumi = -4.*pi*rr_sph**2*hcond_global*dTdr_global
!
!  Normalize hcond_global such that max radiative luminosity matches the
!  stellar luminosity
!
        K0 = star_luminosity/maxval(lumi)
        hcond_global = max(K0*hcond_global,tini) ! avoid values smaller than tini
!
!  Compute derivative of lnK, for the time being just 2nd order finite difference
!
        do j=2,nr-1
           ghcond_global(j) = (log(hcond_global(j+1)) - log(hcond_global(j-1)))/(2.*drr)
        enddo
        ghcond_global(1)=0.
        ghcond_global(nr)=0.
!
!  Re-compute luminosity with normalized hcond_global and calculate the
!  corresponding heating function.
!
        lumi = -4.*pi*rr_sph**2*hcond_global*dTdr_global
        do j=2,nr-1
           heat(j) = (lumi(j+1)-lumi(j-1))/(2.*drr)
        enddo
        heat = heat/(4*pi*(rr_sph+tini)**2)
!
!  Assume that heating and radiation are in equilibrium in the core where
!  the luminosity increases as a function of radius and take care of
!  boundary values.
!
        where (heat<0.) heat=0.
        heat(1)=0.
        heat(nr)=0.
!
!  Divide by luminosity because the profile is multiplied by it in entropy.f90
!
        heat = heat/star_luminosity
!
     endif
!
!  Use interpolation to map the radial profiles to the Cartesian grid.
!
     do ix=1,mx
        do iy=1,my
           do iz=1,mz
              rr_tmp = sqrt(x(ix)**2 + y(iy)**2 + z(iz)**2)

              if (rr_tmp <= rr_sph(1)) then
                 f(ix,iy,iz,ilnrho) = lnrho_global(1)
                 f(ix,iy,iz,iss) = ss_global(1)

                 if (lhcond_global) then
                    hcond_tmp  = hcond_global(1)
                    ghcond_tmp = ghcond_global(1)
                 endif

                 if (lcool_prof_as_var) f(ix,iy,iz,icool_prof) = heat(1)

              elseif (rr_tmp >= rr_sph(nr)) then
                 f(ix,iy,iz,ilnrho) = lnrho_global(nr)
                 f(ix,iy,iz,iss) = ss_global(nr)

                 if (lhcond_global) then
                    hcond_tmp  = hcond_global(nr)
                    ghcond_tmp = ghcond_global(nr)
                 endif

                 if (lcool_prof_as_var) f(ix,iy,iz,icool_prof) = heat(nr)

              else
                 jr = int((rr_tmp - rr_sph(1))/drr) + 1
                 jr = max(1, min(jr, nr-1))
                 q = (rr_tmp - rr_sph(jr))/drr

                 f(ix,iy,iz,ilnrho) = (1. - q)*lnrho_global(jr) + q*lnrho_global(jr+1)
                 f(ix,iy,iz,iss) = (1. - q)*ss_global(jr) + q*ss_global(jr+1)

                 if (lhcond_global) then
                    hcond_tmp  = (1. - q)*hcond_global(jr)  + q*hcond_global(jr+1)
                    ghcond_tmp = (1. - q)*ghcond_global(jr) + q*ghcond_global(jr+1)
                 endif

                 if (lcool_prof_as_var) &
                    f(ix,iy,iz,icool_prof) = (1. - q)*heat(jr) + q*heat(jr+1)
              endif

              if (lhcond_global) then
                 f(ix,iy,iz,iglobal_hcond) = hcond_tmp

                 rr1 = 1.0/max(rr_tmp,tini)
                 f(ix,iy,iz,iglobal_glhc  ) = ghcond_tmp*x(ix)*rr1
                 f(ix,iy,iz,iglobal_glhc+1) = ghcond_tmp*y(iy)*rr1
                 f(ix,iy,iz,iglobal_glhc+2) = ghcond_tmp*z(iz)*rr1
              endif

              if (lcool_prof_as_var) then
                 ! Ensure manually that no additional heating occurs above max_rheat
                 if (rr_tmp > max_rheat) f(ix,iy,iz,icool_prof)=0.
              endif
           enddo
        enddo
     enddo
!
!  Set ss_run_aver equal to the initial ss.
!
     if (lss_running_aver) then
        f(:,:,:,iss_run_aver) = f(:,:,:,iss)
     endif

    endsubroutine initial_condition_all
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat)
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
