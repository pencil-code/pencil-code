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
  real :: r1_poly, width_npoly
  real :: npoly1=1.5, npoly2=3., npoly_jump=1.0, nad=1.5
  real :: npolyK1=1.5, npolyK2=3., wheat=0.1
  character (len=labellen) :: strat_type='polytropic'
!
  namelist /initial_condition_pars/ &
      star_luminosity, Rstar, nad, npoly1, npoly2, npolyK1, npolyK2, &
      r1_poly, width_npoly, npoly_jump, strat_type, wheat
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
      use SharedVariables, only: get_shared_variable !, put_shared_variable
      use EquationOfState, only: rho0, cs20, get_gamma_etc
      use General, only: safe_character_assign
      use FArrayManager
      use Sub, only: poly, register_report_aux
!
      real, dimension (mx,my,mz,mfarray), optional, intent(inout):: f
      real, dimension (nx,*),             optional, intent(out)  :: profiles
!
      real, dimension (2*nxgrid) :: rr_sph, npoly, gnpoly, gg_r, dTdr_global, dlnTdr_global, prof
      real, dimension (2*nxgrid) :: rho_global, TT_global, heat, lumi
      real, dimension (2*nxgrid) :: hcond_global, ghcond_global
      real, dimension (2*nxgrid) :: FF, GG, FF_prime, GG_prime, ggg_r
      real, dimension (2*nxgrid) :: lnrho_global, dlnrhodr_global, cs2_global, ss_global
      real, dimension (2*nxgrid) :: tmp1, tmp2
      real :: rmax, drr, rr_tmp, delr=0.01, gamma, cv
      real, target :: cs2cool=0.0
      real, dimension (:,:), pointer :: cpot
      real, dimension (:,:), pointer :: cpot2
      real, dimension (:), pointer :: g_r
      real, pointer :: r1_pot1, r0_pot, n_pot, g0
      integer :: ir, i, j, n, m, ix, ierr
      integer :: iglobal_hcond, iglobal_glhc
      integer, parameter :: unit=1
      character (len=labellen), dimension(:), pointer :: ipotential
!
!     Retrieve cv and stuff from gravity_r
!
      call get_shared_variable('cpot', cpot)
      call get_shared_variable('cpot2', cpot2)
      call get_shared_variable('g_r', g_r)
      call get_shared_variable('r1_pot1',r1_pot1)
      call get_shared_variable('r0_pot',r0_pot)
      call get_shared_variable('n_pot',n_pot)
      call get_shared_variable('g0',g0)
      call get_gamma_etc(gamma,cv=cv)
!
!  Compute rr_sph (i.e. half-diagonal with radius going from zero to maximum value at the corners)
!
      rmax = sqrt((3./4.)*Lx**2) + 0.1*Lx ! Assuming the domain is a cube!
      do ir=1,2*nxgrid
         rr_sph(ir) = rmax*(ir-1)/float(2*nxgrid)
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
!      print*,'gg_r =', gg_r
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
     do j=2,2*nxgrid
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
     do j=2,2*nxgrid
        lnrho_global(j)=lnrho_global(j-1)+0.5*(dlnrhodr_global(j)+dlnrhodr_global(j-1))*drr
     enddo
!
!  Compute speed of sound and specific entropy
!
     cs2_global=TT_global*cv*gamma*(gamma-1.)
     ss_global=log(cs2_global/cs20)/gamma - & 
             (gamma-1.)/(gamma)*(lnrho_global-log(rho0))
!
!  Use a Gaussian kernel to set lnrho and ss at the positions of the grid points.
!
     do ix=1,mx
        do iy=1,my
           do iz=1,mz
              rr_tmp = sqrt(x(ix)**2+y(iy)**2+z(iz)**2)
              prof = exp(-(rr_sph-rr_tmp)**2/(2*delr**2))/(delr*sqrt(2*pi))
              f(ix,iy,iz,ilnrho) = rmax*sum(prof*lnrho_global)/(2.*nxgrid)
              f(ix,iy,iz,iss) = rmax*sum(prof*ss_global)/(2.*nxgrid)
           enddo
        enddo
     enddo
!
     cs2cool=TT_global(2*nxgrid)*cv*gamma*(gamma-1.)
     if (lroot) print*,'cs2cool = ',cs2cool
!
!  Integrate luminosity from the heating function
!
     heat = 4*pi*rr_sph**2*(star_luminosity/(2*pi*wheat**2)**1.5)*exp(-rr_sph**2/(2.*wheat**2))
     lumi(1)=0.
     do j=2,2*nxgrid
        lumi(j)=lumi(j-1)+0.5*(heat(j)+heat(j-1))*drr
     enddo
!
!  Compute heat conductivity and its gradient (re-use npoly here)
!
      npoly = npolyK1 + (npolyK2-npolyK1)*step(rr_sph,r1_poly,width_npoly)
      gnpoly = (npolyK2-npolyK1)*der_step(rr_sph,r1_poly,width_npoly)
!
      hcond_global = -lumi/(4*pi*(rr_sph+1e-15)**2*(gg_r+1e-15))*cv*(gamma-1)*(npoly+1)
      hcond_global(1) = hcond_global(2)
      hcond_global(2) = hcond_global(3) + .5*(hcond_global(1)-hcond_global(3))
      
!      do j=1,2*nxgrid
!         print*, 'rr_sph(j), hcond_global(j) = ',rr_sph(j),hcond_global(j)
!      enddo

      ghcond_global = -cv*(gamma-1.)/(4*pi)* &
           ( (4*pi*rr_sph**2*heat*(npoly+1) + lumi*gnpoly)*rr_sph*gg_r - lumi*(npoly+1)*(2*gg_r + rr_sph*ggg_r) ) &
           / (rr_sph**3*gg_r**2)
      ghcond_global(1) = 0.
      ghcond_global(2) = 0.5*ghcond_global(3)

!      do j=1,2*nxgrid
!         print*, 'rr_sph(j), hcond_global(j),ghcond_global(j) = ',rr_sph(j),hcond_global(j),ghcond_global(j)
!         print*, 'rr_sph(j), gnpoly(j) = ',rr_sph(j),gnpoly(j)
!      enddo
!
!  Put hcond and ghcond to f-array
!
!  Heat conductivity and its gradient.
!
!      do ix=1,mx
!         do iy=1,my
!            do iz=1,mz
!               rr_tmp = sqrt(x(ix)**2+y(iy)**2+z(iz)**2)
!               prof = exp(-(rr_sph-rr_tmp)**2/(2*delr**2))/(delr*sqrt(2*pi))
!               f(ix,iy,iz,iglobal_hcond) = rmax*sum(prof*hcond_global)/(2.*nxgrid)
! !              print*,'f(ix,iy,iz,iglobal_hcond) =', f(ix,iy,iz,iglobal_hcond)
! !              f(ix,iy,iz,iglobal_glhc) = rmax*sum(prof*ghcond_global)/(2.*nxgrid)
!            enddo
!         enddo
!      enddo
!      print*,'initial condition, f(5,5,5,iglobal_hcond) =', f(5,5,5,iglobal_hcond)
!
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
