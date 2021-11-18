! $Id$
!
!  This module provide a way for users to specify custom initial
!  conditions.
!
!  The module provides a set of standard hooks into the Pencil Code
!  and currently allows the following customizations:
!
!   Description                               | Relevant function call
!  ------------------------------------------------------------------------
!   Initial condition registration            | register_special
!     (pre parameter read)                    |
!   Initial condition initialization          | initialize_special
!     (post parameter read)                   |
!                                             |
!   Initial condition for momentum            | initial_condition_uu
!   Initial condition for density             | initial_condition_lnrho
!   Initial condition for entropy             | initial_condition_ss
!   Initial condition for magnetic potential  | initial_condition_aa
!
!   And a similar subroutine for each module with an "init_XXX" call.
!   The subroutines are organized IN THE SAME ORDER THAT THEY ARE CALLED.
!   First uu, then lnrho, then ss, then aa, and so on.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
! HOW TO USE THIS FILE
! --------------------
!
! Change the line above to
!    linitial_condition = .true.
! to enable use of custom initial conditions.
!
! The rest of this file may be used as a template for your own initial
! conditions. Simply fill out the prototypes for the features you want
! to use.
!
! Save the file with a meaningful name, e.g. mhs_equilibrium.f90, and
! place it in the $PENCIL_HOME/src/initial_condition directory. This
! path has been created to allow users to optionally check their
! contributions in to the Pencil Code SVN repository. This may be
! useful if you are working on/using an initial condition with
! somebody else or may require some assistance from one from the main
! Pencil Code team. HOWEVER, less general initial conditions should
! not go here (see below).
!
! You can also place initial condition files directly in the run
! directory. Simply create the folder 'initial_condition' at the same
! level as the *.in files and place an initial condition file there.
! With pc_setupsrc this file is linked automatically into the local
! src directory. This is the preferred method for initial conditions
! that are not very general.
!
! To use your additional initial condition code, edit the
! Makefile.local in the src directory under the run directory in which
! you wish to use your initial condition. Add a line that says e.g.
!
!    INITIAL_CONDITION =   initial_condition/mhs_equilibrium
!
! Here mhs_equilibrium is replaced by the filename of your new file,
! not including the .f90 extension.
!
! This module is based on Tony's special module.
!
module InitialCondition
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use EquationOfState
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: g0=1.,density_power_law=0.,temperature_power_law=1., plasma_beta=25
  real :: truncation_degree=2.,truncation_scale=1.
  logical :: lexponential_smooth=.false.
  real :: radial_percent_smooth=10.0,rshift=0.0
  real :: gravitational_const=0.
  real :: magnetic_power_law=impossible
  real :: dustdensity_powerlaw=1.5,edtog=0.01
!
! For the magnetic field
!
  logical :: ladd_field=.true.
  character (len=labellen) :: initcond_aa=''
  real, dimension(3) :: B_ext=(/0.0,0.0,0.0/)
  real :: rmode_mag=4.,zmode_mag=16.
  real :: rm_int=-impossible,rm_ext=impossible
  real :: Bz_const=8d-3
  real, dimension(0:6) :: coeff_cs2=(/0.,0.,0.,0.,0.,0.,0./)
!
! For the noise
! 
  real :: rho_rms=0.05
  integer :: xmodes=10,ymodes=10,zmodes=0
  logical :: llowk_noise=.false.,lgaussian_distributed_noise=.true.
  real :: xmid=1.5,rborder_int=0.,rborder_ext=0.,Lxn=1.
!
  real :: r0_pot=0.,qgshear=1.5
  integer :: n_pot=10
!
! Correct extra forces in the centrifugal balance condition
!
  logical :: lcorrect_pressuregradient=.true.
  logical :: lcorrect_selfgravity=.false.
  logical :: lcorrect_lorentzforce=.false.
  logical :: lpolynomial_fit_cs2=.false.
  logical :: ladd_noise_propto_cs=.false.
  logical :: lselfgravity_logspirals=.false.
  real :: ampluu_cs_factor=1d-3
  real :: widthbb1=0.0,widthbb2=0.0
  real :: OOcorot
!
  real :: bump_radius = 1.,bump_ampl = 0.4, bump_width = 0.1
  character (len=labellen) :: ipressurebump='nobump'
  character (len=labellen) :: imidplane='power-law'
!
! Magnetic spiral
  real :: B0_spiral=0.012042837784031205
  real :: etamu0_spiral = 1.0
  real :: Omega0_spiral=1.0, r0_spiral=1.0
!
  namelist /initial_condition_pars/ g0,density_power_law,&
       truncation_degree,truncation_scale,temperature_power_law,&
       lexponential_smooth,radial_percent_smooth,rshift,&
       lcorrect_selfgravity,gravitational_const,xmodes,ymodes,zmodes,&
       rho_rms,llowk_noise,xmid,lgaussian_distributed_noise,&
       rborder_int,rborder_ext,plasma_beta,ladd_field,initcond_aa,B_ext,&
       zmode_mag,rmode_mag,rm_int,rm_ext,Bz_const,r0_pot,qgshear,n_pot,&
       magnetic_power_law,lcorrect_lorentzforce,lcorrect_pressuregradient,&
       lpolynomial_fit_cs2,ladd_noise_propto_cs,ampluu_cs_factor,widthbb1,&
       widthbb2,bump_radius,bump_ampl,bump_width,ipressurebump,imidplane,&
       lselfgravity_logspirals,dustdensity_powerlaw,edtog,&
       B0_spiral,etamu0_spiral,Omega0_spiral,r0_spiral
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
      use SharedVariables, only: put_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      Lxn=Lx-2*(rborder_int+rborder_ext)
!
      if (lcorotational_frame) then
        if (lpointmasses) then
          if (lroot) then
            print*,'initialize_initial_condition: lcorotational frame not coded for pointmasses'
            print*,'switch to gravity_r'
            call fatal_error("","")
          endif
        endif
        OOcorot=rcorot**(-1.5)
      else
        OOcorot=0.
      endif
!
      if (llocal_iso.and.lparticles_blocks) &
           call put_shared_variable('itemperature_power_law',&
             temperature_power_law,&
             caller='initialize_particles')
! for magnetic spiral
      call put_shared_variable('B0_spiral',B0_spiral)
      call put_shared_variable('etamu0_spiral',etamu0_spiral)
      call put_shared_variable('Omega0_spiral',Omega0_spiral)
      call put_shared_variable('r0_spiral',r0_spiral)
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  This subroutine is a general routine that takes
!  the gravity acceleration and adds the centrifugal force
!  that numerically balances it.
!
!  Pressure corrections to ensure centrifugal equilibrium are
!  added in the respective modules
!
!  24-feb-05/wlad: coded
!  04-jul-07/wlad: generalized for any shear
!  08-sep-07/wlad: moved here from initcond
!
      use Gravity, only: acceleration
      use Sub,     only: get_radial_distance,power_law
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: rr_cyl,rr_sph,OO,g_r,tmp
      integer :: i
!
      if (lroot) &
           print*,'centrifugal_balance: initializing velocity field'
!
      if ((rsmooth/=0.).or.(r0_pot/=0)) then
        if (rsmooth/=r0_pot) &
             call fatal_error("centrifugal_balance","rsmooth and r0_pot must be equal")
        if (n_pot<=2) &
             call fatal_error("centrifugal_balance","don't you dare using less smoothing than n_pot=2")
      endif
!
      if (.not.lread_oldsnap) then
      do m=m1,m2
      do n=n1,n2
!
        call get_radial_distance(rr_sph,rr_cyl)
!
        if (lgrav) then
!
! Gravity of a static central body
!
          call acceleration(g_r)
!
! Sanity check
!
          if (any(g_r > 0.)) then
            do i=l1,l2
              if (g_r(i-l1+1) > 0) then
                print*,"centrifugal_balance: gravity at physical point ",&
                     x(i),y(m),z(n),"is directed outwards"
                call fatal_error("","")
              endif
            enddo
          else
            if ( (coord_system=='cylindric')  .or.&
                 (coord_system=='cartesian')) then
              OO=sqrt(max(-g_r/rr_cyl,0.))
            else if (coord_system=='spherical') then
              OO=sqrt(max(-g_r/rr_sph,0.))
            endif
          endif
!
        elseif (lpointmasses) then
!
! Gravity from dynamical point masses with a dominating central body
!
          call power_law(sqrt(g0),rr_sph,qgshear,tmp)
!
          if (lcartesian_coords.or.&
               lcylindrical_coords) then
            OO=tmp
            if (lcylindrical_gravity) &
                 OO=tmp*sqrt(rr_sph/rr_cyl)
          elseif (lspherical_coords) then
            OO=tmp
          endif
!
        else
!
          print*,"both gravity and pointmasses are switched off"
          print*,"there is no gravity to determine the azimuthal velocity"
          call fatal_error("initial_condition_uu","")
!
        endif
!
        if (coord_system=='cartesian') then
          f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) - y(  m  )*OO
          f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + x(l1:l2)*OO
          f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + 0.
        elseif (coord_system=='cylindric') then
          f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) + 0.
          f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + (OO-OOcorot)*rr_cyl
          f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + 0.
        elseif (coord_system=='spherical') then
          f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) + 0.
          f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + 0.
          f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + (OO-OOcorot)*rr_sph
        endif
!
      enddo
      enddo
    endif
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine add_noise(f)
!
      use FArrayManager, only: farray_use_global
      use EquationOfState, only: get_cv1,cs20,gamma_m1,lnrho0
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: cs2
      real :: cv1,cp1
      integer, pointer :: iglobal_cs2
!
      if (llocal_iso) then
        call farray_use_global('cs2',iglobal_cs2)
      elseif (lentropy) then
        call get_cv1(cv1)
      elseif (ltemperature) then
        call get_cp1(cp1)
      endif
!
      do n=n1,n2; do m=m1,m2
        if (llocal_iso) then
          cs2=f(l1:l2,m,n,iglobal_cs2)
        elseif (lentropy) then
          !even with ldensity_nolog=T, this rho is in log
          cs2=cs20*exp(cv1*f(l1:l2,m,n,iss)+ & 
               gamma_m1*(f(l1:l2,m,n,ilnrho)-lnrho0))
        elseif (ltemperature) then
          if (ltemperature_nolog) then
            cs2=f(l1:l2,m,n,iTT)*gamma_m1/cp1
          else
            cs2=exp(f(l1:l2,m,n,ilnTT))*gamma_m1/cp1
          endif
        endif
        call gaunoise_vect(ampluu_cs_factor*sqrt(cs2),f,iux,iuz)
      enddo; enddo
!
    endsubroutine add_noise
!***********************************************************************
    subroutine gaunoise_vect(ampl,f,i1,i2)
!
!  Add Gaussian noise (= normally distributed) white noise for variables i1:i2
!
!  23-may-02/axel: coded
!  10-sep-03/axel: result only *added* to whatever f array had before
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i1,i2
!
      real, dimension (nx) :: r,p,tmp,ampl
      integer :: i,j
!
      intent(in)    :: ampl,i1,i2
      intent(inout) :: f
      integer, save, dimension(mseed) :: rstate=0
!
!  set gaussian random noise vector
!
      do i=i1,i2
        if (lroot.and.m==1.and.n==1) print*,'gaunoise_vect: variable i=',i
        if (modulo(i-i1,2)==0) then
           do j=1,nx
              r(j)=ran0(rstate(1))
              p(j)=ran0(rstate(1))
          enddo
          tmp=sqrt(-2*log(r))*sin(2*pi*p)
        else
          tmp=sqrt(-2*log(r))*cos(2*pi*p)
        endif
        f(l1:l2,m,n,i)=f(l1:l2,m,n,i)+ampl*tmp
      enddo
!
    endsubroutine gaunoise_vect
!***********************************************************************
    function ran0(dummy)
!
!  The 'Minimal Standard' random number generator
!  by Lewis, Goodman and Miller.
!
!  28.08.02/nils: Adapted from Numerical Recipes
!
      integer, intent(inout) :: dummy
!
      integer :: k
      integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836, &
           mask=123459876
      real, parameter :: am=1./im
      real :: ran0
!
      dummy=ieor(dummy,mask)
      k=dummy/iq
      dummy=ia*(dummy-k*iq)-ir*k
      if (dummy<0) dummy=dummy+im
      ran0=am*dummy
      dummy=ieor(dummy,mask)
!
    endfunction ran0
!***********************************************************************
    subroutine poly_fit(cs2)
!
      real, dimension(0:6) :: c
      real, dimension(mx) :: cs2,lncs2
!
!  Fits for different combinations of cs2 and temperature_power_law
!
      if (cs0 .eq. 0.1 .and. temperature_power_law.eq.2) then
         coeff_cs2=(/-8.33551,27.6856,-57.5702,54.6696,-27.2370,6.87119,-0.690690/)
      else if (cs0 .eq. 0.1 .and. temperature_power_law.eq.1) then
         coeff_cs2=(/-6.47454,13.8181,-28.6687,27.1693,-13.5113,3.40305,-0.341599/)
      else
         call fatal_error("poly_fit",&
              "fit not calculated for choice of cs0 and Teff power law")
      endif
      c=coeff_cs2
      lncs2 = c(0) + c(1) * x    + c(2) * x**2 + c(3) * x**3 + &
                     c(4) * x**4 + c(5) * x**5 + c(6) * x**6
!
      cs2=exp(lncs2)
!
    endsubroutine poly_fit
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  07-may-09/wlad: coded
!
      use FArrayManager
      use Gravity,         only: potential,acceleration
      use Sub,             only: get_radial_distance,grad,power_law
      use EquationOfState, only: rho0
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx)   :: strat,tmp1,tmp2,cs2
      real, dimension (mx)   :: rr_sph,rr,rr_cyl,lnrhomid,rho_bump
      real                   :: rmid,lat
      integer, pointer       :: iglobal_cs2,iglobal_glnTT
      integer                :: ics2
      logical                :: lheader,lpresent_zed
!
      if (lroot) print*,&
           'initial_condition_lnrho: locally isothermal approximation'
!
      if (lenergy.and.llocal_iso) then
        if (lroot) then
          print*,'You switched on entropy or temperature evolution,'
          print*,'but you are still using llocal_iso in start.in.'
          print*,'Use one or the other instead.'
        endif
        call fatal_error("initial_condition_lnrho","")
      endif
!
!  Set the sound speed
!
      do m=1,my; do n=1,mz
        lheader=((m==1).and.(n==1).and.lroot)
        call get_radial_distance(rr_sph,rr_cyl)
        if (lcylindrical_gravity.or.lcylinder_in_a_box.or.lcylindrical_coords) then
          rr=rr_cyl 
        elseif (lsphere_in_a_box.or.lspherical_coords) then
          rr=rr_sph
        else
          call fatal_error("initial_condition_lnrho",&
               "no valid coordinate system")
        endif
!
        if (llocal_iso.or.lenergy) then
          if (.not.lpolynomial_fit_cs2) then
            call power_law(cs20,rr,temperature_power_law,cs2,r_ref)
          else 
            call poly_fit(cs2)
          endif  

!
!  Store cs2 in one of the free slots of the f-array
!
          if (llocal_iso) then
            nullify(iglobal_cs2)
            call farray_use_global('cs2',iglobal_cs2)
            ics2=iglobal_cs2
            f(:,m,n,ics2)=cs2
          elseif (.not.lread_oldsnap) then
            if (ltemperature) then
              if (ltemperature_nolog) then
                ics2=iTT
              else
                ics2=ilnTT
              endif
              f(:,m,n,ics2)=cs2
            elseif (lentropy) then
              ics2=iss
              f(:,m,n,ics2)=cs2
            endif
          endif
        else
          ics2=impossible_int 
        endif
      enddo;enddo
!
!  Stratification is only coded for 3D runs. But as
!  cylindrical and spherical coordinates store the
!  vertical direction in different slots, one has to
!  do this trick below to decide whether this run is
!  2D or 3D.
!
      if (.not.lread_oldsnap) then
      lpresent_zed=.false.
      if (lspherical_coords) then
        if (nygrid/=1) lpresent_zed=.true.
      else
        if (nzgrid/=1) lpresent_zed=.true.
      endif
!
!  Pencilize the density allocation.
!
      do n=1,mz
      do m=1,my
!
        lheader=lroot.and.(m==1).and.(n==1)
!
!  Midplane density
!
        call get_radial_distance(rr_sph,rr_cyl)
        if (lcylindrical_gravity.or.lcylinder_in_a_box.or.lcylindrical_coords) then
          rr=rr_cyl
        elseif (lsphere_in_a_box.or.lspherical_coords) then
          rr=rr_sph
        endif
!
        select case (ipressurebump)
        case ('gaussian')
           rho_bump = 1 + bump_ampl*exp(-(rr_cyl - bump_radius)**2/(2*bump_width**2))
        case ('step')
           rho_bump = 1 + .5*bump_ampl*(tanh((rr_cyl-bump_radius)/bump_width) + 1)
        case ('nobump')    
           rho_bump = 1.
        case default
           if (lroot) print*, 'No such value for ipressurebump: ', trim(ipressurebump)
           call fatal_error("initial_condition_lnrho","")
        endselect
!
        select case (imidplane)
        case ('power-law')
           if (lheader) print*,'Radial stratification with power law=',&
                density_power_law
          if (lexponential_smooth) then
            !radial_percent_smooth = percentage of the grid
            !that the smoothing is applied
            rmid=rshift+(xyz1(1)-xyz0(1))/radial_percent_smooth
            lnrhomid=log(rho0) &
                 + density_power_law*log((1-exp( -((rr-rshift)/rmid)**2 ))/rr)
          else
            lnrhomid=log(rho0) -.5*density_power_law*log((rr/r_ref)**2+rsmooth**2)
          endif
        case ('exponential')
          if (lheader) print*,'Exponential disk'
          lnrhomid=log(rho0) - rr/r_ref
        case ('truncated')
          if (lheader) print*,'Truncated disk'
          lnrhomid=log(rho0) - density_power_law*log((rr/r_ref)) &
                   - (rr/(truncation_scale*r_ref))**(truncation_degree-density_power_law)
       case default
           if (lroot) print*, 'No such value for imidplane: ', trim(ipressurebump)
           call fatal_error("initial_condition_lnrho","")
        endselect
!
!  Vertical stratification, if needed
!
        if (.not.lcylindrical_gravity.and.lpresent_zed) then
          if (lheader) &
               print*,"Adding vertical stratification with "//&
               "scale height h/r=",cs0
!
!  Get the sound speed
!
          if (lenergy.or.llocal_iso) then
            cs2=f(:,m,n,ics2)
          else
            cs2=cs20
          endif
!
          if (lspherical_coords.or.lsphere_in_a_box) then
            ! uphi2/r = -gr + dp/dr
            if (lgrav) then
              call acceleration(tmp1)
            elseif (lpointmasses) then
              !call get_totalmass(g0) 
              tmp1=-g0/rr_sph**2
            else
              print*,"both gravity and pointmasses are switched off"
              print*,"there is no gravity to determine the stratification"
              call fatal_error("local_isothermal_density","")
            endif
!
            tmp2=-tmp1*rr_sph - &
                 cs2*(density_power_law + temperature_power_law)/gamma
            lat=pi/2-y(m)
            strat=(tmp2*gamma/cs2) * log(cos(lat))
          else
!
!  The subroutine "potential" yields the whole gradient.
!  I want the function that partially derived in
!  z gives g0/r^3 * z. This is NOT -g0/r
!  The second call takes care of normalizing it
!  i.e., there should be no correction at midplane
!
            if (lgrav) then
              call potential(POT=tmp1,RMN=rr_sph)
              call potential(POT=tmp2,RMN=rr_cyl)
            elseif (lpointmasses) then
              !call get_totalmass(g0)
              tmp1=-g0/rr_sph 
              tmp2=-g0/rr_cyl
            else
              print*,"both gravity and pointmasses are switched off"
              print*,"there is no gravity to determine the stratification"
              call fatal_error("local_isothermal_density","")
            endif
            strat=-(tmp1-tmp2)/cs2
            if (lenergy) strat=gamma*strat
          endif
!
        else
!  No stratification
          strat=0.
        endif
        f(:,m,n,ilnrho) = f(:,m,n,ilnrho) + lnrhomid + strat + log(rho_bump)
      enddo
      enddo
!
!  Correct the velocities by this pressure gradient
!
      if (lcorrect_pressuregradient) &
           call correct_pressure_gradient(f,ics2,temperature_power_law)
!
!  Correct the velocities for self-gravity
!
      if (lcorrect_selfgravity) then
        !if (lselfgravity_logspirals) then
        !  call correct_selfgravity_logspirals(f)
        !else
          call correct_selfgravity(f)
        !endif
      endif
!
!  Set noise in low wavelengths only
!
      if (llowk_noise) call lowk_noise_gaussian_rprof(f)
!
      endif
!
!  Set the thermodynamical variable
!
      if (llocal_iso) then
        call set_thermodynamical_quantities&
             (f,temperature_power_law,ics2,iglobal_cs2,iglobal_glnTT)
      else if (lenergy.and.(.not.lread_oldsnap)) then
        call set_thermodynamical_quantities(f,temperature_power_law,ics2)
      endif
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
    subroutine initial_condition_nd(f)
!
!  Initialize logarithmic density. init_lnrho 
!  will take care of converting it to linear 
!  density if you use ldensity_nolog
!
!  20-sep-17/wlad: coded
!
      use EquationOfState, only: rho0
      use Sub,             only: get_radial_distance
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: rr_sph,rr_cyl
      real, dimension (mx) :: lnrhomid !,strat
      !real, dimension (mx) :: cs2,tmp1,tmp2
      real    :: p !,q
      integer :: k
      !integer :: ics2
      !integer, pointer :: iglobal_cs2
!
      p=-dustdensity_powerlaw 
      if (lroot) print*,'Radial dust density stratification with power law=',p
!
!  Pencilize the density allocation.
!
      do n=1,mz
        do m=1,my
!
!  Midplane density
!
          call get_radial_distance(rr_sph,rr_cyl)
          lnrhomid=log(edtog*rho0)+p*log(rr_cyl/r_ref)
          do k=1,ndustspec
             f(:,m,n,ind(k)) = f(:,m,n,ind(k))+exp(lnrhomid)
          enddo
!
!  Dust stratification          
!
          !if (lgrav) then
          !  call potential(POT=tmp1,RMN=rr_sph)
          !  call potential(POT=tmp2,RMN=rr_cyl)
          !elseif (lpointmasses) then
          !  tmp1=-g0/rr_sph 
          !  tmp2=-g0/rr_cyl
          !endif
          !if (lcylindrical_gravity) then 
          !  strat=0.
          !else
          !  strat=-gamma*(tmp1-tmp2)/vrms2
          !endif
          !f(:,m,n,ilnrho) = f(:,m,n,ilnrho)+strat
!
        enddo
      enddo
!
    endsubroutine initial_condition_nd
!***********************************************************************
    subroutine initial_condition_uud(f)
!
!  Initialize the velocity field.
!
!  This subroutine is a general routine that takes
!  the gravity acceleration and adds the centrifugal force
!  that numerically balances it.
!
!  Pressure corrections to ensure centrifugal equilibrium are
!  added in the respective modules
!
!  24-feb-05/wlad: coded
!  04-jul-07/wlad: generalized for any shear
!  08-sep-07/wlad: moved here from initcond
!
      use Gravity, only: acceleration
      use Sub,     only: get_radial_distance,power_law
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: rr_cyl,rr_sph,OO,g_r,tmp
      integer :: i,k
!
      if (lroot) &
           print*,'centrifugal_balance: initializing velocity field'
!
      if ((rsmooth/=0.).or.(r0_pot/=0)) then
        if (rsmooth/=r0_pot) &
             call fatal_error("centrifugal_balance","rsmooth and r0_pot must be equal")
        if (n_pot<=2) &
             call fatal_error("centrifugal_balance","don't you dare using less smoothing than n_pot=2")
      endif
!
      do m=m1,m2
      do n=n1,n2
!
        call get_radial_distance(rr_sph,rr_cyl)
!
        if (lgrav) then
!
! Gravity of a static central body
!
          call acceleration(g_r)
!
! Sanity check
!
          if (any(g_r > 0.)) then
            do i=l1,l2
              if (g_r(i-l1+1) > 0) then
                print*,"centrifugal_balance: gravity at physical point ",&
                     x(i),y(m),z(n),"is directed outwards"
                call fatal_error("","")
              endif
            enddo
          else
            if ( (coord_system=='cylindric')  .or.&
                 (coord_system=='cartesian')) then
              OO=sqrt(max(-g_r/rr_cyl,0.))
            else if (coord_system=='spherical') then
              OO=sqrt(max(-g_r/rr_sph,0.))
            endif
          endif
!
        elseif (lpointmasses) then
!
! Gravity from dynamical point masses with a dominating central body
!
          call power_law(sqrt(g0),rr_sph,qgshear,tmp)
!
          if (lcartesian_coords.or.&
               lcylindrical_coords) then
            OO=tmp
            if (lcylindrical_gravity) &
                 OO=tmp*sqrt(rr_sph/rr_cyl)
          elseif (lspherical_coords) then
            OO=tmp
          endif
!
        else
!
          print*,"both gravity and pointmasses are switched off"
          print*,"there is no gravity to determine the azimuthal velocity"
          call fatal_error("initial_condition_uud","")
!
        endif
!
        do k=1,ndustspec 
        if (coord_system=='cartesian') then
          f(l1:l2,m,n,iudx(k)) = f(l1:l2,m,n,iudx(k)) - y(  m  )*OO
          f(l1:l2,m,n,iudy(k)) = f(l1:l2,m,n,iudy(k)) + x(l1:l2)*OO
          f(l1:l2,m,n,iudz(k)) = f(l1:l2,m,n,iudz(k)) + 0.
        elseif (coord_system=='cylindric') then
          f(l1:l2,m,n,iudx(k)) = f(l1:l2,m,n,iudx(k)) + 0.
          f(l1:l2,m,n,iudy(k)) = f(l1:l2,m,n,iudy(k)) + (OO-OOcorot)*rr_cyl
          f(l1:l2,m,n,iudz(k)) = f(l1:l2,m,n,iudz(k)) + 0.
        elseif (coord_system=='spherical') then
          f(l1:l2,m,n,iudx(k)) = f(l1:l2,m,n,iudx(k)) + 0.
          f(l1:l2,m,n,iudy(k)) = f(l1:l2,m,n,iudy(k)) + 0.
          f(l1:l2,m,n,iudz(k)) = f(l1:l2,m,n,iudz(k)) + (OO-OOcorot)*rr_sph
       endif
       enddo
!
      enddo
      enddo
!
    endsubroutine initial_condition_uud
!***********************************************************************
    subroutine set_thermodynamical_quantities&
         (f,temperature_power_law,ics2,iglobal_cs2,iglobal_glnTT)
!
!  Subroutine that sets the thermodynamical quantities
!   - static sound speed, temperature or entropy -
!  based on a sound speed which is given as input.
!  This routine is not general. For llocal_iso (locally
!  isothermal approximation, the temperature gradient is
!  stored as a static array, as the (analytical) derivative
!  of an assumed power-law profile for the sound speed
!  (hence the parameter temperature_power_law)
!
!  05-jul-07/wlad: coded
!  16-dec-08/wlad: moved pressure gradient correction to
!                  the density module (correct_pressure_gradient)
!                  Now this subroutine really only sets the thermo
!                  variables.
!
      use FArrayManager
      use EquationOfState, only: gamma,gamma_m1,get_cp1,&
           cs20,cs2bot,cs2top,lnrho0
      use Sub,             only: power_law,get_radial_distance
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: rr,rr_sph,rr_cyl,cs2,lnrho
      real, dimension (nx) :: gslnTT
      real :: cp1,temperature_power_law
      integer, pointer, optional :: iglobal_cs2,iglobal_glnTT
      integer :: ics2
!
      intent(in)    :: temperature_power_law
      intent(inout) :: f
!
!  Break if llocal_iso is used with entropy or temperature
!
      if (lenergy.and.llocal_iso) &
           call fatal_error('set_thermodynamical_quantities','You are '//&
           'evolving the energy, but llocal_iso is switched '//&
           ' on in start.in. Better stop and change it')
!
!  Break if gamma=1.0 and energy is solved
!
      if ((gamma==1.0).and.lenergy) then
        if (lroot) then
          print*,""
          print*,"set_thermodynamical_quantities: gamma=1.0 means "
          print*,"an isothermal disk. You don't need entropy or "
          print*,"temperature for that. Switch to noentropy instead, "
          print*,"which is a better way of having isothermality. "
          print*,"If you do not want isothermality but wants to keep a "
          print*,"static temperature gradient through the simulation, use "
          print*,"noentropy with the switch llocal_iso in init_pars "
          print*,"(start.in file) and add the following line "
          print*,""
          print*,"! MGLOBAL CONTRIBUTION 4"
          print*,""
          print*,"(containing the '!') to the header of the "//&
               "src/cparam.local file"
          print*,""
          call fatal_error('','')
        endif
      endif
!
      if (lroot) print*,'Temperature gradient with power law=',temperature_power_law
!
!  Get the pointers to the global arrays if needed
!
      if (llocal_iso) then
        nullify(iglobal_glnTT)
        call farray_use_global('glnTT',iglobal_glnTT)
      endif
!
      if (lenergy) call get_cp1(cp1)
!
      do m=m1,m2
      do n=n1,n2
!
!  Put in the global arrays if they are to be static
!
        cs2=f(l1:l2,m,n,ics2)
        if (llocal_iso) then
!
          f(l1:l2,m,n,iglobal_cs2) = cs2
!
          call get_radial_distance(rr_sph,rr_cyl);   rr=rr_cyl
          if (lspherical_coords.or.lsphere_in_a_box) rr=rr_sph
!
          gslnTT=-temperature_power_law/((rr/r_ref)**2+rsmooth**2)*rr/r_ref**2
!
          if (lcartesian_coords) then
            f(l1:l2,m,n,iglobal_glnTT  )=gslnTT*x(l1:l2)/rr_cyl
            f(l1:l2,m,n,iglobal_glnTT+1)=gslnTT*y(m)    /rr_cyl
            f(l1:l2,m,n,iglobal_glnTT+2)=0.
          else! (lcylindrical_coords.or.lspherical_coords) then
            f(l1:l2,m,n,iglobal_glnTT  )=gslnTT
            f(l1:l2,m,n,iglobal_glnTT+1)=0.
            f(l1:l2,m,n,iglobal_glnTT+2)=0.
          endif
        elseif (ltemperature) then
!  else do it as temperature ...
          if (ltemperature_nolog) then
            f(l1:l2,m,n,iTT)=cs2*cp1/gamma_m1
          else
            f(l1:l2,m,n,ilnTT)=log(cs2*cp1/gamma_m1)
          endif
        elseif (lentropy) then
!  ... or entropy
          lnrho=f(l1:l2,m,n,ilnrho) ! initial condition, always log
          f(l1:l2,m,n,iss)=1./(gamma*cp1)*(log(cs2/cs20)-gamma_m1*(lnrho-lnrho0))
       else
!
          call fatal_error('set_thermodynamical_quantities', &
               'No thermodynamical variable. Choose if you want '//&
               'a local thermodynamical approximation '//&
               '(switch llocal_iso=T init_pars and entropy=noentropy on '//&
               'Makefile.local), or if you want to compute the '//&
               'temperature directly and evolve it in time.')
        endif
      enddo
      enddo
!
!  Word of warning...
!
      if (lroot.and.llocal_iso) then
        if (associated(iglobal_cs2)) then
          print*,"Max global cs2 = ",&
               maxval(f(l1:l2,m1:m2,n1:n2,iglobal_cs2))
          print*,"Sum global cs2 = ",&
               sum(f(l1:l2,m1:m2,n1:n2,iglobal_cs2))
        endif
        if (associated(iglobal_glnTT)) then
          print*,"Max global glnTT(1) = ",&
               maxval(f(l1:l2,m1:m2,n1:n2,iglobal_glnTT))
          print*,"Sum global glnTT(1) = ",&
               sum(f(l1:l2,m1:m2,n1:n2,iglobal_glnTT))
        endif
      endif
!
      cs2bot=cs20
      cs2top=cs20
!
      if (lroot) &
           print*,"thermodynamical quantities successfully set"
!
!  Add noise if needed.
!
      if (ladd_noise_propto_cs) call add_noise(f)
!
    endsubroutine set_thermodynamical_quantities
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  07-may-09/wlad: coded
!
      use EquationOfState, only: cs20,get_cv1
      use Sub, only: step, power_law
      use FArrayManager, only: farray_use_global
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: amplbb,pblaw
      real, dimension(nx) :: Aphi, rr, rrcyl
      real, dimension(mx) :: bz,aphi_mx,Breal,pp
!
      real :: const,k1,kr
      integer :: i
!
      select case (initcond_aa)
!
!  Constant plasma beta in spherical coordinates. Sets Bphi
!
        case ('plasma-beta') 
          call set_field_constbeta(f)
        case ('Alfven-zconst')
!
!  Radially variable field pointing in the z direction
!  4 Balbus-Hawley wavelengths in the vertical direction
!
!  Bz=Lz/(8pi)*Omega      ==> Aphi = Lz/(8pi) Omega*r/(2-q)
!
!  The smoothed case should be general, since it reduces
!  to the non-smoothed for r0_pot=0.
!
!  B=C*(r2+r02)^-q ==> Aphi=C/(r*(2-q))*(r2+r02)^(1-q/2)
!
!  04-oct-06/wlad: coded
!
          if (lcartesian_coords) then
            amplbb = Lxyz(3)/(2*zmode_mag*pi)
            do n=n1,n2; do m=m1,m2
              rr=sqrt(x(l1:l2)**2+y(m)**2)
              Aphi=amplbb/(rr*(2-qgshear))*(rr**2+r0_pot**2)**(1-qgshear/2.)
              f(l1:l2,m,n,iax) =  -Aphi*y(m)/rr
              f(l1:l2,m,n,iay) =   Aphi*x(l1:l2)/rr
            enddo; enddo
          elseif (lcylindrical_coords) then
            call fatal_error("initial_condition_aa","alfven_zconst: "//&
                 "not implemented for cylindrical coordinates")
          elseif (lspherical_coords) then
            amplbb=Lxyz(2)/(2*zmode_mag*pi)
            pblaw=1-qgshear!-plaw/2.
            do n=n1,n2; do m=m1,m2
              rr=x(l1:l2)
              Aphi=-amplbb/(pblaw+2)*rr**(pblaw+1)*1
              f(l1:l2,m,n,iax)=0.
              f(l1:l2,m,n,iay)=0.
              f(l1:l2,m,n,iaz)=Aphi/sin(y(m))
            enddo; enddo
!
          endif
!
!--------------------------------------------------------------
!
        case ('alfven_rz')
!
!  Alfven wave propagating on radial direction with
!  field pointing to the z direction.
!
!  Bz = B0 cos(k r) ==> Aphi = B0/k sin(k r) + B0/(k2*r)*cos(k r)
!
!  04-oct-06/wlad: coded
!
          if (headtt) print*,'radial alfven wave propagating on z direction'
          if (.not.lcylinder_in_a_box) &
               call fatal_error("alfven_rz, this initial condition",&
               "works only for embedded cylinders")
!
! Choose between cases. The non-smoothed case is singular 
! in r=0 for the potential, thus can only be used if 
! freezing is used with a non-zero r_int. The smoothed case
! has a linear component that prevents singularities and a exponential 
! that prevents the field from growing large in the outer disk
!
          do n=n1,n2; do m=m1,m2
            rrcyl = max(sqrt(x(l1:l2)**2 + y(m)**2),tini)
            if (r_int>0.) then
              if (lroot .and. m==m1 .and. n==n1) then
                print*,'freezing is being used, ok to use singular potentials'
                print*,'Bz=B0cos(k.r) ==> Aphi=B0/k*sin(k.r)+B0/(k^2*r)*cos(k r)'
              endif
              kr = 2*pi*rmode_mag/(r_ext-r_int)
              Aphi =  Bz_const/kr * sin(kr*(rrcyl-r_int)) + &
                   Bz_const/(kr**2*rrcyl)*cos(kr*(rrcyl-r_int))
            else
              if (lroot .and. m==m1 .and. n==n1) &
                   print*,'Softened magnetic field in the center'
              if (rmode_mag < 5) call fatal_error("initial_condition_aa",&
                   "put more wavelengths in the field")
              kr = 2*pi*rmode_mag/r_ext
              k1 = 1. !not tested for other values
              const=Bz_const*exp(1.)*k1/cos(kr/k1)
              Aphi=const/kr*rrcyl*exp(-k1*rrcyl)*sin(kr*rrcyl)
            endif
!
            f(l1:l2,m,n,iax) = Aphi * (-    y(m)/rrcyl)
            f(l1:l2,m,n,iay) = Aphi * ( x(l1:l2)/rrcyl)
!
          enddo; enddo
!
!--------------------------------------------------------------
!
       case ("alfven_rphi")
         if (.not.lcylinder_in_a_box) &
              call fatal_error("alfven_rz, this initial condition",&
              "works only for embedded cylinders")
         do n=n1,n2; do m=m1,m2
           kr = 2*pi*rmode_mag/(r_ext-r_int)
           rrcyl = sqrt(x(l1:l2)**2 + y(m)**2)
           f(l1:l2,m,n,iaz) =  -Bz_const/kr*sin(kr*(rrcyl-r_int))
         enddo; enddo
!
!-------------------------------------------------------------
!
       case ("sine-avoid-boundary")
!
! Sine field in cylindrical coordinates, used in Armitage 1998
!
!   Bz=B0/r * sin(kr*(r-r0))
!
! And 0 outside of the interval r0-rn
! Code the field and find Aphi through solving the
! tridiagonal system for
!
!  Bz= d/dr Aphi + Aphi/r
!
!  -A_(i-1) + A_(i+1) + 2*A_i*dr/r = 2*dr*Bz
!
!  05-apr-08/wlad : coded
!
         if (.not.lcylindrical_coords) &
              call fatal_error("initial_condition_aa",&
              "this IC assumes cylindrical coordinates")
!
         do i=1,mx
           if ((rcyl_mn(i)>=rm_int).and.(rcyl_mn(i)<=rm_ext)) then
             kr = 2*pi*rmode_mag/(rm_ext-rm_int)
             bz(i)=Bz_const/rcyl_mn(i) * sin(kr*(rcyl_mn(i)-rm_int))
           else
             bz(i)=0.
           endif
         enddo
!
         call integrate_field(bz,aphi_mx)
         do m=1,my; do n=1,mz
           f(:,m,n,iay) = aphi_mx
         enddo;enddo
!
!-------------------------------------------------------------
!
       case ("bz-const")
!
         Breal = bz_const
         call set_field(f,Breal)
!
       case ("lambda_over_h_cte") 
         if (zmode_mag==0) &
              call fatal_error("initcond_aa","zmode_mag is zero")
         call power_law(rho0*cs20,x,&
              temperature_power_law+density_power_law,pp,r_ref)
         Breal = sqrt(pp/gamma)/(zmode_mag*2*pi)
         call set_field(f,Breal)
!
       case ("lambda_over_Lz_cte") 
         if (zmode_mag==0) &
              call fatal_error("initcond_aa","zmode_mag is zero")
         if (magnetic_power_law/=impossible) then
           pblaw=magnetic_power_law
         else
           pblaw=1.5+0.5*density_power_law !alfven
         endif
         call power_law(Lxyz(3)/(zmode_mag*2*pi),x,pblaw,Breal,r_ref)
         call set_field(f,Breal)
!
        case ('','nothing')
!
!  Do nothing
!           
        case default
!
!  Catch unknown values.
!
          call fatal_error('initcond_aa', &
              'initcond_aa value "' // trim(initcond_aa) // '" not recognised')
!
        endselect
!
    if (lcorrect_lorentzforce) call correct_lorentz_numerical(f)
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine set_field(f,Breal)
!
!  This routine sets a radially-dependent vertical field, by numerically 
!  integrating it to find the corresponding azimuthal magnetic potential. 
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx) :: bz,aphi_mx,Breal
!
      call cap_field(Breal,bz)
      call integrate_field(bz*x,aphi_mx)
!
      do n=1,mz;do m=1,my
        if (lcylindrical_coords) then
          f(:,m,n,iay) = aphi_mx/x 
        elseif (lspherical_coords) then
          f(:,m,n,iaz) = sin(y(m))*aphi_mx/x 
        else 
          call fatal_error("initial_condition_aa",&
               "bz-const not implemented for the chosen coordinated system")
        endif
      enddo;enddo
!
    endsubroutine set_field
!***********************************************************************
    subroutine set_field_constbeta(f)
!
      use FArrayManager,   only: farray_use_global
      use Sub,             only: get_radial_distance
      use EquationOfState, only: cs20,get_cv1,lnrho0,gamma_m1
      use Boundcond,       only: update_ghosts
      use Messages,        only: fatal_error
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx) :: cs2,pressure,Bphi,BB,tmp
      real, dimension(mx) :: rr_sph,rr_cyl
      integer, pointer :: iglobal_cs2
      real :: cv1
!
      if (.not.lspherical_coords) &
           call fatal_error("set_field_constbeta",&
           "coded for spherical coordinates only")
!
      call update_ghosts(f)
!
      if (llocal_iso) then
        call farray_use_global('cs2',iglobal_cs2)
      elseif (lentropy) then
        call get_cv1(cv1)
      endif

      do m=m1,m2 
        if (llocal_iso) then
          cs2=f(:,m,npoint,iglobal_cs2)
        elseif (lentropy) then
          cs2=cs20*exp(cv1*f(:,m,npoint,iss)+ & 
               gamma_m1*(log(f(:,m,npoint,irho))-lnrho0))
        endif
        pressure = f(:,m,npoint,irho)*cs2 /gamma
        Bphi = sqrt(2*pressure/plasma_beta)
        
        call cap_field(Bphi,BB)
!
!  Bphi = 1/r*d/dr(r*Atheta), so integrate: Atheta=1/r*Int(B*r)dr
!
        call get_radial_distance(rr_sph,rr_cyl)
        call integrate_field(BB*rr_sph,tmp)
        do n=n1,n2
          f(l1:l2,m,n,iay)=f(l1:l2,m,n,iay) + tmp(l1:l2)/rr_sph(l1:l2)
        enddo
      enddo
!
    endsubroutine set_field_constbeta
!***********************************************************************
    subroutine cap_field(Bin,Bout)
!
      use Sub, only: step
!
      real, dimension(mx), intent(in) :: Bin
      real, dimension(mx), intent(out) :: Bout
      integer :: i
!      
      if (rm_int==-impossible.and.rm_ext==impossible) then
        Bout=Bin
      else
        do i=1,mx
          if (widthbb1==0.0) widthbb1=5./dx_1(i)
          if (widthbb2==0.0) widthbb2=5./dx_1(i)
          Bout(i) = Bin(i) * &
               (step(x(i),rm_int,widthbb1)-&
                step(x(i),rm_ext,widthbb2))
        enddo
      endif
!
    endsubroutine cap_field
!***********************************************************************
    subroutine integrate_field(bz,aphi_mx)
!
!  Given a value of the magnetic field B , integrate it numerically to find the 
!  corresponding magnetic potential A, in cylindrical or spherical coordinates. 
!  In cylindrical coordinates, the field is 
!
!    r*Bz = d/dr(r*Aphi) 
!
!  whereas in spherical, it is 
!
!    r*Bphi = d/dr(r*Atheta)
!
!  Use the composite Simpson rule to integrate. 
!
!  26-feb-13/wlad: coded
!
      use FArrayManager, only: farray_use_global
      use Sub, only: gij,curl_mn,get_radial_distance
      use Mpicomm, only: mpibcast_real,mpisend_real,mpirecv_real
!
      real, dimension(mx) :: aphi_mx,bz
      real, dimension(nx) :: tmp
      real, dimension(0:nprocx-1) :: proc_store
      real :: out,in,psum
      integer :: partner, px
!
      call integrate(bz,tmp)
!
!  If the run is serial in x, we're done. Otherwise, take into account that 
!  the contribution of previous x-processors should be summed up. 
!
      if (nprocx/=1) then
!
!  Store the last value of the integral, which should be the starting point 
!  for the next x-processor.
!
         out=tmp(nx)
!
!  Prepare the communication in this yz row.
!
         do px=0,nprocx-1
            partner = px + nprocx*ipy + nprocxy*ipz
            if (iproc/=partner) then
               !Send to all processors in this row.
               call mpisend_real(out,partner,111)
               !Receive from all processors in the same row.
               call mpirecv_real(in,partner,111)
               proc_store(px)=in
            else !data is local
               proc_store(px) = out
            endif
         enddo
!
!  Sum the contributions of the x-processors in this y-row.
!
         psum=sum(proc_store(0:ipx-1))
         tmp=tmp+psum
      endif
!
      aphi_mx(l1:l2)=tmp
      aphi_mx(1:l1-1) =0.;aphi_mx(l2+1:mx)=0.
!
    endsubroutine integrate_field
!***********************************************************************
    subroutine integrate(bb,aa)
!
!  Integrate using Simpson's composite rule. 
!
      real, dimension(mx) :: bb
      real, dimension(0:nx) :: tmp
      real, dimension(nx) :: aa
      integer :: i,ig
!
      tmp(0)=0.
      do i=l1,l2
        ig=i-l1+1
        tmp(ig)=tmp(ig-1)+1./(6.*dx_1(i))*(&
                bb(i-3)   +   bb(i+3)   +&
             4*(bb(i-2)+bb(i)+bb(i+2))  +&
             2*(bb(i-1)   +   bb(i+1)))/3.
      enddo
      aa=tmp(1:nx)
!
    endsubroutine integrate
!***********************************************************************
    subroutine correct_pressure_gradient(f,ics2,temperature_power_law)
!
!  Correct for pressure gradient term in the centrifugal force.
!  For now, it only works for flat (isothermal) or power-law
!  sound speed profiles, because the temperature gradient is
!  constructed analytically.
!
!  21-aug-07/wlad : coded
!
      use FArrayManager
      use Sub,    only: get_radial_distance,grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: glnrho
      real, dimension (nx)   :: rr,rr_cyl,rr_sph
      real, dimension (nx)   :: cs2,fpres_thermal,gslnrho,gslnTT
      integer                :: ics2
      logical                :: lheader
      real :: temperature_power_law
!
      real, dimension(nx) :: xi
      real, dimension(0:6) :: c
!
      if (lroot) print*,'Correcting density gradient on the '//&
           'centrifugal force'
!
      do n=n1,n2 ; do m=m1,m2
        lheader=(lfirstpoint.and.lroot)
!
!  Get the density gradient
!
        call get_radial_distance(rr_sph,rr_cyl)
        call grad(f,ilnrho,glnrho)
        if (lcartesian_coords) then
          gslnrho=(glnrho(:,1)*x(l1:l2) + glnrho(:,2)*y(m))/rr_cyl
        else if (lcylindrical_coords) then
          gslnrho=glnrho(:,1)
        else if (lspherical_coords) then
          gslnrho=glnrho(:,1)
        endif
!
        if (lspherical_coords.or.lsphere_in_a_box) then
          rr=rr_sph
        else
          rr=rr_cyl
        endif
!
!  Get sound speed and calculate the temperature gradient
!
        if (llocal_iso.or.lenergy) then
          cs2=f(l1:l2,m,n,ics2)
          if (.not.lpolynomial_fit_cs2) then
             gslnTT=-temperature_power_law/((rr/r_ref)**2+rsmooth**2)*rr/r_ref**2
          else
             xi=x(l1:l2)
             c=coeff_cs2
             gslnTT = c(1)  + c(2) * 2*xi + c(3) * 3*xi**2 + &
                  c(4) * 4*xi**3 + c(5) * 5*xi**4 + c(6) * 6*xi**5
          endif   
       else
          cs2=cs20
          gslnTT=0.
       endif
!
!  Correct for cartesian or spherical
!
        fpres_thermal=(gslnrho+gslnTT)*cs2/gamma
!
        call correct_azimuthal_velocity(f,fpres_thermal,'pressure gradient')
!
      enddo;enddo
!
    endsubroutine correct_pressure_gradient
!***********************************************************************
    subroutine correct_selfgravity(f)
!
!  Correct for the fluid's self-gravity in the centrifugal force.
!  This routine incorporates selfgravity via the logspirals method
!  which is slightly different in that it uses the accelerations
!  directly instead of taking the gradient of the potential.
!
!  03-dec-07/wlad: coded
!  ??-???-??/vince: included logspirals
!
      use Sub,         only:get_radial_distance,grad
      use Poisson,     only:inverse_laplacian,get_acceleration
      use Boundcond,   only:update_ghosts
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      !real, dimension (nx,ny,nz) :: rho
      real, dimension (nx,ny,nz,3) :: acceleration
!
      real, dimension (mx,my,mz) :: selfpotential
      real, dimension (nx,3) :: gpotself
      real, dimension (nx) :: gspotself,rr_cyl,rr_sph
      real :: rhs_poisson_const
!
      if (lroot) print*,'Correcting for self-gravity on the '//&
           'centrifugal force'
      if (.not.lpoisson) then
        print*,"You want to correct for selfgravity but you "
        print*,"are using POISSON=nopoisson in src/Makefile.local. "
        print*,"Please use a poisson solver."
        call fatal_error("correct_selfgravity","")
      endif
!
!  Poisson constant is 4piG, this has to be consistent with the 
!  constant in poisson_init_pars
!
      if (lselfgravity_logspirals) then
        rhs_poisson_const=1.0
      else
        rhs_poisson_const=4*pi*gravitational_const
      endif
!
!  feed linear density into the poisson solver
!
      selfpotential(l1:l2,m1:m2,n1:n2)=&
           rhs_poisson_const*exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
!
      call inverse_laplacian(selfpotential(l1:l2,m1:m2,n1:n2))
!
      if (lselfgravity_logspirals) then
!
!  Get gravity directly
!
        call get_acceleration(acceleration)
      else
!
!  Define the selfpotential and update the boundaries to take the gradient
!
        f(l1:l2,m1:m2,n1:n2,ipotself)=selfpotential(l1:l2,m1:m2,n1:n2)
        call update_ghosts(f)
      endif
!
      do n=n1,n2 ; do m=m1,m2
!
        if (lselfgravity_logspirals) then
          gspotself = -acceleration(:,m-m1+1,n-n1+1,1)
        else
!
!  Get the potential gradient
!
          call get_radial_distance(rr_sph,rr_cyl)
          call grad(f,ipotself,gpotself)
!
!  correct the angular frequency phidot^2
!
          if (lcartesian_coords) then
            gspotself=(gpotself(:,1)*x(l1:l2) + gpotself(:,2)*y(m))/rr_cyl
          elseif (lcylindrical_coords) then
            gspotself=gpotself(:,1)
          elseif (lspherical_coords) then
            gspotself=gpotself(:,1)*sinth(m) + gpotself(:,2)*costh(m)
          endif
        endif
!
        call correct_azimuthal_velocity(f,gspotself,'self gravity')
!
      enddo;enddo
!
    endsubroutine correct_selfgravity
!***********************************************************************
    subroutine correct_lorentz_numerical(f)
!
      use Sub, only: get_radial_distance,&
                     gij,curl_mn,gij_etc,cross_mn,multsv_mn
      use Boundcond,only: update_ghosts
!
      real, dimension (mx,my,mz,mfarray) :: f      
      real, dimension (nx,3,3) :: aij,bij
      real, dimension (nx,3) :: aa,bb,jj,jxb,jxbr
      real, dimension (nx) :: rho1, fpres_magnetic
!
      if (lroot) print*,'Correcting Lorentz force on centrifugal balance'
!
      call update_ghosts(f)
!
      do n=n1,n2 ; do m=m1,m2
!
        aa=f(l1:l2,m,n,iax:iaz)         !aa
        call gij(f,iaa,aij,1)           !aij
        call curl_mn(aij,bb,aa)         !bb
        call gij_etc(f,iaa,aa,aij,bij)  !bij
        call curl_mn(bij,jj,bb)         !jj
        call cross_mn(jj,bb,jxb)        !jxb
        rho1=1./f(l1:l2,m,n,irho)       !1/rho
        call multsv_mn(rho1,jxb,jxbr)   !jxb/rho
        fpres_magnetic=-jxbr(:,1)
!
        call correct_azimuthal_velocity(f,fpres_magnetic,'magnetic pressure')  
!
      enddo; enddo
!
    endsubroutine correct_lorentz_numerical
!***********************************************************************
    subroutine correct_azimuthal_velocity(f,corr,label)
!
      use Sub, only: get_radial_distance
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nx), intent(in) :: corr
      real, dimension(nx) :: rr_sph, rr_cyl, rr, tmp1, tmp2
      character (len=*) :: label
!
      call get_radial_distance(rr_sph,rr_cyl)
!
      if (lcartesian_coords) then
        tmp1=(f(l1:l2,m,n,iux)**2+f(l1:l2,m,n,iuy)**2)/rr_cyl**2
        tmp2=tmp1 + corr/rr_cyl
      elseif (lcylindrical_coords) then
        tmp1=(f(l1:l2,m,n,iuy)/rr_cyl+OOcorot)**2
        tmp2=tmp1 + corr/rr_cyl
      elseif (lspherical_coords) then
        tmp1=(f(l1:l2,m,n,iuz)/rr_sph)**2
        tmp2=tmp1 + corr/rr_sph
      endif
!
!  Make sure the correction does not impede centrifugal equilibrium
!
      if (lcylindrical_coords.or.lcylinder_in_a_box) then
        rr=rr_cyl
      else
        rr=rr_sph
      endif
      call reality_check(tmp2,rr,label)
!
!  Correct the velocities
!
      if (lcartesian_coords) then
        f(l1:l2,m,n,iux)=-sqrt(tmp2)*y(  m  )
        f(l1:l2,m,n,iuy)= sqrt(tmp2)*x(l1:l2)
      elseif (lcylindrical_coords) then
        f(l1:l2,m,n,iuy)= (sqrt(tmp2)-OOcorot)*rr_cyl
      elseif (lspherical_coords) then
        f(l1:l2,m,n,iuz)= sqrt(tmp2)*rr_sph
      endif
!
    endsubroutine correct_azimuthal_velocity
!***********************************************************************
    subroutine reality_check(tmp,rr,label)
!
!  Catches unphysical negative values of phidot^2, i.e., impossibility 
!  of centrifugal equilibrium.
!
!  18-feb-13/wlad: moved to a separate subroutine because too many 
!                  subroutines called it. 
!
      use Messages, only: warning, fatal_error
!
      real, dimension(nx) :: tmp,rr 
      logical :: lheader
      integer :: i
      character (len=*) :: label
!      
      lheader=lroot.and.lfirstpoint
!
      do i=1,nx
        if (tmp(i)<0.) then
          if (rr(i) < r_int) then
            !it's inside the frozen zone, so
            !just set tmp2 to zero and emit a warning
             tmp(i)=0.
             if ((ip<=10).and.lheader) then
               print*,'reality_check: ',&
                    'cannot have centrifugal equilibrium in the inner ',&
                    'domain. ', label, ' is too strong at ',&
                    'x,y,z=',x(i+nghost),y(m),z(n)
               print*,'the angular frequency here is',tmp(i)
               call warning("","")
             endif
          else
            print*,'reality_check: ',&
                 'cannot have centrifugal equilibrium in the inner ',&
                 'domain. ', label, ' is too strong at ',&
                 'x,y,z=',x(i+nghost),y(m),z(n)
            print*,'the angular frequency here is',tmp(i)
            call fatal_error("","")
          endif
        endif
      enddo
!
    endsubroutine reality_check
!***********************************************************************
    subroutine lowk_noise_gaussian_rprof(f)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  07-may-09/wlad: coded
!
      use General, only: random_number_wrapper
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx,ny,nz) :: lump_of_sines
      real :: Lx,Ly,Lz,d0,phase,xi,yi,zi,nw1,fac
      real :: fmeantmp_rho2,fmeantmp_rho
      real :: fmean_rho2,fmean_rho
      real :: unnormalized_rho_rms
      real :: normalization_factor
      integer :: i,mm,nn,ll,mm1,ll1,nn1,itmp,irho
!
      Lx=Lxyz(1) ; Ly=Lxyz(2) ; Lz=Lxyz(3)
!
      d0=0.2*Lx
!
      if (lroot) print*,'domain size=',Lx,Ly,Lz
!
!  save the original linear density 
!  on the free spot of shock, which isn't set yet
!
      if (lshock) then
        itmp=ishock
        f(l1:l2,m1:m2,n1:n2,itmp)=f(l1:l2,m1:m2,n1:n2,ilnrho)
      else 
        if (lroot) then
          print*,'This is a baroclinic run. You expect '
          print*,'many vortices to form in this run, do you not?'
          print*,'These beasts usually get supersonic, which is '
          print*,'actually what prevents them from growing too much,'
          print*,'as then they shock and dissipate. Yet, you are'
          print*,'NOT using shock viscosity. I recommend you to stop'
          print*,'and switch SHOCK=shock_highorder in src/Makefile.local'
        endif
        call fatal_error("","")     
      endif
!
!  Work with LINEAR density. As this is start time, there is no confusion
!  linear-log. All initial conditions are to be coded in log. We will 
!  construct the density here as linear, then pass to log
!
      irho=ilnrho
!
!  Have to set it to something first, otherwise the sines 
!  can lead to negative densities. The density is saved in the 
!  slot reserved to shock. 
!
      f(l1:l2,m1:m2,n1:n2,irho)=1.
      lump_of_sines=0.
!      
      do ll=-xmodes,xmodes ; do mm=0,ymodes ; do nn=-zmodes,zmodes
!
        if (lroot) call random_number_wrapper(phase)
        call mpibcast_real(phase)
!
        do i=1,nx ; do m=1,ny ; do n=1,nz
          ll1=i+l1-1 ; xi=x(ll1)
          mm1=m+m1-1 ; yi=y(mm1)
          nn1=n+n1-1 ; zi=z(nn1)
!
! Exclude the border
!
          if ((xi.gt.xyz0(1)+2*rborder_int).and.(xi.lt.xyz1(1)-2*rborder_ext)) then
            lump_of_sines(i,m,n)=lump_of_sines(i,m,n) + &
                 sin(2*pi*(ll*xi/Lxn + mm*yi/Ly + nn*zi/Lz+ phase))
          endif
!
        enddo;enddo;enddo !end grid loop
      enddo;enddo;enddo !end modes loop
!
!  Now construct the density and fill in the normalization 
!  constants needed to get a rms equal to the input rho_rms.
!
      fmeantmp_rho2=0.
      fmeantmp_rho=0.
      nw1=1./(nxgrid*nygrid*nzgrid*1.0)
!      
      do n=1,nz
        nn1=n+n1-1
        do m=1,ny
          mm1=m+m1-1
          do i=1,nx
            ll1=i+l1-1 ; xi=x(ll1)
!
            if (lgaussian_distributed_noise) then
              fac=exp(-(.5*(xi-xmid)/d0)**2)
            else
              fac=1
            endif
!
            f(ll1,mm1,nn1,irho)=f(ll1,mm1,nn1,irho) + &
                 lump_of_sines(i,m,n)*fac
          enddo
          fmeantmp_rho2=fmeantmp_rho2+nw1*sum(f(l1:l2,mm1,nn1,irho)**2)
          fmeantmp_rho =fmeantmp_rho +nw1*sum(f(l1:l2,mm1,nn1,irho))
        enddo
      enddo !end grid loop
!
!  Sum the normalization constants over processors, and perform
!  the normalization.
!
      call mpireduce_sum(fmeantmp_rho2,fmean_rho2)
      call mpireduce_sum(fmeantmp_rho,fmean_rho)
      call mpibcast_real(fmean_rho2)
      call mpibcast_real(fmean_rho)
!
      unnormalized_rho_rms=sqrt(fmean_rho2-fmean_rho**2)
      normalization_factor=rho_rms/unnormalized_rho_rms
!
!  Assumes rho0=1.
!
      f(l1:l2,m1:m2,n1:n2,irho)=1.+&
          normalization_factor*(f(l1:l2,m1:m2,n1:n2,irho)-1.)
!
      if (lroot) then
        print*,'max density (linear): ',maxval(f(l1:l2,m1:m2,n1:n2,irho))
        print*,'min density (linear): ',minval(f(l1:l2,m1:m2,n1:n2,irho))
        print*,'rms density (linear): ',&
            unnormalized_rho_rms*normalization_factor
      endif
!
!  convert to log before finishing. 
!
      f(l1:l2,m1:m2,n1:n2,ilnrho)=&
          f(l1:l2,m1:m2,n1:n2,itmp)+log(f(l1:l2,m1:m2,n1:n2,irho))
!
!  keep the original stratification in the ishock slot since it 
!  will be needed when setting the entropy
!
    endsubroutine lowk_noise_gaussian_rprof
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
