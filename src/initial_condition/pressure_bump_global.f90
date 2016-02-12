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
  real :: g0=1.,density_power_law=0.,temperature_power_law=1.
  real :: radial_percent_smooth=10.0,rshift=0.0
  real :: gravitational_const=0.
!
  real, dimension(0:6) :: coeff_cs2=(/0.,0.,0.,0.,0.,0.,0./)
!
  real :: r0_pot=0.,qgshear=1.5
  integer :: n_pot=10
!
! Correct extra forces in the centrifugal balance condition
!
  logical :: lcorrect_pressuregradient=.true.
  logical :: lpolynomial_fit_cs2=.false.
  logical :: ladd_noise_propto_cs=.false.
  real :: ampluu_cs_factor=1d-3
!
!  Quantities for the pressure bump
!
  real :: bump_radius = 1.,bump_ampl = 0.4, bump_width = 0.1
!
  namelist /initial_condition_pars/ g0,density_power_law,&
       temperature_power_law,radial_percent_smooth,rshift,&
       gravitational_const,r0_pot,qgshear,n_pot,&
       lcorrect_pressuregradient,lpolynomial_fit_cs2,&
       ladd_noise_propto_cs,ampluu_cs_factor,&
       bump_radius,bump_ampl,bump_width
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
        elseif (lparticles_nbody) then
!
! Nbody gravity with a dominating but dynamical central body
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
        endif
!
        if (coord_system=='cartesian') then
          f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) - y(  m  )*OO
          f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + x(l1:l2)*OO
          f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + 0.
        elseif (coord_system=='cylindric') then
          f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) + 0.
          f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + OO*rr_cyl
          f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + 0.
        elseif (coord_system=='spherical') then
          f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) + 0.
          f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + 0.
          f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + OO*rr_sph
        endif
!
      enddo
      enddo
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
      use General, only: random_number_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i1,i2
!
      real, dimension (nx) :: r,p,tmp,ampl
      integer :: i
!
      intent(in)    :: ampl,i1,i2
      intent(inout) :: f
!
!  set gaussian random noise vector
!
      do i=i1,i2
        if (lroot.and.m==1.and.n==1) print*,'gaunoise_vect: variable i=',i
        if (modulo(i-i1,2)==0) then
          call random_number_wrapper(r)
          call random_number_wrapper(p)
          tmp=sqrt(-2*log(r))*sin(2*pi*p)
        else
          tmp=sqrt(-2*log(r))*cos(2*pi*p)
        endif
        f(l1:l2,m,n,i)=f(l1:l2,m,n,i)+ampl*tmp
      enddo
!
    endsubroutine gaunoise_vect
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
      real, dimension (mx)   :: rr_sph,rr,rr_cyl,lnrhomid,rho_tmp
      real                   :: lat
      integer, pointer       :: iglobal_cs2,iglobal_glnTT
      integer                :: ics2
      logical                :: lheader,lpresent_zed
!
      if (lroot) print*,&
           'initial_condition_lnrho: locally isothermal approximation'
      if (lroot) print*,'Radial stratification with power law=',&
           density_power_law
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
          elseif (ltemperature) then
            if (ltemperature_nolog) then 
              ics2=iTT
            else
              ics2=ilnTT
            endif
          elseif (lentropy) then
            ics2=iss
          endif
          f(:,m,n,ics2)=cs2
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
        rho_tmp = (1 + bump_ampl*exp(-(rr_cyl - bump_radius)**2/(2*bump_width**2))) * rr_cyl**(-density_power_law)
        lnrhomid = log(rho_tmp)
        !lnrhomid=log(rho0) -.5*density_power_law*log((rr/r_ref)**2+rsmooth**2)
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
            elseif (lparticles_nbody) then
              !call get_totalmass(g0) 
              tmp1=-g0/rr_sph**2
            else
              print*,"both gravity and particles_nbody are switched off"
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
            elseif (lparticles_nbody) then
              !call get_totalmass(g0)
              tmp1=-g0/rr_sph 
              tmp2=-g0/rr_cyl
            else
              print*,"both gravity and particles_nbody are switched off"
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
        f(:,m,n,ilnrho) = f(:,m,n,ilnrho) + lnrhomid + strat
      enddo
      enddo
!
!  Correct the velocities by this pressure gradient
!
      if (lcorrect_pressuregradient) &
           call correct_pressure_gradient(f,ics2,temperature_power_law)
!
!  Set the thermodynamical variable
!
      if (llocal_iso) then
        call set_thermodynamical_quantities&
             (f,temperature_power_law,ics2,iglobal_cs2,iglobal_glnTT)
      else if (lenergy) then 
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
        call correct_azimuthal_velocity(f,fpres_thermal)
!
      enddo;enddo
!
    endsubroutine correct_pressure_gradient
!***********************************************************************
    subroutine correct_azimuthal_velocity(f,corr)
!
      use Sub, only: get_radial_distance
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nx), intent(in) :: corr
      real, dimension(nx) :: rr_sph, rr_cyl, rr, tmp1, tmp2
!
      call get_radial_distance(rr_sph,rr_cyl)
!
      if (lcartesian_coords) then
        tmp1=(f(l1:l2,m,n,iux)**2+f(l1:l2,m,n,iuy)**2)/rr_cyl**2
        tmp2=tmp1 + corr/rr_cyl
      elseif (lcylindrical_coords) then
        tmp1=(f(l1:l2,m,n,iuy)/rr_cyl)**2
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
      call reality_check(tmp2,rr)
!
!  Correct the velocities
!
      if (lcartesian_coords) then
        f(l1:l2,m,n,iux)=-sqrt(tmp2)*y(  m  )
        f(l1:l2,m,n,iuy)= sqrt(tmp2)*x(l1:l2)
      elseif (lcylindrical_coords) then
        f(l1:l2,m,n,iuy)= sqrt(tmp2)*rr_cyl
      elseif (lspherical_coords) then
        f(l1:l2,m,n,iuz)= sqrt(tmp2)*rr_sph
      endif
!
    endsubroutine correct_azimuthal_velocity
!***********************************************************************
    subroutine reality_check(tmp,rr)
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
!      
      lheader=lroot.and.lfirstpoint
!
      do i=1,nx
        if (tmp(i)<0.) then
          if (rr(i) < r_int) then
            !it's inside the frozen zone, so
            !just set tmp2 to zero and emit a warning
             tmp(i)=0.
             if ((ip<=10).and.lheader) &
                  call warning('reality_check','Cannot '//&
                  'have centrifugal equilibrium in the inner '//&
                  'domain. The pressure gradient is too steep.')
          else
            print*,'reality_check: ',&
                 'cannot have centrifugal equilibrium in the inner ',&
                 'domain. The pressure gradient is too steep at ',&
                 'x,y,z=',x(i+nghost),y(m),z(n)
            print*,'the angular frequency here is',tmp(i)
            call fatal_error("","")
          endif
        endif
      enddo
!
    endsubroutine reality_check
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
