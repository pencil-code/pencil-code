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
!!  integer :: dummy
!
  real :: g0=1.,density_power_law=0.,temperature_power_law=1., plasma_beta=25
  logical :: lexponential_smooth=.false.
  real :: radial_percent_smooth=10.0,rshift=0.0
  logical :: lcorrect_selfgravity=.false.
  real :: gravitational_const=0.
  real :: magnetic_power_law=impossible
!
! For the magnetic field
!
  logical :: ladd_field=.true.
  character (len=labellen) :: initcond_aa=''
  real, dimension(3) :: B_ext=(/0.0,0.0,0.0/)
  real :: rmode_mag=4.,zmode_mag=16.,amplbb=1.
  real :: rm_int=0.0,rm_ext=impossible
  real :: Bz_const=8d-3
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
  namelist /initial_condition_pars/ g0,density_power_law,&
       temperature_power_law,lexponential_smooth,&
       radial_percent_smooth,rshift,lcorrect_selfgravity,&
       gravitational_const,xmodes,ymodes,zmodes,rho_rms,&
       llowk_noise,xmid,lgaussian_distributed_noise,rborder_int,&
       rborder_ext,plasma_beta,ladd_field,initcond_aa,B_ext,&
       zmode_mag,rmode_mag,rm_int,rm_ext,amplbb,Bz_const, &
       r0_pot,qgshear,n_pot,magnetic_power_law
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
      Lxn=Lx-2*(rborder_int+rborder_ext)
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
      real, dimension (mx)   :: rr_sph,rr,rr_cyl,lnrhomid
      real                   :: rmid,lat
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
        if (lsphere_in_a_box.or.lspherical_coords) then
          rr=rr_sph
        elseif (lcylinder_in_a_box.or.lcylindrical_coords) then
          rr=rr_cyl
        else
          call fatal_error("initial_condition_lnrho",&
               "no valid coordinate system")
        endif
!
        call power_law(cs20,rr,temperature_power_law,cs2,r_ref)
!
!  Store cs2 in one of the free slots of the f-array
!
        if (llocal_iso) then
          nullify(iglobal_cs2)
          call farray_use_global('cs2',iglobal_cs2)
          ics2=iglobal_cs2
        elseif (ltemperature) then
          ics2=ilnTT
        elseif (lentropy) then
          ics2=iss
        endif
        f(:,m,n,ics2)=cs2
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
        if (lsphere_in_a_box.or.lspherical_coords) then
          rr=rr_sph
        elseif (lcylinder_in_a_box.or.lcylindrical_coords) then
          rr=rr_cyl
        endif
!
        if (lexponential_smooth) then
          !radial_percent_smooth = percentage of the grid
          !that the smoothing is applied
          rmid=rshift+(xyz1(1)-xyz0(1))/radial_percent_smooth
          lnrhomid=log(rho0) &
               + density_power_law*log((1-exp( -((rr-rshift)/rmid)**2 ))/rr)
        else
          lnrhomid=log(rho0) -.5*density_power_law*log((rr/r_ref)**2+rsmooth**2)
        endif
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
          cs2=f(:,m,n,ics2)
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
      call correct_pressure_gradient(f,ics2,temperature_power_law)
!
!  Correct the velocities for self-gravity
!
      call correct_for_selfgravity(f)
!
!  Set noise in low wavelengths only
!
      if (llowk_noise) call lowk_noise_gaussian_rprof(f)
!
!  Set the thermodynamical variable
!
      if (llocal_iso) then
        call set_thermodynamical_quantities&
             (f,temperature_power_law,ics2,iglobal_cs2,iglobal_glnTT)
      else
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
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  07-may-09/wlad: coded
!
      use EquationOfState, only: cs20
      use Sub, only: step_scalar, power_law
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: amplbb,pblaw
      real, dimension(nx) :: Aphi, rr, rrcyl
      real, dimension(mx) :: bz,aphi_mx,Breal,pp
      character (len=labellen) :: intlabel
!
      real :: const,k1,kr
      integer :: i
!
      select case (initcond_aa)
!
!  Add force-free field of constant plasma beta. 
!
        case ('plasma-beta')
          if (ladd_field) then 
            amplbb = 2*cs20*rho0/plasma_beta
!
            do m=m1,m2 
              do n=n1,n2
                Aphi = -amplbb * log(x(l1:l2))
!
                if (lspherical_coords) then 
                  f(l1:l2,m,n,iaz) = Aphi
                elseif (lcylindrical_coords) then
                  f(l1:l2,m,n,iay) = Aphi
                else
                  print*,'case plasma-beta not coded for Cartesian '
                  print*,'coordinates'
                  call fatal_error("initial_condition_aa","")
                endif
              enddo
            enddo
          endif
!
!--------------------------------------------------------------
!
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
            call correct_lorentz_force(f,amplbb,pblaw)
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
              Aphi =  amplbb/kr * sin(kr*(rrcyl-r_int)) + &
                   amplbb/(kr**2*rrcyl)*cos(kr*(rrcyl-r_int))
            else
              if (lroot .and. m==m1 .and. n==n1) &
                   print*,'Softened magnetic field in the center'
              if (rmode_mag < 5) call fatal_error("initial_condition_aa",&
                   "put more wavelengths in the field")
              kr = 2*pi*rmode_mag/r_ext
              k1 = 1. !not tested for other values
              const=amplbb*exp(1.)*k1/cos(kr/k1)
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
           f(l1:l2,m,n,iaz) =  -amplbb/kr*sin(kr*(rrcyl-r_int))
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
             bz(i)=amplbb/rcyl_mn(i) * sin(kr*(rcyl_mn(i)-rm_int))
           else
             bz(i)=0.
           endif
         enddo
!
         intlabel=trim('tridag')
         call integrate_field(bz,aphi_mx,intlabel)
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
        !if (lspherical_coords) call correct_lorentz_numerical(f)
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine set_field(f,Breal)
!
      use Sub, only: step_scalar
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: width
      real, dimension(mx) :: bz,aphi_mx,Breal
      integer :: i
!
      do i=1,mx
!
        width=5./dx_1(i)
!
        bz(i) = Breal(i) * &
             (step_scalar(x(i),rm_int+width,width)-&
              step_scalar(x(i),rm_ext-width,width))
      enddo
!
      call integrate_field_parallel(bz*x,aphi_mx)
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
    subroutine correct_lorentz_numerical(f)
!
      use Sub, only: get_radial_distance,&
                     gij,curl_mn,gij_etc,cross_mn,multsv_mn
!
      real, dimension (mx,my,mz,mfarray) :: f      
      real, dimension (nx,3,3) :: aij,bij
      real, dimension (nx,3) :: aa,bb,jj,jxb,jxbr
      real, dimension (nx) :: rr_sph,rr_cyl,rho1
      real, dimension (nx) :: fpres_magnetic
!
      if (.not.lspherical_coords) then
        call fatal_error("correct_lorentz_numerical",&
            "only implemented for spherical coordinates")
        if ((B_ext(1)/=0).or.(B_ext(3)/=0)) then
          call fatal_error("correct_lorentz_numerical",&
              "only implemented for polar fields")
        endif
      endif
!
      do m=m1,m2
        do n=n1,n2
!
          if (lmagnetic) then
            aa=f(l1:l2,m,n,iax:iaz)         !aa
            call gij(f,iaa,aij,1)           !aij
            call curl_mn(aij,bb,aa)         !bb
            call gij_etc(f,iaa,aa,aij,bij)  !bij
            call curl_mn(bij,jj,bb)         !jj
            call cross_mn(jj,bb,jxb)        !jxb
            rho1=1./f(l1:l2,m,n,irho)       !1/rho
            call multsv_mn(rho1,jxb,jxbr)   !jxb/rho
            fpres_magnetic=-jxbr(:,2)
          else
            fpres_magnetic=0.
          endif
!         
          call get_radial_distance(rr_sph,rr_cyl)
!     
! Not valid for all signs of the pressure gradient
!
          f(l1:l2,m,n,iuz)=f(l1:l2,m,n,iuz)+&
              sqrt(rr_sph*fpres_magnetic/cotth(m))
!
        enddo
      enddo
!
    endsubroutine correct_lorentz_numerical
!***********************************************************************
    subroutine integrate_field_parallel(bz,aphi_mx)
!
      use FArrayManager, only: farray_use_global
      use Sub, only: gij,curl_mn,get_radial_distance
      use EquationOfState, only: cs20,rho0
      use Mpicomm, only: mpibcast_real,mpisend_real,mpirecv_real
!
      real, dimension(mx) :: aphi_mx,bz
      real, dimension(0:nx) :: tmp2
      integer :: i,ig,iprocx,iprocy,iserial_xy
      real, dimension(nprocx*nprocy) :: procsum
      real :: procsum_loc,tmpy,psum

!
      iserial_xy=ipx+nprocx*ipy+1
      tmp2(0)=0.
      do i=l1,l2
        ig=i-l1+1
        tmp2(ig)=tmp2(ig-1)+1./(6.*dx_1(i))*(&
                bz(i-3)   +   bz(i+3)   +&
             4*(bz(i-2)+bz(i)+bz(i+2))  +&
             2*(bz(i-1)   +   bz(i+1)))/3.
      enddo
      procsum_loc=tmp2(nx)
!
      if (ip<=9) print*,'send',ipx+nprocx*ipy+1,procsum_loc
!
      if (lroot) then 
        procsum(1)=procsum_loc
        do iprocx=0,nprocx-1; do iprocy=0,nprocy-1
          iserial_xy=iprocx+nprocx*iprocy+1
          if (iserial_xy/=1) then
            call mpirecv_real(tmpy,1,iserial_xy-1,111)
            procsum(iserial_xy)=tmpy
          endif
          if (ip<=9) print*,'recv',iserial_xy,procsum(iserial_xy)
        enddo; enddo
      else
        call mpisend_real(procsum_loc,1,0,111)
      endif
!
      call mpibcast_real(procsum,nprocx*nprocy)
!
      do n=n1,n2
!
        if (nprocx==1) then 
          tmp2(0)=0.
          do i=l1,l2
            ig=i-l1+1
            tmp2(ig)=tmp2(ig-1)+1./(6.*dx_1(i))*(&
                    bz(i-3)   +   bz(i+3)   +&
                 4*(bz(i-2)+bz(i)+bz(i+2))  +&
                 2*(bz(i-1)   +   bz(i+1)))/3.
          enddo
        else
          if (lfirst_proc_x) then 
            tmp2(0)=0.
            do i=l1,l2
              ig=i-l1+1
              tmp2(ig)=tmp2(ig-1)+1./(6.*dx_1(i))*(&
                   bz(i-3)   +   bz(i+3)   +&
                   4*(bz(i-2)+bz(i)+bz(i+2))  +&
                   2*(bz(i-1)   +   bz(i+1)))/3.
            enddo
          else 
            psum=0.
            do iprocx=0,ipx-1
              iserial_xy=iprocx+nprocx*ipy+1
              psum=psum+procsum(iserial_xy)
            enddo
            tmp2(0)=psum
            do i=l1,l2
              ig=i-l1+1
              tmp2(ig)=tmp2(ig-1)+1./(6.*dx_1(i))*(&
                   bz(i-3)   +   bz(i+3)   +&
                   4*(bz(i-2)+bz(i)+bz(i+2))  +&
                   2*(bz(i-1)   +   bz(i+1)))/3.
            enddo
          endif
        endif
        aphi_mx(l1:l2)=tmp2(1:nx)
        aphi_mx(1:l1-1) =0.;aphi_mx(l2+1:mx)=0.
      enddo
! 
    endsubroutine integrate_field_parallel
!***********************************************************************
    subroutine integrate_field(bz,aphi_mx,integration_method)
!
      use General, only: tridag
!
      real, dimension(mx) :: a_tri,b_tri,c_tri,bz,rhs,aphi_mx
      character (len=labellen) :: integration_method
      integer :: i
!
      select case (integration_method)
!
      case ('tridag') 

        a_tri=-1. ; b_tri=2/(x*dx_1) ; c_tri=1. ; rhs=bz*2/dx_1
!
        a_tri(1) =0.;c_tri(1 )=0.
        a_tri(mx)=0.;c_tri(mx)=0.
!
        call tridag(a_tri,b_tri,c_tri,rhs,aphi_mx)
!
      case ('simpson')  
!
        aphi_mx(1:l1-1)=0.
        do i=l1,l2
          aphi_mx(i) = aphi_mx(i-1)+1./(6.*dx_1(i))*(&
                  bz(i-3)   +   bz(i+3)   +&
               4*(bz(i-2)+bz(i)+bz(i+2))  +&
               2*(bz(i-1)   +   bz(i+1)))/3.
          !aphi_mx(i) = aphi_mx(i-1)+dx/2.*(bz(i-1)+4*bz(i)+bz(i+1))/2.
        enddo
        !aphi_mx(1:l1-1)=0.; aphi_mx(l2+1:mx)=0.
!
      case default
!
!  Catch unknown values.
!
        call fatal_error('integrate_field', &
             'integration_method "' // trim(integration_method) // '"'//&
             ' not recognised')
!
      endselect
!
    endsubroutine integrate_field
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
      real, dimension (nx)   :: cs2,tmp1,tmp2,corr,gslnrho,gslnTT
      integer                :: i,ics2
      logical                :: lheader
      real :: temperature_power_law
!
      if (lroot) print*,'Correcting density gradient on the '//&
           'centrifugal force'
!
      do m=m1,m2
      do n=n1,n2
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
!  Get sound speed and calculate the temperature gradient
!
        cs2=f(l1:l2,m,n,ics2);rr=rr_cyl
        if (lspherical_coords.or.lsphere_in_a_box) rr=rr_sph
        gslnTT=-temperature_power_law/((rr/r_ref)**2+rsmooth**2)*rr/r_ref**2
!
!  Correct for cartesian or spherical
!
        corr=(gslnrho+gslnTT)*cs2
        if (lenergy) corr=corr/gamma
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
        do i=1,nx
          if (tmp2(i)<0.) then
            if (rr(i) < r_int) then
              !it's inside the frozen zone, so
              !just set tmp2 to zero and emit a warning
              tmp2(i)=0.
              if ((ip<=10).and.lheader) &
                   call warning('correct_density_gradient','Cannot '//&
                   'have centrifugal equilibrium in the inner '//&
                   'domain. The pressure gradient is too steep.')
            else
              print*,'correct_density_gradient: ',&
                   'cannot have centrifugal equilibrium in the inner ',&
                   'domain. The pressure gradient is too steep at ',&
                   'x,y,z=',x(i+nghost),y(m),z(n)
              print*,'the angular frequency here is',tmp2(i)
              call fatal_error("","")
            endif
          endif
        enddo
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
      enddo
      enddo
!
    endsubroutine correct_pressure_gradient
!***********************************************************************
    subroutine exponential_fall(f)
!
!  Exponentially falling radial density profile.
!
!  21-aug-07/wlad: coded
!
      use Sub, only: get_radial_distance
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: rr_sph,rr_cyl,arg
      real :: rmid,fac
!
      if (lroot) print*,'setting exponential falling '//&
           'density with e-fold=',r_ref
!
      fac=pi/2
      do m=1,my
      do n=1,mz
        call get_radial_distance(rr_sph,rr_cyl)
        f(:,m,n,ilnrho) = f(:,m,n,ilnrho) + lnrho0 - rr_cyl/r_ref
        if (lexponential_smooth) then
          rmid=rshift+(xyz1(1)-xyz0(1))/radial_percent_smooth
          arg=(rr_cyl-rshift)/rmid
          f(:,m,n,ilnrho) = f(:,m,n,ilnrho)*&
               .5*(1+atan(arg)/fac)
        endif
      enddo
      enddo
!
!  Add self-gravity's contribution to the centrifugal force
!
      call correct_for_selfgravity(f)
!
    endsubroutine exponential_fall
!***********************************************************************
    subroutine correct_for_selfgravity(f)
!
!  Correct for the fluid's self-gravity in the
!  centrifugal force
!
!  03-dec-07/wlad: coded
!
      use Sub,         only:get_radial_distance,grad
      use Poisson,     only:inverse_laplacian
      use Boundcond,   only:update_ghosts
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (mx,my,mz) :: selfpotential
      real, dimension (nx,3) :: gpotself
      real, dimension (nx) :: tmp1,tmp2
      real, dimension (nx) :: gspotself,rr_cyl,rr_sph
      real :: rhs_poisson_const
      logical :: lheader
      integer :: i
!
!  Do nothing if self-gravity is not called
!
      if (lcorrect_selfgravity) then
!
        if (lroot) print*,'Correcting for self-gravity on the '//&
             'centrifugal force'
        if (.not.lpoisson) then 
          print*,"You want to correct for selfgravity but you "
          print*,"are using POISSON=nopoisson in src/Makefile.local. "
          print*,"Please use a poisson solver."
          call fatal_error("correct_for_selfgravity","")
        endif
!
!  Poisson constant is 4piG, this has to be consistent with the 
!  constant in poisson_init_pars
!
        rhs_poisson_const=4*pi*gravitational_const
!
!  feed linear density into the poisson solver
!
        selfpotential(l1:l2,m1:m2,n1:n2)=&
             rhs_poisson_const*exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
        call inverse_laplacian(f,&
             selfpotential(l1:l2,m1:m2,n1:n2))
!
        f(l1:l2,m1:m2,n1:n2,ipotself)=selfpotential(l1:l2,m1:m2,n1:n2)
        call update_ghosts(f)
!
!  update the boundaries for the self-potential
!
        do n=n1,n2
        do m=m1,m2
!
          lheader=(lfirstpoint.and.lroot)
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
            tmp1=(f(l1:l2,m,n,iux)**2+f(l1:l2,m,n,iuy)**2)/rr_cyl**2
            tmp2=tmp1+gspotself/rr_cyl
          elseif (lcylindrical_coords) then
            gspotself=gpotself(:,1)
            tmp1=(f(l1:l2,m,n,iuy)/rr_cyl)**2
            tmp2=tmp1+gspotself/rr_cyl
          elseif (lspherical_coords) then
            gspotself=gpotself(:,1)*sinth(m) + gpotself(:,2)*costh(m)
            tmp1=(f(l1:l2,m,n,iuz)/(rr_sph*sinth(m)))**2
            tmp2=tmp1 + gspotself/(rr_sph*sinth(m)**2)
          endif
!
!  Catch negative values of phidot^2
!
          do i=1,nx
            if (tmp2(i)<0.) then
              if (rr_cyl(i) < r_int) then
                !it's inside the frozen zone, so
                !just set tmp2 to zero and emit a warning
                tmp2(i)=0.
                if ((ip<=10).and.lheader) &
                     call warning('correct_for_selfgravity','Cannot '//&
                     'have centrifugal equilibrium in the inner '//&
                     'domain. Just warning...')
              else
                print*,'correct_for_selfgravity: ',&
                     'cannot have centrifugal equilibrium in the inner ',&
                     'domain. The offending point is ',&
                     'x,y,z=',x(i+nghost),y(m),z(n)
                print*,'the angular frequency here is ',tmp2(i)
                call fatal_error("correct_for_selfgravity","")
              endif
            endif
          enddo
!
!  Correct the velocities
!
          if (lcartesian_coords) then
            f(l1:l2,m,n,iux)=-sqrt(tmp2)*y(  m  )
            f(l1:l2,m,n,iuy)= sqrt(tmp2)*x(l1:l2)
          elseif (lcylindrical_coords) then
            f(l1:l2,m,n,iuy)= sqrt(tmp2)*rr_cyl
          elseif (lspherical_coords) then
            f(l1:l2,m,n,iuz)= sqrt(tmp2)*rr_sph*sinth(m)
          endif
        enddo
        enddo
      endif ! if (lcorrect_selfgravity)
!
    endsubroutine correct_for_selfgravity
!***********************************************************************
    subroutine correct_lorentz_force(f,const,pblaw)
!
!  Correct for the magnetic term in the centrifugal force. The
!  pressure gradient was already corrected in the density and temperature
!  modules
!
!  13-nov-08/wlad : coded
!
      use Sub,      only: get_radial_distance
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx)   :: rr_cyl,rr_sph,Btheta
      real, dimension (mx)   :: tmp,uu2,va2,va2r
      real, dimension (mx)   :: rho1,rho1_jxb
      real                   :: const,pblaw
      logical                :: lheader
      integer                :: i
!
      if (lroot) print*,'Correcting magnetic terms on the '//&
           'centrifugal force'
!
      if (.not.lspherical_coords) then
        call fatal_error("correct_lorentz_force",&
            "only implemented for spherical coordinates")
        if ((B_ext(1)/=0).or.(B_ext(3)/=0)) then
          call fatal_error("correct_lorentz_force",&
              "only implemented for polar fields")
        endif
      endif
!
      do m=1,my
        do n=1,mz
!
          call get_radial_distance(rr_sph,rr_cyl)
!
          lheader=((m==1).and.(n==1).and.lroot)
!
          !this field also has a magnetic pressure gradient
          Btheta=const*rr_sph**pblaw/sin(y(m))
!
          rho1=1./f(:,m,n,ilnrho)
!
          uu2=f(:,m,n,iuz)**2
          va2=rho1*(Btheta+B_ext(2))**2
!
          rho1_jxb=rho1
!
!  set rhomin_jxb>0 in order to limit the jxb term at very low densities.
!
          !if (rhomin_jxb>0) rho1_jxb=min(rho1_jxb,1/rhomin_jxb)
!
!  set va2max_jxb>0 in order to limit the jxb term at very high Alfven speeds.
!  set va2power_jxb to an integer value in order to specify the power
!  of the limiting term,
!
          !if (va2max_jxb>0) then
          !  rho1_jxb = rho1_jxb &
          !      * (1+(va2/va2max_jxb)**va2power_jxb)**(-1.0/va2power_jxb)
          !endif
          va2r=rho1_jxb*va2
!
!  The second term is the magnetic pressure gradient
!
 
          tmp=uu2+va2r*(1+2*pblaw/rr_sph)
!
!  The polar pressure should be
!  -2*cot(theta)/rr * va2r, but I will ignore it
!  for now. It feedbacks on the density,
!  so the initial condition for density and field
!  should be solved iteratively. But as uu>>va, I
!  will just let the system relax to equilibrium in
!  runtime.
!
!  Make sure the correction does not impede centrifugal equilibrium
!
          do i=1,nx
            if (tmp(i)<0.) then
              if (rr_sph(i) < r_int) then
                !it's inside the frozen zone, so
                !just set tmp to zero and emit a warning
                tmp(i)=0.
                if ((ip<=10).and.lheader) &
                    call warning('correct_lorentz_force','Cannot '//&
                    'have centrifugal equilibrium in the inner '//&
                    'domain. The lorentz force is too strong.')
              else
                print*,'correct_lorentz_force: ',&
                    'cannot have centrifugal equilibrium in the inner ',&
                    'domain. The lorentz force is too strong at ',&
                    'x,y,z=',x(i+nghost),y(m),z(n)
                print*,'the angular frequency here is',tmp(i)
                call fatal_error("","")
              endif
            endif
          enddo
!
!  Correct the velocities
!
          f(:,m,n,iuz)= sqrt(tmp)
!
        enddo
      enddo
!
    endsubroutine correct_lorentz_force
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
                                 cs20,cs2bot,cs2top,lnrho0,rho0
      use Sub,             only: power_law,get_radial_distance
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: rr,rr_sph,rr_cyl,cs2,lnrho
      real, dimension (nx) :: gslnTT
      real :: cp1,temperature_power_law
      integer, pointer, optional :: iglobal_cs2,iglobal_glnTT
      integer :: ics2
!
      intent(in)  :: temperature_power_law
      intent(out) :: f
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
          f(l1:l2,m,n,ilnTT)=log(cs2*cp1/gamma_m1)
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
    endsubroutine set_thermodynamical_quantities
!***********************************************************************
    subroutine lowk_noise_gaussian_rprof(f)
!
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
        call mpibcast_real(phase,1)
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
      call mpibcast_real(fmean_rho2,1)
      call mpibcast_real(fmean_rho,1)
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
    subroutine read_initial_condition_pars(unit,iostat)
!
!  07-may-09/wlad: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=initial_condition_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=initial_condition_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!     
      integer, intent(in) :: unit
!
      write(unit,NML=initial_condition_pars)
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
