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
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
  use EquationOfState
!
  implicit none
!
  include '../initial_condition.h'
!
!!  integer :: dummy
!
  real :: g0=1.,plaw=0.,ptlaw=1.
  logical :: lexponential_smooth=.false.
  real :: radial_percent_smooth=10.0,rshift=0.0
  logical :: lcorrect_selfgravity=.false.
  real :: gravitational_const
!
! For the noise
! 
  real :: rho_rms=0.05
  integer :: xmodes=10,ymodes=10,zmodes=0
  logical :: llowk_noise=.false.,lgaussian_distributed_noise=.true.
  real :: xmid=1.5
!
  namelist /initial_condition_pars/ g0,plaw,ptlaw,lexponential_smooth,&
       radial_percent_smooth,rshift,lcorrect_selfgravity,gravitational_const,&
       xmodes,ymodes,zmodes,rho_rms,llowk_noise,xmid,&
       lgaussian_distributed_noise
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
    subroutine initial_condition_all(f)
!
!  Initializes all the f arrays in one call.  This subroutine is called last.
!
!  21-dec-10/ccyang: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_all
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
      use Gravity, only: r0_pot,n_pot,acceleration,qgshear
      use Sub,     only: get_radial_distance,power_law
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: rr_cyl,rr_sph,OO,g_r,tmp
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
      do m=1,my
      do n=1,mz
!
        call get_radial_distance(rr_sph,rr_cyl)
!
        if (lgrav) then
!
! Gravity of a static central body
!
          call acceleration(g_r)
!
          if (any(g_r > 0.)) then
            do i=1,mx
              if (g_r(i) > 0) then
                if ((i <= nghost).or.(i >= mx-nghost)) then
                  !ghost zones, just emit a warning
                  if (ip<=7) then
                    print*,"centrifugal_balance: gravity at ghost point ",&
                         x(i),y(m),z(n),"is directed outwards"
                    call warning("","")
                  endif
                  OO(i)=0
                else
                  !physical point. Break!
                  print*,"centrifugal_balance: gravity at physical point ",&
                       x(i),y(m),z(n),"is directed outwards"
                  call fatal_error("","")
                endif
              else !g_r ne zero
                OO(i)=sqrt(-g_r(i)/rr_cyl(i))
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
          f(:,m,n,iux) = f(:,m,n,iux) - y(m)*OO
          f(:,m,n,iuy) = f(:,m,n,iuy) + x   *OO
          f(:,m,n,iuz) = f(:,m,n,iuz) + 0.
        elseif (coord_system=='cylindric') then
          f(:,m,n,iux) = f(:,m,n,iux) + 0.
          f(:,m,n,iuy) = f(:,m,n,iuy) + OO*rr_cyl
          f(:,m,n,iuz) = f(:,m,n,iuz) + 0.
        elseif (coord_system=='spherical') then
          f(:,m,n,iux) = f(:,m,n,iux) + 0.
          f(:,m,n,iuy) = f(:,m,n,iuy) + 0.
          f(:,m,n,iuz) = f(:,m,n,iuz) + OO*rr_sph
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
      use Mpicomm,         only: stop_it
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
      logical                :: lheader,lenergy,lpresent_zed
!
      if (lroot) print*,&
           'initial_condition_lnrho: locally isothermal approximation'
      if (lroot) print*,'Radial stratification with power law=',plaw
!
      lenergy=ltemperature.or.lentropy
!
      if (lenergy.and.llocal_iso) then
        if (lroot) then
          print*,'You switched on entropy or temperature evolution,'
          print*,'but you are still using llocal_iso in start.in.'
          print*,'Use one or the other instead.'
        endif
        call stop_it("")
      endif
!
!  Set the sound speed
!
      do m=1,my
      do n=1,mz
        lheader=((m==1).and.(n==1).and.lroot)
        call get_radial_distance(rr_sph,rr_cyl)
        if (lsphere_in_a_box.or.lspherical_coords) then
          rr=rr_sph
        elseif (lcylinder_in_a_box.or.lcylindrical_coords) then
          rr=rr_cyl
        else
          call stop_it("initial_condition_lnrho: "//&
               "no valid coordinate system")
        endif
!
        call power_law(cs20,rr,ptlaw,cs2,r_ref)
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
      enddo
      enddo
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
               + plaw*log((1-exp( -((rr-rshift)/rmid)**2 ))/rr)
        else
          lnrhomid=log(rho0)-.5*plaw*log((rr/r_ref)**2+rsmooth**2)
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
              call stop_it("local_isothermal_density")
            endif
!
            tmp2=-tmp1*rr_sph - cs2*(plaw + ptlaw)/gamma
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
              call stop_it("local_isothermal_density")
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
      call correct_pressure_gradient(f,ics2,ptlaw)
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
             (f,ptlaw,ics2,iglobal_cs2,iglobal_glnTT)
      else
        call set_thermodynamical_quantities(f,ptlaw,ics2)
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
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine correct_pressure_gradient(f,ics2,ptlaw)
!
!  Correct for pressure gradient term in the centrifugal force.
!  For now, it only works for flat (isothermal) or power-law
!  sound speed profiles, because the temperature gradient is
!  constructed analytically.
!
!  21-aug-07/wlad : coded
!
      use FArrayManager
      use Mpicomm,only: stop_it
      use Sub,    only: get_radial_distance,grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: glnrho
      real, dimension (nx)   :: rr,rr_cyl,rr_sph
      real, dimension (nx)   :: cs2,tmp1,tmp2,corr,gslnrho,gslnTT
      integer                :: i,ics2
      logical                :: lheader,lenergy
      real :: ptlaw
!
      if (lroot) print*,'Correcting density gradient on the '//&
           'centrifugal force'
!
      lenergy=ltemperature.or.lentropy
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
        gslnTT=-ptlaw/((rr/r_ref)**2+rsmooth**2)*rr/r_ref**2
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
              call stop_it("")
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
      use Mpicomm,     only:stop_it
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
          call stop_it("")
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
                call stop_it("")
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
    subroutine set_thermodynamical_quantities&
         (f,ptlaw,ics2,iglobal_cs2,iglobal_glnTT)
!
!  Subroutine that sets the thermodynamical quantities
!   - static sound speed, temperature or entropy -
!  based on a sound speed which is given as input.
!  This routine is not general. For llocal_iso (locally
!  isothermal approximation, the temperature gradient is
!  stored as a static array, as the (analytical) derivative
!  of an assumed power-law profile for the sound speed
!  (hence the parameter ptlaw)
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
      real :: cp1,ptlaw
      integer, pointer, optional :: iglobal_cs2,iglobal_glnTT
      integer :: ics2
      logical :: lenergy
!
      intent(in)  :: ptlaw
      intent(out) :: f
!
!  Break if llocal_iso is used with entropy or temperature
!
      lenergy=ltemperature.or.lentropy
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
      if (lroot) print*,'Temperature gradient with power law=',ptlaw
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
          gslnTT=-ptlaw/((rr/r_ref)**2+rsmooth**2)*rr/r_ref**2
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
      if (llocal_iso) then
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
      real :: Lx,Ly,Lz,d0,phase,xi,yi,zi,nw1
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
          lump_of_sines(i,m,n)=lump_of_sines(i,m,n) + &
             sin(2*pi*(ll*xi/Lx + mm*yi/Ly + nn*zi/Lz+ phase))
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
            f(ll1,mm1,nn1,irho)=f(ll1,mm1,nn1,irho) + &
                lump_of_sines(i,m,n)*exp(-(.5*(xi-xmid)/d0)**2)
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
