! $Id: baroclinic_run.f90 19193 2012-06-30 12:55:46Z wdobler $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
!
!***************************************************************

!-------------------------------------------------------------------
!
! HOW TO USE THIS FILE
! --------------------
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the 
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users ot optionally check their contributions
! in to the Pencil-Code CVS repository.  This may be useful if you
! are working on/using the additional physics with somebodyelse or 
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module 
! selections to say something like:
!
!    SPECIAL=special/nstar
!
! Where nstar it replaced by the filename of your new module
! upto and not including the .f90
!
!--------------------------------------------------------------------

!
!  This file adds a 2D turbulent potential to emulate the effects of
!  turbulence in a simulation, in a more realistic way than simply 
!  using alpha viscosity. The description can be found in the following papers
!  
!    Laughlin, G., Steinacker, A., & Adams, F. C. 2004, ApJ, 608, 489
!    Ogihara, M., Ida S., & Morbidelli, A. 2007, Icarus, 188, 522
!    Baruteau, C., & Lin, D. 2010, ApJ, 709, 759
!    Horn, R. B., Lyra, W., Mac Low, M.-M., & Sandor, Zs. 2012, ApJ, 750, 34
!  
!  03-oct-12/wlad: coded
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../special.h'
!
  real :: alpha=0.01, temperature_power_law=1.0
  logical :: lturbulent_force=.true.,lcalc_potturb=.true.
!
  namelist /special_init_pars/ alpha,lcalc_potturb
!   
  namelist /special_run_pars/ alpha,lcalc_potturb,lturbulent_force
!
  integer, parameter :: nmode_max = nygrid/2, nmode_min = 1
  real :: logmode_min,logmode_max
!
  type InternalPencils
     real, dimension(nx,3) :: gpotturb
     real, dimension(nx)   :: potturb
  endtype InternalPencils
!
  type (InternalPencils) :: q
!
  real, dimension(nmode_max) :: gauss_ampl, rcenter, phicenter
  real, dimension(nmode_max) :: omega_mode, tmode_init, radial_sigma_inv
  real, dimension(nmode_max) :: tmode_lifetime, tmode_lifetime_inv, tmode_age
  integer, dimension(nmode_max) :: mode_wnumber 
!
  real, dimension(mx) :: rad,amplitude_scaled
  real, dimension(my) :: phi
!
  integer :: ipotturb
  integer :: idiag_potturbm=0,idiag_potturbmax=0,idiag_potturbmin=0
  integer :: idiag_gpotturbx2m=0,idiag_gpotturby2m=0,idiag_gpotturbz2m=0
!
  contains

!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables 
!  which should be know to be able to evaluate
!
!  14-jul-09/wlad: coded
!
      use Cdata
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
           "$Id: baroclinic_run.f90 19193 2012-06-30 12:55:46Z wdobler $")
!
      call farray_register_auxiliary('potturb',ipotturb)
!
      if (lroot) then
        open(4,file=trim(datadir)//'/index_special.pro',status='replace')
        write(4,*) 'potturb  ',ipotturb
        close(4)
      endif
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  14-jul-09/wlad: coded
!
      use Cdata 
      use Mpicomm
      use EquationOfState, only: cs0
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(mx) :: Omega2
      real :: aspect_ratio, amplitude
      logical :: lstarting
!
!  Initialize the turbulent potential to zero. 
!
      f(:,:,:,ipotturb)=0.0
!
!  For readability, use rad and phi instead of x and y. Break if not cylindrical
!
      if (.not.lcylindrical_coords) then
        call fatal_error("initialize_special",&
             "turbulent potential coded for cylindrical coordinates only")
      else
        rad=x
        phi=y
      endif
!
!  Useful constants
!
      logmode_min=log(1.*nmode_min)
      logmode_max=log(1.*nmode_max)
      aspect_ratio = cs0
!
!  Use the scaling from Baruteau & Lin 2010,  that defines the amplitude of the 
!  potential as a function of the Shakura-Sunyayev alpha value. 
!
      amplitude = 8.5d-2 * aspect_ratio * sqrt(alpha)
!
!  Scale the amplitude by r**2 * Omega**2
!
      !Omega2 = g0/rad**3
      Omega2 = 1./rad**3
      amplitude_scaled = rad**2*Omega2 * amplitude
!
!  Set the potential at start time, for debugging purposes.
!
      if (lstarting) call special_before_boundary(f)
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!   14-jul-09/wlad: coded
!
      use Mpicomm
      use Sub, only: grad
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
!  The turbulence force from the potential, for the momentum equation.
!
      if (ldiagnos .and. (idiag_potturbm   /=0 .or. &
                          idiag_potturbmax /=0 .or. & 
                          idiag_potturbmin /=0) &
          ) q%potturb = f(l1:l2,m,n,ipotturb)
!
      if (lturbulent_force .or. (ldiagnos .and. (idiag_gpotturbx2m/=0  .or. &
                                                 idiag_gpotturby2m/=0  .or. &
                                                 idiag_gpotturbz2m/=0) &
                                ) &
          ) call grad(f,ipotturb,q%gpotturb)
!
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine special_before_boundary(f)
!
!  This subroutine calculates the full potential due to the turbulence.
!
!  03-oct-12/wlad: coded
!
      use Mpicomm
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx) :: lambda 
      integer :: k
!
!  This main loop sets the modes.
!
      if (lcalc_potturb) then 
!
!  Generate or update the full list of modes. 
!
        call update_mode_list
!
!  Loop the grid **locally** to sum the contribution of each mode. 
!
        do n=n1,n2;do m=1,my 
!
!  Calculate the mode structure. The modes dimensions are: 
!  Radius: Centred exponentially
!  Azimuth: Pure mode with integer number of wavelengths in 2pi. 
!           Azimuthal position moves according to the Keplerian rate.
!  Time: Grows and fades following a sine curve, of duration equal 
!        to its sound crossing time.
!
          do k=1,nmode_max
            lambda = gauss_ampl(k) * &
!radial gaussian dimension of mode
                 exp(-((rad-rcenter(k))*radial_sigma_inv(k))**2)* &
!azimuthal extent of mode..
                 cos(mode_wnumber(k)*phi(m) - phicenter(k) - &
!..with Keplerianly-shifting center
                 omega_mode(k)*tmode_age(k))*&
!mode grows and fades in time following a sine curve
                 sin(pi*tmode_age(k)*tmode_lifetime_inv(k))
!
!  Sum all the modes.
!
            f(:,m,n,ipotturb) = f(:,m,n,ipotturb) + amplitude_scaled*lambda
!
          enddo
        enddo;enddo
      endif
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine update_mode_list
!
      use Mpicomm, only: mpibcast_real,mpibcast_int
!
      integer :: k
!
!  Create and update the list only on one processor. In the 
!  first timestep it generates the whole mode list, but 
!  afterwards, modes are updated only sporadically. Therefore, 
!  the idle time is minimal, as the list does not need to be 
!  updated too often. Were that not the case, the mode list 
!  would need to be calculated in a load-balanced way. 
!
      if (lroot) then 
        do k=1,nmode_max 
!
          if (it == 1) then 
            call set_mode(k)
            tmode_age(k) = 0.
          else
!
!  Test if the mode has exceeded its lifetime. 
!  If so, replace it by a new random mode.
!
            tmode_age(k) = t-tmode_init(k)
            if (tmode_age(k) > tmode_lifetime(k)) then
!
              if (ldebug .or. ip<=9) print*,'mode ',k,' at r=',&
                   rcenter(k),', of m=',mode_wnumber(k),' gone'
!
              call set_mode(k)
!
              if (ldebug .or. ip<=9) print*,'   replaced by mode at r=',&
                   rcenter(k),', of m=',mode_wnumber(k)
            endif
          endif
        enddo
      endif
!
!  Unfortunately it cannot broadcast a structure. Do one by one. 
!  Also, look at a way to flag individual slots in the arrays and 
!  only broadcast those that were changed. 
!
      call mpibcast_real(gauss_ampl,nmode_max)
      call mpibcast_real(rcenter,nmode_max)
      call mpibcast_real(phicenter,nmode_max)
      call mpibcast_real(radial_sigma_inv,nmode_max)
      call mpibcast_real(tmode_init,nmode_max)
      call mpibcast_real(tmode_lifetime,nmode_max)
      call mpibcast_real(omega_mode,nmode_max)
      call mpibcast_real(tmode_lifetime_inv,nmode_max)
      call mpibcast_real(tmode_age,nmode_max)
      call mpibcast_int(mode_wnumber,nmode_max)
!
    endsubroutine update_mode_list
!***********************************************************************
    subroutine set_mode(k)
!
!  This sets all the needed parameters of the mode, and stores them in 
!  a "structure".
!
      use General, only: random_number_wrapper
!
      real :: gauss_ampl_scl,phicenter_scl
      real :: tmode_lifetime_scl,rcenter_scl,inv_radial_sigma_scl
      real :: tmode_init_scl,omega_mode_scl
      integer :: mode_wnumber_scl
!
      real :: aux1,aux2
      real :: cs_mode
      integer :: k
!
!  Mode's Gaussian-distributed random amplitude.
!
      call random_number_wrapper(aux1)
      call random_number_wrapper(aux2)
      gauss_ampl_scl   = sqrt(-2*log(aux1))*cos(2*pi*aux2)
!
!  Mode's random radial and azimuthal location.
!
      call random_number_wrapper(aux1)
      rcenter_scl    = aux1*(r_ext-r_int) + r_int
      call random_number_wrapper(aux1)
      phicenter_scl  = aux1*2*pi
!
!  Mode's azimuthal wavenumber. 
!
      call random_number_wrapper(aux1)
      aux1 = aux1*(logmode_max-logmode_min)+logmode_min
      mode_wnumber_scl = nint(exp(aux1))
!
!  Radial extent (sigma variance) of mode.
!
      inv_radial_sigma_scl = 4.*mode_wnumber_scl/(pi*rcenter_scl)
!
!  Keplerian frequency at mode center.
!
      !omega_mode_scl = g0*rcenter_scl**(-1.5)
      omega_mode_scl = rcenter_scl**(-1.5)
!
!  Mode's time of appearance.
!
      tmode_init_scl=t
!
!  Mode's duration, which is the sound crossing time.
!
      cs_mode = rcenter_scl**(-0.5*temperature_power_law)
      tmode_lifetime_scl = 2*pi*rcenter_scl/(mode_wnumber_scl*cs_mode)
!
!  Store them in pan-mode array. 
!
      gauss_ampl(k)       = gauss_ampl_scl
      rcenter(k)          = rcenter_scl
      phicenter(k)        = phicenter_scl
      mode_wnumber(k)     = mode_wnumber_scl
      radial_sigma_inv(k) = inv_radial_sigma_scl
      tmode_init(k)       = tmode_init_scl
      tmode_lifetime(k)   = tmode_lifetime_scl
      omega_mode(k)       = omega_mode_scl
!
!  Shortcut for inverse mode lifetime
!
      tmode_lifetime_inv(k) = 1./tmode_lifetime_scl
!
    endsubroutine set_mode
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=special_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=special_init_pars)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif
!
99    return
endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
      use Diagnostics
!
!  reads and registers print parameters relevant to special
!
!   14-jul-09/wlad: coded
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Write information to index.pro
!
      if (lreset) then
        idiag_potturbm=0;idiag_potturbmax=0;idiag_potturbmin=0 
        idiag_gpotturbx2m=0;idiag_gpotturby2m=0;idiag_gpotturbz2m=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'potturbm',idiag_potturbm)
        call parse_name(iname,cname(iname),cform(iname),'potturbmax',idiag_potturbmax)
        call parse_name(iname,cname(iname),cform(iname),'potturbmin',idiag_potturbmin)
        call parse_name(iname,cname(iname),cform(iname),'gpotturbx2m',idiag_gpotturbx2m)
        call parse_name(iname,cname(iname),cform(iname),'gpotturby2m',idiag_gpotturby2m)
        call parse_name(iname,cname(iname),cform(iname),'gpotturbz2m',idiag_gpotturbz2m)
      enddo
!
      if (lwr) then
        write(3,*) 'i_potturbm=',idiag_potturbm
        write(3,*) 'i_potturbmax=',idiag_potturbmax
        write(3,*) 'i_potturbmin=',idiag_potturbmin
        write(3,*) 'i_gpotturbx2m=',idiag_gpotturbx2m
        write(3,*) 'i_gpotturby2m=',idiag_gpotturby2m
        write(3,*) 'i_gpotturbz2m=',idiag_gpotturbz2m
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   momentum equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
      use Cdata
      use Diagnostics
!      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      integer :: j,ju
!
!  Modified momentum equation
!
      if (lturbulent_force) then
        do j=1,3 
          ju=j+iuu-1
          df(l1:l2,m,n,ju) = df(l1:l2,m,n,ju) - q%gpotturb(:,j)
        enddo
      endif
!
      if (ldiagnos) then 
        if (idiag_potturbm/=0)   call sum_mn_name(q%potturb,idiag_potturbm)
        if (idiag_potturbmax/=0) call max_mn_name(q%potturb,idiag_potturbmax)
        if (idiag_potturbmin/=0) call max_mn_name(-q%potturb,idiag_potturbmin,lneg=.true.)
        if (idiag_gpotturbx2m/=0) call sum_mn_name(q%gpotturb(:,1)**2,idiag_gpotturbx2m)
        if (idiag_gpotturby2m/=0) call sum_mn_name(q%gpotturb(:,2)**2,idiag_gpotturby2m)
        if (idiag_gpotturbz2m/=0) call sum_mn_name(q%gpotturb(:,3)**2,idiag_gpotturbz2m)
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
!***********************************************************************
!********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special

