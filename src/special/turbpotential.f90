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
  logical :: lcalc_potturb=.true.              ! Calculate the turbulent potential
  logical :: lturbulent_force=.true.           ! Add the turbulent potential to the momentum equation
  logical :: ltime_dependant_amplitude=.true.  ! Term that makes the modes come in and out of existence smoothly
  logical :: lgravitational_turbulence=.false. ! Choose between asp ratio 4 for the modes (GI) or 1/h (MRI)
  logical :: lcap_modes_at_m6=.false.          ! Set amplitude of modes higher than m>6 to zero. 
  logical :: lupdate_as_var=.false.            ! Writes modes.dat and MODESN like var.dat and VARN
                                               !  though symmetric, it is not that useful, since the 
                                               !  modelist changes only very sporadically (the modes are 
                                               !  long lived).  
!
  real :: rmodes_int=impossible
  real :: rmodes_ext=impossible                ! emulate a dead zone
!
  namelist /special_init_pars/ alpha,lcalc_potturb,lturbulent_force,&
       ltime_dependant_amplitude,lgravitational_turbulence,&
       lcap_modes_at_m6,lupdate_as_var,rmodes_int,rmodes_ext
!
  namelist /special_run_pars/ alpha,lcalc_potturb,lturbulent_force,&
       ltime_dependant_amplitude,lgravitational_turbulence,&
       lcap_modes_at_m6,lupdate_as_var,rmodes_int,rmodes_ext
!
  integer, parameter :: nmode_max = 50
  integer :: mmode_max=nygrid/8, mmode_min = 1
  real :: logmode_min,logmode_max,cs01
  logical :: lcompute_all_modes
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
  real, dimension(nmode_max) :: tmode_lifetime, tmode_lifetime_inv
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
      logmode_min=log(1.*mmode_min)
      logmode_max=log(1.*mmode_max)
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
!  Inverse sound speed
!
      cs01=1./cs0
!
!  Switch for a more intuitive label. 
!
      if (lcap_modes_at_m6) then 
        lcompute_all_modes=.false.
      else
        lcompute_all_modes=.true.
      endif
!
!  Radii to bracket the modes
!
      if (rmodes_int==impossible) rmodes_int=r_int
      if (rmodes_ext==impossible) rmodes_ext=r_ext
!
!  Set the potential at start time, for debugging purposes.
!
      if (lstarting) call special_before_boundary(f,lstarting)
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
    subroutine special_before_boundary(f,lstarting)
!
!  This subroutine calculates the full potential due to the turbulence.
!
!  03-oct-12/wlad: coded
!
      use Mpicomm
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx) :: lambda, time_dependant_amplitude
      real ::  tmode_age
      integer :: k,icount
      logical, optional :: lstarting
      logical :: lstart
!
      if (present(lstarting)) then 
        lstart = lstarting
      else
        lstart=.false.
      endif
!
!  This main loop sets the modes.
!
      if (lcalc_potturb) then 
!
!  Generate or update the full list of modes. 
!
        call update_modes(lstart)
!
!  Loop the grid **locally** to sum the contribution of each mode. 
!
        icount=0
        do k=1,nmode_max
          if (lcompute_all_modes .or. mode_wnumber(k) <= 6) then  
!
            icount=icount+1
!
            tmode_age = t-tmode_init(k)
            do n=n1,n2;do m=1,my 
!
!  Calculate the mode structure. The modes dimensions are: 
!  Radius: Centred exponentially
!  Azimuth: Pure mode with integer number of wavelengths in 2pi. 
!           Azimuthal position moves according to the Keplerian rate.
!  Time: Grows and fades following a sine curve, of duration equal 
!        to its sound crossing time.
!
              if (ltime_dependant_amplitude) then 
                time_dependant_amplitude = &
                     sin(pi*tmode_age*tmode_lifetime_inv(k))
              else
                time_dependant_amplitude = 1.
              endif
!
              lambda = gauss_ampl(k) * &
!radial gaussian dimension of mode
                   exp(-((rad-rcenter(k))*radial_sigma_inv(k))**2)* &
!azimuthal extent of mode..
                   cos(mode_wnumber(k)*phi(m) - phicenter(k) - &
!..with Keplerianly-shifting center
                   omega_mode(k)*tmode_age)* &
!mode grows and fades in time following a sine curve
                   time_dependant_amplitude
!
!  Sum all the modes.
!
              if (icount==1) then
                f(:,m,n,ipotturb) = amplitude_scaled*lambda
              else
                f(:,m,n,ipotturb) = f(:,m,n,ipotturb) + amplitude_scaled*lambda
              endif
!
            enddo;enddo
          endif
        enddo
      endif
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine update_modes(lstart)
!
      use Mpicomm, only: mpibcast_logical,mpibcast_real,mpibcast_int
!
!  Auxiliary array for the broadcast only. 
!
      real, dimension(8) :: g
      integer  :: h
!
      real :: tmode_age
      logical :: flag,lstart
      integer :: k
!
!  For file inquiry
!
      logical :: exist
      integer :: stat
!
!  Create and update the list only on one processor. In the 
!  first timestep it generates the whole mode list, but 
!  afterwards, modes are updated only sporadically. Therefore, 
!  the idle time is minimal, as the list does not need to be 
!  updated too often. Were that not the case, the mode list 
!  would need to be calculated in a load-balanced way. 
!
      if (lstart) then
        do k=1,nmode_max 
          if (lroot) then
            if (k==1) then 
              print*,''
              print*,'it=',it,' t=',t
              print*,'Start of the simulation. Setting the modes'
              print*,''
            endif
            call get_mode(g,h)
          endif
          call mpibcast_real(g,8)
          call mpibcast_int(h)
          call set_values(g,h,k)
        enddo
!
        if (lroot) then 
          call output_modes(trim(datadir)//'/modes.dat')
          if (lwrite_ic.and.lupdate_as_var) &
               call output_modes(trim(datadir)//'/MODES0')
        endif
      else
!
!  Restarting a simulation. Read the list of modes. All procs 
!  read the file and set their modes (would it be more efficient
!  to only have the root read, then broadcast?).
!
        if (it==1.and.itsub==1) then 
!
          if (lroot) then
            print*,''
            print*,'it=',it,' t=',t
            print*,'This is a simulation that is '
            print*,'restarting from a snapshot. The'
            print*,'list of modes will be read from file, '
            print*,'to regenerate the mode structure.'
            print*,''
            inquire(file=trim(datadir)//'/modes.dat',exist=exist)
            if (exist) then
              open(19,file=trim(datadir)//'/modes.dat')
            else
              call fatal_error(trim(datadir)//'/modes.dat','no input file')
            endif
          endif
!
          do k=1,nmode_max
            if (lroot) read(19,*,iostat=stat)    &
                   g(1),g(2),g(3),g(4),g(5),g(6),g(7),g(8),h
            call mpibcast_real(g,8)
            call mpibcast_int(h)
            call set_values(g,h,k)
          enddo
          if (lroot) close(19)
!
        else 
!
!  Test if the mode has exceeded its lifetime. 
!  If so, replace it by a new random mode.
!
          do k=1,nmode_max
            if (lroot) then
              flag=.false.
              tmode_age = t-tmode_init(k)
              if (tmode_age > tmode_lifetime(k)) then
                if (ldebug .or. ip<=9) print*,'mode ',k,' at r=',&
                     rcenter(k),', of m=',mode_wnumber(k),' gone'
!
                call get_mode(g,h)
                flag=.true.
!
                if (ldebug .or. ip<=9) print*,'   replaced by mode at r=',&
                     rcenter(k),', of m=',mode_wnumber(k)
              endif
            endif
            call mpibcast_logical(flag)
            if (flag) then 
!
              if (lroot) then
                if (ldebug .or. ip<=9) &
                     print*,'update_modes: saving the current modes before updating'
                call output_modes(trim(datadir)//'/modes_old.dat')
              endif
!
              call mpibcast_real(g,8)
              call mpibcast_int(h)
              call set_values(g,h,k)
!
!  Update the list of modes as well. Save the previous one. 
!
              if (lroot) then
                if (ldebug .or. ip<=9) &
                     print*,'update_modes: updating the mode output list '
                call output_modes(trim(datadir)//'/modes.dat')
              endif
!
            endif
          enddo
!
        endif
      endif
!
    endsubroutine update_modes
!***********************************************************************
    subroutine wsnap_mode
!
      use General, only: itoa
      use Sub, only: read_snaptime
!
      real, save :: tsnap
      integer, save :: nsnap
      logical, save :: lfirstcall=.true.
      character (len=intlen) :: insnap
!
!  The usual var.dat
!
      if (mod(it,isave)==0) &
           call output_modes(trim(datadir)//'/modes.dat')
!
!  The enumerated VARN
!
      if (lfirstcall) then
        nsnap=floor(t/dsnap)
        tsnap=dsnap*(nsnap+1)
        lfirstcall=.false.
      endif
      if (t >= tsnap) then 
        tsnap = tsnap + dsnap
        nsnap = nsnap + 1
        insnap=itoa(nsnap)
        call output_modes(trim(datadir)//'/MODES'//trim(insnap))
      endif
!
    endsubroutine wsnap_mode
!***********************************************************************
    subroutine output_modes(file)
!
      character (len=*) :: file
      integer :: k
!
      open(18,file=file)
      do k=1,nmode_max 
        write(18,*) gauss_ampl(k),rcenter(k),phicenter(k),&
             radial_sigma_inv(k),tmode_init(k),tmode_lifetime(k),&
             omega_mode(k),tmode_lifetime_inv(k),mode_wnumber(k)
      enddo
      close(18)
!
    endsubroutine output_modes
!***********************************************************************
    subroutine get_mode(g,h)
!
!  This sets all the needed parameters of the mode, and stores them in 
!  a "structure".
!
      !use General, only: random_number_wrapper
!
      real, dimension(8), intent(out) :: g
      integer, intent(out) :: h
!
      real :: gauss_ampl_scl,rcenter_scl,phicenter_scl
      real :: tmode_lifetime_scl,inv_radial_sigma_scl,cs1_mode
      real :: tmode_init_scl,omega_mode_scl,mode_aspect_ratio
      integer :: mode_wnumber_scl
!
      real :: aux1,aux2
!
!  Mode's azimuthal wavenumber. 
!
      call random_number(aux1)
      aux1 = aux1*(logmode_max-logmode_min)+logmode_min
      mode_wnumber_scl = nint(exp(aux1))
!
!  If the mode is smaller than m=6, we set the amplitude to zero; still 
!  we need to calculate the lifetime to know when to replace a mode. 
!
!  Mode's radial location
!
      call random_number(aux1)
      rcenter_scl    = aux1*(rmodes_ext-rmodes_int) + rmodes_int
!
!  Sound speed at mode location
!
      cs1_mode = cs01*rcenter_scl**(0.5*temperature_power_law)
!
!  Mode's time of appearance and duration, which is the sound crossing time.
!
      tmode_init_scl=t
      tmode_lifetime_scl = 2*pi*rcenter_scl*cs1_mode/mode_wnumber_scl
!
!  The rest only needs be calculated if the mode is computed, i.e., m<=6.
!
      if (lcompute_all_modes.or.mode_wnumber_scl <= 6) then 
!
!  Mode's Gaussian-distributed random amplitude.
!
        call random_number(aux1)
        call random_number(aux2)
        gauss_ampl_scl   = sqrt(-2*log(aux1))*cos(2*pi*aux2)
!
!  Mode's random radial and azimuthal location.
!
        call random_number(aux1)
        phicenter_scl  = aux1*2*pi
!
!  Keplerian frequency at mode location.
!
        omega_mode_scl = rcenter_scl**(-1.5)
!
!  Radial extent (sigma variance) of mode. In the original Laughlin
!  paper, the factor 4 was so that the modes had a 4:1 aspect ratio. 
!  This leads to a turbulence that in many ways ressembles that 
!  of self-gravity, with well defined spiral arms and low-m swinging modes. 
!  This perhaps comes from his background of massive disks. MRI turbulence
!  on the other hand, being subsonic, will give rise to structure that is 
!  more localized, with modes bound by the sonic radii. It is better to 
!  have the aspect ratio defined by the sound speed itself. 
!     
!  Aspect ratio at mode center.
!
        if (lgravitational_turbulence) then 
          mode_aspect_ratio = 4.
        else
          !inverse aspect ratio of the disk itself, 1/h
          mode_aspect_ratio = rcenter_scl*omega_mode_scl*cs1_mode
        endif
!
        inv_radial_sigma_scl = &
             mode_aspect_ratio*mode_wnumber_scl/(pi*rcenter_scl)
          
      else
!
!  Modes of m>6. Set amplitude at zero. Whatever for the rest.
!
        gauss_ampl_scl       = 0.
        phicenter_scl        = 1.
        omega_mode_scl       = 1.
        inv_radial_sigma_scl = 1.
!
      endif
!
!  Store them in the "broadcastable" array. 
!
      g(1) = gauss_ampl_scl
      g(2) = rcenter_scl
      g(3) = phicenter_scl
      g(4) = inv_radial_sigma_scl
      g(5) = tmode_init_scl
      g(6) = tmode_lifetime_scl
      g(7) = omega_mode_scl
      g(8) = 1./tmode_lifetime_scl
      h    = mode_wnumber_scl
!
    endsubroutine get_mode
!***********************************************************************
    subroutine set_values(g,h,k)
!
      real, dimension(8) :: g
      integer :: h,k
!
      gauss_ampl(k)          = g(1) 
      rcenter(k)             = g(2) 
      phicenter(k)           = g(3) 
      radial_sigma_inv(k)    = g(4) 
      tmode_init(k)          = g(5) 
      tmode_lifetime(k)      = g(6) 
      omega_mode(k)          = g(7) 
      tmode_lifetime_inv(k)  = g(8) 
      mode_wnumber(k)        = h    
      
    endsubroutine set_values
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
    subroutine special_after_timestep(f,df,dt_)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dt_
!
      if (lupdate_as_var.and.lroot) call wsnap_mode
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dt_)
!
    endsubroutine special_after_timestep
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Potential
!
        case('potturb')
          slices%yz=f(ix_loc,m1:m2,n1:n2,ipotturb)
          slices%xz=f(l1:l2,iy_loc,n1:n2,ipotturb)
          slices%xy=f(l1:l2,m1:m2,iz_loc,ipotturb)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,ipotturb)
          slices%ready = .true.
        endselect
!
    endsubroutine get_slices_special
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

