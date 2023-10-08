! $Id$
!
!  This module includes physics related to planet atmospheres.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
! variables global to this module
!
  real, dimension(my,mz) :: mu_ss=0.
  real, dimension(my) :: lat
  real, dimension(mz) :: lon
  real :: pp2Pa=1., TT2K=1.  !  convert pc units to [Pa] and [K]
!
! variables in the reference profile
!
  real :: dlogp_ref, logp_ref_min, logp_ref_max
  real, dimension(:), allocatable :: p_temp_ref,tau_rad,temp_ref,logp_ref,dteq,Teq_night
  integer :: nref
!
! Init parameters
!
  real :: lon_ss=0., lat_ss=0.  ! unit: [degree]
  real :: Tbot=0., Ttop=100., dTeqbot=0., dTeqtop=100.  ! unit: [K]
  real :: peqtop=1.d2, peqbot=1.d6  ! unit: [Pa]
  real :: tauradtop=1.d4, tauradbot=1.d7  ! unit: [s]
  real :: pradtop=1.d3, pradbot=1.d6 ! unit:[ Pa]
!
!
!
  namelist /special_init_pars/ &
      lon_ss,lat_ss,Tbot,Ttop,peqtop,peqbot,tauradtop,tauradbot,&
      pradtop,pradbot
!
!
! Declare index of new variables in f array (if any).
!
!!   integer :: ispecial=0
!!   integer :: ispecaux=0
!
!! Diagnostic variables (needs to be consistent with reset list below).
!
!!   integer :: idiag_POSSIBLEDIAGNOSTIC=0
!
  contains
!****************************************************************************
  subroutine initialize_mult_special
!
! Dummy routine.
!
  endsubroutine initialize_mult_special
!***********************************************************************
    subroutine register_special
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
      if (lroot) call svn_id( &
           "$Id$")
!
!!      call farray_register_pde('special',ispecial)
!!      call farray_register_auxiliary('specaux',ispecaux)
!!      call farray_register_auxiliary('specaux',ispecaux,communicated=.true.)
!
    endsubroutine register_special
!***********************************************************************
    subroutine register_particles_special(npvar)
!
!  Set up indices for particle variables in special modules.
!
!  4-jan-14/tony: coded
!
      integer :: npvar
!
      if (lroot) call svn_id( &
           "$Id$")
      call keep_compiler_quiet(npvar)
!
!
!!      iqp=npvar+1
!!      npvar=npvar+1
!
    endsubroutine register_particles_special
!***********************************************************************
    subroutine initialize_special(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (.not.ltemperature_nolog) call fatal_error('initialize_special', &
          'special/planet_atmosphere is formulated in TT only')
!
!  convert y and z to latitude and longitude in rad
!
      lat=0.5*pi-y
      lon=z-pi
!
!  unit conversion
!
      call prepare_unit_conversion
!
!  read in the reference P-T profile, and calculate Teq and tau_rad etc.
!
      call prepare_Tref_and_tau
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine finalize_special(f)
!
!  Called right before exiting.
!
!  14-aug-2011/Bourdin.KIS: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine finalize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!!
!!  SAMPLE IMPLEMENTATION
!!
!!      select case (initspecial)
!!        case ('nothing'); if (lroot) print*,'init_special: nothing'
!!        case ('zero', '0'); f(:,:,:,iSPECIAL_VARIABLE_INDEX) = 0.
!!        case default
!!          call fatal_error("init_special: No such value for initspecial:" &
!!              ,trim(initspecial))
!!      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
    lpenc_requested(i_rho)=.true.
    lpenc_requested(i_pp)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils), intent(inout) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Special pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-nov-04/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt_ode

    endsubroutine dspecial_dt_ode
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  Several precalculated Pencils of information are passed for
!  efficiency.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!!      if (headtt) call identify_bcs('special',ispecial)
!
!!
!! SAMPLE DIAGNOSTIC IMPLEMENTATION
!!
!!      if (ldiagnos) then
!!        if (idiag_SPECIAL_DIAGNOSTIC/=0) then
!!          call sum_mn_name(MATHEMATICAL EXPRESSION,idiag_SPECIAL_DIAGNOSTIC)
!!! see also integrate_mn_name
!!        endif
!!      endif
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
!!      use FArrayManager, only: farray_index_append
!
!!      integer :: iname
      logical :: lreset,lwrite
!!!
!!!  reset everything in case of reset
!!!  (this needs to be consistent with what is defined above!)
!!!
      if (lreset) then
!!        idiag_SPECIAL_DIAGNOSTIC=0
      endif
!!
!!      do iname=1,nname
!!        call parse_name(iname,cname(iname),cform(iname),&
!!            'NAMEOFSPECIALDIAGNOSTIC',idiag_SPECIAL_DIAGNOSTIC)
!!      enddo
!!
!!!  write column where which magnetic variable is stored
!!      if (lwr) then
!!        call farray_index_append('idiag_SPECIAL_DIAGNOSTIC',idiag_SPECIAL_DIAGNOSTIC)
!!      endif
!!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of Special variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  momentum equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  continuity equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_dustdensity(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  continuity equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_dustdensity
!***********************************************************************
    subroutine special_calc_energy(f,df,p)
!
!  08-sep-23/hongzhe: specific things for planet atmospheres.
!                     Waiting Kuan to fill in.
!                     Reference: Rogers & Komacek (2014)
!  27-sep-23/hongzhe: outsourced from temperature_idealgas.f90
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx) :: Teq_x,tau_rad_x,log10pp
      real ,dimension(:), allocatable :: Teq_local
      integer :: ix,index
!
!  The local equilibrium T; still in pressure coordinate
!
      allocate(Teq_local(nref))
      Teq_local = Teq_night + dteq*max(0.,mu_ss(nghost+m,nghost+n))
!
!  interpolation for Teq_local at height coordinate (could have arbitrary pressure)
!
!!! %HZ: need to convert p%pp to [Pa]
        log10pp = log10(p%pp*pp2Pa)
        do ix=l1,l2
          ! index of the logp_ref that is just smaller than log10(pressure)
          index = ceiling(log10pp(ix)-logp_ref_min)/dlogp_ref
          if (index>=nref) then
            Teq_x(ix) = Teq_local(nref)
            tau_rad_x(ix) = tau_rad(nref)
          elseif (index<= 1) then 
            Teq_x(ix) = Teq_local(1)
            tau_rad_x(ix) = tau_rad(1)
          else
            ! interpolate T in log10(p) space
            Teq_x(ix) = Teq_local(index)+(Teq_local(index+1)-Teq_local(index))*   &
                            (log10pp(ix)-logp_ref(index))/   &
                            (logp_ref(index+1)-logp_ref(index))   ! unit of K
            ! interpolate log10(tau_rad) in log10 p space
            tau_rad_x(ix) = log10(tau_rad(index))+  &
                                (log10(tau_rad(index+1))-log10(tau_rad(index)))*   &
                                (log10pp(ix)-logp_ref(index))/   &
                                (logp_ref(index+1)-logp_ref(index))
            tau_rad_x(ix) = 10**tau_rad_x(ix)  ! unit of s
          endif
        enddo
!
      df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) - (p%TT-Teq_x)/tau_rad_x
!
      deallocate(Teq_local)
!
    endsubroutine special_calc_energy
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  induction equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_magnetic
!***********************************************************************
    subroutine special_calc_pscalar(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  passive scalar equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array.
!
!  15-jun-09/anders: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,ilncc) = df(l1:l2,m,n,ilncc) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_pscalar
!***********************************************************************
    subroutine special_particles_bfre_bdary(f,fp,ineargrid)
!
!  Called before the loop, in case some particle value is needed
!  for the special density/hydro/magnetic/entropy.
!
!  20-nov-08/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (:,:), intent(in) :: fp
      integer, dimension(:,:) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine special_particles_bfre_bdary
!***********************************************************************
        subroutine special_calc_particles(f,df,fp,dfp,ineargrid)
!
!  Called before the loop, in case some particle value is needed
!  for the special density/hydro/magnetic/entropy.
!
!  20-nov-08/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(in) :: df
      real, dimension (:,:), intent(in) :: fp,dfp
      integer, dimension(:,:) :: ineargrid
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(fp,dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine special_calc_particles
!***********************************************************************
    subroutine special_calc_chemistry(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  induction equation.
!
!
!  15-sep-10/natalia: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_chemistry
!***********************************************************************
    subroutine special_calc_spectra(f,spectrum,spectrumhel,lfirstcall,kind)

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:) :: spectrum,spectrumhel
      logical :: lfirstcall
      character(LEN=3) :: kind

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(spectrum,spectrumhel)
      call keep_compiler_quiet(lfirstcall)
      call keep_compiler_quiet(kind)
  
    endsubroutine special_calc_spectra
!***********************************************************************
    subroutine special_calc_spectra_byte(f,spectrum,spectrumhel,lfirstcall,kind,len)

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:) :: spectrum,spectrumhel
      logical :: lfirstcall
      integer(KIND=ikind1), dimension(3) :: kind
      integer :: len

      !call keep_compiler_quiet(char(kind))
      call keep_compiler_quiet(len)
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(spectrum,spectrumhel)
      call keep_compiler_quiet(lfirstcall)
      
    endsubroutine special_calc_spectra_byte
!***********************************************************************
    subroutine special_before_boundary(f)
!
!  Possibility to modify the f array before the boundaries are
!  communicated.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
!  compute cos(angle between the substellar point)
!
      call get_mu_ss(mu_ss,lon_ss,lat_ss)
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine special_boundconds(f,bc)
!
!  Some precalculated pencils of data are passed in for efficiency,
!  others may be calculated directly from the f array.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (boundary_condition), intent(in) :: bc
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(bc)
!
    endsubroutine special_boundconds
!***********************************************************************
    subroutine special_after_timestep(f,df,dt_,llast)
!
!  Possibility to modify the f and df after df is updated.
!  Used for the Fargo shift, for instance.
!
!  27-nov-08/wlad: coded
!
      logical, intent(in) :: llast
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, intent(in) :: dt_
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dt_)
      call keep_compiler_quiet(llast)
!
    endsubroutine  special_after_timestep
!***********************************************************************
    subroutine special_particles_after_dtsub(f, dtsub, fp, dfp, ineargrid)
!
!  Possibility to modify fp in the end of a sub-time-step.
!
!  28-aug-18/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, intent(in) :: dtsub
      real, dimension(:,:), intent(in) :: fp, dfp
      integer, dimension(:,:), intent(in) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(dtsub)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine special_particles_after_dtsub
!***********************************************************************
    subroutine set_init_parameters(Ntot,dsize,init_distr,init_distr2)
!
!  Possibility to modify the f and df after df is updated.
!  Used for the Fargo shift, for instance.
!
!  27-nov-08/wlad: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(ndustspec) :: dsize,init_distr,init_distr2
      real :: Ntot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(dsize,init_distr,init_distr2)
      call keep_compiler_quiet(Ntot)
!
    endsubroutine  set_init_parameters
!***********************************************************************
    subroutine prepare_unit_conversion
!
!  Constants that convert code units to SI.
!
!  28-sep-23/hongzhe: coded
!
      use EquationOfState, only: getmu
!
      real :: mu  !  mean molecular weight
!
      call getmu(mu_tmp=mu)
!
      if (unit_system=='SI') then
        ! %HZ: need to code correctly!!
        ! m_u is atomic mass unit, not mean molecular weight
        pp2Pa=k_B_cgs/m_u_cgs*1.0e-4/mu*unit_density*unit_temperature
        TT2K=1.*unit_temperature
      else
        call fatal_error('prepare_unit_conversion','please use SI system')
      endif
!
    endsubroutine  prepare_unit_conversion
!***********************************************************************
    subroutine prepare_Tref_and_tau
!
!   Read the reference profile and calculated tau_rad.
!   All quantities in this subroutine are in SI units.
!
!   28-sep-23/xianyu,hongzhe: coded
!
    real :: alpha
    integer :: i
    logical :: lTref_file_exists
!
!  read in Tref, in physical unit.
!
      inquire(FILE='iro-teq-tint100K-regrid-Pa.txt', EXIST=lTref_file_exists)
      if (.not.lTref_file_exists) call fatal_error('initialize_special', &
          'Must provide a Tref file')
!
      open(1,file='iro-teq-tint100K-regrid-Pa.txt')
      read(1,*)
      read(1,*) dlogp_ref, logp_ref_min, logp_ref_max  ! in log10(bar)
      read(1,*) nref
      if(allocated(logp_ref)) deallocate(logp_ref)
      if(allocated(temp_ref)) deallocate(temp_ref)
      allocate(logp_ref(nref),temp_ref(nref))
      read(1,*) logp_ref,temp_ref
      if (lroot) then
        print*, 'nref=',nref
        print *,'Here is the baseline radiative equil T profile'
        print *,'used in the Newtonian relaxation:'
        print *,'p [Pa] and T_eq[K]:'
        do i=1,nref
          print*, 'logp_ref,temp_ref=',logp_ref(i),temp_ref(i)
        enddo
      endif
      close(1)
!
!  convert units, and calculate tau_rad, dteq, and Teq_night
!
      if(allocated(p_temp_ref)) deallocate(p_temp_ref)
      if(allocated(tau_rad))    deallocate(tau_rad)
      if(allocated(dteq))       deallocate(dteq)
      if(allocated(Teq_night))  deallocate(Teq_night)
      allocate( p_temp_ref(nref), tau_rad(nref), dteq(nref),Teq_night(nref) )
      p_temp_ref = 10.**logp_ref
!
      alpha=log(tauradtop/tauradbot)/log(pradtop/pradbot)
      where (p_temp_ref>=pradtop .and. p_temp_ref<=pradbot)
        tau_rad = tauradbot * (p_temp_ref/pradbot)**alpha
        dteq = dTeqbot + (dTeqtop - dTeqbot) * log(p_temp_ref/peqbot)/log(peqtop/peqbot)
      elsewhere (p_temp_ref < pradtop)
        tau_rad = tauradtop
        dteq = dTeqtop
      elsewhere
        tau_rad = tauradbot
        dteq = dteqbot
      endwhere
!
      Teq_night = temp_ref - 0.5*dteq
!
      deallocate(p_temp_ref)
!
    endsubroutine  prepare_Tref_and_tau
!***********************************************************************
    subroutine get_mu_ss(mu_ss,lonss,latss)
!
!  Given the substellar longitude lonss and latitude latss (both in
!  degree), compute the cosine of the angle of the sun from overhead.
!  This will be used in future RT calculation.
!  It also tells you whether you're on the dayside, since
!  mu_ss>0 for dayside, mu_ss=0 on the terminator, and mu_ss<0 on the nightside.
!
!  28-sep-23/xianyu,hongzhe: coded
!
      real, dimension(my,mz), intent(out) :: mu_ss
      real, intent (in) :: lonss, latss
!
      real, PARAMETER :: deg2rad=pi/180.
!
      mu_ss = spread(cos(lon),1,my) * spread(cos(lat),2,mz) &
          * cos(lonss*deg2rad) * cos(latss*deg2rad) + &
           spread(sin(lon),1,my)*spread(cos(lat),2,mz) &
          * sin(lonss*deg2rad) * cos(latss*deg2rad) + &
           spread(sin(lat),2,mz)*sin(latss*deg2rad)
!
    endsubroutine  get_mu_ss
!***********************************************************************
endmodule Special
