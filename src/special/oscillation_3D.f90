! $Id$
!
!  This module provide a way for users to specify custom
!  (i.e. not in the standard Pencil Code) physics, diagnostics etc.
!
!  The module provides a set of standard hooks into the Pencil-Code and
!  currently allows the following customizations:
!
!  Description                                     | Relevant function call
!  ---------------------------------------------------------------------------
!  Special variable registration                   | register_special
!    (pre parameter read)                          |
!  Special variable initialization                 | initialize_special
!    (post parameter read)                         |
!  Special variable finalization                   | finalize_special
!    (deallocation, etc.)                          |
!                                                  |
!  Special initial condition                       | init_special
!   this is called last so may be used to modify   |
!   the mvar variables declared by this module     |
!   or optionally modify any of the other f array  |
!   variables.  The latter, however, should be     |
!   avoided where ever possible.                   |
!                                                  |
!  Special term in the mass (density) equation     | special_calc_density
!  Special term in the momentum (hydro) equation   | special_calc_hydro
!  Special term in the energy equation             | special_calc_energy
!  Special term in the induction (magnetic)        | special_calc_magnetic
!     equation                                     |
!                                                  |
!  Special equation                                | dspecial_dt
!    NOT IMPLEMENTED FULLY YET - HOOKS NOT PLACED INTO THE PENCIL-CODE
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of oscillation_3D_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 2
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
! HOW TO USE THIS FILE
! --------------------
!
! Change the line above to
!   lspecial = .true.
! to enable use of special hooks.
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users ot optionally check their contributions
! in to the Pencil-Code SVN repository.  This may be useful if you
! are working on/using the additional physics with somebodyelse or
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module
! selections to say something like:
!
!   SPECIAL=special/geo_kws
!
! Where geo_kws it replaced by the filename of your new module
! upto and not including the .f90
!
module oscillation_3D
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
! Declare index of new variables in f array (if any).
!
  real :: ampl_psi=0., initpower_psi=0.,  initpower2_psi=0., kpeak_psi=0.
  real :: ampl_dpsi=0., initpower_dpsi=0.,  initpower2_dpsi=0., kpeak_dpsi=0.
  real :: relhel_psi=0., cutoff_psi=0.,  ncutoff_psi=1.
  real :: kgaussian_psi=0.,kgaussian_dpsi=0.
  logical :: lno_noise_psi=.false.
  logical :: lskip_projection_psi=.false., lvectorpotential=.false.
  logical :: lscale_tobox_psi=.true., lscale_tobox_dpsi=.true.
!
  real :: modulation_fact=1.
  character (len=labellen), dimension(ninit) :: init_psi='nothing'
  namelist /oscillation_3D_init_pars/  init_psi, &
    ampl_psi, initpower_psi, initpower2_psi, kpeak_psi, &
    ampl_dpsi, initpower_dpsi, initpower2_dpsi, kpeak_dpsi, &
    lscale_tobox_psi
!
  ! run parameters
  namelist /oscillation_3D_run_pars/ &
    modulation_fact
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_psirms=0,idiag_dpsirms=0
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
      use FArrayManager
!
!  6-oct-03/tony: coded
!
      if (lroot) call svn_id( &
           "$Id$")
!
      call farray_register_pde('ispecial',ispecialvar,array=2)
      ispecialvar2=ispecialvar+1
!
    endsubroutine register_special
!***********************************************************************
    subroutine register_particles_special(npvar)
!
!  Set up indices for particle variables in special modules.
!
      use FArrayManager
!
!  4-jan-14/tony: coded
!
      integer :: npvar
!
      if (lroot) call svn_id( &
           "$Id$")
!
      call farray_register_pde('ispecial',ispecialvar,array=2)
      ispecialvar2=ispecialvar+1
!
    endsubroutine register_particles_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
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
      use Initcond, only: power_randomphase_hel
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: kx
      integer :: j
!
      intent(inout) :: f
!
!  initial condition
!
      do j=1,ninit
        select case (init_psi(j))
          case ('nothing'); if (lroot) print*,'init_psi: nothing'
          case ('set')
            f(:,:,:,ispecialvar)=f(:,:,:,ispecialvar)+ampl_psi
            f(:,:,:,ispecialvar2)=f(:,:,:,ispecialvar2)+ampl_dpsi
          case ('sincos')
            kx=2*pi/Lx
            do m=1,my
            do n=1,mz
              f(:,m,n,ispecialvar)=ampl_psi*sin(kx*x)
              f(:,m,n,ispecialvar2)=ampl_dpsi*cos(kx*x)
            enddo
            enddo
!
!  spectrum
!
            case ('psi_power_randomphase')
              call power_randomphase_hel(ampl_psi, initpower_psi, initpower2_psi, &
                cutoff_psi, ncutoff_psi, kpeak_psi, f, ispecialvar, ispecialvar, &
                relhel_psi, kgaussian_psi, lskip_projection_psi, lvectorpotential, &
                lscale_tobox_psi, lpower_profile_file=.false., lno_noise=lno_noise_psi)
            case ('dpsi_power_randomphase')
              call power_randomphase_hel(ampl_dpsi, initpower_dpsi, initpower2_dpsi, &
                cutoff_psi, ncutoff_psi, kpeak_dpsi, f, ispecialvar2, ispecialvar2, &
                relhel_psi, kgaussian_psi, lskip_projection_psi, lvectorpotential, &
                lscale_tobox_psi, lpower_profile_file=.false., lno_noise=lno_noise_psi)
!
          case default
            !
            !  Catch unknown values
            !
            if (lroot) print*,'init_psi: No such value for init_psi: ', trim(init_psi(j))
            call fatal_error("init_psi","init value not defined")
        endselect
      enddo
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
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
!  Solve wave equation in 3-D with wave speed altered by
!  the strain tensor components.
!
      use Diagnostics
      use Mpicomm
      use Deriv, only: der2,derij
      use FArrayManager, only: farray_index_by_name
      use Sub
!
!  05-jan-25/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: tmp
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
      integer :: ih11_realspace, ih22_realspace, ih33_realspace
      integer :: ih12_realspace, ih23_realspace, ih31_realspace
      integer :: i, j, ij
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
!  Determine indices from other module
!
        ih11_realspace=farray_index_by_name('h11_realspace')
        ih22_realspace=farray_index_by_name('h22_realspace')
        ih33_realspace=farray_index_by_name('h33_realspace')
        ih12_realspace=farray_index_by_name('h12_realspace')
        ih23_realspace=farray_index_by_name('h23_realspace')
        ih31_realspace=farray_index_by_name('h31_realspace')
!
!  Solve wave equation d2psi/dt = del2 psi.
!  Begin with dpsi/dt = dpsi and solve further below ddpsi/dt = del2 psi.
!
        df(l1:l2,m,n,ispecialvar)=df(l1:l2,m,n,ispecialvar)+f(l1:l2,m,n,ispecialvar2)
!
!  Do the following only if at least the first of 6 indices exists.
!
        if (ih11_realspace>0) then
          do ij=1,6
!
!  Select i and j, and the target location ihij in the f-array.
!
            select case (ij)
              case (1); i=1; j=1; ihij=ih11_realspace
              case (2); i=2; j=2; ihij=ih22_realspace
              case (3); i=3; j=3; ihij=ih33_realspace
              case (4); i=1; j=2; ihij=ih12_realspace
              case (5); i=2; j=3; ihij=ih23_realspace
              case (6); i=3; j=1; ihij=ih31_realspace
            endselect
!
!  For diagonal components, add (delta_ij+h_ij)*d2f/dx_i dx_j
!  For each pair of off-diagonal components, add 2h_ij*d2f/dx_i dx_j.
!
            if (i==j) then
              call der2(f,ispecialvar,tmp,i)
              df(l1:l2,m,n,ispecialvar2)=df(l1:l2,m,n,ispecialvar2)+(1.+modulation_fact*f(l1:l2,m,n,ihij))*tmp
            else
              call derij(f,ispecialvar,tmp,i,j)
              df(l1:l2,m,n,ispecialvar2)=df(l1:l2,m,n,ispecialvar2)+2.*modulation_fact*f(l1:l2,m,n,ihij)*tmp
            endif
          enddo
        endif
!
!  diagnostics
!
      if (ldiagnos) then
        call sum_mn_name(f(l1:l2,m,n,ispecialvar)**2,idiag_psirms,lsqrt=.true.)
        call sum_mn_name(f(l1:l2,m,n,ispecialvar2)**2,idiag_dpsirms,lsqrt=.true.)
      endif
!
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=oscillation_3D_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=oscillation_3D_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=oscillation_3D_run_pars, IOSTAT=iostat)
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
      use Diagnostics
      use FArrayManager, only: farray_index_append
      use Sub
!
!   SAMPLE IMPLEMENTATION
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_psirms=0; idiag_dpsirms=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'psirms',idiag_psirms)
        call parse_name(iname,cname(iname),cform(iname),'dpsirms',idiag_dpsirms)
      enddo
!
!  write column where which variable is stored
!
!     if (lwr) then
!       call farray_index_append('i_psirms',idiag_psirms)
!       call farray_index_append('i_dpsirms',idiag_dpsirms)
!     endif
!
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
!  Calculate an additional 'special' term on the right hand side of the
!  energy equation.
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
!!  df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
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
    subroutine special_calc_spectra(f,spectrum,spectrumhel, &
      spectrum_2d,spectrum_2d_hel,&
      lfirstcall,kind)

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:) :: spectrum,spectrumhel
      real, dimension (:,:) :: spectrum_2d,spectrum_2d_hel
      logical :: lfirstcall
      character(LEN=3) :: kind

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(spectrum,spectrumhel)
      call keep_compiler_quiet(spectrum_2d,spectrum_2d_hel)
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
    subroutine input_persist_special_id(id,done)
!
      integer :: id
      logical :: done

      call keep_compiler_quiet(id)
      call keep_compiler_quiet(done)

    endsubroutine input_persist_special_id
!***********************************************************************
    subroutine input_persist_special
            
    endsubroutine input_persist_special
!***********************************************************************
    logical function output_persistent_special()

      output_persistent_special=.false.

    endfunction output_persistent_special
!***********************************************************************
endmodule oscillation_3D
