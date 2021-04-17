! $Id$

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lmagnetic = .false.
! CPARAM logical, parameter :: lbfield = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED bb(3); bbb(3); bij(3,3); jxbr(3); ss12; b2; uxb(3); jj(3)
! PENCILS PROVIDED aa(3) ; diva; del2a(3); aij(3,3), bunit(3); va2
!
!***************************************************************
module Magnetic
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include 'magnetic.h'
!
  real, dimension(3) :: B_ext_inv=(/0.0,0.0,0.0/)
  real, dimension (mz,3) :: aamz
  real, dimension (nz,3) :: bbmz,jjmz
  real :: inertial_length=0.,linertial_2
  logical :: lelectron_inertia=.false.
  logical :: lcalc_aameanz=.false., lcalc_aamean=.false.
  logical, dimension(7) :: lresi_dep=.false. 
  logical :: lcovariant_magnetic=.false.
  integer :: pushpars2c, pushdiags2c  ! should be procedure pointer (F2003)
!
  integer :: idiag_axmz=0,idiag_aymz=0
  integer :: idiag_bxmz=0,idiag_bymz=0
  integer :: idiag_bsinphz=0,idiag_bcosphz=0
!
  contains
!***********************************************************************
    subroutine register_magnetic
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaa, etc; increase nvar accordingly
!
!  3-may-2002/wolf: dummy routine
!
!  identify version number
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_magnetic
!***********************************************************************
    subroutine initialize_magnetic(f)
!
!  Perform any post-parameter-read initialization
!
!  24-nov-2002/tony: dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Precalculate 1/mu (moved here from register.f90)
!
      mu01=1./mu0
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_magnetic
!***********************************************************************
    subroutine init_aa(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_aa
!***********************************************************************
    subroutine pencil_criteria_magnetic
!
!  All pencils that the Magnetic module depends on are specified here.
!
!  20-11-04/anders: coded
!
    endsubroutine pencil_criteria_magnetic
!***********************************************************************
    subroutine pencil_interdep_magnetic(lpencil_in)
!
!  Interdependency among pencils provided by the Magnetic module
!  is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_magnetic
!***********************************************************************
    subroutine magnetic_before_boundary(f)
!
!  Conduct pre-processing required before boundary conditions and pencil
!  calculations.
!
!  29-may-14/ccyang: dummy
!
      real, dimension(mx,my,mz,mfarray), intent(in):: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine magnetic_before_boundary
!***********************************************************************
    subroutine calc_pencils_magnetic_std(f,p)
!
!  Standard version (_std): global variable lpencil contains information about needed pencils.
!
      real, dimension (mx,my,mz,mfarray), intent(inout):: f
      type (pencil_case),                 intent(out)  :: p
!
      call calc_pencils_magnetic_pencpar(f,p,lpencil)
!
    endsubroutine calc_pencils_magnetic_std
!***********************************************************************
    subroutine calc_pencils_magnetic_pencpar(f,p,lpenc_loc)
!
!  Calculate Magnetic pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      logical, dimension(:) :: lpenc_loc
!
      intent(in)  :: f, lpenc_loc
      intent(inout) :: p
!
      if (lpenc_loc(i_aa)) p%aa=0.0
      if (lpenc_loc(i_bb)) p%bb=0.0
      if (lpenc_loc(i_bbb)) p%bbb=0.0
      if (lpenc_loc(i_bunit)) p%bunit=0.0
      if (lpenc_loc(i_b2)) p%b2=0.0
      if (lpenc_loc(i_jxbr)) p%jxbr=0.0
      if (lpenc_loc(i_bij)) p%bij=0.0
      if (lpenc_loc(i_uxb)) p%uxb=0.0
      if (lpenc_loc(i_jj)) p%jj=0.0
      if (lpenc_loc(i_va2)) p%va2=0.0
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_magnetic_pencpar
!***********************************************************************
    subroutine update_char_vel_magnetic(f)
!
! Dummy 
!
!  25-sep-15/MR+joern: coded
!
      real, dimension(mx,my,mz,mfarray), intent(INOUT) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine update_char_vel_magnetic
!***********************************************************************
    subroutine daa_dt(f,df,p)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f, df, p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine daa_dt
!***********************************************************************
    subroutine calc_diagnostics_magnetic(f,p)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f, p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_diagnostics_magnetic
!***********************************************************************
    subroutine time_integrals_magnetic(f,p)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine time_integrals_magnetic
!***********************************************************************
    subroutine df_diagnos_magnetic(df,p)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)  ::  df, p
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine df_diagnos_magnetic
!***********************************************************************
    subroutine magnetic_after_boundary(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine magnetic_after_boundary
!***********************************************************************
    subroutine magnetic_calc_spectra(f,spectrum,spectrum_hel,lfirstcall,kind)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:) :: spectrum,spectrum_hel
      logical :: lfirstcall
      character(LEN=3) :: kind
!
      call fatal_error("magnetic_calc_spectra","impossible: iaakim=0, ieekim=0")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(spectrum)
      call keep_compiler_quiet(spectrum_hel)
      call keep_compiler_quiet(lfirstcall)
      call keep_compiler_quiet(kind)
!
    endsubroutine magnetic_calc_spectra
!***********************************************************************
    subroutine rescaling_magnetic(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine rescaling_magnetic
!***********************************************************************
    subroutine read_magnetic_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_magnetic_init_pars
!***********************************************************************
    subroutine write_magnetic_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_magnetic_init_pars
!***********************************************************************
    subroutine read_magnetic_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_magnetic_run_pars
!***********************************************************************
    subroutine write_magnetic_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_magnetic_run_pars
!***********************************************************************
    subroutine rprint_magnetic(lreset,lwrite)
!
!  reads and registers print parameters relevant for magnetic fields
!  dummy routine
!
!   3-may-02/axel: coded
!  26-aug-13/MR: unneeded output of idiag* removed
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)

    endsubroutine rprint_magnetic
!***********************************************************************
    subroutine get_slices_magnetic(f,slices)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_magnetic
!***********************************************************************
    subroutine bdry_magnetic(f,quench,task)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (nx) :: quench
      character (len=*), intent(in) :: task

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(quench)
      call keep_compiler_quiet(task)
      call fatal_error('bdry_magnetic','not to be called w/o B-field')
!
    endsubroutine bdry_magnetic
!***********************************************************************
    subroutine calc_mfield
!
!  Dummy routine
!
    endsubroutine calc_mfield
!***********************************************************************
    subroutine bb_unitvec_shock(f,bb_hat)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,3), intent (out) :: bb_hat
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(bb_hat)
!
    endsubroutine bb_unitvec_shock
!***********************************************************************
    subroutine input_persistent_magnetic(id,done)
!
!  Dummy routine
!
      integer, optional :: id
      logical, optional :: done
!
      if (present (id)) call keep_compiler_quiet(id)
      if (present (done)) call keep_compiler_quiet(done)
!
    endsubroutine input_persistent_magnetic
!***********************************************************************
    logical function output_persistent_magnetic()
!
!  Dummy routine
!
      output_persistent_magnetic = .false.
!
    endfunction output_persistent_magnetic
!***********************************************************************
    subroutine dynamical_resistivity(uc)
!
!  dummy
!
      real, intent(in) :: uc
!
      call keep_compiler_quiet(uc)
!
    endsubroutine dynamical_resistivity
!***********************************************************************
    subroutine split_update_magnetic(f)
!
!  dummy
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine split_update_magnetic
!***********************************************************************
    subroutine magnetic_after_timestep(f,df,dtsub)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dtsub
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dtsub)
!
    endsubroutine magnetic_after_timestep
!***********************************************************************
    subroutine expand_shands_magnetic
!
!  Dummy
!
    endsubroutine expand_shands_magnetic
!***********************************************************************
    subroutine get_bext(B_ext_out)
!
!  Dummy
!
      real, dimension(3), intent(out) :: B_ext_out
!
      B_ext_out = 0.0
!
    endsubroutine get_bext
!***********************************************************************
endmodule Magnetic
