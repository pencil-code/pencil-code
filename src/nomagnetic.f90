! $Id$

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lmagnetic = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED bb(3); bij(3,3); jxbr(3); ss12; b2; uxb(3); jj(3)
! PENCILS PROVIDED aa(3) ; diva; del2a(3); aij(3,3), bunit(3)
!
!***************************************************************
module Magnetic
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id
!
  implicit none
!
  include 'magnetic.h'
!
  real, dimension(3) :: B_ext_inv=(/0.0,0.0,0.0/)
  real, dimension (ninit) :: amplaa=0.0,kx_aa=1.,ky_aa=1.,kz_aa=1.
  real :: kx=1.,ky=1.,kz=1.
  real :: brms=0., bmz_beltrami_phase=0.
  real, dimension (mz,3) :: aamz
  real, dimension (nz,3) :: bbmz,jjmz
  real :: inertial_length=0.,linertial_2
  logical :: lelectron_inertia=.false.
  logical :: lcalc_aamean=.false.
!
  integer :: idiag_b2m=0,idiag_bm2=0,idiag_j2m=0,idiag_jm2=0,idiag_abm=0
  integer :: idiag_jbm=0,idiag_epsM=0,idiag_vArms=0,idiag_vAmax=0
  integer :: idiag_brms=0,idiag_bmax=0,idiag_jrms=0,idiag_jmax=0
  integer :: idiag_bx2m=0, idiag_by2m=0, idiag_bz2m=0,idiag_bmz=0
  integer :: idiag_axmz=0,idiag_aymz=0
  integer :: idiag_bxmz=0,idiag_bymz=0,idiag_bzmz=0,idiag_bmx=0,idiag_bmy=0
  integer :: idiag_bxmxy=0,idiag_bymxy=0,idiag_bzmxy=0
  integer :: idiag_uxbm=0,idiag_oxuxbm=0,idiag_jxbxbm=0,idiag_uxDxuxbm=0
  integer :: idiag_b2mphi=0
  integer :: idiag_bmxy_rms=0
  integer :: idiag_bsinphz=0
  integer :: idiag_bcosphz=0
  integer :: idiag_magfricmax=0
!
  contains
!***********************************************************************
    subroutine register_magnetic()
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaa, etc; increase nvar accordingly
!  3-may-2002/wolf: dummy routine
!
!  identify version number
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_magnetic
!***********************************************************************
    subroutine initialize_magnetic(f,lstarting)
!
!  Perform any post-parameter-read initialization
!
!  24-nov-2002/tony: dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
!  Precalculate 1/mu (moved here from register.f90)
!
      mu01=1./mu0
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
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
    subroutine pencil_criteria_magnetic()
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
    subroutine calc_pencils_magnetic(f,p)
!
!  Calculate Magnetic pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in)  :: f
      intent(inout) :: p
!
      if (lpencil(i_aa)) p%aa=0.0
      if (lpencil(i_bb)) p%bb=0.0
      if (lpencil(i_bunit)) p%bunit=0.0
      if (lpencil(i_b2)) p%b2=0.0
      if (lpencil(i_jxbr)) p%jxbr=0.0
      if (lpencil(i_bij)) p%bij=0.0
      if (lpencil(i_uxb)) p%uxb=0.0
      if (lpencil(i_jj)) p%jj=0.0
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_magnetic
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
    subroutine calc_lmagnetic_pars(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lmagnetic_pars
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
    subroutine read_magnetic_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      call keep_compiler_quiet(present(iostat))
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
    subroutine read_magnetic_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      call keep_compiler_quiet(present(iostat))
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
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!
      if (lwr) then
        write(3,*) 'nname=',nname
        write(3,*) 'nnamexy=',nnamexy
        write(3,*) 'nnamez=',nnamez
        write(3,*) 'iaa=',iaa
        write(3,*) 'iax=',iax
        write(3,*) 'iay=',iay
        write(3,*) 'iaz=',iaz
        write(3,*) 'idiag_bmxy_rms=',idiag_bmxy_rms
      endif
!
      call keep_compiler_quiet(lreset)

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
      integer :: id
      logical :: done
!
      call keep_compiler_quiet(id)
      call keep_compiler_quiet(done)
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
    subroutine dynamical_resistivity(umax)
!
!  dummy
!
      real, intent(in) :: umax
!     
      call keep_compiler_quiet(umax)
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
    subroutine expand_shands_magnetic()
!
!  Dummy
!
    endsubroutine expand_shands_magnetic
!***********************************************************************
endmodule Magnetic
