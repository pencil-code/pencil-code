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
! PENCILS PROVIDED bb(3); bij(3,3); jxbr(3); ss12; b2; uxb(3)
! PENCILS PROVIDED diva; del2a(3)
!
!***************************************************************

module Magnetic

  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet

  implicit none

  include 'magnetic.h'

  real, dimension (ninit) :: amplaa=0.0,kx_aa=1.,ky_aa=1.,kz_aa=1.
  real :: kx=1.,ky=1.,kz=1.,ABC_A=1.,ABC_B=1.,ABC_C=1.
  real :: brms=0., bmz_beltrami_phase=0.
  real, dimension(nx) :: meanfield_EMFdotB
  real :: inertial_length=0.,linertial_2
  logical :: lelectron_inertia=.false.

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_b2m=0,idiag_bm2=0,idiag_j2m=0,idiag_jm2=0,idiag_abm=0
  integer :: idiag_jbm=0,idiag_epsM=0,idiag_vArms=0,idiag_vAmax=0
  integer :: idiag_brms=0,idiag_bmax=0,idiag_jrms=0,idiag_jmax=0
  integer :: idiag_bx2m=0, idiag_by2m=0, idiag_bz2m=0,idiag_bmz=0
  integer :: idiag_bxmz=0,idiag_bymz=0,idiag_bzmz=0,idiag_bmx=0,idiag_bmy=0
  integer :: idiag_bxmxy=0,idiag_bymxy=0,idiag_bzmxy=0
  integer :: idiag_uxbm=0,idiag_oxuxbm=0,idiag_jxbxbm=0,idiag_uxDxuxbm=0
  integer :: idiag_b2mphi=0
  integer :: idiag_bmxy_rms=0
  integer :: idiag_bsinphz=0
  integer :: idiag_bcosphz=0

  contains

!***********************************************************************
    subroutine register_magnetic()
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaa, etc; increase nvar accordingly
!  3-may-2002/wolf: dummy routine
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_aa called twice')
      first = .false.
!
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id$")
!
    endsubroutine register_magnetic
!***********************************************************************
    subroutine initialize_magnetic(f,lstarting)
!
!  Perform any post-parameter-read initialization
!
!  24-nov-2002/tony: dummy routine
      use Cdata
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
    subroutine init_aa(f,xx,yy,zz)
!
!  initialise magnetic field; called from start.f90
!  3-may-2002/wolf: dummy routine
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(xx,yy,zz)
!
    endsubroutine init_aa
!***********************************************************************
    subroutine pert_aa(f)
!
!   perturb magnetic field when reading old NON-magnetic snapshot
!   called from run.f90
!   30-july-2004/dave: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine pert_aa
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
      if (lpencil(i_bb)) p%bb=0.
      if (lpencil(i_b2)) p%b2=0.
      if (lpencil(i_jxbr)) p%jxbr=0.
      if (lpencil(i_bij)) p%bij=0.
      if (lpencil(i_uxb)) p%uxb=0.
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_magnetic
!***********************************************************************
    subroutine daa_dt(f,df,p)
!
!  magnetic field evolution
!  3-may-2002/wolf: dummy routine
!
      use Cdata
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
!   1-jul-08/axel: dummy
!
      use Cdata
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
    subroutine df_diagnos_magnetic(f,df,p)
!
!  11-oct-06/axel: dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)  :: f, df, p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine df_diagnos_magnetic
!***********************************************************************
    subroutine rescaling_magnetic(f)
!
!  Dummy routine
!
!  22-feb-05/axel: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine rescaling_magnetic
!***********************************************************************
    subroutine read_magnetic_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) call keep_compiler_quiet(iostat)
      call keep_compiler_quiet(unit)

    endsubroutine read_magnetic_init_pars
!***********************************************************************
    subroutine write_magnetic_init_pars(unit)
      integer, intent(in) :: unit

      call keep_compiler_quiet(unit)

    endsubroutine write_magnetic_init_pars
!***********************************************************************
    subroutine read_magnetic_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) call keep_compiler_quiet(iostat)
      call keep_compiler_quiet(unit)

    endsubroutine read_magnetic_run_pars
!***********************************************************************
    subroutine write_magnetic_run_pars(unit)
      integer, intent(in) :: unit

      call keep_compiler_quiet(unit)

    endsubroutine write_magnetic_run_pars
!***********************************************************************
    subroutine rprint_magnetic(lreset,lwrite)
!
!  reads and registers print parameters relevant for magnetic fields
!  dummy routine
!
!   3-may-02/axel: coded
!
      use Cdata
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!  idl needs this even if everything is zero
!
      if (lwr) then
        write(3,*) 'i_abm=',idiag_abm
        write(3,*) 'i_jbm=',idiag_jbm
        write(3,*) 'i_b2m=',idiag_b2m
        write(3,*) 'i_bm2=',idiag_bm2
        write(3,*) 'i_j2m=',idiag_j2m
        write(3,*) 'i_jm2=',idiag_jm2
        write(3,*) 'i_epsM=',idiag_epsM
        write(3,*) 'i_brms=',idiag_brms
        write(3,*) 'i_bmax=',idiag_bmax
        write(3,*) 'i_jrms=',idiag_jrms
        write(3,*) 'i_jmax=',idiag_jmax
        write(3,*) 'i_vArms=',idiag_vArms
        write(3,*) 'i_vAmax=',idiag_vAmax
        write(3,*) 'i_bx2m=',idiag_bx2m
        write(3,*) 'i_by2m=',idiag_by2m
        write(3,*) 'i_bz2m=',idiag_bz2m
        write(3,*) 'i_uxbm=',idiag_uxbm
        write(3,*) 'i_oxuxbm=',idiag_oxuxbm
        write(3,*) 'i_jxbxbm=',idiag_jxbxbm
        write(3,*) 'i_uxDxuxbm=',idiag_uxDxuxbm
        write(3,*) 'i_bxmz=',idiag_bxmz
        write(3,*) 'i_bymz=',idiag_bymz
        write(3,*) 'i_bzmz=',idiag_bzmz
        write(3,*) 'i_bmx=',idiag_bmx
        write(3,*) 'i_bmy=',idiag_bmy
        write(3,*) 'i_bmz=',idiag_bmz
        write(3,*) 'i_bxmxy=',idiag_bxmxy
        write(3,*) 'i_bymxy=',idiag_bymxy
        write(3,*) 'i_bzmxy=',idiag_bzmxy
        write(3,*) 'i_b2mphi=',idiag_b2mphi
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
!  Write slices for animation of magnetic variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices)
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
!  18-aug-2006/tobi: coded
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,3), intent (out) :: bb_hat

      if (NO_WARN) then
        print *,f
        bb_hat = 0.
      endif

    endsubroutine bb_unitvec_shock
!***********************************************************************
    subroutine input_persistent_magnetic(id,lun,done)
! 
!  Dummy routine
! 
      use Cdata, only: lroot
!
      integer :: id,lun
      logical :: done
!
      if (NO_WARN) print*,id,lun,done
!
    endsubroutine input_persistent_magnetic
!***********************************************************************
    subroutine output_persistent_magnetic(lun)
!
!  Dummy routine
!
      integer :: lun
!
      if (NO_WARN) print*,lun
!
    endsubroutine output_persistent_magnetic
!***********************************************************************
    subroutine remove_mean_emf(f,df)
!
!  Dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,my,mz,mvar), intent (inout) :: df

    endsubroutine remove_mean_emf
!***********************************************************************
endmodule Magnetic
