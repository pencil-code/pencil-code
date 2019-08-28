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
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MAUX CONTRIBUTION 18
!
! PENCILS PROVIDED stress_ij(6)
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
module Special
!
  use Cparam
  use Cdata
  use Initcond
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
! Declare index of new variables in f array (if any).
!
  character (len=labellen) :: inithij='nothing'
  character (len=labellen) :: initgij='nothing'
  character (len=labellen) :: ctrace_factor='1/3'
  character (len=labellen) :: cstress_prefactor='6'
  character (len=labellen) :: fourthird_in_stress='4/3'
  character (len=labellen) :: cc_light='1'
  character (len=labellen) :: aux_stress='stress'
  real :: amplhij=0., amplgij=0., dummy=0.
  real :: trace_factor=0., stress_prefactor, fourthird_factor, EGWpref
  real :: nscale_factor_conformal=1., tshift=0.
  logical :: lno_transverse_part=.false.
  logical :: lswitch_sign_e_X=.true., ldebug_print=.false., lkinGW=.true.
  logical :: lStress_as_aux=.true., lreynolds=.false.
  logical :: lggTX_as_aux=.true., lhhTX_as_aux=.true.
  logical :: lremove_mean_hij=.false., lremove_mean_gij=.false.
  real, dimension(3,3) :: ij_table
  real :: c_light2=1.
!
!  Do this here because shared variables for this array doesn't work on Beskow.
!
  integer, parameter :: nk=nxgrid/2
  real, dimension(nk) :: specSCL, specVCT, specTpq
  real, dimension(nk) :: specGWs   ,specGWh   ,specGWm   ,specStr
  real, dimension(nk) :: specGWshel,specGWhhel,specGWmhel,specStrhel
  public :: specGWs, specGWshel, specGWh, specGWhhel, specGWm, specGWmhel
  public :: specStr, specStrhel, specSCL, specVCT, specTpq
!
! input parameters
  namelist /special_init_pars/ &
    ctrace_factor, cstress_prefactor, fourthird_in_stress, lno_transverse_part, &
    inithij, initgij, amplhij, amplgij, lStress_as_aux, &
    lggTX_as_aux, lhhTX_as_aux
!
! run parameters
  namelist /special_run_pars/ &
    ctrace_factor, cstress_prefactor, fourthird_in_stress, lno_transverse_part, &
    ldebug_print, lswitch_sign_e_X, &
    nscale_factor_conformal, tshift, cc_light, &
    lStress_as_aux, lkinGW, aux_stress, &
    lggTX_as_aux, lhhTX_as_aux, lremove_mean_hij, lremove_mean_gij
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_g11pt=0       ! DIAG_DOC: $g_{11}(x_1,y_1,z_1,t)$
  integer :: idiag_g22pt=0       ! DIAG_DOC: $g_{22}(x_1,y_1,z_1,t)$
  integer :: idiag_g33pt=0       ! DIAG_DOC: $g_{33}(x_1,y_1,z_1,t)$
  integer :: idiag_g12pt=0       ! DIAG_DOC: $g_{12}(x_1,y_1,z_1,t)$
  integer :: idiag_g23pt=0       ! DIAG_DOC: $g_{23}(x_1,y_1,z_1,t)$
  integer :: idiag_g31pt=0       ! DIAG_DOC: $g_{31}(x_1,y_1,z_1,t)$
  integer :: idiag_hhTpt=0       ! DIAG_DOC: $h_{T}(x_1,y_1,z_1,t)$
  integer :: idiag_hhXpt=0       ! DIAG_DOC: $h_{X}(x_1,y_1,z_1,t)$
  integer :: idiag_ggTpt=0       ! DIAG_DOC: $\dot{h}_{T}(x_1,y_1,z_1,t)$
  integer :: idiag_ggXpt=0       ! DIAG_DOC: $\dot{h}_{X}(x_1,y_1,z_1,t)$
  integer :: idiag_hhTp2=0       ! DIAG_DOC: $h_{T}(x_1,y_1,z_1,t)$
  integer :: idiag_hhXp2=0       ! DIAG_DOC: $h_{X}(x_1,y_1,z_1,t)$
  integer :: idiag_ggTp2=0       ! DIAG_DOC: $\dot{h}_{T}(x_1,y_1,z_1,t)$
  integer :: idiag_ggXp2=0       ! DIAG_DOC: $\dot{h}_{X}(x_1,y_1,z_1,t)$
  integer :: idiag_hrms=0        ! DIAG_DOC: $\bra{h_T^2+h_X^2}^{1/2}$
  integer :: idiag_EEGW=0        ! DIAG_DOC: $\bra{g_T^2+g_X^2}\,c^2/(32\pi G)$
  integer :: idiag_gg2m=0        ! DIAG_DOC: $\bra{g_T^2+g_X^2}$
  integer :: idiag_hhT2m=0       ! DIAG_DOC: $\bra{h_T^2}$
  integer :: idiag_hhX2m=0       ! DIAG_DOC: $\bra{h_X^2}$
  integer :: idiag_hhTXm=0       ! DIAG_DOC: $\bra{h_T h_X}$
  integer :: idiag_ggT2m=0       ! DIAG_DOC: $\bra{g_T^2}$
  integer :: idiag_ggX2m=0       ! DIAG_DOC: $\bra{g_X^2}$
  integer :: idiag_ggTXm=0       ! DIAG_DOC: $\bra{g_T g_X}$
!
  contains
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
!   3-aug-17/axel: coded
!
      use Sub, only: register_report_aux
      use FArrayManager
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Register ggT and ggX as auxiliary arrays
!  May want to do this only when Fourier transform is enabled.
!
      if (lggTX_as_aux) then
        call farray_register_auxiliary('ggT',iggT)
        call farray_register_auxiliary('ggX',iggX)
        call farray_register_auxiliary('ggTim',iggTim)
        call farray_register_auxiliary('ggXim',iggXim)
      endif
!
      if (lhhTX_as_aux) then
        call farray_register_auxiliary('hhT',ihhT)
        call farray_register_auxiliary('hhX',ihhX)
        call farray_register_auxiliary('hhTim',ihhTim)
        call farray_register_auxiliary('hhXim',ihhXim)
      endif
!
      if (lStress_as_aux) then
        call farray_register_auxiliary('StT',iStressT)
        call farray_register_auxiliary('StX',iStressX)
        call farray_register_auxiliary('StTim',iStressTim)
        call farray_register_auxiliary('StXim',iStressXim)
        call farray_register_auxiliary('Str',iStress_ij,vector=6)
      endif
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      use EquationOfState, only: cs0
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  set index table
!
      ij_table(1,1)=1
      ij_table(2,2)=2
      ij_table(3,3)=3
      ij_table(1,2)=4
      ij_table(2,3)=5
      ij_table(3,1)=6
      ij_table(2,1)=4
      ij_table(3,2)=5
      ij_table(1,3)=6
!
!  determine trace factor
!
      select case (ctrace_factor)
        case ('0'); trace_factor=.0
        case ('1/2'); trace_factor=.5
        case ('1/3'); trace_factor=onethird
        case default
          call fatal_error("initialize_special: No such value for ctrace_factor:" &
              ,trim(ctrace_factor))
      endselect
!
!  Determine fourthird_in_stress. This factor is normally 4/3, but it can be
!  set to unity in case we want to pretend that the kinematic Beltrami field
!  has the same prefactor as a magnetic one.
!
      select case (fourthird_in_stress)
        case ('1'); fourthird_factor=1.
        case ('4/3'); fourthird_factor=fourthird
        case default
          call fatal_error("initialize_special: No such value for fourthird_in_stress:" &
              ,trim(fourthird_in_stress))
      endselect
!
!  determine stress_prefactor and GW energy prefactor,
!  which is EGWpref=.5*16.*pi/stress_prefactor**2
!  At the moment, only the case stress_prefactor=6 is to be used.
!  The other cases are kept for backward compatibility.
!
      select case (cstress_prefactor)
        case ('1'); stress_prefactor=1.; EGWpref=8.*pi
        case ('6'); stress_prefactor=6.; EGWpref=1./6.
        case ('6old'); stress_prefactor=6.; EGWpref=1./(32.*pi)
        case ('16pi'); stress_prefactor=16.*pi; EGWpref=1./(32.*pi)
        case ('16pi_corr'); stress_prefactor=16.*pi; EGWpref=1./(16.*pi)
        case ('16piG/c^2'); stress_prefactor=16.*pi*G_Newton_cgs/c_light_cgs**2;
          EGWpref=c_light_cgs**2/(32.*pi*G_Newton_cgs)
        case default
          call fatal_error("initialize_special: No such value for ctrace_factor:" &
              ,trim(ctrace_factor))
      endselect
      if (headt) print*,'stress_prefactor=',stress_prefactor
      if (headt) print*,'EGWpref=',EGWpref
!
!  set speed of light
!
      select case (cc_light)
        case ('1'); c_light2=1.
        case ('cgs'); c_light2=c_light_cgs**2
        case default
          call fatal_error("initialize_special: No such value for cc_light:" &
              ,trim(ctrace_factor))
      endselect
      if (headt) print*,'c_light2=',c_light2
!
!  Determine whether Reynolds stress needs to be computed:
!  Compute Reynolds stress when lhydro or lhydro_kinematic,
!  provided lkinGW=T
!
      lreynolds=(lhydro.or.lhydro_kinematic).and.lkinGW
!
!  give a warning if cs0**2=1
!
!     if (cs0==1.) call fatal_error('gravitational_waves_hij6', &
!         'cs0 should probably not be unity')
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
!
!  initial condition for hij
!
      select case (inithij)
        case ('nothing')
          if (lroot) print*,'init_special: nothing'
          f(:,:,:,ihhT:ihhXim)=0.
        case default
          call fatal_error("init_special: No such value for inithij:" &
              ,trim(inithij))
      endselect
!
!  initial condition for gij
!
      select case (initgij)
        case ('nothing')
          if (lroot) print*,'init_special: nothing'
          f(:,:,:,iggT:iggXim)=0.
        case default
          call fatal_error("init_special: No such value for initgij:" &
              ,trim(initgij))
      endselect
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!   1-apr-06/axel: coded
!
!  Velocity field needed for Reynolds stress
!
      if (lreynolds) then
        lpenc_requested(i_uu)=.true.
        lpenc_requested(i_rho)=.true.
        if (trace_factor/=0.) lpenc_requested(i_u2)=.true.
      endif
!
!  Magnetic field needed for Maxwell stress
!
      if (lmagnetic) then
        lpenc_requested(i_bb)=.true.
        if (trace_factor/=0.) lpenc_requested(i_b2)=.true.
      endif
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-jul-06/tony: coded
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
!  24-aug-17/axel: coded
!   7-jun-18/axel: included 4/3*rho factor
!
      use Deriv, only: derij
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: tmp
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
      integer :: i, j, ij
!
!  Construct stress tensor; notice opposite signs for u and b.
!
      p%stress_ij=0.0
      do j=1,3
      do i=1,j
        ij=ij_table(i,j)
        if (lreynolds) p%stress_ij(:,ij)=p%stress_ij(:,ij)+p%uu(:,i)*p%uu(:,j)*fourthird_factor*p%rho
        if (lmagnetic) p%stress_ij(:,ij)=p%stress_ij(:,ij)-p%bb(:,i)*p%bb(:,j)
!
!  Remove trace.
!
        if (i==j) then
          if (lreynolds) p%stress_ij(:,ij)=p%stress_ij(:,ij)-trace_factor*p%u2*fourthird_factor*p%rho
          if (lmagnetic) p%stress_ij(:,ij)=p%stress_ij(:,ij)+trace_factor*p%b2
        endif
      enddo
      enddo
!
    endsubroutine calc_pencils_special
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
!  This routine computes various diagnostic quantities, but those
!  would be more easily done further below when sign_switch is known.
!  But this only affects the helicities. The values of EEGW and hrms are
!  however correct and agree with those of gravitational_waves_hij6.
!
!  06-oct-03/tony: coded
!  07-feb-18/axel: added nscale_factor=0 (no expansion), =.5 (radiation era)
!
      use Diagnostics
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: scale_factor, stress_prefactor2, sign_switch=0
      type (pencil_case) :: p
!
      integer :: ij
!
      intent(in) :: p
      intent(inout) :: f,df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
!  Compute scale factor.
!  Note: to prevent division by zero, it is best to put tstart=1. in start.in.
!  If that is not done, one can put here tshift=1., for example.
!  If that is not the case either, we put scale_factor=1.
!  At the next timestep, this will no longer be a problem.
!
      if (t+tshift==0.) then
        scale_factor=1.
      else
        scale_factor=(t+tshift)**nscale_factor_conformal
      endif
      stress_prefactor2=stress_prefactor/scale_factor
!
!  Assemble rhs of GW equations.
!
      do ij=1,6
        f(l1:l2,m,n,iStress_ij+ij-1)=stress_prefactor2*p%stress_ij(:,ij)
      enddo
!
!  diagnostics
!
       if (ldiagnos) then
         if (lggTX_as_aux) then
           if (idiag_EEGW/=0) call sum_mn_name((f(l1:l2,m,n,iggT)**2+f(l1:l2,m,n,iggTim)**2 &
                                               +f(l1:l2,m,n,iggX)**2+f(l1:l2,m,n,iggXim)**2 &
                                               )*nwgrid*EGWpref,idiag_EEGW)
           if (idiag_gg2m/=0) call sum_mn_name((f(l1:l2,m,n,iggT)**2+f(l1:l2,m,n,iggTim)**2 &
                                               +f(l1:l2,m,n,iggX)**2+f(l1:l2,m,n,iggXim)**2 &
                                               )*nwgrid,idiag_gg2m)
           if (idiag_ggT2m/=0) call sum_mn_name((f(l1:l2,m,n,iggT)**2+f(l1:l2,m,n,iggTim)**2 &
                                                )*nwgrid,idiag_ggT2m)
           if (idiag_ggX2m/=0) call sum_mn_name((f(l1:l2,m,n,iggX)**2+f(l1:l2,m,n,iggXim)**2 &
                                                )*nwgrid,idiag_ggX2m)
           if (idiag_ggTXm/=0) call sum_mn_name((f(l1:l2,m,n,iggT  )*f(l1:l2,m,n,iggXim) &
                                                -f(l1:l2,m,n,iggTim)*f(l1:l2,m,n,iggX  ) &
                                                )*nwgrid*sign_switch,idiag_ggTXm)
         endif
         if (lhhTX_as_aux) then
           if (idiag_hrms/=0) call sum_mn_name((f(l1:l2,m,n,ihhT)**2+f(l1:l2,m,n,ihhTim)**2 &
                                               +f(l1:l2,m,n,ihhX)**2+f(l1:l2,m,n,ihhXim)**2 &
                                               )*nwgrid,idiag_hrms,lsqrt=.true.)
           if (idiag_hhT2m/=0) call sum_mn_name((f(l1:l2,m,n,ihhT)**2+f(l1:l2,m,n,ihhTim)**2 &
                                                )*nwgrid,idiag_hhT2m)
           if (idiag_hhX2m/=0) call sum_mn_name((f(l1:l2,m,n,ihhX)**2+f(l1:l2,m,n,ihhXim)**2 &
                                                )*nwgrid,idiag_hhX2m)
           if (idiag_hhTXm/=0) call sum_mn_name((f(l1:l2,m,n,ihhT  )*f(l1:l2,m,n,ihhXim) &
                                                -f(l1:l2,m,n,ihhTim)*f(l1:l2,m,n,ihhX  ) &
                                                )*nwgrid*sign_switch,idiag_hhTXm)
         endif
!
         if (lroot.and.m==mpoint.and.n==npoint) then
           if (idiag_g11pt/=0) call save_name(f(lpoint,m,n,igij+1-1),idiag_g11pt)
           if (idiag_g22pt/=0) call save_name(f(lpoint,m,n,igij+2-1),idiag_g22pt)
           if (idiag_g33pt/=0) call save_name(f(lpoint,m,n,igij+3-1),idiag_g33pt)
           if (idiag_g12pt/=0) call save_name(f(lpoint,m,n,igij+4-1),idiag_g12pt)
           if (idiag_g23pt/=0) call save_name(f(lpoint,m,n,igij+5-1),idiag_g23pt)
           if (idiag_g31pt/=0) call save_name(f(lpoint,m,n,igij+6-1),idiag_g31pt)
           if (lhhTX_as_aux) then
             if (idiag_hhTpt/=0) call save_name(f(lpoint,m,n,ihhT),idiag_hhTpt)
             if (idiag_hhXpt/=0) call save_name(f(lpoint,m,n,ihhX),idiag_hhXpt)
           endif
           if (lggTX_as_aux) then
             if (idiag_ggTpt/=0) call save_name(f(lpoint,m,n,iggT),idiag_ggTpt)
             if (idiag_ggXpt/=0) call save_name(f(lpoint,m,n,iggX),idiag_ggXpt)
           endif
         endif
!
         if (lroot.and.m==mpoint2.and.n==npoint2) then
           if (lhhTX_as_aux) then
             if (idiag_hhTp2/=0) call save_name(f(lpoint2,m,n,ihhT),idiag_hhTp2)
             if (idiag_hhXp2/=0) call save_name(f(lpoint2,m,n,ihhX),idiag_hhXp2)
           endif
           if (lggTX_as_aux) then
             if (idiag_ggTp2/=0) call save_name(f(lpoint2,m,n,iggT),idiag_ggTp2)
             if (idiag_ggXp2/=0) call save_name(f(lpoint2,m,n,iggX),idiag_ggXp2)
           endif
         endif
       endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine special_before_boundary(f)
!
!  Possibility to modify the f array before the boundaries are
!  communicated.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  13-may-18/axel: added remove_mean_value for hij and gij
!
      use Sub, only: remove_mean_value
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  07-aug-17/axel: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    endsubroutine special_after_boundary
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
!  Compute the transverse part of the stress tensor by going into Fourier space.
!
      if (llast) call compute_gT_and_gX_from_gij(f,'St')
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dt_)
!
    endsubroutine special_after_timestep
!***********************************************************************
    subroutine compute_gT_and_gX_from_gij(f,label)
!
!  Compute the transverse part of the stress tensor by going into Fourier space.
!
!  07-aug-17/axel: coded
!
      use Fourier, only: fourier_transform, fft_xyz_parallel
      use SharedVariables, only: put_shared_variable
!
      real, dimension (:,:,:,:), allocatable :: Tpq_re, Tpq_im
      real, dimension (:,:,:), allocatable :: S_T_re, S_T_im, S_X_re, S_X_im
      real, dimension (:), allocatable :: kx, ky, kz
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (6) :: Pij=0., e_T, e_X, Sij_re, Sij_im, delij=0.
      real, dimension (3) :: e1, e2, kvec, VCT_re, VCT_im
!     real, dimension(nk) :: specGWs   ,specGWh   ,specStr
!     real, dimension(nk) :: specGWshel,specGWhhel,specStrhel
      integer :: i,j,p,q,ik,ikx,iky,ikz,stat,ij,pq,ip,jq,jStress_ij
      logical :: lscale_tobox1=.true.
      real :: kscale_factor, fact, sign_switch
      real :: ksqr, one_over_k2, k1, k2, k3, k1sqr, k2sqr, k3sqr
      real :: hhTre, hhTim, hhXre, hhXim, coefAre, coefAim
      real :: ggTre, ggTim, ggXre, ggXim, coefBre, coefBim
      real :: cosot, sinot, om12, om1, om, omt1
      real :: SCLij_re, SCLij_im
      intent(inout) :: f
      character (len=2) :: label
!
!  Check that the relevant arrays are registered
!
      if (.not.lStress_as_aux.and.label=='St') call fatal_error('compute_gT_and_gX_from_gij','lStress_as_aux must be true')
!
!  For testing purposes, if lno_transverse_part=T, we would not need to
!  compute the Fourier transform, so we would skip the rest.
!
!  Allocate memory for arrays.
!
      allocate(S_T_re(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for S_T_re')
!
      allocate(S_T_im(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for S_T_im')
!
      allocate(S_X_re(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for S_X_re')
!
      allocate(S_X_im(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for S_X_im')
!
      allocate(Tpq_re(nx,ny,nz,6),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for Tpq_re')
!
      allocate(Tpq_im(nx,ny,nz,6),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for Tpq_im')
!
      allocate(kx(nxgrid),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij', &
          'Could not allocate memory for kx')
      allocate(ky(nygrid),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij', &
          'Could not allocate memory for ky')
      allocate(kz(nzgrid),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij', &
          'Could not allocate memory for kz')
!
!  calculate k^2
!
      kscale_factor=1
      if (lscale_tobox1) kscale_factor=2*pi/Lx
      kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2)*kscale_factor
!
      kscale_factor=1
      if (lscale_tobox1) kscale_factor=2*pi/Ly
      ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2)*kscale_factor
!
      kscale_factor=1
      if (lscale_tobox1) kscale_factor=2*pi/Lz
      kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2)*kscale_factor
!
!  set delta_ij
!
      delij(1)=1.
      delij(2)=1.
      delij(3)=1.
      delij(4)=0.
      delij(5)=0.
      delij(6)=0.
!
!  Set k^2 array. Note that in Fourier space, kz is the fastest index and has
!  the full nx extent (which, currently, must be equal to nxgrid).
!  But call it one_over_k2.
!
!  Assemble stress, Tpq
!
      Tpq_re=0.0
      Tpq_im=0.0
      do ij=1,6
        if (label=='St') then
          jStress_ij=iStress_ij-1+ij
          Tpq_re(:,:,:,ij)=f(l1:l2,m1:m2,n1:n2,jStress_ij)
        endif
      enddo
!
!  Fourier transform all 6 components
!
      do ij=1,6
        call fourier_transform(Tpq_re(:,:,:,ij),Tpq_im(:,:,:,ij))
!--     call fft_xyz_parallel(Tpq_re(:,:,:,ij),Tpq_im(:,:,:,ij))
      enddo
!
!  Set ST=SX=0 and reset all spectra.
!
      S_T_re=0. ; S_T_im=0.
      S_X_re=0. ; S_X_im=0.
      specGWs=0.; specGWshel=0.
      specGWh=0.; specGWhhel=0.
      specGWm=0.; specGWmhel=0.
      specStr=0.; specStrhel=0.
      specSCL=0.; specVCT=0.; specTpq=0.
!
!  P11, P22, P33, P12, P23, P31
!
      do iky=1,nz
        do ikx=1,ny
          do ikz=1,nx
!
!  compute e_T and e_X; determine first preferred direction,
!  which is a component with the smallest component by modulus.
!
            k1=kx(ikx+ipy*ny)
            k2=ky(iky+ipz*nz)
            k3=kz(ikz+ipx*nx)
            k1sqr=k1**2
            k2sqr=k2**2
            k3sqr=k3**2
            ksqr=k1sqr+k2sqr+k3sqr
!
!  find two vectors e1 and e2 to compute e_T and e_X
!
            if (lroot.and.ikx==1.and.iky==1.and.ikz==1) then
              e1=0.
              e2=0.
              Pij(1)=0.
              Pij(2)=0.
              Pij(3)=0.
              Pij(4)=0.
              Pij(5)=0.
              Pij(6)=0.
              one_over_k2=0.
            else
              one_over_k2=1./ksqr
              if(abs(k1)<abs(k2)) then
                if(abs(k1)<abs(k3)) then !(k1 is pref dir)
                  e1=(/0.,-k3,+k2/)
                  e2=(/k2sqr+k3sqr,-k2*k1,-k3*k1/)
                else !(k3 is pref dir)
                  e1=(/k2,-k1,0./)
                  e2=(/k1*k3,k2*k3,-(k1sqr+k2sqr)/)
                endif
              else !(k2 smaller than k1)
                if(abs(k2)<abs(k3)) then !(k2 is pref dir)
                  e1=(/-k3,0.,+k1/)
                  e2=(/+k1*k2,-(k1sqr+k3sqr),+k3*k2/)
                else !(k3 is pref dir)
                  e1=(/k2,-k1,0./)
                  e2=(/k1*k3,k2*k3,-(k1sqr+k2sqr)/)
                endif
              endif
              e1=e1/sqrt(e1(1)**2+e1(2)**2+e1(3)**2)
              e2=e2/sqrt(e2(1)**2+e2(2)**2+e2(3)**2)
              Pij(1)=1.-k1sqr*one_over_k2
              Pij(2)=1.-k2sqr*one_over_k2
              Pij(3)=1.-k3sqr*one_over_k2
              Pij(4)=-k1*k2*one_over_k2
              Pij(5)=-k2*k3*one_over_k2
              Pij(6)=-k3*k1*one_over_k2
            endif
!
!  set k vector
!
            kvec(1)=k1
            kvec(2)=k2
            kvec(3)=k3
!
!  compute e_T and e_X
!
            do j=1,3
            do i=1,3
              ij=ij_table(i,j)
              e_T(ij)=e1(i)*e1(j)-e2(i)*e2(j)
              e_X(ij)=e1(i)*e2(j)+e2(i)*e1(j)
            enddo
            enddo
!
!  possibility of swapping the sign
!
            sign_switch=1.
            if (lswitch_sign_e_X) then
              if (k3<0.) then
                sign_switch=-1.
                e_X=-e_X
              elseif (k3==0.) then
                if (k2<0.) then
                  sign_switch=-1.
                  e_X=-e_X
                elseif (k2==0.) then
                  if (k1<0.) then
                    sign_switch=-1.
                    e_X=-e_X
                  endif
                endif
              endif
            endif
!
!  Traceless-tansverse projection:
!  Sij = (Pip*Pjq-.5*Pij*Ppq)*Tpq
!
            Sij_re=0.
            Sij_im=0.
            do j=1,3
            do i=1,j
            do q=1,3
            do p=1,3
              ij=ij_table(i,j)
              pq=ij_table(p,q)
              ip=ij_table(i,p)
              jq=ij_table(j,q)
              Sij_re(ij)=Sij_re(ij)+(Pij(ip)*Pij(jq)-.5*Pij(ij)*Pij(pq))*Tpq_re(ikz,ikx,iky,pq)
              Sij_im(ij)=Sij_im(ij)+(Pij(ip)*Pij(jq)-.5*Pij(ij)*Pij(pq))*Tpq_im(ikz,ikx,iky,pq)
            enddo
            enddo
            enddo
            enddo
!
            SCLij_re=0.
            SCLij_im=0.
            do q=1,3
            do p=1,3
              pq=ij_table(p,q)
              SCLij_re=SCLij_re+(Pij(pq)-delij(pq))*Tpq_re(ikz,ikx,iky,pq)
              SCLij_im=SCLij_im+(Pij(pq)-delij(pq))*Tpq_im(ikz,ikx,iky,pq)
            enddo
            enddo
!
            do q=1,3
              VCT_re(q)=+twothird*kvec(q)*SCLij_im
              VCT_im(q)=-twothird*kvec(q)*SCLij_re
              do p=1,3
                pq=ij_table(p,q)
                VCT_re(q)=VCT_re(q)+2.*kvec(p)*one_over_k2*Tpq_im(ikz,ikx,iky,pq)
                VCT_im(q)=VCT_im(q)-2.*kvec(p)*one_over_k2*Tpq_re(ikz,ikx,iky,pq)
              enddo
            enddo
!
!  Compute S_T and S_X. Loop over all i and j for simplicity to avoid
!  special treatments of diagonal and off-diagonal terms,
!  AB: it looks like one could reuse part of Tpq_re and Tpq_im for S_T_re, ...,
!  AB: S_X_im to save memory
!
            do j=1,3
            do i=1,3
              ij=ij_table(i,j)
              S_T_re(ikz,ikx,iky)=S_T_re(ikz,ikx,iky)+.5*e_T(ij)*Sij_re(ij)
              S_T_im(ikz,ikx,iky)=S_T_im(ikz,ikx,iky)+.5*e_T(ij)*Sij_im(ij)
              S_X_re(ikz,ikx,iky)=S_X_re(ikz,ikx,iky)+.5*e_X(ij)*Sij_re(ij)
              S_X_im(ikz,ikx,iky)=S_X_im(ikz,ikx,iky)+.5*e_X(ij)*Sij_im(ij)
            enddo
            enddo
!
!  Compute exact solution for hT, hX, gT, and gX in Fourier space.
!
            hhTre=f(nghost+ikz,nghost+ikx,nghost+iky,ihhT  )
            hhXre=f(nghost+ikz,nghost+ikx,nghost+iky,ihhX  )
            hhTim=f(nghost+ikz,nghost+ikx,nghost+iky,ihhTim)
            hhXim=f(nghost+ikz,nghost+ikx,nghost+iky,ihhXim)
!
            ggTre=f(nghost+ikz,nghost+ikx,nghost+iky,iggT  )
            ggXre=f(nghost+ikz,nghost+ikx,nghost+iky,iggX  )
            ggTim=f(nghost+ikz,nghost+ikx,nghost+iky,iggTim)
            ggXim=f(nghost+ikz,nghost+ikx,nghost+iky,iggXim)
!
!  compute omega (but assume c=1), omega*t, etc.
!
            om12=one_over_k2
            om1=sqrt(om12)
            om=1./om1
            omt1=1./(om*t)
!
!  compute cos(om*dt) and sin(om*dt) to get from one timestep to the next.
!
            cosot=cos(om*dt)
            sinot=sin(om*dt)
!
!  Solve wave equation for hT and gT from one timestep to the next.
!
            coefAre=(hhTre-om12*S_T_re(ikz,ikx,iky))
            coefAim=(hhTim-om12*S_T_im(ikz,ikx,iky))
            coefBre=ggTre*om1 !+omt1*om12*S_T_re(ikz,ikx,iky)
            coefBim=ggTim*om1 !+omt1*om12*S_T_im(ikz,ikx,iky)
            f(nghost+ikz,nghost+ikx,nghost+iky,ihhT  )=coefAre*cosot+coefBre*sinot+om12*S_T_re(ikz,ikx,iky)
            f(nghost+ikz,nghost+ikx,nghost+iky,ihhTim)=coefAim*cosot+coefBim*sinot+om12*S_T_im(ikz,ikx,iky)
            f(nghost+ikz,nghost+ikx,nghost+iky,iggT  )=coefBre*cosot*om-coefAre*om*sinot !-omt1*om1*S_T_re(ikz,ikx,iky)
            f(nghost+ikz,nghost+ikx,nghost+iky,iggTim)=coefBim*cosot*om-coefAim*om*sinot !-omt1*om1*S_T_im(ikz,ikx,iky)
!
!  Debug output
!
            if (ldebug_print) then
              if (nint(k1)==2.and.nint(k2)==0.and.nint(k3)==0) then
                print*,'AXEL0: ',coefAre,coefBre,hhTre,ggTre
              endif
            endif
!
!  Solve wave equation for hX and gX from one timestep to the next.
!
            coefAre=(hhXre-om12*S_X_re(ikz,ikx,iky))
            coefAim=(hhXim-om12*S_X_im(ikz,ikx,iky))
            coefBre=ggXre*om1 !+omt1*om12*S_X_re(ikz,ikx,iky)
            coefBim=ggXim*om1 !+omt1*om12*S_X_im(ikz,ikx,iky)
            f(nghost+ikz,nghost+ikx,nghost+iky,ihhX  )=coefAre*cosot+coefBre*sinot+om12*S_X_re(ikz,ikx,iky)
            f(nghost+ikz,nghost+ikx,nghost+iky,ihhXim)=coefAim*cosot+coefBim*sinot+om12*S_X_im(ikz,ikx,iky)
            f(nghost+ikz,nghost+ikx,nghost+iky,iggX  )=coefBre*cosot*om-coefAre*om*sinot !-omt1*om1*S_X_re(ikz,ikx,iky)
            f(nghost+ikz,nghost+ikx,nghost+iky,iggXim)=coefBim*cosot*om-coefAim*om*sinot !-omt1*om1*S_X_im(ikz,ikx,iky)
!
!  Set origin to zero. It is given by (1,1,1) on root processor.
!
            if (lroot) f(nghost+1,nghost+1,nghost+1,ihhT  ) = 0.
            if (lroot) f(nghost+1,nghost+1,nghost+1,ihhTim) = 0.
            if (lroot) f(nghost+1,nghost+1,nghost+1,iggT  ) = 0.
            if (lroot) f(nghost+1,nghost+1,nghost+1,iggTim) = 0.
!
            if (lroot) f(nghost+1,nghost+1,nghost+1,ihhX  ) = 0.
            if (lroot) f(nghost+1,nghost+1,nghost+1,ihhXim) = 0.
            if (lroot) f(nghost+1,nghost+1,nghost+1,iggX  ) = 0.
            if (lroot) f(nghost+1,nghost+1,nghost+1,iggXim) = 0.
!
!  Set stress components in f-array.
!
            f(nghost+ikz,nghost+ikx,nghost+iky,iStressT  )=S_T_re(ikz,ikx,iky)
            f(nghost+ikz,nghost+ikx,nghost+iky,iStressTim)=S_T_im(ikz,ikx,iky)
            f(nghost+ikz,nghost+ikx,nghost+iky,iStressX  )=S_X_re(ikz,ikx,iky)
            f(nghost+ikz,nghost+ikx,nghost+iky,iStressXim)=S_X_im(ikz,ikx,iky)
!
!  Sum up energy and helicity spectra. Divide by kscale_factor to have integers
!  for the Fortran index ik. Note, however, that 
!
            if (lspec) then
              ik=1+nint(sqrt(ksqr)/kscale_factor)
!
!  Debug output
!
              if (ldebug_print) then
                if (ik <= 5) write(*,1000) iproc,ik,k1,k2,k3,f(nghost+ikz,nghost+ikx,nghost+iky,iggX  )
                1000 format(2i5,1p,4e11.2)
              endif
!
              if (ik <= nk) then
!
!  Gravitational wave energy spectrum computed from hdot (=g)
!
                if (GWs_spec) then
                  specGWs(ik)=specGWs(ik) &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,iggX  )**2 &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,iggXim)**2 &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,iggT  )**2 &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,iggTim)**2
                  specGWshel(ik)=specGWshel(ik)+2*sign_switch*( &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,iggXim) &
                     *f(nghost+ikz,nghost+ikx,nghost+iky,iggT  ) &
                     -f(nghost+ikz,nghost+ikx,nghost+iky,iggX  ) &
                     *f(nghost+ikz,nghost+ikx,nghost+iky,iggTim) )
                endif
!
!  Gravitational wave strain spectrum computed from h
!
                if (GWh_spec) then
                  specGWh(ik)=specGWh(ik) &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,ihhX  )**2 &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,ihhXim)**2 &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,ihhT  )**2 &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,ihhTim)**2
                  specGWhhel(ik)=specGWhhel(ik)+2*sign_switch*( &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,ihhXim) &
                     *f(nghost+ikz,nghost+ikx,nghost+iky,ihhT  ) &
                     -f(nghost+ikz,nghost+ikx,nghost+iky,ihhX  ) &
                     *f(nghost+ikz,nghost+ikx,nghost+iky,ihhTim) )
                endif
!
!  Gravitational wave mixed spectrum computed from h and g
!
                if (GWm_spec) then
                  specGWm(ik)=specGWm(ik) &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,ihhX)  *f(nghost+ikz,nghost+ikx,nghost+iky,iggX  ) &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,ihhXim)*f(nghost+ikz,nghost+ikx,nghost+iky,iggXim) &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,ihhT)  *f(nghost+ikz,nghost+ikx,nghost+iky,iggT  ) &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,ihhTim)*f(nghost+ikz,nghost+ikx,nghost+iky,iggTim)
                  specGWmhel(ik)=specGWmhel(ik)-sign_switch*( &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,ihhXim) &
                     *f(nghost+ikz,nghost+ikx,nghost+iky,iggT  ) &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,iggXim) &
                     *f(nghost+ikz,nghost+ikx,nghost+iky,ihhT  ) &
                     -f(nghost+ikz,nghost+ikx,nghost+iky,ihhX  ) &
                     *f(nghost+ikz,nghost+ikx,nghost+iky,iggTim) &
                     -f(nghost+ikz,nghost+ikx,nghost+iky,iggX  ) &
                     *f(nghost+ikz,nghost+ikx,nghost+iky,ihhTim) )
                endif
!
!  Stress spectrum computed from Str
!  ?not used currently
!
                if (Str_spec) then
                  specStr(ik)=specStr(ik) &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,iStressX  )**2 &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,iStressXim)**2 &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,iStressT  )**2 &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,iStressTim)**2
                  specStrhel(ik)=specStrhel(ik)+2*sign_switch*( &
                     +f(nghost+ikz,nghost+ikx,nghost+iky,iStressXim) &
                     *f(nghost+ikz,nghost+ikx,nghost+iky,iStressT  ) &
                     -f(nghost+ikz,nghost+ikx,nghost+iky,iStressX  ) &
                     *f(nghost+ikz,nghost+ikx,nghost+iky,iStressTim) )
                endif
!
!  Stress spectrum computed from the scalar mode, SCL.
!
                if (SCL_spec) then
                  specSCL(ik)=specSCL(ik)+.5*(SCLij_re**2+SCLij_im**2)
                endif
!
!  Spectrum computed from the vector modes...
!
                if (VCT_spec) then
                  do q=1,3
                    specVCT(ik)=specVCT(ik)+.5*(VCT_re(q)**2+VCT_im(q)**2)
                  enddo
                endif
!
!  Spectrum computed from the total unprojected stress
!
                if (Tpq_spec) then
                  do q=1,3
                  do p=1,3
                    pq=ij_table(p,q)
                    specTpq(ik)=specTpq(ik)+.5*(Tpq_re(ikz,ikx,iky,pq)**2+Tpq_im(ikz,ikx,iky,pq)**2)
                  enddo
                  enddo
                endif
!
              endif
!
!  end of lspec
!
            endif
!
!  end of ikx, iky, and ikz loops
!
          enddo
        enddo
      enddo
!
!  Communicate as shared variables
!
!     if (lspec) then
!if (ldebug_print) print*,'AXEL5: iproc,specGWs=',iproc,specGWs
!       if (GWs_spec) then
!         call put_shared_variable('specGWs   ',specGWs   )
!         call put_shared_variable('specGWshel',specGWshel)
!       endif
!       if (GWh_spec) then
!         call put_shared_variable('specGWh   ',specGWh   )
!         call put_shared_variable('specGWhhel',specGWhhel)
!       endif
!       if (Str_spec) then
!         call put_shared_variable('specStr   ',specStr   )
!         call put_shared_variable('specStrhel',specStrhel)
!       endif
!     endif
!
!  back to real space
!
  !   call fourier_transform(S_T_re,S_T_im,linv=.true.)
  !   call fourier_transform(S_X_re,S_X_im,linv=.true.)
!--   call fft_xyz_parallel(S_T_re,S_T_im,linv=.true.)
!--   call fft_xyz_parallel(S_X_re,S_X_im,linv=.true.)
!
!  Deallocate arrays.
!
      if (allocated(S_T_re)) deallocate(S_T_re)
      if (allocated(S_X_im)) deallocate(S_X_im)
      if (allocated(Tpq_re)) deallocate(Tpq_re)
      if (allocated(Tpq_im)) deallocate(Tpq_im)
      if (allocated(kx)) deallocate(kx)
      if (allocated(ky)) deallocate(ky)
      if (allocated(kz)) deallocate(kz)
!
    endsubroutine compute_gT_and_gX_from_gij
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
      use Diagnostics
!!      use FArrayManager, only: farray_index_append
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!!!
!!!  reset everything in case of reset
!!!  (this needs to be consistent with what is defined above!)
!!!
      if (lreset) then
        idiag_g11pt=0; idiag_g22pt=0; idiag_g33pt=0
        idiag_g12pt=0; idiag_g23pt=0; idiag_g31pt=0
        idiag_hhTpt=0; idiag_hhXpt=0; idiag_ggTpt=0; idiag_ggXpt=0
        idiag_hhTp2=0; idiag_hhXp2=0; idiag_ggTp2=0; idiag_ggXp2=0
        idiag_hhT2m=0; idiag_hhX2m=0; idiag_hhTXm=0; idiag_hrms=0
        idiag_ggT2m=0; idiag_ggX2m=0; idiag_ggTXm=0; idiag_gg2m=0
        idiag_EEGW=0
        cformv=''
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'g11pt',idiag_g11pt)
        call parse_name(iname,cname(iname),cform(iname),'g22pt',idiag_g22pt)
        call parse_name(iname,cname(iname),cform(iname),'g33pt',idiag_g33pt)
        call parse_name(iname,cname(iname),cform(iname),'g12pt',idiag_g12pt)
        call parse_name(iname,cname(iname),cform(iname),'g23pt',idiag_g23pt)
        call parse_name(iname,cname(iname),cform(iname),'g31pt',idiag_g31pt)
        if (lhhTX_as_aux) then
          call parse_name(iname,cname(iname),cform(iname),'hhTpt',idiag_hhTpt)
          call parse_name(iname,cname(iname),cform(iname),'hhXpt',idiag_hhXpt)
          call parse_name(iname,cname(iname),cform(iname),'hhTp2',idiag_hhTp2)
          call parse_name(iname,cname(iname),cform(iname),'hhXp2',idiag_hhXp2)
          call parse_name(iname,cname(iname),cform(iname),'hrms',idiag_hrms)
          call parse_name(iname,cname(iname),cform(iname),'hhT2m',idiag_hhT2m)
          call parse_name(iname,cname(iname),cform(iname),'hhX2m',idiag_hhX2m)
          call parse_name(iname,cname(iname),cform(iname),'hhTXm',idiag_hhTXm)
        endif
        if (lggTX_as_aux) then
          call parse_name(iname,cname(iname),cform(iname),'ggTpt',idiag_ggTpt)
          call parse_name(iname,cname(iname),cform(iname),'ggXpt',idiag_ggXpt)
          call parse_name(iname,cname(iname),cform(iname),'ggTp2',idiag_ggTp2)
          call parse_name(iname,cname(iname),cform(iname),'ggXp2',idiag_ggXp2)
          call parse_name(iname,cname(iname),cform(iname),'EEGW',idiag_EEGW)
          call parse_name(iname,cname(iname),cform(iname),'gg2m',idiag_gg2m)
          call parse_name(iname,cname(iname),cform(iname),'ggT2m',idiag_ggT2m)
          call parse_name(iname,cname(iname),cform(iname),'ggX2m',idiag_ggX2m)
          call parse_name(iname,cname(iname),cform(iname),'ggTXm',idiag_ggTXm)
        endif
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='hhT'.or.cnamev=='hhX'.or.cnamev=='ggT'.or.cnamev=='ggX'.or. &
              cnamev=='h22'.or.cnamev=='h33'.or.cnamev=='h23') cformv='DEFINED'
      endif
!
!!!  write column where which magnetic variable is stored
!!      if (lwr) then
!!        call farray_index_append('i_SPECIAL_DIAGNOSTIC',i_SPECIAL_DIAGNOSTIC)
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
      use Slices_methods, only: assign_slices_scal
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  hhT
!
        case ('hhT'); call assign_slices_scal(slices,f,ihhT)
!
!  hhX
!
        case ('hhX'); call assign_slices_scal(slices,f,ihhX)
!
!  ggT
!
        case ('ggT'); call assign_slices_scal(slices,f,iggT)
!
!  ggX
!
        case ('ggX'); call assign_slices_scal(slices,f,iggX)
!
!  hh22
!
        case ('h22'); call assign_slices_scal(slices,f,ihij-1+2)
!
!  hh33
!
        case ('h33'); call assign_slices_scal(slices,f,ihij-1+3)
!
!  hh23
!
        case ('h23'); call assign_slices_scal(slices,f,ihij-1+5)
!
      endselect
!
!  The following is just a comment to remind ourselves how
!  the remaining 3 offdiagonal terms are being accessed.
!
      !ij_table(1,2)=4
      !ij_table(2,3)=5
      !ij_table(3,1)=6
      !ij_table(2,1)=4
      !ij_table(3,2)=5
      !ij_table(1,3)=6
!
    endsubroutine get_slices_special
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
