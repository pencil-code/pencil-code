! $Id$
!
!  This module provide a way for users to specify custom
!  (i.e. not in the standard Pencil Code) physics, diagnostics etc.
!
!  The module provides a set of standard hooks into the Pencil-Code and
!  currently allows the following customizations:
!
!  Description                                     | Relevant function
!  call
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
!  Special term in the mass (density) equation     |
!  special_calc_density
!  Special term in the momentum (hydro) equation   | special_calc_hydro
!  Special term in the energy equation             | special_calc_energy
!  Special term in the induction (magnetic)        |
!  special_calc_magnetic
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
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED mu5; gmu5(3); del2mu5
! PENCILS PROVIDED ugmu5
!! PENCILS PROVIDED ucurlb
!! PENCILS PROVIDED curlj(3); jjij(3,3)
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
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
! Declare index of new variables in f array (if any).
!
   real :: diffmu5, lambda5, mu5_const=0.
   real :: meanmu5
   real, pointer :: eta
   real :: cdtchiral=1.
   real, dimension (nx) :: diffus_mu5_1, diffus_mu5_2, diffus_mu5_3
   real, dimension (nx) :: diffus_bb_1, diffus_bb_2, diffus_uu_1,diffus_special
   real, dimension (nx) :: uxbjrms, mu5bjrms
   real, dimension (nx) :: uxbj, mu5bj
   integer :: imu5
!
  character (len=labellen) :: initspecial='nothing'
!
  namelist /special_init_pars/ &
      initspecial, mu5_const
!
  namelist /special_run_pars/ &
      diffmu5, lambda5, cdtchiral
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_mu5m=0      ! DIAG_DOC: $\left<\mu_5\right>$
  integer :: idiag_mu5rms=0    ! DIAG_DOC: $\left<\mu_5^2\right>^{1/2}$
  integer :: idiag_bgmu5rms=0  ! DIAG_DOC:$\left<(\Bv\cdot\nabla\mu_5)^2\right>^{1/2}$
  integer :: idiag_diffus_mu5_1=0
  integer :: idiag_diffus_mu5_2=0
  integer :: idiag_diffus_mu5_3=0
  integer :: idiag_diffus_bb_1=0
  integer :: idiag_diffus_bb_2=0
  integer :: idiag_diffus_uu_1=0
  integer :: idiag_uxbjm=0
  integer :: idiag_mu5bjm=0
!  integer :: idiag_bcurljm=0
  integer :: idiag_uxbjrms=0
  integer :: idiag_mu5bjrms=0
!  integer :: idiag_bcurljrms=0
!
  contains
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
      use FArrayManager, only: farray_register_pde
!
      if (lroot) call svn_id( &
           "$Id$")
!
      call farray_register_pde('mu5',imu5)
!!      call farray_register_auxiliary('specaux',ispecaux)
!!      call
!farray_register_auxiliary('specaux',ispecaux,communicated=.true.)
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
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      use SharedVariables, only : get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ierr
!
      call keep_compiler_quiet(f)
!
      if (lmagnetic.and.lrun) then
          call get_shared_variable('eta',eta,ierr)
          if (ierr/=0) call fatal_error("initialize_special: ", &
              "cannot get shared var eta")
      endif
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
      use Initcond
!
      real, dimension (mx,my,mz,mfarray) :: f,df
!
      intent(inout) :: f
!
!  initial conditions
!
      select case (initspecial)
!
        case ('nothing'); if (lroot) print*,'init_special: nothing'
!
        case ('zero')
          f(:,:,:,imu5) = 0.
!
        case ('const')
          f(:,:,:,imu5) = mu5_const
!
        case default
          call fatal_error("init_special: No such value for initspecial:" &
              ,trim(initspecial))
      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
      lpenc_requested(i_b2)=.true.
      lpenc_requested(i_mu5)=.true.
!      lpenc_requested(i_curlj)=.true.
!      lpenc_requested(i_ucurlb)=.true.
      lpenc_requested(i_ugmu5)=.true.
      if (ldt) lpenc_requested(i_rho1)=.true.
!      lpenc_requested(i_jjij)=.true.
      if (diffmu5/=0.) lpenc_requested(i_del2mu5)=.true.
      if (lhydro.or.lhydro_kinematic) lpenc_requested(i_uu)=.true.
      if (lmagnetic) lpenc_requested(i_bb)=.true.
!      if (lmagnetic) lpenc_requested(i_jij)=.true.
!      if (lmagnetic.and.lhydro) lpenc_requested(i_ub)=.true.
      if (lmagnetic.and.lhydro) lpenc_requested(i_jb)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified
!  here.
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
      use Sub, only: del2, dot2_mn, del2v_etc, grad, dot, u_dot_grad, gij
      use Sub, only: multsv, curl, curl_mn
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      if (lpencil(i_mu5)) p%mu5=f(l1:l2,m,n,imu5)
!      if (lpencil(i_ucurlb)) call dot(p%uu,p%jj,p%ucurlb)
      if (lpencil(i_ugmu5)) call dot(p%uu,p%gmu5,p%ugmu5)
      if (lpencil(i_del2mu5)) call del2(f,imu5,p%del2mu5)
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
!  06-oct-03/tony: coded
!
      use Sub, only: multsv, dot_mn, dot2_mn, dot_mn_vm_trans, dot
      use Diagnostics, only: sum_mn_name
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: chiraldiffusion
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
      real, dimension (nx) :: bgmu5, EB, uujj, bbjj!, bcurlj
      real, dimension (nx,3) :: mu5bb
      real, parameter :: alpha_fine_structure=1./137.
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!!    if (headtt) call identify_bcs('mu5',imu5)
!
!  Compute E.B
!
      EB=eta*(p%jb-p%mu5*p%b2)
!
!  Evolution of mu5
!
      df(l1:l2,m,n,imu5) = df(l1:l2,m,n,imu5)-p%ugmu5 &
      +diffmu5*p%del2mu5+lambda5*EB
      diffus_mu5_1 = lambda5*eta*p%b2/(p%mu5)*sqrt(dxyz_2)
      diffus_mu5_2 = lambda5*eta*p%b2
      diffus_mu5_3 = diffmu5*dxyz_2
!                          
!  Additions to evolution of bb
!
      if (lmagnetic) then
        call multsv(p%mu5,p%bb,mu5bb)
         df(l1:l2,m,n,iax:iaz) = df(l1:l2,m,n,iax:iaz)+eta*(mu5bb)
      endif
      diffus_bb_1 = eta*p%mu5*sqrt(dxyz_2)
!
!  Additions to evolution of uu
!
      if (lhydro) then
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz)
      endif  
!
!  Contribution to the time-step
!
      if (lfirst.and.ldt) then
        diffus_special = cdtchiral*max(diffus_mu5_1, diffus_mu5_2, diffus_mu5_3, &
        diffus_bb_1, diffus_bb_2, diffus_uu_1)
      endif
!
!  diagnostics
!
      if (ldiagnos) then
        if (idiag_mu5m/=0) call sum_mn_name(p%mu5,idiag_mu5m)
        if (idiag_mu5rms/=0) call sum_mn_name(p%mu5**2,idiag_mu5rms,lsqrt=.true.)
        endif
        if (idiag_bgmu5rms/=0) then
          call dot_mn(p%bb,p%gmu5,bgmu5)
          call sum_mn_name(bgmu5**2,idiag_bgmu5rms,lsqrt=.true.)
        endif
!        if (idiag_diffus_mu5_1/=0)  call
!        sum_mn_name(diffus_mu5_1,idiag_diffus_mu5_1)
!        if (idiag_diffus_mu5_2/=0)  call
!        sum_mn_name(diffus_mu5_2,idiag_diffus_mu5_2)
!        if (idiag_diffus_mu5_3/=0)  call
!        sum_mn_name(diffus_mu5_3,idiag_diffus_mu5_3)
!        if (idiag_diffus_gtheta5_1/=0)  call
!        sum_mn_name(diffus_gtheta5_1,idiag_diffus_gtheta5_1)
!        if (idiag_diffus_gtheta5_2/=0)  call
!        sum_mn_name(diffus_gtheta5_2,idiag_diffus_gtheta5_2)
!        if (idiag_diffus_bb_1/=0)  call
!        sum_mn_name(diffus_bb_1,idiag_diffus_bb_1)
!        if (idiag_diffus_bb_2/=0)  call
!        sum_mn_name(diffus_bb_2,idiag_diffus_bb_2)
!        if (idiag_diffus_uu_1/=0)  call
!        sum_mn_name(diffus_uu_1,idiag_diffus_uu_1)
        if (idiag_uxbjm/=0)  call sum_mn_name(uxbj,idiag_uxbjm)
        if (idiag_mu5bjm/=0)  call sum_mn_name(mu5bj,idiag_mu5bjm)
!        if (idiag_bcurljm/=0)  call sum_mn_name(bcurlj,idiag_bcurljm)
!        if (idiag_test/=0)  call
!        sum_mn_name(test2,idiag_test,lsqrt=.true.)
        if (idiag_uxbjrms/=0)  call sum_mn_name(uxbj**2,idiag_uxbjrms,lsqrt=.true.)
        if (idiag_mu5bjrms/=0)  call sum_mn_name(mu5bj**2,idiag_mu5bjrms,lsqrt=.true.)
!        if (idiag_bcurljrms/=0)  call
!        sum_mn_name(bcurlj**2,idiag_bcurljrms,lsqrt=.true.)
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
      use Diagnostics, only: parse_name
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
        idiag_mu5m=0; idiag_mu5rms=0; idiag_bgmu5rms=0;
        idiag_diffus_mu5_1=0; idiag_diffus_mu5_2=0; idiag_diffus_mu5_3=0; 
        idiag_diffus_bb_1=0; idiag_diffus_bb_2=0; idiag_diffus_uu_1=0; 
        idiag_uxbjm=0;
        idiag_mu5bjm=0; idiag_uxbjrms=0;
        idiag_mu5bjrms=0!; idiag_bcurljm=0; idiag_bcurljrms=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'mu5m',idiag_mu5m)
        call parse_name(iname,cname(iname),cform(iname),'mu5rms',idiag_mu5rms)
        call parse_name(iname,cname(iname),cform(iname),'bgmu5rms',idiag_bgmu5rms)
        call parse_name(iname,cname(iname),cform(iname),'diffus_mu5_1',idiag_diffus_mu5_1)
        call parse_name(iname,cname(iname),cform(iname),'diffus_mu5_2',idiag_diffus_mu5_2)
        call parse_name(iname,cname(iname),cform(iname),'diffus_mu5_3',idiag_diffus_mu5_3)
        call parse_name(iname,cname(iname),cform(iname),'diffus_bb_1',idiag_diffus_bb_1)
        call parse_name(iname,cname(iname),cform(iname),'diffus_bb_2',idiag_diffus_bb_2)
        call parse_name(iname,cname(iname),cform(iname),'diffus_uu_1',idiag_diffus_uu_1)
        call parse_name(iname,cname(iname),cform(iname),'uxbjm',idiag_uxbjm)
        call parse_name(iname,cname(iname),cform(iname),'mu5bjm',idiag_mu5bjm)
!        call
!        parse_name(iname,cname(iname),cform(iname),'bcurljm',idiag_bcurljm)
        call parse_name(iname,cname(iname),cform(iname),'uxbjrms',idiag_uxbjrms)
        call parse_name(iname,cname(iname),cform(iname),'mu5bjrms',idiag_mu5bjrms)
!        call
!        parse_name(iname,cname(iname),cform(iname),'bcurljrms',idiag_bcurljrms)
      enddo
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
    subroutine calc_lspecial_pars(f)
!
!  Calculate meanmu5, which should be subtracted from mu5 in eqn for
!  gtheta5
!
!  11-oct-15/jenny: coded
!
      use Mpicomm, only: mpiallreduce_sum
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: fact, meanmu5_tmp
      intent(inout) :: f
!
!  compute meanmu5
!
      meanmu5=0.
      do n=n1,n2
      do m=m1,m2
        meanmu5=meanmu5+sum(f(l1:l2,m,n,imu5))
      enddo
      enddo
!
!  communicate and divide by all mesh meshpoints
!
      call mpiallreduce_sum(meanmu5,meanmu5_tmp,1)
      fact=1./(nw*ncpus)
      meanmu5=fact*meanmu5_tmp
!
    endsubroutine calc_lspecial_pars
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
    subroutine special_calc_particles(f,fp,ineargrid)
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
    subroutine special_after_timestep(f,df,dt_)
!
!  Possibility to modify the f and df after df is updated.
!  Used for the Fargo shift, for instance.
!
!  27-nov-08/wlad: coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, intent(in) :: dt_
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dt_)
!
    endsubroutine  special_after_timestep
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
  subroutine initialize_mult_special
          
  endsubroutine initialize_mult_special

!***********************************************************************
  subroutine finalize_mult_special
        

  endsubroutine finalize_mult_special
!*********************************************************************** 
endmodule Special
