! $Id$
!
!  Advecto-resistive gauge, dLamRA/dt = -u.gradLamRA - U.A_resistive + eta*del2A_resistive
!  A_advecto-resistive=A_resistive+grad(LamRA)
!
!  25-feb-07/axel: adapted from nospecial.f90
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
! PENCILS PROVIDED gLamRA(3)
!***************************************************************
!
module Special
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../special.h'
!
!  square of wave speed for gauge field
!
  real, pointer :: eta
  logical, pointer :: lweyl_gauge
  integer :: iaadv=0, iaadvx=0, iaadvy=0, iaadvz=0
!
  ! input parameters
  real :: ampl=1e-3,kx=1.,ky=0.,kz=0.
  logical :: ladvecto_resistive=.true.
  logical :: laa_adv_as_aux=.false.
  character(len=50) :: init='zero'
  namelist /special_init_pars/ &
    ladvecto_resistive,init,ampl,kx,ky,kz, &
    laa_adv_as_aux
!
  ! run parameters
  namelist /special_run_pars/ &
    ladvecto_resistive, &
    laa_adv_as_aux
!
! Declare any index variables necessary for main or
!
  integer :: iLamRA=0
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_LamRAm=0     ! DIAG_DOC: $\left<\Lambda_{r\to a}\right>$
  integer :: idiag_LamRApt=0    ! DIAG_DOC: $\Lambda_{r\to a}(x1,y1,z1)$
  integer :: idiag_LamRAp2=0    ! DIAG_DOC: $\Lambda_{r\to a}(x2,y2,z2)$
  integer :: idiag_LamRArms=0   ! DIAG_DOC: $\left<\Lambda_{r\to a}^2\right>^{1/2}$
  integer :: idiag_LamRAbzm=0   ! DIAG_DOC: $\left<\Lambda_{r\to a} B_z\right>$
  integer :: idiag_LamRAbzmz=0  ! DIAG_DOC: $\left<\Lambda_{r\to a} B_z\right>_{xy}$
  integer :: idiag_gLamRAbm=0     ! DIAG_DOC: $\left<\Lambda_{r\to a}\Bv\right>$
  integer :: idiag_apbrms=0     ! DIAG_DOC: $\left<(\Av'\Bv)^2\right>^{1/2}$
  integer :: idiag_jxarms=0     ! DIAG_DOC: $\left<(\Jv\times\Av)^2\right>^{1/2}$
  integer :: idiag_jxaprms=0    ! DIAG_DOC: $\left<(\Jv\times\Av')^2\right>^{1/2}$
  integer :: idiag_jxgLamRArms=0  ! DIAG_DOC: $\left<(\Jv\times\nabla\Lambda_{r\to a})^2\right>^{1/2}$
  integer :: idiag_gLamRArms=0    ! DIAG_DOC: $\left<(\nabla\Lambda_{r\to a})^2\right>^{1/2}$
  integer :: idiag_divabrms=0   ! DIAG_DOC: $\left<[(\nabla\cdot\Av)\Bv]^2\right>^{1/2}$
  integer :: idiag_divapbrms=0  ! DIAG_DOC: $\left<[(\nabla\cdot\Av')\Bv]^2\right>^{1/2}$
  integer :: idiag_d2LamRAbrms=0  ! DIAG_DOC: $\left<[(\nabla^2\Lambda_{r\to a})\Bv]^2\right>^{1/2}$
  integer :: idiag_d2LamRArms=0   ! DIAG_DOC: $\left<[\nabla^2\Lambda_{r\to a}]^2\right>^{1/2}$
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/tony: coded
!
      use FArrayManager
      use Sub, only: register_report_aux
!
      call farray_register_pde('LamRA',iLamRA)
      ispecialvar=iLamRA
!
!  Compute vector potential in advective gauge as auxiliary array.
!
      if (laa_adv_as_aux) call register_report_aux('aadv',iaadv,iaadvx,iaadvy,iaadvz)
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      use SharedVariables, only : get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ierr
!
!  Initialize module variables which are parameter dependent
!  wave speed of gauge potential
!
      if (lrun) then
        call get_shared_variable('lweyl_gauge',lweyl_gauge,ierr)
        if (ierr/=0) &
            call fatal_error("initialize_special: ", "cannot get lweyl_gauge")
        if (.not.lweyl_gauge) then
          call get_shared_variable('eta',eta,ierr)
          if (ierr/=0) &
              call fatal_error("initialize_special: ", "cannot get shared eta")
        endif
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Initcond
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      select case (init)
        case ('nothing'); if (lroot) print*,'init_special: nothing'
        case ('zero'); f(:,:,:,iLamRA)=0.
        case ('sinwave-x'); call sinwave(ampl,f,iLamRA,kx=kx)
        case ('sinwave-y'); call sinwave(ampl,f,iLamRA,ky=ky)
        case ('sinwave-z'); call sinwave(ampl,f,iLamRA,kz=kz)
!
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_special: No such value for init: ', trim(init)
          call stop_it("")
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
!  25-feb-07/axel: adapted
!
      if (.not.lweyl_gauge) then
        lpenc_requested(i_diva)=.true.
      endif
      lpenc_requested(i_aa)=.true.
      lpenc_requested(i_uu)=.true.
      lpenc_requested(i_gLamRA)=.true.
!
!  Diagnostic pencils
!
      if (idiag_apbrms/=0) lpenc_diagnos(i_ab)=.true.
      if (idiag_gLamRAbm/=0) lpenc_diagnos(i_bb)=.true.
      if (idiag_jxarms/=0 .or. idiag_jxaprms/=0 .or. idiag_jxgLamRArms/=0) lpenc_diagnos(i_jj)=.true.
      if (idiag_divabrms/=0 .or. idiag_divapbrms/=0 .or. idiag_d2LamRAbrms/=0) then
        lpenc_diagnos(i_bb)=.true.
        lpenc_diagnos(i_diva)=.true.
      endif
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate pencils for advective gauge.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-nov-04/tony: coded
!   7-mar-26/axel: made p%gLamRA a pencil
!
      use Sub, only: grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      call grad(f,iLamRA,p%gLamRA)
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
!  several precalculated Pencils of information are passed if for
!  efficiency.
!
!   06-oct-03/tony: coded
!
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: jxa, divab
      real, dimension (nx) :: LamRA, del2LamRA, ua, ugLamRA, gLamRAb, jxa2, divab2
!
      intent(in) :: p
      intent(inout) :: f,df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
!
      if (headtt) call identify_bcs('LamRA',iLamRA)
!
!  solve gauge condition
!
      if (lhydro.or.lhydro_kinematic) then
        call dot(p%uu,p%gLamRA,ugLamRA)
!
!  U.A term
!
        if (lmagnetic) then
          call dot(p%uu,p%aa,ua)
        else
          ua=0.
        endif
!
! based on resistive gauge
!
        if (ladvecto_resistive) then
          call del2(f,iLamRA,del2LamRA)
          df(l1:l2,m,n,iLamRA)=df(l1:l2,m,n,iLamRA)-ugLamRA-ua+eta*del2LamRA
        else
          df(l1:l2,m,n,iLamRA)=df(l1:l2,m,n,iLamRA)-ugLamRA-ua-eta*p%diva
        endif
      else
        call fatal_error('dspecial_dt','need hydro or hydro_kinematic for advective gauge')
      endif
!
!  Possibility of vector potential in advective gauge as auxiliary output.
!  Compute Aadv=AWeyl+gLamRA, where gLam is evolved in time.
!  By contrast to the Coulomb gauge, there is here a plus sign.
!
      if (laa_adv_as_aux) f(l1:l2,m,n,iaadvx:iaadvz)=p%aa+p%gLamRA
!
!  diagnostics
!
      if (ldiagnos) then
        LamRA=f(l1:l2,m,n,iLamRA)
        if (idiag_LamRAm/=0) call sum_mn_name(LamRA,idiag_LamRAm)
        if (idiag_LamRArms/=0) call sum_mn_name(LamRA**2,idiag_LamRArms,lsqrt=.true.)
        if (idiag_LamRAbzm/=0) call sum_mn_name(LamRA*p%bb(:,3),idiag_LamRAbzm)
        if (idiag_gLamRAbm/=0.or.idiag_apbrms/=0) then
          call dot(p%gLamRA,p%bb,gLamRAb)
          if (idiag_gLamRAbm/=0) call sum_mn_name(gLamRAb,idiag_gLamRAbm)
          if (idiag_apbrms/=0) call sum_mn_name((p%ab+gLamRAb)**2,idiag_apbrms,lsqrt=.true.)
        endif
!
!  (JxA^r)_rms
!
        if (idiag_jxarms/=0) then
          call cross(p%jj,p%aa,jxa)
          call dot2(jxa,jxa2)
          if (idiag_jxarms/=0) call sum_mn_name(jxa2,idiag_jxarms,lsqrt=.true.)
        endif
!
!  (JxA^ar)_rms
!
        if (idiag_jxarms/=0) then
          call cross(p%jj,p%aa+p%gLamRA,jxa)
          call dot2(jxa,jxa2)
          if (idiag_jxaprms/=0) call sum_mn_name(jxa2,idiag_jxaprms,lsqrt=.true.)
        endif
!
!  (JxgLamRA)_rms
!
        if (idiag_jxgLamRArms/=0) then
          call cross(p%jj,p%gLamRA,jxa)
          call dot2(jxa,jxa2)
          if (idiag_jxgLamRArms/=0) call sum_mn_name(jxa2,idiag_jxgLamRArms,lsqrt=.true.)
        endif
!
!  (gLamRA)_rms
!
        if (idiag_jxgLamRArms/=0) then
          call dot2(p%gLamRA,jxa2)
          if (idiag_gLamRArms/=0) call sum_mn_name(jxa2,idiag_gLamRArms,lsqrt=.true.)
        endif
!
!  [(divA^r)*B]_rms
!
        if (idiag_divabrms/=0) then
          call multsv(p%diva,p%bb,divab)
          call dot2(divab,divab2)
          if (idiag_divabrms/=0) call sum_mn_name(divab2,idiag_divabrms,lsqrt=.true.)
        endif
!
!  [(divA^ar)*B]_rms
!
        if (idiag_divapbrms/=0) then
          call multsv(p%diva+del2LamRA,p%bb,divab)
          call dot2(divab,divab2)
          if (idiag_divapbrms/=0) call sum_mn_name(divab2,idiag_divapbrms,lsqrt=.true.)
        endif
!
!  (del2LamRA*B)_rms
!
        if (idiag_d2LamRAbrms/=0) then
          call multsv(del2LamRA,p%bb,divab)
          call dot2(divab,divab2)
          if (idiag_d2LamRAbrms/=0) call sum_mn_name(divab2,idiag_d2LamRAbrms,lsqrt=.true.)
        endif
!
!  (del2LamRA)_rms
!
        if (idiag_d2LamRArms/=0) then
          if (idiag_d2LamRArms/=0) call sum_mn_name(del2LamRA**2,idiag_d2LamRArms,lsqrt=.true.)
        endif
!
!  check for point 1
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_LamRApt/=0) call save_name(LamRA(lpoint-nghost),idiag_LamRApt)
        endif
!
!  check for point 2
!
        if (lroot.and.m==mpoint2.and.n==npoint2) then
          if (idiag_LamRAp2/=0) call save_name(LamRA(lpoint2-nghost),idiag_LamRAp2)
        endif
!
      endif
!
      if (l1davgfirst .or. (ldiagnos .and. ldiagnos_need_zaverages)) then
        if (idiag_LamRAbzmz/=0)   call xysum_mn_name_z(p%bb(:,3),idiag_LamRAbzmz)
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
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!   06-oct-03/tony: coded
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
      use Sub
!
!  define counters
!
      integer :: iname,inamez
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
        idiag_LamRAm=0; idiag_LamRApt=0; idiag_LamRAp2=0; idiag_LamRArms=0
        idiag_gLamRAbm=0; idiag_apbrms=0
        idiag_LamRAbzm=0; idiag_LamRAbzmz=0
        idiag_jxarms=0; idiag_divabrms=0
        idiag_jxaprms=0; idiag_divapbrms=0
        idiag_jxgLamRArms=0; idiag_d2LamRAbrms=0
        idiag_gLamRArms=0; idiag_d2LamRArms=0
        cformv=''
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'LamRAm',idiag_LamRAm)
        call parse_name(iname,cname(iname),cform(iname),'LamRArms',idiag_LamRArms)
        call parse_name(iname,cname(iname),cform(iname),'LamRAbzm',idiag_LamRAbzm)
        call parse_name(iname,cname(iname),cform(iname),'LamRApt',idiag_LamRApt)
        call parse_name(iname,cname(iname),cform(iname),'LamRAp2',idiag_LamRAp2)
        call parse_name(iname,cname(iname),cform(iname),'gLamRAbm',idiag_gLamRAbm)
        call parse_name(iname,cname(iname),cform(iname),'apbrms',idiag_apbrms)
        call parse_name(iname,cname(iname),cform(iname),'jxarms',idiag_jxarms)
        call parse_name(iname,cname(iname),cform(iname),'jxaprms',idiag_jxaprms)
        call parse_name(iname,cname(iname),cform(iname),'jxgLamRArms',idiag_jxgLamRArms)
        call parse_name(iname,cname(iname),cform(iname),'gLamRArms',idiag_gLamRArms)
        call parse_name(iname,cname(iname),cform(iname),'divabrms',idiag_divabrms)
        call parse_name(iname,cname(iname),cform(iname),'divapbrms',idiag_divapbrms)
        call parse_name(iname,cname(iname),cform(iname),'d2LamRAbrms',idiag_d2LamRAbrms)
        call parse_name(iname,cname(iname),cform(iname),'d2LamRArms',idiag_d2LamRArms)
      enddo
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'LamRAbzmz',idiag_LamRAbzmz)
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='LamRA') cformv='DEFINED'
      endif
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        call farray_index_append('i_LamRAm',idiag_LamRAm)
        call farray_index_append('i_LamRArms',idiag_LamRArms)
        call farray_index_append('i_LamRAbzm',idiag_LamRAbzm)
        call farray_index_append('i_LamRAbzmz',idiag_LamRAbzmz)
        call farray_index_append('i_LamRApt',idiag_LamRApt)
        call farray_index_append('i_LamRAp2',idiag_LamRAp2)
        call farray_index_append('i_gLamRAbm',idiag_gLamRAbm)
        call farray_index_append('i_apbrms',idiag_apbrms)
        call farray_index_append('iLamRA',iLamRA)
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of electric potential
!
!  26-feb-07/axel: adapted from gross_pitaevskii
!
      use Slices_methods, only: assign_slices_scal

      real, dimension (mx,my,mz,mvar+maux) :: f
      type (slice_data) :: slices
!
      integer :: inamev
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  LamRA
!
        case ('LamRA'); call assign_slices_scal(slices,f,iLamRA)
!
      endselect
!
    endsubroutine get_slices_special
!***********************************************************************
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
!
endmodule Special
