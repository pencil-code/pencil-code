! $Id$
!
!  Advecto-resistive gauge, dLam/dt = -u.gradLambda - U.A_resistive + eta*del2A_resistive
!  A_advecto-resistive=A_resistive+grad Lambda
!
!  25-feb-07/axel: adapted from nospecial.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Special

  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet

  implicit none

  include '../special.h'
!
!  square of wave speed for gauge field
!
  real, pointer :: eta
  logical, pointer :: lweyl_gauge

  ! input parameters
  real :: ampl=1e-3,kx=1.,ky=0.,kz=0.
  logical :: ladvecto_resistive=.true.
  character(len=50) :: init='zero'
  namelist /special_init_pars/ &
    ladvecto_resistive,init,ampl,kx,ky,kz

  ! run parameters
  namelist /special_run_pars/ &
    ladvecto_resistive
!
! Declare any index variables necessary for main or 
! 
   integer :: iLam=0
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_Lamm=0       ! DIAG_DOC: $\left<\Lam\right>$
  integer :: idiag_Lampt=0      ! DIAG_DOC: $\Lam(x1,y1,z1)>$
  integer :: idiag_Lamp2=0      ! DIAG_DOC: $\Lam(x2,y2,z2)>$
  integer :: idiag_Lamrms=0     ! DIAG_DOC: $\left<\Lam^2\right>^{1/2}$
  integer :: idiag_Lambzm=0     ! DIAG_DOC: $\left<\Lam B_z\right>$
  integer :: idiag_Lambzmz=0    ! DIAG_DOC: $\left<\Lam B_z\right>_{xy}$
  integer :: idiag_gLambm=0     ! DIAG_DOC: $\left<\Lam\Bv\right>$
  integer :: idiag_apbrms=0     ! DIAG_DOC: $\left<(\Av'\Bv)^2\right>^{1/2}$
  integer :: idiag_jxarms=0     ! DIAG_DOC: $\left<(\Jv\times\Av)^2\right>^{1/2}$
  integer :: idiag_jxaprms=0    ! DIAG_DOC: $\left<(\Jv\times\Av')^2\right>^{1/2}$
  integer :: idiag_jxgLamrms=0  ! DIAG_DOC: $\left<(\Jv\times\nabla\Lambda)^2\right>^{1/2}$
  integer :: idiag_gLamrms=0    ! DIAG_DOC: $\left<(\nabla\Lambda)^2\right>^{1/2}$
  integer :: idiag_divabrms=0   ! DIAG_DOC: $\left<[(\nabla\cdot\Av)\Bv]^2\right>^{1/2}$
  integer :: idiag_divapbrms=0  ! DIAG_DOC: $\left<[(\nabla\cdot\Av')\Bv]^2\right>^{1/2}$
  integer :: idiag_d2Lambrms=0  ! DIAG_DOC: $\left<[(\nabla^2\Lambda)\Bv]^2\right>^{1/2}$
  integer :: idiag_d2Lamrms=0   ! DIAG_DOC: $\left<[\nabla^2\Lambda)^2\right>^{1/2}$
!
  contains

!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/tony: coded
!
      use FArrayManager
!
      call farray_register_pde('Lam',iLam)
      ispecialvar=iLam
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      use SharedVariables, only : get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
      integer :: ierr
!
!  Initialize module variables which are parameter dependent
!  wave speed of gauge potential
!
      if (.not.lstarting) then
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
      call keep_compiler_quiet(lstarting)
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
        case ('zero'); f(:,:,:,iLam)=0.
        case ('sinwave-x'); call sinwave(ampl,f,iLam,kx=kx)
        case ('sinwave-y'); call sinwave(ampl,f,iLam,ky=ky)
        case ('sinwave-z'); call sinwave(ampl,f,iLam,kz=kz)

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
      if (idiag_apbrms/=0) lpenc_diagnos(i_ab)=.true.
      if (idiag_gLambm/=0) lpenc_diagnos(i_bb)=.true.
      if (idiag_jxarms/=0 .or. idiag_jxaprms/=0 .or. idiag_jxgLamrms/=0) lpenc_diagnos(i_jj)=.true.
      if (idiag_divabrms/=0 .or. idiag_divapbrms/=0 .or. idiag_d2Lambrms/=0) then
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
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   24-nov-04/tony: coded
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
      real, dimension (nx,3) :: gLam,jxa,divab
      real, dimension (nx) :: Lam,del2Lam,ua,ugLam,gLamb,jxa2,divab2
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
!
      if (headtt) call identify_bcs('Lam',iLam)
!
!  solve gauge condition
!
      if (lhydro.or.lhydro_kinematic) then
        call grad(f,iLam,gLam)
        call dot(p%uu,gLam,ugLam)
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
          call del2(f,iLam,del2Lam)
          df(l1:l2,m,n,iLam)=df(l1:l2,m,n,iLam)-ugLam-ua+eta*del2Lam
        else
          df(l1:l2,m,n,iLam)=df(l1:l2,m,n,iLam)-ugLam-ua-eta*p%diva
        endif
      else
        call fatal_error('dspecial_dt','no advective if no hydro')
      endif
!
!  diagnostics
!
      if (ldiagnos) then
        Lam=f(l1:l2,m,n,iLam)
        if (idiag_Lamm/=0) call sum_mn_name(Lam,idiag_Lamm)
        if (idiag_Lamrms/=0) call sum_mn_name(Lam**2,idiag_Lamrms,lsqrt=.true.)
        if (idiag_Lambzm/=0) call sum_mn_name(Lam*p%bb(:,3),idiag_Lambzm)
        if (idiag_gLambm/=0.or.idiag_apbrms/=0) then
          call dot(gLam,p%bb,gLamb)
          if (idiag_gLambm/=0) call sum_mn_name(gLamb,idiag_gLambm)
          if (idiag_apbrms/=0) call sum_mn_name((p%ab+gLamb)**2,idiag_apbrms,lsqrt=.true.)
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
          call cross(p%jj,p%aa+gLam,jxa)
          call dot2(jxa,jxa2)
          if (idiag_jxaprms/=0) call sum_mn_name(jxa2,idiag_jxaprms,lsqrt=.true.)
        endif
!
!  (JxgLam)_rms
!
        if (idiag_jxgLamrms/=0) then
          call cross(p%jj,gLam,jxa)
          call dot2(jxa,jxa2)
          if (idiag_jxgLamrms/=0) call sum_mn_name(jxa2,idiag_jxgLamrms,lsqrt=.true.)
        endif
!
!  (gLam)_rms
!
        if (idiag_jxgLamrms/=0) then
          call dot2(gLam,jxa2)
          if (idiag_gLamrms/=0) call sum_mn_name(jxa2,idiag_gLamrms,lsqrt=.true.)
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
          call multsv(p%diva+del2Lam,p%bb,divab)
          call dot2(divab,divab2)
          if (idiag_divapbrms/=0) call sum_mn_name(divab2,idiag_divapbrms,lsqrt=.true.)
        endif
!
!  (del2Lambda*B)_rms
!
        if (idiag_d2Lambrms/=0) then
          call multsv(del2Lam,p%bb,divab)
          call dot2(divab,divab2)
          if (idiag_d2Lambrms/=0) call sum_mn_name(divab2,idiag_d2Lambrms,lsqrt=.true.)
        endif
!
!  (del2Lambda)_rms
!
        if (idiag_d2Lamrms/=0) then
          if (idiag_d2Lamrms/=0) call sum_mn_name(del2Lam**2,idiag_d2Lamrms,lsqrt=.true.)
        endif
!
!  check for point 1
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_Lampt/=0) call save_name(Lam(lpoint-nghost),idiag_Lampt)
        endif
!
!  check for point 2
!
        if (lroot.and.m==mpoint2.and.n==npoint2) then
          if (idiag_Lamp2/=0) call save_name(Lam(lpoint2-nghost),idiag_Lamp2)
        endif
!
      endif
!
      if (l1davgfirst .or. (ldiagnos .and. ldiagnos_need_zaverages)) then
        if (idiag_Lambzmz/=0)   call xysum_mn_name_z(p%bb(:,3),idiag_Lambzmz)
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
!  read namelist
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
!
!  write name list
!
      write(unit,NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
!  write name list
!
      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif
!
99  endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
!
!  write name list
!
      write(unit,NML=special_run_pars)
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
        idiag_Lamm=0; idiag_Lampt=0; idiag_Lamp2=0; idiag_Lamrms=0
        idiag_gLambm=0; idiag_apbrms=0
        idiag_Lambzm=0; idiag_Lambzmz=0
        idiag_jxarms=0; idiag_divabrms=0
        idiag_jxaprms=0; idiag_divapbrms=0
        idiag_jxgLamrms=0; idiag_d2Lambrms=0
        idiag_gLamrms=0; idiag_d2Lamrms=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'Lamm',idiag_Lamm)
        call parse_name(iname,cname(iname),cform(iname),'Lamrms',idiag_Lamrms)
        call parse_name(iname,cname(iname),cform(iname),'Lambzm',idiag_Lambzm)
        call parse_name(iname,cname(iname),cform(iname),'Lampt',idiag_Lampt)
        call parse_name(iname,cname(iname),cform(iname),'Lamp2',idiag_Lamp2)
        call parse_name(iname,cname(iname),cform(iname),'gLambm',idiag_gLambm)
        call parse_name(iname,cname(iname),cform(iname),'apbrms',idiag_apbrms)
        call parse_name(iname,cname(iname),cform(iname),'jxarms',idiag_jxarms)
        call parse_name(iname,cname(iname),cform(iname),'jxaprms',idiag_jxaprms)
        call parse_name(iname,cname(iname),cform(iname),'jxgLamrms',idiag_jxgLamrms)
        call parse_name(iname,cname(iname),cform(iname),'gLamrms',idiag_gLamrms)
        call parse_name(iname,cname(iname),cform(iname),'divabrms',idiag_divabrms)
        call parse_name(iname,cname(iname),cform(iname),'divapbrms',idiag_divapbrms)
        call parse_name(iname,cname(iname),cform(iname),'d2Lambrms',idiag_d2Lambrms)
        call parse_name(iname,cname(iname),cform(iname),'d2Lamrms',idiag_d2Lamrms)
      enddo
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Lambzmz',idiag_Lambzmz)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_Lamm=',idiag_Lamm
        write(3,*) 'i_Lamrms=',idiag_Lamrms
        write(3,*) 'i_Lambzm=',idiag_Lambzm
        write(3,*) 'i_Lambzmz=',idiag_Lambzmz
        write(3,*) 'i_Lampt=',idiag_Lampt
        write(3,*) 'i_Lamp2=',idiag_Lamp2
        write(3,*) 'i_gLambm=',idiag_gLambm
        write(3,*) 'i_apbrms=',idiag_apbrms
        write(3,*) 'iLam=',iLam
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (slice_data) :: slices
!
      integer :: inamev
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Lam
!
        case ('Lam')
          slices%yz=f(ix_loc,m1:m2,n1:n2,iLam)
          slices%xz=f(l1:l2,iy_loc,n1:n2,iLam)
          slices%xy=f(l1:l2,m1:m2,iz_loc,iLam)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,iLam)
          slices%ready = .true.
!
      endselect
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine calc_lspecial_pars(f)
!
!  dummy routine
!
!  15-jan-08/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lspecial_pars
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_magnetic
!!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine special_after_timestep(f,df,dt_)
!
!   Possibility to modify the f and df after df is updated
!   Used for the fargo shift, for instance.
!
!   27-nov-08/wlad: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dt_
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dt_)
!
    endsubroutine  special_after_timestep
!********************************************************************
    subroutine special_calc_particles(fp)
!
!   Called before the loop, in case some particle value is needed 
!   for the special density/hydro/magnetic/entropy
!
!   20-nov-08/wlad: coded
!
      real, dimension (:,:), intent(in) :: fp
!
      call keep_compiler_quiet(fp)
!
    endsubroutine special_calc_particles
!***********************************************************************
    subroutine special_calc_particles_nbody(fsp)
!
!   Called before the loop, in case some massive particles value 
!   is needed for the special density/hydro/magnetic/entropy
!
!   20-nov-08/wlad: coded
!
      real, dimension (:,:), intent(in) :: fsp
!
      call keep_compiler_quiet(fsp)
!
!
    endsubroutine special_calc_particles_nbody
!***********************************************************************
    subroutine special_boundconds(f,bc)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (boundary_condition) :: bc
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(bc)
!
    endsubroutine special_boundconds
!***********************************************************************
    subroutine special_before_boundary(f)
!
!   Possibility to modify the f array before the boundaries are
!   communicated.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_before_boundary
!***********************************************************************

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

