! $Id$
!
!  This modules solves the passive scalar advection equation.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lpscalar = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED cc; cc1; lncc; glncc(3); uglncc; del2lncc
! PENCILS PROVIDED hlncc(3,3)
!
!***************************************************************
module Pscalar
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'pscalar.h'
!
  character (len=labellen) :: initlncc='zero', initlncc2='zero'
  character (len=40) :: tensor_pscalar_file
  logical :: nopscalar=.false.
  real, dimension(3) :: gradC0=(/0.0,0.0,0.0/)
  real :: ampllncc=0.1, widthlncc=0.5, cc_min=0.0, lncc_min
  real :: ampllncc2=0.0, kx_lncc=1.0, ky_lncc=1.0, kz_lncc=1.0, radius_lncc=0.0
  real :: epsilon_lncc=0.0, cc_const=0.0
  logical :: lupw_lncc=.false.
  logical :: ldustdrift=.false.
!
  namelist /pscalar_init_pars/ &
      initlncc,initlncc2,ampllncc,ampllncc2,kx_lncc,ky_lncc,kz_lncc, &
      radius_lncc,epsilon_lncc,widthlncc,cc_min,cc_const,lupw_lncc, ldustdrift
!
  real :: pscalar_diff=0.0, tensor_pscalar_diff=0.0
  real :: rhoccm=0.0, cc2m=0.0, gcc2m=0.0
!
  namelist /pscalar_run_pars/ &
      pscalar_diff,nopscalar,tensor_pscalar_diff,gradC0,lupw_lncc
!
  integer :: idiag_rhoccm=0     ! DIAG_DOC: $\left<\varrho c\right>$
  integer :: idiag_ccmax=0      ! DIAG_DOC: $\max(c)$
  integer :: idiag_ccmin=0      ! DIAG_DOC:
  integer :: idiag_lnccm=0      ! DIAG_DOC:
  integer :: idiag_mcct=0       ! DIAG_DOC:
  integer :: idiag_gcc5m=0      ! DIAG_DOC:
  integer :: idiag_gcc10m=0     ! DIAG_DOC:
  integer :: idiag_ucm=0        ! DIAG_DOC:
  integer :: idiag_uudcm=0      ! DIAG_DOC:
  integer :: idiag_Cz2m=0       ! DIAG_DOC:
  integer :: idiag_Cz4m=0       ! DIAG_DOC:
  integer :: idiag_cc1m=0       ! DIAG_DOC:
  integer :: idiag_cc2m=0       ! DIAG_DOC:
  integer :: idiag_cc3m=0       ! DIAG_DOC:
  integer :: idiag_cc4m=0       ! DIAG_DOC:
  integer :: idiag_cc5m=0       ! DIAG_DOC:
  integer :: idiag_cc6m=0       ! DIAG_DOC:
  integer :: idiag_cc7m=0       ! DIAG_DOC:
  integer :: idiag_cc8m=0       ! DIAG_DOC:
  integer :: idiag_cc9m=0       ! DIAG_DOC:
  integer :: idiag_cc10m=0      ! DIAG_DOC:
  integer :: idiag_Crmsm=0      ! DIAG_DOC:
  integer :: idiag_gcc1m=0      ! DIAG_DOC:
  integer :: idiag_gcc2m=0      ! DIAG_DOC:
  integer :: idiag_gcc3m=0      ! DIAG_DOC:
  integer :: idiag_gcc4m=0      ! DIAG_DOC:
  integer :: idiag_gcc6m=0      ! DIAG_DOC:
  integer :: idiag_gcc7m=0      ! DIAG_DOC:
  integer :: idiag_gcc8m=0      ! DIAG_DOC:
  integer :: idiag_gcc9m=0      ! DIAG_DOC:
  integer :: idiag_lnccmz=0     ! DIAG_DOC:
  integer :: idiag_lnccmy=0     ! DIAG_DOC:
  integer :: idiag_lnccmx=0     ! DIAG_DOC:
  integer :: idiag_ccmz=0       ! DIAG_DOC:
  integer :: idiag_ccmy=0       ! DIAG_DOC:
  integer :: idiag_ccmx=0       ! DIAG_DOC:
  integer :: idiag_ccglnrm=0    ! DIAG_DOC: $\left<c\nabla_z\varrho\right>$
!
  contains
!***********************************************************************
    subroutine register_pscalar()
!
!  Initialise variables which should know that we solve for passive
!  scalar: ilncc; increase nvar accordingly.
!
!  6-jul-02/axel: coded
!
      use FArrayManager
!
      call farray_register_pde('lncc',ilncc)
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',lncc $'
          if (nvar == mvar) write(4,*) ',lncc'
        else
          write(4,*) ',lncc $'
        endif
        write(15,*) 'lncc = fltarr(mx,my,mz)*one'
      endif
!
    endsubroutine register_pscalar
!***********************************************************************
    subroutine initialize_pscalar(f)
!
!  Perform any necessary post-parameter read initialization.
!  Dummy routine
!
!  24-nov-02/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_pscalar
!***********************************************************************
    subroutine init_lncc(f)
!
!  Initialise passive scalar field; called from start.f90.
!
!   6-jul-2001/axel: coded
!
      use Sub
      use Initcond
      use InitialCondition, only: initial_condition_lncc
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      select case (initlncc)
        case ('zero'); f(:,:,:,ilncc)=0.
        case ('constant'); f(:,:,:,ilncc)=log(cc_const)
        case ('gaussian-x'); call gaussian(ampllncc,f,ilncc,kx=kx_lncc)
        case ('gaussian-y'); call gaussian(ampllncc,f,ilncc,ky=ky_lncc)
        case ('gaussian-z'); call gaussian(ampllncc,f,ilncc,kz=kz_lncc)
        case ('parabola-x'); call parabola(ampllncc,f,ilncc,kx=kx_lncc)
        case ('parabola-y'); call parabola(ampllncc,f,ilncc,ky=ky_lncc)
        case ('parabola-z'); call parabola(ampllncc,f,ilncc,kz=kz_lncc)
        case ('gaussian-noise'); call gaunoise(ampllncc,f,ilncc,ilncc)
        case ('wave-x'); call wave(ampllncc,f,ilncc,kx=kx_lncc)
        case ('wave-y'); call wave(ampllncc,f,ilncc,ky=ky_lncc)
        case ('wave-z'); call wave(ampllncc,f,ilncc,kz=kz_lncc)
        case ('propto-ux'); call wave_uu(ampllncc,f,ilncc,kx=kx_lncc)
        case ('propto-uy'); call wave_uu(ampllncc,f,ilncc,ky=ky_lncc)
        case ('propto-uz'); call wave_uu(ampllncc,f,ilncc,kz=kz_lncc)
        case ('tang-discont-z')
          print*,'init_lncc: widthlncc=',widthlncc
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilncc)=-1.0+2*0.5*(1.+tanh(z(n)/widthlncc))
          enddo; enddo
        case ('hor-tube'); call htube2(ampllncc,f,ilncc,ilncc,radius_lncc,epsilon_lncc)
        case default; call fatal_error('init_lncc','bad initlncc='//trim(initlncc))
      endselect
!
!  Superimpose something else.
!
      select case (initlncc2)
        case ('wave-x'); call wave(ampllncc2,f,ilncc,ky=5.)
      endselect
!
!  Interface for user's own initial condition.
!
      if (linitial_condition) call initial_condition_lncc(f)
!
!  Add floor value if cc_min is set
!
      if (cc_min/=0.) then
        lncc_min=log(cc_min)
        if (lroot) print*,'set floor value for cc; cc_min=',cc_min
        f(:,:,:,ilncc)=max(lncc_min,f(:,:,:,ilncc))
      endif
!
    endsubroutine init_lncc
!***********************************************************************
    subroutine pencil_criteria_pscalar()
!
!  All pencils that the Pscalar module depends on are specified here.
!
!  20-11-04/anders: coded
!
      integer :: i
!
      if (lhydro) lpenc_requested(i_uglncc)=.true.
      if (pscalar_diff/=0.) then
        lpenc_requested(i_glncc)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_del2lncc)=.true.
      endif
      do i=1,3
        if (gradC0(i)/=0.) lpenc_requested(i_uu)=.true.
      enddo
      if (tensor_pscalar_diff/=0.0) then
        lpenc_requested(i_hlncc)=.true.
        lpenc_requested(i_glncc)=.true.
      endif
!
      if (idiag_rhoccm/=0 .or. idiag_ccmax/=0 .or. idiag_ccmin/=0 .or. &
          idiag_ucm/=0 .or. idiag_uudcm/=0 .or. idiag_Cz2m/=0 .or. &
          idiag_Cz4m/=0 .or. idiag_Crmsm/=0 .or. idiag_rhoccm/=0 .or. &
          idiag_mcct/=0) lpenc_diagnos(i_cc)=.true.
      if (idiag_rhoccm/=0 .or. idiag_Cz2m/=0 .or. idiag_Cz4m/=0 .or. &
          idiag_Crmsm/=0 .or. idiag_mcct/=0) lpenc_diagnos(i_rho)=.true.
      if (idiag_ucm/=0 .or. idiag_uudcm/=0) lpenc_diagnos(i_uu)=.true.
      if (idiag_uudcm/=0) lpenc_diagnos(i_uglncc)=.true.
      if (idiag_lnccmz/=0 .or. idiag_lnccmy/=0 .or. idiag_lnccmx/=0) &
          lpenc_diagnos(i_lncc)=.true.
      if (idiag_ccmz/=0 .or. idiag_ccmy/=0 .or. idiag_ccmx/=0) &
          lpenc_diagnos(i_cc)=.true.
      if (idiag_ccglnrm/=0) lpenc_requested(i_glnrho)=.true.
!
    endsubroutine pencil_criteria_pscalar
!***********************************************************************
    subroutine pencil_interdep_pscalar(lpencil_in)
!
!  Interdependency among pencils provided by the Pscalar module
!  is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_cc1)) lpencil_in(i_cc)=.true.
      if (lpencil_in(i_uglncc)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_glncc)=.true.
      endif
!
    endsubroutine pencil_interdep_pscalar
!***********************************************************************
    subroutine calc_pencils_pscalar(f,p)
!
!  Calculate Pscalar pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
! cc
      if (lpencil(i_cc)) p%cc=exp(f(l1:l2,m,n,ilncc))
! cc1
      if (lpencil(i_cc1)) p%cc1=1/p%cc
! glncc
      if (lpencil(i_glncc)) call grad(f,ilncc,p%glncc)
! uglncc
      if (lpencil(i_uglncc)) then
        if (ldustdrift) then
          call u_dot_grad(f,ilncc,p%glncc,p%udropav,p%uglncc,UPWIND=lupw_lncc)
        else
          call u_dot_grad(f,ilncc,p%glncc,p%uu,p%uglncc,UPWIND=lupw_lncc)
        endif
      endif
! del2lncc
      if (lpencil(i_del2lncc)) call del2(f,ilncc,p%del2lncc)
! hlncc
      if (lpencil(i_hlncc)) call g2ij(f,ilncc,p%hlncc)
!
    endsubroutine calc_pencils_pscalar
!***********************************************************************
    subroutine dlncc_dt(f,df,p)
!
!  Passive scalar evolution.
!
!  Calculate dc/dt=-uu.glncc + pscaler_diff*[del2lncc + (glncc+glnrho).glncc]
!
!   6-jul-02/axel: coded
!
      use Diagnostics
      use Special, only: special_calc_pscalar
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: diff_op
      integer :: j
!
      intent(in)  :: f
      intent(out) :: df
!
!  Identify module and boundary conditions.
!
      if (nopscalar) then
        if (headtt.or.ldebug) print*,'not SOLVED: dlncc_dt'
      else
        if (headtt.or.ldebug) print*,'SOLVE dlncc_dt'
      endif
      if (headtt) call identify_bcs('cc',ilncc)
!
!  Gradient of passive scalar.
!  Allow for possibility to turn off passive scalar
!  without changing file size and recompiling everything.
!
      if (.not. nopscalar) then ! i.e. if (pscalar)
!
!  Passive scalar equation.
!
        if (lhydro) df(l1:l2,m,n,ilncc) = df(l1:l2,m,n,ilncc) - p%uglncc
!
!  Diffusion operator.
!
        if (pscalar_diff/=0.) then
          if (headtt) print*,'dlncc_dt: pscalar_diff=',pscalar_diff
          call dot_mn(p%glncc+p%glnrho,p%glncc,diff_op)
          diff_op=diff_op+p%del2lncc
          df(l1:l2,m,n,ilncc) = df(l1:l2,m,n,ilncc) + pscalar_diff*diff_op
        endif
!
!  Add advection of imposed constant gradient of lncc (called gradC0).
!  Makes sense really only for periodic boundary conditions.
!  This gradient can have arbitary direction.
!
        do j=1,3
          if (gradC0(j)/=0.) then
            df(l1:l2,m,n,ilncc) = df(l1:l2,m,n,ilncc) - gradC0(j)*p%uu(:,j)
          endif
        enddo
!
!  Tensor diffusion (but keep the isotropic one).
!
        if (tensor_pscalar_diff/=0.) &
            call tensor_diff(df,p,tensor_pscalar_diff)
!
        if (lspecial) call special_calc_pscalar(f,df,p)
!
      endif
!
!  Diagnostics.
!
!  output for double and triple correlators (assume z-gradient of cc)
!  <u_k u_j d_j c> = <u_k c uu.gradlncc>
!
      if (ldiagnos) then
        if (idiag_mcct/=0)   call integrate_mn_name(p%rho*p%cc,idiag_mcct)
        if (idiag_rhoccm/=0) call sum_mn_name(p%rho*p%cc,idiag_rhoccm)
        if (idiag_ccmax/=0)  call max_mn_name(p%cc,idiag_ccmax)
        if (idiag_ccmin/=0)  call max_mn_name(-p%cc,idiag_ccmin,lneg=.true.)
        if (idiag_ucm/=0)    call sum_mn_name(p%uu(:,3)*p%cc,idiag_ucm)
        if (idiag_uudcm/=0)  &
            call sum_mn_name(p%uu(:,3)*p%cc*p%uglncc,idiag_uudcm)
        if (idiag_Cz2m/=0)   call sum_mn_name(p%rho*p%cc*z(n)**2,idiag_Cz2m)
        if (idiag_Cz4m/=0)   call sum_mn_name(p%rho*p%cc*z(n)**4,idiag_Cz4m)
        if (idiag_Crmsm/=0)  &
            call sum_mn_name((p%rho*p%cc)**2,idiag_Crmsm,lsqrt=.true.)
        if (idiag_ccglnrm/=0) call sum_mn_name(p%cc*p%glnrho(:,3),idiag_ccglnrm)
      endif
!
      if (l1davgfirst) then
         if (idiag_lnccmz/=0) call xysum_mn_name_z(p%lncc,idiag_lnccmz)
         if (idiag_lnccmy/=0) call xzsum_mn_name_y(p%lncc,idiag_lnccmy)
         if (idiag_lnccmx/=0) call yzsum_mn_name_x(p%lncc,idiag_lnccmx)
         if (idiag_ccmz/=0)   call xysum_mn_name_z(p%cc,idiag_ccmz)
         if (idiag_ccmy/=0)   call xzsum_mn_name_y(p%cc,idiag_ccmy)
         if (idiag_ccmx/=0)   call yzsum_mn_name_x(p%cc,idiag_ccmx)
      endif
!
    endsubroutine dlncc_dt
!***********************************************************************
    subroutine read_pscalar_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=pscalar_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=pscalar_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_pscalar_init_pars
!***********************************************************************
    subroutine write_pscalar_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=pscalar_init_pars)
!
    endsubroutine write_pscalar_init_pars
!***********************************************************************
    subroutine read_pscalar_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=pscalar_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=pscalar_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_pscalar_run_pars
!***********************************************************************
    subroutine write_pscalar_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=pscalar_run_pars)
!
    endsubroutine write_pscalar_run_pars
!***********************************************************************
    subroutine rprint_pscalar(lreset,lwrite)
!
!  Reads and registers print parameters relevant for passive scalar.
!
!   6-jul-02/axel: coded
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamez, inamey, inamex
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_rhoccm=0; idiag_ccmax=0; idiag_ccmin=0.; idiag_lnccm=0
        idiag_ucm=0; idiag_uudcm=0; idiag_Cz2m=0; idiag_Cz4m=0
        idiag_Crmsm=0; idiag_mcct=0
        idiag_lnccmz=0; idiag_lnccmy=0; idiag_lnccmx=0
        idiag_ccmz=0; idiag_ccmy=0; idiag_ccmx=0; idiag_ccglnrm=0
      endif
!
!  Check for those quantities that we want to evaluate online.
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'mcct',idiag_mcct)
        call parse_name(iname,cname(iname),cform(iname),'rhoccm',idiag_rhoccm)
        call parse_name(iname,cname(iname),cform(iname),'ccmax',idiag_ccmax)
        call parse_name(iname,cname(iname),cform(iname),'ccmin',idiag_ccmin)
        call parse_name(iname,cname(iname),cform(iname),'lnccm',idiag_lnccm)
        call parse_name(iname,cname(iname),cform(iname),'ucm',idiag_ucm)
        call parse_name(iname,cname(iname),cform(iname),'uudcm',idiag_uudcm)
        call parse_name(iname,cname(iname),cform(iname),'Cz2m',idiag_Cz2m)
        call parse_name(iname,cname(iname),cform(iname),'Cz4m',idiag_Cz4m)
        call parse_name(iname,cname(iname),cform(iname),'Crmsm',idiag_Crmsm)
        call parse_name(iname,cname(iname),cform(iname),'ccglnrm',idiag_ccglnrm)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'lnccmz',idiag_lnccmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'ccmz',idiag_ccmz)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'lnccmy',idiag_lnccmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'ccmy',idiag_ccmy)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'lnccmx',idiag_lnccmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'ccmx',idiag_ccmx)
      enddo
!
!  Write column where which passive scalar variable is stored.
!
      if (lwr) then
        write(3,*) 'ilncc=',ilncc
        write(3,*) 'icc=0'
      endif
!
    endsubroutine rprint_pscalar
!***********************************************************************
    subroutine get_slices_pscalar(f,slices)
!
!  Write slices for animation of pscalar variables.
!
!  26-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices.
!
      select case (trim(slices%name))
!
!  Passive scalar.
!
        case ('cc')
          slices%yz =exp(f(ix_loc,m1:m2,n1:n2,ilncc))
          slices%xz =exp(f(l1:l2,iy_loc,n1:n2,ilncc))
          slices%xy =exp(f(l1:l2,m1:m2,iz_loc,ilncc))
          slices%xy2=exp(f(l1:l2,m1:m2,iz2_loc,ilncc))
          if (lwrite_slice_xy3) slices%xy3=exp(f(l1:l2,m1:m2,iz3_loc,ilncc))
          if (lwrite_slice_xy4) slices%xy4=exp(f(l1:l2,m1:m2,iz4_loc,ilncc))
          slices%ready=.true.
!
!  Logarithmic passive scalar.
!
        case ('lncc')
          slices%yz =f(ix_loc,m1:m2,n1:n2,ilncc)
          slices%xz =f(l1:l2,iy_loc,n1:n2,ilncc)
          slices%xy =f(l1:l2,m1:m2,iz_loc,ilncc)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,ilncc)
          if (lwrite_slice_xy3) slices%xy3=f(l1:l2,m1:m2,iz3_loc,ilncc)
          if (lwrite_slice_xy4) slices%xy4=f(l1:l2,m1:m2,iz4_loc,ilncc)
          slices%ready=.true.
!
      endselect
!
    endsubroutine get_slices_pscalar
!***********************************************************************
    subroutine calc_mpscalar
!
!  calculate mean magnetic field from xy- or z-averages
!
!  14-apr-03/axel: adaped from calc_mfield
!
      use Diagnostics
!
      logical,save :: first=.true.
      real :: lnccm
!
!  Magnetic energy in horizontally averaged field
!  The bxmz and bymz must have been calculated,
!  so they are present on the root processor.
!
      if (idiag_lnccm/=0) then
        if (idiag_lnccmz==0) then
          if (first) print*
          if (first) print*,"NOTE: to get lnccm, lnccmz must also be set in xyaver"
          if (first) print*,"      We proceed, but you'll get lnccm=0"
          lnccm=0.
        else
          lnccm=sqrt(sum(fnamez(:,:,idiag_lnccmz)**2)/(nz*nprocz))
        endif
        call save_name(lnccm,idiag_lnccm)
      endif
!
    endsubroutine calc_mpscalar
!***********************************************************************
    subroutine tensor_diff(df,p,tensor_pscalar_diff)
!
!  reads file
!
!  11-jul-02/axel: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: tensor_pscalar_diff
!
      real, save, dimension (nx,ny,nz,3) :: bunit,hhh
      real, dimension (nx) :: tmp,scr
      integer :: iy,iz,i,j
      logical, save :: first=.true.
!
!  read H and Bunit arrays and keep them in memory
!
      if (first) then
        open(1,file=trim(directory)//'/bunit.dat',form='unformatted')
        print*,'read bunit.dat with dimension: ',nx,ny,nz,3
        read(1) bunit,hhh
        close(1)
        print*,'read bunit.dat; bunit=',bunit
      endif
!
!  tmp = (Bunit.G)^2 + H.G + Bi*Bj*Gij
!  for details, see tex/mhd/thcond/tensor_der.tex
!
      call dot_mn(bunit,p%glncc,scr)
      call dot_mn(hhh,p%glncc,tmp)
      tmp=tmp+scr**2
!
!  dot with bi*bj
!
      iy=m-m1+1
      iz=n-n1+1
      do j=1,3
      do i=1,3
        tmp=tmp+bunit(:,iy,iz,i)*bunit(:,iy,iz,j)*p%hlncc(:,i,j)
      enddo
      enddo
!
!  and add result to the dlncc/dt equation
!
      df(l1:l2,m,n,ilncc)=df(l1:l2,m,n,ilncc)+tensor_pscalar_diff*tmp
!
      first=.false.
!
    endsubroutine tensor_diff
!***********************************************************************
endmodule Pscalar
