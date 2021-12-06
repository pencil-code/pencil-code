! $Id$
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
  use Messages
!
  implicit none
!
  include '../special.h'
!
  !
  ! Slice precalculation buffers
  !
  real, target, dimension(:,:,:), allocatable :: oo_xy_meanx,oo_xy2_meanx,oo_xy3_meanx,oo_xy4_meanx
  real, target, dimension(:,:,:), allocatable :: oo_xz_meanx,oo_xz2_meanx,oo_yz_meanx
  real, target, dimension(:,:,:,:,:,:), allocatable :: oo_r_meanx
  real, target, dimension(:,:,:), allocatable :: uu_xy_meanx,uu_xy2_meanx,uu_xy3_meanx,uu_xy4_meanx
  real, target, dimension(:,:,:), allocatable :: uu_xz_meanx,uu_xz2_meanx,uu_yz_meanx
  real, target, dimension(:,:,:,:,:,:), allocatable :: uu_r_meanx

  real, dimension(nygrid,3) :: mean_u
!
  character(len=24) :: initspecial='nothing'
  real :: central_vel=0,ampluu_spec=0,Re_tau=180
!
!!  character, len(50) :: initcustom
!
! input parameters
  namelist /internal_flow_init_pars/ &
       initspecial,central_vel,ampluu_spec,Re_tau
!
  integer :: idiag_turbint=0
  integer :: idiag_uxm_central=0,idiag_tau_w=0
  integer :: ivid_uu_meanx=0, ivid_oo_meanx=0
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
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      use Slices_methods, only: alloc_slice_buffers
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Initialize any module variables which are parameter dependent
!
      call keep_compiler_quiet(f)
!
      if (ivid_uu_meanx/=0) &
        call alloc_slice_buffers(uu_xy_meanx,uu_xz_meanx,uu_yz_meanx, &
                                 uu_xy2_meanx,uu_xy3_meanx,uu_xy4_meanx,uu_xz2_meanx,uu_r_meanx)

      if (ivid_oo_meanx/=0) &
        call alloc_slice_buffers(oo_xy_meanx,oo_xz_meanx,oo_yz_meanx, &
                                 oo_xy2_meanx,oo_xy3_meanx,oo_xy4_meanx,oo_xz2_meanx,oo_r_meanx)

    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Mpicomm
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j
      real :: height,h2
!
      intent(inout) :: f
      !
      ! Select case
      !
      select case (initspecial)
      case ('nothing'); if (lroot) print*,'init_special: nothing'
      case ('poiseulle_xy')
        f(:,:,:,iux:iuz)=0
        call poiseulle_flowx_wally(f,central_vel)
      case ('poiseulle_xy_noise')
        f(:,:,:,iux:iuz)=0
        height=Lxyz(2)/2
        h2=height**2
        call gaunoise(ampluu_spec,f,iux,iuz)
        do j=m1,m2
          do i=iux,iuz
            f(l1:l2,j,n1:n2,i)=f(l1:l2,j,n1:n2,i)&
               *(1-(y(j)-xyz0(2)-height)**2/h2)
          enddo
        enddo
        call poiseulle_flowx_wally(f,central_vel)
      case ('velocity_defect_xy')
        call velocity_defect_flowx_wally(f,central_vel,Re_tau)
      case ('log_law_xy')
        call log_law_flowx_wally(f,central_vel,Re_tau)
      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*,'init_special: No such value for initspecial: ', &
             trim(initspecial)
        call stop_it("")
      endselect
      !
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
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
      use Deriv, only: der_pencil
      use Slices_methods, only: store_slices
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: ufluct
      real, dimension (nx) :: ufluct2
      real, dimension (3) :: meanx_uu
      type (pencil_case) :: p
      integer :: i,j
      real, dimension (my) :: tmp,du_mean_dy
      real :: tau_tmp
      real, pointer :: nu
!
      intent(in) :: f,p
      intent(inout) :: df
!
      call get_shared_variable('nu',nu,caller='dspecial_dt')
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
      !
      ! Write video slices
      !
      if (lvideo.and.lfirst) then
        if (ivid_uu_meanx/=0) then
          do j=1,3
            !ufluct(:,j)=p%uu(:,j)-mean_u(m+ny*ipy_meanx-nghost,j)   !!!MR: corrected, please check!
            ufluct(:,j)=p%uu(:,j)-mean_u(m+ny*ipy-nghost,j)
          enddo
          call store_slices(ufluct,uu_xy_meanx,uu_xz_meanx,uu_yz_meanx,uu_xy2_meanx,uu_xy3_meanx, &
                            uu_xy4_meanx,uu_xz2_meanx,uu_r_meanx)
        endif
        if (ivid_oo_meanx/=0) then
          do j=1,3
            ufluct(:,j)=p%oo(:,j)-sum(p%oo(:,j))/nx
          enddo
          call store_slices(ufluct,oo_xy_meanx,oo_xz_meanx,oo_yz_meanx,oo_xy2_meanx,oo_xy3_meanx, &
                            oo_xy4_meanx,oo_xz2_meanx,oo_r_meanx)
        endif
      endif
!
! Write diagnostics
!
      if (ldiagnos) then
        if (idiag_turbint/=0) then
          do j=1,3
            meanx_uu(j)=mean_u(m+ny*ipy-nghost,j)
            ufluct(:,j)=p%uu(:,j)-meanx_uu(j)
          enddo
          call dot2(ufluct,ufluct2)
          call sum_mn_name(ufluct2,idiag_turbint)
        endif
        if (idiag_uxm_central/=0) then
          if (m==m2.and.n==n1) then
            fname(idiag_uxm_central)=maxval(mean_u)
            itype_name(idiag_uxm_central)=ilabel_max
          endif
!          call max_mn_name(meanx_uu(1),idiag_uxm_central)
        endif
        if (idiag_tau_w/=0) then
          if (m==m1.and.n==n1) then
            if (ipy==0) then
              tmp=0
              tmp(m1:m2)=mean_u(1:ny,1)
              call der_pencil(2,tmp,du_mean_dy)
              if (lroot) print*,'rhom is hardcoded in internal_flow.f90: dspecial_dt'
              tau_tmp=du_mean_dy(m1+3)*nu*1.2
!              tau_tmp=-(mean_u(2,1)-mean_u(1,1))/(y(l1+1)-y(l1+0))
            else
              tau_tmp=0
            endif
            !call max_mn_name(tau_tmp,idiag_tau_w)
            itype_name(idiag_tau_w)=ilabel_max
            fname(idiag_tau_w)=tau_tmp
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
      read(parallel_unit, NML=internal_flow_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=internal_flow_init_pars)
!
    endsubroutine write_special_init_pars
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
        idiag_turbint=0
        idiag_tau_w=0
        idiag_uxm_central=0
        ivid_uu_meanx=0; ivid_oo_meanx=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'turbint',idiag_turbint)
        call parse_name(iname,cname(iname),cform(iname),'tau_w',idiag_tau_w)
        call parse_name(iname,cname(iname),cform(iname),'uxm_central',idiag_uxm_central)
      enddo
!
!  check for those quantities for which we want video slices
!
      do iname=1,nnamev
        call parse_name(iname,cnamev(iname),cformv(iname),'uu_meanx', ivid_uu_meanx)
        call parse_name(iname,cnamev(iname),cformv(iname),'oo_meanx', ivid_oo_meanx)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        call farray_index_append('i_turbint',idiag_turbint)
        call farray_index_append('i_tau_w',idiag_tau_w)
        call farray_index_append('i_uxm_central',idiag_uxm_central)
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of special variables.
!
!  26-jun-06/tony: dummy
!
      use Slices_methods, only: assign_slices_vec

      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
      !
      !  Loop over slices
      !
      select case (trim(slices%name))
        !
        !  Mean vorticity (derived variable)
        !
      case ('oo_meanx')
        call assign_slices_vec(slices,oo_xy_meanx,oo_xz_meanx,oo_yz_meanx,oo_xy2_meanx, &
                                      oo_xy3_meanx,oo_xy4_meanx,oo_xz2_meanx,oo_r_meanx)
      case ('uu_meanx')
        call assign_slices_vec(slices,uu_xy_meanx,uu_xz_meanx,uu_yz_meanx,uu_xy2_meanx, &
                                      uu_xy3_meanx,uu_xy4_meanx,uu_xz2_meanx,uu_r_meanx)
      endselect
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Mean flow velocitites
!
!  14-mar-08/nils: coded
!
      use Sub
      use Mpicomm, only: mpireduce_sum, mpibcast_real
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nygrid,3) :: mean_u_tmp
      real :: faq
      integer :: j,k
!
!  calculate mean of velocity in xz planes
!
      if (lvideo.and.lfirst .or. ldiagnos) then
        mean_u_tmp=0
        faq=nxgrid*nzgrid
        do j=m1,m2
          do k=1,3
            mean_u_tmp(j+ny*ipy-nghost,k)=sum(f(l1:l2,j,n1:n2,k+iux-1))/faq
          enddo
        enddo
        do k=1,3
          call mpireduce_sum(mean_u_tmp(:,k),mean_u(:,k),nygrid)
          call mpibcast_real(mean_u(:,k),nygrid)
        enddo
      endif
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine special_boundconds(f,bc)
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   2008-06-19/nils: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (boundary_condition) :: bc
!
      select case (bc%bcname)
      case ('poi')
        select case (bc%location)
        case (iBC_X_TOP)
          call bc_poi_x(f,-1,'top',iux,REL=.true.,val=bc%value1)
        case (iBC_X_BOT)
          call bc_poi_x(f,-1,'bot',iux,REL=.true.,val=bc%value1)
        end select
      end select
!
    endsubroutine special_boundconds
!***********************************************************************
    subroutine poiseulle_flowx_wally(f,central_vel)
      !
      ! Set initial Poiseulle flow in x-direction.
      ! The walls are in the y-direction
      !
      ! 2008.02.18: Nils Erland (Coded)
      !
      real, dimension (mx,my,mz,mfarray) :: f
      real :: central_vel, height, h2
      integer :: i,j,k
      !
      height=Lxyz(2)/2
      h2=height**2
      !
      do j=m1,m2
        f(l1:l2,j,n1:n2,iux)=f(l1:l2,j,n1:n2,iux)+central_vel*&
            (1-(y(j)-xyz0(2)-height)**2/h2)
    enddo
      !
    endsubroutine poiseulle_flowx_wally
!***********************************************************************
    subroutine velocity_defect_flowx_wally(f,central_vel,Re_tau)
      !
      ! Set initial turbulent flow in x-direction.
      ! The walls are in the y-direction.
      ! This method is based on the velocity defect law of Milikan (1938).
      ! More data on this can be found in e.g. Pope's book on turbulent flows.
      !
      ! 2008.03.22: Nils Erland (Coded)
      !
      use Mpicomm, only: ipy
      use SharedVariables, only: get_shared_variable
      !
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, intent(in) :: central_vel,Re_tau
      real :: B1,kappa,utau,height, h2
      real, pointer :: nu
      integer :: i,j,k
!
      real :: y_pluss, defect,def_y,log_y,lw,defect_min,log_max
!
      call get_shared_variable('nu',nu,caller='dspecial_dt')
!
if (nu==0) then
  print*,'Nu is zero, setting it to 1.5e-5 for now, but this should be fixed!'
  nu=1.5e-5
endif
      height=Lxyz(2)/2
      h2=height**2
      B1=0.2
      kappa=0.41
      utau=Re_tau*nu/height
      lw=nu/utau
      log_y=4
      def_y=30
!print*,'Re_tau=',Re_tau
!print*,'height=',height
!print*,'u_tau=',utau
!print*,'lw=',lw
      !
      ! Set some interpolation parameters
      !
      defect_min=central_vel&
           +utau*log(def_y*lw/height)/kappa&
           -B1*utau
      log_max=log_y
      !
      ! Add mean turbulent velocity profile to the possibly already
      ! existing turnulent velocity field. As there should be less turbulence
      ! close to the walls we scale the existing turbulence field with the
      ! velocity profile.
      !
      do j=m1,m2
        if (y(j)<xyz0(2)+height) then
          y_pluss=(y(j)-xyz0(2))/lw
          defect=central_vel&
               +utau*log(y_pluss*lw/height)/kappa&
               -B1*utau
          if (y_pluss .le. log_y) then
            f(l1:l2,j,n1:n2,iux)=&
                 (f(l1:l2,j,n1:n2,iux)/central_vel+1)*y_pluss*utau
          elseif (y_pluss .ge. def_y) then
            f(l1:l2,j,n1:n2,iux)=(f(l1:l2,j,n1:n2,iux)/central_vel+1)*defect
          else
            f(l1:l2,j,n1:n2,iux)=(f(l1:l2,j,n1:n2,iux)/central_vel+1)*&
                 (defect_min*(y_pluss-log_y)/(def_y-log_y) &
                 +log_max*utau*(def_y-y_pluss)/(def_y-log_y))
          endif
        else
          y_pluss=-(y(j)-xyz0(2)-Lxyz(2))/lw
          defect=central_vel&
               +utau*log(y_pluss*lw/height)/kappa&
               -B1*utau
          if (y_pluss .le. log_y) then
            f(l1:l2,j,n1:n2,iux)=(f(l1:l2,j,n1:n2,iux)/central_vel+1)*y_pluss*utau
          elseif (y_pluss .ge. def_y) then
            f(l1:l2,j,n1:n2,iux)=(f(l1:l2,j,n1:n2,iux)/central_vel+1)*defect
          else
            f(l1:l2,j,n1:n2,iux)=(f(l1:l2,j,n1:n2,iux)/central_vel+1)*&
                 (defect_min*(y_pluss-log_y)/(def_y-log_y) &
                 +log_max*utau*(def_y-y_pluss)/(def_y-log_y))
          endif
        endif
        !
      enddo
      if (lfirst_proc_y)  f(l1:l2,m1,n1:n2,iux)=0.0
      if (llast_proc_y) f(l1:l2,m2,n1:n2,iux)=0.0
      !
    endsubroutine velocity_defect_flowx_wally
!***********************************************************************
    subroutine log_law_flowx_wally(f,central_vel,Re_tau)
      !
      ! Set initial turbulent flow in x-direction.
      ! The walls are in the y-direction.
      ! This method is based on log-law.
      ! More data on this can be found in e.g. Pope's book on turbulent flows.
      !
      ! 2008.05.20: Nils Erland (Coded)
      !
      use Mpicomm, only: ipy
      use SharedVariables, only: get_shared_variable
      !
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, intent(in) :: central_vel,Re_tau
      real :: B1,kappa,utau,height, h2
      real, pointer :: nu
      integer :: i,j,k
!
      real :: y_pluss, lw, u_log, u_lam
      !
      call get_shared_variable('nu',nu,caller='dspecial_dt')
!
if (nu==0) then
  print*,'WARNING: Nu is zero, setting it to 1.5e-5 for now, but this should be fixed!'
  nu=1.5e-5
endif
!
!
!
      height=Lxyz(2)/2
      h2=height**2
      B1=5.5
      kappa=2.5
      utau=Re_tau*nu/height
      lw=nu/utau
      !
      ! Add mean turbulent velocity profile to the possibly already
      ! existing turnulent velocity field. As there should be less turbulence
      ! close to the walls we scale the existing turbulence field with the
      ! velocity profile.
      !
      do j=m1,m2
        if (y(j)<xyz0(2)+height) then
          y_pluss=(y(j)-xyz0(2))/lw
          u_log=(kappa*log(y_pluss+tini)+B1)*utau
          u_lam=y_pluss*utau
          f(l1:l2,j,n1:n2,iux)=&
               (f(l1:l2,j,n1:n2,iux)/central_vel+1)*min(u_log,u_lam)
        else
          y_pluss=-(y(j)-xyz0(2)-Lxyz(2))/lw
          u_log=(kappa*log(y_pluss+tini)+B1)*utau
          u_lam=y_pluss*utau
          f(l1:l2,j,n1:n2,iux)=&
               (f(l1:l2,j,n1:n2,iux)/central_vel+1)*min(u_log,u_lam)
        endif
      enddo
      if (lfirst_proc_y) then
        f(l1:l2,m1,n1:n2,iux)=0.0
        f(l1:l2,m1+1,n1:n2,iux)=(y(m1+1)-xyz0(2))/lw*utau
      endif
      if (llast_proc_y) then
        f(l1:l2,m2,n1:n2,iux)=0.0
        f(l1:l2,m2-1,n1:n2,iux)=-(y(m2-1)-xyz0(2)-Lxyz(2))/lw*utau
      endif
      !
    endsubroutine log_law_flowx_wally
!***********************************************************************
    subroutine bc_poi_x(f,sgn,topbot,j,rel,val)
!
! Poiseulle inflow
!
!  03-jan-08/nils: Coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, optional :: val
      real :: umax,y2,height,h2
      integer :: sgn,i,j,jj
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) then
          ! Multiply by three halfs to get max velocity from mean velocity
          umax=val!*3/2
        else
          umax=0
        endif
!
        height=Lxyz(2)/2
        h2=height**2
!
        do jj=m1,m2
          y2=(y(jj)-xyz0(2)-height)**2
          f(l1,jj,n1:n2,j)=umax*(1-y2/h2)
        enddo
!
!
        if (relative) then
          do i=1,nghost; f(l1-i,:,:,j)=2*f(l1,:,:,j)+sgn*f(l1+i,:,:,j); enddo
        else
          do i=1,nghost; f(l1-i,:,:,j)=              sgn*f(l1+i,:,:,j); enddo
          f(l1,:,:,j)=(4.*f(l1+1,:,:,j)-f(l1+2,:,:,j))/3.
        endif
!
      case ('top')               ! top boundary
        if (present(val)) then
          umax=val
        else
          umax=0
        endif
!
!
        height=Lxyz(2)/2
        h2=height**2
!
        do jj=m1,m2
          y2=(y(jj)-xyz0(2)-height)**2
          f(l2,jj,n1:n2,j)=umax*(1-y2/h2)
        enddo
!
        if (relative) then
          do i=1,nghost; f(l2+i,:,:,j)=2*f(l2,:,:,j)+sgn*f(l2-i,:,:,j); enddo
        else
          do i=1,nghost; f(l2+i,:,:,j)=              sgn*f(l2-i,:,:,j); enddo
          f(l2,:,:,j)=(4.*f(l2-1,:,:,j)-f(l2-2,:,:,j))/3.
        endif
!
      case default
        print*, "bc_poi_x: ", topbot, " should be `top' or `bot'"
!
      endselect
!
    endsubroutine bc_poi_x
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
