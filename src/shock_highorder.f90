! $Id$
!
!  This modules implements viscous heating and diffusion terms
!  here for shock viscosity
!    nu_total = nu + nu_shock*dx^2*smooth(max5(-(div u))))
!
!  NOTE: this works and has been tested for periodic boundaries.
!  With the current version, if your shock fronts reach a non-periodic
!  boundary, unexpected things may happen, so you should monitor the
!  behavior on the boundaries in this case.
!
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lshock = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
! COMMUNICATED AUXILIARIES 1
!
! PENCILS PROVIDED shock; gshock(3); shock_perp; gshock_perp(3)
!
!***************************************************************
module Shock
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include 'shock.h'
!
  integer :: ishock_max=1
  logical :: lshock_linear = .false.
  logical :: lgaussian_smooth=.false.
  logical :: lforce_periodic_shockviscosity=.false.
  logical :: lfix_Re_mesh=.false.
  real    :: div_threshold=0.0
  real    :: shock_linear = 0.01
  real    :: shock_div_pow = 1.
  logical :: lrewrite_shock_boundary=.false., lconvergence_only=.true.
!
  namelist /shock_run_pars/ &
      ishock_max, lgaussian_smooth, lforce_periodic_shockviscosity, &
      div_threshold, lrewrite_shock_boundary, lfix_Re_mesh, lshock_linear, shock_linear, &
      shock_div_pow, lconvergence_only
!
!  Diagnostic variables for print.in
! (needs to be consistent with reset list below)
!
  integer :: idiag_shockm=0        ! DIAG_DOC:
  integer :: idiag_shockmin=0      ! DIAG_DOC:
  integer :: idiag_shockmax=0      ! DIAG_DOC:
  integer :: idiag_gshockmax = 0   ! DIAG_DOC: $\max\left|\nabla\nu_{shock}\right|$
!
! xy averaged diagnostics given in xyaver.in written every it1d timestep
!
  integer :: idiag_shockmz=0       ! XYAVG_DOC:
!
! xz averaged diagnostics given in xzaver.in
!
  integer :: idiag_shockmy=0       ! XZAVG_DOC:
!
! yz averaged diagnostics given in yzaver.in
!
  integer :: idiag_shockmx=0       ! YZAVG_DOC:
!
  real, dimension (-3:3,-3:3,-3:3) :: smooth_factor
!
  interface shock_divu_perp
    module procedure shock_divu_perp_pencil
  endinterface
!
  contains
!***********************************************************************
    subroutine register_shock()
!
!  19-nov-02/tony: coded
!  24-jan-05/tony: modified from visc_shock.f90
!
      use FArrayManager
      use Messages, only: svn_id
!
      call farray_register_auxiliary('shock',ishock,communicated=.true.)
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Writing files for use with IDL
!
      if (naux+naux_com <  maux+maux_com) aux_var(aux_count)=',shock $'
      if (naux+naux_com  == maux+maux_com) aux_var(aux_count)=',shock'
      aux_count=aux_count+1
      if (lroot) write(15,*) 'shock = fltarr(mx,my,mz)*one'
!
    endsubroutine register_shock
!***********************************************************************
    subroutine initialize_shock(f)
!
!  20-nov-02/tony: coded
!
      use Messages, only: fatal_error
      use Sub, only: register_report_aux, smoothing_kernel
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (-3:3) :: weights
      integer :: i,j,k
!
!  Initialize shock profile to zero
!
      f(:,:,:,ishock)=0.0
      if (ldivu_perp) then
        call register_report_aux('shock_perp',ishock_perp,communicated=.true.)
        if (.not. lmagnetic) call fatal_error('initialize_shock', &
                    'ldivu_perp requires magnetic field module for B_perp')
        f(:,:,:,ishock_perp)=0.0
      endif
!
!  Calculate the smoothing factors
!
      call smoothing_kernel(smooth_factor,lgaussian_smooth)
!
!      if (lgaussian_smooth) then
!        weights = (/1.,9.,45.,70.,45.,9.,1./)
!      else
!        weights = (/1.,6.,15.,20.,15.,6.,1./)
!      endif
!!
!      if (nxgrid > 1) then
!        do i = -3,3
!          smooth_factor(i,:,:) = smooth_factor(i,:,:)*weights(i)
!        enddo
!      else
!        smooth_factor(-3:-1,:,:) = 0.
!        smooth_factor(+1:+3,:,:) = 0.
!      endif
!!
!      if (nygrid > 1) then
!        do j = -3,3
!          smooth_factor(:,j,:) = smooth_factor(:,j,:)*weights(j)
!        enddo
!      else
!        smooth_factor(:,-3:-1,:) = 0.
!        smooth_factor(:,+1:+3,:) = 0.
!      endif
!!
!      if (nzgrid > 1) then
!        do k = -3,3
!          smooth_factor(:,:,k) = smooth_factor(:,:,k)*weights(k)
!        enddo
!      else
!        smooth_factor(:,:,-3:-1) = 0.
!        smooth_factor(:,:,+1:+3) = 0.
!      endif
!!
!      smooth_factor = smooth_factor / sum(smooth_factor)
!
!  Check that smooth order is within bounds
!
      if (lroot) then
        if (ishock_max < 1.or.ishock_max > 3) then
          call fatal_error('initialize_shock', &
                           'ishock_max needs to be between 1 and 3.')
        endif
      endif
!
!  Die if periodic boundary condition for shock viscosity, but not for
!  velocity field. It can lead to subtle numerical errors near non-
!  periodic boundaries if the shock viscosity is assumed periodic.
!
      if (lrun .and. .not. lforce_periodic_shockviscosity) then
        if ((bcx(ishock) == 'p') .and. .not. all(bcx(iux:iuz) == 'p')) then
          if (lroot) then
            print*, 'initialize_shock: shock viscosity has bcx=''p'', but the velocity field is not'
            print*, '                  periodic! (you must set a proper boundary condition for the'
            print*, '                  shock viscosity)'
            print*, 'initialize_shock: bcx=', bcx
            print*, 'initialize_shock: to suppress this error,'
            print*, '                  set lforce_periodic_shockviscosity=T in &shock_run_pars'
          endif
          call fatal_error('initialize_shock','')
        endif
        if (bcy(ishock)=='p' .and. .not. all(bcy(iux:iuz)=='p')) then
          if (lroot) then
            print*, 'initialize_shock: shock viscosity has bcy=''p'', but the velocity field is not'
            print*, '                  periodic! (you must set a proper boundary condition for the'
            print*, '                  shock viscosity)'
            print*, 'initialize_shock: bcy=', bcy
            print*, 'initialize_shock: to suppress this error,'
            print*, '                  set lforce_periodic_shockviscosity=T in &shock_run_pars'
          endif
          call fatal_error('initialize_shock','')
        endif
        if (bcz(ishock)=='p' .and. .not. all(bcz(iux:iuz)=='p')) then
          if (lroot) then
            print*, 'initialize_shock: shock viscosity has bcz=''p'', but the velocity field is not'
            print*, '                  periodic! (you must set a proper boundary condition for the'
            print*, '                  shock viscosity)'
            print*, 'initialize_shock: bcz=', bcz
            print*, 'initialize_shock: to suppress this error,'
            print*, '                  set lforce_periodic_shockviscosity=T in &shock_run_pars'
          endif
          call fatal_error('initialize_shock','')
        endif
      endif
!
    endsubroutine initialize_shock
!***********************************************************************
    subroutine read_shock_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=shock_run_pars, IOSTAT=iostat)
!
    endsubroutine read_shock_run_pars
!***********************************************************************
    subroutine write_shock_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=shock_run_pars)
!
    endsubroutine write_shock_run_pars
!***********************************************************************
    subroutine rprint_shock(lreset,lwrite)
!
!  Writes ishock to index.pro file
!
!  24-nov-03/tony: adapted from rprint_ionization
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamex, inamey, inamez
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_shockm=0; idiag_shockmin=0; idiag_shockmax=0
        idiag_gshockmax = 0
        idiag_shockmx=0; idiag_shockmy=0; idiag_shockmz=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_shock: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),&
            'shockm',idiag_shockm)
        call parse_name(iname,cname(iname),cform(iname),&
            'shockmin',idiag_shockmin)
        call parse_name(iname,cname(iname),cform(iname),&
            'shockmax',idiag_shockmax)
        call parse_name(iname, cname(iname), cform(iname), 'gshockmax', idiag_gshockmax)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'shockmx',idiag_shockmx)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'shockmy',idiag_shockmy)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'shockmz',idiag_shockmz)
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then 
        where(cnamev=='shock') cformv='DEFINED'
      endif
!
!  Write column where which shock variable is stored.
!
      if (present(lwrite)) then
        if (lwrite) then
          call farray_index_append('ishock',ishock)
          call farray_index_append('ishock_perp',ishock_perp)
        endif
      endif
!
    endsubroutine rprint_shock
!***********************************************************************
    subroutine get_slices_shock(f,slices)
!
!  Write slices for animation of shock variable.
!
!  26-jul-06/tony: coded
!
      use Slices_methods, only: assign_slices_scal

      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Shock profile
!
        case ('shock'); call assign_slices_scal(slices,f,ishock)
!
      endselect
!
    endsubroutine get_slices_shock
!***********************************************************************
    subroutine pencil_criteria_shock()
!
!  All pencils that the Viscosity module depends on are specified here.
!
!  20-11-04/anders: coded
!
      if (idiag_shockm/=0 .or. idiag_shockmin/=0 .or. idiag_shockmax/=0 .or. &
          idiag_shockmx/=0 .or. idiag_shockmy/=0 .or. idiag_shockmz/=0) then
        lpenc_diagnos(i_shock)=.true.
      endif
!
      if (idiag_gshockmax /= 0) lpenc_diagnos(i_gshock) = .true.
!
    endsubroutine pencil_criteria_shock
!***********************************************************************
    subroutine pencil_interdep_shock(lpencil_in)
!
!  Interdependency among pencils from the Viscosity module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension (npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_shock
!***********************************************************************
    subroutine calc_pencils_shock(f,p)
!
!  Calculate Viscosity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      use Diagnostics
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
! shock
      if (lpencil(i_shock)) p%shock=f(l1:l2,m,n,ishock)
! gshock
      if (lpencil(i_gshock)) call grad(f,ishock,p%gshock)
      if (ldivu_perp) then
! shock_perp
        if (lpencil(i_shock_perp)) p%shock_perp=f(l1:l2,m,n,ishock_perp)
! gshock_perp
        if (lpencil(i_gshock_perp)) call grad(f,ishock_perp,p%gshock_perp)
      endif
!
      if (ldiagnos) then
        if (idiag_shockm/=0)   call sum_mn_name( p%shock,idiag_shockm)
        if (idiag_shockmin/=0) call max_mn_name(-p%shock,idiag_shockmin,lneg=.true.)
        if (idiag_shockmax/=0) call max_mn_name( p%shock,idiag_shockmax)
        if (idiag_gshockmax /= 0) call max_mn_name(sqrt(sum(p%gshock**2, dim=2)), idiag_gshockmax)
      endif
!
      if (l1davgfirst) then
        call yzsum_mn_name_x(p%shock,idiag_shockmx)
        call xzsum_mn_name_y(p%shock,idiag_shockmy)
        call xysum_mn_name_z(p%shock,idiag_shockmz)
      endif
!
    endsubroutine calc_pencils_shock
!***********************************************************************
    subroutine calc_shock_profile(f)
!
!  Calculate divu based shock profile to be used in viscosity and
!  diffusion type terms.
!
!  23-nov-02/tony: coded
!  17-dec-08/ccyang: add divergence threshold
!
      use Boundcond, only: boundconds_x, boundconds_y, boundconds_z
      use Mpicomm, only: initiate_isendrcv_bdry, finalize_isendrcv_bdry, &
                         mpiallreduce_max, MPI_COMM_WORLD
      use Magnetic, only: bb_unitvec_shock
      use Sub, only: div
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
!
      real, dimension (mx,my,mz) :: tmp
      real, dimension (nx) :: penc, penc1, penc_perp
      real, dimension (mx,3) :: bb_hat
      integer :: imn
      integer :: i,j,k
      integer :: ni,nj,nk
      real :: shock_max, a=0., a1
!
!  Initialize shock to impossibly large number to force code crash in case of
!  inconsistencies in boundary conditions.
!
      f(:,:,:,ishock)=impossible
!
!  Compute divergence of velocity field.
!
      call boundconds_x(f,iux,iux)
!
      call initiate_isendrcv_bdry(f,iux,iuz)
!
      if (lshock_linear) a1 = shock_linear / dxmax
!
      do imn=1,ny*nz
!
        n = nn(imn)
        m = mm(imn)
!
        if (necessary(imn)) then
          call finalize_isendrcv_bdry(f,iux,iuz)
          call boundconds_y(f,iuy,iuy)
          call boundconds_z(f,iuz,iuz)
        endif
!
! The following will calculate div u for any coordinate system.
!
        call div(f,iuu,penc)
        if (lconvergence_only) then
          f(l1:l2,m,n,ishock) = max(0.0,-penc)
        else
          f(l1:l2,m,n,ishock) = abs(penc)
        endif
        if (shock_div_pow /= 1.) f(l1:l2,m,n,ishock)=f(l1:l2,m,n,ishock)**shock_div_pow
!
!  Add the linear term if requested
!
        if (lshock_linear) f(l1:l2,m,n,ishock) = f(l1:l2,m,n,ishock) + a1 * wave_speed(f)
!
      enddo
!
!  Cut off small divergence if requested.
!
      if (div_threshold > 0.0) &
          where(abs(f(:,:,:,ishock)) < div_threshold) f(:,:,:,ishock) = 0.0
!
!  Take maximum over a number of grid cells
!
      ni = merge(ishock_max,0,nxgrid > 1)
      nj = merge(ishock_max,0,nygrid > 1)
      nk = merge(ishock_max,0,nzgrid > 1)
!
!  Because of a bug in the shearing boundary conditions we must first manually
!  set the y boundary conditions on the shock profile.
!
      if (lshear) then
        call boundconds_y(f,ishock,ishock)
        call initiate_isendrcv_bdry(f,ishock,ishock)
        call finalize_isendrcv_bdry(f,ishock,ishock)
      endif
!
      call boundconds_x(f,ishock,ishock)
!
      call initiate_isendrcv_bdry(f,ishock,ishock)
!
      tmp = 0.0
!
      do imn=1,ny*nz
!
        n = nn(imn)
        m = mm(imn)
!
        if (necessary(imn)) then
          call finalize_isendrcv_bdry(f,ishock,ishock)
          call boundconds_y(f,ishock,ishock)
          call boundconds_z(f,ishock,ishock)
        endif
!
        penc = 0.0
!
        do k=-nk,nk
        do j=-nj,nj
        do i=-ni,ni
          penc = max(penc,f(l1+i:l2+i,m+j,n+k,ishock))
        enddo
        enddo
        enddo
!
        tmp(l1:l2,m,n) = penc
!
      enddo
!
      f(:,:,:,ishock) = tmp
!
!  Smooth with a Gaussian profile
!
      ni = merge(3,0,nxgrid > 1)
      nj = merge(3,0,nygrid > 1)
      nk = merge(3,0,nzgrid > 1)
!
!  Because of a bug in the shearing boundary conditions we must first manually
!  set the y boundary conditions on the shock profile.
!
      if (lshear) then
        call boundconds_y(f,ishock,ishock)
        call initiate_isendrcv_bdry(f,ishock,ishock)
        call finalize_isendrcv_bdry(f,ishock,ishock)
      endif
!
      call boundconds_x(f,ishock,ishock)
      call initiate_isendrcv_bdry(f,ishock,ishock)
!
      tmp = 0.0
!
      do imn=1,ny*nz
!
        n = nn(imn)
        m = mm(imn)
!
        if (necessary(imn)) then
          call finalize_isendrcv_bdry(f,ishock,ishock)
          call boundconds_y(f,ishock,ishock)
          call boundconds_z(f,ishock,ishock)
        endif
!
        penc = 0.0
!
        do k=-nk,nk
        do j=-nj,nj
        do i=-ni,ni
          penc = penc + smooth_factor(i,j,k)*f(l1+i:l2+i,m+j,n+k,ishock)
        enddo
        enddo
        enddo
!
        tmp(l1:l2,m,n) = penc
!
      enddo
!
      fix_Re: if (lfix_Re_mesh) then
!
!  Scale given a fixed mesh Reynolds number.
!
        if (headtt) print *, 'Shock: fix mesh Reynolds number'
        call mpiallreduce_max(maxval(tmp(l1:l2,m1:m2,n1:n2)), shock_max, comm=MPI_COMM_WORLD)
        shock: if (shock_max > 0.) then
          a1 = 0.
          scan_z: do n = n1, n2
            scan_y: do m = m1, m2
              penc = 0.
              penc1 = 0.
              xdir: if (nxgrid > 1) then
                penc = penc + abs(f(l1:l2,m,n,iux) * dx_1(l1:l2))
                penc1 = penc1 + dx_1(l1:l2)**2
              endif xdir
              ydir: if (nygrid > 1) then
                penc = penc + abs(f(l1:l2,m,n,iuy) * dy_1(m))
                penc1 = penc1 + dy_1(m)**2
              endif ydir
              zdir: if (nzgrid > 1) then
                penc = penc + abs(f(l1:l2,m,n,iuz) * dz_1(n))
                penc1 = penc1 + dz_1(n)**2
              endif zdir
              a1 = max(a1, maxval(penc / penc1))
            enddo scan_y
          enddo scan_z
          call mpiallreduce_max(a1, a, comm=MPI_COMM_WORLD)
          a = a / (re_mesh * pi * shock_max)
        else shock
          a = dxmin**2
        endif shock
        f(:,:,:,ishock) = a * tmp
      else fix_Re
!
!  Scale by dxmin**2.
!
!  The shearing boundary conditions have a bug that can cause nonsensical
!  numbers in the corners (bug #61). We can overwrite this bug by defining the
!  shock viscosity in the ghost zones too. THIS WOULD NOT BE NECESSARY
!  IF THE BUG IN THE SHEARING BOUNDARY CONDITIONS WOULD BE FIXED.
!
        if (.not.lrewrite_shock_boundary) then
          f(l1:l2,m1:m2,n1:n2,ishock) = tmp(l1:l2,m1:m2,n1:n2)*dxmin**2
        else
          f(:,:,:,ishock) = tmp*dxmin**2
        endif
      endif fix_Re
!
      if (ldivu_perp) then
        call boundconds_x(f,ishock_perp,ishock_perp)
        call initiate_isendrcv_bdry(f,ishock_perp,ishock_perp)
!
        do imn=1,ny*nz
!
          n = nn(imn)
          m = mm(imn)
!
          if (necessary(imn)) then
            call finalize_isendrcv_bdry(f,ishock_perp,ishock_perp)
            call boundconds_y(f,ishock_perp,ishock_perp)
            call boundconds_z(f,ishock_perp,ishock_perp)
          endif
!
          call div(f,iuu,penc)
          call bb_unitvec_shock(f,bb_hat)
          call shock_divu_perp(f,bb_hat,penc,penc_perp)
          f(l1:l2,m,n,ishock_perp)=max(0.,-penc_perp)
!
        enddo
        tmp = 0.0
        do imn=1,ny*nz
!
          n = nn(imn)
          m = mm(imn)
!
          if (necessary(imn)) then
            call finalize_isendrcv_bdry(f,ishock_perp,ishock_perp)
            call boundconds_y(f,ishock_perp,ishock_perp)
            call boundconds_z(f,ishock_perp,ishock_perp)
          endif
!
          penc = 0.0
!
          do k=-nk,nk
          do j=-nj,nj
          do i=-ni,ni
            penc = penc + smooth_factor(i,j,k)*f(l1+i:l2+i,m+j,n+k,ishock_perp)
          enddo
          enddo
          enddo
!
          tmp(l1:l2,m,n) = penc
!
        enddo
        if (.not.lrewrite_shock_boundary) then
          f(l1:l2,m1:m2,n1:n2,ishock_perp) = tmp(l1:l2,m1:m2,n1:n2) * dxmin**2
        else
          f(:,:,:,ishock_perp) = tmp * dxmin**2
        endif
!
      endif
!
    endsubroutine calc_shock_profile
!***********************************************************************
    subroutine shock_divu_perp_pencil(f,bb_hat,divu,divu_perp)
!
!  Calculate `perpendicular divergence' of u.
!  nabla_perp.uu = nabla.uu - (1/b2)*bb.(bb.nabla)*uu
!
!  16-aug-06/tobi: coded
!  07-jun-18/fred: revised to include higher order gradient u
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,3), intent (in) :: bb_hat
      real, dimension (nx), intent (in) :: divu
      real, dimension (nx), intent (out) :: divu_perp
!
      real, dimension (nx) :: fac
!
      divu_perp = divu
!
      if (nxgrid/=1) then
         fac=bb_hat(l1:l2,1)/(60*dx)
         divu_perp(:) = divu_perp(:)                          &
           - fac*(bb_hat(l1:l2,1)*( f(l1+3:l2+3,m  ,n  ,iux)   &
                                   + 45.0*(f(l1+1:l2+1,m  ,n  ,iux)   &
                                         - f(l1-1:l2-1,m  ,n  ,iux))  &
                                   -  9.0*(f(l1+2:l2+2,m  ,n  ,iux)   &
                                         - f(l1-2:l2-2,m  ,n  ,iux))  &
                                         - f(l1-3:l2-3,m  ,n  ,iux) ) &
                       + bb_hat(l1:l2,2)*( f(l1+3:l2+3,m  ,n  ,iuy)   &
                                   + 45.0*(f(l1+1:l2+1,m  ,n  ,iuy)   &
                                         - f(l1-1:l2-1,m  ,n  ,iuy))  &
                                   -  9.0*(f(l1+2:l2+2,m  ,n  ,iuy)   &
                                         - f(l1-2:l2-2,m  ,n  ,iuy))  &
                                         - f(l1-3:l2-3,m  ,n  ,iuy) ) &
                       + bb_hat(l1:l2,3)*( f(l1+3:l2+3,m  ,n  ,iuz)   &
                                   + 45.0*(f(l1+1:l2+1,m  ,n  ,iuz)   &
                                         - f(l1-1:l2-1,m  ,n  ,iuz))  &
                                   -  9.0*(f(l1+2:l2+2,m  ,n  ,iuz)   &
                                         - f(l1-2:l2-2,m  ,n  ,iuz))  &
                                         - f(l1-3:l2-3,m  ,n  ,iuz) ) )
      endif
!
      if (nygrid/=1) then
         fac=bb_hat(l1:l2,2)/(60*dy)
         divu_perp(:) = divu_perp(:)                          &
           - fac*(bb_hat(l1:l2,1)*( f(  l1:l2  ,m+3,n  ,iux)   &
                                   + 45.0*(f(  l1:l2  ,m+1,n  ,iux)   &
                                         - f(  l1:l2  ,m-1,n  ,iux))  &
                                   -  9.0*(f(  l1:l2  ,m+2,n  ,iux)   &
                                         - f(  l1:l2  ,m-2,n  ,iux))  &
                                         - f(  l1:l2  ,m-3,n  ,iux) ) &
                       + bb_hat(l1:l2,2)*( f(  l1:l2  ,m+3,n  ,iuy)   &
                                   + 45.0*(f(  l1:l2  ,m+1,n  ,iuy)   &
                                         - f(  l1:l2  ,m-1,n  ,iuy))  &
                                   -  9.0*(f(  l1:l2  ,m+2,n  ,iuy)   &
                                         - f(  l1:l2  ,m-2,n  ,iuy))  &
                                         - f(  l1:l2  ,m-3,n  ,iuy) ) &
                       + bb_hat(l1:l2,3)*( f(  l1:l2  ,m+3,n  ,iuz)   &
                                   + 45.0*(f(  l1:l2  ,m+1,n  ,iuz)   &
                                         - f(  l1:l2  ,m-1,n  ,iuz))  &
                                   -  9.0*(f(  l1:l2  ,m+2,n  ,iuz)   &
                                         - f(  l1:l2  ,m-2,n  ,iuz))  &
                                         - f(  l1:l2  ,m-3,n  ,iuz) ) )
      endif
!
      if (nzgrid/=1) then
         fac=bb_hat(l1:l2,3)/(60*dz)
         divu_perp(:) = divu_perp(:)                          &
           - fac*(bb_hat(l1:l2,1)*( f(  l1:l2  ,m  ,n+3,iux)   &
                                   + 45.0*(f(  l1:l2  ,m,  n+1,iux)   &
                                         - f(  l1:l2  ,m,  n-1,iux))  &
                                   -  9.0*(f(  l1:l2  ,m,  n+2,iux)   &
                                         - f(  l1:l2  ,m,  n-2,iux))  &
                                         - f(  l1:l2  ,m  ,n-3,iux) ) &
                       + bb_hat(l1:l2,2)*( f(  l1:l2  ,m  ,n+3,iuy)   &
                                   + 45.0*(f(  l1:l2  ,m,  n+1,iuy)   &
                                         - f(  l1:l2  ,m,  n-1,iuy))  &
                                   -  9.0*(f(  l1:l2  ,m,  n+2,iuy)   &
                                         - f(  l1:l2  ,m,  n-2,iuy))  &
                                         - f(  l1:l2  ,m  ,n-3,iuy) ) &
                       + bb_hat(l1:l2,3)*( f(  l1:l2  ,m  ,n+3,iuz)   &
                                   + 45.0*(f(  l1:l2  ,m,  n+1,iuz)   &
                                         - f(  l1:l2  ,m,  n-1,iuz))  &
                                   -  9.0*(f(  l1:l2  ,m,  n+2,iuz)   &
                                         - f(  l1:l2  ,m,  n-2,iuz))  &
                                         - f(  l1:l2  ,m  ,n-3,iuz) ) )
      endif
!
    endsubroutine shock_divu_perp_pencil
!***********************************************************************
    subroutine calc_shock_profile_simple(f)
!
!  Calculate divu based shock profile to be used in viscosity and
!  diffusion type terms.
!
!  12-apr-05/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_shock_profile_simple
!***********************************************************************
    function wave_speed(f) result(speed)
!
!  Calculate the wave speeds along one pencil.
!
!  23-nov-13/ccyang: coded.
!
      use EquationOfState, only: eoscalc, rho0
      use Magnetic, only: get_bext
      use Sub, only: gij, curl_mn
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx) :: speed
!
      real, dimension(nx,3,3) :: aij
      real, dimension(nx,3) :: bb
      real, dimension(nx) :: b2
      real, dimension(3) :: B_ext
!
!  Get the sound speed.
!
      call eoscalc(f, nx, cs2=speed)
!
!  Add the Alfven speed.
!
      alfven: if (lmagnetic) then
!
        bfield: if (lbfield) then
          bb = f(l1:l2,m,n,ibx:ibz)
        else bfield
          call gij(f, iaa, aij, 1)
          call curl_mn(aij, bb, f(:,m,n,iax:iaz))
        endif bfield
        call get_bext(B_ext)
        b2 = sum((bb + spread(B_ext,1,nx))**2, dim=2)
!
        density: if (ldensity) then
          if (ldensity_nolog) then
            speed = speed + mu01 * b2 / f(l1:l2,m,n,irho)
          else
            speed = speed + mu01 * b2 / exp(f(l1:l2,m,n,ilnrho))
          endif
        else density
          speed = speed + mu01 / rho0 * b2
        endif density
!
      endif alfven
!
      speed = sqrt(speed)
!
    endfunction wave_speed
!***********************************************************************
endmodule Shock
