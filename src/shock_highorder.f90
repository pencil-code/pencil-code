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
  use Cdata
  use Cparam
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'shock.h'
!
  integer :: ishock_max=1
  logical :: lgaussian_smooth=.false.
  logical :: lforce_periodic_shockviscosity=.false.
  real    :: div_threshold=0.
  logical :: lrewrite_shock_boundary=.false.
!
  namelist /shock_run_pars/ &
      ishock_max,lgaussian_smooth,lforce_periodic_shockviscosity,&
      div_threshold,lrewrite_shock_boundary
!
  integer :: idiag_shockm=0, idiag_shockmin=0, idiag_shockmax=0
  integer :: idiag_shockmx=0, idiag_shockmy=0, idiag_shockmz=0
!
  real, dimension (-3:3,-3:3,-3:3) :: smooth_factor
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
    subroutine initialize_shock(f,lstarting)
!
!  20-nov-02/tony: coded
!
       use Messages, only: fatal_error
!
       real, dimension (mx,my,mz,mfarray) :: f
       logical, intent(in) :: lstarting
!
       real, dimension (-3:3) :: weights
       integer :: i,j,k
!
!  Initialize shock profile to zero
!
      f(:,:,:,ishock)=0.0
!
!  Calculate the smoothing factors
!
      smooth_factor = 1.
!
      if (lgaussian_smooth) then
        weights = (/1.,9.,45.,70.,45.,9.,1./)
      else
        weights = (/1.,6.,15.,20.,15.,6.,1./)
      endif
!
      if (nxgrid > 1) then
        do i = -3,3
          smooth_factor(i,:,:) = smooth_factor(i,:,:)*weights(i)
        enddo
      else
        smooth_factor(-3:-1,:,:) = 0.
        smooth_factor(+1:+3,:,:) = 0.
      endif
!
      if (nygrid > 1) then
        do j = -3,3
          smooth_factor(:,j,:) = smooth_factor(:,j,:)*weights(j)
        enddo
      else
        smooth_factor(:,-3:-1,:) = 0.
        smooth_factor(:,+1:+3,:) = 0.
      endif
!
      if (nzgrid > 1) then
        do k = -3,3
          smooth_factor(:,:,k) = smooth_factor(:,:,k)*weights(k)
        enddo
      else
        smooth_factor(:,:,-3:-1) = 0.
        smooth_factor(:,:,+1:+3) = 0.
      endif
!
      smooth_factor = smooth_factor / sum(smooth_factor)
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
      if (.not. lstarting .and. .not. lforce_periodic_shockviscosity) then
        if (bcx(ishock)=='p' .and. .not. all(bcx(iux:iuz)=='p')) then
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
    subroutine read_shock_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_shock_init_pars
!***********************************************************************
    subroutine write_shock_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_shock_init_pars
!***********************************************************************
    subroutine read_shock_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=shock_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=shock_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_shock_run_pars
!***********************************************************************
    subroutine write_shock_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=shock_run_pars)
!
    endsubroutine write_shock_run_pars
!*******************************************************************
    subroutine rprint_shock(lreset,lwrite)
!
!  Writes ishock to index.pro file
!
!  24-nov-03/tony: adapted from rprint_ionization
!
      use Diagnostics
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
!  Write column where which shock variable is stored.
!
      if (present(lwrite)) then
        if (lwrite) write(3,*) 'ishock=',ishock
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
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Shock profile
!
        case ('shock')
          slices%yz= f(ix_loc,m1:m2 ,n1:n2  ,ishock)
          slices%xz= f(l1:l2 ,iy_loc,n1:n2  ,ishock)
          slices%xy= f(l1:l2 ,m1:m2 ,iz_loc ,ishock)
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ishock)
          if (lwrite_slice_xy3) &
              slices%xy3=f(l1:l2,m1:m2,iz3_loc,ishock)
          if (lwrite_slice_xy4) &
              slices%xy4=f(l1:l2,m1:m2,iz4_loc,ishock)
          slices%ready=.true.
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
!
      if (ldiagnos) then
        if (idiag_shockm/=0)   call sum_mn_name( p%shock,idiag_shockm)
        if (idiag_shockmin/=0) call max_mn_name(-p%shock,idiag_shockmin,lneg=.true.)
        if (idiag_shockmax/=0) call max_mn_name( p%shock,idiag_shockmax)
      endif
!
      if (l1davgfirst) then
        if (idiag_shockmx/=0)  call yzsum_mn_name_x(p%shock,idiag_shockmx)
        if (idiag_shockmy/=0)  call xzsum_mn_name_y(p%shock,idiag_shockmy)
        if (idiag_shockmz/=0)  call xysum_mn_name_z(p%shock,idiag_shockmz)
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
                         mpiallreduce_max
      use Sub, only: div
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
!
      integer, dimension (3) :: max_loc
      real, dimension (mx,my,mz) :: tmp
      real, dimension (nx) :: penc
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
! The following will calculate div u for any coordinate system.
        call div(f,iuu,penc)
        f(l1:l2,m,n,ishock) = max(0.0,-penc)
!
      enddo
!
!  Cut off small divergence if requested.
!
      if (div_threshold > 0.0) &
          where(f(:,:,:,ishock) < div_threshold) f(:,:,:,ishock) = 0.0
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
!  Scale given a fixed mesh Reynolds number or
!
      if (ldynamical_diffusion) then
        if (headtt) print *, 'Shock: fix mesh Reynolds number at ', re_mesh
        if (lfirst) then
          max_loc = (/ l1-1, m1-1, n1-1 /) + maxloc(tmp(l1:l2,m1:m2,n1:n2))
          a1 = tmp(max_loc(1),max_loc(2),max_loc(3))
          call mpiallreduce_max(a1, shock_max)
          if (shock_max > 0.) then
            if (shock_max == a1) then
              a1 = sqrt((f(max_loc(1),max_loc(2),max_loc(3),iux) * xprim(max_loc(1)))**2 &
                      + (f(max_loc(1),max_loc(2),max_loc(3),iuy) * yprim(max_loc(2)))**2 &
                      + (f(max_loc(1),max_loc(2),max_loc(3),iuz) * zprim(max_loc(3)))**2) / (pi * re_mesh * shock_max)
            else
              a1 = 0.
            endif
            call mpiallreduce_max(a1, a)
          endif
        endif
        f(:,:,:,ishock) = a * tmp
      else
!
!  Scale by dxmax**2
!
        if (.not.lrewrite_shock_boundary) then
          f(l1:l2,m1:m2,n1:n2,ishock) = tmp(l1:l2,m1:m2,n1:n2)*dxmax**2
        else
          f(:,:,:,ishock) = tmp*dxmax**2
        endif
      endif
!
    endsubroutine calc_shock_profile
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
endmodule Shock
