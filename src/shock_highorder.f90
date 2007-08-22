! $Id: shock_highorder.f90,v 1.6 2007-08-22 18:59:04 dhruba Exp $

!  This modules implements viscous heating and diffusion terms
!  here for shock viscosity
!    nu_total = nu + nu_shock*dx^2*smooth(max5(-(div u))))
!
!  NOTE: this works and has been tested for periodic boundaries.
!  With the current version, if your shock fronts reach a non-periodic
!  boundary, unexpected things may happen, so you should monitor the
!  behavior on the boundaries in this case.
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
! PENCILS PROVIDED shock,gshock,shock_perp,gshock_perp
!
!***************************************************************

module Shock

  use Cparam

  implicit none

  include 'shock.h'

  integer :: ishock_max = 1
  logical :: lgaussian_smooth = .false.

  ! run parameters
  namelist /shock_run_pars/ ishock_max,lgaussian_smooth

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_shockmax=0

  real, dimension (-3:3,-3:3,-3:3) :: smooth_factor

  contains

!***********************************************************************
    subroutine register_shock()
!
!  19-nov-02/tony: coded
!  24-jan-05/tony: modified from visc_shock.f90
!
      use Cdata
      use Mpicomm, only: stop_it
      use Messages, only: cvs_id
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_shock called twice')
      first = .false.
!
      ishock = mvar + naux_com + 1
      naux = naux + 1
      naux_com = naux_com + 1
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_shock: shock viscosity nvar = ', nvar
        print*, 'ishock = ', ishock
      endif
!
!  Put variable name in array
!
      varname(ishock) = 'shock'
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: shock_highorder.f90,v 1.6 2007-08-22 18:59:04 dhruba Exp $")
!
! Check we aren't registering too many auxiliary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux= ', maux
        call stop_it('register_shock: naux > maux')
      endif
      if (naux_com > maux_com) then
        if (lroot) write(0,*) 'naux_com = ', naux_com, ', maux_com = ', maux_com
        call stop_it('register_shock: naux_com > maux_com')
      endif
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
       use Cdata, only: ishock,lroot
       use Messages, only: fatal_error

       real, dimension (mx,my,mz,mfarray) :: f
       logical, intent(in) :: lstarting

       real, dimension (-3:3) :: weights
       integer :: i,j,k
!
!  Initialize shock profile to zero
!
      f(:,:,:,ishock)=0.0

!
!  Calculate factors for polynomial smoothing
!
      smooth_factor = 1.

      if (lgaussian_smooth) then
        weights = (/1.,9.,45.,70.,45.,9.,1./)
      else
        weights = (/1.,6.,15.,20.,15.,6.,1./)
      endif

      if (nxgrid > 1) then
        do i = -3,3
          smooth_factor(i,:,:) = smooth_factor(i,:,:)*weights(i)
        enddo
      else
        smooth_factor(-3:-1,:,:) = 0.
        smooth_factor(+1:+3,:,:) = 0.
      endif

      if (nygrid > 1) then
        do j = -3,3
          smooth_factor(:,j,:) = smooth_factor(:,j,:)*weights(j)
        enddo
      else
        smooth_factor(:,-3:-1,:) = 0.
        smooth_factor(:,+1:+3,:) = 0.
      endif

      if (nzgrid > 1) then
        do k = -3,3
          smooth_factor(:,:,k) = smooth_factor(:,:,k)*weights(k)
        enddo
      else
        smooth_factor(:,:,-3:-1) = 0.
        smooth_factor(:,:,+1:+3) = 0.
      endif

      smooth_factor = smooth_factor/sum(smooth_factor)

!
!  Check that smooth order is within bounds
!
      if (lroot) then
        if (ishock_max < 1.or.ishock_max > 3) then
          call fatal_error('initialize_shock', &
                           'ishock_max needs to be between 1 and 3.')
        endif
      endif

    endsubroutine initialize_shock
!***********************************************************************
    subroutine read_shock_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit
!
    endsubroutine read_shock_init_pars
!***********************************************************************
    subroutine write_shock_init_pars(unit)
      integer, intent(in) :: unit
!
      if (NO_WARN) print*,unit
!
    endsubroutine write_shock_init_pars
!***********************************************************************
    subroutine read_shock_run_pars(unit,iostat)
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
      integer, intent(in) :: unit

      write(unit,NML=shock_run_pars)

    endsubroutine write_shock_run_pars
!*******************************************************************
    subroutine rprint_shock(lreset,lwrite)
!
!  Writes ishock to index.pro file
!
!  24-nov-03/tony: adapted from rprint_ionization
!
      use Cdata, only: lroot,ip
      use Cdata, only: nname,cname,cform
      use Cdata, only: ishock
      use Sub, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
      integer :: iname
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_shockmax=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_shock: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),&
            'shockmax',idiag_shockmax)
      enddo
!
!  write column where which shock variable is stored
!
      if (present(lwrite)) then
        if (lwrite) then
          write(3,*) 'i_shockmax=',idiag_shockmax
          write(3,*) 'ishock=',ishock
        endif
      endif
!
      if (NO_WARN) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_shock
!***********************************************************************
    subroutine get_slices_shock(f,slices)
!
!  Write slices for animation of shock variable.
!
!  26-jul-06/tony: coded
!
      use Cdata, only: ishock

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
          slices%yz= f(slices%ix,m1:m2    ,n1:n2     ,ishock)
          slices%xz= f(l1:l2    ,slices%iy,n1:n2     ,ishock)
          slices%xy= f(l1:l2    ,m1:m2    ,slices%iz ,ishock)
          slices%xy2=f(l1:l2    ,m1:m2    ,slices%iz2,ishock)
          slices%ready = .true.

      endselect

    endsubroutine get_slices_shock
!!***********************************************************************
    subroutine pencil_criteria_shock()
!
!  All pencils that the Viscosity module depends on are specified here.
!
!  20-11-04/anders: coded
!
!   dummy
!
      if (idiag_shockmax/=0) then
          lpenc_diagnos(i_shock)=.true.
      endif

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
      if (NO_WARN) print*, lpencil_in !(keep compiler quiet)
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
      use Cdata, only: m,n,ishock
      use Cdata, only: ldiagnos
      use Sub, only: grad,max_mn_name

      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p

      intent(in) :: f
      intent(inout) :: p

      ! shock
      if (lpencil(i_shock)) p%shock=f(l1:l2,m,n,ishock)

      ! gshock
      if (lpencil(i_gshock)) call grad(f,ishock,p%gshock)

      if (ldiagnos) then
        if (idiag_shockmax/=0) call max_mn_name(p%shock,idiag_shockmax)
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
!
      use Cdata, only: m,n,mm,nn,necessary
      use Cdata, only: iux,iuy,iuz,ishock
      use Cdata, only: dxmax
! For spherical polar coordinate system 
      use Cdata, only: r1_mn,cotth,lspherical_coords
! For cylindrical coordinate system
      use Cdata, only: rcyl_mn1,lcylindrical_coords
      use Boundcond, only: boundconds_x,boundconds_y,boundconds_z
      use Mpicomm, only: initiate_isendrcv_bdry,finalize_isendrcv_bdry
      use Deriv, only: der

      real, dimension (mx,my,mz,mfarray), intent (inout) :: f

      real, dimension (mx,my,mz) :: tmp
      real, dimension (nx) :: penc
      integer :: imn
      integer :: i,j,k
      integer :: ni,nj,nk

!
!  Compute divergence
!
      call boundconds_x(f,iux,iux)

      call initiate_isendrcv_bdry(f,iuy,iuz)

      do imn=1,ny*nz

        n = nn(imn)
        m = mm(imn)

        if (necessary(imn)) then
          call finalize_isendrcv_bdry(f,iuy,iuz)
          call boundconds_y(f,iuy,iuy)
          call boundconds_z(f,iuz,iuz)
        endif

        f(l1:l2,m,n,ishock) = 0.

        if (nxgrid > 1) then
          call der(f,iux,penc,1)
          f(l1:l2,m,n,ishock) = f(l1:l2,m,n,ishock) + max(0.,-penc)
        endif

        if (nygrid > 1) then
          call der(f,iuy,penc,2)
         
          f(l1:l2,m,n,ishock) = f(l1:l2,m,n,ishock) + max(0.,-penc)
        endif

        if (nzgrid > 1) then
          call der(f,iuz,penc,3)
          f(l1:l2,m,n,ishock) = f(l1:l2,m,n,ishock) + max(0.,-penc)
        endif
        if(lspherical_coords) then
          penc= r1_mn*f(l1:l2,m,n,iux) + r1_mn*cotth(m)*f(l1:l2,m,n,iuy)
          f(l1:l2,m,n,ishock) = f(l1:l2,m,n,ishock) + max(0.,-penc)
        endif
        if(lcylindrical_coords) then
          penc= r1_mn*f(l1:l2,m,n,iux) 
          f(l1:l2,m,n,ishock) = f(l1:l2,m,n,ishock) + max(0.,-penc)
        endif
      enddo

!
!  Take maximum over a number of grid cells
!
      ni = merge(ishock_max,0,nxgrid > 1)
      nj = merge(ishock_max,0,nygrid > 1)
      nk = merge(ishock_max,0,nzgrid > 1)

      call boundconds_x(f,ishock,ishock)
      call initiate_isendrcv_bdry(f,ishock,ishock)

      tmp = 0.

      do imn=1,ny*nz

        n = nn(imn)
        m = mm(imn)

        if (necessary(imn)) then
          call finalize_isendrcv_bdry(f,ishock,ishock)
          call boundconds_y(f,ishock,ishock)
          call boundconds_z(f,ishock,ishock)
        endif

        penc = 0.

        do k=-nk,nk
        do j=-nj,nj
        do i=-ni,ni
          penc = max(penc,f(l1+i:l2+i,m+j,n+k,ishock))
        enddo
        enddo
        enddo

        tmp(l1:l2,m,n) = penc

      enddo

      f(:,:,:,ishock) = tmp

!
!  Smooth with a Gaussian profile
!
      ni = merge(3,0,nxgrid > 1)
      nj = merge(3,0,nygrid > 1)
      nk = merge(3,0,nzgrid > 1)

      call boundconds_x(f,ishock,ishock)
      call initiate_isendrcv_bdry(f,ishock,ishock)

      tmp = 0.

      do imn=1,ny*nz

        n = nn(imn)
        m = mm(imn)

        if (necessary(imn)) then
          call finalize_isendrcv_bdry(f,ishock,ishock)
          call boundconds_y(f,ishock,ishock)
          call boundconds_z(f,ishock,ishock)
        endif

        penc = 0.

        do k=-nk,nk
        do j=-nj,nj
        do i=-ni,ni
          penc = penc + smooth_factor(i,j,k)*f(l1+i:l2+i,m+j,n+k,ishock)
        enddo
        enddo
        enddo

        tmp(l1:l2,m,n) = penc

      enddo

!
!  Scale by dxmax**2
!
      f(:,:,:,ishock) = tmp*dxmax**2

    endsubroutine calc_shock_profile
!!***********************************************************************
    subroutine calc_shock_profile_simple(f)
!
!  Calculate divu based shock profile to be used in viscosity and
!  diffusion type terms.
!
!  12-apr-05/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f

    endsubroutine calc_shock_profile_simple
!***********************************************************************
endmodule Shock
