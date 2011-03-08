! $Id$
!
!  Module for boundary conditions. Extracted from (no)mpicomm, since
!  all non-periodic (external) boundary conditions require the same
!  code for serial and parallel runs.
!
module Boundcond
!
  use Cdata
  use Cparam
  use Messages
  use Mpicomm
!
  implicit none
!
  private
!
  public :: update_ghosts
  public :: boundconds, boundconds_x, boundconds_y, boundconds_z
  public :: bc_per_x, bc_per_y, bc_per_z
!
  contains
!***********************************************************************
    subroutine update_ghosts(a)
!
!  Update all ghost zones of a.
!
!  21-sep-02/wolf: extracted from wsnaps
!
      real, dimension (mx,my,mz,mfarray) :: a
!
      call boundconds_x(a)
      call initiate_isendrcv_bdry(a)
      call finalize_isendrcv_bdry(a)
      call boundconds_y(a)
      call boundconds_z(a)
!
    endsubroutine update_ghosts
!***********************************************************************
    subroutine boundconds(f,ivar1_opt,ivar2_opt)
!
!  Apply boundary conditions in all three directions.
!  Note that we _must_ call boundconds_{x,y,z} in this order, or edges and
!  corners will not be OK.
!
!  10-oct-02/wolf: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
      call boundconds_x(f,ivar1,ivar2)
      call boundconds_y(f,ivar1,ivar2)
      call boundconds_z(f,ivar1,ivar2)
!
    endsubroutine boundconds
!***********************************************************************
    subroutine boundconds_x(f,ivar1_opt,ivar2_opt)
!
!  Boundary conditions in x, except for periodic part handled by communication.
!  Remark: boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-'corners').
!
!   8-jul-02/axel: split up into different routines for x,y and z directions
!  11-nov-02/wolf: unified bot/top, now handled by loop
!  15-dec-06/wolf: Replaced "if (bcx1(1)=='she') then" by "any" command
!
      use EquationOfState
      use Shear
      use Special, only: special_boundconds
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      real, dimension (mcom) :: fbcx12
      real, dimension (mcom) :: fbcx2_12
      integer :: ivar1, ivar2, j, k, one
      logical :: ip_ok
      character (len=bclen), dimension(mcom) :: bc12
      character (len=3) :: topbot
      type (boundary_condition) :: bc
!
      if (ldebug) print*, 'boundconds_x: ENTER: boundconds_x'
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
      select case (nxgrid)
!
      case (1)
        if (ldebug) print*, 'boundconds_x: no x-boundary'
!
!  Boundary conditions in x.
!
      case default
!
!  Use the following construct to keep compiler from complaining if
!  we have no variables (and boundconds) at all (samples/no-modules):
!
        one = min(1,mcom)
        if (any(bcx1(1:one)=='she')) then
          call boundcond_shear(f,ivar1,ivar2)
        else
          do k=1,2
            ! loop over 'bot','top'
            if (k==1) then
              topbot='bot'; bc12=bcx1; fbcx12=fbcx1; fbcx2_12=fbcx1_2; ip_ok=lfirst_proc_x
            else
              topbot='top'; bc12=bcx2; fbcx12=fbcx2; fbcx2_12=fbcx2_2; ip_ok=llast_proc_x
            endif
!
            do j=ivar1,ivar2
!
! Natalia: the next line is for the dustdensity case.
! If ndustspec is large, it is stupid to set bc for all dust species
! in start.in. But if one does not set them, they becomes 'p' by default
! Since this problem is crutial  only for aerosol + chemistry
! the following condition is used. But this place should be modifyed somehow
! Any ideas?
!
!            if ((bc12(j)=='p') .and. lchemistry .and. ldustdensity) bc12(j)=''
!
              if (ldebug) write(*,'(A,I1,A,I2,A,A)') ' bcx',k,'(',j,')=',bc12(j)
              if (ip_ok) then
                select case (bc12(j))
                case ('0')
                  ! BCX_DOC: zero value in ghost zones, free value on boundary
                  call bc_zero_x(f,topbot,j)
                case ('p')
                  ! BCX_DOC: periodic
                  call bc_per_x(f,topbot,j)
                case ('s')
                  ! BCX_DOC: symmetry, $f_{N+i}=f_{N-i}$;
                  ! BCX_DOC: implies $f'(x_N)=f'''(x_0)=0$
                  call bc_sym_x(f,+1,topbot,j)
                case ('ss')
                  ! BCX_DOC: symmetry, plus function value given
                  call bc_symset_x(f,+1,topbot,j)
                case ('s0d')
                  ! BCX_DOC: symmetry, function value such that df/dx=0
                  call bc_symset0der_x(f,topbot,j)
                case ('a')
                  ! BCX_DOC: antisymmetry, $f_{N+i}=-f_{N-i}$;
                  ! BCX_DOC: implies $f(x_N)=f''(x_0)=0$
                  call bc_sym_x(f,-1,topbot,j)
                case ('a2')
                  ! BCX_DOC: antisymmetry relative to boundary value,
                  ! BCX_DOC: $f_{N+i}=2 f_{N}-f_{N-i}$;
                  ! BCX_DOC: implies $f''(x_0)=0$
                  call bc_sym_x(f,-1,topbot,j,REL=.true.)
                case ('cpc')
                  ! BCX_DOC: cylindrical perfect conductor
                  ! BCX_DOC: implies $f''+f'/R=0$
                  call bc_cpc_x(f,topbot,j)
                case ('cpp')
                  ! BCX_DOC: cylindrical perfect conductor
                  ! BCX_DOC: implies $f''+f'/R=0$
                  call bc_cpp_x(f,topbot,j)
                case ('cpz')
                  ! BCX_DOC: cylindrical perfect conductor
                  ! BCX_DOC: implies $f''+f'/R=0$
                  call bc_cpz_x(f,topbot,j)
                case ('v')
                  ! BCX_DOC: vanishing third derivative
                  call bc_van_x(f,topbot,j)
                case ('cop')
                  ! BCX_DOC: copy value of last physical point to all ghost cells
                  call bc_copy_x(f,topbot,j)
                case ('1s')
                  ! BCX_DOC: onesided
                  call bc_onesided_x(f,topbot,j)
                case ('1so')
                  ! BCX_DOC: onesided
                  call bc_onesided_x_old(f,topbot,j)
                case ('cT')
                  ! BCX_DOC: constant temperature (implemented as
                  ! BCX_DOC: condition for entropy $s$ or temperature $T$)
                  call bc_ss_temp_x(f,topbot)
                case ('c1')
                  ! BCX_DOC: constant temperature (or maybe rather constant
                  ! BCX_DOC: conductive flux??)
                  if (j==iss)   call bc_ss_flux_x(f,topbot)
                  if (j==ilnTT) call bc_lnTT_flux_x(f,topbot)
                case ('Fgs')
                  ! BCX_DOC: Fconv = - chi_t*rho*T*grad(s)
                  if (j==iss) call bc_ss_flux_turb_x(f,topbot)
                case ('sT')
                  ! BCX_DOC: symmetric temperature, $T_{N-i}=T_{N+i}$;
                  ! BCX_DOC: implies $T'(x_N)=T'''(x_0)=0$
                  if (j==iss) call bc_ss_stemp_x(f,topbot)
                case ('asT')
                  ! BCX_DOC: select entropy for uniform ghost temperature
                  ! BCX_DOC: matching fluctuating boundary value,
                  ! BCX_DOC: $T_{N-i}=T_{N}=$;
                  ! BCX_DOC: implies $T'(x_N)=T'(x_0)=0$
                  if (j==iss) call bc_ss_a2stemp_x(f,topbot)
                case ('db')
                  ! BCX_DOC:
                  call bc_db_x(f,topbot,j)
                case ('f')
                  ! BCX_DOC: ``freeze'' value, i.e. maintain initial
                  !  value at boundary
                  call bc_freeze_var_x(topbot,j)
                  call bc_sym_x(f,-1,topbot,j,REL=.true.)
                  ! antisymm wrt boundary
                case ('fg')
                  ! BCX_DOC: ``freeze'' value, i.e. maintain initial
                  !  value at boundary, also mantaining the
                  !  ghost zones at the initial coded value, i.e.,
                  !  keep the gradient frozen as well
                  call bc_freeze_var_x(topbot,j)
                case ('1')
                  ! BCX_DOC: $f=1$ (for debugging)
                  call bc_one_x(f,topbot,j)
                case ('set')
                  ! BCX_DOC: set boundary value to \var{fbcx12}
                  call bc_sym_x(f,-1,topbot,j,REL=.true.,val=fbcx12)
                case ('der')
                  ! BCX_DOC: set derivative on boundary to \var{fbcx12}
                  call bc_set_der_x(f,topbot,j,fbcx12(j))
                case ('slo')
                  ! BCX_DOC: set slope at the boundary = \var{fbcx12}
                  call bc_slope_x(f,fbcx12,topbot,j)
                case ('dr0')
                  ! BCX_DOC: set boundary value [really??]
                  call bc_dr0_x(f,fbcx12,topbot,j)
                case ('ovr')
                  ! BCX_DOC: overshoot boundary condition
                  ! BCX_DOC:  ie $(d/dx-1/\mathrm{dist}) f = 0.$
                  call bc_overshoot_x(f,fbcx12,topbot,j)
                case ('ant')
                  ! BCX_DOC: stops and prompts for adding documentation
                  call bc_antis_x(f,fbcx12,topbot,j)
                case ('e1')
                  ! BCX_DOC: extrapolation [describe]
                  call bcx_extrap_2_1(f,topbot,j)
                case ('e2')
                  ! BCX_DOC: extrapolation [describe]
                  call bcx_extrap_2_2(f,topbot,j)
                case ('e3')
                  ! BCX_DOC: extrapolation in log [maintain a power law]
                  call bcx_extrap_2_3(f,topbot,j)
                case ('hat')
                  ! BCX_DOC: top hat jet profile in spherical coordinate.
                  !Defined only for the bottom boundary
                  call bc_set_jethat_x(f,j,topbot,fbcx12,fbcx2_12)
                case ('spd')
                  ! BCX_DOC:  sets $d(rA_{\alpha})/dr = \mathtt{fbcx12(j)}$
                  call bc_set_spder_x(f,topbot,j,fbcx12(j))
                case ('sfr')
                  ! BCX_DOC: stress-free boundary condition
                  ! BCX_DOC: for spherical coordinate system.
                  call bc_set_sfree_x(f,topbot,j)
                case ('nfr')
                  ! BCX_DOC: Normal-field bc for spherical coordinate system.
                  ! BCX_DOC: Some people call this the ``(angry) hedgehog bc''.
                  call bc_set_nfr_x(f,topbot,j)
                case ('sa2')
                  ! BCX_DOC: $(d/dr)(r B_{\phi}) = 0$ imposes
                  ! BCX_DOC: boundary condition on 2nd derivative of
                  ! BCX_DOC: $r A_{\phi}$. Same applies to $\theta$ component.
                  call bc_set_sa2_x(f,topbot,j)
                case ('pfc')
                  ! BCX_DOC: perfect-conductor in spherical
                  ! BCX_DOC: coordinate: $d/dr( A_r) + 2/r = 0$ .
                  call bc_set_pfc_x(f,topbot,j)
                case ('fix')
                  ! BCX_DOC: set boundary value [really??]
                  call bc_fix_x(f,topbot,j,fbcx12(j))
                case ('fil')
                  ! BCX_DOC: set boundary value from a file
                  call bc_file_x(f,topbot,j)
                case ('cfb')
                  ! BCZ_DOC: radial centrifugal balance
                  if (lcylindrical_coords) then
                    call bc_lnrho_cfb_r_iso(f,topbot)
                  else
                    print*,'not implemented for other than cylindrical'
                    stop
                  endif
                case ('g')
                  ! BCX_DOC: set to given value(s) or function
                  call bc_force_x(f, -1, topbot, j)
                case ('nil')
                  ! BCX_DOC: do nothing; assume that everything is set
                case ('ioc')
                  ! BCX_DOC: inlet/outlet on western/eastern hemisphere
                  ! BCX_DOC: in cylindrical coordinates
                  call bc_inlet_outlet_cyl(f,topbot,j,fbcx12)
                case ('')
                  ! BCX_DOC: do nothing; assume that everything is set
                case default
                  bc%bcname=bc12(j)
                  bc%ivar=j
                  bc%location=(((k-1)*2)-1)   ! -1/1 for x bot/top
                  bc%value1=fbcx12(j)
                  bc%value2=fbcx12(j)
                  bc%done=.false.
!
                  call special_boundconds(f,bc)
!
                  if (.not.bc%done) then
                    write(unit=errormsg,fmt='(A,A4,A,I3)') &
                         "No such boundary condition bcx1/2 = ", &
                         bc12(j), " for j=", j
                    call stop_it_if_any(.true.,"boundconds_x: "//trim(errormsg))
                  endif
                endselect
              endif
            enddo
          enddo
        endif
      endselect
!
      ! Catch any 'stop_it_if_any' calls from single MPI ranks that may
      ! have occured inside the above select statement. This final call
      ! for all MPI ranks is necessary to prevent dead-lock situations.
      call stop_it_if_any(.false.,"")
!
    endsubroutine boundconds_x
!***********************************************************************
    subroutine boundconds_y(f,ivar1_opt,ivar2_opt)
!
!  Boundary conditions in y, except for periodic part handled by communication.
!  Remark: boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-'corners').
!
!   8-jul-02/axel: split up into different routines for x,y and z directions
!  11-nov-02/wolf: unified bot/top, now handled by loop
!
      use Special, only: special_boundconds
      use EquationOfState
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      real, dimension (mcom) :: fbcy12
      integer :: ivar1, ivar2, j, k
      logical :: ip_ok
      character (len=bclen), dimension(mcom) :: bc12
      character (len=3) :: topbot
      type (boundary_condition) :: bc
!
      if (ldebug) print*,'boundconds_y: ENTER: boundconds_y'
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
      select case (nygrid)
!
      case (1)
        if (ldebug) print*,'boundconds_y: no y-boundary'
!
!  Boundary conditions in y
!
      case default
        do k=1,2                ! loop over 'bot','top'
          if (k==1) then
            topbot='bot'; bc12=bcy1; fbcy12=fbcy1; ip_ok=lfirst_proc_y
          else
            topbot='top'; bc12=bcy2; fbcy12=fbcy2; ip_ok=llast_proc_y
          endif
!
          do j=ivar1,ivar2
!
! Natalia: the next line is for the dustdensity case.
! If ndustspec is large, it is stupid to set bc for all dust species
! in start.in. But if one does not set them, they becomes 'p' by default
! Since this problem is crutial  only for aerosol + chemistry
! the following condition is used. But this place should be modifyed somehow
! Any ideas?
!
!            if ((bc12(j)=='p') .and. lchemistry .and. ldustdensity) bc12(j)=''
!
            if (ldebug) write(*,'(A,I1,A,I2,A,A)') ' bcy',k,'(',j,')=',bc12(j)
            if (ip_ok) then
              select case (bc12(j))
              case ('p')
                ! BCY_DOC: periodic
                call bc_per_y(f,topbot,j)
              case ('pp')
                ! BCY_DOC: periodic across the pole
                call bc_pper_y(f,topbot,j)
              case ('ap')
                ! BCY_DOC: anti-periodic across the pole
                call bc_aper_y(f,topbot,j)
              case ('s')
                ! BCY_DOC: symmetry symmetry, $f_{N+i}=f_{N-i}$;
                  ! BCX_DOC: implies $f'(y_N)=f'''(y_0)=0$
                call bc_sym_y(f,+1,topbot,j)
              case ('ss')
                ! BCY_DOC: symmetry, plus function value given
                call bc_symset_y(f,+1,topbot,j)
              case ('s0d')
                ! BCY_DOC: symmetry, function value such that df/dy=0
                call bc_symset0der_y(f,topbot,j)
              case ('a')
                ! BCY_DOC: antisymmetry
                call bc_sym_y(f,-1,topbot,j)
              case ('a2')
                ! BCY_DOC: antisymmetry relative to boundary value
                call bc_sym_y(f,-1,topbot,j,REL=.true.)
              case ('v')
                ! BCY_DOC: vanishing third derivative
                call bc_van_y(f,topbot,j)
              case ('1s')
                ! BCY_DOC: onesided
                call bc_onesided_y(f,topbot,j)
              case ('cT')
                ! BCY_DOC: constant temp.
                if (j==iss) call bc_ss_temp_y(f,topbot)
              case ('sT')
                ! BCY_DOC: symmetric temp.
                if (j==iss) call bc_ss_stemp_y(f,topbot)
              case ('asT')
                ! BCY_DOC: select entropy for uniform ghost temperature
                ! BCY_DOC: matching fluctuating boundary value,
                ! BCY_DOC: $T_{N-i}=T_{N}=$;
                ! BCY_DOC: implies $T'(x_N)=T'(x_0)=0$
                if (j==iss) call bc_ss_a2stemp_y(f,topbot)
              case ('f')
                ! BCY_DOC: freeze value
                ! tell other modules not to change boundary value
                call bc_freeze_var_y(topbot,j)
                call bc_sym_y(f,-1,topbot,j,REL=.true.) ! antisymm wrt boundary
              case ('s+f')
                ! BCY_DOC: freeze value
                ! tell other modules not to change boundary value
                call bc_freeze_var_y(topbot,j)
                call bc_sym_y(f,+1,topbot,j) ! symm wrt boundary
              case ('fg')
                ! BCY_DOC: ``freeze'' value, i.e. maintain initial
                !  value at boundary, also mantaining the
                !  ghost zones at the initial coded value, i.e.,
                !  keep the gradient frozen as well
                call bc_freeze_var_y(topbot,j)
              case ('fBs')
                ! BCZ_DOC: frozen-in B-field (s)
                call bc_frozen_in_bb(topbot,j)
                call bc_sym_y(f,+1,topbot,j) ! symmetry
              case ('1')
                ! BCY_DOC: f=1 (for debugging)
                call bc_one_y(f,topbot,j)
              case ('set')
                ! BCY_DOC: set boundary value
                call bc_sym_y(f,-1,topbot,j,REL=.true.,val=fbcy12)
              case ('e1')
                ! BCY_DOC: extrapolation
                call bcy_extrap_2_1(f,topbot,j)
              case ('e2')
                ! BCY_DOC: extrapolation
                call bcy_extrap_2_2(f,topbot,j)
              case ('e3')
                ! BCY_DOC: extrapolation in log [maintain a power law]
                call bcy_extrap_2_3(f,topbot,j)
              case ('der')
                ! BCY_DOC: set derivative on the boundary
                call bc_set_der_y(f,topbot,j,fbcy12(j))
              case ('cop')
                ! BCY_DOC: outflow: copy value of last physical point to
                ! BCY_DOC: all ghost cells
                call bc_copy_y(f,topbot,j)
              case ('c+k')
                ! BCY_DOC: no-inflow: copy value of last physical point
                ! BCY_DOC: to all ghost cells, but suppressing any inflow
                call bc_copy_y_noinflow(f,topbot,j)
              case ('sfr')
                ! BCY_DOC: stress-free boundary condition for spherical
                ! BCY_DOC: coordinate system.
                call bc_set_sfree_y(f,topbot,j)
              case ('nfr')
                ! BCY_DOC: Normal-field bc for spherical coordinate system.
                ! BCY_DOC: Some people call this the ``(angry) hedgehog bc''.
                call bc_set_nfr_y(f,topbot,j)
              case ('pfc')
                ! BCY_DOC: perfect conducting boundary condition
                ! BCY_DOC: along $\theta$ boundary
                call bc_set_pfc_y(f,topbot,j)
              case ('twi')
                call bc_twist_xz(f,topbot,j)
              case ('nil','')
                ! BCY_DOC: do nothing; assume that everything is set
              case default
                bc%bcname=bc12(j)
                bc%ivar=j
                bc%value1=fbcy12(j)
                bc%value2=fbcy12(j)
                bc%location=(((k-1)*4)-2)   ! -2/2 for y bot/top
                bc%done=.false.
!
                if (lspecial) call special_boundconds(f,bc)
!
                if (.not.bc%done) then
                  write(unit=errormsg,fmt='(A,A4,A,I3)') "No such boundary condition bcy1/2 = ", &
                       bc12(j), " for j=", j
                  call stop_it_if_any(.true.,"boundconds_y: "//trim(errormsg))
                endif
              endselect
            endif
          enddo
        enddo
      endselect
!
      ! Catch any 'stop_it_if_any' calls from single MPI ranks that may
      ! have occured inside the above select statement. This final call
      ! for all MPI ranks is necessary to prevent dead-lock situations.
      call stop_it_if_any(.false.,"")
!
    endsubroutine boundconds_y
!***********************************************************************
    subroutine boundconds_z(f,ivar1_opt,ivar2_opt)
!
!  Boundary conditions in z, except for periodic part handled by communication.
!  Remark: boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-'corners').
!
!   8-jul-02/axel: split up into different routines for x,y and z directions
!  11-nov-02/wolf: unified bot/top, now handled by loop
!
      use Special, only: special_boundconds
      use EquationOfState
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
      real, dimension (mcom) :: fbcz12, fbcz12_1, fbcz12_2, fbcz_zero=0.
      integer :: ivar1, ivar2, j, k
      logical :: ip_ok
      character (len=bclen), dimension(mcom) :: bc12
      character (len=3) :: topbot
      type (boundary_condition) :: bc
!
      if (ldebug) print*,'boundconds_z: ENTER: boundconds_z'
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
      select case (nzgrid)
!
      case (1)
        if (ldebug) print*,'boundconds_z: no z-boundary'
!
!  Boundary conditions in z
!
      case default
        do k=1,2                ! loop over 'bot','top'
          if (k==1) then
            topbot='bot'
            bc12=bcz1
            fbcz12=fbcz1
            fbcz12_1=fbcz1_1
            fbcz12_2=fbcz1_2
            ip_ok=lfirst_proc_z
          else
            topbot='top'
            bc12=bcz2
            fbcz12=fbcz2
            fbcz12_1=fbcz2_1
            fbcz12_2=fbcz2_2
            ip_ok=llast_proc_z
          endif
!
          do j=ivar1,ivar2
            if (ldebug) write(*,'(A,I1,A,I2,A,A)') ' bcz',k,'(',j,')=',bc12(j)
            if (ip_ok) then
              select case (bc12(j))
              case ('0')
                ! BCZ_DOC: zero value in ghost zones, free value on boundary
                call bc_zero_z(f,topbot,j)
              case ('p')
                ! BCZ_DOC: periodic
                call bc_per_z(f,topbot,j)
              case ('s')
                ! BCZ_DOC: symmetry
                call bc_sym_z(f,+1,topbot,j)
              case ('s0d')
                ! BCZ_DOC: symmetry, function value such that df/dz=0
                call bc_symset0der_z(f,topbot,j)
              case ('0ds')
                ! BCZ_DOC: symmetry, function value such that df/dz=0
                call bc_symset0der_z_v2(f,topbot,j)
              case ('a')
                ! BCZ_DOC: antisymmetry
                call bc_sym_z(f,-1,topbot,j)
              case ('a2')
                ! BCZ_DOC: antisymmetry relative to boundary value
                call bc_sym_z(f,-1,topbot,j,REL=.true.)
              case ('a0d')
                ! BCZ_DOC: antisymmetry with zero derivative
                call bc_sym_z(f,+1,topbot,j,VAL=fbcz_zero)
              case ('v')
                ! BCZ_DOC: vanishing third derivative
                call bc_van_z(f,topbot,j)
              case ('v3')
                ! BCZ_DOC: vanishing third derivative
                call bc_van3rd_z(f,topbot,j)
              case ('1s')
                ! BCZ_DOC: one-sided
                call bc_onesided_z(f,topbot,j)
              case ('fg')
                ! BCZ_DOC: ``freeze'' value, i.e. maintain initial
                !  value at boundary, also mantaining the
                !  ghost zones at the initial coded value, i.e.,
                !  keep the gradient frozen as well
                call bc_freeze_var_z(topbot,j)
              case ('c1')
                ! BCZ_DOC: complex
                if (j==iss) call bc_ss_flux(f,topbot)
                if (j==iaa) call bc_aa_pot(f,topbot)
                if (j==ilnTT) call bc_lnTT_flux_z(f,topbot)
              case ('Fgs')
                ! BCZ_DOC: Fconv = - chi_t*rho*T*grad(s)
                if (j==iss) call bc_ss_flux_turb(f,topbot)
              case ('c3')
                ! BCZ_DOC: constant flux at the bottom with a variable hcond
                if (j==ilnTT) call bc_ADI_flux_z(f,topbot)
              case ('pfe')
                ! BCZ_DOC: potential field extrapolation
                if (j==iaa) call bc_aa_pot_field_extrapol(f,topbot)
              case ('p1D')
                ! BCZ_DOC: potential field extrapolation in 1D
                if (j==iay) call bc_aa_pot_1D(f,topbot)
              case ('pot')
                ! BCZ_DOC: potential magnetic field
                if (j==iaa) call bc_aa_pot2(f,topbot)
              case ('pwd')
                ! BCZ_DOC: a variant of 'pot' for nprocx=1
                if (j==iaa) call bc_aa_pot3(f,topbot)
              case ('d2z')
                ! BCZ_DOC:
                call bc_del2zero(f,topbot,j)
              case ('hds')
                ! BCZ_DOC: hydrostatic equilibrium with a high-frequency filter
                call bc_lnrho_hdss_z_iso(f,topbot)
              case ('cT')
                ! BCZ_DOC: constant temp.
                if (j==ilnrho) call bc_lnrho_temp_z(f,topbot)
                call bc_ss_temp_z(f,topbot)
              case ('cT2')
                ! BCZ_DOC: constant temp. (keep lnrho)
                if (j==iss)   call bc_ss_temp2_z(f,topbot)
              case ('hs')
                ! BCZ_DOC: hydrostatic equilibrium
                if (.not.lgrav) call fatal_error('boundconds_z', &
                  'hs boundary condition requires gravity')
                if (j==ilnrho) call bc_lnrho_hds_z_iso(f,topbot)
                if (j==irho_b) call bc_lnrho_hds_z_iso(f,topbot)
                if (j==ipp)    call bc_pp_hds_z_iso(f,topbot)
              case ('hse')
                ! BCZ_DOC: hydrostatic extrapolation
                ! BCZ_DOC: rho or lnrho is extrapolated linearily and the
                ! BCZ_DOC: temperature is calculated in hydrostatic equilibrium.
                if (.not. lgrav) &
                    call fatal_error ('boundconds_z', "'hse' requires gravity")
                if (.not. leos) call fatal_error ('boundconds_z', &
                    "'hse' requires an eos module")
                if ((ilnrho == 0) .or. (ilnTT == 0)) &
                    call fatal_error ('boundconds_z', "'hse' requires lnrho and lnTT")
                if (j /= ilnTT) call fatal_error ('boundconds_z', &
                    "'hse' works only in lnTT")
                call bcz_hydrostatic_temp(f,topbot)
              case ('cp')
                ! BCZ_DOC: constant pressure
                ! BCZ_DOC:
                if (j==ilnrho) call bc_lnrho_pressure_z(f,topbot)
              case ('sT')
                ! BCZ_DOC: symmetric temp.
                ! BCZ_DOC:
                if (j==iss) call bc_ss_stemp_z(f,topbot)
              case ('ism')
                ! BCZ_DOC: special for interstellar runs
                if (j==iuz) call bc_steady_z(f,topbot,iuz)
                if (j==irho.or.j==ilnrho) call bc_onesided_z(f,topbot,j)
                if (j==iss) call bc_ctz(f,topbot,iss)
              case ('asT')
                ! BCZ_DOC: select entropy for uniform ghost temperature
                ! BCZ_DOC: matching fluctuating boundary value,
                ! BCZ_DOC: $T_{N-i}=T_{N}=$;
                ! BCZ_DOC: implies $T'(x_N)=T'(x_0)=0$
                if (j==iss) call bc_ss_a2stemp_z(f,topbot)
              case ('c2')
                ! BCZ_DOC: complex
                ! BCZ_DOC:
                if (j==iss) call bc_ss_temp_old(f,topbot)
              case ('db')
                ! BCZ_DOC: complex
                ! BCZ_DOC:
                call bc_db_z(f,topbot,j)
              case ('ce')
                ! BCZ_DOC: complex
                ! BCZ_DOC:
                if (j==iss) call bc_ss_energy(f,topbot)
              case ('e1')
                ! BCZ_DOC: extrapolation
                call bc_extrap_2_1(f,topbot,j)
              case ('e2')
                ! BCZ_DOC: extrapolation
                call bc_extrap_2_2(f,topbot,j)
              case ('ex')
                ! BCZ_DOC: simple linear extrapolation in first order
                call bcz_extrapol(f,topbot,j)
              case ('exf')
                ! BCZ_DOC: simple linear extrapolation in first order
                !  with a fixed value in the first ghost cell
                call bcz_extrapol_fixed(f,topbot,j)
              case ('exd')
                ! BCZ_DOC: simple linear extrapolation in first order
                !  with an included damping to zero (useful for velocities)
                call bcz_extrapol_damped(f,topbot,j)
              case ('b1')
                ! BCZ_DOC: extrapolation with zero value (improved 'a')
                call bc_extrap0_2_0(f,topbot,j)
              case ('b2')
                ! BCZ_DOC: extrapolation with zero value (improved 'a')
                call bc_extrap0_2_1(f,topbot,j)
              case ('b3')
                ! BCZ_DOC: extrapolation with zero value (improved 'a')
                call bc_extrap0_2_2(f,topbot,j)
              case ('f')
                ! BCZ_DOC: freeze value
                ! tell other modules not to change boundary value
                call bc_freeze_var_z(topbot,j)
                call bc_sym_z(f,-1,topbot,j,REL=.true.) ! antisymm wrt boundary
              case ('fBs')
                ! BCZ_DOC: frozen-in B-field (s)
                call bc_frozen_in_bb(topbot,j)
                call bc_sym_z(f,+1,topbot,j) ! symmetry
              case ('fB')
                ! BCZ_DOC: frozen-in B-field (a2)
                call bc_frozen_in_bb(topbot,j)
                call bc_sym_z(f,-1,topbot,j,REL=.true.) ! antisymm wrt boundary
              case ('g')
                ! BCZ_DOC: set to given value(s) or function
                 call bc_force_z(f,-1,topbot,j)
              case ('gs')
                ! BCZ_DOC:
                 call bc_force_z(f,+1,topbot,j)
              case ('1')
                ! BCZ_DOC: f=1 (for debugging)
                call bc_one_z(f,topbot,j)
              case ('StS')
                ! BCZ_DOC: solar surface boundary conditions
                if (j==ilnrho) call bc_stellar_surface(f,topbot)
              case ('set')
                ! BCZ_DOC: set boundary value
                call bc_sym_z(f,-1,topbot,j,REL=.true.,val=fbcz12)
              case ('der')
                ! BCZ_DOC: set derivative on the boundary
                call bc_set_der_z(f,topbot,j,fbcz12(j))
              case ('div')
                ! BCZ_DOC: set the divergence of $\uv$ to a given value
                ! BCZ_DOC: use bc = 'div' for iuz
                call bc_set_div_z(f,topbot,j,fbcz12(j))
              case ('ovr')
                ! BCZ_DOC: set boundary value
                call bc_overshoot_z(f,fbcz12,topbot,j)
              case ('inf')
                ! BCZ_DOC: allow inflow, but no outflow
                call bc_inflow_z(f,topbot,j)
              case ('ouf')
                ! BCZ_DOC: allow outflow, but no inflow
                call bc_outflow_z(f,topbot,j)
              case ('in')
                ! BCZ_DOC: allow inflow, but no outflow
                ! BCZ_DOC: forces ghost cells and boundary to not point outwards
                call bc_inflow_z(f,topbot,j,.true.)
              case ('out')
                ! BCZ_DOC: allow outflow, but no inflow
                ! BCZ_DOC: forces ghost cells and boundary to not point inwards
                call bc_outflow_z(f,topbot,j,.true.)
              case ('in0')
                ! BCZ_DOC: allow inflow, but no outflow
                ! BCZ_DOC: forces ghost cells and boundary to not point outwards
                ! BCZ_DOC: relaxes to vanishing 1st derivative at boundary
                call bc_inflow_zero_deriv_z(f,topbot,j)
              case ('ou0')
                ! BCZ_DOC: allow outflow, but no inflow
                ! BCZ_DOC: forces ghost cells and boundary to not point inwards
                ! BCZ_DOC: relaxes to vanishing 1st derivative at boundary
                call bc_outflow_zero_deriv_z(f,topbot,j)
              case ('ubs')
                ! BCZ_DOC: symmetric outflow,
                ! but match boundary inflow (experimental)
                call bc_steady_z(f,topbot,j)
              case ('win')
                ! BCZ_DOC: forces massflux given as
                ! BCZ_DOC: $\Sigma \rho_i ( u_i + u_0)=\textrm{fbcz1/2}(\rho)$
                if (j==ilnrho) then
                   call bc_wind_z(f,topbot,fbcz12(j))     !
                   call bc_sym_z(f,+1,topbot,j)           !  's'
                   call bc_sym_z(f,+1,topbot,iuz)         !  's'
                endif
              case ('cop')
                ! BCZ_DOC: copy value of last physical point to all ghost cells
                call bc_copy_z(f,topbot,j)
              case ('nil')
                ! BCZ_DOC: do nothing; assume that everything is set
              case default
                bc%bcname=bc12(j)
                bc%ivar=j
                bc%location=(((k-1)*6)-3)   ! -3/3 for z bot/top
                bc%value1=fbcz12_1(j)
                bc%value2=fbcz12_2(j)
                bc%done=.false.
!
                if (lspecial) call special_boundconds(f,bc)
!
                if (.not.bc%done) then
                  write(unit=errormsg,fmt='(A,A4,A,I3)') "No such boundary condition bcz1/2 = ", &
                       bc12(j), " for j=", j
                  call stop_it_if_any(.true.,"boundconds_z: "//trim(errormsg))
                endif
              endselect
            endif
          enddo
        enddo
      endselect
!
! Catch any 'stop_it_if_any' calls from single MPI ranks that may
! have occured inside the above select statement. This final call
! for all MPI ranks is necessary to prevent dead-lock situations.
!
      call stop_it_if_any(.false.,"")
!
    endsubroutine boundconds_z
!***********************************************************************
    subroutine bc_per_x(f,topbot,j)
!
!  Periodic boundary condition
!
!  11-nov-02/wolf: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
      character (len=3) :: topbot
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (nprocx==1) f(1:l1-1,:,:,j) = f(l2i:l2,:,:,j)
!
      case ('top')               ! top boundary
        if (nprocx==1) f(l2+1:mx,:,:,j) = f(l1:l1i,:,:,j)
!
      case default
        print*, "bc_per_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_per_x
!***********************************************************************
    subroutine bc_per_y(f,topbot,j)
!
!  Periodic boundary condition
!
!  11-nov-02/wolf: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
      character (len=3) :: topbot
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (nprocy==1) f(:,1:m1-1,:,j) = f(:,m2i:m2,:,j)
!
      case ('top')               ! top boundary
        if (nprocy==1) f(:,m2+1:my,:,j) = f(:,m1:m1i,:,j)
!
      case default
        print*, "bc_per_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_per_y
!***********************************************************************
    subroutine bc_pper_y(f,topbot,j)
!
!  Periodic boundary condition across the pole
!
!  15-june-10/dhruba: aped
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j,nhalf
      character (len=3) :: topbot
!
      nhalf=(n1+n2)/2
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (nprocz==1) then
           f(:,1,n1:nhalf,j) = f(:,m1+2,nhalf+1:n2,j)
           f(:,2,n1:nhalf,j) = f(:,m1+1,nhalf+1:n2,j)
           f(:,3,n1:nhalf,j) = f(:,m1,nhalf+1:n2,j)
!
           f(:,1,nhalf+1:n2,j) = f(:,m1+2,n1:nhalf,j)
           f(:,2,nhalf+1:n2,j) = f(:,m1+1,n1:nhalf,j)
           f(:,3,nhalf+1:n2,j) = f(:,m1,n1:nhalf,j)
        endif
      case ('top')               ! top boundary
        if (nprocz==1) then
           f(:,m2+1,n1:nhalf,j) = f(:,m2,nhalf+1:n2,j)
           f(:,m2+2,n1:nhalf,j) = f(:,m2-1,nhalf+1:n2,j)
           f(:,m2+3,n1:nhalf,j) = f(:,m2-2,nhalf+1:n2,j)
!
           f(:,m2+1,nhalf+1:n2,j) = f(:,m2,n1:nhalf,j)
           f(:,m2+2,nhalf+1:n2,j) = f(:,m2-1,n1:nhalf,j)
           f(:,m2+3,nhalf+1:n2,j) = f(:,m2-2,n1:nhalf,j)
        endif
      case default
        print*, "bc_pper_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_pper_y
!***********************************************************************
    subroutine bc_aper_y(f,topbot,j)
!
!  Anti-periodic boundary condition across the pole
!
!  15-june-10/dhruba: aped
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j,nhalf
      character (len=3) :: topbot
!
      nhalf=(n1+n2)/2
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (nprocz==1) then
           f(:,1:m1-1,n1:nhalf,j) = -f(:,m1:m1i,nhalf+1:n2,j)
           f(:,1:m1-1,nhalf+1:n2,j) = -f(:,m1:m1i,n1:nhalf,j)
        endif
      case ('top')               ! top boundary
        if (nprocz==1) then
           f(:,m2+1:my,n1:nhalf,j) = -f(:,m2i:m2,nhalf+1:n2,j)
           f(:,m2+1:my,nhalf+1:n2,j) = -f(:,m2i:m2,n1:nhalf,j)
        endif
      case default
        print*, "bc_aper_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_aper_y
!***********************************************************************
    subroutine bc_per_z(f,topbot,j)
!
!  Periodic boundary condition
!
!  11-nov-02/wolf: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
      character (len=3) :: topbot
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (nprocz==1) f(:,:,1:n1-1,j) = f(:,:,n2i:n2,j)
!
      case ('top')               ! top boundary
        if (nprocz==1) f(:,:,n2+1:mz,j) = f(:,:,n1:n1i,j)
!
      case default
        print*, "bc_per_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_per_z
!***********************************************************************
    subroutine bc_sym_x(f,sgn,topbot,j,rel,val)
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  11-nov-02/wolf: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mcom), optional :: val
      integer :: sgn,i,j
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) f(l1,:,:,j)=val(j)
        if (relative) then
          do i=1,nghost; f(l1-i,:,:,j)=2*f(l1,:,:,j)+sgn*f(l1+i,:,:,j); enddo
        else
          do i=1,nghost; f(l1-i,:,:,j)=              sgn*f(l1+i,:,:,j); enddo
          if (sgn<0) f(l1,:,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(l2,:,:,j)=val(j)
        if (relative) then
          do i=1,nghost; f(l2+i,:,:,j)=2*f(l2,:,:,j)+sgn*f(l2-i,:,:,j); enddo
        else
          do i=1,nghost; f(l2+i,:,:,j)=              sgn*f(l2-i,:,:,j); enddo
          if (sgn<0) f(l2,:,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case default
        print*, "bc_sym_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_sym_x
!***********************************************************************
    subroutine bc_cpc_x(f,topbot,j)
!
!  This condition gives A"+A'/R=0.
!  We compute the A1 point using a 2nd-order formula,
!  i.e. A1 = - (1-dx/2R)*A_(-1)/(1+x/2R).
!  Next, we compute A2 using a 4th-order formula,
!  and finally A3 using a 6th-order formula.
!  this can not be used in the setup for -a ..a with cpc on both sides,
!  for both sides A=0 on the boundary does for example not allow a constant Bz
!  removed this restriction in cpp
!  note that for A!=0 boundary conditions for Aphi and Az are not the same,
!  hence cpz
!
!  11-nov-09/axel+koen: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (my,mz) :: extra1,extra2
      integer :: i,j
      real :: dxR
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        dxR=-dx/x(l1)
        i=-0; f(l2+i,:,:,j)=0.
        i=-1; f(l2+i,:,:,j)=-(1.-.5*dxR)*f(l2-i,:,:,j)/(1.+.5*dxR)
        extra1=(1.+.5*dxR)*f(l2+i,:,:,j)+(1.-.5*dxR)*f(l2-i,:,:,j)
        i=-2; f(l2+i,:,:,j)=(-(1.-   dxR)*f(l2-i,:,:,j)+16.*extra1)/(1.+dxR)
        extra2=(1.+dxR)*f(l2+i,:,:,j)+(1.-dxR)*f(l2-i,:,:,j)-10.*extra1
        i=-3; f(l2+i,:,:,j)=(-(2.-3.*dxR)*f(l2-i,:,:,j)+27.*extra2)/(2.+3.*dxR)
!
      case ('top')               ! top boundary
        dxR=-dx/x(l2)
        i=0; f(l2+i,:,:,j)=0.
        i=1; f(l2+i,:,:,j)=-(1.-.5*dxR)*f(l2-i,:,:,j)/(1.+.5*dxR)
        extra1=(1.+.5*dxR)*f(l2+i,:,:,j)+(1.-.5*dxR)*f(l2-i,:,:,j)
        i=2; f(l2+i,:,:,j)=(-(1.-   dxR)*f(l2-i,:,:,j)+16.*extra1)/(1.+dxR)
        extra2=(1.+dxR)*f(l2+i,:,:,j)+(1.-dxR)*f(l2-i,:,:,j)-10.*extra1
        i=3; f(l2+i,:,:,j)=(-(2.-3.*dxR)*f(l2-i,:,:,j)+27.*extra2)/(2.+3.*dxR)
!
      case default
        print*, "bc_cpc_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_cpc_x
!***********************************************************************
    subroutine bc_cpz_x(f,topbot,j)
!
!  This condition gives R(RA)"-(RA)'=0, i e perfect conductor condition
!  for Az in cylindrical coordinates.
!  We compute the A1 point using a 2nd-order formula,
!  Next, we compute A2 using a 4th-order formula,
!  and finally A3 using a 6th-order formula.
!
!  28-feb-11/koen: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (my,mz) :: f1_co,f2_co
      integer :: i,j
      real :: dxR
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        dxR=dx/x(l1)
        i=-1; f(l1+i,:,:,j)=(f(l1,:,:,j)*2+f(l1-i,:,:,j)*(dxR/2-1))/(dxR/2+1)
        f1_co=(1+dxR/2)*f(l1+i,:,:,j)+(1-dxR/2)*f(l1-i,:,:,j)
        i=-2; f(l1+i,:,:,j)=(30*f(l1,:,:,j)+16*f1_co+(dxR-1)*f(l1-i,:,:,j))/(dxR+1)
        f2_co=(1+dxR)*f(l1+i,:,:,j)+(1-dxR)*f(l1-i,:,:,j)
        i=-3; f(l1+i,:,:,j)=(490*f(l1,:,:,j)-270*f1_co+27*f2_co+(3*dxR-2)*f(l1-i,:,:,j))/(3*dxR+2)
!
      case ('top')               ! top boundary
        dxR=dx/x(l2)
        i=1; f(l2+i,:,:,j)=(f(l2,:,:,j)*(2+dxR**2)+f(l2-i,:,:,j)*(dxR/2-1))/(dxR/2+1)
        f1_co=(1+dxR/2)*f(l2+i,:,:,j)+(1-dxR/2)*f(l2-i,:,:,j)
        i=2; f(l2+i,:,:,j)=(30*f(l2,:,:,j)+16*f1_co+(dxR-1)*f(l2-i,:,:,j))/(dxR+1)
        f2_co=(1+dxR)*f(l2+i,:,:,j)+(1-dxR)*f(l2-i,:,:,j)
        i=3; f(l2+i,:,:,j)=(490*f(l2,:,:,j)-270*f1_co+27*f2_co+(3*dxR-2)*f(l2-i,:,:,j))/(3*dxR+2)
!
      case default
        print*, "bc_cpz_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_cpz_x
!***********************************************************************
    subroutine bc_cpp_x(f,topbot,j)
!
!  This condition gives RA"+A'=0, i e perfect conductor condition
!  for Aphi in cylindrical coordinates.
!  We compute the A1 point using a 2nd-order formula,
!  i.e. A1 = - (1-dx/2R)*A_(-1)/(1+x/2R).
!  Next, we compute A2 using a 4th-order formula,
!  and finally A3 using a 6th-order formula.
!
!  28-feb-11/koen: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (my,mz) :: f1_co,f2_co
      integer :: i,j
      real :: dxR
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        dxR=dx/x(l1)
        i=-1; f(l1+i,:,:,j)=(f(l1,:,:,j)*(2+dxR**2)+f(l1-i,:,:,j)*(dxR/2-1))/(dxR/2+1)
        f1_co=(1+dxR/2)*f(l1+i,:,:,j)+(1-dxR/2)*f(l1-i,:,:,j)
        i=-2; f(l1+i,:,:,j)=((30+12*dxR**2)*f(l1,:,:,j)+16*f1_co+(dxR-1)*f(l1-i,:,:,j))/(dxR+1)
        f2_co=(1+dxR)*f(l1+i,:,:,j)+(1-dxR)*f(l1-i,:,:,j)
        i=-3; f(l1+i,:,:,j)=((490+180*dxR**2)*f(l1,:,:,j)-270*f1_co+27*f2_co+(3*dxR-2)*f(l1-i,:,:,j))/(3*dxR+2)
!
      case ('top')               ! top boundary
        dxR=dx/x(l2)
        i=1; f(l2+i,:,:,j)=(f(l2,:,:,j)*(2+dxR**2)+f(l2-i,:,:,j)*(dxR/2-1))/(dxR/2+1)
        f1_co=(1+dxR/2)*f(l2+i,:,:,j)+(1-dxR/2)*f(l2-i,:,:,j)
        i=2; f(l2+i,:,:,j)=((30+12*dxR**2)*f(l2,:,:,j)+16*f1_co+(dxR-1)*f(l2-i,:,:,j))/(dxR+1)
        f2_co=(1+dxR)*f(l2+i,:,:,j)+(1-dxR)*f(l2-i,:,:,j)
        i=3; f(l2+i,:,:,j)=((490+180*dxR**2)*f(l2,:,:,j)-270*f1_co+27*f2_co+(3*dxR-2)*f(l2-i,:,:,j))/(3*dxR+2)
!
      case default
        print*, "bc_cpp_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_cpp_x
!***********************************************************************
    subroutine bc_symset_x(f,sgn,topbot,j,rel,val)
!
!  This routine works like bc_sym_x, but sets the function value to val
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  11-nov-02/wolf: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mcom), optional :: val
      integer :: sgn,i,j
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) f(l1,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(l1-i,:,:,j)=2*f(l1,:,:,j)+sgn*f(l1+i,:,:,j); enddo
        else
          do i=1,nghost; f(l1-i,:,:,j)=              sgn*f(l1+i,:,:,j); enddo
          f(l1,:,:,j)=(4.*f(l1+1,:,:,j)-f(l1+2,:,:,j))/3.
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(l2,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(l2+i,:,:,j)=2*f(l2,:,:,j)+sgn*f(l2-i,:,:,j); enddo
        else
          do i=1,nghost; f(l2+i,:,:,j)=              sgn*f(l2-i,:,:,j); enddo
          f(l2,:,:,j)=(4.*f(l2-1,:,:,j)-f(l2-2,:,:,j))/3.
        endif
!
      case default
        print*, "bc_symset_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_symset_x
!***********************************************************************
    subroutine bc_symset0der_x(f,topbot,j)
!
!  This routine works like bc_sym_x, but sets the function value to what
!  it should be for vanishing one-sided derivative.
!  This is the routine to be used as regularity condition on the axis.
!
!  12-nov-09/axel+koen: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j,i1=1,i2=2,i3=3,i4=4,i5=5,i6=6
!
      select case (topbot)
!
!  bottom (left end of the domain)
!
      case ('bot')               ! bottom boundary
        f(l1,m1:m2,n1:n2,j)=(360.*f(l1+i1,m1:m2,n1:n2,j) &
                            -450.*f(l1+i2,m1:m2,n1:n2,j) &
                            +400.*f(l1+i3,m1:m2,n1:n2,j) &
                            -225.*f(l1+i4,m1:m2,n1:n2,j) &
                             +72.*f(l1+i5,m1:m2,n1:n2,j) &
                             -10.*f(l1+i6,m1:m2,n1:n2,j))/147.
        do i=1,nghost; f(l1-i,:,:,j)=f(l1+i,:,:,j); enddo
!
!  top (right end of the domain)
!
      case ('top')               ! top boundary
        f(l2,m1:m2,n1:n2,j)=(360.*f(l2-i1,m1:m2,n1:n2,j) &
                            -450.*f(l2-i2,m1:m2,n1:n2,j) &
                            +400.*f(l2-i3,m1:m2,n1:n2,j) &
                            -225.*f(l2-i4,m1:m2,n1:n2,j) &
                             +72.*f(l2-i5,m1:m2,n1:n2,j) &
                             -10.*f(l2-i6,m1:m2,n1:n2,j))/147.
        do i=1,nghost; f(l2+i,:,:,j)=f(l2-i,:,:,j); enddo
!
      case default
        print*, "bc_symset0der_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_symset0der_x
!***********************************************************************
    subroutine bc_slope_x(f,slope,topbot,j,rel,val)
!
! FIXME: This documentation is almost certainly wrong
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  25-feb-07/axel: adapted from bc_sym_x
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mcom), optional :: val
      real, dimension (mcom) :: slope
      integer :: i,j
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) f(l1,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost
            f(l1-i,:,:,j)=2*f(l1,:,:,j)+slope(j)*f(l1+i,:,:,j)*x(l1+i)/x(l1-i)
          enddo
        else
          do i=1,nghost
            f(l1-i,:,:,j)=f(l1+i,:,:,j)*(x(l1+i)/x(l1-i))**slope(j)
          enddo
!         f(l1,:,:,j)=(2.*x(l1+1)*f(l1+1,:,:,j)-&
!          .5*x(l1+2)*f(l1+2,:,:,j))/(1.5*x(l1))
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(l2,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost
            f(l2+i,:,:,j)=2*f(l2,:,:,j)+slope(j)*f(l2-i,:,:,j)
          enddo
        else
          do i=1,nghost
            f(l2+i,:,:,j)=f(l2-i,:,:,j)*(x(l2-i)/x(l2+i))**slope(j)
          enddo
!         f(l2,:,:,j)=(2.*x(l2-1)*f(l2-1,:,:,j)-&
!           .5*x(l2-2)*f(l2-2,:,:,j))/(1.5*x(l2))
        endif
!
      case default
        print*, "bc_slope_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_slope_x
!***********************************************************************
    subroutine bc_dr0_x(f,slope,topbot,j,rel,val)
!
! FIXME: This documentation is almost certainly wrong
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  25-feb-07/axel: adapted from bc_sym_x
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mcom), optional :: val
      real, dimension (mcom) :: slope
      integer :: i,j
      ! Abbreviations to keep compiler from complaining in 1-d or 2-d:
      integer :: l1_4, l1_5, l1_6
      integer :: l2_4, l2_5, l2_6
      logical, optional :: rel
      logical :: relative
!
      l1_4=l1+4; l1_5=l1+5; l1_6=l1+6
      l2_4=l2-4; l2_5=l2-5; l2_6=l2-6
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) f(l1,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost
            f(l1-i,:,:,j)=2*f(l1,:,:,j)+slope(j)*f(l1+i,:,:,j)*x(l1+i)/x(l1-i)
          enddo
        else
          f(l1,:,:,j)=(360.*x(l1+1)*f(l1+1,:,:,j)-450.*x(l1+2)*f(l1+2,:,:,j) &
                      +400.*x(l1+3)*f(l1+3,:,:,j)-225.*x(l1_4)*f(l1_4,:,:,j) &
                       +72.*x(l1_5)*f(l1_5,:,:,j)- 10.*x(l1_6)*f(l1_6,:,:,j) &
                      )/(147.*x(l1))
          do i=1,nghost
            f(l1-i,:,:,j)=f(l1+i,:,:,j)+(2.*dx/x(l1))*i*f(l1,:,:,j)
          enddo
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(l2,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost
            f(l2+i,:,:,j)=2*f(l2,:,:,j)+slope(j)*f(l2-i,:,:,j)
          enddo
        else
          f(l2,:,:,j)=(360.*x(l2-1)*f(l2-1,:,:,j)-450.*x(l2-2)*f(l2-2,:,:,j) &
                      +400.*x(l2-3)*f(l2-3,:,:,j)-225.*x(l2_4)*f(l2_4,:,:,j) &
                       +72.*x(l2_5)*f(l2_5,:,:,j)- 10.*x(l2_6)*f(l2_6,:,:,j) &
                      )/(147.*x(l2))
          do i=1,nghost
            f(l2+i,:,:,j)=f(l2-i,:,:,j)-(2.*dx/x(l2))*i*f(l2,:,:,j)
          enddo
        endif
!
      case default
        print*, "bc_slope_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_dr0_x
!***********************************************************************
    subroutine bc_overshoot_x(f,dist,topbot,j)
!
!  Overshoot boundary conditions, ie (d/dx-1/dist) f = 0.
!  Is implemented as d/dx [ f*exp(-x/dist) ] = 0,
!  so f(l1-i)*exp[-x(l1-i)/dist] = f(l1+i)*exp[-x(l1+i)/dist],
!  or f(l1-i) = f(l1+i)*exp{[x(l1-i)-x(l1+i)]/dist}.
!
!  25-feb-07/axel: adapted from bc_sym_x
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mcom) :: dist
      integer :: i,j
!
      select case (topbot)
!
!  bottom
!
      case ('bot')               ! bottom boundary
        do i=1,nghost
          f(l1-i,:,:,j)=f(l1+i,:,:,j)*exp((x(l1-i)-x(l1+i))/dist(j))
        enddo
!
!  top
!
      case ('top')               ! top boundary
        do i=1,nghost
          f(l2+i,:,:,j)=f(l2-i,:,:,j)*exp((x(l2+i)-x(l2-i))/dist(j))
        enddo
!
!  default
!
      case default
        print*, "bc_overshoot_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_overshoot_x
!***********************************************************************
    subroutine bc_overshoot_z(f,dist,topbot,j)
!
!  Overshoot boundary conditions, ie (d/dz-1/dist) f = 0.
!  Is implemented as d/dz [ f*exp(-z/dist) ] = 0,
!  so f(n1-i)*exp[-z(n1-i)/dist] = f(n1+i)*exp[-z(n1+i)/dist],
!  or f(n1-i) = f(n1+i)*exp{[z(n1-i)-z(n1+i)]/dist}.
!
!  25-feb-07/axel: adapted from bc_sym_z
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mcom) :: dist
      integer :: i,j
!
      select case (topbot)
!
!  bottom
!
      case ('bot')               ! bottom boundary
        do i=1,nghost
          f(:,:,n1-i,j)=f(:,:,n1+i,j)*exp((z(n1-i)-z(n1+i))/dist(j))
        enddo
!
!  top
!
      case ('top')               ! top boundary
        do i=1,nghost
          f(:,:,n2+i,j)=f(:,:,n2-i,j)*exp((z(n2+i)-z(n2-i))/dist(j))
        enddo
!
!  default
!
      case default
        print*, "bc_overshoot_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_overshoot_z
!***********************************************************************
    subroutine bc_antis_x(f,slope,topbot,j,rel,val)
!
!  Print a warning to prompt potential users to document this.
!  This routine seems an experimental one to me (Axel)
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  25-feb-07/axel: adapted from bc_slope_x
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mcom), optional :: val
      real, dimension (mcom) :: slope
      integer :: i,j
      logical, optional :: rel
      logical :: relative
!
!  Print a warning to prompt potential users to document this.
!
      call fatal_error('bc_antis_x','outdated/invalid? Document if needed')
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) f(l1,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost
            f(l1-i,:,:,j)=2*f(l1,:,:,j)+slope(j)*f(l1+i,:,:,j)*x(l1+i)/x(l1-i)
          enddo
        else
          f(l1,:,:,j)=0.
          do i=1,nghost
            f(l1-i,:,:,j)=-f(l1+i,:,:,j)*(x(l1+i)/x(l1-i))**slope(j)
          enddo
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(l2,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost
            f(l2+i,:,:,j)=2*f(l2,:,:,j)+slope(j)*f(l2-i,:,:,j)
          enddo
        else
          f(l2,:,:,j)=0.
          do i=1,nghost
            f(l2+i,:,:,j)=-f(l2-i,:,:,j)*(x(l2-i)/x(l2+i))**slope(j)
          enddo
        endif
!
      case default
        print*, "bc_antis_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_antis_x
!***********************************************************************
    subroutine bc_sym_y(f,sgn,topbot,j,rel,val)
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  11-nov-02/wolf: coded
!  10-apr-05/axel: added val argument
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mcom), optional :: val
      integer :: sgn,i,j
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) f(:,m1,:,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,m1-i,:,j)=2*f(:,m1,:,j)+sgn*f(:,m1+i,:,j); enddo
        else
          do i=1,nghost; f(:,m1-i,:,j)=              sgn*f(:,m1+i,:,j); enddo
          if (sgn<0) f(:,m1,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(:,m2,:,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,m2+i,:,j)=2*f(:,m2,:,j)+sgn*f(:,m2-i,:,j); enddo
        else
          do i=1,nghost; f(:,m2+i,:,j)=              sgn*f(:,m2-i,:,j); enddo
          if (sgn<0) f(:,m2,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case default
        print*, "bc_sym_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_sym_y
!***********************************************************************
    subroutine bc_symset_y(f,sgn,topbot,j,rel,val)
!
!  This routine works like bc_sym_x, but sets the function value to what
!  it should be for vanishing one-sided derivative.
!  At the moment the derivative is only 2nd order accurate.
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  11-nov-02/wolf: coded
!  10-apr-05/axel: added val argument
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mcom), optional :: val
      integer :: sgn,i,j
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) f(l1:l2,m1,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,m1-i,:,j)=2*f(:,m1,:,j)+sgn*f(:,m1+i,:,j); enddo
        else
          do i=1,nghost; f(:,m1-i,:,j)=              sgn*f(:,m1+i,:,j); enddo
          f(:,m1,:,j)=(4.*f(:,m1+1,:,j)-f(:,m1+2,:,j))/3.
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(l1:l2,m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,m2+i,:,j)=2*f(:,m2,:,j)+sgn*f(:,m2-i,:,j); enddo
        else
          do i=1,nghost; f(:,m2+i,:,j)=              sgn*f(:,m2-i,:,j); enddo
          f(:,m2,:,j)=(4.*f(:,m2-1,:,j)-f(:,m2-2,:,j))/3.
        endif
!
      case default
        print*, "bc_symset_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_symset_y
!***********************************************************************
    subroutine bc_symset0der_y(f,topbot,j)
!
!  This routine works like bc_sym_y, but sets the function value to what
!  it should be for vanishing one-sided derivative.
!  This is the routine to be used as regularity condition on the axis.
!
!  19-nov-09/axel: adapted from bc_symset0der_x
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j,i1=1,i2=2,i3=3,i4=4,i5=5,i6=6
!
      select case (topbot)
!
!  bottom (left end of the domain)
!
      case ('bot')               ! bottom boundary
        f(:,m1,:,j)=(360.*f(:,m1+i1,:,j) &
                    -450.*f(:,m1+i2,:,j) &
                    +400.*f(:,m1+i3,:,j) &
                    -225.*f(:,m1+i4,:,j) &
                     +72.*f(:,m1+i5,:,j) &
                     -10.*f(:,m1+i6,:,j))/147.
        do i=1,nghost; f(:,m1-i,:,j)=f(:,m1+i,:,j); enddo
!
!  top (right end of the domain)
!
      case ('top')               ! top boundary
        f(:,m2,:,j)=(360.*f(:,m2-i1,:,j) &
                    -450.*f(:,m2-i2,:,j) &
                    +400.*f(:,m2-i3,:,j) &
                    -225.*f(:,m2-i4,:,j) &
                     +72.*f(:,m2-i5,:,j) &
                     -10.*f(:,m2-i6,:,j))/147.
        do i=1,nghost; f(:,m2+i,:,j)=f(:,m2-i,:,j); enddo
!
      case default
        print*, "bc_symset0der_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_symset0der_y
!***********************************************************************
    subroutine bc_sym_z(f,sgn,topbot,j,rel,val)
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  11-nov-02/wolf: coded
!  10-apr-05/axel: added val argument
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mcom), optional :: val
      integer :: sgn,i,j
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) f(:,:,n1,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,:,n1-i,j)=2*f(:,:,n1,j)+sgn*f(:,:,n1+i,j); enddo
        else
          do i=1,nghost; f(:,:,n1-i,j)=              sgn*f(:,:,n1+i,j); enddo
          if (sgn<0) f(:,:,n1,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(:,:,n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,:,n2+i,j)=2*f(:,:,n2,j)+sgn*f(:,:,n2-i,j); enddo
        else
          do i=1,nghost; f(:,:,n2+i,j)=              sgn*f(:,:,n2-i,j); enddo
          if (sgn<0) f(:,:,n2,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case default
        print*, "bc_sym_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_sym_z
!***********************************************************************
    subroutine bc_symset0der_z(f,topbot,j)
!
!  This routine works like bc_sym_z, but sets the function value to what
!  it should be for vanishing one-sided derivative.
!  This is the routine to be used as regularity condition on the axis.
!
!  22-nov-09/axel: adapted from bc_symset0der_y
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j,i1=1,i2=2,i3=3,i4=4,i5=5,i6=6
!
      select case (topbot)
!
!  bottom (left end of the domain)
!
      case ('bot')               ! bottom boundary
        f(:,:,n1,j)=(360.*f(:,:,n1+i1,j) &
                    -450.*f(:,:,n1+i2,j) &
                    +400.*f(:,:,n1+i3,j) &
                    -225.*f(:,:,n1+i4,j) &
                     +72.*f(:,:,n1+i5,j) &
                     -10.*f(:,:,n1+i6,j))/147.
        do i=1,nghost; f(:,:,n1-i,j)=f(:,:,n1+i,j); enddo
!
!  top (right end of the domain)
!
      case ('top')               ! top boundary
        f(:,:,n2,j)=(360.*f(:,:,n2-i1,j) &
                    -450.*f(:,:,n2-i2,j) &
                    +400.*f(:,:,n2-i3,j) &
                    -225.*f(:,:,n2-i4,j) &
                     +72.*f(:,:,n2-i5,j) &
                     -10.*f(:,:,n2-i6,j))/147.
        do i=1,nghost; f(:,:,n2+i,j)=f(:,:,n2-i,j); enddo
!
      case default
        print*, "bc_symset0der_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_symset0der_z
!***********************************************************************
    subroutine bc_set_der_x(f,topbot,j,val)
!
!  Sets the derivative on the boundary to a given value.
!
!  14-may-2006/tobi: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
      real, intent (in) :: val
!
      integer :: i
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost; f(l1-i,:,:,j) = f(l1+i,:,:,j) - 2*i*dx*val; enddo
!
      case ('top')               ! top boundary
        do i=1,nghost; f(l2+i,:,:,j) = f(l2-i,:,:,j) + 2*i*dx*val; enddo
!
      case default
        call warning('bc_set_der_x',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_der_x
!***********************************************************************
    subroutine bc_fix_x(f,topbot,j,val)
!
!  Sets the value of f, particularly:
!    A_{\alpha}= <val>
!  on the boundary to a given value
!
!  27-apr-2007/dhruba: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
!
      real, intent (in) :: val
      integer :: i
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost; f(l1-i,:,:,j)=val; enddo
      case ('top')               ! top boundary
        do i=1,nghost; f(l2+i,:,:,j)=val; enddo
      case default
        call warning('bc_fix_x',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_fix_x
!***********************************************************************
    subroutine bc_file_x(f,topbot,j)
!
!  Sets the value of f from a file
!
!   9-jan-2008/axel+nils+natalia: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
!
      real, dimension (:,:,:,:), allocatable :: bc_file_x_array
      integer :: i,lbc0,lbc1,lbc2,stat
      real :: lbc,frac
      logical, save :: lbc_file_x=.true.
!
!  Allocate memory for large array.
!
      allocate(bc_file_x_array(mx,my,mz,mvar),stat=stat)
      if (stat>0) call fatal_error('bc_file_x', &
          'Could not allocate memory for bc_file_x_array')
!
      if (lbc_file_x) then
        if (lroot) then
          print*,'opening bc_file_x.dat'
          open(9,file=trim(directory_snap)//'/bc_file_x.dat',form='unformatted')
          read(9,end=99) bc_file_x_array
          close(9)
        endif
        lbc_file_x=.false.
      endif
!
      select case (topbot)
!
!  x - Udrift_bc*t = dx * (ix - Udrift_bc*t/dx)
!
      case ('bot')               ! bottom boundary
        lbc=Udrift_bc*t*dx_1(1)+1.
        lbc0=int(lbc)
        frac=mod(lbc,real(lbc0))
        lbc1=mx+mod(-lbc0,mx)
        lbc2=mx+mod(-lbc0-1,mx)
        do i=1,nghost
          f(l1-i,:,:,j)=(1-frac)*bc_file_x_array(lbc1,:,:,j) &
                           +frac*bc_file_x_array(lbc2,:,:,j)
        enddo
      case ('top')               ! top boundary
!
!  note: this "top" thing hasn't been adapted or tested yet.
!  The -lbc0-1 has been changed to +lbc0+1, but has not been tested yet.
!
        lbc=Udrift_bc*t*dx_1(1)+1.
        lbc0=int(lbc)
        frac=mod(lbc,real(lbc0))
        lbc1=mx+mod(+lbc0,mx)
        lbc2=mx+mod(+lbc0+1,mx)
        do i=1,nghost
          f(l2+i,:,:,j)=(1-frac)*bc_file_x_array(lbc1,:,:,j) &
                           +frac*bc_file_x_array(lbc2,:,:,j)
        enddo
      case default
        call warning('bc_fix_x',topbot//" should be 'top' or 'bot'")
!
      endselect
!
      goto 98
99    continue
      if (lroot) print*,'need file with dimension: ',mx,my,mz,mvar
!
      call stop_it("boundary file bc_file_x.dat not found")
!
!  Deallocate array
!
98    if (allocated(bc_file_x_array)) deallocate(bc_file_x_array)
!
    endsubroutine bc_file_x
!***********************************************************************
    subroutine bc_set_spder_x(f,topbot,j,val)
!
!  Sets the derivative, particularly:
!    d(rA_{\alpha})/dr = <val>
!  on the boundary to a given value
!
!  27-apr-2007/dhruba: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
!
      real, intent (in) :: val
      integer :: i
!
      if (lspherical_coords)then
        select case (topbot)
        case ('bot')               ! bottom boundary
        do i=1,nghost
          f(l1-i,:,:,j)=f(l1+i,:,:,j)-2*i*dx*(val-f(l1,:,:,j)*r1_mn(1))
        enddo
      case ('top')               ! top boundary
        do i=1,nghost
          f(l2+i,:,:,j)=f(l2-i,:,:,j)+2*i*dx*(val-f(l2,:,:,j)*r1_mn(nx))
        enddo
!
      case default
        call warning('bc_set_spder_x',topbot//" should be 'top' or 'bot'")
!
      endselect
    else
      call stop_it('bc_set_spder_x valid only in spherical coordinate system')
    endif
!
    endsubroutine bc_set_spder_x
! **********************************************************************
    subroutine bc_set_pfc_x(f,topbot,j)
!
!  In spherical polar coordinate system,
!  at a radial boundary set : $A_{\theta} = 0$ and $A_{phi} = 0$,
!  and demand $div A = 0$ gives the condition on $A_r$ to be
!  $d/dr( A_r) + 2/r = 0$ . This subroutine sets this condition of
!  $j$ the component of f. As this is related to setting the
!  perfect conducting boundary condition we call this "pfc".
!
!  25-Aug-2007/dhruba: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
! The coding assumes we are using 6-th order centered finite difference for our
! derivatives.
        f(l1-1,:,:,j)= f(l1+1,:,:,j) +  2.*60.*f(l1,:,:,j)*dx/(45.*x(l1))
        f(l1-2,:,:,j)= f(l1+2,:,:,j) +  2.*60.*f(l1,:,:,j)*dx/(9.*x(l1))
        f(l1-3,:,:,j)= f(l1+3,:,:,j) +  2.*60.*f(l1,:,:,j)*dx/x(l1)
      case ('top')               ! top boundary
        f(l2+1,:,:,j)= f(l2-1,:,:,j) -  2.*60.*f(l2,:,:,j)*dx/(45.*x(l2))
        f(l2+2,:,:,j)= f(l2-2,:,:,j) -  2.*60.*f(l2,:,:,j)*dx/(9.*x(l2))
        f(l2+3,:,:,j)= f(l2-3,:,:,j) -  2.*60.*f(l2,:,:,j)*dx/(x(l2))
!
      case default
        call warning('bc_set_pfc_x',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_pfc_x
!***********************************************************************
    subroutine bc_set_nfr_x(f,topbot,j)
!
!  Normal-field (or angry-hedgehog) boundary condition for spherical
!  coordinate system.
!  d_r(A_{\theta}) = -A_{\theta}/r  with A_r = 0 sets B_{r} to zero
!  in spherical coordinate system.
!  (compare with next subroutine sfree )
!
!  25-Aug-2007/dhruba: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
      integer :: k
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do k=1,nghost
          f(l1-k,:,:,j)= f(l1+k,:,:,j)*(x(l1+k)/x(l1-k))
        enddo
!
     case ('top')               ! top boundary
       do k=1,nghost
         f(l2+k,:,:,j)= f(l2-k,:,:,j)*(x(l2-k)/x(l2+k))
       enddo
!
      case default
        call warning('bc_set_nfr_x',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_nfr_x
! **********************************************************************
    subroutine bc_set_sa2_x(f,topbot,j)
!
!  To set the boundary condition:
!  d_r(r B_{\phi} = 0 we need to se
!  (d_r)^2(r A_{\theta}) = 0 which sets the condition 'a2'
!  on r A_{\theta} and vice-versa for A_{\phi}
!
!  03-Dec-2009/dhruba: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
      integer :: k
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do k=1,nghost
          f(l1-k,:,:,j)= f(l1,:,:,j)*2.*(x(l1)/x(l1-k))&
                         -f(l1+k,:,:,j)*(x(l1+k)/x(l1-k))
        enddo
!
     case ('top')               ! top boundary
       do k=1,nghost
         f(l2+k,:,:,j)= f(l2,:,:,j)*2.*(x(l2)/x(l2+k))&
                        -f(l2-k,:,:,j)*(x(l2-k)/x(l2+k))
       enddo
!
      case default
        call warning('bc_set_sa2_x',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_sa2_x
! **********************************************************************
    subroutine bc_set_sfree_x(f,topbot,j)
!
!  Details are given in an appendix in the manual.
!  Lambda effect : stresses due to Lambda effect are added to the stress-tensor.
!  For rotation along the z direction and also for not very strong rotation such
!  that the breaking of rotational symmetry is only due to gravity, the only
!  new term is appears in the r-phi component. This implies that this term
!  affects only the boundary condition of u_{\phi} for the radial boundary.
!
!  25-Aug-2007/dhruba: coded
!  21-Mar-2009/axel: get llambda_effect using get_shared_variable
!
      use SharedVariables, only : get_shared_variable
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
!
      real, pointer :: nu,Lambda_V0t,Lambda_V0b,Lambda_V1t,Lambda_V1b
      logical, pointer :: llambda_effect
      integer :: ierr,k
      real :: lambda_exp,lambda_exp_sinth
! -------- Either case get the lambda variables first -----------
!
      call get_shared_variable('nu',nu,ierr)
      if (ierr/=0) call stop_it("bc_set_sfree_x: "//&
          "there was a problem when getting nu")
      call get_shared_variable('llambda_effect',llambda_effect,ierr)
      if (ierr/=0) call stop_it("bc_set_sfree_x: "//&
          "there was a problem when getting llambda_effect")
      if (llambda_effect) then
         call get_shared_variable('Lambda_V0t',Lambda_V0t,ierr)
         if (ierr/=0) call stop_it("bc_set_sfree_x: " &
              // "problem getting shared var Lambda_V0t")
         call get_shared_variable('Lambda_V1t',Lambda_V1t,ierr)
         if (ierr/=0) call stop_it("bc_set_sfree_x: " &
              // "problem getting shared var Lambda_V1t")
         call get_shared_variable('Lambda_V0b',Lambda_V0b,ierr)
         if (ierr/=0) call stop_it("bc_set_sfree_x: " &
              // "problem getting shared var Lambda_V0b")
         call get_shared_variable('Lambda_V1b',Lambda_V1b,ierr)
         if (ierr/=0) call stop_it("bc_set_sfree_x: " &
              // "problem getting shared var Lambda_V1b")
      endif
!
      select case (topbot)
!
! Bottom boundary
!
      case ('bot')
!
        if ((llambda_effect).and.(j==iuz)) then
          do iy=1,my
            lambda_exp=-(Lambda_V0b+Lambda_V1b*sinth(iy)*sinth(iy))
            lambda_exp_sinth=lambda_exp*sinth(iy)
            do k=1,nghost
               if (Omega==0) then
                 f(l1-k,iy,:,j)= f(l1+k,iy,:,j)*(x(l1-k)/x(l1+k)) &
                     **(1-(lambda_exp/nu))
               else
                 f(l1-k,iy,:,j)= (f(l1+k,iy,:,j)+Omega*x(l1+k)*sinth(iy)) &
                      *(x(l1-k)/x(l1+k))**(1-(lambda_exp/nu)) &
                      -Omega*x(l1-k)*sinth(iy)
               endif
            enddo
          enddo
        else
          do k=1,nghost
            f(l1-k,:,:,j)= f(l1+k,:,:,j)*(x(l1-k)/x(l1+k))
          enddo
        endif
!
! Top boundary
!
      case ('top')
        if ((llambda_effect).and.(j==iuz)) then
          do iy=1,my
            lambda_exp=-(Lambda_V0t+Lambda_V1t*sinth(iy)*sinth(iy))
            lambda_exp_sinth=lambda_exp*sinth(iy)
            do k=1,nghost
              if (Omega==0) then
                f(l2+k,iy,:,j)= f(l2-k,iy,:,j)*((x(l2+k)/x(l2-k)) &
                    **(1-(lambda_exp/nu)))
              else
                f(l2+k,iy,:,j)= (f(l2-k,iy,:,j)+Omega*x(l2-k)*sinth(iy)) &
                     *(x(l2+k)/x(l2-k))**(1-(lambda_exp/nu)) &
                     -Omega*x(l2+k)*sinth(iy)
              endif
            enddo
          enddo
        else
          do k=1,nghost
            f(l2+k,:,:,j)= f(l2-k,:,:,j)*(x(l2+k)/x(l2-k))
          enddo
        endif
!
      case default
        call warning('bc_set_sfree_x',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_sfree_x
! **********************************************************************
    subroutine bc_set_jethat_x(f,jj,topbot,fracall,uzeroall)
!
!  Sets tophat velocity profile at the inner (bot) boundary
!
!  03-jan-2008/dhruba: coded
!
      use Sub, only: step
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent(in) :: jj
      integer :: i,j,k
      real, dimension(mcom),intent(in) :: fracall,uzeroall
      real :: frac,uzero,ylim,ymid,y1,zlim,zmid,z1
      real :: yhat_min,yhat_max,zhat_min,zhat_max
      real, parameter :: width_hat=0.01
      real, dimension (ny) :: hatprofy
      real, dimension (nz) :: hatprofz
!
      y1 = xyz1(2)
      z1 = xyz1(3)
      frac = fracall(jj)
      uzero = uzeroall(jj)
!
      if (lspherical_coords)then
!
        select case (topbot)
        case ('bot')               ! bottom boundary
          ylim = (y1-y0)*frac
          ymid = y0+(y1-y0)/2.
          yhat_min=ymid-ylim/2.
          yhat_max=ymid+ylim/2
          hatprofy=step(y(m1:m2),yhat_min,width_hat)*(1.-step(y(m1:m2),yhat_max,width_hat))
          zlim = (z1-z0)*frac
          zmid = z0+(z1-z0)/2.
          zhat_min=zmid-zlim/2.
          zhat_max=zmid+zlim/2
          hatprofz=step(z(n1:n2),zhat_min,width_hat)*(1.-step(z(n1:n2),zhat_max,width_hat))
          do j=m1,m2
            do k=n1,n2
                f(l1,j,k,iux)= uzero*hatprofy(j)*hatprofz(k)
                do i=1,nghost
                  f(l1-i,j,k,iux)= uzero*hatprofy(j)*hatprofz(k)
                enddo
            enddo
          enddo
!
        case ('top')               ! top boundary
          call warning('bc_set_jethat_x','Jet flowing out of the exit boundary ?')
          do i=1,nghost
            f(l2+i,:,:,jj)=0.
          enddo
!
        case default
          call warning('bc_set_jethat_x',topbot//" should be 'top' or 'bot'")
        endselect
!
      else
        call stop_it('Boundary condition jethat is valid only in spherical coordinate system')
      endif
!
    endsubroutine bc_set_jethat_x
! **********************************************************************
    subroutine bc_set_nfr_y(f,topbot,j)
!
!  Stress-free boundary condition for spherical coordinate system.
!  d_{\theta}(A_{\phi}) = -A_{\phi}cot(\theta)/r  with A_{\theta} = 0 sets
!  B_{\theta}=0 in spherical polar
!  coordinate system. This subroutine sets only the first part of this
!  boundary condition for 'j'-th component of f.
!
!  25-Aug-2007/dhruba: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
      integer :: k
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do k=1,nghost
          f(:,m1-k,:,j)= f(:,m1+k,:,j)*sinth(m1+k)*sin1th(m1-k)
        enddo
       case ('top')               ! top boundary
         do k=1,nghost
           f(:,m2+k,:,j)= f(:,m2-k,:,j)*sinth(m2-k)*sin1th(m2+k)
         enddo
!
      case default
        call warning('bc_set_nfr_y',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_nfr_y
! **********************************************************************
    subroutine bc_set_sfree_y(f,topbot,j)
!
!  Stress-free boundary condition for spherical coordinate system.
!  d_{\theta}(u_{\phi}) = u_{\phi}cot(\theta)  with u_{\theta} = 0 sets
!  S_{\theta \phi} component of the strain matrix to be zero in spherical
!  coordinate system. This subroutine sets only the first part of this
!  boundary condition for 'j'-th component of f.
!
!  25-Aug-2007/dhruba: coded
!
      use SharedVariables, only : get_shared_variable
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
      real, pointer :: Lambda_H1,nu
      real, pointer :: LH1_rprof(:)
      logical, pointer :: llambda_effect
      integer :: ierr,k,ix
      real :: cos2thm_k,cos2thmpk,somega
      real,dimension(mx):: LH1
! -------- Either case get the lambda variables first -----------
!
      call get_shared_variable('nu',nu,ierr)
      if (ierr/=0) call stop_it("bc_set_sfree_y: "//&
          "there was a problem when getting nu")
      call get_shared_variable('llambda_effect',llambda_effect,ierr)
      if (ierr/=0) call stop_it("bc_set_sfree_y: "//&
          "there was a problem when getting llambda_effect")
      if (llambda_effect) then
         call get_shared_variable('Lambda_H1',Lambda_H1,ierr)
         if (ierr/=0) call stop_it("bc_set_sfree_y: " &
              // "problem getting shared var Lambda_H1")
         call get_shared_variable('LH1_rprof',LH1_rprof,ierr)
         if (ierr/=0) call stop_it("bc_set_sfree_y: " &
              // "problem getting shared var LH1_rprof")
         LH1=Lambda_H1*LH1_rprof
      endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if ((llambda_effect).and.(j==iuz).and.(Lambda_H1/=0.)) then
          do k=1,nghost
              cos2thm_k= costh(m1-k)**2-sinth(m1-k)**2
              cos2thmpk= costh(m1+k)**2-sinth(m1+k)**2
            if (Omega==0) then
               do ix=1,mx
                  f(ix,m1-k,:,j)= f(ix,m1+k,:,j)* &
                       (exp(LH1(ix)*cos2thmpk/(4.*nu))*sin1th(m1+k)) &
                       *(exp(-LH1(ix)*cos2thm_k/(4.*nu))*sinth(m1-k))
               enddo
            else
              do ix=1,mx
                somega=x(ix)*Omega*sinth(m1-k)*( &
                   exp(2*cos2thm_k*LH1(ix)/(4.*nu))&
                        -exp((cos2thmpk+cos2thm_k)*LH1(ix)/(4.*nu)) )
                f(ix,m1-k,:,j)= f(ix,m1+k,:,j)* &
                   (exp(LH1(ix)*cos2thmpk/(4.*nu))*sin1th(m1+k)) &
                   *(exp(-LH1(ix)*cos2thm_k/(4.*nu))*sinth(m1-k)) &
                      +somega
              enddo
            endif
          enddo
        else
          do k=1,nghost
            f(:,m1-k,:,j)= f(:,m1+k,:,j)*sinth(m1-k)*sin1th(m1+k)
          enddo
        endif
      case ('top')               ! top boundary
        if ((llambda_effect).and.(j==iuz).and.(Lambda_H1/=0)) then
          do k=1,nghost
            cos2thm_k= costh(m2-k)**2-sinth(m2-k)**2
            cos2thmpk= costh(m2+k)**2-sinth(m2+k)**2
            if (Omega==0)then
               do ix=1,mx
                  f(ix,m2+k,:,j)= f(ix,m2-k,:,j)* &
                   (exp(LH1(ix)*cos2thm_k/(4.*nu))*sin1th(m2-k)) &
                  *(exp(-LH1(ix)*cos2thmpk/(4.*nu))*sinth(m2+k))
               enddo
             else
              do ix=1,mx
                somega=x(ix)*Omega*sinth(m2+k)*( &
                   exp(2*cos2thmpk*LH1(ix)/(4.*nu))&
                        -exp((cos2thmpk+cos2thm_k)*LH1(ix)/(4.*nu)) )
                f(ix,m2+k,:,j)= f(ix,m2-k,:,j)* &
                     (exp(LH1(ix)*cos2thm_k/(4.*nu))*sin1th(m2-k)) &
                    *(exp(-LH1(ix)*cos2thmpk/(4.*nu))*sinth(m2+k)) &
                      +somega
              enddo
             endif
          enddo
        else
          do k=1,nghost
            f(:,m2+k,:,j)= f(:,m2-k,:,j)*sinth(m2+k)*sin1th(m2-k)
          enddo
        endif
!
     case default
        call warning('bc_set_sfree_y',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_sfree_y
! **********************************************************************
    subroutine bc_set_pfc_y(f,topbot,j)
!
!  In spherical polar coordinate system,
!  at a theta boundary set : $A_{r} = 0$ and $A_{\phi} = 0$,
!  and demand $div A = 0$ gives the condition on $A_{\theta}$ to be
!  $d/d{\theta}( A_{\theta}) + \cot(\theta)A_{\theta} = 0$ .
!  This subroutine sets this condition on
!  $j$ the component of f. As this is related to setting the
!  perfect conducting boundary condition we call this "pfc".
!
!  25-Aug-2007/dhruba: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
      real :: cottheta
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
!
!  The coding assumes we are using 6-th order centered finite difference for our
!  derivatives.
!
        cottheta= cotth(m1)
        f(:,m1-1,:,j)= f(:,m1+1,:,j) +  60.*dy*cottheta*f(:,m1,:,j)/45.
        f(:,m1-2,:,j)= f(:,m1+2,:,j) -  60.*dy*cottheta*f(:,m1,:,j)/9.
        f(:,m1-3,:,j)= f(:,m1+3,:,j) +  60.*dy*cottheta*f(:,m1,:,j)
      case ('top')               ! top boundary
        cottheta= cotth(m2)
        f(:,m2+1,:,j)= f(:,m2-1,:,j) -  60.*dy*cottheta*f(:,m2,:,j)/45.
        f(:,m2+2,:,j)= f(:,m2-2,:,j) +  60.*dy*cottheta*f(:,m2,:,j)/9.
        f(:,m2+3,:,j)= f(:,m2-3,:,j) -  60.*dy*cottheta*f(:,m2,:,j)
!
      case default
        call warning('bc_set_pfc_y',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_pfc_y
!***********************************************************************
    subroutine bc_set_der_y(f,topbot,j,val)
!
!  Sets the derivative on the boundary to a given value.
!
!  14-may-2006/tobi: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
      real, intent (in) :: val
!
      integer :: i
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost; f(:,m1-i,:,j) = f(:,m1+i,:,j) - 2*i*dy*val; enddo
!
      case ('top')               ! top boundary
        do i=1,nghost; f(:,m2+i,:,j) = f(:,m2-i,:,j) + 2*i*dy*val; enddo
!
      case default
        call warning('bc_set_der_y',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_der_y
!***********************************************************************
    subroutine bc_set_der_z(f,topbot,j,val)
!
!  Sets the derivative on the boundary to a given value.
!
!  14-may-2006/tobi: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
      real, intent (in) :: val
!
      integer :: i
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost; f(:,:,n1-i,j) = f(:,:,n1+i,j) - 2*i*dz*val; enddo
!
      case ('top')               ! top boundary
        do i=1,nghost; f(:,:,n2+i,j) = f(:,:,n2-i,j) + 2*i*dz*val; enddo
!
      case default
        call warning('bc_set_der_z',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_der_z
!***********************************************************************
    subroutine bc_set_div_z(f,topbot,j,val)
!
!  Sets the derivative on the boundary to a given value
!
!  17-may-2010/bing: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real, dimension (nx,ny) :: fac,duz_dz
      real, intent(in) :: val
!
      integer, intent (in) :: j
!
      integer :: iref=-1,pos
!
      if (j/=iuz) call fatal_error('bc_set_div_z','only implemented for div(u)=0')
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        iref = n1
!
      case ('top')               ! top boundary
        iref = n2
!
      case default
        call warning('bc_set_der_x',topbot//" should be 'top' or 'bot'")
!
      endselect
!
! take the x derivative of ux
      if (nxgrid/=1) then
        fac=(1./60)*spread(dx_1(l1:l2),2,ny)
        duz_dz= fac*(+45.0*(f(l1+1:l2+1,m1:m2,iref,iux)-f(l1-1:l2-1,m1:m2,iref,iux)) &
            -  9.0*(f(l1+2:l2+2,m1:m2,iref,iux)-f(l1-2:l2-2,m1:m2,iref,iux)) &
            +      (f(l1+3:l2+3,m1:m2,iref,iux)-f(l1-3:l2-3,m1:m2,iref,iux)))
      else
        if (ip<=5) print*, 'bc_set_div_z: Degenerate case in x-direction'
      endif
!
! take the y derivative of uy and add to dux/dx
      if (nygrid/=1) then
        fac=(1./60)*spread(dy_1(m1:m2),1,nx)
        duz_dz=duz_dz + fac*(+45.0*(f(l1:l2,m1+1:m2+1,iref,iuy)-f(l1:l2,m1-1:m2-1,iref,iuy)) &
            -  9.0*(f(l1:l2,m1+2:m2+2,iref,iuy)-f(l1:l2,m1-2:m2-2,iref,iuy)) &
            +      (f(l1:l2,m1+3:m2+3,iref,iuy)-f(l1:l2,m1-3:m2-3,iref,iuy)))
      else
        if (ip<=5) print*, 'bc_set_div_z: Degenerate case in y-direction'
      endif
!
! add given number to set div(u)=val; default val=0
! duz/dz = val - dux/dx - duy/dy
!
      duz_dz = val - duz_dz
!
! set the derivative of uz at the boundary
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do pos=1,nghost
          f(l1:l2,m1:m2,n1-pos,j) = f(l1:l2,m1:m2,n1+pos,j) - 2*pos*dz*duz_dz
        enddo
!
      case ('top')               ! top boundary
        do pos=1,nghost
          f(l1:l2,m1:m2,n2+pos,j) = f(l1:l2,m1:m2,n2-pos,j) + 2*pos*dz*duz_dz
        enddo
!
      case default
        call warning('bc_set_div_z',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_div_z
!***********************************************************************
    subroutine bc_van_x(f,topbot,j)
!
!  Vanishing boundary conditions.
!  (TODO: clarify what this means)
!
!  26-apr-06/tobi: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
          do i=1,nghost
            f(l1-i,:,:,j)=((nghost+1-i)*f(l1,:,:,j))/(nghost+1)
          enddo
!
      case ('top')               ! top boundary
          do i=1,nghost
            f(l2+i,:,:,j)=((nghost+1-i)*f(l2,:,:,j))/(nghost+1)
          enddo
!
      case default
        print*, "bc_van_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_van_x
!***********************************************************************
    subroutine bc_van_y(f,topbot,j)
!
!  Vanishing boundary conditions.
!  (TODO: clarify what this means)
!
!  26-apr-06/tobi: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
          do i=1,nghost
            f(:,m1-i,:,j)=((nghost+1-i)*f(:,m1,:,j))/(nghost+1)
          enddo
!
      case ('top')               ! top boundary
          do i=1,nghost
            f(:,m2+i,:,j)=((nghost+1-i)*f(:,m2,:,j))/(nghost+1)
          enddo
!
      case default
        print*, "bc_van_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_van_y
!***********************************************************************
    subroutine bc_van_z(f,topbot,j)
!
!  Vanishing boundary conditions.
!  (TODO: clarify what this means)
!
!  26-apr-06/tobi: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
          do i=1,nghost
            f(:,:,n1-i,j)=((nghost+1-i)*f(:,:,n1,j))/(nghost+1)
          enddo
!
      case ('top')               ! top boundary
          do i=1,nghost
            f(:,:,n2+i,j)=((nghost+1-i)*f(:,:,n2,j))/(nghost+1)
          enddo
!
      case default
        print*, "bc_van_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_van_z
!***********************************************************************
    subroutine bc_van3rd_z(f,topbot,j)
!
!  Boundary condition with vanishing 3rd derivative
!  (useful for vertical hydrostatic equilibrium in discs)
!
!  19-aug-03/anders: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      real, dimension (:,:), allocatable :: cpoly0,cpoly1,cpoly2
      integer :: i,stat
!
!  Allocate memory for large arrays.
!
      allocate(cpoly0(mx,my),stat=stat)
      if (stat>0) call fatal_error('bc_van3rd_z', &
          'Could not allocate memory for cpoly0')
      allocate(cpoly1(mx,my),stat=stat)
      if (stat>0) call fatal_error('bc_van3rd_z', &
          'Could not allocate memory for cpoly1')
      allocate(cpoly2(mx,my),stat=stat)
      if (stat>0) call fatal_error('bc_van3rd_z', &
          'Could not allocate memory for cpoly2')
!
      select case (topbot)
!
      case ('bot')
        cpoly0(:,:)=f(:,:,n1,j)
        cpoly1(:,:)=-(3*f(:,:,n1,j)-4*f(:,:,n1+1,j)+f(:,:,n1+2,j))/(2*dz)
        cpoly2(:,:)=-(-f(:,:,n1,j)+2*f(:,:,n1+1,j)-f(:,:,n1+2,j)) /(2*dz**2)
        do i=1,nghost
          f(:,:,n1-i,j) = cpoly0(:,:) - cpoly1(:,:)*i*dz + cpoly2(:,:)*(i*dz)**2
        enddo
!
      case ('top')
        cpoly0(:,:)=f(:,:,n2,j)
        cpoly1(:,:)=-(-3*f(:,:,n2,j)+4*f(:,:,n2-1,j)-f(:,:,n2-2,j))/(2*dz)
        cpoly2(:,:)=-(-f(:,:,n2,j)+2*f(:,:,n2-1,j)-f(:,:,n2-2,j))/(2*dz**2)
        do i=1,nghost
          f(:,:,n2+i,j) = cpoly0(:,:) + cpoly1(:,:)*i*dz + cpoly2(:,:)*(i*dz)**2
        enddo
!
      endselect
!
!  Deallocate arrays.
!
      if (allocated(cpoly0)) deallocate(cpoly0)
      if (allocated(cpoly1)) deallocate(cpoly1)
      if (allocated(cpoly2)) deallocate(cpoly2)
!
    endsubroutine bc_van3rd_z
!***********************************************************************
    subroutine bc_onesided_x(f,topbot,j)
!
!  One-sided conditions.
!  These expressions result from combining Eqs(207)-(210), astro-ph/0109497,
!  corresponding to (9.207)-(9.210) in Ferriz-Mas proceedings.
!
!  05-apr-03/axel: coded
!  07-jan-09/axel: corrected
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j,k
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
          k=l1-1
          f(k,:,:,j)=7*f(k+1,:,:,j) &
                   -21*f(k+2,:,:,j) &
                   +35*f(k+3,:,:,j) &
                   -35*f(k+4,:,:,j) &
                   +21*f(k+5,:,:,j) &
                    -7*f(k+6,:,:,j) &
                      +f(k+7,:,:,j)
          k=l1-2
          f(k,:,:,j)=9*f(k+1,:,:,j) &
                   -35*f(k+2,:,:,j) &
                   +77*f(k+3,:,:,j) &
                  -105*f(k+4,:,:,j) &
                   +91*f(k+5,:,:,j) &
                   -49*f(k+6,:,:,j) &
                   +15*f(k+7,:,:,j) &
                    -2*f(k+8,:,:,j)
          k=l1-3
          f(k,:,:,j)=9*f(k+1,:,:,j) &
                   -45*f(k+2,:,:,j) &
                  +147*f(k+3,:,:,j) &
                  -315*f(k+4,:,:,j) &
                  +441*f(k+5,:,:,j) &
                  -399*f(k+6,:,:,j) &
                  +225*f(k+7,:,:,j) &
                   -72*f(k+8,:,:,j) &
                   +10*f(k+9,:,:,j)
!
      case ('top')               ! top boundary
          k=l2+1
          f(k,:,:,j)=7*f(k-1,:,:,j) &
                   -21*f(k-2,:,:,j) &
                   +35*f(k-3,:,:,j) &
                   -35*f(k-4,:,:,j) &
                   +21*f(k-5,:,:,j) &
                    -7*f(k-6,:,:,j) &
                      +f(k-7,:,:,j)
          k=l2+2
          f(k,:,:,j)=9*f(k-1,:,:,j) &
                   -35*f(k-2,:,:,j) &
                   +77*f(k-3,:,:,j) &
                  -105*f(k-4,:,:,j) &
                   +91*f(k-5,:,:,j) &
                   -49*f(k-6,:,:,j) &
                   +15*f(k-7,:,:,j) &
                    -2*f(k-8,:,:,j)
          k=l2+3
          f(k,:,:,j)=9*f(k-1,:,:,j) &
                   -45*f(k-2,:,:,j) &
                  +147*f(k-3,:,:,j) &
                  -315*f(k-4,:,:,j) &
                  +441*f(k-5,:,:,j) &
                  -399*f(k-6,:,:,j) &
                  +225*f(k-7,:,:,j) &
                   -72*f(k-8,:,:,j) &
                   +10*f(k-9,:,:,j)
!
      case default
        print*, "bc_onesided_x ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_onesided_x
!***********************************************************************
    subroutine bc_onesided_x_old(f,topbot,j)
!
!  One-sided conditions.
!  These expressions result from combining Eqs(207)-(210), astro-ph/0109497,
!  corresponding to (9.207)-(9.210) in Ferriz-Mas proceedings.
!
!  05-apr-03/axel: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j,k
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost
          k=l1-i
          f(k,:,:,j)=7*f(k+1,:,:,j) &
                   -21*f(k+2,:,:,j) &
                   +35*f(k+3,:,:,j) &
                   -35*f(k+4,:,:,j) &
                   +21*f(k+5,:,:,j) &
                    -7*f(k+6,:,:,j) &
                      +f(k+7,:,:,j)
        enddo
!
      case ('top')               ! top boundary
        do i=1,nghost
          k=l2+i
          f(k,:,:,j)=7*f(k-1,:,:,j) &
                   -21*f(k-2,:,:,j) &
                   +35*f(k-3,:,:,j) &
                   -35*f(k-4,:,:,j) &
                   +21*f(k-5,:,:,j) &
                    -7*f(k-6,:,:,j) &
                      +f(k-7,:,:,j)
        enddo
!
      case default
        print*, "bc_onesided_x_old ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_onesided_x_old
!***********************************************************************
    subroutine bc_onesided_y(f,topbot,j)
!
!  One-sided conditions.
!  These expressions result from combining Eqs(207)-(210), astro-ph/0109497,
!  corresponding to (9.207)-(9.210) in Ferriz-Mas proceedings.
!
!  05-apr-03/axel: coded
!  07-jan-09/axel: corrected
!  26-jan-09/nils: adapted from bc_onesided_x
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j,k
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
          k=m1-1
          f(:,k,:,j)=7*f(:,k+1,:,j) &
                   -21*f(:,k+2,:,j) &
                   +35*f(:,k+3,:,j) &
                   -35*f(:,k+4,:,j) &
                   +21*f(:,k+5,:,j) &
                    -7*f(:,k+6,:,j) &
                      +f(:,k+7,:,j)
          k=m1-2
          f(:,k,:,j)=9*f(:,k+1,:,j) &
                   -35*f(:,k+2,:,j) &
                   +77*f(:,k+3,:,j) &
                  -105*f(:,k+4,:,j) &
                   +91*f(:,k+5,:,j) &
                   -49*f(:,k+6,:,j) &
                   +15*f(:,k+7,:,j) &
                    -2*f(:,k+8,:,j)
          k=m1-3
          f(:,k,:,j)=9*f(:,k+1,:,j) &
                   -45*f(:,k+2,:,j) &
                  +147*f(:,k+3,:,j) &
                  -315*f(:,k+4,:,j) &
                  +441*f(:,k+5,:,j) &
                  -399*f(:,k+6,:,j) &
                  +225*f(:,k+7,:,j) &
                   -72*f(:,k+8,:,j) &
                   +10*f(:,k+9,:,j)
!
      case ('top')               ! top boundary
          k=m2+1
          f(:,k,:,j)=7*f(:,k-1,:,j) &
                   -21*f(:,k-2,:,j) &
                   +35*f(:,k-3,:,j) &
                   -35*f(:,k-4,:,j) &
                   +21*f(:,k-5,:,j) &
                    -7*f(:,k-6,:,j) &
                      +f(:,k-7,:,j)
          k=m2+2
          f(:,k,:,j)=9*f(:,k-1,:,j) &
                   -35*f(:,k-2,:,j) &
                   +77*f(:,k-3,:,j) &
                  -105*f(:,k-4,:,j) &
                   +91*f(:,k-5,:,j) &
                   -49*f(:,k-6,:,j) &
                   +15*f(:,k-7,:,j) &
                    -2*f(:,k-8,:,j)
          k=m2+3
          f(:,k,:,j)=9*f(:,k-1,:,j) &
                   -45*f(:,k-2,:,j) &
                  +147*f(:,k-3,:,j) &
                  -315*f(:,k-4,:,j) &
                  +441*f(:,k-5,:,j) &
                  -399*f(:,k-6,:,j) &
                  +225*f(:,k-7,:,j) &
                   -72*f(:,k-8,:,j) &
                   +10*f(:,k-9,:,j)
!
      case default
        print*, "bc_onesided_7 ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_onesided_y
!***********************************************************************
    subroutine bc_onesided_z_orig(f,topbot,j)
!
!  One-sided conditions.
!  These expressions result from combining Eqs(207)-(210), astro-ph/0109497,
!  corresponding to (9.207)-(9.210) in Ferriz-Mas proceedings.
!
!  05-apr-03/axel: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j,k
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost
          k=n1-i
          f(:,:,k,j)=7*f(:,:,k+1,j) &
                   -21*f(:,:,k+2,j) &
                   +35*f(:,:,k+3,j) &
                   -35*f(:,:,k+4,j) &
                   +21*f(:,:,k+5,j) &
                    -7*f(:,:,k+6,j) &
                      +f(:,:,k+7,j)
        enddo
!
      case ('top')               ! top boundary
        do i=1,nghost
          k=n2+i
          f(:,:,k,j)=7*f(:,:,k-1,j) &
                   -21*f(:,:,k-2,j) &
                   +35*f(:,:,k-3,j) &
                   -35*f(:,:,k-4,j) &
                   +21*f(:,:,k-5,j) &
                    -7*f(:,:,k-6,j) &
                      +f(:,:,k-7,j)
        enddo
!
      case default
        print*, "bc_onesided_z ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_onesided_z_orig
!***********************************************************************
    subroutine bc_onesided_z(f,topbot,j)
!
!  One-sided conditions.
!  These expressions result from combining Eqs(207)-(210), astro-ph/0109497,
!  corresponding to (9.207)-(9.210) in Ferriz-Mas proceedings.
!
!  05-apr-03/axel: coded
!  10-mar-09/axel: corrected
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j,k
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
          k=n1-1
          f(:,:,k,j)=7*f(:,:,k+1,j) &
                   -21*f(:,:,k+2,j) &
                   +35*f(:,:,k+3,j) &
                   -35*f(:,:,k+4,j) &
                   +21*f(:,:,k+5,j) &
                    -7*f(:,:,k+6,j) &
                      +f(:,:,k+7,j)
          k=n1-2
          f(:,:,k,j)=9*f(:,:,k+1,j) &
                   -35*f(:,:,k+2,j) &
                   +77*f(:,:,k+3,j) &
                  -105*f(:,:,k+4,j) &
                   +91*f(:,:,k+5,j) &
                   -49*f(:,:,k+6,j) &
                   +15*f(:,:,k+7,j) &
                    -2*f(:,:,k+8,j)
          k=n1-3
          f(:,:,k,j)=9*f(:,:,k+1,j) &
                   -45*f(:,:,k+2,j) &
                  +147*f(:,:,k+3,j) &
                  -315*f(:,:,k+4,j) &
                  +441*f(:,:,k+5,j) &
                  -399*f(:,:,k+6,j) &
                  +225*f(:,:,k+7,j) &
                   -72*f(:,:,k+8,j) &
                   +10*f(:,:,k+9,j)
!
      case ('top')               ! top boundary
          k=n2+1
          f(:,:,k,j)=7*f(:,:,k-1,j) &
                   -21*f(:,:,k-2,j) &
                   +35*f(:,:,k-3,j) &
                   -35*f(:,:,k-4,j) &
                   +21*f(:,:,k-5,j) &
                    -7*f(:,:,k-6,j) &
                      +f(:,:,k-7,j)
          k=n2+2
          f(:,:,k,j)=9*f(:,:,k-1,j) &
                   -35*f(:,:,k-2,j) &
                   +77*f(:,:,k-3,j) &
                  -105*f(:,:,k-4,j) &
                   +91*f(:,:,k-5,j) &
                   -49*f(:,:,k-6,j) &
                   +15*f(:,:,k-7,j) &
                    -2*f(:,:,k-8,j)
          k=n2+3
          f(:,:,k,j)=9*f(:,:,k-1,j) &
                   -45*f(:,:,k-2,j) &
                  +147*f(:,:,k-3,j) &
                  -315*f(:,:,k-4,j) &
                  +441*f(:,:,k-5,j) &
                  -399*f(:,:,k-6,j) &
                  +225*f(:,:,k-7,j) &
                   -72*f(:,:,k-8,j) &
                   +10*f(:,:,k-9,j)
!
      case default
        print*, "bc_onesided_z ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_onesided_z
!***********************************************************************
    subroutine bc_extrap_2_1(f,topbot,j)
!
!  Extrapolation boundary condition.
!  Correct for polynomials up to 2nd order, determined 1 further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  19-jun-03/wolf: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        f(:,:,n1-1,j)=0.25*(  9*f(:,:,n1,j)- 3*f(:,:,n1+1,j)- 5*f(:,:,n1+2,j)+ 3*f(:,:,n1+3,j))
        f(:,:,n1-2,j)=0.05*( 81*f(:,:,n1,j)-43*f(:,:,n1+1,j)-57*f(:,:,n1+2,j)+39*f(:,:,n1+3,j))
        f(:,:,n1-3,j)=0.05*(127*f(:,:,n1,j)-81*f(:,:,n1+1,j)-99*f(:,:,n1+2,j)+73*f(:,:,n1+3,j))
!
      case ('top')               ! top boundary
        f(:,:,n2+1,j)=0.25*(  9*f(:,:,n2,j)- 3*f(:,:,n2-1,j)- 5*f(:,:,n2-2,j)+ 3*f(:,:,n2-3,j))
        f(:,:,n2+2,j)=0.05*( 81*f(:,:,n2,j)-43*f(:,:,n2-1,j)-57*f(:,:,n2-2,j)+39*f(:,:,n2-3,j))
        f(:,:,n2+3,j)=0.05*(127*f(:,:,n2,j)-81*f(:,:,n2-1,j)-99*f(:,:,n2-2,j)+73*f(:,:,n2-3,j))
!
      case default
        print*, "bc_extrap_2_1: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_extrap_2_1
!***********************************************************************
    subroutine bcx_extrap_2_1(f,topbot,j)
!
!  Extrapolation boundary condition for x.
!  Correct for polynomials up to 2nd order, determined 1 further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  19-jun-03/wolf: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        f(l1-1,:,:,j)=0.25*(  9*f(l1,:,:,j)- 3*f(l1+1,:,:,j)- 5*f(l1+2,:,:,j)+ 3*f(l1+3,:,:,j))
        f(l1-2,:,:,j)=0.05*( 81*f(l1,:,:,j)-43*f(l1+1,:,:,j)-57*f(l1+2,:,:,j)+39*f(l1+3,:,:,j))
        f(l1-3,:,:,j)=0.05*(127*f(l1,:,:,j)-81*f(l1+1,:,:,j)-99*f(l1+2,:,:,j)+73*f(l1+3,:,:,j))
!
      case ('top')               ! top boundary
        f(l2+1,:,:,j)=0.25*(  9*f(l2,:,:,j)- 3*f(l2-1,:,:,j)- 5*f(l2-2,:,:,j)+ 3*f(l2-3,:,:,j))
        f(l2+2,:,:,j)=0.05*( 81*f(l2,:,:,j)-43*f(l2-1,:,:,j)-57*f(l2-2,:,:,j)+39*f(l2-3,:,:,j))
        f(l2+3,:,:,j)=0.05*(127*f(l2,:,:,j)-81*f(l2-1,:,:,j)-99*f(l2-2,:,:,j)+73*f(l2-3,:,:,j))
!
      case default
        print*, "bcx_extrap_2_1: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bcx_extrap_2_1
!***********************************************************************
    subroutine bcy_extrap_2_1(f,topbot,j)
!
!  Extrapolation boundary condition for y.
!  Correct for polynomials up to 2nd order, determined 1 further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  19-jun-03/wolf: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        f(:,m1-1,:,j)=0.25*(  9*f(:,m1,:,j)- 3*f(:,m1+1,:,j)- 5*f(:,m1+2,:,j)+ 3*f(:,m1+3,:,j))
        f(:,m1-2,:,j)=0.05*( 81*f(:,m1,:,j)-43*f(:,m1+1,:,j)-57*f(:,m1+2,:,j)+39*f(:,m1+3,:,j))
        f(:,m1-3,:,j)=0.05*(127*f(:,m1,:,j)-81*f(:,m1+1,:,j)-99*f(:,m1+2,:,j)+73*f(:,m1+3,:,j))
!
      case ('top')               ! top boundary
        f(:,m2+1,:,j)=0.25*(  9*f(:,m2,:,j)- 3*f(:,m2-1,:,j)- 5*f(:,m2-2,:,j)+ 3*f(:,m2-3,:,j))
        f(:,m2+2,:,j)=0.05*( 81*f(:,m2,:,j)-43*f(:,m2-1,:,j)-57*f(:,m2-2,:,j)+39*f(:,m2-3,:,j))
        f(:,m2+3,:,j)=0.05*(127*f(:,m2,:,j)-81*f(:,m2-1,:,j)-99*f(:,m2-2,:,j)+73*f(:,m2-3,:,j))
!
      case default
        print*, "bcy_extrap_2_1: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bcy_extrap_2_1
!***********************************************************************
    subroutine bc_extrap_2_2(f,topbot,j)
!
!  Extrapolation boundary condition.
!  Correct for polynomials up to 2nd order, determined 2 further degrees
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  19-jun-03/wolf: coded
!  01-jul-03/axel: introduced abbreviations n1p4,n2m4
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j,n1p4,n2m4
!
!  abbreviations, because otherwise the ifc compiler complains
!  for 1-D runs without vertical extent
!
      n1p4=n1+4
      n2m4=n2-4
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        f(:,:,n1-1,j)=0.2   *(  9*f(:,:,n1,j)                 -  4*f(:,:,n1+2,j)- 3*f(:,:,n1+3,j)+ 3*f(:,:,n1p4,j))
        f(:,:,n1-2,j)=0.2   *( 15*f(:,:,n1,j)- 2*f(:,:,n1+1,j)-  9*f(:,:,n1+2,j)- 6*f(:,:,n1+3,j)+ 7*f(:,:,n1p4,j))
        f(:,:,n1-3,j)=1./35.*(157*f(:,:,n1,j)-33*f(:,:,n1+1,j)-108*f(:,:,n1+2,j)-68*f(:,:,n1+3,j)+87*f(:,:,n1p4,j))
!
      case ('top')               ! top boundary
        f(:,:,n2+1,j)=0.2   *(  9*f(:,:,n2,j)                 -  4*f(:,:,n2-2,j)- 3*f(:,:,n2-3,j)+ 3*f(:,:,n2m4,j))
        f(:,:,n2+2,j)=0.2   *( 15*f(:,:,n2,j)- 2*f(:,:,n2-1,j)-  9*f(:,:,n2-2,j)- 6*f(:,:,n2-3,j)+ 7*f(:,:,n2m4,j))
        f(:,:,n2+3,j)=1./35.*(157*f(:,:,n2,j)-33*f(:,:,n2-1,j)-108*f(:,:,n2-2,j)-68*f(:,:,n2-3,j)+87*f(:,:,n2m4,j))
!
      case default
        print*, "bc_extrap_2_2: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_extrap_2_2
!***********************************************************************
    subroutine bcx_extrap_2_2(f,topbot,j)
!
!  Extrapolation boundary condition.
!  Correct for polynomials up to 2nd order, determined 2 further degrees
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  19-jun-03/wolf: coded
!  01-jul-03/axel: introduced abbreviations n1p4,n2m4
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j,l1p4,l2m4
!
!  abbreviations, because otherwise the ifc compiler complains
!  for 1-D runs without vertical extent
!
      l1p4=l1+4
      l2m4=l2-4
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        f(l1-1,:,:,j)=0.2   *(  9*f(l1,:,:,j)                 -  4*f(l1+2,:,:,j)- 3*f(l1+3,:,:,j)+ 3*f(l1p4,:,:,j))
        f(l1-2,:,:,j)=0.2   *( 15*f(l1,:,:,j)- 2*f(l1+1,:,:,j)-  9*f(l1+2,:,:,j)- 6*f(l1+3,:,:,j)+ 7*f(l1p4,:,:,j))
        f(l1-3,:,:,j)=1./35.*(157*f(l1,:,:,j)-33*f(l1+1,:,:,j)-108*f(l1+2,:,:,j)-68*f(l1+3,:,:,j)+87*f(l1p4,:,:,j))
!
      case ('top')               ! top boundary
        f(l2+1,:,:,j)=0.2   *(  9*f(l2,:,:,j)                 -  4*f(l2-2,:,:,j)- 3*f(l2-3,:,:,j)+ 3*f(l2m4,:,:,j))
        f(l2+2,:,:,j)=0.2   *( 15*f(l2,:,:,j)- 2*f(l2-1,:,:,j)-  9*f(l2-2,:,:,j)- 6*f(l2-3,:,:,j)+ 7*f(l2m4,:,:,j))
        f(l2+3,:,:,j)=1./35.*(157*f(l2,:,:,j)-33*f(l2-1,:,:,j)-108*f(l2-2,:,:,j)-68*f(l2-3,:,:,j)+87*f(l2m4,:,:,j))
!
      case default
        print*, "bcx_extrap_2_2: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bcx_extrap_2_2
!***********************************************************************
    subroutine bcy_extrap_2_2(f,topbot,j)
!
!  Extrapolation boundary condition.
!  Correct for polynomials up to 2nd order, determined 2 further degrees
!  of freedom by minimizing L2 norm of coefficient vector.
!
!   19-jun-03/wolf: coded
!   01-jul-03/axel: introduced abbreviations n1p4,n2m4
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j,m1p4,m2m4
!
!  abbreviations, because otherwise the ifc compiler complains
!  for 1-D runs without vertical extent
!
      m1p4=m1+4
      m2m4=m2-4
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        f(:,m1-1,:,j)=0.2   *(  9*f(:,m1,:,j)                 -  4*f(:,m1+2,:,j)- 3*f(:,m1+3,:,j)+ 3*f(:,m1p4,:,j))
        f(:,m1-2,:,j)=0.2   *( 15*f(:,m1,:,j)- 2*f(:,m1+1,:,j)-  9*f(:,m1+2,:,j)- 6*f(:,m1+3,:,j)+ 7*f(:,m1p4,:,j))
        f(:,m1-3,:,j)=1./35.*(157*f(:,m1,:,j)-33*f(:,m1+1,:,j)-108*f(:,m1+2,:,j)-68*f(:,m1+3,:,j)+87*f(:,m1p4,:,j))
!
      case ('top')               ! top boundary
        f(:,m2+1,:,j)=0.2   *(  9*f(:,m2,:,j)                 -  4*f(:,m2-2,:,j)- 3*f(:,m2-3,:,j)+ 3*f(:,m2m4,:,j))
        f(:,m2+2,:,j)=0.2   *( 15*f(:,m2,:,j)- 2*f(:,m2-1,:,j)-  9*f(:,m2-2,:,j)- 6*f(:,m2-3,:,j)+ 7*f(:,m2m4,:,j))
        f(:,m2+3,:,j)=1./35.*(157*f(:,m2,:,j)-33*f(:,m2-1,:,j)-108*f(:,m2-2,:,j)-68*f(:,m2-3,:,j)+87*f(:,m2m4,:,j))
!
      case default
        print*, "bcy_extrap_2_2: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bcy_extrap_2_2
!***********************************************************************
    subroutine bcy_extrap_2_3(f,topbot,j)
!
!  Extrapolation boundary condition in logarithm:
!  It maintains a power law
!
!  18-dec-08/wlad: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j,l,i
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost
          do n=1,mz
            do l=1,mx
              if (f(l,m1+i,n,j)/=0.) then
                f(l,m1-i,n,j)=f(l,m1,n,j)**2/f(l,m1+i,n,j)
              else
                f(l,m1-i,n,j)=0.
              endif
            enddo
          enddo
        enddo
!
      case ('top')               ! top boundary
        do i=1,nghost
          do n=1,mz
            do l=1,mx
              if (f(l,m2-i,n,j)/=0.) then
                f(l,m2+i,n,j)=f(l,m2,n,j)**2/f(l,m2-i,n,j)
              else
                f(l,m2+i,n,j)=0.
              endif
            enddo
          enddo
        enddo
!
      case default
        print*, "bcy_extrap_2_3: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bcy_extrap_2_3
!***********************************************************************
    subroutine bc_extrap0_2_0(f,topbot,j)
!
!  Extrapolation boundary condition for f(bdry)=0.
!  Correct for polynomials up to 2nd order, determined no further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  09-oct-03/wolf: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      select case (topbot)
!
! Nyquist-filtering
!
      case ('bot')
! bottom boundary
        f(:,:,n1  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n1-1,j)=(1/11.)*&
             (-17*f(:,:,n1+1,j)- 9*f(:,:,n1+2,j)+ 8*f(:,:,n1+3,j))
        f(:,:,n1-2,j)=      2*&
             (- 2*f(:,:,n1+1,j)-   f(:,:,n1+2,j)+   f(:,:,n1+3,j))
        f(:,:,n1-3,j)=(3/11.)*&
             (-27*f(:,:,n1+1,j)-13*f(:,:,n1+2,j)+14*f(:,:,n1+3,j))
!
      case ('top')
! top boundary
        f(:,:,n2  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n2+1,j)=(1/11.)*&
             (-17*f(:,:,n2-1,j)- 9*f(:,:,n2-2,j)+ 8*f(:,:,n2-3,j))
        f(:,:,n2+2,j)=      2*&
             (- 2*f(:,:,n2-1,j)-   f(:,:,n2-2,j)+   f(:,:,n2-3,j))
        f(:,:,n2+3,j)=(3/11.)*&
             (-27*f(:,:,n2-1,j)-13*f(:,:,n2-2,j)+14*f(:,:,n2-3,j))
!
      case default
        print*, "bc_extrap0_2_0: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_extrap0_2_0
!***********************************************************************
    subroutine bc_extrap0_2_1(f,topbot,j)
!
!  Extrapolation boundary condition for f(bdry)=0.
!  Correct for polynomials up to 2nd order, determined 1 further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  NOTE: This is not the final formula, but just bc_extrap_2_1() with f(bdry)=0
!
!  09-oct-03/wolf: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        f(:,:,n1  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n1-1,j)=0.25*(- 3*f(:,:,n1+1,j)- 5*f(:,:,n1+2,j)+ 3*f(:,:,n1+3,j))
        f(:,:,n1-2,j)=0.05*(-43*f(:,:,n1+1,j)-57*f(:,:,n1+2,j)+39*f(:,:,n1+3,j))
        f(:,:,n1-3,j)=0.05*(-81*f(:,:,n1+1,j)-99*f(:,:,n1+2,j)+73*f(:,:,n1+3,j))
!
      case ('top')               ! top boundary
        f(:,:,n2  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n2+1,j)=0.25*(- 3*f(:,:,n2-1,j)- 5*f(:,:,n2-2,j)+ 3*f(:,:,n2-3,j))
        f(:,:,n2+2,j)=0.05*(-43*f(:,:,n2-1,j)-57*f(:,:,n2-2,j)+39*f(:,:,n2-3,j))
        f(:,:,n2+3,j)=0.05*(-81*f(:,:,n2-1,j)-99*f(:,:,n2-2,j)+73*f(:,:,n2-3,j))
!
      case default
        print*, "bc_extrap0_2_1: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_extrap0_2_1
!***********************************************************************
    subroutine bc_extrap0_2_2(f,topbot,j)
!
!  Extrapolation boundary condition for f(bdry)=0.
!  Correct for polynomials up to 2nd order, determined 1 further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  NOTE: This is not the final formula, but just bc_extrap_2_2() with f(bdry)=0
!
!  09-oct-03/wolf: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j,n1p4,n2m4
!
!  abbreviations, because otherwise the ifc compiler complains
!  for 1-D runs without vertical extent
!
      n1p4=n1+4
      n2m4=n2-4
!
      select case (topbot)
!
      case ('bot')
! bottom boundary
        f(:,:,n1  ,j)= 0.       ! set bdry value=0 (indep of initcond)
        f(:,:,n1-1,j)=0.2   *(                 -  4*f(:,:,n1+2,j)- 3*f(:,:,n1+3,j)+ 3*f(:,:,n1p4,j))
        f(:,:,n1-2,j)=0.2   *(- 2*f(:,:,n1+1,j)-  9*f(:,:,n1+2,j)- 6*f(:,:,n1+3,j)+ 7*f(:,:,n1p4,j))
        f(:,:,n1-3,j)=1./35.*(-33*f(:,:,n1+1,j)-108*f(:,:,n1+2,j)-68*f(:,:,n1+3,j)+87*f(:,:,n1p4,j))
!
      case ('top')
! top boundary
        f(:,:,n2  ,j)= 0.       ! set bdry value=0 (indep of initcond)
        f(:,:,n2+1,j)=0.2   *(                 -  4*f(:,:,n2-2,j)- 3*f(:,:,n2-3,j)+ 3*f(:,:,n2m4,j))
        f(:,:,n2+2,j)=0.2   *(- 2*f(:,:,n2-1,j)-  9*f(:,:,n2-2,j)- 6*f(:,:,n2-3,j)+ 7*f(:,:,n2m4,j))
        f(:,:,n2+3,j)=1./35.*(-33*f(:,:,n2-1,j)-108*f(:,:,n2-2,j)-68*f(:,:,n2-3,j)+87*f(:,:,n2m4,j))
!
      case default
        print*, "bc_extrap0_2_2: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_extrap0_2_2
!***********************************************************************
    subroutine bcx_extrap_2_3(f,topbot,j)
!
!  Extrapolation boundary condition in logarithm:
!  It maintains a power law
!
!  18-dec-08/wlad: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j,i
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost
          do n=1,mz
            do m=1,my
              if (f(l1+i,m,n,j)/=0.) then
                f(l1-i,m,n,j)=f(l1,m,n,j)**2/f(l1+i,m,n,j)
              else
                f(l1-i,m,n,j)=0.
              endif
            enddo
          enddo
        enddo
!
      case ('top')               ! top boundary
        do i=1,nghost
          do n=1,mz
            do m=1,my
              if (f(l2-i,m,n,j)/=0.) then
                f(l2+i,m,n,j)=f(l2,m,n,j)**2/f(l2-i,m,n,j)
              else
                f(l2+i,m,n,j)=0.
              endif
            enddo
          enddo
        enddo
!
      case default
        print*, "bcx_extrap_2_3: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bcx_extrap_2_3
!***********************************************************************
    subroutine bcz_extrapol(f,topbot,j)
!
!  Simple linear extrapolation in first order.
!  The last two grid points are used to determine the slope.
!
!  23-nov-10/Bourdin.KIS: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i, j
!
      real, dimension (mx,my) :: m
!
!
      select case (topbot)
      case ('bot')
        ! bottom (left end of the domain)
        m = (f(:,:,n1+1,j) - f(:,:,n1,j)) / (z(n1+1) - z(n1))
        do i = 1, nghost
          f(:,:,n1-i,j) = f(:,:,n1,j) + m * (z(n1-i) - z(n1))
        enddo
      case ('top')
        ! top (right end of the domain)
        m = (f(:,:,n2,j) - f(:,:,n2-1,j)) / (z(n2) - z(n2-1))
        do i = 1, nghost
          f(:,:,n2+i,j) = f(:,:,n2,j) + m * (z(n2+i) - z(n2))
        enddo
      case default
        call fatal_error ('bcz_extrapol', 'invalid argument', lfirst_proc_xy)
      endselect
!
    endsubroutine bcz_extrapol
!***********************************************************************
    subroutine bcz_extrapol_fixed(f,topbot,j)
!
!  Simple linear extrapolation in first order
!  with a fixed value in the first ghost cell.
!  The last two grid points are used to determine the slope.
!
!  23-nov-10/Bourdin.KIS: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i, j
!
      real, dimension (mx,my) :: m
!
!
      select case (topbot)
      case ('bot')
        ! bottom (left end of the domain)
        m = (f(:,:,n1+1,j) - f(:,:,n1,j)) / (z(n1+1) - z(n1))
        do i = 2, nghost
          f(:,:,n1-i,j) = f(:,:,n1-1,j) + m * (z(n1-i) - z(n1-1))
        enddo
      case ('top')
        ! top (right end of the domain)
        m = (f(:,:,n2,j) - f(:,:,n2-1,j)) / (z(n2) - z(n2-1))
        do i = 2, nghost
          f(:,:,n2+i,j) = f(:,:,n2+1,j) + m * (z(n2+i) - z(n2+1))
        enddo
      case default
        call fatal_error ('bcz_extrapol_fixed', 'invalid argument', lfirst_proc_xy)
      endselect
!
    endsubroutine bcz_extrapol_fixed
!***********************************************************************
    subroutine bcz_extrapol_damped(f,topbot,j)
!
!  Simple linear extrapolation in first order
!  with an included damping to zero (useful for velocities).
!  The last two grid points are used to determine the slope.
!  The parameters 'fbcz_bot' and 'fbcz_top' are used as damping factors,
!  which should have values between 0.0 (no damping) and 1.0 (full damping).
!  A typical value of 0.001 corresponds to a half-value time of ~1000 timesteps.
!  Negative values let the damping be persistent, even if 'lfade_damp' is true.
!
!  23-nov-10/Bourdin.KIS: coded
!
      use SharedVariables, only: get_shared_variable
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i, j
!
      real, dimension (mx,my) :: m
      real :: gamma_bot, gamma_top, tau, fade_fact
      real, pointer :: tdamp, tfade_start
      logical, pointer :: ldamp_fade
!
!
      ! bottom and top damping factors:
      gamma_bot = 1.0 - abs (fbcz_bot(j))
      gamma_top = 1.0 - abs (fbcz_top(j))
!
      call get_shared_variable ('ldamp_fade', ldamp_fade)
      if (ldamp_fade) then
        ! fading of damping is active
        call get_shared_variable ('tdamp', tdamp)
        call get_shared_variable ('tfade_start', tfade_start)
        if (t > tfade_start) then
          if (t < tdamp) then
            ! tau is a normalized t, the transition interval is [-0.5, 0.5]:
            tau = (t-tfade_start) / (tdamp-tfade_start) - 0.5
            fade_fact = 0.5 * (1 - tau * (3 - 4*tau**2))
            ! apply damping with fading:
            gamma_bot = 1.0 - abs (fbcz_bot(j)) * fade_fact
            gamma_top = 1.0 - abs (fbcz_top(j)) * fade_fact
          else
            ! damping has already ended (t >= tdamp)
            gamma_bot = 1.0
            gamma_top = 1.0
          endif
          ! apply fading-persistent damping:
          if (fbcz_bot(j) < 0.0) gamma_bot = 1.0 - abs (fbcz_bot(j))
          if (fbcz_top(j) < 0.0) gamma_top = 1.0 - abs (fbcz_top(j))
        endif
      endif
!
      select case (topbot)
      case ('bot')
        ! bottom (left end of the domain)
        m = (f(:,:,n1+1,j) - f(:,:,n1,j)) / (z(n1+1) - z(n1))
        do i = 1, nghost
          f(:,:,n1-i,j) = (f(:,:,n1,j) + m * (z(n1-i) - z(n1))) * gamma_bot**i
        enddo
        m = (f(:,:,n1+1,j) - f(:,:,n1-1,j)) / (z(n1+1) - z(n1-1))
        f(:,:,n1,j) = f(:,:,n1+1,j) + m * (z(n1) - z(n1+1))
      case ('top')
        ! top (right end of the domain)
        m = (f(:,:,n2,j) - f(:,:,n2-1,j)) / (z(n2) - z(n2-1))
        do i = 1, nghost
          f(:,:,n2+i,j) = (f(:,:,n2,j) + m * (z(n2+i) - z(n2))) * gamma_top**i
        enddo
        m = (f(:,:,n2+1,j) - f(:,:,n2-1,j)) / (z(n2+1) - z(n2-1))
        f(:,:,n2,j) = f(:,:,n2-1,j) + m * (z(n2) - z(n2-1))
      case default
        call fatal_error ('bcz_extrapol_damped', 'invalid argument', lfirst_proc_xy)
      endselect
!
    endsubroutine bcz_extrapol_damped
!***********************************************************************
    subroutine bc_db_z(f,topbot,j)
!
!  "One-sided" boundary condition for density.
!  Set ghost zone to reproduce one-sided boundary condition
!  (2nd order):
!  Finding the derivatives on the boundary using a one
!  sided final difference method. This derivative is being
!  used to calculate the boundary points. This will probably
!  only be used for ln(rho)
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the density module
!  13-aug-2002/nils: moved into boundcond
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      real, dimension (:,:), allocatable :: fder
      integer :: i, stat
!
!  Allocate memory for large array.
!
      allocate(fder(mx,my),stat=stat)
      if (stat>0) call fatal_error('bc_db_z', &
          'Could not allocate memory for fder')
!
      select case (topbot)
!
! Bottom boundary
!
      case ('bot')
        do i=1,nghost
          fder=(-3*f(:,:,n1-i+1,j)+4*f(:,:,n1-i+2,j)&
               -f(:,:,n1-i+3,j))/(2*dz)
          f(:,:,n1-i,j)=f(:,:,n1-i+2,j)-2*dz*fder
        enddo
      case ('top')
        do i=1,nghost
          fder=(3*f(:,:,n2+i-1,j)-4*f(:,:,n2+i-2,j)&
               +f(:,:,n2+i-3,j))/(2*dz)
          f(:,:,n2+i,j)=f(:,:,n2+i-2,j)+2*dz*fder
        enddo
      case default
        print*,"bc_db_z: invalid argument for 'bc_db_z'"
      endselect
!
!  Deallocate array.
!
      if (allocated(fder)) deallocate(fder)
!
    endsubroutine bc_db_z
!***********************************************************************
    subroutine bc_db_x(f,topbot,j)
!
!  "One-sided" boundary condition for density.
!  Set ghost zone to reproduce one-sided boundary condition
!  (2nd order):
!  Finding the derivatives on the boundary using a one
!  sided final difference method. This derivative is being
!  used to calculate the boundary points. This will probably
!  only be used for ln(rho)
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the density module
!  13-aug-2002/nils: moved into boundcond
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      real, dimension (:,:), allocatable :: fder
      integer :: i,stat
!
!  Allocate memory for large array.
!
      allocate(fder(my,mz),stat=stat)
      if (stat>0) call fatal_error('bc_db_x', &
          'Could not allocate memory for fder')
!
      select case (topbot)
!
! Bottom boundary
!
      case ('bot')
        do i=1,nghost
          fder=(-3*f(l1-i+1,:,:,j)+4*f(l1-i+2,:,:,j)&
               -f(l1-i+3,:,:,j))/(2*dx)
          f(l1-i,:,:,j)=f(l1-i+2,:,:,j)-2*dx*fder
        enddo
      case ('top')
        do i=1,nghost
          fder=(3*f(l2+i-1,:,:,j)-4*f(l2+i-2,:,:,j)&
               +f(l2+i-3,:,:,j))/(2*dx)
          f(l2+i,:,:,j)=f(l2+i-2,:,:,j)+2*dx*fder
        enddo
      case default
        print*,"bc_db_x: invalid argument for 'bc_db_x'"
      endselect
!
!  Deallocate array.
!
      if (allocated(fder)) deallocate(fder)
!
!
    endsubroutine bc_db_x
!***********************************************************************
    subroutine bc_force_z(f,sgn,topbot,j)
!
!  Force values of j-th variable on vertical boundary topbot.
!  This can either be used for freezing variables at the boundary, or for
!  enforcing a certain time-dependent function of (x,y).
!
!  Currently this is hard-coded for velocity components (ux,uy) and quite
!  useless. Plan is to read time-dependent velocity field from disc and
!  apply it as boundary condition here.
!
!  26-apr-2004/wolf: coded
!
      use EquationOfState, only: gamma_m1, cs2top, cs2bot
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: sgn,i,j
!
      select case (topbot)
!
!  lower boundary
!
      case ('bot')
         select case (force_lower_bound)
         case ('uxy_sin-cos')
            call bc_force_uxy_sin_cos(f,n1,j)
         case ('axy_sin-cos')
            call bc_force_axy_sin_cos(f,n1,j)
         case ('uxy_convection')
            call uu_driver(f)
         !case ('kepler')
         !   call bc_force_kepler(f,n1,j)
         case ('mag_time')
            call bc_force_aa_time(f)
         case ('mag_convection')
            call bc_force_aa_time(f)
            call uu_driver(f)
         case ('cT')
            f(:,:,n1,j) = log(cs2bot/gamma_m1)
         case ('vel_time')
            call bc_force_ux_time(f,n1,j)
         case default
            if (lroot) print*, "No such value for force_lower_bound: <", &
                 trim(force_lower_bound),">"
            call stop_it("")
         endselect
         !
         !  Now fill ghost zones imposing antisymmetry w.r.t. the values just set:
         !
         do i=1,nghost; f(:,:,n1-i,j)=2*f(:,:,n1,j)+sgn*f(:,:,n1+i,j); enddo
!
!  upper boundary
!
      case ('top')
         select case (force_upper_bound)
         case ('uxy_sin-cos')
            call bc_force_uxy_sin_cos(f,n2,j)
         case ('axy_sin-cos')
            call bc_force_axy_sin_cos(f,n2,j)
         case ('uxy_convection')
            call uu_driver(f)
         !case ('kepler')
         !   call bc_force_kepler(f,n2,j)
         case ('cT')
            f(:,:,n2,j) = log(cs2top/gamma_m1)
         case ('vel_time')
            call bc_force_ux_time(f,n2,j)
         case default
            if (lroot) print*, "No such value for force_upper_bound: <", &
                 trim(force_upper_bound),">"
            call stop_it("")
         endselect
         !
         !  Now fill ghost zones imposing antisymmetry w.r.t. the values just set:
         !
         do i=1,nghost; f(:,:,n2+i,j)=2*f(:,:,n2,j)+sgn*f(:,:,n2-i,j); enddo
      case default
        print*,"bc_force_z: invalid argument topbot=",topbot
      endselect
!
    endsubroutine bc_force_z
!***********************************************************************
    subroutine bc_force_x(f, sgn, topbot, j)
!
!  Force values of j-th variable on x-boundaries topbot.
!
!  09-mar-2007/dintrans: coded
!
      use SharedVariables, only : get_shared_variable
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, pointer :: ampl_forc, k_forc, w_forc
      integer :: sgn, i, j, ierr
!
      select case (topbot)
!
!  lower boundary
!
      case ('bot')
         select case (force_lower_bound)
         case ('vel_time')
           if (j /= iuy) call stop_it("BC_FORCE_X: only valid for uy")
           call get_shared_variable('ampl_forc', ampl_forc, ierr)
           if (ierr/=0) call stop_it("BC_FORCE_X: "//&
               "there was a problem when getting ampl_forc")
           call get_shared_variable('k_forc', k_forc, ierr)
           if (ierr/=0) call stop_it("BC_FORCE_X: "//&
               "there was a problem when getting k_forc")
           call get_shared_variable('w_forc', w_forc, ierr)
           if (ierr/=0) call stop_it("BC_FORCE_X: "//&
               "there was a problem when getting w_forc")
           if (headtt) print*, 'BC_FORCE_X: ampl_forc, k_forc, w_forc=',&
               ampl_forc, k_forc, w_forc
           f(l1,:,:,iuy) = spread(ampl_forc*sin(k_forc*y)*cos(w_forc*t), 2, mz)
         case default
            if (lroot) print*, "No such value for force_lower_bound: <", &
                 trim(force_lower_bound),">"
            call stop_it("")
         endselect
         !
         !  Now fill ghost zones imposing antisymmetry w.r.t. the values just set:
         !
         do i=1,nghost; f(l1-i,:,:,j)=2*f(l1,:,:,j)+sgn*f(l1+i,:,:,j); enddo
!
!  upper boundary
!
      case ('top')
         select case (force_upper_bound)
         case ('vel_time')
            if (j /= iuy) call stop_it("BC_FORCE_X: only valid for uy")
            call get_shared_variable('ampl_forc', ampl_forc, ierr)
            if (ierr/=0) call stop_it("BC_FORCE_X: "//&
                   "there was a problem when getting ampl_forc")
            call get_shared_variable('k_forc', k_forc, ierr)
            if (ierr/=0) call stop_it("BC_FORCE_X: "//&
                   "there was a problem when getting k_forc")
            call get_shared_variable('w_forc', w_forc, ierr)
            if (ierr/=0) call stop_it("BC_FORCE_X: "//&
                   "there was a problem when getting w_forc")
            if (headtt) print*, 'BC_FORCE_X: ampl_forc, k_forc, w_forc=',&
                   ampl_forc, k_forc, w_forc
            f(l2,:,:,iuy) = spread(ampl_forc*sin(k_forc*y)*cos(w_forc*t), 2, mz)
         case default
            if (lroot) print*, "No such value for force_upper_bound: <", &
                 trim(force_upper_bound),">"
            call stop_it("")
         endselect
         !
         !  Now fill ghost zones imposing antisymmetry w.r.t. the values just set:
         !
         do i=1,nghost; f(l2+i,:,:,j)=2*f(l2,:,:,j)+sgn*f(l2-i,:,:,j); enddo
      case default
        print*,"bc_force_x: invalid argument topbot=",topbot
      endselect
!
    endsubroutine bc_force_x
!***********************************************************************
    subroutine bc_force_uxy_sin_cos(f,idz,j)
!
!  Set (ux, uy) = (cos y, sin x) in vertical layer
!
!  26-apr-2004/wolf: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: idz,j
      real :: kx,ky
!
      if (iuz == 0) call stop_it("BC_FORCE_UXY_SIN_COS: Bad idea...")
!
      if (j==iux) then
        if (Ly>0) then; ky=2*pi/Ly; else; ky=0.; endif
        f(:,:,idz,j) = spread(cos(ky*y),1,mx)
      elseif (j==iuy) then
        if (Lx>0) then; kx=2*pi/Lx; else; kx=0.; endif
        f(:,:,idz,j) = spread(sin(kx*x),2,my)
      elseif (j==iuz) then
        f(:,:,idz,j) = 0.
      endif
!
    endsubroutine bc_force_uxy_sin_cos
!***********************************************************************
    subroutine bc_force_axy_sin_cos(f,idz,j)
!
!  Set (ax, ay) = (cos y, sin x) in vertical layer
!
!  26-apr-2004/wolf: coded
!  10-apr-2005/axel: adapted for A
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: idz,j
      real :: kx,ky
!
      if (iaz == 0) call stop_it("BC_FORCE_AXY_SIN_COS: Bad idea...")
!
      if (j==iax) then
        if (Ly>0) then; ky=2*pi/Ly; else; ky=0.; endif
        f(:,:,idz,j) = spread(cos(ky*y),1,mx)
      elseif (j==iay) then
        if (Lx>0) then; kx=2*pi/Lx; else; kx=0.; endif
        f(:,:,idz,j) = spread(sin(kx*x),2,my)
      elseif (j==iaz) then
        f(:,:,idz,j) = 0.
      endif
!
    endsubroutine bc_force_axy_sin_cos
!!***********************************************************************
    subroutine bc_one_x(f,topbot,j)
!
!  Set bdry values to 1 for debugging purposes
!
!  11-jul-02/wolf: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
      character (len=3) :: topbot
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
          f(1:l1-1,:,:,j)=1.
!
      case ('top')               ! top boundary
          f(l2+1:mx,:,:,j)=1.
!
      case default
        print*, "bc_one_x: ",topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_one_x
!***********************************************************************
    subroutine bc_one_y(f,topbot,j)
!
!  Set bdry values to 1 for debugging purposes
!
!  11-jul-02/wolf: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
      character (len=3) :: topbot
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
          f(:,1:m1-1,:,j)=1.
!
      case ('top')               ! top boundary
          f(:,m2+1:my,:,j)=1.
!
      case default
        print*, "bc_one_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_one_y
!***********************************************************************
    subroutine bc_one_z(f,topbot,j)
!
!  Set bdry values to 1 for debugging purposes
!
!  11-jul-02/wolf: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
      character (len=3) :: topbot
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
          f(:,:,1:n1-1,j)=1.
!
      case ('top')               ! top boundary
          f(:,:,n2+1:mz,j)=1.
!
      case default
        print*, "bc_one_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_one_z
!***********************************************************************
    subroutine bc_freeze_var_x(topbot,j)
!
!  Tell other modules that variable with slot j is to be frozen in on
!  given boundary
!
      integer :: j
      character (len=3) :: topbot
!
      lfrozen_bcs_x = .true.    ! set flag
!
      select case (topbot)
      case ('bot')               ! bottom boundary
        lfrozen_bot_var_x(j) = .true.
      case ('top')               ! top boundary
        lfrozen_top_var_x(j) = .true.
      case default
        print*, "bc_freeze_var_x: ", topbot, " should be 'top' or 'bot'"
      endselect
!
    endsubroutine bc_freeze_var_x
!***********************************************************************
    subroutine bc_freeze_var_y(topbot,j)
!
!  Tell other modules that variable with slot j is to be frozen in on
!  given boundary
!
      integer :: j
      character (len=3) :: topbot
!
      lfrozen_bcs_y = .true.    ! set flag
!
      select case (topbot)
      case ('bot')               ! bottom boundary
        lfrozen_bot_var_y(j) = .true.
      case ('top')               ! top boundary
        lfrozen_top_var_y(j) = .true.
      case default
        print*, "bc_freeze_var_y: ", topbot, " should be 'top' or 'bot'"
      endselect
!
    endsubroutine bc_freeze_var_y
!***********************************************************************
    subroutine bc_freeze_var_z(topbot,j)
!
!  Tell other modules that variable with slot j is to be frozen in on
!  given boundary
!
      integer :: j
      character (len=3) :: topbot
!
      lfrozen_bcs_z = .true.    ! set flag
!
      select case (topbot)
      case ('bot')               ! bottom boundary
        lfrozen_bot_var_z(j) = .true.
      case ('top')               ! top boundary
        lfrozen_top_var_z(j) = .true.
      case default
        print*, "bc_freeze_var_z: ", topbot, " should be 'top' or 'bot'"
      endselect
!
    endsubroutine bc_freeze_var_z
!***********************************************************************
     subroutine uu_driver(f,quenching)
!
!  Simulated velocity field used as photospherec motions
!  Use of velocity field produced by Boris Gudiksen
!
!  27-mai-04/bing: coded
!  11-aug-06/axel: make it compile with nprocx>0, renamed quenching -> quen
!  18-jun-08/bing: quenching depends on B^2, not only Bz^2
!
       use EquationOfState, only : gamma,gamma_m1,gamma_inv,cs20,lnrho0
       use Mpicomm, only : mpisend_real, mpirecv_real
       use Syscalls, only : file_exists
!
       real, dimension (mx,my,mz,mfarray) :: f
!
       real, dimension (:,:), save, allocatable :: uxl,uxr,uyl,uyr
       real, dimension (:,:), allocatable :: uxd,uyd,quen,pp,betaq,fac
       real, dimension (:,:), allocatable :: bbx,bby,bbz,bb2,tmp
       integer :: tag_xl=321,tag_yl=322,tag_xr=323,tag_yr=324
       integer :: tag_tl=345,tag_tr=346,tag_dt=347
       integer :: lend=0,ierr,frame=0,stat,pos,iref,px,py
       real, save :: tl=0.,tr=0.,delta_t=0.
       real  :: zmin
       logical, optional :: quenching
       logical :: quench
!
       character (len=*), parameter :: vel_times_dat = 'driver/vel_times.dat'
       character (len=*), parameter :: vel_field_dat = 'driver/vel_field.dat'
       integer :: unit=1
!
       intent (inout) :: f
!
       if (lroot .and. .not. file_exists(vel_times_dat)) &
           call stop_it_if_any(.true.,'uu_driver: Could not find file "'//trim(vel_times_dat)//'"')
       if (lroot .and. .not. file_exists(vel_field_dat)) &
           call stop_it_if_any(.true.,'uu_driver: Could not find file "'//trim(vel_field_dat)//'"')
!
       ierr = 0
       stat = 0
       if (.not.allocated(uxl))  allocate(uxl(nx,ny),stat=ierr)
       if (.not.allocated(uxr))  allocate(uxr(nx,ny),stat=stat)
       ierr = max(stat,ierr)
       if (.not.allocated(uyl))  allocate(uyl(nx,ny),stat=stat)
       ierr = max(stat,ierr)
       if (.not.allocated(uyr))  allocate(uyr(nx,ny),stat=stat)
       ierr = max(stat,ierr)
       allocate(uxd(nx,ny),stat=stat);         ierr = max(stat,ierr)
       allocate(uyd(nx,ny),stat=stat);         ierr = max(stat,ierr)
       allocate(quen(nx,ny),stat=stat);        ierr = max(stat,ierr)
       allocate(pp(nx,ny),stat=stat);          ierr = max(stat,ierr)
       allocate(betaq(nx,ny),stat=stat);       ierr = max(stat,ierr)
       allocate(fac(nx,ny),stat=stat);         ierr = max(stat,ierr)
       allocate(bbx(nx,ny),stat=stat);         ierr = max(stat,ierr)
       allocate(bby(nx,ny),stat=stat);         ierr = max(stat,ierr)
       allocate(bbz(nx,ny),stat=stat);         ierr = max(stat,ierr)
       allocate(bb2(nx,ny),stat=stat);         ierr = max(stat,ierr)
       allocate(tmp(nxgrid,nygrid),stat=stat); ierr = max(stat,ierr)
!
       if (ierr>0) call stop_it_if_any(.true.,'uu_driver: '// &
           'Could not allocate memory for all variable, please check')
!
       if (present(quenching)) then
         quench = quenching
       else
         ! Right now quenching is per default active
         quench=.true.
       endif
!
!  Read the time table
!
       if ((t*unit_time<tl+delta_t) .or. (t*unit_time>=tr+delta_t)) then
!
         if (lroot) then
           inquire(IOLENGTH=lend) tl
           open (unit,file=vel_times_dat,form='unformatted',status='unknown',recl=lend,access='direct')
!
           ierr = 0
           frame = 0
           do while (ierr == 0)
             frame=frame+1
             read (unit,rec=frame,iostat=ierr) tl
             read (unit,rec=frame+1,iostat=ierr) tr
             if (ierr /= 0) then
               frame=1
               delta_t = t*unit_time                  ! EOF is reached => read again
               read (unit,rec=frame,iostat=ierr) tl
               read (unit,rec=frame+1,iostat=ierr) tr
               ierr=-1
             else
               if (t*unit_time>=tl+delta_t .and. t*unit_time<tr+delta_t) ierr=-1
               ! correct time step is reached
             endif
           enddo
           close (unit)
!
           do px=0, nprocx-1
             do py=0, nprocy-1
               if ((px /= 0) .or. (py /= 0)) then
                 call mpisend_real (tl, 1, px+py*nprocx, tag_tl)
                 call mpisend_real (tr, 1, px+py*nprocx, tag_tr)
                 call mpisend_real (delta_t, 1, px+py*nprocx, tag_dt)
               endif
             enddo
           enddo
         else
           call mpirecv_real (tl, 1, 0, tag_tl)
           call mpirecv_real (tr, 1, 0, tag_tr)
           call mpirecv_real (delta_t, 1, 0, tag_dt)
         endif
!
! Read velocity field
!
         if (lroot) then
           open (unit,file=vel_field_dat,form='unformatted',status='unknown',recl=lend*nxgrid*nygrid,access='direct')
!
           read (unit,rec=2*frame-1) tmp
           do px=0, nprocx-1
             do py=0, nprocy-1
               if ((px /= 0) .or. (py /= 0)) then
                 uxl = tmp(px*nx+1:(px+1)*nx,py*ny+1:(py+1)*ny)
                 call mpisend_real (uxl, (/ nx, ny /), px+py*nprocx, tag_xl)
               endif
             enddo
           enddo
           uxl = tmp(1:nx,1:ny)
!
           read (unit,rec=2*frame) tmp
           do px=0, nprocx-1
             do py=0, nprocy-1
               if ((px /= 0) .or. (py /= 0)) then
                 uyl = tmp(px*nx+1:(px+1)*nx,py*ny+1:(py+1)*ny)
                 call mpisend_real (uyl, (/ nx, ny /), px+py*nprocx, tag_yl)
               endif
             enddo
           enddo
           uyl = tmp(1:nx,1:ny)
!
           read (unit,rec=2*frame+1) tmp
           do px=0, nprocx-1
             do py=0, nprocy-1
               if ((px /= 0) .or. (py /= 0)) then
                 uxr = tmp(px*nx+1:(px+1)*nx,py*ny+1:(py+1)*ny)
                 call mpisend_real (tmp(px*nx+1:(px+1)*nx,py*ny+1:(py+1)*ny), (/ nx, ny /), px+py*nprocx, tag_xr)
               endif
             enddo
           enddo
           uxr = tmp(1:nx,1:ny)
!
           read (unit,rec=2*frame+2) tmp
           uyr = tmp(1:nx,1:ny)
           do px=0, nprocx-1
             do py=0, nprocy-1
               if ((px /= 0) .or. (py /= 0)) then
                 uyr = tmp(px*nx+1:(px+1)*nx,py*ny+1:(py+1)*ny)
                 call mpisend_real (tmp(px*nx+1:(px+1)*nx,py*ny+1:(py+1)*ny), (/ nx, ny /), px+py*nprocx, tag_yr)
               endif
             enddo
           enddo
           uyr = tmp(1:nx,1:ny)
!
           close (unit)
         else
           call mpirecv_real (uxl, (/ nx, ny /), 0, tag_xl)
           call mpirecv_real (uyl, (/ nx, ny /), 0, tag_yl)
           call mpirecv_real (uxr, (/ nx, ny /), 0, tag_xr)
           call mpirecv_real (uyr, (/ nx, ny /), 0, tag_yr)
         endif
!
         uxl = uxl / 10. / unit_velocity
         uxr = uxr / 10. / unit_velocity
         uyl = uyl / 10. / unit_velocity
         uyr = uyr / 10. / unit_velocity
!
       endif
!
!   simple linear interploation between timesteps
!
       if (tr /= tl) then
         uxd  = (t*unit_time - (tl+delta_t)) * (uxr - uxl) / (tr - tl) + uxl
         uyd  = (t*unit_time - (tl+delta_t)) * (uyr - uyl) / (tr - tl) + uyl
       else
         uxd = uxl
         uyd = uyl
       endif
!
!   suppress footpoint motion at low plasma beta
!
       zmin = minval(abs(z(n1:n2)))
       iref = n1
       do pos=n1,n2
         if (abs(z(pos))==zmin) iref=pos; exit
       enddo
!
!   Calculate B^2 for plasma beta
!
       if (quench) then
!-----------------------------------------------------------------------
         if (nygrid/=1) then
           fac=(1./60)*spread(dy_1(m1:m2),1,nx)
           bbx= fac*(+ 45.0*(f(l1:l2,m1+1:m2+1,iref,iaz)-f(l1:l2,m1-1:m2-1,iref,iaz)) &
               -  9.0*(f(l1:l2,m1+2:m2+2,iref,iaz)-f(l1:l2,m1-2:m2-2,iref,iaz)) &
               +      (f(l1:l2,m1+3:m2+3,iref,iaz)-f(l1:l2,m1-3:m2-3,iref,iaz)))
         else
           if (ip<=5) print*, 'uu_driver: Degenerate case in y-direction'
         endif
         if (nzgrid/=1) then
           fac=(1./60)*spread(spread(dz_1(iref),1,nx),2,ny)
           bbx= bbx -fac*(+ 45.0*(f(l1:l2,m1:m2,iref+1,iay)-f(l1:l2,m1:m2,iref-1,iay)) &
               -  9.0*(f(l1:l2,m1:m2,iref+2,iay)-f(l1:l2,m1:m2,iref-2,iay)) &
               +      (f(l1:l2,m1:m2,iref+3,iay)-f(l1:l2,m1:m2,iref-2,iay)))
         else
           if (ip<=5) print*, 'uu_driver: Degenerate case in z-direction'
         endif
!-----------------------------------------------------------------------
         if (nzgrid/=1) then
           fac=(1./60)*spread(spread(dz_1(iref),1,nx),2,ny)
           bby= fac*(+ 45.0*(f(l1:l2,m1:m2,iref+1,iax)-f(l1:l2,m1:m2,iref-1,iax)) &
               -  9.0*(f(l1:l2,m1:m2,iref+2,iax)-f(l1:l2,m1:m2,iref-2,iax)) &
               +      (f(l1:l2,m1:m2,iref+3,iax)-f(l1:l2,m1:m2,iref-3,iax)))
         else
           if (ip<=5) print*, 'uu_driver: Degenerate case in z-direction'
         endif
         if (nxgrid/=1) then
           fac=(1./60)*spread(dx_1(l1:l2),2,ny)
           bby=bby-fac*(+45.0*(f(l1+1:l2+1,m1:m2,iref,iaz)-f(l1-1:l2-1,m1:m2,iref,iaz)) &
               -  9.0*(f(l1+2:l2+2,m1:m2,iref,iaz)-f(l1-2:l2-2,m1:m2,iref,iaz)) &
               +      (f(l1+3:l2+3,m1:m2,iref,iaz)-f(l1-3:l2-3,m1:m2,iref,iaz)))
         else
           if (ip<=5) print*, 'uu_driver: Degenerate case in x-direction'
         endif
!-----------------------------------------------------------------------
         if (nxgrid/=1) then
           fac=(1./60)*spread(dx_1(l1:l2),2,ny)
           bbz= fac*(+ 45.0*(f(l1+1:l2+1,m1:m2,iref,iay)-f(l1-1:l2-1,m1:m2,iref,iay)) &
               -  9.0*(f(l1+2:l2+2,m1:m2,iref,iay)-f(l1-2:l2-2,m1:m2,iref,iay)) &
               +      (f(l1+3:l2+3,m1:m2,iref,iay)-f(l1-3:l2-3,m1:m2,iref,iay)))
         else
           if (ip<=5) print*, 'uu_driver: Degenerate case in x-direction'
         endif
         if (nygrid/=1) then
           fac=(1./60)*spread(dy_1(m1:m2),1,nx)
           bbz=bbz-fac*(+45.0*(f(l1:l2,m1+1:m2+1,iref,iax)-f(l1:l2,m1-1:m2-1,iref,iax)) &
               -  9.0*(f(l1:l2,m1+2:m2+2,iref,iax)-f(l1:l2,m1-2:m2-2,iref,iax)) &
               +      (f(l1:l2,m1+3:m2+3,iref,iax)-f(l1:l2,m1-3:m2-3,iref,iax)))
         else
           if (ip<=5) print*, 'uu_driver: Degenerate case in y-direction'
         endif
!-----------------------------------------------------------------------
!
         bb2 = bbx*bbx + bby*bby + bbz*bbz
         bb2 = bb2/(2.*mu0)
!
         if (ltemperature) then
           pp=gamma_m1*gamma_inv*exp(f(l1:l2,m1:m2,iref,ilnrho)+f(l1:l2,m1:m2,iref,ilnTT))
         else if (lentropy) then
           if (pretend_lnTT) then
             pp=gamma_m1*gamma_inv*exp(f(l1:l2,m1:m2,iref,ilnrho)+f(l1:l2,m1:m2,iref,iss))
           else
             pp=gamma*(f(l1:l2,m1:m2,iref,iss)+ &
                 f(l1:l2,m1:m2,iref,ilnrho))-gamma_m1*lnrho0
             pp=exp(pp) * cs20*gamma_inv
           endif
         else
           pp=gamma_inv*cs20*exp(lnrho0)
         endif
!
!   limit plasma beta
!
         betaq = pp / max(tini,bb2)*1e-3
!
         quen=(1.+betaq**2)/(10.+betaq**2)
       else
         quen(:,:)=1.
       endif
!
!   Fill z=0 layer with velocity field
!
       f(l1:l2,m1:m2,iref,iux)=uxd*quen
       f(l1:l2,m1:m2,iref,iuy)=uyd*quen
       if (iref/=n1) f(l1:l2,m1:m2,n1,iux:iuz)=0.
!
       if (allocated(uxd)) deallocate(uxd)
       if (allocated(uyd)) deallocate(uyd)
       if (allocated(quen)) deallocate(quen)
       if (allocated(pp)) deallocate(pp)
       if (allocated(betaq)) deallocate(betaq)
       if (allocated(fac)) deallocate(fac)
       if (allocated(bbx)) deallocate(bbx)
       if (allocated(bby)) deallocate(bby)
       if (allocated(bbz)) deallocate(bbz)
       if (allocated(bb2)) deallocate(bb2)
       if (allocated(tmp)) deallocate(tmp)
!
     endsubroutine uu_driver
!***********************************************************************
    subroutine bc_force_aa_time(f)
!
!  Reads in time series of magnetograms
!
!  17-feb-10/bing: coded
!  25-jul-10/Bourdin.KIS: parallelized
!
      use Fourier, only : field_extrapol_z_parallel
      use Mpicomm, only : mpisend_real, mpirecv_real, &
                          mpisend_logical, mpirecv_logical
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, save :: t_l=0., t_r=0., delta_t=0.
      integer :: ierr, lend, frame, stat, rec_l, rec_r
      integer :: rec_vxl, rec_vxr, rec_vyl, rec_vyr ! l- and r-record position if file
      integer, parameter :: bnx=nxgrid, bny=ny/nprocx ! data in pencil shape
      integer, parameter :: enx=nygrid, eny=nx/nprocy ! transposed data in pencil shape
      integer :: px, py, partner
      integer, parameter :: tag_l=208, tag_r=209, tag_dt=210
      logical, save :: luse_vel_field = .false., first_run = .true.
      logical :: ex
!
      ! temporal storage for frames before (l) and after (r) current time step:
      real, dimension (:,:), allocatable, save :: Bz0_l, Bz0_r
      real, dimension (:,:), allocatable, save :: vx_l, vx_r, vy_l, vy_r
      real, dimension (:,:), allocatable :: vx_tmp, vy_tmp
      ! current magnetic field z-component (interpolated between l and r):
      real, dimension (:,:), allocatable :: Bz0
      ! current velocity x- and y-component (interpolated between l and r):
      real, dimension (:,:), allocatable, save :: vx, vy
!
      real, dimension (:,:,:), allocatable, save :: exp_fact ! exponential factor
      real, dimension (:,:), allocatable :: k_2 ! wave vector length
      integer :: kx_start, pos_z, i
      real :: delta_z
      real, parameter :: reduce_factor=0.25, enhance_factor=1.0
!
      real :: mu0_SI, u_b, time_SI
!
      character (len=*), parameter :: mag_field_dat = 'driver/mag_field.dat'
      character (len=*), parameter :: mag_times_dat = 'driver/mag_times.dat'
      character (len=*), parameter :: mag_vel_field_dat = 'driver/mag_vel_field.dat'
!
      if (first_run) then
!
!  Check for consistency:
        if ((.not. lequidist(1)) .or. (.not. lequidist(2))) &
            call fatal_error ('bc_force_aa_time', 'not yet implemented for non-equidistant grids', lfirst_proc_xy)
        if (mod (nygrid, nprocxy) /= 0) &
            call fatal_error ('bc_force_aa_time', 'nygrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
        if (mod (nxgrid, nprocxy) /= 0) &
            call fatal_error ('bc_force_aa_time', 'nxgrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
!
!  Check for existence of necessary driver files:
        if (lfirst_proc_xy) then
          inquire (file=mag_field_dat, exist=ex)
          if (.not. ex) call fatal_error ('bc_force_aa_time', 'File does not exists: '//trim(mag_field_dat), .true.)
          inquire (file=mag_times_dat, exist=ex)
          if (.not. ex) call fatal_error ('bc_force_aa_time', 'File does not exists: '//trim(mag_times_dat), .true.)
          inquire (file=mag_vel_field_dat, exist=ex)
          if (ex) then
            luse_vel_field = .true.
            print *, 'bc_force_aa_time: using time dependant magnetogram _with_ corresponding horizontal velocities.'
          else
            print *, 'bc_force_aa_time: using time dependant magnetogram _without_ corresponding horizontal velocities.'
          endif
        endif
!
        do i=1,nprocxy-1
          if (lroot) then
            call mpisend_logical(luse_vel_field,1,i,i)
          elseif (iproc==i) then
            call mpirecv_logical(luse_vel_field,1,0,iproc)
          endif
        enddo
!
!  Initialization of magnetograms and velocity fields.
        allocate(Bz0_l(bnx,bny),stat=stat)
        if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for Bz0_l',.true.)
        allocate(Bz0_r(bnx,bny),stat=stat)
        if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for Bz0_r',.true.)
        if (luse_vel_field) then
          allocate(vx_l(nx,ny),stat=stat)
          if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for vx_l',.true.)
          allocate(vx_r(nx,ny),stat=stat)
          if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for vx_r',.true.)
          allocate(vy_l(nx,ny),stat=stat)
          if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for vy_l',.true.)
          allocate(vy_r(nx,ny),stat=stat)
          if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for vy_r',.true.)
          allocate(vx(nx,ny),stat=stat)
          if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for vx',.true.)
          allocate(vy(nx,ny),stat=stat)
          if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for vy',.true.)
        endif
!
        first_run = .false.
!
      endif
!
      allocate(Bz0(bnx,bny),stat=stat)
      if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for Bz0',.true.)
!
      time_SI = t*unit_time
!
      if (t_r+delta_t <= time_SI) then
!
        if (lfirst_proc_xy) then
!  Read and distribute Bz data (in pencil shape)
!
          inquire (IOLENGTH=lend) t_l
          open (10,file=mag_times_dat,form='unformatted',status='unknown', &
              recl=lend,access='direct')
!
          ierr = 0
          t_l = 0.
          frame = 0
          do while (ierr == 0)
            frame=frame+1
            read (10,rec=frame,iostat=ierr) t_l
            read (10,rec=frame+1,iostat=ierr) t_r
            if (ierr /= 0) then
              ! EOF is reached => read again
              delta_t = time_SI
              frame=1
              read (10,rec=frame,iostat=ierr) t_l
              read (10,rec=frame+1,iostat=ierr) t_r
              ierr=-1
            else
              ! test, if correct time step is reached
              if ((t_l+delta_t < time_SI) .and. (t_r+delta_t > time_SI)) ierr = -1
            endif
          enddo
          close (10)
!
          if (luse_vel_field) then
            allocate (vx_tmp(nxgrid,nygrid), stat=stat)
            if (stat>0) call fatal_error ('bc_force_aa_time', 'Could not allocate memory for vx_tmp', .true.)
            allocate (vy_tmp(nxgrid,nygrid), stat=stat)
            if (stat>0) call fatal_error ('bc_force_aa_time', 'Could not allocate memory for vy_tmp', .true.)
            open (10, file=mag_vel_field_dat, form='unformatted', status='unknown', &
                recl=lend*nxgrid*nygrid, access='direct')
!
            rec_vxl = 1 + (frame-1)*nprocxy*2
            rec_vxr = 1 + frame*nprocxy*2
            rec_vyl = rec_vxl + frame*nprocxy
            rec_vyr = rec_vxr + frame*nprocxy
!
            ! read _l data in the order of occurence in file
            read (10,rec=rec_vxl) vx_l
            read (10,rec=rec_vyl) vy_l
!
            ! send _l data to remote
            do py = 1, nprocy
              do px = 1, nprocx
                partner = px + py*nprocx + ipz*nprocxy
                if (partner == iproc) cycle
                vx_l = vx_tmp(1+(px-1)*nprocx:px*nprocx,1+(py-1)*nprocy:py*nprocy)
                vy_l = vy_tmp(1+(px-1)*nprocx:px*nprocx,1+(py-1)*nprocy:py*nprocy)
                call mpisend_real (vx_l, (/ nx, ny /), partner, tag_l)
                call mpisend_real (vy_l, (/ nx, ny /), partner, tag_r)
              enddo
            enddo
            ! read local _l data
            vx_l = vx_tmp(1:nprocx,1:nprocy)
            vy_l = vy_tmp(1:nprocx,1:nprocy)
!
            ! read _r data in the order of occurence in file
            read (10,rec=rec_vxr) vx_r
            read (10,rec=rec_vyr) vy_r
!
            ! send _r data to remote
            do py = 1, nprocy
              do px = 1, nprocx
                partner = px + py*nprocx + ipz*nprocxy
                if (partner == iproc) cycle
                vx_r = vx_tmp(1+(px-1)*nprocx:px*nprocx,1+(py-1)*nprocy:py*nprocy)
                vy_r = vy_tmp(1+(px-1)*nprocx:px*nprocx,1+(py-1)*nprocy:py*nprocy)
                call mpisend_real (vx_r, (/ nx, ny /), partner, tag_l)
                call mpisend_real (vy_r, (/ nx, ny /), partner, tag_r)
              enddo
            enddo
            ! read local _r data
            vx_r = vx_tmp(1:nprocx,1:nprocy)
            vy_r = vy_tmp(1:nprocx,1:nprocy)
!
            close (10)
            if (allocated (vx_tmp)) deallocate (vx_tmp)
            if (allocated (vy_tmp)) deallocate (vy_tmp)
          endif
!
          open (10,file=mag_field_dat,form='unformatted',status='unknown', &
              recl=lend*bnx*bny,access='direct')
          rec_l = 1 + (frame-1)*nprocxy
          rec_r = 1 + frame*nprocxy
          do py=1, nprocxy-1
            partner = py + ipz*nprocxy
            ! read Bz data for remote processors
            read (10,rec=rec_l+py) Bz0_l
            read (10,rec=rec_r+py) Bz0_r
            ! send Bz data to remote
            call mpisend_real (Bz0_l, (/ bnx, bny /), partner, tag_l)
            call mpisend_real (Bz0_r, (/ bnx, bny /), partner, tag_r)
            call mpisend_real (t_l, 1, partner, tag_l)
            call mpisend_real (t_r, 1, partner, tag_r)
            call mpisend_real (delta_t, 1, partner, tag_dt)
          enddo
          ! read local Bz data
          read (10,rec=rec_l) Bz0_l
          read (10,rec=rec_r) Bz0_r
          close (10)
!
        else
!
          if (luse_vel_field) then
            ! wait for vx and vy data from root processor
            call mpirecv_real (vx_l, (/ nx, ny /), ipz*nprocxy, tag_l)
            call mpirecv_real (vy_l, (/ nx, ny /), ipz*nprocxy, tag_r)
            call mpirecv_real (vx_r, (/ nx, ny /), ipz*nprocxy, tag_l)
            call mpirecv_real (vy_r, (/ nx, ny /), ipz*nprocxy, tag_r)
          endif
!
          ! wait for Bz data from root processor
          call mpirecv_real (Bz0_l, (/ bnx, bny /), ipz*nprocxy, tag_l)
          call mpirecv_real (Bz0_r, (/ bnx, bny /), ipz*nprocxy, tag_r)
          call mpirecv_real (t_l, 1, ipz*nprocxy, tag_l)
          call mpirecv_real (t_r, 1, ipz*nprocxy, tag_r)
          call mpirecv_real (delta_t, 1, ipz*nprocxy, tag_dt)
!
        endif
!
        ! convert units
        mu0_SI = 4.*pi*1.e-7
        u_b = unit_velocity*sqrt(mu0_SI/mu0*unit_density)
!
        Bz0_l = Bz0_l * 1e-4 / u_b
        Bz0_r = Bz0_r * 1e-4 / u_b
!
        if (luse_vel_field) then
          vx_l = vx_l / unit_velocity
          vy_l = vy_l / unit_velocity
          vx_r = vx_r / unit_velocity
          vy_r = vy_r / unit_velocity
        endif
!
      endif
!
      Bz0 = (time_SI - (t_l+delta_t)) * (Bz0_r - Bz0_l) / (t_r - t_l) + Bz0_l
!
      if (luse_vel_field) then
        vx = (time_SI - (t_l+delta_t)) * (vx_r - vx_l) / (t_r - t_l) + vx_l
        vy = (time_SI - (t_l+delta_t)) * (vy_r - vy_l) / (t_r - t_l) + vy_l
      endif
!
!  Fourier Transform of Bz0:
!
      if (.not. allocated (exp_fact)) then
!
        ! Setup exponential factor for bottom boundary
        allocate (exp_fact(enx,eny,nghost+1), stat=stat)
        if (stat > 0) call fatal_error ('bc_force_aa_time', 'Could not allocate memory for exp_fact', .true.)
        allocate (k_2(enx,eny), stat=stat)
        if (stat > 0) call fatal_error ('bc_force_aa_time', 'Could not allocate memory for k_2', .true.)
!
        ! Get wave numbers already in transposed pencil shape and calculate exp(|k|)
        kx_start = (ipx+ipy*nprocx)*eny
        k_2 = spread (ky_fft(1:enx), 2, eny)**2 + spread (kx_fft(kx_start+1:kx_start+eny), 1, enx)**2
        if (kx_start == 0) k_2(1,1) = 1.0 ! dummy value to avoid division by zero
        exp_fact = spread (exp (sqrt (k_2)), 3, nghost+1)
!
        ! Setup increase of fourrier coefficients for bottom boundary
        do pos_z = 1, nghost
          delta_z = z(n1) - z(n1-nghost+pos_z-1) ! dz is positive => increase
          ! Enhance (delta_z<0) or reduce (delta_z>0) extrapolation
          if (delta_z < 0.0) delta_z = delta_z * enhance_factor
          if (delta_z > 0.0) delta_z = delta_z * reduce_factor
          ! Include normalization factor for fourier transform: 1/(nxgrid*nygrid)
          exp_fact(:,:,pos_z) = exp_fact(:,:,pos_z) ** delta_z / (k_2 * nxgrid*nygrid)
        enddo
        exp_fact(:,:,nghost+1) = 1.0 / (k_2 * nxgrid*nygrid)
        if (kx_start == 0) exp_fact(1,1,:) = 1.0 / (nxgrid*nygrid)
!
        if (allocated (k_2)) deallocate (k_2)
      endif
!
      call field_extrapol_z_parallel (Bz0, f(l1:l2,m1:m2,n1-nghost:n1,iax:iay), exp_fact)
      call communicate_vect_field_ghosts (f, 'bot')
      if (luse_vel_field) call communicate_vect_field_ghosts (f, 'bot', iux)
!
      if (allocated(Bz0)) deallocate(Bz0)
!
    endsubroutine bc_force_aa_time
!***********************************************************************
    subroutine bc_lnTT_flux_x(f,topbot)
!
!  Constant flux boundary condition for temperature (called when bcx='c1')
!
!  12-Mar-2007/dintrans: coded
!
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=3) :: topbot
!
      real, pointer :: hcond0, hcond1, Fbot
      real, dimension (:,:), allocatable :: tmp_yz
      integer :: i,ierr,stat
!
!  Allocate memory for large array.
!
      allocate(tmp_yz(my,mz),stat=stat)
      if (stat>0) call fatal_error('bc_lnTT_flux_x', &
          'Could not allocate memory for tmp_yz')
!
!  Do the 'c1' boundary condition (constant heat flux) for lnTT.
!  check whether we want to do top or bottom (this is processor dependent)
!
      call get_shared_variable('hcond0',hcond0,ierr)
      if (ierr/=0) call stop_it("bc_lnTT_flux_x: "//&
           "there was a problem when getting hcond0")
      call get_shared_variable('hcond1',hcond1,ierr)
      if (ierr/=0) call stop_it("bc_lnTT_flux_x: "//&
           "there was a problem when getting hcond1")
      call get_shared_variable('Fbot',Fbot,ierr)
      if (ierr/=0) call stop_it("bc_lnTT_flux_x: "//&
           "there was a problem when getting Fbot")
!
      if (headtt) print*,'bc_lnTT_flux_x: Fbot,hcond,dx=',Fbot,hcond0*hcond1,dx
!
      select case (topbot)
!
!  bottom boundary
!  ===============
!
      case ('bot')
        tmp_yz=-Fbot/(hcond0*hcond1)/exp(f(l1,:,:,ilnTT))
!
!  enforce dlnT/dx = - Fbot/(K*T)
!
        do i=1,nghost
          f(l1-i,:,:,ilnTT)=f(l1+i,:,:,ilnTT)-2*i*dx*tmp_yz
        enddo
!
      case default
        call fatal_error('bc_lnTT_flux_x','invalid argument')
!
      endselect
!
!  Deallocate large array.
!
      if (allocated(tmp_yz)) deallocate(tmp_yz)
!
    endsubroutine bc_lnTT_flux_x
!***********************************************************************
    subroutine bc_lnTT_flux_z(f,topbot)
!
!  Constant flux boundary condition for temperature
!  (called when bcz='c1')
!
!  12-May-07/dintrans: coded
!
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=3) :: topbot
!
      real, dimension (:,:), allocatable :: tmp_xy
      real, pointer :: hcond0, Fbot
      integer :: i,ierr,stat
!
!  Allocate memory for large array.
!
      allocate(tmp_xy(mx,my),stat=stat)
      if (stat>0) call fatal_error('bc_lnTT_flux_x', &
          'Could not allocate memory for tmp_xy')
!
!  Do the 'c1' boundary condition (constant heat flux) for lnTT or TT (if
!  ltemperature_nolog=.true.) at the bottom _only_.
!  lnTT version: enforce dlnT/dz = - Fbot/(K*T)
!    TT version: enforce   dT/dz = - Fbot/K
!
      call get_shared_variable('hcond0',hcond0,ierr)
      if (ierr/=0) call stop_it("bc_lnTT_flux_z: "//&
           "there was a problem when getting hcond0")
      call get_shared_variable('Fbot',Fbot,ierr)
      if (ierr/=0) call stop_it("bc_lnTT_flux_z: "//&
           "there was a problem when getting Fbot")
!
      if (headtt) print*,'bc_lnTT_flux_z: Fbot,hcond,dz=',Fbot,hcond0,dz
!
      select case (topbot)
      case ('bot')
        if (ltemperature_nolog) then
          tmp_xy=-Fbot/hcond0
        else
          tmp_xy=-Fbot/hcond0/exp(f(:,:,n1,ilnTT))
        endif
        do i=1,nghost
          f(:,:,n1-i,ilnTT)=f(:,:,n1+i,ilnTT)-2.*i*dz*tmp_xy
        enddo
!
      case default
        call fatal_error('bc_lnTT_flux_z','invalid argument')
!
      endselect
!
!  Deallocate large array.
!
      if (allocated(tmp_xy)) deallocate(tmp_xy)
!
    endsubroutine bc_lnTT_flux_z
!***********************************************************************
    subroutine bc_ss_flux_x(f,topbot)
!
!  Constant flux boundary condition for entropy (called when bcx='c1')
!
!  17-mar-07/dintrans: coded
!
      use EquationOfState, only: gamma, gamma_m1, lnrho0, cs20
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=3) :: topbot
!
      real, dimension (:,:), allocatable :: tmp_yz,cs2_yz
      real, pointer :: FbotKbot, FtopKtop
      integer :: i,ierr,stat
!
!  Allocate memory for large arrays.
!
      allocate(tmp_yz(my,mz),stat=stat)
      if (stat>0) call fatal_error('bc_ss_flux_x', &
          'Could not allocate memory for tmp_yz')
      allocate(cs2_yz(my,mz),stat=stat)
      if (stat>0) call fatal_error('bc_ss_flux_x', &
          'Could not allocate memory for cs2_yz')
!
!  Do the 'c1' boundary condition (constant heat flux) for entropy and
!  pretend_lnTT. Check whether we want to do top or bottom (this is
!  processor dependent)
!
      call get_shared_variable('FbotKbot',FbotKbot,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux_x: "//&
           "there was a problem when getting FbotKbot")
!
      select case (topbot)
!
!  bottom boundary
!  ===============
!
      case ('bot')
        if (headtt) print*,'bc_ss_flux_x: FbotKbot=',FbotKbot
!
!  Deal with the simpler pretend_lnTT=T case first. Now ss is actually
!  lnTT and the boundary condition reads glnTT=FbotKbot/T
!
        if (pretend_lnTT) then
          do i=1,nghost
            f(l1-i,:,:,iss)=f(l1+i,:,:,iss)+2*i*dx*FbotKbot/exp(f(l1,:,:,iss))
          enddo
        else
!
!  calculate Fbot/(K*cs2)
!
!       cs2_yz=cs20*exp(gamma_m1*(f(l1,:,:,ilnrho)-lnrho0)+cv1*f(l1,:,:,iss))
!
!  Both, bottom and top boundary conditions are corrected for linear density
!
          if (ldensity_nolog) then
            cs2_yz=cs20*exp(gamma_m1*(log(f(l1,:,:,irho))-lnrho0)+gamma*f(l1,:,:,iss))
          else
            cs2_yz=cs20*exp(gamma_m1*(f(l1,:,:,ilnrho)-lnrho0)+gamma*f(l1,:,:,iss))
          endif
          tmp_yz=FbotKbot/cs2_yz
!
!  enforce ds/dx + gamma_m1/gamma*dlnrho/dx = - gamma_m1/gamma*Fbot/(K*cs2)
!
          do i=1,nghost
!           f(l1-i,:,:,iss)=f(l1+i,:,:,iss)+(cp-cv)* &
            if (ldensity_nolog) then
              f(l1-i,:,:,iss)=f(l1+i,:,:,iss)+gamma_m1/gamma* &
                  (log(f(l1+i,:,:,irho))-log(f(l1-i,:,:,irho))+2*i*dx*tmp_yz)
            else
              f(l1-i,:,:,iss)=f(l1+i,:,:,iss)+gamma_m1/gamma* &
                  (f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho)+2*i*dx*tmp_yz)
            endif
          enddo
        endif
!
!  top boundary
!  ============
!
      case ('top')
!
        call get_shared_variable('FtopKtop',FtopKtop,ierr)
        if (ierr/=0) call stop_it("bc_ss_flux_x: "//&
             "there was a problem when getting FtopKtop")
!
        if (headtt) print*,'bc_ss_flux_x: FtopKtop=',FtopKtop
!
!  Deal with the simpler pretend_lnTT=T case first. Now ss is actually
!  lnTT and the boundary condition reads glnTT=FtopKtop/T
!
        if (pretend_lnTT) then
          do i=1,nghost
            f(l2+i,:,:,iss)=f(l2-i,:,:,iss)-2*i*dx*FtopKtop/exp(f(l2,:,:,iss))
          enddo
        else
!
!  calculate Ftop/(K*cs2)
!
          if (ldensity_nolog) then
            cs2_yz=cs20*exp(gamma_m1*(log(f(l2,:,:,irho))-lnrho0)+gamma*f(l2,:,:,iss))
          else
            cs2_yz=cs20*exp(gamma_m1*(f(l2,:,:,ilnrho)-lnrho0)+gamma*f(l2,:,:,iss))
          endif
          tmp_yz=FtopKtop/cs2_yz
!
!  enforce ds/dx + gamma_m1/gamma*dlnrho/dx = - gamma_m1/gamma*Ftop/(K*cs2)
!
          do i=1,nghost
            if (ldensity_nolog) then
              f(l2+i,:,:,iss)=f(l2-i,:,:,iss)+gamma_m1/gamma* &
                  (log(f(l2-i,:,:,irho))-log(f(l2+i,:,:,irho))-2*i*dx*tmp_yz)
            else
              f(l2+i,:,:,iss)=f(l2-i,:,:,iss)+gamma_m1/gamma* &
                  (f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho)-2*i*dx*tmp_yz)
            endif
          enddo
!
        endif
!
        case default
          call fatal_error('bc_ss_flux_x','invalid argument')
!
        endselect
!
!  Deallocate large arrays.
!
      if (allocated(tmp_yz)) deallocate(tmp_yz)
      if (allocated(cs2_yz)) deallocate(cs2_yz)
!
    endsubroutine bc_ss_flux_x
!***********************************************************************
    subroutine bc_del2zero(f,topbot,j)
!
!  Potential field boundary condition
!
!  11-oct-06/wolf: Adapted from Tobi's bc_aa_pot2
!
      use Fourier, only: fourier_transform_xy_xy
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
      integer, intent (in) :: j
!
      real, dimension (:,:), allocatable :: kx,ky,kappa,exp_fact,tmp_re,tmp_im
      integer :: i,stat
!
!  Allocate memory for large arrays.
!
      allocate(kx(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_del2zero', &
          'Could not allocate memory for kx')
      allocate(ky(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_del2zero', &
          'Could not allocate memory for ky')
      allocate(kappa(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_del2zero', &
          'Could not allocate memory for kappa')
      allocate(exp_fact(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_del2zero', &
          'Could not allocate memory for exp_fact')
      allocate(tmp_re(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_del2zero', &
          'Could not allocate memory for tmp_im')
!
!  Get local wave numbers
!
      kx = spread(kx_fft(ipx*nx+1:ipx*nx+nx),2,ny)
      ky = spread(ky_fft(ipy*ny+1:ipy*ny+ny),1,nx)
!
!  Calculate 1/k^2, zero mean
!
      if (lshear) then
        kappa = sqrt((kx+deltay*ky/Lx)**2+ky**2)
      else
        kappa = sqrt(kx**2 + ky**2)
      endif
!
!  Check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  Potential field condition at the bottom
!
      case ('bot')
!
        do i=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          exp_fact = exp(-kappa*(z(n1+i)-z(n1-i)))
!
!  Determine potential field in ghost zones
!
          !  Fourier transforms of x- and y-components on the boundary
          tmp_re = f(l1:l2,m1:m2,n1+i,j)
          tmp_im = 0.0
          call fourier_transform_xy_xy(tmp_re,tmp_im)
          tmp_re = tmp_re*exp_fact
          tmp_im = tmp_im*exp_fact
          ! Transform back
          call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
          f(l1:l2,m1:m2,n1-i,j) = tmp_re
!
        enddo
!
!  Potential field condition at the top
!
      case ('top')
!
        do i=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          exp_fact = exp(-kappa*(z(n2+i)-z(n2-i)))
!
!  Determine potential field in ghost zones
!
          !  Fourier transforms of x- and y-components on the boundary
          tmp_re = f(l1:l2,m1:m2,n2-i,j)
          tmp_im = 0.0
          call fourier_transform_xy_xy(tmp_re,tmp_im)
          tmp_re = tmp_re*exp_fact
          tmp_im = tmp_im*exp_fact
          ! Transform back
          call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
          f(l1:l2,m1:m2,n2+i,j) = tmp_re
!
        enddo
!
      case default
!
        if (lroot) print*,"bc_del2zero: invalid argument"
!
      endselect
!
!  Deallocate large arrays.
!
      if (allocated(kx)) deallocate(kx)
      if (allocated(ky)) deallocate(ky)
      if (allocated(kappa)) deallocate(kappa)
      if (allocated(exp_fact)) deallocate(exp_fact)
      if (allocated(tmp_re)) deallocate(tmp_re)
      if (allocated(tmp_im)) deallocate(tmp_im)
!
    endsubroutine bc_del2zero
!***********************************************************************
    subroutine bc_zero_x(f,topbot,j)
!
!  Zero value in the ghost zones.
!
!  11-aug-2009/anders: implemented
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=3) :: topbot
      integer :: j
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        f(1:l1-1,:,:,j)=0.0
!
!  Top boundary.
!
      case ('top')
        f(l2+1:mx,:,:,j)=0.0
!
!  Default.
!
      case default
        print*, "bc_zero_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_zero_x
!***********************************************************************
    subroutine bc_zero_z(f,topbot,j)
!
!  Zero value in the ghost zones.
!
!  13-aug-2007/anders: implemented
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=3) :: topbot
      integer :: j
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        f(:,:,1:n1-1,j)=0.0
!
!  Top boundary.
!
      case ('top')
        f(:,:,n2+1:mz,j)=0.0
!
!  Default.
!
      case default
        print*, "bc_zero_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_zero_z
!***********************************************************************
    subroutine bc_inflow_z(f,topbot,j,lforce_ghost)
!
!  Inflow boundary conditions.
!
!  If the velocity vector points out of the box, the velocity boundary
!  condition is set to 's', otherwise it is set to 'a'.
!  If 'lforce_ghost' is true, the boundary and ghost cell values are forced
!  to not point outwards. Otherwise the boundary value is forced to be 0.
!
!  25-dec-2010/Bourdin.KIS: adapted from 'bc_outflow_z'
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
      logical, optional :: lforce_ghost
!
      integer :: i, ix, iy
      logical :: lforce
!
      lforce = .false.
      if (present (lforce_ghost)) lforce = lforce_ghost
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do iy=1,my; do ix=1,mx
          if (f(ix,iy,n1,j)>0.0) then  ! 's'
            do i=1,nghost; f(ix,iy,n1-i,j)=+f(ix,iy,n1+i,j); enddo
          else                         ! 'a'
            do i=1,nghost; f(ix,iy,n1-i,j)=-f(ix,iy,n1+i,j); enddo
            f(ix,iy,n1,j)=0.0
          endif
          if (lforce) then
            do i = 0, nghost
              if (f(ix,iy,n1-i,j) < 0.0) f(ix,iy,n1-i,j) = 0.0
            enddo
          endif
        enddo; enddo
!
!  Top boundary.
!
      case ('top')
        do iy=1,my; do ix=1,mx
          if (f(ix,iy,n2,j)<0.0) then  ! 's'
            do i=1,nghost; f(ix,iy,n2+i,j)=+f(ix,iy,n2-i,j); enddo
          else                         ! 'a'
            do i=1,nghost; f(ix,iy,n2+i,j)=-f(ix,iy,n2-i,j); enddo
            f(ix,iy,n2,j)=0.0
          endif
          if (lforce) then
            do i = 0, nghost
              if (f(ix,iy,n2+i,j) > 0.0) f(ix,iy,n2+i,j) = 0.0
            enddo
          endif
        enddo; enddo
!
!  Default.
!
      case default
        print*, "bc_inflow_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_inflow_z
!***********************************************************************
    subroutine bc_outflow_z(f,topbot,j,lforce_ghost)
!
!  Outflow boundary conditions.
!
!  If the velocity vector points out of the box, the velocity boundary
!  condition is set to 's', otherwise it is set to 'a'.
!  If 'lforce_ghost' is true, the boundary and ghost cell values are forced
!  to not point inwards. Otherwise the boundary value is forced to be 0.
!
!  12-aug-2007/anders: implemented
!  25-dec-2010/Bourdin.KIS: added forcing of boundary and ghost cell values
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
      logical, optional :: lforce_ghost
!
      integer :: i, ix, iy
      logical :: lforce
!
      lforce = .false.
      if (present (lforce_ghost)) lforce = lforce_ghost
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do iy=1,my; do ix=1,mx
          if (f(ix,iy,n1,j)<0.0) then  ! 's'
            do i=1,nghost; f(ix,iy,n1-i,j)=+f(ix,iy,n1+i,j); enddo
          else                         ! 'a'
            do i=1,nghost; f(ix,iy,n1-i,j)=-f(ix,iy,n1+i,j); enddo
            f(ix,iy,n1,j)=0.0
          endif
          if (lforce) then
            do i = 0, nghost
              if (f(ix,iy,n1-i,j) > 0.0) f(ix,iy,n1-i,j) = 0.0
            enddo
          endif
        enddo; enddo
!
!  Top boundary.
!
      case ('top')
        do iy=1,my; do ix=1,mx
          if (f(ix,iy,n2,j)>0.0) then  ! 's'
            do i=1,nghost; f(ix,iy,n2+i,j)=+f(ix,iy,n2-i,j); enddo
          else                         ! 'a'
            do i=1,nghost; f(ix,iy,n2+i,j)=-f(ix,iy,n2-i,j); enddo
            f(ix,iy,n2,j)=0.0
          endif
          if (lforce) then
            do i = 0, nghost
              if (f(ix,iy,n2+i,j) < 0.0) f(ix,iy,n2+i,j) = 0.0
            enddo
          endif
        enddo; enddo
!
!  Default.
!
      case default
        print*, "bc_outflow_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_outflow_z
!***********************************************************************
    subroutine bc_inflow_zero_deriv_z(f,topbot,j)
!
!  Inflow boundary condition, tries to create zero 1st derivative at boundary.
!
!  If the velocity vector points out of the box, the velocity boundary
!  condition is set to 's', otherwise it is set to relax to zero derivative.
!  The boundary and ghost cell values are forced to not point outwards.
!
!  27-dec-2010/Bourdin.KIS: adapted from 'bc_outflow_const_deriv_z'
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      integer :: i, ix, iy
!
      select case (topbot)
!
      ! bottom boundary
      case ('bot')
        do iy = 1, my
          do ix = 1, mx
            ! 's' boundary condition
            do i = 1, nghost
              f(ix,iy,n1-i,j) = f(ix,iy,n1+i,j)
            enddo
            ! force inflow
            do i = 0, nghost
              if (f(ix,iy,n1-i,j) < 0.0) f(ix,iy,n1-i,j) = 0.0
            enddo
          enddo
        enddo
!
      ! top boundary
      case ('top')
        do iy = 1, my
          do ix = 1, mx
            ! 's' boundary condition
            do i = 1, nghost
              f(ix,iy,n2+i,j) = f(ix,iy,n2-i,j)
            enddo
            ! force inflow
            do i = 0, nghost
              if (f(ix,iy,n2+i,j) > 0.0) f(ix,iy,n2+i,j) = 0.0
            enddo
          enddo
        enddo
!
      case default
        print*, "bc_inflow_zero_deriv_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_inflow_zero_deriv_z
!***********************************************************************
    subroutine bc_outflow_zero_deriv_z(f,topbot,j)
!
!  Outflow boundary condition, tries to create zero 1st derivative at boundary.
!
!  If the velocity vector points out of the box, the velocity boundary
!  condition is set to 's', otherwise it is set to relax to zero derivative.
!  The boundary and ghost cell values are forced to not point inwards.
!
!  27-dec-2010/Bourdin.KIS: adapted from 'bc_outflow_z'
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      integer :: i, ix, iy
!
      select case (topbot)
!
      ! bottom boundary
      case ('bot')
        do iy = 1, my
          do ix = 1, mx
            ! 's' boundary condition
            do i = 1, nghost
              f(ix,iy,n1-i,j) = f(ix,iy,n1+i,j)
            enddo
            ! force outflow
            do i = 0, nghost
              if (f(ix,iy,n1-i,j) > 0.0) f(ix,iy,n1-i,j) = 0.0
            enddo
          enddo
        enddo
!
      ! top boundary
      case ('top')
        do iy = 1, my
          do ix = 1, mx
            ! 's' boundary condition
            do i = 1, nghost
              f(ix,iy,n2+i,j) = f(ix,iy,n2-i,j)
            enddo
            ! force outflow
            do i = 0, nghost
              if (f(ix,iy,n2+i,j) < 0.0) f(ix,iy,n2+i,j) = 0.0
            enddo
          enddo
        enddo
!
      case default
        print*, "bc_outflow_zero_deriv_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_outflow_zero_deriv_z
!***********************************************************************
    subroutine bc_steady_z(f,topbot,j)
!
!  Steady in/outflow boundary conditions.
!
!  If the velocity vector points out of the box, the velocity boundary
!  condition is set to 's', otherwise it is set to boundary value.
!
!  06-nov-2010/fred: implemented
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      integer :: i, ix, iy
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do iy=1,my; do ix=1,mx
          if (f(ix,iy,n1,j)<0.0) then  !'s'
            do i=1,nghost; f(ix,iy,n1-i,j)=f(ix,iy,n1+i,j); enddo
          else                         !'u(n1)'
            do i=1,nghost; f(ix,iy,n1-i,j)=f(ix,iy,n1,j); enddo
          endif
        enddo; enddo
!
!  Top boundary.
!
      case ('top')
        do iy=1,my; do ix=1,mx
          if (f(ix,iy,n2,j)>0.0) then  !'s'
            do i=1,nghost; f(ix,iy,n2+i,j)=f(ix,iy,n2-i,j); enddo
          else                         !'u(n2)'
            do i=1,nghost; f(ix,iy,n2+i,j)=f(ix,iy,n2,j); enddo
          endif
        enddo; enddo
!
!  Default.
!
      case default
        print*, "bc_steady_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_steady_z
!***********************************************************************
    subroutine bc_copy_x(f,topbot,j)
!
!  Copy value in last grid point to all ghost cells.
!
!  11-aug-2009/anders: implemented
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      integer :: i
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do i=1,nghost; f(l1-i,:,:,j)=f(l1,:,:,j); enddo
!
!  Top boundary.
!
      case ('top')
        do i=1,nghost; f(l2+i,:,:,j)=f(l2,:,:,j); enddo
!
!  Default.
!
      case default
        print*, "bc_copy_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_copy_x
!***********************************************************************
    subroutine bc_copy_y(f,topbot,j)
!
!  Copy value in last grid point to all ghost cells.
!
!  08-june-2010/wlyra: implemented
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      integer :: i
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do i=1,nghost; f(:,m1-i,:,j)=f(:,m1,:,j); enddo
!
!  Top boundary.
!
      case ('top')
        do i=1,nghost; f(:,m2+i,:,j)=f(:,m2,:,j); enddo
!
!  Default.
!
      case default
        print*, "bc_copy_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_copy_y
!***********************************************************************
    subroutine bc_copy_y_noinflow(f,topbot,j)
!
!  Copy value in last grid point to all ghost cells. Set to zero if
!  the sign is wrong. This bc is different from outflow (cop). Outflow
!  just copies the last point to the ghost cells, thus permitting both
!  outflow (uy pointing out of the box) and inflow (uy pointing back to
!  the box). 'c+k' is a no-inflow, purely outflow boundary. It sets the
!  velocity to zero if that was pointing back to the box. The 'k' means
!  "kill". "copy if outflow, kill if inflow".
!
!  08-june-2010/wlyra: implemented
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real :: value
      integer :: j,l,n
!
      integer :: i
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do l=1,mx; do n=1,mz
          value=0.
          if (f(l,m1,n,j)<0) value=f(l,m1,n,j)
          do i=1,nghost
            f(l,m1-i,n,j)=value
          enddo
        enddo;enddo
!
!  Top boundary.
!
      case ('top')
        do l=1,mx; do n=1,mz
          value=0.
          if (f(l,m2,n,j) > 0) value=f(l,m2,n,j)
          do i=1,nghost
            f(l,m2+i,n,j)=value
          enddo
        enddo; enddo
!
!  Default.
!
      case default
        print*, "bc_copy_y_noinflow: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_copy_y_noinflow
!***********************************************************************
    subroutine bc_copy_z(f,topbot,j)
!
!  Copy value in last grid point to all ghost cells.
!
!  15-aug-2007/anders: implemented
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      integer :: i
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do i=1,nghost; f(:,:,n1-i,j)=f(:,:,n1,j); enddo
!
!  Top boundary.
!
      case ('top')
        do i=1,nghost; f(:,:,n2+i,j)=f(:,:,n2,j); enddo
!
!  Default.
!
      case default
        print*, "bc_copy_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_copy_z
!***********************************************************************
    subroutine bc_frozen_in_bb(topbot,j)
!
!  Set flags to indicate that magnetic flux is frozen-in at the
!  z boundary. The implementation occurs in daa_dt where magnetic
!  diffusion is switched off in that layer.
!
      use SharedVariables, only: get_shared_variable
!
      character (len=3) :: topbot
      integer :: j
!
      logical, save :: lfirstcall=.true.
      logical, pointer, save, dimension (:) :: lfrozen_bb_bot, lfrozen_bb_top
      integer :: ierr
!
      if (lfirstcall) then
        call get_shared_variable('lfrozen_bb_bot',lfrozen_bb_bot,ierr)
        if (ierr/=0) call fatal_error('bc_frozen_in_bb', &
            'there was a problem getting lfrozen_bb_bot')
        call get_shared_variable('lfrozen_bb_top',lfrozen_bb_top,ierr)
        if (ierr/=0) call fatal_error('bc_frozen_in_bb', &
            'there was a problem getting lfrozen_bb_top')
      endif
!
      select case (topbot)
      case ('bot')               ! bottom boundary
        lfrozen_bb_bot(j-iax+1) = .true.    ! set flag
      case ('top')               ! top boundary
        lfrozen_bb_top(j-iax+1) = .true.    ! set flag
      case default
        print*, "bc_frozen_in_bb: ", topbot, " should be 'top' or 'bot'"
      endselect
!
      lfirstcall=.false.
!
    endsubroutine bc_frozen_in_bb
!***********************************************************************
    subroutine bcz_hydrostatic_temp(f,topbot)
!
!  The logarithmic density in the ghost cells is used to calculate the
!  logarithmic tempterature under the asumption of a hydrostatic equilibrium.
!
!  19-nov-2010/Bourdin.KIS: coded
!
      use EquationOfState, only: gamma, gamma_m1, get_cp1
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
!
      integer :: i
      real, dimension (mx,my) :: T_inv
      real :: g_ref, delta_z, inv_cp_cv, cp_inv
      real, dimension (:), pointer :: gravz_zpencil
!
!
      call get_shared_variable ('gravz_zpencil', gravz_zpencil)
      call get_cp1 (cp_inv)
      inv_cp_cv = gamma / gamma_m1 * cp_inv
!
      select case (topbot)
      case ('bot')
        ! bottom (left end of the domain)
        do i = 1, nghost
          delta_z = z(n1-i) - z(n1-i+1)
          g_ref = gravz_zpencil(n1-i+1)
          T_inv = exp (-f(:,:,n1-i+1,ilnTT))
          f(:,:,n1-i,ilnTT) = f(:,:,n1-i+1,ilnTT) + f(:,:,n1-i+1,ilnrho) - f(:,:,n1-i,ilnrho) + g_ref*delta_z*inv_cp_cv*T_inv
        enddo
      case ('top')
        ! top (right end of the domain)
        do i = 1, nghost
          delta_z = z(n2+i) - z(n2+i-1)
          g_ref = gravz_zpencil(n2+i-1)
          T_inv = exp (-f(:,:,n2+i-1,ilnTT))
          f(:,:,n2+i,ilnTT) = f(:,:,n2+i-1,ilnTT) + f(:,:,n2+i-1,ilnrho) - f(:,:,n2+i,ilnrho) + g_ref*delta_z*inv_cp_cv*T_inv
        enddo
      case default
        call fatal_error ('bcz_hydrostatic_temp', 'invalid argument', lfirst_proc_xy)
      endselect
!
    endsubroutine bcz_hydrostatic_temp
!***********************************************************************
    subroutine bc_aa_pot_field_extrapol(f,topbot)
!
!  Potential field extrapolation in z-direction for the ghost cells.
!
!  9-jul-2010/Bourdin.KIS: coded
!
      use Fourier, only: vect_pot_extrapol_z_parallel
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
!
      real, dimension (:,:,:), allocatable, save :: exp_fact_top, exp_fact_bot
      integer, parameter :: bnx=nygrid, bny=nx/nprocy
      integer :: kx_start, stat, pos_z
      real :: delta_z, reduce_factor=1.
!
      ! reduce_factor reduces the structure increase at the bottom boundary
      ! to help numerically resolving the strong gradients in the ghost cells.
      ! Set reduce_factor to [0,1] by using fbcz_bot(iaa) in run_pars.
      ! A value of 1 (default) switches this reducing-mechanism off.
      ! A value of 0 just clones the vector field A from the n1-layer.
      if (fbcz_bot(iaa) /= 0.) reduce_factor = fbcz_bot(iaa)
!
      if (.not. ((lfirst_proc_z .and. (topbot == 'bot')) .or. (llast_proc_z .and. (topbot == 'top')))) &
          call fatal_error ('bc_aa_pot_field_extrapol', 'Only implemented for topmost or downmost z-layer.', lfirst_proc_xy)
!
      if (mod (nx, nprocy) /= 0) &
          call fatal_error ('bc_aa_pot_field_extrapol', 'nx needs to be an integer multiple of nprocy.', lfirst_proc_xy)
!
!  Check whether we want to do top or bottom z boundary
!
      select case (topbot)
      case ('bot')
        if (.not. allocated (exp_fact_bot)) then
          ! Setup exponential factor for bottom boundary
          allocate (exp_fact_bot(bnx,bny,nghost), stat=stat)
          if (stat > 0) call fatal_error ('bc_aa_pot_field_extrapol', 'Could not allocate memory for exp_fact_bot', .true.)
          ! Get wave numbers already in transposed pencil shape and calculate exp(|k|)
          kx_start = (ipx+ipy*nprocx)*bny
          exp_fact_bot = spread (exp (sqrt (spread (ky_fft(1:bnx), 2, bny) ** 2 + &
                                            spread (kx_fft(kx_start+1:kx_start+bny), 1, bnx) ** 2)), 3, nghost)
          do pos_z = 1, nghost
            ! dz is positive => enhance structures or contrast
            delta_z = reduce_factor * (z(n1) - z(n1-nghost+pos_z-1))
            ! Include normalization factor for fourier transform: 1/(nxgrid*nygrid)
            exp_fact_bot(:,:,pos_z) = exp_fact_bot(:,:,pos_z) ** delta_z / (nxgrid*nygrid)
          enddo
        endif
        call vect_pot_extrapol_z_parallel &
             (f(l1:l2,m1:m2,n1,iax:iaz), f(l1:l2,m1:m2,n1-nghost:n1-1,iax:iaz), exp_fact_bot)
      case ('top')
        if (.not. allocated (exp_fact_top)) then
          ! Setup exponential factor for top boundary
          allocate (exp_fact_top(bnx,bny,nghost), stat=stat)
          if (stat > 0) call fatal_error ('bc_aa_pot_field_extrapol', 'Could not allocate memory for exp_fact_top', .true.)
          ! Get wave numbers already in transposed pencil shape and calculate exp(|k|)
          kx_start = (ipx+ipy*nprocx)*bny
          exp_fact_top = spread (exp (sqrt (spread (ky_fft(1:bnx), 2, bny) ** 2 + &
                                            spread (kx_fft(kx_start+1:kx_start+bny), 1, bnx) ** 2)), 3, nghost)
          do pos_z = 1, nghost
            ! dz is negative => decay of structures
            delta_z = z(n2) - z(n2+pos_z)
            ! Include normalization factor for fourier transform: 1/(nxgrid*nygrid)
            exp_fact_top(:,:,pos_z) = exp_fact_top(:,:,pos_z) ** delta_z / (nxgrid*nygrid)
          enddo
        endif
        call vect_pot_extrapol_z_parallel &
             (f(l1:l2,m1:m2,n2,iax:iaz), f(l1:l2,m1:m2,n2+1:n2+nghost,iax:iaz), exp_fact_top)
      case default
        call fatal_error ('bc_aa_pot_field_extrapol', 'invalid argument', lfirst_proc_xy)
      endselect
!
!  The vector potential needs to be known outside of (l1:l2,m1:m2) as well
!
      call communicate_vect_field_ghosts(f,topbot)
!
    endsubroutine bc_aa_pot_field_extrapol
!***********************************************************************
    subroutine bc_aa_pot3(f,topbot)
!
!  Potential field boundary condition
!
!  11-oct-06/wolf: Adapted from Tobi's bc_aa_pot2
!
      use Fourier, only: fourier_transform_xy_xy
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
!
      real, dimension (:,:,:), allocatable :: aa_re,aa_im
      real, dimension (:,:), allocatable :: kx,ky,kappa,exp_fact
      real, dimension (:,:), allocatable :: tmp_re,tmp_im
      real    :: delta_z
      integer :: i,j,stat
!
!  Allocate memory for large arrays.
!
      allocate(aa_re(nx,ny,iax:iaz),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot3', &
          'Could not allocate memory for aa_re',.true.)
      allocate(aa_im(nx,ny,iax:iaz),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot3', &
          'Could not allocate memory for aa_im',.true.)
      allocate(kx(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot3', &
          'Could not allocate memory for kx',.true.)
      allocate(ky(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot3', &
          'Could not allocate memory for ky',.true.)
      allocate(kappa(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot3', &
          'Could not allocate memory for kappa',.true.)
      allocate(exp_fact(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot3', &
          'Could not allocate memory for exp_fact',.true.)
      allocate(tmp_re(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot3', &
          'Could not allocate memory for tmp_re',.true.)
      allocate(tmp_im(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot3', &
          'Could not allocate memory for tmp_im',.true.)
!
!  Get local wave numbers
!
      kx = spread(kx_fft(ipx*nx+1:ipx*nx+nx),2,ny)
      ky = spread(ky_fft(ipy*ny+1:ipy*ny+ny),1,nx)
!
!  Calculate 1/k^2, zero mean
!
      kappa = sqrt(kx**2 + ky**2)
!
!  Fourier transforms of x- and y-components on the boundary
!  Check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
      case ('bot')
        ! Potential field condition at the bottom
        do j=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          delta_z  = z(n1+j) - z(n1-j)
          exp_fact = exp(-kappa*delta_z)
          ! Determine potential field in ghost zones
          do i=iax,iaz
            tmp_re = f(l1:l2,m1:m2,n1+j,i)
            tmp_im = 0.0
            call fourier_transform_xy_xy(tmp_re,tmp_im)
            aa_re(:,:,i) = tmp_re*exp_fact
            aa_im(:,:,i) = tmp_im*exp_fact
          enddo
         ! Transform back
          do i=iax,iaz
            tmp_re = aa_re(:,:,i)
            tmp_im = aa_im(:,:,i)
            call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
            f(l1:l2,m1:m2,n1-j,i) = tmp_re
          enddo
        enddo
!
      case ('top')
        ! Potential field condition at the top
        do j=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          delta_z  = z(n2+j) - z(n2-j)
          exp_fact = exp(-kappa*delta_z)
          ! Determine potential field in ghost zones
          do i=iax,iaz
            tmp_re = f(l1:l2,m1:m2,n2-j,i)
            tmp_im = 0.0
            call fourier_transform_xy_xy(tmp_re,tmp_im)
            aa_re(:,:,i) = tmp_re*exp_fact
            aa_im(:,:,i) = tmp_im*exp_fact
          enddo
          ! Transform back
          do i=iax,iaz
            tmp_re = aa_re(:,:,i)
            tmp_im = aa_im(:,:,i)
            call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
            f(l1:l2,m1:m2,n2+j,i) = tmp_re
          enddo
        enddo
!
      case default
        call fatal_error('bc_aa_pot3', 'invalid argument', lfirst_proc_xy)
!
      endselect
!
!  The vector potential needs to be known outside of (l1:l2,m1:m2) as well
!
      call communicate_vect_field_ghosts(f,topbot)
!
!  Deallocate large arrays.
!
      if (allocated(aa_re)) deallocate(aa_re)
      if (allocated(aa_im)) deallocate(aa_im)
      if (allocated(kx)) deallocate(kx)
      if (allocated(ky)) deallocate(ky)
      if (allocated(kappa)) deallocate(kappa)
      if (allocated(exp_fact)) deallocate(exp_fact)
      if (allocated(tmp_re)) deallocate(tmp_re)
      if (allocated(tmp_im)) deallocate(tmp_im)
!
    endsubroutine bc_aa_pot3
!***********************************************************************
    subroutine bc_aa_pot2(f,topbot)
!
!  Potential field boundary condition
!
!  10-oct-06/tobi: Coded
!
      use Fourier, only: fourier_transform_xy_xy, fourier_transform_y_y
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
!
      real, dimension (:,:,:), allocatable :: aa_re,aa_im
      real, dimension (:,:), allocatable :: kx,ky,kappa
      real, dimension (:,:), allocatable :: tmp_re,tmp_im,fac
      integer :: i,j,stat
!
!  Allocate memory for large arrays.
!
      allocate(aa_re(nx,ny,iax:iaz),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot2', &
          'Could not allocate memory for aa_re',.true.)
      allocate(aa_im(nx,ny,iax:iaz),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot2', &
          'Could not allocate memory for aa_im',.true.)
      allocate(kx(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot2', &
          'Could not allocate memory for kx',.true.)
      allocate(ky(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot2', &
          'Could not allocate memory for ky',.true.)
      allocate(kappa(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot2', &
          'Could not allocate memory for kappa',.true.)
      allocate(tmp_re(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot2', &
          'Could not allocate memory for tmp_re',.true.)
      allocate(tmp_im(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot2', &
          'Could not allocate memory for tmp_im',.true.)
      allocate(fac(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot2', &
          'Could not allocate memory for fac',.true.)
!
!  Get local wave numbers
!
      if (nxgrid>1) then
        kx = spread(kx_fft(ipx*nx+1:ipx*nx+nx),2,ny)
        ky = spread(ky_fft(ipy*ny+1:ipy*ny+ny),1,nx)
      else
        kx(1,:) = 0.0
        ky(1,:) = ky_fft(ipy*ny+1:ipy*ny+ny)
      endif
!
!  Calculate 1/k^2, zero mean
!
      kappa = sqrt(kx**2 + ky**2)
!
!  Fourier transforms of x- and y-components on the boundary
!  Check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
      case ('bot')
        ! Potential field condition at the bottom
        do i=iax,iaz
          tmp_re = f(l1:l2,m1:m2,n1,i)
          tmp_im = 0.0
          if (nxgrid>1) then
            call fourier_transform_xy_xy(tmp_re,tmp_im)
          else
            call fourier_transform_y_y(tmp_re,tmp_im)
          endif
          aa_re(:,:,i) = tmp_re
          aa_im(:,:,i) = tmp_im
        enddo
        ! Determine potential field in ghost zones
        do j=1,nghost
          fac = exp(-j*kappa*dz)
          do i=iax,iaz
            tmp_re = fac*aa_re(:,:,i)
            tmp_im = fac*aa_im(:,:,i)
            if (nxgrid>1) then
              call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
            else
              call fourier_transform_y_y(tmp_re,tmp_im,linv=.true.)
            endif
            f(l1:l2,m1:m2,n1-j,i) = tmp_re
          enddo
        enddo
!
      case ('top')
        ! Potential field condition at the top
        do i=iax,iaz
          tmp_re = f(l1:l2,m1:m2,n2,i)
          tmp_im = 0.0
          if (nxgrid>1) then
            call fourier_transform_xy_xy(tmp_re,tmp_im)
          else
            call fourier_transform_y_y(tmp_re,tmp_im)
          endif
          aa_re(:,:,i) = tmp_re
          aa_im(:,:,i) = tmp_im
        enddo
        ! Determine potential field in ghost zones
        do j=1,nghost
          fac = exp(-j*kappa*dz)
          do i=iax,iaz
            tmp_re = fac*aa_re(:,:,i)
            tmp_im = fac*aa_im(:,:,i)
            if (nxgrid>1) then
              call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
            else
              call fourier_transform_y_y(tmp_re,tmp_im,linv=.true.)
            endif
            f(l1:l2,m1:m2,n2+j,i) = tmp_re
          enddo
        enddo
!
      case default
        call fatal_error('bc_aa_pot2', 'invalid argument', lfirst_proc_xy)
!
      endselect
!
!  The vector potential needs to be known outside of (l1:l2,m1:m2) as well
!
        call communicate_vect_field_ghosts(f,topbot)
!
!  Deallocate large arrays.
!
      if (allocated(aa_re)) deallocate(aa_re)
      if (allocated(aa_im)) deallocate(aa_im)
      if (allocated(kx)) deallocate(kx)
      if (allocated(ky)) deallocate(ky)
      if (allocated(kappa)) deallocate(kappa)
      if (allocated(tmp_re)) deallocate(tmp_re)
      if (allocated(tmp_im)) deallocate(tmp_im)
      if (allocated(fac)) deallocate(fac)
!
    endsubroutine bc_aa_pot2
!***********************************************************************
      subroutine bc_aa_pot(f,topbot)
!
!  Potential field boundary condition for magnetic vector potential at
!  bottom or top boundary (in z).
!
!  14-jun-2002/axel: adapted from similar
!   8-jul-2002/axel: introduced topbot argument
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=3) :: topbot
!
      real, dimension (:,:), allocatable :: f2,f3
      real, dimension (:,:,:), allocatable :: fz
      integer :: j, stat
!
!  Allocate memory for large arrays.
!
      allocate(f2(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot', &
          'Could not allocate memory for f2',.true.)
      allocate(f3(nx,ny),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot', &
          'Could not allocate memory for f3',.true.)
      allocate(fz(nx,ny,nghost+1),stat=stat)
      if (stat>0) call fatal_error('bc_aa_pot', &
          'Could not allocate memory for fz',.true.)
!
!  potential field condition
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  potential field condition at the bottom
!
      case ('bot')
        if (headtt) print*,'bc_aa_pot: pot-field bdry cond at bottom'
        if (mod(nxgrid,nygrid)/=0) &
             call fatal_error("bc_aa_pot", "pot-field doesn't work "//&
                 "with mod(nxgrid,nygrid)/=1", lfirst_proc_xy)
        do j=0,1
          f2=f(l1:l2,m1:m2,n1+1,iax+j)
          f3=f(l1:l2,m1:m2,n1+2,iax+j)
          call potential_field(fz,f2,f3,-1)
          f(l1:l2,m1:m2,1:n1,iax+j)=fz
        enddo
!
        f2=f(l1:l2,m1:m2,n1,iax)
        f3=f(l1:l2,m1:m2,n1,iay)
        call potentdiv(fz,f2,f3,-1)
        f(l1:l2,m1:m2,1:n1,iaz)=-fz
!
!  potential field condition at the top
!
      case ('top')
        if (headtt) print*,'bc_aa_pot: pot-field bdry cond at top'
        if (mod(nxgrid,nygrid)/=0) &
             call fatal_error("bc_aa_pot", "pot-field doesn't work "//&
                 "with mod(nxgrid,nygrid)/=1", lfirst_proc_xy)
        do j=0,1
          f2=f(l1:l2,m1:m2,n2-1,iax+j)
          f3=f(l1:l2,m1:m2,n2-2,iax+j)
          call potential_field(fz,f2,f3,+1)
          f(l1:l2,m1:m2,n2:mz,iax+j)=fz
        enddo
!
        f2=f(l1:l2,m1:m2,n2,iax)
        f3=f(l1:l2,m1:m2,n2,iay)
        call potentdiv(fz,f2,f3,+1)
        f(l1:l2,m1:m2,n2:mz,iaz)=-fz
      case default
        call fatal_error('bc_aa_pot', 'invalid argument', lfirst_proc_xy)
      endselect
!
      call communicate_vect_field_ghosts(f,topbot)
!
!  Deallocate large arrays.
!
      if (allocated(f2)) deallocate(f2)
      if (allocated(f3)) deallocate(f3)
      if (allocated(fz)) deallocate(fz)
!
      endsubroutine bc_aa_pot
!***********************************************************************
      subroutine potential_field(fz,f2,f3,irev)
!
!  solves the potential field boundary condition;
!  fz is the boundary layer, and f2 and f3 are the next layers inwards.
!  The condition is the same on the two sides.
!
!  20-jan-00/axel+wolf: coded
!  22-mar-00/axel: corrected sign (it is the same on both sides)
!  29-sep-06/axel: removed multiple calls, removed normalization, non-para
!
      use Fourier, only: fourier_transform_xy_xy
!
      real, dimension (:,:,:) :: fz
      real, dimension (:,:) :: f2,f3
      integer :: irev
!
      real, dimension (:,:), allocatable :: fac,kk,f1r,f1i,g1r,g1i
      real, dimension (:,:), allocatable :: f2r,f2i,f3r,f3i
      real :: delz
      integer :: i,stat
!
!  Allocate memory for large arrays.
!
      allocate(fac(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potential_field', &
          'Could not allocate memory for fac',.true.)
      allocate(kk(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potential_field', &
          'Could not allocate memory for kk',.true.)
      allocate(f1r(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potential_field', &
          'Could not allocate memory for f1r',.true.)
      allocate(f1i(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potential_field', &
          'Could not allocate memory for f1i',.true.)
      allocate(g1r(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potential_field', &
          'Could not allocate memory for g1r',.true.)
      allocate(g1i(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potential_field', &
          'Could not allocate memory for g1i',.true.)
      allocate(f2r(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potential_field', &
          'Could not allocate memory for f2r',.true.)
      allocate(f2i(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potential_field', &
          'Could not allocate memory for f2i',.true.)
      allocate(f3r(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potential_field', &
          'Could not allocate memory for f3r',.true.)
      allocate(f3i(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potential_field', &
          'Could not allocate memory for f3i',.true.)
!
!  initialize workspace
!
      f2r=f2; f2i=0
      f3r=f3; f3i=0
!
!  Transform; real and imaginary parts
!
      call fourier_transform_xy_xy(f2r,f2i)
      call fourier_transform_xy_xy(f3r,f3i)
!
!  define wave vector
!  calculate sqrt(k^2)
!
      kk=sqrt(spread(kx_fft(ipx*nx+1:ipx*nx+nx)**2,2,ny)+spread(ky_fft(ipy*ny+1:ipy*ny+ny)**2,1,nx))
!
!  one-sided derivative
!
      fac=1./(3.+2.*dz*kk)
      f1r=fac*(4.*f2r-f3r)
      f1i=fac*(4.*f2i-f3i)
!
!  set ghost zones
!
      do i=0,nghost
        delz=i*dz
        fac=exp(-kk*delz)
        g1r=fac*f1r
        g1i=fac*f1i
!
!  Transform back
!
        call fourier_transform_xy_xy(g1r,g1i,linv=.true.)
!
!  reverse order if irev=-1 (if we are at the bottom)
!
        if (irev==+1) fz(:,:,       i+1) = g1r
        if (irev==-1) fz(:,:,nghost-i+1) = g1r
      enddo
!
!  Deallocate large arrays.
!
      if (allocated(fac)) deallocate(fac)
      if (allocated(kk)) deallocate(kk)
      if (allocated(f1r)) deallocate(f1r)
      if (allocated(f1i)) deallocate(f1i)
      if (allocated(g1r)) deallocate(g1i)
      if (allocated(f2r)) deallocate(f2r)
      if (allocated(f2i)) deallocate(f2i)
      if (allocated(f3r)) deallocate(f3r)
      if (allocated(f3i)) deallocate(f3i)
!
    endsubroutine potential_field
!***********************************************************************
    subroutine potentdiv(fz,f2,f3,irev)
!
!  solves the divA=0 for potential field boundary condition;
!  f2 and f3 correspond to Ax and Ay (input) and fz corresponds to Ax (out)
!  In principle we could save some ffts, by combining with the potential
!  subroutine above, but this is now easier
!
!  22-mar-02/axel: coded
!  29-sep-06/axel: removed multiple calls, removed normalization, non-para
!   7-oct-06/axel: corrected sign for irev==+1.
!
      use Fourier, only: fourier_transform_xy_xy
!
      real, dimension (:,:,:) :: fz
      real, dimension (:,:) :: f2,f3
      integer :: irev
!
      real, dimension (:,:), allocatable :: fac,kk,kkkx,kkky
      real, dimension (:,:), allocatable :: f1r,f1i,g1r,g1i,f2r,f2i,f3r,f3i
      real, dimension (:), allocatable :: ky
      real, dimension (nx) :: kx
      real :: delz
      integer :: i,stat
!
!  Allocate memory for large arrays.
!
      allocate(fac(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potentdiv', &
          'Could not allocate memory for fac',.true.)
      allocate(kk(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potentdiv', &
          'Could not allocate memory for kk',.true.)
      allocate(kkkx(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potentdiv', &
          'Could not allocate memory for kkkx',.true.)
      allocate(kkky(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potentdiv', &
          'Could not allocate memory for kkky',.true.)
      allocate(f1r(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potentdiv', &
          'Could not allocate memory for f1r',.true.)
      allocate(f1i(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potentdiv', &
          'Could not allocate memory for f1i',.true.)
      allocate(g1r(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potentdiv', &
          'Could not allocate memory for g1r',.true.)
      allocate(g1i(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potentdiv', &
          'Could not allocate memory for g1i',.true.)
      allocate(f2r(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potentdiv', &
          'Could not allocate memory for f2r',.true.)
      allocate(f2i(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potentdiv', &
          'Could not allocate memory for f2i',.true.)
      allocate(f3r(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potentdiv', &
          'Could not allocate memory for f3r',.true.)
      allocate(f3i(nx,ny),stat=stat)
      if (stat>0) call fatal_error('potentdiv', &
          'Could not allocate memory for f3i',.true.)
      allocate(ky(nygrid),stat=stat)
      if (stat>0) call fatal_error('potentdiv', &
          'Could not allocate memory for ky',.true.)
!
      f2r=f2; f2i=0
      f3r=f3; f3i=0
!
!  Transform
!
      call fourier_transform_xy_xy(f2r,f2i)
      call fourier_transform_xy_xy(f3r,f3i)
!
!  define wave vector
!
      kx=cshift((/(i-nx/2,i=0,nx-1)/),+nx/2)*2*pi/Lx
      ky=cshift((/(i-nygrid/2,i=0,nygrid-1)/),+nygrid/2)*2*pi/Ly
!
!  calculate 1/k^2, zero mean
!
      kk=sqrt(spread(kx**2,2,ny)+spread(ky(ipy*ny+1:(ipy+1)*ny)**2,1,nx))
      kkkx=spread(kx,2,ny)
      kkky=spread(ky(ipy*ny+1:(ipy+1)*ny),1,nx)
!
!  calculate 1/kk
!
      kk(1,1)=1.
      fac=1./kk
      fac(1,1)=0.
!
      f1r=fac*(-kkkx*f2i-kkky*f3i)
      f1i=fac*(+kkkx*f2r+kkky*f3r)
!
!  set ghost zones
!
      do i=0,nghost
        delz=i*dz
        fac=exp(-kk*delz)
        g1r=fac*f1r
        g1i=fac*f1i
!
!  Transform back
!
        call fourier_transform_xy_xy(g1r,g1i,linv=.true.)
!
!  reverse order if irev=-1 (if we are at the bottom)
!  but reverse sign if irev=+1 (if we are at the top)
!
        if (irev==+1) fz(:,:,       i+1) = -g1r
        if (irev==-1) fz(:,:,nghost-i+1) = +g1r
      enddo
!
!  Deallocate large arrays.
!
      if (allocated(fac)) deallocate(fac)
      if (allocated(kk)) deallocate(kk)
      if (allocated(kkkx)) deallocate(kkkx)
      if (allocated(kkky)) deallocate(kkky)
      if (allocated(f1r)) deallocate(f1r)
      if (allocated(f1i)) deallocate(f1i)
      if (allocated(g1r)) deallocate(g1i)
      if (allocated(f2r)) deallocate(f2r)
      if (allocated(f2i)) deallocate(f2i)
      if (allocated(f3r)) deallocate(f3r)
      if (allocated(f3i)) deallocate(f3i)
      if (allocated(ky)) deallocate(ky)
!
    endsubroutine potentdiv
!***********************************************************************
    subroutine bc_wind_z(f,topbot,massflux)
!
!  Calculates u_0 so that rho*(u+u_0)=massflux.
!  Set 'win' for rho and
!  massflux can be set as fbcz1/2(rho) in run.in.
!
!  18-06-2008/bing: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j,ipt,ntb=-1
      real :: massflux,u_add
      real :: local_flux,local_mass
      real :: total_flux,total_mass
      real :: get_lf,get_lm
      integer :: nroot
!
      if (headtt) then
        print*,'bc_wind: Massflux',massflux
!
!   check wether routine can be implied
!
        if (.not.(lequidist(1) .and. lequidist(2))) &
            call fatal_error('bc_wind_z', &
            'non equidistant grid in x and y not implemented')
!
!   check for warnings
!
        if (.not. ldensity)  &
            call warning('bc_wind',"no defined density, using rho=1 ?")
      endif
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        ntb = n1
        nroot = 0
!
!  Top boundary.
!
      case ('top')
        ntb = n2
        nroot = ipz*nprocx*nprocy
!
!  Default.
!
      case default
        print*, "bc_wind: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
      local_flux=sum(exp(f(l1:l2,m1:m2,ntb,ilnrho))*f(l1:l2,m1:m2,ntb,iuz))
      local_mass=sum(exp(f(l1:l2,m1:m2,ntb,ilnrho)))
!
!  One  processor has to collect the data
!
      if (iproc/=nroot) then
        ! send to first processor at given height
        !
        call mpisend_real(local_flux,1,nroot,111+iproc)
        call mpisend_real(local_mass,1,nroot,211+iproc)
      else
        total_flux=local_flux
        total_mass=local_mass
        do i=0,nprocx-1
          do j=0,nprocy-1
            ipt = i+nprocx*j+ipz*nprocx*nprocy
            if (ipt/=nroot) then
              call mpirecv_real(get_lf,1,ipt,111+ipt)
              call mpirecv_real(get_lm,1,ipt,211+ipt)
              total_flux=total_flux+get_lf
              total_mass=total_mass+get_lm
            endif
          enddo
        enddo
!
!  Get u0 addition rho*(u+u0) = wind
!  rho*u + u0 *rho =wind
!  u0 = (wind-rho*u)/rho
!
        u_add = (massflux-total_flux) / total_mass
      endif
!
!  now distribute u_add
!
      if (iproc/=nroot) then
        call mpirecv_real(u_add,1,nroot,311+iproc)
      else
        do i=0,nprocx-1
          do j=0,nprocy-1
            ipt = i+nprocx*j+ipz*nprocx*nprocy
            if (ipt/=nroot) then
              call mpisend_real(u_add,1,ipt,311+ipt)
            endif
          enddo
        enddo
      endif
!
!  Set boundary
!
      f(l1:l2,m1:m2,ntb,iuz) =  f(l1:l2,m1:m2,ntb,iuz)+u_add
!
    endsubroutine bc_wind_z
!***********************************************************************
    subroutine bc_ADI_flux_z(f,topbot)
!
!  Constant flux boundary condition for temperature (called when bcz='c3')
!  at the bottom _only_ in the ADI case where hcond(n1)=hcond(x)
!  TT version: enforce dT/dz = - Fbot/K
!  30-jan-2009/dintrans: coded
!
      use SharedVariables, only: get_shared_variable
!
      real, pointer :: Fbot
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: tmp_x
      integer :: i, ierr
!
      call get_shared_variable('Fbot', Fbot, ierr)
      if (ierr/=0) call stop_it("bc_lnTT_flux_z: "//&
           "there was a problem when getting Fbot")
!
      if (headtt) print*,'bc_ADI_flux_z: Fbot, hcondADI, dz=', &
           Fbot, hcondADI, dz
!
      if (topbot=='bot') then
        tmp_x=-Fbot/hcondADI
        do i=1,nghost
          f(:,4,n1-i,ilnTT)=f(:,4,n1+i,ilnTT)-2.*i*dz*tmp_x
        enddo
      else
        call fatal_error('bc_ADI_flux_z', 'invalid argument')
      endif
!
    endsubroutine bc_ADI_flux_z
!***********************************************************************
    subroutine bc_force_ux_time(f, idz, j)
!
!  Set ux = ampl_forc*sin(k_forc*x)*cos(w_forc*t)
!
!  05-jun-2009/dintrans: coded from bc_force_uxy_sin_cos
!  Note: the ampl_forc, k_forc & w_forc run parameters are set in
!  'hydro' and shared using the 'shared_variables' module
!
      use SharedVariables, only : get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: idz, j, ierr
      real    :: kx
      real, pointer :: ampl_forc, k_forc, w_forc, x_forc, dx_forc
!
      if (headtt) then
        if (iuz == 0) call stop_it("BC_FORCE_UX_TIME: Bad idea...")
        if (Lx  == 0) call stop_it("BC_FORCE_UX_TIME: Lx cannot be 0")
        if (j /= iux) call stop_it("BC_FORCE_UX_TIME: only valid for ux")
      endif
      call get_shared_variable('ampl_forc', ampl_forc, ierr)
      if (ierr/=0) call stop_it("BC_FORCE_UX_TIME: "//&
           "there was a problem when getting ampl_forc")
      call get_shared_variable('k_forc', k_forc, ierr)
      if (ierr/=0) call stop_it("BC_FORCE_UX_TIME: "//&
           "there was a problem when getting k_forc")
      call get_shared_variable('w_forc', w_forc, ierr)
      if (ierr/=0) call stop_it("BC_FORCE_UX_TIME: "//&
           "there was a problem when getting w_forc")
      call get_shared_variable('x_forc', x_forc, ierr)
      if (ierr/=0) call stop_it("BC_FORCE_UX_TIME: "//&
           "there was a problem when getting x_forc")
      call get_shared_variable('dx_forc', dx_forc, ierr)
      if (ierr/=0) call stop_it("BC_FORCE_UX_TIME: "//&
           "there was a problem when getting dx_forc")
      if (headtt) print*, 'bc_force_ux_time: ampl_forc, k_forc, '//&
           'w_forc, x_forc, dx_forc=', ampl_forc, k_forc, w_forc, &
           x_forc, dx_forc
!
      if (k_forc /= impossible) then
        kx=2*pi/Lx*k_forc
        f(:,:,idz,j) = spread(ampl_forc*sin(kx*x)*cos(w_forc*t), 2, my)
      else
        f(:,:,idz,j) = spread(ampl_forc*exp(-((x-x_forc)/dx_forc)**2)*cos(w_forc*t), 2, my)
      endif
!
    endsubroutine bc_force_ux_time
!***********************************************************************
    subroutine bc_inlet_outlet_cyl(f,topbot,j,val)
!
! For pi/2 < y < 3pi/4,
! set r and theta velocity corresponding to a constant x-velocity
! and symmetric for lnrho/rho.
!
! Otherwise, set symmetric for velocities, and constant
! for lnrho/rho.
!
! NB! Assumes y to have the range 0 < y < 2pi
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j,i
      real, dimension(mcom) :: val
!
      select case (topbot)
      case ('bot')
        call fatal_error('bc_inlet_outlet_cyl', &
          'this boundary condition is not allowed for bottom boundary')
      case ('top')
        do m=m1,m2
          if (      (y(m)>=xyz0(2) +   Lxyz(2)/4)&
              .and. (y(m)<=xyz0(2) + 3*Lxyz(2)/4)) then
            if (j==iux) then
              f(l2,m,:,j) = cos(y(m))*val(j)
              do i=1,nghost; f(l2+i,m,:,j) = 2*f(l2,m,:,j) - f(l2-i,m,:,j); enddo
            elseif (j==iuy) then
              f(l2,m,:,j) = -sin(y(m))*val(j)
              do i=1,nghost; f(l2+i,m,:,j) = 2*f(l2,m,:,j) - f(l2-i,m,:,j); enddo
            elseif ((j==ilnrho) .or. (j==irho)) then
              do i=1,nghost; f(l2+i,m,:,j) = f(l2-i,m,:,j); enddo
            endif
!
          else
            if (j==iux) then
              do i=1,nghost; f(l2+i,m,:,j) = f(l2-i,m,:,j); enddo
            elseif (j==iuy) then
              do i=1,nghost; f(l2+i,m,:,j) = f(l2-i,m,:,j); enddo
            elseif ((j==ilnrho) .or. (j==irho)) then
              f(l2,m,:,j) = val(j)
              do i=1,nghost; f(l2+i,m,:,j) = 2*f(l2,m,:,j) - f(l2-i,m,:,j); enddo
            endif
          endif
        enddo
      endselect
!
    endsubroutine bc_inlet_outlet_cyl
!***********************************************************************
    subroutine bc_pp_hds_z_iso(f,topbot)
!
!  Boundary condition for pressure
!
!  This sets \partial_{z} p = \rho g_{z},
!  i.e. it enforces hydrostatic equlibrium at the boundary for the
!  pressure with an isothermal EOS.
!
!  16-dec-2009/dintrans: coded
!
      use Gravity, only: gravz
      use EquationOfState, only : cs20
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
      real    :: haut
      integer :: i
!
      haut=cs20/gravz
      if (topbot=='bot') then
        do i=1,nghost
          f(:,:,n1-i,ipp) = f(:,:,n1+i,ipp)-2.0*i*dz*f(:,:,n1,ipp)/haut
        enddo
      else
        do i=1,nghost
          f(:,:,n2+i,ipp) = f(:,:,n2-i,ipp)+2.0*i*dz*f(:,:,n2,ipp)/haut
        enddo
      endif
!
    endsubroutine bc_pp_hds_z_iso
!***********************************************************************
    subroutine bc_symset0der_z_v2(f,topbot,j)
!
!  This routine modified from bc_sym_z, but to a lower order.
!  Only available for z axis, activate with "0ds"
!  This is the routine to be used as regularity condition on the axis.
!
!  25-Oct-10/tijmen & bing: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j
!
      select case (topbot)
!
!  bottom (left end of the domain)
      case ('bot')
        f(:,:,n1,j)=(-18.*f(:,:,n1+1,j) &
                     +9.*f(:,:,n1+2,j) &
                     -2.*f(:,:,n1+3,j))/11.
!
        do i=1,nghost; f(:,:,n1-i,j)=f(:,:,n1+i,j); enddo
!
!  top (right end of the domain)
      case ('top')
        f(:,:,n2,j)=(+18.*f(:,:,n2-1,j) &
                     -9.*f(:,:,n2-2,j) &
                     +2.*f(:,:,n2-3,j))/11.
!
        do i=1,nghost; f(:,:,n2+i,j)=f(:,:,n2-i,j); enddo
!
      case default
        print*, "bc_symset0der_z_v2: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_symset0der_z_v2
!***********************************************************************
    subroutine bc_aa_pot_1D(f,topbot)
!
!  Computes a potential field extrapolation for a
!  1D magnetic field boundary with nprocx >= 1
!
!  27-Oct-10/bing: coded
!
      use Fourier, only: fourier_transform_other
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
      real, dimension (nxgrid) :: fft_az_r,fft_az_i,A_r,A_i,kx,exp_fact
      real, dimension (nxgrid) :: iay_global
      integer :: i,j,ipos,dir
!
      kx = abs(kx_fft(ipx*nx+1:ipx*nx+nx))
!
      select case (topbot)
!
!  bottom (left end of the domain)
      case ('bot')
        ipos = n1
        dir = -1
!
!  top (right end of the domain)
      case ('top')
        ipos = n2
        dir = 1
!
      case default
        print*, "bc_aa_pot_1D: ", topbot, " should be 'top' or 'bot'"
        ipos=1
        dir=0
!
      endselect
!
      if (nygrid>1) call fatal_error('bc_aa_pot_1D','only for nygrid=1')
!
      if (iproc==0) then
        iay_global(1:nx) = f(l1:l2,m1,ipos,iay)
        if (nprocx>1) then
          do j=1,nprocx-1
            call mpirecv_real(iay_global(j*nx+1:(j+1)*nx),nx,j,j*100)
          enddo
        endif
        fft_az_r=iay_global
        call fourier_transform_other(fft_az_r,fft_az_i)
      else
        if (nprocx>1) call mpisend_real(f(l1:l2,m1,ipos,iay),nx,0,iproc*100)
      endif
!
      do i=1,nghost
        if (iproc==0) then
!
          exp_fact = exp(-kx*(z(ipos+dir*i)-z(ipos)))
!
          A_r = exp_fact*fft_az_r
          A_i = exp_fact*fft_az_i
!
          call fourier_transform_other(A_r,A_i,linv=.true.)
!
          f(l1:l2,m1,ipos+dir*i,iay) = A_r(1:nx)
!
          if (nprocx>1) then
            do j=1,nprocx-1
              call mpisend_real(A_r(j*nx+1:(j+1)*nx),nx,j,j*100)
            enddo
          endif
        else
          if (nprocx>1) call mpirecv_real(f(l1:l2,m1,ipos+dir*i,iay),nx,0,iproc*100)
        endif
!
      enddo
!
    endsubroutine bc_aa_pot_1D
!***********************************************************************
    subroutine bc_ctz(f,topbot,j)
!
!  Set entropy to match temperature in the ghost zones to boundary value
!  value. Density ghost zones need to be calculated again here and corners
!  must be included to avoid NAN's.
!
!  13-feb-11/fred: check that 'ism' or 'a2' also set for bcz density.
!
      use EquationOfState, only: get_cv1,get_cp1
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j,k
      real :: cv1,cp1,cv,cp
!
      call get_cv1(cv1);cv=1./cv1
      call get_cp1(cp1);cp=1./cp1
!
      call bc_onesided_z(f,topbot,j-1)
!
      if (ldensity_nolog) then
        where (f(:,:,:,j-1)<=0.0) f(:,:,:,j-1)=tini*15.0
      endif
!
      if (.not.ldensity_nolog) &
          f(:,:,:,j-1)=exp(f(:,:,:,j-1))
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do k=1,3
          f(:,:,n1-k,j)=f(:,:,n1,j) &
              +(cp-cv)*(log(f(:,:,n1,j-1))-log(f(:,:,n1-k,j-1)))
        enddo
!
      case ('top')               ! top boundary
        do k=1,3
          f(:,:,n2+k,j)=f(:,:,n2,j) &
              +(cp-cv)*(log(f(:,:,n2,j-1))-log(f(:,:,n2+k,j-1)))
        enddo
!
      case default
        print*, "bc_ctz ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
      if (.not.ldensity_nolog) &
          f(:,:,:,j-1)=log(f(:,:,:,j-1))
!
    endsubroutine bc_ctz
!***********************************************************************
    subroutine bc_twist_xz(f,topbot,j)
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j,ix,iz
      real :: cx,cz,rad,radm,amp,amp0
!
      cx = (l1+l2)/2.
      cz = (n1+n2)/2.
!
      radm = (l1+l2)/16.
      amp0 = 1.
!
      if ((j < iux) .or. (j > iuz)) call fatal_error('bc_twist_xz', &
          'only for iux to iuz')
!
      select case (topbot)
      case ('bot')
        do ix=l1,l2
          do iz=n1,n2
            rad = sqrt( (cx-ix)**2. + (cz-iz)**2.)
            amp = amp0*rad/radm * exp(-rad/radm+1.)
            f(ix,1:m1,iz,iux) = -amp * cos( real(iz-nghost)/(1.*nz)*pi)
            f(ix,1:m1,iz,iuz) = amp * cos( real(ix-nghost)/(1.*nx)*pi)
            f(ix,1:m1,iz,iuy) = 0.
          enddo
        enddo
      case default
        call fatal_error('bc_twist_xz','only for bottom up to now')
      endselect
!
    endsubroutine  bc_twist_xz
!***********************************************************************
endmodule Boundcond
