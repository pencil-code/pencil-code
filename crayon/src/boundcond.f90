! $Id: boundcond.f90 13597 2010-04-06 10:25:38Z tavo.buk $
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
!  Boundary conditions in x except for periodic part handled by communication.
!  boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-`corners').
!
!   8-jul-02/axel: split up into different routines for x,y and z directions
!  11-nov-02/wolf: unified bot/top, now handled by loop
!  15-dec-06/wolf: Replaced "if (bcx1(1)=='she') then" by "any" command
!
      use EquationOfState
      use Shear
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      real, dimension (mcom) :: fbcx12
      real, dimension (mcom) :: fbcx2_12
      integer :: ivar1, ivar2, j, k, ip_ok, one
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
        if (headtt) print*, 'boundconds_x: no x-boundary'
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
          do k=1,2                ! loop over 'bot','top'
            if (k==1) then
              topbot='bot'; bc12=bcx1; fbcx12=fbcx1; fbcx2_12=fbcx1_2; ip_ok=0
            else
              topbot='top'; bc12=bcx2; fbcx12=fbcx2; fbcx2_12=fbcx2_2; ip_ok=nprocx-1
            endif
!
            do j=ivar1,ivar2
              if (ldebug) write(*,'(A,I1,A,I2,A,A)') ' bcx',k,'(',j,')=',bc12(j)
              if (ipx==ip_ok) then
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
                  ! BCX_DOC: implies $f"+f'/R=0$
                  call bc_cpc_x(f,topbot,j)
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
                case ('c1')
                  ! BCX_DOC: constant temperature (or maybe rather constant
                  ! BCX_DOC: conductive flux??)
                  if (j==iss)   call bc_ss_flux_x(f,topbot)
                case ('db')
                  ! BCX_DOC:
                  call bc_db_x(f,topbot,j)
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
                case ('sa2')
                  ! BCX_DOC: $(d/dr)(r B_{\phi}) = 0$ imposes
                  ! BCX_DOC: boundary condition on 2nd derivative of
                  ! BCX_DOC: $r A_{\phi}$. Same applies to $\theta$ component.
                  call bc_set_sa2_x(f,topbot,j)
                case ('pfc')
                  !BCX_DOC: perfect-conductor in spherical coordinate: $d/dr( A_r) + 2/r = 0$ .
                  call bc_set_pfc_x(f,topbot,j)
                case ('fix')
                  ! BCX_DOC: set boundary value [really??]
                  call bc_fix_x(f,topbot,j,fbcx12(j))
                case ('fil')
                  ! BCX_DOC: set boundary value from a file
                  call bc_file_x(f,topbot,j)
                case ('nil')
                case ('')
                  ! BCX_DOC: do nothing; assume that everything is set
                case default
                  bc%bcname=bc12(j)
                  bc%ivar=j
                  bc%location=(((k-1)*2)-1)   ! -1/1 for x bot/top
                  bc%value1=fbcx12(j)
                  bc%done=.false.
!
                  if (.not.bc%done) then
                    write(unit=errormsg,fmt='(A,A4,A,I3)') &
                         "No such boundary condition bcx1/2 = ", &
                         bc12(j), " for j=", j
                    call fatal_error("boundconds_x",errormsg)
                  endif
                endselect
              endif
            enddo
          enddo
        endif
      endselect
!
    endsubroutine boundconds_x
!***********************************************************************
    subroutine boundconds_y(f,ivar1_opt,ivar2_opt)
!
!  Boundary conditions in x except for periodic part handled by communication.
!  boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-`corners').
!
!   8-jul-02/axel: split up into different routines for x,y and z directions
!  11-nov-02/wolf: unified bot/top, now handled by loop
!
      use EquationOfState
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      real, dimension (mcom) :: fbcy12
      integer :: ivar1, ivar2, j, k, ip_ok
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
        if (headtt) print*,'boundconds_y: no y-boundary'
!
!  Boundary conditions in y
!
      case default
        do k=1,2                ! loop over 'bot','top'
          if (k==1) then
            topbot='bot'; bc12=bcy1; fbcy12=fbcy1; ip_ok=0
          else
            topbot='top'; bc12=bcy2; fbcy12=fbcy2; ip_ok=nprocy-1
          endif
!
          do j=ivar1,ivar2
            if (ldebug) write(*,'(A,I1,A,I2,A,A)') ' bcy',k,'(',j,')=',bc12(j)
            if (ipy==ip_ok) then
              select case (bc12(j))
              case ('p')
                ! BCY_DOC: periodic
                call bc_per_y(f,topbot,j)
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
                ! BCX_DOC: extrapolation in log [maintain a power law]
                call bcy_extrap_2_3(f,topbot,j)
              case ('der')
                ! BCY_DOC: set derivative on the boundary
                call bc_set_der_y(f,topbot,j,fbcy12(j))
              case ('')
               ! do nothing; assume that everything is set
              case default
                bc%bcname=bc12(j)
                bc%ivar=j
                bc%location=(((k-1)*4)-2)   ! -2/2 for y bot/top
                bc%done=.false.
!
                if (.not.bc%done) then
                  write(unit=errormsg,fmt='(A,A4,A,I3)') "No such boundary condition bcy1/2 = ", &
                       bc12(j), " for j=", j
                  call fatal_error("boundconds_y",errormsg)
                endif
              endselect
            endif
          enddo
        enddo
      endselect
!
    endsubroutine boundconds_y
!***********************************************************************
    subroutine boundconds_z(f,ivar1_opt,ivar2_opt)
!
!  Boundary conditions in x except for periodic part handled by communication.
!  boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-`corners').
!
!   8-jul-02/axel: split up into different routines for x,y and z directions
!  11-nov-02/wolf: unified bot/top, now handled by loop
!
      use EquationOfState
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      real, dimension (mcom) :: fbcz12, fbcz12_1, fbcz12_2, fbcz_zero=0.
      integer :: ivar1, ivar2, j, k, ip_ok
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
        if (headtt) print*,'boundconds_z: no z-boundary'
!
!  Boundary conditions in z
!
      case default
        !call get_shared_variable('Fbot',Fbot,ierr)
        !if (ierr/=0) call stop_it("boundconds_z: "//&
        !     "there was a problem when getting Fbot")
        !call get_shared_variable('Ftop',Ftop,ierr)
        !if (ierr/=0) call stop_it("boundconds_z: "//&
        !     "there was a problem when getting Fbot")
        !call get_shared_variable('FbotKbot',FbotKbot,ierr)
        !if (ierr/=0) call stop_it("boundconds_z: "//&
        !     "there was a problem when getting FbotKbot")
        !call get_shared_variable('FtopKtop',FtopKtop,ierr)
        !if (ierr/=0) call stop_it("boundconds_z: "//&
        !     "there was a problem when getting FtopKtop")
        do k=1,2                ! loop over 'bot','top'
          if (k==1) then
            topbot='bot'
            bc12=bcz1
            fbcz12=fbcz1
            fbcz12_1=fbcz1_1
            fbcz12_2=fbcz1_2
            ip_ok=0
            !Ftopbot=Fbot
            !FtopbotK=FbotKbot
          else
            topbot='top'
            bc12=bcz2
            fbcz12=fbcz2
            fbcz12_1=fbcz2_1
            fbcz12_2=fbcz2_2
            ip_ok=nprocz-1
            !Ftopbot=Ftop
            !FtopbotK=FtopKtop
          endif
!
          do j=ivar1,ivar2
            if (ldebug) write(*,'(A,I1,A,I2,A,A)') ' bcz',k,'(',j,')=',bc12(j)
            if (ipz==ip_ok) then
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
              case ('db')
                ! BCZ_DOC: complex
                ! BCZ_DOC:
                call bc_db_z(f,topbot,j)
              case ('e1')
                ! BCZ_DOC: extrapolation
                call bc_extrap_2_1(f,topbot,j)
              case ('e2')
                ! BCZ_DOC: extrapolation
                call bc_extrap_2_2(f,topbot,j)
              case ('b1')
                ! BCZ_DOC: extrapolation with zero value (improved 'a')
                call bc_extrap0_2_0(f,topbot,j)
              case ('b2')
                ! BCZ_DOC: extrapolation with zero value (improved 'a')
                call bc_extrap0_2_1(f,topbot,j)
              case ('b3')
                ! BCZ_DOC: extrapolation with zero value (improved 'a')
                call bc_extrap0_2_2(f,topbot,j)
              case ('fBs')
                ! BCZ_DOC: frozen-in B-field (s)
                call bc_frozen_in_bb(topbot,j)
                call bc_sym_z(f,+1,topbot,j) ! symmetry
              case ('fB')
                ! BCZ_DOC: frozen-in B-field (a2)
                call bc_frozen_in_bb(topbot,j)
                call bc_sym_z(f,-1,topbot,j,REL=.true.) ! antisymm wrt boundary
              case ('1')
                ! BCZ_DOC: f=1 (for debugging)
                call bc_one_z(f,topbot,j)
              case ('set')
                ! BCZ_DOC: set boundary value
                call bc_sym_z(f,-1,topbot,j,REL=.true.,val=fbcz12)
              case ('der')
                ! BCZ_DOC: set derivative on the boundary
                call bc_set_der_z(f,topbot,j,fbcz12(j))
              case ('ovr')
                ! BCZ_DOC: set boundary value
                call bc_overshoot_z(f,fbcz12,topbot,j)
              case ('ouf')
                ! BCZ_DOC: allow outflow, but no inflow (experimental)
                call bc_outflow_z(f,topbot,j)
              case ('win')
                ! BCZ_DOC: forces massflux given as
                ! BCZ_DOC: $\Sigma \rho_i ( u_i + u_0) = \textrm{fbcz1/2}(\rho)$
                if (j==ilnrho) then
                   call bc_wind_z(f,topbot,fbcz12(j))     !
                   call bc_sym_z(f,-1,topbot,j,REL=.true.)!  'a2'
                   call bc_sym_z(f,+1,topbot,iuz)         !  's'
                endif
              case ('cop')
                ! BCZ_DOC: copy value of last physical point to all ghost cells
                call bc_copy_z(f,topbot,j)
              case ('nil')
                ! do nothing; assume that everything is set
              case default
                bc%bcname=bc12(j)
                bc%ivar=j
                bc%location=(((k-1)*6)-3)   ! -3/3 for z bot/top
                bc%value1=fbcz12_1(j)
                bc%value2=fbcz12_2(j)
                bc%done=.false.
!
                if (.not.bc%done) then
                  write(unit=errormsg,fmt='(A,A4,A,I3)') "No such boundary condition bcz1/2 = ", &
                       bc12(j), " for j=", j
                  call fatal_error("boundconds_z",errormsg)
                endif
              endselect
            endif
          enddo
        enddo
      endselect
!
    endsubroutine boundconds_z
!***********************************************************************
    subroutine bc_per_x(f,topbot,j)
!
!  periodic boundary condition
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
        print*, "bc_per_x: ", topbot, " should be `top' or `bot'"
!
      endselect
!
    endsubroutine bc_per_x
!***********************************************************************
    subroutine bc_per_y(f,topbot,j)
!
!  periodic boundary condition
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
        print*, "bc_per_y: ", topbot, " should be `top' or `bot'"
!
      endselect
!
    endsubroutine bc_per_y
!***********************************************************************
    subroutine bc_per_z(f,topbot,j)
!
!  periodic boundary condition
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
        print*, "bc_per_z: ", topbot, " should be `top' or `bot'"
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
        if (present(val)) f(l1,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(l1-i,:,:,j)=2*f(l1,:,:,j)+sgn*f(l1+i,:,:,j); enddo
        else
          do i=1,nghost; f(l1-i,:,:,j)=              sgn*f(l1+i,:,:,j); enddo
          if (sgn<0) f(l1,:,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(l2,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(l2+i,:,:,j)=2*f(l2,:,:,j)+sgn*f(l2-i,:,:,j); enddo
        else
          do i=1,nghost; f(l2+i,:,:,j)=              sgn*f(l2-i,:,:,j); enddo
          if (sgn<0) f(l2,:,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case default
        print*, "bc_sym_x: ", topbot, " should be `top' or `bot'"
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
        dxR=-dx/x(l2)
        i=-0; f(l2+i,:,:,j)=0.
        i=-1; f(l2+i,:,:,j)=-(1.-.5*dxR)*f(l2-i,:,:,j)/(1.+.5*dxR)
        extra1=(1.+.5*dxR)*f(l2+i,:,:,j)+(1.-.5*dxR)*f(l2-i,:,:,j)
        i=-2; f(l2+i,:,:,j)=(-(1.-   dxR)*f(l2-i,:,:,j)+16.*extra1)/(1.+dxR)
        extra2=(1.+dxR)*f(l2+i,:,:,j)+(1.-dxR)*f(l2-i,:,:,j)-10.*extra1
        i=-3; f(l2+i,:,:,j)=(-(2.-3.*dxR)*f(l2-i,:,:,j)+27.*extra2)/(2.+3.*dxR)
!
      case ('top')               ! top boundary
        dxR=dx/x(l2)
        i=0; f(l2+i,:,:,j)=0.
        i=1; f(l2+i,:,:,j)=-(1.-.5*dxR)*f(l2-i,:,:,j)/(1.+.5*dxR)
        extra1=(1.+.5*dxR)*f(l2+i,:,:,j)+(1.-.5*dxR)*f(l2-i,:,:,j)
        i=2; f(l2+i,:,:,j)=(-(1.-   dxR)*f(l2-i,:,:,j)+16.*extra1)/(1.+dxR)
        extra2=(1.+dxR)*f(l2+i,:,:,j)+(1.-dxR)*f(l2-i,:,:,j)-10.*extra1
        i=3; f(l2+i,:,:,j)=(-(2.-3.*dxR)*f(l2-i,:,:,j)+27.*extra2)/(2.+3.*dxR)
!
      case default
        print*, "bc_cpc_x: ", topbot, " should be `top' or `bot'"
!
      endselect
!
    endsubroutine bc_cpc_x
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
        print*, "bc_symset_x: ", topbot, " should be `top' or `bot'"
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
        print*, "bc_symset0der_x: ", topbot, " should be `top' or `bot'"
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
!         f(l1,:,:,j)=(2.*x(l1+1)*f(l1+1,:,:,j)-.5*x(l1+2)*f(l1+2,:,:,j))/(1.5*x(l1))
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
!         f(l2,:,:,j)=(2.*x(l2-1)*f(l2-1,:,:,j)-.5*x(l2-2)*f(l2-2,:,:,j))/(1.5*x(l2))
        endif
!
      case default
        print*, "bc_slope_x: ", topbot, " should be `top' or `bot'"
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
        print*, "bc_slope_x: ", topbot, " should be `top' or `bot'"
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
        print*, "bc_overshoot_x: ", topbot, " should be `top' or `bot'"
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
        print*, "bc_overshoot_z: ", topbot, " should be `top' or `bot'"
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
        print*, "bc_antis_x: ", topbot, " should be `top' or `bot'"
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
        if (present(val)) f(l1:l2,m1,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,m1-i,:,j)=2*f(:,m1,:,j)+sgn*f(:,m1+i,:,j); enddo
        else
          do i=1,nghost; f(:,m1-i,:,j)=              sgn*f(:,m1+i,:,j); enddo
          if (sgn<0) f(:,m1,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(l1:l2,m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,m2+i,:,j)=2*f(:,m2,:,j)+sgn*f(:,m2-i,:,j); enddo
        else
          do i=1,nghost; f(:,m2+i,:,j)=              sgn*f(:,m2-i,:,j); enddo
          if (sgn<0) f(:,m2,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case default
        print*, "bc_sym_y: ", topbot, " should be `top' or `bot'"
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
        print*, "bc_symset_y: ", topbot, " should be `top' or `bot'"
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
        print*, "bc_symset0der_y: ", topbot, " should be `top' or `bot'"
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
        if (present(val)) f(l1:l2,m1:m2,n1,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,:,n1-i,j)=2*f(:,:,n1,j)+sgn*f(:,:,n1+i,j); enddo
        else
          do i=1,nghost; f(:,:,n1-i,j)=              sgn*f(:,:,n1+i,j); enddo
          if (sgn<0) f(:,:,n1,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(l1:l2,m1:m2,n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,:,n2+i,j)=2*f(:,:,n2,j)+sgn*f(:,:,n2-i,j); enddo
        else
          do i=1,nghost; f(:,:,n2+i,j)=              sgn*f(:,:,n2-i,j); enddo
          if (sgn<0) f(:,:,n2,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case default
        print*, "bc_sym_z: ", topbot, " should be `top' or `bot'"
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
        print*, "bc_symset0der_z: ", topbot, " should be `top' or `bot'"
!
      endselect
!
    endsubroutine bc_symset0der_z
!***********************************************************************
    subroutine bc_set_der_x(f,topbot,j,val)
!
!  Sets the derivative on the boundary to a given value
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
        call warning('bc_set_der_x',topbot//" should be `top' or `bot'")
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
        do i=1,nghost;f(l1-i,:,:,j)=val; enddo
      case ('top')               ! top boundary
        do i=1,nghost; f(l2+i,:,:,j)=val; enddo
      case default
        call warning('bc_fix_x',topbot//" should be `top' or `bot'")
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
        call warning('bc_fix_x',topbot//" should be `top' or `bot'")
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
    subroutine bc_set_pfc_x(f,topbot,j)
!
! In spherical polar coordinate system,
! at a radial boundary set : $A_{\theta} = 0$ and $A_{phi} = 0$,
! and demand $div A = 0$ gives the condition on $A_r$ to be
! $d/dr( A_r) + 2/r = 0$ . This subroutine sets this condition of
! $j$ the component of f. As this is related to setting the
! perfect conducting boundary condition we call this "pfc".
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
        call warning('bc_set_pfc_x',topbot//" should be `top' or `bot'")
!
      endselect
!
    endsubroutine bc_set_pfc_x
!***********************************************************************
    subroutine bc_set_sa2_x(f,topbot,j)
! To set the boundary condition:
! d_r(r B_{\phi} = 0 we need to se
! (d_r)^2(r A_{\theta}) = 0 which sets the condition 'a2'
! on r A_{\theta} and vice-versa for A_{\phi}
!
!  3-Dec-2009/dhruba: coded
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
        call warning('bc_set_sa2_x',topbot//" should be `top' or `bot'")
!
      endselect
!
    endsubroutine bc_set_sa2_x
! **********************************************************************
    subroutine bc_set_der_y(f,topbot,j,val)
!
!  Sets the derivative on the boundary to a given value
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
        call warning('bc_set_der_y',topbot//" should be `top' or `bot'")
!
      endselect
!
    endsubroutine bc_set_der_y
!***********************************************************************
    subroutine bc_set_der_z(f,topbot,j,val)
!
!  Sets the derivative on the boundary to a given value
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
        call warning('bc_set_der_z',topbot//" should be `top' or `bot'")
!
      endselect
!
    endsubroutine bc_set_der_z
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
        print*, "bc_van_x: ", topbot, " should be `top' or `bot'"
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
        print*, "bc_van_y: ", topbot, " should be `top' or `bot'"
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
        print*, "bc_van_z: ", topbot, " should be `top' or `bot'"
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
!   5-apr-03/axel: coded
!   7-jan-09/axel: corrected
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
        print*, "bc_onesided_x ", topbot, " should be `top' or `bot'"
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
!   5-apr-03/axel: coded
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
        print*, "bc_onesided_x_old ", topbot, " should be `top' or `bot'"
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
!   5-apr-03/axel: coded
!   7-jan-09/axel: corrected
!   26-jan-09/nils: adapted from bc_onesided_x
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
        print*, "bc_onesided_7 ", topbot, " should be `top' or `bot'"
!
      endselect
!
    endsubroutine bc_onesided_y
!***********************************************************************
    subroutine bc_onesided_z(f,topbot,j)
!
!  One-sided conditions.
!  These expressions result from combining Eqs(207)-(210), astro-ph/0109497,
!  corresponding to (9.207)-(9.210) in Ferriz-Mas proceedings.
!
!   5-apr-03/axel: coded
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
        print*, "bc_onesided_z ", topbot, " should be `top' or `bot'"
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
!   19-jun-03/wolf: coded
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
        print*, "bc_extrap_2_1: ", topbot, " should be `top' or `bot'"
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
!   19-jun-03/wolf: coded
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
        print*, "bcx_extrap_2_1: ", topbot, " should be `top' or `bot'"
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
!   19-jun-03/wolf: coded
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
        print*, "bcy_extrap_2_1: ", topbot, " should be `top' or `bot'"
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
!   19-jun-03/wolf: coded
!    1-jul-03/axel: introduced abbreviations n1p4,n2m4
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
        print*, "bc_extrap_2_2: ", topbot, " should be `top' or `bot'"
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
!   19-jun-03/wolf: coded
!    1-jul-03/axel: introduced abbreviations n1p4,n2m4
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
        print*, "bcx_extrap_2_2: ", topbot, " should be `top' or `bot'"
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
!    1-jul-03/axel: introduced abbreviations n1p4,n2m4
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
        print*, "bcy_extrap_2_2: ", topbot, " should be `top' or `bot'"
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
!   18-dec-08/wlad: coded
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
        print*, "bcy_extrap_2_3: ", topbot, " should be `top' or `bot'"
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
!    9-oct-03/wolf: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      select case (topbot)
!
!       case ('bot')               ! bottom boundary
!         f(:,:,n1  ,j)= 0.       ! set bdry value=0 (indep of initcond)
!         f(:,:,n1-1,j)=- 3*f(:,:,n1+1,j)+  f(:,:,n1+2,j)
!         f(:,:,n1-2,j)=- 8*f(:,:,n1+1,j)+3*f(:,:,n1+2,j)
!         f(:,:,n1-3,j)=-15*f(:,:,n1+1,j)+6*f(:,:,n1+2,j)
!
!       case ('top')               ! top boundary
!         f(:,:,n2  ,j)= 0.       ! set bdry value=0 (indep of initcond)
!         f(:,:,n2+1,j)=- 3*f(:,:,n2-1,j)+  f(:,:,n2-2,j)
!         f(:,:,n2+2,j)=- 8*f(:,:,n2-1,j)+3*f(:,:,n2-2,j)
!         f(:,:,n2+3,j)=-15*f(:,:,n2-1,j)+6*f(:,:,n2-2,j)
!
!! Nyquist-filtering
      case ('bot')               ! bottom boundary
        f(:,:,n1  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n1-1,j)=(1/11.)*(-17*f(:,:,n1+1,j)- 9*f(:,:,n1+2,j)+ 8*f(:,:,n1+3,j))
        f(:,:,n1-2,j)=      2*(- 2*f(:,:,n1+1,j)-   f(:,:,n1+2,j)+   f(:,:,n1+3,j))
        f(:,:,n1-3,j)=(3/11.)*(-27*f(:,:,n1+1,j)-13*f(:,:,n1+2,j)+14*f(:,:,n1+3,j))
!
      case ('top')               ! top boundary
        f(:,:,n2  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n2+1,j)=(1/11.)*(-17*f(:,:,n2-1,j)- 9*f(:,:,n2-2,j)+ 8*f(:,:,n2-3,j))
        f(:,:,n2+2,j)=      2*(- 2*f(:,:,n2-1,j)-   f(:,:,n2-2,j)+   f(:,:,n2-3,j))
        f(:,:,n2+3,j)=(3/11.)*(-27*f(:,:,n2-1,j)-13*f(:,:,n2-2,j)+14*f(:,:,n2-3,j))
!
! !! Nyquist-transparent
!       case ('bot')               ! bottom boundary
!         f(:,:,n1  ,j)=0.        ! set bdry value=0 (indep of initcond)
!         f(:,:,n1-1,j)=(1/11.)*(-13*f(:,:,n1+1,j)-14*f(:,:,n1+2,j)+10*f(:,:,n1+3,j))
!         f(:,:,n1-2,j)=(1/11.)*(-48*f(:,:,n1+1,j)-17*f(:,:,n1+2,j)+20*f(:,:,n1+3,j))
!         f(:,:,n1-3,j)=         - 7*f(:,:,n1+1,j)- 4*f(:,:,n1+2,j)+ 4*f(:,:,n1+3,j)
!
!       case ('top')               ! top boundary
!         f(:,:,n2  ,j)=0.        ! set bdry value=0 (indep of initcond)
!         f(:,:,n2+1,j)=(1/11.)*(-13*f(:,:,n2-1,j)-14*f(:,:,n2-2,j)+10*f(:,:,n2-3,j))
!         f(:,:,n2+2,j)=(1/11.)*(-48*f(:,:,n2-1,j)-17*f(:,:,n2-2,j)+20*f(:,:,n2-3,j))
!         f(:,:,n2+3,j)=         - 7*f(:,:,n2-1,j)- 4*f(:,:,n2-2,j)+ 4*f(:,:,n2-3,j)
!
      case default
        print*, "bc_extrap0_2_0: ", topbot, " should be `top' or `bot'"
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
!    9-oct-03/wolf: coded
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
        print*, "bc_extrap0_2_1: ", topbot, " should be `top' or `bot'"
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
!    9-oct-03/wolf: coded
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
        f(:,:,n1  ,j)= 0.       ! set bdry value=0 (indep of initcond)
        f(:,:,n1-1,j)=0.2   *(                 -  4*f(:,:,n1+2,j)- 3*f(:,:,n1+3,j)+ 3*f(:,:,n1p4,j))
        f(:,:,n1-2,j)=0.2   *(- 2*f(:,:,n1+1,j)-  9*f(:,:,n1+2,j)- 6*f(:,:,n1+3,j)+ 7*f(:,:,n1p4,j))
        f(:,:,n1-3,j)=1./35.*(-33*f(:,:,n1+1,j)-108*f(:,:,n1+2,j)-68*f(:,:,n1+3,j)+87*f(:,:,n1p4,j))
!
      case ('top')               ! top boundary
        f(:,:,n2  ,j)= 0.       ! set bdry value=0 (indep of initcond)
        f(:,:,n2+1,j)=0.2   *(                 -  4*f(:,:,n2-2,j)- 3*f(:,:,n2-3,j)+ 3*f(:,:,n2m4,j))
        f(:,:,n2+2,j)=0.2   *(- 2*f(:,:,n2-1,j)-  9*f(:,:,n2-2,j)- 6*f(:,:,n2-3,j)+ 7*f(:,:,n2m4,j))
        f(:,:,n2+3,j)=1./35.*(-33*f(:,:,n2-1,j)-108*f(:,:,n2-2,j)-68*f(:,:,n2-3,j)+87*f(:,:,n2m4,j))
!
      case default
        print*, "bc_extrap0_2_2: ", topbot, " should be `top' or `bot'"
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
!   18-dec-08/wlad: coded
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
        print*, "bcx_extrap_2_3: ", topbot, " should be `top' or `bot'"
!
      endselect
!
    endsubroutine bcx_extrap_2_3
!***********************************************************************
    subroutine bc_db_z(f,topbot,j)
!
!  ``One-sided'' boundary condition for density.
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
!  ``One-sided'' boundary condition for density.
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
        print*, "bc_one_x: ",topbot, " should be `top' or `bot'"
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
        print*, "bc_one_y: ", topbot, " should be `top' or `bot'"
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
        print*, "bc_one_z: ", topbot, " should be `top' or `bot'"
!
      endselect
!
    endsubroutine bc_one_z
!***********************************************************************
    subroutine bc_ss_flux_x(f,topbot)
!
!  constant flux boundary condition for entropy (called when bcx='c1')
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
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is processor dependent)
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
!  calculate Fbot/(K*cs2)
!
!       cs2_yz=cs20*exp(gamma_m1*(f(l1,:,:,ilnrho)-lnrho0)+cv1*f(l1,:,:,iss))

! Both, bottom and tom boundary conditions are corrected for linear density

        if (ldensity_nolog) then
           cs2_yz=cs20*exp(gamma_m1*(log(f(l1,:,:,ilnrho))-lnrho0)+gamma*f(l1,:,:,iss))
        else
           cs2_yz=cs20*exp(gamma_m1*(f(l1,:,:,ilnrho)-lnrho0)+gamma*f(l1,:,:,iss))
        endif
        tmp_yz=FbotKbot/cs2_yz
!
!  enforce ds/dx + gamma_m1/gamma*dlnrho/dx = - gamma_m1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
!         f(l1-i,:,:,iss)=f(l1+i,:,:,iss)+(cp-cv)* &
           if (ldensity_nolog) then
              f(l1-i,:,:,iss)=f(l1+i,:,:,iss)+gamma_m1/gamma* &
                   (log(f(l1+i,:,:,ilnrho))-log(f(l1-i,:,:,ilnrho))+2*i*dx*tmp_yz)
           else
              f(l1-i,:,:,iss)=f(l1+i,:,:,iss)+gamma_m1/gamma* &
                   (f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho)+2*i*dx*tmp_yz)
           endif
        enddo
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
        if (headtt) print*,'bc_ss_flux_z: FtopKtop=',FtopKtop
!
!  calculate Ftop/(K*cs2)
!
        if (ldensity_nolog) then
           cs2_yz=cs20*exp(gamma_m1*(log(f(l2,:,:,ilnrho))-lnrho0)+gamma*f(l2,:,:,iss))
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
                   (f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho)-2*i*dx*tmp_yz)
           else
              f(l2+i,:,:,iss)=f(l2-i,:,:,iss)+gamma_m1/gamma* &
                   (log(f(l2-i,:,:,ilnrho))-log(f(l2+i,:,:,ilnrho))-2*i*dx*tmp_yz)
           endif
!          f(l1-i,:,:,iss)=f(l1+i,:,:,iss)+gamma_m1/gamma* &
!              (f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho)+2*i*dx*tmp_yz)
        enddo
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
        print*, "bc_zero_x: ", topbot, " should be `top' or `bot'"
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
        print*, "bc_zero_z: ", topbot, " should be `top' or `bot'"
!
      endselect
!
    endsubroutine bc_zero_z
!***********************************************************************
    subroutine bc_outflow_z(f,topbot,j)
!
!  Outflow boundary conditions.
!
!  If the velocity vector points out of the box, the velocity boundary
!  condition is set to 's', otherwise it is set to 'a'.
!
!  12-aug-2007/anders: implemented
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
          if (f(ix,iy,n1,j)<=0.0) then  ! 's'
            do i=1,nghost; f(ix,iy,n1-i,j)=+f(ix,iy,n1+i,j); enddo
          else                          ! 'a'
            do i=1,nghost; f(ix,iy,n1-i,j)=-f(ix,iy,n1+i,j); enddo
            f(ix,iy,n1,j)=0.0
          endif
        enddo; enddo
!
!  Top boundary.
!
      case ('top')
        do iy=1,my; do ix=1,mx
          if (f(ix,iy,n2,j)>=0.0) then  ! 's'
            do i=1,nghost; f(ix,iy,n2+i,j)=+f(ix,iy,n2-i,j); enddo
          else                          ! 'a'
            do i=1,nghost; f(ix,iy,n2+i,j)=-f(ix,iy,n2-i,j); enddo
            f(ix,iy,n2,j)=0.0
          endif
        enddo; enddo
!
!  Default.
!
      case default
        print*, "bc_outflow_z: ", topbot, " should be `top' or `bot'"
!
      endselect
!
    endsubroutine bc_outflow_z
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
        print*, "bc_copy_z: ", topbot, " should be `top' or `bot'"
!
      endselect
!
    endsubroutine bc_copy_x
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
        print*, "bc_copy_z: ", topbot, " should be `top' or `bot'"
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
      use SharedVariables
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
        print*, "bc_frozen_in_bb: ", topbot, " should be `top' or `bot'"
      endselect
!
      lfirstcall=.false.
!
    endsubroutine bc_frozen_in_bb
!***********************************************************************
    subroutine bc_wind_z(f,topbot,massflux)
!
!  Calculates u_0 so that rho*(u+u_0)=massflux
!  massflux can be set as fbcz1/2(rho) in run.in
!
!  18-06-2008/bing: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,ipt,ntb
      real :: massflux,u_add
      real :: local_flux,local_mass
      real :: total_flux,total_mass
      real :: get_lf,get_lm
!
      if (headtt) then
         print*,'bc_wind: Massflux',massflux
!
!   check wether routine can be implied
!
         if (nprocx .gt. 1)  &
              call stop_it('bc_wind: nprocx > 1 not yet implemented')
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
!
!  Top boundary.
!
       case ('top')
         ntb = n2
!
!  Default.
!
      case default
        print*, "bc_wind: ", topbot, " should be `top' or `bot'"
!
      endselect
!
      local_flux=sum(exp(f(l1:l2,m1:m2,ntb,ilnrho))*f(l1:l2,m1:m2,ntb,iuz))
      local_mass=sum(exp(f(l1:l2,m1:m2,ntb,ilnrho)))
!
!  One  processor has to collect the data
!
      if (ipy .ne. 0) then
         ! send to first processor at given height
         !
         call mpisend_real(local_flux,1,ipz*nprocy,111+iproc)
         call mpisend_real(local_mass,1,ipz*nprocy,211+iproc)
      else
         do i=1,nprocy-1
            ipt=ipz*nprocy+i
            call mpirecv_real(get_lf,1,ipt,111+ipt)
            call mpirecv_real(get_lm,1,ipt,211+ipt)
            total_flux=total_flux+get_lf
            total_mass=total_mass+get_lm
         enddo
         total_flux=total_flux+local_flux
         total_mass=total_mass+local_mass
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
      if (ipy .eq. 0) then
         do i=1,nprocy-1
            ipt=ipz*nprocy+i
            call mpisend_real(u_add,1,ipt,311+ipt)
         enddo
      else
         call mpirecv_real(u_add,1,ipz*nprocy,311+iproc)
      endif
!
!  Set boundary
!
      f(l1:l2,m1:m2,ntb,iuz) =  f(l1:l2,m1:m2,ntb,iuz)+u_add
!
     endsubroutine bc_wind_z
!***********************************************************************
endmodule Boundcond
