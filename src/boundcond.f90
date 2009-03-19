! $iD: BOundcond.f90 10428 2009-02-24 10:14:46Z sven.bingert $

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   boundcond.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!

!!!  Module for boundary conditions. Extracted from (no)mpicomm, since
!!!  all non-periodic (external) boundary conditions require the same
!!!  code for serial and parallel runs.

module Boundcond

  use Cdata
  use Messages
  use Mpicomm
  use Viscosity, Only:llambda_effect,Lambda_V0,Lambda_Omega

  implicit none

  private

  public :: boundconds, boundconds_x, boundconds_y, boundconds_z
  public :: bc_per_x, bc_per_y, bc_per_z
  public :: update_ghosts
  public :: nscbc_boundtreat

  integer, pointer :: iglobal_gg

  contains

!***********************************************************************
    subroutine boundconds(f,ivar1_opt,ivar2_opt)
!
!  Apply boundary conditions in all three directions.
!  Note that we _must_ call boundconds_{x,y,z} in this order, or edges and
!  corners will not be OK.
!
!  10-oct-02/wolf: coded
!
      use Cparam
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt

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
      use Special, only: special_boundconds
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
      select case(nxgrid)
!
      case(1)
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
          if (.not. lstart) then
            call boundcond_shear(f,ivar1,ivar2)
          endif
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
                select case(bc12(j))
                case ('p')
                  ! BCX_DOC: periodic
                  call bc_per_x(f,topbot,j)
                case ('s')
                  ! BCX_DOC: symmetry, $f_{N+i}=f_{N-i}$;
                  ! BCX_DOC: implies $f'(x_N)=f'''(x_0)=0$
                  call bc_sym_x(f,+1,topbot,j)
                case ('ss')
                  ! BCX_DOC: symmetry [???]
                  call bc_symset_x(f,+1,topbot,j)
                case ('a')
                  ! BCX_DOC: antisymmetry, $f_{N+i}=-f_{N-i}$;
                  ! BCX_DOC: implies $f(x_N)=f''(x_0)=0$
                  call bc_sym_x(f,-1,topbot,j)
                case ('a2')
                  ! BCX_DOC: antisymmetry relative to boundary value,
                  ! BCX_DOC: $f_{N+i}=2 f_{N}-f_{N-i}$;
                  ! BCX_DOC: implies $f''(x_0)=0$
                  call bc_sym_x(f,-1,topbot,j,REL=.true.)
                case ('v')
                  ! BCX_DOC: vanishing third derivative
                  call bc_van_x(f,topbot,j)
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
                  !if (j==iss) 
                  !if (j==ilnTT)  then
                  !  force_lower_bound='cT'
                  !  force_upper_bound='cT'
                  !  call bc_force_x(f,-1,topbot,j)
                  !endif
                case ('c1')
                  ! BCX_DOC: constant temperature (or maybe rather constant
                  ! BCX_DOC: conductive flux??)
                  if (j==iss)   call bc_ss_flux_x(f,topbot)
                  if (j==ilnTT) call bc_lnTT_flux_x(f,topbot)
                case ('sT')
                  ! BCX_DOC: symmetric temperature, $T_{N-i}=T_{N+i}$;
                  ! BCX_DOC: implies $T'(x_N)=T'''(x_0)=0$
                  if (j==iss) call bc_ss_stemp_x(f,topbot)
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
                  ! BCX_DOC: set boundary value [really??]
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
                  !BCX_DOC: top hat jet profile in spherical coordinate. 
                  !Defined only for the bottom boundary 
                  call bc_set_jethat_x(f,j,topbot,fbcx12,fbcx2_12)
                case ('spd')
                  ! BCX_DOC:  sets $d(rA_{\alpha})/dr = \mathtt{fbcx12(j)}$
                  call bc_set_spder_x(f,topbot,j,fbcx12(j))
                case('sfr')
                  ! BCX_DOC: stress-free boundary condition for spherical coordinate system. 
                  call bc_set_sfree_x(f,topbot,j)
                case('nfr')
                  ! BCX_DOC: Normal-field bc for spherical coordinate system.
                  ! BCX_DOC: Some people call this the ``(angry) hedgehog bc''.
                  call bc_set_nfr_x(f,topbot,j)
                case('pfc')
                  !BCX_DOC: perfect-conductor in spherical coordinate: $d/dr( A_r) + 2/r = 0$ . 
                  call bc_set_pfc_x(f,topbot,j)
                 case ('fix')
                  ! BCX_DOC: set boundary value [really??]
                  call bc_fix_x(f,topbot,j,fbcx12(j))
                 case ('fil')
                  ! BCX_DOC: set boundary value from a file
                  call bc_file_x(f,topbot,j)
                case('cfb')
                  ! BCZ_DOC: radial centrifugal balance 
                  if (lcylindrical_coords) then
                    call bc_lnrho_cfb_r_iso(f,topbot,j)
                  else
                    print*,'not implemented for other than cylindrical'
                    stop
                  endif
                case ('')
                  ! BCX_DOC: do nothing; assume that everything is set
                case default
                  bc%bcname=bc12(j)
                  bc%ivar=j
                  bc%location=(((k-1)*2)-1)   ! -1/1 for x bot/top
                  bc%value1=fbcx12(j)
                  bc%done=.false.

                  call special_boundconds(f,bc)

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
      use Special, only: special_boundconds
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
      select case(nygrid)
!
      case(1)
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
              select case(bc12(j))
              case ('p')
                ! BCY_DOC: periodic
                call bc_per_y(f,topbot,j)
              case ('s')
                ! BCY_DOC: symmetry symmetry, $f_{N+i}=f_{N-i}$;
                  ! BCX_DOC: implies $f'(y_N)=f'''(y_0)=0$
                call bc_sym_y(f,+1,topbot,j)
              case ('ss')
                ! BCY_DOC: symmetry [???]
                call bc_symset_y(f,+1,topbot,j)
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
              case ('f')
                ! BCY_DOC: freeze value
                ! tell other modules not to change boundary value
                call bc_freeze_var_y(topbot,j)
                call bc_sym_y(f,-1,topbot,j,REL=.true.) ! antisymm wrt boundary
              case ('fg')
                ! BCY_DOC: ``freeze'' value, i.e. maintain initial
                !  value at boundary, also mantaining the 
                !  ghost zones at the initial coded value, i.e., 
                !  keep the gradient frozen as well
                call bc_freeze_var_y(topbot,j)
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
              case('sfr')
                  ! BCY_DOC: stress-free boundary condition for spherical coordinate system. 
                call bc_set_sfree_y(f,topbot,j)
              case('nfr')
                  ! BCY_DOC: Normal-field bc for spherical coordinate system.
                  ! BCY_DOC: Some people call this the ``(angry) hedgehog bc''.
                call bc_set_nfr_y(f,topbot,j)
              case('pfc')
                  !BCY_DOC: perfect conducting boundary condition along $\theta$ boundary  
                call bc_set_pfc_y(f,topbot,j)
              case ('')
               ! do nothing; assume that everything is set
              case default
                bc%bcname=bc12(j)
                bc%ivar=j
                bc%location=(((k-1)*4)-2)   ! -2/2 for y bot/top
                bc%done=.false.

                if (lspecial) call special_boundconds(f,bc)

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
!!      use Entropy, only: hcond0,hcond1,Fbot,FbotKbot,Ftop,FtopKtop,chi, &
!!                         lmultilayer,lheatc_chiconst
      use Special, only: special_boundconds
      !use Density
      use EquationOfState
      !use SharedVariables, only : get_shared_variable
      !use Mpicomm,         only : stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
      !real, pointer :: Fbot, Ftop, FbotKbot, FtopKtop
!
      real, dimension (mcom) :: fbcz12, fbcz12_1, fbcz12_2
      !real :: Ftopbot,FtopbotK
      integer :: ivar1, ivar2, j, k, ip_ok, ierr
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
      select case(nzgrid)
!
      case(1)
        if (headtt) print*,'boundconds_z: no z-boundary'
!
!  Boundary conditions in z
!
      case default
        !call get_shared_variable('Fbot',Fbot,ierr)
        !if (ierr/=0) call stop_it("boundcond_z: "//&
        !     "there was a problem when getting Fbot")
        !call get_shared_variable('Ftop',Ftop,ierr)
        !if (ierr/=0) call stop_it("boundcond_z: "//&
        !     "there was a problem when getting Fbot")
        !call get_shared_variable('FbotKbot',FbotKbot,ierr)
        !if (ierr/=0) call stop_it("boundcond_z: "//&
        !     "there was a problem when getting FbotKbot")
        !call get_shared_variable('FtopKtop',FtopKtop,ierr)
        !if (ierr/=0) call stop_it("boundcond_z: "//&
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
              select case(bc12(j))
              case ('0')
                ! BCZ_DOC: zero value in ghost zones, free value on boundary
                call bc_zero_z(f,topbot,j)
              case ('p')
                ! BCZ_DOC: periodic
                call bc_per_z(f,topbot,j)
              case ('s')
                ! BCZ_DOC: symmetry
                call bc_sym_z(f,+1,topbot,j)
              case ('a')
                ! BCZ_DOC: antisymmetry
                call bc_sym_z(f,-1,topbot,j)
              case ('a2')
                ! BCZ_DOC: antisymmetry relative to boundary value
                call bc_sym_z(f,-1,topbot,j,REL=.true.)
              case ('a3')
                ! BCZ_DOC: a2 - wiggles
                call bc_asym3(f,topbot,j)
              case ('v')
                ! BCZ_DOC: vanishing third derivative
                call bc_van_z(f,topbot,j)
              case ('v3')
                ! BCZ_DOC: vanishing third derivative
                call bc_van3rd_z(f,topbot,j)
              case ('1s')
                ! BCZ_DOC: one-sided
                call bc_onesided_z(f,topbot,j)
              case ('c1')
                ! BCZ_DOC: complex
                if (j==iss) call bc_ss_flux(f,topbot)
                if (j==iaa) call bc_aa_pot(f,topbot)
                if (j==ilnTT) call bc_lnTT_flux_z(f,topbot)
              case ('c3')
                ! BCZ_DOC: constant flux at the bottom with a variable hcond
                if (j==ilnTT) call bc_ADI_flux_z(f,topbot)
              case ('pot')
                ! BCZ_DOC: potential magnetic field
                if (j==iaa) call bc_aa_pot2(f,topbot)
              case ('pwd')
                ! BCZ_DOC: a variant of `pot'
                if (j==iaa) call bc_aa_pot3(f,topbot)
              case ('d2z')
                ! BCZ_DOC: 
                call bc_del2zero(f,topbot,j)
              case ('hds')
                ! BCZ_DOC: hydrostatic equilibrium with 
                !          a high-frequency filter
                if (llocal_iso) then 
                  call bc_lnrho_hdss_z_liso(f,topbot)
                else
                  call bc_lnrho_hdss_z_iso(f,topbot)
                endif
              case ('cT')
                ! BCZ_DOC: constant temp.
                ! BCZ_DOC: 
                if (j==ilnrho) call bc_lnrho_temp_z(f,topbot)
                call bc_ss_temp_z(f,topbot)
                !if (j==iss) then
                !   if (pretend_lnTT) then
                !       force_lower_bound='cT'
                !       force_upper_bound='cT'
                !      call bc_force_z(f,-1,topbot,j)                      
                !   else
                ! endif
                !endif
                !if (j==ilnTT)  then
                !   force_lower_bound='cT'
                !   force_upper_bound='cT'
                !  call bc_force_z(f,-1,topbot,j)
                !endif
              case ('cT2')
                ! BCZ_DOC: constant temp. (keep lnrho)
                ! BCZ_DOC: 
                if (j==iss)   call bc_ss_temp2_z(f,topbot)
              case ('hs')
                ! BCZ_DOC: hydrostatic equilibrium
                if (llocal_iso) then !non local
                  if (j==ilnrho) call bc_lnrho_hds_z_liso(f,topbot)
!                 if (j==iss)    call bc_lnrho_hydrostatic_z(f,topbot)
                else
                  if (j==ilnrho) call bc_lnrho_hds_z_iso(f,topbot)
!                 if (j==iss)    call bc_lnrho_hydrostatic_z(f,topbot)
                endif
              case ('cp')
                ! BCZ_DOC: constant pressure
                ! BCZ_DOC: 
                if (j==ilnrho) call bc_lnrho_pressure_z(f,topbot)
              case ('sT')
                ! BCZ_DOC: symmetric temp.
                ! BCZ_DOC: 
                if (j==iss) call bc_ss_stemp_z(f,topbot)
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

                if (lspecial) call special_boundconds(f,bc)

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
    subroutine nscbc_boundtreat(f,df)
!
!  Boundary treatment of the df-array. 
!
!  This is a way to impose (time-
!  dependent) boundary conditions by solving a so-called characteristic
!  form of the fluid equations on the boundaries, as opposed to setting 
!  actual values of the variables in the f-array. The method is called 
!  Navier-Stokes characteristic boundary conditions (NSCBC).
!  The current implementation only solves a simplified version of the
!  equations, namely a set of so-called local one-dimensional inviscid
!  (LODI) relations. This means that transversal and viscous terms are
!  dropped on the boundaries.
!
!  The treatment should be done after the y-z-loop, but before the Runge-
!  Kutta solver adds to the f-array.
!
!   7-jul-08/arne: coded.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df

      intent(inout) :: f
      intent(inout) :: df
!
      if (nscbc1(1) /= '' .or. nscbc2(1) /= '') &
          call nscbc_boundtreat_xyz(f,df,1)
      if (nscbc1(2) /= '' .or. nscbc2(2) /= '') &
          call nscbc_boundtreat_xyz(f,df,2)
      if (nscbc1(3) /= '' .or. nscbc2(3) /= '') &
          call nscbc_boundtreat_xyz(f,df,3)
     
    endsubroutine nscbc_boundtreat
!***********************************************************************
    subroutine nscbc_boundtreat_xyz(f,df,j)
!
!   NSCBC boundary treatment.
!   j = 1, 2 or 3 for x, y or z-boundaries respectively.
!
!   7-jul-08/arne: coded.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=nscbc_len), dimension(3) :: bc12
      character (len=3) :: topbot
      integer j,k,direction,ip_ok,ip_test
      real, dimension(mcom) :: valx,valy,valz

      intent(inout) :: f
      intent(inout) :: df
      intent(in)    :: j
!
      do k=1,2                ! loop over 'bot','top'
        if (k==1) then
          topbot='bot'; bc12(j)=nscbc1(j);!val=bt_val1(j)
          valx=fbcx1; valy=fbcy1; valz=fbcz1; ip_ok=0
        else
          topbot='top'; bc12(j)=nscbc2(j);!val=bt_val2(j)
          valx=fbcx2; valy=fbcy2; valz=fbcz2
          if (j==1) ip_ok=nprocx-1
          if (j==2) ip_ok=nprocy-1
          if (j==3) ip_ok=nprocz-1
        endif
        if (j==1) ip_test=ipx
        if (j==2) ip_test=ipy
        if (j==3) ip_test=ipz
!
!  Check if this is a physical boundary
!
        if (ip_test==ip_ok) then
          select case(bc12(j))
          case('part_ref_outlet')
!   Partially reflecting outlet.
            if (j==1) then 
              call bc_nscbc_prf_x(f,df,topbot,.false.)
            elseif (j==2) then 
              call bc_nscbc_prf_y(f,df,topbot,.false.)
            elseif (j==3) then 
              call bc_nscbc_prf_z(f,df,topbot,.false.)
            endif
          case('part_ref_inlet')
!   Partially reflecting inlet, ie. impose a velocity u_t.
            if (j==1) then 
              direction = 1
              call bc_nscbc_prf_x(f,df,topbot,.true.,linlet=.true.,&
                  u_t=valx(direction))
            elseif (j==2) then 
              direction = 2
              call bc_nscbc_prf_y(f,df,topbot,.true.,linlet=.true.,&
                  u_t=valy(direction))
            elseif (j==3) then 
              direction = 3
              call bc_nscbc_prf_z(f,df,topbot,.true.,linlet=.true.,&
                  u_t=valz(direction))
            endif
          case('ref_inlet')
!   Partially reflecting inlet, ie. impose a velocity u_t.
            if (j==1) then 
              direction = 1
              call bc_nscbc_prf_x(f,df,topbot,.false.,linlet=.true.,&
                  u_t=valx(direction))
            elseif (j==2) then 
              direction = 2
              call bc_nscbc_prf_y(f,df,topbot,.false.,linlet=.true.,&
                  u_t=valy(direction))
            elseif (j==3) then 
              direction = 3
              call bc_nscbc_prf_z(f,df,topbot,.false.,linlet=.true.,&
                  u_t=valz(direction))
            endif
          case('subsonic_inflow')
! Subsonic inflow 
            if (j==1) then
              call bc_nscbc_subin_x(f,df,topbot,val=valx)
            elseif (j==2) then
            endif
          case('subson_nref_outflow')
            if (j==1) then
              call bc_nscbc_nref_subout_x(f,df,topbot)
            elseif (j==2) then
              call bc_nscbc_nref_subout_y(f,df,topbot)
            endif
          case('')
!   Do nothing.
          case('none')
            print*,'nscbc_boundtreat_xyz: doing nothing!'
!   Do nothing.
          case default
            call fatal_error("nscbc_boundtreat_xyz",'You must specify nscbc bouncond!')
          endselect
        endif
      end do
    end subroutine
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
      select case(topbot)

      case('bot')               ! bottom boundary
        if (nprocx==1) f(1:l1-1,:,:,j) = f(l2i:l2,:,:,j)

      case('top')               ! top boundary
        if (nprocx==1) f(l2+1:mx,:,:,j) = f(l1:l1i,:,:,j)

      case default
        print*, "bc_per_x: ", topbot, " should be `top' or `bot'"

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
      select case(topbot)

      case('bot')               ! bottom boundary
        if (nprocy==1) f(:,1:m1-1,:,j) = f(:,m2i:m2,:,j)

      case('top')               ! top boundary
        if (nprocy==1) f(:,m2+1:my,:,j) = f(:,m1:m1i,:,j)

      case default
        print*, "bc_per_y: ", topbot, " should be `top' or `bot'"

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
      select case(topbot)

      case('bot')               ! bottom boundary
        if (nprocz==1) f(:,:,1:n1-1,j) = f(:,:,n2i:n2,j)

      case('top')               ! top boundary
        if (nprocz==1) f(:,:,n2+1:mz,j) = f(:,:,n1:n1i,j)

      case default
        print*, "bc_per_z: ", topbot, " should be `top' or `bot'"

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

      select case(topbot)

      case('bot')               ! bottom boundary
        if (present(val)) f(l1,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(l1-i,:,:,j)=2*f(l1,:,:,j)+sgn*f(l1+i,:,:,j); enddo
        else
          do i=1,nghost; f(l1-i,:,:,j)=              sgn*f(l1+i,:,:,j); enddo
          if (sgn<0) f(l1,:,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif

      case('top')               ! top boundary
        if (present(val)) f(l2,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(l2+i,:,:,j)=2*f(l2,:,:,j)+sgn*f(l2-i,:,:,j); enddo
        else
          do i=1,nghost; f(l2+i,:,:,j)=              sgn*f(l2-i,:,:,j); enddo
          if (sgn<0) f(l2,:,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif

      case default
        print*, "bc_sym_x: ", topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_sym_x
!***********************************************************************
    subroutine bc_symset_x(f,sgn,topbot,j,rel,val)
!
! FIXME: Get documentation right
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

      select case(topbot)

      case('bot')               ! bottom boundary
        if (present(val)) f(l1,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(l1-i,:,:,j)=2*f(l1,:,:,j)+sgn*f(l1+i,:,:,j); enddo
        else
          do i=1,nghost; f(l1-i,:,:,j)=              sgn*f(l1+i,:,:,j); enddo
          f(l1,:,:,j)=(4.*f(l1+1,:,:,j)-f(l1+2,:,:,j))/3.
        endif

      case('top')               ! top boundary
        if (present(val)) f(l2,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(l2+i,:,:,j)=2*f(l2,:,:,j)+sgn*f(l2-i,:,:,j); enddo
        else
          do i=1,nghost; f(l2+i,:,:,j)=              sgn*f(l2-i,:,:,j); enddo
          f(l2,:,:,j)=(4.*f(l2-1,:,:,j)-f(l2-2,:,:,j))/3.
        endif

      case default
        print*, "bc_symset_x: ", topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_symset_x
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

      select case(topbot)

      case('bot')               ! bottom boundary
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

      case('top')               ! top boundary
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

      case default
        print*, "bc_slope_x: ", topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_slope_x
!***********************************************************************
    subroutine bc_dr0_x(f,slope,topbot,j,rel,val)

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
      integer :: l1_4=l1+4, l1_5=l1+5, l1_6=l1+6
      integer :: l2_4=l2-4, l2_5=l2-5, l2_6=l2-6
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif

      select case(topbot)

      case('bot')               ! bottom boundary
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

      case('top')               ! top boundary
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

      case default
        print*, "bc_slope_x: ", topbot, " should be `top' or `bot'"

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
      select case(topbot)
!
!  bottom
!
      case('bot')               ! bottom boundary
        do i=1,nghost
          f(l1-i,:,:,j)=f(l1+i,:,:,j)*exp((x(l1-i)-x(l1+i))/dist(j))
        enddo
!
!  top
!
      case('top')               ! top boundary
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
      select case(topbot)
!
!  bottom
!
      case('bot')               ! bottom boundary
        do i=1,nghost
          f(:,:,n1-i,j)=f(:,:,n1+i,j)*exp((z(n1-i)-z(n1+i))/dist(j))
        enddo
!
!  top
!
      case('top')               ! top boundary
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
! FIXME: Get documentation right
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
      if (present(rel)) then; relative=rel; else; relative=.false.; endif

      select case(topbot)

      case('bot')               ! bottom boundary
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

      case('top')               ! top boundary
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

      case default
        print*, "bc_antis_x: ", topbot, " should be `top' or `bot'"

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

      select case(topbot)

      case('bot')               ! bottom boundary
        if (present(val)) f(l1:l2,m1,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,m1-i,:,j)=2*f(:,m1,:,j)+sgn*f(:,m1+i,:,j); enddo
        else
          do i=1,nghost; f(:,m1-i,:,j)=              sgn*f(:,m1+i,:,j); enddo
          if (sgn<0) f(:,m1,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif

      case('top')               ! top boundary
        if (present(val)) f(l1:l2,m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,m2+i,:,j)=2*f(:,m2,:,j)+sgn*f(:,m2-i,:,j); enddo
        else
          do i=1,nghost; f(:,m2+i,:,j)=              sgn*f(:,m2-i,:,j); enddo
          if (sgn<0) f(:,m2,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif

      case default
        print*, "bc_sym_y: ", topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_sym_y
!***********************************************************************
    subroutine bc_symset_y(f,sgn,topbot,j,rel,val)
!
! FIXME: Get documentation right
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

      select case(topbot)

      case('bot')               ! bottom boundary
        if (present(val)) f(l1:l2,m1,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,m1-i,:,j)=2*f(:,m1,:,j)+sgn*f(:,m1+i,:,j); enddo
        else
          do i=1,nghost; f(:,m1-i,:,j)=              sgn*f(:,m1+i,:,j); enddo
          f(:,m1,:,j)=(4.*f(:,m1+1,:,j)-f(:,m1+2,:,j))/3.
        endif

      case('top')               ! top boundary
        if (present(val)) f(l1:l2,m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,m2+i,:,j)=2*f(:,m2,:,j)+sgn*f(:,m2-i,:,j); enddo
        else
          do i=1,nghost; f(:,m2+i,:,j)=              sgn*f(:,m2-i,:,j); enddo
          f(:,m2,:,j)=(4.*f(:,m2-1,:,j)-f(:,m2-2,:,j))/3.
        endif

      case default
        print*, "bc_symset_y: ", topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_symset_y
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

      select case(topbot)

      case('bot')               ! bottom boundary
        if (present(val)) f(l1:l2,m1:m2,n1,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,:,n1-i,j)=2*f(:,:,n1,j)+sgn*f(:,:,n1+i,j); enddo
        else
          do i=1,nghost; f(:,:,n1-i,j)=              sgn*f(:,:,n1+i,j); enddo
          if (sgn<0) f(:,:,n1,j) = 0. ! set bdry value=0 (indep of initcond)
        endif

      case('top')               ! top boundary
        if (present(val)) f(l1:l2,m1:m2,n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,:,n2+i,j)=2*f(:,:,n2,j)+sgn*f(:,:,n2-i,j); enddo
        else
          do i=1,nghost; f(:,:,n2+i,j)=              sgn*f(:,:,n2-i,j); enddo
          if (sgn<0) f(:,:,n2,j) = 0. ! set bdry value=0 (indep of initcond)
        endif

      case default
        print*, "bc_sym_z: ", topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_sym_z
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

      integer :: i

      select case(topbot)

      case('bot')               ! bottom boundary
        do i=1,nghost; f(l1-i,:,:,j) = f(l1+i,:,:,j) - 2*i*dx*val; enddo

      case('top')               ! top boundary
        do i=1,nghost; f(l2+i,:,:,j) = f(l2-i,:,:,j) + 2*i*dx*val; enddo

      case default
        call warning('bc_set_der_x',topbot//" should be `top' or `bot'")

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

      real, intent (in) :: val
      integer :: i

      select case(topbot)

      case('bot')               ! bottom boundary
        do i=1,nghost;f(l1-i,:,:,j)=val; enddo
      case('top')               ! top boundary
        do i=1,nghost; f(l2+i,:,:,j)=val; enddo
      case default
        call warning('bc_fix_x',topbot//" should be `top' or `bot'")

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
      real, dimension (mx,my,mz,mvar) :: bc_file_x_array
      integer, intent (in) :: j
      integer :: i,lbc0,lbc1,lbc2
      real :: lbc,frac
      logical, save :: lbc_file_x=.true.

      if (lbc_file_x) then
        if (lroot) then
          print*,'opening bc_file_x.dat'
          open(9,file=trim(directory_snap)//'/bc_file_x.dat',form='unformatted')
          read(9,end=99) bc_file_x_array
          close(9)
        endif
        lbc_file_x=.false.
      endif

      select case(topbot)
!
!  x - Udrift_bc*t = dx * (ix - Udrift_bc*t/dx)
!
      case('bot')               ! bottom boundary
        lbc=Udrift_bc*t*dx_1(1)+1.
        lbc0=int(lbc)
        frac=mod(lbc,real(lbc0))
        lbc1=mx+mod(-lbc0,mx)
        lbc2=mx+mod(-lbc0-1,mx)
        do i=1,nghost
          f(l1-i,:,:,j)=(1-frac)*bc_file_x_array(lbc1,:,:,j) &
                           +frac*bc_file_x_array(lbc2,:,:,j)
        enddo
      case('top')               ! top boundary
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

      endselect
!
      goto 98
99    continue
      if (lroot) print*,'need file with dimension: ',mx,my,mz,mvar
      call stop_it("boundary file bc_file_x.dat not found")
98  endsubroutine bc_file_x
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

      real, intent (in) :: val
      integer :: i

      if (lspherical_coords)then
        select case(topbot)
        case('bot')               ! bottom boundary
        do i=1,nghost
          f(l1-i,:,:,j)=f(l1+i,:,:,j)-2*i*dx*(val-f(l1,:,:,j)*r1_mn(1))
        enddo
      case('top')               ! top boundary
        do i=1,nghost
          f(l2+i,:,:,j)=f(l2-i,:,:,j)+2*i*dx*(val-f(l2,:,:,j)*r1_mn(nx))
        enddo

      case default
        call warning('bc_set_spder_x',topbot//" should be `top' or `bot'")

      endselect
    else
      call stop_it('Boundary condition spder is valid only in spherical coordinate system')
    endif
!
    endsubroutine bc_set_spder_x
! **********************************************************************
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
  

      select case(topbot)

      case('bot')               ! bottom boundary
! The coding assumes we are using 6-th order centered finite difference for our
! derivatives. 
        f(l1-1,:,:,j)= f(l1+1,:,:,j) +  2.*60.*f(l1,:,:,j)*dx/(45.*x(l1))
        f(l1-2,:,:,j)= f(l1+2,:,:,j) +  2.*60.*f(l1,:,:,j)*dx/(9.*x(l1))
        f(l1-3,:,:,j)= f(l1+3,:,:,j) +  2.*60.*f(l1,:,:,j)*dx/x(l1)
      case('top')               ! top boundary
        f(l2+1,:,:,j)= f(l2-1,:,:,j) -  2.*60.*f(l2,:,:,j)*dx/(45.*x(l2))
        f(l2+2,:,:,j)= f(l2-2,:,:,j) -  2.*60.*f(l2,:,:,j)*dx/(9.*x(l2))
        f(l2+3,:,:,j)= f(l2-3,:,:,j) -  2.*60.*f(l2,:,:,j)*dx/(x(l2))

      case default
        call warning('bc_set_pfc_x',topbot//" should be `top' or `bot'")

      endselect
!
    endsubroutine bc_set_pfc_x
!***********************************************************************
    subroutine bc_set_nfr_x(f,topbot,j)
!
! Normal-field (or angry-hedgehog) boundary condition for spherical
! coordinate system. 
! d_r(A_{\theta}) = -A_{\theta}/r  with A_r = 0 sets B_{r} to zero
! in spherical coordinate system. 
! (compare with next subroutine sfree )
!
!  25-Aug-2007/dhruba: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j


      select case(topbot)

      case('bot')               ! bottom boundary
! The coding assumes we are using 6-th order centered finite difference for our
! derivatives. 
        f(l1-1,:,:,j)= f(l1+1,:,:,j) +  60.*f(l1,:,:,j)*dx/(45.*x(l1))
        f(l1-2,:,:,j)= f(l1+2,:,:,j) +  60.*f(l1,:,:,j)*dx/(9.*x(l1))
        f(l1-3,:,:,j)= f(l1+3,:,:,j) +  60.*f(l1,:,:,j)*dx/(x(l1))
      case('top')               ! top boundary
        f(l2+1,:,:,j)= f(l2-1,:,:,j) -  60.*f(l1,:,:,j)*dx/(45.*x(l2))
        f(l2+2,:,:,j)= f(l2-2,:,:,j) -  60.*f(l1,:,:,j)*dx/(9.*x(l2))
        f(l2+3,:,:,j)= f(l2-3,:,:,j) -  60.*f(l1,:,:,j)*dx/(x(l2))

      case default
        call warning('bc_set_nfr_x',topbot//" should be `top' or `bot'")

      endselect
!
    endsubroutine bc_set_nfr_x
! **********************************************************************
    subroutine bc_set_sfree_x(f,topbot,j)
!
! Stress-free boundary condition for spherical coordinate system. 
! d_r(u_{\theta}) = u_{\theta}/r  with u_r = 0 sets S_{r \theta}
! component of the strain matrix to be zero in spherical coordinate system. 
! This subroutine sets only the first part of this boundary condition for 'j'-th
! component of f. 
! Lambda effect : stresses due to Lambda effect are added to the stress-tensor. 
! For rotation along the z direction and also for not very strong rotation such
! that the breaking of rotational symmetry is only due to gravity, the only 
! new term is appears in the r-phi component. This implies that this term
! affects only the boundary condition of u_{\phi} for the radial boundary. 

!  25-Aug-2007/dhruba: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real, dimension (my,mz) :: boundary_value
      integer, intent (in) :: j

      select case(topbot)
!
! bottom boundary
!
      case('bot')
!
! The coding assumes we are using 6-th order centered finite difference for our
! derivatives. 
!
        if(llambda_effect) then
          boundary_value(:,:)=f(l1,:,:,j)/x(l1)-Lambda_V0*(f(l1,:,:,iuz)/x(l1)+Lambda_Omega)
        else
          boundary_value(:,:)=f(l1,:,:,j)/x(l1)
        endif
        f(l1-1,:,:,j)= f(l1+1,:,:,j) -  60.*boundary_value(:,:)*dx/45.
        f(l1-2,:,:,j)= f(l1+2,:,:,j) -  60.*boundary_value(:,:)*dx/9.
        f(l1-3,:,:,j)= f(l1+3,:,:,j) -  60.*boundary_value(:,:)*dx
!
! top boundary
!
      case('top')
        if(llambda_effect) then
          boundary_value(:,:)=f(l2,:,:,j)/x(l2)-Lambda_V0*f(l2,:,:,iuz)/x(l2)
        else
          boundary_value(:,:)=f(l2,:,:,j)/x(l2)  
        endif
        f(l2+1,:,:,j)= f(l2-1,:,:,j) +  60.*boundary_value(:,:)*dx/45.
        f(l2+2,:,:,j)= f(l2-2,:,:,j) +  60.*boundary_value(:,:)*dx/9.
        f(l2+3,:,:,j)= f(l2-3,:,:,j) +  60.*boundary_value(:,:)*dx

      case default
        call warning('bc_set_sfree_x',topbot//" should be `top' or `bot'")
!
      endselect
!
    endsubroutine bc_set_sfree_x
! **********************************************************************
    subroutine bc_set_jethat_x(f,jj,topbot,fracall,uzeroall)
!
! Sets tophat velocity profile at the inner (bot) boundary
!
!  3-jan-2008/dhruba: coded
!
      use Sub
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent(in) :: jj
      integer :: i,j,k
      real, dimension(mcom),intent(in) :: fracall,uzeroall
      real :: frac,uzero,ylim,ymid,ydif,y1,zlim,zmid,zdif,z1
      real :: yhat_min,yhat_max,zhat_min,zhat_max
      real :: width_hat=0.01
      real, dimension (ny) :: hatprofy
      real, dimension (nz) :: hatprofz
      y1 = xyz1(2)
      z1 = xyz1(3)
      frac = fracall(jj)
      uzero = uzeroall(jj)
!      write(*,*) frac,uzero,y0,z0,y1,z1
     if (lspherical_coords)then
!
        select case(topbot)
        case('bot')               ! bottom boundary
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
        case('top')               ! top boundary
          call warning('bc_set_jethat_x','Jet flowing out of the exit boundary ?')
          do i=1,nghost
            f(l2+i,:,:,j)=0.
          enddo
!
        case default
          call warning('bc_set_jethat_x',topbot//" should be `top' or `bot'")
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
! Stress-free boundary condition for spherical coordinate system. 
! d_{\theta}(A_{\phi}) = -A_{\phi}cot(\theta)/r  with A_{\theta} = 0 sets 
! B_{\theta}=0 in spherical polar
! coordinate system. This subroutine sets only the first part of this 
! boundary condition for 'j'-th component of f. 
!
!  25-Aug-2007/dhruba: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
      real :: cottheta
!
      select case(topbot)

      case('bot')               ! bottom boundary
!
! The coding assumes we are using 6-th order centered finite difference for our
! derivatives. 
!
        cottheta= cotth(m1)
        f(:,m1-1,:,j)= f(:,m1+1,:,j) +  60.*dy*cottheta*f(:,m1,:,j)/45.
        f(:,m1-2,:,j)= f(:,m1+2,:,j) +  60.*dy*cottheta*f(:,m1,:,j)/9.
        f(:,m1-3,:,j)= f(:,m1+3,:,j) +  60.*dy*cottheta*f(:,m1,:,j)
      case('top')               ! top boundary
        cottheta= cotth(m2)
        f(:,m2+1,:,j)= f(:,m2-1,:,j) -  60.*dy*cottheta*f(:,m2,:,j)/45.
        f(:,m2+2,:,j)= f(:,m2-2,:,j) -  60.*dy*cottheta*f(:,m2,:,j)/9.
        f(:,m2+3,:,j)= f(:,m2-3,:,j) -  60.*dy*cottheta*f(:,m2,:,j)
!
      case default
        call warning('bc_set_nfr_y',topbot//" should be `top' or `bot'")
!
      endselect
!
    endsubroutine bc_set_nfr_y
! **********************************************************************
    subroutine bc_set_sfree_y(f,topbot,j)
!
! Stress-free boundary condition for spherical coordinate system. 
! d_{\theta}(u_{\phi}) = u_{\phi}cot(\theta)/r  with u_{\theta} = 0 sets 
! S_{\theta \phi} component of the strain matrix to be zero in spherical 
! coordinate system. This subroutine sets only the first part of this 
! boundary condition for 'j'-th component of f. 
!
!  25-Aug-2007/dhruba: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
      real :: cottheta
!
      select case(topbot)

      case('bot')               ! bottom boundary
!
! The coding assumes we are using 6-th order centered finite difference for our
! derivatives. 
!
        cottheta= cotth(m1)
        f(:,m1-1,:,j)= f(:,m1+1,:,j) -  60.*dy*cottheta*f(:,m1,:,j)/45.
        f(:,m1-2,:,j)= f(:,m1+2,:,j) -  60.*dy*cottheta*f(:,m1,:,j)/9.
        f(:,m1-3,:,j)= f(:,m1+3,:,j) -  60.*dy*cottheta*f(:,m1,:,j)
      case('top')               ! top boundary
        cottheta= cotth(m2)
        f(:,m2+1,:,j)= f(:,m2-1,:,j) +  60.*dy*cottheta*f(:,m2,:,j)/45.
        f(:,m2+2,:,j)= f(:,m2-2,:,j) +  60.*dy*cottheta*f(:,m2,:,j)/9.
        f(:,m2+3,:,j)= f(:,m2-3,:,j) +  60.*dy*cottheta*f(:,m2,:,j)

      case default
        call warning('bc_set_sfree_y',topbot//" should be `top' or `bot'")
!
      endselect
!
    endsubroutine bc_set_sfree_y
! **********************************************************************
    subroutine bc_set_pfc_y(f,topbot,j)
!
! In spherical polar coordinate system,
! at a theta boundary set : $A_{r} = 0$ and $A_{\phi} = 0$,
! and demand $div A = 0$ gives the condition on $A_{\theta}$ to be
! $d/d{\theta}( A_{\theta}) + \cot(\theta)A_{\theta} = 0$ . 
! This subroutine sets this condition on 
! $j$ the component of f. As this is related to setting the
! perfect conducting boundary condition we call this "pfc". 
!
!  25-Aug-2007/dhruba: coded
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
      real :: cottheta
!
      select case(topbot)
!
      case('bot')               ! bottom boundary
!
! The coding assumes we are using 6-th order centered finite difference for our
! derivatives. 
!
        cottheta= cotth(m1)
        f(:,m1-1,:,j)= f(:,m1+1,:,j) +  2.*60.*dy*cottheta*f(:,m1,:,j)/45.
        f(:,m1-2,:,j)= f(:,m1+2,:,j) -  2.*60.*dy*cottheta*f(:,m1,:,j)/9.
        f(:,m1-3,:,j)= f(:,m1+3,:,j) +  2.*60.*dy*cottheta*f(:,m1,:,j)
      case('top')               ! top boundary
        cottheta= cotth(m2)
        f(:,m2+1,:,j)= f(:,m2-1,:,j) -  2.*60.*dy*cottheta*f(:,m2,:,j)/45.
        f(:,m2+2,:,j)= f(:,m2-2,:,j) +  2.*60.*dy*cottheta*f(:,m2,:,j)/9.
        f(:,m2+3,:,j)= f(:,m2-3,:,j) -  2.*60.*dy*cottheta*f(:,m2,:,j)

      case default
        call warning('bc_set_pfc_y',topbot//" should be `top' or `bot'")
!
      endselect
!
    endsubroutine bc_set_pfc_y
!***********************************************************************
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

      integer :: i

      select case(topbot)

      case('bot')               ! bottom boundary
        do i=1,nghost; f(:,m1-i,:,j) = f(:,m1+i,:,j) - 2*i*dy*val; enddo

      case('top')               ! top boundary
        do i=1,nghost; f(:,m2+i,:,j) = f(:,m2-i,:,j) + 2*i*dy*val; enddo

      case default
        call warning('bc_set_der_y',topbot//" should be `top' or `bot'")

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

      integer :: i

      select case(topbot)

      case('bot')               ! bottom boundary
        do i=1,nghost; f(:,:,n1-i,j) = f(:,:,n1+i,j) - 2*i*dz*val; enddo

      case('top')               ! top boundary
        do i=1,nghost; f(:,:,n2+i,j) = f(:,:,n2-i,j) + 2*i*dz*val; enddo

      case default
        call warning('bc_set_der_z',topbot//" should be `top' or `bot'")

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

      select case(topbot)

      case('bot')               ! bottom boundary
          do i=1,nghost
            f(l1-i,:,:,j)=((nghost+1-i)*f(l1,:,:,j))/(nghost+1)
          enddo

      case('top')               ! top boundary
          do i=1,nghost
            f(l2+i,:,:,j)=((nghost+1-i)*f(l2,:,:,j))/(nghost+1)
          enddo

      case default
        print*, "bc_van_x: ", topbot, " should be `top' or `bot'"

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

      select case(topbot)

      case('bot')               ! bottom boundary
          do i=1,nghost
            f(:,m1-i,:,j)=((nghost+1-i)*f(:,m1,:,j))/(nghost+1)
          enddo

      case('top')               ! top boundary
          do i=1,nghost
            f(:,m2+i,:,j)=((nghost+1-i)*f(:,m2,:,j))/(nghost+1)
          enddo

      case default
        print*, "bc_van_y: ", topbot, " should be `top' or `bot'"

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

      select case(topbot)

      case('bot')               ! bottom boundary
          do i=1,nghost
            f(:,:,n1-i,j)=((nghost+1-i)*f(:,:,n1,j))/(nghost+1)
          enddo

      case('top')               ! top boundary
          do i=1,nghost
            f(:,:,n2+i,j)=((nghost+1-i)*f(:,:,n2,j))/(nghost+1)
          enddo

      case default
        print*, "bc_van_z: ", topbot, " should be `top' or `bot'"

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
    real, dimension (mx,my) :: cpoly0,cpoly1,cpoly2
    integer :: i,j

    select case(topbot)

    case('bot')
      cpoly0(:,:)=f(:,:,n1,j)
      cpoly1(:,:)=-(3*f(:,:,n1,j)-4*f(:,:,n1+1,j)+f(:,:,n1+2,j))/(2*dz)
      cpoly2(:,:)=-(-f(:,:,n1,j)+2*f(:,:,n1+1,j)-f(:,:,n1+2,j)) /(2*dz**2)
      do i=1,nghost
        f(:,:,n1-i,j) = cpoly0(:,:) - cpoly1(:,:)*i*dz + cpoly2(:,:)*(i*dz)**2
      enddo

    case('top')
      cpoly0(:,:)=f(:,:,n2,j)
      cpoly1(:,:)=-(-3*f(:,:,n2,j)+4*f(:,:,n2-1,j)-f(:,:,n2-2,j))/(2*dz)
      cpoly2(:,:)=-(-f(:,:,n2,j)+2*f(:,:,n2-1,j)-f(:,:,n2-2,j))/(2*dz**2)
      do i=1,nghost
        f(:,:,n2+i,j) = cpoly0(:,:) + cpoly1(:,:)*i*dz + cpoly2(:,:)*(i*dz)**2
      enddo

    endselect

    endsubroutine bc_van3rd_z
!***********************************************************************
    subroutine bc_asym3(f,topbot,j)
!
!  Generalized antisymmetric bc (a al `a2') with removal of Nyquist wiggles
!  Does not seem to help against wiggles -- use upwinding instead
!
!  TEMPORARY HACK: Commented out calculation of Nyquist, as this creates
!  problems for some 2D runs and this boundary condition was not really
!  helpful so far. Will either have to find a better solution or remove
!  this altogether. wd, 21-jun-2003
!
!  17-jun-03/wolf: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: Nyquist=impossible
      integer :: j
!
      select case(topbot)

      case('bot')               ! bottom boundary
        ! Nyquist = 0.25*(f(:,:,n1,j)-2*f(:,:,n1+1,j)+f(:,:,n1+2,j))
        ! Nyquist = 0.0625*(     f(:,:,n1  ,j)+f(:,:,n1+4,j) &
        !                   - 4*(f(:,:,n1+1,j)+f(:,:,n1+3,j)) &
        !                   + 6* f(:,:,n1+2,j) )
        f(:,:,n1-1,j) = 2*f(:,:,n1,j) - f(:,:,n1+1,j) -4*Nyquist
        f(:,:,n1-2,j) = 2*f(:,:,n1,j) - f(:,:,n1+2,j)
        f(:,:,n1-3,j) = 2*f(:,:,n1,j) - f(:,:,n1+3,j) -4*Nyquist

      case('top')               ! top boundary
        ! Nyquist = 0.25*(f(:,:,n2,j)-2*f(:,:,n2-1,j)+f(:,:,n2-2,j))
        ! Nyquist = 0.0625*(     f(:,:,n2  ,j)+f(:,:,n2-4,j) &
        !                   - 4*(f(:,:,n2-1,j)+f(:,:,n2-3,j)) &
        !                   + 6* f(:,:,n2-2,j) )
        f(:,:,n2+1,j) = 2*f(:,:,n2,j) - f(:,:,n2-1,j) -4*Nyquist
        f(:,:,n2+2,j) = 2*f(:,:,n2,j) - f(:,:,n2-2,j)
        f(:,:,n2+3,j) = 2*f(:,:,n2,j) - f(:,:,n2-3,j) -4*Nyquist

      case default
        print*, "bc_asym3: ", topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_asym3
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
      integer :: i,j,k
!
      select case(topbot)

      case('bot')               ! bottom boundary
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

      case('top')               ! top boundary
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

      case default
        print*, "bc_onesided_x ", topbot, " should be `top' or `bot'"

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
      select case(topbot)

      case('bot')               ! bottom boundary
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

      case('top')               ! top boundary
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

      case default
        print*, "bc_onesided_x_old ", topbot, " should be `top' or `bot'"

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
      integer :: i,j,k
!
      select case(topbot)

      case('bot')               ! bottom boundary
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
                               
      case('top')               ! top boundary
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

      case default
        print*, "bc_onesided_7 ", topbot, " should be `top' or `bot'"

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
!   5-apr-03/axel: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j,k
!
      select case(topbot)

      case('bot')               ! bottom boundary
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

      case('top')               ! top boundary
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

      case default
        print*, "bc_onesided_z ", topbot, " should be `top' or `bot'"

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
!   5-apr-03/axel: coded
!  10-mar-09/axel: corrected
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j,k
!
      select case(topbot)

      case('bot')               ! bottom boundary
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
                               
      case('top')               ! top boundary
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

      case default
        print*, "bc_onesided_z ", topbot, " should be `top' or `bot'"

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
      select case(topbot)

      case('bot')               ! bottom boundary
        f(:,:,n1-1,j)=0.25*(  9*f(:,:,n1,j)- 3*f(:,:,n1+1,j)- 5*f(:,:,n1+2,j)+ 3*f(:,:,n1+3,j))
        f(:,:,n1-2,j)=0.05*( 81*f(:,:,n1,j)-43*f(:,:,n1+1,j)-57*f(:,:,n1+2,j)+39*f(:,:,n1+3,j))
        f(:,:,n1-3,j)=0.05*(127*f(:,:,n1,j)-81*f(:,:,n1+1,j)-99*f(:,:,n1+2,j)+73*f(:,:,n1+3,j))

      case('top')               ! top boundary
        f(:,:,n2+1,j)=0.25*(  9*f(:,:,n2,j)- 3*f(:,:,n2-1,j)- 5*f(:,:,n2-2,j)+ 3*f(:,:,n2-3,j))
        f(:,:,n2+2,j)=0.05*( 81*f(:,:,n2,j)-43*f(:,:,n2-1,j)-57*f(:,:,n2-2,j)+39*f(:,:,n2-3,j))
        f(:,:,n2+3,j)=0.05*(127*f(:,:,n2,j)-81*f(:,:,n2-1,j)-99*f(:,:,n2-2,j)+73*f(:,:,n2-3,j))

      case default
        print*, "bc_extrap_2_1: ", topbot, " should be `top' or `bot'"

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
      select case(topbot)

      case('bot')               ! bottom boundary
        f(l1-1,:,:,j)=0.25*(  9*f(l1,:,:,j)- 3*f(l1+1,:,:,j)- 5*f(l1+2,:,:,j)+ 3*f(l1+3,:,:,j))
        f(l1-2,:,:,j)=0.05*( 81*f(l1,:,:,j)-43*f(l1+1,:,:,j)-57*f(l1+2,:,:,j)+39*f(l1+3,:,:,j))
        f(l1-3,:,:,j)=0.05*(127*f(l1,:,:,j)-81*f(l1+1,:,:,j)-99*f(l1+2,:,:,j)+73*f(l1+3,:,:,j))

      case('top')               ! top boundary
        f(l2+1,:,:,j)=0.25*(  9*f(l2,:,:,j)- 3*f(l2-1,:,:,j)- 5*f(l2-2,:,:,j)+ 3*f(l2-3,:,:,j))
        f(l2+2,:,:,j)=0.05*( 81*f(l2,:,:,j)-43*f(l2-1,:,:,j)-57*f(l2-2,:,:,j)+39*f(l2-3,:,:,j))
        f(l2+3,:,:,j)=0.05*(127*f(l2,:,:,j)-81*f(l2-1,:,:,j)-99*f(l2-2,:,:,j)+73*f(l2-3,:,:,j))

      case default
        print*, "bcx_extrap_2_1: ", topbot, " should be `top' or `bot'"

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
      select case(topbot)

      case('bot')               ! bottom boundary
        f(:,m1-1,:,j)=0.25*(  9*f(:,m1,:,j)- 3*f(:,m1+1,:,j)- 5*f(:,m1+2,:,j)+ 3*f(:,m1+3,:,j))
        f(:,m1-2,:,j)=0.05*( 81*f(:,m1,:,j)-43*f(:,m1+1,:,j)-57*f(:,m1+2,:,j)+39*f(:,m1+3,:,j))
        f(:,m1-3,:,j)=0.05*(127*f(:,m1,:,j)-81*f(:,m1+1,:,j)-99*f(:,m1+2,:,j)+73*f(:,m1+3,:,j))

      case('top')               ! top boundary
        f(:,m2+1,:,j)=0.25*(  9*f(:,m2,:,j)- 3*f(:,m2-1,:,j)- 5*f(:,m2-2,:,j)+ 3*f(:,m2-3,:,j))
        f(:,m2+2,:,j)=0.05*( 81*f(:,m2,:,j)-43*f(:,m2-1,:,j)-57*f(:,m2-2,:,j)+39*f(:,m2-3,:,j))
        f(:,m2+3,:,j)=0.05*(127*f(:,m2,:,j)-81*f(:,m2-1,:,j)-99*f(:,m2-2,:,j)+73*f(:,m2-3,:,j))

      case default
        print*, "bcy_extrap_2_1: ", topbot, " should be `top' or `bot'"

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
      select case(topbot)

      case('bot')               ! bottom boundary
        f(:,:,n1-1,j)=0.2   *(  9*f(:,:,n1,j)                 -  4*f(:,:,n1+2,j)- 3*f(:,:,n1+3,j)+ 3*f(:,:,n1p4,j))
        f(:,:,n1-2,j)=0.2   *( 15*f(:,:,n1,j)- 2*f(:,:,n1+1,j)-  9*f(:,:,n1+2,j)- 6*f(:,:,n1+3,j)+ 7*f(:,:,n1p4,j))
        f(:,:,n1-3,j)=1./35.*(157*f(:,:,n1,j)-33*f(:,:,n1+1,j)-108*f(:,:,n1+2,j)-68*f(:,:,n1+3,j)+87*f(:,:,n1p4,j))

      case('top')               ! top boundary
        f(:,:,n2+1,j)=0.2   *(  9*f(:,:,n2,j)                 -  4*f(:,:,n2-2,j)- 3*f(:,:,n2-3,j)+ 3*f(:,:,n2m4,j))
        f(:,:,n2+2,j)=0.2   *( 15*f(:,:,n2,j)- 2*f(:,:,n2-1,j)-  9*f(:,:,n2-2,j)- 6*f(:,:,n2-3,j)+ 7*f(:,:,n2m4,j))
        f(:,:,n2+3,j)=1./35.*(157*f(:,:,n2,j)-33*f(:,:,n2-1,j)-108*f(:,:,n2-2,j)-68*f(:,:,n2-3,j)+87*f(:,:,n2m4,j))

      case default
        print*, "bc_extrap_2_2: ", topbot, " should be `top' or `bot'"

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
      select case(topbot)

      case('bot')               ! bottom boundary
        f(l1-1,:,:,j)=0.2   *(  9*f(l1,:,:,j)                 -  4*f(l1+2,:,:,j)- 3*f(l1+3,:,:,j)+ 3*f(l1p4,:,:,j))
        f(l1-2,:,:,j)=0.2   *( 15*f(l1,:,:,j)- 2*f(l1+1,:,:,j)-  9*f(l1+2,:,:,j)- 6*f(l1+3,:,:,j)+ 7*f(l1p4,:,:,j))
        f(l1-3,:,:,j)=1./35.*(157*f(l1,:,:,j)-33*f(l1+1,:,:,j)-108*f(l1+2,:,:,j)-68*f(l1+3,:,:,j)+87*f(l1p4,:,:,j))

      case('top')               ! top boundary
        f(l2+1,:,:,j)=0.2   *(  9*f(l2,:,:,j)                 -  4*f(l2-2,:,:,j)- 3*f(l2-3,:,:,j)+ 3*f(l2m4,:,:,j))
        f(l2+2,:,:,j)=0.2   *( 15*f(l2,:,:,j)- 2*f(l2-1,:,:,j)-  9*f(l2-2,:,:,j)- 6*f(l2-3,:,:,j)+ 7*f(l2m4,:,:,j))
        f(l2+3,:,:,j)=1./35.*(157*f(l2,:,:,j)-33*f(l2-1,:,:,j)-108*f(l2-2,:,:,j)-68*f(l2-3,:,:,j)+87*f(l2m4,:,:,j))

      case default
        print*, "bcx_extrap_2_2: ", topbot, " should be `top' or `bot'"

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
      select case(topbot)

      case('bot')               ! bottom boundary
        f(:,m1-1,:,j)=0.2   *(  9*f(:,m1,:,j)                 -  4*f(:,m1+2,:,j)- 3*f(:,m1+3,:,j)+ 3*f(:,m1p4,:,j))
        f(:,m1-2,:,j)=0.2   *( 15*f(:,m1,:,j)- 2*f(:,m1+1,:,j)-  9*f(:,m1+2,:,j)- 6*f(:,m1+3,:,j)+ 7*f(:,m1p4,:,j))
        f(:,m1-3,:,j)=1./35.*(157*f(:,m1,:,j)-33*f(:,m1+1,:,j)-108*f(:,m1+2,:,j)-68*f(:,m1+3,:,j)+87*f(:,m1p4,:,j))

      case('top')               ! top boundary
        f(:,m2+1,:,j)=0.2   *(  9*f(:,m2,:,j)                 -  4*f(:,m2-2,:,j)- 3*f(:,m2-3,:,j)+ 3*f(:,m2m4,:,j))
        f(:,m2+2,:,j)=0.2   *( 15*f(:,m2,:,j)- 2*f(:,m2-1,:,j)-  9*f(:,m2-2,:,j)- 6*f(:,m2-3,:,j)+ 7*f(:,m2m4,:,j))
        f(:,m2+3,:,j)=1./35.*(157*f(:,m2,:,j)-33*f(:,m2-1,:,j)-108*f(:,m2-2,:,j)-68*f(:,m2-3,:,j)+87*f(:,m2m4,:,j))

      case default
        print*, "bcy_extrap_2_2: ", topbot, " should be `top' or `bot'"

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
      select case(topbot)
!
      case('bot')               ! bottom boundary
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
      case('top')               ! top boundary
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
      select case(topbot)

!       case('bot')               ! bottom boundary
!         f(:,:,n1  ,j)= 0.       ! set bdry value=0 (indep of initcond)
!         f(:,:,n1-1,j)=- 3*f(:,:,n1+1,j)+  f(:,:,n1+2,j)
!         f(:,:,n1-2,j)=- 8*f(:,:,n1+1,j)+3*f(:,:,n1+2,j)
!         f(:,:,n1-3,j)=-15*f(:,:,n1+1,j)+6*f(:,:,n1+2,j)

!       case('top')               ! top boundary
!         f(:,:,n2  ,j)= 0.       ! set bdry value=0 (indep of initcond)
!         f(:,:,n2+1,j)=- 3*f(:,:,n2-1,j)+  f(:,:,n2-2,j)
!         f(:,:,n2+2,j)=- 8*f(:,:,n2-1,j)+3*f(:,:,n2-2,j)
!         f(:,:,n2+3,j)=-15*f(:,:,n2-1,j)+6*f(:,:,n2-2,j)

!! Nyquist-filtering
      case('bot')               ! bottom boundary
        f(:,:,n1  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n1-1,j)=(1/11.)*(-17*f(:,:,n1+1,j)- 9*f(:,:,n1+2,j)+ 8*f(:,:,n1+3,j))
        f(:,:,n1-2,j)=      2*(- 2*f(:,:,n1+1,j)-   f(:,:,n1+2,j)+   f(:,:,n1+3,j))
        f(:,:,n1-3,j)=(3/11.)*(-27*f(:,:,n1+1,j)-13*f(:,:,n1+2,j)+14*f(:,:,n1+3,j))

      case('top')               ! top boundary
        f(:,:,n2  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n2+1,j)=(1/11.)*(-17*f(:,:,n2-1,j)- 9*f(:,:,n2-2,j)+ 8*f(:,:,n2-3,j))
        f(:,:,n2+2,j)=      2*(- 2*f(:,:,n2-1,j)-   f(:,:,n2-2,j)+   f(:,:,n2-3,j))
        f(:,:,n2+3,j)=(3/11.)*(-27*f(:,:,n2-1,j)-13*f(:,:,n2-2,j)+14*f(:,:,n2-3,j))

! !! Nyquist-transparent
!       case('bot')               ! bottom boundary
!         f(:,:,n1  ,j)=0.        ! set bdry value=0 (indep of initcond)
!         f(:,:,n1-1,j)=(1/11.)*(-13*f(:,:,n1+1,j)-14*f(:,:,n1+2,j)+10*f(:,:,n1+3,j))
!         f(:,:,n1-2,j)=(1/11.)*(-48*f(:,:,n1+1,j)-17*f(:,:,n1+2,j)+20*f(:,:,n1+3,j))
!         f(:,:,n1-3,j)=         - 7*f(:,:,n1+1,j)- 4*f(:,:,n1+2,j)+ 4*f(:,:,n1+3,j)

!       case('top')               ! top boundary
!         f(:,:,n2  ,j)=0.        ! set bdry value=0 (indep of initcond)
!         f(:,:,n2+1,j)=(1/11.)*(-13*f(:,:,n2-1,j)-14*f(:,:,n2-2,j)+10*f(:,:,n2-3,j))
!         f(:,:,n2+2,j)=(1/11.)*(-48*f(:,:,n2-1,j)-17*f(:,:,n2-2,j)+20*f(:,:,n2-3,j))
!         f(:,:,n2+3,j)=         - 7*f(:,:,n2-1,j)- 4*f(:,:,n2-2,j)+ 4*f(:,:,n2-3,j)

      case default
        print*, "bc_extrap0_2_0: ", topbot, " should be `top' or `bot'"

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
      select case(topbot)

      case('bot')               ! bottom boundary
        f(:,:,n1  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n1-1,j)=0.25*(- 3*f(:,:,n1+1,j)- 5*f(:,:,n1+2,j)+ 3*f(:,:,n1+3,j))
        f(:,:,n1-2,j)=0.05*(-43*f(:,:,n1+1,j)-57*f(:,:,n1+2,j)+39*f(:,:,n1+3,j))
        f(:,:,n1-3,j)=0.05*(-81*f(:,:,n1+1,j)-99*f(:,:,n1+2,j)+73*f(:,:,n1+3,j))

      case('top')               ! top boundary
        f(:,:,n2  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n2+1,j)=0.25*(- 3*f(:,:,n2-1,j)- 5*f(:,:,n2-2,j)+ 3*f(:,:,n2-3,j))
        f(:,:,n2+2,j)=0.05*(-43*f(:,:,n2-1,j)-57*f(:,:,n2-2,j)+39*f(:,:,n2-3,j))
        f(:,:,n2+3,j)=0.05*(-81*f(:,:,n2-1,j)-99*f(:,:,n2-2,j)+73*f(:,:,n2-3,j))

      case default
        print*, "bc_extrap0_2_1: ", topbot, " should be `top' or `bot'"

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
      select case(topbot)

      case('bot')               ! bottom boundary
        f(:,:,n1  ,j)= 0.       ! set bdry value=0 (indep of initcond)
        f(:,:,n1-1,j)=0.2   *(                 -  4*f(:,:,n1+2,j)- 3*f(:,:,n1+3,j)+ 3*f(:,:,n1p4,j))
        f(:,:,n1-2,j)=0.2   *(- 2*f(:,:,n1+1,j)-  9*f(:,:,n1+2,j)- 6*f(:,:,n1+3,j)+ 7*f(:,:,n1p4,j))
        f(:,:,n1-3,j)=1./35.*(-33*f(:,:,n1+1,j)-108*f(:,:,n1+2,j)-68*f(:,:,n1+3,j)+87*f(:,:,n1p4,j))

      case('top')               ! top boundary
        f(:,:,n2  ,j)= 0.       ! set bdry value=0 (indep of initcond)
        f(:,:,n2+1,j)=0.2   *(                 -  4*f(:,:,n2-2,j)- 3*f(:,:,n2-3,j)+ 3*f(:,:,n2m4,j))
        f(:,:,n2+2,j)=0.2   *(- 2*f(:,:,n2-1,j)-  9*f(:,:,n2-2,j)- 6*f(:,:,n2-3,j)+ 7*f(:,:,n2m4,j))
        f(:,:,n2+3,j)=1./35.*(-33*f(:,:,n2-1,j)-108*f(:,:,n2-2,j)-68*f(:,:,n2-3,j)+87*f(:,:,n2m4,j))

      case default
        print*, "bc_extrap0_2_2: ", topbot, " should be `top' or `bot'"

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
      select case(topbot)
!
      case('bot')               ! bottom boundary
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
      case('top')               ! top boundary
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
      real, dimension (mx,my) :: fder
      integer :: i,j
!
      select case(topbot)
!
! Bottom boundary
!
      case('bot')
        do i=1,nghost
          fder=(-3*f(:,:,n1-i+1,j)+4*f(:,:,n1-i+2,j)&
               -f(:,:,n1-i+3,j))/(2*dz)
          f(:,:,n1-i,j)=f(:,:,n1-i+2,j)-2*dz*fder
        enddo
      case('top')
        do i=1,nghost
          fder=(3*f(:,:,n2+i-1,j)-4*f(:,:,n2+i-2,j)&
               +f(:,:,n2+i-3,j))/(2*dz)
          f(:,:,n2+i,j)=f(:,:,n2+i-2,j)+2*dz*fder
        enddo
      case default
        print*,"bc_db_z: invalid argument for 'bc_db_z'"
      endselect
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
      real, dimension (my,mz) :: fder
      integer :: i,j
!
      select case(topbot)
!
! Bottom boundary
!
      case('bot')
        do i=1,nghost
          fder=(-3*f(l1-i+1,:,:,j)+4*f(l1-i+2,:,:,j)&
               -f(l1-i+3,:,:,j))/(2*dx)
          f(l1-i,:,:,j)=f(l1-i+2,:,:,j)-2*dx*fder
        enddo
      case('top')
        do i=1,nghost
          fder=(3*f(l2+i-1,:,:,j)-4*f(l2+i-2,:,:,j)&
               +f(l2+i-3,:,:,j))/(2*dx)
          f(l2+i,:,:,j)=f(l2+i-2,:,:,j)+2*dx*fder
        enddo
      case default
        print*,"bc_db_x: invalid argument for 'bc_db_x'"
      endselect
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
      use EquationOfState, only: gamma1, cs2top, cs2bot
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: sgn,i,j
!
      select case(topbot)
!
!  lower boundary
!
      case('bot')
         select case (force_lower_bound)
         case ('uxy_sin-cos')
            call bc_force_uxy_sin_cos(f,n1,j)
         case ('axy_sin-cos')
            call bc_force_axy_sin_cos(f,n1,j)
         case ('uxy_convection')
            call uu_driver(f)
         !case ('kepler')
         !   call bc_force_kepler(f,n1,j)
         case ('cT')
            f(:,:,n1,j) = log(cs2bot/gamma1)
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
      case('top')
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
            f(:,:,n2,j) = log(cs2top/gamma1)
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
    subroutine bc_force_x(f,sgn,topbot,j)
!
!  Force values of j-th variable on x-boundaries topbot.
!
!  09-mar-2007/dintrans: coded
!
      use EquationOfState, only: gamma1, cs2top, cs2bot
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: sgn,i,j
!
      select case(topbot)
!
!  lower boundary
!
      case('bot')
         select case (force_lower_bound)
         case ('cT')
            f(l1,:,:,ilnTT) = log(cs2bot/gamma1)
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
      case('top')
         select case (force_upper_bound)
         case ('cT')
            f(l2,:,:,ilnTT) = log(cs2top/gamma1)
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
      select case(topbot)

      case('bot')               ! bottom boundary
          f(1:l1-1,:,:,j)=1.

      case('top')               ! top boundary
          f(l2+1:mx,:,:,j)=1.

      case default
        print*, "bc_one_x: ",topbot, " should be `top' or `bot'"

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
      select case(topbot)

      case('bot')               ! bottom boundary
          f(:,1:m1-1,:,j)=1.

      case('top')               ! top boundary
          f(:,m2+1:my,:,j)=1.

      case default
        print*, "bc_one_y: ", topbot, " should be `top' or `bot'"

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
      select case(topbot)

      case('bot')               ! bottom boundary
          f(:,:,1:n1-1,j)=1.

      case('top')               ! top boundary
          f(:,:,n2+1:mz,j)=1.

      case default
        print*, "bc_one_z: ", topbot, " should be `top' or `bot'"

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

      select case(topbot)
      case('bot')               ! bottom boundary
        lfrozen_bot_var_x(j) = .true.
      case('top')               ! top boundary
        lfrozen_top_var_x(j) = .true.
      case default
        print*, "bc_freeze_var_x: ", topbot, " should be `top' or `bot'"
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

      select case(topbot)
      case('bot')               ! bottom boundary
        lfrozen_bot_var_y(j) = .true.
      case('top')               ! top boundary
        lfrozen_top_var_y(j) = .true.
      case default
        print*, "bc_freeze_var_y: ", topbot, " should be `top' or `bot'"
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

      select case(topbot)
      case('bot')               ! bottom boundary
        lfrozen_bot_var_z(j) = .true.
      case('top')               ! top boundary
        lfrozen_top_var_z(j) = .true.
      case default
        print*, "bc_freeze_var_z: ", topbot, " should be `top' or `bot'"
      endselect
!
    endsubroutine bc_freeze_var_z
!***********************************************************************
    subroutine update_ghosts(a)
!
!  update all ghost zones of a
!  21-sep-02/wolf: extracted from wsnaps
!
      use Cparam
      use Mpicomm
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
     subroutine uu_driver(f)
!
!  Simulated velocity field used as photospherec motions
!  Use of velocity field produced by Boris Gudiksen
!
!  27-mai-04/bing: coded
!  11-aug-06/axel: make it compile with nprocx>0, renamed quenching -> quen
!  18-jun-08/bing: quenching depends on B^2, not only Bz^2
!
       use EquationOfState, only : gamma,gamma1,gamma11,cs20,lnrho0

       real, dimension (mx,my,mz,mfarray) :: f
       real, dimension (nxgrid,nygrid),save :: uxl,uxr,uyl,uyr
       real, dimension (nxgrid,nygrid) :: uxd,uyd
       real, dimension (nx,ny) :: quen,pp,betaq,fac
       real, dimension (nx,ny) :: bbx,bby,bbz,bb2
       integer :: lend,iostat=0,i=0,j
       real,save :: tl=0.,tr=0.,delta_t=0.

       intent (inout) :: f
!
!     Read the time table
!
       if (t*unit_time < tl+delta_t .or. t*unit_time>=tr+delta_t .and. iostat /= -2) then
!
          inquire(IOLENGTH=lend) tl
          close (10)
          open (10,file='driver/time_k',form='unformatted',status='unknown',recl=lend,access='direct')
!
          iostat = 0
          i=0
          do while (iostat == 0)
            i=i+1
            read (10,rec=i,iostat=iostat) tl
            read (10,rec=i+1,iostat=iostat) tr
            if (iostat /= 0) then
              i=1
              delta_t = t*unit_time                  ! EOF is reached => read again
              read (10,rec=i,iostat=iostat) tl
              read (10,rec=i+1,iostat=iostat) tr
              iostat=-1
            else
              if (t*unit_time>=tl+delta_t .and. t*unit_time<tr+delta_t)  iostat=-1 ! correct time step is reached
            endif
          enddo
          close (10)
!
! Read velocity field
!
          open (10,file='driver/vel_k.dat',form='unformatted',status='unknown',recl=lend*nxgrid*nygrid,access='direct')
          read (10,rec=(2*i-1)) uxl
          read (10,rec=2*i)     uyl

          read (10,rec=2*i+1)   uxr
          read (10,rec=2*i+2)   uyr
          close (10)

          uxl = uxl / 10. / unit_velocity
          uxr = uxr / 10. / unit_velocity
          uyl = uyl / 10. / unit_velocity
          uyr = uyr / 10. / unit_velocity

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
!   Calculate B^2 for plasma beta
!
!----------------------------------------------------------------------------------------
       if (nygrid/=1) then
          fac=(1./60)*spread(dy_1(m1:m2),1,nx)
          bbx= fac*(+ 45.0*(f(l1:l2,m1+1:m2+1,n1,iaz)-f(l1:l2,m1-1:m2-1,n1,iaz)) &
                -  9.0*(f(l1:l2,m1+2:m2+2,n1,iaz)-f(l1:l2,m1-2:m2-2,n1,iaz)) &
                +      (f(l1:l2,m1+3:m2+3,n1,iaz)-f(l1:l2,m1-3:m2-3,n1,iaz)))
       else
          if (ip<=5) print*, 'uu_driver: Degenerate case in y-direction'
       endif
       if (nzgrid/=1) then
          fac=(1./60)*spread(spread(dz_1(n1),1,nx),2,ny)
          bbx= bbx -fac*(+ 45.0*(f(l1:l2,m1:m2,n1+1,iay)-f(l1:l2,m1:m2,n1-1,iay)) &
               -  9.0*(f(l1:l2,m1:m2,n1+2,iay)-f(l1:l2,m1:m2,n1-2,iay)) &
               +      (f(l1:l2,m1:m2,n1+3,iay)-f(l1:l2,m1:m2,n1-2,iay)))
       else
          if (ip<=5) print*, 'uu_driver: Degenerate case in z-direction'
       endif
!----------------------------------------------------------------------------------------
       if (nzgrid/=1) then
          fac=(1./60)*spread(spread(dz_1(n1),1,nx),2,ny)
          bby= fac*(+ 45.0*(f(l1:l2,m1:m2,n1+1,iax)-f(l1:l2,m1:m2,n1-1,iax)) &
               -  9.0*(f(l1:l2,m1:m2,n1+2,iax)-f(l1:l2,m1:m2,n1-2,iax)) &
               +      (f(l1:l2,m1:m2,n1+3,iax)-f(l1:l2,m1:m2,n1-3,iax)))
       else
          if (ip<=5) print*, 'uu_driver: Degenerate case in z-direction'
       endif
       if (nxgrid/=1) then
          fac=(1./60)*spread(dx_1(l1:l2),2,ny)
          bby= bby -fac*(+ 45.0*(f(l1+1:l2+1,m1:m2,n1,iaz)-f(l1-1:l2-1,m1:m2,n1,iaz)) &
               -  9.0*(f(l1+2:l2+2,m1:m2,n1,iaz)-f(l1-2:l2-2,m1:m2,n1,iaz)) &
               +      (f(l1+3:l2+3,m1:m2,n1,iaz)-f(l1-3:l2-3,m1:m2,n1,iaz)))
       else
          if (ip<=5) print*, 'uu_driver: Degenerate case in x-direction'
       endif
!----------------------------------------------------------------------------------------
       if (nxgrid/=1) then
          fac=(1./60)*spread(dx_1(l1:l2),2,ny)
          bbz= fac*(+ 45.0*(f(l1+1:l2+1,m1:m2,n1,iay)-f(l1-1:l2-1,m1:m2,n1,iay)) &
               -  9.0*(f(l1+2:l2+2,m1:m2,n1,iay)-f(l1-2:l2-2,m1:m2,n1,iay)) &
               +      (f(l1+3:l2+3,m1:m2,n1,iay)-f(l1-3:l2-3,m1:m2,n1,iay)))
       else
          if (ip<=5) print*, 'uu_driver: Degenerate case in x-direction'
       endif
       if (nygrid/=1) then
          fac=(1./60)*spread(dy_1(m1:m2),1,nx)
          bbz= bbz -fac*(+ 45.0*(f(l1:l2,m1+1:m2+1,n1,iax)-f(l1:l2,m1-1:m2-1,n1,iax)) &
               -  9.0*(f(l1:l2,m1+2:m2+2,n1,iax)-f(l1:l2,m1-2:m2-2,n1,iax)) &
               +      (f(l1:l2,m1+3:m2+3,n1,iax)-f(l1:l2,m1-3:m2-3,n1,iax)))
       else
          if (ip<=5) print*, 'uu_driver: Degenerate case in y-direction'
       endif
!----------------------------------------------------------------------------------------
!
       bb2 = bbx*bbx + bby*bby + bbz*bbz
       bb2 = bb2/(2.*mu0)
!
       if (ltemperature) then
          pp = gamma1*gamma11*exp(f(l1:l2,m1:m2,n1,ilnrho)+f(l1:l2,m1:m2,n1,ilnTT))
       else if (lentropy) then          
          if (pretend_lnTT) then
             pp = gamma1*gamma11*exp(f(l1:l2,m1:m2,n1,ilnrho)+f(l1:l2,m1:m2,n1,iss))
          else
             pp = gamma* (f(l1:l2,m1:m2,n1,iss)+f(l1:l2,m1:m2,n1,ilnrho))-gamma1*lnrho0
             pp = exp(pp) * cs20*gamma11
          endif
       else
          pp=gamma11*cs20*exp(lnrho0)
       endif
!
!   limit plasma beta
!
       betaq = pp / max(tini,bb2)
!
       quen=(1.+betaq**2)/(1e3+betaq**2)
!
!   Fill bottom layer with velocity field
!
       f(l1:l2,m1:m2,n1,iux)=uxd(ipx*nx+1:(ipx+1)*nx,ipy*ny+1:(ipy+1)*ny)*quen
       f(l1:l2,m1:m2,n1,iuy)=uyd(ipx*nx+1:(ipx+1)*nx,ipy*ny+1:(ipy+1)*ny)*quen
!
     endsubroutine uu_driver
!***********************************************************************
    subroutine bc_lnTT_flux_x(f,topbot)
!
!  constant flux boundary condition for temperature (called when bcx='c1')
!  12-Mar-2007/dintrans: coded
!
      use SharedVariables, only: get_shared_variable
!
      real, pointer :: hcond0, hcond1, Fbot
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (my,mz) :: tmp_yz
      integer :: i,ierr
!
!  Do the `c1' boundary condition (constant heat flux) for lnTT.
!  check whether we want to do top or bottom (this is precessor dependent)
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

      select case(topbot)
!
!  bottom boundary
!  ===============
!
      case('bot')
        tmp_yz=-Fbot/(hcond0*hcond1)/exp(f(l1,:,:,ilnTT))
!
!  enforce dlnT/dx = - Fbot/(K*T)
!
        do i=1,nghost
          f(l1-i,:,:,ilnTT)=f(l1+i,:,:,ilnTT)-2*i*dx*tmp_yz
        enddo

      case default
        call fatal_error('bc_lnTT_flux_x','invalid argument')

      endselect
!
    endsubroutine bc_lnTT_flux_x
!***********************************************************************
    subroutine bc_lnTT_flux_z(f,topbot)
!
!  constant flux boundary condition for temperature (called when bcz='c1')
!  12-May-07/dintrans: coded
!
      use SharedVariables, only: get_shared_variable
!
      real, pointer :: hcond0, Fbot
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: tmp_xy
      integer :: i,ierr
!
!  Do the `c1' boundary condition (constant heat flux) for lnTT or TT (if
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

      select case(topbot)
      case('bot')
        if (ltemperature_nolog) then
          tmp_xy=-Fbot/hcond0
        else
          tmp_xy=-Fbot/hcond0/exp(f(:,:,n1,ilnTT))
        endif
        do i=1,nghost
          f(:,:,n1-i,ilnTT)=f(:,:,n1+i,ilnTT)-2.*i*dz*tmp_xy
        enddo

      case default
        call fatal_error('bc_lnTT_flux_z','invalid argument')

      endselect
!
    endsubroutine bc_lnTT_flux_z
!***********************************************************************
    subroutine bc_ss_flux_x(f,topbot)
!
!  constant flux boundary condition for entropy (called when bcx='c1')
!  17-mar-07/dintrans: coded
!
      use EquationOfState, only: gamma, gamma1, lnrho0, cs20
      use SharedVariables, only: get_shared_variable
!
      real, pointer :: FbotKbot
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (my,mz) :: tmp_yz,cs2_yz
      integer :: i,ierr
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
      call get_shared_variable('FbotKbot',FbotKbot,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux_x: "//&
           "there was a problem when getting FbotKbot")
!
      select case(topbot)
!
!  bottom boundary
!  ===============
!
      case('bot')
        if (headtt) print*,'bc_ss_flux_x: FbotKbot=',FbotKbot
!
!  calculate Fbot/(K*cs2)
!
!       cs2_yz=cs20*exp(gamma1*(f(l1,:,:,ilnrho)-lnrho0)+cv1*f(l1,:,:,iss))
        cs2_yz=cs20*exp(gamma1*(f(l1,:,:,ilnrho)-lnrho0)+gamma*f(l1,:,:,iss))
        tmp_yz=FbotKbot/cs2_yz
!
!  enforce ds/dx + gamma1/gamma*dlnrho/dx = - gamma1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
!         f(l1-i,:,:,iss)=f(l1+i,:,:,iss)+(cp-cv)* &
          f(l1-i,:,:,iss)=f(l1+i,:,:,iss)+gamma1/gamma* &
              (f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho)+2*i*dx*tmp_yz)
        enddo
!
!  top boundary
!  ============
!
      case('top')
        call fatal_error('bc_ss_flux_x','not implemented for top')

      case default
        call fatal_error('bc_ss_flux_x','invalid argument')

      endselect
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

      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
      integer, intent (in) :: j

      real, dimension (nx,ny) :: kx,ky,kappa,exp_fact
      real, dimension (nx,ny) :: tmp_re,tmp_im
      integer :: i
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
!  Check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  Potential field condition at the bottom
!
      case('bot')

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

        enddo
!
!  Potential field condition at the top
!
      case('top')

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

        enddo

      case default

        if (lroot) print*,"bc_del2zero: invalid argument"

      endselect

    endsubroutine bc_del2zero
!***********************************************************************
    subroutine bc_zero_z(f,topbot,j)
!
!  Zero value in the ghost zones.
!
!  13-aug-2007/anders: implemented
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      select case(topbot)
!
!  Bottom boundary.
!
      case('bot')
        f(:,:,1:n1-1,j)=0.0
!
!  Top boundary.
!
      case('top')
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
      select case(topbot)
!
!  Bottom boundary.
!
      case('bot')
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
      case('top')
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
      select case(topbot)
!
!  Bottom boundary.
!
      case('bot')
        do i=1,nghost; f(:,:,n1-i,j)=f(:,:,n1,j); enddo
!
!  Top boundary.
!
      case('top')
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
      select case(topbot)
      case('bot')               ! bottom boundary
        lfrozen_bb_bot(j-iax+1) = .true.    ! set flag
      case('top')               ! top boundary
        lfrozen_bb_top(j-iax+1) = .true.    ! set flag
      case default
        print*, "bc_frozen_in_bb: ", topbot, " should be `top' or `bot'"
      endselect
!
      lfirstcall=.false.
!
    endsubroutine bc_frozen_in_bb
!***********************************************************************
    subroutine bc_aa_pot3(f,topbot)
!
!  Potential field boundary condition
!
!  11-oct-06/wolf: Adapted from Tobi's bc_aa_pot2
!
      use Fourier, only: fourier_transform_xy_xy
      use Mpicomm, only: communicate_bc_aa_pot
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
!
      real, dimension (nx,ny,iax:iaz) :: aa_re,aa_im
      real, dimension (nx,ny) :: kx,ky,kappa,kappa1,exp_fact
      real, dimension (nx,ny) :: tmp_re,tmp_im
      real    :: delta_z
      integer :: i,j
!
!  Get local wave numbers
!
      kx = spread(kx_fft(ipx*nx+1:ipx*nx+nx),2,ny)
      ky = spread(ky_fft(ipy*ny+1:ipy*ny+ny),1,nx)
!
!  Calculate 1/k^2, zero mean
!
      kappa = sqrt(kx**2 + ky**2)
      where (kappa > 0)
        kappa1 = 1/kappa
      elsewhere
        kappa1 = 0
      endwhere
!
!  Check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  Potential field condition at the bottom
!
      case('bot')

        do j=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          delta_z  = z(n1+j) - z(n1-j)
          exp_fact = exp(-kappa*delta_z)
!
!  Determine potential field in ghost zones
!
          !  Fourier transforms of x- and y-components on the boundary
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
!  The vector potential needs to be known outside of (l1:l2,m1:m2) as well
!
        call communicate_bc_aa_pot(f,topbot)
!
!  Potential field condition at the top
!
      case('top')

        do j=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          delta_z  = z(n2+j) - z(n2-j)
          exp_fact = exp(-kappa*delta_z)
!
!  Determine potential field in ghost zones
!
          !  Fourier transforms of x- and y-components on the boundary
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
!  The vector potential needs to be known outside of (l1:l2,m1:m2) as well
!
        call communicate_bc_aa_pot(f,topbot)

      case default

        if (lroot) print*,"bc_aa_pot2: invalid argument"

      endselect

    endsubroutine bc_aa_pot3
!***********************************************************************
    subroutine bc_aa_pot2(f,topbot)
!
!  Potential field boundary condition
!
!  10-oct-06/tobi: Coded
!
      use Fourier, only: fourier_transform_xy_xy, fourier_transform_y_y
      use Mpicomm, only: communicate_bc_aa_pot

      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot

      real, dimension (nx,ny,iax:iaz) :: aa_re,aa_im
      real, dimension (nx,ny) :: kx,ky,kappa,kappa1
      real, dimension (nx,ny) :: tmp_re,tmp_im
      real, dimension (nx,ny) :: fac
      integer :: i,j
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
      where (kappa > 0)
        kappa1 = 1/kappa
      elsewhere
        kappa1 = 0
      endwhere
!
!  Check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  Potential field condition at the bottom
!
      case('bot')
!
!  Fourier transforms of x- and y-components on the boundary
!
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
!
!  Determine potential field in ghost zones
!
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
!  The vector potential needs to be known outside of (l1:l2,m1:m2) as well
!
        call communicate_bc_aa_pot(f,topbot)
!
!  Potential field condition at the top
!
      case('top')
!
!  Fourier transforms of x- and y-components on the boundary
!
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
!
!  Determine potential field in ghost zones
!
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
!  The vector potential needs to be known outside of (l1:l2,m1:m2) as well
!
        call communicate_bc_aa_pot(f,topbot)

      case default

        if (lroot) print*,"bc_aa_pot2: invalid argument"

      endselect

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
      use Mpicomm, only: stop_it,communicate_bc_aa_pot
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny) :: f2,f3
      real, dimension (nx,ny,nghost+1) :: fz
      integer :: j
!
!  potential field condition
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  potential field condition at the bottom
!
      case('bot')
        if (headtt) print*,'bc_aa_pot: pot-field bdry cond at bottom'
        if (mod(nxgrid,nygrid)/=0) &
             call stop_it("bc_aa_pot: pot-field doesn't work "//&
                          "with mod(nxgrid,nygrid)/=1")
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
        call communicate_bc_aa_pot(f,topbot)
!
!  potential field condition at the top
!
      case('top')
        if (headtt) print*,'bc_aa_pot: pot-field bdry cond at top'
        if (mod(nxgrid,nygrid)/=0) &
             call stop_it("bc_aa_pot: pot-field doesn't work "//&
                          "with mod(nxgrid,nygrid)/=1")
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
        call communicate_bc_aa_pot(f,topbot)
      case default
        if (lroot) print*,"bc_aa_pot: invalid argument"
      endselect
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
      use Fourier
!
      real, dimension (nx,ny) :: fac,kk,f1r,f1i,g1r,g1i,f2,f2r,f2i,f3,f3r,f3i
      real, dimension (nx,ny,nghost+1) :: fz
      real, dimension (nx) :: kx
      real, dimension (nygrid) :: ky
      real :: delz
      integer :: i,irev
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
!
      kx=cshift((/(i-(nx-1)/2,i=0,nx-1)/),+(nx-1)/2)*2*pi/Lx
      ky=cshift((/(i-(nygrid-1)/2,i=0,nygrid-1)/),+(nygrid-1)/2)*2*pi/Ly
!
!  calculate 1/k^2, zero mean
!
      kk=sqrt(spread(kx**2,2,ny)+spread(ky(ipy*ny+1:(ipy+1)*ny)**2,1,nx))
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
      use Fourier
!
      real, dimension (nx,ny) :: fac,kk,kkkx,kkky,f1r,f1i,g1r,g1i,f2,f2r,f2i,f3,f3r,f3i
      real, dimension (nx,ny,nghost+1) :: fz
      real, dimension (nx) :: kx
      real, dimension (nygrid) :: ky
      real :: delz
      integer :: i,irev
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
    endsubroutine potentdiv
!***********************************************************************
    subroutine bc_wind_z(f,topbot,massflux)
!
!  Calculates u_0 so that rho*(u+u_0)=massflux 
!  massflux can be set as fbcz1/2(rho) in run.in
!  
!  18-06-2008/bing: coded
!
      use Mpicomm, only: mpisend_real,mpirecv_real,stop_it

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
         if (.not.(lequidist(1) .and. lequidist(2))) &
              call stop_it("bc_wind_z:non equidistant grid in x and y not implemented")
         if (nprocx .gt. 1)  &
              call stop_it('bc_wind: nprocx > 1 not yet implemented')
!
!   check for warnings
!
         if (.not. ldensity)  &
              call warning('bc_wind',"no defined density, using rho=1 ?")
      endif
!
      select case(topbot)
!
!  Bottom boundary.
!
      case('bot')
         ntb = n1
!
!  Top boundary.
!
       case('top')
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
    subroutine bc_nscbc_prf_x(f,df,topbot,non_reflecting_inlet,linlet,u_t)
!
!   Calculate du and dlnrho at a partially reflecting outlet/inlet normal to 
!   x-direction acc. to LODI relations. Uses a one-sided finite diff. stencil.
!
!   7-jul-08/arne: coded.
!  25-nov-08/nils: extended to work in multiple dimensions and with cross terms
!                  i.e. not just the LODI equations.
!
      use MpiComm, only: stop_it
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice, der_pencil, der2_pencil
      use Chemistry
      use Viscosity

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      logical, optional :: linlet
      logical :: llinlet, non_reflecting_inlet
      real, optional :: u_t
      real, dimension(ny,nz) :: dlnrho_dx
      real, dimension(ny,nz) :: rho0, L_1, L_3, L_4, L_5,parallell_term_uz
      real, dimension(ny,nz) :: parallell_term_rho,dlnrho_dy,dlnrho_dz
      real, dimension(ny,nz) :: parallell_term_ux,d2u1_dy2,d2u1_dz2
      real, dimension(ny,nz) :: d2u2_dy2,d2u2_dz2,d2u3_dy2,d2u3_dz2
      real, dimension(ny,nz) :: prefac1, prefac2,parallell_term_uy
      real, dimension(ny,nz,3) :: div_rho
      real, dimension(ny,nz,3,3) :: dui_dxj
      real, dimension (my,mz) :: cs0_ar,cs20_ar
      real, dimension (mx,my,mz) :: cs2_full
      real, dimension (my,mz) :: tmp22,tmp12,tmp2_lnrho,tmp33,tmp13,tmp3_lnrho
      real, dimension (my,mz) :: tmp23,tmp32
      real, dimension (my) :: tmpy
      real, dimension (mz) :: tmpz
      real :: Mach,KK,nu
      integer lll,i
      integer sgn

      intent(inout) :: f
      intent(inout) :: df

      llinlet = .false.
      if (present(linlet)) llinlet = linlet
      if (llinlet.and..not.present(u_t)) call stop_it(&
           'bc_nscbc_prf_x: when using linlet=T, you must also specify u_t)')
      select case(topbot)
      case('bot')
        lll = l1
        sgn = 1
      case('top')
        lll = l2
        sgn = -1
      case default
        print*, "bc_nscbc_prf_x: ", topbot, " should be `top' or `bot'"
      endselect
!
!  Find density
!
      if (ldensity_nolog) then
        rho0 = f(lll,m1:m2,n1:n2,ilnrho)
      else
        rho0 = exp(f(lll,m1:m2,n1:n2,ilnrho))
      endif
!
!  Get viscoity
!
      call getnu(nu)
!
!  Set arrays for the speed of sound and for the speed of sound squared (is it
!  really necessarry to have both arrays?) 
!  Set prefactors to be used later.
!
      if (leos_chemistry) call get_cs2_full(cs2_full)
      if (leos_idealgas) then
        cs20_ar(m1:m2,n1:n2)=cs20
        cs0_ar(m1:m2,n1:n2)=cs0
        prefac1 = -1./(2.*cs20)
        prefac2 = -1./(2.*rho0*cs0)
      elseif (leos_chemistry) then
        cs20_ar=cs2_full(lll,:,:)
        cs0_ar=cs20_ar**0.5
        prefac1 = -1./(2.*cs20_ar(m1:m2,n1:n2))
        prefac2 = -1./(2.*rho0*cs0_ar(m1:m2,n1:n2))
      else
        print*,"bc_nscbc_prf_x: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC boundary treatment only implemented for an ideal gas." 
        print*,"Boundary treatment skipped."
        return
      endif
!
!  Calculate one-sided derivatives in the boundary normal direction
!
      call der_onesided_4_slice(f,sgn,ilnrho,div_rho(:,:,1),lll,1)
      call der_onesided_4_slice(f,sgn,iux,dui_dxj(:,:,1,1),lll,1)
      call der_onesided_4_slice(f,sgn,iuy,dui_dxj(:,:,2,1),lll,1)
      call der_onesided_4_slice(f,sgn,iuz,dui_dxj(:,:,3,1),lll,1)
!
!  Do central differencing in the directions parallell to the boundary 
!  first in the y-direction......
!
      if (nygrid /= 1) then
        do i=n1,n2
          call der_pencil(2,f(lll,:,i,iuy),tmp22(:,i))
          call der_pencil(2,f(lll,:,i,ilnrho),tmp2_lnrho(:,i))
          call der_pencil(2,f(lll,:,i,iux),tmp12(:,i))
          call der_pencil(2,f(lll,:,i,iuz),tmp32(:,i))
          call der2_pencil(2,f(lll,:,i,iux),tmpy)
          d2u1_dy2(:,i)=tmpy(m1:m2)
          call der2_pencil(2,f(lll,:,i,iuy),tmpy)
          d2u2_dy2(:,i)=tmpy(m1:m2)
          call der2_pencil(2,f(lll,:,i,iuz),tmpy)
          d2u3_dy2(:,i)=tmpy(m1:m2)
        enddo
      else
        tmp32=0
        tmp22=0
        tmp12=0
        tmp2_lnrho=0
        d2u1_dy2=0
        d2u2_dy2=0
        d2u3_dy2=0
      endif
      dui_dxj(:,:,1,2)=tmp12(m1:m2,n1:n2)
      dui_dxj(:,:,2,2)=tmp22(m1:m2,n1:n2)
      dui_dxj(:,:,3,2)=tmp32(m1:m2,n1:n2)
      div_rho(:,:,2)=tmp2_lnrho(m1:m2,n1:n2)
!
!  .... then in the z-direction
!
      if (nzgrid /= 1) then
        do i=m1,m2
          call der_pencil(3,f(lll,i,:,iuz),tmp33(i,:))
          call der_pencil(3,f(lll,i,:,ilnrho),tmp3_lnrho(i,:))
          call der_pencil(3,f(lll,i,:,iux),tmp13(i,:))
          call der_pencil(3,f(lll,i,:,iuy),tmp23(i,:))
          call der2_pencil(3,f(lll,i,:,iux),tmpz)
          d2u1_dz2(i,:)=tmpz(n1:n2)
          call der2_pencil(3,f(lll,i,:,iuy),tmpz)
          d2u2_dz2(i,:)=tmpz(n1:n2)
          call der2_pencil(3,f(lll,i,:,iuz),tmpz)
          d2u3_dz2(i,:)=tmpz(n1:n2)
        enddo
      else
        tmp33=0
        tmp23=0
        tmp13=0
        tmp3_lnrho=0
        d2u1_dz2=0
        d2u2_dz2=0
        d2u3_dz2=0
      endif
      dui_dxj(:,:,1,3)=tmp13(m1:m2,n1:n2)
      dui_dxj(:,:,2,3)=tmp23(m1:m2,n1:n2)
      dui_dxj(:,:,3,3)=tmp33(m1:m2,n1:n2)
      div_rho(:,:,3)=tmp3_lnrho(m1:m2,n1:n2)
!
!  Find divergence of rho if we solve for logarithm of rho
!
      if (.not. ldensity_nolog) then
        do i=1,3
          div_rho(:,:,i)=div_rho(:,:,i)*rho0
        enddo
      endif
!
!  Find the L_i's (which really is the Lodi equations)
!
      if (llinlet) then
!  This L_1 is correct only for p=cs^2*rho
        L_1 = (f(lll,m1:m2,n1:n2,iux) - sgn*cs0_ar(m1:m2,n1:n2))&
            *(cs20_ar(m1:m2,n1:n2)*div_rho(:,:,1) - sgn*rho0*cs0_ar(m1:m2,n1:n2)&
            *dui_dxj(:,:,1,1))
        L_3=0
        L_4=0
        if (non_reflecting_inlet) then
!
!  The inlet in non-reflecting only when nscbc_sigma_in is set to 0, this 
!  might however lead to problems as the inlet velocity will tend to drift 
!  away from the target velocity u_t. This problem should be overcome by 
!  setting a small but non-zero nscbc_sigma_in.
!
          L_5 = nscbc_sigma_in*cs20_ar(m1:m2,n1:n2)*rho0&
              *sgn*(f(lll,m1:m2,n1:n2,iux)-u_t)
        else
          L_5 = L_1
        endif
      else
!
!  Find Mach number 
!  (NILS: I do not think this is a good way to determine the Mach
!  number since this is a local Mach number for this processor. Furthermore
!  by determining the Mach number like this we will see the Mach number varying
!  with the phase of an acoustic wave as the wave pass through the boundary.
!  I think that what we really want is a Mach number averaged over the 
!  timescale of several acoustic waves. How could this be done????)
!
        Mach=sum(f(lll,m1:m2,n1:n2,iux)/cs0_ar(m1:m2,n1:n2))/(ny*nz)
!
!  Find the parameter determining 
!
        KK=nscbc_sigma_out*(1-Mach**2)*cs0/Lxyz(1)
!
!  Find the L_i's
!
        L_1 = KK*(rho0*cs20-p_infty)
        L_3 = f(lll,m1:m2,n1:n2,iux)*dui_dxj(:,:,2,1)
        L_4 = f(lll,m1:m2,n1:n2,iux)*dui_dxj(:,:,3,1)
        L_5 = (f(lll,m1:m2,n1:n2,iux) - sgn*cs0_ar(m1:m2,n1:n2))*&
             (cs20_ar(m1:m2,n1:n2)*div_rho(:,:,1)&
             - sgn*rho0*cs0_ar(m1:m2,n1:n2)*dui_dxj(:,:,1,1))
      end if
!
!  Add terms due to derivatives parallell to the boundary
!  NILS: Viscous terms in the x direction are missing!
!
      parallell_term_rho &
           =rho0*dui_dxj(:,:,2,2)+f(lll,m1:m2,n1:n2,iuy)*div_rho(:,:,2)&
           +rho0*dui_dxj(:,:,3,3)+f(lll,m1:m2,n1:n2,iuz)*div_rho(:,:,3)
      parallell_term_ux &
           =f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,1,2)&
           +f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,1,3)&
           +nu*(d2u1_dy2+d2u1_dz2)
      parallell_term_uy &
           =f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,2,2)&
           +f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,2,3)&
           +cs20_ar(m1:m2,n1:n2)*div_rho(:,:,2)/rho0&
           +nu*(d2u2_dy2+d2u2_dz2)
      parallell_term_uz &
           =f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,3,2)&
           +f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,3,3)&
           +cs20_ar(m1:m2,n1:n2)*div_rho(:,:,3)/rho0&
           +nu*(d2u3_dy2+d2u3_dz2)






!
!  NILS: Currently the implementation with the parallell terms does not 
!  NILS: seem to work. Set at parallell terms to zero for now.
!
      parallell_term_rho=0
      parallell_term_ux=0
      parallell_term_uy=0
      parallell_term_uz=0
!
!  Find the evolution equations at the boundary
!
      select case(topbot)
      ! NB: For 'top' L_1 plays the role of L5 and L_5 the role of L1
      case('bot')
        df(lll,m1:m2,n1:n2,ilnrho) = prefac1*(L_5 + L_1)-parallell_term_rho
        df(lll,m1:m2,n1:n2,iux) = prefac2*(L_1 - L_5)-parallell_term_ux
        df(lll,m1:m2,n1:n2,iuy) = -L_3-parallell_term_uy
        df(lll,m1:m2,n1:n2,iuz) = -L_4-parallell_term_uz
      case('top')
        df(lll,m1:m2,n1:n2,ilnrho) = prefac1*(L_1 + L_5)-parallell_term_rho
        df(lll,m1:m2,n1:n2,iux) = prefac2*(L_5 - L_1)-parallell_term_ux
        df(lll,m1:m2,n1:n2,iuy) = -L_3-parallell_term_uy
        df(lll,m1:m2,n1:n2,iuz) = -L_4-parallell_term_uz
      endselect
!
!  Check if we are solving for logrho or rho
!
      if (.not. ldensity_nolog) then
        df(lll,m1:m2,n1:n2,ilnrho)=df(lll,m1:m2,n1:n2,ilnrho)/rho0
      endif
!
! Impose required variables at the boundary
!
      if (llinlet) then
        if (.not. non_reflecting_inlet) f(lll,m1:m2,n1:n2,iux) = u_t
        f(lll,m1:m2,n1:n2,iuy) = 0
        f(lll,m1:m2,n1:n2,iuz) = 0
      endif
!
    endsubroutine bc_nscbc_prf_x
!***********************************************************************
    subroutine bc_nscbc_prf_y(f,df,topbot,non_reflecting_inlet,linlet,u_t)
!
!   Calculate du and dlnrho at a partially reflecting outlet/inlet normal to 
!   y-direction acc. to LODI relations. Uses a one-sided finite diff. stencil.
!
!   7-jul-08/arne: coded.
!  25-nov-08/nils: extended to work in multiple dimensions and with cross terms
!                  i.e. not just the LODI equations.
!
      use MpiComm, only: stop_it
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice, der_pencil, der2_pencil
      use Chemistry
      use Viscosity

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      logical, optional :: linlet
      logical :: llinlet, non_reflecting_inlet
      real, optional :: u_t
      !real, parameter :: sigma = 1.
      real, dimension(nx,nz) :: dlnrho_dx,dlnrho_dy,dlnrho_dz
      real, dimension(nx,nz) :: rho0, L_1, L_3, L_4, L_5,parallell_term_uz
      real, dimension(nx,nz) :: parallell_term_rho
      real, dimension(nx,nz) :: parallell_term_ux,d2u1_dx2,d2u1_dz2
      real, dimension(nx,nz) :: d2u2_dx2,d2u2_dz2,d2u3_dx2,d2u3_dz2
      real, dimension(nx,nz) :: prefac1, prefac2,parallell_term_uy
      real, dimension(nx,nz,3) :: div_rho
      real, dimension(nx,nz,3,3) :: dui_dxj
      real, dimension (mx,mz) :: cs0_ar,cs20_ar
      real, dimension (mx,my,mz) :: cs2_full
      real, dimension (mx,mz) :: tmp22,tmp12,tmp2_lnrho,tmp33,tmp13,tmp3_lnrho
      real, dimension (mx,mz) :: tmp23,tmp32
      real, dimension (mx) :: tmpx
      real, dimension (mz) :: tmpz
      real :: Mach,KK,nu
      integer lll,i
      integer sgn

      intent(inout) :: f
      intent(inout) :: df

      llinlet = .false.
      if (present(linlet)) llinlet = linlet
      if (llinlet.and..not.present(u_t)) call stop_it(&
           'bc_nscbc_prf_y: when using linlet=T, you must also specify u_t)')
      select case(topbot)
      case('bot')
        lll = m1
        sgn = 1
      case('top')
        lll = m2
        sgn = -1
      case default
        print*, "bc_nscbc_prf_y: ", topbot, " should be `top' or `bot'"
      endselect
!
!  Find density
!
      if (ldensity_nolog) then
        rho0 = f(l1:l2,lll,n1:n2,ilnrho)
      else
        rho0 = exp(f(l1:l2,lll,n1:n2,ilnrho))
      endif
!
!  Get viscoity
!
      call getnu(nu)
!
!  Set arrays for the speed of sound and for the speed of sound squared (is it
!  really necessarry to have both arrays?) 
!  Set prefactors to be used later.
!
      if (leos_chemistry) call get_cs2_full(cs2_full)
      if (leos_idealgas) then
        cs20_ar(l1:l2,n1:n2)=cs20
        cs0_ar(l1:l2,n1:n2)=cs0
        prefac1 = -1./(2.*cs20)
        prefac2 = -1./(2.*rho0*cs0)
      elseif (leos_chemistry) then
        cs20_ar=cs2_full(:,lll,:)
        cs0_ar=cs20_ar**0.5
        prefac1 = -1./(2.*cs20_ar(l1:l2,n1:n2))
        prefac2 = -1./(2.*rho0*cs0_ar(l1:l2,n1:n2))
      else
        print*,"bc_nscbc_prf_y: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC boundary treatment only implemented for an ideal gas." 
        print*,"Boundary treatment skipped."
        return
      endif
!
!  Calculate one-sided derivatives in the boundary normal direction
!
      call der_onesided_4_slice(f,sgn,ilnrho,div_rho(:,:,2),lll,2)
      call der_onesided_4_slice(f,sgn,iux,dui_dxj(:,:,1,2),lll,2)
      call der_onesided_4_slice(f,sgn,iuy,dui_dxj(:,:,2,2),lll,2)
      call der_onesided_4_slice(f,sgn,iuz,dui_dxj(:,:,3,2),lll,2)
!
!  Do central differencing in the directions parallell to the boundary 
!  first in the x-direction......
!
      if (nxgrid /= 1) then
        do i=n1,n2
          call der_pencil(1,f(:,lll,i,iuy),tmp22(:,i))
          call der_pencil(1,f(:,lll,i,ilnrho),tmp2_lnrho(:,i))
          call der_pencil(1,f(:,lll,i,iux),tmp12(:,i))
          call der_pencil(1,f(:,lll,i,iuz),tmp32(:,i))
          call der2_pencil(1,f(:,lll,i,iux),tmpx)
          d2u1_dx2(:,i)=tmpx(l1:l2)
          call der2_pencil(1,f(:,lll,i,iuy),tmpx)
          d2u2_dx2(:,i)=tmpx(l1:l2)
          call der2_pencil(1,f(:,lll,i,iuz),tmpx)
          d2u3_dx2(:,i)=tmpx(l1:l2)
        enddo
      else
        tmp32=0
        tmp22=0
        tmp12=0
        tmp2_lnrho=0
        d2u1_dx2=0
        d2u2_dx2=0
        d2u3_dx2=0
      endif
      dui_dxj(:,:,3,2)=tmp32(l1:l2,n1:n2)
      dui_dxj(:,:,2,2)=tmp22(l1:l2,n1:n2)
      dui_dxj(:,:,1,2)=tmp12(l1:l2,n1:n2)
      div_rho(:,:,1)=tmp2_lnrho(l1:l2,n1:n2)
!
!  .... then in the z-direction
!
      if (nzgrid /= 1) then
        do i=l1,l2
          call der_pencil(3,f(i,lll,:,iuz),tmp33(i,:))
          call der_pencil(3,f(i,lll,:,ilnrho),tmp3_lnrho(i,:))
          call der_pencil(3,f(i,lll,:,iux),tmp13(i,:))
          call der_pencil(3,f(i,lll,:,iuy),tmp23(i,:))
          call der2_pencil(3,f(i,lll,:,iux),tmpz)
          d2u1_dz2(i,:)=tmpz(n1:n2)
          call der2_pencil(3,f(i,lll,:,iuy),tmpz)
          d2u2_dz2(i,:)=tmpz(n1:n2)
          call der2_pencil(3,f(i,lll,:,iuz),tmpz)
          d2u3_dz2(i,:)=tmpz(n1:n2)
        enddo
      else
        tmp33=0
        tmp23=0
        tmp13=0
        tmp3_lnrho=0
        d2u1_dz2=0
        d2u2_dz2=0
        d2u3_dz2=0
      endif
      dui_dxj(:,:,3,3)=tmp33(l1:l2,n1:n2)
      dui_dxj(:,:,2,3)=tmp23(l1:l2,n1:n2)
      dui_dxj(:,:,1,3)=tmp13(l1:l2,n1:n2)
      div_rho(:,:,3)=tmp3_lnrho(l1:l2,n1:n2)
!
!  Find divergence of rho if we solve for logarithm of rho
!
      if (.not. ldensity_nolog) then
        do i=1,3
          div_rho(:,:,i)=div_rho(:,:,i)*rho0
        enddo
      endif
!
!  Find the L_i's (which really is the Lodi equations)
!
      if (llinlet) then
        L_1 = (f(l1:l2,lll,n1:n2,iuy) - sgn*cs0_ar(l1:l2,n1:n2))*&
            (cs20*div_rho(:,:,2) - sgn*rho0*cs0_ar(l1:l2,n1:n2)*dui_dxj(:,:,2,2))
        L_3=0
        L_4=0
        if (non_reflecting_inlet) then
!
!  The inlet in non-reflecting only when nscbc_sigma_in is set to 0, this 
!  might however lead to problems ..... NILS: Must SPECIFY!!!!!!!!!!
!
          L_5 = nscbc_sigma_in*cs20_ar(l1:l2,n1:n2)*rho0&
              *(sgn*f(l1:l2,lll,n1:n2,iuy)-sgn*u_t)
        else
          L_5 = L_1
        endif
      else
!
!  Find Mach number 
!  (NILS: I do not think this is a good way to determine the Mach
!  number since this is a local Mach number for this processor. Furthermore
!  by determining the Mach number like this we will see the Mach number varying
!  with the phase of an acoustic wave as the wave pass through the boundary.
!  I think that what we really want is a Mach number averaged over the 
!  timescale of several acoustic waves. How could this be done????)
!
        Mach=sum(f(l1:l2,lll,n1:n2,iuy)/cs0_ar(l1:l2,n1:n2))/(nx*nz)
!
!  Find the parameter determining 
!
        KK=nscbc_sigma_out*(1-Mach**2)*cs0/Lxyz(1)
!
!  Find the L_i's
!
        L_1 = KK*(rho0*cs20-p_infty)
        L_3 = f(l1:l2,lll,n1:n2,iuy)*dui_dxj(:,:,1,2)
        L_4 = f(l1:l2,lll,n1:n2,iuy)*dui_dxj(:,:,3,2)
        L_5 = (f(l1:l2,lll,n1:n2,iuy) - sgn*cs0_ar(l1:l2,n1:n2))*&
             (cs20_ar(l1:l2,n1:n2)*div_rho(:,:,2)&
             -sgn*rho0*cs0_ar(l1:l2,n1:n2)*dui_dxj(:,:,2,2))
      end if
!
!  Add terms due to derivatives parallell to the boundary
!

      parallell_term_rho &
          =rho0*dui_dxj(:,:,1,1)+f(l1:l2,lll,n1:n2,iux)*div_rho(:,:,1)&
          +rho0*dui_dxj(:,:,3,3)+f(l1:l2,lll,n1:n2,iuz)*div_rho(:,:,3)
      parallell_term_ux &
           =f(l1:l2,lll,n1:n2,iux)*dui_dxj(:,:,1,1)&
           +f(l1:l2,lll,n1:n2,iuz)*dui_dxj(:,:,1,3)&
           +cs20_ar(l1:l2,n1:n2)*div_rho(:,:,1)/rho0&
           +nu*(d2u1_dx2+d2u1_dz2)
      parallell_term_uy &
           =f(l1:l2,lll,n1:n2,iux)*dui_dxj(:,:,2,1)&
           +f(l1:l2,lll,n1:n2,iuz)*dui_dxj(:,:,2,3)&
           +nu*(d2u2_dx2+d2u2_dz2)
      parallell_term_uz &
           =f(l1:l2,lll,n1:n2,iux)*dui_dxj(:,:,3,1)&
           +f(l1:l2,lll,n1:n2,iuz)*dui_dxj(:,:,3,3)&
           +cs20_ar(l1:l2,n1:n2)*div_rho(:,:,3)/rho0&
           +nu*(d2u3_dx2+d2u3_dz2)



      parallell_term_rho=0
      parallell_term_ux=0
      parallell_term_uy=0
      parallell_term_uz=0
!
!  Find the evolution equations at the boundary
!
      select case(topbot)
      ! NB: For 'top' L_1 plays the role of L5 and L_5 the role of L1
      case('bot')
        df(l1:l2,lll,n1:n2,ilnrho) = prefac1*(L_5 + L_1)-parallell_term_rho
        df(l1:l2,lll,n1:n2,iuy) = prefac2*(L_1 - L_5)-parallell_term_uy
        df(l1:l2,lll,n1:n2,iux) = -L_3-parallell_term_ux
        df(l1:l2,lll,n1:n2,iuz) = -L_4-parallell_term_uz
      case('top')
        df(l1:l2,lll,n1:n2,ilnrho) = prefac1*(L_1 + L_5)-parallell_term_rho
        df(l1:l2,lll,n1:n2,iuy) = prefac2*(L_5 - L_1)-parallell_term_uy
        df(l1:l2,lll,n1:n2,iux) = -L_3-parallell_term_ux
        df(l1:l2,lll,n1:n2,iuz) = -L_4-parallell_term_uz
      endselect
!
!  Check if we are solving for logrho or rho
!
      if (.not. ldensity_nolog) then
        df(l1:l2,lll,n1:n2,ilnrho)=df(l1:l2,lll,n1:n2,ilnrho)/rho0
      endif
!
! Impose required variables at the boundary
!
      if (llinlet) then
        f(l1:l2,lll,n1:n2,iux) = 0
        if (.not. non_reflecting_inlet) f(l1:l2,lll,n1:n2,iuy) = u_t
        f(l1:l2,lll,n1:n2,iuz) = 0
      endif
!
    endsubroutine bc_nscbc_prf_y
!***********************************************************************
    subroutine bc_nscbc_prf_z(f,df,topbot,non_reflecting_inlet,linlet,u_t)
!
!   Calculate du and dlnrho at a partially reflecting outlet/inlet normal to 
!   z-direction acc. to LODI relations. Uses a one-sided finite diff. stencil.
!
!   7-jul-08/arne: coded.
!
      use MpiComm, only: stop_it
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      logical, optional :: linlet
      logical :: llinlet,non_reflecting_inlet
      real, optional :: u_t
      !real, parameter :: sigma = 1.
      real, dimension(nx,ny) :: du_dz, dlnrho_dz, rho0, L_1, L_5
      real, dimension(nx,ny) :: dp_prefac, prefac1, prefac2
      integer lll,sgn

      intent(in) :: f
      intent(out) :: df

      call fatal_error('bc_nscbc_prf_z',&
           'This sub is not properly implemented yet! Adapt from bc_nscbc_prf_x')

      llinlet = .false.
      if (present(linlet)) llinlet = linlet
      if (llinlet.and..not.present(u_t)) call stop_it(&
           'bc_nscbc_prf_z: when using linlet=T, you must also specify u_t)')

      select case(topbot)
      case('bot')
        lll = n1
        sgn = 1
      case('top')
        lll = n2
        sgn = -1
      case default
        print*, "bc_nscbc_prf_z: ", topbot, " should be `top' or `bot'"
      endselect
      if (leos_idealgas) then
        if (ldensity_nolog) then
          rho0 = f(l1:l2,m1:m2,lll,ilnrho)
          ! ``dp = cs20*drho''
          dp_prefac = cs20
        else
          rho0 = exp(f(l1:l2,m1:m2,lll,ilnrho))
          ! ``dp = cs20*rho0*dlnrho''
          dp_prefac = cs20*rho0
        endif
        prefac1 = -1./(2*rho0*cs20)
        prefac2 = -1./(2*rho0*cs0)
      else
        print*,"bc_nscbc_prf_y: leos_idealgas=",leos_idealgas,"." 
        print*,"NSCBC boundary treatment only implemented for an ideal gas."
        print*,"Boundary treatment skipped."
        return
      endif
      call der_onesided_4_slice(f,sgn,ilnrho,dlnrho_dz,lll,3)
      call der_onesided_4_slice(f,sgn,iuz,du_dz,lll,3)
      L_1 = (f(l1:l2,m1:m2,lll,iuz) - sgn*cs0)*&
            (dp_prefac*dlnrho_dz - sgn*rho0*cs0*du_dz)
      if (llinlet) then
        L_5 = nscbc_sigma_in*cs20*rho0*(sgn*f(l1:l2,m1:m2,lll,iuz)-sgn*u_t)
      else
        L_5 = 0
      end if
      select case(topbot)
      ! NB: For 'top' L_1 plays the role of L5 and L_5 the role of L1
      case('bot')
        df(l1:l2,m1:m2,lll,ilnrho) = prefac1*(L_5 + L_1)
        df(l1:l2,m1:m2,lll,iuz) = prefac2*(L_5 - L_1)
      case('top')
        df(l1:l2,m1:m2,lll,ilnrho) = prefac1*(L_1 + L_5)
        df(l1:l2,m1:m2,lll,iuz) = prefac2*(L_1 - L_5)
      endselect
!
    endsubroutine bc_nscbc_prf_z
!***********************************************************************
    subroutine bc_nscbc_subin_x(f,df,topbot,val)
!
!   nscbc case 
!   subsonic inflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    ux,uy,T are fixed, drho  is calculated 
!
!   16-nov-08/natalia: coded.
!
      use MpiComm, only: stop_it
      !use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice, der_pencil
      use Chemistry
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: mom2
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension(ny,nz) :: dux_dx, L_1, L_2, L_5, dpp_dx
      real, dimension(my,mz) :: rho0, gamma0, dmom2_dy, TT0 
      real, dimension (my,mz) :: cs2x,cs0_ar,cs20_ar
      real, dimension (mx,my,mz) :: cs2_full, gamma_full, rho_full, pp
      real, dimension(nchemspec) :: YYi
      real, dimension (mcom), optional :: val
      integer :: lll,i, sgn,k
      real :: Mach_num, u_t, T_t
!
      intent(inout) :: f
      intent(out) :: df
!

     if (.not.present(val)) call stop_it(&
           'bc_nscbc_subin_x: you must specify fbcx)')

      u_t=val(iux)
      T_t=val(ilnTT)
      do k=1,nchemspec
       YYi(k)=val(ichemspec(k))
      enddo

      if (leos_chemistry) then
        call get_cs2_full(cs2_full)
        call get_gamma_full(gamma_full)
      endif
!
      select case(topbot)
      case('bot')
        lll = l1
        sgn = 1
      case('top')
        lll = l2
        sgn = -1
      case default
        print*, "bc_nscbc_subin_x: ", topbot, " should be `top' or `bot'"
      endselect
!
      if (leos_chemistry) then
         cs20_ar=cs2_full(lll,:,:)
         cs0_ar=cs2_full(lll,:,:)**0.5
         gamma0=gamma_full(lll,:,:)
         TT0=exp(f(lll,:,:,ilnTT))
         rho_full=exp(f(:,:,:,ilnrho))
         rho0(:,:) = rho_full(lll,:,:)
         mom2(lll,:,:)=rho0(:,:)*f(lll,:,:,iuy)
         do i=1,my
         do k=1,mz
          if (minval(gamma_full(:,i,k))<=0.) then
           pp(:,i,k)=0.
          else
           pp(:,i,k)=cs2_full(:,i,k)*rho_full(:,i,k)/gamma_full(:,i,k)
          endif
         enddo
         enddo
      else
        print*,"bc_nscbc_subin_x: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC subsonic inflos is only implemented for "//&
            "the chemistry case." 
        print*,"Boundary treatment skipped."
        return
      endif
!
!
      !  call der_onesided_4_slice(f,sgn,ilnrho,dlnrho_dx,lll,1)
        call der_onesided_4_slice(f,sgn,iux,dux_dx,lll,1)
        call der_onesided_4_slice(pp,sgn,dpp_dx,lll,1)
!
        do i=1,mz
          call der_pencil(2,mom2(lll,:,i),dmom2_dy(:,i))
        enddo

        select case(topbot)
        case('bot')
          L_1 = (f(lll,m1:m2,n1:n2,iux) - cs0_ar(m1:m2,n1:n2))*&
              (dpp_dx - rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dux_dx)
          L_5 =L_1-2.*rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*&
              df(lll,m1:m2,n1:n2,iux)
        case('top')
          L_5 = (f(lll,m1:m2,n1:n2,iux) + cs0_ar(m1:m2,n1:n2))*&
              (dpp_dx + rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dux_dx)
          L_1 = L_5+2.*rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*&
              df(lll,m1:m2,n1:n2,iux)
        endselect
        L_2 = 0.5*(gamma0(m1:m2,n1:n2)-1.)*(L_5+L_1) &
            +rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2)*df(lll,m1:m2,n1:n2,ilnTT)
        if (ldensity_nolog) then
          df(lll,m1:m2,n1:n2,ilnrho) = -1./cs20_ar(m1:m2,n1:n2)*&
              (L_2+0.5*(L_5 + L_1)) ! -dmom2_dy(m1:m2,n1:n2)
        else
          df(lll,m1:m2,n1:n2,ilnrho) = &
              -1./rho0(m1:m2,n1:n2)/cs20_ar(m1:m2,n1:n2) &
              *(L_2+0.5*(L_5 + L_1)) !&
          ! -1./rho0(m1:m2,n1:n2)*dmom2_dy(m1:m2,n1:n2)
        endif

!
! Natalia: this subroutine is still under construction
! Please, do not remove this commented part !!!!!
!
  
!
!  this conditions can be important! 
!  check withour them
!

!        do k=1,nchemspec
!         f(lll,m1:m2,n1:n2,ichemspec(k))=YYi(k)  
!        enddo
   !      f(lll,m1:m2,n1:n2,iux) = u_t
   !      f(lll,m1:m2,n1:n2,ilnTT) = T_t
   !      df(lll,m1:m2,n1:n2,ilnTT)=0.


    endsubroutine bc_nscbc_subin_x
!***********************************************************************
    subroutine bc_nscbc_nref_subout_x(f,df,topbot)
!
!   nscbc case 
!   subsonic non-reflecting outflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    drho. dT, dux, duy  are calculated, p_inf can be 
!   fixed (if nscbc_sigma <>0)
!
!   16-nov-08/natalia: coded.
!
      use MpiComm, only: stop_it
      !use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice,der_pencil
      use Chemistry
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension (my,mz) :: rho0,gamma0
      real, dimension (mx,my,mz) :: mom2, rho_ux2, rho_uy2
      real, dimension (mx,my,mz) :: rho_gamma, rhoE_p, pp
      real, dimension (ny,nz) :: dux_dx, drho_dx, dpp_dx,dYk_dx
      real, dimension (ny,nz) :: drho_prefac,p_infx, KK, L_1, L_2, L_5
      real, dimension (my,mz) :: cs2x,cs0_ar,cs20_ar,dmom2_dy
      real, dimension (my,mz) :: drhoE_p_dy,dYk_dy, dux_dy
      real, dimension (ny,nz,nchemspec) :: bound_rhs_Y
      real, dimension (ny,nz) :: bound_rhs_T
      real, dimension (mx,my,mz) :: cs2_full, gamma_full, rho_full
!
      integer :: lll, sgn,i,j,k
      real :: Mach_num
!
      intent(inout) :: f
      intent(out) :: df
!
      if (leos_chemistry) then
        call get_cs2_full(cs2_full)
        call get_gamma_full(gamma_full)
      endif
!
      select case(topbot)
      case('bot')
        lll = l1
        sgn = 1
        if (leos_chemistry) then 
          call get_p_infx(p_infx,'bot')
        endif
      case('top')
        lll = l2
        sgn = -1
        if (leos_chemistry) then
          call get_p_infx(p_infx,'top')
        endif
      case default
        print*, "bc_nscbc_subin_x: ", topbot, " should be `top' or `bot'"
      endselect
!
      if (leos_chemistry) then
         cs20_ar=cs2_full(lll,:,:)
         cs0_ar=cs2_full(lll,:,:)**0.5
         gamma0=gamma_full(lll,:,:)
!
        if (ldensity_nolog) then
          rho_full = f(:,:,:,ilnrho)
          rho0(:,:) = f(lll,:,:,ilnrho)
          drho_prefac=-1./cs20_ar(m1:m2,n1:n2)
        else
          rho_full = exp(f(:,:,:,ilnrho))
          rho0(:,:) = rho_full(lll,:,:)
          drho_prefac=-1./rho0(m1:m2,n1:n2)/cs20_ar(m1:m2,n1:n2)
        endif
         
         do i=1,my
         do k=1,mz
          if (minval(gamma_full(:,i,k))<=0.) then
           pp(:,i,k)=0.
          else
           pp(:,i,k)=cs2_full(:,i,k)*rho_full(:,i,k)/gamma_full(:,i,k)
          endif
         enddo
         enddo

         mom2(lll,:,:)=rho0(:,:)*f(lll,:,:,iuy)
         rho_ux2(lll,:,:)=rho0(:,:)*f(lll,:,:,iux)*f(lll,:,:,iux)
         rho_uy2(lll,:,:)=rho0(:,:)*f(lll,:,:,iuy)*f(lll,:,:,iuy)
         do i=1,my
         do k=1,mz
          if (gamma0(i,k)<=0.) then
           rho_gamma(lll,i,k)=0.
           rhoE_p(lll,i,k)=0.
          else
           rho_gamma(lll,i,k)=rho0(i,k)/gamma0(i,k)
           rhoE_p(lll,i,k)=0.5*rho_ux2(lll,i,k)+0.5*rho_uy2(lll,i,k) &
                         +cs20_ar(i,k)/(gamma0(i,k)-1.)*rho0(i,k)
          endif
         enddo
         enddo
      else
        print*,"bc_nscbc_subin_x: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC subsonic inflos is only implemented "//&
            "for the chemistry case." 
        print*,"Boundary treatment skipped."
        return
      endif
!
      call der_onesided_4_slice(rho_full,sgn,drho_dx,lll,1)
      call der_onesided_4_slice(pp,sgn,dpp_dx,lll,1)
      call der_onesided_4_slice(f,sgn,iux,dux_dx,lll,1)
!
      do i=1,mz
        call der_pencil(2,mom2(lll,:,i),dmom2_dy(:,i))
        call der_pencil(2,f(lll,:,i,iux),dux_dy(:,i))
        call der_pencil(2,rhoE_p(lll,:,i),drhoE_p_dy(:,i))
      enddo
!
      Mach_num=maxval(f(lll,m1:m2,n1:n2,iux)/cs0_ar(m1:m2,n1:n2))
      KK=nscbc_sigma_out*(1.-Mach_num*Mach_num)*cs0_ar(m1:m2,n1:n2)/Lxyz(1)
!
      select case(topbot)
      case('bot')
        L_5=KK*(cs20_ar(m1:m2,n1:n2)/gamma0(m1:m2,n1:n2)*&
            rho0(m1:m2,n1:n2)-p_infx)
        L_1 = (f(lll,m1:m2,n1:n2,iux) - cs0_ar(m1:m2,n1:n2))*&
            (dpp_dx- rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dux_dx)
        call get_rhs_Y('bot',1,bound_rhs_Y)
        call get_rhs_T('bot',1,bound_rhs_T)
     case('top')
        L_1=KK*(cs20_ar(m1:m2,n1:n2)/gamma0(m1:m2,n1:n2)*&
            rho0(m1:m2,n1:n2)-p_infx)
        L_5 = (f(lll,m1:m2,n1:n2,iux) + cs0_ar(m1:m2,n1:n2))*&
            ( dpp_dx+ rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dux_dx)
        call get_rhs_Y('top',1,bound_rhs_Y)
        call get_rhs_T('top',1,bound_rhs_T)
      endselect
!
      L_2 = f(lll,m1:m2,n1:n2,iux)*(cs20_ar(m1:m2,n1:n2)*drho_dx-dpp_dx)
!
      if (ldensity_nolog) then
        df(lll,m1:m2,n1:n2,ilnrho) = &
            drho_prefac*(L_2+0.5*(L_5 + L_1))!-dmom2_dy(m1:m2,n1:n2)
      else
        df(lll,m1:m2,n1:n2,ilnrho) = &
            drho_prefac*(L_2+0.5*(L_5 + L_1)) !&
      !   -1./rho0(m1:m2,n1:n2)*dmom2_dy(m1:m2,n1:n2)
      endif
      df(lll,m1:m2,n1:n2,iux) = -1./&
          (2.*rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2))*(L_5 - L_1) !&
      !      -f(lll,m1:m2,n1:n2,iux)*dux_dy(m1:m2,n1:n2)
      df(lll,m1:m2,n1:n2,ilnTT) = -1./&
          (rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2))*(-L_2 &
          +0.5*(gamma0(m1:m2,n1:n2)-1.)*(L_5+L_1)) 
!        -1./(rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2))* &
!        (gamma0(m1:m2,n1:n2)-1.) &
!        *gamma0(m1:m2,n1:n2)*drhoE_p_dy(m1:m2,n1:n2) !&
    !    +bound_rhs_T(:,:)
!
     
!
      if (nchemspec>1) then
       do k=1,nchemspec
          call der_onesided_4_slice(f,sgn,ichemspec(k),dYk_dx,lll,1)
          do i=1,mz
            call der_pencil(2,f(lll,:,i,ichemspec(k)),dYk_dy(:,i))
          enddo
          df(lll,m1:m2,n1:n2,ichemspec(k))=&
              -f(lll,m1:m2,n1:n2,iux)*dYk_dx &
              -f(lll,m1:m2,n1:n2,iuy)*dYk_dy(m1:m2,n1:n2) &
              +bound_rhs_Y(:,:,k)
      
       ! if (lfilter) then
         do i=m1,m2
         do j=n1,n2
           if ((f(lll,i,j,ichemspec(k))+df(lll,i,j,ichemspec(k))*dt)<-1e-25 ) df(lll,i,j,ichemspec(k))=-1e-25*dt
           if ((f(lll,i,j,ichemspec(k))+df(lll,i,j,ichemspec(k))*dt)>1.) df(lll,i,j,ichemspec(k))=1.*dt
         enddo
         enddo
       ! endif

       enddo

      endif
!

      select case(topbot)
      case('bot')
       do i=1,nghost; f(l1-i,:,:,ilnTT)=2*f(l1,:,:,ilnTT)-f(l1+i,:,:,ilnTT); enddo 
      case('top')
        do i=1,nghost; f(l2+i,:,:,ilnTT)=2*f(l2,:,:,ilnTT)-f(l2-i,:,:,ilnTT); enddo
      case default
        print*, "bc_nscbc_subin_x: ", topbot, " should be `top' or `bot'"
      endselect



    endsubroutine bc_nscbc_nref_subout_x
!***********************************************************************
    subroutine bc_nscbc_nref_subout_y(f,df,topbot)
!
!   nscbc case 
!   subsonic non-reflecting outflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    drho. dT, dux, duy  are calculated, p_inf can be fixed 
!   (if nscbc_sigma <>0)
!
!   16-nov-08/natalia: coded.
!
      use MpiComm, only: stop_it
      !use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice, der_pencil
      use Chemistry
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension (mx,mz) :: rho0,gamma0
      real, dimension (mx,my,mz) :: mom2, rho_ux2, rho_uy2, rho_gamma, rhoE_p
      real, dimension (nx,nz) :: duy_dy, dlnrho_dy, L_1, L_2, L_5
      real, dimension (nx,nz) :: dp_prefac,drho_prefac,p_infy, KK, dpp_dy
      real, dimension (mx,mz) :: cs2y,cs0_ar,cs20_ar,dmom2_dx,drhoE_p_dx,duy_dx
      real, dimension (mx,my,mz) :: cs2_full,gamma_full,pp, rho_full
      integer :: mmm, sgn,i
      real :: Mach_num
!
      intent(in) :: f
      intent(out) :: df
!
      if (leos_chemistry) then 
        call get_cs2_full(cs2_full)
        call get_gamma_full(gamma_full)
      endif
!
      select case(topbot)
      case('bot')
        mmm = m1
        sgn = 1
        if (leos_chemistry) then 
!          call calc_cs2y(cs2y,'bot',f)
!          call get_gammay(gamma0,'bot')
          call get_p_infy(p_infy,'bot')
        endif
      case('top')
        mmm = m2
        sgn = -1
        if (leos_chemistry) then
!          call calc_cs2y(cs2y,'top',f)
!          call get_gammay(gamma0,'top')
          call get_p_infy(p_infy,'top')
        endif
      case default
        print*, "bc_nscbc_subin_y: ", topbot, " should be `top' or `bot'"
      endselect
!
      if (leos_chemistry) then
        cs20_ar=cs2_full(:,mmm,:)
        cs0_ar=cs2_full(:,mmm,:)**0.5
        gamma0=gamma_full(:,mmm,:)
        if (ldensity_nolog) then
          rho0(:,:) = f(:,mmm,:,ilnrho)
!          dp_prefac = cs20_ar(l1:l2,n1:n2)/gamma0(l1:l2,n1:n2)
          drho_prefac=-1./cs20_ar(l1:l2,n1:n2)
        else
          rho0(:,:) = exp(f(:,mmm,:,ilnrho))
!          dp_prefac = cs20_ar(l1:l2,n1:n2)*&
!          rho0(l1:l2,n1:n2)/gamma0(l1:l2,n1:n2)
          drho_prefac=-1./rho0(l1:l2,n1:n2)/cs20_ar(l1:l2,n1:n2)
        endif
!
        rho_full=exp(f(:,:,:,ilnrho))
        pp=cs2_full*rho_full/gamma_full
        mom2(:,mmm,:)=rho0(:,:)*f(:,mmm,:,iuy)
        rho_ux2(:,mmm,:)=rho0(:,:)*f(:,mmm,:,iux)*f(:,mmm,:,iux)
        rho_uy2(:,mmm,:)=rho0(:,:)*f(:,mmm,:,iuy)*f(:,mmm,:,iuy)
        rho_gamma(:,mmm,:)=rho0(:,:)/gamma0(:,:)
        rhoE_p(:,mmm,:)=0.5*rho_ux2(:,mmm,:)+0.5*rho_uy2(:,mmm,:) &
            +cs20_ar(:,:)/(gamma0(:,:)-1.)*rho0(:,:)
      else
        print*,"bc_nscbc_subin_y: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC subsonic inflow is only implemented "//&
            "for the chemistry case." 
        print*,"Boundary treatment skipped."
        return
      endif
!
      call der_onesided_4_slice(f,sgn,ilnrho,dlnrho_dy,mmm,2)
      call der_onesided_4_slice(f,sgn,iux,duy_dy,mmm,2)
      call der_onesided_4_slice(pp,sgn,dpp_dy,mmm,2)
!
      do i=1,mz
        call der_pencil(1,mom2(:,mmm,i),dmom2_dx(:,i))
        call der_pencil(1,f(:,mmm,i,iuy),duy_dx(:,i))
        call der_pencil(1,rhoE_p(:,mmm,i),drhoE_p_dx(:,i))
      enddo
!
      Mach_num=maxval(f(l1:l2,mmm,n1:n2,iuy)/cs0_ar(l1:l2,n1:n2))
      KK=nscbc_sigma_out*(1.-Mach_num*Mach_num)*cs0_ar(l1:l2,n1:n2)/Lxyz(2)
!
      select case(topbot)
      case('bot')
        L_5 = KK*(cs20_ar(l1:l2,n1:n2)/&
            gamma0(l1:l2,n1:n2)*rho0(l1:l2,n1:n2)-p_infy)
        L_1 = (f(l1:l2,mmm,n1:n2,iuy) - cs0_ar(l1:l2,n1:n2))*&
            (dpp_dy - rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2)*duy_dy)
      case('top')
        L_1 = KK*(cs20_ar(l1:l2,n1:n2)/&
            gamma0(l1:l2,n1:n2)*rho0(l1:l2,n1:n2)-p_infy)
        L_5 = (f(l1:l2,mmm,n1:n2,iuy) + cs0_ar(l1:l2,n1:n2))*&
            (dpp_dy + rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2)*duy_dy)
      endselect
      L_2 = f(l1:l2,mmm,n1:n2,iuy)*(cs20_ar(l1:l2,n1:n2)*dlnrho_dy-dpp_dy)
!
      if (ldensity_nolog) then
        df(l1:l2,mmm,n1:n2,ilnrho) = &
            drho_prefac*(L_2+0.5*(L_5 + L_1))-dmom2_dx(l1:l2,n1:n2)
      else
        df(l1:l2,mmm,n1:n2,ilnrho) = drho_prefac*(L_2+0.5*(L_5 + L_1)) &
            -1./rho0(l1:l2,n1:n2)*dmom2_dx(l1:l2,n1:n2)
      endif
!
      df(l1:l2,mmm,n1:n2,iuy) = -1./&
          (2.*rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2))*(L_5 - L_1) &
          -f(l1:l2,mmm,n1:n2,iuy)*duy_dx(l1:l2,n1:n2)
      df(l1:l2,mmm,n1:n2,ilnTT) = -1./&
          (rho0(l1:l2,n1:n2)*cs20_ar(l1:l2,n1:n2))*(-L_2 &
          +0.5*(gamma0(l1:l2,n1:n2)-1.)*(L_5+L_1)) &
          -1./(rho0(l1:l2,n1:n2)*cs20_ar(l1:l2,n1:n2))*&
          (gamma0(l1:l2,n1:n2)-1.) &
          *gamma0(l1:l2,n1:n2)*drhoE_p_dx(l1:l2,n1:n2)
!
    endsubroutine bc_nscbc_nref_subout_y
!***********************************************************************
    subroutine bc_ADI_flux_z(f,topbot)
!
!  constant flux boundary condition for temperature (called when bcz='c3')
!  30-jan-2009/dintrans: coded 
!
      use SharedVariables, only: get_shared_variable
!
      real, pointer :: hcond(:),Fbot
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: tmp_x
      integer :: i,ierr
!
!  Do the `c3' boundary condition (constant heat flux) for TT
!  at the bottom _only_ in the ADI case where hcond is variable.
!  TT version: enforce dT/dz = - Fbot/K
!      
      call get_shared_variable('hcond0',hcond,ierr)
      if (ierr/=0) call stop_it("bc_lnTT_flux_z: "//&
           "there was a problem when getting hcond")      
      call get_shared_variable('Fbot',Fbot,ierr)
      if (ierr/=0) call stop_it("bc_lnTT_flux_z: "//&
           "there was a problem when getting Fbot")      
 
      if(headtt) print*,'bc_ADI_flux_z: Fbot,hcond,dz=',Fbot,hcond,dz

      if (topbot.eq.'bot') then
        tmp_x=-Fbot/hcond
        do i=1,nghost
          f(:,4,n1-i,ilnTT)=f(:,4,n1+i,ilnTT)-2.*i*dz*tmp_x
        enddo
      else
        call fatal_error('bc_ADI_flux_z','invalid argument')
      endif
!
    endsubroutine bc_ADI_flux_z
!***********************************************************************
endmodule Boundcond
