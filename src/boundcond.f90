! $Id: boundcond.f90,v 1.191 2007-11-20 21:21:22 wlyra Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   boundcond.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!

!!!  Module for boundary conditions. Extracted from (no)mpicomm, since
!!!  all non-periodic (external) boundary conditions require the same
!!!  code for serial and parallel runs.

module Boundcond

  use Mpicomm
  use Messages

  implicit none

  private

  public :: boundconds, boundconds_x, boundconds_y, boundconds_z
  public :: bc_per_x, bc_per_y, bc_per_z
  public :: update_ghosts

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
      use Cdata
      use Entropy
      use EquationOfState
      use Magnetic
      use Radiation
      use Shear
      use Special, only: special_boundconds
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      real, dimension (mcom) :: fbcx12
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
          call boundcond_shear(f,ivar1,ivar2)
        else
          do k=1,2                ! loop over 'bot','top'
            if (k==1) then
              topbot='bot'; bc12=bcx1; fbcx12=fbcx1; ip_ok=0
            else
              topbot='top'; bc12=bcx2; fbcx12=fbcx2; ip_ok=nprocx-1
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
                  if (j==iss)   call bc_ss_flux_x(f,topbot,FbotKbot)
                  if (j==ilnTT) call bc_lnTT_flux_x(f,topbot,hcond0,hcond1,Fbot)
                case ('sT')
                  ! BCX_DOC: symmetric temperature, $T_{N-i}=T_{N+i}$;
                  ! BCX_DOC: implies $T'(x_N)=T'''(x_0)=0$
                  if (j==iss) call bc_ss_stemp_x(f,topbot)
                case ('in')
                  ! BCX_DOC: inflow [Need details!]
                  if (j==ie) call bc_ee_inflow_x(f,topbot)
                case ('out')
                  ! BCX_DOC: outflow [Need details!]
                  if (j==ie) call bc_ee_outflow_x(f,topbot)
                case ('db')
                  ! BCX_DOC:
                  call bc_db_x(f,topbot,j)
                case ('f')
                  ! BCX_DOC: ``freeze'' value, i.e. maintain initial
                  !  value at boundary
                  call bc_freeze_var_x(topbot,j)
                  call bc_sym_x(f,-1,topbot,j,REL=.true.) 
                  ! antisymm wrt boundary
                case ('1')
                  ! BCX_DOC: $f=1$ (for debugging)
                  call bc_one_x(f,topbot,j)
                case ('set')
                  ! BCX_DOC: set boundary value to |fbcx12|
                  call bc_sym_x(f,-1,topbot,j,REL=.true.,val=fbcx12)
                case ('der')
                  ! BCX_DOC: set derivative on boundary to |fbcx12|
                  call bc_set_der_x(f,topbot,j,fbcx12(j))
                case ('slo')
                  ! BCX_DOC: set slope at the boundary = fbcx12
                  call bc_slope_x(f,fbcx12,topbot,j)
                case ('dr0')
                  ! BCX_DOC: set boundary value [really??]
                  call bc_dr0_x(f,fbcx12,topbot,j)
                case ('ovr')
                  ! BCX_DOC: overshoot boundary condition
                  ! BCX_DOC:  ie (d/dx-1/dist) f = 0. 
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
                case ('spd')
                  ! BCX_DOC: Only for spherical polar coordinate system, 
                  ! BCX_DOC:  sets  d(rA_{\alpha})/dr = fbcx12(j)
                  call bc_set_spder_x(f,topbot,j,fbcx12(j))
                case('sfr')
                  ! BCX_DOC: stress-free boundary condition for spherical
                  ! BCX_DOC: coordinates.
                  ! BCX_DOC: $d_r(u_{\theta}) = u_{\theta}/r$  with $u_r = 0$
                  ! BCX_DOC: sets $S_{r \theta}$ component of the strain
                  ! BCX_DOC: matrix to be zero in spherical coordinate  
                  ! BCX_DOC: system.
                  ! BCX_DOC: This subroutine sets only the first part of
                  ! BCX_DOC: this boundary condition for 'j'-th component
                  ! BCX_DOC: of f. 
                  ! BCX_DOC: \nFIXME: this description is completely
                  ! BCX_DOC: incomprehensible [wd, 12-nov-2007] 
                  call bc_set_sfree_x(f,topbot,j)
                case('pfc')
                  ! BCX_DOC: In spherical polar coordinate system,
                  ! BCX_DOC: at a radial boundary set : $A_{\theta} = 0$ and 
                  ! BCX_DOC: $A_{phi} = 0$, and demand $div A = 0$ gives the 
                  ! BCX_DOC: condition on $A_r$ to be  $d/dr( A_r) + 2/r = 0$ . 
                  ! BCX_DOC: This subroutine sets this condition on  $j$-th 
                  ! BCX_DOC: component of f. As this is related to setting the
                  ! BCX_DOC: perfect conducting boundary condition we call 
                  ! BCX_DOC: this "pfc".  
                  call bc_set_pfc_x(f,topbot,j)
                 case ('fix')
                  ! BCX_DOC: set boundary value [really??]
                  call bc_fix_x(f,topbot,j,fbcx12(j))
                case('cfb')
                  ! BCZ_DOC: radial centrifugal balance 
                  if (lcylindrical_coords) then
                    call bc_lnrho_cfb_r_iso(f,topbot,j)
                  else
                    print*,'not implemented for other than cylindrical'
                    stop
                  endif
                case ('ouf')
                  ! BCX_DOC: allow outflow, but no inflow (experimental)
                  call bc_outflow_x(f,topbot,j)
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
      use Cdata
      use Entropy
      use Magnetic
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
              case ('der')
                ! BCY_DOC: set derivative on the boundary
                call bc_set_der_y(f,topbot,j,fbcy12(j))
              case('sfr')
                  ! BCY_DOC: stress-free boundary condition for spherical
                  ! BCY_DOC: coordinates.
                  ! BCY_DOC: $(1/r)d_\theta(u_{\phi}) =\cot(\theta) u_{\phi}$
                  ! BCY_DOC: with $u_{\theta} = 0$ sets $S_{\phi \theta}$
                  ! BCY_DOC: component of the strain matrix to be zero in
                  ! BCY_DOC: spherical coordinate system.
                  ! BCY_DOC: This subroutine sets only the first part of this
                  ! BCY_DOC: boundary condition for 'j'-th component of f. 
                  ! BCX_DOC: \nFIXME: this description is completely
                  ! BCX_DOC: incomprehensible [wd, 12-nov-2007] 
                  call bc_set_sfree_y(f,topbot,j)
              case('pfc')
                  ! BCY_DOC: In spherical polar coordinate system,
                  ! BCY_DOC: at a theta boundary set : $A_r = 0$ and 
                  ! BCY_DOC: $A_{\phi} = 0$, and demand $div A = 0$ gives the 
                  ! BCY_DOC: condition on $A_r$ to be  
                  ! BCY_DOC: $d/d\theta( A_\theta) + cot(\theta)A_{\theta} = 0$.
                  ! BCY_DOC: This subroutine sets this condition on  $j$-th 
                  ! BCY_DOC: component of f. As this is related to setting the
                  ! BCY_DOC: perfect conducting boundary condition we call 
                  ! BCY_DOC: this "pfc".  
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
      use Cdata
      use Entropy, only: hcond0,hcond1,Fbot,FbotKbot,Ftop,FtopKtop,chi, &
                         lmultilayer,lheatc_chiconst
      use Magnetic
      use Special, only: special_boundconds
      !use Density
      use EquationOfState
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      real, dimension (mcom) :: fbcz12, fbcz12_1, fbcz12_2
      real :: Ftopbot,FtopbotK
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
      select case(nzgrid)
!
      case(1)
        if (headtt) print*,'boundconds_z: no z-boundary'
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
            ip_ok=0
            Ftopbot=Fbot
            FtopbotK=FbotKbot
          else
            topbot='top'
            bc12=bcz2
            fbcz12=fbcz2
            fbcz12_1=fbcz2_1
            fbcz12_2=fbcz2_2
            ip_ok=nprocz-1
            Ftopbot=Ftop
            FtopbotK=FtopKtop
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
                if (j==iss) call bc_ss_flux(f,topbot,hcond0,hcond1,Ftopbot,FtopbotK,chi, &
                                  lmultilayer,lheatc_chiconst)
                if (j==iaa) call bc_aa_pot(f,topbot)
                if (j==ilnTT) call bc_lnTT_flux_z(f,topbot,hcond0,Fbot)
              case ('pot')
                ! BCZ_DOC: 
                if (j==iaa) call bc_aa_pot2(f,topbot)
              case ('pwd')
                ! BCZ_DOC: 
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
              case ('fB')
                ! BCZ_DOC: frozen-in B-field
                ! tell other modules not to change boundary value
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
    subroutine bc_per_x(f,topbot,j)
!
!  periodic boundary condition
!  11-nov-02/wolf: coded
!
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
    subroutine bc_set_spder_x(f,topbot,j,val)
!
!  Sets the derivative, particularly: 
!    d(rA_{\alpha})/dr = <val>
!  on the boundary to a given value
!
!  27-apr-2007/dhruba: coded
!
      use Cdata
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j

      real, intent (in) :: val
      integer :: i

      if(lspherical_coords)then
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
      use Cdata
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
    subroutine bc_set_sfree_x(f,topbot,j)
! "stress-free" boundary condition for spherical coordinate system. 
! d_r(u_{\theta}) = u_{\theta}/r  with u_r = 0 sets S_{r \theta}
! component of the strain matrix to be zero in spherical coordinate system. 
! This subroutine sets only the first part of this boundary condition for 'j'-th
! component of f. 
!
!  25-Aug-2007/dhruba: coded
!
      use Cdata
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j


      select case(topbot)

      case('bot')               ! bottom boundary
! The coding assumes we are using 6-th order centered finite difference for our
! derivatives. 
        f(l1-1,:,:,j)= f(l1+1,:,:,j) -  60.*f(l1,:,:,j)*dx/(45.*x(l1))
        f(l1-2,:,:,j)= f(l1+2,:,:,j) -  60.*f(l1,:,:,j)*dx/(9.*x(l1))
        f(l1-3,:,:,j)= f(l1+3,:,:,j) -  60.*f(l1,:,:,j)*dx/(x(l1))
      case('top')               ! top boundary
        f(l2+1,:,:,j)= f(l2-1,:,:,j) +  60.*f(l1,:,:,j)*dx/(45.*x(l2))
        f(l2+2,:,:,j)= f(l2-2,:,:,j) +  60.*f(l1,:,:,j)*dx/(9.*x(l2))
        f(l2+3,:,:,j)= f(l2-3,:,:,j) +  60.*f(l1,:,:,j)*dx/(x(l2))

      case default
        call warning('bc_set_sfree_x',topbot//" should be `top' or `bot'")

      endselect
!
    endsubroutine bc_set_sfree_x
! **********************************************************************
    subroutine bc_set_sfree_y(f,topbot,j)
! "stress-free" boundary condition for spherical coordinate system. 
! d_r(u_{\theta}) = u_{\theta}/r  with u_r = 0 sets S_{r \theta}
! component of the strain matrix to be zero in spherical coordinate system. 
! This subroutine sets only the first part of this boundary condition for 'j'-th
! component of f. 
!
!  25-Aug-2007/dhruba: coded
!
      use Cdata
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
      real :: cottheta

      select case(topbot)

      case('bot')               ! bottom boundary
! The coding assumes we are using 6-th order centered finite difference for our
! derivatives. 
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

      endselect
!
    endsubroutine bc_set_sfree_y
! **********************************************************************
    subroutine bc_set_pfc_y(f,topbot,j)
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
      use Cdata
!
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer, intent (in) :: j
      real :: cottheta

      select case(topbot)

      case('bot')               ! bottom boundary
! The coding assumes we are using 6-th order centered finite difference for our
! derivatives. 
        cottheta= cotth(m1)
        f(:,m1-1,:,j)= f(:,m1+1,:,j) +  2.*60.*dy*cottheta*f(:,m1,:,j)/45.
        f(:,m1-2,:,j)= f(:,m1+2,:,j) +  2.*60.*dy*cottheta*f(:,m1,:,j)/9.
        f(:,m1-3,:,j)= f(:,m1+3,:,j) +  2.*60.*dy*cottheta*f(:,m1,:,j)
      case('top')               ! top boundary
        cottheta= cotth(m2)
        f(:,m2+1,:,j)= f(:,m2-1,:,j) -  2.*60.*dy*cottheta*f(:,m2,:,j)/45.
        f(:,m2+2,:,j)= f(:,m2-2,:,j) -  2.*60.*dy*cottheta*f(:,m2,:,j)/9.
        f(:,m2+3,:,j)= f(:,m2-3,:,j) -  2.*60.*dy*cottheta*f(:,m2,:,j)

      case default
        call warning('bc_set_pfc_y',topbot//" should be `top' or `bot'")

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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
    use Cdata
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
      use Cdata
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
    subroutine bc_onesided_z(f,topbot,j)
!
!  One-sided conditions.
!  These expressions result from combining Eqs(207)-(210), astro-ph/0109497,
!  corresponding to (9.207)-(9.210) in Ferriz-Mas proceedings.
!
!   5-apr-03/axel: coded
!
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
    subroutine bc_extrap0_2_0(f,topbot,j)
!
!  Extrapolation boundary condition for f(bdry)=0.
!  Correct for polynomials up to 2nd order, determined no further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!    9-oct-03/wolf: coded
!
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
!***********************************************************************
!    subroutine bc_force_kepler(f,idz,j)
!!
!! DOCUMENT ME!!
!!
!      use Cdata, only: x,y,iux,iuy,iuz,pi, &
!                       mx,my,mz,m1,m2,l1,l2,mfarray,nx
!      use Mpicomm, only: stop_it
!      use Hydro, only: kep_cutoff_pos_ext,kep_cutoff_width_ext
!      use Hydro, only: kep_cutoff_pos_int,kep_cutoff_width_int
!      use Hydro, only: u_out_kep
!      use FArrayManager
!
!      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!      integer, intent(in) :: idz,j
!      real, dimension (nx,3) :: gg
!      real, dimension (nx) :: r,g_r
!      real, dimension (mx,my), save :: ux,uy,uz
!      real :: r1_ext,r2_ext
!      real :: r1_int,r2_int
!      integer :: m
!      logical, save :: initialize=.true.
!
!      if (initialize) then
!
!        r1_int=kep_cutoff_pos_int-kep_cutoff_width_int/2
!        r2_int=kep_cutoff_pos_int+kep_cutoff_width_int/2
!
!        r1_ext=kep_cutoff_pos_ext-kep_cutoff_width_ext/2
!        r2_ext=kep_cutoff_pos_ext+kep_cutoff_width_ext/2
!
!        do m=m1,m2
!
!          r=sqrt(x(l1:l2)**2+y(m)**2)
!
!          !
!          ! First time round, look up gg in the f-array
!          ! and get it's index.
!          !
!          if (.not.associated(iglobal_gg)) then
!            call farray_use_global('gg',iglobal_gg,vector=3)
!            if (.not.associated(iglobal_gg)) then
!              call fatal_error("bc_force_kepler", &
!               "Could not get global gg from an f-array slot")
!            endif
!          endif
!
!          gg=f(l1:l2,m,idz,iglobal_gg:iglobal_gg+2)
!
!          g_r=sqrt(gg(:,1)**2+gg(:,2)**2)
!
!          ux(l1:l2,m)=-y(  m  )*sqrt(g_r/(r+10*tiny(r)))
!          uy(l1:l2,m)= x(l1:l2)*sqrt(g_r/(r+10*tiny(r)))
!
!          where (r>r1_int.and.r<=r2_int)
!            ux(l1:l2,m)=ux(l1:l2,m)*sin((pi/2)*(r-r1_int)/(r2_int-r1_int))**2
!            uy(l1:l2,m)=uy(l1:l2,m)*sin((pi/2)*(r-r1_int)/(r2_int-r1_int))**2
!          endwhere
!
!          where (r>r1_ext.and.r<=r2_ext)
!            ux(l1:l2,m)=ux(l1:l2,m)*cos((pi/2)*(r-r1_ext)/(r2_ext-r1_ext))**2
!            uy(l1:l2,m)=uy(l1:l2,m)*cos((pi/2)*(r-r1_ext)/(r2_ext-r1_ext))**2
!          endwhere
!
!          where (r<=r1_int.or.r>r2_ext)
!            ux(l1:l2,m)=0
!            uy(l1:l2,m)=0
!          endwhere
!
!          ! outflow velocity = u_out_kep * u_kep
!          uz(l1:l2,m)=u_out_kep*sqrt(ux(l1:l2,m)**2+uy(l1:l2,m)**2)
!
!          initialize=.false.
!
!        enddo
!
!      endif
!
!      if     (j==iux) then
!        f(l1:l2,m1:m2,idz,j)=ux(l1:l2,m1:m2)
!      elseif (j==iuy) then
!        f(l1:l2,m1:m2,idz,j)=uy(l1:l2,m1:m2)
!      elseif (j==iuz) then
!        f(l1:l2,m1:m2,idz,j)=uz(l1:l2,m1:m2)
!      else
!        call stop_it('BC_FORCE_KEPLER: only implemented for uu')
!      endif
!
!    endsubroutine bc_force_kepler
!!***********************************************************************
    subroutine bc_one_x(f,topbot,j)
!
!  Set bdry values to 1 for debugging purposes
!
!  11-jul-02/wolf: coded
!
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
      use Cdata
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
!
       use Cdata
       use EquationOfState, only : gamma,gamma1,gamma11,cs20,lnrho0

       real, dimension (mx,my,mz,mfarray) :: f
       real, dimension (nxgrid,nygrid),save :: uxl,uxr,uyl,uyr
       real, dimension (nxgrid,nygrid) :: uxd,uyd
       real, dimension (nx,ny) :: quen,pp,betaq,fac
       real, dimension (nx,ny) :: bbz,bb2
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
!   First get Bz component:
!
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

       bb2 = bbz*bbz
       bb2 = bb2/(2*mu0)*300.
!
       if (ltemperature) pp = gamma1*gamma11*exp(f(l1:l2,m1:m2,n1,ilnrho)+f(l1:l2,m1:m2,n1,ilnTT))
!
       if (lentropy) then
          if (pretend_lnTT) then
             pp = gamma1*gamma11*exp(f(l1:l2,m1:m2,n1,ilnrho)+f(l1:l2,m1:m2,n1,iss))
          else
             pp = gamma* (f(l1:l2,m1:m2,n1,iss)+f(l1:l2,m1:m2,n1,ilnrho))-gamma1*lnrho0
             pp = exp(pp) * cs20*gamma11
          endif
       endif
!
!   limit plasma beta
!
       where (bb2 .gt. sqrt(tini))
          betaq =  pp / bb2
       elsewhere
          betaq = pp * sqrt(tini)
       endwhere
!
       quen=(1.+betaq**2)/(3.+betaq**2)
!
!   Fill the ghost cells and the bottom layer with velocity field
!
       do j=1,n1
         f(l1:l2,m1:m2,j,iux)=uxd(ipx*nx+1:(ipx+1)*nx,ipy*ny+1:(ipy+1)*ny)*quen
         f(l1:l2,m1:m2,j,iuy)=uyd(ipx*nx+1:(ipx+1)*nx,ipy*ny+1:(ipy+1)*ny)*quen
       enddo
!
     endsubroutine uu_driver
!***********************************************************************
    subroutine bc_lnTT_flux_x(f,topbot,hcond0,hcond1,Fbot)
!
!  constant flux boundary condition for temperature (called when bcx='c1')
!  12-Mar-2007/dintrans: coded
!
      use Cdata
!
      real, intent(in) :: hcond0, hcond1, Fbot
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (my,mz) :: tmp_yz
      integer :: i
!
!  Do the `c1' boundary condition (constant heat flux) for lnTT.
!  check whether we want to do top or bottom (this is precessor dependent)
!
      if(headtt) print*,'bc_lnTT_flux_x: Fbot,hcond,dx=',Fbot,hcond0*hcond1,dx

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
    subroutine bc_lnTT_flux_z(f,topbot,hcond0,Fbot)
!
!  constant flux boundary condition for temperature (called when bcz='c1')
!  12-May-07/dintrans: coded
!
      use Cdata
!
      real, intent(in) :: hcond0, Fbot
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: tmp_xy
      integer :: i
!
!  Do the `c1' boundary condition (constant heat flux) for lnTT or TT (if
!  ltemperature_nolog=.true.) at the bottom _only_.
!  lnTT version: enforce dlnT/dz = - Fbot/(K*T)
!    TT version: enforce   dT/dz = - Fbot/K
!
      if(headtt) print*,'bc_lnTT_flux_z: Fbot,hcond,dz=',Fbot,hcond0,dz

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
    subroutine bc_ss_flux_x(f,topbot,FbotKbot)
!
!  constant flux boundary condition for entropy (called when bcx='c1')
!  17-mar-07/dintrans: coded
!
      use Cdata
      use EquationOfState, only: gamma, gamma1, lnrho0, cs20
!
      real, intent(in) :: FbotKbot
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (my,mz) :: tmp_yz,cs2_yz
      integer :: i
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!  ===============
!
      case('bot')
        if(headtt) print*,'bc_ss_flux_x: FbotKbot=',FbotKbot
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
!  Pontential field boundary condition
!
!  11-oct-06/wolf: Adapted from Tobi's bc_aa_pot2
!
      use Cdata
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
      use Cdata, only: mx, my, mz, mfarray, n1, n2, nghost
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
    subroutine bc_outflow_x(f,topbot,j)
!
!  Outflow boundary conditions.
!
!  If velocity vector points out of the box, the velocity in the
!  ghost cells is copied from the last grid point.
!
!  For inwards pointing velocity vector, the velocity is set to zero
!  in the ghost cells.
!
!  12-aug-2007/anders: implemented
!
      use Cdata, only: mx, my, mz, mfarray, l1, l2, nghost
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      integer :: i, iz, iy
!
      select case(topbot)
!
!  Bottom boundary.
!
      case('bot')
        do iz=1,mz; do iy=1,my; 
          if (f(l1,iy,iz,j)<=0.0) then  ! simply copy to ghost cells
            do i=1,nghost; f(l1-i,iy,iz,j)=f(l1,iy,iz,j); enddo
          else                          ! zero, suppressing inflow
            do i=1,nghost; f(l1-i,iy,iz,j)=0.0; enddo
          endif
        enddo; enddo
!
!  Top boundary.
!
      case('top')
        do iz=1,mz;do iy=1,my; 
          if (f(l2,iy,iz,j)>=0.0) then
            do i=1,nghost; f(l2+i,iy,iz,j)=f(l2,iy,iz,j); enddo
          else
            do i=1,nghost; f(l2+i,iy,iz,j)=0.0; enddo
          endif
        enddo; enddo
!
!  Default.
!
      case default
        print*, "bc_outflow_x: ", topbot, " should be `top' or `bot'"
!
      endselect
!
    endsubroutine bc_outflow_x
!***********************************************************************
    subroutine bc_outflow_z(f,topbot,j)
!
!  Outflow boundary conditions.
!
!  If velocity vector points out of the box, the velocity in the
!  ghost cells is copied from the last grid point.
!
!  For inwards pointing velocity vector, the velocity is set to zero
!  in the ghost cells.
!
!  12-aug-2007/anders: implemented
!
      use Cdata, only: mx, my, mz, mfarray, n1, n2, nghost
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
          if (f(ix,iy,n1,j)<=0.0) then  ! simply copy to ghost cells
            do i=1,nghost; f(ix,iy,n1-i,j)=f(ix,iy,n1,j); enddo
          else                          ! zero, suppressing inflow
            do i=1,nghost; f(ix,iy,n1-i,j)=0.0; enddo
          endif
        enddo; enddo
!
!  Top boundary.
!
      case('top')
        do iy=1,my; do ix=1,mx
          if (f(ix,iy,n2,j)>=0.0) then
            do i=1,nghost; f(ix,iy,n2+i,j)=f(ix,iy,n2,j); enddo
          else
            do i=1,nghost; f(ix,iy,n2+i,j)=0.0; enddo
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
      use Cdata
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
endmodule Boundcond
