!***********************************************************************
    subroutine boundconds(f,ivar1_opt,ivar2_opt)
!
!  Apply boundary conditions in all three directions.
!  Note that we _must_ call boundconds_{x,y,z} in this order, or edges and
!  corners will not be OK.
!
!  10-oct-02/wolf: coded
!
      real, dimension (:,:,:,:) :: f
      integer, optional :: ivar1_opt, ivar2_opt
      integer :: ivar1, ivar2
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
      call boundconds_x(f,ivar1,ivar2)
      call boundconds_y(f,ivar1,ivar2)
      call boundconds_z(f,ivar1,ivar2)
!
    endsubroutine boundconds
!***********************************************************************
    subroutine boundconds_x_c(f,ivar1_opt,ivar2_opt)
!
!  Envelope for being called from C code.
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt

      call boundconds_x(f,ivar1_opt,ivar2_opt)

    endsubroutine boundconds_x_c
!***********************************************************************
    subroutine check_boundconds_x(f,ivar1_opt,ivar2_opt)
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
!  30-sep-16/MR: new BCs 'n1s' = Neumann + 1-sided, "c1s" = heatflux + 1-sided (only for z boundaries)
!
      use General, only: itoa
      use Shear
      use Special, only: special_boundconds
!
      real, dimension (:,:,:,:) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, j, topbot
      type (boundary_condition) :: bc
      character(LEN=3) :: bc_code, cjvar 
!
      if (nghost<=0 .or. nxgrid<=1) return

      ivar1=1; ivar2=min(mcom,size(f,4))
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  Use the following construct to keep compiler from complaining if
!  we have no variables (and boundconds) at all (samples/no-modules):
!
      if (all(bcx12(ivar1:ivar2,:)=='she')) then
        !call boundcond_shear(f,ivar1,ivar2)
      else
        do topbot=BOT,TOP
!
          do j=ivar1,ivar2

            cjvar=itoa(j)
            bc_code=bcx12(j,topbot)
!
            if (bc_code == 'she') then
              if (bcx12(j,1) /= bcx12(j,2)) then
                errmsg='generalize '//trim(bc_code)//' to have sheared periodic boundary on only one end'; goto 10
              endif
            else
              select case (bc_code)
              case ('cT')
                if (.not.(j==iss .or. j==ilnTT)) then; errmsg=' not allowed for variable no. '//trim(cjvar)); goto 10; endif
              case ('c1')
                if (.not.(j==iss .or. j==ilnTT)) then; errmsg=' not allowed for variable no. '//trim(cjvar)); goto 10; endif
              case ('Fgs')
                if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 10; endif
              case ('Fct')
                if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 10; endif
              case ('Fcm')
                if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 10; endif
              case ('sT')
                if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 10; endif
              case ('asT')
                if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 10; endif
              case ('hat')
                if (topbot==TOP) then; errmsg=' only defined at x-bottom'; goto 10; endif
              case ('jet')
                if (topbot==TOP) then; errmsg=' only defined at x-bottom'; goto 10; endif
              case ('sfr')
                if (.not.lspherical_coords) then; errmsg=' only defined for spherical coordinates'); goto 10; endif
                ! BCX_DOC: for spherical coordinate system.
              case ('sr1')
                if (.not.lspherical_coords) then; errmsg=' only defined for spherical coordinates'); goto 10; endif
                ! BCX_DOC: Stress-free bc for spherical coordinate system.
              case ('nfr')
                if (.not.lspherical_coords) then; errmsg=' only defined for spherical coordinates'); goto 10; endif
                ! BCX_DOC: Normal-field bc for spherical coordinate system.
              case ('nr1')
                if (.not.lspherical_coords) then; errmsg=' only defined for spherical coordinates'); goto 10; endif
                ! BCX_DOC: Normal-field bc for spherical coordinate system.
              case ('sa2')
                ! BCX_DOC: $(d/dr)(r B_{\phi}) = 0$ imposes
                ! BCX_DOC: boundary condition on 2nd derivative of
                ! BCX_DOC: $r A_{\phi}$. Same applies to $\theta$ component.
                !call bc_set_sa2_x(f,topbot,j)
              case ('pfc')
                if (.not.lspherical_coords) then; errmsg=' only defined for spherical coordinates'); goto 10; endif
                ! BCX_DOC: perfect-conductor in spherical
              case ('cfb')
                ! BCX_DOC: radial centrifugal balance
                !call bc_lnrho_cfb_r_iso(f,topbot)
              case ('ioc')
                if (.not.lcylindrical_coords) then; errmsg=' only defined for spherical coordinates'); goto 10; endif
                ! BCX_DOC: in cylindrical coordinates
                !call bc_inlet_outlet_cyl(f,topbot,j,fbcx(j,topbot))
              case default
                bc%bcname=bc_code
                bc%ivar=j
                bc%location=(((topbot-1)*2)-1)   ! -1/1 for x bot/top
                bc%value1=fbcx(j,topbot)
                bc%value2=fbcx(j,topbot)
                bc%done=.false.
!
                call special_boundconds(f,bc)
!
                if (.not.bc%done) then; errmsg=" no such boundary condition for j="//trim(cjvar); goto 10; endif
              endselect
  10          call fatal_error('check_bounconds_x',trim(bc_code)//errmsg)
            endif
          enddo
        enddo
      endif
!
    endsubroutine check_boundconds_x
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
      use General, only: itoa
      use Special, only: special_boundconds
!
      real, dimension (:,:,:,:) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, j, topbot
      type (boundary_condition) :: bc
      character(LEN=3) :: bc_code, cjvar 
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
      if (nghost<=0.or.nygrid<=1) return
!
!  Boundary conditions in y
!
      do topbot=BOT,TOP              ! loop over 'bot','top'
!
        do j=ivar1,ivar2
!
          select case (bcy12(j,topbot))
          case ('pp')
            if (.not.lspherical_coords) then; errmsg=' only defined for spherical coordinates'); goto 20; endif
            ! BCY_DOC: periodic across the pole
            call bc_pper_y(f,topbot,j,+1)
          case ('yy')
            if (.not.lspherical_coords) then; errmsg=' only defined for spherical coordinates'); goto 20; endif
            ! BCY_DOC: Yin-Yang grid
            call bc_yy_y(f,topbot,j)
          case ('ap')
            if (.not.lspherical_coords) then; errmsg=' only defined for spherical coordinates'); goto 20; endif
            ! BCY_DOC: anti-periodic across the pole
            call bc_pper_y(f,topbot,j,-1)
          case ('cT')
            ! BCY_DOC: constant temp.
            if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 20; endif
          case ('sT')
            ! BCY_DOC: symmetric temp.
            if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 20; endif
          case ('asT')
            ! BCY_DOC: select entropy for uniform ghost temperature
            ! BCY_DOC: matching fluctuating boundary value,
            ! BCY_DOC: $T_{N-i}=T_{N}=$;
            ! BCY_DOC: implies $T'(x_N)=T'(x_0)=0$
            if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 20; endif
          case ('sfr')
            ! BCY_DOC: stress-free boundary condition for spherical
            ! BCY_DOC: coordinate system.
            if (.not.lspherical_coords) then; errmsg=' only defined for spherical coordinates'); goto 20; endif
            if (j==iux.or.j==iuy) then; errmsg='stress-free BC at theta boundary only allowed for uu_phi'; goto 20; endif
          case ('nfr')
            ! BCY_DOC: Normal-field bc for spherical coordinate system.
            if (.not.lspherical_coords) then; errmsg=' only defined for spherical coordinates'); goto 20; endif
            !call bc_set_nfr_y(f,topbot,j)
          case ('spt')
            ! BCY_DOC: spherical perfect conducting boundary condition
            if (.not.lspherical_coords) then; errmsg=' only defined for spherical coordinates'); goto 20; endif
            !call bc_spt_y(f,topbot,j)
          case ('pfc')
            ! BCY_DOC: perfect conducting boundary condition
            ! BCY_DOC: along $\theta$ boundary
            if (.not.lspherical_coords) then; errmsg=' only defined for spherical coordinates'); goto 20; endif
!joern: WARNING, this bc will NOT give a perfect-conductor boundary condition
            !call bc_set_pfc_y(f,topbot,j)
          case default
            bc%bcname=bcy12(j,topbot)
            bc%ivar=j
            bc%value1=fbcy(j,topbot)
            bc%value2=fbcy(j,topbot)
            bc%location=(((topbot-1)*4)-2)   ! -2/2 for y bot/top
            bc%done=.false.
!
            if (lspecial) call special_boundconds(f,bc)
!
            if (.not.bc%done) then; errmsg=" no such boundary condition for j="//trim(cjvar); goto 20; endif
          endselect
  20      call fatal_error('check_bounconds_x',trim(bc_code)//errmsg)
        enddo
      enddo
!
    endsubroutine check_boundconds_y
!***********************************************************************
    subroutine check_boundconds_z(f,ivar1_opt,ivar2_opt)
!
!  Boundary conditions in z, except for periodic part handled by communication.
!  Remark: boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-'corners').
!
!   8-jul-02/axel: split up into different routines for x,y and z directions
!  11-nov-02/wolf: unified bot/top, now handled by loop
!  02-apr-13/MR  : added new boundary condition 'fs' = frozen boundary value
!                  + symmetry about boundary; added 'fa' for alternative reference to
!                  already existing freezing condition (includes antisymmetry)
!  30-dec-16/MR: added BC 'a1s' for constant alpha mean-field model in one dimension
!
      use General, only: itoa
      use Gravity, only: gravz_profile
      use Special, only: special_boundconds
      use EquationOfState
!
      real, dimension (:,:,:,:) :: f
      integer, optional :: ivar1_opt, ivar2_opt
      integer :: ivar1, ivar2, j, topbot
      type (boundary_condition) :: bc
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
      if (nghost<=0.or.nzgrid<=1) return
!
!  Boundary conditions in z
!
      do topbot=BOT,TOP                ! loop over 'bot','top'
!
        do j=ivar1,ivar2

          select case (bcz12(j,topbot))
          case ('yy')
            ! BCZ_DOC: Yin-Yang grid
            if (.not.lspherical_coords) then; errmsg=' only defined for spherical coordinates'); goto 30; endif
          case ('a1s')
            ! BCZ_DOC: special for perfect conductor with const alpha and etaT when A considered as B; one-sided for 1st and 2nd derivative in two first inner points
            !call pc_aasb_const_alpha(f,topbot,j)
            !call set_ghosts_for_onesided_ders(f,topbot,j,3,.true.)
          case ('c1')
            if (.not.(j==iss.or.j==iaa.or.j==ilnTT)) then 
              errmsg=' not allowed for variable no. '//trim(cjvar)
              goto 30
            endif
          case ('c1s')
            if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 30; endif
          case ('Fgs')
            if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 30; endif
          case ('Fct')
            if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 30; endif
          case ('c3')
            if (j/=ilnTT) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 30; endif
          case ('pfe')
            if (j/=iaa) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 30; endif
          case ('p1D')
            if (j/=iay) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 30; endif
          case ('pot')
            if (j/=iaa) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 30; endif
          case ('pwd')
            if (j/=iaa) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 30; endif
          case ('hds')
            ! BCZ_DOC: hydrostatic equilibrium with a high-frequency filter
            !call bc_lnrho_hdss_z_iso(f,topbot)
          case ('cT')
            if (.not.(j==ilnrho.or.j==iss.or.j==iTT.or.j==ilnTT)) then
              errmsg=' not allowed for variable no. '//trim(cjvar)
              goto 30
            endif
          case ('cT1')
            ! BCZ_DOC: constant temperature using one-sided derivatives
            !call bc_ss_temp_z(f,topbot,.true.)
          case ('cT2')
            if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 30; endif
          case ('cT3')
            if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 30; endif
          case ('hs')
            if (.not. lgrav) then 
              errmsg=' requires gravity'
              goto 30
            endif
            if ((.not. ltemperature .or. ltemperature_nolog) .and. (gravz_profile /= 'const')) then
              errmsg=' requires a constant gravity profile'
              goto 30
            endif
            if (.not. lequidist(3)) then
              errmsg=' requires symmetric grid distances on the z boundary')
              goto 30
            endif
            if (.not.((j==ilnrho) .or. (j==irho_b) .or. (j==iss) .or. (j==ipp))) then
              errmsg=" requires density or pressure"
              goto 30
            endif
          case ('hse')
            if (.not. lgrav) then; errmsg="'hse' requires gravity"; goto 30; endif
            if (.not. leos)  then; errmsg="'hse' requires an eos module"; goto 30; endif
            if ((ilnrho == 0) .or. (ilnTT == 0)) then
              errmsg=" requires lnrho and lnTT"
              goto 30
            endif
            if (.not.(j==ilnTT .or. j == ilnrho)) then
              errmsg=' not allowed for variable no. '//trim(cjvar)
              goto 30
            endif
          case ('cp')
            if (j/=ilnrho) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 10; endif
          case ('sT')
            if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 10; endif
          case ('ctz')
            if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 10; endif
          case ('cdz')
            ! BCZ_DOC: for interstellar runs limit rho
            !call bc_cdz(f,topbot,j)
          case ('ism')
            ! BCZ_DOC: exponential decay/growth in rho/T by scale height
            !call bc_ism(f,topbot,j)
          case ('asT')
            if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 10; endif
          case ('c2')
            if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 10; endif
          case ('ce')
            if (j/=iss) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 10; endif
          case ('f','fa')
            ! BCZ_DOC: freeze value + antisymmetry
            ! tell other modules not to change boundary value
            !call bc_freeze_var_z(topbot,j)
            !call bc_sym_z(f,topbot,j,-1,REL=.true.) ! antisymm wrt boundary
          case ('fs')
            ! BCZ_DOC: freeze value + symmetry
            ! tell other modules not to change boundary value
            !call bc_freeze_var_z(topbot,j)
            !call bc_sym_z(f,topbot,j,+1) ! symmetric wrt boundary
          case ('fBs')
            ! BCZ_DOC: frozen-in B-field (s)
            !call bc_frozen_in_bb(topbot,j)
            !call bc_sym_z(f,topbot,j,+1) ! symmetry
          case ('fB')
            ! BCZ_DOC: frozen-in B-field (a2)
            !call bc_frozen_in_bb(topbot,j)
            !call bc_sym_z(f,topbot,j,-1,REL=.true.) ! antisymm wrt boundary
          case ('StS')
            if (j/=ilnrho) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 30; endif
          case ('div')
            ! BCZ_DOC: set the divergence of $\uv$ to a given value
            ! BCZ_DOC: use bc = 'div' for iuz
            if (j/=iuz) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 30; endif
          case ('win')
            if (j/=ilnrho) then; errmsg=' not allowed for variable no. '//trim(cjvar); goto 30; endif
          case default
            bc%bcname=bcz12(j,topbot)
            bc%ivar=j
            bc%location=(((topbot-1)*6)-3)   ! -3/3 for z bot/top
            bc%value1=fbcz_1(j,topbot)
            bc%value2=fbcz_2(j,topbot)
            bc%done=.false.
!
            if (lspecial) call special_boundconds(f,bc)
!
            if (.not.bc%done) then; errmsg=" no such boundary condition for j="//trim(cjvar); goto 10; endif
          endselect
  30      call fatal_error('check_bounconds_x',trim(bc_code)//errmsg)
        enddo
      enddo
!
    endsubroutine check_boundconds_z
