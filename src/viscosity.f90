
! $Id: viscosity.f90,v 1.34 2006-11-03 13:50:41 brandenb Exp $

!  This modules implements viscous heating and diffusion terms
!  here for cases 1) nu constant, 2) mu = rho.nu 3) constant and 

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lviscosity = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
! PENCILS PROVIDED fvisc, diffus_total, visc_heat
!
!***************************************************************

module Viscosity

  use Cparam
  use Cdata
  use Messages

  implicit none

  include 'viscosity.h'

  integer, parameter :: nvisc_max = 4
  character (len=labellen), dimension(nvisc_max) :: ivisc=''
  real :: nu=0., nu_mol=0., nu_hyper2=0., nu_hyper3=0., nu_shock=0.
  real, dimension(3) :: nuvec_hyper3=0.

  ! dummy logical
  logical :: lvisc_first=.false.

  logical :: lvisc_simplified=.false.
  logical :: lvisc_rho_nu_const=.false.
  logical :: lvisc_nu_const=.false.
  logical :: lvisc_nu_shock=.false.
  logical :: lvisc_hyper2_simplified=.false.
  logical :: lvisc_hyper3_simplified=.false.
  logical :: lvisc_hyper3_rho_nu_const=.false.
  logical :: lvisc_hyper3_nu_vector=.false.
  logical :: lvisc_hyper3_rho_nu_const_bulk=.false.
  logical :: lvisc_hyper3_nu_const=.false.
  logical :: lvisc_smag_simplified=.false.
  logical :: lvisc_smag_cross_simplified=.false.
  logical :: lvisc_snr_damp=.false.

  ! input parameters
  !integer :: dummy1
  !namelist /viscosity_init_pars/ dummy1

  ! run parameters
  namelist /viscosity_run_pars/ &
      nu, nu_hyper2, nu_hyper3, ivisc, nu_mol, C_smag, nu_shock, nuvec_hyper3
 
  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_epsK=0,idiag_epsK2=0,idiag_epsK_LES=0
  integer :: idiag_dtnu=0
  integer :: idiag_meshRemax=0
  integer :: idiag_nuD2uxbxm=0, idiag_nuD2uxbym=0, idiag_nuD2uxbzm=0

  contains

!***********************************************************************
    subroutine register_viscosity()
!
!  19-nov-02/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_viscosity called twice')
      first = .false.
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_viscosity: constant viscosity'
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: viscosity.f90,v 1.34 2006-11-03 13:50:41 brandenb Exp $")

      ivisc(1)='nu-const'

! Following test unnecessary as no extra variable is evolved
!
!      if (nvar > mvar) then
!        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
!        call stop_it('Register_viscosity: nvar > mvar')
!      endif
!
    endsubroutine register_viscosity
!***********************************************************************
    subroutine initialize_viscosity(lstarting)
!
!  20-nov-02/tony: coded
!
      use Cdata
      use Mpicomm, only: stop_it

      logical, intent(in) :: lstarting
      integer :: i
!
!  Some viscosity types need the rate-of-strain tensor and grad(lnrho)
!
      lvisc_simplified=.false.
      lvisc_rho_nu_const=.false.
      lvisc_nu_const=.false.
      lvisc_nu_shock=.false.
      lvisc_hyper2_simplified=.false.
      lvisc_hyper3_simplified=.false.
      lvisc_hyper3_rho_nu_const=.false.
      lvisc_hyper3_nu_vector=.false.
      lvisc_hyper3_rho_nu_const_bulk=.false.
      lvisc_hyper3_nu_const=.false.
      lvisc_smag_simplified=.false.
      lvisc_smag_cross_simplified=.false.
      lvisc_snr_damp=.false.

      do i=1,nvisc_max
        select case (ivisc(i))
        case ('simplified', '0')
          if (lroot) print*,'viscous force: nu*del2v'
          lvisc_simplified=.true.
        case('rho_nu-const', '1')
          if (lroot) print*,'viscous force: mu/rho*(del2u+graddivu/3)'
          lvisc_rho_nu_const=.true.
        case('nu-const')
          if (lroot) print*,'viscous force: nu*(del2u+graddivu/3+2S.glnrho)'
          if (nu/=0.) lpenc_requested(i_sij)=.true.
          lvisc_nu_const=.true.
        case('nu-shock')
          if (lroot) print*,'viscous force: nu_shock*(XXXXXXXXXXX)'
          lvisc_nu_shock=.true.
          if (.not. lshock) &
           call stop_it('initialize_viscosity: shock viscosity'// &
                           ' but module setting SHOCK=noshock')
        case ('hyper2_simplified', 'hyper4')
          if (lroot) print*,'viscous force: nu_hyper*del4v'
          lvisc_hyper2_simplified=.true.
        case ('hyper3_simplified', 'hyper6')
          if (lroot) print*,'viscous force: nu_hyper*del6v'
          lvisc_hyper3_simplified=.true.
        case ('hyper3_rho_nu-const')
          if (lroot) print*,'viscous force: nu_hyper/rho*del6v'
          lvisc_hyper3_rho_nu_const=.true.
       case ('hyper3_nu-vector')
          if (lroot) print*,&
               'viscous force(i): 1/rho*(nu_x*d^6/dx^6 + nu_y*d^6/dy^6 + nu_z*d^6/dz^6)v'
          lvisc_hyper3_nu_vector=.true.
        case ('hyper3_rho_nu-const_bulk')
          if (lroot) print*,'viscous force: duj/dt = nu_hyper/rho*d6uj/dxj6'
          lvisc_hyper3_rho_nu_const_bulk=.true.
        case ('hyper3_nu-const')
          if (lroot) print*,'viscous force: nu*(del6u+S.glnrho)'
          lpenc_requested(i_uij5)=.true.
          lvisc_hyper3_nu_const=.true.
        case ('smagorinsky_simplified')
          if (lroot) print*,'viscous force: Smagorinsky_simplified'
          if (lroot) lvisc_LES=.true.
          lpenc_requested(i_sij)=.true.
          lvisc_smag_simplified=.true.
        case ('smagorinsky_cross_simplified')
          if (lroot) print*,'viscous force: Smagorinsky_simplified'
          if (lroot) lvisc_LES=.true.
          lvisc_smag_cross_simplified=.true.
        case ('snr_damp')
          if (lroot) print*,'viscous force: SNR damping'
          lvisc_snr_damp=.true.
        case ('')
          ! do nothing
        case default
          if (lroot) print*, 'No such such value for ivisc(',i,'): ', trim(ivisc(i))
          call stop_it('calc_viscous_forcing')
        endselect
      enddo
!
!  If we're timestepping, die or warn if the viscosity coefficient that
!  corresponds to the chosen viscosity type is not set.
!
      if (lrun) then
        if ( (lvisc_simplified.or.lvisc_rho_nu_const.or.lvisc_nu_const) &
            .and.nu==0.0) &
            call warning('initialize_viscosity', &
            'Viscosity coefficient nu is zero!')
        if (lvisc_hyper2_simplified.and.nu_hyper2==0.0) &
            call fatal_error('initialize_viscosity', &
            'Viscosity coefficient nu_hyper2 is zero!')
        if ( (lvisc_hyper3_simplified.or.lvisc_hyper3_rho_nu_const.or. &
              lvisc_hyper3_rho_nu_const_bulk.or.lvisc_hyper3_nu_const).and. &
              nu_hyper3==0.0 ) &
            call fatal_error('initialize_viscosity', &
            'Viscosity coefficient nu_hyper3 is zero!')
        if ( (lvisc_hyper3_nu_vector) .and.  &         
             ((nuvec_hyper3(1)==0. .and. nxgrid/=1 ).or. &
              (nuvec_hyper3(2)==0. .and. nygrid/=1 ).or. &
              (nuvec_hyper3(3)==0. .and. nzgrid/=1 )) ) & 
            call fatal_error('initialize_viscosity', &
             'A viscosity coefficient of nuvec_hyper3 is zero!')
        if ( (lvisc_smag_simplified.or.lvisc_smag_cross_simplified).and. &
             C_smag==0.0 ) &
            call fatal_error('initialize_viscosity', &
            'Viscosity coefficient C_smag is zero!')
        if (lvisc_nu_shock.and.nu_shock==0.0) &
            call fatal_error('initialize_viscosity', &
            'Viscosity coefficient nu_shock is zero!')
      endif
!
      if (NO_WARN) print*, lstarting
!
    endsubroutine initialize_viscosity
!***********************************************************************
    subroutine read_viscosity_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit
!
    endsubroutine read_viscosity_init_pars
!***********************************************************************
    subroutine write_viscosity_init_pars(unit)
      integer, intent(in) :: unit
!
      if (NO_WARN) print*,unit
!
    endsubroutine write_viscosity_init_pars
!***********************************************************************
    subroutine read_viscosity_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=viscosity_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=viscosity_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_viscosity_run_pars
!***********************************************************************
    subroutine write_viscosity_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=viscosity_run_pars)
!
    endsubroutine write_viscosity_run_pars
!*******************************************************************
    subroutine rprint_viscosity(lreset,lwrite)
!
!  Writes ishock to index.pro file
!
!  24-nov-03/tony: adapted from rprint_ionization
!
      use Cdata
      use Sub
! 
      logical :: lreset
      logical, optional :: lwrite
      integer :: iname
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtnu=0
        idiag_nu_LES=0
        idiag_epsK=0
        idiag_epsK2=0
        idiag_epsK_LES=0
        idiag_meshRemax=0
        idiag_nuD2uxbxm=0; idiag_nuD2uxbym=0; idiag_nuD2uxbzm=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'rprint_viscosity: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtnu',idiag_dtnu)
        call parse_name(iname,cname(iname),cform(iname),'nu_LES',idiag_nu_LES)
        call parse_name(iname,cname(iname),cform(iname),'epsK',idiag_epsK)
        call parse_name(iname,cname(iname),cform(iname),'epsK2',idiag_epsK2)
        call parse_name(iname,cname(iname),cform(iname),&
            'epsK_LES',idiag_epsK_LES)
        call parse_name(iname,cname(iname),cform(iname),&
            'meshRemax',idiag_meshRemax)
      enddo
!
!  write column where which ionization variable is stored
!
      if (present(lwrite)) then
        if (lwrite) then
          write(3,*) 'i_dtnu=',idiag_dtnu
          write(3,*) 'i_nu_LES=',idiag_nu_LES
          write(3,*) 'i_epsK=',idiag_epsK
          write(3,*) 'i_epsK2=',idiag_epsK2
          write(3,*) 'i_epsK_LES=',idiag_epsK_LES
          write(3,*) 'i_meshRemax=',idiag_meshRemax
          write(3,*) 'i_nuD2uxbxm=',idiag_nuD2uxbxm
          write(3,*) 'i_nuD2uxbym=',idiag_nuD2uxbym
          write(3,*) 'i_nuD2uxbzm=',idiag_nuD2uxbzm
          write(3,*) 'ihyper=',ihyper
          write(3,*) 'itest=',0
        endif
      endif
!   
      if(NO_WARN) print*,lreset  !(to keep compiler quiet)
!        
    endsubroutine rprint_viscosity
!***********************************************************************
    subroutine pencil_criteria_viscosity()
!    
!  All pencils that the Viscosity module depends on are specified here.
!
!  20-11-04/anders: coded
!
      if ((lentropy.or.ltemperature) .and. &
          (lvisc_simplified .or. lvisc_rho_nu_const .or. &
           lvisc_nu_const .or. lvisc_nu_shock)) lpenc_requested(i_TT1)=.true.
      if (lvisc_rho_nu_const .or. lvisc_nu_const .or. &
          lvisc_smag_cross_simplified) then
        if (lentropy.or.ltemperature) lpenc_requested(i_sij2)=.true.
        lpenc_requested(i_graddivu)=.true.
      endif
      if (lvisc_smag_simplified) lpenc_requested(i_ss12)=.true.
      if (lvisc_simplified .or. lvisc_rho_nu_const .or. lvisc_nu_const .or. &
          lvisc_smag_simplified .or. lvisc_smag_cross_simplified) &
          lpenc_requested(i_del2u)=.true.
      if (lvisc_hyper3_simplified .or. lvisc_hyper3_rho_nu_const .or. &
          lvisc_hyper3_nu_const) &
          lpenc_requested(i_del6u)=.true.
      if (lvisc_hyper3_rho_nu_const_bulk) lpenc_requested(i_del6u_bulk)=.true.
      if (lvisc_hyper2_simplified) lpenc_requested(i_del4u)=.true.
      if (lvisc_rho_nu_const .or. lvisc_hyper3_rho_nu_const .or. &
          lvisc_hyper3_rho_nu_const_bulk .or.  lvisc_hyper3_nu_vector .or. &
          lvisc_smag_simplified .or. lvisc_smag_cross_simplified) &
          lpenc_requested(i_rho1)=.true.
      if (lvisc_nu_const .or. &
          lvisc_smag_simplified .or. lvisc_smag_cross_simplified)  &
          lpenc_requested(i_sglnrho)=.true.
      if (lvisc_hyper3_nu_const) lpenc_requested(i_uij5glnrho)=.true.
      if (ldensity.and.lvisc_nu_shock) then
        lpenc_requested(i_graddivu)=.true.
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_gshock)=.true.
        lpenc_requested(i_divu)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
!
      if (idiag_meshRemax/=0) lpenc_diagnos(i_u2)=.true.
      if (idiag_epsK/=0.or.idiag_epsK_LES/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_sij2)=.true.
      endif
      if (idiag_epsK2/=0) then
        lpenc_diagnos(i_visc_heat)=.true.
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_uu)=.true.
      endif
      if (lvisc_nu_shock.and.idiag_epsK/=0) then
        lpenc_diagnos(i_fvisc)=.true.
        lpenc_diagnos(i_diffus_total)=.true.
        lpenc_diagnos(i_shock)=.true.
        lpenc_diagnos(i_divu)=.true.
        lpenc_diagnos(i_rho)=.true.
      endif
      if ( (idiag_meshRemax/=0 .or. idiag_dtnu/=0) .and. lvisc_nu_shock) &
          lpenc_diagnos(i_shock)=.true.
!
    endsubroutine pencil_criteria_viscosity
!***********************************************************************
    subroutine pencil_interdep_viscosity(lpencil_in)
!
!  Interdependency among pencils from the Viscosity module is specified here.
!
!  20-11-04/anders: coded
!
      use Cdata
!
      logical, dimension (npencils) :: lpencil_in
!      
      if (NO_WARN) print*, lpencil_in !(keep compiler quiet)

      if (lpencil_in(i_visc_heat)) lpencil_in(i_rho)=.true.
!
    endsubroutine pencil_interdep_viscosity
!***********************************************************************
    subroutine calc_pencils_viscosity(f,p)
!     
!  Calculate Viscosity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      use Cdata
      use Sub
      use Interstellar, only: calc_snr_damping
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (nx,3) :: tmp,tmp2,tmp3
      real, dimension (nx) :: murho1,nu_smag
!
      integer :: i
!
      intent(in) :: f
      intent(inout) :: p
!
!  Viscosity operator
!
      p%fvisc=0.0
      p%diffus_total=0.0
      p%visc_heat=0.0
!!
!! Calculate the heating term due to the viscous force
!!
!      if (lpencil(i_visc_heat)) then
!        call dot_mn(p%uu,p%fvisc,ufvisc)
!        p%visc_heat = -ufvisc
!      endif
!
      if (lvisc_simplified) then
!
!  viscous force: nu*del2v
!  -- not physically correct (no momentum conservation), but
!  numerically easy and in most cases qualitatively OK
!
        p%fvisc=p%fvisc+nu*p%del2u
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          call fatal_error('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_simplified')
        endif
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+nu
      endif
!
      if (lvisc_rho_nu_const) then
!
!  viscous force: mu/rho*(del2u+graddivu/3)
!  -- the correct expression for rho*nu=const
!
        murho1=nu*p%rho1  !(=mu/rho)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i)+murho1*(p%del2u(:,i)+1./3.*p%graddivu(:,i))
        enddo  
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat + 2*nu*p%sij2*p%rho1
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+murho1
      endif
!
      if (lvisc_nu_const) then
!
!  viscous force: nu*(del2u+graddivu/3+2S.glnrho)
!  -- the correct expression for nu=const
!
        if(ldensity) then
          p%fvisc=p%fvisc+2*nu*p%sglnrho+nu*(p%del2u+1./3.*p%graddivu)
        else
          p%fvisc=p%fvisc+nu*(p%del2u+1./3.*p%graddivu)
        endif
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat + 2*nu*p%sij2
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+nu
     endif
!
!  viscous force: nu_shock
!
      if (lvisc_nu_shock) then
        if (ldensity) then
          tmp=nu_shock*(p%shock*(p%divu*p%glnrho+p%graddivu)+p%divu*p%gshock)
          !call multsv(p%divu,p%glnrho,tmp2)
          !tmp=tmp2 + p%graddivu
          !call multsv(nu_shock*p%shock,tmp,tmp2)
          !call multsv_add(tmp2,nu_shock*p%divu,p%gshock,tmp)
          if (lfirst.and.ldt) p%diffus_total=p%diffus_total+(nu_shock*p%shock)
          if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat + &
                                                 nu_shock*p%shock*p%divu**2
          p%fvisc=p%fvisc+tmp
        endif
      endif
!
      if (lvisc_hyper2_simplified) then
!
!  viscous force: nu_hyper2*de46v (not momentum-conserving)
!
        p%fvisc=p%fvisc+nu_hyper2*p%del4u
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          call fatal_error('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper2_simplified')
        endif
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+nu_hyper2*dxyz_4/dxyz_2
      endif
!
      if (lvisc_hyper3_simplified) then
!
!  viscous force: nu_hyper3*del6v (not momentum-conserving)
!
        p%fvisc=p%fvisc+nu_hyper3*p%del6u
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          call fatal_error('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper3_simplified')
        endif
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+nu_hyper3*dxyz_6/dxyz_2
      endif
!
      if (lvisc_hyper3_rho_nu_const) then
!
!  viscous force: mu/rho*del6u
!
        murho1=nu_hyper3*p%rho1  ! (=mu_hyper3/rho)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i)+murho1*p%del6u(:,i)
        enddo
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          call fatal_error('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper3_rho_nu_const')
        endif
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+nu_hyper3*dxyz_6/dxyz_2
      endif
!
      if (lvisc_hyper3_nu_vector) then
!
!  viscous force: f_i = mu_i/rho*del6u
!  Used for non-cubic cells
!         
         call del6fjv(f,nuvec_hyper3,iuu,tmp3)
         do i=1,3
            p%fvisc(:,i)=p%fvisc(:,i)+tmp3(:,i)*p%rho1
         enddo
         if (lpencil(i_visc_heat)) then  ! Heating term not implemented
            call fatal_error('calc_pencils_viscosity', 'viscous heating term '// &
                 'is not implemented for lvisc_hyper3_nu_vector')
         endif
!
! WL: why is diffus due to viscosity prop to dx_1**4 while
! the one due to mass diffusion and resistivity is prop to dx_1**6 ?
!
         if (lfirst.and.ldt) p%diffus_total=p%diffus_total+&
                 nuvec_hyper3(1)*dx_1(l1:l2)**4 + &
                 nuvec_hyper3(2)*dy_1(m)**4 + &
                 nuvec_hyper3(3)*dz_1(n)**4
!
      endif
!
      if (lvisc_hyper3_rho_nu_const_bulk) then
!
!  viscous force: mu/rho*d6uj/dx6
!
        murho1=nu_hyper3*p%rho1  ! (=mu_hyper3/rho)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i)+murho1*p%del6u_bulk(:,i)
        enddo
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          call fatal_error('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper3_rho_nu_const_bulk')
        endif
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+nu_hyper3*dxyz_6/dxyz_2
      endif
!
      if (lvisc_hyper3_nu_const) then
!
!  viscous force: nu_hyper3*(del6u+S.glnrho), where S_ij=d^5 u_i/dx_j^5
!
        p%fvisc=p%fvisc+nu_hyper3*(p%del6u+p%uij5glnrho)
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          call fatal_error('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper3_nu_const')
        endif
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+nu_hyper3*dxyz_6/dxyz_2
      endif
!
!  viscous force: Handle damping at the core of SNRs
!
      if (linterstellar.and.lvisc_snr_damp) then
        call calc_snr_damping(p)
      endif
!
!  viscous force: nu_hyper3*(del6u+S.glnrho), where S_ij=d^5 u_i/dx_j^5
!
      if (lvisc_smag_simplified) then
!
!  viscous force: nu_smag*(del2u+graddivu/3+2S.glnrho)
!  where nu_smag=(C_smag*dxmax)**2*sqrt(2*SS)
!
        if (ldensity) then
!
! Find nu_smag
!
          nu_smag=(C_smag*dxmax)**2.*sqrt(2*p%sij2)
!
! Calculate viscous force
!
          call multsv_mn(nu_smag,p%sglnrho,tmp2)
          call multsv_mn(nu_smag,p%del2u+1./3.*p%graddivu,tmp)
          p%fvisc=p%fvisc+2*tmp2+tmp
          if (lpencil(i_visc_heat)) then  ! Heating term not implemented
            call fatal_error('calc_pencils_viscosity','viscous heating term '//&
                'is not implemented for lvisc_smag_simplified')
          endif
          if (lfirst.and.ldt) p%diffus_total=p%diffus_total+nu_smag
        else
          if (lfirstpoint) print*, 'calc_viscous_force: '// &
              "ldensity better be .true. for ivisc='smagorinsky'"
        endif
      endif
!
      if (lvisc_smag_cross_simplified) then
!
!  viscous force: nu_smag*(del2u+graddivu/3+2S.glnrho)
!  where nu_smag=(C_smag*dxmax)**2*sqrt(S:J)
!
        if (ldensity) then
          nu_smag=(C_smag*dxmax)**2.*p%ss12
!
! Calculate viscous force
!
          call multsv_mn(nu_smag,p%sglnrho,tmp2)
          call multsv_mn(nu_smag,p%del2u+1./3.*p%graddivu,tmp)
          p%fvisc=p%fvisc+2*tmp2+tmp
          if (lpencil(i_visc_heat)) then  ! Heating term not implemented
            call fatal_error('calc_pencils_viscosity','viscous heating term '//&
                'is not implemented for lvisc_smag_cross_simplified')
          endif
          if (lfirst.and.ldt) p%diffus_total=p%diffus_total+nu_smag
        else
          if (lfirstpoint) print*, 'calc_viscous_force: '// &
              "ldensity better be .true. for ivisc='smagorinsky'"
        endif
     endif
!
     if (NO_WARN) print*, f    !(to keep compiler quiet)
!
    endsubroutine calc_pencils_viscosity
!***********************************************************************
    subroutine calc_viscosity(f)
!    
      real, dimension (mx,my,mz,mfarray) :: f
!
      if(NO_WARN) print*,f  !(to keep compiler quiet)
!
    endsubroutine calc_viscosity
!***********************************************************************
    subroutine calc_viscous_heat(df,p,Hmax)
!
!  calculate viscous heating term for right hand side of entropy equation
!
!  20-nov-02/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!      
      real, dimension (nx) :: Hmax
!
!  Add viscous heat (which has units of energy/mass) to the RHS
!  of the entropy ...
!
      if (lentropy) then
         df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) +  p%TT1*p%visc_heat
      endif
!
!  ... or temperature equation.
!
      if (ltemperature) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + p%cv1*p%TT1*p%visc_heat
      endif
!      
      if (lfirst .and. ldt) Hmax=Hmax+p%visc_heat
!
      if(NO_WARN) print*,p  !(keep compiler quiet)
!        
    endsubroutine calc_viscous_heat
!***********************************************************************
    subroutine calc_viscous_force(df,p)
!
!  calculate viscous force term for right hand side of  equation
!
!  20-nov-02/tony: coded
!   9-jul-04/nils: added Smagorinsky viscosity

      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: nu_smag
      real, dimension (nx,3) :: nuD2uxb
      type (pencil_case) :: p
!      
      intent (in) :: p
      intent (inout) :: df

!
!  Add viscosity to equation of motion
!
      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + p%fvisc
!
!  Calculate max total diffusion coefficient for timestep calculation etc.
!
      diffus_nu=max(diffus_nu,p%diffus_total*dxyz_2)
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (lvisc_smag_simplified) then
          if (ldensity) then
            nu_smag=(C_smag*dxmax)**2.*sqrt(2*p%sij2)
          endif
        endif
        if (lvisc_smag_cross_simplified) then
          if (ldensity) then
            nu_smag=(C_smag*dxmax)**2.*p%ss12
          endif
        endif
        if (idiag_dtnu/=0) &
            call max_mn_name(diffus_nu/cdtv,idiag_dtnu,l_dt=.true.)
        if (idiag_nu_LES /= 0) call sum_mn_name(nu_smag,idiag_nu_LES)
        if (idiag_meshRemax/=0) &
           call max_mn_name(sqrt(p%u2(:))*dxmax/p%diffus_total,idiag_meshRemax)
!  Viscous heating as explicit analytical term.
        if (idiag_epsK/=0) then
!           call dot2(p%uu,fvisc,visc_heat)
!           call sum_mn_name(visc_heat,idiag_epsK)
          if (lvisc_nu_const) call sum_mn_name(2*nu*p%rho*p%sij2,idiag_epsK)
          if (lvisc_rho_nu_const) call sum_mn_name(2*nu*p%sij2,idiag_epsK)
          if (lvisc_nu_shock) &  ! Heating from shock viscosity.
              call sum_mn_name((nu_shock*p%shock*p%divu**2)*p%rho,idiag_epsK)
        endif
!  Viscosity power (per volume):
!    P = u_i tau_ij,j = d/dx_j(u_i tau_ij) - u_i,j tau_ij
!  The first term is a kinetic energy flux and the second an energy loss which
!  is the viscous heating. This can be rewritten as the analytical heating
!  term in the manual (for nu-const viscosity).
!  The average of P should, for periodic boundary conditions, equal the
!  analytical heating term, so that epsK and epsK2 are equals.
!  However, in strongly random flow, there can be significant difference
!  between epsK and epsK2, even for periodic boundaries.
        if (idiag_epsK2/=0) then
          call sum_mn_name(p%visc_heat*p%rho,idiag_epsK2)
        endif
!  Viscous heating for Smagorinsky viscosity.
        if (idiag_epsK_LES/=0) then
          if (lvisc_smag_simplified) then
            call sum_mn_name(2*nu_smag*p%rho*p%sij2,idiag_epsK_LES)
          else if (lvisc_smag_cross_simplified) then
            call sum_mn_name(2*nu_smag*p%rho*p%sij2,idiag_epsK_LES)
          endif
        endif
!
!  correlation of viscous term with b-field; relevant for MTA
!  (MTA: minimal tau approximation)
!
        if (lmagnetic) then
          if (idiag_nuD2uxbxm/=0.or.idiag_nuD2uxbym/=0.or.idiag_nuD2uxbzm/=0) then
!           call curl(f,iaa,bb)
            call cross(p%fvisc,p%bb,nuD2uxb)
            call sum_mn_name(nuD2uxb(:,1),idiag_nuD2uxbxm)
            call sum_mn_name(nuD2uxb(:,2),idiag_nuD2uxbym)
            call sum_mn_name(nuD2uxb(:,3),idiag_nuD2uxbzm)
          endif
        endif
      endif
!
    end subroutine calc_viscous_force
!***********************************************************************

endmodule Viscosity
