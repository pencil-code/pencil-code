! $Id: viscosity.f90 13444 2010-03-12 11:34:41Z dhruba.mitra $
!
!  This modules implements viscous heating and diffusion terms
!  here for cases 1) nu constant, 2) mu = rho.nu 3) constant and
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lviscosity = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED fvisc(3); diffus_total; diffus_total2; diffus_total3
! PENCILS PROVIDED visc_heat; nu; gradnu(3); sgnu(3)
!
!***************************************************************
module Viscosity
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'viscosity.h'
!
  integer, parameter :: nvisc_max=4
  character (len=labellen), dimension(nvisc_max) :: ivisc=''
  real :: nu=0.0, nu_shock=0.0
!
  logical :: lvisc_first=.false.
  logical :: lvisc_simplified=.false.
  logical :: lvisc_rho_nu_const=.false.
  logical :: lvisc_nu_const=.false.
  logical :: lvisc_nu_shock=.false.
  logical :: lvisc_heat_as_aux=.false.
  logical, pointer:: lviscosity_heat
!
  namelist /viscosity_run_pars/ &
      nu, ivisc, nu_shock, lvisc_heat_as_aux
!
! other variables (needs to be consistent with reset list below)
  integer :: idiag_fviscm=0     ! DIAG_DOC: Mean value of viscous acceleration
  integer :: idiag_fviscmin=0   ! DIAG_DOC: Min value of viscous acceleration
  integer :: idiag_fviscmax=0   ! DIAG_DOC: Max value of viscous acceleration
  integer :: idiag_epsK=0  ! DIAG_DOC: $\left<2\nu\varrho\Strain^2\right>$
  integer :: idiag_epsK2=0      ! DIAG_DOC:
  integer :: idiag_epsK_LES=0   ! DIAG_DOC:
  integer :: idiag_dtnu=0       ! DIAG_DOC: $\delta t/[c_{\delta t,{\rm v}}\,
                                ! DIAG_DOC:   \delta x^2/\nu_{\rm max}]$
                                ! DIAG_DOC:   \quad(time step relative to
                                ! DIAG_DOC:   viscous time step;
                                ! DIAG_DOC:  see \S~\ref{time-step})
  integer :: idiag_meshRemax=0  ! DIAG_DOC:
  integer :: idiag_nuD2uxbxm=0  ! DIAG_DOC:
  integer :: idiag_nuD2uxbym=0  ! DIAG_DOC:
  integer :: idiag_nuD2uxbzm=0  ! DIAG_DOC:

  contains
!***********************************************************************
    subroutine register_viscosity()
!
!  19-nov-02/tony: coded
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id: viscosity.f90 13444 2010-03-12 11:34:41Z dhruba.mitra $")
!
    endsubroutine register_viscosity
!***********************************************************************
    subroutine initialize_viscosity(lstarting)
!
!  20-nov-02/tony: coded
!
      use FArrayManager, only: farray_register_auxiliary
      use Mpicomm, only: stop_it
      use SharedVariables, only: put_shared_variable,get_shared_variable
!
      logical, intent(in) :: lstarting
!
      integer :: i, ierr
!
!  Default viscosity.
!
      if ( (nu/=0.0).and.(ivisc(1)=='') ) ivisc(1)='nu-const'
!
!  Some viscosity types need the rate-of-strain tensor and grad(lnrho)
!
      lvisc_simplified=.false.
      lvisc_rho_nu_const=.false.
      lvisc_nu_const=.false.
      lvisc_nu_shock=.false.
!
      do i=1,nvisc_max
        select case (ivisc(i))
        case ('simplified', '0')
          if (lroot) print*,'viscous force: nu*del2v'
          lvisc_simplified=.true.
        case ('rho-nu-const','rho_nu-const', '1')
          if (lroot) print*,'viscous force: mu/rho*(del2u+graddivu/3)'
          lvisc_rho_nu_const=.true.
        case ('nu-const')
          if (lroot) print*,'viscous force: nu*(del2u+graddivu/3+2S.glnrho)'
          if (nu/=0.) lpenc_requested(i_sij)=.true.
          lvisc_nu_const=.true.
        case ('nu-shock','shock')
          if (lroot) print*,'viscous force: nu_shock*(XXXXXXXXXXX)'
          lvisc_nu_shock=.true.
          if (.not. lshock) &
           call stop_it('initialize_viscosity: shock viscosity'// &
                           ' but module setting SHOCK=noshock')
        case ('none')
          ! do nothing
        case ('')
          ! do nothing
        case default
          if (lroot) print*, 'No such value for ivisc(',i,'): ', trim(ivisc(i))
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
        if (lvisc_nu_shock.and.nu_shock==0.0) &
            call fatal_error('initialize_viscosity', &
            'Viscosity coefficient nu_shock is zero!')
      endif
!
!  Register an extra aux slot for dissipation rate if requested (so
!  visc_heat is written sto snapshots and can be easily analyzed later).
!    NB: We are doing this here, rather than in register_viscosity, as the
!  register_XXX routines are called before read_{start,run}pars, so
!  lvisc_heat_as_aux isn't known there. This implies that we need to
!  append the ivisc_heat line to index.pro manually.
!
      if (lvisc_heat_as_aux) then
        call farray_register_auxiliary('visc_heat',ivisc_heat)
!
        if (lroot) then
          open(3,file=trim(datadir)//'/index.pro', POSITION='append')
          write(3,*) 'ivisc_heat=',ivisc_heat
          close(3)
        endif
      endif
!
!  Shared variables.
!
      call get_shared_variable('lviscosity_heat',lviscosity_heat,ierr)
      if (ierr/=0) call stop_it("initialize_viscosity: " &
          // "problem getting shared var lviscosity_heat")

!
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_viscosity
!***********************************************************************
    subroutine read_viscosity_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      call keep_compiler_quiet(present(iostat))
!
    endsubroutine read_viscosity_init_pars
!***********************************************************************
    subroutine write_viscosity_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_viscosity_init_pars
!***********************************************************************
    subroutine read_viscosity_run_pars(unit,iostat)
!
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
!
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
      use Diagnostics, only: parse_name
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
      if (lroot.and.ip<14) print*,'rprint_viscosity: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'fviscm',idiag_fviscm)
        call parse_name(iname,cname(iname),cform(iname),'fviscmin',idiag_fviscmin)
        call parse_name(iname,cname(iname),cform(iname),'fviscmax',idiag_fviscmax)
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
!  write column where which viscosity variable is stored
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
          write(3,*) 'itest=',0
        endif
      endif
!
    endsubroutine rprint_viscosity
!***********************************************************************
    subroutine pencil_criteria_viscosity()
!
!  All pencils that the Viscosity module depends on are specified here.
!
!  20-11-04/anders: coded
!
      if (lentropy .and. &
          (lvisc_simplified .or. lvisc_rho_nu_const .or. &
           lvisc_nu_const .or. lvisc_nu_shock)) &
            lpenc_requested(i_TT1)=.true.
      if (lvisc_rho_nu_const .or. lvisc_nu_const) then 
        if (lentropy.and.lviscosity_heat) &
             lpenc_requested(i_sij2)=.true.
        lpenc_requested(i_graddivu)=.true.
      endif
      if (lvisc_simplified .or. lvisc_rho_nu_const .or. lvisc_nu_const) &
          lpenc_requested(i_del2u)=.true.
      if (lvisc_rho_nu_const) lpenc_requested(i_rho1)=.true.
      if (lvisc_nu_const) lpenc_requested(i_sglnrho)=.true.
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
      if (lvisc_heat_as_aux) then
        lpenc_diagnos(i_visc_heat)=.true.
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_sij2)=.true.
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
      logical, dimension (npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
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
      use Diagnostics, only: max_mn_name, sum_mn_name
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (nx,3) :: tmp,tmp2
      real, dimension (nx) :: murho1
!
      integer :: i
!
      intent(inout) :: f,p
!
!  Viscous force and viscous heating are calculated here (for later use).
!
      p%fvisc=0.0
      if (lpencil(i_visc_heat)) p%visc_heat=0.0
      if (lfirst.and.ldt) p%diffus_total=0.0
!
      if (lvisc_simplified) then
!
!  viscous force: nu*del2v
!  -- not physically correct (no momentum conservation), but
!  numerically easy and in most cases qualitatively OK
!
        p%fvisc=p%fvisc+nu*p%del2u
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          if (headtt) then
            call warning('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_simplified')
          endif
        endif
!
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
          p%fvisc(:,i)=p%fvisc(:,i) + &
              murho1*(p%del2u(:,i)+1.0/3.0*p%graddivu(:,i))
        enddo
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*nu*p%sij2*p%rho1
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+murho1
      endif
!
      if (lvisc_nu_const) then
!
!  viscous force: nu*(del2u+graddivu/3+2S.glnrho)
!  -- the correct expression for nu=const
!
        if (ldensity) &
          p%fvisc = p%fvisc + 2*nu*p%sglnrho+nu*(p%del2u + 1./3.*p%graddivu)
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*nu*p%sij2
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+nu
      endif
!
!  viscous force: nu_shock
!
      if (lvisc_nu_shock) then
        if (ldensity) then
          call multsv(p%divu,p%glnrho,tmp2)
          tmp=tmp2 + p%graddivu
          call multsv(nu_shock*p%shock,tmp,tmp2)
          call multsv_add(tmp2,nu_shock*p%divu,p%gshock,tmp)
          p%fvisc=p%fvisc+tmp
          if (lfirst.and.ldt) p%diffus_total=p%diffus_total+(nu_shock*p%shock)
          if (lpencil(i_visc_heat)) &
              p%visc_heat=p%visc_heat+nu_shock*p%shock*p%divu**2
        endif
      endif
!
     if (lvisc_heat_as_aux) f(l1:l2,m,n,ivisc_heat) = p%visc_heat
!
    endsubroutine calc_pencils_viscosity
!***********************************************************************
    subroutine calc_viscosity(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_viscosity
!***********************************************************************
    subroutine calc_viscous_heat(df,p,Hmax)
!
!  calculate viscous heating term for right hand side of entropy equation
!
!  20-nov-02/tony: coded
!
      real, dimension (mx,my,mz,mvar)    :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: Hmax
!
!  Add viscous heat (which has units of energy/mass) to the RHS
!  of the entropy (both with and without pretend_lnTT), or of
!  the temperature equation. Divide by cv if pretend_lnTT.
!
      if (lentropy) &
        df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + p%TT1*p%visc_heat
!
!  calculate maximum heating (for time step constraint), so it is
!  only done on the first of the 3 substeps.
!
      if (lfirst .and. ldt) Hmax=Hmax+p%visc_heat
!
    endsubroutine calc_viscous_heat
!***********************************************************************
    subroutine calc_viscous_force(df,p)
!
!  Calculate viscous force term for right hand side of equation of motion.
!
!  20-nov-02/tony: coded
!   9-jul-04/nils: added Smagorinsky viscosity
!
      use Diagnostics, only: sum_mn_name, max_mn_name
      use Sub, only: cross
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
      if (lfirst.and.ldt) diffus_nu =p%diffus_total *dxyz_2
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_fviscm/=0)   call sum_mn_name(p%fvisc,idiag_fviscm)
        if (idiag_fviscmin/=0) call max_mn_name(-p%fvisc,idiag_fviscmin,lneg=.true.)
        if (idiag_fviscmax/=0) call max_mn_name(p%fvisc,idiag_fviscmax)
        if (idiag_dtnu/=0) &
            call max_mn_name(diffus_nu/cdtv,idiag_dtnu,l_dt=.true.)
        if (idiag_nu_LES /= 0) call sum_mn_name(nu_smag,idiag_nu_LES)
        if (idiag_meshRemax/=0) &
           call max_mn_name(sqrt(p%u2(:))*dxmax/p%diffus_total,idiag_meshRemax)
!  Viscous heating as explicit analytical term.
        if (idiag_epsK/=0) then
          if (lvisc_nu_const)     call sum_mn_name(2*nu*p%rho*p%sij2,idiag_epsK)
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
        if (idiag_epsK2/=0) call sum_mn_name(p%visc_heat*p%rho,idiag_epsK2)
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
    endsubroutine calc_viscous_force
!***********************************************************************
endmodule Viscosity
