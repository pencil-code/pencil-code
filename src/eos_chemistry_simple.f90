! $Id$
!
!  Equation of state for an ideal gas without ionization.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: leos = .true., leos_ionization=.false., leos_temperature_ionization=.false.
! CPARAM logical, parameter :: leos_idealgas = .false., leos_chemistry = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED lnTT;  glnTT(3); TT; TT1; gTT(3)
! PENCILS PROVIDED mu1; pp; del2pp; glnmu(3)
! PENCILS PROVIDED rho1gpp(3); glnpp(3); del2lnTT
! PENCILS PROVIDED glnRR(3), RRmix
!
! PENCILS PROVIDED hss(3,3); hlnTT(3,3); del2ss; del6ss; del6lnTT
! PENCILS PROVIDED yH; ee; ss; delta; glnmumol(3); ppvap; csvap2; cs2
! PENCILS PROVIDED cp1tilde; cp; gamma; gamma_m1
! PENCILS PROVIDED rho_anel; gradcp(3)
!
!
!***************************************************************
module EquationOfState
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'eos.h'
  include 'eos_params.h'
!
  integer :: iglobal_cs2, iglobal_glnTT
!
  real :: mu=1.
  real :: cs0=1., rho0=1.
  real :: cs20=1., lnrho0=0.
  real :: gamma=impossible
  real :: Rgas_cgs=0., Rgas, Rgas_unit_sys=1., error_cp=1e-6, scale_Rgas=1.
  real :: cp=impossible
  real :: cs2bot=1., cs2top=1.
  real :: lnTT0=impossible, TT0=impossible
  integer :: ieosvars=-1, ieosvar1=-1, ieosvar2=-1, ieosvar_count=0
  integer :: ll1,ll2,mm1,mm2,nn1,nn2
  logical :: leos_isothermal=.false., leos_isentropic=.false.
  logical :: leos_isochoric=.false., leos_isobaric=.false.
  logical :: leos_localisothermal=.false.
  character (len=20) :: input_file
  logical ::  lcheminp_eos=.false.
  integer :: imass=1!, iTemp1=2,iTemp2=3,iTemp3=4
  real :: Cp_const=impossible
  real :: Pr_number=0.7
  logical :: lpres_grad = .false.
  logical :: linterp_pressure=.false. !This is only used when ogrid to interpolate pressure instead of temperature
!
 real, dimension(:,:), pointer :: species_constants
!
  namelist /eos_init_pars/ mu, cp, cs0, rho0, gamma, error_cp, Cp_const, lpres_grad, linterp_pressure, scale_Rgas
!
  namelist /eos_run_pars/  mu, cs0, rho0, Cp_const, Pr_number
!
  contains
!***********************************************************************
    subroutine register_eos
!
!  14-jun-03/axel: adapted from register_eos
!
      use FArrayManager
      use SharedVariables, only: put_shared_variable
!
      ilnTT = 0
!
      if (lpres_grad) then
!
!  Register gradient of pressure
!
        call farray_register_auxiliary('gpx',igpx)
        call farray_register_auxiliary('gpy',igpy)
!
      endif
!
      if (linterp_pressure) call farray_register_auxiliary('pp',ipp)
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          '$Id$')
!
      call put_shared_variable('gamma',gamma,caller='register_eos')
      call put_shared_variable('linterp_pressure',linterp_pressure)
      call put_shared_variable('scale_Rgas',scale_Rgas)

    endsubroutine register_eos
!***********************************************************************
    subroutine units_eos
!
!  This routine calculates things related to units and must be called
!  before the rest of the units are being calculated.
!
!  22-jun-06/axel: adapted from initialize_eos
!  16-mar-10/Natalia
!
      logical :: chemin=.false.,cheminp=.false.
!
! Initialize variable selection code (needed for RELOADing)
!
      ieosvars=-1
      ieosvar_count=0
!
      if (unit_system == 'cgs') then
         Rgas_unit_sys = k_B_cgs/m_u_cgs
      elseif (unit_system == 'SI') then
         Rgas_unit_sys = k_B_cgs/m_u_cgs*1.e-4
      endif
!
      if (unit_temperature == impossible) then
        call fatal_error('units_eos','unit_temperature is not found')
      else
        Rgas=Rgas_unit_sys*unit_temperature*scale_Rgas/unit_velocity**2
      endif
!
      inquire(file='chem.inp',exist=cheminp)
      inquire(file='chem.in',exist=chemin)
      if (chemin .and. cheminp) &
        call fatal_error('units_eos','chem.inp and chem.in found. Please decide for one')
!
      if (cheminp) input_file='chem.inp'
      if (chemin) input_file='chem.in'
      lcheminp_eos = cheminp .or. chemin
!
    endsubroutine units_eos
!***********************************************************************
    subroutine initialize_eos(f)
!
! Initialize variable selection code (needed for RELOADing)
!
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f

      ieosvars=-1
      ieosvar_count=0
!
!  write constants to disk. In future we may want to deal with this
!  using an include file or another module.
!
      if (pretend_lnTT) then
        call warning('initialize_eos','pretend_lnTT is not used with eos_chemistry_simple')
        pretend_lnTT=.false.
      endif
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro',position="append")
        write (1,'(a,1pd26.16)') 'k_B=',k_B
        write (1,'(a,1pd26.16)') 'm_H=',m_H
        write (1,*) 'lnTTO=',lnTT0
        write (1,*) 'cp=',cp
        close (1)
      endif
!
      if ((nxgrid==1) .and. (nygrid==1) .and. (nzgrid==1)) then
        ll1=1; ll2=mx; mm1=m1; mm2=m2; nn1=n1; nn2=n2
      elseif (nxgrid==1) then
        ll1=l1; ll2=l2
      else
        ll1=1; ll2=mx
      endif
!
      if (nygrid==1) then
        mm1=m1; mm2=m2
      else
        mm1=1; mm2=my
      endif
!
      if (nzgrid==1) then
        nn1=n1; nn2=n2
      else
        nn1=1;  nn2=mz
      endif
!
      call get_shared_variable('species_constants',species_constants,caller='initialize_eos')
!
    endsubroutine initialize_eos
!***********************************************************************
    subroutine select_eos_variable(variable,findex)
!
!  Select eos variable
!
!   02-apr-06/tony: implemented
!
      use FArrayManager
      use General, only: itoa
!
      character (len=*), intent(in) :: variable
      integer, intent(in) :: findex
      integer :: this_var=0
      integer, save :: ieosvar_selected=0
      integer, parameter :: ieosvar_lnrho = 2**0
      integer, parameter :: ieosvar_rho   = 2**1
      integer, parameter :: ieosvar_ss    = 2**2
      integer, parameter :: ieosvar_lnTT  = 2**3
      integer, parameter :: ieosvar_TT    = 2**4
      integer, parameter :: ieosvar_cs2   = 2**5
      integer, parameter :: ieosvar_pp    = 2**6
!
      if (ieosvar_count==0) ieosvar_selected=0
!
      if (ieosvar_count>=2) call fatal_error("select_eos_variable", &
        "2 thermodynamic quantities have already been defined while attempting to add a 3rd") !//variable)
!
      ieosvar_count=ieosvar_count+1
!
!      select case (variable)
      if (variable=='ss') then
        this_var=ieosvar_ss
        if (findex<0) leos_isentropic=.true.
      elseif (variable=='cs2') then
        this_var=ieosvar_cs2
        if (findex==-2) then
          leos_localisothermal=.true.
!
          call farray_register_global('cs2',iglobal_cs2)
          call farray_register_global('glnTT',iglobal_glnTT,vector=3)
!
        elseif (findex<0) then
          leos_isothermal=.true.
        endif
      elseif (variable=='lnTT') then
        this_var=ieosvar_lnTT
        if (findex<0) leos_isothermal=.true.
      elseif (variable=='TT') then
        this_var=ieosvar_TT
      elseif (variable=='lnrho') then
        this_var=ieosvar_lnrho
        if (findex<0) leos_isochoric=.true.
      elseif (variable=='rho') then
        this_var=ieosvar_rho
        if (findex<0) leos_isochoric=.true.
      elseif (variable=='pp') then
        this_var=ieosvar_pp
        if (findex<0) leos_isobaric=.true.
      else
        call fatal_error("select_eos_variable","no such thermodynamic variable: "//trim(variable))
      endif
      if (ieosvar_count==1) then
        ieosvar1=findex
        ieosvar_selected=ieosvar_selected+this_var
        return
      endif
!
! Ensure the indexes are in the correct order.
!
      if (this_var<ieosvar_selected) then
        ieosvar2=ieosvar1
        ieosvar1=findex
      else
        ieosvar2=findex
      endif
      ieosvar_selected=ieosvar_selected+this_var
      select case (ieosvar_selected)
        case (ieosvar_lnrho+ieosvar_ss)
          if (lroot) print*,"select_eos_variable: Using lnrho and ss"
          ieosvars=ilnrho_ss
        case (ieosvar_rho+ieosvar_ss)
          if (lroot) print*,"select_eos_variable: Using rho and ss"
          ieosvars=irho_ss
        case (ieosvar_lnrho+ieosvar_lnTT)
          if (lroot) print*,"select_eos_variable: Using lnrho and lnTT"
          ieosvars=ilnrho_lnTT
        case (ieosvar_lnrho+ieosvar_TT)
          if (lroot) print*,"select_eos_variable: Using lnrho and TT"
          ieosvars=ilnrho_TT
        case (ieosvar_rho+ieosvar_lnTT)
          if (lroot) print*,"select_eos_variable: Using rho and lnTT"
          ieosvars=irho_lnTT
        case (ieosvar_lnrho+ieosvar_cs2)
          if (lroot) print*,"select_eos_variable: Using lnrho and cs2"
          ieosvars=ilnrho_cs2
        case (ieosvar_rho+ieosvar_cs2)
          if (lroot) print*,"select_eos_variable: Using rho and cs2",iproc
          ieosvars=irho_cs2
        case (ieosvar_rho+ieosvar_TT)
          if (lroot) print*, 'select_eos_variable: Using rho and TT'
          ieosvars=irho_TT
        case default
          call not_implemented("select_eos_variable", &
               "thermodynamic variable combination ieosvar_selected="//trim(itoa(ieosvar_selected)))
      endselect
!
    endsubroutine select_eos_variable
!***********************************************************************
   subroutine getmu(f,mu_tmp)
!
!  Calculate average particle mass in the gas relative to
!
!   12-aug-03/tony: implemented
!   23 may-10/nils: fleshed it up
!
      real, dimension (mx,my,mz,mfarray), optional :: f
      real, optional, intent(out) :: mu_tmp
!
      call keep_compiler_quiet(mu_tmp)
      call keep_compiler_quiet(present(f))

      call not_implemented("getmu","in eos_chemistry_simple")
!
    endsubroutine getmu
!***********************************************************************
    subroutine rprint_eos(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      if (lwrite_slices) then
        if (lpres_grad) where(cnamev=='gpx'.or.cnamev=='gpy') cformv='DEFINED'
        where(cnamev=='pp') cformv='DEFINED'
      endif
!
    endsubroutine rprint_eos
!***********************************************************************
    subroutine get_slices_eos(f,slices)
!
!  Write slices for animation of Eos variables.
!
!  26-jul-06/tony: coded

      use Slices_methods, only: assign_slices_scal
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices.
!
      select case (trim(slices%name))
!
        case ('lnTT'); call assign_slices_scal(slices,f,ilnTT)
        case ('cp'); call assign_slices_scal(slices,f,icp)
        case ('viscosity'); call assign_slices_scal(slices,f,iviscosity)
        case ('gpx'); if (lpres_grad) call assign_slices_scal(slices,f,igpx)
        case ('gpy'); if (lpres_grad) call assign_slices_scal(slices,f,igpy)
        case ('pp')
          if (linterp_pressure) then
            call assign_slices_scal(slices,f,ipp)
          elseif (ldensity_nolog .or. ltemperature_nolog) then
            if (lwrite_slice_yz) slices%yz=(f(ix_loc,m1:m2,n1:n2,iTT)*f(ix_loc,m1:m2,n1:n2,irho)) &
                                           *f(ix_loc,m1:m2,n1:n2,iRR)
            if (lwrite_slice_xz) slices%xz=(f(l1:l2,iy_loc,n1:n2,iTT)*f(l1:l2,iy_loc,n1:n2,irho)) &
                                           *f(l1:l2,iy_loc,n1:n2,iRR)
            if (lwrite_slice_xz2) slices%xz=(f(l1:l2,iy2_loc,n1:n2,iTT)*f(l1:l2,iy2_loc,n1:n2,irho)) &
                                            *f(l1:l2,iy2_loc,n1:n2,iRR)
            if (lwrite_slice_xy) slices%xy=(f(l1:l2,m1:m2,iz_loc,iTT)*f(l1:l2,m1:m2,iz_loc,irho)) &
                                           *f(l1:l2,m1:m2,iz_loc,iRR)
            if (lwrite_slice_xy2) slices%xy2=(f(l1:l2,m1:m2,iz2_loc,iTT)*f(l1:l2,m1:m2,iz2_loc,irho)) &
                                             *f(l1:l2,m1:m2,iz2_loc,iRR)
            if (lwrite_slice_xy3) slices%xy3=(f(l1:l2,m1:m2,iz3_loc,iTT)*f(l1:l2,m1:m2,iz3_loc,irho)) &
                                             *f(l1:l2,m1:m2,iz3_loc,iRR)
            if (lwrite_slice_xy4) slices%xy4=(f(l1:l2,m1:m2,iz4_loc,iTT)*f(l1:l2,m1:m2,iz4_loc,irho)) &
                                             *f(l1:l2,m1:m2,iz4_loc,iRR)
          else
            if (lwrite_slice_yz) slices%yz=exp(f(ix_loc,m1:m2,n1:n2,ilnTT)+f(ix_loc,m1:m2,n1:n2,ilnrho)) &
                                              *f(ix_loc,m1:m2,n1:n2,iRR)
            if (lwrite_slice_xz) slices%xz=exp(f(l1:l2,iy_loc,n1:n2,ilnTT)+f(l1:l2,iy_loc,n1:n2,ilnrho)) &
                                              *f(l1:l2,iy_loc,n1:n2,iRR)
            if (lwrite_slice_xz2) slices%xz=exp(f(l1:l2,iy2_loc,n1:n2,ilnTT)+f(l1:l2,iy2_loc,n1:n2,ilnrho)) &
                                               *f(l1:l2,iy2_loc,n1:n2,iRR)
            if (lwrite_slice_xy) slices%xy=exp(f(l1:l2,m1:m2,iz_loc,ilnTT)+f(l1:l2,m1:m2,iz_loc,ilnrho)) &
                                              *f(l1:l2,m1:m2,iz_loc,iRR)
            if (lwrite_slice_xy2) slices%xy2=exp(f(l1:l2,m1:m2,iz2_loc,ilnTT)+f(l1:l2,m1:m2,iz2_loc,ilnrho)) &
                                                *f(l1:l2,m1:m2,iz2_loc,iRR)
            if (lwrite_slice_xy3) slices%xy3=exp(f(l1:l2,m1:m2,iz3_loc,ilnTT)+f(l1:l2,m1:m2,iz3_loc,ilnrho)) &
                                                *f(l1:l2,m1:m2,iz3_loc,iRR)
            if (lwrite_slice_xy4) slices%xy4=exp(f(l1:l2,m1:m2,iz4_loc,ilnTT)+f(l1:l2,m1:m2,iz4_loc,ilnrho)) &
                                                *f(l1:l2,m1:m2,iz4_loc,iRR)
          endif
        slices%ready=.true.
!
      endselect
!
    endsubroutine get_slices_eos
!***********************************************************************
   subroutine pencil_criteria_eos
!
!  All pencils that the EquationOfState module depends on are specified here.
!
!  02-04-06/tony: coded
!
!  EOS is a pencil provider but evolves nothing so it is unlokely that
!  it will require any pencils for it's own use.
!
      lpenc_requested(i_cv) = .true.
      lpenc_requested(i_cp) = .true.
      lpenc_requested(i_cv1) = .true.
      lpenc_requested(i_cp1) = .true.
!
      lpenc_requested(i_nu) = .true.
      lpenc_requested(i_gradnu) = .true.
      lpenc_requested(i_lambda) = .true.
      lpenc_requested(i_glambda) = .true.
      lpenc_requested(i_glncp) = .true.

      lpenc_requested(i_TT)=.true.
      lpenc_requested(i_TT1)=.true.
      if (.not. ldustdensity) then
        lpenc_requested(i_lnTT)=.true.
      endif
      lpenc_requested(i_del2lnTT)=.true.
!
      if (ltemperature_nolog) then
        lpenc_requested(i_gTT)=.true.
      endif
      lpenc_requested(i_glnrho)=.true.
!
      lpenc_requested(i_rho1gpp)=.true.
      lpenc_requested(i_RRmix)=.true.
      lpenc_requested(i_glnRR)=.true.
      lpenc_requested(i_pp)=.true.
      lpenc_requested(i_glnTT)=.true.
!
      lpenc_requested(i_cs2) = .true.
!
    endsubroutine pencil_criteria_eos
!***********************************************************************
    subroutine pencil_interdep_eos(lpencil_in)
!
!  Interdependency among pencils from the Entropy module is specified here.
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_cv1)) lpencil_in(i_cv) = .true.
      if (lpencil_in(i_cp1)) lpencil_in(i_cp) = .true.
!
    endsubroutine pencil_interdep_eos
!***********************************************************************
    subroutine calc_pencils_eos_std(f,p)
!
! Envelope adjusting calc_pencils_eos_pencpar to the standard use with
! lpenc_loc=lpencil
!
!  9-oct-15/MR: coded
!
      real, dimension (mx,my,mz,mfarray),intent(INOUT):: f
      type (pencil_case),                intent(INOUT):: p
!
      call calc_pencils_eos_pencpar(f,p,lpencil)
!
    endsubroutine calc_pencils_eos_std
!***********************************************************************
    subroutine calc_pencils_eos_pencpar(f,p,lpenc_loc)
!
!  Calculate Entropy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  02-apr-06/tony: coded
!  09-oct-15/MR: added mask parameter lpenc_loc.
!
      use Sub, only: grad, del2, dot2
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc
!
      intent(in) :: lpenc_loc
      intent(inout) :: f
      intent(out) :: p

      real, dimension(nx,3) :: glnDiff_full_add, glncp
      real, dimension(nx) :: D_th, R_mix
      integer :: i,k,j2,j3
!
! Cp/Cv pencils
!
          if (lpencil(i_cp)) then
            p%cp = f(l1:l2,m,n,icp)
            if (lpencil(i_cp1)) p%cp1 = 1./p%cp
          endif
!
          if (lpencil(i_cv))  then
            p%cv = 0.
            if (Cp_const < impossible) then
              do k = 1,nchemspec
                p%cv = p%cv + (Cp_const-Rgas)/species_constants(k,imass)*f(l1:l2,m,n,ichemspec(k))
              enddo
            else
              R_mix = 0.
              do k = 1,nchemspec
                R_mix = R_mix + Rgas/species_constants(k,imass)*f(l1:l2,m,n,ichemspec(k))
              enddo
              p%cv = p%cp - R_mix
            endif
            if (lpencil(i_cv1)) p%cv1 = 1./p%cv
          endif
!
!  Temperature
!
       if (lpenc_loc(i_lnTT)) then
         if (ltemperature_nolog) then
          p%lnTT=log(f(l1:l2,m,n,iTT))
         else
          p%lnTT=f(l1:l2,m,n,ilnTT)
         endif
       endif
!
       if (lpenc_loc(i_TT))  then
         if (ltemperature_nolog) then
           p%TT=f(l1:l2,m,n,iTT)
         else
           p%TT=exp(f(l1:l2,m,n,ilnTT))
         endif
!
         if (minval(p%TT)==0.) call fatal_error('calc_pencils_eos','p%TT=0')
       endif
!
       if (lpenc_loc(i_TT1)) then
         if (ltemperature_nolog) then
           p%TT1=1./f(l1:l2,m,n,iTT)
         else
           p%TT1=1./exp(f(l1:l2,m,n,ilnTT))
         endif
       endif
!
!  Temperature laplacian and gradient
!
      if (lpenc_loc(i_gTT)) call grad(f,iTT,p%gTT)
      if (lpenc_loc(i_glnTT)) then
        if (ltemperature_nolog) then
          p%glnTT(:,1)=p%gTT(:,1)/p%TT(:)
          p%glnTT(:,2)=p%gTT(:,2)/p%TT(:)
          p%glnTT(:,3)=p%gTT(:,3)/p%TT(:)
        else
          call grad(f,ilnTT,p%glnTT)
        endif
      endif
!
      if (lpenc_loc(i_del2lnTT)) then
        if (.not. ltemperature_nolog) call del2(f,ilnTT,p%del2lnTT)
      endif
!
! Viscosity of a mixture
!
      if (lpencil(i_nu)) then
        p%nu = f(l1:l2,m,n,iviscosity)
        if (lpencil(i_gradnu)) call grad(f(:,:,:,iviscosity),p%gradnu)
      endif
!
! Calculate thermal conductivity & diffusivity
!
      D_th = f(l1:l2,m,n,iviscosity)/Pr_number
      if (lpencil(i_lambda)) p%lambda = p%cp*p%rho*D_th
      if (lpencil(i_glambda)) then
         call grad(f(:,:,:,icp),p%glncp)
         do i = 1,3
           p%glncp(:,i) = p%glncp(:,i)*p%cp1
           glnDiff_full_add(:,i) = p%gradnu(:,i)/p%nu
           p%glambda(:,i) = p%lambda*(p%glnrho(:,i)+glnDiff_full_add(:,i) + p%glncp(:,i))
         enddo
      endif
!
!  Mean molecular weight
!
        if (lpenc_loc(i_RRmix)) p%RRmix = f(l1:l2,m,n,iRR)
        if (lpenc_loc(i_glnRR)) then
          call grad(f(:,:,:,iRR),p%glnRR)
          do i=1,3
            p%glnRR(:,i) = p%glnRR(:,i)/p%RRmix
          enddo
        endif
!
!  Pressure
!
      if (lpenc_loc(i_pp)) p%pp = p%TT*p%rho*p%RRmix
      if (linterp_pressure) then
        f(l1:l2,m,n,ipp) = p%pp
      endif
!
!  Logarithmic pressure gradient
!
      if (lpenc_loc(i_rho1gpp)) then
        do i=1,3
          p%rho1gpp(:,i) = p%pp*p%rho1(:)*(p%glnrho(:,i)+p%glnTT(:,i)+p%glnRR(:,i))
        enddo
      endif
      if (lpres_grad) then
        f(l1:l2,m,n,igpx) = p%rho1gpp(:,1)*f(l1:l2,m,n,irho)
        f(l1:l2,m,n,igpy) = p%rho1gpp(:,2)*f(l1:l2,m,n,irho)
      endif
!
!  Energy per unit mass (this is not used now)
!
      if (lpencil(i_ee)) p%ee = p%cv*p%TT
!
! Gradient of lnpp (removed)
!
! This is not used now since advec_cs2 is not computed but will probably be needed
!
      if (lpencil(i_cs2)) then
        if (any(p%cv1 == 0.0)) then
        else
          p%cs2 = p%cp*p%cv1*p%TT*p%RRmix
        endif
      endif
!
    endsubroutine calc_pencils_eos_pencpar
!***********************************************************************
    subroutine read_eos_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=eos_init_pars, IOSTAT=iostat)
!
    endsubroutine read_eos_init_pars
!***********************************************************************
    subroutine write_eos_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=eos_init_pars)
!
    endsubroutine write_eos_init_pars
!***********************************************************************
    subroutine read_eos_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=eos_run_pars, IOSTAT=iostat)
!
    endsubroutine read_eos_run_pars
!***********************************************************************
    subroutine write_eos_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=eos_run_pars)
!
    endsubroutine write_eos_run_pars
!***********************************************************************
!**************ONLY DUMMY ROUTINES BELOW********************************
!***********************************************************************
    subroutine isothermal_entropy(lnrho,T0,ss)
!
      real, dimension (mx,my,mz), intent(in ) :: lnrho
      real, dimension (mx,my,mz), intent(out) :: ss
      real, intent(in) :: T0

      call not_implemented("isothermal_entropy","in eos_chemistry_simple")

      call keep_compiler_quiet(T0)
      call keep_compiler_quiet(lnrho,ss)

    endsubroutine isothermal_entropy
!***********************************************************************
    subroutine bc_ss_flux(f,topbot,lone_sided)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  26-aug-2003/tony: distributed across ionization modules
!   3-oct-16/MR: added new optional switch lone_sided
!
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      logical, optional :: lone_sided
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux
!***********************************************************************
    subroutine bc_ss_flux_turb(f,topbot)
!
!  dummy routine
!
!   4-may-2009/axel: dummy routine
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_turb
!***********************************************************************
    subroutine bc_ss_flux_turb_x(f,topbot)
!
!  dummy routine
!
!   31-may-2010/pete: dummy routine
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_turb_x
!***********************************************************************
    subroutine bc_ss_flux_condturb_x(f,topbot)
!
!   23-apr-2014/pete: dummy
!
      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_condturb_x
!***********************************************************************
    subroutine bc_ss_flux_condturb_mean_x(f,topbot)
!
!   07-jan-2015/pete: dummy
!
      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_condturb_mean_x
!***********************************************************************
    subroutine bc_ss_flux_condturb_z(f,topbot)
!
!   15-jul-2014/pete: dummy
!
      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_condturb_z
!***********************************************************************
    subroutine bc_ss_temp_old(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  23-jun-2003/tony: implemented for leos_fixed_ionization
!  26-aug-2003/tony: distributed across ionization modules
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_old
!***********************************************************************
 subroutine bc_ss_temp_x(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_x
!***********************************************************************
    subroutine bc_ss_temp_y(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_y
!***********************************************************************
    subroutine bc_ss_temp_z(f,topbot,lone_sided)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
      logical, optional :: lone_sided
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_z
!***********************************************************************
    subroutine bc_lnrho_temp_z(f,topbot)
!
!  boundary condition for lnrho *and* ss: constant temperature
!
!  27-sep-2002/axel: coded
!  19-aug-2005/tobi: distributed across ionization modules
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_temp_z
!***********************************************************************
    subroutine bc_lnrho_pressure_z(f,topbot)
!
!  boundary condition for lnrho: constant pressure
!
!   4-apr-2003/axel: coded
!   1-may-2003/axel: added the same for top boundary
!  19-aug-2005/tobi: distributed across ionization modules
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_pressure_z
!***********************************************************************
    subroutine bc_ss_temp2_z(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp2_z
!***********************************************************************
    subroutine bc_ss_temp3_z(f,topbot)
!
!  31-jan-2013/axel: coded to impose cs2bot and dcs2bot at bottom
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_temp3_z","in eos_chemistry_simple")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp3_z
!***********************************************************************
    subroutine bc_ss_stemp_x(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_stemp_x
!***********************************************************************
   subroutine bc_ss_stemp_y(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_stemp_y
!***********************************************************************
    subroutine bc_ss_stemp_z(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_stemp_z
!***********************************************************************
    subroutine bc_ss_a2stemp_x(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_a2stemp_x
!***********************************************************************
   subroutine bc_ss_a2stemp_y(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_a2stemp_y
!***********************************************************************
    subroutine bc_ss_a2stemp_z(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_a2stemp_z
!***********************************************************************
    subroutine bc_ss_energy(f,topbot)
!
!  boundary condition for entropy
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the entropy module
!  26-aug-2003/tony: distributed across ionization modules
!
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_energy
!***********************************************************************
    subroutine bc_stellar_surface(f,topbot)
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_stellar_surface
!***********************************************************************
    subroutine bc_lnrho_cfb_r_iso(f,topbot)
!
!  Boundary condition for radial centrifugal balance
!
!  21-aug-2006/wlad: coded
!
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_cfb_r_iso
!***********************************************************************
    subroutine bc_lnrho_hds_z_iso(f,topbot)
!
!  Boundary condition for density *and* entropy.
!
!  12-Juil-2006/dintrans: coded
!
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hds_z_iso
!***********************************************************************
    subroutine bc_lnrho_hdss_z_iso(f,topbot)
!
!  Smooth out density perturbations with respect to hydrostatic
!  stratification in Fourier space.
!  05-jul-07/tobi: Adapted from bc_aa_pot3
!
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent(IN) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hdss_z_iso
!***********************************************************************
    subroutine bc_ism(f,topbot,j)
!
!  30-nov-15/fred: Replaced bc_ctz and bc_cdz.
!  Apply observed scale height locally from Reynolds 1991, Manchester & Taylor
!  1981 for warm ionized gas - dominant scale height above 500 parsecs.
!  Apply constant local temperature across boundary for entropy.
!  Motivation to prevent numerical spikes in shock fronts, which cannot be
!  absorbed in only three ghost cells, but boundary thermodynamics still
!  responsive to interior dynamics.
!  06-jun-22/fred update to allow setting scale height in start.in or run.in
!  default is density_scale_factor=impossible so that scale_factor is 0.9, assuming
!  unit_length = 1 kpc and scale is 400 pc. To change scale height add to
!  start_pars or run_pars density_scale_factor=... in dimensionless units
!  Copied from eos_ionization written for entropy - may need revision
!  Currently not correct for energy variable
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j,k
      real :: density_scale1, density_scale
!
      if (density_scale_factor==impossible) then
        density_scale=density_scale_cgs/unit_length
      else
        density_scale=density_scale_factor
      endif
      density_scale1=1./density_scale
!
      select case (topbot)
!
      case(BOT)               ! bottom boundary
        do k=1,nghost
          if (j==irho .or. j==ilnrho) then
            if (ldensity_nolog) then
              f(:,:,k,j)=f(:,:,n1,j)*exp(-(z(n1)-z(k))*density_scale1)
            else
              f(:,:,k,j)=f(:,:,n1,j) - (z(n1)-z(k))*density_scale1
            endif
          else if (j==iss) then
            f(:,:,k,j)=f(:,:,n1,j) + (z(n1)-z(k))*density_scale1
            !if (ldensity_nolog) then
            !  f(:,:,n1-k,j)=f(:,:,n1,j)+(cp-cv) * &
            !      (log(f(:,:,n1,j-1))-log(f(:,:,n1-k,j-1))) + &
            !      cv*log((z(n1)-z(n1-k))*density_scale+1.)
            !else
            !  f(:,:,n1-k,j)=f(:,:,n1,j)+(cp-cv)*&
            !      (f(:,:,n1,j-1)-f(:,:,n1-k,j-1))+&
            !      cv*log((z(n1)-z(n1-k))*density_scale+1.)
            !endif
          else
            call fatal_error('bc_ism','only for irho, ilnrho or iss')
          endif
        enddo
!
      case(TOP)               ! top boundary
        do k=1,nghost
          if (j==irho .or. j==ilnrho) then
            if (ldensity_nolog) then
              f(:,:,n2+k,j)=f(:,:,n2,j)*exp(-(z(n2+k)-z(n2))*density_scale1)
            else
              f(:,:,n2+k,j)=f(:,:,n2,j) - (z(n2+k)-z(n2))*density_scale1
            endif
          else if (j==iss) then
            f(:,:,n2+k,j)=f(:,:,n2,j) + (z(n2+k)-z(n2))*density_scale1
            !if (ldensity_nolog) then
            !  f(:,:,n2+k,j)=f(:,:,n2,j)+(cp-cv)*&
            !      (log(f(:,:,n2,j-1))-log(f(:,:,n2+k,j-1)))+&
            !      cv*log((z(n2+k)-z(n2))*density_scale+1.)
            !else
            !  f(:,:,n2+k,j)=f(:,:,n2,j)+(cp-cv)*&
            !      (f(:,:,n2,j-1)-f(:,:,n2+k,j-1))+&
            !      cv*log((z(n2+k)-z(n2))*density_scale+1.)
            !endif
          else
            call fatal_error('bc_ism','only for irho, ilnrho or iss')
          endif
        enddo
!
      case default
        print*, "bc_ism: topbot should be BOT or TOP"
!
      endselect
!
    endsubroutine bc_ism
!***********************************************************************
    subroutine get_gamma_etc(gamma,cp,cv)
!
      real, optional, intent(OUT) :: gamma, cp,cv
!
      if (headt) call warning('get_gamma_etc','gamma, cp, and cv are not constant in eos_chemistry_simple.'// &
                              achar(10)//'The values provided are for one-atomic ideal gas. Use at own risk')
      if (present(gamma)) gamma=5./3.
      if (present(cp)) cp=1.
      if (present(cv)) cv=3./5.

    endsubroutine get_gamma_etc
!***********************************************************************
    subroutine pushpars2c(p_par)

    integer, parameter :: n_pars=0
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call keep_compiler_quiet(p_par)

    endsubroutine pushpars2c
!***********************************************************************
!********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING        *************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include 'eos_dummies.inc'
!***********************************************************************
endmodule EquationOfState
