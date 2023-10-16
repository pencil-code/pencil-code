! $Id$
!
!  Equation of state for an ideal gas without ionization.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: leos = .true., leos_ionization=.false.
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
  real :: lnTT0=impossible
!
  real :: mu=1.
  real :: cs0=1., rho0=1.
  real :: cs20=1., lnrho0=0.
  real :: gamma=impossible
  real :: Rgas_cgs=0., Rgas, Rgas_unit_sys=1., error_cp=1e-6, scale_Rgas=1.
  real :: cp=impossible
  real :: cs2bot=1., cs2top=1.
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
 real, dimension(nchemspec,18) :: species_constants
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
    subroutine initialize_eos
!
      use SharedVariables, only: put_shared_variable
!
! Initialize variable selection code (needed for RELOADing)
!
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

      call put_shared_variable('linterp_pressure',linterp_pressure)
      call put_shared_variable('scale_Rgas',scale_Rgas)
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
        "2 thermodynamic quantities have already been defined while attempting to add a 3rd: ") !//variable)
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
        call fatal_error("select_eos_variable","unknown thermodynamic variable")
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
          if (lroot) print*,"Thermodynamic variable combination, ieosvar_selected= ",ieosvar_selected
          call fatal_error("select_eos_variable", &
             "This thermodynamic variable combination is not implemented: ")
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
    subroutine find_mass(element_name,MolMass)
!
!  Find mass of element
!
!  05-feb-08/nils: coded
!
      character (len=*), intent(in) :: element_name
      real, intent(out) :: MolMass
!
      select case (element_name)
      case ('H')
        MolMass=1.00794
      case ('C')
        MolMass=12.0107
      case ('N')
        MolMass=14.00674
      case ('O')
        MolMass=15.9994
      case ('Ar','AR')
        MolMass=39.948
      case ('He','HE')
        MolMass=4.0026
      case ('S')
        MolMass=32.0655
      case ('CLOUD')
        MolMass=0.
      case default
        if (lroot) print*,'element_name=',element_name
        call fatal_error('find_mass','no such element: '//trim(element_name))
      end select
!
    endsubroutine find_mass
!***********************************************************************
   subroutine find_species_index(species_name,ind_glob,ind_chem,found_specie)
!
!  Find index in the f array for specie
!
!  05-feb-08/nils: coded
!
      integer, intent(out) :: ind_glob
      integer, intent(inout) :: ind_chem
      character (len=*), intent(in) :: species_name
      integer :: k
      logical, intent(out) :: found_specie
!
      ind_glob=0
    !  ind_chem=0
      do k=1,nchemspec
        if (trim(varname(ichemspec(k)))==species_name) then
          ind_glob=k+ichemspec(1)-1
          ind_chem=k
          exit
        endif
!print*, trim(varname(ichemspec(k))),(species_name)
      enddo
!
!  Check if the species was really found
!
      if ((ind_glob==0)) then
        found_specie=.false.
     !  if (lroot) print*,' no species has been found  ',' species index= ', ind_glob,ind_chem,species_name
     !   call fatal_error('find_species_index','no species has been found')
      else
        found_specie=.true.
    !    if (lroot) print*,species_name,'   species index= ',ind_chem
      endif
!
    endsubroutine find_species_index
!***********************************************************************
    subroutine read_species(input_file)
!
!  This subroutine reads all species information from chem.inp
!  See the chemkin manual for more information on
!  the syntax of chem.inp.
!
!  06-mar-08/nils: coded
!
      logical :: IsSpecie=.false., emptyfile
      integer :: k,file_id=123, StartInd, StopInd
      character (len=80) :: ChemInpLine
      character (len=*) :: input_file
!
      emptyFile=.true.
      k=1
      open(file_id,file=input_file)
      dataloop: do
        read(file_id,'(80A)',end=1000) ChemInpLine(1:80)
        emptyFile=.false.
!
!  Check if we are reading a line within the species section
!
        if (ChemInpLine(1:7)=="SPECIES")            IsSpecie=.true.
        if (ChemInpLine(1:3)=="END" .and. IsSpecie) IsSpecie=.false.
!
!  Read in species
!
        if (IsSpecie) then
          if (ChemInpLine(1:7) /= "SPECIES") then
            StartInd=1; StopInd =0
            stringloop: do
              StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
              if (StopInd==StartInd) then
                StartInd=StartInd+1
              else
                if (k>nchemspec) then
                  print*,'nchemspec=',nchemspec
                  call fatal_error('read_species',"there were too many species, increase nchemspec")
                endif
                varname(ichemspec(k))=trim(ChemInpLine(StartInd:StopInd-1))
                StartInd=StopInd
                k=k+1
              endif
              if (StartInd==80) exit
            enddo stringloop
          endif
        endif
      enddo dataloop
!
!  Stop if chem.inp is empty
!
1000  if (emptyFile)  call fatal_error('read_species','input file chem.inp was empty')
!
!  Check if nchemspec where not too large
!
      if (k<nchemspec-1) then
        print*,'nchemspec=',nchemspec
        call fatal_error("read_species","there were too few species, decrease nchemspec")
      endif
!
      close(file_id)
!
    endsubroutine read_species
!***********************************************************************
   subroutine read_thermodyn(input_file)
!
!  This subroutine reads the thermodynamical data for all species
!  from chem.inp. See the chemkin manual for more information on
!  the syntax of chem.inp.
!
!  06-mar-08/nils: coded
!
      character (len=*), intent(in) :: input_file
      integer :: file_id=123, ind_glob, ind_chem
      character (len=80) :: ChemInpLine
      integer :: In1,In2,In3,In4,In5,iElement,iTemperature,StopInd
      integer :: NumberOfElement_i
      logical :: IsThermo=.false., found_specie, existing_specie
      real, dimension(4) :: MolMass
      real, dimension(3) :: tmp_temp
      character (len=5) :: NumberOfElement_string,element_string
      character (len=10) :: specie_string,TemperatureNr_i
      real :: nne
      integer, dimension(7) :: iaa1,iaa2
!
      integer :: iTemp1=2,iTemp2=3,iTemp3=4
!
      ind_chem=0
!
!  Initialize some index pointers
!
      iaa1(1)=5;iaa1(2)=6;iaa1(3)=7;iaa1(4)=8
      iaa1(5)=9;iaa1(6)=10;iaa1(7)=11
!
      iaa2(1)=12;iaa2(2)=13;iaa2(3)=14;iaa2(4)=15
      iaa2(5)=16;iaa2(6)=17;iaa2(7)=18
!
      open(file_id,file=input_file)
      dataloop2: do
        read(file_id,'(80A)',end=1001) ChemInpLine(1:80)
!
! Check if we are reading a line within the thermo section
!
        if (ChemInpLine(1:6)=="THERMO") IsThermo=.true.
        if (ChemInpLine(1:3)=="END" .and. IsThermo) IsThermo=.false.
!
! Read in thermo data
!
        if (IsThermo) then
          if (ChemInpLine(1:7) /= "THERMO") then
            StopInd=index(ChemInpLine,' ')
            specie_string=trim(ChemInpLine(1:StopInd-1))
!
            call find_species_index(specie_string,ind_glob,ind_chem,found_specie)
!
! Check if we are working with a specie that was found under the SPECIES
! section of chem.inp.
!
            if (ChemInpLine(80:80)=="1") then
              if (found_specie) then
                existing_specie=.true.
              else
                existing_specie=.false.
              endif
            endif
!
! What problems are in the case of  ind_chem=0?
!
            if (ind_chem>0 .and. ind_chem<=nchemspec) then
!
            if (existing_specie) then
            if (found_specie) then
!
! Find molar mass
!
              MolMass=0
              do iElement=1,4
                In1=25+(iElement-1)*5
                In2=26+(iElement-1)*5
                In3=27+(iElement-1)*5
                In4=29+(iElement-1)*5
                if (ChemInpLine(In1:In1)==' ') then
                  MolMass(iElement)=0
                else
                  element_string=trim(ChemInpLine(In1:In2))
                  call find_mass(element_string,MolMass(iElement))
                  In5=verify(ChemInpLine(In3:In4),' ')+In3-1
                  NumberOfElement_string=trim(ChemInpLine(In5:In4))
                  read (unit=NumberOfElement_string,fmt='(I5)') NumberOfElement_i
                  MolMass(iElement)=MolMass(iElement)*NumberOfElement_i
                endif
              enddo
              species_constants(ind_chem,imass)=sum(MolMass)
!
! Find temperature-ranges for low and high temperature fitting
!
              do iTemperature=1,3
                In1=46+(iTemperature-1)*10
                In2=55+(iTemperature-1)*10
                if (iTemperature==3) In2=73
                In3=verify(ChemInpLine(In1:In2),' ')+In1-1
                TemperatureNr_i=trim(ChemInpLine(In3:In2))
                read (unit=TemperatureNr_i,fmt='(F10.1)') nne
                tmp_temp(iTemperature)=nne
              enddo
              species_constants(ind_chem,iTemp1)=tmp_temp(1)
              species_constants(ind_chem,iTemp2)=tmp_temp(3)
              species_constants(ind_chem,iTemp3)=tmp_temp(2)
!
            elseif (ChemInpLine(80:80)=="2") then
              ! Read iaa1(1):iaa1(5)
              read (unit=ChemInpLine(1:75),fmt='(5E15.8)') species_constants(ind_chem,iaa1(1):iaa1(5))
            elseif (ChemInpLine(80:80)=="3") then
              ! Read iaa1(6):iaa5(3)
              read (unit=ChemInpLine(1:75),fmt='(5E15.8)') species_constants(ind_chem,iaa1(6):iaa2(3))
            elseif (ChemInpLine(80:80)=="4") then
              ! Read iaa2(4):iaa2(7)
              read (unit=ChemInpLine(1:75),fmt='(4E15.8)') species_constants(ind_chem,iaa2(4):iaa2(7))
            endif
!
          endif
          endif
          endif !(from ind_chem>0 query)
        endif
      enddo dataloop2
1001  continue
      close(file_id)
!
   endsubroutine read_thermodyn
!***********************************************************************
    subroutine write_thermodyn
!
!  This subroutine writes the thermodynamical data for every specie
!  to ./data/chem.out.
!
!  06-mar-08/nils: coded
!
      use General
!
      character (len=fnlen) :: input_file="./data/chem.out"
      character (len=intlen) :: ispec
      integer :: file_id=123,k
      integer, dimension(7) :: iaa1,iaa2
      integer :: iTemp1=2,iTemp3=4
!
!      Initialize some index pointers
!
      iaa1(1)=5;iaa1(2)=6;iaa1(3)=7;iaa1(4)=8
      iaa1(5)=9;iaa1(6)=10;iaa1(7)=11
!
      iaa2(1)=12;iaa2(2)=13;iaa2(3)=14;iaa2(4)=15
      iaa2(5)=16;iaa2(6)=17;iaa2(7)=18
!
      open(file_id,file=input_file)
      write(file_id,*) 'Specie'
      write(file_id,*) 'MolMass Temp1 Temp2 Temp3'
      write(file_id,*) 'a1(1)  a1(2)  a1(3)  a1(4)  a1(5)  a1(6)  a1(7)'
      write(file_id,*) 'a2(1)  a2(2)  a2(3)  a2(4)  a2(5)  a2(6)  a2(7)'
      write(file_id,*) '***********************************************'
      dataloop2: do k=1,nchemspec
        write(file_id,*) varname(ichemspec(k))
        write(file_id,'(F10.2,3F10.2)') species_constants(k,imass),species_constants(k,iTemp1:iTemp3)
        write(file_id,'(7E12.5)') species_constants(k,iaa1)
        write(file_id,'(7E12.5)') species_constants(k,iaa2)
      enddo dataloop2
!
      close(file_id)
!
      if (lroot) then
        print*,'Write pc_constants.pro in chemistry.f90'
        open (143,FILE=trim(datadir)//'/pc_constants.pro',POSITION="append")
        write (143,*) 'specname=strarr(',nchemspec,')'
        write (143,*) 'specmass=fltarr(',nchemspec,')'
        do k=1,nchemspec
          ispec=itoa(k-1)
          write (143,*) 'specname[',trim(ispec),']=',"'",trim(varname(ichemspec(k))),"'"
          write (143,*) 'specmass[',trim(ispec),']=',species_constants(k,imass)
        enddo
        close (143)
      endif
!
    endsubroutine write_thermodyn
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
   subroutine ioninit(f)
!
!  the ionization fraction has to be set to a value yH0 < yH < yHmax before
!  rtsafe is called for the first time
!
!  12-jul-03/tobi: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine ioninit
!***********************************************************************
    subroutine ioncalc(f)
!
!   calculate degree of ionization and temperature
!   This routine is called from equ.f90 and operates on the full 3-D array.
!
!   13-jun-03/tobi: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine ioncalc
!***********************************************************************
    subroutine pressure_gradient_farray(f,cs2,cp1tilde)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(out) :: cs2,cp1tilde
!
      cs2=impossible
      cp1tilde=impossible
      call not_implemented("pressure_gradient_farray","in eos_chemistry_simple")
!
      call keep_compiler_quiet(f)
!
    endsubroutine pressure_gradient_farray
!***********************************************************************
    subroutine pressure_gradient_point(lnrho,ss,cs2,cp1tilde)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: cs2,cp1tilde
!
      call not_implemented("pressure_gradient_point","in eos_chemistry_simple")
!
      call keep_compiler_quiet(cs2,cp1tilde,ss,lnrho)
!
    endsubroutine pressure_gradient_point
!***********************************************************************
    subroutine temperature_gradient(f,glnrho,gss,glnTT)
!
!   Calculate thermodynamical quantities
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,3), intent(in) :: glnrho,gss
      real, dimension(nx,3), intent(out) :: glnTT
!
      call not_implemented("temperature_gradient","in eos_chemistry_simple")
!
!  given that we just stopped, it cannot become worse by setting
!  cp1tilde to impossible, which allows the compiler to compile.
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(glnrho,gss,glnTT)
!
    endsubroutine temperature_gradient
!***********************************************************************
    subroutine temperature_laplacian(f,del2lnrho,del2ss,del2lnTT)
!
!   Calculate thermodynamical quantities
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(in) :: del2lnrho,del2ss
      real, dimension(nx), intent(out) :: del2lnTT
!
      call not_implemented("temperature_laplacian","in eos_chemistry_simple")
!
     call keep_compiler_quiet(f)
     call keep_compiler_quiet(del2lnrho,del2ss,del2lnTT)
!
    endsubroutine temperature_laplacian
!***********************************************************************
    subroutine temperature_hessian(f,hlnrho,hss,hlnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1
!   and optionally hlnPP and hlnTT
!   hP/rho=cs2*(hlnrho+cp1*hss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,3,3), intent(in) :: hlnrho,hss
      real, dimension(nx,3,3), intent(out) :: hlnTT
!
      call not_implemented("temperature_hessian","in eos_chemistry_simple")
!
      hlnTT(:,:,:)=0
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(hss)
      call keep_compiler_quiet(hlnrho)
!
    endsubroutine temperature_hessian
!***********************************************************************
    subroutine eoscalc_farray(f,psize,lnrho,yH,lnTT,ee,pp,cs2,kapparho)
!
!   dummy routine to calculate thermodynamical quantities
!   copied from eo_idealgas
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: psize
      real, dimension(psize), optional :: lnrho,lnTT
      real, dimension(psize), optional :: yH,ee,pp,cs2,kapparho
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(present(lnrho),present(lnTT))
      call keep_compiler_quiet(present(yH),present(ee))
      call keep_compiler_quiet(present(pp),present(kapparho))
      call keep_compiler_quiet(present(cs2))
!
      call not_implemented("eoscalc_farray","in eos_chemistry_simple")
!
    endsubroutine eoscalc_farray
!***********************************************************************
    subroutine eoscalc_point(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)
!
!   Calculate thermodynamical quantities
!
!   22-jun-06/axel: reinstated cp,cp1,cv,cv1 in hopefully all the places.
!
      integer, intent(in) :: ivars
      real, intent(in) :: var1,var2
      real, optional :: lnrho,ss
      real, optional :: yH,lnTT
      real, optional :: ee,pp,cs2
!
      call not_implemented("eoscalc_point","in eos_chemistry_simple")
!
      call keep_compiler_quiet(present(lnrho),present(lnTT))
      call keep_compiler_quiet(present(pp),present(ee))
      call keep_compiler_quiet(present(yH))
      call keep_compiler_quiet(present(ss),present(cs2))
      call keep_compiler_quiet(ivars)
      call keep_compiler_quiet(var1,var2)
!
    endsubroutine eoscalc_point
!***********************************************************************
    subroutine eoscalc_pencil(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)
!
!   Calculate thermodynamical quantities
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1 to
!                   subroutine pressure_gradient
!   27-mar-06/tony: Introduces cv, cv1, gamma1 to make faster
!                   + more explicit
!   31-mar-06/tony: I removed messy lcalc_cp stuff completely. cp=1.
!                   is just fine.
!   22-jun-06/axel: reinstated cp,cp1,cv,cv1 in hopefully all the places.
!
      integer, intent(in) :: ivars
      real, dimension(nx), intent(in) :: var1,var2
      real, dimension(nx), optional :: lnrho,ss
      real, dimension(nx), optional :: yH,lnTT
      real, dimension(nx), optional :: ee,pp,cs2
!
      call not_implemented("eoscalc_pencil","in eos_chemistry_simple")
!
      call keep_compiler_quiet(present(lnrho),present(lnTT))
      call keep_compiler_quiet(present(pp),present(ee))
      call keep_compiler_quiet(present(yH))
      call keep_compiler_quiet(present(ss),present(cs2))
      call keep_compiler_quiet(ivars)
      call keep_compiler_quiet(var1,var2)
!
    endsubroutine eoscalc_pencil
!***********************************************************************
    subroutine get_soundspeed(TT,cs2)
!
!  Calculate sound speed for given temperature
!
!  20-Oct-03/tobi: Coded
!
      real, intent(in)  :: TT
      real, intent(out) :: cs2
!
      cs2=impossible
      call not_implemented("get_soundspeed","in eos_chemistry_simple")
!
      call keep_compiler_quiet(TT)
!
    endsubroutine get_soundspeed
!***********************************************************************
    subroutine isothermal_entropy(f,T0)
!
!  Isothermal stratification (for lnrho and ss)
!  This routine should be independent of the gravity module used.
!  When entropy is present, this module also initializes entropy.
!
!  Sound speed (and hence Temperature), is
!  initialised to the reference value:
!           sound speed: cs^2_0            from start.in
!           density: rho0 = exp(lnrho0)
!
!  11-jun-03/tony: extracted from isothermal routine in Density module
!                  to allow isothermal condition for arbitrary density
!  17-oct-03/nils: works also with leos_ionization=T
!  18-oct-03/tobi: distributed across ionization modules
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, intent(in) :: T0
!
      cs2top=cs2bot
      call not_implemented("isothermal_entropy","in eos_chemistry_simple")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(T0)
!
    endsubroutine isothermal_entropy
!***********************************************************************
    subroutine isothermal_lnrho_ss(f,T0,rho0)
!
!  Isothermal stratification for lnrho and ss (for yH=0!)
!
!  Currently only implemented for ionization_fixed.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, intent(in) :: T0,rho0
!
      call not_implemented("isothermal_lnrho_ss","in eos_chemistry_simple")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(T0)
      call keep_compiler_quiet(rho0)
!
    endsubroutine isothermal_lnrho_ss
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
    subroutine get_stratz(z, rho0z, dlnrho0dz, eth0z)
!
!  Get background stratification in z direction.
!
!  13-oct-14/ccyang: dummy
!
      real, dimension(:), intent(in) :: z
      real, dimension(:), intent(out), optional :: rho0z, dlnrho0dz, eth0z
!
      call not_implemented('get_stratz', 'stratification for eos_chemistry_simple')
!
      call keep_compiler_quiet(z)
      if (present(rho0z)) call keep_compiler_quiet(rho0z)
      if (present(dlnrho0dz)) call keep_compiler_quiet(dlnrho0dz)
      if (present(eth0z)) call keep_compiler_quiet(eth0z)
!
    endsubroutine get_stratz
!***********************************************************************
   subroutine getdensity(f,EE,TT,yH,rho_full)
!
     real, dimension (mx,my,mz,mfarray) :: f
     real, dimension (mx,my,mz), intent(out) :: rho_full
     real, intent(in), optional :: EE,TT,yH
!
      call keep_compiler_quiet(yH,EE,TT)
      call keep_compiler_quiet(rho_full)
      call keep_compiler_quiet(f)
!
   endsubroutine getdensity
!***********************************************************************
   subroutine gettemperature(f,TT_full)
!
     real, dimension (mx,my,mz,mfarray) :: f
     real, dimension (mx,my,mz), intent(out) :: TT_full
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(TT_full)
!
   endsubroutine gettemperature
!***********************************************************************
  subroutine getpressure(pp,TT,rho,mu1)
!
     real, dimension (nx), intent(out) :: pp
     real, dimension (nx), intent(in) :: TT, rho, mu1
!
      call keep_compiler_quiet(rho,mu1,TT)
      call keep_compiler_quiet(pp)
!
   endsubroutine getpressure
!***********************************************************************
     subroutine get_average_pressure(average_density,average_pressure)
!
!   01-dec-2009/piyali+dhrube: coded
!
      real, intent(in):: average_density
      real, intent(out):: average_pressure
      call keep_compiler_quiet(average_density)
      call keep_compiler_quiet(average_pressure)

    endsubroutine get_average_pressure
!***********************************************************************
    subroutine get_gamma_etc(gamma,cp,cv)
!
      real, intent(OUT) :: gamma
      real, optional, intent(OUT) :: cp,cv
!
      if (headt) call warning('get_gamma_etc','gamma, cp, and cv are not constant in eos_chemistry_simple.'// &
                              achar(10)//'The values provided are for one-atomic ideal gas. Use at own risk')
      gamma=5./3.
      if (present(cp)) cp=1.
      if (present(cv)) cv=3./5.

    endsubroutine get_gamma_etc
!***********************************************************************
    subroutine eosperturb(f,psize,ee,pp,ss)
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(in), optional :: ee,pp,ss
!
      call not_implemented("eosperturb","in eos_chemistry_simple")
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(present(ee),present(pp),present(ss))
!
    endsubroutine eosperturb
!***********************************************************************
    subroutine pushpars2c(p_par)

    integer, parameter :: n_pars=0
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call keep_compiler_quiet(p_par)

    endsubroutine pushpars2c
!***********************************************************************
endmodule EquationOfState
