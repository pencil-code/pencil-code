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
! PENCILS PROVIDED lnTT; gss(3); glnTT(3); TT; TT1; gTT(3)
! PENCILS PROVIDED pp; del2pp; mu1; gmu1(3); glnmu(3)
! PENCILS PROVIDED rho1gpp(3); glnpp(3); del2lnTT
!
! PENCILS PROVIDED hss(3,3); hlnTT(3,3); del2ss; del6ss; del6lnTT
! PENCILS PROVIDED yH; ee; ss; delta; glnmumol(3); ppvap; csvap2
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
  logical :: lpp_as_aux=.false.
  real :: gamma=impossible
  real :: Rgas_cgs=0., Rgas, Rgas_unit_sys=1., error_cp=1e-6
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
!
  real :: Cp_const=impossible
  real :: Pr_number=0.7
  logical :: lpres_grad = .false.
!
!NILS: Why do we spend a lot of memory allocating these variables here????
!MR: Is now allocated only once.
 real, dimension(mx,my,mz), target :: mu1_full
 integer :: imass=1
!
  namelist /eos_init_pars/ mu, cp, cs0, rho0, gamma, error_cp, lpp_as_aux
!
  namelist /eos_run_pars/  mu, cs0, rho0, lpp_as_aux
!
  contains
!***********************************************************************
    subroutine register_eos
!
!  14-jun-03/axel: adapted from register_eos
!
      use SharedVariables, only: put_shared_variable
      use Sub, only: register_report_aux
!
      ilnTT = 0
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          '$Id$')
!
!  pressure as optional auxiliary variable
!
      if (lpp_as_aux) call register_report_aux('pp',ipp)
!
      call put_shared_variable('gamma',gamma,caller='register_eos')

      if (.not.ldensity) then
        call put_shared_variable('rho0',rho0)
        call put_shared_variable('lnrho0',lnrho0)
      endif
      if (lchemistry) call put_shared_variable('mu1_full',mu1_full)

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
        call fatal_error('units_eos','unit_temperature not found')
      else
        Rgas=Rgas_unit_sys*unit_temperature/unit_velocity**2
      endif
!
      inquire(file='chem.inp',exist=cheminp)
      inquire(file='chem.in',exist=chemin)
      if (chemin .and. cheminp) call fatal_error('eos_chemistry', &
          'both chem.inp and chem.in found. Please decide for one')
!
      if (cheminp) input_file='chem.inp'
      if (chemin) input_file='chem.in'
      lcheminp_eos = cheminp .or. chemin
!      inquire(FILE=input_file, EXIST=lcheminp_eos)
!
      if (.not. lcheminp_eos ) then
        call fatal_error('units_eos','file chem.imp not found')
      elseif (lroot) then
        print*,'units_eos: chem.imp is found! Now cp, cv, gamma, mu are pencils ONLY!'
      endif
!
    endsubroutine units_eos
!***********************************************************************
    subroutine initialize_eos(f)
!
! Initialize variable selection code (needed for RELOADing)
!
      real, dimension (mx,my,mz,mfarray) :: f

      ieosvars=-1
      ieosvar_count=0
!
!  write constants to disk. In future we may want to deal with this
!  using an include file or another module.
!
      if (pretend_lnTT) then
        call warning('initialize_eos','pretend_lnTT is not used with ionization')
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
      if (dimensionality==0) then
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
           "2 thermodynamic quantities have already been defined while attempting to add a 3rd: "// &
           trim(variable))
!
      ieosvar_count=ieosvar_count+1
!
!      select case (variable)
      if (variable=='ss') then
          this_var=ieosvar_ss
          if (findex<0) then
            leos_isentropic=.true.
          endif
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
          if (findex<0) then
            leos_isothermal=.true.
          endif
      elseif (variable=='TT') then
          this_var=ieosvar_TT
      elseif (variable=='lnrho') then
          this_var=ieosvar_lnrho
          if (findex<0) then
            leos_isochoric=.true.
          endif
      elseif (variable=='rho') then
          this_var=ieosvar_rho
          if (findex<0) then
            leos_isochoric=.true.
          endif
      elseif (variable=='pp') then
          this_var=ieosvar_pp
          if (findex<0) then
            leos_isobaric=.true.
          endif
      else
        call fatal_error("select_eos_variable", &
             "unknown thermodynamic variable")
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
        case default
          call not_implemented("select_eos_variable", &
             "thermodynamic variable combination ieosvar_selected="//trim(itoa(ieosvar_selected)))
      endselect
!
    endsubroutine select_eos_variable
!***********************************************************************
    subroutine rprint_eos(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(present(lwrite))
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='lnTT'.or.cnamev=='pp') cformv='DEFINED'
      endif
!
    endsubroutine rprint_eos
!***********************************************************************
    subroutine get_slices_eos(f,slices)
!
!  Write slices for animation of Eos variables.
!
!  26-jul-06/tony: coded
!
      use Slices_methods, only: assign_slices_scal
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices.
!
      select case (trim(slices%name))
!
!  Temperature.
!
        case ('lnTT'); call assign_slices_scal(slices,f,ilnTT)
!
!  Pressure.
!
        case ('pp')
          if (ldensity_nolog .or. ltemperature_nolog) &
              call not_implemented('get_slices_eos','pressure slice for ldensity_nolog or ltemperature_nolog')
          if (lwrite_slice_yz) slices%yz=Rgas*exp(f(ix_loc,m1:m2,n1:n2,ilnTT)+f(ix_loc,m1:m2,n1:n2,ilnrho)) &
                                             *mu1_full(ix_loc,m1:m2,n1:n2)
          if (lwrite_slice_xz) slices%xz=Rgas*exp(f(l1:l2,iy_loc,n1:n2,ilnTT)+f(l1:l2,iy_loc,n1:n2,ilnrho)) &
                                             *mu1_full(l1:l2,iy_loc,n1:n2)
          if (lwrite_slice_xz2) slices%xz=Rgas*exp(f(l1:l2,iy2_loc,n1:n2,ilnTT)+f(l1:l2,iy2_loc,n1:n2,ilnrho)) &
                                              *mu1_full(l1:l2,iy2_loc,n1:n2)
          if (lwrite_slice_xy) slices%xy=Rgas*exp(f(l1:l2,m1:m2,iz_loc,ilnTT)+f(l1:l2,m1:m2,iz_loc,ilnrho)) &
                                             *mu1_full(l1:l2,m1:m2,iz_loc)
          if (lwrite_slice_xy2) slices%xy2=Rgas*exp(f(l1:l2,m1:m2,iz2_loc,ilnTT)+f(l1:l2,m1:m2,iz2_loc,ilnrho)) &
                                               *mu1_full(l1:l2,m1:m2,iz2_loc)
          if (lwrite_slice_xy3) slices%xy3=Rgas*exp(f(l1:l2,m1:m2,iz3_loc,ilnTT)+f(l1:l2,m1:m2,iz3_loc,ilnrho)) &
                                               *mu1_full(l1:l2,m1:m2,iz3_loc)
          if (lwrite_slice_xy4) slices%xy4=Rgas*exp(f(l1:l2,m1:m2,iz4_loc,ilnTT)+f(l1:l2,m1:m2,iz4_loc,ilnrho)) &
                                               *mu1_full(l1:l2,m1:m2,iz4_loc)
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
      lpenc_requested(i_glnTT)=.true.
      lpenc_requested(i_glnrho)=.true.
      lpenc_requested(i_glnrho2)=.true.
      lpenc_requested(i_del2lnrho)=.true.
!
      lpenc_requested(i_glnpp)=.true.
      lpenc_requested(i_del2pp)=.true.
      lpenc_requested(i_mu1)=.true.
      lpenc_requested(i_gmu1)=.true.
      lpenc_requested(i_pp)=.true.
      lpenc_requested(i_rho1gpp)=.true.
!
!  pp pencil if lpp_as_aux
!
      if (lpp_as_aux) lpenc_requested(i_pp)=.true.
!
    endsubroutine pencil_criteria_eos
!***********************************************************************
    subroutine pencil_interdep_eos(lpencil_in)
!
!  Interdependency among pencils from the Entropy module is specified here.
!
!  20-11-04/anders: coded
!
! Modified by Natalia. Taken from  eos_temperature_ionization module
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_TT1))    lpencil_in(i_TT)=.true.
!
      if (lpencil_in(i_pp))  then
         lpencil_in(i_mu1)=.true.
         lpencil_in(i_rho)=.true.
         lpencil_in(i_TT)=.true.
       endif
       if (lpencil_in(i_rho1gpp))  then
         lpencil_in(i_mu1)=.true.
         lpencil_in(i_gmu1)=.true.
         lpencil_in(i_pp)=.true.
         lpencil_in(i_rho)=.true.
         lpencil_in(i_glnTT)=.true.
       endif
       if (lpencil_in(i_glnpp))  then
         lpencil_in(i_pp)=.true.
         lpencil_in(i_rho)=.true.
         lpencil_in(i_rho1gpp)=.true.
       endif
!
       if (lpencil_in(i_ee)) then
         lpencil_in(i_mu1)=.true.
         lpencil_in(i_TT)=.true.
       endif
!
      call keep_compiler_quiet(lpencil_in)
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
      real, dimension(nx) :: glnTT2,TT1del2TT,del2lnrho
      real, dimension(nx) :: rho1del2rho,del2mu1,gmu12,glnpp2
      real, dimension(nx) :: del2TT, gradTgradT
!
      intent(in) :: lpenc_loc
      intent(inout) :: f,p
!
      integer :: i
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
        if (lpenc_loc(i_glnTT)) then
          if (ltemperature_nolog) then
            call grad(f,iTT,p%glnTT)
            p%glnTT(:,1)=p%glnTT(:,1)/p%TT(:)
            p%glnTT(:,2)=p%glnTT(:,2)/p%TT(:)
            p%glnTT(:,3)=p%glnTT(:,3)/p%TT(:)
          else
            call grad(f,ilnTT,p%glnTT)
          endif
        endif
!
        if (ltemperature_nolog) then
          if (lpenc_loc(i_gTT)) call grad(f,iTT,p%gTT)
            call dot2(p%gTT,gradTgradT)
            call del2(f,iTT,del2TT)
            p%del2lnTT = -p%TT1*p%TT1*gradTgradT+p%TT1*del2TT
         !NILS: The call below does not yield del2lnTT but rather del2TT,
         !NILS: this should be fixed before used. One should also look
         !NILS: through the chemistry module to make sure del2lnTT is used
         !NILS: corectly.
 !         if (lpenc_loc(i_del2lnTT)) call del2(f,iTT,p%del2lnTT)
 !         call fatal_error('calc_pencils_eos',&
 !             'del2lnTT is not correctly implemented - this must be fixed!')
        else
          if (lpenc_loc(i_del2lnTT)) call del2(f,ilnTT,p%del2lnTT)
        endif
!
!  Mean molecular weight
!
        if (lpenc_loc(i_mu1)) p%mu1=mu1_full(l1:l2,m,n)
        if (lpenc_loc(i_gmu1)) call grad(mu1_full,p%gmu1)
        if (lpenc_loc(i_glnmu)) then
          do i = 1,3
            p%glnmu(:,i) = -p%gmu1(:,i)/p%mu1(:)
          enddo
        endif
!
!  Pressure
!
        if (lpenc_loc(i_pp)) p%pp = Rgas*p%TT*p%mu1*p%rho
!
!  Logarithmic pressure gradient
!
        if (lpenc_loc(i_rho1gpp)) then
          do i=1,3
            p%rho1gpp(:,i) = p%pp/p%rho(:)*(p%glnrho(:,i)+p%glnTT(:,i)+p%gmu1(:,i)/p%mu1(:))
          enddo
        endif
!
! Gradient of lnpp
!
       if (lpenc_loc(i_glnpp)) then
         do i=1,3
           p%glnpp(:,i)=p%rho1gpp(:,i)*p%rho(:)/p%pp(:)
         enddo
       endif
!
! Laplacian of pressure
!
       if (lpenc_loc(i_del2pp)) then
         call dot2(p%glnTT,glnTT2)
         TT1del2TT=p%del2lnTT+glnTT2
         rho1del2rho=p%del2lnrho+p%glnrho2
         call dot2(p%gmu1,gmu12)
         call dot2(p%glnpp,glnpp2)
         call del2(mu1_full,del2mu1)
         p%del2pp = p%pp*glnpp2 &
                   +p%pp*(rho1del2rho+TT1del2TT+del2mu1/p%mu1) &
                   -p%pp*(p%glnrho2+glnTT2+gmu12/p%mu1**2)
       endif
!
!  Energy per unit mass (this has been moved to chemistry.f90 in order
!  to get the correct cv.
!
      !if (lpenc_loc(i_ee)) p%ee = p%cv*p%TT
!
!  pressure and cp as optional auxiliary pencils
!
      if (lpp_as_aux) f(l1:l2,m,n,ipp)=p%pp
!
    endsubroutine calc_pencils_eos_pencpar
!***********************************************************************
   subroutine getdensity(f,EE,TT,yH,rho_full)
!
     real, dimension (mx,my,mz,mfarray) :: f
     real, dimension (mx,my,mz), intent(out) :: rho_full
     real, intent(in), optional :: EE,TT,yH
!
      if (ldensity_nolog) then
        rho_full=f(:,:,:,ilnrho)
      else
        rho_full=exp(f(:,:,:,ilnrho))
      endif
!
      call keep_compiler_quiet(present(yH))
      call keep_compiler_quiet(present(EE))
      call keep_compiler_quiet(present(TT))
!
   endsubroutine getdensity
!***********************************************************************
   subroutine gettemperature(f,TT_full)
!
     real, dimension (mx,my,mz,mfarray) :: f
     real, dimension (mx,my,mz), intent(out) :: TT_full
!
      if (ltemperature_nolog) then
        TT_full=f(:,:,:,ilnTT)
      else
        TT_full=exp(f(:,:,:,ilnTT))
      endif
!
   endsubroutine gettemperature
!***********************************************************************
  subroutine getpressure(pp,TT,rho,mu1)
!
     real, dimension (nx), intent(out) :: pp
     real, dimension (nx), intent(in) :: TT, rho, mu1
     integer :: j2,j3
!
     pp=Rgas*mu1*rho*TT
!
   endsubroutine getpressure
!***********************************************************************
    subroutine get_gamma_etc(gamma,cp,cv)
!
      real, optional, intent(OUT) :: gamma, cp,cv
!
      if (headt) call warning('get_gamma_etc','gamma, cp, and cv are not constant in eos_chemistry.'// &
                              achar(10)//'The values provided are for one-atomic ideal gas. Use at own risk')
      if (present(gamma)) gamma=5./3.
      if (present(cp)) cp=1.
      if (present(cv)) cv=3./5.

    endsubroutine get_gamma_etc
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
    subroutine isothermal_entropy(lnrho,T0,ss)
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
      real, dimension (mx,my,mz), intent(in) :: lnrho
      real, dimension (mx,my,mz), intent(out):: ss
      real, intent(in) :: T0
!
      cs2top=cs2bot
      call not_implemented('isothermal_entropy','in eos_chemistry')
!
      call keep_compiler_quiet(lnrho,ss)
      call keep_compiler_quiet(T0)
!
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
      call not_implemented('bc_ss_temp3_z','in eos_chemistry')
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
      call not_implemented("bc_stellar_surface","in eos_chemistry")
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
!  unit_length = 1 kpc and scale is 403 pc. To change scale height add to
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
        call fatal_error('bc_ss_flux','topbot should be BOT or TOP')
!
      endselect
!
    endsubroutine bc_ism
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr

    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(cs20,p_par(1))

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
