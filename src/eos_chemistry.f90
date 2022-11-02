! $Id$
!
!  Equation of state for an ideal gas without ionization.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: leos = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED lnTT;  glnTT(3); TT; TT1; gTT(3)
! PENCILS PROVIDED pp; del2pp; mu1; gmu1(3); glnmu(3)
! PENCILS PROVIDED rho1gpp(3); glnpp(3); del2lnTT
!
! PENCILS PROVIDED hss(3,3); hlnTT(3,3); del2ss; del6ss; del6lnTT
! PENCILS PROVIDED yH; ee; ss; delta; glnmumol(3); ppvap; csvap2; cs2
! PENCILS PROVIDED cp1tilde; cp; gamma_m1; gamma
! PENCILS PROVIDED rho_anel; gradcp(3)
!
!
!***************************************************************
module EquationOfState
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'eos.h'
!
  integer, parameter :: ilnrho_ss=1,ilnrho_ee=2,ilnrho_pp=3
  integer, parameter :: ilnrho_lnTT=4,ilnrho_cs2=5
  integer, parameter :: irho_cs2=6, irho_ss=7, irho_lnTT=8, ilnrho_TT=9
  integer, parameter :: ipp_ss=11, irho_TT=10, ipp_cs2=12
  integer, parameter :: irho_eth=13, ilnrho_eth=14
!
  integer :: iglobal_cs2, iglobal_glnTT, ics
!
  real :: lnTT0=impossible
!
  real :: mu=1.
  real :: cs0=1., rho0=1.
  real :: cs20=1., lnrho0=0.
  logical :: lpp_as_aux=.false.
  real :: gamma=5./3.
  real :: Rgas_cgs=0., Rgas, Rgas_unit_sys=1.,  error_cp=1e-6
  real :: gamma_m1    !(=gamma-1)
  real :: gamma1   !(=1/gamma)
  real :: cp=impossible, cp1=impossible, cv=impossible, cv1=impossible
  real :: cs2bot=1., cs2top=1.
  integer :: ieosvars=-1, ieosvar1=-1, ieosvar2=-1, ieosvar_count=0
  integer :: ll1,ll2,mm1,mm2,nn1,nn2
  logical :: leos_isothermal=.false., leos_isentropic=.false.
  logical :: leos_isochoric=.false., leos_isobaric=.false.
  logical :: leos_localisothermal=.false.
  character (len=20) :: input_file
  logical, SAVE ::  lcheminp_eos=.false.
  logical :: l_gamma_m1=.false.
  logical :: l_gamma=.false.
  logical :: l_cp=.false.
  integer :: imass=1!, iTemp1=2,iTemp2=3,iTemp3=4
!
  real, dimension(nchemspec,18) :: species_constants
!
  real :: Cp_const=impossible
  real :: Pr_number=0.7
  logical :: lpres_grad = .false.
!
!NILS: Why do we spend a lot of memory allocating these variables here????
!MR: Is now allocated only once.
 real, dimension(mx,my,mz), target :: mu1_full
!
  namelist /eos_init_pars/ mu, cp, cs0, rho0, gamma, error_cp, lpp_as_aux
!
  namelist /eos_run_pars/  mu, cp, cs0, rho0, gamma, error_cp, lpp_as_aux
!
  contains
!***********************************************************************
    subroutine register_eos
!
!  14-jun-03/axel: adapted from register_eos
!
      use Sub, only: register_report_aux
!
      leos_chemistry=.true.
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
      use Mpicomm, only: stop_it
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
        call stop_it('unit_temperature is not found!')
      else
        Rgas=Rgas_unit_sys*unit_temperature/unit_velocity**2
      endif
!
      inquire(file='chem.inp',exist=cheminp)
      inquire(file='chem.in',exist=chemin)
      if(chemin .and. cheminp) call fatal_error('eos_chemistry',&
          'chem.inp and chem.in found. Please decide for one')
!
      if (cheminp) input_file='chem.inp'
      if (chemin) input_file='chem.in'
      lcheminp_eos = cheminp .or. chemin
!      inquire(FILE=input_file, EXIST=lcheminp_eos)
!
      if (lroot) then
!
       if (.not. lcheminp_eos ) then
        call fatal_error('initialize_eos',&
                        'chem.imp is not found!')
       else
        print*,'units_eos: chem.imp is found! Now cp, cv, gamma, mu are pencils ONLY!'
       endif
      endif
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

      if (.not.ldensity) then
        call put_shared_variable('rho0',rho0,caller='initialize_eos')
        call put_shared_variable('lnrho0',lnrho0)
      endif
      if (lchemistry) call put_shared_variable('mu1_full',mu1_full,caller='initialize_eos')
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
      if (ieosvar_count>=2) &
        call fatal_error("select_eos_variable", &
             "2 thermodynamic quantities have already been defined while attempting to add a 3rd: ") !//variable)
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
          if (lroot) print*,"select_eos_variable: Thermodynamic variable combination, ieosvar_selected= ",ieosvar_selected
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
!
    endsubroutine getmu
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
        case ('pp')
          if (ldensity_nolog .or. ltemperature_nolog) &
              call fatal_error('get_slices_eos',&
              'pp not implemented for ldensity_nolog .or. ltemperature_nolog')
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
         if (minval(p%TT)==0.) then
           call fatal_error('calc_pencils_eos','p%TT=0!')
         endif         
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
            p%rho1gpp(:,i) = p%pp/p%rho(:) &
               *(p%glnrho(:,i)+p%glnTT(:,i)+p%gmu1(:,i)/p%mu1(:))
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
! Laplasian of pressure
!
       if (lpenc_loc(i_del2pp)) then
         call dot2(p%glnTT,glnTT2)
         TT1del2TT=p%del2lnTT+glnTT2
         rho1del2rho=p%del2lnrho+p%glnrho2
         call dot2(p%gmu1,gmu12)
         call dot2(p%glnpp,glnpp2)
         call del2(mu1_full,del2mu1)
         p%del2pp&
             =p%pp*glnpp2&
             +p%pp*(rho1del2rho+TT1del2TT+del2mu1/p%mu1)&
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
    subroutine get_cp1(cp1_)
!
!  04-nov-06/axel: added to alleviate spurious use of pressure_gradient
!
!  return the value of cp1 to outside modules
!
      real, intent(out) :: cp1_
      call fatal_error('get_cp1','SHOULD NOT BE CALLED WITH eos_chemistry')
      cp1_=impossible
!
    endsubroutine get_cp1
!***********************************************************************
    subroutine get_cv1(cv1_)
!
!  22-dec-10/PJK: adapted from get_cp1
!
!  return the value of cv1 to outside modules
!
      real, intent(out) :: cv1_
      call fatal_error('get_cv1','SHOULD NOT BE CALLED WITH eos_chemistry')
      cv1_=impossible
!
    endsubroutine get_cv1
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
      call fatal_error('pressure_gradient_farray','SHOULD NOT BE CALLED WITH eos_chemistry')
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
      call fatal_error('pressure_gradient_farray','SHOULD NOT BE CALLED WITH eos_chemistry')
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
     call fatal_error('temperature_gradien','SHOULD NOT BE CALLED WITH eos_chemistry')
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
      call fatal_error('temperature_laplacian','SHould not be called!')
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
      call fatal_error('temperature_hessian', &
        'This routine is not coded for eos_chemistry')
!
      hlnTT(:,:,:)=0
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(hss)
      call keep_compiler_quiet(hlnrho)
!
    endsubroutine temperature_hessian
!***********************************************************************
    subroutine eosperturb(f,psize,ee,pp,ss)
!
!  Set f(l1:l2,m,n,iss), depending on the valyes of ee and pp
!  Adding pressure perturbations is not implemented
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(in), optional :: ee, pp, ss
      real, dimension(psize) :: lnrho_
!
      if (psize==nx) then
        lnrho_=f(l1:l2,m,n,ilnrho)
      elseif (psize==mx) then
        lnrho_=f(:,m,n,ilnrho)
      else
        call not_implemented("eosperturb")
      endif
!
     call fatal_error('eosperturb', &
        'This routine is not coded for eos_chemistry')
!
      call keep_compiler_quiet(present(ee))
      call keep_compiler_quiet(present(pp))
      call keep_compiler_quiet(present(ss))
!
    endsubroutine eosperturb
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
      call fatal_error('eoscalc_farray', &
          'This routine is not coded for eos_chemistry')
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
      call fatal_error('eoscalc_point', &
        'This routine is not coded for eos_chemistry')
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
      call fatal_error('eoscalc_pencil', &
        'This routine is not coded for eos_chemistry')
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
      call fatal_error('get_soundspeed', &
        'This routine is not coded for eos_chemistry')
!
      call keep_compiler_quiet(TT)
!
    endsubroutine get_soundspeed
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
      call fatal_error('isothermal_entropy', &
          'This routine is not coded for eos_chemistry')
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
      call fatal_error('isothermal_lnrho_ss', &
        'This routine is not coded for eos_chemistry')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(T0)
      call keep_compiler_quiet(rho0)
!
    endsubroutine isothermal_lnrho_ss
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call fatal_error('bc_ss_temp3_z', &
          'not implemented in eos_chemistry.f90')
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
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
      character (len=3) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_energy
!***********************************************************************
    subroutine bc_stellar_surface(f,topbot)
!
      real, dimension (:,:,:,:) :: f
      character (len=3) :: topbot
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
      character (len=3), intent (in) :: topbot
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
      character (len=3), intent (in) :: topbot
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
      character (len=3), intent (in) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hdss_z_iso
!***********************************************************************
    subroutine find_mass(element_name,MolMass)
!
!  Find mass of element
!
!  05-feb-08/nils: coded
!
      use Mpicomm, only: stop_it
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
        call stop_it('find_mass: Element not found!')
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
     !   call fatal_error('find_species_index',&
      !                 'no species has been found')
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
      use Mpicomm, only: stop_it
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
                  call stop_it("There were too many species, "//&
                      "please increase nchemspec!")
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
1000  if (emptyFile)  call stop_it('The input file chem.inp was empty!')
!
!  Check if nchemspec where not too large
!
      if (k<nchemspec-1) then
        print*,'nchemspec=',nchemspec
        call stop_it("There were too few species, "//&
            "please decrease nchemspec!")
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
            call find_species_index(specie_string,ind_glob,ind_chem,&
                found_specie)
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
                  read (unit=NumberOfElement_string,fmt='(I5)') &
                      NumberOfElement_i
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
              read (unit=ChemInpLine(1:75),fmt='(5E15.8)')  &
                  species_constants(ind_chem,iaa1(1):iaa1(5))
!
            elseif (ChemInpLine(80:80)=="3") then
              ! Read iaa1(6):iaa5(3)
              read (unit=ChemInpLine(1:75),fmt='(5E15.8)')  &
                  species_constants(ind_chem,iaa1(6):iaa2(3))
            elseif (ChemInpLine(80:80)=="4") then
              ! Read iaa2(4):iaa2(7)
              read (unit=ChemInpLine(1:75),fmt='(4E15.8)')  &
                  species_constants(ind_chem,iaa2(4):iaa2(7))
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
        write(file_id,'(F10.2,3F10.2)') species_constants(k,imass),&
            species_constants(k,iTemp1:iTemp3)
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
          write (143,*) 'specname[',trim(ispec),']=',"'",&
              trim(varname(ichemspec(k))),"'"
          write (143,*) 'specmass[',trim(ispec),']=',species_constants(k,imass)
        enddo
        close (143)
      endif
!
    endsubroutine write_thermodyn
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
      call fatal_error('get_stratz', 'Stratification for this EOS is not implemented. ')
!
      call keep_compiler_quiet(z)
      if (present(rho0z)) call keep_compiler_quiet(rho0z)
      if (present(dlnrho0dz)) call keep_compiler_quiet(dlnrho0dz)
      if (present(eth0z)) call keep_compiler_quiet(eth0z)
!
    endsubroutine get_stratz
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr

    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(cs20,p_par(1))

    endsubroutine pushpars2c
!***********************************************************************
endmodule EquationOfState
