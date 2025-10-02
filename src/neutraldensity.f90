! $Id$
!
!  This module is used both for the initial condition and during run time.
!  It contains dlnrhon_dt and init_lnrhon, among other auxiliary routines.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lneutraldensity = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED lnrhon; rhon; rhon1; glnrhon(3); grhon(3);
! PENCILS PROVIDED unglnrhon; ungrhon; del2rhon; glnrhon2
! PENCILS PROVIDED del6lnrhon; del6rhon; snglnrhon(3); alpha; zeta
!
!***************************************************************
module NeutralDensity
!
  use Cdata
  use Messages
!
  implicit none
!
  include 'neutraldensity.h'
!
  real :: kx_lnrhon=1.,ky_lnrhon=1.,kz_lnrhon=1.
  real :: ampllnrhon=0.,rhon_left=1.,rhon_right=1.
  real :: diffrhon=0.,diffrhon_hyper3=0.,diffrhon_shock=0.
  real :: lnrhon_const=0., rhon_const=1.
  real :: lnrhon_int=0.,lnrhon_ext=0.
  real :: lnrhon0,lnrhon_left,lnrhon_right,alpha,zeta
  real, dimension(3) :: diffrhon_hyper3_aniso=0.
  integer, parameter :: ndiff_max=4
  logical :: lmass_source=.false.,lcontinuity_neutral=.true.
  logical :: lupw_lnrhon=.false.,lupw_rhon=.false.
  logical :: ldiffn_normal=.false.,ldiffn_hyper3=.false.,ldiffn_shock=.false.
  logical :: ldiffn_hyper3lnrhon=.false.,ldiffn_hyper3_aniso=.false.
  logical :: lfreeze_lnrhonint=.false.,lfreeze_lnrhonext=.false.
  logical :: ldiffn_hyper3_polar=.false., luse_as_ionization=.false.
  logical :: lpretend_star,lramp_up
  real :: star_form_threshold=1.,star_form_exponent=1.5
  character (len=labellen), dimension(ninit) :: initlnrhon='nothing'
  character (len=labellen), dimension(ndiff_max) :: idiffn=''
  character (len=labellen) :: borderlnrhon='nothing', alpha_prescription='const'
  character (len=intlen) :: iinit_str
!
  namelist /neutraldensity_init_pars/ &
       ampllnrhon,initlnrhon,    &
       rhon_left,rhon_right,lnrhon_const,rhon_const, &
       idiffn,lneutraldensity_nolog,    &
       lcontinuity_neutral,lnrhon0,lnrhon_left,lnrhon_right, &
       alpha,zeta,kx_lnrhon,ky_lnrhon,kz_lnrhon,lpretend_star,&
       star_form_threshold,lramp_up,star_form_exponent, &
       alpha_prescription, luse_as_ionization
!
  namelist /neutraldensity_run_pars/ &
       diffrhon,diffrhon_hyper3,diffrhon_shock,   &
       lupw_lnrhon,lupw_rhon,idiffn,     &
       lnrhon_int,lnrhon_ext, &
       lfreeze_lnrhonint,lfreeze_lnrhonext,         &
       lnrhon_const,lcontinuity_neutral,borderlnrhon,    &
       diffrhon_hyper3_aniso,alpha,zeta,lpretend_star, &
       star_form_threshold,lramp_up,star_form_exponent, &
       alpha_prescription, luse_as_ionization
!
!  Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_rhonm=0,idiag_rhon2m=0,idiag_lnrhon2m=0
  integer :: idiag_rhonmin=0,idiag_rhonmax=0,idiag_unglnrhonm=0
  integer :: idiag_lnrhonmphi=0,idiag_rhonmphi=0,idiag_dtnd=0
  integer :: idiag_rhonmz=0, idiag_rhonmy=0, idiag_rhonmx=0
  integer :: idiag_rhonmxy=0, idiag_rhonmr=0
  integer :: idiag_neutralmass=0, idiag_alprec=0
!
! Auxiliaries
!
  real :: alpha_time
  real, dimension(nx) :: diffus_diffrhon

  integer :: enum_alpha_prescription = 0
  integer :: enum_borderlnrhon = 0
  contains
!***********************************************************************
    subroutine register_neutraldensity
!
!  Initialise variables which should know that we solve the
!  compressible hydro equations: ilnrhon; increase nvar accordingly.
!
!  28-feb-07/wlad: adapted from density
!
      use FArrayManager
!
      if (.not.lcartesian_coords) call not_implemented('register_neutraldensity','non-Cartesian coords')
!
      if (lneutraldensity_nolog) then
        call farray_register_pde('rhon',irhon)
        ilnrhon=irhon
      else
        call farray_register_pde('lnrhon',ilnrhon)
      endif
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_neutraldensity
!***********************************************************************
    subroutine initialize_neutraldensity
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  For compatibility with other applications, we keep the possibility
!  of giving diffrhon units of dxmin*cs0, but cs0 is not well defined general
!
!  28-feb-07/wlad: adapted
!
      use General, only: itoa
      use BorderProfiles, only: request_border_driving
      use FArrayManager
!
      integer :: i
      logical :: lnothing
!
!  Turn off continuity equation term for 0-D runs.
!
      if (dimensionality==0) then
        lcontinuity_neutral=.false.
        print*, 'initialize_neutraldensity: 0-D run, turned off continuity equation'
      endif
!
!  Initialize dust diffusion
!
      ldiffn_normal=.false.
      ldiffn_shock=.false.
      ldiffn_hyper3=.false.
      ldiffn_hyper3lnrhon=.false.
      ldiffn_hyper3_aniso=.false.
      ldiffn_hyper3_polar=.false.
!
      lnothing=.false.
!
      do i=1,ndiff_max
        select case (idiffn(i))
        case ('normal')
          if (lroot) print*,'diffusion: div(D*grad(rhon))'
          ldiffn_normal=.true.
        case ('hyper3')
          if (lroot) print*,'diffusion: (d^6/dx^6+d^6/dy^6+d^6/dz^6)rhon'
          ldiffn_hyper3=.true.
        case ('hyper3lnrhon')
          if (lneutraldensity_nolog) call fatal_error('initialize_neutraldensity', &
                             "'hyper3lnrhon' diffusion not meaningful for linear density")
          if (lroot) print*,'diffusion: (d^6/dx^6+d^6/dy^6+d^6/dz^6)lnrhon'
          ldiffn_hyper3lnrhon=.true.
       case ('hyper3_aniso')
          if (.not.lneutraldensity_nolog) &
            call not_implemented('initialize_neutraldensity',"anisotropic hyperdiffusion for log density")
          if (lroot) print*,'diffusion: (Dx*d^6/dx^6 + Dy*d^6/dy^6 + Dz*d^6/dz^6)rhon'
          ldiffn_hyper3_aniso=.true.
        case ('hyper3_cyl','hyper3-cyl','hyper3_sph','hyper3-sph')
          if (lroot) print*,'diffusion: Dhyper/pi^4 *(Delta(rhon))^6/Deltaq^2'
          ldiffn_hyper3_polar=.true.
        case ('shock')
          if (lroot) print*,'diffusion: shock diffusion'
          ldiffn_shock=.true.
          call not_implemented("initialize_neutraldensity", &
                               "shock diffusion assumes ions velocity for neutrals")
        case ('')
          if (lroot .and. (.not. lnothing)) print*,'diffusion: nothing'
        case default
          call fatal_error('initialize_neutraldensity', &
               'No such value for idiff('//trim(itoa(i))//'): '//trim(idiffn(i)))
        endselect
        lnothing=.true.
      enddo
!
!  If we're timestepping, die or warn if the the diffusion coefficient that
!  corresponds to the chosen diffusion type is not set.
!
      if (lrun) then
        if (ldiffn_normal.and.diffrhon==0.0) &
            call warning('initialize_neutraldensity','diffusion coefficient diffrhon is zero')
        if ( (ldiffn_hyper3 .or. ldiffn_hyper3lnrhon) .and. diffrhon_hyper3==0.0) &
          call fatal_error('initialize_neutraldensity','diffusion coefficient diffrhon_hyper3 is zero')
        if ( (ldiffn_hyper3_aniso) .and.  &
            ((diffrhon_hyper3_aniso(1)==0. .and. nxgrid/=1 ).or. &
             (diffrhon_hyper3_aniso(2)==0. .and. nygrid/=1 ).or. &
             (diffrhon_hyper3_aniso(3)==0. .and. nzgrid/=1 )) ) &
            call fatal_error('initialize_neutraldensity', &
                             'A diffusion coefficient of diffrhon_hyper3 is zero')
        if (ldiffn_shock .and. diffrhon_shock==0.0) &
          call fatal_error('initialize_neutraldensity','diffusion coefficient diffrhon_shock is zero')
      endif

      if (.not.lneutraldensity_nolog) then         ! for logarithmic density

        if (ldiffn_hyper3.or.ldiffn_hyper3_aniso) then
!
!  To obtain del6rhon with log density, one must work with a COMAUX slot in the f-array for linear density.
!
          if (lroot) print*,"initialize_neutraldensity: Creating comm. aux variable rhon for del6rhon ", &
                            "and setting exponentiation BC for rhon"
          if (irhon==0) call farray_register_auxiliary('rhon',irhon,communicated=.true.)

          if (.not.lperi(1)) then
            bcx12(irhon,:)='exp'
            fbcx(irhon,:) = float(ilnrhon)
          endif
          if (.not.lperi(2)) then
            bcy12(irhon,:)='exp'
            fbcy(irhon,:) = float(ilnrhon)
          endif
          if (.not.lperi(3)) then
            bcz12(irhon,:)='exp'
            fbcz(irhon,:) = float(ilnrhon)
          endif

        endif
      endif

      if (lpretend_star.and..not.lcylindrical_coords) &
        call not_implemented('initialize_neutraldensity', &
                             "lpretend_star for other than cylindrical coordinates")
!
      if (lfreeze_lnrhonint) lfreeze_varint(ilnrhon) = .true.
      if (lfreeze_lnrhonext) lfreeze_varext(ilnrhon) = .true.
!
!  Tell the BorderProfiles module if we intend to use border driving, so
!  that the module can request the right pencils.
!
      select case (borderlnrhon)
!
      case ('zero','0','constant','stratification')
         call request_border_driving(borderlnrhon)
      case ('nothing')
        if (lroot.and.ip<=5) print*,"initialize_neutraldensity: borderlnrhon='nothing'"
      case default
         call fatal_error('initialize_neutraldensity','no such borderlnrhon: '//trim(borderlnrhon))
      endselect
!
    endsubroutine initialize_neutraldensity
!***********************************************************************
    subroutine read_neutraldensity_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=neutraldensity_init_pars, IOSTAT=iostat)
!
    endsubroutine read_neutraldensity_init_pars
!***********************************************************************
    subroutine write_neutraldensity_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=neutraldensity_init_pars)
!
    endsubroutine write_neutraldensity_init_pars
!***********************************************************************
    subroutine read_neutraldensity_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=neutraldensity_run_pars, IOSTAT=iostat)
!
    endsubroutine read_neutraldensity_run_pars
!***********************************************************************
    subroutine write_neutraldensity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=neutraldensity_run_pars)
!
    endsubroutine write_neutraldensity_run_pars
!***********************************************************************
    subroutine init_lnrhon(f)
!
!  initialise lnrhon; called from start.f90
!
!  28-feb-07/wlad: adapted
!
      use General, only: itoa,complex_phase
      use Initcond
      use InitialCondition, only: initial_condition_lnrhon
      use General, only: notanumber
      use EquationOfState, only: cs20, cs2bot,cs2top
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      logical :: lnothing
      integer :: j
!
!  Set default values for sound speed at top and bottom.
!  These may be updated in one of the following initialization routines.
!
      cs2top=cs20; cs2bot=cs20
!
!  different initializations of lnrhon (called from start).
!  If initrhon does't match, f=0 is assumed (default).
!
      lnothing=.true.

      do j=1,ninit

         if (initlnrhon(j)/='nothing') then

            lnothing=.false.

            iinit_str=itoa(j)

            select case (initlnrhon(j))
!
! some one-liners from density
!
            case ('zero', '0'); f(:,:,:,ilnrhon)=0.
            case ('const_lnrhon'); f(:,:,:,ilnrhon)=lnrhon_const
            case ('const_rhon'); f(:,:,:,ilnrhon)=log(rhon_const)
            case ('constant'); f(:,:,:,ilnrhon)=log(rhon_left)
            case ('scale-ions')
              if (ldensity_nolog) then
                f(:,:,:,ilnrhon)=log(rhon_const)+log(f(:,:,:,ilnrho))
              else
                f(:,:,:,ilnrhon)=log(rhon_const)+f(:,:,:,ilnrho)
              endif
            case ('sinwave-z'); call sinwave(ampllnrhon,f,ilnrhon,kz=kz_lnrhon)
            case ('gaussian-noise')
               if (lnrhon_left /= 0.) f(:,:,:,ilnrhon)=lnrhon_left
               call gaunoise(ampllnrhon,f,ilnrhon,ilnrhon)
            case default
               call fatal_error('init_lnrhon','no such initlnrhon('// &
                                trim(iinit_str)//'): '//trim(initlnrhon(j)))
!
            endselect

            if (lroot) print*,'init_lnrhon: initlnrhon('//trim(iinit_str)//') = ',trim(initlnrhon(j))
         endif

      enddo
!
!  Interface for user's own initial conditon
!
      if (linitial_condition) call initial_condition_lnrhon(f)
!
      if (lnothing.and.lroot) print*,'init_lnrhon: nothing'
!
!  If unlogarithmic density considered, take exp of lnrhon resulting from
!  initialisation above. (irhon=ilnrhon!)
!
      if (lneutraldensity_nolog) f(:,:,:,irhon)=exp(f(:,:,:,ilnrhon))
!
!  sanity check
!
      if (notanumber(f(l1:l2,m1:m2,n1:n2,ilnrhon))) call error('init_lnrhon','imaginary density values')
!
    endsubroutine init_lnrhon
!***********************************************************************
    subroutine pencil_criteria_neutraldensity
!
!  All pencils that the NeutralDensity module depends on are specified here.
!
!  28-feb-07/wlad: adapted
!
!  always needed for ionization and recombination
!
      lpenc_requested(i_rho)  =.true.
      lpenc_requested(i_rhon) =.true.
      lpenc_requested(i_alpha)=.true.
      lpenc_requested(i_zeta)=.true.
!
!  needed for temperature prescription
!
      if (alpha_prescription=='Temp_dep') then
        lpenc_requested(i_TT)=.true.
      endif

      if (.not.lneutraldensity_nolog) then
        lpenc_requested(i_rho1) =.true.
        lpenc_requested(i_rhon1)=.true.
      endif
!
      if (lpretend_star) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_uu)=.true.
        lpenc_requested(i_rcyl_mn1)=.true.
      endif
!
      if (lcontinuity_neutral) then
        lpenc_requested(i_divun)=.true.
        if (lneutraldensity_nolog) then
          lpenc_requested(i_ungrhon)=.true.
        else
          lpenc_requested(i_unglnrhon)=.true.
        endif
      endif
      if (ldiffn_shock) then
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_gshock)=.true.
        if (lneutraldensity_nolog) then
          lpenc_requested(i_grhon)=.true.
          lpenc_requested(i_del2rhon)=.true.
        else
          lpenc_requested(i_glnrhon)=.true.
          lpenc_requested(i_glnrhon2)=.true.
        endif
      endif
      if (ldiffn_normal) then
        if (lneutraldensity_nolog) then
          lpenc_requested(i_del2rhon)=.true.
        else
          lpenc_requested(i_glnrhon2)=.true.
        endif
      endif
      if (ldiffn_hyper3.or.ldiffn_hyper3_aniso) lpenc_requested(i_del6rhon)=.true.
      if (ldiffn_hyper3.and..not.lneutraldensity_nolog) lpenc_requested(i_rhon)=.true.
      if (ldiffn_hyper3lnrhon) lpenc_requested(i_del6lnrhon)=.true.
      if (ldiffn_hyper3_polar.and..not.lneutraldensity_nolog) lpenc_requested(i_rhon1)=.true.
!
      if (lmass_source) lpenc_requested(i_rcyl_mn)=.true.
!
      lpenc_diagnos2d(i_lnrhon)=.true.
      lpenc_diagnos2d(i_rhon)=.true.
!
      if (idiag_rhonm/=0 .or. idiag_rhonmz/=0 .or. idiag_rhonmy/=0 .or. &
          idiag_rhonmx/=0 .or. idiag_rhon2m/=0 .or. idiag_rhonmin/=0 .or. &
          idiag_rhonmax/=0 .or. idiag_rhonmxy/=0 .or. idiag_neutralmass/=0) &
        lpenc_diagnos(i_rhon)=.true.
      if (idiag_lnrhon2m/=0) lpenc_diagnos(i_lnrhon)=.true.
      if (idiag_unglnrhonm/=0) lpenc_diagnos(i_unglnrhon)=.true.
      if (idiag_alprec/=0) lpenc_diagnos(i_alpha)=.true.
!
    endsubroutine pencil_criteria_neutraldensity
!***********************************************************************
    subroutine pencil_interdep_neutraldensity(lpencil_in)
!
!  Interdependency among pencils from the NeutralDensity module is
!    specified here.
!
!  28-feb-07/wlad: adapted
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lneutraldensity_nolog) then
        if (lpencil_in(i_rhon1)) lpencil_in(i_rhon)=.true.
      else
        if (lpencil_in(i_rhon)) lpencil_in(i_rhon1)=.true.
      endif
      if (lpencil_in(i_unglnrhon)) then
        lpencil_in(i_uun)=.true.
        lpencil_in(i_glnrhon)=.true.
      endif
      if (lpencil_in(i_ungrhon)) then
        lpencil_in(i_uun)=.true.
        lpencil_in(i_grhon)=.true.
      endif
      if (lpencil_in(i_glnrhon2)) lpencil_in(i_glnrhon)=.true.
      if (lpencil_in(i_snglnrhon)) then
        lpencil_in(i_snij)=.true.
        lpencil_in(i_glnrhon)=.true.
      endif
!  The pencils glnrhon and grhon come in a bundle.
      if (lpencil_in(i_glnrhon) .and. lpencil_in(i_grhon)) then
        if (lneutraldensity_nolog) then
          lpencil_in(i_grhon)=.false.
        else
          lpencil_in(i_glnrhon)=.false.
        endif
      endif
!
    endsubroutine pencil_interdep_neutraldensity
!***********************************************************************
    subroutine calc_pencils_neutraldensity(f,p)
!
!  Calculate NeutralDensity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  28-feb-07/wlad: adapted
!
      use Sub, only: grad,dot,dot2,u_dot_grad,del2,del6,multmv,g2ij
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      intent(inout) :: f,p
!
      real, dimension (nx) :: tmp,smooth_step_threshold
      real :: unit_alpha=impossible
      integer :: i
!
! lnrhon
      if (lpencil(i_lnrhon)) then
        if (lneutraldensity_nolog) then
          p%lnrhon=log(f(l1:l2,m,n,irhon))
        else
          p%lnrhon=f(l1:l2,m,n,ilnrhon)
        endif
      endif
! rhon1 and rhon
      if (lneutraldensity_nolog) then
        if (lpencil(i_rhon)) p%rhon=f(l1:l2,m,n,irhon)
        if (lpencil(i_rhon1)) p%rhon1=1.0/p%rhon
      else
        if (lpencil(i_rhon1)) p%rhon1=exp(-f(l1:l2,m,n,ilnrhon))
        if (lpencil(i_rhon)) p%rhon=1.0/p%rhon1
      endif
! glnrhon and grhon
      if (lpencil(i_glnrhon).or.lpencil(i_grhon)) then
        if (lneutraldensity_nolog) then
          call grad(f,ilnrhon,p%grhon)
          if (lpencil(i_glnrhon)) then
            do i=1,3
              p%glnrhon(:,i)=p%grhon(:,i)/p%rhon
            enddo
          endif
        else
          call grad(f,ilnrhon,p%glnrhon)
          if (lpencil(i_grhon)) then
            do i=1,3
              p%grhon(:,i)=p%rhon*p%glnrhon(:,i)
            enddo
          endif
        endif
      endif
! unglnrhon
      if (lpencil(i_unglnrhon)) then
        if (lneutraldensity_nolog) then
          call dot(p%uun,p%glnrhon,p%unglnrhon)
        else
          call u_dot_grad(f,ilnrhon,p%glnrhon,p%uun,p%unglnrhon,UPWIND=lupw_lnrhon)
        endif
      endif
! ungrhon
      if (lpencil(i_ungrhon)) then
        if (lneutraldensity_nolog) then
          call u_dot_grad(f,ilnrhon,p%grhon,p%uun,p%ungrhon,UPWIND=lupw_rhon)
        else
          call dot(p%uun,p%grhon,p%ungrhon)
        endif
      endif
! glnrhon2
      if (lpencil(i_glnrhon2)) call dot2(p%glnrhon,p%glnrhon2)
! del2rhon
      if (lpencil(i_del2rhon)) then
        call del2(f,ilnrhon,p%del2rhon)    ! can be either del2rhon or del2lnrhon
        if (.not.lneutraldensity_nolog) p%del2rhon=p%del2rhon+p%glnrhon2  ! here del2rhon/rhon
      endif
!
! del6rhon
! Note that irhon either points to the linear neutraldensity as a PDE variable
! for lneutraldensity_nolog, or as a COMAUX variable otherwise.
!
      if (lpencil(i_del6rhon).and.irhon>0) call del6(f,irhon,p%del6rhon)
! del6lnrhon
      if (lpencil(i_del6lnrhon)) call del6(f,ilnrhon,p%del6lnrhon)
! snglnrhon
      if (lpencil(i_snglnrhon)) call multmv(p%snij,p%glnrhon,p%snglnrhon)
!
! ionization and recombination pencils
!
      p%zeta=zeta
      if (lpretend_star) then
!
! Star formation rate
! These lines below recover d/dt(rho_star)=sfr_const*omega*rho_gas**1.5
!
! There is a threshold that has to be smoothed. The star formation
! rate falls drastically after the threshold of 5-10 solar masses
! per cubic parsec. An arctan smoothing over a tenth of this value
! is okay to avoid numerical disasters.
!
        tmp=(p%rho-star_form_threshold)/(.1*star_form_threshold)
        smooth_step_threshold=.5*(1+atan(tmp)*2*pi_1)

        !OO=p%uu(*,2)*p%rcyl_mn1
        p%alpha=(alpha_time*smooth_step_threshold)*(p%uu(:,2)*p%rcyl_mn1)*p%rho1**(2-star_form_exponent)
!
      else
!
!  Choice of different alpha prescriptions
!
        select case (alpha_prescription)
        case ('const'); p%alpha=alpha
        case ('Temp_dep')
          if (lentropy) then
            unit_alpha=1./(unit_time*unit_density)
            p%alpha=unit_alpha*1.14e-19*4.309e6*p%TT**(-.6166)/(1.+.6703*p%TT**.53)
  !         if (ip<10) print*,'AXEL: p%alpha=',p%alpha
  !         if (ip<10) print*,'AXEL: unit_alpha=',unit_alpha
          else
            call fatal_error('calc_pencils_neutraldensity', &
              'no energy equation is used')
          endif
        case default
          call fatal_error('calc_pencils_neutraldensity', &
            'No such value for alpha_prescription')
        endselect
      endif
!
    endsubroutine calc_pencils_neutraldensity
!***********************************************************************
    subroutine neutraldensity_after_boundary(f)

      use General, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray) :: f
      real :: ramping_period

      if (lpretend_star) then
!
! Smooth transition over a ramping period of 5 orbits.
!
        if (lramp_up) then              ! MR: strange use of lpoint!
          ramping_period=2*pi*x(lpoint) !omega=v/r; v=1; 1/omega=r
          if (t <= ramping_period) then
            alpha_time=alpha*(sin((.5*pi)*(t/ramping_period))**2)
          else
            alpha_time=alpha
          endif
        else
          alpha_time=alpha
        endif

      endif

      call keep_compiler_quiet(f)

    endsubroutine neutraldensity_after_boundary
!***********************************************************************
    subroutine neutraldensity_before_boundary(f)

      use General, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray) :: f
!
! Fill global rhon array using the ilnrhon data.
!
      if (.not.lneutraldensity_nolog.and.irhon/=0) &
        f(l1:l2,m1:m2,n1:n2,irhon) = exp(f(l1:l2,m1:m2,n1:n2,ilnrhon))

    endsubroutine neutraldensity_before_boundary
!***********************************************************************
    subroutine dlnrhon_dt(f,df,p)
!
!  continuity equation
!  calculate dlnrhon/dt = - un.gradlnrhon - divun
!
!  28-feb-07/wlad: adapted
!
      use Deriv, only: der6
      use Sub, only: identify_bcs, del6fj, dot_mn
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)  :: f,p
      intent(out) :: df
!
      real, dimension(nx) :: fdiff,gshockglnrhon,gshockgrhon,tmp,diffus_diffrhon3
      integer :: j
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dlnrhon_dt: SOLVE dlnrhon_dt'
      if (headtt) then
        call identify_bcs('lnrhon',ilnrhon)
        if (.not.lneutraldensity_nolog.and.irhon/=0) call identify_bcs('rhon',irhon)
      endif
!
!  continuity equation
!
      if (lcontinuity_neutral) then
        if (lneutraldensity_nolog) then
          df(l1:l2,m,n,irhon)   = df(l1:l2,m,n,irhon)   - p%ungrhon   - p%rhon*p%divun
        else
          df(l1:l2,m,n,ilnrhon) = df(l1:l2,m,n,ilnrhon) - p%unglnrhon - p%divun
        endif
      endif
!
!  Ionization and recombination
!  If the helium fraction Y_He=0.079 is to be included, we can redefine alpha correspondingly.
!
      if (lneutraldensity_nolog) then
         df(l1:l2,m,n,irhon)   = df(l1:l2,m,n,irhon)   - p%zeta*p%rhon        + p%alpha*p%rho**2
         df(l1:l2,m,n,ilnrho ) = df(l1:l2,m,n,ilnrho ) + p%zeta*p%rhon        - p%alpha*p%rho**2
      else
         if (luse_as_ionization) then
           df(l1:l2,m,n,ilnrhon) = df(l1:l2,m,n,ilnrhon) - p%alpha*p%rhon1
         else
           df(l1:l2,m,n,ilnrhon) = df(l1:l2,m,n,ilnrhon) - p%zeta               + p%alpha*p%rho**2*p%rhon1
           df(l1:l2,m,n,ilnrho ) = df(l1:l2,m,n,ilnrho ) + p%zeta*p%rhon*p%rho1 - p%alpha*p%rho
         endif
      endif
!
!  Normal mass diffusion
!
      if (dimensionality>0) then

        fdiff=0.0
        diffus_diffrhon=0.; diffus_diffrhon3=0.
!
        if (ldiffn_normal) then  ! Normal diffusion operator
          fdiff = fdiff + diffrhon*p%del2rhon
          if (lupdate_courant_dt) diffus_diffrhon=diffus_diffrhon+diffrhon*dxyz_2
          if (headtt) print*,'dlnrhon_dt: diffrhon=', diffrhon
        endif
!
!  Hyper diffusion
!
        if (ldiffn_hyper3) then
          if (lneutraldensity_nolog) then
            fdiff = fdiff + diffrhon_hyper3*p%del6rhon
          else
            fdiff = fdiff + 1/p%rhon*diffrhon_hyper3*p%del6rhon
          endif
          if (lupdate_courant_dt) diffus_diffrhon3=diffus_diffrhon3+diffrhon_hyper3*dxyz_6
          if (headtt) print*,'dlnrhon_dt: diffrhon_hyper3=', diffrhon_hyper3
        endif
!
        if (ldiffn_hyper3_aniso) then
          if (lneutraldensity_nolog) then
            call del6fj(f,diffrhon_hyper3_aniso,ilnrhon,tmp)
            fdiff = fdiff + tmp
            if (lupdate_courant_dt) diffus_diffrhon3=diffus_diffrhon3 + &
                                                 diffrhon_hyper3_aniso(1)*dline_1(:,1)**6 + &
                                                 diffrhon_hyper3_aniso(2)*dline_1(:,2)**6 + &
                                                 diffrhon_hyper3_aniso(3)*dline_1(:,3)**6
            if (headtt) print*,'dlnrhon_dt: diffrhon_hyper3=(Dx,Dy,Dz)=',diffrhon_hyper3_aniso
          endif
        endif
!
        if (ldiffn_hyper3_polar) then
          do j=1,3
            call der6(f,ilnrhon,tmp,j,IGNOREDX=.true.)
            if (.not.lneutraldensity_nolog) tmp=tmp*p%rhon1
            fdiff = fdiff + diffrhon_hyper3*pi4_1*tmp*dline_1(:,j)**2
          enddo
          if (lupdate_courant_dt) diffus_diffrhon3=diffus_diffrhon3+diffrhon_hyper3*pi4_1*dxmin_pencil**4
          if (headtt) print*,'dlnrhon_dt: diffrhon_hyper3=', diffrhon_hyper3
        endif
!
        if (ldiffn_hyper3lnrhon) then
          if (.not. lneutraldensity_nolog) fdiff = fdiff + diffrhon_hyper3*p%del6lnrhon
          if (lupdate_courant_dt) diffus_diffrhon3=diffus_diffrhon3+diffrhon_hyper3*dxyz_6
          if (headtt) print*,'dlnrhon_dt: diffrhon_hyper3=', diffrhon_hyper3
        endif
!
!  Shock diffusion
!
        if (ldiffn_shock) then
          if (lneutraldensity_nolog) then
            call dot_mn(p%gshock,p%grhon,gshockgrhon)
            df(l1:l2,m,n,ilnrhon) = df(l1:l2,m,n,ilnrhon) + &
                diffrhon_shock*p%shock*p%del2rhon + diffrhon_shock*gshockgrhon
          else
            call dot_mn(p%gshock,p%glnrhon,gshockglnrhon)
            df(l1:l2,m,n,ilnrhon) = df(l1:l2,m,n,ilnrhon) + &
                diffrhon_shock*p%shock*p%del2rhon + diffrhon_shock*gshockglnrhon
          endif
          if (lupdate_courant_dt) diffus_diffrhon=diffus_diffrhon+diffrhon_shock*p%shock*dxyz_2
          if (headtt) print*,'dlnrhon_dt: diffrhon_shock=', diffrhon_shock
        endif
!
!  Add diffusion term to continuity equation
!
        if (lneutraldensity_nolog) then
          df(l1:l2,m,n,irhon)   = df(l1:l2,m,n,irhon)   + fdiff
        else
          df(l1:l2,m,n,ilnrhon) = df(l1:l2,m,n,ilnrhon) + fdiff
        endif
!
        if (lupdate_courant_dt) then
          if (headtt.or.ldebug) print*,'dlnrhon_dt: max(diffus_diffrhon) =',maxval(diffus_diffrhon)
          maxdiffus=max(maxdiffus,diffus_diffrhon)
          maxdiffus3=max(maxdiffus3,diffus_diffrhon3)
        endif
      endif
!
!  Apply border profile
!
      if (lborder_profiles) call set_border_neutraldensity(f,df,p)
!
      call calc_diagnostics_neutraldens(p)

    endsubroutine dlnrhon_dt
!***********************************************************************
    subroutine calc_diagnostics_neutraldens(p)
! 
      use Diagnostics

      type (pencil_case) :: p
!
!  2d-averages
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      if (l2davgfirst) then
        call zsum_mn_name_xy(p%rhon,idiag_rhonmxy)
        call phisum_mn_name_rz(p%lnrhon,idiag_lnrhonmphi)
        call phisum_mn_name_rz(p%rhon,idiag_rhonmphi)
      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1
!
      if (l1davgfirst) then
         call phizsum_mn_name_r(p%rhon,idiag_rhonmr)
         call xysum_mn_name_z(p%rhon,idiag_rhonmz)
         call yzsum_mn_name_x(p%rhon,idiag_rhonmx)
         call xzsum_mn_name_y(p%rhon,idiag_rhonmy)
      endif
!
!  Calculate density diagnostics
!
      if (ldiagnos) then
        call sum_mn_name(p%rhon,idiag_rhonm)
        call sum_mn_name(p%alpha,idiag_alprec)
        call sum_lim_mn_name(p%rhon,idiag_neutralmass,p)
        if (idiag_rhonmin/=0)  call max_mn_name(-p%rhon,idiag_rhonmin,lneg=.true.)
        call max_mn_name(p%rhon,idiag_rhonmax)
        if (idiag_rhon2m/=0)   call sum_mn_name(p%rhon**2,idiag_rhon2m)
        if (idiag_lnrhon2m/=0) call sum_mn_name(p%lnrhon**2,idiag_lnrhon2m)
        call sum_mn_name(p%unglnrhon,idiag_unglnrhonm)
        if (dimensionality>0) then
          if (idiag_dtnd/=0)   call max_mn_name(diffus_diffrhon/cdtv,idiag_dtnd,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_diagnostics_neutraldens
!***********************************************************************
    subroutine set_border_neutraldensity(f,df,p)
!
!  Calculates the driving term for the border profile
!  of the lnrhon variable.
!
!  28-jul-06/wlad: coded
!
      use BorderProfiles, only: border_driving
      use EquationOfState, only: cs0,cs20
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension(mx,my,mz,mvar) :: df

      real, dimension(nx) :: f_target !,OO_sph,OO_cyl,cs,theta
!      real :: r0_pot=0.1
!
      select case (borderlnrhon)
!
      case ('zero','0')
         f_target=0.
         if (lneutraldensity_nolog) f_target=exp(f_target)
         call border_driving(f,df,p,f_target,ilnrhon)
      case ('constant')
         f_target=lnrhon_const
         if (lneutraldensity_nolog) f_target=exp(f_target)
         call border_driving(f,df,p,f_target,ilnrhon)
      case ('stratification')
         !OO_sph = sqrt((r_mn**2 + r0_pot**2)**(-1.5))
         !OO_cyl = sqrt((rcyl_mn**2 + r0_pot**2)**(-1.5))
         !cs = OO_cyl*rcyl_mn*cs0
         !f_target=lnrhon_const - 0.5*(theta/cs0)**2
         f_target=(p%rcyl_mn-p%r_mn)/(cs20*p%r_mn)
         if (lneutraldensity_nolog) f_target=exp(f_target)
         call border_driving(f,df,p,f_target,ilnrhon)
      case ('nothing')
      endselect
!
!
    endsubroutine set_border_neutraldensity
!***********************************************************************
    subroutine rprint_neutraldensity(lreset,lwrite)
!
!  reads and registers print parameters relevant for compressible part
!
!  28-feb-07/wlad: adapted
!
      use Diagnostics, only: parse_name
      use FArrayManager, only: farray_index_append
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamex, inamey, inamez, inamexy, irz, inamer
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_rhonm=0; idiag_rhon2m=0; idiag_lnrhon2m=0; idiag_unglnrhonm=0
        idiag_rhonmin=0; idiag_rhonmax=0; idiag_dtnd=0
        idiag_lnrhonmphi=0; idiag_rhonmphi=0
        idiag_rhonmz=0; idiag_rhonmy=0; idiag_rhonmx=0
        idiag_rhonmxy=0; idiag_rhonmr=0; idiag_neutralmass=0
        diffrhon=0.; idiag_alprec=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_neutraldensity: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'rhonm',idiag_rhonm)
        call parse_name(iname,cname(iname),cform(iname),'rhon2m',idiag_rhon2m)
        call parse_name(iname,cname(iname),cform(iname),'rhonmin',idiag_rhonmin)
        call parse_name(iname,cname(iname),cform(iname),'rhonmax',idiag_rhonmax)
        call parse_name(iname,cname(iname),cform(iname),'lnrhon2m',idiag_lnrhon2m)
        call parse_name(iname,cname(iname),cform(iname),'unglnrhonm',idiag_unglnrhonm)
        call parse_name(iname,cname(iname),cform(iname),'dtnd',idiag_dtnd)
        call parse_name(iname,cname(iname),cform(iname),'neutralmass',idiag_neutralmass)
        call parse_name(iname,cname(iname),cform(iname),'alprec',idiag_alprec)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhonmz',idiag_rhonmz)
      enddo
!
!  check for those quantities for which we want xz-averages
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'rhonmy',idiag_rhonmy)
      enddo
!
!  check for those quantities for which we want yz-averages
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'rhonmx',idiag_rhonmx)
      enddo
!
!  check for those quantities for which we want phiz-averages
!
      do inamer=1,nnamer
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'rhonmr',idiag_rhonmr)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'rhonmxy',idiag_rhonmxy)
      enddo
!
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'lnrhonmphi',idiag_lnrhonmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'rhonmphi',idiag_rhonmphi)
      enddo
!
!  write column where which density variable is stored
!
      if (lwr) then
        call farray_index_append('ilnrhon',ilnrhon)
      endif
!
    endsubroutine rprint_neutraldensity
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr
    use General , only: string_to_enum

    integer, parameter :: n_pars=30
    integer(KIND=ikind8), dimension(n_pars) :: p_par


    call copy_addr(diffrhon,p_par(1))
    call copy_addr(diffrhon_hyper3,p_par(2))
    call copy_addr(diffrhon_shock,p_par(3))
    call copy_addr(lnrhon_const,p_par(4))
    call copy_addr(rhon_const,p_par(5))
    call copy_addr(alpha,p_par(6))
    call copy_addr(zeta,p_par(7))
    call copy_addr(lcontinuity_neutral,p_par(8)) ! bool
    call copy_addr(lupw_lnrhon,p_par(9)) ! bool
    call copy_addr(lupw_rhon,p_par(10)) ! bool
    call copy_addr(ldiffn_normal,p_par(11)) ! bool
    call copy_addr(ldiffn_hyper3,p_par(12)) ! bool
    call copy_addr(ldiffn_shock,p_par(13)) ! bool
    call copy_addr(ldiffn_hyper3lnrhon,p_par(14)) ! bool
    call copy_addr(ldiffn_hyper3_aniso,p_par(15)) ! bool
    call copy_addr(ldiffn_hyper3_polar,p_par(16)) ! bool
    call copy_addr(luse_as_ionization,p_par(17)) ! bool
    call copy_addr(lpretend_star,p_par(18)) ! bool
    call copy_addr(star_form_threshold,p_par(19))
    call copy_addr(star_form_exponent,p_par(20))
    call copy_addr(alpha_time,p_par(21))
    call copy_addr(diffrhon_hyper3_aniso,p_par(22)) ! real3
    call string_to_enum(enum_alpha_prescription,alpha_prescription)
    call copy_addr(enum_alpha_prescription,p_par(23)) ! int
    call string_to_enum(enum_borderlnrhon,borderlnrhon)
    call copy_addr(enum_borderlnrhon,p_par(24)) ! int
    endsubroutine pushpars2c
!***********************************************************************
endmodule NeutralDensity
