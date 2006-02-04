! $Id: dustvelocity_onefluid.f90,v 1.5 2006-02-04 10:10:31 ajohan Exp $
!
!  This module takes care of everything related to dust velocity
!  in the one-fluid limit (valid for vanishing friction time).
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED divud,ood,od2,oud,ud2,udij,sdij,udgud,uud
! PENCILS PROVIDED del2ud,del6ud,graddivud
!
!***************************************************************

module Dustvelocity

!  Note that Omega is already defined in cdata.

  use Cdata
  use Messages
  use Hydro

  implicit none

  include 'dustvelocity.h'

  public :: dust_geometry, dimd1, rhods, surfd, mdplus, mdminus
  public :: ad, scolld, ustcst, tausd1, tausd
  public :: unit_md, dust_chemistry, mumon, mmon, mi, md
  character (len=labellen) :: dust_geometry='sphere', dust_chemistry='nothing'
  real, dimension(ndustspec,ndustspec) :: scolld
  real, dimension(nx,ndustspec) :: tausd1
  real, dimension(ndustspec) :: md=1.0, mdplus=0.0, mdminus=0.0, surfd=0.0
  real, dimension(ndustspec) :: mi=0.0, ad=1.0, tausd=0.0
  real :: dimd1=0.0, rhods=0.0, ustcst=0.0, unit_md=0.0
  real :: mumon=0.0, mmon=0.0

  ! init parameters
  real :: nud=0.0
  real :: kx_uud=1.0, ky_uud=1.0, kz_uud=1.0
  real :: ampluud=0.0, ampl_udx=0.0, ampl_udy=0.0, ampl_udz=0.0
  real :: phase_udx=0.0, phase_udy=0.0, phase_udz=0.0
  real :: nd0=0.0, rhod0=0.0
  logical :: ldustcoagulation=.false.,ldustcondensation=.false.
  logical :: ladvection_dust=.true.,lcoriolisforce_dust=.true.
  logical :: lpressure_gradient=.true.
  logical :: lviscosity_dust=.true.,lvshear_dust_global_eps=.false.
  character (len=labellen), dimension(ninit) :: inituud='nothing'
  character (len=labellen) :: iviscd='simplified'

  namelist /dustvelocity_init_pars/ &
      ampl_udx, ampl_udy, ampl_udz, phase_udx, phase_udy, phase_udz, &
      ampluud, inituud, Omega, &
      kx_uud, ky_uud, kz_uud, lvshear_dust_global_eps

  ! run parameters
  namelist /dustvelocity_run_pars/ &
       nud, iviscd, ladvection_dust, lcoriolisforce_dust, &
       lpressure_gradient

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_ud2m=0
  integer :: idiag_udxm=0,idiag_udym=0,idiag_udzm=0
  integer :: idiag_udx2m=0,idiag_udy2m=0,idiag_udz2m=0
  integer :: idiag_udm2=0,idiag_oudm=0,idiag_od2m=0
  integer :: idiag_udxpt=0,idiag_udypt=0,idiag_udzpt=0
  integer :: idiag_udrms=0,idiag_udmax=0,idiag_odrms=0
  integer :: idiag_odmax=0,idiag_rdudmax=0
  integer :: idiag_udxmz=0,idiag_udymz=0,idiag_udzmz=0
  integer :: idiag_udx2mz=0,idiag_udy2mz=0,idiag_udz2mz=0
  integer :: idiag_udmx=0,idiag_udmy=0,idiag_udmz=0
  integer :: idiag_udxmxy=0,idiag_udymxy=0,idiag_udzmxy=0
  integer :: idiag_divud2m=0,idiag_epsKd=0
  integer :: idiag_dtud=0,idiag_dtnud=0
  integer :: idiag_rdudxm=0,idiag_rdudym=0,idiag_rdudzm=0
  integer :: idiag_rdudx2m=0


  contains

!***********************************************************************
    subroutine register_dustvelocity()
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  18-mar-03/axel+anders: adapted from hydro
!
      use Cdata
      use Sub
      use General, only: chn
!
      logical, save :: first=.true.
      integer :: k
      character(len=4) :: sdust
!
      if (.not. first) call fatal_error('register_dustvelocity','module registration called twice')
      first = .false.
!
      ldustvelocity = .true.
!
      iuud = nvar+1
      iudx = nvar+1
      iudy = nvar+2
      iudz = nvar+3
      nvar = nvar+3
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_dustvelocity: nvar = ', nvar
        print*, 'register_dustvelocity: iudx,iudy,iudz = ', &
            iudx,iudy,iudz
      endif
!
!  Put variable name in array
!
      varname(iudx) = 'udx'
      varname(iudy) = 'udy'
      varname(iudz) = 'udz'
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: dustvelocity_onefluid.f90,v 1.5 2006-02-04 10:10:31 ajohan Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call fatal_error('register_dustvelocity','nvar > mvar')
      endif
!
    endsubroutine register_dustvelocity
!***********************************************************************
    subroutine initialize_dustvelocity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  02-feb-06/anders: adapted
!
!  Turn off dust viscosity if zero viscosity
!
      if (nud==0.0) lviscosity_dust=.false.
      if (lroot) print*, &
          'initialize_dustvelocity: lviscosity_dust=', lviscosity_dust
!
    endsubroutine initialize_dustvelocity
!***********************************************************************
    subroutine copy_bcs_dust
!
!  Copy boundary conditions on first dust species to all others
!    
!  02-feb-06/anders: dummy
!
    endsubroutine copy_bcs_dust
!***********************************************************************
    subroutine init_uud(f)
!
!  initialise uud; called from start.f90
!
!  02-feb-06/anders: adapted
!
      use Cdata
      use EquationOfState, only: gamma, beta_glnrho_global, beta_glnrho_scaled
      use Sub
      use Global
      use Gravity
      use Initcond
      use EquationOfState, only: pressure_gradient,cs20
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: lnrho,rho,cs2,rhod,cp1tilde
      real :: eps,cs,eta_glnrho,v_Kepler
      integer :: j,l
      logical :: lnothing
!
!  inituud corresponds to different initializations of uud (called from start).
!
      lnothing=.false.
      do j=1,ninit
        select case(inituud(j))

        case('nothing')
          if (lroot .and. .not. lnothing) print*, 'init_uud: nothing'
          lnothing=.true.
        case('zero', '0')
          f(:,:,:,iudx(1):iudz(1))=0.0
          if (lroot) print*,'init_uud: zero dust velocity'
        case('gaussian-noise')
          call gaunoise(ampluud,f,iudx(1),iudz(1))
        case('sinwave-phase')
          call sinwave_phase(f,iudx(1),ampl_udx,kx_uud,ky_uud,kz_uud,phase_udx)
          call sinwave_phase(f,iudy(1),ampl_udy,kx_uud,ky_uud,kz_uud,phase_udy)
          call sinwave_phase(f,iudz(1),ampl_udz,kx_uud,ky_uud,kz_uud,phase_udz)
        case('udx_sinx')
          do l=1,mx; f(l,:,:,iudx(1)) = ampluud*sin(kx_uud*x(l)); enddo
        case('udy_siny')
          do m=1,my; f(:,m,:,iudy(1)) = ampluud*sin(ky_uud*y(m)); enddo
        case('sinwave-z-x')
          if (lroot) print*, 'init_uud: sinwave-z-x, ampluud=', ampluud
          call sinwave(ampluud,f,iudz(1),kx=kx_uud)
        case('udz_sinz')
          do n=1,mz; f(:,:,n,iudz(1)) = ampluud*sin(kz_uud*z(n)); enddo
        case('udz_siny')
          do m=m1,m2
            f(:,m,:,iudz(1)) = f(:,m,:,iudz(1)) + ampluud*sin(ky_uud*y(m))
          enddo
        case('udx_sinxsinysinz')
          do l=1,mx; do m=1,my; do n=1,mz
            f(l,m,n,iudx(1)) = &
                ampluud*sin(kx_uud*x(l))*sin(ky_uud*y(m))*sin(kz_uud*z(n))
          enddo; enddo; enddo
        case('udy_sinxsinysinz')
          do l=1,mx; do m=1,my; do n=1,mz
            f(l,m,n,iudy(1)) = &
                ampluud*sin(kx_uud*x(l))*sin(ky_uud*y(m))*sin(kz_uud*z(n))
          enddo; enddo; enddo
        case('udz_sinxsinysinz')
          do l=1,mx; do m=1,my; do n=1,mz
            f(l,m,n,iudz(1)) = &
                ampluud*sin(kx_uud*x(l))*sin(ky_uud*y(m))*sin(kz_uud*z(n))
          enddo; enddo; enddo

        case('dragforce_equilibrium')
!
!  Equilibrium between drag force and other forces.
!
          if (lroot) then
            print*, 'init_uud: drag force equilibrium'
            print*, 'init_uud: beta_glnrho_scaled=', beta_glnrho_scaled
          endif

          if (ldensity_nolog) then
            if (ldustdensity_log) then
              eps=sum(exp(f(l1:l2,m1:m2,n1:n2,ind(1))))/sum(f(l1:l2,m1:m2,n1:n2,ilnrho))
            else
              eps=sum(f(l1:l2,m1:m2,n1:n2,ind(1)))/sum(f(l1:l2,m1:m2,n1:n2,ilnrho))
            endif
          else
            if (ldustdensity_log) then
              eps=sum(exp(f(l1:l2,m1:m2,n1:n2,ind(1))))/sum(exp(f(l1:l2,m1:m2,n1:n2,ilnrho)))
            else
              eps=sum(f(l1:l2,m1:m2,n1:n2,ind(1)))/sum(exp(f(l1:l2,m1:m2,n1:n2,ilnrho)))
            endif
          endif

          if (lroot) print*, 'init_uud: average dust-to-gas ratio=', eps
          
          do l=l1,l2; do m=m1,m2; do n=n1,n2
            cs=sqrt(cs20)

            if (.not. lvshear_dust_global_eps) then
              if (ldensity_nolog) then
                if (ldustdensity_log) then
                  eps=exp(f(l,m,n,ind(1)))/f(l,m,n,ilnrho)
                else
                  eps=f(l,m,n,ind(1))/f(l,m,n,ilnrho)
                endif
              else
                if (ldustdensity_log) then
                  eps=exp(f(l,m,n,ind(1)))/exp(f(l,m,n,ilnrho))
                else
                  eps=f(l,m,n,ind(1))/exp(f(l,m,n,ilnrho))
                endif
              endif
            endif

            if (beta_glnrho_scaled(1)/=0.0) then
              f(l,m,n,iudx(1)) = 0.0
              f(l,m,n,iudy(1)) = f(l,m,n,iudy(1)) + &
                  1/gamma*cs20*beta_glnrho_scaled(1)*(1+eps)/ &
                  (2*Omega*(1.0+2*eps+eps**2))
            endif
          enddo; enddo; enddo
!
!  Catch unknown values
!
        case default
          write (unit=errormsg,fmt=*) 'No such such value for inituu: ', trim(inituud(j))
          call fatal_error('init_uud',errormsg)

        endselect
!
!  End loop over initial conditions
!        
      enddo
!
    endsubroutine init_uud
!***********************************************************************
    subroutine pencil_criteria_dustvelocity()
! 
!  All pencils that the Dustvelocity module depends on are specified here.
! 
!  20-11-04/anders: coded
!
      integer :: i
!
      lpenc_requested(i_uud)=.true.
      if (ladvection_dust) lpenc_requested(i_udgud)=.true.
      if (lpressure_gradient) then
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_epsd)=.true.
        if (lentropy) lpenc_requested(i_gss)=.true.
        if (lentropy) lpenc_requested(i_cp1tilde)=.true.
      endif
      if (lviscosity_dust) then
        if ((iviscd=='nud-const' .or. iviscd=='hyper3_nud-const') &
            .and. ldustdensity) then
          lpenc_requested(i_sdij)=.true.
          lpenc_requested(i_glnnd)=.true.
        endif
        if (iviscd=='simplified' .or. iviscd=='nud-const') &
            lpenc_requested(i_del2ud)=.true.
        if (iviscd=='hyper3_simplified' .or. iviscd=='hyper3_nud-const' .or. &
            iviscd=='hyper3_rhod_nud-const') &
            lpenc_requested(i_del6ud)=.true.
        if (iviscd=='nud-const' .or. iviscd=='hyper3_nud-const') &
            lpenc_requested(i_sdglnnd)=.true.
        if (iviscd=='nud-const') lpenc_requested(i_graddivud)=.true.
        if (iviscd=='hyper3_rhod_nud-const') lpenc_requested(i_rhod)=.true.
      endif
!
      lpenc_diagnos(i_uud)=.true.
      if (idiag_divud2m/=0) lpenc_diagnos(i_divud)=.true.
      if (idiag_rdudmax/=0 .or. idiag_rdudxm/=0 .or. &
          idiag_rdudym/=0 .or. idiag_rdudzm/=0 .or. &
          idiag_rdudx2m/=0) &
          lpenc_diagnos(i_rhod)=.true.
      if (idiag_udrms/=0 .or. idiag_udmax/=0 .or. &
          idiag_rdudmax/=0 .or. idiag_ud2m/=0 .or. &
          idiag_udm2/=0) &
          lpenc_diagnos(i_ud2)=.true.
      if (idiag_odrms/=0 .or. idiag_odmax/=0 .or. &
          idiag_od2m/=0) lpenc_diagnos(i_od2)=.true.
      if (idiag_oudm/=0) lpenc_diagnos(i_oud)=.true.
!
    endsubroutine pencil_criteria_dustvelocity
!***********************************************************************
    subroutine pencil_interdep_dustvelocity(lpencil_in)
!
!  Interdependency among pencils provided by the Dustvelocity module
!  is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_ud2)) lpencil_in(i_uud)=.true.
      if (lpencil_in(i_divud)) lpencil_in(i_udij)=.true.
      if (lpencil_in(i_udgud)) then
        lpencil_in(i_uud)=.true.
        lpencil_in(i_udij)=.true.
      endif
      if (lpencil_in(i_ood)) lpencil_in(i_udij)=.true.
      if (lpencil_in(i_od2)) lpencil_in(i_ood)=.true.
      if (lpencil_in(i_oud)) then
        lpencil_in(i_uud)=.true.
        lpencil_in(i_ood)=.true.
      endif
      if (lpencil_in(i_sdij)) then
        if (iviscd=='nud-const') then
          lpencil_in(i_udij)=.true.
          lpencil_in(i_divud)=.true.
        endif
      endif
!
    endsubroutine pencil_interdep_dustvelocity
!***********************************************************************
    subroutine calc_pencils_dustvelocity(f,p)
!
!  Calculate Dustvelocity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  13-nov-04/anders: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
!
      real, dimension (nx,3,3) :: tmp_pencil_3x3
      integer :: i,j
!
      intent(in) :: f
      intent(inout) :: p
! uud
      if (lpencil(i_uud)) p%uud(:,:,1)=f(l1:l2,m,n,iudx(1):iudz(1))
! ud2
      if (lpencil(i_ud2)) call dot2_mn(p%uud(:,:,1),p%ud2(:,1))
! udij
      if (lpencil(i_udij)) call gij(f,iuud(1),p%udij,1)
! divud
      if (lpencil(i_divud)) &
          p%divud(:,1) = p%udij(:,1,1,1) + p%udij(:,2,2,1) + p%udij(:,3,3,1)
! udgud      
      if (lpencil(i_udgud)) call multmv_mn(p%udij,p%uud(:,:,1),p%udgud)
! ood
      if (lpencil(i_ood)) then
        p%ood(:,1,1)=p%udij(:,3,2,1)-p%udij(:,2,3,1)
        p%ood(:,2,1)=p%udij(:,1,3,1)-p%udij(:,3,1,1)
        p%ood(:,3,1)=p%udij(:,2,1,1)-p%udij(:,1,2,1)
      endif
! od2
      if (lpencil(i_od2)) call dot2_mn(p%ood(:,:,1),p%od2(:,1))
! oud
      if (lpencil(i_oud)) call dot_mn(p%ood(:,:,1),p%uud(:,:,1),p%oud(:,1))
! sdij
      if (lpencil(i_sdij)) then
        select case (iviscd)
        case ('nud-const')
          do j=1,3
            do i=1,3
              p%sdij(:,i,j,1)=.5*(p%udij(:,i,j,1)+p%udij(:,j,i,1))
            enddo
            p%sdij(:,j,j,1)=p%sdij(:,j,j,1)-.333333*p%divud(:,1)
          enddo
        case ('hyper3_nud-const')
          call gij(f,iuud(1),tmp_pencil_3x3,5)
          do i=1,3
            do j=1,3
              p%sdij(:,i,j,1)=tmp_pencil_3x3(:,i,j)
            enddo
          enddo
        case default
          if (headtt) then
            write (unit=errormsg,fmt=*) 'No rate-of-strain tensor matches iviscd=', iviscd
            call warning('calc_pencils_dustvelocity',errormsg)
          endif
        endselect
      endif
! del2ud
      if (lpencil(i_del2ud)) call del2v(f,iuud(1),p%del2ud(:,:,1))
! del6ud
      if (lpencil(i_del6ud)) call del6v(f,iuud(1),p%del6ud(:,:,1))
! graddivud          
      if (lpencil(i_graddivud)) &
          call del2v_etc(f,iuud(1),GRADDIV=p%graddivud(:,:,1))
!
      if (.not.lhydro) p%uu=p%uud(:,:,1)
!
    endsubroutine calc_pencils_dustvelocity
!***********************************************************************
    subroutine duud_dt(f,df,p)
!
!  Dust velocity evolution
!  Calculate duud/dt = - uud.graduud - 2Omega x uud - 1/tausd*(uud-uu)
!
!  18-mar-03/axel+anders: adapted from hydro
!
      use Cdata
      use General
      use Sub
      use EquationOfState, only: gamma, beta_glnrho_global, beta_glnrho_scaled
      use IO
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!      
      real, dimension (nx,3) :: fviscd,tausd13,tausg13,del2ud,AA_sfta,BB_sfta
      real, dimension (nx) :: csrho,tausg1,mudrhod1
      real :: c2,s2 !(coefs for Coriolis force with inclined Omega)
      integer :: i,j,l
!
      intent(in) :: f,p
      intent(out) :: df
!
!  Identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'duud_dt: SOLVE duud_dt'
      if (headtt) then
        call identify_bcs('udx',iudx(1))
        call identify_bcs('udy',iudy(1))
        call identify_bcs('udz',iudz(1))
      endif
!
!  Advection term
!
      if (ladvection_dust) df(l1:l2,m,n,iudx(1):iudz(1)) = &
          df(l1:l2,m,n,iudx(1):iudz(1)) - p%udgud(:,:,1)
!
!  Coriolis force, -2*Omega x ud
!  Omega=(-sin_theta, 0, cos_theta)
!  theta corresponds to latitude
!
      if (Omega/=0. .and. lcoriolisforce_dust) then
        if (theta==0) then
          if (headtt) &
              print*,'duud_dt: add Coriolis force; Omega=',Omega
          c2=2*Omega
          df(l1:l2,m,n,iudx(1)) = df(l1:l2,m,n,iudx(1)) + c2*p%uud(:,2,1)
          df(l1:l2,m,n,iudy(1)) = df(l1:l2,m,n,iudy(1)) - c2*p%uud(:,1,1)
        else
          if (headtt) print*, &
              'duud_dt: Coriolis force; Omega,theta=',Omega,theta
          c2=2*Omega*cos(theta*pi/180.)
          s2=2*Omega*sin(theta*pi/180.)
          df(l1:l2,m,n,iudx(1)) = &
              df(l1:l2,m,n,iudx(1)) + c2*p%uud(:,2,1)
          df(l1:l2,m,n,iudy(1)) = &
              df(l1:l2,m,n,iudy(1)) - c2*p%uud(:,1,1) + s2*p%uud(:,3,1)
          df(l1:l2,m,n,iudz(1)) = &
              df(l1:l2,m,n,iudz(1))                   + s2*p%uud(:,2,1)
        endif
      endif
!
!  Add pressure gradient force
!
      if (lpressure_gradient) then
        do j=1,3
          if (lentropy) then
            df(l1:l2,m,n,iudx(1)-1+j) = df(l1:l2,m,n,iudx(1)-1+j) &
                - 1/(1+p%epsd(:,1))*p%cs2* &
                (p%glnrho(:,j) + p%cp1tilde*p%gss(:,j))
          else
            df(l1:l2,m,n,iudx(1)-1+j) = df(l1:l2,m,n,iudx(1)-1+j) &
                - 1/(1+p%epsd(:,1))*p%cs2*p%glnrho(:,j)/gamma
          endif
        enddo
!
!  Add pressure force from global density gradient.
!  
        if (maxval(abs(beta_glnrho_global))/=0.0) then
          if (headtt) print*, 'duud_dt: adding global pressure gradient force'
          do j=1,3
            df(l1:l2,m,n,(iudx(1)-1)+j) = df(l1:l2,m,n,(iudx(1)-1)+j) &
                - 1/(1+p%epsd(:,1))*1/gamma*p%cs2*beta_glnrho_scaled(j)
          enddo
        endif
      endif
!
!  Add viscosity on dust
!
      if (lviscosity_dust) then
!
        fviscd=0.0
        diffus_nud=0.0  ! Do not sum viscosity from all dust species
!
        select case (iviscd)
!
!  Viscous force: nud*del2ud
!     -- not physically correct (no momentum conservation)
!
        case('simplified')
          if (headtt) print*, 'Viscous force (dust): nud*del2ud'
          fviscd = fviscd + nud*p%del2ud(:,:,1)
          if (lfirst.and.ldt) diffus_nud=diffus_nud+nud*dxyz_2
!
!  Viscous force: nud*(del2ud+graddivud/3+2Sd.glnnd)
!    -- the correct expression for nud=const
!
        case('nud-const')
          if (headtt) print*, &
              'Viscous force (dust): nud*(del2ud+graddivud/3+2Sd.glnnd)'
          if (ldustdensity) then
            fviscd = fviscd + 2*nud*p%sdglnnd(:,:,1) + &
                nud*(p%del2ud(:,:,1)+1/3.*p%graddivud(:,:,1))
          else
            fviscd = fviscd + nud*(p%del2ud(:,:,1)+1/3.*p%graddivud(:,:,1))
          endif
          if (lfirst.and.ldt) diffus_nud=diffus_nud+nud*dxyz_2
!
!  Viscous force: nud*del6ud (not momentum-conserving)
!
        case('hyper3_simplified')
          if (headtt) print*, 'Viscous force (dust): nud*del6ud'
          fviscd = fviscd + nud*p%del6ud(:,:,1)
          if (lfirst.and.ldt) diffus_nud=diffus_nud+nud*dxyz_6

        case('hyper3_rhod_nud-const')
!
!  Viscous force: mud/rhod*del6ud
!
          if (headtt) print*, 'Viscous force (dust): mud/rhod*del6ud'
          mudrhod1=nud/p%rhod(:,1)   ! = mud/rhod
          do i=1,3
            fviscd(:,i) = fviscd(:,i) + mudrhod1*p%del6ud(:,i,1)
          enddo
          if (lfirst.and.ldt) diffus_nud=diffus_nud+nud*dxyz_6

        case('hyper3_nud-const')
!
!  Viscous force: nud*(del6ud+S.glnnd), where S_ij=d^5 ud_i/dx_j^5
!
          if (headtt) print*, 'Viscous force (dust): nud*(del6ud+S.glnnd)'
          fviscd = fviscd + nud*(p%del6ud(:,:,1)+p%sdglnnd(:,:,1))
          if (lfirst.and.ldt) diffus_nud=diffus_nud+nud*dxyz_6

        case default

          write (unit=errormsg,fmt=*) 'No such value for iviscd: ', trim(iviscd)
          call fatal_error('duud_dt',errormsg)

        endselect

      df(l1:l2,m,n,iudx(1):iudz(1)) = df(l1:l2,m,n,iudx(1):iudz(1)) + fviscd

      endif
!
!  ``uud/dx'' for timestep
!
      if (lfirst .and. ldt) then
        advec_uud=max(advec_uud,abs(p%uud(:,1,1))*dx_1(l1:l2)+ &
                                abs(p%uud(:,2,1))*dy_1(  m  )+ &
                                abs(p%uud(:,3,1))*dz_1(  n  ))
        if (idiag_dtud/=0) &
            call max_mn_name(advec_uud/cdt,idiag_dtud,l_dt=.true.)
        if (idiag_dtnud/=0) &
            call max_mn_name(diffus_nud/cdtv,idiag_dtnud,l_dt=.true.)
      endif
      if (headtt.or.ldebug) then
        print*,'duud_dt: max(advec_uud) =',maxval(advec_uud)
        print*,'duud_dt: max(diffus_nud) =',maxval(diffus_nud)
      endif
!
!  Calculate diagnostic variables
!
      if (ldiagnos) then
        if ((headtt.or.ldebug) .and. (ip<6)) &
            print*, 'duud_dt: Calculate diagnostic values...'
        if (idiag_udrms/=0) &
            call sum_mn_name(p%ud2(:,1),idiag_udrms,lsqrt=.true.)
        if (idiag_udmax/=0) &
            call max_mn_name(p%ud2(:,1),idiag_udmax,lsqrt=.true.)
        if (idiag_rdudmax/=0) &
            call max_mn_name(p%rhod(:,1)**2*p%ud2(:,1),idiag_rdudmax, &
            lsqrt=.true.)
        if (idiag_ud2m/=0) call sum_mn_name(p%ud2(:,1),idiag_ud2m)
        if (idiag_udxm/=0) call sum_mn_name(p%uud(:,1,1),idiag_udxm)
        if (idiag_udym/=0) call sum_mn_name(p%uud(:,2,1),idiag_udym)
        if (idiag_udzm/=0) call sum_mn_name(p%uud(:,3,1),idiag_udzm)
        if (idiag_udx2m/=0) &
            call sum_mn_name(p%uud(:,1,1)**2,idiag_udx2m)
        if (idiag_udy2m/=0) &
            call sum_mn_name(p%uud(:,2,1)**2,idiag_udy2m)
        if (idiag_udz2m/=0) &
            call sum_mn_name(p%uud(:,3,1)**2,idiag_udz2m)
        if (idiag_udm2/=0) call max_mn_name(p%ud2(:,1),idiag_udm2)
        if (idiag_divud2m/=0) &
            call sum_mn_name(p%divud(:,1)**2,idiag_divud2m)
        if (idiag_rdudxm/=0) &
            call sum_mn_name(p%rhod(:,1)*p%uud(:,1,1),idiag_rdudxm)
        if (idiag_rdudym/=0) &
            call sum_mn_name(p%rhod(:,1)*p%uud(:,2,1),idiag_rdudym)
        if (idiag_rdudzm/=0) &
            call sum_mn_name(p%rhod(:,1)*p%uud(:,3,1),idiag_rdudzm)
        if (idiag_rdudx2m/=0) &
            call sum_mn_name((p%rhod(:,1)*p%uud(:,1,1))**2,idiag_rdudx2m)
!
!  xy-averages
!
        if (idiag_udxmz/=0) &
            call xysum_mn_name_z(p%uud(:,1,1),idiag_udxmz)
        if (idiag_udymz/=0) &
            call xysum_mn_name_z(p%uud(:,2,1),idiag_udymz)
        if (idiag_udzmz/=0) &
            call xysum_mn_name_z(p%uud(:,3,1),idiag_udzmz)
        if (idiag_udx2mz/=0) &
            call xysum_mn_name_z(p%uud(:,1,1)**2,idiag_udx2mz)
        if (idiag_udy2mz/=0) &
            call xysum_mn_name_z(p%uud(:,2,1)**2,idiag_udy2mz)
        if (idiag_udz2mz/=0) &
            call xysum_mn_name_z(p%uud(:,3,1)**2,idiag_udz2mz)
!
!  z-averages
!
        if (idiag_udxmxy/=0) &
            call zsum_mn_name_xy(p%uud(:,1,1),idiag_udxmxy)
        if (idiag_udymxy/=0) &
            call zsum_mn_name_xy(p%uud(:,2,1),idiag_udymxy)
        if (idiag_udzmxy/=0) &
            call zsum_mn_name_xy(p%uud(:,3,1),idiag_udzmxy)
!
!  kinetic field components at one point (=pt)
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_udxpt/=0) &
              call save_name(p%uud(lpoint-nghost,1,1),idiag_udxpt)
          if (idiag_udypt/=0) &
              call save_name(p%uud(lpoint-nghost,2,1),idiag_udypt)
          if (idiag_udzpt/=0) &
              call save_name(p%uud(lpoint-nghost,3,1),idiag_udzpt)
        endif
!
!  Things related to vorticity and helicity
!
        if (idiag_odrms/=0) &
            call sum_mn_name(p%od2,idiag_odrms,lsqrt=.true.)
        if (idiag_odmax/=0) &
            call max_mn_name(p%od2,idiag_odmax,lsqrt=.true.)
        if (idiag_od2m/=0) call sum_mn_name(p%od2,idiag_od2m)
        if (idiag_oudm/=0) call sum_mn_name(p%oud,idiag_oudm)
!          
      endif
!
    endsubroutine duud_dt
!***********************************************************************
    subroutine read_dustvelocity_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
                                                                                                   
      if (present(iostat)) then
        read(unit,NML=dustvelocity_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=dustvelocity_init_pars,ERR=99)
      endif
                                                                                                   
                                                                                                   
99    return
    endsubroutine read_dustvelocity_init_pars
!***********************************************************************
    subroutine write_dustvelocity_init_pars(unit)
      integer, intent(in) :: unit
                                                                                                   
      write(unit,NML=dustvelocity_init_pars)
                                                                                                   
    endsubroutine write_dustvelocity_init_pars
!***********************************************************************
    subroutine read_dustvelocity_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
                                                                                                   
      if (present(iostat)) then
        read(unit,NML=dustvelocity_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=dustvelocity_run_pars,ERR=99)
      endif
                                                                                                   
                                                                                                   
99    return
    endsubroutine read_dustvelocity_run_pars
!***********************************************************************
    subroutine write_dustvelocity_run_pars(unit)
      integer, intent(in) :: unit
                                                                                                   
      write(unit,NML=dustvelocity_run_pars)
                                                                                                   
    endsubroutine write_dustvelocity_run_pars
!***********************************************************************
    subroutine rprint_dustvelocity(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Cdata
      use Sub
      use General, only: chn
!
      integer :: iname,inamez
      logical :: lreset,lwr
      logical, optional :: lwrite
!
!  Write information to index.pro that should not be repeated for i
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      
      if (lwr) then
        write(3,*) 'ndustspec=',ndustspec
        write(3,*) 'nname=',nname
      endif
!
!  reset everything in case of reset
!
      if (lreset) then
        idiag_dtud=0; idiag_dtnud=0; idiag_ud2m=0; idiag_udx2m=0
        idiag_udxm=0; idiag_udym=0; idiag_udzm=0
        idiag_udy2m=0; idiag_udz2m=0; idiag_udm2=0; idiag_oudm=0; idiag_od2m=0
        idiag_udxpt=0; idiag_udypt=0; idiag_udzpt=0; idiag_udrms=0
        idiag_udmax=0; idiag_odrms=0; idiag_odmax=0; idiag_rdudmax=0
        idiag_udmx=0; idiag_udmy=0; idiag_udmz=0; idiag_divud2m=0
        idiag_epsKd=0; idiag_rdudxm=0;idiag_rdudym=0; idiag_rdudzm=0;
        idiag_rdudx2m=0; idiag_udx2mz=0; idiag_udy2mz=0; idiag_udz2mz=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'rprint_dustvelocity: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname), &
            'dtud',idiag_dtud)
        call parse_name(iname,cname(iname),cform(iname), &
            'dtnud',idiag_dtnud)
        call parse_name(iname,cname(iname),cform(iname), &
            'udxm',idiag_udxm)
        call parse_name(iname,cname(iname),cform(iname), &
            'udym',idiag_udym)
        call parse_name(iname,cname(iname),cform(iname), &
            'udzm',idiag_udzm)
        call parse_name(iname,cname(iname),cform(iname), &
            'ud2m',idiag_ud2m)
        call parse_name(iname,cname(iname),cform(iname), &
            'udx2m',idiag_udx2m)
        call parse_name(iname,cname(iname),cform(iname), &
            'udy2m',idiag_udy2m)
        call parse_name(iname,cname(iname),cform(iname), &
            'udz2m',idiag_udz2m)
        call parse_name(iname,cname(iname),cform(iname), &
            'udm2',idiag_udm2)
        call parse_name(iname,cname(iname),cform(iname), &
            'od2m',idiag_od2m)
        call parse_name(iname,cname(iname),cform(iname), &
            'oudm',idiag_oudm)
        call parse_name(iname,cname(iname),cform(iname), &
            'udrms',idiag_udrms)
        call parse_name(iname,cname(iname),cform(iname), &
            'udmax',idiag_udmax)
        call parse_name(iname,cname(iname),cform(iname), &
            'rdudmax',idiag_rdudmax)
        call parse_name(iname,cname(iname),cform(iname), &
            'rdudxm',idiag_rdudxm)
        call parse_name(iname,cname(iname),cform(iname), &
            'rdudym',idiag_rdudym)
        call parse_name(iname,cname(iname),cform(iname), &
            'rdudzm',idiag_rdudzm)
        call parse_name(iname,cname(iname),cform(iname), &
            'rdudx2m',idiag_rdudx2m)
        call parse_name(iname,cname(iname),cform(iname), &
            'odrms',idiag_odrms)
        call parse_name(iname,cname(iname),cform(iname), &
            'odmax',idiag_odmax)
        call parse_name(iname,cname(iname),cform(iname), &
            'udmx',idiag_udmx)
        call parse_name(iname,cname(iname),cform(iname), &
            'udmy',idiag_udmy)
        call parse_name(iname,cname(iname),cform(iname), &
            'udmz',idiag_udmz)
        call parse_name(iname,cname(iname),cform(iname), &
            'divud2m',idiag_divud2m)
        call parse_name(iname,cname(iname),cform(iname), &
            'epsKd',idiag_epsKd)
        call parse_name(iname,cname(iname),cform(iname), &
            'udxpt',idiag_udxpt)
        call parse_name(iname,cname(iname),cform(iname), &
            'udypt',idiag_udypt)
        call parse_name(iname,cname(iname),cform(iname), &
            'udzpt',idiag_udzpt)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'udxmz',idiag_udxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'udymz',idiag_udymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'udzmz',idiag_udzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'udx2mz',idiag_udx2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'udy2mz',idiag_udy2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'udz2mz',idiag_udz2mz)
      enddo
!
!  Write dust index in short notation
!
      if (lwr) then
        write(3,*) 'iuud=', iudx(1)
        write(3,*) 'iudx=', iudx(1)
        write(3,*) 'iudy=', iudy(1)
        write(3,*) 'iudz=', iudz(1)
      endif
!
    endsubroutine rprint_dustvelocity
!***********************************************************************

endmodule Dustvelocity
