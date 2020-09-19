! $Id: streamfunction.f90 19193 2012-06-30 12:55:46Z wdobler $
!
!  This module solves for a streamfunction for the 2D velocity field.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 1 
! COMMUNICATED AUXILIARIES 0
!
!***************************************************************
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
  use Deriv
!
  implicit none
!
  include '../special.h'
!
!  Global to this method
!
  real, dimension(mx,my,mz) :: psi
!
!  Code constants
!
  real :: amplpsi=1d-5   !amplitude of initial perturbation for psi
  real :: tolerance=1d-6,tolerance_coarsest=impossible
  real :: alpha_sor=impossible
!
!  Variables for calculating the Rayleigh number. Default are Europa values.
!
  real ::  gravity_z = 1.3  ! gravity in m/s^2
  real ::  rho0_bq = 917.   ! density in Kg/m^3, Boussinesq value
  real ::  alpha_thermal = 1.65d-4  ! Thermal expansivity in K^-1
  real ::  kappa = 1d-6     ! Thermal diffusivity in m^2/s
  !real ::  cp = 2000.      ! Specific heat in J Kg^-1 K^-1
  real ::  eta_0 = 1d13     ! Viscosity at melting temperature in Pa s
  !real ::  dslab_km = 20.  ! Ice shell thikness in km
  real :: Tbot=270, Tupp=100
!
!  Variables for variable viscosity
!
  real ::  TT_melt = 270.   ! Melting tempertature in K
  real ::  Avisc = 4.00     ! constants for viscosity prescriptions
  real ::  Bvisc = 0.00
  real ::  Cvisc = 0.00
  real ::  deltaT1, Lz1
!
!  Variables for tidal heating
!
  real :: mu_ice = 4d9     ! Rigidity of ice in Pa s
  real :: epsi_0 = 2.1d-5  ! Amplitude of original tidal flexing
  real :: EuropaPeriod = 3.0682204d5 ! Period of Europa, in seconds  
  real :: mu1_ice, OmegaEuropa, kappa1
!
!  These are the needed internal "pencils".
!
  type InternalPencils
     real, dimension(nx,3)   :: uu
     real, dimension(nx) :: ugTT,u2,qtidal,eta,devsigzz
  endtype InternalPencils
!
  type (InternalPencils) :: q
!
  character (len=labellen) :: initpsi='nothing'
  character (len=labellen) :: iconv_viscosity='Newtonian'
  character (len=labellen) :: iorder_sor_solver='high_order'
!
  logical :: lprint_residual=.false.,ltidal_heating=.true.
  logical :: ltemperature_advection=.true.,ltemperature_diffusion=.true.
  logical :: lmultigrid=.true.
!  
  integer :: maxit=1000
!
  real :: Ra=impossible
  logical :: lsave_residual=.true.
  real :: kx_TT=2*pi, kz_TT=pi, ampltt=0.01
  logical :: lpoisson_test=.false.
  logical :: ldirect_solver=.true.
!
  integer :: npost=5,npre=5
  integer :: igamma=1,n_vcycles=1
  integer :: nx_coarsest=3
!
  logical :: lsave_residual_grid=.false.
  logical :: ladimensional=.true.
  logical :: lsplit_temperature=.false.
!
  namelist /special_init_pars/ amplpsi,alpha_sor,lprint_residual,&
       tolerance,maxit,gravity_z,rho0_bq,alpha_thermal,kappa,eta_0,&
       iconv_viscosity,Avisc,Bvisc,Cvisc,&
       Tbot,Tupp,Ra,iorder_sor_solver,lsave_residual,&
       kx_TT,kz_TT,ampltt,initpsi,lmultigrid,lpoisson_test,npost,npre,igamma,n_vcycles,&
       ldirect_solver,nx_coarsest,lsave_residual_grid,ladimensional,lsplit_temperature
!
  namelist /special_run_pars/ amplpsi,alpha_sor,Avisc,lprint_residual,&
       tolerance,maxit,gravity_z,rho0_bq,alpha_thermal,kappa,eta_0,&
       iconv_viscosity,Avisc,Bvisc,Cvisc,&
       Tbot,Tupp,Ra,iorder_sor_solver,lsave_residual,&
       ltidal_heating,ltemperature_advection,ltemperature_diffusion,lmultigrid,&
       npost,npre,igamma,n_vcycles,ldirect_solver,nx_coarsest,lsave_residual_grid,&
       lsplit_temperature
!
!  Declare index of new variables in f array. Surface density, midplane
!  temperature, and mass accretion rate.
!
  integer :: ipsi=0, iresidual=0
!
!  Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_uqxmin=0, idiag_uqxmax=0, idiag_uqxrms=0, idiag_uqxm=0, idiag_uqx2m=0
  integer :: idiag_uqzmin=0, idiag_uqzmax=0, idiag_uqzrms=0, idiag_uqzm=0, idiag_uqz2m=0
  integer :: idiag_uq2m=0, idiag_uqrms=0, idiag_uqmax=0
  integer :: idiag_qtidalmin=0, idiag_qtidalmax=0, idiag_qtidalm=0
  integer :: idiag_dTdz1=0, idiag_dTdz2=0, idiag_dTdz3=0, idiag_dTdz4=0 
  integer :: idiag_nusselt_num=0, idiag_nusselt_den=0
  integer :: idiag_TTmax_cline=0, idiag_TTmin_cline=0
  integer :: idiag_devsigzz1=0,idiag_devsigzz2=0,idiag_devsigzz3=0,idiag_devsigzz4=0
  integer :: idiag_icount=0,idiag_residual=0
!
  logical :: lsolver_highorder, lsolver_loworder
!
  integer :: icount_save
  real :: residual_save
!  
  logical :: lviscosity_const
  logical :: lviscosity_var_newtonian
  logical :: lviscosity_var_blankenbach
!
  real, dimension(mz) :: Tbar = 0.0, etabar=0.0
  real :: gTT_conductive = 0.0, eta_alpha2=0.0
!
  contains
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!  01-aug-11/wlad: adapted
!
      use FArrayManager, only: farray_register_pde,farray_register_auxiliary
!
      if (lroot) call svn_id( &
          "$Id: alphadisk.f90 19193 2012-06-30 12:55:46Z wdobler $")
!
!  Register the streamfunction. The streamfunction needs to be a pde so it
!  can have boundaries. 
!
      call farray_register_pde('psi',ipsi)
      call farray_register_auxiliary('residual',iresidual)
!
! Write to read special variables with pc_varcontent
!
      if (lroot) then
        open(4,file=trim(datadir)//'/index_special.pro',status='replace')
        write(4,*) 'ipsi= ',ipsi
        write(4,*) 'iresidual= ',iresidual
        close(4)
      endif
!
      if (lsave_residual) then
         open(9,file='residuals.dat',status='replace')
      endif
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!  01-aug-11/wlad: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: dslab,delta_T,delta
!
      if (Lxyz(1)/(nxgrid-1) .ne. Lxyz(3)/(nzgrid-1)) then
        print*,Lxyz(1)/(nxgrid-1),Lxyz(3)/(nzgrid-1)
        call fatal_error("initialize_special","hx ne hz")
      endif
      !if (Lxyz(1)/nxgrid .ne. Lxyz(3)/nzgrid) then
      !   print*,Lxyz(1)/nxgrid,Lxyz(3)/nzgrid
      !  call fatal_error("initialize_special","dx ne dz")
      !endif
!
!  Pre-calculate Rayleigh number in case it is not given
!
      delta_T=Tbot-Tupp
      deltaT1=1./delta_T
      dslab = Lxyz(3)
      Lz1 = 1./Lxyz(3)
!
!  Conductive (linear) Temperature Gradient, used if the T=Tbar(z) + theta(x,y,z,t)
!  split is used. It similarly defines the conductive viscosity eta = etabar(z) * nu(x,y,z,t).
!
      if (lsplit_temperature) then
         do n=n1,n2
            Tbar(n) = Tbot + (n-n1)*(Tupp-Tbot)/(n2-n1)
         enddo
         etabar = exp(-Bvisc*Tbar*deltaT1)
         gTT_conductive = -delta_T*Lz1
         eta_alpha2 = (Bvisc*deltaT1*gTT_conductive)**2
      endif
!
       if (tolerance_coarsest==impossible) tolerance_coarsest=tolerance
!
!  Stuff for the coefficient of SOR. Should be between 0 and 2 for stability.
!     
      if (alpha_sor == impossible) then
        delta=min(dx,dy,dz)/dslab
        alpha_sor= 2./(1+sin(pi*delta))
      endif
!
!  Europa-specific stuff
!
      mu1_ice=1./mu_ice
      OmegaEuropa=1./EuropaPeriod
!
      kappa1=1./kappa
!
      if (ladimensional) then
         if (eta_0 /= 1.0)   call fatal_error("initialize_special",&
              "You are using dimensionless variables, set eta_0=1 and vary Ra")
         if (kappa /= 1.0)   call fatal_error("initialize_special",&
              "You are using dimensionless variables, set kappa=1 and vary Ra")
         if (delta_T /= 1.0) call fatal_error("initialize_special",&
              "You are using dimensionless variables, set Tupp=0 and Tbot=1 and vary Ra")
         if (dslab /= 1.0)   call fatal_error("initialize_special",&
              "You are using dimensionless variables, set Lxyz(3)=1 and vary Ra")
      else
         if (lstart.and.Ra/=impossible) call fatal_error("initialize_special",&
             "You are using dimensioned variables but you set the Rayleigh number. Uncomment it from your start.in")
      endif
!
      if (Ra==impossible) then
        Ra = (gravity_z*alpha_thermal*rho0_bq*delta_T*dslab**3)/(kappa*eta_0)
      endif
!
      if (lroot) print*,'Rayleigh number=',Ra
!
      select case (iorder_sor_solver)
!
      case ('low_order')
         lsolver_loworder=.true.
         lsolver_highorder=.false.
      case ('high_order') 
         lsolver_loworder=.false.
         lsolver_highorder=.true.
      case default  
        write(unit=errormsg,fmt=*) &
             'initialize_special: No such value for iorder_sor_solver: ', &
             trim(iorder_sor_solver)
        call fatal_error('initialize_special',errormsg)
      endselect
!
      select case (iconv_viscosity)
!
      case ('constant')
         lviscosity_const=.true.
         lviscosity_var_newtonian=.false.
         lviscosity_var_blankenbach=.false.
      case ('Newtonian') 
         lviscosity_const=.false.
         lviscosity_var_newtonian=.true.
         lviscosity_var_blankenbach=.false.
      case ('Blankenbach-variable') 
         lviscosity_const=.false.
         lviscosity_var_newtonian=.false.
         lviscosity_var_blankenbach=.true.
      case default
        write(unit=errormsg,fmt=*) &
             'initialize_special: No such value for iconv_viscosity: ', &
             trim(iconv_viscosity)
        call fatal_error('initialize_special',errormsg)
      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  Initialise special condition; called from start.f90.
!
!  06-oct-2003/tony: coded
!  01-aug-11/wlad: adapted
!
      use Initcond, only: gaunoise
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: amplpsi_
!
!  Give an initial guess for psi
!      
      select case (initpsi)
!
      case ('noise')
         print*,'init_special: gaussian noise initialization'
         call gaunoise(amplpsi,f,ipsi)
      case ('single-mode')
         if (.not.ladimensional) call fatal_error("init_special",&
              "not yet coded for dimensioned variables. Use initpsi='noise' instead")
         print*,'init_special: single-mode initialization for constant viscosity'
         print*,'derived from the solution of the dispersion relation'
         amplpsi_ = -ampltt * Ra*kx_TT/(kz_TT**2 + kx_TT**2)**2
         do n=n1,n2
           do m=m1,m2
             f(l1:l2,m,n,ipsi) = f(l1:l2,m,n,ipsi) + &
                  amplpsi_*sin(kx_TT*x(l1:l2))*sin(kz_TT*z(n))
           enddo
         enddo
         print*,'init_special: done single-mode initialization'
!
      case ('zero')
         f(l1:l2,m1:m2,n1:n2,ipsi)=0.
!         
      case ('nothing')
                   
      case default
!
!  Catch unknown values.
!
          write(unit=errormsg,fmt=*) 'No such value initpsi = ',trim(initpsi)
          call fatal_error('init_special',errormsg)
!
      endselect
      !
      call special_after_boundary(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
      if (ltemperature_advection) lpenc_requested(i_gTT)=.true.
      if (ltemperature_diffusion) lpenc_requested(i_del2TT)=.true.
      if (ltidal_heating) lpenc_requested(i_TT1)=.true.
!
      if (idiag_dTdz1/=0 .or. &
          idiag_dTdz2/=0 .or. &
          idiag_dTdz3/=0 .or. &
          idiag_dTdz4/=0 .or. &
          idiag_nusselt_num/=0) lpenc_diagnos(i_gTT)=.true.
!
      if (idiag_nusselt_den/=0 .or. &
          idiag_TTmax_cline/=0 .or. &
          idiag_TTmin_cline/=0) lpenc_diagnos(i_TT)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils), intent(inout) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Special pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-nov-04/tony: coded
!  01-aug-11/wlad: adapted
!
      use Deriv, only: der,derij
      use Sub, only: u_dot_grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (nx) :: derxzpsi
!
      intent(inout) :: f
      intent(inout) :: p
!
      if (ltidal_heating .or. &
          idiag_devsigzz1 /=0 .or. &
          idiag_devsigzz2 /=0 .or. &
          idiag_devsigzz3 /=0 .or. &
          idiag_devsigzz4 /=0) then 
         if (lviscosity_const) then 
            q%eta=eta_0
         else if (lviscosity_var_newtonian) then
            q%eta=eta_0*exp(Avisc*(TT_melt*p%TT1 - 1.))
         else if (lviscosity_var_blankenbach) then
            q%eta=eta_0*exp(-Bvisc*p%TT*deltaT1 + Cvisc*(1-z(n))*Lz1)
         else   
            call fatal_error("calc_pencils_special",&
                 "The world is flat, and we never got here.")
         endif
      endif
!
!  Calculate velocities via streamfunction. 
!
      if (ltemperature_advection) then 
         call der( f,ipsi,q%uu(:,1),3)
         q%uu(:,2)=0.
         call der(-f,ipsi,q%uu(:,3),1)
!
!  Calculate temperature advection term. 
!  
         call u_dot_grad(f,iTT,p%gTT,q%uu,q%ugTT)
      endif
!
      if (ltidal_heating) then
        call fatal_error("calc_pencils_special","these are not yet normalized")
        q%qtidal = .5*q%eta*(epsi_0*OmegaEuropa)**2 / (1+(OmegaEuropa*q%eta*mu1_ice)**2)
      endif
!
      if (idiag_uq2m/=0.or.idiag_uqrms/=0.or.idiag_uqmax/=0) &
           q%u2=q%uu(:,1)**2 + q%uu(:,3)**2
!
!  duz/dz for topography calculation. duz/dz = d/dz (-dpsi/dx) = -d2dxdz psi = derij(psi,x,z)
!      
      if (idiag_devsigzz1/=0.or.idiag_devsigzz2/=0.or.idiag_devsigzz3/=0.or.idiag_devsigzz4/=0) then 
         call derij(f,ipsi,derxzpsi,1,3)
         q%devsigzz = q%eta * derxzpsi
      endif
!
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine special_after_boundary(f)
!
      !use Boundcond, only: update_ghosts
!
      real, dimension (mx,my,mz,mfarray) :: f   
!    
      call solve_for_psi(f)
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine solve_for_psi(f)
!
      real, dimension (mx,my,mz,mfarray) :: f   
      real, dimension (mx,mz) :: psi,psi_old,eta
!
      real, dimension (nx,nz) :: rhs,tmp
      real, dimension (nx,nz) :: alpha_factor,beta_factor
!
      real :: residual,aout,dTTdx
      real :: alpha=impossible,beta=impossible
      real :: Rayleigh_factor
!
      integer :: icount,i
!
      if (ladimensional) then
         Rayleigh_factor=Ra
      else
         Rayleigh_factor=alpha_thermal*rho0_bq*gravity_z
      endif
!
!  Define r.h.s.
!
      if (lpoisson_test) then
         psi=0. 
         rhs=1.
      else
!
         call calc_viscosity(f,eta)
!
         do n=n1,n2; do m=m1,m2
            do i=l1,l2
               call f_function(eta,aout,i,n) ; alpha = aout/eta(i,n)
               call g_function(eta,aout,i,n) ; beta  = aout/eta(i,n)
!
               dTTdx=dx_1(i)/60*(+ 45.0*(f(i+1,m,n,iTT)-f(i-1,m,n,iTT)) &
                                 -  9.0*(f(i+2,m,n,iTT)-f(i-2,m,n,iTT)) &
                                 +      (f(i+3,m,n,iTT)-f(i-3,m,n,iTT)))            
!
               if (lsplit_temperature) then 
                 rhs(i-l1+1,n-n1+1) = Rayleigh_factor*dTTdx/(eta(i,n)*etabar(n))
                 alpha_factor(i-l1+1,n-n1+1)=alpha + eta_alpha2
               else
                 rhs(i-l1+1,n-n1+1) = Rayleigh_factor*dTTdx/eta(i,n)
                 alpha_factor(i-l1+1,n-n1+1)=alpha
               endif
               beta_factor(i-l1+1,n-n1+1)=beta
!
            enddo;enddo
         enddo
!         
         psi(l1:l2,n1:n2)=f(l1:l2,mpoint,n1:n2,ipsi)
      
         call update_bounds_psi(psi)
      endif
!
      if (lmultigrid) then
         call multigrid(f,psi,rhs,alpha_factor,beta_factor)
         call update_bounds_psi(psi)
      else
!
!  Initial residual
!
         residual=1e33
!
!  Start counter
!      
         icount=0
!         
         do while (residual > tolerance)
!
!  Calculate psi via SOR
!
            psi_old=psi
!
            call successive_over_relaxation(psi,rhs,alpha_factor,beta_factor)
!
            call update_bounds_psi(psi)
!
            tmp=(psi(l1:l2,n1:n2) - psi_old(l1:l2,n1:n2))**2
            residual = sqrt(sum(tmp)/sum(psi(l1:l2,n1:n2)**2))
!
            if (lsave_residual_grid) then
               do m=1,my
                  f(:,m,:,iresidual) = psi-psi_old
               enddo
            endif
!
! Increase counter. 
!
            icount=icount+1
            if (lprint_residual) &
                 print*, icount,residual,tolerance
            if (lsave_residual) then
               write(9,*) icount,residual
            endif
!
            if (icount >= maxit) then
               call fatal_error("solve_for_psi","reached limit of iterations: did not converge.")
               if (lsave_residual) close(9)
               goto 200
            endif
!
         enddo
!
200      continue
         
         if (lsave_residual) close(9)
!
         icount_save=icount
         residual_save=residual
      endif
!
      do m=m1,m2
         f(:,m,:,ipsi)=psi
      enddo
!
    endsubroutine solve_for_psi
!***********************************************************************
    subroutine successive_over_relaxation(u,r,alp,bet,h)
!
      real, dimension (:,:) :: u
      real, dimension (:,:) :: r
      real, dimension (:,:) :: alp,bet
      real :: cc,ufactor,vterm,u_new,alpha_sor_
      real :: alpha=impossible,beta=impossible
      real, optional :: h
      integer :: ic,nc,k,l
      integer :: nux,nuz,nrx,nrz,nn1,nn2,ll1,ll2,lcc1,lcc2
!
      nux=size(u,1); nuz=size(u,2)
      nrx=size(r,1); nrz=size(r,2)
      ll1=nghost+1; ll2=nux-nghost
      nn1=nghost+1; nn2=nuz-nghost
!
      if (present(h)) then 
         alpha_sor_=2./(1+sin(pi*h))
      else
         alpha_sor_=alpha_sor
      endif
!
      if (lperi(1)) then
         lcc1=ll1
         lcc2=ll2
      else
         lcc1=ll1+1
         lcc2=ll2-1
      endif
!
      do nc=nn1+1,nn2-1
         k=nc-nn1+1
         do ic=lcc1,lcc2
           l=ic-ll1+1
!     
           cc    =   r(l,k)
           alpha = alp(l,k)
           beta  = bet(l,k)
!     
           if (lsolver_highorder) then
              if (present(h)) then
                 call solve_highorder(u,alpha,beta,ufactor,vterm,ic,nc)
              else
                 call solve_highorder(u,alpha,beta,ufactor,vterm,ic,nc,h)
              endif
           else
              if (present(h)) then
                 call solve_loworder(u,alpha,beta,ufactor,vterm,ic,nc)
              else
                 call solve_loworder(u,alpha,beta,ufactor,vterm,ic,nc,h)
              endif
           endif
!     
           u_new=(cc-vterm)/ufactor
           u(ic,nc) = u(ic,nc)*(1-alpha_sor_) +  alpha_sor_*u_new
!     
        enddo
     enddo
!
    endsubroutine successive_over_relaxation
!***********************************************************************
    subroutine solve_highorder(a,alpha,beta,u,v,i,k,h)
!
!   Solve the momentum equation of mantle convection for the streamfunction
!      
!   F2(psi,i,k,h) + G2(psi,i,k,h) + A*F(psi,i,k,h) + B*G(psi,i,k,h) = C
!         
!   d4/dx4 + d4/dz4 + 2*d4dx2dz2 + A*d2/dz2 - A*d2dx2 + B*d2dxdz = C 
!         
      real, dimension (:,:), intent(in) :: a
      real, intent(in) :: alpha,beta
      real, intent(out) :: u,v
      real :: dx1,dz1,fac,fac2
      real :: d4dx4,d4dz4,d2dx2,d2dz2,d2dxdz,d4dx2dz2
!
      integer :: i,k
!
      real, optional :: h
!
      if (present(h)) then
         dx1 = 1./h
         dz1 = 1./h
      else
         dx1=dx_1(i)
         dz1=dz_1(k)
      endif
!      
!  The quantity u collects all terms that multiply psi[i,k]
!  v collects everything else. 
!      
      u = 0.
      v = 0.
!           
! d4dx4
!
      fac=(1.0/6)*dx1**4
      d4dx4=fac*(+ 56.0* a(i,k)               &
                 - 39.0*(a(i+1,k)+a(i-1,k)) &
                 + 12.0*(a(i+2,k)+a(i-2,k)) &
                 -      (a(i+3,k)+a(i-3,k)))
!
! fac2 is the factor by which this term is multiplied in the momentum equation. 
!
      fac2 = 1.
      u = u + fac2*fac*56.
      v = v + fac2*d4dx4 
!
! d4dz4
!
      fac=(1.0/6)*dz1**4
      d4dz4=fac*(+ 56.0* a(i,k)               &
                 - 39.0*(a(i,k+1)+a(i,k-1)) &
                 + 12.0*(a(i,k+2)+a(i,k-2)) &
                 -      (a(i,k+3)+a(i,k-3)))

      fac2 = 1.
      u = u + fac2*fac*56.
      v = v + fac2*d4dz4
!
! d2dx2         
!
      fac=(1./180)*dx1**2
      d2dx2=fac*(-490.0* a(i,k)               &
                 +270.0*(a(i+1,k)+a(i-1,k)) &
                 - 27.0*(a(i+2,k)+a(i-2,k)) &
                 +  2.0*(a(i+3,k)+a(i-3,k)))
          
      fac2 = -alpha
      u = u + fac2*fac*(-490)
      v = v + fac2*d2dx2
!
! d2dz2
!          
      fac=(1./180)*dz1**2
      d2dz2=fac*(-490.0* a(i,k) &
                 +270.0*(a(i,k+1)+a(i,k-1)) &
                 - 27.0*(a(i,k+2)+a(i,k-2)) &
                 +  2.0*(a(i,k+3)+a(i,k-3)))
          
      fac2 = alpha 
      u = u + fac2*fac*(-490)
      v = v + fac2*d2dz2
!                    
! d2dxdz
!          
      fac=(1./60.**2)*dx1*dz1
      d2dxdz=fac*( &
           45.*((45.*(a(i+1,k+1)-a(i-1,k+1))   &
                 -9.*(a(i+2,k+1)-a(i-2,k+1))   &
                    +(a(i+3,k+1)-a(i-3,k+1)))  &
               -(45.*(a(i+1,k-1)-a(i-1,k-1))   &
                 -9.*(a(i+2,k-1)-a(i-2,k-1))   &
                    +(a(i+3,k-1)-a(i-3,k-1)))) &
           -9.*((45.*(a(i+1,k+2)-a(i-1,k+2))   &
                 -9.*(a(i+2,k+2)-a(i-2,k+2))   &
                    +(a(i+3,k+2)-a(i-3,k+2)))  &
               -(45.*(a(i+1,k-2)-a(i-1,k-2))   &
                 -9.*(a(i+2,k-2)-a(i-2,k-2))   &
                    +(a(i+3,k-2)-a(i-3,k-2)))) &
              +((45.*(a(i+1,k+3)-a(i-1,k+3))   &
                 -9.*(a(i+2,k+3)-a(i-2,k+3))   &
                    +(a(i+3,k+3)-a(i-3,k+3)))  &
               -(45.*(a(i+1,k-3)-a(i-1,k-3))   &
                 -9.*(a(i+2,k-3)-a(i-2,k-3))   &
                    +(a(i+3,k-3)-a(i-3,k-3)))) &
           )         
          
      fac2 = beta
      u = u + fac2*fac*0.
      v = v + fac2*d2dxdz
!
! d4dx2dz2       
!
      fac=(1./180.**2)*dx1**2*dz1**2
      d4dx2dz2=fac*( &
           - 490*(  -490.0* a(i,k)                   &
                    +270.0*(a(i+1,k)+a(i-1,k))       &
                    - 27.0*(a(i+2,k)+a(i-2,k))       &
                    +  2.0*(a(i+3,k)+a(i-3,k)))      &
           + 270.*((-490.0* a(i,k+1)                 &
                    +270.0*(a(i+1,k+1)+a(i-1,k+1))   &
                    - 27.0*(a(i+2,k+1)+a(i-2,k+1))   &
                    +  2.0*(a(i+3,k+1)+a(i-3,k+1)))  &
                  +(-490.0* a(i,k-1)                 &
                    +270.0*(a(i+1,k-1)+a(i-1,k-1))   &
                    - 27.0*(a(i+2,k-1)+a(i-2,k-1))   &
                    +  2.0*(a(i+3,k-1)+a(i-3,k-1)))) &
             -27.*((-490.0* a(i,k+2)                 &
                    +270.0*(a(i+1,k+2)+a(i-1,k+2))   &
                    - 27.0*(a(i+2,k+2)+a(i-2,k+2))   &
                    +  2.0*(a(i+3,k+2)+a(i-3,k+2)))  &
                  +(-490.0* a(i,k-2)                 &
                    +270.0*(a(i+1,k-2)+a(i-1,k-2))   &
                    - 27.0*(a(i+2,k-2)+a(i-2,k-2))   &
                    +  2.0*(a(i+3,k-2)+a(i-3,k-2)))) &
              +2.*((-490.0* a(i,k+3)                 &
                    +270.0*(a(i+1,k+3)+a(i-1,k+3))   &
                    - 27.0*(a(i+2,k+3)+a(i-2,k+3))   &
                    +  2.0*(a(i+3,k+3)+a(i-3,k+3)))  &
                  +(-490.0* a(i,k-3)                 &
                    +270.0*(a(i+1,k-3)+a(i-1,k-3))   &
                    - 27.0*(a(i+2,k-3)+a(i-2,k-3))   &
                    +  2.0*(a(i+3,k-3)+a(i-3,k-3)))) &
           )
!
      fac2 = 2.
      u = u + fac2*fac*(490.**2)
      v = v + fac2*d4dx2dz2
!
! Must remove u*a(i,k) from v 
!
      v = v - u*a(i,k)
!
    endsubroutine solve_highorder
!***********************************************************************
    subroutine solve_loworder(a,alpha,beta,u,v,i,k,h)
!
      real, dimension (:,:), intent(in) :: a
      real, intent(in) :: alpha,beta
      real, intent(out) :: u,v
!
      real :: v_1, v_2, v_3, v_4, v_5
      real :: v51,v52,v53,dx1,dz1
      integer :: i,k
!
      real, optional :: h
!
      if (present(h)) then
         dx1 = 1./h
         dz1 = 1./h
      else
         dx1=dx_1(i)
         dz1=dz_1(k)
      endif
!
! a is updated on the fly. So, do it swiping, so that a(i-1), that is updated, gets used for a(i)
!
      u = 20.*(dx1**2 * dz1**2)
!
      v_1 = alpha*(a(i,k+1) + a(i,k-1) - a(i+1,k) - a(i-1,k))*(dx1*dz1)
!
      v_2 = ((a(i,k+2) + a(i,k-2)) + (a(i+2,k) + a(i-2,k)) - 2.*(a(i+1,k+1) + a(i+1,k-1) + a(i-1,k+1) + a(i-1,k-1)))&
              *(dx1**2 * dz1**2)
!
      v_3 = beta*(a(i+1,k+1) - a(i+1,k-1) - a(i-1,k+1) + a(i-1,k-1))*(dx1*dz1)
!      
      v_4 = 4*(a(i+1,k+1) - 2.*a(i+1,k) + a(i+1,k-1) - 2.*a(i,k+1) - 2.*a(i,k-1) + a(i-1,k+1) - 2.*a(i-1,k) + a(i-1,k-1))&
              *(dx1**2 * dz1**2)

      v = v_1 + v_2 + v_3 + v_4

!      u = 2*alpha*(dx1**2 - dz1**2) + 6.*(dx1**4 + dz1**4) + &
!              4.*beta*dx1*dz1 + 8.*dx1**2*dz1**2
!      v_1 = (a(i,k+1) + a(i,k-1))*dz1**2 - (a(i+1,k) + a(i-1,k)) * dx1**2
!      v_2 = ( a(i,k+2) - 4.*a(i,k+1) - 4.*a(i,k-1) + a(i,k-2))   * dz1**4
!      v_3 = ( a(i+2,k) - 4.*a(i+1,k) - 4.*a(i-1,k) + a(i-2,k))   * dx1**4
!      v_4 = (-a(i-1,k) -    a(i,k-1) +    a(i-1,k-1)) * dx1*dz1
!!
!      v51 = a(i+1,k+1) - 2.*a(i,k+1) + a(i-1,k+1)
!      v52 = a(i+1,k  )               + a(i-1,k  )
!      v53 = a(i+1,k-1) - 2.*a(i,k-1) + a(i-1,k-1)
!      v_5 = (v51 - 2.*v52 + v53) * dx1**2*dz1**2
!!
!      v = alpha*v_1 + v_2 + v_3 + 4.*beta*v_4 + 2.*v_5
!
    endsubroutine solve_loworder
!***********************************************************************      
    subroutine calc_viscosity(f,eta)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,mz), intent(out) :: eta
!
!  Viscosities normalized by eta_0
!
      select case (iconv_viscosity)
!
      case ('constant')
        eta = eta_0
!
      case ('Netwonian') 
        eta = eta_0*exp(Avisc * (TT_melt/f(:,mpoint,:,iTT) - 1.))
!
      case ('Blankenbach-variable')
         do n=1,mz
            eta(:,n) = eta_0*exp(-Bvisc * f(:,mpoint,n,iTT)*deltaT1 + &
                                  Cvisc * (1-z(n)*Lz1)*Lz1 )
         enddo
!
      case default  
        write(unit=errormsg,fmt=*) &
             'calc_viscosity: No such value for iconv_viscosity: ', &
             trim(iconv_viscosity)
        call fatal_error('calc_viscosity',errormsg)
      endselect
!     
    endsubroutine calc_viscosity
!***********************************************************************
    subroutine f_function(a,aout,i,n)
!
!  Calculate F = dz2 - dx2, 6th order derivatives. 
!
      real, dimension(mx,mz), intent(in) :: a
      integer, intent(in) :: i,n
      real :: df2x,df2z,aout,fac
!
      fac=(1./180)*dx_1(i)**2
      df2x=fac*(-490.0*a(i,n) &
                +270.0*(a(i+1,n)+a(i-1,n)) &
                - 27.0*(a(i+2,n)+a(i-2,n)) &
                +  2.0*(a(i+3,n)+a(i-3,n)))

      fac=(1./180)*dz_1(n)**2
      df2z=fac*(-490.0*a(i,n) &
                +270.0*(a(i,n+1)+a(i,n-1)) &
                - 27.0*(a(i,n+2)+a(i,n-2)) &
                +  2.0*(a(i,n+3)+a(i,n-3)))

      aout = df2z-df2x
!
    endsubroutine f_function
!***********************************************************************
    subroutine g_function(a,aout,i,n)
!
!   Calculate G = dzdx, 6th order derivatives. 
!
      real, dimension (mx,mz), intent(in) :: a
      integer, intent(in) :: i,n
      real :: aout,fac
!
      fac=(1./60.**2)*dz_1(n)*dx_1(i)
      aout=fac*( &
              45.*((45.*(a(i+1,n+1)-a(i-1,n+1))  &
                    -9.*(a(i+2,n+1)-a(i-2,n+1))  &
                       +(a(i+3,n+1)-a(i-3,n+1))) &
                  -(45.*(a(i+1,n-1)-a(i-1,n-1))  &
                    -9.*(a(i+2,n-1)-a(i-2,n-1))  &
                       +(a(i+3,n-1)-a(i-3,n-1))))&
              -9.*((45.*(a(i+1,n+2)-a(i-1,n+2))  &
                    -9.*(a(i+2,n+2)-a(i-2,n+2))  &
                       +(a(i+3,n+2)-a(i-3,n+2))) &
                  -(45.*(a(i+1,n-2)-a(i-1,n-2))  &
                    -9.*(a(i+2,n-2)-a(i-2,n-2))  &
                       +(a(i+3,n-2)-a(i-3,n-2))))&
                 +((45.*(a(i+1,n+3)-a(i-1,n+3))  &
                    -9.*(a(i+2,n+3)-a(i-2,n+3))  &
                       +(a(i+3,n+3)-a(i-3,n+3))) &
                  -(45.*(a(i+1,n-3)-a(i-1,n-3))  &
                    -9.*(a(i+2,n-3)-a(i-2,n-3))  &
                       +(a(i+3,n-3)-a(i-3,n-3))))&
                   )
!
    endsubroutine g_function
!***********************************************************************
    subroutine update_bounds_psi(a)
!
      real, dimension(:,:) :: a
      integer :: i,nn1,nn2,ll1,ll2,mxa,mza,ll1i,ll2i
!
      mxa=size(a,1); mza=size(a,2)
      ll1=nghost+1; ll2=mxa-nghost
      nn1=nghost+1; nn2=mza-nghost
!
!  Set boundary of psi - vertical, zero
!
      a(:,nn1)=0.
      a(:,nn2)=0.
!
!  Zero also the second derivative
!
      do i=1,nghost
        a(:,nn1-i) = 2*a(:,nn1) - a(:,nn1+i)
      enddo
      do i=1,nghost
        a(:,nn2+i) = 2*a(:,nn2) - a(:,nn2-i)
      enddo
!
      if (lperi(1)) then 
!
!  Periodic in the lateral 
!
         ll1i=ll1+nghost-1;ll2i=ll2-nghost+1
         a(1   :ll1-1,:) = a(ll2i:ll2,:)
         a(ll2+1:mxa ,:) = a(ll1:ll1i,:)
      else
         a(ll1,:)=0.
         a(ll2,:)=0.
!
         do i=1,nghost
            a(ll1-i,:) = 2*a(ll1,:) - a(ll1+i,:)
         enddo
         do i=1,nghost
            a(ll2+i,:) = 2*a(ll2,:) - a(ll2-i,:)
         enddo
      endif
!
    endsubroutine update_bounds_psi
!***********************************************************************
    subroutine special_calc_energy(f,df,p)
!
      use Diagnostics, only: max_mn_name,sum_mn_name,integrate_mn_name
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: dTdz1,dTdz2,dTdz3,dTdz4,nusselt_num,nusselt_den,TTmin_cline,TTmax_cline
      real, dimension (nx) :: devsigzz1,devsigzz2,devsigzz3,devsigzz4
      real, dimension (nx) :: diffus_special,advec_special
      type (pencil_case) :: p
!      
!  Advection
!
      if (ltemperature_advection) then
         df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) - q%ugTT
         if (lsplit_temperature) then
            df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) - q%uu(:,3)*gTT_conductive
         endif
      endif
!
!  Conduction (diffusion)
!
      if (ltemperature_diffusion) & 
           df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + kappa*p%del2TT
!
!  Tidal heating 
!
      if (ltidal_heating) &
           df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + q%qtidal
!
      if (lfirst.and.ldt) then 
        advec_special=sum(abs(q%uu)*dline_1,2)
        maxadvec=maxadvec+advec_special
!
        diffus_special = kappa*dxyz_2
        maxdiffus=max(maxdiffus,diffus_special)
!
        if (headtt.or.ldebug) then
          print*,'special_calc_energy: max(advec_special)  =',maxval(advec_special)
          print*,'special_calc_energy: max(diffus_special) =',maxval(diffus_special)
        endif
      endif
!
      if (ldiagnos) then 
         if (idiag_uqxmin/=0) call max_mn_name(-q%uu(:,1)   ,idiag_uqxmin,lneg=.true.)
         if (idiag_uqxmax/=0) call max_mn_name( q%uu(:,1)   ,idiag_uqxmax)
         if (idiag_uqxm/=0)   call sum_mn_name( q%uu(:,1)   ,idiag_uqxm)
         if (idiag_uqx2m/=0)  call sum_mn_name( q%uu(:,1)**2,idiag_uqx2m)              
         if (idiag_uqxrms/=0) call sum_mn_name( q%uu(:,1)**2,idiag_uqxrms,lsqrt=.true.)
!
         if (idiag_uqzmin/=0) call max_mn_name(-q%uu(:,3)   ,idiag_uqzmin,lneg=.true.)
         if (idiag_uqzmax/=0) call max_mn_name( q%uu(:,3)   ,idiag_uqzmax)
         if (idiag_uqzm/=0)   call sum_mn_name( q%uu(:,3)   ,idiag_uqzm)
         if (idiag_uqz2m/=0)  call sum_mn_name( q%uu(:,3)**2,idiag_uqz2m)
         if (idiag_uqzrms/=0) call sum_mn_name( q%uu(:,3)**2,idiag_uqzrms,lsqrt=.true.)
!
         if (idiag_uq2m/=0)   call sum_mn_name(q%u2,idiag_uq2m)
         if (idiag_uqrms/=0)  call sum_mn_name(q%u2,idiag_uqrms,lsqrt=.true.)
         if (idiag_uqmax/=0)  call max_mn_name(q%u2,idiag_uqmax,lsqrt=.true.)
!
         if (idiag_qtidalmin/=0) call max_mn_name(-q%qtidal,idiag_qtidalmin,lneg=.true.)
         if (idiag_qtidalmax/=0) call max_mn_name( q%qtidal,idiag_qtidalmax)
         if (idiag_qtidalm/=0)   call sum_mn_name( q%qtidal,idiag_qtidalm)
         if (idiag_icount/=0)    call max_mn_name(0*x(l1:l2)+icount_save,idiag_icount)
         if (idiag_residual/=0)  call max_mn_name(0*x(l1:l2)+residual_save,idiag_residual)
!
!  Calculate for benchmark diagnostic the (negative of the) temperature gradient at the
!  four corners of the grid, labeled as below: 
!         
! 1 .______. 2
!   |      |
!   |      |
!   .______. 
! 4          3         
!         
!  I.e., 1 is (x,z)=( 0,Lz), top     left corner, above upwelling
!        2 is (x,z)=(Lx,Lz), top    right corner, above downwelling
!        3 is (x,z)=(Lx, 0), bottom right corner, below downwelling 
!        4 is (x,z)=( 0, 0), bottom  left corner, below upwelling
!
         if (idiag_dTdz1/=0) then
           if (n == n2) then
             dTdz1 = -p%gTT(1,3) - gTT_conductive
           else
             dTdz1 = -impossible
           endif
           call max_mn_name(dTdz1,idiag_dTdz1)
         endif
!
         if (idiag_dTdz2/=0) then
           if (n == n2) then
             dTdz2 = -p%gTT(nx,3) - gTT_conductive
           else
             dTdz2 = -impossible
           endif
           call max_mn_name(dTdz2,idiag_dTdz2)
         endif
!
         if (idiag_dTdz3/=0) then
           if (n == n1) then
             dTdz3 = -p%gTT(nx,3) - gTT_conductive
           else
             dTdz3 = -impossible
           endif
           call max_mn_name(dTdz3,idiag_dTdz3)
         endif
!
         if (idiag_dTdz4/=0) then
           if (n == n1) then
             dTdz4 = -p%gTT(1,3) - gTT_conductive
           else
             dTdz4 = -impossible
           endif
           call max_mn_name(dTdz4,idiag_dTdz4)
         endif
!
         if (idiag_nusselt_num/=0) then
           if (n == n2) then
             nusselt_num=-p%gTT(:,3) - gTT_conductive
           else
             nusselt_num=0.
           endif
           call integrate_mn_name(nusselt_num,idiag_nusselt_num)
         endif
!         
         if (idiag_nusselt_den/=0) then
           if (n == n1) then
             nusselt_den=p%TT + Tbar(n)
           else
             nusselt_den=0.
           endif
           call integrate_mn_name(nusselt_den,idiag_nusselt_den)
         endif
!         
!  temperature extreme at center line
!  min(abs(dTdz)) and z le 0.5
!
!  min(abs(dTdz)) and z ge 0.5
!
!  record position as well

         if (idiag_TTmax_cline/=0) then
            if (z(n) .ge. 0.5*Lxyz(3)) then
               TTmax_cline=p%TT(nx/2) + Tbar(n)
            else
               TTmax_cline=-impossible
            endif
            call max_mn_name(TTmax_cline,idiag_TTmax_cline)
         endif

         if (idiag_TTmin_cline/=0) then
            if (z(n) .lt. 0.5*Lxyz(3)) then
               TTmin_cline=p%TT(nx/2) + Tbar(n)
            else
               TTmin_cline=impossible
            endif
            call max_mn_name(-TTmin_cline,idiag_TTmin_cline,lneg=.true.)
         endif
!         
         if (idiag_devsigzz1/=0) then
            if (n == n2) then
               devsigzz1 = q%devsigzz(1)
            else
               devsigzz1 = -impossible
            endif
            call max_mn_name(devsigzz1,idiag_devsigzz1)
         endif
!         
         if (idiag_devsigzz2/=0) then
            if (n == n2) then
               devsigzz2 = q%devsigzz(nx)
            else
               devsigzz2 = -impossible
            endif
            call max_mn_name(devsigzz2,idiag_devsigzz2)
         endif
!
         if (idiag_devsigzz3/=0) then
            if (n == n1) then
               devsigzz3 = q%devsigzz(nx)
            else
               devsigzz3 = -impossible
            endif
            call max_mn_name(devsigzz3,idiag_devsigzz3)
         endif
!
         if (idiag_devsigzz4/=0) then
            if (n == n1) then
               devsigzz4 = q%devsigzz(1)
            else
               devsigzz4 = -impossible
            endif
            call max_mn_name(devsigzz4,idiag_devsigzz4)
         endif         
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_energy
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
      use Diagnostics, only: parse_name
      use FArrayManager, only: farray_index_append
!
      logical :: lreset
      logical, optional :: lwrite
!
      logical :: lwr
      integer :: iname, inamex
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_uqxmin=0; idiag_uqxmax=0; idiag_uqxrms=0; idiag_uqxm=0; idiag_uqx2m=0
        idiag_uqzmin=0; idiag_uqzmax=0; idiag_uqzrms=0; idiag_uqzm=0; idiag_uqz2m=0
        idiag_uq2m=0; idiag_uqrms=0; idiag_uqmax=0; idiag_qtidalmin=0
        idiag_qtidalmax=0; idiag_qtidalm=0
        idiag_dTdz1=0; idiag_dTdz2=0; idiag_dTdz3=0; idiag_dTdz4=0
        idiag_nusselt_num=0; idiag_nusselt_den=0
        idiag_TTmax_cline=0; idiag_TTmin_cline=0
        idiag_devsigzz1=0; idiag_devsigzz2=0; idiag_devsigzz3=0; idiag_devsigzz4=0
        idiag_icount=0; idiag_residual=0
        cformv=''
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'uqxmin',idiag_uqxmin)
        call parse_name(iname,cname(iname),cform(iname),'uqxmax',idiag_uqxmax)
        call parse_name(iname,cname(iname),cform(iname),'uqxrms',idiag_uqxrms)
        call parse_name(iname,cname(iname),cform(iname),'uqxm',idiag_uqxm)
        call parse_name(iname,cname(iname),cform(iname),'uqx2m',idiag_uqx2m)
!
        call parse_name(iname,cname(iname),cform(iname),'uqzmin',idiag_uqzmin)
        call parse_name(iname,cname(iname),cform(iname),'uqzmax',idiag_uqzmax)
        call parse_name(iname,cname(iname),cform(iname),'uqzrms',idiag_uqzrms)
        call parse_name(iname,cname(iname),cform(iname),'uqzm',idiag_uqzm)
        call parse_name(iname,cname(iname),cform(iname),'uqz2m',idiag_uqz2m)
!
        call parse_name(iname,cname(iname),cform(iname),'uq2m',idiag_uq2m)
        call parse_name(iname,cname(iname),cform(iname),'uqrms',idiag_uqrms)
        call parse_name(iname,cname(iname),cform(iname),'uqmax',idiag_uqmax)
!
        call parse_name(iname,cname(iname),cform(iname),'qtidalmin',idiag_qtidalmin)
        call parse_name(iname,cname(iname),cform(iname),'qtidalmax',idiag_qtidalmax)
        call parse_name(iname,cname(iname),cform(iname),'qtidalm',idiag_qtidalm)
!
        call parse_name(iname,cname(iname),cform(iname),'dTdz1',idiag_dTdz1)
        call parse_name(iname,cname(iname),cform(iname),'dTdz2',idiag_dTdz2)
        call parse_name(iname,cname(iname),cform(iname),'dTdz3',idiag_dTdz3)
        call parse_name(iname,cname(iname),cform(iname),'dTdz4',idiag_dTdz4)
!
        call parse_name(iname,cname(iname),cform(iname),'nusselt_num',idiag_nusselt_num)
        call parse_name(iname,cname(iname),cform(iname),'nusselt_den',idiag_nusselt_den)  
!
        call parse_name(iname,cname(iname),cform(iname),'TTmax_cline',idiag_TTmax_cline)
        call parse_name(iname,cname(iname),cform(iname),'TTmin_cline',idiag_TTmin_cline)
!
        call parse_name(iname,cname(iname),cform(iname),'devsigzz1',idiag_devsigzz1)
        call parse_name(iname,cname(iname),cform(iname),'devsigzz2',idiag_devsigzz2)
        call parse_name(iname,cname(iname),cform(iname),'devsigzz3',idiag_devsigzz3)
        call parse_name(iname,cname(iname),cform(iname),'devsigzz4',idiag_devsigzz4)
!        
        call parse_name(iname,cname(iname),cform(iname),'icount',idiag_icount)
        call parse_name(iname,cname(iname),cform(iname),'residual',idiag_residual)
     enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        !call parse_name(inamex,cnamex(inamex),cformx(inamex),'sigmamx', &
        !    idiag_sigmamx)
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='psi') cformv='DEFINED'
      endif
!
      if (lwr) then
        call farray_index_append('uqxmin',idiag_uqxmin)
        call farray_index_append('uqxmax',idiag_uqxmax)
        call farray_index_append('uqxrms',idiag_uqxrms)
        call farray_index_append('uqxm',idiag_uqxm)
        call farray_index_append('uqx2m',idiag_uqx2m)
!
        call farray_index_append('uqzmin',idiag_uqzmin)
        call farray_index_append('uqzmax',idiag_uqzmax)
        call farray_index_append('uqzrms',idiag_uqzrms)
        call farray_index_append('uqzm',idiag_uqzm)
        call farray_index_append('uqz2m',idiag_uqz2m)
!
        call farray_index_append('uq2m',idiag_uq2m)
        call farray_index_append('uqrms',idiag_uqrms)
        call farray_index_append('uqmax',idiag_uqmax)
!
        call farray_index_append('qtidalmin',idiag_qtidalmin)
        call farray_index_append('qtidalmax',idiag_qtidalmax)
        call farray_index_append('qtidalm',idiag_qtidalm)
!
        call farray_index_append('dTdz1',idiag_dTdz1)
        call farray_index_append('dTdz2',idiag_dTdz2)
        call farray_index_append('dTdz3',idiag_dTdz3)
        call farray_index_append('dTdz4',idiag_dTdz4)
!
        call farray_index_append('nusselt_num',idiag_nusselt_num)
        call farray_index_append('nusselt_den',idiag_nusselt_den)
!
        call farray_index_append('TTmax_cline',idiag_TTmax_cline)
        call farray_index_append('TTmin_cline',idiag_TTmin_cline)
!
        call farray_index_append('devsigzz1',idiag_devsigzz1)
        call farray_index_append('devsigzz2',idiag_devsigzz2)
        call farray_index_append('devsigzz3',idiag_devsigzz3)
        call farray_index_append('devsigzz4',idiag_devsigzz4)
!
        call farray_index_append('icount',idiag_icount)
        call farray_index_append('residual',idiag_residual)
     endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of Special variables.
!
!  26-jun-06/tony: dummy
!
      use Slices_methods, only: assign_slices_scal
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Streamfunction
!
        case('psi'); call assign_slices_scal(slices,f,ipsi)

      endselect
!
    endsubroutine get_slices_special
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!  Stuff for multigrid
!
!***********************************************************************
!***********************************************************************
!***********************************************************************    
    subroutine multigrid(f,psi,rhs,alpha_factor,beta_factor)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,mz) :: psi
      real, dimension(nx,nz) :: rhs
      real, dimension(nx,nz) :: alpha_factor,beta_factor
!
      call full_multigrid(f,psi,rhs,alpha_factor,beta_factor,n_vcycles)
!
    endsubroutine multigrid
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!  Stuff for multigrid, adapted from numerical recipes. 
!
!***********************************************************************
!***********************************************************************
!***********************************************************************    
    subroutine full_multigrid(f,psi,r,alp,bet,ncycle)
!
!  Full Multigrid Algorithm for solution of linear elliptic equations.
!  On input u contains the right-hand side ρ in an N × N array,
!  while on output it returns the solution. The dimension N is related
!  to the number of grid levels used in the solution, ng below, by N = 2**ng+1.
!  ncycle is the number of V-cycles to be used at each level.
!
!  13-jan-16/wlad: from numerical recipes.
!
      implicit none
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,mz), intent(inout) :: psi
      real, dimension(nx,nz), intent(in) :: r
      real, dimension(nx,nz), intent(in) :: alp,bet
      integer, intent(in) :: ncycle
      integer :: j,jcycle,ng,ngrid
      integer :: mvx,mvz,nvx,nvz,nmin,nxc,nzc,mxc,mzc
      integer:: ll1,ll2,ll1_,ll2_,nn1,nn2,nn1_,nn2_
!
      real :: residuals
      integer :: icount
      integer :: nux,nuz
!
!  Define a type so we can have an array of pointers to arrays of grid variables.
!
      type ptr2d
         real, pointer :: rhs(:,:)
         real, pointer :: alpha(:,:)
         real, pointer :: beta(:,:)
      endtype ptr2d
      type(ptr2d), allocatable :: grid(:)
      real, dimension(:,:), pointer :: uj,uj_1
      real, dimension(:,:), allocatable :: uj_old
!
      mvx=size(psi,1); mvz=size(psi,2)
      nvx=size(r  ,1); nvz=size(r  ,2)
      nmin=min(nvx,nvz)

      ngrid=nint(log(nmin-1.0)/log(2.0))
      if (nmin /= 2**ngrid+1) call fatal_error("full_multigrid",&
           "nxgrid-1 or nzgrid-1 must be a power of 2 in full_multigrid")
!
      allocate(grid(ngrid))
!
      nxc=nvx; nzc=nvz 
      ng=ngrid
!
!  Allocate storage for r.h.s. on grid ng, and fill it with the input r.h.s.
!
      allocate(grid(ng)%rhs(nxc,nzc))
      grid(ng)%rhs=r

      allocate(grid(ng)%alpha(nxc,nzc))
      grid(ng)%alpha=alp

      allocate(grid(ng)%beta(nxc,nzc))
      grid(ng)%beta=bet
      
!
!  Similarly allocate storage and fill r.h.s. on all coarse grids by restricting from finer grids.
!
      do
         if ((nxc <= 3).or.(nzc <= 3)) exit
         nxc=nxc/2+1
         nzc=nzc/2+1
         ng=ng-1
!         
         allocate(grid(ng)%rhs(nxc,nzc))
         grid(ng)%rhs=restrict(grid(ng+1)%rhs)
!
         allocate(grid(ng)%alpha(nxc,nzc))
         grid(ng)%alpha=restrict(grid(ng+1)%alpha)
!         
         allocate(grid(ng)%beta(nxc,nzc))
         grid(ng)%beta=restrict(grid(ng+1)%beta)
      enddo
!
!  Initial solution on coarsest grid.
!
      !nn=3
      mxc=nxc+2*nghost     
      mzc=nzc+2*nghost
      allocate(uj(mxc,mzc))
      uj=0.0
      call solve_coarsest(uj,grid(1)%rhs,grid(1)%alpha,grid(1)%beta)
!
!  Nested iteration loop. 
!
      do j=2,ngrid
         ll1_=nghost+1; ll2_=mxc-nghost
         nn1_=nghost+1; nn2_=mzc-nghost
!         
         nxc=2*nxc-1; nzc=2*nzc-1
         mxc=nxc+2*nghost; mzc=nzc+2*nghost
!
         ll1=nghost+1; ll2=mxc-nghost
         nn1=nghost+1; nn2=mzc-nghost
!
         uj_1=>uj
         allocate(uj(mxc,mzc))
         uj=0.0
!
!  Prolongate (interpolate) from grid j-1 to next finer grid j. 
!
         uj(ll1:ll2,nn1:nn2)=prolongate(uj_1(ll1_:ll2_,nn1_:nn2_))
         deallocate(uj_1)
!
!  V-cycle loop. Iterate until convergence. 
!
         residuals=1e33
         icount=0
         do while (residuals > tolerance)
            nux=size(uj,1); nuz=size(uj,2)
            allocate(uj_old(nux,nuz))
            uj_old=uj
            !do jcycle=1,ncycle
            call vcycle(j,uj,grid(j)%rhs,grid(j)%alpha,grid(j)%beta)
            !enddo
            if (j==ngrid.and.lsave_residual_grid) then
              do m=1,my
                f(:,m,:,iresidual) = uj-uj_old
              enddo
            endif
            residuals=sqrt(sum((uj-uj_old)**2)/(nux*nuz))
            icount=icount+1
            if (lprint_residual) print*,icount,residuals,j
            if (lsave_residual) write(9,*) icount,residuals,j
            deallocate(uj_old)
            if (icount >= maxit) goto 100
            !call fatal_error("full multigrid",&
            !     "reached limit of iterations: did not converge.")
         enddo
100      continue 
      enddo
      if (lsave_residual) close(9)
      icount_save=icount
      residual_save=residuals
!      
!  Return solution in psi.
!
      psi=uj
!
!  Deallocate all grids
!
      deallocate(uj)
      do j=1,ngrid
         deallocate(grid(j)%rhs)
         deallocate(grid(j)%alpha)
         deallocate(grid(j)%beta)
!         
      enddo
      deallocate(grid)
!
      contains
!--------------------------------------------------------------------
      recursive subroutine vcycle(jgrid,u,rrhs,alpp,bett)
!
!  Recursive multigrid iteration. On input, jgrid is the current level,
!  u is the current value of the solution, and rhs is the right-hand
!  side. On output u contains the improved solution at the current level.
!  Parameters: npre and npost are the number of relaxation sweeps before
!  and after the coarse-grid correction is computed.
!
!  13-jan-16/wlad: from numerical recipes.
!        
        implicit none
        integer, intent(in) :: jgrid
        real, dimension(:,:), intent(inout) :: u
        real, dimension(:,:), intent(in) :: rrhs
!
        real, dimension(:,:), intent(in) :: alpp,bett
!        
        integer :: jpost,jpre
        integer :: mmvx,mmvz,mmux,mmuz
        integer :: lv1,lv2,lu1,lu2,nv1,nv2,nu1,nu2
!        
        integer :: icycle
        real, dimension((size(rrhs,1)+1)/2,(size(rrhs,2)+1)/2) :: res
        real, dimension((size(rrhs,1)+1)/2,(size(rrhs,2)+1)/2) :: alpp_,bett_
        real, dimension((size(rrhs,1)+1)/2+2*nghost,(size(rrhs,2)+1)/2+2*nghost) :: v
!        
        if (jgrid == 1) then
!
!  Bottom of V: Solve on coarsest grid.
!
           call solve_coarsest(u,rrhs,alpp,bett)
        else
!
!  On downward stroke of the V.
!
!  Pre-smoothing.              
!
           do jpre=1,npre
              call relaxation(u,rrhs,alpp,bett)
           enddo
!
!  Restriction of the residual is the next r.h.s.
!
           res=restrict(resid(u,rrhs,alpp,bett))
!
           alpp_=restrict(alpp)
           bett_=restrict(bett)
!
!  Zero for initial guess in next relaxation.
!
           v=0.0
!
!  Recursive call for the coarse grid correction.
!
           do icycle=1,merge(igamma,1,jgrid/=2)
              call vcycle(jgrid-1,v,res,alpp_,bett_)
           enddo
!
!  On upward stroke of V.
!
           mmvx=size(v,1); mmvz=size(v,2)
           lv1=nghost+1; lv2=mmvx-nghost
           nv1=nghost+1; nv2=mmvz-nghost
!           
           mmux=size(u,1); mmuz=size(u,2)
           lu1=nghost+1; lu2=mmux-nghost
           nu1=nghost+1; nu2=mmuz-nghost
!
           u(lu1:lu2,nu1:nu2)=u(lu1:lu2,nu1:nu2)+prolongate(v(lv1:lv2,nv1:nv2))
!
!  Post-smoothing.
!
           do jpost=1,npost
              call relaxation(u,rrhs,alpp,bett)
           enddo
        endif
!
      endsubroutine vcycle
!--------------------------------------------------------------------
    endsubroutine full_multigrid
!********************************************************************                
    function restrict(uf)
!
!  Half-weighting restriction. If Nc is the coarse-grid dimension, the
!  fine-grid solution is input in the (2Nc − 1) × (2Nc − 1) array uf,
!  the coarse-grid solution is returned in the Nc × Nc array restrict.
!
!  13-jan-16/wlad: from numerical recipes.
!
      implicit none
      real, dimension(:,:), intent(in) :: uf
      real, dimension((size(uf,1)+1)/2,(size(uf,2)+1)/2) :: restrict
      integer :: ncx,ncz,nfx,nfz
!
      nfx=size(uf,1); nfz=size(uf,2)
      ncx=(nfx+1)/2; ncz=(nfz+1)/2
!
      !if (nfx == nfz) then 
!
!  Interior points
!      
         restrict(2:ncx-1,2:ncz-1) =   0.5 *uf(3:nfx-2:2,3:nfz-2:2) + &
                                     0.125*(uf(4:nfx-1:2,3:nfz-2:2) + uf(2:nfx-3:2,3:nfz-2:2)+ &
                                            uf(3:nfx-2:2,4:nfz-1:2) + uf(3:nfx-2:2,2:nfz-3:2))
!
!  Boundary points      
!
         restrict(1:ncx,1    ) = uf(1:nfx:2,1      )
         restrict(1:ncx,ncz  ) = uf(1:nfx:2,nfz    )
         restrict(1    ,1:ncz) = uf(1      ,1:nfz:2)
         restrict(ncx  ,1:ncz) = uf(nfx    ,1:nfz:2)
      !else
      !   restrict=1.
      !   print*,'restricted, go on'
      !endif
!      
    endfunction restrict
!********************************************************************
    function prolongate(uc)
!      
!  Coarse-to-fine prolongation by bilinear interpolation. If Nf is the
!  fine-grid dimension and Nc the coarse-grid dimension, then
!  Nf = 2Nc − 1. The coarse-grid solution is input as uc, the fine-grid
!  solution is returned in prolongate.
!
!  13-jan-16/wlad: from numerical recipes.
!
      implicit none
      real, dimension(:,:), intent(in) :: uc
      real, dimension(2*size(uc,1)-1,2*size(uc,2)-1) :: prolongate
      integer :: ncx,ncz,nfx,nfz
!
      ncx=size(uc,1); ncz=size(uc,2)
      nfx=2*ncx-1; nfz=2*ncz-1

      !if (ncx == ncz) then 
!
! Do elements that are copies.
!
         prolongate(1:nfx:2,1:nfz:2)=uc(1:ncx,1:ncz)
!
! Do odd-numbered columns, interpolating vertically.
!
         prolongate(2:nfx-1:2,1:nfz:2) = 0.5*(prolongate(3:nfx:2,1:nfz:2) + prolongate(1:nfx-2:2,1:nfz:2))
!
! Do even-numbered columns, interpolating horizontally.
!
         prolongate(1:nfx,2:nfz-1:2) = 0.5*(prolongate(1:nfx,3:nfz:2) + prolongate(1:nfx,1:nfz-2:2))
      !else
      !   prolongate=1.
      !   print*,'prolongated, go on'
      !endif
!
    endfunction prolongate
!********************************************************************
    subroutine solve_coarsest(u,rhs,alp,bet)
!
!  Solution of the model problem on the coarsest grid. 
!  Input in rhs(1:3,1:3) and the solution is returned in u(1:3,1:3). 
!
!  13-jan-16/wlad: from numerical recipes.
!
      implicit none
      real, dimension(:,:), intent(inout) :: u
      real, dimension(:,:), intent(in) :: rhs

      real, dimension(:,:), intent(in) :: alp,bet

      real, dimension(size(rhs,1),size(rhs,2)) :: u_old
      real :: h,hx,hz
      real :: alpha=impossible,beta=impossible
!
      real :: res,ufactor,vterm
      integer :: icount,iu
      integer :: nu,nux,nuz,nrx,nrz,lc1,lc2,nc1,nc2
!
      nux=size(u,1); nuz=size(u,2)
      nrx=size(rhs,1); nrz=size(rhs,2)
!
      u = 0.0
      hx = Lxyz(1)/(nrx-1)
      hz = Lxyz(3)/(nrz-1)
      if (hx==hz) then
        h = hx
      else
        print*,"hx=",hx,"hz=",hz
        call fatal_error("solve_coarsest","grid cells are not square")
      endif
!
      if (ldirect_solver) then
         if (nx_coarsest /= 3) call fatal_error("solve_coarsest",&
              "direct solver only for nx=3 in the coarsest grid")
         if (nux/=nuz) call fatal_error("solve_coarsest",&
              "direct solver only for square grids")
         nu=nux
         iu=(nu+1)/2
         if (lpoisson_test) then 
            u(iu,iu) = -h**2 * rhs(2,2)/4.0
         else
            alpha=alp(2,2)
            beta=bet(2,2)
            if (lsolver_highorder) then 
               call solve_highorder(u,alpha,beta,ufactor,vterm,iu,iu,h)
            else
               call solve_loworder(u,alpha,beta,ufactor,vterm,iu,iu,h)
            endif
            u(iu,iu) =  rhs(2,2)/ufactor 
         endif
      else
!
!  No multigrid, do SOR. 
!         
         lc1=nghost+1; lc2=nux-nghost
         nc1=nghost+1; nc2=nuz-nghost      
!      
         res=1e33
         icount=0
!
         do while (res > tolerance_coarsest)
            icount=icount+1
            u_old=u(lc1:lc2,nc1:nc2)
            !call successive_over_relaxation(u,rhs,alp,bet,h)
            
            call relaxation(u,rhs,alp,bet)
            res = sqrt(sum((u(lc1:lc2,nc1:nc2)-u_old)**2)/(nrx*nrz))
            if (ldebug) print*,'solve_coarsest:',icount,res,tolerance_coarsest
         enddo
!
      endif
!
    endsubroutine solve_coarsest
!********************************************************************
    subroutine relaxation(u,rhs,alp,bet)
!
!  Red-black Gauss-Seidel relaxation for model problem. The current
!  value of the solution u is updated, using the right-hand-side
!  function rhs. u and rhs are square arrays of the same odd dimension.
!
!  13-jan-16/wlad: from numerical recipes.
!
      implicit none
      real, dimension(:,:), intent(inout) :: u
      real, dimension(:,:), intent(in) :: rhs

      real, dimension(:,:), intent(in) :: alp,bet

      integer :: l,k
      real :: h,h2
!
      real :: cc,ufactor,vterm
      real :: alpha=impossible,beta=impossible
!
      integer :: mvx,mvz,nvx,nvz
      integer :: ic,nc,icc1,icc2
      real :: hx,hz
!
      mvx=size(u  ,1); mvz=size(u  ,2)
      nvx=size(rhs,1); nvz=size(rhs,2)
!
      hx=Lxyz(1)/(nvx-1)
      hz=Lxyz(3)/(nvz-1)
      if (hx==hz) then
        h = hx
      else
        print*,"hx=",hx,"hz=",hz
        call fatal_error("relaxation","grid cells are not square")
      endif

!
! Convection
!
      if (lperi(1)) then
         icc1=1
         icc2=nvx
      else
         icc1=2
         icc2=nvx-1
      endif
!         
      if (.not.lpoisson_test) then 
         call update_bounds_psi(u)
         alpha=0.0
         beta=0.0
         do ic=icc1,icc2
            l=ic+nghost
            do nc=2,nvz-1
              k=nc+nghost
               cc = rhs(ic,nc)
               alpha=alp(ic,nc)
               beta=bet(ic,nc)
               if (lsolver_highorder) then
                  call solve_highorder(u,alpha,beta,ufactor,vterm,l,k,h)
               else
                  call solve_loworder(u,alpha,beta,ufactor,vterm,l,k,h)
               endif
               u(l,k) = (cc-vterm)/ufactor
            enddo
         enddo
         call update_bounds_psi(u)
      else
!
! Poisson
!     
         h2=h**2
!
         do ic=2,nvx-1
            l=ic+nghost
            do nc=2,nvz-1
               k=nc+nghost
               u(l,k) = 0.25*(u(l+1,k) + u(l-1,k) + &
                              u(l,k+1) + u(l,k-1) - h2*rhs(ic,nc))
            enddo
         enddo
      endif
!
    endsubroutine relaxation
!********************************************************************
    function resid(u,rhs,alp,bet)
!
!  Returns minus the residual for the model problem. Input quantities are u and rhs,
!  while the residual is returned in resid. All three quantities are square arrays
!  with the same odd dimension.
!
!  13-jan-16/wlad: from numerical recipes.
!
      implicit none
      real, dimension(:,:), intent(in) :: u
      real, dimension(:,:), intent(in) :: rhs

      real, dimension(:,:), intent(in) :: alp,bet

      real, dimension(size(rhs,1),size(rhs,2)) :: resid
      integer :: l,k
      real :: h,h2i
!
      real :: ufactor,vterm,lhs
      real :: alpha=impossible,beta=impossible
!
      integer :: mvx,mvz,nvx,nvz
      integer :: ic,nc
      real :: hx,hz
!
      mvx=size(u  ,1); mvz=size(u  ,2)
      nvx=size(rhs,1); nvz=size(rhs,2)
!
      hx=Lxyz(1)/(nvx-1)
      hz=Lxyz(3)/(nvz-1)
      if (hx==hz) then
        h = hx
      else
        print*,"hx=",hx,"hz=",hz
        call fatal_error("resid","grid cells are not square")
      endif
!
!  Interior points.
!
      if (lpoisson_test) then 
         h2i=1.0/(h**2)
         do ic=2,nvx-1
            l=ic+nghost
            do nc=2,nvz-1
               k=nc+nghost
               resid(ic,nc) = -h2i*(u(l+1,k  ) + u(l-1,k  ) + &
                                    u(l  ,k+1) + u(l  ,k-1) - &
                                    4.0*u(l,k))+rhs(ic,nc)
            enddo
         enddo
      else
         do ic=2,nvx-1
            l=ic+nghost
            do nc=2,nvz-1
               k=nc+nghost
               alpha=alp(ic,nc)
               beta=bet(ic,nc)
               if (lsolver_highorder) then
                  call solve_highorder(u,alpha,beta,ufactor,vterm,l,k,h)
               else
                  call solve_loworder(u,alpha,beta,ufactor,vterm,l,k,h)
               endif
               lhs = ufactor*u(l,k) + vterm
               resid(ic,nc) = rhs(ic,nc) - lhs
            enddo
         enddo
      endif
!
!  Boundary points.
!      
      resid(1:nvx,1    )=0.0
      resid(1:nvx,nvz  )=0.0
      resid(1    ,1:nvz)=0.0
      resid(nvx  ,1:nvz)=0.0
!
    endfunction resid
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
