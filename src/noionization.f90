! $Id: noionization.f90,v 1.59 2003-08-26 16:40:36 mee Exp $

!  Dummy routine for noionization

module Ionization

  use Cparam
  use Cdata
  use Density

  implicit none

  interface thermodynamics              ! Overload the `thermodynamics' function
    module procedure thermodynamics_pencil   ! explicit f implicit m,n
    module procedure thermodynamics_point    ! explocit lnrho, ss
  end interface

  interface ionget
    module procedure ionget_pencil
    module procedure ionget_point
  end interface

  interface perturb_energy              ! Overload subroutine perturb_energy
    module procedure perturb_energy_pencil
    module procedure perturb_energy_point
  end interface

  ! secondary parameters calculated in initialize
  real :: TT_ion,TT_ion_,ss_ion,kappa0
  real :: lnrho_H,lnrho_e,lnrho_e_,lnrho_p,lnrho_He
  real :: yHmin,yHmax
  real :: xHe_term,yH_term,one_yH_term

  !  lionization initialized to .false.
  !  cannot currently be reset to .true. in namelist
  !  because the namelist is now not even read
  logical :: lionization=.false.,lionization_fixed=.false.
  real :: xHe=0.1

  ! input parameters
  integer :: dummy_ni 
  namelist /ionization_init_pars/ xHe

  ! run parameters
  namelist /ionization_run_pars/ xHe

  contains

!***********************************************************************
    subroutine register_ionization()
!
!  14-jun-03/axel: adapted from register_ionization
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_ionization called twice')
      first = .false.
!
      iyH = 0
      iTT = 0

      if ((ip<=8) .and. lroot) then
        print*, 'register_ionization: ionization nvar = ', nvar
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: noionization.f90,v 1.59 2003-08-26 16:40:36 mee Exp $")
!
!  Check we aren't registering too many auxiliary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('Register_ionization: naux > maux')
      endif
!
    endsubroutine register_ionization
!***********************************************************************
    subroutine initialize_ionization()
!
      if(ip==0) print*,'initialize_ionization: keeping compiler quiet'
!
    endsubroutine initialize_ionization
!*******************************************************************
    subroutine getmu(mu)
      real, intent(out) :: mu
      mu=1.+3.97153*xHe  
    endsubroutine getmu
!*******************************************************************
    subroutine rprint_ionization(lreset)
!
!  Writes iyH and iTT to index.pro file
!
!  14-jun-03/axel: adapted from rprint_radiation
!
      use Cdata
      use Sub
! 
      logical :: lreset
!
!  write column where which ionization variable is stored
!
      write(3,*) 'nname=',nname
      write(3,*) 'iyH=',iyH
      write(3,*) 'iTT=',iTT
!   
      if(ip==0) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_ionization
!***********************************************************************
    subroutine ioninit(f)
!   
      real, dimension (mx,my,mz,mvar+maux), intent(inout) :: f
!   
      if(ip==0) print*,f(1,1,1,1)  !(keep compiler quiet)
!   
    endsubroutine ioninit
!***********************************************************************
    subroutine ioncalc(f)
!
    real, dimension (mx,my,mz,mvar+maux) :: f
!
    if(ip==0) print*,f(1,1,1,1)  !(keep compiler quiet)
!
    endsubroutine ioncalc
!***********************************************************************
    subroutine perturb_energy_point(lnrho,EE,ss,TT)
      use Mpicomm, only: stop_it
      
      real, intent(in) :: lnrho,EE
      real, intent(out) :: ss,TT

      call stop_it("perturb_energy_point: NOT IMPLEMENTED IN NO IONIZATION")
      ss=0.
      TT=0.
      if (ip==0) print*,lnrho,EE
    end subroutine perturb_energy_point
!***********************************************************************
    subroutine perturb_energy_pencil(lnrho,EE,ss,TT)
      use Mpicomm, only: stop_it
      
      real, dimension(nx), intent(in) :: lnrho,EE
      real, dimension(nx), intent(out) :: ss,TT

      call stop_it("perturb_energy_pencil: NOT IMPLEMENTED IN NO IONIZATION")
      ss=0.
      TT=0.
      if (ip==0) print*,lnrho,EE
    end subroutine perturb_energy_pencil
!***********************************************************************
    subroutine getdensity(EE,TT,yH,rho)
      use Mpicomm, only: stop_it
      
      real, intent(in) :: EE,TT,yH
      real, intent(out) :: rho

      call stop_it("getdensity: NOT IMPLEMENTED IN NO IONIZATION")
      !rho = EE / ((1.5*(1+yH+xHe)*TT + yH*TT_ion) * ss_ion)
      if (ip==0) print*,EE,TT,yH
      rho=0.

    end subroutine getdensity
!***********************************************************************
    subroutine isothermal_density_ion(pot,tmp)
!
      real, dimension (nx), intent(in) :: pot
      real, dimension (nx), intent(out) :: tmp
!
      tmp=pot
!
    end subroutine isothermal_density_ion
!***********************************************************************
    subroutine ionget_pencil(f,yH,TT)
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(:), intent(inout) :: yH,TT
      real, dimension(size(yH)) :: lnrho,ss
!
      if (size(lnrho)==nx) lnrho=f(l1:l2,m,n,ilnrho)
      if (size(ss)==nx) ss=f(l1:l2,m,n,iss)
!
      if (size(lnrho)==mx) lnrho=f(:,m,n,ilnrho)
      if (size(ss)==mx) ss=f(:,m,n,iss)
!
      TT=exp(gamma*ss+gamma1*(lnrho-lnrho0)+alog(cs20/gamma1))
!
    endsubroutine ionget_pencil
!***********************************************************************
    subroutine ionget_point(lnrho,ss,yH,TT)
!
      use Cdata
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: yH,TT
!
      TT=exp(gamma*ss+gamma1*(lnrho-lnrho0)+alog(cs20/gamma1))
      yH=0.
!
    endsubroutine ionget_point
!***********************************************************************
    subroutine thermodynamics_pencil(lnrho,ss,yH,TT,cs2,cp1tilde,ee)
!
!  Calculate thermodynamical quantities, cs2, 1/T, and cp1tilde
!  cs2=(dp/drho)_s is the adiabatic sound speed
!  TT1=1/T is the inverse temperature
!  neutral gas: cp1tilde=kB_over_mp/cp=0.4 ("nabla_ad" maybe better name)
!  in general: cp1tilde=dlnPdS/dlnPdlnrho
!
!   2-feb-03/axel: simple example coded
!  15-jun-03/axel: made compatible with current ionization routine
!
      use Cdata
      use General
      use Sub
!
      real, dimension(nx), intent(in) :: lnrho,ss,yH,TT
      real, dimension(nx), optional :: cs2,cp1tilde,ee
!
      if (present(cs2)) cs2=gamma1*exp(gamma*ss+gamma1*(lnrho-lnrho0)+alog(cs20/gamma1))
      if (present(cp1tilde)) cp1tilde=1.
      if (present(ee)) ee=exp(gamma*ss+gamma1*(lnrho-lnrho0)+alog(cs20/gamma1))/gamma
!
      if (ip==0) print*,yH,TT
    endsubroutine thermodynamics_pencil
!***********************************************************************
    subroutine thermodynamics_point(lnrho,ss,yH,TT,cs2,cp1tilde,ee)
!
!  Calculate thermodynamical quantities, cs2, 1/T, and cp1tilde
!  cs2=(dp/drho)_s is the adiabatic sound speed
!  TT1=1/T is the inverse temperature
!  neutral gas: cp1tilde=kB_over_mp/cp=0.4 ("nabla_ad" maybe better name)
!  in general: cp1tilde=dlnPdS/dlnPdlnrho
!
!   2-feb-03/axel: simple example coded
!  15-jun-03/axel: made compatible with current ionization routine
!
      use Cdata
      use General
      use Sub
!
      real, intent(in) :: lnrho,ss,yH,TT
      real, optional :: cs2,cp1tilde,ee
!
      if (present(cs2)) cs2=gamma1*exp(gamma*ss+gamma1*(lnrho-lnrho0)+alog(cs20/gamma1))
      if (present(cp1tilde)) cp1tilde=1.
      if (present(ee)) ee=exp(gamma*ss+gamma1*(lnrho-lnrho0)+alog(cs20/gamma1))/gamma
!
      if (ip==0) print*,yH,TT
    endsubroutine thermodynamics_point
!***********************************************************************
    subroutine bc_ss_flux(f,topbot,hcond0,hcond1,Fbot, FbotKbot, chi, &
                lmultilayer,lcalc_heatcond_constchi)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  26-aug-2003/tony: distributed across ionization modules
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      real, intent(in) :: Fbot, FbotKbot, hcond0, hcond1, chi
      logical, intent(in) :: lmultilayer, lcalc_heatcond_constchi
      
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my) :: tmp_xy,cs2_xy,rho_xy
      integer :: i
      
!
      if(ldebug) print*,'bc_ss_flux: ENTER - cs20,cs0=',cs20,cs0
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('strange-bot')
        if(headtt) print*,'bc_ss_flux: hcond0,hcond1=',hcond0,hcond1
        if ((bcz1(ilnrho) /= "a2") .and. (bcz1(ilnrho) /= "a3"))&
             call stop_it("bc_ss_flux: Inconsistent boundary conditions 1.")
        tmp_xy = gamma1/cs20 & ! 1/T_0 (i.e. 1/T at boundary)
                 * exp(-gamma*f(:,:,n1,iss) &
                       - gamma1*(f(:,:,n1,ilnrho)-lnrho0))
        tmp_xy = Fbot/(hcond0*hcond1) * tmp_xy ! F_heat/(hcond T_0)
        do i=1,nghost
          f(:,:,n1-i,iss) = &
               (2*i*dz*tmp_xy &
                + 2*gamma1*(f(:,:,n1+i,ilnrho)-f(:,:,n1,ilnrho)) &
               )/gamma &
               + f(:,:,n1+i,iss)
        enddo
!
!  bottom boundary
!  ===============
!
      case('bot')
        if (lmultilayer) then
          if(headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0*hcond1
        else
          if(headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0
        endif
!       if(bcz1(ilnrho)/="a2") call stop_it("bc_ss_flux: bad lnrho bc")
!
!  calculate Fbot/(K*cs2)
!
        rho_xy=exp(f(:,:,n1,ilnrho))
        cs2_xy=cs20*exp(gamma1*(f(:,:,n1,ilnrho)-lnrho0)+gamma*f(:,:,n1,iss))
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy. 
!
        if(lcalc_heatcond_constchi) then
          tmp_xy=Fbot/(rho_xy*chi*cs2_xy)
        else
          tmp_xy=FbotKbot/cs2_xy
        endif
!
!  enforce ds/dz + gamma1/gamma*dlnrho/dz = - gamma1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+gamma1/gamma* &
              (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)+2*i*dz*tmp_xy)
        enddo
!
!  top boundary
!  ============
!
      case('top')
        if(headtt) print*,'bc_ss_flux: hcond0=',hcond0
        if ((bcz2(ilnrho) /= "a2") .and. (bcz2(ilnrho) /= "a3")) &
             call stop_it("bc_ss_flux: Inconsistent boundary conditions 2.")
        tmp_xy = gamma1/cs20 & ! 1/T_0 (i.e. 1/T at boundary)
                 * exp(-gamma*f(:,:,n2,iss) &
                       - gamma1*(f(:,:,n2,ilnrho)-lnrho0))
        tmp_xy = FbotKbot * tmp_xy ! F_heat/(hcond T_0)
        do i=1,nghost
          f(:,:,n2+i,iss) = &
               (-2*i*dz*tmp_xy &
                + 2*gamma1*(f(:,:,n2-i,ilnrho)-f(:,:,n2,ilnrho)) &
               )/gamma &
               + f(:,:,n2-i,iss)
        enddo
      case default
        print*,"bc_ss_flux: invalid argument"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_flux
!***********************************************************************
    subroutine bc_ss_temp_old(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  23-jun-2003/tony: implemented for lionization_fixed
!  26-aug-2003/tony: distributed across ionization modules
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my) :: tmp_xy
      integer :: i
!
      if(ldebug) print*,'bc_ss_temp_old: ENTER - cs20,cs0=',cs20,cs0
!
!  Do the `c2' boundary condition (fixed temperature/sound speed) for entropy.
!  This assumes that the density is already set (ie density must register
!  first!)
!  tmp_xy = s(x,y) on the boundary.
!  gamma*s/cp = [ln(cs2/cs20)-(gamma-1)ln(rho/rho0)]
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if ((bcz1(ilnrho) /= "a2") .and. (bcz1(ilnrho) /= "a3")) &
          call stop_it("bc_ss_temp_old: Inconsistent boundary conditions 3.")
        if (ldebug) print*, &
                'bc_ss_temp_old: set bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) &
              print*,'bc_ss_temp_old: cannot have cs2bot<=0'
        tmp_xy = (-gamma1*(f(:,:,n1,ilnrho)-lnrho0) &
             + alog(cs2bot/cs20)) / gamma
        f(:,:,n1,iss) = tmp_xy
        do i=1,nghost
          f(:,:,n1-i,iss) = 2*tmp_xy - f(:,:,n1+i,iss)
        enddo
 
!
!  top boundary
!
      case('top')
        if ((bcz1(ilnrho) /= "a2") .and. (bcz1(ilnrho) /= "a3")) &
          call stop_it("bc_ss_temp_old: Inconsistent boundary conditions 3.")
        if (ldebug) print*, &
                   'bc_ss_temp_old: set top temperature - cs2top=',cs2top
        if (cs2top<=0.) print*, &
                   'bc_ss_temp_old: cannot have cs2top<=0'
  !     if (bcz1(ilnrho) /= "a2") &
  !          call stop_it("BOUNDCONDS: Inconsistent boundary conditions 4.")
        tmp_xy = (-gamma1*(f(:,:,n2,ilnrho)-lnrho0) &
                 + alog(cs2top/cs20)) / gamma
        f(:,:,n2,iss) = tmp_xy
        do i=1,nghost
          f(:,:,n2+i,iss) = 2*tmp_xy - f(:,:,n2-i,iss)
        enddo 
      case default
        print*,"bc_ss_temp_old: invalid argument"
        call stop_it("")
      endselect
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
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: tmp
      integer :: i
!
      if(ldebug) print*,'bc_ss_temp_x: cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (ldebug) print*, &
                   'bc_ss_temp_x: set x bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp_x: cannot have cs2bot<=0'
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(l1,:,:,iss) = 0.5*tmp - gamma1/gamma*(f(l1,:,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(l1-i,:,:,iss) = -f(l1+i,:,:,iss) + tmp &
               - gamma1/gamma*(f(l1+i,:,:,ilnrho)+f(l1-i,:,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                       'bc_ss_temp_x: set x top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*, &
                       'bc_ss_temp_x: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(l2,:,:,iss) = 0.5*tmp - gamma1/gamma*(f(l2,:,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(l2+i,:,:,iss) = -f(l2-i,:,:,iss) + tmp &
               - gamma1/gamma*(f(l2-i,:,:,ilnrho)+f(l2+i,:,:,ilnrho)-2*lnrho0)
        enddo

      case default
        print*,"bc_ss_temp_x: invalid argument"
        call stop_it("")
      endselect
      

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
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: tmp
      integer :: i
!
      if(ldebug) print*,'bc_ss_temp_y: cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (ldebug) print*, &
                   'bc_ss_temp_y: set y bottom temperature - cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp_y: cannot have cs2bot<=0'
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(:,m1,:,iss) = 0.5*tmp - gamma1/gamma*(f(:,m1,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m1-i,:,iss) = -f(:,m1+i,:,iss) + tmp &
               - gamma1/gamma*(f(:,m1+i,:,ilnrho)+f(:,m1-i,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                     'bc_ss_temp_y: set y top temperature - cs2top=',cs2top
        if (cs2top<=0.) print*, &
                     'bc_ss_temp_y: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(:,m2,:,iss) = 0.5*tmp - gamma1/gamma*(f(:,m2,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m2+i,:,iss) = -f(:,m2-i,:,iss) + tmp &
               - gamma1/gamma*(f(:,m2-i,:,ilnrho)+f(:,m2+i,:,ilnrho)-2*lnrho0)
        enddo

      case default
        print*,"bc_ss_temp_y: invalid argument"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_temp_y
!***********************************************************************
    subroutine bc_ss_temp_z(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: tmp
      integer :: i
!
      if(ldebug) print*,'bc_ss_temp_z: cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (ldebug) print*, &
                   'bc_ss_temp_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp_z: cannot have cs2bot<=0'
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(:,:,n1,iss) = 0.5*tmp - gamma1/gamma*(f(:,:,n1,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,:,n1-i,iss) = -f(:,:,n1+i,iss) + tmp &
               - gamma1/gamma*(f(:,:,n1+i,ilnrho)+f(:,:,n1-i,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                     'bc_ss_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*,'bc_ss_temp_z: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(:,:,n2,iss) = 0.5*tmp - gamma1/gamma*(f(:,:,n2,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,:,n2+i,iss) = -f(:,:,n2-i,iss) + tmp &
               - gamma1/gamma*(f(:,:,n2-i,ilnrho)+f(:,:,n2+i,ilnrho)-2*lnrho0)
        enddo
      case default
        print*,"bc_ss_temp_z: invalid argument"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_temp_z
!***********************************************************************
    subroutine bc_ss_stemp_x(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: i
!
      if(ldebug) print*,'bc_ss_stemp_x: cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (cs2bot<=0.) print*, &
                        'bc_ss_stemp_x: cannot have cs2bot<=0'
        do i=1,nghost
          f(l1-i,:,:,iss) = f(l1+i,:,:,iss) &
               + gamma1/gamma*(f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho))
        enddo
!
!  top boundary
!
      case('top')
        if (cs2top<=0.) print*, &
                        'bc_ss_stemp_x: cannot have cs2top<=0'
        do i=1,nghost
          f(l2+i,:,:,iss) = f(l2-i,:,:,iss) &
               + gamma1/gamma*(f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho))
        enddo

      case default
        print*,"bc_ss_stemp_x: invalid argument"
        call stop_it("")
      endselect
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
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: i
!
      if(ldebug) print*,'bc_ss_stemp_y: cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (cs2bot<=0.) print*, &
                       'bc_ss_stemp_y: cannot have cs2bot<=0'
        do i=1,nghost
          f(:,m1-i,:,iss) = f(:,m1+i,:,iss) &
               + gamma1/gamma*(f(:,m1+i,:,ilnrho)-f(:,m1-i,:,ilnrho))
        enddo
!
!  top boundary
!
      case('top')
        if (cs2top<=0.) print*, &
                       'bc_ss_stemp_y: cannot have cs2top<=0'
        do i=1,nghost
          f(:,m2+i,:,iss) = f(:,m2-i,:,iss) &
               + gamma1/gamma*(f(:,m2-i,:,ilnrho)-f(:,m2+i,:,ilnrho))
        enddo

      case default
        print*,"bc_ss_stemp_y: invalid argument"
        call stop_it("")
      endselect
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
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: i
!
      if(ldebug) print*,'bc_ss_stemp_z: cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
          if (cs2bot<=0.) print*, &
                                  'bc_ss_stemp_z: cannot have cs2bot<=0'
          do i=1,nghost
             f(:,:,n1-i,iss) = f(:,:,n1+i,iss) &
                  + gamma1/gamma*(f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho))
          enddo
!
!  top boundary
!
      case('top')
        if (cs2top<=0.) print*, &
                 'bc_ss_stemp_z: cannot have cs2top<=0'
         do i=1,nghost
           f(:,:,n2+i,iss) = f(:,:,n2-i,iss) &
                + gamma1/gamma*(f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho))
         enddo
      case default
        print*,"bc_ss_stemp_z: invalid argument"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_stemp_z
!***********************************************************************
    subroutine bc_ss_energy(f,topbot)
!
!  boundary condition for entropy
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the entropy module
!  26-aug-2003/tony: distributed across ionization modules
!
      use Mpicomm, only: stop_it
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my) :: cs2_2d
      integer :: i
!
!  The 'ce' boundary condition for entropy makes the energy constant at
!  the boundaries.
!  This assumes that the density is already set (ie density must register
!  first!)
!
    select case(topbot)
!
! Bottom boundary
!
    case('bot')
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      cs2_2d=cs20*exp(gamma1*f(:,:,n1,ilnrho)+gamma*f(:,:,n1,iss))
      do i=1,nghost
         f(:,:,n1-i,iss)=1./gamma*(-gamma1*f(:,:,n1-i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo

!
! Top boundary
!
    case('top')
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      cs2_2d=cs20*exp(gamma1*f(:,:,n2,ilnrho)+gamma*f(:,:,n2,iss))
      do i=1,nghost
         f(:,:,n2+i,iss)=1./gamma*(-gamma1*f(:,:,n2+i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo
    case default
      print*,"bc_ss_energy: invalid argument"
      call stop_it("")
    endselect

    end subroutine bc_ss_energy
!***********************************************************************
endmodule ionization
