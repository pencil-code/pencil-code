! $Id: eos_temperature_ionization.f90,v 1.8 2006-04-05 13:57:39 theine Exp $

!  Dummy routine for ideal gas

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ss,ee,pp,lnTT,cs2,glnTT,TT,TT1
! PENCILS PROVIDED yH,del2lnTT,cv,cv1,cp,cp1,mu1
! PENCILS PROVIDED glnTT,dppdlnTT,dppdlnrho,rho1gpp
!
!***************************************************************

module EquationOfState

  use Cparam
  use Cdata
  use Messages

  implicit none

  include 'eos.h'

  interface eoscalc ! Overload subroutine `eoscalc' function
    module procedure eoscalc_pencil   ! explicit f implicit m,n
    module procedure eoscalc_point    ! explicit lnrho, ss
    module procedure eoscalc_farray   ! explicit lnrho, ss
  end interface

  interface pressure_gradient ! Overload subroutine `pressure_gradient'
    module procedure pressure_gradient_farray  ! explicit f implicit m,n
    module procedure pressure_gradient_point   ! explicit lnrho, ss
  end interface

  ! integers specifying which independent variables to use in eoscalc
  integer, parameter :: ilnrho_ss=1,ilnrho_ee=2,ilnrho_pp=3,ilnrho_lnTT=4

  !  secondary parameters calculated in initialize
  real :: mu1_0,Rgas
  real :: TT_ion,lnTT_ion,TT_ion_,lnTT_ion_
  real :: ss_ion,kappa0
  real :: rho_H,rho_e,rho_e_,rho_He
  real :: lnrho_H,lnrho_e,lnrho_e_,lnrho_He

  real :: xHe=0.1,yH_const=0.0
  logical :: lconst_yH=.false.

  ! input parameters
  namelist /eos_init_pars/ xHe,lconst_yH,yH_const

  ! run parameters
  namelist /eos_run_pars/ xHe,lconst_yH,yH_const

  real :: cs0=impossible, rho0=impossible, cp=impossible
  real :: cs20=impossible, lnrho0=impossible
  logical :: lcalc_cp=.false. 
!  real :: gamma=impossible, gamma1=impossible, gamma11=impossible
  real :: gamma=5./3., gamma1=impossible, gamma11=impossible
  real :: cs2bot=impossible, cs2top=impossible 
  real :: cs2cool=impossible
  real :: mpoly=impossible, mpoly0=impossible
  real :: mpoly1=impossible, mpoly2=impossible
  integer :: isothtop=0
  real, dimension (3) :: beta_glnrho_global=impossible
  real, dimension (3) :: beta_glnrho_scaled=impossible

  contains

!***********************************************************************
    subroutine register_eos()
!
!  14-jun-03/axel: adapted from register_eos
!
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call fatal_error('register_eos','module registration called twice')
      first = .false.
!
      leos=.true.
      leos_temperature_ionization=.true.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           '$Id: eos_temperature_ionization.f90,v 1.8 2006-04-05 13:57:39 theine Exp $')
!
    endsubroutine register_eos
!***********************************************************************
    subroutine initialize_eos()
!
!  called by run.f90 after reading parameters, but before the time loop
!
      if (lroot) print*,'initialize_eos: ENTER'
!
!  Useful constants for ionization
!
      mu1_0 = 1/(1 + 4*xHe)
      Rgas = k_B/m_H  !(Tobi, this doesn't agree with the literature...)
      !Rgas = k_B/m_u
      TT_ion = chiH/k_B
      lnTT_ion = log(TT_ion)
      TT_ion_ = chiH_/k_B
      lnTT_ion_ = log(chiH_/k_B)
      rho_H = (1/mu1_0)*m_H*((m_H/hbar)*(chiH/hbar)/(2*pi))**(1.5)
      lnrho_H = log(rho_H)
      rho_e = (1/mu1_0)*m_H*((m_e/hbar)*(chiH/hbar)/(2*pi))**(1.5)
      lnrho_e = log(rho_e)
      rho_He = (1/mu1_0)*m_H*((4*m_H/hbar)*(chiH/hbar)/(2*pi))**(1.5)
      lnrho_He = log(rho_He)
      rho_e_ = (1/mu1_0)*m_H*((m_e/hbar)*(chiH_/hbar)/(2*pi))**(1.5)
      lnrho_e_ = log(rho_e_)
      kappa0 = sigmaH_*mu1_0/m_H/4.0
!
!  write scale non-free constants to file; to be read by idl
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro')
        write (1,*) 'mu1_0=',mu1_0
        write (1,*) 'Rgas=',Rgas
        write (1,*) 'TT_ion=',TT_ion
        write (1,*) 'lnTT_ion=',lnTT_ion
        write (1,*) 'TT_ion_=',TT_ion_
        write (1,*) 'lnTT_ion_=',lnTT_ion_
        write (1,*) 'lnrho_H=',lnrho_H
        write (1,*) 'lnrho_e=',lnrho_e
        write (1,*) 'lnrho_He=',lnrho_He
        write (1,*) 'lnrho_e_=',lnrho_e_
        write (1,*) 'kappa0=',kappa0
        close (1)
      endif

    endsubroutine initialize_eos
!*******************************************************************
    subroutine select_eos_variable(variable,findex)
!
!  dummy (but to be changed)
!
      character (len=*), intent(in) :: variable
      integer, intent(in) :: findex

      if (NO_WARN) print *,variable,findex

    endsubroutine select_eos_variable
!***********************************************************************
    subroutine pencil_criteria_eos()
!
!  dummy (but to be changed)
!
    endsubroutine pencil_criteria_eos
!***********************************************************************
    subroutine pencil_interdep_eos(lpencil_in)
!
!  dummy (but to be changed)
!
      logical, dimension(npencils) :: lpencil_in

      if (lpencil_in(i_cs2)) then
        lpencil_in(i_rho1)=.true.
        lpencil_in(i_cv1)=.true.
        lpencil_in(i_TT1)=.true.
        lpencil_in(i_dppdlnrho)=.true.
        lpencil_in(i_dppdlnTT)=.true.
      endif

      if (lpencil_in(i_rho1gpp)) then
        lpencil_in(i_rho1)=.true.
        lpencil_in(i_dppdlnrho)=.true.
        lpencil_in(i_glnrho)=.true.
        lpencil_in(i_dppdlnTT)=.true.
        lpencil_in(i_glnTT)=.true.
      endif

      if (lpencil_in(i_cv1)) lpencil_in(i_cv)=.true.

      if (lpencil_in(i_cv)) then
        lpencil_in(i_yH)=.true.
        lpencil_in(i_TT1)=.true.
        lpencil_in(i_mu1)=.true.
      endif

      if (lpencil_in(i_cp1)) lpencil_in(i_cp)=.true.

      if (lpencil_in(i_cp)) then
        lpencil_in(i_yH)=.true.
        lpencil_in(i_TT1)=.true.
        lpencil_in(i_mu1)=.true.
      endif

      if (lpencil_in(i_pp)) then
        lpencil_in(i_mu1)=.true.
        lpencil_in(i_rho)=.true.
        lpencil_in(i_TT)=.true.
      endif

      if (lpencil_in(i_yH).and.not(lconst_yH)) then
        lpencil_in(i_rho1)=.true.
        lpencil_in(i_TT1)=.true.
      endif

      if (lpencil_in(i_dppdlnTT)) then
        lpencil_in(i_yH)=.true.
        lpencil_in(i_TT1)=.true.
        lpencil_in(i_pp)=.true.
      endif

      if (lpencil_in(i_dppdlnrho)) then
        lpencil_in(i_yH)=.true.
        lpencil_in(i_pp)=.true.
      endif

      if (lpencil_in(i_mu1)) lpencil_in(i_yH)=.true.

      if (lpencil_in(i_ee)) then
        lpencil_in(i_mu1)=.true.
        lpencil_in(i_TT)=.true.
        lpencil_in(i_yH)=.true.
      endif

      if (lpencil_in(i_ss)) then
        lpencil_in(i_yH)=.true.
        lpencil_in(i_lnrho)=.true.
        lpencil_in(i_lnTT)=.true.
      endif

      if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
      if (lpencil_in(i_TT1)) lpencil_in(i_TT)=.true.

      if (NO_WARN) print *,lpencil_in

    endsubroutine pencil_interdep_eos
!*******************************************************************
    subroutine calc_pencils_eos(f,p)
!
!  dummy (but to be changed)
!
      use Sub, only: grad,del2

      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p

      real, dimension (nx) :: rhs,yH_term,tmp,sqrtrhs
      integer :: i

      if (NO_WARN) print *,f,p

!
!  Temperature
!
      if (lpencil(i_lnTT)) p%lnTT=f(l1:l2,m,n,ilnTT)
      if (lpencil(i_TT)) p%TT=exp(p%lnTT)
      if (lpencil(i_TT1)) p%TT1=1/p%TT

!
!  Temperature laplacian and gradient
!
      if (lpencil(i_glnTT)) call grad(f,ilnTT,p%glnTT)
      if (lpencil(i_del2lnTT)) call del2(f,ilnTT,p%del2lnTT)

!
!  Ionization fraction
!
      if (lpencil(i_yH)) then
        if (lconst_yH) then
          p%yH = yH_const
        else
          where (TT_ion*p%TT1 < -log(tiny(TT_ion)))
            rhs = rho_e*p%rho1*(p%TT1*TT_ion)**(-1.5)*exp(-TT_ion*p%TT1)
            sqrtrhs = sqrt(rhs)
            p%yH = 2*sqrtrhs/(sqrtrhs+sqrt(4+rhs))
          elsewhere
            p%yH = 0
          endwhere
        endif
      endif

!
!  Mean molecular weight
!
      if (lpencil(i_mu1)) p%mu1 = mu1_0*(1 + p%yH + xHe)

!
!  Pressure
!
      if (lpencil(i_pp)) p%pp = Rgas*p%mu1*p%rho*p%TT

!
!  Common term involving the ionization fraction
!
      if (lpencil(i_cv).or.lpencil(i_dppdlnTT).or.lpencil(i_dppdlnrho)) then
        yH_term = p%yH*(1-p%yH)/((2-p%yH)*(1+p%yH+xHe))
      endif

!
!  Specific heat at constant volume (i.e. density)
!
      if (lpencil(i_cv)) then
        tmp = 1.5 + TT_ion*p%TT1
        p%cv = Rgas*p%mu1*(1.5 + yH_term*tmp**2)
      endif

      if (lpencil(i_cv1)) p%cv1=1/p%cv

!
!  Specific heat at constant pressure
!
      if (lpencil(i_cp)) then
        tmp = 2.5+TT_ion*p%TT1
        p%cp = Rgas*p%mu1*(2.5 + yH_term*tmp**2)
      endif

      if (lpencil(i_cp1)) p%cp1=1/p%cp

!
!  Coefficients for the pressure gradient
!
      if (lpencil(i_dppdlnTT)) then
        tmp = 1.5+TT_ion*p%TT1
        p%dppdlnTT = p%pp*(1 + yH_term*tmp)
      endif

      if (lpencil(i_dppdlnrho)) p%dppdlnrho = p%pp*(1 - yH_term)

!
!  Logarithmic pressure gradient
!
      if (lpencil(i_rho1gpp)) then
        do i=1,3
          p%rho1gpp(:,i)=p%rho1*(p%dppdlnrho*p%glnrho(:,i) + &
                                 p%dppdlnTT*p%glnTT(:,i))
        enddo
      endif

!
!  Sound speed
!
      if (lpencil(i_cs2)) then
        p%cs2 = p%rho1*(p%rho1*p%cv1*p%TT1*(p%dppdlnTT)**2 + p%dppdlnrho)
      endif

!
!  Energy per unit mass
!AB: Tobi, is this correct?
!
      if (lpencil(i_ee)) p%ee = 1.5*Rgas*p%mu1*p%TT + p%yH*Rgas*mu1_0*TT_ion

!
!  Entropy per unit mass
!  The contributions from each particle species contain the mixing entropy
!
      if (lpencil(i_ss)) then
        tmp = 2.5 - 1.5*(lnTT_ion-p%lnTT)
        p%ss = tmp + lnrho_H - p%lnrho
        where (p%yH > 0)
          p%ss = p%ss + p%yH*(tmp + lnrho_e - p%lnrho - log(p%yH))
        endwhere
        if (xHe > 0) then
          p%ss = p%ss + xHe*(tmp + lnrho_He - p%lnrho - log(xHe))
        endif
        p%ss = Rgas*mu1_0*p%ss
      endif

    endsubroutine calc_pencils_eos
!*******************************************************************
    subroutine getmu(mu)
!
!  Calculate average particle mass in the gas relative to
!
!   12-aug-03/tony: implemented dummy
!
      real, intent(out) :: mu

      mu=0.
    endsubroutine getmu
!*******************************************************************
    subroutine rprint_eos(lreset,lwrite)
!
!  Writes iyH and ilnTT to index.pro file
!
!  02-apr-03/tony: implemented dummy
!
! 
      logical :: lreset
      logical, optional :: lwrite
!   
      if (NO_WARN) print*,lreset,present(lwrite)  !(keep compiler quiet)
!   
    endsubroutine rprint_eos
!***********************************************************************
    subroutine ioninit(f)
!   
      real, dimension (mx,my,mz,mvar+maux), intent(inout) :: f
!   
      if(NO_WARN) print*,f  !(keep compiler quiet)
!
    endsubroutine ioninit
!***********************************************************************
    subroutine ioncalc(f)
!
    real, dimension (mx,my,mz,mvar+maux) :: f
!
    if(NO_WARN) print*,f  !(keep compiler quiet)
!
    endsubroutine ioncalc
!***********************************************************************
    subroutine eosperturb(f,psize,ee,pp)

      real, dimension(mx,my,mz,mvar+maux), intent(inout) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(in), optional :: ee,pp

      call not_implemented("eosperturb")

    end subroutine eosperturb
!***********************************************************************
    subroutine getdensity(EE,TT,yH,rho)
      
      real, intent(in) :: EE,TT,yH
      real, intent(inout) :: rho

      call fatal_error('getdensity','SHOULD NOT BE CALLED WITH NOEOS')
      rho = 0. 
      if (NO_WARN) print*,yH,EE,TT  !(keep compiler quiet)

    end subroutine getdensity
!***********************************************************************
    subroutine pressure_gradient_farray(f,cs2,cp1tilde)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   02-apr-04/tony: implemented dummy
!
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(nx), intent(out) :: cs2,cp1tilde
!
!
      cs2=impossible
      cp1tilde=impossible
      call fatal_error('pressure_gradient_farray','SHOULD NOT BE CALLED WITH NOEOS')
      if (NO_WARN) print*,f  !(keep compiler quiet)
!
    endsubroutine pressure_gradient_farray
!***********************************************************************
    subroutine pressure_gradient_point(lnrho,ss,cs2,cp1tilde)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   02-apr-04/tony: implemented dummy
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: cs2,cp1tilde
!
      cs2=impossible
      cp1tilde=impossible
      call fatal_error('pressure_gradient_point','SHOULD NOT BE CALLED WITH NOEOS')
      if (NO_WARN) print*,lnrho,ss  !(keep compiler quiet)
!
    endsubroutine pressure_gradient_point
!***********************************************************************
    subroutine temperature_gradient(f,glnrho,gss,glnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   02-apr-04/tony: implemented dummy
!
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(nx,3), intent(in) :: glnrho,gss
      real, dimension(nx,3), intent(out) :: glnTT
!
!
      glnTT=0.
      call fatal_error('temperature_gradient','SHOULD NOT BE CALLED WITH NOEOS')
      if (NO_WARN) print*,f,glnrho,gss  !(keep compiler quiet)
    endsubroutine temperature_gradient
!***********************************************************************
    subroutine temperature_laplacian(f,del2lnrho,del2ss,del2lnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   12-dec-05/tony: adapted from subroutine temperature_gradient
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(nx), intent(in) :: del2lnrho,del2ss
      real, dimension(nx), intent(out) :: del2lnTT
!
      call fatal_error('temperature_laplacian','SHOULD NOT BE CALLED WITH NOEOS')
!
      del2lnTT=0.
      if (NO_WARN) print*,f,del2lnrho,del2ss !(keep compiler quiet)
    endsubroutine temperature_laplacian
!***********************************************************************
    subroutine temperature_hessian(f,hlnrho,hss,hlnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally hlnPP and hlnTT
!   hP/rho=cs2*(hlnrho+cp1tilde*hss)
!
!   13-may-04/tony: adapted from idealgas dummy
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(nx,3), intent(in) :: hlnrho,hss
      real, dimension(nx,3) :: hlnTT
!
      call fatal_error('temperature_hessian','now I do not believe you'// &
                           ' intended to call this!')
!
      if (NO_WARN) print*,f,hlnrho,hss !(keep compiler quiet)
!
    endsubroutine temperature_hessian
!***********************************************************************
    subroutine eoscalc_farray(f,psize,lnrho,ss,yH,mu1,lnTT,ee,pp,kapparho)
!
!   Calculate thermodynamical quantities
!
!   04-apr-06/tobi: Adapted for this EOS module
!
      use Mpicomm, only: stop_it

      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(out), optional :: lnrho,ss
      real, dimension(psize), intent(out), optional :: yH,lnTT,mu1
      real, dimension(psize), intent(out), optional :: ee,pp,kapparho

      real, dimension(psize) :: lnrho_,lnTT_
      real, dimension(psize) :: rho1,TT1,rhs,tmp

      select case (psize)

      case (nx)
        lnrho_=f(l1:l2,m,n,ilnrho)
        lnTT_=f(l1:l2,m,n,ilnTT)

      case (mx)
        lnrho_=f(:,m,n,ilnrho)
        lnTT_=f(:,m,n,ilnTT)

      case default
        call stop_it("eoscalc: no such pencil size")

      end select


      if (present(lnrho)) lnrho=lnrho_

      if (present(lnTT)) lnTT=lnTT_

      if (present(yH).or.present(ss).or.present(ee).or.present(pp).or. &
          present(mu1).or.present(kapparho)) then
        rho1 = exp(-lnrho_)
        TT1 = exp(-lnTT_)
        if (lconst_yH) then
          yH = yH_const
        else
          where (TT_ion*TT1 < -log(tiny(TT_ion)))
            rhs = rho_e*rho1*(TT1*TT_ion)**(-1.5)*exp(-TT_ion*TT1)
            yH = 2*sqrt(rhs)/(sqrt(rhs)+sqrt(4+rhs))
          elsewhere
            yH = 0
          endwhere
        endif
      endif

      if (present(mu1).or.present(ss).or.present(ee).or.present(pp)) then
        mu1 = mu1_0*(1 + yH + xHe)
      endif

      if (present(ee)) ee = 1.5*Rgas*mu1/TT1 + yH*Rgas*mu1_0*TT_ion

      if (present(pp)) pp = Rgas*mu1/(rho1*TT1)

      if (present(ss)) then
        tmp = 2.5 - 1.5*(lnTT_ion-lnTT_)
        ss = tmp + lnrho_H - lnrho_
        where (yH > 0)
          ss = ss + yH*(tmp + lnrho_e - lnrho - log(yH))
        endwhere
        if (xHe > 0) then
          ss = ss + xHe*(tmp + lnrho_He - lnrho - log(xHe))
        endif
        ss = Rgas*mu1_0*ss
      endif

      if (present(kapparho)) kapparho=0

    endsubroutine eoscalc_farray
!***********************************************************************
    subroutine eoscalc_point(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp)
!
!   Calculate thermodynamical quantities
!
!   02-apr-04/tony: implemented dummy
!
      integer, intent(in) :: ivars
      real, intent(in) :: var1,var2
      real, intent(out), optional :: lnrho,ss
      real, intent(out), optional :: yH,lnTT
      real, intent(out), optional :: ee,pp
!
!
      call fatal_error('eoscalc_point','SHOULD NOT BE CALLED WITH NOEOS')

      if (present(lnrho)) lnrho=0.
      if (present(ss)) ss=0.
      if (present(yH)) yH=impossible
      if (present(lnTT)) lnTT=0.
      if (present(ee)) ee=0.
      if (present(pp)) pp=0.
      if (NO_WARN) print*,ivars,var1,var2  !(keep compiler quiet)
!
    endsubroutine eoscalc_point
!***********************************************************************
    subroutine eoscalc_pencil(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp)
!
!   Calculate thermodynamical quantities
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1tilde to
!                   subroutine pressure_gradient
!
      integer, intent(in) :: ivars
      real, dimension(nx), intent(in) :: var1,var2
      real, dimension(nx), intent(out), optional :: lnrho,ss
      real, dimension(nx), intent(out), optional :: yH,lnTT
      real, dimension(nx), intent(out), optional :: ee,pp
!
!
      call fatal_error('eoscalc_pencil','SHOULD NOT BE CALLED WITH NOEOS')

      if (present(lnrho)) lnrho=0.
      if (present(ss)) ss=0.
      if (present(yH)) yH=impossible
      if (present(lnTT)) lnTT=0.
      if (present(ee)) ee=0.
      if (present(pp)) pp=0.
      if (NO_WARN) print*,ivars,var1,var2  !(keep compiler quiet)
!
    endsubroutine eoscalc_pencil
!***********************************************************************
    subroutine get_soundspeed(lnTT,cs2)
!
!  Calculate sound speed for given temperature
!
!  02-apr-04/tony: dummy coded
!
      real, intent(in)  :: lnTT
      real, intent(out) :: cs2
!
      call not_implemented('get_soundspeed')
      cs2=0.
      if (NO_WARN) print*,lnTT
!
    end subroutine get_soundspeed
!***********************************************************************
    subroutine read_eos_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit !(keep compiler quiet)
!
    endsubroutine read_eos_init_pars
!***********************************************************************
    subroutine write_eos_init_pars(unit)
      integer, intent(in) :: unit
!
      if (NO_WARN) print*,unit !(keep compiler quiet)
!
    endsubroutine write_eos_init_pars
!***********************************************************************
    subroutine read_eos_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit !(keep compiler quiet)
!
    endsubroutine read_eos_run_pars
!***********************************************************************
    subroutine write_eos_run_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit !(keep compiler quiet)

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
      real, dimension(mx,my,mz,mvar+maux), intent(inout) :: f
      real, intent(in) :: T0
      real, dimension(nx) :: lnrho,ss
      real :: ss_offset=0.
!
!  if T0 is different from unity, we interpret
!  ss_offset = ln(T0)/gamma as an additive offset of ss
!
      if (T0/=1.) ss_offset=alog(T0)/gamma
!
      do n=n1,n2
      do m=m1,m2
        lnrho=f(l1:l2,m,n,ilnrho)
        ss=-gamma1*(lnrho-lnrho0)/gamma
          !+ other terms for sound speed not equal to cs_0
        f(l1:l2,m,n,iss)=ss+ss_offset
      enddo
      enddo
!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  The values calculated here may be revised in the entropy module.
!
      cs2bot=cs20
      cs2top=cs20
!
    endsubroutine isothermal_entropy
!***********************************************************************
    subroutine isothermal_lnrho_ss(f,T0,rho0)
!
!  Isothermal stratification for lnrho and ss (for yH=0!)
!
!  Currently only implemented for ionization_fixed.
!
      real, dimension(mx,my,mz,mvar+maux), intent(inout) :: f
      real, intent(in) :: T0,rho0
!
      if (NO_WARN) print*,f,T0,rho0
!
    endsubroutine isothermal_lnrho_ss
!***********************************************************************
    subroutine bc_ss_flux(f,topbot,hcond0,hcond1,Fheat,FheatK,chi, &
                lmultilayer,lcalc_heatcond_constchi)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  26-aug-2003/tony: distributed across ionization modules
!
      use Gravity
!
      real, intent(in) :: Fheat, FheatK, hcond0, hcond1, chi
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
!  ===============
!
      case('bot')
        if (lmultilayer) then
          if(headtt) print*,'bc_ss_flux: Fbot,hcond=',Fheat,hcond0*hcond1
        else
          if(headtt) print*,'bc_ss_flux: Fbot,hcond=',Fheat,hcond0
        endif
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
          tmp_xy=Fheat/(rho_xy*chi*cs2_xy)
        else
          tmp_xy=FheatK/cs2_xy
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
        if (lmultilayer) then
          if(headtt) print*,'bc_ss_flux: Ftop,hcond=',Fheat,hcond0*hcond1
        else
          if(headtt) print*,'bc_ss_flux: Ftop,hcond=',Fheat,hcond0
        endif
!
!  calculate Ftop/(K*cs2)
!
        rho_xy=exp(f(:,:,n2,ilnrho))
        cs2_xy=cs20*exp(gamma1*(f(:,:,n2,ilnrho)-lnrho0)+gamma*f(:,:,n2,iss))
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy. 
!
        if(lcalc_heatcond_constchi) then
          tmp_xy=Fheat/(rho_xy*chi*cs2_xy)
        else
          tmp_xy=FheatK/cs2_xy
        endif
!
!  enforce ds/dz + gamma1/gamma*dlnrho/dz = - gamma1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n2+i,iss)=f(:,:,n2-i,iss)+gamma1/gamma* &
              (f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho)-2*i*dz*tmp_xy)
        enddo
      case default
        call fatal_error('bc_ss_flux','invalid argument')
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
!  23-jun-2003/tony: implemented for leos_fixed_ionization
!  26-aug-2003/tony: distributed across ionization modules
!
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
        if ((bcz1(ilnrho) /= 'a2') .and. (bcz1(ilnrho) /= 'a3')) &
          call fatal_error('bc_ss_temp_old','Inconsistent boundary conditions 3.')
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
        if ((bcz1(ilnrho) /= 'a2') .and. (bcz1(ilnrho) /= 'a3')) &
          call fatal_error('bc_ss_temp_old','Inconsistent boundary conditions 3.')
        if (ldebug) print*, &
                   'bc_ss_temp_old: set top temperature - cs2top=',cs2top
        if (cs2top<=0.) print*, &
                   'bc_ss_temp_old: cannot have cs2top<=0'
  !     if (bcz1(ilnrho) /= 'a2') &
  !          call fatal_error(bc_ss_temp_old','Inconsistent boundary conditions 4.')
        tmp_xy = (-gamma1*(f(:,:,n2,ilnrho)-lnrho0) &
                 + alog(cs2top/cs20)) / gamma
        f(:,:,n2,iss) = tmp_xy
        do i=1,nghost
          f(:,:,n2+i,iss) = 2*tmp_xy - f(:,:,n2-i,iss)
        enddo 
      case default
        call fatal_error('bc_ss_temp_old','invalid argument')
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
        call fatal_error('bc_ss_temp_x','invalid argument')
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
        call fatal_error('bc_ss_temp_y','invalid argument')
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
        call fatal_error('bc_ss_temp_z','invalid argument')
      endselect
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
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: tmp
      integer :: i
!
      if(ldebug) print*,'bc_lnrho_temp_z: cs20,cs0=',cs20,cs0
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
                 'bc_lnrho_temp_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0. .and. lroot) print*, &
                 'bc_lnrho_temp_z: cannot have cs2bot<=0'
        tmp = 2/gamma*log(cs2bot/cs20)
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        f(:,:,n1,iss) = 0.5*tmp - gamma1/gamma*(f(:,:,n1,ilnrho)-lnrho0)
        do i=1,nghost; f(:,:,n1-i,iss) = 2*f(:,:,n1,iss)-f(:,:,n1+i,iss); enddo
!
!  set density in the ghost zones so that dlnrho/dz + ds/dz = gz/cs2bot
!  for the time being, we don't worry about lnrho0 (assuming that it is 0)
!
        tmp=-gravz/cs2bot
        do i=1,nghost
          f(:,:,n1-i,ilnrho) = f(:,:,n1+i,ilnrho) +f(:,:,n1+i,iss) &
                                                  -f(:,:,n1-i,iss) +2*i*dz*tmp
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                    'bc_lnrho_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0. .and. lroot) print*, &
                    'bc_lnrho_temp_z: cannot have cs2top<=0'
        tmp = 2/gamma*log(cs2top/cs20)
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        f(:,:,n2,iss) = 0.5*tmp - gamma1/gamma*(f(:,:,n2,ilnrho)-lnrho0)
        do i=1,nghost; f(:,:,n2+i,iss) = 2*f(:,:,n2,iss)-f(:,:,n2-i,iss); enddo
!
!  set density in the ghost zones so that dlnrho/dz + ds/dz = gz/cs2top
!  for the time being, we don't worry about lnrho0 (assuming that it is 0)
!
        tmp=gravz/cs2top
        do i=1,nghost
          f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho) +f(:,:,n2-i,iss) &
                                                  -f(:,:,n2+i,iss) +2*i*dz*tmp
        enddo

      case default
        call fatal_error('bc_lnrho_temp_z','invalid argument')
      endselect
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
      use Cdata
      use Gravity, only: lnrho_bot,lnrho_top,ss_bot,ss_top
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: i
!
      if(ldebug) print*,'bc_lnrho_pressure_z: cs20,cs0=',cs20,cs0
!
!  Constant pressure, i.e. antisymmetric
!  This assumes that the entropy is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('top')
        if (ldebug) print*,'bc_lnrho_pressure_z: lnrho_top,ss_top=',lnrho_top,ss_top
!
!  fix entropy if inflow (uz>0); otherwise leave s unchanged
!  afterwards set s antisymmetrically about boundary value
!
        if(lentropy) then
!         do m=m1,m2
!         do l=l1,l2
!           if (f(l,m,n1,iuz)>=0) then
!             f(l,m,n1,iss)=ss_bot
!           else
!             f(l,m,n1,iss)=f(l,m,n1+1,iss)
!           endif
!         enddo
!         enddo
          f(:,:,n2,iss)=ss_top
          do i=1,nghost; f(:,:,n2+i,iss)=2*f(:,:,n2,iss)-f(:,:,n2-i,iss); enddo
!
!  set density value such that pressure is constant at the bottom
!
          f(:,:,n2,ilnrho)=lnrho_top+ss_top-f(:,:,n2,iss)
        else
          f(:,:,n2,ilnrho)=lnrho_top
        endif
!
!  make density antisymmetric about boundary
!  another possibility might be to enforce hydrostatics
!  ie to set dlnrho/dz=-g/cs^2, assuming zero entropy gradient
!
        do i=1,nghost
          f(:,:,n2+i,ilnrho)=2*f(:,:,n2,ilnrho)-f(:,:,n2-i,ilnrho)
        enddo
!
!  top boundary
!
      case('bot')
        if (ldebug) print*,'bc_lnrho_pressure_z: lnrho_bot,ss_bot=',lnrho_bot,ss_bot
!
!  fix entropy if inflow (uz>0); otherwise leave s unchanged
!  afterwards set s antisymmetrically about boundary value
!
        if(lentropy) then
!         do m=m1,m2
!         do l=l1,l2
!           if (f(l,m,n1,iuz)>=0) then
!             f(l,m,n1,iss)=ss_bot
!           else
!             f(l,m,n1,iss)=f(l,m,n1+1,iss)
!           endif
!         enddo
!         enddo
          f(:,:,n1,iss)=ss_bot
          do i=1,nghost; f(:,:,n1-i,iss)=2*f(:,:,n1,iss)-f(:,:,n1+i,iss); enddo
!
!  set density value such that pressure is constant at the bottom
!
          f(:,:,n1,ilnrho)=lnrho_bot+ss_bot-f(:,:,n1,iss)
        else
          f(:,:,n1,ilnrho)=lnrho_bot
        endif
!
!  make density antisymmetric about boundary
!  another possibility might be to enforce hydrostatics
!  ie to set dlnrho/dz=-g/cs^2, assuming zero entropy gradient
!
        do i=1,nghost
          f(:,:,n1-i,ilnrho)=2*f(:,:,n1,ilnrho)-f(:,:,n1+i,ilnrho)
        enddo
!
      case default
        call fatal_error('bc_lnrho_pressure_z','invalid argument')
      endselect
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: tmp
      integer :: i
!
      if(ldebug) print*,'bc_ss_temp2_z: cs20,cs0=',cs20,cs0
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
                   'bc_ss_temp2_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp2_z: cannot have cs2bot<=0'
        tmp = 1/gamma*alog(cs2bot/cs20)
        do i=0,nghost
          f(:,:,n1-i,iss) = tmp &
               - gamma1/gamma*(f(:,:,n1-i,ilnrho)-lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                     'bc_ss_temp2_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*,'bc_ss_temp2_z: cannot have cs2top<=0'
        tmp = 1/gamma*alog(cs2top/cs20)
        do i=0,nghost
          f(:,:,n2+i,iss) = tmp &
               - gamma1/gamma*(f(:,:,n2+i,ilnrho)-lnrho0)
        enddo
      case default
        call fatal_error('bc_ss_temp2_z','invalid argument')
      endselect
!
    endsubroutine bc_ss_temp2_z
!***********************************************************************
    subroutine bc_ss_stemp_x(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
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
        call fatal_error('bc_ss_stemp_x','invalid argument')
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
        call fatal_error('bc_ss_stemp_y','invalid argument')
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
        call fatal_error('bc_ss_stemp_z','invalid argument')
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
      call fatal_error('bc_ss_energy','invalid argument')
    endselect

    end subroutine bc_ss_energy
!***********************************************************************
endmodule EquationOfState
