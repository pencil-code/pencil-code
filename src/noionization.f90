! $Id: noionization.f90,v 1.54 2003-08-09 19:47:27 mee Exp $

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

  interface ioncalc_ss                  ! Overload the 'ioncalc_ss' function
    module procedure ioncalc_ss_penc
    module procedure ioncalc_ss_point
  end interface

  interface ionget
    module procedure ionget_pencil
    module procedure ionget_point
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
  real :: yH0=0.5,yHacc=1e-5,xHe=0.1

  ! input parameters
  integer :: dummy_ni 
  namelist /ionization_init_pars/ yH0,yHacc,xHe

  ! run parameters
  namelist /ionization_run_pars/ yH0,yHacc,xHe

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
           "$Id: noionization.f90,v 1.54 2003-08-09 19:47:27 mee Exp $")
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
    subroutine ioncalc_ss_point(lnrho,TT,ss)
!
      real,intent(in) :: lnrho,TT
      real, intent(out) :: ss
!
      ss=(log(TT)-gamma1*(lnrho-lnrho0)-alog(cs20/gamma1))/gamma
!
    end subroutine ioncalc_ss_point
!***********************************************************************
    subroutine ioncalc_ss_penc(lnrho,TT,ss)
!
      real, dimension(nx), intent(in) :: lnrho,TT
      real, dimension(nx), intent(out) :: ss
!
      ss=(log(TT)-gamma1*(lnrho-lnrho0)-alog(cs20/gamma1))/gamma
!
    end subroutine ioncalc_ss_penc
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
    endsubroutine thermodynamics_point
!***********************************************************************
endmodule ionization
