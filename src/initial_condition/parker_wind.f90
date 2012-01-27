! $Id: parker_wind.f90,v 1.0 2012-01-24 09:08:03 Dhruba, joern Exp $
!
! Initial condition for the parker wind in spherical coordinates.  
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
!
!!!!!!!
!! NB: This is a very early stage of the initial condition. DO NOT USE IT!!!
!! dhruba+ joern: very much a work in progress. 
!!!!!!!

module InitialCondition
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet

!
  implicit none
!
  include '../initial_condition.h'
!
  real :: rcrit=1., mdot=0.
  real, dimension(mx) :: vv=0.
  real, dimension(mx) :: den=0.
  
  namelist /initial_condition_pars/ &
       rcrit, mdot
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  26-jan-12/joern: coded
!
      if (lroot) call svn_id( &
           "$Id: parker_wind.f90, v 1.0 2012-01-24 09:08:03 Dhruba, joern Exp$")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  24-jan-12/joern: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer :: iy, iz
!
      if (.not. lanelastic) then
        call parker_wind_iteration(vv,den)
        do iy=m1,m2;do iz=n1,n2
          f(:,iy,iy,ilnrho)=log(den)
          f(:,iy,iz,iuu)=vv
          f(:,iy,iz,iuu+1)=0.
          f(:,iy,iz,iuu+2)=0.
        enddo; enddo
      endif

!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho 
!  will take care of converting it to linear 
!  density if you use ldensity_nolog
!
!  24-jan-12/joern: coded
!
      use SharedVariables
      use EquationOfState
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f      
!      
      if (lanelastic) then
        call parker_wind_iteration(vv,den)
      do iy=m1,m2;do iz=n1,n2
          f(:,iy,iy,ilnrho)=log(den)
          f(:,iy,iz,iux)=vv
          f(:,iy,iz,iuy)=0.
          f(:,iy,iz,iuz)=0.
        enddo; enddo
      endif

   call keep_compiler_quiet(f)

    endsubroutine initial_condition_lnrho
!********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize entropy.
!
!  07-sep-11/simon: coded
!

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f      
!
      call keep_compiler_quiet(f)
!      
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-sep-11/simon: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f      
!      
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine initial_condition_lncc(f)
!
!  Initialize entropy.
!
!  07-sep-11/simon: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f      
      integer :: l, m, n
!      
    endsubroutine initial_condition_lncc
!***********************************************************************
    subroutine read_initial_condition_pars(unit,iostat)
!
!  07-may-09/wlad: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=initial_condition_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=initial_condition_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!     
      integer, intent(in) :: unit
!
      write(unit,NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
    subroutine parker_wind_iteration(vel,rho)
!
      use SharedVariables
      use EquationOfState
!     
      real, dimension (mx), intent(out) :: vel
      real, dimension (mx), intent(out) :: rho
      real, dimension (mx) :: cs20logx,GM_r
      real :: Ecrit
      integer :: j
!
      Ecrit=0.5*cs20-cs20*log(cs0)-2*cs20*log(rcrit)-2*cs20
      cs20logx=2*cs20*log(x)
      GM_r=2.*rcrit*cs20/x
!     
      vel=x/rcrit 
      do j=1,2000
        where (x>=rcrit)
          vel=sqrt(2.*(Ecrit+cs20logx &
          +GM_r+cs20*log(vel)))
        elsewhere
          vel=exp(0.5*vel**2-Ecrit &
          -cs20logx-GM_r)
        endwhere
      enddo
!
      rho=Mdot/(4*pi*x**2*vel)
print*,'Ecrit=',Ecrit
print*, 'GM=', 2.*rcrit*cs20
print*,'Rho0=',log(rho(l1))
!
    endsubroutine parker_wind_iteration 
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
