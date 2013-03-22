! $Id: centrifugal_balance.f90 19193 2012-06-30 12:55:46Z wdobler $
!
!  This module provide a way for users to specify custom initial
!  conditions.
!
!  The module provides a set of standard hooks into the Pencil Code
!  and currently allows the following customizations:
!
!   Description                               | Relevant function call
!  ------------------------------------------------------------------------
!   Initial condition registration            | register_special
!     (pre parameter read)                    |
!   Initial condition initialization          | initialize_special
!     (post parameter read)                   |
!                                             |
!   Initial condition for momentum            | initial_condition_uu
!   Initial condition for density             | initial_condition_lnrho
!   Initial condition for entropy             | initial_condition_ss
!   Initial condition for magnetic potential  | initial_condition_aa
!
!   And a similar subroutine for each module with an "init_XXX" call.
!   The subroutines are organized IN THE SAME ORDER THAT THEY ARE CALLED.
!   First uu, then lnrho, then ss, then aa, and so on.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
! HOW TO USE THIS FILE
! --------------------
!
! Change the line above to
!    linitial_condition = .true.
! to enable use of custom initial conditions.
!
! The rest of this file may be used as a template for your own initial
! conditions. Simply fill out the prototypes for the features you want
! to use.
!
! Save the file with a meaningful name, e.g. mhs_equilibrium.f90, and
! place it in the $PENCIL_HOME/src/initial_condition directory. This
! path has been created to allow users to optionally check their
! contributions in to the Pencil Code SVN repository. This may be
! useful if you are working on/using an initial condition with
! somebody else or may require some assistance from one from the main
! Pencil Code team. HOWEVER, less general initial conditions should
! not go here (see below).
!
! You can also place initial condition files directly in the run
! directory. Simply create the folder 'initial_condition' at the same
! level as the *.in files and place an initial condition file there.
! With pc_setupsrc this file is linked automatically into the local
! src directory. This is the preferred method for initial conditions
! that are not very general.
!
! To use your additional initial condition code, edit the
! Makefile.local in the src directory under the run directory in which
! you wish to use your initial condition. Add a line that says e.g.
!
!    INITIAL_CONDITION =   initial_condition/mhs_equilibrium
!
! Here mhs_equilibrium is replaced by the filename of your new file,
! not including the .f90 extension.
!
! This module is based on Tony's special module.
!
module InitialCondition
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use EquationOfState
!
  implicit none
!
  include '../initial_condition.h'
!
!!  integer :: dummy
!
!  real    :: Hd=0.5, eps_dtog=1.0, vdampl_dust=1.0
!  real    :: OOg=1.0, input_fac=1.0, tau=0.5
!  logical :: ldragforce_gas=.true.
  real    :: Hd=1., eps_dtog=1., vdampl_dust=1.
  real    :: OOg=1., input_fac=1., tau=1.
  logical :: ldragforce_gas=.true.
  real    :: cb20=0.0
!
  namelist /initial_condition_pars/ Hd,OOg,eps_dtog,vdampl_dust,tau, &
       ldragforce_gas,input_fac, cb20
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  07-may-09/wlad: coded
!
      if (lroot) call svn_id( &
         "$Id: centrifugal_balance.f90 19193 2012-06-30 12:55:46Z wdobler $")
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
!  This subroutine is a general routine that takes
!  the gravity acceleration and adds the centrifugal force
!  that numerically balances it.
!
!  Pressure corrections to ensure centrifugal equilibrium are
!  added in the respective modules
!
!  24-feb-05/wlad: coded
!  04-jul-07/wlad: generalized for any shear
!  08-sep-07/wlad: moved here from initcond
!
      real, dimension (mx,my,mz,mfarray) :: f
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  07-may-09/wlad: coded
!
      use EquationOfState, only: cs20,rho0
      use Sub, only: hypergeometric2F1
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: lntmp,tmp1,tmp2
!
      real :: a1,b1,c1,d1,f1,chi
      real :: expp,expm,fac
      integer :: i,j
!
      a1=cb20 ; b1=cs20*eps_dtog ; c1=Hd**2
      d1=eps_dtog*vdampl_dust/tau    ; f1=OOg**2
!
      if (a1/=0) then
        chi=f1*c1/a1 ! (Hd/Hg)**2
      else
        chi=0.
      endif

!
      do m=m1,m2
        do n=n1,n2
!
          if (a1/=0) then 
            expp=exp( z(n)**2/(2*c1))
            expm=exp(-z(n)**2/(2*c1))
!
!  Reduces to isothermal when b=0 (no photoelectric effect)
!
            lntmp=z(n)**2/(2*c1) - (1+chi)*log(a1*expp+b1)
            tmp1=rho0*exp(lntmp)
!
!  Effect of drag backreaction with eps=eps(z)
!
            if (b1/=0 .and. ldragforce_gas) then 
              fac=d1*c1/(b1*(1-chi)) * (b1/a1*expm+1)**(-chi)/(a1/b1*expp+1)
              do i=l1,l2
                j=i-l1+1
                tmp2(j)=input_fac*fac*&
                     hypergeometric2F1(-chi,1-chi,2-chi,-b1/a1*expm,tol=1e-2)
              enddo
            else
              tmp2=0.
            endif
!
          else
            !only photoelectric
            tmp1=1.
          endif
!
          print*,z(n),c1,chi,a1,expp,b1

            f(l1:l2,m,n,ilnrho) = log(tmp1+tmp2)
!
        enddo
      enddo
!
      if (lroot) then
        print*,"min(lnrho)=",minval(f(l1:l2,m1:m2,n1:n2,ilnrho))
        print*,"max(lnrho)=",maxval(f(l1:l2,m1:m2,n1:n2,ilnrho))
      endif
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_nd(f)
!
      real, dimension (mx,my,mz,mfarray) :: f      
      real :: rho00, rhod00
!
      rho00=1.
      rhod00 = eps_dtog*rho00 !eps_dtog*Hrho/Hnd*rho00
      do n=n1,n2
         f(:,:,n,ind) = rhod00*exp(-z(n)**2/(2*Hd**2))
      enddo
!
      endsubroutine initial_condition_nd
!***********************************************************************
    subroutine initial_condition_uud(f)
!
!  Initial condition for dust velocity in the fluid approximation.
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      do n=1,mz
         f(:,:,n,iudz(1)) = -vdampl_dust*z(n)
      enddo
!
    endsubroutine initial_condition_uud
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
