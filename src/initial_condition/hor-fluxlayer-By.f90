! $Id: hor-fluxlayer-By.f90 23094 2015-02-10 17:54:43Z mreinhardt@nordita.org $
!
!  Horizontal flux layer in B_y, with tanh profiles switching on (off) at z0 (z1)
!  Density and entropy also set appropriately.
!  To reproduce initial condition in Catteneo & Hughes, JFM, v. 196, p. 323, 1988
!
!  Written as initial_condition/hor-fluxlayer-By.f90,  04-mar-15/ngrs.
!  (From code originally added to magnetic.f90, initcond.f90 in apr-14, 
!   but never checked in.)
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
module InitialCondition
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use EquationOfState
  use Sub, only: step, der_step
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: amplaa=1.74110, widthaa=0.0333333, z0aa=-1.3, z1aa=-1.0
  !logical :: lpress_equil=.false., lpress_equil_via_ss=.false.
  !logical :: lpress_equil_alt=.false.
!
  namelist /initial_condition_pars/ &
      amplaa,widthaa,z0aa,z1aa !, &
      !lpress_equil,lpress_equil_via_ss,lpress_equil_alt
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
           "$Id: hor-fluxlayer-By.f90 23094 2015-02-10 17:54:43Z mreinhardt@nordita.org $")
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
    subroutine initial_condition_all(f,profiles)
!
!  Initializes all the f arrays in one call. This subroutine is called last.
!
      real, dimension (mx,my,mz,mfarray), optional, intent(inout):: f
      real, dimension (:,:),              optional, intent(out)  :: profiles
!
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_all
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
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
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
!
      use EquationOfState, only: gamma,gamma_m1,gamma1,cs20,rho0,lnrho0
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the vector potential.
!
!  04-mar-15/ngrs: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
!     call hfluxlayer_By(amplaa,f,iaa,z0aa,z1aa,widthaa)  ! grs: add CH88 B layer -- now done here:
!
      if (amplaa==0) then
        f(:,:,:,iaa:iaa+2)=0
        if (lroot) print*,'hor-fluxlayer-By: set variable to zero; iaa=',iaa
      else
        if (lroot) print*,'hor-fluxlayer-By: horizontal flux layer; iaa=',iaa
        if ((ip<=16).and.lroot) print*,'hfluxlayer_By: ampl,widthaa,zflayer0,zflayer1=',amplaa,widthaa,z0aa,z1aa
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,iaa  )=amplaa*widthaa/2.*( log(cosh((z(n)-z0aa)/widthaa)) - log(cosh((z(n)-z1aa)/widthaa)) )
          f(l1:l2,m,n,iaa+1)=0.0
          f(l1:l2,m,n,iaa+2)=0.0
        enddo; enddo
      endif
!
    endsubroutine initial_condition_aa
!********************************************************************
    subroutine read_initial_condition_pars(unit,iostat)
!
!  07-may-09/wlad: coded
!
      include '../unit.h'
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
