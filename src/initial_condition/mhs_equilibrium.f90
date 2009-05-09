!  $Id$
!
!  Initial condition (density, magnetic field, velocity) 
!  for magnetohydrostatical equilibrium in a global accretion
!  disk with an imposed (cylindrically symmetric) sound speed 
!  profile in spherical coordinates. 
!
!  07-may-09/axel: adapted from noinitial_condition.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
module InitialCondition
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'initial_condition.h'
!
  real :: density_power_law,temperature_power_law
!
  namelist /initial_condition_init_pars/ &
      density_power_law,temperature_power_law
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/wlad: coded
!
!  Identify CVS/SVN version information.
!
      if (lroot) call cvs_id( &
           "$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Initialize any module variables which are parameter dependent
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
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
!  Initialize logarithmic density. init_lnrho 
!  will take care of converting it to linear 
!  density if you use ldensity_nolog
!
!  07-may-09/wlad: coded
!
      use EquationOfState, only: cs20,rho0
      use FArrayManager, only: farray_use_global
      use Gravity, only: potential
      use Sub, only: get_radial_distance,power_law,grad
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension (mx) :: rr_sph,rr_cyl,rr
      real, dimension (mx) :: lnrhomid,strat
      real, dimension (mx) :: cs2,tmp1,tmp2
!
      real, dimension (nx,3) :: grho_mn,glnrho_mn,glnTT_mn
      real, dimension (nx)   :: cs2_mn,rr_sph_mn,rr_cyl_mn
!
      real :: p,q
      integer :: ics2,irho,j
      integer, pointer :: iglobal_cs2,iglobal_glnTT
      logical :: lheader
!
      p=-density_power_law 
      q=-temperature_power_law
!
      if (lroot) print*,&
           'initial_condition_lnrho: locally isothermal approximation'
      if (lroot) print*,'Radial density stratification with power law=',p
      if (lroot) print*,'Radial temperature stratification with power law=',q
!
!  Set the sound speed - a power law in cylindrical radius.
!
      do m=1,my
        do n=1,mz
          lheader=((m==1).and.(n==1).and.lroot)
          call get_radial_distance(rr_sph,rr_cyl)
          rr=rr_cyl
          call power_law(cs20,rr,-q,cs2,r_ref)
!
!  Store cs2 in one of the free slots of the f-array
!
          nullify(iglobal_cs2)
          call farray_use_global('cs2',iglobal_cs2)
          ics2=iglobal_cs2
          f(:,m,n,ics2)=cs2
          
          nullify(iglobal_glnTT)
          call farray_use_global('glnTT',iglobal_glnTT)
          f(:,m,n,iglobal_glnTT  )=q/rr_sph
          f(:,m,n,iglobal_glnTT+1)=q/rr_sph*cotth(m)
          f(:,m,n,iglobal_glnTT+2)=0.
!
        enddo
      enddo
!
!  Pencilize the density allocation.
!
      do n=1,mz
        do m=1,my
!
          lheader=lroot.and.(m==1).and.(n==1)
!
!  Midplane density
!
          call get_radial_distance(rr_sph,rr_cyl)
          lnrhomid=log(rho0)+p*log(rr_cyl) 
!
!  Vertical stratification
!
          cs2=f(:,m,n,ics2)
          call potential(POT=tmp1,RMN=rr_sph)
          call potential(POT=tmp2,RMN=rr_cyl)
          strat=-(tmp1-tmp2)/cs2
          f(:,m,n,ilnrho) = lnrhomid+strat
!
!  Set the azimuthal velocities
! 
!  Commented out part with the analytical derivation, 
!  since the numerical one (used below) yields better 
!  cancelation.
!
!          call acceleration(g_r)
!          OOK2=max(-g_r/(rr_sph*sinth(m)**3),0.)
!          H2=cs2/OOK2

!          OO2=OOK2* (1 + H2/rr_cyl**2*(p+q) + q*(1-sinth(m)))
!          OO=sqrt(OO2)

!          f(:,m,n,iuz) = rr_cyl*OO
 
        enddo
      enddo
!
!  Use ilnrho as rho - works only for ldensity_nolog (which isn't passed 
!  yet, but, hey, it's MY own custom initial condition file. I know what
!  it does. :-) 
!
      irho=iglobal_glnTT+2  !take an empty slot of f to put irho
      f(:,:,:,irho)=exp(f(:,:,:,ilnrho))
!
!  Azimuthal speed that perfectly balances the pressure gradient. 
!
      do m=m1,m2
        do n=n1,n2

          call grad(f,irho,grho_mn)
          do j=1,3
            glnrho_mn(:,j)=grho_mn(:,j)/f(l1:l2,m,n,irho)
          enddo
          cs2_mn=f(l1:l2,m,n,iglobal_cs2)
          
          glnTT_mn(:,1:2)=f(l1:l2,m,n,iglobal_glnTT:iglobal_glnTT+1)
          
          call get_radial_distance(rr_sph_mn,rr_cyl_mn)
     
          f(l1:l2,m,n,iuz)=sqrt(rr_sph_mn*cs2_mn*&
              (glnrho_mn(:,2)+glnTT_mn(:,2))/cotth(m))
!
        enddo
      enddo
!
! Revert the free slot used for irho to its original value.
!
      f(:,:,:,iglobal_glnTT+2)=0.
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ss
!********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic potential.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine read_initial_condition_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=initial_condition_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=initial_condition_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_initial_condition_init_pars
!***********************************************************************
    subroutine write_initial_condition_init_pars(unit)
!     
      integer, intent(in) :: unit
!
      write(unit,NML=initial_condition_init_pars)
!
    endsubroutine write_initial_condition_init_pars
!***********************************************************************
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
