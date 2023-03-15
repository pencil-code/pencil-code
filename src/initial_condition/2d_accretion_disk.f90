! $Id$
!
!  21-feb-23/hongzhe: adapted from initial_condition/centrifugal_balance.f90
! 
!  This module sets up the initial condition for an axisymmetric accretion
!  disk in spherical coordinates, using the following structure:
!  (denoting the radial distance in spherical coordinates in r)
!  (and that in cylindrical coordinates in R)
!
!  gravity:  potential propto r**(-1)
!  
!  In the R direction:
!  velocity: keplerian R**(-0.5), plus noise propto cs, i.e. also R**(-0.5)
!  sound speed: cs propto R**(-0.5)
!  density: propto R**(-1.5)
!  pressure: propto R**(-2.5)
!
!  In the z direction:
!  All the quantities depend on a single variable l=z/R=1/tan(theta) and
!  have some analytical form.
!  
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
! global variables to this module
!
  integer :: mcrit  !  m<mcrit is corona region
  real, dimension (my) :: frho=1.,fu=1.
!
  character (len=labellen) :: inituu='0',initlnrho='uniform'
  real :: uukep_inner=1.,ampluu_noise=0.
  real :: rmin_lnrho=1., rmax_lnrho=10.
  logical :: lsimple_density=.false., luniform_density=.false.
!
  namelist /initial_condition_pars/ uukep_inner, &
       ampluu_noise,lsimple_density,luniform_density,&
       inituu,initlnrho,rmin_lnrho,rmax_lnrho
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
         "$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  21-feb-23/hongzhe: coded
!
      use EquationOfState, only: gamma,cs20
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: tmp
      integer :: mm
!
      do m=1,my
        tmp = sin(y(m))-1.+cs20*gamma/(gamma-1.)
        if (tmp<=0.) then
          frho(m) = 0.001
          mcrit = m+1
        else
          frho(m) = ( tmp*(gamma-1.)/cs20/gamma )**(1./(gamma-1.))
          if (frho(m)<=0.001) then
            frho(m) = 0.001
            mcrit = m+1
          endif
        endif
        tmp = sin(y(m))-cs20*gamma/(gamma-1.)*frho(m)**(gamma-1.)
        fu(m) = sqrt( max(0.,tmp) )
      enddo
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  21-feb-23/hongzhe: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: fact
!
      select case (inituu)
      case ('0')
        f(:,:,:,iux:iuz) = 0.
      case ('kep-cyl')  !  Kepelerian velocity using cylindrical radius
        do m=1,my
        do n=1,mz
          f(:,m,n,iux) = f(:,m,n,iux) + 0.
          f(:,m,n,iuy) = f(:,m,n,iuy) + 0.
          fact = (x*sin(y(m)))**(-0.5)*fu(m)
          f(:,m,n,iuz) = f(:,m,n,iuz) + uukep_inner*fact
          call gaunoise_vect(ampluu_noise*fact,f,iux,iuz)
        enddo
        enddo
      case ('kep-cyl-noz')  !  same as 'kep_cyl' but no vertical structure
        do m=1,my
        do n=1,mz
          f(:,m,n,iux) = f(:,m,n,iux) + 0.
          f(:,m,n,iuy) = f(:,m,n,iuy) + 0.
          fact = (x*sin(y(m)))**(-0.5)
          f(:,m,n,iuz) = f(:,m,n,iuz) + uukep_inner*fact
          call gaunoise_vect(ampluu_noise*fact,f,iux,iuz)
        enddo
        enddo
      case default
        if (lroot) print*,'initial_condition_uu: No profile of inituu="'//trim(inituu)//'"'
      endselect
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine gaunoise_vect(ampl,f,i1,i2)
!
!  Add Gaussian noise (= normally distributed) white noise for variables i1:i2
!
!  23-may-02/axel: coded
!  10-sep-03/axel: result only *added* to whatever f array had before
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i1,i2
!
      real, dimension (nx) :: r,p,tmp,ampl
      integer :: i,j
!
      intent(in)    :: ampl,i1,i2
      intent(inout) :: f
      integer, save, dimension(mseed) :: rstate=0
!
!  set gaussian random noise vector
!
      do i=i1,i2
        if (lroot.and.m==1.and.n==1) print*,'gaunoise_vect: variable i=',i
        if (modulo(i-i1,2)==0) then
           do j=1,nx
              r(j)=ran0(rstate(1))
              p(j)=ran0(rstate(1))
          enddo
          tmp=sqrt(-2*log(r))*sin(2*pi*p)
        else
          tmp=sqrt(-2*log(r))*cos(2*pi*p)
        endif
        f(l1:l2,m,n,i)=f(l1:l2,m,n,i)+ampl*tmp
      enddo
!
    endsubroutine gaunoise_vect
!***********************************************************************
    function ran0(dummy)
!
!  The 'Minimal Standard' random number generator
!  by Lewis, Goodman and Miller.
!
!  28.08.02/nils: Adapted from Numerical Recipes
!
      integer, intent(inout) :: dummy
!
      integer :: k
      integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836, &
           mask=123459876
      real, parameter :: am=1./im
      real :: ran0
!
      dummy=ieor(dummy,mask)
      k=dummy/iq
      dummy=ia*(dummy-k*iq)-ir*k
      if (dummy<0) dummy=dummy+im
      ran0=am*dummy
      dummy=ieor(dummy,mask)
!
    endfunction ran0
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  21-feb-23/hongzhe: adapted
!
      use FArrayManager
      use EquationOfState, only: gamma,lnrho0,cs20
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: fact
!
      select case (initlnrho)
      case ('uniform')
        f(:,:,:,ilnrho) = lnrho0 + log(0.001)
      case ('hydrostatic')
        do n=1,mz
        do m=1,mcrit-1  !  corona
          fact=x**(-1/(gamma-1.))*frho(m)
          f(:,m,n,ilnrho) = f(:,m,n,ilnrho) + lnrho0 + log(fact)
        enddo
        do m=mcrit,my  !  disk
          fact=(x*sin(y(m)))**(-1/(gamma-1.))*frho(m)
          f(:,m,n,ilnrho) = f(:,m,n,ilnrho) + lnrho0 + log(fact)
        enddo
        enddo
      case ('thin-disk')  !  powerlaw in r_sph and exponential in z
        do m=1,my
        do n=1,mz
          f(:,m,n,ilnrho) = f(:,m,n,ilnrho) - 1.5*log(x) - (cos(y(m)))**2/cs20
        enddo
        enddo
      case ('thin-disk-gaur')  !  powerlaw+gaussian in r_sph and exponential in z
        fact = (x-(rmin_lnrho+rmax_lnrho)/2.)/2./(rmax_lnrho-rmin_lnrho)*5.
        do m=1,my
        do n=1,mz
          f(:,m,n,ilnrho) = f(:,m,n,ilnrho) - 1.5*log(x) - fact**2 - (cos(y(m))/2.)**2/cs20
        enddo
        enddo
      case default
        if (lroot) print*,'initial_condition_uu: No profile of initlnrho="'//trim(initlnrho)//'"'
      endselect
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
!***********************************************************************
    subroutine initial_condition_aa(f)
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
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat)
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=initial_condition_pars)
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
