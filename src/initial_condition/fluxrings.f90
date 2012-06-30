! $Id$
!
!  Magnetic flux rings. Constructed from a canonical ring which is the
!  rotated and translated; see Pencil Code manual, Section C.3.
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
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: fring1=0.0, Iring1=0.0, Rring1=1.0, wr1=0.3
  real :: fring2=0.0, Iring2=0.0, Rring2=1.0, wr2=0.3
  real :: fring3=0.0, Iring3=0.0, Rring3=1.0, wr3=0.3
  real :: radius=0.1, epsilonaa=0.01, widthaa=0.5, x0aa=0.0, z0aa=0.0
  real, dimension(3) :: axisr1=(/0,0,1/), dispr1=(/0.0,0.5,0.0/)
  real, dimension(3) :: axisr2=(/1,0,0/), dispr2=(/0.0,-0.5,0.0/)
  real, dimension(3) :: axisr3=(/1,0,0/), dispr3=(/0.0,-0.5,0.0/)
  integer :: nrings=2
  character (len=labellen) :: fring_profile='tanh'
  real, dimension(ninit) :: amplaa
  character (len=labellen), dimension(ninit) :: initring='nothing'
!
  namelist /initial_condition_pars/  fring1, Iring1, Rring1, wr1, &
       axisr1, dispr1, fring2, Iring2, Rring2, wr2, axisr2, dispr2, &
       fring3, Iring3, Rring3, wr3,  axisr3, dispr3, nrings, radius, &
       fring_profile, amplaa, initring
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
!  Initialize magnetic field.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer :: j
!
      do j=1,ninit
!
        select case (initring(j))
          case ('fluxrings', '4')
            call fluxrings(amplaa(j),f,iaa,iaa,fring_profile)
          case ('fluxrings_WB')
            call fluxrings(amplaa(j),f,iuu,iaa,fring_profile)
          case ('fluxrings_BW')
            call fluxrings(amplaa(j),f,iaa,iuu,fring_profile)
          case ('fluxrings_WW')
            call fluxrings(amplaa(j),f,iuu,iuu,fring_profile)
          case ('nothing')
!
          case default
!
!  Catch unknown values
!
            call fatal_error('initring', &
                 'initring "' // trim(initring(j)) // '" not recognised')
!
        endselect
!
      enddo
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine fluxrings(ampl,f,ivar1,ivar2,prof)
!
!  Initialize the magnetic vector potential.
!
!  Magnetic flux rings. Constructed from a canonical ring which is the
!  rotated and translated (see Pencil Code manual, Section C.3):
!
!    AA(xxx) = D*AA0(D^(-1)*(xxx-xxx_disp)) ,
!
!  where AA0(xxx) is the canonical ring and D the rotation matrix
!  corresponding to a rotation by phi around z, followed by a
!  rotation by theta around y.
!  The array was already initialized to zero before calling this
!  routine.
!
!  Optional argument `profile' allows to choose a different profile (see
!  norm_ring())
!
!  07-may-09/wlad: coded
!
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: tmpv
      real, dimension (nx) :: xx1,yy1,zz1
      real, dimension(3) :: axis,disp
      real :: phi,theta,ct,st,cp,sp
      real :: fring,Iring,R0,width
      integer :: i,ivar3,ivar=-1,ivar1,ivar2
      real :: ampl
      character (len=labellen) :: prof
!
!  fix ivar3=ivar1 (for now)
!
      ivar3=ivar1
!
!  initialize each ring
!
      if (any((/fring1,fring2,Iring1,Iring2/) /= 0.)) then
        ! fringX is the magnetic flux, IringX the current
        if (lroot) then
          print*, 'fluxrings: Initialising magnetic flux rings'
        endif
        do i=1,nrings
          if (i==1) then
            fring = fring1      ! magnetic flux along ring
            Iring = Iring1      ! current along ring (for twisted flux tube)
            R0    = Rring1      ! radius of ring
            width = wr1         ! ring thickness
            axis  = axisr1      ! orientation
            disp  = dispr1      ! position
            ivar  = ivar1
          elseif (i==2) then
            fring = fring2
            Iring = Iring2
            R0    = Rring2
            width = wr2
            axis  = axisr2
            disp  = dispr2
            ivar  = ivar2
          elseif (i==3) then
            fring = fring3
            Iring = Iring3
            R0    = Rring3
            width = wr3
            axis  = axisr3
            disp  = dispr3
            ivar  = ivar3
          else
            call stop_it('fluxrings: nrings is too big')
          endif
          phi   = atan2(axis(2),axis(1)+epsi)
          theta = atan2(sqrt(axis(1)**2+axis(2)**2)+epsi,axis(3))
          ct = cos(theta); st = sin(theta)
          cp = cos(phi)  ; sp = sin(phi)
          ! Calculate D^(-1)*(xxx-disp)
          do n=n1,n2; do m=m1,m2
            xx1= ct*cp*(x(l1:l2)-disp(1))+ct*sp*(y(m)-disp(2))-st*(z(n)-disp(3))
            yy1=-   sp*(x(l1:l2)-disp(1))+   cp*(y(m)-disp(2))
            zz1= st*cp*(x(l1:l2)-disp(1))+st*sp*(y(m)-disp(2))+ct*(z(n)-disp(3))
            call norm_ring(xx1,yy1,zz1,fring,Iring,R0,width,tmpv,PROFILE=prof)
            ! calculate D*tmpv
            f(l1:l2,m,n,ivar  ) = f(l1:l2,m,n,ivar  ) + ampl*( &
                 + ct*cp*tmpv(:,1) - sp*tmpv(:,2) + st*cp*tmpv(:,3))
            f(l1:l2,m,n,ivar+1) = f(l1:l2,m,n,ivar+1) + ampl*( &
                 + ct*sp*tmpv(:,1) + cp*tmpv(:,2) + st*sp*tmpv(:,3))
            f(l1:l2,m,n,ivar+2) = f(l1:l2,m,n,ivar+2) + ampl*( &
                 - st   *tmpv(:,1)                + ct   *tmpv(:,3))
          enddo; enddo
        enddo
      endif
      if (lroot) print*, 'fluxrings: Magnetic flux rings initialized'
!
    endsubroutine fluxrings
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
    subroutine norm_ring(xx1,yy1,zz1,fring,Iring,R0,width,vv,profile)
!
!  Generate vector potential for a flux ring of magnetic flux FRING,
!  current Iring (not correctly normalized), radius R0 and thickness
!  WIDTH in normal orientation (lying in the x-y plane, centred at (0,0,0)).
!
!   1-may-02/wolf: coded (see the manual, Section C.3)
!   7-jun-09/axel: added gaussian and constant (or box) profiles
!
      use Mpicomm, only: stop_it
      use Sub, only: erfunc
!
      real, dimension (nx,3) :: vv
      real, dimension (nx) :: xx1,yy1,zz1,phi,tmp
      real :: fring,Iring,R0,width
      character (len=*) :: profile
!
!  magnetic ring, define r-R
!
      tmp = sqrt(xx1**2+yy1**2)-R0
!
!  choice of different profile functions
!
      select case (profile)
!
!  gaussian profile, exp(-.5*(x/w)^2)/(sqrt(2*pi)*eps),
!  so its derivative is .5*(1.+erf(-x/(sqrt(2)*eps))
!
      case ('gaussian')
        vv(:,3) = - fring * .5*(1.+erfunc(tmp/(sqrt(2.)*width))) &
                          * exp(-.5*(zz1/width)**2)/(sqrt(2.*pi)*width)
!
!  tanh profile, so the delta function is approximated by 1/cosh^2.
!  The name tanh is misleading, because the actual B frofile is
!  1./cosh^2, but this is harder to write.
!
      case ('tanh')
        vv(:,3) = - fring * 0.5*(1+tanh(tmp/width)) &
                          * 0.5/width/cosh(zz1/width)**2
!
!  constant profile, so the delta function is approximated by the function
!  delta(x) = 1/2w, if -w < x < w.
!
      case ('const')
        vv(:,3) = - fring * 0.5*(1.+max(-1.,min(tmp/width,1.))) &
                          * 0.25/width*(1.-sign(1.,abs(zz1)-width))
!
!  there is no default option here
!
      case default
        call stop_it('norm_ring: No such fluxtube profile')
      endselect
!
!  current ring (to twist the B-lines)
!
      tmp = width - sqrt(tmp**2 + zz1**2)
      tmp = Iring*0.5*(1+tanh(tmp/width))     ! Now the A_phi component
      phi = atan2(yy1,xx1)
      vv(:,1) = - tmp*sin(phi)
      vv(:,2) =   tmp*cos(phi)
!
    endsubroutine norm_ring
!********************************************************************
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
