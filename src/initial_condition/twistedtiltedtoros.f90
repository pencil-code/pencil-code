! $Id: bipole.f90
! 
!
!  17-jul-12/piyali: adapted from fluxrings.f90 but for spherical
!  geometry. Also some major changes, hence a new initial condition
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
!!  integer :: dummy
!
  real :: fring=1e-3,r0=0.2,Iring=0.0,rho0,tilt=0.0,width=0.02, &
          posx=0.98,dIring=0.0
  real :: scale_aa=1.0
  character (len=8) :: extfieldfile
  logical :: luse_only_botbdry=.false.
!
  namelist /initial_condition_pars/ fring, r0, Iring, posx,&
  tilt,extfieldfile,rho0,width,luse_only_botbdry,scale_aa, &
  dIring
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
         "$Id: fieldloop.f90 17229 2011-07-19 20:14:00Z sven.bingert $")
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
    subroutine initial_condition_lnrho(f)
!
!  Initialize the magnetic vector potential.
!
!  02-may-12/piyali: coded
!
      use Sub, only: curl
      use IO,  only: input_snap,input_snap_finalize
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx,3) :: bb
      real :: xi
      integer :: l
!
      call initial_condition_aa(f) 
      do m=m1,m2
        do n=n1,n2
          call curl(f,iaa,bb)
          f(l1:l2,m,n,ilnrho)=log(rho0+bb(:,1)**2+bb(:,2)**2+bb(:,3)**2)
        end do
      end do
    endsubroutine initial_condition_lnrho  
!
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  02-may-12/piyali: coded
!
      use Boundcond, only: update_ghosts 
      use Deriv, only: der
      use Hydro, only: find_umax
      use Sub, only: curl_mn,gij
      use IO,  only: input_snap,input_snap_finalize
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx,3,3) :: aij
      real, dimension (nx,3) :: tmpv,aa,bb
      real, dimension (nx) :: xx0,yy0,zz0,xx1,yy1,zz1,dist,distxy
      real, dimension (nx) :: tmpx,tmpy,tmpz,dfy,dfz
      real :: xi,umax
      integer :: l
!
!  IMPLEMENTATION OF INSERTION OF BIPOLES (Non-Potential part) 
!  (Yeates, Mackay and van Ballegooijen 2008, Sol. Phys, 247, 103)
!
      do n=n1,n2
        do m=m1,m2
          do l=l1,l2
!
            if (lcartesian_coords) then 
              xi=((x(l)**2+z(n)**2)/2.+y(m)**2)/r0**2
              f(l,m,n,iax) = fring*Iring*z(n)*exp(-2*xi)
              f(l,m,n,iay) = fring*r0*exp(-xi)
              f(l,m,n,iaz) = -fring*Iring*x(l)*exp(-2*xi)
            else if (lcylindrical_coords) then 
              call fatal_error('initial_condition_aa','Bipoles not coded  &
              for cylindrical coordinates')
            endif
          enddo
!
!  Then set up the helical field
!
          if (lspherical_coords) then 
              xx0=x(l1:l2)*sinth(m)*cos(z(n))
              yy0=x(l1:l2)*sinth(m)*sin(z(n))
              zz0=x(l1:l2)*costh(m)
          ! Calculate D^(-1)*(xxx-disp)
              xx1=xx0-posx
              yy1=cos(tilt*pi/180.0)*yy0+sin(tilt*pi/180.0)*zz0
              zz1=-sin(tilt*pi/180.0)*yy0+cos(tilt*pi/180.0)*zz0
              call norm_ring(xx1,yy1,zz1,fring,Iring,r0,width,posx,tmpv,PROFILE='gaussian')
            ! calculate D*tmpv
              tmpx=tmpv(:,1)
              tmpy=cos(tilt*pi/180.0)*tmpv(:,2)-sin(tilt*pi/180.0)*tmpv(:,3)
              tmpz=sin(tilt*pi/180.0)*tmpv(:,2)+cos(tilt*pi/180.0)*tmpv(:,3)
              dist=sqrt(xx0**2+yy0**2+zz0**2)
              distxy=sqrt(xx0**2+yy0**2)
              f(l1:l2,m,n,iax) =  scale_aa*f(l1:l2,m,n,iax)+& 
                     xx0*tmpx/dist+yy0*tmpy/dist+zz0*tmpz/dist
              f(l1:l2,m,n,iay) = scale_aa*f(l1:l2,m,n,iay)+ &
                     xx0*zz0*tmpx/(dist*distxy)+yy0*zz0*tmpy/(dist*distxy) &
                      -distxy*tmpz/dist
              f(l1:l2,m,n,iaz) = scale_aa*f(l1:l2,m,n,iaz)-yy0*tmpx/distxy+xx0*tmpy/distxy               
          endif
        enddo
      enddo
      call update_ghosts(f)
!
! The arrays iux,iuy,iuz is used to store the surface radial
! magnetic field, theta-comp of vel and phi-comp of velocity at the surface
! respectively if required by the problem
!
        do n=n1,n2
          do m=m1,m2
            aa=f(l1:l2,m,n,iax:iaz)
            call gij(f,iaa,aij,1)
            call curl_mn(aij,bb,aa)
            f(l1:l2,m,n,iux)=bb(:,1)**2
          enddo
        enddo  
        call find_umax(f,umax) 
        call update_ghosts(f)
        do n=n1,n2
          do m=m1,m2
            call der(f,iux,dfz,3)
            call der(f,iux,dfy,2)
!            
! Normalize to unity.          
!
            f(l1,m,n,iuy)=dIring*dfz(1)/umax
            f(l1,m,n,iuz)=-dIring*dfy(1)/umax
          enddo
        enddo  
      f(l1:l2,m1:m2,n1:n2,iux)=0.0

      if (luse_only_botbdry) then
        f(l1+1:l2,:,:,:)=0.0
      endif
!
    endsubroutine initial_condition_aa
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
    subroutine norm_ring(xx1,yy1,zz1,fring,Iring,r0,width,posx,vv,profile)
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
      real, dimension (nx) :: xx1,yy1,zz1,phi,tmp,cs,sn
      real :: fring,Iring,r0,width,posx,width2
      character (len=*) :: profile
      integer :: l
!
!  magnetic ring, define r-R
!
      tmp = sqrt(xx1**2+yy1**2)-r0
!
!  choice of different profile functions
!
      select case (profile)
!
!  gaussian profile, exp(-.5*(x/w)^2)/(sqrt(2*pi)*eps),
!  so its derivative is .5*(1.+erf(-x/(sqrt(2)*eps))
!
      case ('gaussian')
!        width2=20.0*width
!        do l=1, nx
!        if (abs(zz1(l)).le.width2) then
!          if (sqrt(tmp(l)**2+zz1(l)**2).le.width2) then
!            vv(l,3) = + fring * .5*(erfunc(sqrt(width2**2-zz1(l)**2)/width/sqrt(2.))- &
!                      erfunc(tmp(l)/(sqrt(2.)*width))) &
!                     *exp(-.5*(zz1(l)/width)**2)/(sqrt(2.*pi)*width)
!          else 
!            if (tmp(l).le.0) then
!              vv(l,3) = + fring *(erfunc(sqrt(width2**2-zz1(l)**2)/width/sqrt(2.))) &
!                     *exp(-.5*(zz1(l)/width)**2)/(sqrt(2.*pi)*width)
!            else
!              vv(l,3)=0.0
!            endif         
!          endif
!        else
!          vv(l,3)=0.0
!        endif  
!        enddo
          vv(:,3) = + fring * .5*(1.-erfunc(tmp/(sqrt(2.)*width))) &
                   *exp(-.5*(zz1/width)**2)/(sqrt(2.*pi)*width)
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
      case default
        call stop_it('norm_ring: No such fluxtube profile')
      endselect
!
!  current ring (to twist the B-lines)
!
      phi = atan2(yy1,xx1)
      tmp = sqrt((xx1-r0*cos(phi))**2 + (yy1-r0*sin(phi))**2+zz1**2)
      tmp = Iring*fring*exp(-0.5*(tmp/width)**2)/(sqrt(2.*pi)*width)    ! Now the A_phi component
      vv(:,1) = - tmp*sin(phi)
      vv(:,2) =   tmp*cos(phi)
!
    endsubroutine norm_ring
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
