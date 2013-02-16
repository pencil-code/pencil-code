!  $Id$
!
!  Initial condition (density, magnetic field, velocity) 
!  for magnetohydrostatical equilibrium in a global accretion
!  disk with an imposed (cylindrically symmetric) sound speed 
!  profile in spherical coordinates. 
!
!  07-may-09/wlad: adapted from noinitial_condition.f90
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
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: density_power_law,temperature_power_law,plasma_beta
  logical :: lnumerical_mhsequilibrium=.true.
  logical :: lintegrate_potential=.true.
!
  namelist /initial_condition_pars/ &
      density_power_law,temperature_power_law,plasma_beta,&
      lnumerical_mhsequilibrium,lintegrate_potential
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  07-oct-09/wlad: coded
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      use Gravity,       only: acceleration
      use Sub,           only: get_radial_distance
      use FArrayManager, only: farray_use_global

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: rr_sph,rr_cyl,g_r
      real, dimension (mx) :: OOK2,OO2,tmp1,H2
      real :: p,q,tmp2,ksi
      integer, pointer :: iglobal_cs2
      integer :: ics2
!
!  Set the sound speed
!
      call set_sound_speed(f)
!
!  Get the sound speed
!
      nullify(iglobal_cs2)
      call farray_use_global('cs2',iglobal_cs2)
      ics2=iglobal_cs2 
!
!  Analytical expression that leads to an equilibrium configuration. 
!  Commented out because the numerical equilibrium enforced below 
!  is a lot better when it comes to reproduce what the code actually
!  solves. 
!
      if (.not.lnumerical_mhsequilibrium) then 
!
        if (lmagnetic) then 
          ksi=(1.+plasma_beta)/plasma_beta
        else
          ksi=1.
        endif
!
        p=-density_power_law
        q=-temperature_power_law
!
        do m=1,my;do n=1,mz
          call acceleration(g_r)
          call get_radial_distance(rr_sph,rr_cyl)
          OOK2=max(-g_r/(rr_sph*sinth(m)**3),0.)
!
          H2=f(:,m,n,ics2)/OOK2
!
          tmp1=H2/rr_cyl**2*(ksi*(p+q-2.) + 2.)
          tmp2=ksi*(q+1.)*(1.-sinth(m))
!
          OO2=OOK2*(tmp1+tmp2+sinth(m))
!
          f(:,m,n,iuz) = f(:,m,n,iuz) + rr_cyl*sqrt(OO2)
        enddo;enddo
!
      endif
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
      use EquationOfState, only: rho0
      use FArrayManager,   only: farray_use_global
      use Gravity,         only: potential
      use Sub,             only: get_radial_distance
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: rr_sph,rr_cyl
      real, dimension (mx) :: lnrhomid,strat
      real, dimension (mx) :: cs2,tmp1,tmp2
      real    :: p,q
      integer :: ics2
      integer, pointer :: iglobal_cs2
!
      p=-density_power_law 
      q=-temperature_power_law
!
      if (lroot) print*,&
           'initial_condition_lnrho: locally isothermal approximation'
      if (lroot) print*,'Radial density stratification with power law=',p
      if (lroot) print*,'Radial temperature stratification with power law=',q
!
!  Get the sound speed globals
!
      nullify(iglobal_cs2)
      call farray_use_global('cs2',iglobal_cs2)
      ics2=iglobal_cs2 
!
!  Pencilize the density allocation.
!
      do n=1,mz
        do m=1,my
!
!  Midplane density
!
          call get_radial_distance(rr_sph,rr_cyl)
          lnrhomid=log(rho0)+p*log(rr_cyl) 
          f(:,m,n,ilnrho) = f(:,m,n,ilnrho)+lnrhomid
!
!  Vertical stratification
!
          cs2=f(:,m,n,ics2)
          call potential(POT=tmp1,RMN=rr_sph)
          call potential(POT=tmp2,RMN=rr_cyl)
          strat=-(tmp1-tmp2)/cs2
          f(:,m,n,ilnrho) = f(:,m,n,ilnrho)+strat
!
        enddo
      enddo
!
!  If the run is non-magnetic, this is the last routine called. 
!  Enforce numerical equilibrium then
!
      if ((.not.lmagnetic).and.lnumerical_mhsequilibrium) &
          call enforce_numerical_equilibrium(f,lhd=.true.)
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential. Constant plasma 
!  beta magnetic field. 
!
!  07-may-09/wlad: coded
!
      use FArrayManager, only: farray_use_global
      use Sub, only: gij,curl_mn,get_radial_distance
      use EquationOfState, only: cs20,rho0
      use Mpicomm, only: mpibcast_real,mpisend_real,mpirecv_real
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(nx) :: pressure,Bphi,Atheta
      real, dimension(mx) :: rr_sph,rr_cyl
      real, dimension(nx) :: tmp
      real, dimension(0:nx) :: tmp2
      real :: dr,psum
      integer, pointer :: iglobal_cs2
      integer :: i,ics2,irho,iprocx,iprocy,iserial_xy
      real, dimension(ny,nprocx*nprocy) :: procsum
      real, dimension(ny) :: procsum_loc,tmpy
!
!  Get the sound speed globals 
!
      nullify(iglobal_cs2)
      call farray_use_global('cs2',iglobal_cs2)
      ics2=iglobal_cs2 
!
!  Density is already in linear after init_lnrho, so
!
      irho=ilnrho
!
      do m=m1,m2
!
        pressure=f(l1:l2,m,npoint,irho)*f(l1:l2,m,npoint,ics2)
!
!  The following line assumes mu0=1
!
        Bphi = sqrt(2*pressure/plasma_beta)
!
!  Bphi = 1/r*d/dr(r*Atheta), so integrate: Atheta=1/r*Int(B*r)dr
!
        call get_radial_distance(rr_sph,rr_cyl)
        !dr=rr_sph(2)-rr_sph(1)
        tmp=Bphi*rr_sph(l1:l2)
!
        iserial_xy=ipx+nprocx*ipy+1
        tmp2(0)=0.
        do i=1,nx
          dr=rr_sph(i+l1-1)-rr_sph(i+l1-2) 
          tmp2(i)=tmp2(i-1)+tmp(i)*dr
        enddo
        procsum_loc(m-m1+1)=tmp2(nx)
      enddo
!
      if (ip<=9) print*,'send',ipx+nprocx*ipy+1,procsum_loc(mpoint)
!
      if (lroot) then 
        procsum(:,1)=procsum_loc
        do iprocx=0,nprocx-1; do iprocy=0,nprocy-1
          iserial_xy=iprocx+nprocx*iprocy+1
          if (iserial_xy/=1) then
            call mpirecv_real(tmpy,ny,iserial_xy-1,111)
            procsum(:,iserial_xy)=tmpy
          endif
          if (ip<=9) print*,'recv',iserial_xy,procsum(mpoint,iserial_xy)
        enddo; enddo
      else
        call mpisend_real(procsum_loc,ny,0,111)
      endif
!
      call mpibcast_real(procsum,(/nprocx*nprocy,ny/))
!
      do m=m1,m2;do n=n1,n2
!
        pressure=f(l1:l2,m,n,irho)*f(l1:l2,m,n,ics2)
!
!  The following line assumes mu0=1
!
        Bphi = sqrt(2*pressure/plasma_beta)
!
!  Bphi = 1/r*d/dr(r*Atheta), so integrate: Atheta=1/r*Int(B*r)dr
!
        call get_radial_distance(rr_sph,rr_cyl)
        !dr=rr_sph(2)-rr_sph(1)
        tmp=Bphi*rr_sph(l1:l2)
!
        if (nprocx==1) then 
          tmp2(0)=0.
          do i=1,nx
            dr=rr_sph(i+l1-1)-rr_sph(i+l1-2)
            tmp2(i)=tmp2(i-1)+tmp(i)*dr
          enddo
        else
          if (lfirst_proc_x) then 
            tmp2(0)=0.
            do i=1,nx
              dr=rr_sph(i+l1-1)-rr_sph(i+l1-2) 
              tmp2(i)=tmp2(i-1)+tmp(i)*dr
            enddo
          else 
            psum=0.
            do iprocx=0,ipx-1
              iserial_xy=iprocx+nprocx*ipy+1
              psum=psum+procsum(m-m1+1,iserial_xy)
            enddo
            tmp2(0)=psum
            do i=1,nx
              dr=rr_sph(i+l1-1)-rr_sph(i+l1-2) 
              tmp2(i)=tmp2(i-1)+tmp(i)*dr
            enddo
          endif
        endif
        Atheta=tmp2(1:nx)/rr_sph(l1:l2)
        f(l1:l2,m,n,iay)=f(l1:l2,m,n,iay)+Atheta
      enddo;enddo
!
!  All quantities are set. Enforce numerical equilibrium.
!
    if (lnumerical_mhsequilibrium) &
         call enforce_numerical_equilibrium(f,lhd=.false.)
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine set_sound_speed(f)
!
!  Set the thermo-related quantities. Illustrates that 
!  the user can define as many internal routines as wanted.
!
!  10-may-09/wlad : moved from initial_condition_lnrho
!
      use FArrayManager,   only: farray_use_global
      use EquationOfState, only: cs20
      use Sub,             only: get_radial_distance
!
      real, dimension (mx,my,mz,mfarray) :: f      
      real, dimension (mx) :: rr_sph,rr_cyl,cs2
      integer, pointer :: iglobal_cs2,iglobal_glnTT
      real    :: q
      integer :: ics2,iglnTT
!
      q=-temperature_power_law
!
!  Get the globals needed to store sound speed and temperature gradient
!
      nullify(iglobal_cs2)
      call farray_use_global('cs2',iglobal_cs2);ics2=iglobal_cs2
      nullify(iglobal_glnTT)
      call farray_use_global('glnTT',iglobal_glnTT);iglnTT=iglobal_glnTT
!
!  Set the sound speed - a power law in cylindrical radius.
!
      do m=1,my
        do n=1,mz
          call get_radial_distance(rr_sph,rr_cyl)
          cs2=cs20*rr_cyl**q
!
!  Store cs2 in one of the free slots of the f-array
!
          f(:,m,n,ics2)=cs2
!
!  Same for the temperature gradient
!          
          f(:,m,n,iglnTT  )=q/rr_sph
          f(:,m,n,iglnTT+1)=q/rr_sph*cotth(m)
          f(:,m,n,iglnTT+2)=0.
!
        enddo
      enddo
!
    endsubroutine set_sound_speed
!***********************************************************************
    subroutine enforce_numerical_equilibrium(f,lhd)
!
!  This subroutine does exactly what's done in runtime to the 
!  velocity field, in order to ensure precise numerical 
!  magnetohydrostatical equilibrium. 
!
      use Sub, only: grad,get_radial_distance,&
                     gij,curl_mn,gij_etc,cross_mn,multsv_mn
      use FArrayManager, only: farray_use_global
!      use BoundCond, only: update_ghosts
!
      real, dimension (mx,my,mz,mfarray) :: f      
      real, dimension (nx,3,3) :: aij,bij
      real, dimension (nx,3) :: glnTT,grho,glnrho
      real, dimension (nx,3) :: aa,bb,jj,jxb,jxbr
      real, dimension (nx) :: cs2,rr_sph,rr_cyl,rho1
      real, dimension (nx) :: fpres_thermal,fpres_magnetic
      integer, pointer :: iglobal_cs2, iglobal_glnTT
      integer :: j,irho,ics2,iglnTT
      logical :: lhd
!
!  Get the temperature globals
!
      nullify(iglobal_cs2)
      call farray_use_global('cs2',iglobal_cs2);ics2=iglobal_cs2
      nullify(iglobal_glnTT)
      call farray_use_global('glnTT',iglobal_glnTT);iglnTT=iglobal_glnTT
!
!  Use ilnrho as rho - works only for ldensity_nolog (which isn't passed 
!  yet, but, hey, it's MY own custom initial condition file. I know what
!  it does. :-) 
!
      if (lhd) then 
        irho=iglnTT+2  !take an empty slot of f to put irho
        f(l1:l2,m1:m2,n1:n2,irho)=exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
      else
        irho=ilnrho
      endif
!
      !call update_ghosts(f)
!
!  Azimuthal speed that perfectly balances the pressure gradient. 
!
      do m=m1,m2
        do n=n1,n2
!
          call grad(f,irho,grho)
          do j=1,3
            glnrho(:,j)=grho(:,j)/f(l1:l2,m,n,irho)
          enddo
          cs2=f(l1:l2,m,n,ics2) 
          glnTT(:,1:2)=f(l1:l2,m,n,iglnTT:iglnTT+1)
!
          fpres_thermal=cs2*(glnrho(:,2)+glnTT(:,2))
!
          if (lmagnetic) then
            aa=f(l1:l2,m,n,iax:iaz)         !aa
            call gij(f,iaa,aij,1)           !aij
            call curl_mn(aij,bb,aa)         !bb
            call gij_etc(f,iaa,aa,aij,bij)  !bij
            call curl_mn(bij,jj,bb)         !jj
            call cross_mn(jj,bb,jxb)        !jxb
            rho1=1./f(l1:l2,m,n,irho)       !1/rho
            call multsv_mn(rho1,jxb,jxbr)   !jxb/rho
            fpres_magnetic=-jxbr(:,2)
          else
            fpres_magnetic=0.
          endif
!         
          call get_radial_distance(rr_sph,rr_cyl)
     
          f(l1:l2,m,n,iuz)=f(l1:l2,m,n,iuz)+&
              sqrt(rr_sph*(fpres_thermal+fpres_magnetic)/cotth(m))
!
        enddo
      enddo
!
!  Revert the free slot used for irho to its original value.
!
      if (lhd) f(:,:,:,iglnTT+2)=0.
!
    endsubroutine enforce_numerical_equilibrium
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
