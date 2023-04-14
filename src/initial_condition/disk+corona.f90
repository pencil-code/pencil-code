! $Id$
!
!  This module provide a way for users to specify custom initial
!  conditions.
!
!  The module provides a set of standard hooks into the Pencil Code
!  and currently allows the following customizations:
!
!   Description                               | Relevant function call
!  ------------------------------------------------------------------------
!   Initial condition registration            | register_initial_condition
!     (pre parameter read)                    |
!   Initial condition initialization          | initialize_initial_condition
!     (post parameter read)                   |
!                                             |
!   Initial condition for momentum            | initial_condition_uu
!   Initial condition for density             | initial_condition_lnrho
!   Initial condition for entropy             | initial_condition_ss
!   Initial condition for magnetic potential  | initial_condition_aa
!                                             |
!   Initial condition for all in one call     | initial_condition_all
!     (called last)                           |
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
  real :: l0_d, l_d1, psipn0, psipn1, psi0, psi1, lnrho0_c
  real :: m0, r0_d, rs=1.0, apara, h_d, dsteep, ngamma, &
          rho0_c, cs0_c, Tc, rpos,psi0d,pbeta,amplaa=1.0
  character (len=labellen) :: initaa_dc='nothing'
!
  namelist /initial_condition_pars/ r0_d, apara,amplaa, & 
  h_d, rho0_c, dsteep, ngamma, cs0_c, Tc, rpos, pbeta,initaa_dc
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

!  /mayank
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!      
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uu
!***********************************************************************

    subroutine initial_condition_all(f,profiles)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  07-may-09/wlad: coded

!/mayank

      use Mpicomm, only: mpibcast
      use EquationOfState, only: get_cp1, cs0, cs20, cs2bot, cs2top, rho0, lnrho0, &
                             gamma, gamma1, gamma_m1
      use FArrayManager,   only: farray_use_global
      use Sub,             only: get_radial_distance, location_in_proc
      use SharedVariables, only: get_shared_variable
      
      real, dimension (mx,my,mz,mfarray), optional, intent(inout):: f
      real, dimension (:,:),              optional, intent(out)  :: profiles
      real, dimension (mx, my, mz) :: psi, rcut, masked, xmesh, tmasked, masked3
      real, dimension (mx) :: rr_sph, rr_cyl
      real, dimension (nx) :: rho_d, lnrho_c, psipn, l_d, vphi_d, lnT_d, aphi
      real, dimension (nx) :: dzphi_d, drphi_d, pg, dzrho_d, drrho_d, kphi, pg_c, pbeta_c
      real, pointer :: g0
      real :: cp1, m00, rpn=1.2, pg0, drphi_d0, drrho_d0
      integer :: lpos, mpos, npos, procno,procnomax

! Some constants

      call get_shared_variable('g0',g0)
      lnrho0_c=log(rho0_c)
      rs=2*g0/c_light**2
      l0_d=sqrt(g0*r0_d**3.0)/(r0_d-rs)
      l_d1=l0_d*(2.0*rs/r0_d)**apara
      psipn0=-g0/(r0_d-rs)
      psipn1=-g0/rs
      psi0d=psipn0+1.0/(2.0-2.0*apara)*(l0_d/r0_d)**2.0
      psi1=psipn1+1.0/(2.0-2.0*apara)*(l_d1/(2.0*rs))**2.0
!    
      call get_cp1(cp1) 
      m00=(g0/G_Newton)
      do n=n1,n2
        do m=1,my
          xmesh(:,m,n)=x
          call get_radial_distance(rr_sph,rr_cyl)
          rcut(:,m,n)=rr_sph          
! Specific angular momentum
!        
          l_d=l0_d*(rr_cyl(l1:l2)/r0_d)**apara
!       
! Pseudo-Newtonian potential
!        
        psipn=-g0/(rr_sph(l1:l2)-rs)

! Gravitational Potential for the system
        
        psi(l1:l2,m,n)=psipn+1.0/(2.0-2.0*apara)*(l_d/rr_cyl(l1:l2))**2.0
        enddo
      enddo
      if (location_in_proc((/rpos,xyz0(2), 0.0/),lpos,mpos,npos)) then
          psi0=psi(lpos,mpos,npos)
          procno=iproc
          print*, psi0, x(lpos),procno
      else
          procno=-100 
      endif
!
! The following is a hack as procno is not set for other processors
!
      call find_procno(procno,procnomax)
      call mpibcast(psi0,procnomax)
      print*, iproc, procnomax, psi0
      mask: where ((-(psi-psi0) .gt. 0.0) .and. (abs(xmesh) .gt. rpos))
               masked=1.0
             elsewhere
               masked=0.0
            end where mask
!
      mask2: where ((-(psi-psi0) .gt. 0.0) .and. (abs(xmesh) .gt. rpos))
               tmasked=0.0
             elsewhere
               tmasked=1.0
            end where mask2
      mask3: where (rcut .gt. rpn)
               masked3=1.0
            elsewhere
               masked3=0.0
            end where mask3
!
      do n=n1,n2
        do m=m1,m2
          call get_radial_distance(rr_sph,rr_cyl)
       
! Specific angular momentum
!
          l_d=l0_d*(rr_cyl(l1:l2)/r0_d)**apara
!
! Pseudo-Newtonian potential
!
        psipn=-g0/(rr_sph(l1:l2)-rs)
        
!  Disk Density
                  
        rho_d=masked(l1:l2,m,n)*exp(lnrho0)*(1.0-gamma*(psi(l1:l2,m,n)-psi0)/(cs0**2*(ngamma+1.0)))**ngamma
!
! Corona Density
        lnrho_c=lnrho0_c-masked3(l1:l2,m,n)*(psipn-psipn1)*gamma/cs0_c**2
        
        
        pg=exp(lnrho0)**(-1.0/ngamma)*cs0**2.0/gamma*(rho_d)**(1.0+1.0/ngamma)
        pg_c=exp(lnrho_c)*cs0_c**2.0/gamma
        f(l1:l2,m,n,ilnrho) = log(rho_d+exp(lnrho_c))
        vphi_d=l_d*x(l1:l2)/rr_cyl(l1:l2)**2*masked(l1:l2,m,n)
        f(l1:l2,m,n,iuy) = vphi_d  
        lnT_d=(log(rho_d)-lnrho0)/ngamma+log(cs0**2.0/(gamma_m1/cp1))
        f(l1:l2,m,n,ilnTT)=log(exp(lnT_d)+Tc*tmasked(l1:l2,m,n))
        enddo
      enddo
      select case (initaa_dc)
        case ('kato')
          do n=n1,n2
            do m=m1,m2
              call get_radial_distance(rr_sph,rr_cyl)
        
              drphi_d0=g0/((r0_d*rs-rs)**2+(1e-3*rs)**2)*40.0*rs/sqrt((r0_d*rs)**2.0+&
                     (1e-3*rs)**2)-l0_d**2.0/((1-apara)*sqrt((r0_d*rs)**2+(1e-3*rs)**2)**3.0)
!
              pg0=exp(lnrho0)**(-1.0/ngamma)*cs0**2.0/gamma*(exp(lnrho0))**(1.0+1.0/ngamma)
!
              drrho_d0=-gamma*drphi_d0/(cs0**2.0*(ngamma+1.0))*exp(lnrho0)*ngamma*(1.0/ &
                  (cs0**2*(ngamma+1.0)))**(ngamma-1.0)
!        
              kphi=sqrt(2.0*mu0*pg0/(pbeta*(drrho_d0)**2.0))
              f(l1:l2,m,n,iay)=amplaa*x(l1:l2)*exp(f(l1:l2,m,n,ilnrho)) !*masked(l1:l2,m,n)
            enddo
          enddo
        case ('loop')
          do n=n1,n2
            do m=m1,m2
              call get_radial_distance(rr_sph,rr_cyl)
              aphi=amplaa*rpos**2.0/(rpos**2.0+rr_sph(l1:l2)**2)**1.5*(1+15.0*rpos**2.0*&
                x(l1:l2)**2.0/(8.0*(rpos**2.0+rr_sph(l1:l2)**2)**2.0))
              f(l1:l2,m,n,iay)=aphi*x(l1:l2)*masked(l1:l2,m,n)
            enddo
          enddo
        case ('nothing')
!         Do nothing
        case default
          call fatal_error('initialize_magnetic','No such initaa_dc')
      endselect 
!
!
    endsubroutine initial_condition_all
!***********************************************************************
    subroutine find_procno(procno,procnomax)
!
!  Find the absolute maximum of the velocity.
!
!  13-sep-2023/piyali: adapted from find_umax
!
      use Mpicomm, only: mpiallreduce_max, MPI_COMM_WORLD
!
      integer, intent(in) :: procno
      integer, intent(out) :: procnomax
!
      integer :: procno1
!
!  Find the maximum.
!
      procno1 = procno
      call mpiallreduce_max(procno1, procnomax, comm=MPI_COMM_WORLD)
!
    endsubroutine find_procno
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
! 
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
! 
         call keep_compiler_quiet(f)
! 
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f

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
!********************************************************************!
endmodule InitialCondition
