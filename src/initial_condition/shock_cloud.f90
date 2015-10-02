! $Id$
!
!  Shock-Cloud Initial Condition
!    sets up a 2D version similar to
!    The Magnetohydrodynamics of Shock-Cloud Interaction in Three Dimensions
!    Min-Su Shin, James M Stone, and Gregory F Snyder
!    http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2008ApJ...680..336S&link_type=ABSTRACT
! 
!  19-aug-13/mcnallcp: Initial commit, created for an out-of-date version of thermal_energy.f90
!    which is missing the call to the _ss init routine - need to remove the extra call after the school
!
! -------------------------------------
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
!
  integer :: dummy
!
  real :: xpos, uu_left, uu_right, widthuu
  real :: eth_left, eth_right
  real :: rho_left, rho_right
  real :: nval = 8.0, chib = 10.0, rcoreb = 0.62, rboundb = 1.77d0, machs = 10.0 
  real :: xblobi, yblobi, zblobi

  namelist /initial_condition_pars/ nval, chib, rcoreb, rboundb, machs, xblobi, yblobi, zblobi
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
      use EquationOfState, only: getmu, gamma, gamma_m1
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: pleft, pright, alfa, prat, cleft
      real :: tempr
!
      xpos = -2.66d0 ! shock position

!      nval = 8.0
!      chi = 10.0
!      rcore = 0.62
!      rbound = 1.77d0
!      machs = 10.0
!
      print*,'rho_max', 1.0 + (chib - 1.0)/( 1.0 + (0.0/rcoreb)**nval)
      pleft = 1.0
      !rholeft=one+(eqpar(chi_)-one)/(one+(eqpar(rbound_)/eqpar(rcore_))**eqpar(nval_))
      rho_left = 1.0 + (chib - 1.0)/( 1.0 + (rboundb/rcoreb)**nval)
      uu_left = 0.0
!
      ! compute the RH related states
      !Prat=one/(one+(eqpar(machs_)**2-one)*two*eqpar(gamma_)/(eqpar(gamma_)+one))
      prat = 1.0/(1.0 + (machs**2 -1.0)*2.0*gamma/(gamma+1.0))
      !alfa=(eqpar(gamma_)+1)/(eqpar(gamma_)-one)
      alfa = (gamma +1.0)/(gamma_m1)
      !cleft=dsqrt(eqpar(gamma_)*pleft/rholeft)
      cleft = sqrt(gamma * pleft/rho_left)
      !rhoright=rholeft*(alfa+Prat)/(alfa*Prat+one)
      rho_right = rho_left *(alfa + prat)/(alfa*prat + 1.0)
      !pright=pleft/Prat
      pright = pleft/prat
      !vright=cleft*eqpar(machs_)*(one-(alfa*Prat+one)/(alfa+Prat))
      uu_right = cleft * machs *(1.0-(alfa*prat+1.0)/(alfa+prat))
!
      widthuu = 3e-2
      eth_left = pleft/gamma_m1
      eth_right = pright/gamma_m1 
!
      if (lroot) then
        print*,' init Pressure ratio prat',prat,'alfa',alfa,'cleft',cleft
        print*,' init rho_left',rho_left,' rho_right', rho_right
        print*,'Values to use for B_ext'
!       p%beta1=0.5*p%b2*mu01/p%pp  is inverse beta
        print*,' assuming permittivity mu0 = 1.0'
        print*,'Beta = 0.5 B_ext = ', sqrt(1.0/0.5*pleft/1.0*2.0)
        print*,'Beta = 1   B_ext = ', sqrt(pleft/1.0*2.0)
        print*,'Beta = 10  B_ext = ', sqrt(1.0/10.0*pleft/1.0*2.0)
      endif
!     Rony's defs are backwards - swap
      tempr = rho_left
      rho_left = rho_right
      rho_right = tempr
      tempr = uu_left
      uu_left = uu_right
      uu_right = tempr
      tempr = eth_left
      eth_left = eth_right
      eth_right = tempr
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_all(f,profiles)
!
!  Initializes all the f arrays in one call.  This subroutine is called last.
!
!  21-dec-10/ccyang: coded
!  10-feb-15/MR    : added optional parameter 'profiles' (intended to replace f)
!
      real, dimension (mx,my,mz,mfarray), optional, intent(inout):: f
      real, dimension (:,:),              optional, intent(out)  :: profiles
!
!  SAMPLE IMPLEMENTATION
!
      !workaround for missing call in thermal_energy
      call initial_condition_ss(f)
      call keep_compiler_quiet(f)
!
     if (present(profiles)) then
       call fatal_error('initial_condition_all', &
                        'returning of profiles not implemented')
       call keep_compiler_quiet(profiles)
     endif
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
      call pjump(f,iux,xpos,uu_left,uu_right,widthuu,'x')
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
     
      call pjump(f,irho,xpos,rho_left,rho_right,widthuu,'x')
      !call sharpblob(ampli, f , irho, radiusi, xblobi, yblobi, zblobi)
      call ronyblob(chib, nval, f , irho, rcoreb, rboundb, xblobi, yblobi, zblobi)
      f(l1:l2,m1:m2,n1:n2,irho) = log(f(l1:l2,m1:m2,n1:n2,irho))
!
    endsubroutine initial_condition_lnrho

!***********************************************************************
    subroutine sharpblob(ampl,f,i,radius,xblob,yblob,zblob)
!
!  Single sharpblob.
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real, intent(in) :: xblob,yblob,zblob
      real, intent(in) :: ampl,radius
      real :: radius21
!
!  Single  blob. - specified in logrho
!
      radius21=1./radius**2
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,i)=f(l1:l2,m,n,i)+ampl*exp(-(((x(l1:l2)-xblob)**2+(y(m)-yblob)**2+(z(n)-zblob)**2)*radius21)**32.0)
      enddo; enddo
!
    endsubroutine sharpblob
!***********************************************************************
    subroutine ronyblob(chib,nval,f,i,rcoreb,rboundb,xblob,yblob,zblob)
!
!  Single RonyBlob.
!
      integer, intent(in) :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real, intent(in) :: xblob,yblob,zblob
      real, intent(in) :: chib,nval,rcoreb,rboundb
      integer :: l
      real :: radiuspoint, rhopoint
!
!  Single  blob. - specified in logrho
!
!the AMRVAC code is:
!rval(ix^S)=dsqrt(^D&x(ix^S,^D)**2+)
!where(rval(ix^S)<eqpar(rbound_))
!  w(ix^S,rho_)=one+(eqpar(chi_)-one)/(one+(rval(ix^S)/eqpar(rcore_))**eqpar(nval_))
      do n=n1,n2; do m=m1,m2; do l=l1,l2
        radiuspoint = sqrt((x(l)-xblob)**2+(y(m)-yblob)**2)
        if (radiuspoint < rbound) then
          rhopoint =  1.0 + (chib - 1.0) /(1.0 + (radiuspoint/rcoreb)**nval )
          f(l,m,n,i) = rhopoint
        endif
      enddo; enddo; enddo
!
    endsubroutine ronyblob
!***********************************************************************
    subroutine pjump(f,i,xpos,fleft,fright,width,dir)
!
!  jump
!
!  19-sep-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: profxy
      real, dimension (mx) :: profx
      real, dimension (my) :: profy
      real, dimension (mz) :: profz
      real :: fleft,fright,width,xpos
      character(len=*) :: dir
      integer :: l,m
!
!  jump; check direction
!
      select case (dir)
!
      case ('x')
        profx=fleft+(fright-fleft)*.5*(1.+tanh((x-xpos)/width))
        f(:,:,:,i)=f(:,:,:,i)+spread(spread(profx,2,my),3,mz)
      case default
        print*, 'pjump: not implemented'
      end select
    endsubroutine pjump
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call pjump(f,ieth,xpos,eth_left,eth_right,widthuu,'x')
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
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine initial_condition_aatest(f)
!
!  Initialize testfield.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_aatest
!***********************************************************************
    subroutine initial_condition_uutest(f)
!
!  Initialize testflow.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uutest
!***********************************************************************
    subroutine initial_condition_lncc(f)
!
!  Initialize passive scalar.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lncc
!***********************************************************************
    subroutine initial_condition_chiral(f)
!
!  Initialize chiral.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_chiral
!***********************************************************************
    subroutine initial_condition_chemistry(f)
!
!  Initialize chemistry.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_chemistry
!***********************************************************************
    subroutine initial_condition_uud(f)
!
!  Initialize dust fluid velocity.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uud
!***********************************************************************
    subroutine initial_condition_nd(f)
!
!  Initialize dust fluid density.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_nd
!***********************************************************************
    subroutine initial_condition_uun(f)
!
!  Initialize neutral fluid velocity.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uun
!***********************************************************************
    subroutine initial_condition_lnrhon(f)
!
!  Initialize neutral fluid density.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrhon
!***********************************************************************
    subroutine initial_condition_ecr(f)
!
!  Initialize cosmic rays.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ecr
!***********************************************************************
    subroutine initial_condition_fcr(f)
!
!  Initialize cosmic ray flux.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_fcr
!***********************************************************************
    subroutine initial_condition_solid_cells(f)
!
!  Initialize solid cells.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_solid_cells
!***********************************************************************
    subroutine initial_condition_cctest(f)
!
!  Initialize testscalar.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_cctest
!***********************************************************************
    subroutine initial_condition_xxp(f,fp)
!
!  Initialize particles' positions.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (:,:), intent(inout) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine initial_condition_xxp
!***********************************************************************
    subroutine initial_condition_vvp(f,fp)
!
!  Initialize particles' velocities.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (:,:), intent(inout) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine initial_condition_vvp
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
    subroutine initial_condition_clean_up
!
!  04-may-11/dhruba: coded
! dummy
!      
    endsubroutine initial_condition_clean_up
!***********************************************************************
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
