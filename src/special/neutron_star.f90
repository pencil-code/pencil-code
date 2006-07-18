! $Id: neutron_star.f90,v 1.4 2006-07-18 11:18:36 mee Exp $
!
!  This module incorporates all the modules used for Natalia's
!  neutron star -- disk coupling simulations (referred to as nstar)
!
!  This sample modules solves a special set of problems related
!  to computing the accretion through a thin disk onto a rigid surface
!  and the spreading of accreted material along the neutron star's surface.
!  One-dimensional problems along the disc and perpendicular to it
!  can also be considered.
!

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

!-------------------------------------------------------------------
!
! HOW TO USE THIS FILE
! --------------------
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the 
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users ot optionally check their contributions
! in to the Pencil-Code CVS repository.  This may be useful if you
! are working on/using the additional physics with somebodyelse or 
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module 
! selections to say something like:
!
!    SPECIAL=special/nstar
!
! Where nstar it replaced by the filename of your new module
! upto and not including the .f90
!
!--------------------------------------------------------------------

module Special

  use Cparam
  use Cdata
  use Messages
!  use Density, only: rho_up
  use EquationOfState

  implicit none

  include 'special.h'
  
  ! input parameters 
 ! logical :: sharp=.false., smooth=.false.

  logical :: lmass_source_NS=.false. 
  logical :: leffective_gravity=.false.

  character (len=labellen), dimension(ninit) :: initnstar='default'
  real :: rho_star=1.,rho_disk=1., rho_surf=1.
!
! Keep some over used pencils
!
  real, dimension(nx) :: z_2

  namelist /neutron_star_init_pars/ &
      initnstar,lmass_source_NS,leffective_gravity, rho_star,rho_disk,rho_surf!,sharp

  ! run parameters

  namelist /neutron_star_run_pars/ &
      lmass_source_NS,leffective_gravity, rho_star,rho_disk,rho_surf

!!
!! Declare any index variables necessary for main or 
!! 
!!   integer :: iSPECIAL_VARIABLE_INDEX=0
!!  
!! other variables (needs to be consistent with reset list below)
!!
!!   integer :: i_POSSIBLEDIAGNOSTIC=0
!!

  contains

!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables 
!  which should be know to be able to evaluate
! 
!
!  6-oct-03/tony: coded
!
      use Cdata
   !   use Density
      use EquationOfState
      use Mpicomm
!
      logical, save :: first=.true.
!
! A quick sanity check
!
      if (.not. first) call stop_it('register_special called twice')
      first = .false.

!!
!! MUST SET lspecial = .true. to enable use of special hooks in the Pencil-Code 
!!   THIS IS NOW DONE IN THE HEADER ABOVE
!
!
!
!! 
!! Set any required f-array indexes to the next available slot 
!!  
!!
!      iSPECIAL_VARIABLE_INDEX = nvar+1             ! index to access entropy
!      nvar = nvar+1
!
!      iSPECIAL_AUXILLIARY_VARIABLE_INDEX = naux+1             ! index to access entropy
!      naux = naux+1
!
!
!  identify CVS version information (if checked in to a CVS repository!)
!  CVS should automatically update everything between $Id: neutron_star.f90,v 1.4 2006-07-18 11:18:36 mee Exp $ 
!  when the file in committed to a CVS repository.
!
      if (lroot) call cvs_id( &
           "$Id: neutron_star.f90,v 1.4 2006-07-18 11:18:36 mee Exp $")
!
!
!  Perform some sanity checks (may be meaningless if certain things haven't 
!  been configured in a custom module but they do no harm)
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('register_special: naux > maux')
      endif
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_special: nvar > mvar')
      endif
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      use Cdata 
   !   use Density
      use EquationOfState

!
      real, dimension (mx,my,mz,mvar+maux) :: f
!!
!!  Initialize any module variables which are parameter dependant  
!!
!
! DO NOTHING
      if(NO_WARN) print*,f  !(keep compiler quiet)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f,xx,yy,zz)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Cdata
   !   use Density
      use EquationOfState
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      intent(in) :: xx,yy,zz
      intent(inout) :: f

!!
      select case(initnstar)
        case('default')
          if(lroot) print*,'init_special: Default neutron star setup'
          call density_step(f,xx,zz)
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_special: No such value for initnstar: ', trim(initnstar)
          call stop_it("")
      endselect
!
      if(NO_WARN) print*,f,xx,yy,zz  !(keep compiler quiet)
!
    endsubroutine init_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   24-nov-04/tony: coded
!
      use Cdata  
  !    use Density
      use EquationOfState
      
!
      real, dimension (mx,my,mz,mvar+maux) :: f       
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!     
      if(NO_WARN) print*,f(1,1,1,1),p   !(keep compiler quiet)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  several precalculated Pencils of information are passed if for
!  efficiency.
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Global
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df     
      type (pencil_case) :: p

!
      intent(in) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
!!      if (headtt) call identify_bcs('ss',iss)
!
!!
!! SAMPLE DIAGNOSTIC IMPLEMENTATION
!!
!!      if(ldiagnos) then
!!        if(i_SPECIAL_DIAGNOSTIC/=0) then
!!          call sum_mn_name(SOME MATHEMATICAL EXPRESSION,i_SPECIAL_DIAGNOSTIC)
!!! see also integrate_mn_name
!!        endif
!!      endif

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,f,df,p

    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=special_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_init_pars,ERR=99)
      endif

99    return
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=special_init_pars)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
    
      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif

99    return
endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
                                                                                                   
      if (NO_WARN) print*,unit
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!   06-oct-03/tony: coded
!
!AB:  use Cdata   !(not needed, right?)
      use Sub

      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite

!!
!!!   SAMPLE IMPLEMENTATION
!!
!!      integer :: iname
!!!
!!!  reset everything in case of reset
!!!  (this needs to be consistent with what is defined above!)
!!!
      if (lreset) then
!!        i_SPECIAL_DIAGNOSTIC=0
      endif
!!
!!      do iname=1,nname
!!        call parse_name(iname,cname(iname),cform(iname),'NAMEOFSPECIALDIAGNOSTIC',i_SPECIAL_DIAGNOSTIC)
!!      enddo
!!
!!!  write column where which magnetic variable is stored
!!      if (lwr) then
!!        write(3,*) 'i_SPECIAL_DIAGNOSTIC=',i_SPECIAL_DIAGNOSTIC
!!      endif
!!

    endsubroutine rprint_special
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
     ! use Density    
     ! use EquationOfState
    
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      integer :: i, l_sz, tmp_int
!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!  
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
!!


!  mass sources and sinks for the boundary layer on NS in 1D approximation
!
      if (lmass_source_NS) call mass_source_NS(f,df,p%rho)
!
! Natalia
! deceleration zone in a case of a Keplerian disk

       if (laccelerat_zone) then
         if (n .gt. nzgrid-ac_dc_size .and. dt .gt. 0.) then
            if (lnstar_entropy) then   
              if (nxgrid /= 1) then       
!natalia  .... REASON FOR BLOCK COMMENT ...
!natalia                df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)&
!natalia                  -1./(5.*dt)*(f(l1:l2,m,n,ilnrho)-f(l1:l2,m,n-1,ilnrho))
!natalia
!natalia
!natalia                df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)&
!natalia                  -1./(5.*dt)*(f(l1:l2,m,n,ilnrho)-f(l1:l2,m,nzgrid-ac_dc_size,ilnrho))
!natalia        
!natalia                df(l1:H_disk_point+4,m,n,ilnrho)=df(l1:H_disk_point+4,m,n,ilnrho) &
!natalia                  -1./(5.*dt)*(f(l1:H_disk_point+4,m,n,ilnrho) &
!natalia                 -log(rho_surf)-(1.-(x(1:H_disk_point)/H_disk)**2))
!natalia
!natalia                df(l1:H_disk_point+4,m,n,ilnrho)=df(l1:H_disk_point+4,m,n,ilnrho) &
!natalia                  -1./(5.*dt)*(f(l1:H_disk_point+4,m,n,ilnrho) &
!natalia                  -log(rho_surf)-(1.-(x(1:H_disk_point)/H_disk)))
!natalia     
!natalia                do i=1,H_disk_point+4
!natalia                  df(i,m,n,ilnrho)=df(i,m,n,ilnrho)&
!natalia                    -1./(5.*dt)*(f(i,m,n,ilnrho)-f(i+1,m,n,ilnrho))
!natalia                enddo    
!natalia                df(H_disk_point+5:l2,m,n,ilnrho)=df(H_disk_point+5:l2,m,n,ilnrho) &
!natalia                  -1./(5.*dt)*(f(H_disk_point+5:l2,m,n,ilnrho)
!natalia                  -log(rho_surf)-(1.-(x(H_disk_point+1)/H_disk)))
!natalia       
!natalia                df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) & 
!natalia                  -1./(5.*dt)*(f(l1:l2,m,n,ilnrho) &
!natalia                  -log(rho_surf)-(1.-M_star/2./z(n)**3*x(l1:l2)**2*mu/Rgas/T_star))
!natalia
!natalia                df(H_disk_point+5:l2,m,n,ilnrho)=df(H_disk_point+5:l2,m,n,ilnrho) &
!natalia                  -1./(5.*dt)*(f(H_disk_point+5:l2,m,n,ilnrho)-f(H_disk_point+5:l2,m,n-1,ilnrho))
!natalia
                if (H_disk_point .ge. nxgrid) then
                  df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) &
                    -1./(5.*dt)*(f(l1:l2,m,n,ilnrho) &
                    -log(rho_surf)-(1.-M_star/2./z(n)**3*x(l1:l2)**2*mu/Rgas/T_star))
                else
                  df(l1:H_disk_point+4,m,n,ilnrho)= &
                    df(l1:H_disk_point+4,m,n,ilnrho) &
                    -1./(5.*dt)*(f(l1:H_disk_point+4,m,n,ilnrho) &
                    -log(rho_surf)-(1.-M_star/2./z(n)**3 &
                                    * x(l1:H_disk_point+4)**2*mu/Rgas/T_star))
                  
                  df(H_disk_point+5:l2,m,n,ilnrho)=i &
                    df(H_disk_point+5:l2,m,n,ilnrho) &
                    -1./(5.*dt)*(f(H_disk_point+5:l2,m,n,ilnrho) &
                    -log(rho_surf)-(1.-M_star/2./z(n)**3 &
                                    * x(H_disk_point+4)**2*mu/Rgas/T_star))
                endif ! (H_disk_point .lt. nxgrid)
              else
!  if nxgrid==1
                if (lnstar_T_const) then
                else            
                endif    
              endif 
            endif ! (lnstar_entropy)
         endif  ! (
!
       if (ldecelerat_zone) then
!
         if (n .le. ac_dc_size+4 .and. dt .gt. 0.) then
!
            if (lnstar_entropy) then          
              if (lnstar_T_const) then
                df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) &
                  -0.2*p%rho1(:)/dt &
                  *(p%rho(:)-rho_star*exp(-M_star/R_star/cs0**2*gamma*(1.-R_star/z(n))))
              else            
!
!natalia   .... REASON FOR BLOCK COMMENT ...
!natalia                 df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) &
!natalia                   -1./p%rho(:)/(5.*dt) &
!natalia                   *(p%rho(:)-rho_star*exp(-M_star/R_star/p%cs2(:)*(1.-R_star/z(n))))
!natalia                 df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) &
!natalia                   -1./p%rho(:)/(5.*dt) &
!natalia                   *(p%rho(:)-rho_star*exp(-M_star/R_star/(gamma1*T_star)*gamma*(1.-R_star/z(n))))
!natalia                
!natalia                 df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) &
!natalia                   -1./(5.*dt)*(p%rho(:)-rho_star)/p%rho(:)
              endif    
            else
              df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) &
                 -1./p%rho(:)/(5.*dt) *(p%rho(:) &
                 -rho_star*exp(-M_star/R_star/p%cs2(:)*(1.-R_star/z(n))))
            endif
          endif
        endif
      endif  
!
! surface zone in a case of a Keplerian disk
!
      if (lsurface_zone) then
          if ( dt .gt.0.) then
            l_sz=l2-5
!
!            df(l_sz:l2,m,n,ilnrho)=df(l_sz:l2,m,n,ilnrho)&
!                 -1./(5.*dt)*(1.-rho_surf/exp(f(l_sz:l2,m,n,ilnrho)))
!
           do i=l_sz,l2   
             df(i,m,n,ilnrho)=df(i,m,n,ilnrho)&
                   -0.2/dt*(f(i,m,n,ilnrho)-f(i-1,m,n,ilnrho))
           enddo
         endif ! ()
      endif
!
! Keep compiler quiet by ensuring every parameter is used
!
      if (NO_WARN) print*,df,p
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   16-jul-06/natalia: coded
!
      use Cdata
!      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      integer :: j,l_sz
!
! add effective gravity term = -Fgrav+Fcentrifugal
! Natalia
!
      if (leffective_gravity) then
        if (headtt) &
          print*,'duu_dt: Effectiv gravity; Omega, Rstar=', Omega, R_star, M_star

          df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)- &
            M_star/z(n)**2*(1.-p%uu(:,2)*p%uu(:,2)*z(n)/M_star)
          if (nxgrid /= 1) then  
            df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)- &
              M_star/z(n)**2/sqrt(z(n)**2+x(l1:l2)**2)*x(l1:l2)*(z(n)-R_star)/(Lxyz(1)*0.5)
          endif
      endif
!
! acceleration zone in a case of a Keplerian disk
!
      if (laccelerat_zone) then
        if (n .ge. nzgrid-ac_dc_size  .and. dt .gt.0.) then
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)&
             -1./(5.*dt)*(p%uu(:,2)-sqrt(M_star/z(n)))
       
          if (nxgrid == 1) then  
            df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)&
              -1./(5.*dt)*(p%uu(:,3)+accretion_flux/p%rho(:))
          else
            df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-1./(5.*dt)*(p%uu(:,1)-0.)
            df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)&
            -1./(5.*dt)*(f(l1:l2,m,n,iuz)-f(l1:l2,m,n-1,iuz))

!natalia  .... REASON FOR BLOCK QUOTE ...
!natalia            df(l1:H_disk_point+4,m,n,iuz)=df(l1:H_disk_point+4,m,n,iuz)&
!natalia              -1./(5.*dt)*(f(l1:H_disk_point+4,m,n,iuz)-f(l1:H_disk_point+4,m,n-1,iuz))
!natalia
!natalia            df(l1:H_disk_point+4,m,n,iuz)=df(l1:H_disk_point+4,m,n,iuz)&
!natalia              -1./(5.*dt)*(p%uu(1:H_disk_point,3)+accretion_flux/p%rho(1:H_disk_point))
!natalia
!natalia            df(H_disk_point+5:l2,m,n,iuz)=df(H_disk_point+5:l2,m,n,iuz)&
!natalia              -1./(5.*dt)*(f(H_disk_point+5:l2,m,n,iuz)-f(H_disk_point+5:l2,m,nzgrid-ac_dc_size,iuz))
           endif 
         endif
       endif
!
! deceleration zone in a case of a Keplerian disk
!
      if (ldecelerat_zone) then
        if (n .le. ac_dc_size+4  .and. dt .gt.0.) then
          if (lnstar_entropy) then  
!            if (lnstar_T_const) then
            df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)&
!               -1./(5.*dt)*(p%uu(:,2)-f(l1:l2,m,n-1,iuy))
               -1./(5.*dt)*(p%uu(:,2)-0.)   
                 
            df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)&
!               -1./(5.*dt)*(p%uu(:,3)-f(l1:l2,m,n-1,iuz)) 
               -1./(5.*dt)*(p%uu(:,3)-0.)   
!            endif  
          else 
            df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)&
                           -1./(5.*dt)*(p%uu(:,2)-0.)
!
            df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)&
                          -1./(5.*dt)*(p%uu(:,3)-0.)
!
          endif   
        endif
      endif

! surface zone in a case of a Keplerian disk

      if (lsurface_zone) then
        if ( dt .gt.0.) then
!
          l_sz=l2-5
!

!natalia  ... REASON FOR BLOCK QUOTE ...
!          do j=l_sz,l2   
!            df(j,m,n,iux)=df(j,m,n,iux)&
!                 -1./(3.*dt)*(-f(j-1,m,n,iux)+f(j,m,n,iux))
!            df(j,m,n,iux)=df(j,m,n,iux)&
!                 -1./(10.*dt)*(f(j,m,n,iux)-f(j+1,m,n,iux))
!          enddo
!
        if (lnstar_1D) then
             df(l_sz:l2,m,n,iux)=df(l_sz:l2,m,n,iux)&
                   -1./(2.*dt)*(f(l_sz:l2,m,n,iux)-0.)
!
!
            df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)&
                  -1./(5.*dt)*(f(l1:l2,m,n,iuy)-sqrt(M_star/xyz0(3)))
!
            df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)&
             -1./(5.*dt)*(f(l1:l2,m,n,iuz)+accretion_flux/p%rho(:))
!
         else
           do j=l_sz,l2   
             df(j,m,n,iux)=df(j,m,n,iux)&
                  -1./(5.*dt)*(f(j,m,n,iux)-f(j-1,m,n,iux))
             df(j,m,n,iuy)=df(j,m,n,iuy)&
                  -1./(5.*dt)*(f(j,m,n,iuy)-f(j-1,m,n,iuy))
             df(j,m,n,iuz)=df(j,m,n,iuz)&
                  -1./(5.*dt)*(f(j,m,n,iuz)-f(j-1,m,n,iuz))
           enddo
!
!natalia ... REASON FOR BLOCK QUOTE ...
!           df(l_sz:l2,m,n,iux)=df(l_sz:l2,m,n,iux)&
!                  -1./(2.*dt)*(f(l_sz:l2,m,n,iux)-0.)
!
!
!           df(l_sz:l2,m,n,iuy)=df(l_sz:l2,m,n,iuy)&
!                  -1./(5.*dt)*(f(l_sz:l2,m,n,iuy)-sqrt(M_star/z(n)))
!
!           df(l_sz:l2,m,n,iuz)=df(l_sz:l2,m,n,iuz)&
!                  -1./(5.*dt)*(f(l_sz:l2,m,n,iuz)+accretion_flux/p%rho(:))
 
     
         endif 
    
        
         endif

      endif

    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!  
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
!!

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,df,p

    endsubroutine special_calc_magnetic
!!***********************************************************************
    subroutine special_calc_entropy(df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!  
!!
!!  df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + SOME NEW TERM
!!
!!

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,df,p

    endsubroutine special_calc_entropy
!*************************************************************************
    subroutine mass_source_NS(f,df,rho)
!
!  add mass sources and sinks
!
!  2006/Natalia
!
      use Cdata
!     
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
!     real, dimension(nx) :: fint,fext,pdamp
      real  ::  sink_area, V_acc, V_0, rho_0, ksi, integral_rho=0., flux
      integer :: sink_area_points=50, i
      integer :: idxz  
      real, dimension (nx), intent(in) :: rho 

       sink_area=Lxyz(3)/(nzgrid-1.)*sink_area_points

! 
!  No clue what this index is good for, but nzgrid-30 is not a
!  valid index for e.g. 2-d runs, so sanitize it to avoid
!  `Array reference at (1) is out of bounds' with g95 -Wall
! 
       idxz = min(nzgrid-30,n2)

       V_0=f(4,4,idxz,iuz)
       rho_0=exp(f(4,4,idxz,ilnrho))
       flux=accretion_flux
   
       flux=V_0*rho_0     

!       V_0=rho_0*V_acc*(sink_area_points+1)/integral_rho
!
!
!       ksi=2.*((Lxyz(3)/(nzgrid-1.)*(sink_area_points+4-n))/sink_area)/sink_area
       ksi=1./sink_area

       if ( 25 .gt. n .and. n .lt. sink_area_points+25) then 
         df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)-flux*ksi/rho(:)
       endif
 
       if ( n .eq. 25 .or. n .eq. sink_area_points+25) then
         df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)-0.5*flux*ksi/rho(:)
       endif
!
       if (headtt) print*,'dlnrho_dt: mass source*rho = ', flux/sink_area
!
    endsubroutine mass_source_NS
!***********************************************************************
    subroutine density_step(f,xx,zz)
!
!Natalia
!Initialization of density in a case of the step-like distribution
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (my,mz) :: lnrho_2d
      real, dimension (mx,my,mz) :: xx, zz
!
      integer :: step_width, step_length,i
      real :: H_disk_min, L_disk_min, hdisk, ldisk, ll,  ln_ro_l, ln_ro_r, ln_ro_u
!
      hdisk=H_disk 
      ldisk=L_disk
!
      H_disk_min=Lxyz(1)/(nxgrid-1)
      L_disk_min=Lxyz(3)/(nzgrid-1)
!
      if (H_disk .gt. Lxyz(1)-H_disk_min) hdisk=Lxyz(1)
      if (H_disk .lt. H_disk_min) hdisk=0.
!
      if (L_disk .gt. Lxyz(3)-L_disk_min) ldisk=Lxyz(3)
      if (L_disk .lt. L_disk_min) ldisk=0.
!
      step_width=nint((nxgrid-1)*hdisk/Lxyz(1))
      step_length=nint((nzgrid-1)*(Lxyz(3)-ldisk)/Lxyz(3))
!
      if (hdisk .EQ. Lxyz(1) .AND. ldisk .EQ. Lxyz(3))  f(:,:,:,ilnrho)=log(rho_star)
!
      if (hdisk .EQ. 0. .AND. ldisk .EQ. 0.) f(:,:,:,ilnrho)=log(rho_disk)
!
      if (hdisk .EQ. Lxyz(1) .AND. ldisk .LT. Lxyz(3)) then
        f(:,:,1:step_length+3,ilnrho)=log(rho_disk)
        f(:,:,step_length+3+1:mz,ilnrho)=log(rho_star)
      endif
!
      if (sharp) then
        if (hdisk .LT. Lxyz(1) .AND. ldisk .EQ. Lxyz(3)) then
          f(1:step_width+3,:,:,ilnrho)=log(rho_star)
          f(step_width+3+1:mx,:,:,ilnrho)=log(rho_disk)
        endif
!
        if (hdisk .GT. 0.  .AND. hdisk .LT. Lxyz(1) ) then
          if (ldisk .GT. 0.  .AND. ldisk .LT. Lxyz(3)) then
            f(1:step_width+3,:,step_length+3+1:mz,ilnrho)=log(rho_star)
            f(step_width+3+1:mx,:,step_length+3+1:mz,ilnrho)=log(rho_disk)
            f(:,:,1:step_length+3,ilnrho)=log(rho_disk)
          end if
        end if
      end if 
!
      if (smooth) then
        ln_ro_r=log(rho_disk)
        ln_ro_l=log(rho_star)
        ln_ro_u=log(rho_surf)
!
        ll=Lxyz(3)-ldisk
!
        if (nxgrid/=1.and.nzgrid/=1) then
!natalia  ... REASON FOR BLOCK COMMENT ...
!natalia          lnrho_2d(:,:)=(zz(l1,:,:)-R_star)/Lxyz(3)*(ln_ro_r-ln_ro_l)+ln_ro_l
!natalia
!natalia          do i=1,l1
!natalia           f(i,:,:,ilnrho)=lnrho_2d(:,:)
!natalia          enddo 
!natalia       
!natalia          do i=l1+1,mx
!natalia           f(i,:,:,ilnrho)=(xx(i,:,:)-xx(l1,:,:))/Lxyz(1)*(ln_ro_l-lnrho_2d(:,:))+ln_ro_l 
!natalia          enddo        
!natalia
!natalia          do i=1,H_disk_point+4
!natalia           f(i,:,:,ilnrho)=ln_ro_r
!natalia          enddo 
!natalia       
!natalia          do i=H_disk_point+5,mx
!natalia            f(i,:,:,ilnrho)=ln_ro_u 
!natalia          enddo    
!natalia
!natalia
!natalia        if (H_disk.GT.0.) f(:,:,:,ilnrho)=f(:,:,:,ilnrho)+(1.-(xx(:,:,:)/H_disk)**2)
!natalia          if (H_disk.GT.0.) 
!natalia        f(:,:,:,ilnrho)=ln_ro_u+(1.-(xx(:,:,:)/H_disk)**2)
!
          do i=1,H_disk_point_int+4
!
            f(i,:,:,ilnrho)=ln_ro_u+(1.-(xx(i,:,:)/H_disk)**2)
!natalia           f(i,:,:,ilnrho)=ln_ro_u+(1.-M_star/2./zz(i,:,:)**3*x(i)**2*mu/Rgas/T_star)
 
          enddo 

          do i=H_disk_point_int+5,mx
 
            f(i,:,:,ilnrho)=f(H_disk_point_int+4,:,:,ilnrho)

!natalia            f(i,:,:,ilnrho)=f(H_disk_point+4,:,:,ilnrho)
!natalia            f(i,:,:,ilnrho)=ln_ro_u+(1.-M_star/2./zz(i,:,:)**3*x(H_disk_point+4)**2*mu/Rgas/T_star)
          enddo    
        else 
          if (nzgrid .GT. 1) then 
            f(:,:,:,ilnrho)=(zz(:,:,:)-R_star)/Lxyz(3)*(ln_ro_r-ln_ro_l)+ln_ro_l
          else
            if (H_disk.GT.0.) f(:,:,:,ilnrho)=(xx(:,:,:)-0.)/Lxyz(1)*(ln_ro_u-ln_ro_r)+ln_ro_r
          endif 
        endif     
      endif

    endsubroutine density_step
!***********************************************************************
endmodule Special

