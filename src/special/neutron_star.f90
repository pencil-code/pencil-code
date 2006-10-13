! $Id: neutron_star.f90,v 1.16 2006-10-13 15:14:46 nbabkovs Exp $
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
  logical :: lsharp=.false., lsmooth=.false.

  logical :: lmass_source_NS=.false. 
  logical :: leffective_gravity=.false.

  character (len=labellen) :: initnstar='default'
  real :: rho_star=1.,rho_disk=1., rho_surf=1.

  real :: uu_left=0.
  real :: uy_left=0.,uy_right=0.
 
  real :: H_disk=0.
  real :: L_disk=0.
  real :: R_star=0.
  real :: M_star=0. 
  real :: T_star=0.
  real :: T_disk=0.
  real :: accretion_flux=0.

  logical :: lextrapolate_bot_density=.false.
  logical :: ltop_velocity_kep=.false.
  logical :: laccelerat_zone=.false.
  logical :: ldecelerat_zone=.false.
  logical :: lsurface_zone=.false.
  logical :: lnstar_T_const=.false.
  logical :: lnstar_entropy=.false.
  logical :: lnstar_1D=.false.
  integer :: ac_dc_size=5
  integer :: H_disk_point=0

  real :: beta_hand=1.
  real :: nu_for_1D=1. 
  real :: mu_local=0.
  logical :: l1D_cooling=.false.,l1D_heating=.false.
  logical :: l1D_cool_heat=.false.
  logical :: lgrav_x_mdf=.false. 
!
! Keep some over used pencils
!
  real, dimension(nx) :: z_2
  integer :: H_disk_point_int=0

! start parameters
  namelist /neutron_star_init_pars/ &
      initnstar,lmass_source_NS,leffective_gravity, &
      laccelerat_zone, ldecelerat_zone, lsurface_zone, &
      rho_star,rho_disk,rho_surf, &
      H_disk_point_int, &
      H_disk, H_disk_point, &
      L_disk, R_star, M_star, &
      T_star,accretion_flux, T_disk, &
      uu_left, uy_left, uy_right, &
      l1D_cooling,l1D_heating,beta_hand, &
      nu_for_1D, ltop_velocity_kep, lextrapolate_bot_density, &
      lnstar_entropy, lnstar_T_const, lnstar_1D, &
      lgrav_x_mdf

! run parameters
  namelist /neutron_star_run_pars/ &
      lmass_source_NS,leffective_gravity, rho_star,rho_disk,rho_surf, &
      laccelerat_zone, ldecelerat_zone, lsurface_zone, &
       H_disk, H_disk_point, &
       L_disk, R_star, M_star, T_star, &
       accretion_flux, lnstar_entropy, &
       lnstar_T_const,lnstar_1D, &
       l1D_cooling,l1D_heating, lgrav_x_mdf
!!
!! Declare any index variables necessary for main or 
!! 
!!   integer :: iSPECIAL_VARIABLE_INDEX=0
!  
!! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_dtcrad=0
  integer :: idiag_dtchi=0
!
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
!  CVS should automatically update everything between $Id: neutron_star.f90,v 1.16 2006-10-13 15:14:46 nbabkovs Exp $ 
!  when the file in committed to a CVS repository.
!
      if (lroot) call cvs_id( &
           "$Id: neutron_star.f90,v 1.16 2006-10-13 15:14:46 nbabkovs Exp $")
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
    l1D_cool_heat=l1D_cooling.or.l1D_heating
    if (l1D_cool_heat.and.lroot) &
          print*, 'neutron_star: 1D cooling or heating'
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
          call entropy_step(f,xx,zz,T_star)
          call velocity_step(f)
        case('sharp')
          if(lroot) print*,'init_special: Sharp neutron star setup'
          lsharp=.true.
          call density_step(f,xx,zz)
          call entropy_step(f,xx,zz,T_star)
          call velocity_kep_disk(f,zz)
        case('smooth')
          if(lroot) print*,'init_special: Sharp neutron star setup'
          lsmooth=.true.
          call density_step(f,xx,zz)
          call entropy_step(f,xx,zz,T_star)
          call velocity_kep_disk(f,zz)
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
    subroutine pencil_criteria_special()
! 
!  All pencils that this special module depends on are specified here.
! 
!  18-07-06/tony: coded
!
      use Cdata
!
      if (laccelerat_zone)  lpenc_requested(i_rho)=.true.
    !  if (lmass_source_NS)  lpenc_requested(i_rho)=.true.
!Natalia (accretion on a NS)
       if (lnstar_entropy) then
         lpenc_requested(i_TT)=.true.
          lpenc_requested(i_lnTT)=.true.
         lpenc_requested(i_cs2)=.true.
         lpenc_requested(i_ss)=.true.
         lpenc_requested(i_rho)=.true.
      endif
!
!
    endsubroutine pencil_criteria_special
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
      if(ldiagnos) then
        if (idiag_dtcrad/=0) &
          call max_mn_name(sqrt(advec_crad2)/cdt,idiag_dtcrad,l_dt=.true.)
        if (idiag_dtchi/=0) &
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
      endif

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,f,df,p

    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=neutron_star_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=neutron_star_init_pars,ERR=99)
      endif

99    return
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=neutron_star_init_pars)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
    
      if (present(iostat)) then
        read(unit,NML=neutron_star_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=neutron_star_run_pars,ERR=99)
      endif

99    return
endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
                                                                                                   
      write(unit,NML=neutron_star_run_pars)

    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!   06-oct-03/tony: coded
!
      use Sub
!
!  define diagnostics variable
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite

!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtcrad=0
        idiag_dtchi=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtcrad',idiag_dtcrad)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
      enddo
!
!  write column where which magnetic variable is stored
      if (lwr) then
        write(3,*) 'i_dtcrad=',idiag_dtcrad
        write(3,*) 'i_dtchi=',idiag_dtchi
      endif
!
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
      use EquationOfState
    
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      integer :: i, l_sz, tmp_int
      real :: cs2_star
!


    
       call eoscalc(ilnrho_lnTT,log(rho_star),log(T_star), cs2=cs2_star)

!  mass sources and sinks for the boundary layer on NS in 1D approximation
!
      if (lmass_source_NS) call mass_source_NS(f,df,p%rho)
!
! Natalia
! deceleration zone in a case of a Keplerian disk

 

       if (laccelerat_zone) then
     
         
         if (n .GE. nzgrid-ac_dc_size .AND. dt .GT. 0.) then
       
            if (lnstar_entropy) then   
             
             if (nxgrid .LE. 1) then       
              if (lnstar_T_const) then
              else            
              endif    

             else

         ! df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)&
         !     -1./(5.*dt)*(f(l1:l2,m,n,ilnrho)-f(l1:l2,m,n-1,ilnrho))


         !    df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)&
         !     -1./(5.*dt)*(f(l1:l2,m,n,ilnrho)-f(l1:l2,m,nzgrid-ac_dc_size,ilnrho))
        
        !df(l1:H_disk_point+4,m,n,ilnrho)=df(l1:H_disk_point+4,m,n,ilnrho) &
        !  -1./(5.*dt)*(f(l1:H_disk_point+4,m,n,ilnrho) &
        !  -log(rho_surf)-(1.-(x(1:H_disk_point)/H_disk)**2))

        !  df(l1:H_disk_point+4,m,n,ilnrho)=df(l1:H_disk_point+4,m,n,ilnrho) &
        !  -1./(5.*dt)*(f(l1:H_disk_point+4,m,n,ilnrho) &
        !  -log(rho_surf)-(1.-(x(1:H_disk_point)/H_disk)))
     
    ! do i=1,H_disk_point+4
    !         df(i,m,n,ilnrho)=df(i,m,n,ilnrho)&
    !               -1./(5.*dt)*(f(i,m,n,ilnrho)-f(i+1,m,n,ilnrho))
    ! enddo    
    !     df(H_disk_point+5:l2,m,n,ilnrho)=df(H_disk_point+5:l2,m,n,ilnrho) &
    !-1./(5.*dt)*(f(H_disk_point+5:l2,m,n,ilnrho)
    !-log(rho_surf)-(1.-(x(H_disk_point+1)/H_disk)))
       
        ! df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) & 
        !   -1./(5.*dt)*(f(l1:l2,m,n,ilnrho) &
	!   -log(rho_surf)-(1.-M_star/2./z(n)**3*x(l1:l2)**2*gamma/cs2_star))
			     



        ! df(H_disk_point+5:l2,m,n,ilnrho)=df(H_disk_point+5:l2,m,n,ilnrho) &
        ! -1./(5.*dt)*(f(H_disk_point+5:l2,m,n,ilnrho)-f(H_disk_point+5:l2,m,n-1,ilnrho))

           if (H_disk_point .GE. nxgrid) then
	  
	     df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) &
	     -1./(5.*dt)*(f(l1:l2,m,n,ilnrho) &
             -log(rho_surf)-(1.-M_star/2./z(n)**3*x(l1:l2)**2*gamma/cs2_star))
			   
			   
	   
	   else
	   
	 	   
             df(l1:H_disk_point+4,m,n,ilnrho)=df(l1:H_disk_point+4,m,n,ilnrho) &
             -1./(5.*dt)*(f(l1:H_disk_point+4,m,n,ilnrho) &
     	     -log(rho_surf)-(1.-M_star/2./z(n)**3*x(l1:H_disk_point+4)**2*gamma/cs2_star))
			
	  
			      
             df(H_disk_point+5:l2,m,n,ilnrho)=df(H_disk_point+5:l2,m,n,ilnrho) &
             -1./(5.*dt)*(f(H_disk_point+5:l2,m,n,ilnrho) &
	     -log(rho_surf)-(1.-M_star/2./z(n)**3*x(H_disk_point+4)**2*gamma/cs2_star))
	    endif
						    
             endif 
           endif
           


         endif 

       if (ldecelerat_zone) then
        
         if (n .LE. ac_dc_size+4 .AND. dt .GT. 0.) then
       
            if (lnstar_entropy) then          
              if (lnstar_T_const) then
               df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) &
                -1./p%rho(:)/(5.*dt) &
                *(p%rho(:)-rho_star*exp(-M_star/R_star/cs0**2*gamma*(1.-R_star/z(n))))
               else            

              !   df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) &
              !           -1./p%rho(:)/(5.*dt) &
              !   *(p%rho(:)-rho_star*exp(-M_star/R_star/p%cs2(:)*(1.-R_star/z(n))))
              ! df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) &
              !  -1./p%rho(:)/(5.*dt) &
              !  *(p%rho(:)-rho_star*exp(-M_star/R_star/(gamma1*T_star)*gamma*(1.-R_star/z(n))))
              
           !  df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) &
           !   -1./(5.*dt)*(p%rho(:)-rho_star)/p%rho(:)
          
            endif    

            else
               df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) &
                         -1./p%rho(:)/(5.*dt) &
                 *(p%rho(:)-rho_star*exp(-M_star/R_star/p%cs2(:)*(1.-R_star/z(n))))
            endif
      
          endif
       endif
           
     endif  

! surface zone in a case of a Keplerian disk

      if (lsurface_zone) then


          if ( dt .GT.0.) then
            l_sz=l2-5

          !  df(l_sz:l2,m,n,ilnrho)=df(l_sz:l2,m,n,ilnrho)&
          !       -1./(5.*dt)*(1.-rho_surf/exp(f(l_sz:l2,m,n,ilnrho)))

           if (lnstar_1D) then
 
           else 
            do i=l_sz,l2   
              df(i,m,n,ilnrho)=df(i,m,n,ilnrho)&
              -1./(5.*dt)*(f(i,m,n,ilnrho)-f(i-1,m,n,ilnrho) &
	      +M_star/z(n)**3*(x(i)-x(i-1))*x(i-1)*gamma/cs2_star)
          		    
		!  df(i,m,n,ilnrho)=df(i,m,n,ilnrho)&
                !   -1./(5.*dt)*(f(i,m,n,ilnrho)-f(i-1,m,n,ilnrho))		    
		   
           enddo
           endif 
		


         endif
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

          if (lgrav_x_mdf) then 
            df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)- &
              M_star/z(n)**2/sqrt(z(n)**2+x(l1:l2)**2)*x(l1:l2)*(z(n)-R_star)/(Lxyz(1)*0.5)
           else
            df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-M_star/z(n)**3*x(l1:l2)
           endif
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
            df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)&
                           -1./(5.*dt)*(p%uu(:,1)-0.)

            
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
    subroutine special_calc_magnetic(f,df,p)
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
      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
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
    subroutine special_calc_entropy(f,df,p)
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
      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      integer :: j, l_sz, l_sz_1

      if (l1D_cool_heat) call rad_cool_heat_1D(f,df,p)

!f (accretion on NS)
!
    if (lnstar_entropy) then
   if (T_disk.EQ.0) then
     T_disk=cs0**2/gamma1
   endif 
  
 
       if ( dt .GT. 0..AND. n .GT. 24 .AND. n .LT. nzgrid-20) then
   
         if (lnstar_T_const) then
    
          df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)-1./(dt)*(p%TT(:)-T_disk)/T_disk
           !    df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
           !   -1./(5.*dt)*(f(l1:l2,m,n,iss)*gamma+gamma1*f(l1:l2,m,n,ilnrho))/p%rho(:)/p%TT(:)


        else
       
        endif
 
    
       endif 


      if (ldecelerat_zone) then
    
         if ( dt .GT. 0..AND. n .LE. ac_dc_size+4 ) then
          if (lnstar_T_const) then
          df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
           -1./(2.*dt)*(f(l1:l2,m,n,iss)*gamma+gamma1*f(l1:l2,m,n,ilnrho))/p%rho(:)/T_disk    
          !  df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
          ! -1./(5.*dt)*(f(l1:l2,m,n,iss)-log(TT_cs0)/gamma)/p%rho(:)/T_disk
          !  df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
          ! -1./(5.*dt)*(p%TT(:)-T_disk)/T_disk

          else
              df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
           -1./(5.*dt)*(p%TT(:)-T_star)/T_star
         !   df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
         !  -1./(5.*dt)*(f(l1:l2,m,n,iss)*gamma+gamma1*f(l1:l2,m,n,ilnrho))/p%rho(:)/T_star!p%TT(:)

          !    df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
          !     -1./(2.*dt)*(f(l1:l2,m,n,iss)*gamma+gamma1*f(l1:l2,m,n,ilnrho))/p%rho(:)/T_star   
 
          !  df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
          !  -1./(5.*dt)*(f(l1:l2,m,n,iss)-log(T_star)/gamma)/p%rho(:)/T_star
   
          endif
        end if  


    !     endif 
     endif  
   
     if (laccelerat_zone) then
         if (n .GE. nzgrid-ac_dc_size  .AND. dt .GT.0.) then
                   
          if (nxgrid .LE.1) then
              if (lnstar_T_const) then   
                 df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
                 -1./(5.*dt)*(p%TT(:)-T_disk)/T_disk
              else  

              !    df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
              !   -1./(5.*dt)*(p%TT(:)-T_disk)/T_disk
              !    df(l1:H_disk_point+4,m,n,iss)=df(l1:H_disk_point+4,m,n,iss) &
               !   -1./(5.*dt)*(p%TT(1:H_disk_point)-T_disk)/T_disk
               ! df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
               !  -1./(5.*dt)*(f(l1:l2,m,n,iss)-log(T_disk)/gamma)
              ! df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
              ! -1./(5.*dt)*(f(l1:l2,m,n,iss)*gamma+gamma1*f(l1:l2,m,n-1,ilnrho))/p%rho(:)/T_disk
      

              endif
         else
         endif     
 

     endif 

     endif  
    endif

       if (lsurface_zone) then
          if ( dt .GT.0.) then
            l_sz=l2-5
            l_sz_1=nxgrid-5

          if (lnstar_1D) then   
     !       df(l_sz:l2,m,n,iss)=df(l_sz:l2,m,n,iss) &
     !       -1./(5.*dt)*(f(l_sz:l2,m,n,iss)-log(T_disk)/gamma) &
     !       /p%rho(l_sz_1:nxgrid)/p%TT(l_sz_1:nxgrid)  
          else

            do j=l_sz,l2   
             df(j,m,n,iss)=df(j,m,n,iss)&
               -1./(5.*dt)*(f(j,m,n,iss)-f(j-1,m,n,iss))
            enddo 

          endif

      !        df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
      !       -1./(5.*dt)*(p%TT(:)-T_disk)/T_disk
        
          !    df(l_sz:l2,m,n,iss)=df(l1:l2,m,n,iss) &
       !     -1./(1.*dt)*(f(l_sz:l2,m,n,iss)*gamma+gamma1*f(l_sz:l2,m,n,ilnrho))/ &
       !     p%rho(l_sz_1:nxgrid)/T_disk!p%TT(l_sz_1:nxgrid) 

      !   df(l_sz:l2,m,n,iss)=df(l_sz:l2,m,n,iss) &
      !      -1./(5.*dt)*(p%TT(l_sz_1:nxgrid)-T_disk)/T_disk

       
         endif
      endif

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,df,p

    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine special_boundconds(f,bc)
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
!      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      type (boundary_condition) :: bc
!
      select case (bc%bcname)
       case ('stp')
         select case (bc%location)
         case (iBC_X_TOP)
           call bc_BL_x(f,-1, bc)
         case (iBC_X_BOT)
           call bc_BL_x(f,-1, bc)
         case (iBC_Z_TOP)
           call bc_BL_z(f,-1, bc)
         case (iBC_Z_BOT)
           call bc_BL_z(f,-1, bc)
         endselect
         bc%done=.true.
      endselect

      if (NO_WARN) print*,f(1,1,1,1),bc%bcname
!
    endsubroutine special_boundconds
!***********************************************************************
!
!  PRIVATE UTITLITY ROUTINES
!
!***********************************************************************
    subroutine rad_cool_heat_1D(f,df,p)
!
!  heat conduction
!  Natalia (NS)
!   12-apr-06/axel: adapted from Wolfgang's more complex version
!

!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
      real, dimension (mx,my,mz,mvar) :: df
       real, dimension (nx) :: diffus_chi1
      real, dimension (nx) :: thdiff_1D
      real ::  beta
      integer :: l_sz, l_sz_1,j 

      intent(in) :: f,p
      intent(out) :: df
 
 


!   cooling in 1D case 
!
    if (l1D_cooling) then

      beta=beta_hand  !1e6

      thdiff_1D =-16./3.*sigmaSB/kappa_es*p%TT**4 &
                 *p%rho1*beta

      l_sz=l2-10
      l_sz_1=nxgrid-10 

      if (ldecelerat_zone) then
        if (n .GT. ac_dc_size+4)   df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff_1D
        if (headtt) print*,'calc_heatcond_diffusion: added thdiff_1D'
      else
        df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff_1D
       if (headtt) print*,'calc_heatcond_diffusion: added thdiff_1D'
   
       endif
    endif 
 
!   heating in 1D case 
!
    if (l1D_heating) then
!
!  commented out reference to nu for the time being.
!  This line should really be in the viscosity module,
!  but this is currently not used anyway.
!
 !    thdiff_1D =p%rho*nu*(1.5*f(l1:l2,m,n,iuy)/xyz0(3))**2
      thdiff_1D =p%rho*nu_for_1D*(1.5*f(l1:l2,m,n,iuy)/xyz0(3))**2
      
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff_1D

      if (headtt) print*,'calc_heatcond_diffusion: added thdiff_1D'
    endif 
 

!
    endsubroutine rad_cool_heat_1D
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
      real :: cs2_star
!

       call eoscalc(ilnrho_lnTT,log(rho_star),log(T_star), cs2=cs2_star)

      
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
      if (lsharp) then
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
      if (lsmooth) then
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

 !   do i=1,44!H_disk_point+4
        do i=1,H_disk_point_int+4
    !    f(i,:,:,ilnrho)=ln_ro_u+(1.-(xx(i,:,:)/H_disk)**2)
     f(i,:,:,ilnrho)=ln_ro_u+(1.-M_star/2./zz(i,:,:)**3*x(i)**2*gamma/cs2_star)
 
       enddo 

		   
		   
       do i=H_disk_point_int+5,mx
    !     do i=45,mx
	 
     !	f(i,:,:,ilnrho)=f(44,:,:,ilnrho)
	
    ! f(i,:,:,ilnrho)=f(H_disk_point+4,:,:,ilnrho)
     f(i,:,:,ilnrho)=ln_ro_u+(1.-M_star/2./zz(i,:,:)**3*x(H_disk_point_int+4)**2*gamma/cs2_star)
     
     
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
!***************************************************************
    subroutine entropy_step(f,xx,zz,T_star)
!Natalia
!Initialization of entropy in a case of the step-like distribution
 use EquationOfState

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx, zz
      real, dimension (nx) ::  lnrho, lnTT,ss
      integer :: step_width, step_length, mi,ni, li,  decel_zone
      real :: H_disk_min, L_disk_min, hdisk, ldisk, ll, T_star, const_tmp


      decel_zone=ac_dc_size+4

      hdisk=H_disk 
      ldisk=L_disk

      H_disk_min=Lxyz(1)/(nxgrid-1)
      L_disk_min=Lxyz(3)/(nzgrid-1)

      if (H_disk .GT. Lxyz(1)-H_disk_min) hdisk=Lxyz(1)
      if (H_disk .LT. H_disk_min) hdisk=0.

      if (L_disk .GT. Lxyz(3)-L_disk_min) ldisk=Lxyz(3)
      if (L_disk .LT. L_disk_min) ldisk=0.

      step_width=nint((nxgrid-1)*hdisk/Lxyz(1))
      step_length=nint((nzgrid-1)*(Lxyz(3)-ldisk)/Lxyz(3))


      lnTT=log(T_star)!  log(T0)
 if (T_disk.EQ.0) then
     T_disk=cs0**2/gamma1
   endif 
    

      print*,'T_star=',T_star
      do ni=n1,n2;
       do mi=m1,m2;
     if (lnstar_T_const) then
       f(l1:l2,mi,ni,iss)=-f(l1:l2,mi,ni,ilnrho)*gamma1/gamma
      else
      ! lnrho=f(l1:l2,mi,ni,ilnrho)
      ! const_tmp=M_star/sigmaSB*c_light*3./4.

     !  lnTT=0.25*log(T_star**4+const_tmp*exp(f(l1:l2,mi,ni,ilnrho))*(1./zz(l1:l2,mi,ni)-1./R_star))
    
     !  call eoscalc(4,lnrho,lnTT,ss=ss)
  
      ! f(l1:l2,mi,ni,iss)=ss  

     if (nxgrid .LE. 1) then
        f(l1:l2,mi,ni,iss)=-f(l1:l2,mi,ni,ilnrho)*gamma1/gamma
     else
         f(l1:l2,mi,ni,iss)=-f(l1:l2,mi,ni,ilnrho)*gamma1/gamma

     !  lnrho=f(l1:l2,mi,ni,ilnrho)
    
   
     !  lnTT=log(T_disk)
   
     !  call eoscalc(4,lnrho,lnTT,ss=ss)
  
     !  f(l1:l2,mi,ni,iss)=ss  


   !    f(l1,mi,ni,iss)=-f(l1,mi,ni,ilnrho)*gamma1/gamma
   
   !     do li=l1+1,l2;
   !         lnrho=f(li,mi,ni,ilnrho)
   !         lnTT=(xx(li,mi,ni)-0.)/Lxyz(1)*(log(TT_cs0)-log(TT_cs0*0.01))+log(TT_cs0*0.01)
   !         call eoscalc(4,lnrho,lnTT,ss=ss)
   !         f(l1:l2,mi,ni,iss)=ss  
   !    enddo

     endif 



     endif

       end do 
     end do   

     !   f(l1:l2,:,:,iss)=lnTT0/gamma

    !    f(:,:,:,iss)=(zz(:,:,:)-R_star)/Lxyz(3)*(lnTT0-lnTT0/10.)/gamma+lnTT0/10./gamma

    endsubroutine entropy_step
!***********************************************************************
    subroutine velocity_step(f)
!Natalia
!Initialization of velocity in a case of the step-like distribution

      use Cdata

      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: step_width, step_length
      real ::  H_disk_min, L_disk_min, hdisk, ldisk
    
      hdisk=H_disk 
      ldisk=L_disk
    
      H_disk_min=Lxyz(1)/(nxgrid-1)
      L_disk_min=Lxyz(3)/(nzgrid-1)

      if (H_disk .GT. Lxyz(1)-H_disk_min) hdisk=Lxyz(1)
      if (H_disk .LT. H_disk_min) hdisk=0.

      if (L_disk .GT. Lxyz(3)-L_disk_min) ldisk=Lxyz(3)
      if (L_disk .LT. L_disk_min) ldisk=0.


       
      step_width=nint((nxgrid-1)*hdisk/Lxyz(1))
      step_length=nint((nzgrid-1)*(Lxyz(3)-ldisk)/Lxyz(3))


      if (hdisk .EQ. Lxyz(1) .AND. ldisk .EQ. Lxyz(3))  then
        f(:,:,:,iuz)=uu_left
        f(:,:,:,iuy)=uy_left
      end if

      if (hdisk .EQ. 0. .AND. ldisk .EQ. 0.) f(:,:,:,iuy)=uy_right
      if (hdisk .EQ. 0. .OR. ldisk .EQ. 0.) f(:,:,:,iuy)=uy_right
      
       
      if (hdisk .EQ. Lxyz(1) .AND. ldisk .LT. Lxyz(3)) then
        f(:,:,step_length+3+1:mz,iuz)=uu_left
        f(:,:,step_length+3+1:mz,iuy)=uy_left
        f(:,:,1:step_length+3,iuy)=uy_right
      endif

      if (hdisk .LT. Lxyz(1) .AND. ldisk .EQ. Lxyz(3)) then
        f(1:step_width+3,:,:,iuz)=uu_left
        f(1:step_width+3,:,:,iuy)=uy_left
        f(step_width+3+1:mx,:,:,iuy)=uy_right
      endif


      if (hdisk .GT. 0.  .AND. hdisk .LT. Lxyz(1) ) then
        if (ldisk .GT. 0.  .AND. ldisk .LT. Lxyz(3)) then
          f(1:step_width+3,:,step_length+3+1:mz,iuz)=uu_left
          f(1:step_width+3,:,step_length+3+1:mz,iuy)=uy_left
          f(step_width+3+1:mx,:,step_length+3+1:mz,iuy)=uy_right
          f(:,:,1:step_length+3,iuy)=uy_right
        end if
      end if

    endsubroutine  velocity_step
!***********************************************************************
    subroutine velocity_kep_disk(f,zz)
!Natalia
!Initialization of velocity in a case of the step-like distribution

      use Cdata

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: zz
      integer :: step_length, decel_zone
      real ::   L_disk_min,  ldisk,   ll

      decel_zone=ac_dc_size+4
       

      ll=Lxyz(3)-L_disk
      ldisk=L_disk

   if (nzgrid .GT. 1) then
      L_disk_min=Lxyz(3)/(nzgrid-1)

      if (L_disk .GT. Lxyz(3)-L_disk_min) ldisk=Lxyz(3)
      if (L_disk .LT. L_disk_min) ldisk=0.
     
      step_length=nint((nzgrid-1)*ll/Lxyz(3))

      if (ldisk .LT. L_disk_min) then
        f(:,:,:,iuz)=uu_left
        f(:,:,:,iuy)=(zz-R_star)/ll*sqrt(M_star/(ll+R_star))
      endif

      if (ldisk .GE. Lxyz(3)) then
        f(:,:,:,iuz)=uu_left
        f(:,:,:,iuy)=sqrt(M_star/zz)
      endif
!
      if (ldisk .GT. 0.  .AND. ldisk .LT. Lxyz(3)) then

        f(:,:,:,iuz)=uu_left
       if (ldecelerat_zone .AND. decel_zone .LT. nzgrid) then 
 
        f(:,:,step_length+3+1:mz,iuy)=sqrt(M_star/zz(:,:,step_length+3+1:mz))

        f(:,:,decel_zone+1:step_length+3,iuy)= &
           (zz(:,:,decel_zone+1:step_length+3)-R_star-(decel_zone-4)*L_disk_min) &
           /(ll-(decel_zone-4)*L_disk_min)*sqrt(M_star/(ll+R_star))
        f(:,:,1:decel_zone,iuy)=0.

       else
         f(:,:,step_length+3+1:mz,iuy)=sqrt(M_star/zz(:,:,step_length+3+1:mz))
         f(:,:,1:step_length+3,iuy)= &
           (zz(:,:,1:step_length+3)-R_star)/ll*sqrt(M_star/(ll+R_star))
       end if

      end if
   else
    f(:,:,:,iux)=uu_left
    f(:,:,:,iuz)=uu_left
    f(:,:,:,iuy)=sqrt(M_star/xyz0(3))
   endif


    endsubroutine  
!***********************************************************************
    subroutine bc_BL_x(f,sgn,bc)
!
! Natalia
!  11-may-06
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: sgn
      type (boundary_condition) :: bc
      integer :: i,j

      j=bc%ivar
      if (bc%location==iBC_X_BOT) then
      ! bottom boundary
        if (j == 1) then 
            f(l1,:,:,j) = 0.
            do i=1,nghost; f(l2+i,:,:,j)=2*f(l2,:,:,j)+sgn*f(l2-i,:,:,j); enddo
        else
           do i=1,nghost; f(l1-i,:,:,j)= f(l1+i,:,:,j); enddo
        endif
          !   f(l1,:,:,j) = 0. ! set bdry value=0 (indep of initcond)

      elseif (bc%location==iBC_X_TOP) then
      ! top boundary
        if (nxgrid <= 1) then
          if (j == 1) then
            f(l2,m1:m2,n1:n2,j)=bc%value1
            do i=1,nghost; f(l2+i,:,:,j)=2*f(l2,:,:,j)+sgn*f(l2-i,:,:,j); enddo
          else
            f(l2+1,:,:,j)=0.25*(  9*f(l2,:,:,j)- 3*f(l2-1,:,:,j)- 5*f(l2-2,:,:,j)+ 3*f(l2-3,:,:,j))
            f(l2+2,:,:,j)=0.05*( 81*f(l2,:,:,j)-43*f(l2-1,:,:,j)-57*f(l2-2,:,:,j)+39*f(l2-3,:,:,j))
            f(l2+3,:,:,j)=0.05*(127*f(l2,:,:,j)-81*f(l2-1,:,:,j)-99*f(l2-2,:,:,j)+73*f(l2-3,:,:,j))
          endif
        else
!            f(l2+1,:,:,j)=0.25*(  9*f(l2,:,:,j)- 3*f(l2-1,:,:,j)- 5*f(l2-2,:,:,j)+ 3*f(l2-3,:,:,j))
!            f(l2+2,:,:,j)=0.05*( 81*f(l2,:,:,j)-43*f(l2-1,:,:,j)-57*f(l2-2,:,:,j)+39*f(l2-3,:,:,j))
!            f(l2+3,:,:,j)=0.05*(127*f(l2,:,:,j)-81*f(l2-1,:,:,j)-99*f(l2-2,:,:,j)+73*f(l2-3,:,:,j))
!
          do i=1,nghost; f(l1+i,:,:,j)=f(l1-i,:,:,j); enddo
!
        endif
      else
        print*, "bc_BL_x: ", bc%location, " should be `top(", &
                        iBC_X_TOP,")' or `bot(",iBC_X_BOT,")'"
      endif
!
    endsubroutine bc_BL_x
 !***********************************************************************
    subroutine bc_BL_z(f,sgn,bc)
!
!  Step boundary conditions.
!
!  11-feb-06/nbabkovs
!
      use Cdata
      use EquationOfState
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: value1,value2 
      type (boundary_condition) :: bc
      real, dimension(nx) :: lnrho,lnTT,ss
      integer :: sgn,i,j, step_width, n1p4,n2m4, i_tmp
      real :: H_disk_min, L_disk_min, ddz, ddx
    !  integer, parameter :: ilnrho_lnTT=4

      j=bc%ivar
      H_disk_min=Lxyz(1)/(nxgrid-1)
      step_width=nint((nxgrid-1)*H_disk/Lxyz(1))
      ddx=H_disk_min

      L_disk_min=Lxyz(3)/(nzgrid-1)
      ddz=L_disk_min

      if (j == 4 .or. j==5) then
        value1=log(bc%value1)
        value2=log(bc%value2)
      else
        value1=bc%value1
        value2=bc%value2
      endif

      if (bc%location==iBC_Z_BOT) then
      ! bottom boundary
        if (lextrapolate_bot_density .and. j>=4) then
          n1p4=n1+4

          f(:,:,n1-1,j)=0.2   *(  9*f(:,:,n1,j)                 -  4*f(:,:,n1+2,j)- 3*f(:,:,n1+3,j)+ 3*f(:,:,n1p4,j))
          f(:,:,n1-2,j)=0.2   *( 15*f(:,:,n1,j)- 2*f(:,:,n1+1,j)-  9*f(:,:,n1+2,j)- 6*f(:,:,n1+3,j)+ 7*f(:,:,n1p4,j))
          f(:,:,n1-3,j)=1./35.*(157*f(:,:,n1,j)-33*f(:,:,n1+1,j)-108*f(:,:,n1+2,j)-68*f(:,:,n1+3,j)+87*f(:,:,n1p4,j))
        else
          if (j==5) then
            lnrho=f(l1:l2,m1,n1,ilnrho)
            if (lnstar_T_const) then 
              lnTT=log(cs0**2/(gamma1))
            else     
              lnTT=log(T_star)
            endif
            !+ other terms for sound speed not equal to cs_0
            call eoscalc(4,lnrho,lnTT,ss=ss)
            f(l1:l2,m1,n1,iss)=ss 
            !  print*, 'boundary entropy ', ss
            !ss=exp(ss-(-log(cs0**2/(gamma1))-gamma1*lnrho)/gamma)
            !   ss=exp(log(cs0**2/(gamma1))+gamma*ss+gamma1*lnrho)
            !print*, 'boundary entropy ', ss
          else
            if (H_disk >= H_disk_min .and. H_disk <= Lxyz(1)-H_disk_min) then
              f(1:step_width+3,:,n1,j)=value1
              f(step_width+3+1:mx,:,n1,j)=value2
            endif
            if (H_disk < H_disk_min)    f(:,:,n1,j)=value2
            if (H_disk > Lxyz(1)-H_disk_min)    f(:,:,n1,j)=value1
          endif
          do i=1,nghost; f(:,:,n1-i,j)=2*f(:,:,n1,j)+sgn*f(:,:,n1+i,j); enddo
        endif
      elseif (bc%location==iBC_Z_TOP) then
      ! top boundary

        if (ltop_velocity_kep .and. j==2) then 
          f(:,:,n2,j)=sqrt(M_star/(R_star+Lxyz(3)))
        else
          if (nxgrid <= 1) then
            if (j==5) then
              lnrho=f(l1:l2,m2,n2,ilnrho)
              if (T_disk.EQ.0) then    
              lnTT=log(cs0**2/(gamma1))
              else
              lnTT=log(T_disk)    
              endif           
              call eoscalc(4,lnrho,lnTT,ss=ss)
              f(l1:l2,m2,n2,iss)=ss
            else 
              if (H_disk >= H_disk_min .and. H_disk <= Lxyz(1)-H_disk_min) then
                f(1:step_width+3,:,n2,j)=value1
                f(step_width+3+1:mx,:,n2,j)=value2
              endif
              if (H_disk < H_disk_min)    f(:,:,n2,j)=value2
              if (H_disk > Lxyz(1)-H_disk_min)    f(:,:,n2,j)=value1
            endif
          else

!           if (j==4) then
!             f(1:H_disk_point+4,:,n2,j)=value1
!             f(H_disk_point+5:mx,:,n2,j)=value2
!
!             do i=1,H_disk_point+4
!               f(i,:,n2,ilnrho)=log(5.)+(1.-(ddx*i/H_disk)**2)
!             enddo 
!
!             do i=H_disk_point+5,mx
!               f(i,:,n2,ilnrho)=f(H_disk_point+4,:,n2,ilnrho)
!             enddo    
!
!           else 
            n2m4=n2-4
            i_tmp=H_disk_point+5

            f(i_tmp:mx,:,n2+1,j)=0.2   *(   9*f(i_tmp:mx,:,n2  ,j) &
                                         -  4*f(i_tmp:mx,:,n2-2,j) &
                                         -  3*f(i_tmp:mx,:,n2-3,j) &
                                         +  3*f(i_tmp:mx,:,n2m4,j))
            f(i_tmp:mx,:,n2+2,j)=0.2   *(  15*f(i_tmp:mx,:,n2  ,j) &
                                        -   2*f(i_tmp:mx,:,n2-1,j) &
                                        -   9*f(i_tmp:mx,:,n2-2,j) &
                                        -   6*f(i_tmp:mx,:,n2-3,j) &
                                        +   7*f(i_tmp:mx,:,n2m4,j))
            f(i_tmp:mx,:,n2+3,j)=1./35.*( 157*f(i_tmp:mx,:,n2  ,j) &
                                         - 33*f(i_tmp:mx,:,n2-1,j) &
                                         -108*f(i_tmp:mx,:,n2-2,j) &
                                         - 68*f(i_tmp:mx,:,n2-3,j) &
                                         + 87*f(i_tmp:mx,:,n2m4,j))
!           endif

          endif
        endif

        do i=1,nghost; f(:,:,n2+i,j)=2*f(:,:,n2,j)+sgn*f(:,:,n2-i,j); enddo

      else 
        print*, "bc_BL_z: ", bc%location, " should be `top(", &
                        iBC_X_TOP,")' or `bot(",iBC_X_BOT,")'"
      endif
!
    endsubroutine bc_BL_z
!***********************************************************************
    subroutine special_before_boundary(f)
!
!   Possibility to modify the f array before the boundaries are 
!   communicated.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-jul-06/tony: coded
!
      use Cdata
!      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
!
      if (NO_WARN) print*,f(1,1,1,1)
!
    endsubroutine special_before_boundary
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include 'special_dummies.inc'
!********************************************************************
endmodule Special

