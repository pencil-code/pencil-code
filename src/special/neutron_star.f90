! $Id: neutron_star.f90,v 1.41 2006-12-11 17:09:30 nbabkovs Exp $
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
  use Sub, only: keep_compiler_quiet
!  use Density, only: rho_up
  use EquationOfState


  implicit none

  include 'special.h'

  ! input parameters
  logical :: lmass_source_NS=.false.
  logical :: leffective_gravity=.false.

  character (len=labellen) :: initnstar='default'
  real :: rho_star=1.,rho_disk=1., rho_surf=1.

  real :: uu_init=0.


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
  integer :: L_disk_point=0

  real :: beta_hand=1.
  real :: nu_for_1D=1.
  real :: star_rot=0.
  logical :: hot_star=.false.
  logical :: new_T_profile=.false.

 logical :: l1D_cooling=.false.,l1D_heating=.false.

  logical :: lraddif_local=.false.



  logical :: l1D_cool_heat=.false.
  logical :: lgrav_x_mdf=.false.
!
! Keep some over used pencils
!
  real, dimension(nx) :: z_2


! start parameters
  namelist /neutron_star_init_pars/ &
      initnstar,lmass_source_NS,leffective_gravity, &
      laccelerat_zone, ldecelerat_zone, lsurface_zone, &
      rho_star,rho_disk, rho_surf,&
      lraddif_local,&
      L_disk_point, R_star, M_star, &
      T_star,accretion_flux, T_disk, &
      uu_init, &
      l1D_cooling,l1D_heating,beta_hand, &
      nu_for_1D, ltop_velocity_kep, lextrapolate_bot_density, &
      lnstar_entropy, lnstar_T_const, lnstar_1D, &
      lgrav_x_mdf, star_rot, hot_star, new_T_profile

! run parameters
  namelist /neutron_star_run_pars/ &
      lmass_source_NS,leffective_gravity, rho_star,rho_disk, &
      laccelerat_zone, ldecelerat_zone, lsurface_zone,lraddif_local, &
       accretion_flux, lnstar_entropy, &
       lnstar_T_const,lnstar_1D, &
       l1D_cooling,l1D_heating, lgrav_x_mdf, star_rot, hot_star, new_T_profile
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
!  CVS should automatically update everything between $Id: neutron_star.f90,v 1.41 2006-12-11 17:09:30 nbabkovs Exp $
!  when the file in committed to a CVS repository.
!
      if (lroot) call cvs_id( &
           "$Id: neutron_star.f90,v 1.41 2006-12-11 17:09:30 nbabkovs Exp $")
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
          call density_init(f,xx,zz)
          call entropy_init(f,xx,zz)
          call velocity_init(f,zz)
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
     if (lraddif_local) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_gss)=.true.
        lpenc_requested(i_del2lnrho)=.true.
        lpenc_requested(i_del2ss)=.true.
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
!   06-oct-03/tony: coded
!
      use Cdata
      ! use Viscosity
      use EquationOfState

      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      integer :: i, l_sz, tmp_int
      real :: cs2_star, p_gas,p_rad, Sigma_rho, grad_rho

       call eoscalc(ilnrho_lnTT,log(rho_disk),log(T_disk), cs2=cs2_star)

!  mass sources and sinks for the boundary layer on NS in 1D approximation
!
      if (lmass_source_NS) call mass_source_NS(f,df,p%rho)

       if (laccelerat_zone) then
         if (n .GE. nzgrid-ac_dc_size .AND. dt .GT. 0.) then
          if (lnstar_entropy) then
            if (nxgrid == 1) then
              if (lnstar_T_const) then
              endif
            else
             l_sz=l2-10
	if (hot_star) then
	 df(l1:l_sz,m,n,ilnrho)=df(l1:l_sz,m,n,ilnrho) &
	 -1./(5.*dt)*(f(l1:l_sz,m,n,ilnrho) &
	  -log(rho_disk)-(-M_star/2./z(n)**3 &
	     *x(l1:l_sz)**2*gamma/cs2_star))
	     
!	df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) &
!	 -1./(5.*dt)*(f(l1:l2,m,n,ilnrho) -f(l1:l2,m,n-1,ilnrho)) 
	 
	else   
               df(l1:l_sz,m,n,ilnrho)=df(l1:l_sz,m,n,ilnrho) &
                  -1./(5.*dt)*(f(l1:l_sz,m,n,ilnrho) &
                  -log(rho_disk)-(-M_star/2./z(n)**3 &
                  *x(l1:l_sz)**2*gamma/(p%cs2(l1:l_sz))))
        endif
 

               df(l_sz+1:l2,m,n,ilnrho)=df(l_sz+1:l2,m,n,ilnrho) &
	         -1./(5.*dt)*(f(l_sz+1:l2,m,n,ilnrho)-f(l_sz,m,n,ilnrho))				 
  
            endif
          endif
         endif
       endif

       if (ldecelerat_zone) then
         if (n .LE. ac_dc_size+4 .AND. dt .GT. 0.) then

            if (lnstar_entropy) then
              if (lnstar_T_const) then
               df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) &
                -1./p%rho(:)/(5.*dt)*(p%rho(:)-rho_star &
                *exp(-M_star/R_star/cs0**2*gamma*(1.-R_star/z(n))))
              endif

            if (hot_star) then
	     df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) &
	       -1./(5.*dt)*(f(l1:l2,m,n,ilnrho)-f(l1:l2,m,n+1,ilnrho))
	    endif
	    
            else


            endif
         endif
       endif


! surface zone in a case of a Keplerian disk

      if (lsurface_zone) then
         if ( dt .GT.0.) then
            l_sz=l2-5
          if (lnstar_1D) then
           else
           do i=l_sz,l2
            if (n .LT. nzgrid-ac_dc_size .AND. dt .GT. 0.) then

            if (hot_star) then
	     df(i,m,n,ilnrho)=df(i,m,n,ilnrho)&
	     -1./(5.*dt)*(f(i,m,n,ilnrho)-f(i-1,m,n,ilnrho))
	    else
             df(i,m,n,ilnrho)=df(i,m,n,ilnrho)&
             -1./(5.*dt)*(f(i,m,n,ilnrho)-f(i-1,m,n,ilnrho) &
             +M_star/z(n)**3*(x(i)-x(i-1))*x(i-1)*gamma/p%cs2(i-1))
            endif
            else
            df(i,m,n,ilnrho)=df(i,m,n,ilnrho)&
             -1./(5.*dt)*(f(i,m,n,ilnrho)-f(i-1,m,n,ilnrho))
            endif
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

       l_sz=l2-5*0

      if (leffective_gravity) then
        if (headtt) &
          print*,'duu_dt: Effectiv gravity; Omega, Rstar=', Omega, R_star, M_star

          df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)- &
            M_star/z(n)**2*(1.-p%uu(:,2)*p%uu(:,2)*z(n)/M_star)

           df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+ &
	               +p%uu(:,1)*p%uu(:,2)/z(n)
		       

        if (nxgrid /= 1) then

          if (lgrav_x_mdf) then
            df(l1:l_sz,m,n,iux)=df(l1:l_sz,m,n,iux)- &
              M_star/z(n)**2/sqrt(z(n)**2+x(l1:l_sz)**2)*x(l1:l_sz)*(z(n)-R_star)/(Lxyz(1)*0.5)
           else
            df(l1:l_sz,m,n,iux)=df(l1:l_sz,m,n,iux)-M_star/z(n)**3*x(l1:l_sz)
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

         l_sz=l2-5
	 
	 if (hot_star) then
	  df(l1:l_sz,m,n,iux)=df(l1:l_sz,m,n,iux)&
	    -1./(5.*dt)*(f(l1:l_sz,m,n,iux)-f(l1:l_sz,m,n-1,iux))
			 
	 else
             df(l1:l_sz,m,n,iux)=df(l1:l_sz,m,n,iux) &
                              -1./(5.*dt)*(f(l1:l_sz,m,n,iux)-0.)
         endif

        !    df(l_sz+1:l2,m,n,iux)=df(l_sz+1:l2,m,n,iux) &
        !     -1./(5.*dt)*(f(l_sz+1:l2,m,n,iux)-f(l_sz+1:l2,m,n-1,iux))


         l_sz=l2-5

        ! if (hot_star) then
	! else
             df(l1:l_sz,m,n,iuz)=df(l1:l_sz,m,n,iuz)&
              -1./(5.*dt)*(f(l1:l_sz,m,n,iuz)-f(l1:l_sz,m,n-1,iuz))
        !    df(l1:l_sz,m,n,iux)=df(l1:l_sz,m,n,iux)&
	!       -1./(5.*dt)*(f(l1:l_sz,m,n,iux)-f(l1:l_sz,m,n-1,iux))
			  

        ! endif
           endif
        endif
      endif
!
! deceleration zone in a case of a Keplerian disk
!
      if (ldecelerat_zone) then
        if (n .le. ac_dc_size+4  .and. dt .gt.0.) then
          if (lnstar_entropy) then

            df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)&
                -1./(5.*dt)*(p%uu(:,2)-sqrt(M_star/z(n))*star_rot)
		
!	if (hot_star) then
!	else	
            df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)&
                -1./(5.*dt)*(p%uu(:,3)-0.)
		
	if (hot_star) then	
	      df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)&
	        -1./(5.*dt)*(p%uu(:,1)-0.)
	endif	      	
		
!	endif
		
	!if (hot_star) then	
	!     df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)&
	!     -1./(5.*dt)*(f(l1:l2,m,n,iux)-f(l1:l2,m,n-1,iux))	
        ! endif
          else
           df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux) &
                    -1./(5.*dt)*(p%uu(:,1)-0.)


           df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)&
                    -1./(5.*dt)*(p%uu(:,2)-0.)

            df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)&
                    -1./(5.*dt)*(p%uu(:,3)-0.)
         endif
        endif
      endif

! surface zone in a case of a Keplerian disk

      if (lsurface_zone) then
        if ( dt .gt.0.) then
!
          l_sz=l2-5

         if (lnstar_1D) then
             df(l_sz:l2,m,n,iux)=df(l_sz:l2,m,n,iux)&
                   -1./(2.*dt)*(f(l_sz:l2,m,n,iux)-0.)
            df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy) &
                  -1./(5.*dt)*(f(l1:l2,m,n,iuy) &
                  -sqrt(M_star/xyz0(3)))
            df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)&
                  -1./(5.*dt)*(f(l1:l2,m,n,iuz)&
                  +accretion_flux/p%rho(:))
         else
	 
	 
	
!	if (hot_star) then
	
!	if (n.gt.ac_dc_size+4) then
!	      do j=l_sz,l2
	
!	        df(j,m,n,iux)=df(j,m,n,iux)&
!	          -1./(5.*dt)*(f(j,m,n,iux)-f(j-1,m,n,iux))
!	        df(j,m,n,iuy)=df(j,m,n,iuy)&
!	          -1./(5.*dt)*(f(j,m,n,iuy)-f(j-1,m,n,iuy))
!         	df(j,m,n,iuz)=df(j,m,n,iuz)&
!	          -1./(5.*dt)*(f(j,m,n,iuz)-f(j-1,m,n,iuz))
!          	enddo
!        endif			      
!	else
           do j=l_sz,l2
             df(j,m,n,iux)=df(j,m,n,iux)&
                  -1./(5.*dt)*(f(j,m,n,iux)-f(j-1,m,n,iux))
             df(j,m,n,iuy)=df(j,m,n,iuy)&
                  -1./(5.*dt)*(f(j,m,n,iuy)-f(j-1,m,n,iuy))
             df(j,m,n,iuz)=df(j,m,n,iuz)&
                  -1./(5.*dt)*(f(j,m,n,iuz)-f(j-1,m,n,iuz))
           enddo
!        endif 

         endif

       endif
      endif

    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
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

      use Cdata

      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      integer :: j, l_sz, l_sz_1, li
      real :: dT_dx_i1

        if (l1D_cool_heat) call rad_cool_heat_1D(f,df,p)

        if (lraddif_local) call raddif_local(f,df,p)

    if (lnstar_entropy) then


      if (T_disk.EQ.0) then
         T_disk=cs0**2/gamma1
      endif


      if ( dt .GT. 0..AND. n .GT. 24 .AND. n .LT. nzgrid-20) then
         if (lnstar_T_const) then
           df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
                   -1./(dt)*(p%TT(:)-T_disk)/T_disk
         endif
      endif


      if (ldecelerat_zone) then

         if ( dt .GT. 0..AND. n .LE. ac_dc_size+4 ) then
          if (lnstar_T_const) then
           df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
            -1./(2.*dt)*(f(l1:l2,m,n,iss) &
            *gamma+gamma1*f(l1:l2,m,n,ilnrho))/p%rho(:)/T_disk

          else


          endif
         endif
      endif

     if (laccelerat_zone) then
         if (n .GE. nzgrid-ac_dc_size  .AND. dt .GT.0.) then
           if (nxgrid .LE.1) then
              if (lnstar_T_const) then
                 df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
                 -1./(5.*dt)*(p%TT(:)-T_disk)/T_disk
              endif
           else
     
          if (hot_star .and. new_T_profile) then
               
            df(l1,m,n,iss)=df(l1,m,n,iss) &
                 -1./(5.*dt)*(p%TT(1)-T_disk)/p%TT(1)
   
            do li=l1+1,l2      
  
            dT_dx_i1=-M_star/z(n)**3*x(li-l1) &
              *3./16./sigmaSB*c_light*p%rho(li-l1)/p%TT(li-l1)**3

            df(li,m,n,iss)=df(li,m,n,iss) &
            -1./(5.*dt)*(p%TT(li-l1+1)-p%TT(li-l1)-dT_dx_i1*dx)/p%TT(li-l1+1)
            enddo 

          endif 


           endif

          endif

     endif


       if (lsurface_zone) then
          if ( dt .GT.0.) then
            l_sz=l2-5
            l_sz_1=nxgrid-5

            if (lnstar_1D) then
             else
              if (hot_star .and. new_T_profile) then
                
                df(l1,m,n,iss)=df(l1,m,n,iss) &
                   -1./(5.*dt)*(p%TT(1)-T_disk)/p%TT(1)
   
                do li=l_sz,l2      
  
                 dT_dx_i1=-M_star/z(n)**3*x(li-l1) &
                          *3./16./sigmaSB*c_light  &
                          *p%rho(li-l1)/p%TT(li-l1)**3

                 df(li,m,n,iss)=df(li,m,n,iss) &
                   -1./(5.*dt)*(p%TT(li-l1+1)-p%TT(li-l1)-dT_dx_i1*dx)/p%TT(li-l1+1)
                enddo 
              else       
               do j=l_sz,l2
                 df(j,m,n,iss)=df(j,m,n,iss)*0.&
                 -1./(5.*dt)*(f(j,m,n,iss)-f(j-1,m,n,iss))
               enddo
              endif  
            endif
          endif
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

  !     print*,' cooling  ',thdiff_1D
      endif
    endif

!   heating in 1D case
!
    if (l1D_heating) then

      thdiff_1D =p%rho*nu_for_1D*(1.5*f(l1:l2,m,n,iuy)/xyz0(3))**2

      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff_1D


    !   if (lsurface_zone) then
       !    df(l2-l_sz:l2,m,n,iss) = df(l2-l_sz:l2,m,n,iss) + thdiff_1D(l_sz_1:nxgrid)

    !     df(l1:l2-l_sz,m,n,iss) = df(l1:l2-l_sz,m,n,iss) + thdiff_1D(1:l_sz_1)
    !   else
    !       df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff_1D
    !   endif

  ! print*,'heating',thdiff_1D

  endif

!
    endsubroutine rad_cool_heat_1D
!*************************************************************************
 subroutine raddif_local(f,df,p)
!
!  heat conduction
!  Natalia (NS)
!   12-apr-06/axel: adapted from Wolfgang's more complex version
!
      use Sub, only: max_mn_name,dot
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: glnT,glnThcond !,glhc
      real, dimension (nx) :: chix,diffus_chi1
      real, dimension (nx) :: thdiff,g2,thdiff_1D
      real, dimension (nx) :: hcond
      real ::  beta
      integer :: l_sz, l_sz_1

      intent(in) :: f,p
      intent(out) :: df

!
!  Heat conduction
!
      chix = p%rho1*p%rho1*p%TT**3*16./3.*sigmaSB/kappa_es!hcond
      glnT = gamma*p%gss + gamma1*p%glnrho ! grad ln(T)
      glnThcond = glnT !... + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
      call dot(glnT,glnThcond,g2)
!
!AB:  derivs of chix missing??
!
      thdiff = chix * (gamma*p%del2ss+gamma1*p%del2lnrho + g2)
!print*,' p%del2ss(10)'   ,p%del2ss(10)

   !  add heat conduction to entropy equation
    !


!All calculations till 1.10.2006 were made with this option
! Check it
        if (ldecelerat_zone) then
          if (nxgrid == 1) then
           if (n .gt. ac_dc_size+4) df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
          else
           if (n .gt. ac_dc_size+4) df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
          endif
         if (headtt) print*,'calc_heatcond_diffusion: added thdiff'
        else
         df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
         if (headtt) print*,'calc_heatcond_diffusion: added thdiff'
        endif

!


      if (headtt) print*,'calc_heatcond_diffusion: added thdiff_1D'



        df(l1:l2,m,n,iuz) = &
         df(l1:l2,m,n,iuz)-p%rho1*16./3.*sigmaSB/c_light*p%TT**4*glnT(:,3)

        df(l1:l2,m,n,iuy) = &
          df(l1:l2,m,n,iuy)-p%rho1*16./3.*sigmaSB/c_light*p%TT**4*glnT(:,2)

        df(l1:l2,m,n,iux) = &
         df(l1:l2,m,n,iux)-p%rho1*16./3.*sigmaSB/c_light*p%TT**4*glnT(:,1)
!
!  include constraint from radiative time step
!
      if (lfirst.and.ldt) then
        advec_crad2=p%rho1*16./3.*sigmaSB/c_light*p%TT**4
      endif

     if (headtt) print*,'calc_radiation_pressure: added to z-component'
!
!  check maximum diffusion from thermal diffusion
!  With heat conduction, the second-order term for entropy is
!  gamma*chix*del2ss
!
      if (lfirst.and.ldt) then
! Calculate timestep limitation
        diffus_chi=max(diffus_chi,gamma*chix*dxyz_2)

     !   diffus_chi1=min(gamma*chix*dxyz_2, &
      !              real(sigmaSB*kappa_es*p%TT**3*4.*p%cp1tilde))
     !   diffus_chi=max(diffus_chi,diffus_chi1)
      endif
    endsubroutine raddif_local
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

       sink_area=dz*sink_area_points

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
!      ksi=2.*((Lxyz(3)/(nzgrid-1.)*(sink_area_points+4-n))/sink_area)/sink_area
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
!********************************************************************
     subroutine density_init(f,xx,zz)
     !
     !Natalia
     !Initialization of density in a case of the step-like distribution
     !
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx, zz
      real ::  ln_ro_l, ln_ro_r, ln_ro_u, cs2_disk
      integer :: i

       call eoscalc(ilnrho_lnTT,log(rho_disk),log(T_disk), cs2=cs2_disk)


        ln_ro_r=log(rho_disk)
        ln_ro_l=log(rho_star)
        ln_ro_u=log(rho_surf)

       if (nxgrid/=1.and.nzgrid/=1) then

          f(:,:,:,ilnrho)= &
            -M_star/2./zz(:,:,:)**3*xx(:,:,:)**2*gamma/cs2_disk


        f(:,:,:,ilnrho)= f(:,:,:,ilnrho) &
                       +((zz(:,:,:)-R_star)/Lxyz(3))**0.25 &
                       *(ln_ro_r-ln_ro_l)+ln_ro_l
       else
         if (nzgrid .GT. 1) then
            f(:,:,:,ilnrho)=((zz(:,:,:)-R_star)/Lxyz(3))**0.25 &
                           *(ln_ro_r-ln_ro_l)+ln_ro_l
         else
            f(:,:,:,ilnrho)=(xx(:,:,:)-0.)/Lxyz(1) &
                           *(ln_ro_u-ln_ro_r)+ln_ro_r
         endif
       endif

      endsubroutine density_init
!***************************************************************

      subroutine entropy_init(f,xx,zz)

     ! use EquationOfState

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx, zz
      real, dimension (nx) ::  lnrho, lnTT,ss, TT_0
      integer ::  mi,ni,li,i
      real ::  ll, cs2_star, dT_dx_i1


      if (T_star.GT.0)  lnTT=log(T_star)
        print*,'T_star=',T_star

      if (T_disk.EQ.0) then
         T_disk=cs0**2/gamma1
      endif

         !cs2_star=sqrt(T_star*gamma1)
         call eoscalc(ilnrho_lnTT,log(rho_disk),log(T_disk), cs2=cs2_star)

      do ni=n1,n2;
       do mi=m1,m2;

         if (lnstar_T_const) then
           f(l1:l2,mi,ni,iss)=-f(l1:l2,mi,ni,ilnrho)*gamma1/gamma
         else

           if (nxgrid .LE. 1) then

             lnrho=f(l1:l2,mi,ni,ilnrho)
             lnTT=(zz(l1:l2,mi,ni)-R_star)/Lxyz(3) &
                 *(log(T_disk)-log(T_star))+log(T_star)

             call eoscalc(4,lnrho,lnTT,ss=ss)
               f(l1:l2,mi,ni,iss)=ss

           else

            lnrho=f(l1:l2,mi,ni,ilnrho)

            lnTT=(zz(l1:l2,mi,ni)-R_star)/Lxyz(3) &
                *(log(T_disk)-log(T_star))+log(T_star)

        !  if (hot_star) then
	  
	!  for li=l1,l2 do
	!    lnTT(li-3)=(xx(li,mi,ni)-0)/Lxyz(1) &
	!       *(log(T_disk*0.7)-lnTT(li-3))+lnTT(li-3)
        !  enddo			  
	  
	!  endif

        if (hot_star .and. new_T_profile) then
    
         
              TT_0(1)=exp(lnTT(l1))

                do li=l1+1,l2
                dT_dx_i1=-M_star/(zz(li-1,mi,ni))**3*xx(li-1,mi,ni) &
                *3./16./sigmaSB*c_light*exp(f(li-1,mi,ni,ilnrho)) &
                /TT_0(li-l1)**3

                TT_0(li-l1+1)=TT_0(li-l1)+dT_dx_i1*dx           

                lnTT(li-l1+1)=log(TT_0(li-l1+1))

   
                enddo
               
        endif
            

            call eoscalc(4,lnrho,lnTT,ss=ss)

            f(l1:l2,mi,ni,iss)=ss

             !  f(l1,mi,ni,iss)=-f(l1,mi,ni,ilnrho)*gamma1/gamma

           endif
         endif

       end do
      end do


      endsubroutine entropy_init
!*********************************************************************
!***********************************************************************
      subroutine velocity_init(f,zz)
      !Natalia
      !Initialization of velocity in a case of the step-like distribution

      use Cdata

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: zz
      integer :: decel_zone
      real ::   ll
      integer :: L_disk_point=46

      decel_zone=ac_dc_size+4
      ll=L_disk_point*dz

     !ll=Lxyz(3)-L_disk
     if (nzgrid .GT. 1) then

      f(:,:,:,iux)=uu_init
      f(:,:,:,iuz)=uu_init
      
      if (hot_star) then

         f(:,:,:,iuy)=sqrt(M_star/xyz0(3))
	 

      else
       if (ldecelerat_zone .AND. decel_zone .LT. nzgrid) then

        f(:,:,L_disk_point+4:mz,iuy)= &
          sqrt(M_star/zz(:,:,L_disk_point+4:mz))

        f(:,:,decel_zone+1:L_disk_point+3,iuy)= &
           (zz(:,:,decel_zone+1:L_disk_point+3)-R_star-(decel_zone-4)*dz) &
           /(ll-(decel_zone-4)*dz)*sqrt(M_star/(ll+R_star))
        f(:,:,1:decel_zone,iuy)=0.

       else
         f(:,:,L_disk_point+4:mz,iuy)= &
           sqrt(M_star/zz(:,:,L_disk_point+3+1:mz))
         f(:,:,1:L_disk_point+3,iuy)= &
           (zz(:,:,1:L_disk_point+3)-R_star)/ll*sqrt(M_star/(ll+R_star))
       endif

     endif
     else
      f(:,:,:,iux)=uu_init
      f(:,:,:,iuz)=uu_init
      f(:,:,:,iuy)=sqrt(M_star/xyz0(3))
     endif


    endsubroutine  velocity_init
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
          do i=1,nghost; f(l1-i,:,:,j)=-f(l1+i,:,:,j); enddo
              f(l1,:,:,j) = 0.
        !  do i=1,nghost; f(l1-i,:,:,j)=2*f(l1,:,:,j)+sgn*f(l1+i,:,:,j); enddo

        else
          do i=1,nghost; f(l1-i,:,:,j)= f(l1+i,:,:,j); enddo

          ! if (j==5) then
        !     f(l1,:,n1,j)=(f(l1+1,:,n1,j)+f(l1-1,:,n1,j))/2
        !   endif

       endif

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


           f(l2+1,:,:,j)=0.25*(  9*f(l2,:,:,j)- 3*f(l2-1,:,:,j)- 5*f(l2-2,:,:,j)+ 3*f(l2-3,:,:,j))
            f(l2+2,:,:,j)=0.05*( 81*f(l2,:,:,j)-43*f(l2-1,:,:,j)-57*f(l2-2,:,:,j)+39*f(l2-3,:,:,j))
            f(l2+3,:,:,j)=0.05*(127*f(l2,:,:,j)-81*f(l2-1,:,:,j)-99*f(l2-2,:,:,j)+73*f(l2-3,:,:,j))

        endif
      else
        print*, "bc_BL_x: ", bc%location, " should be `top(", &
                        iBC_X_TOP,")' or `bot(",iBC_X_BOT,")'"
      endif
!
    endsubroutine bc_BL_x
 !********************************************************************
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
      real :: value1,value2, dT_dx_i1
      type (boundary_condition) :: bc
      real, dimension(nx) :: lnrho,lnTT,ss,TT_0
      integer :: sgn,i,j,  n1p4,n2m4, i_tmp,li

    j=bc%ivar


      if (j == 4 .or. j==5) then
        value1=log(bc%value1)
        value2=log(bc%value2)
      else
        value1=bc%value1
        value2=bc%value2
      endif

      if (bc%location==iBC_Z_BOT) then
      ! bottom boundary
    !    if (lextrapolate_bot_density .and. j>=4) then
!          n1p4=n1+4

!          f(:,:,n1-1,j)=0.2   *(  9*f(:,:,n1,j)  -  4*f(:,:,n1+2,j)- 3*f(:,:,n1+3,j)+ 3*f(:,:,n1p4,j))
!          f(:,:,n1-2,j)=0.2   *( 15*f(:,:,n1,j)- 2*f(:,:,n1+1,j)-  9*f(:,:,n1+2,j)- 6*f(:,:,n1+3,j)+ 7*f(:,:,n1p4,j))
!          f(:,:,n1-3,j)=1./35.*(157*f(:,:,n1,j)-33*f(:,:,n1+1,j)-108*f(:,:,n1+2,j)-68*f(:,:,n1+3,j)+87*f(:,:,n1p4,j))
!        else



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
           ! f(l1:l2,m1,n1,iss)=log(T_star)/gamma

          else
            if ( j==4 ) then
            else
	    if (j==1 .and. hot_star) then
	    else
	      if (j==2) then
               f(:,:,n1,j)=sqrt(M_star/R_star)*star_rot
	      else
              f(:,:,n1,j)=value1
	      endif
	    endif  
           endif
          endif


          do i=1,nghost; f(:,:,n1-i,j)=2*f(:,:,n1,j)+sgn*f(:,:,n1+i,j); enddo
!endif
      elseif (bc%location==iBC_Z_TOP) then
      ! top boundary

        if (ltop_velocity_kep .and. j==2) then
          f(:,:,n2,j)=sqrt(M_star/(R_star+Lxyz(3)))
        else
         ! if (nxgrid <= 1) then
            if (j==5) then
              lnrho=f(l1:l2,m2,n2,ilnrho)
              if (T_disk.EQ.0) then
              lnTT=log(cs0**2/(gamma1))
              else

               if (hot_star .and. new_T_profile) then
                                  
                TT_0(1)=T_disk

                do li=l1+1,l2
                dT_dx_i1=-M_star/(R_star+Lxyz(3))**3*(li-l1)*dx &
                *3./16./sigmaSB*c_light*exp(f(li-1,m2,n2,4)) &
                /TT_0(li-l1)**3

                TT_0(li-l1+1)=TT_0(li-l1)+dT_dx_i1*dx

                enddo
                lnTT=log(TT_0)    
               !  lnTT=log(T_disk)
               else   
                lnTT=log(T_disk)
               endif 
              endif
              call eoscalc(4,lnrho,lnTT,ss=ss)
              f(l1:l2,m2,n2,iss)=ss
            else
              if (j==4) then
              else
                f(:,:,n2,j)=value1
              endif
            endif
         ! else

         ! endif
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

