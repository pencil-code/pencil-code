! $Id$
!
!  This module includes physics related to planet atmospheres.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
! variables global to this module
!
  real, dimension(my,mz) :: mu_ss=0.
  real, dimension(my) :: lat
  real, dimension(mz) :: lon
  real :: pp2Pa=1., TT2K=1., tt2s=1.  !  convert pc units to [Pa] and [K] and [s]
!
! variables in the reference profile
!
  real :: dlogp_ref, logp_ref_min, logp_ref_max
  real, dimension(:), allocatable :: p_temp_ref,tau_rad,temp_ref,logp_ref,dTeq,Teq_night
  integer :: nref
!
! Init parameters
!
  real :: lon_ss=0., lat_ss=0.  ! unit: [degree]
  real :: dTeqbot=0., dTeqtop=100.  ! unit: [K]
  real :: peqtop=1.d2, peqbot=1.d6  ! unit: [Pa]
  real :: tauradtop=1.d4, tauradbot=1.d7  ! unit: [s]
  real :: pradtop=1.d3, pradbot=1.d6 ! unit:[ Pa]
!
! Run parameters
!
  real :: tau_slow_heating=-1.
  real :: Bext_dipole=0.
!
!
!
  namelist /special_init_pars/ &
      lon_ss,lat_ss,peqtop,peqbot,tauradtop,tauradbot,&
      pradtop,pradbot,dTeqbot,dTeqtop,pp2Pa,TT2K,tt2s
!
  namelist /special_run_pars/ &
      tau_slow_heating,Bext_dipole
!
!
! Declare index of new variables in f array (if any).
!
!!   integer :: ispecial=0
!!   integer :: ispecaux=0
!
!! Diagnostic variables (needs to be consistent with reset list below).
!
!!   integer :: idiag_POSSIBLEDIAGNOSTIC=0
!
  contains
!****************************************************************************
    subroutine initialize_special(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (.not.ltemperature_nolog) call fatal_error('initialize_special', &
          'special/planet_atmosphere is formulated in TT only')
!
!  convert y and z to latitude and longitude in rad
!
      lat=0.5*pi-y
      lon=z-pi
!
!  unit conversion
!
      call prepare_unit_conversion
!
!  read in the reference P-T profile, and calculate Teq and tau_rad etc.
!
      call prepare_Tref_and_tau
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
    lpenc_requested(i_rho)=.true.
    lpenc_requested(i_pp)=.true.
    lpenc_requested(i_uu)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine special_calc_energy(f,df,p)
!
!  08-sep-23/hongzhe: specific things for planet atmospheres.
!                     Reference: Rogers & Komacek (2014)
!  27-sep-23/hongzhe: outsourced from temperature_idealgas.f90
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx) :: Teq_x,tau_rad_x,log10pp
      real ,dimension(:), allocatable :: Teq_local
      real :: f_slow_heating
      integer :: ix,index
!
!  The local equilibrium T; still in pressure coordinate
!
      allocate(Teq_local(nref))
      Teq_local = Teq_night + dTeq*max(0.,mu_ss(m,n))
!
!  interpolation for Teq_local at height coordinate (could have arbitrary pressure)
!
      log10pp = log10(p%pp*pp2Pa)
      do ix=1,nx
        ! index of the logp_ref that is just smaller than log10(pressure)
        index = 1+floor((log10pp(ix)-logp_ref_min)/dlogp_ref)
        if (index>=nref) then
          Teq_x(ix) = Teq_local(nref)
          tau_rad_x(ix) = tau_rad(nref)
        elseif (index<= 1) then
          Teq_x(ix) = Teq_local(1)
          tau_rad_x(ix) = tau_rad(1)
        else
          ! interpolate T in log10(p) space
          Teq_x(ix) = Teq_local(index)+(Teq_local(index+1)-Teq_local(index))*   &
                          (log10pp(ix)-logp_ref(index))/   &
                          (logp_ref(index+1)-logp_ref(index))   ! unit of K
          ! interpolate log10(tau_rad) in log10 p space
          tau_rad_x(ix) = log10(tau_rad(index))+  &
                              (log10(tau_rad(index+1))-log10(tau_rad(index)))*   &
                              (log10pp(ix)-logp_ref(index))/   &
                              (logp_ref(index+1)-logp_ref(index))
          tau_rad_x(ix) = 10**tau_rad_x(ix)  ! unit of s
        endif
      enddo
!
      if (tau_slow_heating>0) then
        f_slow_heating = min(1.d0,t/tau_slow_heating)
      else
        f_slow_heating = 1.
      endif
  !
      df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) - f_slow_heating*(p%TT-Teq_x/TT2K)/(tau_rad_x/tt2s)
!
      deallocate(Teq_local)
!
    endsubroutine special_calc_energy
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!  09-oct-23/hongzhe: add a background dipole field
!
      use Sub, only: cross_mn
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx,3) :: uxb_ext,Bdipole
!
!  the r,theta,phi components of the dipole filed
!
      Bdipole(:,1) = Bext_dipole / (x(l1:l2)**3.)
      Bdipole(:,2) = Bext_dipole / (x(l1:l2)**3.)
      Bdipole(:,3) = Bext_dipole / (x(l1:l2)**3.)
!
      call cross_mn(p%uu,Bdipole,uxb_ext)
      df(l1:l2,m,n,iax:iaz) = df(l1:l2,m,n,iax:iaz) + uxb_ext
!
    endsubroutine special_calc_magnetic
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
!  compute cos(angle between the substellar point)
!
      call get_mu_ss(mu_ss,lon_ss,lat_ss)
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine prepare_unit_conversion
!
!  Constants that convert code units to SI.
!
!  28-sep-23/hongzhe: coded
!
 !     use EquationOfState, only: getmu
!
!      real :: mu  !  mean molecular weight
!
!      call getmu(mu_tmp=mu)
!
      if (unit_system=='SI') then
        ! %HZ: need to code correctly!!
        ! m_u is atomic mass unit, not mean molecular weight
        !pp2Pa=k_B_cgs/m_u_cgs*1.0e-4/mu*unit_density*unit_temperature
        !TT2K=1.*unit_temperature
      else
        call fatal_error('prepare_unit_conversion','please use SI system')
      endif
!
    endsubroutine  prepare_unit_conversion
!***********************************************************************
    subroutine prepare_Tref_and_tau
!
!   Read the reference profile and calculated tau_rad.
!   All quantities in this subroutine are in SI units.
!
!   28-sep-23/xianyu,hongzhe: coded
!
    real :: alpha
    integer :: i
    logical :: lTref_file_exists
!
!  read in Tref, in physical unit.
!
!!HZ: need use a less specific file name
      inquire(FILE='iro-teq-tint100K-regrid-Pa.txt', EXIST=lTref_file_exists)
      if (.not.lTref_file_exists) call fatal_error('initialize_special', &
          'Must provide a Tref file')
!
      open(1,file='iro-teq-tint100K-regrid-Pa.txt')
      read(1,*)
      read(1,*) dlogp_ref, logp_ref_min, logp_ref_max  ! in log10(bar)
      read(1,*) nref
      if(allocated(logp_ref)) deallocate(logp_ref)
      if(allocated(temp_ref)) deallocate(temp_ref)
      allocate(logp_ref(nref),temp_ref(nref))
      read(1,*) logp_ref,temp_ref
      if (lroot) then
        print*, 'nref=',nref
        print *,'Here is the baseline radiative equil T profile'
        print *,'used in the Newtonian relaxation:'
        print *,'p [Pa] and T_eq[K]:'
        do i=1,nref
          print*, 'logp_ref,temp_ref=',logp_ref(i),temp_ref(i)
        enddo
      endif
      close(1)
!
!  convert units, and calculate tau_rad, dTeq, and Teq_night
!
      if(allocated(p_temp_ref)) deallocate(p_temp_ref)
      if(allocated(tau_rad))    deallocate(tau_rad)
      if(allocated(dTeq))       deallocate(dTeq)
      if(allocated(Teq_night))  deallocate(Teq_night)
      allocate( p_temp_ref(nref), tau_rad(nref), dTeq(nref),Teq_night(nref) )
      p_temp_ref = 10.**logp_ref
  !
      where (p_temp_ref>=peqtop .and. p_temp_ref<=peqbot)
        dTeq = dTeqbot + (dTeqtop - dTeqbot) * log(p_temp_ref/peqbot)/log(peqtop/peqbot)
      elsewhere (p_temp_ref < peqtop)
        dTeq = dTeqtop
      elsewhere
        dTeq = dTeqbot
      endwhere
!
      alpha=log(tauradtop/tauradbot)/log(pradtop/pradbot)
      where (p_temp_ref>=pradtop .and. p_temp_ref<=pradbot)
        tau_rad = tauradbot * (p_temp_ref/pradbot)**alpha
      elsewhere (p_temp_ref < pradtop)
        tau_rad = tauradtop
      elsewhere
        tau_rad = tauradbot
      endwhere
!
      Teq_night = temp_ref - 0.5*dTeq
!
      deallocate(temp_ref,p_temp_ref)
!
    endsubroutine  prepare_Tref_and_tau
!***********************************************************************
    subroutine get_mu_ss(mu_ss,lonss,latss)
!
!  Given the substellar longitude lonss and latitude latss (both in
!  degree), compute the cosine of the angle of the sun from overhead.
!  This will be used in future RT calculation.
!  It also tells you whether you're on the dayside, since
!  mu_ss>0 for dayside, mu_ss=0 on the terminator, and mu_ss<0 on the nightside.
!
!  28-sep-23/xianyu,hongzhe: coded
!
      real, dimension(my,mz), intent(out) :: mu_ss
      real, intent (in) :: lonss, latss
!
      real, PARAMETER :: deg2rad=pi/180.
!
      mu_ss = spread(cos(lon),1,my) * spread(cos(lat),2,mz) &
             * cos(lonss*deg2rad) * cos(latss*deg2rad) + &
              spread(sin(lon),1,my)*spread(cos(lat),2,mz) &
             * sin(lonss*deg2rad) * cos(latss*deg2rad) + &
              spread(sin(lat),2,mz)*sin(latss*deg2rad)
!
    endsubroutine  get_mu_ss
!***********************************************************************
!********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING        *************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copied dummy routines from nospecial.f90 for any              **
!**  routine not implemented in this file.                         **
!**                                                                **
    include '../special_dummies.inc'
!*********************************************************************** 
endmodule Special
