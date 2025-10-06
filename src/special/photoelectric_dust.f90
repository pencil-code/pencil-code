! $Id: baroclinic_run.f90 19193 2012-06-30 12:55:46Z wdobler $

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
!
!   Equation of state from Lyra & Kuchner (2013), due to photoelectric
!   heating of the gas by the dust. A polytropic term is added for
!   completeness and to provide support against gravity in stratified
!   simulations. The pressure gradient itself is
!
!    gradP = cs20/(gamma*rho0)*(p%rho*p%grhop(:,j) + p%rhop*p%grho(:,j))
!
!   17-may-12/wlad: coded
!   20-dec-12/wlad: streamlined, and moved here to special
!
module Special
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../special.h'
!
  real :: mu=1.0, kappa=0.0, factor_localiso=0.0, factor_photoelectric=1.0
  real, dimension(ndustspec) :: const_pr=0.0
  logical :: ldust_pressureforce=.true.
  logical :: lradiation_PRdrag=.false.
  logical :: ldivrhop_by_vol=.false.
!
  namelist /special_init_pars/ mu, kappa, factor_localiso, factor_photoelectric, &
       ldust_pressureforce, lradiation_PRdrag, const_pr, ldivrhop_by_vol
!
  namelist /special_run_pars/ mu, kappa, factor_localiso, factor_photoelectric, &
       ldust_pressureforce, lradiation_PRdrag, const_pr, ldivrhop_by_vol
!
!  integer, parameter :: nmode_max = 50
!  real, dimension(nmode_max) :: gauss_ampl, rcenter, phicenter
!
  type InternalPencils
    real, dimension(nx,3) :: fpres
    real, dimension(nx,3) :: fpres_photoelectric
    real, dimension(nx,3) :: fpres_polytropic
    real, dimension(nx,3) :: fpres_localisotropic
  endtype InternalPencils
  type (InternalPencils) :: q
  !$omp threadprivate(q)
!
  real :: const1,const2,const3
  integer :: idiag_photom=0,idiag_photomax=0
  integer :: idiag_photomin=0,idiag_polym=0
  integer :: idiag_polymax=0,idiag_polymin=0
  integer :: idiag_dtcp=0
!
  logical, pointer :: lpressuregradient_gas
  real :: gamma, gamma1
!
  contains
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables 
!  which should be know to be able to evaluate
!
!  14-jul-09/wlad: coded
!
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  14-jul-09/wlad: coded
!
      use Cdata 
      use EquationOfState, only: cs20, get_gamma_etc
      use Mpicomm
      use SharedVariables
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (lentropy) call fatal_error('initialize_special','This code should be used with noentropy')
!
      call get_gamma_etc(gamma)
      gamma1=1./gamma
!
!  Initializing the pressures. Seeting mu=1 makes the polytrope a linear barotrope.
!
      const1=kappa*mu
      const2=factor_photoelectric * cs20 * gamma1
      const3=factor_localiso
!
      call get_shared_variable('lpressuregradient_gas',lpressuregradient_gas)
      if (lrun) then 
        if (lpressuregradient_gas) &
          call fatal_error('initialize_special','switch lpressuregradient_gas=F in hydro_run_pars') 
      endif  
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
      use Mpicomm
!
      if (const1 /= 0.0) then 
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
!
      if (const2 /= 0.0) then 
        lpenc_requested(i_glnrho)=.true.
        if (lparticles) then 
          lpenc_requested(i_rhop)=.true.
          lpenc_requested(i_grhop)=.true.
        elseif (ldustdensity) then 
          lpenc_requested(i_rhodsum)=.true.
          lpenc_requested(i_glnrhodsum)=.true.
        else          
          call fatal_error("pencil_criteria_special",&
               "the world is flat, and we never got here")
        endif
      endif  
!
      if (const3 /= 0.0) then
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
!
      lpenc_requested(i_cs2)=.true.
!
      if (lradiation_PRdrag) then
        lpenc_requested(i_uud)=.true.
      endif
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate the pressure gradient. For the polytropic, P=S*rho**mu, 
!  so fpres = - 1./rho * grad(P) = - S*mu * rho**(mu-2) * grad(rho) 
!  
!  For the photoelectric pressure, P=cs20/(gamma*/rho0) * rho*rhod
!
      !use Particles_sub, only: find_grid_volume
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nx) :: areas, rhop_tmp
      real, dimension(nx,3) :: grhop_tmp
      type (pencil_case) :: p
      integer :: j, l
!
      if (ldivrhop_by_vol) then
        areas = x(l1:l2)*y(m)/dy_1(m)
        grhop_tmp(:,1) = p%grhop(:,1)/areas
        grhop_tmp(:,2) = p%grhop(:,2)/areas
        grhop_tmp(:,3) = p%grhop(:,3)/areas
        rhop_tmp = p%rhop(:)/areas
      else
        grhop_tmp = p%grhop
        rhop_tmp = p%rhop
      endif
!
      do j=1,3
        if (const1/=0.0) then 
          q%fpres_polytropic(:,j) = -const1 * p%rho**(mu-1) * p%glnrho(:,j)
        else
          q%fpres_polytropic(:,j) = 0.
        endif
!
        if (const2/=0.0) then
          if (lparticles) then
            q%fpres_photoelectric(:,j) = -const2 * (grhop_tmp(:,j) + rhop_tmp*p%glnrho(:,j))
          else
            q%fpres_photoelectric(:,j) = -const2 * p%rhodsum * (p%glnrhodsum(:,j) + p%glnrho(:,j))
          endif   
        else
          q%fpres_photoelectric(:,j) = 0.
        endif
!
        if (const3/=0.0) then
          q%fpres_localisotropic(:,j) = -const3*p%cs2*(p%glnrho(:,j)+p%glnTT(:,j))
        endif
      enddo
!
      q%fpres = q%fpres_localisotropic + q%fpres_photoelectric + q%fpres_polytropic
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine special_before_boundary(f)
!
!  This subroutine calculates the full potential due to the turbulence.
!
!  03-oct-12/wlad: coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_before_boundary
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
    subroutine rprint_special(lreset,lwrite)
!
      use Diagnostics
!
!  reads and registers print parameters relevant to special
!
!   14-jul-09/wlad: coded
!
      use FArrayManager, only: farray_index_append
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Write information to index.pro
!
      if (lreset) then
        idiag_photom=0;idiag_photomax=0;idiag_photomin=0
        idiag_polym=0;idiag_polymax=0;idiag_polymin=0
        idiag_dtcp=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'photom',idiag_photom)
        call parse_name(iname,cname(iname),cform(iname),'photomax',idiag_photomax) 
        call parse_name(iname,cname(iname),cform(iname),'photomin',idiag_photomin) 
        call parse_name(iname,cname(iname),cform(iname),'polym',idiag_polym)
        call parse_name(iname,cname(iname),cform(iname),'polymax',idiag_polymax)
        call parse_name(iname,cname(iname),cform(iname),'polymin',idiag_polymin)
        call parse_name(iname,cname(iname),cform(iname),'dtcp',idiag_dtcp)
      enddo
!
      if (lwr) then
        call farray_index_append('i_photom',idiag_photom)
        call farray_index_append('i_photomax',idiag_photomax)
        call farray_index_append('i_photomin',idiag_photomin)
        call farray_index_append('i_polym',idiag_polym)
        call farray_index_append('i_polymax',idiag_polymax)
        call farray_index_append('i_polymin',idiag_polymin)
        call farray_index_append('i_dtcp',idiag_dtcp)
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   momentum equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
      use Cdata
      use Diagnostics
!      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      integer :: j,ju,k
!
      if (lfirst.and.ldt) advec_cs2 = (const3*p%cs2 + const2*gamma1 + const1) * dxyz_2
      if (headtt.or.ldebug) print*, 'dss_dt: max(advec_cs2) =', maxval(advec_cs2)
!
!  Modified momentum equation
!
      if (ldust_pressureforce) then
        do j=1,3 
          ju=j+iuu-1
          df(l1:l2,m,n,ju) = df(l1:l2,m,n,ju) + q%fpres(:,j)
        enddo
      endif
!
      if (lradiation_PRdrag) then
         do k=1,ndustspec 
            df(l1:l2,m,n,iudx(k)) = df(l1:l2,m,n,iudx(k)) - 2*const_pr(k)*p%uud(:,1,k)
            df(l1:l2,m,n,iudy(k)) = df(l1:l2,m,n,iudy(k)) -   const_pr(k)*p%uud(:,2,k)
         enddo
      endif
!
      if (ldiagnos) then 
        if (idiag_photom/=0)   call sum_mn_name( q%fpres_photoelectric(:,1),idiag_photom)
        if (idiag_photomax/=0) call max_mn_name( q%fpres_photoelectric,idiag_photomax)
        if (idiag_photomin/=0) call max_mn_name(-q%fpres_photoelectric,idiag_photomin,lneg=.true.)
        if (idiag_polym/=0)    call sum_mn_name( q%fpres_polytropic(:,1),idiag_polym)
        if (idiag_polymax/=0)  call max_mn_name( q%fpres_polytropic,idiag_polymax)
        if (idiag_polymin/=0)  call max_mn_name(-q%fpres_polytropic,idiag_polymin,lneg=.true.)
        if (idiag_dtcp/=0)     call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtcp,l_dt=.true.)
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine dspecial_dt_ode

    endsubroutine dspecial_dt_ode
!***********************************************************************
    subroutine pushpars2c

      use Syscalls, only: copy_addr
      use General , only: string_to_enum

      integer, parameter :: n_pars=20
      integer(KIND=ikind8), dimension(n_pars) :: p_par

      call copy_addr(mu,p_par(1))
      call copy_addr(ldust_pressureforce,p_par(2)) ! bool
      call copy_addr(lradiation_prdrag,p_par(3)) ! bool
      call copy_addr(ldivrhop_by_vol,p_par(4)) ! bool
      call copy_addr(const1,p_par(5))
      call copy_addr(const2,p_par(6))
      call copy_addr(const3,p_par(7))
      call copy_addr(gamma1,p_par(8))
      call copy_addr(const_pr,p_par(9)) ! (ndustspec)

    endsubroutine pushpars2c
!***********************************************************************
!
!***********************************************************************
!***********************************************************************
!********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special

