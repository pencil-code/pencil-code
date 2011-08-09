!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 2
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
!!    SPECIAL=special/nstar
!! Where nstar it replaced by the filename of your new module
! upto and not including the .f90
!
!--------------------------------------------------------------------
!
module Special
!
  use Cdata
  use Cparam
  use Messages, only: svn_id, fatal_error
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../special.h'
!
  real :: mtemp=21,mrho=41,mzeta=61
  !real, dimension (mtemp,mrho,mzeta) :: eta_table
  !real, dimension (mtemp) :: TT_table
  !real, dimension (mrho)  :: rho_table
  !real, dimension (mzeta) :: zeta_table
  real, dimension (21,41,61) :: eta_table
  real, dimension (21) :: TT_table
  real, dimension (41)  :: rho_table
  real, dimension (61) :: zeta_table
  real, dimension (nx) :: deltaz
  real :: dummy,minval_zeta_table=1.1e-24 
  real :: unit_density_cgs=1.,unit_length_cgs=1.
  real :: unit_time_cgs=1.,unit_velocity_cgs
  real :: zeta_radionuclides_ref=4d-17,zeta_radionuclides
  real, dimension(nx,nz,0:ncpus-1) :: density_column
  logical :: lzeta_cosmicray,lzeta_xray,lzeta_nuclides
  real :: dust_to_gas=1d-5
!
  namelist /special_init_pars/ dummy
!
  namelist /special_run_pars/ unit_density_cgs,&       
       unit_length_cgs,unit_time_cgs,&
       lzeta_cosmicray,lzeta_xray,lzeta_nuclides, &
       minval_zeta_table,dust_to_gas
!
! Declare index of new variables in f array (if any).
!
!!   integer :: ispecial=0
!!   integer :: ispecaux=0
  integer :: ieta=0
  integer :: izeta=0
!
!! Diagnostic variables (needs to be consistent with reset list below).
!
!!   integer :: idiag_POSSIBLEDIAGNOSTIC=0
!
  contains
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
           "$Id$")
!
!!      call farray_register_pde('special',ispecial)
      call farray_register_auxiliary('eta',ieta)
      call farray_register_auxiliary('zeta',izeta)
!!      call farray_register_auxiliary('specaux',ispecaux,communicated=.true.)
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      use General, only: itoa
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting,exist
      real :: TT,rho,zeta,eta,teqm
      integer :: mtable,itable,irho,itemp,izet
      integer :: iitemp,iirho,iizet,stat
      real :: unit_eta_cgs,unit_eta1_cgs
      character (len=fnlen)  :: tablefile
      character (len=intlen) :: sdust 
!
!  Work out the units
!
      unit_velocity_cgs=unit_length_cgs/unit_time_cgs
      unit_eta_cgs=unit_length_cgs**2/unit_time_cgs
      unit_eta1_cgs=1./unit_eta_cgs
!
      if (lroot) then 
        print*,'unit_velocity=',unit_velocity_cgs,' cm/s'
        print*,'unit_time=',unit_time_cgs,' s'
        print*,'unit_resistivity=',unit_eta_cgs,' cm2/s'
      endif
!
!  Scale the ionization due to radionuclides with the dust to gas ratio
!
      zeta_radionuclides=zeta_radionuclides_ref*dust_to_gas
!
!  Read the lookup table into a local 3D array. All processors should have it 
!
      sdust=itoa(nint(-log10(dust_to_gas)))
      tablefile=trim('./table/table.eta.a0p1um.d2gmr1em'//trim(sdust)//&
           '.mgeps1em4.within2.dat')
      print*,'reading table ',tablefile
!     
      inquire(file=tablefile,exist=exist)
      if (exist) then
        open(19,file=tablefile)
      else
        call fatal_error('resistivity','no input file')
      endif

      mtable=52521
!
      do itable=1,mtable
        read(19,*,iostat=stat) itemp,irho,izet,TT,rho,zeta,eta,teqm
        iitemp=itemp+1;iirho=irho+1;iizet=izet+1
!
        TT_table(iitemp)=TT
        rho_table(iirho)=rho
        zeta_table(iizet)=zeta
        eta_table(iitemp,iirho,iizet)=eta*unit_eta1_cgs
      enddo
!
      close(19)
!
      deltaz=abs(x(l1:l2)*cos(y(2))-x(l1:l2)*cos(y(1)))
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
!      select case (initspecial)
!        case ('nothing'); if (lroot) print*,'init_special: nothing'
!        case ('zero', '0'); f(:,:,:,iSPECIAL_VARIABLE_INDEX) = 0.
!        case default
!          call fatal_error("init_special: No such value for initspecial:" &
!              ,trim(initspecial))
!      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
      lpenc_requested(i_jj)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Special pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-nov-04/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
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
!  Several precalculated Pencils of information are passed for
!  efficiency.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!!      if (headtt) call identify_bcs('special',ispecial)
!
!!
!! SAMPLE DIAGNOSTIC IMPLEMENTATION
!!
!!      if (ldiagnos) then
!!        if (idiag_SPECIAL_DIAGNOSTIC/=0) then
!!          call sum_mn_name(MATHEMATICAL EXPRESSION,idiag_SPECIAL_DIAGNOSTIC)
!!! see also integrate_mn_name
!!        endif
!!      endif
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=special_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_init_pars,ERR=99)
      endif
!
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
!
      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
!!      integer :: iname
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
!        i_SPECIAL_DIAGNOSTIC=0
      endif
!
!      do iname=1,nname
!        call parse_name(iname,cname(iname),cform(iname),&
!            'NAMEOFSPECIALDIAGNOSTIC',i_SPECIAL_DIAGNOSTIC)
!      enddo
!
!  write column where which magnetic variable is stored
      if (lwr) then
        write(3,*) 'ieta=' ,ieta
        write(3,*) 'izeta=',izeta
      endif

    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of Special variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine calc_lspecial_pars(f)
!
!  Dummy routine.
!
!  15-jan-08/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lspecial_pars
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  momentum equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  continuity equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  entropy equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  induction equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension(nx) :: eta
      integer :: j,ju
!
!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!
      eta = f(l1:l2,m,n,ieta)
!
      do j=1,3
        ju=iaa+j-1
        df(l1:l2,m,n,ju) = df(l1:l2,m,n,ju) - eta*p%jj(:,j)
      enddo
!
!  Contribution to the timestep
!
      if (lfirst.and.ldt) diffus_eta=diffus_eta+eta
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_magnetic
!***********************************************************************
    subroutine special_before_boundary(f)
!
      use FArrayManager, only: farray_use_global
!
!  Possibility to modify the f array before the boundaries are
!  communicated.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: TT,rho,zeta
      integer :: it1,it2,ir1,ir2,iz1,iz2
      integer :: i,m,n
!
      real :: deta,eta_zeta1,eta_zeta2
      real :: eta_rho1_zeta1,eta_rho2_zeta1
      real :: eta_rho1_zeta2,eta_rho2_zeta2
      real :: dtemp,delta_temp
      real :: drho,delta_rho
      real :: dzeta,delta_zeta,mu_gas
      real :: gamma_loc,cp_loc,Rgas_cgs,cs2
      integer, pointer :: iglobal_cs2
      integer :: ics2
!
!  Pre-calculate the column densities for the parallel case
!
      if (nprocy/=1) call calc_column_densities(f)
!
      gamma_loc=1.4
      Rgas_cgs=8.3d7
      mu_gas=2.3
      cp_loc=gamma_loc*Rgas_cgs/(mu_gas*(gamma_loc-1))
!
      nullify(iglobal_cs2)
      call farray_use_global('cs2',iglobal_cs2)
      ics2=iglobal_cs2 
!
      do n=n1,n2; do m=m1,m2; do i=l1,l2
        rho=f(i,m,n,ilnrho)*unit_density_cgs
        cs2=f(i,m,n,ics2)*unit_velocity_cgs**2
        TT =cs2/(cp_loc*(gamma_loc-1))
        call calc_ionization_rate(f,zeta,i,m,n)
        f(i,m,n,izeta)=zeta
!
! Get the closest two values of temperature, density, and ionization rate
!
        call get_closest_value(TT  ,TT_table  ,it1,it2)
        call get_closest_value(rho ,rho_table ,ir1,ir2)
        call get_closest_value(zeta,zeta_table,iz1,iz2)
!
        if (ip <= 7) then 
          print*,'Temperature closest cells    ',it1,it2,TT  
          print*,'Density closest cells        ',ir1,ir2,rho
          print*,'Ionization rate closest cells',iz1,iz2,zeta
        endif
! 
! Do the trilinear interpolation
!
! First with temperature
!        
        deta=eta_table(it2,ir1,iz1)-eta_table(it1,ir1,iz1)
        dtemp=TT_table(it2)-TT_table(it1)
        if (it2 == it1) dtemp = 1.
        delta_temp = TT-TT_table(it1)
        eta_rho1_zeta1=eta_table(it1,ir1,iz1)  + deta/dtemp * delta_temp
!
        deta=eta_table(it2,ir2,iz1)-eta_table(it1,ir2,iz1)
        eta_rho2_zeta1=eta_table(it1,ir2,iz1) + deta/dtemp * delta_temp
!
! Interpolate in density. Get eta for lower face of the cube
!
        deta=eta_rho2_zeta1-eta_rho1_zeta1
        drho=rho_table(ir2)-rho_table(ir1)
        if (ir2 == ir1) drho = 1.
        delta_rho = rho-rho_table(ir1)
        eta_zeta1=eta_rho1_zeta1 + deta/drho * delta_rho
!
! Now do the upper face of the cube
!
        deta=eta_table(it2,ir1,iz2)-eta_table(it1,ir1,iz2)
        eta_rho1_zeta2=eta_table(it1,ir1,iz2) + deta/dtemp * delta_temp
        deta=eta_table(it2,ir2,iz2)-eta_table(it1,ir2,iz2)
        eta_rho2_zeta2=eta_table(it1,ir2,iz2) + deta/dtemp * delta_temp
!
        deta=eta_rho2_zeta2-eta_rho1_zeta2
        eta_zeta2=eta_rho1_zeta2 + deta/drho * delta_rho
!
! Now interpolate in zeta
!
        deta=eta_zeta2-eta_zeta1
        dzeta=zeta_table(iz2)-zeta_table(iz1)
        if (iz1 == iz2) dzeta=1.
        delta_zeta = zeta - zeta_table(iz1)
!
        f(i,m,n,ieta) = eta_zeta1 + deta/dzeta*delta_zeta
!
      enddo;enddo;enddo
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine get_closest_value(array,array_table,itmp1,itmp2)

      real :: array,tmp
      real, dimension(:), intent(in) :: array_table
      integer, intent(out) :: itmp1,itmp2
      integer :: i,itmp,msize

      tmp=minval(abs(array-array_table))
      itmp=0
      msize=size(array_table,1)

      do i=1,msize
        if (tmp .eq. abs(array-array_table(i))) then 
          itmp=i
        endif
      enddo
      if (array .gt. array_table(itmp)) then
        itmp1=itmp
        itmp2=itmp+1
      else
        itmp2=itmp
        itmp1=itmp-1
      endif
!
      if (itmp .eq. 0) then 
        itmp1=1
        itmp2=1
      endif
      if (itmp .eq. msize) then
        itmp1=msize
        itmp2=msize
      endif
!
    endsubroutine get_closest_value
!***********************************************************************
    subroutine calc_column_densities(f)
!
      use Mpicomm,only: mpibcast_real,mpisend_real,mpirecv_real
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i,n,j,ii,nn
      real :: step
      real, dimension(ny) :: density
      real, dimension(nx,nz) :: density_column_local
!
      do i=l1,l2
        ii=i-l1+1
        step=deltaz(ii)*unit_length_cgs
        do n=n1,n2
          nn=n-n1+1
          density=f(i,m1:m2,n,ilnrho)*unit_density_cgs
          density_column_local(ii,nn)=sum(density)*step
        enddo
      enddo
!
      if (.not.lroot) then 
        call mpisend_real(density_column_local,(/nx,nz/),0,111)
      else
        density_column(:,:,0)=density_column_local
        do j=1,ncpus-1 
          call mpirecv_real(density_column_local,(/nx,nz/),j,111)
          density_column(:,:,j)=density_column_local
        enddo
      endif
      call mpibcast_real(density_column,(/nx,nz,ncpus/))
!
    endsubroutine calc_column_densities
!***********************************************************************
    subroutine calc_ionization_rate(f,zeta,ilp,mlp,nlp)
!
!  This only works for serial runs right now
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, intent(out) :: zeta
      real :: sigma_above,sigma_below
      real :: sigma_above_loc,sigma_below_loc
      integer, intent(in) :: ilp,mlp,nlp
      integer :: m,mm,ii
      real :: cray_ionization,xray_ionization
      real :: cray_depth,xray_depth
      real :: zeta_xray,zeta_cray
      real :: step,tab,tbl,gamcr
      real, dimension(ny) :: density
      integer :: iserial,ipy_loc
!
!  Constants
!
      cray_ionization=5.0d-18 !s-1, half of 1d-17 
      xray_ionization=5.2d-15 !s-1, twice 2.6d-15, which is for H2 molecule
!                             !     X-ray luminosity in young stars
      gamcr=3./4 !Umebayashi fit
!
!  Penetration depth
!
      cray_depth=96. !g/cm2
      xray_depth=8.
!     
!  Cosmic ray ionization
!
      ii=ilp-l1+1
      mm=mlp-m1+1
      step=deltaz(ii)*unit_length_cgs
      density=f(ilp,m1:m2,nlp,ilnrho)*unit_density_cgs
!
! Get stuff above in the local processor
!
      sigma_above_loc=0.
      do m=mm,ny
        sigma_above_loc=sigma_above_loc+density(m)*step
      enddo
!
! Get stuff from the processors above
!
      sigma_above=sigma_above_loc
      if (lmpicomm) then
        if (ipy/=nprocy-1) then 
          do ipy_loc=ipy+1,nprocy-1
            iserial=ipx+nprocx*ipy_loc+nprocx*nprocy*ipz
            sigma_above=sigma_above+density_column(ii,mm,iserial)
          enddo
        endif
      endif
!
! Get stuff below in the local processor
!
      sigma_below_loc=0.
      do m=1,mm
        sigma_below_loc=sigma_below_loc+density(m)*step
      enddo
!
! Get stuff from the processors below
!
      sigma_below=sigma_below_loc
      if (lmpicomm) then 
        if (ipy/=0) then 
          do ipy_loc=0,ipy-1
            iserial=ipx+nprocx*ipy_loc+nprocx*nprocy*ipz
            sigma_below=sigma_below+density_column(ii,mm,iserial)
          enddo
        endif
      endif
!
      tab=sigma_above/cray_depth
      tbl=sigma_below/cray_depth
!
      zeta_cray=cray_ionization*exp(-tab)*(1+tab**gamcr)**(-1./gamcr)
      zeta_cray=zeta_cray+&
                cray_ionization*exp(-tbl)*(1+tbl**gamcr)**(-1./gamcr)
!
      zeta_xray=xray_ionization*r1_mn(ilp-l1+1)**2*&
           (exp(-sigma_above/xray_depth) + exp(-sigma_below/xray_depth))
!
      zeta=0.
!      if (ilp==l1) print*,x(ilp),x(ilp)*cos(y(mlp)),zeta_xray
      if (lzeta_cosmicray) zeta=zeta+zeta_cray
      if (lzeta_xray)      zeta=zeta+zeta_xray
      if (lzeta_nuclides)  zeta=zeta+zeta_radionuclides
!
      if (zeta.lt.minval_zeta_table) zeta=minval_zeta_table 
!
    endsubroutine calc_ionization_rate
!***********************************************************************
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
