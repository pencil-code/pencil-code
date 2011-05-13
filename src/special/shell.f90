! $Id: nospecial.f90 13879 2010-05-13 15:28:09Z sven.bingert $
!
! Runs a shell model (only GOY implemented) for turbulence
! Copied from Dhruba's code
! Uses a slaved, Adams-Bashforth scheme run in parallel with the PC
! time stepping
!
! Not production level
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

module Special

  use Cdata
  use Cparam
  use Messages, only: svn_id, fatal_error
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'special.h'
!
! constants
!
   complex*16 :: zi=cmplx(0.,1.), zone=cmplx(1.,1.)
!
! run parameters
!
  integer :: nshell=20                                        !number of shells
  real*8  :: lambda=2.                                        !shell separation
  real*8  :: visc=1e-7                                        !viscosity
  real*8  :: f_amp=1.0                                        !forcing amplitude
  real*8  :: deltG=1e-5                                       !internal time step
  real*8  :: k_amp=1.                                         !k normalization
  real*8  :: u_amp=1.                                         !u normalization
  real    :: nstability=2e6                                   !timesteps during start
  complex*16 :: zforce
!
! output parameters
!
  logical :: lstructure=.true.                                !output structure functions
  integer :: nord=1                                           !order of structure function output
!
! real-space parameters
!
  logical :: l_needspecialuu=.false.                          !Need real-space velocity
  logical :: l_quenched=.true.                                !NOT update real-space velocity vectors
!
! runtime variables
!
! aval,bval,cval= k_n*a_n, k_n*b_n/k, k_n*c_n/k^2 from Jensen 1999
!
  real*8, allocatable, dimension(:) :: aval, bval, cval       !evolution equation coefficients
  real*8, allocatable, dimension(:) :: exfac, exrat           !evolution equation variables
  real*8, allocatable, dimension(:) :: umod, eddy_time        !velocity magnitude and turnover times
  real*8, allocatable, dimension(:) :: kval, ksquare          !k vectors
  complex*16, allocatable, dimension(:) :: zu, zgprev         !fourier space velocities
!
! real-space variables
!
  real, allocatable, dimension(:,:) :: kvec, uvec             !to generate real-space velocities
  real, allocatable, dimension(:) :: tnext_vel                !if l_quenched=F, update kvec,uvec
                                                              !at eddy turnover times
  character (len=labellen) :: time_step_type='nothing'        !Shell contribution to time-step
  integer :: minshell=-1, maxshell=-1                         !bounds of physically relevant
                                                              !shell (time-step contraints)
!
  namelist /special_init_pars/ &
      nshell, visc, f_amp, lambda, deltG, k_amp, u_amp, nstability,&
      l_needspecialuu
!
  namelist /special_run_pars/ &
      lstructure, l_quenched, nord,time_step_type, minshell, maxshell
!
!! Diagnostic variables (needs to be consistent with reset list below).
!
   integer :: idiag_uu4=0, idiag_deltM=0, idiag_eddy_time=0
!
  contains
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
      if (lroot) call svn_id( &
           "$Id: nospecial.f90 13879 2010-05-13 15:28:09Z sven.bingert $")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
      Use Mpicomm, only: mpibcast_cmplx_arr
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
      integer :: ns, ierr
!
! allocate arrays (init_special calls initialize_special)
!
      if(.not.allocated(zu)) then
        allocate(zu(nshell), zgprev(nshell), eddy_time(nshell))
        allocate(kval(nshell),umod(nshell),aval(nshell),bval(nshell))
        allocate(cval(nshell),exfac(nshell), exrat(nshell),ksquare(nshell))
        if (l_needspecialuu) allocate(kvec(nshell,3),uvec(nshell,3))
        if (l_needspecialuu .and. .not. l_quenched) allocate(tnext_vel(nshell))
      endif
!
      call keep_compiler_quiet(f)
!
! set up k vector and evolution coefficients
!
      do ns=1,nshell
        kval(ns)=k_amp*(lambda**dfloat(ns))
        ksquare(ns)=kval(ns)*kval(ns)
        exfac(ns) = dexp(-visc*ksquare(ns)*deltG)
        exrat(ns) = (1.0d0 - exfac(ns))/(visc*ksquare(ns))
!
        if (ns .ge. (nshell-1)) then
          aval(ns)=0.
        else
          aval(ns)=kval(ns)
        endif
        if ((ns .eq. 1) .or. (ns .eq. nshell)) then
          bval(ns)=0.
        else
          bval(ns)=-kval(ns-1)/2.
        endif
        if (ns .le. 2) then
          cval(ns)=0.
        else
          cval(ns)=-kval(ns-2)/2.
        endif
      enddo
!      
      zforce= zone*f_amp
!
! If not called by init_special, read in goyvar.dat
! calculate umod, turnover times
!
      if (.not.lstarting) then
        call GOY_read_snapshot(trim(directory_snap)//'/goyvar.dat')
        umod=abs(zu)
        eddy_time=1./kval/umod
      endif
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!
      Use General
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(nshell) :: random_phase_uu_init
      real*8, dimension(nshell) ::randomd
      integer :: ns, count=1
      logical :: count_out
      real :: noneoverthree=-1./3.
!
!  set up random phase
!
      call random_number_wrapper(random_phase_uu_init)
      random_phase_uu_init=random_phase_uu_init*2.*pi
      randomd=dble(random_phase_uu_init)
!
!  set up initial velocity (k^-1/3 powerlaw); only on root
!
      if (lroot) then
        call initialize_special(f,lstarting=.true.)
        do ns=1,nshell
          umod(ns)=u_amp*kval(ns)**noneoverthree
          zu(ns)=dcmplx(umod(ns)*dcos(randomd(ns)), &
              umod(ns)*dsin(randomd(ns)))
          zgprev=dcmplx(0.,0.)
        enddo
        umod=abs(zu)
        call write_structure()
!
! stabilize velocities
!
        do while (count .lt. nstability)
          call GOY_advance_timestep()
          count=count+1
          if (mod(count, 100000)==0) then
            print*, 'stabilizing, count=',count
            umod=abs(zu)
            call write_structure()
          endif
        enddo
      endif
!
! if we need real-space velocities, 
!
      if (l_needspecialuu) then
        if (lroot) call get_vectors(kvec,uvec)
      endif
!
!  write vars
!
      call GOY_write_snapshot(trim(directory_snap)//'/goyvar.dat',ENUM=.false.)
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
!  does not advance the shell model, that is handled
!  in special_after_timestep
!
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: deltm1, minturn, deltmbase
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  calculate timestep for GOY
!
      if (visc/=0. .and. lfirstpoint) then
        if (lroot) then
          select case (time_step_type)
            case ('smallest')
              minturn=minval(eddy_time(3:nshell-2)) !minimum turnover largest k
              deltmbase=min(0.4*minturn,0.08/visc/ksquare(nshell-2))
            case ('largest')
              minturn=1./umod(3)/kval(nshell-2)  !crossing time of smallest eddies
             deltmbase=min(0.4*minturn,0.08/visc/ksquare(nshell-2))
            case ('nothing')
              deltmbase=dtmax
          endselect
        deltm1=1./deltmbase
        endif
        call mpibcast_real(deltm1,1)
        dt1_special=deltm1
      endif
!
      if (ldiagnos) then
!
! output structure functions to file
!
        if (lstructure.and.lfirstpoint.and.lroot) call write_structure()
!
! standard diagnostics
!
        if (lroot) then
          if (idiag_uu4/=0)   call save_name(sngl(umod(4)),idiag_uu4) 
          if (idiag_deltm/=0) call save_name(deltmbase,idiag_deltm)       
          if (idiag_eddy_time/=0) call save_name(sngl(eddy_time(nshell-2)),idiag_eddy_time)   
        endif
      endif
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
!
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
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
!
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
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
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
      use Diagnostics
      use Sub
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!!!
!!!  reset everything in case of reset
!!!  (this needs to be consistent with what is defined above!)
!!!
      if (lreset) then
        idiag_uu4=0
        idiag_deltm=0
        idiag_eddy_time=0
      endif
!!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'uu4',idiag_uu4)
        call parse_name(iname,cname(iname),cform(iname),'deltm',idiag_deltm)
        call parse_name(iname,cname(iname),cform(iname),'eddy_time',idiag_eddy_time)
      enddo
!!
!!!  write column where which magnetic variable is stored
      if (lwr) then
        write(3,*) 'i_uu4=',idiag_uu4
        write(3,*) 'i_deltm=',idiag_deltm
        write(3,*) 'i_eddy_time=',idiag_eddy_time
      endif
!!
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
    endsubroutine special_calc_magnetic
!***********************************************************************
    subroutine special_calc_pscalar(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  passive scalar equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array.
!
!  15-jun-09/anders: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,ilncc) = df(l1:l2,m,n,ilncc) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_pscalar
!***********************************************************************
    subroutine special_calc_particles(fp)
!
!  Called before the loop, in case some particle value is needed
!  for the special density/hydro/magnetic/entropy.
!
!  20-nov-08/wlad: coded
!
      real, dimension (:,:), intent(in) :: fp
!
      call keep_compiler_quiet(fp)
!
    endsubroutine special_calc_particles
!***********************************************************************
    subroutine special_calc_particles_nbody(fsp)
!
!  Called before the loop, in case some massive particles value
!  is needed for the special density/hydro/magnetic/entropy.
!
!  20-nov-08/wlad: coded
!
      real, dimension (:,:), intent(in) :: fsp
!
      call keep_compiler_quiet(fsp)
!
    endsubroutine special_calc_particles_nbody
!***********************************************************************
    subroutine special_before_boundary(f)
!
!  Possibility to modify the f array before the boundaries are
!  communicated.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_boundconds(f,bc)
!
!  Some precalculated pencils of data are passed in for efficiency,
!  others may be calculated directly from the f array.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (boundary_condition) :: bc
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(bc)
!
    endsubroutine special_boundconds
!***********************************************************************
    subroutine special_after_timestep(f,df,dt_)
!
!  Used to advance velocity
!  Advance the shell model, with dt<deltm
!  The physical system is advancing dt_
!
      Use Mpicomm, only: mpibcast_cmplx_arr, mpibcast_real,stop_it
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dt_                  !dt_ passed by timestep
      integer :: count, advances   !number of GOY timesteps to call
      logical, dimension(nshell) :: shell
      logical :: update_vecs
!
      umod=abs(zu)
      eddy_time(:)=1./kval(:)/umod(:)
!
!     note: need to set a flag for more options than eddy_time for the
!     next if statement   
!
!  Is the physical system advancing too slowly?
!
      if ((minshell .ne. -1) .and. (dt_ .gt. eddy_time(minshell))) then
        print*, 'Shell turn-over times unresolved (dt_>turnover)'
        print*, 'dt_=', dt_,' tturnover(nshell-2)=', eddy_time(minshell)
        call stop_it("")
      endif
!
!  Number of advances to call
!
      count=floor(dt_/deltG)
      if (count .le. 5) then
        print*, 'Shell timestepping too long (decrease deltaG)'
        print*, 'dt_/deltaG=',count
        call stop_it("")
      endif
!
!  Call advances
!
      do advances=1,count
        call GOY_advance_timestep()
      enddo
!
!  Ancillary stuff
!
      if (llast) then
!
!  Make certain all processors are on the same page
!
        call mpibcast_cmplx_arr(zu,nshell)
        if (lout) call &
            GOY_write_snapshot(trim(directory_snap)//'/GOYVAR',ENUM=.true., &
                FLIST='goyvarN.list')
      endif
!
      call keep_compiler_quiet(f,df)
!
!  Do we need to update real-space velocity vectors?
!
      if (.not. l_quenched) then
        call check_time(shell, update_vecs)
        if (update_vecs) call get_vectors(kvec, uvec, shell)
        call mpibcast_real(kvec,(/nshell,3/))
        call mpibcast_real(uvec,(/nshell,3/))
      endif
!
    endsubroutine  special_after_timestep
!***********************************************************************
    subroutine GOY_advance_timestep()
!
!  advances the shell model 1 internal timestep
!
      complex*16:: zt1,zt2,zt3,zgt
      integer :: ns
!
!!	  SHELL 1
      zt1 = exfac(1)*zu(1)
      zt2 = aval(1)*zu(2)*zu(3) + zforce
      zgt = zi*dconjg(zt2)
      zt3 = 0.50d0*(3.0D0*zgt - zgprev(1))
      zu(1) = zt1 + exrat(1)*zt3
      umod(1)=zabs(zu(1))
      zgprev(1) = zgt
!!	  SHELL 2
      zt1 = exfac(2)*zu(2)
      zt2 = aval(2)*zu(3)*zu(4) + bval(2)*zu(1)*zu(3)
      zgt = zi*dconjg(zt2)
      zt3 = 0.50d0*(3.0d0*zgt - zgprev(2))
      zu(2) = zt1 + exrat(2)*zt3
      umod(2)=zabs(zu(2))
      zgprev(2) = zgt
!!	  INNER (NSHELL-4) SHELLS
      do ns=3,nshell-2
        zt1 = exfac(ns)*zu(ns)
        zt2 =  aval(ns)*zu(ns+1)*zu(ns+2) + & 
                  bval(ns)*zu(ns-1)*zu(ns+1) + &
                  cval(ns)*zu(ns-1)*zu(ns-2)
        zgt = zi*dconjg(zt2)
        zt3 = 0.50d0*(3.0D0*zgt - zgprev(ns))
        zu(ns) = zt1 + exrat(ns)*zt3
        umod(ns)=zabs(zu(ns))
        zgprev(ns) = zgt
      enddo
!!	  (NSHELL-1)TH SHELL
      zt1 = exfac(nshell-1)*zu(nshell-1)
      zt2 = bval(nshell-1)*zu(nshell-2)*zu(nshell) + &
               cval(nshell-1)*zu(nshell-2)*zu(nshell-3)
      zgt = zi*dconjg(zt2)
      zt3 = 0.50d0*(3.0D0*zgt - zgprev(nshell-1))
      zu(nshell-1) = zt1 + exrat(nshell-1)*zt3
      umod(nshell-1)=zabs(zu(nshell-1))
      zgprev(nshell-1) = zgt
!!	  LAST SHELL
      zt1 = exfac(nshell)*zu(nshell)
      zt2 = cval(nshell)*zu(nshell-1)*zu(nshell-2)
      zgt = zi*dconjg(zt2)
      zt3 = 0.50d0*(3.0D0*zgt - zgprev(nshell))
      zu(nshell) = zt1 + exrat(nshell)*zt3
      umod(nshell)=zabs(zu(nshell))
      zgprev(nshell) = zgt
!
    endsubroutine GOY_advance_timestep
!***********************************************************************
    subroutine GOY_write_snapshot(chsnap,enum,flist)
!
!  Write uu snapshot to file.
!
      character (len=*) :: chsnap
      logical :: enum
      character (len=*), optional :: flist
!
      logical :: lsnap
!
      if (present(flist)) then
        call wsnap_GOY(chsnap,enum,lsnap,flist)
      else
        call wsnap_GOY(chsnap,enum,lsnap)
      endif
!
    endsubroutine GOY_write_snapshot
!***********************************************************************
    subroutine GOY_read_snapshot(filename)
!
!  Read uu snapshot from file.
!
      character (len=*) :: filename
      real :: t_goy
!
      call input_GOY(filename,t_goy)
!
    endsubroutine GOY_read_snapshot
!***********************************************************************
    subroutine wsnap_GOY(snapbase,enum,lsnap,flist)
!
      use General
      use Io
      use Sub
!
     character (len=*) :: snapbase, flist
     logical :: enum, lsnap
!
     integer, save :: ifirst=0, nsnap
     real, save :: tsnap, tsnap_goy
     character (len=fnlen), save :: fgoy
     character (len=fnlen) :: snapname
     character (len=5) :: nsnap_ch
!
     optional :: flist
!
     if (enum) then
!
       if (ifirst==0) then
          call safe_character_assign(fgoy,trim(datadir)//'/tsnap.dat')
          call read_snaptime(fgoy,tsnap,nsnap,dsnap,t)
          ifirst=1
        endif
!
        call update_snaptime(fgoy,tsnap,nsnap,dsnap,t,lsnap,nsnap_ch, ENUM=.true.)
        if (lsnap) then
          snapname=snapbase//nsnap_ch
          call output_GOY(snapname)
          if (ip<=10 .and. lroot) &
              print*,'wsnap_GOY: written snapshot ', snapname
          if (present(flist)) call log_filename_to_file(snapname,flist)
        endif
      else
        snapname=snapbase
        call output_GOY(snapname)
        if (ip<=10 .and. lroot) &
             print*,'wsnap_GOY: written snapshot ', snapname
        if (present(flist)) call log_filename_to_file(snapname,flist)
      endif
!
    endsubroutine wsnap_GOY
!***********************************************************************
    subroutine rsnap_GOY(filename)
!
      Use Mpicomm, only: mpibcast_cmplx_arr, mpibcast_real
!
      character (len=*) :: filename
      real :: t_goy
!
      if (lroot) then
        call input_GOY(filename,t_goy)
      endif
!
      call mpibcast_cmplx_arr(zu,nshell)
      call mpibcast_cmplx_arr(zgprev,nshell)
      if (l_needspecialuu) then
        call mpibcast_real(kvec, (/nshell,3/))
        call mpibcast_real(uvec, (/nshell,3/)) 
      endif
!
    endsubroutine rsnap_GOY
!***********************************************************************
    subroutine output_GOY(filename)
!
      integer :: lun_output=91
      real :: t_goy
      character (len=*) :: filename
!
      if (lroot) then
        t_goy=t
        open(lun_output,FILE=filename,FORM='unformatted')
        if (l_needspecialuu) then
          write(lun_output) t_goy,zu, zgprev, kvec, uvec
        else       
          write(lun_output) t_goy,zu, zgprev
        endif
        close(lun_output)
      endif
!
    endsubroutine output_GOY
!***********************************************************************
    subroutine input_GOY(filename,t_goy)
!
      Use Cdata
!
      integer :: lun_input=92
      real :: t_goy
      character (len=*) :: filename
!
      if (lroot) then
        open(lun_input,FILE=filename,FORM='unformatted')
        if (l_needspecialuu) then
          read(lun_input) t_goy, zu, zgprev, kvec, uvec
        else
          read(lun_input) t_goy, zu, zgprev
        endif
        close(lun_input)
      endif
!
    endsubroutine input_GOY
!***********************************************************************
    subroutine write_structure()
!
      integer :: lun_input=92, iord
      real :: t_goy
      character (120) :: filename
      real, dimension(nord,nshell) :: sf
!
      if (lroot) then
!
! calc new sf
!
        sf=0
        call calc_structure(sf)
!
! write sf
!
        do iord=1,nord
          call get_structure_filenames(filename,iord)
          t_goy=t
          open(UNIT=lun_input,FILE=trim(filename),&
              POSITION='append',FORM='unformatted')
          write(UNIT=lun_input), t_goy, sf(iord,:)
          close(UNIT=lun_input)
        enddo
      endif
!
    endsubroutine write_structure
!***********************************************************************
    subroutine calc_structure(sf)
!
      real, dimension(nord,nshell) :: sf
      integer :: ns, iord
      real :: sftemp
!
      do ns=1,nshell
        sftemp=1.
        do iord=1,nord
          sftemp=sftemp*umod(ns)
          sf(iord,ns)=sftemp
        enddo
      enddo
!
    endsubroutine calc_structure
!***********************************************************************
    subroutine get_structure_filenames(filename,iord)
!
      use General, only: chn, safe_character_assign
!
      character (len=*) :: filename
      character (len=5) :: ch
      integer :: iord
!
      call chn(iord,ch)
      call safe_character_assign(filename,trim(datadir)//&
          '/t_structure_'//trim(ch)//'.dat')
!
    endsubroutine get_structure_filenames
!***********************************************************************
    subroutine get_vectors(vec_array1, vec_array2, shell)
!
      real,dimension(nshell,3) :: vec_array1, vec_array2
      logical, dimension(nshell), optional :: shell
      integer :: ns
      real, dimension(3) :: a,b
!
! accumulate perpendicular pair/pairs of vectors,
! vec_array1=k-vectors, weighted by kval
! vec_array2=u-vectors, unit-length
!
      if (.not.present(shell)) then
        do ns=1, nshell
          call get_vector_pair(a,b)
          vec_array1(ns,:)=a*kval(ns)
          vec_array2(ns,:)=b
        enddo
      else
        do ns=1, nshell
          if (shell(ns)) then
            call get_vector_pair(a,b)
            vec_array1(ns,:)=a*kval(ns)
            vec_array2(ns,:)=b
          endif
        enddo
      endif
!
    endsubroutine get_vectors
!***********************************************************************
    subroutine get_vector_pair(vec1, vec2)
!
      real,dimension(3) :: vec1, vec2, a
      real :: norm
!
!  general perpendicular pair of unit vectors
!
!  get first unit vector
!
      call get_unit_vector(vec1)
!
! generate 2nd vector, if reasonable, project onto plane perpendicular
! to first and renormalize, otherwise try again
!
      do
        call get_unit_vector(a)
        vec2=a-dot_product(a,vec1)*vec1
        norm=sum(vec2**2)
        if (norm.gt.0.01) exit
      enddo
      vec2=vec2/(norm**(0.5)) 
!
    endsubroutine get_vector_pair
!***********************************************************************
    subroutine get_unit_vector(vec1)
!
! generate random unit vector
!
      use General, only: random_number_wrapper
!
      real,dimension(3) :: vec1
      real :: norm
!
      do
        call random_number_wrapper(vec1)
        vec1=vec1-0.5                        !generate random vector in unit cube
        norm=sum(vec1**2)
        if ((norm.le.0.25).and.(norm.gt.0.01)) exit !if in half-unit sphere, keep
      enddo                                                  !otherwise try again
      vec1=vec1/(norm**(0.5)) !normalize
!
    endsubroutine get_unit_vector
!***********************************************************************
    subroutine calc_gas_velocity_GOY(xpos, uup, shell1,shell2)
!
! provide velocity at point xpos
!
      real, dimension(3) :: xpos, uup
      integer, optional :: shell1, shell2
!
      integer :: ns, nstart, nend
!
      uup=0
      if (present(shell1)) then; nstart=shell1; else; nstart=3;  endif
      if (present(shell2)) then; nend=shell2;   else; nend=ns-2; endif
      do ns=nstart,nend
        uup=uup+uvec(ns,:)*2*(&
            real(zu(ns))*cos(dot_product(xpos,kvec(ns,:)))-&
            aimag(zu(ns))*sin(dot_product(xpos,kvec(ns,:))))
      enddo

    endsubroutine calc_gas_velocity_GOY
!***********************************************************************
    subroutine check_time(shell, update_vecs)
!
! determine whether we need to update real-space velocity vectors
!
      logical :: update_vecs
      logical, dimension(nshell) :: shell
!
      integer :: ns
!
! when do we need to update?  First set up update times
!
      if (headt) call new_update_times(tnext_vel)
!
      update_vecs=.false.
      shell(:)=.false.
      do ns=3,nshell-2       
        if (t .gt. tnext_vel(ns)) then
          shell(ns)=.true.   !this shell update needed
          update_vecs=.true. !any updates needed
        endif
      enddo
!
!  If we will be updating, also update the time for next update
!
      if (update_vecs) call new_update_times(tnext_vel, shell)
!
    endsubroutine check_time
!***********************************************************************
    subroutine new_update_times(tnext, shell)

      real, dimension(nshell) :: tnext
      logical, dimension(nshell), optional :: shell
      integer :: ns
!
! get new real-space velocity vector update times
!
      if (.not.present(shell)) then
        do ns=3,nshell-2
          tnext(ns)=t+eddy_time(ns)
        enddo
      else
        do ns=3, nshell-2
          if (shell(ns)) tnext(ns)=t+eddy_time(ns)
        enddo
      endif
!
    endsubroutine new_update_times
!***********************************************************************
    subroutine get_largest_turnover(t_largestturnover)
!
! share largest meaningful turnover time with the world if they ask for it
!
      Use Mpicomm, only:stop_it
!
      real :: t_largestturnover
!
      if (maxshell.eq. -1) then
        print*, 'get_largest_turnover called without maxshell set'
        call stop_it("")
      else
        t_largestturnover=eddy_time(maxshell)
      endif
!
    endsubroutine get_largest_turnover
!***********************************************************************
    subroutine special_clean_up()
!
      deallocate(zu)
      deallocate(kval,aval,bval, exfac, exrat)
      deallocate(cval,ksquare, umod, eddy_time, zgprev)
      if (l_needspecialuu) deallocate(kvec,uvec)
      if (l_needspecialuu .and. .not. l_quenched) deallocate(tnext_vel)
!
    endsubroutine special_clean_up
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
