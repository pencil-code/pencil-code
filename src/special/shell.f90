! $Id$
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
  real    :: nstability1=1e6,nstability2=1e6,nstability3=1e6  !timesteps during start
  complex*16 :: zforce
  integer :: stability_output=100000
  integer :: i_starting_local=0
!
! output parameters
!
  logical :: lstructure=.true.                                !output structure functions
  integer :: nord=1                                           !order of structure function output
  logical :: l_altu=.false.
!
! real-space parameters
!
  logical :: l_needspecialuu=.false.                          !Need real-space velocity
  logical :: l_quenched=.true.                                !NOT update real-space velocity vectors
  logical :: l_nogoy=.false.
!
! runtime variables
!
! aval,bval,cval= k_n*a_n, k_n*b_n/k, k_n*c_n/k^2 from Jensen 1999
!
  real*8, allocatable, dimension(:) :: aval, bval, cval       !evolution equation coefficients
  real*8, allocatable, dimension(:) :: exfac, exrat           !evolution equation variables
  real*8, allocatable, dimension(:) :: umod, umalt,eddy_time  !velocity magnitude and turnover times
  real, allocatable, dimension(:) :: uav
  real*8, allocatable, dimension(:) :: kval, ksquare          !k vectors
  complex*16, allocatable, dimension(:) :: zu, zgprev         !fourier space velocities
!
! real-space variables
!
  real, allocatable, dimension(:,:,:) :: kvec, uvec      !to generate real-space velocities
  real, allocatable, dimension(:,:)   :: vec_split, vec_split0
!
  character (len=labellen) :: time_step_type='nothing'        !Shell contribution to time-step
  integer :: minshell=-1, maxshell=-1, minshell2=-1           !bounds of physically relevant
                                                              !shell (time-step contraints)
!
  real, dimension(:),pointer:: uup_shared                            !shared variable (gas velocity)
  real,pointer :: turnover_shared
  logical,pointer :: vel_call, turnover_call
!
  namelist /special_init_pars/ &
      nshell, visc, f_amp, lambda, deltG, k_amp, u_amp,&
      l_needspecialuu, l_altu, l_quenched, stability_output,&
      nstability1, nstability2, nstability3, minshell, maxshell,&
      minshell2, l_nogoy
!
  namelist /special_run_pars/ &
      lstructure, nord,time_step_type 
!
!! Diagnostic variables (needs to be consistent with reset list below).
!
   integer :: idiag_uuGOY=0, idiag_deltM=0, idiag_eddy_time=0
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
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
      Use Mpicomm, only: mpibcast_cmplx_arr, stop_it, mpifinalize
      use General
      use Sub, only: cross
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting, kdat_exists
      integer :: ns, ierr, nk,nk2, ik, count
      character (len=50) :: filename
      character (len=5) :: ch
      real :: kav
      real, dimension(1) :: fran
      real, dimension(2000) :: kkx, kky, kkz
!
! allocate arrays (init_special calls initialize_special)
!
      if(.not.allocated(zu)) then
        allocate(zu(nshell),    zgprev(nshell), eddy_time(nshell))
        allocate(kval(nshell),  umod(nshell),   aval(nshell), bval(nshell))
        allocate(cval(nshell),  exfac(nshell),  exrat(nshell),ksquare(nshell))
        allocate(uav(nshell+1))
        if (l_altu) allocate(umalt(nshell))
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
! set minshell2
!
      if (minshell2 .lt. 0) minshell2=minshell
!
!  generate/allocate k/u vectors (if needed, and this is called by start)
!

      if (l_needspecialuu) then
        if ((maxshell .lt. 0) .or. (minshell .lt. 0)) then
          print*, 'Set maxshell and minshell if you want l_needspecialuu'
          call stop_it("")
        else
!
!  prep k/u vectors
!
          if (.not. allocated(kvec)) then
            allocate(kvec(3,nshell,3),uvec(3,nshell,3), vec_split(nshell,3))
            if (.not. l_quenched) &
                allocate(vec_split0(nshell,3))
          endif
          if (i_starting_local .eq. 1) then; do ns=maxshell, minshell2
!
            call get_unit_vector(vec_split(ns,:))
!
!  read in possible vectors
!
            call chn(int(kval(ns)),ch)
            call safe_character_assign(filename,'k'//trim(ch)//'.dat')
            inquire(FILE=filename, EXIST=kdat_exists)
            if (.not. kdat_exists) then              
              print*, 'k.dat file ',filename,' does not exist, exiting'
              call stop_it("")
            else
              open(9,file=filename,status='old')
              read(9,*) nk,kav
              if (nk>2000) then
                print*,'too many kvectors (shell.f90)'
                call mpifinalize
              endif
              read(9,*) (kkx(ik),ik=1,nk)
              read(9,*) (kky(ik),ik=1,nk)
              read(9,*) (kkz(ik),ik=1,nk)
              close(9)
            endif
!
!  choose k vectors, u vectors
!  1st pair
!
            call random_number_wrapper(fran)
            ik=nk*(.9999*fran(1))+1
            kvec(1,ns,:)=(/kkx(ik),kky(ik),kkz(ik)/)
            call get_perp_vector(kvec(1, ns,:), uvec(1, ns,:))
!
! 2nd/3rd pairs
!
            do count=2,3
              nk2=1
              do
                call random_number_wrapper(fran)
                ik=nk*(.9999*fran(1))+1
                kvec(count,ns,:)=(/kkx(ik),kky(ik),kkz(ik)/)
                if (abs(dot_product(uvec(count-1,ns,:),kvec(count,ns,:))) .gt. 0.86) exit
                if (nk2 .gt. nk) then
                  call get_perp_vector(kvec(count-1, ns,:), uvec(count-1, ns,:))
                  nk2=0
                endif
                nk2=nk2+1
              enddo 
              call cross(kvec(count-1,ns,:),kvec(count,ns,:), uvec(count,ns,:))
              uvec(count,ns,:)=uvec(count,ns,:)/(sum(uvec(count,ns,:)**2))**(0.5)
            enddo
          enddo; endif
        endif
!
      endif
!
!  end l_needspecialuu section
!
!  read goydata, bcast etc... (if not starting)
!
      if (.not. lstarting) then
        call GOY_read_snapshot(trim(directory_snap)//'/goyvar.dat')
        umod=abs(zu); if (l_altu)  call calc_altu(); call calc_eddy_time()
        call GOY_write_snapshot(trim(directory_snap)//'/GOYVAR',ENUM=.true., &
            FLIST='goyvarN.list',snapnum=0)
        call GOY_write_snapshot(trim(directory_snap)//'/goyvar.dat',ENUM=.false.)
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
      real :: tstarting
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
        i_starting_local=1
        call initialize_special(f,lstarting=.true.)
        i_starting_local=2
        do ns=1,nshell
          umod(ns)=u_amp*kval(ns)**noneoverthree
          zu(ns)=dcmplx(umod(ns)*dcos(randomd(ns)), &
              umod(ns)*dsin(randomd(ns)))
          zgprev=dcmplx(0.,0.)
        enddo
        umod=abs(zu)
        if (l_altu) call calc_altu()
        uav=0
        eddy_time=0
!
! if we need real-space velocities, 
!
        if (l_needspecialuu .and. .not. l_quenched) vec_split0=vec_split
!
        call write_structure(tstarting=0.)
!
! stabilize velocities
!
        count=0
        do while (count .lt. nstability1)
          call GOY_advance_timestep()
          count=count+1
        enddo
!
! start tracking the velocities to get eddy time estimates
!
        count=0
        do while (count .lt. nstability2)
          call GOY_advance_timestep()
          if (mod(count, stability_output)==0) then
            if (l_altu) call calc_altu()
            call calc_eddy_time()
          endif
          count=count+1
        enddo
!
! more stabilization, with some data accumulation thrown in
!
        count=0
        do while (count .lt. nstability3)
          call GOY_advance_timestep()
          if (.not. l_quenched) call update_vec(vec_split,sngl(deltG))
          count=count+1
          if (mod(count, stability_output)==0) then
            print*, 'stabilizing, count=',count
            if (l_altu) call calc_altu()
            call calc_eddy_time()
            tstarting=deltG*count
            call write_structure(tstarting=tstarting)
          endif
        enddo
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
      if (lfirstpoint .and. lfirst) then
        if (lroot) then
          select case (time_step_type)
            case ('smallest')
              minturn=minval(eddy_time(maxshell:minshell))
              deltmbase=min(cdt*minturn,0.08/visc/ksquare(minshell))
            case ('largest')
              if (l_altu) then
                 minturn=1./maxval(umalt(maxshell:minshell))/kval(minshell)
              else
                minturn=1./umod(maxshell)/kval(minshell)  !crossing time of smallest eddies
              endif
              deltmbase=min(cdt*minturn,0.08/visc/ksquare(minshell))
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
        if (idiag_deltM/=0)     call max_name(deltmbase,idiag_deltM)
        if (idiag_uuGOY/=0) then
          if (l_altu) then
            call save_name(sngl(umalt(maxshell)),idiag_uuGOY) 
          else
            call save_name(sngl(umod(maxshell)),idiag_uuGOY)
          endif
        endif
        if (idiag_deltM/=0)     call max_name(deltmbase,idiag_deltM)       
        if (idiag_eddy_time/=0) call save_name(sngl(eddy_time(minshell)),idiag_eddy_time)   
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
        idiag_uuGOY=0
        idiag_deltM=0
        idiag_eddy_time=0
      endif
!!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'uuGOY',idiag_uuGOY)
        call parse_name(iname,cname(iname),cform(iname),'deltM',idiag_deltM)
        call parse_name(iname,cname(iname),cform(iname),'eddy_time',idiag_eddy_time)
      enddo
!!
!!!  write column where which magnetic variable is stored
      if (lwr) then
        write(3,*) 'i_uuGOY=',idiag_uuGOY
        write(3,*) 'i_deltM=',idiag_deltM
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
      use SharedVariables, only: get_shared_variable
!
!  Called before the loop, in case some particle value is needed
!  for the special density/hydro/magnetic/entropy.
!
!  Also hijacked to get gas velocity using shared variable.
!  shared logical vel_call triggers the hijack
!
      real, dimension (:,:), intent(in) :: fp
!
      real, dimension(3) :: xpos
      logical, save :: first_call=.true.
      integer :: ierr
!
      if (first_call .and. l_needspecialuu) then
        call get_shared_variable('uup_shared',uup_shared,ierr)
        call get_shared_variable('vel_call',vel_call,ierr)
        call get_shared_variable('turnover_call',turnover_call,ierr)
        call get_shared_variable('turnover_shared',turnover_shared,ierr)
        first_call=.false.
print*,'special_calc_particles, shared variable,iproc=',iproc
      endif
!
      if (vel_call) then
        xpos=uup_shared
        call calc_gas_velocity(xpos, uup_shared)
        vel_call=.false.
      endif
!
      if (turnover_call) then
        call get_largest_turnover(turnover_shared)
        turnover_call=.false.
      endif
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
!  Is the physical system advancing too slowly?
!
      if ((minshell .ne. -1) .and. (dt_ .gt. eddy_time(minshell))) then
        print*, 'Shell turn-over times unresolved (dt_>turnover)'
        print*, 'dt_=', dt_,' tturnover(minshell)=', eddy_time(minshell)
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
      if (lroot) then
        do advances=1,count
          if (.not. l_nogoy) &
              call GOY_advance_timestep()
          if (.not. l_quenched) then
            call update_vec(vec_split, sngl(deltG))
          endif
        enddo
      endif
!
      call mpibcast_cmplx_arr(zu,nshell)
      call mpibcast_cmplx_arr(zgprev,nshell)
      umod=zabs(zu)
      if (l_altu) call calc_altu()
      call calc_eddy_time()
      if (l_needspecialuu) then
        call mpibcast_real(vec_split, (/nshell,3/))
        if (.not. l_quenched) call mpibcast_real(vec_split0, (/nshell,3/))
      endif
!
!  Ancillary stuff
!
      if (llast) then
!
        if (lout) call &
            GOY_write_snapshot(trim(directory_snap)//'/GOYVAR',ENUM=.true., &
                FLIST='goyvarN.list')
!
        if (mod(it,isave)==0) call &
            GOY_write_snapshot(trim(directory_snap)//'/goyvar.dat',ENUM=.false.)
      endif
!
      call keep_compiler_quiet(f,df)
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
    subroutine GOY_write_snapshot(chsnap,enum,flist,snapnum)
!
!  Write uu snapshot to file.
!
      character (len=*) :: chsnap
      logical :: enum
      character (len=*), optional :: flist
      integer, optional :: snapnum
!
      logical :: lsnap
!
      if (present(snapnum)) then
        if (present(flist)) then
          call wsnap_GOY(chsnap,enum,lsnap,flist, snapnum=snapnum)
        else
          call wsnap_GOY(chsnap,enum,lsnap,snapnum=snapnum)
        endif
      else
        if (present(flist)) then
          call wsnap_GOY(chsnap,enum,lsnap,flist)
        else
          call wsnap_GOY(chsnap,enum,lsnap)
        endif
      endif
!
    endsubroutine GOY_write_snapshot
!***********************************************************************
    subroutine GOY_read_snapshot(filename)
!
!  Read uu snapshot from file.
!
      character (len=*) :: filename
!
      call rsnap_GOY(filename)
!
    endsubroutine GOY_read_snapshot
!***********************************************************************
    subroutine wsnap_GOY(snapbase,enum,lsnap,flist, snapnum)
!
      use General
      use Io
      use Sub
!
     character (len=*) :: snapbase, flist
     logical :: enum, lsnap
     integer, optional :: snapnum
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
       if (present(snapnum)) then
         call chn(snapnum,nsnap_ch)
         lsnap=.true.
       else
         if (ifirst==0) then
            call safe_character_assign(fgoy,trim(datadir)//'/tsnap.dat')
            call read_snaptime(fgoy,tsnap,nsnap,dsnap,t)
            ifirst=1
         endif
!
         call update_snaptime(fgoy,tsnap,nsnap,dsnap,t,lsnap,nsnap_ch, ENUM=.true.)
       endif
!
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
      if (lroot) call input_GOY(filename,t_goy)
!
      call mpibcast_cmplx_arr(zu,nshell)
      call mpibcast_cmplx_arr(zgprev,nshell)
      call mpibcast_real(uav, nshell+1)
      if (l_needspecialuu) then
        call mpibcast_real(kvec, (/3,nshell,3/))
        call mpibcast_real(uvec, (/3,nshell,3/)) 
        call mpibcast_real(vec_split, (/nshell,3/))
        if (.not. l_quenched) call mpibcast_real(vec_split0, (/nshell,3/))
      endif
      umod=abs(zu)
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
        if (l_needspecialuu .and. .not. l_quenched) then
          write(lun_output) t_goy,zu, zgprev, kvec, uvec, uav, vec_split, vec_split0
        else if (l_needspecialuu) then
          write(lun_output) t_goy,zu, zgprev, kvec, uvec, uav, vec_split
        else
          write(lun_output) t_goy,zu, zgprev, uav
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
        if (l_needspecialuu .and. .not. l_quenched) then
          read(lun_input)   t_goy,zu, zgprev, kvec, uvec, uav, vec_split, vec_split0
        else if (l_needspecialuu) then
          read(lun_input)   t_goy,zu, zgprev, kvec, uvec, uav, vec_split
        else
          read(lun_input)   t_goy,zu, zgprev, uav
        endif
        close(lun_input)
      endif
!
    endsubroutine input_GOY
!***********************************************************************
    subroutine write_structure(tstarting)
!
      real, optional :: tstarting
!
      integer :: lun_input=92, iord
      real :: t_goy
      character (120) :: filename
      real, dimension(nord,nshell) :: sf
      real, dimension(nshell) :: phase, vecdot
!
      if (lroot) then
!
! calc new sf
!
        sf=0
        call calc_structure(sf, phase, vecdot)
        iord=1
!
! write sf
!
        t_goy=t
        if (present(tstarting)) t_goy=tstarting
        call get_structure_filenames(filename,iord)
        open(UNIT=lun_input,FILE=trim(filename),&
            POSITION='append',FORM='unformatted')
        if (.not. l_quenched) then
          write(UNIT=lun_input), t_goy, sf(1,:), phase(:),vecdot(:), eddy_time(:)
        else
          write(UNIT=lun_input), t_goy, sf(1,:), phase(:)
        endif
        close(UNIT=lun_input)
        if (nord .gt. 1) then
          do iord=2,nord
            call get_structure_filenames(filename,iord)
            open(UNIT=lun_input,FILE=trim(filename),&
                POSITION='append',FORM='unformatted')
            write(UNIT=lun_input), t_goy, sf(iord,:)
            close(UNIT=lun_input)
          enddo
        endif
      endif
!
    endsubroutine write_structure
!***********************************************************************
    subroutine calc_structure(sf, phase, vecdot)
!
      real, dimension(nord,nshell) :: sf
      real, dimension(nshell) :: phase, vecdot
      integer :: ns, iord
      real :: sftemp
!
      do ns=1,nshell
        sftemp=1.
        do iord=1,nord
          if (l_altu) then
            sftemp=sftemp*umalt(ns)
          else
            sftemp=sftemp*umod(ns)
          endif
          sf(iord,ns)=sftemp
        enddo
        phase(ns)=real(zu(ns))/umod(ns)
        if (.not. l_quenched) then
          vecdot(ns)=dot_product(vec_split0(ns,:),vec_split(ns,:))
        endif
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
    subroutine get_perp_vector(vec1, vec2)
!
      real,dimension(3) :: vec1, vec2, a
      real :: norm
!
!  general vec2 unit, perpendicular to vec1
!
! generate 2nd vector, if reasonable, project onto plane perpendicular
! to first and renormalize, otherwise try again
!
      do
        call get_unit_vector(a)
        vec2=a-dot_product(a,vec1)*vec1/sum(vec1**2)
        norm=sum(vec2**2)
        if (norm.gt.0.01) exit
      enddo
      vec2=vec2/(norm**(0.5)) 
!
    endsubroutine get_perp_vector
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
    subroutine update_vec(vec_array, dt_)
!
!  update vectors
!
      use Sub, only: cross
!
      real,dimension(nshell,3) :: vec_array
      real :: dt_
!
      integer :: ns
      real, dimension(3) :: a, uv
      real :: norm, scale
!
      do ns=maxshell,minshell2
!
        scale=pi*((dt_/eddy_time(ns))**(0.5))
        do
          uv=vec_array(ns,:)
          call get_unit_vector(a)
!
! scale a for random walk (|a| = pi sqrt(dt_/eddy_time))
!
          scale=pi*((dt_/eddy_time(ns))**(0.5))
          a=a*scale
!
          uv=uv+a
          norm=sum(uv**2);
          if (norm .gt. 0.01) exit
        enddo
!
        vec_array(ns,:)=uv/(norm**(0.5)) 
!
      enddo
!
    endsubroutine update_vec
!***********************************************************************
    subroutine calc_gas_velocity(xpos, uup, shell1,shell2)
!
! provide velocity at point xpos
!
      real, dimension(3) :: xpos, uup
      integer, optional :: shell1, shell2
!
      integer :: ns, nstart, nend, count
!
      uup=0
      if (present(shell1)) then; nstart=shell1; else; nstart=maxshell;  endif
      if (present(shell2)) then; nend=shell2
      else if (minshell2 .gt. 0) then
        nend=minshell2
      else
        nend=minshell
      endif
!
      do ns=nstart,nend
        if (l_altu) then; do count=1,3
          uup=uup+vec_split(ns,count)*uvec(count,ns,:)*2**(0.5)*(&
              real(zu(ns))*cos(dot_product(xpos,kvec(count,ns,:)))-&
              aimag(zu(ns))*sin(dot_product(xpos,kvec(count,ns,:))))*&
              umalt(ns)/umod(ns) 
        enddo; else; do count=1,3
          uup=uup+vec_split(ns,count)*uvec(count,ns,:)*2**(0.5)*(&
              real(zu(ns))*cos(dot_product(xpos,kvec(count,ns,:)))-&
              aimag(zu(ns))*sin(dot_product(xpos,kvec(count,ns,:))))
        enddo; endif
      enddo
!
    endsubroutine calc_gas_velocity
!***********************************************************************
    subroutine get_largest_turnover(t_largestturnover)
!
! share largest meaningful turnover time with the world if they ask for it
!
      Use Mpicomm, only: stop_it
!
      real :: t_largestturnover
!
      if (maxshell.eq. -1) then
        print*, 'get_largest_turnover called without maxshell set'
        call stop_it("")
      else
        t_largestturnover=maxval(eddy_time(maxshell:minshell))
      endif
!
    endsubroutine get_largest_turnover
!***********************************************************************
    subroutine calc_altu()
!
      integer :: ns
!
      umalt(1)=umod(1); umalt(nshell-1:nshell)=umod(nshell-1:nshell)
      do ns=2,nshell-2
        umalt(ns)=abs(aimag(zu(ns+2)*zu(ns+1)*zu(ns)-(1./4.)* &
            zu(ns-1)*zu(ns)*zu(ns+1)))**(1./3.)
      enddo
    endsubroutine calc_altu
!***********************************************************************
    subroutine calc_eddy_time()
!
      integer :: ns
!
      uav(nshell+1)=uav(nshell+1)+1
      if (l_altu) then
        uav(1:nshell)=(uav(1:nshell)*(uav(nshell+1)-1)+umalt(:))/uav(nshell+1)
      else
        uav(1:nshell)=(uav(1:nshell)*(uav(nshell+1)-1)+umod(:))/uav(nshell+1)
      endif
      eddy_time=1./kval/uav(1:nshell)
!
    endsubroutine calc_eddy_time
!***********************************************************************
    subroutine special_clean_up()
!
      call GOY_write_snapshot(trim(directory_snap)//'/goyvar.dat',ENUM=.false.)
!
      deallocate(zu)
      deallocate(kval,aval,bval, exfac, exrat)
      deallocate(cval,ksquare, umod, eddy_time, zgprev)
      if (l_needspecialuu) then
        deallocate(kvec,uvec)
      endif
      if (l_altu) deallocate(umalt)
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
