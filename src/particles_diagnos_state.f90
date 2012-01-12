! $Id$
!
!  This module tracks the evolution of particles in terms of
!  user defined "states".  Currently the only implemented state is in terms
!  of the local gas velocity.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 5
! MAUX CONTRIBUTION 0
!
! CPARAM logical, parameter :: lparticles_diagnos_state=.true.
!
!
!***************************************************************
module Particles_diagnos_state
!
  use Cdata
  use Cparam
  use Messages
  use Particles_cdata
  use Particles_mpicomm
  use Particles_sub
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_diagnos_state.h'
!
  character (len=labellen) :: state_def ='velz'
  real :: slow_vel_trigger=0.,slow_vel_trigger2=0.
!
  namelist /particles_diagnos_state_run_pars/ &
      slow_vel_trigger, state_def
!
  integer :: idiag_upm=0, idiag_vprms=0, idiag_uzplus=0, idiag_uzminus=0
  integer :: idiag_uzminuscount=0, idiag_uzpluscount=0
!
  contains
!***********************************************************************
    subroutine register_pars_diagnos_state()
!
!  Set up indices for access to the fp and dfp arrays
!
      if (lroot) call svn_id( &
           "$Id$")
!
      ipss=npvar+1  !particle state
      ipst=npvar+2  !time particle entered current state
!
      ipxx=npvar+3  !x position particle entered current state
      ipyy=npvar+4  !y position particle entered current state
      ipzz=npvar+5  !z position particle entered current state
!
!     individual particle labels are tracked through ipar,
!
      npvar=npvar+5
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_pars_diagnos_state', &
            'npvar > mpvar')
      endif
!
    endsubroutine register_pars_diagnos_state
!***********************************************************************
    subroutine initialize_pars_diagnos_state(f,lstarting)
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      slow_vel_trigger2=slow_vel_trigger**2
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_pars_diagnos_state
!***********************************************************************
    subroutine init_particles_diagnos_state(fp)
!
      real, dimension (mpar_loc,mpvar) :: fp
!
      integer :: k
!
      intent (out) :: fp
!
      do k=1,npar_loc
        fp(k,ipst)=t
        fp(k,ipss)=-1.
        fp(k,ipxx:ipzz)=fp(k,ixp:izp)
      enddo
!
    endsubroutine init_particles_diagnos_state
!***********************************************************************
    subroutine insert_particles_diagnos_state(fp, npar_loc_old)
!
! Allow for the insertion of particles
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: npar_loc_old
!
      integer :: k
!
      intent (inout) :: fp
!
      do k=npar_loc_old+1,npar_loc
        fp(k,ipst)=t
        fp(k,ipss)=-1.
        fp(k,ipxx:ipzz)=fp(k,ixp:izp)
      enddo
!
    endsubroutine insert_particles_diagnos_state
!***********************************************************************
    subroutine read_pars_diag_state_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_diagnos_state_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_diagnos_state_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_pars_diag_state_run_pars
!***********************************************************************
    subroutine write_pars_diag_state_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_diagnos_state_run_pars)
!
    endsubroutine write_pars_diag_state_run_pars
!***********************************************************************
    subroutine rprint_particles_diagnos_state(lreset,lwrite)
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname,inamez,inamey,inamex,inamexy,inamexz,inamer,inamerz
      logical :: lwr
!
!  Write information to index.pro.
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if (lwr) then
        write(3,*) 'ipss=', ipss
        write(3,*) 'ipst=', ipst
        write(3,*) 'ipxx=', ipxx
        write(3,*) 'ipyy=', ipyy
        write(3,*) 'ipzz=', ipzz
      endif
!
!  Reset everything in case of reset.
!
      if (lreset) then
        idiag_upm=0; idiag_uzminus=0; idiag_uzplus=0
        idiag_uzminuscount=0; idiag_uzpluscount=0
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'upm',idiag_upm)
        call parse_name(iname,cname(iname),cform(iname),'uzminus',idiag_uzminus)
        call parse_name(iname,cname(iname),cform(iname),'uzplus',idiag_uzplus)
        call parse_name(iname,cname(iname),cform(iname),'uzminuscount',idiag_uzminuscount)
        call parse_name(iname,cname(iname),cform(iname),'uzpluscount',idiag_uzpluscount)
      enddo
!
    endsubroutine rprint_particles_diagnos_state
!***********************************************************************
    subroutine persistence_check(fp, k, uup)
!
      real, dimension (mpar_loc,mpvar) :: fp
      real :: r_new
      integer :: k
      real, dimension(3) :: uup
      logical :: l_swap
!
!  check for state change
!
      l_swap=.false.
      call evaluate_state(fp, k, uup, l_swap, r_new)
!
!  execute state change, write state change to file
!
      if (l_swap) then
        if (.not. (fp(k,ipss) < 0)) call &
            data_store_persistence(fp, k, r_new, .false.)
        fp(k,ipss)=r_new
        fp(k,ipst)=t
        fp(k,ipzz)=fp(k,izp)
      endif
!
!  if final particle in pencil, close file
!
      if (k==k2_imn(imn)) &
          call data_store_persistence(fp, -1, -1., .true.)
!
!  some state diagnostics
!
      if (ldiagnos) then
        if (idiag_upm/=0) &
            call sum_par_name((/sum((fp(k,ivpx:ivpz)-uup)**2)/),idiag_upm,lsqrt=.true.)
        if (idiag_uzminus/=0) then
          if (uup(3) < 0) then
            call sum_par_name((/ uup(3) /),idiag_uzminus)
            call sum_par_name((/ 1. /),idiag_uzminuscount)
          else
            call sum_par_name((/ 0. /),idiag_uzminus)
            call sum_par_name((/ 0. /),idiag_uzminuscount)
          endif
        endif
        if (idiag_uzplus/=0) then
          if (uup(3) >= 0) then
            call sum_par_name((/ uup(3) /),idiag_uzplus)
            call sum_par_name((/ 1. /),idiag_uzpluscount)
          else
            call sum_par_name((/ 0. /),idiag_uzplus)
            call sum_par_name((/ 0. /),idiag_uzpluscount)
          endif
        endif
      endif
!
     endsubroutine persistence_check
!***********************************************************************
    subroutine evaluate_state(fp, k, uup, l_swap, r_new)
!
      real, dimension (mpar_loc,mpvar) :: fp
      real :: r_new
      integer :: k
      real, dimension(3) :: uup
      logical :: l_swap
!
! as more cases are added, the call to
! evaluate_state may need to be moved/changed/elaborated on
!
      select case(state_def)
      case ('velz') !gas z velocity up or down in lab frame (particle settling)
        if (uup(3) > 0) then
          r_new=2
        else
          r_new=1
        endif
      case default
        call fatal_error('evaluate_state','No such such value for '// &
          'state_def: '//trim(state_def))
      endselect
!
      if (abs(fp(k, ipss)-r_new) > 0.1) then
        l_swap=.true.
      elseif (abs(fp(k, ipss)-r_new) <= 0.1) then
        l_swap=.false.
      endif
!
    endsubroutine evaluate_state
!***********************************************************************
    subroutine data_store_persistence(fp, k, r_new, l_close)
!
      use General, only: safe_character_assign, itoa
!
      real, dimension (mpar_loc,mpvar) :: fp
      real :: r_new
      integer :: k
      logical :: l_close                   !do we close the output file?
!
      logical, save       :: l_opened=.false. !is state output file open
!
      integer             :: lun_input=91.
      real                :: time_span
      real, dimension(3)  :: travel
      character (len=128) :: filename
!
!  Close file if told to and it is opened
!  Open file if it is closed, and not told to close it
!
      if (l_close .and. l_opened) then
        close(UNIT=lun_input)
        l_opened=.false.
      elseif ((.not. l_close) .and. (.not. l_opened)) then
        call safe_character_assign(filename,trim(datadir)//&
            '/persistence'//trim(itoa(iproc))//'.dat')
        open(UNIT=lun_input,FILE=trim(filename),&
            POSITION='append',FORM='formatted')
        l_opened=.true.
      endif
!
      if (l_opened) then
        time_span=t-fp(k,ipst)
        travel=fp(k,ixp:izp)-fp(k,ipxx:ipzz)
        write(UNIT=lun_input,FMT= &
            '(f12.5, i9, f4.1, f12.5, f4.1, f12.5, f12.5, f12.5)') &
            t, ipar(k), fp(k,ipss),&
            time_span, r_new, travel(1), travel(2), travel(3)
!
!  ie: current time, particle label, old state, time spend in old state,
!      new state, distance traveled in old state (x,y,z)
!
      endif
!
    endsubroutine data_store_persistence
!***********************************************************************
endmodule Particles_diagnos_state
