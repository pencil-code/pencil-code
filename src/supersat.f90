! $Id$
!
! This module is used to solve the equation of supersaturation
! for either the Smoluchowski approach or the swarm model.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lsupersat = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ssat
! PENCILS PROVIDED gssat(3); ugssat
! PENCILS PROVIDED del2ssat
!***************************************************************
module Supersat
!
  use Cdata
  use Cparam
  use Messages
!
  implicit none
!
  include 'supersat.h'
!
!  Init parameters.
!
  real :: ssat_const=0.0
  logical :: nosupersat=.false., reinitialize_ssat=.false.
  character (len=labellen) :: initssat='impossible'
!
  namelist /supersat_init_pars/ &
           initssat, ssat_const
!
!  Run parameters.
!
  real :: supersat_diff=0.0
  real :: supersat_sink=0.0, Rsupersat_sink=0.5
  real, dimension(3) :: gradC0=(/0.0,0.0,0.0/)
  logical :: lsupersat_sink=.false.
  logical :: lupw_ssat=.false.

  namelist /supersat_run_pars/ &
      lupw_ssat, lsupersat_sink, Rsupersat_sink, supersat_sink, &
      supersat_diff, gradC0
!
! Declare index of new variables in f array
!
!  integer :: issat=0
! XY: type of "issat" is defined in "cdata.f90"
!
!  Diagnostics variables
!
  integer :: idiag_ssatrms=0
!
  contains
!***********************************************************************
    subroutine register_supersat()
!
!  Initialise the ssat variable and increase nvar accordingly
!
!   3-jun-16/xiangyu: adapted from pscalar_nolog
!
      use FArrayManager
!
      call farray_register_pde('ssat', issat)
      issat = 0                 ! needed for idl
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_supersat
!***********************************************************************
    subroutine initialize_supersat(f)
!
!  Perform any necessary post-parameter read initialization
!  Since the passive scalar is often used for diagnostic purposes
!  one may want to reinitialize it to its initial distribution.
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!
      if (lroot) print*, 'Supersaturation routine'
!
!  set to zero and then call the same initial condition
!  that was used in start.csh
!
      if (reinitialize_ssat) then
        f(:,:,:,issat)=0.
        call init_ssat(f)
      endif
!
    endsubroutine initialize_supersat
!
!***********************************************************************
    subroutine init_ssat(f)
!  initialise passive scalar field; called from start.f90
      use Sub
      use Initcond
!??AXEL use InitialCondition, only: initial_condition_ssat
!
      real, dimension (mx,my,mz,mfarray) :: f
       select case (initssat)
        case ('nothing')
        case ('zero'); f(:,:,:,issat)=0.0
        case ('constant'); f(:,:,:,issat)=ssat_const
       endselect
    endsubroutine init_ssat
   
!***********************************************************************
    subroutine pencil_criteria_supersat()
      lpenc_requested(i_ssat)=.true.
            
      if (lsupersat_sink) lpenc_requested(i_ssat)=.true.
      if (supersat_diff/=0.) lpenc_requested(i_del2ssat)=.true.
 
      lpenc_diagnos(i_ssat)=.true.
    endsubroutine pencil_criteria_supersat
!***********************************************************************
    subroutine pencil_interdep_supersat(lpencil_in)
!  Interdependency among pencils provided by the Pscalar module
!  is specified here.
      logical, dimension(npencils) :: lpencil_in
!
      lpencil_in(i_ssat)=.true.
      if (lpencil_in(i_ugssat)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_gssat)=.true.
      endif
    endsubroutine pencil_interdep_supersat
!**********************************************************************
    subroutine calc_pencils_supersat(f,p)
!
!  Calculate supersat Pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      
! ssat
      if (lpencil(i_ssat)) p%ssat=f(l1:l2,m,n,issat)
! gssat
      if (lpencil(i_gssat)) then
        call grad(f,issat,p%gssat)
      endif
! ugssat
      if (lpencil(i_ugssat)) then
        call u_dot_grad(f,issat,p%gssat,p%uu,p%ugssat,UPWIND=lupw_ssat)
      endif
! del2ssat
      if (lpencil(i_del2ssat)) then
        call del2(f,issat,p%del2ssat)
      endif
!
    endsubroutine calc_pencils_supersat
!***********************************************************************
    subroutine dssat_dt(f,df,p)
      use Diagnostics
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: diff_op,diff_op2,bump,gcgu
      real :: ssat_xyaver
      real :: tau=10., A1=5e-4
      real :: lam_gradC_fact=1., om_gradC_fact=1., gradC_fact=1.
      integer :: j, k
      integer, parameter :: nxy=nxgrid*nygrid
!
      intent(in)  :: f
      intent(out) :: df
!
      character(len=2) :: id
!  Identify module and boundary conditions.
!
      if (nosupersat) then
        if (headtt.or.ldebug) print*,'not SOLVED: dssat_dt'
      else
        if (headtt.or.ldebug) print*,'SOLVE dssat_dt'
      endif
      if (headtt) then
          write(id,'(i0)')
          call identify_bcs('ssat'//trim(id),issat)
      endif
!  Passive scalar equation.
!
        df(l1:l2,m,n,issat)=df(l1:l2,m,n,issat)-p%ugssat &
                +supersat_diff*p%del2ssat
!
!  Passive scalar sink/source.
!
!        if (lsupersat_sink) then
!          if (Rsupersat_sink==0) then
!            bump=supersat_sink
!          else
!            bump=supersat_sink*exp(-0.5*(x(l1:l2)**2+y(m)**2+z(n)**2)/Rsupersat_sink**2)
!          endif
!          df(l1:l2,m,n,issat)=df(l1:l2,m,n,issat)-spread(bump,2)*p%ssat
!        endif
! 1-June-16/XY coded: to be completed 
         if (lsupersat_sink) then
                 if (Rsupersat_sink==0) then
                    bump=-f(l1:l2,m,n,issat)/tau
                  else
                    bump=-f(l1:l2,m,n,issat)/tau+ &
                    !A1*fp(k,ivpz)
!AB: this fp doesn't exist, so I remove it for now, so it compiles
                    A1
                 endif
                 df(l1:l2,m,n,issat)=df(l1:l2,m,n,issat)-p%ugssat+bump 
         endif
!       
        if (idiag_ssatrms/=0) & 
            call sum_mn_name(p%ssat**2,idiag_ssatrms,lsqrt=.true.)
 
    endsubroutine dssat_dt
!
!***********************************************************************
    subroutine read_supersat_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=supersat_init_pars, IOSTAT=iostat)
!
    endsubroutine read_supersat_init_pars
!***********************************************************************
    subroutine write_supersat_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=supersat_init_pars)
!
    endsubroutine write_supersat_init_pars
!***********************************************************************
    subroutine read_supersat_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=supersat_run_pars, IOSTAT=iostat)
!
    endsubroutine read_supersat_run_pars
!***********************************************************************
    subroutine write_supersat_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=supersat_run_pars)
!
    endsubroutine write_supersat_run_pars
!
!***********************************************************************
    subroutine rprint_supersat(lreset,lwrite)
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
      logical :: lwr
!
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!      
      if (lreset) then
        idiag_ssatrms=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ssatrms',idiag_ssatrms)
      enddo
!
      if (lwr) then 
        write(3,*) 'issat = ', issat
      endif
    endsubroutine rprint_supersat 
endmodule Supersat
