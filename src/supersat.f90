! $Id$
!
! This module is used to solve the passive scalar 
! equation of supersaturation for swarm model. 
!
! **************************************************
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
! PENCILS PROVIDED cc; cc1
! PENCILS PROVIDED gcc(3); ugcc
! PENCILS PROVIDED del2cc
!***************************************************************
module Superstat
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
  logical :: nosupersat=.false., reinitialize_cc=.false.
  logical :: reinitialize_lncc=.false.
  character (len=labellen) :: initlncc='impossible'
 !
  namelist /supersat_init_pars/ &
           initlncc
!
!  Run parameters.
!
  real :: supersat_diff=0.0
  real :: supersat_sink=0.0, Rsupersat_sink=0.5
  real, dimension(3) :: gradC0=(/0.0,0.0,0.0/)
  logical :: lsupersat_sink=.false.

  namelist /supersat_run_pars/ &
      lsupersat_sink, Rsupersat_sink, supersat_sink &
      supersat_diff, gradC0
!
!Diagnostics variables
integer :: idiag_ccrms=0
!**********************************
    subroutine register_supersat()
      use FArrayManager
!
      lsupersat = .true.
!
      call farray_register_pde('cc', icc)
      ilncc = 0                 ! needed for idl
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
      if (reinitialize_cc) then
        f(:,:,:,icc)=0.
        call init_lncc(f)
      endif
!
    endsubroutine initialize_supersat
!
!***********************************************************************
    subroutine init_lncc(f)
!  initialise passive scalar field; called from start.f90
      use Sub
      use Initcond
      use InitialCondition, only: initial_condition_lncc
!
      real, dimension (mx,my,mz,mfarray) :: f
       select case (initcc)
        case ('nothing')
        case ('zero'); f(:,:,:,icc)=0.0
        case ('constant'); f(:,:,:,icc)=cc_const
       endselect
    endsubroutine init_lncc
   
!***********************************************************************
    subroutine pencil_criteria_supersat()
      lpenc_requested(i_cc)=.true.
            
      if (lsupersat_sink) lpenc_requested(i_cc)=.true.
      if (supersat_diff/=0.) lpenc_requested(i_del2cc)=.true.
 
      lpenc_diagnos(i_cc)=.true.
    endsubroutine pencil_criteria_supersat
!***********************************************************************
    subroutine pencil_interdep_supersat(lpencil_in)
!  Interdependency among pencils provided by the Pscalar module
!  is specified here.
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_cc1)) lpencil_in(i_cc)=.true.
      if (lpencil_in(i_ugcc)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_gcc)=.true.
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
      
! cc
      if (lpencil(i_cc)) p%cc=f(l1:l2,m,n,icc)
! cc1
      if (lpencil(i_cc1)) p%cc1=1./p%cc
! gcc
      if (lpencil(i_gcc)) then
        call grad(f,icc,p%gcc)
      endif
! ugcc
      if (lpencil(i_ugcc)) then
        call u_dot_grad(f,icc,p%gcc,p%uu,p%ugcc,UPWIND=lupw_cc)
      endif
! del2cc
      if (lpencil(i_del2cc)) then
          call del2(f,icc+i-1,p%del2cc(:,i))
      endif
!
    endsubroutine calc_pencils_supersat
!***********************************************************************
    subroutine dlncc_dt(f,df,p)
      use Diagnostics
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: diff_op,diff_op2,bump,gcgu
      real :: cc_xyaver
      real :: tau=10., A1=5.*e-4
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
        if (headtt.or.ldebug) print*,'not SOLVED: dlncc_dt'
      else
        if (headtt.or.ldebug) print*,'SOLVE dlncc_dt'
      endif
      if (headtt) then
          write(id,'(i0)')
          call identify_bcs('cc'//trim(id),icc)
      endif
!  Passive scalar equation.
!
        df(l1:l2,m,n,icc)=df(l1:l2,m,n,icc)-p%ugcc &
                +supersat_diff*p%del2cc
!
!  Passive scalar sink/source.
!
!        if (supersat_sink) then
!          if (Rsupersat_sink==0) then
!            bump=supersat_sink
!          else
!            bump=supersat_sink*exp(-0.5*(x(l1:l2)**2+y(m)**2+z(n)**2)/Rsupersat_sink**2)
!          endif
!          df(l1:l2,m,n,icc)=df(l1:l2,m,n,icc)-spread(bump,2)*p%cc
!        endif
! 1-June-16/XY coded: to be completed 
         if (supersat_sink) then
                print*,"XY" 
                 if (Rsupersat_sink==0) then
                    bump=-f(l1:l2,m,n,icc)/tau
                  else
                    bump=-f(l1:l2,m,n,icc)/tau+ &
                    A1*fp(k,ivpz)
                 endif
                 df(l1:l2,m,n,icc)=df(l1:l2,m,n,icc)-p%ugcc+bump 
         endif
!       
        if (idiag_ccrms/=0) & 
            call sum_mn_name(p%cc(:,1)**2,idiag_ccrms,lsqrt=.true.)
 
    endsubroutine dlncc_dt
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
        idiag_ccrms=0,idiag_uzcmz=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ccrms',idiag_ccrms)
      enddo
!
      if (lwr) then 
        write(3,*) 'ilncc=0'
        write(3,*) 'icc = ', icc
      endif
    endsubroutine rprint_supersat 
endmodule Supersat
