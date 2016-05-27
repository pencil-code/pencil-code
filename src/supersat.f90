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
! CPARAM logical, parameter :: lpscalar = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED cc; cc1
! PENCILS PROVIDED gcc; ugcc
! PENCILS PROVIDED gcc2; gcc1
! PENCILS PROVIDED del2cc; del6cc
! PENCILS PROVIDED g5cc; g5ccglnrho
! PENCILS PROVIDED hcc
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
  real, dimension(3) :: gradC0=(/0.0,0.0,0.0/)
  real :: ampllncc=impossible, widthlncc=impossible, lncc_min
  real :: ampllncc2=impossible, radius_lncc=impossible
  real :: kx_lncc=impossible, ky_lncc=impossible,kz_lncc=impossible
  real :: epsilon_lncc=impossible
  real :: cc_left=1., cc_right=0.
  real :: amplcc=0.1, widthcc=0.5, cc_min=0.0
  real :: amplcc2=0.0, kx_cc=1.0, ky_cc=1.0, kz_cc=1.0, radius_cc=0.0
  real :: kxx_cc=0.0, kyy_cc=0.0, kzz_cc=0.0
  real :: epsilon_cc=0.0, cc_const=1.0
  real :: zoverh=1.0, hoverr=0.05, powerlr=3.0
  logical :: nosupersat=.false., reinitialize_cc=.false.
  logical :: reinitialize_lncc=.false.
  character (len=labellen) :: initlncc='impossible', initlncc2='impossible'
  character (len=labellen) :: initcc='nothing', initcc2='zero'
  character (len=40) :: tensor_pscalar_file
!
  namelist /pscalar_init_pars/ &
      initcc, initcc2,amplcc, amplcc2, kx_cc, ky_cc, kz_cc, radius_cc, &
      cc_left, cc_right, &
      epsilon_cc, widthcc, cc_min, cc_const, initlncc, initlncc2, ampllncc, &
      ampllncc2, kx_lncc, ky_lncc, kz_lncc, radius_lncc, epsilon_lncc, &
      widthlncc, kxx_cc, kyy_cc, kzz_cc, hoverr, powerlr, zoverh
!
!  Run parameters.
!
  real :: pscalar_diff=0.0, tensor_pscalar_diff=0.0, soret_diff=0.0
  real :: diffcc_shock = 0.
  real :: pscalar_diff_hyper3=0.0
  real :: rhoccm=0.0, cc2m=0.0, gcc2m=0.0
  real :: pscalar_sink=0.0, Rpscalar_sink=0.5
  real :: lam_gradC=0.0, om_gradC=0.0, lambda_cc=0.0
  real :: scalaracc=0.0
  real :: LLambda_cc=0.0
  logical :: lpscalar_sink=.false., lgradC_profile=.false., lreactions=.false.
!**********************************
    subroutine register_pscalar()
      use FArrayManager
!
      !lpscalar_nolog = .true.
      lpscalar = .true.
!
      call farray_register_pde('cc', icc)
      ilncc = 0                 ! needed for idl
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_pscalar
!***********************************************************************
    subroutine initialize_pscalar(f)
!
!  Perform any necessary post-parameter read initialization
!  Since the passive scalar is often used for diagnostic purposes
!  one may want to reinitialize it to its initial distribution.
!
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!
!  set to zero and then call the same initial condition
!  that was used in start.csh
!
      if (reinitialize_cc) then
        f(:,:,:,icc)=0.
        call init_lncc(f)
      endif
!
      if (lroot .and. diffcc_shock /= 0.) print*, 'initialize_pscalar: shock diffusion, diffcc_shock = ', diffcc_shock
!
      if (lnotpassive) scalaracc=3./5./hoverr**2
!
    endsubroutine initialize_pscalar
!
!**********************************************************************
    subroutine calc_pencils_pscalar(f,p)
!
!  Calculate pscalar Pencils.
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
      real, dimension(nx) :: dot2_tmp
      integer :: i
! cc
      if (lpencil(i_cc)) p%cc=f(l1:l2,m,n,icc)
! cc1
      if (lpencil(i_cc1)) p%cc1=1./p%cc
! gcc
      if (lpencil(i_gcc)) then
          call grad(f,icc+i-1,p%gcc)
      endif
! ugcc
      if (lpencil(i_ugcc)) then
          call u_dot_grad(f,icc,p%gcc,p%uu,p%ugcc,UPWIND=lupw_cc)
      endif
!
    endsubroutine calc_pencils_pscalar
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
      real :: lam_gradC_fact=1., om_gradC_fact=1., gradC_fact=1.
      integer :: j, k
      integer, parameter :: nxy=nxgrid*nygrid
!
      intent(in)  :: f
      intent(out) :: df
!
      character(len=2) :: id
!  Passive scalar equation.
!
        df(l1:l2,m,n)=df(l1:l2,m,n)-p%ugcc
!
!  Passive scalar sink.
!
        if (pscalar_sink) then
          if (Rpscalar_sink==0) then
            bump=pscalar_sink
          else
            bump=pscalar_sink*exp(-0.5*(x(l1:l2)**2+y(m)**2+z(n)**2)/Rpscalar_sink**2)
          endif
          df(l1:l2,m,n)=df(l1:l2,m,n)-spread(bump,2)*p%cc
        endif
    endsubroutine dlncc_dt
!
!***********************************************************************
    subroutine read_pscalar_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=pscalar_init_pars, IOSTAT=iostat)
!
    endsubroutine read_pscalar_init_pars
!***********************************************************************
    subroutine write_pscalar_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=pscalar_init_pars)
!
    endsubroutine write_pscalar_init_pars
!***********************************************************************
    subroutine read_pscalar_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=pscalar_run_pars, IOSTAT=iostat)
!
    endsubroutine read_pscalar_run_pars
!***********************************************************************
    subroutine write_pscalar_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=pscalar_run_pars)
!
    endsubroutine write_pscalar_run_pars
!
endmodule Supersat
