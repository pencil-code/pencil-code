! $Id$

!  Cosmic Ray Flux
!

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!
! CPARAM logical, parameter :: lcosmicrayflux = .true.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
!
!***************************************************************

module Cosmicrayflux

  use Cparam
  use Messages

  implicit none

  include 'cosmicrayflux.h'

  character (len=labellen) :: initfcr='zero'
  ! input parameters
  real :: amplfcr=0.
  namelist /cosmicrayflux_init_pars/ &
       amplfcr
  ! run parameters
  real :: tau=0.,kpara=0.,kperp=0.,tau1=0.
  namelist /cosmicrayflux_run_pars/ &
       tau, kpara, kperp

  ! other variables (needs to be consistent with reset list below)
!  integer :: idiag_b2m=0,idiag_bm2=0,idiag_j2m=0,idiag_jm2=0

  contains

!***********************************************************************
    subroutine register_cosmicrayflux()
!
!  Initialise variables which should know that we solve for the vector
!  potential: ifcr, etc; increase nvar accordingly
!
!  1-may-02/wolf: coded
!
      use Cdata
      use FArrayManager
!
      call farray_register_pde('fcr',ifcr,vector=3)
      ifcrx = ifcr; ifcry = ifcr+1; ifcrz = ifcr+2
!
!  Identify version number.
!
      if (lroot) call cvs_id( &
          "$Id$")
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',fcr $'
          if (nvar == mvar) write(4,*) ',fcr'
        else
          write(4,*) ',fcr $'
        endif
        write(15,*) 'fcr = fltarr(mx,my,mz,3)*one'
      endif
!
    endsubroutine register_cosmicrayflux
!***********************************************************************
    subroutine initialize_cosmicrayflux(f)
!
!  Perform any post-parameter-read initialization
!
!  24-nov-02/tony: dummy routine - nothing to do at present
!  20-may-03/axel: reinitalize_aa added
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f

      if (tau /= 0.) tau1=1./tau

!
    endsubroutine initialize_cosmicrayflux
!***********************************************************************
    subroutine init_fcr(f)
!
!  initialise magnetic field; called from start.f90
!  AB: maybe we should here call different routines (such as rings)
!  AB: and others, instead of accummulating all this in a huge routine.
!  We have an init parameter (initaa) to stear magnetic i.c. independently.
!
!   7-nov-2001/wolf: coded
!
      use Cdata
      use Mpicomm
!      use EquationOfState
!      use Gravity, only: gravz
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mfarray) :: f
!,tmp,prof
!      real, dimension (nx,3) :: bb
!      real, dimension (nx) :: b2,fact
!      real :: beq2
!
      select case(initfcr)

      case('zero', '0'); f(:,:,:,ifcrx:ifcrz) = 0.
! probably no more cases needed for fcr
      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'init_fcr: No such such value for initfcr: ', trim(initfcr)
        call stop_it(" ")

      endselect
!
!
    endsubroutine init_fcr
!***********************************************************************
    subroutine pencil_criteria_cosmicrayflux()
!
!   All pencils that the Magnetic module depends on are specified here.
!
!  19-11-04/anders: coded
!
      use Cdata
!
      lpenc_requested(i_gecr)=.true.
      lpenc_requested(i_bb)=.true.

!      lpenc_requested(i_fcr)=.true.
!
    endsubroutine pencil_criteria_cosmicrayflux
!***********************************************************************
    subroutine pencil_interdep_cosmicrayflux(lpencil_in)
!
!
      logical, dimension(npencils) :: lpencil_in
!
!
    endsubroutine pencil_interdep_cosmicrayflux
!***********************************************************************
    subroutine calc_pencils_cosmicrayflux(f,p)
!
!  Calculate Cosmicray Flux pencils - to be done
!
!
      use Cdata
      use Sub
      use Deriv
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
!      real, dimension (nx,3) :: bb_ext,bb_ext_pot,ee_ext,jj_ext
!      real, dimension (nx) :: rho1_jxb,quenching_factor,alpha_total
!      real :: B2_ext,c,s
!      integer :: i,j
!
      intent(in)  :: f
      intent(inout) :: p
! fcr
!      if (lpencil(i_fcr)) p%fcr=f(l1:l2,m,n,ifcrx:ifcrz)
!
    endsubroutine calc_pencils_cosmicrayflux
!***********************************************************************
    subroutine dfcr_dt(f,df,p)
!
!  Cosmicray Flux evolution
!
!  08-mar-05/snod: adapted from daa_dt
!
      use Cdata
      use Sub
      use Slices
      use IO, only: output_pencil
      use Mpicomm, only: stop_it
!      use EquationOfState, only: eoscalc,gamma1
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: BuiBujgecr, bunit
      real, dimension (nx)   :: b2, b21
      real, dimension (nx)   :: tmp
      integer :: i,j
      type (pencil_case) :: p
!
!      real, dimension (nx) :: uxb_dotB0,oxuxb_dotB0,jxbxb_dotB0,uxDxuxb_dotB0
!      real, dimension (nx) :: gpxb_dotB0,uxj_dotB0,hall_ueff2
!      real, dimension (nx) :: b2b13,sign_jo
!      real, dimension (nx) :: eta_mn,eta_tot
!      real, dimension (nx) :: eta_smag,etatotal,fres2
!      integer :: i,j
!
      intent(in)     :: f
      intent(inout)  :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dfcr_dt: SOLVE'
      if (headtt) then
        call identify_bcs('Fecx',ifcrx)
        call identify_bcs('Fecy',ifcry)
        call identify_bcs('Fecz',ifcrz)
      endif
!
!
      call dot2_mn(p%bb,b2)
!  with frequency omega_Bz_ext
      b21=1./max(tini,b2)
      call multsv_mn(sqrt(b21),p%bb,bunit)
!
!

      do i=1,3
        tmp=0.
        do j=1,3
          tmp=tmp+bunit(:,i)*bunit(:,j)*p%gecr(:,j)
        enddo
        BuiBujgecr(:,i)=tmp
      enddo
!
!  Cosmic Ray Flux equation
!
      df(l1:l2,m,n,ifcrx:ifcrz) = df(l1:l2,m,n,ifcrx:ifcrz) &
      - tau1*f(l1:l2,m,n,ifcrx:ifcrz)                       &
      - kperp*p%gecr                                        &
      - (kpara - kperp)*BuiBujgecr
!
!
!  Calculate diagnostic quantities
!
      if (ldiagnos) then
!
!  cosmicrayflux components at one point (=pt)
!
!        if (lroot.and.m==mpoint.and.n==npoint) then
!          if (idiag_fcrxpt/=0) call save_name(p%fcr(lpoint-nghost,1),idiag_fcrxpt)
!          if (idiag_fcrypt/=0) call save_name(p%fcr(lpoint-nghost,2),idiag_fcrypt)
!         if (idiag_fcrzpt/=0) call save_name(p%fcr(lpoint-nghost,3),idiag_fcrzpt)
!        endif
!
!  v_A = |B|/sqrt(rho); in units where "4pi"=1
!
!
      endif ! endif (ldiagnos)
!
!  debug output
!
!
!  write B-slices for output in wvid in run.f90
!  Note: ix is the index with respect to array with ghost zones.
!
!      if (lvideo.and.lfirst) then
!        do j=1,3
!          bb_yz(m-m1+1,n-n1+1,j)=p%bb(ix_loc-l1+1,j)
!          if (m==iy_loc)  bb_xz(:,n-n1+1,j)=p%bb(:,j)
!          if (n==iz_loc)  bb_xy(:,m-m1+1,j)=p%bb(:,j)
!          if (n==iz2_loc) bb_xy2(:,m-m1+1,j)=p%bb(:,j)
!        enddo
!        b2_yz(m-m1+1,n-n1+1)=p%b2(ix_loc-l1+1)
!        if (m==iy_loc)  b2_xz(:,n-n1+1)=p%b2
!        if (n==iz_loc)  b2_xy(:,m-m1+1)=p%b2
!        if (n==iz2_loc) b2_xy2(:,m-m1+1)=p%b2
!        jb_yz(m-m1+1,n-n1+1)=p%jb(ix_loc-l1+1)
!        if (m==iy_loc)  jb_xz(:,n-n1+1)=p%jb
!        if (n==iz_loc)  jb_xy(:,m-m1+1)=p%jb
!        if (n==iz2_loc) jb_xy2(:,m-m1+1)=p%jb
!        if (bthresh_per_brms/=0) call calc_bthresh
!        call vecout(41,trim(directory)//'/bvec',p%bb,bthresh,nbvec)
!      endif
!
    endsubroutine dfcr_dt
!***********************************************************************
    subroutine read_cosmicrayflux_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=cosmicrayflux_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=cosmicrayflux_init_pars,ERR=99)
      endif

99    return
    endsubroutine read_cosmicrayflux_init_pars
!***********************************************************************
    subroutine write_cosmicrayflux_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=cosmicrayflux_init_pars)

    endsubroutine write_cosmicrayflux_init_pars
!***********************************************************************
    subroutine read_cosmicrayflux_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=cosmicrayflux_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=cosmicrayflux_run_pars,ERR=99)
      endif

99    return
    endsubroutine read_cosmicrayflux_run_pars
!***********************************************************************
    subroutine write_cosmicrayflux_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=cosmicrayflux_run_pars)

    endsubroutine write_cosmicrayflux_run_pars
!***********************************************************************
    subroutine rprint_cosmicrayflux(lreset,lwrite)
!
!  reads and registers print parameters relevant for cosmicrayflux
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Cdata
      use Sub
!
      integer :: iname,inamez,ixy,irz
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of RELOAD
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
!        idiag_b2m=0; idiag_bm2=0; idiag_j2m=0; idiag_jm2=0; idiag_abm=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
!        call parse_name(iname,cname(iname),cform(iname),'dteta',idiag_dteta)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
!        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bxmz',idiag_bxmz)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do ixy=1,nnamexy
!        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bxmxy',idiag_bxmxy)
      enddo
!
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
!        call parse_name(irz,cnamerz(irz),cformrz(irz),'brmphi'  ,idiag_brmphi)
      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!
      if (lwr) then
        write(3,*) 'ifcr=',ifcr
      endif
!
    endsubroutine rprint_cosmicrayflux
!***********************************************************************

endmodule Cosmicrayflux
