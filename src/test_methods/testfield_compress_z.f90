! $Id$

!  This modules deals with all aspects of testfield fields; if no
!  testfield fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testfield relevant subroutines listed in here.

!  Note: this routine requires that MVAR and MAUX contributions
!  together with njtest are set correctly in the cparam.local file.
!  njtest must be set at the end of the file such that 6*njtest=MVAR.
!
!  Example:
!  ! MVAR CONTRIBUTION 35
!  ! MAUX CONTRIBUTION 35
!  integer, parameter :: njtest=5

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ltestfield = .true.
! CPARAM logical, parameter :: ltestfield_z = .true.
! CPARAM logical, parameter :: ltestfield_xy = .false.
! CPARAM logical, parameter :: ltestfield_xz = .false.
!
!***************************************************************

module Testfield

  use Cparam
  use Cdata
  use Messages

  implicit none

  include '../testfield.h'
!
! Slice precalculation buffers
!
  real, target, dimension(:,:,:), allocatable :: bb11_xy, bb11_xy2, bb11_xy3, bb11_xy4
  real, target, dimension(:,:,:), allocatable :: bb11_xz, bb11_xz2, bb11_yz
!
!  cosine and sine function for setting test fields and analysis
!
  real, dimension(mz) :: cz,sz,c2z,csz,s2z,c2kz,s2kz
  real :: phase_testfield=.0,cs0test=1.,cs20test=1.
!
  character (len=labellen), dimension(ninit) :: initaatest='nothing'
  character (len=labellen), dimension(ninit) :: inituutest='nothing'
  character (len=labellen), dimension(ninit) :: inithhtest='nothing'
  real, dimension (ninit) :: kx_aatest=1.,ky_aatest=1.,kz_aatest=1.
  real, dimension (ninit) :: phasex_aatest=0.,phasez_aatest=0.
  real, dimension (ninit) :: amplaatest=0.,ampluutest=0.,amplhhtest=0.
  integer :: iuxtest=0,iuytest=0,iuztest=0
  integer :: iu0xtest=0,iu0ztest=0,ihxtest=0,ihhtest=0
  integer :: iE0=0

  ! input parameters
  real, dimension(3) :: Btest_ext=(/0.,0.,0./)
  real :: taainit=0.,daainit=0.,taainit_previous=0.
  logical :: reinitialize_aatest=.false.
  logical :: reinitialize_from_mainrun=.false.
  logical :: zextent=.true.,lsoca=.false.,lset_bbtest2=.false.
  logical :: luxb_as_aux=.true.,ljxb_as_aux=.true.
  logical :: lugu_as_aux=.true.,lugh_as_aux=.true., lSgh_as_aux=.true.
  logical :: linit_aatest=.false.
  logical :: lignore_uxbtestm=.false., lignore_jxbtestm=.false.
  logical :: lignore_ugutestm=.false., lignore_ughtestm=.false., lignore_Sghtestm=.false.
  logical :: lugutest=.true., lfprestest=.true., lSghtest=.true.
  logical :: lphase_adjust=.false.
  logical :: luse_main_run=.true., lvisc_simplified_testfield=.false.
  logical :: lremove_meanaa0x_test=.false., lremove_meanaa0y_test=.false., &
             lremove_meanuu0x_test=.false., lremove_meanuu0y_test=.false., &
             lzero_only=.false., lremove_E0=.false., lremove_F0=.false., lremove_Q0=.false.
  character (len=labellen) :: itestfield='B11-B21',itestfield_method='(i)'
  real :: ktestfield=1., ktestfield1=1.
  real :: lin_testfield=0.,lam_testfield=0.,om_testfield=0.,delta_testfield=0.
  integer :: naainit
  real :: bamp=1.,bamp1=1.,bamp12=1.
  real :: damp_uxb=0.

  namelist /testfield_init_pars/ &
       Btest_ext,zextent,initaatest, &
       initaatest,inituutest,inithhtest, &
       amplaatest,ampluutest,amplhhtest, &
       kx_aatest,ky_aatest,kz_aatest, &
       cs0test,phasex_aatest,phasez_aatest, &
       luxb_as_aux,ljxb_as_aux,lugu_as_aux,lugh_as_aux,lSgh_as_aux, &
       luse_main_run, &
       lugutest, lfprestest, lSghtest

  ! run parameters
  real :: etatest=0.,etatest1=0.,nutest=0.,nutest1=0.,eta_hyper3_test=0.
  real :: ampl_fcont_aatest=1.,ampl_fcont_uutest=1.
  real, dimension(njtest) :: rescale_aatest=0.,rescale_uutest=0.,rescale_hhtest=0.
  logical :: ltestfield_newz=.true.,leta_rank2=.true.
  logical :: lforcing_cont_aatest=.false.,lforcing_cont_uutest=.false.
  logical :: lupw_uutest=.false., lupw_hhtest=.false.
  real :: rho0test=1., rho0test1

  namelist /testfield_run_pars/ &
       reinitialize_aatest,reinitialize_from_mainrun, &
       Btest_ext,zextent,lsoca, &
       lset_bbtest2,itestfield,ktestfield,itestfield_method, &
       etatest,etatest1,nutest,nutest1, &
       lin_testfield,lam_testfield,om_testfield,delta_testfield, &
       ltestfield_newz,leta_rank2,lphase_adjust,phase_testfield, &
       cs0test, &
       luxb_as_aux,ljxb_as_aux,lugu_as_aux,lugh_as_aux,lSgh_as_aux, &
       lignore_uxbtestm,lignore_jxbtestm, &
       lignore_ugutestm,lignore_ughtestm, lignore_Sghtestm, &
       lforcing_cont_aatest,ampl_fcont_aatest, &
       lforcing_cont_uutest,ampl_fcont_uutest, &
       daainit,linit_aatest,bamp, &
       rescale_aatest,rescale_uutest, rescale_hhtest, rho0test, &
       lupw_uutest, lupw_hhtest, luse_main_run, lvisc_simplified_testfield, Omega, &
       lugutest, lfprestest, lSghtest, &
       lremove_meanaa0x_test, lremove_meanaa0y_test, &
       lremove_meanuu0x_test, lremove_meanuu0y_test, &
       damp_uxb, lzero_only, lremove_E0, lremove_F0, lremove_Q0

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_alp11=0      ! DIAG_DOC: $\alpha_{11}$
  integer :: idiag_alp21=0      ! DIAG_DOC: $\alpha_{21}$
  integer :: idiag_alp31=0      ! DIAG_DOC: $\alpha_{31}$
  integer :: idiag_alp12=0      ! DIAG_DOC: $\alpha_{12}$
  integer :: idiag_alp22=0      ! DIAG_DOC: $\alpha_{22}$
  integer :: idiag_alp32=0      ! DIAG_DOC: $\alpha_{32}$
  integer :: idiag_eta11=0      ! DIAG_DOC: $\eta_{11}k$
  integer :: idiag_eta21=0      ! DIAG_DOC: $\eta_{21}k$
  integer :: idiag_eta12=0      ! DIAG_DOC: $\eta_{12}k$
  integer :: idiag_eta22=0      ! DIAG_DOC: $\eta_{22}k$
  integer :: idiag_alpK=0       ! DIAG_DOC: $\alpha^K$
  integer :: idiag_alpM=0       ! DIAG_DOC: $\alpha^M$
  integer :: idiag_alpMK=0      ! DIAG_DOC: $\alpha^{MK}$
  integer :: idiag_phi11=0      ! DIAG_DOC: $\phi_{11}$
  integer :: idiag_phi21=0      ! DIAG_DOC: $\phi_{21}$
  integer :: idiag_phi12=0      ! DIAG_DOC: $\phi_{12}$
  integer :: idiag_phi22=0      ! DIAG_DOC: $\phi_{22}$
  integer :: idiag_phi32=0      ! DIAG_DOC: $\phi_{32}$
  integer :: idiag_psi11=0      ! DIAG_DOC: $\psi_{11}k$
  integer :: idiag_psi21=0      ! DIAG_DOC: $\psi_{21}k$
  integer :: idiag_psi12=0      ! DIAG_DOC: $\psi_{12}k$
  integer :: idiag_psi22=0      ! DIAG_DOC: $\psi_{22}k$
  integer :: idiag_sig1=0       ! DIAG_DOC: $\sigma_1$
  integer :: idiag_sig2=0       ! DIAG_DOC: $\sigma_2$
  integer :: idiag_sig3=0       ! DIAG_DOC: $\sigma_3$
  integer :: idiag_tau1=0       ! DIAG_DOC: $\tau_1$
  integer :: idiag_tau2=0       ! DIAG_DOC: $\tau_2$
  integer :: idiag_phiK=0       ! DIAG_DOC: $\phi^K$
  integer :: idiag_phiM=0       ! DIAG_DOC: $\phi^M$
  integer :: idiag_phiMK=0      ! DIAG_DOC: $\phi^{MK}$
  integer :: idiag_alp11cc=0    ! DIAG_DOC: $\alpha_{11}\cos^2 kz$
  integer :: idiag_alp21sc=0    ! DIAG_DOC: $\alpha_{21}\sin kz\cos kz$
  integer :: idiag_alp12cs=0    ! DIAG_DOC: $\alpha_{12}\cos kz\sin kz$
  integer :: idiag_alp22ss=0    ! DIAG_DOC: $\alpha_{22}\sin^2 kz$
  integer :: idiag_eta11cc=0    ! DIAG_DOC: $\eta_{11}\cos^2 kz$
  integer :: idiag_eta21sc=0    ! DIAG_DOC: $\eta_{21}\sin kz\cos kz$
  integer :: idiag_eta12cs=0    ! DIAG_DOC: $\eta_{12}\cos kz\sin kz$
  integer :: idiag_eta22ss=0    ! DIAG_DOC: $\eta_{22}\sin^2 kz$
  integer :: idiag_s2kzDFm=0    ! DIAG_DOC: $\left<\sin2kz\nabla\cdot F\right>$
  integer :: idiag_M11=0        ! DIAG_DOC: ${\cal M}_{11}$
  integer :: idiag_M22=0        ! DIAG_DOC: ${\cal M}_{22}$
  integer :: idiag_M33=0        ! DIAG_DOC: ${\cal M}_{33}$
  integer :: idiag_M11cc=0      ! DIAG_DOC: ${\cal M}_{11}\cos^2 kz$
  integer :: idiag_M11ss=0      ! DIAG_DOC: ${\cal M}_{11}\sin^2 kz$
  integer :: idiag_M22cc=0      ! DIAG_DOC: ${\cal M}_{22}\cos^2 kz$
  integer :: idiag_M22ss=0      ! DIAG_DOC: ${\cal M}_{22}\sin^2 kz$
  integer :: idiag_M12cs=0      ! DIAG_DOC: ${\cal M}_{12}\cos kz\sin kz$
  integer :: idiag_bx11pt=0     ! DIAG_DOC: $b_x^{11}$
  integer :: idiag_bx21pt=0     ! DIAG_DOC: $b_x^{21}$
  integer :: idiag_bx12pt=0     ! DIAG_DOC: $b_x^{12}$
  integer :: idiag_bx22pt=0     ! DIAG_DOC: $b_x^{22}$
  integer :: idiag_bx0pt=0      ! DIAG_DOC: $b_x^{0}$
  integer :: idiag_by11pt=0     ! DIAG_DOC: $b_y^{11}$
  integer :: idiag_by21pt=0     ! DIAG_DOC: $b_y^{21}$
  integer :: idiag_by12pt=0     ! DIAG_DOC: $b_y^{12}$
  integer :: idiag_by22pt=0     ! DIAG_DOC: $b_y^{22}$
  integer :: idiag_by0pt=0      ! DIAG_DOC: $b_y^{0}$
  integer :: idiag_u11rms=0     ! DIAG_DOC: $\left<u_{11}^2\right>^{1/2}$
  integer :: idiag_u21rms=0     ! DIAG_DOC: $\left<u_{21}^2\right>^{1/2}$
  integer :: idiag_u12rms=0     ! DIAG_DOC: $\left<u_{12}^2\right>^{1/2}$
  integer :: idiag_u22rms=0     ! DIAG_DOC: $\left<u_{22}^2\right>^{1/2}$
  integer :: idiag_h11rms=0     ! DIAG_DOC: $\left<h_{11}^2\right>^{1/2}$
  integer :: idiag_h21rms=0     ! DIAG_DOC: $\left<h_{21}^2\right>^{1/2}$
  integer :: idiag_h12rms=0     ! DIAG_DOC: $\left<h_{12}^2\right>^{1/2}$
  integer :: idiag_h22rms=0     ! DIAG_DOC: $\left<h_{22}^2\right>^{1/2}$
  integer :: idiag_j11rms=0     ! DIAG_DOC: $\left<j_{11}^2\right>^{1/2}$
  integer :: idiag_b11rms=0     ! DIAG_DOC: $\left<b_{11}^2\right>^{1/2}$
  integer :: idiag_b21rms=0     ! DIAG_DOC: $\left<b_{21}^2\right>^{1/2}$
  integer :: idiag_b12rms=0     ! DIAG_DOC: $\left<b_{12}^2\right>^{1/2}$
  integer :: idiag_b22rms=0     ! DIAG_DOC: $\left<b_{22}^2\right>^{1/2}$
  integer :: idiag_ux0m=0       ! DIAG_DOC: $\left<u_{0_x}\right>$
  integer :: idiag_uy0m=0       ! DIAG_DOC: $\left<u_{0_y}\right>$
  integer :: idiag_ux11m=0      ! DIAG_DOC: $\left<u_{11_x}\right>$
  integer :: idiag_uy11m=0      ! DIAG_DOC: $\left<u_{11_y}\right>$
  integer :: idiag_u0rms=0      ! DIAG_DOC: $\left<u_{0}^2\right>^{1/2}$
  integer :: idiag_b0rms=0      ! DIAG_DOC: $\left<b_{0}^2\right>^{1/2}$
  integer :: idiag_h0rms=0      ! DIAG_DOC: $\left<h_{0}^2\right>^{1/2}$
  integer :: idiag_rho0m=0      ! DIAG_DOC: $\left<\exp h_{0}\right>$
  integer :: idiag_u0max=0      ! DIAG_DOC: $\operatorname{max}\left|\boldsymbol{u}_{0}\right|$
  integer :: idiag_b0max=0      ! DIAG_DOC: $\operatorname{max}\left|\boldsymbol{b}_{0}\right|$
  integer :: idiag_h0max=0      ! DIAG_DOC: $\operatorname{max}h_{0}$
  integer :: idiag_bhrms=0      ! DIAG_DOC: $\left<b_{h}^2\right>^{1/2}$
  integer :: idiag_jb0m=0       ! DIAG_DOC: $\left<jb_{0}\right>$
  integer :: idiag_E11rms=0     ! DIAG_DOC: $\left<{\cal E}_{11}^2\right>^{1/2}$
  integer :: idiag_E21rms=0     ! DIAG_DOC: $\left<{\cal E}_{21}^2\right>^{1/2}$
  integer :: idiag_E12rms=0     ! DIAG_DOC: $\left<{\cal E}_{12}^2\right>^{1/2}$
  integer :: idiag_E22rms=0     ! DIAG_DOC: $\left<{\cal E}_{22}^2\right>^{1/2}$
  integer :: idiag_E0rms=0      ! DIAG_DOC: $\left<{\cal E}_{0}^2\right>^{1/2}$
  integer :: idiag_E0mrms=0     ! DIAG_DOC: $\left<\left<{\cal E}_{0}\right>^2\right>^{1/2}$
  integer :: idiag_Ex11pt=0     ! DIAG_DOC: ${\cal E}_x^{11}$
  integer :: idiag_Ex21pt=0     ! DIAG_DOC: ${\cal E}_x^{21}$
  integer :: idiag_Ex12pt=0     ! DIAG_DOC: ${\cal E}_x^{12}$
  integer :: idiag_Ex22pt=0     ! DIAG_DOC: ${\cal E}_x^{22}$
  integer :: idiag_Ex0pt=0      ! DIAG_DOC: ${\cal E}_x^{0}$
  integer :: idiag_Ey11pt=0     ! DIAG_DOC: ${\cal E}_y^{11}$
  integer :: idiag_Ey21pt=0     ! DIAG_DOC: ${\cal E}_y^{21}$
  integer :: idiag_Ey12pt=0     ! DIAG_DOC: ${\cal E}_y^{12}$
  integer :: idiag_Ey22pt=0     ! DIAG_DOC: ${\cal E}_y^{22}$
  integer :: idiag_Ey0pt=0      ! DIAG_DOC: ${\cal E}_y^{0}$
  integer :: idiag_bamp=0       ! DIAG_DOC: bamp
  integer :: idiag_E111z=0      ! DIAG_DOC: ${\cal E}_1^{11}$
  integer :: idiag_E211z=0      ! DIAG_DOC: ${\cal E}_2^{11}$
  integer :: idiag_E311z=0      ! DIAG_DOC: ${\cal E}_3^{11}$
  integer :: idiag_E121z=0      ! DIAG_DOC: ${\cal E}_1^{21}$
  integer :: idiag_E221z=0      ! DIAG_DOC: ${\cal E}_2^{21}$
  integer :: idiag_E321z=0      ! DIAG_DOC: ${\cal E}_3^{21}$
  integer :: idiag_E112z=0      ! DIAG_DOC: ${\cal E}_1^{12}$
  integer :: idiag_E212z=0      ! DIAG_DOC: ${\cal E}_2^{12}$
  integer :: idiag_E312z=0      ! DIAG_DOC: ${\cal E}_3^{12}$
  integer :: idiag_E122z=0      ! DIAG_DOC: ${\cal E}_1^{22}$
  integer :: idiag_E222z=0      ! DIAG_DOC: ${\cal E}_2^{22}$
  integer :: idiag_E322z=0      ! DIAG_DOC: ${\cal E}_3^{22}$
  integer :: idiag_E10z=0       ! DIAG_DOC: ${\cal E}_1^{0}$
  integer :: idiag_E20z=0       ! DIAG_DOC: ${\cal E}_2^{0}$
  integer :: idiag_E30z=0       ! DIAG_DOC: ${\cal E}_3^{0}$
  integer :: idiag_EBpq=0       ! DIAG_DOC: ${\cal E}\cdot\Bv^{pq}$
  integer :: idiag_E0Um=0       ! DIAG_DOC: ${\cal E}^0\cdot\Uv$
  integer :: idiag_E0Wm=0       ! DIAG_DOC: ${\cal E}^0\cdot\Wv$
  integer :: idiag_bx0mz=0      ! DIAG_DOC: $\left<b_{x}\right>_{xy}$
  integer :: idiag_by0mz=0      ! DIAG_DOC: $\left<b_{y}\right>_{xy}$
  integer :: idiag_bz0mz=0      ! DIAG_DOC: $\left<b_{z}\right>_{xy}$
  integer :: idiag_M11z=0       ! DIAG_DOC: $\left<{\cal M}_{11}\right>_{xy}$
  integer :: idiag_M22z=0       ! DIAG_DOC: $\left<{\cal M}_{22}\right>_{xy}$
  integer :: idiag_M33z=0       ! DIAG_DOC: $\left<{\cal M}_{33}\right>_{xy}$
  integer :: ivid_bb11=0
!
!  arrays for horizontally averaged uxb and jxb
!
  real, dimension (nz,3,njtest) :: uxbtestmz,jxbtestmz,ugutestmz,Sghtestmz
  real, dimension (nz,  njtest) :: ughtestmz
!
!  auxiliaries
!
  logical :: lBext
  integer :: jtest_start

  contains

!***********************************************************************
    subroutine register_testfield
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaatest, etc; increase nvar accordingly
!
!   3-jun-05/axel: adapted from register_magnetic
!
      use Cdata
      use FArrayManager
      use Mpicomm
      use Sub
!
!  Register test field.
!
      ntestfield = 3*njtest
      ntestflow = 3*njtest
      ntestlnrho = njtest
      call farray_register_pde('aatest',iaatest,vector=3,array=-njtest)
      call farray_index_append('ntestfield',ntestfield)
      call farray_register_pde('uutest',iuutest,vector=3,array=-njtest)
      call farray_index_append('ntestflow',ntestflow)
      call farray_register_pde('lnrhotest',ihhtest,array=-njtest)
      call farray_index_append('ntestlnrho',njtest)
!
!  Set first and last index of test solution.
!  These values are used in this form in start, but later overwritten.
!
      iaztestpq=iaatest+ntestfield-1
      iuztestpq=iuutest+ntestflow-1
      ihhtestpq=ihhtest+ntestlnrho-1
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',aatest,uutest,lnrhotest $'
          if (nvar == mvar) write(4,*) ',aatest,uutest,lnrhotest'
        else
          write(4,*) ',aatest,uutest,lnrhotest $'
        endif
        write(15,*) 'aatest = fltarr(mx,my,mz,ntestfield)*one'
        write(15,*) 'uutest = fltarr(mx,my,mz,ntestflow)*one'
        write(15,*) 'lnrhotest = fltarr(mx,my,mz,ntestlnrho)*one'
      endif
!
    endsubroutine register_testfield
!***********************************************************************
    subroutine initialize_testfield(f)
!
!  Perform any post-parameter-read initialization
!
!   2-jun-05/axel: adapted from magnetic
!
      use Cdata
      use FArrayManager
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(mz) :: ztestfield
      real :: ktestfield_effective
      integer :: jtest, nscal
      real, pointer :: rhoref
!
!  Precalculate etatest if 1/etatest (==etatest1) is given instead
!
      if (etatest1/=0.) etatest=1./etatest1
!
!  Precalculate nutest if 1/nutest (==nutest1) is given instead
!
      if (nutest1/=0.) nutest=1./nutest1
!
!  square of sound speed
!
      cs20test=cs0test**2
!
!  inverse of reference density
!
      if (lmagnetic) then
        call get_shared_variable('rhoref', rhoref, caller='initialize_testfield')
        if (rhoref==impossible) &
          call warning('initialize_testfield','rhoref in magnetic = impossible')
      else
        allocate(rhoref); rhoref=impossible
      endif
!
      if (rho0test==0.) then
        if (rhoref/=impossible.and.rhoref/=0.) then
          rho0test1=1./rhoref
          call warning('initialize_testfield','rho0test=0, set to rhoref from magnetic')
        else
          rho0test1=1.
          call warning('initialize_testfield','rho0test=0, set to 1')
        endif
      else
        if (lmagnetic.and.rho0test/=rhoref) then
          if (lroot) print*, 'Warning - initialize_testfield: rho0test=', rho0test, 'but rhoref from magnetic =', rhoref
        endif
        rho0test1=1./rho0test
      endif
!
!  set cosine and sine function for setting test fields and analysis
!  Choice of using rescaled z-array or original z-array
!  Define ktestfield_effective to deal with boxes bigger than 2pi.
!
      if (ltestfield_newz) then
        ktestfield_effective=ktestfield*(2.*pi/Lz)
        ztestfield=ktestfield_effective*(z-z0)-pi
      else
        ktestfield_effective=ktestfield
        ztestfield=z*ktestfield_effective
      endif
      cz=cos(ztestfield)
      sz=sin(ztestfield)
      c2kz=cos(2*ztestfield)
      s2kz=sin(2*ztestfield)
!
!  calculate cosz*sinz, cos^2, and sinz^2, to take moments with
!  of alpij and etaij. This is useful if there is a mean Beltrami field
!  in the main calculation (lmagnetic=.true.) with phase zero.
!  Optionally, one can determine the phase in the actual field
!  and modify the following calculations in testfield_after_boundary.
!  They should also be with respect to k*z, not just z.
!
      c2z=cz**2
      s2z=sz**2
      csz=cz*sz
!
!  debug output
!
      if (lroot.and.ip<14) then
        print*,'cz=',cz
        print*,'sz=',sz
      endif
!
!  Also inverse of ktestfield_effective, but only if different from zero
!
      if (ktestfield==0) then
        ktestfield1=1.
      else
        ktestfield1=1./ktestfield_effective
      endif
!
!  calculate inverse testfield amplitude (unless it is set to zero)
!
      if (bamp==0.) then
        bamp1=1.
        bamp12=1.
      else
        bamp1=1./bamp
        bamp12=1./bamp**2
      endif
!
!  calculate iE0
!
      if (lrun) then
        select case (itestfield)
          case ('B11-B22'); iE0=0
        case default
          call fatal_error('initialize_testfield','undefined itestfield value')
        endselect
      endif
!
!  Override iE0 if njtest is big enough.
!  This method of identifying the location of the reference field is not very elegant.
!
      if (njtest>=5) iE0=njtest
!
!  Suppress integration of the test problems proper for lzero_only=.true.
!
      if (lzero_only.and.iE0>0) then
        jtest_start=iE0
      else
        jtest_start=1
      endif
!
!  rescale the testfield
!  (in future, could call something like init_aa_simple)
!
      if (reinitialize_aatest) then
!
!  Don't rescale the reference field, so if existent rescale only to njtest-1
!  on njtest-2.
!
        if (iE0>=5) then
          nscal=4
        else
          nscal=njtest
        endif
        do jtest=jtest_start,nscal
          iaxtest=iaatest+3*(jtest-1); iaztest=iaxtest+2
          iuxtest=iuutest+3*(jtest-1); iuztest=iuxtest+2
          ihxtest=ihhtest+  (jtest-1)
          f(:,:,:,iaxtest:iaztest)=rescale_aatest(jtest)*f(:,:,:,iaxtest:iaztest)
          f(:,:,:,iuxtest:iuztest)=rescale_uutest(jtest)*f(:,:,:,iuxtest:iuztest)
          f(:,:,:,ihxtest        )=rescale_hhtest(jtest)*f(:,:,:,ihxtest)
        enddo
      endif
!
!  set lrescaling_testfield=T if linit_aatest=T
!
      if (linit_aatest) lrescaling_testfield=.true.

      if (luse_main_run) then
        if (.not.(lmagnetic.and.lhydro)) then
          luse_main_run=.false.
          call warning('initialize_testfield',  &
                       'Main run missing or incomplete -> switching to kinematic mode')
        elseif (.not.ldensity) then
          call warning('initialize_testfield', &
                       "Main run doesn't contain density -> fluctuating grad(h) is set to zero")
        endif
      endif
!
!  Register an extra aux slot for uxb if requested (so uxb is written
!  to snapshots and can be easily analyzed later). For this to work you
!  must reserve enough auxiliary workspace by setting, for example,
!     ! MAUX CONTRIBUTION 9
!  in the beginning of your src/cparam.local file, *before* setting
!  ncpus, nprocy, etc.
!
!  After a reload, we need to rewrite index.pro, but the auxiliary
!  arrays are already allocated and must not be allocated again.
!
      if (luxb_as_aux) then
        if (iuxbtest==0) then
          call farray_register_auxiliary('uxb',iuxbtest,vector=3,array=njtest)
        else
          if (lroot) print*, 'initialize_testfield: iuxbtest = ', iuxbtest
          call farray_index_append('iuxbtest',iuxbtest)
        endif
      endif
!
!  possibility of using jxb as auxiliary array
!
      if (ljxb_as_aux) then
        if (ijxbtest==0) then
          call farray_register_auxiliary('jxb',ijxbtest,vector=3,array=njtest)
        else
          if (lroot) print*, 'initialize_testfield: ijxbtest = ', ijxbtest
          call farray_index_append('ijxbtest',ijxbtest)
        endif
      endif
!
!  possibility of using ugu as auxiliary array
!
      if (lugu_as_aux) then
        if (iugutest==0) then
          call farray_register_auxiliary('ugu',iugutest,vector=3,array=njtest)
        else
          if (lroot) print*, 'initialize_testfield: iugutest = ', iugutest
          call farray_index_append('iugutest',iugutest)
        endif
      endif
!
!  possibility of using ugh as auxiliary array
!
      if (lugh_as_aux) then
        if (iughtest==0) then
          call farray_register_auxiliary('ugh',iughtest,array=njtest)
        else
          if (lroot) print*, 'initialize_testfield: iughtest = ', iughtest
          call farray_index_append('iughtest',iughtest)
        endif
      endif
!
!  possibility of using Sgh as auxiliary array
!
      if (lSgh_as_aux) then
        if (iSghtest==0) then
          call farray_register_auxiliary('Sgh',iSghtest,vector=3,array=njtest)
        else
          if (lroot) print*, 'initialize_testfield: iSghtest = ', iSghtest
          call farray_index_append('iSghtest',iSghtest)
        endif
      endif
!
      lBext=any(Btest_ext /= 0.)
!
      lremove_E0 = lremove_E0 .and. lmagnetic
      lremove_F0 = lremove_F0 .and. lhydro
      lremove_Q0 = lremove_Q0 .and. ldensity
!
      if (ivid_bb11/=0) then
        !call alloc_slice_buffers(bb11_xy,bb11_xz,bb11_yz,bb11_xy2,bb11_xy3,bb11_xy4,bb11_xz2)
        if (lwrite_slice_xy .and..not.allocated(bb11_xy) ) allocate(bb11_xy (nx,ny,3))
        if (lwrite_slice_xz .and..not.allocated(bb11_xz) ) allocate(bb11_xz (nx,nz,3))
        if (lwrite_slice_yz .and..not.allocated(bb11_yz) ) allocate(bb11_yz (ny,nz,3))
        if (lwrite_slice_xy2.and..not.allocated(bb11_xy2)) allocate(bb11_xy2(nx,ny,3))
        if (lwrite_slice_xy3.and..not.allocated(bb11_xy3)) allocate(bb11_xy3(nx,ny,3))
        if (lwrite_slice_xy4.and..not.allocated(bb11_xy4)) allocate(bb11_xy4(nx,ny,3))
        if (lwrite_slice_xz2.and..not.allocated(bb11_xz2)) allocate(bb11_xz2(nx,nz,3))
      endif
!
!  write testfield information to a file (for convenient post-processing)
!
      if (lroot) then
        open(1,file=trim(datadir)//'/testfield_info.dat',STATUS='unknown')
        write(1,'(a,i1)') 'zextent=',merge(1,0,zextent)
        write(1,'(a,i1)') 'lsoca='  ,merge(1,0,lsoca)
        write(1,'(3a)') "itestfield='",trim(itestfield)//"'"
        write(1,'(3a)') "itestfield_method='",trim(itestfield_method)//"'"
        write(1,'(a,f5.2)') 'ktestfield=',ktestfield
        write(1,'(a,f7.4)') 'lin_testfield=',lin_testfield
        write(1,'(a,f7.4)') 'lam_testfield=',lam_testfield
        write(1,'(a,f7.4)') 'om_testfield=', om_testfield
        write(1,'(a,f7.4)') 'delta_testfield=',delta_testfield
        close(1)
      endif
!
    endsubroutine initialize_testfield
!***********************************************************************
    subroutine init_aatest(f)
!
!  initialise testfield; called from start.f90
!
!  27-nov-09/axel: adapted from init_aatest in testfield_z
!
      use Cdata
      use Mpicomm
      use Initcond
      use Sub
      use InitialCondition, only: initial_condition_aatest
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      iaxtest=iaatest; iaztest=iaxtest+2
      iuxtest=iuutest; iuztest=iuxtest+2

      do j=1,ninit
!
!  select an initial condition for aatest
!
        select case (initaatest(j))
        case ('zero'); f(:,:,:,iaatest:iaatest+3*njtest-1)=0.
        case ('gaussian-noise-1'); call gaunoise(amplaatest(j),f,iaxtest+0,iaztest+0)
        case ('gaussian-noise-2'); call gaunoise(amplaatest(j),f,iaxtest+3,iaztest+3)
        case ('gaussian-noise-3'); call gaunoise(amplaatest(j),f,iaxtest+6,iaztest+6)
        case ('gaussian-noise-4'); call gaunoise(amplaatest(j),f,iaxtest+9,iaztest+9)
        case ('gaussian-noise-5'); call gaunoise(amplaatest(j),f,iaxtest+12,iaztest+12)
        case ('sinwave-x-1'); call sinwave(amplaatest(j),f,iaxtest+0+1,kx=kx_aatest(j))
        case ('sinwave-x-2'); call sinwave(amplaatest(j),f,iaxtest+3+1,kx=kx_aatest(j))
        case ('sinwave-x-3'); call sinwave(amplaatest(j),f,iaxtest+6+1,kx=kx_aatest(j))
        case ('sinwave-x-4'); call sinwave(amplaatest(j),f,iaxtest+9+1,kx=kx_aatest(j))
        case ('sinwave-x-5'); call sinwave(amplaatest(j),f,iaxtest+12+1,kx=kx_aatest(j))
        case ('sinwave-x-6'); if (njtest==6) call sinwave(amplaatest(j),f,iaxtest+15+1,kx=kx_aatest(j))
        case ('Beltrami-x-1'); call beltrami(amplaatest(j),f,iaxtest+0,kx=-kx_aatest(j),phase=phasex_aatest(j))
        case ('Beltrami-z-1'); call beltrami(amplaatest(j),f,iaxtest+0,kz=-kz_aatest(j),phase=phasez_aatest(j))
        case ('Beltrami-z-2'); call beltrami(amplaatest(j),f,iaxtest+3,kz=-kz_aatest(j),phase=phasez_aatest(j))
        case ('Beltrami-z-3'); call beltrami(amplaatest(j),f,iaxtest+6,kz=-kz_aatest(j),phase=phasez_aatest(j))
        case ('Beltrami-z-4'); call beltrami(amplaatest(j),f,iaxtest+9,kz=-kz_aatest(j),phase=phasez_aatest(j))
        case ('Beltrami-z-5'); call beltrami(amplaatest(j),f,iaxtest+12,kz=-kz_aatest(j),phase=phasez_aatest(j))
        case ('nothing'); !(do nothing)
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*, 'init_aatest: initaatest=', trim(initaatest(j))
          call stop_it("")
        endselect
!
!  select an initial condition for uutest
!
      select case (inituutest(j))
        case ('zero'); f(:,:,:,iuutest:iuutest+3*njtest-1)=0.
        case ('sinwave-x-1'); call sinwave(ampluutest(j),f,iuxtest+0,kx=kx_aatest(j))
        case ('sinwave-x-2'); call sinwave(ampluutest(j),f,iuxtest+3,kx=kx_aatest(j))
        case ('sinwave-x-3'); call sinwave(ampluutest(j),f,iuxtest+6,kx=kx_aatest(j))
        case ('sinwave-x-4'); call sinwave(ampluutest(j),f,iuxtest+9,kx=kx_aatest(j))
        case ('sinwave-x-5'); call sinwave(ampluutest(j),f,iuxtest+12,kx=kx_aatest(j))
        case ('sinwave-x-6'); if (njtest==6) call sinwave(ampluutest(j),f,iuxtest+15,kx=kx_aatest(j))
        case ('nothing'); !(do nothing)
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*, 'init_aatest: inituutest=', trim(inituutest(j))
          call stop_it("")
        endselect
!
!  select on initial condition for hhtest
!
      select case (inithhtest(j))
        case ('zero'); f(:,:,:,ihhtest:ihhtest+njtest-1)=0.
        case ('sinwave-x-1'); call sinwave(amplhhtest(j),f,ihhtest+0,kx=kx_aatest(j))
        case ('sinwave-x-2'); call sinwave(amplhhtest(j),f,ihhtest+1,kx=kx_aatest(j))
        case ('sinwave-x-3'); call sinwave(amplhhtest(j),f,ihhtest+2,kx=kx_aatest(j))
        case ('sinwave-x-4'); call sinwave(amplhhtest(j),f,ihhtest+3,kx=kx_aatest(j))
        case ('sinwave-x-5'); call sinwave(amplhhtest(j),f,ihhtest+4,kx=kx_aatest(j))
        case ('sinwave-x-6'); if (njtest==6) call sinwave(amplhhtest(j),f,ihhtest+5,kx=kx_aatest(j))
        case ('nothing'); !(do nothing)
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*, 'init_aatest: inithhtest=', trim(inithhtest(j))
          call stop_it("")
        endselect
!
      enddo
!
!  Interface for user's own subroutine
!
      if (linitial_condition) call initial_condition_aatest(f)
!
    endsubroutine init_aatest
!***********************************************************************
    subroutine pencil_criteria_testfield
!
!   All pencils that the Testfield module depends on are specified here.
!
!  26-jun-05/anders: adapted from magnetic
!
      use Cdata
!
      if (luse_main_run) then
        lpenc_requested(i_bbb)=.true.
        lpenc_requested(i_jj)=.true.
        lpenc_requested(i_uu)=.true.
        if (ldensity) lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_sij)=.true.
!
        lpenc_diagnos(i_oo)=.true.
      endif
!
      if (lforcing_cont_uutest.or.lforcing_cont_aatest) lpenc_requested(i_fcont)=.true.

    endsubroutine pencil_criteria_testfield
!***********************************************************************
    subroutine pencil_interdep_testfield(lpencil_in)
!
!  Interdependency among pencils from the Testfield module is specified here.
!
!  26-jun-05/anders: adapted from magnetic
!
      use General, only: keep_compiler_quiet
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_testfield
!***********************************************************************
    subroutine read_testfield_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=testfield_init_pars, IOSTAT=iostat)
!
    endsubroutine read_testfield_init_pars
!***********************************************************************
    subroutine write_testfield_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=testfield_init_pars)
!
    endsubroutine write_testfield_init_pars
!***********************************************************************
    subroutine read_testfield_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=testfield_run_pars, IOSTAT=iostat)
!
    endsubroutine read_testfield_run_pars
!***********************************************************************
    subroutine write_testfield_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=testfield_run_pars)
!
    endsubroutine write_testfield_run_pars
!***********************************************************************
    subroutine daatest_dt(f,df,p)
!
!  testfield evolution:
!
!  calculate da^(pq)/dt=Uxb^(pq)+uxB^(pq)+uxb-<uxb>+eta*del2A^(pq),
!  and du^(pq)/dt=Jxb^(pq)+jxB^(pq)+jxb-<jxb>+eta*del2U^(pq)-Ugu-ugU-gh,
!  and dh^(pq)/dt=-Ugh-ugH+G'-cs2*divu^(pq)
!    where p=1,2 and q=1 (if B11-B21) and optionally q=2 (if B11-B22),
!  and  da^(0)/dt=(uxb)'+Eext+eta*del2A^(0),
!  with du^(0)/dt=(jxb)'+Fext+eta*del2U^(0)-(ugu)'-cs2*gradh^(0),
!   and dh^(0)/dt=G'-cs2*divu^(0)-(ugh)'
!
!  also calculate corresponding Lorentz force in connection with
!  testflow method
!
!   3-jun-05/axel: coded
!  16-mar-08/axel: Lorentz force added for testfield method
!  25-jan-09/axel: added Maxwell stress tensor calculation
!  27-nov-09/axel: adapted from testfield_z, and added velocity equation
!  15-feb-10/axel: adapted from testfield_nonlinear_z, and added enthalpy
!  27-sep-13/MR  : changes due to uxbtestmz(mz,...  --> uxbtestmz(nz,...
!
      use Cdata
      use Diagnostics
      use Hydro, only: uumz,guumz,lcalc_uumeanz,coriolis_cartesian
      use Density, only: lnrhomz,glnrhomz,lcalc_glnrhomean
      use Magnetic, only: aamz,bbmz,jjmz,lcalc_aameanz,B_ext_inv
      use Mpicomm, only: stop_it
      use Sub
      use Slices_methods, only: store_slices
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      real, dimension (nx,3) :: uxB,bbtest,B0_imposed,uum,umxbtest
      real, dimension (nx,3) :: B0test=0,J0test=0
      real, dimension (nx,3) :: duxbtest,djxbtest,dugutest,dSghtest,eetest
      real, dimension (nx,3) :: uxbtest,uxbtestK,uxbtestM,uxbtestMK
      real, dimension (nx,3) :: jxbtest,jxbtestK,jxbtestM,jxbtestMK
      real, dimension (nx,3,3,njtest) :: Mijpq
      real, dimension (nx,njtest) :: hpq
      real, dimension (nx,3,njtest) :: Fipq,Eipq,upq,bpq,jpq
      real, dimension (nx,njtest) :: Rpq
      real, dimension (nx) :: alpK,alpM,alpMK,phiK,phiM,phiMK
      real, dimension (nx) :: advec_uu, advec_va2
      real, dimension (nx,3) :: del2Atest,uufluct,bbfluct,jjfluct
      real, dimension (nx,3) :: graddivAtest,aatest,jjtest
      real, dimension (nx,3) :: jxbrtest,jxbtest1,jxbtest2
      real, dimension (nx,3) :: del2Utest,uutest,ugutest,Sghtest
      real, dimension (nx,3) :: ugum,ghhtest,graddivutest,sglnrho
      real, dimension (nx,3) :: u0ref,b0ref,j0ref
      real, dimension (nx,3,3) :: aijtest,bijtest,uijtest,sijtest,Mijtest,uijm
      real, dimension (nx) :: jbpq,upq2,jpq2,bpq2,Epq2,s2kzDF1,s2kzDF2
      real, dimension (nx), parameter :: unity=1.
      real, dimension (nx) :: ughm,ughtest,dughtest,divutest,diffus_eta
      real, dimension (nz) :: rhomean1
      integer :: jtest,j,nl
      integer, parameter :: i1=1, i2=2, i3=3, i4=4
      logical,save :: ltest_uxb=.false.,ltest_jxb=.false.
      logical,save :: ltest_ugu=.false.,ltest_ugh=.false., ltest_Sgh=.false.
      real :: fac
!
      intent(in)   :: f,p
      intent(inout):: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'daatest_dt: SOLVE'
      if (headtt) then
        if (iaxtest /= 0) call identify_bcs('Axtest',iaxtest)
        if (iaytest /= 0) call identify_bcs('Aytest',iaytest)
        if (iaztest /= 0) call identify_bcs('Aztest',iaztest)
      endif
!
      !rhomean1=exp(-lnrhomz)

      nl=n-n1+1
!
!  Calculate uufluct=U-Umean.
!  Note that uumz has dimensions (mz,3), not (nz,3).
!
      if (luse_main_run) then
        if (lcalc_uumeanz) then
          do j=1,3
            uufluct(:,j)=p%uu(:,j)-uumz(n,j)
          enddo

          uijm=0.
          do j=1,3
            uum(:,j)=uumz(n,j)
            uijm(:,j,3)=guumz(nl,j)
          enddo

        else
          uufluct=p%uu
        endif
!
!  Calculate bbfluct=B-Bmean and jjfluct=J-Jmean.
!  Note that, unlike uumz, bbmz and jjmz have dimensions nz*3.
!
        if (lcalc_aameanz) then
          do j=1,3
            bbfluct(:,j)=p%bbb(:,j)-bbmz(nl,j)
            jjfluct(:,j)=p%jj(:,j)-jjmz(nl,j)
          enddo
        else
          bbfluct=p%bbb
          jjfluct=p%jj
        endif
      endif
!
!  loop over all fields, but do it backwards,
!  so we compute the zero field first
!
      advec_uu=0.; advec_va2=0.
!
!testfields:
      do jtest=njtest,jtest_start,-1

        iaxtest=iaatest+3*(jtest-1); iaztest=iaxtest+2
        iuxtest=iuutest+3*(jtest-1); iuztest=iuxtest+2
        ihxtest=ihhtest+   jtest-1 
!
!  compute bbtest et al.
!
        aatest=f(l1:l2,m,n,iaxtest:iaztest)
        uutest=f(l1:l2,m,n,iuxtest:iuztest)
        call gij(f,iuxtest,uijtest,1)
        call gij_etc(f,iuxtest,uutest,uijtest,DEL2=del2Utest,GRADDIV=graddivutest)
        call gij(f,iaxtest,aijtest,1)
        call gij_etc(f,iaxtest,aatest,aijtest,bijtest,DEL2=del2Atest)
        call curl_mn(aijtest,bbtest,aatest)
        call curl_mn(bijtest,jjtest,bbtest)
        call div_mn(uijtest,divutest,uutest)
        call grad(f,ihxtest,ghhtest)
!if (lroot.and.n==n1.and.m==m1.and.jtest==5) print*, 'ghhtest=',maxval(ghhtest)
!
!  Set u0ref, b0ref, and j0ref (if iE0=5 or 6).
!
        if (jtest==iE0) then
          u0ref=uutest
          b0ref=bbtest
          j0ref=jjtest
          if (.not.luse_main_run) then
            uufluct=u0ref
            bbfluct=b0ref
            jjfluct=j0ref
          endif
        endif
!
!  Magnetic diffusion term, valid for eta=const.
!
        !df(l1:l2,m,n,iaxtest:iaztest)= df(l1:l2,m,n,iaxtest:iaztest)-etatest*jjtest      ! Weyl gauge!!!
        df(l1:l2,m,n,iaxtest:iaztest)= df(l1:l2,m,n,iaxtest:iaztest)+etatest*del2Atest
!
!  Viscous force, valid for nu=const.
!
        if (lvisc_simplified_testfield) then
          df(l1:l2,m,n,iuxtest:iuztest) = df(l1:l2,m,n,iuxtest:iuztest)+nutest*del2Utest
        else
          df(l1:l2,m,n,iuxtest:iuztest) = df(l1:l2,m,n,iuxtest:iuztest) &
                                         +nutest*(del2Utest+1./3.*graddivutest)
          if (luse_main_run.and.lSghtest) then
!
!  Contribution 2*nu(stest.grad Hmean + Smean.grad htest).
!
            sglnrho=0.
            if (ldensity.and.lcalc_glnrhomean) then
              call traceless_strain(uijtest,divutest,sijtest,uutest,.true.)
              sglnrho = sglnrho+sijtest(:,:,3)*glnrhomz(nl)                             ! sij^T.grad(Hbar)
            endif
            if (lcalc_uumeanz) sglnrho(:,3) = sglnrho(:,3) + guumz(nl,3)*ghhtest(:,3)   ! Sijbar.grad(hh^T)

            df(l1:l2,m,n,iuxtest:iuztest) = df(l1:l2,m,n,iuxtest:iuztest)+nutest*2.*sglnrho
          endif
        endif
!
!  Centrifugal and Coriolis forces.
!
        if (lhydro) then
          call coriolis_cartesian(df,uutest,iuxtest)
        elseif (Omega/=0.) then
          df(l1:l2,m,n,iuxtest  )=df(l1:l2,m,n,iuxtest  )+2*Omega*uutest(:,2)
          df(l1:l2,m,n,iuxtest+1)=df(l1:l2,m,n,iuxtest+1)-2*Omega*uutest(:,1)
        endif
!
!  With imposed field, calculate uutest x B0 and jjtest x B0 terms.
!  This applies to all terms, including the reference fields.
!
        if (lBext) then
          do j=1,3
            B0_imposed(:,j)=Btest_ext(j)
          enddo
          call cross_mn(uutest,B0_imposed,uxbtest)
          call cross_mn(jjtest,B0_imposed,jxbtest)
          df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest)+uxbtest
          df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest)+jxbtest*rho0test1
        endif
!
!  hh-dependent terms
!
        if (lfprestest) &
          df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest)-cs20test*ghhtest ! ... -c_s^2 grad h^T; h=ln(rho)

        df(l1:l2,m,n,ihxtest)=df(l1:l2,m,n,ihxtest)-divutest                           ! ... -div h^T
!
!  add Ubar x b^0 and Ubar x b^T terms to induction equation
!  -Ubar.gradu^0 and -Ubar.gradu^T, as well as
!  -u^0.gradUbar and -u^T.gradUbar to momentum equation
!  and -Ubar.gradh^0, -Ubar.gradh^T, -u^0.gradHbar, -u^T.gradHbar to continuity equation
!
        if (luse_main_run) then
          if (lcalc_uumeanz) then    ! lcalc_uumeanz should be enforced
!
!  terms in induction equation
!
            call cross_mn(uum,bbtest,umxbtest)                                         ! Ubar x b^T
            df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest)+umxbtest
!
!  terms in momentum equation
!
            if (lugutest) then
              call u_dot_grad(f,iuxtest,uijtest,uum   ,ugum,UPWIND=lupw_uutest        )  ! Ubar.grad u^T
              call u_dot_grad(f,iuxtest,uijm   ,uutest,ugum,UPWIND=.false.,LADD=.true.)  ! u^T.grad Ubar
              df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest)-ugum
            endif
!
!  terms in enthalpy equation
!
            call dot(uum,ghhtest,ughm)                                                   ! Ubar.grad h^T
          else
            ughm=0.
          endif
!
!  mean density terms
!
!if (lroot.and.m==m1.and.n==n1) print*, 'glnrhomz,ughm=', maxval(abs(glnrhomz)), maxval(abs(ughm))
          if (ldensity.and.lcalc_glnrhomean) &
            ughm = ughm + uutest(:,3)*glnrhomz(nl)                                     ! u^T.grad Hbar

          df(l1:l2,m,n,ihxtest)=df(l1:l2,m,n,ihxtest)-ughm
        endif
!
!  possibility of non-soca terms (always used for reference problem)
!
        if (.not.lsoca .or. jtest==iE0) then
!
!  subtract average emf, unless we ignore the (uxb)' term (lignore_uxbtestm=T)
!
          if (iuxbtest/=0.and..not.ltest_uxb) then                                ! here uxb is aux!
            uxbtest=f(l1:l2,m,n,iuxbtest+3*(jtest-1):iuxbtest+3*jtest-1)
            if (lignore_uxbtestm) then
              duxbtest=uxbtest
            else
              do j=1,3
                duxbtest(:,j)=uxbtest(:,j)-uxbtestmz(nl,j,jtest)                  ! (utest x btest)'
              enddo
            endif
            if (jtest/=iE0.and.damp_uxb>0.) duxbtest=damp_uxb*duxbtest
            df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest) + duxbtest
          endif
!
!  subtract average jxb, unless we ignore the (jxb)' term (lignore_jxbtestm=T)
!
          if (ijxbtest/=0.and..not.ltest_jxb) then                                ! here jxb is aux!
            jxbtest=f(l1:l2,m,n,ijxbtest+3*(jtest-1):ijxbtest+3*jtest-1)
            if (lignore_jxbtestm) then
              djxbtest=jxbtest
            else
              do j=1,3
                djxbtest(:,j)=jxbtest(:,j)-jxbtestmz(nl,j,jtest)                  ! (jtest x btest)'
              enddo
            endif
            df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest) + djxbtest*rho0test1    ! (jxb)'/rhobar
          endif
!
!  subtract average ugu, unless we ignore the (ugu)' term (lignore_ugutestm=T)
!
          if (lugutest) then
            if (iugutest/=0.and..not.ltest_ugu) then                            ! here ugu is aux!
              ugutest=f(l1:l2,m,n,iugutest+3*(jtest-1):iugutest+3*jtest-1)
              if (lignore_ugutestm) then
                dugutest=ugutest
              else
                do j=1,3
                  dugutest(:,j)=ugutest(:,j)-ugutestmz(nl,j,jtest)              ! (utest.grad utest)'
                enddo
              endif
              df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest) - dugutest
            endif
          endif
!
!  subtract average Sgh, unless we ignore the (Sgh)' term (lignore_Sghtestm=T)
!
          if (.not.lvisc_simplified_testfield.and.lSghtest) then
            if (iSghtest/=0.and..not.ltest_Sgh) then                            ! here Sgh is aux!
              Sghtest=f(l1:l2,m,n,iSghtest+3*(jtest-1):iSghtest+3*jtest-1)
              if (lignore_Sghtestm) then
                dSghtest=Sghtest
              else
                do j=1,3
                  dSghtest(:,j)=Sghtest(:,j)-Sghtestmz(nl,j,jtest)              ! (Stest.grad htest)'
                enddo
              endif
              df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest) + 2.*nutest*dSghtest
            endif
          endif
!
!  subtract average ugh, unless we ignore the (ugh)' term (lignore_ughtestm=T)
!
          if (iughtest/=0.and..not.ltest_ugh) then                            ! here ugh is aux!
            ughtest=f(l1:l2,m,n,iughtest+(jtest-1))
            if (lignore_ughtestm) then
              dughtest=ughtest
            else
              dughtest=ughtest-ughtestmz(nl,jtest)                            ! (utest.grad htest)'
            endif

!if (lroot.and.n==n1.and.m==m1.and.jtest==iE0) print*, 'ughtest,ughtestmz,dughtest=',maxval(ughtest),ughtestmz(nl,jtest),maxval(dughtest)
            df(l1:l2,m,n,ihxtest)=df(l1:l2,m,n,ihxtest)-dughtest
          endif
!
!  end of .not.lsoca 
!
        endif
!
        if (jtest==iE0) then
!
!  Add possibility of forcing that is not delta-correlated in time.
!  (no pseudo-enthalpy forcing at the moment)
!
          if (lforcing_cont_uutest) &
            df(l1:l2,m,n,iuxtest:iuztest)= df(l1:l2,m,n,iuxtest:iuztest) &
                                          +ampl_fcont_uutest*p%fcont(:,:,1)    ! Indices 1,2 universal?
          if (lforcing_cont_aatest) &
            df(l1:l2,m,n,iaxtest:iaztest)= df(l1:l2,m,n,iaxtest:iaztest) &
                                          +ampl_fcont_aatest*p%fcont(:,:,2)
!
        elseif (jtest<=4) then                               ! exclude homogeneous test problem, generalize?
!
!  do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!
          select case (itestfield)
          case ('B11-B22')
            call set_bbtest_B11_B22(B0test,jtest)
            call set_J0test_B11_B22(J0test,jtest)
          case default
            call fatal_error('daatest_dt','undefined itestfield value')
          endselect
!
!  u x B^T
!
          call cross_mn(uufluct,B0test,uxB)
!
!  jxB^T + J^Txb
!
          call cross_mn(jjfluct,B0test,jxbtest1)             ! jxB^T
          call cross_mn(J0test,bbfluct,jxbtest2)             ! J^Txb
          !call multsv_mn(p%rho1,jxbtest1+jxbtest2,jxbrtest)
          jxbrtest=jxbtest1+jxbtest2
          df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest)+uxB
          df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest)+jxbrtest*rho0test1
!
!  r.h.s. of test problems completed
!
        endif
!
!  Restore uxbtest, jxbtest, ugutest and ughtest from f-array (if not already done above).
!
        if (lsoca .and. jtest/=iE0) then
          if (luxb_as_aux) uxbtest=f(l1:l2,m,n,iuxbtest+3*(jtest-1):iuxbtest+3*jtest-1)
          if (ljxb_as_aux) jxbtest=f(l1:l2,m,n,ijxbtest+3*(jtest-1):ijxbtest+3*jtest-1)
          if (lugu_as_aux) ugutest=f(l1:l2,m,n,iugutest+3*(jtest-1):iugutest+3*jtest-1)
          if (lugh_as_aux) ughtest=f(l1:l2,m,n,iughtest+jtest-1)
        endif
!
        if ((ldiagnos.or.l1davgfirst).and. &
          (lsoca.or.ltest_uxb.or.idiag_b0rms/=0.or.idiag_bhrms/=0.or. &
           idiag_j11rms/=0.or.idiag_b11rms/=0.or.idiag_b21rms/=0.or. &
           idiag_b12rms/=0.or.idiag_b22rms/=0.or. &
           idiag_s2kzDFm/=0.or. &
           idiag_M11cc/=0.or.idiag_M11ss/=0.or. &
           idiag_M22cc/=0.or.idiag_M22ss/=0.or. &
           idiag_M12cs/=0.or. &
           idiag_M11/=0.or.idiag_M22/=0.or.idiag_M33/=0.or. &
           idiag_M11z/=0.or.idiag_M22z/=0.or.idiag_M33z/=0)) then
!
!  Settings for quadratic quantities (magnetic stress)
!
          if (idiag_M11cc/=0.or.idiag_M11ss/=0.or. &
              idiag_M22cc/=0.or.idiag_M22ss/=0.or. &
              idiag_M12cs/=0.or. &
              idiag_M11/=0.or.idiag_M22/=0.or.idiag_M33/=0.or. &
              idiag_M11z/=0.or.idiag_M22z/=0.or.idiag_M33z/=0) then
            call dyadic2(bbtest,Mijtest)
            Mijpq(:,:,:,jtest)=Mijtest*bamp12
          endif
        endif
!
        if (ldiagnos.or.l1davgfirst.or.lvideo) then

          bpq(:,:,jtest)=bbtest
          upq(:,:,jtest)=uutest
          if (idiag_jb0m/=0.or.idiag_j11rms/=0) jpq(:,:,jtest)=jjtest
          hpq(:,jtest)=f(l1:l2,m,n,ihxtest)
!
!  evaluate different contributions to <uxb>, <jxb>/rho0 + <u.gradu> + 2*nu*<S'*grad h> and u.grad h.
!
          if (jtest==iE0) then
            fac=1.
          else
            fac=bamp1
          endif
          Eipq(:,:,jtest)=uxbtest*fac
          Fipq(:,:,jtest)=rho0test1*jxbtest
          if (lugutest) Fipq(:,:,jtest)=Fipq(:,:,jtest)-ugutest
          if (.not.lvisc_simplified_testfield.and.lSghtest) Fipq(:,:,jtest)=Fipq(:,:,jtest)+2.*nutest*Sghtest
          Fipq(:,:,jtest)=Fipq(:,:,jtest)*fac
          Rpq(:,jtest)=-ughtest*fac
        endif
!
        if (lfirst.and.ldt) then
!
! advective and Alfven timestep
!
          advec_uu=max(advec_uu,sum(abs(uutest)*dline_1,2))
          advec_va2=max(advec_va2,sum((bbtest*dline_1)**2,2)*(mu01*rho0test1)) 
        endif
      enddo   ! jtest loop = loop over test problems
!
!  compute kinetic, magnetic, and magneto-kinetic contributions
!
      if (luse_main_run.and.any(B_ext_inv/=0.)) then
        call cross_mn(p%uu,p%bbb-b0ref,uxbtestK)
        call cross_mn(p%uu-u0ref,p%bbb,uxbtestM)
        call cross_mn(p%uu-u0ref,p%bbb-b0ref,uxbtestMK)
        call cross_mn(p%jj,p%bbb-b0ref,jxbtestK)
        call cross_mn(p%jj-j0ref,p%bbb,jxbtestM)
        call cross_mn(p%jj-j0ref,p%bbb-b0ref,jxbtestMK)
        call dot(B_ext_inv,uxbtestK,alpK)
        call dot(B_ext_inv,uxbtestM,alpM)
        call dot(B_ext_inv,uxbtestMK,alpMK)
        call dot(B_ext_inv,jxbtestK,phiK)
        call dot(B_ext_inv,jxbtestM,phiM)
        call dot(B_ext_inv,jxbtestMK,phiMK)
      endif
!
      if (lremove_E0) then
        do j=1,3 
          df(l1:l2,m,n,iax+j-1) = df(l1:l2,m,n,iax+j-1)-uxbtestmz(nl,j,iE0)
        enddo
      endif
      if (lremove_F0) then
!if (ldiagnos.and.m==m1.and.lroot) print*, 'maxmin(ugutestmz etc.)=', &
!maxval(abs(ugutestmz(nl,:,iE0)-jxbtestmz(nl,:,iE0)*rho0test1-2.*nutest*Sghtestmz(nl,:,iE0)))
        do j=1,3 
          df(l1:l2,m,n,iux+j-1) = df(l1:l2,m,n,iux+j-1)+ugutestmz(nl,j,iE0)-jxbtestmz(nl,j,iE0)*rho0test1-2.*nutest*Sghtestmz(nl,j,iE0)
        enddo
      endif
      if (lremove_Q0) then
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho)+ughtestmz(nl,iE0)
      endif
!
      if (lfirst.and.ldt) then
        maxadvec=maxadvec+advec_uu
        advec2=advec2+advec_va2
!
        if (headtt.or.ldebug) print*,'daatest_dt: max(advec_uu) =',maxval(advec_uu)
!
!  diffusive time step, just take the max of diffus_eta (if existent)
!  and whatever is calculated here
!
        diffus_eta=max(etatest,nutest)*dxyz_2
        maxdiffus=max(maxdiffus,diffus_eta)
      endif
!
!  in the following block, we have already swapped the 4-6 entries with 7-9
!  The g95 compiler doesn't like to see an index that is out of bounds,
!  so prevent this warning by writing i3=3 and i4=4
!
      if (l1davgfirst) then
!
        call xysum_mn_name_z(bpq(:,1,iE0),idiag_bx0mz)
        call xysum_mn_name_z(bpq(:,2,iE0),idiag_by0mz)
        call xysum_mn_name_z(bpq(:,3,iE0),idiag_bz0mz)
        call xysum_mn_name_z(Eipq(:,1,i1),idiag_E111z)
        call xysum_mn_name_z(Eipq(:,2,i1),idiag_E211z)
        call xysum_mn_name_z(Eipq(:,3,i1),idiag_E311z)
        call xysum_mn_name_z(Eipq(:,1,i2),idiag_E121z)
        call xysum_mn_name_z(Eipq(:,2,i2),idiag_E221z)
        call xysum_mn_name_z(Eipq(:,3,i2),idiag_E321z)
        call xysum_mn_name_z(Eipq(:,1,i3),idiag_E112z)
        call xysum_mn_name_z(Eipq(:,2,i3),idiag_E212z)
        call xysum_mn_name_z(Eipq(:,3,i3),idiag_E312z)
        call xysum_mn_name_z(Eipq(:,1,i4),idiag_E122z)
        call xysum_mn_name_z(Eipq(:,2,i4),idiag_E222z)
        call xysum_mn_name_z(Eipq(:,3,i4),idiag_E322z)
        call xysum_mn_name_z(Eipq(:,1,iE0),idiag_E10z)
        call xysum_mn_name_z(Eipq(:,2,iE0),idiag_E20z)
        call xysum_mn_name_z(Eipq(:,3,iE0),idiag_E30z)
        call xysum_mn_name_z(Mijpq(:,1,1,i1),idiag_M11z)
        call xysum_mn_name_z(Mijpq(:,2,2,i1),idiag_M22z)
        call xysum_mn_name_z(Mijpq(:,3,3,i1),idiag_M33z)
!
      endif
!
      if (ldiagnos) then
!
!  Volume averages for alpha and eta.
!
        if (idiag_alp11/=0) call sum_mn_name(+cz(n)*Eipq(:,1,1)+sz(n)*Eipq(:,1,i2),idiag_alp11)
        if (idiag_alp21/=0) call sum_mn_name(+cz(n)*Eipq(:,2,1)+sz(n)*Eipq(:,2,i2),idiag_alp21)
        if (idiag_alp31/=0) call sum_mn_name(+cz(n)*Eipq(:,3,1)+sz(n)*Eipq(:,3,i2),idiag_alp31)
        if (idiag_eta12/=0) call sum_mn_name(-(-sz(n)*Eipq(:,1,i1)+cz(n)*Eipq(:,1,i2))*ktestfield1,idiag_eta12)
        if (idiag_eta22/=0) call sum_mn_name(-(-sz(n)*Eipq(:,2,i1)+cz(n)*Eipq(:,2,i2))*ktestfield1,idiag_eta22)
!
!  same for jxb
!
        if (idiag_phi11/=0) call sum_mn_name(+cz(n)*Fipq(:,1,1)+sz(n)*Fipq(:,1,i2),idiag_phi11)
        if (idiag_phi21/=0) call sum_mn_name(+cz(n)*Fipq(:,2,1)+sz(n)*Fipq(:,2,i2),idiag_phi21)
        if (idiag_psi12/=0) call sum_mn_name(-(-sz(n)*Fipq(:,1,i1)+cz(n)*Fipq(:,1,i2))*ktestfield1,idiag_psi12)
        if (idiag_psi22/=0) call sum_mn_name(-(-sz(n)*Fipq(:,2,i1)+cz(n)*Fipq(:,2,i2))*ktestfield1,idiag_psi22)
!
!  weighted averages alpha and eta
!  Still need to do this for iE0 /= 0 case.
!
        if (idiag_alp11cc/=0) call sum_mn_name(c2z(n)*(+cz(n)*Eipq(:,1,1)+sz(n)*Eipq(:,1,i2)),idiag_alp11cc)
        if (idiag_alp21sc/=0) call sum_mn_name(csz(n)*(+cz(n)*Eipq(:,2,1)+sz(n)*Eipq(:,2,i2)),idiag_alp21sc)
        if (idiag_eta12cs/=0) call sum_mn_name(-csz(n)*(-sz(n)*Eipq(:,1,i1)+cz(n)*Eipq(:,1,i2))*ktestfield1,idiag_eta12cs)
        if (idiag_eta22ss/=0) call sum_mn_name(-s2z(n)*(-sz(n)*Eipq(:,2,i1)+cz(n)*Eipq(:,2,i2))*ktestfield1,idiag_eta22ss)
!
!  Compute kinetic, magnetic, and magneto-kinetic contributions with
!  imposed-field method
!
        if (idiag_alpK/=0) call sum_mn_name(alpK,idiag_alpK)
        if (idiag_alpM/=0) call sum_mn_name(alpM,idiag_alpM)
        if (idiag_alpMK/=0) call sum_mn_name(alpMK,idiag_alpMK)
        if (idiag_phiK/=0) call sum_mn_name(phiK,idiag_phiK)
        if (idiag_phiM/=0) call sum_mn_name(phiM,idiag_phiM)
        if (idiag_phiMK/=0) call sum_mn_name(phiMK,idiag_phiMK)
!
!  Divergence of current helicity flux
!
        if (idiag_s2kzDFm/=0) then
          eetest=etatest*jjtest-duxbtest
          s2kzDF1=2*ktestfield*c2kz(n)*(&
            bijtest(:,1,3)*eetest(:,1)+&
            bijtest(:,2,3)*eetest(:,2)+&
            bijtest(:,3,3)*eetest(:,3))
          s2kzDF2=4*ktestfield**2*s2kz(n)*(&
            bbtest(:,1)*eetest(:,1)+&
            bbtest(:,2)*eetest(:,2))
          call sum_mn_name(s2kzDF1-s2kzDF2,idiag_s2kzDFm)
        endif
!
!  Maxwell tensor and its weighted averages
!
        if (idiag_M11/=0)   call sum_mn_name(       Mijpq(:,1,1,i1),idiag_M11)
        if (idiag_M22/=0)   call sum_mn_name(       Mijpq(:,2,2,i1),idiag_M22)
        if (idiag_M33/=0)   call sum_mn_name(       Mijpq(:,3,3,i1),idiag_M33)
        if (idiag_M11cc/=0) call sum_mn_name(c2z(n)*Mijpq(:,1,1,i1),idiag_M11cc)
        if (idiag_M11ss/=0) call sum_mn_name(s2z(n)*Mijpq(:,1,1,i1),idiag_M11ss)
        if (idiag_M22cc/=0) call sum_mn_name(c2z(n)*Mijpq(:,2,2,i1),idiag_M22cc)
        if (idiag_M22ss/=0) call sum_mn_name(s2z(n)*Mijpq(:,2,2,i1),idiag_M22ss)
        if (idiag_M12cs/=0) call sum_mn_name(cz(n)*sz(n)*Mijpq(:,1,2,i1),idiag_M12cs)
!
!  Projection of EMF from testfield against testfield itself
!
        if (idiag_EBpq/=0) call sum_mn_name(cz(n)*Eipq(:,1,1) &
                                           +sz(n)*Eipq(:,2,1),idiag_EBpq)
!
!  print warning if alp12 and alp12 are needed, but njtest is too small XX
!
        if (njtest<=2 .and. &
            (idiag_alp12/=0.or.idiag_alp22/=0.or.idiag_alp32/=0.or. &
            (leta_rank2.and.(idiag_eta11/=0.or.idiag_eta21/=0)).or. &
            (.not.leta_rank2.and.(idiag_eta12/=0.or.idiag_eta22/=0)).or. &
            (leta_rank2.and.(idiag_eta11cc/=0.or.idiag_eta21sc/=0)))) then
          call stop_it('njtest is too small if alpi2 or etai2 for i=1,2,3 are needed')
        else
!
!  Remaining coefficients
!
          if (idiag_alp12/=0) call sum_mn_name(+cz(n)*Eipq(:,1,i3)+sz(n)*Eipq(:,1,i4),idiag_alp12)
          if (idiag_alp22/=0) call sum_mn_name(+cz(n)*Eipq(:,2,i3)+sz(n)*Eipq(:,2,i4),idiag_alp22)
          if (idiag_alp32/=0) call sum_mn_name(+cz(n)*Eipq(:,3,i3)+sz(n)*Eipq(:,3,i4),idiag_alp32)
          if (idiag_alp12cs/=0) call sum_mn_name(csz(n)*(+cz(n)*Eipq(:,1,i3)+sz(n)*Eipq(:,1,i4)),idiag_alp12cs)
          if (idiag_alp22ss/=0) call sum_mn_name(s2z(n)*(+cz(n)*Eipq(:,2,i3)+sz(n)*Eipq(:,2,i4)),idiag_alp22ss)
          if (idiag_eta11/=0) call sum_mn_name((-sz(n)*Eipq(:,1,i3)+cz(n)*Eipq(:,1,i4))*ktestfield1,idiag_eta11)
          if (idiag_eta21/=0) call sum_mn_name((-sz(n)*Eipq(:,2,i3)+cz(n)*Eipq(:,2,i4))*ktestfield1,idiag_eta21)
          if (idiag_eta11cc/=0) call sum_mn_name(c2z(n)*(-sz(n)*Eipq(:,1,i3)+cz(n)*Eipq(:,1,i4))*ktestfield1,idiag_eta11cc)
          if (idiag_eta21sc/=0) call sum_mn_name(csz(n)*(-sz(n)*Eipq(:,2,i3)+cz(n)*Eipq(:,2,i4))*ktestfield1,idiag_eta21sc)
!
!  same for jxb/rho0 + u.grad u.
!
          if (idiag_phi12/=0) call sum_mn_name(+cz(n)*Fipq(:,1,i3)+sz(n)*Fipq(:,1,i4),idiag_phi12)
          if (idiag_phi22/=0) call sum_mn_name(+cz(n)*Fipq(:,2,i3)+sz(n)*Fipq(:,2,i4),idiag_phi22)
          if (idiag_phi32/=0) call sum_mn_name(+cz(n)*Fipq(:,3,i3)+sz(n)*Fipq(:,3,i4),idiag_phi32)
          if (idiag_psi11/=0) call sum_mn_name((-sz(n)*Fipq(:,1,i3)+cz(n)*Fipq(:,1,i4))*ktestfield1,idiag_psi11)
          if (idiag_psi21/=0) call sum_mn_name((-sz(n)*Fipq(:,2,i3)+cz(n)*Fipq(:,2,i4))*ktestfield1,idiag_psi21)
!
!  and u.grad h
!

          if (idiag_sig1/=0) call sum_mn_name(+cz(n)*Rpq(:,i1)+sz(n)*Rpq(:,i2),idiag_sig1)
          if (idiag_sig2/=0) call sum_mn_name(+cz(n)*Rpq(:,i3)+sz(n)*Rpq(:,i4),idiag_sig2)
          if (idiag_tau1/=0) call sum_mn_name((-sz(n)*Rpq(:,i1)+cz(n)*Rpq(:,i2))*ktestfield1,idiag_tau1)
          if (idiag_tau2/=0) call sum_mn_name((-sz(n)*Rpq(:,i3)+cz(n)*Rpq(:,i4))*ktestfield1,idiag_tau2)
!
        endif
!
!  Volume-averaged dot products of mean emf and velocity and of mean emf and vorticity
!
        if (luse_main_run.and.iE0/=0) then
          if (idiag_E0Um/=0) call sum_mn_name(uxbtestmz(nl,1,iE0)*p%uu(:,1) &
                                             +uxbtestmz(nl,2,iE0)*p%uu(:,2) &
                                             +uxbtestmz(nl,3,iE0)*p%uu(:,3),idiag_E0Um)
          if (idiag_E0Wm/=0) call sum_mn_name(uxbtestmz(nl,1,iE0)*p%oo(:,1) &
                                             +uxbtestmz(nl,2,iE0)*p%oo(:,2) &
                                             +uxbtestmz(nl,3,iE0)*p%oo(:,3),idiag_E0Wm)
        endif
        if (iE0/=0.and.idiag_E0mrms/=0) call sum_mn_name(sum(uxbtestmz(nl,:,iE0)**2)*unity,idiag_E0mrms,lsqrt=.true.)
!
!  diagnostics for delta function driving, but doesn't seem to work
!
        if (idiag_bamp/=0) call sum_mn_name(bamp*unity,idiag_bamp)
!
!  diagnostics for single points
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_bx0pt/=0)  call save_name(bpq(lpoint-nghost,1,iE0),idiag_bx0pt)
          if (idiag_bx11pt/=0) call save_name(bpq(lpoint-nghost,1,i1),idiag_bx11pt)
          if (idiag_bx21pt/=0) call save_name(bpq(lpoint-nghost,1,i2),idiag_bx21pt)
          if (idiag_bx12pt/=0) call save_name(bpq(lpoint-nghost,1,i3),idiag_bx12pt)
          if (idiag_bx22pt/=0) call save_name(bpq(lpoint-nghost,1,i4),idiag_bx22pt)
          if (idiag_by0pt/=0)  call save_name(bpq(lpoint-nghost,2,iE0),idiag_by0pt)
          if (idiag_by11pt/=0) call save_name(bpq(lpoint-nghost,2,i1),idiag_by11pt)
          if (idiag_by21pt/=0) call save_name(bpq(lpoint-nghost,2,i2),idiag_by21pt)
          if (idiag_by12pt/=0) call save_name(bpq(lpoint-nghost,2,i3),idiag_by12pt)
          if (idiag_by22pt/=0) call save_name(bpq(lpoint-nghost,2,i4),idiag_by22pt)
          if (idiag_Ex0pt/=0)  call save_name(Eipq(lpoint-nghost,1,iE0),idiag_Ex0pt)
          if (idiag_Ex11pt/=0) call save_name(Eipq(lpoint-nghost,1,i1),idiag_Ex11pt)
          if (idiag_Ex21pt/=0) call save_name(Eipq(lpoint-nghost,1,i2),idiag_Ex21pt)
          if (idiag_Ex12pt/=0) call save_name(Eipq(lpoint-nghost,1,i3),idiag_Ex12pt)
          if (idiag_Ex22pt/=0) call save_name(Eipq(lpoint-nghost,1,i4),idiag_Ex22pt)
          if (idiag_Ey0pt/=0)  call save_name(Eipq(lpoint-nghost,2,iE0),idiag_Ey0pt)
!         if (idiag_bamp/=0)   call save_name(bamp,idiag_bamp)
          if (idiag_Ey11pt/=0) call save_name(Eipq(lpoint-nghost,2,i1),idiag_Ey11pt)
          if (idiag_Ey21pt/=0) call save_name(Eipq(lpoint-nghost,2,i2),idiag_Ey21pt)
          if (idiag_Ey12pt/=0) call save_name(Eipq(lpoint-nghost,2,i3),idiag_Ey12pt)
          if (idiag_Ey22pt/=0) call save_name(Eipq(lpoint-nghost,2,i4),idiag_Ey22pt)
        endif
!
!  rms values of small scales fields bpq in response to the test fields Bpq
!  Obviously idiag_b0rms and idiag_b12rms cannot both be invoked!
!  Needs modification!
!
        if (idiag_jb0m/=0) then
          call dot(jpq(:,:,iE0),bpq(:,:,iE0),jbpq)
          call sum_mn_name(jbpq,idiag_jb0m)
        endif
!
        if (idiag_ux0m/=0) call sum_mn_name(upq(:,1,iE0),idiag_ux0m)
        if (idiag_uy0m/=0) call sum_mn_name(upq(:,2,iE0),idiag_uy0m)
        if (idiag_ux11m/=0) call sum_mn_name(upq(:,1,i1),idiag_ux11m)
        if (idiag_uy11m/=0) call sum_mn_name(upq(:,2,i1),idiag_uy11m)
!
        if (idiag_u0rms/=0) then
          call dot2(upq(:,:,iE0),bpq2)
          call sum_mn_name(bpq2,idiag_u0rms,lsqrt=.true.)
        endif
        if (idiag_u0max/=0) then
          call dot2(upq(:,:,iE0),bpq2)
          call max_mn_name(bpq2,idiag_u0max,lsqrt=.true.)
        endif
        if (idiag_h0rms/=0) call sum_mn_name(hpq(:,iE0)**2,idiag_h0rms,lsqrt=.true.)
        call max_mn_name(hpq(:,iE0),idiag_h0max)
        if (idiag_rho0m/=0) call sum_mn_name(exp(hpq(:,iE0)),idiag_rho0m)
!
        if (idiag_b0rms/=0) then
          call dot2(bpq(:,:,iE0),bpq2)
          call sum_mn_name(bpq2,idiag_b0rms,lsqrt=.true.)
        endif
        if (idiag_b0max/=0) then
          call dot2(bpq(:,:,iE0),bpq2)
          call max_mn_name(bpq2,idiag_b0max,lsqrt=.true.)
        endif
        if (idiag_bhrms/=0.and.njtest==6) then
          call dot2(bpq(:,:,iE0-1),bpq2)
          call sum_mn_name(bpq2,idiag_bhrms,lsqrt=.true.)
        endif
!
        if (idiag_u11rms/=0) then
          call dot2(upq(:,:,1),upq2)
          call sum_mn_name(upq2,idiag_u11rms,lsqrt=.true.)
        endif
        if (idiag_h11rms/=0) call sum_mn_name(hpq(:,1)**2,idiag_h11rms,lsqrt=.true.)
!
        if (idiag_u21rms/=0) then
          call dot2(upq(:,:,i2),upq2)
          call sum_mn_name(upq2,idiag_u21rms,lsqrt=.true.)
        endif
        if (idiag_h21rms/=0) call sum_mn_name(hpq(:,i2)**2,idiag_h21rms,lsqrt=.true.)
!
        if (idiag_u12rms/=0) then
          call dot2(upq(:,:,i3),upq2)
          call sum_mn_name(upq2,idiag_u12rms,lsqrt=.true.)
        endif
        if (idiag_h12rms/=0) call sum_mn_name(hpq(:,i3)**2,idiag_h12rms,lsqrt=.true.)
!
        if (idiag_u22rms/=0) then
          call dot2(upq(:,:,i4),upq2)
          call sum_mn_name(upq2,idiag_u22rms,lsqrt=.true.)
        endif
        if (idiag_h22rms/=0) call sum_mn_name(hpq(:,i4)**2,idiag_h22rms,lsqrt=.true.)
!
        if (idiag_b11rms/=0) then
          call dot2(bpq(:,:,1),bpq2)
          call sum_mn_name(bpq2,idiag_b11rms,lsqrt=.true.)
        endif
!
        if (idiag_j11rms/=0) then
          call dot2(jpq(:,:,1),jpq2)
          call sum_mn_name(jpq2,idiag_j11rms,lsqrt=.true.)
          !-test- call sum_mn_name(jpq(:,1,1)**2,idiag_j11rms,lsqrt=.true.)
        endif
!
        if (idiag_b21rms/=0) then
          call dot2(bpq(:,:,i2),bpq2)
          call sum_mn_name(bpq2,idiag_b21rms,lsqrt=.true.)
        endif
!
        if (idiag_b12rms/=0) then
          call dot2(bpq(:,:,i3),bpq2)
          call sum_mn_name(bpq2,idiag_b12rms,lsqrt=.true.)
        endif
!
        if (idiag_b22rms/=0) then
          call dot2(bpq(:,:,i4),bpq2)
          call sum_mn_name(bpq2,idiag_b22rms,lsqrt=.true.)
        endif
!
        if (idiag_E0rms/=0) then
          call dot2(Eipq(:,:,iE0),Epq2)
          call sum_mn_name(Epq2,idiag_E0rms,lsqrt=.true.)
        endif
!
        if (idiag_E11rms/=0) then
          call dot2(Eipq(:,:,i1),Epq2)
          call sum_mn_name(Epq2,idiag_E11rms,lsqrt=.true.)
        endif
!
        if (idiag_E21rms/=0) then
          call dot2(Eipq(:,:,i2),Epq2)
          call sum_mn_name(Epq2,idiag_E21rms,lsqrt=.true.)
        endif
!
        if (idiag_E12rms/=0) then
          call dot2(Eipq(:,:,i3),Epq2)
          call sum_mn_name(Epq2,idiag_E12rms,lsqrt=.true.)
        endif
!
        if (idiag_E22rms/=0) then
          call dot2(Eipq(:,:,i4),Epq2)
          call sum_mn_name(Epq2,idiag_E22rms,lsqrt=.true.)
        endif
!
      endif
!
!  write B-slices for output in wvid in run.f90
!
      if (lvideo.and.lfirst.and.ivid_bb11/=0) &
        call store_slices(bpq(:,:,1),bb11_xy,bb11_xz,bb11_yz,bb11_xy2,bb11_xy3,bb11_xy4,bb11_xz2)
!
    endsubroutine daatest_dt
!***********************************************************************
    subroutine get_slices_testfield(f,slices)
!
!  Write slices for animation of magnetic variables.
!
!  12-sep-09/axel: adapted from the corresponding magnetic routine
!
      use General, only: keep_compiler_quiet
      use Slices_methods, only: assign_slices_vec
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Magnetic field
!
        case ('bb11')
          call assign_slices_vec(slices,bb11_xy,bb11_xz,bb11_yz,bb11_xy2,bb11_xy3,bb11_xy4,bb11_xz2)

      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine get_slices_testfield
!***********************************************************************
    subroutine testfield_before_boundary(f)
!
!  Actions to take before boundary conditions are set.
!
!    4-oct-18/axel+nishant: adapted from testflow
!
      use General, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine testfield_before_boundary
!***********************************************************************
    subroutine testfield_after_boundary(f)
!
!  Calculate <uxb^T> + <u^Txb>, which is needed when lsoca=.false.
!  Also calculate <jxb^T> + <j^Txb>, which is needed when 
!
!  30-nov-09/axel: adapted from testfield_z.f90
!  25-sep-13/MR  : removed parameter p, restricted calculation of pencil case
!  27-sep-13/MR  : changes due to uxbtestmz(mz,...  --> uxbtestmz(nz,...;
!                  pencil calculation corrected; communication simplified
!
      use Cdata
      use Sub
      use Density, only: glnrhomz,lcalc_glnrhomean,calc_pencils_density
      use Hydro, only: calc_pencils_hydro,uumz,guumz,lcalc_uumeanz
      use Magnetic, only: calc_pencils_magnetic, idiag_bcosphz, idiag_bsinphz, &
                          aamz,bbmz,jjmz,lcalc_aameanz
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mz) :: c,s
!
      real, dimension (nx,3,3) :: aijtest,bijtest,uijtest,sijtest,uij0ref
      real, dimension (nx,3) :: aatest,bbtest,jjtest,uutest,uxbtest,jxbtest,Sghtest
      real, dimension (nx,3) :: uxbtest1,jxbtest1
      real, dimension (nx,3) :: uxbtest2,jxbtest2
      real, dimension (nx,3) :: ghhtest,ugutest
      real, dimension (nx,3) :: u0ref,b0ref,j0ref,gh0ref
      real, dimension (nx,3) :: uufluct,ghhfluct,bbfluct,jjfluct
      real, dimension (nx,3,3) :: sijfluct, sij0ref
      real, dimension (nx) :: divutest,ughtest
      real, dimension (mx,3) :: aatestmx,uutestmx
      real, dimension (my,3) :: aatestmy,uutestmy
      real, dimension (3) :: uutestm0
      integer :: jtest,j,juxb,jjxb,jugu,jugh,jSgh,nl,lll,mmm
      logical :: headtt_save
      real :: fac, bcosphz, bsinphz
!
      type (pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc
!
!  In this routine we will reset headtt after the first pencil,
!  so we need to reset it afterwards.
!
      headtt_save=headtt
      fac=1./nxygrid
! 
!  Stop if iE0 is too small.
!
      if (iE0<5) &
        call fatal_error('testfield_after_boundary','need njtest>=5 for u0ref')
!
!  initialize buffer for mean fields
!
      uxbtestmz=0.
      jxbtestmz=0.
      ugutestmz=0.
      ughtestmz=0.
      Sghtestmz=0.
!
      if (luse_main_run) then
        lpenc_loc = .false.
        lpenc_loc((/i_uu,i_uij,i_divu,i_sij,i_bbb,i_bb,i_bij,i_jj,i_lnrho,i_glnrho/))=.true.
      endif
!
!  Start mn loop
!
mn:   do n=n1,n2
        nl=n-n1+1
        do m=m1,m2
!
        if (luse_main_run) then
!
!  Begin by getting/computing fields from main run.
!
          call calc_pencils_hydro(f,p,lpenc_loc)
          call calc_pencils_density(f,p,lpenc_loc)
          call calc_pencils_magnetic(f,p,lpenc_loc)
!
!  Calculate uufluct=U-Umean.
!  Note that uumz has dimensions mz*3, not nz*3.
!
          if (lcalc_uumeanz) then     ! should be enforced
            do j=1,3
              uufluct(:,j)=f(l1:l2,m,n,iux+j-1)-uumz(n,j)
            enddo
          else
            uufluct=f(l1:l2,m,n,iux:iuz)
          endif

          if (ldensity) then
!
!  Calculate ghhfluct = grad(H) - d_z Hmean.
!
            ghhfluct=p%glnrho
            if (lcalc_glnrhomean) ghhfluct(:,3)=ghhfluct(:,3)-glnrhomz(nl)
          else
            ghhfluct=0.
          endif
!
!  Calculate Sfluct = S(U) - S(Umean).
!  Note that glnrhomz has dimensions mz*3, not nz*3.
!
          sijfluct=p%sij
          if (lcalc_uumeanz) sijfluct(:,3,3) = sijfluct(:,3,3) - guumz(nl,3)
!
!
!  Calculate bbfluct=B-Bmean and jjfluct=J-Jmean.
!  Note that, unlike uumz, bbmz and jjmz have dimensions nz*3.
!
          if (lcalc_aameanz) then
            do j=1,3
              bbfluct(:,j)=p%bbb(:,j)-bbmz(nl,j)
              jjfluct(:,j)=p%jj(:,j)-jjmz(nl,j)
            enddo
          else
            bbfluct=p%bbb
            jjfluct=p%jj
          endif
        endif
!
!  Count jtest backward, so we have already access to the reference fields.
!
        do jtest=njtest,jtest_start,-1
!
!  Compute test solutions aatest, bbtest, jjtest, and uutest,
!  and set bbref, jjref, uuref if jtest=iE0 (first loop iteration)
!
          iaxtest=iaatest+3*(jtest-1); iaztest=iaxtest+2
          iuxtest=iuutest+3*(jtest-1); iuztest=iuxtest+2
          ihxtest=ihhtest+  (jtest-1)
          aatest=f(l1:l2,m,n,iaxtest:iaztest)
          uutest=f(l1:l2,m,n,iuxtest:iuztest)
!
!  calculate uutest, bbtest, ghhtest
!
          call gij(f,iaxtest,aijtest,1)
          call gij(f,iuxtest,uijtest,1)
          call gij_etc(f,iaxtest,aatest,aijtest,bijtest)
          call curl_mn(aijtest,bbtest,aatest)
          call curl_mn(bijtest,jjtest,bbtest)
          call grad(f,ihxtest,ghhtest)
          call div_mn(uijtest,divutest,uutest)
          call traceless_strain(uijtest,divutest,sijtest,uutest,.true.)
!
          if (jtest==iE0) then

!  Set u0ref, b0ref, and j0ref (if test points to "zero-problem").
!  Also compute u0 x b0, j0 x b0, u0.grad u0, and u0.grad h0 and put into corresponding arrays.
!  They continue to exist throughout the jtest loop.
!
            u0ref=uutest
            b0ref=bbtest
            j0ref=jjtest
            uij0ref=uijtest
            gh0ref=ghhtest
            sij0ref=sijtest
            call cross_mn(u0ref,b0ref,uxbtest)        ! u0 x b0
            call cross_mn(j0ref,b0ref,jxbtest)        ! j0 x b0
            call u_dot_grad(f,iuxtest,uij0ref,u0ref,ugutest,UPWIND=lupw_uutest)  ! u0.grad u0
            call u_dot_grad(f,ihxtest,gh0ref,u0ref,ughtest,UPWIND=lupw_hhtest)   ! u0.grad h0
            call multmv(sij0ref,gh0ref,Sghtest)       ! S(u0).grad h0
            if (.not.luse_main_run) then
              uufluct=u0ref
              bbfluct=b0ref
              jjfluct=j0ref
              sijfluct=sij0ref
              ghhfluct=gh0ref
            endif
          else
!
!  Calculate uxb and jxb, depending on whether we use
!  Testfield Method (i) or (ii), or mixed ones of either (iii) or (iv).
!
            select case (itestfield_method)
            case ('ju', '(i)')
              call cross_mn(jjfluct,bbtest,jxbtest1)    ! j x btest
              call cross_mn(jjtest,b0ref,jxbtest2)      ! + jtest x b0
              call cross_mn(uufluct,bbtest,uxbtest1)    ! u x btest
              call cross_mn(uutest,b0ref,uxbtest2)      ! + utest x b0
            case ('bb', '(ii)')
              call cross_mn(j0ref,bbtest,jxbtest1)      ! j0 x btest
              call cross_mn(jjtest,bbfluct,jxbtest2)    ! + jtest x b
              call cross_mn(u0ref,bbtest,uxbtest1)      ! u0 x btest
              call cross_mn(uutest,bbfluct,uxbtest2)    ! + utest x b
            case ('bu', '(iii)')
              call cross_mn(j0ref,bbtest,jxbtest1)      ! j0 x btest
              call cross_mn(jjtest,bbfluct,jxbtest2)    ! + jtest x b
              call cross_mn(uufluct,bbtest,uxbtest1)    ! u x btest
              call cross_mn(uutest,b0ref,uxbtest2)      ! + utest x b0
            case ('jb', '(iv)')
              call cross_mn(jjfluct,bbtest,jxbtest1)    ! j x btest
              call cross_mn(jjtest,b0ref,jxbtest2)      ! + jtest x b0
              call cross_mn(u0ref,bbtest,uxbtest1)      ! u0 x btest
              call cross_mn(uutest,bbfluct,uxbtest2)    ! + utest x b
            case default
              call fatal_error('testfield_after_boundary','unknown itestfield_method '//trim(itestfield_method))
            endselect
!
            uxbtest=uxbtest1+uxbtest2
            jxbtest=jxbtest1+jxbtest2
!
!  Only one of the two possibilities for ugu, ugh and sgh, respectively, implemented here.
!
            call u_dot_grad(f,iuxtest,uijtest,uufluct,ugutest,UPWIND=lupw_uutest)                         ! (u.grad)(utest)   tbc
            call u_dot_grad(f,iuutest+3*(njtest-1),uij0ref,uutest,ugutest,UPWIND=lupw_uutest,LADD=.true.) ! (utest.grad)(u0)  tbc
!
            call u_dot_grad(f,ihxtest,ghhtest,uufluct,ughtest,UPWIND=lupw_hhtest)                         ! u.grad(htest)     tbc
            call u_dot_grad(f,ihhtest+3*(njtest-1),gh0ref,uutest,ughtest,UPWIND=lupw_hhtest,LADD=.true.)  ! utest.grad(h0)    tbc
!
            call multmv(sijfluct,ghhtest,Sghtest)                                                         ! S'.grad(htest)
            call multmv(sijtest,gh0ref,Sghtest,ladd=.true.)                                               ! Stest.grad(h0)
!
          endif
!
!  Put uxb, jxb, ugu, ugh and sgh into f-array
!
          if (iuxbtest/=0) then
            juxb=iuxbtest+3*(jtest-1)
            f(l1:l2,m,n,juxb:juxb+2)=uxbtest
          endif
          if (ijxbtest/=0) then
            jjxb=ijxbtest+3*(jtest-1)
            f(l1:l2,m,n,jjxb:jjxb+2)=jxbtest
          endif
          if (iugutest/=0) then 
            jugu=iugutest+3*(jtest-1)
            f(l1:l2,m,n,jugu:jugu+2)=ugutest
          endif
          if (iughtest/=0) then 
            jugh=iughtest+(jtest-1)
            f(l1:l2,m,n,jugh)=ughtest
          endif
          if (iSghtest/=0) then
            jSgh=iSghtest+3*(jtest-1)
            f(l1:l2,m,n,jSgh:jSgh+2)=Sghtest
          endif
!
!  Add corresponding contribution into averaged arrays, uxbtestmz, jxbtestmz.
!  Do the same for ugutestmz, Sghtestmz and ughtestmz.
!
          uxbtestmz(nl,:,jtest)=uxbtestmz(nl,:,jtest)+fac*sum(uxbtest,1)
          jxbtestmz(nl,:,jtest)=jxbtestmz(nl,:,jtest)+fac*sum(jxbtest,1)
          ugutestmz(nl,:,jtest)=ugutestmz(nl,:,jtest)+fac*sum(ugutest,1)
          Sghtestmz(nl,:,jtest)=Sghtestmz(nl,:,jtest)+fac*sum(Sghtest,1)
          ughtestmz(nl,jtest)  =ughtestmz(nl,  jtest)+fac*sum(ughtest)
!
          headtt=.false.
!
!  finish jtest and mn loops
!
        enddo  ! jtest loop
        enddo
      enddo mn
!
!  do communication
!
      call finalize_aver(nprocxy,12,uxbtestmz)
      call finalize_aver(nprocxy,12,jxbtestmz)
      call finalize_aver(nprocxy,12,ugutestmz)
      call finalize_aver(nprocxy,12,ughtestmz)
      call finalize_aver(nprocxy,12,Sghtestmz)

      iaxtest=iaatest+3*(iE0-1)

      if (lremove_meanaa0x_test) then

        fac=1./nyzgrid
        do j=1,3
          do lll=1,mx
            aatestmx(lll,j)=fac*sum(f(lll,m1:m2,n1:n2,iaxtest+j-1))
          enddo
        enddo
        call finalize_aver(nprocyz,23,aatestmx)
!
        do j=1,3
          do lll=1,mx
            f(lll,:,:,iaxtest+j-1) = f(lll,:,:,iaxtest+j-1)-aatestmx(lll,j)
          enddo
        enddo

      endif
      
      if (lremove_meanaa0y_test) then

        fac=1./nxzgrid
        do j=1,3
          do mmm=1,my
            aatestmy(mmm,j)=fac*sum(f(l1:l2,mmm,n1:n2,iaxtest+j-1))
          enddo
        enddo
        call finalize_aver(nprocxz,13,aatestmy)
!
        do j=1,3
          do mmm=1,my
            f(:,mmm,:,iaxtest+j-1) = f(:,mmm,:,iaxtest+j-1)-aatestmy(mmm,j)
          enddo
        enddo

      endif
!
      if (lremove_meanuu0x_test) then

        fac=1./nyzgrid
        do j=1,3
          do lll=1,mx
            uutestmx(lll,j)=fac*sum(f(lll,m1:m2,n1:n2,iuxtest+j-1))
          enddo
        enddo
        call finalize_aver(nprocyz,23,uutestmx)
!
!  Remove volume average from uutestmx as this is already removed by the z-averaging.
!
        uutestm0=sum(uutestmx(l1:l2,:),1)/nxgrid
        call finalize_aver(nprocx,1,uutestm0)
        do j=1,3
          uutestmx(:,j)=uutestmx(:,j)-uutestm0(j)
        enddo
!
        do j=1,3
          do lll=1,mx
            f(lll,:,:,iuxtest+j-1) = f(lll,:,:,iuxtest+j-1)-uutestmx(lll,j)
          enddo
        enddo

      endif
      
      if (lremove_meanuu0y_test) then

        fac=1./nxzgrid
        do j=1,3
          do mmm=1,my
            uutestmy(mmm,j)=fac*sum(f(l1:l2,mmm,n1:n2,iuxtest+j-1))
          enddo
        enddo
        call finalize_aver(nprocxz,13,uutestmy)
!
!  Remove volume average from uutestmy as this is already removed by the z-averaging.
!
        if (.not.lremove_meanuu0x_test) then
          uutestm0=sum(uutestmy(m1:m2,:),1)/nygrid
          call finalize_aver(nprocy,2,uutestm0)
        endif
        do j=1,3
          uutestmy(:,j)=uutestmy(:,j)-uutestm0(j)
        enddo

        do j=1,3
          do mmm=1,my
            f(:,mmm,:,iuxtest+j-1) = f(:,mmm,:,iuxtest+j-1)-uutestmy(mmm,j)
          enddo
        enddo

      endif
!
!  calculate cosz*sinz, cos^2, and sinz^2, to take moments with
!  of alpij and etaij. This is useful if there is a mean Beltrami field
!  in the main calculation (lmagnetic=.true.) with phase zero.
!  Here we modify the calculations depending on the phase of the
!  actual field.
!
!  Calculate phase_testfield (for Beltrami fields)
!
      if (lphase_adjust) then
        if (lroot) then
          if (idiag_bcosphz/=0.and.idiag_bsinphz/=0) then
            bcosphz=fname(idiag_bcosphz)
            bsinphz=fname(idiag_bsinphz)
            phase_testfield=atan2(bsinphz,bcosphz)
          else
            call fatal_error('testfield_after_boundary', &
            'need bcosphz, bsinphz in print.in for lphase_adjust=T')
          endif
        endif
        call mpibcast_real(phase_testfield)
        c=cos(z+phase_testfield)
        s=sin(z+phase_testfield)
        c2z=c**2
        s2z=s**2
        csz=c*s
      endif
      if (ip<9) print*,'iproc,phase_testfield=',iproc,phase_testfield
!
!  reset headtt
!
      headtt=headtt_save
!
    endsubroutine testfield_after_boundary
!***********************************************************************
    subroutine rescaling_testfield(f)
!
!  Rescale testfield by factor rescale_aatest(jtest),
!  which could be different for different testfields
!
!  18-may-08/axel: rewrite from rescaling as used in magnetic
!
      use Cdata
      use Sub
      use Hydro, only: uumz,lcalc_uumeanz
      use Density, only: lnrhomz,lcalc_lnrhomean
      use Magnetic, only: aamz,lcalc_aameanz
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=fnlen) :: file
      logical :: ltestfield_out
      integer,save :: ifirst=0
      integer :: j,jtest,jaatest,juutest,jaa,juu
!
      intent(inout) :: f
!
!  reinitialize aatest periodically if requested
!
      if (linit_aatest) then

        file=trim(datadir)//'/tinit_aatest.dat'
        if (ifirst==0) then
          call read_snaptime(trim(file),taainit,naainit,daainit,t)
          if (taainit==0 .or. taainit < t-daainit) then
            taainit=t+daainit
          endif
          ifirst=1
        endif
!
        if (t >= taainit) then
          do jtest=jtest_start,njtest
            iaxtest=iaatest+3*(jtest-1); iaztest=iaxtest+2
            iuxtest=iuutest+3*(jtest-1); iuztest=iuxtest+2
            ihxtest=ihhtest+  (jtest-1)
            if (jtest/=iE0) then
              f(l1:l2,m1:m2,n1:n2,iaxtest:iaztest)=rescale_aatest(jtest)*f(l1:l2,m1:m2,n1:n2,iaxtest:iaztest)
              f(l1:l2,m1:m2,n1:n2,iuxtest:iuztest)=rescale_uutest(jtest)*f(l1:l2,m1:m2,n1:n2,iuxtest:iuztest)
              f(l1:l2,m1:m2,n1:n2,ihxtest)=rescale_hhtest(jtest)*f(l1:l2,m1:m2,n1:n2,ihxtest)
            endif
          enddo
!
!  Reinitialize reference fields with fluctuations of main run.
!
          if (reinitialize_from_mainrun.and.iE0/=0) then
            if (lcalc_aameanz.and.lcalc_uumeanz.and.lcalc_lnrhomean) then
              jtest=iE0
              iaxtest=iaatest+3*(jtest-1)
              iuxtest=iuutest+3*(jtest-1)
              ihxtest=ihhtest+  (jtest-1)
              do n=n1,n2
                do j=1,3
                  jaatest=iaxtest+j-1;  jaa=iax+j-1
                  juutest=iuxtest+j-1;  juu=iux+j-1
                  f(l1:l2,m1:m2,n,jaatest)=f(l1:l2,m1:m2,n,jaa)-aamz(n,j)
                  f(l1:l2,m1:m2,n,juutest)=f(l1:l2,m1:m2,n,juu)-uumz(n,j)
                enddo
                f(l1:l2,m1:m2,n,ihxtest)=f(l1:l2,m1:m2,n,ihxtest)-lnrhomz(n-n1+1)
              enddo
            else
              call fatal_error('rescaling_testfield', &
                               'reinitialize_from_mainrun needs lcalc_aameanz and lcalc_uumeanz and lcalc_lnrhomean')
            endif
          endif
!
!  Update next time for rescaling
!
          call update_snaptime(file,taainit,naainit,daainit,t,ltestfield_out)
        endif
      endif
!
    endsubroutine rescaling_testfield
!***********************************************************************
    subroutine set_bbtest_Beltrami(B0test,jtest)
!
!  set testfield
!
!  29-mar-08/axel: coded
!
      use Cdata
!
      real, dimension (nx,3) :: B0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: B0test
!
!  set B0test for each of the 9 cases
!
      select case (jtest)
      case (1); B0test(:,1)=bamp*cz(n); B0test(:,2)=bamp*sz(n); B0test(:,3)=0.
      case default; B0test=0.
      endselect
!
    endsubroutine set_bbtest_Beltrami
!***********************************************************************
    subroutine set_bbtest_B11_B21(B0test,jtest)
!
!  set testfield
!
!   3-jun-05/axel: coded
!
      use Cdata
!
      real, dimension (nx,3) :: B0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: B0test
!
!  set B0test for each of the 9 cases
!
      select case (jtest)
      case (1); B0test(:,1)=bamp*cz(n); B0test(:,2)=0.; B0test(:,3)=0.
      case (2); B0test(:,1)=bamp*sz(n); B0test(:,2)=0.; B0test(:,3)=0.
      case default; B0test=0.
      endselect
!
    endsubroutine set_bbtest_B11_B21
!***********************************************************************
    subroutine set_J0test_B11_B21(J0test,jtest)
!
!  set testfield
!
!   3-jun-05/axel: coded
!
      use Cdata
!
      real, dimension (nx,3) :: J0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: J0test
!
!  set J0test for each of the 9 cases
!
      select case (jtest)
      case (1); J0test(:,1)=0.; J0test(:,2)=-bamp*ktestfield*sz(n); J0test(:,3)=0.
      case (2); J0test(:,1)=0.; J0test(:,2)=+bamp*ktestfield*cz(n); J0test(:,3)=0.
      case default; J0test=0.
      endselect
!
    endsubroutine set_J0test_B11_B21
!***********************************************************************
    subroutine set_J0test_B11_B22(J0test,jtest)
!
!  set testfield
!
!   3-jun-05/axel: coded
!
      use Cdata
!
      real, dimension (nx,3) :: J0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: J0test
!
!  set J0test for each of the 9 cases
!
      select case (jtest)
      case (1); J0test(:,1)=0.; J0test(:,2)=-bamp*ktestfield*sz(n); J0test(:,3)=0.
      case (2); J0test(:,1)=0.; J0test(:,2)=+bamp*ktestfield*cz(n); J0test(:,3)=0.
      case (3); J0test(:,1)=+bamp*ktestfield*sz(n); J0test(:,2)=0.; J0test(:,3)=0.
      case (4); J0test(:,1)=-bamp*ktestfield*cz(n); J0test(:,2)=0.; J0test(:,3)=0.
      case default; J0test=0.
      endselect
!
    endsubroutine set_J0test_B11_B22
!***********************************************************************
    subroutine set_bbtest_B11_B22 (B0test,jtest)
!
!  set testfield
!
!   3-jun-05/axel: coded
!
      use Cdata
!
      real, dimension (nx,3) :: B0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: B0test
!
!  set B0test for each of the 4 cases
!
      select case (jtest)
      case (1); B0test(:,1)=bamp*cz(n); B0test(:,2)=0.; B0test(:,3)=0.
      case (2); B0test(:,1)=bamp*sz(n); B0test(:,2)=0.; B0test(:,3)=0.
      case (3); B0test(:,1)=0.; B0test(:,2)=bamp*cz(n); B0test(:,3)=0.
      case (4); B0test(:,1)=0.; B0test(:,2)=bamp*sz(n); B0test(:,3)=0.
      case default; B0test=0.
      endselect
!
    endsubroutine set_bbtest_B11_B22
!***********************************************************************
    subroutine set_bbtest_B11_B22_lin (B0test,jtest)
!
!  set testfield
!
!  25-Mar-09/axel: adapted from set_bbtest_B11_B22
!
      use Cdata
!
      real, dimension (nx,3) :: B0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: B0test
!
!  set B0test for each of the 9 cases
!
      select case (jtest)
      case (1); B0test(:,1)=bamp     ; B0test(:,2)=0.; B0test(:,3)=0.
      case (2); B0test(:,1)=bamp*z(n); B0test(:,2)=0.; B0test(:,3)=0.
      case (3); B0test(:,1)=0.; B0test(:,2)=bamp     ; B0test(:,3)=0.
      case (4); B0test(:,1)=0.; B0test(:,2)=bamp*z(n); B0test(:,3)=0.
      case default; B0test=0.
      endselect
!
    endsubroutine set_bbtest_B11_B22_lin
!***********************************************************************
    subroutine rprint_testfield(lreset,lwrite)
!
!  reads and registers print parameters relevant for testfield fields
!
!   3-jun-05/axel: adapted from rprint_magnetic
!
      use Cdata
      use Diagnostics
!
      integer :: iname,inamez
      logical :: lreset
      logical, optional :: lwrite
!
!  reset everything in case of RELOAD
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_bx0mz=0; idiag_by0mz=0; idiag_bz0mz=0
        idiag_E111z=0; idiag_E211z=0; idiag_E311z=0
        idiag_E121z=0; idiag_E221z=0; idiag_E321z=0
        idiag_E10z=0; idiag_E20z=0; idiag_E30z=0
        idiag_EBpq=0; idiag_E0Um=0; idiag_E0Wm=0
        idiag_alp11=0; idiag_alp21=0; idiag_alp31=0
        idiag_alp12=0; idiag_alp22=0; idiag_alp32=0
        idiag_eta11=0; idiag_eta21=0
        idiag_eta12=0; idiag_eta22=0
        idiag_phi11=0; idiag_phi21=0
        idiag_phi12=0; idiag_phi22=0; idiag_phi32=0
        idiag_psi11=0; idiag_psi21=0
        idiag_psi12=0; idiag_psi22=0
        idiag_sig1=0; idiag_sig2=0; idiag_sig3=0
        idiag_tau1=0; idiag_tau2=0
        idiag_alp11cc=0; idiag_alp21sc=0; idiag_alp12cs=0; idiag_alp22ss=0
        idiag_eta11cc=0; idiag_eta21sc=0; idiag_eta12cs=0; idiag_eta22ss=0
        idiag_alpK=0; idiag_alpM=0; idiag_alpMK=0
        idiag_phiK=0; idiag_phiM=0; idiag_phiMK=0
        idiag_s2kzDFm=0
        idiag_M11=0; idiag_M22=0; idiag_M33=0
        idiag_M11cc=0; idiag_M11ss=0; idiag_M22cc=0; idiag_M22ss=0
        idiag_M12cs=0
        idiag_M11z=0; idiag_M22z=0; idiag_M33z=0; idiag_bamp=0
        idiag_jb0m=0; idiag_u0rms=0; idiag_b0rms=0; idiag_h0rms=0; idiag_E0rms=0; idiag_E0mrms=0
        idiag_u0max=0; idiag_b0max=0; idiag_h0max=0; idiag_rho0m=0; idiag_bhrms=0
        idiag_ux0m=0; idiag_uy0m=0
        idiag_ux11m=0; idiag_uy11m=0
        idiag_u11rms=0; idiag_u21rms=0; idiag_u12rms=0; idiag_u22rms=0
        idiag_h11rms=0; idiag_h21rms=0; idiag_h12rms=0; idiag_h22rms=0
        idiag_j11rms=0; idiag_b11rms=0; idiag_b21rms=0; idiag_b12rms=0; idiag_b22rms=0
        idiag_E11rms=0; idiag_E21rms=0; idiag_E12rms=0; idiag_E22rms=0
        idiag_bx0pt=0; idiag_bx11pt=0; idiag_bx21pt=0; idiag_bx12pt=0; idiag_bx22pt=0
        idiag_by0pt=0; idiag_by11pt=0; idiag_by21pt=0; idiag_by12pt=0; idiag_by22pt=0
        idiag_Ex0pt=0; idiag_Ex11pt=0; idiag_Ex21pt=0; idiag_Ex12pt=0; idiag_Ex22pt=0
        idiag_Ey0pt=0; idiag_Ey11pt=0; idiag_Ey21pt=0; idiag_Ey12pt=0; idiag_Ey22pt=0
        ivid_bb11=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'alp11',idiag_alp11)
        call parse_name(iname,cname(iname),cform(iname),'alp21',idiag_alp21)
        call parse_name(iname,cname(iname),cform(iname),'alp31',idiag_alp31)
        call parse_name(iname,cname(iname),cform(iname),'alp12',idiag_alp12)
        call parse_name(iname,cname(iname),cform(iname),'alp22',idiag_alp22)
        call parse_name(iname,cname(iname),cform(iname),'alp32',idiag_alp32)
        call parse_name(iname,cname(iname),cform(iname),'eta11',idiag_eta11)
        call parse_name(iname,cname(iname),cform(iname),'eta21',idiag_eta21)
        call parse_name(iname,cname(iname),cform(iname),'eta12',idiag_eta12)
        call parse_name(iname,cname(iname),cform(iname),'eta22',idiag_eta22)
        call parse_name(iname,cname(iname),cform(iname),'alpK',idiag_alpK)
        call parse_name(iname,cname(iname),cform(iname),'alpM',idiag_alpM)
        call parse_name(iname,cname(iname),cform(iname),'alpMK',idiag_alpMK)
        call parse_name(iname,cname(iname),cform(iname),'phi11',idiag_phi11)
        call parse_name(iname,cname(iname),cform(iname),'phi21',idiag_phi21)
        call parse_name(iname,cname(iname),cform(iname),'phi12',idiag_phi12)
        call parse_name(iname,cname(iname),cform(iname),'phi22',idiag_phi22)
        call parse_name(iname,cname(iname),cform(iname),'phi32',idiag_phi32)
        call parse_name(iname,cname(iname),cform(iname),'psi11',idiag_psi11)
        call parse_name(iname,cname(iname),cform(iname),'psi21',idiag_psi21)
        call parse_name(iname,cname(iname),cform(iname),'psi12',idiag_psi12)
        call parse_name(iname,cname(iname),cform(iname),'psi22',idiag_psi22)
        call parse_name(iname,cname(iname),cform(iname),'sig1',idiag_sig1)
        call parse_name(iname,cname(iname),cform(iname),'sig2',idiag_sig2)
        call parse_name(iname,cname(iname),cform(iname),'sig3',idiag_sig3)
        call parse_name(iname,cname(iname),cform(iname),'tau1',idiag_tau1)
        call parse_name(iname,cname(iname),cform(iname),'tau2',idiag_tau2)
        call parse_name(iname,cname(iname),cform(iname),'phiK',idiag_phiK)
        call parse_name(iname,cname(iname),cform(iname),'phiM',idiag_phiM)
        call parse_name(iname,cname(iname),cform(iname),'phiMK',idiag_phiMK)
        call parse_name(iname,cname(iname),cform(iname),'alp11cc',idiag_alp11cc)
        call parse_name(iname,cname(iname),cform(iname),'alp21sc',idiag_alp21sc)
        call parse_name(iname,cname(iname),cform(iname),'alp12cs',idiag_alp12cs)
        call parse_name(iname,cname(iname),cform(iname),'alp22ss',idiag_alp22ss)
        call parse_name(iname,cname(iname),cform(iname),'eta11cc',idiag_eta11cc)
        call parse_name(iname,cname(iname),cform(iname),'eta21sc',idiag_eta21sc)
        call parse_name(iname,cname(iname),cform(iname),'eta12cs',idiag_eta12cs)
        call parse_name(iname,cname(iname),cform(iname),'eta22ss',idiag_eta22ss)
        call parse_name(iname,cname(iname),cform(iname),'s2kzDFm',idiag_s2kzDFm)
        call parse_name(iname,cname(iname),cform(iname),'M11',idiag_M11)
        call parse_name(iname,cname(iname),cform(iname),'M22',idiag_M22)
        call parse_name(iname,cname(iname),cform(iname),'M33',idiag_M33)
        call parse_name(iname,cname(iname),cform(iname),'M11cc',idiag_M11cc)
        call parse_name(iname,cname(iname),cform(iname),'M11ss',idiag_M11ss)
        call parse_name(iname,cname(iname),cform(iname),'M22cc',idiag_M22cc)
        call parse_name(iname,cname(iname),cform(iname),'M22ss',idiag_M22ss)
        call parse_name(iname,cname(iname),cform(iname),'M12cs',idiag_M12cs)
        call parse_name(iname,cname(iname),cform(iname),'bx11pt',idiag_bx11pt)
        call parse_name(iname,cname(iname),cform(iname),'bx21pt',idiag_bx21pt)
        call parse_name(iname,cname(iname),cform(iname),'bx12pt',idiag_bx12pt)
        call parse_name(iname,cname(iname),cform(iname),'bx22pt',idiag_bx22pt)
        call parse_name(iname,cname(iname),cform(iname),'bx0pt',idiag_bx0pt)
        call parse_name(iname,cname(iname),cform(iname),'by11pt',idiag_by11pt)
        call parse_name(iname,cname(iname),cform(iname),'by21pt',idiag_by21pt)
        call parse_name(iname,cname(iname),cform(iname),'by12pt',idiag_by12pt)
        call parse_name(iname,cname(iname),cform(iname),'by22pt',idiag_by22pt)
        call parse_name(iname,cname(iname),cform(iname),'by0pt',idiag_by0pt)
        call parse_name(iname,cname(iname),cform(iname),'Ex11pt',idiag_Ex11pt)
        call parse_name(iname,cname(iname),cform(iname),'Ex21pt',idiag_Ex21pt)
        call parse_name(iname,cname(iname),cform(iname),'Ex12pt',idiag_Ex12pt)
        call parse_name(iname,cname(iname),cform(iname),'Ex22pt',idiag_Ex22pt)
        call parse_name(iname,cname(iname),cform(iname),'Ex0pt',idiag_Ex0pt)
        call parse_name(iname,cname(iname),cform(iname),'Ey11pt',idiag_Ey11pt)
        call parse_name(iname,cname(iname),cform(iname),'Ey21pt',idiag_Ey21pt)
        call parse_name(iname,cname(iname),cform(iname),'Ey12pt',idiag_Ey12pt)
        call parse_name(iname,cname(iname),cform(iname),'Ey22pt',idiag_Ey22pt)
        call parse_name(iname,cname(iname),cform(iname),'Ey0pt',idiag_Ey0pt)
        call parse_name(iname,cname(iname),cform(iname),'u11rms',idiag_u11rms)
        call parse_name(iname,cname(iname),cform(iname),'u21rms',idiag_u21rms)
        call parse_name(iname,cname(iname),cform(iname),'u12rms',idiag_u12rms)
        call parse_name(iname,cname(iname),cform(iname),'u22rms',idiag_u22rms)
        call parse_name(iname,cname(iname),cform(iname),'h11rms',idiag_h11rms)
        call parse_name(iname,cname(iname),cform(iname),'h21rms',idiag_h21rms)
        call parse_name(iname,cname(iname),cform(iname),'h12rms',idiag_h12rms)
        call parse_name(iname,cname(iname),cform(iname),'h22rms',idiag_h22rms)
        call parse_name(iname,cname(iname),cform(iname),'j11rms',idiag_j11rms)
        call parse_name(iname,cname(iname),cform(iname),'b11rms',idiag_b11rms)
        call parse_name(iname,cname(iname),cform(iname),'b21rms',idiag_b21rms)
        call parse_name(iname,cname(iname),cform(iname),'b12rms',idiag_b12rms)
        call parse_name(iname,cname(iname),cform(iname),'b22rms',idiag_b22rms)
        call parse_name(iname,cname(iname),cform(iname),'jb0m',idiag_jb0m)
        call parse_name(iname,cname(iname),cform(iname),'bamp',idiag_bamp)
        call parse_name(iname,cname(iname),cform(iname),'ux0m',idiag_ux0m)
        call parse_name(iname,cname(iname),cform(iname),'uy0m',idiag_uy0m)
        call parse_name(iname,cname(iname),cform(iname),'ux11m',idiag_ux11m)
        call parse_name(iname,cname(iname),cform(iname),'uy11m',idiag_uy11m)
        call parse_name(iname,cname(iname),cform(iname),'u0rms',idiag_u0rms)
        call parse_name(iname,cname(iname),cform(iname),'b0rms',idiag_b0rms)
        call parse_name(iname,cname(iname),cform(iname),'h0rms',idiag_h0rms)
        call parse_name(iname,cname(iname),cform(iname),'rho0m',idiag_rho0m)
        call parse_name(iname,cname(iname),cform(iname),'u0max',idiag_u0max)
        call parse_name(iname,cname(iname),cform(iname),'b0max',idiag_b0max)
        call parse_name(iname,cname(iname),cform(iname),'h0max',idiag_h0max)
        call parse_name(iname,cname(iname),cform(iname),'bhrms',idiag_bhrms)
        call parse_name(iname,cname(iname),cform(iname),'E11rms',idiag_E11rms)
        call parse_name(iname,cname(iname),cform(iname),'E21rms',idiag_E21rms)
        call parse_name(iname,cname(iname),cform(iname),'E12rms',idiag_E12rms)
        call parse_name(iname,cname(iname),cform(iname),'E22rms',idiag_E22rms)
        call parse_name(iname,cname(iname),cform(iname),'E0rms',idiag_E0rms)
        call parse_name(iname,cname(iname),cform(iname),'E0mrms',idiag_E0mrms)
        call parse_name(iname,cname(iname),cform(iname),'EBpq',idiag_EBpq)
        call parse_name(iname,cname(iname),cform(iname),'E0Um',idiag_E0Um)
        call parse_name(iname,cname(iname),cform(iname),'E0Wm',idiag_E0Wm)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bx0mz',idiag_bx0mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'by0mz',idiag_by0mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bz0mz',idiag_bz0mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E111z',idiag_E111z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E211z',idiag_E211z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E311z',idiag_E311z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E121z',idiag_E121z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E221z',idiag_E221z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E321z',idiag_E321z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E112z',idiag_E112z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E212z',idiag_E212z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E312z',idiag_E312z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E122z',idiag_E122z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E222z',idiag_E222z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E322z',idiag_E322z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E10z',idiag_E10z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E20z',idiag_E20z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E30z',idiag_E30z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'M11z',idiag_M11z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'M22z',idiag_M22z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'M33z',idiag_M33z)
      enddo
!
!  check for those quantities for which we want video slices
!
      do iname=1,nnamev
        call parse_name(iname,cnamev(iname),cformv(iname),'bb11', ivid_bb11)
      enddo
!
    endsubroutine rprint_testfield

endmodule Testfield
