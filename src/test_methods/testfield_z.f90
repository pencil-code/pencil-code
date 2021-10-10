! $Id$
!  This modules deals with all aspects of testfield fields; if no
!  testfield fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testfield relevant subroutines listed in here.

!  Note: this routine requires that MVAR and MAUX contributions
!  together with njtest are set correctly in the cparam.local file.
!  njtest must be set at the end of the file such that 3*njtest=MVAR.
!
!  Example:
!  ! MVAR CONTRIBUTION 12
!  ! MAUX CONTRIBUTION 12
!  integer, parameter :: njtest=4
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ltestfield = .true.
! CPARAM logical, parameter :: ltestfield_z = .true.
! CPARAM logical, parameter :: ltestfield_xy = .false.
! CPARAM logical, parameter :: ltestfield_xz  = .false.
!
!***************************************************************

module Testfield

  use Cparam
  use Messages
  use Testfield_general

  implicit none

  include '../testfield.h'
!
! Slice precalculation buffers
!
  real, target, dimension(:,:,:), allocatable :: bb11_xy, bb11_xy2, &
                                                 bb11_xz, bb11_yz
!
!  cosine and sine function for setting test fields and analysis
!
  real, dimension(mz) :: cz,sz,c2z,csz,s2z,c2kz,s2kz
  real :: phase_testfield=.0
!
  integer :: iE0=0
!
! run parameters
!
  logical :: lphase_adjust=.false.
  real :: ktestfield=1., ktestfield1=1.
  real :: kdamp_2ndord=0., kdamp_iter=0., dt_iter=0., reduce_iter=1.
  real :: chiraltest=0.
  real :: alpha_incoherent=0.
  real, dimension(nx) :: alpha_tmp
  logical :: ltestfield_newz=.true.
  logical :: llorentzforce_testfield=.false.
  logical :: ltest_uxb=.false.,ltest_jxb=.false.
  logical :: lalpha_incoherent=.false.

  !!!! new input pars
  namelist /testfield_run_pars/ &
       B_ext,reinitialize_aatest,lsoca,lsoca_jxb, &
       etatest,etatest1,etatest_hyper3,iresistivity_test, &
       chiraltest, itestfield,ktestfield, &
       lin_testfield,lam_testfield,om_testfield,delta_testfield, &
       ltestfield_newz,leta_rank2,lphase_adjust, &
       ltestfield_taver,llorentzforce_testfield, &
       ltestfield_profile_eta_z, &
       luxb_as_aux,ljxb_as_aux,lignore_uxbtestm, &
       ltest_uxb,ltest_jxb, &
       lforcing_cont_aatest,ampl_fcont_aatest, &
       daainit,linit_aatest,bamp, &
       rescale_aatest,tau_aatest, &
       lalpha_incoherent, alpha_incoherent, &
!
!                                         the following parameter relevant for artificically introduced 2nd order in time equation for a_test for suppressing
!                                         unstable eigenmodes of the homogeneous equations
!
       kdamp_2ndord, &                  ! damping factor in front of \dot{a_test} in artificial 2nd order equation (if om_testfield=0)
!                                         \dot{\dot{a_test}} + kdamp_2ndord*\dot{a_test} = \etatest\nabla^2 a_test + u x b_test (default=0);
!
!                                         or for testfields with harmonic time dependence (if om_testfield/=0), factor in front of rhs of the artificial
!                                         evolution equation \dot{Re{a_test}} = kdamp_2ndord*(Re{rhs} + omega*Im{a_test})
!
!                                         the following parameters relevant for iterative procedure (see Raedler & Rheinhardt GAFD 2007)
! 
       kdamp_iter, &                    ! controls to which extent the uxb term formed with a_test^(i) is used in iteration step i (default=0)
       dt_iter, &                       ! integration time for a single step of iterative procedure, if =0 no iteration is performed (default=0) 
       reduce_iter                      ! factor by which kdamp_iter is reduced after each iteration (default=1 - no reduction)
!
! diagnostic variables (needs to be consistent with reset list below)
  integer :: idiag_alp11=0      ! DIAG_DOC: $\alpha_{11}$
  integer :: idiag_alp21=0      ! DIAG_DOC: $\alpha_{21}$
  integer :: idiag_alp31=0      ! DIAG_DOC: $\alpha_{31}$
  integer :: idiag_alp12=0      ! DIAG_DOC: $\alpha_{12}$
  integer :: idiag_alp22=0      ! DIAG_DOC: $\alpha_{22}$
  integer :: idiag_alp32=0      ! DIAG_DOC: $\alpha_{32}$
  integer :: idiag_alp13=0      ! DIAG_DOC: $\alpha_{13}$
  integer :: idiag_alp23=0      ! DIAG_DOC: $\alpha_{23}$
  integer :: idiag_eta11=0      ! DIAG_DOC: $\eta_{113}k$ or $\eta_{11}k$ if leta_rank2=T
  integer :: idiag_eta21=0      ! DIAG_DOC: $\eta_{213}k$ or $\eta_{21}k$ if leta_rank2=T
  integer :: idiag_eta31=0      ! DIAG_DOC: $\eta_{313}k$
  integer :: idiag_eta12=0      ! DIAG_DOC: $\eta_{123}k$ or $\eta_{12}k$ if leta_rank2=T
  integer :: idiag_eta22=0      ! DIAG_DOC: $\eta_{223}k$ or $\eta_{22}k$ if leta_rank2=T
  integer :: idiag_eta32=0      ! DIAG_DOC: $\eta_{323}k$
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
  integer :: idiag_b11rms=0     ! DIAG_DOC: $\left<b_{11}^2\right>^{1/2}$
  integer :: idiag_b21rms=0     ! DIAG_DOC: $\left<b_{21}^2\right>^{1/2}$
  integer :: idiag_b12rms=0     ! DIAG_DOC: $\left<b_{12}^2\right>^{1/2}$
  integer :: idiag_b22rms=0     ! DIAG_DOC: $\left<b_{22}^2\right>^{1/2}$
  integer :: idiag_b0rms=0      ! DIAG_DOC: $\left<b_{0}^2\right>^{1/2}$
  integer :: idiag_jb0m=0       ! DIAG_DOC: $\left<jb_{0}\right>$
  integer :: idiag_E11rms=0     ! DIAG_DOC: $\left<{\cal E}_{11}^2\right>^{1/2}$
  integer :: idiag_E21rms=0     ! DIAG_DOC: $\left<{\cal E}_{21}^2\right>^{1/2}$
  integer :: idiag_E12rms=0     ! DIAG_DOC: $\left<{\cal E}_{12}^2\right>^{1/2}$
  integer :: idiag_E22rms=0     ! DIAG_DOC: $\left<{\cal E}_{22}^2\right>^{1/2}$
  integer :: idiag_E0rms=0      ! DIAG_DOC: $\left<{\cal E}_{0}^2\right>^{1/2}$
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
  integer :: idiag_alp11z=0     ! DIAG_DOC: $\alpha_{11}(z,t)$
  integer :: idiag_alp21z=0     ! DIAG_DOC: $\alpha_{21}(z,t)$
  integer :: idiag_alp12z=0     ! DIAG_DOC: $\alpha_{12}(z,t)$
  integer :: idiag_alp22z=0     ! DIAG_DOC: $\alpha_{22}(z,t)$
  integer :: idiag_alp13z=0     ! DIAG_DOC: $\alpha_{13}(z,t)$
  integer :: idiag_alp23z=0     ! DIAG_DOC: $\alpha_{23}(z,t)$
  integer :: idiag_eta11z=0     ! DIAG_DOC: $\eta_{11}(z,t)$
  integer :: idiag_eta21z=0     ! DIAG_DOC: $\eta_{21}(z,t)$
  integer :: idiag_eta12z=0     ! DIAG_DOC: $\eta_{12}(z,t)$
  integer :: idiag_eta22z=0     ! DIAG_DOC: $\eta_{22}(z,t)$
  integer :: idiag_uzjx1z=0     ! DIAG_DOC: $u_z j^{11}_x$
  integer :: idiag_uzjy1z=0     ! DIAG_DOC: $u_z j^{11}_y$
  integer :: idiag_uzjz1z=0     ! DIAG_DOC: $u_z j^{11}_z$
  integer :: idiag_uzjx2z=0     ! DIAG_DOC: $u_z j^{21}_x$
  integer :: idiag_uzjy2z=0     ! DIAG_DOC: $u_z j^{21}_y$
  integer :: idiag_uzjz2z=0     ! DIAG_DOC: $u_z j^{21}_z$
  integer :: idiag_uzjx3z=0     ! DIAG_DOC: $u_z j^{12}_x$
  integer :: idiag_uzjy3z=0     ! DIAG_DOC: $u_z j^{12}_y$
  integer :: idiag_uzjz3z=0     ! DIAG_DOC: $u_z j^{12}_z$
  integer :: idiag_uzjx4z=0     ! DIAG_DOC: $u_z j^{22}_x$
  integer :: idiag_uzjy4z=0     ! DIAG_DOC: $u_z j^{22}_y$
  integer :: idiag_uzjz4z=0     ! DIAG_DOC: $u_z j^{22}_z$
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

  logical :: lcomplex=.false., liter=.false., lfirst_iter=.true., ltestfield_linear=.false.
  real    :: t_iter_last=0., taainit_previous
  integer :: njtestl
  integer, parameter :: nkeep=16
  integer, dimension(nkeep) :: idiag_keep=0
!
!  arrays for horizontally averaged uxb and jxb
!
  real, dimension(nz,3,njtest) :: uxbtestm,jxbtestm    ! TB improved: declare smaller (njtestl) if possible, requires tb allocatable

  contains
!
!***********************************************************************
    subroutine initialize_testfield(f)
!
!  Perform any post-parameter-read initialization
!
!   2-jun-05/axel: adapted from magnetic
!   6-sep-3 /MR: outsourced unspecific stuff to testfield_general
!  19-nov-13/MR: slice buffers dynamically allocated
!
      use Diagnostics, only: gen_form_legend
      use Cdata
      use FarrayManager, only: farray_register_auxiliary, farray_index_append
      use General, only: operator(.in.)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(mz) :: ztestfield, c, s
      real :: ktestfield_effective
      integer :: i, jtest, ierr, iaatest2, iostat
      character(LEN=640) :: fform
      real, dimension(nname) :: buffer

      call initialize_testfield_general(f)
!
!  set cosine and sine function for setting test fields and analysis
!  Choice of using rescaled z-array or original z-array
!  Define ktestfield_effective to deal with boxes bigger than 2pi.
!
      if (ltestfield_newz) then
        ktestfield_effective=ktestfield*(2.*pi/Lz)
        ztestfield=ktestfield_effective*(z-z0-Lz/2)
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
      c=cz
      s=sz
      c2z=c**2
      s2z=s**2
      csz=c*s
!
!  debug output
!
      if (lroot.and.ip<14) then
        print*,'cz=',cz
        print*,'sz=',sz
      endif
!
!  Also calculate its inverse, but only if different from zero
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
!  calculate iE0; set ltestfield_linear in specific cases
!
      ltestfield_linear=.false.
      if (lrun) then
        select case (itestfield)
        case ('Beltrami'); iE0=1
        case ('B11-B22_lin','linear'); iE0=0; ltestfield_linear=.true.
        case ('B11-B21+B=0'); iE0=3
        case ('B11-B22+B=0'); iE0=5
        case ('B11-B21'); iE0=0
        case ('B11-B22'); iE0=0
        case ('B11'); iE0=1
        case ('B12'); iE0=1
        case ('B=0'); iE0=1
        case default
          call fatal_error('initialize_testfield','undefined itestfield value')
        endselect
      endif
!
      liter = dt_iter/=0.                       ! iterative procedure enabled

      if (kdamp_2ndord==0..and..not.liter) then
        njtestl=njtest
      else
        njtestl=njtest/2                        ! TB improved for liter=T: auxiliary variables instead of doubling njtest !!!
      endif

      if (liter) then
!
        if (kdamp_iter==0.) then
          lsoca=.true.
        else
          kdamp_iter = min(kdamp_iter,.9)       ! kdamp_iter limited to 0.9
        endif
!
        t_iter_last=t

        iaatest2 = iaatest+3*njtestl            ! start index of the right hand sides of the iterated test problems
!
! in case of restart: check for non-zero right hand sides, if present it is no longer the first iteration 
!
        if (any(f(:,:,:,iaatest2:iaatest2+3*njtestl-1)/=0.)) lfirst_iter=.false.
!
! indices of those diagnostics whose values accumulate - 
!
        idiag_keep = (/ idiag_alp11,idiag_alp12,idiag_alp21,idiag_alp22, &
                        idiag_eta11,idiag_eta12,idiag_eta21,idiag_eta22, &
                        idiag_alp11cc,idiag_alp21sc,idiag_alp12cs,idiag_alp22ss, &
                        idiag_eta11cc,idiag_eta21sc,idiag_eta12cs,idiag_eta22ss   /)
!
        if (lfirst_iter) then
!
! mark accumulating diagnostics to prevent them from being reset immediately after diagnostic output
!
          fname_keep(idiag_keep)=impossible
!
        else
!
! in case of restart: read the accumulated values from the end of the time series 
!
          call gen_form_legend(fform)
          open(1,file=trim(datadir)//'/time_series.dat',position='append',IOSTAT=iostat)
          backspace(1)
          read(1,trim(fform)) buffer
          close(1)
          fname_keep(idiag_keep) = buffer(idiag_keep)
!
        endif
!
      endif
!
! enable complex calculation for harmonic testfields 
!
      lcomplex = kdamp_2ndord/=0..and.(om_testfield/=0. .or. lam_testfield/=0.)
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
          call farray_register_auxiliary('uxbtest',iuxbtest,vector=3,array=-njtestl)
        else
          if (lroot) print*, 'initialize_testfield: iuxbtest = ', iuxbtest
          call farray_index_append('iuxbtest',iuxbtest)
        endif
      endif
!
!  possibility of using jxb as auxiliary array (is intended to be
!  used in connection with testflow method)
!
      if (ljxb_as_aux) then
        if (ijxbtest==0) then
          call farray_register_auxiliary('jxbtest',ijxbtest,vector=3,array=-njtestl)
        else
          if (lroot) print*, 'initialize_testfield: ijxbtest = ', ijxbtest
          call farray_index_append('ijxbtest',ijxbtest)
        endif
      endif
!
!  incoherent alpha effect
!
      if (lalpha_incoherent) then
        alpha_tmp=sqrt2*alpha_incoherent*cos(x(l1:l2))
      endif
!
!  allocate slice buffers
!
      if (lwrite_slices) then
        if ('bb11' .in. cnamev) &
          allocate(bb11_xy(nx,ny,3), bb11_xy2(nx,ny,3), &
                   bb11_xz(nx,nz,3), bb11_yz(ny,nz,3) )
      endif
!
!  write testfield information to a file (for convenient post-processing)
!
      if (lroot) then
        open(1,file=trim(datadir)//'/testfield_info.dat',STATUS='unknown')
        write(1,'(a,i1)') 'lsoca='  ,merge(1,0,lsoca)
        write(1,'(a,i1)') 'lsoca_jxb='  ,merge(1,0,lsoca_jxb)
        write(1,'(3a)') "itestfield='",trim(itestfield)//"'"
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
!    where p=1,2 and q=1 (if B11-B21) and optionally q=2 (if B11-B22)
!
!  also calculate corresponding Lorentz force in connection with
!  testflow method
!
!   3-jun-05/axel: coded
!  16-mar-08/axel: Lorentz force added for testfield method
!  25-jan-09/axel: added Maxwell stress tensor calculation
!   5-jun-13/MR  : second-order-in-time fake test equations introduced 
!                  (for damping of unwanted unstable solutions)
!   5-jun-13/axel+MR: correction of the former; df(l1:l2,m,n,iaxtest:iaztest) --> daatest
!   6-jun-13/MR: further corrected, alternative formulation added
!  16-aug-13/MR: iterative procedure and complex treatment for harmonic testfields added
!  20-aug-13/MR: calc_uxb and calc_diffusive_part introduced
!  27-sep-13/MR: changes due to uxbtestm(mz,...  -->  uxbtestm(nz,...
!  19-nov-13/MR: complex p=(lam_testfield,om_testfield) in complex calculation branch enabled
!  21-nov-13/MR: suppressed time-dependence of testfield in complex calculation for lam_testfield/=0
!  10-oct-21/axel: added possibility of incoherent alpha effect
!
      use Diagnostics
      use Cdata
      use Hydro, only: uumz,lcalc_uumeanz
      use Mpicomm, only: stop_it
      use Sub
      use General, only: operator(.in.)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)     :: p
      intent(inout)  :: f,df
!
      real, dimension (nx,3) :: uxb,B0test,bbtest,duxbbtest
      real, dimension (nx,3) :: jxbtest,djxbrtest,eetest
      real, dimension (nx,3) :: J0test=0,jxB0rtest,J0xbrtest
      real, dimension (nx,3,3,njtestl) :: Mijpq
      real, dimension (nx,3,njtestl) :: bpq,jpq
      real, dimension (nx,3) :: uufluct,daatest
      real, dimension (nx,3) :: del2Atest2,graddivatest,aatest,jjtest,jxbrtest
      real, dimension (nx,3,3) :: aijtest,bijtest,Mijtest
      real, dimension (nx) :: jbpq,bpq2,Epq2,s2kzDF1,s2kzDF2,diffus_eta
      real, dimension (nx), parameter :: unity=1.
      real, dimension (:,:,:,:), allocatable :: Eipq
!
!  auxiliary arrays for imaginary parts
!
      real, dimension (:,:),     allocatable :: daatest2,uxb2,bbtest2,duxbbtest2
!
      integer :: jtest, j, iuxtest, iuztest
      integer :: i1=1, i2=2, i3=3, i4=4, i5=5, iaxtest2, iaztest2
      integer :: iswitch_iter=0, nl
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'daatest_dt: SOLVE'
      if (headtt) then
        if (iaxtest /= 0) call identify_bcs('Axtest',iaxtest)
        if (iaytest /= 0) call identify_bcs('Aytest',iaytest)
        if (iaztest /= 0) call identify_bcs('Aztest',iaztest)
      endif

      if (lcomplex) then
!
!  allocate as fake complex with Eipq(1,...) - real and Eipq(2,...) - imaginary part
!
        allocate(Eipq(2,nx,3,njtestl))
      else
        allocate(Eipq(1,nx,3,njtestl))
      endif
!
!  calculate uufluct=U-Umean
!
      if (lcalc_uumeanz) then
        do j=1,3
          uufluct(:,j)=p%uu(:,j)-uumz(n,j)
        enddo
      else
        uufluct=p%uu
      endif
!
!  Multiply by exponential factor if lam_testfield is different from zero.
!  Allow also for linearly increasing testfields
!  Keep bamp1=1 for oscillatory test fields.
!
      if (lam_testfield/=0..or.lin_testfield/=0. .or. &
          om_testfield/=0..or.delta_testfield/=0.) then
        bamp=1.
        if (lam_testfield/=0.) then
          if (.not.lcomplex) then
            taainit_previous=taainit-daainit
            bamp=bamp*exp(lam_testfield*(t-taainit_previous))
          endif
          bamp1=1./bamp
        endif
        if (lin_testfield/=0.) then
          taainit_previous=taainit-daainit
          bamp=bamp*lin_testfield*(t-taainit_previous-daainit/2.)
          bamp1=1./bamp
        endif
        if (om_testfield/=0.) then
          if (.not.lcomplex) bamp=bamp*cos(om_testfield*t)
          bamp1=1.
        endif
        if (delta_testfield/=0.) then
          if (t>=delta_testfield_next.and.t<=(delta_testfield_next+dt)) then
            delta_testfield_time=t
            delta_testfield_next=t+delta_testfield
          endif
          if (t>=delta_testfield_time-.1*dt.and.t<=delta_testfield_time+.1*dt) then
            bamp=1.
          else
            bamp=0.
          endif
          bamp1=1.
        endif
      endif
!
! if complex formulation for om_testfield/=0 used, allocate arrays for imaginary parts of auxiliary quantities
!
      if (lcomplex) allocate(daatest2(nx,3),uxb2(nx,3),bbtest2(nx,3),duxbbtest2(nx,3))
!
!  do each of the njtestl testfields at a time,
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further down in the file.
!
      nl=n-n1+1
!
      do jtest=1,njtestl
!
        iaxtest=iaatest+3*(jtest-1)
        iaztest=iaxtest+2
!
!  calculate indices for "second part of a_test" - can be \dot{a_test} or Im{a_test} or rhs of problem (i) of iterative procedure = u x b_test^(i-1)
!
        if (kdamp_2ndord/=0..or.liter) then
          iaxtest2=iaxtest+3*njtestl
          iaztest2=iaxtest2+2
        endif
!
        select case (itestfield)
          case ('Beltrami');    call set_bbtest_Beltrami(B0test,jtest)
          case ('B11-B22_lin','linear'); call set_bbtest_B11_B22_lin(B0test,jtest)
          case ('B11-B21+B=0'); call set_bbtest_B11_B21(B0test,jtest)
          case ('B11-B22+B=0'); call set_bbtest_B11_B22(B0test,jtest)
          case ('B11-B21');     call set_bbtest_B11_B21(B0test,jtest)
          case ('B11-B22');     call set_bbtest_B11_B22(B0test,jtest)
          case ('B11');         call set_bbtest_B11_B22(B0test,jtest)
          case ('B12');         call set_bbtest_B11_B22(B0test,3)
          case ('B=0')                                                    !(dont do anything)
          case default
            call fatal_error('daatest_dt','undefined itestfield value')
        endselect
!
!  add an external field to the testfield, if present
!
        if (B_ext(1)/=0.) B0test(:,1)=B0test(:,1)+B_ext(1)
        if (B_ext(2)/=0.) B0test(:,2)=B0test(:,2)+B_ext(2)
        if (B_ext(3)/=0.) B0test(:,3)=B0test(:,3)+B_ext(3)
!
!  put diffusion into daatest
!
        call calc_diffusive_part(f,p,iaxtest,daatest)
        if (lcomplex) call calc_diffusive_part(f,p,iaxtest2,daatest2)
!
!  add u' \times Btest (in iterative procedure only for first problem)
!
        if (lfirst_iter) then 
          call cross_mn(uufluct,B0test,uxb)
          if (lalpha_incoherent) call multsv_mn_add(alpha_tmp,B0test,uxb)
          daatest=daatest+uxb
        endif
!
!  Add non-SOCA terms:
!  Use f-array for uxbbtest (if space has been allocated for this) and
!  if we don't test (i.e. if ltest_uxb=.false.).
!
        if (.not.lsoca) then
!
          if ((iuxbtest/=0.and..not.ltest_uxb).and.chiraltest==0.) then
            uxb=f(l1:l2,m,n,iuxbtest+3*(jtest-1):iuxbtest+3*jtest-1)
          else
            call calc_uxb(f,p,iaxtest,uxb,bbtest)
          endif
          if (lalpha_incoherent) call multsv_mn_add(alpha_tmp,bbtest,uxb)
!
!  subtract average emf, unless we ignore the <uxb> term (lignore_uxbtestm=T)
!
          if (lignore_uxbtestm) then
            duxbbtest=uxb
          else
            do j=1,3
              duxbbtest(:,j)=uxb(:,j)-uxbtestm(nl,j,jtest)
            enddo
          endif
!
!  in iterative procedure damp u x b_test term with kdamp_iter
!
          if (liter.and.kdamp_iter/=0.) duxbbtest=kdamp_iter*duxbbtest
!
!  add to rhs of test-field equation
!
          daatest=daatest+duxbbtest
!
!  the same for imaginary part
!
          if (lcomplex) then
!
            call calc_uxb(f,p,iaxtest2,uxb2,bbtest2)
!
            if (lignore_uxbtestm) then
              duxbbtest2=uxb2
            else
              do j=1,3
                duxbbtest2(:,j)=uxb2(:,j)-uxbtestm(nl,j,jtest+njtestl)
              enddo
            endif
!
            daatest2=daatest2+duxbbtest2
!
          endif
        endif
!
!  add chiral effect term
!  (if lsoca/=T, it may not work!)
!
        if (chiraltest/=0.) daatest=daatest+chiraltest*bbtest
!
!  add possibility of forcing that is not delta-correlated in time
!
        if (lforcing_cont_aatest) &
          daatest=daatest+ampl_fcont_aatest*p%fcont(:,:,2)    ! if lcomplex=T: forcing assumed real!
!
!  add possibility of artificial friction
!
        if (ltestfield_artifric) then
          daatest=daatest-tau1_aatest*f(l1:l2,m,n,iaxtest:iaztest)
          if (lcomplex) &
            daatest2=daatest2-tau1_aatest*f(l1:l2,m,n,iaxtest2:iaztest2)
        endif
!
!  modification for 2nd order time derivative and complex formulation in harmonically time-dependent case
!
        if (kdamp_2ndord/=0.) then
!
          if (.not.lcomplex) then   ! for 2nd order time derivative
!        
            df(l1:l2,m,n,iaxtest2:iaztest2)=df(l1:l2,m,n,iaxtest2:iaztest2)+daatest                               ! \dot{f2} = rhs
            df(l1:l2,m,n,iaxtest :iaztest )=df(l1:l2,m,n,iaxtest :iaztest )+ &                                    ! \dot{b_test} = f2 - kdamp_2ndord*b_test
                                             f(l1:l2,m,n,iaxtest2:iaztest2)-kdamp_2ndord*f(l1:l2,m,n,iaxtest:iaztest)
!
!  alternatively
!
!            df(l1:l2,m,n,iaxtest2:iaztest2)= df(l1:l2,m,n,iaxtest2:iaztest2)+daatest &                           ! \dot{f2} = \dot{\dot{b_test}} 
!                                                                                                                 !          = rhs - kdamp_2ndord*f2 
!                                            -kdamp_2ndord*f(l1:l2,m,n,iaxtest2:iaztest2)                         !          = rhs - kdamp_2ndord*\dot{b_test}
! 
!            df(l1:l2,m,n,iaxtest :iaztest )= df(l1:l2,m,n,iaxtest :iaztest )+f(l1:l2,m,n,iaxtest2:iaztest2)      ! \dot{b_test} = f2
 
          else                      ! complex formulation for harmonically time-dependent testfields (om_testfield/=0, constant velocity!) 
             
            df(l1:l2,m,n,iaxtest2:iaztest2)= df(l1:l2,m,n,iaxtest2:iaztest2)+daatest2 &                           ! \dot{Im{a_test}} = 
                                            -om_testfield*f(l1:l2,m,n,iaxtest:iaztest)                            ! Im{rhs} - omega*Re{a_test}
            if (lam_testfield/=0.) &
              df(l1:l2,m,n,iaxtest2:iaztest2)= df(l1:l2,m,n,iaxtest2:iaztest2) &                                  ! \dot{Im{a_test}} = 
                                              -lam_testfield*f(l1:l2,m,n,iaxtest2:iaztest2)                       ! \dot{Im{a_test}} - lambda*Im{a_test}
!
            df(l1:l2,m,n,iaxtest:iaztest )=df(l1:l2,m,n,iaxtest:iaztest ) &                                       ! \dot{Re{a_test}} = 
                                            +kdamp_2ndord*(daatest+om_testfield*f(l1:l2,m,n,iaxtest2:iaztest2))   ! kdamp_2ndord*(Re{rhs} + omega*Im{a_test})
            if (lam_testfield/=0.) &
              df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest) &                                       ! \dot{Re{a_test}} = 
                                              -(kdamp_2ndord*lam_testfield)*f(l1:l2,m,n,iaxtest:iaztest)          ! \dot{Re{a_test}} - kdamp_2ndord*lambda*Re{a_test})
          endif
!
        else        ! unmodified test equations
!
          df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest)+daatest
!
!  for iterative procedure
!
          if (liter) then
!
!  for iterative step i>1 add stored right hand side  (1-kdamp_iter)*(u x b_test^(i-1))
! 
            if (.not.lfirst_iter) &
              df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest)+f(l1:l2,m,n,iaxtest2:iaztest2)
!
!  if integration time >= dt_iter
!
            if ( t-t_iter_last >= dt_iter ) then
              
              if (lsoca) call calc_uxb(f,p,iaxtest,uxb,bbtest)
!
!  subtract average emf, unless we ignore the <uxb> term (lignore_uxbtestm=T)
!
              if (lignore_uxbtestm) then
                duxbbtest=uxb
              else
                do j=1,3
                  duxbbtest(:,j)=uxb(:,j)-uxbtestm(nl,j,jtest)
                enddo
              endif
!
!  store (1.-kdamp_iter)*(u x b_test^(i)) for use in next iterative step
!  
              f(l1:l2,m,n,iaxtest2:iaztest2) = (1.-kdamp_iter)*duxbbtest
!
            endif
          endif
        endif
!
!  Calculate Lorentz force for sinlge B11 testfield and add to duu
!
        if (llorentzforce_testfield.and.lhydro) then
!
          aatest=f(l1:l2,m,n,iaxtest:iaztest)
          call gij(f,iaxtest,aijtest,1)
          call gij_etc(f,iaxtest,aatest,aijtest,bijtest,del2Atest2,graddivatest)
!
!  calculate jpq and bpq
!
          call curl_mn(aijtest,bbtest,aatest)
          call curl_mn(bijtest,jjtest,bbtest)
!
!  calculate jpq x B0pq
!
          select case (itestfield)
            case ('B11'); call set_J0test_B11_B21(J0test,jtest)
            case ('B12'); call set_J0test_B11_B21(J0test,jtest)
            case ('B=0') !(dont do anything)
            case default
              call fatal_error('daatest_dt','undefined itestfield value')
          endselect
          call cross_mn(J0test+jjtest,B0test+bbtest,jxbrtest)
          call multsv_mn(p%rho1,jxbrtest,jxbrtest)
!
!  add them all together
!
          df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+jxbrtest
        endif
!
!  Calculate Lorentz force
!
        if (ltestflow.or.(ldiagnos.and.idiag_jb0m/=0)) then
!
          iuxtest=iuutest+4*(jtest-1)                       !!!MR: only correct if jtest refers to number of testflow
          iuztest=iuxtest+2
          aatest=f(l1:l2,m,n,iaxtest:iaztest)
          call gij(f,iaxtest,aijtest,1)
          call gij_etc(f,iaxtest,aatest,aijtest,bijtest,del2Atest2,graddivatest)
!
!  calculate jpq x bpq
!
          call curl_mn(aijtest,bbtest,aatest)
          call curl_mn(bijtest,jjtest,bbtest)
        endif
!
!  calculate jpq x B0pq
!
        if (ltestflow) then
!
          call cross_mn(jjtest,B0test,jxB0rtest)
!
!  calculate J0pq x bpq
!
          select case (itestfield)
!           case ('B11-B21+B=0'); call set_J0test(J0test,jtest)
            case ('B11-B21'); call set_J0test_B11_B21(J0test,jtest)
!           case ('B11-B22'); call set_J0test_B11_B22(J0test,jtest)
            case ('B=0') !(dont do anything)
          case default
            call fatal_error('daatest_dt','undefined itestfield value')
          endselect
          call cross_mn(J0test,bbtest,J0xbrtest)
!
!  add them all together
!
          if (lsoca_jxb) then
            df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest) &
                                          +jxB0rtest+J0xbrtest
          else
!
!  use f-array for uxb (if space has been allocated for this) and
!  if we don't test (i.e. if ltest_jxb=.false.)
!
            if (ijxbtest/=0.and..not.ltest_jxb) then
              jxbtest=f(l1:l2,m,n,ijxbtest+3*(jtest-1):ijxbtest+3*jtest-1)
            else
              call cross_mn(jjtest,bbtest,jxbrtest)
            endif
!
!  subtract average jxb
!
            do j=1,3
              djxbrtest(:,j)=jxbtest(:,j)-jxbtestm(nl,j,jtest)
            enddo
            df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest) &
              +jxB0rtest+J0xbrtest+djxbrtest
!
          endif
        endif
!
!  calculate alpha, begin by calculating uxbbtest (if not already done above)
!
        if ((ldiagnos.or.l1davgfirst).and. &
          (lsoca.or.ltest_uxb.or.idiag_b0rms/=0.or. &
           idiag_b11rms/=0.or.idiag_b21rms/=0.or. &
           idiag_b12rms/=0.or.idiag_b22rms/=0.or. &
           idiag_s2kzDFm/=0.or. &
           idiag_M11cc/=0.or.idiag_M11ss/=0.or. &
           idiag_M22cc/=0.or.idiag_M22ss/=0.or. &
           idiag_M12cs/=0.or. &
           idiag_M11/=0.or.idiag_M22/=0.or.idiag_M33/=0.or. &
           idiag_M11z/=0.or.idiag_M22z/=0.or.idiag_M33z/=0)) then

          call calc_uxb(f,p,iaxtest,uxb,bbtest)
          if (lalpha_incoherent) call multsv_mn_add(alpha_tmp,bbtest,uxb)
          if (lcomplex) call calc_uxb(f,p,iaxtest2,uxb2,bbtest2)

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
        if (ldiagnos) then
!
          Eipq(1,:,:,jtest)=uxb*bamp1
          if (lcomplex) Eipq(2,:,:,jtest)=uxb2*bamp1
!
!  only real part saved here (imaginary is in bbtest2)
!
          bpq(:,:,jtest)=bbtest
!
          if (idiag_jb0m/=0) jpq(:,:,jtest)=jjtest
!
          if (idiag_uzjx1z/=0.or.idiag_uzjy1z/=0.or.idiag_uzjz1z/=0.or. &
              idiag_uzjx2z/=0.or.idiag_uzjy2z/=0.or.idiag_uzjz2z/=0.or. &
              idiag_uzjx3z/=0.or.idiag_uzjy3z/=0.or.idiag_uzjz3z/=0) then
!
            call gij_etc(f,iaxtest,aatest,aijtest,bijtest,del2Atest2,graddivatest)
            call curl_mn(aijtest,bbtest,aatest)
            call curl_mn(bijtest,jjtest,bbtest)
!
            jpq(:,:,jtest)=jjtest
!
          endif
!
        endif
!
      enddo    ! end loop over njtestl testfields
!
      if ( liter .and. n==n2 .and. m==m2 ) then
!
!  for iterative procedure at end of all loops: if integration time >= dt_iter update start time t_iter_last, reduce kdamp_iter
!                                               signal that iteration level has been changed in iswitch_iter
        if (t-t_iter_last >= dt_iter) then
          t_iter_last = t
          if (lfirst_iter) then 
            iswitch_iter=-1
          else
            iswitch_iter=1
          endif
          lfirst_iter = .false.
          kdamp_iter = reduce_iter*kdamp_iter
        endif
!
      endif
!
!  diffusive time step, just take the max of diffus_eta (if existent)
!  and whatever is calculated here
!
      if (lfirst.and.ldt) then
        diffus_eta=etatest*dxyz_2
        maxdiffus=max(maxdiffus,diffus_eta)
      endif
!
      if (ldiagnos) then
!
        if (liter .and. n==n1 .and. m==m1) then
!
!  for iterative procedure at end of first run of all loops:
!  if iteration level has been changed after last diagnostic output, use the values from the preceding iteration level
!
          if (iswitch_iter/=0) then
            if (iswitch_iter==-1) then
!
!  to overwrite impossible values in fname_keep after first iteration
!
              fname_keep(idiag_keep) = fname(idiag_keep)
            else
!
!  to accumulate otherwise
!
              fname_keep(idiag_keep) = fname_keep(idiag_keep)+fname(idiag_keep)
            endif
            iswitch_iter=0
          endif
!
!  reset entries for accumulating values in fname (which are not reset in prints)
!
          fname(idiag_keep) = 0.
!
        endif
!
!  in the following block, we have already swapped the 4-6 entries with 7-9
!  The g95 compiler doesn't like to see an index that is out of bounds,
!  so prevent this warning by writing i3=3 and i4=4
!
!  only real parts of testsolutions and mean EMF taken into account here
!
        if (iE0>0) call xysum_mn_name_z(bpq(:,1,iE0),idiag_bx0mz)
        if (iE0>0) call xysum_mn_name_z(bpq(:,2,iE0),idiag_by0mz)
        if (iE0>0) call xysum_mn_name_z(bpq(:,3,iE0),idiag_bz0mz)
        call xysum_mn_name_z(p%uu(:,3)*jpq(:,1,i1),idiag_uzjx1z)
        call xysum_mn_name_z(p%uu(:,3)*jpq(:,2,i1),idiag_uzjy1z)
        call xysum_mn_name_z(p%uu(:,3)*jpq(:,3,i1),idiag_uzjz1z)
        call xysum_mn_name_z(p%uu(:,3)*jpq(:,1,i2),idiag_uzjx2z)
        call xysum_mn_name_z(p%uu(:,3)*jpq(:,2,i2),idiag_uzjy2z)
        call xysum_mn_name_z(p%uu(:,3)*jpq(:,3,i2),idiag_uzjz2z)
        call xysum_mn_name_z(p%uu(:,3)*jpq(:,1,i3),idiag_uzjx3z)
        call xysum_mn_name_z(p%uu(:,3)*jpq(:,2,i3),idiag_uzjy3z)
        call xysum_mn_name_z(p%uu(:,3)*jpq(:,3,i3),idiag_uzjz3z)
        call xysum_mn_name_z(p%uu(:,3)*jpq(:,1,i4),idiag_uzjx4z)
        call xysum_mn_name_z(p%uu(:,3)*jpq(:,2,i4),idiag_uzjy4z)
        call xysum_mn_name_z(p%uu(:,3)*jpq(:,3,i4),idiag_uzjz4z)
        call xysum_mn_name_z(Eipq(1,:,1,i1),idiag_E111z)
        call xysum_mn_name_z(Eipq(1,:,2,i1),idiag_E211z)
        call xysum_mn_name_z(Eipq(1,:,3,i1),idiag_E311z)
        call xysum_mn_name_z(Eipq(1,:,1,i2),idiag_E121z)
        call xysum_mn_name_z(Eipq(1,:,2,i2),idiag_E221z)
        call xysum_mn_name_z(Eipq(1,:,3,i2),idiag_E321z)
        call xysum_mn_name_z(Eipq(1,:,1,i3),idiag_E112z)
        call xysum_mn_name_z(Eipq(1,:,2,i3),idiag_E212z)
        call xysum_mn_name_z(Eipq(1,:,3,i3),idiag_E312z)
        call xysum_mn_name_z(Eipq(1,:,1,i4),idiag_E122z)
        call xysum_mn_name_z(Eipq(1,:,2,i4),idiag_E222z)
        call xysum_mn_name_z(Eipq(1,:,3,i4),idiag_E322z)
        if (iE0>0) call xysum_mn_name_z(Eipq(1,:,1,iE0),idiag_E10z)
        if (iE0>0) call xysum_mn_name_z(Eipq(1,:,2,iE0),idiag_E20z)
        if (iE0>0) call xysum_mn_name_z(Eipq(1,:,3,iE0),idiag_E30z)
        call xysum_mn_name_z(Mijpq(:,1,1,i1),idiag_M11z)
        call xysum_mn_name_z(Mijpq(:,2,2,i1),idiag_M22z)
        call xysum_mn_name_z(Mijpq(:,3,3,i1),idiag_M33z)
!
        if (ltestfield_linear) then
          call xysum_mn_name_z(Eipq(1,:,1,i1),idiag_alp11z)
          call xysum_mn_name_z(Eipq(1,:,2,i1),idiag_alp21z)
          call xysum_mn_name_z(Eipq(1,:,1,i3),idiag_alp12z)
          call xysum_mn_name_z(Eipq(1,:,2,i3),idiag_alp22z)
          if (njtestl >= 5) then
            call xysum_mn_name_z(Eipq(1,:,1,i5),idiag_alp13z)
            call xysum_mn_name_z(Eipq(1,:,2,i5),idiag_alp23z)
          endif
          call xysum_mn_name_z(+(-z(n)*Eipq(1,:,1,i3)+Eipq(1,:,1,i4)),idiag_eta11z)
          call xysum_mn_name_z(+(-z(n)*Eipq(1,:,2,i3)+Eipq(1,:,2,i4)),idiag_eta21z)
          call xysum_mn_name_z(-(-z(n)*Eipq(1,:,1,i1)+Eipq(1,:,1,i2)),idiag_eta12z)
          call xysum_mn_name_z(-(-z(n)*Eipq(1,:,2,i1)+Eipq(1,:,2,i2)),idiag_eta22z)
        else
          call xysum_mn_name_z(cz(n)*Eipq(1,:,1,i1)+sz(n)*Eipq(1,:,1,i2),idiag_alp11z)
          call xysum_mn_name_z(cz(n)*Eipq(1,:,2,i1)+sz(n)*Eipq(1,:,2,i2),idiag_alp21z)
          call xysum_mn_name_z(cz(n)*Eipq(1,:,1,i3)+sz(n)*Eipq(1,:,1,i4),idiag_alp12z)
          call xysum_mn_name_z(cz(n)*Eipq(1,:,2,i3)+sz(n)*Eipq(1,:,2,i4),idiag_alp22z)
!          if (idiag_alp11z/=0) call fatal_error('daatest_dt','works only for ltestfield_linear')
          call xysum_mn_name_z(-(sz(n)*Eipq(1,:,1,i3)-cz(n)*Eipq(1,:,1,i4))*ktestfield1,idiag_eta11z)
          call xysum_mn_name_z(+(sz(n)*Eipq(1,:,1,i1)-cz(n)*Eipq(1,:,1,i2))*ktestfield1,idiag_eta12z)
          call xysum_mn_name_z(-(sz(n)*Eipq(1,:,2,i3)-cz(n)*Eipq(1,:,2,i4))*ktestfield1,idiag_eta21z)
          call xysum_mn_name_z(+(sz(n)*Eipq(1,:,2,i1)-cz(n)*Eipq(1,:,2,i2))*ktestfield1,idiag_eta22z)
        endif
!
!  averages of alpha and eta
!  Don't subtract E0 field if iE=0
!
        if (iE0==0) then
          if (ltestfield_linear) then
            if (idiag_alp11/=0) call sum_mn_name(Eipq(:,:,1,1),idiag_alp11)
            if (idiag_alp21/=0) call sum_mn_name(Eipq(:,:,2,1),idiag_alp21)
            if (idiag_alp31/=0) call sum_mn_name(Eipq(:,:,3,1),idiag_alp31)
            if (njtestl >= 5) then
              if (idiag_alp13/=0) call sum_mn_name(Eipq(:,:,1,i5),idiag_alp13)
              if (idiag_alp23/=0) call sum_mn_name(Eipq(:,:,2,i5),idiag_alp23)
            endif
            if (idiag_eta12/=0) call sum_mn_name(-(-z(n)*Eipq(:,:,1,i1)+Eipq(:,:,1,i2)),idiag_eta12)
            if (idiag_eta22/=0) call sum_mn_name(-(-z(n)*Eipq(:,:,2,i1)+Eipq(:,:,2,i2)),idiag_eta22)
          else
            if (idiag_alp11/=0) call sum_mn_name(+cz(n)*Eipq(:,:,1,1)+sz(n)*Eipq(:,:,1,i2),idiag_alp11)
            if (idiag_alp21/=0) call sum_mn_name(+cz(n)*Eipq(:,:,2,1)+sz(n)*Eipq(:,:,2,i2),idiag_alp21)
            if (idiag_alp31/=0) call sum_mn_name(+cz(n)*Eipq(:,:,3,1)+sz(n)*Eipq(:,:,3,i2),idiag_alp31)
            if (leta_rank2) then
              if (idiag_eta12/=0) call sum_mn_name(-(-sz(n)*Eipq(:,:,1,i1)+cz(n)*Eipq(:,:,1,i2))*ktestfield1,idiag_eta12)
              if (idiag_eta22/=0) call sum_mn_name(-(-sz(n)*Eipq(:,:,2,i1)+cz(n)*Eipq(:,:,2,i2))*ktestfield1,idiag_eta22)
            else
              if (idiag_eta11/=0) call sum_mn_name((-sz(n)*Eipq(:,:,1,1)+cz(n)*Eipq(:,:,1,i2))*ktestfield1,idiag_eta11)
              if (idiag_eta21/=0) call sum_mn_name((-sz(n)*Eipq(:,:,2,1)+cz(n)*Eipq(:,:,2,i2))*ktestfield1,idiag_eta21)
              if (idiag_eta31/=0) call sum_mn_name((-sz(n)*Eipq(:,:,3,1)+cz(n)*Eipq(:,:,3,i2))*ktestfield1,idiag_eta31)
            endif
          endif
!
!  Subtract E0 field if it has been calculated.
!  Do this only for the leta_rank2=T option.
!
        else
          if (idiag_alp11/=0) call sum_mn_name(  +cz(n)*(Eipq(:,:,1,i1)-Eipq(:,:,1,iE0)) &
                                                 +sz(n)*(Eipq(:,:,1,i2)-Eipq(:,:,1,iE0)),idiag_alp11)
          if (idiag_alp21/=0) call sum_mn_name(  +cz(n)*(Eipq(:,:,2,i1)-Eipq(:,:,2,iE0)) &
                                                 +sz(n)*(Eipq(:,:,2,i2)-Eipq(:,:,2,iE0)),idiag_alp21)
          if (idiag_alp31/=0) call sum_mn_name(  +cz(n)*(Eipq(:,:,3,i1)-Eipq(:,:,3,iE0)) &
                                                 +sz(n)*(Eipq(:,:,3,i2)-Eipq(:,:,3,iE0)),idiag_alp31)
          if (idiag_eta12/=0) call sum_mn_name(-(-sz(n)*(Eipq(:,:,1,i1)-Eipq(:,:,1,iE0)) &
                                                 +cz(n)*(Eipq(:,:,1,i2)-Eipq(:,:,1,iE0)))*ktestfield1,idiag_eta12)
          if (idiag_eta22/=0) call sum_mn_name(-(-sz(n)*(Eipq(:,:,2,i1)-Eipq(:,:,2,iE0)) &
                                                 +cz(n)*(Eipq(:,:,2,i2)-Eipq(:,:,2,iE0)))*ktestfield1,idiag_eta22)
        endif
!
!  weighted averages alpha and eta
!  Still need to do this for iE0 /= 0 case.
!
        if (idiag_alp11cc/=0) call sum_mn_name(c2z(n)*(+cz(n)*Eipq(:,:,1,1)+sz(n)*Eipq(:,:,1,i2)),idiag_alp11cc)
        if (idiag_alp21sc/=0) call sum_mn_name(csz(n)*(+cz(n)*Eipq(:,:,2,1)+sz(n)*Eipq(:,:,2,i2)),idiag_alp21sc)
        if (leta_rank2) then
          if (idiag_eta12cs/=0) call sum_mn_name(-csz(n)*(-sz(n)*Eipq(:,:,1,i1)+cz(n)*Eipq(:,:,1,i2))*ktestfield1,idiag_eta12cs)
          if (idiag_eta22ss/=0) call sum_mn_name(-s2z(n)*(-sz(n)*Eipq(:,:,2,i1)+cz(n)*Eipq(:,:,2,i2))*ktestfield1,idiag_eta22ss)
        endif
!
!  Divergence of current helicity flux
!
        if (idiag_s2kzDFm/=0) then
          eetest=etatest*jjtest-uxb
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
        if (idiag_EBpq/=0) call sum_mn_name(cz(n)*Eipq(1,:,1,1) &
                                           +sz(n)*Eipq(1,:,2,1),idiag_EBpq)
!
!  print warning if alpi2 or etai2 or alpi3 are needed, but njtest is too small
!
        if (njtestl<=2 .and. &
            (idiag_alp12/=0.or.idiag_alp22/=0.or.idiag_alp32/=0.or. &
            (leta_rank2.and.(idiag_eta11/=0.or.idiag_eta21/=0)).or. &
            (.not.leta_rank2.and.(idiag_eta12/=0.or.idiag_eta22/=0.or.idiag_eta32/=0)).or. &
            (leta_rank2.and.(idiag_eta11cc/=0.or.idiag_eta21sc/=0)))) then
          call stop_it('njtest is too small if alpi2 or etai2 for i=1,2,3 are needed')
        else if (njtestl < 5 .and. (idiag_alp13/=0 .or. idiag_alp23/=0)) then
          call stop_it('njtest is too small if alpi3 for i=1,2,3 are needed')
        else
!
!  Don't subtract E0 field if iE0==0
!
          if (iE0==0) then
            if (ltestfield_linear) then
              if (idiag_alp12/=0) call sum_mn_name(Eipq(:,:,1,i3),idiag_alp12)
              if (idiag_alp22/=0) call sum_mn_name(Eipq(:,:,2,i3),idiag_alp22)
              if (idiag_alp32/=0) call sum_mn_name(Eipq(:,:,3,i3),idiag_alp32)
              if (idiag_eta11/=0) call sum_mn_name(+(-z(n)*Eipq(:,:,1,i3)+Eipq(:,:,1,i4)),idiag_eta11)
              if (idiag_eta21/=0) call sum_mn_name(+(-z(n)*Eipq(:,:,2,i3)+Eipq(:,:,2,i4)),idiag_eta21)
            else
              if (idiag_alp12/=0) call sum_mn_name(+cz(n)*Eipq(:,:,1,i3)+sz(n)*Eipq(:,:,1,i4),idiag_alp12)
              if (idiag_alp22/=0) call sum_mn_name(+cz(n)*Eipq(:,:,2,i3)+sz(n)*Eipq(:,:,2,i4),idiag_alp22)
              if (idiag_alp32/=0) call sum_mn_name(+cz(n)*Eipq(:,:,3,i3)+sz(n)*Eipq(:,:,3,i4),idiag_alp32)
              if (idiag_alp12cs/=0) call sum_mn_name(csz(n)*(+cz(n)*Eipq(:,:,1,i3)+sz(n)*Eipq(:,:,1,i4)),idiag_alp12cs)
              if (idiag_alp22ss/=0) call sum_mn_name(s2z(n)*(+cz(n)*Eipq(:,:,2,i3)+sz(n)*Eipq(:,:,2,i4)),idiag_alp22ss)
              if (leta_rank2) then
                if (idiag_eta11/=0) call sum_mn_name((-sz(n)*Eipq(:,:,1,i3)+cz(n)*Eipq(:,:,1,i4))*ktestfield1,idiag_eta11)
                if (idiag_eta21/=0) call sum_mn_name((-sz(n)*Eipq(:,:,2,i3)+cz(n)*Eipq(:,:,2,i4))*ktestfield1,idiag_eta21)
                if (idiag_eta11cc/=0) & 
                  call sum_mn_name(c2z(n)*(-sz(n)*Eipq(:,:,1,i3)+cz(n)*Eipq(:,:,1,i4))*ktestfield1,idiag_eta11cc)
                if (idiag_eta21sc/=0) &
                  call sum_mn_name(csz(n)*(-sz(n)*Eipq(:,:,2,i3)+cz(n)*Eipq(:,:,2,i4))*ktestfield1,idiag_eta21sc)
              else
                if (idiag_eta12/=0) call sum_mn_name((-sz(n)*Eipq(:,:,1,i3)+cz(n)*Eipq(:,:,1,i4))*ktestfield1,idiag_eta12)
                if (idiag_eta22/=0) call sum_mn_name((-sz(n)*Eipq(:,:,2,i3)+cz(n)*Eipq(:,:,2,i4))*ktestfield1,idiag_eta22)
                if (idiag_eta32/=0) call sum_mn_name((-sz(n)*Eipq(:,:,3,i3)+cz(n)*Eipq(:,:,3,i4))*ktestfield1,idiag_eta32)
                if (idiag_eta12cs/=0) &
                  call sum_mn_name(csz(n)*(-sz(n)*Eipq(:,:,1,i3)+cz(n)*Eipq(:,:,1,i4))*ktestfield1,idiag_eta12cs)
                if (idiag_eta22ss/=0) &
                  call sum_mn_name(s2z(n)*(-sz(n)*Eipq(:,:,2,i3)+cz(n)*Eipq(:,:,2,i4))*ktestfield1,idiag_eta22ss)
              endif
            endif
          else
            if (idiag_alp12/=0) call sum_mn_name( +cz(n)*(Eipq(:,:,1,i3)-Eipq(:,:,1,iE0)) &
                                                  +sz(n)*(Eipq(:,:,1,i4)-Eipq(:,:,1,iE0)),idiag_alp12)
            if (idiag_alp22/=0) call sum_mn_name( +cz(n)*(Eipq(:,:,2,i3)-Eipq(:,:,2,iE0)) &
                                                  +sz(n)*(Eipq(:,:,2,i4)-Eipq(:,:,2,iE0)),idiag_alp22)
            if (idiag_alp32/=0) call sum_mn_name( +cz(n)*(Eipq(:,:,3,i3)-Eipq(:,:,3,iE0)) &
                                                  +sz(n)*(Eipq(:,:,3,i4)-Eipq(:,:,3,iE0)),idiag_alp32)
            if (idiag_eta11/=0) call sum_mn_name((-sz(n)*(Eipq(:,:,1,i3)-Eipq(:,:,1,iE0)) &
                                                  +cz(n)*(Eipq(:,:,1,i4)-Eipq(:,:,1,iE0)))*ktestfield1,idiag_eta11)
            if (idiag_eta21/=0) call sum_mn_name((-sz(n)*(Eipq(:,:,2,i3)-Eipq(:,:,2,iE0)) &
                                                  +cz(n)*(Eipq(:,:,2,i4)-Eipq(:,:,2,iE0)))*ktestfield1,idiag_eta21)
          endif
        endif
!
!  Volume-averaged dot products of mean emf and velocity and of mean emf and vorticity
!
        if (iE0/=0) then
          if (idiag_E0Um/=0) call sum_mn_name(uxbtestm(nl,1,iE0)*p%uu(:,1) &
                                             +uxbtestm(nl,2,iE0)*p%uu(:,2) &
                                             +uxbtestm(nl,3,iE0)*p%uu(:,3),idiag_E0Um)
          if (idiag_E0Wm/=0) call sum_mn_name(uxbtestm(nl,1,iE0)*p%oo(:,1) &
                                             +uxbtestm(nl,2,iE0)*p%oo(:,2) &
                                             +uxbtestm(nl,3,iE0)*p%oo(:,3),idiag_E0Wm)
        endif
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
          if (idiag_Ex0pt/=0)  call save_name(Eipq(1,lpoint-nghost,1,iE0),idiag_Ex0pt)
          if (idiag_Ex11pt/=0) call save_name(Eipq(1,lpoint-nghost,1,i1),idiag_Ex11pt)
          if (idiag_Ex21pt/=0) call save_name(Eipq(1,lpoint-nghost,1,i2),idiag_Ex21pt)
          if (idiag_Ex12pt/=0) call save_name(Eipq(1,lpoint-nghost,1,i3),idiag_Ex12pt)
          if (idiag_Ex22pt/=0) call save_name(Eipq(1,lpoint-nghost,1,i4),idiag_Ex22pt)
          if (idiag_Ey0pt/=0)  call save_name(Eipq(1,lpoint-nghost,2,iE0),idiag_Ey0pt)
!         if (idiag_bamp/=0)   call save_name(bamp,idiag_bamp)
          if (idiag_Ey11pt/=0) call save_name(Eipq(1,lpoint-nghost,2,i1),idiag_Ey11pt)
          if (idiag_Ey21pt/=0) call save_name(Eipq(1,lpoint-nghost,2,i2),idiag_Ey21pt)
          if (idiag_Ey12pt/=0) call save_name(Eipq(1,lpoint-nghost,2,i3),idiag_Ey12pt)
          if (idiag_Ey22pt/=0) call save_name(Eipq(1,lpoint-nghost,2,i4),idiag_Ey22pt)
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
        if (idiag_b0rms/=0) then
          call dot2(bpq(:,:,iE0),bpq2)
          call sum_mn_name(bpq2,idiag_b0rms,lsqrt=.true.)
        endif
!
        if (idiag_b11rms/=0) then
          call dot2(bpq(:,:,1),bpq2)
          call sum_mn_name(bpq2,idiag_b11rms,lsqrt=.true.)
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
          call dot2(Eipq(1,:,:,iE0),Epq2)
          call sum_mn_name(Epq2,idiag_E0rms,lsqrt=.true.)
        endif
!
        if (idiag_E11rms/=0) then
          call dot2(Eipq(1,:,:,i1),Epq2)
          call sum_mn_name(Epq2,idiag_E11rms,lsqrt=.true.)
        endif
!
        if (idiag_E21rms/=0) then
          call dot2(Eipq(1,:,:,i2),Epq2)
          call sum_mn_name(Epq2,idiag_E21rms,lsqrt=.true.)
        endif
!
        if (idiag_E12rms/=0) then
          call dot2(Eipq(1,:,:,i3),Epq2)
          call sum_mn_name(Epq2,idiag_E12rms,lsqrt=.true.)
        endif
!
        if (idiag_E22rms/=0) then
          call dot2(Eipq(1,:,:,i4),Epq2)
          call sum_mn_name(Epq2,idiag_E22rms,lsqrt=.true.)
        endif
!
      endif
!
!  write B-slices for output in wvid in run.f90
!  Note: ix is the index with respect to array with ghost zones.
!
      if (lvideo.and.lfirst) then
!
        if ('bb11' .in. cnamev) then
!
!  first test solution
!
          do j=1,3
            bb11_yz(m-m1+1,n-n1+1,j)=bpq(ix_loc-l1+1,j,1)
            if (m==iy_loc)  bb11_xz(:,n-n1+1,j)=bpq(:,j,1)
            if (n==iz_loc)  bb11_xy(:,m-m1+1,j)=bpq(:,j,1)
            if (n==iz2_loc) bb11_xy2(:,m-m1+1,j)=bpq(:,j,1)
          enddo
        endif
      endif
!
!      if (lcomplex) deallocate(daatest2,uxb2,bbtest2,duxbbtest2)

    endsubroutine daatest_dt
!***********************************************************************
    subroutine get_slices_testfield(f,slices)
!
!  Write slices for animation of magnetic variables.
!
!  12-sep-09/axel: adapted from the corresponding magnetic routine
!
      use General, only: keep_compiler_quiet
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
          if (slices%index>=3) then
            slices%ready=.false.
          else
            slices%index=slices%index+1
            slices%yz =>bb11_yz(:,:,slices%index)
            slices%xz =>bb11_xz(:,:,slices%index)
            slices%xy =>bb11_xy(:,:,slices%index)
            slices%xy2=>bb11_xy2(:,:,slices%index)
            if (slices%index<=3) slices%ready=.true.
          endif
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
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine testfield_before_boundary
!***********************************************************************
    subroutine testfield_after_boundary(f)
!
!  calculate <uxb>, which is needed when lsoca=.false.
!
!  21-jan-06/axel: coded
!  16-aug-13/MR: MPI communication simplified; changes for iterative procedure
!  20-aug-13/MR: changes for complex calculation: testfield loop then over njtest instead of njtestl
!  27-sep-13/MR: changes due to uxbtestm(mz,...  -->  uxbtestm(nz,...; removed p from parameter list
!  restricted pencil calculation; simplified communication
!
      use Sub
      use Cdata
      use Hydro, only: calc_pencils_hydro
      use Magnetic, only: idiag_bcosphz, idiag_bsinphz
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension (mz) :: c,s
      real, dimension(nx,3,3) :: aijtest,bijtest
      real, dimension(nx,3) :: aatest,bbtest,jjtest,uxbtest,jxbtest
      real, dimension(nx,3) :: del2Atest2,graddivatest
      integer :: jtest,j,juxb,jjxb, njtest_loc,nl
      logical :: headtt_save
      real :: fac, bcosphz, bsinphz, fac1=0., fac2=1.
      type(pencil_case),dimension(:), allocatable :: p          ! vector as scalar quantities not allocatable
      logical, dimension(:), allocatable :: lpenc_loc
!
      uxbtestm=0.; jxbtestm=0.
!
!  In this routine we will reset headtt after the first pencil,
!  so we need to reset it afterwards.
!
      headtt_save=headtt
      fac=1./nxygrid
!
!  do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further up in the file.
!
      if ( .not.lsoca .or. .not.lsoca_jxb .or. liter.and.(t-t_iter_last >= dt_iter) ) then
        if (lcomplex) then
          njtest_loc = njtest
        else
          njtest_loc = njtestl
        endif
      endif
!
!  calculate uxb for nonSOCA or if in iterative procedure integration time for present iteration level 
!  is reached: then uxb needed for rhs of next level
!
      if ( .not.lsoca .or. liter.and.(t-t_iter_last >= dt_iter) ) then

        allocate(p(1),lpenc_loc(npencils))
        lpenc_loc = .false.; lpenc_loc(i_uu)=.true.

        do jtest=1,njtest_loc
!
          iaxtest=iaatest+3*(jtest-1)
          iaztest=iaxtest+2
!
!  prepare coefficients for possible time integration of uxb,
!  But do this only when we are on the last time step
!
          if (ltestfield_taver) then
            if (llast) then
              fac2=1./(nuxb(jtest)+1.)
              fac1=1.-fac2
            endif
          endif
!
!  calculate uxb, and put it into f-array
!  If ltestfield_taver is turned on, the time evolution becomes unstable,
!  so that does not currently work. It was introduced in revision 10236.
!
          do n=n1,n2
            nl=n-n1+1
            do m=m1,m2
!
              call calc_pencils_hydro(f,p(1),lpenc_loc)
              call calc_uxb(f,p(1),iaxtest,uxbtest,bbtest)
              if (lalpha_incoherent) call multsv_mn_add(alpha_tmp,bbtest,uxbtest)
!
              juxb=iuxbtest+3*(jtest-1)
              if (ltestfield_taver) then
                if (llast) then
                  if (iuxbtest/=0) f(l1:l2,m,n,juxb:juxb+2)= &
                          fac1*f(l1:l2,m,n,juxb:juxb+2)+fac2*uxbtest
                endif
              else
                if (iuxbtest/=0) f(l1:l2,m,n,juxb:juxb+2)=uxbtest
              endif
              uxbtestm(nl,:,jtest)=uxbtestm(nl,:,jtest)+fac*sum(uxbtest,1)
              headtt=.false.
            enddo
          enddo
!
!  update counter, but only when we are on the last substep
!
          if (ltestfield_taver) then
            if (llast) then
              nuxb(jtest)=nuxb(jtest)+1
            endif
          endif
!
        enddo
        call finalize_aver(nprocxy,12,uxbtestm)
!
      endif
!
!  Do the same for jxb; do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further up in the file.
!
      if (.not.lsoca_jxb) then
!
!  this not valid for complex calculation.
!
        do jtest=1,njtest_loc
!
          iaxtest=iaatest+3*(jtest-1)
          iaztest=iaxtest+2
!
          do n=n1,n2
            nl=n-n1+1
            do m=m1,m2
              aatest=f(l1:l2,m,n,iaxtest:iaztest)
              call gij(f,iaxtest,aijtest,1)
              call gij_etc(f,iaxtest,aatest,aijtest,bijtest,del2Atest2,graddivatest)
              call curl_mn(aijtest,bbtest,aatest)
              call curl_mn(bijtest,jjtest,bbtest)
              call cross_mn(jjtest,bbtest,jxbtest)
              jjxb=ijxbtest+3*(jtest-1)
              if (ijxbtest/=0) f(l1:l2,m,n,jjxb:jjxb+2)=jxbtest
              jxbtestm(nl,:,jtest)=jxbtestm(nl,:,jtest)+fac*sum(jxbtest,1)
              headtt=.false.
            enddo
          enddo
        enddo
        call finalize_aver(nprocxy,12,jxbtestm)
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
      if (ip<13) print*,'iproc,phase_testfield=',iproc,phase_testfield
!
!  reset headtt
!
      headtt=headtt_save
!
    endsubroutine testfield_after_boundary
!***********************************************************************
    subroutine set_bbtest_Beltrami(B0test,jtest)
!
!  set testfield
!
!  29-mar-08/axel: coded
!
      use Cdata, only: nx,n
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
      use Cdata, only: nx,n
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
      use Cdata, only: nx,n
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
    subroutine set_bbtest_B11_B22(B0test,jtest)
!
!  set testfield
!
!   3-jun-05/axel: coded
!
      use Cdata, only: nx,n
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
      case (3); B0test(:,1)=0.; B0test(:,2)=bamp*cz(n); B0test(:,3)=0.
      case (4); B0test(:,1)=0.; B0test(:,2)=bamp*sz(n); B0test(:,3)=0.
      case default; B0test=0.
      endselect
!
    endsubroutine set_bbtest_B11_B22
!***********************************************************************
    subroutine set_bbtest_B11_B22_lin(B0test,jtest)
!
!  set testfield
!
!  25-Mar-09/axel: adapted from set_bbtest_B11_B22
!  22-Juil-10/emeric: added the fifth testfield to measure alpi3
!                     done only for linear testfields, should also be done for the other cases
!
      use Cdata, only: nx,n,z
!
      real, dimension (nx,3) :: B0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: B0test
!
!  set B0test for each of the 5 cases
!
      select case (jtest)
      case (1); B0test(:,1)=bamp     ; B0test(:,2)=0.; B0test(:,3)=0.
      case (2); B0test(:,1)=bamp*z(n); B0test(:,2)=0.; B0test(:,3)=0.
      case (3); B0test(:,1)=0.; B0test(:,2)=bamp     ; B0test(:,3)=0.
      case (4); B0test(:,1)=0.; B0test(:,2)=bamp*z(n); B0test(:,3)=0.
      case (5); B0test(:,1)=0. ; B0test(:,2)=0 ; B0test(:,3)=bamp
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
!  14-aug-13/MR  : removed unnecessary output of idiags into index.pro
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
        idiag_uzjx1z=0; idiag_uzjy1z=0; idiag_uzjz1z=0
        idiag_uzjx2z=0; idiag_uzjy2z=0; idiag_uzjz2z=0
        idiag_uzjx3z=0; idiag_uzjy3z=0; idiag_uzjz3z=0
        idiag_uzjx4z=0; idiag_uzjy4z=0; idiag_uzjz4z=0
        idiag_E111z=0; idiag_E211z=0; idiag_E311z=0
        idiag_E121z=0; idiag_E221z=0; idiag_E321z=0
        idiag_alp11z=0; idiag_alp21z=0; idiag_alp12z=0; idiag_alp22z=0; idiag_alp13z=0; idiag_alp23z=0
        idiag_eta11z=0; idiag_eta21z=0; idiag_eta12z=0; idiag_eta22z=0
        idiag_E10z=0; idiag_E20z=0; idiag_E30z=0
        idiag_EBpq=0; idiag_E0Um=0; idiag_E0Wm=0
        idiag_alp11=0; idiag_alp21=0; idiag_alp31=0
        idiag_alp12=0; idiag_alp22=0; idiag_alp32=0
        idiag_alp13=0; idiag_alp23=0
        idiag_eta11=0; idiag_eta21=0; idiag_eta31=0
        idiag_eta12=0; idiag_eta22=0; idiag_eta32=0
        idiag_alp11cc=0; idiag_alp21sc=0; idiag_alp12cs=0; idiag_alp22ss=0
        idiag_eta11cc=0; idiag_eta21sc=0; idiag_eta12cs=0; idiag_eta22ss=0
        idiag_s2kzDFm=0
        idiag_M11=0; idiag_M22=0; idiag_M33=0
        idiag_M11cc=0; idiag_M11ss=0; idiag_M22cc=0; idiag_M22ss=0
        idiag_M12cs=0
        idiag_M11z=0; idiag_M22z=0; idiag_M33z=0; idiag_bamp=0
        idiag_jb0m=0; idiag_b0rms=0; idiag_E0rms=0
        idiag_b11rms=0; idiag_b21rms=0; idiag_b12rms=0; idiag_b22rms=0
        idiag_E11rms=0; idiag_E21rms=0; idiag_E12rms=0; idiag_E22rms=0
        idiag_bx0pt=0; idiag_bx11pt=0; idiag_bx21pt=0; idiag_bx12pt=0; idiag_bx22pt=0
        idiag_by0pt=0; idiag_by11pt=0; idiag_by21pt=0; idiag_by12pt=0; idiag_by22pt=0
        idiag_Ex0pt=0; idiag_Ex11pt=0; idiag_Ex21pt=0; idiag_Ex12pt=0; idiag_Ex22pt=0
        idiag_Ey0pt=0; idiag_Ey11pt=0; idiag_Ey21pt=0; idiag_Ey12pt=0; idiag_Ey22pt=0
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
        call parse_name(iname,cname(iname),cform(iname),'alp13',idiag_alp13)
        call parse_name(iname,cname(iname),cform(iname),'alp23',idiag_alp23)
        call parse_name(iname,cname(iname),cform(iname),'eta11',idiag_eta11)
        call parse_name(iname,cname(iname),cform(iname),'eta21',idiag_eta21)
        call parse_name(iname,cname(iname),cform(iname),'eta31',idiag_eta31)
        call parse_name(iname,cname(iname),cform(iname),'eta12',idiag_eta12)
        call parse_name(iname,cname(iname),cform(iname),'eta22',idiag_eta22)
        call parse_name(iname,cname(iname),cform(iname),'eta32',idiag_eta32)
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
        call parse_name(iname,cname(iname),cform(iname),'b11rms',idiag_b11rms)
        call parse_name(iname,cname(iname),cform(iname),'b21rms',idiag_b21rms)
        call parse_name(iname,cname(iname),cform(iname),'b12rms',idiag_b12rms)
        call parse_name(iname,cname(iname),cform(iname),'b22rms',idiag_b22rms)
        call parse_name(iname,cname(iname),cform(iname),'jb0m',idiag_jb0m)
        call parse_name(iname,cname(iname),cform(iname),'bamp',idiag_bamp)
        call parse_name(iname,cname(iname),cform(iname),'b0rms',idiag_b0rms)
        call parse_name(iname,cname(iname),cform(iname),'E11rms',idiag_E11rms)
        call parse_name(iname,cname(iname),cform(iname),'E21rms',idiag_E21rms)
        call parse_name(iname,cname(iname),cform(iname),'E12rms',idiag_E12rms)
        call parse_name(iname,cname(iname),cform(iname),'E22rms',idiag_E22rms)
        call parse_name(iname,cname(iname),cform(iname),'E0rms',idiag_E0rms)
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
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzjx1z',idiag_uzjx1z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzjy1z',idiag_uzjy1z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzjz1z',idiag_uzjz1z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzjx2z',idiag_uzjx2z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzjy2z',idiag_uzjy2z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzjz2z',idiag_uzjz2z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzjx3z',idiag_uzjx3z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzjy3z',idiag_uzjy3z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzjz3z',idiag_uzjz3z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzjx4z',idiag_uzjx4z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzjy4z',idiag_uzjy4z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzjz4z',idiag_uzjz4z)
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
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp11z',idiag_alp11z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp21z',idiag_alp21z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp12z',idiag_alp12z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp22z',idiag_alp22z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp13z',idiag_alp13z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alp23z',idiag_alp23z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta11z',idiag_eta11z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta21z',idiag_eta21z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta12z',idiag_eta12z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eta22z',idiag_eta22z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E10z',idiag_E10z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E20z',idiag_E20z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E30z',idiag_E30z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'M11z',idiag_M11z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'M22z',idiag_M22z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'M33z',idiag_M33z)
      enddo
!
    endsubroutine rprint_testfield

endmodule Testfield
