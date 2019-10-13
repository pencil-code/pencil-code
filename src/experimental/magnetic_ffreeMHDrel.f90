! $Id$

!  Relativistic treatment of force-free magnetic fields.
!  Still quite experimental.
!  This modules deals with all aspects of magnetic fields; if no
!  magnetic fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  magnetically relevant subroutines listed in here.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lmagnetic = .true.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Magnetic
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'magnetic.h'
!
  character (len=labellen) :: initaa='zero',initaa2='zero'

  ! input parameters
  real, dimension(3) :: A0=(/0.,0.,1./),B_ext=(/0.,0.,0./),k_aa=(/0.,0.,0./)
  real, dimension(3) :: axisr1=(/0,0,1/),dispr1=(/0.,0.5,0./)
  real, dimension(3) :: axisr2=(/1,0,0/),dispr2=(/0.,-0.5,0./)
  real, dimension (nx,3) :: bbb
  real :: fring1=0.,Iring1=0.,Rring1=1.,wr1=0.3
  real :: fring2=0.,Iring2=0.,Rring2=1.,wr2=0.3
  real :: amplaa=0.
  real :: radius=.1,epsilonaa=1e-2,widthaa=.5,z0aa=0.
  real :: by_left=0.,by_right=0.
  real :: ABC_A=1.,ABC_B=1.,ABC_C=1.
  real :: amplaa2=0.,kx_aa2=impossible,ky_aa2=impossible,kz_aa2=impossible
  real :: bthresh=0.
  logical :: lpress_equil=.false.
  character (len=40) :: kinflow=''
!
! Slice precalculation buffers
!
  real, target, dimension (:,:,:), allocatable :: bb_xy
  real, target, dimension (:,:,:), allocatable :: bb_xy2
  real, target, dimension (:,:,:), allocatable :: bb_xy3
  real, target, dimension (:,:,:), allocatable :: bb_xy4
  real, target, dimension (:,:,:), allocatable :: bb_xz
  real, target, dimension (:,:,:), allocatable :: bb_yz
  real, target, dimension (:,:,:), allocatable :: bb_xz2

  namelist /magnetic_init_pars/ &
       eta,A0,B_ext,k_aa, &
       fring1,Iring1,Rring1,wr1,axisr1,dispr1, &
       fring2,Iring2,Rring2,wr2,axisr2,dispr2, &
       radius,epsilonaa,z0aa,widthaa,by_left,by_right, &
       initaa,initaa2,amplaa,amplaa2, &
       kx_aa2,ky_aa2,kz_aa2, lpress_equil

  ! run parameters
  real :: eta=0.,B2min=0,height_eta=0.,eta_out=0.
  real :: tau_aa_exterior=0.
  real :: inertial_length=0.,linertial_2
  logical :: lelectron_inertia=.false.

  namelist /magnetic_run_pars/ &
       eta,B_ext,B2min,k_aa, &
       height_eta,eta_out,tau_aa_exterior, &
       kinflow,ABC_A,ABC_B,ABC_C, &
       bthresh

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_dive2m=0,idiag_divee2m=0
  integer :: idiag_e2m=0,idiag_b2m=0,idiag_bm2=0,idiag_j2m=0,idiag_jm2=0
  integer :: idiag_abm=0,idiag_jbm=0,idiag_epsM=0
  integer :: idiag_brms=0,idiag_bmax=0,idiag_jrms=0,idiag_jmax=0
  integer :: idiag_vArms=0,idiag_vAmax=0,idiag_b2mphi=0
  integer :: idiag_bx2m=0, idiag_by2m=0, idiag_bz2m=0
  integer :: idiag_bxmz=0,idiag_bymz=0,idiag_bzmz=0,idiag_bmx=0
  integer :: idiag_bmy=0,idiag_bmz=0,idiag_bxmxy=0,idiag_bymxy=0,idiag_bzmxy=0
  integer :: idiag_uxbm=0,idiag_oxuxbm=0,idiag_jxbxbm=0,idiag_uxDxuxbm=0
  integer :: ivid_bb=0

  contains

!***********************************************************************
    subroutine register_magnetic()
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaa, etc; increase nvar accordingly
!
!  1-may-02/wolf: coded
!
      use FArrayManager
!
      call farray_register_pde('aa',iaa,vector=3)
      iax = iaa; iay = iaa+1; iaz = iaa+2
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
          if (nvar < mvar) write(4,*) ',aa $'
          if (nvar == mvar) write(4,*) ',aa'
        else
          write(4,*) ',aa $'
        endif
        write(15,*) 'aa = fltarr(mx,my,mz,3)*one'
      endif
!
    endsubroutine register_magnetic
!***********************************************************************
    subroutine initialize_magnetic()
!
!  Perform any post-parameter-read initialization
!
!  24-nov-2002/tony: dummy routine - nothing to do at present
      mu01=1./mu0

      if (ivid_bb/=0) then
        !call alloc_slice_buffers(bb_xy,bb_xz,bb_yz,bb_xy2,bb_xy3,bb_xy4,bb_xz2)
        if (lwrite_slice_xy .and..not.allocated(bb_xy) ) allocate(bb_xy (nx,ny,3))
        if (lwrite_slice_xz .and..not.allocated(bb_xz) ) allocate(bb_xz (nx,nz,3))
        if (lwrite_slice_yz .and..not.allocated(bb_yz) ) allocate(bb_yz (ny,nz,3))
        if (lwrite_slice_xy2.and..not.allocated(bb_xy2)) allocate(bb_xy2(nx,ny,3))
        if (lwrite_slice_xy3.and..not.allocated(bb_xy3)) allocate(bb_xy3(nx,ny,3))
        if (lwrite_slice_xy4.and..not.allocated(bb_xy4)) allocate(bb_xy4(nx,ny,3))
        if (lwrite_slice_xz2.and..not.allocated(bb_xz2)) allocate(bb_xz2(nx,nz,3))
      endif

    endsubroutine initialize_magnetic
!***********************************************************************
    subroutine init_aa(f)
!
!  initialise magnetic field; called from start.f90
!  AB: maybe we should here call different routines (such as rings)
!  AB: and others, instead of accummulating all this in a huge routine.
!  We have an init parameter (initaa) to stear magnetic i.c. independently.
!
!   7-nov-2001/wolf: coded
!
      use Mpicomm
      use Gravity, only: gravz
      use Sub
      use Initcond
      use InitialCondition, only: initial_condition_aa
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: cosphi,sinphi
      real, dimension (nx,3) :: bb
      real, dimension (nx) :: b2
      real, dimension (3) :: A0xB0,kxA0,A0xkxA0
      real :: omega_aa,kx_aa,ky_aa,kz_aa
      integer :: j
!
      kx_aa=k_aa(1)
      ky_aa=k_aa(2)
      kz_aa=k_aa(3)
      select case (initaa)

      case ('zero'); f(:,:,:,iax:iaz) = 0.
      case ('wave')
!
!  wave
!
        if (lroot) print*,'init_aa: B_ext=',B_ext
        if (lroot) print*,'init_aa: k_aa=',k_aa
        if (lroot) print*,'init_aa: A0=',A0
!
!  calculate frequency, sin(phi), cos(phi)
!
        !omega_aa=sqrt(kx_aa**2+ky_aa**2+kz_aa**2)
        omega_aa=kx_aa**2+ky_aa**2+kz_aa**2
        do n=n1,n2; do m=m1,m2
          sinphi=amplaa*sin(kx_aa*x(l1:l2)+ky_aa*y(m)+kz_aa*z(n))
          cosphi=amplaa*cos(kx_aa*x(l1:l2)+ky_aa*y(m)+kz_aa*z(n))
!
!  calculate amplitudes
!
          call cross(A0,B_ext,A0xB0)
          call cross(k_aa,A0,kxA0)
          call cross(A0,kxA0,A0xkxA0)
print*,'init_aa: A0=',A0
print*,'init_aa: A0xB0=',A0xB0
print*,'init_aa: A0xkxA0=',A0xkxA0
!
!  Calculate A and S
!
          do j=1,3
            f(l1:l2,m,n,iaa-1+j)=A0(j)*sinphi
            f(l1:l2,m,n,iuu-1+j)=(A0xB0(j)+A0xkxA0(j)*cosphi)*omega_aa*cosphi*omega_aa
          enddo
        enddo; enddo
        !
      case ('current')
!
!  wave
!
        if (lroot) print*,'init_aa: B_ext=',B_ext
        if (lroot) print*,'init_aa: k_aa=',kx_aa,ky_aa,kz_aa
!
        do n=n1,n2; do m=m1,m2
          sinphi=amplaa*sin(kx_aa*x(l1:l2)+ky_aa*y(m)+kz_aa*z(n))
          cosphi=amplaa*cos(kx_aa*x(l1:l2)+ky_aa*y(m)+kz_aa*z(n))
          f(l1:l2,m,n,iay)=cosphi
          f(l1:l2,m,n,iuy)=sinphi**2
        enddo; enddo
!
      case default
!
!  Catch unknown values
!
        if (lroot) print*, 'init_aa: No such value for initaa: ', trim(initaa)
        call stop_it("")

      endselect
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_aa(f)
!
    endsubroutine init_aa
!***********************************************************************
    subroutine daa_dt(f,df,uu,rho1,TT1,uij,bb,va2,shock,gshock)
!
!  solve relativistic force-free MHD equations
!
!  21-jul-03/axel: turned to ffreeMHDrel, adapted from magnetic
!
      use Diagnostics
      use IO
      use Slices_methods, only: store_slices
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: uij
      real, dimension (nx,3) :: aa,jj=0,uxB,JxB,JxBr,oxuxb,jxbxb
      real, dimension (nx,3) :: oo,oxu,uxDxuxb,gshock
      real, dimension (nx) :: rho1,J2,TT1,ab,jb,bx,by,bz,va2,shock
      real, dimension (nx) :: uxb_dotB0,oxuxb_dotB0,jxbxb_dotB0,uxDxuxb_dotB0
      real, dimension (nx) :: bx2, by2, bz2  ! bx^2, by^2 and bz^2
      real :: tmp,eta_out1
      integer :: j,nbthresh
!
      real, dimension (nx,3,3) :: Bij
      real, dimension (nx,3) :: uu,SS,BB,CC,EE
      real, dimension (nx,3) :: divS,curlS,curlB,curlE,del2S,del2A,graddivA
      real, dimension (nx,3) :: SgB,BgS,BdivS,CxE,curlBxB,curlExE,divEE
      real, dimension (nx,3) :: glnrho
      real, dimension (nx) :: B2,B21,E2,divE,divE2,divEE2,ou,o2,sij2
      real, dimension (nx) :: ux,uy,uz,ux2,uy2,uz2
      real, dimension (nx) :: diffus_eta,diffus_nu,advec_va2
      real :: c2=1
!
!
      intent(in)     :: f,uu,rho1,TT1,uij,bb,shock,gshock
      intent(out)    :: bij,va2
      intent(inout)  :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'daa_dt: SOLVE'
      if (headtt) then
        call identify_bcs('Ax',iax)
        call identify_bcs('Ay',iay)
        call identify_bcs('Az',iaz)
      endif
!
!  abbreviations
!
      SS=f(l1:l2,m,n,iux:iuz)
!
!  calculate div and curl of S
!
      call gij(f,iuu,Sij,1)
      call curl_mn(Sij,curlS)
      call trace_mn(Sij,divS)
!
!  calculate del2S only if nu is finite
!
      if (nu/=0.) then
        call del2v(f,iuu,del2S)
      else
        del2S=0.
      endif
!
!  calculate 1/bb^2

      call dot2_mn(BB,B2)
      B21=1./max(B2,B2min)
!
!  "stretching" terms
!
      call multmv_mn(Sij,BB,BgS)
      call multmv_mn(Bij,SS,SgB)
!
!  calculate C = 2*B_k B_k,i
!
      call multmv_mn_transp(Bij,2*BB,CC)
!
!  calculate E = (BxS)/B2, and then E2
!
      call cross_mn(BB,SS,EE)
      call multvs_mn(EE,B21,EE)
      call dot2_mn(EE,E2)
!
!  calculate divE = (S.curlB - B.curlS - C.E)/B2
!
      call dot_mn    (SS,curlB,divE)
      call dot_mn_sub(BB,curlS,divE)
      call dot_mn_sub(CC,EE,divE)
      divE=B21*divE
!
!  calculate curlE = (SgB-BgS+BdivS-CxE)/B2
!
      call multvs_mn(BB,divS,BdivS)
      call cross_mn(CC,EE,CxE)
      call multvs_mn(SgB-BgS+BdivS-CxE,B21,curlE)
!
!  calculate curlB x B, curlE x E, and divE*E
!
      call cross_mn(curlB,BB,curlBxB)
      call cross_mn(curlE,EE,curlExE)
      call multvs_mn(EE,divE,divEE)
!
!  advance Poynting flux and magnetic field
!
      df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+nu*del2S+curlBxB+curlExE+divEE
      df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+eta*del2A-EE
!
!  ``va^2/dx^2'' and ``eta/dx^2'' for timestep
!
      if (lfirst.and.ldt) then
        advec_va2=B2*dxyz_2
        advec2=advec2+advec_va2

        diffus_nu=nu*dxyz_2   ! isn't this done elsewhere ?
        diffus_eta=eta*dxyz_2
        maxdiffus=max(maxdiffus,diffus_eta,diffus_nu)
!
        if (headtt.or.ldebug) then
          print*,'daa_dt: max(advec_va2) =',maxval(advec_va2)
          print*,'daa_dt: max(diffus_nu) =',maxval(diffus_nu)
          print*,'daa_dt: max(diffus_eta) =',maxval(diffus_eta)
        endif
      endif
!
!  calculate B-field, and then max and mean (w/o imposed field, if any)
!
      if (ldiagnos) then
        aa=f(l1:l2,m,n,iax:iaz)
        call dot_mn(aa,bbb,ab)
        call dot2_mn(bbb,b2) !(squared field w/o imposed field)
        call dot2_mn(divEE,divEE2)
        if (idiag_abm/=0)    call sum_mn_name(ab,idiag_abm)
        if (idiag_b2m/=0)    call sum_mn_name(b2,idiag_b2m)
        if (idiag_e2m/=0)    call sum_mn_name(e2,idiag_e2m)
        if (idiag_dive2m/=0) call sum_mn_name(dive**2,idiag_dive2m)
        if (idiag_divee2m/=0) call sum_mn_name(divee2,idiag_divee2m)
        if (idiag_b2mphi/=0) call phisum_mn_name_rz(b2,idiag_b2mphi)
        if (idiag_bm2/=0) call max_mn_name(b2,idiag_bm2)
        if (idiag_brms/=0) call sum_mn_name(b2,idiag_brms,lsqrt=.true.)
        if (idiag_bmax/=0) call max_mn_name(b2,idiag_bmax,lsqrt=.true.)
        if (idiag_bx2m/=0) then
           bx2 = bb(:,1)*bb(:,1)
           call sum_mn_name(bx2,idiag_bx2m)
        endif
        if (idiag_by2m/=0) then
           by2 = bb(:,2)*bb(:,2)
           call sum_mn_name(by2,idiag_by2m)
        endif
        if (idiag_bz2m/=0) then
           bz2 = bb(:,3)*bb(:,3)
           call sum_mn_name(bz2,idiag_bz2m)
        endif
      endif
!
!  1-D averages.
!
      if (l1ddiagnos) then
        bx=bb(:,1); by=bb(:,2); bz=bb(:,3)
        if (idiag_bxmz/=0) call xysum_mn_name_z(bx,idiag_bxmz)
        if (idiag_bymz/=0) call xysum_mn_name_z(by,idiag_bymz)
        if (idiag_bzmz/=0) call xysum_mn_name_z(bz,idiag_bzmz)
      endif
!
!  2-D averages.
!
      if (l2davgfirst) then
        bx=bb(:,1); by=bb(:,2); bz=bb(:,3)
        if (idiag_bxmxy/=0) call zsum_mn_name_xy(bx,idiag_bxmxy)
        if (idiag_bymxy/=0) call zsum_mn_name_xy(by,idiag_bymxy)
        if (idiag_bzmxy/=0) call zsum_mn_name_xy(bz,idiag_bzmxy)
      endif
!
!  write B-slices for output in wvid in run.f90
!  Note: ix is the index with respect to array with ghost zones.
!
        if (lvideo.and.lfirst.and.ivid_bb) then
          call store_slices(p%bb,bb_xy,bb_xz,bb_yz,bb_xy2,bb_xy3,bb_xy4,bb_xz2)
          call vecout(41,trim(directory_snap)//'/bvec.dat',bbb,bthresh,nbthresh)
        endif
!
!  set alven speed to zero for proper time steppeing
!
      va2=0
!
!  calculate max and rms current density
!  at the moment (and in future?) calculate max(b^2) and mean(b^2).
!
      if (ldiagnos) then
        !
        !  v_A = |B|/sqrt(rho); in units where "4pi"=1
        !
        if (idiag_vArms/=0) call sum_mn_name(va2,idiag_vArms,lsqrt=.true.)
        if (idiag_vAmax/=0) call max_mn_name(va2,idiag_vAmax,lsqrt=.true.)
        !
        ! <J.B>
        !
        if (idiag_jbm/=0) then
          call dot_mn(jj,bb,jb)
          call sum_mn_name(jb,idiag_jbm)
        endif
        !
        ! <J^2> and J^2|max
        !
        if (idiag_jrms/=0 .or. idiag_jmax/=0 .or. idiag_j2m/=0 &
            .or.idiag_jm2/=0.or.idiag_epsM/=0) then
          call dot2_mn(jj,j2)
          if (idiag_j2m/=0) call sum_mn_name(j2,idiag_j2m)
          if (idiag_jm2/=0) call max_mn_name(j2,idiag_jm2)
          if (idiag_jrms/=0) call sum_mn_name(j2,idiag_jrms,lsqrt=.true.)
          if (idiag_jmax/=0) call max_mn_name(j2,idiag_jmax,lsqrt=.true.)
          if (idiag_epsM/=0) call sum_mn_name(eta*j2,idiag_epsM)
        endif
        !
        if (idiag_uxbm/=0) then
          call cross_mn(uu,bbb,uxb)
          uxb_dotB0=B_ext(1)*uxb(:,1)+B_ext(2)*uxb(:,2)+B_ext(3)*uxb(:,3)
          uxb_dotB0=uxb_dotB0/(B_ext(1)**2+B_ext(2)**2+B_ext(3)**2)
          call sum_mn_name(uxb_dotB0,idiag_uxbm)
        endif
        !
        if (idiag_jxbxbm/=0) then
          call cross_mn(jj,bbb,jxb)
          call cross_mn(jxb,bbb,jxbxb)
          jxbxb_dotB0=B_ext(1)*jxbxb(:,1)+B_ext(2)*jxbxb(:,2)+B_ext(3)*jxbxb(:,3)
          jxbxb_dotB0=jxbxb_dotB0/(B_ext(1)**2+B_ext(2)**2+B_ext(3)**2)
          call sum_mn_name(jxbxb_dotB0,idiag_jxbxbm)
        endif
        !
        if (idiag_oxuxbm/=0) then
          oo(:,1)=uij(:,3,2)-uij(:,2,3)
          oo(:,2)=uij(:,1,3)-uij(:,3,1)
          oo(:,3)=uij(:,2,1)-uij(:,1,2)
          call cross_mn(oo,uu,oxu)
          call cross_mn(oxu,bbb,oxuxb)
          oxuxb_dotB0=B_ext(1)*oxuxb(:,1)+B_ext(2)*oxuxb(:,2)+B_ext(3)*oxuxb(:,3)
          oxuxb_dotB0=oxuxb_dotB0/(B_ext(1)**2+B_ext(2)**2+B_ext(3)**2)
          call sum_mn_name(oxuxb_dotB0,idiag_oxuxbm)
        endif
        !
        !  < u x curl(uxB) > = < E_i u_{j,j} - E_j u_{j,i} >
        !   ( < E_1 u2,2 + E1 u3,3 - E2 u2,1 - E3 u3,1 >
        !   ( < E_2 u1,1 + E2 u3,3 - E1 u2,1 - E3 u3,2 >
        !   ( < E_3 u1,1 + E3 u2,2 - E1 u3,1 - E2 u2,3 >
        !
        if (idiag_uxDxuxbm/=0) then
          call cross_mn(uu,bbb,uxb)
          uxDxuxb(:,1)=uxb(:,1)*(uij(:,2,2)+uij(:,3,3))-uxb(:,2)*uij(:,2,1)-uxb(:,3)*uij(:,3,1)
          uxDxuxb(:,2)=uxb(:,2)*(uij(:,1,1)+uij(:,3,3))-uxb(:,1)*uij(:,1,2)-uxb(:,3)*uij(:,3,2)
          uxDxuxb(:,3)=uxb(:,3)*(uij(:,1,1)+uij(:,2,2))-uxb(:,1)*uij(:,1,3)-uxb(:,2)*uij(:,2,3)
          uxDxuxb_dotB0=B_ext(1)*uxDxuxb(:,1)+B_ext(2)*uxDxuxb(:,2)+B_ext(3)*uxDxuxb(:,3)
          uxDxuxb_dotB0=uxDxuxb_dotB0/(B_ext(1)**2+B_ext(2)**2+B_ext(3)**2)
          call sum_mn_name(uxDxuxb_dotB0,idiag_uxDxuxbm)
        endif
        !
      endif
!
!  debug output
!
      if (headtt .and. lfirst .and. ip<=4) then
        call output_pencil('aa.dat',aa,3)
        call output_pencil('bb.dat',bb,3)
        call output_pencil('jj.dat',jj,3)
        call output_pencil('del2A.dat',del2A,3)
        call output_pencil('JxBr.dat',JxBr,3)
        call output_pencil('JxB.dat',JxB,3)
        call output_pencil('df.dat',df(l1:l2,m,n,:),mvar)
      endif
!
if (ip<3.and.m==4.and.n==4) write(61) ss,Sij,curlS,divS,del2A,curlB
if (ip<3.and.m==4.and.n==4) write(61) BB,B2,BgS,SgB,Bij,CC,EE,B21
if (ip<3.and.m==4.and.n==4) write(61) divE,BdivS,CxE,curlBxB,curlE,curlExE,divEE
!
      call keep_compiler_quiet(shock,gshock(:,1))
!
    endsubroutine daa_dt
!***********************************************************************
    subroutine calculate_vars_magnetic(f,bb,bij,aij,del2A,graddivA)
!
!  calculate bb
!  possibility to add external field
!  Note; for diagnostics purposes keep copy bbb of original field
!
!  06-feb-04/bing: coded
!  27-mar-07/axel: this routine needs to be replaced by calc_pencils_magnetic
!
      use Sub

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3,3) :: bij,aij
      real, dimension (nx,3) :: bb,del2A,graddivA

      intent(in)  :: f
      intent(out) :: bb,bij,del2A,graddivA

      call gij(f,iaa,aij,1)
      call curl_mn(aij,bb,aa)
      call gij_etc(f,iaa,aa,aij,bij,del2A,graddivA)
!
!  possibility to add external field
!  Note; for diagnostics purposes keep copy of original field
!
       if (ldiagnos) bbb=bb
       if (B_ext(1)/=0.) bb(:,1)=bb(:,1)+B_ext(1)
       if (B_ext(2)/=0.) bb(:,2)=bb(:,2)+B_ext(2)
       if (B_ext(3)/=0.) bb(:,3)=bb(:,3)+B_ext(3)
       if (headtt) print*,'calculate_vars_magnetic: B_ext=',B_ext

    endsubroutine calculate_vars_magnetic
!***********************************************************************
    subroutine read_magnetic_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=magnetic_init_pars, IOSTAT=iostat)
!
    endsubroutine read_magnetic_init_pars
!***********************************************************************
    subroutine write_magnetic_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=magnetic_init_pars)
!
    endsubroutine write_magnetic_init_pars
!***********************************************************************
    subroutine read_magnetic_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=magnetic_run_pars, IOSTAT=iostat)
!
    endsubroutine read_magnetic_run_pars
!***********************************************************************
    subroutine write_magnetic_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=magnetic_run_pars)
!
    endsubroutine write_magnetic_run_pars
!***********************************************************************
    subroutine input_persistent_magnetic(id,done)
!
      integer, intent(in), optional :: id
      logical, intent(inout), optional :: done
!
      if (present (id)) call keep_compiler_quiet(id)
      if (present (done)) call keep_compiler_quiet(done)
!
    endsubroutine input_persistent_magnetic
!***********************************************************************
    subroutine rprint_magnetic(lreset,lwrite)
!
!  reads and registers print parameters relevant for magnetic fields
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
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
        idiag_dive2m=0; idiag_divee2m=0
        idiag_e2m=0; idiag_b2m=0; idiag_bm2=0; idiag_j2m=0; idiag_jm2=0
        idiag_abm=0; idiag_jbm=0; idiag_epsM=0; idiag_b2mphi=0
        idiag_brms=0; idiag_bmax=0; idiag_jrms=0; idiag_jmax=0; idiag_vArms=0
        idiag_vAmax=0; idiag_bx2m=0; idiag_by2m=0; idiag_bz2m=0
        idiag_bxmz=0; idiag_bymz=0; idiag_bzmz=0; idiag_bmx=0; idiag_bmy=0
        idiag_bmz=0; idiag_bxmxy=0; idiag_bymxy=0; idiag_bzmxy=0
        idiag_uxbm=0; idiag_oxuxbm=0; idiag_jxbxbm=0.; idiag_uxDxuxbm=0.
        ivid_bb=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'abm',idiag_abm)
        call parse_name(iname,cname(iname),cform(iname),'jbm',idiag_jbm)
        call parse_name(iname,cname(iname),cform(iname),'b2m',idiag_b2m)
        call parse_name(iname,cname(iname),cform(iname),'e2m',idiag_e2m)
        call parse_name(iname,cname(iname),cform(iname),'dive2m',idiag_dive2m)
        call parse_name(iname,cname(iname),cform(iname),'divee2m',idiag_divee2m)
        call parse_name(iname,cname(iname),cform(iname),'bm2',idiag_bm2)
        call parse_name(iname,cname(iname),cform(iname),'j2m',idiag_j2m)
        call parse_name(iname,cname(iname),cform(iname),'jm2',idiag_jm2)
        call parse_name(iname,cname(iname),cform(iname),'epsM',idiag_epsM)
        call parse_name(iname,cname(iname),cform(iname),'brms',idiag_brms)
        call parse_name(iname,cname(iname),cform(iname),'bmax',idiag_bmax)
        call parse_name(iname,cname(iname),cform(iname),'jrms',idiag_jrms)
        call parse_name(iname,cname(iname),cform(iname),'jmax',idiag_jmax)
        call parse_name(iname,cname(iname),cform(iname),'vArms',idiag_vArms)
        call parse_name(iname,cname(iname),cform(iname),'vAmax',idiag_vAmax)
        call parse_name(iname,cname(iname),cform(iname),'bx2m',idiag_bx2m)
        call parse_name(iname,cname(iname),cform(iname),'by2m',idiag_by2m)
        call parse_name(iname,cname(iname),cform(iname),'bz2m',idiag_bz2m)
        call parse_name(iname,cname(iname),cform(iname),'uxbm',idiag_uxbm)
        call parse_name(iname,cname(iname),cform(iname),'jxbxbm',idiag_jxbxbm)
        call parse_name(iname,cname(iname),cform(iname),'oxuxbm',idiag_oxuxbm)
        call parse_name(iname,cname(iname),cform(iname),&
            'uxDxuxbm',idiag_uxDxuxbm)
        call parse_name(iname,cname(iname),cform(iname),'bmx',idiag_bmx)
        call parse_name(iname,cname(iname),cform(iname),'bmy',idiag_bmy)
        call parse_name(iname,cname(iname),cform(iname),'bmz',idiag_bmz)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bxmz',idiag_bxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bymz',idiag_bymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bzmz',idiag_bzmz)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do ixy=1,nnamexy
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bxmxy',idiag_bxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bymxy',idiag_bymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bzmxy',idiag_bzmxy)
      enddo
!
! Currently need to force zaverage calculation at every lout step for
! bmx and bmy.
!
      if ((idiag_bmx+idiag_bmy)>0) ldiagnos_need_zaverages=.true.
!
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'b2mphi',idiag_b2mphi)
      enddo
!
!  check for those quantities for which we want video slices
!
      do iname=1,nnamev
        call parse_name(iname,cnamev(iname),cformv(iname),'bb',ivid_bb)
      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!
      if (lwr) then
        call farray_index_append('i_abm',idiag_abm)
        call farray_index_append('i_jbm',idiag_jbm)
        call farray_index_append('i_b2m',idiag_b2m)
        call farray_index_append('i_e2m',idiag_e2m)
        call farray_index_append('i_dive2m',idiag_dive2m)
        call farray_index_append('i_divee2m',idiag_divee2m)
        call farray_index_append('i_bm2',idiag_bm2)
        call farray_index_append('i_j2m',idiag_j2m)
        call farray_index_append('i_jm2',idiag_jm2)
        call farray_index_append('i_epsM',idiag_epsM)
        call farray_index_append('i_brms',idiag_brms)
        call farray_index_append('i_bmax',idiag_bmax)
        call farray_index_append('i_jrms',idiag_jrms)
        call farray_index_append('i_jmax',idiag_jmax)
        call farray_index_append('i_vArms',idiag_vArms)
        call farray_index_append('i_vAmax',idiag_vAmax)
        call farray_index_append('i_bx2m',idiag_bx2m)
        call farray_index_append('i_by2m',idiag_by2m)
        call farray_index_append('i_bz2m',idiag_bz2m)
        call farray_index_append('i_uxbm',idiag_uxbm)
        call farray_index_append('i_oxuxbm',idiag_oxuxbm)
        call farray_index_append('i_jxbxbm',idiag_jxbxbm)
        call farray_index_append('i_uxDxuxbm',idiag_uxDxuxbm)
        call farray_index_append('i_bxmz',idiag_bxmz)
        call farray_index_append('i_bymz',idiag_bymz)
        call farray_index_append('i_bzmz',idiag_bzmz)
        call farray_index_append('i_bmx',idiag_bmx)
        call farray_index_append('i_bmy',idiag_bmy)
        call farray_index_append('i_bmz',idiag_bmz)
        call farray_index_append('i_bxmxy',idiag_bxmxy)
        call farray_index_append('i_bymxy',idiag_bymxy)
        call farray_index_append('i_bzmxy',idiag_bzmxy)
        call farray_index_append('i_b2mphi',idiag_b2mphi)
        call farray_index_append('nname',nname)
        call farray_index_append('nnamexy',nnamexy)
        call farray_index_append('nnamez',nnamez)
        call farray_index_append('iaa',iaa)
        call farray_index_append('iax',iax)
        call farray_index_append('iay',iay)
        call farray_index_append('iaz',iaz)
      endif
!
    endsubroutine rprint_magnetic
!***********************************************************************
    subroutine calc_mfield
!
!  calculate mean magnetic field from xy- or z-averages
!
!  19-jun-02/axel: moved from print to here
!   9-nov-02/axel: corrected bxmy(m,j); it used bzmy instead!
!
      use Sub
!
      logical,save :: first=.true.
      real, dimension(nx) :: bymx,bzmx
      real, dimension(ny,nprocy) :: bxmy,bzmy
      real :: bmx,bmy,bmz
      integer :: l,j
!
!  For vector output (of bb vectors) we need brms
!  on all processors. It suffices to have this for times when lout=.true.,
!  but we need to broadcast the result to all procs.
!
!  calculate brms (this requires that brms is set in print.in)
!  broadcast result to other processors
!
      if (idiag_brms/=0) then
        if (iproc==0) brms=fname(idiag_brms)
        call mpibcast_real(brms)
      endif

      if (.not.lroot) return
!
!  Magnetic energy in vertically averaged field
!  The bymxy and bzmxy must have been calculated,
!  so they are present on the root processor.
!
      if (idiag_bmx/=0) then
        if (idiag_bymxy==0.or.idiag_bzmxy==0) then
          if (first) print*,"calc_mfield:                  WARNING"
          if (first) print*, &
                  "calc_mfield: NOTE: to get bmx, bymxy and bzmxy must also be set in zaver"
          if (first) print*, &
                  "calc_mfield:       We proceed, but you'll get bmx=0"
          bmx=0.
        else
          do l=1,nx
            bymx(l)=sum(fnamexy(l,:,:,idiag_bymxy))/(ny*nprocy)
            bzmx(l)=sum(fnamexy(l,:,:,idiag_bzmxy))/(ny*nprocy)
          enddo
          bmx=sqrt(sum(bymx**2+bzmx**2)/nx)
        endif
        call save_name(bmx,idiag_bmx)
      endif
!
!  similarly for bmy
!
      if (idiag_bmy/=0) then
        if (idiag_bxmxy==0.or.idiag_bzmxy==0) then
          if (first) print*,"calc_mfield:                  WARNING"
          if (first) print*, &
                  "calc_mfield: NOTE: to get bmy, bxmxy and bzmxy must also be set in zaver"
          if (first) print*, &
                  "calc_mfield:       We proceed, but you'll get bmy=0"
          bmy=0.
        else
          do j=1,nprocy
          do m=1,ny
            bxmy(m,j)=sum(fnamexy(:,m,j,idiag_bxmxy))/nx
            bzmy(m,j)=sum(fnamexy(:,m,j,idiag_bzmxy))/nx
          enddo
          enddo
          bmy=sqrt(sum(bxmy**2+bzmy**2)/(ny*nprocy))
        endif
        call save_name(bmy,idiag_bmy)
      endif
!
!  Magnetic energy in horizontally averaged field
!  The bxmz and bymz must have been calculated,
!  so they are present on the root processor.
!
      if (idiag_bmz/=0) then
        if (idiag_bxmz==0.or.idiag_bymz==0) then
          if (first) print*,"calc_mfield:                  WARNING"
          if (first) print*, &
                  "calc_mfield: NOTE: to get bmz, bxmz and bymz must also be set in xyaver"
          if (first) print*, &
                  "calc_mfield:       This may be because we renamed zaver.in into xyaver.in"
          if (first) print*, &
                  "calc_mfield:       We proceed, but you'll get bmz=0"
          bmz=0.
        else
          bmz=sqrt(sum(fnamez(:,:,idiag_bxmz)**2+fnamez(:,:,idiag_bymz)**2)/(nz*nprocz))
        endif
        call save_name(bmz,idiag_bmz)
      endif
!
      first = .false.
!
    endsubroutine calc_mfield
!***********************************************************************
    subroutine bc_frozen_in_bb_z(topbot)
!
!  Dummy routine for frozen-in flux at boundary
!
      character (len=3) :: topbot
!

      print*, 'WARNING:'
      print*, '  bc_frozen_in_bb_z not implemented for magnetic_ffreeMHDrel !!'
!
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_frozen_in_bb_z
!***********************************************************************
    subroutine bc_frozen_in_bb_z(topbot)
!
!  Set flags to indicate that magnetic flux is frozen-in at the
!  z boundary. The implementation occurs in daa_dt where magnetic
!  diffusion is switched off in that layer.
!
      character (len=3) :: topbot
!
      select case (topbot)
      case ('bot')               ! bottom boundary
        lfrozen_bz_z_bot = .true.    ! set flag
      case ('top')               ! top boundary
        lfrozen_bz_z_top = .true.    ! set flag
      case default
        print*, "bc_frozen_in_bb_z: ", topbot, " should be `top' or `bot'"
      endselect
!
    endsubroutine bc_frozen_in_bb_z
!***********************************************************************
      subroutine bc_aa_pot(f,topbot)
!
!  Potential field boundary condition for magnetic vector potential
!
!  14-jun-2002/axel: adapted from similar
!   8-jul-2002/axel: introduced topbot argument
!
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny) :: f2,f3
      real, dimension (nx,ny,nghost+1) :: fz
      integer :: j
!
!  potential field condition
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  potential field condition at the bottom
!
      case ('bot')
        if (headtt) print*,'bc_aa_pot: potential field boundary condition at the bottom'
        if (nprocy/=1) &
             call stop_it("bc_aa_pot: potential field doesn't work yet with nprocy/=1")
        do j=0,1
          f2=f(l1:l2,m1:m2,n1+1,iax+j)
          f3=f(l1:l2,m1:m2,n1+2,iax+j)
          call potential_field(fz,f2,f3,-1)
          f(l1:l2,m1:m2,1:n1,iax+j)=fz
        enddo
        !
        f2=f(l1:l2,m1:m2,n1,iax)
        f3=f(l1:l2,m1:m2,n1,iay)
        call potentdiv(fz,f2,f3,-1)
        f(l1:l2,m1:m2,1:n1,iaz)=-fz
!
!  potential field condition at the top
!
      case ('top')
        if (headtt) print*,'bc_aa_pot: potential field boundary condition at the top'
        if (nprocy/=1) &
             call stop_it("bc_aa_pot: potential field doesn't work yet with nprocy/=1")
        do j=0,1
          f2=f(l1:l2,m1:m2,n2-1,iax+j)
          f3=f(l1:l2,m1:m2,n2-2,iax+j)
          call potential_field(fz,f2,f3,+1)
          f(l1:l2,m1:m2,n2:mz,iax+j)=fz
        enddo
        !
        f2=f(l1:l2,m1:m2,n2,iax)
        f3=f(l1:l2,m1:m2,n2,iay)
        call potentdiv(fz,f2,f3,+1)
        f(l1:l2,m1:m2,n2:mz,iaz)=-fz
      case default
        if (lroot) print*,"bc_aa_pot: invalid argument"
      endselect
!
      endsubroutine bc_aa_pot
!***********************************************************************
      subroutine potential_field(fz,f2,f3,irev)
!
!  solves the potential field boundary condition;
!  fz is the boundary layer, and f2 and f3 are the next layers inwards.
!  The condition is the same on the two sides.
!
!  20-jan-00/axel+wolf: coded
!  22-mar-00/axel: corrected sign (it is the same on both sides)
!
      use Fourier
!
      real, dimension (nx,ny) :: fac,kk,f1r,f1i,g1r,g1i,f2,f2r,f2i,f3,f3r,f3i
      real, dimension (nx,ny,nghost+1) :: fz
      real, dimension (nx) :: kx
      real, dimension (ny) :: ky
      real :: delz
      integer :: i,irev
!
      f2r=f2; f2i=0
      f3r=f3; f3i=0
!
!  Transform
!
      call fourier_transform(f2r,f2i) ! x-direction
      call fourier_transform(f2r,f2i) ! y-direction
!
      call fourier_transform(f3r,f3i) ! x-direction
      call fourier_transform(f3r,f3i) ! y-direction
!
!  define wave vector
!
      kx=cshift((/(i-(nx-1)/2,i=0,nx-1)/),+(nx-1)/2)*2*pi/Lx
      ky=cshift((/(i-(ny-1)/2,i=0,ny-1)/),+(ny-1)/2)*2*pi/Ly
!
!  calculate 1/k^2, zero mean
!
      kk=sqrt(spread(kx**2,2,ny)+spread(ky**2,1,nx))
!
!  one-sided derivative
!
      fac=1./(3.+2.*dz*kk)
      f1r=fac*(4.*f2r-f3r)
      f1i=fac*(4.*f2i-f3i)
!
!  set ghost zones
!
      do i=0,nghost
        delz=i*dz
        fac=exp(-kk*delz)
        g1r=fac*f1r
        g1i=fac*f1i
!
!  Transform back
!
        call fourier_transform(g1r,g1i,linv=.true.) ! x-direction
        call fourier_transform(g1r,g1i,linv=.true.) ! y-direction
!
!  reverse order if irev=-1 (if we are at the bottom)
!
        if (irev==+1) fz(:,:,       i+1) = g1r/(nx*ny)  ! Renormalize
        if (irev==-1) fz(:,:,nghost-i+1) = g1r/(nx*ny)  ! Renormalize
      enddo
!
    endsubroutine potential_field
!***********************************************************************
      subroutine potentdiv(fz,f2,f3,irev)
!
!  solves the divA=0 for potential field boundary condition;
!  f2 and f3 correspond to Ax and Ay (input) and fz corresponds to Ax (out)
!  In principle we could save some ffts, by combining with the potential
!  subroutine above, but this is now easier
!
!  22-mar-02/axel: coded
!
      use Fourier
!
      real, dimension (nx,ny) :: fac,kk,kkkx,kkky,f1r,f1i,g1r,g1i,f2,f2r,f2i,f3,f3r,f3i
      real, dimension (nx,ny,nghost+1) :: fz
      real, dimension (nx) :: kx
      real, dimension (ny) :: ky
      real :: delz
      integer :: i,irev
!
      f2r=f2; f2i=0
      f3r=f3; f3i=0
!
!  Transform
!
      call fourier_transform(f2r,f2i) ! x-direction
      call fourier_transform(f2r,f2i) ! y-direction
!
      call fourier_transform(f3r,f3i) ! x-direction
      call fourier_transform(f3r,f3i) ! y-direction
!
!  define wave vector
!
      kx=cshift((/(i-nx/2,i=0,nx-1)/),+nx/2)
      ky=cshift((/(i-ny/2,i=0,ny-1)/),+ny/2)
!
!  calculate 1/k^2, zero mean
!
      kk=sqrt(spread(kx**2,2,ny)+spread(ky**2,1,nx))
      kkkx=spread(kx,2,ny)
      kkky=spread(ky,1,nx)
!
!  calculate 1/kk
!
      kk(1,1)=1.
      fac=1./kk
      fac(1,1)=0.
!
      f1r=fac*(-kkkx*f2i-kkky*f3i)
      f1i=fac*(+kkkx*f2r+kkky*f3r)
!
!  set ghost zones
!
      do i=0,nghost
        delz=i*dz
        fac=exp(-kk*delz)
        g1r=fac*f1r
        g1i=fac*f1i
!
!  Transform back
!
        call fourier_transform(g1r,g1i,linv=.true.) ! x-direction
        call fourier_transform(g1r,g1i,linv=.true.) ! y-direction
!
!  reverse order if irev=-1 (if we are at the bottom)
!
        if (irev==+1) fz(:,:,       i+1) = g1r/(nx*ny)  ! Renormalize
        if (irev==-1) fz(:,:,nghost-i+1) = g1r/(nx*ny)  ! Renormalize
      enddo
!
    endsubroutine potentdiv
!***********************************************************************
    subroutine bb_unitvec_shock(f,bb_hat)
!
!  Dummy routine
!
!  18-aug-2006/tobi: coded
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,3), intent (out) :: bb_hat
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(bb_hat)
!
    endsubroutine bb_unitvec_shock
!***********************************************************************
    subroutine split_update_magnetic(f)
!
!  dummy
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine split_update_magnetic
!***********************************************************************
endmodule Magnetic
