! $Id: magnetic.f90,v 1.150 2003-11-23 21:59:37 brandenb Exp $

!  This modules deals with all aspects of magnetic fields; if no
!  magnetic fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  magnetically relevant subroutines listed in here.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Magnetic

  use Cparam

  implicit none

  character (len=labellen) :: initaa='zero',initaa2='zero'
  character (len=labellen) :: ires='eta-const'
  ! input parameters
  real, dimension(3) :: axisr1=(/0,0,1/),dispr1=(/0.,0.5,0./)
  real, dimension(3) :: axisr2=(/1,0,0/),dispr2=(/0.,-0.5,0./)
  real :: fring1=0.,Iring1=0.,Rring1=1.,wr1=0.3
  real :: fring2=0.,Iring2=0.,Rring2=1.,wr2=0.3
  real :: amplaa=0., kx_aa=1.,ky_aa=1.,kz_aa=1.
  real :: radius=.1,epsilonaa=1e-2,widthaa=.5,z0aa=0.
  real :: by_left=0.,by_right=0.
  real :: ABC_A=1.,ABC_B=1.,ABC_C=1.
  real :: amplaa2=0.,kx_aa2=impossible,ky_aa2=impossible,kz_aa2=impossible
  real :: bthresh=0.,bthresh_per_brms=0.,brms=0.,bthresh_scl=1.
  integer :: nbvec,nbvecmax=nx*ny*nz/4
  logical :: lpress_equil=.false.
  character (len=40) :: kinflow=''
  real :: alpha_effect
  complex, dimension(3) :: coefaa=(/0.,0.,0./), coefbb=(/0.,0.,0./)

  namelist /magnetic_init_pars/ &
       fring1,Iring1,Rring1,wr1,axisr1,dispr1, &
       fring2,Iring2,Rring2,wr2,axisr2,dispr2, &
       radius,epsilonaa,z0aa,widthaa,by_left,by_right, &
       initaa,initaa2,amplaa,amplaa2,kx_aa,ky_aa,kz_aa,coefaa,coefbb, &
       kx_aa2,ky_aa2,kz_aa2, lpress_equil

  ! run parameters
  real, dimension(3) :: B_ext=(/0.,0.,0./)
  real :: eta=0.,height_eta=0.,eta_out=0.
  real :: tau_aa_exterior=0.

  namelist /magnetic_run_pars/ &
       eta,B_ext,alpha_effect, &
       height_eta,eta_out,tau_aa_exterior, &
       kinflow,kx_aa,ky_aa,kz_aa,ABC_A,ABC_B,ABC_C, &
       bthresh,bthresh_per_brms,ires

  ! other variables (needs to be consistent with reset list below)
  integer :: i_b2m=0,i_bm2=0,i_j2m=0,i_jm2=0,i_abm=0,i_jbm=0,i_ubm,i_epsM=0
  integer :: i_bxpt=0,i_bypt=0,i_bzpt=0
  integer :: i_aybym2=0,i_exaym2=0
  integer :: i_brms=0,i_bmax=0,i_jrms=0,i_jmax=0,i_vArms=0,i_vAmax=0
  integer :: i_bx2m=0, i_by2m=0, i_bz2m=0
  integer :: i_bxbym=0, i_bxbzm=0, i_bybzm=0
  integer :: i_bxmz=0,i_bymz=0,i_bzmz=0,i_bmx=0,i_bmy=0,i_bmz=0
  integer :: i_bxmxy=0,i_bymxy=0,i_bzmxy=0
  integer :: i_uxbm=0,i_oxuxbm=0,i_jxbxbm=0,i_gpxbm=0,i_uxDxuxbm=0
  integer :: i_uxbmx=0,i_uxbmy=0,i_uxbmz=0,i_uxjm=0,i_ujxbm
  integer :: i_brmphi=0,i_bpmphi=0,i_bzmphi=0,i_b2mphi=0,i_jbmphi=0

  contains

!***********************************************************************
    subroutine register_magnetic()
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaa, etc; increase nvar accordingly
!
!  1-may-02/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_aa called twice')
      first = .false.
!
      lmagnetic = .true.
      iaa = nvar+1              ! indices to access aa
      iax = iaa
      iay = iaa+1
      iaz = iaa+2
      nvar = nvar+3             ! added 3 variables
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_magnetic: nvar = ', nvar
        print*, 'register_magnetic: iaa,iax,iay,iaz = ', iaa,iax,iay,iaz
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: magnetic.f90,v 1.150 2003-11-23 21:59:37 brandenb Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_magnetic: nvar > mvar')
      endif
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

    endsubroutine initialize_magnetic
!***********************************************************************
    subroutine init_aa(f,xx,yy,zz)
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
      use Density
      use Gravity, only: gravz
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
      real, dimension (nx,3) :: bb
      real, dimension (nx) :: b2
!
      select case(initaa)

      case('zero', '0'); f(:,:,:,iax:iaz) = 0.
      case('mode'); call modev(amplaa,coefaa,f,iaa,kx_aa,ky_aa,kz_aa,xx,yy,zz)
      case('modeb'); call modeb(amplaa,coefbb,f,iaa,kx_aa,ky_aa,kz_aa,xx,yy,zz)
      case('gaussian-noise'); call gaunoise(amplaa,f,iax,iaz)
      case('Beltrami-x', '11'); call beltrami(amplaa,f,iaa,KX=kx_aa)
      case('Beltrami-y', '12'); call beltrami(amplaa,f,iaa,KY=ky_aa)
      case('Beltrami-z', '1');  call beltrami(amplaa,f,iaa,KZ=kz_aa)
      case('propto-ux'); call wave_uu(amplaa,f,iaa,kx=kx_aa)
      case('propto-uy'); call wave_uu(amplaa,f,iaa,ky=ky_aa)
      case('propto-uz'); call wave_uu(amplaa,f,iaa,kz=kz_aa)
      case('diffrot'); call diffrot(amplaa,f,iay,xx,yy,zz)
      case('hor-tube'); call htube(amplaa,f,iax,iaz,xx,yy,zz,radius,epsilonaa)
      case('hor-fluxlayer'); call hfluxlayer(amplaa,f,iaa,xx,yy,zz,z0aa,widthaa)
      case('mag-support'); call magsupport(amplaa,f,zz,gravz,cs0,rho0)
      case('halfcos-Bx'); call halfcos_x(amplaa,f,iaa,xx,yy,zz)
      case('uniform-Bx'); call uniform_x(amplaa,f,iaa,xx,yy,zz)
      case('uniform-By'); call uniform_y(amplaa,f,iaa,xx,yy,zz)
      case('uniform-Bz'); call uniform_z(amplaa,f,iaa,xx,yy,zz)
      case('Bz(x)', '3'); call vfield(amplaa,f,iaa,xx)
      case('xjump'); call bjump(f,iaz,by_left,by_right,widthaa,'x')
      case('fluxrings', '4'); call fluxrings(f,iaa,xx,yy,zz)
      case('sinxsinz'); call sinxsinz(amplaa,f,iaa,kx_aa,ky_aa,kz_aa)
      case('crazy', '5'); call crazy(amplaa,f,iaa)
      case('Alfven-x'); call alfven_x(amplaa,f,iuu,iaa,ilnrho,xx,kx_aa)
      case('Alfven-z'); call alfven_z(amplaa,f,iuu,iaa,zz,kz_aa)
      case('Alfvenz-rot'); call alfvenz_rot(amplaa,f,iuu,iaa,zz,kz_aa,Omega)
      case('Alfven-circ-x')
        !
        !  circularly polarised Alfven wave in x direction
        !
        if (lroot) print*,'init_aa: circular Alfven wave -> x'
        f(:,:,:,iay) = amplaa/kx_aa*sin(kx_aa*xx)
        f(:,:,:,iaz) = amplaa/kx_aa*cos(kx_aa*xx)

      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'init_aa: No such such value for initaa: ', trim(initaa)
        call stop_it("")

      endselect
!
!    If not already used in initaa one can still use kx_aa etc. 
!    to define the wavenumber of the 2nd field. (For old runs!)
!
       if (kx_aa2==impossible) kx_aa2 = kx_aa
       if (ky_aa2==impossible) ky_aa2 = ky_aa
       if (kz_aa2==impossible) kz_aa2 = kz_aa
!
!  superimpose something else
!
      select case(initaa2)
        case('Beltrami-x'); call beltrami(amplaa2,f,iaa,KX=kx_aa2)
        case('Beltrami-y'); call beltrami(amplaa2,f,iaa,KY=ky_aa2)
        case('Beltrami-z'); call beltrami(amplaa2,f,iaa,KZ=kz_aa2)
      endselect
!
!  allow for pressure equilibrium (for isothermal tube)
!  assume that ghost zones have already been set.
!
      if (lpress_equil) then
        if(lroot) print*,'init_aa: adjust lnrho to have pressure equilib; cs0=',cs0
        do n=n1,n2
        do m=m1,m2
          call curl(f,iaa,bb)
          call dot2_mn(bb,b2)
          f(l1:l2,m,n,ilnrho)=f(l1:l2,m,n,ilnrho)-b2/(2.*cs0**2)
        enddo
        enddo
      endif
!
    endsubroutine init_aa
!***********************************************************************
    subroutine daa_dt(f,df,uu,rho1,TT1,uij,bij,bb)
!
!  magnetic field evolution
!
!  calculate dA/dt=uxB+3/2 Omega_0 A_y x_dir -eta mu_0 J +alpha*bb
!  add JxB/rho to momentum equation
!  add eta mu_0 J2/rho to entropy equation
!
!  22-nov-01/nils: coded
!   1-may-02/wolf: adapted for pencil_modular
!  17-jun-03/ulf:  added bx^2, by^2 and bz^2 as separate diagnostics
!   8-aug-03/axel: introduced B_ext21=1./B_ext**2, and set to 1 to prevent division by 0.
!  12-aug-03/christer: added alpha effect (alpha in the equation above)
!
      use Cdata
      use Sub
      use Slices
      use IO, only: output_pencil
      use Special, only: special_calc_magnetic
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: uij,bij
      real, dimension (nx,3) :: bb,aa,jj,uxB,uu,JxB,JxBr,oxuxb,jxbxb
      real, dimension (nx,3) :: gpxb,glnrho,uxj
      real, dimension (nx,3) :: del2A,oo,oxu,bbb,uxDxuxb,del6A,fres
      real, dimension (nx) :: rho1,J2,TT1,b2,b2tot,ab,jb,ub,bx,by,bz,va2
      real, dimension (nx) :: uxb_dotB0,oxuxb_dotB0,jxbxb_dotB0,uxDxuxb_dotB0
      real, dimension (nx) :: gpxb_dotB0,uxj_dotB0,ujxb
      real, dimension (nx) :: bx2, by2, bz2  ! bx^2, by^2 and bz^2
      real, dimension (nx) :: bxby, bxbz, bybz
      real :: tmp,eta_out1,B_ext21=1.
      integer :: j
!
      intent(in)  :: f,uu,rho1,TT1,uij
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
!  calculate B-field, and then max and mean (w/o imposed field, if any)
!
      call curl(f,iaa,bb)
      call dot2_mn(bb,b2)
      if (ldiagnos) then
        aa=f(l1:l2,m,n,iax:iaz)
        if (i_b2m/=0)    call sum_mn_name(b2,i_b2m)
        if (i_bm2/=0) call max_mn_name(b2,i_bm2)
        if (i_brms/=0) call sum_mn_name(b2,i_brms,lsqrt=.true.)
        if (i_bmax/=0) call max_mn_name(b2,i_bmax,lsqrt=.true.)
        if (i_aybym2/=0) call sum_mn_name(2*aa(:,2)*bb(:,2),i_aybym2)
        if (i_abm/=0) then
           call dot_mn(aa,bb,ab)
           call sum_mn_name(ab,i_abm)
        endif
        if (i_ubm/=0) then
           call dot_mn(uu,bb,ub)
           call sum_mn_name(ub,i_ubm)
        endif
        if (i_bx2m/=0) then
           bx2 = bb(:,1)*bb(:,1)
           call sum_mn_name(bx2,i_bx2m)
        endif
        if (i_by2m/=0) then
           by2 = bb(:,2)*bb(:,2)
           call sum_mn_name(by2,i_by2m)
        endif
        if (i_bz2m/=0) then
           bz2 = bb(:,3)*bb(:,3)
           call sum_mn_name(bz2,i_bz2m)
        endif
        if (i_bxbym/=0) then
           bxby = bb(:,1)*bb(:,2)
           call sum_mn_name(bxby,i_bxbym)
        endif
        if (i_bxbzm/=0) then
           bxbz = bb(:,1)*bb(:,3)
           call sum_mn_name(bxbz,i_bxbzm)
        endif
        if (i_bybzm/=0) then
           bybz = bb(:,2)*bb(:,3)
           call sum_mn_name(bybz,i_bybzm)
        endif
!
!  this doesn't need to be as frequent (check later)
!
        if (i_bxmz/=0.or.i_bxmxy/=0) bx=bb(:,1)
        if (i_bymz/=0.or.i_bymxy/=0) by=bb(:,2)
        if (i_bzmz/=0.or.i_bzmxy/=0) bz=bb(:,3)
        if (i_bxmz/=0) call xysum_mn_name_z(bx,i_bxmz)
        if (i_bymz/=0) call xysum_mn_name_z(by,i_bymz)
        if (i_bzmz/=0) call xysum_mn_name_z(bz,i_bzmz)
        if (i_bxmxy/=0) call zsum_mn_name_xy(bx,i_bxmxy)
        if (i_bymxy/=0) call zsum_mn_name_xy(by,i_bymxy)
        if (i_bzmxy/=0) call zsum_mn_name_xy(bz,i_bzmxy)
      endif
!
!  phi-averages
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      if (l2davgfirst) then
        if (i_b2mphi/=0) call phisum_mn_name_rz(b2,i_b2mphi)
        bx=bb(:,1)
        by=bb(:,2)
        bz=bb(:,3)
        if (i_brmphi/=0) call phisum_mn_name_rz(bx*pomx+by*pomy,i_brmphi)
        if (i_bpmphi/=0) call phisum_mn_name_rz(bx*phix+by*phiy,i_bpmphi)
        if (i_bzmphi/=0) call phisum_mn_name_rz(bz,i_bzmphi)
        if (i_b2mphi/=0) call phisum_mn_name_rz(b2,i_b2mphi)
        if (i_jbmphi/=0) then
          call dot_mn(jj,bb,jb)
          call phisum_mn_name_rz(jb,i_jbmphi)
        endif
      endif
!
!  write B-slices for output in wvid in run.f90
!  Note: ix is the index with respect to array with ghost zones.
!
        if(lvid.and.lfirst) then
          do j=1,3
            bb_yz(m-m1+1,n-n1+1,j)=bb(ix-l1+1,j)
            if (m.eq.iy)  bb_xz(:,n-n1+1,j)=bb(:,j)
            if (n.eq.iz)  bb_xy(:,m-m1+1,j)=bb(:,j)
            if (n.eq.iz2) bb_xy2(:,m-m1+1,j)=bb(:,j)
          enddo
          b2_yz(m-m1+1,n-n1+1)=b2(ix-l1+1)
          if (m.eq.iy)  b2_xz(:,n-n1+1)=b2
          if (n.eq.iz)  b2_xy(:,m-m1+1)=b2
          if (n.eq.iz2) b2_xy2(:,m-m1+1)=b2
          if(bthresh_per_brms/=0) call calc_bthresh
          call vecout(41,trim(directory)//'/bvec',bb,bthresh,nbvec)
        endif
!
!  possibility to add external field
!  Note; for diagnostics purposes keep copy of original field
!
      if (ldiagnos) bbb=bb
      if (B_ext(1)/=0.) bb(:,1)=bb(:,1)+B_ext(1)
      if (B_ext(2)/=0.) bb(:,2)=bb(:,2)+B_ext(2)
      if (B_ext(3)/=0.) bb(:,3)=bb(:,3)+B_ext(3)
      if (headtt) print*,'daa_dt: B_ext=',B_ext
!
!  calculating the current jj, and simultaneously del2A.
!
!  --old--   call del2v_etc(f,iaa,del2A,curlcurl=jj)
!  The following two routines also calculate del2A and jj, but they
!  also produce the magnetic field gradient matrix, B_{i,j}, which is
!  needed for cosmic ray evolution. If the overhead in calculating
!  Bij becomes noticeable (even though no extra derivatives are calculated)
!  one should make it switchable between this one and del2v_etc.
!
      call bij_etc(f,iaa,bij,del2A)
      call curl_mn(bij,jj)
!
!  calculate JxB/rho (when hydro is on) and J^2 (when entropy is on)
!  add JxB/rho to momentum equation, and eta mu_0 J2/rho to entropy equation
!
      if (lhydro) then
        call cross_mn(jj,bb,JxB)
        call multsv_mn(rho1,JxB,JxBr)
        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+JxBr
        if(lentropy) then
          call dot2_mn(jj,J2)
          df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+eta*J2*rho1*TT1
        endif
      endif
!
!  calculate uxB+eta*del2A and add to dA/dt
!  (Note: the linear shear term is added later)
!
      call cross_mn(uu,bb,uxB)
!
!  calculate restive term
!  (either normal resisitivity or sixth order hyper resistivity)
!
      if (ires .eq. 'eta-const') then
        fres=eta*del2A
      elseif (ires .eq. 'hyper6') then
        call del6v(f,iaa,del6A)
        fres=eta*del6A
      else
        if (lroot) print*,'daa_dt: no such ires:',ires
        call stop_it("")
      endif
!
!  add to dA/dt
!
      df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+uxB+fres
!
!  add alpha effect if alpha_effect /= 0
!      
      if(alpha_effect/=0.) then
         df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+alpha_effect*bb
      endif
!
!  Possibility of adding extra diffusivity in some halo of given geometry:
!  Note that eta_out is total eta in halo (not eta_out+eta)
!
      if(height_eta/=0.) then
        if (headtt) print*,'daa_dt: height_eta,eta_out=',height_eta,eta_out
        tmp=(z(n)/height_eta)**2
        eta_out1=eta_out*(1.-exp(-tmp**5/amax1(1.-tmp,1e-5)))-eta
        df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)-eta_out1*jj
      endif
!
!  possibility of relaxation of A in exterior region
!
      if (tau_aa_exterior/=0.) call calc_tau_aa_exterior(f,df)
!
!  For the timestep calculation, need maximum Alfven speed.
!  and maximum diffusion (for timestep)
!  This must include the imposed field (if there is any)
!  The b2 calculated above for only updated when diagnos=.true.
!
        if (lfirst.and.ldt) then
          call dot2_mn(bb,b2tot)
          va2=b2tot*rho1
          maxadvec2=amax1(maxadvec2,va2)
          maxdiffus=amax1(maxdiffus,eta)
        endif

      if (lspecial) call special_calc_magnetic(f,df,uu,rho1,TT1,uij)
!
!  calculate max and rms current density
!  at the moment (and in future?) calculate max(b^2) and mean(b^2).
!
      if (ldiagnos) then
        !
        !  magnetic field components at one point (=pt)
        !
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (i_bxpt/=0) call save_name(bbb(lpoint-nghost,1),i_bxpt)
          if (i_bypt/=0) call save_name(bbb(lpoint-nghost,2),i_bypt)
          if (i_bzpt/=0) call save_name(bbb(lpoint-nghost,3),i_bzpt)
        endif
        !
        !  v_A = |B|/sqrt(rho); in units where "4pi"=1
        !
        if (i_vArms/=0) call sum_mn_name(va2,i_vArms,lsqrt=.true.)
        if (i_vAmax/=0) call max_mn_name(va2,i_vAmax,lsqrt=.true.)
        !
        ! <J.B>
        !
        if (i_jbm/=0) then
          call dot_mn(jj,bb,jb)
          call sum_mn_name(jb,i_jbm)
        endif
        !
        ! <J^2> and J^2|max
        !
        if (i_jrms/=0.or.i_jmax/=0.or.i_j2m/=0.or.i_jm2/=0.or.i_epsM/=0) then
          call dot2_mn(jj,j2)
          if (i_j2m/=0) call sum_mn_name(j2,i_j2m)
          if (i_jm2/=0) call max_mn_name(j2,i_jm2)
          if (i_jrms/=0) call sum_mn_name(j2,i_jrms,lsqrt=.true.)
          if (i_jmax/=0) call max_mn_name(j2,i_jmax,lsqrt=.true.)
          if (i_epsM/=0) call sum_mn_name(eta*j2,i_epsM)
        endif
        !
        !  calculate surface integral <2ExA>*dS
        !
        if (i_exaym2/=0) call helflux(aa,uxb,jj)
        !
        !  calculate B_ext21
        !
        B_ext21=B_ext(1)**2+B_ext(2)**2+B_ext(3)**2
        if(B_ext21/=0.) then
          B_ext21=1./B_ext21
        else
          B_ext21=1.
        endif
        !
        !  calculate emf for alpha effect (for imposed field)
        !
        if (i_uxbm/=0.or.i_uxbmx/=0.or.i_uxbmy/=0.or.i_uxbmz/=0) then
          call cross_mn(uu,bbb,uxb)
          uxb_dotB0=B_ext(1)*uxb(:,1)+B_ext(2)*uxb(:,2)+B_ext(3)*uxb(:,3)
          uxb_dotB0=uxb_dotB0*B_ext21
          call sum_mn_name(uxb_dotB0,i_uxbm)
          if (i_uxbmx/=0) call sum_mn_name(uxb(:,1),i_uxbmx)
          if (i_uxbmy/=0) call sum_mn_name(uxb(:,2),i_uxbmy)
          if (i_uxbmz/=0) call sum_mn_name(uxb(:,3),i_uxbmz)
        endif
        !
        !  calculate <uxj>.B0/B0^2
        !
        if (i_uxjm/=0) then
          call cross_mn(jj,uu,uxj)
          uxj_dotB0=B_ext(1)*uxj(:,1)+B_ext(2)*uxj(:,2)+B_ext(3)*uxj(:,3)
          uxj_dotB0=uxj_dotB0*B_ext21
          call sum_mn_name(uxj_dotB0,i_uxjm)
        endif
        !
        !  calculate <u.(jxb)>
        !
        if (i_ujxbm/=0) then
          call cross_mn(jj,bbb,jxb)
          call dot_mn(uu,jxb,ujxb)
          call sum_mn_name(ujxb,i_ujxbm)
        endif
        !
        !  magnetic triple correlation term (for imposed field)
        !
        if (i_jxbxbm/=0) then
          call cross_mn(jj,bbb,jxb)
          call cross_mn(jxb,bbb,jxbxb)
          jxbxb_dotB0=B_ext(1)*jxbxb(:,1)+B_ext(2)*jxbxb(:,2)+B_ext(3)*jxbxb(:,3)
          jxbxb_dotB0=jxbxb_dotB0*B_ext21
          call sum_mn_name(jxbxb_dotB0,i_jxbxbm)
        endif
        !
        !  triple correlation from Reynolds tensor (for imposed field)
        !
        if (i_oxuxbm/=0) then
          oo(:,1)=uij(:,3,2)-uij(:,2,3)
          oo(:,2)=uij(:,1,3)-uij(:,3,1)
          oo(:,3)=uij(:,2,1)-uij(:,1,2)
          call cross_mn(oo,uu,oxu)
          call cross_mn(oxu,bbb,oxuxb)
          oxuxb_dotB0=B_ext(1)*oxuxb(:,1)+B_ext(2)*oxuxb(:,2)+B_ext(3)*oxuxb(:,3)
          oxuxb_dotB0=oxuxb_dotB0*B_ext21
          call sum_mn_name(oxuxb_dotB0,i_oxuxbm)
        endif
        !
        !  triple correlation from pressure gradient (for imposed field)
        !  (assume cs2=1, and that no entropy evolution is included)
        !
        if (i_gpxbm/=0) then
          call grad(f,ilnrho,glnrho)
          call cross_mn(glnrho,bbb,gpxb)
          gpxb_dotB0=B_ext(1)*gpxb(:,1)+B_ext(2)*gpxb(:,2)+B_ext(3)*gpxb(:,3)
          gpxb_dotB0=gpxb_dotB0*B_ext21
          call sum_mn_name(oxuxb_dotB0,i_gpxbm)
        endif
        !
        !  < u x curl(uxB) > = < E_i u_{j,j} - E_j u_{j,i} >
        !   ( < E_1 u2,2 + E1 u3,3 - E2 u2,1 - E3 u3,1 >
        !   ( < E_2 u1,1 + E2 u3,3 - E1 u2,1 - E3 u3,2 >
        !   ( < E_3 u1,1 + E3 u2,2 - E1 u3,1 - E2 u2,3 >
        !
        if (i_uxDxuxbm/=0) then
          call cross_mn(uu,bbb,uxb)
          uxDxuxb(:,1)=uxb(:,1)*(uij(:,2,2)+uij(:,3,3))-uxb(:,2)*uij(:,2,1)-uxb(:,3)*uij(:,3,1)
          uxDxuxb(:,2)=uxb(:,2)*(uij(:,1,1)+uij(:,3,3))-uxb(:,1)*uij(:,1,2)-uxb(:,3)*uij(:,3,2)
          uxDxuxb(:,3)=uxb(:,3)*(uij(:,1,1)+uij(:,2,2))-uxb(:,1)*uij(:,1,3)-uxb(:,2)*uij(:,2,3)
          uxDxuxb_dotB0=B_ext(1)*uxDxuxb(:,1)+B_ext(2)*uxDxuxb(:,2)+B_ext(3)*uxDxuxb(:,3)
          uxDxuxb_dotB0=uxDxuxb_dotB0*B_ext21
          call sum_mn_name(uxDxuxb_dotB0,i_uxDxuxbm)
        endif
        !
      endif
!
!  debug output
!
      if (headtt .and. lfirst .and. ip<=4) then
        call output_pencil(trim(directory)//'/aa.dat',aa,3)
        call output_pencil(trim(directory)//'/bb.dat',bb,3)
        call output_pencil(trim(directory)//'/jj.dat',jj,3)
        call output_pencil(trim(directory)//'/del2A.dat',del2A,3)
        call output_pencil(trim(directory)//'/JxBr.dat',JxBr,3)
        call output_pencil(trim(directory)//'/JxB.dat',JxB,3)
        call output_pencil(trim(directory)//'/df.dat',df(l1:l2,m,n,:),mvar)
      endif
!     
    endsubroutine daa_dt
!***********************************************************************
    subroutine calc_bthresh()
!
!  calculate bthresh from brms, give warnings if there are problems
!
!   6-aug-03/axel: coded
!
      use Cdata
!
!  give warning if brms is not set in prints.in
!
      if(i_brms==0) then
        if(lroot.and.lfirstpoint) then
          print*,'calc_bthresh: need to set brms in print.in to get bthresh'
        endif
      endif
!
!  if nvec exceeds nbvecmax (=1/4) of points per processor, then begin to
!  increase scaling factor on bthresh. These settings will stay in place
!  until the next restart
!
      if(nbvec>nbvecmax.and.lfirstpoint) then
        print*,'calc_bthresh: processor ',iproc,': bthresh_scl,nbvec,nbvecmax=', &
                                                   bthresh_scl,nbvec,nbvecmax
        bthresh_scl=bthresh_scl*1.2
      endif
!
!  calculate bthresh as a certain fraction of brms
!
      bthresh=bthresh_scl*bthresh_per_brms*brms
!
    endsubroutine calc_bthresh
!***********************************************************************
    subroutine calc_tau_aa_exterior(f,df)
!
!  magnetic field relaxation to zero on time scale tau_aa_exterior within
!  exterior region. For the time being this means z > zgrav.
!
!  29-jul-02/axel: coded
!
      use Cdata
      use Gravity
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: scl
      integer :: j
!
      intent(in) :: f
      intent(out) :: df
!
      if (headtt) print*,'calc_tau_aa_exterior: tau=',tau_aa_exterior
      if(z(n)>zgrav) then
        scl=1./tau_aa_exterior
        do j=iax,iaz
          df(l1:l2,m,n,j)=df(l1:l2,m,n,j)-scl*f(l1:l2,m,n,j)
        enddo
      endif
!
    endsubroutine calc_tau_aa_exterior
!***********************************************************************
    subroutine helflux(aa,uxb,jj)
!
!  magnetic helicity flux (preliminary)
!
!  14-aug-03/axel: coded
!
      use Cdata
      use Sub
!
      real, dimension (nx,3), intent(in) :: aa,uxb,jj
      real, dimension (nx,3) :: ee
      real, dimension (nx) :: FHx,FHz
      real :: FH
!
      ee=eta*jj-uxb
!
!  calculate magnetic helicity flux in the X and Z directions
!
      FHx=-2*ee(:,3)*aa(:,2)*dsurfyz
      FHz=+2*ee(:,1)*aa(:,2)*dsurfxy
!
!  sum up contribution per pencil
!  and then stuff result into surf_mn_name for summing up all processors.
!
      FH=FHx(nx)-FHx(1)
      if(ipz==0       .and.n==n1) FH=FH-sum(FHz)
      if(ipz==nprocz-1.and.n==n2) FH=FH+sum(FHz)
      call surf_mn_name(FH,i_exaym2)
!
    endsubroutine helflux
!***********************************************************************
    subroutine rprint_magnetic(lreset,lwrite)
!
!  reads and registers print parameters relevant for magnetic fields
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
        i_b2m=0; i_bm2=0; i_j2m=0; i_jm2=0; i_abm=0; i_jbm=0; i_ubm=0; i_epsM=0
        i_bxpt=0; i_bypt=0; i_bzpt=0
        i_aybym2=0; i_exaym2=0
        i_brms=0; i_bmax=0; i_jrms=0; i_jmax=0; i_vArms=0; i_vAmax=0
        i_bx2m=0; i_by2m=0; i_bz2m=0
        i_bxbym=0; i_bxbzm=0; i_bybzm=0
        i_bxmz=0; i_bymz=0; i_bzmz=0; i_bmx=0; i_bmy=0; i_bmz=0
        i_bxmxy=0; i_bymxy=0; i_bzmxy=0
        i_uxbm=0; i_oxuxbm=0; i_jxbxbm=0.; i_gpxbm=0.; i_uxDxuxbm=0.
        i_uxbmx=0; i_uxbmy=0; i_uxbmz=0
        i_uxjm=0; i_ujxbm=0
        i_brmphi=0; i_bpmphi=0; i_bzmphi=0; i_b2mphi=0; i_jbmphi=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'aybym2',i_aybym2)
        call parse_name(iname,cname(iname),cform(iname),'exaym2',i_exaym2)
        call parse_name(iname,cname(iname),cform(iname),'abm',i_abm)
        call parse_name(iname,cname(iname),cform(iname),'jbm',i_jbm)
        call parse_name(iname,cname(iname),cform(iname),'ubm',i_ubm)
        call parse_name(iname,cname(iname),cform(iname),'b2m',i_b2m)
        call parse_name(iname,cname(iname),cform(iname),'bm2',i_bm2)
        call parse_name(iname,cname(iname),cform(iname),'j2m',i_j2m)
        call parse_name(iname,cname(iname),cform(iname),'jm2',i_jm2)
        call parse_name(iname,cname(iname),cform(iname),'epsM',i_epsM)
        call parse_name(iname,cname(iname),cform(iname),'brms',i_brms)
        call parse_name(iname,cname(iname),cform(iname),'bmax',i_bmax)
        call parse_name(iname,cname(iname),cform(iname),'jrms',i_jrms)
        call parse_name(iname,cname(iname),cform(iname),'jmax',i_jmax)
        call parse_name(iname,cname(iname),cform(iname),'vArms',i_vArms)
        call parse_name(iname,cname(iname),cform(iname),'vAmax',i_vAmax)
        call parse_name(iname,cname(iname),cform(iname),'bx2m',i_bx2m)
        call parse_name(iname,cname(iname),cform(iname),'by2m',i_by2m)
        call parse_name(iname,cname(iname),cform(iname),'bz2m',i_bz2m)
        call parse_name(iname,cname(iname),cform(iname),'bxbym',i_bxbym)
        call parse_name(iname,cname(iname),cform(iname),'bxbzm',i_bxbzm)
        call parse_name(iname,cname(iname),cform(iname),'bybzm',i_bybzm)
        call parse_name(iname,cname(iname),cform(iname),'uxbm',i_uxbm)
        call parse_name(iname,cname(iname),cform(iname),'uxbmx',i_uxbmx)
        call parse_name(iname,cname(iname),cform(iname),'uxbmy',i_uxbmy)
        call parse_name(iname,cname(iname),cform(iname),'uxbmz',i_uxbmz)
        call parse_name(iname,cname(iname),cform(iname),'uxjm',i_uxjm)
        call parse_name(iname,cname(iname),cform(iname),'ujxbm',i_ujxbm)
        call parse_name(iname,cname(iname),cform(iname),'jxbxbm',i_jxbxbm)
        call parse_name(iname,cname(iname),cform(iname),'oxuxbm',i_oxuxbm)
        call parse_name(iname,cname(iname),cform(iname),'gpxbm',i_gpxbm)
        call parse_name(iname,cname(iname),cform(iname),'uxDxuxbm',i_uxDxuxbm)
        call parse_name(iname,cname(iname),cform(iname),'bmx',i_bmx)
        call parse_name(iname,cname(iname),cform(iname),'bmy',i_bmy)
        call parse_name(iname,cname(iname),cform(iname),'bmz',i_bmz)
        call parse_name(iname,cname(iname),cform(iname),'bxpt',i_bxpt)
        call parse_name(iname,cname(iname),cform(iname),'bypt',i_bypt)
        call parse_name(iname,cname(iname),cform(iname),'bzpt',i_bzpt)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bxmz',i_bxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bymz',i_bymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bzmz',i_bzmz)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do ixy=1,nnamexy
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bxmxy',i_bxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bymxy',i_bymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bzmxy',i_bzmxy)
      enddo
!
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'brmphi',i_brmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'bpmphi',i_bpmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'bzmphi',i_bzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'b2mphi',i_b2mphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'jbmphi',i_jbmphi)
      enddo
!
!  write column, i_XYZ, where our variable XYZ is stored
!
      if (lwr) then
        write(3,*) 'i_aybym2=',i_aybym2
        write(3,*) 'i_exaym2=',i_exaym2
        write(3,*) 'i_abm=',i_abm
        write(3,*) 'i_jbm=',i_jbm
        write(3,*) 'i_ubm=',i_ubm
        write(3,*) 'i_b2m=',i_b2m
        write(3,*) 'i_bm2=',i_bm2
        write(3,*) 'i_j2m=',i_j2m
        write(3,*) 'i_jm2=',i_jm2
        write(3,*) 'i_epsM=',i_epsM
        write(3,*) 'i_brms=',i_brms
        write(3,*) 'i_bmax=',i_bmax
        write(3,*) 'i_jrms=',i_jrms
        write(3,*) 'i_jmax=',i_jmax
        write(3,*) 'i_vArms=',i_vArms
        write(3,*) 'i_vAmax=',i_vAmax
        write(3,*) 'i_bx2m=',i_bx2m
        write(3,*) 'i_by2m=',i_by2m
        write(3,*) 'i_bz2m=',i_bz2m
        write(3,*) 'i_bxbym=',i_bxbym
        write(3,*) 'i_bxbzm=',i_bxbzm
        write(3,*) 'i_bybzm=',i_bybzm
        write(3,*) 'i_uxbm=',i_uxbm
        write(3,*) 'i_uxbmx=',i_uxbmx
        write(3,*) 'i_uxbmy=',i_uxbmy
        write(3,*) 'i_uxbmz=',i_uxbmz
        write(3,*) 'i_uxjm=',i_uxjm
        write(3,*) 'i_ujxbm=',i_ujxbm
        write(3,*) 'i_oxuxbm=',i_oxuxbm
        write(3,*) 'i_jxbxbm=',i_jxbxbm
        write(3,*) 'i_gpxbm=',i_gpxbm
        write(3,*) 'i_uxDxuxbm=',i_uxDxuxbm
        write(3,*) 'nname=',nname
        write(3,*) 'iaa=',iaa
        write(3,*) 'iax=',iax
        write(3,*) 'iay=',iay
        write(3,*) 'iaz=',iaz
        write(3,*) 'nnamez=',nnamez
        write(3,*) 'i_bxmz=',i_bxmz
        write(3,*) 'i_bymz=',i_bymz
        write(3,*) 'i_bzmz=',i_bzmz
        write(3,*) 'i_bmx=',i_bmx
        write(3,*) 'i_bmy=',i_bmy
        write(3,*) 'i_bmz=',i_bmz
        write(3,*) 'i_bxpt=',i_bxpt
        write(3,*) 'i_bypt=',i_bypt
        write(3,*) 'i_bzpt=',i_bzpt
        write(3,*) 'nnamexy=',nnamexy
        write(3,*) 'i_bxmxy=',i_bxmxy
        write(3,*) 'i_bymxy=',i_bymxy
        write(3,*) 'i_bzmxy=',i_bzmxy
        write(3,*) 'i_brmphi=',i_brmphi
        write(3,*) 'i_bpmphi=',i_bpmphi
        write(3,*) 'i_bzmphi=',i_bzmphi
        write(3,*) 'i_b2mphi=',i_b2mphi
        write(3,*) 'i_jbmphi=',i_jbmphi
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
      use Cdata
      use Sub
!
      logical,save :: first=.true.
      real, dimension(nx) :: bymx,bzmx
      real, dimension(ny,nprocy) :: bxmy,bzmy
      real :: bmx,bmy,bmz
      integer :: l,j
!
!  Magnetic energy in vertically averaged field
!  The bymxy and bzmxy must have been calculated,
!  so they are present on the root processor.
!
        if (i_bmx/=0) then
          if(i_bymxy==0.or.i_bzmxy==0) then
            if(first) print*,"calc_mfield:                  WARNING"
            if(first) print*, &
                    "calc_mfield: NOTE: to get bmx, bymxy and bzmxy must also be set in zaver"
            if(first) print*, &
                    "calc_mfield:       We proceed, but you'll get bmx=0"
            bmx=0.
          else
            do l=1,nx
              bymx(l)=sum(fnamexy(l,:,:,i_bymxy))/(ny*nprocy)
              bzmx(l)=sum(fnamexy(l,:,:,i_bzmxy))/(ny*nprocy)
            enddo
            bmx=sqrt(sum(bymx**2+bzmx**2)/nx)
          endif
          call save_name(bmx,i_bmx)
        endif
!
!  similarly for bmy
!
        if (i_bmy/=0) then
          if(i_bxmxy==0.or.i_bzmxy==0) then
            if(first) print*,"calc_mfield:                  WARNING"
            if(first) print*, &
                    "calc_mfield: NOTE: to get bmy, bxmxy and bzmxy must also be set in zaver"
            if(first) print*, &
                    "calc_mfield:       We proceed, but you'll get bmy=0"
            bmy=0.
          else
            do j=1,nprocy
            do m=1,ny
              bxmy(m,j)=sum(fnamexy(:,m,j,i_bxmxy))/nx
              bzmy(m,j)=sum(fnamexy(:,m,j,i_bzmxy))/nx
            enddo
            enddo
            bmy=sqrt(sum(bxmy**2+bzmy**2)/(ny*nprocy))
          endif
          call save_name(bmy,i_bmy)
        endif
!
!  Magnetic energy in horizontally averaged field
!  The bxmz and bymz must have been calculated,
!  so they are present on the root processor.
!
        if (i_bmz/=0) then
          if(i_bxmz==0.or.i_bymz==0) then
            if(first) print*,"calc_mfield:                  WARNING"
            if(first) print*, &
                    "calc_mfield: NOTE: to get bmz, bxmz and bymz must also be set in xyaver"
            if(first) print*, &
                    "calc_mfield:       This may be because we renamed zaver.in into xyaver.in"
            if(first) print*, &
                    "calc_mfield:       We proceed, but you'll get bmz=0"
            bmz=0.
          else
            bmz=sqrt(sum(fnamez(:,:,i_bxmz)**2+fnamez(:,:,i_bymz)**2)/(nz*nprocz))
          endif
          call save_name(bmz,i_bmz)
        endif
!
      first = .false.
    endsubroutine calc_mfield
!***********************************************************************
    subroutine alfven_x(ampl,f,iuu,iaa,ilnrho,xx,kx)
!
!  Alfven wave propagating in the z-direction
!  ux = cos(kz-ot), for B0z=1 and rho=1.
!  Ay = sin(kz-ot), ie Bx=-cos(kz-ot)
!
!  satisfies the equations
!  dlnrho/dt = -ux'
!  dux/dt = -cs2*(lnrho)'
!  duy/dt = B0*By'  ==>  dux/dt = B0*Ay''
!  dBy/dt = B0*uy'  ==>  dAy/dt = -B0*ux
!
!   8-nov-03/axel: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx
      real :: ampl,kx
      integer :: iuu,iaa,ilnrho
!
!  ux and Ay
!
      f(:,:,:,ilnrho)=ampl*sin(kx*xx)
      f(:,:,:,iuu+0)=+ampl*sin(kx*xx)
      f(:,:,:,iuu+1)=+ampl*sin(kx*xx)
      f(:,:,:,iaa+2)=-ampl*cos(kx*xx)
!
    endsubroutine alfven_x
!***********************************************************************
    subroutine alfven_z(ampl,f,iuu,iaa,zz,kz)
!
!  Alfven wave propagating in the z-direction
!  ux = cos(kz-ot), for B0z=1 and rho=1.
!  Ay = sin(kz-ot), ie Bx=-cos(kz-ot)
!
!  satisfies the equations
!  dux/dt = Bx'  ==>  dux/dt = -Ay''
!  dBx/dt = ux'  ==>  dAy/dt = -ux.
!
!  18-aug-02/axel: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: zz
      real :: ampl,kz
      integer :: iuu,iaa
!
!  ux and Ay
!
      f(:,:,:,iuu+0)=+ampl*cos(kz*zz)
      f(:,:,:,iaa+1)=+ampl*sin(kz*zz)
!
    endsubroutine alfven_z
!***********************************************************************
    subroutine alfvenz_rot(ampl,f,iuu,iaa,zz,kz,O)
!
!  Alfven wave propagating in the z-direction
!  ux = cos(kz-ot), for B0z=1 and rho=1.
!  Ay = sin(kz-ot), ie Bx=-cos(kz-ot)
!
!  satisfies the equations
!  dux/dt - 2Omega*uy = -Ay''
!  duy/dt + 2Omega*ux = +Ax''
!  dAx/dt = +uy
!  dAy/dt = -ux
!
!  18-aug-02/axel: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: zz
      real :: ampl,kz,O,fac
      integer :: iuu,iaa
!
!  ux and Ay
!
      print*,'alfvenz_rot: Alfven wave with rotation; O,kz=',O,kz
      fac=-O+sqrt(O**2+kz**2)
      f(:,:,:,iuu+0)=+ampl*cos(kz*zz)*fac/kz
      f(:,:,:,iuu+1)=-ampl*sin(kz*zz)*fac/kz
      f(:,:,:,iaa+0)=-ampl*cos(kz*zz)
      f(:,:,:,iaa+1)=+ampl*sin(kz*zz)
!
    endsubroutine alfvenz_rot
!***********************************************************************
    subroutine fluxrings(f,ivar,xx,yy,zz,profile)
!
!  Magnetic flux rings. Constructed from a canonical ring which is the
!  rotated and translated:
!    AA(xxx) = D*AA0(D^(-1)*(xxx-xxx_disp)) ,
!  where AA0(xxx) is the canonical ring and D the rotation matrix
!  corresponding to a rotation by phi around z, followed by a
!  rotation by theta around y.
!  The array was already initialized to zero before calling this
!  routine.
!  Optional argument `profile' allows to choose a different profile (see
!  norm_ring())
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,3)    :: tmpv
      real, dimension (mx,my,mz)      :: xx,yy,zz,xx1,yy1,zz1
      real, dimension(3) :: axis,disp
      real    :: phi,theta,ct,st,cp,sp
      real    :: fring,Iring,R0,width
      integer :: i,ivar
      character (len=*), optional :: profile
      character (len=labellen) :: prof
!
      if (present(profile)) then
        prof = profile
      else
        prof = 'tanh'
      endif

      if (any((/fring1,fring2,Iring1,Iring2/) /= 0.)) then
        ! fringX is the magnetic flux, IringX the current
        if (lroot) then
          print*, 'fluxrings: Initialising magnetic flux rings'
        endif
        do i=1,2
          if (i==1) then
            fring = fring1      ! magnetic flux along ring
            Iring = Iring1      ! current along ring (for twisted flux tube)
            R0    = Rring1      ! radius of ring
            width = wr1         ! ring thickness
            axis  = axisr1 ! orientation
            disp  = dispr1    ! position
          else
            fring = fring2
            Iring = Iring2
            R0    = Rring2
            width = wr2
            axis  = axisr2
            disp  = dispr2
          endif
          phi   = atan2(axis(2),axis(1)+epsi)
          theta = atan2(sqrt(axis(1)**2+axis(2)**2)+epsi,axis(3))
          ct = cos(theta); st = sin(theta)
          cp = cos(phi)  ; sp = sin(phi)
          ! Calculate D^(-1)*(xxx-disp)
          xx1 =  ct*cp*(xx-disp(1)) + ct*sp*(yy-disp(2)) - st*(zz-disp(3))
          yy1 = -   sp*(xx-disp(1)) +    cp*(yy-disp(2))
          zz1 =  st*cp*(xx-disp(1)) + st*sp*(yy-disp(2)) + ct*(zz-disp(3))
          call norm_ring(xx1,yy1,zz1,fring,Iring,R0,width,tmpv,PROFILE=prof)
          ! calculate D*tmpv
          f(:,:,:,ivar  ) = f(:,:,:,ivar  ) + amplaa*( &
               + ct*cp*tmpv(:,:,:,1) - sp*tmpv(:,:,:,2) + st*cp*tmpv(:,:,:,3))
          f(:,:,:,ivar+1) = f(:,:,:,ivar+1) + amplaa*( &
               + ct*sp*tmpv(:,:,:,1) + cp*tmpv(:,:,:,2) + st*sp*tmpv(:,:,:,3))
          f(:,:,:,ivar+2) = f(:,:,:,ivar+2) + amplaa*( &
               - st   *tmpv(:,:,:,1)                    + ct   *tmpv(:,:,:,3))
        enddo
      endif
      if (lroot) print*, 'fluxrings: Magnetic flux rings initialized'
!
    endsubroutine fluxrings
!***********************************************************************
    subroutine norm_ring(xx,yy,zz,fring,Iring,R0,width,vv,profile)
!
!  Generate vector potential for a flux ring of magnetic flux FRING,
!  current Iring (not correctly normalized), radius R0 and thickness
!  WIDTH in normal orientation (lying in the x-y plane, centred at (0,0,0)).
!
!   1-may-02/wolf: coded
!
      use Cdata, only: mx,my,mz
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,3) :: vv
      real, dimension (mx,my,mz)   :: xx,yy,zz,phi,tmp
      real :: fring,Iring,R0,width
      character (len=*) :: profile
!
      vv = 0.
!
!  magnetic ring
!
      tmp = sqrt(xx**2+yy**2)-R0

      select case(profile)

      case('tanh')
        vv(:,:,:,3) = - fring * 0.5*(1+tanh(tmp/width)) &
                              * 0.5/width/cosh(zz/width)**2

      case default
        call stop_it('norm_ring: No such fluxtube profile')
      endselect
!
!  current ring (to twist the B-lines)
!
!      tmp = tmp**2 + zz**2 + width**2  ! need periodic analog of this
      tmp = width - sqrt(tmp**2 + zz**2)
      tmp = Iring*0.5*(1+tanh(tmp/width))     ! Now the A_phi component
      phi = atan2(yy,xx)
      vv(:,:,:,1) = - tmp*sin(phi)
      vv(:,:,:,2) =   tmp*cos(phi)
!
    endsubroutine norm_ring
!***********************************************************************
      subroutine bc_aa_pot(f,topbot)
!
!  Potential field boundary condition for magnetic vector potential
!
!  14-jun-2002/axel: adapted from similar 
!   8-jul-2002/axel: introduced topbot argument
!
      use Cdata
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx,ny) :: f2,f3
      real, dimension (nx,ny,nghost+1) :: fz
      integer :: j
!
!  pontential field condition
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  pontential field condition at the bottom
!
      case('bot')
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
!  pontential field condition at the top
!
      case('top')
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
        if(lroot) print*,"bc_aa_pot: invalid argument"
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
     use Cdata
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
      call fft(f2r, f2i, nx*ny, nx,    nx,-1) ! x-direction
      call fft(f2r, f2i, nx*ny, ny, nx*ny,-1) ! y-direction
!
      call fft(f3r, f3i, nx*ny, nx,    nx,-1) ! x-direction
      call fft(f3r, f3i, nx*ny, ny, nx*ny,-1) ! y-direction
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
        call fft(g1r, g1i, nx*ny, nx,    nx,+1) ! x-direction
        call fft(g1r, g1i, nx*ny, ny, nx*ny,+1) ! y-direction
!
!  reverse order if irev=-1 (if we are at the bottom)
!
        if(irev==+1) fz(:,:,       i+1) = g1r/(nx*ny)  ! Renormalize
        if(irev==-1) fz(:,:,nghost-i+1) = g1r/(nx*ny)  ! Renormalize
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
     use Cdata
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
      call fft(f2r, f2i, nx*ny, nx,    nx,-1) ! x-direction
      call fft(f2r, f2i, nx*ny, ny, nx*ny,-1) ! y-direction
!
      call fft(f3r, f3i, nx*ny, nx,    nx,-1) ! x-direction
      call fft(f3r, f3i, nx*ny, ny, nx*ny,-1) ! y-direction
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
        call fft(g1r, g1i, nx*ny, nx,    nx,+1) ! x-direction
        call fft(g1r, g1i, nx*ny, ny, nx*ny,+1) ! y-direction
!
!  reverse order if irev=-1 (if we are at the bottom)
!
        if(irev==+1) fz(:,:,       i+1) = g1r/(nx*ny)  ! Renormalize
        if(irev==-1) fz(:,:,nghost-i+1) = g1r/(nx*ny)  ! Renormalize
      enddo
!
    endsubroutine potentdiv
!***********************************************************************

endmodule Magnetic
