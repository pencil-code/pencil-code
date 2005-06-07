! $Id: magnetic.f90,v 1.236 2005-06-07 21:21:28 brandenb Exp $

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
  character (len=labellen) :: iresistivity='eta-const'
  ! input parameters
  real, dimension(3) :: B_ext=(/0.,0.,0./),B_ext_tmp
  real, dimension(3) :: axisr1=(/0,0,1/),dispr1=(/0.,0.5,0./)
  real, dimension(3) :: axisr2=(/1,0,0/),dispr2=(/0.,-0.5,0./)
  real, dimension (nx,3) :: bbb
  real :: fring1=0.,Iring1=0.,Rring1=1.,wr1=0.3
  real :: fring2=0.,Iring2=0.,Rring2=1.,wr2=0.3
  real :: amplaa=0., kx_aa=1.,ky_aa=1.,kz_aa=1.
  real :: radius=.1,epsilonaa=1e-2,widthaa=.5,x0aa=0.,z0aa=0.
  real :: by_left=0.,by_right=0.,bz_left=0.,bz_right=0.
  real :: ABC_A=1.,ABC_B=1.,ABC_C=1.
  real :: amplaa2=0.,kx_aa2=impossible,ky_aa2=impossible,kz_aa2=impossible
  real :: bthresh=0.,bthresh_per_brms=0.,brms=0.,bthresh_scl=1.
  real :: eta_shock=0.
  real :: rhomin_JxB=0.,va2max_JxB=0.
  real :: omega_Bz_ext
  real :: mu_r=-0.5 !(still needed for backwards compatibility)
  real :: mu_ext_pot=-0.5
  real :: rescale_aa=1.
  real :: ampl_B0=0.,D_smag=0.17
  integer :: nbvec,nbvecmax=nx*ny*nz/4,va2power_JxB=5
  logical :: lpress_equil=.false., lpress_equil_via_ss=.false.
  logical :: llorentzforce=.true.,linduction=.true.
  ! dgm: for hyper diffusion in any spatial variation of eta
  logical :: lresistivity_hyper=.false.,leta_const=.true.
  logical :: lfrozen_bz_z_bot=.false.,lfrozen_bz_z_top=.false.
  logical :: reinitalize_aa=.false.
  logical :: lB_ext_pot=.false.
  logical :: lee_ext=.false.,lbb_ext=.false.,ljj_ext=.false.
  logical :: lforce_free_test=.false.
  character (len=40) :: kinflow=''
  real :: nu_ni,nu_ni1,hall_term,alpha_effect
  real :: displacement_gun=0.
  complex, dimension(3) :: coefaa=(/0.,0.,0./), coefbb=(/0.,0.,0./)
  ! dgm: for perturbing magnetic field when reading NON-magnetic snapshot
  real :: pertamplaa=0.
  real :: initpower_aa=0.,cutoff_aa=0.
  character (len=labellen) :: pertaa='zero'
  integer :: N_modes_aa=1

  namelist /magnetic_init_pars/ &
       B_ext, &
       fring1,Iring1,Rring1,wr1,axisr1,dispr1, &
       fring2,Iring2,Rring2,wr2,axisr2,dispr2, &
       radius,epsilonaa,x0aa,z0aa,widthaa, &
       by_left,by_right,bz_left,bz_right, &
       initaa,initaa2,amplaa,amplaa2,kx_aa,ky_aa,kz_aa,coefaa,coefbb, &
       kx_aa2,ky_aa2,kz_aa2,lpress_equil,lpress_equil_via_ss,mu_r, &
       mu_ext_pot,lB_ext_pot,lforce_free_test, &
       ampl_B0,initpower_aa,cutoff_aa,N_modes_aa

  ! run parameters
  real :: eta=0.,height_eta=0.,eta_out=0.
  real :: eta_int=0.,eta_ext=0.,wresistivity=.01
  real :: tau_aa_exterior=0.

  namelist /magnetic_run_pars/ &
       eta,B_ext,omega_Bz_ext,alpha_effect,nu_ni,hall_term, &
       height_eta,eta_out,tau_aa_exterior, &
       kinflow,kx_aa,ky_aa,kz_aa,ABC_A,ABC_B,ABC_C, &
       bthresh,bthresh_per_brms, &
       iresistivity,lresistivity_hyper, &
       eta_int,eta_ext,eta_shock,wresistivity, &
       rhomin_JxB,va2max_JxB,va2power_JxB,llorentzforce,linduction, &
       reinitalize_aa,rescale_aa,lB_ext_pot, &
       lee_ext,lbb_ext,ljj_ext,displacement_gun, &
       pertaa,pertamplaa,D_smag

  ! other variables (needs to be consistent with reset list below)
  integer :: i_b2m=0,i_bm2=0,i_j2m=0,i_jm2=0,i_abm=0,i_jbm=0,i_ubm,i_epsM=0
  integer :: i_bxpt=0,i_bypt=0,i_bzpt=0,i_epsM_LES=0
  integer :: i_aybym2=0,i_exaym2=0,i_exjm2=0
  integer :: i_brms=0,i_bmax=0,i_jrms=0,i_jmax=0,i_vArms=0,i_vAmax=0,i_dtb=0
  integer :: i_beta1m=0,i_beta1max=0
  integer :: i_bx2m=0, i_by2m=0, i_bz2m=0
  integer :: i_bxbym=0, i_bxbzm=0, i_bybzm=0,i_djuidjbim
  integer :: i_bxmz=0,i_bymz=0,i_bzmz=0,i_bmx=0,i_bmy=0,i_bmz=0
  integer :: i_bxmxy=0,i_bymxy=0,i_bzmxy=0
  integer :: i_bxmxz=0,i_bymxz=0,i_bzmxz=0
  integer :: i_uxbm=0,i_oxuxbm=0,i_jxbxbm=0,i_gpxbm=0,i_uxDxuxbm=0
  integer :: i_uxbmx=0,i_uxbmy=0,i_uxbmz=0,i_uxjm=0,i_ujxbm
  integer :: i_b2b13m=0
  integer :: i_brmphi=0,i_bpmphi=0,i_bzmphi=0,i_b2mphi=0,i_jbmphi=0
  integer :: i_uxbrmphi,i_uxbpmphi,i_uxbzmphi
  integer :: i_dteta=0

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
!  Put variable names in array
!
      varname(iax) = 'ax'
      varname(iay) = 'ay'
      varname(iaz) = 'az'
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: magnetic.f90,v 1.236 2005-06-07 21:21:28 brandenb Exp $")
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
    subroutine initialize_magnetic(f)
!
!  Perform any post-parameter-read initialization
!
!  24-nov-02/tony: dummy routine - nothing to do at present
!  20-may-03/axel: reinitalize_aa added
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
!  set to zero and then rescale the magnetic field
!  (in future, could call something like init_aa_simple)
!
      if (reinitalize_aa) then
        f(:,:,:,iax:iaz)=rescale_aa*f(:,:,:,iax:iaz)
      endif
!
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
      real, dimension (mx,my,mz)      :: xx,yy,zz,tmp,prof
      real, dimension (nx,3) :: bb
      real, dimension (nx) :: b2,fact
      real :: beq2
!
      select case(initaa)

      case('zero', '0'); f(:,:,:,iax:iaz) = 0.
      case('rescale'); f(:,:,:,iax:iaz)=amplaa*f(:,:,:,iax:iaz)
      case('mode'); call modev(amplaa,coefaa,f,iaa,kx_aa,ky_aa,kz_aa,xx,yy,zz)
      case('modeb'); call modeb(amplaa,coefbb,f,iaa,kx_aa,ky_aa,kz_aa,xx,yy,zz)
      case('power_randomphase')
         call power_randomphase(amplaa,initpower_aa,cutoff_aa,f,iax,iaz)
      case('random-isotropic-KS')
         call random_isotropic_KS(amplaa,initpower_aa,cutoff_aa,f,iax,iaz,N_modes_aa)
      case('gaussian-noise'); call gaunoise(amplaa,f,iax,iaz)
      case('gaussian-noise-rprof')
        tmp=sqrt(xx**2+yy**2+zz**2)
        call gaunoise_rprof(amplaa,tmp,prof,f,iax,iaz)
      case('Beltrami-x', '11'); call beltrami(amplaa,f,iaa,KX=kx_aa)
      case('Beltrami-y', '12'); call beltrami(amplaa,f,iaa,KY=ky_aa)
      case('Beltrami-z', '1');  call beltrami(amplaa,f,iaa,KZ=kz_aa)
      case('propto-ux'); call wave_uu(amplaa,f,iaa,kx=kx_aa)
      case('propto-uy'); call wave_uu(amplaa,f,iaa,ky=ky_aa)
      case('propto-uz'); call wave_uu(amplaa,f,iaa,kz=kz_aa)
      case('diffrot'); call diffrot(amplaa,f,iay,xx,yy,zz)
      case('hor-tube'); call htube(amplaa,f,iax,iaz,xx,yy,zz,radius,epsilonaa)
      case('hor-fluxlayer'); call hfluxlayer(amplaa,f,iaa,xx,yy,zz,z0aa,widthaa)
      case('ver-fluxlayer'); call vfluxlayer(amplaa,f,iaa,xx,yy,zz,x0aa,widthaa)
      case('mag-support'); call magsupport(amplaa,f,zz,gravz,cs0,rho0)
      case('arcade-x'); call arcade_x(amplaa,f,iaa,xx,yy,zz,kx_aa,kz_aa)
      case('halfcos-Bx'); call halfcos_x(amplaa,f,iaa,xx,yy,zz)
      case('uniform-Bx'); call uniform_x(amplaa,f,iaa,xx,yy,zz)
      case('uniform-By'); call uniform_y(amplaa,f,iaa,xx,yy,zz)
      case('uniform-Bz'); call uniform_z(amplaa,f,iaa,xx,yy,zz)
      case('Bz(x)', '3'); call vfield(amplaa,f,iaa,xx)
      case('vfield2'); call vfield2(amplaa,f,iaa,xx)
      case('xjump'); call bjump(f,iaa,by_left,by_right,bz_left,bz_right,widthaa,'x')
      case('fluxrings', '4'); call fluxrings(f,iaa,xx,yy,zz)
      case('sinxsinz'); call sinxsinz(amplaa,f,iaa,kx_aa,ky_aa,kz_aa)
      case('cosxcosy'); call cosx_cosy_cosz(amplaa,f,iaz,kx_aa,ky_aa,0.)
      case('sinxsiny'); call sinx_siny_cosz(amplaa,f,iaz,kx_aa,ky_aa,0.)
      case('cosxcoscosy'); call cosx_coscosy_cosz(amplaa,f,iaz,kx_aa,ky_aa,0.)
      case('crazy', '5'); call crazy(amplaa,f,iaa)
      case('Alfven-x'); call alfven_x(amplaa,f,iuu,iaa,ilnrho,xx,kx_aa)
      case('Alfven-z'); call alfven_z(amplaa,f,iuu,iaa,zz,kz_aa,mu0)
      case('Alfvenz-rot'); call alfvenz_rot(amplaa,f,iuu,iaa,zz,kz_aa,Omega)
      case('Alfvenz-rot-shear'); call alfvenz_rot_shear(amplaa,f,iuu,iaa,zz,kz_aa,Omega)
      case('force-free-jet')
        lB_ext_pot=.true.
        call force_free_jet(mu_ext_pot,xx,yy,zz)
      case('Alfven-circ-x')
        !
        !  circularly polarised Alfven wave in x direction
        !
        if (lroot) print*,'init_aa: circular Alfven wave -> x'
        f(:,:,:,iay) = amplaa/kx_aa*sin(kx_aa*xx)
        f(:,:,:,iaz) = amplaa/kx_aa*cos(kx_aa*xx)
      case('geo-benchmark-case1','geo-benchmark-case2'); call geo_benchmark_B(f)

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
        case('gaussian-noise'); call gaunoise(amplaa2,f,iax,iaz)
      endselect
!
!  allow for pressure equilibrium (for isothermal tube)
!  assume that ghost zones have already been set.
!  corrected expression below for gamma /= 1 case.
!  The beq2 expression for 2*mu0*p is not general yet.
!
      if (lpress_equil.or.lpress_equil_via_ss) then
        if(lroot) print*,'init_aa: adjust lnrho to have pressure equilib; cs0=',cs0
        do n=n1,n2
        do m=m1,m2
          call curl(f,iaa,bb)
          call dot2_mn(bb,b2)
          if (gamma==1.) then
            f(l1:l2,m,n,ilnrho)=f(l1:l2,m,n,ilnrho)-b2/(2.*cs0**2)
          else
            beq2=2.*rho0*cs0**2
            fact=max(1e-6,1.-b2/beq2)
            if (lentropy.and.lpress_equil_via_ss) then
              f(l1:l2,m,n,iss)=f(l1:l2,m,n,iss)+fact/gamma
            else
              f(l1:l2,m,n,ilnrho)=f(l1:l2,m,n,ilnrho)+fact/gamma1
            endif
          endif
        enddo
        enddo
      endif
!
    endsubroutine init_aa
!***********************************************************************
    subroutine pert_aa(f)
!
!   perturb magnetic field when reading old NON-magnetic snapshot
!   called from run.f90
!   30-july-2004/dave: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      xx=spread(spread(x,2,my),3,mz)
      yy=spread(spread(y,1,mx),3,mz)
      zz=spread(spread(z,1,mx),2,my)
      initaa=pertaa
      amplaa=pertamplaa
      call init_aa(f,xx,yy,zz)
!
    endsubroutine pert_aa
!***********************************************************************
    subroutine daa_dt(f,df,uu,uij,rho1,TT1,bb,bij,aij,jj,JxBr,del2A,graddivA,va2,shock,gshock)
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
!   8-aug-03/axel: introduced B_ext21=1./B_ext**2, and set =1 to avoid div. by 0
!  12-aug-03/christer: added alpha effect (alpha in the equation above)
!  26-may-04/axel: ambipolar diffusion added
!  18-jun-04/axel: Hall term added
!
      use Cdata
      use Sub
      use Slices
      use Global, only: get_global
      use IO, only: output_pencil
      use Special, only: special_calc_magnetic
      use Mpicomm, only: stop_it
      use Ionization, only: eoscalc,gamma1
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: uij,bij,aij
      real, dimension (nx,3) :: bb,aa,jj,uxB,uu,JxB,JxBr,oxuxb,jxbxb,JxBrxB
      real, dimension (nx,3) :: gpxb,glnrho,uxj,gshock,geta
      real, dimension (nx,3) :: del2A,graddivA,oo,oxu,uxDxuxb,del6A,fres,del4A
      real, dimension (nx,3) :: ee_ext
      real, dimension (nx) :: rho1,J2,TT1,b2,b2tot,ab,jb,ub,bx,by,bz,va2
      real, dimension (nx) :: uxb_dotB0,oxuxb_dotB0,jxbxb_dotB0,uxDxuxb_dotB0
      real, dimension (nx) :: gpxb_dotB0,uxj_dotB0,ujxb,shock,rho1_JxB
      real, dimension (nx) :: hall_ueff2
      real, dimension (nx) :: bx2, by2, bz2  ! bx^2, by^2 and bz^2
      real, dimension (nx) :: bxby, bxbz, bybz
      real, dimension (nx) :: b2b13,jo,sign_jo
      real, dimension (nx) :: eta_mn,divA,eta_tot,del4A2        ! dgm: 
      real, dimension (nx) :: etatotal,pp,djuidjbi
      real :: tmp,eta_out1,B_ext21=1.
      integer :: j,i
      real, dimension (nx) :: eta_smag
!
      intent(in)     :: f,uu,rho1,TT1,uij,bij,aij,bb,jj,del2A,shock,gshock
      intent(out)    :: va2
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
!  Note that the B field is now calculated in calculate_vars_magnetic
!
      call dot2_mn(bb,b2)
!
!  Calculate some diagnostic quantities
!
      if (ldiagnos) then

        if (gamma1/=0.) then
          call eoscalc(f,nx,pp=pp)
          if (i_beta1m/=0) call sum_mn_name(0.5*b2/pp,i_beta1m)
          if (i_beta1max/=0) call max_mn_name(0.5*b2/pp,i_beta1max)
        endif

        aa=f(l1:l2,m,n,iax:iaz)

        if (i_b2m/=0) call sum_mn_name(b2,i_b2m)
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

        if (i_djuidjbim/=0) then
           call multmm_sc(uij,bij,djuidjbi)
           call sum_mn_name(djuidjbi,i_djuidjbim)
        endif 

!
!  this doesn't need to be as frequent (check later)
!
        if (i_bxmz/=0.or.i_bxmxy/=0.or.i_bxmxz) bx=bb(:,1)
        if (i_bymz/=0.or.i_bymxy/=0.or.i_bymxz) by=bb(:,2)
        if (i_bzmz/=0.or.i_bzmxy/=0.or.i_bzmxz) bz=bb(:,3)
        if (i_bxmz/=0) call xysum_mn_name_z(bx,i_bxmz)
        if (i_bymz/=0) call xysum_mn_name_z(by,i_bymz)
        if (i_bzmz/=0) call xysum_mn_name_z(bz,i_bzmz)
        if (i_bxmxy/=0) call zsum_mn_name_xy(bx,i_bxmxy)
        if (i_bymxy/=0) call zsum_mn_name_xy(by,i_bymxy)
        if (i_bzmxy/=0) call zsum_mn_name_xy(bz,i_bzmxy)
        if (i_bxmxz/=0) call ysum_mn_name_xz(bx,i_bxmxz)
        if (i_bymxz/=0) call ysum_mn_name_xz(by,i_bymxz)
        if (i_bzmxz/=0) call ysum_mn_name_xz(bz,i_bzmxz)
      endif
!
!  calculate Alfven speed
!  This must include the imposed field (if there is any)
!  The b2 calculated above for only updated when diagnos=.true.
!
      call dot2_mn(bb,b2tot)
      va2=b2tot*mu01*rho1
!
!  calculate JxB/rho (when hydro is on) and J^2 (when entropy is on)
!  add JxB/rho to momentum equation, and eta mu_0 J2/rho to entropy equation
!  set rhomin_JxB>0 in order to limit the JxB term at very low densities.
!  set va2max_JxB>0 in order to limit the JxB term at very high alven speeds.
!  set va2power_JxB to an integer value in order to specify the power
!  of the limiting term,
!
      if (lhydro) then
        call cross_mn(jj,bb,JxB)
        rho1_JxB=rho1
        if (rhomin_JxB>0) rho1_JxB=min(rho1_JxB,1/rhomin_JxB)
        if (va2max_JxB>0) rho1_JxB=rho1_JxB/(1+(va2/va2max_JxB)**va2power_JxB)
        call multsv_mn(rho1_JxB,JxB,JxBr)
        if(llorentzforce) df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+JxBr
        if(lentropy) then
          call dot2_mn(jj,J2)
          df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+(eta*mu0)*J2*rho1*TT1
        endif
      endif
!
!  calculate uxB+eta*del2A and add to dA/dt
!  (Note: the linear shear term is added later)
!  Read external electric field (currently for spheromak experiments)
!
      if (linduction) then
        call cross_mn(uu,bb,uxB)
        if (lee_ext) then
          call get_global(ee_ext,m,n,'ee_ext')
          uxB=uxB+ee_ext
        endif
      else
        uxB=0.
      endif
!
!  calculate restive term
!
      select case (iresistivity)

      case ('eta-const')
        fres=eta*del2A
        etatotal=eta
      case ('hyper3')
        call del6v(f,iaa,del6A)
        fres=eta*del6A
        etatotal=eta
      case ('hyper2')
        call del4v(f,iaa,del4A)
        fres=eta*del4A
        etatotal=eta  
      case ('shell')
        call eta_shell(eta_mn,geta)
        call div(f,iaa,divA)
        do j=1,3; fres(:,j)=eta_mn*del2A(:,j)+geta(:,j)*divA; enddo
        etatotal=eta_mn
      case ('shock')
        if (eta_shock/=0) then
          call div(f,iaa,divA)
          eta_tot=eta+eta_shock*shock
          geta=eta_shock*gshock
          do j=1,3; fres(:,j)=eta_tot*del2A(:,j)+geta(:,j)*divA; enddo
          etatotal=eta+eta_shock*shock
        else
          fres=eta*del2A
          etatotal=eta
        endif
      case ('Smagorinsky')
        call dot2_mn(jj,J2)
        eta_smag=(D_smag*dxmax)**2.*sqrt(J2)
        call multsv(eta_smag+eta,del2A,fres)
        etatotal=eta_smag+eta
      case ('Smagorinsky_cross')        
        oo(:,1)=uij(:,3,2)-uij(:,2,3)
        oo(:,2)=uij(:,1,3)-uij(:,3,1)
        oo(:,3)=uij(:,2,1)-uij(:,1,2)
        call dot(jj,oo,jo)
        sign_jo=1.
        do i=1,nx 
          if (jo(i) .lt. 0) sign_jo(i)=-1.
        enddo
        eta_smag=(D_smag*dxmax)**2.*sign_jo*sqrt(jo*sign_jo)
        call multsv(eta_smag+eta,del2A,fres)
        etatotal=eta_smag+eta
      case default
        if (lroot) print*,'daa_dt: no such ires:',iresistivity
        call stop_it("")
      end select
      if (headtt) print*,'daa_dt: iresistivity=',iresistivity
!
!  Switch off diffusion of horizontal components in boundary slice if
!  requested by boundconds
!
      if (lfrozen_bz_z_bot) then
        !
        ! Only need to do this for nonperiodic z direction, on bottommost
        ! processor and in bottommost pencils
        !
        if ((.not. lperi(3)) .and. (ipz == 0) .and. (n == n1)) then
          fres(:,1) = 0.
          fres(:,2) = 0.
        endif
      endif
!
!  Add to dA/dt
!
      df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+uxB+fres
!
!  Ambipolar diffusion in the strong coupling approximation
!
      if (nu_ni/=0.) then
        nu_ni1=1./nu_ni
        call cross_mn(JxBr,bb,JxBrxB)
        df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+nu_ni1*JxBrxB
        etatotal=etatotal+nu_ni1*va2
      endif
!
!  Hall term
!
      if (hall_term/=0.) then
        if (headtt) print*,'daa_dt: hall_term=',hall_term
        df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)-hall_term*JxB
        if (lfirst.and.ldt) then
          advec_hall=abs(uu(:,1)-hall_term*jj(:,1))*dx_1(l1:l2)+ &
                     abs(uu(:,2)-hall_term*jj(:,2))*dy_1(  m  )+ &
                     abs(uu(:,3)-hall_term*jj(:,3))*dz_1(  n  )
        endif
        if (headtt.or.ldebug) print*,'duu_dt: max(advec_hall) =',&
                                     maxval(advec_hall)
      endif
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
        eta_out1=eta_out*(1.-exp(-tmp**5/max(1.-tmp,1e-5)))-eta
        df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)-(eta_out1*mu0)*jj
      endif
!
!  possibility of relaxation of A in exterior region
!
      if (tau_aa_exterior/=0.) call calc_tau_aa_exterior(f,df)
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
          call dot_mn(jj,bb,jb)
          jb_yz(m-m1+1,n-n1+1)=jb(ix-l1+1)
          if (m.eq.iy)  jb_xz(:,n-n1+1)=jb
          if (n.eq.iz)  jb_xy(:,m-m1+1)=jb
          if (n.eq.iz2) jb_xy2(:,m-m1+1)=jb
          if(bthresh_per_brms/=0) call calc_bthresh
          call vecout(41,trim(directory)//'/bvec',bb,bthresh,nbvec)
        endif
!
!  phi-averages
!  Note that this does not necessarily happen with ldiagnos=.true.
!  This block must be done after jj has been calculated.
!
      if (l2davgfirst) then
        bx=bb(:,1)
        by=bb(:,2)
        bz=bb(:,3)
        call phisum_mn_name_rz(bx*pomx+by*pomy,i_brmphi)
        call phisum_mn_name_rz(bx*phix+by*phiy,i_bpmphi)
        call phisum_mn_name_rz(bz,i_bzmphi)
        call phisum_mn_name_rz(b2,i_b2mphi)
        if (i_jbmphi/=0) then
          call dot_mn(jj,bb,jb)
          call phisum_mn_name_rz(jb,i_jbmphi)
        endif
        if (any((/i_uxbrmphi,i_uxbpmphi,i_uxbzmphi/) /= 0)) then
          call cross_mn(uu,bb,uxb)
          call phisum_mn_name_rz(uxb(:,1)*pomx+uxb(:,2)*pomy,i_uxbrmphi)
          call phisum_mn_name_rz(uxb(:,1)*phix+uxb(:,2)*phiy,i_uxbpmphi)
          call phisum_mn_name_rz(uxb(:,3)                   ,i_uxbzmphi)
        endif
      endif
!
!  ``va^2/dx^2'' and ``eta/dx^2'' for timestep
!
      if (lfirst.and.ldt) then
        advec_va2=((bb(:,1)*dx_1(l1:l2))**2+ &
                   (bb(:,2)*dy_1(  m  ))**2+ &
                   (bb(:,3)*dz_1(  n  ))**2)*mu01*rho1
        diffus_eta=etatotal*dxyz_2
        if (ldiagnos.and.i_dteta/=0) then
          call max_mn_name(diffus_eta/cdtv,i_dteta,l_dt=.true.)
        endif
      endif
      if (headtt.or.ldebug) then
        print*,'duu_dt: max(advec_va2) =',maxval(advec_va2)
        print*,'duu_dt: max(diffus_eta) =',maxval(diffus_eta)
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
        if (i_dtb/=0) call max_mn_name(sqrt(advec_va2)/cdt,i_dtb,l_dt=.true.)
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
        if (i_jrms/=0 .or. i_jmax/=0 .or. i_j2m/=0 .or. i_jm2/=0 &
            .or. i_epsM/=0 .or. i_epsM_LES/=0) then
          call dot2_mn(jj,j2)
          if (i_j2m/=0) call sum_mn_name(j2,i_j2m)
          if (i_jm2/=0) call max_mn_name(j2,i_jm2)
          if (i_jrms/=0) call sum_mn_name(j2,i_jrms,lsqrt=.true.)
          if (i_jmax/=0) call max_mn_name(j2,i_jmax,lsqrt=.true.)
          if (i_epsM/=0 .and. iresistivity /= 'hyper3') &
              call sum_mn_name(eta*j2,i_epsM)
          if (i_epsM_LES/=0) call sum_mn_name(eta_smag*j2,i_epsM_LES)
        endif
        !
        ! epsM need del4A in cases with hyperresistivity
        !
        if (i_epsM/=0 .and. iresistivity == 'hyper3') then
          call del4v(f,iaa,del4A)
          call dot2_mn(del4A,del4A2)
          if (i_epsM/=0) call sum_mn_name(eta*del4A2,i_epsM)
        endif
        !
        !  calculate surface integral <2ExA>*dS
        !
        if (i_exaym2/=0) call helflux(aa,uxb,jj)
        !
        !  calculate surface integral <2ExJ>*dS
        !
        if (i_exjm2/=0) call curflux(uxb,jj)
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
        !     < E_2 u1,1 + E2 u3,3 - E1 u2,1 - E3 u3,2 >
        !     < E_3 u1,1 + E3 u2,2 - E1 u3,1 - E2 u2,3 > )
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
        !  < b2 b1,3 >
        !
        if (i_b2b13m/=0) then
          b2b13=bb(:,2)*bij(:,1,3)
          call sum_mn_name(b2b13,i_b2b13m)
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
    subroutine calculate_vars_magnetic(f,bb,jj,bij,aij,del2A,graddivA)
!
!  Calculation of bb
!
!  06-feb-04/bing: coded
!  14-mar-04/axel: allow external magnetic field to precess about z-axis
!  11-sep-04/axel: calculate bb from aij matrix
!
      use Cdata
      use Sub
      use Deriv
      use Global, only: get_global

      real, dimension (mx,my,mz,mvar+maux) :: f       
      real, dimension (nx,3,3) :: aij,bij
      real, dimension (nx,3) :: aa,bb,jj,del2A,graddivA
      real, dimension (nx,3) :: bb_ext,bb_ext_pot,ee_ext,jj_ext
      real :: B2_ext,c,s

      intent(in)  :: f
      intent(out) :: bb,jj,bij,aij,del2A,graddivA
!
!  The following routines calculates the gradient matrix of
!  the vector potential (needed to calculate bb) and bij
!  (needed for cosmic ray evolution and thermal conduction).
!
      call gij(f,iaa,aij,1)
      call bij_etc(f,iaa,bij,del2A,graddivA)
!
!  use aij to calculate bb, and
!  use bij to calculate jj
!
      aa=f(l1:l2,m,n,iax:iaz)
      call curl_mn(aij,bb,aa)
      call curl_mn(bij,jj,bb)
!
!  in spherical geometry, del2A is best written as graddivA-jj.
!  After that we can rescale jj by mu01.
!
      if (lspherical) del2A=graddivA-jj
      if (mu01/=1.) jj=mu01*jj
!
!  Note; for diagnostics purposes keep copy of original field
!
      if (ldiagnos) bbb=bb
!
!  possibility to add external field
!
      B2_ext=B_ext(1)**2+B_ext(2)**2+B_ext(3)**2
!
!  allow external field to precess about z-axis
!  with frequency omega_Bz_ext
!
      if (B2_ext/=0.) then
        if (omega_Bz_ext==0.) then
          B_ext_tmp=B_ext
        elseif (omega_Bz_ext/=0.) then
          c=cos(omega_Bz_ext*t)
          s=sin(omega_Bz_ext*t)
          B_ext_tmp(1)=B_ext(1)*c-B_ext(2)*s
          B_ext_tmp(2)=B_ext(1)*s+B_ext(2)*c
          B_ext_tmp(3)=B_ext(3)
        endif
!
!  add the external field
!
        if (B_ext(1)/=0.) bb(:,1)=bb(:,1)+B_ext_tmp(1)
        if (B_ext(2)/=0.) bb(:,2)=bb(:,2)+B_ext_tmp(2)
        if (B_ext(3)/=0.) bb(:,3)=bb(:,3)+B_ext_tmp(3)
        if (headtt) print*,'calculate_vars_magnetic: B_ext=',B_ext
        if (headtt) print*,'calculate_vars_magnetic: B_ext_tmp=',B_ext_tmp
      endif
!
!  add the external potential field
!
      if (lB_ext_pot) then
        call get_global(bb_ext_pot,m,n,'B_ext_pot')
        bb=bb+bb_ext_pot
      endif
!
!  add external B-field (currently for spheromak experiments)
!
      if (lbb_ext) then
        call get_global(bb_ext,m,n,'bb_ext')
        bb=bb+bb_ext
      endif
!
!  external current (currently for spheromak experiments)
!
      if (ljj_ext) then
        call get_global(ee_ext,m,n,'ee_ext')
        !call get_global(jj_ext,m,n,'jj_ext')
        !jj=jj+jj_ext
        jj=jj-ee_ext*displacement_gun
      endif

    endsubroutine calculate_vars_magnetic
!***********************************************************************
    subroutine eta_shell(eta_mn,geta)
!
!   24-nov-03/dave: coded 
!
      use Cdata
      use Sub, only: step, der_step
!
      real, dimension (nx) :: eta_mn
      real, dimension (nx) :: prof,eta_r
      real, dimension (nx,3) :: geta
      real :: d_int=0.,d_ext=0.
!
      eta_r=0.
!
      if (eta_int > 0.) d_int=eta_int-eta
      if (eta_ext > 0.) d_ext=eta_ext-eta
!
!     calculate steps in resistivity
!
      prof=step(r_mn,r_int,wresistivity)
      eta_mn=d_int*(1-prof)
      prof=step(r_mn,r_ext,wresistivity)
      eta_mn=eta+eta_mn+d_ext*prof
!
!     calculate radial derivative of steps and gradient of eta
!
      prof=der_step(r_mn,r_int,wresistivity)
      eta_r=-d_int*prof
      prof=der_step(r_mn,r_ext,wresistivity)
      eta_r=eta_r+d_ext*prof
      geta=evr*spread(eta_r,2,3)
!
    endsubroutine eta_shell
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
    subroutine curflux(uxb,jj)
!
!  current helicity flux (preliminary)
!
!  27-nov-03/axel: adapted from helflux
!
      use Cdata
      use Sub
!
      real, dimension (nx,3), intent(in) :: uxb,jj
      real, dimension (nx,3) :: ee
      real, dimension (nx) :: FCx,FCz
      real :: FC
!
      ee=eta*jj-uxb
!
!  calculate current helicity flux in the X and Z directions
!  Could speed up by only calculating here boundary points!
!
      FCx=2*(ee(:,2)*jj(:,3)-ee(:,3)*jj(:,2))*dsurfyz
      FCz=2*(ee(:,1)*jj(:,2)-ee(:,2)*jj(:,1))*dsurfxy
!
!  sum up contribution per pencil
!  and then stuff result into surf_mn_name for summing up all processors.
!
      FC=FCx(nx)-FCx(1)
      if(ipz==0       .and.n==n1) FC=FC-sum(FCz)
      if(ipz==nprocz-1.and.n==n2) FC=FC+sum(FCz)
      call surf_mn_name(FC,i_exjm2)
!
    endsubroutine curflux
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
      integer :: iname,inamez,ixy,ixz,irz
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
        i_bxpt=0; i_bypt=0; i_bzpt=0; i_epsM_LES=0
        i_aybym2=0; i_exaym2=0; i_exjm2=0
        i_brms=0; i_bmax=0; i_jrms=0; i_jmax=0; i_vArms=0; i_vAmax=0; i_dtb=0
        i_beta1m=0; i_beta1max=0
        i_bx2m=0; i_by2m=0; i_bz2m=0
        i_bxbym=0; i_bxbzm=0; i_bybzm=0; i_djuidjbim=0
        i_bxmz=0; i_bymz=0; i_bzmz=0; i_bmx=0; i_bmy=0; i_bmz=0
        i_bxmxy=0; i_bymxy=0; i_bzmxy=0
        i_bxmxz=0; i_bymxz=0; i_bzmxz=0
        i_uxbm=0; i_oxuxbm=0; i_jxbxbm=0.; i_gpxbm=0.; i_uxDxuxbm=0.
        i_uxbmx=0; i_uxbmy=0; i_uxbmz=0
        i_uxjm=0; i_ujxbm=0
        i_b2b13m=0
        i_brmphi=0; i_bpmphi=0; i_bzmphi=0; i_b2mphi=0; i_jbmphi=0
        i_uxbrmphi=0; i_uxbpmphi=0; i_uxbzmphi=0;
        i_dteta=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dteta',i_dteta)
        call parse_name(iname,cname(iname),cform(iname),'aybym2',i_aybym2)
        call parse_name(iname,cname(iname),cform(iname),'exaym2',i_exaym2)
        call parse_name(iname,cname(iname),cform(iname),'exjm2',i_exjm2)
        call parse_name(iname,cname(iname),cform(iname),'abm',i_abm)
        call parse_name(iname,cname(iname),cform(iname),'jbm',i_jbm)
        call parse_name(iname,cname(iname),cform(iname),'ubm',i_ubm)
        call parse_name(iname,cname(iname),cform(iname),'b2m',i_b2m)
        call parse_name(iname,cname(iname),cform(iname),'bm2',i_bm2)
        call parse_name(iname,cname(iname),cform(iname),'j2m',i_j2m)
        call parse_name(iname,cname(iname),cform(iname),'jm2',i_jm2)
        call parse_name(iname,cname(iname),cform(iname),'epsM',i_epsM)
        call parse_name(iname,cname(iname),cform(iname),'epsM_LES',i_epsM_LES)
        call parse_name(iname,cname(iname),cform(iname),'brms',i_brms)
        call parse_name(iname,cname(iname),cform(iname),'bmax',i_bmax)
        call parse_name(iname,cname(iname),cform(iname),'jrms',i_jrms)
        call parse_name(iname,cname(iname),cform(iname),'jmax',i_jmax)
        call parse_name(iname,cname(iname),cform(iname),'vArms',i_vArms)
        call parse_name(iname,cname(iname),cform(iname),'vAmax',i_vAmax)
        call parse_name(iname,cname(iname),cform(iname),'beta1m',i_beta1m)
        call parse_name(iname,cname(iname),cform(iname),'beta1max',i_beta1max)
        call parse_name(iname,cname(iname),cform(iname),'dtb',i_dtb)
        call parse_name(iname,cname(iname),cform(iname),'bx2m',i_bx2m)
        call parse_name(iname,cname(iname),cform(iname),'by2m',i_by2m)
        call parse_name(iname,cname(iname),cform(iname),'bz2m',i_bz2m)
        call parse_name(iname,cname(iname),cform(iname),'bxbym',i_bxbym)
        call parse_name(iname,cname(iname),cform(iname),'bxbzm',i_bxbzm)
        call parse_name(iname,cname(iname),cform(iname),'bybzm',i_bybzm)
        call parse_name(iname,cname(iname),cform(iname),'djuidjbim',i_djuidjbim)
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
        call parse_name(iname,cname(iname),cform(iname),'b2b13m',i_b2b13m)
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
!  check for those quantities for which we want y-averages
!
      do ixz=1,nnamexz
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bxmxz',i_bxmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bymxz',i_bymxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bzmxz',i_bzmxz)
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
        call parse_name(irz,cnamerz(irz),cformrz(irz),'brmphi'  ,i_brmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'bpmphi'  ,i_bpmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'bzmphi'  ,i_bzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'b2mphi'  ,i_b2mphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'jbmphi'  ,i_jbmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'uxbrmphi',i_uxbrmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'uxbpmphi',i_uxbpmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'uxbzmphi',i_uxbzmphi)
      enddo
!
!  write column, i_XYZ, where our variable XYZ is stored
!
      if (lwr) then
        write(3,*) 'i_dteta=',i_dteta
        write(3,*) 'i_aybym2=',i_aybym2
        write(3,*) 'i_exaym2=',i_exaym2
        write(3,*) 'i_exjm2=',i_exjm2
        write(3,*) 'i_abm=',i_abm
        write(3,*) 'i_jbm=',i_jbm
        write(3,*) 'i_ubm=',i_ubm
        write(3,*) 'i_b2m=',i_b2m
        write(3,*) 'i_bm2=',i_bm2
        write(3,*) 'i_j2m=',i_j2m
        write(3,*) 'i_jm2=',i_jm2
        write(3,*) 'i_epsM=',i_epsM
        write(3,*) 'i_epsM_LES=',i_epsM_LES
        write(3,*) 'i_brms=',i_brms
        write(3,*) 'i_bmax=',i_bmax
        write(3,*) 'i_jrms=',i_jrms
        write(3,*) 'i_jmax=',i_jmax
        write(3,*) 'i_vArms=',i_vArms
        write(3,*) 'i_vAmax=',i_vAmax
        write(3,*) 'i_beta1m=',i_beta1m
        write(3,*) 'i_beta1max=',i_beta1max
        write(3,*) 'i_dtb=',i_dtb
        write(3,*) 'i_bx2m=',i_bx2m
        write(3,*) 'i_by2m=',i_by2m
        write(3,*) 'i_bz2m=',i_bz2m
        write(3,*) 'i_bxbym=',i_bxbym
        write(3,*) 'i_bxbzm=',i_bxbzm
        write(3,*) 'i_bybzm=',i_bybzm
        write(3,*) 'i_djuidjbim=',i_djuidjbim
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
        write(3,*) 'i_b2b13m=',i_b2b13m
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
        write(3,*) 'nnamexz=',nnamexz
        write(3,*) 'i_bxmxy=',i_bxmxy
        write(3,*) 'i_bymxy=',i_bymxy
        write(3,*) 'i_bzmxy=',i_bzmxy
        write(3,*) 'i_bxmxz=',i_bxmxz
        write(3,*) 'i_bymxz=',i_bymxz
        write(3,*) 'i_bzmxz=',i_bzmxz
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
!  ux and Ay.
!  Don't overwrite the density, just add to the log of it.
!
      f(:,:,:,ilnrho)=ampl*sin(kx*xx)+f(:,:,:,ilnrho)
      f(:,:,:,iuu+0)=+ampl*sin(kx*xx)
      f(:,:,:,iuu+1)=+ampl*sin(kx*xx)
      f(:,:,:,iaa+2)=-ampl*cos(kx*xx)
!
    endsubroutine alfven_x
!***********************************************************************
    subroutine alfven_z(ampl,f,iuu,iaa,zz,kz,mu0)
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
      real :: ampl,kz,mu0
      integer :: iuu,iaa
!
!  ux and Ay
!
      f(:,:,:,iuu+0)=+ampl*cos(kz*zz)
      f(:,:,:,iaa+1)=+ampl*sin(kz*zz)*sqrt(mu0)
!
    endsubroutine alfven_z
!***********************************************************************
    subroutine alfvenz_rot(ampl,f,iuu,iaa,zz,kz,O)
!
!  Alfven wave propagating in the z-direction (with Coriolis force)
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
      use Cdata, only: lroot
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: zz
      real :: ampl,kz,O,fac
      integer :: iuu,iaa
!
!  ux, uy, Ax and Ay
!
      if (lroot) print*,'alfvenz_rot: Alfven wave with rotation; O,kz=',O,kz
      fac=-O+sqrt(O**2+kz**2)
      f(:,:,:,iuu+0)=-ampl*sin(kz*zz)*fac/kz
      f(:,:,:,iuu+1)=-ampl*cos(kz*zz)*fac/kz
      f(:,:,:,iaa+0)=+ampl*sin(kz*zz)/kz
      f(:,:,:,iaa+1)=+ampl*cos(kz*zz)/kz
!
    endsubroutine alfvenz_rot
!***********************************************************************
    subroutine alfvenz_rot_shear(ampl,f,iuu,iaa,zz,kz,O)
!
!  Alfven wave propagating in the z-direction (with Coriolis force and shear)
!
!  satisfies the equations
!  dux/dt - 2*Omega*uy = -Ay''
!  duy/dt + 1/2*Omega*ux = +Ax''
!  dAx/dt = 3/2*Omega*Ay + uy
!  dAy/dt = -ux
!
!  Assume B0=rho0=mu0=1
!
!  28-june-04/anders: coded
!
      use Cdata, only: lroot
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: zz
      real :: ampl,kz,O
      complex :: fac
      integer :: iuu,iaa
!
!  ux, uy, Ax and Ay
!
      if (lroot) print*,'alfvenz_rot_shear: '// &
          'Alfven wave with rotation and shear; O,kz=',O,kz
      fac=cmplx(O-sqrt(16*kz**2+O**2),0.)
      f(:,:,:,iuu+0)=f(:,:,:,iuu+0) + ampl*fac/(4*kz)*sin(kz*zz)
      f(:,:,:,iuu+1)=f(:,:,:,iuu+1) + ampl*real(exp(cmplx(0,zz*kz))* &
          fac*sqrt(2*kz**2+O*fac)/(sqrt(2.)*kz*(-6*O-fac)))
      f(:,:,:,iaa+0)=ampl*sin(kz*zz)/kz
      f(:,:,:,iaa+1)=-ampl*2*sqrt(2.)*aimag(exp(cmplx(0,zz*kz))* &
          sqrt(2*kz**2+O*fac)/(-6*O-fac)/(cmplx(0,kz)))
!
    endsubroutine alfvenz_rot_shear
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
    subroutine force_free_jet(mu,xx,yy,zz)
!
!  Force free magnetic field configuration for jet simulations
!  with a fixed accretion disk at the bottom boundary.
!
!  The input parameter mu specifies the radial dependency of
!  the magnetic field in the disk.
!
!  Solves the laplace equation in cylindrical coordinates for the
!  phi-component of the vector potential. A_r and A_z are taken to
!  be zero.
!
!    nabla**2 A_phi - A_phi / r**2 = 0
!
!  For the desired boundary condition in the accretion disk
!
!    B_r=B0*r**(mu-1)  (z == 0)
!
!  the solution is
!
!    A_phi = Hypergeometric2F1( (1-mu)/2, (2+mu)/2, 2, xi**2 )
!            *xi*(r**2+z**2)**(mu/2)
!
!  where xi = sqrt(r**2/(r**2+z**2))
!
!
!  30-may-04/tobi: coded
!
      use Cdata, only: x,y,z,lroot,directory,ip,m,n,pi,r_ref
      use Sub, only: hypergeometric2F1,gamma_function
      use Global, only: set_global
      use Deriv, only: der
      use Io, only: output

      real, intent(in) :: mu
      real, dimension(mx,my,mz), intent(in) :: xx,yy,zz
      real :: xi2,A_phi
      real :: a,b,c,r2
      real :: B1r_,B1z_,B1
      real, parameter :: tol=10*epsilon(1.0)
      integer :: l,ierr
      integer, dimension(1) :: ll,mm,nn
      real, dimension(mx,my,mz) :: Ax_ext,Ay_ext
      real, dimension(nx,3) :: bb_ext_pot
      real, dimension(nx) :: bb_x,bb_y,bb_z
!
!  calculate un-normalized |B| at r=r_ref and z=0 for later normalization
!
      if (lroot.and.ip<=5) print*,'FORCE_FREE_JET: calculating normalization'

      B1r_=sin(pi*mu/2)*gamma_function(   abs(mu) /2) / &
                        gamma_function((1+abs(mu))/2)

      B1z_=cos(pi*mu/2)*gamma_function((1+abs(mu))/2) / &
                        gamma_function((2+abs(mu))/2)

      B1=sqrt(4/pi)*r_ref**(mu-1)*sqrt(B1r_**2+B1z_**2)
!
!  calculate external vector potential
!
      if (lroot) print*,'FORCE_FREE_JET: calculating external vector potential'

      if (lforce_free_test) then

        if (lroot) print*,'FORCE_FREE_JET: using analytic solution for mu=-1'
        Ax_ext=-2*yy*(1-zz/sqrt(xx**2+yy**2+zz**2))/(xx**2+yy**2)/B1
        Ay_ext= 2*xx*(1-zz/sqrt(xx**2+yy**2+zz**2))/(xx**2+yy**2)/B1

      else

        do l=1,mx
        do m=1,my
        do n=1,mz

          r2=x(l)**2+y(m)**2
          xi2=r2/(r2+z(n)**2)
          A_phi=hypergeometric2F1((1-mu)/2,(2+mu)/2,2.0,xi2,tol) &
               *sqrt(xi2)*sqrt(r2+z(n)**2)**mu/B1

          Ax_ext(l,m,n)=-y(m)*A_phi/sqrt(r2)
          Ay_ext(l,m,n)= x(l)*A_phi/sqrt(r2)

        enddo
        enddo
        enddo

      endif

!
!  calculate external magnetic field
!
      if (lroot.and.ip<=5) &
        print*,'FORCE_FREE_JET: calculating the external magnetic field'

      do n=n1,n2
      do m=m1,m2
        call der(Ay_ext,bb_x,3)
        bb_ext_pot(:,1)=-bb_x
        call der(Ax_ext,bb_y,3)
        bb_ext_pot(:,2)= bb_y
        call der(Ay_ext,bb_z,1)
        bb_ext_pot(:,3)= bb_z
        call der(Ax_ext,bb_z,2)
        bb_ext_pot(:,3)=bb_ext_pot(:,3)-bb_z
        call set_global(bb_ext_pot,m,n,'B_ext_pot')
      enddo
      enddo

      if (ip<=5) then
        call output(trim(directory)//'/Ax_ext.dat',Ax_ext,1)
        call output(trim(directory)//'/Ay_ext.dat',Ay_ext,1)
      endif

    endsubroutine force_free_jet
!***********************************************************************
    subroutine geo_benchmark_B(f)
!
!  30-june-04/grs: coded
!
      use Cdata
      use Sub, only: calc_unitvects_sphere
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mvar+maux), intent(inout) :: f     
      real, dimension(nx) :: theta_mn,ar,atheta,aphi
      real :: C_int,C_ext,A_int,A_ext

      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        call calc_unitvects_sphere()
        theta_mn=acos(z_mn/r_mn)
        phi_mn=atan2(y_mn,x_mn)

! calculate ax,ay,az (via ar,atheta,aphi) inside shell (& leave zero outside shell)
          select case(initaa) 
            case('geo-benchmark-case1')
              if (lroot .and. imn==1) print*, 'geo_benchmark_B: geo-benchmark-case1'
              C_int=-( -1./63.*r_int**4 + 11./84.*r_int**3*r_ext             &
                     + 317./1050.*r_int**2*r_ext**2                         &
                     - 1./5.*r_int**2*r_ext**2*log(r_int) )
              C_ext=-( -1./63.*r_ext**9 + 11./84.*r_ext**8*r_int             &
                     + 317./1050.*r_ext**7*r_int**2                         &
                     - 1./5.*r_ext**7*r_int**2*log(r_ext) )
              A_int=5./2.*(r_ext-r_int)
              A_ext=5./8.*(r_ext**4-r_int**4)

              where (r_mn < r_int)
                ar=C_int*ampl_B0*80.*2.*(3.*sin(theta_mn)**2-2.)*r_mn
                atheta=3.*C_int*ampl_B0*80.*sin(2.*theta_mn)*r_mn 
                aphi=ampl_B0*A_int*r_mn*sin(theta_mn)
              endwhere

              where (r_mn <= r_ext .and. r_mn >= r_int) 
                ar=ampl_B0*80.*2.*(3.*sin(theta_mn)**2-2.)*                 &
                   (   1./36.*r_mn**5 - 1./12.*(r_int+r_ext)*r_mn**4        &
                     + 1./14.*(r_int**2+4.*r_int*r_ext+r_ext**2)*r_mn**3    &
                     - 1./3.*(r_int**2*r_ext+r_int*r_ext**2)*r_mn**2        &
                     - 1./25.*r_int**2*r_ext**2*r_mn                        &
                     + 1./5.*r_int**2*r_ext**2*r_mn*log(r_mn) )
                atheta=-ampl_B0*80.*sin(2.*theta_mn)*                        &
                   (   7./36.*r_mn**5 - 1./2.*(r_int+r_ext)*r_mn**4         &
                     + 5./14.*(r_int**2+4.*r_int*r_ext+r_ext**2)*r_mn**3    &
                     - 4./3.*(r_int**2*r_ext+r_int*r_ext**2)*r_mn**2        &
                     + 2./25.*r_int**2*r_ext**2*r_mn                        &
                     + 3./5.*r_int**2*r_ext**2*r_mn*log(r_mn) )
                aphi=ampl_B0*5./8.*sin(theta_mn)*                           &
                   ( 4.*r_ext*r_mn - 3.*r_mn**2 - r_int**4/r_mn**2 ) 
              endwhere

              where (r_mn > r_ext)
                ar=C_ext*ampl_B0*80.*2.*(3.*sin(theta_mn)**2-2.)/r_mn**4
                atheta=-2.*C_ext*ampl_B0*80.*sin(2.*theta_mn)/r_mn**4
                aphi=ampl_B0*A_ext/r_mn**2*sin(theta_mn)
              endwhere
  
          ! debug checks -- look at a pencil near the centre...
          if (ip<=4 .and. imn==(ny+1)*nz/2) then
            print*,'r_int,r_ext',r_int,r_ext
            print'(a45,2i6,2f15.7)','geo_benchmark_B: minmax(r_mn), imn, iproc:', iproc, imn, minval(r_mn), maxval(r_mn)
            print'(a45,2i6,2f15.7)','geo_benchmark_B: minmax(theta_mn), imn, iproc:', iproc, imn, minval(theta_mn), maxval(theta_mn)
            print'(a45,2i6,2f15.7)','geo_benchmark_B: minmax(phi_mn), imn, iproc:', iproc, imn, minval(phi_mn), maxval(phi_mn)
            print'(a45,2i6,2f15.7)','geo_benchmark_B: minmax(ar), imn, iproc:', iproc, imn, minval(ar), maxval(ar)
            print'(a45,2i6,2f15.7)','geo_benchmark_B: minmax(atheta), imn, iproc:', iproc, imn, minval(atheta), maxval(atheta)
            print'(a45,2i6,2f15.7)','geo_benchmark_B: minmax(aphi), imn, iproc:', iproc, imn, minval(aphi), maxval(aphi)
          endif

            case('geo-benchmark-case2')
              if (lroot .and. imn==1) print*, 'geo_benchmark_B: geo-benchmark-case2 not yet coded.'

            case default
              if (lroot .and. imn==1) print*,'geo_benchmark_B: case not defined!'
              call stop_it("")
          endselect

          f(l1:l2,m,n,iax)=sin(theta_mn)*cos(phi_mn)*ar + cos(theta_mn)*cos(phi_mn)*atheta - sin(phi_mn)*aphi
          f(l1:l2,m,n,iay)=sin(theta_mn)*sin(phi_mn)*ar + cos(theta_mn)*sin(phi_mn)*atheta + cos(phi_mn)*aphi
          f(l1:l2,m,n,iaz)=cos(theta_mn)*ar - sin(theta_mn)*atheta
      enddo

      if (ip<=14) then
        print*,'geo_benchmark_B: minmax(ax) on iproc:', iproc, minval(f(l1:l2,m1:m2,n1:n2,iax)),maxval(f(l1:l2,m1:m2,n1:n2,iax))
        print*,'geo_benchmark_B: minmax(ay) on iproc:', iproc, minval(f(l1:l2,m1:m2,n1:n2,iay)),maxval(f(l1:l2,m1:m2,n1:n2,iay))
        print*,'geo_benchmark_B: minmax(az) on iproc:', iproc, minval(f(l1:l2,m1:m2,n1:n2,iaz)),maxval(f(l1:l2,m1:m2,n1:n2,iaz))
      endif

    endsubroutine geo_benchmark_B
    
!***********************************************************************
    subroutine bc_frozen_in_bb_z(topbot)
!
!  Set flags to indicate that magnetic flux is frozen-in at the
!  z boundary. The implementation occurs in daa_dt where magnetic
!  diffusion is switched off in that layer.
!
      use Cdata
!
      character (len=3) :: topbot
!
      select case(topbot)
      case('bot')               ! bottom boundary
        lfrozen_bz_z_bot = .true.    ! set flag
      case('top')               ! top boundary
        lfrozen_bz_z_top = .true.    ! set flag
      case default
        print*, "bc_frozen_in_bb_z: ", topbot, " should be `top' or `bot'"
      endselect
!
    endsubroutine bc_frozen_in_bb_z
!***********************************************************************
      subroutine bc_aa_pot(f,topbot)
!
!  Potential field boundary condition for magnetic vector potential at
!  bottom or top boundary (in z).
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
