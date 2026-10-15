! $Id$
!
!  Initial data for the relativistic fluid and magnetic field in curved
!  spacetime (hydro with lgr_metric=T) together with the Z4c metric of
!  special/z4c, which are overwritten here:
!
!  'flrw':         homogeneous fluid in a spatially flat FLRW universe,
!                  synchronous gauge (alpha=1, beta=0), a=1 and H=H0 at t=0;
!                  with a uniform magnetic field B_ext a Bianchi I universe.
!  'flrw_perturb': as 'flrw', plus a linear adiabatic plane-wave mode along x
!                  in synchronous gauge with the amplitudes delta0, theta0,
!                  eta0 of Ma & Bertschinger (1995, ApJ 455, 7) and h=0 at t=0;
!                  hdot and etadot follow from the constraints.
!  'flrw_singlemode': linearly perturbed dust (or p=cs20*e) FLRW universe of
!                  FLRWSolver (Macpherson, Lasky & Price 2017, PRD 95, 064028),
!                  as in its flrw_singlemode.F90: longitudinal-gauge growing
!                  mode with Phi = phi0 sin(k x) (perturb_dir='x') or
!                  phi0 sum_i sin(k x^i) ('all'), k = 2 pi/Lx, H = H0, a=1,
!                  lapse alpha0*sqrt(1+2 Phi).
!  'alfven':       large-amplitude circularly polarized Alfven wave along x
!                  (exact solution, e.g. Del Zanna et al. 2007, A&A 473, 11)
!                  in flat spacetime with constant lapse alpha0, scale factor
!                  a0 and shift beta0, or, if H0/=0, in an FLRW universe.
!                  The guide field is taken from B_ext(1) = sqrt(gamma) B^x.
!  'gw_standing':  linear standing gravitational wave in TT gauge in an
!                  initially homogeneous, uniformly magnetized fluid at rest,
!                  test D of Duez, Liu, Shapiro & Stephens (2005, PRD 72,
!                  024028) with the analytic solution of Paper II (PRD 72,
!                  024029). Their z axis is x here (x,y,z -> y,z,x):
!                  gamma_yy = 1 + h+, gamma_zz = 1 - h+, gamma_yz = hx with
!                  h+,x = (hplus0, hcross0) sin(k x) cos(k t), K_ij = 0 at t=0;
!                  fluid energy density e0, field sqrt(gamma) B^i = B_ext.
!
!  13-oct-26: coded
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
module InitialCondition
!
  use Cdata
  use Quiet
  use Messages
!
  implicit none
!
  include '../initial_condition.h'
!
  character(len=1), dimension(3), parameter :: cxyz=(/'x','y','z'/)
!
  character(len=labellen) :: initgr='flrw'
  real :: H0=0., e0=1., alpha0=1., a0=1., delta0=0., theta0=0., eta0=0., eta_alfven=0.
  real, dimension(3) :: beta0=0.
  real :: phi0=0., hplus0=0., hcross0=0.
  integer :: nkgr=1
  character(len=labellen) :: perturb_dir='x'
!
  namelist /initial_condition_pars/ initgr, H0, e0, alpha0, a0, beta0, nkgr, &
      delta0, theta0, eta0, eta_alfven, phi0, perturb_dir, hplus0, hcross0
!
  integer :: ichi=0, iKh=0, ith=0, ialp=0
  integer, dimension(3) :: iGam=0, ibeta=0
  integer, dimension(6) :: igt=0, iAt=0
  real, dimension(:), pointer :: B_ext
  real :: kgr
!
  contains
!***********************************************************************
    subroutine register_initial_condition
!
      call svn_id("$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  13-oct-26/TP: coded
!
      use FArrayManager, only: farray_index_by_name
      use SharedVariables, only: get_shared_variable
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i, j
!
      ichi=farray_index_by_name('chi')
      iKh=farray_index_by_name('Kh')
      ith=farray_index_by_name('theta')
      ialp=farray_index_by_name('alp')
      do i=1,3
        iGam(i)=farray_index_by_name('Gam'//cxyz(i))
        ibeta(i)=farray_index_by_name('beta'//cxyz(i))
      enddo
      do j=1,3; do i=1,j
        igt(ij_table(i,j))=farray_index_by_name('gt'//cxyz(i)//cxyz(j))
        iAt(ij_table(i,j))=farray_index_by_name('At'//cxyz(i)//cxyz(j))
      enddo; enddo
      if (any((/ichi,iKh,ith,ialp,iGam,ibeta,igt,iAt/)<=0)) &
        call fatal_error('initialize_initial_condition','gr_matter needs special/z4c')
      if (lmagnetic) then
        call get_shared_variable('B_ext',B_ext,caller='initialize_initial_condition')
      else
        allocate(B_ext(3)); B_ext=0.
      endif
!
      kgr=2*pi*nkgr/Lxyz(1)
      if (H0/=0.) e0=3*H0**2/(8*pi)
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_all(f,profiles)
!
!  13-oct-26/TP: coded
!
      use EquationOfState, only: cs20
      use Sub, only: sym_inverse
!
      real, dimension(mx,my,mz,mfarray), optional, intent(inout) :: f
      real, dimension(:,:), optional, intent(out) :: profiles
!
      real, parameter :: eps=1e-5
      real, dimension(mx,6) :: g, Kij, gu, gp, gm
      real, dimension(mx,3) :: vu, Bu, Atil
      real, dimension(mx) :: e, detg, ph, chi, gtup, gtum
      real :: hdot, etadot, cA, ebg, Etot, uboost, etaL, Sperp
      integer :: i, it
!
      if (present(profiles)) call fatal_error('initial_condition_all','profiles not implemented')
!
      if (initgr=='flrw_singlemode') then
        call flrw_singlemode(f)
        return
      endif
!
      ph=kgr*x
      g=0.; Kij=0.; vu=0.; Atil=0.
      do i=1,3
        g(:,ij_table(i,i))=a0**2
        Kij(:,ij_table(i,i))=-a0**2*H0
      enddo
      e=e0
      etaL=eta_alfven
!
      select case (initgr)
!
!
!  With a uniform field B_ext, the fluid energy density is reduced by the
!  magnetic energy density B^2/2 so that the Hamiltonian constraint holds;
!  for B_ext/=0 the universe is then of Bianchi type I.
!
      case ('flrw')
        if (lmagnetic) e=e0-.5*sum(B_ext**2)/a0**4
!
!  Linear adiabatic mode, metric gamma_ij = delta_ij + h_ij with
!  h_ij = cos(kx) diag(h+4eta,-2eta,-2eta) and h(0)=0;
!  K_ij = -H gamma_ij - hdot_ij/2 (synchronous gauge, a=1).
!
      case ('flrw_perturb')
        if (H0==0. .or. a0/=1. .or. alpha0/=1. .or. any(beta0/=0.)) call fatal_error( &
            'initial_condition_all','flrw_perturb needs H0/=0, a0=1, alpha0=1, beta0=0')
        hdot=2*(kgr**2*eta0+4*pi*e0*delta0)/H0
        etadot=4*pi*(1+cs20)*e0*theta0/kgr**2
        g(:,i11)=g(:,i11)+4*eta0*cos(ph)
        g(:,i22)=g(:,i22)-2*eta0*cos(ph)
        g(:,i33)=g(:,i33)-2*eta0*cos(ph)
        Kij(:,i11)=-H0*g(:,i11)-.5*(hdot+4*etadot)*cos(ph)
        Kij(:,i22)=-H0*g(:,i22)+etadot*cos(ph)
        Kij(:,i33)=-H0*g(:,i33)+etadot*cos(ph)
        e=e0*(1+delta0*cos(ph))
        if (lmagnetic) e=e-.5*sum(B_ext**2)/a0**4
        vu(:,1)=theta0/kgr*sin(ph)
!
!  Circularly polarized Alfven wave. In the fluid frame (where v^x=0),
!  B = B0 (1, eta cos(kx), eta sin(kx)) and v = -c eta (0, cos(kx), sin(kx)),
!  with c^2 = [Q - sqrt(Q^2-4 eta^2 B0^4)]/(2 eta^2 B0^2), Q = (e+p) + B0^2 (1+eta^2).
!  B0 is B_ext(1)/a0^2, the physical guide field.
!
!  In FLRW (H0/=0) the uniform momentum density S_x = c eta^2 B0^2 of this
!  state would violate the momentum constraint in a periodic box, so the
!  solution is boosted along x with the velocity u for which S_x=0; it then
!  has v^x=-u, transverse amplitude eta*gamma_u*(1-u c) and phase speed
!  (c-u)/(1-u c). The momentum constraint is satisfied by K_xy and K_xz,
!  and the background fluid energy density is chosen such that the
!  Hamiltonian constraint holds, including the uniform A_ij A^ij of K_xy, K_xz.
!
      case ('alfven')
        if (.not.lmagnetic) call fatal_error('initial_condition_all','alfven needs magnetic')
        if (H0/=0. .and. (a0/=1. .or. alpha0/=1. .or. any(beta0/=0.))) call fatal_error( &
            'initial_condition_all','alfven with H0/=0 needs a0=1, alpha0=1, beta0=0')
        ebg=e0
        do it=1,100
          call alfven_state(ebg,cs20,B_ext(1)/a0**2,eta_alfven,H0/=0.,cA,uboost,etaL,Etot,Sperp)
          if (H0==0.) exit
          ebg=ebg+((6*H0**2-2*(8*pi*Sperp/kgr)**2)/(16*pi)-Etot)/(1+cs20)
        enddo
        if (lroot) print*,'initial_condition_all: alfven: c, u, (c-u)/(1-uc), e=', &
            cA, uboost, (cA-uboost)/(1-uboost*cA), ebg
        e=ebg
        vu(:,1)=-uboost/a0
        vu(:,2)=-cA*eta_alfven*sqrt(1-uboost**2)*cos(ph)/a0
        vu(:,3)=-cA*eta_alfven*sqrt(1-uboost**2)*sin(ph)/a0
        Atil(:,2)=-B_ext(1)*etaL*cos(ph)/kgr
        Atil(:,3)=-B_ext(1)*etaL*sin(ph)/kgr
!
!  Standing gravitational wave at t=0 (K_ij=0) in a fluid at rest.
!
      case ('gw_standing')
        if (H0/=0. .or. a0/=1. .or. alpha0/=1. .or. any(beta0/=0.)) call fatal_error( &
            'initial_condition_all','gw_standing needs H0=0, a0=1, alpha0=1, beta0=0')
        call metric_at(x,g)
!
      case default
        call fatal_error('initial_condition_all','no such initgr: '//trim(initgr))
      endselect
!
!  Z4c variables; Gam^i = -d_j gt^ij by centered differences of the
!  analytic gt^xx (only x dependence).
!
      call sym_inverse(g,gu,detg)
      chi=detg**(-1./3.)
      call metric_at(x+eps,gp); call metric_at(x-eps,gm)
      call sym_inverse(gp,gu,detg)
      gtup=gu(:,i11)*detg**(1./3.)
      call sym_inverse(gm,gu,detg)
      gtum=gu(:,i11)*detg**(1./3.)
      call sym_inverse(g,gu,detg)
!
      do n=1,mz; do m=1,my
        f(:,m,n,ichi)=chi
        f(:,m,n,iKh)=contract_sym(gu,Kij)
        f(:,m,n,ith)=0.
        f(:,m,n,ialp)=alpha0
        do i=1,6
          f(:,m,n,igt(i))=chi*g(:,i)
          f(:,m,n,iAt(i))=chi*(Kij(:,i)-g(:,i)*contract_sym(gu,Kij)/3.)
        enddo
        f(:,m,n,iGam(1))=-(gtup-gtum)/(2*eps)
        f(:,m,n,iGam(2))=0.
        f(:,m,n,iGam(3))=0.
        do i=1,3
          f(:,m,n,ibeta(i))=beta0(i)
        enddo
!
!  Vector potential and densitized field sqrt(gamma) B^i = B_ext + curl A.
!
        if (lmagnetic) then
          f(:,m,n,iax:iaz)=Atil
          Bu(:,1)=B_ext(1)
          Bu(:,2)=B_ext(2)+B_ext(1)*etaL*cos(ph)
          Bu(:,3)=B_ext(3)+B_ext(1)*etaL*sin(ph)
          if (initgr/='alfven') Bu(:,2:3)=spread(B_ext(2:3),1,mx)
          Bu=Bu/spread(sqrt(1./chi**3),2,3)
        else
          Bu=0.
        endif
        call prim2con(g,gu,e,vu,Bu,cs20,f(:,m,n,irho),f(:,m,n,iux:iuz))
!
!  Momentum constraint d_x K_xy = 8 pi S_y, d_x K_xz = 8 pi S_z for the
!  transverse momentum ~ (cos(kx), sin(kx)) of the Alfven wave in FLRW (a=1).
!
        if (initgr=='alfven' .and. H0/=0.) then
          f(:,m,n,iAt(i12))=8*pi*f(:,m,n,iuz)/kgr
          f(:,m,n,iAt(i13))=-8*pi*f(:,m,n,iuy)/kgr
        endif
      enddo; enddo
!
    endsubroutine initial_condition_all
!***********************************************************************
    subroutine flrw_singlemode(f)
!
!  Single-mode initial data of FLRWSolver (flrw_singlemode.F90):
!  gamma_ij = a^2 (1-2 Phi) delta_ij, alpha = alpha0 sqrt(1+2 Phi),
!  K_ij = -a adot (1-2 Phi)/alpha delta_ij with adot = H a,
!  e = rho0 (1 + delta), delta = (-2 k^2/(3 H^2) - 2) Phi,
!  v^i = -2/(3 a H) d_i Phi, with rho0 = 3 H^2/(8 pi a^2).
!  The conformal metric is flat, so Gam^i=0.
!
!  14-oct-26/TP: coded
!
      use EquationOfState, only: cs20
      use Sub, only: sym_inverse
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      real, dimension(mx,6) :: g, gu
      real, dimension(mx,3) :: vu, Bu
      real, dimension(mx) :: phi, alp, e, detg, Kd
      real :: rho0, pert_rho, pert_v
      integer :: i
!
      if (H0==0.) call fatal_error('flrw_singlemode','need H0/=0')
      rho0=3*H0**2/(8*pi*a0**2)
      pert_rho=-2*kgr**2/(3*H0**2)-2
      pert_v=-2/(3*a0*H0)
      Bu=0.
      if (lmagnetic) call fatal_error('flrw_singlemode','no magnetic field implemented')
!
      do n=1,mz; do m=1,my
        vu=0.
        select case (perturb_dir)
        case ('x')
          phi=phi0*sin(kgr*x)
          vu(:,1)=pert_v*phi0*kgr*cos(kgr*x)
        case ('all')
          phi=phi0*(sin(kgr*x)+sin(kgr*y(m))+sin(kgr*z(n)))
          vu(:,1)=pert_v*phi0*kgr*cos(kgr*x)
          vu(:,2)=pert_v*phi0*kgr*cos(kgr*y(m))
          vu(:,3)=pert_v*phi0*kgr*cos(kgr*z(n))
        case default
          call fatal_error('flrw_singlemode','no such perturb_dir: '//trim(perturb_dir))
        endselect
        alp=alpha0*sqrt(1+2*phi)
        g=0.
        do i=1,3
          g(:,ij_table(i,i))=a0**2*(1-2*phi)
        enddo
        Kd=-H0*a0**2*(1-2*phi)/alp
        e=rho0*(1+pert_rho*phi)
        call sym_inverse(g,gu,detg)
        f(:,m,n,ichi)=detg**(-1./3.)
        f(:,m,n,iKh)=3*Kd/(a0**2*(1-2*phi))
        f(:,m,n,ith)=0.
        f(:,m,n,ialp)=alp
        do i=1,6
          f(:,m,n,igt(i))=0.
          f(:,m,n,iAt(i))=0.
        enddo
        do i=1,3
          f(:,m,n,igt(ij_table(i,i)))=1.
          f(:,m,n,iGam(i))=0.
          f(:,m,n,ibeta(i))=beta0(i)
        enddo
        call prim2con(g,gu,e,vu,Bu,cs20,f(:,m,n,irho),f(:,m,n,iux:iuz))
      enddo; enddo
!
    endsubroutine flrw_singlemode
!***********************************************************************
    subroutine metric_at(xx,g)
!
!  Spatial metric as a function of x (for the derivative of gt^ij).
!
      real, dimension(:) :: xx
      real, dimension(:,:) :: g
      integer :: i
!
      g=0.
      do i=1,3
        g(:,ij_table(i,i))=a0**2
      enddo
      if (initgr=='flrw_perturb') then
        g(:,i11)=g(:,i11)+4*eta0*cos(kgr*xx)
        g(:,i22)=g(:,i22)-2*eta0*cos(kgr*xx)
        g(:,i33)=g(:,i33)-2*eta0*cos(kgr*xx)
      elseif (initgr=='gw_standing') then
        g(:,i22)=g(:,i22)+hplus0*sin(kgr*xx)
        g(:,i33)=g(:,i33)-hplus0*sin(kgr*xx)
        g(:,i23)=hcross0*sin(kgr*xx)
      endif
!
    endsubroutine metric_at
!***********************************************************************
    subroutine prim2con(g,gu,e,vu,Bu,w,Econ,Scon)
!
!  Densitized conserved variables sqrt(gamma) E, sqrt(gamma) S_i from the
!  rest-frame energy density e, Eulerian velocity v^i and field B^i,
!  for p = w e and ideal MHD (see hydro_after_boundary_gr).
!
!  13-oct-26/TP: coded
!
      use Sub, only: raise_index,lower_index,contract_vec
      real, dimension(:,:) :: g, gu, vu, Bu, Scon
      real, dimension(:) :: e, Econ
      real :: w
!
      real, dimension(size(e),3) :: vd, Bd, Ed, Eu
      real, dimension(size(e)) :: sg, v2, zz, B2, E2
!
      sg=sqrt(det_sym(g))
      call raise_index(g,vu,vd)
      call lower_index(g,Bu,Bd)
      v2=sum(vu*vd,2); B2=sum(Bu*Bd,2)
      zz=(1+w)*e/(1-v2)
      Ed(:,1)=-sg*(vu(:,2)*Bu(:,3)-vu(:,3)*Bu(:,2))
      Ed(:,2)=-sg*(vu(:,3)*Bu(:,1)-vu(:,1)*Bu(:,3))
      Ed(:,3)=-sg*(vu(:,1)*Bu(:,2)-vu(:,2)*Bu(:,1))
      call lower_index(gu,Ed,Eu)
      E2=contract_vec(Eu,Ed)
      Econ=sg*(zz-w*e+.5*(E2+B2))
      Scon(:,1)=sg*(zz*vd(:,1)+sg*(Eu(:,2)*Bu(:,3)-Eu(:,3)*Bu(:,2)))
      Scon(:,2)=sg*(zz*vd(:,2)+sg*(Eu(:,3)*Bu(:,1)-Eu(:,1)*Bu(:,3)))
      Scon(:,3)=sg*(zz*vd(:,3)+sg*(Eu(:,1)*Bu(:,2)-Eu(:,2)*Bu(:,1)))
!
    endsubroutine prim2con
!***********************************************************************
    subroutine alfven_state(e,w,B0,eta,lboost,c,u,etaL,Etot,Sperp)
!
!  Phase speed c of the circularly polarized Alfven wave in the fluid frame,
!  and, if lboost, the boost velocity u that removes the total momentum
!  along x; also the transverse field amplitude etaL, the total energy
!  density Etot and the transverse momentum amplitude Sperp in the boosted frame.
!
!  14-oct-26/TP: coded
!
      real :: e, w, B0, eta, c, u, etaL, Etot, Sperp
      logical :: lboost
!
      real :: zz, T00, T0x, Txx, T0p, Txp, r, gam2
!
      c=alfven_speed(e,w,B0,eta)
      zz=(1+w)*e/(1-(c*eta)**2)
      T00=zz-w*e+.5*B0**2*(1+eta**2+(c*eta)**2)
      T0x=c*eta**2*B0**2
      Txx=w*e+.5*B0**2*(1+eta**2+(c*eta)**2)-B0**2
      T0p=-(zz+B0**2)*c*eta
      Txp=-eta*B0**2
      u=0.
      if (lboost) then
        r=T0x/(T00+Txx)
        u=2*r/(1+sqrt(1-4*r**2))
      endif
      gam2=1/(1-u**2)
      etaL=eta*sqrt(gam2)*(1-u*c)
      Etot=gam2*(T00-2*u*T0x+u**2*Txx)
      Sperp=sqrt(gam2)*(T0p-u*Txp)
!
    endsubroutine alfven_state
!***********************************************************************
    real function alfven_speed(e,w,B0,eta)
!
!  Phase speed of the circularly polarized Alfven wave.
!
      real :: e, w, B0, eta, Q
!
      Q=(1+w)*e+B0**2*(1+eta**2)
      if (eta==0.) then
        alfven_speed=sqrt(B0**2/Q)
      else
        alfven_speed=sqrt((Q-sqrt(Q**2-4*eta**2*B0**4))/(2*eta**2*B0**2))
      endif
!
    endfunction alfven_speed
!***********************************************************************
    function det_sym(a)
!
      real, dimension(:,:) :: a
      real, dimension(size(a,1)) :: det_sym
!
      det_sym=a(:,i11)*(a(:,i22)*a(:,i33)-a(:,i23)**2) &
             -a(:,i12)*(a(:,i12)*a(:,i33)-a(:,i23)*a(:,i13)) &
             +a(:,i13)*(a(:,i12)*a(:,i23)-a(:,i22)*a(:,i13))
!
    endfunction det_sym
!***********************************************************************
    function contract_sym(au,b)
!
      real, dimension(:,:) :: au, b
      real, dimension(size(au,1)) :: contract_sym
!
      contract_sym=au(:,i11)*b(:,i11)+au(:,i22)*b(:,i22)+au(:,i33)*b(:,i33) &
                  +2*(au(:,i12)*b(:,i12)+au(:,i13)*b(:,i13)+au(:,i23)*b(:,i23))
!
    endfunction contract_sym
!***********************************************************************
    subroutine read_initial_condition_pars(iomsg)
!
      use File_io, only: parallel_unit
!
      character(len=*), intent(out) :: iomsg
      integer :: iostat
!
      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat, IOMSG=iomsg)
      if (iostat==0) iomsg=""
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
