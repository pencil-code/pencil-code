! $Id: interstellar.f90,v 1.2 2002-11-19 14:58:58 ngrs Exp $

!  This modules solves contains ISM and SNe 

module Interstellar

  use Cparam
  use Cdata
  use Density

  implicit none

  real :: t_next_SNI=0e0,t_interval_SNI=3.64e-3,h_SNI=0.325,ampl_SN=5.0
  real :: x_SN,y_SN,z_SN

  ! input parameters
  integer :: dummy
  namelist /interstellar_init_pars/ dummy

  ! run parameters
  namelist /interstellar_run_pars/ &
      t_next_SNI,t_interval_SNI,h_SNI,ampl_SN,x_SN,y_SN,z_SN

 
  contains

!***********************************************************************
    subroutine register_interstellar()
!
!  19-nov-02/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_interstellar called twice')
      first = .false.
!
      linterstellar = .true.
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_lncc'
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: interstellar.f90,v 1.2 2002-11-19 14:58:58 ngrs Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_interstellar: nvar > mvar')
      endif
!
    endsubroutine register_interstellar
!***********************************************************************
    subroutine calc_heat_cool_interstellar(df,rho1,TT1)
!
!  adapted from calc_heat_cool
!
      use Cdata
      use Mpicomm
      use Density, only : rho0
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df  !,f
      real, dimension (nx) :: rho1,TT1       !,cs2
      real, dimension (nx) :: rho,TT,heat,cool
!  cp1=1/cp used to convert TT (and ss) into interstellar code units
!  (useful, as many conditions conveniently expressed in terms of TT)
!  code units based on:
!    [length]  = 1kpc  = 3.09 10^21 cm
!    [time]    = 1Gyr  = 3.15 10^16 s
!    [rho]     =       = 1.00 10^-24 g/cm^3
!  Lambdaunits converts coolH into interstellar code units.
!   (this should really just be incorporated into coolH coefficients)
      real :: cp1=27.8,Lambdaunits=3.29e-18
      real :: rhoUV=0.1,TUV=7000.,T0UV=12000.,cUV=5.e-4
      real, dimension(6) ::                                                 &
        coolT=(/ 500.,     2000.,    8000.,    1.e5,    4.e7,     1.e9 /),  &
        coolH=(/ 5.595e15, 2.503e17, 1.156e12, 4.45e29, 8.054e20, 0.   /),  &
        coolB=(/ 2.,       1.5,      2.867,    -.65,    0.5,      0.   /)
      integer :: i,l
!
      intent(in) :: rho1,TT1   !,f,cs2
      intent(inout) :: df
!
!  identifier
!
      if(headtt) print*,'calc_heat_cool_interstellar'
!
!  define T in K, for calculation of both UV heating and radiative cooling
!
      TT=cp1/TT1 
      rho=1.0/rho1
!
!  add T-dept radiative cooling, from Rosen et al., ApJ, 413, 137, 1993
!  cooling is Lambda*rho^2, with (eq 7)
!     Lambda=coolH(i)*TT*coolB(i),   for coolT(i) <= T < coolT(i+1)
!  nb: our coefficients coolH(i) differ from those in Rosen et al. by
!   factor (mu mp)^2, with mu=1.2, since Rosen works in number density, n.
!   (their cooling = Lambda*n^2,  rho=mu mp n.)
!  The factor Lambdaunits converts from cgs units to code units.
!  We've increased coolT(1), to avoid creating gas too cold to resolve.
!
      cool=0.0
      do l=l1,l2
        do i=1,5
          if (coolT(i) <= TT(l) .and. TT(l) < coolT(i+1)) then
            cool(l)=cool(l) + Lambdaunits*rho(l)*coolH(i)*TT(l)**coolB(i)
          endif
        enddo
      enddo
!
!  add UV heating, cf. Wolfire et al., ApJ, 443, 152, 1995
!  with the values above, this gives about 0.012 erg/g/s (T < ~1.e4 K)
!  nb: need rho0 from density_[init/run]_pars, if i want to implement
!      the the arm/interarm scaling.
!
      heat=rhoUV*(rho0/1.38)**1.4*Lambdaunits*coolH(3)*TUV**coolB(3)*   &
                               0.5*(1.0+tanh(cUV*(T0UV-TT)))
!      heat=rhoUV*Lambdaunits*coolH(3)*TUV**coolB(3)*         &
!                               0.5*(1.0+tanh(cUV*(T0UV-TT)))
!
!  heat and cool were calculated in terms of de/dt [~ erg/g/s], 
!  so just multiply by TT1 to get ds/dt [~ erg/g/s/K]:
!
      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + TT1*(heat - cool)
!
    endsubroutine calc_heat_cool_interstellar
!***********************************************************************
    subroutine check_SN(f,df)
!
!  Checks for SNe, and implements appropriately:
!   relevant subroutines in entropy.f90
!
    use Cdata
!
    real, dimension(mx,my,mz,mvar) :: f,df
    logical :: l_SNI=.false.   !only allow SNII if no SNI this step
!
!  Do separately for SNI (simple scheme) and SNII (Boris' scheme)
!
    call check_SNI(f,df,l_SNI)
    call check_SNII(f,df,l_SNI)
!
    endsubroutine check_SN
!***********************************************************************
    subroutine check_SNI(f,df,l_SNI)
!
!  If time for next SNI, then implement, and calculate time of subsequent SNI
!
    use Cdata
    use Mpicomm
    use General
!
    real, dimension(mx,my,mz,mvar) :: f,df
    real, dimension(1) :: fran
    real :: rnd
    logical :: l_SNI
!
    l_SNI=.false.
    if (t >= t_next_SNI) then
      call position_SNI
      !call explode_SN(f,df,1)
!  pre-determine time for next SNI
      if (lroot) then
        call random_number_wrapper(fran)   
        rnd=(fran(1)-0.5)*0.4
        t_next_SNI=t_next_SNI + rnd*t_interval_SNI
        print*,'Next SNI at time: ',t_next_SNI
      endif
      call mpibcast_real(t_next_SNI,1)
      l_SNI=.true.
    endif
!
    endsubroutine check_SNI
!***********************************************************************
    subroutine check_SNII(f,df,l_SNI)
!
!   check if time for next SNII.  based on current density/T structure...
!
    use Cdata
    use General
! 
    real, dimension(mx,my,mz,mvar) :: f,df
    logical :: l_SNI
!
    if (.not. l_SNI) then
      if (ip==0) print*, 'check_SNII: still to write'
    endif
!
    endsubroutine check_SNII
!***********************************************************************
    subroutine position_SNI()
!
!   determine position for next SNI (w/ fixed scale-height)
!
    use Cdata
    use Mpicomm
    use General
!
    real, dimension(nzgrid) :: cum_prob_SNI
    real :: rnd,zn
    real, dimension(3) :: fran,fmax
    integer :: i, nzskip=20   !prevent SNI from being too close to boundaries
!
!  Lx,Ly,Lz,x0,y0,z0 not necessarily (correctly) set in cdata module. 
!  (dx,dy,dz OK.  Lx,Ly,Lz assume periodic domains, not true for z.
!   x0,y0,z0 not set (their value also depends on periodicity.) )
!  The following assumes x, y periodic, z non-periodic (see start.f90)
    Lx=Lxyz(1);       Ly=Lxyz(2);       Lz=Lxyz(3)
    !dx=Lx/nxgrid;     dy=Ly/nygrid;     dz=Lz/(nzgrid-1)    !already OK
    x0=xyz0(1)+.5*dx; y0=xyz0(2)+.5*dy; z0=xyz0(3)
    !if (lroot) print*, 'dx,dy,dz',dx,dy,dz
    !if (lroot) print*, 'x0,y0,z0,Lx,Ly,Lz',x0,y0,z0,Lx,Ly,Lz
!
!  Cumulative probability funcion in z currently calculated each time.
!  It constant, and could be stored (and calculated in init)
!
    cum_prob_SNI(1:nzskip)=0.0
    do n=nzskip+1,nzgrid-nzskip
      zn=z0+(n-1)*dz
      cum_prob_SNI(n)=cum_prob_SNI(n-1)+exp(-(zn/h_SNI)**2)
    enddo
    cum_prob_SNI=cum_prob_SNI/cum_prob_SNI(nzgrid-nzskip)
!  The following should never be needed, but just in case floating point 
!  errors ever lead to cum_prob_SNI(nzgrid-nzskip) < rnd < 1.
    cum_prob_SNI(nzgrid-nzskip+1:nzgrid)=1.0   
    !if (lroot) print*, 'cum_prob',cum_prob_SNI
!
    call random_number(fran)
    rnd=fran(1) ! call random_number(rnd)   
    do n=nzskip+1,nzgrid-nzskip
      if (cum_prob_SNI(n-1) <= rnd .and. rnd < cum_prob_SNI(n)) &
        z_SN=z0+(n-1)*dz
    enddo
    !if (lroot) print*, 'z',rnd,z_SN
!
    rnd=fran(2) ! call random_number(rnd)   
    i=int(rnd*nxgrid)+1
    x_SN=x0+(i-1)*dx
    !if (lroot) print*, 'x',rnd,x_SN
!
    rnd=fran(3) ! call random_number(rnd)   
    i=int(rnd*nygrid)+1
    y_SN=y0+(i-1)*dy
    !if (lroot) print*, 'y',rnd,y_SN
!
    fmax=(/ x_SN, y_SN, z_SN /)
    call mpibcast_real(fmax,3)
    x_SN=fmax(1); y_SN=fmax(2); z_SN=fmax(3)
    if (lroot) print*, 'position_SNI:',x_SN,y_SN,z_SN
!
    endsubroutine position_SNI

!***********************************************************************
    subroutine position_SNII(f,df)
!
!   determine position for next SNII (using Boris' scheme)
!
    use Cdata
    use General
!
    real, dimension(mx,my,mz,mvar) :: f,df
!
    if (ip==0) print*, 'check_SNII: still to write'
!
    endsubroutine position_SNII
!***********************************************************************
    subroutine explode_SN(f,df,i_SN)
!
!   implement SN (of either type), at pre-calculated position
!
    use Cdata
    use Mpicomm
!
    real, dimension(mx,my,mz,mvar) :: f,df
    real :: dx_SN_in,dx_SN_out_x0,dx_SN_out_x1,dy_SN_in,dy_SN_out_y
    real :: dy_SN_out_x0a,dy_SN_out_x0b,dy_SN_out_x1a,dy_SN_out_x1b
    real :: dz_SN,yshift
    real, dimension(nx,ny,nz) :: dr2_SN
    real :: width_SN,c_SN,rho_SN,TT_SN,lnrhol,ssl
!  The following arrays can be streamlined, after debugging.
    real, dimension(nx) :: lnrho,rho,dd,TT,ss,dss,profile
    real :: cnorm_SN=1.5484             ! (int exp(-r^6) 4\pi r^2 dr)^{1/3}
    real :: cp1=27.8,TT_limit=1.e7       
    real, dimension(2) :: fmax,fmax_tmp 
    integer :: i_SN,l,mshift
    integer :: point_width=4
    logical :: lfound
!
!  Need rho and TT at the SN location.  At present done very basicly.
!   (Be careful that if-condition is actually satisfied, in f.p...)
!  Could be done in position_SN, to save on communication -- but then
!   need access to f in position_SN.  (And also need do similar loop
!   conditional anyway -- or else be able to broadcast from non-root.)
! 
    lfound=.false.
    rho_SN=0.0; TT_SN=0.0
    do n=n1,n2
      do m=m1,m2
        do l=l1,l2
          if (x(l)==x_SN .and. y(m)==y_SN .and. z(n)==z_SN) then
            lnrhol=f(l,m,n,ilnrho)
            ssl=f(l,m,n,ient)
            rho_SN=exp(lnrhol)
! TT_SN=cs2_SN/gamma1;  *cp1 to give in K
            TT_SN=cs20*exp(gamma1*(lnrhol-lnrho0) + gamma*ssl)/gamma1*cp1
            lfound=.true.
          endif
        enddo 
      enddo 
    enddo
    if (.not. lfound) print*,'explode_SN: SN not found!'
!
!  Use reduce_max to get the true SN values on all processors
!
    fmax_tmp=(/ rho_SN, TT_SN /)
    call mpireduce_max(fmax_tmp,fmax,2)
    call mpibcast_real(fmax,2)
    rho_SN=fmax(1); TT_SN=fmax(2)
    if (lroot) print*,'explode_SN: pre-SN, rho_SN, TT_SN:',rho_SN,TT_SN
!
!  Obtain distance to SN.  If incorporated with explosion itself,
!   this could be done on pencils to save memory.
!   (w/ real, dimension(nx) :: dr2_SN,dx_SN_in,dx_SN_out_x0,dx_SN_out_x1)
!  Make sure this will work w/ mass relocation, though.
!
    do n=n1,n2
      dz_SN=abs(z(n)-z_SN)
      do m=m1,m2
!  consider all possible positions in xy plane, to get shortest.
!  can be made more efficient later.
!  y-separations with no boundaries crossed in x
        dy_SN_in=abs(y(m)-y_SN)                 !dyi
        dy_SN_out_y=Ly-dy_SN_in                 !dyii
!  y-separations across sliding periodic boundary at x=x0
        yshift=y(m)-deltay
        mshift=0
        if (yshift < y0) then 
          mshift=int((y0-yshift)/Ly)+1
          yshift=yshift+mshift*Ly
        elseif (yshift > y0+Ly) then      
          mshift=int((yshift-(y0+Ly))/Ly)+1
          yshift=yshift-mshift*Ly
        endif
        dy_SN_out_x0a=abs(yshift-y_SN)          !dyiii
        dy_SN_out_x0b=Ly-dy_SN_out_x0a          !dyiv
!  y-separations across sliding periodic boundary at x=x1
        yshift=y(m)+deltay
        mshift=0
        if (yshift < y0) then 
          mshift=int((y0-yshift)/Ly)+1
          yshift=yshift+mshift*Ly
        elseif (yshift > y0+Ly) then
          mshift=int((yshift-(y0+Ly))/Ly)+1
          yshift=yshift-mshift*Ly
        endif
        dy_SN_out_x1a=abs(yshift-y_SN)          !dyv
        dy_SN_out_x1b=Ly-dy_SN_out_x1a          !dyvi
        do l=l1,l2
!  x-separations associated with each of the above
          dx_SN_in=abs(x(l)-x_SN)               !dxi=dxii
          dx_SN_out_x0=Lx+(x(l)-x_SN)           !dxiii=dxiv
          dx_SN_out_x1=Lx-(x(l)-x_SN)           !dxv=dxvi
          dr2_SN(l,m,n)=min(                                                 &
       dx_SN_in**2+dy_SN_in**2, dx_SN_in**2+dy_SN_out_y**2,                  &
       dx_SN_out_x0**2+dy_SN_out_x0a**2, dx_SN_out_x0**2+dy_SN_out_x0b**2,   &
       dx_SN_out_x1**2+dy_SN_out_x1a**2, dx_SN_out_x1**2+dy_SN_out_x1b**2)   &
                  + dz**2
        enddo
      enddo
    enddo
!
!  Can now add energy, and, if necessary, relocate mass
!
    width_SN=point_width*dx
    c_SN=ampl_SN/(cnorm_SN*width_SN)**3  !normalision for profile below
!
!  If nec., create shell in density, as well as adding energy
!!    if (TT_SN < TT_limit) then
!!      print*, 'explode_SN, mass relocation: still to write'
!  No need to move mass;  simply inject energy
!!    else     ! if TT_SN >= TT_limit
!  The following is v. inefficient -- but should help with debugging
!  Put in summation checks to ensure energy is added properly..
      do n=n1,n2
        do m=m1,m2
          profile=exp(-min( (dr2_SN(:,m,n)/width_SN**2)**3, 75. ))
          lnrho=f(l1:l2,m,n,ilnrho)
          rho=exp(lnrho)
          ss=f(l1:l2,m,n,ient)
          TT=cs20*exp(gamma1*(lnrho-lnrho0) + gamma*ss)/gamma1*cp1
          dss=c_SN*profile/rho/TT           ! s in dimensional units
          df(l1:l2,m,n,ient)=dss*cp1        ! back to s/cp
        enddo
      enddo
!!    endif 
    
    endsubroutine explode_SN
!***********************************************************************
endmodule interstellar
