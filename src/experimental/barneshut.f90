!
! This module implements gas/dust self-gravity using a Barnes-Hut octree
! algorithm. The main function, do_barneshut, takes the variable 'phi' as input,
! which is the gas/dust density grid from the selfgravity module, and modifies
! 'phi' in place, becoming the gravitational potential.
!
! To use this module:
! 1. In Makefile.local, set POISSON=experimental/barneshut (FOURIER is
!    subsequently not needed)
! 2. Update start.in with &poisson_init_pars/appropriate parameters.
!
! The variable 'themap' contains a complete list of the region/point pairings
! for the BH map for each processor (where 'regions' are being integrated to
! calculate the potential at each 'point'). Each such pairing is given by ten
! integers: the indices of the local point coordinates (i,j,k), the upper/lower
! boundaries of the region in each dimension, and the processor ID of the
! region. The subroutine 'mkmap' is run twice: first, a "dummy" run of the
! subroutine calculates the total number of regions, after which 'themap' is
! allocated. A second pass then populates 'themap' with the needed information.
! The memory required per core for the map for a reasonably-sized run in
! spherical coordinates is of order 10^2 MB. The runtime then consists of
! MPI-based exchange of density fields between processors, then a single loop
! over the point-region pairings in 'themap' takes care of the integration.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lpoisson=.true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Poisson
!
  use Cdata
  use Cparam
  use Messages
!
  implicit none
!
  real :: octree_theta = 0.5
  logical :: lshowtime = .false. ! Print BH integration at each time step for
                                 ! root processor. Debugging/testing only.
  logical :: lsquareregions = .true. ! .true. = make regions as square-shaped as
                                     ! possible. Results in fewer regions.
  logical :: lprecalcdists = .false. ! Pre-calculate point-region distances.
                            ! Not technically correct, but may work as an
                            ! approximation for sufficiently small octree_theta,
                            ! and provides a significant speedup.
  real :: octree_maxdist = huge(1.0)
  real :: octree_smoothdist = 0.15
  real, dimension (nx,ny,nz) :: phi, vols
  integer :: pp, i, j, k, xs, ys, zs, ii, jj, kk
  integer, dimension (0:ncpus-1) :: sx, sy, sz
  real, dimension (nx) :: xc
  real, dimension (ny) :: yc
  real, dimension (nz) :: zc
  real, dimension (nx,ny,nz) :: xmesh, ymesh, zmesh
  real, dimension (nx,0:ncpus-1) :: xrecv
  real, dimension (ny,0:ncpus-1) :: yrecv
  real, dimension (nz,0:ncpus-1) :: zrecv
  integer, allocatable :: themap(:,:), themap_sort(:,:)
  real, allocatable :: regdist1(:)
  real, allocatable :: regsmooth(:)
  integer :: nregions, iregion, irl, iru, jrl, jru, krl, kru, nprecalc=0
  integer, parameter :: nlt = 1e7 ! 3*80MB for 1e7
  real :: lkmin, lkmax, dxlt
  real, dimension(nlt) :: xlt, sinlt, coslt
!
  namelist /poisson_init_pars/ &
      octree_theta, lshowtime, lsquareregions, lprecalcdists, octree_maxdist, &
      octree_smoothdist
!
  namelist /poisson_run_pars/ &
      octree_theta, lshowtime, lsquareregions, lprecalcdists, octree_maxdist, &
      octree_smoothdist
!
  contains
!***********************************************************************
    subroutine inverse_laplacian(phi)
!
      use General, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: phi
!
      intent(inout) :: phi
!
      if (lcylindrical_coords) then 
        if (lroot) print*,'You are using cylindrical coordinates. '//&
             'Use poisson_cyl.f90 instead'
        call fatal_error("barneshut","")
      endif
      if (modulo(log(real(nx))/log(2.0),1.0) .gt. tini .or. &
        modulo(log(real(ny))/log(2.0),1.0) .gt. tini .or. &
        modulo(log(real(nz))/log(2.0),1.0) .gt. tini) then
        if (lroot) print*,'Grid zones per processor in each axis '//&
            'must be a power of 2.'
        call fatal_error("barneshut","")
      endif
!
      call do_barneshut(phi)
      call keep_compiler_quiet(f)
!
    endsubroutine inverse_laplacian
!***********************************************************************
    subroutine initialize_poisson
!
    use Mpicomm
!
    real :: ldx, xw, yw, zw
    integer, dimension (nx) :: xind
    integer, dimension (ny) :: yind
    integer, dimension (nz) :: zind
    integer :: iprecalc
    logical :: lprecalc = .false.
!
    ! Trimming ghost zones
    xc = x(l1:l2)
    yc = y(m1:m2)
    zc = z(n1:n2)
!
    do i=1,nx
      do j=1,ny
        do k=1,nz
          if (lspherical_coords) then
            xmesh(i,j,k) = xc(i)*sin(yc(j))*cos(zc(k))
            ymesh(i,j,k) = xc(i)*sin(yc(j))*sin(zc(k))
            zmesh(i,j,k) = xc(i)*cos(yc(j))
          else
            xmesh(i,j,k) = xc(i)
            ymesh(i,j,k) = yc(j)
            zmesh(i,j,k) = zc(k)
          endif
        enddo
      enddo
    enddo
!
    do pp=0,ncpus-1
      if (pp/=iproc) then
        call mpisendrecv_real(xc,nx,pp,118,xrecv(:,pp),pp,118)
        call mpisendrecv_real(yc,ny,pp,119,yrecv(:,pp),pp,119)
        call mpisendrecv_real(zc,nz,pp,120,zrecv(:,pp),pp,120)
      endif
    enddo
!
    xrecv(:,iproc) = xc
    yrecv(:,iproc) = yc
    zrecv(:,iproc) = zc
!
    ! Calculate volumes for converting density to mass
    if (lspherical_coords) then
      do j=1,ny
        do i=1,nx
          if (grid_func(1) .eq. 'power-law') then
            ! Not 100% sure formula is correct, right idea though:
            ldx = (((xc(i)-xyz0(1))/xyz1(1))**(1.0-coeff_grid(1))) &
              *xyz1(1)+xyz0(1)
          else
            ldx = dx
          endif
          vols(i,j,:) = ldx*dy*dz*sin(yc(j))*xc(i)**2.0
        enddo
      enddo
    else
      vols = dx*dy*dz
    endif
!
    if (lsquareregions) then
      do pp=0,ncpus-1
        if (lspherical_coords) then
          xw = abs(xrecv(nx,pp)-xrecv(1,pp))
          yw = abs(yrecv(ny,pp)-yrecv(1,pp))*xrecv(nx/2,pp)
          zw = abs(zrecv(nz,pp)-zrecv(1,pp))*xrecv(nx/2,pp)*sin(yrecv(ny/2,pp))
        else
          xw = abs(xrecv(nx,pp)-xrecv(1,pp))
          yw = abs(yrecv(ny,pp)-yrecv(1,pp))
          zw = abs(zrecv(nz,pp)-zrecv(1,pp))
        endif
        ! If a region is longer along one axis (or two), it will be split along
        ! that axis (or axes) to make the resulting sub-regions roughly cubical.
        sx(pp) = max(nx/max(roundtwo(xw/yw),roundtwo(xw/zw),1),1)
        sy(pp) = max(ny/max(roundtwo(yw/xw),roundtwo(yw/zw),1),1)
        sz(pp) = max(nz/max(roundtwo(zw/xw),roundtwo(zw/yw),1),1)
      enddo
    endif
!
    ! Lookup table setup:
    if (lspherical_coords) then
      lkmin = -pi-dy-dz
      lkmax = pi+dy+dz
      dxlt = (lkmax-lkmin)/real(nlt)
      do i=1,nlt
        xlt(i) = lkmin+dxlt*real(i-1)
        sinlt(i) = sin(xlt(i))
        coslt(i) = cos(xlt(i))
      enddo
    endif
!
    do i=1,nx
      xind(i) = i
    enddo
    do j=1,ny
      yind(j) = j
    enddo
    do k=1,nz
      zind(k) = k
    enddo
!
    ! First pass only counts regions, second pass actually populates 'themap'
    do iprecalc=1,2
      if (lprecalc) then
        if (lroot) print*,"# regions on proc 0:",nregions
        allocate(themap(10,nregions))
        allocate(themap_sort(10,nregions))
        allocate(regsmooth(nregions))
        regsmooth = 1.0
        allocate(regdist1(nregions))
      endif
      nregions = 0 ! Advanced inside 'mkmap'
      do pp=0,ncpus-1
        if (lsquareregions) then
          do kk=1,nz/sz(pp)
            do jj=1,ny/sy(pp)
              do ii=1,nx/sx(pp)
                do k=1,nz
                  do j=1,ny
                    do i=1,nx
                      irl = (ii-1)*sx(pp)+1 ; iru = ii*sx(pp)
                      jrl = (jj-1)*sy(pp)+1 ; jru = jj*sy(pp)
                      krl = (kk-1)*sz(pp)+1 ; kru = kk*sz(pp)
                      call mkmap((/i,j,k/),(/sx(pp),sy(pp),sz(pp)/), &
                        xrecv(irl:iru,pp),yrecv(jrl:jru,pp),zrecv(krl:kru,pp), &
                        xind(irl:iru),yind(jrl:jru),zind(krl:kru),pp,lprecalc)
                    enddo !i
                  enddo !j
                enddo !k
              enddo !ii
            enddo !jj
          enddo !kk
        else
          do k=1,nz
            do j=1,ny
              do i=1,nx
                call mkmap((/i,j,k/),(/nx,ny,nz/), xrecv(:,pp),yrecv(:,pp), &
                  zrecv(:,pp),xind,yind,zind,pp,lprecalc)
              enddo !i
            enddo !j
          enddo !k
        endif !lsquareregions
      enddo !pp
      lprecalc = .true.
    enddo
!
    if (lprecalcdists) then
      nprecalc = nregions
    else
      ! sorting 'themap' so that all precalculated regions are at the 
      ! beginning of the array, i.e.,
      ! precalc: themap(1:nprecalc,:)
      ! non-precalc: themap(nprecalc+1:nregions)
      do iregion=1,nregions
        if (themap(4,iregion) .eq. themap(5,iregion) .and. &
        themap(6,iregion) .eq. themap(7,iregion) .and. &
        themap(8,iregion) .eq. themap(9,iregion)) then
          nprecalc = nprecalc+1
          themap_sort(:,nprecalc) = themap(:,iregion)
        else
          themap_sort(:,nregions-(iregion-nprecalc)+1) = themap(:,iregion)
        endif
      enddo
      themap = themap_sort
      deallocate(themap_sort)
    endif
!
    endsubroutine initialize_poisson
!***********************************************************************
    subroutine do_barneshut(phi) ! 'phi' is density from selfgravity module,
                                 !  will be transformed into potential.
!
    use Mpicomm
!
    real, dimension (nx,ny,nz,0:ncpus-1) :: phirecv
    real, dimension (nx,ny,nz) :: phi
    real :: xreg, yreg, zreg, summreg, summreg1
    real :: tstart, tstop
!
    if (lshowtime .and. lroot) call cpu_time(tstart)
!
    phirecv = 0.0
    phi = phi*vols ! 'phi' now in mass units
!
    ! Send masses to other processors
    do pp=0,ncpus-1
      if (pp/=iproc) then
        call mpisendrecv_real(phi,(/nx,ny,nz/),pp,117, &
                              phirecv(:,:,:,pp),pp,117)
      endif
    enddo
!
    phirecv(:,:,:,iproc) = phi
!
    phi = 0.0
!
    do iregion=1,nprecalc
      i   = themap(1, iregion) ; irl = themap(4,iregion) ; iru = themap(5,iregion)
      j   = themap(2, iregion) ; jrl = themap(6,iregion) ; jru = themap(7,iregion)
      k   = themap(3, iregion) ; krl = themap(8,iregion) ; kru = themap(9,iregion)
      pp  = themap(10,iregion)
      phi(i,j,k) = phi(i,j,k) - regsmooth(iregion)* &
        sum(phirecv(irl:iru,jrl:jru,krl:kru,pp))*regdist1(iregion)
    enddo
    do iregion=nprecalc+1,nregions
      i   = themap(1, iregion) ; irl = themap(4,iregion) ; iru = themap(5,iregion)
      j   = themap(2, iregion) ; jrl = themap(6,iregion) ; jru = themap(7,iregion)
      k   = themap(3, iregion) ; krl = themap(8,iregion) ; kru = themap(9,iregion)
      pp  = themap(10,iregion)
      summreg = sum(phirecv(irl:iru,jrl:jru,krl:kru,pp))*regsmooth(iregion)
      summreg1 = 1.0/summreg
      xreg = sum(xrecv(irl:iru,pp) &
        *sum(sum(phirecv(irl:iru,jrl:jru,krl:kru,pp),3),2))*summreg1
      yreg = sum(yrecv(jrl:jru,pp) &
        *sum(sum(phirecv(irl:iru,jrl:jru,krl:kru,pp),3),1))*summreg1
      zreg = sum(zrecv(krl:kru,pp) &
        *sum(sum(phirecv(irl:iru,jrl:jru,krl:kru,pp),2),1))*summreg1
      phi(i,j,k) = phi(i,j,k) - summreg/get_dist((/i,j,k/),(/xreg,yreg,zreg/))
    enddo
!
    phi = phi/(4.0*pi) ! The selfgravity module is going to multiply by a factor
                       ! of 4pi that we don't want (I think).
!
    if (lshowtime .and. lroot) then
      call cpu_time(tstop)
      print '("Octree integration time = ",f6.3," seconds.")',tstop-tstart
    endif
!
    endsubroutine do_barneshut
!***********************************************************************
    recursive subroutine mkmap(ipos,dimi,xi,yi,zi,xsind,ysind,zsind,ppi,laddreg)
      integer, dimension(3) :: ipos, dimi ! Local position; dimensions of region
      real, dimension(dimi(1)) :: xi ! Region coordinates
      real, dimension(dimi(2)) :: yi !
      real, dimension(dimi(3)) :: zi !
      integer, dimension(dimi(1)) :: xsind ! Region indices (for building map)
      integer, dimension(dimi(2)) :: ysind !
      integer, dimension(dimi(3)) :: zsind !
      integer :: sx, sy, sz, ii, ji, ki, il, iu, jl, ju, kl, ku, ppi
      real :: lxi, lyi, lzi, dist, xwi, ywi, zwi, width
      logical :: laddreg, lsplit
!
      ! If the point where the potential is being calculated is inside the
      ! region in question, then we MUST subdivide further (unless the region is
      ! a single cell-- this is taken care of in the 'if' statement below).
      lsplit = ((ipos(1) .ge. xsind(1) .and. ipos(1) .le. xsind(dimi(1))) &
          .and. (ipos(2) .ge. ysind(1) .and. ipos(2) .le. ysind(dimi(2))) &
          .and. (ipos(3) .ge. zsind(1) .and. ipos(3) .le. zsind(dimi(3))) &
          .and. (iproc .eq. ppi))
!
      lxi = sum(xi)/real(dimi(1))
      lyi = sum(yi)/real(dimi(2))
      lzi = sum(zi)/real(dimi(3))
      dist = get_dist(ipos,(/lxi,lyi,lzi/))
      xwi = abs(xi(dimi(1))-xi(1))
      ywi = abs(yi(dimi(2))-yi(1))
      zwi = abs(zi(dimi(3))-zi(1))
      if (lspherical_coords) then
        width = max(xwi,lxi*ywi,lxi*zwi*sin(lyi))
      else
        width = max(xwi,ywi,zwi)
      endif
      if ((width/dist > octree_theta .or. lsplit) .and. product(dimi)/=1) then
        sx = max(dimi(1)/2,1)
        sy = max(dimi(2)/2,1)
        sz = max(dimi(3)/2,1)
        do ki=1,dimi(3)/sz
          do ji=1,dimi(2)/sy
            do ii=1,dimi(1)/sx
              il = sx*(ii-1)+1; iu = sx*ii
              jl = sy*(ji-1)+1; ju = sy*ji
              kl = sz*(ki-1)+1; ku = sz*ki
              call mkmap(ipos,(/sx,sy,sz/),xi(il:iu),yi(jl:ju),zi(kl:ku), &
                xsind(il:iu),ysind(jl:ju),zsind(kl:ku),ppi,laddreg)
            enddo
          enddo
        enddo
      elseif (.not. lsplit .and. dist .lt. octree_maxdist) then
        nregions = nregions+1
        if (laddreg) then
          themap(:,nregions) = (/ipos(1),ipos(2),ipos(3), xsind(1),xsind(dimi(1)), &
            ysind(1),ysind(dimi(2)),zsind(1),zsind(dimi(3)),ppi/)
          if (lprecalcdists) regdist1(nregions) = 1.0/dist
          if (dist .gt. (octree_maxdist-octree_smoothdist)) then
           regsmooth(nregions) = 0.5* &
            (cos(pi*((octree_maxdist-dist)/octree_smoothdist+1.0))+1.0)
          endif
        endif
      endif
    endsubroutine mkmap
!***********************************************************************
    function get_dist(p1,p2) result(rdist)
    integer, dimension(3) :: p1
    real, dimension(3) :: p2
    real :: x2, y2, z2, rdist
    real, dimension(2) :: sincosy, sincosz
!
    if (lspherical_coords) then
      sincosy = sincoslf(p2(2))
      sincosz = sincoslf(p2(3))
      x2 = p2(1)*sincosy(1)*sincosz(2)
      y2 = p2(1)*sincosy(1)*sincosz(1)
      z2 = p2(1)*sincosy(2)
    else
      x2 = p2(1)
      y2 = p2(2)
      z2 = p2(3)
    endif
      rdist = sqrt((xmesh(p1(1),p1(2),p1(3))-x2)**2+ &
        (ymesh(p1(1),p1(2),p1(3))-y2)**2+(zmesh(p1(1),p1(2),p1(3))-z2)**2)
!
    endfunction get_dist
!***********************************************************************
    function sincoslf(ang)
    real,dimension(2) :: sincoslf
    real :: ang, sinx, cosx !, a0, a1, a2
    integer :: ilk
!
    ilk = nint((ang-lkmin)/dxlt)+1
    sinx = sinlt(ilk)
    cosx = coslt(ilk)
    ! Quadratic interpolation (probably not needed; slower but uses less memory)
    !a0 = (ang-xlt(ilk))*(ang-xlt(ilk+1))/ &
    !  ((xlt(ilk-1)-xlt(ilk))*(xlt(ilk-1)-xlt(ilk+1)))
    !a1 = (ang-xlt(ilk-1))*(ang-xlt(ilk+1))/ &
    !  ((xlt(ilk)-xlt(ilk-1))*(xlt(ilk)-xlt(ilk+1)))
    !a2 = (ang-xlt(ilk-1))*(ang-xlt(ilk))/ &
    !  ((xlt(ilk+1)-xlt(ilk-1))*(xlt(ilk+1)-xlt(ilk)))
    !sinx = sinlt(ilk-1)*a0+sinlt(ilk)*a1+sinlt(ilk+1)*a2
    !cosx = coslt(ilk-1)*a0+coslt(ilk)*a1+coslt(ilk+1)*a2
    sincoslf = (/sinx,cosx/)
!
    endfunction sincoslf
!***********************************************************************
    function roundtwo(rin) result(iout)
    real :: rin
    integer :: iout
!
    iout = 2**nint(log(rin)/log(2.0))
    endfunction roundtwo
!***********************************************************************
    subroutine read_poisson_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=poisson_init_pars, IOSTAT=iostat)
!
    endsubroutine read_poisson_init_pars
!***********************************************************************
    subroutine write_poisson_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=poisson_init_pars)
!
    endsubroutine write_poisson_init_pars
!***********************************************************************
    subroutine read_poisson_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=poisson_run_pars, IOSTAT=iostat)
!
    endsubroutine read_poisson_run_pars
!***********************************************************************
    subroutine write_poisson_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=poisson_run_pars)
!
    endsubroutine write_poisson_run_pars
!***********************************************************************
    subroutine inverse_laplacian_semispectral(f,phi)
!
!  Solve the Poisson equation by Fourier transforming on a periodic grid.
!
!  15-may-2006/anders+jeff: dummy
!
      use General, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: phi
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(phi)
      call keep_compiler_quiet(f)
!
    endsubroutine inverse_laplacian_semispectral
!***********************************************************************
endmodule Poisson
