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
  logical :: lsquareregions = .false. ! .true. = make regions as square-shaped as
                                      ! possible. Results in fewer regions.
  logical :: lprecalcdists = .false. ! Pre-calculate point-region distances.
                            ! Not technically correct, but may work as an
                            ! approximation for sufficiently small octree_theta,
                            ! and provides a significant speedup.
  logical :: lnorepeatsumming = .true. ! (very) slow to start, quick to run
  logical :: lmakecartoon = .false. ! print some map data to stdout for tree illustration purposes
  logical :: ltreestatus = .false.
  logical :: lreadoctree = .false.
  logical :: lwriteoctree = .false.
  real :: octree_maxdist = 1e308
  real :: octree_smoothdist = 0.15
  real, dimension (nx,ny,nz) :: vols
  integer, dimension (0:ncpus-1) :: sx, sy, sz
  real, dimension (nx) :: xc
  real, dimension (ny) :: yc
  real, dimension (nz) :: zc
  real, dimension (nx,ny,nz) :: xmesh, ymesh, zmesh
  real, dimension (nx,0:ncpus-1) :: xrecv
  real, dimension (ny,0:ncpus-1) :: yrecv
  real, dimension (nz,0:ncpus-1) :: zrecv
  integer, allocatable :: themap_group(:,:), themap_single(:,:)
  real, allocatable :: regdist1_single(:), regsmooth_single(:)
  real, allocatable :: regdist1_group(:), regsmooth_group(:)
  logical, allocatable :: luseprevioussum(:)
  integer :: nsingle, ngroup, iregion, irl, iru, jrl, jru, krl, kru, nprecalc=0
  integer, parameter :: nlt = 1e7 ! 3*80MB for 1e7
  real :: lkmin, lkmax, dxlt
  real, dimension(nlt) :: xlt, sinlt, coslt
!
  namelist /poisson_init_pars/ &
      octree_theta, lshowtime, lsquareregions, lprecalcdists, octree_maxdist, &
      octree_smoothdist, lnorepeatsumming, ltreestatus, lwriteoctree, lreadoctree, &
      lmakecartoon
!
  namelist /poisson_run_pars/ &
      octree_theta, lshowtime, lsquareregions, lprecalcdists, octree_maxdist, &
      octree_smoothdist, lnorepeatsumming, ltreestatus, lwriteoctree, lreadoctree
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
    integer :: pp, i, j, k, xs, ys, zs, ii, jj, kk
    intent(inout) :: phi
!
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
    real :: xw, yw, zw, dx1, dy1, dz1
    integer, dimension (nx) :: xind
    integer, dimension (ny) :: yind
    integer, dimension (nz) :: zind
    real, dimension (nx) :: dxc_1
    real, dimension (ny) :: dyc_1
    real, dimension (nz) :: dzc_1
    integer :: pp, i, j, k, xs, ys, zs, ii, jj, kk
    integer :: iprecalc, iregtmp, info, itmp
    logical :: lprecalc = .false., cartooncond, doingsort
    logical, allocatable :: tmpmask(:)
    character :: treestatus
!
    ! Trimming ghost zones
    xc = x(l1:l2)
    yc = y(m1:m2)
    zc = z(n1:n2)
    dxc_1 = dx_1(l1:l2)
    dyc_1 = dy_1(m1:m2)
    dzc_1 = dz_1(n1:n2)
!
    do i=1,nx
      do j=1,ny
        do k=1,nz
          if (lspherical_coords) then
            xmesh(i,j,k) = xc(i)*sin(yc(j))*cos(zc(k))
            ymesh(i,j,k) = xc(i)*sin(yc(j))*sin(zc(k))
            zmesh(i,j,k) = xc(i)*cos(yc(j))
          elseif (lcylindrical_coords) then
            xmesh(i,j,k) = xc(i)*cos(yc(j))
            ymesh(i,j,k) = xc(i)*sin(yc(j))
            zmesh(i,j,k) = zc(k)
          elseif (lcartesian_coords) then
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
        call mpisendrecv_real(xc,nx,pp,281,xrecv(:,pp),pp,281)
        call mpisendrecv_real(yc,ny,pp,282,yrecv(:,pp),pp,282)
        call mpisendrecv_real(zc,nz,pp,283,zrecv(:,pp),pp,283)
      endif
    enddo
!
    xrecv(:,iproc) = xc
    yrecv(:,iproc) = yc
    zrecv(:,iproc) = zc
!
    ! Calculate volumes for converting density to mass
    do i=1,nx
      do j=1,ny
        do k=1,nz
          if (nxgrid ==1) then
            dx1=1.0!/Lxyz(1)
          else
            dx1=dxc_1(i)
          endif
          if (nygrid ==1) then
            dy1=1.0!/Lxyz(2)
          else
            dy1=dyc_1(j)
          endif
          if (nzgrid ==1) then
            dz1=1.0!/Lxyz(3)
          else
            dz1=dzc_1(k)
          endif
          vols(i,j,k) = 1.0/(dx1*dy1*dz1)
          if (lcylindrical_coords .and. nygrid/=1) then
            vols(i,j,k) = vols(i,j,k)*xc(i)
          elseif (lspherical_coords) then
            if (nygrid/=1) vols(i,j,k) = vols(i,j,k)*xc(i)
            if (nzgrid/=1) vols(i,j,k) = vols(i,j,k)*sin(yc(j))*xc(i)
          endif
        enddo
      enddo
    enddo
!
    if (lsquareregions) then
      do pp=0,ncpus-1
        if (lspherical_coords) then
          xw = abs(xrecv(nx,pp)-xrecv(1,pp))
          yw = abs(yrecv(ny,pp)-yrecv(1,pp))*xrecv(nx/2,pp)
          zw = abs(zrecv(nz,pp)-zrecv(1,pp))*xrecv(nx/2,pp)*sin(yrecv(ny/2,pp)) ! this needs to get fixed for 2D case (where ny/2=0)
        elseif (lcylindrical_coords) then
          xw = abs(xrecv(nx,pp)-xrecv(1,pp))
          yw = abs(yrecv(ny,pp)-yrecv(1,pp))*xrecv(nx/2,pp)
          zw = abs(zrecv(nz,pp)-zrecv(1,pp))
        elseif (lcartesian_coords) then
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
    if (lspherical_coords .or. lcylindrical_coords) then
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
!        if (lroot) print*,"barneshut: 1x1x1 ~regions/proc:",nsingle,"; >1x1x1 ~regions/proc:",ngroup
        if (lroot) print '("barneshut: ",I0," 1x1x1 and ",I0," >1x1x1 region pairings on proc 0")', nsingle, ngroup
        allocate(themap_group(10,ngroup))
        allocate(themap_single(10,nsingle))
        allocate(regsmooth_single(nsingle))
        allocate(regsmooth_group(ngroup))
        regsmooth_single = 1.0 ; regsmooth_group = 1.0
        allocate(regdist1_single(nsingle))
        if (lprecalcdists) allocate(regdist1_group(ngroup))
        if (lreadoctree) then
          call read_octree
          exit
        endif
      endif
      nsingle = 0 ! Advanced inside 'mkmap'
      ngroup = 0 ! Advanced inside 'mkmap'
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
                if (ltreestatus.and.lprecalc) then
                  write(*,char(nsingle+ngroup))
                endif
              enddo !i
            enddo !j
          enddo !k
        endif !lsquareregions
      enddo !pp
      lprecalc = .true.
    enddo !iprecalc
!
    allocate(luseprevioussum(ngroup))
    luseprevioussum = .false.
    if (lnorepeatsumming) then
      do iregion=2,ngroup
        if (all(themap_group(4:,iregion).eq.themap_group(4:,iregion-1))) &
          luseprevioussum(iregion) = .true.
      enddo
    endif
!
    if (lwriteoctree) call write_octree
!
!    if (lmakecartoon) then
!      do iregion=1,nregions
!        if (lcylindrical_coords.or.lcartesian_coords) then
!          cartooncond = (themap(2,iregion).eq.(ny/2))
!        else
!          cartooncond = (themap(3,iregion).eq.(nz/2))
!        endif
!        if ((themap(1,iregion).eq.(nx/2)).and.cartooncond) then
!          print*,themap(:,iregion)
!        endif
!      enddo
!    endif
!
    endsubroutine initialize_poisson
!***********************************************************************
    subroutine do_barneshut(phi) ! 'phi' is density from selfgravity module,
                                 !  will be transformed into potential.
!
    use Mpicomm
!
    real, dimension (nx,ny,nz,0:ncpus-1) :: phirecv=0.0
    real, dimension (nx,ny,nz) :: phi
    real :: xreg, yreg, zreg, summreg, summreg1
    real :: tstart_loop1, tstop_loop1, tstart_loop2, tstop_loop2, tstart_mpi, tstop_mpi
    integer :: pp, i, j, k, xs, ys, zs, ii, jj, kk
!
    phi = phi*vols ! 'phi' now in mass units
!
    ! Send masses to other processors
    if (lroot .and. lshowtime) call cpu_time(tstart_mpi)
    do pp=0,ncpus-1
      if (pp/=iproc) then
        call mpisendrecv_real(phi,(/nx,ny,nz/),pp,284,phirecv(:,:,:,pp),pp,284)
      endif
    enddo
!
    if (lroot .and. lshowtime) call cpu_time(tstop_mpi)
!
    phirecv(:,:,:,iproc) = phi
!
    phi = 0.0
!
    if (lroot .and. lshowtime) call cpu_time(tstart_loop1)
    do iregion=1,nsingle
      i   = themap_single(1, iregion) ; irl = themap_single(5,iregion) ; iru = themap_single(6,iregion)
      j   = themap_single(2, iregion) ; jrl = themap_single(7,iregion) ; jru = themap_single(8,iregion)
      k   = themap_single(3, iregion) ; krl = themap_single(9,iregion) ; kru = themap_single(10,iregion)
      pp  = themap_single(4, iregion)
      phi(i,j,k) = phi(i,j,k) - regsmooth_single(iregion)* &
        sum(phirecv(irl:iru,jrl:jru,krl:kru,pp))*regdist1_single(iregion)
    enddo
    if (lroot .and. lshowtime) call cpu_time(tstop_loop1)
    if (lroot .and. lshowtime) call cpu_time(tstart_loop2)
    if (.not.lprecalcdists) then
      do iregion=1,ngroup
        i   = themap_group(1,iregion) 
        j   = themap_group(2,iregion) 
        k   = themap_group(3,iregion) 
        if (.not.luseprevioussum(iregion)) then
          pp  = themap_group(4,iregion)
          irl = themap_group(5,iregion) ; iru = themap_group(6,iregion)
          jrl = themap_group(7,iregion) ; jru = themap_group(8,iregion)
          krl = themap_group(9,iregion) ; kru = themap_group(10,iregion)
          summreg = sum(phirecv(irl:iru,jrl:jru,krl:kru,pp))
          summreg1 = 1.0/summreg
          xreg = sum(xrecv(irl:iru,pp) &
            *sum(sum(phirecv(irl:iru,jrl:jru,krl:kru,pp),3),2))*summreg1
          yreg = sum(yrecv(jrl:jru,pp) &
            *sum(sum(phirecv(irl:iru,jrl:jru,krl:kru,pp),3),1))*summreg1
          zreg = sum(zrecv(krl:kru,pp) &
            *sum(sum(phirecv(irl:iru,jrl:jru,krl:kru,pp),2),1))*summreg1
        endif
        phi(i,j,k) = phi(i,j,k) - summreg*regsmooth_group(iregion)/get_dist((/i,j,k/),(/xreg,yreg,zreg/))
      enddo
    else
      do iregion=1,ngroup
        i   = themap_group(1,iregion) 
        j   = themap_group(2,iregion) 
        k   = themap_group(3,iregion) 
        pp  = themap_group(4,iregion)
        irl = themap_group(5,iregion) ; iru = themap_group(6,iregion)
        jrl = themap_group(7,iregion) ; jru = themap_group(8,iregion)
        krl = themap_group(9,iregion) ; kru = themap_group(10,iregion)
        if (.not.luseprevioussum(iregion)) then
          summreg = sum(phirecv(irl:iru,jrl:jru,krl:kru,pp))
        endif
        phi(i,j,k) = phi(i,j,k) - summreg*regsmooth_group(iregion)*regdist1_group(iregion)
      enddo
    endif
!
    if (lroot .and. lshowtime) call cpu_time(tstop_loop2)
!
    phi = phi/(4.0*pi) ! The selfgravity module is going to multiply by a factor
                       ! of 4pi that we don't want (I think).
!
    if (lroot .and. lshowtime) then
      print '("barneshut: MPI time = ",f10.6," seconds.")',tstop_mpi-tstart_mpi
      print '("barneshut: Loop1 time = ",f10.6," seconds.")',tstop_loop1-tstart_loop1
      print '("barneshut: Loop2 time = ",f10.6," seconds.")',tstop_loop2-tstop_loop1
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
    integer, dimension(10) :: newregion
    integer :: sx, sy, sz, ii, ji, ki, il, iu, jl, ju, kl, ku, ppi, newindex
    real :: lxi, lyi, lzi, dist, xwi, ywi, zwi, width
    logical :: laddreg, lsplit, lregionmatched, lnotself
!
    ! If the point where the potential is being calculated is inside the
    ! region in question, then we MUST subdivide further (unless the region is
    ! a single cell-- this is taken care of in the 'if' statement below).
    lsplit = ((ipos(1) .ge. xsind(1) .and. ipos(1) .le. xsind(dimi(1))) &
        .and. (ipos(2) .ge. ysind(1) .and. ipos(2) .le. ysind(dimi(2))) &
        .and. (ipos(3) .ge. zsind(1) .and. ipos(3) .le. zsind(dimi(3))) &
        .and. (iproc .eq. ppi))
    lnotself = (((ipos(1) .ne. xsind(1) .or. ipos(1) .ne. xsind(dimi(1))) &
            .or. (ipos(2) .ne. ysind(1) .or. ipos(2) .ne. ysind(dimi(2))) &
            .or. (ipos(3) .ne. zsind(1) .or. ipos(3) .ne. zsind(dimi(3))) &
            .or. (iproc .ne. ppi)))
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
    elseif (lcylindrical_coords) then
      width = max(xwi,lxi*ywi,zwi)
    elseif (lcartesian_coords) then
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
    elseif (((.not. lsplit).or.(product(dimi).eq.1)) .and. (dist .lt. octree_maxdist) .and. lnotself) then
      newregion = (/ipos(1),ipos(2),ipos(3), ppi, xsind(1), &
       xsind(dimi(1)), ysind(1),ysind(dimi(2)),zsind(1),zsind(dimi(3))/)
      if ((newregion(5).ne.newregion(6)).or.(newregion(7).ne.newregion(8)).or. &
       (newregion(9).ne.newregion(10))) then
        if (laddreg) then
          lregionmatched = .false.
          if (lnorepeatsumming) then
            do iregion=1,ngroup
              if (all(newregion(4:10).eq.themap_group(4:10,iregion))) then
                themap_group(:,iregion+2:ngroup+1) = themap_group(:,iregion+1:ngroup)
                themap_group(:,iregion+1) = newregion
                lregionmatched = .true.
                exit
              endif
            enddo
          endif
          if (.not.lregionmatched) themap_group(:,ngroup+1) = newregion
        endif
        ngroup = ngroup+1
        if (lprecalcdists) regdist1_group(ngroup) = 1.0/dist
        if (dist .gt. (octree_maxdist-octree_smoothdist)) then
         regsmooth_group(ngroup) = 0.5*(cos(pi*((octree_maxdist-dist)/octree_smoothdist+1.0))+1.0)
        endif
      else
        nsingle = nsingle+1
        if (laddreg) then
          themap_single(1:10,nsingle) = newregion(1:10)
          regdist1_single(nsingle) = 1.0/dist
        endif
        if (dist .gt. (octree_maxdist-octree_smoothdist)) then
         regsmooth_single(nsingle) = 0.5*(cos(pi*((octree_maxdist-dist)/octree_smoothdist+1.0))+1.0)
        endif
      endif
    endif
!
    endsubroutine mkmap
!***********************************************************************
    function get_dist(p1,p2) result(rdist)
    integer, dimension(3) :: p1
    real, dimension(3) :: p2
    real :: x2, y2, z2, rdist
    real, dimension(2) :: sincosy, sincosz
!
    if (lspherical_coords) then
      sincosy = sincoslf(p2(2)) ! returns (/sin, cos/)
      sincosz = sincoslf(p2(3))
      x2 = p2(1)*sincosy(1)*sincosz(2)
      y2 = p2(1)*sincosy(1)*sincosz(1)
      z2 = p2(1)*sincosy(2)
    elseif (lcylindrical_coords) then
      sincosy = sincoslf(p2(2))
      x2 = p2(1)*sincosy(2)
      y2 = p2(1)*sincosy(1)
      z2 = p2(3)
    elseif (lcartesian_coords) then
      x2 = p2(1)
      y2 = p2(2)
      z2 = p2(3)
    endif
      rdist = sqrt((xmesh(p1(1),p1(2),p1(3))-x2)**2+ &
        (ymesh(p1(1),p1(2),p1(3))-y2)**2+(zmesh(p1(1),p1(2),p1(3))-z2)**2)
!
    endfunction get_dist
!***********************************************************************
    function sincoslf(ang) ! returns (/sin,cos/)
    real,dimension(2) :: sincoslf
    real :: ang, sinx, cosx !, a0, a1, a2
    integer :: ilk
!
    ilk = nint((ang-lkmin)/dxlt)+1
    sinx = sinlt(ilk)
    cosx = coslt(ilk)
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
    subroutine write_octree
      call system('mkdir -p '//trim(directory_snap)//'/bhmap')
      open(unit = 10, status='replace', &
        file=trim(directory_snap)//trim('/bhmap/themap_single.dat'),form='unformatted')
      write(10) themap_single
      close(10)
      open(unit = 10, status='replace', &
        file=trim(directory_snap)//trim('/bhmap/themap_group.dat'),form='unformatted')
      write(10) themap_group
      close(10)
      open(unit = 10, status='replace', &
        file=trim(directory_snap)//trim('/bhmap/regsmooth_single.dat'),form='unformatted')
      write(10) regsmooth_single
      close(10)
      open(unit = 10, status='replace', &
        file=trim(directory_snap)//trim('/bhmap/regsmooth_group.dat'),form='unformatted')
      write(10) regsmooth_group
      close(10)
      open(unit = 10, status='replace', &
        file=trim(directory_snap)//trim('/bhmap/regdist1_single.dat'),form='unformatted')
      write(10) regdist1_single
      close(10)
      if (lprecalcdists) then
        open(unit = 10, status='replace', &
          file=trim(directory_snap)//trim('/bhmap/regdist1_group.dat'),form='unformatted')
        write(10) regdist1_group
        close(10)
      endif
    endsubroutine write_octree
!***********************************************************************
    subroutine read_octree
      open(unit = 10, status='old', &
        file=trim(directory_snap)//trim('/bhmap/themap_single.dat'),form='unformatted')
      read(10) themap_single
      close(10)
      open(unit = 10, status='old', &
        file=trim(directory_snap)//trim('/bhmap/themap_group.dat'),form='unformatted')
      read(10) themap_group
      close(10)
      open(unit = 10, status='old', &
        file=trim(directory_snap)//trim('/bhmap/regsmooth_single.dat'),form='unformatted')
      read(10) regsmooth_single
      close(10)
      open(unit = 10, status='old', &
        file=trim(directory_snap)//trim('/bhmap/regsmooth_group.dat'),form='unformatted')
      read(10) regsmooth_single
      close(10)
      open(unit = 10, status='old', &
        file=trim(directory_snap)//trim('/bhmap/regdist1_single.dat'),form='unformatted')
      read(10) regdist1_single
      close(10)
      if (lprecalcdists) then
        open(unit = 10, status='old', &
          file=trim(directory_snap)//trim('/bhmap/regdist1_group.dat'),form='unformatted')
        read(10) regdist1_group
        close(10)
      endif
    endsubroutine read_octree
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
    subroutine get_acceleration(acceleration)
!
      use General, only: keep_compiler_quiet
!
      real, dimension(nx,ny,nz,3), intent(out) :: acceleration           !should I (CAN I?) make this allocatable?
!
      call keep_compiler_quiet(acceleration)
!
    endsubroutine get_acceleration
!***********************************************************************
endmodule Poisson
