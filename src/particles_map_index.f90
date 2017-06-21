! $Id$
!
!  This module contains subroutines for particle-mesh related
!  operations.
!
!** AUTOMATIC CPARAM.INC GENERATION ************************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_blocks = .false.
!
!***********************************************************************
module Particles_map
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_mpicomm
!
  implicit none
!
  include "particles_map.h"
!
  interface interpolate_linear
    module procedure interpolate_linear_range
    module procedure interpolate_linear_scalar
  endinterface
!
  public :: distribute_particles, collect_particles
  public :: pm_interpolation, pm_assignment
  public :: map_particles
!
!  Data types for particles in cell.
!
  type, public :: particle
    integer :: proc     ! Parent process rank if a ghost particle; -1, otherwise.
    integer :: id       ! Index in the fp array.
    real :: weight      ! Weight contributed to a cell in code mass unit.
    real :: eps         ! Particle-to-gas density ratio.
    real, dimension(3) :: x, xi ! Particle position in real space and in index space.
    real, dimension(3) :: v, dv ! Particle velocity and its change.
  endtype particle
!
  type, public :: pic
    real :: rho         ! Gas density.
    real, dimension(3) :: u, du ! Gas velocity and its change.
    integer :: np       ! Number of particles in the cell.
    type(particle), dimension(:), pointer :: p  ! Particles in the cell.
  endtype pic
!
!  Module variables
!
  real, dimension(:,:), allocatable :: xi
  real, dimension(3) :: dxi_diag = 0.0
  real :: rinf = 0.0
!
  contains
!***********************************************************************
!***********************************************************************
!  PUBLIC PROCEDURES GO BELOW HERE.
!***********************************************************************
!***********************************************************************
    subroutine initialize_particles_map()
!
!  Perform any post-parameter-read initialization.
!
!  28-apr-16/ccyang: coded.
!
!  Note: Currently, this subroutine is called after modules
!    Particles_mpicomm and Particles.
!
!  Check the particle-mesh interpolation method.
!
      lparticlemesh_gab = .false.
!
      pm: select case (particle_mesh)
!
      case ('ngp', 'NGP') pm
!       Nearest-Grid-Point
        rinf = 0.0
        lparticlemesh_cic = .false.
        lparticlemesh_tsc = .false.
        if (lroot) print *, 'initialize_particles_map: selected nearest-grid-point for particle-mesh method. '
!
      case ('cic', 'CIC') pm
!       Cloud-In-Cell
        rinf = 0.5
        lparticlemesh_cic = .true.
        lparticlemesh_tsc = .false.
        if (lroot) print *, 'initialize_particles_map: selected cloud-in-cell for particle-mesh method. '
!
      case ('tsc', 'TSC') pm
!       Triangular-Shaped-Cloud
        rinf = 1.0
        lparticlemesh_cic = .false.
        lparticlemesh_tsc = .true.
        if (lroot) print *, 'initialize_particles_map: selected triangular-shaped-cloud for particle-mesh method. '
!
      case ('3rd') pm
!       Third-order kernel
        rinf = 1.5
        lparticlemesh_cic = .false.
        lparticlemesh_tsc = .false.
        if (lroot) print *, 'initialize_particles_map: selected third-order kernel for particle-mesh method. '
      case ('6th') pm
!       Third-order kernel
        rinf = 3.0
        lparticlemesh_cic = .false.
        lparticlemesh_tsc = .false.
        if (lroot) print *, 'initialize_particles_map: selected sixth-order kernel for particle-mesh method. '
      case ('etsc', 'ETSC') pm
!       Extended TSC with radius 2
        rinf = 2.0
        lparticlemesh_cic = .false.
        lparticlemesh_tsc = .false.
        if (lroot) print *, 'initialize_particles_map: ', &
                            'selected extended triangular-shaped-cloud with radius 2 for particle-mesh method. '
      case ('etsc2', 'ETSC2') pm
!       Extended TSC with radius 3
        rinf = 3.0
        lparticlemesh_cic = .false.
        lparticlemesh_tsc = .false.
        if (lroot) print *, 'initialize_particles_map: ', &
                            'selected extended triangular-shaped-cloud with radius 3 for particle-mesh method. '
      case default pm
        call fatal_error('initialize_particles_map', "unknown particle-mesh type '" // trim(particle_mesh) // "'")
!
      endselect pm
!
      dxi_diag = merge((/rinf, rinf, rinf/), (/0.0, 0.0, 0.0/), lactive_dimension)
!
    endsubroutine initialize_particles_map
!***********************************************************************
    subroutine collect_particles(f, fp, npsend, sendlist, ghost, cell, ngp_send, ngp_recv, lupdate_par, lupdate_gas, lpmbr)
!
!  Collect the change in velocities and deallocate the working arrays.
!
!  08-may-16/ccyang: coded.
!  21-jun-17/ccyang: accommodated Particles_mass.
!
!  Input/Output Arguments
!      f, fp
!          The f and fp arrays.
!  Input Arguments
!      npsend
!          Number of particles that are required to be communicated.
!      sendlist
!          Indices of the particles to be communicated; deallocated when
!          done.
!      ghost
!          Received particles from other processes; deallocated when
!          done.
!      cell
!          Complete final paricle-in-cell information.
!      ngp_send
!          The i-th element is the number of particles to be sent from
!          this process to process pointed to by iproc_comm(i).
!          The zeroth is that to be sent to itself.
!      ngp_recv
!          The i-th element is the number of particles to be received
!          from process pointed to by iproc(i).
!          The zeroth is that to be received from itself.
!      lupdate_par
!          Update the particle velocities or not.
!      lupdate_gas
!          Update the gas velocity or not.
!      lpmbr
!          Use particle-mesh back reaction or not.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
      type(particle), dimension(:), pointer :: ghost
      integer, dimension(:), pointer :: sendlist
      type(pic), dimension(nx,ny,nz), intent(inout) :: cell
      integer, dimension(0:nproc_comm), intent(in) :: ngp_send, ngp_recv
      integer, intent(in) :: npsend
      logical, intent(in) :: lupdate_par, lupdate_gas, lpmbr
!
      type(particle), dimension(:), allocatable :: packet
      real, dimension(npar_loc,3) :: dmv
      integer :: stat, j
!
!  Collect change in particle velocities.
!
      par: if (lupdate_par) then
        call pic_unset_particles(cell, ghost, dmv)
        call ghost_particles_collect(ghost, ngp_send, ngp_recv, dmv)
        if (lparticles_mass) then
          forall(j = 1:3) fp(1:npar_loc,ivpx+j-1) = fp(1:npar_loc,ivpx+j-1) + dmv(:,j) / fp(1:npar_loc,imp)
        else
          fp(1:npar_loc,ivpx:ivpz) = fp(1:npar_loc,ivpx:ivpz) + dmv / mp_swarm
        endif
      endif par
!
!  Collect change in gas velocity.
!
      gas: if (lupdate_gas) then
        pmbr: if (lpmbr) then
          allocate(packet(npsend), stat=stat)
          if (stat /= 0) call fatal_error_local('collect_particles', 'unable to allocate working array packet.')
          call pack_particles(fp, sendlist, packet)
          forall(j = 1:3) packet%v(j) = dmv(sendlist,j)  ! Use v to send dmv.
          call ghost_particles_send(npsend, packet, xi(sendlist,:), ngp_send, ngp_recv, ghost)
          call back_reaction(f, cell, npar_loc + size(ghost), (/(dmv(:,j), ghost%v(j), j = 1, 3)/), &
                                                              (/(xi(:,j), ghost%xi(j), j = 1, 3)/))
          deallocate(packet, stat=stat)
          if (stat /= 0) call warning('collect_particles', 'unable to deallocate working array packet.')
        else pmbr
          call pic_unset_gas(f, cell)
        endif pmbr
      endif gas
!
!  Deallocate the particles.
!
      call pic_deallocate(int(nw), cell)
      deallocate(sendlist, ghost, stat=stat)
      if (stat /= 0) call warning('collect_particles', 'unable to deallocate ghost particles.')
!
    endsubroutine collect_particles
!***********************************************************************
    subroutine distribute_particles(f, fp, npsend, sendlist, ghost, cell, ngp_send, ngp_recv)
!
!  Distributes weighted particles into each cell.
!
!  08-may-16/ccyang: coded.
!  21-jun-17/ccyang: accommodated Particles_mass.
!
!  Input Arguments
!      f, fp
!          The f and fp arrays.
!  Output Arguments
!      npsend
!          Number of particles that are required to be communicated.
!      sendlist
!          Indices of the particles to be communicated; needs to be
!          deallocated when done.
!      ghost
!          Received particles from other processes; needs to be
!          deallocated when done.
!      cell
!          Complete initial paricle-in-cell information.
!      ngp_send
!          The i-th element is the number of particles to be sent from
!          this process to process pointed to by iproc_comm(i).
!          The zeroth is that to be sent to itself.
!      ngp_recv
!          The i-th element is the number of particles to be received
!          from process pointed to by iproc(i).
!          The zeroth is that to be received from itself.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, intent(out) :: npsend
      integer, dimension(:), pointer :: sendlist
      type(particle), dimension(:), pointer :: ghost
      type(pic), dimension(nx,ny,nz), intent(out) :: cell
      integer, dimension(0:nproc_comm), intent(out) :: ngp_send, ngp_recv
!
      type(particle), dimension(:), allocatable :: packet
      integer :: ngp, stat
!
!  Assign gas properties.
!
      call pic_set_gas(f, cell)
!
!  Send ghost particles.
!
      call ghost_particles_count(xi, sendlist, ngp_send, ngp_recv)
      npsend = size(sendlist)
      ngp = sum(ngp_recv)
      allocate(packet(npsend), ghost(ngp), stat=stat)
      if (stat /= 0) call fatal_error_local('distribute_particles', 'unable to allocate ghost particles.')
      call pack_particles(fp, sendlist, packet)
      call ghost_particles_send(npsend, packet, xi(sendlist,:), ngp_send, ngp_recv, ghost)
!
!  Count particles in each cell.
!
      cell%np = 0
      call pic_count_particles(npar_loc, xi, cell)
      call pic_count_particles(ngp, (/ ghost%xi(1), ghost%xi(2), ghost%xi(3) /), cell)
      call pic_allocate(int(nw), cell)
!
!  Distribute particles into cells.
!
      cell%np = 0
      setpar: if (lparticles_mass) then
        call pic_set_particles(npar_loc, spread(-1,1,npar_loc), xi, fp(1:npar_loc,ixp:izp), fp(1:npar_loc,ivpx:ivpz), cell, &
                fp(1:npar_loc,imp))
        call pic_set_particles(ngp, ghost%proc, (/ ghost%xi(1), ghost%xi(2), ghost%xi(3) /), &
                (/ ghost%x(1), ghost%x(2), ghost%x(3) /), (/ ghost%v(1), ghost%v(2), ghost%v(3) /), cell, ghost%weight)
      else setpar
        call pic_set_particles(npar_loc, spread(-1,1,npar_loc), xi, fp(1:npar_loc,ixp:izp), fp(1:npar_loc,ivpx:ivpz), cell)
        call pic_set_particles(ngp, ghost%proc, (/ ghost%xi(1), ghost%xi(2), ghost%xi(3) /), &
                (/ ghost%x(1), ghost%x(2), ghost%x(3) /), (/ ghost%v(1), ghost%v(2), ghost%v(3) /), cell)
      endif setpar
      call pic_set_eps(cell)
!
      deallocate(packet, stat=stat)
      if (stat /= 0) call warning('distribute_particles', 'unable to deallocate the working array.')
!
    endsubroutine distribute_particles
!***********************************************************************
    subroutine map_nearest_grid(fp, ineargrid)
!
!  Update the coordinates in index space, and round them into ineargrid.
!
!  27-apr-16/ccyang: coded
!
      use Grid, only: real_to_index
!
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, dimension(mpar_loc,3), intent(inout) :: ineargrid
!
      logical :: lflag
      integer :: istat
!
!  Check the status and size of xi.
!
      lflag = .true.
      xisize: if (allocated(xi)) then
        if (size(xi, 1) == npar_loc) then
          lflag = .false.
        else
          deallocate (xi)
        endif
      endif xisize
!
!  (Re)allocate xi as necessary.
!
      alloc: if (lflag) then
        allocate(xi(npar_loc,3), stat=istat)
        if (istat > 0) then
          call fatal_error_local('map_nearest_grid', 'unable to allocate xi. ')
          return
        endif
      endif alloc
!
!  Find the coordinates in index space.
!
      call real_to_index(npar_loc, fp(1:npar_loc,ixp:izp), xi)
      ineargrid(1:npar_loc,:) = nint(xi)
!
    endsubroutine map_nearest_grid
!***********************************************************************
    subroutine map_particles(f, cell)
!
!  Assign the properties of the particles onto the grid for diagnostics.
!
!  08-may-16/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      type(pic), dimension(nx,ny,nz), intent(in) :: cell
!
      integer :: ix, iy, iz, j, l, m, n
      real :: dz1, dyz1, dv1, total_weight
!
!  Loop over each cell.
!
      zscan: do n = n1, n2
        iz = n - nghost
        dz1 = merge(dz_1(n), 1.0, nzgrid > 1)
!
        yscan: do m = m1, m2
          iy = m - nghost
          dyz1 = dz1 * merge(dy_1(m), 1.0, nygrid > 1)
!
          xscan: do l = l1, l2
            ix = l - nghost
            total_weight = sum(cell(ix,iy,iz)%p%weight)
!
!  Assign rhop
!
            rhop: if (irhop > 0) then
              if (total_weight > 0.0) then
                dv1 = dyz1 * merge(dx_1(l), 1.0, nxgrid > 1)
                f(l,m,n,irhop) = dv1 * total_weight
              else
                f(l,m,n,irhop) = 0.0
              endif
            endif rhop
!
!  Assign uup
!
            uup: if (iuup > 0) then
              if (total_weight > 0.0) then
                forall(j = 1:3) &
                    f(l,m,n,iuup+j-1) = sum(cell(ix,iy,iz)%p%weight * (cell(ix,iy,iz)%p%v(j) + cell(ix,iy,iz)%p%dv(j))) &
                                      / total_weight
              else
                f(l,m,n,iupx:iupz) = 0.0
              endif
            endif uup
!
          enddo xscan
        enddo yscan
      enddo zscan
!
    endsubroutine map_particles
!***********************************************************************
    subroutine map_xxp_grid(f, fp, ineargrid, lmapsink_opt)
!
!  Assign particles to rhop and uup.
!
!  08-may-16/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
      logical, intent(in), optional :: lmapsink_opt
!
      type(pic), dimension(nx,ny,nz) :: cell
      type(particle), dimension(:), pointer :: ghost
      integer, dimension(:), pointer :: sendlist
      integer, dimension(0:nproc_comm) :: ngp_send, ngp_recv
      integer :: npsend, stat
!
!  Sanity check.
!
      call keep_compiler_quiet(ineargrid)
      if (irhop == 0 .and. iuup == 0) return
      if (present(lmapsink_opt)) call fatal_error('map_xxp_grid', 'lmapsink_opt is not implemented. ')
      if (lparticles_blocks) call fatal_error('map_xxp_grid', 'particles_blocks version is not implemented yet. ')
!
!  Find np.
!
      if (inp /= 0) call find_np(f, ineargrid)
!
!  Module Particles_drag has done this on the fly.
!
      if (lrun .and. lparticles_drag .and. it > 1) return
!
!  Distribute the particles.
!
      nullify(ghost)
      call distribute_particles(f, fp, npsend, sendlist, ghost, cell, ngp_send, ngp_recv)
!
!  Delegate to map_particles.
!
      call map_particles(f, cell)
!
!  Release the working arrays.
!
      call pic_deallocate(int(nw), cell)
      deallocate(sendlist, ghost, stat=stat)
      if (stat /= 0) call warning('map_xxp_grid', 'unable to deallocate ghost particles.')
!
    endsubroutine map_xxp_grid
!***********************************************************************
    subroutine sort_particles_imn(fp, ineargrid, ipar, dfp, f)
!
!  Sort the particles by pencils.
!
!  27-apr-16/ccyang: coded.
!
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
      integer, dimension(mpar_loc,3), intent(inout) :: ineargrid
      integer, dimension(mpar_loc), intent(inout) :: ipar
      real, dimension(mx,my,mz,mfarray), intent(in), optional :: f
      real, dimension(mpar_loc,mpvar), intent(inout), optional :: dfp
!
      integer, dimension(npar_loc) :: k_sorted
      integer :: i, k
!
      if (present(f)) call keep_compiler_quiet(f)
!
!  Sanity check.
!
      if (lparticles_potential) call fatal_error('sort_particles_imn', 'particles_potential is not supported. ')
!
!  Count particles in each pencil.
!
      npar_imn = 0
      cnt: do k = 1, npar_loc
        i = imn_array(ineargrid(k,2), ineargrid(k,3))
        npar_imn(i) = npar_imn(i) + 1
      enddo cnt
!
!  Calculate the beginning index for the particles in each pencil.
!
      k1_imn(1) = 1
      do i = 2, ny * nz
        k1_imn(i) = k1_imn(i-1) + npar_imn(i-1)
      enddo
      k2_imn = k1_imn - 1
!
!  Distribute the particles in each pencil.
!
      dist: do k = 1, npar_loc
        i = imn_array(ineargrid(k,2), ineargrid(k,3))
        k2_imn(i) = k2_imn(i) + 1
        k_sorted(k2_imn(i)) = k
      enddo dist
!
!  Reset k1_imn and k2_imn for pencils with no particles.
!
      reset: where(npar_imn == 0)
        k1_imn = 0
        k2_imn = 0
      endwhere reset
!
!  Reorganize the particle data.
!
      fp(1:npar_loc,:) = fp(k_sorted,:)
      if (present(dfp)) dfp(1:npar_loc,:) = dfp(k_sorted,:)
      xi = xi(k_sorted,:)
      ineargrid(1:npar_loc,:) = ineargrid(k_sorted,:)
      ipar(1:npar_loc) = ipar(k_sorted)
!
   endsubroutine sort_particles_imn
!***********************************************************************
!***********************************************************************
!  LOCAL SUBROUTINES GO BELOW HERE.
!***********************************************************************
!***********************************************************************
    subroutine back_reaction(f, cell, np, dmv, xi)
!
!  Particle-mesh the momentum changes in particles back to the gas and
!  update the gas velocity.
!
!  20-sep-15/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      type(pic), dimension(nx,ny,nz), intent(in) :: cell
      integer, intent(in) :: np
      real, dimension(np,3), intent(in) :: dmv, xi
!
      real, dimension(mx,my,mz) :: dp
      integer :: i
!
!  Operate on each component.
!
      comp: do i = 1, 3
        call pm_assignment(np, dmv(:,i), xi, dp)
        f(l1:l2,m1:m2,n1:n2,iuu+i-1) = f(l1:l2,m1:m2,n1:n2,iuu+i-1) - dp(l1:l2,m1:m2,n1:n2) / cell%rho + cell%du(i)
      enddo comp
!
    endsubroutine back_reaction
!***********************************************************************
    pure subroutine block_of_influence(xi, lower_corner, upper_corner, prune, domain)
!
!  Check the block of influence from a particle and return the two
!  diagonal corners of the block in terms of the cell indices.
!
!  If the optional argument prune is present and is .true., the portion
!  outside of the domain will be trimmed off the block.
!
!  If the optional argument domain is present and is .true., the corners
!  of the block are in terms of the domain sides: The value -1 denote
!  outside of the lower boundary, 0 inside the domain, and +1 outside of
!  the upper boundary.
!
!  10-sep-15/ccyang: coded.
!
      real, dimension(3), intent(in) :: xi
      integer, dimension(3), intent(out) :: lower_corner, upper_corner
      logical, intent(in), optional :: prune, domain
!
!##ccyang
!      integer, dimension(3), parameter :: xi1 = (/ l1, m1, n1 /)
!      integer, dimension(3), parameter :: xi2 = (/ l2, m2, n2 /)
      integer, dimension(3) :: xi1, xi2
      xi1 = (/ l1, m1, n1 /)
      xi2 = (/ l2, m2, n2 /)
!##ccyang
!
!  Find the indices of the corners.
!
      lower_corner = nint(xi - dxi_diag)
      upper_corner = nint(xi + dxi_diag)
!
!  Trim the block, if requested.
!
      opt1: if (present(prune)) then
        clip: if (prune) then
          lower_corner = min(max(lower_corner, xi1), xi2)
          upper_corner = min(max(upper_corner, xi1), xi2)
        endif clip
      endif opt1
!
!  Convert to domain sides, if requested.
!
      opt2: if (present(domain)) then
        conv: if (domain) then
          lower_corner = check_side(lower_corner, xi1, xi2)
          upper_corner = check_side(upper_corner, xi1, xi2)
        endif conv
      endif opt2
!
    endsubroutine block_of_influence
!***********************************************************************
    pure subroutine find_np(f, ineargrid)
!
!  Count particles in each cell.
!
!  10-may-16/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
!
      integer :: ix0, iy0, iz0, k
!
      f(:,:,:,inp) = 0.0
      incr: do k = 1, npar_loc
        ix0 = ineargrid(k,1)
        iy0 = ineargrid(k,2)
        iz0 = ineargrid(k,3)
        f(ix0,iy0,iz0,inp) = f(ix0,iy0,iz0,inp) + 1.0
      enddo incr
!
    endsubroutine find_np
!***********************************************************************
    subroutine ghost_particles_collect(ghost, ngp_send, ngp_recv, dmv)
!
!  Communicate the momentum change of the ghost particles and add it to
!  dmv.
!
!  28-feb-16/ccyang: coded.
!
      use Mpicomm, only: mpisend_int, mpirecv_int, mpisend_real, mpirecv_real
!
      type(particle), dimension(:), intent(in) :: ghost
      integer, dimension(0:nproc_comm), intent(in) :: ngp_send, ngp_recv
      real, dimension(npar_loc,3), intent(inout) :: dmv
!
      integer, parameter :: tag = 100
      integer, dimension(maxval(ngp_send)) :: id
      integer, dimension(maxval(ngp_recv)) :: ibuf
      real, dimension(3*maxval(ngp_send)) :: dp
      real, dimension(3*maxval(ngp_recv)) :: rbuf
      integer :: ip, iproc_target, nsend, nrecv
      integer :: i, j, k, n3s, n3r
!
!  Shake hands with each process.
!
      j = 0
      proc: do ip = 0, nproc_comm
        if (ip > 0) then
          iproc_target = iproc_comm(ip)
        else
          iproc_target = iproc
        endif
        nsend = ngp_recv(ip)
        nrecv = ngp_send(ip)
        k = j + nsend
        j = j + 1
        n3s = 3 * nsend
        n3r = 3 * nrecv
        if (any(ghost(j:k)%proc /= iproc_target)) &
            call fatal_error_local('ghost_particles_collect', 'inconsistant number of particles to be sent. ')
        ibuf(1:nsend) = ghost(j:k)%id
        rbuf(1:n3s) = (/ ghost(j:k)%dv(1), ghost(j:k)%dv(2), ghost(j:k)%dv(3) /)
!
!  Send the particles.
!
        comm: if (iproc_target > iproc) then
          call mpisend_int(ibuf, nsend, iproc_target, tag+iproc_target)
          call mpisend_real(rbuf, n3s, iproc_target, tag+iproc_target)
          call mpirecv_int(id, nrecv, iproc_target, tag+iproc)
          call mpirecv_real(dp, n3r, iproc_target, tag+iproc)
        elseif (iproc_target < iproc) then comm
          call mpirecv_int(id, nrecv, iproc_target, tag+iproc)
          call mpirecv_real(dp, n3r, iproc_target, tag+iproc)
          call mpisend_int(ibuf, nsend, iproc_target, tag+iproc_target)
          call mpisend_real(rbuf, n3s, iproc_target, tag+iproc_target)
        else comm
          id(1:nrecv) = ibuf(1:nsend)
          dp(1:n3r) = rbuf(1:n3s)
        endif comm
!
!  Collect the change of particle momenta.
!
        do i = 1, nrecv
          dmv(id(i),:) = dmv(id(i),:) + dp(i:n3r:nrecv)
        enddo
        j = k
      enddo proc
!
    endsubroutine ghost_particles_collect
!***********************************************************************
    subroutine ghost_particles_count(xi, ighost, ngp_send, ngp_recv)
!
!  Count the ghost particles to be communicated.
!
!  28-feb-16/ccyang: coded.
!
      use Mpicomm, only: mpisend_int, mpirecv_int
      use Sub, only: get_where
!
      real, dimension(npar_loc,3), intent(in) :: xi
      integer, dimension(:), pointer :: ighost
      integer, dimension(0:nproc_comm), intent(out) :: ngp_send, ngp_recv
!
      integer, parameter :: tag = 100
      logical, dimension(npar_loc) :: lghost
      integer, dimension(-1:1,-1:1,-1:1) :: neighbor_send
      integer :: ip, iproc_target, j
!
!  Tag the particles near the domain boundary and count the directions to be sent.
!
      lghost = .false.
      ngp_send = 0
      par: do ip = 1, npar_loc
        call tag_send_directions(xi(ip,:), neighbor_send)
        cnt: if (any(neighbor_send >= 0)) then
          lghost(ip) = .true.
          forall(j = 0:nproc_comm) ngp_send(j) = ngp_send(j) + count(neighbor_send == j)
        endif cnt
      enddo par
!
!  Get the indices of the would-be ghost particles.
!
      call get_where(lghost, ighost)
!
!  Communicate the number of ghost particles.
!
      ngp_recv(0) = ngp_send(0)
      proc: do ip = 1, nproc_comm
        iproc_target = iproc_comm(ip)
        comm: if (iproc_target > iproc) then
          call mpisend_int(ngp_send(ip), iproc_target, tag + iproc_target)
          call mpirecv_int(ngp_recv(ip), iproc_target, tag + iproc)
        else comm
          call mpirecv_int(ngp_recv(ip), iproc_target, tag + iproc)
          call mpisend_int(ngp_send(ip), iproc_target, tag + iproc_target)
        endif comm
      enddo proc
!
    endsubroutine ghost_particles_count
!***********************************************************************
    subroutine ghost_particles_send(np, packet, xi, ngp_send, ngp_recv, ghost)
!
!  Send particles near the boundary.
!
!  28-feb-16/ccyang: coded.
!  21-jun-17/ccyang: accommodated Particles_mass.
!
      use Grid, only: real_to_index
      use Mpicomm, only: mpisend_int, mpirecv_int, mpisend_real, mpirecv_real
!
      integer, intent(in) :: np
      type(particle), dimension(np), intent(in) :: packet
      real, dimension(np,3), intent(in) :: xi
      integer, dimension(0:nproc_comm), intent(in) :: ngp_send, ngp_recv
      type(particle), dimension(:), intent(out) :: ghost
!
      type :: buffer
        integer, dimension(:), pointer :: ibuf
        real, dimension(:), pointer :: rbuf
      endtype buffer
!
      integer, parameter :: tag = 100
      real, dimension(:,:), allocatable :: xigp
      integer, dimension(-1:1,-1:1,-1:1) :: neighbor_send
      type(buffer), dimension(0:nproc_comm) :: sendbuf
      integer, dimension(maxval(ngp_recv)) :: ibuf
      integer, dimension(0:nproc_comm) :: ngp
      real, dimension(7*maxval(ngp_recv)) :: rbuf
      real, dimension(3) :: xp
      integer :: nattr, neighbor, nsend, nrecv, stat
      integer :: ip, iproc_target, i, j, k, m, n
!
!  Set number of attributes to send.
!
      nattr = 6
      if (lparticles_mass) nattr = nattr + 1
!
!  Allocate buffers.
!
      alloc: do ip = 0, nproc_comm
        allocate(sendbuf(ip)%ibuf(ngp_send(ip)), sendbuf(ip)%rbuf(nattr*ngp_send(ip)), stat=stat)
        if (stat /= 0) call fatal_error_local('ghost_particles_send', 'cannot allocate the send buffers. ')
      enddo alloc
!
!  Prepare the buffers.
!
      ngp = 0
      par: do ip = 1, np
        call tag_send_directions(xi(ip,:), neighbor_send)
        zscan: do k = -1, 1
          xscan: do i = -1, 1
            yscan: do j = -1, 1
              neighbor = neighbor_send(i,j,k)
              valid: if (neighbor >= 0) then
                n = ngp(neighbor)
                sendbuf(neighbor)%ibuf(n+1) = packet(ip)%id
                xp = packet(ip)%x
                call wrap_particle_position(xp, (/i,j,k/))
                if (lparticles_mass) then
                  sendbuf(neighbor)%rbuf(nattr*n+1:nattr*(n+1)) = (/ xp, packet(ip)%v, packet(ip)%weight /)
                else
                  sendbuf(neighbor)%rbuf(nattr*n+1:nattr*(n+1)) = (/ xp, packet(ip)%v /)
                endif
                ngp(neighbor) = n + 1
              endif valid
            enddo yscan
          enddo xscan
        enddo zscan
      enddo par
!
      n = 0
      proc: do ip = 0, nproc_comm
        nsend = ngp(ip)
        nrecv = ngp_recv(ip)
        m = nattr * nrecv
        if (nsend /= ngp_send(ip)) call fatal_error_local('ghost_particles_send', 'inconsistant number of particles to be sent. ')
!
!  Send the particles.
!
        if (ip > 0) then
          iproc_target = iproc_comm(ip)
        else
          iproc_target = iproc
        endif
!
        comm: if (ip == 0) then
          ibuf(1:nrecv) = sendbuf(ip)%ibuf
          rbuf(1:m) = sendbuf(ip)%rbuf
        elseif (iproc_target > iproc) then comm
          call mpisend_int(sendbuf(ip)%ibuf, nsend, iproc_target, tag+iproc_target)
          call mpisend_real(sendbuf(ip)%rbuf, nattr*nsend, iproc_target, tag+iproc_target)
          call mpirecv_int(ibuf, nrecv, iproc_target, tag+iproc)
          call mpirecv_real(rbuf, m, iproc_target, tag+iproc)
        elseif (iproc_target < iproc) then comm
          call mpirecv_int(ibuf, nrecv, iproc_target, tag+iproc)
          call mpirecv_real(rbuf, m, iproc_target, tag+iproc)
          call mpisend_int(sendbuf(ip)%ibuf, nsend, iproc_target, tag+iproc_target)
          call mpisend_real(sendbuf(ip)%rbuf, nattr*nsend, iproc_target, tag+iproc_target)
        endif comm
!
!  Assemble the particle data.
!
        j = n + 1
        k = n + nrecv
        ghost(j:k)%proc = iproc_target
        ghost(j:k)%id = ibuf(1:nrecv)
        if (lparticles_mass) ghost(j:k)%weight = rbuf(7:m:nattr)
        comp: forall(i = 1:3)
          ghost(j:k)%x(i) = rbuf(i:m:nattr)
          ghost(j:k)%v(i) = rbuf(i+3:m:nattr)
        endforall comp
        n = k
      enddo proc
!
!  Convert the particle positions in real space to those in index space.
!
      allocate(xigp(n,3))
      call real_to_index(n, (/ ghost%x(1), ghost%x(2), ghost%x(3) /), xigp)
      forall(i = 1:3) ghost%xi(i) = xigp(:,i)
!
!  Deallocate working arrays.
!
      deallocate(xigp, stat=stat)
      if (stat /= 0) call warning('ghost_particles_send', 'cannot deallocate xigp. ')
      dealloc: do ip = 0, nproc_comm
        deallocate(sendbuf(ip)%ibuf, sendbuf(ip)%rbuf, stat=stat)
        if (stat /= 0) call warning('ghost_particles_send', 'cannot deallocate the send buffers. ')
      enddo dealloc
!
    endsubroutine ghost_particles_send
!***********************************************************************
    subroutine pack_particles(fp, list, packet)
!
!  Pack the information of selected particles into a packet.
!
!  18-sep-15/ccyang: coded.
!  21-jun-17/ccyang: pack particle mass if present.
!
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, dimension(:), intent(in) :: list
      type(particle), dimension(size(list)), intent(out) :: packet
!
      integer :: i
!
      packet%proc = iproc
      packet%id = list
      if (lparticles_mass) then
        packet%weight = fp(list, imp)
      else
        packet%weight = 0.0
      endif
      packet%eps = 0.0
      comp: forall(i = 1:3)
        packet%x(i) = fp(list, ixp+i-1)
        packet%v(i) = fp(list, ivpx+i-1)
        packet%dv(i) = 0.0
      endforall comp
!
    endsubroutine pack_particles
!***********************************************************************
    subroutine pic_allocate(n, cell)
!
!  Allocate space for particles given known number of particles.
!
!  08-feb-15/ccyang: coded.
!
      integer, intent(in) :: n
      type(pic), dimension(n), intent(inout) :: cell
!
      integer :: i, istat
!
!  Loop over each cell.
!
      loop: do i = 1, n
        allocate(cell(i)%p(cell(i)%np), stat=istat)
        if (istat /= 0) call fatal_error_local('pic_allocate', 'unable to allocate particle array.')
      enddo loop
!
    endsubroutine pic_allocate
!***********************************************************************
    pure subroutine pic_count_particles(npar, xi, cell)
!
!  Count and accumulate the number of particles that have influence in
!  each cell.
!
!  11-feb-15/ccyang: coded.
!
      integer, intent(in) :: npar
      real, dimension(npar,3), intent(in) :: xi
      type(pic), dimension(nx,ny,nz), intent(inout) :: cell
!
      integer, dimension(3) :: xi1, xi2
      integer :: ip, ix, iy, iz, l, m, n
!
      par: do ip = 1, npar
        call block_of_influence(xi(ip,:), xi1, xi2, prune=.true.)
        zscan: do n = xi1(3), xi2(3)
          iz = n - nghost
          yscan: do m = xi1(2), xi2(2)
            iy = m - nghost
            xscan: do l = xi1(1), xi2(1)
              ix = l - nghost
              cell(ix,iy,iz)%np = cell(ix,iy,iz)%np + 1
            enddo xscan
          enddo yscan
        enddo zscan
      enddo par
!
    endsubroutine pic_count_particles
!***********************************************************************
    subroutine pic_deallocate(n, cell)
!
!  Deallocate space for particles.
!
!  09-feb-15/ccyang: coded.
!
      integer, intent(in) :: n
      type(pic), dimension(n), intent(inout) :: cell
!
      integer :: i, istat
!
!  Loop over each cell.
!
      loop: do i = 1, n
        dealloc: if (associated(cell(i)%p)) then
          deallocate(cell(i)%p, stat=istat)
          if (istat /= 0) call warning('pic_deallocate', 'unable to deallocate particle array.')
          nullify(cell(i)%p)
        endif dealloc
      enddo loop
!
    endsubroutine pic_deallocate
!***********************************************************************
    pure subroutine pic_set_eps(cell)
!
!  Find and set the particle-to-gas density ratios.
!
!  11-feb-15/ccyang: coded
!
      type(pic), dimension(nx,ny,nz), intent(inout) :: cell
!
      integer :: ix, iy, iz, l, m, n
      real :: dz1, dyz1, dv1
!
!  Loop over each cell.
!
      zscan: do n = n1, n2
        iz = n - nghost
        dz1 = merge(dz_1(n), 1.0, nzgrid > 1)
        yscan: do m = m1, m2
          iy = m - nghost
          dyz1 = dz1 * merge(dy_1(m), 1.0, nygrid > 1)
          xscan: do l = l1, l2
            ix = l - nghost
            dv1 = dyz1 * merge(dx_1(l), 1.0, nxgrid > 1)
            cell(ix,iy,iz)%p%eps = dv1 / cell(ix,iy,iz)%rho * cell(ix,iy,iz)%p%weight
          enddo xscan
        enddo yscan
      enddo zscan
!
    endsubroutine pic_set_eps
!***********************************************************************
    subroutine pic_set_gas(f, cell)
!
!  Assigns gas properties to each cell.
!
!  08-feb-15/ccyang: coded.
!
      use Particles_sub, only: get_gas_density
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      type(pic), dimension(nx,ny,nz), intent(inout) :: cell
!
      integer :: ix, iy, iz, l, m, n
!
!  Loop over each cell.
!
      zscan: do n = n1, n2
        iz = n - nghost
        yscan: do m = m1, m2
          iy = m - nghost
          xscan: do l = l1, l2
            ix = l - nghost
            cell(ix,iy,iz)%rho = get_gas_density(f, l, m, n)
            if (iuu /= 0) then
              cell(ix,iy,iz)%u = f(l,m,n,iux:iuz)
            else
              cell(ix,iy,iz)%u = 0.0
            endif
            cell(ix,iy,iz)%du = 0.0
          enddo xscan
        enddo yscan
      enddo zscan
!
    endsubroutine pic_set_gas
!***********************************************************************
    pure subroutine pic_set_particles(npar, proc, xi, xp, vp, cell, mass)
!
!  Set the properties of the particles in each cell.
!
!  20-sep-15/ccyang: coded.
!  21-jun-17/ccyang: added optional argument mass.
!
      integer, intent(in) :: npar
      integer, dimension(npar), intent(in) :: proc
      real, dimension(npar,3), intent(in) :: xi, xp, vp
      type(pic), dimension(nx,ny,nz), intent(inout) :: cell
      real, dimension(npar), intent(in), optional :: mass
!
      real, dimension(3) :: dxi
      integer, dimension(3) :: xi1, xi2
      integer :: np, ip, ix, iy, iz, l, m, n
      real :: mp
!
!  Weigh and distribute particles to the cells.
!
      par: do ip = 1, npar
        if (present(mass)) then
          mp = mass(ip)
        else
          mp = mp_swarm
        endif
!
        call block_of_influence(xi(ip,:), xi1, xi2, prune=.true.)
!
        zscan: do n = xi1(3), xi2(3)
          iz = n - nghost
          dxi(3) = xi(ip,3) - real(n)
!
          yscan: do m = xi1(2), xi2(2)
            iy = m - nghost
            dxi(2) = xi(ip,2) - real(m)
!
            xscan: do l = xi1(1), xi2(1)
              ix = l - nghost
              dxi(1) = xi(ip,1) - real(l)
!
              np = cell(ix,iy,iz)%np + 1
              cell(ix,iy,iz)%p(np)%proc = proc(ip)
              cell(ix,iy,iz)%p(np)%id = ip
              cell(ix,iy,iz)%p(np)%weight = mp * weigh_particle(dxi(1), dxi(2), dxi(3))
              cell(ix,iy,iz)%p(np)%xi = xi(ip,:)
              cell(ix,iy,iz)%p(np)%x = xp(ip,:)
              cell(ix,iy,iz)%p(np)%v = vp(ip,:)
              cell(ix,iy,iz)%p(np)%dv = 0.0
              cell(ix,iy,iz)%np = np
            enddo xscan
          enddo yscan
        enddo zscan
      enddo par
!
    endsubroutine pic_set_particles
!***********************************************************************
    subroutine pic_unset_gas(f, cell)
!
!  Add the change of the gas velocity back to the f-array.
!
!  12-feb-15/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      type(pic), dimension(nx,ny,nz), intent(in) :: cell
!
      integer :: k
!
!  Loop over each compnent.
!
      forall(k = 1:3) f(l1:l2,m1:m2,n1:n2,iuu+k-1) = f(l1:l2,m1:m2,n1:n2,iuu+k-1) + cell%du(k)
!
    endsubroutine pic_unset_gas
!***********************************************************************
    pure subroutine pic_unset_particles(cell, ghost, dmv)
!
!  Collect the change of particle momenta over each cell into dmv and
!  ghost%dv.
!
!  10-jun-15/ccyang: coded.
!
      type(pic), dimension(nx,ny,nz), intent(in) :: cell
      type(particle), dimension(:), intent(inout) :: ghost
      real, dimension(npar_loc,3), intent(out) :: dmv
!
      real, dimension(3) :: dp
      integer :: id, ix, iy, iz, j
!
!  Initialization.
!
      dmv = 0.0
      forall(j = 1:3) ghost%dv(j) = 0.0
!
!  Process each cell.
!
      zscan: do iz = 1, nz
        yscan: do iy = 1, ny
          xscan: do ix = 1, nx
            par: do j = 1, cell(ix,iy,iz)%np
              id = cell(ix,iy,iz)%p(j)%id
              dp = cell(ix,iy,iz)%p(j)%weight * cell(ix,iy,iz)%p(j)%dv
              if (cell(ix,iy,iz)%p(j)%proc < 0) then
                dmv(id,:) = dmv(id,:) + dp
              else
                ghost(id)%dv = ghost(id)%dv + dp
              endif
            enddo par
          enddo xscan
        enddo yscan
      enddo zscan
!
    endsubroutine pic_unset_particles
!***********************************************************************
    pure subroutine pm_assignment(np, fp, xi, fg)
!
!  Assign the property fp at (particle) positions xi in index space onto
!  fg on the grid.  Note that fp is converted to density in fg, being
!  divided by the cell volume.  Also, ghost particles should be
!  included in the total number np.
!
!  20-sep-15/ccyang: coded.
!
      integer, intent(in) :: np
      real, dimension(mx,my,mz), intent(out) :: fg
      real, dimension(np,3), intent(in) :: xi
      real, dimension(np), intent(in) :: fp
!
      integer, dimension(3) :: xi1, xi2
      real, dimension(3) :: dxi
      integer :: ip, l, m, n
      real :: w, dz1, dyz1, dv1
!
!  Make the assignment.
!
      fg = 0.0
      loop: do ip = 1, np
        call block_of_influence(xi(ip,:), xi1, xi2)
        zscan: do n = xi1(3), xi2(3)
          dxi(3) = xi(ip,3) - real(n)
          dz1 = merge(dz_1(n), 1.0, nzgrid > 1)
          yscan: do m = xi1(2), xi2(2)
            dxi(2) = xi(ip,2) - real(m)
            dyz1 = dz1 * merge(dy_1(m), 1.0, nygrid > 1)
            xscan: do l = xi1(1), xi2(1)
              dxi(1) = xi(ip,1) - real(l)
              dv1 = dyz1 * merge(dx_1(l), 1.0, nxgrid > 1)
              w = weigh_particle(dxi(1), dxi(2), dxi(3))
              fg(l,m,n) = fg(l,m,n) + dv1 * w * fp(ip)
            enddo xscan
          enddo yscan
        enddo zscan
      enddo loop
!
    endsubroutine pm_assignment
!***********************************************************************
    pure subroutine pm_interpolation(np, fg, xi, fp)
!
!  Interpolate the field fg on the grid to fp at the (particle)
!  positions xi in index space.  Ghost cells are assumed updated.
!
!  05-jun-15/ccyang: coded.
!
      integer, intent(in) :: np
      real, dimension(mx,my,mz), intent(in) :: fg
      real, dimension(np,3), intent(in) :: xi
      real, dimension(np), intent(out) :: fp
!
      integer, dimension(3) :: xi1, xi2
      real, dimension(3) :: dxi
      integer :: ip, l, m, n
      real :: w, a
!
!  Make the interpolation.
!
      loop: do ip = 1, np
        a = 0.0
        call block_of_influence(xi(ip,:), xi1, xi2)
        zscan: do n = xi1(3), xi2(3)
          dxi(3) = xi(ip,3) - real(n)
          yscan: do m = xi1(2), xi2(2)
            dxi(2) = xi(ip,2) - real(m)
            xscan: do l = xi1(1), xi2(1)
              dxi(1) = xi(ip,1) - real(l)
              w = weigh_particle(dxi(1), dxi(2), dxi(3))
              a = a + w * fg(l,m,n)
            enddo xscan
          enddo yscan
        enddo zscan
        fp(ip) = a
      enddo loop
!
    endsubroutine pm_interpolation
!***********************************************************************
    subroutine tag_send_directions(xi, neighbor_send)
!
!  Specifies the directions and the process ranks a particle needs to be
!  sent.
!
!  28-feb-16/ccyang: coded.
!
!  Input Argument
!      xi
!          The position of the particle in local index space.
!  Output Argument
!      neighbor_send
!          The element (i,j,k) contains the index to iproc_comm in the
!          direction (i,j,k), 0 if send to itself, -1 if not a valid
!          direction.
!
      use General, only: find_proc
      use Mpicomm, only: index_to_iproc_comm
!
      real, dimension(3), intent(in) :: xi
      integer, dimension(-1:1,-1:1,-1:1), intent(out) :: neighbor_send
!
      real, parameter :: xi0 = nghost + 0.5
      logical, dimension(-1:1,-1:1,-1:1) :: lsend
      integer, dimension(3) :: side1, side2
      logical :: lflag
      integer :: ipxn, ipyn1, ipyn2, ipzn
      integer :: ix, iy, iz, n
      real :: eta
!
!  Disable all directions initially.
!
      neighbor_send = -1
!
!  Find the block of influence from the particle.
!
      call block_of_influence(xi, side1, side2, domain=.true.)
!
!  Check if the block boundary is outside of the local domain in each direction.
!
      xscan: do ix = side1(1), side2(1)
        shear: if (is_shear_boundary(ix)) then
!         Is a shear boundary.
          if (.not. lequidist(2)) &
              call fatal_error_local('tag_send_directions', 'not implemented for shear with nonequidistant y grid. ')
          if (llast_proc_x) then
            ipxn = 0
          else
            ipxn = nprocx - 1
          endif
          eta = modulo(xi(2) - xi0 + real(ny * ipy) + real(ix) * deltay / dy, real(nygrid)) - 0.5
          ipyn1 = find_ipy(eta - rinf)
          ipyn2 = find_ipy(eta + rinf)
          lflag = ipyn1 /= ipyn2
          ipyn1 = modulo(ipyn1, nprocy)
          ipyn2 = modulo(ipyn2, nprocy)
          zscan1: do iz = side1(3), side2(3)
            ipzn = modulo(ipz + iz, nprocz)
            setn: if (lflag) then
              neighbor_send(ix,-1,iz) = find_proc(ipxn,ipyn1,ipzn)
              neighbor_send(ix,+1,iz) = find_proc(ipxn,ipyn2,ipzn)
            else setn
              neighbor_send(ix,0,iz) = find_proc(ipxn,ipyn1,ipzn)
            endif setn
          enddo zscan1
        else shear
!         Is not a shear boundary.
          zscan2: do iz = side1(3), side2(3)
            yscan: do iy = side1(2), side2(2)
              if (ix == 0 .and. iz == 0 .and. iy == 0) cycle
              n = neighbors_par(ix,iy,iz)
              if (n >= 0) neighbor_send(ix,iy,iz) = n
            enddo yscan
          enddo zscan2
        endif shear
      enddo xscan
!
!  Transform the process ranks to indices to iproc_comm.
!
      lsend = neighbor_send >= 0
      neighbor_send = index_to_iproc_comm(neighbor_send, lsend)
      if (any(lsend .and. neighbor_send < 0)) &
          call fatal_error_local('tag_send_directions', 'sending particles to non-neighbor process')
!
    endsubroutine tag_send_directions
!***********************************************************************
    subroutine wrap_particle_position(xp, sides)
!
!  Updates the position of a particle if it will be sent over to the
!  opposite side of the computational domain.
!
!  14-sep-15/ccyang: coded.
!
!  Input/Output Argument
!      xp(i)
!          The i-th component of the particle position.
!  Input Argument
!      sides(i)
!          The side in the i-th dimension the particle is sent over,
!          where only +1, 0, and -1 are accepted.
!
      use Grid, only: inverse_grid
!
      real, dimension(3), intent(inout) :: xp
      integer, dimension(3), intent(in) :: sides
!
      logical :: lflag
      real, dimension(1) :: xi
      integer :: n
!
!  Directions in x.
!
      lflag = .false.
      xdir: if (lfirst_proc_x .and. sides(1) == -1) then
        xp(1) = xp(1) + Lxyz(1)
        lflag = .true.
      elseif (llast_proc_x .and. sides(1) == +1) then xdir
        xp(1) = xp(1) - Lxyz(1)
        lflag = .true.
      endif xdir
!
!  Directions in y with shear
!
      ydir: if (lshear .and. nygrid > 1 .and. lflag) then
        xp(2) = xyz0(2) + modulo(xp(2) + real(sides(1)) * deltay - xyz0(2), lxyz(2))
        call inverse_grid(2, (/xp(2)/), xi)
        n = nint(xi(1) + real(sides(2)) * rinf) - nghost
        if (n < 1) then
          xp(2) = xp(2) + lxyz(2)
        elseif (n > nygrid) then
          xp(2) = xp(2) - lxyz(2)
        endif
!
!  or without shear.
!
      elseif (lfirst_proc_y .and. sides(2) == -1) then ydir
        xp(2) = xp(2) + Lxyz(2)
      elseif (llast_proc_y .and. sides(2) == +1) then ydir
        xp(2) = xp(2) - Lxyz(2)
      endif ydir
!
!  Directions in z.
!
      if (lfirst_proc_z .and. sides(3) == -1) then
        xp(3) = xp(3) + Lxyz(3)
      elseif (llast_proc_z .and. sides(3) == +1) then
        xp(3) = xp(3) - Lxyz(3)
      endif
!
    endsubroutine wrap_particle_position
!***********************************************************************
!***********************************************************************
!  LOCAL FUNCTIONS GO BELOW HERE.
!***********************************************************************
!***********************************************************************
    elemental integer function check_side(n, n1, n2)
!
!  Return -1, if n < n1; 0, if n1 <= n <= n2; +1, otherwise.
!
!  08-feb-15/ccyang: coded
!
      integer, intent(in) :: n, n1, n2
!
      if (n < n1) then
        check_side = -1
      else if (n <= n2) then
        check_side = 0
      else
        check_side = 1
      endif
!
    endfunction check_side
!***********************************************************************
    pure integer function find_ipy(eta)
!
!  Helper function for subroutine tag_send_directions.
!
!  01-mar-16/ccyang: coded.
!
      real, intent(in) :: eta
!
      integer :: k
!
      k = floor(eta + 0.5)
      find_ipy = (k - modulo(k, ny)) / ny
!
    endfunction find_ipy
!***********************************************************************
    logical function is_shear_boundary(xside)
!
!  Returns .true. if the x-boundary in the xside side is a shear
!  boundary, where xside can only be +1, 0, or -1.
!
!  13-sep-15/ccyang: coded.
!
      integer, intent(in) :: xside
!
      is_shear_boundary = lshear .and. nygrid > 1 .and. &
          (lfirst_proc_x .and. xside == -1 .or. llast_proc_x .and. xside == +1)
!
    endfunction is_shear_boundary
!***********************************************************************
    elemental real function particlemesh_weighting(dxi) result(weight)
!
!  Returns a normalized weight of a particle according to its distance
!  from the cell center in index space in one dimension.
!
!  13-apr-15/ccyang: coded.
!
      real, intent(in) :: dxi
!
      real :: x
!
      pm: select case (particle_mesh)
!
!  Nearest-Grid-Point scheme.
!
      case ('ngp', 'NGP') pm
        weight = merge(1.0, 0.0, -0.5 <= dxi .and. dxi < 0.5)
!
!  Cloud-In-Cell scheme
!
      case ('cic', 'CIC') pm
        weight = max(1.0 - abs(dxi), 0.0)
!
!  Triangular-Shaped-Cloud scheme
!
      case ('tsc', 'TSC') pm
        x = abs(dxi)
        if (x < 0.5) then
          weight = 0.75 - x**2
        elseif (x < 1.5) then
          weight = 0.5 * (1.5 - x)**2
        else
          weight = 0.0
        endif
!
!  Third-order
!
      case ('3rd') pm
        x = abs(dxi)
        if (x < 1.0) then
          weight = (4.0 - 3.0 * (2.0 - x) * x**2) / 6.0
        elseif (x < 2.0) then
          weight = (2.0 - x)**3 / 6.0
        else
          weight = 0.0
        endif
!
!  Sixth-order
!
      case ('6th') pm
        x = abs(dxi)
        if (x < 0.5) then
          x = x**2
          weight = (5.887e3 - 20.0 * x * (2.31e2 - x * (84.0 - 16.0 * x))) / 1.152e4
        elseif (x < 1.5) then
          weight = (2.3583e4 + x * (-420.0 + x * (-1.638e4 + x * (-5.6e3 + x * (1.512e4 + x * (-6.72e3 + 960.0 * x)))))) / 4.608e4
        elseif (x < 2.5) then
          weight = (4.137e3 + x*(3.0408e4 + x*(-5.922e4 + x*(4.256e4 + x*(-1.512e4 + x*(2.688e3 - 192.0 * x)))))) / 2.304e4
        elseif (x < 3.5) then
          weight = (7.0 - 2.0 * x)**6 / 4.608e4
        else
          weight = 0.0
        endif
!
!  Extended TSC with radius 2.
!
      case ('etsc', 'ETSC') pm
        x = abs(dxi)
        if (x < 0.5) then
          weight = 0.4375 - 0.25 * x**2
        elseif (x < 1.5) then
          weight = 0.5 - 0.25 * x
        elseif (x < 2.5) then
          weight = 0.125 * (2.5 - x)**2
        else
          weight = 0.0
        endif
!
!  Extended TSC with radius 3.
!
      case ('etsc2', 'ETSC2') pm
        x = abs(dxi)
        if (x < 0.5) then
          weight = (2.75 - x**2) / 9.0
        elseif (x < 2.5) then
          weight = (3.0 - x) / 9.0
        elseif (x < 3.5) then
          weight = (3.5 - x)**2 / 18.0
        else
          weight = 0.0
        endif
!
!  Unknown weight function; give insensible value.
!
      case default pm
!
        weight = huge(1.0)
!
      endselect pm
!
    endfunction particlemesh_weighting
!***********************************************************************
    elemental real function weigh_particle(dxi1, dxi2, dxi3) result(weight)
!
!  Weighs a particle according to the particle-mesh method.
!
!  03-feb-15/ccyang: coded.
!
      real, intent(in) :: dxi1, dxi2, dxi3
!
!  Apply the weighting function in each direction.
!
      if (nxgrid > 1) then
        weight = particlemesh_weighting(dxi1)
      else
        weight = 1.0
      endif
!
      if (nygrid > 1) weight = weight * particlemesh_weighting(dxi2)
!
      if (nzgrid > 1) weight = weight * particlemesh_weighting(dxi3)
!
    endfunction weigh_particle
!***********************************************************************
!***********************************************************************
!  DUMMY PROCEDURES GO BELOW HERE.
!***********************************************************************
!***********************************************************************
    subroutine boundcond_neighbour_list()
!
!  26-apr-16/ccyang: dummy.
!
      call fatal_error('boundcond_neighbour_list', 'not implemented. ')
!
    endsubroutine boundcond_neighbour_list
!***********************************************************************
    subroutine cleanup_interpolated_quantities()
!
!  27-apr-16/ccyang: dummy.
!
      if (interp%luu .or. interp%loo .or. interp%lTT .or. interp%lrho .or. interp%lgradTT .or. interp%lbb .or. interp%lee .or. &
          interp%lpp .or. interp%lspecies .or. interp%lnu) &
          call fatal_error('cleanup_interpolated_quantities', 'interp is not implemented. ')
!
    endsubroutine cleanup_interpolated_quantities
!***********************************************************************
    subroutine fill_bricks_with_blocks(a, ab, marray, ivar)
!
!  28-apr-16/ccyang: dummy.
!
      integer, intent(in) :: marray, ivar
      real, dimension(mxb,myb,mzb,marray,0:nblockmax-1), intent(in) :: ab
      real, dimension(mx,my,mz,marray), intent(in) :: a
!
      call keep_compiler_quiet(a)
      call keep_compiler_quiet(ab(:,:,:,:,0))
      call keep_compiler_quiet(marray)
      call keep_compiler_quiet(ivar)
!
      call fatal_error('fill_bricks_with_blocks', 'only implemented for block domain decomposition. ')
!
    endsubroutine fill_bricks_with_blocks
!***********************************************************************
    subroutine fill_blocks_with_bricks(a, ab, marray, ivar)
!
!  28-apr-16/ccyang: dummy.
!
      integer, intent(in) :: marray, ivar
      real, dimension(mxb,myb,mzb,marray,0:nblockmax-1), intent(in) :: ab
      real, dimension(mx,my,mz,marray), intent(in) :: a
!
      call keep_compiler_quiet(a)
      call keep_compiler_quiet(ab(:,:,:,:,0))
      call keep_compiler_quiet(marray)
      call keep_compiler_quiet(ivar)
!
      call fatal_error('fill_blocks_with_bricks', 'only implemented for block domain decomposition. ')
!
    endsubroutine fill_blocks_with_bricks
!***********************************************************************
    subroutine interpolate_linear_range(f, ivar1, ivar2, xxp, gp, inear, iblock, ipar)
!
!  29-mar-16/ccyang: dummy.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, dimension(3), intent(in) :: inear
      real, dimension(3), intent(in) :: xxp
      integer, intent(in) :: ivar1, ivar2, iblock, ipar
      real, dimension(ivar2-ivar1+1), intent(out) :: gp
!
      gp = 0.0
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ivar1)
      call keep_compiler_quiet(ivar2)
      call keep_compiler_quiet(xxp)
      call keep_compiler_quiet(inear)
      call keep_compiler_quiet(iblock)
      call keep_compiler_quiet(ipar)
!
      call fatal_error_local('interpolate_linear_range', 'not implemented. ')
!
    endsubroutine interpolate_linear_range
!***********************************************************************
    subroutine interpolate_linear_scalar(f, ivar, xxp, gp, inear, iblock, ipar)
!
!  29-mar-16/ccyang: dummy.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, dimension(3), intent(in) :: inear
      real, dimension(3), intent(in) :: xxp
      integer, intent(in) :: ivar, iblock, ipar
      real, intent(out) :: gp
!
      gp = 0.0
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ivar)
      call keep_compiler_quiet(xxp)
      call keep_compiler_quiet(gp)
      call keep_compiler_quiet(inear)
      call keep_compiler_quiet(iblock)
      call keep_compiler_quiet(ipar)
!
      call fatal_error_local('interpolate_linear_scalar', 'not implemented. ')
!
    endsubroutine interpolate_linear_scalar
!***********************************************************************
    subroutine interpolate_quadratic(f, ivar1, ivar2, xxp, gp, inear, iblock, ipar)
!
!  30-mar-16/ccyang: dummy.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, dimension(3), intent(in) :: inear
      real, dimension(3), intent(in) :: xxp
      integer, intent(in) :: ivar1, ivar2, iblock, ipar
      real, dimension(ivar2-ivar1+1), intent(in) :: gp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ivar1)
      call keep_compiler_quiet(ivar2)
      call keep_compiler_quiet(xxp)
      call keep_compiler_quiet(gp)
      call keep_compiler_quiet(inear)
      call keep_compiler_quiet(iblock)
      call keep_compiler_quiet(ipar)
!
      call fatal_error_local('interpolate_quadratic', 'not implemented. ')
!
    endsubroutine interpolate_quadratic
!***********************************************************************
    subroutine interpolate_quadratic_spline(f, ivar1, ivar2, xxp, gp, inear, iblock, ipar)
!
!  30-mar-16/ccyang: dummy.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, dimension(3), intent(in) :: inear
      real, dimension(3), intent(in) :: xxp
      integer, intent(in) :: ivar1, ivar2, iblock, ipar
      real, dimension(ivar2-ivar1+1), intent(in) :: gp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ivar1)
      call keep_compiler_quiet(ivar2)
      call keep_compiler_quiet(xxp)
      call keep_compiler_quiet(gp)
      call keep_compiler_quiet(inear)
      call keep_compiler_quiet(iblock)
      call keep_compiler_quiet(ipar)
!
      call fatal_error_local('interpolate_quadratic_spline', 'not implemented. ')
!
    endsubroutine interpolate_quadratic_spline
!***********************************************************************
    subroutine interpolate_quantities(f, fp, p, ineargrid)
!
!  27-apr-16/ccyang: dummy.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
      type(pencil_case), intent(in) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
      call keep_compiler_quiet(p)
!
      if (interp%luu .or. interp%loo .or. interp%lTT .or. interp%lrho .or. interp%lgradTT .or. interp%lbb .or. interp%lee .or. &
          interp%lpp .or. interp%lspecies .or. interp%lnu) &
          call fatal_error('interpolate_quantities', 'interp is not implemented. ')
!
    endsubroutine interpolate_quantities
!***********************************************************************
    subroutine interpolation_consistency_check()
!
!  27-apr-16/ccyang: dummy.
!
      if (interp%luu .or. interp%loo .or. interp%lTT .or. interp%lrho .or. interp%lgradTT .or. interp%lbb .or. interp%lee .or. &
          interp%lpp .or. interp%lspecies .or. interp%lnu) &
          call fatal_error('interpolation_consistency_check', 'interp is not implemented. ')
!
    endsubroutine interpolation_consistency_check
!***********************************************************************
    subroutine map_vvp_grid(f, fp, ineargrid)
!
!  30-mar-16/ccyang: dummy.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine map_vvp_grid
!***********************************************************************
    subroutine shepherd_neighbour_block(fp, ineargrid, kshepherdb, kneighbour, iblock)
!
!  26-apr-16/ccyang: dummy.
!
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
      integer, dimension(nxb,nyb,nzb), intent(in) :: kshepherdb
      integer, dimension(:), intent(in) :: kneighbour
      integer, intent(in) :: iblock
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
      call keep_compiler_quiet(kshepherdb)
      call keep_compiler_quiet(kneighbour)
      call keep_compiler_quiet(iblock)
!
      call fatal_error('shepherd_neighbour_block', 'not implemented. ')
!
    endsubroutine shepherd_neighbour_block
!***********************************************************************
    subroutine shepherd_neighbour_pencil(fp, ineargrid, kshepherd, kneighbour)
!
!  26-apr-16/ccyang: dummy.
!
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
      integer, dimension(nx), intent(in) :: kshepherd
      integer, dimension(:), intent(in) :: kneighbour
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
      call keep_compiler_quiet(kshepherd)
      call keep_compiler_quiet(kneighbour)
!
      call fatal_error('shepherd_neighbour_pencil', 'not implemented. ')
!
    endsubroutine shepherd_neighbour_pencil
!***********************************************************************
    subroutine shepherd_neighbour_pencil3d(fp, ineargrid_c, kshepherd_c, kneighbour_c)
!
!  26-apr-16/ccyang: dummy.
!
      real, dimension(mpar_loc,mparray) :: fp
      integer, dimension(:,:,:) :: kshepherd_c
      integer, dimension(mpar_loc,3) :: ineargrid_c
      integer, dimension(mpar_loc) :: kneighbour_c
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(kshepherd_c)
      call keep_compiler_quiet(ineargrid_c)
      call keep_compiler_quiet(kneighbour_c)
!
      call fatal_error('shepherd_neighbour_pencil3d', 'not implemented. ')
!
    endsubroutine shepherd_neighbour_pencil3d
!***********************************************************************
    subroutine sort_particles_iblock(fp, ineargrid, ipar, dfp)
!
!  26-apr-16/ccyang: dummy.
!
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
      integer, dimension(mpar_loc), intent(in) :: ipar
      real, dimension(mpar_loc,mpvar), intent(in), optional :: dfp
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
      call keep_compiler_quiet(ipar)
      if (present(dfp)) call keep_compiler_quiet(dfp)
!
      call fatal_error('sort_particles_iblock', 'only implemented for block domain decomposition. ')
!
    endsubroutine sort_particles_iblock
!***********************************************************************
endmodule Particles_map
