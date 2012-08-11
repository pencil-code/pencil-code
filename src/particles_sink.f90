! $Id: particles_dust.f90 19206 2012-06-30 21:40:24Z sven.bingert $
!
!  This module takes care of everything related to sink particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_sink=.true.
!
! MPVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
! COMMUNICATED AUXILIARIES 1
!
!***************************************************************
module Particles_sink
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
  use Particles_sub
!
  implicit none
!
  include 'particles_sink.h'
!
  real :: sink_radius=0.0, rhop_sink_create=1.0
  logical :: lsink_radius_dx_unit=.false., lrhop_roche_unit=.false.
!
  integer :: idiag_nparsink=0
!
  namelist /particles_sink_init_pars/ &
      sink_radius, lsink_radius_dx_unit, rhop_sink_create, lrhop_roche_unit
!
  namelist /particles_sink_run_pars/ &
      sink_radius, lsink_radius_dx_unit, rhop_sink_create, lrhop_roche_unit
!
  contains
!***********************************************************************
    subroutine register_particles_sink()
!
!  Set up indices for access to the fp and dfp arrays
!
!  07-aug-12/anders: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
           "$Id: particles_dust.f90 19206 2012-06-30 21:40:24Z sven.bingert $")
!
!  Indices for sink particle radius.
!
      israd=npvar+1
      pvarname(npvar+1)='israd'
      npvar=npvar+1
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_sink_particles','npvar > mpvar')
      endif
!
    endsubroutine register_particles_sink
!***********************************************************************
    subroutine initialize_particles_sink(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  07-aug-12/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      if (lsink_radius_dx_unit) sink_radius=sink_radius*dx
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_sink
!***********************************************************************
    subroutine init_particles_sink(f,fp)
!
!  Initial sink particle radii.
!
!  07-aug-12/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
!
!  Set all particles to have zero sink radius.
!
      fp(1:npar_loc,israd)=0.0
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_particles_sink
!***********************************************************************
    subroutine create_particles_sink(f,fp,dfp,ineargrid)
!
!  Create sink particles based on local particle density.
!
!  07-aug-12/anders: coded
!
      use Boundcond
      use Diagnostics
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mpvar) :: fp, dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      real, dimension(1) :: rhop_interp
      integer :: k, ix0, iy0, iz0, npar_sink
!
      if (lparticles_blocks) then

      else
        do imn=1,ny*nz
          if (npar_imn(imn)/=0) then
            do k=k1_imn(imn),k2_imn(imn)
              if (fp(k,israd)==0.0) then
                ix0=ineargrid(k,1)
                iy0=ineargrid(k,2)
                iz0=ineargrid(k,3)
                if (lparticlemesh_cic) then
                  call interpolate_linear(f,irhop,irhop, &
                      fp(k,ixp:izp),rhop_interp,ineargrid(k,:),0,ipar(k))
                elseif (lparticlemesh_tsc) then
                  if (linterpolate_spline) then
                    call interpolate_quadratic_spline(f,irhop,irhop, &
                        fp(k,ixp:izp),rhop_interp,ineargrid(k,:),0,ipar(k))
                  else
                    call interpolate_quadratic(f,irhop,irhop, &
                        fp(k,ixp:izp),rhop_interp,ineargrid(k,:),0,ipar(k))
                  endif
                else
                  rhop_interp=f(ix0,iy0,iz0,irhop:irhop)
                endif
                if (rhop_interp(1)>=rhop_sink_create) then
                  fp(k,israd)=sink_radius
                endif
              endif
            enddo
          endif
        enddo
      endif
!
      if (ldiagnos) then
        if (idiag_nparsink/=0) then
          npar_sink=0
          do k=1,npar_loc
            if (fp(k,israd)/=0.0) npar_sink=npar_sink+1
          enddo
          call save_name(float(npar_sink),idiag_nparsink)
        endif
      endif
!
    endsubroutine create_particles_sink
!***********************************************************************
    subroutine remove_particles_sink(f,fp,dfp,ineargrid)
!
!  Remove particles in the vicinity of sink particles.
!
!  07-aug-12/anders: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mpvar) :: fp, dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      real, dimension(3) :: xxps, vvps
      real :: rhops, rads, rads2, dist2
      integer, dimension(3) ::  dis=(/-1,0,+1/)
      integer :: j, k, jsink

!
      if (lmpicomm) then

      else
!
!  Loop over sink particles, starting at npar_loc so that the sink particle
!  position in fp does not change when particles are removed.
!
        j=npar_loc
        do while (j>=1)
          if (fp(j,israd)/=0.0) then
!
!  Store sink particle information in separate variables, as this might give
!  better cache efficiency.
!
            xxps=fp(j,ixp:izp)
            vvps=fp(j,ivpx:ivpz)
            rhops=fp(j,irhopswarm)
            rads=fp(j,israd)
            rads2=fp(j,israd)**2
            k=npar_loc
            jsink=j
!
!  Loop over all particles to see which ones get removed by the sink particle.
!
            do while (k>=1)
              if (j/=k) then
!
!  Find minimum distance by directional splitting. This makes it easier to
!  incorporate periodic boundary conditions.
!
                if (minval(abs(fp(k,ixp)-(xxps(1)+Lx*dis)))<=rads) then
                  if (minval(abs(fp(k,iyp)-(xxps(2)+Ly*dis)))<=rads) then
                    if (minval(abs(fp(k,izp)-(xxps(3)+Lz*dis)))<=rads) then
!
!  Particle is constrained to be within cube of size rads. Estimate whether
!  the particle is also within the sphere of size rads.
!
                      dist2=minval(((fp(k,ixp)-(xxps(1)+Lx*dis))**2+ &
                          (fp(k,iyp)-(xxps(2)+Ly*dis))**2+ &
                          (fp(k,izp)-(xxps(3)+Lz*dis))**2))
!
!  Do not allow sink particle to accrete itself. This is already excluded
!  above by demanding j/=k, but in pathological cases where j=npar_loc, the
!  sink particle may be moved down to replace a removed particle.
!
                      if (dist2<=rads2) then
                        if (lparticles_mass) then
                          vvps=(rhops*vvps+fp(k,irhopswarm)* &
                              fp(k,ivpx:ivpz))/(rhops+fp(k,irhopswarm))
                          rhops=rhops+fp(k,irhopswarm)
                        endif
!
!  If the sink particle has index npar_loc, then it will be moved to k.
!
                        if (j==npar_loc) then
                          jsink=k
                        endif
!
!  Remove particle, replacing it in fp with the particle of index npar_loc.
!
                        call remove_particle(fp,ipar,k,dfp,ineargrid)
                      endif
                    endif
                  endif
                endif
              endif
              k=k-1
            enddo
!
!  Give accreted mass and velocity to the sink particle. This conserves mass
!  and momentum, but angular momentum is not conserved - this is assumed to
!  be stored in rotation of the sink particle.
!
            fp(jsink,ivpx:ivpz)=vvps
            fp(jsink,irhopswarm)=rhops
          endif
          j=j-1
        enddo
      endif
!
    endsubroutine remove_particles_sink
!***********************************************************************
    subroutine read_particles_sink_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_sink_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_sink_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_sink_init_pars
!***********************************************************************
    subroutine write_particles_sink_init_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_sink_init_pars)
!
    endsubroutine write_particles_sink_init_pars
!***********************************************************************
    subroutine read_particles_sink_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_sink_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_sink_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_sink_run_pars
!***********************************************************************
    subroutine write_particles_sink_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_sink_run_pars)
!
    endsubroutine write_particles_sink_run_pars
!***********************************************************************
    subroutine rprint_particles_sink(lreset,lwrite)
!
!  Read and register print parameters relevant for particles sink radius.
!
!  11-aug-12/anders: coded
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip<14) &
          print*, 'rprint_particles_sink: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'nparsink', &
            idiag_nparsink)
      enddo
!
    endsubroutine rprint_particles_sink
!***********************************************************************
endmodule Particles_sink
