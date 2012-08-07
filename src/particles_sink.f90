! $Id: particles_dust.f90 19206 2012-06-30 21:40:24Z sven.bingert $
!
!  This module takes care of everything related to sink particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
! CPARAM logical, parameter :: lparticles_sink=.true.
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
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mpvar) :: fp, dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      real, dimension(1) :: rhop_interp
      integer :: k, ix0, iy0, iz0
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
                  if (rhop_interp(1)>=rhop_sink_create) then
                    fp(k,israd)=sink_radius
                  endif
                endif
              endif
            enddo
          endif
        enddo
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
    endsubroutine remove_particles_sink
!***********************************************************************
endmodule Particles_sink
