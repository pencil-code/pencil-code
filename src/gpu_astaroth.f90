! $Id$
!
! MODULE_DOC: This module contains GPU related types and functions to be used with the ASTAROTH nucleus.
!
! CPARAM logical, parameter :: lgpu = .true.
!
!**************************************************************************
!
module GPU
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Mpicomm, only: stop_it

  implicit none

  external initialize_gpu_c
  external finalize_gpu_c
  external rhs_gpu_c

  include 'gpu.h'

contains

!***********************************************************************
    subroutine initialize_GPU
!
      character(LEN=30) :: str
!
      str=''
      if (lanelastic) str='anelastic'
      if (lboussinesq) str='boussinesq'
      if (lenergy) str='energy'
      if (lentropy) str='entropy'
      if (ltemperature) str='temperature'
      if (lshock) str='shock'
      if (lmagnetic) str='magnetic'
      if (lforcing) str='forcing'
      if (lgrav) str='gravity'
      if (lheatflux) str='heatflux'
      if (lhyperresistivity_strict) str='hyperresi_strict'
      if (lhyperviscosity_strict) str='hypervisc_strict'
      if (lADI) str='implicit_physics'
      if (llorenz_gauge) str='lorenz_gauge'
      if (ldustvelocity) str='dustvelocity'
      if (ldustdensity) str='dustdensity'
      if (ltestscalar) str='testscalar'
      if (ltestfield) str='testfield'
      if (ltestflow) str='testflow'
      if (linterstellar) str='interstellar'
      if (lcosmicray) str='cosmicray'
      if (lcosmicrayflux) str='cosmicrayflux'
      if (lshear) str='shear'
      if (lpscalar) str='pscalar'
      if (lsupersat) str='supersat'
      if (lradiation) str='radiation'
      if (leos) str='eos'
      if (lchemistry) str='chemistry'
      if (lchiral) str='chiral'
      if (ldetonate) str='detonate'
      if (lneutralvelocity) str='neutralvelocity'
      if (lneutraldensity) str='neutraldensity'
      if (lopacity) str='opacity'
      if (lpolymer) str='polymer'
      if (lpointmasses) str='pointmasses'
      if (lpoisson) str='poisson'
      if (lselfgravity) str='selfgravity'
      if (lsolid_cells) str='solid_cells'
      if (lspecial) str='special'
      if (lviscosity) str='viscosity'
      if (lpower_spectrum) str='power_spectrum'
      if (lparticles) str='particles'

      if (str/='') &
        call stop_it('No GPU implementation for module "'//trim(str)//'"')
!
      call initialize_gpu_c(0)

    endsubroutine initialize_GPU
!**************************************************************************
    subroutine finalize_GPU
!
      call finalize_gpu_c(0)
!
    endsubroutine finalize_GPU
!**************************************************************************
    subroutine rhs_GPU(f,df)
!
      real, dimension (mx,my,mz,mfarray), intent(INOUT) :: f
      real, dimension (mx,my,mz,mvar) :: df
!
      call rhs_gpu_c(f(1,1,1,iuu),f(1,1,1,ilnrho),df(1,1,1,iuu),df(1,1,1,ilnrho))
!
    endsubroutine rhs_GPU
!**************************************************************************
endmodule  GPU
