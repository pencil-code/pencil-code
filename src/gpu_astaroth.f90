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
      !if (lviscosity) str='viscosity'
      if (lpower_spectrum) str='power_spectrum'
      if (lparticles) str='particles'

      if (str/='') &
        call stop_it('No GPU implementation for module "'//trim(str)//'"')
!
      call initialize_gpu_c

    endsubroutine initialize_GPU
!**************************************************************************
    subroutine finalize_GPU
!
      call finalize_gpu_c
!
    endsubroutine finalize_GPU
!**************************************************************************
    subroutine rhs_GPU(f,isubstep)
!
      real, dimension (mx,my,mz,mfarray), intent(INOUT) :: f
      integer,                            intent(IN)    :: isubstep
!
      integer :: ll, mm, nn
      real :: val
      logical, save :: lvery_first=.true.
      logical :: lfull_inner

      goto 1
      val=1.
      do nn=1,mz
        do mm=1,my
          do ll=1,mx
            f(ll,mm,nn,iux)=val; val=val+1.
      enddo; enddo; enddo
      f(1,1,1,iuy)=-1.; f(1,1,1,iuz)=-1.; f(1,1,1,ilnrho)=-1.

1     lfull_inner = (lsnap.or.lsnap_down.or.lspec) .and. isubstep==itorder

      call rhs_gpu_c(f(1,1,1,iux),f(1,1,1,iuy),f(1,1,1,iuz),f(1,1,1,ilnrho), &
                     isubstep,lfull_inner,lvery_first)
!
      lvery_first=.false.

      return
      if (.not.lroot) return
      do nn=1,mz   !  nghost+1,mz-nghost   !1,mz
        print*, 'nn=', nn
        do mm=1,my
          print'(14(1x,f7.0))',f(:,mm,nn,iux)
      enddo; enddo

    endsubroutine rhs_GPU
!**************************************************************************
endmodule GPU
