! $Id$
!
!  This modules takes care of viscosity of inertial particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_viscosity = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 6
! COMMUNICATED AUXILIARIES 3
!
!***************************************************************
module Particles_viscosity
!
  use Cparam
  use Cdata
  use Messages
  use Particles_cdata
  use Particles_sub
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_viscosity.h'
!
  integer, parameter :: nviscp_max = 4
  character (len=labellen), dimension(nviscp_max) :: iviscp=''
  real :: nup=0.0
  logical :: lviscp_simplified=.false.
  logical :: lviscp_rhop_nup_const=.false.
  logical :: lviscp_nup_const=.false.
!
  namelist /particles_visc_run_pars/ &
      nup, iviscp
!
  contains
!***********************************************************************
    subroutine register_particles_viscosity()
!
!  07-oct-08/anders: coded
!
      use FArrayManager
      use Sub
!
!  Identify version number.
!
      if (lroot) call cvs_id( &
          '$Id$')
!
!  Set indices for auxiliary variables.
!
      call farray_register_auxiliary('upviscx',ipviscx,vector=3)
      ipviscy=ipviscx+1; ipviscz=ipviscx+2
      call farray_register_auxiliary('upx',iupx,communicated=.true.,vector=3)
      iupy=iupx+1; iupz=iupx+2
!
    endsubroutine register_particles_viscosity
!***********************************************************************
    subroutine initialize_particles_viscosity(lstarting)
!
!  07-oct-08/anders: coded
!
      use FArrayManager
      use Mpicomm, only: stop_it
      use SharedVariables
!
      logical, intent(in) :: lstarting
!
      integer :: i, ierr
!
!  Some viscosity types need the rate-of-strain tensor and grad(lnrho)
!
      lviscp_simplified=.false.
      lviscp_rhop_nup_const=.false.
      lviscp_nup_const=.false.
!
      do i=1,nviscp_max
        select case (iviscp(i))
        case ('simplified', '0')
          if (lroot) print*,'particle viscosity: nup*del2v'
          lviscp_simplified=.true.
        case('rhop_nup-const', '1')
          if (lroot) print*,'particle viscosity: mup/rhop*(del2v+graddivv/3)'
          lviscp_rhop_nup_const=.true.
        case('nup-const')
          if (lroot) &
              print*,'particle viscosity: nup*(del2v+graddivv/3+2S.glnrhop)'
          lviscp_nup_const=.true.
        case ('none','')
          ! do nothing
        case default
          if (lroot) print*, 'No such value for iviscp(',i,'): ', trim(iviscp(i))
          call fatal_error('initialize_particles_viscosity','')
        endselect
      enddo
!
!  If we're timestepping, die or warn if the viscosity coefficient that
!  corresponds to the chosen viscosity type is not set.
!
      if (lrun) then
        if ( (lviscp_simplified.or.lviscp_rhop_nup_const.or.lviscp_nup_const) &
            .and.nup==0.0) &
            call warning('initialize_particles_viscosity', &
            'Viscosity coefficient nup is zero!')
      endif
!
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_viscosity
!***********************************************************************
    subroutine calc_particles_viscosity(f,fp,ineargrid)
!
!  Calculate the Laplacian of the particle velocity field, for use in
!  particle viscosity later.
!
!  07-feb-09/anders: coded
!
      use Boundcond
      use Mpicomm
      use Particles_sub, only: map_vvp_grid
!
      real, dimension (mx,my,mz,mfarray) :: f 
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (nx,3) :: del2uup
      real, save :: dx_2, dy_2, dz_2
      logical, save :: lfirstcall=.true.
!
      if (lfirstcall) then
        dx_2=1/dx**2
        dy_2=1/dy**2
        dz_2=1/dz**2
      endif
!
!  Map the particle velocities as a vector field on the grid.
!
      call map_vvp_grid(f,fp,ineargrid)
!
!  Put boundary conditions on mapped velocity field.
!
      call boundconds_x(f,iupx,iupz)
      call initiate_isendrcv_bdry(f,iupx,iupz)
      call finalize_isendrcv_bdry(f,iupx,iupz)
      call boundconds_y(f,iupx,iupz)
      call boundconds_z(f,iupx,iupz)
!
!  We use a second order discrete Laplacian.
!
      if (lviscp_simplified .or. lviscp_rhop_nup_const) then
        do n=n1,n2; do m=m1,m2
          if (nxgrid/=1) f(l1:l2,m,n,ipviscx:ipviscz) = &
              (f(l1+1:l2+1,m,n,iupx:iupz) - 2*f(l1:l2,m,n,iupx:iupz) + &
               f(l1-1:l2-1,m,n,iupx:iupz))*dx_2
          if (nygrid/=1) then
            if (nxgrid==1) then
              f(l1:l2,m,n,ipviscx:ipviscz) = &
                  (f(l1:l2,m+1,n,iupx:iupz) - 2*f(l1:l2,m,n,iupx:iupz) + &
                   f(l1:l2,m-1,n,iupx:iupz))*dy_2
            else
              f(l1:l2,m,n,ipviscx:ipviscz) = f(l1:l2,m,n,ipviscx:ipviscz) + &
                  (f(l1:l2,m+1,n,iupx:iupz) - 2*f(l1:l2,m,n,iupx:iupz) + &
                   f(l1:l2,m-1,n,iupx:iupz))*dy_2
            endif
          endif
          if (nzgrid/=1) then
            if (nxgrid==1 .and. nygrid==1) then
              f(l1:l2,m,n,ipviscx:ipviscz) = &
                  (f(l1:l2,m,n+1,iupx:iupz) - 2*f(l1:l2,m,n,iupx:iupz) + &
                   f(l1:l2,m,n-1,iupx:iupz))*dz_2
            else
              f(l1:l2,m,n,ipviscx:ipviscz) = f(l1:l2,m,n,ipviscx:ipviscz) + &
                  (f(l1:l2,m,n+1,iupx:iupz) - 2*f(l1:l2,m,n,iupx:iupz) + &
                   f(l1:l2,m,n-1,iupx:iupz))*dz_2
            endif
          endif
        enddo; enddo
      endif
!
      lfirstcall=.false.
!
    endsubroutine calc_particles_viscosity
!***********************************************************************
    subroutine dvvp_dt_viscosity_pencil(f,df,fp,dfp,ineargrid)
!
!  Calculate viscous force term for particles.
!
!  07-feb-09/anders: coded
!
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f 
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: f, df, fp, ineargrid
      intent (out) :: dfp
!
      real, dimension (3) :: fviscp
      integer :: k, ix0, iy0, iz0
!
      if (npar_imn(imn)/=0) then
!
!  Loop over all particles in current pencil.
!
        do k=k1_imn(imn),k2_imn(imn)
          ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
          if (lparticlemesh_cic) then
            call interpolate_linear(f,ipviscx,ipviscz, &
                fp(k,ixp:izp),fviscp,ineargrid(k,:),ipar(k) )
          elseif (lparticlemesh_tsc) then
            if (linterpolate_spline) then
              call interpolate_quadratic_spline(f,ipviscx,ipviscz, &
                  fp(k,ixp:izp),fviscp,ineargrid(k,:),ipar(k) )
            else
              call interpolate_quadratic(f,ipviscx,ipviscz, &
                  fp(k,ixp:izp),fviscp,ineargrid(k,:),ipar(k) )
            endif
          else
            fviscp=f(ix0,iy0,iz0,ipviscx:ipviscz)
          endif
!
!  Viscous force dvp/dt = nu*del2(uup). Does not conserve momentum, but
!  should be okay if the particle density is nearly constant.
!
          if (lviscp_simplified) then
            dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + nup*fviscp
          endif
!
!  Viscous force dvp/dt = (mu/rho)*del2(uup). Conserves momentum and is
!  a reasonable description of viscosity since nup=cs*lambda yields
!  nup*np=constant.
!
          if (lviscp_rhop_nup_const) then
            dfp(k,ivpx:ivpz)=dfp(k,ivpx:ivpz)+nup*fviscp/f(ix0,iy0,iz0,irhop)
          endif
!
        enddo
!
      endif
!
!  Calculate viscosity time-step.
!
      if (lfirst.and.ldt) dt1_max=max(dt1_max,0.25*dxyz_2*nup)
!
    endsubroutine dvvp_dt_viscosity_pencil
!***********************************************************************
    subroutine read_particles_visc_init_pars(unit,iostat)
!
!
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit
!
    endsubroutine read_particles_visc_init_pars
!***********************************************************************
    subroutine write_particles_visc_init_pars(unit)
!
!
!
      integer, intent(in) :: unit
!
      if (NO_WARN) print*,unit
!
    endsubroutine write_particles_visc_init_pars
!***********************************************************************
    subroutine read_particles_visc_run_pars(unit,iostat)
!
!
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_visc_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_visc_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_particles_visc_run_pars
!***********************************************************************
    subroutine write_particles_visc_run_pars(unit)
!
!
!
      integer, intent(in) :: unit
!
      write(unit,NML=particles_visc_run_pars)
!
    endsubroutine write_particles_visc_run_pars
!***********************************************************************
    subroutine calc_viscosity(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (NO_WARN) print*,f  !(to keep compiler quiet)
!
    endsubroutine calc_viscosity
!*******************************************************************
    subroutine rprint_particles_viscosity(lreset,lwrite)
!
!  Writes ishock to index.pro file
!
!  07-oct-08/anders: adapted
!
      use Sub
!
      logical :: lreset
      logical, optional :: lwrite
      integer :: iname
!
      if (present(lwrite)) then
        if (lwrite) then
          write(3,*) 'ipviscx=', ipviscx
          write(3,*) 'ipviscy=', ipviscy
          write(3,*) 'ipviscz=', ipviscz
          write(3,*) 'iupx=', iupx
          write(3,*) 'iupy=', iupy
          write(3,*) 'iupz=', iupz
        endif
      endif
!
      if (NO_WARN) print*,lreset  !(to keep compiler quiet)
!
    endsubroutine rprint_particles_viscosity
!***********************************************************************
endmodule Particles_viscosity
