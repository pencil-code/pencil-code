! $Id$
!
! MODULE_DOC: Runge-Kutta time advance, accurate to order itorder.
! MODULE_DOC: At the moment, itorder can be 1, 2, or 3.
!
module Timestep
!
  use Cparam
  use Cdata
!
  implicit none
!
  private
!
  public :: time_step
!
  contains
!***********************************************************************
    subroutine time_step(f,df,p)
!
!   2-apr-01/axel: coded
!  14-sep-01/axel: moved itorder to cdata
!
      use BorderProfiles, only: border_quenching
      use Equ, only: pde
      use Interstellar, only: calc_snr_damp_int
      use Mpicomm
      use Particles_main, only: particles_timestep_first, &
          particles_timestep_second
      use Shear, only: advance_shear
      use Snapshot, only: shift_dt
      use Entropy
      use Special
      use Boundcond
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: ds
      real :: dt1, dt1_local, dt1_last=0.0
      integer :: j,ienergy, Nsub !, N_max
      real, dimension (mx,my,mz,mfarray) :: ftmp
      real, dimension (nx) :: dt1_hcond_max
      real :: dt1_energy_local,dt1_energy,dt_energy
 !
      if (ltemperature) then
        if (ltemperature_nolog) then
          ienergy = iTT
        else
          ienergy = ilnTT
        endif
      endif
      if (lentropy) then
        if (pretend_lnTT) then
          ienergy = ilnTT
        else
          ienergy = iss
        endif
      endif
      if (lthermal_energy) ienergy=ieth
!
!  Coefficients for up to order 3.
!
      if (itorder==1) then
        alpha_ts=(/ 0.0, 0.0, 0.0 /)
        beta_ts =(/ 1.0, 0.0, 0.0 /)
      elseif (itorder==2) then
        alpha_ts=(/   0.0, -1/2.0, 0.0 /)
        beta_ts =(/ 1/2.0,    1.0, 0.0 /)
      elseif (itorder==3) then
        !alpha_ts=(/   0.0, -2/3.0, -1.0   /)
        !beta_ts =(/ 1/3.0,    1.0,  1/2.0 /)
        !  use coefficients of Williamson (1980)
        alpha_ts=(/   0.0, -5/9.0 , -153/128.0 /)
        beta_ts =(/ 1/3.0, 15/16.0,    8/15.0  /)
      else
        if (lroot) print*,'Not implemented: itorder=',itorder
        call mpifinalize
      endif
!
!  dt_beta_ts may be needed in other modules (like Dustdensity) for fixed dt.
!
      if (.not. ldt) dt_beta_ts=dt*beta_ts
!
!  Set up df and ds for each time sub.
!
      do itsub=1,itorder
        lfirst=(itsub==1)
        llast=(itsub==itorder)
        if (lfirst) then
          df=0.0
          ds=0.0
        else
          df=alpha_ts(itsub)*df !(could be subsumed into pde, but is dangerous!)
          ds=alpha_ts(itsub)*ds
        endif
!
!  Change df according to the chosen physics modules.
!
        call pde(f,df,p)
!
        ds=ds+1.0
!
!  If we are in the first time substep we need to calculate timestep dt.
!  Only do it on the root processor, then broadcast dt to all others.
!
        if (lfirst.and.ldt) then
          dt1_local=maxval(dt1_max(1:nx))
          ! Timestep growth limiter
          if (real(ddt) > 0.) dt1_local=max(dt1_local,dt1_last)
!
!  Get time step for heat conduction in sub step, when sub_step_hcond
!  in coronae module is F, dt1_energy will be 0, and sub iteration below
!  is turned off. This allows turn on/off sub iteration without recompile.
          call calc_hcond_timestep(f,p,dt1_hcond_max)
          dt1_energy_local=maxval(dt1_hcond_max(1:nx))
!
          ! Get global time step limite
          call mpiallreduce_max(dt1_local,dt1)
          call mpiallreduce_max(dt1_energy_local,dt1_energy)
          dt=1.0/dt1
!  calc the number of sub iteration
          Nsub = int(dt1_energy/dt1)+1
! in pdf(f,df,p) ghost cells of f-array are set
          ftmp = f
          if (loutput_varn_at_exact_tsnap) call shift_dt(dt)
          ! Timestep growth limiter
          if (ddt/=0.) dt1_last=dt1_local/ddt
        endif
!
!  Calculate dt_beta_ts.
!
        if (ldt) dt_beta_ts=dt*beta_ts
        if (ip<=6) print*, 'time_step: iproc, dt=', iproc, dt  !(all have same dt?)
!
!  Time evolution of grid variables.
!  (do this loop in pencils, for cache efficiency)
!
        do j=1,mvar; do n=n1,n2; do m=m1,m2
!ajwm Note to self... Just how much overhead is there in calling
!ajwm a sub this often...
            if (lborder_profiles) call border_quenching(f,df,j,dt_beta_ts(itsub))
            f(l1:l2,m,n,j)=f(l1:l2,m,n,j)+dt_beta_ts(itsub)*df(l1:l2,m,n,j)
        enddo; enddo; enddo
!
        if (lspecial) call special_after_timestep(f,df,dt_beta_ts(itsub)*ds)
!
!  Increase time.
!
        t = t + dt_beta_ts(itsub)*ds
!
      enddo
!!
!!  Now do the substeps for the energy (thermal conduction) only ----------------------------------
!!
      dt_energy = dt / dble(Nsub)
      if ((dt_energy /= dt) .and. lroot .and. (mod(it,isave) == 0)) &
      write(*,'(a,I3,a,ES10.2)') 'Nsub = ',Nsub,' | dt_energy = ',dt_energy
!
      IF (Nsub > 1) THEN ! Turn of the subcircle when dt1_energy actully = 0.
!
      do n=n1,n2; do m=m1,m2
        ftmp(l1:l2,m,n,ienergy)=f(l1:l2,m,n,ienergy)
      enddo; enddo
!
      do j=1,Nsub
        do itsub=1,itorder
          lfirst=(itsub==1)
          llast=(itsub==itorder)
          if (lfirst) then
            df=0.0
            ds=0.0
          else
            df=alpha_ts(itsub)*df !(could be subsumed into pde, but is dangerous!)
            ds=alpha_ts(itsub)*ds
          endif
!
          call pde_energy_only(ftmp,df,p)
!
          ds=ds+1.0
!
          if (ldt) dt_beta_ts=dt_energy*beta_ts
          if (ip<=6) print*, 'time_step: iproc, dt=', iproc, dt_energy  !(all have same dt?)
!!!!
!!!!  Time evolution of grid variables.
!!!!  (do this loop in pencils, for cache efficiency)
!!!!
          do n=n1,n2; do m=m1,m2
!ajwm Note to self... Just how much overhead is there in calling
!ajwm a sub this often...
            ftmp(l1:l2,m,n,ienergy)=ftmp(l1:l2,m,n,ienergy)+dt_beta_ts(itsub)*df(l1:l2,m,n,ienergy)
          enddo; enddo
!
        enddo
      enddo
!
       do n=n1,n2; do m=m1,m2
         f(l1:l2,m,n,ienergy)=ftmp(l1:l2,m,n,ienergy)
       enddo; enddo
!
       ENDIF
!
    endsubroutine time_step
!***********************************************************************
    subroutine pde_energy_only(f,df,p)
!
!  Calculate thermal conduction for sub_circle time step.
!
!  22-jun-12/feng/bingert: coded
!
      use Boundcond
      use BorderProfiles, only: calc_pencils_borderprofiles
      use Density
      use Diagnostics
      use Entropy
      use EquationOfState
      use Grid, only: calc_pencils_grid
      use Hydro
      use Lorenz_gauge
      use Magnetic
      use Mpicomm
      use Special
      use Sub
!
      logical :: early_finalize
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(inout)  :: f       ! inout due to  lshift_datacube_x,
                                ! density floor, or velocity ceiling
      intent(out)    :: df, p
!
!  Print statements when they are first executed.
!
      headtt = headt .and. lfirst .and. lroot
!
      if (headtt.or.ldebug) print*,'pde: ENTER'
!
!  For chemistry with LSODE
!
      lchemonly=.false.
!
!  Grid spacing. For non equidistant grid or non-cartesian coordinates
!  the grid spacing is calculated in the (m,n) loop below.
!
      if (lcartesian_coords .and. all(lequidist)) then
        if (old_cdtv) then
          dxyz_2 = max(dx_1(l1:l2)**2,dy_1(m1)**2,dz_1(n1)**2)
        else
          dline_1(:,1)=dx_1(l1:l2)
          dline_1(:,2)=dy_1(m1)
          dline_1(:,3)=dz_1(n1)
          dxyz_2 = dline_1(:,1)**2+dline_1(:,2)**2+dline_1(:,3)**2
          dxyz_4 = dline_1(:,1)**4+dline_1(:,2)**4+dline_1(:,3)**4
          dxyz_6 = dline_1(:,1)**6+dline_1(:,2)**6+dline_1(:,3)**6
        endif
      endif
!
!  Shift entire data cube by one grid point at the beginning of each
!  time-step. Useful for smearing out possible x-dependent numerical
!  diffusion, e.g. in a linear shear flow.
!
      if (lfirst .and. lshift_datacube_x) then
        call boundconds_x(f)
        do  n=n1,n2; do m=m1,m2
          f(:,m,n,:)=cshift(f(:,m,n,:),1,1)
        enddo; enddo
      endif
!
!  Need to finalize communication early either for test purposes, or
!  when radiation transfer of global ionization is calculated.
!  This could in principle be avoided (but it not worth it now)
!
      early_finalize=test_nonblocking.or. &
                     leos_ionization.or.lradiation_ray.or. &
                     lhyperviscosity_strict.or.lhyperresistivity_strict.or. &
                     ltestscalar.or.ltestfield.or.ltestflow.or. &
                     lparticles_prepencil_calc.or.lsolid_cells.or. &
                     lchemistry.or.lweno_transport
!
!  Write crash snapshots to the hard disc if the time-step is very low.
!  The user must have set crash_file_dtmin_factor>0.0 in &run_pars for
!  this to be done.
!
!
!  For debugging purposes impose minimum or maximum value on certain variables.
!
      call impose_energy_floor(f)
!
!  Call "before_boundary" hooks (for f array precalculation)
!
      if (lspecial)      call special_before_boundary(f)
!
!  Fetch fp to the special module
!
!  Prepare x-ghost zones; required before f-array communication
!  AND shock calculation
!
      call boundconds_x(f,ilnTT,ilnTT)
!
!  Initiate (non-blocking) communication and do boundary conditions.
!  Required order:
!  1. x-boundaries (x-ghost zones will be communicated) - done above
!  2. communication
!  3. y- and z-boundaries
!
      if (ldebug) print*,'pde: bef. initiate_isendrcv_bdry'
      call initiate_isendrcv_bdry(f,ilnTT,ilnTT)
     if (early_finalize) then
       call finalize_isendrcv_bdry(f,ilnTT,ilnTT)
       call boundconds_y(f,ilnTT,ilnTT)
       call boundconds_z(f,ilnTT,ilnTT)
     endif
!
!DM I suggest the following lhydro_pars, lmagnetic_pars be renamed to
! hydro_after_boundary etc.
!AB: yes, we should rename these step by step
!AB: so calc_polymer_after_boundary -> polymer_after_boundary
!!!      if (lhydro)                 call calc_lhydro_pars(f)
!!!      if (lmagnetic)              call calc_lmagnetic_pars(f)
!!!   if (lmagnetic)              call magnetic_after_boundary(f)
!!!      if (lentropy)               call calc_lentropy_pars(f)
!!!      if (ldensity)               call calc_ldensity_pars(f)
!!!      if (lspecial)               call calc_lspecial_pars(f)
!AB: could be renamed to special_after_boundary etc
      !if (lspecial)               call special_after_boundary(f)
!
!
!------------------------------------------------------------------------------
!  Do loop over m and n.
!
      mn_loop: do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        lfirstpoint=(imn==1)      ! true for very first m-n loop
        llastpoint=(imn==(ny*nz)) ! true for very last m-n loop
!
!  Make sure all ghost points are set.
!
       if (.not.early_finalize.and.necessary(imn)) then
         call finalize_isendrcv_bdry(f,ilnTT,ilnTT)
         call boundconds_y(f,ilnTT,ilnTT)
         call boundconds_z(f,ilnTT,ilnTT)
       endif
!
!  Grid spacing. In case of equidistant grid and cartesian coordinates
!  this is calculated before the (m,n) loop.
!
        if (lspherical_coords.or.lcylindrical_coords.or. &
            .not.all(lequidist)) then
          if (old_cdtv) then
!
!  The following is only kept for backwards compatibility. Will be deleted in
!  the future.
!
            dxyz_2 = max(dx_1(l1:l2)**2,dy_1(m)**2,dz_1(n)**2)
          else
            if (lspherical_coords) then
              dline_1(:,1)=dx_1(l1:l2)
              dline_1(:,2)=r1_mn*dy_1(m)
              dline_1(:,3)=r1_mn*sin1th(m)*dz_1(n)
            else if (lcylindrical_coords) then
              dline_1(:,1)=dx_1(l1:l2)
              dline_1(:,2)=rcyl_mn1*dy_1(m)
              dline_1(:,3)=dz_1(n)
            else if (lcartesian_coords) then
              dline_1(:,1)=dx_1(l1:l2)
              dline_1(:,2)=dy_1(m)
              dline_1(:,3)=dz_1(n)
            endif
            dxmax_pencil = 0.
            if (nxgrid /= 1) dxmax_pencil=1./dx_1(l1:l2)
            if (nygrid /= 1) dxmax_pencil=max(1./dy_1(m),dxmax_pencil)
            if (nzgrid /= 1) dxmax_pencil=max(1./dz_1(n),dxmax_pencil)
            dxmin_pencil = 0.
            if (nxgrid /= 1) dxmin_pencil=1./dx_1(l1:l2)
            if (nygrid /= 1) dxmin_pencil=min(1./dy_1(m),dxmin_pencil)
            if (nzgrid /= 1) dxmin_pencil=min(1./dz_1(n),dxmin_pencil)
!
            dxyz_2 = dline_1(:,1)**2+dline_1(:,2)**2+dline_1(:,3)**2
            dxyz_4 = dline_1(:,1)**4+dline_1(:,2)**4+dline_1(:,3)**4
            dxyz_6 = dline_1(:,1)**6+dline_1(:,2)**6+dline_1(:,3)**6
          endif
        endif
!
!  Calculate grid/geometry related pencils.
!
        call calc_pencils_grid(f,p)
!
!  Calculate pencils for the pencil_case.
!  Note: some no-modules (e.g. nohydro) also calculate some pencils,
!  so it would be wrong to check for lhydro etc in such cases.
!
! To check ghost cell consistency, please uncomment the following 2 lines:
!       if (.not. lpencil_check_at_work .and. necessary(imn)) &
!       call check_ghosts_consistency (f, 'before calc_pencils_*')
!!!                              call calc_pencils_hydro(f,p)
call calc_pencils_sub_circle(f,p)
!!!                              call calc_pencils_density(f,p)
call calc_pencils_eos(f,p)
!!!                              call calc_pencils_entropy(f,p)
!!!        if (llorenz_gauge)    call calc_pencils_lorenz_gauge(f,p)
!!!        if (lmagnetic)        call calc_pencils_magnetic(f,p)
!!!        if (lspecial)         call calc_pencils_special(f,p)
!
!  --------------------------------------------------------
!  NO CALLS MODIFYING PENCIL_CASE PENCILS BEYOND THIS POINT
!  --------------------------------------------------------
!
!  hydro, density, and entropy evolution
!  Note that pressure gradient is added in dss_dt to momentum,
!  even if lentropy=.false.
!
!!!        call dss_dt(f,df,p)
call calc_heatcond_tensor(f,df,p)
call calc_heatcond_glnTT(df,p)
call calc_heatcond_glnTT_iso(df,p)
!
!  End of loops over m and n.
!
        headtt=.false.
      enddo mn_loop
!
!
!  -------------------------------------------------------------
!  NO CALLS MODIFYING DF BEYOND THIS POINT (APART FROM FREEZING)
!  -------------------------------------------------------------
!
!  Reset lwrite_prof.
!
      lwrite_prof=.false.
!
    endsubroutine pde_energy_only
!***********************************************************************
    subroutine calc_pencils_sub_circle(f,p)
!
! Calculate pencil for spitzer heat conduction and
! hcond_grad and hcond_grad_iso
!
! 22 Jun F. Chen
!
    use Sub
!
    real, dimension (mx,my,mz,mfarray) :: f
    type (pencil_case) :: p
!
!    integer :: i
    real :: B2_ext
    real, dimension (nx) :: quench
    real, dimension (3) :: B_ext=(/0.0,0.0,1d-8/)
    logical :: luse_Bext_in_b2=.true.
!
    intent(inout) :: f,p
! lnrho
    p%lnrho=f(l1:l2,m,n,ilnrho)
! glnrho
    call grad(f,ilnrho,p%glnrho)
! aa
    p%aa=f(l1:l2,m,n,iax:iaz)
! aij
    call gij(f,iaa,p%aij,1)
! bb
    call curl_mn(p%aij,p%bb,p%aa)
    p%bbb=p%bb
    B2_ext=B_ext(1)**2+B_ext(2)**2+B_ext(3)**2
!
    if (B2_ext/=0.0) then
!
!  Add the external field.
!
      if (B_ext(1)/=0.0) p%bb(:,1)=p%bb(:,1)+B_ext(1)
      if (B_ext(2)/=0.0) p%bb(:,2)=p%bb(:,2)+B_ext(2)
      if (B_ext(3)/=0.0) p%bb(:,3)=p%bb(:,3)+B_ext(3)
    endif
!
! b2 (default is that B_ext is not included), but this can be changed
! by setting luse_Bext_in_b2=.true.
!
    if (luse_Bext_in_b2) then
    if (lpencil(i_b2)) call dot2_mn(p%bb,p%b2)
    else
    if (lpencil(i_b2)) call dot2_mn(p%bbb,p%b2)
    endif
! bunit
    quench = 1.0/max(tini,sqrt(p%b2))
    if (luse_Bext_in_b2) then
        p%bunit(:,1) = p%bb(:,1)*quench
        p%bunit(:,2) = p%bb(:,2)*quench
        p%bunit(:,3) = p%bb(:,3)*quench
    else
        p%bunit(:,1) = p%bbb(:,1)*quench
        p%bunit(:,2) = p%bbb(:,2)*quench
        p%bunit(:,3) = p%bbb(:,3)*quench
    endif
!
!  bij, del2a, graddiva
!  For non-cartesian coordinates jj is always required for del2a=graddiva-jj
!
    call gij_etc(f,iaa,p%aa,p%aij,p%bij)
!
    endsubroutine calc_pencils_sub_circle
!***********************************************************************
endmodule Timestep
