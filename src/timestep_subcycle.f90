! $Id$
!
! MODULE_DOC: This is a highly specified timestep module currently only working
! MODULE_DOC: together with the special module coronae.f90.
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
  include 'timestep.h'
!
  contains
!***********************************************************************
    subroutine initialize_timestep
! 
!  Coefficients for order 3.
!     
      use Messages, only: fatal_error
      use General, only: itoa
!
      if (itorder==3) then
!
!  Coefficients for order 3.  (use coefficients of Williamson (1980))
!
        alpha_ts=(/   0.0, -5/9.0 , -153/128.0 /)
        beta_ts =(/ 1/3.0, 15/16.0,    8/15.0  /)
      else
        call fatal_error('initialize_timestep','Not implemented: itorder= '// &
                         trim(itoa(itorder)))
      endif
      
    endsubroutine initialize_timestep
!***********************************************************************
    subroutine time_step(f,df,p)
!
!   2-apr-01/axel: coded
!  14-sep-01/axel: moved itorder to cdata
!
      use Boundcond
      use BorderProfiles, only: border_quenching
      use Equ, only: pde
      use Mpicomm, only: mpiallreduce_max
      use Particles_main, only: particles_timestep_first, &
          particles_timestep_second
      use Shear, only: advance_shear
      use Sub, only: shift_dt
      use Energy
      use Special
      use Boundcond
      use General, only: notanumber
!
      real, dimension (mx,my,mz,mfarray) :: f,fsub,fm1
      real, dimension (mx,my,mz) :: fm2, fj
      real, dimension (mx,my,mz,mvar) :: df,dfsub,dfm1
      type (pencil_case) :: p
      real :: ds
      real :: dt1, dt1_local, dt1_last=0.0
      integer :: ienergy, Nsub, itRKL, s
      real, dimension (nx) :: dt1_hcond_max
      real, dimension (itorder) :: coeff_fm1, coeff_fm2, coeff_fsub, coeff_dfm1, coeff_dfsub
      real :: dt1_energy_local,dt1_energy,dt_energy,dt_RKL,dt_sub_RKL
!
      ienergy = ilnTT
!
!  dt_beta_ts may be needed in other modules (like Dustdensity) for fixed dt.
!
      if (.not. ldt) dt_beta_ts=dt*beta_ts
!
!  Set up df and ds for each time sub.
!
      do itsub=1,3
        lfirst=(itsub==1)
        llast=(itsub==3)
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
          ! Get global time step limite
          call mpiallreduce_max(dt1_local,dt1,MPI_COMM_WORLD)
          dt=1.0/dt1
! in pde(f,df,p) ghost cells of f-array are set
          fsub = f
          if (loutput_varn_at_exact_tsnap) call shift_dt(dt)
          ! Timestep growth limiter
          if (ddt/=0.) dt1_last=dt1_local/ddt
        endif
!
!  Calculate dt_beta_ts.
!
        if (ldt) dt_beta_ts=dt*beta_ts
        if (ip<=6) print*, 'time_step: iproc, dt=', iproc_world, dt  !(all have same dt?)
!
!  Apply border quenching.
!
        if (lborder_profiles) call border_quenching(f,df,dt_beta_ts(itsub))
!
!  Time evolution of grid variables.
!
        f(l1:l2,m1:m2,n1:n2,1:mvar) = f(l1:l2,m1:m2,n1:n2,1:mvar) + dt_beta_ts(itsub)*df(l1:l2,m1:m2,n1:n2,1:mvar)
!
        if (lspecial) call special_after_timestep(f,df,dt_beta_ts(itsub)*ds,llast)
!
!  Increase time.
!
        t = t + dt_beta_ts(itsub)*ds
!
      enddo
!--------------------------------------------------------------------------------*
!  Now do the substeps for the energy (thermal conduction) only -----------------*
!--------------------------------------------------------------------------------*
!
      if (sub_step_hcond) then
!
!  initialized fsub
        fsub(l1:l2,m1:m2,n1:n2,ienergy)=f(l1:l2,m1:m2,n1:n2,ienergy)
!  update boundary is necessary for calculate time steps
        call boundconds_x(fsub,ilnTT,ilnTT)
        call initiate_isendrcv_bdry(fsub,ilnTT,ilnTT)
        call finalize_isendrcv_bdry(fsub,ilnTT,ilnTT)
        call boundconds_y(fsub,ilnTT,ilnTT)
        call boundconds_z(fsub,ilnTT,ilnTT)
!
!  This should be a better position to determin dt_energy, changes in lnTT
!  by other terms are taken into account.
!  Get time step for heat conduction in sub step
        call calc_hcond_timestep(fsub,p,dt1_hcond_max)
        dt1_energy_local=maxval(dt1_hcond_max(1:nx))
        call mpiallreduce_max(dt1_energy_local,dt1_energy,MPI_COMM_WORLD)
!
        dt_energy = 1d0/dt1_energy
!  Set time step to the super-timestep
        itRKL = itorder
        dt_RKL = (itRKL*itRKL+itRKL-2d0)/4d0*dt_energy
!
!  calc the number of sub iteration
        Nsub = int(dt / dt_RKL)+1
        dt_sub_RKL = dt / dble(Nsub)
!
!  calc timestep coefficient
        call RKL_coeff(itRKL, coeff_fm1, coeff_fm2, coeff_fsub, coeff_dfm1, coeff_dfsub)
!
! initalize fsub, fm1, fm2 for sub cycle with RKL steps
        fm1=fsub
!
        do j=1,Nsub
!
! one RKL step contains itRKL sub steps
!
        do s=1,itRKL
          if (s == 1) then
            dfsub=0.
            call pde_energy_only(fsub,dfsub,p,dt_sub_RKL)
            fm1(l1:l2,m1:m2,n1:n2,ienergy)=fsub(l1:l2,m1:m2,n1:n2,ienergy) &
                +coeff_dfm1(s)*dt_sub_RKL*dfsub(l1:l2,m1:m2,n1:n2,ienergy)
                fm2(l1:l2,m1:m2,n1:n2) = fsub(l1:l2,m1:m2,n1:n2,ienergy)
          else
            dfm1=0.
            call pde_energy_only(fm1,dfm1,p,dt_sub_RKL)
            fj(l1:l2,m1:m2,n1:n2)=coeff_fm1(s) * fm1(l1:l2,m1:m2,n1:n2,ienergy) &
                +coeff_fm2(s) * fm2(l1:l2,m1:m2,n1:n2)         &
                +coeff_fsub(s)*fsub(l1:l2,m1:m2,n1:n2,ienergy) &
                +coeff_dfm1(s) *dt_sub_RKL* dfm1(l1:l2,m1:m2,n1:n2,ienergy) &
                +coeff_dfsub(s)*dt_sub_RKL*dfsub(l1:l2,m1:m2,n1:n2,ienergy)
!
! set Yj-1 and Yj-2 for the next sub step
            fm2(l1:l2,m1:m2,n1:n2)=fm1(l1:l2,m1:m2,n1:n2,ienergy)
            fm1(l1:l2,m1:m2,n1:n2,ienergy)=fj(l1:l2,m1:m2,n1:n2)
          endif
        enddo ! end of itRKL sub steps
!
! set inital value for the next RKL step
        fsub(l1:l2,m1:m2,n1:n2,ienergy)=fj(l1:l2,m1:m2,n1:n2)
!
        if (notanumber(fsub(:,:,:,ienergy))) then
           print*, 'fsub contains NaN in proc',iproc_world, 'in No.',j,'subcycle'
           STOP
        endif
!
        enddo ! end of sub cycle
!
! set temperature after heat conduction back to f-array
        f(l1:l2,m1:m2,n1:n2,ienergy)=fsub(l1:l2,m1:m2,n1:n2,ienergy)
!
      ENDIF
!
    endsubroutine time_step
!***********************************************************************
    subroutine RKL_coeff(itRKL, coeff_fm1, coeff_fm2, coeff_fsub, coeff_dfm1, coeff_dfsub)
!
    real, dimension (itorder) :: a,b,coeff_fm1, coeff_fm2, coeff_fsub, coeff_dfm1, coeff_dfsub
    integer :: j,itRKL
!
    b(1) = 1d0/3d0
    b(2) = 1d0/3d0
    a(1) = 1d0-b(1)
    a(2) = 1d0-b(2)
!
    coeff_dfm1(1) = 4d0/3d0/(itRKL*itRKL+itRKL-2d0)
!
    coeff_fm1(2) = 3d0/2d0
    coeff_fm2(2) = -1d0/2d0
    coeff_fsub(2) = 1d0-coeff_fm1(2)-coeff_fm2(2)
    coeff_dfm1(2) = 6d0/(itRKL*itRKL+itRKL-2d0)
    coeff_dfsub(2) = -a(1)*coeff_dfm1(2)
!
    do j=3, itRKL
!
      b(j) = (j*j+j-2d0)/(2d0*j)/(j+1d0)
      a(j) = 1d0-b(j)
!
      coeff_fm1(j) = (2d0*j-1d0)/j*b(j)/b(j-1)
      coeff_fm2(j) = -(j-1d0)/j*b(j)/b(j-2)
      coeff_fsub(j) = 1d0-coeff_fm1(j)-coeff_fm2(j)
      coeff_dfm1(j) = 4d0*(2d0*j-1d0)/j/(itRKL*itRKL+itRKL-2d0)*b(j)/b(j-2)
      coeff_dfsub(j) = -a(j-1)*coeff_dfm1(j)
    enddo
!
    endsubroutine RKL_coeff
!***********************************************************************
    subroutine pde_energy_only(f,df,p,dt_in)
!
!  Calculate thermal conduction for sub_cycle time step.
!
!  22-jun-12/feng/bingert: coded
!
      use Boundcond
      use BorderProfiles, only: calc_pencils_borderprofiles
      use Density
      use Diagnostics
      use Energy
      use EquationOfState
      use Grid, only: calc_pencils_grid
      use Hydro
      use Lorenz_gauge
      use Magnetic
      use Special
      use Sub
!
      logical :: early_finalize
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: dt_in
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
!  Shift entire data cube by one grid point at the beginning of each
!  time-step. Useful for smearing out possible x-dependent numerical
!  diffusion, e.g. in a linear shear flow.
!
      if (lfirst .and. lshift_datacube_x) then
        call boundconds_x(f)
        f=cshift(f,1,1)
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
                     lparticles_spin.or.lsolid_cells.or. &
                     lchemistry.or.lweno_transport
!
!  Write crash snapshots to the hard disc if the time-step is very low.
!  The user must have set crash_file_dtmin_factor>0.0 in &run_pars for
!  this to be done.
!
!
!  For debugging purposes impose minimum or maximum value on certain variables.
!
! Sven.Bingert: routine not existing      
!  call impose_entropy_floor(f)
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
       if (.not. lcartesian_coords .or. .not.all(lequidist)) call get_grid_mn
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
       call calc_pencils_sub_cycle(f,p)
       call calc_pencils_eos(f,p)
!
!  --------------------------------------------------------
!  NO CALLS MODIFYING PENCIL_CASE PENCILS BEYOND THIS POINT
!  --------------------------------------------------------
!
!  hydro, density, and entropy evolution
!  Note that pressure gradient is added in denergy_dt to momentum,
!  even if lentropy=.false.
!
!!!        call denergy_dt(f,df,p)
        if (Kpara /= 0.)          call calc_heatcond_tensor(f,df,p)
        if (hcond_grad /= 0.)     call calc_heatcond_glnTT(df,p)
        if (hcond_grad_iso /= 0.) call calc_heatcond_glnTT_iso(df,p)
        if (lfilter_farray)       call filter_lnTT(f,df,dt_in)
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
    subroutine calc_pencils_sub_cycle(f,p)
!
! Calculate pencil for spitzer heat conduction and
! hcond_grad and hcond_grad_iso
!
! 22 Jun F. Chen
!
    use Magnetic, only: get_bext
    use Messages, only: fatal_error
    use Sub, only: grad, gij, curl_mn, dot2_mn, gij_etc
!
    real, dimension (mx,my,mz,mfarray) :: f
    type (pencil_case) :: p
!
    real :: B2_ext
    real, dimension (nx) :: quench
    logical :: luse_Bext_in_b2=.true.
    real, dimension(3) :: B_ext
!
    intent(inout) :: f,p
!
!  Get the external magnetic field if exists.
!
    if (lmagnetic) call get_bext(B_ext)
!
    p%lnrho=f(l1:l2,m,n,ilnrho)
    call grad(f,ilnrho,p%glnrho)
    p%aa=f(l1:l2,m,n,iax:iaz)
    call gij(f,iaa,p%aij,1)
    call curl_mn(p%aij,p%bb,p%aa)
    p%bbb=p%bb
    B2_ext=B_ext(1)**2+B_ext(2)**2+B_ext(3)**2
!
    if (B2_ext/=0.0) then
      if (B_ext(1)/=0.0) p%bb(:,1)=p%bb(:,1)+B_ext(1)
      if (B_ext(2)/=0.0) p%bb(:,2)=p%bb(:,2)+B_ext(2)
      if (B_ext(3)/=0.0) p%bb(:,3)=p%bb(:,3)+B_ext(3)
    endif
!
!  b2 now (since 18 June 2013) includes B_ext by default.
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
    endsubroutine calc_pencils_sub_cycle
!***********************************************************************
    subroutine pushpars2c(p_par)

    integer, parameter :: n_pars=2
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr_c(alpha_ts,p_par(1))  ! (3)
    call copy_addr_c(beta_ts ,p_par(2))  ! (3)

    endsubroutine pushpars2c
!***********************************************************************
endmodule Timestep
