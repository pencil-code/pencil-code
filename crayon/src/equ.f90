! $Id: equ.f90 13579 2010-03-31 11:58:14Z AxelBrandenburg $
!
!  This module adds evolution terms to the dynamical equations, calling
!  subroutines in the chosen set of physics modules.  
!
module Equ
!
  use Cdata
  use Messages
!
  implicit none
!
  private
!
  public :: pde, debug_imn_arrays, initialize_pencils
!
  contains
!***********************************************************************
    include 'pencil_init.inc' ! defines subroutine initialize_pencils()
!***********************************************************************
    subroutine pde(f,df,p)
!
!  Call the different evolution equations.
!
!  10-sep-01/axel: coded
!
      use Boundcond
      use Density
      use Diagnostics
      use Entropy
      use EquationOfState
      use Gravity
      use Grid, only: calc_pencils_grid
      use Hydro
      use Magnetic
      use Mpicomm
      use Shear
      use Shock, only: calc_pencils_shock, calc_shock_profile, &
                       calc_shock_profile_simple
      use Sub
      use Viscosity, only: calc_viscosity, calc_pencils_viscosity
!
      logical :: early_finalize
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: maxadvec,advec2,maxdiffus,maxdiffus3
      intent(inout)  :: f       ! inout due to  lshift_datacube_x,
                                ! density floor, or velocity ceiling
      intent(out)    :: df, p
!
!  Print statements when they are first executed.
!
      headtt = headt .and. lfirst .and. lroot
!
      if (headtt.or.ldebug) print*,'pde: ENTER'
      if (headtt) call svn_id( &
           "$Id: equ.f90 13579 2010-03-31 11:58:14Z AxelBrandenburg $")
!
!  Initialize counter for calculating and communicating print results.
!  Do diagnostics only in the first of the 3 (=itorder) substeps.
!
      ldiagnos   =lfirst.and.lout
      l1davgfirst=lfirst.and.l1davg
      l2davgfirst=lfirst.and.l2davg
!
!  Record times for diagnostic and 2d average output.
!
      if (ldiagnos)    tdiagnos=t    ! (diagnostics are for THIS time)
      if (l1davgfirst) t1ddiagnos=t  ! (1-D averages are for THIS time)
      if (l2davgfirst) t2davgfirst=t ! (2-D averages are for THIS time)
!
!  Grid spacing. 
!
      dxyz_2 = dx_1(l1:l2)**2+dy_1(m1)**2+dz_1(n1)**2
      dxyz_4 = dx_1(l1:l2)**4+dy_1(m1)**4+dz_1(n1)**4
      dxyz_6 = dx_1(l1:l2)**6+dy_1(m1)**6+dz_1(n1)**6
!
!  Write crash snapshots to the hard disc if the time-step is very low.
!  The user must have set crash_file_dtmin_factor>0.0 in &run_pars for
!  this to be done.
!
      if (crash_file_dtmin_factor > 0.0) call output_crash_files(f)      
!
!  Call "before_boundary" hooks (for f array precalculation)
!
      if (ldensity)      call density_before_boundary(f)
      if (lshear)        call shear_before_boundary(f)
!
!  Initiate shock profile calculation and use asynchronous to handle
!  communication along processor/periodic boundaries.
!
      if (lshock) call calc_shock_profile(f)
!
!  Prepare x-ghost zones; required before f-array communication
!  AND shock calculation
!
      call boundconds_x(f)
!
!  Initiate (non-blocking) communication and do boundary conditions.
!  Required order:
!  1. x-boundaries (x-ghost zones will be communicated) - done above
!  2. communication
!  3. y- and z-boundaries
!
      if (ldebug) print*,'pde: bef. initiate_isendrcv_bdry'
      call initiate_isendrcv_bdry(f)
      if (early_finalize) then
        call finalize_isendrcv_bdry(f)
        call boundconds_y(f)
        call boundconds_z(f)
      endif
!
!  Set inverse timestep to zero before entering loop over m and n.
!
      if (lfirst.and.ldt) then
        if (dtmax/=0.0) then
          dt1_max=1/dtmax
        else
          dt1_max=0.0
        endif
      endif
!
!  Calculate shock profile (simple).
!
      if (lshock) call calc_shock_profile_simple(f)
!
!  Calculate averages, currently only required for certain settings
!  in hydro of the testfield procedure (only when lsoca=.false.)
!
      call timing('pde','before mn-loop')
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
!  Store the velocity part of df array in a temporary array 
!  while solving the anelastic case.
!
        call timing('pde','before ldensity_anelastic',mnloop=.true.)
!
!  Make sure all ghost points are set.
!
        if (.not.early_finalize.and.necessary(imn)) then
          call finalize_isendrcv_bdry(f)
          call boundconds_y(f)
          call boundconds_z(f)
        endif
        call timing('pde','finished boundconds_z',mnloop=.true.)
!
!  For each pencil, accumulate through the different modules
!  advec_XX and diffus_XX, which are essentially the inverse
!  advective and diffusive timestep for that module.
!  (note: advec_cs2 and advec_va2 are inverse _squared_ timesteps)
!
        if (lfirst.and.ldt.and.(.not.ldt_paronly)) then
          if (lhydro) then
            advec_uu=0.0
          endif
          if (ldensity) then
            diffus_diffrho=0.0; advec_lnrho=0.0
            if (leos) advec_cs2=0.0
          endif
          if (lentropy) then
            diffus_chi=0.0
          endif
          if (lmagnetic) then
            advec_va2=0.0
            diffus_eta=0.0
          endif
          if (lviscosity) then
            diffus_nu=0.0
          endif
          if (lshear) then
            advec_shear=0.0
          endif
        endif
!
!  Grid spacing.
!
        dxyz_2 = dx_1(l1:l2)**2+dy_1(m1)**2+dz_1(n1)**2
        dxyz_4 = dx_1(l1:l2)**4+dy_1(m1)**4+dz_1(n1)**4
        dxyz_6 = dx_1(l1:l2)**6+dy_1(m1)**6+dz_1(n1)**6
!
!  Calculate grid/geometry related pencils.
!
        call calc_pencils_grid(f,p)
!
!  Calculate pencils for the pencil_case.
!  Note: some no-modules (e.g. nohydro) also calculate some pencils,
!  so it would be wrong to check for lhydro etc in such cases.
!
                              call calc_pencils_hydro(f,p)
                              call calc_pencils_density(f,p)
                              call calc_pencils_eos(f,p)
        if (lshock)           call calc_pencils_shock(f,p)
        if (lviscosity)       call calc_pencils_viscosity(f,p)
        if (lmagnetic)        call calc_pencils_magnetic(f,p)
        if (lgrav)            call calc_pencils_gravity(f,p)
        if (lshear)           call calc_pencils_shear(f,p)
!
!  --------------------------------------------------------
!  NO CALLS MODIFYING PENCIL_CASE PENCILS BEYOND THIS POINT
!  --------------------------------------------------------
!
!  hydro, density, entropy, and magnetic evolution
!  Note that pressure gradient is added in dss_dt to momentum,
!  even if lentropy=.false.
!
        call duu_dt(df,p)
        call dlnrho_dt(df,p)
        call dss_dt(df,p)
        call daa_dt(df,p)
!
!  Add gravity, if present
!
        if (lgrav) then
          if (lhydro) call duu_dt_grav(f,df,p)
        endif
!
!  Add shear if present
!
        if (lshear) call shearing(f,df,p)
!
!  In max_mn maximum values of u^2 (etc) are determined sucessively
!  va2 is set in magnetic (or nomagnetic)
!  In rms_mn sum of all u^2 (etc) is accumulated
!  Calculate maximum advection speed for timestep; needs to be done at
!  the first substep of each time step
!  Note that we are (currently) accumulating the maximum value,
!  not the maximum squared!
!
!  The dimension of the run ndim (=0, 1, 2, or 3) enters the viscous time step.
!  This has to do with the term on the diagonal, cdtv depends on order of scheme
!
        if (lfirst.and.ldt.and.(.not.ldt_paronly)) then
!
!  sum or maximum of the advection terms?
!  (lmaxadvec_sum=.false. by default)
!
          maxadvec=0.0
          maxdiffus=0.0
          if (lhydro) maxadvec=maxadvec+advec_uu
          if (lshear) maxadvec=maxadvec+advec_shear
          if (ldensity.or.lmagnetic) then
            advec2=0.0
            if (ldensity) advec2=advec2+advec_cs2+advec_lnrho**2
            if (lmagnetic) advec2=advec2+advec_va2
            maxadvec=maxadvec+sqrt(advec2)
          endif
!
!  Time step constraints from each module.
!  (At the moment, magnetic and testfield use the same variable.)
!
          if (lviscosity) maxdiffus=max(diffus_nu,maxdiffus)
          if (ldensity)   maxdiffus=max(diffus_diffrho,maxdiffus)
          if (lentropy)   maxdiffus=max(diffus_chi,maxdiffus)
          if (lmagnetic)  maxdiffus=max(diffus_eta,maxdiffus)
!
!  cdt, cdtv, and cdtc are empirical non-dimensional coefficients
!
          dt1_advec  = maxadvec/cdt
          dt1_diffus = maxdiffus/cdtv
          dt1_max    = max(dt1_max,sqrt(dt1_advec**2+dt1_diffus**2))
!
!  Diagnostics showing how close to advective and diffusive time steps we are
!
          if (ldiagnos.and.idiag_dtv/=0) then
            call max_mn_name(maxadvec/cdt,idiag_dtv,l_dt=.true.)
          endif
          if (ldiagnos.and.idiag_dtdiffus/=0) then
            call max_mn_name(maxdiffus/cdtv,idiag_dtdiffus,l_dt=.true.)
          endif
!
!  Regular and hyperdiffusive mesh Reynolds numbers
!
          if (ldiagnos) then
            if (idiag_Rmesh/=0) &
                call max_mn_name(pi_1*maxadvec/(maxdiffus+tini),idiag_Rmesh)
            if (idiag_Rmesh3/=0) &
                call max_mn_name(pi5_1*maxadvec/(maxdiffus3+tini),idiag_Rmesh3)
            if (idiag_maxadvec/=0) &
                call max_mn_name(maxadvec,idiag_maxadvec)
          endif
        endif
!
        call timing('pde','end of mn loop',mnloop=.true.)
!
!  End of loops over m and n.
!
        headtt=.false.
      enddo mn_loop
      call timing('pde','at the end of the mn_loop')
!
!  ----------------------------------------
!  NO CALLS MODIFYING DF BEYOND THIS POINT 
!  ----------------------------------------
!
!  Check for NaNs in the advection time-step.
!
      if (notanumber(dt1_advec)) then
        print*, 'pde: dt1_advec contains a NaN at iproc=', iproc
        if (lhydro)           print*, 'advec_uu   =',advec_uu
        if (lshear)           print*, 'advec_shear=',advec_shear
        if (lentropy)         print*, 'advec_cs2  =',advec_cs2
        if (lmagnetic)        print*, 'advec_va2  =',advec_va2
        call fatal_error_local('pde','')
      endif
!
!  Diagnostics.
!
      if (ldiagnos) call diagnostic
!
!  1-D diagnostics.
!
      if (l1davgfirst) then
        call xyaverages_z
        call xzaverages_y
        call yzaverages_x
      endif
!
!  2-D averages.
!
      if (l2davgfirst) then
        if (lwrite_yaverages)   call yaverages_xz
        if (lwrite_zaverages)   call zaverages_xy
      endif
!
!  Note: zaverages_xy are also needed if bmx and bmy are to be calculated
!  (of course, yaverages_xz does not need to be calculated for that).
!
      if (.not.l2davgfirst.and.ldiagnos.and.ldiagnos_need_zaverages) then
        if (lwrite_zaverages) call zaverages_xy
      endif
!
!  Reset lwrite_prof.
!
      lwrite_prof=.false.
!
    endsubroutine pde
!***********************************************************************
    subroutine debug_imn_arrays
!
!  For debug purposes: writes out the mm, nn, and necessary arrays.
!
!  23-nov-02/axel: coded
!
      open(1,file=trim(directory)//'/imn_arrays.dat')
      do imn=1,ny*nz
        if (necessary(imn)) write(1,'(a)') '----necessary=.true.----'
        write(1,'(4i6)') imn,mm(imn),nn(imn)
      enddo
      close(1)
!
    endsubroutine debug_imn_arrays
!***********************************************************************
    subroutine output_crash_files(f)
!
!  Write crash snapshots when time-step is low.
!
!  15-aug-2007/anders: coded
!
      use Snapshot
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      integer, save :: icrash=0
      character (len=10) :: filename
      character (len=1) :: icrash_string
!
      if ( (it>1) .and. (itsub==1) .and. (dt<=crash_file_dtmin_factor*dtmin) ) then
        write(icrash_string, fmt='(i1)') icrash
        filename='crash'//icrash_string//'.dat'
        call wsnap(trim(directory_snap)//'/'//filename,f,mvar_io,.false.)
        if (lroot) then
          print*, 'Time-step is very low - writing '//trim(filename)
          print*, '(it, itsub=', it, itsub, ')'
          print*, '(t, dt=', t, dt, ')'
        endif
!
!  Next crash index, cycling from 0-9 to avoid excessive writing of
!  snapshots to the hard disc.
!        
        icrash=icrash+1
        icrash=mod(icrash,10)
      endif
!
    endsubroutine output_crash_files
!***********************************************************************
endmodule Equ
