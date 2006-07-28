
! $Id: equ.f90,v 1.313 2006-07-28 13:23:30 wlyra Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!
!***************************************************************

module Equ
!
  use Cdata
  use Messages
!
  implicit none
!
  private
!
  public :: pde, debug_imn_arrays
  public :: initialize_pencils,pencil_consistency_check
!
  contains

!***********************************************************************

      include 'pencil_init.inc' ! defines subroutine initialize_pencils()

!***********************************************************************
      subroutine collect_UUmax
!
!  Calculate the maximum effective advection velocity in the domain;
!  needed for determining dt at each timestep
!
!   2-sep-01/axel: coded
!
      use Mpicomm
      use Cdata
      use Sub
!
      real, dimension(1) :: fmax_tmp,fmax
!
!  communicate over all processors
!  the result is then present only on the root processor
!  reassemble using old names
!
      fmax_tmp(1)=UUmax
      call mpireduce_max(fmax_tmp,fmax,1)
      if (lroot) UUmax=fmax(1)
!
      endsubroutine collect_UUmax
!***********************************************************************
    subroutine diagnostic
!
!  calculate diagnostic quantities
!   2-sep-01/axel: coded
!  14-aug-03/axel: began adding surface integrals
!
      use Mpicomm
      use Cdata
      use Sub
!
      real, dimension (mname) :: fmax_tmp, fsum_tmp, fmax, fsum, fweight_tmp
      real :: dv
      integer :: iname,imax_count,isum_count,nmax_count,nsum_count
      logical :: lweight_comm
!
!  go through all print names, and sort into communicators
!  corresponding to their type
!
      imax_count=0
      isum_count=0
      lweight_comm=.false.
      do iname=1,nname
        if (itype_name(iname)<0) then
          imax_count=imax_count+1
          fmax_tmp(imax_count)=fname(iname)
        elseif (itype_name(iname)>0) then
          isum_count=isum_count+1
          fsum_tmp(isum_count)=fname(iname)
          if (itype_name(iname)==ilabel_sum_weighted .or. &
              itype_name(iname)==ilabel_sum_weighted_sqrt) then
            fweight_tmp(isum_count)=fweight(iname)
            lweight_comm=.true.
          endif
        endif
      enddo
      nmax_count=imax_count
      nsum_count=isum_count
!
!  communicate over all processors
!
      call mpireduce_max(fmax_tmp,fmax,nmax_count)
      call mpireduce_sum(fsum_tmp,fsum,nsum_count)
      if (lweight_comm) call mpireduce_sum(fweight_tmp,fweight,nsum_count)
!
!  the result is present only on the root processor
!
      if (lroot) then
!
!  sort back into original array
!
         imax_count=0
         isum_count=0
         do iname=1,nname
           if (itype_name(iname)<0) then ! max
             imax_count=imax_count+1

             if (itype_name(iname)==ilabel_max)            &
                 fname(iname)=fmax(imax_count)

             if (itype_name(iname)==ilabel_max_sqrt)       &
                 fname(iname)=sqrt(fmax(imax_count))

             if (itype_name(iname)==ilabel_max_dt)         &
                 fname(iname)=fmax(imax_count)

             if (itype_name(iname)==ilabel_max_neg)        &
                 fname(iname)=-fmax(imax_count)

             if (itype_name(iname)==ilabel_max_reciprocal) &
                 fname(iname)=1./fmax(imax_count)

           elseif (itype_name(iname)>0) then ! sum
             isum_count=isum_count+1

             if (itype_name(iname)==ilabel_sum)            &
                 fname(iname)=fsum(isum_count)/(nw*ncpus)

             if (itype_name(iname)==ilabel_sum_sqrt)       &
                 fname(iname)=sqrt(fsum(isum_count)/(nw*ncpus))

             if (itype_name(iname)==ilabel_sum_par)        &
                 fname(iname)=fsum(isum_count)/npar

             if (itype_name(iname)==ilabel_integrate) then
               dv=1.
               if (nxgrid/=1) dv=dv*dx
               if (nygrid/=1) dv=dv*dy
               if (nzgrid/=1) dv=dv*dz
               fname(iname)=fsum(isum_count)*dv
              endif

              if (itype_name(iname)==ilabel_surf)          & 
                  fname(iname)=fsum(isum_count)

             if (itype_name(iname)==ilabel_sum_weighted) then
               if (fweight(isum_count)/=0.0) then
                 fname(iname)=fsum(isum_count)/fweight(isum_count)
               else
                 fname(iname)=0.0
               endif
             endif

             if (itype_name(iname)==ilabel_sum_weighted_sqrt) then
               if (fweight(isum_count)/=0.0) then
                 fname(iname)=sqrt(fsum(isum_count)/fweight(isum_count))
               else
                 fname(iname)=0.0
               endif
             endif

           endif

         enddo
!
      endif
!
    endsubroutine diagnostic
!***********************************************************************
    subroutine xyaverages_z()
!
!  Calculate xy-averages (still depending on z)
!  NOTE: these averages depend on z, so after summation in x and y they
!  are still distributed over nprocz CPUs; hence the dimensions of fsumz
!  (and fnamez).
!  In other words: the whole xy-average is present in one and the same fsumz,
!  but the result is not complete on any of the processors before
!  mpireduce_sum has been called. This is simpler than collecting results
!  first in processors with the same ipz and different ipy, and then
!  assemble result from the subset of ipz processors which have ipy=0
!  back on the root processor.
!
!   6-jun-02/axel: coded
!
      use Mpicomm
      use Cdata
!
      real, dimension (nz,nprocz,mnamez) :: fsumz
!
!  communicate over all processors
!  the result is only present on the root processor
!
      if (nnamez>0) then
        call mpireduce_sum(fnamez,fsumz,nz*nprocz*nnamez)
        if (lroot) &
            fnamez(:,:,1:nnamez)=fsumz(:,:,1:nnamez)/(nx*ny*nprocy)
      endif
!
    endsubroutine xyaverages_z
!***********************************************************************
    subroutine xzaverages_y()
!
!  Calculate xz-averages (still depending on y)
!
!  12-oct-05/anders: adapted from xyaverages_z
!
      use Mpicomm
      use Cdata
!
      real, dimension (ny,nprocy,mnamey) :: fsumy
!
!  communicate over all processors
!  the result is only present on the root processor
!
      if (nnamey>0) then
        call mpireduce_sum(fnamey,fsumy,ny*nprocy*nnamey)
        if (lroot) &
            fnamey(:,:,1:nnamey)=fsumy(:,:,1:nnamey)/(nx*nz*nprocz)
      endif
!
    endsubroutine xzaverages_y
!***********************************************************************
    subroutine yzaverages_x()
!
!  Calculate yz-averages (still depending on x)
!
!   2-oct-05/anders: adapted from xyaverages_z
!
      use Mpicomm
      use Cdata
!
      real, dimension (nx,mnamex) :: fsumx
!
!  communicate over all processors
!  the result is only present on the root processor
!
      if (nnamex>0) then
        call mpireduce_sum(fnamex,fsumx,nx*nnamex)
        if (lroot) &
            fnamex(:,1:nnamex)=fsumx(:,1:nnamex)/(ny*nprocy*nz*nprocz)
      endif
!
    endsubroutine yzaverages_x
!***********************************************************************
    subroutine yaverages_xz()
!
!  Calculate y-averages (still depending on x and z)
!  NOTE: these averages depend on x and z, so after summation in y they
!  are still distributed over nprocy CPUs; hence the dimensions of fsumxz
!  (and fnamexz).
!
!   7-jun-05/axel: adapted from zaverages_xy
!
      use Mpicomm
      use Cdata
      use Sub
!
      real, dimension (nx,nz,nprocz,mnamexz) :: fsumxz
!
!  communicate over all processors
!  the result is only present on the root processor
!
      if (nnamexz>0) then
        call mpireduce_sum(fnamexz,fsumxz,nnamexz*nx*nz*nprocz)
        if (lroot) &
            fnamexz(:,:,:,1:nnamez)=fsumxz(:,:,:,1:nnamez)/(ny*nprocy)
      endif
!
    endsubroutine yaverages_xz
!***********************************************************************
    subroutine zaverages_xy()
!
!  Calculate z-averages (still depending on x and y)
!  NOTE: these averages depend on x and y, so after summation in z they
!  are still distributed over nprocy CPUs; hence the dimensions of fsumxy
!  (and fnamexy).
!
!  19-jun-02/axel: coded
!
      use Mpicomm
      use Cdata
      use Sub
!
      real, dimension (nx,ny,nprocy,mnamexy) :: fsumxy
!
!  communicate over all processors
!  the result is only present on the root processor
!
      if (nnamexy>0) then
        call mpireduce_sum(fnamexy,fsumxy,nnamexy*nx*ny*nprocy)
        if (lroot) &
            fnamexy(:,:,:,1:nnamexy)=fsumxy(:,:,:,1:nnamexy)/(nz*nprocz)
      endif
!
    endsubroutine zaverages_xy
!***********************************************************************
    subroutine phiaverages_rz()
!
!  calculate azimuthal averages (as functions of r_cyl,z)
!  NOTE: these averages depend on (r and) z, so after summation they
!  are still distributed over nprocz CPUs; hence the dimensions of fsumrz
!  (and fnamerz).
!
!  9-dec-02/wolf: coded
!
      use Mpicomm
      use Cdata
      use Sub
!
      integer :: i
      real, dimension (nrcyl,0:nz,nprocz,mnamerz) :: fsumrz
!
!  communicate over all processors
!  the result is only present on the root processor
!  normalize by sum of unity which is accumulated in fnamerz(:,0,:,1)
!
      if (nnamerz>0) then
        call mpireduce_sum(fnamerz,fsumrz,mnamerz*nrcyl*(nz+1)*nprocz)
        if (lroot) then
          do i=1,nnamerz
            fnamerz(:,1:nz,:,i)=fsumrz(:,1:nz,:,i)/spread(fsumrz(:,0,:,1),2,nz)
          enddo
        endif
      endif
!
    endsubroutine phiaverages_rz
!***********************************************************************
    subroutine pde(f,df,p)
!
!  call the different evolution equations (now all in their own modules)
!
!  10-sep-01/axel: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Global
      use Hydro
      use Gravity
      use Selfgravity
      use Entropy
      use Magnetic
      use Testfield
      use Radiation
      use EquationOfState
      use Pscalar
      use Chiral
      use Dustvelocity
      use Dustdensity
      use Cosmicray
      use CosmicrayFlux
      use Special
      use Boundcond
      use Shear
      use Density
      use Shock, only: calc_pencils_shock, calc_shock_profile, calc_shock_profile_simple
      use Viscosity, only: calc_viscosity, calc_pencils_viscosity, &
                           lvisc_first, idiag_epsK
      use Particles_main
      use Planet
      use BorderProfiles
!
      logical :: early_finalize
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: maxadvec,maxdiffus,pfreeze,pfreeze_int,pfreeze_ext
      integer :: iv
!
!  print statements when they are first executed
!
      headtt = headt .and. lfirst .and. lroot
!
      if (headtt.or.ldebug) print*,'pde: ENTER'
      if (headtt) call cvs_id( &
           "$Id: equ.f90,v 1.313 2006-07-28 13:23:30 wlyra Exp $")
!
!  initialize counter for calculating and communicating print results
!
      ldiagnos=lfirst.and.lout
      l2davgfirst=lfirst.and.l2davg
!
!  record times for diagnostic and 2d average output
!
      if (ldiagnos) tdiagnos=t !(diagnostics are for THIS time)
      if (l2davgfirst) t2davgfirst=t !(2-D averages are for THIS time)
!
!  need to finalize communication early either for test purposes, or
!  when radiation transfer of global ionization is calculatearsd.
!  This could in principle be avoided (but it not worth it now)
!
      early_finalize=test_nonblocking.or.leos_ionization.or.lradiation_ray
!
!  Calculate the potential of the self gravity. Must be done before
!  communication in order to be able to take the gradient of the potential
!  later.
!
      call calc_selfpotential(f)
!
!  Check for dust grain mass interval overflows
!  (should consider having possibility for all modules to fiddle with the
!   f array before boundary conditions are sent)
!
      if (ldustdensity .and. ldustnulling) call null_dust_vars(f)
      if (ldustdensity .and. lmdvar .and. itsub==1) call redist_mdbins(f)
!
! Prepare x-ghost zones required before f-array communication AND shock calculation
!
      call boundconds_x(f)
!
!  Initiate shock profile calculation and use asynchronous to handle
!  communication along processor/periodic boundaries.
!
      if (lshock) call calc_shock_profile(f)
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
        call boundconds_z(f,df=df)
      endif
!
!  Apply global boundary conditions to particle positions and communiate
!  migrating particles between the processors.
!
!
      if (lparticles) call particles_boundconds(f)
!
!  set inverse timestep to zero before entering loop over m and n
!
      if (lfirst.and.ldt) then
        if (dtmax/=0) then 
          dt1_max=1./dtmax
        else
          dt1_max=0.
        endif
      endif
!
!  Calculate ionization degree (needed for thermodynamics)
!  Radiation transport along rays
!
      if (leos_ionization.or.leos_temperature_ionization) call ioncalc(f)
      if (lradiation_ray)  call radtransfer(f)
      if (lshock) call calc_shock_profile_simple(f)
      if (lvisc_hyper.or.lvisc_smagorinsky) then
        if ((lvisc_first.and.lfirst).or..not.lvisc_first) call calc_viscosity(f)
      endif
!
!  Calculate averages for testfield procedure (only when lsoca=.false.)
!
      if (ltestfield) call calc_ltestfield_pars(f)
!
!  do loop over y and z
!  set indices and check whether communication must now be completed
!  if test_nonblocking=.true., we communicate immediately as a test.
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        lfirstpoint=(imn==1)      ! true for very first m-n loop
        llastpoint=(imn==(ny*nz)) ! true for very last m-n loop

!        if (loptimise_ders) der_call_count=0 !DERCOUNT
!
! make sure all ghost points are set
!
        if (.not.early_finalize.and.necessary(imn)) then
          call finalize_isendrcv_bdry(f)
          call boundconds_y(f)
          call boundconds_z(f,df=df)
        endif
!
!  coordinates are needed frequently
!  --- but not for isotropic turbulence; and there are many other
!  circumstances where this is not needed.
!  Note: cylindrical radius currently only needed for phi-averages.
!
        call calc_unitvects_sphere()
!
!  calculate profile for phi-averages if needed
!  Note that rcyl_mn is also needed for Couette flow experiments,
!  so let's hope that everybody remembers to do averages as well...
!
        if (l2davgfirst.and.lwrite_phiaverages) then
          call calc_phiavg_general()
          call calc_phiavg_profile()
          call calc_phiavg_unitvects()
        elseif (lcylindrical) then
          call calc_phiavg_general()
        endif
!
!  general phiaverage quantities -- useful for debugging
!
        if (l2davgfirst) then
          call phisum_mn_name_rz(rcyl_mn,idiag_rcylmphi)
          call phisum_mn_name_rz(phi_mn,idiag_phimphi)
          call phisum_mn_name_rz(z_mn,idiag_zmphi)
          call phisum_mn_name_rz(r_mn,idiag_rmphi)
        endif
!
!  For each pencil, accumulate through the different modules
!  advec_XX and diffus_XX, which are essentially the inverse
!  advective and diffusive timestep for that module.
!  (note: advec_cs2 and advec_va2 are inverse _squared_ timesteps)
!  
        if (lfirst) then
          advec_uu=0.; advec_shear=0.; advec_hall=0.
          advec_cs2=0.; advec_va2=0.; advec_crad2=0.; advec_uud=0;
          diffus_pscalar=0.
          diffus_chiral=0.; diffus_diffrho=0.; diffus_cr=0.
          diffus_eta=0.; diffus_nu=0.; diffus_chi=0.
          diffus_nud=0.; diffus_diffnd=0.
        endif
!
!  The following is only kept for backwards compatibility.
!  Will be deleted in the future.
!
        if (old_cdtv) then
          dxyz_2 = max(dx_1(l1:l2)**2,dy_1(m)**2,dz_1(n)**2)
        else
          dxyz_2 = dx_1(l1:l2)**2+dy_1(m)**2+dz_1(n)**2
          dxyz_4 = dx_1(l1:l2)**4+dy_1(m)**4+dz_1(n)**4
          dxyz_6 = dx_1(l1:l2)**6+dy_1(m)**6+dz_1(n)**6
        endif
!
!  Calculate pencils for the pencil_case
!
                            call calc_pencils_hydro(f,p)
                            call calc_pencils_density(f,p)
                            call calc_pencils_eos(f,p)
        if (lshock)         call calc_pencils_shock(f,p)
        if (lviscosity)     call calc_pencils_viscosity(f,p)
                            call calc_pencils_entropy(f,p)
                            call calc_pencils_magnetic(f,p)
        if (lgrav)          call calc_pencils_gravity(f,p)
        if (lselfgravity)   call calc_pencils_selfgravity(f,p)
        if (lpscalar)       call calc_pencils_pscalar(f,p)
        if (ldustvelocity)  call calc_pencils_dustvelocity(f,p)
        if (ldustdensity)   call calc_pencils_dustdensity(f,p)
        if (lcosmicray)     call calc_pencils_cosmicray(f,p)
        if (lcosmicrayflux) call calc_pencils_cosmicrayflux(f,p)
        if (lchiral)        call calc_pencils_chiral(f,p)
        if (lradiation)     call calc_pencils_radiation(f,p)
        if (lspecial)       call calc_pencils_special(f,p)
        if (lplanet)        call calc_pencils_planet(f,p)
        if (lparticles)     call particles_calc_pencils(f,p)
!
!  --------------------------------------------------------
!  NO CALLS MODIFYING PENCIL_CASE PENCILS BEYOND THIS POINT
!  --------------------------------------------------------
!
!  Average of phi-variables to remove mean flow
!  and calculate stresses
!
        if ((lplanet).and.(lfirst)) &
             call planet_phiavg(p)
!
!  hydro, density, and entropy evolution
!
        call duu_dt(f,df,p)
        call dlnrho_dt(f,df,p)
        call dss_dt(f,df,p)
!
!  Magnetic field evolution
!
        call daa_dt(f,df,p)
!
!  Testfield evolution
!
        if (ltestfield) call daatest_dt(f,df,p)
!
!  Passive scalar evolution
!
        call dlncc_dt(f,df,p)
!
!  Dust evolution
!
        call duud_dt(f,df,p)
        call dndmd_dt(f,df,p)
!
!  Add gravity, if present
!
        if (lgrav) then
          if (lhydro.or.ldustvelocity) then 
             call duu_dt_grav(f,df,p)             
          endif
        endif
!
!  Self-gravity
!
        if (lselfgravity) call duu_dt_selfgrav(f,df,p)
!
!  cosmic ray energy density
!
        if (lcosmicray) call decr_dt(f,df,p)
!
!  cosmic ray flux
!
        if (lcosmicrayflux) call dfcr_dt(f,df,p)
!
!  chirality of left and right handed aminoacids
!
        if (lchiral) call dXY_chiral_dt(f,df,p)
!
!  Evolution of radiative energy
!
        if (lradiation_fld) call de_dt(f,df,p,gamma)
!
!  Add and extra 'special' physics 
!
        if (lspecial)                    call dspecial_dt(f,df,p)
!
!  Add radiative cooling and radiative pressure (for ray method)
!
        if (lradiation_ray.and.(lentropy.or.ltemperature)) then
          call radiative_cooling(f,df,p)
          call radiative_pressure(f,df,p)
        endif
!
!  Add shear if present
!
        if (lshear)                      call shearing(f,df)
!
!  Add planet gravity
!
        if (lplanet) then
           call auxcall_gravcomp(g0,r0_pot,n_pot,p)
        endif
!
        if (lparticles) call particles_pde_pencil(f,df,p)
!
! Apply border profile before freezing
!
        if (lborder_profiles) then
           do iv=1,nvar
              call border_driving(f,df,iv)
          enddo
        endif
!
!  -------------------------------------------------------------
!  NO CALLS MODIFYING DF BEYOND THIS POINT (APART FROM FREEZING)
!  -------------------------------------------------------------
!
!  Similar to velocity damping, but more radical: set some df=0 for r_mn<r_int:

        if (any(lfreeze_varint)) then
          if (headtt) &
              print*, 'pde: freezing variables for r < ', rfreeze_int, &
                  ' : ', lfreeze_varint
          if (lcylindrical) then
            pfreeze_int = &
                quintic_step(rcyl_mn,rfreeze_int,wfreeze_int,SHIFT=fshift_int)
          else
            pfreeze_int = &
                quintic_step(r_mn   ,rfreeze_int,wfreeze_int,SHIFT=fshift_int)
          endif
!          
          do iv=1,nvar
            if (lfreeze_varint(iv)) &
                df(l1:l2,m,n,iv) = pfreeze_int*df(l1:l2,m,n,iv)
          enddo
        endif

        if (any(lfreeze_varext)) then
          if (headtt) &
              print*, 'pde: freezing variables for r > ', rfreeze_ext, &
                  ' : ', lfreeze_varext
          if (lcylindrical) then
            pfreeze_ext = &
                1-quintic_step(rcyl_mn,rfreeze_ext,wfreeze_ext,SHIFT=fshift_ext)
          else
            pfreeze_ext = &
                1-quintic_step(r_mn   ,rfreeze_ext,wfreeze_ext,SHIFT=fshift_ext)
          endif
!
          do iv=1,nvar
            if (lfreeze_varext(iv)) &
                df(l1:l2,m,n,iv) = pfreeze_ext*df(l1:l2,m,n,iv)
          enddo
        endif

        if (any(lfreeze_varsquare)) then
          if (headtt) &
              print*, 'pde: freezing variables inside square : ', &
                  lfreeze_varsquare
          pfreeze = 1. - quintic_step(x_mn,xfreeze_square, wfreeze,SHIFT=-1.) &
                       * quintic_step(y_mn,yfreeze_square,-wfreeze,SHIFT=-1.)
!
          do iv=1,nvar          
            if (lfreeze_varsquare(iv)) &
                df(l1:l2,m,n,iv) = pfreeze*df(l1:l2,m,n,iv)
          enddo
        endif

!
!  Freeze components of variables in boundary slice if specified by boundary
!  condition 'f'
!
!  freezing boundary conditions in x
!
        if (lfrozen_bcs_x) then ! are there any frozen vars at all?
!
! Only need to do this for nonperiodic x direction, on left/right-most
! processor and in left/right--most pencils
!
          if (.not. lperi(1)) then
            if (ipx == 0) then
              do iv=1,nvar
                if (lfrozen_bot_var_x(iv)) df(l1,m,n,iv) = 0.
              enddo
            endif
            if (ipx == nprocx-1) then
              do iv=1,nvar
                if (lfrozen_top_var_x(iv)) df(l2,m,n,iv) = 0.
              enddo
            endif
          endif
!
        endif
!
!  freezing boundary conditions in y
!
        if (lfrozen_bcs_y) then ! are there any frozen vars at all?
!
! Only need to do this for nonperiodic y direction, on bottom/top-most
! processor and in bottom/top-most pencils
!
          if (.not. lperi(2)) then
            if ((ipy == 0) .and. (m == m1)) then
              do iv=1,nvar
                if (lfrozen_bot_var_y(iv)) df(l1:l2,m,n,iv) = 0.
              enddo
            endif
            if ((ipy == nprocy-1) .and. (m == m2)) then
              do iv=1,nvar
                if (lfrozen_top_var_y(iv)) df(l1:l2,m,n,iv) = 0.
              enddo
            endif
          endif
        endif
!
!  freezing boundary conditions in z
!
        if (lfrozen_bcs_z) then ! are there any frozen vars at all?
!
! Only need to do this for nonperiodic z direction, on bottom/top-most
! processor and in bottom/top-most pencils
!
          if (.not. lperi(3)) then
            if ((ipz == 0) .and. (n == n1)) then
              do iv=1,nvar
                if (lfrozen_bot_var_z(iv)) df(l1:l2,m,n,iv) = 0.
              enddo
            endif
            if ((ipz == nprocz-1) .and. (n == n2)) then
              do iv=1,nvar
                if (lfrozen_top_var_z(iv)) df(l1:l2,m,n,iv) = 0.
              enddo
            endif
          endif
        endif
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
        if (lfirst.and.ldt) then
!
!  sum or maximum of the advection terms?
!  (lmaxadvec_sum=.false. by default)
!
          maxadvec=advec_uu+advec_shear+advec_hall+ &
              sqrt(advec_cs2+advec_va2+advec_crad2)
          maxdiffus=max(diffus_nu,diffus_chi,diffus_eta,diffus_diffrho, &
              diffus_pscalar,diffus_cr,diffus_nud,diffus_diffnd,diffus_chiral)
          if (nxgrid==1.and.nygrid==1.and.nzgrid==1) then
            maxadvec=0.
            maxdiffus=0.
          endif
          dt1_advec=maxadvec/cdt
          dt1_diffus=maxdiffus/cdtv
          dt1_max=max(dt1_max,sqrt(dt1_advec**2+dt1_diffus**2))
!
          if (ldiagnos.and.idiag_dtv/=0) then
            call max_mn_name(maxadvec/cdt,idiag_dtv,l_dt=.true.)
          endif
        endif
!
!  Display derivative info
!
!debug   if (loptimise_ders.and.lout) then                         !DERCOUNT
!debug     do iv=1,nvar                                            !DERCOUNT
!debug     do ider=1,8                                             !DERCOUNT
!debug     do j=1,3                                                !DERCOUNT
!debug     do k=1,3                                                !DERCOUNT
!debug       if (der_call_count(iv,ider,j,k) .gt. 1) then          !DERCOUNT
!debug         print*,'DERCOUNT: '//varname(iv)//' derivative ', & !DERCOUNT
!debug                                                 ider,j,k, & !DERCOUNT
!debug                                               ' called ', & !DERCOUNT
!debug                              der_call_count(iv,ider,j,k), & !DERCOUNT
!debug                                                  'times!'   !DERCOUNT
!debug       endif                                                 !DERCOUNT
!debug     enddo                                                   !DERCOUNT
!debug     enddo                                                   !DERCOUNT
!debug     enddo                                                   !DERCOUNT
!debug     enddo                                                   !DERCOUNT
!debug     if (maxval(der_call_count).gt.1) call fatal_error( &        !DERCOUNT
!debug      'pde','ONE OR MORE DERIVATIVES HAS BEEN DOUBLE CALLED') !DERCOUNT
!debug   endif
!
!  end of loops over m and n
!
        headtt=.false.
     enddo
!
!  Check for NaNs in the advection time-step.
!
     if (notanumber(dt1_advec)) then
       print*, 'pde: dt1_advec contains a NaN at iproc=', iproc
       call fatal_error_local('pde','')
     endif
!
      if (lradiation_fld) f(:,:,:,idd)=DFF_new
!
!  Calculate the gradient of the potential if there is room allocated in the
!  f-array.
!
      if (igpotselfx/=0) then
        call initiate_isendrcv_bdry(f,igpotselfx,igpotselfz)
        call finalize_isendrcv_bdry(f,igpotselfx,igpotselfz)
        call boundconds_x(f,igpotselfx,igpotselfz)
        call boundconds_y(f,igpotselfx,igpotselfz)
        call boundconds_z(f,igpotselfx,igpotselfz)
      endif
!
!  Change dfp according to the chosen particle modules
!
      if (lparticles) call particles_pde(f,df)
!
!  in case of lvisc_hyper=true epsK is calculated for the whole array 
!  at not just for one pencil, it must therefore be added outside the
!  m,n loop.
!      
!ajwm idiag_epsK needs close inspection... and requires tidying up
!ajwm to be consistent in the viscosity.f90 routine.
!      if (lvisc_hyper .and. ldiagnos) fname(idiag_epsK)=epsK_hyper

!
!  diagnostic quantities
!  collect from different processors UUmax for the time step
!
      if (lfirst.and.ldt) call collect_UUmax
      if (ldiagnos) then
        call diagnostic
        call xyaverages_z
        call xzaverages_y
        call yzaverages_x
      endif
!
!  2-D averages
!
      if (l2davgfirst) then
        if (lwrite_yaverages) call yaverages_xz
        if (lwrite_zaverages) call zaverages_xy
        if (lwrite_phiaverages) call phiaverages_rz
      endif
      if (.not.l2davgfirst.and.ldiagnos.and.ldiagnos_need_zaverages) then
!
!  Note: zaverages_xy are also needed if bmx and bmy are to be calculated
!  (Of course, yaverages_xz does not need to be calculated for that.)
!
        if (lwrite_zaverages) call zaverages_xy
      endif
!
!  Force reinitializing of dust variables if certain criteria are fulfilled
!
      call reinit_criteria_dust
!
!  reset lwrite_prof
!
      lwrite_prof=.false.
!
!  Fold df from first ghost zone into main df.
!  Currently only needed for smoothed out particle drag force.
!
      if (lhydro .and. lfold_df) call fold_df(df,iux,iuz)
!
    endsubroutine pde
!***********************************************************************
    subroutine debug_imn_arrays
!
!  for debug purposes: writes out the mm, nn, and necessary arrays
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
    subroutine pencil_consistency_check(f,df,p)
!
!  This subroutine checks the run for missing and for superfluous pencils.
!  First a reference df is calculated with all the requested pencils. Then
!  the pencil request is flipped one by one (a pencil that is requested
!  is not calculated, a pencil that is not requested is calculated). A
!  complete set of pencils should fulfil
!    - Calculating a not requested pencil should not change df
!    - Not calculating a requested pencil should change df
!  The run has a problem when
!    - Calculating a not requested pencil changes df
!      (program dies with error message)
!    - Not calculating a requested pencil does not change df
!      (program gives a warning)
!  If there are missing pencils, the programmer should go into the 
!  pencil_criteria_XXX subs and request the proper pencils (based cleverly
!  on run parameters).
!
!  18-apr-05/tony: coded
!
      use Cdata
      use General, only: random_number_wrapper, random_seed_wrapper
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, allocatable, dimension(:,:,:,:) :: df_ref, f_other
      real, allocatable, dimension(:) :: fname_ref
      real, dimension (nx) :: dt1_max_ref
      integer :: i,j,k,penc,iv
      integer, dimension (mseed) :: iseed_org
      logical :: lconsistent=.true., ldie=.false.
!
      if (lroot) print*, &
          'pencil_consistency_check: checking the pencil case'      
!
! Prevent code from dying due to any errors...
!
      call life_support_on
!
!  Allocate memory for alternative df, fname
!
      allocate(df_ref(mx,my,mz,mvar))
      allocate(fname_ref(mname))
      allocate(f_other(mx,my,mz,mvar+maux))
!
!  Check requested pencils
!
      lfirst=.true.
      headt=.false.
      itsub=1                   ! some modules like dustvelocity.f90
                                ! reference dt_beta(itsub)
      call random_seed_wrapper(get=iseed_org)
      call random_seed_wrapper(put=iseed_org)
      do i=1,mvar+maux
        call random_number_wrapper(f_other(:,:,:,i))
      enddo
      df_ref=0.0
      call initialize_pencils(p,penc0)
!
!  Calculate reference results with all requested pencils on
!
      lpencil=lpenc_requested
      call pde(f_other,df_ref,p)
      dt1_max_ref=dt1_max
!
      do penc=1,npencils 
        df=0.0
        call random_seed_wrapper(put=iseed_org)
        do i=1,mvar+maux
          call random_number_wrapper(f_other(:,:,:,i))
        enddo
        call initialize_pencils(p,penc0) 
!
!  Calculate results with one pencil swapped
!
        lpencil=lpenc_requested
        lpencil(penc)=(.not. lpencil(penc))
        call pde(f_other,df,p)
!
!  Compare results...
!
        lconsistent=.true.
        do i=1,nx
          if (dt1_max(i)/=dt1_max_ref(i)) then
            lconsistent=.false.
            exit
          endif
        enddo
        if (lconsistent) then
f_loop:   do iv=1,mvar
            do k=n1,n2; do j=m1,m2; do i=l1,l2
              lconsistent=(df(i,j,k,iv)==df_ref(i,j,k,iv))
              if (.not. lconsistent) exit f_loop
            enddo; enddo; enddo
          enddo f_loop
        endif
! 
        if (lconsistent .and. lpenc_requested(penc)) then
          if (lroot) print '(a,i4,a)', &
              'pencil_consistency_check: OPTIMISATION POTENTIAL... pencil '// &
              trim(pencil_names(penc))//' (',penc,')', &
              'is requested, but does not appear to be required!'
        elseif ( (.not. lconsistent) .and. (.not. lpenc_requested(penc)) ) then
          if (lroot) print '(a,i4,a)', &
              'pencil_consistency_check: MISSING PENCIL... pencil '// &
              trim(pencil_names(penc))//' (',penc,')', &
              'is not requested, but calculating it changes the results!'
          ldie=.true.
        endif
      enddo
!
!  Check diagnostic pencils
!
      lout=.true.
      lfirst=.true.
      df=0.0
      call random_seed_wrapper(put=iseed_org)
      do i=1,mvar
        call random_number_wrapper(f_other(:,:,:,i))
      enddo
      fname=0.0
      call initialize_pencils(p,penc0) 
!
!  Calculate reference diagnostics with all diagnostic pencils on
!
      lpencil=(lpenc_diagnos.or.lpenc_requested)
      ldiagnos=.true.
      call pde(f_other,df,p)
      fname_ref=fname

      do penc=1,npencils 
        df=0.0
        call random_seed_wrapper(put=iseed_org)
        do i=1,mvar
          call random_number_wrapper(f_other(:,:,:,i))
        enddo
        fname=0.0
        call initialize_pencils(p,penc0) 
!
!  Calculate diagnostics with one pencil swapped
!
        lpencil=(lpenc_diagnos.or.lpenc_requested)
        lpencil(penc)=(.not. lpencil(penc))
        call pde(f_other,df,p)
!
!  Compare results...
!
        lconsistent=.true.
        do k=1,mname
          lconsistent=(fname(k)==fname_ref(k))
          if (.not.lconsistent) exit
        enddo
!        
!  ref = result same as "correct" reference result
!    d = swapped pencil set as diagnostic
!    r = swapped pencil set as requested (can take over for diagnostic pencil) 
!
!   ref +  d +  r = d not needed but set, r not needed but set; optimize d
!   ref +  d + !r = d not needed but set, r not needed and not set; optimize d
!  !ref +  d +  r = d needed and set, r needed and set; d superfluous, but OK
!  !ref +  d + !r = d needed and set; all OK
!   ref + !d +  r = d not needed and not set; all OK
!  !ref + !d +  r = d needed and not set, r needed and set; all OK
!  !ref + !d + !r = d needed and not set, r needed and not set; missing d
!
        if (lpencil_check_diagnos_opti .and. &
            lconsistent .and. lpenc_diagnos(penc) ) then
          if (lroot) print '(a,i4,a)', &
              'pencil_consistency_check: OPTIMISATION POTENTIAL... pencil '// &
              trim(pencil_names(penc))//' (',penc,')', &
              'is requested for diagnostics, '// &
              'but does not appear to be required!'
        elseif ( (.not. lconsistent) .and. (.not. lpenc_diagnos(penc)) .and. &
            (.not. lpenc_requested(penc)) ) then
          if (lroot) print '(a,i4,a)', &
              'pencil_consistency_check: MISSING PENCIL... pencil '// &
              trim(pencil_names(penc))//' (',penc,')', &
              'is not requested for diagnostics, '// &
              'but calculating it changes the diagnostics!'
          ldie=.true.
        endif
      enddo
!
!  Clean up
!
      call random_seed_wrapper(put=iseed_org)
      headt=.true.
      lout=.false.
      lfirst=.false.
      df=0.0      
      deallocate(df_ref)
      deallocate(fname_ref)
      deallocate(f_other)
!
!  Return the code to its former mortal state
!
      call life_support_off
!
      if (ldie) call fatal_error('pencil_consistency_check','DYING')
!
      if (NO_WARN) print*, f(1,1,1,1) !(keep compiler quiet)
!
    endsubroutine pencil_consistency_check
!***********************************************************************

endmodule Equ
