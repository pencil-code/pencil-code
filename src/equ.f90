! $Id: equ.f90,v 1.231 2005-06-07 21:21:28 brandenb Exp $

module Equ

  use Cdata

  implicit none

  contains

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
      if(lroot) UUmax=fmax(1)
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
      integer :: iname,imax_count,isum_count,nmax_count,nsum_count
      real :: dv
      real, dimension (mname) :: fmax_tmp,fsum_tmp,fmax,fsum
!
!  go through all print names, and sort into communicators
!  corresponding to their type
!
      imax_count=0
      isum_count=0
      do iname=1,nname
        if(itype_name(iname)<0) then
          imax_count=imax_count+1
          fmax_tmp(imax_count)=fname(iname)
        elseif(itype_name(iname)>0) then
          isum_count=isum_count+1
          fsum_tmp(isum_count)=fname(iname)
        endif
      enddo
      nmax_count=imax_count
      nsum_count=isum_count
!
!  communicate over all processors
!
      call mpireduce_max(fmax_tmp,fmax,nmax_count)
      call mpireduce_sum(fsum_tmp,fsum,nsum_count)
!


!
!  the result is present only on the root processor
!
      if(lroot) then
!        fsum=fsum/(nw*ncpus)
!
!  sort back into original array
!  need to take sqare root if |itype|=2
!  (in current version, don't need itype=2 anymore)
!
         imax_count=0
         isum_count=0
         do iname=1,nname
           if(itype_name(iname)<0) then ! max
             imax_count=imax_count+1

             if(itype_name(iname)==ilabel_max)            &
                 fname(iname)=fmax(imax_count)

             if(itype_name(iname)==ilabel_max_sqrt)       &
                 fname(iname)=sqrt(fmax(imax_count))

             if(itype_name(iname)==ilabel_max_dt)         &
                 fname(iname)=fmax(imax_count)

             if(itype_name(iname)==ilabel_max_neg)        &
                 fname(iname)=-fmax(imax_count)

             if(itype_name(iname)==ilabel_max_reciprocal) &
                 fname(iname)=1./fmax(imax_count)

           elseif(itype_name(iname)>0) then ! sum
             isum_count=isum_count+1

             if(itype_name(iname)==ilabel_sum)            &
                 fname(iname)=fsum(isum_count)/(nw*ncpus)

             if(itype_name(iname)==ilabel_sum_sqrt)       &
                 fname(iname)=sqrt(fsum(isum_count)/(nw*ncpus))

             if(itype_name(iname)==ilabel_integrate) then
               dv=1.
               if (nxgrid/=1) dv=dv*dx
               if (nygrid/=1) dv=dv*dy
               if (nzgrid/=1) dv=dv*dz
               fname(iname)=fsum(isum_count)*dv
              endif

              if(itype_name(iname)==ilabel_surf)          &
                  fname(iname)=fsum(isum_count)
           endif

         enddo
         !nmax_count=imax_count
         !nsum_count=isum_count
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
      use Sub
!
      real, dimension (nz,nprocz,mnamez) :: fsumz
!
!  communicate over all processors
!  the result is only present on the root processor
!
      if(nnamez>0) then
        call mpireduce_sum(fnamez,fsumz,nnamez*nz*nprocz)
        if(lroot) fnamez=fsumz/(nx*ny*nprocy)
      endif
!
    endsubroutine xyaverages_z
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
      if (nnamexy>0) then
        call mpireduce_sum(fnamexz,fsumxz,nnamexz*nx*nz*nprocz)
        if(lroot) fnamexz=fsumxz/(ny*nprocy)
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
        if(lroot) fnamexy=fsumxy/(nz*nprocz)
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
      if(nnamerz>0) then
        call mpireduce_sum(fnamerz,fsumrz,nnamerz*nrcyl*(nz+1)*nprocz)
        if(lroot) then
          do i=1,nnamerz
            fnamerz(:,1:nz,:,i)=fsumrz(:,1:nz,:,i)/spread(fsumrz(:,0,:,1),2,nz)
          enddo
        endif
      endif
!
    endsubroutine phiaverages_rz
!***********************************************************************
    subroutine pde(f,df)
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
      use Entropy
      use Magnetic
      use Testfield
      use Radiation
      use Ionization
      use Pscalar
      use Chiral
      use Dustvelocity
      use Dustdensity
      use CosmicRay
      use Boundcond
      use Shear
      use Density
      use Viscosity, only: calc_viscosity,lvisc_first
!
      logical :: early_finalize
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: uij,udij,bij,aij
      real, dimension (nx,3) :: uu,glnrho,bb,jj,JxBr,gshock,del2A,graddivA
      real, dimension (nx,3,ndustspec) :: uud,gnd
      real, dimension (nx,ndustspec) :: divud,ud2
      real, dimension (nx) :: lnrho,divu,u2,rho,rho1
      real, dimension (nx) :: cs2,va2,TT1,cc,cc1,shock
      real, dimension (nx) :: maxadvec,maxdiffus
      integer :: iv
!
!  print statements when they are first executed
!
      headtt = headt .and. lfirst .and. lroot

      if (headtt.or.ldebug) print*,'pde: ENTER'
      if (headtt) call cvs_id( &
           "$Id: equ.f90,v 1.231 2005-06-07 21:21:28 brandenb Exp $")
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
      early_finalize=test_nonblocking.or.lionization.or.lradiation_ray
!
!  Check for dust grain mass interval overflows
!  (should consider having possibility for all modules to fiddle with the
!   f array before boundary conditions are sent)
!
      if (ldustdensity .and. ldustnulling) call null_dust_vars(f)
      if (ldustdensity .and. lmdvar .and. itsub == 1) call redist_mdbins(f)
!
!  Initiate (non-blocking) communication and do boundary conditions.
!  Required order:
!  1. x-boundaries (x-ghost zones will be communicated)
!  2. communication
!  3. y- and z-boundaries
!
      call boundconds_x(f)
      if (ldebug) print*,'pde: bef. initiate_isendrcv_bdry'
      call initiate_isendrcv_bdry(f)
      if (early_finalize) call finalise_isendrcv_bdry(f)
!
!  Calculate ionization degree (needed for thermodynamics)
!  Radiation transport along rays
!
      if (lionization)    call ioncalc(f)
      if (lradiation_ray) call radtransfer(f)
      if (lvisc_shock.or.lvisc_hyper.or.lvisc_smagorinsky) then
        if ((lvisc_first.and.lfirst).or..not.lvisc_first) call calc_viscosity(f)
      endif
!  Turbulence parameters (alpha, scale height, etc.)      
      if (lcalc_turbulence_pars) call calc_turbulence_pars(f)
!
!  set inverse timestep to zero before entering loop over m and n
!
      if (lfirst.and.ldt) dt1_max=0.0
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
        if (necessary(imn)) then  ! make sure all ghost points are set
          if (.not.early_finalize) call finalise_isendrcv_bdry(f)
          call boundconds_y(f)
          call boundconds_z(f)
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
!
        if (l2davgfirst.and.lwrite_phiaverages) then
          call calc_phiavg_general()
          call calc_phiavg_profile()
          call calc_phiavg_unitvects()
        endif
!
!  general phiaverage quantities -- useful for debugging
!
        if (l2davgfirst) then
          call phisum_mn_name_rz(rcyl_mn,i_rcylmphi)
          call phisum_mn_name_rz(phi_mn ,i_phimphi)
          call phisum_mn_name_rz(z_mn   ,i_zmphi)
          call phisum_mn_name_rz(r_mn   ,i_rmphi)
        endif
!
!  For each pencil, accumulate through the different modules
!  advec_XX and diffus_XX, which are essentially the inverse
!  advective and diffusive timestep for that module.
!  (note: advec_cs2 and advec_va2 are inverse _squared_ timesteps)
!  
        advec_uu=0.; advec_shear=0.; advec_hall=0.
        advec_cs2=0.; advec_va2=0.; advec_uud=0;
        diffus_pscalar=0.
        diffus_chiral=0.; diffus_diffrho=0.; diffus_cr=0.
        diffus_eta=0.; diffus_nu=0.; diffus_chi=0.; diffus_nud=0.
!
!  The following is only kept for backwards compatibility.
!  Will be deleted in the future.
!
        if (old_cdtv) then
          dxyz_2 = max(dx_1(l1:l2)**2,dy_1(m)**2,dz_1(n)**2)
        else
          dxyz_2 = dx_1(l1:l2)**2+dy_1(m)**2+dz_1(n)**2
        endif
!
!  Calculate inverse density and magnetic field
!  WD: Also needed with heat conduction, so we better calculate it in all
!  cases. Could alternatively have a switch lrho1known and check for it,
!  or initialise to 1e35.
!
        call calculate_some_vars(f,lnrho,rho,rho1,bb,jj,bij,aij,del2A,graddivA)
!
!  hydro, density, and entropy evolution
!  They all are needed for setting some variables even
!  if their evolution is turned off.
!
        call duu_dt   (f,df,uu,u2,divu,rho,rho1,glnrho,uij,bij,shock,gshock)
        call dlnrho_dt(f,df,uu,divu,lnrho,rho,glnrho,shock,gshock)
!
!  Entropy evolution
!
        call dss_dt(f,df,uu,divu,lnrho,rho,rho1,glnrho,cs2,TT1,shock,gshock,bb,bij)
!
!  Magnetic field evolution
!
        if (lmagnetic) call daa_dt(f,df,uu,uij,rho1,TT1,bb,bij,aij,jj,JxBr,del2A,graddivA,va2,shock,gshock)
!
!  Testfield evolution
!
        if (ltestfield) call daatest_dt(f,df,uu)
!
!  Passive scalar evolution
!
        call dlncc_dt(f,df,uu,rho,glnrho,cc,cc1)
!
!  Dust evolution
!
        call duud_dt (f,df,uu,rho,rho1,glnrho,cs2,JxBr,uud,ud2,divud,udij)
        call dndmd_dt(f,df,rho,rho1,TT1,cs2,cc,cc1,uud,divud,gnd)
!
!  Add gravity, if present
!  Shouldn't we call this one in hydro itself?
!  WD: there is some virtue in calling all of the dXX_dt in equ.f90
!  AB: but it is not really a new dXX_dt, because XX=uu.
!  duu_dt_grav now also takes care of dust velocity
!
        if (lgrav) then
          if (lhydro) call duu_dt_grav(f,df,uu,rho)
        endif
!
!  cosmic ray energy density
!
        if (lcosmicray) call decr_dt(f,df,uu,rho1,divu,bij,bb)
!
!  chirality of left and right handed aminoacids
!
        if (lchiral) call dXY_chiral_dt(f,df,uu)
!
!  Evolution of radiative energy
!
        if (lradiation_fld) call de_dt(f,df,rho1,divu,uu,uij,TT1,gamma)
!
!  Add radiative cooling (for ray method)
!
        if (lradiation_ray.and.lentropy) call radiative_cooling(f,df,lnrho,TT1)
!
!  Add shear if present
!
        if (lshear)                     call shearing(f,df)
!
!  =======================================
!  NO CALLS MODIFYING DF BEYOND THIS POINT
!  =======================================
!
!  Freeze components of variables in boundary slice if specified by boundary
!  condition 'f'
!
        if (lfrozen_bcs_z) then ! are there any frozen vars at all?
          !
          ! Only need to do this for nonperiodic z direction, on bottommost
          ! processor and in bottommost pencils
          !
          if ((.not. lperi(3)) .and. (ipz == 0) .and. (n == n1)) then
            do iv=1,nvar
              if (lfrozen_bot_var_z(iv)) df(l1:l2,m,n,iv) = 0.
              if (lfrozen_top_var_z(iv)) df(l1:l2,m,n,iv) = 0.
            enddo
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
          maxadvec=advec_uu+advec_shear+advec_hall+sqrt(advec_cs2+advec_va2)
          maxdiffus=max(diffus_nu,diffus_chi,diffus_eta,diffus_diffrho, &
                        diffus_pscalar,diffus_cr,diffus_nud,diffus_chiral)
          dt1_advec=maxadvec/cdt
          dt1_diffus=maxdiffus/cdtv

          dt1_max=max(dt1_max,sqrt(dt1_advec**2+dt1_diffus**2))

          if (ldiagnos.and.i_dtv/=0) then
            call max_mn_name(maxadvec/cdt,i_dtv,l_dt=.true.)
            !call max_mn_name(maxdiffus/cdtv,i_dtv,l_dt=.true.)
          endif
        endif
!
!  calculate density diagnostics: mean density
!  Note that p/rho = gamma1*e = cs2/gamma, so e = cs2/(gamma1*gamma).
!
        if (ldiagnos) then
          if (ldensity .or. ldensity_fixed) then
            ! Nothing seems to depend on lhydro here:
            ! if(lhydro) then
            if (i_ekin/=0) call sum_mn_name(.5*rho*u2,i_ekin)
            if (i_ekintot/=0) call integrate_mn_name(.5*rho*u2,i_ekintot)
            if (i_rhom/=0) call sum_mn_name(rho,i_rhom)
            if (i_rhomin/=0) call max_mn_name(-rho,i_rhomin,lneg=.true.)
            if (i_rhomax/=0) call max_mn_name(rho,i_rhomax)
          endif
          !
          !  Mach number, rms and max
          !
          if (i_Marms/=0) call sum_mn_name(u2/cs2,i_Marms,lsqrt=.true.)
          if (i_Mamax/=0) call max_mn_name(u2/cs2,i_Mamax,lsqrt=.true.)
        endif
!
!  end of loops over m and n
!
        headtt=.false.
      enddo
!
      if (lradiation_fld) f(:,:,:,idd)=DFF_new
!
!  in case of lvisc_hyper=true epsK is calculated for the whole array 
!  at not just for one pencil, it must therefore be added outside the
!  m,n loop.
!      
      if (lvisc_hyper .and. ldiagnos) fname(i_epsK)=epsK_hyper
!
!  diagnostic quantities
!  collect from different processors UUmax for the time step
!
      if (lfirst.and.ldt) call collect_UUmax
      if (ldiagnos) then
        call diagnostic
        call xyaverages_z
      endif
!
!  2-D averages
!
      if (l2davgfirst) then
        if (lwrite_yaverages) call yaverages_xz
        if (lwrite_zaverages) call zaverages_xy
        if (lwrite_phiaverages) call phiaverages_rz
      endif
!
!  Note: zaverages_xy are also needed if bmx and bmy are to be calculated
!
      if (.not.l2davgfirst.and.(i_bmx+i_bmy)>0) then
        if (lwrite_yaverages) call yaverages_xz
        if (lwrite_zaverages) call zaverages_xy
      endif
!
    endsubroutine pde
!***********************************************************************
    subroutine debug_imn_arrays
!
!  for debug purposes: writes out the mm, nn, and necessary arrays
!
!  23-nov-02/axel: coded
!
      use Mpicomm
!
      open(1,file=trim(directory)//'/imn_arrays.dat')
      do imn=1,ny*nz
        if(necessary(imn)) write(1,'(a)') '----necessary=.true.----'
        write(1,'(4i6)') imn,mm(imn),nn(imn)
      enddo
      close(1)
!
    endsubroutine debug_imn_arrays
!***********************************************************************
     subroutine calculate_some_vars(f,lnrho,rho,rho1,bb,jj,bij,aij,del2A,graddivA)
!
!   Calculation of some variables used by routines later at time
!
!   06-febr-04/bing: coded
!
       use Magnetic
       use Density

       real, dimension (mx,my,mz,mvar+maux) :: f
       real, dimension (nx) :: lnrho,rho,rho1
       real, dimension (nx,3) :: bb,jj,del2A,graddivA
       real, dimension (nx,3,3) :: bij,aij

       intent(in)  :: f
       intent(out) :: lnrho,rho,rho1,bb,jj,bij,aij,del2A

       if (ldensity .or. ldensity_fixed) then
          call calculate_vars_rho(f,lnrho,rho,rho1)
       else
          lnrho=0.                ! Default for nodensity.f90
          rho=1.
          rho1=1.
       endif

       if (lmagnetic) then
          call calculate_vars_magnetic(f,bb,jj,bij,aij,del2A,graddivA)
       else
          bb=0.                 ! Default for nomagnetic.f90
       endif

     endsubroutine calculate_some_vars
!***********************************************************************
     
endmodule Equ
