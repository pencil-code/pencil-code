! $Id: equ.f90,v 1.175 2003-11-26 17:34:40 theine Exp $

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
            if(itype_name(iname)<0) then
               imax_count=imax_count+1
               if(itype_name(iname)==-1) fname(iname)=fmax(imax_count)
               if(itype_name(iname)==-2) fname(iname)=sqrt(fmax(imax_count))
               if(itype_name(iname)==-3) fname(iname)=fmax(imax_count)
            elseif(itype_name(iname)>0) then
               isum_count=isum_count+1
               if(itype_name(iname)==+1) fname(iname)=fsum(isum_count)/(nw*ncpus)
               if(itype_name(iname)==+2) fname(iname)=sqrt(fsum(isum_count)/(nw*ncpus))
               if(itype_name(iname)==+3) then
                  dv=1.
                  if (nxgrid/=1) dv=dv*dx
                  if (nygrid/=1) dv=dv*dy
                  if (nzgrid/=1) dv=dv*dz
                  fname(iname)=fsum(isum_count)*dv
               endif
               if(itype_name(iname)==+4) fname(iname)=fsum(isum_count)
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
      use Radiation
      use Ionization
      use Pscalar
      use Dustvelocity
      use Dustdensity
      use CosmicRay
      use Boundcond
      use Shear
      use Density
!
      logical :: early_finalize
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: uij,udij,bij
      real, dimension (nx,3) :: uu,uud,glnrho,glnrhod,bb,gshock
      real, dimension (nx) :: lnrho,lnrhod,divu,divud,u2,ud2,rho,rho1
      real, dimension (nx) :: cs2, TT1, shock, UUtemp
      real :: facdiffus,facdss,facdlnrho
!
!  print statements when they are first executed
!
      headtt = headt .and. lfirst .and. lroot

      if (headtt.or.ldebug) print*,'pde: ENTER'
      if (headtt) call cvs_id( &
           "$Id: equ.f90,v 1.175 2003-11-26 17:34:40 theine Exp $")
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
!  when radiation transfer of global ionization is calculated.
!  This could in principle be avoided (but it not worth it now)
!
      early_finalize=test_nonblocking.or.lionization.or.lradiation_ray
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
      if(lionization) call ioncalc(f)
      if(lradiation_ray) call radtransfer(f)
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
        if (i_rcylmphi/=0) call phisum_mn_name_rz(rcyl_mn,i_rcylmphi)
        if (i_phimphi/=0)  call phisum_mn_name_rz(phi_mn ,i_phimphi)
        if (i_zmphi/=0)    call phisum_mn_name_rz(z_mn   ,i_zmphi)
        if (i_rmphi/=0)    call phisum_mn_name_rz(r_mn   ,i_rmphi)
      endif
!
!  for each pencil, accumulate through the different routines
!  maximum diffusion, maximum advection (keep as nx-array)
!  and maximum heating
!  
        maxdiffus  = 0.
        maxadvec2  = 0.
!
!  calculate inverse density
!  WD: Also needed with heat conduction, so we better calculate it in all
!  cases. Could alternatively have a switch lrho1known and check for it,
!  or initialise to 1e35.
!
        if (ldensity) then
          rho1=exp(-f(l1:l2,m,n,ilnrho))
        else
          rho1=1.               ! for all the modules that use it
        endif
!
!  hydro, density, and entropy evolution
!  They all are needed for setting some variables even
!  if their evolution is turned off.
!
        call duu_dt   (f,df,uu,glnrho,divu,rho1,u2,uij,shock,gshock)
        call dlnrho_dt(f,df,uu,glnrho,divu,lnrho,shock,gshock)
        call dss_dt   (f,df,uu,glnrho,divu,rho1,lnrho,cs2,TT1,shock,gshock)
        call dlncc_dt (f,df,uu,glnrho)
!
!  dust equations
!
        call duud_dt   (f,df,uu,uud,divud,ud2,udij)
        call dlnrhod_dt(f,df,uud,glnrhod,divud,lnrhod)
!
!  Add gravity, if present
!  Shouldn't we call this one in hydro itself?
!  WD: there is some virtue in calling all of the dXX_dt in equ.f90
!  AB: but it is not really a new dXX_dt, because XX=uu.
!  duu_dt_grav now also takes care of dust velocity
!
        if (lhydro.or.ldustvelocity) then
          if(lgrav) call duu_dt_grav(f,df,uu,rho1)
        endif
!
!  Magnetic field evolution
!
        if (lmagnetic) call daa_dt(f,df,uu,rho1,TT1,uij,bij,bb)
!
!  cosmic ray energy density
        if (lcosmicray) call decr_dt(f,df,uu,rho1,divu,bij,bb)
!
!  Evolution of radiative energy
!
        if (lradiation_fld) call de_dt(f,df,rho1,divu,uu,uij,TT1,gamma)
!
!  Add radiative cooling (for ray method)
!
        if (lradiation_ray.and.lentropy) call radiative_cooling(f,df)
!
!  Add shear if present
!
        if (lshear) call shearing(f,df)
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
          !  timestep limited by diffusion
          !
          if (old_cdtv) then
            facdiffus=cdt/(cdtv*dxmin)
          else
            facdiffus=cdt/(cdtv*dxmin)*dimensionality
          endif
          !
          !  sum or maximum of the above terms?
          !  (old_UUmax=.true. by default)
          !
          if (old_UUmax) then
            call max_mn(sqrt(maxadvec2)+(facdiffus*maxdiffus),UUmax)
          else
            UUtemp=amax1(sqrt(maxadvec2),facdiffus*maxdiffus)
            !
            !  timestep limited by entropy variation
            !
            !if (lentropy) then
              !maxdss=maxval(abs(df(l1:l2,m,n,iss)))
              !facdss=cdt*dxmin/cdts
              !UUtemp=amax1(facdss*maxdss,UUtemp)
            !endif
            !
            !  timestep limited by density variation
            !
            !if (ldensity) then
              !maxdlnrho=maxval(abs(df(l1:l2,m,n,ilnrho)))
              !facdlnrho=cdt*dxmin/cdtr
              !UUtemp=amax1(facdlnrho*maxdlnrho,UUtemp)
            !endif
            call max_mn(UUtemp,UUmax)
          endif
        endif
!
!  calculate density diagnostics: mean density
!  Note that p/rho = gamma1*e = cs2/gamma, so e = cs2/(gamma1*gamma).
!
        if (ldiagnos) then
          if (ldensity) then
            ! Nothing seems to depend on lhydro here:
            ! if(lhydro) then
            rho=exp(f(l1:l2,m,n,ilnrho))
            if (i_ekin/=0) call sum_mn_name(.5*rho*u2,i_ekin)
            if (i_ekintot/=0) call integrate_mn_name(.5*rho*u2,i_ekintot)
            if (i_rhom/=0) call sum_mn_name(rho,i_rhom)
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
        if (lwrite_zaverages) call zaverages_xy
        if (lwrite_phiaverages) call phiaverages_rz
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

endmodule Equ
