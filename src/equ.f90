! $Id: equ.f90,v 1.121 2003-02-03 14:49:13 dobler Exp $

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
!
      use Mpicomm
      use Cdata
      use Sub
!
      integer :: iname,imax_count,isum_count,nmax_count,nsum_count
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
        fsum=fsum/(nw*ncpus)
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
        elseif(itype_name(iname)>0) then
          isum_count=isum_count+1
          if(itype_name(iname)==+1) fname(iname)=fsum(isum_count)
          if(itype_name(iname)==+2) fname(iname)=sqrt(fsum(isum_count))
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
      real, dimension (nrcyl,nz,nprocz,mnamerz) :: fsumrz
!
!  communicate over all processors
!  the result is only present on the root processor
!
      if(nnamerz>0) then
        call mpireduce_sum(fnamerz,fsumrz,nnamerz*nrcyl*nz*nprocz)
!        if(lroot) fnamerz=fsumrz/(nx*ny*nprocy)
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
      use Slices
      use Sub
      use Global
      use Hydro
      use Gravity
      use Entropy
      use Magnetic
      use Radiation
      use Pscalar
      use Boundcond
      use IO
      use Shear
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3,3) :: uij
      real, dimension (nx,3) :: uu,glnrho
      real, dimension (nx) :: lnrho,divu,u2,rho,ee=0.,rho1
      real :: fac, facheat
!
!  print statements when they are first executed
!
      headtt = headt .and. lfirst .and. lroot

      if (headtt.or.ldebug) print*,'ENTER: pde'
      if (headtt) call cvs_id( &
           "$Id: equ.f90,v 1.121 2003-02-03 14:49:13 dobler Exp $")
!
!  initialize counter for calculating and communicating print results
!
      ldiagnos=lfirst.and.lout
      if (ldiagnos) tdiagnos=t !(diagnostics are for THIS time)
!
!  Initiate (non-blocking) communication and do boundary conditions.
!  Required order:
!  1. x-boundaries (x-ghost zones will be communicated)
!  2. communication
!  3. y- and z-boundaries
!
      call boundconds_x(f)
      if (ldebug) print*,'PDE: bef. initiate_isendrcv_bdry'
      call initiate_isendrcv_bdry(f)
      if (test_nonblocking) call finalise_isendrcv_bdry(f)
!
!  do loop over y and z
!  set indices and check whether communication must now be completed
!  if test_nonblocking=.true., we communicate immediately as a test.
!
      lfirstpoint=.true.        ! true for very first m-n loop
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        if (necessary(imn)) then ! make sure all ghost points are set
          if (.not.test_nonblocking) call finalise_isendrcv_bdry(f)
          call boundconds_y(f)
          call boundconds_z(f)
        endif
!
!  coordinates are needed frequently
!  --- but not for isotropic turbulence; and there are many other
!  circumstances where this is not needed.
!
        x_mn = x(l1:l2)
        y_mn = spread(y(m),1,nx)
        z_mn = spread(z(n),1,nx)
        rcyl_mn = sqrt(x_mn**2+y_mn**2) ! Needed for phi-averages
        r_mn    = sqrt(x_mn**2+y_mn**2+z_mn**2)
!
!  calculate profile for phi-averages if needed
!
        if (ldiagnos) call calc_phiavg_profile()
!
!  for each pencil, accummulate through the different routines
!  maximum diffusion, maximum advection (keep as nx-array)
!  and maximum heating
!  
        maxdiffus  = 0.
        maxadvec2  = 0.
        maxheating = 0.
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
        call duu_dt   (f,df,uu,glnrho,divu,rho1,u2,uij)
        call dlnrho_dt(f,df,uu,glnrho,divu,lnrho)
        call dss_dt   (f,df,uu,glnrho,divu,rho1,lnrho,cs2,TT1)
        call dlncc_dt (f,df,uu,glnrho)
!
!  Add gravity, if present
!  Shouldn't we call this one in hydro itself?
!  WD: there is some virtue in calling all of the dXX_dt in equ.f90
!  AB: but it is not really a new dXX_dt, because XX=uu.
!
        if (lhydro) then
          if(lgrav) call duu_dt_grav(f,df)
        endif
!
!  Magnetic field evolution
!
        if (lmagnetic) call daa_dt(f,df,uu,rho1,TT1)
!
!  Evolution of radiative energy
!
        if (lradiation) call de_dt(f,df,rho1,divu,uu,uij,TT1,gamma)
!
!  Add shear if precent
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
        if (lfirst.and.ldt) then
          fac=cdt/(cdtv*dxmin)
          facheat=dxmin/cdt

          call max_mn((facheat*maxheating)+ &
               sqrt(maxadvec2)+(fac*maxdiffus),UUmax)
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
            if (gamma1/=0.) ee=cs2/(gamma*gamma1)
            if (i_eth/=0)  call sum_mn_name(rho*ee,i_eth)
            if (i_ekin/=0) call sum_mn_name(.5*rho*u2,i_ekin)
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
        lfirstpoint=.false.
      enddo
      if (lradiation) f(:,:,:,idd)=DFF_new
!
!  diagnostic quantities
!  collect from different processors UUmax for the time step
!
      if (lfirst.and.ldt) call collect_UUmax
      if (ldiagnos) then
        call diagnostic
        call xyaverages_z
        call zaverages_xy
      endif
!
    endsubroutine pde
!***********************************************************************
    subroutine rmwig_xyaverage(f,ivar)
!
!  Removes wiggles from the xyaverage of variable ivar.
!  This routine works currently only on one processor, which
!  may not be too bad an approximation even several procs.
!
!  28-Sep-02/axel: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mz) :: xyaver,xyaver_smooth
      real :: del_average,rhom1,rhom2
      integer :: ivar
!
!  calculate horizontal average and smooth in the vertical direction
!
      do n=1,mz
        xyaver(n)=sum(exp(f(l1:l2,m1:m2,n,ivar)))/(nx*ny)
      enddo
      call smooth_4th(xyaver,xyaver_smooth,.true.)
rhom1=sum(xyaver(n1:n2))/nz
rhom2=sum(xyaver_smooth(n1:n2))/nz
!
      do n=1,mz
        del_average=xyaver(n)-xyaver_smooth(n)
        f(l1:l2,m1:m2,n,ivar)=alog(exp(f(l1:l2,m1:m2,n,ivar))-del_average)
      enddo
!
!  print identifier
!
      if (lroot) then
        if (ivar == ilnrho) then
          print*,'RMWIG: removing wiggles in xyaverage of rho, t=',t,rhom1,rhom2
        else
          write(*,'(" ",A,I3,A,G12.5)') &
          'RMWIG: removing wiggles in xyaverage of variable ', ivar, 't=', t
        endif
      endif
!
    endsubroutine rmwig_xyaverage
!***********************************************************************
    subroutine rmwig_lnxyaverage(f,ivar)
!
!  Removes wiggles from the xyaverage of variable ivar.
!  This routine works currently only on one processor, which
!  may not be too bad an approximation even several procs.
!  This routine operates only on the log of rho, so its not mass conserving
!
!  28-Sep-02/axel: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mz) :: xyaver,xyaver_smooth
      real :: del_average
      integer :: ivar
!
!  print identifier
!
      if (lroot) then
        if (ivar == ilnrho) then
          print*,'RMWIG: removing wiggles in xyaverage of lnrho, t=',t
        else
          write(*,'(" ",A,I3,A,G12.5)') &
          'RMWIG: removing wiggles in xyaverage of variable ', ivar, 't=', t
        endif
      endif
!
!  calculate horizontal average and smooth in the vertical direction
!
      do n=1,mz
        xyaver(n)=sum(f(l1:l2,m1:m2,n,ivar))/(nx*ny)
      enddo
      call smooth_4th(xyaver,xyaver_smooth,.true.)
!
      do n=1,mz
        del_average=xyaver(n)-xyaver_smooth(n)
        f(l1:l2,m1:m2,n,ivar)=f(l1:l2,m1:m2,n,ivar)-del_average
      enddo
!
    endsubroutine rmwig_lnxyaverage
!***********************************************************************
      subroutine smooth_4th(a,b,lbdry)
!
!  28-Sep-02/axel: adapted from f77 version
!
      use Cdata
!
      real, dimension (mz) :: a(mz),b(mz)
      real :: am1,a0,a1,a2,a3,a4
      logical :: lbdry
!
!  smooth a vertical average
!  Note: a and b must be different
!  smoothing on the boundaries only when lbdry=.true.
!
      a2=-1./16.
      a1=1./4.
      a0=5./8.
      do 100 n=3,mz-2
        b(n)=a0*a(n)+a1*(a(n-1)+a(n+1))+a2*(a(n-2)+a(n+2))
100   continue
!
!  next nearest points to the boundaries
!
      a3=1./16.
      a2=-1./4.
      a1=3./8.
      a0=3./4.
      am1=1./16.
      b(   2)=am1*a( 1)+a0*a(   2)+a1*a(   3)+a2*a(   4)+a3*a(   5)
      b(mz-1)=am1*a(mz)+a0*a(mz-1)+a1*a(mz-2)+a2*a(mz-3)+a3*a(mz-4)
!
!  points at the boundaries
!
      if (lbdry) then
        a4=-1./16.
        a3=1./4.
        a2=-3./8.
        a1=1./4.
        a0=15./16.
        b( 1)=a0*a( 1)+a1*a(   2)+a2*a(   3)+a3*a(   4)+a4*a(   5)
        b(mz)=a0*a(mz)+a1*a(mz-1)+a2*a(mz-2)+a3*a(mz-3)+a4*a(mz-4)
      else
        b( 1)=a( 1)
        b(mz)=a(mz)
      endif
!
    endsubroutine smooth_4th
!***********************************************************************
    subroutine rmwig(f,df,ivar,awig,explog)
!
!  Remove small scale oscillations (`wiggles') from a component of f,
!  normally from lnrho. Sometimes necessary since Nyquist oscillations in
!  lnrho do not affect the equation of motion at all (dlnrho=0); thus, in
!  order to keep lnrho smooth one needs to smooth lnrho in sporadic time
!  intervals.
!    Since this is a global operation, we need to do a full loop through
!  imn here (including communication and boundary conditions[?]) and
!  collect the result in df, from where it is applied to f after the
!  loop.
!    This version removes in the three directions consecutively (damping
!  each Nyquist frequency fully). This implies three sets of
!  communication, but that seems to be the way to do it.
!
!  30-Aug-02/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Boundcond
!
      real, dimension (mx,my,mz,mvar) :: f,df
      logical, optional :: explog
      integer :: ivar
      real :: awig
!
      if (lroot) then
        if (ivar == ilnrho) then
          print*,'RMWIG: removing wiggles in lnrho, t=',t
        else
          write(*,'(" ",A,I3,A,G12.5)') &
               'RMWIG: removing wiggles in variable ', ivar, 't=', t
        endif
      endif
!
!  Check whether we want to smooth the actual variable f, or exp(f)
!  The latter can be useful if the variable is lnrho or lncc.
!
      if (present(explog)) then
        f(:,:,:,ivar)=exp(f(:,:,:,ivar))
        if(lroot) print*,'RMWIG: turn f into exp(f), ivar=',ivar
      endif
!
!  Apply boundconds and smooth in all three directions consecutively
!  Note 1: we _can't_ do all boundconds first and then call all rmwig_1d.
!  That's an experimental fact which I [wd] don't exactly understand.
!  Note2: this will clearly cause trouble for explog=T and most boundary
!  conditions, because these are now applied to exp(f), rather than to f.
!  Could mend this by exp-ing and log-ing three times, but am reluctant
!  to do so.
!
!  NB: Need to stick to the order x-y-z, or boundconds will be wrong
!
      call boundconds_x(f)
      call rmwig_1d(f,df,ivar,awig,1) ! x direction
      call boundconds_y(f)
      call rmwig_1d(f,df,ivar,awig,2) ! y direction
      call boundconds_z(f)
      call rmwig_1d(f,df,ivar,awig,3) ! z direction
!
!  Revert back to f if we have been working on exp(f)
!
      if (present(explog)) then
        f(l1:l2,m1:m2,n1:n2,ivar)=alog(f(l1:l2,m1:m2,n1:n2,ivar))
        if(lroot) print*,'RMWIG: turn f back into alog(f), ivar=',ivar
      endif
!
    endsubroutine rmwig
!***********************************************************************
    subroutine rmwig_1d(f,df,ivar,awig,idir)
!
!  Remove small scale oscillations (`wiggles') in the x direction from
!  the ivar component of f (normally lnrho).
!
!  30-Aug-02/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Deriv
!-- use Wsnaps
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx) :: tmp
      real :: awig
      integer :: ivar,idir
!
!  Initiate communication of ghost zones
!  At the moment we communicate all variables, even though we only
!  need to worry about one. The extra overhead is insignificant, because
!  this routine is used only sporadically.
!
      call initiate_isendrcv_bdry(f)
!
!  do loop over y and z
!  set indices and check whether communication must now be completed
!
      lfirstpoint=.true.        ! true for very first m-n loop
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        if (necessary(imn)) then 
          call finalise_isendrcv_bdry(f)
        endif
        call der6(f,ivar,tmp,idir,IGNOREDX=.true.)
        df(l1:l2,m,n,ivar) = 1./64.*tmp
      enddo
!
!  Not necessary to do this in a (cache-efficient) loop updating in
!  timestep, since this routine will be applied only infrequently.
!
      f(l1:l2,m1:m2,n1:n2,ivar) = &
           f(l1:l2,m1:m2,n1:n2,ivar) + awig*df(l1:l2,m1:m2,n1:n2,ivar)
!
!-- print*,'WRITE df (from der6) for testing; idir=',idir
!-- if (idir==1) call wsnap(trim(directory)//'/XX',df,.false.)
!-- if (idir==2) call wsnap(trim(directory)//'/YY',df,.false.)
!-- if (idir==3) call wsnap(trim(directory)//'/ZZ',df,.false.)
!
    endsubroutine rmwig_1d
!***********************************************************************
      subroutine rmwig_old(f,df,ivar,explog)
!
!  Remove small scale oscillations (`wiggles') from a component of f,
!  normally from lnrho. Sometimes necessary since Nyquist oscillations in
!  lnrho do not affect the equation of motion at all (dlnrho=0); thus, in
!  order to keep lnrho smooth one needs to smooth lnrho in sporadic time
!  intervals.
!    Since this is a global operation, we need to do a full loop through
!  imn here (including communication and boundary conditions[?]) and
!  collect the result in df, from where it is applied to f after the
!  loop.
!
!  This version uses an additive operator D_x^6+D_y^6+D_z^6, and fully
!  damps only the diagonal Nyquist wave number. Multiplicative damping in
!  the three directions is almost certainly better; see the new rmwig().
!  The current routine is kept for testing only and should be reomved at
!  some point.
!
!  8-Jul-02/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Boundcond
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx) :: tmp
      logical, optional :: explog
      integer :: ivar
!
      if (lroot) then
        if (ivar == ilnrho) then
          print*,'RMWIG: removing wiggles in lnrho, t=',t
        else
          write(*,'(" ",A,I3,A,G12.5)') &
               'RMWIG: removing wiggles in variable ', ivar, 't=', t
        endif
      endif
!
!  initiate communication and do boundary conditions
!  need to call boundconds, because it also deals with x-boundaries!
!AB: We don't need to communicate all variables though; just lnrho
!WD: I wouldn't care, since this should be applied quite infrequently
!
      if (ldebug) print*,'RMWIG: bef. initiate_isendrcv_bdry'
      call boundconds(f)
      call initiate_isendrcv_bdry(f)
!
!  Check whether we want to smooth on the actual variable, or on exp(f)
!  The latter can be useful if the variable is lnrho or lncc.
!
      if (present(explog)) then
        f(:,:,:,ivar)=exp(f(:,:,:,ivar))
        if(lroot) print*,'RMWIG: turn f into exp(f), ivar=',ivar
      endif
!
!  do loop over y and z
!  set indices and check whether communication must now be completed
!
      lfirstpoint=.true.        ! true for very first m-n loop
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        if (necessary(imn)) then 
          call finalise_isendrcv_bdry(f)
        endif
        call del6_nodx(f,ivar,tmp)
! 1/64 would be fine for 1d runs, but in 3d we have higher wave numbers
!        df(l1:l2,m,n,ivar) = 1./64.*tmp
        df(l1:l2,m,n,ivar) = 1./192.*tmp
      enddo
!
!  Not necessary to do this in a (cache-efficient) loop updating in
!  timestep, since this routine will be applied only infrequently.
!
      f(l1:l2,m1:m2,n1:n2,ivar) = &
           f(l1:l2,m1:m2,n1:n2,ivar) + df(l1:l2,m1:m2,n1:n2,ivar)
!
!  Check whether we want to smooth on the actual variable, or on exp(f)
!  The latter can be useful if the variable is lnrho or lncc.
!
      if (present(explog)) then
        f(l1:l2,m1:m2,n1:n2,ivar)=alog(f(l1:l2,m1:m2,n1:n2,ivar))
        if(lroot) print*,'RMWIG: turn f back into alog(f), ivar=',ivar
      endif
!
    endsubroutine rmwig_old
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
