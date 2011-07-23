! $Id$
!
!  Filters for smoothing, removing trends (like spurious build-up of
!  horizontal momentum) and similar.
!
module Filter
!
  use Cdata
  use Mpicomm
  use Sub
!
  implicit none
!
  private
!
  public :: rmwig, rmwig_xyaverage
!
  contains
!***********************************************************************
    subroutine rmwig(f,df,ivar1,ivar2,awigg,explog)
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
!  each Nyquist frequency fully). This implies three instances of
!  communication, but that seems to be the way to do it.
!    Note: Wolfgang believes that the necessity for using rmwig on lnrho
!  is no longer given, because upwinding solves the problem.
!
!  30-Aug-02/wolf: coded
!  28-jul-03/axel: moved to own module, allowed range ivar1-ivar2
!
      use Boundcond
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      logical, optional :: explog
      integer :: ivar1,ivar2
      real :: awigg
!
!  print identifier
!
      if (lroot.and.ip<14) then
        write(*,'(" ",A,I2," -",I2,A,G12.5)') &
               'rmwig: removing wiggles in variable ', ivar1,ivar2, ', t=', t
      endif
!
!  Check whether we want to smooth the actual variable f, or exp(f)
!  The latter can be useful if the variable is lnrho or lncc.
!
      if (present(explog)) then
        f(:,:,:,ilnrho)=exp(f(:,:,:,ilnrho))
        if (lroot) print*,'rmwig: turn f into exp(f), ilnrho=',ilnrho
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
      call rmwig_1d(f,df,ivar1,ivar2,awigg,1) ! x directiong
      call boundconds_y(f)
      call rmwig_1d(f,df,ivar1,ivar2,awigg,2) ! y direction
      call boundconds_z(f)
      call rmwig_1d(f,df,ivar1,ivar2,awigg,3) ! z direction
!
!  Revert back to f if we have been working on exp(f)
!
      if (present(explog)) then
        f(l1:l2,m1:m2,n1:n2,ilnrho)=log(f(l1:l2,m1:m2,n1:n2,ilnrho))
        if (lroot) print*,'rmwig: turn f back into log(f), ilnrho=',ilnrho
      endif
!
    endsubroutine rmwig
!***********************************************************************
    subroutine rmwig_1d(f,df,ivar1,ivar2,awigg,idir)
!
!  Remove small scale oscillations (`wiggles') in the x direction from
!  the ivar component of f (normally lnrho).
!  df is only used as work array
!
!  30-Aug-02/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Deriv
!-- use Wsnaps
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: tmp
      real :: awigg
      integer :: ivar,ivar1,ivar2,idir
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
          call finalize_isendrcv_bdry(f)
        endif
        do ivar=ivar1,ivar2
          call der6(f,ivar,tmp,idir,IGNOREDX=.true.)
          df(l1:l2,m,n,ivar) = 1./64.*tmp
        enddo
      enddo
!
!  Not necessary to do this in a (cache-efficient) loop updating in
!  timestep, since this routine will be applied only infrequently.
!
      do ivar=ivar1,ivar2
        f(l1:l2,m1:m2,n1:n2,ivar) = &
             f(l1:l2,m1:m2,n1:n2,ivar) + awigg*df(l1:l2,m1:m2,n1:n2,ivar)
      enddo
!
!-- print*,'rmwig_1d: WRITE df (from der6) for testing; idir=',idir
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
      use Boundcond
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: tmp
      logical, optional :: explog
      integer :: ivar
!
      if (lroot.and.ip<14) then
        if (ivar == ilnrho) then
          print*,'rmwig_old: removing wiggles in lnrho, t=',t
        else
          write(*,'(" ",A,I2,A,G12.5)') &
               'rmwig_old: removing wiggles in variable', ivar, ', t=', t
        endif
      endif
!
!  initiate communication and do boundary conditions
!  need to call boundconds, because it also deals with x-boundaries!
!AB: We don't need to communicate all variables though; just lnrho
!WD: I wouldn't care, since this should be applied quite infrequently
!
      if (ldebug) print*,'rmwig_old: bef. initiate_isendrcv_bdry'
      call boundconds(f)
      call initiate_isendrcv_bdry(f)
!
!  Check whether we want to smooth on the actual variable, or on exp(f)
!  The latter can be useful if the variable is lnrho or lncc.
!
      if (present(explog)) then
        f(:,:,:,ivar)=exp(f(:,:,:,ivar))
        if (lroot) print*,'rmwig_old: turn f into exp(f), ivar=',ivar
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
          call finalize_isendrcv_bdry(f)
        endif
        call del6(f,ivar,tmp,ignoredx=.true.)
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
        f(l1:l2,m1:m2,n1:n2,ivar)=log(f(l1:l2,m1:m2,n1:n2,ivar))
        if (lroot) print*,'rmwig_old: turn f back into log(f), ivar=',ivar
      endif
!
    endsubroutine rmwig_old
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
      real, dimension (mx,my,mz,mfarray) :: f
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
        f(l1:l2,m1:m2,n,ivar)=log(exp(f(l1:l2,m1:m2,n,ivar))-del_average)
      enddo
!
!  print identifier
!
      if (lroot) then
        if (ivar == ilnrho) then
          print*,'rmwig_xyaverage: removing wiggles in xyaverage of rho, t=',t,rhom1,rhom2
        else
          write(*,'(" ",A,I3,A,G12.5)') &
          'rmwig_xyaverage: removing wiggles in xyaverage of variable ', ivar, 't=', t
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
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mz) :: xyaver,xyaver_smooth
      real :: del_average
      integer :: ivar
!
!  print identifier
!
      if (lroot) then
        if (ivar == ilnrho) then
          print*,'rmwig_lnxyaverage: removing wiggles in xyaverage of lnrho, t=',t
        else
          write(*,'(" ",A,I3,A,G12.5)') &
          'rmwig_lnxyaverage: removing wiggles in xyaverage of variable ', ivar, 't=', t
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
endmodule Filter
