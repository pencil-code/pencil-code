! $Id: dustdensity.f90,v 1.18 2003-12-09 10:35:30 ajohan Exp $

!  This module is used both for the initial condition and during run time.
!  It contains dlnrhod_dt and init_lnrhod, among other auxiliary routines.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Dustdensity

  use Cparam
  use Cdata
  use Dustvelocity

  implicit none

  real, dimension(ndustspec) :: ampllnrhod=0.,amplrhod=0.,cdiffrhod=0.
  real, dimension(ndustspec) :: lnrhod_const=0.,rhod_const=1.
  real, dimension(ndustspec) :: dust_to_gas_ratio=0.
  real, dimension(ndustspec) :: kx_lnrhod,ky_lnrhod,kz_lnrhod
  real :: rhod0=1.,lnrhod0,ampllnrhod_all,dust_to_gas_ratio_all
  real :: kx_lnrhod_all,ky_lnrhod_all,kz_lnrhod_all,amplrhod_all
  real :: rhod_const_all,lnrhod_const_all,cdiffrhod_all
  character (len=labellen), dimension(ndustspec) :: initlnrhod='zero'
  character (len=labellen) :: initlnrhod_all=''

  namelist /dustdensity_init_pars/ &
       rhod0, ampllnrhod, ampllnrhod_all, initlnrhod, initlnrhod_all, &
       dust_to_gas_ratio, dust_to_gas_ratio_all, &
       kx_lnrhod, kx_lnrhod_all, ky_lnrhod, ky_lnrhod_all, &
       kz_lnrhod, kz_lnrhod_all, amplrhod, amplrhod_all, &
       rhod_const, rhod_const_all, lnrhod_const, lnrhod_const_all

  namelist /dustdensity_run_pars/ &
       rhod0, cdiffrhod, cdiffrhod_all

  ! diagnostic variables (needs to be consistent with reset list below)
  integer, dimension(ndustspec) :: i_rhodm=0

  contains

!***********************************************************************
    subroutine register_dustdensity()
!
!  Initialise variables which should know that we solve the
!  compressible hydro equations: ilnrhod; increase nvar accordingly.
!
!   4-jun-02/axel: adapted from hydro
!
      use Mpicomm, only: stop_it
      use Sub
      use General, only: chn
!
      logical, save :: first=.true.
      integer :: idust
      character (len=4) :: sidust
!
      if (.not. first) call stop_it('register_dustdensity: called twice')
      first = .false.
!
      ldustdensity = .true.
!
      do idust=1,ndustspec
        if (idust .eq. 1) then
          ilnrhod(1) = iuud(1)+3         ! indix to access lam
        else
          ilnrhod(idust) = ilnrhod(idust-1)+4
        endif  
        nvar = nvar+1                 ! add 1 variable pr. dust layer
!
        if ((ip<=8) .and. lroot) then
          print*, 'register_dustdensity: idust = ', idust
          print*, 'register_dustdensity: nvar = ', nvar
          print*, 'register_dustdensity: ilnrhod = ', ilnrhod(idust)
        endif
      enddo
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: dustdensity.f90,v 1.18 2003-12-09 10:35:30 ajohan Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_dustdensity: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      do idust=1,ndustspec
        call chn(idust,sidust)
        if (ndustspec .eq. 1) sidust = ''
        if (lroot) then
          if (maux == 0) then
            if (nvar < mvar) write(4,*) ',lnrhod'//trim(sidust)//' $'
            if (nvar == mvar) write(4,*) ',lnrhod'//trim(sidust)
          else
            write(4,*) ',lnrhod'//trim(sidust)//' $'
          endif
          write(15,*) 'lnrhod'//trim(sidust)//' = fltarr(mx,my,mz,1)*one'
        endif
      enddo
!
    endsubroutine register_dustdensity
!***********************************************************************
    subroutine initialize_dustdensity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded 
!  08-dec-03/anders: Copy *_all parameters to whole array
      integer :: idust
!
!  If *_all set, make all empty *(:) = *_all
!
      if (initlnrhod_all .ne. '') then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustdensity: initlnrhod_all=',initlnrhod_all
        do idust=1,ndustspec
          if (initlnrhod(idust) .eq. 'zero') initlnrhod(idust)=initlnrhod_all
        enddo
      endif
!           
      if (dust_to_gas_ratio_all .ne. 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustdensity: dust_to_gas_ratio_all=', &
            dust_to_gas_ratio_all
        do idust=1,ndustspec
          if (dust_to_gas_ratio(idust) .eq. 0.) &
              dust_to_gas_ratio(idust)=dust_to_gas_ratio_all
        enddo
      endif
!           
      if (kx_lnrhod_all .ne. 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustdensity: kx_lnrhod_all=',kx_lnrhod_all
        do idust=1,ndustspec
          if (kx_lnrhod(idust) .eq. 0.) kx_lnrhod(idust)=kx_lnrhod_all
        enddo
      endif
!           
      if (ky_lnrhod_all .ne. 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustdensity: ky_lnrhod_all=',ky_lnrhod_all
        do idust=1,ndustspec
          if (ky_lnrhod(idust) .eq. 0.) ky_lnrhod(idust)=ky_lnrhod_all
        enddo
      endif
!           
      if (kz_lnrhod_all .ne. 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustdensity: kz_lnrhod_all=',kz_lnrhod_all
        do idust=1,ndustspec
          if (kz_lnrhod(idust) .eq. 0.) kz_lnrhod(idust)=kz_lnrhod_all
        enddo
      endif
!           
      if (amplrhod_all .ne. 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustdensity: amplrhod_all=',amplrhod_all
        do idust=1,ndustspec
          if (amplrhod(idust) .eq. 0.) amplrhod(idust)=amplrhod_all
        enddo
      endif
!           
      if (ampllnrhod_all .ne. 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustdensity: ampllnrhod_all=',ampllnrhod_all
        do idust=1,ndustspec
          if (ampllnrhod(idust) .eq. 0.) ampllnrhod(idust)=ampllnrhod_all
        enddo
      endif
!           
      if (rhod_const_all .ne. 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustdensity: rhod_const_all=',rhod_const_all
        do idust=1,ndustspec
          if (rhod_const(idust) .eq. 1.) rhod_const(idust)=rhod_const_all
        enddo
      endif
!           
      if (lnrhod_const_all .ne. 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustdensity: lnrhod_const_all=',lnrhod_const_all
        do idust=1,ndustspec
          if (lnrhod_const(idust) .eq. 0.) lnrhod_const(idust)=lnrhod_const_all
        enddo
      endif
!           
      if (cdiffrhod_all .ne. 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustdensity: cdiffrhod_all=',cdiffrhod_all
        do idust=1,ndustspec
          if (cdiffrhod(idust) .eq. 0.) cdiffrhod(idust)=cdiffrhod_all
        enddo
      endif
!
    endsubroutine initialize_dustdensity
!***********************************************************************
    subroutine init_lnrhod(f,xx,yy,zz)
!
!  initialise lnrhod; called from start.f90
!
!  7-nov-01/wolf: coded
! 28-jun-02/axel: added isothermal
!
      use Mpicomm
      use Sub
      use IO
      use Global
      use Gravity
      use Initcond
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      integer :: i,j,idust
!
!  Loop over dust layers
!
      do idust=1,ndustspec
!
!  different initializations of lnrhod (called from start).
!
        lnrhod0=alog(rhod0)

        select case(initlnrhod(idust))
 
        case('zero'); if(lroot) print*,'init_lnrhod: zero lnrhod'
        case('const_rhod'); f(:,:,:,ilnrhod(idust))=alog(rhod_const(idust))
        case('const_lnrhod'); f(:,:,:,ilnrhod(idust))=lnrhod_const(idust)
        case('frac_of_gas_loc')
          if (dust_to_gas_ratio(idust) .lt. 0.) &
              call stop_it("init_lnrhod: Negative dust_to_gas_ratio!")
          f(:,:,:,ilnrhod(idust))=alog(dust_to_gas_ratio(idust))+f(:,:,:,ilnrho)
        case('frac_of_gas_glo')
          if (dust_to_gas_ratio(idust) .lt. 0.) &
              call stop_it("init_lnrhod: Negative dust_to_gas_ratio!")
          do i=1,mx
            do j=1,my
              f(i,j,:,ilnrhod(idust)) = &
                  alog(dust_to_gas_ratio(idust))+f(4,4,:,ilnrho)
            enddo
          enddo
        case default
!
!  Catch unknown values
!
          if (lroot) print*, 'init_lnrhod: No such value for initlnrhod: ', &
              trim(initlnrhod(idust))
          call stop_it("")

      endselect
!
!  sanity check
!
        if ( notanumber(f(:,:,:,ilnrhod(idust))) ) then
          STOP "init_lnrhod: Imaginary dustdensity values"
        endif
!
!  End loop over dust layers
!
      enddo
!
      if(ip==0) print*,xx,yy,zz ! keep compiler quiet
!
    endsubroutine init_lnrhod
!***********************************************************************
    subroutine dlnrhod_dt(f,df,uud,glnrhod,divud,lnrhod)
!
!  continuity equation
!  calculate dlnrhod/dt = - u.gradlnrhod - divud
!
!   7-jun-02/axel: incoporated from subroutine pde
!
      use Sub
      use Density, only: cs0
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: uud,glnrhod
      real, dimension (nx) :: lnrhod,divud,uglnrhod,glnrhod2
      real, dimension (nx) :: del2lnrhod,rhod
      real :: diffrhod
      integer :: idust
!
      intent(in)  :: f,uud,divud
      intent(out) :: df,glnrhod,lnrhod
!
!  Loop over dust layers
!
      do idust=1,ndustspec
!
!  identify module and boundary conditions
!
        if (headtt.or.ldebug) print*,'dlnrhod_dt: SOLVE dlnrhod_dt'
        if (headtt) call identify_bcs('lnrhod',ilnrhod(idust))
!
!  define lnrhod; calculate dustdensity gradient and avection term
!
        lnrhod=f(l1:l2,m,n,ilnrhod(idust))
        call grad(f,ilnrhod(idust),glnrhod)
        call dot_mn(uud,glnrhod,uglnrhod)
!
!  continuity equation
!
        df(l1:l2,m,n,ilnrhod(idust))=df(l1:l2,m,n,ilnrhod(idust))-uglnrhod-divud
!
!  mass diffusion, in units of dxmin*cs0
!
        if (cdiffrhod(idust) /= 0.) then
          diffrhod=cdiffrhod(idust)*dxmin*cs0
          call del2(f,ilnrhod(idust),del2lnrhod)
          call dot2_mn(glnrhod,glnrhod2)
          df(l1:l2,m,n,ilnrhod(idust)) = &
              df(l1:l2,m,n,ilnrhod(idust))+diffrhod*(del2lnrhod+glnrhod2)
          maxdiffus=amax1(maxdiffus,diffrhod)
        endif
!
        if (ldiagnos) then
          rhod=exp(f(l1:l2,m,n,ilnrhod(idust)))
          if (i_rhodm(idust)/=0) call sum_mn_name(rhod,i_rhodm(idust))
        endif
!
!  End loop over dust layers
!
      enddo
!
    endsubroutine dlnrhod_dt
!***********************************************************************
    subroutine rprint_dustdensity(lreset,lwrite)
!
!  reads and registers print parameters relevant for compressible part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Sub
      use General, only: chn
!
      integer :: iname,idust
      logical :: lreset,lwr
      logical, optional :: lwrite
      character (len=4) :: sidust
!
!  Write information to index.pro that should not be repeated for idust
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite

      if (lwr) then
        write(3,*) 'ndustspec=',ndustspec
        write(3,*) 'nname=',nname
      endif
!
!  Loop over dust layers
!
      do idust=1,ndustspec
        call chn(idust,sidust)
        if (ndustspec .eq. 1) sidust=''
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
        if (lreset) then
          i_rhodm(idust)=0
        endif
!
!  iname runs through all possible names that may be listed in print.in
!
        if(lroot.and.ip<14) print*,'rprint_dustdensity: run through parse list'
        do iname=1,nname
          call parse_name(iname,cname(iname),cform(iname), &
              'rhodm'//trim(sidust),i_rhodm(idust))
        enddo
!
!  write column where which magnetic variable is stored
!
        if (lwr) then
          write(3,*) 'i_rhodm'//trim(sidust)//'=',i_rhodm(idust)
          write(3,*) 'ilnrhod'//trim(sidust)//'=',ilnrhod(idust)
        endif
!
!  End loop over dust layers
!
      enddo
!
    endsubroutine rprint_dustdensity
!***********************************************************************

endmodule Dustdensity
