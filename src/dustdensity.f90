! $Id: dustdensity.f90,v 1.16 2003-12-06 13:52:21 ajohan Exp $

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

  real, dimension(dustlayers) :: rhod0=1.,lnrhod0
  real, dimension(dustlayers) :: ampllnrhod=0.,amplrhod=0.,cdiffrhod=0.
  real, dimension(dustlayers) :: lnrhod_const=0.,rhod_const=1.
  real, dimension(dustlayers) :: dust_to_gas_ratio=0.
  real, dimension(dustlayers) :: kx_lnrhod,ky_lnrhod,kz_lnrhod
  character (len=labellen), dimension(dustlayers) :: initlnrhod='zero'

  namelist /dustdensity_init_pars/ &
       rhod0,ampllnrhod,initlnrhod,dust_to_gas_ratio, &
       kx_lnrhod,ky_lnrhod,kz_lnrhod,amplrhod,rhod_const

  namelist /dustdensity_run_pars/ &
       rhod0,cdiffrhod

  ! diagnostic variables (needs to be consistent with reset list below)
  integer, dimension(dustlayers) :: i_rhodm=0

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
      do idust=1,dustlayers
        if (idust .eq. 1) then
          ilnrhod(1) = nvar+1         ! indix to access lam
        else
          ilnrhod(idust) = ilnrhod(idust-1)+1
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
           "$Id: dustdensity.f90,v 1.16 2003-12-06 13:52:21 ajohan Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_dustdensity: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      do idust=1,dustlayers
        call chn(idust,sidust)
        if (dustlayers .eq. 1) sidust = ''
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
!
!  do nothing
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
      do idust=1,dustlayers
!
!  different initializations of lnrhod (called from start).
!
        lnrhod0(idust)=alog(rhod0(idust))

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
      do idust=1,dustlayers
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
        write(3,*) 'dustlayers=',dustlayers
        write(3,*) 'nname=',nname
      endif
!
!  Loop over dust layers
!
      do idust=1,dustlayers
        call chn(idust,sidust)
        if (dustlayers .eq. 1) sidust=''
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
