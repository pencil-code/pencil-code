! $Id: feautrier.f90,v 1.26 2003-06-15 09:28:10 brandenb Exp $

!!!  NOTE: this routine will perhaps be renamed to radiation_feautrier
!!!  or it may be combined with radiation_ray.

module Radiation

!  Radiation (solves transfer equation along rays)
!  The direction of the ray is given by the vector (lrad,mrad,nrad),
!  and the parameters radx0,rady0,radz0 gives the maximum number of
!  steps of the direction vector in the corresponding direction.

  use Cparam

  implicit none

  real, dimension (mx,my,mz) :: Srad,kaprho
  logical :: nocooling=.false.,test_radiation=.false.,output_Qrad=.false.
!
!  default values for one pair of vertical rays
!
  integer :: radx=0,rady=0,radz=1,rad2max=1
!
!  definition of dummy variables for FLD routine
!
  real :: DFF_new=0.  !(dum)
  integer :: i_frms=0,i_fmax=0,i_Erad_rms=0,i_Erad_max=0
  integer :: i_Egas_rms=0,i_Egas_max=0

  namelist /radiation_init_pars/ &
       radx,rady,radz,rad2max,output_Qrad,test_radiation

  namelist /radiation_run_pars/ &
       radx,rady,radz,rad2max,output_Qrad,test_radiation,nocooling

  contains

!***********************************************************************
    subroutine register_radiation()
!
!  initialise radiation flags
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if(.not. first) call stop_it('register_radiation called twice')
      first = .false.
!
      lradiation = .true.
      lradiation_ray = .true.
!
!  set indices for auxiliary variables
!
      iQrad = mvar + naux +1; naux = naux + 1
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_radiation: radiation naux = ', naux
        print*, 'iQrad = ', iQrad
      endif
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: feautrier.f90,v 1.26 2003-06-15 09:28:10 brandenb Exp $")
!
! Check we aren't registering too many auxiliary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('register_feautrier: naux > maux')
      endif
!
    endsubroutine register_radiation
!***********************************************************************
    subroutine initialize_radiation()
!
!  nothing to be done here
!
    endsubroutine initialize_radiation
!***********************************************************************
    subroutine radcalc(f)
!
!  calculate source function and opacity
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
      use Ionization
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(nx) :: lnrho,ss,yH,TT,kappa
!
!  test
!
print*,'radcalc: test_radiation=',test_radiation
      if(test_radiation) then
        if(lroot) print*,'radcalc: put Srad=kaprho=1 (as a test)'
        Srad=1.
        kaprho=1.
        return
      endif
!
!  Use the ionization module to calculate temperature
!  At the moment we don't calculate ghost zones (ok for vertical arrays)  
!
      do n=n1,n2
      do m=m1,m2
         lnrho=f(l1:l2,m,n,ilnrho)
         ss=f(l1:l2,m,n,ient)
         yH=f(l1:l2,m,n,iyH)
         TT=f(l1:l2,m,n,iTT)
         Srad(l1:l2,m,n)=sigmaSB*TT**4/pi
         kappa=.25*exp(lnrho-lnrho_ion_)*(TT_ion_/TT)**1.5 &
               *exp(TT_ion_/TT)*yH*(1.-yH)*kappa0
         kaprho(l1:l2,m,n)=kappa*exp(lnrho)
      enddo
      enddo
!
    endsubroutine radcalc
!***********************************************************************
    function feautrier(f)
!
!  Solves the transfer equation using Feautrier's method
!  At the moment for vertical rays only
!
!  01-apr/tobi: coded
!
      use Cdata
      use General
      use Ionization
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(mx,my,mz) :: feautrier
      real, dimension(nz) :: kaprho,tau,Srad_,Prad_
      real, dimension(nz) :: a,b,c
      integer :: lrad,mrad,nrad
      logical :: err=.false.
!
!  loop fastest over first index --> could do this for l-vectors
!
      do mrad=m1,m2
      do lrad=l1,l2
         kaprho=f(lrad,mrad,n1:n2,ikappa)*exp(f(lrad,mrad,n1:n2,ilnrho))
         tau=spline_integral(z,kaprho)
         Srad_=Srad(lrad,mrad,n1:n2)
!
!  top boundary: P'=P, together with P1-P0=dtau*P'+.5*dtau^2*P" and P"=P-S
!  lower boundary: P=S
!
         b(1)=1.+2./(tau(2)-tau(1))+2./(tau(2)-tau(1))**2
         c(1)=-2./(tau(2)-tau(1))**2
         a(nz)=0.
         b(nz)=1.
!
!  interior, P"=P-S
!
         a(2:nz-1)=  -2./(tau(3:nz)-tau(1:nz-2))/(tau(2:nz-1)-tau(1:nz-2))
         b(2:nz-1)=1.+2./(tau(3:nz)-tau(1:nz-2))/(tau(2:nz-1)-tau(1:nz-2)) &
                     +2./(tau(3:nz)-tau(1:nz-2))/(tau(3:nz)-tau(2:nz-1))
         c(2:nz-1)=  -2./(tau(3:nz)-tau(1:nz-2))/(tau(3:nz)-tau(2:nz-1))
!
!  solve tri-diagonal matrix, and give detailed error output if problems
!
         call tridag(a,b,c,Srad_,Prad_,err=err)
         if (err) then
            print*,'lnrho=',f(lrad,mrad,n1:n1+5,ilnrho),'...', &
                            f(lrad,mrad,n2-5:n2,ilnrho)
            print*,'ss=',f(lrad,mrad,n1:n1+5,ient),'...', &
                         f(lrad,mrad,n2-5:n2,ient)
            print*,'tau=',tau(1:6),'...',tau(n2-n1-5:n2-n1)
            print*,'kappa=',f(lrad,mrad,n1:n1+5,ikappa),'...', &
                            f(lrad,mrad,n2-5:n2,ikappa)
            stop
         endif
!
         feautrier(lrad,mrad,n1:n2)=Prad_
!         print*,'Prad',Prad_
      enddo
      enddo
    endfunction feautrier
!***********************************************************************
    function feautrier_double(f)
!
!  Solves the transfer equation using Feautrier's method
!  At the moment for vertical rays only
!
!  01-apr/tobi: coded
!  11-apr/axel: turned some matrix stuff into double precision
!
      use Cdata
      use General
      use Ionization
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(mx,my,mz) :: feautrier_double
      double precision, dimension(nz) :: kaprho,tau,Srad_,Prad_
      double precision, dimension(nz) :: a,b,c
      integer :: lrad,mrad,nrad
      logical :: err=.false.
!
!  loop fastest over first index --> could do this for l-vectors
!
      do mrad=m1,m2
      do lrad=l1,l2
         kaprho=f(lrad,mrad,n1:n2,ikappa)*exp(f(lrad,mrad,n1:n2,ilnrho))
         tau=spline_integral_double(z,kaprho)
         Srad_=Srad(lrad,mrad,n1:n2)
!
!  top boundary: P'=P, together with P1-P0=dtau*P'+.5*dtau^2*P" and P"=P-S
!  lower boundary: P=S
!
         b(1)=1.+2./(tau(2)-tau(1))+2./(tau(2)-tau(1))**2
         c(1)=-2./(tau(2)-tau(1))**2
         a(nz)=0.
         b(nz)=1.
!
!  interior, P"=P-S
!
         a(2:nz-1)=  -2./(tau(3:nz)-tau(1:nz-2))/(tau(2:nz-1)-tau(1:nz-2))
         b(2:nz-1)=1.+2./(tau(3:nz)-tau(1:nz-2))/(tau(2:nz-1)-tau(1:nz-2)) &
                     +2./(tau(3:nz)-tau(1:nz-2))/(tau(3:nz)-tau(2:nz-1))
         c(2:nz-1)=  -2./(tau(3:nz)-tau(1:nz-2))/(tau(3:nz)-tau(2:nz-1))
!
!  solve tri-diagonal matrix, and give detailed error output if problems
!
         call tridag_double(a,b,c,Srad_,Prad_,err=err)
         if (err) then
            print*,'lnrho=',f(lrad,mrad,n1:n1+5,ilnrho),'...', &
                            f(lrad,mrad,n2-5:n2,ilnrho)
            print*,'ss=',f(lrad,mrad,n1:n1+5,ient),'...', &
                         f(lrad,mrad,n2-5:n2,ient)
            print*,'tau=',tau(1:6),'...',tau(n2-n1-5:n2-n1)
            print*,'kappa=',f(lrad,mrad,n1:n1+5,ikappa),'...', &
                            f(lrad,mrad,n2-5:n2,ikappa)
            stop
         endif
!
         feautrier_double(lrad,mrad,n1:n2)=Prad_
!         print*,'Prad',Prad_
      enddo
      enddo
!
    endfunction feautrier_double
!***********************************************************************
    subroutine radtransfer(f)
!
!  Integration radiation transfer equation along rays
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
      use Sub
!
      real, dimension(mx,my,mz,mvar+maux) :: f
!
!  identifier
!
      if(lroot.and.headt) print*,'radtransfer'
!
      call radcalc(f)
      f(:,:,:,iQrad)=-Srad+feautrier_double(f)
!
    endsubroutine radtransfer
!***********************************************************************
    subroutine radiative_cooling(f,df)
!
!  calculate source function
!
!  25-mar-03/axel+tobi: coded
!
      use Cdata
      use Ionization
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: formfactor=0.5
!
!  Add radiative cooling
!
      do n=n1,n2
      do m=m1,m2
         if(.not. nocooling) then
            df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient) &
                              +4.*pi*f(l1:l2,m,n,ikappa) &
                               *f(l1:l2,m,n,iQrad) &
                               /f(l1:l2,m,n,iTT)*formfactor
         endif
      enddo
      enddo
!
    endsubroutine radiative_cooling
!***********************************************************************
    subroutine output_radiation(lun)
!
!  Optional output of derived quantities along with VAR-file
!  Called from wsnap
!
!   5-apr-03/axel: coded
!
      use Cdata
      use Ionization
!
      integer, intent(in) :: lun
!
!  identifier
!
!!$
!!$print*,'output_Qrad,lun=',output_Qrad,lun,size(Qrad)
!!$
!!$      if(lroot.and.headt) print*,'output_radiation',Qrad(4,4,4)
!!$      if(output_Qrad) write(2) 'help'
!!$      if(output_Qrad) write(90) 'help'
!!$!
    endsubroutine output_radiation
!***********************************************************************
    subroutine init_rad(f,xx,yy,zz)
!
!  initialise radiation; called from start.f90
!  We have an init parameter (initrad) to stear radiation i.c. independently.
!
!   15-jul-2002/nils: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
      real :: nr1,nr2
      integer :: l12
!
      if(ip==0) print*,yy !(keep compiler quiet)
    endsubroutine init_rad
!***********************************************************************
   subroutine de_dt(f,df,rho1,divu,uu,uij,TT1,gamma)
!
!  Dummy routine for Flux Limited Diffusion routine
!
!  15-jul-2002/nils: dummy routine
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: uu
      real, dimension (nx) :: rho1,TT1
      real, dimension (nx,3,3) :: uij
      real, dimension (nx) :: divu
      real :: gamma
!
      if(ip==0) print*,f,df,rho1,divu,uu,uij,TT1,gamma !(keep compiler quiet)
    endsubroutine de_dt
!*******************************************************************
    subroutine rprint_radiation(lreset)
!
!  Writes iQrad, etc, to index.pro file
!  Also: dummy routine for Flux Limited Diffusion routine
!  reads and registers print parameters relevant for radiative part
!
!  16-jul-02/nils: adapted from rprint_hydro
!  14-jun-03/axel: moved iTT to rprint_ionization
!
      use Cdata
      use Sub
!  
      logical :: lreset
!
!  write column where which radiative variable is stored
!
      write(3,*) 'i_frms=',i_frms
      write(3,*) 'i_fmax=',i_fmax
      write(3,*) 'i_Erad_rms=',i_Erad_rms
      write(3,*) 'i_Erad_max=',i_Erad_max
      write(3,*) 'i_Egas_rms=',i_Egas_rms
      write(3,*) 'i_Egas_max=',i_Egas_max
      write(3,*) 'nname=',nname
      write(3,*) 'ie=',ie
      write(3,*) 'ifx=',ifx
      write(3,*) 'ify=',ify
      write(3,*) 'ifz=',ifz
      write(3,*) 'iQrad=',iQrad
!   
      if(ip==0) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_radiation
!***********************************************************************
    subroutine  bc_ee_inflow_x(f,topbot)
!
!  Dummy routine for Flux Limited Diffusion routine
!
!  8-aug-02/nils: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
!
    end subroutine bc_ee_inflow_x
!***********************************************************************
    subroutine  bc_ee_outflow_x(f,topbot)
!
!  Dummy routine for Flux Limited Diffusion routine
!
!  8-aug-02/nils: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
!
    end subroutine bc_ee_outflow_x
!***********************************************************************
    subroutine init_equil(f)
!
!  Routine for calculating equilibrium solution of radiation
!
!  18-jul-02/nils: coded
!
      use Cdata
      use Ionization
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(nx) :: lnrho,ss,source,TT1,cs2,cp1tilde
!
!  Use the ionization module to calculate temperature
!  At the moment we don't calculate ghost zones (ok for vertical arrays)  
!
      do n=n1,n2
      do m=m1,m2
         lnrho=f(l1:l2,m,n,ilnrho)
         ss=f(l1:l2,m,n,ient)
         call thermodynamics(lnrho,ss,cs2,TT1,cp1tilde)
         source=sigmaSB/(pi*TT1**4)
         f(l1:l2,m,n,iQrad) = source
      enddo
      enddo
!
    end subroutine init_equil
!***********************************************************************
endmodule Radiation
