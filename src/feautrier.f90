! $Id: feautrier.f90,v 1.12 2003-04-05 21:22:33 brandenb Exp $

!!!  NOTE: this routine will perhaps be renamed to radiation_feautrier
!!!  or it may be combined with radiation_ray.

module Radiation

!  Radiation (solves transfer equation along rays)
!  The direction of the ray is given by the vector (lrad,mrad,nrad),
!  and the parameters radx0,rady0,radz0 gives the maximum number of
!  steps of the direction vector in the corresponding direction.

  use Cparam

  implicit none

  real, dimension (mx,my,mz) :: Qrad,Srad,kappa,TT
  logical :: nocooling=.false.,output_Qrad=.false.
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
       radx,rady,radz,rad2max

  namelist /radiation_run_pars/ &
       radx,rady,radz,rad2max,nocooling,output_Qrad

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
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: feautrier.f90,v 1.12 2003-04-05 21:22:33 brandenb Exp $")
!
    endsubroutine register_radiation
!***********************************************************************
    subroutine initialize_radiation()
!
!  nothing to be done here
!
    endsubroutine initialize_radiation
!***********************************************************************
    subroutine source_function(f)
!
!  calculate source function
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
      use Ionization
!
      real, dimension(mx,my,mz,mvar), intent(in) :: f
      real, dimension(nx) :: lnrho,ss,yH,TT_,kappa_
!
!  Use the ionization module to calculate temperature
!  At the moment we don't calculate ghost zones (ok for vertical arrays)  
!
      do n=n1,n2
      do m=m1,m2
         lnrho=f(l1:l2,m,n,ilnrho)
         ss=f(l1:l2,m,n,ient)
         yH=ionfrac(lnrho,ss)
         call ioncalc(lnrho,ss,yH,TT=TT_,kappa=kappa_)
         Srad(l1:l2,m,n)=sigmaSB*TT_**4/pi
         TT(l1:l2,m,n)=TT_
         kappa(l1:l2,m,n)=kappa_
      enddo
      enddo
!
    endsubroutine source_function
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
!
      real, dimension(mx,my,mz,mvar) :: f
      real, dimension(mx,my,mz) :: feautrier
      real, dimension(nz) :: kaprho,tau,Srad_,Prad_
      real, dimension(nz) :: a,b,c
      integer :: lrad,mrad,nrad
      logical :: err=.false.
!
      do lrad=l1,l2
      do mrad=m1,m2
         kaprho=kappa(lrad,mrad,n1:n2)*exp(f(lrad,mrad,n1:n2,ilnrho))
         tau=spline_integral(z,kaprho)
         Srad_=Srad(lrad,mrad,n1:n2)
         !print*,'kappa=',kappa(lrad,mrad,n1:n2)
         !print*,'tau=',tau
         !print*,'Srad=',Srad_
!
         b(1)=1.+2./(tau(2)-tau(1))+2./(tau(2)-tau(1))**2
         c(1)=-2./(tau(2)-tau(1))**2
         a(nz)=0.
         b(nz)=1.
!
         a(2:nz-1)=  -2./(tau(3:nz)-tau(1:nz-2))/(tau(2:nz-1)-tau(1:nz-2))
         b(2:nz-1)=1.+2./(tau(3:nz)-tau(1:nz-2))/(tau(2:nz-1)-tau(1:nz-2)) &
                     +2./(tau(3:nz)-tau(1:nz-2))/(tau(3:nz)-tau(2:nz-1))
         c(2:nz-1)=  -2./(tau(3:nz)-tau(1:nz-2))/(tau(3:nz)-tau(2:nz-1))
!
         call tridag(a,b,c,Srad_,Prad_,err=err)
         if (err) then
            print*,'tau=',tau
            stop
         endif
!
         feautrier(lrad,mrad,n1:n2)=Prad_
         !print*,'Prad',Prad_
      enddo
      enddo
    endfunction feautrier
!***********************************************************************
    subroutine radtransfer(f)
!
!  Integration radiation transfer equation along rays
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
      use Sub
      use Ionization
!
      real, dimension(mx,my,mz,mvar) :: f
!
!  identifier
!
      if(lroot.and.headt) print*,'radtransfer'
!
      call source_function(f)
      Qrad=-Srad
      Qrad=Qrad+feautrier(f)
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
!
      real, dimension (mx,my,mz,mvar) :: f,df
!
!  Add radiative cooling
!
      do n=n1,n2
      do m=m1,m2
         if(.not. nocooling) then
            df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient) &
                              +4.*pi*kappa(l1:l2,m,n) &
                               *Qrad(l1:l2,m,n) &
                               /TT(l1:l2,m,n)
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
!
      integer, intent(in) :: lun
!
!  identifier
!
      if(lroot.and.headt) print*,'output_radiation',Qrad(4,4,4)
      if(output_Qrad) write(lun) Qrad
!
    endsubroutine output_radiation
!***********************************************************************
    subroutine init_rad(f,xx,yy,zz)
!
!  Dummy routine for Flux Limited Diffusion routine
!  initialise radiation; called from start.f90
!
!  15-jul-2002/nils: dummy routine
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
!
      if(ip==0) print*,f,xx,yy,zz !(keep compiler quiet)
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
!  Dummy routine for Flux Limited Diffusion routine
!  reads and registers print parameters relevant for radiative part
!
!  16-jul-02/nils: adapted from rprint_hydro
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

endmodule Radiation
