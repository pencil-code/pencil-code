! $Id: radiation_exp.f90,v 1.1 2003-06-12 22:50:30 theine Exp $

!!!  NOTE: this routine will perhaps be renamed to radiation_feautrier
!!!  or it may be combined with radiation_ray.

module Radiation

!  Radiation (solves transfer equation along rays)
!  The direction of the ray is given by the vector (lrad,mrad,nrad),
!  and the parameters radx0,rady0,radz0 gives the maximum number of
!  steps of the direction vector in the corresponding direction.

  use Cparam

  implicit none

  real, dimension (mx,my,mz) :: Qrad,Srad,chi,TT
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
       radx,rady,radz,rad2max,output_Qrad

  namelist /radiation_run_pars/ &
       radx,rady,radz,rad2max,output_Qrad,nocooling

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
           "$Id: radiation_exp.f90,v 1.1 2003-06-12 22:50:30 theine Exp $")
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
      real, dimension(nx) :: lnrho,ss,yH,TT_,kappa
!
!  Use the ionization module to calculate temperature
!  At the moment we don't calculate ghost zones (ok for vertical arrays)  
!
      do n=n1,n2
      do m=m1,m2
         lnrho=f(l1:l2,m,n,ilnrho)
         ss=f(l1:l2,m,n,ient)
         call ioncalc(lnrho,ss,yH,TT=TT_,kappa=kappa)
         Srad(l1:l2,m,n)=sigmaSB*TT_**4/pi
         chi(l1:l2,m,n)=kappa*exp(lnrho)
         TT(l1:l2,m,n)=TT_
      enddo
      enddo
!
    endsubroutine source_function
!***********************************************************************
    function mean_intensity(f)
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
      real, dimension(mx,my,mz,mvar) :: f
      real, dimension(mx,my,mz) :: mean_intensity,Iup,Idown
      real :: dtau,emdtau
      integer :: lr,mr,nr
!
!  loop fastest over first index --> could do this for l-vectors
!
      do mr=m1,m2
      do lr=l1,l2
        Iup(lr,mr,n1)=Srad(lr,mr,n1)
        do nr=n1,n2
          dtau=.5*(chi(lr,mr,nr)+chi(lr,mr,nr+1))*dz
          emdtau=exp(-dtau)
          Iup(lr,mr,nr+1)=Iup(lr,mr,nr)*emdtau &
                          +(1.-emdtau)*Srad(lr,mr,nr) &
                          +(emdtau-1+dtau) &
                           *(Srad(lr,mr,nr+1)-Srad(lr,mr,nr))/dtau
        enddo
        Idown(lr,mr,n2)=0
        do nr=n2,n1,-1
          dtau=.5*(chi(lr,mr,nr)+chi(lr,mr,nr-1))*dz
          emdtau=exp(-dtau)
          Idown(lr,mr,nr-1)=Idown(lr,mr,nr)*emdtau &
                            +(1.-emdtau)*Srad(lr,mr,nr) &
                            +(emdtau-1+dtau) &
                             *(Srad(lr,mr,nr-1)-Srad(lr,mr,nr))/dtau
        enddo
      enddo
      enddo
      mean_intensity=.5*(Iup+Idown)
    endfunction mean_intensity
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
      Qrad=-Srad+mean_intensity(f)
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
      real, dimension (mx,my,mz,mvar) :: f,df
      real :: formfactor=0.5
!
!  Add radiative cooling
!
      do n=n1,n2
      do m=m1,m2
         if(.not. nocooling) then
            df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient) &
                              +4.*pi*chi(l1:l2,m,n) &
                               *Qrad(l1:l2,m,n) &
                               /TT(l1:l2,m,n)*formfactor &
                               *exp(-f(l1:l2,m,n,ilnrho))
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
      if(lroot.and.headt) print*,'output_radiation',Qrad(4,4,4)
      if(output_Qrad) write(lun) Qrad,Srad,chi,TT
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
