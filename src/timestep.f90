! $Id: timestep.f90,v 1.38 2006-06-22 08:58:07 mee Exp $

module Timestep

  use Cparam
  use Cdata

  implicit none

  private

  public :: rk_2n, border_profiles, timestep_autopsy
!
!  border_prof_[x-z] could be of size n[x-z], but having the same
!  length as f() (in the given dimension) gives somehow more natural code.
!
  real, dimension(mx) :: border_prof_x=1.0
  real, dimension(my) :: border_prof_y=1.0
  real, dimension(mz) :: border_prof_z=1.0

  contains

!***********************************************************************
    subroutine rk_2n(f,df,p)
!
!  Runge Kutta advance, accurate to order itorder
!  At the moment, itorder can be 1, 2, or 3.
!
!   2-apr-01/axel: coded
!  14-sep-01/axel: moved itorder to cdata
!
      use Mpicomm
      use Cdata
      use Equ
      use Particles_main
      use Interstellar, only: calc_interstellar_snr_damp_int
      use Shear, only: advance_shear
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: ds
      real, dimension(1) :: dt1, dt1_local
      integer :: j
!
!  coefficients for up to order 3
!
      if (itorder==1) then
        alpha=(/ 0., 0., 0. /)
        beta =(/ 1., 0., 0. /)
      elseif (itorder==2) then
        alpha=(/ 0., -1./2., 0. /)
        beta=(/ 1./2.,  1.,  0. /)
      elseif (itorder==3) then
        !alpha=(/0., -2./3., -1./)
        !beta=(/1./3., 1., 1./2./)
        !  use coefficients of Williamson (1980)
        alpha=(/  0. ,  -5./9., -153./128. /)
        beta=(/ 1./3., 15./16.,    8./15.  /)
      else
        if (lroot) print*,'Not implemented: itorder=',itorder
        call mpifinalize
      endif
!
!  dt_beta may be needed in other modules (like Dustdensity) for fixed dt
!
      if (.not. ldt) dt_beta=dt*beta
!
!  Set up df and ds for each time sub
!
      do itsub=1,itorder
        if (itsub==1) then
          lfirst=.true.
          df=0.
          ds=0.
        else
          lfirst=.false.
          df=alpha(itsub)*df  !(could be subsumed into pde, but could be dangerous!)
          ds=alpha(itsub)*ds
        endif
!
!  Set up particle derivative array.
!
        if (lparticles) call particles_timestep_first()
!
!  Change df according to the chosen physics modules
!
        call pde(f,df,p)
!
        ds=ds+1.
!
!  If we are in the first time substep we need to calculate timestep dt.
!  This is done here because it uses UUmax which was calculated in pde.
!  Only do it on the root processor, then broadcast dt to all others.
!
        if (lfirst.and.ldt) then 
          dt1_local=maxval(dt1_max(1:nx))


          !Timestep growth limiter
          if (real(ddt) .gt. 0.) dt1_local=max(dt1_local(1),dt1_last)
          call mpireduce_max(dt1_local,dt1,1)
          if (lroot) dt=1.0/dt1(1)
          !Timestep growth limiter
          if (ddt/=0.) dt1_last=dt1_local(1)/ddt
          call mpibcast_real(dt,1)
        endif

        if (ldt) dt_beta=dt*beta
        if (ip<=6) print*,'TIMESTEP: iproc,dt=',iproc,dt  !(all have same dt?)
!
!
        if (linterstellar) call calc_interstellar_snr_damp_int(dt_beta(itsub))
!
!  Time evolution of grid variables
!  (do this loop in pencils, for cache efficiency)
!
        do j=1,mvar; do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,j)=f(l1:l2,m,n,j)+dt_beta(itsub)*df(l1:l2,m,n,j) &
                        *border_prof_x(l1:l2)*border_prof_y(m)*border_prof_z(n)
        enddo; enddo; enddo
!
!  Time evolution of particle variables
!
        if (lparticles) call particles_timestep_second()
!
!  Advance deltay of the shear
!
        if (lshear) call advance_shear(dt_beta(itsub)*ds)
!
!  Increase time
!
        t=t+dt_beta(itsub)*ds
!
      enddo
!
    endsubroutine rk_2n
!***********************************************************************
    subroutine timestep_autopsy
!
!  After the event, determine where the timestep too short occured 
!  Kinda like playing Cluedo... Just without the dice. 
!
!  25-aug-04/tony: coded
! 
      use Cdata
      use Cparam
      use Mpicomm, only: start_serialize, end_serialize

      real :: dt_local, dt1_max_local, dt1_max_global
      integer :: l
      integer, dimension(1) :: turn
 
      dt1_max_global=1./dt  !Could more accurately calculate this and mpireduce
      dt1_max_local=maxval(dt1_max)
      dt_local=1.0/dt1_max_local

      if (lroot) then
        print*,"-------- General Description of Time Step Failure -----------"
        print*,"  it=",it
        print*,"  t=",t
      endif
!Note: ALL processors will do this.
! Identify the murderer
  
! Procs testify in serial
     call start_serialize
!        if ( dt >= dt_local ) then
          print*,"------------------ START OF CONFESSION (", iproc, ") ----------------------"
!          print*,"  Ok, you got me... I (processor - ", iproc,") did it."
!          print*,"  I handle the calculation for: "
!          maxadvec=advec_uu+advec_shear+advec_hall+sqrt(advec_cs2+advec_va2)
!          maxdiffus=max(diffus_nu,diffus_chi,diffus_eta,diffus_diffrho, &
!              diffus_pscalar,diffus_cr,diffus_nud,diffus_diffnd,diffus_chiral)
!          if (nxgrid==1.and.nygrid==1.and.nzgrid==1) then
!            maxadvec=0.
!            maxdiffus=0.
!          endif
!          dt1_advec=maxadvec/cdt
!          dt1_diffus=maxdiffus/cdtv
!
!          if (nxgrid/=1) print*,"   ",x(l1)," < x < ",x(l2)
!          if (nygrid/=1) print*,"   ",y(m1)," < y < ",y(m2)
!          if (nzgrid/=1) print*,"   ",z(n1)," < z < ",z(n2)

!! In the kitchen?
!          print*,"Can't be more specific about a location(s) than the x grid index (including ghost zone) and coordinate:"
!          do l=1,nx
!            if (dt1_max(l) >= dt1_max_global) then
!               print*,"    f(",l+nghost-1,",?,?,?)"," -> x =",x(l)
!            endif
!          enddo 

! With the lead pipe?
!          if (maxval(sqrt(dt1_advec**2+dt1_diffus**2)) < dt1_max_local) then  
!            print *,"  It appears it is not a CFL advection/diffusion limit."
!            print*,"  Perhaps another limiter eg. the cooling time" 
!          else
!            if (maxval(dt1_advec)>maxval(dt1_diffus)) then
!              print*,"  It appears the dagger was in the form of an advection term."
!              print*,"   Here's the line up, the big guy is the offender:"
!              if (lhydro) &
!                print*,"     Fluid velocity: maxval(advec_uu)        = ", maxval(advec_uu)
!              if (lshear) &
!                print*,"     Shear velocity: maxval(advec_shear)     = ", maxval(advec_shear)
!              if (lmagnetic) &
!                print*,"     Hall effect:    maxval(advec_hall)      = ", maxval(advec_hall)
!              if (leos) &
!                print*,"     Sound speed:    maxval(sqrt(advec_cs2)) = ", sqrt(maxval(advec_cs2))
!              if (lmagnetic) &
!                print*,"     Alfven speed:    maxval(sqrt(advec_va2)) = ", sqrt(maxval(advec_va2))
!            else
!              print*,"  It appears the dagger was in the form of an diffusion term."
!              print*,"   Here's the line up, the big guy is the offender:"
!              if (lhydro) &
!                print*,"     Fluid viscosity: maxval(diffus_nu)               = ", maxval(diffus_nu)
!              if (lentropy) &
!                print*,"     Thermal diffusion: maxval(diffus_chi)            = ", maxval(diffus_chi)
!              if (lmagnetic) &
!                print*,"     Magnetic diffusion: maxval(diffus_eta)           = ", maxval(diffus_eta)
!              if (ldensity) &
!                print*,"     Mass diffusion: maxval(diffus_diffrho)           = ", maxval(diffus_diffrho)
!              if (lcosmicray) &
!                print*,"     Passive scalar diffusion: maxval(diffus_pscalar) = ", maxval(diffus_pscalar)
!              if (lcosmicray) &
!                print*,"     Cosmic ray diffusion: maxval(diffus_cr)          = ", maxval(diffus_cr)
!              if (ldustvelocity) &
!                print*,"     Dust viscosity: maxval(diffus_nud)               = ", maxval(diffus_nud)
!              if (lchiral) &
!                print*,"     Chirality diffusion: maxval(diffus_chiral)       = ", maxval(diffus_chiral)
!            endif
!  
!            print*,"  Also, cdt (advection), cdtv (diffusion) = ",cdt,cdtv
            print*,"     Fluid velocity:       maxval(advec_uu)           = ", maxval(advec_uu)/cdt
            print*,"     Shear velocity:       maxval(advec_shear)        = ", maxval(advec_shear)/cdt
            print*,"     Hall effect:          maxval(advec_hall)         = ", maxval(advec_hall)/cdt
            print*,"     Sound speed:          maxval(sqrt(advec_cs2))    = ", sqrt(maxval(advec_cs2))/cdt
            print*,"     Alfven speed:          maxval(sqrt(advec_va2))    = ", sqrt(maxval(advec_va2))/cdt
            print*,"     Fluid viscosity:      maxval(diffus_nu)          = ", maxval(diffus_nu)/cdtv
            print*,"     Thermal diffusion:    maxval(diffus_chi)         = ", maxval(diffus_chi)/cdtv
            print*,"     Magnetic diffusion:   maxval(diffus_eta)         = ", maxval(diffus_eta)/cdtv
            print*,"     Mass diffusion:       maxval(diffus_diffrho)     = ", maxval(diffus_diffrho)/cdtv
            print*,"     Passive scalar diffusion: maxval(diffus_pscalar) = ", maxval(diffus_pscalar)/cdtv
            print*,"     Cosmic ray diffusion: maxval(diffus_cr)          = ", maxval(diffus_cr)/cdtv
            print*,"     Dust viscosity:       maxval(diffus_nud)         = ", maxval(diffus_nud)/cdtv
            print*,"     Chirality diffusion:  maxval(diffus_chiral)      = ", maxval(diffus_chiral)/cdtv
            print*,"------------------- END OF CONFESSION -----------------------"
  
!          endif
!        endif
     call end_serialize
 
    endsubroutine timestep_autopsy
!***********************************************************************
    subroutine border_profiles()
!
!  Position-dependent quenching factor that multiplies rhs of pde
!  by a factor that goes gradually to zero near the boundaries.
!  border_frac is a 3-D array, separately for all three directions.
!  border_frac=1 would affect everything between center and border.
!
      use Cdata

      real, dimension(nx) :: xi
      real, dimension(ny) :: eta
      real, dimension(nz) :: zeta
      real :: border_width,lborder,uborder
!
!  x-direction
!
      border_prof_x(l1:l2)=1

      if ((border_frac_x(1)>0) .and. (.not. lperi(1))) then
        border_width=border_frac_x(1)*Lxyz(1)/2
        lborder=xyz0(1)+border_width
        xi=1-max(lborder-x(l1:l2),0.0)/border_width
        border_prof_x(l1:l2)=min(border_prof_x(l1:l2),xi**2*(3-2*xi))
      endif

      if ((border_frac_x(2)>0) .and. (.not. lperi(1))) then
        border_width=border_frac_x(2)*Lxyz(1)/2
        uborder=xyz1(1)-border_width
        xi=1-max(x(l1:l2)-uborder,0.0)/border_width
        border_prof_x(l1:l2)=min(border_prof_x(l1:l2),xi**2*(3-2*xi))
      endif
!
!  y-direction
!
      border_prof_y(m1:m2)=1

      if ((border_frac_y(1)>0) .and. (.not. lperi(2))) then
        border_width=border_frac_y(1)*Lxyz(2)/2
        lborder=xyz0(2)+border_width
        eta=1-max(lborder-y(m1:m2),0.0)/border_width
        border_prof_y(m1:m2)=min(border_prof_y(m1:m2),eta**2*(3-2*eta))
      endif

      if ((border_frac_y(2)>0) .and. (.not. lperi(2))) then
        border_width=border_frac_y(2)*Lxyz(2)/2
        uborder=xyz1(2)-border_width
        eta=1-max(y(m1:m2)-uborder,0.0)/border_width
        border_prof_y(m1:m2)=min(border_prof_y(m1:m2),eta**2*(3-2*eta))
      endif
!
!  z-direction
!
      border_prof_z(n1:n2)=1

      if ((border_frac_z(1)>0) .and. (.not. lperi(3))) then
        border_width=border_frac_z(1)*Lxyz(3)/2
        lborder=xyz0(3)+border_width
        zeta=1-max(lborder-z(n1:n2),0.0)/border_width
        border_prof_z(n1:n2)=min(border_prof_z(n1:n2),zeta**2*(3-2*zeta))
      endif

      if ((border_frac_z(2)>0) .and. (.not. lperi(3))) then
        border_width=border_frac_z(2)*Lxyz(3)/2
        uborder=xyz1(3)-border_width
        zeta=1-max(z(n1:n2)-uborder,0.0)/border_width
        border_prof_z(n1:n2)=min(border_prof_z(n1:n2),zeta**2*(3-2*zeta))
      endif
!
    endsubroutine border_profiles
!***********************************************************************

endmodule Timestep
