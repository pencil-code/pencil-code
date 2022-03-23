!  $Id$
!
!  Initial condition (density, magnetic field, velocity) 
!  for magnetohydrostatical equilibrium in a global accretion
!  disk with an imposed (cylindrically symmetric) sound speed 
!  profile in spherical coordinates. 
!
!  07-may-09/wlad: adapted from noinitial_condition.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
module InitialCondition
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: g0=1, qgshear=1.5
  real :: density_power_law=1.5,temperature_power_law=1.0,plasma_beta=25.
  real :: dustdensity_powerlaw=1.5,edtog=0.01
  real :: rm_int=0.0,rm_ext=impossible
  real :: tm_bot=0.0,tm_top=impossible
  real :: ampluu_cs_factor=0.01
  logical :: lnumerical_mhsequilibrium=.true.
  logical :: lintegrate_potential=.true.
  logical :: lcap_field=.false.,lcap_field_radius=.false.,lcap_field_theta=.false.
  logical :: ladd_noise_propto_cs=.false. 
  logical :: ladd_field=.true.,ladd_field_azimuthal=.true.,ladd_field_vertical=.false.
!
  namelist /initial_condition_pars/ &
      g0,qgshear,density_power_law,temperature_power_law,plasma_beta,&
      lnumerical_mhsequilibrium,lintegrate_potential,&
      rm_int,rm_ext,tm_bot,tm_top,lcap_field_radius,lcap_field_theta, &
      ladd_noise_propto_cs, ampluu_cs_factor, ladd_field,&
      dustdensity_powerlaw,edtog,ladd_field_azimuthal,ladd_field_vertical
!
  real :: ksi=1.
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  07-oct-09/wlad: coded
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
   subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
        if (lmagnetic) then 
          ksi=(1.+plasma_beta)/plasma_beta
          if (ladd_field) then 
            if (ladd_field_azimuthal.and.ladd_field_vertical) then 
              call fatal_error("initialize_initial_condition",&
                   "Both ladd_field_azimuthal and ladd_field_vertical are true. Choose one.")
            else if ((.not.ladd_field_azimuthal).and.(.not.ladd_field_vertical)) then 
              call fatal_error("initialize_initial_condition",&
                   "Both ladd_field_azimuthal and ladd_field_vertical are false. Choose one.")
            endif
          endif
        else
          ksi=1.
        endif
!
        lcap_field=lcap_field_radius.or.lcap_field_theta
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      use Gravity,        only: acceleration
      use Sub,            only: get_radial_distance, power_law
      use FArrayManager,  only: farray_use_global
      use EquationofState,only: gamma

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: rr_sph,rr_cyl,g_r
      real, dimension (mx) :: OOK2,OO2,tmp,H2
      real :: p,q,OOcorot
      integer, pointer :: iglobal_cs2
      integer :: ics2
!
      if (.not.lspherical_coords) call fatal_error("initial_condition_uu",&
           "This method is only for spherical coordinates. Use centrifugal_balance for cylindrical.")
!
!  Set the sound speed
!
      call set_sound_speed(f)
!
!  Get the sound speed
!
      if (llocal_iso) then 
        nullify(iglobal_cs2)
        call farray_use_global('cs2',iglobal_cs2)
        ics2=iglobal_cs2 
      elseif (lentropy) then
        ics2=iss
      elseif (ltemperature) then 
        ics2=iTT
      endif
!
!  Analytical expression that leads to an equilibrium configuration. 
!  Commented out because the numerical equilibrium enforced below 
!  is a lot better when it comes to reproduce what the code actually
!  solves. 
!
      if (.not.lnumerical_mhsequilibrium) then 
!
        if (lcorotational_frame) then 
          OOcorot=rcorot**(-1.5)
        else
          OOcorot=0.
        endif
!
        p=-density_power_law
        q=-temperature_power_law
!
        do m=1,my;do n=1,mz
          call get_radial_distance(rr_sph,rr_cyl)
          if (lgrav) then 
            call acceleration(g_r)
          ! Cylindrical Keplerian velocity: GM/rr_cyl**3
            if (.not.lcylindrical_gravity) then 
              OOK2=max(-g_r/(rr_sph*sinth(m)**3),0.)
            else
              OOK2=max(-g_r/rr_cyl,0.)
            endif
          elseif (lpointmasses) then 
            call power_law(g0,rr_cyl,2*qgshear,OOK2)
          endif

          !H2 defined as in Fromang et al. 2011
          H2=f(:,m,n,ics2)/(gamma*OOK2)
!
!  pressure correction
!
          if (.not.lcylindrical_gravity) then 
            tmp=1 + H2/rr_cyl**2*(ksi*(p+q-2.) + 2.) + q*(1-sinth(m)) 
          else
            tmp=1 + H2/rr_cyl**2*(ksi*(p+q-2.) + 2.)
          endif
          OO2=OOK2*tmp
!
          f(:,m,n,iuz) = f(:,m,n,iuz) + rr_cyl*(sqrt(OO2)-OOcorot)
!
        enddo;enddo
!
      endif
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho 
!  will take care of converting it to linear 
!  density if you use ldensity_nolog
!
!  07-may-09/wlad: coded
!
      use EquationOfState, only: rho0,gamma
      use FArrayManager,   only: farray_use_global
      use Gravity,         only: potential
      use Sub,             only: get_radial_distance
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: rr_sph,rr_cyl
      real, dimension (mx) :: lnrhomid,strat
      real, dimension (mx) :: cs2,tmp1,tmp2
      real    :: p,q
      integer :: ics2
      integer, pointer :: iglobal_cs2
!
      if (ldensity_linearstart) call fatal_error("initial_condition_lnrho",&
           "Switch off ldensity_linearstart. This routine assumes lnrho in initial condition")
!
      p=-density_power_law 
      q=-temperature_power_law
!
      if (llocal_iso) then 
        if (lroot) print*,&
             'initial_condition_lnrho: locally isothermal approximation'
      else
        if (lroot) print*,&
             'initial_condition_lnrho: gamma=',gamma
      endif
      if (lroot) print*,'Radial density stratification with power law=',p
      if (lroot) print*,'Radial temperature stratification with power law=',q
!
!  Get the sound speed globals
!
      if (llocal_iso) then 
        nullify(iglobal_cs2)
        call farray_use_global('cs2',iglobal_cs2)
        ics2=iglobal_cs2 
      elseif (lentropy) then 
        ics2=iss
      elseif (ltemperature) then 
        ics2=iTT
      endif
!
!  Pencilize the density allocation.
!
      do n=1,mz
        do m=1,my
!
!  Midplane density
!
          call get_radial_distance(rr_sph,rr_cyl)
          lnrhomid=log(rho0)+p*log(rr_cyl/r_ref) 
          f(:,m,n,ilnrho) = f(:,m,n,ilnrho)+lnrhomid
!
!  Vertical stratification
!
          cs2=f(:,m,n,ics2)
!
          if (lgrav) then
            call potential(POT=tmp1,RMN=rr_sph)
            call potential(POT=tmp2,RMN=rr_cyl)
          elseif (lpointmasses) then
            tmp1=-g0/rr_sph 
            tmp2=-g0/rr_cyl
          endif
          if (lcylindrical_gravity) then 
            strat=0.
          else
            strat=-gamma*(tmp1-tmp2)/(cs2*ksi)
          endif
          f(:,m,n,ilnrho) = f(:,m,n,ilnrho)+strat
!
        enddo
      enddo
!
!  If the run is non-magnetic, this is the last routine called. 
!  Enforce numerical equilibrium then
!
      if ((.not.lmagnetic).and.lnumerical_mhsequilibrium) &
           call enforce_numerical_equilibrium(f,lhd=.true.)
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_nd(f)
!
!  Initialize logarithmic density. init_lnrho 
!  will take care of converting it to linear 
!  density if you use ldensity_nolog
!
!  20-sep-17/wlad: coded
!
      use EquationOfState, only: rho0
      use Sub, only: get_radial_distance
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: rr_sph,rr_cyl
      real, dimension (mx) :: lnrhomid
      real    :: p
      integer :: k
!
      p=-dustdensity_powerlaw 
      if (lroot) print*,'Radial dust density stratification with power law=',p
!
!  Pencilize the density allocation.
!
      do n=1,mz
        do m=1,my
!
!  Midplane density
!
          call get_radial_distance(rr_sph,rr_cyl)
          lnrhomid=log(edtog*rho0)+p*log(rr_cyl/r_ref)
          do k=1,ndustspec
             f(:,m,n,ind(k)) = f(:,m,n,ind(k))+exp(lnrhomid)
          enddo
!
        enddo
      enddo
!
    endsubroutine initial_condition_nd
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential. Constant plasma 
!  beta magnetic field. 
!
!  07-may-09/wlad: coded
!
      use FArrayManager, only: farray_use_global
      use Sub, only: gij,curl_mn,get_radial_distance
      use Mpicomm, only: mpibcast_real,mpisend_real,mpirecv_real
      use EquationOfState, only: cs0,rho0
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(nx) :: pressure,Bphi,Atheta,BB
      real, dimension(mx) :: rr_sph,rr_cyl
      real, dimension(nx) :: tmp
      real, dimension(0:nx) :: tmp2
      real :: dr,psum
      integer, pointer :: iglobal_cs2
      integer :: i,ics2,irho,iprocx,iprocy,iserial_xy
      real, dimension(ny,nprocx*nprocy) :: procsum
      real, dimension(ny) :: procsum_loc,tmpy

      addfield: if (ladd_field) then 
!
!  Get the sound speed globals 
!
         nullify(iglobal_cs2)
         call farray_use_global('cs2',iglobal_cs2)
         ics2=iglobal_cs2 
!
!  Density is already in linear after init_lnrho, so
!
         irho=ilnrho
!
         azimuthal_field: if (ladd_field_azimuthal) then 
!
           do m=m1,m2
!
             pressure=f(l1:l2,m,npoint,irho)*f(l1:l2,m,npoint,ics2)
!
!  The following line assumes mu0=1
!
             BB = sqrt(2*pressure/plasma_beta)
             if (lcap_field) then 
               call cap_field(BB,Bphi)
             else
               Bphi=BB
             endif
!
!  Bphi = 1/r*d/dr(r*Atheta), so integrate: Atheta=1/r*Int(B*r)dr
!
             call get_radial_distance(rr_sph,rr_cyl)
             !dr=rr_sph(2)-rr_sph(1)
             tmp=Bphi*rr_sph(l1:l2)
!
             iserial_xy=ipx+nprocx*ipy+1
             tmp2(0)=0.
             do i=1,nx
               dr=rr_sph(i+l1-1)-rr_sph(i+l1-2) 
               tmp2(i)=tmp2(i-1)+tmp(i)*dr
             enddo
             procsum_loc(m-m1+1)=tmp2(nx)
           enddo
!
            if (ip<=9) print*,'send',ipx+nprocx*ipy+1,procsum_loc(mpoint)
!
            if (lroot) then 
              procsum(:,1)=procsum_loc
              do iprocx=0,nprocx-1
                do iprocy=0,nprocy-1
                  iserial_xy=iprocx+nprocx*iprocy+1
                  if (iserial_xy/=1) then
                    call mpirecv_real(tmpy,ny,iserial_xy-1,111)
                    procsum(:,iserial_xy)=tmpy
                  endif
                  if (ip<=9) print*,'recv',iserial_xy,procsum(mpoint,iserial_xy)
                enddo
              enddo
            else
              call mpisend_real(procsum_loc,ny,0,111)
            endif
!
            call mpibcast_real(procsum,(/nprocx*nprocy,ny/))
!
            do m=m1,m2
              do n=n1,n2
!
                pressure=f(l1:l2,m,n,irho)*f(l1:l2,m,n,ics2)
!
!  The following line assumes mu0=1
!
                BB = sqrt(2*pressure/plasma_beta)
                call cap_field(BB,Bphi)
!
!  Bphi = 1/r*d/dr(r*Atheta), so integrate: Atheta=1/r*Int(B*r)dr
!
                call get_radial_distance(rr_sph,rr_cyl)
                !dr=rr_sph(2)-rr_sph(1)
                tmp=Bphi*rr_sph(l1:l2)
!
                if (nprocx==1) then 
                  tmp2(0)=0.
                  do i=1,nx
                    dr=rr_sph(i+l1-1)-rr_sph(i+l1-2)
                    tmp2(i)=tmp2(i-1)+tmp(i)*dr
                  enddo
                else
                  if (lfirst_proc_x) then 
                    tmp2(0)=0.
                    do i=1,nx
                      dr=rr_sph(i+l1-1)-rr_sph(i+l1-2) 
                      tmp2(i)=tmp2(i-1)+tmp(i)*dr
                    enddo
                  else 
                    psum=0.
                    do iprocx=0,ipx-1
                      iserial_xy=iprocx+nprocx*ipy+1
                      psum=psum+procsum(m-m1+1,iserial_xy)
                    enddo
                    tmp2(0)=psum
                    do i=1,nx
                      dr=rr_sph(i+l1-1)-rr_sph(i+l1-2) 
                      tmp2(i)=tmp2(i-1)+tmp(i)*dr
                    enddo
                  endif
                endif
                Atheta=tmp2(1:nx)/rr_sph(l1:l2)
                f(l1:l2,m,n,iay)=f(l1:l2,m,n,iay)+Atheta
              enddo
            enddo

         else if (ladd_field_vertical) then 
           do m=m1,m2
!
             !pressure=f(l1:l2,mpoint,npoint,irho)*f(l1:l2,mpoint,npoint,ics2)
             call get_radial_distance(rr_sph,rr_cyl)
             pressure=rho0*cs0**2/rr_cyl(l1:l2)**(density_power_law+temperature_power_law)
!
!  The following line assumes mu0=1
!
             BB = sqrt(2*pressure/plasma_beta)
             if (lcap_field) then 
               call cap_field(BB,Bphi)
             else
               Bphi=BB
             endif
!
!  Bphi = 1/r*d/dr(r*Atheta), so integrate: Atheta=1/r*Int(B*r)dr
!
             call get_radial_distance(rr_sph,rr_cyl)
             !dr=rr_sph(2)-rr_sph(1)
             tmp=Bphi*rr_cyl(l1:l2)
!
             iserial_xy=ipx+nprocx*ipy+1
             tmp2(0)=0.
             do i=1,nx
               dr=rr_cyl(i+l1-1)-rr_cyl(i+l1-2) 
               tmp2(i)=tmp2(i-1)+tmp(i)*dr
             enddo
             procsum_loc(m-m1+1)=tmp2(nx)
           enddo
!
           if (ip<=9) print*,'send',ipx+nprocx*ipy+1,procsum_loc(mpoint)
!
           if (lroot) then 
             procsum(:,1)=procsum_loc
             do iprocx=0,nprocx-1
               do iprocy=0,nprocy-1
                 iserial_xy=iprocx+nprocx*iprocy+1
                 if (iserial_xy/=1) then
                   call mpirecv_real(tmpy,ny,iserial_xy-1,111)
                   procsum(:,iserial_xy)=tmpy
                 endif
                 if (ip<=9) print*,'recv',iserial_xy,procsum(mpoint,iserial_xy)
               enddo
             enddo
           else
             call mpisend_real(procsum_loc,ny,0,111)
           endif
!
           call mpibcast_real(procsum,(/nprocx*nprocy,ny/))
!
           do m=m1,m2
             do n=n1,n2
!
               !pressure=f(l1:l2,mpoint,npoint,irho)*f(l1:l2,mpoint,npoint,ics2)
               call get_radial_distance(rr_sph,rr_cyl)
               pressure=rho0*cs0**2/rr_cyl(l1:l2)**(density_power_law+temperature_power_law)
!
!  The following line assumes mu0=1
!
               BB = sqrt(2*pressure/plasma_beta)
               call cap_field(BB,Bphi)
!
!  Bphi = 1/r*d/dr(r*Atheta), so integrate: Atheta=1/r*Int(B*r)dr
!
               call get_radial_distance(rr_sph,rr_cyl)
               !dr=rr_sph(2)-rr_sph(1)
               tmp=Bphi*rr_cyl(l1:l2)
!
               if (nprocx==1) then 
                 tmp2(0)=0.
                 do i=1,nx
                   dr=rr_cyl(i+l1-1)-rr_cyl(i+l1-2)
                   tmp2(i)=tmp2(i-1)+tmp(i)*dr
                 enddo
               else
                 if (lfirst_proc_x) then 
                   tmp2(0)=0.
                   do i=1,nx
                     dr=rr_cyl(i+l1-1)-rr_cyl(i+l1-2) 
                     tmp2(i)=tmp2(i-1)+tmp(i)*dr
                   enddo
                 else 
                   psum=0.
                   do iprocx=0,ipx-1
                     iserial_xy=iprocx+nprocx*ipy+1
                     psum=psum+procsum(m-m1+1,iserial_xy)
                   enddo
                   tmp2(0)=psum
                   do i=1,nx
                     dr=rr_cyl(i+l1-1)-rr_cyl(i+l1-2) 
                     tmp2(i)=tmp2(i-1)+tmp(i)*dr
                   enddo
                 endif
               endif
               Atheta=tmp2(1:nx)/rr_cyl(l1:l2)
               f(l1:l2,m,n,iaz)=f(l1:l2,m,n,iaz)+Atheta
             enddo
           enddo

         else
!
!  No field configuration chosen
!
           call fatal_error("initial_conditon_aa",&
                "The world is flat, and we never got here.")

         endif azimuthal_field
!
!  All quantities are set. Enforce numerical equilibrium.
!
         if (lnumerical_mhsequilibrium) &
              call enforce_numerical_equilibrium(f,lhd=.false.)
!         
       endif addfield
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine cap_field(Bin,Bout)
!
      use Sub, only: step
!
      real, dimension(nx) :: Bin,Brad,Bout
      real :: width
      integer :: i
!
      if (lcap_field_radius) then
        do i=1,nx
          width=5./dx_1(i-1+l1)
          Brad(i) = Bin(i) * &
               (step(x(i),rm_int,width)-&
                step(x(i),rm_ext,width))
        enddo
      else
        Brad=Bin
      endif
!
      if (lcap_field_theta) then
        width=1./dy_1(m-1+m1)
        Bout = Brad * &
             (step(y(m),tm_bot,width)-&
              step(y(m),tm_top,width))
      else
        Bout=Brad
      endif
!
      endsubroutine cap_field
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
!
      use EquationOfState, only: gamma,gamma_m1,get_cp1,cs20,lnrho0
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: cs2,lnrho
      real :: cp1
!
      call get_cp1(cp1)
!
      do m=m1,m2; do n=n1,n2
        if (ldensity_nolog) then 
          lnrho=log(f(l1:l2,m,n,irho))
        else
          lnrho=f(l1:l2,m,n,ilnrho)
        endif     
!
!  The sound speed is stored in the energy slot
!
        if (lentropy) then 
          cs2=f(l1:l2,m,n,iss)
          if (pretend_lnTT) then
            f(l1:l2,m,n,iss)=log(cs2*cp1/gamma_m1)
          else
            f(l1:l2,m,n,iss)=1./(gamma*cp1)*(log(cs2/cs20)-gamma_m1*(lnrho-lnrho0))
          endif
        elseif (ltemperature) then 
          cs2=f(l1:l2,m,n,iTT)
          f(l1:l2,m,n,iTT)=cs2*cp1/gamma_m1
        endif
!
      enddo;enddo
!
    endsubroutine initial_condition_ss
!***********************************************************************     
    subroutine set_sound_speed(f)
!
!  Set the thermo-related quantities. Illustrates that 
!  the user can define as many internal routines as wanted.
!
!  10-may-09/wlad : moved from initial_condition_lnrho
!
      use FArrayManager,   only: farray_use_global
      use EquationOfState, only: cs20
      use Sub,             only: get_radial_distance
!
      real, dimension (mx,my,mz,mfarray) :: f      
      real, dimension (mx) :: rr_sph,rr_cyl,cs2
      integer, pointer :: iglobal_cs2,iglobal_glnTT
      real    :: q
      integer :: ics2,iglnTT
!
      q=-temperature_power_law
!
!  Get the globals needed to store sound speed and temperature gradient
!
      if (llocal_iso) then 
        nullify(iglobal_cs2)
        call farray_use_global('cs2',iglobal_cs2);ics2=iglobal_cs2
        nullify(iglobal_glnTT)
        call farray_use_global('glnTT',iglobal_glnTT);iglnTT=iglobal_glnTT
      elseif (lentropy) then 
          ics2=iss
      elseif (ltemperature) then 
          ics2=iTT
      endif
!
!  Set the sound speed - a power law in cylindrical radius.
!
      do m=1,my
        do n=1,mz
          call get_radial_distance(rr_sph,rr_cyl)
          cs2=cs20*(rr_cyl/r_ref)**q
!
!  Store cs2 in one of the free slots of the f-array
!
          f(:,m,n,ics2)=cs2
!
!  Set the velocity noise if a fraction of sound speed
!          
          if (ladd_noise_propto_cs) then 
            call gaunoise_vect(ampluu_cs_factor*sqrt(cs2),f,iux,iuz)
          endif
!
!  Same for the temperature gradient
!
          if (llocal_iso) then             
            f(:,m,n,iglnTT  )=q/rr_sph
            f(:,m,n,iglnTT+1)=q/rr_sph*cotth(m)
            f(:,m,n,iglnTT+2)=0.
          endif
!
        enddo
      enddo
!
    endsubroutine set_sound_speed
!***********************************************************************
    subroutine gaunoise_vect(ampl,f,i1,i2)
!
!  Add Gaussian noise (= normally distributed) white noise for variables i1:i2
!
!  23-may-02/axel: coded
!  10-sep-03/axel: result only *added* to whatever f array had before
!
      use General, only: random_number_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i1,i2
!
      real, dimension (mx) :: r,p,tmp,ampl
      integer :: i
!
      intent(in)    :: ampl,i1,i2
      intent(inout) :: f
!
!  set gaussian random noise vector
!
      do i=i1,i2
        if (lroot.and.m==1.and.n==1) print*,'gaunoise_vect: variable i=',i
        if (modulo(i-i1,2)==0) then
          call random_number_wrapper(r)
          call random_number_wrapper(p)
          tmp=sqrt(-2*log(r))*sin(2*pi*p)
        else
          tmp=sqrt(-2*log(r))*cos(2*pi*p)
        endif
        f(:,m,n,i)=f(:,m,n,i)+ampl*tmp
      enddo
!
    endsubroutine gaunoise_vect
!***********************************************************************
    subroutine enforce_numerical_equilibrium(f,lhd)
!
!  This subroutine does exactly what's done in runtime to the 
!  velocity field, in order to ensure precise numerical 
!  magnetohydrostatical equilibrium. 
!
      use Sub, only: grad,get_radial_distance,&
                     gij,curl_mn,gij_etc,cross_mn,multsv_mn
      use FArrayManager, only: farray_use_global
!
      real, dimension (mx,my,mz,mfarray) :: f      
      real, dimension (nx,3,3) :: aij,bij
      real, dimension (nx,3) :: glnTT,grho,glnrho
      real, dimension (nx,3) :: aa,bb,jj,jxb,jxbr
      real, dimension (nx) :: cs2,rr_sph,rr_cyl,rho1
      real, dimension (nx) :: fpres_thermal,fpres_magnetic
      integer, pointer :: iglobal_cs2, iglobal_glnTT
      integer :: j,ics2,iglnTT
      logical :: lhd
!
      if (.not.llocal_iso) call fatal_error("enforce_numerical_equilibrium",&
           "not coded for other than the local isothermal approximation."//&
           "Switch lnumerical_mhsequilibrium=F in start.in")
!
      if (lcylindrical_gravity) call fatal_error("enforce_numerical_equilibrium",&
           "switch lenforce_mhsequilibrium=F in initial_condition_pars")
!
!  Get the temperature globals
!
      nullify(iglobal_cs2)
      call farray_use_global('cs2',iglobal_cs2);ics2=iglobal_cs2
      nullify(iglobal_glnTT)
      call farray_use_global('glnTT',iglobal_glnTT);iglnTT=iglobal_glnTT
!
!  Azimuthal speed that perfectly balances the pressure gradient. 
!
      do m=m1,m2
        do n=n1,n2
          !the initial condition is always in log density,
          !even when using ldensity_nolog
          call grad(f,ilnrho,glnrho)
          cs2=f(l1:l2,m,n,ics2) 
          glnTT(:,1:2)=f(l1:l2,m,n,iglnTT:iglnTT+1)
!
          fpres_thermal=cs2*(glnrho(:,2)+glnTT(:,2))
!
          if (lmagnetic) then
            aa=f(l1:l2,m,n,iax:iaz)         !aa
            call gij(f,iaa,aij,1)           !aij
            call curl_mn(aij,bb,aa)         !bb
            call gij_etc(f,iaa,aa,aij,bij)  !bij
            call curl_mn(bij,jj,bb)         !jj
            call cross_mn(jj,bb,jxb)        !jxb
            rho1=1./f(l1:l2,m,n,irho)       !1/rho
            call multsv_mn(rho1,jxb,jxbr)   !jxb/rho
            fpres_magnetic=-jxbr(:,2)
          else
            fpres_magnetic=0.
          endif
!         
          call get_radial_distance(rr_sph,rr_cyl)
!
          f(l1:l2,m,n,iuz)=f(l1:l2,m,n,iuz)+&
              sqrt(rr_sph*(fpres_thermal+fpres_magnetic)/cotth(m))
!
        enddo
      enddo
!
    endsubroutine enforce_numerical_equilibrium
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat)
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
  endmodule InitialCondition
