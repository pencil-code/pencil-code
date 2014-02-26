! $Id: aerosol_init.f90 21410 2013-12-25 15:44:43Z Nbabkovskaia 
!
!  This module provide a way for users to specify custom initial
!  conditions.
!
!  The module provides a set of standard hooks into the Pencil Code
!  and currently allows the following customizations:
!
!   Description                               | Relevant function call
!  ------------------------------------------------------------------------
!   Initial condition registration            | register_initial_condition
!     (pre parameter read)                    |
!   Initial condition initialization          | initialize_initial_condition
!     (post parameter read)                   |
!                                             |
!   Initial condition for momentum            | initial_condition_uu
!   Initial condition for density             | initial_condition_lnrho
!   Initial condition for entropy             | initial_condition_ss
!   Initial condition for magnetic potential  | initial_condition_aa
!
!   And a similar subroutine for each module with an "init_XXX" call.
!   The subroutines are organized IN THE SAME ORDER THAT THEY ARE CALLED.
!   First uu, then lnrho, then ss, then aa, and so on.
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
  use EquationOfState
!
  implicit none
!
  include '../initial_condition.h'
!
     real :: init_ux=impossible,init_uy=impossible,init_uz=impossible
     integer :: imass=1, spot_number=10
     integer :: index_H2O=3
     integer :: index_N2=4, i_point=1
     integer :: Ndata=10
     real :: dYw=1.,dYw1=1.,dYw2=1., init_water1=0., init_water2=0.
     real :: init_x1=0.,init_x2=0.,init_TT1, init_TT2
     real, dimension(nchemspec) :: init_Yk_1, init_Yk_2
     real :: X_wind=impossible, spot_size=10.
     real :: AA=0.66e-4, d0=2.4e-6 , BB0=1.5*1e-16
     real :: dsize_min=0., dsize_max=0., r0=0., r02=0.,  Period=2. 
     real, dimension(ndustspec) :: dsize, dsize0
     real, dimension(2000) :: Ntot_data
     logical :: lreinit_water=.false.,lwet_spots=.false.
     logical :: linit_temperature=.false., lcurved_xz=.false.
     logical :: ltanh_prof_xy=.false.,ltanh_prof_xz=.false.
     logical :: llog_distribution=.true., lcurved_xy=.false.
     logical :: lACTOS=.false.,lACTOS_read=.true., lACTOS_write=.true., lsinhron=.false.

!
    namelist /initial_condition_pars/ &
     init_ux, init_uy,init_uz,init_x1,init_x2, init_water1, init_water2, &
     lreinit_water, dYw,dYw1, dYw2, X_wind, spot_number, spot_size, lwet_spots, &
     linit_temperature, init_TT1, init_TT2, dsize_min, dsize_max, r0, r02, d0, lcurved_xz, lcurved_xy, &
     ltanh_prof_xz,ltanh_prof_xy, Period, BB0, index_N2, index_H2O, lACTOS, lACTOS_read, lACTOS_write, &
     i_point,Ndata, lsinhron
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  07-may-09/wlad: coded
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
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer :: i,j
      real :: del=10.,bs
!
        if ((init_ux /=impossible) .and. (nygrid>1)) then
         do i=1,my
           f(:,i,:,iux)=cos(Period*PI*y(i)/Lxyz(2))*init_ux
         enddo
        endif
        if ((init_uz /=impossible) .and. (nzgrid>1)) then
         bs=9.81e2*(293.-290.)/293.
         do i=1,mx
            f(i,:,:,iuz)=-sqrt(Lxyz(1)*bs) &
                        *((exp(2.*x(i)/Lxyz(1))+exp(-2.*x(i)/Lxyz(1)))/2.)**(-2)
         enddo
        endif
!
        if ((init_uy /=impossible) .and. (X_wind /= impossible)) then
          do j=1,mx
             f(j,:,:,iuy)=f(j,:,:,iuy) &
              +(init_uy+0.)*0.5+((init_uy-0.)*0.5)  &
              *(exp((x(j)+X_wind)/del)-exp(-(x(j)+X_wind)/del)) &
              /(exp((x(j)+X_wind)/del)+exp(-(x(j)+X_wind)/del))
          enddo
        endif
!
        if ((init_uz /=impossible) .and. (X_wind /= impossible)) then
          do j=1,mx
             f(j,:,:,iuz)=f(j,:,:,iuz) &
              +(init_uz+0.)*0.5+((init_uz-0.)*0.5)  &
              *(exp((x(j)+X_wind)/del)-exp(-(x(j)+X_wind)/del)) &
              /(exp((x(j)+X_wind)/del)+exp(-(x(j)+X_wind)/del))
          enddo
         endif
!
        if ((init_uy /=impossible) .and. (X_wind == impossible)) then
          f(:,:,:,iuy)=f(:,:,:,iuy)+init_uy
        endif
        if ((init_uz /=impossible) .and. (X_wind == impossible)) then
          f(:,:,:,iuz)=f(:,:,:,iuz)+init_uz
        endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine initial_condition_chemistry(f)
!
!  Initialize chemistry.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (my,mz) :: init_water1,init_water2, &
               init_water1_tmp,init_water2_tmp
      !real, dimension (mx,my,mz), intent(inout) :: f
      real, dimension (ndustspec) ::  lnds
      real :: ddsize, ddsize0, del, air_mass, PP
      integer :: i, ii
      logical ::  lstart1=.false., lstart2=.false.
!
      if (llog_distribution) then
        ddsize=(alog(dsize_max)-alog(dsize_min))/(max(ndustspec,2)-1)
      else
        ddsize=(dsize_max-dsize_min)/(max(ndustspec,2)-1) 
      endif
!
      do i=0,(ndustspec-1)
        if (llog_distribution) then
          lnds(i+1)=alog(dsize_min)+i*ddsize
          dsize(i+1)=exp(lnds(i+1))
        else
          lnds(i+1)=dsize_min+i*ddsize
          dsize(i+1)=lnds(i+1)
        endif
      enddo
!
      if (lACTOS) then
       call ACTOS_data(f)
      else
       call air_field_local(f, air_mass, PP)
      endif
!
      call reinitialization(f, air_mass, PP)
!
    endsubroutine initial_condition_chemistry
!***********************************************************************
    subroutine initial_condition_uud(f)
!
!  Initialize dust fluid velocity.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uud
!***********************************************************************
    subroutine initial_condition_nd(f)
!
!  Initialize dust fluid density.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (1300,6) :: coeff_loc
      real, dimension (mx,nchemspec) :: init_distr_loc
      real, dimension (6) :: ctmp
      integer :: i,k,ll1
      real :: r0, delta, tmp
!
      open(143,file="coeff_part.dat")
        do i=1,1300
          read(143,'(f15.6,f15.6,f15.6,f15.6,f15.6,f15.6)'),ctmp
          coeff_loc(i,:)=ctmp
        enddo
      close(143)
!          
       r0=2e-5
       delta=1.2
! 
       ll1=int((x(l1)-xyz0(1))/dx)
!
       do i=l1,l2
       do k=1,ndustspec
         f(i,:,:,ind(k)) = Ntot_data(ll1+i-3)/(2.*pi)**0.5/dsize(k)/alog(delta) &
             * exp(-(alog(2.*dsize(k))-alog(2.*r0))**2/(2.*(alog(delta))**2))  &
             /exp(f(i,:,:,ilnrho))
! 
         tmp=dsize(k)*1e4
          init_distr_loc(i,k)=(coeff_loc(ll1+i-3,1) &
                    *exp(-0.5*( (tmp-coeff_loc(ll1+i-3,2)) /coeff_loc(ll1+i-3,3) )**2)  &
                    +coeff_loc(ll1+i-3,4)+coeff_loc(ll1+i-3,5)*tmp+coeff_loc(ll1+i-3,6)*tmp**2)/1e-4                !1e-4=dsize
          if (init_distr_loc(i,k) .le. 1e-10) init_distr_loc(i,k)=1e-10
!
          f(i,:,:,ind(k)) = f(i,:,:,ind(k)) + init_distr_loc(i,k)/exp(f(i,:,:,ilnrho))
        enddo
!                     
        enddo
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_nd
!***********************************************************************
    subroutine initial_condition_uun(f)
!
!  Initialize neutral fluid velocity.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uun
!***********************************************************************
    subroutine read_initial_condition_pars(unit,iostat)
!
!  07-may-09/wlad: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

!
      if (present(iostat)) then
        read(unit,NML=initial_condition_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=initial_condition_pars,ERR=99)
      endif
!
      99    return
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
!  07-may-09/wlad: coded
!
      integer, intent(in) :: unit
!
     write(unit,NML=initial_condition_pars)
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
    subroutine air_field_local(f, air_mass, PP)
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: sum_Y, tmp
!
      logical :: emptyfile=.true.
      logical :: found_specie
      integer :: file_id=123, ind_glob, ind_chem
      character (len=80) :: ChemInpLine
      character (len=10) :: specie_string
      character (len=1)  :: tmp_string
      integer :: i,j,k=1,index_YY, j1,j2,j3, iter
      real :: YY_k, air_mass, TT=300.
      real, intent(out) :: PP ! (in dynes = 1atm)
      real, dimension(nchemspec)    :: stor2
      integer, dimension(nchemspec) :: stor1
!
      integer :: StartInd,StopInd,StartInd_1,StopInd_1
      integer :: iostat, i1,i2,i3
!
      air_mass=0.
      StartInd_1=1; StopInd_1 =0
      open(file_id,file="air.dat")
!
      if (lroot) print*, 'the following parameters and '//&
          'species are found in air.dat (volume fraction fraction in %): '
!
      dataloop: do
!
        read(file_id,'(80A)',IOSTAT=iostat) ChemInpLine
        if (iostat < 0) exit dataloop
        emptyFile=.false.
        StartInd_1=1; StopInd_1=0
        StopInd_1=index(ChemInpLine,' ')
        specie_string=trim(ChemInpLine(1:StopInd_1-1))
        tmp_string=trim(ChemInpLine(1:1))
!
        if (tmp_string == '!' .or. tmp_string == ' ') then
        elseif (tmp_string == 'T') then
          StartInd=1; StopInd =0
!
          StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
          StartInd=verify(ChemInpLine(StopInd:),' ')+StopInd-1
          StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
!
          read (unit=ChemInpLine(StartInd:StopInd),fmt='(E14.7)') TT
          if (lroot) print*, ' Temperature, K   ', TT
!
        elseif (tmp_string == 'P') then
!
          StartInd=1; StopInd =0
!
          StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
          StartInd=verify(ChemInpLine(StopInd:),' ')+StopInd-1
          StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
!
          read (unit=ChemInpLine(StartInd:StopInd),fmt='(E14.7)') PP
          if (lroot) print*, ' Pressure, Pa   ', PP
!
        else
!
          call find_species_index(specie_string,ind_glob,ind_chem,found_specie)
!
          if (found_specie) then
!
            StartInd=1; StopInd =0
!
            StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
            StartInd=verify(ChemInpLine(StopInd:),' ')+StopInd-1
            StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
            read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)') YY_k
            if (lroot) print*, ' volume fraction, %,    ', YY_k, &
                species_constants(ind_chem,imass)
!
            if (species_constants(ind_chem,imass)>0.) then
             air_mass=air_mass+YY_k*0.01/species_constants(ind_chem,imass)
            endif
!
            if (StartInd==80) exit
!
            stor1(k)=ind_chem
            stor2(k)=YY_k
            k=k+1
          endif
!
        endif
      enddo dataloop
!
!  Stop if air.dat is empty
!
      if (emptyFile)  call fatal_error("air_field", "I can only set existing fields")
      air_mass=1./air_mass
!
      do j=1,k-1
        f(:,:,:,ichemspec(stor1(j)))=stor2(j)*0.01
      enddo
!
      sum_Y=0.
!
      do j=1,nchemspec
        sum_Y=sum_Y+f(:,:,:,ichemspec(j))
      enddo
      do j=1,nchemspec
        f(:,:,:,ichemspec(j))=f(:,:,:,ichemspec(j))/sum_Y
      enddo
!
!
      do j=1,nchemspec
       init_Yk_1(j)=f(l1,m1,n1,ichemspec(j))
       init_Yk_2(j)=f(l1,m1,n1,ichemspec(j))
      enddo
!
      if (mvar < 5) then
        call fatal_error("air_field", "I can only set existing fields")
      endif
        if (ltemperature_nolog) then
          f(:,:,:,iTT)=TT
        else
          f(:,:,:,ilnTT)=alog(TT)!+f(:,:,:,ilnTT)
        endif
        if (ldensity_nolog) then
          f(:,:,:,ilnrho)=(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/TT)/unit_mass*unit_length**3
        else
          tmp=(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/TT)/unit_mass*unit_length**3
          f(:,:,:,ilnrho)=alog(tmp)
        endif
!
        if (ltemperature_nolog) then
          f(:,:,:,iTT)=TT
        else
          f(:,:,:,ilnTT)=alog(TT)!+f(:,:,:,ilnTT)
        endif
        if (ldensity_nolog) then
          f(:,:,:,ilnrho)=(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/TT)/unit_mass*unit_length**3
        else
          tmp=(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/TT)/unit_mass*unit_length**3
          f(:,:,:,ilnrho)=alog(tmp)
        endif
!
      if (lroot) print*, 'local:Air temperature, K', TT
      if (lroot) print*, 'local:Air pressure, dyn', PP
      if (lroot) print*, 'local:Air density, g/cm^3:'
      if (lroot) print '(E10.3)',  PP/(k_B_cgs/m_u_cgs)*air_mass/TT
      if (lroot) print*, 'local:Air mean weight, g/mol', air_mass
      if (lroot) print*, 'local:R', k_B_cgs/m_u_cgs
!
      close(file_id)
!!
    endsubroutine air_field_local
!***********************************************************************
    subroutine ACTOS_data(f)
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: sum_Y, tmp, air_mass
      real, dimension (20000) ::  PP_data, rhow_data, TT_data
      real, dimension (20000) ::  ux_data, uy_data, uz_data, ttime2
      real, dimension (1340) ::  ttime
      real, dimension (6) ::  coeff
!
      logical :: emptyfile=.true., lfind
      logical :: found_specie
      integer :: file_id=123, ind_glob, ind_chem,jj
      character (len=800) :: ChemInpLine
      integer :: i,j,k=1,index_YY, j1,j2,j3, iter, ll1, mm1, nn1
      real ::  TT=300., tmp2
!      real, intent(out) :: PP ! (in dynes = 1atm)
      real, dimension(nchemspec)    :: stor2
      integer, dimension(nchemspec) :: stor1
      real, dimension(29)    :: input_data, input_data2
!
      integer :: StartInd,StopInd,StartInd_1,StopInd_1
      integer :: iostat, i1,i2,i3




      if (lACTOS_write) then

       if (lsinhron) then
         open(143,file="coeff_part.dat")
         do i=1,1300
            read(143,'(f15.6,f15.6,f15.6,f15.6,f15.6,f15.6,f15.6)'),coeff,ttime(i)
         enddo
         close(143)


       endif

!
!      air_mass=0.
      StartInd_1=1; StopInd_1 =0 
     open(file_id,file="ACTOS_data.out")
      open(143,file="ACTOS_new.out")
!      open(file_id,file="ACTOS_xyz_data.out")
!      open(143,file="ACTOS_xyz_new.out")
!
!      dataloop: do
       j=1
       i=1
       do  while (j<780000) 
!
        read(file_id,'(80A)',IOSTAT=iostat) ChemInpLine
         StartInd=1; StopInd =0
         StopInd=index(ChemInpLine(StartInd:),'	')+StartInd-1
!
        if (i_point==1) then
          i=1
        endif
        
        if (i==1) then
        k=1
        do  while (k<30) 
!
         if (StopInd==StartInd) then
            StartInd=StartInd+1
          else
           if (k<29) then
            StopInd=index(ChemInpLine(StartInd:),'	')+StartInd-1
           else
            StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
           endif
            read (unit=ChemInpLine(StartInd:StopInd-1),fmt='(f12.6)') tmp2
             input_data(k)=tmp2

            StartInd=StopInd
            k=k+1
          endif
        enddo
!
         if (lsinhron) then
!
           ttime2(j)=input_data(23)
           lfind=.false.
!
           do jj=1,1300 
!print*, jj, input_data(23), ttime(jj), lfind,  abs(ttime(jj)-input_data(23))
            if (abs(ttime(jj)-input_data(23))<.05) then 
              lfind=.true.
            else
              lfind=.false.
            endif
!
            if ((ttime2(j) == ttime2(j-1)) .or. (j==1)) lfind=.false.
            if (lfind) then 
              write(143,'(29f15.6)'),input_data
             print*, jj, input_data(23), ttime(jj), lfind,  abs(ttime(jj)-input_data(23))
            endif
           enddo
!
         else
!
!          if ((input_data(1)>3545.53) .and. (input_data(1)<3660.53)) then
          if ((input_data(1)>3130.) .and. (input_data(1)<6347.)) then
           write(143,'(29f15.6)'),input_data
          endif
         endif

           i=i+1
        elseif (i>1) then
          i=i+1
          if (i==i_point) i=1
        endif
        j=j+1
        enddo
!
      close(file_id)
      close(143)
!
     endif

!
!        emptyFile=.false.
!      enddo dataloop
!
    if (lACTOS_read) then
        open(143,file="ACTOS_new.out")
        do i=1,Ndata 
          read(143,'(29f15.6)'),input_data
          TT_data(i)=input_data(10)+272.15
!
          PP_data(i)=input_data(7)*1e3   !dyn
!
          rhow_data(i)=input_data(16)*1e-6 !g/cm3
!          uvert(i)=input_data(26)*1e2

          ux_data(i)=input_data(2)*1e2 
          uy_data(i)=input_data(3)*1e2 
          uz_data(i)=input_data(4)*1e2
          Ntot_data(i)=input_data(11)

        enddo
      close(143)
!
!      open(143,file="ACTOS_xyz_new.out")
!        do i=1,Ndata
!          read(143,'(29f15.6)'),input_data2
!          ux_data(i)=uvert(i)/cos(input_data2(19))/cos(input_data2(20))
!          uy_data(i)=uvert(i)/cos(input_data2(19))/sin(input_data2(20))
!          uz_data(i)=uvert(i)/sin(input_data2(19))
          
         
!
!print*,ux_data(i),uy_data(i),i, sin(input_data2(20)),input_data2(20),input_data2(20)*180./3.1415
!!     
!        enddo
!      close(143)
!
        ll1=int((x(l1)-xyz0(1))/dx)
        
        do i=l1,l2
          f(i,:,:,ilnTT)=alog(TT_data(ll1+i-3))
          f(i,:,:,ichemspec(index_H2O))=rhow_data(ll1+i-3)/1e-2  !g/cm3
!
          if (init_ux /= impossible) then
            f(i,:,:,iux)=init_ux
          else
            f(i,:,:,iux)=ux_data(ll1+i-3)
          endif
!          
        enddo

        if (nygrid>1) then
          mm1=int((y(m1)-xyz0(2))/dy)
!          
          if (init_uy /= impossible) then
            f(:,:,:,iuy)=init_uy
          else
           do i=m1,m2
            f(:,i,:,iuy)=uy_data(mm1+i-3)
           enddo
          endif
!            
        endif

        if (nzgrid>1) then
          nn1=int((z(m1)-xyz0(3))/dz)
!
          if (init_uy /= impossible) then
           f(:,:,:,iuz)=init_uz
          else
           do i=n1,n2
            f(:,:,i,iuz)=uz_data(nn1+i-3)
           enddo
          endif
!
        endif

!
       f(:,:,:,ichemspec(index_N2))=0.7
       f(:,:,:,ichemspec(1))=1.-f(:,:,:,ichemspec(index_N2))-f(:,:,:,ichemspec(index_H2O))

!  Stop if air.dat is empty
!
!      if (emptyFile)  call fatal_error("ACTOS data", "I can only set existing fields")
!      air_mass=1./air_mass
!
!         
       do iter=1,4
!   
       sum_Y=0.
       do k=1,nchemspec
         sum_Y=sum_Y + f(:,:,:,ichemspec(k))/species_constants(k,imass)
       enddo
       air_mass=1./sum_Y
!
         do i=l1,l2
           f(i,:,:,ilnrho)=alog(PP_data(ll1+i-3)/(k_B_cgs/m_u_cgs)*&
           air_mass(i,:,:)/exp(f(i,:,:,ilnTT)))/unit_mass*unit_length**3

         enddo
!
       if (iter<4) then
         do i=l1,l2
           f(i,:,:,ichemspec(index_H2O))=rhow_data(ll1+i-3)/exp(f(i,:,:,ilnrho))
         enddo
           f(:,:,:,ichemspec(1))=1.-f(:,:,:,ichemspec(index_N2))-f(:,:,:,ichemspec(index_H2O))
       endif
!
       enddo
!
      if (lroot) print*, 'local:Air temperature, K', maxval(exp(f(l1:l2,m1:m2,n1:n2,ilnTT))), &
                                                     minval(exp(f(l1:l2,m1:m2,n1:n2,ilnTT)))
      if (lroot) print*, 'local:Air pressure, dyn', maxval(PP_data), minval(PP_data)
      if (lroot) print*, 'local:Air density, g/cm^3:'
      if (lroot) print '(E10.3)',  maxval(exp(f(:,:,:,ilnrho))), minval(exp(f(:,:,:,ilnrho)))
      if (lroot) print*, 'local:Air mean weight, g/mol', maxval(air_mass),minval(air_mass)
      if (lroot) print*, 'local:R', k_B_cgs/m_u_cgs
!
       endif
!
    endsubroutine ACTOS_data

!***********************************************************************
    subroutine reinitialization(f, air_mass, PP)
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: sum_Y, air_mass_ar, tmp
      real, dimension (mx,my,mz) :: init_water1_,init_water2_
      real, dimension (my,mz) :: init_water1_min,init_water2_max
      real , dimension (my) :: init_x1_ary, init_x2_ary, del_ary, del_ar1y, del_ar2y
      real , dimension (mz) :: init_x1_arz, init_x2_arz, del_arz, del_ar1z, del_ar2z
!
      integer :: i,j,k, j1,j2,j3, iter
      real :: YY_k, air_mass,  PP, del, psat1, psat2, psf_1, psf_2 
      real :: air_mass_1, air_mass_2, sum1, sum2, init_water_1, init_water_2 
      logical :: spot_exist=.true., lmake_spot, lline_profile=.false.
      real ::  Rgas_loc=8.314472688702992E+7, T_tmp
      real :: aa0= 6.107799961, aa1= 4.436518521e-1
      real :: aa2= 1.428945805e-2, aa3= 2.650648471e-4
      real :: aa4= 3.031240396e-6, aa5= 2.034080948e-8, aa6= 6.136820929e-11

      intent(in) :: air_mass
!

!  Reinitialization of T, water => rho
!
      if (linit_temperature) then
        del=(init_x2-init_x1)*0.2
        if (lcurved_xz) then
          do j=n1,n2         
            init_x1_arz(j)=init_x1*(1-0.6*sin(4.*PI*z(j)/Lxyz(3)))
            init_x2_arz(j)=init_x2*(1+0.6*sin(4.*PI*z(j)/Lxyz(3)))
          enddo
          del_ar1z(:)=del*(1-0.6*sin(4.*PI*z(:)/Lxyz(3)))
          del_ar2z(:)=del*(1+0.6*sin(4.*PI*z(:)/Lxyz(3)))
        elseif (lcurved_xy) then
          do j=m1,m2         
            init_x1_ary(j)=init_x1*(1-0.1*sin(4.*PI*y(j)/Lxyz(2)))
            init_x2_ary(j)=init_x2*(1+0.1*sin(4.*PI*y(j)/Lxyz(2)))
          enddo
          del_ar1y(:)=del*(1-0.1*sin(4.*PI*y(:)/Lxyz(2)))
          del_ar2y(:)=del*(1+0.1*sin(4.*PI*y(:)/Lxyz(2)))
        else
          init_x1_ary=init_x1
          init_x2_ary=init_x2
          init_x1_arz=init_x1
          init_x2_arz=init_x2
          del_ar1y(:)=del
          del_ar2y(:)=del
          del_ar1z(:)=del
          del_ar2z(:)=del
        endif
!
          
        do i=l1,l2
          if (x(i)<0) then
            del_ary=del_ar1y
            del_arz=del_ar1z
          else
            del_ary=del_ar2y 
            del_arz=del_ar2z
          endif
!        
          if (ltanh_prof_xy) then
            do j=m1,m2
              f(i,j,:,ilnTT)=log((init_TT2+init_TT1)*0.5  &
                             +((init_TT2-init_TT1)*0.5)  &
                  *(exp(x(i)/del_ary(j))-exp(-x(i)/del_ary(j))) &
                  /(exp(x(i)/del_ary(j))+exp(-x(i)/del_ary(j))))
            enddo
          elseif (ltanh_prof_xz) then
            do j=n1,n2
            f(i,:,j,ilnTT)=log((init_TT2+init_TT1)*0.5  &
                             +((init_TT2-init_TT1)*0.5)  &
              *(1.-exp(-2.*x(i)/del_arz(j))) &
              /(1.+exp(-2.*x(i)/del_arz(j))))
            enddo
          else
          do j=m1,m2
            if (x(i)<=init_x1_ary(j)) then
              f(i,j,:,ilnTT)=alog(init_TT1)
            endif
            if (x(i)>=init_x2_ary(j)) then
              f(i,j,:,ilnTT)=alog(init_TT2)
            endif
            if (x(i)>init_x1_ary(j) .and. x(i)<init_x2_ary(j)) then
              if (init_x1_ary(j) /= init_x2_ary(j)) then
                f(i,j,:,ilnTT)=&
                   alog((x(i)-init_x1_ary(j))/(init_x2_ary(j)-init_x1_ary(j)) &
                   *(init_TT2-init_TT1)+init_TT1)
              endif
            endif
          enddo
          endif
        enddo
!        
        if (ldensity_nolog) then
!          f(:,:,:,ilnrho)=(PP/(k_B_cgs/m_u_cgs)*&
!            air_mass/exp(f(:,:,:,ilnTT)))/unit_mass*unit_length**3
        else
          tmp=(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/exp(f(:,:,:,ilnTT)))/unit_mass*unit_length**3
          f(:,:,:,ilnrho)=alog(tmp)
        endif
      endif
!
       if (lreinit_water) then
!
       if (init_TT2==0) init_TT2=init_TT1
!
        T_tmp=init_TT1-273.15
        psat1=(aa0 + aa1*T_tmp    + aa2*T_tmp**2  &
                   + aa3*T_tmp**3 + aa4*T_tmp**4  &
                   + aa5*T_tmp**5 + aa6*T_tmp**6)*1e3
        T_tmp=init_TT2-273.15
        psat2=(aa0 + aa1*T_tmp    + aa2*T_tmp**2  &
                   + aa3*T_tmp**3 + aa4*T_tmp**4  &
                   + aa5*T_tmp**5 + aa6*T_tmp**6)*1e3
!
        psf_1=psat1 
        if (r02 /= 0) then
          psf_2=psat2
        else
          psf_2=psat2
        endif
!
! Recalculation of the air_mass for different boundary conditions
!
        air_mass_1=air_mass
        air_mass_2=air_mass

        do iter=1,3
!
          init_Yk_1(index_H2O)=psf_1/(PP*air_mass_1/18.)*dYw1
          init_Yk_2(index_H2O)=psf_2/(PP*air_mass_2/18.)*dYw2
!
           sum1=0.
           sum2=0.
           do k=1,nchemspec
            if (ichemspec(k)/=ichemspec(index_N2)) then
              sum1=sum1+init_Yk_1(k)
              sum2=sum2+init_Yk_2(k)
            endif
           enddo
!
           init_Yk_1(index_N2)=1.-sum1
           init_Yk_2(index_N2)=1.-sum2
!
!  Recalculation of air_mass 
!
             sum1=0.
             sum2=0.
             do k=1,nchemspec
               sum1=sum1 + init_Yk_1(k)/species_constants(k,imass)
               sum2=sum2 + init_Yk_2(k)/species_constants(k,imass)
             enddo
               air_mass_1=1./sum1
               air_mass_2=1./sum2
        enddo

           init_water_1=init_Yk_1(index_H2O)
           init_water_2=init_Yk_2(index_H2O)
!
! End of Recalculation of the air_mass for different boundary conditions
!
!  Different profiles
!
           if (ltanh_prof_xz .or. ltanh_prof_xy) then
             do i=l1,l2
               f(i,:,:,ichemspec(index_H2O))= &
                   (init_water_2+init_water_1)*0.5  &
                 +((init_water_2-init_water_1)*0.5)  &
                   *(exp(x(i)/del)-exp(-x(i)/del)) &
                   /(exp(x(i)/del)+exp(-x(i)/del))
             enddo
           elseif (lwet_spots) then
              f(:,:,:,ilnTT)=log(init_TT1)
              f(:,:,:,ichemspec(index_H2O))=init_water_1
              call spot_init(f,init_TT2,init_water_2)
           elseif (.not. lwet_spots) then
! Initial conditions for the  0dcase: cond_evap
             f(:,:,:,ichemspec(index_H2O))=init_Yk_1(index_H2O)
           endif
!
!
           sum_Y=0.
           do k=1,nchemspec
             if (ichemspec(k)/=ichemspec(index_N2)) &
               sum_Y=sum_Y+f(:,:,:,ichemspec(k))
           enddo
           f(:,:,:,ichemspec(index_N2))=1.-sum_Y
!
           sum_Y=0.
           do k=1,nchemspec
             sum_Y=sum_Y+f(:,:,:,ichemspec(k)) &
               /species_constants(k,imass)
           enddo
           air_mass_ar=1./sum_Y
!
! end of loot do iter=1,2
!         enddo
!
         if (ldensity_nolog) then
           f(:,:,:,ilnrho)=(PP/(k_B_cgs/m_u_cgs)&
            *air_mass_ar/exp(f(:,:,:,ilnTT)))/unit_mass*unit_length**3
         else
           tmp=(PP/(k_B_cgs/m_u_cgs) &
            *air_mass_ar/exp(f(:,:,:,ilnTT)))/unit_mass*unit_length**3
           f(:,:,:,ilnrho)=alog(tmp) 
         endif
!
         if ((nxgrid>1) .and. (nygrid==1).and. (nzgrid==1)) then
            f(:,:,:,iux)=f(:,:,:,iux)+init_ux
         endif
!       
         if (lroot) print*, ' Saturation Pressure, Pa   ', psf_1, psf_2
         if (lroot) print*, ' psf, Pa   ',  psf_1, psf_2
         if (lroot) print*, ' pw1, Pa   ', (exp(f(l1,m1,n1,ilnrho))*Rgas_loc*init_TT1/18.), PP*air_mass_1/18.
         if (lroot) print*, ' pw2, Pa   ', (exp(f(l2,m1,n1,ilnrho))*Rgas_loc*init_TT2/18.), PP*air_mass_2/18.
         if (lroot) print*, ' saturated water mass fraction', psat1/PP, psat2/PP
!         if (lroot) print*, 'New Air density, g/cm^3:'
!         if (lroot) print '(E10.3)',  PP/(k_B_cgs/m_u_cgs)*maxval(air_mass_ar)/TT
         if (lroot) print*, 'New Air mean weight, g/mol', maxval(air_mass_ar)
         if  ((lroot) .and. (nx >1 )) then
            print*, 'density', exp(f(l1,4,4,ilnrho)), exp(f(l2,4,4,ilnrho))
          endif
         if  ((lroot) .and. (nx >1 )) then
            print*, 'temperature', exp(f(l1,4,4,ilnTT)), exp(f(l2,4,4,ilnTT))
          endif
       endif
!
    endsubroutine reinitialization
!*********************************************************************************************************
    subroutine spot_init(f,init_TT_2,init_H2O_2)
!
!  Initialization of the dust spot positions and dust distribution
!
!  10-may-10/Natalia: coded
!
      use General, only: random_number_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: k, j, j1,j2,j3, lx=0,ly=0,lz=0
      real ::  RR, init_TT_2, init_H2O_2
      real, dimension (3,spot_number) :: spot_posit
      logical :: spot_exist=.true.
! 
      spot_posit(:,:)=0.0
      do j=1,spot_number
        spot_exist=.true.
        lx=0;ly=0; lz=0
        if (nxgrid/=1) then
          lx=1
            call random_number_wrapper(spot_posit(1,j))
            spot_posit(1,j)=spot_posit(1,j)*Lxyz(1)
          if ((spot_posit(1,j)-1.5*spot_size<xyz0(1)) .or. &
            (spot_posit(1,j)+1.5*spot_size>xyz0(1)+Lxyz(1)))  &
            spot_exist=.false.
            print*,'positx',spot_posit(1,j),spot_exist
!          if ((spot_posit(1,j)-1.5*spot_size<xyz0(1)) )  &
!            spot_exist=.false.
!            print*,'positx',spot_posit(1,j),spot_exist
        endif
        if (nygrid/=1) then
          ly=1
            call random_number_wrapper(spot_posit(2,j))
            spot_posit(2,j)=spot_posit(2,j)*Lxyz(2)
          if ((spot_posit(2,j)-1.5*spot_size<xyz0(2)) .or. &
           (spot_posit(2,j)+1.5*spot_size>xyz0(2)+Lxyz(2)))  &
          spot_exist=.false.
            print*,'posity',spot_posit(2,j),spot_exist
        endif
        if (nzgrid/=1) then
          lz=1
            call random_number_wrapper(spot_posit(3,j))
            spot_posit(3,j)=spot_posit(3,j)*Lxyz(3)
          if ((spot_posit(3,j)-1.5*spot_size<xyz0(3)) .or. &
           (spot_posit(3,j)+1.5*spot_size>xyz0(3)+Lxyz(3)))  &
           spot_exist=.false.
        endif
             do j1=1,mx; do j2=1,my; do j3=1,mz
               RR= (lx*x(j1)-spot_posit(1,j))**2 &
                   +ly*(y(j2)-spot_posit(2,j))**2 &
                   +lz*(z(j3)-spot_posit(3,j))**2
               RR=sqrt(RR)
!
               if ((RR<spot_size) .and. (spot_exist)) then
                f(j1,j2,j3,ichemspec(index_H2O)) = init_H2O_2
                f(j1,j2,j3,ilnTT)=log(init_TT_2)
               endif
             enddo; enddo; enddo

      enddo
!
    endsubroutine spot_init
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
