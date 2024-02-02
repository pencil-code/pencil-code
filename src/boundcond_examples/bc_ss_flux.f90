!
!  Example of preprocessor output after unrolling and inlining functions calls
!  Has declarations for GPU parallelization, but not for starting ending target registration
!
!***********************************************************************
subroutine bc_ss_flux(f,topbot,lone_sided)
use General
use DensityMethods
use Deriv
integer :: i
real :: facs
integer :: l1_,l2_,sdf
real, dimension (:,:,:,:) :: f
integer :: explicit_index_1, explicit_index_0, explicit_index_2, explicit_index_3
real, dimension (size(f,1),size(f,2)) :: tmp_xy,cs2_xy,rho_xy
real, pointer :: hcond0_kramers, nkramers
real, dimension (size(rho_xy(:,m))) :: fac
integer :: k,off
integer :: n_save, m_save
integer, intent(IN) :: topbot
real, pointer :: Fbot,Ftop,FtopKtop,FbotKbot,chi
logical, pointer :: lmultilayer, lheatc_chiconst, lheatc_kramers
logical, optional :: lone_sided
if (ldebug) print*,'bc_ss_flux: ENTER - cs20,cs0=',cs20,cs0
call get_shared_variable('Fbot',Fbot,caller='bc_ss_flux')
call get_shared_variable('Ftop',Ftop)
call get_shared_variable('FbotKbot',FbotKbot)
call get_shared_variable('FtopKtop',FtopKtop)
call get_shared_variable('chi',chi)
call get_shared_variable('lmultilayer',lmultilayer)
call get_shared_variable('lheatc_chiconst',lheatc_chiconst)
call get_shared_variable('lheatc_kramers',lheatc_kramers)
if (lheatc_kramers) then
call get_shared_variable('hcond0_kramers',hcond0_kramers)
call get_shared_variable('nkramers',nkramers)
endif
select case (topbot)
case(BOT)
if (pretend_lnTT) then
!$omp teams distribute parallel do
do explicit_index_0=1,size(tmp_xy,dim=1)
do explicit_index_1=1,size(tmp_xy,dim=2)
tmp_xy(explicit_index_0,explicit_index_1)=-FbotKbot/exp(f(explicit_index_0,explicit_index_1,n1,iss))
end do
end do
!$omp end teams distribute parallel do

do i=1,nghost
!$omp teams distribute parallel do
do explicit_index_0=1,size(f,dim=1)
do explicit_index_1=1,size(f,dim=2)
f(explicit_index_0,explicit_index_1,n1-i,iss)=f(explicit_index_0,explicit_index_1,n1+i,iss)-dz2_bound(-i)*tmp_xy(explicit_index_0,explicit_index_1)
end do
end do
!$omp end teams distribute parallel do

enddo
else
if (ldensity_nolog) then
if (lreference_state) then
!$omp teams distribute parallel do
do explicit_index_0=1,l2-l1+1
do explicit_index_1=1,size(rho_xy,dim=2)
rho_xy(l1 + explicit_index_0,explicit_index_1)= f(l1 + explicit_index_0,explicit_index_1,n1,ilnrho)  +reference_state(explicit_index_0,iref_rho)
end do
end do
!$omp end teams distribute parallel do

else
!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
do explicit_index_1=1,size(rho_xy,dim=2)
rho_xy(explicit_index_0,explicit_index_1)=f(explicit_index_0,explicit_index_1,n1,ilnrho)
end do
end do
!$omp end teams distribute parallel do

endif
else
!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
do explicit_index_1=1,size(rho_xy,dim=2)
rho_xy(explicit_index_0,explicit_index_1)=exp(f(explicit_index_0,explicit_index_1,n1,ilnrho))
end do
end do
!$omp end teams distribute parallel do

endif
!$omp teams distribute parallel do
do explicit_index_0=1,size(cs2_xy,dim=1)
do explicit_index_1=1,size(cs2_xy,dim=2)
cs2_xy(explicit_index_0,explicit_index_1) = f(explicit_index_0,explicit_index_1,n1,iss)
end do
end do
!$omp end teams distribute parallel do

!$omp teams distribute parallel do
do explicit_index_0=1,l2-l1+1
do explicit_index_1=1,size(cs2_xy,dim=2)
if (lreference_state)  cs2_xy(l1 + explicit_index_0,explicit_index_1) = cs2_xy(l1 + explicit_index_0,explicit_index_1) + reference_state(explicit_index_0,iref_s)
end do
end do
!$omp end teams distribute parallel do

if (ldensity_nolog) then
!$omp teams distribute parallel do
do explicit_index_0=1,size(cs2_xy,dim=1)
do explicit_index_1=1,size(cs2_xy,dim=2)
cs2_xy(explicit_index_0,explicit_index_1)=cs20*exp(gamma_m1*(log(rho_xy(explicit_index_0,explicit_index_1))-lnrho0)+cv1*cs2_xy(explicit_index_0,explicit_index_1))
end do
end do
!$omp end teams distribute parallel do

else
!$omp teams distribute parallel do
do explicit_index_0=1,size(cs2_xy,dim=1)
do explicit_index_1=1,size(cs2_xy,dim=2)
cs2_xy(explicit_index_0,explicit_index_1)=cs20*exp(gamma_m1*(f(explicit_index_0,explicit_index_1,n1,ilnrho)-lnrho0)+cv1*cs2_xy(explicit_index_0,explicit_index_1))
end do
end do
!$omp end teams distribute parallel do

endif
if (lheatc_chiconst) then
!$omp teams distribute parallel do
do explicit_index_0=1,size(tmp_xy,dim=1)
do explicit_index_1=1,size(tmp_xy,dim=2)
tmp_xy(explicit_index_0,explicit_index_1)=Fbot/(rho_xy(explicit_index_0,explicit_index_1)*chi*cs2_xy(explicit_index_0,explicit_index_1))
end do
end do
!$omp end teams distribute parallel do

else if (lheatc_kramers) then
!$omp teams distribute parallel do
do explicit_index_0=1,size(tmp_xy,dim=1)
do explicit_index_1=1,size(tmp_xy,dim=2)
tmp_xy(explicit_index_0,explicit_index_1)=Fbot*rho_xy(explicit_index_0,explicit_index_1)**(2*nkramers)*(cp*gamma_m1)**(6.5*nkramers)  /(hcond0_kramers*cs2_xy(explicit_index_0,explicit_index_1)**(6.5*nkramers+1.))
end do
end do
!$omp end teams distribute parallel do

else
!$omp teams distribute parallel do
do explicit_index_0=1,size(tmp_xy,dim=1)
do explicit_index_1=1,size(tmp_xy,dim=2)
tmp_xy(explicit_index_0,explicit_index_1)=FbotKbot/cs2_xy(explicit_index_0,explicit_index_1)
end do
end do
!$omp end teams distribute parallel do

endif
if (loptest(lone_sided)) then
call not_implemented('bc_ss_flux', 'one-sided BC')
n_save=n
m_save=m
n=iz
do m=1,size(f,2)
if (3==1) then
if (nxgrid/=1) then
!$omp teams distribute parallel do
do explicit_index_0=1,size(fac,dim=1)
fac(explicit_index_0)=(1./60)*dx_1(l1 + explicit_index_0)
end do
!$omp end teams distribute parallel do

!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
rho_xy(explicit_index_0,m)=fac(explicit_index_0)*(+ 45.0*(f(l1+1 + explicit_index_0,m,n,ilnrho)-f(l1-1 + explicit_index_0,m,n,ilnrho))  -      9.0*(f(l1+2 + explicit_index_0,m,n,ilnrho)-f(l1-2 + explicit_index_0,m,n,ilnrho))  +      (f(l1+3 + explicit_index_0,m,n,ilnrho)-f(l1-3 + explicit_index_0,m,n,ilnrho)))
end do
!$omp end teams distribute parallel do

else
!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
rho_xy(explicit_index_0,m)=0.
end do
!$omp end teams distribute parallel do

if (ip<=5) print*, 'der_other: Degenerate case in x-direction'
endif
else
sdf=size(rho_xy(:,m))
if (sdf>nx) then
l1_=1
l2_=sdf
else
l1_=l1
l2_=l2
endif
if (3==2) then
if (nygrid/=1) then
facs=(1./60.)*dy_1(m)
!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
rho_xy(explicit_index_0,m)=facs*(+ 45.0*(f(l1_ + explicit_index_0,m+1,n,ilnrho)-f(l1_ + explicit_index_0,m-1,n,ilnrho))  -     9.0*(f(l1_ + explicit_index_0,m+2,n,ilnrho)-f(l1_ + explicit_index_0,m-2,n,ilnrho))  +      (f(l1_ + explicit_index_0,m+3,n,ilnrho)-f(l1_ + explicit_index_0,m-3,n,ilnrho)))
end do
!$omp end teams distribute parallel do

!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
if (lspherical_coords)   rho_xy(explicit_index_0,m)=rho_xy(explicit_index_0,m)*r1_mn(explicit_index_0)
end do
!$omp end teams distribute parallel do

!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
if (lcylindrical_coords) rho_xy(explicit_index_0,m)=rho_xy(explicit_index_0,m)*rcyl_mn1(explicit_index_0)
end do
!$omp end teams distribute parallel do

else
!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
rho_xy(explicit_index_0,m)=0.
end do
!$omp end teams distribute parallel do

if (ip<=5) print*, 'der_other: Degenerate case in y-direction'
endif
elseif (3==3) then
if (nzgrid/=1) then
if (lcoarse_mn) then
facs = (1./60.) * dz_1(n) * nphis1(m)
!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
rho_xy(explicit_index_0,m)=facs*(+ 45.0*(f(l1_ + explicit_index_0,m,+1,ilnrho)-f(l1_ + explicit_index_0,m,-1,ilnrho))  -  9.0*(f(l1_ + explicit_index_0,m,+2,ilnrho)-f(l1_ + explicit_index_0,m,-2,ilnrho))       +      (f(l1_ + explicit_index_0,m,+3,ilnrho)-f(l1_ + explicit_index_0,m,-3,ilnrho)) )
end do
!$omp end teams distribute parallel do

else
facs = (1./60.) * dz_1(n)
!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
rho_xy(explicit_index_0,m)=facs*(+ 45.0*(f(l1_ + explicit_index_0,m,n+1,ilnrho)-f(l1_ + explicit_index_0,m,n-1,ilnrho))  -     9.0*(f(l1_ + explicit_index_0,m,n+2,ilnrho)-f(l1_ + explicit_index_0,m,n-2,ilnrho))  +      (f(l1_ + explicit_index_0,m,n+3,ilnrho)-f(l1_ + explicit_index_0,m,n-3,ilnrho)))
end do
!$omp end teams distribute parallel do

endif
!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
if (lspherical_coords) rho_xy(explicit_index_0,m)=rho_xy(explicit_index_0,m)*r1_mn(explicit_index_0)*sin1th(m)
end do
!$omp end teams distribute parallel do

else
!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
rho_xy(explicit_index_0,m)=0.
end do
!$omp end teams distribute parallel do

if (ip<=5) print*, 'der_other: Degenerate case in z-direction'
endif
endif
endif
enddo
n=n_save
m=m_save
if (ldensity_nolog) then
if (lreference_state) then
!$omp teams distribute parallel do
do explicit_index_0=1,l2-l1+1
do explicit_index_1=1,size(rho_xy,dim=2)
rho_xy(l1 + explicit_index_0,explicit_index_1) = &
rho_xy(l1 + explicit_index_0,explicit_index_1)/(f(l1 + explicit_index_0,explicit_index_1,n1,ilnrho) &
 +reference_state(explicit_index_0,iref_rho))
end do
end do
!$omp end teams distribute parallel do

else
!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
do explicit_index_1=1,size(rho_xy,dim=2)
rho_xy(explicit_index_0,explicit_index_1) = rho_xy(explicit_index_0,explicit_index_1)/f(explicit_index_0,explicit_index_1,n1,ilnrho)
end do
end do
!$omp end teams distribute parallel do

endif
endif
call bval_from_neumann(f,topbot,iss,3,rho_xy)
if (loptest(.true.)) then
off=nghost-1
else
off=nghost
endif
if (topbot==BOT) then
if (3==1) then
do k=l1-1,l1-off,-1
!$omp teams distribute parallel do
do explicit_index_0=1,size(f,dim=2)
do explicit_index_1=1,size(f,dim=3)
f(k,explicit_index_0,explicit_index_1,iss)=   7.*(f(k+1,explicit_index_0,explicit_index_1,iss)-f(k+6,explicit_index_0,explicit_index_1,iss))  -21.*(f(k+2,explicit_index_0,explicit_index_1,iss)-f(k+5,explicit_index_0,explicit_index_1,iss))      +35.*(f(k+3,explicit_index_0,explicit_index_1,iss)-f(k+4,explicit_index_0,explicit_index_1,iss))  + f(k+7,explicit_index_0,explicit_index_1,iss)
end do
end do
!$omp end teams distribute parallel do

enddo
elseif (3==2) then
do k=m1-1,m1-off,-1
!$omp teams distribute parallel do
do explicit_index_0=1,size(f,dim=1)
do explicit_index_1=1,size(f,dim=3)
f(explicit_index_0,k,explicit_index_1,iss)=   7.*(f(explicit_index_0,k+1,explicit_index_1,iss)-f(explicit_index_0,k+6,explicit_index_1,iss))  -21.*(f(explicit_index_0,k+2,explicit_index_1,iss)-f(explicit_index_0,k+5,explicit_index_1,iss))      +35.*(f(explicit_index_0,k+3,explicit_index_1,iss)-f(explicit_index_0,k+4,explicit_index_1,iss))  + f(explicit_index_0,k+7,explicit_index_1,iss)
end do
end do
!$omp end teams distribute parallel do

enddo
elseif (3==3) then
do k=n1-1,n1-off,-1
!$omp teams distribute parallel do
do explicit_index_0=1,size(f,dim=1)
do explicit_index_1=1,size(f,dim=2)
f(explicit_index_0,explicit_index_1,k,iss)=   7.*(f(explicit_index_0,explicit_index_1,k+1,iss)-f(explicit_index_0,explicit_index_1,k+6,iss))  -21.*(f(explicit_index_0,explicit_index_1,k+2,iss)-f(explicit_index_0,explicit_index_1,k+5,iss))      +35.*(f(explicit_index_0,explicit_index_1,k+3,iss)-f(explicit_index_0,explicit_index_1,k+4,iss))  + f(explicit_index_0,explicit_index_1,k+7,iss)
end do
end do
!$omp end teams distribute parallel do

enddo
endif
else
if (3==1) then
do k=l2+1,l2+off
!$omp teams distribute parallel do
do explicit_index_0=1,size(f,dim=2)
do explicit_index_1=1,size(f,dim=3)
f(k,explicit_index_0,explicit_index_1,iss)=   7.*(f(k-1,explicit_index_0,explicit_index_1,iss)-f(k-6,explicit_index_0,explicit_index_1,iss))  -21.*(f(k-2,explicit_index_0,explicit_index_1,iss)-f(k-5,explicit_index_0,explicit_index_1,iss))      +35.*(f(k-3,explicit_index_0,explicit_index_1,iss)-f(k-4,explicit_index_0,explicit_index_1,iss))  + f(k-7,explicit_index_0,explicit_index_1,iss)
end do
end do
!$omp end teams distribute parallel do

enddo
elseif (3==2) then
do k=m2+1,m2+off
!$omp teams distribute parallel do
do explicit_index_0=1,size(f,dim=1)
do explicit_index_1=1,size(f,dim=3)
f(explicit_index_0,k,explicit_index_1,iss)=   7.*(f(explicit_index_0,k-1,explicit_index_1,iss)-f(explicit_index_0,k-6,explicit_index_1,iss))  -21.*(f(explicit_index_0,k-2,explicit_index_1,iss)-f(explicit_index_0,k-5,explicit_index_1,iss))      +35.*(f(explicit_index_0,k-3,explicit_index_1,iss)-f(explicit_index_0,k-4,explicit_index_1,iss))  + f(explicit_index_0,k-7,explicit_index_1,iss)
end do
end do
!$omp end teams distribute parallel do

enddo
elseif (3==3) then
do k=n2+1,n2+off
!$omp teams distribute parallel do
do explicit_index_0=1,size(f,dim=1)
do explicit_index_1=1,size(f,dim=2)
f(explicit_index_0,explicit_index_1,k,iss)=   7.*(f(explicit_index_0,explicit_index_1,k-1,iss)-f(explicit_index_0,explicit_index_1,k-6,iss))  -21.*(f(explicit_index_0,explicit_index_1,k-2,iss)-f(explicit_index_0,explicit_index_1,k-5,iss))      +35.*(f(explicit_index_0,explicit_index_1,k-3,iss)-f(explicit_index_0,explicit_index_1,k-4,iss))  + f(explicit_index_0,explicit_index_1,k-7,iss)
end do
end do
!$omp end teams distribute parallel do

enddo
endif
endif
else
do i=1,nghost
!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
do explicit_index_1=1,size(rho_xy,dim=2)
rho_xy(explicit_index_0,explicit_index_1) = f(explicit_index_0,explicit_index_1,n1+i,ilnrho)-f(explicit_index_0,explicit_index_1,n1-i,ilnrho)
end do
end do
!$omp end teams distribute parallel do

if (ldensity_nolog) then
if (lreference_state) then
!$omp teams distribute parallel do
do explicit_index_0=1,l2-l1+1
do explicit_index_1=1,size(rho_xy,dim=2)
rho_xy(l1 + explicit_index_0,explicit_index_1) = rho_xy(l1 + explicit_index_0,explicit_index_1)/(f(l1 + explicit_index_0,explicit_index_1,n1,ilnrho)  +reference_state(explicit_index_0,iref_rho))
end do
end do
!$omp end teams distribute parallel do

else
!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
do explicit_index_1=1,size(rho_xy,dim=2)
rho_xy(explicit_index_0,explicit_index_1) = rho_xy(explicit_index_0,explicit_index_1)/f(explicit_index_0,explicit_index_1,n1,ilnrho)
end do
end do
!$omp end teams distribute parallel do

endif
endif
!$omp teams distribute parallel do
do explicit_index_0=1,size(f,dim=1)
do explicit_index_1=1,size(f,dim=2)
f(explicit_index_0,explicit_index_1,n1-i,iss)=f(explicit_index_0,explicit_index_1,n1+i,iss)+cp*(cp-cv)*(rho_xy(explicit_index_0,explicit_index_1)+dz2_bound(-i)*tmp_xy(explicit_index_0,explicit_index_1))
end do
end do
!$omp end teams distribute parallel do

enddo
endif
endif
case(TOP)
if (pretend_lnTT) then
!$omp teams distribute parallel do
do explicit_index_0=1,size(tmp_xy,dim=1)
do explicit_index_1=1,size(tmp_xy,dim=2)
tmp_xy(explicit_index_0,explicit_index_1)=-FtopKtop/exp(f(explicit_index_0,explicit_index_1,n2,iss))
end do
end do
!$omp end teams distribute parallel do

do i=1,nghost
!$omp teams distribute parallel do
do explicit_index_0=1,size(f,dim=1)
do explicit_index_1=1,size(f,dim=2)
f(explicit_index_0,explicit_index_1,n2-i,iss)=f(explicit_index_0,explicit_index_1,n2+i,iss)-dz2_bound(i)*tmp_xy(explicit_index_0,explicit_index_1)
end do
end do
!$omp end teams distribute parallel do

enddo
else
call getrho(f(:,:,n2,ilnrho),rho_xy)
!$omp teams distribute parallel do
do explicit_index_0=1,size(cs2_xy,dim=1)
do explicit_index_1=1,size(cs2_xy,dim=2)
cs2_xy(explicit_index_0,explicit_index_1) = f(explicit_index_0,explicit_index_1,n2,iss)
end do
end do
!$omp end teams distribute parallel do

!$omp teams distribute parallel do
do explicit_index_0=1,l2-l1+1
do explicit_index_1=1,size(cs2_xy,dim=2)
if (lreference_state)  cs2_xy(l1 + explicit_index_0,explicit_index_1) = cs2_xy(l1 + explicit_index_0,explicit_index_1) + reference_state(explicit_index_0,iref_s)
end do
end do
!$omp end teams distribute parallel do

if (ldensity_nolog) then
!$omp teams distribute parallel do
do explicit_index_0=1,size(cs2_xy,dim=1)
do explicit_index_1=1,size(cs2_xy,dim=2)
cs2_xy(explicit_index_0,explicit_index_1)=cs20*exp(gamma_m1*(log(rho_xy(explicit_index_0,explicit_index_1))-lnrho0)+cv1*cs2_xy(explicit_index_0,explicit_index_1))
end do
end do
!$omp end teams distribute parallel do

else
!$omp teams distribute parallel do
do explicit_index_0=1,size(cs2_xy,dim=1)
do explicit_index_1=1,size(cs2_xy,dim=2)
cs2_xy(explicit_index_0,explicit_index_1)=cs20*exp(gamma_m1*(f(explicit_index_0,explicit_index_1,n2,ilnrho)-lnrho0)+cv1*cs2_xy(explicit_index_0,explicit_index_1))
end do
end do
!$omp end teams distribute parallel do

endif
if (lheatc_chiconst) then
!$omp teams distribute parallel do
do explicit_index_0=1,size(tmp_xy,dim=1)
do explicit_index_1=1,size(tmp_xy,dim=2)
tmp_xy(explicit_index_0,explicit_index_1)=Ftop/(rho_xy(explicit_index_0,explicit_index_1)*chi*cs2_xy(explicit_index_0,explicit_index_1))
end do
end do
!$omp end teams distribute parallel do

else if (lheatc_kramers) then
!$omp teams distribute parallel do
do explicit_index_0=1,size(tmp_xy,dim=1)
do explicit_index_1=1,size(tmp_xy,dim=2)
tmp_xy(explicit_index_0,explicit_index_1)=Ftop*rho_xy(explicit_index_0,explicit_index_1)**(2*nkramers)*(cp*gamma_m1)**(6.5*nkramers)  /(hcond0_kramers*cs2_xy(explicit_index_0,explicit_index_1)**(6.5*nkramers+1.))
end do
end do
!$omp end teams distribute parallel do

else
!$omp teams distribute parallel do
do explicit_index_0=1,size(tmp_xy,dim=1)
do explicit_index_1=1,size(tmp_xy,dim=2)
tmp_xy(explicit_index_0,explicit_index_1)=FtopKtop/cs2_xy(explicit_index_0,explicit_index_1)
end do
end do
!$omp end teams distribute parallel do

endif
if (loptest(lone_sided)) then
call not_implemented('bc_ss_flux', 'one-sided BC')
call getderlnrho_z(f,n2,rho_xy)
call bval_from_neumann(f,topbot,iss,3,rho_xy)
call set_ghosts_for_onesided_ders(f,topbot,iss,3,.true.)
else
do i=1,nghost
!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
do explicit_index_1=1,size(rho_xy,dim=2)
rho_xy(explicit_index_0,explicit_index_1) = f(explicit_index_0,explicit_index_1,n2+i,ilnrho)-f(explicit_index_0,explicit_index_1,n2-i,ilnrho)
end do
end do
!$omp end teams distribute parallel do

if (ldensity_nolog) then
if (lreference_state) then
!$omp teams distribute parallel do
do explicit_index_0=1,l2-l1+1
do explicit_index_1=1,size(rho_xy,dim=2)
rho_xy(l1 + explicit_index_0,explicit_index_1) = rho_xy(l1 + explicit_index_0,explicit_index_1)/(f(l1 + explicit_index_0,explicit_index_1,n2,ilnrho)  +reference_state(explicit_index_0,iref_rho))
end do
end do
!$omp end teams distribute parallel do

else
!$omp teams distribute parallel do
do explicit_index_0=1,size(rho_xy,dim=1)
do explicit_index_1=1,size(rho_xy,dim=2)
rho_xy(explicit_index_0,explicit_index_1) = rho_xy(explicit_index_0,explicit_index_1)/f(explicit_index_0,explicit_index_1,n2,ilnrho)
end do
end do
!$omp end teams distribute parallel do

endif
endif
!$omp teams distribute parallel do
do explicit_index_0=1,size(f,dim=1)
do explicit_index_1=1,size(f,dim=2)
f(explicit_index_0,explicit_index_1,n2+i,iss)=f(explicit_index_0,explicit_index_1,n2-i,iss)+cp*(cp-cv)*(-rho_xy(explicit_index_0,explicit_index_1)-dz2_bound(i)*tmp_xy(explicit_index_0,explicit_index_1))
end do
end do
!$omp end teams distribute parallel do

enddo
endif
endif
case default
call fatal_error('bc_ss_flux','invalid argument')
endselect
endsubroutine bc_ss_flux
!***********************************************************************