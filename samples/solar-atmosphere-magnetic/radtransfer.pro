@data/pc_constants.pro
arad = 1.8049444e-14

dtau_thresh_min = 2.9730177E-02
dtau_thresh_max = 87.33655

Srad = arad*exp(4*lnTT)

nrad = 0
lrad = 1
mrad = 1

dlength = sqrt((dx*lrad)^2+(dy*mrad)^2+(dz*nrad)^2)

tau = fltarr(mx,my,mz)
Qrad = fltarr(mx,my,mz)

for n = n1,n2 do begin
for m = m1,m2 do begin
for l = l1,l2 do begin

  dtau_m = sqrt(kapparho[l-lrad,m-mrad,n-nrad]*kapparho[l,m,n])*dlength
  dtau_p = sqrt(kapparho[l,m,n]*kapparho[l+lrad,m+mrad,n+nrad])*dlength

  dSdtau_m=(Srad[l,m,n]-Srad[l-lrad,m-mrad,n-nrad])/dtau_m
  dSdtau_p=(Srad[l+lrad,m+mrad,n+nrad]-Srad[l,m,n])/dtau_p

  Srad1st=(dSdtau_p*dtau_m+dSdtau_m*dtau_p)/(dtau_m+dtau_p)
  Srad2nd=2*(dSdtau_p-dSdtau_m)/(dtau_m+dtau_p)

  if (dtau_m gt dtau_thresh_max) then begin
    emdtau=0.0
    emdtau1=1.0
    emdtau2=-1.0
  endif else if (dtau_m lt dtau_thresh_min) then begin
    emdtau1=dtau_m*(1-0.5*dtau_m*(1-dtau_m/3))
    emdtau=1-emdtau1
    emdtau2=-dtau_m^2*(0.5-dtau_m/3)
  endif else begin
    emdtau=exp(-dtau_m)
    emdtau1=1-emdtau
    emdtau2=emdtau*(1+dtau_m)-1
  endelse
  tau[l,m,n]=tau[l-lrad,m-mrad,n-nrad]+dtau_m
  Qrad[l,m,n]=Qrad[l-lrad,m-mrad,n-nrad]*emdtau $
             -Srad1st*emdtau1-Srad2nd*emdtau2

endfor
endfor
endfor

if (lrad ne 0 or mrad ne 0) then begin

  Qrad0 = fltarr(mx,my,mz)

  Qrad0[l1-lrad,m1-mrad:m2,n1-nrad:n2] = $
   Qrad[l2     ,m1-mrad:m2,n1-nrad:n2]
  Qrad0[l1-lrad:l2,m1-mrad,n1-nrad:n2] = $
   Qrad[l1-lrad:l2,m2     ,n1-nrad:n2]
  ;Qrad0[l1-lrad,m1-mrad,*] = Qrad[l2,m2,*]

  for n = n1,n2 do begin
  for m = m1,m2 do begin
  for l = l1,l2 do begin

    Qrad0[l,m,n] = Qrad0[l-lrad,m-mrad,n-nrad]
    Qrad[l,m,n] = Qrad[l,m,n]+Qrad0[l,m,n]*exp(-tau[l,m,n])
  
  endfor
  endfor
  endfor

endif

end
