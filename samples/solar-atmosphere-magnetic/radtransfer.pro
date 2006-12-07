@data/pc_constants.pro
arad = 1.8049444E-16

dtau_thresh_min = 2.9730177E-02
dtau_thresh_max = 87.33655

Srad = arad*exp(4*lnTT)

Isurf = fltarr(nx,ny,3,3)

nrad = 1
nnstart=n1 & nnstop=n2 & nsign=+1

for lrad = -1,1 do begin
for mrad = -1,1 do begin

    llstart=l1 & llstop=l2 & ll1=l1 & ll2=l2 & lsign=+1
    mmstart=m1 & mmstop=m2 & mm1=m1 & mm2=m2 & msign=+1

  if (lrad>0) then begin
    llstart=l1 & llstop=l2 & ll1=l1-lrad & lsign=+1
  endif
  if (lrad<0) then begin
    llstart=l2 & llstop=l1 & ll2=l2-lrad & lsign=-1
  endif
  if (mrad>0) then begin
    mmstart=m1 & mmstop=m2 & mm1=m1-mrad &  msign=+1
  endif                                    
  if (mrad<0) then begin                   
    mmstart=m2 & mmstop=m1 & mm2=m2-mrad &  msign=-1
  endif
  if (nrad>0) then begin
    nnstart=n1 & nnstop=n2 & nn1=n1-nrad &  nsign=+1
  endif                                    
  if (nrad<0) then begin                   
    nnstart=n2 & nnstop=n1 & nn2=n2-nrad &  nsign=-1
  endif

  dlength = sqrt((dx*lrad)^2+(dy*mrad)^2+(dz*nrad)^2)
  
  tau = fltarr(mx,my,mz)
  Qrad = fltarr(mx,my,mz)
  
  for irep = 0,2 do begin
  
    for n = nnstart,nnstop,nsign do begin
    for m = mmstart,mmstop,msign do begin
    for l = llstart,llstop,lsign do begin
    
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
  
    if (lrad ne 0) then begin
      Qrad[llstart-lrad    ,mm1:mm2,n1:n2] = Qrad[    llstop,mm1:mm2,n1:n2]
       tau[llstart-lrad    ,mm1:mm2,n1:n2] =  tau[    llstop,mm1:mm2,n1:n2]
    endif
    if (mrad ne 0) then begin
      Qrad[ll1:ll2,mmstart-mrad    ,n1:n2] = Qrad[ll1:ll2,    mmstop,n1:n2]
       tau[ll1:ll2,mmstart-mrad    ,n1:n2] =  tau[ll1:ll2,    mmstop,n1:n2]
    endif
  
  endfor

  slice = reform((Qrad+Srad)[l1:l2,m1:m2,n2])

  print,lrad,mrad,min(slice),max(slice)
  contour,slice,/fill,nlevels=64

  Isurf[*,*,lrad+1,mrad+1] = slice

endfor
endfor

end
