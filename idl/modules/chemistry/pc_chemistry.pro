;
; Calculate useful parameters related to chemistry
;
pro pc_chemistry,f,molfrac=molfrac,massfrac=massfrac,concentration=concentration, tot_concentration=tot_concentration, PP=PP, mean_molar_mass=mean_molar_mass, tot_massfrac=tot_massfrac, debug=debug
;
@./data/pc_constants.pro
@./data/index.pro
  pc_read_dim,obj=dim,/quiet
  l1=dim.l1 & l2=dim.l2
  m1=dim.m1 & m2=dim.m2
  n1=dim.n1 & n2=dim.n2
  pc_read_param,obj=start,/quiet
  dimx=size(f.x)
  dimy=size(f.y)
  dimz=size(f.z)
  mx=dimx[1]
  my=dimy[1]
  mz=dimz[1]
  chem1=min(ichemspec)
  chem2=max(ichemspec)
  nchemspec=chem2-chem1+1
  if (keyword_set(debug)) then begin
     print,"l1,l2,m1,m2,n1,n2=",l1,l2,m1,m2,n1,n2
     print,specname
     print,specmass
     print,"mx,my,mz=,",mx,my,mz
     print,"chem1,chem2,nchemspec=",chem1,chem2,nchemspec
     print,"is_defined(PP)=",is_defined(PP)
     print,"keyword_set(PP)=",keyword_set(PP)
  endif
;
; Set dependencies
;
lmolfrac=0 & lconcentration=0 & ltot_concentration=0
lPP=0 & lmean_molar_mass=0
if keyword_set(molfrac) then begin
   lmolfrac=1 & lmean_molar_mass=1
endif
if keyword_set(mean_molar_mass) then begin
   lmean_molar_mass=1
endif
if keyword_set(concentration) then begin
   lconcentration=1 & lmolfrac=1 & lmean_molar_mass=1
endif
if is_defined(PP) then begin
;   if keyword_set(PP) then begin
   lPP=1 & lconcentration=1 & lmolfrac=1 & lmean_molar_mass=1
endif
print,"lPP=",lPP
;
; Put mass fractions in appropriate array
;
massfrac=fltarr(mx,my,mz,nchemspec)
for i=0,nchemspec-1 do begin
   com="massfrac(*,*,*,"+str(i)+")=f.yy"+str(1+i)
   ret=execute(com)
end
print,"min,max(massfrac)=",min(massfrac(l1:l2,m1:m2,n1:n2,*)),max(massfrac(l1:l2,m1:m2,n1:n2,*))
tot_massfrac=fltarr(mx,my,mz)
tot_massfrac=total(massfrac,4)
print,"min,max(tot_massfrac)=",min(tot_massfrac(l1:l2,m1:m2,n1:n2)),max(tot_massfrac(l1:l2,m1:m2,n1:n2))
;
; Calculate mean molar mass
;
if (lmean_molar_mass) then begin
   mean_molar_mass=fltarr(mx,my,mz)*0.0
   for i=0,mx-1 do begin
      for j=0,my-1 do begin
         for k=0,mz-1 do begin
            for ichem=0,nchemspec-1 do begin
               mean_molar_mass(i,j,k) = mean_molar_mass(i,j,k)+massfrac(i,j,k,ichem)/specmass(ichem)
            end
         end
      end
   end
   mean_molar_mass=1./mean_molar_mass
   print,"min,max(mean_molar_mass)=",min(mean_molar_mass(l1:l2,m1:m2,n1:n2)),max(mean_molar_mass(l1:l2,m1:m2,n1:n2))
endif
;
; Calculate mole fraction
;
if (lmolfrac) then begin
   molfrac=fltarr(mx,my,mz,nchemspec)*0.0
   for i=0,mx-1 do begin
      for j=0,my-1 do begin
         for k=0,mz-1 do begin
            for ichem=0,nchemspec-1 do begin
               molfrac(i,j,k,ichem) = massfrac(i,j,k,ichem)*mean_molar_mass(i,j,k)/specmass(ichem)
            end
         end
      end
   end
   print,"min,max(molfrac)=",min(molfrac(l1:l2,m1:m2,n1:n2,*)),max(molfrac(l1:l2,m1:m2,n1:n2,*))
endif
;
; Calculate concentration for each species
;
if (lconcentration) then begin
   if start.ldensity_nolog then begin
      print,"linear density"
      rho=f.rho
   endif else begin
      print,"logarithmic density"
      rho=exp(f.lnrho)
   end
   print,"min,max(rho)=",min(rho(l1:l2,m1:m2,n1:n2)),max(rho(l1:l2,m1:m2,n1:n2))
   concentration=fltarr(mx,my,mz,nchemspec)*0.0
   for i=0,mx-1 do begin
      for j=0,my-1 do begin
         for k=0,mz-1 do begin
            for ichem=0,nchemspec-1 do begin
               concentration(i,j,k,ichem) = molfrac(i,j,k,ichem)*rho(i,j,k)/mean_molar_mass(i,j,k)
            end
         end
      end
   end
   print,"min,max(concentration)=",min(concentration(l1:l2,m1:m2,n1:n2,*)),max(concentration(l1:l2,m1:m2,n1:n2,*))
;
; Calculate total concentration
;
   tot_concentration=fltarr(mx,my,mz)
   tot_concentration=total(concentration,4)
   print,"min,max(tot_concentration)=",min(tot_concentration(l1:l2,m1:m2,n1:n2)),max(tot_concentration(l1:l2,m1:m2,n1:n2))
endif
;
; Calculate pressure
;
if (lPP) then begin
   if start.ltemperature_nolog then begin
      print,"linear temperature"
      TT=f.tt
   endif else begin
      print,"logarithmic temperature"
      TT=exp(f.lntt)
   end
   print,"min,max(TT)=",min(TT(l1:l2,m1:m2,n1:n2)),max(TT(l1:l2,m1:m2,n1:n2))
   if (start.unit_system eq "cgs") then begin
      Rgas=8.314e7              ;erg/mol/K 
   endif else if (start.unit_system eq "si") then begin
      Rgas=8.314e7              ;J/mol/K
   endif else begin
      print,"No such unit system!"
      return
   end
   
   PP=fltarr(mx,my,mz,nchemspec)*0.0
   for i=0,mx-1 do begin
      for j=0,my-1 do begin
         for k=0,mz-1 do begin
            PP(i,j,k) = Rgas*tot_concentration(i,j,k)*TT(i,j,k)
         end
      end
   end
   print,"min,max(PP)=",min(PP(l1:l2,m1:m2,n1:n2)),max(PP(l1:l2,m1:m2,n1:n2))
endif


END



