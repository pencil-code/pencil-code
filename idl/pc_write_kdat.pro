;
;  $Id: pc_write_kdat.pro,v 1.1 2007-11-13 07:09:49 ajohan Exp $
;
;  Program to calculate all possible vector k that have a specified
;  length (in 1, 2 or 3 dimensions). The output file, k.dat, can be used
;  directly as input to forced turbulence simulations with the Pencil
;  Code.
;
pro pc_write_kdat, kmean, deltak, ndim
;
;  Chosen interval in k.
;
k0=kmean-deltak
k1=kmean+deltak
;
;  Initialise mean length of k vector.
;
klen=0.0
;
;  Go through all possible k and save combinations within chosen length
;  interval.
;
if (ndim eq 1) then begin
  for ik=-ceil(k1),ceil(k1) do begin
    k=sqrt(ik^2)
    if ( (k ge k0) and (k le k1) ) then begin
      if (n_elements(kk) eq 0) then kk=ik else kk=[kk,ik]
      klen=klen+k
    endif
  endfor
endif
;
if (ndim eq 2) then begin
  for ik=-ceil(k1),ceil(k1) do begin
    for ip=-ceil(k1),ceil(k1) do begin
      k=sqrt(ik^2+ip^2)
      if ( (k ge k0) and (k le k1) ) then begin
        if (n_elements(kk) eq 0) then kk=ik else kk=[kk,ik]
        if (n_elements(pp) eq 0) then pp=ip else pp=[pp,ip]
        klen=klen+k
      endif
    endfor
  endfor
endif
;
if (ndim eq 3) then begin
  for ik=-ceil(k1),ceil(k1) do begin
    for ip=-ceil(k1),ceil(k1) do begin
      for iq=-ceil(k1),ceil(k1) do begin
        k=sqrt(ik^2+ip^2+iq^2)
        if ( (k ge k0) and (k le k1) ) then begin
          if (n_elements(kk) eq 0) then kk=ik else kk=[kk,ik]
          if (n_elements(pp) eq 0) then pp=ip else pp=[pp,ip]
          if (n_elements(qq) eq 0) then qq=iq else qq=[qq,iq]
          klen=klen+k
        endif
      endfor
    endfor
  endfor
endif
;
;  Mean length of k vectors.
;
klen=klen/n_elements(kk)
;
;  Write to file.
;
if (ndim le 2) then qq=fltarr(n_elements(kk))
if (ndim le 1) then pp=fltarr(n_elements(kk))
close, 1
openw, 1, 'k.dat'
  printf, 1, n_elements(kk), klen
  printf, 1, kk, format='(9i8)'
  printf, 1, pp, format='(9i8)'
  printf, 1, qq, format='(9i8)'
close, 1
;
end
