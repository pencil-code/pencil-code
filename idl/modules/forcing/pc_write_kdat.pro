;
;  $Id$
;
;  Program to calculate all possible vectors k that have a specified
;  length (in 1, 2 or 3 dimensions). The output file, k.dat, can be used
;  directly as input to forced turbulence simulations with the Pencil
;  Code.
;
pro pc_write_kdat, kmean, deltak, ndim, kmin=kmin, format=format, $
    screen=screen, file=file
;
;  Print to screen and file by default.
;
default, screen, 1
default, file, 1
;
;  Format of the output.
;
default, format, '(9f8.1)'
;
;  Smallest wavenumber is kmin. If box size is 2*pi, then kmin=1.0 is fine.
;
default, kmin, 1.0
;
;  If kmin is not an integer, one must be careful and get enough significant
;  digits on the k vectors.
;
if (kmin ne fix(kmin)) then begin
  print, 'Warning: kmin is not an integer, so please check that the k vectors have enough significant digits'
  if (format eq '(9f8.1)') then format='(9f13.7)'
endif
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
;  Information on how many k's were found
;
print, 'Found '+strtrim(n_elements(kk),2)+' vectors in the specified range.'
print, n_elements(kk), klen*kmin
;
;  Print to screen.
;
if (ndim le 2) then qq=fltarr(n_elements(kk))
if (ndim le 1) then pp=fltarr(n_elements(kk))
if (screen) then begin
  print, kk*kmin, format=format
  print, pp*kmin, format=format
  print, qq*kmin, format=format
endif
;
;  Write to file.
;
if (file) then begin
  close, 1
  openw, 1, 'k.dat'
    printf, 1, n_elements(kk), klen*kmin
    printf, 1, kk*kmin, format=format
    printf, 1, pp*kmin, format=format
    printf, 1, qq*kmin, format=format
  close, 1
endif
;
end
