;  $Id$
;
;  Calculate auxiliary particle variables such as particle velocity.
;
;  Author: Anders Johansen
;
function pc_particles_aux, np=np, vvpsum=vvpsum, var=var, $
    dim=dim, datadir=datadir

  if (n_elements(dim) eq 0) then pc_read_dim, obj=dim, datadir=datadir


  if (var eq 'vvp') then begin

    result=fltarr(dim.mx,dim.my,dim.mz,3)

    for l=dim.l1,dim.l2 do begin
      for m=dim.m1,dim.m2 do begin
        for n=dim.n1,dim.n2 do begin
          if (np[l,m,n] ne 0.0) then result[l,m,n,*]=vvpsum[l,m,n,*]/np[l,m,n]
       endfor
      endfor
    endfor

  endif else begin
; Make sure and missing cases say so...
    message, 'pc_particles_aux: CASE NOT IMPLEMENTED'
  endelse

  return, result

end

