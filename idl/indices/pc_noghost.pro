;; $Id: pc_noghost.pro,v 1.5 2007-12-10 07:49:50 ajohan Exp $
;;
;;  Trim variable of ghost zones and empty dimensions
;;
;;  Author: Anders Johansen (ajohan@astro.ku.dk)
;;  $Date: 2007-12-10 07:49:50 $
;;  $Revision: 1.5 $
;;
;;  07-may-04/anders: coded
;;
function pc_noghost,var,dim=dim,proc=proc,run2D=run2D,quiet=quiet
;
; Size of dimensions in input array
;
if (n_elements(dim) ne 1) then pc_read_dim, obj=dim, proc=proc, /quiet
;
; For 2-D runs it is possible to write data without ghost zones in the
; missing direction.
;
run2Dxy=0 & run2Dxz=0
if (keyword_set(run2D)) then begin
  if (dim.nz eq 1) then run2Dxy=1
  if (dim.ny eq 1) then run2Dxz=1
endif
;
varsize=size(var)
;
; Array must be given with ghost zones, or there is not much point in
; removing them.
;
if (run2Dxy) then begin
  if ((varsize[1] ne dim.mx) or (varsize[2] ne dim.my)) then begin
    message, 'Attempted to remove ghosts from an array not of '+ $
             'size [mx,my,*,*]... Skipping.', /info
    return, var
  endif
endif else if (run2Dxz) then begin
  if ((varsize[1] ne dim.mx) or (varsize[2] ne dim.mz)) then begin
    message, 'Attempted to remove ghosts from an array not of '+ $
             'size [mx,mz,*,*]... Skipping.', /info
    return, var
  endif
endif else begin
  if ((varsize[1] ne dim.mx) or (varsize[2] ne dim.my) or $
      (varsize[3] ne dim.mz)) then begin
    message, 'Attempted to remove ghosts from an array not of '+ $
             'size [mx,my,mz,*,*]... Skipping.', /info
    return, var
  endif
endelse
;
; Remove ghost zones and empty dimensions.
;
if (run2Dxy) then begin
  return, reform(var[dim.l1:dim.l2,dim.m1:dim.m2,*,*,*])
endif else if (run2Dxz) then begin
  return, reform(var[dim.l1:dim.l2,dim.n1:dim.n2,*,*,*])
endif else begin
  return, reform(var[dim.l1:dim.l2,dim.m1:dim.m2,dim.n1:dim.n2,*,*,*])
endelse
;
end
