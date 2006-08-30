; $Id: pc_noghost.pro,v 1.3 2006-08-30 13:28:37 dintrans Exp $
;
;  Trim variable of ghost zones and empty dimensions
;
;  Author: Anders Johansen (ajohan@astro.ku.dk)
;  $Date: 2006-08-30 13:28:37 $
;  $Revision: 1.3 $
;
;  07-may-04/anders: coded
;
;
function pc_noghost,var,dim=dim,proc=proc,quiet=quiet,run2D=run2D
;
; Size of dimensions in input array
;
if (n_elements(dim) ne 1) then pc_read_dim,obj=dim,proc=proc,quiet=quiet

varsize=size(var)

if (not keyword_set(run2D)) then begin
  if (varsize[1] ne dim.mx) or (varsize[2] ne dim.my) or (varsize[3] ne dim.mz) then begin
    message,'Attempted to remove ghosts from an array not of size [mx,my,mz,*,*]... Skipping.',/info
    return,var
  endif
endif

if (not keyword_set(run2D)) then begin
  return, reform(var[dim.l1:dim.l2,dim.m1:dim.m2,dim.n1:dim.n2,*,*])
endif else begin
  if (dim.ny eq 1) then $
    return, reform(var[dim.l1:dim.l2,dim.n1:dim.n2,*,*]) else $
    return, reform(var[dim.l1:dim.l2,dim.m1:dim.m2,*,*])
endelse

end
