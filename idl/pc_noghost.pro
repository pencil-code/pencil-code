; $Id: pc_noghost.pro,v 1.2 2004-05-11 15:29:19 mee Exp $
;
;  Trim variable of ghost zones and empty dimensions
;
;  Author: Anders Johansen (ajohan@astro.ku.dk)
;  $Date: 2004-05-11 15:29:19 $
;  $Revision: 1.2 $
;
;  07-may-04/anders: coded
;
;
function pc_noghost,var,dim=dim,proc=proc,quiet=quiet
;
; Size of dimensions in input array
;
if (n_elements(dim) ne 1) then pc_read_dim,obj=dim,proc=proc,quiet=quiet

varsize=size(var)

if (varsize[1] ne dim.mx) or (varsize[2] ne dim.mx) or (varsize[3] ne dim.mz) then begin
  message,'Attempted to remove ghosts from an array not of size [mx,my,mz,*,*]... Skipping.',/info
  return,var
endif

return, reform(var[dim.l1:dim.l2,dim.m1:dim.m2,dim.n1:dim.n2,*,*])

end
