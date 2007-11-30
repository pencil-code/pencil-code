;; $Id: pc_noghost.pro,v 1.4 2007-11-30 13:59:22 ajohan Exp $
;;
;;  Trim variable of ghost zones and empty dimensions
;;
;;  Author: Anders Johansen (ajohan@astro.ku.dk)
;;  $Date: 2007-11-30 13:59:22 $
;;  $Revision: 1.4 $
;;
;;  07-may-04/anders: coded
;;
function pc_noghost,var,dim=dim,proc=proc,quiet=quiet
;
; Size of dimensions in input array
;
if (n_elements(dim) ne 1) then pc_read_dim,obj=dim,proc=proc,/quiet
;
varsize=size(var)
;
; Array must be given with ghost zones, or there is not much point in
; removing them.
;
if (varsize[1] ne dim.mx) or (varsize[2] ne dim.my) or $
   (varsize[3] ne dim.mz) then begin
  message, 'Attempted to remove ghosts from an array not of '+ $
           'size [mx,my,mz,*,*]... Skipping.', /info
  return,var
endif
;
; Remove ghost zones and empty dimensions.
;
return, reform(var[dim.l1:dim.l2,dim.m1:dim.m2,dim.n1:dim.n2,*,*,*])
;
end
