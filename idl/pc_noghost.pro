; $Id: pc_noghost.pro,v 1.1 2004-05-07 14:38:29 mee Exp $
;
;  Trim variable of ghost zones and empty dimensions
;
;  Author: Anders Johansen (ajohan@astro.ku.dk)
;  $Date: 2004-05-07 14:38:29 $
;  $Revision: 1.1 $
;
;  07-may-04/anders: coded
;
;
function pc_noghost,var,dim=dim,proc=proc
;
; Size of dimensions in input array
;
if (n_elements(dim) ne 1) then pc_read_dim,obj=dim,proc=proc

return, reform(var[dim.l1:dim.l2,dim.m1:dim.m2,dim.n1:dim.n2])

end
