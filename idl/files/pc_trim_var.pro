; $Id: pc_trim_var.pro,v 1.3 2004-05-11 08:45:06 ajohan Exp $
;
;  Trim variable of ghost zones and empty dimensions
;
;  Author: Anders Johansen (ajohan@astro.ku.dk)
;  $Date: 2004-05-11 08:45:06 $
;  $Revision: 1.3 $
;
;  07-may-04/anders: coded
;
;
function pc_trim_var, var
;
; Size of dimensions in input array
;
  mx_arr = n_elements(var(*,0,0,0))
  my_arr = n_elements(var(0,*,0,0))
  mz_arr = n_elements(var(0,0,*,0))
;
; Must consider number of dimensions of the input array
;
  if (my_arr eq 1 and mz_arr eq 1) then return, reform(var(3:mx_arr-4))
  if (my_arr ne 1 and mz_arr eq 1) then return, reform(var(3:mx_arr-4),var(3:my_arr-4))
;
; 3+D array (the fourth dimension can be velocity component, etc.)
;
  return, reform(var(3:mx_arr-4,3:my_arr-4,3:mz_arr-4,*,*))

end
