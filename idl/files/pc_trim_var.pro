; $Id: pc_trim_var.pro,v 1.1 2004-05-07 08:57:07 ajohan Exp $
;
;  Trim variable of ghost zones and empty dimensions
;
;  Author: Anders Johansen (ajohan@astro.ku.dk)
;  $Date: 2004-05-07 08:57:07 $
;  $Revision: 1.1 $
;
;  07-may-04/anders: coded
;
;
function pc_trim_var, var

  pc_read_dim, l1=l1, l2=l2, m1=m1, m2=m2, n1=n1, n2=n2

  return, reform(var(l1:l2,m1:m2,n1:n2,*))

end
