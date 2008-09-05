function levi_civita,i,j,k,debug=debug
if keyword_set(debug) then print,'$Id$'
;
;  totally antisymmetric tensor
;  20-mar-04/axel: adapted from fortran code
;
if (i eq 0 and j eq 1 and k eq 2) or $
   (i eq 1 and j eq 2 and k eq 0) or $
   (i eq 2 and j eq 0 and k eq 1) then begin
    value=+1.
end else if $
   (i eq 2 and j eq 1 and k eq 0) or $
   (i eq 0 and j eq 2 and k eq 1) or $
   (i eq 1 and j eq 0 and k eq 2) then begin
    value=-1.
end else begin
    value=0.
endelse
;
return,value
end
