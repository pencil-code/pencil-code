;$Id$
function means,f,i
if n_elements(i) eq 0 then return,total(f)/n_elements(f)
;
s=size(f)
return,total(f,i)/s[i]
end
