;$Id: means.pro,v 1.1 2004-06-08 13:18:23 brandenb Exp $
function means,f,i
if n_elements(i) eq 0 then return,total(f)/n_elements(f)
;
s=size(f)
return,total(f,i)/s(i)
end
