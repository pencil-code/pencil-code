;$Id: means.pro,v 1.2 2004-06-08 13:20:35 brandenb Exp $
function means,f,i
if n_elements(i) eq 0 then return,total(f)/n_elements(f)
;
s=size(f)
return,total(f,i)/s[i]
end
