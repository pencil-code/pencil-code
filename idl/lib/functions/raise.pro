FUNCTION raise,a,b
;
;  compute a^b such that its sign is negative for negative values of a
;
x=abs(a)^b
i=where(a lt 0.)
if i(0) ne -1 then x(i)=-x(i)
;
return,x
END
