function div,f
COMPILE_OPT IDL2,HIDDEN
return,xder(f[*,*,*,0])+yder(f[*,*,*,1])+zder(f[*,*,*,2])
end
