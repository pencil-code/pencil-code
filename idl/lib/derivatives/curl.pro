;;; Calculate the curl of a 3-d vector field
function curlx,f
COMPILE_OPT IDL2,HIDDEN
 return,yder(f[*,*,*,2])-zder(f[*,*,*,1])
end
function curly,f
COMPILE_OPT IDL2,HIDDEN
 return,zder(f[*,*,*,0])-xder(f[*,*,*,2])
end
function curlz,f
COMPILE_OPT IDL2,HIDDEN
 return,xder(f[*,*,*,1])-yder(f[*,*,*,0])
end
function curl,f
COMPILE_OPT IDL2,HIDDEN
 w=make_array(size=size(f),/nozero)
 w[*,*,*,0]=curlx(f)
 w[*,*,*,1]=curly(f)
 w[*,*,*,2]=curlz(f)
 return,w
end
