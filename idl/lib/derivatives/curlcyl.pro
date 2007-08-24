;;; Calculate the curl of a 3-d vector field
function curlr,f,xx
COMPILE_OPT IDL2,HIDDEN
 return,yder(f[*,*,*,2])/xx -zder(f[*,*,*,1])
end
function curlphi,f,xx
COMPILE_OPT IDL2,HIDDEN
 return,zder(f[*,*,*,0])-xder(f[*,*,*,2])
end
function curlz,f,xx
COMPILE_OPT IDL2,HIDDEN
 return,xder(f[*,*,*,1])+f[*,*,*,1]/xx-yder(f[*,*,*,0])/xx
end
function curlcyl,f,x
COMPILE_OPT IDL2,HIDDEN
 w=make_array(size=size(f),/nozero)
 tmp=size(f)
 my=tmp[2]
 mz=tmp[3]
 xx=spread(x,[1,2],[my,mz])
 w[*,*,*,0]=curlr(f,xx)
 w[*,*,*,1]=curlphi(f,xx)
 w[*,*,*,2]=curlz(f,xx)
 return,w
end
