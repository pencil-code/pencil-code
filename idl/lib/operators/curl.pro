function curl,f,debug=debug,ibcx=ibcx,xfull=xfull
compile_opt IDL2,HIDDEN
if keyword_set(debug) then print,'$Id: curl.pro,v 1.2 2004-05-05 17:10:32 mee Exp $'
;
s=size(f)
if keyword_set(xfull) then s[1]=n_elements(xfull)
w=make_array(size=s,/nozero)
;
if n_elements(ibcx) eq 0 then begin
  w[*,*,*,0]=yder(f[*,*,*,2])-zder(f[*,*,*,1])
  w[*,*,*,1]=zder(f[*,*,*,0])-xder(f[*,*,*,2])
  w[*,*,*,2]=xder(f[*,*,*,1])-yder(f[*,*,*,0])
end else begin
  if keyword_set(xfull) then begin
    w[*,*,*,0]=yder(f[*,*,*,2])-zder(f[*,*,*,1])
    w[*,*,*,1]=zder(f[*,*,*,0])-xder(f[*,*,*,2],ibcx=ibcx[2],xfull=xfull)
    w[*,*,*,2]=xder(f[*,*,*,1],ibcx=ibcx[1],xfull=xfull)-yder(f[*,*,*,0])
  end else begin
    w[*,*,*,0]=yder(f[*,*,*,2])-zder(f[*,*,*,1])
    w[*,*,*,1]=zder(f[*,*,*,0])-xder(f[*,*,*,2],ibcx=ibcx[2])
    w[*,*,*,2]=xder(f[*,*,*,1],ibcx=ibcx[1])-yder(f[*,*,*,0])
  end
end
;
return,w
end
