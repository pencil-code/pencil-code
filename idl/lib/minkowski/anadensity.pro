
th=[1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,20.]
s = size(th)
Nth = s(1)
T=findgen(Nth)
W=findgen(Nth)
L=findgen(Nth)
pc_read_var,obj=obj,varfile='VAR535',variables=['oo'],/magic,/trimall
for ith=0,Nth-1 do begin
  print,th[ith]
  isostruct,obj.oo,th[ith],object=M
  T[ith]=M.V0/(2.*M.V1)
  W[ith]=2.*M.V1/(!Pi*M.V2)
  L[ith]=3.*M.V2/(4.*M.V3)
endfor
WINDOW,0
plot, th,T,PSYM=2
WINDOW,1
plot, th,W,PSYM=2
WINDOW,2
plot, th,L,PSYM=2
END
