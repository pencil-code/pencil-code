nchemspec=0
nreactions=0
;
openr,1,'data/net_reactions.dat'
readf,1,nchemspec,nreactions
print,nchemspec,nreactions
net_react_p=fltarr(nchemspec,nreactions)
net_react_m=fltarr(nchemspec,nreactions)
fo='(8e10.2)'
;
it=0
while not eof(1) do begin
  readf,1,t
  readf,1,net_react_p,net_react_m,fo=fo
  print,t
  if it eq 0 then begin
    tt=t
    nnet_react_p=net_react_p
    nnet_react_m=net_react_m
  endif else begin
    tt=[tt,t]
    nnet_react_p=[nnet_react_p,net_react_p]
    nnet_react_m=[nnet_react_m,net_react_m]
  endelse
  it=it+1
endwhile
nt=it
;
net_react_p=reform(nnet_react_p,nchemspec,nt,nreactions)
net_react_m=reform(nnet_react_m,nchemspec,nt,nreactions)
help,net_react_p,net_react_m
close,1
end
