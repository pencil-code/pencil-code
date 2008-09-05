;$Id$
function pc_accumulate,f,n1
;
;  Computes the accumulative mean of f(t)
;  Used to show the convergence of Lyapunov exponents, etc.
;  The onset point for averaging can be choosen (n1)
;  23-mar-91
;
if n_elements(n1) eq 0 then n1=0
n2=n_elements(f)-1
a=f(n1:n2)
for i=n1,n2 do a(i-n1)=total(f(n1:i))/(i-n1+1)
return,a
end
