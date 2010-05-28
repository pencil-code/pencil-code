;$Id$
FUNCTION pc_round,x,d
;
;  rounds off a number to the accuracy of d
;
default,d,1.
return,fix(x/d+.5)*d
END
