;***********************************************************************
function zder2,f
COMPILE_OPT IDL2,HIDDEN
;
print,'You have not forced compilation of a particular set of derivative routines.'
print,'This is accomplished by manually compiling routines from lib/derivatives.'
print,"Or automatically by running 'pc_init' or the '.r start' method depending upon your preference."
message,"ABORTING",/INFO
;
return,0.
end
