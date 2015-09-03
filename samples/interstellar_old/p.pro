; IDL Version 6.0 (linux x86 m32)
; Journal File for mee@wormhole
; Working directory: /home/data/mee/cvs-src/branches/eos/samples/interstellar
; Date: Thu Mar 24 17:42:43 2005
savep=!p 
pc_read_var,obj=data,variables=['tt','pp','cs2'],/add,/magic,/trimall
;  name               minval         maxval          mean           rms
; uu_x     -->       -15.6781        15.6781    2.89508E-11        1.58480
; uu_y     -->       -15.6783        15.6780    1.72630E-07        1.58480
; uu_z     -->       -15.1612        17.2703      0.0151254        1.61486
; lnrho    -->      -0.374380       0.309380      0.0481443       0.214412
; ss       -->        2.35894        5.97678        2.72103        2.80595
; shock    -->        0.00000      0.0244843     0.00197302     0.00594523
; tt       -->        7001.07    2.76456E+06        84998.7        383710.
; pp       -->        35.9513        18690.5        570.575        2620.37
; cs2      -->        62.1057        24524.1        754.014        3403.85
; t =   7.41688e-05
; % Program caused arithmetic error: Floating underflow
loadct,3
contour,data.lnrho[16,*,*],data.y,data.z,/fill,nlevels=256,/iso,title='Density (colorscale) and Velocity (arrows)'
wdvelovect,data.uu[16,*,*,1],data.uu[16,*,*,2],data.y,data.z,/over

!p=savep
end
