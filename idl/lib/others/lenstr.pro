 function lenstr,str
COMPILE_OPT IDL2,HIDDEN
;+
; ROUTINE:         lenstr
; USEAGE:          result=lenstr(str)
;
; input:
;  str             a single string or string array.
; 
; output:
;  result          length of the string(s) in normalized units
;                  the number of elements of RESULT matches the number of 
;                  elements of STRING. 
;
; procedure:
;                  This function returns the physical length of the
;                  string on the output device, not the number of
;                  characters.  This is done by first switching to 'X'
;                  and writing the string(s) with XYOUTS in graphics
;                  mode 5, which disables display to the screen but
;                  does not interfere with operation of XYOUTS.  The
;                  WIDTH keyword parameter of XYOUTS is used to
;                  retrieve the physical length of the string(s).
;
;  Example:
;;; 
;     for i=(byte(' '))[0],(byte('z'))[0] do print,string(i),lenstr(string(i))
;;;
;
;  author:  Paul Ricchiazzi                            7apr93
;           Institute for Computational Earth System Science
;           University of California, Santa Barbara
;-

dsave=!d.name
psave=!p
!p.charsize=1.

set_plot,'X'

device,get_graphics=oldg,set_graphics=5
nn=n_elements(str)


case nn of

     0:w=0

     1:xyouts,0,0,/device,str,width=w

     else:begin
       w=fltarr(nn)
       for i=0,nn-1 do begin
         xyouts,0,0,/device,str[i],width=ww
         w[i]=ww
       endfor
     end
endcase

fac1=float(!d.x_ch_size)/!d.x_vsize   ; ratio of char width to device1 width

device,set_graphics=oldg
set_plot,dsave
!p=psave

fac2=float(!d.x_ch_size)/!d.x_vsize   ; ratio of char width to device2 width

return,w*fac2/fac1                    ; string width adjusted for device width
end


