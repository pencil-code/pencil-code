pro esrg_legend,labels,pos=pos,linestyle=linestyle,psym=psym,clrbox=clrbox,$
           color=color,thick=thick,fg_color=fg_color,bg_color=bg_color,$
           box=box,numsym=numsym,ystretch=ystretch,silent=silent,$
           center=center,fill=fill,right=right,norm=norm,spos=spos, $
           lenfac=lenfac,charsize=charsize,xstretch=xstretch
;+
; ROUTINE: legend
;
; PURPOSE: draws a plot legend 
;
; USEAGE:  legend,labels,pos=pos,linestyle=linestyle,psym=psym,clrbox=clrbox,$
;                 color=color,thick=thick,fg_color=fg_color,bg_color=bg_color,$
;                 box=box,numsym=numsym,ystretch=ystretch,silent=silent,$
;                 center=center,fill=fill,right=right,norm=norm,spos=spos,$
;                 charsize=charsize,xstretch=xstretch
;               
; INPUT:
;
; labels	
;     an array of labels (required). Elements of LABELS which do not
;     begin with the characters ".l",".f",".c", or ".r" are treated as
;     legend labels.  These elements are drawn to the right of a line
;     or symbol of given LINESTYLE, PSYM, COLOR or THICK type. The
;     total number of labels in LABELS (i.e., not counting legend
;     titles) should correspond to the number of elements in those
;     input parameters.
;     
;     If any of the label strings start with the characters
;     ".l",".f",".c", or ".r" that element of the LABELS array is
;     treated as a legend title.  Legend titles appear as left
;     justified (".l"), filled (".f") centered (".c") or right
;     justified (".r") strings in the legend box and do not correspond
;     to elements of LINESTYLE, PSYM, COLOR or THICK.  Null strings
;     ("") or blank strings (" ") are also treated as legend
;     titles. Lines or symbols are not drawn to the right of legend
;     titles.
;     
;     If none of the keywords, LINESTYPE, PSYM, COLOR or THICK are
;     set, then all elements of LABELS will be treated as legend
;     titles and will be left justified unless prepended by ".f", ".c"
;     or ".r".  or unless one of the keywords keywords "center",
;     "fill" "right" are set.  These keywords set the default, string
;     justification but may always be overridden by one of the
;     justification flags.
;     
;     Consider,
;     
;            labels=['.ctitle','1','2','.lsecond title,'a','b']
;            linestyle= [1,2,3,4]
;     
;     In this example LABELS contains 4 legend labels and 2 legend
;     titles. The correspondence with the linestyle input parameter is
;     as follows:
;     
;                      labels(1) &#60;==> linestyle(0)
;                      labels(2) &#60;==> linestyle(1)
;                      labels(4) &#60;==> linestyle(2)
;                      labels(5) &#60;==> linestyle(3)
;     
;     To simplify input, LABELS may also be specified as a single
;     string with the different LABEL elements separated by
;     backslashes (\) (E.g., labels='.cCloud height\\1\2\3\4\5'
;
; KEYWORD INPUT:
;
; pos
;     position and size of legend area in normalized data 
;     coordinates.  Format of this four element vector is
;     [x0,y0,x1,y1] where x0,y0 is the position of lower
;     left corner and x1,y1 is the position of the upper
;     right corner.  For example pos=[.5,.5,.9,.9] specifies
;     a legend in the upper right quadrant of the data 
;     window.  If the POS parameter is not set or if POS is
;     set to zero, then the CURBOX procedure is called to
;     set the legend position interactively.  
;     
;     NOTE: the value of POS can be retrieved by setting POS
;     to a named variable which has an initial value of zero.
;
; linestyle   
;     an array of linestyle types, one for each label in LABELS
;
; psym	
;     an array of symbol types, one for each label in LABELS. 
;     PSYM can be either a integer or string array.  If PSYM is a
;     string array then the library routine USERSYMBOL is used to
;     create the symbols.  In this case if the symbol specifier is
;     appended with '_f' then a filled symbol will be generated
;
; clrbox
;     if set, a color-filled box appears to the left of each legend
;     label.  The number of boxes drawn can be changed with NUMSYM.
;     If clrbox is set to 2, a color-filled circle appears to the left
;     of each legend label.
;
; color       
;     an array of color indices, one for each label in LABELS.
;     Any elements of COLOR set to -1 causes the default color,
;     !P.COLOR, to be used.
;
; thick       
;     an array of line thicknesses, one for each label in LABELS
;     (default=!p.thick)
;
; numsym      
;     number of symbol nodes used to indicate linestyle or symbol
;     type.  The length of the line is measured in units of LENFAC
;     character widths so that the length of the line
;     = LENFAC*(NUMSYM-1) * X_CHARSIZE
;     (default=2 when linestyle, color or thick set, otherwise default=1)
;
; lenfac
;     factor which controls length of line between symbols (see numsym)
;     (default=5)
;
; fg_color    
;     color of box and legend titles (default=!P.COLOR)
;
; bg_color    
;     background color. Setting BG_COLOR erases the area covered by
;     the legend (filling it with color BG_COLOR) prior to writing the
;     legend.  If both BG_COLOR and !p.color are zero then the
;     background color is reset to 255 to gaurantee a readability.
;               
; box         
;     if set draw a box around the legend text
;
; silent      
;     if not set, print box position string to the terminal
;     and show popup help on the CURBOX cursor routine
;     
;     silent=1 don't print position string
;     
;     silent=2 don't print position string, don't show help widget
;
; center
;     if set, default action is to center legend labels
;
; fill
;    if set, default action is to fill legend labels
;
; right
;    if set, default action is to right justify labels
;
; norm
;    if set, normalized coordinates are used to specify POS both
;    on input and printed output.  This option allows placement
;    of the legend outside of a plot region. 
;
; ystretch    
;     factor by which to increase vertical spacing of legend 
;     entries. This parameter is particularly useful when used
;     with the SPOS keyword.  (default = 1)  
;     
;     NOTE: the aspect ratio of the legend box can be
;     modified on the fly by pushing the cursor box against
;     one of the window boundaries and pressing the middle
;     button.  
;
; spos
;    string value which specifies legend position to one of a set of
;    standard positions.  if spos is set, it supercedes the value of pos
;
;                'bl' = bottom left
;                'br' = bottom right
;                'tl' = top left
;                'tr' = top right
;                'top'= top center
;                'bot'= bottom center
;
; charsize
;    character size used when legend position specified with SPOS.  
;    CHARSIZE has no effect when POS is used to set position.
;
; xstretch
;    Adjust length of labels by this factor -- necessary with
;    postscript output, where lenstr() does not return the correct result.
;
; OUTPUT:       none
;
; PROCEDURE:
;    When called without the POS keyword, the legend position and size
;    is set with the CURBOX routine.  The legend is positioned by
;    dragging the box where you want the legend to appear.  The size
;    of the legend area can be decreased/increased by the left/middle
;    mouse buttons.  When the right mouse button is pressed the legend
;    is drawn and a numerical positioning string giving the current
;    value of the POS keyword is written to the terminal (nothing
;    printed if SILENT is set).  You can run LEGEND in batch mode by
;    pasting this value of POS into the LEGEND command line.  The best
;    way to get good-looking legends on your hardcopies is to size
;    your graphics window to the approximate shape of the output
;    media.  For example a plot which is to be printed in landscape
;    mode should be previewed on a window which is approximately 1.4
;    times as wide as it is tall.
;
;    NOTE:    The values returned for the POS keyword are based on a
;    computation of the length of the text strings in your legend.  If
;    you change the contents of the legend titles or if you change the
;    default text font, you must rerun LEGEND in interactive mode to
;    determine new values for the POS paramter.
;                     
;
;
;; EXAMPLE       interactive mode (put the box anywhere you want and press
;;               the right mouse button)
;
;   dcolors
;   plot,6*randf(20,3),/nodata 
;   for i=1,6 do oplot,i*randf(20,3),li=i-1,color=i
;   lb='.cFirst bunch\\First\Second\Third\\.cSecond bunch\\forth\fifth\sixth'
;   legend,lb,li=[0,1,2,3,4,5],/box,bg=0,color=[1,2,3,4,5,6]
;
;; EXAMPLE       interactive mode. retrieve the value of POS for later calls:
;
;   !p.multi=[0,1,2]
;   plot,[0,.4],yrange=[0,1],li=0 & oplot,[0,.6],li=2 & oplot,[0,.9],li=3
;   legpos=0
;   lb=['.cLegend Demo','','one','two','three']
;   legend,lb,pos=legpos,li=[0,2,3]
;   plot,[0,.4],yrange=[0,1],li=0 & oplot,[0,.6],li=2 & oplot,[0,.9],li=3
;   legend,lb,pos=legpos,li=[0,2,3]
;   !p.multi=0
;
;
;; EXAMPLE       use fill mode to print a figure caption
;
;  w11x8
;  !p.multi=[0,1,2]
;  plot,dist(20)
;  t=$
;   'When called without the POS keyword, the legend position and size is\'+$
;   'set with the CURBOX routine.  The legend is positioned by dragging the\'+$
;   'box where you want the legend to appear.  The size of the legend area\'+$
;   'can be decreased/increased by the left/middle mouse buttons.  When the\'+$
;   'right mouse button is pressed the legend is drawn and a numerical\'+$
;   'positioning string giving the current value of the POS keyword is\'+$
;   'written to the terminal (nothing printed if SILENT is set).'
;
;   legend,t,bg=0,pos=[0.00,-0.52,0.47,-0.18]         ; default left justified
;   legend,t,bg=0,/right,pos=[0.53,-0.52,1.00,-0.18]  ; right justified
;   legend,t,bg=0,/center,pos=[0.27,-1.00,0.74,-0.66] ; centered
;
;; NOTE: procedure PFILL provides more elaborate text formatting capability
;   
;; EXAMPLE       batch mode:
;
;   plot,[-50,50],[0,1.5],/nodata & x=findgen(101)-50.
;   li=indgen(6) & pos=[0.66,0.54,0.91,0.89] 
;   for i=0,5 do oplot,x,1./(1.+(x/(i+2))^2),li=li(i) 
;   labels='.cPlot Key\\First\Second\Third\Fourth\Fifth\Sixth'
;   legend,' ',bg=80,pos=pos+.02,/box
;   legend,labels,li=li,pos=pos,bg=40,/box
;                
;; EXAMPLE       batch mode with symbols generated by USERSYMBOL:
;
;   plot,[-50,50],[0,1.5],/nodata & x=findgen(101)-50.
;   psym=['TRIANGLE','DIAMOND','PENTAGON','CIRCLE','SQUARE','SPADE']
;   pos=[0.66,0.54,0.91,0.89] 
;   for i=0,5 do begin &$
;     usersymbol,psym(i) &$
;     oplot,x,1./(1.+(x/(i+2))^2),psym=8,symsize=2 &$
;   endfor
;   labels='.cPlot Key\\First\Second\Third\Fourth\Fifth\Sixth'
;   legend,' ',bg=4,pos=pos+.02,/box
;   legend,labels,psym=psym,pos=pos,bg=10,/box
;                
;
;  author:  Paul Ricchiazzi                            4dec92
;           Institute for Computational Earth System Science
;           University of California, Santa Barbara
;
; REVISIONS
; 15dec92: added legend titles, and CHARSIZE parameter
; 11jan93: added numsym parameter
; 20jan93: added thick parameter
;  2feb93: fixed positioning bug  legend titles
; 25mar93: removed the NOBOX option, now you set BOX to get a box
; 29mar93: added the bg_color option to blank out legend area before writing
;  8apr93: use CURBOX for cursor control and LENSTR for exact string size
; 27apr93: improved alignment symbol and label (ylab=yy-.3*charsize*cnvrty)
;  9jun93: center legend labels when legend titles are longer (see dxcen)
; 17jun93: added ystretch option to increase vertical spacing
; 17jun93: default line thickness read from !p.thick
; 30sep93: .l implied when LINESTYLE, PSYM, COLOR, THICK  not set. see NOLINES
; 28Jun94: LABELS is now normal input param, not a keyword.
; 28Jun94: legend "titles" don't correspond to LINESTYLE, PSYM, COLOR, THICK
;          vector elements; no need to put in dummy values. see examples
; 18Aug94: added USYM option
; 28Jun95: added the .f format option
;  8Sep95: added CLRBOX option
;  5oct95: added charsize adjustment for Y axis
;   sep96: added spos keyword, numsym default = 1 unless linestyle or color set
;   nov96: CLRBOX=2 yields filled circles instead of boxes
;   may03: wd: added xstretch keyword
;-
on_error,2

n=n_elements(labels)
if n eq 0 then begin & xhelp,'legend' & return & end
if n eq 1 then labls=str_sep(labels,'\') else labls=labels
n=n_elements(labls)
jblnks=where(labls eq ' ',nblnks) 

if (n_elements(xstretch) eq 0) then xstretch=1

case 1 of
  keyword_set(center):djust='c'
  keyword_set(right): djust='r'
  keyword_set(fill):  djust='f'
  else:               djust='l'
end
if nblnks gt 0 then labls(jblnks)='.'+djust+' '
if keyword_set(silent) eq 0 then silent=0
;
nls=n_elements(linestyle)
nps=n_elements(psym)
ncl=n_elements(color)
ntk=n_elements(thick) 

nolines=nls eq 0 and nps eq 0 and ncl eq 0 and ntk eq 0
;

if n_elements(fg_color) eq 0 then fgc=!p.color else fgc=fg_color
if n_elements(ystretch) eq 0 then ystretch=1.

titles=strarr(n)
just=strarr(n)

; real label length = slenmx - 2*(number of exclamation marks)

nlbls=0
for i=0,n-1 do begin
  ss=labls(i)
  if strmid(ss,0,1) eq '.' then begin
    just(i)=strmid(ss,1,1)
    tlen=strlen(ss)-2
    ttt=strmid(ss,2,tlen)
    if ttt eq '' then titles(i)=' ' else titles(i)=ttt
  endif else begin
    if nolines or ss eq '' then begin
      just(i)=djust
      if ss eq '' then titles(i)=' ' else titles(i)=ss
    endif else begin
      nlbls=nlbls+1
    endelse
  endelse
endfor

if n_elements(numsym) eq 0 then begin
  if keyword_set(clrbox) or (nls eq 0 and ncl eq 0) then numsym=1 else numsym=2
endif 

if nlbls gt 0 then begin

  drawline= ( nls gt 0 or ncl gt 0 or ntk gt 0 ) and numsym gt 1
  
  if nls eq 0 then begin & linestyle=replicate(-1,nlbls)   & nls=nlbls & end
  if nps eq 0 then begin & psym=intarr(nlbls)              & nps=nlbls & end 
  if ncl eq 0 then begin & color=replicate(-1,nlbls)       & ncl=nlbls & end 
  if ntk eq 0 then begin & thick=replicate(!p.thick,nlbls) & ntk=nlbls & end

  nof='Number of '
  mess=' does not match number of labels'
  if nls ne nlbls then  message,nof+'linestyles'+mess,/continue
  if nps ne nlbls then  message,nof+'psym types'+mess,/continue
  if ncl ne nlbls then  message,nof+'colors'+mess,/continue
  if ntk ne nlbls then  message,nof+'thicknesses'+mess,/continue
  if nls ne nlbls or nps ne nlbls or ncl ne nlbls or ntk ne nlbls then return
endif

cnvrtx=float(!d.x_ch_size)/!d.x_vsize   ; conversion factor from
cnvrty=float(!d.y_ch_size)/!d.y_vsize   ; characters to normalized coordinates

iti=where(titles eq '',niti)            ; these labels are not legend titles

if not keyword_set(lenfac) then lenfac=4

if niti gt 0 then $
  slenmx=max(xstretch*lenstr(labls(iti)))+(lenfac*(numsym-1)+2)*cnvrtx else slenmx=0

tlenmx=max(xstretch*lenstr(titles))
xsize=(slenmx > tlenmx) + 2*cnvrtx       ; add two character widths
dxcen=.5*(xsize-2*cnvrtx-slenmx)         ; offset to center legend labels

ysize=1.1*(n+1)*cnvrty

if max([!x.window,!y.window],min=mn) eq mn then $
  plot,[0,1],[0,1],xstyle=4,ystyle=4,/nodata

sx0=!x.window(0)
sx1=!x.window(1)-sx0
sy0=!y.window(0)
sy1=!y.window(1)-sy0

if not ( keyword_set(pos) or keyword_set(spos)) then begin
;  xsz=xsize*!d.x_size
;  ysz=ysize*!d.y_size*ystretch
;  if silent le 1 then curbox,x,y,xsz,ysz,/message else $
;                      curbox,x,y,xsz,ysz

  x0=.5*(1.-xsize)
  x1=.5*(1.+xsize)
  y0=.5*(1.-ysize*ystretch)
  y1=.5*(1.+ysize*ystretch)
  if silent le 1 then curbox,x0,x1,y0,y1,/init,/message else $
                      curbox,x0,x1,y0,y1,/init
  if keyword_set(norm) then begin
    pos=[x0,y0,x1,y1]
  endif else begin
    pos=[(x0-sx0)/sx1,(y0-sy0)/sy1,(x1-sx0)/sx1,(y1-sy0)/sy1]
  endelse
  posstring=string(form='(a,4(f10.2,a))',$
           ',pos=[',pos(0),',',pos(1),',',pos(2),',',pos(3),']')

  if keyword_set(silent) eq 0 then print,strcompress(posstring,/remove_all)

endif else begin
  
  if keyword_set(spos) then begin
    if not keyword_set(norm) then begin
      xsz=xsize/sx1
      ysz=ysize/sy1
    endif else begin
      xsz=xsize
      ysz=ysize
    endelse
    if not keyword_set(charsize) then charsize=!p.charsize
    if charsize eq 0 then charsize=1
    xsz=xsz*charsize
    ysz=ysz*charsize
    case spos of
      "bl" : begin & x0=.05 & x1=x0+xsz & y0=.05 & y1=y0+ysz*ystretch & end
      "br" : begin & x1=.95 & x0=x1-xsz & y0=.05 & y1=y0+ysz*ystretch & end
      "tl" : begin & x0=.05 & x1=x0+xsz & y1=.95 & y0=y1-ysz*ystretch & end
      "tr" : begin & x1=.95 & x0=x1-xsz & y1=.95 & y0=y1-ysz*ystretch & end
      "top": begin 
              x0=.5*(1-xsz) & x1=.5*(1+xsz)
              y1=.98        & y0=y1-ysz
             end
      "bot": begin 
              x0=.5*(1-xsz) & x1=.5*(1+xsz)
              y0=.02        & y1=y0+ysz
             end
      else : message,'check value of spos'
    endcase

    if not keyword_set(norm) then begin
      x0=sx0+x0*sx1         ;                    
      x1=sx0+x1*sx1         ; normalized coordinates of
      y0=sy0+y0*sy1         ; legend box         
      y1=sy0+y1*sy1         ;                    
    endif

  endif else begin
    if n_elements(pos) ne 4 then message,'Incorrect POS specification'
    if keyword_set(norm) then begin
      x0=pos(0)
      x1=pos(2)
      y0=pos(1)
      y1=pos(3)
    endif else begin
      x0=sx0+pos(0)*sx1         ;                    
      x1=sx0+pos(2)*sx1         ; normalized coordinates of
      y0=sy0+pos(1)*sy1         ; legend box         
      y1=sy0+pos(3)*sy1         ;                    
    endelse
  endelse

endelse

; blank out legend area

if n_elements(bg_color) ne 0 then begin
  if !p.color eq 0 and bg_color eq 0 then bgc=255 else bgc=bg_color
  polyfill,[x0,x1,x1,x0],[y0,y0,y1,y1],color=bgc,/normal
endif
if keyword_set(box) then begin
  xbox=[x0,x1,x1,x0,x0]
  ybox=[y0,y0,y1,y1,y0]
  plots,xbox,ybox,color=fgc,/norm
endif
;

xcharsize=(x1-x0)/xsize       ; find the largest charsize that will fit
ycharsize=(y1-y0)/ysize       ; in the legend box
if( xcharsize lt ycharsize ) then begin
  chrsz=xcharsize
endif else begin
  chrsz=ycharsize
endelse

symsize=chrsz
dx=chrsz*cnvrtx
dy=1.1*cnvrty*(y1-y0)/ysize ; dy=1.1*chrsz*cnvrty

;
if numsym eq 1 then begin
  xb=x0+2*dx+dxcen
  xbar=[xb,xb]
  xlab=x0+2*dx+dxcen+(lenfac*(numsym-1)+2)*dx
endif else begin
  xbar=x0+dx+dxcen+lenfac*(numsym-1)*dx*findgen(numsym)/(numsym-1)
  xlab=x0+dx+dxcen+(lenfac*(numsym-1)+2)*dx
endelse

yy=y1-dy

jlbl=0
stdspc=xstretch*lenstr(' ')*chrsz

for i=0,n-1 do begin
  ylab=yy-.3*chrsz*cnvrty

  if titles(i) eq '' then begin

    clr=color(jlbl)
    if clr eq -1 then clr=!p.color
    ls=linestyle(jlbl)
    ps=psym(jlbl)
    
    if (size(ps))(1) eq 7 then begin
      usersymbol,ps
      ps=8
    endif    
    thk=thick(jlbl)
    ybar=replicate(yy,numsym)

    if keyword_set(clrbox) then begin
      if clrbox eq 2 then begin
        xcbox=1.3*cos(findrng(0,360,37)*!dtor)
        ycbox=1.3*sin(findrng(0,360,37)*!dtor)
      endif else begin
        xcbox=1.3*[-1,1,1,-1,-1]+1
        ycbox=1.3*[-1,-1,1,1,-1]
      endelse
      usersym,xcbox,ycbox,/fill,color=clr
      plots,xbar,ybar,psym=8,linestyle=ls,color=clr,/norm,symsize=symsize
      usersym,xcbox,ycbox
      plots,xbar,ybar,psym=8,linestyle=ls,/norm,symsize=symsize,thick=thk
    endif else begin
      if drawline then begin
        plots,xbar,ybar,psym=ps,linestyle=ls,color=clr,/norm,$
            symsize=symsize,thick=thk
      endif else begin
        if ps ne 0 then plots,max(xbar),max(ybar),psym=ps,color=clr,/norm,$
            symsize=symsize,thick=thk
      endelse
    endelse

    xyouts,xlab,ylab,labls(i),color=fgc,charsize=chrsz,/norm
    jlbl=jlbl+1

  endif else begin

    case just(i) of
      'l':begin & align=0.0 & xtitle=x0+dx      & end
      'f':begin & align=0.0 & xtitle=x0+dx      & end
      'c':begin & align=0.5 & xtitle=.5*(x0+x1) & end
      'r':begin & align=1.0 & xtitle=x1-dx      & end
    endcase
;    stop
    if just(i) eq 'f' then begin
      tstr=str_sep(titles(i),' ')
      lngth=xstretch*lenstr(tstr)*chrsz
      nwords=n_elements(tstr)
      
      dspc=(x1-x0-2*dx-total(lngth))/((nwords-1)>1)
      if dspc gt 4*stdspc then dspc=2*stdspc
      for j=0,nwords-1 do begin 
        xyouts,xtitle,ylab,tstr(j),color=fgc,charsize=chrsz,/norm 
        xtitle=xtitle+dspc+lngth(j)
      endfor
    endif else begin  
      xyouts,xtitle,ylab,titles(i),color=fgc,charsize=chrsz,$
           alignment=align,/norm
    endelse
  endelse
  yy=yy-dy
endfor

return
end
