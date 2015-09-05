PRO SAVEIMAGE, FILE, BMP=BMP, PNG=PNG, PICT=PICT, JPEG=JPEG, TIFF=TIFF, $
  QUALITY=QUALITY, DITHER=DITHER, CUBE=CUBE, QUIET=QUIET, OUTSIZE=OUTSIZE

;+
; NAME:
;    SAVEIMAGE
;
; PURPOSE:
;    Save the current graphics window to an output file (GIF by default).
;
;    The output formats supported are:
;    GIF   8-bit with color table,
;    BMP   8-bit with color table,
;    PNG   8-bit with color table,
;    PICT  8-bit with color table,
;    JPEG 24-bit true color,
;    TIFF 24-bit true-color.
;
;    Any conversions necessary to convert 8-bit or 24-bit images onscreen to
;    8-bit or 24-bit output files are done automatically.
;
; CATEGORY:
;    Input/Output.
;
; CALLING SEQUENCE:
;    SAVEIMAGE, FILE
;
; INPUTS:
;    FILE     Name of the output file (GIF format by default).
;
; OPTIONAL INPUTS:
;    None.
;
; KEYWORD PARAMETERS:
;    BMP      Set this keyword to create BMP format (8-bit with color table).
;    PNG      Set this keyword to create PNG format (8-bit with color table).
;    PICT     Set this keyword to create PICT format (8-bit with color table).
;    JPEG     Set this keyword to create JPEG format (24-bit true color).
;    TIFF     Set this keyword to create TIFF format (24-bit true color).
;    QUALITY  If set to a named variable, specifies the quality for
;             JPEG output (default 75). Ranges from 0 ("terrible") to
;             100 ("excellent"). Smaller quality values yield higher
;             compression ratios and smaller output files.
;    DITHER   If set, dither the output image when creating 8-bit output
;             which is read from a 24-bit display (default is no dithering).
;    CUBE     If set, use the color cube method to quantize colors when
;             creating 8-bit output which is read from a 24-bit display
;             (default is to use the statistical method). This may improve
;             the accuracy of colors in the output image, especially white.
;    QUIET    Set this keyword to suppress the information message
;             (default is to print an information message).
;    OUTSIZE  Size of output image file in pixels
;
; OUTPUTS:
;    None.
;
; OPTIONAL OUTPUTS:
;    None
;
; COMMON BLOCKS:
;    None
;
; SIDE EFFECTS:
;    The output file is overwritten if it exists.
;
; RESTRICTIONS:
;    Requires IDL 5.0 or higher (square bracket array syntax).
;
; EXAMPLE:
;
;openr, lun, filepath('hurric.dat', subdir='examples/data'), /get_lun
;image = bytarr(440, 330)
;readu, lun, image
;free_lun, lun
;loadct, 13
;tvscl, image
;saveimage, 'hurric.gif'
;
; MODIFICATION HISTORY:
; Liam.Gumley@ssec.wisc.edu
; http://cimss.ssec.wisc.edu/~gumley
; $Id: saveimage.pro,v 1.17 2000/02/08 19:39:38 gumley Exp $
; Modified: 29-Oct-2000, zarro@smmdac.nascom.nasa.gov - added outsize keyword
; Version 3, 13-Aug-2003, William Thompson
;       Use SSW_WRITE_GIF instead of WRITE_GIF
; Version 4, 28-June-2006, Zarro - set DECOMPOSED=0 if undefined, and
;           reset to original value when completed.
;
; Copyright (C) 1999 Liam E. Gumley
;
; This program is free software; you can redistribute it and/or
; modify it under the terms of the GNU General Public License
; as published by the Free Software Foundation; either version 2
; of the License, or (at your option) any later version.
;
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.
;
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
;-

rcs_id = '$Id: saveimage.pro,v 1.17 2000/02/08 19:39:38 gumley Exp $'

;-------------------------------------------------------------------------------
;- CHECK INPUT
;-------------------------------------------------------------------------------

;- Check arguments
if (n_params() ne 1) then message, 'Usage: SAVEIMAGE, FILE'
if (n_elements(file) eq 0) then message, 'Argument FILE is undefined'
if (n_elements(file) gt 1) then message, 'Argument FILE must be a scalar string'

;- Check keywords
output = 'GIF'
if keyword_set(bmp)  then output = 'BMP'
if keyword_Set(png)  then output = 'PNG'
if keyword_set(pict) then output = 'PICT'
if keyword_set(jpeg) then output = 'JPEG'
if keyword_set(tiff) then output = 'TIFF'
if (n_elements(quality) eq 0) then quality = 100

;- Check for TVRD capable device
if ((!d.flags and 128)) eq 0 then message, 'Unsupported graphics device'

;- Check for open window
if (!d.flags and 256) ne 0 then begin
  if (!d.window lt 0) then message, 'No graphics windows are open'
endif

;- Get display depth
depth = 8
if (!d.n_colors gt 256) then depth = 24
entry_depth=depth

;-------------------------------------------------------------------------------
;- GET CONTENTS OF GRAPHICS WINDOW
;-------------------------------------------------------------------------------

;- Handle window devices (other than the Z buffer)
if (!d.flags and 256) ne 0 then begin

  ;- Copy the contents of the current display to a pixmap
  current_window = !d.window
  xsize = !d.x_size
  ysize = !d.y_size
  window, /free, /pixmap, xsize=xsize, ysize=ysize, retain=2
  device, copy=[0, 0, xsize, ysize, 0, 0, current_window]

  ;- Set decomposed color mode for 24-bit displays
  version = float(!version.release)
  if (depth gt 8) then begin
    if (version gt 5.1) then begin
     device, get_decomposed=entry_decomposed
     if n_elements(entry_decomposed) eq 0 then entry_decomposed=0
    endif
    device, decomposed=1
  endif

endif

;- Read the pixmap contents into an array

depth=8
if (depth gt 8) then begin
  image = tvrd(order=0, true=1)
endif else begin
  image = tvrd(order=0)
endelse

;-- resize output?

resize=0
if n_elements(outsize) ne 0 then begin
 sz=size(image)
 if sz(0) eq 3 then begin
  nx=sz(2) & ny=sz(3)
 endif else begin
  nx=sz(1) & ny=sz(2)
 endelse 
 new_nx=float(outsize(0))
 new_ny=new_nx
 if n_elements(outsize) gt 1 then new_ny=float(outsize(1))
 resize=((nx ne new_nx) or (ny ne new_ny)) and (new_nx gt 1) and (new_ny gt 1)
endif
 
;- Handle window devices (other than the Z buffer)
if (!d.flags and 256) ne 0 then begin

  ;- Restore decomposed color mode for 24-bit displays
  if (entry_depth gt 8) then begin
    if (version gt 5.1) then begin
      device, decomposed=entry_decomposed
    endif else begin
      device, decomposed=0
      if (keyword_set(quiet) eq 0) then $
        print, 'Decomposed color was turned off'
    endelse
  endif

  ;- Delete the pixmap
  wdelete, !d.window
  wset, current_window

endif

;- Get the current color table
tvlct, r, g, b, /get

;- If an 8-bit image was read, reduce the number of colors
if (depth le 8) then begin
  reduce_colors, image, index
  r = r[index]
  g = g[index]
  b = b[index]
endif

;-------------------------------------------------------------------------------
;- WRITE OUTPUT FILE
;-------------------------------------------------------------------------------

case 1 of

  ;- Save the image in 8-bit output format
  (output eq 'GIF')  or (output eq 'BMP') or $
  (output eq 'PICT') or (output eq 'PNG') : begin

    if (depth gt 8) then begin

      ;- Convert 24-bit image to 8-bit
      case keyword_set(cube) of
        0 : image = color_quan(image, 1, r, g, b, colors=256, $
              dither=keyword_set(dither))
        1 : image = color_quan(image, 1, r, g, b, cube=6)
      endcase

      ;- Sort the color table from darkest to brightest
      table_sum = total([[long(r)], [long(g)], [long(b)]], 2)
      table_index = sort(table_sum)
      image_index = sort(table_index)
      r = r[table_index]
      g = g[table_index]
      b = b[table_index]
      oldimage = image
      image[*] = image_index[temporary(oldimage)]

    endif

    ;- Save the image

    if resize then image=congrid(temporary(image),new_nx,new_ny)

    case output of
      'GIF'  : ssw_write_gif,  file, image, r, g, b
      'BMP'  : write_bmp,  file, image, r, g, b
      'PNG'  : write_png,  file, image, r, g, b
      'PICT' : write_pict, file, image, r, g, b
    endcase

  end

  ;- Save the image in 24-bit output format
  (output eq 'JPEG') or (output eq 'TIFF') : begin

    ;- Convert 8-bit image to 24-bit

    if (depth le 8) then begin
      info = size(image)
      nx = info[1]
      ny = info[2]
      true = bytarr(3, nx, ny)
      true[0, *, *] = r[image]
      true[1, *, *] = g[image]
      true[2, *, *] = b[image]
      image = temporary(true)
    endif

    ;- If TIFF format output, reverse image top to bottom
    if (output eq 'TIFF') then image = reverse(temporary(image), 3)

    ;- Write the image

    if resize then image=congrid(temporary(image),3,new_nx,new_ny)             
    case output of
      'JPEG' : write_jpeg, file, image, true=1, quality=quality
      'TIFF' : write_tiff, file, image, 1
    endcase

  end

endcase

;- Print information for the user
if (keyword_set(quiet) eq 0) then $
  print, file, output, format='("Created ",a," in ",a," format")'

END
