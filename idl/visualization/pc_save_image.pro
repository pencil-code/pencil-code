;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_save_image.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Routine to save a copy of a given window ID or of the current window
;;;   into a given graphics file format. The image data can also be retured.
;;;
;;;  Parameters:
;;;   file       filename including suffix
;;;   window     window ID, otherwise the current window is taken
;;;   image      a buffer with the image data is returned
;;;   quality    quality setting for lossy formats, like JPEG
;;;   grey       keyword that chooses greyscale conversion
;;;   crop       array of crop coordinates for cropping the image data


; Saves an image from a given window ID of from the current window
pro pc_save_image, file, window=window, image=image, quality=quality, grey=grey, crop=crop

	if (keyword_set (window)) then win = window else win = !D.window
	if (win lt 0) then begin
		print, 'Window number is invalid or no window is open.'
		stop
	end

	if (keyword_set (grey)) then begin
		true_color = 1 - grey
	end else begin
		true_color = 1
	end

	; try to determine image format from file extension
	img_type = strlowcase (strmid (file, strlen (file)-3, 3))

	if (img_type eq "jpg") then begin
		default, quality, 98
		image = tvrd (true=true_color)
		if (keyword_set (crop)) then image = image[*,crop[0]:crop[1],crop[2]:crop[3]]
		WRITE_JPEG, file, image, /true, /progressive, quality=quality
	end else if (img_type eq "gif") then begin
		image = tvrd ()
		if (keyword_set (crop)) then image = image[*,crop[0]:crop[1],crop[2]:crop[3]]
		tvlct, r, g, b, /get
		WRITE_GIF, file, image, r, g, b
	end else if (img_type eq "ppm") then begin
		image = tvrd (true=true_color, /order)
		if (keyword_set (crop)) then image = image[*,crop[0]:crop[1],crop[2]:crop[3]]
		WRITE_PPM, file, image
	end else if (img_type eq "png") then begin
		image = tvrd (true=true_color)
		if (keyword_set (crop)) then image = image[*,crop[0]:crop[1],crop[2]:crop[3]]
		WRITE_PNG, file, image
	end else begin
		print, 'Image file type "'+img_type+'" is unknown.'
		stop
	end

end

