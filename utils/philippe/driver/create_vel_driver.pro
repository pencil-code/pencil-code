N_times = 128
delta_t = 2.0
switch_on_time = 0.0
data_size_x = 96
data_size_y = 96
image_filename = "../../utils/philippe/driver/KFU-logo.png"
driver_filename = "vel-z"

read_png, image_filename, image
if (size (image, /n_dimensions) eq 3) then image = reform (image[0,*,*])
image = reform (image[*,*])
image -= min (image)
image /= double (max (image))
image -= 0.5

image_size = size (image, /dimensions)
image_size_x = image_size[0]
image_size_y = image_size[1]

data = dblarr (data_size_x, data_size_y)
data[0:image_size_x-1,0:image_size_y-1] = image
data = spread (data, 2, N_times)

diff_x = data_size_x - image_size_x
diff_y = data_size_y - image_size_y
diff = min ([diff_x, diff_y])
for pos = 0, diff do begin
	data[*,*,pos] = 0.0
	data[pos:pos+(image_size_x-1), pos:pos+(image_size_y-1), pos] = image
	tvscl, data[*,*,pos]
	wait, 0.1
end
for pos = diff+1, 2*diff do begin
	space = 2*diff - pos
	data[*,*,pos] = 0.0
	data[space:space+(image_size_x-1), space:space+(image_size_y-1), pos] = image
	tvscl, data[*,*,pos]
	wait, 0.1
end

times = dblarr (N_times)
for pos = 0, N_times-1 do begin
	times[pos] = switch_on_time + delta_t * pos
end

openw, lun_data, driver_filename+".dat", /get_lun
openw, lun_times, driver_filename+"_times.dat", /get_lun
for pos = 0, N_times-1 do begin
	writeu, lun_data, reform (data[*,*,pos])
	writeu, lun_times, times[pos]
end
close, lun_data
close, lun_times
free_lun, lun_data
free_lun, lun_times

END
