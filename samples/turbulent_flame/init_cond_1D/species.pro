; General settings:
default, zoom, 4
default, specie, 'H2'

; Process videofiles
default, videofiles_processed, 0
if (not videofiles_processed) then begin
	print, "Reading all videofiles."
	spawn, 'src/read_all_videofiles.x'
	videofiles_processed = 1
end

; Visualize selected species
@pc_visualize_species

