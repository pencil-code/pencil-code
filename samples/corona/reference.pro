; Settings:

default, tolerance, 1d-10
default, reference_file, 'reference.xdr'

HR_ohm = pc_get_quantity ('HR_ohm', 'var.dat')

test_mean = mean (HR_ohm)
test_min = min (HR_ohm)
test_max = max (HR_ohm)

restore, reference_file

success = 1

if (abs (test_mean / ref_mean - 1.0) gt tolerance) then begin
	print, "ERROR: mean differs: ", test_mean, " != ", ref_mean
	success = 0
end

if (abs (test_min / ref_min - 1.0) gt tolerance) then begin
	print, "ERROR: min differs: ", test_min, " != ", ref_min
	success = 0
end

if (abs (test_max / ref_max - 1.0) gt tolerance) then begin
	print, "ERROR: max differs: ", test_max, " != ", ref_max
	success = 0
end

if (success) then print, "IDL autotest: >SUCCESS<"

END

