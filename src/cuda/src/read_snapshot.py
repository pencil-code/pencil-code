import numpy as np

MOST_RECENT_SNAPSHOT = 7423174123783277124123809

def read_grid_info(datafile):
	#READ NX NY NZ
	nninf = open(datafile, "r")

	line = nninf.readline()
	line = line.split()

	NX = int(line[0])
	NY = int(line[1])
	NZ = int(line[2])
	GRID_SIZE = int(line[3])
	COMP_DOMAIN_SIZE_X = int(line[4])   
	COMP_DOMAIN_SIZE_Y = int(line[5])   
	COMP_DOMAIN_SIZE_Z = int(line[6])   
	PAD_SIZE = int(line[7])   
	BOUND_SIZE = int(line[8])   
	DX = float(line[9])   
	DY = float(line[10])   
	DZ = float(line[11])   
	time = float(line[12])

	nninf.close()

        #TODO: Should also contain the location of coordinate origin point. 
	print " NX =", NX, " NY =", NY, " NZ =", NZ, "GRID_SIZE =", GRID_SIZE
	print " COMP_DOMAIN_SIZE_X =", COMP_DOMAIN_SIZE_X, " COMP_DOMAIN_SIZE_Y =",COMP_DOMAIN_SIZE_Y, " COMP_DOMAIN_SIZE_Z =",COMP_DOMAIN_SIZE_Z
	print " PAD_SIZE =",PAD_SIZE, " BOUND_SIZE =",BOUND_SIZE, " DX =",DX, " DY =",DY, " DZ =",DZ, " time =",time
	print " NX = COMP_DOMAIN_SIZE_X + 2*PAD_SIZE"
	print " NY = COMP_DOMAIN_SIZE_Y + 2*BOUND_SIZE"
	print " NZ = COMP_DOMAIN_SIZE_Z + 2*BOUND_SIZE"
	print " GRID_SIZE = NX*NY*NZ"

	return NX, NY, NZ, GRID_SIZE, COMP_DOMAIN_SIZE_X, COMP_DOMAIN_SIZE_Y, COMP_DOMAIN_SIZE_Z, PAD_SIZE, BOUND_SIZE, DX, DY, DZ, time

def read_grid_dat(filename, NX, NY, NZ, PAD_SIZE, BOUND_SIZE, COMP_DOMAIN_SIZE_X, COMP_DOMAIN_SIZE_Y, COMP_DOMAIN_SIZE_Z, include_ghost, ldouble = 0):
	#Read and reconstruct the binary data

	print "Reading " + filename + " from memory"

	datafile = open(filename, "rb")

        if ldouble:
        	dtype = np.dtype(np.double) #double precision
	else:
		dtype = np.dtype(np.single) #C-type floating point number
	all_numbers = np.fromfile(datafile, dtype) #Read the binary file

	datafile.close()

	#Copy data into cube form
	
	data_cube = np.zeros((NX, NY, NZ))

	for i in range(0,NX): 
		for j in range(0,NY): 
			for k in range(0,NZ): 
				ind = i + j*NX + k*NX*NY
				data_cube[i,j,k] = all_numbers[ind]

	#Clear useless array to free memory
	all_numbers = []

	print "File " + filename + " Read"

	if include_ghost==0:
		#Remove ghost and padding zones
		grid = np.zeros((COMP_DOMAIN_SIZE_X, COMP_DOMAIN_SIZE_Y, COMP_DOMAIN_SIZE_Z))

		xtop = NX-PAD_SIZE
		ytop = NY-BOUND_SIZE
		ztop = NZ-BOUND_SIZE

		xbot = PAD_SIZE
		ybot = BOUND_SIZE
		zbot = BOUND_SIZE

		grid = data_cube[xbot:xtop,ybot:ytop,zbot:ztop]
		print "Ghost zones and padding removed."

		#Clear useless array to free memory
		data_cube = []
	else:
		#Leaves the ghost and padding zones untouched for debugging purposes
		grid = data_cube

	return grid

def read_whole_grid(data_dir="data/", snapshot=MOST_RECENT_SNAPSHOT, include_ghost=0, 
                    include_lnrho=1, include_ux=1, include_uy=1, include_uz=1, ldouble = 0):
	#Reads all grids from the prefered datafiles

	#Convert the snapshot number into a string
	if snapshot == MOST_RECENT_SNAPSHOT:
		snapshot_str = ""		
	else:
		snapshot_str = str(snapshot)

	#Read grid info
	grid_info_path = data_dir + "grid_info.ac" #TODO: time variable needs to be handled differently if we have multiple snapshots
	NX, NY, NZ, GRID_SIZE, COMP_DOMAIN_SIZE_X, COMP_DOMAIN_SIZE_Y, COMP_DOMAIN_SIZE_Z, PAD_SIZE, BOUND_SIZE, DX, DY, DZ, time = read_grid_info(grid_info_path)



	#Create the axis vectors 
	if include_ghost==0:
		xx = np.arange((COMP_DOMAIN_SIZE_X))
		yy = np.arange((COMP_DOMAIN_SIZE_Y))
		zz = np.arange((COMP_DOMAIN_SIZE_Z))
	else:
		CX_TOP = NX - PAD_SIZE
		CY_TOP = NY - BOUND_SIZE
		CZ_TOP = NZ - BOUND_SIZE
		CX_BOT = PAD_SIZE
		CY_BOT = BOUND_SIZE
		CZ_BOT = BOUND_SIZE

		xx = np.zeros((NX))
		yy = np.zeros((NY))
		zz = np.zeros((NZ))

		xx[CX_BOT:CX_TOP] = np.arange((COMP_DOMAIN_SIZE_X))
		yy[CY_BOT:CY_TOP] = np.arange((COMP_DOMAIN_SIZE_Y))
		zz[CZ_BOT:CZ_TOP] = np.arange((COMP_DOMAIN_SIZE_Z))

	#Scale the axis vectors TODO: Place the origin!
	xx = xx*DX
	yy = yy*DY
	zz = zz*DZ

	#Define grid paths
	lnrho_path = data_dir + "lnrho" + snapshot_str + ".dat"
	uu_x_path = data_dir + "uu_x" + snapshot_str + ".dat"
	uu_y_path = data_dir + "uu_y" + snapshot_str + ".dat"
	uu_z_path = data_dir + "uu_z" + snapshot_str + ".dat"

	#Read grid data
	if include_lnrho:
		lnrho = read_grid_dat(lnrho_path, NX, NY, NZ, PAD_SIZE, BOUND_SIZE, COMP_DOMAIN_SIZE_X, COMP_DOMAIN_SIZE_Y, COMP_DOMAIN_SIZE_Z, include_ghost, ldouble = ldouble)
	else:
		lnrho = np.zeros((xx.size, yy.size, zz.size))

	if include_ux:
		uu_x  = read_grid_dat(uu_x_path, NX, NY, NZ, PAD_SIZE, BOUND_SIZE, COMP_DOMAIN_SIZE_X, COMP_DOMAIN_SIZE_Y, COMP_DOMAIN_SIZE_Z, include_ghost, ldouble = ldouble)
        else:
                uu_x = np.zeros((xx.size, yy.size, zz.size))

	if include_uy:
		uu_y  = read_grid_dat(uu_y_path, NX, NY, NZ, PAD_SIZE, BOUND_SIZE, COMP_DOMAIN_SIZE_X, COMP_DOMAIN_SIZE_Y, COMP_DOMAIN_SIZE_Z, include_ghost, ldouble = ldouble)
        else:
                uu_y = np.zeros((xx.size, yy.size, zz.size))

	if include_uz:
		uu_z  = read_grid_dat(uu_z_path, NX, NY, NZ, PAD_SIZE, BOUND_SIZE, COMP_DOMAIN_SIZE_X, COMP_DOMAIN_SIZE_Y, COMP_DOMAIN_SIZE_Z, include_ghost, ldouble = ldouble)
        else:
                uu_z = np.zeros((xx.size, yy.size, zz.size))

	uu_tot = np.sqrt( np.power(uu_x,2.0) + np.power(uu_y,2.0) + np.power(uu_z,2.0) )

	return xx, yy, zz, lnrho, uu_x, uu_y, uu_z, uu_tot, time

def check_inf_nan(lnrho, uu_x, uu_y, uu_z, uu_tot, snapshot=MOST_RECENT_SNAPSHOT):
	#Calculate some simple diagnostics and check for inf or nan in the array

	print "\nCalculating diagnostics..."

	rho_min = np.nanmin(np.exp(lnrho))
	rho_max = np.nanmax(np.exp(lnrho))
	rho_rms = np.sqrt(np.mean(np.square(np.exp(lnrho))))

	lnrho_min = np.nanmin(lnrho)
	lnrho_max = np.nanmax(lnrho)
	lnrho_rms = np.sqrt(np.mean(np.square(lnrho)))

	uu_tot_min = np.nanmin(uu_tot)
	uu_tot_max = np.nanmax(uu_tot)
	uu_tot_rms = np.sqrt(np.mean(np.square(uu_tot)))

	uu_x_min = np.nanmin(uu_x)
	uu_x_max = np.nanmax(uu_x)
	uu_x_rms = np.sqrt(np.mean(np.square(uu_x)))

	uu_y_min = np.nanmin(uu_y)
	uu_y_max = np.nanmax(uu_y)
	uu_y_rms = np.sqrt(np.mean(np.square(uu_y)))

	uu_z_min = np.nanmin(uu_z)
	uu_z_max = np.nanmax(uu_z)
	uu_z_rms = np.sqrt(np.mean(np.square(uu_z)))
	
	if snapshot == MOST_RECENT_SNAPSHOT:
		print "The most recent snapshot: "
	else:
		print "For snapshot number: " + str(snapshot)

	print "rho_min = %e ; rho_max = %e ; rho_rms = %e" % (rho_min, rho_max, rho_rms)
	print "lnrho_min = %e ; lnrho_max = %e ; lnrho_rms = %e" % (lnrho_min, lnrho_max, lnrho_rms)
	print "uu_tot_min = %e ; uu_tot_max = %e ; uu_tot_rms = %e" % (uu_tot_min, uu_tot_max, uu_tot_rms)
	print "uu_x_min = %e ; uu_x_max = %e ; uu_x_rms = %e" % (uu_x_min, uu_x_max, uu_x_rms)
	print "uu_y_min = %e ; uu_y_max = %e ; uu_y_rms = %e" % (uu_y_min, uu_y_max, uu_y_rms)
	print "uu_z_min = %e ; uu_z_max = %e ; uu_z_rms = %e" % (uu_z_min, uu_z_max, uu_z_rms)

	print "\nChecking the presence of nan or inf.."

	lnrho_nan = np.nanmax(np.isnan(lnrho))
	lnrho_inf = np.nanmax(np.isinf(lnrho))
	rho_nan = np.nanmax(np.isnan(np.exp(lnrho)))
	rho_inf = np.nanmax(np.isinf(np.exp(lnrho)))
	uu_tot_nan = np.nanmax(np.isnan(uu_tot))
	uu_tot_inf = np.nanmax(np.isinf(uu_tot))
	uu_x_nan = np.nanmax(np.isnan(uu_x))
	uu_x_inf = np.nanmax(np.isinf(uu_x))
	uu_y_nan = np.nanmax(np.isnan(uu_y))
	uu_y_inf = np.nanmax(np.isinf(uu_y))
	uu_z_nan = np.nanmax(np.isnan(uu_z))
	uu_z_inf = np.nanmax(np.isinf(uu_z))

	print "Number 1 signifies the presence of nan or inf:"
	print "lnrho_nan = %i ; lnrho_inf = %i ; rho_nan = %i ; rho_inf = %i ; uu_tot_nan = %i ; uu_tot_inf = %i" % (lnrho_nan, lnrho_inf, rho_nan, rho_inf, uu_tot_nan, uu_tot_inf)
	print "uu_x_nan = %i ; uu_x_inf = %i ; uu_y_nan = %i ; uu_y_inf = %i ; uu_z_nan = %i ; uu_z_inf = %i" % (uu_x_nan, uu_x_inf, uu_y_nan, uu_y_inf, uu_z_nan, uu_z_inf)



