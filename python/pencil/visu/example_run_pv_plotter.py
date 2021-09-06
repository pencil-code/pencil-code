#
# Example running script for the pv_plotter.py to create visualizations from slices
#
# Author: Leevi Veneranta (leevi.veneranta@aalto.fi)
#
import pencil as pc
import numpy as np

# Coordinate system
CARTESIAN = False
SPHERICAL = True
CYLINDER  = False

# Read in vector data or not
VECTORS = True

# Choose vector data
if VECTORS:
    vectors = ['uu1', 'uu2', 'uu3']
    # vectors = ['bb1', 'bb2', 'bb3']
else:
    vectors = None

# Which scalars to plot
FIELDS = ['lnrho']
# Normalization applied on the scalar
NORM   = 'linear'

# Choose time interval to plot
### 1. Start and end indeces (between 0 <= t <= slices.t[-1])
ISTART = 30
IEND   = 60
### 2. Start and end time instants i.e. the time values found in slices.t
TSTART = None
TEND   = None

################################################################################
##### PARAMETERS ###############################################################
################################################################################
## TALLBOX CARTESIAN PARAMETERS
if CARTESIAN:
    coordinates   = 'cartesian'

    ##### TALL BOX - these work for tallbox ####################################
    xn, yn, zn    = 128, 256, 768                   
    DIR           = '../4pc2OloqIM_data'
    CAMERA_CENTRE = (-850,-850, 1460)
    FOCAL_POINT   = (2*xn, 1.3*yn, zn/5)
    COLOR_RANGE   = None #[-50, 50] 
        
################################################################################
## SPHERICAL MILLENIUM PARAMETERS
elif SPHERICAL:
    DIR           = '/home/leevi/Desktop/MOUNT/testfield_millennium/data'
    coordinates   = 'spherical'
    CAMERA_CENTRE = (-2.9,-2.5, 2.3)
    #CAMERA_CENTRE = (-10,-10,8)
    TEND          = 40095
    FOCAL_POINT   = None
    offset        = None
    window_size   = None
    COLOR_RANGE   = None

################################################################################
else: # CYLINDER
    DIR             = './data'
    coordinates     = 'cylinder'
    CAMERA_CENTRE   = None
    TEND            = 1e30
    FOCAL_POINT     = None
    offset          = None
    window_size     = None
    COLOR_RANGE     = None

################################################################################

# Set the PlotSettings object - see pv_plotter.py PlotSetting for further documentation
# on the possible parameters
settings = pc.visu.pv_plotter.PlotSettings(
    ### General settings
    off_screen      = True, 
    preview         = False, 
    window_size     = (656, 1008), 
    progress_bar    = True, 
    videoformat     = None, 
    imageformat     = 'png', 
    framerate       = 7,
    figdir          = './images/',
    moviedir        = './movies/', 
    bg_color        = 'white',
    timestamp       = True,
    
    ### Plot parameters
    norm            = NORM, 
    coordinates     = coordinates,
    opacities       = {'xy': 1, 'xy2': 1, 'xz': 1, 'yz': 1},
    
    ### Coordinates specific parameters
    offset          = 6., 
    spherical_xz_pos= -1, 
    
    ### Vector plot parameters
    n_vectors       = 350,
    vector_size     = 3,
    vector_method   = 'every_nth', # 'random' | 'every_nth'
    vector_scaling  = 'magnitude',
    vector_max      = 11,
    surface_vectors = True, 
    
    ### Streamline parameters
    streamlines         = True,
    stream_tube_radius  = 0.5,
    stream_variable_radius = 'Magnitude',
    stream_radius_factor = 8,
    stream_src_points   = {'xy': 40, 'xy2': 40, 'xz': 150, 'yz': 150},
    stream_show_source  = False,
    stream_src_radius   = None, #!TODO! MISSING
    stream_log_scale = False,
    stream_params       = {
        'max_steps'             : 2000,
        'max_time'              : 1e60,
        'terminal_speed'        : 1e-60,
        'integration_direction' : 'both',
        'compute_vorticity'     : False,
        'integrator_type'       : 4,
        }, 

    ### Random sampling parameters 
    set_seed_1          = False,
    constant_seed       = True,
    
    ### Camera parameters
    camera_centre   = CAMERA_CENTRE, 
    focal_point     = FOCAL_POINT,
    
    ### Axes parameters
    show_axes       = True, 
    axes_labels     = ('','',''),   #('x','y','z'), 
    axes_font       = 9, 
    
    ##### Colorbar properties
    ### Field specific
    show_field_sbar = True,
    field_cmap      = 'spring',
    field_sbar_title = 'bb\n',
    field_sbar_pos_x = 0.05,
    field_sbar_pos_y = None,
    
    ### Scalar specific
    show_scalar_sbar= True,
    scalar_cmap     = 'jet', 
    scalar_sbar_title = 'log[rho]\n',
    scalar_sbar_pos_x = 0.89,
    scalar_sbar_pos_y = None, 
    
    ### General colorbar properties
    cbar_width = 0.06,
    cbar_height = 0.65,
    cbar_title_font = 10,
    cbar_label_font = 8,
    n_colors        = 256, 
    title_position  = 'upper_right', 
    title_font      = 14,
    str_unit        = ' UNIT', 
    tscale          = 1., 
    time_precision  = 4,   
    )

# Call plotter
pc.visu.pv_plotter.plot(
    # Data
    slice_obj=None, datadir=DIR, precision='f',
    fields=FIELDS,
    #xyzplane=xyzplane,
    vectors=vectors,
    # Plot interval
    tstart=TSTART, tend=TEND, islice=-1,
    istart=ISTART, iend=IEND,
    # Colors & scaling
    color_range=COLOR_RANGE, color_levels=None,
    unit='unit_velocity', rescale=1., par=list(),
    # Extras
    debug=True, settings=settings
    )

################################################################################