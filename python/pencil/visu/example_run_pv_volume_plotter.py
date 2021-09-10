#
# Example running script for the pv_volume_plotter.py to create visualizations
#   from pencil code VAR data
#
# Author: Leevi Veneranta (leevi.veneranta@aalto.fi)
#
from pencil.visu.pv_volume_plotter import Pyvista3DPlot, Plot3DSettings


def runPyvistaVolumePlotter():    
    COORDINATES = 'cartesian'
    # COORDINATES = 'spherical'

    # DATADIR = '/home/leevi/Desktop/MOUNT/testfield_millennium/data'
    # DATADIR = '/home/leevi/Desktop/MOUNT/4pc20loqsMG3.5mc/data'
    DATADIR = '/home/leevi/Desktop/MOUNT/4pc2OloqlM/data'

    # SCALAR_KEY  = 'lnrho'
    SCALAR_KEY  = 'tt'
    # SCALAR_KEY  = 'rho'
    # SCALAR_KEY  = 'ss'
    # SCALAR_KEY = 'cooling'
    
    VECTOR_KEY  = 'bb'
    # VECTOR_KEY = ['ux', 'uy', 'uz']
    # VECTOR_KEY  = ['ax', 'ay', 'az']
    
    IVAR        = -1
    METHOD      = 'points'
    
    settings = Plot3DSettings(
        ##############################
        ### General
        title       = '   ',
        title_font_size = 14,
        n_points    = 2500,
        set_seed_1    = False,
        imageformat = 'png',
        off_screen  = False, 
        # window_size = (1000, 1000),
        # window_size = (1920,1080),
        method      = METHOD,
        show_axes = False,
        show_grid = False,
        scalars_norm = 'log',
        background_color = 'white',
        ##############################
        ### Widgets
        add_volume_opacity_slider = False,
        widget_type = None,
        # widget_type = 'clip slice',
        # widget_type = 'clip box',
        # widget_type = 'plane vectors',
        # widget_type = 'plane streamlines', 
        ##############################
        ### Camera
        camera_centre = None,
        focal_point = None,
        view_up = (0,0,1),
        zoom_factor = 1.5,
        ##############################
        ### Streamlines
        tube_radius         = 0.0005,
        show_source         = False,
        surface_streamlines = False,
        stream_params       = {
                'max_steps': 600, # 5000
                'max_time': 1e60,
                'terminal_speed': 1e-60,       
                'integration_direction': 'both',
                'compute_vorticity': False,
                'integrator_type': 45,   # 2 | 4 | 45
            },    
        source_radius = None,
        source_center = None,
        stream_variable_radius = 'vector magnitude',
        stream_radius_factor = 6,
        ###############################
        ### Vectors
        vector_scaling  = 'magnitude',
        vector_factor   = 1.2e-2,
        vector_max      = 2,
        ###############################
        ### MESH
        show_mesh = True,
        mesh_type = "surface",
        slice_pos = (None,None,None),
        volume_opacity = 1/15,
        culling = None,
        ###############################
        # Scalarbar GENERAL
        vertical_sbar   = True,
        cbar_width      = 0.06,
        cbar_height     = 0.8,
        cbar_title_font = 15,
        cbar_label_font = 12,
        n_colors = 256,
        annotation_color = 'black',
        ###############################
        # Scalarbar: Mesh
        scalar_sbar_title = 'scalars', # 'log[TT]\n',
        show_mesh_sbar  = False,
        mesh_opacity    = 0,  
        mesh_cmap       = 'seismic',
        scalar_sbar_pos_x = 0.91,
        scalar_sbar_pos_y = None,
        ###############################
        # Scalarbar: Field
        field_sbar_title = 'vectors',
        show_field_sbar = False,
        field_opacity   = 1,
        field_cmap      = 'jet',
        field_sbar_pos_x = 0.03,
        field_sbar_pos_y = None,
        ################################
        # Orbit gif
        orbit_gif = False,
        orbit_points = 35,
        orbit_step = 2.0,
        ################################
        # Slices
        normal = 'x',
        origin = None,
        ################################
    )

    plot = Pyvista3DPlot(
                vector_key=VECTOR_KEY, scalar_key=SCALAR_KEY, datadir=DATADIR,
                precision='f', magic=['bb', 'tt'], ivar=IVAR, coordinates=COORDINATES, 
                settings=settings, debug=True, outputdir='/home/leevi/Desktop/stream_output')     
    # plot.preview()
    # plot.vectors()
    plot.streamlines()
    # plot.scalars()
    # plot.contour(isosurfaces=10)
    # plot.movingIsovalue(isosurfaces=100)
    # plot.writeSettings('settings.json')


runPyvistaVolumePlotter()