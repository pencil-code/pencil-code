# pv_plotter.py
#
# 3D Visualization routine for creating plots from slices in cartesian | spherical |
# cylinder coordinates using PyVista.
#
# Authors:  L. Veneranta (leevi.veneranta@aalto.fi)
#
"""
Requirements
------------
These are generic requirements for all PyVista plotter tools including 
``pv_plotter.py``, ``pv_volume_plotter.py`` and ``pv_plotter_utils.py``.

* pyvista
    -> latest version (version >= 0.31.3) otherwise some streamline methods may 
        not work
* imageio-ffmpeg
* tqdm
* numpy
* sklearn --- necessary for pv_volume_plotter.py parameter tests. 

All pyvista plotting tools are tested to work with versions: VTK==9.0.3 and 
pyvista==0.31.3.

Furthermore, in order to save videos, PyVista requires ``imageio-ffmpeg``. Note that
saving only images is faster and should be used instead in combination with e.g.
imagemagick on command line. If using progress_bar parameter, ``tqdm`` is needed. 

Pyvista plot generation does not work on CSC Puhti without loading the module
``mesa-settings``.

General
-------
This file defines a routine for generating 3D visualisation from video slices. It 
is able to generate both images (e.g. png) and videos (mp4/gif). This uses PyVista 
in order to generate plots in cartesian | spherical | cylindrical coordinates.

User should call the function `plot()` to generate the visualizations. See 
documentation of plot() for further usage examples and parameters.


Plot availability for each coordinate system
--------------------------------------------

| Plot type             | Cartesian | Spherical | Cylinder |
|-----------------------|-----------|-----------|----------|
| Scalars               |     x     |     x     |     x    |
| + surface vectors     |     x     |     -     |     -    |
| + surface streamlines |     x     |     -     |     -    |


Things to note
--------------
STREAMLINES: increasing number of source points, steps in integration or even field 
itself can have a large effect on the used memory|time that takes to create one 
frame.


Settings
--------
Class ``PlotSettings()`` defines all possible settings in addition to the parameters 
function plot() takes in. For further information about all the possible parameters, 
see documentation of ``PlotSettings``. 

The class can be initialized in at least three ways: 

1. Pass in the parameters:
>>> settings = PlotSettings(off_screen=True, preview=False, ...)

2. Define a dictionary and pass that in using asterisk notation:
>>> params = {
    'off_screen': True,
    'preview': False,
    ...
}
>>> settings = PlotSettings(**params)

3. Initialize class and change the values in more 'object-oriented' way:
>>> settings = PlotSettings() 
>>> settings.off_screen = True
>>> settings.preview = False

"""
# Python extrenal libraries
import pyvista as pv
from pyvista.utilities.features import transform_vectors_sph_to_cart
import pencil as pc
import numpy as np

# Python internal libraries
from dataclasses import dataclass
from pathlib import Path
import time
import os

# Some generic utility functions are defined in .../visu/pv_plotter_utils.py
from pencil.visu.pv_plotter_utils import *


from icecream import ic

# For sanity checking input arguments. See plot() function documentation for details.
XYZPLANE_KEYS = ['xy', 'xy2', 'xz', 'yz']

# Global variable for constant random seed. This is used when PlotSettings.constant_seed 
# is set to True. The seed is chosen randomly once in the beginning. This ensures same
# randomized points.
CONSTANT_SEED = np.random.randint(0, 2**16) # maximal value is 2**32 - 1


@dataclass
class PlotSettings:
    """
    Defines all possible settings and defaults for all plots. Defining
    this as a class has the advantage of keeping all default and parameter 
    settings in one maintainable place and removes the hassle of e.g. using
    kwargs.
    
    GENERAL PARAMETERS
    ------------------
    off_screen: bool
        Whether plotting should be done on screen or off screen
    preview: bool
        Interactive preview, specific slice can be set by islice parameter in plot()
    window_size: tuple
        should be multiple of macro_block_size=16 for imageio-ffmpeg
    progress_bar: bool
        Enable tqdm progress bar 
    videoformat: str
        Options: None | 'mp4' | 'gif'
    imageformat: str
        Options: None | 'png' | 'jpeg' | etc., any Pillow compatible imageformat
    framerate: int
        Movie framerate. Does not affect gif.
    figdir: str
        Path to save images
    moviedir: str
        Path to save movies
    bg_color: str or 3 item list
        Background color. Either string or 3 tuple in RGB format.
    timestamp: bool
        Add timestamp to saved output filename.
    
    PLOT PARAMETERS
    ---------------
    norm: str
        Normalization applied. Options: 'linear' | 'log'
    coordinates: str
        Options: 'cartesian' | 'spherical' | 'cylinder'
    opacities: dict
        Opacities per slice. Should be dict with keys 'xy', 'xy2', 'yz' and 'xz'.
        #TODO! at the moment, only works for cartesian box plots
    
    COORDINATE SPECIFIC PARAMETERS
    ------------------------------
    offset: float
        Offset for the bottom slice. Z-coordinate of the bottom mesh is set to
        -zn / offset, where zn is the number of points in z-direction.
    spherical_xz_pos: float
        Options: None | -1 | <angle in radians>. Position for xz slice in spherical
        coordinates. If None, position inferred from slice data. If -1, slice set to
        close the gap between meridional slices in one end.
    
    VECTOR PLOT PARAMS
    ------------------
    n_vectors: int
        Depending on chosen vector_method:
            1. method=='every_nth': Plot only every n_vectors vectors
            2. method=='random': Plot a total of n_vectors, source point chosen 
                randomly
    vector_size: float
        Affects the size of plotted vectors.
    vector_method: str
        Options:
            1. 'random' = random sampled points
            2. 'every_nth' = plot every nth vector (defined by n_vectors)
    vector_scaling: str
        Scaling of vectors. Options: 
            1. None, no scaling applied
            2. 'magnitude', vectors scaled by their magnitude (to max size defined 
            vector_max)
            3. 'scalars', vectors scaled based on the scalar data on the surface.
    vector_max: float
        Maximum value for vector size. Useful, e.g. if vectors scaled by their magnitude
        that they're not arbitrarily large.
    surface_vectors: bool
        If True, always one vector component set to zero (depending on which mesh
        were on). If False, the plotted vectors have also their 3rd component, i.e.
        not constrained on to the 2D surface.
    
    STREAMLINE PARAMS
    -----------------
    streamlines: bool
        Enable streamlines. If False, and vectors are supplied vectors are plotted
        instead.
    stream_tube_radius: float
        Radius for the streamline tubes.
    stream_variable_radius: str
        Enable scaling of stream tube radius. Input should be name of the scalar
        array used for scaling. Available arrays include 'Magnitude' and 'scalars'.
        If None, streamline tube radius stays constant
        
        Extra: PyVista includes some internal scalar arrays if stream_params has
        'compute_vorticity'=True. This includes e.g. 'Vorticity' array. The names
        of these are the easiest found by modifying __addSurfaceStreamlines function, 
        see output of: `print(stream.point_arrays)` for available arrays.
        
    stream_radius_factor: float
        Maximal radius of stream tube as a multiple of stream_tube_radius.
    stream_src_points: int or dict
        Number of source points given either as an integer (all meshes have same
        number of src points) or as a dictionary containing keys 'xy', 'xy2', 'xz' 
        and 'yz' defining number of source points for each mesh.
    stream_show_source: bool
        Show source points as spheres in the plot.
    stream_src_radius: float
        #TODO! MISSING TOTALLY, CURRENTLY USELESS PARAMETER
    stream_log_scale: bool
        !WARNING! Not tested yet well, enables pyvista.add_mesh log_scale option
    stream_params: dict
        Any parameters pyvista streamlines_from_source() takes in OTHER THAN
        surface_streamlines, vectors. These parameters include:
            * integrator_type
            * integration_direction
            * initial_step_length
            * step_unit
            * min_step_length
            * max_step_length
            * max_steps
            * terminal_speed
            * max_error
            * max_time
            * compute_vorticity
            * interpolator_type
            * rotation_scale
        See https://docs.pyvista.org/core/filters.html?highlight=streamlines_from_source#pyvista.DataSetFilters.streamlines_from_source
        for more details on the parameters.
        
    RANDOM SAMPLING PARAMS
    ----------------------
    set_seed_1: bool
        For debugging, always sets random seed to 1, thus all sampled points are
        the same. Makes comparing streamlines easier for debugging
    constant_seed: bool
        Sets random seed once in the beginning and resets the random generator always
        with this seed. This has the consequence of randomizing initial points, 
        keeping them constant, e.g. throught the movie. Makes streamlines jump 
        less around the meshes.
        
    CAMERA PARAMS
    -------------
    camera_centre: tuple of floats
        Coordinates for camera centre in form (x, y, z).
    focal_point: tuple of floats
        Focal point for the camera in form (x, y, z).
    
    AXES PARAMS
    -----------
    show_axes: bool
        Show axes or not.
    axes_labels: tuple of str
        Labels for the axes. Should be tuple of length 3 containing strings.
    axes_font: int
        Font for the axes labels 
    
    COLORBAR PROPERTIES: field specific (vectors | streamlines)
    -----------------------------------------------------------
    show_field_sbar: bool
        Show scalarbar for vectors / streamlines
    field_cmap: str
        Matplotlib compatible colormap for vectors / streamlines
    field_sbar_title: str
        Title for fields scalarbar.
    
    pos_x and pos_y are given as a float between 0 and 1 (percentage). That is 
    it defines the distance from leftmost edge (pos_x) and bottom most edge (pos_y). 
    E.g. pos_x = 0 would mean colorbar is set at the very leftmost edge. NOTE! if
    both None, colorbars set automatically, also works for multiple colorbars set
    nicely next to each other.
    
    field_sbar_pos_x: float
        If None set automatically. Else float in range [0,1]. Note, 0.03 is pretty
        good value (fairly close to left edge)
    field_sbar_pos_y: float
        None or float between [0,1]. None works well
    
    COLORBAR PROPERTIES: scalar specific
    ------------------------------------
    show_scalar_sbar: bool
        Show scalarbar for scalars
    scalar_cmap: str
        Matplotlib compatible colormap for scalars.
    scalar_sbar_title: str
        Title for scalar scalarbar. 
    scalar_sbar_pos_x: float    
        None or float between [0,1]. Note, 0.88 is pretty good value (enough 
        close to right edge). See notes above for detailed explanation of pos_x
    scalar_sbar_pos_y: float
        None or float between [0,1]. Note, None works well. See notes above for 
        detailed explanation of pos_y
    
    COLORBAR PROPERTIES: Generic properties
    ---------------------------------------
    Following parameters apply both to field and scalar colorbars.
    
    cbar_width: float
        Width as a percentage, between 0 and 1
    cbar_height: float
        Height as a percentage, between 0 and 1
    cbar_title_font: int
        Title font size for scalarbar
    cbar_label_font: int
        Label font size for scalarbar
    n_colors: int
        Number of colors for the scalarbar.
    _sbar_args: dict
        INTERNAL VALUE. SHOULD NOT BE USED.
    
    TITLE ANNOTATIONS
    -----------------    
    title_position: str
        Options: 'lower_left', 'lower_right', 'upper_left', 'upper_right', 
        'lower_edge', 'upper_edge', 'right_edge', and 'left_edge'. 
    title_font: int
        Font size for the title.
    str_unit: str
        Unit added behind the timestamp in the title. I.e. title is form 
        '<current time step><str_unit>'
    tscale: float
        Multiplicative scaling for the timestamp in title.
    time_precision: int
        Number of decimals shown for the timestamp.
    """
    ### General settings #######################################################
    off_screen: bool   = True
    preview: bool      = False
    window_size: tuple = (1024, 768)  
    progress_bar: bool = True
    videoformat: str   = None
    imageformat: str   = 'png'
    framerate: int     = 15
    figdir: str        = './images/'
    moviedir: str      = './movies/'
    bg_color: str      = 'white'
    timestamp: bool    = True
    
    ### Plot parameters ########################################################
    norm: str          = 'linear'
    coordinates: str   = 'cartesian'
    opacities: dict    = None
    
    ### Coordinate specific parameters #########################################
    offset: float           = 2.0  
    spherical_xz_pos: float = -1
    
    ### Vector plot params #####################################################
    n_vectors: int          = 500
    vector_size: float      = 7
    vector_method: str      = 'random'
    vector_scaling: str     = 'magnitude'
    vector_max: float       = None
    surface_vectors: bool   = True    
          
    ### Streamline parameters ##################################################
    streamlines: bool           = False
    stream_tube_radius: float   = 0.05
    stream_variable_radius: str = 'Magnitude'
    stream_radius_factor: float = 10
    stream_src_points: int      = 50
    stream_show_source: bool    = False
    stream_src_radius: float    = None  #TODO! MISSING TOTALLY, CURRENTLY USELESS PARAMETER
    stream_log_scale: bool      = False # !WARNING! not tested yet too well, enables add_mesh log_scale
    stream_params: dict         = None
    
    ### Random sampling params #################################################
    set_seed_1: bool    = False
    constant_seed: bool = False
    
    ### Camera params ##########################################################
    camera_centre: tuple = None
    focal_point: tuple   = None
    
    ### Axes params ############################################################
    show_axes: bool    = True
    axes_labels: tuple = ('x', 'y', 'z')
    axes_font: int     = 16
    
    ### Colorbar properties ####################################################
    ##### Field (vectors / streamlines) colorbar specific properties ###########
    show_field_sbar: bool   = True
    field_cmap: str         = 'bwr' 
    field_sbar_title: str   = ''
    field_sbar_pos_x: float = None
    field_sbar_pos_y: float = None
    
    ##### Scalar colorbar specific properties ##################################
    show_scalar_sbar: bool   = True
    scalar_cmap: str         = 'bwr'
    scalar_sbar_title: str   = ''
    scalar_sbar_pos_x: float = None
    scalar_sbar_pos_y: float = None
    
    ##### General colorbar properties ##########################################
    cbar_width: float    = 0.06
    cbar_height: float   = 0.65
    cbar_title_font: int = 10
    cbar_label_font: int = 8
    n_colors: int        = 256
    _sbar_args: dict     = None     # !Internally used value, should not be used!
    
    ### Title annotations ######################################################
    title_position: str = 'upper_left'
    title_font: int     = 20
    str_unit: str       = ''
    tscale: float       = 1
    time_precision: int = 3
    
    ################################################################################################

    def __post_init__(self):
        """
        Initializes Paths outside the __init__ function automatically created by
        the dataclass. Especially, this converts the string paths to Path objects.
        """
        self.figdir   = Path(self.figdir) 
        self.moviedir = Path(self.moviedir)

        # self.stream_src_points = {'xy': 50, 'xy2': 50, 'xz': 100, 'yz': 100}
        
        if self.opacities == None:
            self.opacities = {'xy': 1.0, 'xy2': 1.0, 'xz': 1.0, 'yz': 1.0}

        if self.stream_params == None:
            # Set defaults for streamline parameters, can't be set above since
            # only immutable types can be as defaults
            self.stream_params = {
                'max_steps': 1000,
                'max_time': 1e60,
                'terminal_speed': 1e-60,
                'integration_direction': 'both',
                'compute_vorticity': False,
                'integrator_type': 45,
            }
            
        self._sbar_args = {
            'vertical': True,
            'width': self.cbar_width, 
            'height': self.cbar_height, 
            'color': 'black',
            'title_font_size': self.cbar_title_font,
            'label_font_size': self.cbar_label_font,
            'n_colors': self.n_colors
            }


################################################################################
############ INTERNAL FUNCTIONS ################################################
################################################################################

################################################################
### HELPER FUNCTIONS ###########################################
################################################################


def __addVectorsToMesh(
    # PyVista visualization specific
    plotter: pv.Plotter, mesh, key, settings):
    """
    Adds vectors (glyph) to mesh. This should work for any kind of a mesh, e.g. 
    capable of creating vectors on a 2D surface, but also to 3D volumes. 
    
    NOTE! Vectors should already be added to the mesh as an array named 'vectors'
    before calling this function!
    
    Parameters
    ----------
    plotter : pv.Plotter
        Instance of pv.Plotter
    meshes : pyvista mesh
        Mesh to add vectors to.
    key : str
        Name of the mesh, used to name the added vectors on to the plotter. If 
        they'd have same names some would be overwritten.
    settings : PlotSettings, optional
        PlotSettings object
    
    Last updated: 24.6.2021
    """            
    assert settings.n_vectors > 0, "Parameter n_vectors should be greater than 0."    
    vecs = mesh['vectors']
    magnitude = np.sqrt(vecs[:,0]**2 + vecs[:,1]**2 + vecs[:,2]**2)
    mesh['Magnitude'] = magnitude
    
    if settings.vector_scaling == 'magnitude':
        # If vector max not set, use mean + one STD of magnitudes as a reasonable value
        if settings.vector_max is None:
            settings.vector_max = magnitude.mean() + magnitude.std()
        mesh['scale'] = np.array([settings.vector_max if a > settings.vector_max 
                                  else a for a in magnitude])
        scale = 'scale'
    elif settings.vector_scaling == 'scalars':
        scale = 'scalars'
    elif settings.vector_scaling is None:
        scale = None
    else:
        raise ValueError(f'[__add_vectors] Unknown scaling method scaling = {settings.vector_scaling}')
    
    if settings.vector_method == 'every_nth':
        mask = np.zeros_like(vecs)
        mask[::settings.n_vectors,:] = 1
        mesh['mask'] = mask 
        glyph = mesh.glyph(orient='vectors', scale=scale, geom=pv.Arrow(),
                            factor=settings.vector_size).threshold((0.1,1.1), scalars='mask')
    elif settings.vector_method == 'random':
        sampled = randomSampleMeshPoints(settings.n_vectors, mesh, 
                                           set_seed_1=settings.set_seed_1, 
                                           constant_seed=settings.constant_seed,
                                           get_arrays=True,
                                           seed=CONSTANT_SEED)
        glyph = sampled.glyph(orient='vectors', scale=scale, geom=pv.Arrow(),
                              factor=settings.vector_size)
    else:
        raise ValueError(f'[__add_vectors] Unknown sampling method = {settings.vector_method}.' \
            'should be either "every_nth" or "random".')
    
    mesh.set_active_vectors('vectors')
    mesh.set_active_scalars('scalars')
    
    plotter.add_mesh(glyph, name=f'{key} vectors', 
                    show_scalar_bar=settings.show_field_sbar,
                    cmap=settings.field_cmap,
                    scalar_bar_args=create_sbar_args(
                                        settings,
                                        settings.field_sbar_title,
                                        settings.field_sbar_pos_x,
                                        settings.field_sbar_pos_y))


def __addSurfaceStreamlines(plotter: pv.Plotter, mesh: pv.StructuredGrid, key, settings):
    """
    Adds streamlines onto a mesh.
    """    
    vecs = mesh['vectors']
    magnitude = np.sqrt(vecs[:,0]**2 + vecs[:,1]**2 + vecs[:,2]**2)
    mesh['Magnitude'] = magnitude
    
    if isinstance(settings.stream_src_points, dict):
        assert key in settings.stream_src_points.keys(), f'settings.stream_src_points does not contain key {key}' 
        src_points = settings.stream_src_points[key]
    else:
        src_points = settings.stream_src_points
    
    src = randomSampleMeshPoints(src_points, mesh=mesh, 
                                   set_seed_1=settings.set_seed_1,
                                   constant_seed=settings.constant_seed,
                                   seed=CONSTANT_SEED)
    
    stream = mesh.streamlines_from_source(src, surface_streamlines=True, 
                            vectors='vectors', **settings.stream_params)
    
    tube = stream.tube(radius=settings.stream_tube_radius, 
                       scalars=settings.stream_variable_radius,
                       radius_factor=settings.stream_radius_factor,)
    
    if settings.stream_show_source:
        plotter.add_mesh(src, render_points_as_spheres=True, 
                         show_scalar_bar=False)    
    
    if tube.n_points == 0:  
        print(f'WARNING! stream.tube has zero points! Skipping this mesh!')
    else:        
        # scalarbar will use the active scalars to show values
        tube.set_active_scalars('Magnitude')
        
        plotter.add_mesh(tube, name=f'{key} streamlines', 
                        cmap=settings.field_cmap,
                        clim=[magnitude.min(), magnitude.max()],
                        log_scale=settings.stream_log_scale,
                        show_scalar_bar=settings.show_field_sbar,
                        scalar_bar_args=create_sbar_args(settings,
                                                         settings.field_sbar_title,
                                                         settings.field_sbar_pos_x,
                                                         settings.field_sbar_pos_y))


def __normScalars(scalars, cmin, cmax, field, norm='linear', dtype='f') -> tuple:
    """
    Applies required normalization on the 2D scalar slices.
    Parameters
    ----------
    scalars : dictionary
        Dictionary for box slices that has contains the values for 
        the 2D slices.
    cmin : int
        Colorbar minimum
    cmax : int 
        Colorbar maximum
    field : string
        Name of the current field.
    norm : string, optinonal
        String key for the normalization applied on the scalars. Options:
        'log' or 'linear'. Default: 'linear'
    dtype : string, optional
        Data type. Default: 'f'.
    
    Returns
    -------
    scalars : dictionary
        Scalars with the normalized slices
    cmin : int
        Colorbar minimum
    cmax : int
        Colorbar maximum
    Notes
    -----
    Last updated 31.5.2021 9:05
    """
    # Possible data normalizations
    # TODO! type conversions redundant?
    def ln(x): return np.log10(np.exp(np.float32(x))).astype(dtype)
    def fun(x): return np.log10(x).astype(dtype)
        
    for key in scalars.keys():        
        if norm == 'log':
            if 'ln' in field:
                scalars[key] = ln(scalars[key])
            else:
                scalars[key] = fun(scalars[key])
            
    # Apply necessary normalizations for the colorbar [cmin,cmax]
    if norm == 'log':
        if 'ln' in field:
            cmax = ln(cmax)
            cmin = ln(cmin)
        else:
            cmax = np.log10(cmax)
            cmin = np.log10(cmin)
    elif norm == 'linear':
        cmax = max(-cmin, cmax)
        cmin = -cmax
    else:
        print("WARNING: 'norm' undefined, applying 'linear'")
        cmax = max(-cmin, cmax)
        cmin = -cmax

    return scalars, cmin, cmax


def __vectorsFromSlices(slice_obj, fields, datadir, vectors_unit_length=False, 
                        coordinates='cartesian', surface_vectors=True):
    """
    Function to generate vectors from slice data. Fields defines the vector components,
    i.e. fields[0] defines what field is the first vector component, fields[1] the 
    second etc.
    Parameters
    ----------
    slice_obj : pencil.read.allslices.SliceSeries
        Pencil slice data.
    fields : list of length 3
        Defines the fields used for vector components. Should be of length 3.
    vectors_unit_length : bool, optional
        If True, normalizes the vectors to unit length.
    coordinates : string, optional
        Coordinate system used
        
    Returns
    -------
    vectors : dict
        Dictionary of vectors. Vectors are of shape [time, points, 3]
    """
    vectors = {}
    assert len(fields) == 3, "Fields should have exactly 3 values"
    
    def normalize(v):
        shape = v.shape
        norm = np.linalg.norm(v, axis=-1) # Axis 1 since ==> [time, vectors, components]
        return v / norm.reshape(shape[0], shape[1], 1) 
    
    time_shape = slice_obj.t.shape[0]
    
    for slice in ['xy', 'xy2', 'xz', 'yz']:
        if coordinates == 'cartesian':
            v1 = slice_obj.__getattribute__(slice).__getattribute__(
                                    fields[0]).swapaxes(1,2).reshape(time_shape, -1)
            v2 = slice_obj.__getattribute__(slice).__getattribute__(
                                    fields[1]).swapaxes(1,2).reshape(time_shape, -1)
            v3 = slice_obj.__getattribute__(slice).__getattribute__(
                                    fields[2]).swapaxes(1,2).reshape(time_shape, -1)

        elif coordinates == 'spherical':
            v1 = slice_obj.__getattribute__(slice).__getattribute__(fields[0])
            v2 = slice_obj.__getattribute__(slice).__getattribute__(fields[1])
            v3 = slice_obj.__getattribute__(slice).__getattribute__(fields[2])
            
            var = pc.read.var(datadir=datadir, trimall=True, precision='f')
            r, theta, phi = var.x, (var.y*180)/np.pi, (var.z/180)/np.pi
            nr, ntheta, nphi = r.shape[0], theta.shape[0], phi.shape[0]

            # NOTE! vectors now have one extra dimension (time) compared to var 
            # plots. The transfrom fucntion would either need to be run once per
            # time instant or something else (calculate one conversion matrix 
            # and just use that + matrix product). Check the source of pyvista for
            # the conversion --> make it into a matrix?
            shape = (time_shape, nr, ntheta, nphi)
            u = np.zeros(shape)
            v = np.zeros(shape)
            w = np.zeros(shape)
            
            ic(v1.shape, v2.shape, v3.shape)
            ic(shape)

            if slice == 'xy' or slice =='xy2':
                u[:,]
            elif slice == 'xz':
                pass
            elif slice == 'yz':
                pass
            

            # The transformed vectors 
            v1_t, v2_t, v3_t = np.zeros_like(v1), np.zeros_like(v2), np.zeros_like(v3)
            
            # Stupid solution: loop over all time instants converting each time instant individually.
            # Better solution would be to just calculate the matrix once defined in pv.transform_vectors_sph_to_cart
            #   and apply it straight to the v1,v2,v3 without the looping and recalculations. 
            for t in range(time_shape):
                v1_t[t], v2_t[t], v3_t[t] = pv.transform_vectors_sph_to_cart(phi, theta, r, v1[t], v2[t], v3[t])

        elif coordinates == 'cylinder':
            raise NotImplementedError("Currently __vectorsFromSlices() cannot convert vectors from"
                " cylinder coordinates to cartesian which is for the vectors to be added to the"
                " Pyvista plot!")
        
        if surface_vectors:
            # Set correct vector component to zero
            if slice == 'xy' or slice == 'xy2':
                v3 = np.zeros_like(v3)
            elif slice == 'xz':
                v2 = np.zeros_like(v2)
            elif slice == 'yz':
                v1 = np.zeros_like(v1)
        
        n_vectors  = v1.shape[1]
        vecs = np.empty((time_shape, n_vectors, 3))
        vecs[:,:,0] = v1
        vecs[:,:,1] = v2
        vecs[:,:,2] = v3
        
        if vectors_unit_length:
            vectors[slice] = normalize(vecs)
        else:
            vectors[slice] = vecs
        
    return vectors


def __updateVectorData(mesh, vectors, t=0):
    """
    Updates the vector data on a mesh. Does not generate the vector visualization
    itself but initializes / updates the vector field.
    """        
    mesh['vectors'] = vectors[t]
    mesh.set_active_vectors('vectors')
    mesh.set_active_scalars('scalars')


def __getFilename(type, field, settings, idx=-1, timestamp=False):
    """
    Returns a filename given type. Keeps the filenaming logic in a single place.
    
    Parameters 
    ----------
    type: str
        TODO!
    field: str
        
    settings: PlotSettings
    
    idx: int, optional
    timestamp: bool, optional
    
    
    Returns
    -------
        path: pathlib.Path
            Path object pointing to the correct directory and output file of given 
            type.
    """
    if type in ['mp4', 'gif']:
        if timestamp:
            name = f'{field}_movie_{int(time.time())}.{type}'
        else:
            name = f'{field}_movie.{type}'
        path = settings.moviedir / name
    # Filetype is an image
    else:
        if timestamp: 
            name = f'{field}_{idx:04d}_{int(time.time())}.{type}'
        else:
            name = f'{field}_{idx:04d}.{type}'
        path = settings.figdir / name
    
    return path


def __plot_annotations(p: pv.Plotter, s: PlotSettings, slice_obj, itt) -> None:
    """
    Applies annotations to the pyvista.Plotter instance. This includes
    scalarbar, text, axes, etc.
    Parameters
    ----------
    p : pv.Plotter
        Instance of pyvista.Plotter to apply annotations on.
    s : PlotSettings
        PlotSettings object
        
    Notes
    -----
    Last updated: 13.7.2021
    """
    __updateTitleText(p, s, slice_obj, itt)
    
    # Set camera position, sets [centre, focalpoint, viewup]
    if s.camera_centre is not None and s.focal_point is not None:
        p.camera_position = (s.camera_centre, s.focal_point, (0, 0, 1))
    elif s.camera_centre is not None and s.focal_point is None:
        # Now camera is defined but use the focal point given by PyVista
        _, focal, _ = p.camera_position
        p.camera_position = (s.camera_centre, focal, (0,0,1))
        
    # Set other properties of the plot
    if s.show_axes:
        p.show_bounds(xlabel=s.axes_labels[0], ylabel=s.axes_labels[1], zlabel=s.axes_labels[2],
                      color='black', location='outer', font_size=s.axes_font)
    p.background_color = s.bg_color


def __updateTitleText(plotter: pv.Plotter, s: PlotSettings, slice_obj, itt) -> None:
    """
    Updates plot title text. Title is of form <TIME><UNIT> where TIME is the current
    time instant rounded to wanted precision and possibly scaled. UNIT is a unit
    (a string) for the time.
    """
    # First part controls timestamp, scaled by tscale and rounded to time_precision
    # decimal places. NOTE! .{}f notation makes sure there are enough zeros 
    # at the end if it rounds down to fewer decimal places. ==> otherwise the 
    # time unit would start jumping left and right if the number of printed decimals
    # change all the time.
    text = "{:.{}f}".format(slice_obj.t[itt] * s.tscale, s.time_precision) + f' {s.str_unit}'
    
    plotter.add_text(
        text=text,
        position=s.title_position, 
        font_size=s.title_font, 
        color='black', 
        # name needs to be title so that previous title actor is overwritten
        name='title'    
        )


##################################################################
### PLOTTING FUNCTIONS ###########################################
##################################################################


def __pv_cylinder_plot(p: pv.Plotter, slice_obj, scalars, t, lims, settings,
                        dtype='f', datadir='./data',):
    """
    Generate pyvista.StructuredGrid straight from cylinder coordinates.
    
    Parameters
    ----------
    slice_obj : pencil.read.allslices.SliceSeries
        Pencil slice data
    scalars : dict
        Dictionary of numpy.ndarrays
    t : int
        Index for the time instant
    lims : list, optional
        Limits for the colorbar, by default [-1,1]
    settings : PlotSettings
        PlotSettings object
    datadir : str, optional
        Data directory, by default './data'
    
    Returns
    -------
    dict
        Dictionary of meshes on the plotter that can be used to update the scalars.
    
    Notes
    -----
    Last updated: 18.6.2021
    """    
    
    grid = pc.read.grid(datadir=datadir, trim=True, precision=dtype)
    r, theta, z = grid.x, grid.y, grid.z
    
    # Create arrays of grid cell boundaries, which have shape of (x.shape[0] + 1)
    r_bounds   = cellBounds(r)       # Radial distance
    t_bounds   = cellBounds(theta)     # Polar angle theta 
    z_bounds   = cellBounds(z)       # Height z
    
    tz_slice_pos = slice_obj.position.yz[0]
    rz_slice_pos = slice_obj.position.xz[0]
    ###########################################################################
    # MESHES
    meshes = {
        'xy2' : gridFromCylCoords(r_bounds, t_bounds, z[0]),
        'xy'  : gridFromCylCoords(r_bounds, t_bounds, z[-1]),
        'yz'  : gridFromCylCoords(tz_slice_pos, t_bounds, z_bounds), # tz_slice_pos = radial position of slice
        'xz'  : gridFromCylCoords(r_bounds, rz_slice_pos, z_bounds), # rz_slice_pos = theta position of slice
    }

    for key in scalars.keys():
        meshes[key]['scalars'] = scalars[key][t].reshape(-1)
        p.add_mesh(meshes[key], name=f'{key}', #opacity=opacities[key] #TODO! Opacities
                   show_scalar_bar=settings.show_scalar_sbar,
                   scalar_bar_args=create_sbar_args(settings,
                                                      settings.scalar_sbar_title,
                                                      settings.scalar_sbar_pos_x,
                                                      settings.scalar_sbar_pos_y),
                   clim=lims,
                   cmap=settings.scalar_cmap)

    return meshes


def __pv_box_plot(p: pv.Plotter, scalars, t, lims, settings,
                  opacities={'xy': 1.0, 'xy2': 1.0, 'xz': 1.0, 'yz': 1.0}):
    """
    Standalone function that wraps PyVista calls to create a box plot consisting of 4 2D slices with
    the xy slice offsetted below the cube.
    
    Parameters
    -----------
    scalars : dict
        Dictionary of 2D scalar slices
    t : int
        Index of time instant to plot
    lims : list, optional
        Limits for the colorbar, by default [0, 1]
    settings : PlotSettings
        PlotSettings object
    
    Returns
    --------
    meshes : dict
        Dictionary with keys ['z2', 'z1', 'xz', 'yz'] and as values
        PyVista meshes for the corresponding sides. Needed for updating the scalars
        on the meshes.
    
    Notes
    -----
    Last updated: 18.6.2021
    """
    xn = scalars['xy2'][0].shape[0]     # Order is [x,y]
    yn = scalars['xy2'][0].shape[1]     # Order is [x,y]
    zn = scalars['yz'][0].shape[1]  # order is [y,z]
    # Coordinates
    xx = np.linspace(0, xn - 1, xn)
    yy = np.linspace(0, yn - 1, yn)
    zz = np.linspace(0, zn - 1, zn)
    x_z, z_x = np.meshgrid(xx, zz)  # Coordinates for xz plane
    y_z, z_y = np.meshgrid(yy, zz)  # Coordinates for yz plane
    x_y, y_x = np.meshgrid(xx, yy)  # Coordinates for xy plane

    # Offsets for z2, z1 and side slices.
    z1_offset = (-zn/settings.offset) * np.ones_like(x_y)
    z2_offset = (zn-1) * np.ones_like(x_y)
    y1_offset = np.zeros_like(z_x)
    x1_offset = np.zeros_like(z_y)
    
    # Generate meshes
    meshes = {
            'xy'  : pv.StructuredGrid(x_y, y_x, z1_offset),
            'xy2' : pv.StructuredGrid(x_y, y_x, z2_offset),
            'xz'  : pv.StructuredGrid(x_z, y1_offset, z_x),
            'yz'  : pv.StructuredGrid(x1_offset, y_z, z_y),
            }
    
    
    
    # pyvista only adds scalarbar if the title does not exist. If title is not give,
    # it is automatically the point array's name
    for key in scalars.keys():        
        meshes[key]['scalars'] = scalars[key][t].reshape(-1)
        p.add_mesh(meshes[key], name=f'{key}', opacity=opacities[key], 
                   show_scalar_bar=settings.show_scalar_sbar,
                   clim=lims,
                   cmap=settings.scalar_cmap,
                   scalar_bar_args=create_sbar_args(settings,
                                                      settings.scalar_sbar_title,
                                                      settings.scalar_sbar_pos_x,
                                                      settings.scalar_sbar_pos_y))
    return meshes


def __pv_spherical_plot(p: pv.Plotter, slice_obj, scalars, t, lims, settings,
                            dtype='f', datadir='./data'):  
    """
    Initialize Pyvista plot for spherical coordinates.
    
    Parameters
    ----------
    slice_obj : pencil.read.allslices.SliceSeries
        Pencil slice data
    scalars : dict
        Scalars in dictionary
    t : int
        Index of time instant
    dtype : str, optional
        Datatype, by default 'f'
    datadir : str, optional
        Data directoyr, by default './data'
    opacities : list, optional
        Opacities for the meshes, by default [1.0, 1.0, 0.7, 0.7]
    cmap : str, optional
        Matplotlib compatible colormap, by default 'bwr'
    lims : list, optional
        Colorbar limits, by default [-1,1]
    spherical_xz_pos : None | -1 | <angle in radians>
        Theta angle for xz slice. If None, angle inferred from slice data. If -1, set to
        "close" the meridional slices. Otherwise supplied angle is used for the slice positioning
        (in radians).
    off_screen : bool, optional
        Whether plotting should be done on screen or off screen, by default True
    window_size : list, optional
        Window size, by default [1024,768]
    
    Notes
    -----
    - Pencil code gives radians, PyVista uses degrees
    - Note that data is assumed to be in spherical coordinates
    - i.e. (x,y,z) --> (r, theta, phi) where
        - r radial distance
        - theta polar angle / latitude, [0,180)
        - phi azimuthal / longitudinal angle [0,360)
        - phi2 = half of full sphere i.e. azimuthal angle from 0-180
        
    Last updated: 18.6.2021
    """
    grid = pc.read.grid(datadir=datadir, trim=True, precision=dtype)
    r, theta, phi = grid.x, grid.y, grid.z
    
    # Dimensions
    nphi = phi.shape[0]

    # ntiles will tell you how many times you need to replicate the data to
    # cover a full sphere
    ntiles = int(2*np.pi / phi[nphi-1])

    # number of wedge tiles left open in the largest onion layer, between [1, ntiles)
    ntileso = ntiles-1

    # full longitudinal range i.e. [0,360), note original input values for phi
    # overwritten here, we don't need them
    phi = 360 * np.arange(ntiles*nphi) / (ntiles*nphi - 1)

    # Removing ntileso tiles to make the next layer of onion peeling
    phi2 = ((ntiles-ntileso)/ntiles * 360.0) * np.arange((ntiles-ntileso)* nphi) / ((ntiles-ntileso) * nphi-1)
    
    theta = theta * 180.0 / np.pi       # Polar angle [0,180), i.e. colatitude
    # TODO! TEST if the angles go wrong
    theta = np.abs(180 - theta)
    
    # Use wanted theta value and convert from radians to degrees
    if settings.spherical_xz_pos is None:
        xz_slice_pos = slice_obj.position.xz[0] * 180 / np.pi
    elif settings.spherical_xz_pos == -1:
        xz_slice_pos = theta[-1]     # Note, theta is in degrees already
    else:
        # Now spherical_xz_pos should be angle in radians
        xz_slice_pos = settings.spherical_xz_pos * 180 / np.pi
        
    # Create arrays of grid cell boundaries, which have shape of (x.shape[0] + 1)
    # phi_bounds = _cell_bounds(phi)        # Azimuthal angle [0,360) i.e. longitude
    # phi_bounds2 = cellBounds(phi2)     # Azimuthal angle (half of full sphere)
    # theta_bounds = cellBounds(theta)   # Latitude [0,180) i.e. polar angle
    # r_bounds = cellBounds(r)           # Radial distance
    
    phi_bounds2 = phi2
    theta_bounds = theta
    r_bounds = r

    # Phi positions for the meridional slices
    mer1_phi = 0.
    mer2_phi = -ntileso * 360 / ntiles 
        
    ###########################################################################
    # MESHES
    meshes = {
        'yz'  : pv.grid_from_sph_coords(phi_bounds2, theta_bounds, r[-1]),
        'xz'  : pv.grid_from_sph_coords(phi_bounds2, xz_slice_pos, r_bounds),
        'xy'  : pv.grid_from_sph_coords(mer1_phi, theta_bounds, r_bounds),
        'xy2' : pv.grid_from_sph_coords(mer2_phi, theta_bounds, r_bounds),
    } 
    
    for key in scalars.keys():
        # TODO! FIX THIS ==> hack to make the spherical surface data be correct,
        #   otherwise the "bananas" on the surface goes wrong way
        if key == 'yz':
            meshes[key]['scalars'] = scalars[key][t].T.reshape(-1)
        else:
            meshes[key]['scalars'] = scalars[key][t].reshape(-1)
            
        p.add_mesh(meshes[key], name=f'{key}', # opacity=opacities[key] #!TODO! Opacities
                   show_scalar_bar=settings.show_scalar_sbar,
                   scalar_bar_args=create_sbar_args(settings,
                                                      settings.scalar_sbar_title,
                                                      settings.scalar_sbar_pos_x,
                                                      settings.scalar_sbar_pos_y),
                   clim=lims,
                   cmap=settings.scalar_cmap)
    return meshes


def __plot_field(
    slice_obj, scalars, field, it, vectors=None,
    datadir='./data', 
    cmin=None, cmax=None, dtype='f',
    settings=PlotSettings(), debug=False,
    ) -> None:
    """
    Internal function to plot a field at all given time instants.
    Internal function that handles plotting of one field at all given time
    instants. Applies necessary transformations to the 2D slices defined by 
    parameter norm. Creates a visualization based on which coordinates
    one is using.
    Function can save both images and videos. Note, that saving videos requires
    the Python library `imageio-ffmpeg` since PyVista uses that internally.
    Parameters
    ----------
    slice_obj : pencil.read.allslices.SliceSeries
        Pencil slice data.
    scalars : dict
        Dictionary of scalar data. Keys should be compatible with what the plotting functions
        use. Currently allowed keys defined by XYZPLANE_KEYS.
    field : str
        Field that is being plotted.
    it : numpy.ndarray
        List of times to be plotted.
    vectors : numpy.ndarray
        Numpy array of vectors, shape [mesh.n_points, 3]. By default None.
    datadir : str, optional
        Data directory, by default './data'
    cmin : float, optional
        Colorbar minimum. If none, set to min of all slices. By default None.
    cmax : float, optional
        Colorbar maximum. If none, set to max of all slices. By default None.
    dtype : str, optional
        Datatype. Pencil reading routine options are: 'half', 'f' or 'd'. By default 'f'
    settings: PlotSettings, optional
        Dataclass containing plot settings.
    debug : bool, optional
        Enable debugging prints, by default False
        Parameters passed on to plotting functions.
        
    Notes
    -----
    Last updated: 21.6.2021
    """
    # Create figdir and moviedir if it does not exist
    if not settings.figdir.exists() and settings.imageformat is not None:
        os.mkdir(settings.figdir)
    if not settings.moviedir.exists() and settings.videoformat is not None:
        os.mkdir(settings.moviedir)
    # If preview enabled, plotting should be done ON-screen not off
    if settings.preview:
        settings.off_screen = False
    else:
        print(f'Warning! preview=False but off_screen=False. Plotting should be done'
        ' off screen currently. Setting off_screen=True')
        settings.off_screen = True

    # If [cmin,cmax] not given, set automatically to min and max of all 2D slices
    if cmin is None and cmax is None:
        cmin, cmax = 1e38, -1e38
        for key in scalars.keys():       # Loop over ['z2', 'z1', ...]
            cmax = max(cmax, scalars[key].max())
            cmin = min(cmin, scalars[key].min())

    # Apply normalization on scalars and clims
    
    scalars, cmin, cmax = __normScalars(scalars, cmin, cmax, field, norm=settings.norm,
                                        dtype=dtype)

    ########################################################################
    
    # Switch fortran order to Python i.e. [t,y,x] ===> [t,x,y]
    for key in scalars.keys():
        scalars[key] = scalars[key].swapaxes(1,2) # 0th axis is time, switch 1-2 axes around
    
    ########################################################################
    
    # try:
    #     from pyvistaqt import BackgroundPlotter
    #     p = BackgroundPlotter(off_screen=settings.off_screen, 
    #                 window_size=settings.window_size)
    # except Exception as e:
    #     print(f'>>> WARNING! Exception when creating background plotter: {e}')
    #     print(f'---> Creating just a pyvista.Plotter instead!')
    #     p = pv.Plotter(off_screen=settings.off_screen, window_size=settings.window_size)

    p = pv.Plotter(off_screen=settings.off_screen, window_size=settings.window_size)

    # Initialize plotter with the very first time instant i.e. it[0]
    if settings.coordinates == 'cartesian':
        meshes = __pv_box_plot(p, scalars, it[0], [cmin, cmax], settings,
                                  opacities=settings.opacities)
    elif settings.coordinates == 'spherical':
        meshes = __pv_spherical_plot(p, slice_obj, scalars, t=it[0], lims=[cmin,cmax], 
                                    dtype=dtype, datadir=datadir, settings=settings)
    elif settings.coordinates == 'cylinder':
        meshes = __pv_cylinder_plot(p, slice_obj, scalars, t=it[0], lims=[cmin,cmax], 
                                    dtype=dtype, datadir=datadir, settings=settings)
    
    # If vectors ==> create either vectors on meshes or streamlines   
    if vectors is not None:
        for key in vectors.keys():
            # Add vector field data on to the mesh first before plotting
            __updateVectorData(meshes[key], vectors[key], t=0)
            
            # Add vectors visualization to the plotter
            if not settings.streamlines:
                __addVectorsToMesh(plotter=p, mesh=meshes[key], key=key,
                                    settings=settings)
            
            # Add streamline visualization to the plotter
            else:
                __addSurfaceStreamlines(p, meshes[key], key, settings=settings)
            
    # Apply annotations on the plot
    __plot_annotations(p, settings, slice_obj, it[0])
    
    # !TODO! is this an issue or not
    # Set to false, otherwise later in time loop, PyVista recreates scalar bars at
    # each timestep consuming more and more memory ==> crashes
    settings.show_field_sbar = False
    settings.show_scalar_bar = False

    # If preview enabled = show preview window i.e. plot at the first time instant
    # with camera parameters added as title text to the screen
    if settings.preview:
        print(f'Starting preview...')
        plotPreview(p)
        return
    # Open a mp4 file and save first frame
    elif settings.videoformat == 'mp4': 
        path = __getFilename(settings.videoformat, field, settings, 
                             timestamp=settings.timestamp)
        p.open_movie(path, framerate=settings.framerate)
        p.write_frame()
    # Open a gif file and save first frame
    elif settings.videoformat == 'gif':
        path = __getFilename(settings.videoformat, field, settings, 
                             timestamp=settings.timestamp)
        p.open_gif(path)
        p.write_frame()
    
    # Save image of the first time instant    
    if settings.imageformat is not None:
        path = __getFilename(settings.imageformat, field, settings, idx=0,
                             timestamp=settings.timestamp)
        p.screenshot(path)
    
    # Start a TQDM progress bar with maximum of len(it) iterations
    if settings.progress_bar:
        # Create a tqdm progress bar to track the iterations
        from tqdm import tqdm
        pbar = tqdm(total=len(it), desc=f'Field: {field}', unit='frame')
        pbar.update(1)       
    
    # Loop over all except first time instant (done already in initialization)
    for itt in it[1:]:        
        # Update vector / streamline data
        if vectors is not None:              
            for key in vectors.keys():
                # Update the vector values on the mesh
                __updateVectorData(meshes[key], vectors[key], t=itt)
                
                # Add vector visualization to mesh
                if not settings.streamlines:
                    __addVectorsToMesh(plotter=p, mesh=meshes[key], key=key,
                                        settings=settings)
                # Add streamline visualization to mesh
                else:
                    __addSurfaceStreamlines(p, meshes[key], key, settings=settings)
                    
        # Update 2D scalar slices on the box        
        for key in scalars.keys():                    
            meshes[key].set_active_scalars('scalars')
            if key == 'yz' and settings.coordinates == 'spherical':
                # TODO! FIX THIS ==> hack to make the spherical surface data be
                #! correct, otherwise the "bananas" on the surface goes wrong way
                p.update_scalars(scalars[key][itt].T.ravel(), mesh=meshes[key]) 
            else:
                p.update_scalars(scalars[key][itt].ravel(), mesh=meshes[key]) 
                
        # Note this rewrites the previous text actor since it has the same name "title"
        __updateTitleText(p, settings, slice_obj, itt)

        if settings.videoformat is not None:
            p.write_frame()  
        if settings.imageformat is not None:
            path = __getFilename(settings.imageformat, field, settings, idx=itt,
                                 timestamp=settings.timestamp)
            p.screenshot(path)
        if settings.progress_bar:
            pbar.update(1)
        else:
            print(f'> [{itt+1} / {len(it)}] frame done.')

    # Movie done!
    p.close()
    if settings.progress_bar:
        pbar.close()
    if debug:
        print(f'Done rendering images/video for field: {field}.')


################################################################################
############ USER CALLED FUNCTIONS #############################################
################################################################################


def plot(
        slice_obj = None,
        datadir='./data', precision='f',
        # acquire the slice objects
        fields=['uu1', ],
        xyzplane={'xy2': 'xy2', 'xy': 'xy',
                  'yz' : 'yz' , 'xz': 'xz',},
        vectors=None,
        # select data to plot
        tstart=0., tend=1e38, islice=-1,
        istart=None, iend=None,
        # color_range size 2 list cmin and cmax
        color_range=None, color_levels=None,
        # convert data to cgs from code units and rescale to cbar_label
        # par if present is a param object
        unit='unit_velocity', rescale=1., par=list(),
        # Extras
        debug=True,
        # Any parameters __plot() takes
        settings=PlotSettings(),
        ) -> None:    
    """
    Create plot from slice data in cartesian | spherical | cylindrical
    coordinates and saves the output as images and / or videos.
    Creates plots for times [tstart, tend] of all the fields specified in
    fields argument. Internally calls __plot_field() to handle plotting of a 
    given field over all time instants.
    By default this function creates images off screen without creating
    interactive plot window. If plotting done on screen (i.e. off_screen=False)
    and creating a video, one can see the video live and the video angle can
    be rotated, also affecting the output video.
    
    In case of fine tuning the plot parameters, e.g. camera position / focal point, 
    one can turn preview on (in settings) to see an interactive plot window
    that can be rotated around and outputs the camera parameters on both command
    line and the plot window.

    Parameters
    ----------
    slice_obj : pencil.read.allslices.SliceSeries, optional
        Pencil slice data, if supplied no slice data is read. If None, slice data is
        automatically read from datadir. Default: None.
    datadir : str, optional
        Directory of the data. Default: './data'.
    precision : str, optional
        Precision for the data read, parameter e.g. pencil.read.slices() takes in.
        Can be 'half', 'f' or 'd'. Default: 'f'.
    fields : list of strings, optional
        Fields that are plotted. Default: ['uu1',].
    xyzplane : dictionary, optional
        Dictionary having allowed keys defined by XYZPLANE_KEYS. Values should be the
        matching data that is plotted on each of these surfaces. See "Coordinate systems"
        below.
    vectors : list, optional
        List of 3 strings defining the vector components of the vector field. This
        is used to create streamlines / vectors.
    tstart : float, optional
        Start time instant. See parameter islice for more details. Default: 0.
    tend : float, optional
        End time instant. See parameter islice for more details. Default: 1e38
    islice : int, optional
        Sequential integer number of a slice in given period (starting from 0).
        If set to -1, generates all slices in given period [tstart, tend]. 
        Else generates only visualisations of the given slice number islice. 
        Default: -1.
    istart : int, optional
        First index to start plotting from. This index corresponds to time instant
        slice_obj.t[istart]. istart is used if it is not None, otherwise tstart
        is used.
    iend : int, optional
        Similar to istart, end time index corresponding to time slice_obj.t[iend]
    color_range : list of 2, optional
        If supplied, list containing [cmin,cmax] for the colorbar.
    color_levels : string, optional
        If set to 'common', then all times and fields have same colorbar min
        and max values. Default: None.
    unit : string, optional
        TODO!
    rescale : int, optional
        TODO!
    par : list, optional
        TODO!
    debug : boolean, optional
        Enable debugging prints. Default: False.
    **settings : optional
        Dictionary of possible settings for the plot. See docstring at the start of this
        file for an ready dictionary input of all possible parameters settings could
        take in.

    Coordinate systems
    ------------------
    Allowed keys for xyzplane are defined by XYZPLANE_KEYS = ['xy', 'xy2', 'xz', 'yz'].
    In case of each coordinate system, these surfaces correspond to:
    
    - Cartesian: 
        * xy = bottom
        * xy2 = top 
        * xz = right vertical slice
        * yz = left verical slice
    
    - Spherical: 
        * xy = r-theta slice, 1st meridional slice
        * xy2 = r-theta slice, 2nd meridional slice
        * xz = r-phi slice, "plane slice" between meridionals and yz
        * yz = theta-phi slice, spherical surface at the back (or part of it)
        
    - Cylinder:
        * xy = r-theta slice, bottom of the cylinder (circle)
        * xy2 = r-theta slice 2, top of the cylinder (circle)
        * xz = r-z slice, "radial cut", rectangle between top/bottom
        * yz = theta-z slice, shell of the cylinder (or part of it)

    Examples
    --------    
    Minimal usage example, this saves by default 'png' images at all timesteps of field 'uu1'
    to figure directory ./images (creating the directory if it does not exist). Debug should be
    left on so one can see that the script runs properly. By default this creates a cartesian
    plot (i.e. box plot).

    >>> from pencil.visu.pv_plotter import plot, PlotSettings
    >>> plot(debug=True)

    This plotter uses a settings object (dataclass) PlotSettings that contains all
    possible extra settings that can be varied. See documentation on PlotSettings 
    for further information on all of the parameters. This can be used in the 
    following way:
    
    >>> settings = {
                'imageformat': None,
                'videoformat': 'mp4',
                'coordinates': 'cartesian'
                'cbar_title': 'TITLE', 
                'n_colors': 100,
                'camera_centre': (-170, -140, 125),
                'focal_point': (30, 30, 16),
                'offset': 2.5
                }
    >>> settings = PlotSettings(**settings)
    >>> plot(debug=True settings=settings)
    
    Furthermore, 2D slices to be plotted can be modified by passing in xyzplane parameter
    
    >>> xyzplane = {'xy2': 'xy2','xy': 'xy', 'xz': 'xz', 'yz': 'yz'}
    
    where the string values represent the slices to be plotted on the meshes (see "Coordinate
    systems" for what which mesh each key corresponds given coordinate system). These values 
    should be found in slices.__dict__.keys().
    
    Notes
    -----
    - NOTE! The default camera angle is most likely horrible and may not even show
    the data at all. This should be manually modified and e.g. easily found using
    preview=True that shows an interactive plot window that can be rotated around
    and print out the camera parameters at the same time.
    - NOTE! On CSC computers module 'mesa-settings' needs to be loaded (works at
    least on Puhti).    
    """
    ############################################################################
    ### SANITY CHECKS                    
    if not isinstance(xyzplane, dict):
        raise TypeError(f'xyzplane is of type {type(xyzplane)}, should be dict.')
    
    for key in xyzplane.keys():
        if key not in XYZPLANE_KEYS:
            raise KeyError(f'Invalid key {key} in argument xyzplane. Key should' 
                           'be one of: {XYZPLANE_KEYS}')
    
    if settings.coordinates not in ['cartesian', 'spherical', 'cylinder']:
        raise NotImplementedError(f'Unknown coordinate system: {settings.coordinates}!'
                                  f'Should be "cartesian", "spherical" or "cylinder".')
    
    ############################################################################
    
    if slice_obj == None:
        if vectors is not None:
            # Slice object needs to contain both the scalar field and necessary 
            # vector field components
            slice_fields = list(set(fields + vectors))
        else:
            slice_fields = fields 
        if debug:
            print(f'Starting to read slices...')
            print(f'> Reading fields: {slice_fields}')
            T0 = time.time()
        slice_obj = pc.read.slices(datadir=datadir, field=slice_fields, precision=precision)
        if debug:
            T1 = time.time()
            print(f'Reading slices took: {T1-T0} seconds')
    else:
        print("Using supplied slice_obj!")
    if debug:
        print(f'Slice object contains:')
        print(f'> slices = {slice_obj.__dict__.keys()}')
        print(f'> Fields = {slice_obj.xy.__dict__.keys()}')
        
    if vectors is not None:
        assert len(vectors) == 3 and isinstance(vectors, list), f'Vectors should'
        'be a list of 3 elements'
        
        if settings.streamlines and settings.surface_vectors:
            print(f'WARNING! Having surface vectors on zeroes out the third component'
            'but it is useful for streamline variable radius scaling (i.e. scaling would'
            'be based on all 3 components rather than 2).')
            print(f'>>> Automatically settings settings.surface_vectors to False'
            'in order to avoid zeroing out the 3rd component!')

            settings.surface_vectors = False

        vectors = __vectorsFromSlices(slice_obj, vectors, datadir,
                                        vectors_unit_length=False, 
                                        coordinates=settings.coordinates,
                                        surface_vectors=settings.surface_vectors)
    if istart != None and iend != None:
        tstart = slice_obj.t[istart]
        tend   = slice_obj.t[iend]
    ttmp = slice_obj.t[np.where(slice_obj.t <= tend)[0]]
    it = np.where(ttmp >= tstart)[0]    # time indices (boolean) such that t in [tstart, tend]

    if isinstance(xyzplane, dict) and len(xyzplane.values()) < 4:
        raise ValueError(f'xyzplane: plotting requires at least 4 surfaces. Only'
                         f'{len(xyzplane.values)} given in xyzplane.')

    # avoid increasing memory
    dtype = type(slice_obj.__getattribute__(xyzplane[list(xyzplane.keys())[0]])
                          .__getattribute__(fields[0])[0, 0, 0])

    if dtype == np.float16 or dtype == 'half':
        print(f'Caution! dtype {dtype} may cause under/overflow')

    for field in fields:  # fields = ['uu1', ...]
        if not isinstance(par,list) and len(unit) > 0:
            unitscale = par.__getattribute__(unit)*rescale
        else:
            unitscale = 1.

        # Read in slices
        scalars = dict()
        for key in xyzplane.keys():   # ['z2', 'z1', 'yz', 'xz', ...]
            if 'ln' in field:
                scalars[key] = slice_obj.__getattribute__(xyzplane[key]).__getattribute__(
                                            field).astype(dtype) + np.log(unitscale)
            else:
                scalars[key] = slice_obj.__getattribute__(xyzplane[key]).__getattribute__(
                                            field).astype(dtype) * unitscale

        if color_range is None and color_levels == 'common':
            # For current field set common [cmin,cmax] for all time instants
            cmin, cmax = 1e38, -1e38
            for key in scalars.keys():       # Loop over ['z2', 'z1', ...]
                cmax = max(cmax, scalars[key][it].max())
                cmin = min(cmin, scalars[key][it].min())

        elif color_range is None:
            # Color range not supplied, set automatically in __plot_field
            cmin = None
            cmax = None
        else:
            # Color_range is supplied
            cmin = color_range[0]
            cmax = color_range[1]

        # '-1' to generate all slices in given period
        if islice == -1:
            __plot_field(
                slice_obj, scalars, field, it, vectors,
                datadir=datadir, cmin=cmin, cmax=cmax, 
                dtype=dtype, settings=settings, debug=debug)

        # Generate only slice corresponding to index islice
        else:
            __plot_field(
                slice_obj, scalars, field, [islice], vectors,
                datadir=datadir, cmin=cmin, cmax=cmax, 
                dtype=dtype, settings=settings, debug=debug)


################################################################################
