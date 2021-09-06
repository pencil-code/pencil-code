#
# pv_volume_plotter.py
#
# This file defines 3D plotting routines for VAR files using PyVista.
#
# Author: Leevi Veneranta (leevi.veneranta@aalto.fi)
#
"""
`pv_volume_plotter.py` defines multiple different plotting routines for volume 
data from Pencil code VAR files. This includes methods such as:

* scalars --- plots only scalar values on the plotter. Furthermore, scalars can
            be just 2D surfaces, or a volume rendering (only in cartesian coordinates)
            or 3-orthogonal slices defined by the center of the 3 slices.

* vectors --- Vectors can be overlayed on to the scalar data. Vectors can be either
            constant in size or scaled by the magnitude of the vector field or the 
            values on the scalar data.

* streamlines --- Similarly to vectors, streamlines can be overlayed on to the scalar
                data. The radius of the streamlines can be scaled e.g. by the cartesian
                norm of the vector field. In addition, a colormap can be applied on to
                the streamline tubes.

    In addition, there are multiple source point sampling strategies for the vector
    field visualization methods (vectors|streamlines) such as:
    * random points in a volume
    * points on the surface of the volume, useful if streamlines are set to be 
        constrained on to the surface (by setting parameter `surface_streamlines`)
    * points on a plane that slices through the volume
    * points on 3-orthogonal planes that slice through the volume

* isosurfaces --- Plots isosurfaces defined by its isovalues. These can be either
                manually set or calculated automatically between the min max range
                of the supplied scalars. Furthermore, a moving isovalue gif can be
                generated that moves from first isovalue to last and back showing
                an moving representation of the isosurfaces.

For usage examples and tutorials, see Pyvista3DPlot docstring!


Requirements
------------
* sklearn
* tqdm
* numpy
* pyvista (should be at least >= 0.31.3)

All of these can be installed by:
```
pip3 install <Library name>
```


WARNING
-------
I have experienced issues with volume rendering when using linux in VirtualBox due
to missing 3D acceleration. HyperV seems to be able to access system hardware better
and works on that. In my case I have a dedicated NVIDIA GPU - with different system
hardware volume rendering might start throwing VTK warnings/errors and OpenGL related
issues.


NOTES
-----
- Note! In order for `streamlines_from_source()` to work you need to update PyVista
    to the latest version (at the moment 0.31.3)


Known Issues
------------
- Interactive plotting windows have trouble closing sometimes - based on pyvista
    Github this is a known issue in the VTK source.
"""
# External libraries
from sklearn.model_selection import ParameterGrid
from tqdm import tqdm
import pyvista as pv
import pencil as pc
import numpy as np

# Python libraries
from dataclasses import dataclass, asdict
from pathlib import Path
import time
import json
import csv
import os

from pencil.visu.pv_plotter_utils import *

# Libraries used for testing / Debugging
# --> icecream provides nice pretty printing and supports different data structures
# --> memory-profiler has tools for easy python script memory profiling
# from memory_profiler import profile, memory_usage
from icecream import ic


# Constant dictionary defining keys (integer error key) and values (reason for
# terminating streamline integration). VTK internals since Pyvista streamlines
# uses vtk.StreamTracer.
REASONS_FOR_TERMINATION = {
    pv._vtk.vtkStreamTracer.OUT_OF_DOMAIN      : 'OUT_OF_DOMAIN',
    pv._vtk.vtkStreamTracer.NOT_INITIALIZED    : 'NOT_INITIALIZED',
    pv._vtk.vtkStreamTracer.UNEXPECTED_VALUE  : 'UNEXPECTED_VALUE',
    pv._vtk.vtkStreamTracer.OUT_OF_LENGTH      : 'OUT_OF_LENGTH',
    pv._vtk.vtkStreamTracer.OUT_OF_STEPS        : 'OUT_OF_STEPS',
    pv._vtk.vtkStreamTracer.FIXED_REASONS_FOR_TERMINATION_COUNT
                                    : 'FIXED_REASONS_FOR_TERMINATION_COUNT',
}

SAMPLING_METHODS = [
    'surface points', 
    'points', 
    'slice', 
    'slice orthogonal',
    'default'
    ]


def printTerminationReasons(arr, plotter, add_to_image=False):

        def countTerminationReasons(arr):
            counts = dict()
            for i in arr:
                if i not in counts.keys():
                    counts[i] = 1
                else:
                    counts[i] += 1
            return counts

        counts = countTerminationReasons(arr)
        total = sum(counts.values())

        print(f'\nReasons for termination\n---------------------------------')
        title = ''
        
        for key in sorted(counts.keys()):
            reason = REASONS_FOR_TERMINATION[key]
            count  = counts[key]
            print(f'> [{key}] {reason} = {count} ({round(100*count / total, 2)} %)')
            
            if add_to_image:
                header = True
                if header:
                    title += f'Reasons for termination:\n'
                    header = False
                title += f'{reason} = {count} ({round(100*count / total, 2)} %)\n'
            
        print(f'---------------------------------\n')
        
        if add_to_image:
            plotter.add_text(text=title, font_size=8, position='upper_left')
            return plotter
        

def meshSamplingMethods(method, mesh, surface_mesh, n_points=100, normal='x', 
                        origin=None, set_seed_1=False, alwaysReturnList=False,
                        get_arrays=False):
    
    if method == 'surface points':
        src = randomSampleMeshPoints(n_points, surface_mesh, set_seed_1=set_seed_1, 
                                     get_arrays=get_arrays)
    
    elif method == 'points':
        src = randomSampleMeshPoints(n_points, mesh, set_seed_1=set_seed_1, get_arrays=get_arrays)
    
    elif method == 'slice':
        src = mesh.slice(normal=normal, origin=origin)
        src = randomSampleMeshPoints(n_points, src, set_seed_1=set_seed_1, get_arrays=get_arrays)  
    
    elif method == 'slice orthogonal':
        if origin is not None:
            x, y, z = origin[0], origin[1], origin[2]
        else:
            x, y, z = None, None, None
            
        src = mesh.slice_orthogonal(x=x, y=y, z=z)
        for i in range(len(src)):
            src[i] = randomSampleMeshPoints(n_points, src[i], set_seed_1=set_seed_1, 
                                            get_arrays=get_arrays)
        # In this case slice_orthogonal returns a pyvista.MultiBlock which is basically
        # a list of PolyData. Hence no list conversion necessary so return src here!
        return src
    
    else:
        raise ValueError(f'Unknown sampling method {method}. Should be one of'
                         f' the following: {SAMPLING_METHODS}')

    if alwaysReturnList:
        if not isinstance(src, list):
            src = [src]
            
    return src    


@dataclass
class Plot3DSettings:
    """
    Class holding all plot customization parameters for the Pyvista3DPlot class.
    NOTE! See Pyvista3DPlot docstrings for documentation on the parameter tests.
    
    Parameters
    ----------
    title: str
        Title text to be added to the upper left corner. If set to '' title is set
        automatically to contain some information, i.e. datadir, scalar key and 
        vector keys. If you want title to be empty, set it to e.g. '   '
    title_font_size: int
        Font size for the title text
    n_points: int
        Total number of source points (vectors | streamlines)
    set_seed_1: bool
        If True, random seed is set to one, all sampled source points are the same
    imageformat: str
        Any Python Pillow compatible image format, 'png', 'jpeg', etc.
    off_screen: bool
        Plotting done on screen or off screen?
    window_size: tuple or list
        Window size.
    method: str
        Sampling method, how source points for vectors / streamlines are sampled: 
            1.'surface points' --- random sampled points from the surface of the
                                    mesh
            2. 'points' --- random sampled points in the whole mesh
            3. 'slice' --- random sampled points on a slice defined by normal and
                            an origin
            4. 'slice_orthogonal' --- random sampled points on 3 orthogonal slices 
                                        defined by origin
            5. 'default' --- only to be used with streamlines, randomly samples 
                                 mesh in a sphere given radius.
    show_axes: bool
        Show axes for the plot
    show_grid: bool
        Shows a grid on the plotter.
    scalars_norm: str
        Either 'log' or None
    override_log: bool
        If False, and scalar_key contains 'ln' parameter scalars_norm is automatically
        set to None. Else if scalar_key does not contain 'ln', scalar_norm is set
        to 'log'. If True, scalars_norm is kept to whatever user has set it to.
    background_color: str or 3-tuple of ints
        Either a string e.g. 'white', 'black' or 'gray'. Also can be a three int
        RGB tuple e.g (255, 255, 255)

    Widgets
    -------
    add_volume_opacity_slider: bool
        Adds a slider widget controlling parameter `volume_opacity`. Only applies
        when `mesh_type=="volume"`, i.e. volume rendering enabled.
    widget_type: type
        Both of the available widget methods show an "arrow" that can be used to
        choose the angle of the sliced plane. Furthermore the plane itself can be
        moved up and down in the whole volume. 
        
        Available widget types are:
        1) 'clip slice' --- Only shows the current slice, noticeably faster than
                            clip box. 
        2) 'clip box' --- Shows the current slice + the rest of the data "below"
                            the slice.
        3) 'plane vectors' --- Similar to 'clip box' except also vectors are added
                            on to the slice. The vectors can be controlled by changing
                            Plot3DSettings vectors parameters.
    Camera
    ------
    camera_centre: tuple of 3 floats
        Coordinates for the centre of the camera.
    focal_point: tuple of 3 floats
        Coordinates for the camera focal point
    view_up: tuple of 3 floats
        Vector defining up-direction of the camera. By default (0,0,1)
        i.e. "up" (z-direction).
    zoom_factor: float
        Zoom multiplier added to the camera. Larger the value, the more it zooms
        in.

    Streamlines
    -----------
    tube_radius: float
        Radius for stream tubes
    show_source: bool
        Show source points for streamlines
    surface_streamlines: bool
        If true streamlines are confined to a 2D surface
    stream_params: dict
        Any parameters pyvista streamlines_from_source() takes in OTHER THAN
        `surface_streamlines` and `vectors`. The following is the Pyvista's own
        documentation on the pyvista.streamlines_from_source, see link for documentation
        after the list:
            * integrator_type
                The integrator type to be used for streamline generation. The 
                default is Runge-Kutta45. The recognized solvers are: RUNGE_KUTTA2 
                (2), RUNGE_KUTTA4 (4), and RUNGE_KUTTA45 (45). Options are 2, 4, 
                or 45. Default is 45.
            * integration_direction
                Specify whether the streamline is integrated in the upstream or 
                downstream directions (or both). Options are 'both', 'backward', 
                or 'forward'.
            * initial_step_length
                Initial step size used for line integration, expressed in length
                unitsL or cell length units (see step_unit parameter). either the
                starting size for an adaptive integrator, e.g., RK45, or the 
                constant / fixed size for non-adaptive ones, i.e., RK2 and RK4).
            * step_unit
                Uniform integration step unit. The valid unit is now limited to 
                only LENGTH_UNIT ('l') and CELL_LENGTH_UNIT ('cl'). Default is 
                CELL_LENGTH_UNIT: 'cl'.
            * min_step_length
                Minimum step size used for line integration, expressed in length 
                or cell length units. Only valid for an adaptive integrator RK45.
            * max_step_length
                aximum step size used for line integration, expressed in length 
                or cell length units. Only valid for an adaptive integrator RK45.
            * max_steps
                Maximum number of steps for integrating a streamline.
            * terminal_speed
                Terminal speed value, below which integration is terminated.
            * max_error
                Maximum error tolerated throughout streamline integration.
            * max_time
                Specify the maximum length of a streamline expressed in LENGTH_UNIT.
            * compute_vorticity
                Vorticity computation at streamline points (necessary for generating
                proper stream-ribbons using the vtkRibbonFilter.
            * interpolator_type
                Set the type of the velocity field interpolator to locate cells 
                during streamline integration either by points or cells. The cell
                locator is more robust then the point locator. Options are 'point'
                or 'cell' (abbreviations of 'p' and 'c' are also supported).
            * rotation_scale
                This can be used to scale the rate with which the streamribbons 
                twist. The default is 1

        See https://docs.pyvista.org/core/filters.html?highlight=streamlines_from_source#pyvista.DataSetFilters.streamlines_from_source
        for more details on the parameters.

        Note that if stream_params is None the following defaults are set:
        >>> self.stream_params = {
                'max_steps': 1000,
                'max_time': 1e60,
                'terminal_speed': 1e-60,
                'integration_direction': 'both',
                'compute_vorticity': False,
                'integrator_type': 45,
            }
    
    The following parameters are used if method == 'default'. Then streamline source points are 
    initialized from a sphere of radius source_radius and with center source_center.
    
    source_radius: float
        See description above. Radius for the sphere from which source points are 
        sampled.
    source_center: tuple
        Should be a tuple of 3 floats. Center for the sphere where source points
        are sampled from.
    stream_variable_radius: str
        Enable variable radius on the streamlines. Name of the data array used to scale the radius,
        either 'vector magnitude', 'scalars' or None
    stream_radius_factor: float
        Multiplicative factor for the stream variable tube radius scaling.

    Vectors
    -------
    vector_scaling: str
        Scaling method for vectors, options:
            1. None --- no scaling
            2. 'scalars' --- scaled based on scalars no the mesh
            3. 'magnitude' --- scaled based on vector magnitude
    vector_factor: float
        Multiplicative scaling for vectors
    vector_max: float
        Maximum value for vector magnitude scaling. If None, set automatically 
        to mean + one standard deviation
    
    Mesh
    ----
    show_mesh: bool
        Whether to show the mesh or not
    mesh_type: str
        Rendering type for mesh, options: "volume", "orthogonal slices" or "surface"
    slice_pos: tuple
        Applies if mesh_type=='orthogonal slices'. 3-tuple of floats (x,y,z) where
            - x = The X location of the YZ slice
            - y = The Y location of the XZ slice
            - z = The Z location of the XY slice
    volume_opacity: float
        Value between [0,1]. 0 corresponds to "not see through" i.e. opacity=1 and
        1 "see through" i.e. opacity=0. In reality it is not that straight forward
        since this is handled automatically by volume renderer. 

        What this parameter does is that `pyvista.add_volume parameter` `unit_opacity_distance`
        is set to value `mesh.length * volume_opacity`. See pyvista Documentation
        for more information on `unit_opacity_distance`. The default parameter 
        1/25 = 0.04 seems to work fairly well.
    culling: str
        Does not render faces that are culled. Options are 'front' or 'back'. 
        This can be helpful for dense surface meshes, especially when edges are 
        visible, but can cause flat meshes to be partially displayed. By default
        None, i.e. no culling applied.
    
    Scalarbar: General
    ------------------
    vertical_sbar: bool
        Scalarbar vertical or horizontal
    cbar_width: float
        Width as a percentage, between 0 and 1
    cbar_height: float
        Height as a percentage, between 0 and 1
    cbar_title_font: int
        Title font size for scalarbar title.
    cbar_label_font: int
        Label font size for scalarbar labels.
    n_colors: int
        Number of colors for the colorbar.
    _sbar_args: dict
        INTERNAL VALUE, SHOULD NOT BE USED
    annotation_color: str
        Color for the annotations, i.e. scalarbar tick texts etc.
    
    Scalarbar: Mesh
    ---------------
    scalar_sbar_title: str
        Title for the scalarbar for the mesh scalars.
    show_mesh_sbar: bool
        Show scalarbar for the scalars on the mesh
    mesh_opacity: float
        Opacity for the mesh
    mesh_cmap: str
        matplotlib compatible colormap for the mesh
    scalar_sbar_pos_x: float
        The percentage (0 to 1) along the windows’s horizontal direction to place 
        the bottom left corner of the colorbar. If None it is placed automatically.
    scalar_sbar_pos_y: float
        The percentage (0 to 1) along the windows’s vertical direction to place 
        the bottom left corner of the colorbar. If None it is placed automatically.

    Scalarbar: Vectors | Streamlines
    --------------------------------
    field_sbar_title: str
        Title for the scalarbar for the mesh field (vectors|streamlines).
    show_field_sbar: bool
        Show scalarbar for the vectors on the mesh (streamlines/vectors)
    field_opacity: float
        Opacity for the plotted streamlines / vectors
    field_cmap: str
        Colormap for the plotter streamlines / vectors
    field_sbar_pos_x: float
        The percentage (0 to 1) along the windows’s horizontal direction to place 
        the bottom left corner of the colorbar. If None it is placed automatically.
    field_sbar_pos_y: float
        The percentage (0 to 1) along the windows’s vertical direction to place 
        the bottom left corner of the colorbar. If None it is placed automatically.

    Orbit gif
    ---------
    Orbit gif creates a gif that circles around the plotted mesh. The orbital 
    path is always 360 of how many points you have.
    
    orbit_gif: bool
        Whether to create an orbiting gif
    orbit_points: int
        Number of points on the orbiting path
    orbit_step: float   
        Timestep between flying to each camera position
    
    Slices
    ------
    This defines the parameters for the two slicing methods: 'slice orthogonal' 
    and 'slice'. Note that this only affects sampled source points for vectors or 
    streamlines.
    
    normal: str
        tuple(float) or str. Length 3 tuple for the normal vector direction. Can
        also be specified as a string conventional direction as 'x' for (1,0,0) 
        or '-x' for (-1,0,0) etc.
    origin: tuple
        The center (x,y,z) coordinate of the plane on which the slice occurs. Note
        that if method == 'slice orthogonal' then we have that:
            * x corresponds to x location of the YZ slice
            * y corresponds to y location of the XZ slice
            * z corresponds to z location of the XY slice
    """
    ### General
    title: str    = ''
    title_font_size: int = 14
    n_points: int = 100
    set_seed_1: bool = False
    imageformat: str = 'png'
    off_screen: bool = False
    window_size: tuple = (1008,656)
    method: str = 'points'
    show_axes: bool = True
    show_grid: bool = False
    scalars_norm: str = None
    override_log: bool = False
    background_color: str = 'white'

    ### Widgets
    add_volume_opacity_slider: bool = False,
    widget_type: str = None

    ### Camera
    camera_centre: tuple = None
    focal_point: tuple = None
    view_up: tuple = (0, 0, 1)
    zoom_factor: float = None

    ### Streamlines 
    tube_radius: float = 0.005
    show_source: bool = False
    surface_streamlines: bool = False
    stream_params: dict = None
    source_radius: float = None 
    source_center: tuple = None 
    stream_variable_radius: str = None
    stream_radius_factor: float = 1

    ### Vectors
    vector_scaling: str = 'magnitude'
    vector_factor: float = 1
    vector_max: float = None

    ### Mesh
    show_mesh: bool = True
    mesh_type: str = "surface"
    slice_pos: tuple = (None, None, None)
    volume_opacity: float = 0.04
    culling: str = None

    ### Scalarbar: General
    vertical_sbar: bool = True
    cbar_width: float = 0.06
    cbar_height: float = 0.65
    cbar_title_font: int = 10
    cbar_label_font: int = 8
    n_colors: int = 256
    _sbar_args: dict = None
    annotation_color: str = 'black'
    
    ### Scalarbar: mesh
    scalar_sbar_title: str = None
    show_mesh_sbar: bool = True
    mesh_opacity: float = 0.2
    mesh_cmap: str = 'bwr'
    scalar_sbar_pos_x: float = None # TODO! Not used
    scalar_sbar_pos_y: float = None # TODO! Not used

    ### Scalarbar: Vectors | Streamlines
    field_sbar_title: str = None
    show_field_sbar: bool = True
    field_opacity: float = 1.0
    field_cmap: str = 'bwr'
    field_sbar_pos_x: float = None # TODO! Not used
    field_sbar_pos_y: float = None # TODO! Not used

    ### Orbit gif
    orbit_gif: bool = False
    orbit_points: int = 35
    orbit_step: float = 2.0
    
    ### Slices
    normal: str = 'x'
    origin: tuple = None


    def __post_init__(self):
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
            'vertical': self.vertical_sbar,
            'width': self.cbar_width, 
            'height': self.cbar_height, 
            'color': 'black',
            'title_font_size': self.cbar_title_font,
            'label_font_size': self.cbar_label_font,
            'n_colors': self.n_colors,
            'color': self.annotation_color,
        }
        
    
    def update(self, **kwargs):
        """
        For updating any key-value pair defined by the self.settings.
        """
        for key, value in kwargs.items():
                setattr(self, key, value)


class Pyvista3DPlot:
    """
    Implements 3D plotting routines from pencil var files. It is able to generate
    streamlines, vectors and isosurfaces. Optionally one can save the output as an
    image (any Python Pillow compatible image type) or as a gif that orbits (circles)
    around the plot.
    
    Compatible coordinate systems are both cartesian and spherical coordinates.
    Note that there are multiple ways of sampling source points for the vector field
    (i.e. streamlines|vectors) and multiple ways of adding the scalars to the plot
    e.g. volume rendering (only cartesian supported), 3 orthogonal slices or surfaces.
    
    
    Parameter tests
    ---------------
    If at least one of the input settings is given as a dictionary parameter test
    will be enabled automatically. This enables one to input multiple values for 
    settings in a list, and when running streamlines | vectors all these parameter
    combinations are tried automatically. This is implemented using 
    ``sklearn.ParameterGrid``.

    * Note that only methods streamlines() and vectors() support parameter tests.

    * This can be useful for finding good parameters from a large set of options.

    * Log file is written containing an id and all settings matching to the outputted
        pictures (id is added to the output image names too).


    Interface
    ---------
    This class defines the following plotting utilities for the user:

    * scalars --- Adds just the scalar values to the visualization. Scalars can be added
                as surfaces, volume rendering (only cartesian) and 3 orthogonal 
                slices. This is defined by Plot3DSettings.mesh_type.

    * streamlines --- Adds streamlines to visualization. In addition scalars are added
                    to mesh in wanted way, i.e. volume rendering can be used in
                    conjunction with streamlines.

    * vectors --- Adds vectors to visualization. Similarly to streamlines, scalars
                    can be added, e.g. as a volume rendering (only cartesian)

    * contour --- Adds isosurfaces to visualization. These can be either defined 
                manually or the number of isosurfaces can be set generating uniformly
                spaced isosurfaces between minimum and maximum value of scalars.

    * movingIsovalue --- Similar to contour except outputs a "moving" gif visualization
                        that loops from first isosurface to last and backwards.

    Other utility functions:

    * writeSettings --- can be used to write the current Plot3DSettings to a json
                        file.

    * loadSettings --- can be used to load values for Plot3DSettings from a json
                        file.

    * saveToVTK --- saves the current mesh to a VTK file. This can be used later,
                    e.g. loaded in by PyVista.

    * updateSettings --- can be used to update one or multiple values of Plot3DSettings
                        without having to reinitialize the whole Pyvista3DPlot object.
    
    * preview --- Shows an interactive window that can be panned around showing
                the camera parameters. This can be used to find good camera parameters.

    
    Examples
    --------
    1) Very simple usage example
    
    Assuming your data is in `./data` directory:
    >>> from pencil.visu.pv_volume_plotter import Pyvista3DPlot, Plot3DSettings
    >>> plotter = Pyvista3DPlot(vector_key=['ux','uy','uz'], scalar_key='rho', 
                                outputdir='./stream_output', debug=True)
    
    Then plot streamlines 
    >>> plotter.streamlines()

    This should show an interactive plotter window. If not parameter `off_screen`
    might be true. One way to change this is:
    >>> plotter.updateSettings(off_screen=False)

    Now all plots should be done "on screen" that is interactively. When you are
    plotting surfaces, that is `method='surface'` it might be useful to turn mesh
    opacity to full 1.0. But then if you have for example streamlines inside the
    volume these are not visible? This can be fixed by setting parameter
    `culling='back'` from Plot3DSettings. This means that the front meshes in front
    of the camera are not rendered at all, only showing the back meshes.

    How to do this? 
    >>> plotter.updateSettings(culling='back', mesh_opacity=1.0)
    >>> plotter.streamlines()

                            
    2) VOLUME RENDERING SCALARS WITH STREAMLINES OR VECTORS. 
    
    This assumes your data is in cartesian coordinates. First create the settings object:
    >>> settings = {
        'off_screen': False,                          # Plotting should be done on screen
        'method': 'points',                           # Ensure streamline source points are taken from the whole volume
        'scalars_norm': 'log',                        # This might be useful for the scalar data
        'stream_variable_radius': 'vector magnitude', # Enable stream radius scaling
        'stream_radius_factor': 4,                    # Might be necessary to add multiplicative factor if the radius 
                                                      #   doesnt show up (it is so small)
        'mesh_type': 'volume',                        # Enable volume rendering
        #'volume_opacity': 1/25,                      # Change this value if rendered volume is very faint
    }
    
    Now that we have settings in dictionary we can initialize the settings object:
    >>> settings = Plot3DSettings(**settings)
    >>> plotter = Pyvista3DPlot(vector_key=['ux','uy','uz'], scalar_key='rho',
                                outputdir='./stream_output', settings=settings, debug=True)
    
    And finally create the streamlines + volume rendering of scalars
    >>> plotter.streamlines() 
    
    Furthermore you could try adding vectors instead of streamlines:
    >>> plotter.vectors()
    
    If this looks bad, you need to tune the vector parameters, example:
    >>> plotter.updateSettings(vector_scaling='magnitude', 
                                vector_factor=1e-2, # Depending on data you need to adjust this smaller or larger
                                                    #   in order to make all vectors smaller or larger
                                vector_max=13, # If the maximal length of vectors is way too large, tone this down
                                )
    Now try again:
    >>> plotter.vectors()


    3) ISOSURFACES EXAMPLE:
    >>> plotter = Pyvista3DPlot(vector_key=['ux','uy','uz'], scalar_key='rho',
                                    outputdir='./stream_output', debug=True)
    >>> plotter.contour(isosurfaces=10)
    
    Note, contour might sometimes print out warnings that isosurface has no points,
    this is most likely due to the isosurface value doesn't correspond to any existing
    value in the dataset itself. Therefore, it might be more useful to specify manually
    the isovalues by setting the `values` parameter to some known values.
    

    Notes on some missing functionality
    -----------------------------------
    * At the time of development pyvista.add_volume (volume rendering) method supported
        only pyvista.UniformMesh which made it difficult if not impossible to add
        volume renderings of data in spherical coordinates.
    """
    def __init__(self, vector_key, scalar_key, datadir='./data', 
                     precision='f', magic='bb', ivar=-1, coordinates='cartesian',
                     outputdir=f'./stream_output', 
                     settings=Plot3DSettings(), debug=False):
        """
        Initialization for the Pyvista3DPlot object.
        
        Parameters
        ----------
        vector_key: str or list of strings
            Keys defining the vector field. Either a single string e.g. 'bb' or 
            a list of 3 strings e.g. ['uu1', 'uu2', 'uu3']. This is read from the 
            var file.
        scalar_key: str
            Key defining the scalars in the volume. This is read from the var file.
        datadir: str, optional
            Path to the data directory. By default './data'
        precision: str, optional
            Precision for the pencil.read.var routine. By default 'f'.
        magic: str, optional
            Magic for the pencil.read.var routine. By default 'bb'.
        ivar: int, optional
            Index of the var file read. By default -1.
        coordinates: str, optional
            Either 'cartesian' or 'spherical'. By default 'cartesian'
        outputdir: str, optional
            Path to the directory where the output is saved (images / gif). By
            default './stream_output'
        settings: Plot3DSettings, optional
            Plot3DSettings object containing further settings for the plots. See
            source code for possible arguments, explanations and default values.
            By default Plot3DSettings().
        debug: bool, optional
            Enable debug printing. By default False.
        """
        self.datadir = Path(datadir)
        self.outputdir = Path(outputdir)
        if not self.outputdir.exists():
            os.mkdir(self.outputdir)
        assert self.datadir.exists(), f'Given data directory does not exist!' \
        f' Directory given was: {self.datadir}'
        
        self.debug = debug
        
        if self.debug:
            print(f'> Starting to load settings...')
        
        self.loadSettings(settings)

        # Set automatic title containing data directory, scalar and vector keys
        if settings.title == '' or settings.title == None:
            settings.update(title=f'{str(self.datadir)}\nscalar: {scalar_key}\n'
                                f'vector: {vector_key}')
        
        # If override_log is False set normalization on scalars automatically based
        # on whether the key contains 'ln' or not
        if not settings.override_log:
            if 'ln' in scalar_key:
                settings.update(scalars_norm=None)
            else:
                settings.update(scalars_norm='log')

        self.scalar_key = scalar_key
        self.vector_key = vector_key

        self.var = pc.read.var(datadir=datadir, precision=precision, trimall=True, 
                                magic=magic, ivar=ivar)
        self.grid = pc.read.grid(datadir=datadir, trim=True)

        self.coordinates = coordinates
        self.plotter = None
        self.mesh = None
        self.src = None
        
        # id_counter is added to output images names, helps to identify paramtest
        # parameters to images
        self.id_counter = 0

        # Represents either streamlines / or vectors of the vector field. Done this
        # way so we can have only one function call in __handlePlotter() to add
        # this mesh to plotter
        self.field = None
        if self.debug: 
            print("> Starting to load data...")
        self.__3D_data_loader(scalar_key, vector_key)
        self.surface_mesh = self.mesh.extract_surface()

        if self.debug:
            print(f'> var file has keys: \n{self.var.__dict__.keys()}\n')

    
    def writeSettings(self, filename):
        """
        Writes current Plot3DSettings to a json file for later use. These can be
        loaded in by using Pyvista3DPlot.loadSettings
        
        Parameter
        ---------
        filename: str
            Name and path where settings are written to in json format. filename
            SHOULD CONTAIN THE SUFFIX .json
        """
        assert filename.endswith('.json'), 'Filename should contain suffix ".json".'
        path = Path(filename)
        if path.exists():
            print(f'WARNING! Given file exists already, it will be rewritten!')
        with open(path,'w') as f:
            json.dump(asdict(self.settings), f, indent='\t', sort_keys=False)
    

    def saveToVTK(self, filename):
        """
        Saves the current mesh into a .vtk file. File is saved to outputdir and 
        filename should not contain a suffix.
        """
        assert '.' not in filename, "Filename contains '.'. Note filename should not contain a suffix."
        path = self.outputdir / f'{filename}.vtk'
        self.mesh.save(path)


    def loadSettings(self, settings):
        """
        Loads Plot3DSettings to the Pyvista3DPlotter object from the given json 
        file.

        Parameters
        ----------
        settings: str
            Path to json file containing the Plot3DSettings. This should be written
            by the method Pyvist3DPlot.writeSettings.
        """
        # Path given to a file containing settings
        if isinstance(settings, str):
            if self.debug:
                print(f'String given as Plot3DSettings, interpreting as path to settings'
                        'json')
            assert settings.endswith('.json'), 'If parameter settings is a string ' \
            'it should contain suffix ".json"'      

            path = Path(settings)
            assert path.exists(), f'Given path to settings does not exists: {str(path)}'
            with open(path, 'r') as f:
                settings_dict = json.load(f)
                self.settings = Plot3DSettings(**settings_dict)
        else:
            assert isinstance(settings, Plot3DSettings), "Parameter settings should" \
                "either be path to .json file containing settings or a Plot3DSettings" \
                "object!"
            self.settings = settings
        
        # Handle parameter test: i.e. settings of format key -> list.
        settings_dict = asdict(self.settings)
        # This is a problematic parameter since it could be inputted as a list
        window_size = settings_dict.pop('window_size')
        
        self.param_test = False
        # If even one element in settings is a list, then it is a parameter test
        for value in settings_dict.values():
            if isinstance(value, list):
                self.param_test = True
                break
        
        # Handle parameter tests
        if self.param_test:
            # settings may still contain non list elements due to default values,
            # need to change this so that all values are lists
            for key in settings_dict.keys():
                if not isinstance(settings_dict[key], list):
                    settings_dict[key] = [settings_dict[key]]
            
            # Note, it is now a list of list, now ParameterGrid should handle it
            # correctly
            settings_dict['window_size'] = [window_size]
                    
            # Basically a list of settings combinations
            self.settings_grid = ParameterGrid(settings_dict)
            print(f'NOTE! parameter test run detected! Number of settings'
                f' combinations: {len(self.settings_grid)}')

           # Parameter logging
            self.log_path = self.outputdir / 'PARAMETER_TEST.csv'


    def updateSettings(self, **kwargs):
        """
        For updating any key-value pair in the dataclass Plot3DSettings. Multiple
        values can be passed at a time. 
        
        Example usage: 
        >>> plotter = Pyvista3DPlot()
        >>> plotter.updateSettings(n_points=1, mesh_cmap='bwr')
        """
        self.settings.update(**kwargs)


    def streamlines(self):
        """
        Add streamlines to plotter.

        Note that this can be used in conjunction with different source point 
        methods defined in Plot3DSettings.method and also the scalars can be 
        added to the vectors in different ways defined by Plot3DSettings.mesh_type 
        (e.g. just surfaces or volume rendering etc.)

        streamlines support parameter tests, see docstring of Pyvista3DPlot for
        explanation.
        """
        if self.param_test:
            with tqdm(total=len(self.settings_grid), desc='STREAMLINES:') as pbar:
                for settings in self.settings_grid:
                    self.__streamlines(Plot3DSettings(**settings))
                    self.__writeParamtestLog(settings, id=self.id_counter)
                    self.id_counter += 1
                    pbar.update(1)
        else:
            self.__streamlines(self.settings)
            self.id_counter += 1


    def vectors(self):
        """
        Add vectors to plotter. 
        
        Note that this can be used in conjunction with different source point 
        methods defined in Plot3DSettings.method and also the scalars can be 
        added to the vectors in different ways defined by Plot3DSettings.mesh_type 
        (e.g. just surfaces or volume rendering etc.)

        streamlines support parameter tests, see docstring of Pyvista3DPlot for
        explanation.
        """
        if self.param_test:
            with tqdm(total=len(self.settings_grid), desc='VECTORS:') as pbar:
                for settings in self.settings_grid:
                    self.__vectors(Plot3DSettings(**settings))
                    self.__writeParamtestLog(settings, id=self.id_counter)
                    self.id_counter += 1
                    pbar.update(1)
        else:
            self.__vectors(self.settings)
            self.id_counter += 1
    

    def contour(self, values=None, isosurfaces=3,  filename=f'contour_{int(time.time())}',
                cmap=None, use_vector_magnitude_as_scalars=False):
        """
        Similar to movingIsovalue() except an image is saved. Can be used to add
        or multiple isosurfaces to the image. Opens by default an interactive window
        in which camera angle can be adjusted, after which an image is saved.

        values: list of floats
            Specific values for the isosurfaces. Each isosurface generated corresponds
            to a value given in this list. If None, automatically generates `isosurfaces`
            number of isosurfaces between [scalars.min, scalars.max]
        isosurfaces: int
            Number of isosurfaces to be generated between scalars minimum and maximum.
            This applies if parameter `values=None`. Values for the isosurfaces are
            evenly spaced between the interval.
        filename: str
            Output filename, should not contain the image format, it is added 
            automatically based on the given Plot3DSettings.imageformat.
        cmap: str
            Matplotlib compatible colormap
        use_vector_magnitude_as_scalars: bool
            If true, calculates cartesian norm of the current vector field and uses
            that as the scalars.
        """
        assert not self.param_test, 'Only vectors and streamlines support parameter tests!'

        vecs = self.mesh['vectors']
        if use_vector_magnitude_as_scalars:
            magn = np.sqrt(vecs[:,0]**2 + vecs[:,1]**2 + vecs[:,2]**2)
            self.mesh['magnitude'] = magn
            scalars = 'magnitude'
        else:
            scalars = 'scalars'

        if values == None:
            values = np.linspace(self.mesh[scalars].min(), self.mesh[scalars].max(),
                             isosurfaces)
        
        # Apply pyvista countour filter
        surfaces = [self.mesh.contour([v], scalars=scalars) for v in values]

        # Remove empty surfaces if they exist ==> this would throw an error if added
        # to the plotter. THis might be the case sometimes
        to_remove = []
        for i in range(len(surfaces)):
            if surfaces[i].n_points == 0:
                print(f'WARNING: surface has zero points, removing isovalue')
                to_remove.append(i)
        for i in to_remove:
            del surfaces[i]

        lims = [self.mesh[scalars].min(), self.mesh[scalars].max()]
        self.plotter = pv.Plotter(window_size=self.settings.window_size)
        self.__plotterSettings(self.settings)

        for surf in surfaces:
            self.plotter.add_mesh(surf, cmap=cmap, clim=lims)
        
        print('\n--> Pan around the camera to wanted angle, then press "q" to save image!\n')
        self.plotter.show(auto_close=False)
        self.plotter.screenshot(filename=f'{filename}.{self.settings.imageformat}')


    def scalars(self):
        """
        Plot only the scalar values on to the plotter.
        """
        assert not self.param_test, 'Only vectors and streamlines support parameter tests!'
        self.plotter = pv.Plotter(off_screen=self.settings.off_screen, window_size=self.settings.window_size)        
        self.__handlePlotter(self.settings)


    def movingIsovalue(self, values=None, filename=f'moving_isovalue_{time.time()}.gif',
                         isosurfaces=30, show_outline=True, show_mesh_outline=True,
                         cmap=None, use_vector_magnitude_as_scalars=False):
        """
        Creates a moving isovalue gif that cycles showing the isovalues from first
        to last and back making a "moving plot". The isovalues itself can be defined,
        otherwise an isovalue range is created by default of `isosurfaces` amount of 
        isosurfaces between minimum and maximum of scalars.

        values: list of floats
            Specific values for the isosurfaces. Each isosurface generated corresponds
            to a value given in this list. If None, automatically generates `isosurfaces`
            number of isosurfaces between [scalars.min, scalars.max]
        filename: str
            Output filename for the gif. Defaults to 'moving_isovalue_<timestamp>.gif'
        isosurfaces: int
            Number of isosurfaces to be generated between scalars minimum and maximum.
            This applies if parameter `values=None`. Values for the isosurfaces are
            evenly spaced between the interval.
        show_outline: bool
            Shows an outline box in which the mesh is contained in.
        show_mesh_outline: bool
            TODO! 
        cmap: str
            Matplotlib compatible colormap
        use_vector_magnitude_as_scalars: bool
            If true, calculates cartesian norm of the current vector field and uses
            that as the scalars.
        """
        assert not self.param_test, 'Only vectors and streamlines support parameter tests!'
        # scalars = self.mesh['scalars']
        vecs = self.mesh['vectors']
        if use_vector_magnitude_as_scalars:
            magn = np.sqrt(vecs[:,0]**2 + vecs[:,1]**2 + vecs[:,2]**2)
            self.mesh['magnitude'] = magn
            scalars = 'magnitude'
        else:
            scalars = 'scalars'

        if values == None:
            values = np.linspace(self.mesh[scalars].min(), self.mesh[scalars].max(),
                             isosurfaces)
        
        # Apply pyvista countour filter
        surfaces = [self.mesh.contour([v], scalars=scalars) for v in values]

        # Remove empty surfaces if they exist ==> this would throw an error if added
        # to the plotter. THis might be the case sometimes
        to_remove = []
        for i in range(len(surfaces)):
            if surfaces[i].n_points == 0:
                print(f'WARNING: surface has zero points, removing isovalue')
                to_remove.append(i)
        for i in to_remove:
            del surfaces[i]

        surface = surfaces[0].copy()
        plotter = pv.Plotter()
        plotter.open_gif(filename)
        plotter.enable_depth_peeling()
        
        lims = [self.mesh[scalars].min(), self.mesh[scalars].max()]
        plotter.add_mesh(surface, opacity=1.0, clim=lims, cmap=cmap)

        if show_outline:
            # plotter.add_mesh(self.mesh.outline_corners(), color='k')
            plotter.add_mesh(self.mesh.outline(), color='k')

        print('\n--> Pan around the camera to wanted angle, then press "q" to produce the movie!\n')
        plotter.show(auto_close=False)
        print(f'Starting to create the isovalue gif, this might take a moment!')
        with tqdm(total=2*len(surfaces), desc="Moving isovalue rendering:") as pbar:
            for surf in surfaces:
                surface.overwrite(surf)
                plotter.write_frame()
                pbar.update(1)
            for surf in surfaces[::-1]:
                surface.overwrite(surf)
                plotter.write_frame()
                pbar.update(1)
        plotter.close()


    def preview(self):
        """
        !TODO!
        """
        assert not self.param_test, 'Only vectors and streamlines support parameter tests!'
        if self.param_test:
            settings = self.settings_grid[0]
        else:
            settings = self.settings
            
        self.plotter = pv.Plotter(window_size=settings.window_size,
                                  title='Plot Preview')
        scalar_bar_args = {
            'width': settings.sbar_width,
            'height': settings.sbar_height,
            'vertical': settings.vertical_sbar,
        }
        self.plotter.add_mesh(self.surface_mesh, 
                                  opacity=settings.mesh_opacity, 
                                  cmap=settings.mesh_cmap,
                                  show_scalar_bar=settings.show_mesh_sbar,
                                  scalar_bar_args=scalar_bar_args,
                                  )
        self.plotter.show_bounds(color='black', location='outer')
        
        print("--> NOTE! In spreview only the mesh is added by default, not vectors | streamlines are shown!")

        plotPreview(self.plotter)
        
     
    ############################################################################
    ### INTERNALS ##############################################################
    ############################################################################


    def __command_help_print(self):
        help = """
        |-------------------------------------------------------
        | You can pan around the camera and when closing the plot (using 'q')
        | it saves an image at the chosen camera location (assuming you set the imageformat
        | parameter to appropriate file type).
        |-------------------------------------------------------
        |You have interactive plotting on! For controls see:
        |
        |--> https://docs.pyvista.org/plotting/plotting.html?highlight=plotter#pyvista.BasePlotter.update
        |
        | Most important controls:
        |------------------------
        |* q --- Close the rendering window
        |* v --- Isometric camera
        |* w --- Switch all datasets to wireframe presentations
        |* r --- Reset camera
        |* Left-click --- Rotate scene in 3D
        |* CTRL+click --- Rotate scene in 2D
        |* SHIFT+S --- Screenshot (only on BackgroundPlotter)
        |* +/- --- increase point size and line widths
        |---------------------------------------------------------
        """
        print(help)


    def __streamlines(self, settings: Plot3DSettings):
        assert settings.method in SAMPLING_METHODS, f"Unkown sampling method {settings.method}."
        f"Should be one of: {SAMPLING_METHODS}"
        
        if settings.surface_streamlines and settings.method == 'points':
            print(f'\nWARNING! surface_streamlines and sampling method "points" doesnt make sense'
                  'when used together! Rather use "surface points" or one of the slicing methods!\n')
            
        self.plotter = pv.Plotter(off_screen=settings.off_screen, window_size=settings.window_size)
        
        if settings.method == 'default':
            if settings.surface_streamlines:
                raise Warning('Surface streamlines cannot be done with "default"'
                                  'sampling. Use "points".')
            
            stream, self.src = self.mesh.streamlines('vectors', return_source=True, 
                                    n_points=settings.n_points, 
                                    source_radius=settings.source_radius, 
                                    source_center=settings.source_center,
                                    **settings.stream_params)
            self.field = [stream.tube(radius=settings.tube_radius,
                                        scalars=settings.stream_variable_radius,
                                        radius_factor=settings.stream_radius_factor)]
            self.src = [self.src]
        else:
            self.src = meshSamplingMethods(
                    settings.method, self.mesh, self.surface_mesh, settings.n_points, 
                    normal=settings.normal, origin=settings.origin, set_seed_1=settings.set_seed_1, 
                    alwaysReturnList=True)
            
            self.field = []
            for source in self.src:
                if settings.surface_streamlines:
                    mesh = self.surface_mesh
                else:
                    mesh = self.mesh
                
                stream = mesh.streamlines_from_source(source, vectors='vectors', 
                            surface_streamlines=settings.surface_streamlines, 
                            **settings.stream_params)
                
                self.field.append(stream.tube(radius=settings.tube_radius,
                                                scalars=settings.stream_variable_radius,
                                                radius_factor=settings.stream_radius_factor))

        for field in self.field:
            printTerminationReasons(field['ReasonForTermination'], self.plotter)        

        self.__handlePlotter(settings)


    def __vectors(self, settings: Plot3DSettings):
        assert settings.method != 'default', f"Sampling method default is only" \
            "available for streamlines!"
            
        self.plotter = pv.Plotter(off_screen=settings.off_screen, 
                                  window_size=settings.window_size)
        
        if settings.vector_scaling == 'magnitude':
            magnitude = self.mesh['vector magnitude']
            # If vector max not set, use mean + one STD of magnitudes as a reasonable value
            if settings.vector_max is None:
                settings.vector_max = magnitude.mean() + magnitude.std()
            self.mesh['scale'] = np.array([settings.vector_max if a > settings.vector_max 
                                    else a for a in magnitude])
            scale = 'scale'
        elif settings.vector_scaling == 'scalars':
            scale = 'scalars'
        elif settings.vector_scaling is None:
            scale = None
        else:
            raise ValueError(f'[__add_vectors] Unknown scaling method scaling = {settings.vector_scaling}')
        
        self.src = meshSamplingMethods(
            settings.method, self.mesh, self.surface_mesh, n_points=settings.n_points, 
            normal=settings.normal, origin=settings.origin, set_seed_1=settings.set_seed_1, 
            alwaysReturnList=True, get_arrays=True)
        
        self.field = []
        for source in self.src:
            self.field.append(source.glyph(orient='vectors', 
                                           scale=scale, 
                                           factor=settings.vector_factor))
 
        self.__handlePlotter(settings)


    def __3D_data_loader(self, scalar_key, vector_key):
        """
        Reads in a var file, generates a mesh and adds the given vector field and scalar
        field to the mesh. That is, output mesh contains point arrays 'vectors' and 
        'scalars' that corresponds to the given input fields.
        
        vector_key can be either a single string e.g. 'uu' or a list of vector field 
        keys e.g. ['uu1', 'uu2', 'uu3']. 
        """
        if isinstance(vector_key, list):
            assert len(vector_key) == 3, f'Parameter vector_key should contain 3 elements'
            f'if list is given. Given vector_key contains: {vector_key}'
                    
        # Read scalar data
        scalars = self.var.__getattribute__(scalar_key)
        # Or this works too: previously was like this but read order was "C" not "F"
        # scalars = self.var.__getattribute__(scalar_key).swapaxes(0,2).flatten(order='F')
        nz, ny, nx = scalars.shape
        self.nx, self.ny, self.nz = nx, ny, nz

        # Note that the vectors is of shape [points, 3], this is given to the pyvista mesh 
        vectors = np.empty((nx*ny*nz, 3))

        if self.coordinates == 'cartesian':  
            scalars = scalars.flatten()

            # Create the mesh, note only UniformGrid can be used for volume rendering
            origx, origy, origz = self.var.x.min(), self.var.y.min(), self.var.z.min()
            self.mesh = pv.UniformGrid()
            self.mesh.dimensions = (self.nx, self.ny, self.nz)
            self.mesh.origin = (origx, origy, origz)
            self.mesh.spacing = (self.grid.dx, self.grid.dy, self.grid.dz)

            if isinstance(vector_key, list):
                # Vectors given as [z,y,x] thus swap first and last axes and reshape to vector of
                # shape [points,]
                # vectors[:,0] = self.var.__getattribute__(vector_key[0]).swapaxes(0,2).reshape(-1)
                # vectors[:,1] = self.var.__getattribute__(vector_key[1]).swapaxes(0,2).reshape(-1)
                # vectors[:,2] = self.var.__getattribute__(vector_key[2]).swapaxes(0,2).reshape(-1)
                vectors[:,0] = self.var.__getattribute__(vector_key[0]).flatten()
                vectors[:,1] = self.var.__getattribute__(vector_key[1]).flatten()
                vectors[:,2] = self.var.__getattribute__(vector_key[2]).flatten()
            else:
                # Now input is of form [3, z, y, x] ==> output needs to be of shape [points, 3] 
                # vectors = self.var.__getattribute__(vector_key).swapaxes(1,3
                #                                     ).reshape(3,-1).swapaxes(0,1)
                vectors = self.var.__getattribute__(vector_key).reshape(3,-1).swapaxes(0,1)

            # x, y, z = self.var.x, self.var.y, self.var.z
            # xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
            # self.mesh = pv.StructuredGrid(xx, yy, zz) 

        elif self.coordinates == 'spherical':  
            scalars = scalars.T.flatten()

            if isinstance(vector_key, list):
                u = self.var.__getattribute__(vector_key[0]) 
                v = self.var.__getattribute__(vector_key[1])
                w = self.var.__getattribute__(vector_key[2])
            else:
                tmp = self.var.__getattribute__(vector_key)
                u, v, w = tmp[0,:], tmp[1,:], tmp[2,:]
                
            r = self.var.x 
            phi = self.var.z * 180 / np.pi
            theta = self.var.y * 180 / np.pi

            self.mesh = pv.grid_from_sph_coords(phi, theta, r)

            u_t, v_t, w_t = pv.transform_vectors_sph_to_cart(phi, theta, r, u, v, w)
            
            u_t = u_t.swapaxes(0,2).reshape(-1)
            v_t = v_t.swapaxes(0,2).reshape(-1)
            w_t = w_t.swapaxes(0,2).reshape(-1)
            
            vectors[:,0] = u_t
            vectors[:,1] = v_t
            vectors[:,2] = w_t
            
        else:
            raise ValueError(f'Coordinate system {self.coordinates} is not a valid coordinate system!')
        
        if self.settings.scalars_norm == 'log':
            scalars = np.log10(scalars)
        
        # Calculate vector magnitude
        self.mesh['vector magnitude'] = np.sqrt(vectors[:,0]**2 + vectors[:,1]**2 + vectors[:,2]**2)
        
        # TO REMOVE
        self._scalars = scalars.reshape(nx,ny,nz)

        self.mesh['scalars'] = scalars
        self.mesh['vectors'] = vectors
        self.mesh.set_active_scalars('scalars')
        self.mesh.set_active_vectors('vectors')


    def __addWidgets(self, settings):
        print(f'>>> WARNING! with larger datasets some of these widgets might be'
        ' laggy! Try "clip plane" might be generally faster than other widgets.')
        
        field_sbar_args = create_sbar_args(settings, settings.field_sbar_title,
                                settings.field_sbar_pos_x, settings.field_sbar_pos_y) 
        scalar_sbar_args = create_sbar_args(settings, settings.scalar_sbar_title,
                                settings.scalar_sbar_pos_x, settings.scalar_sbar_pos_y)

        if settings.widget_type == 'clip slice':
            self.plotter.add_mesh_slice(self.mesh, scalars='scalars', 
                                        scalar_bar_args=scalar_sbar_args)
        
        elif settings.widget_type == 'clip box':
            self.plotter.add_mesh_clip_plane(self.mesh, scalars='scalars', 
                                            scalar_bar_args=scalar_sbar_args)
        
        elif settings.widget_type == 'plane vectors':
            print(f'> NOTE! "plane vectors" widget can be tuned using Plot3DSettings vector parameters!')
            # Callback function is called by add_plane_widget
            def callback(normal, origin):
                slice = self.mesh.slice(normal=normal, origin=origin)
                
                if settings.vector_scaling == 'magnitude':
                    magnitude = slice['vector magnitude']
                    # If vector max not set, use mean + one STD of magnitudes as a reasonable value
                    if settings.vector_max is None:
                        settings.vector_max = magnitude.mean() + magnitude.std()
                    slice['scale'] = np.array([settings.vector_max if a > settings.vector_max 
                                            else a for a in magnitude])
                    scale = 'scale'
                elif settings.vector_scaling == 'scalars':
                    scale = 'scalars'
                elif settings.vector_scaling is None:
                    scale = None
                else:
                    raise ValueError(f'[__add_vectors] Unknown scaling method scaling = {settings.vector_scaling}')
                
                # Random sample the slice so that not all vectors are plotted
                if self.settings.n_points != None:
                    sampled = randomSampleMeshPoints(settings.n_points, slice, get_arrays=True)
                else: 
                    sampled = slice

                arrows = sampled.glyph(orient='vectors', scale=scale, factor=settings.vector_factor)

                # if self.coordinates == 'cartesian':
                #     origin = 

                self.plotter.add_mesh(self.surface_mesh, opacity=0.2, show_scalar_bar=False)
                self.plotter.add_mesh(slice, opacity=settings.mesh_opacity, name='scalars',
                                    cmap=settings.mesh_cmap, scalar_bar_args=scalar_sbar_args)
                self.plotter.add_mesh(arrows, opacity=settings.field_opacity, name="arrows", 
                                    cmap=settings.field_cmap, scalar_bar_args=field_sbar_args)

            self.plotter.add_plane_widget(callback,)

        else:
            raise ValueError(f'Unknown widget type: {settings.widget_type}.')

        self.__plotterSettings(settings)
        self.plotter.show(window_size=self.settings.window_size, auto_close=False)

        if settings.imageformat is not None:
            self.plotter.screenshot(self.outputdir / f'widget_image_{int(time.time())}.{settings.imageformat}')

        self.plotter.close()


    def __plotterSettings(self, settings):
        self.plotter.add_text(text=settings.title, 
                              position='upper_left',
                              font_size=settings.title_font_size,
                              color=settings.annotation_color)
        
        if self.settings.show_axes:
            self.plotter.show_bounds(color=settings.annotation_color, location='outer')
        if self.settings.show_grid:
            self.plotter.show_grid()

        self.__set_camera(settings)
        self.plotter.background_color = settings.background_color


    def __handlePlotter(self, settings):
        """
        Adds all necessary meshes to the plotter and shows it / saves an image or
        optionally generates an orbit gif. These meshes include streamlines or 
        vectors (called field here), source points for the streamlines / vectors
        and the mesh itself the vector field is in. 
        """
        field_sbar_args = create_sbar_args(settings, settings.field_sbar_title,
                                settings.field_sbar_pos_x, settings.field_sbar_pos_y) 
        scalar_sbar_args = create_sbar_args(settings, settings.scalar_sbar_title,
                                settings.scalar_sbar_pos_x, settings.scalar_sbar_pos_y) 

        if settings.widget_type != None and isinstance(settings.widget_type, str):
            self.__addWidgets(settings)

            # We want to stop execution of __handlePlotter here. Otherwise it'll
            # try to add meshes to an already closed plotter
            return
        
        # Add scalars to plotter
        if settings.show_mesh:
            if settings.mesh_type == "volume":
                assert self.coordinates == 'cartesian', "At the moment only cartesian supports volume rendering!"
                
                # Possible mappers for add_volume are: 'fixed_point', 'gpu', 
                #   'open_gl', 'smart'
                va = self.plotter.add_volume(self.mesh, # self._scalars
                                        scalars=self._scalars,
                                        scalar_bar_args=scalar_sbar_args,
                                        show_scalar_bar=settings.show_mesh_sbar,
                                        cmap=settings.mesh_cmap,
                                        opacity='linear',
                                        mapper='smart',
                                        opacity_unit_distance= self.mesh.length * settings.volume_opacity,
                                        # shade=True
                                        )
                # This adds a slider to the plot moving the opacity_unit_distance parameter
                if settings.add_volume_opacity_slider:
                    f = lambda val: va.GetProperty().SetScalarOpacityUnitDistance(self.mesh.length * val)
                    self.plotter.add_slider_widget(f, [0, 1], 
                                title="Volume opacity")

            elif settings.mesh_type == "surface":
                self.plotter.add_mesh(self.surface_mesh, 
                                    opacity=settings.mesh_opacity, 
                                    cmap=settings.mesh_cmap,
                                    show_scalar_bar=settings.show_mesh_sbar,
                                    scalar_bar_args=scalar_sbar_args,
                                    culling=settings.culling)
            elif settings.mesh_type == "orthogonal slices":
                x, y, z = settings.slice_pos
                slices = self.mesh.slice_orthogonal(x=x, y=y, z=z)
                self.plotter.add_mesh(slices, 
                                    opacity=settings.mesh_opacity, 
                                    cmap=settings.mesh_cmap,
                                    show_scalar_bar=settings.show_mesh_sbar,
                                    scalar_bar_args=scalar_sbar_args)
            else:
                raise ValueError(f"Unknown mesh_type == {settings.mesh_type}!"
                " Options: 'volume', 'surface' or 'orthogonal slices'.")
                
        # Add vectors | streamlines to plotter
        if self.field != None:
            for tmp in self.field:
                if tmp.n_points == 0:
                    print(f'WARNING! This field (vectors | streamlines) has zero points'
                        f' nothing to plot ---> skipping adding field to the plot.')
                else:
                    self.plotter.add_mesh(tmp, 
                                        opacity=settings.field_opacity,
                                        cmap=settings.field_cmap,
                                        show_scalar_bar=settings.show_field_sbar,
                                        scalar_bar_args=field_sbar_args)
                    
            # Show vector | streamline source points
            if settings.show_source:
                for source in self.src:
                    self.plotter.add_mesh(source, render_points_as_spheres=True)
        
        # Apply extra settings on the plotter + annotations
        self.__plotterSettings(settings)

        if settings.orbit_gif:
            print(f'Starting orbit gif generation...')
            orbit = self.plotter.generate_orbital_path(n_points=settings.orbit_points, 
                                           shift=self.mesh.length)
            path = self.outputdir / 'streamline_orbit.gif'
            # open_gif() doesnt work with the Path() object ==> need to give as a string
            # to circumvent this issue
            self.plotter.open_gif(str(path))
            self.plotter.orbit_on_path(orbit, step=settings.orbit_step, write_frames=True)
            self.plotter.close()
            print(f'Orbit gif generation finished!')
            return
        
        if self.debug:
            print(f'Camera angle: {self.plotter.camera_position[0]}')

        # Plotting not off screen ==> show live viewer
        if not settings.off_screen:
            self.__command_help_print()
            self.plotter.show(title=f'Scalar: {self.scalar_key}, vector: {self.vector_key}',
                                auto_close=False)
        
        if settings.imageformat is not None:
            self.plotter.screenshot(
                filename=self.outputdir / f'{self.id_counter:04d}_streamline_{int(time.time())}.{settings.imageformat}')
        
        self.plotter.close()


    def __set_camera(self, settings):
        """
        Set plotter camera using all or some of the possible supplied values. That is:
        - settings.camera_centre
        - settings.focal_point
        - settings.view_up 
        """
        if settings.camera_centre != None and settings.focal_point != None:
            self.plotter.camera_position = (settings.camera_centre, 
                                            settings.focal_point, 
                                            settings.view_up)
        elif settings.camera_centre != None:
            _, focal, _ = self.plotter.camera_position
            self.plotter.camera_position = (settings.camera_centre,
                                            focal,
                                            settings.view_up)
        elif settings.focal_point != None:
            pos, _, _ = self.plotter.camera_position
            self.plotter.camera_position = (pos,
                                            settings.focal_point,
                                            settings.view_up)
        if settings.zoom_factor != None:
            self.plotter.camera.zoom(settings.zoom_factor)


    def __writeParamtestLog(self, settings: dict, id=-1):
        # Note we assume the settings is passed as a dict NOT AS A DATACLASS,
        # sort the dict so that keys-values are always in the same order in the log
        s_dict = dict(sorted(settings.items(), key=lambda x: x[0]))

        writeHeader = False
        if not self.log_path.exists():
            writeHeader = True

        with open(self.log_path, 'a+', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            if writeHeader:
                writer.writerow(['id'] + list(s_dict.keys()))

            writer.writerow([f'{id}'] + list(s_dict.values()))


    ############################################################################