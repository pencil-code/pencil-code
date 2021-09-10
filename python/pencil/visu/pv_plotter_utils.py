#
# pv_plotter_utils.py
# 
# This file defines generic utility routines for PyVista plotters pv_plotter.py 
# and pv_volume_plotter.py.
#
# Author: Leevi Veneranta (leevi.veneranta@aalto.fi)
#
"""
Requirements
------------
See ``pv_plotter.py`` docstring!
"""
# External Python libraries
from tqdm import tqdm
import pyvista as pv
import numpy as np
# import cv2

# Internal Python libraries
from pathlib import Path
import os
import time


def plotPreview(plotter) -> None:
    """
    Routine for interactive plot preview window showing the camera parameters.
    Parameters
    ----------
    plotter : pv.Plotter
        Instance of pv.Plotter
    window_size : list, optional
        Windows size [width, height], by default None
        
    Notes
    -----
    Last updated: 18.6.2021
    """
    raise NotImplementedError("Currently not working properly! Set parameter preview=False")

    plotter.show(title='Plot preview', interactive=True,
                 auto_close=False, interactive_update=True)

    cam_prev = (0,0,0)
    counter = 0
    print(f'WARNING! In order to close window properly, type "CTRL-C" on the'
          'command line!')
    print(f'>>> DO NOT PRESS quit on the plot window!')
    while True:       
        try:
            cam, focal, viewup = getCameraParamsRounded(plotter)

            if cam != cam_prev:
                print(f'[{counter}] CAMERA PARAMS:\n> camera centre = {cam}' 
                    f'\n> focal point = {focal}\n> view up = {viewup}\n')
                plotter.add_text(f'CAMERA PARAMETERS:\n> Camera centre = {cam}' 
                                f'\n> Focal point = {focal}\n> View up = {viewup}', 
                                name='title', font_size=11)
                counter += 1
                cam_prev = cam

            plotter.update()
            
        except KeyboardInterrupt:
            # When CTRL + C pressed
            plotter.close()


def randomSampleMeshPoints(n, mesh, set_seed_1=False, constant_seed=False,
                             get_arrays=False, seed=None):
    """
    Returns randomly (uniformly) sampled points from the given mesh, 
    i.e. a subset of mesh.points
    
    Parameters
    ----------
    n : int
        Number of points to choose from the mesh
    mesh
        Pyvista mesh
    set_seed_1 : bool, optional
        If set_seed_1 set to True, sets seed used by the random generator to 1, 
        i.e. each time.time this is called with the same mesh, function would return
        the same points. For debugging purposes. By default False
    get_arrays : bool, optional
        If True, also samples all the point arrays that the given mesh contains,
        and adds them to the returned pyvista.PolyData
        
    Returns
    -------
    pyvista.PolyData
        Set of n randomly sampled points from the given mesh.
    """
    npoints = mesh.n_points
    assert n <= npoints, f"Number of points chosen must be at most mesh.n_points" \
     f"= {mesh.n_points}"
    
    # For debugging purposes to keep point choices constant  
    if set_seed_1:
        np.random.seed(1)
    elif constant_seed:
        assert seed is not None, "INTERNAL ERROR: seed value not set?"
        np.random.seed(seed)
        
    idx = np.random.randint(0, high=npoints, size=n)
    points = mesh.points[idx,:]
    sampled = pv.PolyData(points)
    
    if get_arrays:
        for key in mesh.array_names:
            sampled[key] = mesh[key][idx]
        
    return sampled


def randomSampleMesh(n, mesh):
    """
    PROBLEM: random samples individual cells which consists of 9 points. Hence ends up being
    individual cells all around the plot each having 9 streamlines
    """
    npoints = mesh.n_points
    assert n <= npoints, f"Number of points chosen must be at most mesh.n_points = {mesh.n_points}"
    idx = np.random.randint(0,high=npoints, size=n)
    mask = np.zeros(npoints)
    mask[idx] = 1
    mesh['mask'] = mask
    return mesh.threshold(value=1, scalars='mask')


def getCameraParamsRounded(plotter: pv.Plotter, precision=3) -> tuple:
    """
    Returns camera parameters rounded to given precision.
    Parameters
    ----------
    plotter : pv.Plotter
        Instance of pv.Plotter
    precision : int, optional
        Precision to round to, by default 3
    Returns
    -------
    tuple
        (camera_centre, focal_point, view_up)
        
    Notes
    -----
    Last updated: 18.6.2021
    """
    cam, focal, viewup = plotter.camera_position
    
    cam = (round(cam[0], precision),
           round(cam[1], precision),
           round(cam[2], precision))
    
    focal = (round(focal[0], precision),
             round(focal[1], precision),
             round(focal[2], precision))
    
    viewup = (round(viewup[0], precision),
              round(viewup[1], precision),
              round(viewup[2],precision))
    
    return cam, focal, viewup


def gridFromCylCoords(r, theta, z) -> pv.StructuredGrid:
    """"
    Create PyVista StructuredGrid out of cylinder coordinates.
    Parameters
    ----------
    r : numpy.ndarray
        Radial distances
    theta : numpy.ndarray
        Polar angles (in radians)
    z : numpy.ndarray
        Z-coordinates (heights)
    Returns
    -------
    pyvista.StructuredGrid
        Mesh given the coordinates
        
    Notes
    -----
    Last updated: 18.6.2021
    """
    x,y,z = np.meshgrid(r, theta, z)
    # Transform grid to cartesian coordinates
    x_cart = x * np.cos(y)
    y_cart = x * np.sin(y)
    z_cart = z
    # Make a grid object
    return pv.StructuredGrid(x_cart, y_cart, z_cart)


def cellBounds(points, bound_position=0.5) -> np.ndarray:
    """
    Calculate coordinate cell boundaries.
    Parameters
    ----------
    points: numpy.array
        One-dimensional array of uniformly spaced values of shape (M,)
    bound_position: bool, optional
        The desired position of the bounds relative to the position
        of the points.
    Returns
    -------
    bounds: numpy.array
        Array of shape (M+1,)
    Examples
    --------
    >>> a = np.arange(-1, 2.5, 0.5)
    >>> a
    array([-1. , -0.5,  0. ,  0.5,  1. ,  1.5,  2. ])
    >>> cell_bounds(a)
    array([-1.25, -0.75, -0.25,  0.25,  0.75,  1.25,  1.75,  2.25])
    """
    assert points.ndim == 1, "Only 1D points are allowed"
    diffs = np.diff(points)
    delta = diffs[0] * bound_position
    bounds = np.concatenate([[points[0] - delta], points + delta])
    return bounds


def create_sbar_args(settings, title, posx, posy):
    """
    Small convenience function to remove duplicate code. Combines all scalarbar
    arguments into one dictionary. Note that PlotSettings __post_init__ creates
    _sbar_args but it is still missing parameters that depend on whether one is
    showing field (stream | vectors) or scalars.
    """
    sbar_args = {
                'title': title,
                'position_x': posx,
                'position_y': posy,
                **settings._sbar_args
                }
    if title is None or title == '':
        # The title works stupidly - if you want pyvista to infer sbar name 
        # automatically from point arrays, you cant pass title parameter at all
        # thus if it is None ==> remove title from sbar_args
        del sbar_args['title']
    
    return sbar_args


# def imgsToVideo(output='test_movie.mp4', imgpath=Path('./images'), outputdir=Path('./'), 
#                 imageformat='.png', fps=10, fourcc=0, debug=False):
#     """
#     NOTE! This seems to have an effect on the image quality, unlike running just
#     ffmpeg on command line.
#     Combine images in directory imgpath to create video named output to outputdir.
#     Images read in are filtered based on imageformat.
#     Note, that the default value for fourcc=0 works at least for creating '.avi'
#     format videos. If '.mp4' are wanted, script uses by default the codec 'mp4v'
#     which should work.
#     Parameters
#     ----------
#     output: str
#         Output filename with videoformat included, e.g. 'test.mp4' or 'test.avi'.
#     imgpath: str or pathlib.Path object
#         Directory to images.
#     outputdir: str or pathlib.Path object
#         Path of output directory where output is saved to.
#     imageformat: str
#         Imageformat, e.g. '.png'
#     fps: int
#         Frames per second for the video.
#     fourcc: int or str
#         From OpenCV doc: 4-character code of codec used to compress the frames. 
#         For example, CV_FOURCC('P','I','M','1') is a MPEG-1 codec, 
#         CV_FOURCC('M','J','P','G') is a motion-jpeg codec etc. List of codes can 
#         be obtained at Video Codecs by FOURCC page.
#         Link: https://www.fourcc.org/codecs.php
#     debug: bool
#         Enable debug printing.
#     Notes
#     -----
#     Seems to be a fairly fast method, writes approx 40 frames/s. Creating a movie
#     out of ~1400 png images takes about 30 seconds.
#     """
#     if debug:
#         T0 = time.time()
#     # Optionally string codec can be passed in insted of an integer.
#     if isinstance(fourcc, str):
#         fourcc = cv2.VideoWriter_fourcc(*fourcc)
#     # If output is mp4 ==> user 'mp4v' codec, seems to work, see:
#     # https://stackoverflow.com/questions/60689660/cant-save-a-video-in-opencv
#     if output.endswith('.mp4'):
#         fourcc = cv2.VideoWriter_fourcc(*'mp4v')
#     if isinstance(outputdir, str):
#         outputdir = Path(outputdir)
#         if not outputdir.exists():
#             print(f'> Output directory does not exist! Creating new directory')
#             os.mkdir(outputdir)
#     if isinstance(imgpath, str):
#         imgpath = Path(imgpath)
#     images = [img for img in os.listdir(imgpath) if img.endswith(imageformat)]  
#     frame = cv2.imread(str(imgpath / images[0]))
#     height, width, layers = frame.shape
#     video = cv2.VideoWriter(str(outputdir / output), fourcc, fps, (width, height))
#     with tqdm(total=len(images), desc='CREATING VIDEO', unit='frame') as pbar:
#         for image in images:
#             video.write(cv2.imread(str(imgpath / image)))
#             pbar.update(1)
#     cv2.destroyAllWindows()
#     video.release()
#     if debug:
#         T1 = time.time()
#         print(f'Creating video took: {T1-T0} seconds.')
#         print(f'Output written to: {str(outputdir / output)}')
#         print(f'Output exists: {(outputdir / output).exists()}')