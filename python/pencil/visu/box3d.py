# box3d.py
#
# Plotting routines for a 3d box and its animation (from slides).
#
# Author: S. Candelaresi (iomsn1@gmail.com).
"""
Contains the classes and methods to read the slice files,
generate box plot and an animation with the simulation box.
"""


def box3d(slices=None):
    """
    Plot the given slices as a 3d box.

    call signature:

    box3d()

    Keyword arguments:

    *slices*:
      A list of 3 or 4 slices, i.e. 2d arrays of shape [m, n].
      The first slices is the top, then left, right and bottom.
    """

    import inspect

    # Assign parameters to the AnimateBox object.
    box3d_return = Box3d()
    argument_dict = inspect.getargvalues(inspect.currentframe()).locals
    for argument in argument_dict:
        setattr(box3d_return, argument, argument_dict[argument])

    box3d_return.plot()

    return box3d_return


class Box3d(object):
    """
    Box plotting class in 3d.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.slices = None
        
        self.axes = True
        self.font_size = 24
        self.font_family = 'arial'
        self.title = ''
        self.colorbar = True
        self.colormap = 'binary'
        self.colorbar_label = r'u_x'
        self.xlabel = r'$x$'
        self.ylabel = r'$y$'
        self.zlabel = r'$z$'
        self.vmin = None
        self.vmax = None


    def plot(self):
        """
        Plot the 3d box.
        """

        import pylab as plt
        from matplotlib.transforms import Affine2D
        import mpl_toolkits.axisartist.floating_axes as floating_axes

        # Prepare the plot.
        width = 6
        height = 5
        plt.rc('text', usetex=True)
        plt.rc('font', family='arial')
        plt.rc("figure.subplot", left=0.15)
        plt.rc("figure.subplot", right=.90)
        plt.rc("figure.subplot", bottom=.15)
        plt.rc("figure.subplot", top=.95)
        fig = plt.figure(figsize=(width, height))

        # Top slice.
        data = self.slices[0]
        transformation = Affine2D().scale(0.5, 0.5).rotate_deg(45).translate(0, 0.5)
        grid_helper = floating_axes.GridHelperCurveLinear(transformation,
                                                          extremes=(0, 1, 0, 1))
        axes = floating_axes.FloatingSubplot(fig, 111, grid_helper=grid_helper)
        fig.add_subplot(axes)
        aux_axes = axes.get_aux_axes(transformation)
        aux_axes.imshow(data, extent=[0, 1, 0, 1])

        # Left slice.
        data = self.slices[1]
        transformation = Affine2D().scale(0.5, 0.5).rotate_deg(45).translate(0, 0.5).skew_deg(20, 40)
        grid_helper = floating_axes.GridHelperCurveLinear(transformation,
                                                          extremes=(0, 1, 0, 1))
        axes = floating_axes.FloatingSubplot(fig, 121, grid_helper=grid_helper)
        fig.add_subplot(axes)
        aux_axes = axes.get_aux_axes(transformation)
        aux_axes.imshow(data, extent=[0, 1, 0, 1])
       
        if len(self.slices) > 3:
            pass
        
#def animate_box(*args, **kwargs):
#    """
#    Read the appropriate slice files and generate and animation.
#    Return object containing the plotted data.
#
#    call signature:
#
#    animate_box()
#
#    Keyword arguments:
#
#    *aaa*:
#      aaa
#
#    *bbb*:
#      bbb
#    """
#
#    import inspect
#
#    # Assign parameters to the AnimateBox object.
#    animate_box_return = AnimateBox()
#    argument_dict = inspect.getargvalues(inspect.currentframe()).locals
#    for argument in argument_dict:
#        setattr(animate_box_return, argument, argument_dict[argument])
#
#    animate_box_return.read()
#    animate_box_return.make_animation()
#
#    return animate_box_return
#    
#
#class AnimateBox(Box3d):
#    """
#    AnimateBox -- holds the parameters, slice data and methods for the box animation.
#    """
#
#    def __init__(self):
#        """
#        Fill members with default values.
#        """
#
#        self.t = None
#        self.slices = None
#        self.datadir = 'data'
#        self.t_min = 0
#        self.t_max = 1
#        self.field = 'uu1'
#        
#        self.axes = True
#        self.font_size = 24
#        self.font_family = 'arial'
#        self.title = ''
#        self.colorbar = True
#        self.colormap = 'binary'
#        self.colorbar_label = r'u_x'
#        self.xlabel = r'$x$'
#        self.ylabel = r'$y$'
#        self.zlabel = r'$z$'
#        self.vmin = None
#        self.vmax = None
#        self.test = False
#
#        self.save_png = False
#        self.save_eps = False
#        self.dpi = 300
#
#
#    def read(self):
#        """
#        Read the slice files.
#        """
#
#        from pencil.read import slices
#
#        self.slices = slices(field=self.field, extension=['xy', 'xz', 'yz', 'xy2'],
#                             datadir=self.datadir)
#
#        return 0
#
#
#    def make_animation(self):
#        """
#        Create the animation.
#        """
#
#        import numpy as np
#        import pylab as plt
#
#        # Create the plotting structure and plot the first slice.
#
#
#        
#        for t_idx in np.arange(self.t_min+1, self.t_max+1):
#            
#            
#            
